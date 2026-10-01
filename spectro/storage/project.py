"""Project file (SQLite) with immutable raw data and a hash-chained audit trail.

Data-integrity rules enforced here and by database triggers:

* Spectral data (x/y arrays, hashes, lineage) can never be updated or deleted;
  records are archived instead.
* The audit table is append-only; each entry stores the SHA-256 of the previous
  entry, so any modification of the log is detected by :meth:`verify`.
* Every mutating call writes its audit entry in the *same transaction* as the
  change, so a change can never exist without its log entry.
"""

from __future__ import annotations

import getpass
import hashlib
import json
import socket
import sqlite3
from contextlib import contextmanager
from dataclasses import dataclass, field
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Iterator

import numpy as np

from spectro import __version__
from spectro.core.io import ParseOptions, file_sha256, load_file
from spectro.core.operations import REGISTRY, apply_pipeline, describe_step
from spectro.core.spectrum import Spectrum

SCHEMA_VERSION = 1
GENESIS = "0" * 64
ROLES = ("", "standard", "calibration", "validation", "mixture", "sample", "blank",
         "divisor", "reference", "unknown")

_SCHEMA = """
CREATE TABLE IF NOT EXISTS meta (key TEXT PRIMARY KEY, value TEXT);
CREATE TABLE IF NOT EXISTS compounds (
    id INTEGER PRIMARY KEY, name TEXT UNIQUE NOT NULL, unit TEXT DEFAULT 'µg/mL',
    mw REAL, notes TEXT DEFAULT '');
CREATE TABLE IF NOT EXISTS trials (
    id INTEGER PRIMARY KEY, name TEXT NOT NULL, description TEXT DEFAULT '',
    status TEXT DEFAULT 'open', metadata TEXT DEFAULT '{}', created_utc TEXT NOT NULL,
    archived INTEGER DEFAULT 0);
CREATE TABLE IF NOT EXISTS spectra (
    id INTEGER PRIMARY KEY, trial_id INTEGER REFERENCES trials(id),
    name TEXT NOT NULL, kind TEXT NOT NULL CHECK (kind IN ('raw','derived')),
    role TEXT DEFAULT '', parent_ids TEXT DEFAULT '[]', pipeline TEXT DEFAULT '[]',
    xdata BLOB NOT NULL, ydata BLOB NOT NULL, npoints INTEGER NOT NULL,
    metadata TEXT DEFAULT '{}', concentrations TEXT DEFAULT '{}', notes TEXT DEFAULT '',
    source_file TEXT DEFAULT '', source_sha256 TEXT DEFAULT '',
    data_sha256 TEXT NOT NULL, created_utc TEXT NOT NULL, archived INTEGER DEFAULT 0);
CREATE TABLE IF NOT EXISTS methods (
    id INTEGER PRIMARY KEY, trial_id INTEGER, name TEXT NOT NULL, type TEXT NOT NULL,
    definition TEXT NOT NULL, created_utc TEXT NOT NULL, archived INTEGER DEFAULT 0);
CREATE TABLE IF NOT EXISTS results (
    id INTEGER PRIMARY KEY, trial_id INTEGER, method_id INTEGER, name TEXT NOT NULL,
    kind TEXT NOT NULL, data TEXT NOT NULL, created_utc TEXT NOT NULL,
    archived INTEGER DEFAULT 0);
CREATE TABLE IF NOT EXISTS audit (
    seq INTEGER PRIMARY KEY AUTOINCREMENT, ts_utc TEXT NOT NULL, user TEXT NOT NULL,
    host TEXT NOT NULL, app_version TEXT NOT NULL, action TEXT NOT NULL,
    entity TEXT NOT NULL, entity_id INTEGER, summary TEXT NOT NULL,
    details TEXT NOT NULL, reason TEXT DEFAULT '', prev_hash TEXT NOT NULL,
    hash TEXT NOT NULL);
CREATE INDEX IF NOT EXISTS audit_entity ON audit(entity, entity_id);
CREATE INDEX IF NOT EXISTS spectra_trial ON spectra(trial_id);

CREATE TRIGGER IF NOT EXISTS audit_no_update BEFORE UPDATE ON audit
BEGIN SELECT RAISE(ABORT, 'audit trail is append-only'); END;
CREATE TRIGGER IF NOT EXISTS audit_no_delete BEFORE DELETE ON audit
BEGIN SELECT RAISE(ABORT, 'audit trail is append-only'); END;
CREATE TRIGGER IF NOT EXISTS spectra_no_delete BEFORE DELETE ON spectra
BEGIN SELECT RAISE(ABORT, 'spectra cannot be deleted; archive them instead'); END;
CREATE TRIGGER IF NOT EXISTS spectra_immutable_data
BEFORE UPDATE OF xdata, ydata, npoints, data_sha256, kind, parent_ids, pipeline,
                 source_file, source_sha256, created_utc ON spectra
BEGIN SELECT RAISE(ABORT, 'spectral data is immutable'); END;
CREATE TRIGGER IF NOT EXISTS methods_no_delete BEFORE DELETE ON methods
BEGIN SELECT RAISE(ABORT, 'methods cannot be deleted; archive them instead'); END;
CREATE TRIGGER IF NOT EXISTS results_no_delete BEFORE DELETE ON results
BEGIN SELECT RAISE(ABORT, 'results cannot be deleted; archive them instead'); END;
CREATE TRIGGER IF NOT EXISTS results_immutable BEFORE UPDATE OF data ON results
BEGIN SELECT RAISE(ABORT, 'results are immutable'); END;
"""


def utc_now() -> str:
    return datetime.now(timezone.utc).isoformat(timespec="milliseconds")


def _json(obj: Any) -> str:
    return json.dumps(obj, sort_keys=True, ensure_ascii=False, default=_default)


def _default(o: Any):
    if isinstance(o, np.ndarray):
        return o.tolist()
    if isinstance(o, (np.floating, np.integer)):
        return o.item()
    if isinstance(o, Spectrum):
        return {"spectrum": o.name}
    raise TypeError(f"not serialisable: {type(o)}")


def _current_user() -> str:
    try:
        return getpass.getuser()
    except Exception:  # pragma: no cover
        return "unknown"


@dataclass
class SpectrumRecord:
    id: int
    trial_id: int | None
    name: str
    kind: str
    role: str
    parent_ids: list[int]
    pipeline: list[dict]
    npoints: int
    metadata: dict
    concentrations: dict[str, float]
    notes: str
    source_file: str
    source_sha256: str
    data_sha256: str
    created_utc: str
    archived: bool


@dataclass
class AuditEntry:
    seq: int
    ts_utc: str
    user: str
    host: str
    app_version: str
    action: str
    entity: str
    entity_id: int | None
    summary: str
    details: dict = field(default_factory=dict)
    reason: str = ""
    prev_hash: str = ""
    hash: str = ""


def entry_hash(prev_hash: str, ts: str, user: str, host: str, version: str, action: str,
               entity: str, entity_id: int | None, summary: str, details: str,
               reason: str) -> str:
    payload = "\x1f".join([prev_hash, ts, user, host, version, action, entity,
                           "" if entity_id is None else str(entity_id), summary,
                           details, reason])
    return hashlib.sha256(payload.encode("utf-8")).hexdigest()


class Project:
    """An open ``.spectro`` project file."""

    def __init__(self, path: str | Path, conn: sqlite3.Connection):
        self.path = Path(path)
        self.conn = conn
        self.user = _current_user()
        self.host = socket.gethostname()
        self._listeners: list = []

    # ------------------------------------------------------------ lifecycle
    @classmethod
    def create(cls, path: str | Path, name: str = "", description: str = "") -> "Project":
        path = Path(path)
        if path.exists():
            raise FileExistsError(f"{path} already exists")
        conn = sqlite3.connect(str(path), isolation_level=None)
        conn.executescript(_SCHEMA)
        proj = cls(path, conn)
        with proj._tx():
            conn.executemany("INSERT INTO meta VALUES (?, ?)", [
                ("schema_version", str(SCHEMA_VERSION)), ("name", name or path.stem),
                ("description", description), ("created_utc", utc_now()),
                ("created_by", proj.user)])
            proj._audit("CREATE", "project", None, f"Created project '{name or path.stem}'",
                        {"path": str(path)})
        return proj

    @classmethod
    def open(cls, path: str | Path) -> "Project":
        path = Path(path)
        if not path.exists():
            raise FileNotFoundError(path)
        conn = sqlite3.connect(str(path), isolation_level=None)
        try:
            conn.executescript(_SCHEMA)
            ver = conn.execute("SELECT value FROM meta WHERE key='schema_version'").fetchone()
        except sqlite3.DatabaseError as exc:
            raise ValueError(f"{path.name} is not a Spectro project") from exc
        if ver is None:
            raise ValueError(f"{path.name} is not a Spectro project")
        proj = cls(path, conn)
        with proj._tx():
            proj._audit("OPEN", "project", None, "Opened project", {"path": str(path)})
        return proj

    def close(self) -> None:
        if self.conn:
            with self._tx():
                self._audit("CLOSE", "project", None, "Closed project", {})
            self.conn.close()
            self.conn = None  # type: ignore[assignment]

    def subscribe(self, callback) -> None:
        """Register ``callback(entity, entity_id)`` called after each change."""
        self._listeners.append(callback)

    def _notify(self, entity: str, entity_id: int | None) -> None:
        for cb in list(self._listeners):
            cb(entity, entity_id)

    @contextmanager
    def _tx(self) -> Iterator[sqlite3.Connection]:
        try:
            self.conn.execute("BEGIN IMMEDIATE")
            yield self.conn
            self.conn.execute("COMMIT")
        except BaseException:
            self.conn.execute("ROLLBACK")
            raise

    # ------------------------------------------------------------ audit
    def _audit(self, action: str, entity: str, entity_id: int | None, summary: str,
               details: dict, reason: str = "") -> int:
        row = self.conn.execute("SELECT hash FROM audit ORDER BY seq DESC LIMIT 1").fetchone()
        prev = row[0] if row else GENESIS
        ts = utc_now()
        det = _json(details)
        h = entry_hash(prev, ts, self.user, self.host, __version__, action, entity,
                       entity_id, summary, det, reason)
        cur = self.conn.execute(
            "INSERT INTO audit (ts_utc, user, host, app_version, action, entity, entity_id,"
            " summary, details, reason, prev_hash, hash) VALUES (?,?,?,?,?,?,?,?,?,?,?,?)",
            (ts, self.user, self.host, __version__, action, entity, entity_id, summary,
             det, reason, prev, h))
        return int(cur.lastrowid)

    def log(self, action: str, summary: str, details: dict | None = None,
            entity: str = "project", entity_id: int | None = None, reason: str = "") -> None:
        """Record an operation that does not change stored data (export,
        calculation, report…)."""
        with self._tx():
            self._audit(action, entity, entity_id, summary, details or {}, reason)
        self._notify("audit", None)

    def audit_entries(self, entity: str | None = None, entity_id: int | None = None,
                      action: str | None = None, text: str | None = None,
                      limit: int | None = None) -> list[AuditEntry]:
        q = "SELECT * FROM audit WHERE 1=1"
        args: list[Any] = []
        if entity:
            q += " AND entity=?"
            args.append(entity)
        if entity_id is not None:
            q += " AND entity_id=?"
            args.append(entity_id)
        if action:
            q += " AND action=?"
            args.append(action)
        if text:
            q += " AND (summary LIKE ? OR details LIKE ? OR reason LIKE ?)"
            args += [f"%{text}%"] * 3
        q += " ORDER BY seq DESC"
        if limit:
            q += f" LIMIT {int(limit)}"
        out = []
        for r in self.conn.execute(q, args):
            e = AuditEntry(*r)
            e.details = json.loads(e.details)
            out.append(e)
        return out

    def spectrum_history(self, spectrum_id: int) -> list[AuditEntry]:
        """Audit entries touching one spectrum (as target, input or output)."""
        rows = self.conn.execute(
            "SELECT * FROM audit WHERE (entity='spectrum' AND entity_id=?) "
            "OR details LIKE ? ORDER BY seq", (spectrum_id, f"%{spectrum_id}%"))
        out = []
        for r in rows:
            e = AuditEntry(*r)
            e.details = json.loads(e.details)
            ids = set(e.details.get("inputs", [])) | set(e.details.get("outputs", []))
            ids |= set(e.details.get("references", []))
            if (e.entity == "spectrum" and e.entity_id == spectrum_id) or spectrum_id in ids:
                out.append(e)
        return out

    def verify(self) -> dict:
        """Check the audit hash chain and every spectrum's data hash."""
        problems: list[str] = []
        prev = GENESIS
        n = 0
        for r in self.conn.execute(
                "SELECT seq, ts_utc, user, host, app_version, action, entity, entity_id,"
                " summary, details, reason, prev_hash, hash FROM audit ORDER BY seq"):
            n += 1
            (seq, ts, user, host, ver, action, entity, eid, summary, det, reason,
             prev_hash, h) = r
            if prev_hash != prev:
                problems.append(f"audit #{seq}: chain broken (previous hash mismatch)")
            if entry_hash(prev_hash, ts, user, host, ver, action, entity, eid, summary,
                          det, reason) != h:
                problems.append(f"audit #{seq}: entry content was modified")
            prev = h
        ns = 0
        for sid, x, y, dh in self.conn.execute("SELECT id, xdata, ydata, data_sha256 FROM spectra"):
            ns += 1
            s = Spectrum(np.frombuffer(x, "<f8"), np.frombuffer(y, "<f8"))
            if s.data_hash() != dh:
                problems.append(f"spectrum #{sid}: data does not match its hash")
        return {"ok": not problems, "problems": problems, "audit_entries": n, "spectra": ns}

    # ------------------------------------------------------------ project meta
    def meta(self, key: str, default: str = "") -> str:
        row = self.conn.execute("SELECT value FROM meta WHERE key=?", (key,)).fetchone()
        return row[0] if row else default

    def set_meta(self, key: str, value: str, reason: str = "") -> None:
        old = self.meta(key)
        if old == value:
            return
        with self._tx():
            self.conn.execute("INSERT OR REPLACE INTO meta VALUES (?, ?)", (key, value))
            self._audit("EDIT", "project", None, f"Changed project {key}",
                        {"field": key, "before": old, "after": value}, reason)
        self._notify("project", None)

    # ------------------------------------------------------------ compounds
    def compounds(self) -> list[dict]:
        return [dict(zip(("id", "name", "unit", "mw", "notes"), r))
                for r in self.conn.execute("SELECT id, name, unit, mw, notes FROM compounds ORDER BY id")]

    def compound_names(self) -> list[str]:
        return [c["name"] for c in self.compounds()]

    def add_compound(self, name: str, unit: str = "µg/mL", mw: float | None = None,
                     notes: str = "") -> int:
        name = name.strip()
        if not name:
            raise ValueError("compound name is empty")
        with self._tx():
            cid = self.conn.execute("INSERT INTO compounds (name, unit, mw, notes) VALUES (?,?,?,?)",
                                    (name, unit, mw, notes)).lastrowid
            self._audit("CREATE", "compound", cid, f"Added compound '{name}'",
                        {"name": name, "unit": unit, "mw": mw})
        self._notify("compound", cid)
        return int(cid)

    def update_compound(self, cid: int, reason: str = "", **fields: Any) -> None:
        cur = dict(zip(("name", "unit", "mw", "notes"), self.conn.execute(
            "SELECT name, unit, mw, notes FROM compounds WHERE id=?", (cid,)).fetchone()))
        changes = {k: v for k, v in fields.items() if k in cur and cur[k] != v}
        if not changes:
            return
        with self._tx():
            for k, v in changes.items():
                self.conn.execute(f"UPDATE compounds SET {k}=? WHERE id=?", (v, cid))
            self._audit("EDIT", "compound", cid, f"Edited compound '{cur['name']}'",
                        {"before": {k: cur[k] for k in changes}, "after": changes}, reason)
        self._notify("compound", cid)

    # ------------------------------------------------------------ trials
    def trials(self, include_archived: bool = False) -> list[dict]:
        q = "SELECT id, name, description, status, metadata, created_utc, archived FROM trials"
        if not include_archived:
            q += " WHERE archived=0"
        return [{"id": r[0], "name": r[1], "description": r[2], "status": r[3],
                 "metadata": json.loads(r[4]), "created_utc": r[5], "archived": bool(r[6])}
                for r in self.conn.execute(q + " ORDER BY id")]

    def trial(self, tid: int) -> dict:
        return next(t for t in self.trials(True) if t["id"] == tid)

    def add_trial(self, name: str, description: str = "", metadata: dict | None = None) -> int:
        with self._tx():
            tid = self.conn.execute(
                "INSERT INTO trials (name, description, metadata, created_utc) VALUES (?,?,?,?)",
                (name, description, _json(metadata or {}), utc_now())).lastrowid
            self._audit("CREATE", "trial", tid, f"Created trial '{name}'",
                        {"name": name, "description": description, "metadata": metadata or {}})
        self._notify("trial", tid)
        return int(tid)

    def update_trial(self, tid: int, reason: str = "", **fields: Any) -> None:
        cur = self.trial(tid)
        changes = {k: v for k, v in fields.items()
                   if k in ("name", "description", "status", "metadata") and cur[k] != v}
        if not changes:
            return
        with self._tx():
            for k, v in changes.items():
                self.conn.execute(f"UPDATE trials SET {k}=? WHERE id=?",
                                  (_json(v) if k == "metadata" else v, tid))
            self._audit("EDIT", "trial", tid, f"Edited trial '{cur['name']}'",
                        {"before": {k: cur[k] for k in changes}, "after": changes}, reason)
        self._notify("trial", tid)

    def archive_trial(self, tid: int, reason: str, archived: bool = True) -> None:
        with self._tx():
            self.conn.execute("UPDATE trials SET archived=? WHERE id=?", (int(archived), tid))
            self._audit("ARCHIVE" if archived else "RESTORE", "trial", tid,
                        f"{'Archived' if archived else 'Restored'} trial #{tid}", {}, reason)
        self._notify("trial", tid)

    # ------------------------------------------------------------ spectra
    def _row_to_record(self, r) -> SpectrumRecord:
        return SpectrumRecord(
            id=r[0], trial_id=r[1], name=r[2], kind=r[3], role=r[4],
            parent_ids=json.loads(r[5]), pipeline=json.loads(r[6]), npoints=r[7],
            metadata=json.loads(r[8]), concentrations=json.loads(r[9]), notes=r[10],
            source_file=r[11], source_sha256=r[12], data_sha256=r[13], created_utc=r[14],
            archived=bool(r[15]))

    _REC_COLS = ("id, trial_id, name, kind, role, parent_ids, pipeline, npoints, metadata, "
                 "concentrations, notes, source_file, source_sha256, data_sha256, "
                 "created_utc, archived")

    def records(self, trial_id: int | None = None, include_archived: bool = False,
                ids: list[int] | None = None) -> list[SpectrumRecord]:
        q = f"SELECT {self._REC_COLS} FROM spectra WHERE 1=1"
        args: list[Any] = []
        if trial_id is not None:
            q += " AND trial_id=?"
            args.append(trial_id)
        if not include_archived:
            q += " AND archived=0"
        if ids is not None:
            q += f" AND id IN ({','.join('?' * len(ids))})"
            args += list(ids)
        return [self._row_to_record(r) for r in self.conn.execute(q + " ORDER BY id", args)]

    def record(self, sid: int) -> SpectrumRecord:
        r = self.conn.execute(f"SELECT {self._REC_COLS} FROM spectra WHERE id=?", (sid,)).fetchone()
        if r is None:
            raise KeyError(f"spectrum #{sid} not found")
        return self._row_to_record(r)

    def spectrum(self, sid: int) -> Spectrum:
        r = self.conn.execute("SELECT name, xdata, ydata, metadata, concentrations, role "
                              "FROM spectra WHERE id=?", (int(sid),)).fetchone()
        if r is None:
            raise KeyError(f"spectrum #{sid} not found")
        meta = json.loads(r[3])
        meta.update({"id": int(sid), "role": r[5]})
        return Spectrum(np.frombuffer(r[1], "<f8").copy(), np.frombuffer(r[2], "<f8").copy(),
                        name=r[0], metadata=meta, concentrations=json.loads(r[4]))

    def spectra(self, ids: list[int]) -> list[Spectrum]:
        return [self.spectrum(i) for i in ids]

    def resolver(self):
        return lambda ref: self.spectrum(int(ref))

    def _insert_spectrum(self, s: Spectrum, trial_id: int | None, kind: str,
                         parents: list[int], pipeline: list[dict], role: str = "",
                         source_file: str = "", source_sha: str = "") -> int:
        meta = {k: v for k, v in s.metadata.items() if k not in ("id", "role")}
        return int(self.conn.execute(
            "INSERT INTO spectra (trial_id, name, kind, role, parent_ids, pipeline, xdata, ydata,"
            " npoints, metadata, concentrations, source_file, source_sha256, data_sha256,"
            " created_utc) VALUES (?,?,?,?,?,?,?,?,?,?,?,?,?,?,?)",
            (trial_id, s.name, kind, role, _json(parents), _json(pipeline),
             s.wavelengths.astype("<f8").tobytes(), s.values.astype("<f8").tobytes(),
             int(s.wavelengths.size), _json(meta), _json(s.concentrations),
             source_file, source_sha, s.data_hash(), utc_now())).lastrowid)

    def import_file(self, path: str | Path, trial_id: int | None,
                    options: ParseOptions | None = None, role: str = "",
                    only: list[int] | None = None,
                    concentrations: dict[int, dict[str, float]] | None = None,
                    names: dict[int, str] | None = None,
                    roles: dict[int, str] | None = None) -> list[int]:
        """Import spectra from any supported file (raw, immutable). ``roles``
        overrides ``role`` per spectrum (e.g. from the Excel sheet name)."""
        path = Path(path)
        result = load_file(path, options)
        sha = file_sha256(path)
        ids: list[int] = []
        with self._tx():
            for i, s in enumerate(result.spectra):
                if only is not None and i not in only:
                    continue
                if names and i in names:
                    s.name = names[i]
                if concentrations and i in concentrations:
                    s.concentrations = dict(concentrations[i])
                s.metadata.update({k: v for k, v in result.metadata.items()
                                   if not k.startswith("_")})
                ids.append(self._insert_spectrum(s, trial_id, "raw", [], [],
                                                 (roles or {}).get(i, role), str(path), sha))
            self._audit("IMPORT", "file", None,
                        f"Imported {len(ids)} spectra from {path.name}",
                        {"file": str(path), "file_sha256": sha, "layout": result.layout,
                         "format": result.metadata, "outputs": ids, "trial_id": trial_id,
                         "data_sha256": {str(i): self.record(i).data_sha256 for i in ids},
                         "warnings": result.warnings})
        self._notify("spectrum", None)
        return ids

    def add_spectrum(self, s: Spectrum, trial_id: int | None, role: str = "",
                     note: str = "manual entry") -> int:
        """Add a raw spectrum created in memory (e.g. pasted data)."""
        with self._tx():
            sid = self._insert_spectrum(s, trial_id, "raw", [], [], role)
            self._audit("CREATE", "spectrum", sid, f"Added spectrum '{s.name}' ({note})",
                        {"outputs": [sid], "data_sha256": s.data_hash()})
        self._notify("spectrum", sid)
        return sid

    def process(self, ids: list[int], steps: list[dict], suffix: str | None = None,
                trial_id: int | None = None, role: str | None = None) -> list[int]:
        """Apply a pipeline to spectra, storing each result as a derived spectrum."""
        if not steps:
            raise ValueError("no processing steps")
        for st in steps:
            if st["op"] not in REGISTRY:
                raise KeyError(f"unknown operation {st['op']}")
        refs = sorted({int(v) for st in steps for p in REGISTRY[st["op"]].params
                       if p.kind == "spectrum"
                       for v in [st.get("params", {}).get(p.name)] if v is not None})
        resolve = self.resolver()
        name_of = lambda ref: self.record(int(ref)).name  # noqa: E731
        label = suffix or " → ".join(REGISTRY[st["op"]].label.split(" (")[0] for st in steps)
        out_ids = []
        with self._tx():
            for sid in ids:
                rec = self.record(sid)
                src = self.spectrum(sid)
                res = apply_pipeline(src, steps, resolve)
                res.name = f"{rec.name} | {label}"
                res.metadata = {k: v for k, v in rec.metadata.items()}
                res.metadata["processing"] = [describe_step(st, name_of) for st in steps]
                new_id = self._insert_spectrum(
                    res, rec.trial_id if trial_id is None else trial_id, "derived",
                    [sid] + refs, steps, rec.role if role is None else role)
                out_ids.append(new_id)
            self._audit("PROCESS", "spectrum", out_ids[0] if len(out_ids) == 1 else None,
                        f"Processed {len(ids)} spectra: {label}",
                        {"inputs": list(ids), "references": refs, "outputs": out_ids,
                         "steps": steps,
                         "description": [describe_step(st, name_of) for st in steps],
                         "input_sha256": {str(i): self.record(i).data_sha256 for i in ids},
                         "output_sha256": {str(i): self.record(i).data_sha256 for i in out_ids}})
        self._notify("spectrum", None)
        return out_ids

    def replay(self, sid: int) -> dict:
        """Recompute a derived spectrum from its parent and compare hashes."""
        rec = self.record(sid)
        if rec.kind != "derived":
            return {"ok": True, "message": "raw spectrum (nothing to replay)"}
        res = apply_pipeline(self.spectrum(rec.parent_ids[0]), rec.pipeline, self.resolver())
        stored = self.spectrum(sid)
        same = res.same_grid(stored) and bool(np.allclose(res.values, stored.values,
                                                          rtol=1e-10, atol=1e-12))
        self.log("VERIFY", f"Replayed processing of spectrum #{sid}: "
                           f"{'reproduced' if same else 'MISMATCH'}",
                 {"inputs": [rec.parent_ids[0]], "outputs": [sid], "reproduced": same},
                 entity="spectrum", entity_id=sid)
        return {"ok": same, "message": "result reproduced exactly" if same
                else "recomputed result differs from stored data"}

    def lineage(self, sid: int) -> list[SpectrumRecord]:
        """Chain of records from the raw ancestor to ``sid``."""
        chain = [self.record(sid)]
        while chain[-1].kind == "derived" and chain[-1].parent_ids:
            chain.append(self.record(chain[-1].parent_ids[0]))
        return list(reversed(chain))

    def update_spectrum(self, sid: int, reason: str = "", **fields: Any) -> None:
        """Edit descriptive fields (name, role, concentrations, metadata, notes,
        trial_id). Numeric data cannot be changed."""
        rec = self.record(sid)
        allowed = ("name", "role", "concentrations", "metadata", "notes", "trial_id")
        changes = {k: v for k, v in fields.items() if k in allowed and getattr(rec, k) != v}
        if not changes:
            return
        with self._tx():
            for k, v in changes.items():
                val = _json(v) if k in ("concentrations", "metadata") else v
                self.conn.execute(f"UPDATE spectra SET {k}=? WHERE id=?", (val, sid))
            self._audit("EDIT", "spectrum", sid, f"Edited spectrum '{rec.name}': "
                        + ", ".join(changes), {"before": {k: getattr(rec, k) for k in changes},
                                               "after": changes}, reason)
        self._notify("spectrum", sid)

    def set_concentrations_bulk(self, values: dict[int, dict[str, float]], reason: str = "") -> None:
        with self._tx():
            before = {}
            for sid, conc in values.items():
                rec = self.record(sid)
                if rec.concentrations != conc:
                    before[sid] = rec.concentrations
                    self.conn.execute("UPDATE spectra SET concentrations=? WHERE id=?",
                                      (_json(conc), sid))
            if before:
                self._audit("EDIT", "spectrum", None,
                            f"Edited concentrations of {len(before)} spectra",
                            {"before": {str(k): v for k, v in before.items()},
                             "after": {str(k): values[k] for k in before},
                             "inputs": list(before)}, reason)
        self._notify("spectrum", None)

    def archive_spectra(self, ids: list[int], reason: str, archived: bool = True) -> None:
        with self._tx():
            for sid in ids:
                self.conn.execute("UPDATE spectra SET archived=? WHERE id=?", (int(archived), sid))
            self._audit("ARCHIVE" if archived else "RESTORE", "spectrum",
                        ids[0] if len(ids) == 1 else None,
                        f"{'Archived' if archived else 'Restored'} {len(ids)} spectra",
                        {"inputs": list(ids)}, reason)
        self._notify("spectrum", None)

    # ------------------------------------------------------------ methods & results
    def save_method(self, name: str, definition: dict, trial_id: int | None = None) -> int:
        with self._tx():
            mid = self.conn.execute(
                "INSERT INTO methods (trial_id, name, type, definition, created_utc) VALUES (?,?,?,?,?)",
                (trial_id, name, definition.get("type", ""), _json(definition), utc_now())).lastrowid
            self._audit("CREATE", "method", mid, f"Saved method '{name}'",
                        {"definition": definition, "trial_id": trial_id})
        self._notify("method", mid)
        return int(mid)

    def methods(self, trial_id: int | None = None) -> list[dict]:
        q = "SELECT id, trial_id, name, type, definition, created_utc FROM methods WHERE archived=0"
        args: list[Any] = []
        if trial_id is not None:
            q += " AND (trial_id=? OR trial_id IS NULL)"
            args.append(trial_id)
        return [{"id": r[0], "trial_id": r[1], "name": r[2], "type": r[3],
                 "definition": json.loads(r[4]), "created_utc": r[5]}
                for r in self.conn.execute(q + " ORDER BY id", args)]

    def archive_method(self, mid: int, reason: str) -> None:
        with self._tx():
            self.conn.execute("UPDATE methods SET archived=1 WHERE id=?", (mid,))
            self._audit("ARCHIVE", "method", mid, f"Archived method #{mid}", {}, reason)
        self._notify("method", mid)

    def save_result(self, name: str, kind: str, data: dict, trial_id: int | None = None,
                    method_id: int | None = None, inputs: list[int] | None = None) -> int:
        with self._tx():
            rid = self.conn.execute(
                "INSERT INTO results (trial_id, method_id, name, kind, data, created_utc)"
                " VALUES (?,?,?,?,?,?)",
                (trial_id, method_id, name, kind, _json(data), utc_now())).lastrowid
            self._audit("CALCULATE", "result", rid, f"Saved result '{name}' ({kind})",
                        {"inputs": inputs or [], "method_id": method_id, "trial_id": trial_id,
                         "data": data})
        self._notify("result", rid)
        return int(rid)

    def results(self, trial_id: int | None = None) -> list[dict]:
        q = "SELECT id, trial_id, method_id, name, kind, data, created_utc FROM results WHERE archived=0"
        args: list[Any] = []
        if trial_id is not None:
            q += " AND trial_id=?"
            args.append(trial_id)
        return [{"id": r[0], "trial_id": r[1], "method_id": r[2], "name": r[3], "kind": r[4],
                 "data": json.loads(r[5]), "created_utc": r[6]}
                for r in self.conn.execute(q + " ORDER BY id", args)]
