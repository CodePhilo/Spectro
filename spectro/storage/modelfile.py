"""Portable model files (.spmodel): a calibrated method plus any spectra it
references (divisors, references), with checksums."""

from __future__ import annotations

import hashlib
import json

from spectro import __version__
from spectro.core.operations import REGISTRY
from spectro.core.spectrum import Spectrum

FORMAT = "spectro-model"


def spectrum_refs(d: dict) -> set[int]:
    """Spectrum ids referenced by a method definition."""
    out: set[int] = set()
    for st in d.get("steps", []):
        op = REGISTRY.get(st.get("op"))
        for p in (op.params if op else []):
            v = st.get("params", {}).get(p.name)
            if p.kind == "spectrum" and v is not None:
                out.add(int(v))
    for sig in [d.get("measurement") or {}] + list(d.get("signals") or []):
        v = sig.get("params", {}).get("reference")
        if v is not None:
            out.add(int(v))
    return out


def remap_refs(d: dict, mapping: dict[int, int]) -> None:
    for st in d.get("steps", []):
        op = REGISTRY.get(st.get("op"))
        for p in (op.params if op else []):
            v = st.get("params", {}).get(p.name)
            if p.kind == "spectrum" and v is not None and int(v) in mapping:
                st["params"][p.name] = mapping[int(v)]
    for sig in [d.get("measurement") or {}] + list(d.get("signals") or []):
        v = sig.get("params", {}).get("reference")
        if v is not None and int(v) in mapping:
            sig["params"]["reference"] = mapping[int(v)]


def _checksum(d: dict) -> str:
    return hashlib.sha256(json.dumps(d, sort_keys=True, ensure_ascii=False).encode()).hexdigest()


def is_calibrated(d: dict) -> bool:
    t = d.get("type")
    if t == "spectral":
        return bool(d.get("fitted"))
    if t == "univariate":
        return d.get("regression") is not None
    if t == "equations":
        return d.get("K") is not None
    return False


def export_model(project, method: dict) -> dict:
    """Model-file document for a saved method (fits spectral models if needed)."""
    from spectro.core.multicomponent import SpectralModel
    from spectro.storage.project import utc_now

    d = json.loads(json.dumps(method["definition"]))
    if d.get("type") == "spectral" and not d.get("fitted"):
        m = SpectralModel.from_dict(d)
        m.fit(project.spectra(d["calibration_ids"]), project.resolver())
        d.update(m.to_dict(include_fit=True))
    if not is_calibrated(d):
        raise ValueError("the method is not calibrated")
    embedded = {}
    for sid in sorted(spectrum_refs(d)):
        sp = project.spectrum(sid)
        embedded[str(sid)] = {"name": sp.name, "wavelengths": sp.wavelengths.tolist(),
                              "values": sp.values.tolist(), "data_sha256": sp.data_hash()}
    d.pop("calibration_ids", None)
    return {"format": FORMAT, "format_version": 1, "app_version": __version__,
            "name": method["name"], "exported_utc": utc_now(),
            "source_project": project.meta("name"), "source_method_id": method["id"],
            "definition": d, "embedded_spectra": embedded, "sha256": _checksum(d)}


def import_model(project, doc: dict, trial_id: int | None, source: str = "") -> int:
    """Validate a model-file document and save it as a method in ``project``."""
    if doc.get("format") != FORMAT:
        raise ValueError("not a Spectro model file")
    d = json.loads(json.dumps(doc["definition"]))
    if _checksum(d) != doc.get("sha256"):
        raise ValueError("the model file was modified after export (checksum mismatch)")
    if not is_calibrated(d):
        raise ValueError("the model file does not contain a calibrated model")
    spectra = {}
    for old, sp in (doc.get("embedded_spectra") or {}).items():
        spec = Spectrum(sp["wavelengths"], sp["values"], name=f"{sp['name']} (from model)")
        if spec.data_hash() != sp.get("data_sha256"):
            raise ValueError(f"embedded spectrum '{sp['name']}' failed its checksum")
        spectra[int(old)] = spec
    mapping = {old: project.add_spectrum(spec, trial_id, role="divisor",
                                         note=f"embedded in model {doc.get('name', '')}")
               for old, spec in spectra.items()}
    remap_refs(d, mapping)
    d["imported_from"] = {"file": source, "sha256": doc["sha256"],
                          "source_project": doc.get("source_project"),
                          "exported_utc": doc.get("exported_utc")}
    return project.save_method(f"{doc.get('name', 'Imported model')} (imported)", d, trial_id)
