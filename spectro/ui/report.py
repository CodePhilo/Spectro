"""Trial report: spectra, processing lineage, methods, results and audit trail."""

from __future__ import annotations

import base64
import html
import json
import tempfile
from pathlib import Path

from PySide6.QtCore import QMarginsF
from PySide6.QtGui import QPageLayout, QPageSize, QPdfWriter, QTextDocument
from PySide6.QtWidgets import QFileDialog, QInputDialog

from spectro import __version__
from spectro.core.operations import describe_step
from spectro.storage.project import utc_now
from spectro.ui.widgets import SpectrumPlot, error, fmt

CSS = """
body { font-family: 'Segoe UI', Arial, sans-serif; font-size: 9pt; color: #0b0b0b; }
h1 { font-size: 16pt; margin-bottom: 2px; } h2 { font-size: 12pt; margin-top: 16px;
border-bottom: 1px solid #c3c2b7; } table { border-collapse: collapse; width: 100%; }
th, td { border: 1px solid #c3c2b7; padding: 3px 5px; text-align: left; vertical-align: top; }
th { background: #f0efec; } .muted { color: #52514e; } td.num { text-align: right; }
"""


def _plot_png(spectra) -> str:
    plot = SpectrumPlot()
    plot.resize(900, 420)
    plot.plot_spectra(spectra)
    from pyqtgraph.exporters import ImageExporter

    exp = ImageExporter(plot.plotItem)
    exp.parameters()["width"] = 1400
    with tempfile.TemporaryDirectory() as d:
        f = Path(d) / "plot.png"
        exp.export(str(f))
        data = base64.b64encode(f.read_bytes()).decode()
    plot.deleteLater()
    return data


def build_html(project, trial_id, include_plot: bool = True) -> str:
    e = html.escape
    trial = project.trial(trial_id) if trial_id is not None else {"name": "All spectra",
                                                                  "description": "", "status": ""}
    recs = project.records(trial_id=trial_id)
    name_of = lambda ref: project.record(int(ref)).name  # noqa: E731
    parts = [f"<html><head><style>{CSS}</style></head><body>",
             f"<h1>{e(trial['name'])}</h1>",
             f"<p class='muted'>Project: {e(project.meta('name'))} ({e(str(project.path))})<br>"
             f"Generated {utc_now()} UTC by {e(project.user)}@{e(project.host)} — "
             f"Spectro {__version__}</p>"]
    if trial.get("description"):
        parts.append(f"<p>{e(trial['description'])}</p>")
    integ = project.verify()
    parts.append(f"<p><b>Data integrity:</b> {'PASSED' if integ['ok'] else 'FAILED'} — "
                 f"{integ['audit_entries']} audit entries, {integ['spectra']} spectra checked.</p>")
    comps = project.compounds()
    if comps:
        parts.append("<h2>Compounds</h2><table><tr><th>Name</th><th>Unit</th><th>MW</th></tr>")
        parts += [f"<tr><td>{e(c['name'])}</td><td>{e(c['unit'] or '')}</td>"
                  f"<td>{'' if c['mw'] is None else c['mw']}</td></tr>" for c in comps]
        parts.append("</table>")
    raw = [r for r in recs if r.kind == "raw"]
    if include_plot and raw:
        png = _plot_png(project.spectra([r.id for r in raw[:40]]))
        parts.append(f"<h2>Raw spectra</h2><img src='data:image/png;base64,{png}' width='640'>")
    parts.append("<h2>Spectra</h2><table><tr><th>#</th><th>Name</th><th>Kind</th><th>Role</th>"
                 "<th>Concentrations</th><th>Source / processing</th><th>Data SHA-256</th></tr>")
    for r in recs:
        conc = ", ".join(f"{k} = {v:g}" for k, v in r.concentrations.items())
        if r.kind == "raw":
            src = e(Path(r.source_file).name) if r.source_file else "manual"
        else:
            src = f"from #{r.parent_ids[0]}: " + "; ".join(e(describe_step(s, name_of))
                                                         for s in r.pipeline)
        parts.append(f"<tr><td>{r.id}</td><td>{e(r.name)}</td><td>{r.kind}</td>"
                     f"<td>{e(r.role)}</td><td>{e(conc)}</td><td>{src}</td>"
                     f"<td class='muted'>{r.data_sha256[:16]}…</td></tr>")
    parts.append("</table>")
    methods = project.methods(trial_id)
    if methods:
        parts.append("<h2>Methods</h2><table><tr><th>#</th><th>Name</th><th>Type</th>"
                     "<th>Definition</th></tr>")
        for m in methods:
            d = dict(m["definition"])
            reg = d.pop("regression", None)
            steps = "; ".join(describe_step(s, name_of) for s in d.get("steps", []))
            info = []
            if steps:
                info.append(f"Steps: {e(steps)}")
            if d.get("measurement"):
                info.append(f"Measurement: {e(json.dumps(d['measurement']))}")
            if reg:
                info.append(f"y = {reg['slope']:.5g}·C + {reg['intercept']:.5g}; r = {reg['r']:.5f}; "
                            f"range {reg['range'][0]:g}–{reg['range'][1]:g}; LOD {fmt(reg['lod'])}; "
                            f"LOQ {fmt(reg['loq'])}")
            for k in ("model_type", "compounds", "ranges", "n_components", "signals"):
                if d.get(k):
                    info.append(f"{k}: {e(json.dumps(d[k]))}")
            parts.append(f"<tr><td>{m['id']}</td><td>{e(m['name'])}</td><td>{m['type']}</td>"
                         f"<td>{'<br>'.join(info)}</td></tr>")
        parts.append("</table>")
    results = project.results(trial_id)
    if results:
        parts.append("<h2>Results</h2>")
        for r in results:
            parts.append(f"<p><b>#{r['id']} {e(r['name'])}</b> ({r['kind']}, {r['created_utc'][:19]})</p>")
            parts.append(_result_table(project, r["data"]))
    ids = {r.id for r in recs}
    entries = [a for a in reversed(project.audit_entries())
               if trial_id is None or a.details.get("trial_id") == trial_id
               or (a.entity == "trial" and a.entity_id == trial_id)
               or (a.entity == "spectrum" and a.entity_id in ids)
               or ids & set(a.details.get("inputs", []) + a.details.get("outputs", []))]
    parts.append(f"<h2>Audit trail ({len(entries)} entries)</h2><table><tr><th>#</th>"
                 "<th>Time (UTC)</th><th>User</th><th>Action</th><th>Summary</th><th>Reason</th>"
                 "</tr>")
    for a in entries:
        parts.append(f"<tr><td>{a.seq}</td><td>{a.ts_utc[:19].replace('T', ' ')}</td>"
                     f"<td>{e(a.user)}</td><td>{a.action}</td><td>{e(a.summary)}</td>"
                     f"<td>{e(a.reason)}</td></tr>")
    parts.append("</table></body></html>")
    return "".join(parts)


def _result_table(project, data: dict) -> str:
    e = html.escape
    ids, found = data.get("ids"), data.get("found")
    if not ids or found is None or not isinstance(found, list):
        return f"<pre>{e(json.dumps(data, indent=1)[:3000])}</pre>"
    comps = data.get("compounds") or [data.get("method", {}).get("compound", "")]
    rows = []
    for sid, f in zip(ids, found):
        try:
            rec = project.record(int(sid))
        except Exception:
            continue
        vals = f if isinstance(f, list) else [f.get("X"), f.get("Y")] if isinstance(f, dict) else [f]
        cells = [f"<td>{e(rec.name)}</td>"]
        for c, v in zip(comps if not isinstance(f, dict) else [data.get("X"), data.get("Y")], vals):
            t = rec.concentrations.get(c)
            cells.append(f"<td class='num'>{fmt(v)}</td><td class='num'>{'' if t is None else fmt(t)}</td>"
                         f"<td class='num'>{'' if not t or v is None else f'{100 * v / t:.2f}'}</td>")
        rows.append("<tr>" + "".join(cells) + "</tr>")
    names = comps if not isinstance(found[0], dict) else [data.get("X"), data.get("Y")]
    head = "<tr><th>Spectrum</th>" + "".join(
        f"<th>{e(str(c))} found</th><th>taken</th><th>rec. %</th>" for c in names) + "</tr>"
    return "<table>" + head + "".join(rows) + "</table>"


def export_report(win) -> None:
    project = win.project
    trials = project.trials()
    labels = [t["name"] for t in trials] + ["(all spectra)"]
    cur = win.current_trial()
    idx = next((i for i, t in enumerate(trials) if t["id"] == cur), 0)
    choice, ok = QInputDialog.getItem(win, "Trial report", "Trial:", labels, idx, False)
    if not ok:
        return
    tid = None if choice == "(all spectra)" else trials[labels.index(choice)]["id"]
    path, _ = QFileDialog.getSaveFileName(win, "Save report", f"{choice}.pdf",
                                          "PDF (*.pdf);;HTML (*.html)")
    if not path:
        return
    try:
        doc_html = build_html(project, tid)
        if path.lower().endswith(".html"):
            Path(path).write_text(doc_html, encoding="utf-8")
        else:
            if not path.lower().endswith(".pdf"):
                path += ".pdf"
            writer = QPdfWriter(path)
            writer.setPageSize(QPageSize(QPageSize.A4))
            writer.setPageLayout(QPageLayout(QPageSize(QPageSize.A4), QPageLayout.Portrait,
                                             QMarginsF(15, 15, 15, 15)))
            writer.setResolution(150)
            doc = QTextDocument()
            doc.setHtml(doc_html)
            doc.setPageSize(writer.pageLayout().paintRectPixels(150).size().toSizeF())
            doc.print_(writer)
    except Exception as exc:
        error(win, exc)
        return
    project.log("EXPORT", f"Exported report for '{choice}'", {"file": path, "trial_id": tid})
    win.statusBar().showMessage(f"Report saved to {path}")
