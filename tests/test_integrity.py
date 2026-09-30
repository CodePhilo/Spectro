"""Import edge cases, project/audit integrity and end-to-end numerical flows."""

import sqlite3
import time

import numpy as np
import pytest

from spectro.core.io import export_csv, export_excel, load_file
from spectro.core.spectrum import Spectrum
from spectro.demo import build_demo_project
from spectro.storage.project import Project
from tests.conftest import GRID, mixture

X = np.arange(200.0, 401.0, 1.0)


def gauss(c=260, w=20, h=0.8):
    return h * np.exp(-0.5 * ((X - c) / w) ** 2)


# --------------------------------------------------------------------------- #
# Import edge cases
# --------------------------------------------------------------------------- #
@pytest.mark.parametrize("content", ["", "\n\n", "Wavelength,Abs\n", "only text\nno numbers",
                                     "1\n2\n3\n4\n5"])
def test_unusable_files_give_clear_errors(tmp_path, content):
    f = tmp_path / "bad.csv"
    f.write_text(content)
    with pytest.raises(ValueError):
        load_file(f)


def test_windows_line_endings_bom_scientific_and_unicode_minus(tmp_path):
    f = tmp_path / "w.csv"
    rows = ["﻿nm,Sample β"] + [f"{a:.1f},{v:.4E}".replace("-", "−")
                                    for a, v in zip(X, gauss() - 0.01)]
    f.write_bytes("\r\n".join(rows).encode("utf-8"))
    s = load_file(f).spectra[0]
    assert s.name == "Sample β" and s.x.size == X.size
    assert np.allclose(s.y, gauss() - 0.01, atol=1e-4)


def test_footer_text_and_blank_lines_are_ignored(tmp_path):
    f = tmp_path / "f.txt"
    f.write_text("Instrument: X\n\nnm\tA\n" + "\n".join(f"{a}\t{v:.5f}" for a, v in
                                                         zip(X, gauss()))
                 + "\n\nEnd of data\nPrinted 2026-01-01\n")
    s = load_file(f).spectra[0]
    assert s.x.size == X.size


def test_duplicate_wavelengths_are_averaged(tmp_path):
    f = tmp_path / "dup.csv"
    xs = np.concatenate([X[:100], X[99:]])        # 299 nm measured twice (overlap)
    ys = np.concatenate([gauss()[:100], gauss()[99:] + np.r_[0.002, np.zeros(101)]])
    f.write_text("nm,A\n" + "\n".join(f"{a},{v:.5f}" for a, v in zip(xs, ys)))
    s = load_file(f).spectra[0]
    assert np.all(np.diff(s.wavelengths) > 0) and s.x.size == X.size
    assert s.value_at(299.0) == pytest.approx(gauss()[99] + 0.001, abs=1e-5)


def test_missing_cells_in_excel_are_handled(tmp_path):
    import openpyxl
    wb = openpyxl.Workbook()
    ws = wb.active
    ws.append(["nm", "A", "B"])
    for i, (a, v) in enumerate(zip(X, gauss())):
        ws.append([a, v, None if i % 50 == 7 else v * 2])
    wb.create_sheet("empty")
    wb.save(tmp_path / "m.xlsx")
    r = load_file(tmp_path / "m.xlsx")
    assert [s.name for s in r.spectra] == ["A", "B"]
    assert np.all(np.isfinite(r.spectra[1].y)) and r.spectra[1].x.size == X.size - 4


def test_large_file_parses_quickly(tmp_path):
    grid = np.arange(190, 1100.5, 0.5)
    spectra = [Spectrum(grid, np.random.default_rng(i).random(grid.size),
                        name=f"S{i}") for i in range(120)]
    export_csv(spectra, tmp_path / "big.csv")
    t = time.perf_counter()
    r = load_file(tmp_path / "big.csv")
    assert len(r.spectra) == 120 and time.perf_counter() - t < 15


def test_unicode_and_special_names_round_trip(tmp_path):
    names = ["Ibuprofène 10 µg/mL", 'Mix "A";1', "Долг/тест", "S:[1]*?"]
    spectra = [Spectrum(X, gauss(250 + i * 10), name=n) for i, n in enumerate(names)]
    export_excel(spectra, tmp_path / "u.xlsx")
    export_csv(spectra, tmp_path / "u.csv")
    for f in ("u.xlsx", "u.csv"):
        assert [s.name for s in load_file(tmp_path / f).spectra] == names, f


# --------------------------------------------------------------------------- #
# Project integrity
# --------------------------------------------------------------------------- #
@pytest.fixture
def proj(tmp_path):
    p = Project.create(tmp_path / "p.spectro")
    yield p
    if p.conn:
        p.close()


def _count(p):
    return p.conn.execute("SELECT COUNT(*) FROM audit").fetchone()[0]


def test_every_mutation_writes_exactly_one_audit_entry(proj, tmp_path):
    f = tmp_path / "d.csv"
    export_csv([mixture({"X": 10}, "A"), mixture({"Y": 5}, "B")], f)
    calls = [
        lambda: proj.add_compound("X"),
        lambda: proj.update_compound(1, "r", unit="mg/L"),
        lambda: proj.add_trial("T"),
        lambda: proj.update_trial(1, "r", description="d"),
        lambda: proj.import_file(f, 1),
        lambda: proj.add_spectrum(mixture({"X": 3}, "C"), 1),
        lambda: proj.process([1], [{"op": "scale", "params": {"factor": 2}}]),
        lambda: proj.update_spectrum(1, "r", name="A2"),
        lambda: proj.set_concentrations_bulk({1: {"X": 10.0}, 2: {"Y": 5.0}}, "r"),
        lambda: proj.archive_spectra([2], "bad"),
        lambda: proj.save_method("m", {"type": "univariate"}),
        lambda: proj.archive_method(1, "old"),
        lambda: proj.save_result("r", "k", {"a": 1}),
        lambda: proj.log("EXPORT", "x", {}),
        lambda: proj.set_meta("name", "New"),
        lambda: proj.archive_trial(1, "done"),
    ]
    for i, call in enumerate(calls):
        before = _count(proj)
        call()
        assert _count(proj) == before + 1, f"call #{i} wrote {_count(proj) - before} entries"
    # no-op edits write nothing
    before = _count(proj)
    proj.update_spectrum(1, "r", name="A2")
    proj.update_compound(1, "r", unit="mg/L")
    proj.set_meta("name", "New")
    assert _count(proj) == before
    assert proj.verify()["ok"]


def test_failed_operations_leave_no_trace(proj, tmp_path):
    f = tmp_path / "d.csv"
    export_csv([mixture({"X": 10}, "A")], f)
    tid = proj.add_trial("T")
    proj.import_file(f, tid)
    n_audit, n_spec = _count(proj), len(proj.records(include_archived=True))
    for bad in ([{"op": "divide", "params": {"reference": 999}}],
                [{"op": "crop", "params": {"start": 500, "end": 600}}],
                [{"op": "nope", "params": {}}]):
        with pytest.raises((KeyError, ValueError)):
            proj.process([1], bad)
    with pytest.raises(Exception):
        proj.import_file(tmp_path / "missing.csv", tid)
    assert _count(proj) == n_audit and len(proj.records(include_archived=True)) == n_spec
    assert proj.verify()["ok"]


def test_tampering_with_spectral_data_or_deleting_log_is_detected(proj, tmp_path):
    f = tmp_path / "d.csv"
    export_csv([mixture({"X": 10}, "A")], f)
    proj.import_file(f, proj.add_trial("T"))
    c = proj.conn
    c.execute("DROP TRIGGER spectra_immutable_data")
    y = np.frombuffer(c.execute("SELECT ydata FROM spectra WHERE id=1").fetchone()[0], "<f8").copy()
    y[100] *= 1.01
    c.execute("UPDATE spectra SET ydata=? WHERE id=1", (y.astype("<f8").tobytes(),))
    rep = proj.verify()
    assert not rep["ok"] and any("spectrum #1" in m for m in rep["problems"])
    c.execute("DROP TRIGGER audit_no_delete")
    c.execute("DELETE FROM audit WHERE seq=2")
    assert any("chain broken" in m for m in proj.verify()["problems"])


def test_two_app_instances_on_one_file_keep_a_valid_chain(tmp_path):
    path = tmp_path / "shared.spectro"
    a = Project.create(path)
    b = Project.open(path)
    try:
        for i in range(10):
            (a if i % 2 else b).add_compound(f"C{i}")
        assert a.verify()["ok"] and b.verify()["ok"]
        assert len(a.compounds()) == 10
    finally:
        a.close()
        b.close()


def test_opening_foreign_files_is_refused(tmp_path):
    (tmp_path / "x.spectro").write_text("not a database")
    with pytest.raises((ValueError, sqlite3.DatabaseError)):
        Project.open(tmp_path / "x.spectro")
    other = sqlite3.connect(tmp_path / "other.spectro")
    other.execute("CREATE TABLE t (a)")
    other.commit()
    other.close()
    with pytest.raises(ValueError):
        Project.open(tmp_path / "other.spectro")


def test_chained_derivations_replay_even_after_archiving_references(proj):
    tid = proj.add_trial("T")
    a = proj.add_spectrum(mixture({"X": 10, "Y": 5}, "mix"), tid)
    y = proj.add_spectrum(mixture({"Y": 10}, "Y10"), tid)
    r1 = proj.process([a], [{"op": "divide", "params": {"reference": y}}])[0]
    dy = proj.process([y], [{"op": "derivative", "params": {"order": 1}}])[0]
    r2 = proj.process([r1], [{"op": "derivative", "params": {"order": 1}},
                             {"op": "crop", "params": {"start": 220, "end": 300}}])[0]
    r3 = proj.process([a], [{"op": "crop", "params": {"start": 220, "end": 390}},
                            {"op": "subtract", "params": {"reference": dy, "factor": 0.1}}])[0]
    proj.archive_spectra([y, dy], "cleanup")
    for sid in (r1, r2, r3):
        assert proj.replay(sid)["ok"], sid
    assert [r.id for r in proj.lineage(r2)] == [a, r1, r2]


# --------------------------------------------------------------------------- #
# End-to-end: stored results equal an independent hand calculation
# --------------------------------------------------------------------------- #
def test_demo_results_match_independent_numpy_calculation(tmp_path):
    p = build_demo_project(tmp_path / "d.spectro")
    try:
        res = next(r for r in p.results() if r["name"].startswith("Ratio difference — PAR"))
        recs = p.records()
        caf10 = next(r for r in recs if r.name == "CAF 10")
        div = p.spectrum(caf10.id)
        std = [r for r in recs if r.kind == "raw" and r.name.startswith("PAR ")
               and r.name.endswith("ug/mL")]   # all 8 standards (PAR 10 is also the divisor)

        def signal(s):  # ratio difference at 243 / 265 nm, written from scratch
            ratio = s.y / np.interp(s.x, div.x, div.y)
            return np.interp(243.0, s.x, ratio) - np.interp(265.0, s.x, ratio)
        c = np.array([r.concentrations["PAR"] for r in std])
        sig = np.array([signal(p.spectrum(r.id)) for r in std])
        b, a = np.polyfit(c, sig, 1)
        for sid, found in zip(res["data"]["ids"], res["data"]["found"]):
            assert found[0] == pytest.approx((signal(p.spectrum(sid)) - a) / b, rel=1e-9)
        # every derived spectrum reproduces from raw data
        for r in p.records():
            if r.kind == "derived":
                assert p.replay(r.id)["ok"]
        assert p.verify()["ok"]
    finally:
        p.close()


def test_excel_export_contains_exactly_the_stored_numbers(tmp_path):
    import openpyxl

    from spectro.storage.tables import export_results_excel
    p = build_demo_project(tmp_path / "d.spectro")
    try:
        export_results_excel(p, tmp_path / "r.xlsx")
        wb = openpyxl.load_workbook(tmp_path / "r.xlsx")
        res = next(r for r in p.results() if r["name"].startswith("PLS2"))
        ws = wb[next(n for n in wb.sheetnames if n.startswith(f"{res['id']} "))]
        found = [[ws.cell(row=5 + i, column=2 + 3 * j).value for j in range(3)]
                 for i in range(len(res["data"]["ids"]))]
        assert np.allclose(found, res["data"]["found"], rtol=1e-12)
    finally:
        p.close()


def test_pdf_and_html_reports(tmp_path, qapp):
    from PySide6.QtGui import QPageSize, QPdfWriter, QTextDocument

    from spectro.ui.report import build_html
    p = build_demo_project(tmp_path / "d.spectro")
    try:
        html = build_html(p, p.trials()[0]["id"])
        assert "PASSED" in html and "Ratio difference" in html and html.count("<tr>") > 40
        w = QPdfWriter(str(tmp_path / "r.pdf"))
        w.setPageSize(QPageSize(QPageSize.A4))
        doc = QTextDocument()
        doc.setHtml(html)
        doc.print_(w)
        del w
        assert (tmp_path / "r.pdf").stat().st_size > 5000
    finally:
        p.close()


@pytest.fixture(scope="module")
def qapp():
    from PySide6.QtWidgets import QApplication
    return QApplication.instance() or QApplication([])


def test_spectrum_rejects_nan_and_mismatched_data():
    with pytest.raises(ValueError):
        Spectrum([1, 2, 3], [1, 2])
    with pytest.raises(ValueError):
        Spectrum([1], [1])
    try:
        s = Spectrum(GRID[:5], [0.1, np.nan, 0.2, 0.3, 0.4])
    except ValueError:
        return
    assert np.all(np.isfinite(s.values)), "NaN accepted into a spectrum"
