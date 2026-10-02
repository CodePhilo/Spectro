import sqlite3

import numpy as np
import pytest

from spectro.core.io import export_csv
from spectro.storage.project import Project
from tests.conftest import mixture


@pytest.fixture
def project(tmp_path):
    p = Project.create(tmp_path / "trial.spectro", "Test")
    yield p
    if p.conn:
        p.close()


def _import_two(project, tmp_path):
    f = tmp_path / "data.csv"
    export_csv([mixture({"X": 10}, "A"), mixture({"Y": 5}, "B")], f)
    tid = project.add_trial("Binary PAR+CAF")
    return tid, project.import_file(f, tid, role="standard")


def test_import_process_and_audit(project, tmp_path):
    tid, ids = _import_two(project, tmp_path)
    assert [r.name for r in project.records(tid)] == ["A", "B"]
    out = project.process([ids[0]], [{"op": "divide", "params": {"reference": ids[1]}},
                                      {"op": "derivative", "params": {"order": 1}}])
    rec = project.record(out[0])
    assert rec.kind == "derived" and rec.parent_ids == [ids[0], ids[1]]
    assert project.replay(out[0])["ok"]
    assert [r.id for r in project.lineage(out[0])] == [ids[0], out[0]]
    actions = [e.action for e in project.audit_entries()]
    assert {"CREATE", "IMPORT", "PROCESS", "VERIFY"} <= set(actions)
    hist = project.spectrum_history(ids[0])
    assert [e.action for e in hist][:2] == ["IMPORT", "PROCESS"]
    assert project.verify()["ok"]


def test_edits_are_logged_with_before_after(project, tmp_path):
    tid, ids = _import_two(project, tmp_path)
    project.update_spectrum(ids[0], reason="typo", name="Std A",
                            concentrations={"X": 10.0, "Y": 0.0})
    e = project.audit_entries(action="EDIT")[0]
    assert e.reason == "typo"
    assert e.details["before"]["name"] == "A" and e.details["after"]["name"] == "Std A"


def test_raw_data_and_audit_are_immutable(project, tmp_path):
    tid, ids = _import_two(project, tmp_path)
    c = project.conn
    with pytest.raises(sqlite3.DatabaseError):
        c.execute("UPDATE spectra SET ydata=? WHERE id=?", (b"x", ids[0]))
    with pytest.raises(sqlite3.DatabaseError):
        c.execute("DELETE FROM spectra WHERE id=?", (ids[0],))
    with pytest.raises(sqlite3.DatabaseError):
        c.execute("UPDATE audit SET summary='x'")
    with pytest.raises(sqlite3.DatabaseError):
        c.execute("DELETE FROM audit")


def test_tampering_is_detected(project, tmp_path):
    _import_two(project, tmp_path)
    c = project.conn
    # A malicious user could drop the trigger and edit the log directly …
    c.execute("DROP TRIGGER audit_no_update")
    c.execute("UPDATE audit SET summary='nothing happened' WHERE action='IMPORT'")
    report = project.verify()
    assert not report["ok"]
    assert any("modified" in p for p in report["problems"])


def test_archive_and_reopen(project, tmp_path):
    tid, ids = _import_two(project, tmp_path)
    project.archive_spectra([ids[1]], reason="bad scan")
    assert [r.id for r in project.records(tid)] == [ids[0]]
    path = project.path
    project.close()
    p2 = Project.open(path)
    try:
        assert len(p2.records(tid, include_archived=True)) == 2
        assert np.allclose(p2.spectrum(ids[0]).values, mixture({"X": 10}).values)
        assert p2.verify()["ok"]
    finally:
        p2.close()


def test_methods_and_results(project):
    mid = project.save_method("RD for X", {"type": "univariate", "name": "RD"})
    rid = project.save_result("Assay", "prediction", {"X": [10.1]}, method_id=mid)
    assert project.methods()[0]["id"] == mid
    assert project.results()[0]["id"] == rid
    with pytest.raises(sqlite3.DatabaseError):
        project.conn.execute("UPDATE results SET data='{}'")


@pytest.mark.parametrize("kind", ["amplitude_centering", "absorption_factor"])
def test_progressive_methods_save_apply_and_travel_as_model_files(tmp_path, kind):
    """A saved progressive method is applied to new spectra after reopening
    (refitting from its calibration ids if needed) and survives export to a
    .spmodel file and import into another project (the divisor is embedded)."""
    from spectro.core import univariate as uv
    from spectro.storage.modelfile import export_model, import_model, is_calibrated
    from tests.test_literature import af, ter

    mk, comps = (ter, ["X", "Y", "Z"]) if kind == "amplitude_centering" else \
        (af, ["MET", "MEB", "DLX"])
    p = Project.create(tmp_path / "a.spectro")
    tid = p.add_trial("T")
    cal = [p.add_spectrum(mk({c: v}, f"{c} {v}"), tid, role="standard")
           for c in comps for v in (4, 8, 12, 16, 20)]
    mixes = [{c: v for c, v in zip(comps, vals)} for vals in ((10, 10, 10), (20, 6, 4))]
    mix_ids = [p.add_spectrum(mk(m, f"mix {i}"), tid, role="mixture") for i, m in enumerate(mixes)]
    if kind == "amplitude_centering":
        div = p.add_spectrum(mk({"Z": 24}, "Z′ 24"), tid, role="divisor")
        m = uv.AmplitudeCentering(275.0, comps, subtract="Y", divisor_compound="Z",
                                  plateau=(340.0, 380.0),
                                  differences={"X": {"w1": 275.0, "w2": 240.0,
                                                     "factor_from": "Y"}}, divisor=div)
    else:
        m = uv.AbsorptionFactorMethod([("MET", 330.0), ("MEB", 280.0), ("DLX", 246.0)])
    m.fit_spectra(p.spectra(cal), p.resolver())
    d = m.to_dict()
    d["calibration_ids"] = cal
    mid = p.save_method("progressive", d, tid)
    true = np.array([[x[c] for c in comps] for x in mixes])

    # a definition without its fit (e.g. edited) is refitted from calibration_ids
    stored = p.methods()[-1]["definition"]
    assert is_calibrated(stored)
    bare = uv.progressive_from_dict({**stored, "regressions": {}})
    assert not bare.is_fitted
    bare.fit_spectra(p.spectra(stored["calibration_ids"]), p.resolver())
    for meth in (uv.progressive_from_dict(stored), bare):
        assert np.allclose(meth.predict_spectra(p.spectra(mix_ids), p.resolver()), true,
                           rtol=1e-3)

    doc = export_model(p, next(x for x in p.methods() if x["id"] == mid))
    if kind == "amplitude_centering":
        assert list(doc["embedded_spectra"]) == [str(div)]
    q = Project.create(tmp_path / "b.spectro")
    t2 = q.add_trial("T")
    new_mid = import_model(q, doc, t2)
    qd = next(x for x in q.methods() if x["id"] == new_mid)["definition"]
    qm = uv.progressive_from_dict(qd)
    qmix = [q.add_spectrum(mk(x, "mix"), t2, role="mixture") for x in mixes]
    assert np.allclose(qm.predict_spectra(q.spectra(qmix), q.resolver()), true, rtol=1e-3)
    p.close()
    q.close()


def test_revising_a_method_creates_an_audited_new_version(project):
    mid = project.save_method("RD for X", {"type": "univariate", "compound": "X",
                                           "steps": [], "measurement": {"kind": "amplitude",
                                                                        "params": {"w1": 250}}})
    with pytest.raises(ValueError, match="reason"):
        project.revise_method(mid, "RD for X", {"type": "univariate"}, "  ")
    d2 = dict(project.method(mid)["definition"])
    d2["measurement"] = {"kind": "amplitude", "params": {"w1": 252}}
    v2 = project.revise_method(mid, "RD for X (252 nm)", d2, "better λ")
    v3 = project.revise_method(v2, "RD for X (252 nm)", project.method(v2)["definition"],
                               "renamed nothing, re-saved")
    live = project.methods()
    assert [m["id"] for m in live] == [v3]
    assert project.method(mid)["archived"] and project.method(v2)["archived"]
    assert project.method(v3)["definition"]["version"] == 3
    assert project.method(v3)["definition"]["revision_of"] == v2
    assert project.method(v2)["definition"]["measurement"]["params"]["w1"] == 252
    assert project.method(mid)["definition"]["measurement"]["params"]["w1"] == 250  # kept
    rev = [e for e in project.audit_entries() if e.action == "REVISE"]
    assert len(rev) == 2 and {e.reason for e in rev} == {"better λ", "renamed nothing, re-saved"}
    assert project.verify()["ok"]
    with pytest.raises(KeyError):
        project.revise_method(mid, "x", d2, "an archived version cannot be revised")


def test_progressive_methods_with_smoothing_steps(tmp_path):
    """Smoothing is applied to standards, divisor and samples alike and is saved
    with the method."""
    from spectro.core import univariate as uv
    from tests.test_literature import TER, ter
    comps = list(TER)
    std = [ter({c: v}, f"{c} {v}", noise=0.002, seed=int(v) + 7 * k)
           for k, c in enumerate(comps) for v in (4, 8, 12, 16, 20, 24)]
    div = ter({"Z": 24}, "Z′", noise=0.0007, seed=99)
    mixes = [ter({"X": 10, "Y": 10, "Z": 10}, "m1", noise=0.002, seed=5),
             ter({"X": 20, "Y": 10, "Z": 10}, "m2", noise=0.002, seed=6)]
    sg = [{"op": "smooth_sg", "params": {"window": 31, "polyorder": 2}}]
    out = {}
    for steps in ([], sg):
        m = uv.AmplitudeCentering(275.0, comps, subtract="Y", divisor_compound="Z",
                                  plateau=(340.0, 380.0),
                                  differences={"X": {"w1": 275.0, "w2": 240.0,
                                                     "factor_from": "Y"}},
                                  divisor=div, steps=steps)
        m.fit_spectra(std)
        out[bool(steps)] = m.predict_spectra(mixes)
        if steps:
            rt = uv.progressive_from_dict({**m.to_dict(), "divisor": None})
            assert rt.steps == sg and "SG smoothing (31 pts" in rt.describe()
    true = np.array([[10, 10, 10], [20, 10, 10]])
    err = {k: np.abs(v / true - 1).max() for k, v in out.items()}
    assert err[True] < 0.5 * err[False] and err[True] < 0.05   # 10.3 % → 3.4 % here
    af = uv.AbsorptionFactorMethod([("X", 340.0)], steps=sg)
    assert uv.progressive_from_dict(af.to_dict()).steps == sg


def test_changing_processing_rebuilds_downstream_and_keeps_history(project, tmp_path):
    tid, ids = _import_two(project, tmp_path)
    a, b = ids
    sm = project.process([b], [{"op": "smooth_sg", "params": {"window": 5, "order": 2}}])[0]
    ratio = project.process([a], [{"op": "divide", "params": {"reference": sm}}])[0]
    d1 = project.process([ratio], [{"op": "derivative", "params": {"order": 1}}])[0]
    with pytest.raises(ValueError):
        project.revise_processing([sm], [{"op": "smooth_sg", "params": {"window": 9}}], " ")
    with pytest.raises(ValueError):
        project.revise_processing([a], [{"op": "derivative", "params": {}}], "raw")
    mapping = project.revise_processing(
        [sm], [{"op": "smooth_sg", "params": {"window": 9, "order": 2}}], "wider window")
    assert set(mapping) == {sm, ratio, d1}
    live = {r.id for r in project.records(tid)}
    assert not {sm, ratio, d1} & live and set(mapping.values()) <= live
    new_ratio = project.record(mapping[ratio])
    assert new_ratio.parent_ids == [a, mapping[sm]]
    assert new_ratio.pipeline[0]["params"]["reference"] == mapping[sm]
    assert project.record(mapping[d1]).parent_ids == [mapping[ratio]]
    assert project.record(mapping[sm]).pipeline[0]["params"]["window"] == 9
    for new in mapping.values():
        assert project.replay(new)["ok"]
    e = project.audit_entries(action="REVISE")[0]
    assert e.reason == "wider window" and e.details["replaced"][str(sm)] == mapping[sm]
    assert project.verify()["ok"]


def test_concentration_edits_reach_processed_versions(project, tmp_path):
    tid, ids = _import_two(project, tmp_path)
    d = project.process([ids[0]], [{"op": "derivative", "params": {"order": 1}}])[0]
    dd = project.process([d], [{"op": "derivative", "params": {"order": 1}}])[0]
    project.update_spectrum(ids[0], reason="weighing", concentrations={"X": 12.5})
    assert project.record(d).concentrations == {"X": 12.5}
    assert project.record(dd).concentrations == {"X": 12.5}
    e = project.audit_entries(action="EDIT")[0]
    assert set(e.details["propagated_to"]) == {d, dd}
    project.set_concentrations_bulk({ids[0]: {"X": 11.0}}, reason="re-weighed")
    assert project.record(dd).concentrations == {"X": 11.0}
    mid = project.save_method("m", {"type": "univariate", "calibration_ids": [d]})
    project.save_result("r", "univariate", {"ids": [ids[1]]})
    use = project.usage([d])
    assert [m["id"] for m in use["methods"]] == [mid] and use["results"] == []
    assert len(project.usage([ids[1]])["results"]) == 1


def test_revising_a_result_creates_an_audited_new_version(project):
    rid = project.save_result("Assay", "univariate", {"rows": [1, 2]})
    with pytest.raises(ValueError):
        project.revise_result(rid, "Assay", {"rows": [1]}, "")
    r2 = project.revise_result(rid, "Assay (corrected)", {"rows": [1], "excluded_ids": [2]},
                               "outlier")
    assert [r["id"] for r in project.results()] == [r2]
    new = project.result(r2)
    assert new["data"]["revision_of"] == rid and new["data"]["version"] == 2
    assert project.result(rid)["archived"] and project.result(rid)["data"] == {"rows": [1, 2]}
    with pytest.raises(KeyError):
        project.revise_result(rid, "x", {}, "old version")
    assert project.audit_entries(action="REVISE")[0].reason == "outlier"
    assert project.verify()["ok"]


def _assay(project):
    """Y standards (one processed by D1), mixtures, a D1 method and a result."""
    from spectro.core import univariate as uv
    tid = project.add_trial("T")
    std = [project.add_spectrum(mixture({"Y": v}, f"std {v}"), tid, role="standard")
           for v in (2, 4, 6, 8, 10)]
    mix = [project.add_spectrum(mixture({"X": 3, "Y": v}, f"mix {v}"), tid, role="mixture")
           for v in (3, 5, 7)]
    m = uv.UnivariateMethod("Y direct", "Y", [], {"kind": "amplitude", "params": {"w1": 350.0}})
    m.calibrate(project.spectra(std), project.resolver())
    d = m.to_dict()
    d["calibration_ids"] = std
    mid = project.save_method("Y direct", d, tid)
    found = [[m.predict(s, project.resolver())] for s in project.spectra(mix)]
    rid = project.save_result("Y assay", "routine", {"ids": mix, "compounds": ["Y"],
                                                     "found": found}, tid, mid, mix)
    return tid, std, mix, mid, rid


def test_results_keep_taken_values_and_confidence_intervals(project):
    from spectro.storage.tables import result_table
    _, std, mix, mid, rid = _assay(project)
    data = project.result(rid)["data"]
    assert data["taken_snapshot"][str(mix[0])] == {"Y": 3}
    reg = project.method(mid)["definition"]["regression"]
    from spectro.core.validation import Regression
    assert data["ci95"][0][0] == pytest.approx(
        Regression.from_dict(reg).ci_at_x(data["found"][0][0]))
    headers, rows, summary = result_table(project, data)
    assert headers[:5] == ["Spectrum", "Y found", "Y ± 95 % CI", "Y taken", "Y recovery %"]
    # a later correction does not silently change the stored recovery …
    project.update_spectrum(mix[0], reason="typo", concentrations={"X": 3, "Y": 3.3})
    _, rows2, summary2 = result_table(project, project.result(rid)["data"])
    assert rows2[0][3] == 3 and rows2[0][4] == pytest.approx(rows[0][4])
    assert any("changed after this result was saved" in x for x in summary2)


def test_recalculate_all_affected_after_corrections(project):
    from spectro.storage.recalc import recalculate_affected, stale_taken
    _, std, mix, mid, rid = _assay(project)
    project.update_spectrum(std[-1], reason="weighing", concentrations={"Y": 12.0})
    out = recalculate_affected(project, [std[-1]], "standard re-weighed")
    assert set(out["methods"]) == {mid} and set(out["results"]) == {rid}, out
    new_m = project.method(out["methods"][mid])["definition"]
    assert new_m["regression"]["x"][-1] == 12.0 and new_m["revision_of"] == mid
    new_r = project.result(out["results"][rid])
    assert new_r["method_id"] == out["methods"][mid] and not stale_taken(project, new_r["data"])
    assert new_r["data"]["found"] != project.result(rid)["data"]["found"]

    # processing correction: spectra are replaced, ids remapped in method and result
    d1 = project.process(std + mix, [{"op": "smooth_sg", "params": {"window": 5, "order": 2}}])
    from spectro.core import univariate as uv
    m = uv.UnivariateMethod("Y sm", "Y", [], {"kind": "amplitude", "params": {"w1": 350.0}})
    m.calibrate(project.spectra(d1[:5]), project.resolver())
    d = m.to_dict()
    d["calibration_ids"] = d1[:5]
    mid2 = project.save_method("Y sm", d)
    rid2 = project.save_result("Y sm assay", "routine",
                               {"ids": d1[5:], "compounds": ["Y"],
                                "found": [[m.predict(s)] for s in project.spectra(d1[5:])]},
                               None, mid2, d1[5:])
    mapping = project.revise_processing(d1, [{"op": "smooth_sg",
                                              "params": {"window": 9, "order": 2}}], "wider")
    out = recalculate_affected(project, [], "re-processed", mapping)
    nm = project.method(out["methods"][mid2])["definition"]
    assert nm["calibration_ids"] == [mapping[i] for i in d1[:5]]
    nr = project.result(out["results"][rid2])["data"]
    assert nr["ids"] == [mapping[i] for i in d1[5:]]


def test_undo_restores_previous_versions(project, tmp_path):
    _, std, mix, mid, rid = _assay(project)
    # concentration edit (propagated) → undo
    sm = project.process([mix[0]], [{"op": "smooth_sg", "params": {"window": 5, "order": 2}}])
    project.update_spectrum(mix[0], reason="typo", concentrations={"X": 3, "Y": 4})
    with pytest.raises(ValueError):
        project.undo(project.undoable()[0].seq, "")
    project.undo(project.undoable()[0].seq, "wrong spectrum")
    assert project.record(mix[0]).concentrations == {"X": 3, "Y": 3}
    assert project.record(sm[0]).concentrations == {"X": 3, "Y": 3}
    # result revision → undo
    r2 = project.revise_result(rid, "renamed", project.result(rid)["data"], "name")
    e = project.undoable()[0]
    assert e.action == "REVISE" and e.entity == "result"
    project.undo(e.seq, "keep the old name")
    assert [r["id"] for r in project.results()] == [rid]
    assert project.result(r2)["archived"]
    with pytest.raises(KeyError):
        project.undo(e.seq, "twice")
    # method revision that a later revision depends on cannot be undone first
    m2 = project.revise_method(mid, "v2", project.method(mid)["definition"], "a")
    m3 = project.revise_method(m2, "v3", project.method(m2)["definition"], "b")
    first = next(x for x in project.undoable() if x.entity == "method"
                 and x.details["after"]["id"] == m2)
    with pytest.raises(ValueError):
        project.undo(first.seq, "out of order")
    project.undo(project.undoable()[0].seq, "undo v3")
    project.undo(first.seq, "undo v2")
    assert [m["id"] for m in project.methods()] == [mid]
    assert project.method(m3)["archived"]
    # re-processing → undo
    mp = project.revise_processing(sm, [{"op": "smooth_sg", "params": {"window": 9}}], "w")
    project.undo(project.undoable()[0].seq, "back")
    assert not project.record(sm[0]).archived and project.record(mp[sm[0]]).archived
    # archive → undo
    project.archive_result(rid, "x")
    project.undo(project.undoable()[0].seq, "oops")
    assert [r["id"] for r in project.results()] == [rid]
    assert project.verify()["ok"]
