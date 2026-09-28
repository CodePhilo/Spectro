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
