from spectro.core.io import load_file
from spectro.demo import build_demo_project, data_dir


def test_every_demo_file_imports():
    files = [f for f in data_dir().rglob("*") if f.is_file() and f.suffix != ".md"]
    assert len(files) == 16
    for f in files:
        r = load_file(f)
        assert r.spectra and all(s.x[0] == 200 and s.x[-1] == 400 for s in r.spectra), f.name


def test_demo_project_recoveries(tmp_path):
    p = build_demo_project(tmp_path / "demo.spectro")
    try:
        assert p.verify()["ok"]
        results = p.results()
        assert len(results) == 4 and len(p.methods()) == 4
        for r in results:
            d = r["data"]
            for sid, found in zip(d["ids"], d["found"]):
                conc = p.record(sid).concentrations
                for c, v in zip(d["compounds"], found):
                    if conc.get(c):
                        assert abs(100 * v / conc[c] - 100) < 3, (r["name"], sid, c)
    finally:
        p.close()
