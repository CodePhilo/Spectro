"""Build a ready-to-explore demo project from the bundled demo data.

Everything is done through the normal project API, so the demo project's audit
trail shows exactly how it was produced (imports, processing, calibrations,
predictions).

    python -m spectro.demo "Spectro demo.spectro"
"""

from __future__ import annotations

import shutil
import sys
from pathlib import Path
from typing import Callable

from spectro.core import univariate as uv
from spectro.core.io import load_file
from spectro.core.multicomponent import SignalEquations, SpectralModel
from spectro.core.naming import concentrations_from_name
from spectro.storage.project import Project

DATA = Path(__file__).resolve().parent / "demo_data"
COMPOUNDS = [("PAR", 151.16, "Paracetamol (acetaminophen)"),
             ("CAF", 194.19, "Caffeine"),
             ("ASA", 180.16, "Aspirin (acetylsalicylic acid)")]
NOTE = ("Simulated from published band positions and absorptivities with realistic "
        "noise, baseline and preparation errors — not measured data.")


def data_dir() -> Path:
    # PyInstaller one-folder builds put package data under sys._MEIPASS
    base = Path(getattr(sys, "_MEIPASS", "")) / "spectro" / "demo_data"
    return base if base.exists() else DATA


def copy_demo_files(target: str | Path) -> Path:
    target = Path(target) / "Spectro demo data"
    shutil.copytree(data_dir(), target, dirs_exist_ok=True)
    return target


def _import(p: Project, rel: str, trial: int, role: str, comps: list[str],
            fixed: dict[str, float] | None = None, names: list[str] | None = None) -> list[int]:
    path = data_dir() / rel
    spectra = load_file(path).spectra
    conc, new_names = {}, {}
    for i, s in enumerate(spectra):
        c = {k: 0.0 for k in comps}
        c.update(concentrations_from_name(s.name, comps))
        if fixed is not None:
            c = dict(fixed)
        conc[i] = c
        if names:
            new_names[i] = names[i]
    return p.import_file(path, trial, role=role, concentrations=conc, names=new_names)


def build_demo_project(path: str | Path,
                       progress: Callable[[str], None] | None = None) -> Project:
    say = progress or (lambda msg: None)
    p = Project.create(path, "Spectro demo — paracetamol / caffeine / aspirin", NOTE)
    for name, mw, full in COMPOUNDS:
        p.add_compound(name, "µg/mL", mw, full)

    # ------------------------------------------------------------------ binary
    say("Importing binary mixture data…")
    tb = p.add_trial("Binary: paracetamol + caffeine (Panadol Extra)",
                     "PAR 2–20 µg/mL, CAF 2–20 µg/mL in methanol, 1 cm. " + NOTE,
                     {"solvent": "methanol", "cell": "1 cm quartz"})
    par = []
    for f in sorted((data_dir() / "binary" / "shimadzu_txt").glob("*.txt")):
        par += _import(p, f"binary/shimadzu_txt/{f.name}", tb, "standard", ["PAR", "CAF"])
    caf = _import(p, "binary/CAF_standards_xy_pairs.csv", tb, "standard", ["PAR", "CAF"])
    mix = _import(p, "binary/lab_mixtures.xlsx", tb, "mixture", ["PAR", "CAF"])
    tabs = _import(p, "binary/panadol_extra_tablets.csv", tb, "sample", ["PAR", "CAF"],
                   fixed={"PAR": 15.3846, "CAF": 2.0},
                   names=[f"Panadol Extra tablet {i} (label claim)" for i in (1, 2, 3)])
    _import(p, "binary/methanol_blank.txt", tb, "blank", [], fixed={})
    for sid in tabs:
        p.update_spectrum(sid, "demo setup",
                          notes="Concentrations are the nominal (label-claim) values, so "
                                "'recovery %' reads as % of label claim.")
    par10 = next(i for i in par if p.record(i).concentrations["PAR"] == 10)
    caf10 = next(i for i in caf if p.record(i).concentrations["CAF"] == 10)
    p.update_spectrum(caf10, "used as divisor", role="divisor")
    p.update_spectrum(par10, "used as divisor", role="divisor")

    say("Processing spectra…")
    p.process(mix, [{"op": "divide", "params": {"reference": caf10}}], "÷ CAF 10 (ratio)")
    p.process(mix, [{"op": "divide", "params": {"reference": par10}}], "÷ PAR 10 (ratio)")
    p.process(par + caf, [{"op": "derivative", "params": {"order": 1, "delta_lambda": 4.0}}],
              "D1 (Δλ 4 nm)")

    resolve = p.resolver()
    unknown = mix + tabs
    say("Calibrating univariate methods…")
    for comp, cal_ids, div, (w1, w2) in (("PAR", par, caf10, (243.0, 265.0)),
                                         ("CAF", caf, par10, (260.0, 275.0))):
        m = uv.UnivariateMethod(f"Ratio difference — {comp}", comp,
                                [{"op": "divide", "params": {"reference": div}}],
                                {"kind": "difference", "params": {"w1": w1, "w2": w2}})
        m.calibrate(p.spectra(cal_ids), resolve)
        d = m.to_dict()
        d["calibration_ids"] = cal_ids
        mid = p.save_method(m.name, d, tb)
        found = [m.predict(s, resolve) for s in p.spectra(unknown)]
        p.save_result(f"{m.name}: lab mixtures and tablets", "univariate",
                      {"method": m.to_dict(), "ids": unknown, "compounds": [comp],
                       "found": [[f] for f in found]}, tb, mid, unknown)

    say("Fitting Vierordt equations…")
    eq = SignalEquations(["PAR", "CAF"], [{"kind": "amplitude", "params": {"w1": 243.0}},
                                          {"kind": "amplitude", "params": {"w1": 273.0}}])
    std_ids = par + caf
    eq.fit(p.spectra(std_ids))
    d = eq.to_dict()
    d["calibration_ids"] = std_ids
    mid = p.save_method("Vierordt simultaneous equations (243 / 273 nm)", d, tb)
    p.save_result("Vierordt: lab mixtures and tablets", "equations",
                  {"ids": unknown, "compounds": ["PAR", "CAF"],
                   "found": eq.predict(p.spectra(unknown)).tolist()}, tb, mid, unknown)

    # ------------------------------------------------------------------ ternary
    say("Importing ternary mixture data…")
    tt = p.add_trial("Ternary: aspirin + paracetamol + caffeine (Excedrin)",
                     "Brereton 5-level design training set (25 mixtures), 5 validation "
                     "mixtures, 3 tablet solutions. " + NOTE)
    comps = ["PAR", "CAF", "ASA"]
    train = _import(p, "ternary/training_set_brereton.jdx", tt, "calibration", comps)
    val = _import(p, "ternary/validation_mixtures.tsv", tt, "validation", comps)
    exc = _import(p, "ternary/excedrin_tablets.spc", tt, "sample", comps,
                  fixed={"PAR": 12.0, "CAF": 3.12, "ASA": 12.0},
                  names=[f"Excedrin tablet {i} (label claim)" for i in (1, 2, 3)])
    _import(p, "ternary/ASA_standards.csv", tt, "standard", comps)

    say("Building PLS model…")
    pls = SpectralModel("PLS2", comps, ranges=[(215.0, 310.0)], n_components=3)
    cv = pls.cross_validate(p.spectra(train), method="venetian", folds=5, max_components=8)
    p.log("CALCULATE", "Cross-validated PLS2 (demo)",
          {"inputs": train, "curve": cv["curve"], "suggested": cv["suggested_components"]})
    pls.n_components = cv["suggested_components"]
    pls.fit(p.spectra(train))
    d = pls.to_dict(include_fit=True)
    d["calibration_ids"] = train
    mid = p.save_method(f"PLS2 ({pls.n_components} LVs, 215–310 nm)", d, tt)
    targets = val + exc
    p.save_result("PLS2: validation mixtures and tablets", "chemometrics",
                  {"ids": targets, "compounds": comps,
                   "found": pls.predict(p.spectra(targets)).tolist()}, tt, mid, targets)
    say("Demo project ready.")
    return p


def main(argv: list[str] | None = None) -> int:
    argv = sys.argv[1:] if argv is None else argv
    path = Path(argv[0] if argv else "Spectro demo.spectro")
    if path.exists():
        print(f"{path} already exists", file=sys.stderr)
        return 1
    p = build_demo_project(path, print)
    rep = p.verify()
    p.close()
    print(f"{path}: {rep['spectra']} spectra, {rep['audit_entries']} audit entries, "
          f"integrity {'OK' if rep['ok'] else 'FAILED'}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
