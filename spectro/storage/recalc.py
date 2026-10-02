"""Recalculating saved results and refitting saved methods on current data.

Used by the result editor, by *Recalculate all affected* after a
concentration or processing correction, and when results are saved (95 %
confidence intervals of the found concentrations, snapshot of the taken
values)."""

from __future__ import annotations

import copy
import math
from typing import Any

from spectro.core import univariate as uv
from spectro.core.multicomponent import SignalEquations, SpectralModel
from spectro.core.validation import Regression

TYPES = ("univariate", "equations", "spectral") + tuple(uv.PROGRESSIVE_TYPES)
_DROP_ON_RECALC = ("X", "Y", "taken", "recovery", "signals", "names", "found_minus_added",
                   "taken_snapshot", "ci95")


def method_from_definition(d: dict):
    t = d.get("type")
    if t in uv.PROGRESSIVE_TYPES:
        return uv.progressive_from_dict(d)
    if t == "univariate":
        return uv.UnivariateMethod.from_dict(d)
    if t == "equations":
        return SignalEquations.from_dict(d)
    if t == "spectral":
        return SpectralModel.from_dict(d)
    raise ValueError(f"unknown method type {t}")


def determine(project, d: dict, ids: list[int], refit: bool = False):
    """Apply a method/model definition to spectra ``ids``.

    Returns (compounds, found rows, method). The method is calibrated on its
    ``calibration_ids`` when it carries no calibration — or always with
    ``refit`` (e.g. after calibration concentrations were corrected)."""
    res = project.resolver()
    m = method_from_definition(d)
    spectra = project.spectra(ids)
    cal = d.get("calibration_ids") or []
    if refit and not cal:
        raise ValueError("this method does not record its calibration spectra, so it cannot "
                         "be refitted — untick 'Refit'")
    if d["type"] == "univariate":
        if m.regression is None or refit:
            m.calibrate(project.spectra(cal), res)
        return [m.compound], [[m.predict(s, res)] for s in spectra], m
    if d["type"] == "equations":
        if m.K is None or refit:
            m.fit(project.spectra(cal), res)
        return m.compounds, m.predict(spectra, res).tolist(), m
    if d["type"] in uv.PROGRESSIVE_TYPES:
        if not m.is_fitted or refit:
            m.fit_spectra(project.spectra(cal), res)
        return m.compounds, m.predict_spectra(spectra, res).tolist(), m
    if not m.is_fitted or refit:
        m.fit(project.spectra(cal), res)
    return m.compounds, m.predict(spectra, res).tolist(), m


def definition_of(m, calibration_ids=None) -> dict:
    d = m.to_dict()
    if calibration_ids:
        d["calibration_ids"] = list(calibration_ids)
    return d


# --------------------------------------------------------------------------- #
# Confidence intervals and taken values stored with a result
# --------------------------------------------------------------------------- #
def _regression(d: dict | None, compound: str) -> Regression | None:
    """The calibration line that turns a signal into ``compound``'s
    concentration (univariate and progressive methods); None otherwise."""
    if not isinstance(d, dict):
        return None
    try:
        if d.get("type") == "univariate" and d.get("regression") and \
                d.get("compound") == compound and not d.get("direct"):
            return Regression.from_dict(d["regression"])
        if d.get("type") in uv.PROGRESSIVE_TYPES and (d.get("regressions") or {}).get(compound):
            return Regression.from_dict(d["regressions"][compound])
    except (TypeError, KeyError, ValueError):
        return None
    return None


def result_definition(project, data: dict, method_id: int | None = None) -> dict | None:
    for key in ("method", "model"):
        d = data.get(key)
        if isinstance(d, dict) and d.get("type") in TYPES:
            return d
    if method_id:
        try:
            return project.method(int(method_id))["definition"]
        except KeyError:
            return None
    return None


def compounds_of(data: dict) -> list[str] | None:
    found = data.get("found")
    if isinstance(found, list) and found and isinstance(found[0], dict):
        return [data.get("X", "X"), data.get("Y", "Y")]
    if data.get("compounds"):
        return list(data["compounds"])
    m = data.get("method")
    if isinstance(m, dict) and m.get("compound"):
        return [m["compound"]]
    return None


def found_rows(data: dict) -> list[list[Any]] | None:
    found = data.get("found")
    comps = compounds_of(data)
    if not (isinstance(found, list) and found and comps):
        return None
    if isinstance(found[0], dict):
        return [[f.get(c) for c in comps] for f in found]
    if isinstance(found[0], list):
        return found
    return None


def confidence_intervals(definition: dict | None, comps: list[str],
                         found: list[list[Any]]) -> list[list[float | None]] | None:
    """± half-widths of the 95 % CI of each found value (calibration line),
    or None when the method has no single calibration line per compound."""
    regs = [_regression(definition, c) for c in comps]
    if not any(regs):
        return None
    out = []
    for row in found:
        out.append([None if (r is None or v is None) else
                    (lambda h: None if math.isnan(h) else h)(r.ci_at_x(float(v)))
                    for r, v in zip(regs, row)])
    return out


def enrich(project, data: dict, method_id: int | None = None) -> dict:
    """Add to a result the values it must keep: the concentrations taken at
    the time of saving (a later correction must not silently change a stored
    recovery) and the 95 % CI of each found value."""
    ids, comps, found = data.get("ids"), compounds_of(data), found_rows(data)
    if not (ids and comps and found):
        return data
    data = dict(data)
    if "taken_snapshot" not in data:
        snap = {}
        for sid in ids:
            try:
                conc = project.record(int(sid)).concentrations
            except (KeyError, TypeError, ValueError):
                continue
            snap[str(int(sid))] = {c: conc.get(c) for c in comps}
        data["taken_snapshot"] = snap
    if "ci95" not in data:
        ci = confidence_intervals(result_definition(project, data, method_id), comps, found)
        if ci is not None:
            data["ci95"] = ci
    return data


def stale_taken(project, data: dict) -> list[tuple[int, str, Any, Any]]:
    """(spectrum, compound, taken when saved, taken now) where they differ."""
    out = []
    for sid, saved in (data.get("taken_snapshot") or {}).items():
        try:
            now = project.record(int(sid)).concentrations
        except KeyError:
            continue
        for c, v in saved.items():
            if now.get(c) != v:
                out.append((int(sid), c, v, now.get(c)))
    return out


# --------------------------------------------------------------------------- #
# Recalculation
# --------------------------------------------------------------------------- #
def recalculated_data(project, data: dict, definition: dict, refit: bool,
                      ids: list[int] | None = None) -> tuple[dict, Any]:
    """A copy of result ``data`` with found values recomputed by
    ``definition`` (refitted if asked) on the current data."""
    ids = [int(i) for i in (ids if ids is not None else data["ids"])]
    comps, found, m = determine(project, definition, ids, refit=refit)
    new = copy.deepcopy(data)
    for k in _DROP_ON_RECALC:
        new.pop(k, None)
    new.update({"ids": ids, "compounds": comps, "found": found})
    if new.get("enrichment_added") and len(comps) == 1:
        new["found_minus_added"] = [f[0] - new["enrichment_added"] for f in found]
    key = "method" if isinstance(new.get("method"), dict) else "model"
    new[key] = definition_of(m, definition.get("calibration_ids"))
    if definition.get("calibration_ids"):
        new["calibration_ids"] = list(definition["calibration_ids"])
    keep = set(ids)
    if new.get("excluded_ids"):
        new["excluded_ids"] = [i for i in new["excluded_ids"] if int(i) in keep]
    return new, m


def remap_ids(obj: Any, mapping: dict[int, int]) -> Any:
    """Deep copy of a method definition / result data with spectrum ids
    replaced (after re-processing created new versions of spectra)."""
    from spectro.storage.modelfile import remap_refs

    obj = copy.deepcopy(obj)

    def walk(o):
        if isinstance(o, dict):
            for k, v in o.items():
                if k in ("calibration_ids", "ids", "excluded_ids", "inputs") and \
                        isinstance(v, list):
                    o[k] = [mapping.get(int(x), int(x)) if isinstance(x, (int, float))
                            and not isinstance(x, bool) else x for x in v]
                elif k == "taken_snapshot" and isinstance(v, dict):
                    o[k] = {str(mapping.get(int(s), int(s))): c for s, c in v.items()}
                else:
                    walk(v)
                    if isinstance(v, dict) and v.get("type") in TYPES:
                        remap_refs(v, mapping)
        elif isinstance(o, list):
            for v in o:
                walk(v)
    walk(obj)
    if isinstance(obj, dict) and obj.get("type") in TYPES:
        remap_refs(obj, mapping)
    return obj


def recalculate_affected(project, ids, reason: str, mapping: dict[int, int] | None = None,
                         progress=None) -> dict:
    """After a correction to spectra ``ids`` (concentrations, or processing
    with ``mapping`` old → new id), bring every saved method and result that
    uses them up to date: methods are remapped and refitted on their
    calibration spectra (new versions), results are recalculated with the
    current version of their method (or their own model, refitted).

    Returns {"methods": {old: new}, "results": {old: new}, "skipped": [...]}."""
    mapping = {int(k): int(v) for k, v in (mapping or {}).items()}
    use = project.usage(list(ids) + list(mapping))
    out = {"methods": {}, "results": {}, "skipped": []}
    for md in use["methods"]:
        d = remap_ids(md["definition"], mapping)
        cal = d.get("calibration_ids")
        if not cal:
            out["skipped"].append(f"method #{md['id']} {md['name']}: calibration spectra not "
                                  "recorded")
            continue
        try:
            m = method_from_definition(d)
            res = project.resolver()
            if d["type"] == "univariate":
                m.calibrate(project.spectra(cal), res)
            elif d["type"] in uv.PROGRESSIVE_TYPES:
                m.fit_spectra(project.spectra(cal), res)
            else:
                m.fit(project.spectra(cal), res)
            nd = definition_of(m, cal)
            for k in ("calibration_ids",):
                nd[k] = cal
            out["methods"][md["id"]] = project.revise_method(md["id"], md["name"], nd, reason)
        except Exception as exc:  # noqa: BLE001 - reported to the user
            out["skipped"].append(f"method #{md['id']} {md['name']}: {exc}")
        if progress:
            progress()
    for r in use["results"]:
        data = remap_ids(r["data"], mapping)
        if not (data.get("ids") and found_rows(data)):
            out["skipped"].append(f"result #{r['id']} {r['name']}: no per-spectrum values")
            continue
        mid = r["method_id"]
        cur = project.current_method_version(mid) if mid else None
        try:
            if cur is not None:
                d, refit, mid = cur["definition"], False, cur["id"]
            else:
                d = result_definition(project, data)
                if d is None:
                    out["skipped"].append(f"result #{r['id']} {r['name']}: no method or model "
                                          "to recalculate with")
                    continue
                d = dict(d)
                if data.get("calibration_ids") and not d.get("calibration_ids"):
                    d["calibration_ids"] = data["calibration_ids"]
                refit = bool(d.get("calibration_ids"))
            new, _ = recalculated_data(project, data, d, refit)
            new["recalculated"] = {"source": "recalculate all affected", "method_id": mid,
                                   "refit": refit}
            out["results"][r["id"]] = project.revise_result(r["id"], r["name"], new, reason,
                                                            mid)
        except Exception as exc:  # noqa: BLE001 - reported to the user
            out["skipped"].append(f"result #{r['id']} {r['name']}: {exc}")
        if progress:
            progress()
    return out
