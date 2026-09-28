"""Extract known concentrations from sample names such as 'PAR 10 CAF 1.3'."""

from __future__ import annotations

import re

_NUM = r"(\d+(?:[.,]\d+)?)"


def concentrations_from_name(name: str, compounds: list[str]) -> dict[str, float]:
    """'PAR 10 CAF 1.3' → {'PAR': 10.0, 'CAF': 1.3} (by compound name); otherwise,
    if the name holds at least one number per compound, numbers are assigned to
    the compounds in order ('mix_10_5' → first=10, second=5)."""
    found: dict[str, float] = {}
    for c in compounds:
        m = re.search(re.escape(c) + r"\s*[_\-:= ]?\s*" + _NUM, name, re.I)
        if m:
            found[c] = float(m.group(1).replace(",", "."))
    if found or not compounds:
        return found
    nums = re.findall(_NUM, name)
    if len(nums) >= len(compounds):
        return {c: float(n.replace(",", ".")) for c, n in zip(compounds, nums)}
    return {}
