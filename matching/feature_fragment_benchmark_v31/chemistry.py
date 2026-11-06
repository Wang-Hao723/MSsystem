"""Graded elemental-balance checks; these checks never create a fixed score."""

from __future__ import annotations

import re
from collections import Counter
from collections.abc import Mapping
from typing import Any


_TOKEN = re.compile(r"([A-Z][a-z]?)(\d*)")


def parse_formula(formula: str | None) -> dict[str, int] | None:
    text = str(formula or "").strip()
    if not text:
        return None
    position = 0
    counts: Counter[str] = Counter()
    for match in _TOKEN.finditer(text):
        if match.start() != position:
            return None
        element, raw_count = match.groups()
        count = int(raw_count) if raw_count else 1
        if count <= 0:
            return None
        counts[element] += count
        position = match.end()
    if position != len(text) or not counts:
        return None
    return dict(sorted(counts.items()))


def _result(status: str, reason: str, **details: Any) -> dict[str, Any]:
    return {
        "status": status,
        "reason": reason,
        "score_contribution": 0.0,
        "likelihood_ratio_modified": False,
        **details,
    }


def grade_product_formula(
    parent_formula: str | None,
    product_formula: str | None,
    *,
    ion_state_complete: bool,
    allowed_adduct_elements: Mapping[str, int] | None = None,
) -> dict[str, Any]:
    parent = parse_formula(parent_formula)
    product = parse_formula(product_formula)
    if parent is None or product is None:
        return _result("ambiguous", "formula_parse_incomplete")
    allowance = Counter(parent)
    for element, count in (allowed_adduct_elements or {}).items():
        if int(count) < 0:
            return _result("ambiguous", "negative_adduct_allowance_invalid")
        allowance[str(element)] += int(count)
    excess = {
        element: count - allowance.get(element, 0)
        for element, count in product.items()
        if count > allowance.get(element, 0)
    }
    if not ion_state_complete:
        return _result(
            "ambiguous",
            "adduct_charge_radical_or_isotope_state_incomplete",
            elemental_excess=excess,
        )
    if excess:
        return _result(
            "incompatible",
            "product_element_count_exceeds_parent_and_explicit_adduct_allowance",
            elemental_excess=excess,
        )
    return _result("compatible", "elemental_subformula_is_possible")


def grade_neutral_loss_balance(
    parent_formula: str | None,
    product_formula: str | None,
    neutral_loss_formula: str | None,
    *,
    ion_state_complete: bool,
) -> dict[str, Any]:
    parent = parse_formula(parent_formula)
    product = parse_formula(product_formula)
    loss = parse_formula(neutral_loss_formula)
    if parent is None or product is None or loss is None:
        return _result("ambiguous", "formula_parse_incomplete")
    if not ion_state_complete:
        return _result(
            "ambiguous", "adduct_charge_radical_or_isotope_state_incomplete"
        )
    elements = sorted(set(parent) | set(product) | set(loss))
    imbalance = {
        element: parent.get(element, 0) - product.get(element, 0) - loss.get(element, 0)
        for element in elements
        if parent.get(element, 0) != product.get(element, 0) + loss.get(element, 0)
    }
    if imbalance:
        return _result(
            "incompatible", "parent_product_loss_elemental_balance_fails", imbalance=imbalance
        )
    return _result("compatible", "parent_product_loss_elemental_balance_matches")
