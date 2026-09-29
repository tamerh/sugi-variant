"""Variant slug contract — MOVED to the shared `sugislug` package (single source of
truth; Sugi Atlas imports the same code to build resolvable variant URLs). This module
re-exports it UNCHANGED, so sugivariant's app, build and tests are untouched and the
slug output (live URLs + the SQLite resolution index) stays byte-identical.

See sugislug/varslug.py for the implementation and history.
"""
from sugislug.varslug import (  # noqa: F401  (re-export for backward compatibility)
    name_gene, parse_hgvs, variant_slugs,
    _norm, _norm_hgvs, _cap, _SLUG_MAX,
)
