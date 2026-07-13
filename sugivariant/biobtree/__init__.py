"""Shared biobtree access layer — single client, shared by gene/drug/disease.

Entity packages import from here for every biobtree call; no entity-specific
logic lives in this module. Adding a new entity = adding `src/atlas/<entity>/`
that calls into `atlas.biobtree.client`.
"""
from sugivariant.biobtree.client import (
    API, CALLS, BiobtreeError,
    search, entry, bbmap, rows, map_targets, map_all, xref_counts,
)

__all__ = ["API", "CALLS", "BiobtreeError", "search", "entry", "bbmap", "rows",
           "map_targets", "map_all", "xref_counts"]
