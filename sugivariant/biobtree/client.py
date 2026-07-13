#!/usr/bin/env python3
"""Biobtree client — REST over HTTP/1.1 with a urllib3 connection pool.

Three primitives + four pure parsers + a per-call reproducibility log:
  search(term, source=None)       -> /api/search dict
  entry(identifier, source)       -> /api/entry dict
  bbmap(ids, chain, page=None)    -> /api/map dict
  rows(resp)                      -> dicts keyed by REST schema cols
  map_targets(resp)               -> flat list of target dicts (escape-aware)
  map_all(ids, chain, cap=60)     -> paginated map, deduped, loop-safe
  xref_counts(entry_resp)         -> {dataset: count} from entry's xref table
  CALLS                           -> per-call reproducibility log
"""
import html
import json
import os
import re
import time

import urllib3

# Two transports, identical response shapes (verified byte-identical for
# map/search/entry):
#   "direct" — Atlas → biobtree Go server `/ws/` endpoints, `mode=lite` injected.
#              ~10x faster for bulk page-gen (no uvicorn hop, no double-JSON
#              encode/decode, no not_found N+1 probes). The default.
#   "gate"   — Atlas → the FastAPI/MCP gate `/api/` (it forwards to Go with
#              mode=lite and adds LLM-UX niceties Atlas doesn't need). Kept as a
#              fallback for environments that only expose the gate.
# Both overridable via env (README documents ATLAS_BIOBTREE).
API_BASE = os.environ.get("ATLAS_BIOBTREE", "http://127.0.0.1:9291").rstrip("/")
TRANSPORT = os.environ.get("ATLAS_BIOBTREE_TRANSPORT", "direct").strip().lower()

# op → (direct path, gate path, inject mode=lite?). entry is hardwired to
# entryLite in Go, so it takes no mode param on either transport.
_OPS = {
    "search": ("/ws/",      "/api/search", True),
    "map":    ("/ws/map/",  "/api/map",    True),
    "entry":  ("/ws/entry/", "/api/entry", False),
}

# Back-compat: some tooling/tests read `API`. Expose the resolved base+prefix.
API = API_BASE + ("" if TRANSPORT == "direct" else "/api")
CALLS = []


class BiobtreeError(RuntimeError):
    """Biobtree REST returned an HTTP error or an unparseable body. A typed
    error so drivers can skip+log an entity rather than crash on a raw
    JSONDecodeError / urllib3 exception."""


# One pool, HTTP keep-alive. Pool-level retries off — we retry in _get() so we
# can distinguish retryable (5xx/429/timeout/truncated body) from terminal
# (4xx) and log. maxsize bounds concurrent connections.
# NOTE: CALLS is a module global; the client is safe for PROCESS-level
# parallelism (each process has its own CALLS) but NOT thread-level.
_POOL = urllib3.PoolManager(
    num_pools=2,
    maxsize=16,
    retries=False,
    # read=120s (was 15): the genome's heaviest genes (TTN, BRCA2) return large
    # map payloads that, under full-regen contention (30 workers) plus the added
    # HPA + variant-annotation calls, exceeded 15s and got skipped — dropping two
    # flagship genes from the v1.2.0 corpus. A generous 2-min read ceiling lets
    # even the biggest payloads finish under load; the 4-retry backoff still
    # bounds a genuinely-stuck request.
    timeout=urllib3.Timeout(connect=5.0, read=120.0),
    block=False,
)

_RETRY_STATUS = {429, 500, 502, 503, 504}
_MAX_ATTEMPTS = 4


_ENTITY_CI = re.compile(r"&(amp|lt|gt|quot|apos|nbsp);", re.I)
_ENTITY_MAP = {"amp": "&", "lt": "<", "gt": ">", "quot": '"', "apos": "'", "nbsp": " "}


def _unescape(obj):
    """Recursively HTML-unescape every string in a biobtree response, once, at
    the single point all responses flow through — so EVERY consumer (map
    targets, search rows, entry attributes, UniProt CC narratives, GO terms)
    renders clean text (3'-end, not 3&apos;-end). Entities never contain a raw
    '|', so this is safe for the pipe-delimited map/search values.

    html.unescape is case-sensitive (only `&amp;`), but source text — notably
    ClinicalTrials.gov titles — carries mixed-case malformed entities like
    `&Amp;`. Normalize the common five case-insensitively first."""
    if isinstance(obj, str):
        return html.unescape(_ENTITY_CI.sub(lambda m: _ENTITY_MAP[m.group(1).lower()], obj))
    if isinstance(obj, list):
        return [_unescape(x) for x in obj]
    if isinstance(obj, dict):
        return {k: _unescape(v) for k, v in obj.items()}
    return obj


def _get(op: str, params: dict) -> dict:
    """GET with bounded exponential-backoff retry on transient failures.

    `op` is the logical operation ("search"/"map"/"entry") — the URL path and
    the `mode=lite` injection are resolved per transport. Raises BiobtreeError
    on a 4xx, on biobtree's inline {"Err": …} body (the direct Go transport
    returns query errors with a 200 status), or after exhausting retries on
    5xx/429/timeout/truncated-body. Guards json.loads so an HTML error page or a
    short read can't surface as a bare JSONDecodeError mid-collection."""
    direct_path, gate_path, inject_mode = _OPS[op]
    if TRANSPORT == "direct":
        url = API_BASE + direct_path
        if inject_mode:
            params = {**params, "mode": "lite"}
    else:
        url = API_BASE + gate_path
    last = None
    for attempt in range(_MAX_ATTEMPTS):
        try:
            r = _POOL.request("GET", url, fields=params)
        except urllib3.exceptions.HTTPError as e:        # connect/read timeout, reset
            last = e
        else:
            if r.status >= 400 and r.status not in _RETRY_STATUS:
                raise BiobtreeError(f"{op} → HTTP {r.status}: {bytes(r.data[:200])!r}")
            if r.status in _RETRY_STATUS:
                last = BiobtreeError(f"{op} → HTTP {r.status}")
            else:
                try:
                    body = json.loads(r.data)
                except (json.JSONDecodeError, ValueError) as e:
                    last = e                              # truncated / non-JSON → retry
                else:
                    # biobtree signals a query-level error (e.g. unknown dataset)
                    # with an {"Err": …} body — the gate turns this into a 503,
                    # but the direct Go server returns it inline with HTTP 200.
                    # Terminal: a bad chain won't fix itself on retry.
                    if isinstance(body, dict) and body.get("Err"):
                        raise BiobtreeError(f"{op} → {body['Err']}")
                    CALLS.append({"path": op, "params": params})
                    return _unescape(body)
        if attempt < _MAX_ATTEMPTS - 1:
            time.sleep(0.5 * (2 ** attempt))              # 0.5s, 1s, 2s
    raise BiobtreeError(f"{op} failed after {_MAX_ATTEMPTS} attempts: {last}")


def search(term: str, source: str = None) -> dict:
    p = {"i": term}
    if source:
        p["s"] = source
    return _get("search", p)


def entry(identifier: str, source: str) -> dict:
    return _get("entry", {"i": identifier, "s": source})


def bbmap(ids: str, chain: str, page: str = None) -> dict:
    params = {"i": ids, "m": chain}
    if page:
        # cursor param is `p` (NOT `page`; FastAPI silently drops unknown
        # query params so the wrong name fails silently).
        params["p"] = page
    return _get("map", params)


def rows(resp: dict) -> list:
    """search/data rows -> list of dicts keyed by schema columns.

    biobtree returns `"data": null` (literal None, not the missing key) when
    a search has zero results — coerce defensively."""
    cols = (resp.get("schema") or "").split("|")
    out = []
    for r in (resp.get("data") or []):
        parts = r.split("|")
        # Length mismatch = an unescaped pipe shifted the columns; zipping it
        # would misalign every later value (the MeSH-row bug class). Skip it.
        if len(parts) != len(cols):
            continue
        out.append(dict(zip(cols, parts)))     # strings already unescaped in _get
    return out


def map_targets(resp: dict) -> list:
    """map response -> flat list of target dicts keyed by schema columns.

    Schema-aware split: some target fields embed an escaped pipe (\\|),
    e.g. bgee_evidence ids ENSG..\\|UBERON:.. — split on real separators only.
    """
    cols = resp.get("schema", "").split("|")
    out = []
    for m in (resp.get("mappings") or []):
        for t in (m.get("targets") or []):
            parts = [p.replace("\x00", "|") for p in t.replace("\\|", "\x00").split("|")]
            if len(parts) != len(cols):     # column-shift guard (see rows())
                continue
            out.append(dict(zip(cols, parts)))   # strings already unescaped in _get
    return out


def map_all(ids, chain, cap=60):
    """All target dicts across pages, deduped. Pagination via the `p=` cursor.
    Dedupe + stop-on-no-new make it loop-safe; `cap` bounds total pages
    (cap*100 rows is the worst case)."""
    out, page, n, seen = [], None, 0, set()
    while True:
        resp = bbmap(ids, chain, page)
        new = 0
        for t in map_targets(resp):
            key = t.get("id") or tuple(t.values())
            if key in seen:
                continue
            seen.add(key)
            out.append(t)
            new += 1
        pg = resp.get("pagination", {}) or {}
        nxt = pg.get("next_token")
        if not pg.get("has_next") or n >= cap or not nxt or new == 0:
            break
        page, n = nxt, n + 1
    return out


def xref_counts(entry_resp):
    """hgnc/ensembl/transcript `entry` carries a dataset|count xref table —
    exact totals for clinvar/spliceai/msigdb/hpo/gwas/collectri/... in 1 call."""
    out = {}
    for r in (entry_resp.get("xrefs", {}) or {}).get("data", []) or []:
        parts = r.split("|")
        if len(parts) < 2:
            continue
        try:
            out[parts[0]] = int(parts[1])
        except ValueError:                  # malformed/empty count → skip, don't crash
            continue
    return out
