"""Union per-shard queue progress into the one file the summary reads.

A sharded refresh runs N workers, each with its own progress file; nothing else
needs to know about the split. This merges them back into the base path so a
status check (or the refresh summary) reads the whole run from one place.

    merge_progress.py --base cache/refresh-progress.json --shards 6
"""
import argparse
import json
import os
import sys


def merge(base, shards):
    done, failed = set(), {}
    for k in range(shards):
        path = base.replace(".json", f"-{k}.json")
        try:
            with open(path, encoding="utf-8") as f:
                d = json.load(f)
        except Exception:
            continue
        done.update(d.get("done") or [])
        for g, e in (d.get("failed") or {}).items():
            failed[g] = e
    # Genes that succeeded in one shard and failed in another are successes.
    for g in done:
        failed.pop(g, None)
    return sorted(done), failed


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--base", required=True)
    ap.add_argument("--shards", type=int, required=True)
    a = ap.parse_args(argv)
    done, failed = merge(a.base, a.shards)
    tmp = f"{a.base}.{os.getpid()}.tmp"
    with open(tmp, "w", encoding="utf-8") as f:
        json.dump({"done": done, "failed": failed}, f, indent=1, sort_keys=True)
    os.replace(tmp, a.base)
    print(f"{len(done)} done, {len(failed)} failed from {a.shards} shards")
    return 0


if __name__ == "__main__":
    sys.exit(main())
