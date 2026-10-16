#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""持久重试下载器：多目标、长间隔重试，直到 gzip 完整。
网络中断（curl rc=35 代理断流）时不会放弃，适合后台挂机。
用法: python dl_persist.py <manifest.tsv> [sleep_sec] [rounds]
manifest.tsv 每行: <url>\t<outpath>
"""
import gzip
import os
import subprocess
import sys
import time

targets = []
with open(sys.argv[1], encoding="utf-8") as f:
    for ln in f:
        parts = ln.strip().lstrip("\ufeff").split()
        # 只接受 "http(s)://...  输出路径" 形式，注释与表头自动跳过
        if len(parts) >= 2 and parts[0].startswith("http"):
            targets.append([parts[0], parts[1]])

sleep_s = int(sys.argv[2]) if len(sys.argv) > 2 else 60
rounds = int(sys.argv[3]) if len(sys.argv) > 3 else 240


def complete(p):
    if not os.path.exists(p):
        return False
    try:
        with gzip.open(p, "rb") as fh:
            while fh.read(1 << 22):
                pass
        return True
    except Exception:
        return False


done = set()
for rnd in range(1, rounds + 1):
    pending = [(u, o) for u, o in targets if o not in done and not complete(o)]
    if not pending:
        print("[ALL DONE] round %d" % rnd, flush=True)
        break
    for url, out in pending:
        os.makedirs(os.path.dirname(out) or ".", exist_ok=True)
        rc = subprocess.run(
            ["curl", "-sL", "--http1.1", "--max-time", "900", "-o", out, url],
            capture_output=True, text=True).returncode
        sz = os.path.getsize(out) if os.path.exists(out) else 0
        ok = complete(out)
        print("[r%d] %s rc=%s bytes=%d %s" % (rnd, os.path.basename(out), rc, sz,
                                              "OK" if ok else "partial/FAIL"), flush=True)
        if ok:
            done.add(out)
    if len(done) == len(targets):
        print("[ALL DONE]", flush=True)
        break
    time.sleep(sleep_s)
else:
    print("[EXIT] rounds exhausted; incomplete:", [o for _, o in targets if o not in done], flush=True)
