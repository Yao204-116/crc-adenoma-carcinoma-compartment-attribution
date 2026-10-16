#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""整包重试下载器：反复完整 GET，直到 gzip 流校验通过为止。
代理隧道不支持 Range/HEAD，只能整包重试。
用法: python dl_loop.py <url> <outpath> [max_attempts]
"""
import gzip
import os
import subprocess
import sys
import time

url = sys.argv[1]
out = sys.argv[2]
max_attempts = int(sys.argv[3]) if len(sys.argv) > 3 else 20


def complete(path):
    """gzip 流能读到 EOF 才算完整。"""
    if not os.path.exists(path):
        return False
    try:
        with gzip.open(path, "rb") as f:
            while f.read(1 << 20):
                pass
        return True
    except Exception:
        return False


for i in range(1, max_attempts + 1):
    if complete(out):
        print("[OK] already complete: %s (%.2f MB)" % (out, os.path.getsize(out) / 1e6))
        break
    print("[try %d/%d] downloading ..." % (i, max_attempts), flush=True)
    r = subprocess.run(
        ["curl", "-sL", "--http1.1", "--max-time", "900", "-o", out, url],
        capture_output=True, text=True,
    )
    sz = os.path.getsize(out) if os.path.exists(out) else 0
    print("   curl rc=%s bytes=%d" % (r.returncode, sz), flush=True)
    if complete(out):
        print("[OK] complete after %d attempt(s): %.2f MB" % (i, sz / 1e6))
        break
    time.sleep(2)
else:
    print("[FAIL] gave up after %d attempts" % max_attempts)
    sys.exit(1)
