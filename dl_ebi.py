#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""EBI(ftp.ebi.ac.uk) 断点续传下载 + zip 完整性校验。
EBI 与 NCBI 走不同网络路径，实测可通；且支持 Range，可真正续传。
用法: python dl_ebi.py <url> <out.zip> [rounds]
"""
import os
import subprocess
import sys
import time
import zipfile

url, out = sys.argv[1], sys.argv[2]
rounds = int(sys.argv[3]) if len(sys.argv) > 3 else 60


def ok_zip(p):
    if not os.path.exists(p):
        return False
    try:
        with zipfile.ZipFile(p) as z:
            return z.testzip() is None
    except Exception:
        return False


for i in range(1, rounds + 1):
    if ok_zip(out):
        print("[OK] zip complete: %s (%.2f MB)" % (out, os.path.getsize(out) / 1e6), flush=True)
        break
    sz0 = os.path.getsize(out) if os.path.exists(out) else 0
    r = subprocess.run(["curl", "-sL", "--http1.1", "--max-time", "300",
                        "-C", "-", "-o", out, url], capture_output=True, text=True)
    sz1 = os.path.getsize(out) if os.path.exists(out) else 0
    print("[try %d] rc=%s  %d -> %d bytes  zip=%s"
          % (i, r.returncode, sz0, sz1, "OK" if ok_zip(out) else "partial"), flush=True)
    if ok_zip(out):
        break
    if sz1 == sz0 and i > 3:
        # 连续多次零增长：退回整包重下
        print("   no progress -> restart from scratch", flush=True)
        try:
            os.remove(out)
        except OSError:
            pass
    time.sleep(3)
else:
    print("[FAIL] incomplete after %d tries" % rounds, flush=True)
    sys.exit(1)
