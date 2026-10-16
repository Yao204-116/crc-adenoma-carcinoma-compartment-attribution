# -*- coding: utf-8 -*-
"""GSE161277 matrix 最终下载器：.part 隔离 + 完整性校验 + rename。

规则：
  - 目标文件 size == HEAD total 且 gzip 头合法 → skip
  - 下载写入 {name}.part（本进程独占），支持从 .part 尺寸续传（校验 206）
  - size == total 且 gzip 头为 \x1f\x8b → rename 原子落盘
用法：python sc161277_dl_final.py 0   # 0/1 = 两进程分片
"""
import os, re, sys, time, urllib.request

UA = {"User-Agent": "Mozilla/5.0 (research; local meta-analysis)"}
BASE = os.path.dirname(os.path.abspath(__file__))
D = os.path.join(BASE, "sc161277")
TXT = os.path.join(BASE, "cohort_hunt2", "gsm_cache", "GSE161277.txt")


def urls():
    txt = open(TXT, encoding="utf-8", errors="replace").read()
    us = sorted(set(u.replace("ftp://", "https://") for u in
                    re.findall(r"!Sample_supplementary_file_\d+ = (ftp://\S+matrix\.mtx\.gz)", txt)))
    return us


def head_total(u):
    for _ in range(5):
        try:
            r = urllib.request.urlopen(urllib.request.Request(u, headers=UA, method="HEAD"), timeout=30)
            return int(r.headers.get("Content-Length") or 0)
        except Exception:
            time.sleep(2)
    return 0


def gzip_head_ok(p):
    try:
        with open(p, "rb") as f:
            return f.read(2) == b"\x1f\x8b"
    except OSError:
        return False


def download(u, dst, part):
    total = head_total(u)
    if os.path.exists(dst) and total and os.path.getsize(dst) == total and gzip_head_ok(dst):
        return "SKIP"
    for rd in range(60):
        have = os.path.getsize(part) if os.path.exists(part) else 0
        if total and have >= total:
            break
        hdr = dict(UA)
        if have:
            hdr["Range"] = "bytes=%d-" % have
        try:
            resp = urllib.request.urlopen(urllib.request.Request(u, headers=hdr), timeout=120)
            mode = "ab"
            if have and (resp.status != 206 or
                         not (resp.headers.get("Content-Range") or "").startswith("bytes %d-" % have)):
                mode = "wb"  # 服务器不支持续传 → 重头
                have = 0
            with resp, open(part, mode) as f:
                while True:
                    chunk = resp.read(1 << 20)
                    if not chunk:
                        break
                    f.write(chunk)
        except Exception:
            time.sleep(2)
        have = os.path.getsize(part) if os.path.exists(part) else 0
        if total and have >= total:
            break
    if total and os.path.exists(part) and os.path.getsize(part) == total and gzip_head_ok(part):
        os.replace(part, dst)
        return "DONE"
    return "FAIL(%s/%s)" % (os.path.getsize(part) if os.path.exists(part) else 0, total)


if __name__ == "__main__":
    shard = int(sys.argv[1]) if len(sys.argv) > 1 else 0
    us = [u for i, u in enumerate(urls()) if i % 2 == shard]
    print("shard %d: %d files" % (shard, len(us)))
    for u in us:
        name = u.rsplit("/", 1)[-1]
        t0 = time.time()
        st = download(u, os.path.join(D, name), os.path.join(D, name + ".p%d" % shard))
        print("%-6s %-58s %.0fs" % (st, name, time.time() - t0), flush=True)
    print("shard %d finished" % shard)
