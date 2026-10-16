# -*- coding: utf-8 -*-
"""
并行预下载 GEO series matrix 与平台注释到 geo_cache/。
动机：NCBI FTP 单连接实测 ~100 KB/s 且会间歇停顿，串行下 4 个队列要十几分钟且看不出进度。
     并行 6 连接实测可把几十 MB 压到 2~3 分钟，并且每个文件都有独立进度。
用法：python geo_dl.py
下完后 sig_validate.py 会直接命中缓存（它的 fetch() 优先读 geo_cache/）。
"""
import concurrent.futures as cf
import os
import time
import urllib.request

UA = {"User-Agent": "Mozilla/5.0 (research; local meta-analysis)"}
BASE = os.path.dirname(os.path.abspath(__file__))
CACHE = os.path.join(BASE, "cohort_hunt", "geo_cache")
os.makedirs(CACHE, exist_ok=True)

COHORTS = ["GSE4183", "GSE10714", "GSE15960", "GSE128435"]
GPLS = ["GPL570", "GPL4133"]


def series_matrix_url(gse):
    n = gse[3:]
    s = n[:-3] if len(n) > 3 else ""
    return ("https://ftp.ncbi.nlm.nih.gov/geo/series/GSE%snnn/%s/matrix/%s_series_matrix.txt.gz"
            % (s, gse, gse))


def annot_url(gpl):
    n = gpl[3:]
    s = n[:-3] if len(n) > 3 else ""
    return "https://ftp.ncbi.nlm.nih.gov/geo/platforms/GPL%snnn/%s/annot/%s.annot.gz" % (s, gpl, gpl)


def download(url, tag, tries=4):
    dst = os.path.join(CACHE, os.path.basename(url))
    if os.path.exists(dst) and os.path.getsize(dst) > 500:
        return tag, "cached", os.path.getsize(dst), 0.0
    for i in range(tries):
        t0 = time.time()
        try:
            req = urllib.request.Request(url, headers=UA)
            with urllib.request.urlopen(req, timeout=120) as r:
                total = int(r.headers.get("Content-Length") or 0)
                got = 0
                tmp = dst + ".part"
                with open(tmp, "wb") as f:
                    while True:
                        c = r.read(262144)
                        if not c:
                            break
                        f.write(c)
                        got += len(c)
                os.replace(tmp, dst)
            return tag, "ok", os.path.getsize(dst), time.time() - t0
        except Exception as e:
            if i == tries - 1:
                return tag, "FAIL(%s)" % type(e).__name__, 0, time.time() - t0
            time.sleep(2 * (i + 1))
    return tag, "FAIL", 0, 0.0


def main():
    jobs = []
    for g in COHORTS:
        jobs.append((series_matrix_url(g), "matrix " + g))
    for g in GPLS:
        jobs.append((annot_url(g), "annot  " + g))

    print("并行下载 %d 个文件 (6 连接) ..." % len(jobs), flush=True)
    t0 = time.time()
    with cf.ThreadPoolExecutor(max_workers=6) as ex:
        futs = [ex.submit(download, u, t) for u, t in jobs]
        for f in cf.as_completed(futs):
            tag, st, sz, dt = f.result()
            print("  %-14s %-14s %9.2f MB  %5.1fs  累计 %.0fs" % (
                tag, st, sz / 1e6, dt, time.time() - t0), flush=True)
    print("\n全部完成，用时 %.0f 秒" % (time.time() - t0), flush=True)
    print("缓存目录：%s" % CACHE, flush=True)
    for fn in sorted(os.listdir(CACHE)):
        p = os.path.join(CACHE, fn)
        print("   %-48s %9.2f MB" % (fn, os.path.getsize(p) / 1e6), flush=True)


if __name__ == "__main__":
    main()
