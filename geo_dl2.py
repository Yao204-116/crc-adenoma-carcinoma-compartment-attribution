# -*- coding: utf-8 -*-
"""
带断点续传的下载器（专治 NCBI FTP 慢且会中途停顿）。
动机：实测同一文件速度在 100 KB/s ~ 7 KB/s 之间剧烈波动，且会把连接"挂住"在低速状态。
     一次下不完就白费，所以必须支持 Range 续传，逐段累积。
用法：python geo_dl2.py GSE4183 GSE128435
"""
import os
import sys
import time
import urllib.request

UA = {"User-Agent": "Mozilla/5.0 (research; local meta-analysis)"}
BASE = os.path.dirname(os.path.abspath(__file__))
CACHE = os.path.join(BASE, "cohort_hunt", "geo_cache")
os.makedirs(CACHE, exist_ok=True)


def url_of(gse):
    """
    支持三种目标：
      GSE#####      → series matrix
      SOFT:GPL##### → 平台 SOFT 文件（!platform_table 里含 GENE_SYMBOL）
      URL:<完整地址> → 任意文件（如 NCBI 的 Homo_sapiens.gene_info.gz）

    重要实测：GSE117606 所在平台 GPL25373 的行 ID 是 `<EntrezID>_at` 形式
    （对照样本：100009676_at / 10000_at / 10001_at），**不是 Agilent 探针名**。
    所以不需要 55.8 MB 的平台 SOFT，只要 Entrez→Symbol 映射即可
    （NCBI 的 Homo_sapiens.gene_info.gz 仅 5.2 MB）。
    另注：GPL17586 (HTA 2.0) 的 soft 有 7.6 GB → 遇到就换队列。
    """
    if gse.startswith("URL:"):
        return gse[4:]
    if gse.startswith("SOFT:"):
        gpl = gse[5:]
        n = gpl[3:]
        s = n[:-3] if len(n) > 3 else ""
        return ("https://ftp.ncbi.nlm.nih.gov/geo/platforms/GPL%snnn/%s/soft/%s_family.soft.gz"
                % (s, gpl, gpl))
    n = gse[3:]
    s = n[:-3] if len(n) > 3 else ""
    return ("https://ftp.ncbi.nlm.nih.gov/geo/series/GSE%snnn/%s/matrix/%s_series_matrix.txt.gz"
            % (s, gse, gse))


def total_size(url):
    try:
        r = urllib.request.urlopen(urllib.request.Request(url, headers=UA, method="HEAD"), timeout=30)
        return int(r.headers.get("Content-Length") or 0)
    except Exception:
        return 0


def grab(gse, rounds=40, per_round=90):
    url = url_of(gse)
    dst = os.path.join(CACHE, os.path.basename(url))
    part = dst + ".part"
    if os.path.exists(dst) and os.path.getsize(dst) > 500:
        print("%-10s 已存在，跳过 (%.2f MB)" % (gse, os.path.getsize(dst) / 1e6), flush=True)
        return True
    total = total_size(url)
    print("%-10s 目标大小 %.2f MB" % (gse, total / 1e6), flush=True)
    last = -1
    for rnd in range(rounds):
        have = os.path.getsize(part) if os.path.exists(part) else 0
        if total and have >= total:
            os.replace(part, dst)
            print("%-10s 完成 %.2f MB" % (gse, have / 1e6), flush=True)
            return True
        hdr = dict(UA)
        if have:
            hdr["Range"] = "bytes=%d-" % have
        t0 = time.time()
        got = 0
        try:
            with urllib.request.urlopen(urllib.request.Request(url, headers=hdr), timeout=per_round) as r:
                if have and r.status == 206:
                    mode = "ab"
                else:
                    have, mode = 0, "wb"     # 服务器不支持 Range → 从头来
                with open(part, mode) as f:
                    while True:
                        c = r.read(262144)
                        if not c:
                            break
                        f.write(c)
                        got += len(c)
                        if time.time() - t0 > per_round:
                            break            # 主动断开，进入下一轮续传
        except Exception as e:
            print("     [轮 %d] 中断 %s，已获取 %.2f MB" % (rnd + 1, type(e).__name__, got / 1e6), flush=True)
        now = os.path.getsize(part) if os.path.exists(part) else 0
        dt = max(time.time() - t0, 0.01)
        print("     [轮 %2d] +%.2f MB 累计 %.2f/%.2f MB  (%.0f KB/s)"
              % (rnd + 1, got / 1e6, now / 1e6, total / 1e6, got / 1024 / dt), flush=True)
        if now == last and got == 0:
            time.sleep(3)
        last = now
    print("%-10s 未完成（已 %.2f MB）" % (gse, (os.path.getsize(part) if os.path.exists(part) else 0) / 1e6), flush=True)
    return False


if __name__ == "__main__":
    todo = sys.argv[1:] or ["GSE4183", "GSE128435"]
    for g in todo:
        grab(g)
    print("\n缓存目录内容：", flush=True)
    for fn in sorted(os.listdir(CACHE)):
        print("   %-50s %9.2f MB" % (fn, os.path.getsize(os.path.join(CACHE, fn)) / 1e6), flush=True)
