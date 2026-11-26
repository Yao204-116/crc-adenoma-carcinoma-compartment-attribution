# -*- coding: utf-8 -*-
"""GSE8671 下载 v5：family.soft.gz（soft/ 卷，健康）HTTP Range 断点续传。
matrix 卷的 series matrix 服务端损坏（82,197 B 后 gzip 流坏死，FTP/HTTP 双协议证实）；
soft 卷前缀解压通过 → 改下 42MB family SOFT，验证脚本改用 SOFT 解析器。
"""
import gzip, os, time, urllib.request

URL = "http://ftp.ncbi.nlm.nih.gov/geo/series/GSE8nnn/GSE8671/soft/GSE8671_family.soft.gz"
EXPECT = 42_051_249
CACHE = os.path.join("cohort_hunt", "geo_cache", "GSE8671_family.soft.gz")
PART = CACHE + ".part"

buf = bytearray()
if os.path.exists(PART):
    with open(PART, "rb") as f:
        buf.extend(f.read())
print("起始已有 %d / %d B" % (len(buf), EXPECT), flush=True)

MAX_ROUNDS = 80
for rd in range(1, MAX_ROUNDS + 1):
    if len(buf) >= EXPECT:
        break
    got0 = len(buf)
    try:
        req = urllib.request.Request(URL, headers={
            "User-Agent": "Mozilla/5.0 (Windows NT 10.0; Win64; x64)",
            "Range": "bytes=%d-" % len(buf)})
        with urllib.request.urlopen(req, timeout=120) as r:
            print("round %d: status %s" % (rd, r.status), flush=True)
            while True:
                chunk = r.read(1 << 16)
                if not chunk:
                    break
                buf.extend(chunk)
                if len(buf) % (1 << 21) < (1 << 16):
                    print("  … %.1f MB" % (len(buf) / 1048576), flush=True)
    except Exception as e:
        print("round %d: +%.2f MB 后异常 %r" % (rd, (len(buf) - got0) / 1048576, e), flush=True)
    with open(PART, "wb") as f:
        f.write(buf)
    print("round %d: 累计 %.2f / %.2f MB" % (rd, len(buf) / 1048576, EXPECT / 1048576), flush=True)
    if len(buf) < EXPECT:
        time.sleep(2)

assert len(buf) >= EXPECT, "下载不完整: %d / %d" % (len(buf), EXPECT)
with open(PART, "wb") as f:
    f.write(buf)
print("全文件 CRC 校验…", flush=True)
with gzip.open(PART, "rb") as g:
    data = g.read()
print("解压成功 %.1f MB，尾部: %s" % (len(data) / 1048576, data[-70:]), flush=True)
os.replace(PART, CACHE)
print("完成 → %s (%d B)" % (CACHE, os.path.getsize(CACHE)), flush=True)
