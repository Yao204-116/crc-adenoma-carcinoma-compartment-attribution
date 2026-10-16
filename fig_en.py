# -*- coding: utf-8 -*-
"""图表英文化工具 + 出图守卫（投稿版专用）。

★ 为什么不直接替换源码里的中文字符串？
  这些出图脚本里，**面板名/区域名列表同时充当「数据字典的键」和「图上显示的文字」**
  （`ORDER`、`panels`、`sets_`、`order` 都是这种双重身份）。把源码字面量整体换成英文，
  键会被一起翻掉，于是 `MAT.get(英文名)` 静默查空 ——
  本仓库实际踩过：fig5 墨迹覆盖率 47.4% → 7.85%，fig6 6.79% → 1.17%，
  **不报错、不崩溃，图只是变空白**。

  正确做法：**数据键一律保持中文，只在「文字上屏」那一行做映射。**

本模块只做三件事：
  1. `panel_en()` —— 显示名（面板 / 子面板 / 模式标签）的**唯一**英文映射，
                     图与补充表共用，避免两处漂移；
  2. `wrap()` / `declutter()` —— 按**实测渲染宽度**折行、把互压的标注错开；
  3. `finish()`   —— 出图收尾守卫：残留汉字/全角标点 → 报错退出；
                     文字横向溢出画布 → 告警；文字互相压叠 → 告警。
                     杜绝「漏翻一张图 / 把字裁掉 / 两块字叠在一起却没人发现」。

英文用词**逐字对齐 BMC 投稿稿的图注**（`Manuscript_BMC_Genomics_投稿版.html`），
不是自由翻译 —— 否则审稿人无法把图对回正文。稿件里的原话：
  "thirteen cell-programme panels"、"the Paneth-cell programme"、
  "the goblet·5, enterocyte-differentiation·8 and metabolic·5 panels"、
  "pure-epithelial-baseline·7 (negative control)"、
  "stroma/activated-CAF (cross-compartment control)"、
  "P = adenoma peak (adenoma ≥ normal and adenoma > carcinoma) · D monotonic decrease ·
   A already decreased at adenoma"。
"""
import os
import re

import matplotlib
import matplotlib.pyplot as plt
from matplotlib import text as mtext
from matplotlib.font_manager import FontProperties
from matplotlib.textpath import TextPath

# ---------------------------------------------------------------------------
# 面板名：中文数据键 → 英文显示名
# 键的清单来自 r3_facts.json / three_layer_pattern.json /
# gse117606_paired_results.json 的实际数据键，**不要凭记忆增删**。
# 未收录的名字原样返回（例如 WGCNA-turquoise 本身就是英文）。
# ---------------------------------------------------------------------------
PANELS_EN = {
    # --- 13 个策划面板（r3_facts.json panel_table）---
    "潘氏细胞·4": "Paneth-cell·4",
    "黏液杯状·5": "goblet·5",
    "上皮分化·核心8": "enterocyte-differentiation·8",
    "相代谢·5": "metabolic·5",
    "纯上皮基线·7": "pure-epithelial-baseline·7",
    "间质·激活/CAF": "stroma/activated-CAF",
    "上皮·分化/功能": "epithelial-differentiation",
    "上皮·基线/结构": "epithelial-structural",
    "泛间质·6": "stromal/activation·6",
    "髓系·激活": "myeloid-activation",
    "间质·基线": "stromal-baseline",
    "上调分量(210)": "up-component·210",
    "下调分量(103)": "down-component·103",
    # --- 配对分析里的三个分量行（gse117606_paired_results.json B_paired）---
    "全签名": "313-gene signature",
    "上调分量(间质/髓系激活)": "stromal/myeloid-activation component",
    "下调分量(上皮功能关闭)": "epithelial-function component",
    # --- fig5 热图里额外出现的行 ---
    "313签名": "signature·313",
    "Hub基因": "hub genes",
    # --- Paneth 敏感性分析的子面板（r3_sens_paneth.json 的键，S3b 表）---
    "全4基因(DEFA5/6+OLFM4+SPINK1)": "all four genes (DEFA5/6 + OLFM4 + SPINK1)",
    "纯潘氏(DEFA5+DEFA6)": "Paneth core only (DEFA5 + DEFA6)",
    "仅OLFM4(干细胞共享)": "OLFM4 only (shared with stem cells)",
    "仅SPINK1": "SPINK1 only",
    "去OLFM4(DEFA5/6+SPINK1)": "without OLFM4 (DEFA5/6 + SPINK1)",
    # --- 三层模式标签（three_layer_pattern.json 的 pattern 字段，S3b 表）---
    "P 腺瘤峰": "P (adenoma peak)",
    "D 单调降": "D (monotone decrease)",
    "A 已降": "A (already reduced)",
}


def panel_en(name):
    """把中文面板名翻成英文；未收录（或已是英文）则原样返回。"""
    return PANELS_EN.get(name.strip(), name)


# ---------------------------------------------------------------------------
# 标题折行：按**实测渲染宽度**断行
# ---------------------------------------------------------------------------
def text_width(s, fontsize, weight="bold"):
    """`s` 在 DejaVu Sans 下的渲染宽度（英寸，不含换行）。"""
    fp = FontProperties(family="DejaVu Sans", size=fontsize, weight=weight)
    return TextPath((0, 0), s, prop=fp).get_extents().width / 72.0


def wrap(s, maxw_in, fontsize, weight="bold"):
    """按词贪心折行，使每行渲染宽度 ≤ `maxw_in` 英寸；原有 `\\n` 视为硬断点。

    ★ 为什么需要：英文化后标题显著变长（一个汉字 ≈ 2 个拉丁字母宽）。
      `bbox_inches="tight"` 的图会被动撑宽画布 → 相对字号变小、留白变多；
      **不带 tight 的图（fig7/fig8/figS1/figS3）则把超出的字直接裁掉**（静默丢字）。
      用「估算字符数」不可靠（i/l 与 m/W 宽度差 3 倍），所以按实测宽度断。

    经验值（DejaVu Sans, 英寸/字符）：粗体 ≈ 0.0079 × fontsize。
    """
    out = []
    for para in s.split("\n"):
        cur = ""
        for wd in para.split(" "):
            trial = wd if not cur else cur + " " + wd
            if cur and text_width(trial, fontsize, weight) > maxw_in:
                out.append(cur)
                cur = wd
            else:
                cur = trial
        out.append(cur)
    return "\n".join(out)


# ---------------------------------------------------------------------------
# 出图守卫
# ---------------------------------------------------------------------------
# 汉字 / CJK 标点 / 全角形式 —— 这些字符 DejaVu Sans 里没有，渲染出来就是豆腐块。
BAD = re.compile(
    "[\u2e80-\u2eff\u3000-\u303f\u3400-\u4dbf\u4e00-\u9fff"
    "\uf900-\ufaff\ufe30-\ufe4f\uff00-\uffef]"
)


def _shown(t):
    """该文字是否**真的会被画出来**（自身与所有祖先都可见）。

    ★ 为什么必须看这个：`fig.findobj(Text)` 会把**已隐藏**的对象也捞出来。
      最典型的是**轴外刻度标签** —— XAxis._update_ticks() 对超出 view 的刻度只把
      `label.set_visible(False)`，对象仍在树里，`get_window_extent()` 依然返回一个
      **落在画布外**的坐标。不过滤就会报假警（fig8 的 '0.30'、'−0.05' 就是这样骗过守卫的）。
    """
    o = t
    while o is not None:
        if not getattr(o, "get_visible", lambda: True)():
            return False
        o = getattr(o, "get_parent", lambda: None)() if hasattr(o, "get_parent") else None
    return True


def _offview_tick_ids(fig):
    """返回「刻度值落在轴 view 区间之外」的刻度标签 id 集合。

    ★ 为什么 `get_visible()` 不够：matplotlib 并不隐藏越界刻度，而是在 `Axis.draw()`
      里**直接跳过不画**（`_update_ticks` → `ticks_to_draw`）。Text 对象依旧在画布树里，
      `get_window_extent()` 照旧返回一个**落在画布外**的坐标 ——
      于是溢出检查会报假警（fig8 的 '0.30'、'−0.05' 就是这么骗过守卫的：
      axes[1] 的 xlim 上限是 0.264，而 0.30 那个刻度的锚点被算到了 1400 px，画布只有 1350 px）。
    """
    out, seen = set(), set()
    for ax in fig.axes:
        for axis, idx in ((ax.xaxis, 0), (ax.yaxis, 1)):
            try:
                lo, hi = sorted(axis.get_view_interval())
            except Exception:
                continue
            for lab in axis.get_ticklabels():
                if id(lab) in seen:
                    continue
                seen.add(id(lab))
                try:
                    v = lab.get_position()[idx]
                except Exception:
                    continue
                if not (lo <= v <= hi):
                    out.add(id(lab))
    return out


def _collect(fig):
    """把画布上所有**会显示**的文字对象收齐（去重）。"""
    seen, out = set(), []
    offview = _offview_tick_ids(fig)

    def add(t):
        if t is None or id(t) in seen or not _shown(t) or id(t) in offview:
            return
        try:
            t.get_text()
        except Exception:
            return
        seen.add(id(t))
        out.append(t)

    for t in fig.findobj(mtext.Text):
        add(t)
    add(getattr(fig, "_suptitle", None))
    for ax in fig.axes:                       # 含 colorbar 的 axes
        add(ax.title)
        add(getattr(ax, "_left_title", None))
        add(getattr(ax, "_right_title", None))
        add(ax.xaxis.label)
        add(ax.yaxis.label)
        for t in ax.get_xticklabels():
            add(t)
        for t in ax.get_yticklabels():
            add(t)
        for t in ax.texts:
            add(t)
        lg = ax.get_legend()
        if lg is not None:
            for t in lg.get_texts():
                add(t)
    return out


def scan(fig):
    """返回画布上所有残留的中文/全角文字（list[str]）。"""
    fig.canvas.draw()                         # 刻度标签在这一步才生成
    return [t.get_text() for t in _collect(fig) if BAD.search(t.get_text())]


def scan_overflow(fig, tol_px=6):
    """返回横向超出画布的 (文字, 左边界, 右边界, 画布宽)。

    ★ 为什么必须查这个：英文串比中文长得多（一个汉字 ≈ 两个拉丁字母宽）。
      用 `bbox_inches="tight"` 保存的图会被动撑宽画布（字号相对变小），
      而**不带 tight 的图（fig7/fig8/figS1/figS3）会把超出的字直接裁掉** —— 静默丢字。
    """
    fig.canvas.draw()
    rend = fig.canvas.get_renderer()
    W = fig.canvas.get_width_height()[0]
    out = []
    for t in _collect(fig):
        s = t.get_text()
        if not s.strip():
            continue
        try:
            bb = t.get_window_extent(renderer=rend)
        except Exception:
            continue
        if bb.x0 < -tol_px or bb.x1 > W + tol_px:
            out.append((s, bb.x0, bb.x1, W))
    return out


def scan_overlap(fig, min_frac=None):
    """返回**互相压叠**的文字对 `(文本A, 文本B, 重叠率)`。

    ★ 为什么需要：溢出检查只管「有没有跑出画布」，管不了「画布内两块字叠在一起」。
      英文化后这类问题集中爆发（一个汉字 ≈ 2 个拉丁字母宽）：
        · fig5  轴外竖排「▲ negative control」压在 y 轴标签上；
        · figS1 面板 A/B 的长注解压到 y 轴标签；
        · fig7  x 轴三组「(epithelial fraction a%→b%)」彼此压叠；
        · figS3 图例盖住条上的分段百分比。
      重叠率按**较小那个框**的面积算（大字压小字要能报出来）。
      阈值可用环境变量 `FIG_OVERLAP_MIN` 临时调严（如 0.15 查「擦边」）。
    """
    if min_frac is None:
        min_frac = float(os.environ.get("FIG_OVERLAP_MIN", "0.32"))
    fig.canvas.draw()
    rend = fig.canvas.get_renderer()
    items = []
    for t in _collect(fig):
        s = t.get_text()
        if not s.strip():
            continue
        try:
            bb = t.get_window_extent(renderer=rend)
        except Exception:
            continue
        if bb.width <= 0 or bb.height <= 0:
            continue
        items.append((s, bb))
    out = []
    for i in range(len(items)):
        sa, ba = items[i]
        for j in range(i + 1, len(items)):
            sb, bb2 = items[j]
            x0, x1 = max(ba.x0, bb2.x0), min(ba.x1, bb2.x1)
            y0, y1 = max(ba.y0, bb2.y0), min(ba.y1, bb2.y1)
            if x1 <= x0 or y1 <= y0:
                continue
            inter = (x1 - x0) * (y1 - y0)
            small = min(ba.width * ba.height, bb2.width * bb2.height)
            if small > 0 and inter / small >= min_frac:
                out.append((sa, sb, inter / small))
    return out


def declutter(fig, artists, step=12, tries=10, min_frac=0.32):
    """把互相压叠的标注沿 y 方向错开，直到不再压叠或达到重试上限。

    返回实际移动次数。改的是 `Annotation` 的 `xytext`（offset points），
    **不动数据坐标** —— 标注的箭头/锚点仍指向原来的数据点。

    ★ 为什么需要：散点上的基因名是固定偏移量，中文名短所以看不出问题；
      换英文后 `FAP` 会压住 `COL11A1`（重叠率 81%）、`CLCA1` 压住 `FCGBP`（39%）。
    """
    moved = 0
    for _ in range(tries):
        fig.canvas.draw()
        rend = fig.canvas.get_renderer()
        boxes = []
        for t in artists:
            try:
                boxes.append((t, t.get_window_extent(renderer=rend)))
            except Exception:
                pass
        hit = None
        for i in range(len(boxes)):
            _, ba = boxes[i]
            for j in range(i + 1, len(boxes)):
                tb, bb2 = boxes[j]
                x0, x1 = max(ba.x0, bb2.x0), min(ba.x1, bb2.x1)
                y0, y1 = max(ba.y0, bb2.y0), min(ba.y1, bb2.y1)
                if x1 <= x0 or y1 <= y0:
                    continue
                inter = (x1 - x0) * (y1 - y0)
                small = min(ba.width * ba.height, bb2.width * bb2.height)
                if small > 0 and inter / small >= min_frac:
                    hit = tb
                    break
            if hit is not None:
                break
        if hit is None:
            break
        x, y = hit.get_position()
        hit.set_position((x, y + step))
        moved += 1
    return moved


def finish(fig, path, dpi=300, **kw):
    """保存前守卫 → `savefig` → 关闭画布。

    只要画布上还有汉字或全角标点，就打印明细并 `SystemExit(1)`，绝不静默出图。
    另外把横向溢出的文字与**互相压叠**的文字对打出来提醒
    （溢出在非 tight 保存下会被裁掉；压叠则是纯目视问题，脚本自己不会报错）。
    """
    left = scan(fig)
    if left:
        print("  !! %s 仍有 %d 处中文/全角文字，拒绝保存：" % (os.path.basename(path), len(left)))
        for s in left[:25]:
            print("       - %r" % s)
        raise SystemExit(1)
    over = scan_overflow(fig)
    if over:
        tight = "bbox_inches" in kw and kw["bbox_inches"] == "tight"
        print("  ~  %s 有 %d 处文字横向溢出画布（%s）："
              % (os.path.basename(path), len(over),
                 "画布会被撑宽" if tight else "★会被裁掉"))
        for s, x0, x1, W in over[:12]:
            print("       [%7.1f .. %7.1f] / %d  %r" % (x0, x1, W, s[:78]))
    ov = scan_overlap(fig)
    if ov:
        print("  ~  %s 有 %d 处文字互相压叠（目视可见，务必修）："
              % (os.path.basename(path), len(ov)))
        for sa, sb, f in sorted(ov, key=lambda x: -x[2])[:12]:
            print("       %3.0f%%  %-46r ×  %r" % (f * 100, sa[:44], sb[:44]))
    fig.savefig(path, dpi=dpi, **kw)
    print("  ok  %-52s %7.0f KB" % (os.path.basename(path), os.path.getsize(path) / 1024.0))
    plt.close(fig)


# ---------------------------------------------------------------------------
# 中文数字字体 → 一次都不需要了：全部英文后固定用 DejaVu Sans
# ---------------------------------------------------------------------------
def use_english_font():
    """统一英文字体。

    中文版脚本原先到处探测 Microsoft YaHei / SimHei；全部英文化之后不再需要，
    而且 YaHei 的 `·`(U+00B7) 与 `×` 字形偏怪，PMC/投稿系统常渲染不出来。
    DejaVu Sans 是 matplotlib 自带字体，`·` `×` `Δ` `−` `→` `∝` `≥` `★` `▲` `§`
    全部有字形（已用 FT2Font.get_charmap() 逐个核过）。
    """
    plt.rcParams["font.family"] = "DejaVu Sans"
    plt.rcParams["axes.unicode_minus"] = True


if __name__ == "__main__":
    print("PANELS_EN: %d 条" % len(PANELS_EN))
    for k, v in PANELS_EN.items():
        print("   %-24s -> %s" % (k, v))
    print("\nBAD 正则自测：")
    for s in ["abc", "Paneth-cell·4", "腺瘤", "（AF2）", "a｜b", "Δ = −0.400"]:
        print("   %-16r -> %s" % (s, "命中" if BAD.search(s) else "干净"))
