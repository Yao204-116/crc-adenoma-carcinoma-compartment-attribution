# -*- coding: utf-8 -*-
"""
build_suppl_tables.py
归档 Supplementary Tables S1-S4：
  S1  313-gene per-gene attribution table
  S2  Curated panel gene membership
  S3  13 panels x 3 cohorts x 2 calibers pattern matrix (+ sensitivity sub-panel)
  S4  GSE161277 single-cell: sample QC table + T=PxU decomposition table
输出：suppl/*.csv + suppl/Supplementary_Tables_S1-S4.xlsx
数字全部来自落盘 JSON/CSV，不重新计算。
"""
import csv, json, os

BASE = os.environ.get("AACA_BASE") or os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(BASE, 'suppl')
os.makedirs(OUT, exist_ok=True)

# 上游 R 产物（WGCNA turquoise 模块 / 富集交集枢纽）。仓库副本放在 upstream_r/ 下；
# 该目录不存在时（本机开发环境）回退到原始位置 D:/2。
_UPSTREAM = os.environ.get("AACA_UPSTREAM_R") or os.path.join(BASE, 'upstream_r')
UP = _UPSTREAM if os.path.isdir(_UPSTREAM) else r"D:/2"

def R(p):
    with open(os.path.join(BASE, p), encoding='utf-8') as f:
        return json.load(f)

def Rd(p):
    with open(os.path.join(BASE, p), encoding='utf-8-sig', newline='') as f:
        return list(csv.reader(f))

ATTR_EN = {
    '上皮内在（LCM 保留）': 'Epithelial-intrinsic (preserved in LCM)',
    '组成来源（LCM 反向）': 'Composition-derived (reversed in LCM)',
    '全组织不一致': 'Whole-tissue inconsistent (unclassified)',
}
TRANS = {('normal', 'adenoma'): 'normal -> adenoma',
         ('adenoma', 'carcinoma'): 'adenoma -> carcinoma'}
SIG = {'paneth': 'Paneth', 'diff_abs': 'Differentiation/absorptive',
       'goblet': 'Goblet', 'stem': 'Stem/progenitor'}

# ---------------------------------------------------------------- S1
rows = Rd('suppl_table1_313genes.csv')
hdr = rows[0]
S1 = [['gene', 'tier', 'mean_logFC', 'q_stouffer',
       'd_GSE4183_whole', 'd_GSE117606_whole', 'd_GSE15960_LCM',
       'n_whole_tissue_cohorts', 'n_agree', 'n_LCM', 'n_agree_LCM',
       'attribution', 'attribution_full']]
for r in rows[1:]:
    d = dict(zip(hdr, r))
    S1.append([d['gene'], d['tier_change'], d['mean_logFC'], d['q_stouffer'],
               d['d_GSE4183'], d['d_GSE117606'], d['d_GSE15960_LCM'],
               d['n_tissue'], d['tissue_agree'], d['n_lcm'], d['lcm_agree'],
               ATTR_EN.get(d['attribution'], d['attribution']).split(' (')[0],
               ATTR_EN.get(d['attribution'], d['attribution'])])
def w_csv(name, rows):
    with open(os.path.join(OUT, name), 'w', encoding='utf-8-sig', newline='') as f:
        csv.writer(f).writerows(rows)

w_csv('S1_313gene_attribution.csv', S1)

# ---------------------------------------------------------------- S2
pd_ = R('panels_def.json')
tier = {d['gene']: d['tier_change'] for d in
        (dict(zip(rows[0], r)) for r in rows[1:])}
S2 = [['panel', 'n_genes', 'order_in_panel', 'gene', 'component_of_313_signature']]
for panel, genes in pd_.items():
    for i, g in enumerate(genes, 1):
        S2.append([panel, len(genes), i, g, tier.get(g, 'not in 313')])
# 两个由 313 签名切出的分量面板
up = [d['gene'] for d in (dict(zip(rows[0], r)) for r in rows[1:]) if d['tier_change'] == 'UP']
dn = [d['gene'] for d in (dict(zip(rows[0], r)) for r in rows[1:]) if d['tier_change'] == 'DOWN']
for panel, genes in [('Up-regulated component', up), ('Down-regulated component', dn)]:
    for i, g in enumerate(genes, 1):
        S2.append([panel, len(genes), i, g, panel.split('-')[0].lower()])

# 用户自有基因集（WGCNA turquoise 模块 + PPI Hub）——R2 ¶3/¶5 引用，必须可查
def genes_from(path, col, filt=None):
    out = []
    if not os.path.exists(path):
        return out
    with open(path, encoding='utf-8-sig', newline='') as f:
        for r in csv.DictReader(f):
            if filt and not filt(r):
                continue
            v = (r.get(col) or '').strip().upper()
            if v:
                out.append(v)
    return out

TURQ = genes_from(os.path.join(UP, '02.WGCNA', 'WGCNA_top_module_genes.csv'), 'Gene',
                  lambda r: (r.get('Module') or '').strip() == 'turquoise'
                            and (r.get('Comparison') or '').strip() == 'Tumor_vs_Polyp')
HUB = genes_from(os.path.join(UP, '03.enrich', 'Hub_Tumor_vs_Polyp.csv'), 'Gene')
UPS, DNS = set(up), set(dn)
for panel, genes in [('WGCNA-turquoise (Tumor vs Polyp)', TURQ),
                     ('turquoise \u2229 up-component', [g for g in TURQ if g in UPS]),
                     ('turquoise \u2229 down-component', [g for g in TURQ if g in DNS]),
                     ('PPI Hub genes (Tumor vs Polyp)', HUB)]:
    for i, g in enumerate(genes, 1):
        S2.append([panel, len(genes), i, g, tier.get(g, 'not in 313')])
print('S2 扩展集：turquoise=%d, turq∩up=%d, turq∩dn=%d, Hub=%d' % (
    len(TURQ), len(UPS & set(TURQ)), len(DNS & set(TURQ)), len(HUB)))
w_csv('S2_panel_gene_membership.csv', S2)

# ---------------------------------------------------------------- S3
r3 = R('r3_facts.json')
pt = r3['panel_table']
COH = ['GSE4183', 'GSE117606', 'GSE15960']
NG = {'上调分量(210)': 210, '下调分量(103)': 103}
LAB = {'P': 'P (adenoma peak)', 'D': 'D (monotone decrease)', 'A': 'A (already reduced)'}
S3 = [['panel', 'n_genes_from_S2',
       'GSE4183_raw', 'GSE4183_drift-corrected',
       'GSE117606_raw', 'GSE117606_drift-corrected',
       'GSE15960_LCM_raw', 'GSE15960_LCM_drift-corrected',
       'pattern_raw', 'pattern_corrected',
       'robust_adenoma_peak_3of3_both_calibers', 'failure_mode']]
for panel, v in pt.items():
    rp, ap = v['raw_patterns'], v['adj_patterns']
    robust = v.get('robust_P_3of3', False)
    if robust:
        fm = ''
    else:
        raw_p = [c for c in COH if rp[c] == 'P']
        adj_p = [c for c in COH if ap[c] == 'P']
        if raw_p and not adj_p:
            fm = 'adenoma peak lost after drift correction'
        elif len(raw_p) < 3:
            fm = 'peak in at most one cohort'
        else:
            fm = 'inconsistent between cohorts / calibers'
    S3.append([panel, len(pd_.get(panel, [])) if panel in pd_ else NG.get(panel, ''),
               rp['GSE4183'], ap['GSE4183'],
               rp['GSE117606'], ap['GSE117606'],
               rp['GSE15960'], ap['GSE15960'],
               v['raw_str'], v['adj_str'],
               'YES' if robust else 'no', fm])
w_csv('S3_panel_pattern_matrix.csv', S3)

sens = R('r3_sens_paneth.json')
S3b = [['sub-panel', 'cohort', 'caliber', 'normal', 'adenoma', 'carcinoma',
        'pattern', 'ordering_normal<ca<ad_holds']]
def walk_sens(o, prefix=''):
    for k, v in o.items():
        if isinstance(v, dict):
            walk_sens(v, prefix + k + ' / ')
        else:
            pass
# sens 结构按子面板 -> 队列 -> {raw,adj} 展开
for sub, cohorts in sens.items():
    if not isinstance(cohorts, dict):
        continue
    for coh, cal in cohorts.items():
        if not isinstance(cal, dict) or 'raw' not in cal:
            continue
        for cal_name in ['raw', 'adj']:
            m = cal.get(cal_name, {})
            if not isinstance(m, dict):
                continue
            no, ad, ca = m.get('normal'), m.get('adenoma'), m.get('carcinoma')
            S3b.append([sub, coh, cal_name,
                        round(no, 3) if no is not None else '',
                        round(ad, 3) if ad is not None else '',
                        round(ca, 3) if ca is not None else '',
                        m.get('pattern', ''),
                        'yes' if (no is not None and ad is not None and ca is not None
                                  and no < ca and ca < ad) else 'no'])
w_csv('S3b_paneth_sensitivity.csv', S3b)

# ---------------------------------------------------------------- S4
scan = R('sc161277_scan.json')
pairs = R('sc161277_pairs.json')
ps = {x['gsm']: x for x in pairs['defa_persample']}

S4a = [['GSM', 'patient', 'stage', 'n_cells_raw', 'n_cells_after_QC',
        'QC_retention_%', 'median_UMI_after_QC', 'n_epithelial',
        'epithelial_%_of_QC', 'included_in_proportion_decomposition', 'note']]
for s in scan:
    gsm = s['gsm']
    n_epi = ps.get(gsm, {}).get('n_epi', '')
    frac = ps.get(gsm, {}).get('epi_frac', s.get('pct_epi_of_qc'))
    pct = round(100 * frac, 1) if isinstance(frac, float) else ''
    inc = 'yes' if (isinstance(pct, float) and pct >= 5 and s['n_cells_qc'] > 500) else 'no (low epithelial yield)'
    note = ''
    if s['stage'] == 'blood':
        inc, note = 'no (negative control)', 'zero epithelial cells under the EPCAM/KRT threshold'
    elif gsm in ('GSM4904242', 'GSM4904243'):
        note = 'second adenoma of Patient 3; pooled by cell-count weighting'
    S4a.append([gsm, s['patient'], s['stage'], s['n_cells_raw'], s['n_cells_qc'],
                round(100 * s['n_cells_qc'] / s['n_cells_raw'], 1),
                s['median_umi'], n_epi, pct, inc, note])
w_csv('S4a_sc161277_sample_QC.csv', S4a)

S4b = [['patient', 'transition', 'signature',
        'T_start', 'T_end', 'U_start', 'U_end', 'P_start', 'P_end',
        'dlogT', 'dlogP', 'dlogU',
        'residual_dlogT-(dlogP+dlogU)',
        'share_composition_%', 'share_state_%', 'decidable', 'state_contribution_dominant']]
for x in pairs['decomp']:
    dT, dP, dU = x['dlogT'], x['dlogP'], x['dlogU']
    denom = dP + dU
    share_p = round(100 * dP / denom, 1) if abs(denom) > 1e-9 else ''
    share_u = round(100 * dU / denom, 1) if abs(denom) > 1e-9 else ''
    dom = ''
    if isinstance(share_u, float) and isinstance(share_p, float):
        dom = 'state (U)' if abs(share_u) >= abs(share_p) else 'composition (P)'
    dec = 'yes' if abs(dT) >= 0.30 else 'no (|dlogT| < 0.3: shares not interpretable)'
    S4b.append([x['patient'], TRANS.get((x['a'], x['b']), f"{x['a']} -> {x['b']}"),
                SIG.get(x['sig'], x['sig']),
                round(x['T_a'], 4), round(x['T_b'], 4),
                round(x['U_a'], 4), round(x['U_b'], 4),
                round(x['P_a'], 4), round(x['P_b'], 4),
                round(dT, 4), round(dP, 4), round(dU, 4),
                round(dT - denom, 4),
                share_p, share_u,
                dec, dom])
w_csv('S4b_sc161277_decomposition.csv', S4b)

# ---------------------------------------------------------------- S5
# 四队列逐队列性能表（R1 ¶4 引用；口径 = GSE117606 as-annotated 69 例）
fc = R('manuscript_facts.json')['four_cohort_auc']
UPK = '上调分量(间质/髓系激活)'
DNK = '下调分量(上皮功能关闭)'
NOTE = {
    'GSE117606': 'adenoma n = 69 as annotated (includes 10 sessile serrated adenomas); '
                 'excluding SSA (n = 59) leaves whole-signature and up-component AUCs within 0.01 '
                 '(0.921 and 0.906) — all downstream analyses use the 59-adenoma definition',
    'GSE128435': 'excluded from all downstream analyses: gene-level concordance 51% (chance) '
                 'points to a technical rather than a biological discrepancy',
    'GSE4183': '15 IBD samples present in the series are not used in any contrast reported here',
}
S5 = [['cohort', 'platform', 'design', 'n_normal', 'n_adenoma', 'n_carcinoma',
       'auc_signature', 'auc_up_component', 'auc_down_component',
       'gene_direction_concordance', 'n_genes_tested',
       'auc_adenoma_vs_normal', 'note']]
for coh, v in fc.items():
    c = v['counts']
    S5.append([coh, v['gpl'],
               '%s normal / %s adenoma / %s carcinoma' % (c.get('normal', ''), c.get('adenoma', ''), c.get('carcinoma', '')),
               c.get('normal', ''), c.get('adenoma', ''), c.get('carcinoma', ''),
               round(v['auc'], 3),
               round(v['component_auc'][UPK], 3), round(v['component_auc'][DNK], 3),
               round(v['gene_dir_concordance'], 3), v['n_tested'],
               round(v['adenoma_vs_normal_auc']['signature'], 3),
               NOTE.get(coh, '')])
w_csv('S5_four_cohort_performance.csv', S5)

# ---------------------------------------------------------------- XLSX
from openpyxl import Workbook
from openpyxl.styles import Font, PatternFill, Alignment, Border, Side
from openpyxl.utils import get_column_letter

wb = Workbook()
wb.remove(wb.active)
HEAD = PatternFill('solid', fgColor='1E293B')
HF = Font(color='FFFFFF', bold=True, size=10)
BODY = Font(size=10)
THIN = Border(*[Side(style='thin', color='D8DEE9')] * 4)

def add_sheet(title, rows, widths=None, freezes='A2'):
    ws = wb.create_sheet(title[:31])
    for r in rows:
        ws.append(r)
    for j, c in enumerate(ws[1], 1):
        c.fill, c.font = HEAD, HF
        c.alignment = Alignment(horizontal='center', vertical='center', wrap_text=True)
    for row in ws.iter_rows(min_row=2):
        for c in row:
            c.font = BODY
            c.border = THIN
    if widths:
        for j, w in enumerate(widths, 1):
            ws.column_dimensions[get_column_letter(j)].width = w
    ws.freeze_panes = freezes
    ws.row_dimensions[1].height = 30
    return ws

add_sheet('S1_313gene_attribution', S1,
          [11, 8, 11, 11, 12, 13, 13, 11, 9, 8, 10, 26, 42])
add_sheet('S2_panel_membership', S2, [26, 9, 14, 12, 22])
add_sheet('S3_panel_pattern_matrix', S3,
          [24, 8, 13, 18, 14, 19, 15, 20, 11, 12, 14, 40])
add_sheet('S3b_paneth_sensitivity', S3b, [22, 12, 9, 10, 10, 11, 10, 16])
add_sheet('S4a_sample_QC', S4a, [13, 11, 13, 11, 12, 11, 12, 11, 12, 20, 40])
add_sheet('S4b_decomposition', S4b,
          [11, 20, 22, 9, 9, 9, 9, 9, 9, 9, 9, 9, 16, 14, 12, 24, 20])
add_sheet('S5_four_cohort_performance', S5, [12, 10, 34, 9, 10, 11, 13, 15, 16, 15, 12, 15, 70])
xlsx = os.path.join(OUT, 'Supplementary_Tables_S1-S5.xlsx')
wb.save(xlsx)

print('CSV rows: S1=%d S2=%d S3=%d S3b=%d S4a=%d S4b=%d S5=%d' %
      (len(S1), len(S2), len(S3), len(S3b), len(S4a), len(S4b), len(S5)))
print('XLSX ->', xlsx, os.path.getsize(xlsx), 'bytes')
# 自检
print('S4a QC 合计 raw=%d qc=%d epi=%d' % (
    sum(r[3] for r in S4a[1:]), sum(r[4] for r in S4a[1:]), sum(r[7] for r in S4a[1:])))
print('S3 robust 面板:', [r[0] for r in S3[1:] if r[10] == 'YES'])
