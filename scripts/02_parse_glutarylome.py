# 解析 PRIDE pepXML（Proteome Discoverer 导出），提取含戊二酰化赖氨酸的 rank-1 PSM
# 戊二酰化 K 残基质量 = 128.09496 (K) + 114.03169 (glutaryl, C5H6O3) = 242.12665
# 可信鉴定：Percolator q-Value <= 0.01（PSM 级 1% FDR，蛋白组学常规阈值）
# 用法：python -I 02_parse_glutarylome.py <pep.xml> <out_prefix> <dataset_label> [alt_K_mass]
# alt_K_mass：PD 导出把戊二酰化 K 的 mod_aminoacid_mass 写错时用（PXD005895 写成 129.095）。
#   接受条件：该 K 位点质量≈alt 且 母离子质量反推该位点增量 = 114.0317±0.01（逐 PSM 校验，不过即弃）
import sys, re, csv, os
import xml.etree.ElementTree as ET
src, out, label = sys.argv[1], sys.argv[2], sys.argv[3]
ALT = float(sys.argv[4]) if len(sys.argv) > 4 else None
R = dict(G=57.02146,A=71.03711,S=87.03203,P=97.05276,V=99.06841,T=101.04768,C=103.00919,L=113.08406,I=113.08406,N=114.04293,D=115.02694,Q=128.05858,K=128.09496,E=129.04259,M=131.04049,H=137.05891,F=147.06841,R=156.10111,Y=163.06333,W=186.07931)
n_alt_rej = 0
GLU_K = 242.12665; TOL = 0.01; QMAX = 0.01
NS = "{http://regis-web.systemsbiology.net/pepXML}"
psm = []; kmass = {}
n_hit = 0
for ev, el in ET.iterparse(src, events=("end",)):
    if el.tag != NS + "search_hit":
        if el.tag == NS + "spectrum_query": el.clear()
        continue
    if el.get("hit_rank") != "1": continue
    n_hit += 1
    pep = el.get("peptide")
    prots = [(el.get("protein"), el.get("protein_descr") or "")]
    for ap in el.findall(NS + "alternative_protein"):
        prots.append((ap.get("protein"), ap.get("protein_descr") or ""))
    q = None
    for s in el.findall(NS + "search_score"):
        if s.get("name") in ("Percolator q-Value", "q-Value", "Percolator q-value"):
            q = float(s.get("value"))
    sites = []; alt_sites = []; mods = {}
    for mi in el.findall(NS + "modification_info"):
        for m in mi.findall(NS + "mod_aminoacid_mass"):
            mods[int(m.get("position"))] = float(m.get("mass"))
    for pos, mass in mods.items():
        aa = pep[pos - 1] if 0 < pos <= len(pep) else "?"
        if aa == "K": kmass[round(mass, 3)] = kmass.get(round(mass, 3), 0) + 1
        if aa == "K" and abs(mass - GLU_K) < TOL: sites.append(pos)
        if aa == "K" and ALT is not None and abs(mass - ALT) < TOL: alt_sites.append(pos)
    if alt_sites:
        # 母离子校验：其余 K 视为 TMT6（229.16293），N 端 TMT6，C 固定 +57.02146，M 氧化 +15.99491
        calc = float(el.get("calc_neutral_pep_mass"))
        base = sum(R.get(a, 0) for a in pep) + 18.01056 + 229.16293
        nC = pep.count("C"); nM = sum(1 for p, ms in mods.items() if pep[p-1] == "M" and abs(ms - 147.0354) < 0.01)
        nT = pep.count("K") - len(alt_sites)
        per = (calc - base - 57.02146 * nC - 15.99491 * nM - 229.16293 * nT) / len(alt_sites)
        if abs(per - 114.03169) < TOL: sites += alt_sites
        else: n_alt_rej += 1
    if sites:
        for acc, d in prots:
            gn = re.search(r"GN=(\S+)", d); os_ = re.search(r"OS=(.+?) (?:GN|PE)=", d)
            psm.append(dict(dataset=label, protein=acc, gene=gn.group(1) if gn else "", organism=os_.group(1) if os_ else "",
                            peptide=pep, glu_pos=";".join(map(str, sites)), q=q))
tmp = out + "_psm.tsv.tmp"
with open(tmp, "w", encoding="utf-8", newline="") as f:
    w = csv.DictWriter(f, fieldnames=["dataset", "protein", "gene", "organism", "peptide", "glu_pos", "q"], delimiter="\t")
    w.writeheader(); w.writerows(psm)
os.replace(tmp, out + "_psm.tsv")
conf = [r for r in psm if r["q"] is not None and r["q"] <= QMAX]
genes = {}
for r in conf:
    g = genes.setdefault((r["gene"], r["organism"]), dict(n_psm=0, peps=set()))
    g["n_psm"] += 1; g["peps"].add(r["peptide"])
tmp = out + "_genes_q01.tsv.tmp"
with open(tmp, "w", encoding="utf-8", newline="") as f:
    f.write("dataset\tgene\torganism\tn_psm\tn_unique_peptides\n")
    for (g, o), v in sorted(genes.items()):
        f.write(f"{label}\t{g}\t{o}\t{v['n_psm']}\t{len(v['peps'])}\n")
os.replace(tmp, out + "_genes_q01.tsv")
top = sorted(kmass.items(), key=lambda x: -x[1])[:8]
print(label, "rank1_hits", n_hit, "| glutaryl-K PSM (any q)", len(psm), "| q<=0.01", len(conf),
      "| genes q<=0.01", len(genes), "| q missing", sum(r["q"] is None for r in psm), "| alt-mass PSM rejected by precursor check", n_alt_rej)
print(" K modified masses (top):", top)
print(" organisms:", sorted({o for (_, o) in genes}))
