# 戊二酰化相关基因集（GRG）：两套小鼠肝 glutarylome（PRIDE）→ MGI 1:1 人同源 + 酶轴
import csv, os, collections
G = os.environ.get("GLUT_GENESET", "geneset")
MGI = os.environ.get("GLUT_MGI", "data/mgi/HOM_MouseHumanSequence.rpt")

def load(f):
    return {r["gene"] for r in csv.DictReader(open(f, encoding="utf-8"), delimiter="\t") if r["organism"] == "Mus musculus" and r["gene"]}
w = load(f"{G}/wagner2017_genes_q01.tsv"); s = load(f"{G}/schmiesing2018_genes_q01.tsv")
noorg = [r for r in csv.DictReader(open(f"{G}/schmiesing2018_genes_q01.tsv", encoding="utf-8"), delimiter="\t") if r["organism"] != "Mus musculus"]

# MGI 同源：按 DB Class Key 分组，只取 1 小鼠 : 1 人 的组
grp = collections.defaultdict(lambda: {"10090": [], "9606": []})
for r in csv.DictReader(open(MGI, encoding="utf-8"), delimiter="\t"):
    if r["NCBI Taxon ID"] in ("10090", "9606"): grp[r["DB Class Key"]][r["NCBI Taxon ID"]].append(r["Symbol"])
m2h = {}; multi = set()
for k, v in grp.items():
    for mg in v["10090"]:
        if len(v["10090"]) == 1 and len(v["9606"]) == 1: m2h[mg] = v["9606"][0]
        elif v["9606"]: multi.add(mg)

rows = []
for mg in sorted(w | s):
    src = ";".join(x for x, ss in (("Wagner2017", w), ("Schmiesing2018", s)) if mg in ss)
    hg = m2h.get(mg, "")
    rows.append(dict(mouse=mg, human=hg, source=src, both=int(mg in w and mg in s),
                     ortholog="1:1" if hg else ("multi" if mg in multi else "none")))
# 酶轴：Reactome R-HSA-71064 赖氨酸分解（12，与 ThemeScreen 同一份）+ 戊二酰化书写/擦除
lys = [l.strip() for l in open(os.environ.get("GLUT_REACTOME", "geneset/reactome_R-HSA-71064.txt"), encoding="utf-8") if l.strip()]
enz = {"SIRT5": "去戊二酰化酶（PMID 24703693）", "KAT2A": "组蛋白戊二酰转移酶（PMID 31542297，Bao 2019 Mol Cell）", "GCDH": "戊二酰辅酶A脱氢酶"}
tmp = f"{G}/GRG_substrate_mouse2human.tsv.tmp"
with open(tmp, "w", encoding="utf-8", newline="") as f:
    wr = csv.DictWriter(f, fieldnames=["mouse", "human", "source", "both", "ortholog"], delimiter="\t"); wr.writeheader(); wr.writerows(rows)
os.replace(tmp, tmp[:-4])
sub_h = sorted({r["human"] for r in rows if r["human"]})
core_h = sorted({r["human"] for r in rows if r["human"] and r["both"]})
grg = sorted(set(sub_h) | set(lys) | set(enz))
for name, lst in (("GRG_full", grg), ("GRG_substrate", sub_h), ("GRG_substrate_both", core_h), ("GRG_enzyme_axis", sorted(set(lys) | set(enz)))):
    open(f"{G}/{name}.txt.tmp", "w", encoding="utf-8").write("\n".join(lst) + "\n"); os.replace(f"{G}/{name}.txt.tmp", f"{G}/{name}.txt")
print("mouse genes: Wagner", len(w), "| Schmiesing", len(s), "| union", len(w | s), "| both", len(w & s))
print("Schmiesing non-mouse rows dropped:", [(r["gene"], r["organism"]) for r in noorg])
print("ortholog:", collections.Counter(r["ortholog"] for r in rows))
print("human substrate", len(sub_h), "| both-datasets", len(core_h), "| enzyme axis", len(set(lys) | set(enz)), "| GRG_full", len(grg))
print("both:", core_h)
