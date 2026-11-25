"""Fixed, independently sourced program scores and IR-Ctrl descriptive comparisons.

Scores are mean log1p(CP10k) of fixed genes minus mean log1p(CP10k) of
expression-bin-matched control genes. They are not measurements of function.
No animal-level p-values are calculated from cells.
"""
from pathlib import Path
import json,hashlib
import numpy as np
import pandas as pd
from scipy import sparse
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

ROOT=Path(__file__).resolve().parents[2]
OUT=Path(__file__).resolve().parent
T=OUT/'tables';F=OUT/'figures';R=OUT/'resources'
genes,cells=json.loads((ROOT/'reanalysis/Result1/tables/seurat_RNA_counts_dimnames.json').read_text())
m=pd.read_csv(ROOT/'reanalysis/Result1/tables/metadata_harmonized.csv',index_col=0).loc[cells]
cts=['HSC','HSC-Cycle','MPP','MkP','PreE','NMP','NMP-Cycle','NP','MP','DC',
     'BP-1','BP-2','BP-3','EP','CLP','ILCP','ILC1','ILC2','T Cell','B Cell']
groups=['Ctrl','IR','DMSO','DMSO-IR']

def gmt(path):
 return {p[0]:list(dict.fromkeys(p[2:])) for line in path.read_text().splitlines() if (p:=line.split('\t'))}
hall=gmt(R/'mh.all.v2025.1.Mm.symbols.gmt');go=gmt(R/'m5.go.bp.v2025.1.Mm.symbols.gmt')
# Freeze these names before examining IR-Ctrl score effects. The full collection is retained.
chosen={
 'p53-associated expression':hall['HALLMARK_P53_PATHWAY'],
 'Apoptosis-associated expression':hall['HALLMARK_APOPTOSIS'],
 'DNA repair':hall['HALLMARK_DNA_REPAIR'],
 'Double-strand break repair':go['GOBP_DOUBLE_STRAND_BREAK_REPAIR'],
 'Ribosome biogenesis':go['GOBP_RIBOSOME_BIOGENESIS'],
 'Cytoplasmic translation':go['GOBP_CYTOPLASMIC_TRANSLATION'],
 'rRNA processing':go['GOBP_RRNA_PROCESSING'],
 'RNA splicing':go['GOBP_RNA_SPLICING'],
 'E2F targets':hall['HALLMARK_E2F_TARGETS'],
 'G2M checkpoint':hall['HALLMARK_G2M_CHECKPOINT'],
 'MYC targets V1':hall['HALLMARK_MYC_TARGETS_V1'],
 'HSC differentiation-associated expression':go['GOBP_HEMATOPOIETIC_STEM_CELL_DIFFERENTIATION'],
}
(R/'fixed_program_gene_lists.json').write_text(json.dumps(chosen,indent=2))
print('Fixed 12 program sets',flush=True)

x=sparse.load_npz(ROOT/'reanalysis/Result1/tables/seurat_RNA_data.npz').tocsc()
assert x.shape==(len(genes),len(cells))
mean=np.asarray(x.mean(axis=1)).ravel()
positive=np.flatnonzero(mean>0)
bins=np.full(len(genes),-1,int)
order=positive[np.argsort(mean[positive],kind='stable')]
bins[order]=np.minimum(np.arange(len(order))*24//len(order),23)
all_targets={g for gs in chosen.values() for g in gs}
gene_to_index={g:i for i,g in enumerate(genes)}
scores={};definitions={};coverage=[]
for k,(name,gs) in enumerate(chosen.items()):
 target=np.array([gene_to_index[g] for g in gs if g in gene_to_index and mean[gene_to_index[g]]>0])
 excluded=set(target)
 # Controls exclude this entire target set, with deterministic per-gene matching.
 rng=np.random.default_rng(61019+k)
 controls=[]
 for i in target:
  pool=np.array([j for j in np.flatnonzero(bins==bins[i]) if j not in excluded])
  controls.extend(rng.choice(pool,size=min(50,len(pool)),replace=False).tolist())
 control=np.unique(controls)
 definitions[name]={'target_genes':[genes[i] for i in target],
                    'control_genes':[genes[i] for i in control], 'seed':61019+k}
 coverage.append([name,len(gs),len(target),len(control)])
 scores[name]=np.asarray(x[target].mean(axis=0)-x[control].mean(axis=0)).ravel()
(R/'score_definitions_and_matched_controls.json').write_text(json.dumps(definitions,indent=2))
pd.DataFrame(coverage,columns=['program','source_genes','measured_expressed_genes','control_genes']).to_csv(T/'program_coverage.csv',index=False)
s=pd.DataFrame(scores,index=cells);s.index.name='cell'
s.to_csv(T/'program_scores.csv.gz')
ms=m[['orig.ident','celltype']].join(s)
med=ms.groupby(['celltype','orig.ident'],observed=True)[list(chosen)].median()
avg=ms.groupby(['celltype','orig.ident'],observed=True)[list(chosen)].mean()
med.to_csv(T/'program_medians_four_conditions.csv');avg.to_csv(T/'program_means_four_conditions.csv')
delta=(med.xs('IR',level=1)-med.xs('Ctrl',level=1)).reindex(cts)
delta.to_csv(T/'program_IR_minus_Ctrl_median.csv')
n=pd.crosstab(m.celltype,m['orig.ident']).reindex(index=cts,columns=groups)
valid=[ct for ct in cts if n.loc[ct,['IR','Ctrl']].min()>=20]
pd.DataFrame({'increases':(delta.loc[valid]>0).sum(),'decreases':(delta.loc[valid]<0).sum(),
 'celltypes_with_at_least_20_cells_per_condition':len(valid)}).to_csv(T/'program_direction_summary.csv')

# Direct marker changes, independent of cell-level significance cutoffs.
av=pd.read_csv(ROOT/'reanalysis/Result1/tables/celltype_group_mean_normalized_expression.csv.gz',index_col=0)
records=[]
for gene in ['Ccng1','Cdkn1a','Phlda3','Trp53inp1','Bbc3','Bax','Fos','Junb','Klf6','Gata2']:
 for ct in cts:
  a=av.loc[gene,f'Ctrl|{ct}'];b=av.loc[gene,f'IR|{ct}']
  records.append([gene,ct,a,b,np.log2((b+1)/(a+1)),n.loc[ct,'Ctrl'],n.loc[ct,'IR']])
mark=pd.DataFrame(records,columns=['gene','celltype','Ctrl_mean_CP10k','IR_mean_CP10k','log2_ratio_plus_1','Ctrl_n','IR_n'])
mark.to_csv(T/'key_gene_expression_changes.csv',index=False)

# Cell-cycle fractions are inferred states; they do not prove cell-cycle arrest.
ph=pd.crosstab([m.celltype,m['orig.ident']],m.Phase,normalize='index')
ph.to_csv(T/'cellcycle_fractions_by_type_and_group.csv')
pd.crosstab(m['orig.ident'],m.Phase,normalize='index').reindex(groups).to_csv(T/'cellcycle_fractions_all_cells.csv')

# Original enrichment is audited, not treated as independently replicated testing.
for label,fn in [('up','2.Conserved_DEG.GO.csv'),('down','2.Conserved_down.GO.csv')]:
 e=pd.read_csv(ROOT/'result/19_HSPC_specific_DEG'/fn)
 e[e.ONTOLOGY=='BP'].head(20).to_csv(T/f'original_enrichment_{label}_top20.csv',index=False)

plt.rcParams.update({'font.family':'DejaVu Sans','font.size':9,'svg.fonttype':'none','savefig.dpi':200})
def save(fig,name):
 for ext in ['png','svg']:fig.savefig(F/f'{name}.{ext}',bbox_inches='tight')
 plt.close(fig)

fig,ax=plt.subplots(figsize=(12.3,6.5),layout='constrained')
v=delta.T.to_numpy();lim=max(np.quantile(abs(v),.95),.05)
im=ax.imshow(v,cmap='RdBu_r',vmin=-lim,vmax=lim,aspect='auto')
ax.set(xticks=range(20),xticklabels=[c+('*' if c not in valid else '') for c in cts],yticks=range(12),yticklabels=list(chosen),title='IR minus Ctrl in independently defined expression programs')
ax.tick_params(axis='x',rotation=60)
fig.colorbar(im,ax=ax,shrink=.7,label='Difference in median matched-control score')
fig.text(.01,-.025,'* Fewer than 20 cells in at least one condition. Colors describe score shifts, not functional activity or animal-level significance.',fontsize=8)
save(fig,'program_response_heatmap')

fig,axs=plt.subplots(1,2,figsize=(11,5),layout='constrained')
for ax,sets in zip(axs,[['p53-associated expression','Apoptosis-associated expression','Double-strand break repair'],['Ribosome biogenesis','Cytoplasmic translation','RNA splicing']]):
 for name in sets:
  ax.plot(range(20),delta[name],marker='o',markersize=3,label=name)
 ax.axhline(0,color='grey',lw=.8);ax.set(xticks=range(20),xticklabels=cts,ylabel='IR minus Ctrl score')
 ax.tick_params(axis='x',rotation=90);ax.legend(fontsize=8,frameon=False)
save(fig,'program_response_profiles')

sub=mark[mark.gene.isin(['Ccng1','Phlda3','Trp53inp1','Fos','Junb','Klf6'])].pivot(index='gene',columns='celltype',values='log2_ratio_plus_1').reindex(columns=cts)
fig,ax=plt.subplots(figsize=(11,3.5),layout='constrained');im=ax.imshow(sub,aspect='auto',cmap='RdBu_r',vmin=-2,vmax=2)
ax.set(xticks=range(20),xticklabels=cts,yticks=range(len(sub)),yticklabels=sub.index,title='Representative RNA expression changes')
ax.tick_params(axis='x',rotation=60);fig.colorbar(im,ax=ax,label='log2[(IR mean CP10k + 1)/(Ctrl mean CP10k + 1)]',shrink=.8)
save(fig,'key_gene_change_heatmap')

print(pd.read_csv(T/'program_direction_summary.csv').to_string(index=False),flush=True)
print('Marker directions in valid populations',mark[mark.celltype.isin(valid)].groupby('gene').log2_ratio_plus_1.agg(lambda z:int((z>0).sum())).to_dict(),flush=True)
