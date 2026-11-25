"""Fixed-depth UMI checks and associations across cell types (not across duplicated AUC cells)."""
from pathlib import Path
import json
import numpy as np
import pandas as pd
from scipy import sparse,stats
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

ROOT=Path(__file__).resolve().parents[2];OUT=Path(__file__).resolve().parent
T=OUT/'tables';F=OUT/'figures';R=OUT/'resources'
genes,cells=json.loads((ROOT/'reanalysis/Result1/tables/seurat_RNA_counts_dimnames.json').read_text())
m=pd.read_csv(ROOT/'reanalysis/Result1/tables/metadata_harmonized.csv',index_col=0).loc[cells]
counts=sparse.load_npz(ROOT/'reanalysis/Result1/tables/seurat_RNA_counts.npz').tocsc()
assert counts.shape==(len(genes),len(cells))
depth=1000
assert np.asarray(counts.sum(axis=0)).min()>=depth
print('Downsampling all 23,006 cells to exactly 1,000 UMIs, retaining every cell',flush=True)
rng=np.random.default_rng(712019)
data=[];indices=[];indptr=[0];rich=[]
for j in range(counts.shape[1]):
 lo,hi=counts.indptr[j:j+2]
 v=rng.multivariate_hypergeometric(counts.data[lo:hi].astype(np.int64),depth,method='marginals')
 keep=v>0
 data.extend(v[keep].tolist());indices.extend(counts.indices[lo:hi][keep].tolist());indptr.append(len(data));rich.append(int(keep.sum()))
 if (j+1)%6000==0:print('Downsampled',j+1,flush=True)
down=sparse.csc_matrix((np.array(data,dtype=np.int32),np.array(indices,dtype=np.int32),np.array(indptr,dtype=np.int64)),shape=counts.shape)
assert (np.asarray(down.sum(axis=0)).ravel()==depth).all()
sparse.save_npz(T/'counts_1000_UMI.npz',down)
m['nFeature_1000UMI']=rich
m[['orig.ident','celltype','Phase','nCount_RNA','nFeature_RNA','nFeature_1000UMI']].to_csv(T/'depth_matched_cell_metrics.csv')
defs=json.loads((R/'score_definitions_and_matched_controls.json').read_text())
gi={g:i for i,g in enumerate(genes)}
logdown=down.astype(float);logdown.data=np.log1p(logdown.data*10)
sc={}
for name,d in defs.items():
 target=[gi[g] for g in d['target_genes']];ctrl=[gi[g] for g in d['control_genes']]
 sc[name]=np.asarray(logdown[target].mean(axis=0)-logdown[ctrl].mean(axis=0)).ravel()
sd=pd.DataFrame(sc,index=cells)
sd.to_csv(T/'program_scores_1000UMI.csv.gz')
med=m[['orig.ident','celltype']].join(sd).groupby(['celltype','orig.ident'],observed=True)[list(defs)].median()
med.to_csv(T/'program_medians_four_conditions_1000UMI.csv')
delta=med.xs('IR',level=1)-med.xs('Ctrl',level=1)
delta.to_csv(T/'program_IR_minus_Ctrl_1000UMI.csv')
n=pd.crosstab(m.celltype,m['orig.ident']);valid=list(n.index[n[['Ctrl','IR']].min(axis=1)>=20])
pd.DataFrame({'increases':(delta.loc[valid]>0).sum(),'decreases':(delta.loc[valid]<0).sum()}).to_csv(T/'program_direction_summary_1000UMI.csv')

# Repeated exact-depth rarefaction of Ctrl for the baseline richness association.
ctypes=sorted(valid)
ctrl_indices=np.flatnonzero(m['orig.ident'].to_numpy()=='Ctrl')
rep_rows=[]
for rep in range(10):
 rng=np.random.default_rng(813019+rep)
 values=[]
 for j in ctrl_indices:
  lo,hi=counts.indptr[j:j+2]
  v=rng.multivariate_hypergeometric(counts.data[lo:hi].astype(np.int64),depth,method='marginals')
  values.append(np.count_nonzero(v))
 rich_by_type=pd.DataFrame({'celltype':m.iloc[ctrl_indices].celltype.to_numpy(),'richness':values}).groupby('celltype').richness.median()
 for ct in ctypes:rep_rows.append([rep,ct,float(rich_by_type.loc[ct])])
 print('Ctrl rarefaction replicate',rep+1,flush=True)
rep=pd.DataFrame(rep_rows,columns=['replicate','celltype','median_detected_genes_1000UMI'])
rep.to_csv(T/'Ctrl_richness_1000UMI_repeats.csv',index=False)

ctmap=m.groupby('integrated_snn_res.0.6').celltype.first()
a=pd.read_csv(T/'original_Augur_IR_Ctrl.csv');a['celltype']=a.cell_type.map(ctmap)
cur=m[m['orig.ident']=='Ctrl'].groupby('celltype').agg(median_UMI=('nCount_RNA','median'),median_genes=('nFeature_RNA','median'))
cur=cur.join(pd.crosstab(m[m['orig.ident']=='Ctrl'].celltype,m[m['orig.ident']=='Ctrl'].Phase,normalize='index').add_prefix('fraction_'))
cur=cur.join(rep.groupby('celltype').median_detected_genes_1000UMI.mean().rename('rarefied_genes'))
tab=a.set_index('celltype').join(cur)
oldD=pd.read_csv(T/'original_Augur_DMSO_Ctrl.csv');oldD['celltype']=oldD.cell_type.map(ctmap)
tab=tab.join(oldD.set_index('celltype').auc.rename('DMSO_vs_Ctrl_AUC'))
tab.to_csv(T/'complexity_celltype_summary.csv')

out=[]
def corr(label,x,y,confound_count=0):
 ok=np.isfinite(x)&np.isfinite(y);x=np.asarray(x)[ok];y=np.asarray(y)[ok]
 r,p=stats.pearsonr(x,y);rho,sp=stats.spearmanr(x,y)
 if confound_count:
  df=len(x)-confound_count-2
  p=2*stats.t.sf(abs(r)*np.sqrt(df/(1-r*r)),df)
  sp=np.nan
 out.append([label,len(x),r,p,rho,sp])
for col in ['median_genes','median_UMI','rarefied_genes']:
 corr(col+'_vs_original_IR_AUC',tab[col],tab.auc)
def resid(y,confound):
 x=np.column_stack([np.ones(len(y)),confound]);return y-x@np.linalg.lstsq(x,y,rcond=None)[0]
q=np.log1p(tab.median_UMI.to_numpy())
corr('partial_log_UMI',resid(tab.median_genes.to_numpy(),q),resid(tab.auc.to_numpy(),q),1)
q2=np.column_stack([q,tab.fraction_S,tab.fraction_G2M])
corr('partial_log_UMI_and_cycle_fractions',resid(tab.median_genes.to_numpy(),q2),resid(tab.auc.to_numpy(),q2),3)
for r,d in rep.groupby('replicate'):
 joined=tab[['auc']].join(d.set_index('celltype'))
 corr(f'rarefaction_repeat_{r}',joined.median_detected_genes_1000UMI,joined.auc)
o=pd.DataFrame(out,columns=['comparison','celltypes','Pearson_r','nominal_Pearson_p','Spearman_rho','nominal_Spearman_p'])
# Partial Pearson correlations use df adjusted for the number of covariates.
# Rank correlations of residuals are descriptive, so their nominal p-values are omitted.
# None of these calculations estimates between-animal treatment variation.
o.to_csv(T/'complexity_correlations.csv',index=False)
loo=[]
for ct in tab.index:
 z=tab.drop(ct);loo.append([ct,stats.spearmanr(z.median_genes,z.auc).statistic,stats.spearmanr(z.rarefied_genes,z.auc).statistic])
pd.DataFrame(loo,columns=['omitted_celltype','raw_richness_rho','rarefied_richness_rho']).to_csv(T/'complexity_leave_one_type_out.csv',index=False)

plt.rcParams.update({'font.family':'DejaVu Sans','font.size':9,'svg.fonttype':'none','savefig.dpi':200})
fig,axs=plt.subplots(1,3,figsize=(13,4),layout='constrained')
for ax,col,title in zip(axs,['median_genes','median_UMI','rarefied_genes'],['Detected genes at observed depth','UMI counts','Detected genes at 1,000 UMIs']):
 ax.scatter(tab[col],tab.auc,s=30,c='#416d91')
 for ct in ['NP','HSC-Cycle','NMP','ILCP']:
  offset={'NP':(-18,5),'HSC-Cycle':(4,-12),'NMP':(3,3),'ILCP':(3,3)}[ct]
  ax.annotate(ct,(tab.loc[ct,col],tab.loc[ct,'auc']),xytext=offset,textcoords='offset points',fontsize=7)
 rho=stats.spearmanr(tab[col],tab.auc).statistic
 ax.set(title=f'{title}\nSpearman rho = {rho:.3f}',xlabel=col,ylabel='Original IR versus Ctrl Augur AUC',ylim=(.65,1.02))
for ext in ['png','svg']:fig.savefig(F/f'complexity_depth_check.{ext}',bbox_inches='tight')
plt.close(fig)
print(o.head(5).to_string(index=False),flush=True)
print('Fixed-depth program directions',pd.read_csv(T/'program_direction_summary_1000UMI.csv').to_string(index=False),flush=True)
