"""Draw sequence-count heatmap; fixed family/member ordering, no clustering or z-score."""
from pathlib import Path
import argparse,csv,json,re
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import Patch
from audit_promoters_v2 import ELEMENTS,read_tsv
p=argparse.ArgumentParser();p.add_argument('--results',required=True,type=Path);p.add_argument('--outdir',required=True,type=Path);a=p.parse_args();a.outdir.mkdir(parents=True,exist_ok=True)
rows=read_tsv(a.results/'cis_element_full_results.tsv');audit={r['member_id']:r for r in read_tsv(a.results/'promoter_extraction_audit.tsv')}
family_order={'Defensin':0,'Snakin_GASA':1,'nsLTP':2}
family={r['member_id']:r['family'] for r in rows}
members=sorted(family,key=lambda m:(family_order.get(family[m],9),int(re.search(r'(\d+)$',m).group())))
names=[r[0] for r in ELEMENTS];lookup={(r['member_id'],r['element']):int(r['nonredundant_genomic_site_count']) for r in rows}
matrix=np.array([[lookup[m,n] for n in names] for m in members]);log=np.log1p(matrix)
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':9,'pdf.fonttype':42,'ps.fonttype':42,'svg.fonttype':'none'})
fig=plt.figure(figsize=(10.4,13.8));ax=fig.add_axes([.14,.14,.66,.77])
im=ax.imshow(log,aspect='auto',cmap='cividis',vmin=0,vmax=float(log.max()),interpolation='nearest')
labels={'MYB_core':'MYB core','E_box':'E-box','GT1GMSCAM4':'GT1GMSCAM4','TGACG_CGTCA_core':'TGACG/CGTCA core','P_box':'P-box','ACGTG_flank_B':'BACGTG (custom)','STRE_like':'STRE-like'}
ax.set_xticks(range(len(names)));ax.set_xticklabels([labels.get(n,n) for n in names],rotation=55,ha='right',fontsize=8.5)
ax.set_yticks(range(len(members)));ax.set_yticklabels(members,fontsize=8)
ax.tick_params(length=0,pad=5);ax.set_xticks(np.arange(-.5,len(names),1),minor=True);ax.set_yticks(np.arange(-.5,len(members),1),minor=True);ax.grid(which='minor',color='white',linewidth=.35);ax.tick_params(which='minor',length=0)
for spine in ax.spines.values():spine.set_visible(False)
for i in range(1,len(members)):
 if family[members[i]]!=family[members[i-1]]:ax.axhline(i-.5,color='white',linewidth=2)
qc=fig.add_axes([.815,.14,.044,.77]);flags=np.array([[int(audit[m]['overlapping_other_gene_count'])>0,int(audit[m]['annotated_five_prime_utr_features'])==0] for m in members])
qc.imshow(flags,aspect='auto',cmap=matplotlib.colors.ListedColormap(['#f2f4f5','#a04d23']),vmin=0,vmax=1,interpolation='nearest');qc.set_yticks([]);qc.set_xticks([0,1]);qc.set_xticklabels(['Gene overlap','No 5′UTR annotation'],rotation=55,ha='right',fontsize=7.5);qc.tick_params(length=0)
for s in qc.spines.values():s.set_visible(False)
cax=fig.add_axes([.89,.65,.018,.23]);cbar=fig.colorbar(im,cax=cax);cbar.set_label('ln(1 + nonredundant site count)',fontsize=8);cbar.ax.tick_params(labelsize=8)
fig.text(.14,.95,'Sequence motifs in annotated upstream regions',fontsize=15,weight='bold')
fig.text(.14,.927,'63 primary AMP genes · 2,000 bp · 15 double-stranded patterns',fontsize=10,color='#4b5962')
fig.text(.14,.035,'Cells show sequence matches. No binding or regulatory activity is inferred.\nBrown flags indicate annotation limitations; complete coordinates and cohort membership are provided in the supplement.',fontsize=8,color='#4b5962',linespacing=1.5)
base=a.outdir/'Figure5_promoter_motif_occurrence'
for ext in ('svg','pdf','png'):fig.savefig(base.with_suffix('.'+ext),dpi=400,facecolor='white')
fig.savefig(a.outdir/'Figure5_preview.png',dpi=105,facecolor='white');plt.close(fig)
(a.outdir/'figure5_qa.json').write_text(json.dumps({'rows':len(members),'columns':len(names),'row_order':members,'column_order':names,'transformation':'natural log(1+count)','clustering':False,'z_score':False,'colormap':'cividis','scale_min':0,'scale_max':float(log.max()),'matplotlib_version':matplotlib.__version__,'numpy_version':np.__version__},indent=2))
print('Figure 5: '+str(matrix.shape))
