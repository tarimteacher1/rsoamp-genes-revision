#!/usr/bin/env python3
"""Versioned, strand-aware upstream motif audit; Python >=3.8, standard library only."""
import argparse, bisect, csv, gzip, hashlib, json, math, re, sys, time
from collections import defaultdict, Counter
from pathlib import Path

LENGTH = 2000
# Exact legacy sequences retained; only double-strand-equivalent TGACG/CGTCA merged.
ELEMENTS = [
 ('ABRE','ACGTG','ABRE','ABRE-associated core'),
 ('G-box','CACGTG','G-box','G-box sequence'),
 ('DRE_CRT','RCCGAC','DRE_CRT','DRE/CRT core consensus'),
 ('MYB_core','CNGTTR','MYB','MYB-associated core consensus'),
 ('MYB_MBS','CAACTG','MYB_MBS','MBS sequence'),
 ('E_box','CANNTG','MYC','E-box consensus; not specific to MYC'),
 ('W-box','TTGACY','W-box','WRKY-associated W-box consensus'),
 ('TC-rich','GTTTTCTTAC','TC-rich','legacy TC-rich sequence'),
 ('LTR','CCGAAA','LTR','legacy LTR sequence'),
 ('GT1GMSCAM4','GAAAAA','GT-1','GT1GMSCAM4 sequence; not generic light GT1'),
 ('TGACG_CGTCA_core','TGACG','as-1_TGACG;CGTCA-motif','single TGACG/CGTCA core; not a complete as-1 element'),
 ('ARE','AAACCA','ARE','legacy ARE sequence'),
 ('P_box','CCTTTTG','P-box_GARE','P-box sequence; not a generic GARE'),
 ('ACGTG_flank_B','BACGTG','ABRE-like','custom legacy ACGTG-core pattern with B flank'),
 ('STRE_like','AGGGG','STRE','STRE-like sequence; yeast functional evidence does not establish plant activity'),
]
IUPAC=dict(zip('ACGTRYSWKMBDHVN',['A','C','G','T','[AG]','[CT]','[GC]','[AT]','[GT]','[AC]','[CGT]','[AGT]','[ACT]','[ACG]','[ACGT]']))
COMP=str.maketrans('ACGTRYSWKMBDHVN','TGCAYRSWMKVHDBN')
def rc(seq):return seq.translate(COMP)[::-1]
def attrs(text):
 return {a.split('=',1)[0]:a.split('=',1)[1] for a in text.strip().strip(';').split(';') if '=' in a}
def read_tsv(path):
 with Path(path).open(encoding='utf-8-sig',newline='') as f:return list(csv.DictReader(f,delimiter='\t'))
def write_tsv(path,rows):
 assert rows, str(path)
 with Path(path).open('w',encoding='utf-8',newline='') as f:
  w=csv.DictWriter(f,fieldnames=list(rows[0]),delimiter='\t');w.writeheader();w.writerows(rows)
def sha(path):
 h=hashlib.sha256()
 with Path(path).open('rb') as f:
  for b in iter(lambda:f.read(4*1024*1024),b''):h.update(b)
 return h.hexdigest()
def bh(ps):
 order=sorted(range(len(ps)),key=ps.__getitem__);out=[1.0]*len(ps);v=1.0
 for j in range(len(order)-1,-1,-1):
  i=order[j];v=min(v,ps[i]*len(ps)/(j+1));out[i]=v
 return out
def fisher(a,n,c,m):
 # Integer hypergeometric sum avoids lgamma cancellation for small target cohorts.
 k=a+c;N=n+m;den=math.comb(N,n)
 lo=max(a,0,n-(N-k));hi=min(n,k)
 return sum(math.comb(k,x)*math.comb(N-k,n-x) for x in range(lo,hi+1))/den
def load_gff(path):
 genes={};alias={};parents={};utr=[]
 with Path(path).open() as f:
  for line in f:
   if line.startswith('#'):continue
   p=line.rstrip('\n').split('\t')
   if len(p)!=9:continue
   a=attrs(p[8]);ident=a.get('ID');parent=a.get('Parent','').split(',')
   if p[2]=='gene':
    name=a.get('Name') or a.get('gene_id') or ident
    if not name:continue
    if name in genes:raise ValueError('Duplicate gene key '+name)
    genes[name]={'gene_id':name,'chrom':p[0],'gene_start':int(p[3]),'gene_end':int(p[4]),'strand':p[6]}
    alias[name]=name
    if ident:alias[ident]=name
   elif ident:parents[ident]=parent
   if p[2].lower()=='five_prime_utr':utr.append(parent)
 def resolve(p,seen=None):
  if p in alias:return {alias[p]}
  seen=set() if seen is None else seen
  if p in seen:return set()
  seen.add(p);out=set()
  for q in parents.get(p,[]):out.update(resolve(q,seen.copy()))
  return out
 utr_count=Counter();unlinked=0
 for ps in utr:
  matched=set()
  for p in ps:matched.update(resolve(p))
  if not matched:unlinked+=1
  for g in matched:utr_count[g]+=1
 return genes,utr_count,unlinked
def indices(genes):
 by=defaultdict(list)
 for k,g in genes.items():by[g['chrom']].append((g['gene_start'],g['gene_end'],k))
 out={}
 for chrom,rows in by.items():
  rows.sort();mx=0;maximum=[]
  for _,e,_ in rows:mx=max(mx,e);maximum.append(mx)
  out[chrom]=(rows,[r[0] for r in rows],maximum)
 return out
def overlapping(idx,chrom,s,e,ident):
 rows,starts,maxima=idx[chrom];j=bisect.bisect_right(starts,e)-1;found=[]
 while j>=0 and maxima[j]>=s+1:
  a,b,k=rows[j]
  if k!=ident and a<=e and b>=s+1:found.append(k)
  j-=1
 return sorted(set(found))
def fasta_contigs(path,hasher):
 with Path(path).open('rb') as f:
  name=None;pieces=[]
  for line in f:
   hasher.update(line)
   if line.startswith(b'>'):
    if name is not None:yield name,b''.join(pieces).decode('ascii').upper()
    name=line[1:].split()[0].decode();pieces=[]
   else:pieces.append(line.strip())
  if name is not None:yield name,b''.join(pieces).decode('ascii').upper()
def main():
 p=argparse.ArgumentParser(description=__doc__)
 for k in ['genome','gff','catalogue','outdir']:p.add_argument('--'+k,required=True,type=Path)
 p.add_argument('--legacy-results',type=Path);args=p.parse_args()
 args.outdir.mkdir(parents=True,exist_ok=False);t=time.time()
 genes,utrs,unlinked=load_gff(args.gff);idx=indices(genes)
 cat=[r for r in read_tsv(args.catalogue) if r.get('catalogue_status','INCLUDED').startswith('INCLUDED')]
 members={r['gene_id']:r for r in cat};assert len(members)==len(cat)
 assert set(members)<=set(genes),'Catalogue gene missing from GFF'
 pats={name:[re.compile('(?=('+''.join(IUPAC[c] for c in seq)+'))') for seq in {s,rc(s)}] for name,s,_,_ in ELEMENTS}
 bychrom=defaultdict(list)
 for k,g in genes.items():bychrom[g['chrom']].append(k)
 audits=[];counts={};genome_hash=hashlib.sha256();seen=set()
 with gzip.open(args.outdir/'upstream_sequences.fa.gz','wt',encoding='ascii') as fout:
  for chrom,sequence in fasta_contigs(args.genome,genome_hash):
   if chrom not in bychrom:continue
   for ident in sorted(bychrom[chrom],key=lambda k:genes[k]['gene_start']):
    g=genes[ident];seen.add(ident)
    assert g['strand'] in ('+','-'),g
    if g['strand']=='+':s=max(0,g['gene_start']-1-LENGTH);e=g['gene_start']-1
    else:s=g['gene_end'];e=min(len(sequence),s+LENGTH)
    assert 0<=s<=e<=len(sequence),(g,s,e)
    seq=sequence[s:e];seq=rc(seq) if g['strand']=='-' else seq
    ov=overlapping(idx,chrom,s,e,ident) if len(seq) else []
    non=sum(seq.count(b) for b in set(seq)-set('ACGT'))
    count={name:len({(hit.start(),hit.end(1)) for pattern in patterns for hit in pattern.finditer(seq)}) for name,patterns in pats.items()}
    counts[ident]=count;full=len(seq)==LENGTH;clean=full and non==0 and not ov
    row={'member_id':members.get(ident,{}).get('member_id',''),'family':members.get(ident,{}).get('family',''),**g,
     'target':int(ident in members),'promoter_start_1based':s+1 if seq else '', 'promoter_end_1based':e,'promoter_length_bp':len(seq),
     'non_acgt_bases':non,'gc_fraction':(seq.count('G')+seq.count('C'))/len(seq) if seq else '',
     'overlapping_other_gene_count':len(ov),'overlapping_other_gene_ids':';'.join(ov),
     'annotated_five_prime_utr_features':utrs[ident],'sequence_sha256':hashlib.sha256(seq.encode()).hexdigest(),
     'primary_full_length':int(full),'S1_clean':int(clean),'S2_clean_utr':int(clean and utrs[ident]>0)}
    audits.append(row);fout.write('>'+ident+'\n'+seq+'\n')
   print('Scanned '+chrom+'; genes='+str(len(counts)),flush=True)
 assert set(genes)==seen, 'GFF genes on absent FASTA contigs'
 audit_by={r['gene_id']:r for r in audits}
 write_tsv(args.outdir/'all_gene_upstream_audit.tsv',audits)
 write_tsv(args.outdir/'promoter_extraction_audit.tsv',[r for r in audits if r['target']])
 wide=[{'gene_id':r['gene_id'],**counts[r['gene_id']]} for r in audits]
 write_tsv(args.outdir/'all_gene_motif_counts.tsv',wide)
 full_rows=[]
 for k in sorted(members):
  for name,consensus,aliases,category in ELEMENTS:
   n=counts[k][name]
   full_rows.append({'member_id':members[k].get('member_id',''),'gene_id':k,'family':members[k].get('family',''),
    'element':name,'iupac_consensus':consensus,'annotation_category':category,'nonredundant_genomic_site_count':n,'present':'YES' if n else 'NO',
    'interpretation':'Sequence match only; no binding or functional activation inferred','legacy_aliases':aliases})
 write_tsv(args.outdir/'cis_element_full_results.tsv',full_rows)
 enrichment=[];summary=[]
 for scenario in ('primary_full_length','S1_clean','S2_clean_utr'):
  tar=[r['gene_id'] for r in audits if r['target'] and r[scenario]]
  bg=[r['gene_id'] for r in audits if not r['target'] and r[scenario]]
  assert tar and bg,(scenario,len(tar),len(bg))
  rows=[]
  for name,consensus,aliases,category in ELEMENTS:
   a=sum(counts[k][name]>0 for k in tar);c=sum(counts[k][name]>0 for k in bg);n=len(tar);m=len(bg);b=n-a;d=m-c
   odds=(a*d)/(b*c) if b*c else ('inf' if a*d else 'undefined')
   rows.append({'scenario':scenario,'element':name,'iupac_consensus':consensus,'annotation_category':category,
    'target_full_length_promoters':n,'target_with_motif':a,'target_fraction':a/n,
    'non_target_background_full_length_promoters':m,'background_with_motif':c,'background_fraction':c/m,
    'prevalence_difference':a/n-c/m,'odds_ratio':odds,'fisher_greater_p':fisher(a,n,c,m),
    'bh_adjusted_p':'','bh_all_45_p':'','enriched_at_bh_fdr_0.05':''})
  for row,q in zip(rows,bh([r['fisher_greater_p'] for r in rows])):
   row['bh_adjusted_p']=q;row['enriched_at_bh_fdr_0.05']='YES' if q<.05 else 'NO'
  enrichment.extend(rows)
  summary.append({'scenario':scenario,'target_n':n,'background_n':m,'minimum_p':min(r['fisher_greater_p'] for r in rows),
   'minimum_bh_q':min(r['bh_adjusted_p'] for r in rows),'significant':[r['element'] for r in rows if r['bh_adjusted_p']<.05],
   'mean_target_gc':sum(audit_by[k]['gc_fraction'] for k in tar)/n,'mean_background_gc':sum(audit_by[k]['gc_fraction'] for k in bg)/m})
 for row,q in zip(enrichment,bh([r['fisher_greater_p'] for r in enrichment])):row['bh_all_45_p']=q
 write_tsv(args.outdir/'cis_element_enrichment.tsv',[r for r in enrichment if r['scenario']=='primary_full_length'])
 write_tsv(args.outdir/'cis_element_sensitivity.tsv',enrichment)
 if args.legacy_results:
  legacy_counts=read_tsv(args.legacy_results/'cis_element_full_results.tsv');legacy_enrich=read_tsv(args.legacy_results/'cis_element_enrichment.tsv')
  alias_map={alias:name for name,seq,aliases,category in ELEMENTS for alias in aliases.split(';')}
  assert all(counts[r['gene_id']][alias_map[r['element']]]==int(r['nonredundant_genomic_site_count']) for r in legacy_counts)
  primary={r['element']:r for r in enrichment if r['scenario']=='primary_full_length'}
  for r in legacy_enrich:
   new=primary[alias_map[r['element']]]
   for k in ('target_full_length_promoters','target_with_motif','non_target_background_full_length_promoters','background_with_motif'):
    assert int(r[k])==new[k],(r['element'],k,r[k],new[k])
  legacy_check={'target_entries_checked':len(legacy_counts),'background_patterns_checked':len(legacy_enrich),'counts_exact_match':True}
 else:legacy_check={'not_requested':True}
 result={'python':sys.version,'genome_sha256':genome_hash.hexdigest(),'gff_sha256':sha(args.gff),'catalogue_sha256':sha(args.catalogue),
  'script_sha256':sha(Path(__file__)),'gene_n':len(genes),'target_n':len(members),'patterns':len(ELEMENTS),'unlinked_utr_features':unlinked,
  'all_target_and_background_sequences_rescanned':True,'legacy_comparison':legacy_check,'scenarios':summary,
  'target_no_utr':[r['member_id'] for r in audits if r['target'] and not r['annotated_five_prime_utr_features']],
  'target_overlaps':[r['member_id'] for r in audits if r['target'] and r['overlapping_other_gene_count']],
  'elapsed_seconds':time.time()-t}
 (args.outdir/'summary.json').write_text(json.dumps(result,indent=2)+'\n')
 (args.outdir/'PROMOTER_V2_COMPLETE.PASS').write_text('PASS: all sequences rescanned; legacy counts reproduced; three predefined scenarios.\n')
 print(json.dumps(result,indent=2),flush=True)
if __name__=='__main__':main()
