"""Independent validation: literal expansion/string search plus normalized hypergeometric weights."""
import csv,gzip,hashlib,itertools,json,math,argparse
from pathlib import Path
p=argparse.ArgumentParser();p.add_argument('--results',required=True,type=Path);p.add_argument('--sources',required=True,type=Path);p.add_argument('--out',required=True,type=Path);a=p.parse_args()
def rows(p):
 with p.open(encoding='utf-8-sig') as f:return list(csv.DictReader(f,delimiter='\t'))
sources=rows(a.sources);wide={r['gene_id']:r for r in rows(a.results/'all_gene_motif_counts.tsv')};audits={r['gene_id']:r for r in rows(a.results/'all_gene_upstream_audit.tsv')}
letters={'A':'A','C':'C','G':'G','T':'T','R':'AG','Y':'CT','N':'ACGT','B':'CGT'};comp=str.maketrans('ACGT','TGCA')
expanded={r['element']:{''.join(t) for t in itertools.product(*(letters[b] for b in r['iupac_consensus']))} for r in sources}
for k,seqs in expanded.items():seqs.update({s.translate(comp)[::-1] for s in seqs.copy()})
checked=0;mismatches=[]
with gzip.open(a.results/'upstream_sequences.fa.gz','rt') as f:
 for header in f:
  assert header.startswith('>');ident=header[1:].strip();seq=next(f).strip();row=audits[ident]
  assert hashlib.sha256(seq.encode()).hexdigest()==row['sequence_sha256']
  assert len(seq)==int(row['promoter_length_bp'])
  assert sum(c not in 'ACGT' for c in seq)==int(row['non_acgt_bases'])
  for k,variants in expanded.items():
   sites=set()
   for motif in variants:
    j=seq.find(motif)
    while j!=-1:sites.add((j,len(motif)));j=seq.find(motif,j+1)
   if len(sites)!=int(wide[ident][k]):mismatches.append((ident,k,len(sites),wide[ident][k]))
   checked+=1
assert not mismatches,mismatches[:10]
stats=rows(a.results/'cis_element_sensitivity.tsv');maxerr=0
for r in stats:
 key=r['scenario'];t=[g for g,d in audits.items() if d['target']=='1' and d[key]=='1'];b=[g for g,d in audits.items() if d['target']=='0' and d[key]=='1']
 x=sum(int(wide[g][r['element']])>0 for g in t);y=sum(int(wide[g][r['element']])>0 for g in b);n=len(t);m=len(b)
 assert (x,y,n,m)==tuple(int(r[k]) for k in ['target_with_motif','background_with_motif','target_full_length_promoters','non_target_background_full_length_promoters'])
 K=x+y;N=n+m;lo=max(0,n-(N-K));hi=min(n,K)
 logs=[math.lgamma(K+1)-math.lgamma(k+1)-math.lgamma(K-k+1)+math.lgamma(N-K+1)-math.lgamma(n-k+1)-math.lgamma(N-K-n+k+1) for k in range(lo,hi+1)]
 shift=max(logs);weights=[math.exp(v-shift) for v in logs];pvalue=sum(w for k,w in zip(range(lo,hi+1),weights) if k>=x)/sum(weights)
 maxerr=max(maxerr,abs(pvalue-float(r['fisher_greater_p'])))
for scenario in {r['scenario'] for r in stats}:
 s=sorted([r for r in stats if r['scenario']==scenario],key=lambda r:float(r['fisher_greater_p']))
 for i,r in enumerate(s):
  q=min(1,min(float(other['fisher_greater_p'])*len(s)/(j+1) for j,other in enumerate(s) if j>=i))
  assert abs(q-float(r['bh_adjusted_p']))<1e-12
assert maxerr<1e-8,maxerr
# Cohort flags are independently reconstructed from reported length/ambiguity/overlap/UTR fields.
for r in audits.values():
 full=int(r['promoter_length_bp'])==2000;clean=full and int(r['non_acgt_bases'])==0 and int(r['overlapping_other_gene_count'])==0
 assert (int(r['primary_full_length']),int(r['S1_clean']),int(r['S2_clean_utr']))==(int(full),int(clean),int(clean and int(r['annotated_five_prime_utr_features'])>0))
result={'all_genes':len(audits),'literal_string_search_counts_checked':checked,'count_mismatches':0,'all_sequence_hashes_checked':True,'scenario_tests_checked':len(stats),'maximum_independent_fisher_p_error':maxerr,'BH_checked':True,'cohort_rules_checked':True,'status':'PASS'}
a.out.write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2))
