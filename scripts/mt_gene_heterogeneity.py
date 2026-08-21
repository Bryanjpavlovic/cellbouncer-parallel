#!/usr/bin/env python3
"""RNA-only gene-level mitochondrial ratio heterogeneity diagnostic.

This is a nonblocking post-hoc QC layer.  It maps used mt sites to genes, fits
gene-specific parental fractions while holding the cell's ambient fraction and
rho at the primary fit, and calibrates the gene-vs-common likelihood gain by a
deterministic parametric bootstrap.  It never edits or rejects the primary
mt_fusion_ratio estimate.
"""
from __future__ import annotations
import argparse, os, math
from datetime import datetime, timezone
import numpy as np
import pandas as pd
from scipy.special import betaln, gammaln

def parse_args():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--ratio-tsv', nargs='+', required=True); p.add_argument('--site-tsv', nargs='+', required=True)
    p.add_argument('--site-gene-tsv', required=True, help='TSV with chrom,pos,gene')
    p.add_argument('--output-prefix', required=True); p.add_argument('--bootstrap-replicates',type=int,default=200)
    p.add_argument('--seed',type=int,default=1); p.add_argument('--alpha',type=float,default=0.05); p.add_argument('--min-genes',type=int,default=2)
    p.add_argument('--grid-step',type=float,default=0.01)
    p.add_argument('--annotated-ratio-out', default=None, help='Optional annotated copy of ratio input; original inputs are never modified')
    return p.parse_args()

def binom_ll(k,n,q):
    q=np.clip(q,1e-12,1-1e-12); return gammaln(n+1)-gammaln(k+1)-gammaln(n-k+1)+k*np.log(q)+(n-k)*np.log1p(-q)
def bb_ll(k,n,q,rho):
    if rho<=1e-12 or n<=1: return binom_ll(k,n,q)
    t=(1-rho)/rho; a=np.clip(q,1e-12,1-1e-12)*t; b=(1-np.clip(q,1e-12,1-1e-12))*t
    return gammaln(n+1)-gammaln(k+1)-gammaln(n-k+1)+betaln(k+a,n-k+b)-betaln(a,b)
def q_for(df,r,c):
    p=(1-r)*df.parent1_alt_probability.to_numpy(float)+r*df.parent2_alt_probability.to_numpy(float)
    amb=df.ambient_alt_fraction.to_numpy(float) if 'ambient_alt_fraction' in df else np.full(len(df),0.5)
    amb=np.where(np.isfinite(amb),amb,0.5); return np.clip((1-c)*p+c*amb,1e-12,1-1e-12)
def ll_grid(df,c,rho,grid):
    k=df.alt_molecules.to_numpy(int); n=(df.ref_molecules+df.alt_molecules).to_numpy(int)
    beta=bool(rho>0)
    vals=[]
    for r in grid:
        q=q_for(df,r,c)
        vals.append(float(np.sum([bb_ll(kk,nn,qq,rho) if beta else binom_ll(kk,nn,qq) for kk,nn,qq in zip(k,n,q)])))
    return np.asarray(vals)
def simulate_counts(df,r,c,rho,rng):
    n=(df.ref_molecules+df.alt_molecules).to_numpy(int); q=q_for(df,r,c); out=[]
    for nn,qq in zip(n,q):
        if rho>1e-12 and nn>1:
            t=(1-rho)/rho; pp=rng.beta(qq*t,(1-qq)*t); out.append(rng.binomial(nn,pp))
        else: out.append(rng.binomial(nn,qq))
    return np.asarray(out,int)
def main():
    a=parse_args(); rng=np.random.default_rng(a.seed); grid=np.arange(0,1+0.5*a.grid_step,a.grid_step)
    ratio=pd.concat([pd.read_csv(p,sep='\t',dtype={'barcode':str}) for p in a.ratio_tsv],ignore_index=True)
    sites=pd.concat([pd.read_csv(p,sep='\t',dtype={'barcode':str}) for p in a.site_tsv],ignore_index=True)
    genes=pd.read_csv(a.site_gene_tsv,sep='\t')
    if not {'chrom','pos','gene'}.issubset(genes): raise SystemExit('site-gene TSV requires chrom,pos,gene')
    if genes.duplicated(['chrom','pos']).any(): raise SystemExit('duplicate chrom+pos in site-gene TSV')
    sites=sites.merge(genes,on=['chrom','pos'],how='left')
    for c in ['ref_molecules','alt_molecules','parent1_alt_probability','parent2_alt_probability','ambient_alt_fraction']:
        sites[c]=pd.to_numeric(sites[c],errors='coerce')
    ratio_idx=ratio.set_index(['library_id','barcode'],drop=False)
    rows=[]; grows=[]
    for key,g in sites.groupby(['library_id','barcode'],sort=True):
        if key not in ratio_idx.index: continue
        rr=ratio_idx.loc[key]
        if isinstance(rr,pd.DataFrame): raise SystemExit(f'duplicate ratio key {key}')
        if str(rr.get('assay_mode','RNA'))!='RNA':
            rows.append({'library_id':key[0],'barcode':key[1],'expression_heterogeneity_status':'NOT_RNA','n_genes':0}); continue
        gg=g[(g.get('used_in_fit',1).astype(str).isin(['1','1.0','True','true'])) & g.gene.notna()].copy()
        gg=gg[(gg.ref_molecules+gg.alt_molecules)>0]
        gene_groups=[(name,x) for name,x in gg.groupby('gene') if len(x)>0]
        if len(gene_groups)<a.min_genes:
            rows.append({'library_id':key[0],'barcode':key[1],'expression_heterogeneity_status':'LOW_GENE_COVERAGE','n_genes':len(gene_groups)}); continue
        r0=float(rr.get('parent2_fraction',np.nan)); c=float(rr.get('ambient_fraction',0) if np.isfinite(pd.to_numeric(rr.get('ambient_fraction',0),errors='coerce')) else 0)
        rho=float(rr.get('overdispersion_rho',0) if np.isfinite(pd.to_numeric(rr.get('overdispersion_rho',0),errors='coerce')) else 0)
        if not np.isfinite(r0): continue
        common=0.; alt=0.; fits=[]
        for gene,x in gene_groups:
            vals=ll_grid(x,c,rho,grid); j=int(np.argmax(vals)); rhat=float(grid[j]); best=float(vals[j]); base=float(ll_grid(x,c,rho,np.array([r0]))[0])
            common+=base; alt+=best; fits.append((gene,x,rhat,best,base))
        stat=max(0.,2*(alt-common)); boot=[]
        for _ in range(a.bootstrap_replicates):
            bcommon=0.; balt=0.
            for gene,x,_,_,_ in fits:
                xb=x.copy(); xb['alt_molecules']=simulate_counts(x,r0,c,rho,rng); xb['ref_molecules']=(x.ref_molecules+x.alt_molecules).to_numpy(int)-xb.alt_molecules.to_numpy(int)
                vals=ll_grid(xb,c,rho,grid); balt+=float(vals.max()); bcommon+=float(ll_grid(xb,c,rho,np.array([r0]))[0])
            boot.append(max(0.,2*(balt-bcommon)))
        p=(1+sum(x>=stat for x in boot))/(len(boot)+1) if boot else np.nan
        status='EXPRESSION_HETEROGENEOUS' if np.isfinite(p) and p<=a.alpha else 'EXPRESSION_HOMOGENEOUS'
        rows.append({'library_id':key[0],'barcode':key[1],'expression_heterogeneity_status':status,'n_genes':len(fits),'heterogeneity_statistic':stat,'bootstrap_p_value':p,'bootstrap_replicates':a.bootstrap_replicates,'primary_parent2_fraction':r0})
        for gene,x,rhat,best,base in fits:
            grows.append({'library_id':key[0],'barcode':key[1],'gene':gene,'n_sites':len(x),'molecules':int((x.ref_molecules+x.alt_molecules).sum()),'gene_parent2_fraction':rhat,'gene_log_likelihood':best,'common_ratio_log_likelihood':base})
    os.makedirs(os.path.dirname(os.path.abspath(a.output_prefix)) or '.',exist_ok=True)
    cell_df=pd.DataFrame(rows)
    cell_df.to_csv(a.output_prefix+'.mt_gene_heterogeneity_cells.tsv',sep='\t',index=False,float_format='%.12g')
    pd.DataFrame(grows).to_csv(a.output_prefix+'.mt_gene_heterogeneity_genes.tsv',sep='\t',index=False,float_format='%.12g')
    annotated_out = a.annotated_ratio_out or (a.output_prefix + '.mt_ratio_gene_annotated.tsv')
    annotation_cols=[c for c in ['library_id','barcode','expression_heterogeneity_status','n_genes','heterogeneity_statistic','bootstrap_p_value'] if c in cell_df.columns]
    annotated=ratio.merge(cell_df[annotation_cols],on=['library_id','barcode'],how='left') if annotation_cols else ratio.copy()
    annotated.to_csv(annotated_out,sep='\t',index=False,na_rep='NA')
    pd.DataFrame([('bootstrap_replicates',a.bootstrap_replicates),('seed',a.seed),('alpha',a.alpha),('grid_step',a.grid_step),('timestamp_utc',datetime.now(timezone.utc).isoformat())],columns=['parameter','value']).to_csv(a.output_prefix+'.mt_gene_heterogeneity_run_parameters.tsv',sep='\t',index=False)
if __name__=='__main__': main()
