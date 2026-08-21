#!/usr/bin/env python3
"""Summarize mt_fusion_ratio diagnostics and optionally learn pooled-rho references."""
from __future__ import annotations
import argparse
import os
from datetime import datetime, timezone
import numpy as np
import pandas as pd
from scipy.stats import spearmanr

def args_parser():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--ratio-tsv', nargs='+', required=True)
    p.add_argument('--site-tsv', nargs='*', default=[])
    p.add_argument('--output-prefix', required=True)
    p.add_argument('--min-rho-cells', type=int, default=20)
    p.add_argument('--min-ratio-molecules', type=int, default=50)
    p.add_argument('--rho-boundary-fraction', type=float, default=0.99)
    p.add_argument('--overdispersion-max', type=float, default=0.25)
    return p.parse_args()

def load(paths):
    return pd.concat([pd.read_csv(p,sep='\t',dtype={'barcode':str}) for p in paths], ignore_index=True)

def main():
    a=args_parser(); r=load(a.ratio_tsv)
    req={'library_id','barcode','canonical_parent1','canonical_parent2','assay_mode','overdispersion_rho','ratio_molecules','ratio_sites_used','status'}
    miss=req-set(r.columns)
    if miss: raise SystemExit('ratio input missing: '+', '.join(sorted(miss)))
    for c in ['overdispersion_rho','overdispersion_se','ratio_molecules','ratio_sites_used','information_condition','raw_parent1_fraction','canonical_parent1_fraction','max_abs_site_influence','fraction_support_from_top_site']:
        if c in r: r[c]=pd.to_numeric(r[c], errors='coerce')
    os.makedirs(os.path.dirname(os.path.abspath(a.output_prefix)) or '.',exist_ok=True)
    # Status counts
    sc=(r.groupby(['assay_mode','status'],dropna=False).size().rename('n_cells').reset_index())
    sc.to_csv(a.output_prefix+'.mt_diagnostic_status.tsv',sep='\t',index=False)
    # Summary by library/pair.
    groups=[]; rho_refs=[]
    gb=['assay_mode','library_id','canonical_parent1','canonical_parent2']
    for key,g in r.groupby(gb,dropna=False,sort=True):
        rho=g['overdispersion_rho'].replace([np.inf,-np.inf],np.nan)
        hi=g[(g['ratio_molecules']>=a.min_ratio_molecules)&rho.notna()].copy()
        corr=np.nan; corr_p=np.nan
        if len(hi)>=max(5,a.min_rho_cells//2) and hi['ratio_molecules'].nunique()>2:
            corr,corr_p=spearmanr(hi['ratio_molecules'],hi['overdispersion_rho'],nan_policy='omit')
        row=dict(zip(gb,key)); row.update({
            'n_cells':len(g),'n_pass':int((g['status']=='PASS').sum()),
            'median_ratio_molecules':float(g['ratio_molecules'].median()),
            'median_ratio_sites_used':float(g['ratio_sites_used'].median()),
            'rho_nonfinite_fraction':float(rho.isna().mean()),
            'rho_boundary_fraction':float((rho>=a.rho_boundary_fraction*a.overdispersion_max).mean()),
            'median_rho':float(rho.median()) if rho.notna().any() else np.nan,
            'median_rho_se':float(g['overdispersion_se'].median()) if 'overdispersion_se' in g else np.nan,
            'rho_vs_molecules_spearman':corr,'rho_vs_molecules_p':corr_p,
            'coverage_trend_supported':int(np.isfinite(corr) and abs(corr)>=0.3 and corr_p<0.01),
            'boundary_fraction_mle':float(((g['canonical_parent1_fraction']<=1e-8)|(g['canonical_parent1_fraction']>=1-1e-8)).mean()) if 'canonical_parent1_fraction' in g else np.nan,
            'median_information_condition':float(g['information_condition'].median()) if 'information_condition' in g else np.nan,
            'median_abs_site_influence':float(g['max_abs_site_influence'].median()) if 'max_abs_site_influence' in g else np.nan,
            'median_top_site_fraction':float(g['fraction_support_from_top_site'].median()) if 'fraction_support_from_top_site' in g else np.nan,
        })
        groups.append(row)
        if str(key[0])!='ATAC' and len(hi)>=a.min_rho_cells:
            pooled=float(hi['overdispersion_rho'].median())
            if np.isfinite(pooled):
                rho_refs.append({'assay_mode':key[0],'library_id':key[1],'canonical_parent1':key[2],
                                 'canonical_parent2':key[3],'rho':pooled,'n_cells':len(hi),
                                 'rho_reference_status':'PASS'})
    pd.DataFrame(groups).to_csv(a.output_prefix+'.mt_diagnostics.tsv',sep='\t',index=False,float_format='%.12g')
    pd.DataFrame(rho_refs,columns=['assay_mode','library_id','canonical_parent1','canonical_parent2','rho','n_cells','rho_reference_status']).to_csv(
        a.output_prefix+'.mt_rho_reference.tsv',sep='\t',index=False,float_format='%.12g')
    # Site depth n>1 diagnostic.
    if a.site_tsv:
        s=load(a.site_tsv)
        if {'ref_molecules','alt_molecules'}.issubset(s):
            n=pd.to_numeric(s.ref_molecules,errors='coerce').fillna(0)+pd.to_numeric(s.alt_molecules,errors='coerce').fillna(0)
            pd.DataFrame([{'n_site_observations':len(s),'fraction_site_observations_n_gt_1':float((n>1).mean()),
                           'median_site_total':float(n.median())}]).to_csv(a.output_prefix+'.mt_site_depth_diagnostics.tsv',sep='\t',index=False)
    pd.DataFrame([
        ('ratio_tsv',';'.join(a.ratio_tsv)),('site_tsv',';'.join(a.site_tsv)),('min_rho_cells',a.min_rho_cells),
        ('min_ratio_molecules',a.min_ratio_molecules),('overdispersion_max',a.overdispersion_max),('timestamp_utc',datetime.now(timezone.utc).isoformat())
    ],columns=['parameter','value']).to_csv(a.output_prefix+'.mt_diagnostic_run_parameters.tsv',sep='\t',index=False)

if __name__=='__main__': main()
