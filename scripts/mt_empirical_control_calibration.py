#!/usr/bin/env python3
"""Empirically calibrate evidence for biparental mitochondrial inheritance.

Negative controls are pure/homotypic calibration-reference cells.  The script
uses likelihood-support statistics already emitted by mt_fusion_ratio, builds
matched empirical nulls, reports empirical p-values and Benjamini-Hochberg
q-values, summarizes positive-control sensitivity when positive controls are
provided, and recommends an operational single-parent epsilon from negative
control profile intervals.  q-values are FDR-adjusted p-values, not posterior
probabilities.
"""
from __future__ import annotations
import argparse, os
from datetime import datetime, timezone
import numpy as np
import pandas as pd


def parse_args():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--ratio-tsv', nargs='+', required=True)
    p.add_argument('--control-metadata', default=None,
                   help='Optional TSV keyed by library_id+barcode with control_type NEGATIVE|POSITIVE|TARGET and optional true_parent')
    p.add_argument('--output-prefix', required=True)
    p.add_argument('--min-null-cells', type=int, default=20)
    p.add_argument('--fdr', type=float, default=0.05)
    p.add_argument('--epsilon-min', type=float, default=0.005)
    p.add_argument('--epsilon-max', type=float, default=0.10)
    p.add_argument('--epsilon-step', type=float, default=0.005)
    p.add_argument('--max-negative-both-rate', type=float, default=0.01)
    return p.parse_args()

def bh(p):
    p=np.asarray(p,float); out=np.full(len(p),np.nan)
    ok=np.isfinite(p); vals=p[ok]
    if not len(vals): return out
    order=np.argsort(vals); ranked=vals[order]; n=len(vals)
    q=ranked*n/np.arange(1,n+1); q=np.minimum.accumulate(q[::-1])[::-1]; q=np.clip(q,0,1)
    tmp=np.empty(n); tmp[order]=q; out[np.where(ok)[0]]=tmp; return out

def bin_depth(x):
    x=pd.to_numeric(x,errors='coerce').fillna(0)
    return pd.cut(x,[-np.inf,19,49,99,199,499,np.inf],labels=['<20','20-49','50-99','100-199','200-499','500+']).astype(str)

def main():
    a=parse_args()
    r=pd.concat([pd.read_csv(p,sep='\t',dtype={'barcode':str}) for p in a.ratio_tsv],ignore_index=True)
    req={'library_id','barcode','assay_mode','canonical_parent1','canonical_parent2','delta_ll_parent1_only','delta_ll_parent2_only','ratio_molecules','ratio_sites_used'}
    miss=req-set(r.columns)
    if miss: raise SystemExit('ratio input missing: '+', '.join(sorted(miss)))
    if r.duplicated(['library_id','barcode']).any(): raise SystemExit('duplicate library_id+barcode in ratio inputs')
    if a.control_metadata:
        m=pd.read_csv(a.control_metadata,sep='\t',dtype={'barcode':str})
        if not {'library_id','barcode','control_type'}.issubset(m): raise SystemExit('control metadata requires library_id, barcode, control_type')
        if m.duplicated(['library_id','barcode']).any(): raise SystemExit('duplicate library_id+barcode in control metadata')
        r=r.merge(m,on=['library_id','barcode'],how='left',suffixes=('','_meta'))
    else:
        # Calibration-reference C++ output is self-describing.
        r['control_type']=np.where(r.get('calibration_true_parent',pd.Series(index=r.index,dtype=object)).fillna('NA').astype(str).isin(['NA','','.']), 'TARGET','NEGATIVE')
        if 'calibration_true_parent' in r: r['true_parent']=r['calibration_true_parent']
    r['control_type']=r['control_type'].fillna('TARGET').astype(str).str.upper()
    if 'true_parent' not in r: r['true_parent']=np.nan
    for c in ['delta_ll_parent1_only','delta_ll_parent2_only','ratio_molecules','ratio_sites_used','parent2_profile_ci_low','parent2_profile_ci_high','canonical_parent2_fraction']:
        if c in r: r[c]=pd.to_numeric(r[c],errors='coerce')
    # Evidence statistic: loss of the correct single-parent region for negatives; minimum loss of either single-parent region otherwise.
    d1=r['delta_ll_parent1_only'].astype(float); d2=r['delta_ll_parent2_only'].astype(float)
    true=r['true_parent'].astype(str)
    is_p1=true.eq(r['canonical_parent1'].astype(str)); is_p2=true.eq(r['canonical_parent2'].astype(str))
    r['both_evidence_statistic']=np.where(is_p1,d1,np.where(is_p2,d2,np.minimum(d1,d2)))
    r['depth_bin']=bin_depth(r['ratio_molecules']); r['site_bin']=pd.cut(pd.to_numeric(r['ratio_sites_used'],errors='coerce').fillna(0),[-np.inf,2,4,7,12,np.inf],labels=['<=2','3-4','5-7','8-12','13+']).astype(str)
    neg=r[(r.control_type=='NEGATIVE') & np.isfinite(r.both_evidence_statistic)].copy()
    if neg.empty: raise SystemExit('No finite NEGATIVE controls available for empirical null calibration')
    strata=['assay_mode','library_id','canonical_parent1','canonical_parent2','depth_bin','site_bin']
    # Pre-index increasingly broad fallbacks.
    null_maps=[]
    for cols in [strata,strata[:4],['assay_mode','canonical_parent1','canonical_parent2'],['assay_mode'],[]]:
        if cols:
            mp={k if isinstance(k,tuple) else (k,): g.both_evidence_statistic.to_numpy(float) for k,g in neg.groupby(cols,dropna=False)}
        else: mp={():neg.both_evidence_statistic.to_numpy(float)}
        null_maps.append((cols,mp))
    ps=[]; ns=[]; levels=[]
    for _,row in r.iterrows():
        stat=row.both_evidence_statistic
        if not np.isfinite(stat): ps.append(np.nan); ns.append(0); levels.append('NO_STATISTIC'); continue
        chosen=None; label=None
        for cols,mp in null_maps:
            key=tuple(row[c] for c in cols)
            vals=mp.get(key)
            if vals is not None and (len(vals)>=a.min_null_cells or not cols):
                chosen=vals; label='+'.join(cols) if cols else 'GLOBAL'; break
        if chosen is None: chosen=neg.both_evidence_statistic.to_numpy(float); label='GLOBAL'
        ps.append((1.0+float(np.sum(chosen>=stat)))/(len(chosen)+1.0)); ns.append(len(chosen)); levels.append(label)
    r['empirical_p_both']=ps; r['empirical_q_both']=bh(ps); r['empirical_null_cells']=ns; r['empirical_null_stratum']=levels
    r['empirical_both_significant']=(r.empirical_q_both<=a.fdr).astype('Int64')

    os.makedirs(os.path.dirname(os.path.abspath(a.output_prefix)) or '.',exist_ok=True)
    out_cols=['library_id','barcode','control_type','true_parent','assay_mode','canonical_parent1','canonical_parent2','both_evidence_statistic','empirical_p_both','empirical_q_both','empirical_null_cells','empirical_null_stratum','empirical_both_significant']
    r[out_cols].to_csv(a.output_prefix+'.mt_empirical_control_cells.tsv',sep='\t',index=False,float_format='%.12g')

    summaries=[]
    for typ,g in r.groupby('control_type',sort=True):
        finite=g.empirical_q_both.notna()
        summaries.append({'control_type':typ,'n_cells':len(g),'n_calibrated':int(finite.sum()),
                          'fraction_empirical_both_significant':float((g.loc[finite,'empirical_q_both']<=a.fdr).mean()) if finite.any() else np.nan,
                          'interpretation':'false_positive_rate' if typ=='NEGATIVE' else ('positive_control_sensitivity' if typ=='POSITIVE' else 'target_call_rate')})
    pd.DataFrame(summaries).to_csv(a.output_prefix+'.mt_empirical_control_summary.tsv',sep='\t',index=False,float_format='%.12g')

    # Calibrate epsilon using false BOTH classification in pure-parent controls.
    eps_grid=np.arange(a.epsilon_min,a.epsilon_max+0.5*a.epsilon_step,a.epsilon_step)
    erows=[]; rec=np.nan
    if {'parent2_profile_ci_low','parent2_profile_ci_high'}.issubset(neg):
        ci=neg[np.isfinite(neg.parent2_profile_ci_low)&np.isfinite(neg.parent2_profile_ci_high)]
        for e in eps_grid:
            both=(ci.parent2_profile_ci_low>e)&(ci.parent2_profile_ci_high<1-e)
            rate=float(both.mean()) if len(ci) else np.nan
            erows.append({'single_parent_epsilon':e,'n_negative_controls':len(ci),'false_both_rate':rate,
                          'meets_max_negative_both_rate':int(np.isfinite(rate) and rate<=a.max_negative_both_rate)})
            if not np.isfinite(rec) and np.isfinite(rate) and rate<=a.max_negative_both_rate: rec=float(e)
    pd.DataFrame(erows).to_csv(a.output_prefix+'.mt_epsilon_calibration.tsv',sep='\t',index=False,float_format='%.12g')
    # Descriptive control SD benchmark, explicitly not a biological posterior.
    sdvals=[]
    if 'canonical_parent2_fraction' in neg:
        for _,g in neg.groupby(['assay_mode','library_id','canonical_parent1','canonical_parent2'],dropna=False):
            x=g.canonical_parent2_fraction.dropna().to_numpy(float)
            if len(x)>=5: sdvals.append(float(np.std(x,ddof=1)))
    sd95=float(np.quantile(sdvals,0.95)) if sdvals else np.nan
    pd.DataFrame([
        ('recommended_single_parent_epsilon',rec),('max_negative_both_rate',a.max_negative_both_rate),
        ('descriptive_negative_control_sd_95',sd95),('fdr_threshold',a.fdr),('min_null_cells',a.min_null_cells),
        ('q_value_interpretation','FDR-adjusted empirical p-value; not posterior probability'),
        ('timestamp_utc',datetime.now(timezone.utc).isoformat())
    ],columns=['parameter','value']).to_csv(a.output_prefix+'.mt_empirical_calibration_recommendations.tsv',sep='\t',index=False)

if __name__=='__main__': main()
