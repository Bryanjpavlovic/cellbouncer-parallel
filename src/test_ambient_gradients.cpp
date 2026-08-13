#include <cmath>
#include <cstdlib>
#include <functional>
#include <iostream>
#include <map>
#include <random>
#include <string>
#include <vector>
#include <optimML/multivar_ml.h>
#include "ambient_rna_three_ap.h"

using namespace std;

static bool close_enough(double a, double b, double rtol=1e-6, double atol=1e-7){
    double scale = max(fabs(a), fabs(b));
    return fabs(a-b) <= atol + rtol*scale;
}

static double finite_diff(const function<double(const vector<double>&)>& f,
    const vector<double>& x, int idx){
    vector<double> xp=x, xm=x;
    double h=1e-6;
    xp[idx]+=h; xm[idx]-=h;
    return (f(xp)-f(xm))/(2.0*h);
}

static bool test_callback(const string& name,
    double (*ll)(const vector<double>&,const map<string,double>&,const map<string,int>&),
    void (*dll)(const vector<double>&,const map<string,double>&,const map<string,int>&,vector<double>&),
    const vector<double>& x, const map<string,double>& dd, const map<string,int>& di){
    vector<double> grad(x.size(),0.0);
    dll(x,dd,di,grad);
    auto fn=[&](const vector<double>& p){return ll(p,dd,di);};
    bool ok=true;
    for(int j=0;j<(int)x.size();++j){
        double num=finite_diff(fn,x,j);
        if(!close_enough(num,grad[j])){
            cerr << name << " gradient mismatch idx=" << j
                 << " analytic=" << grad[j] << " numerical=" << num << "\n";
            ok=false;
        }
    }
    return ok;
}

static bool randomized_gradient_tests(){
    mt19937 rng(1729);
    uniform_real_distribution<double> u(0.08,0.92);
    bool ok=true;
    for(int rep=0;rep<10;++rep){
        double n=20.0+rep;
        double k=max(1.0,min(n-1.0,n*u(rng)));
        map<string,int> di;

        map<string,double> d3{{"n",n},{"k",k},{"p_A",u(rng)},{"p_B",u(rng)},{"p_c",u(rng)}};
        ok &= test_callback("ll_three",ll_three,dll_three,{u(rng),u(rng)},d3,di);

        map<string,double> da{{"n",n},{"k",k},{"p_e",u(rng)},{"c",u(rng)}};
        map<string,int> dia{{"ef_idx",1}};
        vector<double> pa{u(rng),u(rng),u(rng)};
        ok &= test_callback("ll_ambmu",ll_ambmu,dll_ambmu,pa,da,dia);

        map<string,double> de{{"n",n},{"k",k},{"p_e",u(rng)},{"p_c",u(rng)},{"c",u(rng)}};
        ok &= test_callback("ll_err_rates",ll_err_rates,dll_err_rates,{0.01+0.05*u(rng),0.01+0.05*u(rng)},de,di);

        map<string,double> dpf{{"n",n},{"k",k},{"p_e",u(rng)}};
        ok &= test_callback("ll_amb_prof_mixture_fitted_c",ll_amb_prof_mixture,
            dll_amb_prof_mixture,{u(rng),u(rng)},dpf,di);

        map<string,double> dpx{{"n",n},{"k",k},{"p_e",u(rng)},{"c",u(rng)}};
        ok &= test_callback("ll_amb_prof_mixture_fixed_c",ll_amb_prof_mixture,
            dll_amb_prof_mixture,{u(rng)},dpx,di);
    }
    return ok;
}

static bool one_hot_source_mapping_test(){
    vector<double> params{0.55}; // fitted contamination precedes mixture block
    optimML::multivar_ml_solver solver(params,ll_amb_prof_mixture,dll_amb_prof_mixture);
    vector<vector<double>> mix{{0.90,0.10},{0.10,0.90},{0.80,0.20},{0.20,0.80}};
    vector<double> start{0.50,0.50};
    vector<double> n(4,200.0);
    vector<double> pe(4,0.50);
    const double c_true=0.60;
    vector<double> source0{0.90,0.10,0.80,0.20};
    vector<double> k;
    for(double p:source0) k.push_back(200.0*((1.0-c_true)*0.50+c_true*p));
    if(!solver.add_mixcomp(mix) || !solver.add_mixcomp_fracs(start)){
        cerr << "failed to initialize one-hot mixture test\n";
        return false;
    }
    solver.add_data("n",n);
    solver.add_data("k",k);
    solver.add_data("p_e",pe);
    solver.constrain_01(0);
    solver.set_silent(true);
    bool solved=false;
    try { solved=solver.solve(); } catch (...) { solved=false; }
    if(!solved || solver.results_mixcomp.size()!=2 || solver.results.empty()){
        cerr << "one-hot mixture solve failed\n";
        return false;
    }
    if(solver.results_mixcomp[0] < 0.95 || solver.results_mixcomp[1] > 0.05){
        cerr << "source coefficient/label mapping failed: "
             << solver.results_mixcomp[0] << " " << solver.results_mixcomp[1] << "\n";
        return false;
    }
    if(fabs(solver.results[0]-c_true)>0.10){
        cerr << "fitted-c mixture offset test failed: c=" << solver.results[0] << "\n";
        return false;
    }
    return true;
}


static bool source_exclusion_profile_tests(){
    bool ok=true;
    const map<int,double> global{{0,0.50},{1,0.30},{2,0.20}};
    const set<int> excluded{0,2};

    const auto zero=apply_source_exclusion_profile(global,excluded,0.0);
    if(zero.status!="global_profile_no_exclusion" || zero.scoring_profile!=global ||
       !close_enough(zero.scoring_profile_mass_sum,1.0) ||
       !close_enough(zero.effective_removed_source_mass_total,0.0)){
        cerr << "lambda=0 source-exclusion endpoint failed\n";
        ok=false;
    }

    const auto hard=apply_source_exclusion_profile(global,excluded,1.0);
    if(hard.status!="hard_source_exclusion" ||
       !close_enough(hard.scoring_profile.at(0),0.0) ||
       !close_enough(hard.scoring_profile.at(1),1.0) ||
       !close_enough(hard.scoring_profile.at(2),0.0) ||
       !close_enough(hard.scoring_profile_mass_sum,1.0) ||
       !close_enough(hard.excluded_source_mass_raw_total,0.70) ||
       !close_enough(hard.effective_removed_source_mass_total,0.70)){
        cerr << "lambda=1 hard source-exclusion endpoint failed\n";
        ok=false;
    }

    for(double lambda: {0.25,0.50,0.75}){
        const auto partial=apply_source_exclusion_profile(global,excluded,lambda);
        double sum=0.0;
        bool finite=true;
        for(const auto& kv: partial.scoring_profile){
            finite = finite && std::isfinite(kv.second) && kv.second >= 0.0;
            sum += kv.second;
        }
        const double expected_removed=lambda*0.70;
        if(partial.status!="partial_source_exclusion" || !finite ||
           !close_enough(sum,1.0) ||
           !close_enough(partial.scoring_profile_mass_sum,1.0) ||
           !close_enough(partial.effective_removed_source_mass_total,expected_removed) ||
           !close_enough(partial.scoring_profile_renormalization_denominator,1.0-expected_removed)){
            cerr << "interior source exclusion failed at lambda=" << lambda << "\n";
            ok=false;
        }
    }

    const auto fallback=apply_source_exclusion_profile(global,{0,1,2},1.0);
    if(fallback.status!="degenerate_all_mass_excluded_fallback_global" ||
       fallback.scoring_profile!=global ||
       !close_enough(fallback.scoring_profile_mass_sum,1.0) ||
       !close_enough(fallback.effective_removed_source_mass_total,0.0)){
        cerr << "all-mass-excluded fallback failed\n";
        ok=false;
    }

    for(double invalid: {-0.01,1.01,std::numeric_limits<double>::quiet_NaN()}){
        const auto assessment=apply_source_exclusion_profile(global,excluded,invalid);
        if(assessment.status!="invalid_source_exclusion_strength" ||
           assessment.scoring_profile!=global){
            cerr << "invalid source exclusion strength was not rejected\n";
            ok=false;
        }
    }
    return ok;
}

int main(){
    bool ok=randomized_gradient_tests();
    ok &= one_hot_source_mapping_test();
    ok &= source_exclusion_profile_tests();
    if(!ok) return 1;
    cout << "ambient mathematical tests passed\n";
    return 0;
}
