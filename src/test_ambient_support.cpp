#include <cmath>

double logbinom(double n, double k, double p){
    double ll=k*std::log(p)+(n-k)*std::log(1.0-p);
    if(k<n && k!=0){
        double logn=std::log(n), logk=std::log(k), lognk=std::log(n-k);
        ll += n*logn-k*logk-(n-k)*lognk+0.5*(logn-logk-lognk-std::log(2*M_PI));
    }
    return ll;
}

double adjust_p_err(double p, double e_r, double e_a){
    return p - p*e_a + (1.0-p)*e_r;
}
