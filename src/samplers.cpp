// Rcpp kernels for the R reanalysis. R supplies all random-number generation.
// The target, data augmentation and marginal MH moves match truncated_dp.py.
// No Python is called. Counts are doubles because latent richness can exceed 2^31.
#include <Rcpp.h>
#include <algorithm>
#include <cmath>
#include <vector>
using namespace Rcpp;
using std::vector;
// [[Rcpp::plugins(cpp11)]]

double positive_gamma(double shape, double rate, double lower=0) {
  if (lower<=0) {
    double x=R::rgamma(shape,1/rate);
    if (!(x>0) || !R_finite(x)) stop("Gamma draw outside floating-point range");
    return x;
  }
  if (shape>1 && (shape-1)/rate>=lower) {
    for (int i=0;i<100000;++i) {
      double x=R::rgamma(shape,1/rate); if(x>=lower) return x;
    }
    stop("Truncated Gamma rejection did not finish");
  }
  double logtail=R::pgamma(lower,shape,1/rate,false,true);
  double x=R::qgamma(logtail+std::log(R::runif(0,1)),shape,1/rate,false,true);
  if (!R_finite(x) || x<lower) stop("Invalid truncated Gamma quantile");
  return x;
}
void multinomial_count(double n, const vector<double>& p, vector<double>& out) {
  double total=0;for(double x:p) total+=x;
  if (!(total>0) && n>0) stop("Zero allocation probability");
  double left=n;int H=p.size();
  for(int h=0;h<H-1;++h) {
    double q=total>0?std::min(1.,std::max(0.,p[h]/total)):0;
    double x=(left<=0)?0:(q>=1?left:R::rbinom(left,q));
    out[h]=x;left-=x;total-=p[h];
  }
  out[H-1]=left;
}
double poisson_tail_draw(double rate, int cutoff) {
  if(rate>=cutoff) {
    for(int i=0;i<1000000;++i) {double x=R::rpois(rate);if(x>=cutoff)return x;}
    stop("Poisson rejection did not finish");
  }
  double lp=R::ppois(cutoff-1,rate,false,true)+std::log(R::runif(0,1));
  double x=R::qpois(lp,rate,false,true);
  if(x<cutoff || !R_finite(x)) stop("Invalid censored Poisson quantile");
  return x;
}
void column(double rate,const IntegerVector& j,int cv,double cf,vector<double>& p) {
  for(int z=0;z<j.size();++z)p[z]=R::dpois(j[z],rate,false);
  p[j.size()]=cf>0?R::ppois(cv-1,rate,false,false):0;
}
double loglik(const vector<double>& mu,double tail,double detect,
              const NumericVector& f,double cf,double K,bool truncated) {
  double ans=0;
  for(int z=0;z<f.size();++z) {
    if(!(mu[z]>0))return R_NegInf;
    ans+=f[z]*std::log(mu[z]);
  }
  if(cf>0) {if(!(tail>0))return R_NegInf;ans+=cf*std::log(tail);}
  if(truncated){if(!(detect>0))return R_NegInf;ans-=K*std::log(detect);}
  return ans;
}
double stick_logprior(const vector<double>& w,double alpha) {
  int H=w.size();double rem=w[H-1];
  if(!(rem>0))return R_NegInf;
  double lp=(alpha-1)*std::log(rem);
  for(int h=H-2;h>0;--h){rem+=w[h];lp-=std::log(rem);}
  return lp;
}
double logistic(double x) {
  return x>=0?1/(1+std::exp(-x)):std::exp(x)/(1+std::exp(x));
}

// [[Rcpp::export]]
List dp_mcmc_cpp(IntegerVector j, NumericVector f, double shape, double rate,
 double alpha, int H, int sweeps, int burn, int thin=1, bool truncated=true,
 int cv=0, double cf=0, double lower=0, bool mh=true, int mh_scans=1,
 bool swaps=true, NumericVector thresholds=NumericVector::create()) {
  if(j.size()!=f.size() || H<2 || burn>=sweeps || thin<1 || shape<=0 ||
     rate<=0 || alpha<=0 || lower<0 || (sweeps-burn)%thin)
    stop("Invalid sampler arguments");
  double K=sum(f)+cf;if(K<1)stop("At least one observed species is required");
  int J=j.size(), nk=(sweeps-burn)/thin;
  vector<double> theta(H),w(H),mhc(H),sums(H),alloc(H),pr(H),lp(H),
    q(H),w2(H),col(J+1),mu(J),mu2(J),pre(J),muS(J);
  vector<vector<double>> P(H,vector<double>(J+1));
  for(int h=0;h<H;++h)theta[h]=positive_gamma(shape,rate,lower);
  double rest=1;
  for(int h=0;h<H-1;++h){double v=R::rbeta(1,alpha);w[h]=rest*v;rest*=1-v;}
  w[H-1]=rest;
  NumericVector aout(nk),Nout(nk),kout(nk),lastout(nk);
  NumericMatrix small(nk,thresholds.size());
  double sV=1,sT=.5,sW=1,av=0,nv=0,at=0,nt=0,aw=0,nw=0,
    avall=0,nvall=0,atall=0,ntall=0;
  int keep=0;
  for(int it=0;it<sweeps;++it) {
    if(it%1024==0)checkUserInterrupt();
    std::fill(mhc.begin(),mhc.end(),0);std::fill(sums.begin(),sums.end(),0);
    for(int z=0;z<J;++z) {
      double mx=R_NegInf;
      for(int h=0;h<H;++h){
        lp[h]=(w[h]>0)?std::log(w[h])+j[z]*std::log(theta[h])-theta[h]:R_NegInf;
        mx=std::max(mx,lp[h]);
      }
      for(int h=0;h<H;++h)pr[h]=std::exp(lp[h]-mx);
      multinomial_count(f[z],pr,alloc);
      for(int h=0;h<H;++h){mhc[h]+=alloc[h];sums[h]+=alloc[h]*j[z];}
    }
    if(cf>0){
      for(int h=0;h<H;++h)pr[h]=w[h]*R::ppois(cv-1,theta[h],false,false);
      multinomial_count(cf,pr,alloc);
      for(int h=0;h<H;++h)if(alloc[h]>0){
        mhc[h]+=alloc[h];
        for(double z=0;z<alloc[h];++z)sums[h]+=poisson_tail_draw(theta[h],cv);
      }
    }
    double detect=0;
    for(int h=0;h<H;++h){
      pr[h]=w[h]*std::exp(-theta[h]);detect+=w[h]*(-std::expm1(-theta[h]));
    }
    if(truncated) {
      double f0=R::rnbinom(K,detect);
      if(!R_finite(f0) || f0>9007199254740991.)
        stop("Unseen count exceeds exact integer range; do not truncate this draw");
      multinomial_count(f0,pr,alloc);
      for(int h=0;h<H;++h)mhc[h]+=alloc[h];
    }
    double tail=0;for(double x:mhc)tail+=x;rest=1;
    for(int h=0;h<H-1;++h){
      tail=std::max(0.,tail-mhc[h]);
      double v=R::rbeta(1+mhc[h],alpha+tail);w[h]=rest*v;rest*=1-v;
    }
    w[H-1]=rest;int occupied=0;
    for(int h=0;h<H;++h){
      theta[h]=positive_gamma(shape+sums[h],rate+mhc[h],lower);
      if(mhc[h]>0)++occupied;
    }
    for(int scan=0;scan<(mh?mh_scans:0);++scan){
      std::fill(mu.begin(),mu.end(),0);
      double T=0,Z=0;
      for(int h=0;h<H;++h){
        column(theta[h],j,cv,cf,P[h]);q[h]=-std::expm1(-theta[h]);
        for(int z=0;z<J;++z)mu[z]+=w[h]*P[h][z];
        T+=w[h]*P[h][J];Z+=w[h]*q[h];
      }
      double ll=loglik(mu,T,Z,f,cf,K,truncated);
      std::fill(pre.begin(),pre.end(),0);double preT=0,preZ=0,rem=1;
      for(int h=0;h<H-1;++h){
        if(rem>1e-300 && w[h]>0){
          double V=w[h]/rem;
          if(V<1){
            double Vn=logistic(std::log(V)-std::log1p(-V)+sV*R::rnorm(0,1));
            if(Vn>0 && Vn<1){
              double wn=Vn*rem,rr=(1-Vn)/(1-V);
              for(int z=0;z<J;++z)
                mu2[z]=pre[z]+wn*P[h][z]+rr*(mu[z]-pre[z]-w[h]*P[h][z]);
              double Tn=preT+wn*P[h][J]+rr*(T-preT-w[h]*P[h][J]);
              double Zn=preZ+wn*q[h]+rr*(Z-preZ-w[h]*q[h]);
              double lln=loglik(mu2,Tn,Zn,f,cf,K,truncated);
              double lr=lln-ll+alpha*(std::log1p(-Vn)-std::log1p(-V))+std::log(Vn/V);
              ++nv;++nvall;
              if(std::log(R::runif(0,1))<lr){
                ++av;++avall;w[h]=wn;for(int l=h+1;l<H;++l)w[l]*=rr;
                mu=mu2;T=Tn;Z=Zn;ll=lln;
              }
            }
          }
        }
        for(int z=0;z<J;++z)pre[z]+=w[h]*P[h][z];
        preT+=w[h]*P[h][J];preZ+=w[h]*q[h];rem-=w[h];
      }
      for(int h=0;h<H;++h){
        if(w[h]<=0)continue;
        double thn=theta[h]*std::exp(sT*R::rnorm(0,1));++nt;++ntall;
        if(thn<lower || !(thn>0) || !R_finite(thn))continue;
        column(thn,j,cv,cf,col);double qn=-std::expm1(-thn);
        for(int z=0;z<J;++z)mu2[z]=mu[z]+w[h]*(col[z]-P[h][z]);
        double Tn=T+w[h]*(col[J]-P[h][J]),Zn=Z+w[h]*(qn-q[h]);
        double lln=loglik(mu2,Tn,Zn,f,cf,K,truncated);
        double lr=lln-ll+shape*std::log(thn/theta[h])-rate*(thn-theta[h]);
        if(std::log(R::runif(0,1))<lr){
          ++at;++atall;theta[h]=thn;q[h]=qn;P[h]=col;
          mu=mu2;T=Tn;Z=Zn;ll=lln;
        }
      }
      double prior=stick_logprior(w,alpha);
      for(int rep=0;rep<3;++rep){
        if(!R_finite(prior))break;
        double cutoff=std::exp(std::log(1e-3/K)*R::runif(0,1));
        double Ws=0,Wo=0,TS=0,ZS=0;int ns=0;std::fill(muS.begin(),muS.end(),0);
        for(int h=0;h<H;++h){
          if(theta[h]<cutoff){
            Ws+=w[h];++ns;for(int z=0;z<J;++z)muS[z]+=w[h]*P[h][z];
            TS+=w[h]*P[h][J];ZS+=w[h]*q[h];
          }else Wo+=w[h];
        }
        if(ns==0 || ns==H || Ws<=0 || Wo<=0)continue;
        double z=std::log(Ws/Wo)+sW*R::rnorm(0,1);
        double Wsn=logistic(z),Won=logistic(-z);
        double r1=Wsn/Ws,r2=Won/Wo;
        for(int h=0;h<H;++h)w2[h]=w[h]*(theta[h]<cutoff?r1:r2);
        double prior2=stick_logprior(w2,alpha);++nw;
        if(!R_finite(prior2) || !(r1>0) || !(r2>0))continue;
        for(int z=0;z<J;++z)mu2[z]=r1*muS[z]+r2*(mu[z]-muS[z]);
        double Tn=r1*TS+r2*(T-TS),Zn=r1*ZS+r2*(Z-ZS);
        double lln=loglik(mu2,Tn,Zn,f,cf,K,truncated);
        double lr=lln-ll+prior2-prior+ns*std::log(r1)+(H-ns)*std::log(r2);
        if(std::log(R::runif(0,1))<lr){
          ++aw;w=w2;mu=mu2;T=Tn;Z=Zn;ll=lln;prior=prior2;
        }
      }
      if(swaps)for(int rep=0;rep<H;++rep){
        int h1=std::floor(R::runif(0,H)),h2=std::floor(R::runif(0,H));
        if(h1==h2)continue;
        w2=w;std::swap(w2[h1],w2[h2]);double prior2=stick_logprior(w2,alpha);
        if(R_finite(prior2) && std::log(R::runif(0,1))<prior2-prior){
          w=w2;std::swap(theta[h1],theta[h2]);std::swap(q[h1],q[h2]);
          std::swap(P[h1],P[h2]);prior=prior2;
        }
      }
      if(it<burn && (it+1)%200==0){
        if(nv>0)sV*=std::exp(av/nv-.3);if(nt>0)sT*=std::exp(at/nt-.3);
        if(nw>0)sW*=std::exp(aw/nw-.3);
        av=nv=at=nt=aw=nw=0;
      }
    }
    if(it>=burn && (it-burn)%thin==0){
      double Z=0;for(int h=0;h<H;++h)Z+=w[h]*(-std::expm1(-theta[h]));
      aout[keep]=1-Z;kout[keep]=occupied;lastout[keep]=w[H-1];
      Nout[keep]=truncated?K+R::rnbinom(K,Z):K;
      if(!R_finite(Nout[keep]))stop("Nonfinite richness draw");
      for(int b=0;b<thresholds.size();++b){
        double mass=0;for(int h=0;h<H;++h)if(theta[h]<=thresholds[b])mass+=w[h];
        small(keep,b)=mass;
      }
      ++keep;
    }
  }
  return List::create(_["a"]=aout,_["N"]=Nout,_["clusters"]=kout,_["small"]=small,
    _["last_weight"]=lastout,_["accept_v"]=avall/std::max(1.,nvall),
    _["accept_rate"]=atall/std::max(1.,ntall));
}

double predictive(int n,double size,double total,double shape,double rate){
  return R::dnbinom(n,shape+total,(rate+size)/(rate+size+1),true);
}
// Collapsed DP/PY Gibbs sampling, scanning active clusters rather than all K slots.
// [[Rcpp::export]]
List collapsed_cpp(IntegerVector n,double shape,double rate,double theta,double discount,
 int sweeps,int burn,bool singletons=true,bool posterior_draws=false,int residual_terms=80){
  int K=n.size();if(K<1 || theta<=-discount || discount<0 || discount>=1 ||
    burn>=sweeps)stop("Invalid collapsed sampler arguments");
  vector<int> z(K),active,freeids;vector<double> sizes(K,0),totals(K,0),lw(K+1);
  if(singletons){for(int i=0;i<K;++i){z[i]=i;sizes[i]=1;totals[i]=n[i];active.push_back(i);}}
  else{std::fill(z.begin(),z.end(),0);sizes[0]=K;totals[0]=sum(n);active.push_back(0);
    for(int i=K-1;i>=1;--i)freeids.push_back(i);}
  double L0=std::pow(rate/(rate+1),shape);
  NumericVector rb(sweeps-burn),adraw(sweeps-burn),clusters(sweeps-burn);
  for(int it=0;it<sweeps;++it){
    if(it%32==0)checkUserInterrupt();
    for(int i=0;i<K;++i){
      int c=z[i];sizes[c]-=1;totals[c]-=n[i];
      if(sizes[c]==0){
        auto pos=std::find(active.begin(),active.end(),c);active.erase(pos);freeids.push_back(c);
      }
      double mx=R_NegInf;int k=active.size();
      for(int h=0;h<k;++h){
        int c=active[h];lw[h]=std::log(sizes[c]-discount)+predictive(n[i],sizes[c],totals[c],shape,rate);
        mx=std::max(mx,lw[h]);
      }
      lw[k]=std::log(theta+k*discount)+predictive(n[i],0,0,shape,rate);mx=std::max(mx,lw[k]);
      double total=0;for(int h=0;h<=k;++h){lw[h]=std::exp(lw[h]-mx);total+=lw[h];}
      double u=R::runif(0,total);int h=0;while(h<k && u>lw[h])u-=lw[h++];
      if(h==k){c=freeids.back();freeids.pop_back();active.push_back(c);}else c=active[h];
      z[i]=c;sizes[c]+=1;totals[c]+=n[i];
    }
    if(it>=burn){
      int out=it-burn,k=active.size();double val=(theta+k*discount)*L0;
      for(int c:active)val+=(sizes[c]-discount)*std::pow((rate+sizes[c])/(rate+sizes[c]+1),shape+totals[c]);
      rb[out]=val/(theta+K);clusters[out]=k;
      if(posterior_draws){
        double res=R::rgamma(theta+k*discount,1),den=res,num=0;
        for(int c:active){
          double mass=R::rgamma(sizes[c]-discount,1);
          den+=mass;num+=mass*std::exp(-positive_gamma(shape+totals[c],rate+sizes[c]));
        }
        double remaining=res;
        for(int h=1;h<residual_terms;++h){
          double v=R::rbeta(1-discount,theta+k*discount+h*discount);
          num+=remaining*v*std::exp(-positive_gamma(shape,rate));remaining*=1-v;
        }
        num+=remaining*std::exp(-positive_gamma(shape,rate));adraw[out]=num/den;
      }
    }
  }
  return List::create(_["rb"]=rb,_["a"]=adraw,_["clusters"]=clusters);
}

// [[Rcpp::export]]
NumericVector truncated_gamma_check_cpp(int n,double shape,double rate,double lower){
  NumericVector x(n);for(int i=0;i<n;++i)x[i]=positive_gamma(shape,rate,lower);return x;
}

// Independent prior importance sampling; no MCMC update is used here.
// [[Rcpp::export]]
List importance_cpp(IntegerVector j,NumericVector f,double shape,double rate,
 double alpha,int H,int draws,bool truncated=true,double lower=0,int cv=0,double cf=0){
  double K=sum(f)+cf,sw=0,sw2=0;
  vector<double> wx(3,0),w2x(3,0),w2x2(3,0),mu(j.size()),x(3);
  for(int it=0;it<draws;++it){
    if(it%16384==0)checkUserInterrupt();
    std::fill(mu.begin(),mu.end(),0);double rem=1,Z=0,T=0;
    for(int h=0;h<H;++h){
      double v=h==H-1?1:R::rbeta(1,alpha),w=rem*v;rem*=1-v;
      double rate0=positive_gamma(shape,rate,lower);
      Z+=w*(-std::expm1(-rate0));
      for(int z=0;z<j.size();++z)mu[z]+=w*R::dpois(j[z],rate0,false);
      if(cf>0)T+=w*R::ppois(cv-1,rate0,false,false);
    }
    double lw=loglik(mu,T,Z,f,cf,K,truncated),w=std::exp(lw);
    x[0]=1-Z;x[1]=x[0]>.5;x[2]=R::pnbinom(K,K,Z,true,false);
    sw+=w;sw2+=w*w;
    for(int z=0;z<3;++z){wx[z]+=w*x[z];w2x[z]+=w*w*x[z];w2x2[z]+=w*w*x[z]*x[z];}
  }
  NumericVector mean(3),se(3);
  for(int z=0;z<3;++z){
    mean[z]=wx[z]/sw;
    se[z]=std::sqrt(std::max(0.,w2x2[z]-2*mean[z]*w2x[z]+mean[z]*mean[z]*sw2))/sw;
  }
  return List::create(_["mean"]=mean,_["mcse"]=se,_["ess"]=sw*sw/sw2,_["draws"]=draws);
}
