#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::plugins("cpp11")]]


// include our bark beetle code
#include "bbgenerations.h"



// [[Rcpp::export]]
void test() {
  BBGenerations bb;
  bb.setLatitude(45);
}

void dolog(const std::string &s) { 
  Rcpp::Rcout << "bb: " << s << "\n"; 
}

// a simple wrapper class for bb generations:
class BBWrapper {
public:
  
  void setup(double latitude_deg, double lai, bool do_log) {
    bb.setLatitude(latitude_deg);
    bb.setLAI(lai);
    Rcpp::Rcout << "doy with < 14.5 hrs: " << bb.mDay_14_5hrs;
    if (do_log)
      bb.logfunc = &dolog;
    else
      bb.logfunc = nullptr;
  }
  void setClimate(const NumericMatrix &mat) {
     //Rcpp::Rcout << "matrix rows/cols" << mat.nrow() << "/"<< mat.ncol() << ". climate #rows:" << bb.mClimate.size() << "\n";
     bb.mClimate.resize(366);
    strong_frost_days = 0;
     for (int i=0;i<mat.nrow();++i) {
       BBGenerations::SDay &day = bb.mClimate[i];
       day.tmin = mat(i, 0);
       day.tmax = mat(i, 1);
       day.radiation = mat(i, 2);
       if (day.tmin < -15)
         ++strong_frost_days;
     }
  }
  NumericVector run(const NumericMatrix &mat) {
    setClimate(mat);
    double result = bb.calculateGenerations();
    NumericVector v = NumericVector::create(result, strong_frost_days);
    return v;
  }
  
  NumericVector barkTemps() {
    NumericVector r = NumericVector(bb.mEffectiveBarkTemp, bb.mEffectiveBarkTemp+366);
    return r;
  }
  
  void test() {
    for (const auto &d : bb.mClimate) {
      Rcpp::Rcout << d.tmin << ", " << d.tmax << ", " << d.radiation << "\n";
    }
  }
private:
 
  BBGenerations bb; // this is the iLand phenips code with minor mods
  int strong_frost_days;
};


RCPP_MODULE(econ_evaluator_module) {
  class_<BBWrapper>( "BarkBeetleGen" )
  // .field( "min", &Uniform::min )
     .default_constructor()
  
  .method( "setup", &BBWrapper::setup )
  .method( "run", &BBWrapper::run)
  .method( "barkTemps", &BBWrapper::barkTemps)
  .method( "test", &BBWrapper::test);
}



