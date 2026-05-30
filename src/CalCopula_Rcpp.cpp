#include <Rcpp.h>
using namespace Rcpp;

inline int mode_code(const std::string& mode){
  if(mode=="1")   return 1;
  if(mode=="2")   return 2;
  if(mode=="11")  return 11;
  if(mode=="12")  return 12;
  if(mode=="22")  return 22;
  if(mode=="121") return 121;
  if(mode=="122") return 122;
  return -1;
}

inline int out_size(int n1, int n2){
  return std::max(n1, n2);
}

inline bool recycle_compatible(int n1, int n2){
  int n = std::max(n1, n2);
  return (n % n1 == 0) && (n % n2 == 0);
}

inline double get_recycled(const NumericVector& x, int i){
  int n = x.size();
  return x[i % n];
}

// [[Rcpp::export]]
NumericVector archmCopulaLink_cpp(std::string fam, double param, NumericVector p){
  int n = p.size();
  NumericVector out(n);

  if (fam == "amh"){
    for(int i=0;i<n;i++){
      double pi = p[i];
      out[i] = std::log(param + (1.0 - param)/pi);
    }
    return out;
  }
  if (fam == "clayton"){
    for(int i=0;i<n;i++){
      double pi = p[i];
      out[i] = std::pow(pi, -param) - 1.0;
    }
    return out;
  }
  if (fam == "frank"){
    double den = std::exp(-param) - 1.0;
    for(int i=0;i<n;i++){
      double pi = p[i];
      out[i] = -std::log((std::exp(-param*pi) - 1.0)/den);
    }
    return out;
  }
  if (fam == "gumbel"){
    for(int i=0;i<n;i++){
      double pi = p[i];
      out[i] = std::pow(-std::log(pi), param);
    }
    return out;
  }
  if (fam == "joe"){
    for(int i=0;i<n;i++){
      double pi = p[i];
      out[i] = -std::log(1.0 - std::pow(1.0 - pi, param));
    }
    return out;
  }

  stop("Unsupported copulafam in archmCopulaLink_cpp.");
  return out;
}

// [[Rcpp::export]]
NumericVector archmCopulaLink_dev_cpp(std::string fam, double param, NumericVector js){
  int n = js.size();
  NumericVector out(n);

  if (fam == "amh"){
    for(int i=0;i<n;i++){
      double s = js[i];
      out[i] = (param - 1.0)/(param*s*s + (1.0-param)*s);
    }
    return out;
  }
  if (fam == "clayton"){
    for(int i=0;i<n;i++){
      double s = js[i];
      out[i] = -param * std::pow(s, -param-1.0);
    }
    return out;
  }
  if (fam == "frank"){
    for(int i=0;i<n;i++){
      double s = js[i];
      double e = std::exp(-param*s);
      out[i] = param * e/(e - 1.0);
    }
    return out;
  }
  if (fam == "gumbel"){
    for(int i=0;i<n;i++){
      double s = js[i];
      out[i] = -param * std::pow(-std::log(s), param-1.0)/s;
    }
    return out;
  }
  if (fam == "joe"){
    for(int i=0;i<n;i++){
      double s = js[i];
      out[i] = -param * std::pow(1.0 - s, param-1.0)/(1.0 - std::pow(1.0 - s, param));
    }
    return out;
  }

  stop("Unsupported copulafam in archmCopulaLink_dev_cpp.");
  return out;
}

// [[Rcpp::export]]
NumericVector archmCopulaLink_inv_cpp(std::string fam, double param, NumericVector y){
  int n = y.size();
  NumericVector out(n);

  if (fam == "amh"){
    for(int i = 0; i < n; i++){
      double yi = y[i];
      out[i] = (param - 1.0) / (param - std::exp(yi));
    }
    return out;
  }

  if (fam == "clayton"){
    for(int i = 0; i < n; i++){
      double yi = y[i];
      out[i] = std::pow(yi + 1.0, -1.0 / param);
    }
    return out;
  }

  if (fam == "frank"){
    for(int i = 0; i < n; i++){
      double yi = y[i];
      out[i] = -1.0 / param * std::log(1.0 + (std::exp(-param) - 1.0) * std::exp(-yi));
    }
    return out;
  }

  if (fam == "gumbel"){
    for(int i = 0; i < n; i++){
      double yi = y[i];
      out[i] = std::exp(-std::pow(yi, 1.0 / param));
    }
    return out;
  }

  if (fam == "joe"){
    for(int i = 0; i < n; i++){
      double yi = y[i];
      out[i] = 1.0 - std::pow(1.0 - std::exp(-yi), 1.0 / param);
    }
    return out;
  }

  stop("Unsupported copulafam in archmCopulaLink_inv_cpp.");
  return out;
}

// fallback: mimic the R code using acopula::nderive and Copulafn
double dev_copula_fallback_one(std::string fam, double param, double u, double v, std::string mode){
  Environment global_env = Environment::global_env();

  if (!global_env.exists("dev_Copula_fallback_R")) {
    stop("dev_Copula_fallback_R not found in global environment.");
  }

  Function fallback = global_env["dev_Copula_fallback_R"];

  SEXP res = fallback(_["copulafam"] = fam,
                      _["param"] = param,
                      _["p1"] = u,
                      _["p2"] = v,
                      _["mode"] = mode);

  return as<double>(res);
}

// [[Rcpp::export]]
NumericVector dev_Copula_cpp(std::string fam, double param,
                             NumericVector p1, NumericVector p2,
                             std::string mode = "1"){

  int m = mode_code(mode);
  if(m < 0) stop("Unsupported mode in dev_Copula_cpp.");

  int n1 = p1.size(), n2 = p2.size();
  if (!recycle_compatible(n1, n2)) {
    warning("longer object length is not a multiple of shorter object length");
  }

  int n = out_size(n1, n2);
  NumericVector out(n, NA_REAL);

  // ---------------- amh ----------------
  if (fam == "amh") {
    if (!(param >= -1.0 && param <= 1.0)) {
      std::fill(out.begin(), out.end(), NA_REAL);
      return out;
    }

    for (int i = 0; i < n; i++) {
      double u = get_recycled(p1, i);
      double v = get_recycled(p2, i);
      double a = 1.0 - param * (1.0 - u) * (1.0 - v);

      if (m == 1) {
        out[i] = (v * a - u * v * param * (1.0 - v)) / (a * a);
      } else if (m == 2) {
        out[i] = (u * a - u * v * param * (1.0 - u)) / (a * a);
      } else if (m == 12) {
        out[i] = 1.0 - param
          + 2.0 * param * v / (a * a)
          - 2.0 * param * v * (1.0 - param + param * v) * (1.0 - u) / (a * a * a);
      } else {
        out[i] = dev_copula_fallback_one(fam, param, u, v, mode);
      }
    }
    return out;
  }

  // ---------------- clayton ----------------
  if (fam == "clayton") {
    if (!(param >= -1.0)) {
      std::fill(out.begin(), out.end(), NA_REAL);
      return out;
    }

    for (int i = 0; i < n; i++) {
      double u = get_recycled(p1, i);
      double v = get_recycled(p2, i);

      if (m == 1) {
        double val = std::pow(std::pow(u, -param) + std::pow(v, -param) - 1.0,
                              -1.0/param - 1.0) * std::pow(u, -param - 1.0);
        out[i] = R_IsNaN(val) ? 0.0 : val;
      } else if (m == 2) {
        double val = std::pow(std::pow(u, -param) + std::pow(v, -param) - 1.0,
                              -1.0/param - 1.0) * std::pow(v, -param - 1.0);
        out[i] = R_IsNaN(val) ? 0.0 : val;
      } else if (m == 12) {
        out[i] = (1.0 + param) *
          std::pow(std::pow(u, -param) + std::pow(v, -param) - 1.0,
                   -1.0/param - 2.0) *
          std::pow(u, -param - 1.0) *
          std::pow(v, -param - 1.0);
      } else {
        out[i] = dev_copula_fallback_one(fam, param, u, v, mode);
      }
    }
    return out;
  }

  // ---------------- frank ----------------
  if (fam == "frank") {
    if (param == 0.0) {
      std::fill(out.begin(), out.end(), NA_REAL);
      return out;
    }

    for (int i = 0; i < n; i++) {
      double u = get_recycled(p1, i);
      double v = get_recycled(p2, i);

      double out0 = (std::exp(-param*u) - 1.0) * (std::exp(-param*v) - 1.0);
      double out1 = out0 + (std::exp(-param) - 1.0);

      if (m == 1) {
        out[i] = (std::exp(-param*v) - 1.0) * std::exp(-param*u) / out1;
      } else if (m == 2) {
        out[i] = (std::exp(-param*u) - 1.0) * std::exp(-param*v) / out1;
      } else if (m == 11) {
        out[i] = param * (std::exp(-param*v) - 1.0) * std::exp(-param*u) *
          (std::exp(-param*v) - std::exp(-param)) / (out1 * out1);
      } else if (m == 12) {
        out[i] = (1.0 - std::exp(-param)) * param *
          std::exp(-param*u) * std::exp(-param*v) / (out1 * out1);
      } else if (m == 22) {
        out[i] = param * (std::exp(-param*u) - 1.0) * std::exp(-param*v) *
          (std::exp(-param*u) - std::exp(-param)) / (out1 * out1);
      } else if (m == 121) {
        double out2 = (1.0 - std::exp(-param)) * param * std::exp(-param*u) * std::exp(-param*v);
        out[i] = ( out2 * (-param) * (out1 * out1)
          - 2.0 * out2 * out1 * (std::exp(-param*v) - 1.0) * std::exp(-param*u) * (-param) ) /
          std::pow(out1, 4.0);
      } else if (m == 122) {
        double out2 = (1.0 - std::exp(-param)) * param * std::exp(-param*u) * std::exp(-param*v);
        out[i] = ( out2 * (-param) * (out1 * out1)
          - 2.0 * out2 * out1 * (std::exp(-param*u) - 1.0) * std::exp(-param*v) * (-param) ) /
          std::pow(out1, 4.0);
      } else {
        out[i] = dev_copula_fallback_one(fam, param, u, v, mode);
      }
    }
    return out;
  }

  // ---------------- gumbel ----------------
  if (fam == "gumbel") {
    if (!(param >= 1.0)) {
      std::fill(out.begin(), out.end(), NA_REAL);
      return out;
    }

    for (int i = 0; i < n; i++) {
      double u = get_recycled(p1, i);
      double v = get_recycled(p2, i);

      double lu = -std::log(u);
      double lv = -std::log(v);
      double A  = std::pow(lu, param) + std::pow(lv, param);
      double C  = std::exp(-std::pow(A, 1.0/param));

      if (m == 1) {
        out[i] = C * std::pow(A, 1.0/param - 1.0) * std::pow(lu, param-1.0) / u;
      } else if (m == 2) {
        out[i] = C * std::pow(A, 1.0/param - 1.0) * std::pow(lv, param-1.0) / v;
      } else if (m == 12) {
        out[i] = C * (
          std::pow(A, 2.0/param - 2.0) * std::pow(lu, param-1.0)/u * std::pow(lv, param-1.0)/v +
          std::pow(A, 1.0/param - 2.0) * std::pow(lu, param-1.0)/u * std::pow(lv, param-1.0)/v * (param-1.0)
        );
      } else {
        out[i] = dev_copula_fallback_one(fam, param, u, v, mode);
      }
    }
    return out;
  }

  // ---------------- joe ----------------
  if (fam == "joe") {
    if (!(param >= 1.0)) {
      std::fill(out.begin(), out.end(), NA_REAL);
      return out;
    }

    for (int i = 0; i < n; i++) {
      double u = get_recycled(p1, i);
      double v = get_recycled(p2, i);

      double a = std::pow(1.0-u, param) + std::pow(1.0-v, param)
        - std::pow(1.0-u, param) * std::pow(1.0-v, param);

      if (m == 1) {
        out[i] = std::pow(a, 1.0/param - 1.0) *
          ( std::pow(1.0-u, param-1.0) - std::pow(1.0-v, param) * std::pow(1.0-u, param-1.0) );
      } else if (m == 2) {
        out[i] = std::pow(a, 1.0/param - 1.0) *
          ( std::pow(1.0-v, param-1.0) - std::pow(1.0-u, param) * std::pow(1.0-v, param-1.0) );
      } else if (m == 12) {
        out[i] =
          (param - 1.0) * std::pow(a, 1.0/param - 2.0) *
          std::pow(1.0-u, param-1.0) * std::pow(1.0-v, param-1.0) *
          (1.0 - std::pow(1.0-u, param)) * (1.0 - std::pow(1.0-v, param))
          +
          std::pow(a, 1.0/param - 1.0) * param *
          std::pow(1.0-u, param-1.0) * std::pow(1.0-v, param-1.0);
      } else {
        out[i] = dev_copula_fallback_one(fam, param, u, v, mode);
      }
    }
    return out;
  }

  // fully mimic R fallback behavior
  for (int i = 0; i < n; i++) {
    double u = get_recycled(p1, i);
    double v = get_recycled(p2, i);
    out[i] = dev_copula_fallback_one(fam, param, u, v, mode);
  }
  return out;
}