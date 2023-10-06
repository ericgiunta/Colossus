using namespace std;
using namespace Rcpp;

using Eigen::Map;
using Eigen::MatrixXd;
using Eigen::SparseMatrix;
using Eigen::VectorXd;
using Rcpp::as;


bool Check_Risk( IntegerVector Term_n, StringVector tform, NumericVector a_n,NumericMatrix x_all,IntegerVector dfc,int fir,string modelform, bool verbose, bool debugging, IntegerVector KeepConstant, int term_tot, int nthreads, const double gmix_theta, const IntegerVector gmix_term);


List LogLik_Cox_PH_Omnibus( IntegerVector Term_n, StringVector tform, NumericMatrix a_ns,NumericMatrix x_all,IntegerVector dfc,int fir, int der_iden,string modelform, double lr, NumericVector maxiters, int guesses, int halfmax, double epsilon, double dbeta_cap, double abs_max,double dose_abs_max, double deriv_epsilon, NumericMatrix df_groups, NumericVector tu, int double_step ,bool change_all, bool verbose, bool debugging, IntegerVector KeepConstant, int term_tot, string ties_method, int nthreads, NumericVector& STRATA_vals, const VectorXd cens_weight, const double cens_thres, bool strata_bool, bool basic_bool, bool null_bool, bool CR_bool, bool single_bool, bool constraint_bool, const double gmix_theta, const IntegerVector gmix_term, const MatrixXd Lin_Sys, const VectorXd Lin_Res);

List LogLik_Pois_Omnibus(MatrixXd PyrC, IntegerVector Term_n, StringVector tform, NumericMatrix a_ns,NumericMatrix x_all,IntegerVector dfc,int fir, int der_iden,string modelform, double lr, NumericVector maxiters, int guesses, int halfmax, double epsilon, double dbeta_cap, double abs_max,double dose_abs_max, double deriv_epsilon, int double_step ,bool change_all, bool verbose, bool debugging, IntegerVector KeepConstant, int term_tot, int nthreads, const MatrixXd& dfs, bool strata_bool, bool single_bool, bool constraint_bool, const double gmix_theta, const IntegerVector gmix_term, const MatrixXd Lin_Sys, const VectorXd Lin_Res);
