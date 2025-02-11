using namespace std;
using namespace Rcpp;

using Eigen::Map;
using Eigen::MatrixXd;
using Eigen::SparseMatrix;
using Eigen::VectorXd;
using Rcpp::as;

List PLOT_SURV(int reqrdnum, MatrixXd& R, MatrixXd& Rd, NumericVector& a_er, NumericMatrix& df_groups, NumericVector& tu, int verbose, int nthreads);

List PLOT_SURV_Strata(int reqrdnum, MatrixXd& R, MatrixXd& Rd, NumericVector& a_er, NumericMatrix& df_groups, NumericVector& tu, NumericVector& Strata_vals, int verbose, int nthreads);

List Schoenfeld_Calc(int ntime, int totalnum, const  VectorXd& beta_0, const  MatrixXd& df0, const MatrixXd& R, MatrixXd& Lldd_inv, const IntegerMatrix& RiskFail, const vector<vector<int> >& RiskPairs, IntegerVector& dfc, int verbose, IntegerVector KeepConstant, int nthreads);

List Plot_Omnibus(IntegerVector term_n, StringVector tform, NumericVector a_n, NumericMatrix& x_all, IntegerVector dfc, int fir, int der_iden, string modelform, double abs_max, double dose_abs_max, NumericMatrix& df_groups, NumericVector& tu, int verbose, IntegerVector KeepConstant, int term_tot, string ties_method, int nthreads, NumericVector& Strata_vals, const VectorXd& cens_weight, int uniq_v, List model_bool, bool Surv_bool, bool Risk_bool, bool Schoenfeld_bool, bool Risk_Sub_bool, const double gmix_theta, const IntegerVector& gmix_term);

List Assign_Events_Pois(IntegerVector term_n, StringVector tform, NumericVector a_n, NumericMatrix& x_all, IntegerVector dfc, const MatrixXd& PyrC, const MatrixXd& dfs, int fir, string modelform, int verbose, IntegerVector KeepConstant, int term_tot, int nthreads, const double gmix_theta, const IntegerVector gmix_term, List model_bool);

List Poisson_Residuals(const MatrixXd& PyrC, IntegerVector term_n, StringVector tform, NumericVector a_n, NumericMatrix& x_all, IntegerVector dfc, int fir, int der_iden, string modelform, double abs_max, double dose_abs_max, int verbose, IntegerVector KeepConstant, int term_tot, int nthreads, const MatrixXd& dfs, List model_bool, const double gmix_theta, const IntegerVector gmix_term, bool Pearson_bool, bool Deviance_bool);
