using namespace Rcpp;

using std::string;
using std::vector;

using Eigen::Map;
using Eigen::MatrixXd;
using Eigen::SparseMatrix;
using Eigen::VectorXd;
using Rcpp::as;

List LogLik_Cox_PH_Multidose_Omnibus_Serial(IntegerVector term_n, StringVector tform, NumericVector a_n, NumericMatrix& x_all, NumericMatrix& dose_all, IntegerMatrix dose_cols, IntegerVector dose_index, IntegerVector dfc, int fir, string modelform, double lr, List optim_para, int maxiter, int halfmax, double epsilon, double step_max, double thres_step_max, double deriv_epsilon, NumericMatrix df_groups, NumericVector tu, int verbose, IntegerVector KeepConstant, int term_tot, string ties_method, int nthreads, NumericVector& Strata_vals, const VectorXd& cens_weight, List model_bool, const double gmix_theta, const IntegerVector gmix_term, const MatrixXd Lin_Sys, const VectorXd Lin_Res);

List LogLik_Cox_PH_Multidose_Omnibus_Integrated(IntegerVector term_n, StringVector tform, NumericVector a_n, NumericMatrix& x_all, NumericMatrix& dose_all, IntegerMatrix dose_cols, IntegerVector dose_index, IntegerVector dfc, int fir, string modelform, double lr, List optim_para, int maxiter, int halfmax, double epsilon, double step_max, double thres_step_max, double deriv_epsilon, NumericMatrix df_groups, NumericVector tu, int verbose, IntegerVector KeepConstant, int term_tot, string ties_method, int nthreads, NumericVector& Strata_vals, const VectorXd& cens_weight, List model_bool, const double gmix_theta, const IntegerVector gmix_term, const MatrixXd Lin_Sys, const VectorXd Lin_Res);

List LogLik_Pois_PH_Multidose_Omnibus_Serial(const MatrixXd& PyrC, IntegerVector term_n, StringVector tform, NumericVector a_n, NumericMatrix& x_all, NumericMatrix& dose_all, IntegerMatrix dose_cols, IntegerVector dose_index, IntegerVector dfc, int fir, string modelform, double lr, List optim_para, int maxiter, int halfmax, double epsilon, double step_max, double thres_step_max, double deriv_epsilon, const MatrixXd& dfs, int verbose, IntegerVector KeepConstant, int term_tot, int nthreads, List model_bool, const double gmix_theta, const IntegerVector gmix_term, const MatrixXd Lin_Sys, const VectorXd Lin_Res);

List LogLik_Pois_PH_Multidose_Omnibus_Integrated(const MatrixXd& PyrC, IntegerVector term_n, StringVector tform, NumericVector a_n, NumericMatrix& x_all, NumericMatrix& dose_all, IntegerMatrix dose_cols, IntegerVector dose_index, IntegerVector dfc, int fir, string modelform, double lr, List optim_para, int maxiter, int halfmax, double epsilon, double step_max, double thres_step_max, double deriv_epsilon, const MatrixXd& dfs, int verbose, IntegerVector KeepConstant, int term_tot, int nthreads, List model_bool, const double gmix_theta, const IntegerVector gmix_term, const MatrixXd Lin_Sys, const VectorXd Lin_Res);
