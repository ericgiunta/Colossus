using namespace std;
using namespace Rcpp;

using Eigen::Map;
using Eigen::MatrixXd;
using Eigen::SparseMatrix;
using Eigen::VectorXd;
using Rcpp::as;


List LogLik_Cox_PH_Omnibus_Log_Bound(IntegerVector term_n, StringVector tform, NumericVector a_n, NumericMatrix& x_all, IntegerVector dfc, int fir, string modelform, double lr, NumericVector maxiters, int guesses, int halfmax, double epsilon, double abs_max, double dose_abs_max, double deriv_epsilon, NumericMatrix df_groups, NumericVector tu, int verbose, IntegerVector KeepConstant, int term_tot, string ties_method, int nthreads, NumericVector& Strata_vals, const VectorXd& cens_weight, List model_bool, const double gmix_theta, const IntegerVector gmix_term, const MatrixXd Lin_Sys, const VectorXd Lin_Res, double qchi, int para_number, int maxstep, double mult);

List LogLik_Cox_PH_Omnibus_Log_Bound_Search(IntegerVector term_n, StringVector tform, NumericVector a_n, NumericMatrix& x_all, IntegerVector dfc, int fir, string modelform, double lr, NumericVector maxiters, int guesses, int halfmax, double epsilon, double abs_max, double dose_abs_max, double deriv_epsilon, NumericMatrix df_groups, NumericVector tu, int verbose, IntegerVector KeepConstant, int term_tot, string ties_method, int nthreads, NumericVector& Strata_vals, const VectorXd& cens_weight, List model_bool, const double gmix_theta, const IntegerVector gmix_term, const MatrixXd Lin_Sys, const VectorXd Lin_Res, double qchi, int para_number, int maxstep, double mult);

List LogLik_Poisson_Omnibus_Log_Bound(const MatrixXd& PyrC, const MatrixXd& dfs, IntegerVector term_n, StringVector tform, NumericVector a_n, NumericMatrix& x_all, IntegerVector dfc, int fir, string modelform, double lr, NumericVector maxiters, int guesses, int halfmax, double epsilon, double abs_max, double dose_abs_max, double deriv_epsilon, int verbose, IntegerVector KeepConstant, int term_tot, int nthreads, List model_bool, const double gmix_theta, const IntegerVector gmix_term, const MatrixXd Lin_Sys, const VectorXd Lin_Res, double qchi, int para_number, int maxstep, double mult);

List LogLik_Poisson_Omnibus_Log_Bound_Search(const MatrixXd& PyrC, const MatrixXd& dfs, IntegerVector term_n, StringVector tform, NumericVector a_n, NumericMatrix& x_all, IntegerVector dfc, int fir, string modelform, double lr, NumericVector maxiters, int guesses, int halfmax, double epsilon, double abs_max, double dose_abs_max, double deriv_epsilon, int verbose, IntegerVector KeepConstant, int term_tot, int nthreads, List model_bool, const double gmix_theta, const IntegerVector gmix_term, const MatrixXd Lin_Sys, const VectorXd Lin_Res, double qchi, int para_number, int maxstep, double mult);

List LogLik_Cox_PH_Omnibus_Log_Bound_CurveSearch(IntegerVector term_n, StringVector tform, NumericVector a_n, NumericMatrix& x_all, IntegerVector dfc, int fir, string modelform, double lr, List optim_para, int maxiter, int double_step, int halfmax, double epsilon, double abs_max, double dose_abs_max, double deriv_epsilon, NumericMatrix df_groups, NumericVector tu, int verbose, IntegerVector KeepConstant, int term_tot, string ties_method, int nthreads, NumericVector& Strata_vals, const VectorXd& cens_weight, List model_bool, const double gmix_theta, const IntegerVector gmix_term, const MatrixXd Lin_Sys, const VectorXd Lin_Res, double qchi, int para_number, int maxstep, double step_size);

List LogLik_Poisson_Omnibus_Log_Bound_CurveSearch(const MatrixXd& PyrC, const MatrixXd& dfs, IntegerVector term_n, StringVector tform, NumericVector a_n, NumericMatrix& x_all, IntegerVector dfc, int fir, string modelform, double lr, List optim_para, int maxiter, int doublestep, int halfmax, double epsilon, double abs_max, double dose_abs_max, double deriv_epsilon, int verbose, IntegerVector KeepConstant, int term_tot, int nthreads, List model_bool, const double gmix_theta, const IntegerVector gmix_term, const MatrixXd Lin_Sys, const VectorXd Lin_Res, double qchi, int para_number, int maxstep, double step_size);

