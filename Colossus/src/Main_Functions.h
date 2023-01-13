using namespace std;
using namespace Rcpp;

using Eigen::Map;
using Eigen::MatrixXd;
using Eigen::SparseMatrix;
using Eigen::VectorXd;
using Rcpp::as;


List LogLik_Cox_PH( IntegerVector Term_n, StringVector tform, NumericVector a_n,NumericMatrix x_all,IntegerVector dfc,int fir, int der_iden,string modelform, double lr, int maxiter, int halfmax, double epsilon, double dbeta_cap, double abs_max,double dose_abs_max, double deriv_epsilon, NumericMatrix df_groups, NumericVector tu, int double_step ,bool change_all, bool verbose, bool debugging, IntegerVector KeepConstant, int term_tot, string ties_method, int nthreads);

List LogLik_Cox_PH_basic( NumericVector a_n,NumericMatrix x_all,IntegerVector dfc, int der_iden, double lr, int maxiter, int halfmax, double epsilon, double dbeta_cap, double abs_max,double dose_abs_max, double deriv_epsilon, NumericMatrix df_groups, NumericVector tu, int double_step ,bool change_all, bool verbose, bool debugging, IntegerVector KeepConstant, string ties_method, int nthreads);

List LogLik_Cox_PH_STRATA( IntegerVector Term_n, StringVector tform, NumericVector a_n,NumericMatrix x_all,IntegerVector dfc,int fir, int der_iden,string modelform, double lr, int maxiter, int halfmax, double epsilon, double dbeta_cap, double abs_max,double dose_abs_max, double deriv_epsilon, NumericMatrix df_groups, NumericVector tu, int double_step ,bool change_all, bool verbose, bool debugging, IntegerVector KeepConstant, int term_tot, string ties_method, IntegerVector& STRATA_vals, int nthreads);

List Cox_PH_PLOT_SURV(IntegerVector Term_n, StringVector tform, NumericVector a_n, NumericVector a_er,NumericMatrix x_all,IntegerVector dfc,int fir, int der_iden,string modelform, double abs_max,double dose_abs_max, NumericMatrix df_groups, NumericVector tu , bool verbose, bool debugging, IntegerVector KeepConstant, int term_tot, int nthreads);

List Cox_PH_PLOT_RISK(IntegerVector Term_n, StringVector tform, NumericVector a_n,NumericMatrix x_all,IntegerVector dfc,int fir, int der_iden,string modelform, double abs_max,double dose_abs_max, NumericMatrix df_groups, NumericVector tu , bool verbose, bool debugging, IntegerVector KeepConstant, int term_tot, int uniq_v, int nthreads);

NumericMatrix Schoenfeld_Cox_PH( IntegerVector Term_n, StringVector tform, NumericVector a_n,NumericMatrix x_all,IntegerVector dfc,int fir, int der_iden,string modelform, double abs_max,double dose_abs_max, NumericMatrix df_groups, NumericVector tu , bool verbose, bool debugging, IntegerVector KeepConstant, int term_tot,string ties_method, int nthreads);

List LogLik_Poisson( MatrixXd PyrC, IntegerVector Term_n, StringVector tform, NumericVector a_n,NumericMatrix x_all,IntegerVector dfc,int fir, int der_iden,string modelform, double lr, int maxiter, int halfmax, double epsilon, double dbeta_cap, double abs_max,double dose_abs_max, double deriv_epsilon, int double_step,bool change_all, bool verbose, bool debugging, IntegerVector KeepConstant, int term_tot, int nthreads);

List LogLik_Poisson_STRATA( MatrixXd PyrC, IntegerVector Term_n, StringVector tform, NumericVector a_n,NumericMatrix x_all,IntegerVector dfc,int fir, int der_iden,string modelform, double lr, int maxiter, int halfmax, double epsilon, double dbeta_cap, double abs_max,double dose_abs_max, double deriv_epsilon, int double_step,bool change_all, bool verbose, bool debugging, IntegerVector KeepConstant, int term_tot, IntegerVector& STRATA_vals, bool keep_strata, int nthreads);

void Stress_Run( IntegerVector Term_n, StringVector tform, NumericVector a_n,NumericMatrix x_all,IntegerVector dfc,int fir, int der_iden,string modelform, double lr, int maxiter, int halfmax, double epsilon, double dbeta_cap, double abs_max,double dose_abs_max, double deriv_epsilon, NumericMatrix df_groups, NumericVector tu , int double_step,bool change_all, bool verbose, bool debugging, IntegerVector KeepConstant, int term_tot, StringVector debug_checks, string ties_method, int nthreads);

List LogLik_Cox_PH_null( NumericMatrix df_groups, NumericVector tu, bool verbose, string ties_method, int nthreads);

NumericVector RISK_SUBSET(IntegerVector Term_n, StringVector tform, NumericVector a_n,NumericMatrix x_all,IntegerVector dfc,int fir,string modelform, bool verbose, bool debugging, int term_tot, int nthreads);

List LogLik_Cox_PH_Single( IntegerVector Term_n, StringVector tform, NumericVector a_n,NumericMatrix x_all,IntegerVector dfc,int fir,string modelform, NumericMatrix df_groups, NumericVector tu, bool verbose, bool debugging, int term_tot, string ties_method, int nthreads);

List LogLik_Poisson_Single( MatrixXd PyrC, IntegerVector Term_n, StringVector tform, NumericVector a_n,NumericMatrix x_all,IntegerVector dfc,int fir,string modelform, bool verbose, bool debugging, int term_tot, int nthreads);
