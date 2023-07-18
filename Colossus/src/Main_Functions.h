using namespace std;
using namespace Rcpp;

using Eigen::Map;
using Eigen::MatrixXd;
using Eigen::SparseMatrix;
using Eigen::VectorXd;
using Rcpp::as;


List LogLik_Cox_PH( IntegerVector Term_n, StringVector tform, NumericVector a_n,NumericMatrix x_all,IntegerVector dfc,int fir, int der_iden,string modelform, double lr, int maxiter, int halfmax, double epsilon, double dbeta_cap, double abs_max,double dose_abs_max, double deriv_epsilon, NumericMatrix df_groups, NumericVector tu, int double_step ,bool change_all, bool verbose, bool debugging, IntegerVector KeepConstant, int term_tot, string ties_method, int nthreads);

List Cox_PH_PLOT_SURV(IntegerVector Term_n, StringVector tform, NumericVector a_n, NumericVector a_er,NumericMatrix x_all,IntegerVector dfc,int fir, int der_iden,string modelform, double abs_max,double dose_abs_max, NumericMatrix df_groups, NumericVector tu , bool verbose, bool debugging, IntegerVector KeepConstant, int term_tot, int nthreads);

List Cox_PH_PLOT_RISK(IntegerVector Term_n, StringVector tform, NumericVector a_n,NumericMatrix x_all,IntegerVector dfc,int fir, int der_iden,string modelform, double abs_max,double dose_abs_max, NumericMatrix df_groups, NumericVector tu , bool verbose, bool debugging, IntegerVector KeepConstant, int term_tot, int uniq_v, int nthreads);

List Schoenfeld_Cox_PH( IntegerVector Term_n, StringVector tform, NumericVector a_n,NumericMatrix x_all,IntegerVector dfc,int fir, int der_iden,string modelform, double abs_max,double dose_abs_max, NumericMatrix df_groups, NumericVector tu , bool verbose, bool debugging, IntegerVector KeepConstant, int term_tot,string ties_method, int nthreads);

List LogLik_Poisson( MatrixXd PyrC, IntegerVector Term_n, StringVector tform, NumericVector a_n,NumericMatrix x_all,IntegerVector dfc,int fir, int der_iden,string modelform, double lr, int maxiter, int halfmax, double epsilon, double dbeta_cap, double abs_max,double dose_abs_max, double deriv_epsilon, int double_step,bool change_all, bool verbose, bool debugging, IntegerVector KeepConstant, int term_tot, int nthreads);

void Stress_Run( IntegerVector Term_n, StringVector tform, NumericVector a_n,NumericMatrix x_all,IntegerVector dfc,int fir, int der_iden,string modelform, double lr, int maxiter, int halfmax, double epsilon, double dbeta_cap, double abs_max,double dose_abs_max, double deriv_epsilon, NumericMatrix df_groups, NumericVector tu , int double_step,bool change_all, bool verbose, bool debugging, IntegerVector KeepConstant, int term_tot, StringVector debug_checks, string ties_method, int nthreads);

NumericVector RISK_SUBSET(IntegerVector Term_n, StringVector tform, NumericVector a_n,NumericMatrix x_all,IntegerVector dfc,int fir,string modelform, bool verbose, bool debugging, int term_tot, int nthreads);

bool Check_Risk( IntegerVector Term_n, StringVector tform, NumericVector a_n,NumericMatrix x_all,IntegerVector dfc,int fir,string modelform, bool verbose, bool debugging, IntegerVector KeepConstant, int term_tot, int nthreads);


List LogLik_Cox_PH_Omnibus( IntegerVector Term_n, StringVector tform, NumericMatrix a_ns,NumericMatrix x_all,IntegerVector dfc,int fir, int der_iden,string modelform, double lr, NumericVector maxiters, int guesses, int halfmax, double epsilon, double dbeta_cap, double abs_max,double dose_abs_max, double deriv_epsilon, NumericMatrix df_groups, NumericVector tu, int double_step ,bool change_all, bool verbose, bool debugging, IntegerVector KeepConstant, int term_tot, string ties_method, int nthreads, NumericVector& STRATA_vals, const VectorXd cens_weight, const double cens_thres, bool strata_bool, bool basic_bool, bool null_bool, bool CR_bool, bool single_bool);

List LogLik_Pois_Omnibus(MatrixXd PyrC, IntegerVector Term_n, StringVector tform, NumericMatrix a_ns,NumericMatrix x_all,IntegerVector dfc,int fir, int der_iden,string modelform, double lr, NumericVector maxiters, int guesses, int halfmax, double epsilon, double dbeta_cap, double abs_max,double dose_abs_max, double deriv_epsilon, int double_step ,bool change_all, bool verbose, bool debugging, IntegerVector KeepConstant, int term_tot, int nthreads, const MatrixXd& dfs, bool strata_bool, bool single_bool);
