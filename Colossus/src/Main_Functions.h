using namespace std;
using namespace Rcpp;

using Eigen::Map;
using Eigen::MatrixXd;
using Eigen::SparseMatrix;
using Eigen::VectorXd;
using Rcpp::as;


List LogLik_Cox_PH( IntegerVector Term_n, StringVector tform, NumericVector a_n,NumericMatrix x_all,IntegerVector dfc,int fir, int der_iden,string modelform, double lr, int maxiter, int halfmax, double epsilon, double dbeta_cap, double abs_max,double dose_abs_max, double deriv_epsilon, NumericMatrix df_groups, NumericVector tu ,bool change_all, bool verbose, bool debugging, IntegerVector KeepConstant, int term_tot, string ties_method);

List Cox_PH_PLOT_SURV(IntegerVector Term_n, StringVector tform, NumericVector a_n,NumericMatrix x_all,IntegerVector dfc,int fir, int der_iden,string modelform, double abs_max,double dose_abs_max, NumericMatrix df_groups, NumericVector tu , bool verbose, bool debugging, IntegerVector KeepConstant, int term_tot);

List Cox_PH_PLOT_RISK(IntegerVector Term_n, StringVector tform, NumericVector a_n,NumericMatrix x_all,IntegerVector dfc,int fir, int der_iden,string modelform, double abs_max,double dose_abs_max, NumericMatrix df_groups, NumericVector tu , bool verbose, bool debugging, IntegerVector KeepConstant, int term_tot, int uniq_v);

NumericMatrix Schoenfeld_Cox_PH( IntegerVector Term_n, StringVector tform, NumericVector a_n,NumericMatrix x_all,IntegerVector dfc,int fir, int der_iden,string modelform, double abs_max,double dose_abs_max, NumericMatrix df_groups, NumericVector tu , bool verbose, bool debugging, IntegerVector KeepConstant, int term_tot,string ties_method);

List LogLik_Poisson( MatrixXd PyrC, IntegerVector Term_n, StringVector tform, NumericVector a_n,NumericMatrix x_all,IntegerVector dfc,int fir, int der_iden,string modelform, double lr, int maxiter, int halfmax, double epsilon, double dbeta_cap, double abs_max,double dose_abs_max, double deriv_epsilon,bool change_all, bool verbose, bool debugging, IntegerVector KeepConstant, int term_tot);

void Stress_Run( IntegerVector Term_n, StringVector tform, NumericVector a_n,NumericMatrix x_all,IntegerVector dfc,int fir, int der_iden,string modelform, double lr, int maxiter, int halfmax, double epsilon, double dbeta_cap, double abs_max,double dose_abs_max, double deriv_epsilon, NumericMatrix df_groups, NumericVector tu ,bool change_all, bool verbose, bool debugging, IntegerVector KeepConstant, int term_tot, StringVector debug_checks, string ties_method);

List LogLik_Cox_PH_null( int ntime, NumericMatrix df_groups, NumericVector tu, bool verbose, string ties_method);

NumericVector RISK_SUBSET(IntegerVector Term_n, StringVector tform, NumericVector a_n,NumericMatrix x_all,IntegerVector dfc,int fir,string modelform, bool verbose, bool debugging, int term_tot);
