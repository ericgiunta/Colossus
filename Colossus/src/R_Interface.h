using namespace std;
using namespace Rcpp;

using Eigen::Map;
using Eigen::MatrixXd;
using Eigen::SparseMatrix;
using Eigen::VectorXd;
using Rcpp::as;

List cox_ph_Omnibus_transition(IntegerVector Term_n, StringVector tform, NumericMatrix a_ns,IntegerVector dfc,NumericMatrix x_all, int fir, int der_iden,string modelform, List Control, NumericMatrix df_groups, NumericVector tu, IntegerVector KeepConstant, int term_tot, NumericVector STRATA_vals, NumericVector cens_vec, List model_control);

List pois_Omnibus_transition(NumericMatrix dfe, IntegerVector Term_n, StringVector tform, NumericMatrix a_ns,IntegerVector dfc,NumericMatrix x_all, int fir, int der_iden,string modelform, List Control, IntegerVector KeepConstant, int term_tot, NumericMatrix df0, List model_control);

List cox_ph_transition(IntegerVector Term_n, StringVector tform, NumericVector a_n,IntegerVector dfc,NumericMatrix x_all, int fir, int der_iden,string modelform, List Control, NumericMatrix df_groups, NumericVector tu, IntegerVector KeepConstant, int term_tot);

List cox_ph_plot(IntegerVector Term_n, StringVector tform, NumericVector a_n, NumericVector a_er,IntegerVector dfc,NumericMatrix x_all, int fir, int der_iden,string modelform, List Control, NumericMatrix df_groups, NumericVector tu, IntegerVector KeepConstant, int term_tot, vector<string> Plot_Type ,int uniq_v);

List cox_ph_schoenfeld_transition(IntegerVector Term_n, StringVector tform, NumericVector a_n,IntegerVector dfc,NumericMatrix x_all, int fir, int der_iden,string modelform, List Control, NumericMatrix df_groups, NumericVector tu, IntegerVector KeepConstant, int term_tot);

List poisson_transition(NumericMatrix dfe, IntegerVector Term_n, StringVector tform, NumericVector a_n,IntegerVector dfc,NumericMatrix x_all, int fir, int der_iden,string modelform, List Control, IntegerVector KeepConstant, int term_tot);

void Stress_Test(IntegerVector Term_n, StringVector tform, NumericVector a_n,IntegerVector dfc,NumericMatrix x_all, int fir, int der_iden,string modelform, List Control, NumericMatrix df_groups, NumericVector tu, IntegerVector KeepConstant, int term_tot, StringVector test_point);

NumericVector cox_ph_risk_sub(IntegerVector Term_n, StringVector tform, NumericVector a_n,IntegerVector dfc,NumericMatrix x_all, int fir,string modelform, List Control, int term_tot);

void Write_Time_Dep(const NumericMatrix df0_Times, const NumericMatrix df0_dep, const NumericMatrix df0_const, const NumericVector df0_event,double dt, string filename, StringVector tform, NumericVector tu, bool iscox);

NumericMatrix Gen_Fac_Par(const NumericMatrix df0, const NumericVector vals, const NumericVector cols, const int nthreads);

bool risk_check_transition(IntegerVector Term_n, StringVector tform, NumericVector a_n,IntegerVector dfc,NumericMatrix x_all, int fir,string modelform, List Control, IntegerVector KeepConstant, int term_tot);

void Gen_Strat_Weight(string modelform, const MatrixXd& dfs, const MatrixXd& PyrC, VectorXd& s_weights, const int nthreads, const StringVector& tform, const IntegerVector& Term_n, const int& term_tot);
