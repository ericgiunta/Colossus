using namespace std;
using namespace Rcpp;

using Eigen::Map;
using Eigen::MatrixXd;
using Eigen::SparseMatrix;
using Eigen::VectorXd;
using Rcpp::as;


void Cox_Refresh_R_TERM(const int& totalnum, const int& reqrdnum, const int& term_tot, double& dint, double& dslp,double& dose_abs_max, double& abs_max, const MatrixXd& df0, MatrixXd& T0, MatrixXd& Td0, MatrixXd& Tdd0, MatrixXd& Te, MatrixXd& R, MatrixXd& Rd, MatrixXd& Rdd, MatrixXd& Dose, MatrixXd& nonDose,  MatrixXd& TTerm,  MatrixXd& nonDose_LIN, MatrixXd& nonDose_PLIN, MatrixXd& nonDose_LOGLIN, MatrixXd& RdR, MatrixXd& RddR, bool basic_bool, bool single_bool);

void Cox_Refresh_R_SIDES( const int& reqrdnum, const int& ntime, MatrixXd& Rls1, MatrixXd& Rls2, MatrixXd& Rls3, MatrixXd& Lls1, MatrixXd& Lls2, MatrixXd& Lls3, NumericVector& STRATA_vals, bool strata_bool, bool single_boo);

void Cox_Term_Risk_Calc(string modelform, const StringVector& tform, const IntegerVector& Term_n, const int& totalnum, const int& fir, const IntegerVector& dfc, int term_tot, MatrixXd& T0, MatrixXd& Td0, MatrixXd& Tdd0, MatrixXd& Te, MatrixXd& R, MatrixXd& Rd, MatrixXd& Rdd, MatrixXd& Dose, MatrixXd& nonDose, VectorXd beta_0,const  MatrixXd& df0,const double& dint, const double& dslp,  MatrixXd& TTerm,  MatrixXd& nonDose_LIN, MatrixXd& nonDose_PLIN, MatrixXd& nonDose_LOGLIN, MatrixXd& RdR, MatrixXd& RddR, const int& nthreads, bool debugging, const IntegerVector& KeepConstant, bool verbose, bool basic_bool, bool single_bool, int start, const double gmix_theta, const IntegerVector& gmix_term);

void Cox_Side_LL_Calc(const int& reqrdnum, const int& ntime, const IntegerMatrix& RiskFail, const StringMatrix&  RiskGroup_Strata, const vector<string>& RiskGroup, const int& totalnum, const int& fir, MatrixXd& R, MatrixXd& Rd, MatrixXd& Rdd,  MatrixXd& Rls1, MatrixXd& Rls2, MatrixXd& Rls3, MatrixXd& Lls1, MatrixXd& Lls2, MatrixXd& Lls3, const VectorXd& cens_weight, NumericVector& STRATA_vals, VectorXd beta_0 , MatrixXd& RdR, MatrixXd& RddR, vector<double>& Ll, vector<double>& Lld, vector<double>& Lldd, const int& nthreads, bool debugging, const IntegerVector& KeepConstant,string ties_method, bool verbose,bool strata_bool, bool CR_bool, bool basic_bool, bool single_bool, int start, int iter_stop);

void Pois_Term_Risk_Calc(string modelform, const StringVector& tform, const IntegerVector& Term_n, const int& totalnum, const int& fir, const IntegerVector& dfc, int term_tot, MatrixXd& T0, MatrixXd& Td0, MatrixXd& Tdd0, MatrixXd& Te, MatrixXd& R, MatrixXd& Rd, MatrixXd& Rdd, MatrixXd& Dose, MatrixXd& nonDose, VectorXd beta_0,const  MatrixXd& df0,const double& dint, const double& dslp,  MatrixXd& TTerm,  MatrixXd& nonDose_LIN, MatrixXd& nonDose_PLIN, MatrixXd& nonDose_LOGLIN, MatrixXd& RdR, MatrixXd& RddR, const MatrixXd& s_weights, const int& nthreads, bool debugging, const IntegerVector& KeepConstant, bool verbose, bool strata_bool, bool single_bool, int start, const double gmix_theta, const IntegerVector& gmix_term);

void Pois_Dev_LL_Calc(const int& reqrdnum, const int& totalnum, const int& fir, MatrixXd& R, MatrixXd& Rd, MatrixXd& Rdd, VectorXd beta_0 , MatrixXd& RdR, MatrixXd& RddR, vector<double>& Ll, vector<double>& Lld, vector<double>& Lldd, const MatrixXd& PyrC, MatrixXd& dev_temp, const int& nthreads, bool debugging, const IntegerVector& KeepConstant, bool verbose, bool single_bool, int start, int iter_stop, double& dev);

void Log_Bound(const MatrixXd& Lldd_mat, const VectorXd& Lld_vec, const double& Lstar, const double& qchi, const double& L0, const int& para_number, const int& nthreads, const int& totalnum, const int& reqrdnum, IntegerVector KeepConstant, const int& term_tot, const int& step, vector<double>& dbeta, const VectorXd& beta_0, bool upper, bool verbose);

void Simplified_Inform_Matrix(const int& nthreads,const IntegerMatrix& RiskFail, const vector<string>&  RiskGroup, const int& totalnum, const int& ntime, const MatrixXd& R, const MatrixXd& Rd, const MatrixXd& Rdd, const MatrixXd& RdR, const MatrixXd& RddR,const MatrixXd& Rls1,const MatrixXd& Rls2,const MatrixXd& Rls3,const MatrixXd& Lls1,const MatrixXd& Lls2,const MatrixXd& Lls3, vector<double>& InMa, bool debugging,string ties_method, const IntegerVector& KeepConstant);
void Simplified_Inform_Matrix_STRATA(const int& nthreads,const IntegerMatrix& RiskFail, const StringMatrix&  RiskGroup, const int& totalnum, const int& ntime, const MatrixXd& R, const MatrixXd& Rd, const MatrixXd& Rdd, const MatrixXd& RdR, const MatrixXd& RddR,const MatrixXd& Rls1,const MatrixXd& Rls2,const MatrixXd& Rls3,const MatrixXd& Lls1,const MatrixXd& Lls2,const MatrixXd& Lls3, vector<double>& InMa, bool debugging,string ties_method, NumericVector& STRATA_vals, const IntegerVector& KeepConstant);
