using namespace std;
using namespace Rcpp;

using Eigen::Map;
using Eigen::MatrixXd;
using Eigen::SparseMatrix;
using Eigen::VectorXd;
using Rcpp::as;

//----------------------------------------------------------------------------//
List peanut_transition(IntegerVector Term_n, StringVector tform, NumericVector a_n,IntegerVector dfc,NumericMatrix x_all, int fir, int der_iden,string modelform, List Control, NumericMatrix df_groups, NumericVector tu, IntegerVector KeepConstant, int term_tot);

//  The LogLik function is directly called from transition
//  This function converts the R input into a matrix format that c++ is more efficient in
//  Currently return the final beta vector
List LogLik_PEANUT( IntegerVector Term_n, StringVector tform, NumericVector a_n,NumericMatrix x_all,IntegerVector dfc,int fir, int der_iden,string modelform, double lr, int maxiter, int halfmax, double epsilon, double dbeta_cap, double abs_max,double dose_abs_max, double deriv_epsilon, NumericMatrix df_groups, NumericVector tu ,bool change_all, bool verbose, bool debugging, IntegerVector KeepConstant, int term_tot);

template <typename T> int sign(T val);

void removeRow(Eigen::MatrixXd& matrix, unsigned int rowToRemove);
void removeColumn(Eigen::MatrixXd& matrix, unsigned int colToRemove);
//----------------------------------------------------------------------------//
void Make_Subterms(const int& totalnum, const IntegerVector& Term_n,const StringVector&  tform, const IntegerVector& dfc, const int& fir, MatrixXd& T0, MatrixXd& Td0, MatrixXd& Tdd0, MatrixXd& Dose, MatrixXd& nonDose, MatrixXd& nonDose_LIN, MatrixXd& nonDose_PLIN, MatrixXd& nonDose_LOGLIN,const  VectorXd& beta_0,const  MatrixXd& df0,const double& dint, const int& nthreads, bool debugging);

void Make_Risks(string modelform, const StringVector& tform, const IntegerVector& Term_n, const int& totalnum, const int& fir, const MatrixXd& T0, const MatrixXd& Td0, const MatrixXd& Tdd0, MatrixXd& Te, MatrixXd& R, MatrixXd& Rd, MatrixXd& Rdd, MatrixXd& Dose, MatrixXd& nonDose,  MatrixXd& nonDose_LIN, MatrixXd& nonDose_PLIN, MatrixXd& nonDose_LOGLIN, MatrixXd& RdR, MatrixXd& RddR, const int& nthreads, bool debugging);

void Make_Groups(const int& ntime, const MatrixXd& df_m, IntegerMatrix& RiskFail, vector<string>&  RiskGroup,  NumericVector& tu, const int& nthreads, bool debugging );

void Calculate_Sides(const IntegerMatrix& RiskFail, const vector<string>&  RiskGroup, const int& totalnum, const int& ntime, const MatrixXd& R, const MatrixXd& Rd, const MatrixXd& Rdd, MatrixXd& Rls1, MatrixXd& Rls2, MatrixXd& Rls3, MatrixXd& Lls1, MatrixXd& Lls2, MatrixXd& Lls3,const int& nthreads, bool debugging);

void Calc_LogLik(const int& nthreads,const IntegerMatrix& RiskFail, const vector<string>&  RiskGroup, const int& totalnum, const int& ntime, const MatrixXd& R, const MatrixXd& Rd, const MatrixXd& Rdd, const MatrixXd& RdR, const MatrixXd& RddR,const MatrixXd& Rls1,const MatrixXd& Rls2,const MatrixXd& Rls3,const MatrixXd& Lls1,const MatrixXd& Lls2,const MatrixXd& Lls3, vector<double>& Ll, vector<double>& Lld, vector<double>& Lldd, bool debugging);

void Amfit_LogLik(const int& nthreads, const int& totalnum, const MatrixXd& PyrC, const MatrixXd& R, const MatrixXd& Rd, const MatrixXd& Rdd, const MatrixXd& RdR, const MatrixXd& RddR, vector<double>& Ll, vector<double>& Lld, vector<double>& Lldd, bool debugging);

void Calc_Change(const int& nthreads, const int& totalnum,const int& fir, const int& der_iden, const double& dbeta_cap, const double& dose_abs_max, const double& lr, const double& abs_max, const vector<double>& Ll, const vector<double>& Lld, const vector<double>& Lldd, vector<double>& dbeta, const bool change_all,const StringVector&  tform, const double& dint, IntegerVector KeepConstant, bool debugging);
