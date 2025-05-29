using namespace std;
using namespace Rcpp;

using Eigen::Map;
using Eigen::MatrixXd;
using Eigen::SparseMatrix;
using Eigen::VectorXd;
using Rcpp::as;

template <typename T> int sign(T val);



void Intercept_Bound(const int& nthreads, const int& totalnum, const VectorXd& beta_0, vector<double>& dbeta, const IntegerVector& dfc, const  MatrixXd& df0, const IntegerVector& KeepConstant, const StringVector&  tform);

void Calc_Change_Cons(const MatrixXd& Lin_Sys, const VectorXd& Lin_Res, const  VectorXd& beta_0, const int& nthreads, const int& totalnum, const double& dose_abs_max, const double& lr, const double& abs_max, const vector<double>& Ll, const vector<double>& Lld, const vector<double>& Lldd, vector<double>& dbeta, const StringVector&   tform, const double& dint, const double& dslp, IntegerVector KeepConstant);

void Calc_Change(const int& double_step, const int& nthreads, const int& totalnum, const double& dose_abs_max, const double& lr, const double& abs_max, const vector<double>& Ll, const vector<double>& Lld, const vector<double>& Lldd, vector<double>& dbeta, const StringVector& tform, const double& dint, const double& dslp, IntegerVector KeepConstant);

void Calc_Change_Gradient(const int& nthreads, List& model_bool, const int& totalnum, List& optim_para, int& iteration, const double& abs_max, const vector<double>& Lld, NumericVector& m_g_store, NumericVector& v_beta_store, vector<double>& dbeta, IntegerVector KeepConstant);

void Calc_Change_Basic(const int& double_step, const int& nthreads, const int& totalnum, const double& lr, const double& abs_max, const vector<double>& Ll, const vector<double>& Lld, const vector<double>& Lldd, vector<double>& dbeta, IntegerVector KeepConstant);

void Log_Bound(double& deriv_max, const MatrixXd& Lldd_mat, const VectorXd& Lld_vec, const double& Lstar, const double& qchi, const double& L0, const int& para_number, const int& nthreads, const int& totalnum, const int& reqrdnum, IntegerVector KeepConstant, const int& term_tot, const int& step, vector<double>& dbeta, const VectorXd& beta_0, bool upper, bool& trouble, int verbose, double mult);

void Calc_Change_trouble(const int& para_number, const int& nthreads, const int& totalnum, const double& dose_abs_max, const double& lr, const double& abs_max, const vector<double>& Ll, const vector<double>& Lld, const vector<double>& Lldd, vector<double>& dbeta, const StringVector&   tform, const double& dint, const double& dslp, IntegerVector KeepConstant_trouble);

void Calc_Change_Background(const int& double_step, const int& nthreads, const int& totalnum, const int& group_num, const double& dose_abs_max, const double& lr, const double& abs_max, const vector<double>& Ll, const vector<double>& Lld, const vector<double>& Lldd, vector<double>& dbeta, const StringVector& tform, const double& dint, const double& dslp, IntegerVector KeepConstant, vector<int>& strata_cond, vector<double>& LldOdds, vector<double>& LlddOdds, vector<double>& LlddOddsBeta, vector<double>& dstrata);

void Calc_Change_Background_Gradient(const int& nthreads, List& model_bool, const int& totalnum, const int& group_num, List& optim_para, int& iteration, const double& abs_max, const vector<double>& Lld, NumericVector& m_g_store, NumericVector& v_beta_store, vector<double>& dbeta, IntegerVector KeepConstant, vector<int>& strata_cond, vector<double>& LldOdds, vector<double>& dstrata);
