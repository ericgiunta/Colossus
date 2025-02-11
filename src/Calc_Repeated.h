using namespace std;
using namespace Rcpp;

using Eigen::Map;
using Eigen::MatrixXd;
using Eigen::SparseMatrix;
using Eigen::VectorXd;
using Rcpp::as;

template <typename T> int sign(T val);

void removeRow(Eigen::MatrixXd& matrix, unsigned int rowToRemove);

void removeColumn(Eigen::MatrixXd& matrix, unsigned int colToRemove);

void Make_Groups(const int& ntime, const MatrixXd& df_m, IntegerMatrix& RiskFail, vector<vector<int> >& RiskPairs, NumericVector& tu, const int& nthreads);

void Make_Groups_CR(const int& ntime, const MatrixXd& df_m, IntegerMatrix& RiskFail, vector<vector<int> >& RiskPairs, NumericVector& tu, const VectorXd& cens_weight, const int& nthreads);

void Make_Groups_Strata(const int& ntime, const MatrixXd& df_m, IntegerMatrix& RiskFail, vector<vector<vector<int> > >& RiskPairs_Strata, NumericVector& tu, const int& nthreads, NumericVector& Strata_vals);

void Make_Groups_Strata_CR(const int& ntime, const MatrixXd& df_m, IntegerMatrix& RiskFail, vector<vector<vector<int> > >& RiskPairs_Strata, NumericVector& tu, const int& nthreads, NumericVector& Strata_vals, const VectorXd& cens_weight);

void Calculate_Sides(const IntegerMatrix& RiskFail, const vector<vector<int> >& RiskPairs, const int& totalnum, const int& ntime, const MatrixXd& R, const MatrixXd& Rd, const MatrixXd& Rdd, MatrixXd& Rls1, MatrixXd& Rls2, MatrixXd& Rls3, MatrixXd& Lls1, MatrixXd& Lls2, MatrixXd& Lls3, const int& nthreads, const IntegerVector& KeepConstant);

void Calculate_Sides_Gradient(const IntegerMatrix& RiskFail, const vector<vector<int> >& RiskPairs, const int& totalnum, const int& ntime, const MatrixXd& R, const MatrixXd& Rd, MatrixXd& Rls1, MatrixXd& Rls2, MatrixXd& Lls1, MatrixXd& Lls2, const int& nthreads, const IntegerVector& KeepConstant);

void Calculate_Sides_CR(const IntegerMatrix& RiskFail, const vector<vector<int> >& RiskPairs, const int& totalnum, const int& ntime, const MatrixXd& R, const MatrixXd& Rd, const MatrixXd& Rdd, MatrixXd& Rls1, MatrixXd& Rls2, MatrixXd& Rls3, MatrixXd& Lls1, MatrixXd& Lls2, MatrixXd& Lls3, const VectorXd& cens_weight, const int& nthreads, const IntegerVector& KeepConstant);

void Calculate_Sides_CR_Gradient(const IntegerMatrix& RiskFail, const vector<vector<int> >& RiskPairs, const int& totalnum, const int& ntime, const MatrixXd& R, const MatrixXd& Rd, MatrixXd& Rls1, MatrixXd& Rls2, MatrixXd& Lls1, MatrixXd& Lls2, const VectorXd& cens_weight, const int& nthreads, const IntegerVector& KeepConstant);

void Calculate_Sides_CR_SINGLE(const IntegerMatrix& RiskFail, const vector<vector<int> >& RiskPairs, const int& totalnum, const int& ntime, const MatrixXd& R, MatrixXd& Rls1, MatrixXd& Lls1, const VectorXd& cens_weight, const int& nthreads, const IntegerVector& KeepConstant);

void Calculate_Sides_Single(const IntegerMatrix& RiskFail, const vector<vector<int> >& RiskPairs, const int& totalnum, const int& ntime, const MatrixXd& R, MatrixXd& Rls1, MatrixXd& Lls1, const int& nthreads);

void Calculate_Sides_Strata(const IntegerMatrix& RiskFail, const vector<vector<vector<int> > >& RiskPairs_Strata, const int& totalnum, const int& ntime, const MatrixXd& R, const MatrixXd& Rd, const MatrixXd& Rdd, MatrixXd& Rls1, MatrixXd& Rls2, MatrixXd& Rls3, MatrixXd& Lls1, MatrixXd& Lls2, MatrixXd& Lls3, const int& nthreads, NumericVector& Strata_vals, const IntegerVector& KeepConstant);

void Calculate_Sides_Strata_Gradient(const IntegerMatrix& RiskFail, const vector<vector<vector<int> > >& RiskPairs_Strata, const int& totalnum, const int& ntime, const MatrixXd& R, const MatrixXd& Rd, MatrixXd& Rls1, MatrixXd& Rls2, MatrixXd& Lls1, MatrixXd& Lls2, const int& nthreads, NumericVector& Strata_vals, const IntegerVector& KeepConstant);

void Calculate_Sides_Strata_Single(const IntegerMatrix& RiskFail, const vector<vector<vector<int> > >& RiskPairs_Strata, const int& totalnum, const int& ntime, const MatrixXd& R, MatrixXd& Rls1, MatrixXd& Lls1, const int& nthreads, NumericVector& Strata_vals, const IntegerVector& KeepConstant);

void Calculate_Sides_Strata_CR(const IntegerMatrix& RiskFail, const vector<vector<vector<int> > >& RiskPairs_Strata, const int& totalnum, const int& ntime, const MatrixXd& R, const MatrixXd& Rd, const MatrixXd& Rdd, MatrixXd& Rls1, MatrixXd& Rls2, MatrixXd& Rls3, MatrixXd& Lls1, MatrixXd& Lls2, MatrixXd& Lls3, const VectorXd& cens_weight, const int& nthreads, NumericVector& Strata_vals, const IntegerVector& KeepConstant);

void Calculate_Sides_Strata_CR_Gradient(const IntegerMatrix& RiskFail, const vector<vector<vector<int> > >& RiskPairs_Strata, const int& totalnum, const int& ntime, const MatrixXd& R, const MatrixXd& Rd, MatrixXd& Rls1, MatrixXd& Rls2, MatrixXd& Lls1, MatrixXd& Lls2, const VectorXd& cens_weight, const int& nthreads, NumericVector& Strata_vals, const IntegerVector& KeepConstant);

void Calculate_Sides_Strata_Single_CR(const IntegerMatrix& RiskFail, const vector<vector<vector<int> > >& RiskPairs_Strata, const int& totalnum, const int& ntime, const MatrixXd& R, MatrixXd& Rls1, MatrixXd& Lls1, const VectorXd& cens_weight, const int& nthreads, NumericVector& Strata_vals, const IntegerVector& KeepConstant);

void Calc_LogLik(const int& nthreads, const IntegerMatrix& RiskFail, const vector<vector<int> >& RiskPairs, const int& totalnum, const int& ntime, const MatrixXd& R, const MatrixXd& Rd, const MatrixXd& Rdd, const MatrixXd& RdR, const MatrixXd& RddR, const MatrixXd& Rls1, const MatrixXd& Rls2, const MatrixXd& Rls3, const MatrixXd& Lls1, const MatrixXd& Lls2, const MatrixXd& Lls3, vector<double>& Ll, vector<double>& Lld, vector<double>& Lldd, string ties_method, const IntegerVector& KeepConstant);

void Calc_LogLik_Gradient(const int& nthreads, const IntegerMatrix& RiskFail, const vector<vector<int> >& RiskPairs, const int& totalnum, const int& ntime, const MatrixXd& R, const MatrixXd& Rd, const MatrixXd& RdR, const MatrixXd& Rls1, const MatrixXd& Rls2, const MatrixXd& Lls1, const MatrixXd& Lls2, vector<double>& Ll, vector<double>& Lld, string ties_method, const IntegerVector& KeepConstant);

void Calc_LogLik_Single(const int& nthreads, const IntegerMatrix& RiskFail, const vector<vector<int> >& RiskPairs, const int& totalnum, const int& ntime, const MatrixXd& R, const MatrixXd& Rls1, const MatrixXd& Lls1, vector<double>& Ll, string ties_method);

void Calc_LogLik_Strata(const int& nthreads, const IntegerMatrix& RiskFail, const vector<vector<vector<int> > >& RiskPairs_Strata, const int& totalnum, const int& ntime, const MatrixXd& R, const MatrixXd& Rd, const MatrixXd& Rdd, const MatrixXd& RdR, const MatrixXd& RddR, const MatrixXd& Rls1, const MatrixXd& Rls2, const MatrixXd& Rls3, const MatrixXd& Lls1, const MatrixXd& Lls2, const MatrixXd& Lls3, vector<double>& Ll, vector<double>& Lld, vector<double>& Lldd, string ties_method, NumericVector& Strata_vals, const IntegerVector& KeepConstant);

void Calc_LogLik_Strata_Gradient(const int& nthreads, const IntegerMatrix& RiskFail, const vector<vector<vector<int> > >& RiskPairs_Strata, const int& totalnum, const int& ntime, const MatrixXd& R, const MatrixXd& Rd, const MatrixXd& RdR, const MatrixXd& Rls1, const MatrixXd& Rls2, const MatrixXd& Lls1, const MatrixXd& Lls2, vector<double>& Ll, vector<double>& Lld, string ties_method, NumericVector& Strata_vals, const IntegerVector& KeepConstant);

void Calc_LogLik_Strata_SINGLE(const int& nthreads, const IntegerMatrix& RiskFail, const vector<vector<vector<int> > >& RiskPairs_Strata, const int& totalnum, const int& ntime, const MatrixXd& R, const MatrixXd& Rls1, const MatrixXd& Lls1, vector<double>& Ll, string ties_method, NumericVector& Strata_vals, const IntegerVector& KeepConstant);

void Calc_LogLik_Strata_BASIC(const int& nthreads, const IntegerMatrix& RiskFail, const vector<vector<vector<int> > >& RiskPairs_Strata, const int& totalnum, const int& ntime, const MatrixXd& R, const MatrixXd& Rd, const MatrixXd& Rdd, const MatrixXd& RdR, const MatrixXd& Rls1, const MatrixXd& Rls2, const MatrixXd& Rls3, const MatrixXd& Lls1, const MatrixXd& Lls2, const MatrixXd& Lls3, vector<double>& Ll, vector<double>& Lld, vector<double>& Lldd, string ties_method, NumericVector& Strata_vals, const IntegerVector& KeepConstant);

void Calc_LogLik_Basic(const int& nthreads, const IntegerMatrix& RiskFail, const vector<vector<int> >& RiskPairs, const int& totalnum, const int& ntime, const MatrixXd& R, const MatrixXd& Rd, const MatrixXd& Rdd, const MatrixXd& RdR, const MatrixXd& Rls1, const MatrixXd& Rls2, const MatrixXd& Rls3, const MatrixXd& Lls1, const MatrixXd& Lls2, const MatrixXd& Lls3, vector<double>& Ll, vector<double>& Lld, vector<double>& Lldd, string ties_method, const IntegerVector& KeepConstant);

void Calc_LogLik_Linear_ERR(const StringVector& tform, const int& nthreads, const IntegerMatrix& RiskFail, const vector<vector<int> >& RiskPairs, const int& totalnum, const int& ntime, const MatrixXd& R, const MatrixXd& Rd, const MatrixXd& Rdd, const MatrixXd& RdR, const MatrixXd& RddR, const MatrixXd& Rls1, const MatrixXd& Rls2, const MatrixXd& Rls3, const MatrixXd& Lls1, const MatrixXd& Lls2, const MatrixXd& Lls3, vector<double>& Ll, vector<double>& Lld, vector<double>& Lldd, string ties_method, const IntegerVector& KeepConstant);

void Calc_LogLik_Strata_Linear_ERR(const StringVector& tform, const int& nthreads, const IntegerMatrix& RiskFail, const vector<vector<vector<int> > >& RiskPairs_Strata, const int& totalnum, const int& ntime, const MatrixXd& R, const MatrixXd& Rd, const MatrixXd& Rdd, const MatrixXd& RdR, const MatrixXd& RddR, const MatrixXd& Rls1, const MatrixXd& Rls2, const MatrixXd& Rls3, const MatrixXd& Lls1, const MatrixXd& Lls2, const MatrixXd& Lls3, vector<double>& Ll, vector<double>& Lld, vector<double>& Lldd, string ties_method, NumericVector& Strata_vals, const IntegerVector& KeepConstant);

void Poisson_LogLik(const int& nthreads, const int& totalnum, const MatrixXd& PyrC, const MatrixXd& R, const MatrixXd& Rd, const MatrixXd& Rdd, const MatrixXd& RdR, const MatrixXd& RddR, vector<double>& Ll, vector<double>& Lld, vector<double>& Lldd, const IntegerVector& KeepConstant);

void Poisson_LogLik_Gradient(const int& nthreads, const int& totalnum, const MatrixXd& PyrC, const MatrixXd& R, const MatrixXd& Rd, const MatrixXd& RdR, vector<double>& Ll, vector<double>& Lld, const IntegerVector& KeepConstant);

void Poisson_LogLik_Single(const int& nthreads, const int& totalnum, const MatrixXd& PyrC, const MatrixXd& R, vector<double>& Ll);

void Calculate_Null_Sides(const IntegerMatrix& RiskFail, const vector<vector<int> >& RiskPairs, const int& ntime, const MatrixXd& R, MatrixXd& Rls1, MatrixXd& Lls1, const int& nthreads);

void Calc_Null_LogLik(const int& nthreads, const IntegerMatrix& RiskFail, const vector<vector<int> >& RiskPairs, const int& ntime, const MatrixXd& R, const MatrixXd& Rls1, const MatrixXd& Lls1, vector<double>& Ll, string ties_method);

void Calculate_Null_Sides_Strata(const IntegerMatrix& RiskFail, const vector<vector<vector<int> > >& RiskPairs_Strata, const int& ntime, const MatrixXd& R, MatrixXd& Rls1, MatrixXd& Lls1, NumericVector& Strata_vals, const int& nthreads);

void Calc_Null_LogLik_Strata(const int& nthreads, const IntegerMatrix& RiskFail, const vector<vector<vector<int> > >& RiskPairs_Strata, const int& ntime, const MatrixXd& R, const MatrixXd& Rls1, const MatrixXd& Lls1, NumericVector& Strata_vals, vector<double>& Ll, string ties_method);
