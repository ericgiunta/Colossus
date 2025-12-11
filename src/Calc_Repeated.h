//  Copyright 2022 - 2025, Eric Giunta and the project collaborators, Please see main R package for license and usage details

#ifndef SRC_CALC_REPEATED_H_
#define SRC_CALC_REPEATED_H_

#include <string>
#include <vector>

using std::string;
using std::vector;

using Eigen::Map;
using Eigen::Ref;
using Eigen::MatrixXd;
using Eigen::SparseMatrix;
using Eigen::VectorXd;

using Rcpp::as;
using Rcpp::IntegerMatrix;
using Rcpp::IntegerVector;
using Rcpp::NumericVector;
using Rcpp::NumericMatrix;
using Rcpp::StringVector;
using Rcpp::List;

template <typename T> int sign(T val);

void Make_Groups(const int& ntime, const Ref<const MatrixXd>& df_m, IntegerMatrix& RiskFail, vector<vector<int> >& RiskPairs, NumericVector& tu, const int& nthreads);

void Make_Groups_CR(const int& ntime, const Ref<const MatrixXd>& df_m, IntegerMatrix& RiskFail, vector<vector<int> >& RiskPairs, NumericVector& tu, const VectorXd& cens_weight, const int& nthreads);

void Make_Groups_Strata(const int& ntime, const Ref<const MatrixXd>& df_m, IntegerMatrix& RiskFail, vector<vector<vector<int> > >& RiskPairs_Strata, NumericVector& tu, const int& nthreads, NumericVector& Strata_vals);

void Make_Groups_Strata_CR(const int& ntime, const Ref<const MatrixXd>& df_m, IntegerMatrix& RiskFail, vector<vector<vector<int> > >& RiskPairs_Strata, NumericVector& tu, const int& nthreads, NumericVector& Strata_vals, const VectorXd& cens_weight);

void Make_Match(List& model_bool, const Ref<const MatrixXd>& df_m, IntegerMatrix& RiskFail, vector<vector<int> >& RiskPairs, vector<vector<double> >& Recur_Base, vector<vector<vector<double> > >& Recur_First, vector<vector<vector<double> > >& Recur_Second, vector<double>& strata_odds, vector<int>& strata_cond, const int& nthreads);

void Make_Match_Strata(List& model_bool, const Ref<const MatrixXd>& df_m, IntegerMatrix& RiskFail, vector<vector<int> >& RiskPairs, vector<vector<double> >& Recur_Base, vector<vector<vector<double> > >& Recur_First, vector<vector<vector<double> > >& Recur_Second, vector<double>& strata_odds, vector<int>& strata_cond, const int& nthreads, NumericVector& Strata_vals);

void Make_Match_Time(List& model_bool, const int& ntime, const Ref<const MatrixXd>& df_m, IntegerMatrix& RiskFail, vector<vector<int> >& RiskPairs, vector<vector<double> >& Recur_Base, vector<vector<vector<double> > >& Recur_First, vector<vector<vector<double> > >& Recur_Second, vector<double>& strata_odds, vector<int>& strata_cond, const int& nthreads, NumericVector& tu);

void Make_Match_Time_Strata(List& model_bool, const int& ntime, const Ref<const MatrixXd>& df_m, IntegerMatrix& RiskFail, vector<vector<int> >& RiskPairs, vector<vector<double> >& Recur_Base, vector<vector<vector<double> > >& Recur_First, vector<vector<vector<double> > >& Recur_Second, vector<double>& strata_odds, vector<int>& strata_cond, const int& nthreads, NumericVector& tu, NumericVector& Strata_vals);

void Calculate_Sides(List& model_bool, const IntegerMatrix& RiskFail, const vector<vector<int> >& RiskPairs, const int& totalnum, const int& ntime, const MatrixXd& R, const MatrixXd& Rd, const MatrixXd& Rdd, MatrixXd& Rls1, MatrixXd& Rls2, MatrixXd& Rls3, MatrixXd& Lls1, MatrixXd& Lls2, MatrixXd& Lls3, const int& nthreads, const IntegerVector& KeepConstant);

void Calculate_Sides_PO(List& model_bool, const IntegerMatrix& RiskFail, const vector<vector<int> >& RiskPairs, const int& totalnum, const int& ntime, const MatrixXd& R, const MatrixXd& Rd, const MatrixXd& Rdd, MatrixXd& Rls1, MatrixXd& Rls2, MatrixXd& Rls3, MatrixXd& Lls1, MatrixXd& Lls2, MatrixXd& Lls3, const VectorXd& cens_weight, const int& nthreads, const IntegerVector& KeepConstant);

void Calculate_Sides_CR(List& model_bool, const IntegerMatrix& RiskFail, const vector<vector<int> >& RiskPairs, const int& totalnum, const int& ntime, const MatrixXd& R, const MatrixXd& Rd, const MatrixXd& Rdd, MatrixXd& Rls1, MatrixXd& Rls2, MatrixXd& Rls3, MatrixXd& Lls1, MatrixXd& Lls2, MatrixXd& Lls3, const VectorXd& cens_weight, const int& nthreads, const IntegerVector& KeepConstant);

void Calculate_Sides_Strata(List& model_bool, const IntegerMatrix& RiskFail, const vector<vector<vector<int> > >& RiskPairs_Strata, const int& totalnum, const int& ntime, const MatrixXd& R, const MatrixXd& Rd, const MatrixXd& Rdd, MatrixXd& Rls1, MatrixXd& Rls2, MatrixXd& Rls3, MatrixXd& Lls1, MatrixXd& Lls2, MatrixXd& Lls3, const int& nthreads, NumericVector& Strata_vals, const IntegerVector& KeepConstant);

void Calculate_Sides_Strata_CR(List& model_bool, const IntegerMatrix& RiskFail, const vector<vector<vector<int> > >& RiskPairs_Strata, const int& totalnum, const int& ntime, const MatrixXd& R, const MatrixXd& Rd, const MatrixXd& Rdd, MatrixXd& Rls1, MatrixXd& Rls2, MatrixXd& Rls3, MatrixXd& Lls1, MatrixXd& Lls2, MatrixXd& Lls3, const VectorXd& cens_weight, const int& nthreads, NumericVector& Strata_vals, const IntegerVector& KeepConstant);

void Calc_LogLik(List& model_bool, const int& nthreads, const IntegerMatrix& RiskFail, const vector<vector<int> >& RiskPairs, const int& totalnum, const int& ntime, const MatrixXd& R, const MatrixXd& Rd, const MatrixXd& Rdd, const MatrixXd& RdR, const MatrixXd& RddR, const MatrixXd& Rls1, const MatrixXd& Rls2, const MatrixXd& Rls3, const MatrixXd& Lls1, const MatrixXd& Lls2, const MatrixXd& Lls3, vector<double>& Ll, vector<double>& Lld, vector<double>& Lldd, string ties_method, const IntegerVector& KeepConstant);

void Calc_LogLik_PO(List& model_bool, const int& nthreads, const IntegerMatrix& RiskFail, const vector<vector<int> >& RiskPairs, const int& totalnum, const int& ntime, const MatrixXd& R, const MatrixXd& Rd, const MatrixXd& Rdd, const MatrixXd& RdR, const MatrixXd& RddR, const MatrixXd& Rls1, const MatrixXd& Rls2, const MatrixXd& Rls3, const MatrixXd& Lls1, const MatrixXd& Lls2, const MatrixXd& Lls3, const VectorXd& cens_weight, vector<double>& Ll, vector<double>& Lld, vector<double>& Lldd, string ties_method, const IntegerVector& KeepConstant);

void Calc_LogLik_Strata(List& model_bool, const int& nthreads, const IntegerMatrix& RiskFail, const vector<vector<vector<int> > >& RiskPairs_Strata, const int& totalnum, const int& ntime, const MatrixXd& R, const MatrixXd& Rd, const MatrixXd& Rdd, const MatrixXd& RdR, const MatrixXd& RddR, const MatrixXd& Rls1, const MatrixXd& Rls2, const MatrixXd& Rls3, const MatrixXd& Lls1, const MatrixXd& Lls2, const MatrixXd& Lls3, vector<double>& Ll, vector<double>& Lld, vector<double>& Lldd, string ties_method, NumericVector& Strata_vals, const IntegerVector& KeepConstant);

void Calc_LogLik_Strata_Basic(List& model_bool, const int& nthreads, const IntegerMatrix& RiskFail, const vector<vector<vector<int> > >& RiskPairs_Strata, const int& totalnum, const int& ntime, const MatrixXd& R, const MatrixXd& Rd, const MatrixXd& Rdd, const MatrixXd& RdR, const MatrixXd& Rls1, const MatrixXd& Rls2, const MatrixXd& Rls3, const MatrixXd& Lls1, const MatrixXd& Lls2, const MatrixXd& Lls3, vector<double>& Ll, vector<double>& Lld, vector<double>& Lldd, string ties_method, NumericVector& Strata_vals, const IntegerVector& KeepConstant);

void Calc_LogLik_Basic(List& model_bool, const int& nthreads, const IntegerMatrix& RiskFail, const vector<vector<int> >& RiskPairs, const int& totalnum, const int& ntime, const MatrixXd& R, const MatrixXd& Rd, const MatrixXd& Rdd, const MatrixXd& RdR, const MatrixXd& Rls1, const MatrixXd& Rls2, const MatrixXd& Rls3, const MatrixXd& Lls1, const MatrixXd& Lls2, const MatrixXd& Lls3, vector<double>& Ll, vector<double>& Lld, vector<double>& Lldd, string ties_method, const IntegerVector& KeepConstant);

void Calc_LogLik_Linear_ERR(List& model_bool, const StringVector& tform, const int& nthreads, const IntegerMatrix& RiskFail, const vector<vector<int> >& RiskPairs, const int& totalnum, const int& ntime, const MatrixXd& R, const MatrixXd& Rd, const MatrixXd& Rdd, const MatrixXd& RdR, const MatrixXd& RddR, const MatrixXd& Rls1, const MatrixXd& Rls2, const MatrixXd& Rls3, const MatrixXd& Lls1, const MatrixXd& Lls2, const MatrixXd& Lls3, vector<double>& Ll, vector<double>& Lld, vector<double>& Lldd, string ties_method, const IntegerVector& KeepConstant);

void Calc_LogLik_Strata_Linear_ERR(List& model_bool, const StringVector& tform, const int& nthreads, const IntegerMatrix& RiskFail, const vector<vector<vector<int> > >& RiskPairs_Strata, const int& totalnum, const int& ntime, const MatrixXd& R, const MatrixXd& Rd, const MatrixXd& Rdd, const MatrixXd& RdR, const MatrixXd& RddR, const MatrixXd& Rls1, const MatrixXd& Rls2, const MatrixXd& Rls3, const MatrixXd& Lls1, const MatrixXd& Lls2, const MatrixXd& Lls3, vector<double>& Ll, vector<double>& Lld, vector<double>& Lldd, string ties_method, NumericVector& Strata_vals, const IntegerVector& KeepConstant);

void Poisson_LogLik(List& model_bool, const int& nthreads, const int& totalnum, const Ref<const MatrixXd>& PyrC, const MatrixXd& R, const MatrixXd& Rd, const MatrixXd& Rdd, const MatrixXd& RdR, const MatrixXd& RddR, vector<double>& Ll, vector<double>& Lld, vector<double>& Lldd, const IntegerVector& KeepConstant);

void Calculate_Null_Sides(const IntegerMatrix& RiskFail, const vector<vector<int> >& RiskPairs, const int& ntime, const MatrixXd& R, MatrixXd& Rls1, MatrixXd& Lls1, const int& nthreads);

void Calc_Null_LogLik(const int& nthreads, const IntegerMatrix& RiskFail, const vector<vector<int> >& RiskPairs, const int& ntime, const MatrixXd& R, const MatrixXd& Rls1, const MatrixXd& Lls1, vector<double>& Ll, string ties_method);

void Calculate_Null_Sides_Strata(const IntegerMatrix& RiskFail, const vector<vector<vector<int> > >& RiskPairs_Strata, const int& ntime, const MatrixXd& R, MatrixXd& Rls1, MatrixXd& Lls1, NumericVector& Strata_vals, const int& nthreads);

void Calc_Null_LogLik_Strata(const int& nthreads, const IntegerMatrix& RiskFail, const vector<vector<vector<int> > >& RiskPairs_Strata, const int& ntime, const MatrixXd& R, const MatrixXd& Rls1, const MatrixXd& Lls1, NumericVector& Strata_vals, vector<double>& Ll, string ties_method);

void Calculate_Recursive(List& model_bool, const int& group_num, const IntegerMatrix& RiskFail, const vector<vector<int> >& RiskPairs, const int& totalnum, const int& ntime, const MatrixXd& R, const MatrixXd& Rd, const MatrixXd& Rdd, vector<vector<double> >& Recur_Base, vector<vector<vector<double> > >& Recur_First, vector<vector<vector<double> > >& Recur_Second, const int& nthreads, const IntegerVector& KeepConstant);

void Calc_Recur_LogLik(List& model_bool, const int& group_num, const IntegerMatrix& RiskFail, const vector<vector<int> >& RiskPairs, const int& totalnum, const int& ntime, const MatrixXd& R, const MatrixXd& Rd, const MatrixXd& Rdd, const MatrixXd& RdR, const MatrixXd& RddR, double& dev, vector<double>& Ll, vector<double>& Lld, vector<double>& Lldd, vector<vector<double> >& Recur_Base, vector<vector<vector<double> > >& Recur_First, vector<vector<vector<double> > >& Recur_Second, vector<double>& strata_odds, const int& nthreads, const IntegerVector& KeepConstant, vector<int>& strata_cond, vector<double>& LldOdds, vector<double>& LlddOdds, vector<double>& LlddOddsBeta);

void Calc_LogLik_Logist(List& model_bool, const int& nthreads, const int& totalnum, const Ref<const MatrixXd>& CountEvent, const MatrixXd& P, const MatrixXd& Pnot, const MatrixXd& Pd, const MatrixXd& Pdd, const MatrixXd& PdP, const MatrixXd& PnotdP, const MatrixXd& PddP, const MatrixXd& PnotddP, vector<double>& Ll, vector<double>& Lld, vector<double>& Lldd, const IntegerVector& KeepConstant);

#endif  // SRC_CALC_REPEATED_H_
