//  Copyright 2022 - 2025, Eric Giunta and the project collaborators, Please see main R package for license and usage details

#ifndef SRC_GROUPING_H_
#define SRC_GROUPING_H_

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

void Make_Strata(NumericVector& Strata_vals, const Ref<const MatrixXd>& dfs, vector<vector<int> >& RiskPairs_Strata, const int& nthreads);

void Make_Groups_Strata_CR(const int& ntime, const Ref<const MatrixXd>& df_m, IntegerMatrix& RiskFail, vector<vector<vector<int> > >& RiskPairs_Strata, NumericVector& tu, const int& nthreads, NumericVector& Strata_vals, const VectorXd& cens_weight);

void Make_Match(List& model_bool, const Ref<const MatrixXd>& df_m, IntegerMatrix& RiskFail, vector<vector<int> >& RiskPairs, vector<vector<double> >& Recur_Base, vector<vector<vector<double> > >& Recur_First, vector<vector<vector<double> > >& Recur_Second, vector<double>& strata_odds, vector<int>& strata_cond, const int& nthreads);

void Make_Match_Strata(List& model_bool, const Ref<const MatrixXd>& df_m, IntegerMatrix& RiskFail, vector<vector<int> >& RiskPairs, vector<vector<double> >& Recur_Base, vector<vector<vector<double> > >& Recur_First, vector<vector<vector<double> > >& Recur_Second, vector<double>& strata_odds, vector<int>& strata_cond, const int& nthreads, NumericVector& Strata_vals);

void Make_Match_Time(List& model_bool, const int& ntime, const Ref<const MatrixXd>& df_m, IntegerMatrix& RiskFail, vector<vector<int> >& RiskPairs, vector<vector<double> >& Recur_Base, vector<vector<vector<double> > >& Recur_First, vector<vector<vector<double> > >& Recur_Second, vector<double>& strata_odds, vector<int>& strata_cond, const int& nthreads, NumericVector& tu);

void Make_Match_Time_Strata(List& model_bool, const int& ntime, const Ref<const MatrixXd>& df_m, IntegerMatrix& RiskFail, vector<vector<int> >& RiskPairs, vector<vector<double> >& Recur_Base, vector<vector<vector<double> > >& Recur_First, vector<vector<vector<double> > >& Recur_Second, vector<double>& strata_odds, vector<int>& strata_cond, const int& nthreads, NumericVector& tu, NumericVector& Strata_vals);


#endif  // SRC_GROUPING_H_
