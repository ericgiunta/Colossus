//  Copyright 2022 - 2025, Eric Giunta and the project collaborators, Please see main R package for license and usage details

#ifndef SRC_STEP_NEWTON_H_
#define SRC_STEP_NEWTON_H_

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
using Rcpp::IntegerVector;
using Rcpp::NumericVector;
using Rcpp::NumericMatrix;
using Rcpp::StringVector;
using Rcpp::List;

template <typename T> int sign(T val);

void Calc_Change_Cons(const MatrixXd& Lin_Sys, const VectorXd& Lin_Res, const  VectorXd& beta_0, const int& nthreads, const int& totalnum, const double& thres_step_max, const double& lr, const double& step_max, const vector<double>& Ll, const vector<double>& Lld, const vector<double>& Lldd, vector<double>& dbeta, const StringVector&   tform, const double& dint, const double& dslp, IntegerVector KeepConstant);

void Calc_Change(const int& nthreads, const int& totalnum, const double& thres_step_max, const double& lr, const double& step_max, const vector<double>& Ll, const vector<double>& Lld, const vector<double>& Lldd, vector<double>& dbeta, const StringVector& tform, const double& dint, const double& dslp, IntegerVector KeepConstant);

void Calc_Change_Basic(const int& nthreads, const int& totalnum, const double& lr, const double& step_max, const vector<double>& Ll, const vector<double>& Lld, const vector<double>& Lldd, vector<double>& dbeta, IntegerVector KeepConstant);

void Calc_Change_Basic_Cons(const MatrixXd& Lin_Sys, const VectorXd& Lin_Res, const  VectorXd& beta_0, const int& nthreads, const int& totalnum, const double& lr, const double& step_max, const vector<double>& Ll, const vector<double>& Lld, const vector<double>& Lldd, vector<double>& dbeta, IntegerVector KeepConstant);

void Calc_Change_Background(const int& nthreads, const int& totalnum, const int& group_num, const double& thres_step_max, const double& lr, const double& step_max, const vector<double>& Ll, const vector<double>& Lld, const vector<double>& Lldd, vector<double>& dbeta, const StringVector& tform, const double& dint, const double& dslp, IntegerVector KeepConstant, vector<int>& strata_cond, vector<double>& LldOdds, vector<double>& LlddOdds, vector<double>& LlddOddsBeta, vector<double>& dstrata);

#endif  //  SRC_STEP_NEWTON_H_
