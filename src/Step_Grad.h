//  Copyright 2022 - 2025, Eric Giunta and the project collaborators, Please see main R package for license and usage details

#ifndef SRC_STEP_GRAD_H_
#define SRC_STEP_GRAD_H_

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

void Calc_Change_Gradient(const int& nthreads, List& model_bool, const int& totalnum, List& optim_para, int& iteration, const double& step_max, const vector<double>& Lld, NumericVector& m_g_store, NumericVector& v_beta_store, vector<double>& dbeta, IntegerVector KeepConstant);

void Calc_Change_Gradient_Cons(const MatrixXd& Lin_Sys, const VectorXd& Lin_Res, const int& nthreads, List& model_bool, const int& totalnum, List& optim_para, int& iteration, const double& step_max, const vector<double>& Ll,  const vector<double>& Lld, NumericVector& m_g_store, NumericVector& v_beta_store, const VectorXd& beta_0, vector<double>& dbeta, IntegerVector KeepConstant);

void Calc_Change_Background_Gradient(const int& nthreads, List& model_bool, const int& totalnum, const int& group_num, List& optim_para, int& iteration, const double& step_max, const vector<double>& Lld, NumericVector& m_g_store, NumericVector& v_beta_store, vector<double>& dbeta, IntegerVector KeepConstant, vector<int>& strata_cond, vector<double>& LldOdds, vector<double>& dstrata);

#endif  //  SRC_STEP_GRAD_H_
