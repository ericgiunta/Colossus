//  Copyright 2022 - 2025, Eric Giunta and the project collaborators, Please see main R package for license and usage details

#ifndef SRC_STEP_CALC_H_
#define SRC_STEP_CALC_H_

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

void Intercept_Bound(const int& nthreads, const int& totalnum, const VectorXd& beta_0, vector<double>& dbeta, const IntegerVector& dfc, const Ref<const MatrixXd>& df0, const IntegerVector& KeepConstant, const StringVector&  tform);

void Log_Bound(double& deriv_max, const MatrixXd& Lldd_mat, const VectorXd& Lld_vec, const double& Lstar, const double& qchi, const double& L0, const int& para_number, const int& nthreads, const int& totalnum, const int& reqrdnum, IntegerVector KeepConstant, const int& term_tot, const int& step, vector<double>& dbeta, const VectorXd& beta_0, bool upper, bool& trouble, int verbose, double mult);

void Calc_Change_trouble(const int& para_number, const int& nthreads, const int& totalnum, const double& thres_step_max, const double& lr, const double& step_max, const vector<double>& Ll, const vector<double>& Lld, const vector<double>& Lldd, vector<double>& dbeta, const StringVector&   tform, const double& dint, const double& dslp, IntegerVector KeepConstant_trouble);

#endif  //  SRC_STEP_CALC_H_
