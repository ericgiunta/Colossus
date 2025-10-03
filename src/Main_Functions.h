//  Copyright 2022 - 2025, Eric Giunta and the project collaborators, Please see main R package for license and usage details

#ifndef SRC_MAIN_FUNCTIONS_H_
#define SRC_MAIN_FUNCTIONS_H_

#include <string>

using std::string;
using std::vector;

using Eigen::Map;
using Eigen::MatrixXd;
using Eigen::SparseMatrix;
using Eigen::VectorXd;

using Rcpp::as;
using Rcpp::IntegerVector;
using Rcpp::NumericVector;
using Rcpp::NumericMatrix;
using Rcpp::StringVector;
using Rcpp::List;

bool Check_Risk(IntegerVector term_n, StringVector tform, NumericVector a_n, NumericMatrix& x_all, IntegerVector dfc, int fir, string modelform, int verbose, IntegerVector KeepConstant, int term_tot, int nthreads, const double gmix_theta, const IntegerVector gmix_term);


List LogLik_Cox_PH_Omnibus(IntegerVector term_n, StringVector tform, NumericMatrix& a_ns, NumericMatrix& x_all, IntegerVector dfc, int fir, string modelform, double lr, List optim_para, NumericVector maxiters, int guesses, int halfmax, double epsilon, double step_max, double thres_step_max, double deriv_epsilon, NumericMatrix df_groups, NumericVector tu, int verbose, IntegerVector KeepConstant, int term_tot, string ties_method, int nthreads, NumericVector& Strata_vals, const VectorXd& cens_weight, List model_bool, const double gmix_theta, const IntegerVector gmix_term, const MatrixXd Lin_Sys, const VectorXd Lin_Res);

List LogLik_Pois_Omnibus(const MatrixXd& PyrC, IntegerVector term_n, StringVector tform, NumericMatrix& a_ns, NumericMatrix& x_all, IntegerVector dfc, int fir, string modelform, double lr, List optim_para, NumericVector maxiters, int guesses, int halfmax, double epsilon, double step_max, double thres_step_max, double deriv_epsilon, int verbose, IntegerVector KeepConstant, int term_tot, int nthreads, const MatrixXd& dfs, List model_bool, const double gmix_theta, const IntegerVector gmix_term, const MatrixXd Lin_Sys, const VectorXd Lin_Res);

List LogLik_CaseCon_Omnibus(IntegerVector term_n, StringVector tform, NumericMatrix& a_ns, NumericMatrix& x_all, IntegerVector dfc, int fir, string modelform, double lr, List optim_para, NumericVector maxiters, int guesses, int halfmax, double epsilon, double step_max, double thres_step_max, double deriv_epsilon, NumericMatrix df_groups, NumericVector tu, int verbose, IntegerVector KeepConstant, int term_tot, string ties_method, int nthreads, NumericVector& Strata_vals, List model_bool, const double gmix_theta, const IntegerVector gmix_term, const MatrixXd Lin_Sys, const VectorXd Lin_Res);

List LogLik_Logist_Omnibus(const MatrixXd& CountEvent, IntegerVector term_n, StringVector tform, NumericMatrix& a_ns, NumericMatrix& x_all, IntegerVector dfc, int fir, string modelform, double lr, List optim_para, NumericVector maxiters, int guesses, int halfmax, double epsilon, double step_max, double thres_step_max, double deriv_epsilon, int verbose, IntegerVector KeepConstant, int term_tot, int nthreads, List model_bool, const double gmix_theta, const IntegerVector gmix_term, const MatrixXd Lin_Sys, const VectorXd Lin_Res);

#endif  //  SRC_MAIN_FUNCTIONS_H_
