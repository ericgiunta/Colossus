//  Copyright 2022 - 2025, Eric Giunta and the project collaborators, Please see main R package for license and usage details

#ifndef SRC_MAIN_MULTI_H_
#define SRC_MAIN_MULTI_H_

#include <string>

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

List LogLik_Cox_PH_Multidose_Omnibus_Serial(IntegerVector term_n, StringVector tform, Ref<VectorXd> beta_0, Ref<MatrixXd> df0, const Ref<const MatrixXd>& df1, IntegerMatrix dose_cols, IntegerVector dose_index, IntegerVector dfc, int fir, string modelform, double lr, List optim_para, int maxiter, int halfmax, double epsilon, double step_max, double thres_step_max, double deriv_epsilon, const Ref<const MatrixXd>& df_m, NumericVector tu, int verbose, IntegerVector KeepConstant, int term_tot, string ties_method, int nthreads, NumericVector& Strata_vals, const VectorXd& cens_weight, List model_bool, const double gmix_theta, const IntegerVector gmix_term, const Ref<const MatrixXd>& Lin_Sys, const Ref<const VectorXd>& Lin_Res);

List LogLik_Cox_PH_Multidose_Omnibus_Integrated(IntegerVector term_n, StringVector tform, Ref<VectorXd> beta_0, Ref<MatrixXd> df0, const Ref<const MatrixXd>& df1, IntegerMatrix dose_cols, IntegerVector dose_index, IntegerVector dfc, int fir, string modelform, double lr, List optim_para, int maxiter, int halfmax, double epsilon, double step_max, double thres_step_max, double deriv_epsilon, const Ref<const MatrixXd>& df_m, NumericVector tu, int verbose, IntegerVector KeepConstant, int term_tot, string ties_method, int nthreads, NumericVector& Strata_vals, const VectorXd& cens_weight, List model_bool, const double gmix_theta, const IntegerVector gmix_term, const Ref<const MatrixXd>& Lin_Sys, const Ref<const VectorXd>& Lin_Res);

List LogLik_Pois_PH_Multidose_Omnibus_Serial(const Ref<const MatrixXd>& PyrC, IntegerVector term_n, StringVector tform, Ref<VectorXd> beta_0, Ref<MatrixXd> df0, const Ref<const MatrixXd>& df1, IntegerMatrix dose_cols, IntegerVector dose_index, IntegerVector dfc, int fir, string modelform, double lr, List optim_para, int maxiter, int halfmax, double epsilon, double step_max, double thres_step_max, double deriv_epsilon, const Ref<const MatrixXd>& dfs, int verbose, IntegerVector KeepConstant, int term_tot, int nthreads, List model_bool, const double gmix_theta, const IntegerVector gmix_term, const Ref<const MatrixXd>& Lin_Sys, const Ref<const VectorXd>& Lin_Res);

List LogLik_Pois_PH_Multidose_Omnibus_Integrated(const Ref<const MatrixXd>& PyrC, IntegerVector term_n, StringVector tform, Ref<VectorXd> beta_0, Ref<MatrixXd> df0, const Ref<const MatrixXd>& df1, IntegerMatrix dose_cols, IntegerVector dose_index, IntegerVector dfc, int fir, string modelform, double lr, List optim_para, int maxiter, int halfmax, double epsilon, double step_max, double thres_step_max, double deriv_epsilon, const Ref<const MatrixXd>& dfs, int verbose, IntegerVector KeepConstant, int term_tot, int nthreads, List model_bool, const double gmix_theta, const IntegerVector gmix_term, const Ref<const MatrixXd>& Lin_Sys, const Ref<const VectorXd>& Lin_Res);

#endif  //  SRC_MAIN_MULTI_H_
