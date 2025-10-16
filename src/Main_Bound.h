//  Copyright 2022 - 2025, Eric Giunta and the project collaborators, Please see main R package for license and usage details

#ifndef SRC_MAIN_BOUND_H_
#define SRC_MAIN_BOUND_H_

#include <string>

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

List LogLik_Cox_PH_Omnibus_Log_Bound(IntegerVector term_n, StringVector tform, Ref<VectorXd> beta_0, Ref<MatrixXd> df0, IntegerVector dfc, int fir, string modelform, double lr, NumericVector maxiters, int guesses, int halfmax, double epsilon, double step_max, double thres_step_max, double deriv_epsilon, const Ref<const MatrixXd>& df_m, NumericVector tu, int verbose, IntegerVector KeepConstant, int term_tot, string ties_method, int nthreads, NumericVector& Strata_vals, const VectorXd& cens_weight, List model_bool, const double gmix_theta, const IntegerVector gmix_term, const Ref<const MatrixXd>& Lin_Sys, const Ref<const VectorXd>& Lin_Res, double qchi, int para_number, int maxstep, double mult);

List LogLik_Cox_PH_Omnibus_Log_Bound_Search(IntegerVector term_n, StringVector tform, Ref<VectorXd> beta_0, Ref<MatrixXd> df0, IntegerVector dfc, int fir, string modelform, double lr, NumericVector maxiters, int guesses, int halfmax, double epsilon, double step_max, double thres_step_max, double deriv_epsilon, const Ref<const MatrixXd>& df_m, NumericVector tu, int verbose, IntegerVector KeepConstant, int term_tot, string ties_method, int nthreads, NumericVector& Strata_vals, const VectorXd& cens_weight, List model_bool, const double gmix_theta, const IntegerVector gmix_term, const Ref<const MatrixXd>& Lin_Sys, const Ref<const VectorXd>& Lin_Res, double qchi, int para_number, int maxstep, double mult);

List LogLik_Poisson_Omnibus_Log_Bound(const Ref<const MatrixXd>& PyrC, const Ref<const MatrixXd>& dfs, IntegerVector term_n, StringVector tform, Ref<VectorXd> beta_0, Ref<MatrixXd> df0, IntegerVector dfc, int fir, string modelform, double lr, NumericVector maxiters, int guesses, int halfmax, double epsilon, double step_max, double thres_step_max, double deriv_epsilon, int verbose, IntegerVector KeepConstant, int term_tot, int nthreads, List model_bool, const double gmix_theta, const IntegerVector gmix_term, const Ref<const MatrixXd>& Lin_Sys, const Ref<const VectorXd>& Lin_Res, double qchi, int para_number, int maxstep, double mult);

List LogLik_Poisson_Omnibus_Log_Bound_Search(const Ref<const MatrixXd>& PyrC, const Ref<const MatrixXd>& dfs, IntegerVector term_n, StringVector tform, Ref<VectorXd> beta_0, Ref<MatrixXd> df0, IntegerVector dfc, int fir, string modelform, double lr, NumericVector maxiters, int guesses, int halfmax, double epsilon, double step_max, double thres_step_max, double deriv_epsilon, int verbose, IntegerVector KeepConstant, int term_tot, int nthreads, List model_bool, const double gmix_theta, const IntegerVector gmix_term, const Ref<const MatrixXd>& Lin_Sys, const Ref<const VectorXd>& Lin_Res, double qchi, int para_number, int maxstep, double mult);

List LogLik_Cox_PH_Omnibus_Log_Bound_CurveSearch(IntegerVector term_n, StringVector tform, Ref<VectorXd> beta_0, Ref<MatrixXd> df0, IntegerVector dfc, int fir, string modelform, double lr, List optim_para, int maxiter, int halfmax, double epsilon, double step_max, double thres_step_max, double deriv_epsilon, const Ref<const MatrixXd>& df_m, NumericVector tu, int verbose, IntegerVector KeepConstant, int term_tot, string ties_method, int nthreads, NumericVector& Strata_vals, const VectorXd& cens_weight, List model_bool, const double gmix_theta, const IntegerVector gmix_term, const Ref<const MatrixXd>& Lin_Sys, const Ref<const VectorXd>& Lin_Res, double qchi, int para_number, int maxstep, double step_size);

List LogLik_Poisson_Omnibus_Log_Bound_CurveSearch(const Ref<const MatrixXd>& PyrC, const Ref<const MatrixXd>& dfs, IntegerVector term_n, StringVector tform, Ref<VectorXd> beta_0, Ref<MatrixXd> df0, IntegerVector dfc, int fir, string modelform, double lr, List optim_para, int maxiter, int halfmax, double epsilon, double step_max, double thres_step_max, double deriv_epsilon, int verbose, IntegerVector KeepConstant, int term_tot, int nthreads, List model_bool, const double gmix_theta, const IntegerVector gmix_term, const Ref<const MatrixXd>& Lin_Sys, const Ref<const VectorXd>& Lin_Res, double qchi, int para_number, int maxstep, double step_size);

#endif  //  SRC_MAIN_BOUND_H_
