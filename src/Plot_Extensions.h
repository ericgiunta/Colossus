//  Copyright 2022 - 2025, Eric Giunta and the project collaborators, Please see main R package for license and usage details

#ifndef SRC_PLOT_EXTENSIONS_H_
#define SRC_PLOT_EXTENSIONS_H_

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

List PLOT_SURV(int reqrdnum, MatrixXd& R, MatrixXd& Rd, NumericVector& a_er, const Ref<const MatrixXd>& df_m, NumericVector& tu, int verbose, int nthreads);

List PLOT_SURV_Strata(int reqrdnum, MatrixXd& R, MatrixXd& Rd, NumericVector& a_er, const Ref<const MatrixXd>& df_m, NumericVector& tu, NumericVector& Strata_vals, int verbose, int nthreads);

List Schoenfeld_Calc(int ntime, int totalnum, const  VectorXd& beta_0, const Ref<const MatrixXd>& df0, const MatrixXd& R, MatrixXd& Lldd_inv, const IntegerMatrix& RiskFail, const vector<vector<int> >& RiskPairs, IntegerVector& dfc, int verbose, IntegerVector KeepConstant, int nthreads);

List Plot_Omnibus(IntegerVector term_n, StringVector tform, Ref<VectorXd> beta_0, const Ref<const MatrixXd>& df0, IntegerVector dfc, int fir, int der_iden, string modelform, double step_max, double thres_step_max, const Ref<const MatrixXd>& df_m, NumericVector& tu, int verbose, IntegerVector KeepConstant, int term_tot, string ties_method, int nthreads, NumericVector& Strata_vals, const VectorXd& cens_weight, int uniq_v, List model_bool, bool Surv_bool, bool Risk_bool, bool Schoenfeld_bool, bool Risk_Sub_bool, const double gmix_theta, const IntegerVector& gmix_term);

List Assign_Events_Pois(IntegerVector term_n, StringVector tform, Ref<VectorXd> beta_0, Ref<MatrixXd> df0, IntegerVector dfc, const Ref<const MatrixXd>& PyrC, const Ref<const MatrixXd>& dfs, int fir, string modelform, int verbose, IntegerVector KeepConstant, int term_tot, int nthreads, const double gmix_theta, const IntegerVector gmix_term, List model_bool);

List Poisson_Residuals(const Ref<const MatrixXd>& PyrC, IntegerVector term_n, StringVector tform, Ref<VectorXd> beta_0, Ref<MatrixXd> df0, IntegerVector dfc, int fir, string modelform, double step_max, double thres_step_max, int verbose, IntegerVector KeepConstant, int term_tot, int nthreads, const Ref<const MatrixXd>& dfs, List model_bool, const double gmix_theta, const IntegerVector gmix_term, bool Pearson_bool, bool Deviance_bool);

#endif  //  SRC_PLOT_EXTENSIONS_H_
