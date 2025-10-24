//  Copyright 2022 - 2025, Eric Giunta and the project collaborators, Please see main R package for license and usage details

#ifndef SRC_R_INTERFACE_H_
#define SRC_R_INTERFACE_H_

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

List cox_ph_Omnibus_transition(IntegerVector term_n, StringVector tform, NumericMatrix& a_ns, IntegerVector dfc, Ref<MatrixXd> df0, int fir, string modelform, List Control, const MatrixXd df_m, NumericVector tu, IntegerVector KeepConstant, int term_tot, NumericVector Strata_vals, VectorXd cens_weight, List model_control, MatrixXd Lin_Sys, VectorXd Lin_Res);

List pois_Omnibus_transition(MatrixXd PyrC, IntegerVector term_n, StringVector tform, NumericMatrix& a_ns, IntegerVector dfc, Ref<MatrixXd> df0, int fir, string modelform, List Control, IntegerVector KeepConstant, int term_tot, MatrixXd dfs, List model_control, MatrixXd Lin_Sys, VectorXd Lin_Res);

List Plot_Omnibus_transition(IntegerVector term_n, StringVector tform, NumericVector a_n, IntegerVector dfc, Ref<MatrixXd> df0, int fir, int der_iden, string modelform, List Control, const MatrixXd df_m, NumericVector tu, IntegerVector KeepConstant, int term_tot, NumericVector Strata_vals, VectorXd cens_weight, List model_control);

List Assigned_Event_Poisson_transition(MatrixXd PyrC, MatrixXd dfs, IntegerVector term_n, StringVector tform, NumericVector a_n, IntegerVector dfc, Ref<MatrixXd> df0, int fir, string modelform, List Control, IntegerVector KeepConstant, int term_tot, List model_control);

List cox_ph_Omnibus_Bounds_transition(IntegerVector term_n, StringVector tform, NumericVector a_n, IntegerVector dfc, Ref<MatrixXd> df0, int fir, string modelform, List Control, const MatrixXd df_m, NumericVector tu, IntegerVector KeepConstant, int term_tot, NumericVector Strata_vals, VectorXd cens_weight, List model_control, MatrixXd Lin_Sys, VectorXd Lin_Res);

List cox_ph_Omnibus_CurveSearch_transition(IntegerVector term_n, StringVector tform, NumericVector a_n, IntegerVector dfc, Ref<MatrixXd> df0, int fir, string modelform, List Control, const MatrixXd df_m, NumericVector tu, IntegerVector KeepConstant, int term_tot, NumericVector Strata_vals, VectorXd cens_weight, List model_control, MatrixXd Lin_Sys, VectorXd Lin_Res);

List pois_Omnibus_Bounds_transition(MatrixXd PyrC, IntegerVector term_n, StringVector tform, NumericVector a_n, IntegerVector dfc, Ref<MatrixXd> df0, int fir, string modelform, List Control, IntegerVector KeepConstant, int term_tot, MatrixXd dfs, List model_control, MatrixXd Lin_Sys, VectorXd Lin_Res);

List pois_Omnibus_CurveSearch_transition(MatrixXd PyrC, IntegerVector term_n, StringVector tform, NumericVector a_n, IntegerVector dfc, Ref<MatrixXd> df0, int fir, string modelform, List Control, IntegerVector KeepConstant, int term_tot, MatrixXd dfs, List model_control, MatrixXd Lin_Sys, VectorXd Lin_Res);

List pois_Residual_transition(MatrixXd PyrC, IntegerVector term_n, StringVector tform, NumericVector a_n, IntegerVector dfc, Ref<MatrixXd> df0, int fir, string modelform, List Control, IntegerVector KeepConstant, int term_tot, MatrixXd dfs, List model_control);

List cox_ph_multidose_Omnibus_transition(IntegerVector term_n, StringVector tform, NumericVector a_n, IntegerMatrix dose_cols, IntegerVector dose_index, IntegerVector dfc, Ref<MatrixXd> df0, MatrixXd df1, int fir, string modelform, List Control, const MatrixXd df_m, NumericVector tu, IntegerVector KeepConstant, int term_tot, NumericVector Strata_vals, VectorXd cens_weight, List model_control, MatrixXd Lin_Sys, VectorXd Lin_Res);

List pois_multidose_Omnibus_transition(MatrixXd PyrC, IntegerVector term_n, StringVector tform, NumericVector a_n, IntegerMatrix dose_cols, IntegerVector dose_index, IntegerVector dfc, Ref<MatrixXd> df0, MatrixXd df1, int fir, string modelform, List Control, IntegerVector KeepConstant, int term_tot, MatrixXd dfs, List model_control, MatrixXd Lin_Sys, VectorXd Lin_Res);

List caco_Omnibus_transition(IntegerVector term_n, StringVector tform, NumericMatrix& a_ns, IntegerVector dfc, Ref<MatrixXd> df0, int fir, string modelform, List Control, const MatrixXd df_m, NumericVector tu, IntegerVector KeepConstant, int term_tot, NumericVector Strata_vals, VectorXd cens_weight, List model_control, MatrixXd Lin_Sys, VectorXd Lin_Res);

List logist_Omnibus_transition(MatrixXd CountEvent, IntegerVector term_n, StringVector tform, NumericMatrix& a_ns, IntegerVector dfc, Ref<MatrixXd> df0, int fir, string modelform, List Control, IntegerVector KeepConstant, int term_tot, List model_control, MatrixXd Lin_Sys, VectorXd Lin_Res);

void Write_Time_Dep(const MatrixXd dfs_Times, const MatrixXd dfs_dep, const MatrixXd dfs_const, const NumericVector df0_event, double dt, string filename, StringVector tform_tdep, NumericVector tu, bool iscox, int nthreads);

NumericMatrix Gen_Fac_Par(const MatrixXd dfs, const NumericVector vals, const NumericVector cols, const int nthreads);

bool risk_check_transition(IntegerVector term_n, StringVector tform, NumericVector a_n, IntegerVector dfc, Ref<MatrixXd> df0, int fir, string modelform, List Control, List model_control, IntegerVector KeepConstant, int term_tot);

void Gen_Strat_Weight(string modelform, const Ref<const MatrixXd>& dfs, const Ref<const MatrixXd>& PyrC, VectorXd& s_weights, const int nthreads, const StringVector& tform, const IntegerVector& term_n, const int& term_tot);

bool OMP_Check();

#endif  //  SRC_R_INTERFACE_H_
