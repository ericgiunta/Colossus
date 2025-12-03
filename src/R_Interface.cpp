//  Copyright 2022 - 2025, Eric Giunta and the project collaborators, Please see main R package for license and usage details

#include <RcppEigen.h>

#include "R_Interface.h"
#ifdef _OPENMP
#include <omp.h>
#endif

#include <Eigen/Core>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <random>
#include <ctime>

#include "Main_Functions.h"
#include "Main_Bound.h"
#include "Main_Multi.h"
#include "Plot_Extensions.h"


//  [[Rcpp::depends(RcppEigen)]]
//  [[Rcpp::plugins(openmp)]]

using std::string;
using std::vector;
using std::ofstream;
using std::endl;
using std::invalid_argument;

using Eigen::Map;
using Eigen::Ref;
using Eigen::ArrayXd;
using Eigen::MatrixXd;
using Eigen::SparseMatrix;
using Eigen::VectorXd;
using Eigen::VectorXi;
using Eigen::IOFormat;
using Eigen::FullPrecision;
using Eigen::DontAlignCols;

using Rcpp::as;
using Rcpp::wrap;
using Rcpp::IntegerMatrix;
using Rcpp::IntegerVector;
using Rcpp::NumericVector;
using Rcpp::NumericMatrix;
using Rcpp::StringVector;
using Rcpp::List;
using Rcpp::_;
using Rcpp::Rcout;

template <typename T> int sign(T val) {
    return (T(0) < val) - (val < T(0));
}

template<typename Func>
struct lambda_as_visitor_wrapper : Func {
    lambda_as_visitor_wrapper(const Func& f) : Func(f) {}
    template<typename S, typename I>
    void init(const S& v, I i, I j) { return Func::operator()(v, i, j); }
};

template<typename Mat, typename Func>
void visit_lambda(const Mat& m, const Func& f) {
    lambda_as_visitor_wrapper<Func> visitor(f);
    m.visit(visitor);
}

//' Interface between R code and the Cox PH omnibus regression function
//'
//' \code{cox_ph_Omnibus_transition} Called directly from R, Defines the control variables and calls the regression function
//' @inheritParams CPP_template
//'
//' @return LogLik_Cox_PH output : Log-likelihood of optimum, first derivative of log-likelihood, second derivative matrix, parameter list, standard deviation estimate, AIC, model information
//' @noRd
//'
//  [[Rcpp::export]]
List cox_ph_Omnibus_transition(IntegerVector term_n, StringVector tform, NumericMatrix& a_ns, IntegerVector dfc, MatrixXd df0, int fir, string modelform, List Control, const MatrixXd df_m, NumericVector tu, IntegerVector KeepConstant, int term_tot, NumericVector Strata_vals, VectorXd cens_weight, List model_control, MatrixXd Lin_Sys, VectorXd Lin_Res) {
    int verbose = Control["verbose"];
    //
    double lr = Control["lr"];
    NumericVector maxiters = Control["maxiters"];
    int guesses = Control["guesses"];
    int halfmax = Control["halfmax"];
    double epsilon = Control["epsilon"];
    //
    double step_max = Control["step_max"];
    double thres_step_max = Control["thres_step_max"];
    double deriv_epsilon = Control["deriv_epsilon"];
    string ties_method = Control["ties"];
    int nthreads = Control["ncores"];
    //
    double gmix_theta = model_control["gmix_theta"];
    IntegerVector gmix_term = model_control["gmix_term"];
    //
    List model_bool = List::create(
            _["strata"] = model_control["strata"],
            _["basic"] = model_control["basic"],
            _["linear_err"] = model_control["linear_err"],
            _["null"] = model_control["null"],
            _["cr"] = model_control["cr"],
            _["single"] = model_control["single"],
            _["gradient"] = model_control["gradient"],
            _["outcome_prob"] = false,
            _["constraint"] = model_control["constraint"],
            _["observed_info"] = model_control["observed_info"],
            _["log_bound"] = false,
            _["cox"] = true);
    List optim_para = List::create(
            _["lr"] = Control["lr"]);
    if (model_bool["gradient"]) {
        model_bool["momentum"] = model_control["momentum"];
        model_bool["adadelta"] = model_control["adadelta"];
        model_bool["adam"] = model_control["adam"];
        optim_para["momentum_decay"] = model_control["momentum_decay"];
        optim_para["learning_decay"] = model_control["learning_decay"];
        optim_para["epsilon_decay"] = model_control["epsilon_decay"];
        if (model_bool["constraint"]) {
            optim_para["penalty_weight"] = model_control["penalty_weight"];
            optim_para["penalty_method"] = model_control["penalty_method"];
        }
    }
    //
    //  Performs regression
    //----------------------------------------------------------------------------------------------------------------//
    List res = LogLik_Cox_PH_Omnibus(term_n, tform, a_ns, df0, dfc, fir, modelform, lr, optim_para, maxiters, guesses, halfmax, epsilon, step_max, thres_step_max, deriv_epsilon, df_m, tu, verbose, KeepConstant, term_tot, ties_method, nthreads, Strata_vals, cens_weight, model_bool, gmix_theta, gmix_term, Lin_Sys, Lin_Res);
    //----------------------------------------------------------------------------------------------------------------//
    return res;
}

//' Interface between R code and the poisson omnibus regression function
//'
//' \code{pois_Omnibus_transition} Called directly from R, Defines the control variables and calls the regression function
//' @inheritParams CPP_template
//'
//' @return LogLik_Cox_PH output : Log-likelihood of optimum, first derivative of log-likelihood, second derivative matrix, parameter list, standard deviation estimate, AIC, model information
//' @noRd
//'
//  [[Rcpp::export]]
List pois_Omnibus_transition(MatrixXd PyrC, IntegerVector term_n, StringVector tform, NumericMatrix& a_ns, IntegerVector dfc, MatrixXd df0, int fir, string modelform, List Control, IntegerVector KeepConstant, int term_tot, MatrixXd dfs, List model_control, MatrixXd Lin_Sys, VectorXd Lin_Res) {
    //
    int verbose = Control["verbose"];
    //
    double lr = Control["lr"];
    NumericVector maxiters = Control["maxiters"];
    int guesses = Control["guesses"];
    int halfmax = Control["halfmax"];
    double epsilon = Control["epsilon"];
    //
    double step_max = Control["step_max"];
    double thres_step_max = Control["thres_step_max"];
    double deriv_epsilon = Control["deriv_epsilon"];
    int nthreads = Control["ncores"];
    //
    double gmix_theta = model_control["gmix_theta"];
    IntegerVector gmix_term = model_control["gmix_term"];
    //
    List model_bool = List::create(
            _["strata"] = model_control["strata"],
            _["basic"] = false,
            _["linear_err"] = false,
            _["null"] = false,
            _["cr"] = false,
            _["single"] = model_control["single"],
            _["gradient"] = model_control["gradient"],
            _["outcome_prob"] = false,
            _["constraint"] = model_control["constraint"],
            _["observed_info"] = model_control["observed_info"],
            _["log_bound"] = false,
            _["cox"] = false);
    List optim_para = List::create(
            _["lr"] = Control["lr"]);
    if (model_bool["gradient"]) {
        model_bool["momentum"] = model_control["momentum"];
        model_bool["adadelta"] = model_control["adadelta"];
        model_bool["adam"] = model_control["adam"];
        optim_para["momentum_decay"] = model_control["momentum_decay"];
        optim_para["learning_decay"] = model_control["learning_decay"];
        optim_para["epsilon_decay"] = model_control["epsilon_decay"];
        if (model_bool["constraint"]) {
            optim_para["penalty_weight"] = model_control["penalty_weight"];
            optim_para["penalty_method"] = model_control["penalty_method"];
        }
    }
    //
    //  Performs regression
    //----------------------------------------------------------------------------------------------------------------//
    List res = LogLik_Pois_Omnibus(PyrC, term_n, tform, a_ns, df0, dfc, fir, modelform, lr, optim_para, maxiters, guesses, halfmax, epsilon, step_max, thres_step_max, deriv_epsilon, verbose, KeepConstant, term_tot, nthreads, dfs, model_bool, gmix_theta, gmix_term, Lin_Sys, Lin_Res);
    //----------------------------------------------------------------------------------------------------------------//
    return res;
}

//' Interface between R code and the event assignment omnibus function
//'
//' \code{Assigned_Event_transition} Called directly from R, Defines the control variables and calls the assigning functions
//' @inheritParams CPP_template
//'
//' @return list of assigned/predicted background/excess events
//' @noRd
//'
//  [[Rcpp::export]]
List Assigned_Event_Poisson_transition(MatrixXd PyrC, MatrixXd dfs, IntegerVector term_n, StringVector tform, NumericVector a_n, IntegerVector dfc, MatrixXd df0, int fir, string modelform, List Control, IntegerVector KeepConstant, int term_tot, List model_control) {
    int verbose = Control["verbose"];
    //
    string ties_method = Control["ties"];
    int nthreads = Control["ncores"];
    //
    double gmix_theta = model_control["gmix_theta"];
    IntegerVector gmix_term = model_control["gmix_term"];
    //
    Map<VectorXd> beta_0(as<Map<VectorXd> >(a_n));
    //
    List model_bool = List::create(
            _["strata"] = model_control["strata"],
            _["basic"] = false,
            _["linear_err"] = false,
            _["null"] = false,
            _["cr"] = false,
            _["single"] = false,
            _["gradient"] = false,
            _["outcome_prob"] = false,
            _["constraint"] = false,
            _["log_bound"] = false,
            _["cox"] = false);
    //  Performs regression
    List res;
    //----------------------------------------------------------------------------------------------------------------//
    res = Assign_Events_Pois(term_n, tform, beta_0, df0, dfc, PyrC, dfs, fir, modelform, verbose, KeepConstant, term_tot, nthreads, gmix_theta, gmix_term, model_bool);
    //----------------------------------------------------------------------------------------------------------------//
    return res;
}

//' Interface between R code and the plotting omnibus function
//'
//' \code{Plot_Omnibus_transition} Called directly from R, Defines the control variables and calls the plotting functions
//' @inheritParams CPP_template
//'
//' @return LogLik_Cox_PH output : Log-likelihood of optimum, first derivative of log-likelihood, second derivative matrix, parameter list, standard deviation estimate, AIC, model information
//' @noRd
//'
//  [[Rcpp::export]]
List Plot_Omnibus_transition(IntegerVector term_n, StringVector tform, NumericVector a_n, IntegerVector dfc, MatrixXd df0, int fir, int der_iden, string modelform, List Control, const MatrixXd df_m, NumericVector tu, IntegerVector KeepConstant, int term_tot, NumericVector Strata_vals, VectorXd cens_weight, List model_control) {
    int verbose = Control["verbose"];
    //
    double step_max = Control["step_max"];
    double thres_step_max = Control["thres_step_max"];
    string ties_method = Control["ties"];
    int nthreads = Control["ncores"];
    //
    double gmix_theta = model_control["gmix_theta"];
    IntegerVector gmix_term = model_control["gmix_term"];
    //
    List model_bool = List::create(
            _["strata"] = model_control["strata"],
            _["basic"] = model_control["basic"],
            _["linear_err"] = model_control["linear_err"],
            _["observed_info"] = model_control["observed_info"],
            _["null"] = false,
            _["cr"] = model_control["cr"],
            _["single"] = false,
            _["gradient"] = false,
            _["outcome_prob"] = false,
            _["constraint"] = false,
            _["log_bound"] = false,
            _["cox"] = true);
    //
    bool Surv_bool       = model_control["surv"];
    bool Schoenfeld_bool = model_control["schoenfeld"];
    bool Risk_bool       = model_control["risk"];
    bool Risk_Sub_bool   = model_control["risk_subset"];
    int uniq_v           = model_control["unique_values"];
    //
    //  Performs regression
    List res;
    if (uniq_v < 2) {
        res = List::create(_["Failed"] = "Unique_Values too low, expects atleast 2 values");
        return res;
    }
    //
    Map<VectorXd> beta_0(as<Map<VectorXd> >(a_n));
    int ijk_risk = 0;
    vector<float> vv;  //  stores the covariate values
    if (Risk_bool) {
        float dx = 0;
        if (der_iden >=0) {
        } else {
            throw invalid_argument("Incorrect parameter to plot by");
        }
        if (uniq_v > 100) {  //  selects anything above 100 points to be continuous
            vv.resize(100);  //  continuous covariates use 100 steps
        } else {
            vv.resize(uniq_v);  //  factor covariates use the number of factors
        }
        MatrixXd df1 = MatrixXd::Zero(vv.size(), df0.cols());  //  stores memory for the derivative term parameters and columns
        df1 = df1.array();
        ijk_risk = dfc[der_iden] - 1;
        dx = (df0.col(ijk_risk).maxCoeff() - df0.col(ijk_risk).minCoeff())/(vv.size() - 1);  //  varies from max to minimum
        vv[0] = df0.col(ijk_risk).minCoeff();
        generate(vv.begin(), vv.end(), [n = 0, &dx]() mutable { return n++ * dx; });
        #ifdef _OPENMP
        #pragma omp parallel for schedule(dynamic) num_threads(nthreads)
        #endif
        for (vector<float>::size_type ij = 0; ij < vv.size(); ij++) {
            df1(ij, ijk_risk) = vv[ij];  //  fills the column with varying values
        }
        res = Plot_Omnibus(term_n, tform, beta_0, df1, dfc, fir, der_iden, modelform, step_max, thres_step_max, df_m, tu, verbose, KeepConstant, term_tot, ties_method, nthreads, Strata_vals, cens_weight, uniq_v, model_bool, Surv_bool, Risk_bool, Schoenfeld_bool, Risk_Sub_bool, gmix_theta, gmix_term);
    } else {
        res = Plot_Omnibus(term_n, tform, beta_0, df0, dfc, fir, der_iden, modelform, step_max, thres_step_max, df_m, tu, verbose, KeepConstant, term_tot, ties_method, nthreads, Strata_vals, cens_weight, uniq_v, model_bool, Surv_bool, Risk_bool, Schoenfeld_bool, Risk_Sub_bool, gmix_theta, gmix_term);
    }
    //----------------------------------------------------------------------------------------------------------------//
    return res;
}

//' Interface between R code and the Cox PH omnibus bounds regression function
//'
//' \code{cox_ph_Omnibus_Bounds_transition} Called directly from R, Defines the control variables and calls the regression function
//' @inheritParams CPP_template
//'
//' @return LogLik_Cox_PH output : Log-likelihood of optimum, first derivative of log-likelihood, second derivative matrix, parameter list, standard deviation estimate, AIC, model information
//' @noRd
//'
//  [[Rcpp::export]]
List cox_ph_Omnibus_Bounds_transition(IntegerVector term_n, StringVector tform, NumericVector a_n, IntegerVector dfc, MatrixXd df0, int fir, string modelform, List Control, const MatrixXd df_m, NumericVector tu, IntegerVector KeepConstant, int term_tot, NumericVector Strata_vals, VectorXd cens_weight, List model_control, MatrixXd Lin_Sys, VectorXd Lin_Res) {
    int verbose = Control["verbose"];
    //
    double lr = Control["lr"];
    NumericVector maxiters = Control["maxiters"];
    int guesses = Control["guesses"];
    int halfmax = Control["halfmax"];
    double epsilon = Control["epsilon"];
    //
    double step_max = Control["step_max"];
    double thres_step_max = Control["thres_step_max"];
    double deriv_epsilon = Control["deriv_epsilon"];
    string ties_method = Control["ties"];
    int nthreads = Control["ncores"];
    //
    Map<VectorXd> beta_0(as<Map<VectorXd> >(a_n));
    //
    double gmix_theta = model_control["gmix_theta"];
    IntegerVector gmix_term = model_control["gmix_term"];
    double mult = model_control["search_mult"];
    //
    List model_bool = List::create(
            _["strata"] = model_control["strata"],
            _["basic"] = model_control["basic"],
            _["linear_err"] = model_control["linear_err"],
            _["null"] = model_control["null"],
            _["cr"] = model_control["cr"],
            _["single"] = model_control["single"],
            _["gradient"] = model_control["gradient"],
            _["outcome_prob"] = false,
            _["constraint"] = model_control["constraint"],
            _["observed_info"] = model_control["observed_info"],
            _["log_bound"] = true,
            _["cox"] = true);
    List optim_para = List::create(
            _["lr"] = Control["lr"]);
    //
    double qchi    = model_control["qchi"];
    int para_number = model_control["para_number"];
    para_number -= 1;
    //
    int maxstep    = model_control["maxstep"];
    //
    bool manual   = model_control["manual"];
    //  Performs regression
    //----------------------------------------------------------------------------------------------------------------//
    List res;
    if (manual) {
        res = LogLik_Cox_PH_Omnibus_Log_Bound_Search(term_n, tform, beta_0, df0, dfc, fir, modelform, lr, maxiters, guesses, halfmax, epsilon, step_max, thres_step_max, deriv_epsilon, df_m, tu, verbose, KeepConstant, term_tot, ties_method, nthreads, Strata_vals, cens_weight, model_bool, gmix_theta, gmix_term, Lin_Sys, Lin_Res, qchi, para_number, maxstep, mult);
    } else {
        res = LogLik_Cox_PH_Omnibus_Log_Bound(term_n, tform, beta_0, df0, dfc, fir, modelform, lr, maxiters, guesses, halfmax, epsilon, step_max, thres_step_max, deriv_epsilon, df_m, tu, verbose, KeepConstant, term_tot, ties_method, nthreads, Strata_vals, cens_weight, model_bool, gmix_theta, gmix_term, Lin_Sys, Lin_Res, qchi, para_number, maxstep, mult);
    }
    //----------------------------------------------------------------------------------------------------------------//
    return res;
}
//' Interface between R code and the Cox PH omnibus bounds regression function
//'
//' \code{cox_ph_Omnibus_CurveSearch_transition} Called directly from R, Defines the control variables and calls the regression function
//' @inheritParams CPP_template
//'
//' @return LogLik_Cox_PH output : Log-likelihood of optimum, first derivative of log-likelihood, second derivative matrix, parameter list, standard deviation estimate, AIC, model information
//' @noRd
//'
//  [[Rcpp::export]]
List cox_ph_Omnibus_CurveSearch_transition(IntegerVector term_n, StringVector tform, NumericVector a_n, IntegerVector dfc, MatrixXd df0, int fir, string modelform, List Control, const MatrixXd df_m, NumericVector tu, IntegerVector KeepConstant, int term_tot, NumericVector Strata_vals, VectorXd cens_weight, List model_control, MatrixXd Lin_Sys, VectorXd Lin_Res) {
    int verbose = Control["verbose"];
    //
    double lr = Control["lr"];
    int maxiter = Control["maxiter"];
    int halfmax = Control["halfmax"];
    double epsilon = Control["epsilon"];
    //
    double step_max = Control["step_max"];
    double thres_step_max = Control["thres_step_max"];
    double deriv_epsilon = Control["deriv_epsilon"];
    string ties_method = Control["ties"];
    int nthreads = Control["ncores"];
    //
    Map<VectorXd> beta_0(as<Map<VectorXd> >(a_n));
    //
    double gmix_theta = model_control["gmix_theta"];
    IntegerVector gmix_term = model_control["gmix_term"];
    //
    List model_bool = List::create(
            _["strata"] = model_control["strata"],
            _["basic"] = model_control["basic"],
            _["linear_err"] = model_control["linear_err"],
            _["null"] = model_control["null"],
            _["cr"] = model_control["cr"],
            _["single"] = model_control["single"],
            _["gradient"] = model_control["gradient"],
            _["outcome_prob"] = false,
            _["constraint"] = model_control["constraint"],
            _["observed_info"] = model_control["observed_info"],
            _["log_bound"] = true,
            _["cox"] = true);
    List optim_para = List::create(
            _["lr"] = Control["lr"]);
    //
    double qchi    = model_control["qchi"];
    int para_number = model_control["para_number"];
    para_number -= 1;
    //
    int maxstep    = model_control["maxstep"];
    double step_size    = model_control["step_size"];
    //
    //  Performs regression
    //----------------------------------------------------------------------------------------------------------------//
    List res;
    res = LogLik_Cox_PH_Omnibus_Log_Bound_CurveSearch(term_n, tform, beta_0, df0, dfc, fir, modelform, lr, optim_para, maxiter, halfmax, epsilon, step_max, thres_step_max, deriv_epsilon, df_m, tu, verbose, KeepConstant, term_tot, ties_method, nthreads, Strata_vals, cens_weight, model_bool, gmix_theta, gmix_term, Lin_Sys, Lin_Res, qchi, para_number, maxstep, step_size);
    //----------------------------------------------------------------------------------------------------------------//
    return res;
}

//' Interface between R code and the poisson omnibus bounds regression function
//'
//' \code{pois_Omnibus_CurveSearch_transition} Called directly from R, Defines the control variables and calls the regression function
//' @inheritParams CPP_template
//'
//' @return LogLik_Cox_PH output : Log-likelihood of optimum, first derivative of log-likelihood, second derivative matrix, parameter list, standard deviation estimate, AIC, model information
//' @noRd
//'
//  [[Rcpp::export]]
List pois_Omnibus_CurveSearch_transition(MatrixXd PyrC, IntegerVector term_n, StringVector tform, NumericVector a_n, IntegerVector dfc, MatrixXd df0, int fir, string modelform, List Control, IntegerVector KeepConstant, int term_tot, MatrixXd dfs, List model_control, MatrixXd Lin_Sys, VectorXd Lin_Res) {
    // const Map<MatrixXd> PyrC(as<Map<MatrixXd> >(dfe));
    // const Map<MatrixXd> dfs(as<Map<MatrixXd> >(df_strata));
    //
    int verbose = Control["verbose"];
    //
    double lr = Control["lr"];
    int maxiter = Control["maxiter"];
    int halfmax = Control["halfmax"];
    double epsilon = Control["epsilon"];
    //
    double step_max = Control["step_max"];
    double thres_step_max = Control["thres_step_max"];
    double deriv_epsilon = Control["deriv_epsilon"];
    int nthreads = Control["ncores"];
    //
    //
    double gmix_theta = model_control["gmix_theta"];
    IntegerVector gmix_term = model_control["gmix_term"];
    //
    Map<VectorXd> beta_0(as<Map<VectorXd> >(a_n));
    //
    List model_bool = List::create(
            _["strata"] = model_control["strata"],
            _["basic"] = false,
            _["linear_err"] = false,
            _["null"] = false,
            _["cr"] = false,
            _["single"] = model_control["single"],
            _["gradient"] = model_control["gradient"],
            _["outcome_prob"] = false,
            _["constraint"] = model_control["constraint"],
            _["observed_info"] = model_control["observed_info"],
            _["log_bound"] = false,
            _["cox"] = false);
    List optim_para = List::create(
            _["lr"] = Control["lr"]);
    //
    double qchi    = model_control["qchi"];
    int para_number = model_control["para_number"];
    para_number -= 1;
    //
    int maxstep    = model_control["maxstep"];
    double step_size = model_control["step_size"];
    //
    //  Performs regression
    //----------------------------------------------------------------------------------------------------------------//
    List res;
    res = LogLik_Poisson_Omnibus_Log_Bound_CurveSearch(PyrC, dfs, term_n, tform, beta_0, df0, dfc, fir, modelform, lr, optim_para, maxiter, halfmax, epsilon, step_max, thres_step_max, deriv_epsilon, verbose, KeepConstant, term_tot, nthreads, model_bool, gmix_theta, gmix_term, Lin_Sys, Lin_Res, qchi, para_number, maxstep, step_size);
    //----------------------------------------------------------------------------------------------------------------//
    return res;
}

//' Interface between R code and the poisson omnibus bounds regression function
//'
//' \code{pois_Omnibus_Bounds_transition} Called directly from R, Defines the control variables and calls the regression function
//' @inheritParams CPP_template
//'
//' @return LogLik_Cox_PH output : Log-likelihood of optimum, first derivative of log-likelihood, second derivative matrix, parameter list, standard deviation estimate, AIC, model information
//' @noRd
//'
//  [[Rcpp::export]]
List pois_Omnibus_Bounds_transition(MatrixXd PyrC, IntegerVector term_n, StringVector tform, NumericVector a_n, IntegerVector dfc, MatrixXd df0, int fir, string modelform, List Control, IntegerVector KeepConstant, int term_tot, MatrixXd dfs, List model_control, MatrixXd Lin_Sys, VectorXd Lin_Res) {
    //
    int verbose = Control["verbose"];
    //
    double lr = Control["lr"];
    NumericVector maxiters = Control["maxiters"];
    int guesses = Control["guesses"];
    int halfmax = Control["halfmax"];
    double epsilon = Control["epsilon"];
    //
    double step_max = Control["step_max"];
    double thres_step_max = Control["thres_step_max"];
    double deriv_epsilon = Control["deriv_epsilon"];
    int nthreads = Control["ncores"];
    //
    //
    double gmix_theta = model_control["gmix_theta"];
    IntegerVector gmix_term = model_control["gmix_term"];
    double mult = model_control["search_mult"];
    //
    Map<VectorXd> beta_0(as<Map<VectorXd> >(a_n));
    //
    List model_bool = List::create(
            _["strata"] = model_control["strata"],
            _["basic"] = false,
            _["linear_err"] = false,
            _["null"] = false,
            _["cr"] = false,
            _["single"] = model_control["single"],
            _["gradient"] = model_control["gradient"],
            _["outcome_prob"] = false,
            _["constraint"] = model_control["constraint"],
            _["observed_info"] = model_control["observed_info"],
            _["log_bound"] = false,
            _["cox"] = false);
    List optim_para = List::create(
            _["lr"] = Control["lr"]);
    //
    double qchi    = model_control["qchi"];
    int para_number = model_control["para_number"];
    para_number -= 1;
    //
    int maxstep    = model_control["maxstep"];
    //
    bool manual   = model_control["manual"];
    //  Performs regression
    //----------------------------------------------------------------------------------------------------------------//
    List res;
    if (manual) {
        res = LogLik_Poisson_Omnibus_Log_Bound_Search(PyrC, dfs, term_n, tform, beta_0, df0, dfc, fir, modelform, lr, maxiters, guesses, halfmax, epsilon, step_max, thres_step_max, deriv_epsilon, verbose, KeepConstant, term_tot, nthreads, model_bool, gmix_theta, gmix_term, Lin_Sys, Lin_Res, qchi, para_number, maxstep, mult);
    } else {
        res = LogLik_Poisson_Omnibus_Log_Bound(PyrC, dfs, term_n, tform, beta_0, df0, dfc, fir, modelform, lr, maxiters, guesses, halfmax, epsilon, step_max, thres_step_max, deriv_epsilon, verbose, KeepConstant, term_tot, nthreads, model_bool, gmix_theta, gmix_term, Lin_Sys, Lin_Res, qchi, para_number, maxstep, mult);
    }
    //----------------------------------------------------------------------------------------------------------------//
    return res;
}

//' Interface between R code and the poisson residual calculation function
//'
//' \code{pois_Residual_transition} Called directly from R, Defines the control variables and calls the calculation function
//' @inheritParams CPP_template
//'
//' @return Poisson_Residuals output : list of residuals and sum
//' @noRd
//'
//  [[Rcpp::export]]
List pois_Residual_transition(MatrixXd PyrC, IntegerVector term_n, StringVector tform, NumericVector a_n, IntegerVector dfc, MatrixXd df0, int fir, string modelform, List Control, IntegerVector KeepConstant, int term_tot, MatrixXd dfs, List model_control) {
    //
    Map<VectorXd> beta_0(as<Map<VectorXd> >(a_n));
    //
    int verbose = Control["verbose"];
    //
    double step_max = Control["step_max"];
    double thres_step_max = Control["thres_step_max"];
    int nthreads = Control["ncores"];
    //
    double gmix_theta = model_control["gmix_theta"];
    IntegerVector gmix_term = model_control["gmix_term"];
    //
    //
    bool Pearson_bool = model_control["pearson"];
    bool Deviance_bool = model_control["deviance"];
    //
    List model_bool = List::create(
            _["strata"] = model_control["strata"],
            _["basic"] = false,
            _["linear_err"] = false,
            _["null"] = false,
            _["cr"] = false,
            _["single"] = true,
            _["gradient"] = false,
            _["outcome_prob"] = false,
            _["constraint"] = false,
            _["log_bound"] = false,
            _["cox"] = false);
    //  Performs regression
    //----------------------------------------------------------------------------------------------------------------//
    List res = Poisson_Residuals(PyrC, term_n, tform, beta_0, df0, dfc, fir, modelform, step_max, thres_step_max, verbose, KeepConstant, term_tot, nthreads, dfs, model_bool, gmix_theta, gmix_term, Pearson_bool, Deviance_bool);
    //----------------------------------------------------------------------------------------------------------------//
    return res;
}

//' Interface between R code and the Cox PH omnibus regression function
//'
//' \code{cox_ph_multidose_Omnibus_transition} Called directly from R, Defines the control variables and calls the regression function
//' @inheritParams CPP_template
//'
//' @return LogLik_Cox_PH output : Log-likelihood of optimum, first derivative of log-likelihood, second derivative matrix, parameter list, standard deviation estimate, AIC, model information
//' @noRd
//'
//  [[Rcpp::export]]
List cox_ph_multidose_Omnibus_transition(IntegerVector term_n, StringVector tform, NumericVector a_n, IntegerMatrix dose_cols, IntegerVector dose_index, IntegerVector dfc, MatrixXd df0, MatrixXd df1, int fir, string modelform, List Control, const MatrixXd df_m, NumericVector tu, IntegerVector KeepConstant, int term_tot, NumericVector Strata_vals, VectorXd cens_weight, List model_control, MatrixXd Lin_Sys, VectorXd Lin_Res) {
    int verbose = Control["verbose"];
    //
    double lr = Control["lr"];
    int maxiter = Control["maxiter"];
    int halfmax = Control["halfmax"];
    double epsilon = Control["epsilon"];
    //
    double step_max = Control["step_max"];
    double thres_step_max = Control["thres_step_max"];
    double deriv_epsilon = Control["deriv_epsilon"];
    string ties_method = Control["ties"];
    int nthreads = Control["ncores"];
    //
    Map<VectorXd> beta_0(as<Map<VectorXd> >(a_n));
    //
    double gmix_theta = model_control["gmix_theta"];
    IntegerVector gmix_term = model_control["gmix_term"];
    //
    List model_bool = List::create(
            _["strata"] = model_control["strata"],
            _["basic"] = model_control["basic"],
            _["linear_err"] = model_control["linear_err"],
            _["null"] = model_control["null"],
            _["cr"] = model_control["cr"],
            _["single"] = model_control["single"],
            _["gradient"] = model_control["gradient"],
            _["outcome_prob"] = false,
            _["constraint"] = model_control["constraint"],
            _["observed_info"] = model_control["observed_info"],
            _["log_bound"] = false,
            _["cox"] = true);
    List optim_para = List::create(
            _["lr"] = Control["lr"]);
    if (model_bool["gradient"]) {
        model_bool["momentum"] = model_control["momentum"];
        model_bool["adadelta"] = model_control["adadelta"];
        model_bool["adam"] = model_control["adam"];
        optim_para["momentum_decay"] = model_control["momentum_decay"];
        optim_para["learning_decay"] = model_control["learning_decay"];
        optim_para["epsilon_decay"] = model_control["epsilon_decay"];
        if (model_bool["constraint"]) {
            optim_para["penalty_weight"] = model_control["penalty_weight"];
            optim_para["penalty_method"] = model_control["penalty_method"];
        }
    }
    //
    //  Performs regression
    //----------------------------------------------------------------------------------------------------------------//
    List res;
    bool IntegratedSerial = model_control["mcml"];
    if (IntegratedSerial) {
        res = LogLik_Cox_PH_Multidose_Omnibus_Integrated(term_n, tform, beta_0, df0, df1, dose_cols, dose_index, dfc, fir, modelform, lr, optim_para, maxiter, halfmax, epsilon, step_max, thres_step_max, deriv_epsilon, df_m, tu, verbose, KeepConstant, term_tot, ties_method, nthreads, Strata_vals, cens_weight, model_bool, gmix_theta, gmix_term, Lin_Sys, Lin_Res);
    } else {
        res = LogLik_Cox_PH_Multidose_Omnibus_Serial(term_n, tform, beta_0, df0, df1, dose_cols, dose_index, dfc, fir, modelform, lr, optim_para, maxiter, halfmax, epsilon, step_max, thres_step_max, deriv_epsilon, df_m, tu, verbose, KeepConstant, term_tot, ties_method, nthreads, Strata_vals, cens_weight, model_bool, gmix_theta, gmix_term, Lin_Sys, Lin_Res);
    }
    //----------------------------------------------------------------------------------------------------------------//
    return res;
}

//' Interface between R code and the poisson multidose omnibus regression function
//'
//' \code{pois_multidose_Omnibus_transition} Called directly from R, Defines the control variables and calls the regression function
//' @inheritParams CPP_template
//'
//' @return LogLik_Pois output : Log-likelihood of optimum, first derivative of log-likelihood, second derivative matrix, parameter list, standard deviation estimate, AIC, model information
//' @noRd
//'
//  [[Rcpp::export]]
List pois_multidose_Omnibus_transition(MatrixXd PyrC, IntegerVector term_n, StringVector tform, NumericVector a_n, IntegerMatrix dose_cols, IntegerVector dose_index, IntegerVector dfc, MatrixXd df0, MatrixXd df1, int fir, string modelform, List Control, IntegerVector KeepConstant, int term_tot, MatrixXd dfs, List model_control, MatrixXd Lin_Sys, VectorXd Lin_Res) {
    int verbose = Control["verbose"];
    //
    double lr = Control["lr"];
    int maxiter = Control["maxiter"];
    int halfmax = Control["halfmax"];
    double epsilon = Control["epsilon"];
    //
    double step_max = Control["step_max"];
    double thres_step_max = Control["thres_step_max"];
    double deriv_epsilon = Control["deriv_epsilon"];
    int nthreads = Control["ncores"];
    //
    double gmix_theta = model_control["gmix_theta"];
    IntegerVector gmix_term = model_control["gmix_term"];
    Map<VectorXd> beta_0(as<Map<VectorXd> >(a_n));
    //
    List model_bool = List::create(
            _["strata"] = model_control["strata"],
            _["basic"] = false,
            _["linear_err"] = false,
            _["null"] = false,
            _["cr"] = false,
            _["single"] = model_control["single"],
            _["gradient"] = model_control["gradient"],
            _["outcome_prob"] = false,
            _["constraint"] = model_control["constraint"],
            _["observed_info"] = model_control["observed_info"],
            _["log_bound"] = false,
            _["cox"] = false);
    List optim_para = List::create(
            _["lr"] = Control["lr"]);
    if (model_bool["gradient"]) {
        model_bool["momentum"] = model_control["momentum"];
        model_bool["adadelta"] = model_control["adadelta"];
        model_bool["adam"] = model_control["adam"];
        optim_para["momentum_decay"] = model_control["momentum_decay"];
        optim_para["learning_decay"] = model_control["learning_decay"];
        optim_para["epsilon_decay"] = model_control["epsilon_decay"];
        if (model_bool["constraint"]) {
            optim_para["penalty_weight"] = model_control["penalty_weight"];
            optim_para["penalty_method"] = model_control["penalty_method"];
        }
    }
    //
    //  Performs regression
    //----------------------------------------------------------------------------------------------------------------//
    List res;
    bool IntegratedSerial = model_control["mcml"];
    if (IntegratedSerial) {
        res = LogLik_Pois_PH_Multidose_Omnibus_Integrated(PyrC, term_n, tform, beta_0, df0, df1, dose_cols, dose_index, dfc, fir, modelform, lr, optim_para, maxiter, halfmax, epsilon, step_max, thres_step_max, deriv_epsilon, dfs, verbose, KeepConstant, term_tot, nthreads, model_bool, gmix_theta, gmix_term, Lin_Sys, Lin_Res);
    } else {
        res = LogLik_Pois_PH_Multidose_Omnibus_Serial(PyrC, term_n, tform, beta_0, df0, df1, dose_cols, dose_index, dfc, fir, modelform, lr, optim_para, maxiter, halfmax, epsilon, step_max, thres_step_max, deriv_epsilon, dfs, verbose, KeepConstant, term_tot, nthreads, model_bool, gmix_theta, gmix_term, Lin_Sys, Lin_Res);
    }
    //----------------------------------------------------------------------------------------------------------------//
    return res;
}

//' Interface between R code and the matched case-control omnibus regression function
//'
//' \code{caco_Omnibus_transition} Called directly from R, Defines the control variables and calls the regression function
//' @inheritParams CPP_template
//'
//' @return LogLik_CaseCon_Omnibus output : Log-likelihood of optimum, first derivative of log-likelihood, second derivative matrix, parameter list, standard deviation estimate, AIC, model information
//' @noRd
//'
//  [[Rcpp::export]]
List caco_Omnibus_transition(IntegerVector term_n, StringVector tform, NumericMatrix& a_ns, IntegerVector dfc, MatrixXd df0, int fir, string modelform, List Control, const MatrixXd df_m, NumericVector tu, IntegerVector KeepConstant, int term_tot, NumericVector Strata_vals, List model_control, MatrixXd Lin_Sys, VectorXd Lin_Res) {
    int verbose = Control["verbose"];
    //
    double lr = Control["lr"];
    NumericVector maxiters = Control["maxiters"];
    int guesses = Control["guesses"];
    int halfmax = Control["halfmax"];
    double epsilon = Control["epsilon"];
    //
    double step_max = Control["step_max"];
    double thres_step_max = Control["thres_step_max"];
    double deriv_epsilon = Control["deriv_epsilon"];
    string ties_method = Control["ties"];
    int nthreads = Control["ncores"];
    //
    double gmix_theta = model_control["gmix_theta"];
    IntegerVector gmix_term = model_control["gmix_term"];
    //
    List model_bool = List::create(
            _["strata"] = model_control["strata"],
            _["time_risk"] = model_control["time_risk"],
            _["basic"] = model_control["basic"],
            _["linear_err"] = model_control["linear_err"],
            _["null"] = model_control["null"],
            _["cr"] = false,
            _["single"] = model_control["single"],
            _["gradient"] = model_control["gradient"],
            _["outcome_prob"] = false,
            _["cond_thres"] = model_control["conditional_threshold"],
            _["constraint"] = model_control["constraint"],
            _["observed_info"] = model_control["observed_info"],
            _["log_bound"] = false,
            _["cox"] = false);
    List optim_para = List::create(
            _["lr"] = Control["lr"]);
    if (model_bool["gradient"]) {
        model_bool["momentum"] = model_control["momentum"];
        model_bool["adadelta"] = model_control["adadelta"];
        model_bool["adam"] = model_control["adam"];
        optim_para["momentum_decay"] = model_control["momentum_decay"];
        optim_para["learning_decay"] = model_control["learning_decay"];
        optim_para["epsilon_decay"] = model_control["epsilon_decay"];
        if (model_bool["constraint"]) {
            optim_para["penalty_weight"] = model_control["penalty_weight"];
            optim_para["penalty_method"] = model_control["penalty_method"];
        }
    }
    //
    //  Performs regression
    //----------------------------------------------------------------------------------------------------------------//
    List res = LogLik_CaseCon_Omnibus(term_n, tform, a_ns, df0, dfc, fir, modelform, lr, optim_para, maxiters, guesses, halfmax, epsilon, step_max, thres_step_max, deriv_epsilon, df_m, tu, verbose, KeepConstant, term_tot, ties_method, nthreads, Strata_vals, model_bool, gmix_theta, gmix_term, Lin_Sys, Lin_Res);
    //----------------------------------------------------------------------------------------------------------------//
    return res;
}

//' Interface between R code and the logistic omnibus regression function
//'
//' \code{logist_Omnibus_transition} Called directly from R, Defines the control variables and calls the regression function
//' @inheritParams CPP_template
//'
//' @return LogLik_Cox_PH output : Log-likelihood of optimum, first derivative of log-likelihood, second derivative matrix, parameter list, standard deviation estimate, AIC, model information
//' @noRd
//'
//  [[Rcpp::export]]
List logist_Omnibus_transition(MatrixXd CountEvent, IntegerVector term_n, StringVector tform, NumericMatrix& a_ns, IntegerVector dfc, MatrixXd df0, int fir, string modelform, List Control, IntegerVector KeepConstant, int term_tot, List model_control, MatrixXd Lin_Sys, VectorXd Lin_Res) {
    //
    int verbose = Control["verbose"];
    //
    double lr = Control["lr"];
    NumericVector maxiters = Control["maxiters"];
    int guesses = Control["guesses"];
    int halfmax = Control["halfmax"];
    double epsilon = Control["epsilon"];
    //
    double step_max = Control["step_max"];
    double thres_step_max = Control["thres_step_max"];
    double deriv_epsilon = Control["deriv_epsilon"];
    int nthreads = Control["ncores"];
    //
    double gmix_theta = model_control["gmix_theta"];
    IntegerVector gmix_term = model_control["gmix_term"];
    //
    List model_bool = List::create(
            _["strata"] = model_control["strata"],
            _["basic"] = false,
            _["linear_err"] = false,
            _["null"] = false,
            _["cr"] = false,
            _["single"] = model_control["single"],
            _["gradient"] = model_control["gradient"],
            _["outcome_prob"] = false,
            _["constraint"] = model_control["constraint"],
            _["observed_info"] = model_control["observed_info"],
            _["log_bound"] = false,
            _["cox"] = false,
            _["odds"] = model_control["logit_odds"],
            _["ident"] = model_control["logit_ident"],
            _["loglink"] = model_control["logit_loglink"]);
    List optim_para = List::create(
            _["lr"] = Control["lr"]);
    if (model_bool["gradient"]) {
        model_bool["momentum"] = model_control["momentum"];
        model_bool["adadelta"] = model_control["adadelta"];
        model_bool["adam"] = model_control["adam"];
        optim_para["momentum_decay"] = model_control["momentum_decay"];
        optim_para["learning_decay"] = model_control["learning_decay"];
        optim_para["epsilon_decay"] = model_control["epsilon_decay"];
        if (model_bool["constraint"]) {
            optim_para["penalty_weight"] = model_control["penalty_weight"];
            optim_para["penalty_method"] = model_control["penalty_method"];
        }
    }
    //
    //  Performs regression
    //----------------------------------------------------------------------------------------------------------------//
    List res = LogLik_Logist_Omnibus(CountEvent, term_n, tform, a_ns, df0, dfc, fir, modelform, lr, optim_para, maxiters, guesses, halfmax, epsilon, step_max, thres_step_max, deriv_epsilon, verbose, KeepConstant, term_tot, nthreads, model_bool, gmix_theta, gmix_term, Lin_Sys, Lin_Res);
    //----------------------------------------------------------------------------------------------------------------//
    return res;
}


//' Generates csv file with time-dependent columns
//'
//' \code{Write_Time_Dep} Called directly from R, Defines a new matrix which interpolates time-dependent values on a grid
//' @inheritParams CPP_template
//'
//' @return saves a dataframe to be used with time-dependent covariate analysis
//' @noRd
//'
//  [[Rcpp::export]]
void Write_Time_Dep(const NumericMatrix df0_Times, const NumericMatrix df0_dep, const NumericMatrix df0_const, const NumericVector df0_event, double dt, string filename, StringVector tform_tdep, NumericVector tu, bool iscox, int nthreads) {
    const Map<MatrixXd> df_Times(as<Map<MatrixXd> >(df0_Times));
    const Map<MatrixXd> df_dep(as<Map<MatrixXd> >(df0_dep));
    const Map<MatrixXd> df_const(as<Map<MatrixXd> >(df0_const));
    Rcout.precision(10);  //  forces higher precision numbers printed to terminal
    int tot_covs = ceil(2 + df_dep.cols()/2 + df_const.cols() + 1);
    int max_rows = 0;
    if (iscox) {
        max_rows = tu.size();
        dt = tu[1] - tu[0];
        for (int i  = 0; i < tu.size() - 1; i++) {
            if (dt > (tu[i + 1] - tu[i])) {
                dt = tu[i + 1] - tu[i];
            }
        }
    } else {
        max_rows = ceil(((df_Times.col(1).array() - df_Times.col(0).array()).array().abs().maxCoeff()) / dt);
    }
    int True_Rows = 0;
    VectorXd row_store = VectorXd::Zero(tot_covs);
    MatrixXd new_row_store = MatrixXd::Zero(max_rows, tot_covs);
    //
    static const IOFormat CSVFormat(FullPrecision, DontAlignCols, ", ", "\n");
    ofstream file(filename);
    //
    int serial_0 = 0;
    int serial_1 = 0;
    for (int i_row = 0; i_row < df_Times.rows(); i_row++) {
        new_row_store = MatrixXd::Zero(max_rows, tot_covs);
        if (iscox) {
            True_Rows = 0;
            serial_0 = 0;
            serial_1 = 0;
            for (int i = 0; i < tu.size(); i++) {
                if (df_Times.coeff(i_row, 1) >= tu[i]) {
                    serial_1 = i;
                }
                if (df_Times.coeff(i_row, 0) > tu[i]) {
                    serial_0 = i + 1;
                }
            }
            True_Rows = serial_1 - serial_0 + 1;
            #ifdef _OPENMP
            #pragma omp parallel for schedule(dynamic) num_threads(nthreads)
            #endif
            for (int i_inner=serial_0; i_inner < serial_1 + 1; i_inner++) {
                VectorXd dep_temp = VectorXd::Zero(df_dep.cols()/2);
                double t0 = tu[i_inner]- dt/2;
                double t1 = tu[i_inner];
                double ratio = 0;
                if ((df_Times.coeff(i_row, 1) - df_Times.coeff(i_row, 0)) > 0) {
                    ratio = (t1 - df_Times.coeff(i_row, 0))/(df_Times.coeff(i_row, 1) - df_Times.coeff(i_row, 0));
                }
                string func_id = "";
                string delim = "?";
                size_t pos = 0;
                string token = "";
                double temp_tok  = 0;
                int gather_val = 0;
                char tok_char = 'a';
                //  step string is either g, l, a, b for >=, <=, >, <
                for (int i = 0;  i < dep_temp.size(); i++) {
                    func_id = as<string>(tform_tdep[i]);
                    if (func_id == "lin") {
                        dep_temp[i] = ratio * df_dep.coeff(i_row, 2*i + 1) + (1 - ratio) * df_dep.coeff(i_row, 2*i);
                    } else {
                        pos = func_id.find(delim);
                        token = func_id.substr(0, pos);
                        if (token == "step") {
                            func_id.erase(0, pos + delim.length());
                            gather_val = 0;
                            while ((pos = func_id.find(delim)) != string::npos) {
                                token = func_id.substr(0, pos);
                                //
                                tok_char = token[token.length() - 1];
                                if (tok_char == 'g') {
                                    token.pop_back();
                                    temp_tok = stod(token);
                                    if (t1 > temp_tok) {
                                        gather_val = gather_val + 1;
                                    }
                                    if (t1 == temp_tok) {
                                        gather_val = gather_val + 1;
                                    }
                                } else if (tok_char == 'l') {
                                    token.pop_back();
                                    temp_tok = stod(token);
                                    if (t1 < temp_tok) {
                                        gather_val = gather_val + 1;
                                    }
                                    if (t1 == temp_tok) {
                                        gather_val = gather_val + 1;
                                    }
                                } else if (tok_char == 'a') {
                                    token.pop_back();
                                    temp_tok = stod(token);
                                    if (t1 > temp_tok) {
                                        gather_val = gather_val + 1;
                                    }
                                } else if (tok_char == 'b') {
                                    token.pop_back();
                                    temp_tok = stod(token);
                                    if (t1 < temp_tok) {
                                        gather_val = gather_val + 1;
                                    }
                                } else {
                                }
                                //
                                func_id.erase(0, pos + delim.length());
                            }
                            dep_temp[i] = gather_val;
                        } else {
                            Rcout << "C++ Error: " << func_id << " _:_ " << token << endl;
                            throw invalid_argument("time dependent identifier isn't implemented");
                        }
                    }
                }
                int event0 = 0;
                if (i_inner == True_Rows - 1) {
                    event0 = df0_event[i_row];
                }
                new_row_store.row(i_inner) << t0, t1, dep_temp.transpose(), df_const.row(i_row), event0;
            }
        } else {
            True_Rows = ceil((df_Times.coeff(i_row, 1) - df_Times.coeff(i_row, 0))/dt);
            #ifdef _OPENMP
            #pragma omp parallel for schedule(dynamic) num_threads(nthreads)
            #endif
            for (int i_inner = 0; i_inner < True_Rows; i_inner++) {
                VectorXd dep_temp = VectorXd::Zero(df_dep.cols()/2);
                double ratio = (i_inner+0.5) / True_Rows;
                double t0 = df_Times.coeff(i_row, 0) + i_inner * dt;
                double t1 = t0 + dt;
                string func_id = "";
                string delim = "?";
                size_t pos = 0;
                string token = "";
                double temp_tok  = 0;
                int gather_val = 0;
                char tok_char = 'a';
                for (int i = 0;  i < dep_temp.size(); i++) {
                    func_id = as<string>(tform_tdep[i]);
                    if (func_id == "lin") {
                        dep_temp[i] = ratio * df_dep.coeff(i_row, 2*i + 1) + (1 - ratio) * df_dep.coeff(i_row, 2*i);
                    } else {
                        pos = func_id.find(delim);
                        token = func_id.substr(0, pos);
                        if (token == "step") {
                            func_id.erase(0, pos + delim.length());
                            gather_val = 0;
                            while ((pos = func_id.find(delim)) != string::npos) {
                                token = func_id.substr(0, pos);
                                //
                                tok_char = token[token.length() - 1];
                                if (tok_char == 'u') {
                                    token.pop_back();
                                    temp_tok = stod(token);
                                    if (t0 >= temp_tok) {
                                        gather_val = gather_val + 1;
                                    }
                                } else if (tok_char == 'l') {
                                    token.pop_back();
                                    temp_tok = stod(token);
                                    if (t1 <= temp_tok) {
                                        gather_val = gather_val + 1;
                                    }
                                } else {
                                }
                                //
                                func_id.erase(0, pos + delim.length());
                            }
                            dep_temp[i] = gather_val;
                        }
                    }
                }
                int event0 = 0;
                if (i_inner == True_Rows - 1) {
                    t1 = df_Times.coeff(i_row, 1);
                    event0 = df0_event[i_row];
                }
                new_row_store.row(i_inner) << t0 + (t1-t0)*1e-1, t1, dep_temp, df_const.row(i_row), event0;
            }
        }
        if (file.is_open()) {
            file << new_row_store.block(0, 0, True_Rows, tot_covs).format(CSVFormat);
            file << "\n";
        }
    }
    if (file.is_open()) {
        file.close();
    }
}

//' Generates weightings for stratified poisson regression
//'
//' \code{Gen_Strat_Weight} Called from within c++, assigns vector of weights
//' @inheritParams CPP_template
//'
//' @return assigns weight in place and returns nothing
//' @noRd
//'
void Gen_Strat_Weight(string modelform, const Ref<const MatrixXd>& dfs, const Ref<const MatrixXd>& PyrC, VectorXd& s_weights, const int nthreads, const StringVector& tform, const IntegerVector& term_n, const int& term_tot, const double gmix_theta, const IntegerVector& gmix_term) {
    ArrayXd Pyrs  = dfs.transpose() * PyrC.col(0);
    ArrayXd Events = dfs.transpose() * PyrC.col(1);
    ArrayXd weight = Events.array() * Pyrs.array().pow(- 1).array();
    //
    //
    s_weights = dfs * weight.matrix();
    //
    vector<int> lin_count(term_tot, 0);  //  tracking which terms will go to 0 for only being linear
    vector<int> dose_count(term_tot, 0);  // tracking which terms will be a sum of 1s, for being dose non-piecewise
    vector<int> dose_lin_count(term_tot, 0);  // tracking which terms will go to 0 for being dose-piecewise
    for (int ij = 0; ij < (term_n.size()); ij++) {
        int tn = term_n[ij];
        if (as<string>(tform[ij]) == "loglin") {  //  setting parameters to zero makes the subterm 1
        } else if (as<string>(tform[ij]) == "lin") {  //  setting parameters to zero makes the subterm 0
            lin_count[tn] = lin_count[tn] + 1.0;
        } else if (as<string>(tform[ij]) == "plin") {  //  setting parameters to zero makes the subterm 1
        } else if (as<string>(tform[ij]) == "loglin_slope") {  //  the slope paremeter sets the element to 0
        } else if (as<string>(tform[ij]) == "loglin_top") {  //  the top parameter sets the element to 1
            if (ij == 0) {
                dose_count[tn] = dose_count[tn] + 1.0;
            } else if (tform[ij - 1] != "loglin_slope") {
                dose_count[tn] = dose_count[tn] + 1.0;
            } else {}
        } else if (as<string>(tform[ij]) == "lin_slope") {  //  every other dose term sets the elements to 0
            dose_lin_count[tn] = dose_lin_count[tn] + 1;
        } else if (as<string>(tform[ij]) == "lin_int") {
        } else if (as<string>(tform[ij]) == "quad_slope") {
            dose_lin_count[tn] = dose_lin_count[tn] + 1;
        } else if (as<string>(tform[ij]) == "step_slope") {
            dose_lin_count[tn] = dose_lin_count[tn] + 1;
        } else if (as<string>(tform[ij]) == "step_int") {
        } else if (as<string>(tform[ij]) == "lin_quad_slope") {
            dose_lin_count[tn] = dose_lin_count[tn] + 1;
        } else if (as<string>(tform[ij]) == "lin_quad_int") {
        } else if (as<string>(tform[ij]) == "lin_exp_slope") {
            dose_lin_count[tn] = dose_lin_count[tn] + 1;
        } else if (as<string>(tform[ij]) == "lin_exp_int") {
        } else if (as<string>(tform[ij]) == "lin_exp_exp_slope") {
        } else {
            throw invalid_argument("incorrect subterm type");
        }
    }
    //
    vector<double> term_val(term_tot, 0);
    for (int ijk = 0;  ijk < term_tot; ijk++) {
        if (dose_count[ijk] == 0) {  //  If the dose term isn't used
            if (dose_lin_count[ijk] == 0) {  // If no dose terms that default to 0 are used
                dose_count[ijk] = 1.0;  //  the default term value becomes 1
            }
            //  otherwise the default term value is 0
        }
        if (lin_count[ijk] == 0) {  //  if the linear term isn't used, the entire term is 1 times the dose term value, accounting for the piecewise dose values
            term_val[ijk] = dose_count[ijk];
        } else {  //  if the linear term is used, the entire term is 0
            term_val[ijk] = 0;
        }
    }
    double default_val = 0;
    if (modelform == "A") {
        for (int i = 0;  i < term_tot; i++) {
            default_val += term_val[i];
        }
    } else if (modelform == "PA") {
        for (int i=1; i < term_tot; i++) {
            default_val += term_val[i];
        }
        default_val *= term_val[0];
    } else if (modelform == "PAE") {
        for (int i=1; i < term_tot; i++) {
            default_val += term_val[i];
        }
        default_val = (1 + default_val) * term_val[0];
    } else if (modelform == "M") {
        default_val = 1;
        for (int i = 1; i < term_tot; i++) {
            default_val *= term_val[i];
        }
        default_val *= term_val[0];
    } else if (modelform == "ME") {
        default_val = 1;
        for (int i = 1; i < term_tot; i++) {
            default_val *= (1 + term_val[i]);
        }
        default_val *= term_val[0];
    } else if (modelform == "GMIX") {
        double Ta = 1;
        double Tm = 1;
        for (int i = 1; i < term_tot; i++) {
            Ta += (term_val[i] + gmix_term[i] - 1);
            Tm *= (term_val[i] + gmix_term[i]);
        }
        default_val = term_val[0] * pow(Tm, gmix_theta) * pow(Ta, 1-gmix_theta);
    } else {
        throw invalid_argument("Model isn't implemented");
    }
    s_weights = s_weights / default_val;
    return;
}

//' Checks the OMP flag
//'
//' \code{OMP_Check} Called directly from R, checks the omp flag and returns true if omp is enabled
//'
//' @return boolean: True for OMP allowed
//'
//  [[Rcpp::export]]
bool OMP_Check() {
    bool res = false;
    #ifdef _OPENMP
        res = true;
    #endif
    return res;
}
