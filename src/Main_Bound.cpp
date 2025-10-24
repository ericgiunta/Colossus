//  Copyright 2022 - 2025, Eric Giunta and the project collaborators, Please see main R package for license and usage details

#include <RcppEigen.h>

#include "Main_Bound.h"
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

#include "Omnibus_Pieces.h"
#include "Calc_Repeated.h"
#include "Subterms_Risk.h"
#include "Step_Calc.h"
#include "Step_Grad.h"
#include "Step_Newton.h"
#include "Colossus_types.h"


//  [[Rcpp::depends(RcppEigen)]]
//  [[Rcpp::plugins(openmp)]]

using std::endl;
using std::string;
using std::vector;
using std::invalid_argument;

using Eigen::Map;
using Eigen::Ref;
using Eigen::MatrixXd;
using Eigen::SparseMatrix;
using Eigen::VectorXd;

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
using Rcpp::Dimension;

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

//' Primary Cox PH likelihood bounds calcualtion function.
//'
//' \code{LogLik_Cox_PH_Omnibus_Log_Bound} Performs the calls to calculation functions and log-likeihood profile bounds
//'
//' @inheritParams CPP_template
//'
//' @return List of final results: Log-likelihood of optimum, first derivative of log-likelihood, second derivative matrix, parameter list, standard deviation estimate, AIC, model information
//' @noRd
//'
//
List LogLik_Cox_PH_Omnibus_Log_Bound(IntegerVector term_n, StringVector tform, Ref<VectorXd> beta_0, Ref<MatrixXd> df0, IntegerVector dfc, int fir, string modelform, double lr, NumericVector maxiters, int guesses, int halfmax, double epsilon, double step_max, double thres_step_max, double deriv_epsilon, const Ref<const MatrixXd>& df_m, NumericVector tu, int verbose, IntegerVector KeepConstant, int term_tot, string ties_method, int nthreads, NumericVector& Strata_vals, const VectorXd& cens_weight, List model_bool, const double gmix_theta, const IntegerVector gmix_term, const Ref<const MatrixXd>& Lin_Sys, const Ref<const VectorXd>& Lin_Res, double qchi, int para_number, int maxstep, double mult) {
    //
    List temp_list = List::create(_["Status"] = "TEMP");  //  used as a dummy return value for code checking
    //  Time durations are measured from this point on in microseconds
    //
    //  df0: covariate data
    //  ntime: number of event times for Cox PH
    //  totalnum: number of terms used
    //
    //  ------------------------------------------------------------------------- //  initialize
//    const int mat_row = df0.rows();
    int ntime = tu.size();
    int totalnum = term_n.size();
    int reqrdnum = totalnum - sum(KeepConstant);
    //  ------------------------------------------------------------------------- //  initialize
    if (model_bool["null"]) {
        if (verbose >= 1) {
            Rcout << "null model is not compatable with log-based bound calculation" << endl;
        }
        temp_list = List::create(_["Status"] = "FAILED_BAD_MODEL_NULL", _["LogLik"] = R_NaN);
        return temp_list;
    }
    if (model_bool["single"]) {
        if (verbose >= 1) {
            Rcout << "non-derivative model calculation is not compatable with log-based bound calculation" << endl;
        }
        temp_list = List::create(_["Status"] = "FAILED_WITH_BAD_MODEL_SINGLE", _["LogLik"] = R_NaN);
        return temp_list;
    }
    if (model_bool["gradient"]) {
        if (verbose >= 1) {
            Rcout << "gradient descent model calculation is not compatable with log-based bound calculation" << endl;
        }
        temp_list = List::create(_["Status"] = "FAILED_WITH_BAD_MODEL_GRADIENT", _["LogLik"] = R_NaN);
        return temp_list;
    }
    //
    //  cout.precision: controls the number of significant digits printed
    //
    Rcout.precision(7);  //  forces higher precision numbers printed to terminal
    //
    //  Lld_worst: stores the highest magnitude log-likelihood derivative
    //  ---------------------------------------------
    //  To Start, needs to seperate the derivative terms
    //  ---------------------------------------------
    //
    //  ------------------------------------------------------------------------- //  initialize
    // Map<VectorXd> beta_0(as<Map<VectorXd> >(a_n));
    MatrixXd T0;
    MatrixXd Td0;
    MatrixXd Tdd0;
    //
    MatrixXd Te;
    MatrixXd R;
    MatrixXd Rd;
    MatrixXd Rdd;
    //
    MatrixXd Dose;
    MatrixXd nonDose;
    MatrixXd nonDose_LIN;
    MatrixXd nonDose_PLIN;
    MatrixXd nonDose_LOGLIN;
    MatrixXd TTerm;
    double dint = 0.0;  //  the amount of change used to calculate derivatives in threshold paramters
    double dslp = 0.0;
    MatrixXd RdR;
    MatrixXd RddR;
    //  ------------------------------------------------------------------------- //  initialize
    //  ---------------------------------------------
    //  To Start, needs to seperate the derivative terms
    //  ---------------------------------------------
    //
    Cox_Refresh_R_TERM(totalnum, reqrdnum, term_tot, dint, dslp, thres_step_max, step_max, df0, T0, Td0, Tdd0, Te, R, Rd, Rdd, Dose, nonDose, TTerm, nonDose_LIN, nonDose_PLIN, nonDose_LOGLIN, RdR, RddR, model_bool);
    //  ------------------------------------------------------------------------- //  initialize
    IntegerMatrix RiskFail;
    vector<vector<int> > RiskPairs(ntime);
    vector<vector<vector<int> > > RiskPairs_Strata(ntime, vector<vector<int>>(Strata_vals.size()));
    //  ------------------------------------------------------------------------- //  initialize
    if (model_bool["strata"]) {
        RiskFail = IntegerMatrix(ntime, 2*Strata_vals.size());  //  vector giving the event rows
        //
        //  Creates matrices used to identify the event risk groups
        if (model_bool["cr"]) {
            Make_Groups_Strata_CR(ntime, df_m, RiskFail, RiskPairs_Strata, tu, nthreads, Strata_vals, cens_weight);
        } else {
            Make_Groups_Strata(ntime, df_m, RiskFail, RiskPairs_Strata, tu, nthreads, Strata_vals);
        }
    } else {
        RiskFail = IntegerMatrix(ntime, 2);  //  vector giving the event rows
        //
        //  Creates matrices used to identify the event risk groups
        if (model_bool["cr"]) {
            Make_Groups_CR(ntime, df_m, RiskFail, RiskPairs, tu, cens_weight, nthreads);
        } else {
            Make_Groups(ntime, df_m, RiskFail, RiskPairs, tu, nthreads);
        }
    }
    //  ------------------------------------------------------------------------- //  initialize
    MatrixXd Rls1;
    MatrixXd Lls1;
    MatrixXd Rls2;
    MatrixXd Rls3;
    MatrixXd Lls2;
    MatrixXd Lls3;
    vector<double> Ll(reqrdnum, 0.0);  //  log-likelihood values
    vector<double> Lld(reqrdnum, 0.0);  //  log-likelihood derivative values
    vector<double> Lldd(pow(reqrdnum, 2), 0.0);  //  the second derivative matrix has room for every combination, but only the lower triangle is calculated initially
    //  ------------------------------------------------------------------------- //  initialize
    Cox_Refresh_R_SIDES(reqrdnum, ntime, Rls1, Rls2, Rls3, Lls1, Lls2, Lls3, Strata_vals, model_bool);
    //  the log-likelihood is calculated in parallel over the risk groups
    vector <double> Ll_comp(2, Ll[0]);  //  vector to compare values
    double step_max0 = step_max;
    double thres_step_max0 = thres_step_max;
    //
    vector<double> dbeta(totalnum, 0.0);
    //
    //  --------------------------
    //  always starts from initial guess
    //  --------------------------
    vector<double> beta_c(totalnum, 0.0);
    vector<double> beta_a(totalnum, 0.0);
    vector<double> beta_best(totalnum, 0.0);
    vector<double> beta_peak(totalnum, 0.0);
    vector<double> beta_p(totalnum, 0.0);
    VectorXd::Map(&beta_p[0], beta_0.size()) = beta_0;  //  stores previous parameters
    VectorXd::Map(&beta_c[0], beta_0.size()) = beta_0;  //  stores current parameters
    VectorXd::Map(&beta_a[0], beta_0.size()) = beta_0;  //  stores a refrence value for parameters
    VectorXd::Map(&beta_best[0], beta_0.size()) = beta_0;  //  stores the best parameters
    VectorXd::Map(&beta_peak[0], beta_0.size()) = beta_0;  //  stores the best parameters
    int iter_stop  = 0;  //  tracks if the iterations should be stopped for convergence
    //
    //
    vector<double> beta_abs_best(totalnum, 0.0);
    //
    Cox_Refresh_R_SIDES(reqrdnum, ntime, Rls1, Rls2, Rls3, Lls1, Lls2, Lls3, Strata_vals, model_bool);
    fill(Ll.begin(), Ll.end(), 0.0);
    fill(Lld.begin(), Lld.end(), 0.0);
    fill(Lldd.begin(), Lldd.end(), 0.0);
    beta_p = beta_peak;  //
    beta_a = beta_peak;  //
    beta_c = beta_peak;  //
    step_max = step_max0;
    thres_step_max = thres_step_max0;
    //
    for (int i = 0; i < beta_0.size(); i++) {
        beta_0[i] = beta_peak[i];
    }
    Cox_Term_Risk_Calc(modelform, tform, term_n, totalnum, fir, dfc, term_tot, T0, Td0, Tdd0, Te, R, Rd, Rdd, Dose, nonDose, beta_0, df0, thres_step_max, step_max, TTerm, nonDose_LIN, nonDose_PLIN, nonDose_LOGLIN, RdR, RddR, nthreads, KeepConstant, verbose, model_bool, gmix_theta, gmix_term);
    Cox_Side_LL_Calc(reqrdnum, ntime, tform, RiskFail,  RiskPairs, RiskPairs_Strata, totalnum, fir, R, Rd, Rdd, Rls1, Rls2, Rls3, Lls1, Lls2, Lls3, cens_weight, Strata_vals, beta_0, RdR, RddR, Ll, Lld, Lldd, nthreads, KeepConstant, ties_method, verbose, model_bool, iter_stop);
    //
    NumericVector Lldd_vec(reqrdnum * reqrdnum);
    NumericVector Lld_vecc(reqrdnum);
    #ifdef _OPENMP
    #pragma omp parallel for schedule(dynamic) num_threads(nthreads)
    #endif
    for (int ijk = 0; ijk < reqrdnum*(reqrdnum + 1)/2; ijk++) {
        int ij = 0;
        int jk = ijk;
        while (jk > ij) {
            ij++;
            jk -= ij;
        }
        Lldd_vec[ij * reqrdnum + jk] = Lldd[ij*reqrdnum+jk];
        Lldd_vec[jk * reqrdnum + ij] = Lldd_vec[ij * reqrdnum + jk];
        if (ij == jk) {
            Lld_vecc[ij] = Lld[ij];
        }
    }
    Lldd_vec.attr("dim") = Dimension(reqrdnum, reqrdnum);
    Map<MatrixXd> Lldd_mat(as<Map<MatrixXd> >(Lldd_vec));
    Map<VectorXd> Lld_vec(as<Map<VectorXd> >(Lld_vecc));
    double Lstar = Ll[0]-qchi;
    bool upper = true;
    int half_check = 0;
    bool trouble = false;
    IntegerVector KeepConstant_trouble(totalnum);
    for (int i = 0;  i < totalnum; i++) {
        KeepConstant_trouble[i] = KeepConstant[i];
    }
    KeepConstant_trouble[para_number] = 1;
    //
    vector<double> limits(2, 0.0);
    vector<bool>   limit_hit(2, FALSE);
    vector<bool>   limit_converged(2, FALSE);
    vector<double> ll_final(2, 0.0);
    List res_list;
    //
    if (verbose >= 4) {
        Rcout << "C++ Note: STARTING Upper Bound" << endl;
    }
    upper = true;
    int step = -1;
    bool iter_continue = true;
    double max_change = 100;
    double deriv_max = 100;
    ///
    //  variables added for log loop code
    ///
    VectorXd s_weights(1);
    int bound_val = 1;
    //
    while ((step < maxstep) && (iter_continue)) {
        Lldd_mat = as<Map<MatrixXd> >(Lldd_vec);
        Lld_vec = as<Map<VectorXd> >(Lld_vecc);
        step++;
        trouble = false;
        half_check = 0;
        Log_Bound(deriv_max, Lldd_mat, Lld_vec, Lstar, qchi, Ll[0], para_number, nthreads, totalnum, reqrdnum, KeepConstant, term_tot, step, dbeta, beta_0, upper, trouble, verbose, mult);
        if (trouble) {
            Calc_Change_trouble(para_number, nthreads, totalnum, thres_step_max, lr, step_max, Ll, Lld, Lldd, dbeta, tform, thres_step_max, step_max, KeepConstant_trouble);
        }
        beta_p = beta_c;  //
        beta_a = beta_c;  //
        Cox_Refresh_R_SIDES(reqrdnum, ntime, Rls1, Rls2, Rls3, Lls1, Lls2, Lls3, Strata_vals, model_bool);
        for (int ijk = 0; ijk < totalnum; ijk++) {
            if (KeepConstant[ijk] == 0) {
                if ((tform[ijk] == "lin_quad_int") || (tform[ijk] == "lin_exp_int") || (tform[ijk] == "step_int") || (tform[ijk] == "lin_int")) {  //  the threshold values use different maximum deviation values
                    if (abs(dbeta[ijk]) > thres_step_max) {
                        dbeta[ijk] = thres_step_max * sign(dbeta[ijk]);
                    }
                } else {
                    if (abs(dbeta[ijk]) > step_max) {
                        dbeta[ijk] = step_max * sign(dbeta[ijk]);
                    }
                }
            } else {
                dbeta[ijk] = 0;
            }
        }
        max_change = abs(dbeta[0]);
        for (int ij = 0; ij < totalnum; ij++) {
            if (ij == para_number) {
                //  prevent the parameter estimate from crossing the optimum
                //  issue is beta_0[para_number] <= beta_best[para_number]
                if (dbeta[ij] < (beta_peak[para_number] - beta_a[ij])/lr) {
                    dbeta[ij] = (beta_peak[para_number] - beta_a[ij])/lr/2;
                }
            }
            if (abs(dbeta[ij]) > max_change) {
                max_change = abs(dbeta[ij]);
            }
            beta_0[ij] = beta_a[ij] + lr*dbeta[ij];
            beta_c[ij] = beta_0[ij];
        }
        //  ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------//
        //  The same subterm, risk, sides, and log-likelihood calculations are performed every half-step and iteration
        //  ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------//
        Cox_Term_Risk_Calc(modelform, tform, term_n, totalnum, fir, dfc, term_tot, T0, Td0, Tdd0, Te, R, Rd, Rdd, Dose, nonDose, beta_0, df0, thres_step_max, step_max, TTerm, nonDose_LIN, nonDose_PLIN, nonDose_LOGLIN, RdR, RddR, nthreads, KeepConstant, verbose, model_bool, gmix_theta, gmix_term);
        bound_val = 1;
        Cox_Pois_Log_Loop(step_max, model_bool, beta_0, beta_a, beta_c, bound_val, dbeta, df0, dfc, dint, Dose, thres_step_max, dslp, fir, gmix_term, gmix_theta, half_check, halfmax, KeepConstant, limit_hit, lr, modelform, nonDose, nonDose_LIN, nonDose_LOGLIN, nonDose_PLIN, nthreads, R, Rd, Rdd, RddR, RdR, s_weights, T0, Td0, Tdd0, Te, term_n, term_tot, tform, totalnum, TTerm, verbose);
        if (limit_hit[1]) {
            limits[1] = 0;
            ll_final[1] = 0;
            limit_converged[1] = FALSE;
            break;
        }
        Cox_Side_LL_Calc(reqrdnum, ntime, tform, RiskFail,  RiskPairs, RiskPairs_Strata, totalnum, fir, R, Rd, Rdd, Rls1, Rls2, Rls3, Lls1, Lls2, Lls3, cens_weight, Strata_vals, beta_0, RdR, RddR, Ll, Lld, Lldd, nthreads, KeepConstant, ties_method, verbose, model_bool, iter_stop);
        //
        #ifdef _OPENMP
        #pragma omp parallel for schedule(dynamic) num_threads(nthreads)
        #endif
        for (int ijk = 0; ijk < reqrdnum*(reqrdnum + 1)/2; ijk++) {
            int ij = 0;
            int jk = ijk;
            while (jk > ij) {
                ij++;
                jk -= ij;
            }
            Lldd_vec[ij * reqrdnum + jk] = Lldd[ij*reqrdnum+jk];
            Lldd_vec[jk * reqrdnum + ij] = Lldd_vec[ij * reqrdnum + jk];
            if (ij == jk) {
                Lld_vecc[ij] = Lld[ij];
            }
        }
        Lldd_vec.attr("dim") = Dimension(reqrdnum, reqrdnum);
        limits[1] = beta_0[para_number];
        ll_final[1] = Ll[0];
        if ((max_change < epsilon) && (deriv_max < deriv_epsilon)) {
            iter_continue = false;
            limit_converged[1] = TRUE;
        }
    }
    //
    beta_p = beta_peak;  //
    beta_a = beta_peak;  //
    beta_c = beta_peak;  //
    for (int ij = 0; ij < totalnum; ij++) {
        beta_0[ij] = beta_a[ij];
    }
    Cox_Refresh_R_SIDES(reqrdnum, ntime, Rls1, Rls2, Rls3, Lls1, Lls2, Lls3, Strata_vals, model_bool);
    //
    Cox_Term_Risk_Calc(modelform, tform, term_n, totalnum, fir, dfc, term_tot, T0, Td0, Tdd0, Te, R, Rd, Rdd, Dose, nonDose, beta_0, df0, thres_step_max, step_max, TTerm, nonDose_LIN, nonDose_PLIN, nonDose_LOGLIN, RdR, RddR, nthreads, KeepConstant, verbose, model_bool, gmix_theta, gmix_term);
    //
    //  -------------------------------------------------------------------------------------------
    Cox_Side_LL_Calc(reqrdnum, ntime, tform, RiskFail,  RiskPairs, RiskPairs_Strata, totalnum, fir, R, Rd, Rdd, Rls1, Rls2, Rls3, Lls1, Lls2, Lls3, cens_weight, Strata_vals, beta_0, RdR, RddR, Ll, Lld, Lldd, nthreads, KeepConstant, ties_method, verbose, model_bool, iter_stop);
    #ifdef _OPENMP
    #pragma omp parallel for schedule(dynamic) num_threads(nthreads)
    #endif
    for (int ijk = 0; ijk < reqrdnum*(reqrdnum + 1)/2; ijk++) {
        int ij = 0;
        int jk = ijk;
        while (jk > ij) {
            ij++;
            jk -= ij;
        }
        Lldd_vec[ij * reqrdnum + jk] = Lldd[ij*reqrdnum+jk];
        Lldd_vec[jk * reqrdnum + ij] = Lldd_vec[ij * reqrdnum + jk];
        if (ij == jk) {
            Lld_vecc[ij] = Lld[ij];
        }
    }
    Lldd_vec.attr("dim") = Dimension(reqrdnum, reqrdnum);
    Lldd_mat = as<Map<MatrixXd> >(Lldd_vec);
    Lld_vec = as<Map<VectorXd> >(Lld_vecc);
    //
    step_max = step_max0;
    thres_step_max = thres_step_max0;
    //
    upper = false;
    step = -1;
    iter_continue = true;
    while ((step < maxstep) && (iter_continue)) {
        step++;
        trouble = false;
        half_check = 0;
        Log_Bound(deriv_max, Lldd_mat, Lld_vec, Lstar, qchi, Ll[0], para_number, nthreads, totalnum, reqrdnum, KeepConstant, term_tot, step, dbeta, beta_0, upper, trouble, verbose, mult);
        if (trouble) {
            Calc_Change_trouble(para_number, nthreads, totalnum, thres_step_max, lr, step_max, Ll, Lld, Lldd, dbeta, tform, thres_step_max, step_max, KeepConstant_trouble);
        }
        //
        beta_p = beta_c;  //
        beta_a = beta_c;  //
        Cox_Refresh_R_SIDES(reqrdnum, ntime, Rls1, Rls2, Rls3, Lls1, Lls2, Lls3, Strata_vals, model_bool);
        //
        for (int ijk = 0; ijk < totalnum; ijk++) {
            if (KeepConstant[ijk] == 0) {
                if ((tform[ijk] == "lin_quad_int") || (tform[ijk] == "lin_exp_int") || (tform[ijk] == "step_int") || (tform[ijk] == "lin_int")) {  //  the threshold values use different maximum deviation values
                    if (abs(dbeta[ijk]) > thres_step_max) {
                        dbeta[ijk] = thres_step_max * sign(dbeta[ijk]);
                    }
                } else {
                    if (abs(dbeta[ijk]) > step_max) {
                        dbeta[ijk] = step_max * sign(dbeta[ijk]);
                    }
                }
            } else {
                dbeta[ijk] = 0;
            }
        }
        max_change = abs(dbeta[0]);
        for (int ij = 0; ij < totalnum; ij++) {
            if (ij == para_number) {
                //  prevent the parameter estimate from crossing the optimum
                //  issue is beta_0[para_number] <= beta_best[para_number]
                if (dbeta[ij] > (beta_peak[para_number] - beta_a[ij])/lr) {
                    dbeta[ij] = (beta_peak[para_number] - beta_a[ij])/lr/2;
                }
            }
            if (abs(dbeta[ij]) > max_change) {
                max_change = abs(dbeta[ij]);
            }
            beta_0[ij] = beta_a[ij] + lr*dbeta[ij];
            beta_c[ij] = beta_0[ij];
        }
        //  ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------//
        //  The same subterm, risk, sides, and log-likelihood calculations are performed every half-step and iteration
        //  ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------//
        Cox_Term_Risk_Calc(modelform, tform, term_n, totalnum, fir, dfc, term_tot, T0, Td0, Tdd0, Te, R, Rd, Rdd, Dose, nonDose, beta_0, df0, thres_step_max, step_max, TTerm, nonDose_LIN, nonDose_PLIN, nonDose_LOGLIN, RdR, RddR, nthreads, KeepConstant, verbose, model_bool, gmix_theta, gmix_term);
        bound_val = 0;
        Cox_Pois_Log_Loop(step_max, model_bool, beta_0, beta_a, beta_c, bound_val, dbeta, df0, dfc, dint, Dose, thres_step_max, dslp, fir, gmix_term, gmix_theta, half_check, halfmax, KeepConstant, limit_hit, lr, modelform, nonDose, nonDose_LIN, nonDose_LOGLIN, nonDose_PLIN, nthreads, R, Rd, Rdd, RddR, RdR, s_weights, T0, Td0, Tdd0, Te, term_n, term_tot, tform, totalnum, TTerm, verbose);
        if (limit_hit[0]) {
            limits[0] = 0;
            ll_final[0] = 0;
            limit_converged[0] = FALSE;
            break;
        }
        Cox_Side_LL_Calc(reqrdnum, ntime, tform, RiskFail,  RiskPairs, RiskPairs_Strata, totalnum, fir, R, Rd, Rdd, Rls1, Rls2, Rls3, Lls1, Lls2, Lls3, cens_weight, Strata_vals, beta_0, RdR, RddR, Ll, Lld, Lldd, nthreads, KeepConstant, ties_method, verbose, model_bool, iter_stop);
        //
        #ifdef _OPENMP
        #pragma omp parallel for schedule(dynamic) num_threads(nthreads)
        #endif
        for (int ijk = 0; ijk < reqrdnum*(reqrdnum + 1)/2; ijk++) {
            int ij = 0;
            int jk = ijk;
            while (jk > ij) {
                ij++;
                jk -= ij;
            }
            Lldd_vec[ij * reqrdnum + jk] = Lldd[ij*reqrdnum+jk];
            Lldd_vec[jk * reqrdnum + ij] = Lldd_vec[ij * reqrdnum + jk];
            if (ij == jk) {
                Lld_vecc[ij] = Lld[ij];
            }
        }
        Lldd_vec.attr("dim") = Dimension(reqrdnum, reqrdnum);
        Lldd_mat = as<Map<MatrixXd> >(Lldd_vec);
        Lld_vec = as<Map<VectorXd> >(Lld_vecc);
        limits[0] = beta_0[para_number];
        ll_final[0] = Ll[0];
        if ((max_change < epsilon) && (deriv_max < deriv_epsilon)) {
            iter_continue = false;
            limit_converged[0] = TRUE;
        }
    }
    //
    res_list = List::create(_["Parameter_Limits"] = wrap(limits), _["Negative_Risk_Limit_Hit"] = wrap(limit_hit), _["Likelihood_Boundary"] = wrap(ll_final), _["Likelihood_Goal"] = wrap(Lstar), _["Limit_Converged"] = wrap(limit_converged), _["Status"] = "PASSED");
    //  returns a list of results
    return res_list;
}

//' Primary Cox PH likelihood bounds calcualtion function.
//'
//' \code{LogLik_Cox_PH_Omnibus_Log_Bound_Search} Performs the calls to calculation functions and log-likeihood profile bounds
//'
//' @inheritParams CPP_template
//'
//' @return List of final results: Log-likelihood of optimum, first derivative of log-likelihood, second derivative matrix, parameter list, standard deviation estimate, AIC, model information
//' @noRd
//'
//
List LogLik_Cox_PH_Omnibus_Log_Bound_Search(IntegerVector term_n, StringVector tform, Ref<VectorXd> beta_0, Ref<MatrixXd> df0, IntegerVector dfc, int fir, string modelform, double lr, NumericVector maxiters, int guesses, int halfmax, double epsilon, double step_max, double thres_step_max, double deriv_epsilon, const Ref<const MatrixXd>& df_m, NumericVector tu, int verbose, IntegerVector KeepConstant, int term_tot, string ties_method, int nthreads, NumericVector& Strata_vals, const VectorXd& cens_weight, List model_bool, const double gmix_theta, const IntegerVector gmix_term, const Ref<const MatrixXd>& Lin_Sys, const Ref<const VectorXd>& Lin_Res, double qchi, int para_number, int maxstep, double mult) {
    //
    List temp_list = List::create(_["Status"] = "TEMP");  //  used as a dummy return value for code checking
    //  Time durations are measured from this point on in microseconds
    //
    //  df0: covariate data
    //  ntime: number of event times for Cox PH
    //  totalnum: number of terms used
    //
    //  ------------------------------------------------------------------------- //  initialize
    int ntime = tu.size();
    int totalnum = term_n.size();
    int reqrdnum = totalnum - sum(KeepConstant);
    //  ------------------------------------------------------------------------- //  initialize
    if (model_bool["null"]) {
        if (verbose >= 1) {
            Rcout << "null model is not compatable with log-based bound calculation" << endl;
        }
        temp_list = List::create(_["Status"] = "FAILED_WITH_BAD_MODEL_NULL", _["LogLik"] = R_NaN);
        return temp_list;
    }
    if (model_bool["single"]) {
        if (verbose >= 1) {
            Rcout << "non-derivative model calculation is not compatable with log-based bound calculation" << endl;
        }
        temp_list = List::create(_["Status"] = "FAILED_WITH_BAD_MODEL_SINGLE", _["LogLik"] = R_NaN);
        return temp_list;
    }
    if (model_bool["gradient"]) {
        if (verbose >= 1) {
            Rcout << "gradient descent model calculation is not compatable with log-based bound calculation" << endl;
        }
        temp_list = List::create(_["Status"] = "FAILED_WITH_BAD_MODEL_GRADIENT", _["LogLik"] = R_NaN);
        return temp_list;
    }
    //
    //  cout.precision: controls the number of significant digits printed
    //
    Rcout.precision(7);  //  forces higher precision numbers printed to terminal
    double dbeta_max = 0.0;
    //  ------------------------------------------------------------------------- //  initialize
    // Map<VectorXd> beta_0(as<Map<VectorXd> >(a_n));
    VectorXd beta_max = beta_0;
    MatrixXd T0;
    MatrixXd Td0;
    MatrixXd Tdd0;
    //
    MatrixXd Te;
    MatrixXd R;
    MatrixXd Rd;
    MatrixXd Rdd;
    //
    MatrixXd Dose;
    MatrixXd nonDose;
    MatrixXd nonDose_LIN;
    MatrixXd nonDose_PLIN;
    MatrixXd nonDose_LOGLIN;
    MatrixXd TTerm;
    double dint = 0.0;  //  the amount of change used to calculate derivatives in threshold paramters
    double dslp = 0.0;
    MatrixXd RdR;
    MatrixXd RddR;
    //  ------------------------------------------------------------------------- //  initialize
    //  ---------------------------------------------
    //  To Start, needs to seperate the derivative terms
    //  ---------------------------------------------
    //
    Cox_Refresh_R_TERM(totalnum, reqrdnum, term_tot, dint, dslp, thres_step_max, step_max, df0, T0, Td0, Tdd0, Te, R, Rd, Rdd, Dose, nonDose, TTerm, nonDose_LIN, nonDose_PLIN, nonDose_LOGLIN, RdR, RddR, model_bool);
    //  ------------------------------------------------------------------------- //  initialize
    IntegerMatrix RiskFail;
    vector<vector<int> > RiskPairs(ntime);
    vector<vector<vector<int> > > RiskPairs_Strata(ntime, vector<vector<int>>(Strata_vals.size()));
    //  ------------------------------------------------------------------------- //  initialize
    if (model_bool["strata"]) {
        RiskFail = IntegerMatrix(ntime, 2*Strata_vals.size());  //  vector giving the event rows
        //  Creates matrices used to identify the event risk groups
        if (model_bool["cr"]) {
            Make_Groups_Strata_CR(ntime, df_m, RiskFail, RiskPairs_Strata, tu, nthreads, Strata_vals, cens_weight);
        } else {
            Make_Groups_Strata(ntime, df_m, RiskFail, RiskPairs_Strata, tu, nthreads, Strata_vals);
        }
    } else {
        RiskFail = IntegerMatrix(ntime, 2);  //  vector giving the event rows
        //  Creates matrices used to identify the event risk groups
        if (model_bool["cr"]) {
            Make_Groups_CR(ntime, df_m, RiskFail, RiskPairs, tu, cens_weight, nthreads);
        } else {
            Make_Groups(ntime, df_m, RiskFail, RiskPairs, tu, nthreads);
        }
    }
    //  ------------------------------------------------------------------------- //  initialize
    MatrixXd Rls1;
    MatrixXd Lls1;
    MatrixXd Rls2;
    MatrixXd Rls3;
    MatrixXd Lls2;
    MatrixXd Lls3;
    vector<double> Ll(reqrdnum, 0.0);  //  log-likelihood values
    vector<double> Lld(reqrdnum, 0.0);  //  log-likelihood derivative values
    vector<double> Lldd(pow(reqrdnum, 2), 0.0);  //  the second derivative matrix has room for every combination, but only the lower triangle is calculated initially
    //  ------------------------------------------------------------------------- //  initialize
    Cox_Refresh_R_SIDES(reqrdnum, ntime, Rls1, Rls2, Rls3, Lls1, Lls2, Lls3, Strata_vals, model_bool);
    //  the log-likelihood is calculated in parallel over the risk groups
    vector <double> Ll_comp(2, Ll[0]);  //  vector to compare values
    double step_max0 = step_max;
    double thres_step_max0 = thres_step_max;
    //
    vector<double> dbeta(totalnum, 0.0);
    vector<double> dbeta_start(totalnum, 0.0);
    //
    //  --------------------------
    //  always starts from initial guess
    //  --------------------------
    vector<double> beta_c(totalnum, 0.0);
    vector<double> beta_a(totalnum, 0.0);
    vector<double> beta_best(totalnum, 0.0);
    vector<double> beta_p(totalnum, 0.0);
    vector<double> beta_peak(totalnum, 0.0);
    VectorXd::Map(&beta_p[0], beta_0.size()) = beta_0;  //  stores previous parameters
    VectorXd::Map(&beta_c[0], beta_0.size()) = beta_0;  //  stores current parameters
    VectorXd::Map(&beta_a[0], beta_0.size()) = beta_0;  //  stores a refrence value for parameters
    VectorXd::Map(&beta_best[0], beta_0.size()) = beta_0;  //  stores the best parameters
    VectorXd::Map(&beta_peak[0], beta_0.size()) = beta_0;  //  stores the peak parameters
    int iter_stop  = 0;  //  tracks if the iterations should be stopped for convergence
    vector<double> beta_abs_best(totalnum, 0.0);
    //
    Cox_Refresh_R_SIDES(reqrdnum, ntime, Rls1, Rls2, Rls3, Lls1, Lls2, Lls3, Strata_vals, model_bool);
    fill(Ll.begin(), Ll.end(), 0.0);
    fill(Lld.begin(), Lld.end(), 0.0);
    fill(Lldd.begin(), Lldd.end(), 0.0);
    beta_p = beta_peak;  //
    beta_a = beta_peak;  //
    beta_c = beta_peak;  //
    step_max = step_max0;
    thres_step_max = thres_step_max0;
    //
    for (int i = 0; i < beta_0.size(); i++) {
        beta_0[i] = beta_peak[i];
    }
    Cox_Term_Risk_Calc(modelform, tform, term_n, totalnum, fir, dfc, term_tot, T0, Td0, Tdd0, Te, R, Rd, Rdd, Dose, nonDose, beta_0, df0, thres_step_max, step_max, TTerm, nonDose_LIN, nonDose_PLIN, nonDose_LOGLIN, RdR, RddR, nthreads, KeepConstant, verbose, model_bool, gmix_theta, gmix_term);
    Cox_Side_LL_Calc(reqrdnum, ntime, tform, RiskFail,  RiskPairs, RiskPairs_Strata, totalnum, fir, R, Rd, Rdd, Rls1, Rls2, Rls3, Lls1, Lls2, Lls3, cens_weight, Strata_vals, beta_0, RdR, RddR, Ll, Lld, Lldd, nthreads, KeepConstant, ties_method, verbose, model_bool, iter_stop);
    //
    NumericVector Lldd_vec(reqrdnum * reqrdnum);
    NumericVector Lld_vecc(reqrdnum);
    #ifdef _OPENMP
    #pragma omp parallel for schedule(dynamic) num_threads(nthreads)
    #endif
    for (int ijk = 0; ijk < reqrdnum*(reqrdnum + 1)/2; ijk++) {
        int ij = 0;
        int jk = ijk;
        while (jk > ij) {
            ij++;
            jk -= ij;
        }
        Lldd_vec[ij * reqrdnum + jk] = Lldd[ij*reqrdnum+jk];
        Lldd_vec[jk * reqrdnum + ij] = Lldd_vec[ij * reqrdnum + jk];
        if (ij == jk) {
            Lld_vecc[ij] = Lld[ij];
        }
    }
    Lldd_vec.attr("dim") = Dimension(reqrdnum, reqrdnum);
    Map<MatrixXd> Lldd_mat(as<Map<MatrixXd> >(Lldd_vec));
    Map<VectorXd> Lld_vec(as<Map<VectorXd> >(Lld_vecc));
    double Lstar = Ll[0]-qchi;
    double Lmax = Ll[0];
    bool upper = true;
    int half_check = 0;
    bool trouble = false;
    double deriv_max = 100;
    //
    //
    IntegerVector KeepConstant_trouble(totalnum);
    for (int i = 0;  i < totalnum; i++) {
        KeepConstant_trouble[i] = KeepConstant[i];
    }
    KeepConstant_trouble[para_number] = 1;
    //
    vector<double> limits(2, 0.0);
    vector<bool>   limit_hit(2, FALSE);
    vector<bool>   limit_converged(2, FALSE);
    vector<double> ll_final(2, 0.0);
    List res_list;
    if (verbose >= 4) {
        Rcout << "C++ Note: STARTING BOUNDS" << endl;
    }
    //  //
    if (verbose >= 4) {
         Rcout << "C++ Note: STARTING Upper Bound" << endl;
    }
    upper = true;
    //  Now define the list of points to check
    trouble = false;
    Log_Bound(deriv_max, Lldd_mat, Lld_vec, Lstar, qchi, Lmax, para_number, nthreads, totalnum, reqrdnum, KeepConstant, term_tot, 0, dbeta_start, beta_0, upper, trouble, verbose, mult);
    NumericMatrix a_ns(guesses, totalnum);
    NumericVector a_n = a_ns.row(0);
    //  now the dbeta holds the range to check
    //  note that currently the guesses are not bounded
    for (int i = 0; i < guesses; i++) {
        for (int j = 0; j < totalnum; j++) {
            //  use dbeta to assign a_n guess i, parameter j
            //  assign such that i=guesses - 1 gives mult*dbeta
            a_ns(i, j) = beta_0[j] + dbeta_start[j] * (i + 1)/(guesses);
        }
    }
    if (verbose >= 4) {
        for (int i = 0; i < guesses; i++) {
            Rcout << "C++ Note: Initial guess " << i << ": ";
            for (int j = 0; j < totalnum; j++) {
                Rcout << a_ns(i, j) << " ";
            }
            Rcout << " " << endl;
        }
    }
    //  now we have the points to test
    double halves = 0;  //  number of half-steps taken
    int ind0 = fir;  //  used for validations
    int iteration = 0;  //  iteration number
    int maxiter = 0;
    //
    iter_stop  = 0;  //  tracks if the iterations should be stopped for convergence
    //
    NumericMatrix beta_fin(a_ns.rows(), a_ns.cols());
    NumericVector LL_fin(a_ns.rows());
    //
    double Ll_iter_best = 10;
    //  Variables that are used for the risk check function shared across cox, poisson, and log bound functions
    double dev = 0.0;
    MatrixXd dev_temp = MatrixXd::Zero(1, 1);
    MatrixXd PyrC = MatrixXd::Zero(1, 1);
    ///
    for (int guess = 0; guess <guesses; guess++) {
        Cox_Refresh_R_SIDES(reqrdnum, ntime, Rls1, Rls2, Rls3, Lls1, Lls2, Lls3, Strata_vals, model_bool);
        fill(Ll.begin(), Ll.end(), 0.0);
        fill(Lld.begin(), Lld.end(), 0.0);
        fill(Lldd.begin(), Lldd.end(), 0.0);
        beta_p = beta_peak;  //
        beta_a = beta_peak;  //
        beta_c = beta_peak;  //
        step_max = step_max0;
        thres_step_max = thres_step_max0;
        iter_stop = 0;
        halves = 0;
        iteration = 0;
        halves = 0;  //  number of half-steps taken
        ind0 = fir;  //  used for validations
        iteration = 0;  //  iteration number
        Ll_iter_best = 10;
        //
        iter_stop  = 0;  //  tracks if the iterations should be stopped for convergence
        //
        maxiter = maxiters[0];
        a_n = a_ns.row(guess);
        for (int i = 0; i < beta_0.size(); i++) {
            beta_0[i] = a_n[i];
        }
        Cox_Term_Risk_Calc(modelform, tform, term_n, totalnum, fir, dfc, term_tot, T0, Td0, Tdd0, Te, R, Rd, Rdd, Dose, nonDose, beta_0, df0, thres_step_max, step_max, TTerm, nonDose_LIN, nonDose_PLIN, nonDose_LOGLIN, RdR, RddR, nthreads, KeepConstant, verbose, model_bool, gmix_theta, gmix_term);
        if ((R.minCoeff() <= 0) || (R.hasNaN())) {
            iter_stop = 1;
            Ll[0] = 404;  //  note to consider this point too far
        } else {
            Cox_Side_LL_Calc(reqrdnum, ntime, tform, RiskFail,  RiskPairs, RiskPairs_Strata, totalnum, fir, R, Rd, Rdd, Rls1, Rls2, Rls3, Lls1, Lls2, Lls3, cens_weight, Strata_vals, beta_0, RdR, RddR, Ll, Lld, Lldd, nthreads, KeepConstant, ties_method, verbose, model_bool, iter_stop);
            for (int i = 0; i < beta_0.size(); i++) {
                beta_c[i] = beta_0[i];
            }
        }
        while ((iteration < maxiter) && (iter_stop == 0)) {
            iteration++;
            beta_p = beta_c;  //
            beta_a = beta_c;  //
            beta_best = beta_c;  //
            //
            //  calculates the initial change in parameter
            Calc_Change_trouble(para_number, nthreads, totalnum, thres_step_max, lr, step_max, Ll, Lld, Lldd, dbeta, tform, thres_step_max, step_max, KeepConstant_trouble);
            Intercept_Bound(nthreads, totalnum, beta_0, dbeta, dfc, df0, KeepConstant, tform);
            //
            if ((Ll_iter_best > 0) || (Ll_iter_best < Ll[ind0])) {
                Ll_iter_best = Ll[ind0];
                beta_abs_best = beta_c;
            }
            //
            halves = 0;
            while ((Ll[ind0] <= Ll_iter_best) && (halves < halfmax)) {  //  repeats until half-steps maxed or an improvement
                for (int ij = 0; ij < totalnum; ij++) {
                    beta_0[ij] = beta_a[ij] + dbeta[ij];
                    beta_c[ij] = beta_0[ij];
                }
                //  -----------------------------------------------------------------------------------------------------------------------------------------//
                //  The same subterm, risk, sides, and log-likelihood calculations are performed every half-step and iteration
                //  -----------------------------------------------------------------------------------------------------------------------------------------//
                Cox_Term_Risk_Calc(modelform, tform, term_n, totalnum, fir, dfc, term_tot, T0, Td0, Tdd0, Te, R, Rd, Rdd, Dose, nonDose, beta_0, df0, thres_step_max, step_max, TTerm, nonDose_LIN, nonDose_PLIN, nonDose_LOGLIN, RdR, RddR, nthreads, KeepConstant, verbose, model_bool, gmix_theta, gmix_term);
                Cox_Pois_Check_Continue(model_bool, beta_0, beta_best, beta_c, cens_weight, dbeta, dev, dev_temp, fir, halfmax, halves, ind0, iter_stop, KeepConstant, Ll, Ll_iter_best, Lld, Lldd, Lls1, Lls2, Lls3, Lstar, nthreads, ntime, PyrC, R, Rd, Rdd, RddR, RdR, reqrdnum, tform, RiskFail,  RiskPairs, RiskPairs_Strata, Rls1, Rls2, Rls3, Strata_vals, term_n, ties_method, totalnum, TTerm, verbose);
            }
            if (beta_best != beta_c) {  //  if the risk matrices aren't the optimal values, then they must be recalculated
                //  If it goes through every half step without improvement, then the maximum change needs to be decreased
                step_max = step_max*pow(0.5, halfmax);  //  reduces the step sizes
                thres_step_max = thres_step_max*pow(0.5, halfmax);
                //
                beta_p = beta_best;  //
                beta_a = beta_best;  //
                beta_c = beta_best;  //
                for (int ij = 0; ij < totalnum; ij++) {
                    beta_0[ij] = beta_best[ij];
                }
                Cox_Term_Risk_Calc(modelform, tform, term_n, totalnum, fir, dfc, term_tot, T0, Td0, Tdd0, Te, R, Rd, Rdd, Dose, nonDose, beta_0, df0, thres_step_max, step_max, TTerm, nonDose_LIN, nonDose_PLIN, nonDose_LOGLIN, RdR, RddR, nthreads, KeepConstant, verbose, model_bool, gmix_theta, gmix_term);
            }
            dbeta_max = abs(dbeta[0]);
            for (int ij = 1; ij < totalnum; ij++) {
                if (abs(dbeta[ij]) > dbeta_max) {
                    dbeta_max = abs(dbeta[ij]);
                }
            }
            if (step_max < epsilon) {  //  if the maximum change is too low, then it ends
                iter_stop = 1;
            }
            if (dbeta_max < epsilon) {  //  if the maximum change is too low, then it ends
                iter_stop = 1;
            }
        }
        //  -----------------------------------------------
        //  Performing Full Calculation to get full second derivative matrix
        //  -----------------------------------------------
        fill(Ll.begin(), Ll.end(), 0.0);
        if (!model_bool["single"]) {
            fill(Lld.begin(), Lld.end(), 0.0);
            fill(Lldd.begin(), Lldd.end(), 0.0);
        }
        Cox_Side_LL_Calc(reqrdnum, ntime, tform, RiskFail,  RiskPairs, RiskPairs_Strata, totalnum, fir, R, Rd, Rdd, Rls1, Rls2, Rls3, Lls1, Lls2, Lls3, cens_weight, Strata_vals, beta_0, RdR, RddR, Ll, Lld, Lldd, nthreads, KeepConstant, ties_method, verbose, model_bool, iter_stop);
        //
        a_n = beta_0;
        beta_fin(guess, _) = a_n;
        LL_fin[guess] = Ll[0];
        //
    }
    //
    Cox_Refresh_R_SIDES(reqrdnum, ntime, Rls1, Rls2, Rls3, Lls1, Lls2, Lls3, Strata_vals, model_bool);
    fill(Ll.begin(), Ll.end(), 0.0);
    if (!model_bool["single"]) {
        fill(Lld.begin(), Lld.end(), 0.0);
        fill(Lldd.begin(), Lldd.end(), 0.0);
    }
    beta_p = beta_best;  //
    beta_a = beta_best;  //
    beta_c = beta_best;  //
    step_max = step_max0;
    thres_step_max = thres_step_max0;
    iter_stop = 0;
    halves = 0;
    iteration = 0;
    halves = 0;  //  number of half-steps taken
    ind0 = fir;  //  used for validations
    iteration = 0;  //  iteration number
    //
    iter_stop  = 0;  //  tracks if the iterations should be stopped for convergence
    //
    //  next we need to figure out what point to start at
    int best_guess = - 1;
    for (int guess = guesses - 1; guess> - 1; guess--) {
        //  we want to find the closest point at which the loglikelihood is below lstar
        if (LL_fin[guess] < Lstar) {
            best_guess = guess;
        }
    }
    if (verbose >= 3) {
        Rcout << "C++ Note: Upper Guess Results" << endl;
        Rcout << "Guess number, parameter values, Log-Likelihood change" << endl;
        NumericVector beta_temp;
        for (int i = 0;  i < guesses; i++) {
            beta_temp = wrap(beta_fin.row(i));
            if (i == best_guess) {
                Rcout << i << ", " << beta_temp << ", " << LL_fin[i] - Lstar << "<-- Best Guess" << endl;
            } else {
                Rcout << i << ", " << beta_temp << ", " << LL_fin[i] - Lstar << endl;
            }
        }
    }
    NumericVector beta_temp;
    NumericVector beta_temp0;
    if (best_guess == 0) {
        beta_temp = wrap(beta_fin.row(best_guess));  //  the first point was closest, no lower bound
    } else if (best_guess == - 1) {
        beta_temp = wrap(beta_fin.row(guesses - 1));  //  the last point was closest, no upper bound
    } else {
        beta_temp = wrap(beta_fin.row(best_guess));
        beta_temp0 = wrap(beta_fin.row(best_guess - 1));
        for (int i = 0; i < beta_temp.size(); i++) {
            beta_temp[i] = beta_temp[i] + beta_temp0[i];
            beta_temp[i] = beta_temp[i]/2;
        }
    }
    if (verbose >= 4) {
        Rcout << "C++ Note: Initial Guess: ";
        for (int i = 0; i < beta_0.size(); i++) {
            Rcout << beta_temp[i] << " ";
        }
        Rcout << " " << endl;
    }
    for (int i = 0; i < beta_0.size(); i++) {
        a_n[i] = beta_temp[i];
        beta_0[i] = a_n[i];
        beta_c[i] = beta_0[i];
    }
    //
    Cox_Refresh_R_SIDES(reqrdnum, ntime, Rls1, Rls2, Rls3, Lls1, Lls2, Lls3, Strata_vals, model_bool);
    //  Calculates the subterm and term values
    Cox_Term_Risk_Calc(modelform, tform, term_n, totalnum, fir, dfc, term_tot, T0, Td0, Tdd0, Te, R, Rd, Rdd, Dose, nonDose, beta_0, df0, thres_step_max, step_max, TTerm, nonDose_LIN, nonDose_PLIN, nonDose_LOGLIN, RdR, RddR, nthreads, KeepConstant, verbose, model_bool, gmix_theta, gmix_term);
    //
    //  -------------------------------------------------------------------------------------------
    //  Calculates the side sum terms used
    Cox_Side_LL_Calc(reqrdnum, ntime, tform, RiskFail,  RiskPairs, RiskPairs_Strata, totalnum, fir, R, Rd, Rdd, Rls1, Rls2, Rls3, Lls1, Lls2, Lls3, cens_weight, Strata_vals, beta_0, RdR, RddR, Ll, Lld, Lldd, nthreads, KeepConstant, ties_method, verbose, model_bool, iter_stop);
    //
    #ifdef _OPENMP
    #pragma omp parallel for schedule(dynamic) num_threads(nthreads)
    #endif
    for (int ijk = 0; ijk < reqrdnum*(reqrdnum + 1)/2; ijk++) {
        int ij = 0;
        int jk = ijk;
        while (jk > ij) {
            ij++;
            jk -= ij;
        }
        Lldd_vec[ij * reqrdnum + jk] = Lldd[ij*reqrdnum+jk];
        Lldd_vec[jk * reqrdnum + ij] = Lldd_vec[ij * reqrdnum + jk];
        if (ij == jk) {
            Lld_vecc[ij] = Lld[ij];
        }
    }
    Lldd_vec.attr("dim") = Dimension(reqrdnum, reqrdnum);
    Lldd_mat = as<Map<MatrixXd> >(Lldd_vec);
    Lld_vec = as<Map<VectorXd> >(Lld_vecc);
    int step = 1;
    bool iter_continue = true;
    double max_change = 100;
    ///
    //  variables added for log loop code
    ///
    VectorXd s_weights(1);
    int bound_val = 1;
    //
    while ((step < maxstep) && (iter_continue)) {
        Lldd_mat = as<Map<MatrixXd> >(Lldd_vec);
        Lld_vec = as<Map<VectorXd> >(Lld_vecc);
        step++;
        trouble = false;
        Log_Bound(deriv_max, Lldd_mat, Lld_vec, Lstar, qchi, Ll[0], para_number, nthreads, totalnum, reqrdnum, KeepConstant, term_tot, step, dbeta, beta_0, upper, trouble, verbose, mult);
        if (trouble) {
            Calc_Change_trouble(para_number, nthreads, totalnum, thres_step_max, lr, step_max, Ll, Lld, Lldd, dbeta, tform, thres_step_max, step_max, KeepConstant_trouble);
        }
        //
        beta_p = beta_c;  //
        beta_a = beta_c;  //
        Cox_Refresh_R_SIDES(reqrdnum, ntime, Rls1, Rls2, Rls3, Lls1, Lls2, Lls3, Strata_vals, model_bool);
        //
        for (int ijk = 0; ijk < totalnum; ijk++) {
            if (KeepConstant[ijk] == 0) {
                //
                if ((tform[ijk] == "lin_quad_int") || (tform[ijk] == "lin_exp_int") || (tform[ijk] == "step_int") || (tform[ijk] == "lin_int")) {  //  the threshold values use different maximum deviation values
                    if (abs(dbeta[ijk]) > thres_step_max) {
                        dbeta[ijk] = thres_step_max * sign(dbeta[ijk]);
                    }
                } else {
                    if (abs(dbeta[ijk]) > step_max) {
                        dbeta[ijk] = step_max * sign(dbeta[ijk]);
                    }
                }
            } else {
                dbeta[ijk] = 0;
            }
        }
        max_change = abs(dbeta[0]);
        //
        for (int ij = 0; ij < totalnum; ij++) {
            if (ij == para_number) {
                //  prevent the parameter estimate from crossing the optimum
                //  issue is beta_0[para_number] <= beta_best[para_number]
                if (dbeta[ij] < (beta_peak[para_number] - beta_a[ij])/lr) {
                   dbeta[ij] = (beta_peak[para_number] - beta_a[ij])/lr/2;
                }
            }
            if (abs(dbeta[ij]) > max_change) {
                max_change = abs(dbeta[ij]);
            }
            beta_0[ij] = beta_a[ij] + lr*dbeta[ij];
            beta_c[ij] = beta_0[ij];
        }
        //  ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------//
        //  The same subterm, risk, sides, and log-likelihood calculations are performed every half-step and iteration
        //  ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------//
        Cox_Term_Risk_Calc(modelform, tform, term_n, totalnum, fir, dfc, term_tot, T0, Td0, Tdd0, Te, R, Rd, Rdd, Dose, nonDose, beta_0, df0, thres_step_max, step_max, TTerm, nonDose_LIN, nonDose_PLIN, nonDose_LOGLIN, RdR, RddR, nthreads, KeepConstant, verbose, model_bool, gmix_theta, gmix_term);
        bound_val = 1;
        Cox_Pois_Log_Loop(step_max, model_bool, beta_0, beta_a, beta_c, bound_val, dbeta, df0, dfc, dint, Dose, thres_step_max, dslp, fir, gmix_term, gmix_theta, half_check, halfmax, KeepConstant, limit_hit, lr, modelform, nonDose, nonDose_LIN, nonDose_LOGLIN, nonDose_PLIN, nthreads, R, Rd, Rdd, RddR, RdR, s_weights, T0, Td0, Tdd0, Te, term_n, term_tot, tform, totalnum, TTerm, verbose);
        if (limit_hit[1]) {
            limits[1] = 0;
            ll_final[1] = 0;
            limit_converged[1] = FALSE;
            break;
        }
        Cox_Side_LL_Calc(reqrdnum, ntime, tform, RiskFail,  RiskPairs, RiskPairs_Strata, totalnum, fir, R, Rd, Rdd, Rls1, Rls2, Rls3, Lls1, Lls2, Lls3, cens_weight, Strata_vals, beta_0, RdR, RddR, Ll, Lld, Lldd, nthreads, KeepConstant, ties_method, verbose, model_bool, iter_stop);
        //
        #ifdef _OPENMP
        #pragma omp parallel for schedule(dynamic) num_threads(nthreads)
        #endif
        for (int ijk = 0; ijk < reqrdnum*(reqrdnum + 1)/2; ijk++) {
            int ij = 0;
            int jk = ijk;
            while (jk > ij) {
                ij++;
                jk -= ij;
            }
            Lldd_vec[ij * reqrdnum + jk] = Lldd[ij*reqrdnum+jk];
            Lldd_vec[jk * reqrdnum + ij] = Lldd_vec[ij * reqrdnum + jk];
            if (ij == jk) {
                Lld_vecc[ij] = Lld[ij];
            }
        }
        Lldd_vec.attr("dim") = Dimension(reqrdnum, reqrdnum);
        limits[1] = beta_0[para_number];
        ll_final[1] = Ll[0];
        if ((max_change < epsilon) && (deriv_max < deriv_epsilon)) {
           iter_continue = false;
           limit_converged[1] = TRUE;
        }
    }
    //  Now refresh matrices back to the maximum point
    beta_p = beta_peak;  //
    beta_a = beta_peak;  //
    beta_c = beta_peak;  //
    step_max = step_max0;
    thres_step_max = thres_step_max0;
    for (int ij = 0; ij < totalnum; ij++) {
        beta_0[ij] = beta_a[ij];
    }
    Cox_Refresh_R_SIDES(reqrdnum, ntime, Rls1, Rls2, Rls3, Lls1, Lls2, Lls3, Strata_vals, model_bool);
    Cox_Term_Risk_Calc(modelform, tform, term_n, totalnum, fir, dfc, term_tot, T0, Td0, Tdd0, Te, R, Rd, Rdd, Dose, nonDose, beta_0, df0, thres_step_max, step_max, TTerm, nonDose_LIN, nonDose_PLIN, nonDose_LOGLIN, RdR, RddR, nthreads, KeepConstant, verbose, model_bool, gmix_theta, gmix_term);
    Cox_Side_LL_Calc(reqrdnum, ntime, tform, RiskFail,  RiskPairs, RiskPairs_Strata, totalnum, fir, R, Rd, Rdd, Rls1, Rls2, Rls3, Lls1, Lls2, Lls3, cens_weight, Strata_vals, beta_0, RdR, RddR, Ll, Lld, Lldd, nthreads, KeepConstant, ties_method, verbose, model_bool, iter_stop);
    upper = false;
    #ifdef _OPENMP
    #pragma omp parallel for schedule(dynamic) num_threads(nthreads)
    #endif
    for (int ijk = 0; ijk < reqrdnum*(reqrdnum + 1)/2; ijk++) {
        int ij = 0;
        int jk = ijk;
        while (jk > ij) {
            ij++;
            jk -= ij;
        }
        Lldd_vec[ij * reqrdnum + jk] = Lldd[ij*reqrdnum+jk];
        Lldd_vec[jk * reqrdnum + ij] = Lldd_vec[ij * reqrdnum + jk];
        if (ij == jk) {
            Lld_vecc[ij] = Lld[ij];
        }
    }
    Lldd_vec.attr("dim") = Dimension(reqrdnum, reqrdnum);
    Lldd_mat = as<Map<MatrixXd> >(Lldd_vec);
    Lld_vec = as<Map<VectorXd> >(Lld_vecc);
    //
    trouble = false;
    Log_Bound(deriv_max, Lldd_mat, Lld_vec, Lstar, qchi, Lmax, para_number, nthreads, totalnum, reqrdnum, KeepConstant, term_tot, 0, dbeta_start, beta_0, upper, trouble, verbose, mult);
    //  now the dbeta holds the range to check
    //  note that currently the guesses are not bounded
    for (int i = 0; i < guesses; i++) {
        for (int j = 0; j < totalnum; j++) {
            //  use dbeta to assign a_n guess i, parameter j
            //  assign such that i=guesses - 1 gives mult*dbeta
            a_ns(i, j) = beta_0[j] + dbeta_start[j] * (i + 1)/(guesses);
        }
    }
    if (verbose >= 4) {
        for (int i = 0; i < guesses; i++) {
            Rcout << "C++ Note: Initial guess " << i << ": ";
            for (int j = 0; j < totalnum; j++) {
                Rcout << a_ns(i, j) << " ";
            }
            Rcout << " " << endl;
        }
    }
    //  now we have the points to test
    halves = 0;  //  number of half-steps taken
    ind0 = fir;  //  used for validations
    iteration = 0;  //  iteration number
    maxiter = 0;
    //
    iter_stop  = 0;  //  tracks if the iterations should be stopped for convergence
    //
    Ll_iter_best = 10;
    for (int guess = 0; guess <guesses; guess++) {
        Cox_Refresh_R_SIDES(reqrdnum, ntime, Rls1, Rls2, Rls3, Lls1, Lls2, Lls3, Strata_vals, model_bool);
        fill(Ll.begin(), Ll.end(), 0.0);
        fill(Lld.begin(), Lld.end(), 0.0);
        fill(Lldd.begin(), Lldd.end(), 0.0);
        beta_p = beta_peak;  //
        beta_a = beta_peak;  //
        beta_c = beta_peak;  //
        step_max = step_max0;
        thres_step_max = thres_step_max0;
        iter_stop = 0;
        halves = 0;
        iteration = 0;
        Ll_iter_best = 10;
        halves = 0;  //  number of half-steps taken
        ind0 = fir;  //  used for validations
        iteration = 0;  //  iteration number
        iter_stop  = 0;  //  tracks if the iterations should be stopped for convergence
        maxiter = maxiters[0];
        a_n = a_ns.row(guess);
        for (int i = 0; i < beta_0.size(); i++) {
            beta_0[i] = a_n[i];
        }
        Cox_Term_Risk_Calc(modelform, tform, term_n, totalnum, fir, dfc, term_tot, T0, Td0, Tdd0, Te, R, Rd, Rdd, Dose, nonDose, beta_0, df0, thres_step_max, step_max, TTerm, nonDose_LIN, nonDose_PLIN, nonDose_LOGLIN, RdR, RddR, nthreads, KeepConstant, verbose, model_bool, gmix_theta, gmix_term);
        if ((R.minCoeff() <= 0) || (R.hasNaN())) {
            iter_stop = 1;
            Ll[0] = 404;  //  note to consider this point too far
        } else {
            Cox_Side_LL_Calc(reqrdnum, ntime, tform, RiskFail,  RiskPairs, RiskPairs_Strata, totalnum, fir, R, Rd, Rdd, Rls1, Rls2, Rls3, Lls1, Lls2, Lls3, cens_weight, Strata_vals, beta_0, RdR, RddR, Ll, Lld, Lldd, nthreads, KeepConstant, ties_method, verbose, model_bool, iter_stop);
            for (int i = 0; i < beta_0.size(); i++) {
                beta_c[i] = beta_0[i];
            }
        }
        while ((iteration < maxiter) && (iter_stop == 0)) {
            iteration++;
            beta_p = beta_c;  //
            beta_a = beta_c;  //
            beta_best = beta_c;  //
            //  calculates the initial change in parameter
            Calc_Change_trouble(para_number, nthreads, totalnum, thres_step_max, lr, step_max, Ll, Lld, Lldd, dbeta, tform, thres_step_max, step_max, KeepConstant_trouble);
            Intercept_Bound(nthreads, totalnum, beta_0, dbeta, dfc, df0, KeepConstant, tform);
            if ((Ll_iter_best > 0) || (Ll_iter_best < Ll[ind0])) {
                Ll_iter_best = Ll[ind0];
                beta_abs_best = beta_c;
            }
            halves = 0;
            while ((Ll[ind0] <= Ll_iter_best) && (halves < halfmax)) {  //  repeats until half-steps maxed or an improvement
                for (int ij = 0; ij < totalnum; ij++) {
                    beta_0[ij] = beta_a[ij] + dbeta[ij];
                    beta_c[ij] = beta_0[ij];
                }
                //  ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------//
                //  The same subterm, risk, sides, and log-likelihood calculations are performed every half-step and iteration
                //  ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------//
                Cox_Term_Risk_Calc(modelform, tform, term_n, totalnum, fir, dfc, term_tot, T0, Td0, Tdd0, Te, R, Rd, Rdd, Dose, nonDose, beta_0, df0, thres_step_max, step_max, TTerm, nonDose_LIN, nonDose_PLIN, nonDose_LOGLIN, RdR, RddR, nthreads, KeepConstant, verbose, model_bool, gmix_theta, gmix_term);
                Cox_Pois_Check_Continue(model_bool, beta_0, beta_best, beta_c, cens_weight, dbeta, dev, dev_temp, fir, halfmax, halves, ind0, iter_stop, KeepConstant, Ll, Ll_iter_best, Lld, Lldd, Lls1, Lls2, Lls3, Lstar, nthreads, ntime, PyrC, R, Rd, Rdd, RddR, RdR, reqrdnum, tform, RiskFail,  RiskPairs, RiskPairs_Strata, Rls1, Rls2, Rls3, Strata_vals, term_n, ties_method, totalnum, TTerm, verbose);
            }
            if (beta_best != beta_c) {  //  if the risk matrices aren't the optimal values, then they must be recalculated
                //  If it goes through every half step without improvement, then the maximum change needs to be decreased
                step_max = step_max*pow(0.5, halfmax);  //  reduces the step sizes
                thres_step_max = thres_step_max*pow(0.5, halfmax);
                //
                beta_p = beta_best;  //
                beta_a = beta_best;  //
                beta_c = beta_best;  //
                for (int ij = 0; ij < totalnum; ij++) {
                    beta_0[ij] = beta_best[ij];
                }
                Cox_Term_Risk_Calc(modelform, tform, term_n, totalnum, fir, dfc, term_tot, T0, Td0, Tdd0, Te, R, Rd, Rdd, Dose, nonDose, beta_0, df0, thres_step_max, step_max, TTerm, nonDose_LIN, nonDose_PLIN, nonDose_LOGLIN, RdR, RddR, nthreads, KeepConstant, verbose, model_bool, gmix_theta, gmix_term);
            }
            dbeta_max = abs(dbeta[0]);
            for (int ij = 1; ij < totalnum; ij++) {
                if (abs(dbeta[ij]) > dbeta_max) {
                    dbeta_max = abs(dbeta[ij]);
                }
            }
            if (step_max < epsilon) {  //  if the maximum change is too low, then it ends
                iter_stop = 1;
            }
            if (dbeta_max < epsilon) {  //  if the maximum change is too low, then it ends
                iter_stop = 1;
            }
        }
        //  -----------------------------------------------
        //  Performing Full Calculation to get full second derivative matrix
        //  -----------------------------------------------
        fill(Ll.begin(), Ll.end(), 0.0);
        if (!model_bool["single"]) {
            fill(Lld.begin(), Lld.end(), 0.0);
            fill(Lldd.begin(), Lldd.end(), 0.0);
        }
        Cox_Side_LL_Calc(reqrdnum, ntime, tform, RiskFail,  RiskPairs, RiskPairs_Strata, totalnum, fir, R, Rd, Rdd, Rls1, Rls2, Rls3, Lls1, Lls2, Lls3, cens_weight, Strata_vals, beta_0, RdR, RddR, Ll, Lld, Lldd, nthreads, KeepConstant, ties_method, verbose, model_bool, iter_stop);
        //
        a_n = beta_0;
        beta_fin(guess, _) = a_n;
        LL_fin[guess] = Ll[0];
        //
    }
    //
    Cox_Refresh_R_SIDES(reqrdnum, ntime, Rls1, Rls2, Rls3, Lls1, Lls2, Lls3, Strata_vals, model_bool);
    fill(Ll.begin(), Ll.end(), 0.0);
    if (!model_bool["single"]) {
        fill(Lld.begin(), Lld.end(), 0.0);
        fill(Lldd.begin(), Lldd.end(), 0.0);
    }
    beta_p = beta_best;  //
    beta_a = beta_best;  //
    beta_c = beta_best;  //
    step_max = step_max0;
    thres_step_max = thres_step_max0;
    iter_stop = 0;
    halves = 0;
    iteration = 0;
    halves = 0;  //  number of half-steps taken
    ind0 = fir;  //  used for validations
    iteration = 0;  //  iteration number
    iter_stop  = 0;  //  tracks if the iterations should be stopped for convergence
    //  next we need to figure out what point to start at
    best_guess = - 1;
    for (int guess = guesses - 1; guess> - 1; guess--) {
        //  we want to find the closest point at which the loglikelihood is below lstar
        if (LL_fin[guess] < Lstar) {
            best_guess = guess;
        }
    }
    if (verbose >= 3) {
        Rcout << "C++ Note: Lower Guess Results" << endl;
        Rcout << "Guess number, parameter values, Log-Likelihood change" << endl;
        NumericVector beta_temp;
        for (int i = 0;  i < guesses; i++) {
            beta_temp = wrap(beta_fin.row(i));
            if (i == best_guess) {
                Rcout << i << ", " << beta_temp << ", " << LL_fin[i] - Lstar << "<-- Best Guess" << endl;
            } else {
                Rcout << i << ", " << beta_temp << ", " << LL_fin[i] - Lstar << endl;
            }
        }
    }
    if (best_guess == 0) {
        beta_temp = wrap(beta_fin.row(best_guess));  //  the first point was closest, no lower bound
    } else if (best_guess == - 1) {
        beta_temp = wrap(beta_fin.row(guesses - 1));  //  the last point was closest, no upper bound
    } else {
        beta_temp = wrap(beta_fin.row(best_guess));
        beta_temp0 = wrap(beta_fin.row(best_guess - 1));
        for (int i = 0; i < beta_temp.size(); i++) {
            beta_temp[i] = beta_temp[i] + beta_temp0[i];
            beta_temp[i] = beta_temp[i]/2;
        }
    }
    if (verbose >= 4) {
        Rcout << "C++ Note: Initial Guess: ";
        for (int i = 0; i < beta_0.size(); i++) {
            Rcout << beta_temp[i] << " ";
        }
        Rcout << " " << endl;
    }
    for (int i = 0; i < beta_0.size(); i++) {
        a_n[i] = beta_temp[i];
        beta_0[i] = a_n[i];
        beta_c[i] = beta_0[i];
    }
    //  Calculates the subterm and term values
    Cox_Term_Risk_Calc(modelform, tform, term_n, totalnum, fir, dfc, term_tot, T0, Td0, Tdd0, Te, R, Rd, Rdd, Dose, nonDose, beta_0, df0, thres_step_max, step_max, TTerm, nonDose_LIN, nonDose_PLIN, nonDose_LOGLIN, RdR, RddR, nthreads, KeepConstant, verbose, model_bool, gmix_theta, gmix_term);
    //  Calculates the side sum terms used
    Cox_Side_LL_Calc(reqrdnum, ntime, tform, RiskFail,  RiskPairs, RiskPairs_Strata, totalnum, fir, R, Rd, Rdd, Rls1, Rls2, Rls3, Lls1, Lls2, Lls3, cens_weight, Strata_vals, beta_0, RdR, RddR, Ll, Lld, Lldd, nthreads, KeepConstant, ties_method, verbose, model_bool, iter_stop);
    //
    #ifdef _OPENMP
    #pragma omp parallel for schedule(dynamic) num_threads(nthreads)
    #endif
    for (int ijk = 0; ijk < reqrdnum*(reqrdnum + 1)/2; ijk++) {
        int ij = 0;
        int jk = ijk;
        while (jk > ij) {
            ij++;
            jk -= ij;
        }
        Lldd_vec[ij * reqrdnum + jk] = Lldd[ij*reqrdnum+jk];
        Lldd_vec[jk * reqrdnum + ij] = Lldd_vec[ij * reqrdnum + jk];
        if (ij == jk) {
            Lld_vecc[ij] = Lld[ij];
        }
    }
    Lldd_vec.attr("dim") = Dimension(reqrdnum, reqrdnum);
    Lldd_mat = as<Map<MatrixXd> >(Lldd_vec);
    Lld_vec = as<Map<VectorXd> >(Lld_vecc);
    step = 1;
    iter_continue = true;
    max_change = 100;
    while ((step < maxstep) && (iter_continue)) {
        Lldd_mat = as<Map<MatrixXd> >(Lldd_vec);
        Lld_vec = as<Map<VectorXd> >(Lld_vecc);
        step++;
        trouble = false;
        Log_Bound(deriv_max, Lldd_mat, Lld_vec, Lstar, qchi, Ll[0], para_number, nthreads, totalnum, reqrdnum, KeepConstant, term_tot, step, dbeta, beta_0, upper, trouble, verbose, mult);
        if (trouble) {
            Calc_Change_trouble(para_number, nthreads, totalnum, thres_step_max, lr, step_max, Ll, Lld, Lldd, dbeta, tform, thres_step_max, step_max, KeepConstant_trouble);
        }
        beta_p = beta_c;  //
        beta_a = beta_c;  //
        Cox_Refresh_R_SIDES(reqrdnum, ntime, Rls1, Rls2, Rls3, Lls1, Lls2, Lls3, Strata_vals, model_bool);
        //
        for (int ijk = 0; ijk < totalnum; ijk++) {
            if (KeepConstant[ijk] == 0) {
                //
                if ((tform[ijk] == "lin_quad_int") || (tform[ijk] == "lin_exp_int") || (tform[ijk] == "step_int") || (tform[ijk] == "lin_int")) {  //  the threshold values use different maximum deviation values
                    if (abs(dbeta[ijk]) > thres_step_max) {
                        dbeta[ijk] = thres_step_max * sign(dbeta[ijk]);
                    }
                } else {
                    if (abs(dbeta[ijk]) > step_max) {
                        dbeta[ijk] = step_max * sign(dbeta[ijk]);
                    }
                }
            } else {
                dbeta[ijk] = 0;
            }
        }
        //
        max_change = abs(dbeta[0]);
        for (int ij = 0; ij < totalnum; ij++) {
            if (ij == para_number) {
                //  prevent the parameter estimate from crossing the optimum
                //  issue is beta_0[para_number] <= beta_best[para_number]
                if (dbeta[ij] > (beta_peak[para_number] - beta_a[ij])/lr) {
                    dbeta[ij] = (beta_peak[para_number] - beta_a[ij])/lr/2;
                }
            }
            if (abs(dbeta[ij]) > max_change) {
                max_change = abs(dbeta[ij]);
            }
            beta_0[ij] = beta_a[ij] + lr*dbeta[ij];
            beta_c[ij] = beta_0[ij];
        }
        //  ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------//
        //  The same subterm, risk, sides, and log-likelihood calculations are performed every half-step and iteration
        //  ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------//
        Cox_Term_Risk_Calc(modelform, tform, term_n, totalnum, fir, dfc, term_tot, T0, Td0, Tdd0, Te, R, Rd, Rdd, Dose, nonDose, beta_0, df0, thres_step_max, step_max, TTerm, nonDose_LIN, nonDose_PLIN, nonDose_LOGLIN, RdR, RddR, nthreads, KeepConstant, verbose, model_bool, gmix_theta, gmix_term);
        bound_val = 0;
        Cox_Pois_Log_Loop(step_max, model_bool, beta_0, beta_a, beta_c, bound_val, dbeta, df0, dfc, dint, Dose, thres_step_max, dslp, fir, gmix_term, gmix_theta, half_check, halfmax, KeepConstant, limit_hit, lr, modelform, nonDose, nonDose_LIN, nonDose_LOGLIN, nonDose_PLIN, nthreads, R, Rd, Rdd, RddR, RdR, s_weights, T0, Td0, Tdd0, Te, term_n, term_tot, tform, totalnum, TTerm, verbose);
        if (limit_hit[0]) {
            limits[0] = 0;
            ll_final[0] = 0;
            limit_converged[0] = FALSE;
            break;
        }
        Cox_Side_LL_Calc(reqrdnum, ntime, tform, RiskFail,  RiskPairs, RiskPairs_Strata, totalnum, fir, R, Rd, Rdd, Rls1, Rls2, Rls3, Lls1, Lls2, Lls3, cens_weight, Strata_vals, beta_0, RdR, RddR, Ll, Lld, Lldd, nthreads, KeepConstant, ties_method, verbose, model_bool, iter_stop);
        //
        #ifdef _OPENMP
        #pragma omp parallel for schedule(dynamic) num_threads(nthreads)
        #endif
        for (int ijk = 0; ijk < reqrdnum*(reqrdnum + 1)/2; ijk++) {
            int ij = 0;
            int jk = ijk;
            while (jk > ij) {
                ij++;
                jk -= ij;
            }
            Lldd_vec[ij * reqrdnum + jk] = Lldd[ij*reqrdnum+jk];
            Lldd_vec[jk * reqrdnum + ij] = Lldd_vec[ij * reqrdnum + jk];
            if (ij == jk) {
                Lld_vecc[ij] = Lld[ij];
            }
        }
        Lldd_vec.attr("dim") = Dimension(reqrdnum, reqrdnum);
        limits[0] = beta_0[para_number];
        ll_final[0] = Ll[0];
        if ((max_change < epsilon) && (deriv_max < deriv_epsilon)) {
            iter_continue = false;
            limit_converged[0] = TRUE;
        }
    }
    //
    res_list = List::create(_["Parameter_Limits"] = wrap(limits), _["Negative_Risk_Limit_Hit"] = wrap(limit_hit), _["Likelihood_Boundary"] = wrap(ll_final), _["Likelihood_Goal"] = wrap(Lstar), _["Limit_Converged"] = wrap(limit_converged), _["Status"] = "PASSED");
    //  returns a list of results
    return res_list;
}

//' Primary Poisson likelihood bounds calculation function.
//'
//' \code{LogLik_Poisson_Omnibus_Log_Bound} Performs the calls to calculation functions and log-likeihood profile bounds
//'
//' @inheritParams CPP_template
//'
//' @return List of final results: Log-likelihood of optimum, first derivative of log-likelihood, second derivative matrix, parameter list, standard deviation estimate, AIC, model information
//' @noRd
//'
//
List LogLik_Poisson_Omnibus_Log_Bound(const Ref<const MatrixXd>& PyrC, const Ref<const MatrixXd>& dfs, IntegerVector term_n, StringVector tform, Ref<VectorXd> beta_0, Ref<MatrixXd> df0, IntegerVector dfc, int fir, string modelform, double lr, NumericVector maxiters, int guesses, int halfmax, double epsilon, double step_max, double thres_step_max, double deriv_epsilon, int verbose, IntegerVector KeepConstant, int term_tot, int nthreads, List model_bool, const double gmix_theta, const IntegerVector gmix_term, const Ref<const MatrixXd>& Lin_Sys, const Ref<const VectorXd>& Lin_Res, double qchi, int para_number, int maxstep, double mult) {
    //
    List temp_list = List::create(_["Status"] = "TEMP");  //  used as a dummy return value for code checking
    //
    //  Time durations are measured from this point on in microseconds
    //
    //  df0: covariate data
    //  totalnum: number of terms used
    //
    //  ------------------------------------------------------------------------- //  initialize
    const int mat_row = df0.rows();
    int totalnum = term_n.size();
    int reqrdnum = totalnum - sum(KeepConstant);
    //  ------------------------------------------------------------------------- //  initialize
    if (model_bool["single"]) {
        if (verbose >= 1) {
            Rcout << "non-derivative model calculation is not compatable with log-based bound calculation" << endl;
        }
        temp_list = List::create(_["Status"] = "FAILED_WITH_BAD_MODEL_SINGLE", _["LogLik"] = R_NaN);
        return temp_list;
    }
    if (model_bool["gradient"]) {
        if (verbose >= 1) {
            Rcout << "gradient descent model calculation is not compatable with log-based bound calculation" << endl;
        }
        temp_list = List::create(_["Status"] = "FAILED_WITH_BAD_MODEL_GRADIENT", _["LogLik"] = R_NaN);
        return temp_list;
    }
    //
    //  cout.precision: controls the number of significant digits printed
    //
    Rcout.precision(7);  //  forces higher precision numbers printed to terminal
    //  ------------------------------------------------------------------------- //  initialize
    // Map<VectorXd> beta_0(as<Map<VectorXd> >(a_n));
    MatrixXd T0;
    MatrixXd Td0;
    MatrixXd Tdd0;
    MatrixXd Te;
    MatrixXd R;
    MatrixXd Rd;
    MatrixXd Rdd;
    MatrixXd Dose;
    MatrixXd nonDose;
    MatrixXd nonDose_LIN;
    MatrixXd nonDose_PLIN;
    MatrixXd nonDose_LOGLIN;
    MatrixXd TTerm;
    double dint = 0.0;  //  the amount of change used to calculate derivatives in threshold paramters
    double dslp = 0.0;
    MatrixXd RdR;
    MatrixXd RddR;
    //  ------------------------------------------------------------------------- //  initialize
    //  ---------------------------------------------
    //  To Start, needs to seperate the derivative terms
    //  ---------------------------------------------
    //
    T0 = MatrixXd::Zero(mat_row, totalnum);  //  preallocates matrix for Term column
    Te = MatrixXd::Zero(mat_row, 1);  //  preallocates matrix for column terms used for temporary storage
    R = MatrixXd::Zero(mat_row, 1);  //  preallocates matrix for Risks
    //
    Dose = MatrixXd::Constant(mat_row, term_tot, 0.0);  //  matrix of the total dose term values
    nonDose = MatrixXd::Constant(mat_row, term_tot, 1.0);  //  matrix of the total non-dose term values
    nonDose_LIN = MatrixXd::Constant(mat_row, term_tot, 0.0);  //  matrix of Linear subterm values
    nonDose_PLIN = MatrixXd::Constant(mat_row, term_tot, 1.0);  //  matrix of Loglinear subterm values
    nonDose_LOGLIN = MatrixXd::Constant(mat_row, term_tot, 1.0);  //  matrix of Product linear subterm values
    TTerm = MatrixXd::Zero(mat_row, term_tot);  //  matrix of term values
    if (model_bool["single"]) {
    } else {
        Td0 = MatrixXd::Zero(mat_row, reqrdnum);  //  preallocates matrix for Term derivative columns
        Tdd0 = MatrixXd::Zero(mat_row, reqrdnum*(reqrdnum + 1)/2);  //  preallocates matrix for Term second derivative columns
        //
        Rd = MatrixXd::Zero(mat_row, reqrdnum);  //  preallocates matrix for Risk derivatives
        Rdd = MatrixXd::Zero(mat_row, reqrdnum*(reqrdnum + 1)/2);  //  preallocates matrix for Risk second derivatives
        //
        dint = thres_step_max;  //  the amount of change used to calculate derivatives in threshold paramters
        dslp = step_max;
        RdR = MatrixXd::Zero(mat_row, reqrdnum);  //  preallocates matrix for Risk to derivative ratios
        RddR = MatrixXd::Zero(mat_row, reqrdnum*(reqrdnum + 1)/2);  //  preallocates matrix for Risk to second derivative ratios
    }
    VectorXd s_weights;
    if (model_bool["strata"]) {
        s_weights = VectorXd::Zero(mat_row);
        Gen_Strat_Weight(modelform, dfs, PyrC, s_weights, nthreads, tform, term_n, term_tot);
    }
    //  ------------------------------------------------------------------------- //  initialize
    vector<double> Ll(reqrdnum, 0.0);  //  log-likelihood values
    vector<double> Lld(reqrdnum, 0.0);  //  log-likelihood derivative values
    vector<double> Lldd(pow(reqrdnum, 2), 0.0);  //  the second derivative matrix has room for every combination, but only the lower triangle is calculated initially
    MatrixXd dev_temp = MatrixXd::Zero(PyrC.rows(), 2);
    double dev = 0;
    //  the log-likelihood is calculated in parallel over the risk groups
    vector <double> Ll_comp(2, Ll[0]);  //  vector to compare values
    double step_max0 = step_max;
    double thres_step_max0 = thres_step_max;
    //
    vector<double> dbeta(totalnum, 0.0);
    //
    //  --------------------------
    //  always starts from initial guess
    //  --------------------------
    vector<double> beta_c(totalnum, 0.0);
    vector<double> beta_a(totalnum, 0.0);
    vector<double> beta_best(totalnum, 0.0);
    vector<double> beta_peak(totalnum, 0.0);
    vector<double> beta_p(totalnum, 0.0);
    VectorXd::Map(&beta_p[0], beta_0.size()) = beta_0;  //  stores previous parameters
    VectorXd::Map(&beta_c[0], beta_0.size()) = beta_0;  //  stores current parameters
    VectorXd::Map(&beta_a[0], beta_0.size()) = beta_0;  //  stores a refrence value for parameters
    VectorXd::Map(&beta_best[0], beta_0.size()) = beta_0;  //  stores the best parameters
    VectorXd::Map(&beta_peak[0], beta_0.size()) = beta_0;  //  stores the best parameters
    int iter_stop  = 0;  //  tracks if the iterations should be stopped for convergence
    vector<double> beta_abs_best(totalnum, 0.0);
    //
    fill(Ll.begin(), Ll.end(), 0.0);
    fill(Lld.begin(), Lld.end(), 0.0);
    fill(Lldd.begin(), Lldd.end(), 0.0);
    beta_p = beta_peak;  //
    beta_a = beta_peak;  //
    beta_c = beta_peak;  //
    step_max = step_max0;
    thres_step_max = thres_step_max0;
    //
    for (int i = 0; i < beta_0.size(); i++) {
        beta_0[i] = beta_peak[i];
    }
    Pois_Term_Risk_Calc(modelform, tform, term_n, totalnum, fir, dfc, term_tot, T0, Td0, Tdd0, Te, R, Rd, Rdd, Dose, nonDose, beta_0, df0, dint, dslp, TTerm, nonDose_LIN, nonDose_PLIN, nonDose_LOGLIN, RdR, RddR, s_weights, nthreads, KeepConstant, verbose, model_bool, gmix_theta, gmix_term);
    Pois_Dev_LL_Calc(reqrdnum, totalnum, fir, R, Rd, Rdd, beta_0, RdR, RddR, Ll, Lld, Lldd, PyrC, dev_temp, nthreads, KeepConstant, verbose, model_bool, iter_stop, dev);
    NumericVector Lldd_vec(reqrdnum * reqrdnum);
    NumericVector Lld_vecc(reqrdnum);
    #ifdef _OPENMP
    #pragma omp parallel for schedule(dynamic) num_threads(nthreads)
    #endif
    for (int ijk = 0; ijk < reqrdnum*(reqrdnum + 1)/2; ijk++) {
        int ij = 0;
        int jk = ijk;
        while (jk > ij) {
            ij++;
            jk -= ij;
        }
        Lldd_vec[ij * reqrdnum + jk] = Lldd[ij*reqrdnum+jk];
        Lldd_vec[jk * reqrdnum + ij] = Lldd_vec[ij * reqrdnum + jk];
        if (ij == jk) {
            Lld_vecc[ij] = Lld[ij];
        }
    }
    Lldd_vec.attr("dim") = Dimension(reqrdnum, reqrdnum);
    Map<MatrixXd> Lldd_mat(as<Map<MatrixXd> >(Lldd_vec));
    Map<VectorXd> Lld_vec(as<Map<VectorXd> >(Lld_vecc));
    double Lstar = Ll[0]-qchi;
    bool upper = true;
    int half_check = 0;
    bool trouble = false;
    IntegerVector KeepConstant_trouble(totalnum);
    for (int i = 0;  i < totalnum; i++) {
        KeepConstant_trouble[i] = KeepConstant[i];
    }
    KeepConstant_trouble[para_number] = 1;
    //
    vector<double> limits(2, 0.0);
    vector<bool>   limit_hit(2, FALSE);
    vector<bool>   limit_converged(2, FALSE);
    vector<double> ll_final(2, 0.0);
    List res_list;
    //
    if (verbose >= 4) {
        Rcout << "C++ Note: STARTING Upper Bound" << endl;
    }
    upper = true;
    int step = -1;
    bool iter_continue = true;
    double max_change = 100;
    double deriv_max = 100;
    ///
    //  variables added for log loop code
    ///
    int bound_val = 1;
    //
    while ((step < maxstep) && (iter_continue)) {
        step++;
        //
        trouble = false;
        half_check = 0;
        Log_Bound(deriv_max, Lldd_mat, Lld_vec, Lstar, qchi, Ll[0], para_number, nthreads, totalnum, reqrdnum, KeepConstant, term_tot, step, dbeta, beta_0, upper, trouble, verbose, mult);
        if (trouble) {
            Calc_Change_trouble(para_number, nthreads, totalnum, thres_step_max, lr, step_max, Ll, Lld, Lldd, dbeta, tform, thres_step_max, step_max, KeepConstant_trouble);
        }
        //
        beta_p = beta_c;  //
        beta_a = beta_c;  //

        //
        for (int ijk = 0; ijk < totalnum; ijk++) {
            if (KeepConstant[ijk] == 0) {
                //
                if ((tform[ijk] == "lin_quad_int") || (tform[ijk] == "lin_exp_int") || (tform[ijk] == "step_int") || (tform[ijk] == "lin_int")) {  //  the threshold values use different maximum deviation values
                    if (abs(dbeta[ijk]) > thres_step_max) {
                        dbeta[ijk] = thres_step_max * sign(dbeta[ijk]);
                    }
                } else {
                    if (abs(dbeta[ijk]) > step_max) {
                        dbeta[ijk] = step_max * sign(dbeta[ijk]);
                    }
                }
            } else {
                dbeta[ijk] = 0;
            }
        }
        //
        max_change = abs(dbeta[0]);
        for (int ij = 0; ij < totalnum; ij++) {
            if (ij == para_number) {
                //  prevent the parameter estimate from crossing the optimum
                //  issue is beta_0[para_number] <= beta_best[para_number]
                if (dbeta[ij] < (beta_peak[para_number] - beta_a[ij])/lr) {
                    dbeta[ij] = (beta_peak[para_number] - beta_a[ij])/lr/2;
                }
            }
            if (abs(dbeta[ij]) > max_change) {
                max_change = abs(dbeta[ij]);
            }
            beta_0[ij] = beta_a[ij] + lr*dbeta[ij];
            beta_c[ij] = beta_0[ij];
        }
        //  ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------//
        //  The same subterm, risk, sides, and log-likelihood calculations are performed every half-step and iteration
        //  ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------//

        Pois_Term_Risk_Calc(modelform, tform, term_n, totalnum, fir, dfc, term_tot, T0, Td0, Tdd0, Te, R, Rd, Rdd, Dose, nonDose, beta_0, df0, dint, dslp, TTerm, nonDose_LIN, nonDose_PLIN, nonDose_LOGLIN, RdR, RddR, s_weights, nthreads, KeepConstant, verbose, model_bool, gmix_theta, gmix_term);
        bound_val = 1;
        Cox_Pois_Log_Loop(step_max, model_bool, beta_0, beta_a, beta_c, bound_val, dbeta, df0, dfc, dint, Dose, thres_step_max, dslp, fir, gmix_term, gmix_theta, half_check, halfmax, KeepConstant, limit_hit, lr, modelform, nonDose, nonDose_LIN, nonDose_LOGLIN, nonDose_PLIN, nthreads, R, Rd, Rdd, RddR, RdR, s_weights, T0, Td0, Tdd0, Te, term_n, term_tot, tform, totalnum, TTerm, verbose);
        if (limit_hit[1]) {
            limits[1] = 0;
            ll_final[1] = 0;
            limit_converged[1] = FALSE;
            break;
        }
        Pois_Dev_LL_Calc(reqrdnum, totalnum, fir, R, Rd, Rdd, beta_0, RdR, RddR, Ll, Lld, Lldd, PyrC, dev_temp, nthreads, KeepConstant, verbose, model_bool, iter_stop, dev);
        //
        #ifdef _OPENMP
        #pragma omp parallel for schedule(dynamic) num_threads(nthreads)
        #endif
        for (int ijk = 0; ijk < reqrdnum*(reqrdnum + 1)/2; ijk++) {
            int ij = 0;
            int jk = ijk;
            while (jk > ij) {
                ij++;
                jk -= ij;
            }
            Lldd_vec[ij * reqrdnum + jk] = Lldd[ij*reqrdnum+jk];
            Lldd_vec[jk * reqrdnum + ij] = Lldd_vec[ij * reqrdnum + jk];
            if (ij == jk) {
                Lld_vecc[ij] = Lld[ij];
            }
        }
        Lldd_vec.attr("dim") = Dimension(reqrdnum, reqrdnum);
        Lldd_mat = as<Map<MatrixXd> >(Lldd_vec);
        Lld_vec = as<Map<VectorXd> >(Lld_vecc);
        limits[1] = beta_0[para_number];
        ll_final[1] = Ll[0];
        if ((max_change < epsilon) && (deriv_max < deriv_epsilon)) {
            iter_continue = false;
            limit_converged[1] = TRUE;
        }
    }
    if (verbose >= 4) {
        Rcout << "C++ Note: STARTING Lower Bound" << endl;
    }
    beta_p = beta_peak;  //
    beta_a = beta_peak;  //
    beta_c = beta_peak;  //
    for (int ij = 0; ij < totalnum; ij++) {
        beta_0[ij] = beta_a[ij];
    }
    //
    Pois_Term_Risk_Calc(modelform, tform, term_n, totalnum, fir, dfc, term_tot, T0, Td0, Tdd0, Te, R, Rd, Rdd, Dose, nonDose, beta_0, df0, dint, dslp, TTerm, nonDose_LIN, nonDose_PLIN, nonDose_LOGLIN, RdR, RddR, s_weights, nthreads, KeepConstant, verbose, model_bool, gmix_theta, gmix_term);
    Pois_Dev_LL_Calc(reqrdnum, totalnum, fir, R, Rd, Rdd, beta_0, RdR, RddR, Ll, Lld, Lldd, PyrC, dev_temp, nthreads, KeepConstant, verbose, model_bool, iter_stop, dev);
    #ifdef _OPENMP
    #pragma omp parallel for schedule(dynamic) num_threads(nthreads)
    #endif
    for (int ijk = 0; ijk < reqrdnum*(reqrdnum + 1)/2; ijk++) {
        int ij = 0;
        int jk = ijk;
        while (jk > ij) {
            ij++;
            jk -= ij;
        }
        Lldd_vec[ij * reqrdnum + jk] = Lldd[ij*reqrdnum+jk];
        Lldd_vec[jk * reqrdnum + ij] = Lldd_vec[ij * reqrdnum + jk];
        if (ij == jk) {
            Lld_vecc[ij] = Lld[ij];
        }
    }
    Lldd_vec.attr("dim") = Dimension(reqrdnum, reqrdnum);
    Lldd_mat = as<Map<MatrixXd> >(Lldd_vec);
    Lld_vec = as<Map<VectorXd> >(Lld_vecc);
    //
    step_max = step_max0;
    thres_step_max = thres_step_max0;
    //
    upper = false;
    step = -1;
    iter_continue = true;
    while ((step < maxstep) && (iter_continue)) {
        step++;
        trouble = false;
        half_check = 0;
        Log_Bound(deriv_max, Lldd_mat, Lld_vec, Lstar, qchi, Ll[0], para_number, nthreads, totalnum, reqrdnum, KeepConstant, term_tot, step, dbeta, beta_0, upper, trouble, verbose, mult);
        if (trouble) {
            Calc_Change_trouble(para_number, nthreads, totalnum, thres_step_max, lr, step_max, Ll, Lld, Lldd, dbeta, tform, thres_step_max, step_max, KeepConstant_trouble);
        }
        beta_p = beta_c;  //
        beta_a = beta_c;  //
        for (int ijk = 0; ijk < totalnum; ijk++) {
            if (KeepConstant[ijk] == 0) {
                //
                if ((tform[ijk] == "lin_quad_int") || (tform[ijk] == "lin_exp_int") || (tform[ijk] == "step_int") || (tform[ijk] == "lin_int")) {  //  the threshold values use different maximum deviation values
                    if (abs(dbeta[ijk]) > thres_step_max) {
                        dbeta[ijk] = thres_step_max * sign(dbeta[ijk]);
                    }
                } else {
                    if (abs(dbeta[ijk]) > step_max) {
                        dbeta[ijk] = step_max * sign(dbeta[ijk]);
                    }
                }
            } else {
                dbeta[ijk] = 0;
            }
        }
        //
        max_change = abs(dbeta[0]);
        for (int ij = 0; ij < totalnum; ij++) {
            if (ij == para_number) {
                //  prevent the parameter estimate from crossing the optimum
                //  issue is beta_0[para_number] <= beta_best[para_number]
                if (dbeta[ij] > (beta_peak[para_number] - beta_a[ij])/lr) {
                    dbeta[ij] = (beta_peak[para_number] - beta_a[ij])/lr/2;
                }
            }
            if (abs(dbeta[ij]) > max_change) {
                max_change = abs(dbeta[ij]);
            }
            beta_0[ij] = beta_a[ij] + lr*dbeta[ij];
            beta_c[ij] = beta_0[ij];
        }
        //  ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------//
        //  The same subterm, risk, sides, and log-likelihood calculations are performed every half-step and iteration
        //  ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------//
        Pois_Term_Risk_Calc(modelform, tform, term_n, totalnum, fir, dfc, term_tot, T0, Td0, Tdd0, Te, R, Rd, Rdd, Dose, nonDose, beta_0, df0, dint, dslp, TTerm, nonDose_LIN, nonDose_PLIN, nonDose_LOGLIN, RdR, RddR, s_weights, nthreads, KeepConstant, verbose, model_bool, gmix_theta, gmix_term);
        bound_val = 0;
        Cox_Pois_Log_Loop(step_max, model_bool, beta_0, beta_a, beta_c, bound_val, dbeta, df0, dfc, dint, Dose, thres_step_max, dslp, fir, gmix_term, gmix_theta, half_check, halfmax, KeepConstant, limit_hit, lr, modelform, nonDose, nonDose_LIN, nonDose_LOGLIN, nonDose_PLIN, nthreads, R, Rd, Rdd, RddR, RdR, s_weights, T0, Td0, Tdd0, Te, term_n, term_tot, tform, totalnum, TTerm, verbose);
        if (limit_hit[0]) {
            limits[0] = 0;
            ll_final[0] = 0;
            limit_converged[0] = FALSE;
            break;
        }
        Pois_Dev_LL_Calc(reqrdnum, totalnum, fir, R, Rd, Rdd, beta_0, RdR, RddR, Ll, Lld, Lldd, PyrC, dev_temp, nthreads, KeepConstant, verbose, model_bool, iter_stop, dev);
        //
        #ifdef _OPENMP
        #pragma omp parallel for schedule(dynamic) num_threads(nthreads)
        #endif
        for (int ijk = 0; ijk < reqrdnum*(reqrdnum + 1)/2; ijk++) {
            int ij = 0;
            int jk = ijk;
            while (jk > ij) {
                ij++;
                jk -= ij;
            }
            Lldd_vec[ij * reqrdnum + jk] = Lldd[ij*reqrdnum+jk];
            Lldd_vec[jk * reqrdnum + ij] = Lldd_vec[ij * reqrdnum + jk];
            if (ij == jk) {
                Lld_vecc[ij] = Lld[ij];
            }
        }
        Lldd_vec.attr("dim") = Dimension(reqrdnum, reqrdnum);
        Lldd_mat = as<Map<MatrixXd> >(Lldd_vec);
        Lld_vec = as<Map<VectorXd> >(Lld_vecc);
        limits[0] = beta_0[para_number];
        ll_final[0] = Ll[0];
        if ((max_change < epsilon) && (deriv_max < deriv_epsilon)) {
            iter_continue = false;
            limit_converged[0] = TRUE;
        }
    }
    res_list = List::create(_["Parameter_Limits"] = wrap(limits), _["Negative_Risk_Limit_Hit"] = wrap(limit_hit), _["Likelihood_Boundary"] = wrap(ll_final), _["Likelihood_Goal"] = wrap(Lstar), _["Limit_Converged"] = wrap(limit_converged), _["Status"] = "PASSED");
    //  returns a list of results
    return res_list;
}

//' Primary Poisson likelihood bounds calculation function.
//'
//' \code{LogLik_Poisson_Omnibus_Log_Bound_Search} Performs the calls to calculation functions and log-likeihood profile bounds
//'
//' @inheritParams CPP_template
//'
//' @return List of final results: Log-likelihood of optimum, first derivative of log-likelihood, second derivative matrix, parameter list, standard deviation estimate, AIC, model information
//' @noRd
//'
//
List LogLik_Poisson_Omnibus_Log_Bound_Search(const Ref<const MatrixXd>& PyrC, const Ref<const MatrixXd>& dfs, IntegerVector term_n, StringVector tform, Ref<VectorXd> beta_0, Ref<MatrixXd> df0, IntegerVector dfc, int fir, string modelform, double lr, NumericVector maxiters, int guesses, int halfmax, double epsilon, double step_max, double thres_step_max, double deriv_epsilon, int verbose, IntegerVector KeepConstant, int term_tot, int nthreads, List model_bool, const double gmix_theta, const IntegerVector gmix_term, const Ref<const MatrixXd>& Lin_Sys, const Ref<const VectorXd>& Lin_Res, double qchi, int para_number, int maxstep, double mult) {
    //
    List temp_list = List::create(_["Status"] = "TEMP");  //  used as a dummy return value for code checking
    //
    //  Time durations are measured from this point on in microseconds
    //
    //  df0: covariate data
    //  ntime: number of event times for Cox PH
    //  totalnum: number of terms used
    //
    //  ------------------------------------------------------------------------- //  initialize
    const int mat_row = df0.rows();
    int totalnum = term_n.size();
    int reqrdnum = totalnum - sum(KeepConstant);
    //  ------------------------------------------------------------------------- //  initialize
    if (model_bool["single"]) {
        if (verbose >= 1) {
            Rcout << "non-derivative model calculation is not compatable with log-based bound calculation" << endl;
        }
        temp_list = List::create(_["Status"] = "FAILED_WITH_BAD_MODEL_SINGLE", _["LogLik"] = R_NaN);
        return temp_list;
    }
    if (model_bool["gradient"]) {
        if (verbose >= 1) {
            Rcout << "gradient descent model calculation is not compatable with log-based bound calculation" << endl;
        }
        temp_list = List::create(_["Status"] = "FAILED_WITH_BAD_MODEL_GRADIENT", _["LogLik"] = R_NaN);
        return temp_list;
    }
    //  cout.precision: controls the number of significant digits printed
    //
    double dbeta_max = 0.0;
    Rcout.precision(7);  //  forces higher precision numbers printed to terminal
    //  ------------------------------------------------------------------------- //  initialize
    // Map<VectorXd> beta_0(as<Map<VectorXd> >(a_n));
    VectorXd beta_max = beta_0;
    MatrixXd T0;
    MatrixXd Td0;
    MatrixXd Tdd0;
    //
    MatrixXd Te;
    MatrixXd R;
    MatrixXd Rd;
    MatrixXd Rdd;
    //
    MatrixXd Dose;
    MatrixXd nonDose;
    MatrixXd nonDose_LIN;
    MatrixXd nonDose_PLIN;
    MatrixXd nonDose_LOGLIN;
    MatrixXd TTerm;
    double dint = 0.0;  //  the amount of change used to calculate derivatives in threshold paramters
    double dslp = 0.0;
    MatrixXd RdR;
    MatrixXd RddR;
    //  ------------------------------------------------------------------------- //  initialize
    //  ---------------------------------------------
    //  To Start, needs to seperate the derivative terms
    //  ---------------------------------------------
    //
    T0 = MatrixXd::Zero(mat_row, totalnum);  //  preallocates matrix for Term column
    Te = MatrixXd::Zero(mat_row, 1);  //  preallocates matrix for column terms used for temporary storage
    R = MatrixXd::Zero(mat_row, 1);  //  preallocates matrix for Risks
    //
    Dose = MatrixXd::Constant(mat_row, term_tot, 0.0);  //  matrix of the total dose term values
    nonDose = MatrixXd::Constant(mat_row, term_tot, 1.0);  //  matrix of the total non-dose term values
    nonDose_LIN = MatrixXd::Constant(mat_row, term_tot, 0.0);  //  matrix of Linear subterm values
    nonDose_PLIN = MatrixXd::Constant(mat_row, term_tot, 1.0);  //  matrix of Loglinear subterm values
    nonDose_LOGLIN = MatrixXd::Constant(mat_row, term_tot, 1.0);  //  matrix of Product linear subterm values
    TTerm = MatrixXd::Zero(mat_row, term_tot);  //  matrix of term values
    if (model_bool["single"]) {
    } else {
        Td0 = MatrixXd::Zero(mat_row, reqrdnum);  //  preallocates matrix for Term derivative columns
        Tdd0 = MatrixXd::Zero(mat_row, reqrdnum*(reqrdnum + 1)/2);  //  preallocates matrix for Term second derivative columns
        Rd = MatrixXd::Zero(mat_row, reqrdnum);  //  preallocates matrix for Risk derivatives
        Rdd = MatrixXd::Zero(mat_row, reqrdnum*(reqrdnum + 1)/2);  //  preallocates matrix for Risk second derivatives
        dint = thres_step_max;  //  the amount of change used to calculate derivatives in threshold paramters
        dslp = step_max;
        RdR = MatrixXd::Zero(mat_row, reqrdnum);  //  preallocates matrix for Risk to derivative ratios
        RddR = MatrixXd::Zero(mat_row, reqrdnum*(reqrdnum + 1)/2);  //  preallocates matrix for Risk to second derivative ratios
    }
    VectorXd s_weights;
    if (model_bool["strata"]) {
        s_weights = VectorXd::Zero(mat_row);
        Gen_Strat_Weight(modelform, dfs, PyrC, s_weights, nthreads, tform, term_n, term_tot);
    }
    //  ------------------------------------------------------------------------- //  initialize
    MatrixXd Rls1;
    MatrixXd Lls1;
    MatrixXd Rls2;
    MatrixXd Rls3;
    MatrixXd Lls2;
    MatrixXd Lls3;
    vector<double> Ll(reqrdnum, 0.0);  //  log-likelihood values
    vector<double> Lld(reqrdnum, 0.0);  //  log-likelihood derivative values
    vector<double> Lldd(pow(reqrdnum, 2), 0.0);  //  the second derivative matrix has room for every combination, but only the lower triangle is calculated initially
    MatrixXd dev_temp = MatrixXd::Zero(PyrC.rows(), 2);
    double dev = 0;
    //  ------------------------------------------------------------------------- //  initialize
    //  the log-likelihood is calculated in parallel over the risk groups
    vector <double> Ll_comp(2, Ll[0]);  //  vector to compare values
    double step_max0 = step_max;
    double thres_step_max0 = thres_step_max;
    //
    vector<double> dbeta(totalnum, 0.0);
    vector<double> dbeta_start(totalnum, 0.0);
    //
    //  --------------------------
    //  always starts from initial guess
    //  --------------------------
    vector<double> beta_c(totalnum, 0.0);
    vector<double> beta_a(totalnum, 0.0);
    vector<double> beta_best(totalnum, 0.0);
    vector<double> beta_p(totalnum, 0.0);
    vector<double> beta_peak(totalnum, 0.0);
    VectorXd::Map(&beta_p[0], beta_0.size()) = beta_0;  //  stores previous parameters
    VectorXd::Map(&beta_c[0], beta_0.size()) = beta_0;  //  stores current parameters
    VectorXd::Map(&beta_a[0], beta_0.size()) = beta_0;  //  stores a refrence value for parameters
    VectorXd::Map(&beta_best[0], beta_0.size()) = beta_0;  //  stores the best parameters
    VectorXd::Map(&beta_peak[0], beta_0.size()) = beta_0;  //  stores the peak parameters
    int iter_stop  = 0;  //  tracks if the iterations should be stopped for convergence
    vector<double> beta_abs_best(totalnum, 0.0);
    fill(Ll.begin(), Ll.end(), 0.0);
    fill(Lld.begin(), Lld.end(), 0.0);
    fill(Lldd.begin(), Lldd.end(), 0.0);
    beta_p = beta_peak;  //
    beta_a = beta_peak;  //
    beta_c = beta_peak;  //
    step_max = step_max0;
    thres_step_max = thres_step_max0;
    for (int i = 0; i < beta_0.size(); i++) {
        beta_0[i] = beta_peak[i];
    }
    Pois_Term_Risk_Calc(modelform, tform, term_n, totalnum, fir, dfc, term_tot, T0, Td0, Tdd0, Te, R, Rd, Rdd, Dose, nonDose, beta_0, df0, dint, dslp, TTerm, nonDose_LIN, nonDose_PLIN, nonDose_LOGLIN, RdR, RddR, s_weights, nthreads, KeepConstant, verbose, model_bool, gmix_theta, gmix_term);
    Pois_Dev_LL_Calc(reqrdnum, totalnum, fir, R, Rd, Rdd, beta_0, RdR, RddR, Ll, Lld, Lldd, PyrC, dev_temp, nthreads, KeepConstant, verbose, model_bool, iter_stop, dev);
    //
    NumericVector Lldd_vec(reqrdnum * reqrdnum);
    NumericVector Lld_vecc(reqrdnum);
    #ifdef _OPENMP
    #pragma omp parallel for schedule(dynamic) num_threads(nthreads)
    #endif
    for (int ijk = 0; ijk < reqrdnum*(reqrdnum + 1)/2; ijk++) {
        int ij = 0;
        int jk = ijk;
        while (jk > ij) {
            ij++;
            jk -= ij;
        }
        Lldd_vec[ij * reqrdnum + jk] = Lldd[ij*reqrdnum+jk];
        Lldd_vec[jk * reqrdnum + ij] = Lldd_vec[ij * reqrdnum + jk];
        if (ij == jk) {
            Lld_vecc[ij] = Lld[ij];
        }
    }
    Lldd_vec.attr("dim") = Dimension(reqrdnum, reqrdnum);
    Map<MatrixXd> Lldd_mat(as<Map<MatrixXd> >(Lldd_vec));
    Map<VectorXd> Lld_vec(as<Map<VectorXd> >(Lld_vecc));
    double Lstar = Ll[0]-qchi;
    double Lmax = Ll[0];
    bool upper = true;
    int half_check = 0;
    bool trouble = false;
    double deriv_max = 100;
    //
    IntegerVector KeepConstant_trouble(totalnum);
    for (int i = 0;  i < totalnum; i++) {
        KeepConstant_trouble[i] = KeepConstant[i];
    }
    KeepConstant_trouble[para_number] = 1;
    //
    vector<double> limits(2, 0.0);
    vector<bool>   limit_hit(2, FALSE);
    vector<bool>   limit_converged(2, FALSE);
    vector<double> ll_final(2, 0.0);
    List res_list;
    upper = true;
    //  Now define the list of points to check
    trouble = false;
    Log_Bound(deriv_max, Lldd_mat, Lld_vec, Lstar, qchi, Lmax, para_number, nthreads, totalnum, reqrdnum, KeepConstant, term_tot, 0, dbeta_start, beta_0, upper, trouble, verbose, mult);
    NumericMatrix a_ns(guesses, totalnum);
    NumericVector a_n = a_ns.row(0);
    //  now the dbeta holds the range to check
    //  note that currently the guesses are not bounded
    for (int i = 0; i < guesses; i++) {
        for (int j = 0; j < totalnum; j++) {
            //  use dbeta to assign a_n guess i, parameter j
            //  assign such that i=guesses - 1 gives mult*dbeta
            a_ns(i, j) = beta_0[j] + dbeta_start[j] * (i + 1)/(guesses);
        }
    }
    if (verbose >= 4) {
        for (int i = 0; i < guesses; i++) {
            Rcout << "C++ Note: Initial guess " << i << ": ";
            for (int j = 0; j < totalnum; j++) {
                Rcout << a_ns(i, j) << " ";
            }
            Rcout << " " << endl;
        }
    }
    //  now we have the points to test
    double halves = 0;  //  number of half-steps taken
    int ind0 = fir;  //  used for validations
    int iteration = 0;  //  iteration number
    int maxiter = 0;
    iter_stop  = 0;  //  tracks if the iterations should be stopped for convergence
    NumericMatrix beta_fin(a_ns.rows(), a_ns.cols());
    NumericVector LL_fin(a_ns.rows());
    double Ll_iter_best = 10;
    //  Variables that are used for the risk check function shared across cox, poisson, and log bound functions
    VectorXd cens_weight(1);
    int ntime = 1.0;
    IntegerMatrix RiskFail(1);
    vector<vector<int> > RiskPairs;
    vector<vector<vector<int> > > RiskPairs_Strata;
    NumericVector Strata_vals(1);
    string ties_method = "temp";
    //
    for (int guess = 0; guess <guesses; guess++) {
        fill(Ll.begin(), Ll.end(), 0.0);
        fill(Lld.begin(), Lld.end(), 0.0);
        fill(Lldd.begin(), Lldd.end(), 0.0);
        beta_p = beta_peak;  //
        beta_a = beta_peak;  //
        beta_c = beta_peak;  //
        step_max = step_max0;
        thres_step_max = thres_step_max0;
        iter_stop = 0;
        halves = 0;
        iteration = 0;
        halves = 0;  //  number of half-steps taken
        ind0 = fir;  //  used for validations
        iteration = 0;  //  iteration number
        Ll_iter_best = 10;
        iter_stop  = 0;  //  tracks if the iterations should be stopped for convergence
        maxiter = maxiters[0];
        a_n = a_ns.row(guess);
        for (int i = 0; i < beta_0.size(); i++) {
            beta_0[i] = a_n[i];
        }
        Pois_Term_Risk_Calc(modelform, tform, term_n, totalnum, fir, dfc, term_tot, T0, Td0, Tdd0, Te, R, Rd, Rdd, Dose, nonDose, beta_0, df0, dint, dslp, TTerm, nonDose_LIN, nonDose_PLIN, nonDose_LOGLIN, RdR, RddR, s_weights, nthreads, KeepConstant, verbose, model_bool, gmix_theta, gmix_term);
        if ((R.minCoeff() <= 0) || (R.hasNaN())) {
            iter_stop = 1;
            Ll[0] = 404;  //  note to consider this point too far
        } else {
            Pois_Dev_LL_Calc(reqrdnum, totalnum, fir, R, Rd, Rdd, beta_0, RdR, RddR, Ll, Lld, Lldd, PyrC, dev_temp, nthreads, KeepConstant, verbose, model_bool, iter_stop, dev);
            for (int i = 0; i < beta_0.size(); i++) {
                beta_c[i] = beta_0[i];
            }
        }
        while ((iteration < maxiter) && (iter_stop == 0)) {
            iteration++;
            beta_p = beta_c;  //
            beta_a = beta_c;  //
            beta_best = beta_c;  //
            //  calculates the initial change in parameter
            Calc_Change_trouble(para_number, nthreads, totalnum, thres_step_max, lr, step_max, Ll, Lld, Lldd, dbeta, tform, thres_step_max, step_max, KeepConstant_trouble);
            Intercept_Bound(nthreads, totalnum, beta_0, dbeta, dfc, df0, KeepConstant, tform);
            //
            if ((Ll_iter_best > 0) || (Ll_iter_best < Ll[ind0])) {
                Ll_iter_best = Ll[ind0];
                beta_abs_best = beta_c;
            }
            halves = 0;
            while ((Ll[ind0] <= Ll_iter_best) && (halves < halfmax)) {  //  repeats until half-steps maxed or an improvement
                for (int ij = 0; ij < totalnum; ij++) {
                    beta_0[ij] = beta_a[ij] + dbeta[ij];
                    beta_c[ij] = beta_0[ij];
                }
                //  ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------//
                //  The same subterm, risk, sides, and log-likelihood calculations are performed every half-step and iteration
                //  ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------//
                Pois_Term_Risk_Calc(modelform, tform, term_n, totalnum, fir, dfc, term_tot, T0, Td0, Tdd0, Te, R, Rd, Rdd, Dose, nonDose, beta_0, df0, dint, dslp, TTerm, nonDose_LIN, nonDose_PLIN, nonDose_LOGLIN, RdR, RddR, s_weights, nthreads, KeepConstant, verbose, model_bool, gmix_theta, gmix_term);
                Cox_Pois_Check_Continue(model_bool, beta_0, beta_best, beta_c, cens_weight, dbeta, dev, dev_temp, fir, halfmax, halves, ind0, iter_stop, KeepConstant, Ll, Ll_iter_best, Lld, Lldd, Lls1, Lls2, Lls3, Lstar, nthreads, ntime, PyrC, R, Rd, Rdd, RddR, RdR, reqrdnum, tform, RiskFail,  RiskPairs, RiskPairs_Strata, Rls1, Rls2, Rls3, Strata_vals, term_n, ties_method, totalnum, TTerm, verbose);
            }
            if (beta_best != beta_c) {  //  if the risk matrices aren't the optimal values, then they must be recalculated
                //  If it goes through every half step without improvement, then the maximum change needs to be decreased
                step_max = step_max*pow(0.5, halfmax);  //  reduces the step sizes
                thres_step_max = thres_step_max*pow(0.5, halfmax);
                //
                beta_p = beta_best;  //
                beta_a = beta_best;  //
                beta_c = beta_best;  //
                for (int ij = 0; ij < totalnum; ij++) {
                    beta_0[ij] = beta_best[ij];
                }
                Pois_Term_Risk_Calc(modelform, tform, term_n, totalnum, fir, dfc, term_tot, T0, Td0, Tdd0, Te, R, Rd, Rdd, Dose, nonDose, beta_0, df0, dint, dslp, TTerm, nonDose_LIN, nonDose_PLIN, nonDose_LOGLIN, RdR, RddR, s_weights, nthreads, KeepConstant, verbose, model_bool, gmix_theta, gmix_term);
            }
            dbeta_max = abs(dbeta[0]);
            for (int ij = 1; ij < totalnum; ij++) {
                if (abs(dbeta[ij]) > dbeta_max) {
                    dbeta_max = abs(dbeta[ij]);
                }
            }
            if (step_max < epsilon) {  //  if the maximum change is too low, then it ends
                iter_stop = 1;
            }
            if (dbeta_max < epsilon) {  //  if the maximum change is too low, then it ends
                iter_stop = 1;
            }
        }
        //  -----------------------------------------------
        //  Performing Full Calculation to get full second derivative matrix
        //  -----------------------------------------------
        fill(Ll.begin(), Ll.end(), 0.0);
        if (!model_bool["single"]) {
            fill(Lld.begin(), Lld.end(), 0.0);
            fill(Lldd.begin(), Lldd.end(), 0.0);
        }
        Pois_Dev_LL_Calc(reqrdnum, totalnum, fir, R, Rd, Rdd, beta_0, RdR, RddR, Ll, Lld, Lldd, PyrC, dev_temp, nthreads, KeepConstant, verbose, model_bool, iter_stop, dev);
        a_n = beta_0;
        beta_fin(guess, _) = a_n;
        LL_fin[guess] = Ll[0];
    }
    fill(Ll.begin(), Ll.end(), 0.0);
    if (!model_bool["single"]) {
        fill(Lld.begin(), Lld.end(), 0.0);
        fill(Lldd.begin(), Lldd.end(), 0.0);
    }
    beta_p = beta_peak;  //
    beta_a = beta_peak;  //
    beta_c = beta_peak;  //
    step_max = step_max0;
    thres_step_max = thres_step_max0;
    iter_stop = 0;
    halves = 0;
    iteration = 0;
    halves = 0;  //  number of half-steps taken
    ind0 = fir;  //  used for validations
    iteration = 0;  //  iteration number
    //
    iter_stop  = 0;  //  tracks if the iterations should be stopped for convergence
    //
    //  next we need to figure out what point to start at
    int best_guess = - 1;
    for (int guess = guesses - 1; guess> - 1; guess--) {
        //  we want to find the closest point at which the loglikelihood is below lstar
        if (LL_fin[guess] < Lstar) {
            best_guess = guess;
        }
    }
    if (verbose >= 3) {
        Rcout << "C++ Note: Upper Guess Results" << endl;
        Rcout << "Guess number, parameter values, Log-Likelihood change" << endl;
        NumericVector beta_temp;
        for (int i = 0;  i < guesses; i++) {
            beta_temp = wrap(beta_fin.row(i));
            if (i == best_guess) {
                Rcout << i << ", " << beta_temp << ", " << LL_fin[i] - Lstar << "<-- Best Guess" << endl;
            } else {
                Rcout << i << ", " << beta_temp << ", " << LL_fin[i] - Lstar << endl;
            }
        }
    }
    NumericVector beta_temp;
    NumericVector beta_temp0;
    if (best_guess == 0) {
        beta_temp = wrap(beta_fin.row(best_guess));  //  the first point was closest, no lower bound
    } else if (best_guess == - 1) {
        beta_temp = wrap(beta_fin.row(guesses - 1));  //  the last point was closest, no upper bound
    } else {
        beta_temp = wrap(beta_fin.row(best_guess));
        beta_temp0 = wrap(beta_fin.row(best_guess - 1));
        for (int i = 0; i < beta_temp.size(); i++) {
            beta_temp[i] = beta_temp[i] + beta_temp0[i];
            beta_temp[i] = beta_temp[i]/2;
        }
    }
    if (verbose >= 4) {
        Rcout << "C++ Note: Initial Guess: ";
        for (int i = 0; i < beta_0.size(); i++) {
            Rcout << beta_temp[i] << " ";
        }
        Rcout << " " << endl;
    }
    for (int i = 0; i < beta_0.size(); i++) {
        a_n[i] = beta_temp[i];
        beta_0[i] = a_n[i];
        beta_c[i] = beta_0[i];
    }
    //
    //  Calculates the subterm and term values
    Pois_Term_Risk_Calc(modelform, tform, term_n, totalnum, fir, dfc, term_tot, T0, Td0, Tdd0, Te, R, Rd, Rdd, Dose, nonDose, beta_0, df0, dint, dslp, TTerm, nonDose_LIN, nonDose_PLIN, nonDose_LOGLIN, RdR, RddR, s_weights, nthreads, KeepConstant, verbose, model_bool, gmix_theta, gmix_term);
    Pois_Dev_LL_Calc(reqrdnum, totalnum, fir, R, Rd, Rdd, beta_0, RdR, RddR, Ll, Lld, Lldd, PyrC, dev_temp, nthreads, KeepConstant, verbose, model_bool, iter_stop, dev);
    //
    #ifdef _OPENMP
    #pragma omp parallel for schedule(dynamic) num_threads(nthreads)
    #endif
    for (int ijk = 0; ijk < reqrdnum*(reqrdnum + 1)/2; ijk++) {
        int ij = 0;
        int jk = ijk;
        while (jk > ij) {
            ij++;
            jk -= ij;
        }
        Lldd_vec[ij * reqrdnum + jk] = Lldd[ij*reqrdnum+jk];
        Lldd_vec[jk * reqrdnum + ij] = Lldd_vec[ij * reqrdnum + jk];
        if (ij == jk) {
            Lld_vecc[ij] = Lld[ij];
        }
    }
    Lldd_vec.attr("dim") = Dimension(reqrdnum, reqrdnum);
    Lldd_mat = as<Map<MatrixXd> >(Lldd_vec);
    Lld_vec = as<Map<VectorXd> >(Lld_vecc);
    int step = 1;
    bool iter_continue = true;
    double max_change = 100;
    ///
    //  variables added for log loop code
    ///
    int bound_val = 1;
    //
    while ((step < maxstep) && (iter_continue)) {
        step++;
        trouble = false;
        Log_Bound(deriv_max, Lldd_mat, Lld_vec, Lstar, qchi, Ll[0], para_number, nthreads, totalnum, reqrdnum, KeepConstant, term_tot, step, dbeta, beta_0, upper, trouble, verbose, mult);
        if (trouble) {
            Calc_Change_trouble(para_number, nthreads, totalnum, thres_step_max, lr, step_max, Ll, Lld, Lldd, dbeta, tform, thres_step_max, step_max, KeepConstant_trouble);
        }
        beta_p = beta_c;  //
        beta_a = beta_c;  //
        for (int ijk = 0; ijk < totalnum; ijk++) {
            if (KeepConstant[ijk] == 0) {
                //
                if ((tform[ijk] == "lin_quad_int") || (tform[ijk] == "lin_exp_int") || (tform[ijk] == "step_int") || (tform[ijk] == "lin_int")) {  //  the threshold values use different maximum deviation values
                    if (abs(dbeta[ijk]) > thres_step_max) {
                        dbeta[ijk] = thres_step_max * sign(dbeta[ijk]);
                    }
                } else {
                    if (abs(dbeta[ijk]) > step_max) {
                        dbeta[ijk] = step_max * sign(dbeta[ijk]);
                    }
                }
            } else {
                dbeta[ijk] = 0;
            }
        }
        max_change = abs(dbeta[0]);
        //
        for (int ij = 0; ij < totalnum; ij++) {
            if (ij == para_number) {
               //  prevent the parameter estimate from crossing the optimum
               //  issue is beta_0[para_number] <= beta_best[para_number]
               if (dbeta[ij] < (beta_peak[para_number] - beta_a[ij])/lr) {
                   dbeta[ij] = (beta_peak[para_number] - beta_a[ij])/lr/2;
               }
            }
            if (abs(dbeta[ij]) > max_change) {
                max_change = abs(dbeta[ij]);
            }
            beta_0[ij] = beta_a[ij] + lr*dbeta[ij];
            beta_c[ij] = beta_0[ij];
        }
        //  ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------//
        //  The same subterm, risk, sides, and log-likelihood calculations are performed every half-step and iteration
        //  ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------//
        Pois_Term_Risk_Calc(modelform, tform, term_n, totalnum, fir, dfc, term_tot, T0, Td0, Tdd0, Te, R, Rd, Rdd, Dose, nonDose, beta_0, df0, dint, dslp, TTerm, nonDose_LIN, nonDose_PLIN, nonDose_LOGLIN, RdR, RddR, s_weights, nthreads, KeepConstant, verbose, model_bool, gmix_theta, gmix_term);
        bound_val = 1;
        Cox_Pois_Log_Loop(step_max, model_bool, beta_0, beta_a, beta_c, bound_val, dbeta, df0, dfc, dint, Dose, thres_step_max, dslp, fir, gmix_term, gmix_theta, half_check, halfmax, KeepConstant, limit_hit, lr, modelform, nonDose, nonDose_LIN, nonDose_LOGLIN, nonDose_PLIN, nthreads, R, Rd, Rdd, RddR, RdR, s_weights, T0, Td0, Tdd0, Te, term_n, term_tot, tform, totalnum, TTerm, verbose);
        if (limit_hit[1]) {
            limits[1] = 0;
            ll_final[1] = 0;
            limit_converged[1] = FALSE;
            break;
        }
        Pois_Dev_LL_Calc(reqrdnum, totalnum, fir, R, Rd, Rdd, beta_0, RdR, RddR, Ll, Lld, Lldd, PyrC, dev_temp, nthreads, KeepConstant, verbose, model_bool, iter_stop, dev);
        //
        #ifdef _OPENMP
        #pragma omp parallel for schedule(dynamic) num_threads(nthreads)
        #endif
        for (int ijk = 0; ijk < reqrdnum*(reqrdnum + 1)/2; ijk++) {
            int ij = 0;
            int jk = ijk;
            while (jk > ij) {
                ij++;
                jk -= ij;
            }
            Lldd_vec[ij * reqrdnum + jk] = Lldd[ij*reqrdnum+jk];
            Lldd_vec[jk * reqrdnum + ij] = Lldd_vec[ij * reqrdnum + jk];
            if (ij == jk) {
                Lld_vecc[ij] = Lld[ij];
            }
        }
        Lldd_vec.attr("dim") = Dimension(reqrdnum, reqrdnum);
        Lldd_mat = as<Map<MatrixXd> >(Lldd_vec);
        Lld_vec = as<Map<VectorXd> >(Lld_vecc);
        limits[1] = beta_0[para_number];
        ll_final[1] = Ll[0];
        if ((max_change < epsilon) && (deriv_max < deriv_epsilon)) {
           iter_continue = false;
           limit_converged[1] = TRUE;
        }
    }
    //  Now refresh matrices back to the maximum point
    beta_p = beta_peak;  //
    beta_a = beta_peak;  //
    beta_c = beta_peak;  //
    step_max = step_max0;
    thres_step_max = thres_step_max0;
    for (int ij = 0; ij < totalnum; ij++) {
        beta_0[ij] = beta_a[ij];
    }
    Pois_Term_Risk_Calc(modelform, tform, term_n, totalnum, fir, dfc, term_tot, T0, Td0, Tdd0, Te, R, Rd, Rdd, Dose, nonDose, beta_0, df0, dint, dslp, TTerm, nonDose_LIN, nonDose_PLIN, nonDose_LOGLIN, RdR, RddR, s_weights, nthreads, KeepConstant, verbose, model_bool, gmix_theta, gmix_term);
    Pois_Dev_LL_Calc(reqrdnum, totalnum, fir, R, Rd, Rdd, beta_0, RdR, RddR, Ll, Lld, Lldd, PyrC, dev_temp, nthreads, KeepConstant, verbose, model_bool, iter_stop, dev);
    //
    upper = false;
    #ifdef _OPENMP
    #pragma omp parallel for schedule(dynamic) num_threads(nthreads)
    #endif
    for (int ijk = 0; ijk < reqrdnum*(reqrdnum + 1)/2; ijk++) {
        int ij = 0;
        int jk = ijk;
        while (jk > ij) {
            ij++;
            jk -= ij;
        }
        Lldd_vec[ij * reqrdnum + jk] = Lldd[ij*reqrdnum+jk];
        Lldd_vec[jk * reqrdnum + ij] = Lldd_vec[ij * reqrdnum + jk];
        if (ij == jk) {
            Lld_vecc[ij] = Lld[ij];
        }
    }
    Lldd_vec.attr("dim") = Dimension(reqrdnum, reqrdnum);
    Lldd_mat = as<Map<MatrixXd> >(Lldd_vec);
    Lld_vec = as<Map<VectorXd> >(Lld_vecc);
    //
    trouble = false;
    Log_Bound(deriv_max, Lldd_mat, Lld_vec, Lstar, qchi, Lmax, para_number, nthreads, totalnum, reqrdnum, KeepConstant, term_tot, 0, dbeta_start, beta_0, upper, trouble, verbose, mult);
    //  now the dbeta holds the range to check
    //  note that currently the guesses are not bounded
    for (int i = 0; i < guesses; i++) {
        for (int j = 0; j < totalnum; j++) {
            //  use dbeta to assign a_n guess i, parameter j
            //  assign such that i=guesses - 1 gives mult*dbeta
            a_ns(i, j) = beta_0[j] + dbeta_start[j] * (i + 1)/(guesses);
        }
    }
    if (verbose >= 4) {
        for (int i = 0; i < guesses; i++) {
            Rcout << "C++ Note: Initial guess " << i << ": ";
            for (int j = 0; j < totalnum; j++) {
                Rcout << a_ns(i, j) << " ";
            }
            Rcout << " " << endl;
        }
    }
    //  now we have the points to test
    halves = 0;  //  number of half-steps taken
    ind0 = fir;  //  used for validations
    iteration = 0;  //  iteration number
    maxiter = 0;
    iter_stop  = 0;  //  tracks if the iterations should be stopped for convergence
    Ll_iter_best = 10;
    for (int guess = 0; guess <guesses; guess++) {
        fill(Ll.begin(), Ll.end(), 0.0);
        fill(Lld.begin(), Lld.end(), 0.0);
        fill(Lldd.begin(), Lldd.end(), 0.0);
        beta_p = beta_peak;  //
        beta_a = beta_peak;  //
        beta_c = beta_peak;  //
        step_max = step_max0;
        thres_step_max = thres_step_max0;
        iter_stop = 0;
        halves = 0;
        iteration = 0;
        halves = 0;  //  number of half-steps taken
        ind0 = fir;  //  used for validations
        iteration = 0;  //  iteration number
        Ll_iter_best = 10;
        iter_stop  = 0;  //  tracks if the iterations should be stopped for convergence
        maxiter = maxiters[0];
        a_n = a_ns.row(guess);
        for (int i = 0; i < beta_0.size(); i++) {
            beta_0[i] = a_n[i];
        }
        Pois_Term_Risk_Calc(modelform, tform, term_n, totalnum, fir, dfc, term_tot, T0, Td0, Tdd0, Te, R, Rd, Rdd, Dose, nonDose, beta_0, df0, dint, dslp, TTerm, nonDose_LIN, nonDose_PLIN, nonDose_LOGLIN, RdR, RddR, s_weights, nthreads, KeepConstant, verbose, model_bool, gmix_theta, gmix_term);
        if ((R.minCoeff() <= 0) || (R.hasNaN())) {
            iter_stop = 1;
            Ll[0] = 404;  //  note to consider this point too far
        } else {
            Pois_Dev_LL_Calc(reqrdnum, totalnum, fir, R, Rd, Rdd, beta_0, RdR, RddR, Ll, Lld, Lldd, PyrC, dev_temp, nthreads, KeepConstant, verbose, model_bool, iter_stop, dev);
            for (int i = 0; i < beta_0.size(); i++) {
                beta_c[i] = beta_0[i];
            }
        }
        while ((iteration < maxiter) && (iter_stop == 0)) {
            iteration++;
            beta_p = beta_c;  //
            beta_a = beta_c;  //
            beta_best = beta_c;  //
            //  calculates the initial change in parameter
            Calc_Change_trouble(para_number, nthreads, totalnum, thres_step_max, lr, step_max, Ll, Lld, Lldd, dbeta, tform, thres_step_max, step_max, KeepConstant_trouble);
            Intercept_Bound(nthreads, totalnum, beta_0, dbeta, dfc, df0, KeepConstant, tform);
            if ((Ll_iter_best > 0) || (Ll_iter_best < Ll[ind0])) {
                Ll_iter_best = Ll[ind0];
                beta_abs_best = beta_c;
            }
            halves = 0;
            while ((Ll[ind0] <= Ll_iter_best) && (halves < halfmax)) {  //  repeats until half-steps maxed or an improvement
                for (int ij = 0; ij < totalnum; ij++) {
                    beta_0[ij] = beta_a[ij] + dbeta[ij];
                    beta_c[ij] = beta_0[ij];
                }
                //  ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------//
                //  The same subterm, risk, sides, and log-likelihood calculations are performed every half-step and iteration
                //  ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------//
                Pois_Term_Risk_Calc(modelform, tform, term_n, totalnum, fir, dfc, term_tot, T0, Td0, Tdd0, Te, R, Rd, Rdd, Dose, nonDose, beta_0, df0, dint, dslp, TTerm, nonDose_LIN, nonDose_PLIN, nonDose_LOGLIN, RdR, RddR, s_weights, nthreads, KeepConstant, verbose, model_bool, gmix_theta, gmix_term);
                Cox_Pois_Check_Continue(model_bool, beta_0, beta_best, beta_c, cens_weight, dbeta, dev, dev_temp, fir, halfmax, halves, ind0, iter_stop, KeepConstant, Ll, Ll_iter_best, Lld, Lldd, Lls1, Lls2, Lls3, Lstar, nthreads, ntime, PyrC, R, Rd, Rdd, RddR, RdR, reqrdnum, tform, RiskFail,  RiskPairs, RiskPairs_Strata, Rls1, Rls2, Rls3, Strata_vals, term_n, ties_method, totalnum, TTerm, verbose);
            }
            if (beta_best != beta_c) {  //  if the risk matrices aren't the optimal values, then they must be recalculated
                //  If it goes through every half step without improvement, then the maximum change needs to be decreased
                step_max = step_max*pow(0.5, halfmax);  //  reduces the step sizes
                thres_step_max = thres_step_max*pow(0.5, halfmax);
//                iter_check = 1;
                beta_p = beta_best;  //
                beta_a = beta_best;  //
                beta_c = beta_best;  //
                for (int ij = 0; ij < totalnum; ij++) {
                    beta_0[ij] = beta_best[ij];
                }
                Pois_Term_Risk_Calc(modelform, tform, term_n, totalnum, fir, dfc, term_tot, T0, Td0, Tdd0, Te, R, Rd, Rdd, Dose, nonDose, beta_0, df0, dint, dslp, TTerm, nonDose_LIN, nonDose_PLIN, nonDose_LOGLIN, RdR, RddR, s_weights, nthreads, KeepConstant, verbose, model_bool, gmix_theta, gmix_term);
            }
            dbeta_max = abs(dbeta[0]);
            for (int ij = 1; ij < totalnum; ij++) {
                if (abs(dbeta[ij]) > dbeta_max) {
                    dbeta_max = abs(dbeta[ij]);
                }
            }
            if (step_max < epsilon) {  //  if the maximum change is too low, then it ends
                iter_stop = 1;
            }
            if (dbeta_max < epsilon) {  //  if the maximum change is too low, then it ends
                iter_stop = 1;
            }
        }
        //  -----------------------------------------------
        //  Performing Full Calculation to get full second derivative matrix
        //  -----------------------------------------------
        fill(Ll.begin(), Ll.end(), 0.0);
        if (!model_bool["single"]) {
            fill(Lld.begin(), Lld.end(), 0.0);
            fill(Lldd.begin(), Lldd.end(), 0.0);
        }
        Pois_Dev_LL_Calc(reqrdnum, totalnum, fir, R, Rd, Rdd, beta_0, RdR, RddR, Ll, Lld, Lldd, PyrC, dev_temp, nthreads, KeepConstant, verbose, model_bool, iter_stop, dev);
        a_n = beta_0;
        beta_fin(guess, _) = a_n;
        LL_fin[guess] = Ll[0];
    }
    fill(Ll.begin(), Ll.end(), 0.0);
    if (!model_bool["single"]) {
        fill(Lld.begin(), Lld.end(), 0.0);
        fill(Lldd.begin(), Lldd.end(), 0.0);
    }
    beta_p = beta_best;  //
    beta_a = beta_best;  //
    beta_c = beta_best;  //
    step_max = step_max0;
    thres_step_max = thres_step_max0;
    iter_stop = 0;
    halves = 0;
    iteration = 0;
    halves = 0;  //  number of half-steps taken
    ind0 = fir;  //  used for validations
    iteration = 0;  //  iteration number
    iter_stop  = 0;  //  tracks if the iterations should be stopped for convergence
    //  next we need to figure out what point to start at
    best_guess = - 1;
    for (int guess = guesses - 1; guess> - 1; guess--) {
        //  we want to find the closest point at which the loglikelihood is below lstar
        if (LL_fin[guess] < Lstar) {
            best_guess = guess;
        }
    }
    if (verbose >= 3) {
        Rcout << "C++ Note: Lower Guess Results" << endl;
        Rcout << "Guess number, parameter values, Log-Likelihood change" << endl;
        NumericVector beta_temp;
        for (int i = 0;  i < guesses; i++) {
            beta_temp = wrap(beta_fin.row(i));
            if (i == best_guess) {
                Rcout << i << ", " << beta_temp << ", " << LL_fin[i] - Lstar << "<-- Best Guess" << endl;
            } else {
                Rcout << i << ", " << beta_temp << ", " << LL_fin[i] - Lstar << endl;
            }
        }
    }
    if (best_guess == 0) {
        beta_temp = wrap(beta_fin.row(best_guess));  //  the first point was closest, no lower bound
    } else if (best_guess == - 1) {
        beta_temp = wrap(beta_fin.row(guesses - 1));  //  the last point was closest, no upper bound
    } else {
        beta_temp = wrap(beta_fin.row(best_guess));
        beta_temp0 = wrap(beta_fin.row(best_guess - 1));
        for (int i = 0; i < beta_temp.size(); i++) {
            beta_temp[i] = beta_temp[i] + beta_temp0[i];
            beta_temp[i] = beta_temp[i]/2;
        }
    }
    if (verbose >= 4) {
        Rcout << "C++ Note: Initial Guess: ";
        for (int i = 0; i < beta_0.size(); i++) {
            Rcout << beta_temp[i] << " ";
        }
        Rcout << " " << endl;
    }
    for (int i = 0; i < beta_0.size(); i++) {
        a_n[i] = beta_temp[i];
        beta_0[i] = a_n[i];
        beta_c[i] = beta_0[i];
    }
    //  Calculates the subterm and term values
    Pois_Term_Risk_Calc(modelform, tform, term_n, totalnum, fir, dfc, term_tot, T0, Td0, Tdd0, Te, R, Rd, Rdd, Dose, nonDose, beta_0, df0, dint, dslp, TTerm, nonDose_LIN, nonDose_PLIN, nonDose_LOGLIN, RdR, RddR, s_weights, nthreads, KeepConstant, verbose, model_bool, gmix_theta, gmix_term);
    Pois_Dev_LL_Calc(reqrdnum, totalnum, fir, R, Rd, Rdd, beta_0, RdR, RddR, Ll, Lld, Lldd, PyrC, dev_temp, nthreads, KeepConstant, verbose, model_bool, iter_stop, dev);
    //
    #ifdef _OPENMP
    #pragma omp parallel for schedule(dynamic) num_threads(nthreads)
    #endif
    for (int ijk = 0; ijk < reqrdnum*(reqrdnum + 1)/2; ijk++) {
        int ij = 0;
        int jk = ijk;
        while (jk > ij) {
            ij++;
            jk -= ij;
        }
        Lldd_vec[ij * reqrdnum + jk] = Lldd[ij*reqrdnum+jk];
        Lldd_vec[jk * reqrdnum + ij] = Lldd_vec[ij * reqrdnum + jk];
        if (ij == jk) {
            Lld_vecc[ij] = Lld[ij];
        }
    }
    Lldd_vec.attr("dim") = Dimension(reqrdnum, reqrdnum);
    Lldd_mat = as<Map<MatrixXd> >(Lldd_vec);
    Lld_vec = as<Map<VectorXd> >(Lld_vecc);
    step = 1;
    iter_continue = true;
    max_change = 100;
    while ((step < maxstep) && (iter_continue)) {
        step++;
        trouble = false;
        Log_Bound(deriv_max, Lldd_mat, Lld_vec, Lstar, qchi, Ll[0], para_number, nthreads, totalnum, reqrdnum, KeepConstant, term_tot, step, dbeta, beta_0, upper, trouble, verbose, mult);
        if (trouble) {
            Calc_Change_trouble(para_number, nthreads, totalnum, thres_step_max, lr, step_max, Ll, Lld, Lldd, dbeta, tform, thres_step_max, step_max, KeepConstant_trouble);
        }
        //
        beta_p = beta_c;  //
        beta_a = beta_c;  //
        //
        for (int ijk = 0; ijk < totalnum; ijk++) {
            if (KeepConstant[ijk] == 0) {
                //
                if ((tform[ijk] == "lin_quad_int") || (tform[ijk] == "lin_exp_int") || (tform[ijk] == "step_int") || (tform[ijk] == "lin_int")) {  //  the threshold values use different maximum deviation values
                    if (abs(dbeta[ijk]) > thres_step_max) {
                        dbeta[ijk] = thres_step_max * sign(dbeta[ijk]);
                    }
                } else {
                    if (abs(dbeta[ijk]) > step_max) {
                        dbeta[ijk] = step_max * sign(dbeta[ijk]);
                    }
                }
            } else {
                dbeta[ijk] = 0;
            }
        }
        //
        max_change = abs(dbeta[0]);
        for (int ij = 0; ij < totalnum; ij++) {
            if (ij == para_number) {
                //  prevent the parameter estimate from crossing the optimum
                //  issue is beta_0[para_number] <= beta_best[para_number]
                if (dbeta[ij] > (beta_peak[para_number] - beta_a[ij])/lr) {
                    dbeta[ij] = (beta_peak[para_number] - beta_a[ij])/lr/2;
                }
            }
            if (abs(dbeta[ij]) > max_change) {
                max_change = abs(dbeta[ij]);
            }
            beta_0[ij] = beta_a[ij] + lr*dbeta[ij];
            beta_c[ij] = beta_0[ij];
        }
        //  ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------//
        //  The same subterm, risk, sides, and log-likelihood calculations are performed every half-step and iteration
        //  ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------//
        Pois_Term_Risk_Calc(modelform, tform, term_n, totalnum, fir, dfc, term_tot, T0, Td0, Tdd0, Te, R, Rd, Rdd, Dose, nonDose, beta_0, df0, dint, dslp, TTerm, nonDose_LIN, nonDose_PLIN, nonDose_LOGLIN, RdR, RddR, s_weights, nthreads, KeepConstant, verbose, model_bool, gmix_theta, gmix_term);
        bound_val = 0;
        Cox_Pois_Log_Loop(step_max, model_bool, beta_0, beta_a, beta_c, bound_val, dbeta, df0, dfc, dint, Dose, thres_step_max, dslp, fir, gmix_term, gmix_theta, half_check, halfmax, KeepConstant, limit_hit, lr, modelform, nonDose, nonDose_LIN, nonDose_LOGLIN, nonDose_PLIN, nthreads, R, Rd, Rdd, RddR, RdR, s_weights, T0, Td0, Tdd0, Te, term_n, term_tot, tform, totalnum, TTerm, verbose);
        if (limit_hit[0]) {
            limits[0] = 0;
            ll_final[0] = 0;
            limit_converged[0] = FALSE;
            break;
        }
        Pois_Dev_LL_Calc(reqrdnum, totalnum, fir, R, Rd, Rdd, beta_0, RdR, RddR, Ll, Lld, Lldd, PyrC, dev_temp, nthreads, KeepConstant, verbose, model_bool, iter_stop, dev);
        //
        #ifdef _OPENMP
        #pragma omp parallel for schedule(dynamic) num_threads(nthreads)
        #endif
        for (int ijk = 0; ijk < reqrdnum*(reqrdnum + 1)/2; ijk++) {
            int ij = 0;
            int jk = ijk;
            while (jk > ij) {
                ij++;
                jk -= ij;
            }
            Lldd_vec[ij * reqrdnum + jk] = Lldd[ij*reqrdnum+jk];
            Lldd_vec[jk * reqrdnum + ij] = Lldd_vec[ij * reqrdnum + jk];
            if (ij == jk) {
                Lld_vecc[ij] = Lld[ij];
            }
        }
        Lldd_vec.attr("dim") = Dimension(reqrdnum, reqrdnum);
        Lldd_mat = as<Map<MatrixXd> >(Lldd_vec);
        Lld_vec = as<Map<VectorXd> >(Lld_vecc);
        limits[0] = beta_0[para_number];
        ll_final[0] = Ll[0];
        if ((max_change < epsilon) && (deriv_max < deriv_epsilon)) {
           iter_continue = false;
           limit_converged[0] = TRUE;
        }
    }
    //
    res_list = List::create(_["Parameter_Limits"] = wrap(limits), _["Negative_Risk_Limit_Hit"] = wrap(limit_hit), _["Likelihood_Boundary"] = wrap(ll_final), _["Likelihood_Goal"] = wrap(Lstar), _["Limit_Converged"] = wrap(limit_converged), _["Status"] = "PASSED");
    //  returns a list of results
    return res_list;
}

//' Primary Cox PH likelihood bounds calcualtion function.
//'
//' \code{LogLik_Cox_PH_Omnibus_Log_Bound} Performs the calls to calculation functions and log-likeihood profile bounds
//'
//' @inheritParams CPP_template
//'
//' @return List of final results: Log-likelihood of optimum, first derivative of log-likelihood, second derivative matrix, parameter list, standard deviation estimate, AIC, model information
//' @noRd
//'
//
List LogLik_Cox_PH_Omnibus_Log_Bound_CurveSearch(IntegerVector term_n, StringVector tform, Ref<VectorXd> beta_0, Ref<MatrixXd> df0, IntegerVector dfc, int fir, string modelform, double lr, List optim_para, int maxiter, int halfmax, double epsilon, double step_max, double thres_step_max, double deriv_epsilon, const Ref<const MatrixXd>& df_m, NumericVector tu, int verbose, IntegerVector KeepConstant, int term_tot, string ties_method, int nthreads, NumericVector& Strata_vals, const VectorXd& cens_weight, List model_bool, const double gmix_theta, const IntegerVector gmix_term, const Ref<const MatrixXd>& Lin_Sys, const Ref<const VectorXd>& Lin_Res, double qchi, int para_number, int maxstep, double step_size) {
    //
    List temp_list = List::create(_["Status"] = "TEMP");  //  used as a dummy return value for code checking
    //  Time durations are measured from this point on in microseconds
    //
    //  df0: covariate data
    //  ntime: number of event times for Cox PH
    //  totalnum: number of terms used
    //
    //  ------------------------------------------------------------------------- //  initialize
    int ntime = tu.size();
    int totalnum = term_n.size();
    int reqrdnum = totalnum - sum(KeepConstant);
    //  ------------------------------------------------------------------------- //  initialize
    if (model_bool["null"]) {
        if (verbose >= 1) {
            Rcout << "null model is not compatable with log-based bound calculation" << endl;
        }
        temp_list = List::create(_["Status"] = "FAILED_BAD_MODEL_NULL", _["LogLik"] = R_NaN);
        return temp_list;
    }
    if (model_bool["constraint"]) {
        if (verbose >= 1) {
            Rcout << "linear constataints are not compatable with Case-Control model calculation" << endl;
        }
        temp_list = List::create(_["Status"] = "FAILED_WITH_BAD_MODEL_CONSTRAINT", _["LogLik"] = R_NaN);
        return temp_list;
    }
    //
    //  cout.precision: controls the number of significant digits printed
    //  nthreads: number of threads used for parallel operations
    //
    Rcout.precision(7);  //  forces higher precision numbers printed to terminal
    //
    //  Lld_worst: stores the highest magnitude log-likelihood derivative
    //  ---------------------------------------------
    //  To Start, needs to seperate the derivative terms
    //  ---------------------------------------------
    //
    //  ------------------------------------------------------------------------- //  initialize
    // Map<VectorXd> beta_0(as<Map<VectorXd> >(a_n));
    MatrixXd T0;
    MatrixXd Td0;
    MatrixXd Tdd0;
    //
    MatrixXd Te;
    MatrixXd R;
    MatrixXd Rd;
    MatrixXd Rdd;
    //
    MatrixXd Dose;
    MatrixXd nonDose;
    MatrixXd nonDose_LIN;
    MatrixXd nonDose_PLIN;
    MatrixXd nonDose_LOGLIN;
    MatrixXd TTerm;
    double dint = 0.0;  //  the amount of change used to calculate derivatives in threshold paramters
    double dslp = 0.0;
    MatrixXd RdR;
    MatrixXd RddR;
    //  ------------------------------------------------------------------------- //  initialize
    //  ---------------------------------------------
    //  To Start, needs to seperate the derivative terms
    //  ---------------------------------------------
    //
    Cox_Refresh_R_TERM(totalnum, reqrdnum, term_tot, dint, dslp, thres_step_max, step_max, df0, T0, Td0, Tdd0, Te, R, Rd, Rdd, Dose, nonDose, TTerm, nonDose_LIN, nonDose_PLIN, nonDose_LOGLIN, RdR, RddR, model_bool);
    //  ------------------------------------------------------------------------- //  initialize
    IntegerMatrix RiskFail;
    vector<vector<int> > RiskPairs(ntime);
    vector<vector<vector<int> > > RiskPairs_Strata(ntime, vector<vector<int>>(Strata_vals.size()));
    //  ------------------------------------------------------------------------- //  initialize
    if (model_bool["strata"]) {
        RiskFail = IntegerMatrix(ntime, 2*Strata_vals.size());  //  vector giving the event rows
        //
        //  Creates matrices used to identify the event risk groups
        if (model_bool["cr"]) {
            Make_Groups_Strata_CR(ntime, df_m, RiskFail, RiskPairs_Strata, tu, nthreads, Strata_vals, cens_weight);
        } else {
            Make_Groups_Strata(ntime, df_m, RiskFail, RiskPairs_Strata, tu, nthreads, Strata_vals);
        }
    } else {
        RiskFail = IntegerMatrix(ntime, 2);  //  vector giving the event rows
        //
        //  Creates matrices used to identify the event risk groups
        if (model_bool["cr"]) {
            Make_Groups_CR(ntime, df_m, RiskFail, RiskPairs, tu, cens_weight, nthreads);
        } else {
            Make_Groups(ntime, df_m, RiskFail, RiskPairs, tu, nthreads);
        }
    }
    //  ------------------------------------------------------------------------- //  initialize
    MatrixXd Rls1;
    MatrixXd Lls1;
    MatrixXd Rls2;
    MatrixXd Rls3;
    MatrixXd Lls2;
    MatrixXd Lls3;
    vector<double> Ll(reqrdnum, 0.0);  //  log-likelihood values
    vector<double> Lld(reqrdnum, 0.0);  //  log-likelihood derivative values
    vector<double> Lldd(pow(reqrdnum, 2), 0.0);  //  the second derivative matrix has room for every combination, but only the lower triangle is calculated initially
    //  ------------------------------------------------------------------------- //  initialize
    Cox_Refresh_R_SIDES(reqrdnum, ntime, Rls1, Rls2, Rls3, Lls1, Lls2, Lls3, Strata_vals, model_bool);
    //  the log-likelihood is calculated in parallel over the risk groups
    vector <double> Ll_comp(2, Ll[0]);  //  vector to compare values
    double step_max0 = step_max;
    double thres_step_max0 = thres_step_max;
    //
    vector<double> dbeta(totalnum, 0.0);
    //
    //  --------------------------
    //  always starts from initial guess
    //  --------------------------
    vector<double> beta_c(totalnum, 0.0);
    vector<double> beta_a(totalnum, 0.0);
    vector<double> beta_best(totalnum, 0.0);
    vector<double> beta_peak(totalnum, 0.0);
    vector<double> beta_p(totalnum, 0.0);
    VectorXd::Map(&beta_p[0], beta_0.size()) = beta_0;  //  stores previous parameters
    VectorXd::Map(&beta_c[0], beta_0.size()) = beta_0;  //  stores current parameters
    VectorXd::Map(&beta_a[0], beta_0.size()) = beta_0;  //  stores a refrence value for parameters
    VectorXd::Map(&beta_best[0], beta_0.size()) = beta_0;  //  stores the best parameters
    VectorXd::Map(&beta_peak[0], beta_0.size()) = beta_0;  //  stores the best parameters
    int iter_stop  = 0;  //  tracks if the iterations should be stopped for convergence
    //
    //
    vector<double> beta_abs_best(totalnum, 0.0);
    //
    Cox_Refresh_R_SIDES(reqrdnum, ntime, Rls1, Rls2, Rls3, Lls1, Lls2, Lls3, Strata_vals, model_bool);
    fill(Ll.begin(), Ll.end(), 0.0);
    fill(Lld.begin(), Lld.end(), 0.0);
    fill(Lldd.begin(), Lldd.end(), 0.0);
    beta_p = beta_peak;  //
    beta_a = beta_peak;  //
    beta_c = beta_peak;  //
    step_max = step_max0;
    thres_step_max = thres_step_max0;
    //
    for (int i = 0; i < beta_0.size(); i++) {
        beta_0[i] = beta_peak[i];
    }
    Cox_Term_Risk_Calc(modelform, tform, term_n, totalnum, fir, dfc, term_tot, T0, Td0, Tdd0, Te, R, Rd, Rdd, Dose, nonDose, beta_0, df0, thres_step_max, step_max, TTerm, nonDose_LIN, nonDose_PLIN, nonDose_LOGLIN, RdR, RddR, nthreads, KeepConstant, verbose, model_bool, gmix_theta, gmix_term);
    Cox_Side_LL_Calc(reqrdnum, ntime, tform, RiskFail,  RiskPairs, RiskPairs_Strata, totalnum, fir, R, Rd, Rdd, Rls1, Rls2, Rls3, Lls1, Lls2, Lls3, cens_weight, Strata_vals, beta_0, RdR, RddR, Ll, Lld, Lldd, nthreads, KeepConstant, ties_method, verbose, model_bool, iter_stop);
    //
    double Lstar = Ll[0]-qchi;
    double Lpeak = Ll[0];
    IntegerVector KeepConstant_trouble(totalnum);
    for (int i = 0;  i < totalnum; i++) {
        KeepConstant_trouble[i] = KeepConstant[i];
    }
    KeepConstant_trouble[para_number] = 1;
    //
    vector<double> limits(2, 0.0);
    vector<bool>   limit_hit(2, FALSE);
    vector<bool>   limit_converged(2, FALSE);
    vector<double> ll_final(2, 0.0);
    vector<double> width_final(2, 0.0);
    vector<int>    step_final(2, 0.0);
    List res_list;
    //
    if (verbose >= 4) {
        Rcout << "C++ Note: STARTING Upper Bound" << endl;
    }
    bool convgd = false;
    ///
    //  variables added for log loop code
    ///
    VectorXd s_weights(1);
    //  We need the values reserved for the upper, middle, lower estimates and scores
    vector<double> beta_L(totalnum, 0.0);
    vector<double> beta_M(totalnum, 0.0);
    vector<double> beta_H(totalnum, 0.0);
    double L_L = 0.0;
    double L_M = 0.0;
    double L_H = 0.0;
    NumericVector reg_beta(totalnum);
    //  First we need to establish the first interval estimates
    for (int ij = 0; ij < totalnum; ij++) {
        beta_L[ij] = beta_peak[ij];
        beta_H[ij] = beta_peak[ij];
    }
    L_L = Lpeak;
    bool loop_check = true;
    double temp_step = step_size;
    NumericVector temp_L(1);
    List reg_out;
    while ((loop_check) && (temp_step > 1e-3)) {
        //  assign new high point
        beta_H[para_number] = beta_L[para_number] + temp_step;
        for (int ij = 0; ij < totalnum; ij++) {
            beta_0[ij] = beta_H[ij];
        }
        temp_step = temp_step * 0.5;
        //
        reg_out = Cox_Full_Run(reqrdnum, ntime, tform, RiskFail,  RiskPairs, RiskPairs_Strata, totalnum, fir, R, Rd, Rdd, Rls1, Rls2, Rls3, Lls1, Lls2, Lls3, cens_weight, Strata_vals, beta_0, RdR, RddR, Ll, Lld, Lldd, nthreads, KeepConstant, ties_method, verbose, model_bool, iter_stop, term_tot, dint, dslp, thres_step_max, step_max, df0, T0, Td0, Tdd0, Te, Dose, nonDose, TTerm, nonDose_LIN, nonDose_PLIN, nonDose_LOGLIN, modelform, gmix_theta, gmix_term, convgd, lr, optim_para, maxiter, Lin_Sys, Lin_Res, term_n, dfc, halfmax, epsilon, deriv_epsilon);
        temp_L[0] = reg_out["LogLik"];
        if (!is_nan(temp_L)[0]) {
            loop_check = false;
        }
    }
    if (loop_check) {
        limit_hit[1] = true;
        limits[1] = 0;
        ll_final[1] = 0;
        limit_converged[1] = false;
    } else {
        //  Now we can run the actual algorithm
        reg_beta = reg_out["beta_0"];
        for (int ij = 0; ij < totalnum; ij++) {
            beta_H[ij] = reg_beta[ij];
        }
        L_H = reg_out["LogLik"];
        //  now set the mid point value
        for (int ij = 0; ij < totalnum; ij++) {
            beta_M[ij] = (beta_H[ij] + beta_L[ij])/2;
            beta_0[ij] = beta_M[ij];
        }
        reg_out = Cox_Full_Run(reqrdnum, ntime, tform, RiskFail,  RiskPairs, RiskPairs_Strata, totalnum, fir, R, Rd, Rdd, Rls1, Rls2, Rls3, Lls1, Lls2, Lls3, cens_weight, Strata_vals, beta_0, RdR, RddR, Ll, Lld, Lldd, nthreads, KeepConstant, ties_method, verbose, model_bool, iter_stop, term_tot, dint, dslp, thres_step_max, step_max, df0, T0, Td0, Tdd0, Te, Dose, nonDose, TTerm, nonDose_LIN, nonDose_PLIN, nonDose_LOGLIN, modelform, gmix_theta, gmix_term, convgd, lr, optim_para, maxiter, Lin_Sys, Lin_Res, term_n, dfc, halfmax, epsilon, deriv_epsilon);
        reg_beta = reg_out["beta_0"];
        for (int ij = 0; ij < totalnum; ij++) {
            beta_M[ij] = reg_beta[ij];
        }
        L_M = reg_out["LogLik"];
        //
        int step = 0;
        //  now run the bisection until stopping point
        //  while ((step < step_limit) & (abs(beta_low[para_num] - beta_high[para_num]) > control$epsilon) & (!Limit_Hit[2])) {
        while ((step < maxstep) && (abs(beta_L[para_number] - beta_H[para_number]) > epsilon) && (!limit_hit[1])) {
            step = step + 1;
            if (L_L < Lstar) {
                throw invalid_argument("The lower estimate is too high?");
            } else if (L_H < Lstar) {
                //  midpoint is in between the two
                if (L_M < Lstar) {
                    //  the mid point is past the optimum
                    for (int ij = 0; ij < totalnum; ij++) {
                        beta_H[ij] = beta_M[ij];
                    }
                    L_H = L_M;
                } else {
                    //  the mid point is before the optimum
                    for (int ij = 0; ij < totalnum; ij++) {
                        beta_L[ij] = beta_M[ij];
                    }
                    L_L = L_M;
                }
            } else if (L_M < Lstar) {
                for (int ij = 0; ij < totalnum; ij++) {
                    beta_H[ij] = beta_M[ij];
                }
                L_H = L_M;
            } else {
                //  the upper estimate needs to be shifted up
                for (int ij = 0; ij < totalnum; ij++) {
                    beta_L[ij] = beta_H[ij];
                }
                L_L = L_H;
                //  check new high point
                loop_check = true;
                temp_step = step_size;
                while ((loop_check) && (temp_step > 1e-3)) {
                    //  assign new high point
                    beta_H[para_number] = beta_L[para_number] + temp_step;
                    for (int ij = 0; ij < totalnum; ij++) {
                        beta_0[ij] = beta_H[ij];
                    }
                    temp_step = temp_step * 0.5;
                    //
                    reg_out = Cox_Full_Run(reqrdnum, ntime, tform, RiskFail,  RiskPairs, RiskPairs_Strata, totalnum, fir, R, Rd, Rdd, Rls1, Rls2, Rls3, Lls1, Lls2, Lls3, cens_weight, Strata_vals, beta_0, RdR, RddR, Ll, Lld, Lldd, nthreads, KeepConstant, ties_method, verbose, model_bool, iter_stop, term_tot, dint, dslp, thres_step_max, step_max, df0, T0, Td0, Tdd0, Te, Dose, nonDose, TTerm, nonDose_LIN, nonDose_PLIN, nonDose_LOGLIN, modelform, gmix_theta, gmix_term, convgd, lr, optim_para, maxiter, Lin_Sys, Lin_Res, term_n, dfc, halfmax, epsilon, deriv_epsilon);
                    temp_L[0] = reg_out["LogLik"];
                    if (!is_nan(temp_L)[0]) {
                        loop_check = false;
                    }
                }
                if (loop_check) {
                    limit_hit[1] = true;
                    limits[1] = 0;
                    ll_final[1] = 0;
                    limit_converged[1] = false;
                } else {
                    reg_beta = reg_out["beta_0"];
                    for (int ij = 0; ij < totalnum; ij++) {
                        beta_H[ij] = reg_beta[ij];
                    }
                    L_H = reg_out["LogLik"];
                }
            }
            if (!limit_hit[1]) {
                //  now set the mid point value
                for (int ij = 0; ij < totalnum; ij++) {
                    beta_M[ij] = (beta_H[ij] + beta_L[ij])/2;
                    beta_0[ij] = beta_M[ij];
                }
                reg_out = Cox_Full_Run(reqrdnum, ntime, tform, RiskFail,  RiskPairs, RiskPairs_Strata, totalnum, fir, R, Rd, Rdd, Rls1, Rls2, Rls3, Lls1, Lls2, Lls3, cens_weight, Strata_vals, beta_0, RdR, RddR, Ll, Lld, Lldd, nthreads, KeepConstant, ties_method, verbose, model_bool, iter_stop, term_tot, dint, dslp, thres_step_max, step_max, df0, T0, Td0, Tdd0, Te, Dose, nonDose, TTerm, nonDose_LIN, nonDose_PLIN, nonDose_LOGLIN, modelform, gmix_theta, gmix_term, convgd, lr, optim_para, maxiter, Lin_Sys, Lin_Res, term_n, dfc, halfmax, epsilon, deriv_epsilon);
                reg_beta = reg_out["beta_0"];
                for (int ij = 0; ij < totalnum; ij++) {
                    beta_M[ij] = reg_beta[ij];
                }
                L_M = reg_out["LogLik"];
            }
        }
        if ((abs(beta_L[para_number] - beta_H[para_number]) < epsilon) && (!limit_hit[1])) {
            limit_converged[1] = true;
        }
        limits[1] = beta_M[para_number];
        ll_final[1] = L_M;
        width_final[1] = abs(beta_L[para_number] - beta_H[para_number]);
        step_final[1] = step;
    }
    //  upper limit found, now solve lower limit
    if (verbose >= 4) {
        Rcout << "C++ Note: STARTING Lower Bound" << endl;
    }
    for (int ij = 0; ij < totalnum; ij++) {
        beta_L[ij] = beta_peak[ij];
        beta_H[ij] = beta_peak[ij];
    }
    L_H = Lpeak;
    loop_check = true;
    temp_step = step_size;
    while ((loop_check) && (temp_step > 1e-3)) {
        //  assign new high point
        beta_L[para_number] = beta_H[para_number] - temp_step;
        for (int ij = 0; ij < totalnum; ij++) {
            beta_0[ij] = beta_L[ij];
        }
        temp_step = temp_step * 0.5;
        //
        reg_out = Cox_Full_Run(reqrdnum, ntime, tform, RiskFail,  RiskPairs, RiskPairs_Strata, totalnum, fir, R, Rd, Rdd, Rls1, Rls2, Rls3, Lls1, Lls2, Lls3, cens_weight, Strata_vals, beta_0, RdR, RddR, Ll, Lld, Lldd, nthreads, KeepConstant, ties_method, verbose, model_bool, iter_stop, term_tot, dint, dslp, thres_step_max, step_max, df0, T0, Td0, Tdd0, Te, Dose, nonDose, TTerm, nonDose_LIN, nonDose_PLIN, nonDose_LOGLIN, modelform, gmix_theta, gmix_term, convgd, lr, optim_para, maxiter, Lin_Sys, Lin_Res, term_n, dfc, halfmax, epsilon, deriv_epsilon);
        temp_L[0] = reg_out["LogLik"];
        if (!is_nan(temp_L)[0]) {
            loop_check = false;
        }
    }
    if (loop_check) {
        limit_hit[0] = true;
        limits[0] = 0;
        ll_final[0] = 0;
        limit_converged[0] = false;
    } else {
        //  Now we can run the actual algorithm
        reg_beta = reg_out["beta_0"];
        for (int ij = 0; ij < totalnum; ij++) {
            beta_L[ij] = reg_beta[ij];
        }
        L_L = reg_out["LogLik"];
        //  now set the mid point value
        for (int ij = 0; ij < totalnum; ij++) {
            beta_M[ij] = (beta_H[ij] + beta_L[ij])/2;
            beta_0[ij] = beta_M[ij];
        }
        reg_out = Cox_Full_Run(reqrdnum, ntime, tform, RiskFail,  RiskPairs, RiskPairs_Strata, totalnum, fir, R, Rd, Rdd, Rls1, Rls2, Rls3, Lls1, Lls2, Lls3, cens_weight, Strata_vals, beta_0, RdR, RddR, Ll, Lld, Lldd, nthreads, KeepConstant, ties_method, verbose, model_bool, iter_stop, term_tot, dint, dslp, thres_step_max, step_max, df0, T0, Td0, Tdd0, Te, Dose, nonDose, TTerm, nonDose_LIN, nonDose_PLIN, nonDose_LOGLIN, modelform, gmix_theta, gmix_term, convgd, lr, optim_para, maxiter, Lin_Sys, Lin_Res, term_n, dfc, halfmax, epsilon, deriv_epsilon);
        reg_beta = reg_out["beta_0"];
        for (int ij = 0; ij < totalnum; ij++) {
            beta_M[ij] = reg_beta[ij];
        }
        L_M = reg_out["LogLik"];
        //
        int step = 0;
        //  now run the bisection until stopping point
        while ((step < maxstep) && (abs(beta_L[para_number] - beta_H[para_number]) > epsilon) && (!limit_hit[0])) {
            step = step + 1;
            if (L_H < Lstar) {
                throw invalid_argument("The upper estimate is too high?");
            } else if (L_L < Lstar) {
                //  midpoint is in between the two
                if (L_M > Lstar) {
                    //  the mid point is past the optimum
                    for (int ij = 0; ij < totalnum; ij++) {
                        beta_H[ij] = beta_M[ij];
                    }
                    L_H = L_M;
                } else {
                    //  the mid point is before the optimum
                    for (int ij = 0; ij < totalnum; ij++) {
                        beta_L[ij] = beta_M[ij];
                    }
                    L_L = L_M;
                }
            } else if (L_M < Lstar) {
                for (int ij = 0; ij < totalnum; ij++) {
                    beta_L[ij] = beta_M[ij];
                }
                L_L = L_M;
            } else {
                //  the lower estimate needs to be shifted down
                for (int ij = 0; ij < totalnum; ij++) {
                    beta_H[ij] = beta_L[ij];
                }
                L_H = L_L;
                //  check new high point
                loop_check = true;
                temp_step = step_size;
                while ((loop_check) && (temp_step > 1e-3)) {
                    //  assign new high point
                    beta_L[para_number] = beta_H[para_number] - temp_step;
                    for (int ij = 0; ij < totalnum; ij++) {
                        beta_0[ij] = beta_L[ij];
                    }
                    temp_step = temp_step * 0.5;
                    //
                    reg_out = Cox_Full_Run(reqrdnum, ntime, tform, RiskFail,  RiskPairs, RiskPairs_Strata, totalnum, fir, R, Rd, Rdd, Rls1, Rls2, Rls3, Lls1, Lls2, Lls3, cens_weight, Strata_vals, beta_0, RdR, RddR, Ll, Lld, Lldd, nthreads, KeepConstant, ties_method, verbose, model_bool, iter_stop, term_tot, dint, dslp, thres_step_max, step_max, df0, T0, Td0, Tdd0, Te, Dose, nonDose, TTerm, nonDose_LIN, nonDose_PLIN, nonDose_LOGLIN, modelform, gmix_theta, gmix_term, convgd, lr, optim_para, maxiter, Lin_Sys, Lin_Res, term_n, dfc, halfmax, epsilon, deriv_epsilon);
                    temp_L[0] = reg_out["LogLik"];
                    if (!is_nan(temp_L)[0]) {
                        loop_check = false;
                    }
                }
                if (loop_check) {
                    limit_hit[0] = true;
                    limits[0] = 0;
                    ll_final[0] = 0;
                    limit_converged[0] = false;
                } else {
                    reg_beta = reg_out["beta_0"];
                    for (int ij = 0; ij < totalnum; ij++) {
                        beta_L[ij] = reg_beta[ij];
                    }
                    L_L = reg_out["LogLik"];
                }
            }
            if (!limit_hit[0]) {
                //  now set the mid point value
                for (int ij = 0; ij < totalnum; ij++) {
                    beta_M[ij] = (beta_H[ij] + beta_L[ij])/2;
                    beta_0[ij] = beta_M[ij];
                }
                reg_out = Cox_Full_Run(reqrdnum, ntime, tform, RiskFail,  RiskPairs, RiskPairs_Strata, totalnum, fir, R, Rd, Rdd, Rls1, Rls2, Rls3, Lls1, Lls2, Lls3, cens_weight, Strata_vals, beta_0, RdR, RddR, Ll, Lld, Lldd, nthreads, KeepConstant, ties_method, verbose, model_bool, iter_stop, term_tot, dint, dslp, thres_step_max, step_max, df0, T0, Td0, Tdd0, Te, Dose, nonDose, TTerm, nonDose_LIN, nonDose_PLIN, nonDose_LOGLIN, modelform, gmix_theta, gmix_term, convgd, lr, optim_para, maxiter, Lin_Sys, Lin_Res, term_n, dfc, halfmax, epsilon, deriv_epsilon);
                reg_beta = reg_out["beta_0"];
                for (int ij = 0; ij < totalnum; ij++) {
                    beta_M[ij] = reg_beta[ij];
                }
                L_M = reg_out["LogLik"];
            }
        }
        if ((abs(beta_L[para_number] - beta_H[para_number]) < epsilon) && (!limit_hit[0])) {
            limit_converged[0] = true;
        }
        limits[0] = beta_M[para_number];
        ll_final[0] = L_M;
        width_final[0] = abs(beta_L[para_number] - beta_H[para_number]);
        step_final[0] = step;
    }
    res_list = List::create(_["Parameter_Limits"] = wrap(limits), _["Negative_Risk_Limit_Hit"] = wrap(limit_hit), _["Likelihood_Boundary"] = wrap(ll_final), _["Likelihood_Goal"] = wrap(Lstar), _["Limit_Converged"] = wrap(limit_converged), _["Final_Window_Width"] = wrap(width_final), _["Final_Step"] = wrap(step_final), _["Status"] = "PASSED");
    //  returns a list of results
    return res_list;
}

//' Primary Cox PH likelihood bounds calcualtion function.
//'
//' \code{LogLik_Cox_PH_Omnibus_Log_Bound} Performs the calls to calculation functions and log-likeihood profile bounds
//'
//' @inheritParams CPP_template
//'
//' @return List of final results: Log-likelihood of optimum, first derivative of log-likelihood, second derivative matrix, parameter list, standard deviation estimate, AIC, model information
//' @noRd
//'
//
List LogLik_Poisson_Omnibus_Log_Bound_CurveSearch(const Ref<const MatrixXd>& PyrC, const Ref<const MatrixXd>& dfs, IntegerVector term_n, StringVector tform, Ref<VectorXd> beta_0, Ref<MatrixXd> df0, IntegerVector dfc, int fir, string modelform, double lr, List optim_para, int maxiter, int halfmax, double epsilon, double step_max, double thres_step_max, double deriv_epsilon, int verbose, IntegerVector KeepConstant, int term_tot, int nthreads, List model_bool, const double gmix_theta, const IntegerVector gmix_term, const Ref<const MatrixXd>& Lin_Sys, const Ref<const VectorXd>& Lin_Res, double qchi, int para_number, int maxstep, double step_size) {
    //
    List temp_list = List::create(_["Status"] = "TEMP");  //  used as a dummy return value for code checking
    //  Time durations are measured from this point on in microseconds
    //
    //  df0: covariate data
    //  ntime: number of event times for Cox PH
    //  totalnum: number of terms used
    //
    //  ------------------------------------------------------------------------- //  initialize
    const int mat_row = df0.rows();
    //  int ntime = tu.size();
    int totalnum = term_n.size();
    int reqrdnum = totalnum - sum(KeepConstant);
    //  ------------------------------------------------------------------------- //  initialize
    if (model_bool["null"]) {
        if (verbose >= 1) {
            Rcout << "null model is not compatable with log-based bound calculation" << endl;
        }
        temp_list = List::create(_["Status"] = "FAILED_BAD_MODEL_NULL", _["LogLik"] = R_NaN);
        return temp_list;
    }
    if (model_bool["single"]) {
        if (verbose >= 1) {
            Rcout << "non-derivative model calculation is not compatable with log-based bound calculation" << endl;
        }
        temp_list = List::create(_["Status"] = "FAILED_WITH_BAD_MODEL_SINGLE", _["LogLik"] = R_NaN);
        return temp_list;
    }
    if (model_bool["gradient"]) {
        if (verbose >= 1) {
            Rcout << "gradient descent model calculation is not compatable with log-based bound calculation" << endl;
        }
        temp_list = List::create(_["Status"] = "FAILED_WITH_BAD_MODEL_GRADIENT", _["LogLik"] = R_NaN);
        return temp_list;
    }
    if (model_bool["constraint"]) {
        if (verbose >= 1) {
            Rcout << "linear constataints are not compatable with Case-Control model calculation" << endl;
        }
        temp_list = List::create(_["Status"] = "FAILED_WITH_BAD_MODEL_CONSTRAINT", _["LogLik"] = R_NaN);
        return temp_list;
    }
    //
    //  cout.precision: controls the number of significant digits printed
    //  nthreads: number of threads used for parallel operations
    //
    Rcout.precision(7);  //  forces higher precision numbers printed to terminal
    //
    //  Lld_worst: stores the highest magnitude log-likelihood derivative
    //  ---------------------------------------------
    //  To Start, needs to seperate the derivative terms
    //  ---------------------------------------------
    //
    //  ------------------------------------------------------------------------- //  initialize
    // Map<VectorXd> beta_0(as<Map<VectorXd> >(a_n));
    VectorXd beta_max = beta_0;
    MatrixXd T0;
    MatrixXd Td0;
    MatrixXd Tdd0;
    //
    MatrixXd Te;
    MatrixXd R;
    MatrixXd Rd;
    MatrixXd Rdd;
    //
    MatrixXd Dose;
    MatrixXd nonDose;
    MatrixXd nonDose_LIN;
    MatrixXd nonDose_PLIN;
    MatrixXd nonDose_LOGLIN;
    MatrixXd TTerm;
    double dint = 0.0;  //  the amount of change used to calculate derivatives in threshold paramters
    double dslp = 0.0;
    MatrixXd RdR;
    MatrixXd RddR;
    //  ------------------------------------------------------------------------- //  initialize
    //  ---------------------------------------------
    //  To Start, needs to seperate the derivative terms
    //  ---------------------------------------------
    //
    T0 = MatrixXd::Zero(mat_row, totalnum);  //  preallocates matrix for Term column
    Te = MatrixXd::Zero(mat_row, 1);  //  preallocates matrix for column terms used for temporary storage
    R = MatrixXd::Zero(mat_row, 1);  //  preallocates matrix for Risks
    //
    Dose = MatrixXd::Constant(mat_row, term_tot, 0.0);  //  matrix of the total dose term values
    nonDose = MatrixXd::Constant(mat_row, term_tot, 1.0);  //  matrix of the total non-dose term values
    nonDose_LIN = MatrixXd::Constant(mat_row, term_tot, 0.0);  //  matrix of Linear subterm values
    nonDose_PLIN = MatrixXd::Constant(mat_row, term_tot, 1.0);  //  matrix of Loglinear subterm values
    nonDose_LOGLIN = MatrixXd::Constant(mat_row, term_tot, 1.0);  //  matrix of Product linear subterm values
    TTerm = MatrixXd::Zero(mat_row, term_tot);  //  matrix of term values
    if (model_bool["single"]) {
    } else {
        Td0 = MatrixXd::Zero(mat_row, reqrdnum);  //  preallocates matrix for Term derivative columns
        Tdd0 = MatrixXd::Zero(mat_row, reqrdnum*(reqrdnum + 1)/2);  //  preallocates matrix for Term second derivative columns
        Rd = MatrixXd::Zero(mat_row, reqrdnum);  //  preallocates matrix for Risk derivatives
        Rdd = MatrixXd::Zero(mat_row, reqrdnum*(reqrdnum + 1)/2);  //  preallocates matrix for Risk second derivatives
        dint = thres_step_max;  //  the amount of change used to calculate derivatives in threshold paramters
        dslp = step_max;
        RdR = MatrixXd::Zero(mat_row, reqrdnum);  //  preallocates matrix for Risk to derivative ratios
        RddR = MatrixXd::Zero(mat_row, reqrdnum*(reqrdnum + 1)/2);  //  preallocates matrix for Risk to second derivative ratios
    }
    VectorXd s_weights;
    if (model_bool["strata"]) {
        s_weights = VectorXd::Zero(mat_row);
        Gen_Strat_Weight(modelform, dfs, PyrC, s_weights, nthreads, tform, term_n, term_tot);
    }
    //  ------------------------------------------------------------------------- //  initialize
    MatrixXd Rls1;
    MatrixXd Lls1;
    MatrixXd Rls2;
    MatrixXd Rls3;
    MatrixXd Lls2;
    MatrixXd Lls3;
    vector<double> Ll(reqrdnum, 0.0);  //  log-likelihood values
    vector<double> Lld(reqrdnum, 0.0);  //  log-likelihood derivative values
    vector<double> Lldd(pow(reqrdnum, 2), 0.0);  //  the second derivative matrix has room for every combination, but only the lower triangle is calculated initially
    MatrixXd dev_temp = MatrixXd::Zero(PyrC.rows(), 2);
    double dev = 0;
    //  ------------------------------------------------------------------------- //  initialize
    //  the log-likelihood is calculated in parallel over the risk groups
    vector <double> Ll_comp(2, Ll[0]);  //  vector to compare values
    double step_max0 = step_max;
    double thres_step_max0 = thres_step_max;
    //
    vector<double> dbeta(totalnum, 0.0);
    vector<double> dbeta_start(totalnum, 0.0);
    //
    //  --------------------------
    //  always starts from initial guess
    //  --------------------------
    vector<double> beta_c(totalnum, 0.0);
    vector<double> beta_a(totalnum, 0.0);
    vector<double> beta_best(totalnum, 0.0);
    vector<double> beta_p(totalnum, 0.0);
    vector<double> beta_peak(totalnum, 0.0);
    VectorXd::Map(&beta_p[0], beta_0.size()) = beta_0;  //  stores previous parameters
    VectorXd::Map(&beta_c[0], beta_0.size()) = beta_0;  //  stores current parameters
    VectorXd::Map(&beta_a[0], beta_0.size()) = beta_0;  //  stores a refrence value for parameters
    VectorXd::Map(&beta_best[0], beta_0.size()) = beta_0;  //  stores the best parameters
    VectorXd::Map(&beta_peak[0], beta_0.size()) = beta_0;  //  stores the peak parameters
    int iter_stop  = 0;  //  tracks if the iterations should be stopped for convergence
    vector<double> beta_abs_best(totalnum, 0.0);
    fill(Ll.begin(), Ll.end(), 0.0);
    fill(Lld.begin(), Lld.end(), 0.0);
    fill(Lldd.begin(), Lldd.end(), 0.0);
    beta_p = beta_peak;  //
    beta_a = beta_peak;  //
    beta_c = beta_peak;  //
    step_max = step_max0;
    thres_step_max = thres_step_max0;
    for (int i = 0; i < beta_0.size(); i++) {
        beta_0[i] = beta_peak[i];
    }
    Pois_Term_Risk_Calc(modelform, tform, term_n, totalnum, fir, dfc, term_tot, T0, Td0, Tdd0, Te, R, Rd, Rdd, Dose, nonDose, beta_0, df0, dint, dslp, TTerm, nonDose_LIN, nonDose_PLIN, nonDose_LOGLIN, RdR, RddR, s_weights, nthreads, KeepConstant, verbose, model_bool, gmix_theta, gmix_term);
    Pois_Dev_LL_Calc(reqrdnum, totalnum, fir, R, Rd, Rdd, beta_0, RdR, RddR, Ll, Lld, Lldd, PyrC, dev_temp, nthreads, KeepConstant, verbose, model_bool, iter_stop, dev);
    //
    double Lstar = Ll[0]-qchi;
    double Lpeak = Ll[0];
    IntegerVector KeepConstant_trouble(totalnum);
    for (int i = 0;  i < totalnum; i++) {
        KeepConstant_trouble[i] = KeepConstant[i];
    }
    KeepConstant_trouble[para_number] = 1;
    //
    vector<double> limits(2, 0.0);
    vector<bool>   limit_hit(2, FALSE);
    vector<bool>   limit_converged(2, FALSE);
    vector<double> ll_final(2, 0.0);
    vector<double> width_final(2, 0.0);
    vector<int>    step_final(2, 0.0);
    List res_list;
    //
    if (verbose >= 4) {
        Rcout << "C++ Note: STARTING Upper Bound" << endl;
    }
    bool convgd = false;
    ///
    //  We need the values reserved for the upper, middle, lower estimates and scores
    vector<double> beta_L(totalnum, 0.0);
    vector<double> beta_M(totalnum, 0.0);
    vector<double> beta_H(totalnum, 0.0);
    double L_L = 0.0;
    double L_M = 0.0;
    double L_H = 0.0;
    NumericVector reg_beta(totalnum);
    //  First we need to establish the first interval estimates
    for (int ij = 0; ij < totalnum; ij++) {
        beta_L[ij] = beta_peak[ij];
        beta_H[ij] = beta_peak[ij];
    }
    L_L = Lpeak;
    bool loop_check = true;
    double temp_step = step_size;
    NumericVector temp_L(1);
    List reg_out;
    while ((loop_check) && (temp_step > 1e-3)) {
        //  assign new high point
        beta_H[para_number] = beta_L[para_number] + temp_step;
        for (int ij = 0; ij < totalnum; ij++) {
            beta_0[ij] = beta_H[ij];
        }
        temp_step = temp_step * 0.5;
        //
        reg_out = Pois_Full_Run(PyrC, reqrdnum, tform, totalnum, fir, R, Rd, Rdd, s_weights, beta_0, RdR, RddR, Ll, Lld, Lldd, nthreads, KeepConstant, verbose, model_bool, iter_stop, term_tot, dint, dslp, thres_step_max, step_max, df0, T0, Td0, Tdd0, Te, Dose, nonDose, TTerm, nonDose_LIN, nonDose_PLIN, nonDose_LOGLIN, modelform, gmix_theta, gmix_term, convgd, lr, optim_para, maxiter, Lin_Sys, Lin_Res, term_n, dfc, halfmax, epsilon, deriv_epsilon);
        temp_L[0] = reg_out["LogLik"];
        if (!is_nan(temp_L)[0]) {
            loop_check = false;
        }
    }
    if (loop_check) {
        limit_hit[1] = true;
        limits[1] = 0;
        ll_final[1] = 0;
        limit_converged[1] = false;
    } else {
        //  Now we can run the actual algorithm
        reg_beta = reg_out["beta_0"];
        for (int ij = 0; ij < totalnum; ij++) {
            beta_H[ij] = reg_beta[ij];
        }
        L_H = reg_out["LogLik"];
        //  now set the mid point value
        for (int ij = 0; ij < totalnum; ij++) {
            beta_M[ij] = (beta_H[ij] + beta_L[ij])/2;
            beta_0[ij] = beta_M[ij];
        }
        reg_out = Pois_Full_Run(PyrC, reqrdnum, tform, totalnum, fir, R, Rd, Rdd, s_weights, beta_0, RdR, RddR, Ll, Lld, Lldd, nthreads, KeepConstant, verbose, model_bool, iter_stop, term_tot, dint, dslp, thres_step_max, step_max, df0, T0, Td0, Tdd0, Te, Dose, nonDose, TTerm, nonDose_LIN, nonDose_PLIN, nonDose_LOGLIN, modelform, gmix_theta, gmix_term, convgd, lr, optim_para, maxiter, Lin_Sys, Lin_Res, term_n, dfc, halfmax, epsilon, deriv_epsilon);
        reg_beta = reg_out["beta_0"];
        for (int ij = 0; ij < totalnum; ij++) {
            beta_M[ij] = reg_beta[ij];
        }
        L_M = reg_out["LogLik"];
        //
        int step = 0;
        //  now run the bisection until stopping point
        while ((step < maxstep) && (abs(beta_L[para_number] - beta_H[para_number]) > epsilon) && (!limit_hit[1])) {
            step = step + 1;
            if (L_L < Lstar) {
                throw invalid_argument("The lower estimate is too high?");
            } else if (L_H < Lstar) {
                //  midpoint is in between the two
                if (L_M < Lstar) {
                    //  the mid point is past the optimum
                    for (int ij = 0; ij < totalnum; ij++) {
                        beta_H[ij] = beta_M[ij];
                    }
                    L_H = L_M;
                } else {
                    //  the mid point is before the optimum
                    for (int ij = 0; ij < totalnum; ij++) {
                        beta_L[ij] = beta_M[ij];
                    }
                    L_L = L_M;
                }
            } else if (L_M < Lstar) {
                for (int ij = 0; ij < totalnum; ij++) {
                    beta_H[ij] = beta_M[ij];
                }
                L_H = L_M;
            } else {
                //  the upper estimate needs to be shifted up
                for (int ij = 0; ij < totalnum; ij++) {
                    beta_L[ij] = beta_H[ij];
                }
                L_L = L_H;
                //  check new high point
                loop_check = true;
                temp_step = step_size;
                while ((loop_check) && (temp_step > 1e-3)) {
                    //  assign new high point
                    beta_H[para_number] = beta_L[para_number] + temp_step;
                    for (int ij = 0; ij < totalnum; ij++) {
                        beta_0[ij] = beta_H[ij];
                    }
                    temp_step = temp_step * 0.5;
                    //
                    reg_out = Pois_Full_Run(PyrC, reqrdnum, tform, totalnum, fir, R, Rd, Rdd, s_weights, beta_0, RdR, RddR, Ll, Lld, Lldd, nthreads, KeepConstant, verbose, model_bool, iter_stop, term_tot, dint, dslp, thres_step_max, step_max, df0, T0, Td0, Tdd0, Te, Dose, nonDose, TTerm, nonDose_LIN, nonDose_PLIN, nonDose_LOGLIN, modelform, gmix_theta, gmix_term, convgd, lr, optim_para, maxiter, Lin_Sys, Lin_Res, term_n, dfc, halfmax, epsilon, deriv_epsilon);
                    temp_L[0] = reg_out["LogLik"];
                    if (!is_nan(temp_L)[0]) {
                        loop_check = false;
                    }
                }
                if (loop_check) {
                    limit_hit[1] = true;
                    limits[1] = 0;
                    ll_final[1] = 0;
                    limit_converged[1] = false;
                } else {
                    reg_beta = reg_out["beta_0"];
                    for (int ij = 0; ij < totalnum; ij++) {
                        beta_H[ij] = reg_beta[ij];
                    }
                    L_H = reg_out["LogLik"];
                }
            }
            if (!limit_hit[1]) {
                //  now set the mid point value
                for (int ij = 0; ij < totalnum; ij++) {
                    beta_M[ij] = (beta_H[ij] + beta_L[ij])/2;
                    beta_0[ij] = beta_M[ij];
                }
                reg_out = Pois_Full_Run(PyrC, reqrdnum, tform, totalnum, fir, R, Rd, Rdd, s_weights, beta_0, RdR, RddR, Ll, Lld, Lldd, nthreads, KeepConstant, verbose, model_bool, iter_stop, term_tot, dint, dslp, thres_step_max, step_max, df0, T0, Td0, Tdd0, Te, Dose, nonDose, TTerm, nonDose_LIN, nonDose_PLIN, nonDose_LOGLIN, modelform, gmix_theta, gmix_term, convgd, lr, optim_para, maxiter, Lin_Sys, Lin_Res, term_n, dfc, halfmax, epsilon, deriv_epsilon);
                reg_beta = reg_out["beta_0"];
                for (int ij = 0; ij < totalnum; ij++) {
                    beta_M[ij] = reg_beta[ij];
                }
                L_M = reg_out["LogLik"];
            }
        }
        if ((abs(beta_L[para_number] - beta_H[para_number]) < epsilon) && (!limit_hit[1])) {
            limit_converged[1] = true;
        }
        limits[1] = beta_M[para_number];
        ll_final[1] = L_M;
        width_final[1] = abs(beta_L[para_number] - beta_H[para_number]);
        step_final[1] = step;
    }
    //  upper limit found, now solve lower limit
    for (int ij = 0; ij < totalnum; ij++) {
        beta_L[ij] = beta_peak[ij];
        beta_H[ij] = beta_peak[ij];
    }
    L_H = Lpeak;
    loop_check = true;
    temp_step = step_size;
    while ((loop_check) && (temp_step > 1e-3)) {
        //  assign new high point
        beta_L[para_number] = beta_H[para_number] - temp_step;
        for (int ij = 0; ij < totalnum; ij++) {
            beta_0[ij] = beta_L[ij];
        }
        temp_step = temp_step * 0.5;
        //
        reg_out = Pois_Full_Run(PyrC, reqrdnum, tform, totalnum, fir, R, Rd, Rdd, s_weights, beta_0, RdR, RddR, Ll, Lld, Lldd, nthreads, KeepConstant, verbose, model_bool, iter_stop, term_tot, dint, dslp, thres_step_max, step_max, df0, T0, Td0, Tdd0, Te, Dose, nonDose, TTerm, nonDose_LIN, nonDose_PLIN, nonDose_LOGLIN, modelform, gmix_theta, gmix_term, convgd, lr, optim_para, maxiter, Lin_Sys, Lin_Res, term_n, dfc, halfmax, epsilon, deriv_epsilon);
        temp_L[0] = reg_out["LogLik"];
        if (!is_nan(temp_L)[0]) {
            loop_check = false;
        }
    }
    if (loop_check) {
        limit_hit[0] = true;
        limits[0] = 0;
        ll_final[0] = 0;
        limit_converged[0] = false;
    } else {
        //  Now we can run the actual algorithm
        reg_beta = reg_out["beta_0"];
        for (int ij = 0; ij < totalnum; ij++) {
            beta_L[ij] = reg_beta[ij];
        }
        L_L = reg_out["LogLik"];
        //  now set the mid point value
        for (int ij = 0; ij < totalnum; ij++) {
            beta_M[ij] = (beta_H[ij] + beta_L[ij])/2;
            beta_0[ij] = beta_M[ij];
        }
        reg_out = Pois_Full_Run(PyrC, reqrdnum, tform, totalnum, fir, R, Rd, Rdd, s_weights, beta_0, RdR, RddR, Ll, Lld, Lldd, nthreads, KeepConstant, verbose, model_bool, iter_stop, term_tot, dint, dslp, thres_step_max, step_max, df0, T0, Td0, Tdd0, Te, Dose, nonDose, TTerm, nonDose_LIN, nonDose_PLIN, nonDose_LOGLIN, modelform, gmix_theta, gmix_term, convgd, lr, optim_para, maxiter, Lin_Sys, Lin_Res, term_n, dfc, halfmax, epsilon, deriv_epsilon);
        reg_beta = reg_out["beta_0"];
        for (int ij = 0; ij < totalnum; ij++) {
            beta_M[ij] = reg_beta[ij];
        }
        L_M = reg_out["LogLik"];
        //
        int step = 0;
        //  now run the bisection until stopping point
        while ((step < maxstep) && (abs(beta_L[para_number] - beta_H[para_number]) > epsilon) && (!limit_hit[0])) {
            step = step + 1;
            if (L_H < Lstar) {
                throw invalid_argument("The upper estimate is too high?");
            } else if (L_L < Lstar) {
                //  midpoint is in between the two
                if (L_M > Lstar) {
                    //  the mid point is past the optimum
                    for (int ij = 0; ij < totalnum; ij++) {
                        beta_H[ij] = beta_M[ij];
                    }
                    L_H = L_M;
                } else {
                    //  the mid point is before the optimum
                    for (int ij = 0; ij < totalnum; ij++) {
                        beta_L[ij] = beta_M[ij];
                    }
                    L_L = L_M;
                }
            } else if (L_M < Lstar) {
                for (int ij = 0; ij < totalnum; ij++) {
                    beta_L[ij] = beta_M[ij];
                }
                L_L = L_M;
            } else {
                //  the lower estimate needs to be shifted down
                for (int ij = 0; ij < totalnum; ij++) {
                    beta_H[ij] = beta_L[ij];
                }
                L_H = L_L;
                //  check new high point
                loop_check = true;
                temp_step = step_size;
                while ((loop_check) && (temp_step > 1e-3)) {
                    //  assign new high point
                    beta_L[para_number] = beta_H[para_number] - temp_step;
                    for (int ij = 0; ij < totalnum; ij++) {
                        beta_0[ij] = beta_L[ij];
                    }
                    temp_step = temp_step * 0.5;
                    //
                    reg_out = Pois_Full_Run(PyrC, reqrdnum, tform, totalnum, fir, R, Rd, Rdd, s_weights, beta_0, RdR, RddR, Ll, Lld, Lldd, nthreads, KeepConstant, verbose, model_bool, iter_stop, term_tot, dint, dslp, thres_step_max, step_max, df0, T0, Td0, Tdd0, Te, Dose, nonDose, TTerm, nonDose_LIN, nonDose_PLIN, nonDose_LOGLIN, modelform, gmix_theta, gmix_term, convgd, lr, optim_para, maxiter, Lin_Sys, Lin_Res, term_n, dfc, halfmax, epsilon, deriv_epsilon);
                    temp_L[0] = reg_out["LogLik"];
                    if (!is_nan(temp_L)[0]) {
                        loop_check = false;
                    }
                }
                if (loop_check) {
                    limit_hit[0] = true;
                    limits[0] = 0;
                    ll_final[0] = 0;
                    limit_converged[0] = false;
                } else {
                    reg_beta = reg_out["beta_0"];
                    for (int ij = 0; ij < totalnum; ij++) {
                        beta_L[ij] = reg_beta[ij];
                    }
                    L_L = reg_out["LogLik"];
                }
            }
            if (!limit_hit[0]) {
                //  now set the mid point value
                for (int ij = 0; ij < totalnum; ij++) {
                    beta_M[ij] = (beta_H[ij] + beta_L[ij])/2;
                    beta_0[ij] = beta_M[ij];
                }
                reg_out = Pois_Full_Run(PyrC, reqrdnum, tform, totalnum, fir, R, Rd, Rdd, s_weights, beta_0, RdR, RddR, Ll, Lld, Lldd, nthreads, KeepConstant, verbose, model_bool, iter_stop, term_tot, dint, dslp, thres_step_max, step_max, df0, T0, Td0, Tdd0, Te, Dose, nonDose, TTerm, nonDose_LIN, nonDose_PLIN, nonDose_LOGLIN, modelform, gmix_theta, gmix_term, convgd, lr, optim_para, maxiter, Lin_Sys, Lin_Res, term_n, dfc, halfmax, epsilon, deriv_epsilon);
                reg_beta = reg_out["beta_0"];
                for (int ij = 0; ij < totalnum; ij++) {
                    beta_M[ij] = reg_beta[ij];
                }
                L_M = reg_out["LogLik"];
            }
        }
        if ((abs(beta_L[para_number] - beta_H[para_number]) < epsilon) && (!limit_hit[0])) {
            limit_converged[0] = true;
        }
        limits[0] = beta_M[para_number];
        ll_final[0] = L_M;
        width_final[0] = abs(beta_L[para_number] - beta_H[para_number]);
        step_final[0] = step;
    }
    res_list = List::create(_["Parameter_Limits"] = wrap(limits), _["Negative_Risk_Limit_Hit"] = wrap(limit_hit), _["Likelihood_Boundary"] = wrap(ll_final), _["Likelihood_Goal"] = wrap(Lstar), _["Limit_Converged"] = wrap(limit_converged), _["Final_Window_Width"] = wrap(width_final), _["Final_Step"] = wrap(step_final), _["Status"] = "PASSED");
    //  returns a list of results
    return res_list;
}
