//  Copyright 2022 - 2025, Eric Giunta and the project collaborators, Please see main R package for license and usage details

#include <RcppEigen.h>

#include "Main_Functions.h"
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
using std::reduce;
using std::accumulate;
using std::isinf;
using std::isnan;

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

//' checks if the model is viable
//'
//' \code{Check_Risk} Calculates risks and checks for negative values
//'
//' @inheritParams CPP_template
//'
//' @return True for viable point, False for negative error
//' @noRd
//'
//
bool Check_Risk(IntegerVector term_n, StringVector tform, Ref<VectorXd> beta_0, Ref<MatrixXd> df0, IntegerVector dfc, int fir, string modelform, int verbose, IntegerVector KeepConstant, int term_tot, int nthreads, const double gmix_theta, const IntegerVector gmix_term) {
    //
    List temp_list = List::create(_["Status"] = "FAILED");  //  used as a dummy return value for code checking
    if (verbose >= 3) {
        Rcout << "C++ Note: START_RISK_CHECK" << endl;
    }
    //
    //
    int totalnum = term_n.size();
    //
    Rcout.precision(7);  //  forces higher precision numbers printed to terminal
    //  ---------------------------------------------
    //  To Start, needs to seperate the derivative terms
    //  ---------------------------------------------
    //
    // Map<VectorXd> beta_0(as<Map<VectorXd> >(a_n));
    const int mat_row = df0.rows();
    MatrixXd T0 = MatrixXd::Zero(mat_row, totalnum);  //  preallocates matrix for Term column
    //
    MatrixXd Te = MatrixXd::Zero(mat_row, 1);  //  preallocates matrix for column terms used for temporary storage
    MatrixXd R = MatrixXd::Zero(mat_row, 1);  //  preallocates matrix for Risks
    //
    MatrixXd Dose = MatrixXd::Constant(mat_row, term_tot, 0.0);  //  matrix of the total dose term values
    MatrixXd nonDose = MatrixXd::Constant(mat_row, term_tot, 1.0);  //  matrix of the total non-dose term values
    MatrixXd nonDose_LIN = MatrixXd::Constant(mat_row, term_tot, 0.0);  //  matrix of Linear subterm values
    MatrixXd nonDose_PLIN = MatrixXd::Constant(mat_row, term_tot, 1.0);  //  matrix of Loglinear subterm values
    MatrixXd nonDose_LOGLIN = MatrixXd::Constant(mat_row, term_tot, 1.0);  //  matrix of Product linear subterm values
    MatrixXd TTerm = MatrixXd::Zero(mat_row, term_tot);  //  matrix of term values
    //
    //  Calculates the subterm and term values
    Make_subterms_Single(totalnum, term_n, tform, dfc, fir, T0, Dose, nonDose, TTerm, nonDose_LIN, nonDose_PLIN, nonDose_LOGLIN, beta_0, df0, nthreads, KeepConstant);
    //  ---------------------------------------------------------
    //  Prints off a series of calculations to check at what point values are changing
    //  ---------------------------------------------------------
    //  Calculates the risk for each row
    Make_Risks_Single(modelform, tform, term_n, totalnum, fir, T0, Te, R, Dose, nonDose, TTerm, nonDose_LIN, nonDose_PLIN, nonDose_LOGLIN, nthreads, KeepConstant, gmix_theta, gmix_term);
    //
    //  Removes infinite values
    //
    if ((R.minCoeff() <= 0) || (R.hasNaN())) {
        return FALSE;
    }
    return TRUE;
}

//' Primary Cox PH regression with multiple starting points and optional combinations of null, stratification, competing risks, multiplicative log-linear model, and no derivative calculation.
//'
//' \code{LogLik_Cox_PH_Omnibus} Performs the calls to calculation functions, Structures the Cox PH regression, With verbose option prints out time stamps and intermediate sums of terms and derivatives
//'
//' @inheritParams CPP_template
//'
//' @return List of final results: Log-likelihood of optimum, first derivative of log-likelihood, second derivative matrix, parameter list, standard deviation estimate, AIC, model information
//' @noRd
//'
//
List LogLik_Cox_PH_Omnibus(IntegerVector term_n, StringVector tform, NumericMatrix& a_ns, Ref<MatrixXd> df0, IntegerVector dfc, int fir, string modelform, double lr, List optim_para, NumericVector maxiters, int guesses, int halfmax, double epsilon, double step_max, double thres_step_max, double deriv_epsilon, const Ref<const MatrixXd>& df_m, NumericVector tu, int verbose, IntegerVector KeepConstant, int term_tot, string ties_method, int nthreads, NumericVector& Strata_vals, const VectorXd& cens_weight, List model_bool, const double gmix_theta, const IntegerVector gmix_term, const Ref<const MatrixXd>& Lin_Sys, const Ref<const VectorXd>& Lin_Res) {
    //
    List temp_list = List::create(_["Status"] = "TEMP");  //  used as a dummy return value for code checking
    //
    //  df0: covariate data
    //  ntime: number of event times for Cox PH
    //  totalnum: number of terms used
    //
    //  ------------------------------------------------------------------------- //  initialize
    const int mat_row = df0.rows();
    int ntime = tu.size();
    int totalnum;
    int reqrdnum;
    //  ------------------------------------------------------------------------- //  initialize
    if (!model_bool["null"]) {
        totalnum = term_n.size();
        reqrdnum = totalnum - sum(KeepConstant);
    } else {
        totalnum = 1;
        reqrdnum = 1;
    }
    //
    //  cout.precision: controls the number of significant digits printed
    //  nthreads: number of threads used for parallel operations
    //
    Rcout.precision(7);  //  forces higher precision numbers printed to terminal
    //
    //  Lld_worst: stores the highest magnitude log-likelihood derivative
    //
    //
    double Lld_worst = 0.0;  //  stores derivative value used to determine if every parameter is near convergence
    double dbeta_max = 0.0;  //  stores the largest step taken, determines if the step sizes are small enough for convergence
    //
    //  ---------------------------------------------
    //  To Start, needs to seperate the derivative terms
    //  ---------------------------------------------
    //
    //  ------------------------------------------------------------------------- //  initialize
    NumericVector a_n = a_ns.row(0);
    Map<VectorXd> beta_0(as<Map<VectorXd> >(a_n));
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
    if (!model_bool["null"]) {
        //  ---------------------------------------------
        //  To Start, needs to seperate the derivative terms
        //  ---------------------------------------------
        //
        Cox_Refresh_R_TERM(totalnum, reqrdnum, term_tot, dint, dslp, thres_step_max, step_max, df0, T0, Td0, Tdd0, Te, R, Rd, Rdd, Dose, nonDose, TTerm, nonDose_LIN, nonDose_PLIN, nonDose_LOGLIN, RdR, RddR, model_bool);
    } else {
        R = MatrixXd::Constant(mat_row, 1, 1.0);
    }
    //  ------------------------------------------------------------------------- //  initialize
    IntegerMatrix RiskFail;
    vector<vector<int> > RiskPairs(ntime);
    vector<vector<vector<int> > > RiskPairs_Strata(ntime, vector<vector<int>>(Strata_vals.size()));
    //
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
    if (model_bool["null"]) {
        if (model_bool["strata"]) {
            Rls1 = MatrixXd::Zero(ntime, Strata_vals.size());  //  precomputes a series of sums used frequently in the log-liklihood calculations
            Lls1 = MatrixXd::Zero(ntime, Strata_vals.size());
            Calculate_Null_Sides_Strata(RiskFail, RiskPairs_Strata, ntime, R, Rls1, Lls1, Strata_vals, nthreads);
            Calc_Null_LogLik_Strata(nthreads, RiskFail, RiskPairs_Strata, ntime, R, Rls1, Lls1, Strata_vals, Ll, ties_method);
        } else {
            Rls1 = MatrixXd::Zero(ntime, 1);  //  precomputes a series of sums used frequently in the log-liklihood calculations
            Lls1 = MatrixXd::Zero(ntime, 1);
            //  the log-likelihood is calculated in parallel over the risk groups
            //
            Calculate_Null_Sides(RiskFail, RiskPairs, ntime, R, Rls1, Lls1, nthreads);
            Calc_Null_LogLik(nthreads, RiskFail, RiskPairs, ntime, R, Rls1, Lls1, Ll, ties_method);
        }
        //
        List res_list = List::create(_["LogLik"] = wrap(Ll[0]), _["AIC"] = -2*Ll[0], _["BIC"] = -2*Ll[0], _["Status"] = "PASSED");
        //  returns a list of results
        return res_list;
    }
    Cox_Refresh_R_SIDES(reqrdnum, ntime, Rls1, Rls2, Rls3, Lls1, Lls2, Lls3, Strata_vals, model_bool);
    //  the log-likelihood is calculated in parallel over the risk groups
    vector <double> Ll_comp(2, Ll[0]);  //  vector to compare values
    double step_max0 = step_max;
    double thres_step_max0 = thres_step_max;
    //
    vector<double> dbeta(totalnum, 0.0);
    NumericVector m_g_store(reqrdnum);
    NumericVector v_beta_store(reqrdnum);
    //
    //  --------------------------
    //  always starts from initial guess
    //  --------------------------
    vector<double> beta_c(totalnum, 0.0);
    vector<double> beta_a(totalnum, 0.0);
    vector<double> beta_best(totalnum, 0.0);
    vector<double> beta_p(totalnum, 0.0);
    VectorXd::Map(&beta_p[0], beta_0.size()) = beta_0;  //  stores previous parameters
    VectorXd::Map(&beta_c[0], beta_0.size()) = beta_0;  //  stores current parameters
    VectorXd::Map(&beta_a[0], beta_0.size()) = beta_0;  //  stores a refrence value for parameters
    VectorXd::Map(&beta_best[0], beta_0.size()) = beta_0;  //  stores the best parameters
    double halves = 0;  //  number of half-steps taken
    int ind0 = fir;  //  used for validations
    int iteration = 0;  //  iteration number
    int maxiter = 0;
    //
    bool convgd = FALSE;
    int iter_stop  = 0;  //  tracks if the iterations should be stopped for convergence
    //
    NumericMatrix beta_fin(a_ns.rows(), a_ns.cols());
    NumericVector LL_fin(a_ns.rows());
    //
    double Ll_abs_best = 10;
    vector<double> beta_abs_best(totalnum, 0.0);
    int guess_abs_best = - 1;
    double Ll_iter_best = 10;
    ///
    //  Variables that are used for the risk check function shared across cox, poisson, and log bound functions
    double dev = 0.0;
    MatrixXd dev_temp = MatrixXd::Zero(1, 1);
    double Lstar = 0.0;
    MatrixXd PyrC = MatrixXd::Zero(1, 1);
    //
    for (int guess = 0; guess <guesses; guess++) {
        Cox_Refresh_R_TERM(totalnum, reqrdnum, term_tot, dint, dslp, thres_step_max, step_max, df0, T0, Td0, Tdd0, Te, R, Rd, Rdd, Dose, nonDose, TTerm, nonDose_LIN, nonDose_PLIN, nonDose_LOGLIN, RdR, RddR, model_bool);
        Cox_Refresh_R_SIDES(reqrdnum, ntime, Rls1, Rls2, Rls3, Lls1, Lls2, Lls3, Strata_vals, model_bool);
        fill(Ll.begin(), Ll.end(), 0.0);
        if (!model_bool["single"]) {
            fill(Lld.begin(), Lld.end(), 0.0);
            fill(Lldd.begin(), Lldd.end(), 0.0);
        }
        if (model_bool["gradient"]) {
            m_g_store.fill(0);
            v_beta_store.fill(0);
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
        Ll_iter_best = 10;
        //
        convgd = FALSE;
        iter_stop  = 0;  //  tracks if the iterations should be stopped for convergence
        //
        maxiter = maxiters[guess];
        a_n = a_ns.row(guess);
        for (int i = 0; i < beta_0.size(); i++) {
            beta_0[i] = a_n[i];
        }
        if (verbose >= 4) {
            Rcout << "C++ Note: starting guess " << guess << endl;
        }
        Cox_Term_Risk_Calc(modelform, tform, term_n, totalnum, fir, dfc, term_tot, T0, Td0, Tdd0, Te, R, Rd, Rdd, Dose, nonDose, beta_0, df0, thres_step_max, step_max, TTerm, nonDose_LIN, nonDose_PLIN, nonDose_LOGLIN, RdR, RddR, nthreads, KeepConstant, verbose, model_bool, gmix_theta, gmix_term);
        if ((R.minCoeff() <= 0) || (R.hasNaN())) {
            if (verbose >= 1) {
                Rcout << "C++ Error: A non-positive risk was detected: " << R.minCoeff() << endl;
                Rcout << "C++ Warning: final failing values ";
                for (int ijk = 0; ijk < totalnum; ijk++) {
                    Rcout << beta_0[ijk] << " ";
                }
                Rcout << " " << endl;
            }
            //
            temp_list = List::create(_["beta_0"] = wrap(beta_0), _["Deviation"] = R_NaN, _["Status"] = "FAILED_WITH_NEGATIVE_RISK", _["LogLik"] = R_NaN);
            return temp_list;
        }
        //
        //  -------------------------------------------------------------------------------------------
        Cox_Side_LL_Calc(reqrdnum, ntime, tform, RiskFail, RiskPairs, RiskPairs_Strata, totalnum, fir, R, Rd, Rdd, Rls1, Rls2, Rls3, Lls1, Lls2, Lls3, cens_weight, Strata_vals, beta_0, RdR, RddR, Ll, Lld, Lldd, nthreads, KeepConstant, ties_method, verbose, model_bool, iter_stop);
        //
        for (int i = 0; i < beta_0.size(); i++) {
            beta_c[i] = beta_0[i];
        }
        while ((iteration < maxiter) && (iter_stop == 0)) {
            iteration++;
            beta_p = beta_c;  //
            beta_a = beta_c;  //
            beta_best = beta_c;  //
            //
            //  calculates the initial change in parameter
            if (model_bool["gradient"]) {
                if (model_bool["constraint"]) {
                    Calc_Change_Gradient_Cons(Lin_Sys, Lin_Res, nthreads, model_bool, totalnum, optim_para, iteration, step_max, Ll, Lld, m_g_store, v_beta_store, beta_0, dbeta, KeepConstant);
                } else {
                    Calc_Change_Gradient(nthreads, model_bool, totalnum, optim_para, iteration, step_max, Lld, m_g_store, v_beta_store, dbeta, KeepConstant);
                }
            } else if (model_bool["basic"]) {
                if (model_bool["constraint"]) {
                    Calc_Change_Basic_Cons(Lin_Sys, Lin_Res, beta_0, nthreads, totalnum, lr, step_max, Ll, Lld, Lldd, dbeta, KeepConstant);
                } else {
                    Calc_Change_Basic(nthreads, totalnum, lr, step_max, Ll, Lld, Lldd, dbeta, KeepConstant);
                }
            } else {
                if (model_bool["constraint"]) {
                    Calc_Change_Cons(Lin_Sys, Lin_Res, beta_0, nthreads, totalnum, thres_step_max, lr, step_max, Ll, Lld, Lldd, dbeta, tform, thres_step_max, step_max, KeepConstant);
                } else {
                    Calc_Change(nthreads, totalnum, thres_step_max, lr, step_max, Ll, Lld, Lldd, dbeta, tform, thres_step_max, step_max, KeepConstant);
                }
                Intercept_Bound(nthreads, totalnum, beta_0, dbeta, dfc, df0, KeepConstant, tform);
            }
            if ((Ll_iter_best > 0) || (Ll_iter_best < Ll[ind0])) {
                Ll_iter_best = Ll[ind0];
            }
            //
            if (model_bool["gradient"]) {
                //
                for (int ij = 0; ij < totalnum; ij++) {
                    beta_0[ij] = beta_a[ij] + dbeta[ij];
                    beta_c[ij] = beta_0[ij];
                }
                //
                Cox_Term_Risk_Calc(modelform, tform, term_n, totalnum, fir, dfc, term_tot, T0, Td0, Tdd0, Te, R, Rd, Rdd, Dose, nonDose, beta_0, df0, thres_step_max, step_max, TTerm, nonDose_LIN, nonDose_PLIN, nonDose_LOGLIN, RdR, RddR, nthreads, KeepConstant, verbose, model_bool, gmix_theta, gmix_term);
                //
                Cox_Pois_Check_Continue(model_bool, beta_0, beta_best, beta_c, cens_weight, dbeta, dev, dev_temp, fir, halfmax, halves, ind0, iter_stop, KeepConstant, Ll, Ll_iter_best, Lld, Lldd, Lls1, Lls2, Lls3, Lstar, nthreads, ntime, PyrC, R, Rd, Rdd, RddR, RdR, reqrdnum, tform, RiskFail,  RiskPairs, RiskPairs_Strata, Rls1, Rls2, Rls3, Strata_vals, term_n, ties_method, totalnum, TTerm, verbose);
            } else {
                halves = 0;
                while ((Ll[ind0] <= Ll_iter_best) && (halves < halfmax)) {  //  repeats until half-steps maxed or an improvement
                    for (int ij = 0; ij < totalnum; ij++) {
                        beta_0[ij] = beta_a[ij] + dbeta[ij];
                        beta_c[ij] = beta_0[ij];
                    }
                    //  ----------------------------------------------------------------------------------------------------------//
                    //  The same subterm, risk, sides, and log-likelihood calculations are performed every half-step and iteration
                    //  ----------------------------------------------------------------------------------------------------------//
                    Cox_Term_Risk_Calc(modelform, tform, term_n, totalnum, fir, dfc, term_tot, T0, Td0, Tdd0, Te, R, Rd, Rdd, Dose, nonDose, beta_0, df0, thres_step_max, step_max, TTerm, nonDose_LIN, nonDose_PLIN, nonDose_LOGLIN, RdR, RddR, nthreads, KeepConstant, verbose, model_bool, gmix_theta, gmix_term);
                    //
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
            }
            dbeta_max = abs(dbeta[0]);
            for (int ij = 1; ij < totalnum; ij++) {
                if (abs(dbeta[ij]) > dbeta_max) {
                    dbeta_max = abs(dbeta[ij]);
                }
            }
            if (Lld_worst < deriv_epsilon) {  //  ends if the derivatives are low enough
                iter_stop = 1;
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
        a_n = beta_0;
        beta_fin(guess, _) = a_n;
        LL_fin[guess] = Ll[0];
        if ((Ll_abs_best > 0) || (Ll_abs_best < Ll[ind0])) {
            Ll_abs_best = Ll[ind0];
            beta_abs_best = beta_c;
            guess_abs_best = guess;
        }
        //
    }
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
    convgd = FALSE;
    iter_stop  = 0;  //  tracks if the iterations should be stopped for convergence
    //
    int guess_max = guess_abs_best;
    if (verbose >= 3) {
        Rcout << "C++ Note: Guess Results" << endl;
        Rcout << "Guess number, parameter values, Log-Likelihood" << endl;
        NumericVector beta_temp;
        for (int i = 0;  i < guesses; i++) {
            beta_temp = wrap(beta_fin.row(i));
            if (i == guess_max) {
                Rcout << i << ", " << beta_temp << ", " << LL_fin[i] << "<-- Best Guess" << endl;
            } else {
                Rcout << i << ", " << beta_temp << ", " << LL_fin[i] << endl;
            }
        }
    }
    //
    maxiter = maxiters[guesses];
    a_n = beta_abs_best;
    for (int i = 0; i < beta_0.size(); i++) {
        beta_0[i] = a_n[i];
    }
    for (int i = 0; i < beta_0.size(); i++) {
        beta_c[i] = beta_0[i];
    }
    if (model_bool["gradient"]) {
        m_g_store.fill(0);
        v_beta_store.fill(0);
    }
    //
    //  Calculates the subterm and term values
    Cox_Term_Risk_Calc(modelform, tform, term_n, totalnum, fir, dfc, term_tot, T0, Td0, Tdd0, Te, R, Rd, Rdd, Dose, nonDose, beta_0, df0, thres_step_max, step_max, TTerm, nonDose_LIN, nonDose_PLIN, nonDose_LOGLIN, RdR, RddR, nthreads, KeepConstant, verbose, model_bool, gmix_theta, gmix_term);
    //
    //  -------------------------------------------------------------------------------------------
    //  Calculates the side sum terms used
    Cox_Side_LL_Calc(reqrdnum, ntime, tform, RiskFail,  RiskPairs, RiskPairs_Strata, totalnum, fir, R, Rd, Rdd, Rls1, Rls2, Rls3, Lls1, Lls2, Lls3, cens_weight, Strata_vals, beta_0, RdR, RddR, Ll, Lld, Lldd, nthreads, KeepConstant, ties_method, verbose, model_bool, iter_stop);
    List res_list;
    //
    if (model_bool["single"]) {
        res_list = List::create(_["LogLik"] = wrap(Ll[0]), _["beta_0"] = wrap(beta_0), _["AIC"] = 2*(totalnum-accumulate(KeepConstant.begin(), KeepConstant.end(), 0.0))-2*Ll[0], _["BIC"] = (totalnum-accumulate(KeepConstant.begin(), KeepConstant.end(), 0.0))*log(mat_row)-2*Ll[0], _["Status"] = "PASSED");
        //  returns a list of results
        return res_list;
    }
    //
    for (int i = 0; i < beta_0.size(); i++) {
        beta_c[i] = beta_0[i];
    }
    while ((iteration < maxiter) && (iter_stop == 0)) {
        iteration++;
        beta_p = beta_c;  //
        beta_a = beta_c;  //
        beta_best = beta_c;  //
        //
        //  calculates the initial change in parameter
        if (model_bool["gradient"]) {
                if (model_bool["constraint"]) {
                    Calc_Change_Gradient_Cons(Lin_Sys, Lin_Res, nthreads, model_bool, totalnum, optim_para, iteration, step_max, Ll, Lld, m_g_store, v_beta_store, beta_0, dbeta, KeepConstant);
                } else {
                    Calc_Change_Gradient(nthreads, model_bool, totalnum, optim_para, iteration, step_max, Lld, m_g_store, v_beta_store, dbeta, KeepConstant);
                }
        } else if (model_bool["basic"]) {
            if (model_bool["constraint"]) {
                Calc_Change_Basic_Cons(Lin_Sys, Lin_Res, beta_0, nthreads, totalnum, lr, step_max, Ll, Lld, Lldd, dbeta, KeepConstant);
            } else {
                Calc_Change_Basic(nthreads, totalnum, lr, step_max, Ll, Lld, Lldd, dbeta, KeepConstant);
            }
        } else {
            if (model_bool["constraint"]) {
                Calc_Change_Cons(Lin_Sys, Lin_Res, beta_0, nthreads, totalnum, thres_step_max, lr, step_max, Ll, Lld, Lldd, dbeta, tform, thres_step_max, step_max, KeepConstant);
            } else {
                Calc_Change(nthreads, totalnum, thres_step_max, lr, step_max, Ll, Lld, Lldd, dbeta, tform, thres_step_max, step_max, KeepConstant);
            }
            Intercept_Bound(nthreads, totalnum, beta_0, dbeta, dfc, df0, KeepConstant, tform);
        }
        //
        if ((Ll_abs_best > 0) || (Ll_abs_best < Ll[ind0])) {
            Ll_abs_best = Ll[ind0];
            beta_abs_best = beta_c;
        }
        //
        if (model_bool["gradient"]) {
            //
            for (int ij = 0; ij < totalnum; ij++) {
                beta_0[ij] = beta_a[ij] + dbeta[ij];
                beta_c[ij] = beta_0[ij];
            }
            //
            Cox_Term_Risk_Calc(modelform, tform, term_n, totalnum, fir, dfc, term_tot, T0, Td0, Tdd0, Te, R, Rd, Rdd, Dose, nonDose, beta_0, df0, thres_step_max, step_max, TTerm, nonDose_LIN, nonDose_PLIN, nonDose_LOGLIN, RdR, RddR, nthreads, KeepConstant, verbose, model_bool, gmix_theta, gmix_term);
            //
            Cox_Pois_Check_Continue(model_bool, beta_0, beta_best, beta_c, cens_weight, dbeta, dev, dev_temp, fir, halfmax, halves, ind0, iter_stop, KeepConstant, Ll, Ll_abs_best, Lld, Lldd, Lls1, Lls2, Lls3, Lstar, nthreads, ntime, PyrC, R, Rd, Rdd, RddR, RdR, reqrdnum, tform, RiskFail,  RiskPairs, RiskPairs_Strata, Rls1, Rls2, Rls3, Strata_vals, term_n, ties_method, totalnum, TTerm, verbose);
        } else {
            halves = 0;
            while ((Ll[ind0] <= Ll_abs_best) && (halves < halfmax)) {  //  repeats until half-steps maxed or an improvement
                for (int ij = 0; ij < totalnum; ij++) {
                    beta_0[ij] = beta_a[ij] + dbeta[ij];
                    beta_c[ij] = beta_0[ij];
                }
                //  ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------//
                //  The same subterm, risk, sides, and log-likelihood calculations are performed every half-step and iteration
                //  ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------//
                Cox_Term_Risk_Calc(modelform, tform, term_n, totalnum, fir, dfc, term_tot, T0, Td0, Tdd0, Te, R, Rd, Rdd, Dose, nonDose, beta_0, df0, thres_step_max, step_max, TTerm, nonDose_LIN, nonDose_PLIN, nonDose_LOGLIN, RdR, RddR, nthreads, KeepConstant, verbose, model_bool, gmix_theta, gmix_term);
                Cox_Pois_Check_Continue(model_bool, beta_0, beta_best, beta_c, cens_weight, dbeta, dev, dev_temp, fir, halfmax, halves, ind0, iter_stop, KeepConstant, Ll, Ll_abs_best, Lld, Lldd, Lls1, Lls2, Lls3, Lstar, nthreads, ntime, PyrC, R, Rd, Rdd, RddR, RdR, reqrdnum, tform, RiskFail,  RiskPairs, RiskPairs_Strata, Rls1, Rls2, Rls3, Strata_vals, term_n, ties_method, totalnum, TTerm, verbose);
            }
            if (beta_best != beta_c) {  //  if the risk matrices aren't the optimal values, then they must be recalculated
                //  If it goes through every half step without improvement, then the maximum change needs to be decreased
                step_max = step_max*pow(0.5, halfmax);  //  reduces the step sizes
                thres_step_max = thres_step_max*pow(0.5, halfmax);
                beta_p = beta_best;  //
                beta_a = beta_best;  //
                beta_c = beta_best;  //
                for (int ij = 0; ij < totalnum; ij++) {
                    beta_0[ij] = beta_best[ij];
                }
                Cox_Term_Risk_Calc(modelform, tform, term_n, totalnum, fir, dfc, term_tot, T0, Td0, Tdd0, Te, R, Rd, Rdd, Dose, nonDose, beta_0, df0, thres_step_max, step_max, TTerm, nonDose_LIN, nonDose_PLIN, nonDose_LOGLIN, RdR, RddR, nthreads, KeepConstant, verbose, model_bool, gmix_theta, gmix_term);
                //
            }
        }
        Cox_Side_LL_Calc(reqrdnum, ntime, tform, RiskFail,  RiskPairs, RiskPairs_Strata, totalnum, fir, R, Rd, Rdd, Rls1, Rls2, Rls3, Lls1, Lls2, Lls3, cens_weight, Strata_vals, beta_0, RdR, RddR, Ll, Lld, Lldd, nthreads, KeepConstant, ties_method, verbose, model_bool, iter_stop);
        Lld_worst = 0;
        for (int ij = 0; ij < reqrdnum; ij++) {
            if (abs(Lld[ij]) > Lld_worst) {
                Lld_worst = abs(Lld[ij]);
            }
        }
        dbeta_max = abs(dbeta[0]);
        for (int ij = 1; ij < totalnum; ij++) {
            if (abs(dbeta[ij]) > dbeta_max) {
                dbeta_max = abs(dbeta[ij]);
            }
        }
        if (Lld_worst < deriv_epsilon) {  //  ends if the derivatives are low enough
            iter_stop = 1;
            convgd = TRUE;
        }
        Ll_comp[1] = Ll[0];
        if (step_max < epsilon) {  //  if the maximum change is too low, then it ends
            iter_stop = 1;
        }
        if (dbeta_max < epsilon) {  //  if the maximum change is too low, then it ends
            iter_stop = 1;
        }
    }
    if (Lld_worst < deriv_epsilon) {  //  ends if the derivatives are low enough
        iter_stop = 1;
        convgd = TRUE;
    }
    //  -----------------------------------------------
    //  Performing Full Calculation to get full second derivative matrix
    //  -----------------------------------------------
    fill(Ll.begin(), Ll.end(), 0.0);
    if (!model_bool["single"]) {
        fill(Lld.begin(), Lld.end(), 0.0);
        fill(Lldd.begin(), Lldd.end(), 0.0);
    }
    bool true_gradient = model_bool["gradient"];
    model_bool["gradient"] = false;
    Cox_Term_Risk_Calc(modelform, tform, term_n, totalnum, fir, dfc, term_tot, T0, Td0, Tdd0, Te, R, Rd, Rdd, Dose, nonDose, beta_0, df0, thres_step_max, step_max, TTerm, nonDose_LIN, nonDose_PLIN, nonDose_LOGLIN, RdR, RddR, nthreads, KeepConstant, verbose, model_bool, gmix_theta, gmix_term);
    Cox_Side_LL_Calc(reqrdnum, ntime, tform, RiskFail,  RiskPairs, RiskPairs_Strata, totalnum, fir, R, Rd, Rdd, Rls1, Rls2, Rls3, Lls1, Lls2, Lls3, cens_weight, Strata_vals, beta_0, RdR, RddR, Ll, Lld, Lldd, nthreads, KeepConstant, ties_method, verbose, model_bool, iter_stop);
    model_bool["gradient"] = true_gradient;
    //
    if ((Ll_abs_best > 0) || (Ll_abs_best < Ll[ind0])) {
        Ll_abs_best = Ll[ind0];
        beta_abs_best = beta_c;
    }
    //
    List para_list;
    if (!model_bool["basic"]) {
        para_list = List::create(_["term_n"] = term_n, _["tforms"] = tform);  //  stores the term information
    }
    List control_list = List::create(_["Iteration"] = iteration, _["Maximum Step"] = dbeta_max, _["Derivative Limiting"] = Lld_worst);  //  stores the total number of iterations used
    //
    NumericVector Lldd_vec(reqrdnum * reqrdnum);  //  simplfied information matrix
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
    }
    Lldd_vec.attr("dim") = Dimension(reqrdnum, reqrdnum);
    const Map<MatrixXd> Lldd_mat(as<Map<MatrixXd> >(Lldd_vec));
    //
    MatrixXd cov;
    VectorXd stdev = VectorXd::Zero(totalnum);
    if (model_bool["observed_info"]) {
        cov = - 1 * Lldd_mat.inverse().matrix();  //  uses inverse information matrix to calculate the standard deviation
        for (int ij = 0; ij < totalnum; ij++) {
            if (KeepConstant[ij] == 0) {
                int pij_ind = ij - sum(head(KeepConstant, ij));
                stdev(ij) = sqrt(cov(pij_ind, pij_ind));
            }
        }
    } else {
        vector<double> InMa(pow(reqrdnum, 2), 0.0);
        if (model_bool["strata"]) {
            if (model_bool["cr"]) {
                Expected_Inform_Matrix_Cox_Strata_CR(nthreads, RiskFail, RiskPairs_Strata, totalnum, ntime, R, Rd, RdR, cens_weight, InMa, Strata_vals, KeepConstant);
            } else {
                Expected_Inform_Matrix_Cox_Strata(nthreads, RiskFail, RiskPairs_Strata, totalnum, ntime, R, Rd, RdR, InMa, Strata_vals, KeepConstant);
            }
        } else {
            if (model_bool["cr"]) {
                Expected_Inform_Matrix_Cox_CR(nthreads, RiskFail, RiskPairs, totalnum, ntime, R, Rd, RdR, cens_weight, InMa, KeepConstant);
            } else {
                Expected_Inform_Matrix_Cox(nthreads, RiskFail, RiskPairs, totalnum, ntime, R, Rd, RdR, InMa, KeepConstant);
            }
        }
        NumericVector InMa_vec(reqrdnum * reqrdnum);  //  simplfied information matrix
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
            InMa_vec[ij * reqrdnum + jk] = InMa[ij*reqrdnum+jk];
            InMa_vec[jk * reqrdnum + ij] = InMa[ij * reqrdnum + jk];
        }
        InMa_vec.attr("dim") = Dimension(reqrdnum, reqrdnum);
        const Map<MatrixXd> InMa_mat(as<Map<MatrixXd> >(InMa_vec));    //
        cov = InMa_mat.inverse().matrix();  //  uses inverse information matrix to calculate the standard deviation
        for (int ij = 0; ij < totalnum; ij++) {
            if (KeepConstant[ij] == 0) {
                int pij_ind = ij - sum(head(KeepConstant, ij));
                stdev(ij) = sqrt(cov(pij_ind, pij_ind));
            }
        }
    }
    //
    if (model_bool["basic"]) {
        res_list = List::create(_["LogLik"] = wrap(Ll[0]), _["First_Der"] = wrap(Lld), _["Second_Der"] = Lldd_vec, _["beta_0"] = wrap(beta_0), _["Standard_Deviation"] = wrap(stdev), _["Covariance"] = wrap(cov), _["AIC"] = 2*(totalnum-accumulate(KeepConstant.begin(), KeepConstant.end(), 0.0))-2*Ll[0], _["BIC"] = (totalnum-accumulate(KeepConstant.begin(), KeepConstant.end(), 0.0))*log(df0.rows())-2*Ll[0], _["Control_List"] = control_list, _["Converged"] = convgd, _["Status"] = "PASSED");
    } else {
        res_list = List::create(_["LogLik"] = wrap(Ll[0]), _["First_Der"] = wrap(Lld), _["Second_Der"] = Lldd_vec, _["beta_0"] = wrap(beta_0), _["Standard_Deviation"] = wrap(stdev), _["Covariance"] = wrap(cov), _["AIC"] = 2*(totalnum-accumulate(KeepConstant.begin(), KeepConstant.end(), 0.0))-2*Ll[0], _["BIC"] = (totalnum-accumulate(KeepConstant.begin(), KeepConstant.end(), 0.0))*log(df0.rows())-2*Ll[0], _["Parameter_Lists"] = para_list, _["Control_List"] = control_list, _["Converged"] = convgd, _["Status"] = "PASSED");
    }
    //  returns a list of results
    return res_list;
}

//' Primary poisson regression with multiple starting points and optional combinations of stratification and no derivative calculation.
//'
//' \code{LogLik_Pois_Omnibus} Performs the calls to calculation functions, Structures the poisson regression, With verbose option prints out time stamps and intermediate sums of terms and derivatives
//'
//' @inheritParams CPP_template
//'
//' @return List of final results: Log-likelihood of optimum, first derivative of log-likelihood, second derivative matrix, parameter list, standard deviation estimate, AIC, model information
//' @noRd
//'
//
List LogLik_Pois_Omnibus(const Ref<const MatrixXd>& PyrC, IntegerVector term_n, StringVector tform, NumericMatrix& a_ns, Ref<MatrixXd> df0, IntegerVector dfc, int fir, string modelform, double lr, List optim_para, NumericVector maxiters, int guesses, int halfmax, double epsilon, double step_max, double thres_step_max, double deriv_epsilon, int verbose, IntegerVector KeepConstant, int term_tot, int nthreads, const Ref<const MatrixXd>& dfs, List model_bool, const double gmix_theta, const IntegerVector gmix_term, const Ref<const MatrixXd>& Lin_Sys, const Ref<const VectorXd>& Lin_Res) {
    //
    List temp_list = List::create(_["Status"] = "FAILED");  //  used as a dummy return value for code checking
    //
    //  Time durations are measured from this point on in microseconds
    //
    //  df0: covariate data
    //  ntime: number of event times for Cox PH
    //  totalnum: number of terms used
    //
    //  ------------------------------------------------------------------------- //  initialize
    const int mat_row = df0.rows();
    //
    int totalnum = term_n.size();
    int reqrdnum = totalnum - sum(KeepConstant);
    //  cout.precision: controls the number of significant digits printed
    //  nthreads: number of threads used for parallel operations
    //
    Rcout.precision(7);  //  forces higher precision numbers printed to terminal
    //
    //  Lld_worst: stores the highest magnitude log-likelihood derivative
    //
    //
    double Lld_worst = 0.0;  //  stores derivative value used to determine if every parameter is near convergence
    double dbeta_max = 0.0;  //  stores the largest step taken, determines if the step sizes are small enough for convergence
    //
    //  ---------------------------------------------
    //  To Start, needs to seperate the derivative terms
    //  ---------------------------------------------
    //
    //  ------------------------------------------------------------------------- //  initialize
    NumericVector a_n = a_ns.row(0);
    Map<VectorXd> beta_0(as<Map<VectorXd> >(a_n));
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
    double dev = 0.0;
    //  the log-likelihood is calculated in parallel over the risk groups
    vector <double> Ll_comp(2, Ll[0]);  //  vector to compare values
    double step_max0 = step_max;
    double thres_step_max0 = thres_step_max;
    //
    vector<double> dbeta(totalnum, 0.0);
    NumericVector m_g_store(reqrdnum);
    NumericVector v_beta_store(reqrdnum);
    //
    //  --------------------------
    //  always starts from initial guess
    //  --------------------------
    vector<double> beta_c(totalnum, 0.0);
    vector<double> beta_a(totalnum, 0.0);
    vector<double> beta_best(totalnum, 0.0);
    vector<double> beta_p(totalnum, 0.0);
    VectorXd::Map(&beta_p[0], beta_0.size()) = beta_0;  //  stores previous parameters
    VectorXd::Map(&beta_c[0], beta_0.size()) = beta_0;  //  stores current parameters
    VectorXd::Map(&beta_a[0], beta_0.size()) = beta_0;  //  stores a refrence value for parameters
    VectorXd::Map(&beta_best[0], beta_0.size()) = beta_0;  //  stores the best parameters
    double halves = 0;  //  number of half-steps taken
    int ind0 = fir;  //  used for validations
    int iteration = 0;  //  iteration number
    int maxiter = 0;
    //
    bool convgd = FALSE;
    int iter_stop  = 0;  //  tracks if the iterations should be stopped for convergence
    //
    NumericMatrix beta_fin(a_ns.rows(), a_ns.cols());
    NumericVector LL_fin(a_ns.rows());
    //
    double Ll_abs_best = 10;
    vector<double> beta_abs_best(totalnum, 0.0);
    int guess_abs_best = - 1;
    double Ll_iter_best = 10;
    //  Variables that are used for the risk check function shared across cox, poisson, and log bound functions
    VectorXd cens_weight(1);
    MatrixXd Lls1 = MatrixXd::Zero(1, 1);
    MatrixXd Lls2 = MatrixXd::Zero(1, 1);
    MatrixXd Lls3 = MatrixXd::Zero(1, 1);
    double Lstar = 0.0;
    int ntime = 1.0;
    IntegerMatrix RiskFail(1);
    vector<vector<int> > RiskPairs;
    vector<vector<vector<int> > > RiskPairs_Strata;
    MatrixXd Rls1 = MatrixXd::Zero(1, 1);
    MatrixXd Rls2 = MatrixXd::Zero(1, 1);
    MatrixXd Rls3 = MatrixXd::Zero(1, 1);
    NumericVector Strata_vals(1);
    string ties_method = "temp";
    ///
    for (int guess = 0; guess <guesses; guess++) {
        fill(Ll.begin(), Ll.end(), 0.0);
        if (!model_bool["single"]) {
            fill(Lld.begin(), Lld.end(), 0.0);
            fill(Lldd.begin(), Lldd.end(), 0.0);
        }
        if (model_bool["gradient"]) {
            m_g_store.fill(0);
            v_beta_store.fill(0);
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
        Ll_iter_best = 10;
        //
        convgd = FALSE;
        iter_stop  = 0;  //  tracks if the iterations should be stopped for convergence
        //
        maxiter = maxiters[guess];
        a_n = a_ns.row(guess);
        for (int i = 0; i < beta_0.size(); i++) {
            beta_0[i] = a_n[i];
        }
        //
        Pois_Term_Risk_Calc(modelform, tform, term_n, totalnum, fir, dfc, term_tot, T0, Td0, Tdd0, Te, R, Rd, Rdd, Dose, nonDose, beta_0, df0, dint, dslp, TTerm, nonDose_LIN, nonDose_PLIN, nonDose_LOGLIN, RdR, RddR, s_weights, nthreads, KeepConstant, verbose, model_bool, gmix_theta, gmix_term);
        if ((R.minCoeff() <= 0) || (R.hasNaN())) {
            if (verbose >= 1) {
                Rcout << "C++ Error: A non-positive risk was detected: " << R.minCoeff() << endl;
                Rcout << "C++ Warning: final failing values ";
                for (int ijk = 0; ijk < totalnum; ijk++) {
                    Rcout << beta_0[ijk] << " ";
                }
                Rcout << " " << endl;
                //
                Rcout << "C++ Warning: final failing terms ";
                for (int ijk = 0; ijk < totalnum; ijk++) {
                    Rcout << tform[ijk] << " ";
                }
                Rcout << " " << endl;
            }
            temp_list = List::create(_["beta_0"] = wrap(beta_0), _["Deviation"] = R_NaN, _["Status"] = "FAILED_WITH_NEGATIVE_RISK", _["LogLik"] = R_NaN);
            return temp_list;
        }
        //
        //  -------------------------------------------------------------------------------------------
        //
        //  Calculates log-likelihood
        Pois_Dev_LL_Calc(reqrdnum, totalnum, fir, R, Rd, Rdd, beta_0, RdR, RddR, Ll, Lld, Lldd, PyrC, dev_temp, nthreads, KeepConstant, verbose, model_bool, iter_stop, dev);
        //
        for (int i = 0; i < beta_0.size(); i++) {
            beta_c[i] = beta_0[i];
        }
        while ((iteration < maxiter) && (iter_stop == 0)) {
            iteration++;
            beta_p = beta_c;  //
            beta_a = beta_c;  //
            beta_best = beta_c;  //
            //
            //  calculates the initial change in parameter
            if (model_bool["gradient"]) {
                if (model_bool["constraint"]) {
                    Calc_Change_Gradient_Cons(Lin_Sys, Lin_Res, nthreads, model_bool, totalnum, optim_para, iteration, step_max, Ll, Lld, m_g_store, v_beta_store, beta_0, dbeta, KeepConstant);
                } else {
                    Calc_Change_Gradient(nthreads, model_bool, totalnum, optim_para, iteration, step_max, Lld, m_g_store, v_beta_store, dbeta, KeepConstant);
                }
            } else if (model_bool["constraint"]) {
                Calc_Change_Cons(Lin_Sys, Lin_Res, beta_0, nthreads, totalnum, thres_step_max, lr, step_max, Ll, Lld, Lldd, dbeta, tform, thres_step_max, step_max, KeepConstant);
            } else {
                Calc_Change(nthreads, totalnum, thres_step_max, lr, step_max, Ll, Lld, Lldd, dbeta, tform, thres_step_max, step_max, KeepConstant);
            }
            Intercept_Bound(nthreads, totalnum, beta_0, dbeta, dfc, df0, KeepConstant, tform);
            if ((Ll_iter_best > 0) || (Ll_iter_best < Ll[ind0])) {
                Ll_iter_best = Ll[ind0];
            }
            if (model_bool["gradient"]) {
                //
                for (int ij = 0; ij < totalnum; ij++) {
                    beta_0[ij] = beta_a[ij] + dbeta[ij];
                    beta_c[ij] = beta_0[ij];
                }
                //
                Pois_Term_Risk_Calc(modelform, tform, term_n, totalnum, fir, dfc, term_tot, T0, Td0, Tdd0, Te, R, Rd, Rdd, Dose, nonDose, beta_0, df0, dint, dslp, TTerm, nonDose_LIN, nonDose_PLIN, nonDose_LOGLIN, RdR, RddR, s_weights, nthreads, KeepConstant, verbose, model_bool, gmix_theta, gmix_term);
                //
                Cox_Pois_Check_Continue(model_bool, beta_0, beta_best, beta_c, cens_weight, dbeta, dev, dev_temp, fir, halfmax, halves, ind0, iter_stop, KeepConstant, Ll, Ll_iter_best, Lld, Lldd, Lls1, Lls2, Lls3, Lstar, nthreads, ntime, PyrC, R, Rd, Rdd, RddR, RdR, reqrdnum, tform, RiskFail,  RiskPairs, RiskPairs_Strata, Rls1, Rls2, Rls3, Strata_vals, term_n, ties_method, totalnum, TTerm, verbose);
            } else {
                halves = 0;
                while ((Ll[ind0] <= Ll_iter_best) && (halves < halfmax)) {  //  repeats until half-steps maxed or an improvement
                    for (int ij = 0; ij < totalnum; ij++) {
                        beta_0[ij] = beta_a[ij] + dbeta[ij];
                        beta_c[ij] = beta_0[ij];
                    }
                    //  -----------------------------------------------------------------------------------------------------------------------------------------//
                    //  The same subterm, risk, sides, and log-likelihood calculations are performed every half-step and iteration
                    //  -----------------------------------------------------------------------------------------------------------------------------------------//
                    Pois_Term_Risk_Calc(modelform, tform, term_n, totalnum, fir, dfc, term_tot, T0, Td0, Tdd0, Te, R, Rd, Rdd, Dose, nonDose, beta_0, df0, dint, dslp, TTerm, nonDose_LIN, nonDose_PLIN, nonDose_LOGLIN, RdR, RddR, s_weights, nthreads, KeepConstant, verbose, model_bool, gmix_theta, gmix_term);
                    Cox_Pois_Check_Continue(model_bool, beta_0, beta_best, beta_c, cens_weight, dbeta, dev, dev_temp, fir, halfmax, halves, ind0, iter_stop, KeepConstant, Ll, Ll_iter_best, Lld, Lldd, Lls1, Lls2, Lls3, Lstar, nthreads, ntime, PyrC, R, Rd, Rdd, RddR, RdR, reqrdnum, tform, RiskFail,  RiskPairs, RiskPairs_Strata, Rls1, Rls2, Rls3, Strata_vals, term_n, ties_method, totalnum, TTerm, verbose);
                }
                if (beta_best != beta_c) {  //  if the risk matrices aren't the optimal values, then they must be recalculated
                    //  If it goes through every half step without improvement, then the maximum change needs to be decreased
                    step_max = step_max*pow(0.5, halfmax);  //  reduces the step sizes
                    thres_step_max = thres_step_max*pow(0.5, halfmax);
                    beta_p = beta_best;  //
                    beta_a = beta_best;  //
                    beta_c = beta_best;  //
                    for (int ij = 0; ij < totalnum; ij++) {
                        beta_0[ij] = beta_best[ij];
                    }
                    Pois_Term_Risk_Calc(modelform, tform, term_n, totalnum, fir, dfc, term_tot, T0, Td0, Tdd0, Te, R, Rd, Rdd, Dose, nonDose, beta_0, df0, dint, dslp, TTerm, nonDose_LIN, nonDose_PLIN, nonDose_LOGLIN, RdR, RddR, s_weights, nthreads, KeepConstant, verbose, model_bool, gmix_theta, gmix_term);
                    //  Calculates log-likelihood
                    Pois_Dev_LL_Calc(reqrdnum, totalnum, fir, R, Rd, Rdd, beta_0, RdR, RddR, Ll, Lld, Lldd, PyrC, dev_temp, nthreads, KeepConstant, verbose, model_bool, iter_stop, dev);
                    //
                }
            }
            Lld_worst = 0;
            for (int ij = 0; ij < reqrdnum; ij++) {
                if (abs(Lld[ij]) > Lld_worst) {
                    Lld_worst = abs(Lld[ij]);
                }
            }
            dbeta_max = abs(dbeta[0]);
            for (int ij = 1; ij < totalnum; ij++) {
                if (abs(dbeta[ij]) > dbeta_max) {
                    dbeta_max = abs(dbeta[ij]);
                }
            }
            if (Lld_worst < deriv_epsilon) {  //  ends if the derivatives are low enough
                iter_stop = 1;
            }
            Ll_comp[1] = Ll[0];
            if (step_max < epsilon) {  //  if the maximum change is too low, then it ends
                iter_stop = 1;
            }
            if (dbeta_max < epsilon) {  //  if the maximum change is too low, then it ends
                iter_stop = 1;
            }
            if (model_bool["single"]) {
                iter_stop = 1;
            } else {}
        }
        //  -----------------------------------------------
        //  Performing Full Calculation to get full second derivative matrix
        //  -----------------------------------------------
        fill(Ll.begin(), Ll.end(), 0.0);
        if (!model_bool["single"]) {
            fill(Lld.begin(), Lld.end(), 0.0);
            fill(Lldd.begin(), Lldd.end(), 0.0);
        }
        //  Calculates log-likelihood
        Pois_Dev_LL_Calc(reqrdnum, totalnum, fir, R, Rd, Rdd, beta_0, RdR, RddR, Ll, Lld, Lldd, PyrC, dev_temp, nthreads, KeepConstant, verbose, model_bool, iter_stop, dev);
        //
        a_n = beta_0;
        //
        beta_fin(guess, _) = a_n;
        LL_fin[guess] = Ll[0];
        if ((Ll_abs_best > 0) || (Ll_abs_best < Ll[ind0])) {
            Ll_abs_best = Ll[ind0];
            beta_abs_best = beta_c;
            guess_abs_best = guess;
        }
    }
    //
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
    convgd = FALSE;
    iter_stop  = 0;  //  tracks if the iterations should be stopped for convergence
    //
    int guess_max = guess_abs_best;
    if (verbose >= 3) {
        Rcout << "C++ Note: Guess Results" << endl;
        Rcout << "Guess number, parameter values, Log-Likelihood" << endl;
        NumericVector beta_temp;
        for (int i = 0;  i < guesses; i++) {
            beta_temp = wrap(beta_fin.row(i));
            if (i == guess_max) {
                Rcout << i << ", " << beta_temp << ", " << LL_fin[i] << "<-- Best Guess" << endl;
            } else {
                Rcout << i << ", " << beta_temp << ", " << LL_fin[i] << endl;
            }
        }
    }
    //
    maxiter = maxiters[guesses];
    a_n = beta_abs_best;
    for (int i = 0; i < beta_0.size(); i++) {
        beta_0[i] = a_n[i];
    }
    for (int i = 0; i < beta_0.size(); i++) {
        beta_c[i] = beta_0[i];
    }
    if (model_bool["gradient"]) {
        m_g_store.fill(0);
        v_beta_store.fill(0);
    }
    //  Calculates the subterm and term values
    Pois_Term_Risk_Calc(modelform, tform, term_n, totalnum, fir, dfc, term_tot, T0, Td0, Tdd0, Te, R, Rd, Rdd, Dose, nonDose, beta_0, df0, dint, dslp, TTerm, nonDose_LIN, nonDose_PLIN, nonDose_LOGLIN, RdR, RddR, s_weights, nthreads, KeepConstant, verbose, model_bool, gmix_theta, gmix_term);
    //  Calculates log-likelihood
    Pois_Dev_LL_Calc(reqrdnum, totalnum, fir, R, Rd, Rdd, beta_0, RdR, RddR, Ll, Lld, Lldd, PyrC, dev_temp, nthreads, KeepConstant, verbose, model_bool, iter_stop, dev);
    //
    List res_list;
    //
    if (model_bool["single"]) {
        res_list = List::create(_["LogLik"] = wrap(Ll[0]), _["beta_0"] = wrap(beta_0), _["AIC"] = 2*(totalnum-accumulate(KeepConstant.begin(), KeepConstant.end(), 0.0))+dev, _["BIC"] = (totalnum-accumulate(KeepConstant.begin(), KeepConstant.end(), 0.0))*log(mat_row)-2*Ll[0], _["Status"] = "PASSED");
        //  returns a list of results
        return res_list;
    }
    for (int i = 0; i < beta_0.size(); i++) {
        beta_c[i] = beta_0[i];
    }
    while ((iteration < maxiter) && (iter_stop == 0)) {
        iteration++;
        beta_p = beta_c;  //
        beta_a = beta_c;  //
        beta_best = beta_c;  //
        //
        //  calculates the initial change in parameter
        if (model_bool["gradient"]) {
            if (model_bool["constraint"]) {
                Calc_Change_Gradient_Cons(Lin_Sys, Lin_Res, nthreads, model_bool, totalnum, optim_para, iteration, step_max, Ll, Lld, m_g_store, v_beta_store, beta_0, dbeta, KeepConstant);
            } else {
                Calc_Change_Gradient(nthreads, model_bool, totalnum, optim_para, iteration, step_max, Lld, m_g_store, v_beta_store, dbeta, KeepConstant);
            }
        } else if (model_bool["constraint"]) {
            Calc_Change_Cons(Lin_Sys, Lin_Res, beta_0, nthreads, totalnum, thres_step_max, lr, step_max, Ll, Lld, Lldd, dbeta, tform, thres_step_max, step_max, KeepConstant);
        } else {
            Calc_Change(nthreads, totalnum, thres_step_max, lr, step_max, Ll, Lld, Lldd, dbeta, tform, thres_step_max, step_max, KeepConstant);
        }
        Intercept_Bound(nthreads, totalnum, beta_0, dbeta, dfc, df0, KeepConstant, tform);
        //
        if ((Ll_abs_best > 0) || (Ll_abs_best < Ll[ind0])) {
            Ll_abs_best = Ll[ind0];
            beta_abs_best = beta_c;
        }
        //
        if (model_bool["gradient"]) {
            //
            for (int ij = 0; ij < totalnum; ij++) {
                beta_0[ij] = beta_a[ij] + dbeta[ij];
                beta_c[ij] = beta_0[ij];
            }
            //
            Pois_Term_Risk_Calc(modelform, tform, term_n, totalnum, fir, dfc, term_tot, T0, Td0, Tdd0, Te, R, Rd, Rdd, Dose, nonDose, beta_0, df0, dint, dslp, TTerm, nonDose_LIN, nonDose_PLIN, nonDose_LOGLIN, RdR, RddR, s_weights, nthreads, KeepConstant, verbose, model_bool, gmix_theta, gmix_term);
            //
            Cox_Pois_Check_Continue(model_bool, beta_0, beta_best, beta_c, cens_weight, dbeta, dev, dev_temp, fir, halfmax, halves, ind0, iter_stop, KeepConstant, Ll, Ll_abs_best, Lld, Lldd, Lls1, Lls2, Lls3, Lstar, nthreads, ntime, PyrC, R, Rd, Rdd, RddR, RdR, reqrdnum, tform, RiskFail,  RiskPairs, RiskPairs_Strata, Rls1, Rls2, Rls3, Strata_vals, term_n, ties_method, totalnum, TTerm, verbose);
        } else {
            halves = 0;
            while ((Ll[ind0] <= Ll_abs_best) && (halves < halfmax)) {  //  repeats until half-steps maxed or an improvement
                for (int ij = 0; ij < totalnum; ij++) {
                    beta_0[ij] = beta_a[ij] + dbeta[ij];
                    beta_c[ij] = beta_0[ij];
                }
                //  ----------------------------------------------------------------------------------------------------------------------------------------------------------//
                //  The same subterm, risk, sides, and log-likelihood calculations are performed every half-step and iteration
                //  ----------------------------------------------------------------------------------------------------------------------------------------------------------//
                Pois_Term_Risk_Calc(modelform, tform, term_n, totalnum, fir, dfc, term_tot, T0, Td0, Tdd0, Te, R, Rd, Rdd, Dose, nonDose, beta_0, df0, dint, dslp, TTerm, nonDose_LIN, nonDose_PLIN, nonDose_LOGLIN, RdR, RddR, s_weights, nthreads, KeepConstant, verbose, model_bool, gmix_theta, gmix_term);
                Cox_Pois_Check_Continue(model_bool, beta_0, beta_best, beta_c, cens_weight, dbeta, dev, dev_temp, fir, halfmax, halves, ind0, iter_stop, KeepConstant, Ll, Ll_abs_best, Lld, Lldd, Lls1, Lls2, Lls3, Lstar, nthreads, ntime, PyrC, R, Rd, Rdd, RddR, RdR, reqrdnum, tform, RiskFail,  RiskPairs, RiskPairs_Strata, Rls1, Rls2, Rls3, Strata_vals, term_n, ties_method, totalnum, TTerm, verbose);
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
                //  Calculates log-likelihood
                Pois_Dev_LL_Calc(reqrdnum, totalnum, fir, R, Rd, Rdd, beta_0, RdR, RddR, Ll, Lld, Lldd, PyrC, dev_temp, nthreads, KeepConstant, verbose, model_bool, iter_stop, dev);
            }
        }
        Lld_worst = 0;
        for (int ij = 0; ij < reqrdnum; ij++) {
            if (abs(Lld[ij]) > Lld_worst) {
                Lld_worst = abs(Lld[ij]);
            }
        }
        dbeta_max = abs(dbeta[0]);
        for (int ij = 1; ij < totalnum; ij++) {
            if (abs(dbeta[ij]) > dbeta_max) {
                dbeta_max = abs(dbeta[ij]);
            }
        }
        if (Lld_worst < deriv_epsilon) {  //  ends if the derivatives are low enough
            iter_stop = 1;
            convgd = TRUE;
        }
        Ll_comp[1] = Ll[0];
        if (step_max < epsilon) {  //  if the maximum change is too low, then it ends
            iter_stop = 1;
        }
        if (dbeta_max < epsilon) {  //  if the maximum change is too low, then it ends
            iter_stop = 1;
        }
        if (model_bool["single"]) {
            iter_stop = 1;
        }
    }
    if (Lld_worst < deriv_epsilon) {  //  ends if the derivatives are low enough
        iter_stop = 1;
        convgd = TRUE;
    }
    //  -----------------------------------------------
    //  Performing Full Calculation to get full second derivative matrix
    //  -----------------------------------------------
    bool true_gradient = model_bool["gradient"];
    model_bool["gradient"] = false;
    Pois_Term_Risk_Calc(modelform, tform, term_n, totalnum, fir, dfc, term_tot, T0, Td0, Tdd0, Te, R, Rd, Rdd, Dose, nonDose, beta_0, df0, dint, dslp, TTerm, nonDose_LIN, nonDose_PLIN, nonDose_LOGLIN, RdR, RddR, s_weights, nthreads, KeepConstant, verbose, model_bool, gmix_theta, gmix_term);
    Pois_Dev_LL_Calc(reqrdnum, totalnum, fir, R, Rd, Rdd, beta_0, RdR, RddR, Ll, Lld, Lldd, PyrC, dev_temp, nthreads, KeepConstant, verbose, model_bool, iter_stop, dev);
    model_bool["gradient"] = true_gradient;
    //
    if ((Ll_abs_best > 0) || (Ll_abs_best < Ll[ind0])) {
        Ll_abs_best = Ll[ind0];
        beta_abs_best = beta_c;
    }
    //
    List para_list = List::create(_["term_n"] = term_n, _["tforms"] = tform);  //  stores the term information
    List control_list = List::create(_["Iteration"] = iteration, _["Maximum Step"] = dbeta_max, _["Derivative Limiting"] = Lld_worst);  //  stores the total number of iterations used
    //
    NumericVector Lldd_vec(reqrdnum * reqrdnum);
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
    }
    Lldd_vec.attr("dim") = Dimension(reqrdnum, reqrdnum);
    const Map<MatrixXd> Lldd_mat(as<Map<MatrixXd> >(Lldd_vec));
    //
    MatrixXd cov;
    VectorXd stdev = VectorXd::Zero(totalnum);
    if (model_bool["observed_info"]) {
        cov = - 1 * Lldd_mat.inverse().matrix();  //  uses inverse information matrix to calculate the standard deviation
        for (int ij = 0; ij < totalnum; ij++) {
            if (KeepConstant[ij] == 0) {
                int pij_ind = ij - sum(head(KeepConstant, ij));
                stdev(ij) = sqrt(cov(pij_ind, pij_ind));
            }
        }
    } else {
        vector<double> InMa(pow(reqrdnum, 2), 0.0);
        Expected_Inform_Matrix_Poisson(nthreads, totalnum, PyrC, R, Rd, RdR, InMa, KeepConstant);
        NumericVector InMa_vec(reqrdnum * reqrdnum);  //  simplfied information matrix
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
            InMa_vec[ij * reqrdnum + jk] = InMa[ij*reqrdnum+jk];
            InMa_vec[jk * reqrdnum + ij] = InMa[ij * reqrdnum + jk];
        }
        InMa_vec.attr("dim") = Dimension(reqrdnum, reqrdnum);
        const Map<MatrixXd> InMa_mat(as<Map<MatrixXd> >(InMa_vec));    //
        cov = InMa_mat.inverse().matrix();  //  uses inverse information matrix to calculate the standard deviation
        for (int ij = 0; ij < totalnum; ij++) {
            if (KeepConstant[ij] == 0) {
                int pij_ind = ij - sum(head(KeepConstant, ij));
                stdev(ij) = sqrt(cov(pij_ind, pij_ind));
            }
        }
    }
    //
    res_list = List::create(_["LogLik"] = wrap(Ll[0]), _["First_Der"] = wrap(Lld), _["Second_Der"] = Lldd_vec, _["beta_0"] = wrap(beta_0), _["Standard_Deviation"] = wrap(stdev), _["Covariance"] = wrap(cov), _["AIC"] = 2*(totalnum-accumulate(KeepConstant.begin(), KeepConstant.end(), 0.0))+dev, _["BIC"] = (totalnum-accumulate(KeepConstant.begin(), KeepConstant.end(), 0.0))*log(mat_row)-2*Ll[0], _["Deviation"] = dev, _["Parameter_Lists"] = para_list, _["Control_List"] = control_list, _["Converged"] = convgd, _["Status"] = "PASSED");
    //  returns a list of results
    return res_list;
}

//' Primary Matched Case-Control starting point
//'
//' \code{LogLik_CaseCon_Omnibus} Primary Matched Case-Control starting point
//'
//' @inheritParams CPP_template
//'
//' @return only if it works
//' @noRd
//'
//
List LogLik_CaseCon_Omnibus(IntegerVector term_n, StringVector tform, NumericMatrix& a_ns, Ref<MatrixXd> df0, IntegerVector dfc, int fir, string modelform, double lr, List optim_para, NumericVector maxiters, int guesses, int halfmax, double epsilon, double step_max, double thres_step_max, double deriv_epsilon, const Ref<const MatrixXd>& df_m, NumericVector tu, int verbose, IntegerVector KeepConstant, int term_tot, string ties_method, int nthreads, NumericVector& Strata_vals, List model_bool, const double gmix_theta, const IntegerVector gmix_term, const Ref<const MatrixXd>& Lin_Sys, const Ref<const VectorXd>& Lin_Res) {
    //
    List temp_list = List::create(_["Status"] = "TEMP");  //  used as a dummy return value for code checking
    if (model_bool["constraint"]) {
        if (verbose >= 1) {
            Rcout << "linear constataints are currently not compatable with Case-Control model calculation" << endl;
        }
        temp_list = List::create(_["Status"] = "FAILED_WITH_BAD_MODEL_CONSTRAINT", _["LogLik"] = R_NaN);
        return temp_list;
    }
    //
    //  df0: covariate data
    //  ntime: number of event times for Cox PH
    //  totalnum: number of terms used
    //
    //  ------------------------------------------------------------------------- //  initialize
    const int mat_row = df0.rows();
    int ntime = tu.size();
    int totalnum;
    int reqrdnum;
    //  ------------------------------------------------------------------------- //  initialize
    if (!model_bool["null"]) {
        totalnum = term_n.size();
        reqrdnum = totalnum - sum(KeepConstant);
    } else {
        totalnum = 1;
        reqrdnum = 1;
    }
    //
    //  cout.precision: controls the number of significant digits printed
    //  nthreads: number of threads used for parallel operations
    //
    Rcout.precision(7);  //  forces higher precision numbers printed to terminal
    //
    //  Lld_worst: stores the highest magnitude log-likelihood derivative
    //
    //
    double Lld_worst = 0.0;  //  stores derivative value used to determine if every parameter is near convergence
    double dbeta_max = 0.0;  //  stores the largest step taken, determines if the step sizes are small enough for convergence
    //
    //  ---------------------------------------------
    //  To Start, needs to seperate the derivative terms
    //  ---------------------------------------------
    //
    //  ------------------------------------------------------------------------- //  initialize
    NumericVector a_n = a_ns.row(0);
    Map<VectorXd> beta_0(as<Map<VectorXd> >(a_n));
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
    if (!model_bool["null"]) {
        Cox_Refresh_R_TERM(totalnum, reqrdnum, term_tot, dint, dslp, thres_step_max, step_max, df0, T0, Td0, Tdd0, Te, R, Rd, Rdd, Dose, nonDose, TTerm, nonDose_LIN, nonDose_PLIN, nonDose_LOGLIN, RdR, RddR, model_bool);
    } else {
        R = MatrixXd::Constant(mat_row, 1, 1.0);
    }
    //  ------------------------------------------------------------------------- //  initialize
    int group_num = ntime*Strata_vals.size();
    IntegerMatrix RiskFail(group_num, 2);
    vector<vector<int> > RiskPairs(group_num);
    //
    vector<vector<double> > Recur_Base(group_num);
    vector<vector<vector<double> > > Recur_First(group_num, vector<vector<double>>(reqrdnum));
    vector<vector<vector<double> > > Recur_Second(group_num, vector<vector<double>>(reqrdnum*(reqrdnum + 1)/2));
    vector<double> strata_odds(group_num, 0.0);
    vector<double> strata_def(group_num, 0.0);
    vector<int> strata_cond(group_num, 1);
    //  ------------------------------------------------------------------------- //  initialize
    if (model_bool["time_risk"]) {
        if (model_bool["strata"]) {
            Make_Match_Time_Strata(model_bool, ntime, df_m, RiskFail, RiskPairs, Recur_Base, Recur_First, Recur_Second, strata_odds, strata_cond, nthreads, tu, Strata_vals);
        } else {
            Make_Match_Time(model_bool, ntime, df_m, RiskFail, RiskPairs, Recur_Base, Recur_First, Recur_Second, strata_odds, strata_cond, nthreads, tu);
        }
    } else {
        if (model_bool["strata"]) {
            Make_Match_Strata(model_bool, df_m, RiskFail, RiskPairs, Recur_Base, Recur_First, Recur_Second, strata_odds, strata_cond, nthreads, Strata_vals);
        } else {
            Make_Match(model_bool, df_m, RiskFail, RiskPairs, Recur_Base, Recur_First, Recur_Second, strata_odds, strata_cond, nthreads);
        }
    }
    for (int i = 0; i < group_num; i++) {
        strata_def[i] = strata_odds[i];
    }
    int reqrdcond = group_num - reduce(strata_cond.begin(), strata_cond.end());
    int totreqr = reqrdnum+reqrdcond;
    //  ------------------------------------------------------------------------- //  initialize
    vector<double> Ll(reqrdnum, 0.0);  //  log-likelihood values
    vector<double> Lld(reqrdnum, 0.0);  //  log-likelihood derivative values
    vector<double> Lldd(pow(reqrdnum, 2), 0.0);  //  the second derivative matrix has room for every combination, but only the lower triangle is calculated initially
    //
    vector<double> LldOdds(reqrdcond, 0.0);  //  storing the first derivative with respect to each strata odds ratio
    vector<double> LlddOdds(reqrdcond, 0.0);  //  storing the second derivative with respect to each strata odds ratio
    vector<double> LlddOddsBeta(reqrdcond*reqrdnum, 0.0);  //  storing the second derivative with respect to each strata odds ratio and model parameter
    double dev = 0.0;  //  Deviation needs to be calculated seperately, because the conditional and unconditional portions have different formula?
    if (model_bool["null"]) {
        Calculate_Recursive(model_bool, group_num, RiskFail, RiskPairs, totalnum, ntime, R, Rd, Rdd, Recur_Base, Recur_First, Recur_Second, nthreads, KeepConstant);
        Calc_Recur_LogLik(model_bool, group_num, RiskFail, RiskPairs, totalnum, ntime, R, Rd, Rdd, RdR, RddR, dev, Ll, Lld, Lldd, Recur_Base, Recur_First, Recur_Second, strata_odds, nthreads, KeepConstant, strata_cond, LldOdds, LlddOdds, LlddOddsBeta);
        //
        List res_list = List::create(_["LogLik"] = wrap(Ll[0]), _["Deviation"] = wrap(dev), _["Status"] = "PASSED", _["StrataOdds"] = wrap(strata_odds), _["FreeParameters"] = wrap(reqrdnum), _["FreeSets"] = wrap(reqrdcond));
        //  returns a list of results
        return res_list;
    }
    //  ------------------------------------------------------------------------- //  initialize
    //  the log-likelihood is calculated in parallel over the risk groups
    vector <double> Ll_comp(2, Ll[0]);  //  vector to compare values
    double step_max0 = step_max;
    double thres_step_max0 = thres_step_max;
    //
    vector<double> dbeta(totalnum, 0.0);
    vector<double> dstrata(group_num, 0.0);
    NumericVector m_g_store(totreqr);
    NumericVector v_beta_store(totreqr);
    //
    //  --------------------------
    //  always starts from initial guess
    //  --------------------------
    vector<double> beta_c(totalnum, 0.0);
    vector<double> beta_a(totalnum, 0.0);
    vector<double> beta_best(totalnum, 0.0);
    vector<double> beta_p(totalnum, 0.0);
    VectorXd::Map(&beta_p[0], beta_0.size()) = beta_0;  //  stores previous parameters
    VectorXd::Map(&beta_c[0], beta_0.size()) = beta_0;  //  stores current parameters
    VectorXd::Map(&beta_a[0], beta_0.size()) = beta_0;  //  stores a refrence value for parameters
    VectorXd::Map(&beta_best[0], beta_0.size()) = beta_0;  //  stores the best parameters
    //
    vector<double> strata_c(group_num, 0.0);
    vector<double> strata_a(group_num, 0.0);
    vector<double> strata_best(group_num, 0.0);
    vector<double> strata_p(group_num, 0.0);
    for (int i = 0; i < group_num; i++) {
        strata_p[i] = strata_odds[i];
        strata_c[i] = strata_odds[i];
        strata_a[i] = strata_odds[i];
        strata_best[i] = strata_odds[i];
    }
    //
    double halves = 0;  //  number of half-steps taken
    int ind0 = fir;  //  used for validations
    int iteration = 0;  //  iteration number
    int maxiter = 0;
    //
    bool convgd = FALSE;
    int iter_stop  = 0;  //  tracks if the iterations should be stopped for convergence
    //
    NumericMatrix strata_fin(a_ns.rows(), group_num);
    NumericMatrix beta_fin(a_ns.rows(), totalnum);
    NumericVector LL_fin(a_ns.rows());
    //
    double Ll_abs_best = 10;
    vector<double> beta_abs_best(totalnum, 0.0);
    vector<double> strata_abs_best(group_num, 0.0);
    int guess_abs_best = - 1;
    double Ll_iter_best = 10;
    ///
    //  Variables that are used for the risk check function shared across cox, poisson, and log bound functions
    MatrixXd dev_temp = MatrixXd::Zero(1, 1);
    MatrixXd PyrC = MatrixXd::Zero(1, 1);
    //
    //
    for (int guess = 0; guess <guesses; guess++) {
        Cox_Refresh_R_TERM(totalnum, reqrdnum, term_tot, dint, dslp, thres_step_max, step_max, df0, T0, Td0, Tdd0, Te, R, Rd, Rdd, Dose, nonDose, TTerm, nonDose_LIN, nonDose_PLIN, nonDose_LOGLIN, RdR, RddR, model_bool);
        fill(Ll.begin(), Ll.end(), 0.0);
        if (!model_bool["single"]) {
            fill(Lld.begin(), Lld.end(), 0.0);
            fill(Lldd.begin(), Lldd.end(), 0.0);
            fill(LldOdds.begin(), LldOdds.end(), 0.0);
            fill(LlddOdds.begin(), LlddOdds.end(), 0.0);
            fill(LlddOddsBeta.begin(), LlddOddsBeta.end(), 0.0);
        }
        dev = 0.0;
        if (model_bool["gradient"]) {
            m_g_store.fill(0);
            v_beta_store.fill(0);
        }
        beta_p = beta_best;  //
        beta_a = beta_best;  //
        beta_c = beta_best;  //
        strata_p = strata_best;  //
        strata_a = strata_best;  //
        strata_c = strata_best;  //
        //
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
        convgd = FALSE;
        iter_stop  = 0;  //  tracks if the iterations should be stopped for convergence
        //
        maxiter = maxiters[guess];
        a_n = a_ns.row(guess);
        for (int i = 0; i < beta_0.size(); i++) {
            beta_0[i] = a_n[i];
        }
        for (int i = 0; i < group_num; i++) {
            strata_odds[i] = strata_def[i];
        }
        if (verbose >= 4) {
            Rcout << "C++ Note: starting guess " << guess << endl;
        }
        Cox_Term_Risk_Calc(modelform, tform, term_n, totalnum, fir, dfc, term_tot, T0, Td0, Tdd0, Te, R, Rd, Rdd, Dose, nonDose, beta_0, df0, thres_step_max, step_max, TTerm, nonDose_LIN, nonDose_PLIN, nonDose_LOGLIN, RdR, RddR, nthreads, KeepConstant, verbose, model_bool, gmix_theta, gmix_term);
        if ((R.minCoeff() <= 0) || (R.hasNaN())) {
            if (verbose >= 1) {
                Rcout << "C++ Error: A non-positive risk was detected: " << R.minCoeff() << endl;
                Rcout << "C++ Warning: final failing values ";
                for (int ijk = 0; ijk < totalnum; ijk++) {
                    Rcout << beta_0[ijk] << " ";
                }
                Rcout << " " << endl;
            }
            //
            temp_list = List::create(_["beta_0"] = wrap(beta_0), _["StrataOdds"] = wrap(strata_odds), _["Deviation"] = R_NaN, _["Status"] = "FAILED_WITH_NEGATIVE_RISK", _["LogLik"] = R_NaN);
            return temp_list;
        }
        Calculate_Recursive(model_bool, group_num, RiskFail, RiskPairs, totalnum, ntime, R, Rd, Rdd, Recur_Base, Recur_First, Recur_Second , nthreads, KeepConstant);
        Calc_Recur_LogLik(model_bool, group_num, RiskFail, RiskPairs, totalnum, ntime, R, Rd, Rdd, RdR, RddR, dev, Ll, Lld, Lldd, Recur_Base, Recur_First, Recur_Second, strata_odds, nthreads, KeepConstant, strata_cond, LldOdds, LlddOdds, LlddOddsBeta);
        Print_LL(reqrdnum, totalnum, beta_0, Ll, Lld, Lldd, verbose, model_bool);
        Print_LL_Background(reqrdnum, totalnum, group_num, reqrdcond, strata_odds, LldOdds, LlddOdds, LlddOddsBeta, verbose, model_bool);
        //  NOW WE RUN THE ITERATIONS
        for (int i = 0; i < beta_0.size(); i++) {
            beta_c[i] = beta_0[i];
        }
        for (vector<double>::size_type i = 0; i < strata_odds.size(); i++) {
            strata_c[i] = strata_odds[i];
        }
        while ((iteration < maxiter) && (iter_stop == 0)) {
            iteration++;
            beta_p = beta_c;  //
            beta_a = beta_c;  //
            beta_best = beta_c;  //
            strata_p = strata_c;  //
            strata_a = strata_c;  //
            strata_best = strata_c;  //
            if (model_bool["gradient"]) {
                Calc_Change_Background_Gradient(nthreads, model_bool, totalnum, group_num, optim_para, iteration, step_max, Lld, m_g_store, v_beta_store, dbeta, KeepConstant, strata_cond, LldOdds, dstrata);
            } else {
                Calc_Change_Background(nthreads, totalnum, group_num, thres_step_max, lr, step_max, Ll, Lld, Lldd, dbeta, tform, thres_step_max, step_max, KeepConstant, strata_cond, LldOdds, LlddOdds, LlddOddsBeta, dstrata);
                Intercept_Bound(nthreads, totalnum, beta_0, dbeta, dfc, df0, KeepConstant, tform);
            }
            if ((Ll_iter_best > 0) || (Ll_iter_best < Ll[ind0])) {
                Ll_iter_best = Ll[ind0];
            }
            //
            if (model_bool["gradient"]) {
                //
                halves = 0;
                while (halves < halfmax) {
                    for (int ij = 0; ij < totalnum; ij++) {
                        beta_0[ij] = beta_a[ij] + dbeta[ij];
                        beta_c[ij] = beta_0[ij];
                    }
                    for (int ij = 0; ij < group_num; ij++) {
                        strata_odds[ij] = strata_a[ij] + dstrata[ij];
                        strata_c[ij] = strata_odds[ij];
                    }
                    Cox_Term_Risk_Calc(modelform, tform, term_n, totalnum, fir, dfc, term_tot, T0, Td0, Tdd0, Te, R, Rd, Rdd, Dose, nonDose, beta_0, df0, thres_step_max, step_max, TTerm, nonDose_LIN, nonDose_PLIN, nonDose_LOGLIN, RdR, RddR, nthreads, KeepConstant, verbose, model_bool, gmix_theta, gmix_term);
                    //
                    if ((R.minCoeff() <= 0) || (R.hasNaN())) {
                        halves++;
                        for (int ijk = 0; ijk < totalnum; ijk++) {
                            int tij = term_n[ijk];
                            if (TTerm.col(tij).minCoeff() <= 0) {
                                dbeta[ijk] = dbeta[ijk] / 2.0;
                            } else if (isinf(TTerm.col(tij).maxCoeff())) {
                                dbeta[ijk] = dbeta[ijk] / 2.0;
                            } else if (isnan(TTerm.col(tij).minCoeff())) {
                                dbeta[ijk] = dbeta[ijk] / 2.0;
                            }
                        }
                    } else {
                        halves = halfmax + 1;
                        Calculate_Recursive(model_bool, group_num, RiskFail, RiskPairs, totalnum, ntime, R, Rd, Rdd, Recur_Base, Recur_First, Recur_Second , nthreads, KeepConstant);
                        Calc_Recur_LogLik(model_bool, group_num, RiskFail, RiskPairs, totalnum, ntime, R, Rd, Rdd, RdR, RddR, dev, Ll, Lld, Lldd, Recur_Base, Recur_First, Recur_Second, strata_odds, nthreads, KeepConstant, strata_cond, LldOdds, LlddOdds, LlddOddsBeta);
                        Print_LL(reqrdnum, totalnum, beta_0, Ll, Lld, Lldd, verbose, model_bool);
                        Print_LL_Background(reqrdnum, totalnum, group_num, reqrdcond, strata_odds, LldOdds, LlddOdds, LlddOddsBeta, verbose, model_bool);
                        if (Ll[ind0] > Ll_abs_best) {  //  if a better point wasn't found, takes a half-step
                            for (int ijk = 0; ijk < totalnum; ijk++) {
                                beta_best[ijk] = beta_c[ijk];
                            }
                            for (int ijk = 0; ijk < group_num; ijk++) {
                                strata_best[ijk] = strata_c[ijk];  //
                            }
                        }
                        for (int ijk = 0; ijk < totalnum; ijk++) {  //  totalnum*(totalnum + 1)/2
                            beta_0[ijk] = beta_c[ijk];
                        }
                        for (int ijk = 0; ijk < group_num; ijk++) {
                            strata_odds[ijk] = strata_c[ijk];  //
                        }
                    }
                }
            } else {
                halves = 0;
                while ((Ll[ind0] <= Ll_iter_best) && (halves < halfmax)) {  //  repeats until half-steps maxed or an improvement
                    for (int ij = 0; ij < totalnum; ij++) {
                        beta_0[ij] = beta_a[ij] + dbeta[ij];
                        beta_c[ij] = beta_0[ij];
                    }
                    for (int ij = 0; ij < group_num; ij++) {
                        strata_odds[ij] = strata_a[ij] + dstrata[ij];
                        strata_c[ij] = strata_odds[ij];
                    }
                    //  -----------------------------------------------------------------------------------------------------------------------------------------------//
                    //  The same subterm, risk, sides, and log-likelihood calculations are performed every half-step and iteration
                    //  -----------------------------------------------------------------------------------------------------------------------------------------------//
                    Cox_Term_Risk_Calc(modelform, tform, term_n, totalnum, fir, dfc, term_tot, T0, Td0, Tdd0, Te, R, Rd, Rdd, Dose, nonDose, beta_0, df0, thres_step_max, step_max, TTerm, nonDose_LIN, nonDose_PLIN, nonDose_LOGLIN, RdR, RddR, nthreads, KeepConstant, verbose, model_bool, gmix_theta, gmix_term);
                    //
                    if ((R.minCoeff() <= 0) || (R.hasNaN())) {
//                        #ifdef _OPENMP
//                        #pragma omp parallel for schedule(dynamic) num_threads(nthreads)
//                        #endif
                        for (int ijk = 0; ijk < totalnum; ijk++) {
                            int tij = term_n[ijk];
                            if (TTerm.col(tij).minCoeff() <= 0) {
                                dbeta[ijk] = dbeta[ijk] / 2.0;
                            } else if (isinf(TTerm.col(tij).maxCoeff())) {
                                dbeta[ijk] = dbeta[ijk] / 2.0;
                            } else if (isnan(TTerm.col(tij).minCoeff())) {
                                dbeta[ijk] = dbeta[ijk] / 2.0;
                            }
                        }
                        halves+=0.2;
                    } else {
                        halves++;
                        Calculate_Recursive(model_bool, group_num, RiskFail, RiskPairs, totalnum, ntime, R, Rd, Rdd, Recur_Base, Recur_First, Recur_Second , nthreads, KeepConstant);
                        Calc_Recur_LogLik(model_bool, group_num, RiskFail, RiskPairs, totalnum, ntime, R, Rd, Rdd, RdR, RddR, dev, Ll, Lld, Lldd, Recur_Base, Recur_First, Recur_Second, strata_odds, nthreads, KeepConstant, strata_cond, LldOdds, LlddOdds, LlddOddsBeta);
                        Print_LL(reqrdnum, totalnum, beta_0, Ll, Lld, Lldd, verbose, model_bool);
                        Print_LL_Background(reqrdnum, totalnum, group_num, reqrdcond, strata_odds, LldOdds, LlddOdds, LlddOddsBeta, verbose, model_bool);
                        //
                        if (Ll[ind0] <= Ll_abs_best) {  //  if a better point wasn't found, takes a half-step
                            for (int ijk = 0; ijk < totalnum; ijk++) {
                                dbeta[ijk] = dbeta[ijk] * 0.5;  //
                            }
                            for (int ijk = 0; ijk < group_num; ijk++) {
                                dstrata[ijk] = dstrata[ijk] * 0.5;  //
                            }
                        } else {  //  if improved, updates the best vector
                            for (int ijk = 0; ijk < totalnum; ijk++) {
                                beta_best[ijk] = beta_c[ijk];
                            }
                            for (int ijk = 0; ijk < group_num; ijk++) {
                                strata_best[ijk] = strata_odds[ijk];  //
                            }
                        }
                        for (int ijk = 0; ijk < totalnum; ijk++) {  //  totalnum*(totalnum + 1)/2
                            beta_0[ijk] = beta_c[ijk];
                        }
                        for (int ijk = 0; ijk < group_num; ijk++) {
                            strata_odds[ijk] = strata_c[ijk];  //
                        }
                    }
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
                    strata_p = strata_best;  //
                    strata_a = strata_best;  //
                    strata_c = strata_best;  //
                    for (int ij = 0; ij < group_num; ij++) {
                        strata_odds[ij] = strata_best[ij];
                    }
                    Cox_Term_Risk_Calc(modelform, tform, term_n, totalnum, fir, dfc, term_tot, T0, Td0, Tdd0, Te, R, Rd, Rdd, Dose, nonDose, beta_0, df0, thres_step_max, step_max, TTerm, nonDose_LIN, nonDose_PLIN, nonDose_LOGLIN, RdR, RddR, nthreads, KeepConstant, verbose, model_bool, gmix_theta, gmix_term);
                }
            }
            for (int ij = 0; ij < reqrdnum; ij++) {
                if (abs(Lld[ij]) > Lld_worst) {
                    Lld_worst = abs(Lld[ij]);
                }
            }
            dbeta_max = abs(dbeta[0]);
            for (int ij = 1; ij < totalnum; ij++) {
                if (abs(dbeta[ij]) > dbeta_max) {
                    dbeta_max = abs(dbeta[ij]);
                }
            }
            for (int ij = 0; ij < group_num; ij++) {
                if (abs(dstrata[ij]) > dbeta_max) {
                    dbeta_max = abs(dstrata[ij]);
                }
            }
            if (Lld_worst < deriv_epsilon) {  //  ends if the derivatives are low enough
                iter_stop = 1;
            }
            if (step_max < epsilon) {  //  if the maximum change is too low, then it ends
                iter_stop = 1;
            }
            if (dbeta_max < epsilon) {  //  if the maximum change is too low, then it ends
                iter_stop = 1;
            }
        }
        if (model_bool["gradient"]) {
            beta_p = beta_best;  //
            beta_a = beta_best;  //
            beta_c = beta_best;  //
            for (int ij = 0; ij < totalnum; ij++) {
                beta_0[ij] = beta_best[ij];
            }
            strata_p = strata_best;  //
            strata_a = strata_best;  //
            strata_c = strata_best;  //
            for (int ij = 0; ij < group_num; ij++) {
                strata_odds[ij] = strata_best[ij];
            }
            Cox_Term_Risk_Calc(modelform, tform, term_n, totalnum, fir, dfc, term_tot, T0, Td0, Tdd0, Te, R, Rd, Rdd, Dose, nonDose, beta_0, df0, thres_step_max, step_max, TTerm, nonDose_LIN, nonDose_PLIN, nonDose_LOGLIN, RdR, RddR, nthreads, KeepConstant, verbose, model_bool, gmix_theta, gmix_term);
        }
        //  -----------------------------------------------
        //  Performing Full Calculation to get full second derivative matrix
        //  -----------------------------------------------
        fill(Ll.begin(), Ll.end(), 0.0);
        if (!model_bool["single"]) {
            fill(Lld.begin(), Lld.end(), 0.0);
            fill(Lldd.begin(), Lldd.end(), 0.0);
            fill(LldOdds.begin(), LldOdds.end(), 0.0);
            fill(LlddOdds.begin(), LlddOdds.end(), 0.0);
            fill(LlddOddsBeta.begin(), LlddOddsBeta.end(), 0.0);
        }
        dev = 0.0;
        Calculate_Recursive(model_bool, group_num, RiskFail, RiskPairs, totalnum, ntime, R, Rd, Rdd, Recur_Base, Recur_First, Recur_Second , nthreads, KeepConstant);
        Calc_Recur_LogLik(model_bool, group_num, RiskFail, RiskPairs, totalnum, ntime, R, Rd, Rdd, RdR, RddR, dev, Ll, Lld, Lldd, Recur_Base, Recur_First, Recur_Second, strata_odds, nthreads, KeepConstant, strata_cond, LldOdds, LlddOdds, LlddOddsBeta);
        Print_LL(reqrdnum, totalnum, beta_0, Ll, Lld, Lldd, verbose, model_bool);
        Print_LL_Background(reqrdnum, totalnum, group_num, reqrdcond, strata_odds, LldOdds, LlddOdds, LlddOddsBeta, verbose, model_bool);
        a_n = beta_0;
        beta_fin(guess, _) = a_n;
        for (int i = 0; i < group_num; i++) {
            strata_fin(guess, i) = strata_odds[i];
        }
        LL_fin[guess] = Ll[0];
        if ((Ll_abs_best > 0) || (Ll_abs_best < Ll[ind0])) {
            Ll_abs_best = Ll[ind0];
            beta_abs_best = beta_c;
            strata_abs_best = strata_c;
            guess_abs_best = guess;
        }
        //
    }
    //  ----------------------------------------------------------------------------------- //
    //             NOW WE NEED TO FIND THE BEST, ASSIGN, AND CONTINUE REGRESSION
    //  ----------------------------------------------------------------------------------- //
    fill(Ll.begin(), Ll.end(), 0.0);
    if (!model_bool["single"]) {
        fill(Lld.begin(), Lld.end(), 0.0);
        fill(Lldd.begin(), Lldd.end(), 0.0);
        fill(LldOdds.begin(), LldOdds.end(), 0.0);
        fill(LlddOdds.begin(), LlddOdds.end(), 0.0);
        fill(LlddOddsBeta.begin(), LlddOddsBeta.end(), 0.0);
    }
    dev = 0.0;
    beta_p = beta_best;  //
    beta_a = beta_best;  //
    beta_c = beta_best;  //
    strata_p = strata_best;  //
    strata_a = strata_best;  //
    strata_c = strata_best;  //
    step_max = step_max0;
    thres_step_max = thres_step_max0;
    iter_stop = 0;
    halves = 0;
    iteration = 0;
    halves = 0;  //  number of half-steps taken
    ind0 = fir;  //  used for validations
    iteration = 0;  //  iteration number
    //
    convgd = FALSE;
    iter_stop  = 0;  //  tracks if the iterations should be stopped for convergence
    //
    int guess_max = guess_abs_best;
    if (verbose >= 3) {
        Rcout << "C++ Note: Guess Results" << endl;
        Rcout << "Guess number, parameter values, Log-Likelihood" << endl;
        NumericVector beta_temp;
        for (int i = 0;  i < guesses; i++) {
            beta_temp = wrap(beta_fin.row(i));
            if (i == guess_max) {
                Rcout << i << ", " << beta_temp << ", " << LL_fin[i] << "<-- Best Guess" << endl;
            } else {
                Rcout << i << ", " << beta_temp << ", " << LL_fin[i] << endl;
            }
        }
    }
    //
    maxiter = maxiters[guesses];
    a_n = beta_abs_best;
    for (int i = 0; i < beta_0.size(); i++) {
        beta_0[i] = a_n[i];
        beta_c[i] = beta_0[i];
    }
    for (int i = 0; i < group_num; i++) {
        strata_odds[i] = strata_best[i];
        strata_c[i] = strata_odds[i];
    }
    if (model_bool["gradient"]) {
        m_g_store.fill(0);
        v_beta_store.fill(0);
    }
    //
    //  Calculates the subterm and term values
    Cox_Term_Risk_Calc(modelform, tform, term_n, totalnum, fir, dfc, term_tot, T0, Td0, Tdd0, Te, R, Rd, Rdd, Dose, nonDose, beta_0, df0, thres_step_max, step_max, TTerm, nonDose_LIN, nonDose_PLIN, nonDose_LOGLIN, RdR, RddR, nthreads, KeepConstant, verbose, model_bool, gmix_theta, gmix_term);
    //
    //  -------------------------------------------------------------------------------------------
    Calculate_Recursive(model_bool, group_num, RiskFail, RiskPairs, totalnum, ntime, R, Rd, Rdd, Recur_Base, Recur_First, Recur_Second , nthreads, KeepConstant);
    Calc_Recur_LogLik(model_bool, group_num, RiskFail, RiskPairs, totalnum, ntime, R, Rd, Rdd, RdR, RddR, dev, Ll, Lld, Lldd, Recur_Base, Recur_First, Recur_Second, strata_odds, nthreads, KeepConstant, strata_cond, LldOdds, LlddOdds, LlddOddsBeta);
    Print_LL(reqrdnum, totalnum, beta_0, Ll, Lld, Lldd, verbose, model_bool);
    Print_LL_Background(reqrdnum, totalnum, group_num, reqrdcond, strata_odds, LldOdds, LlddOdds, LlddOddsBeta, verbose, model_bool);
    List res_list;
    //
    if (model_bool["single"]) {
        res_list = List::create(_["LogLik"] = wrap(Ll[0]), _["Deviation"] = wrap(dev), _["beta_0"] = wrap(beta_0), _["StrataOdds"] = wrap(strata_odds), _["FreeParameters"] = wrap(reqrdnum), _["FreeSets"] = wrap(reqrdcond), _["Status"] = "PASSED");
        //  returns a list of results
        return res_list;
    }
    //
    while ((iteration < maxiter) && (iter_stop == 0)) {
        iteration++;
        Print_LL(reqrdnum, totalnum, beta_0, Ll, Lld, Lldd, verbose, model_bool);
        Print_LL_Background(reqrdnum, totalnum, group_num, reqrdcond, strata_odds, LldOdds, LlddOdds, LlddOddsBeta, verbose, model_bool);
        beta_p = beta_c;  //
        beta_a = beta_c;  //
        beta_best = beta_c;  //
        strata_p = strata_c;  //
        strata_a = strata_c;  //
        strata_best = strata_c;  //
        if (model_bool["gradient"]) {
            Calc_Change_Background_Gradient(nthreads, model_bool, totalnum, group_num, optim_para, iteration, step_max, Lld, m_g_store, v_beta_store, dbeta, KeepConstant, strata_cond, LldOdds, dstrata);
        } else {
            Calc_Change_Background(nthreads, totalnum, group_num, thres_step_max, lr, step_max, Ll, Lld, Lldd, dbeta, tform, thres_step_max, step_max, KeepConstant, strata_cond, LldOdds, LlddOdds, LlddOddsBeta, dstrata);
            Intercept_Bound(nthreads, totalnum, beta_0, dbeta, dfc, df0, KeepConstant, tform);
        }
        if ((Ll_iter_best > 0) || (Ll_iter_best < Ll[ind0])) {
            Ll_iter_best = Ll[ind0];
        }
        //
        if (model_bool["gradient"]) {
            //
            halves = 0;
            while (halves < halfmax) {
                for (int ij = 0; ij < totalnum; ij++) {
                    beta_0[ij] = beta_a[ij] + dbeta[ij];
                    beta_c[ij] = beta_0[ij];
                }
                for (int ij = 0; ij < group_num; ij++) {
                    strata_odds[ij] = strata_a[ij] + dstrata[ij];
                    strata_c[ij] = strata_odds[ij];
                }
                Cox_Term_Risk_Calc(modelform, tform, term_n, totalnum, fir, dfc, term_tot, T0, Td0, Tdd0, Te, R, Rd, Rdd, Dose, nonDose, beta_0, df0, thres_step_max, step_max, TTerm, nonDose_LIN, nonDose_PLIN, nonDose_LOGLIN, RdR, RddR, nthreads, KeepConstant, verbose, model_bool, gmix_theta, gmix_term);
                //
                if ((R.minCoeff() <= 0) || (R.hasNaN())) {
                    halves++;
                    for (int ijk = 0; ijk < totalnum; ijk++) {
                        int tij = term_n[ijk];
                        if (TTerm.col(tij).minCoeff() <= 0) {
                            dbeta[ijk] = dbeta[ijk] / 2.0;
                        } else if (isinf(TTerm.col(tij).maxCoeff())) {
                            dbeta[ijk] = dbeta[ijk] / 2.0;
                        } else if (isnan(TTerm.col(tij).minCoeff())) {
                            dbeta[ijk] = dbeta[ijk] / 2.0;
                        }
                    }
                } else {
                    halves = halfmax + 1;
                    Calculate_Recursive(model_bool, group_num, RiskFail, RiskPairs, totalnum, ntime, R, Rd, Rdd, Recur_Base, Recur_First, Recur_Second , nthreads, KeepConstant);
                    Calc_Recur_LogLik(model_bool, group_num, RiskFail, RiskPairs, totalnum, ntime, R, Rd, Rdd, RdR, RddR, dev, Ll, Lld, Lldd, Recur_Base, Recur_First, Recur_Second, strata_odds, nthreads, KeepConstant, strata_cond, LldOdds, LlddOdds, LlddOddsBeta);
                    Print_LL(reqrdnum, totalnum, beta_0, Ll, Lld, Lldd, verbose, model_bool);
                    Print_LL_Background(reqrdnum, totalnum, group_num, reqrdcond, strata_odds, LldOdds, LlddOdds, LlddOddsBeta, verbose, model_bool);
                    if (Ll[ind0] > Ll_abs_best) {  //  if a better point wasn't found, takes a half-step
                        for (int ijk = 0; ijk < totalnum; ijk++) {
                            beta_best[ijk] = beta_c[ijk];
                        }
                        for (int ijk = 0; ijk < group_num; ijk++) {
                            strata_best[ijk] = strata_c[ijk];  //
                        }
                    }
                    for (int ijk = 0; ijk < totalnum; ijk++) {  //  totalnum*(totalnum + 1)/2
                        beta_0[ijk] = beta_c[ijk];
                    }
                    for (int ijk = 0; ijk < group_num; ijk++) {
                        strata_odds[ijk] = strata_c[ijk];  //
                    }
                }
            }
        } else {
            halves = 0;
            while ((Ll[ind0] <= Ll_iter_best) && (halves < halfmax)) {  //  repeats until half-steps maxed or an improvement
                for (int ij = 0; ij < totalnum; ij++) {
                    beta_0[ij] = beta_a[ij] + dbeta[ij];
                    beta_c[ij] = beta_0[ij];
                }
                for (int ij = 0; ij < group_num; ij++) {
                    strata_odds[ij] = strata_a[ij] + dstrata[ij];
                    strata_c[ij] = strata_odds[ij];
                }
                //  ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------//
                //  The same subterm, risk, sides, and log-likelihood calculations are performed every half-step and iteration
                //  ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------//
                Cox_Term_Risk_Calc(modelform, tform, term_n, totalnum, fir, dfc, term_tot, T0, Td0, Tdd0, Te, R, Rd, Rdd, Dose, nonDose, beta_0, df0, thres_step_max, step_max, TTerm, nonDose_LIN, nonDose_PLIN, nonDose_LOGLIN, RdR, RddR, nthreads, KeepConstant, verbose, model_bool, gmix_theta, gmix_term);
                //
                if ((R.minCoeff() <= 0) || (R.hasNaN())) {
                    for (int ijk = 0; ijk < totalnum; ijk++) {
                        int tij = term_n[ijk];
                        if (TTerm.col(tij).minCoeff() <= 0) {
                            dbeta[ijk] = dbeta[ijk] / 2.0;
                        } else if (isinf(TTerm.col(tij).maxCoeff())) {
                            dbeta[ijk] = dbeta[ijk] / 2.0;
                        } else if (isnan(TTerm.col(tij).minCoeff())) {
                            dbeta[ijk] = dbeta[ijk] / 2.0;
                        }
                    }
                    halves+=0.2;
                } else {
                    halves++;
                    Calculate_Recursive(model_bool, group_num, RiskFail, RiskPairs, totalnum, ntime, R, Rd, Rdd, Recur_Base, Recur_First, Recur_Second , nthreads, KeepConstant);
                    Calc_Recur_LogLik(model_bool, group_num, RiskFail, RiskPairs, totalnum, ntime, R, Rd, Rdd, RdR, RddR, dev, Ll, Lld, Lldd, Recur_Base, Recur_First, Recur_Second, strata_odds, nthreads, KeepConstant, strata_cond, LldOdds, LlddOdds, LlddOddsBeta);
                    Print_LL(reqrdnum, totalnum, beta_0, Ll, Lld, Lldd, verbose, model_bool);
                    Print_LL_Background(reqrdnum, totalnum, group_num, reqrdcond, strata_odds, LldOdds, LlddOdds, LlddOddsBeta, verbose, model_bool);
                    //
                    if (Ll[ind0] <= Ll_iter_best) {  //  if a better point wasn't found, takes a half-step
                        for (int ijk = 0; ijk < totalnum; ijk++) {
                            dbeta[ijk] = dbeta[ijk] * 0.5;  //
                        }
                        for (int ijk = 0; ijk < group_num; ijk++) {
                            dstrata[ijk] = dstrata[ijk] * 0.5;  //
                        }
                    } else {  //  if improved, updates the best vector
                        for (int ijk = 0; ijk < totalnum; ijk++) {
                            beta_best[ijk] = beta_c[ijk];
                        }
                        for (int ijk = 0; ijk < group_num; ijk++) {
                            strata_best[ijk] = strata_odds[ijk];  //
                        }
                    }
                    for (int ijk = 0; ijk < totalnum; ijk++) {  //  totalnum*(totalnum + 1)/2
                        beta_0[ijk] = beta_c[ijk];
                    }
                    for (int ijk = 0; ijk < group_num; ijk++) {
                        strata_odds[ijk] = strata_c[ijk];  //
                    }
                }
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
                strata_p = strata_best;  //
                strata_a = strata_best;  //
                strata_c = strata_best;  //
                for (int ij = 0; ij < group_num; ij++) {
                    strata_odds[ij] = strata_best[ij];
                }
                Cox_Term_Risk_Calc(modelform, tform, term_n, totalnum, fir, dfc, term_tot, T0, Td0, Tdd0, Te, R, Rd, Rdd, Dose, nonDose, beta_0, df0, thres_step_max, step_max, TTerm, nonDose_LIN, nonDose_PLIN, nonDose_LOGLIN, RdR, RddR, nthreads, KeepConstant, verbose, model_bool, gmix_theta, gmix_term);
            }
        }
        for (int ij = 0; ij < reqrdnum; ij++) {
            if (abs(Lld[ij]) > Lld_worst) {
                Lld_worst = abs(Lld[ij]);
            }
        }
        dbeta_max = abs(dbeta[0]);
        for (int ij = 1; ij < totalnum; ij++) {
            if (abs(dbeta[ij]) > dbeta_max) {
                dbeta_max = abs(dbeta[ij]);
            }
        }
        for (int ij = 0; ij < group_num; ij++) {
            if (abs(dstrata[ij]) > dbeta_max) {
                dbeta_max = abs(dstrata[ij]);
            }
        }
        if (Lld_worst < deriv_epsilon) {  //  ends if the derivatives are low enough
            iter_stop = 1;
        }
        if (step_max < epsilon) {  //  if the maximum change is too low, then it ends
            iter_stop = 1;
        }
        if (dbeta_max < epsilon) {  //  if the maximum change is too low, then it ends
            iter_stop = 1;
        }
    }
    //  ----------------------------------------------------------------------------------- //
    //               NOW WE WRAP UP
    //  ----------------------------------------------------------------------------------- //
    if (model_bool["gradient"]) {
        beta_p = beta_best;  //
        beta_a = beta_best;  //
        beta_c = beta_best;  //
        for (int ij = 0; ij < totalnum; ij++) {
            beta_0[ij] = beta_best[ij];
        }
        strata_p = strata_best;  //
        strata_a = strata_best;  //
        strata_c = strata_best;  //
        for (int ij = 0; ij < group_num; ij++) {
            strata_odds[ij] = strata_best[ij];
        }
    }
    fill(Ll.begin(), Ll.end(), 0.0);
    if (!model_bool["single"]) {
        fill(Lld.begin(), Lld.end(), 0.0);
        fill(Lldd.begin(), Lldd.end(), 0.0);
        fill(LldOdds.begin(), LldOdds.end(), 0.0);
        fill(LlddOdds.begin(), LlddOdds.end(), 0.0);
        fill(LlddOddsBeta.begin(), LlddOddsBeta.end(), 0.0);
    }
    dev = 0.0;
    bool true_gradient = model_bool["gradient"];
    model_bool["gradient"] = false;
    Cox_Term_Risk_Calc(modelform, tform, term_n, totalnum, fir, dfc, term_tot, T0, Td0, Tdd0, Te, R, Rd, Rdd, Dose, nonDose, beta_0, df0, thres_step_max, step_max, TTerm, nonDose_LIN, nonDose_PLIN, nonDose_LOGLIN, RdR, RddR, nthreads, KeepConstant, verbose, model_bool, gmix_theta, gmix_term);
    Calculate_Recursive(model_bool, group_num, RiskFail, RiskPairs, totalnum, ntime, R, Rd, Rdd, Recur_Base, Recur_First, Recur_Second , nthreads, KeepConstant);
    Calc_Recur_LogLik(model_bool, group_num, RiskFail, RiskPairs, totalnum, ntime, R, Rd, Rdd, RdR, RddR, dev, Ll, Lld, Lldd, Recur_Base, Recur_First, Recur_Second, strata_odds, nthreads, KeepConstant, strata_cond, LldOdds, LlddOdds, LlddOddsBeta);
    model_bool["gradient"] = true_gradient;
    Print_LL(reqrdnum, totalnum, beta_0, Ll, Lld, Lldd, verbose, model_bool);
    Print_LL_Background(reqrdnum, totalnum, group_num, reqrdcond, strata_odds, LldOdds, LlddOdds, LlddOddsBeta, verbose, model_bool);
    //
    List para_list;
    if (!model_bool["basic"]) {
        para_list = List::create(_["term_n"] = term_n, _["tforms"] = tform);  //  stores the term information
    }
    List control_list = List::create(_["Iteration"] = iteration, _["Maximum Step"] = dbeta_max, _["Derivative Limiting"] = Lld_worst);  //  stores the total number of iterations used
    //
    //
    NumericVector Lldd_vec(reqrdnum * reqrdnum);  //  simplfied information matrix
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
    }
    Lldd_vec.attr("dim") = Dimension(reqrdnum, reqrdnum);
    const Map<MatrixXd> Lldd_mat(as<Map<MatrixXd> >(Lldd_vec));
    //
    MatrixXd cov;
    VectorXd stdev = VectorXd::Zero(totalnum);
    //
    if (model_bool["observed_info"]) {
        cov = - 1 * Lldd_mat.inverse().matrix();  //  uses inverse information matrix to calculate the standard deviation
        for (int ij = 0; ij < totalnum; ij++) {
            if (KeepConstant[ij] == 0) {
                int pij_ind = ij - sum(head(KeepConstant, ij));
                stdev(ij) = sqrt(cov(pij_ind, pij_ind));
            }
        }
    } else {
        vector<double> InMa(pow(reqrdnum, 2), 0.0);
        //
        Calc_Recursive_Exp(model_bool, group_num, RiskFail, RiskPairs, totalnum, ntime, R, Rd, Recur_Base, Recur_First, Recur_Second, nthreads, KeepConstant);
        Expected_Inform_Matrix_CaseCon(model_bool, group_num, RiskFail, RiskPairs, totalnum, ntime, R, Rd, Rdd, RdR, RddR, InMa, Recur_Base, Recur_First, Recur_Second, strata_odds, nthreads, KeepConstant, strata_cond);
        //
        NumericVector InMa_vec(reqrdnum * reqrdnum);  //  simplfied information matrix
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
            InMa_vec[ij * reqrdnum + jk] = InMa[ij*reqrdnum+jk];
            InMa_vec[jk * reqrdnum + ij] = InMa[ij * reqrdnum + jk];
        }
        InMa_vec.attr("dim") = Dimension(reqrdnum, reqrdnum);
        const Map<MatrixXd> InMa_mat(as<Map<MatrixXd> >(InMa_vec));    //
        cov = InMa_mat.inverse().matrix();  //  uses inverse information matrix to calculate the standard deviation
        for (int ij = 0; ij < totalnum; ij++) {
            if (KeepConstant[ij] == 0) {
                int pij_ind = ij - sum(head(KeepConstant, ij));
                stdev(ij) = sqrt(cov(pij_ind, pij_ind));
            }
        }
    }
    //
    res_list = List::create(_["LogLik"] = wrap(Ll[0]), _["Deviation"] = wrap(dev), _["First_Der"] = wrap(Lld), _["Second_Der"] = Lldd_vec, _["beta_0"] = wrap(beta_0), _["StrataOdds"] = wrap(strata_odds), _["Standard_Deviation"] = wrap(stdev), _["Covariance"] = wrap(cov), _["Parameter_Lists"] = para_list, _["Control_List"] = control_list, _["Converged"] = convgd, _["FreeParameters"] = wrap(reqrdnum), _["FreeSets"] = wrap(reqrdcond), _["Status"] = "PASSED");
    //  returns a list of results
    return res_list;
}

//' Primary logistic regression with multiple starting points and optional combinations of stratification and no derivative calculation.
//'
//' \code{LogLik_Logist_Omnibus} Performs the calls to calculation functions, Structures the logistic regression, With verbose option prints out time stamps and intermediate sums of terms and derivatives
//'
//' @inheritParams CPP_template
//'
//' @return List of final results: Log-likelihood of optimum, first derivative of log-likelihood, second derivative matrix, parameter list, standard deviation estimate, AIC, model information
//' @noRd
//'
//
List LogLik_Logist_Omnibus(const Ref<const MatrixXd>& CountEvent, IntegerVector term_n, StringVector tform, NumericMatrix& a_ns, Ref<MatrixXd> df0, IntegerVector dfc, int fir, string modelform, double lr, List optim_para, NumericVector maxiters, int guesses, int halfmax, double epsilon, double step_max, double thres_step_max, double deriv_epsilon, int verbose, IntegerVector KeepConstant, int term_tot, int nthreads, List model_bool, const double gmix_theta, const IntegerVector gmix_term, const Ref<const MatrixXd>& Lin_Sys, const Ref<const VectorXd>& Lin_Res) {
    //
    List temp_list = List::create(_["Status"] = "FAILED");  //  used as a dummy return value for code checking
    //
    //  Time durations are measured from this point on in microseconds
    //
    //  df0: covariate data
    //  ntime: number of event times for Cox PH
    //  totalnum: number of terms used
    //
    //  ------------------------------------------------------------------------- //  initialize
    const int mat_row = df0.rows();
    //
    int totalnum = term_n.size();
    int reqrdnum = totalnum - sum(KeepConstant);
    //  cout.precision: controls the number of significant digits printed
    //  nthreads: number of threads used for parallel operations
    //
    Rcout.precision(7);  //  forces higher precision numbers printed to terminal
    //
    //  Lld_worst: stores the highest magnitude log-likelihood derivative
    //
    //
    double Lld_worst = 0.0;  //  stores derivative value used to determine if every parameter is near convergence
    double dbeta_max = 0.0;  //  stores the largest step taken, determines if the step sizes are small enough for convergence
    //
    //  ---------------------------------------------
    //  To Start, needs to seperate the derivative terms
    //  ---------------------------------------------
    //
    //  ------------------------------------------------------------------------- //  initialize
    NumericVector a_n = a_ns.row(0);
    Map<VectorXd> beta_0(as<Map<VectorXd> >(a_n));
    MatrixXd T0;
    MatrixXd Td0;
    MatrixXd Tdd0;
    //
    MatrixXd Te;
    MatrixXd R;
    MatrixXd Rd;
    MatrixXd Rdd;
    MatrixXd P;
    MatrixXd Pnot;
    MatrixXd Pd;
    MatrixXd Pdd;
    //
    MatrixXd Dose;
    MatrixXd nonDose;
    MatrixXd nonDose_LIN;
    MatrixXd nonDose_PLIN;
    MatrixXd nonDose_LOGLIN;
    MatrixXd TTerm;
    MatrixXd RdR;
    MatrixXd RddR;
    MatrixXd PdP;
    MatrixXd PddP;
    MatrixXd PnotdP;
    MatrixXd PnotddP;
    //  ------------------------------------------------------------------------- //  initialize
    //  ---------------------------------------------
    //  To Start, needs to seperate the derivative terms
    //  ---------------------------------------------
    //
    T0 = MatrixXd::Zero(mat_row, totalnum);  //  preallocates matrix for Term column
    Te = MatrixXd::Zero(mat_row, 1);  //  preallocates matrix for column terms used for temporary storage
    R = MatrixXd::Zero(mat_row, 1);  //  preallocates matrix for Risks
    P = MatrixXd::Zero(mat_row, 1);  //  preallocates matrix for Probabilities
    Pnot = MatrixXd::Zero(mat_row, 1);  //  preallocates matrix for 1-P
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
        Pd = MatrixXd::Zero(mat_row, reqrdnum);  //  preallocates matrix for probability derivatives
        Pdd = MatrixXd::Zero(mat_row, reqrdnum*(reqrdnum + 1)/2);  //  preallocates matrix for probability second derivatives
        //
        RdR = MatrixXd::Zero(mat_row, reqrdnum);  //  preallocates matrix for Risk to derivative ratios
        RddR = MatrixXd::Zero(mat_row, reqrdnum*(reqrdnum + 1)/2);  //  preallocates matrix for Risk to second derivative ratios
        PdP = MatrixXd::Zero(mat_row, reqrdnum);  //  preallocates matrix for probability to derivative ratios
        PddP = MatrixXd::Zero(mat_row, reqrdnum*(reqrdnum + 1)/2);  //  preallocates matrix for probability to second derivative ratios
        PnotdP = MatrixXd::Zero(mat_row, reqrdnum);  //  preallocates matrix for 1-probability to derivative ratios
        PnotddP = MatrixXd::Zero(mat_row, reqrdnum*(reqrdnum + 1)/2);  //  preallocates matrix for 1-probability to second derivative ratios
    }
    //  ------------------------------------------------------------------------- //  initialize
    vector<double> Ll(reqrdnum, 0.0);  //  log-likelihood values
    vector<double> Lld(reqrdnum, 0.0);  //  log-likelihood derivative values
    vector<double> Lldd(pow(reqrdnum, 2), 0.0);  //  the second derivative matrix has room for every combination, but only the lower triangle is calculated initially
    MatrixXd dev_temp = MatrixXd::Zero(CountEvent.rows(), 3);
    double dev = 0.0;
    //  the log-likelihood is calculated in parallel over the risk groups
    vector <double> Ll_comp(2, Ll[0]);  //  vector to compare values
    double step_max0 = step_max;
    double thres_step_max0 = thres_step_max;
    //
    vector<double> dbeta(totalnum, 0.0);
    NumericVector m_g_store(reqrdnum);
    NumericVector v_beta_store(reqrdnum);
    //
    //  --------------------------
    //  always starts from initial guess
    //  --------------------------
    vector<double> beta_c(totalnum, 0.0);
    vector<double> beta_a(totalnum, 0.0);
    vector<double> beta_best(totalnum, 0.0);
    vector<double> beta_p(totalnum, 0.0);
    VectorXd::Map(&beta_p[0], beta_0.size()) = beta_0;  //  stores previous parameters
    VectorXd::Map(&beta_c[0], beta_0.size()) = beta_0;  //  stores current parameters
    VectorXd::Map(&beta_a[0], beta_0.size()) = beta_0;  //  stores a refrence value for parameters
    VectorXd::Map(&beta_best[0], beta_0.size()) = beta_0;  //  stores the best parameters
    double halves = 0;  //  number of half-steps taken
    int ind0 = fir;  //  used for validations
    int iteration = 0;  //  iteration number
    int maxiter = 0;
    //
    bool convgd = FALSE;
    int iter_stop  = 0;  //  tracks if the iterations should be stopped for convergence
    //
    NumericMatrix beta_fin(a_ns.rows(), a_ns.cols());
    NumericVector LL_fin(a_ns.rows());
    //
    double Ll_abs_best = 10;
    vector<double> beta_abs_best(totalnum, 0.0);
    int guess_abs_best = - 1;
    double Ll_iter_best = 10;
    //  Variables that are used for the risk check function shared across cox, poisson, and log bound functions
    VectorXd cens_weight(1);
    MatrixXd Lls1 = MatrixXd::Zero(1, 1);
    MatrixXd Lls2 = MatrixXd::Zero(1, 1);
    MatrixXd Lls3 = MatrixXd::Zero(1, 1);
    IntegerMatrix RiskFail(1);
    vector<vector<int> > RiskPairs;
    vector<vector<vector<int> > > RiskPairs_Strata;
    MatrixXd Rls1 = MatrixXd::Zero(1, 1);
    MatrixXd Rls2 = MatrixXd::Zero(1, 1);
    MatrixXd Rls3 = MatrixXd::Zero(1, 1);
    NumericVector Strata_vals(1);
    string ties_method = "temp";
    ///
    for (int guess = 0; guess <guesses; guess++) {
        fill(Ll.begin(), Ll.end(), 0.0);
        if (!model_bool["single"]) {
            fill(Lld.begin(), Lld.end(), 0.0);
            fill(Lldd.begin(), Lldd.end(), 0.0);
        }
        if (model_bool["gradient"]) {
            m_g_store.fill(0);
            v_beta_store.fill(0);
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
        Ll_iter_best = 10;
        //
        convgd = FALSE;
        iter_stop  = 0;  //  tracks if the iterations should be stopped for convergence
        //
        maxiter = maxiters[guess];
        a_n = a_ns.row(guess);
        for (int i = 0; i < beta_0.size(); i++) {
            beta_0[i] = a_n[i];
        }
        //
        Cox_Term_Risk_Calc(modelform, tform, term_n, totalnum, fir, dfc, term_tot, T0, Td0, Tdd0, Te, R, Rd, Rdd, Dose, nonDose, beta_0, df0, thres_step_max, step_max, TTerm, nonDose_LIN, nonDose_PLIN, nonDose_LOGLIN, RdR, RddR, nthreads, KeepConstant, verbose, model_bool, gmix_theta, gmix_term);
        LinkCovertRP(model_bool, reqrdnum, R, Rd, Rdd, RdR, RddR, P, Pd, Pdd, Pnot, PdP, PddP, PnotdP, PnotddP);
        //
        //
        if ((P.minCoeff() <= 0) || (P.maxCoeff() >= 1) || (P.hasNaN())) {
            if (verbose >= 1) {
                if (P.minCoeff() <= 0) {
                Rcout << "C++ Error: An invalid probability was detected: " << P.minCoeff() << endl;
                }
                if (P.maxCoeff() >= 1) {
                Rcout << "C++ Error: An invalid probability was detected: " << P.maxCoeff() << endl;
                }
                if (P.hasNaN()) {
                    Rcout << "C++ Error: An invalid probability was detected" << endl;
                }
                Rcout << "C++ Warning: final failing values ";
                for (int ijk = 0; ijk < totalnum; ijk++) {
                    Rcout << beta_0[ijk] << " ";
                }
                Rcout << " " << endl;
                //
                Rcout << "C++ Warning: final failing terms ";
                for (int ijk = 0; ijk < totalnum; ijk++) {
                    Rcout << tform[ijk] << " ";
                }
                Rcout << " " << endl;
            }
            temp_list = List::create(_["beta_0"] = wrap(beta_0), _["Deviation"] = R_NaN, _["Status"] = "FAILED_WITH_NEGATIVE_RISK", _["LogLik"] = R_NaN);
            return temp_list;
        }
        //
        //  -------------------------------------------------------------------------------------------
        //
        //  Calculates log-likelihood
        Calc_LogLik_Logist(model_bool, nthreads, totalnum, CountEvent, P, Pnot, Pd, Pdd, PdP, PnotdP, PddP, PnotddP, Ll, Lld, Lldd, KeepConstant);
        //
        for (int i = 0; i < beta_0.size(); i++) {
            beta_c[i] = beta_0[i];
        }
        while ((iteration < maxiter) && (iter_stop == 0)) {
            iteration++;
            beta_p = beta_c;  //
            beta_a = beta_c;  //
            beta_best = beta_c;  //
            //
            //  calculates the initial change in parameter
            if (model_bool["gradient"]) {
                if (model_bool["constraint"]) {
                    Calc_Change_Gradient_Cons(Lin_Sys, Lin_Res, nthreads, model_bool, totalnum, optim_para, iteration, step_max, Ll, Lld, m_g_store, v_beta_store, beta_0, dbeta, KeepConstant);
                } else {
                    Calc_Change_Gradient(nthreads, model_bool, totalnum, optim_para, iteration, step_max, Lld, m_g_store, v_beta_store, dbeta, KeepConstant);
                }
            } else if (model_bool["constraint"]) {
                Calc_Change_Cons(Lin_Sys, Lin_Res, beta_0, nthreads, totalnum, thres_step_max, lr, step_max, Ll, Lld, Lldd, dbeta, tform, thres_step_max, step_max, KeepConstant);
            } else {
                Calc_Change(nthreads, totalnum, thres_step_max, lr, step_max, Ll, Lld, Lldd, dbeta, tform, thres_step_max, step_max, KeepConstant);
            }
            Intercept_Bound(nthreads, totalnum, beta_0, dbeta, dfc, df0, KeepConstant, tform);
            if ((Ll_iter_best > 0) || (Ll_iter_best < Ll[ind0])) {
                Ll_iter_best = Ll[ind0];
            }
            if (model_bool["gradient"]) {
                //
                for (int ij = 0; ij < totalnum; ij++) {
                    beta_0[ij] = beta_a[ij] + dbeta[ij];
                    beta_c[ij] = beta_0[ij];
                }
                //
                Cox_Term_Risk_Calc(modelform, tform, term_n, totalnum, fir, dfc, term_tot, T0, Td0, Tdd0, Te, R, Rd, Rdd, Dose, nonDose, beta_0, df0, thres_step_max, step_max, TTerm, nonDose_LIN, nonDose_PLIN, nonDose_LOGLIN, RdR, RddR, nthreads, KeepConstant, verbose, model_bool, gmix_theta, gmix_term);
                LinkCovertRP(model_bool, reqrdnum, R, Rd, Rdd, RdR, RddR, P, Pd, Pdd, Pnot, PdP, PddP, PnotdP, PnotddP);
                //
                if ((P.minCoeff() <= 0) || (P.maxCoeff() >= 1) || (R.hasNaN()))  {
                    // #ifdef _OPENMP
                    // #pragma omp parallel for schedule(dynamic) num_threads(nthreads)
                    // #endif
                    for (int ijk = 0; ijk < totalnum; ijk++) {
                        dbeta[ijk] = dbeta[ijk] / 2.0;
                    }
                    halves+=0.2;
                } else {
                    halves++;
                    Calc_LogLik_Logist(model_bool, nthreads, totalnum, CountEvent, P, Pnot, Pd, Pdd, PdP, PnotdP, PddP, PnotddP, Ll, Lld, Lldd, KeepConstant);
                    //
                    if (Ll[ind0] <= Ll_abs_best) {  //  if a better point wasn't found, takes a half-step
                        // #ifdef _OPENMP
                        // #pragma omp parallel for schedule(dynamic) num_threads(nthreads)
                        // #endif
                        for (int ijk = 0; ijk < totalnum; ijk++) {
                            dbeta[ijk] = dbeta[ijk] * 0.5;
                        }
                    } else {  //  if improved, updates the best vector
                        // #ifdef _OPENMP
                        // #pragma omp parallel for schedule(dynamic) num_threads(nthreads)
                        // #endif
                        for (int ijk = 0; ijk < totalnum; ijk++) {
                            beta_best[ijk] = beta_c[ijk];
                        }
                    }
                    // #ifdef _OPENMP
                    // #pragma omp parallel for schedule(dynamic) num_threads(nthreads)
                    // #endif
                    for (int ijk = 0; ijk < totalnum; ijk++) {
                        beta_0[ijk] = beta_c[ijk];
                    }
                }
                //
            } else {
                halves = 0;
                while ((Ll[ind0] <= Ll_iter_best) && (halves < halfmax)) {  //  repeats until half-steps maxed or an improvement
                    for (int ij = 0; ij < totalnum; ij++) {
                        beta_0[ij] = beta_a[ij] + dbeta[ij];
                        beta_c[ij] = beta_0[ij];
                    }
                    //  ---------------------------------------------------------------------------------------------------------------------------//
                    //  The same subterm, risk, sides, and log-likelihood calculations are performed every half-step and iteration
                    //  ---------------------------------------------------------------------------------------------------------------------------//
                    Cox_Term_Risk_Calc(modelform, tform, term_n, totalnum, fir, dfc, term_tot, T0, Td0, Tdd0, Te, R, Rd, Rdd, Dose, nonDose, beta_0, df0, thres_step_max, step_max, TTerm, nonDose_LIN, nonDose_PLIN, nonDose_LOGLIN, RdR, RddR, nthreads, KeepConstant, verbose, model_bool, gmix_theta, gmix_term);
                    LinkCovertRP(model_bool, reqrdnum, R, Rd, Rdd, RdR, RddR, P, Pd, Pdd, Pnot, PdP, PddP, PnotdP, PnotddP);
                    //
                    if ((P.minCoeff() <= 0) || (P.maxCoeff() >= 1) || (R.hasNaN()))  {
                        // #ifdef _OPENMP
                        // #pragma omp parallel for schedule(dynamic) num_threads(nthreads)
                        // #endif
                        for (int ijk = 0; ijk < totalnum; ijk++) {
                            dbeta[ijk] = dbeta[ijk] / 2.0;
                        }
                        halves+=0.2;
                    } else {
                        halves++;
                        Calc_LogLik_Logist(model_bool, nthreads, totalnum, CountEvent, P, Pnot, Pd, Pdd, PdP, PnotdP, PddP, PnotddP, Ll, Lld, Lldd, KeepConstant);
                        //
                        if (Ll[ind0] <= Ll_abs_best) {  //  if a better point wasn't found, takes a half-step
                            // #ifdef _OPENMP
                            // #pragma omp parallel for schedule(dynamic) num_threads(nthreads)
                            // #endif
                            for (int ijk = 0; ijk < totalnum; ijk++) {
                                dbeta[ijk] = dbeta[ijk] * 0.5;
                            }
                        } else {  //  if improved, updates the best vector
                            // #ifdef _OPENMP
                            // #pragma omp parallel for schedule(dynamic) num_threads(nthreads)
                            // #endif
                            for (int ijk = 0; ijk < totalnum; ijk++) {
                                beta_best[ijk] = beta_c[ijk];
                            }
                        }
                        // #ifdef _OPENMP
                        // #pragma omp parallel for schedule(dynamic) num_threads(nthreads)
                        // #endif
                        for (int ijk = 0; ijk < totalnum; ijk++) {
                            beta_0[ijk] = beta_c[ijk];
                        }
                    }
                    //
                }
                if (beta_best != beta_c) {  //  if the risk matrices aren't the optimal values, then they must be recalculated
                    //  If it goes through every half step without improvement, then the maximum change needs to be decreased
                    step_max = step_max*pow(0.5, halfmax);  //  reduces the step sizes
                    thres_step_max = thres_step_max*pow(0.5, halfmax);
                    beta_p = beta_best;  //
                    beta_a = beta_best;  //
                    beta_c = beta_best;  //
                    for (int ij = 0; ij < totalnum; ij++) {
                        beta_0[ij] = beta_best[ij];
                    }
                    Cox_Term_Risk_Calc(modelform, tform, term_n, totalnum, fir, dfc, term_tot, T0, Td0, Tdd0, Te, R, Rd, Rdd, Dose, nonDose, beta_0, df0, thres_step_max, step_max, TTerm, nonDose_LIN, nonDose_PLIN, nonDose_LOGLIN, RdR, RddR, nthreads, KeepConstant, verbose, model_bool, gmix_theta, gmix_term);
                    LinkCovertRP(model_bool, reqrdnum, R, Rd, Rdd, RdR, RddR, P, Pd, Pdd, Pnot, PdP, PddP, PnotdP, PnotddP);
                    //  Calculates log-likelihood
                    Calc_LogLik_Logist(model_bool, nthreads, totalnum, CountEvent, P, Pnot, Pd, Pdd, PdP, PnotdP, PddP, PnotddP, Ll, Lld, Lldd, KeepConstant);
                    //
                }
            }
            Lld_worst = 0;
            for (int ij = 0; ij < reqrdnum; ij++) {
                if (abs(Lld[ij]) > Lld_worst) {
                    Lld_worst = abs(Lld[ij]);
                }
            }
            dbeta_max = abs(dbeta[0]);
            for (int ij = 1; ij < totalnum; ij++) {
                if (abs(dbeta[ij]) > dbeta_max) {
                    dbeta_max = abs(dbeta[ij]);
                }
            }
            if (Lld_worst < deriv_epsilon) {  //  ends if the derivatives are low enough
                iter_stop = 1;
            }
            Ll_comp[1] = Ll[0];
            if (step_max < epsilon) {  //  if the maximum change is too low, then it ends
                iter_stop = 1;
            }
            if (dbeta_max < epsilon) {  //  if the maximum change is too low, then it ends
                iter_stop = 1;
            }
            if (model_bool["single"]) {
                iter_stop = 1;
            } else {}
        }
        //  -----------------------------------------------
        //  Performing Full Calculation to get full second derivative matrix
        //  -----------------------------------------------
        fill(Ll.begin(), Ll.end(), 0.0);
        if (!model_bool["single"]) {
            fill(Lld.begin(), Lld.end(), 0.0);
            fill(Lldd.begin(), Lldd.end(), 0.0);
        }
        //  Calculates log-likelihood
        Calc_LogLik_Logist(model_bool, nthreads, totalnum, CountEvent, P, Pnot, Pd, Pdd, PdP, PnotdP, PddP, PnotddP, Ll, Lld, Lldd, KeepConstant);
        //
        a_n = beta_0;
        //
        beta_fin(guess, _) = a_n;
        LL_fin[guess] = Ll[0];
        if ((Ll_abs_best > 0) || (Ll_abs_best < Ll[ind0])) {
            Ll_abs_best = Ll[ind0];
            beta_abs_best = beta_c;
            guess_abs_best = guess;
        }
    }
    //
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
    convgd = FALSE;
    iter_stop  = 0;  //  tracks if the iterations should be stopped for convergence
    //
    int guess_max = guess_abs_best;
    if (verbose >= 3) {
        Rcout << "C++ Note: Guess Results" << endl;
        Rcout << "Guess number, parameter values, Log-Likelihood" << endl;
        NumericVector beta_temp;
        for (int i = 0;  i < guesses; i++) {
            beta_temp = wrap(beta_fin.row(i));
            if (i == guess_max) {
                Rcout << i << ", " << beta_temp << ", " << LL_fin[i] << "<-- Best Guess" << endl;
            } else {
                Rcout << i << ", " << beta_temp << ", " << LL_fin[i] << endl;
            }
        }
    }
    //
    maxiter = maxiters[guesses];
    a_n = beta_abs_best;
    for (int i = 0; i < beta_0.size(); i++) {
        beta_0[i] = a_n[i];
    }
    for (int i = 0; i < beta_0.size(); i++) {
        beta_c[i] = beta_0[i];
    }
    if (model_bool["gradient"]) {
        m_g_store.fill(0);
        v_beta_store.fill(0);
    }
    //  Calculates the subterm and term values
    Cox_Term_Risk_Calc(modelform, tform, term_n, totalnum, fir, dfc, term_tot, T0, Td0, Tdd0, Te, R, Rd, Rdd, Dose, nonDose, beta_0, df0, thres_step_max, step_max, TTerm, nonDose_LIN, nonDose_PLIN, nonDose_LOGLIN, RdR, RddR, nthreads, KeepConstant, verbose, model_bool, gmix_theta, gmix_term);
    LinkCovertRP(model_bool, reqrdnum, R, Rd, Rdd, RdR, RddR, P, Pd, Pdd, Pnot, PdP, PddP, PnotdP, PnotddP);
    //  Calculates log-likelihood
    Calc_LogLik_Logist(model_bool, nthreads, totalnum, CountEvent, P, Pnot, Pd, Pdd, PdP, PnotdP, PddP, PnotddP, Ll, Lld, Lldd, KeepConstant);
    //
    List res_list;
    //
    if (model_bool["single"]) {
        res_list = List::create(_["LogLik"] = wrap(Ll[0]), _["beta_0"] = wrap(beta_0), _["AIC"] = 2*(totalnum-accumulate(KeepConstant.begin(), KeepConstant.end(), 0.0)), _["BIC"] = (totalnum-accumulate(KeepConstant.begin(), KeepConstant.end(), 0.0))*log(mat_row)-2*Ll[0], _["Status"] = "PASSED");
        //  returns a list of results
        return res_list;
    }
    for (int i = 0; i < beta_0.size(); i++) {
        beta_c[i] = beta_0[i];
    }
    while ((iteration < maxiter) && (iter_stop == 0)) {
        iteration++;
        beta_p = beta_c;  //
        beta_a = beta_c;  //
        beta_best = beta_c;  //
        //  calculates the initial change in parameter
        if (model_bool["gradient"]) {
            if (model_bool["constraint"]) {
                Calc_Change_Gradient_Cons(Lin_Sys, Lin_Res, nthreads, model_bool, totalnum, optim_para, iteration, step_max, Ll, Lld, m_g_store, v_beta_store, beta_0, dbeta, KeepConstant);
            } else {
                Calc_Change_Gradient(nthreads, model_bool, totalnum, optim_para, iteration, step_max, Lld, m_g_store, v_beta_store, dbeta, KeepConstant);
            }
        } else if (model_bool["constraint"]) {
            Calc_Change_Cons(Lin_Sys, Lin_Res, beta_0, nthreads, totalnum, thres_step_max, lr, step_max, Ll, Lld, Lldd, dbeta, tform, thres_step_max, step_max, KeepConstant);
        } else {
            Calc_Change(nthreads, totalnum, thres_step_max, lr, step_max, Ll, Lld, Lldd, dbeta, tform, thres_step_max, step_max, KeepConstant);
        }
        Intercept_Bound(nthreads, totalnum, beta_0, dbeta, dfc, df0, KeepConstant, tform);
        //
        if ((Ll_abs_best > 0) || (Ll_abs_best < Ll[ind0])) {
            Ll_abs_best = Ll[ind0];
            beta_abs_best = beta_c;
        }
        //
        if (model_bool["gradient"]) {
            //
            for (int ij = 0; ij < totalnum; ij++) {
                beta_0[ij] = beta_a[ij] + dbeta[ij];
                beta_c[ij] = beta_0[ij];
            }
            //
            Cox_Term_Risk_Calc(modelform, tform, term_n, totalnum, fir, dfc, term_tot, T0, Td0, Tdd0, Te, R, Rd, Rdd, Dose, nonDose, beta_0, df0, thres_step_max, step_max, TTerm, nonDose_LIN, nonDose_PLIN, nonDose_LOGLIN, RdR, RddR, nthreads, KeepConstant, verbose, model_bool, gmix_theta, gmix_term);
            LinkCovertRP(model_bool, reqrdnum, R, Rd, Rdd, RdR, RddR, P, Pd, Pdd, Pnot, PdP, PddP, PnotdP, PnotddP);
            //
            //
            if ((P.minCoeff() <= 0) || (P.maxCoeff() >= 1) || (R.hasNaN()))  {
                // #ifdef _OPENMP
                // #pragma omp parallel for schedule(dynamic) num_threads(nthreads)
                // #endif
                for (int ijk = 0; ijk < totalnum; ijk++) {
                    dbeta[ijk] = dbeta[ijk] / 2.0;
                }
                halves+=0.2;
            } else {
                halves++;
                Calc_LogLik_Logist(model_bool, nthreads, totalnum, CountEvent, P, Pnot, Pd, Pdd, PdP, PnotdP, PddP, PnotddP, Ll, Lld, Lldd, KeepConstant);
                //
                if (Ll[ind0] <= Ll_abs_best) {  //  if a better point wasn't found, takes a half-step
                    // #ifdef _OPENMP
                    // #pragma omp parallel for schedule(dynamic) num_threads(nthreads)
                    // #endif
                    for (int ijk = 0; ijk < totalnum; ijk++) {
                        dbeta[ijk] = dbeta[ijk] * 0.5;
                    }
                } else {  //  if improved, updates the best vector
                    // #ifdef _OPENMP
                    // #pragma omp parallel for schedule(dynamic) num_threads(nthreads)
                    // #endif
                    for (int ijk = 0; ijk < totalnum; ijk++) {
                        beta_best[ijk] = beta_c[ijk];
                    }
                }
                // #ifdef _OPENMP
                // #pragma omp parallel for schedule(dynamic) num_threads(nthreads)
                // #endif
                for (int ijk = 0; ijk < totalnum; ijk++) {
                    beta_0[ijk] = beta_c[ijk];
                }
            }
            //
        } else {
            halves = 0;
            while ((Ll[ind0] <= Ll_abs_best) && (halves < halfmax)) {  //  repeats until half-steps maxed or an improvement
                for (int ij = 0; ij < totalnum; ij++) {
                    beta_0[ij] = beta_a[ij] + dbeta[ij];
                    beta_c[ij] = beta_0[ij];
                }
                //  -----------------------------------------------------------------------------------------------------------------------------------//
                //  The same subterm, risk, sides, and log-likelihood calculations are performed every half-step and iteration
                //  -----------------------------------------------------------------------------------------------------------------------------------//
                Cox_Term_Risk_Calc(modelform, tform, term_n, totalnum, fir, dfc, term_tot, T0, Td0, Tdd0, Te, R, Rd, Rdd, Dose, nonDose, beta_0, df0, thres_step_max, step_max, TTerm, nonDose_LIN, nonDose_PLIN, nonDose_LOGLIN, RdR, RddR, nthreads, KeepConstant, verbose, model_bool, gmix_theta, gmix_term);
                LinkCovertRP(model_bool, reqrdnum, R, Rd, Rdd, RdR, RddR, P, Pd, Pdd, Pnot, PdP, PddP, PnotdP, PnotddP);
                //
                if ((P.minCoeff() <= 0) || (P.maxCoeff() >= 1) || (R.hasNaN()))  {
                    // #ifdef _OPENMP
                    // #pragma omp parallel for schedule(dynamic) num_threads(nthreads)
                    // #endif
                    for (int ijk = 0; ijk < totalnum; ijk++) {
                        dbeta[ijk] = dbeta[ijk] / 2.0;
                    }
                    halves+=0.2;
                } else {
                    halves++;
                    Calc_LogLik_Logist(model_bool, nthreads, totalnum, CountEvent, P, Pnot, Pd, Pdd, PdP, PnotdP, PddP, PnotddP, Ll, Lld, Lldd, KeepConstant);
                    //
                    if (Ll[ind0] <= Ll_abs_best) {  //  if a better point wasn't found, takes a half-step
                        // #ifdef _OPENMP
                        // #pragma omp parallel for schedule(dynamic) num_threads(nthreads)
                        // #endif
                        for (int ijk = 0; ijk < totalnum; ijk++) {
                            dbeta[ijk] = dbeta[ijk] * 0.5;
                        }
                    } else {  //  if improved, updates the best vector
                        // #ifdef _OPENMP
                        // #pragma omp parallel for schedule(dynamic) num_threads(nthreads)
                        // #endif
                        for (int ijk = 0; ijk < totalnum; ijk++) {
                            beta_best[ijk] = beta_c[ijk];
                        }
                    }
                    // #ifdef _OPENMP
                    // #pragma omp parallel for schedule(dynamic) num_threads(nthreads)
                    // #endif
                    for (int ijk = 0; ijk < totalnum; ijk++) {
                        beta_0[ijk] = beta_c[ijk];
                    }
                }
                //
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
                LinkCovertRP(model_bool, reqrdnum, R, Rd, Rdd, RdR, RddR, P, Pd, Pdd, Pnot, PdP, PddP, PnotdP, PnotddP);
                //  Calculates log-likelihood
                Calc_LogLik_Logist(model_bool, nthreads, totalnum, CountEvent, P, Pnot, Pd, Pdd, PdP, PnotdP, PddP, PnotddP, Ll, Lld, Lldd, KeepConstant);
            }
        }
        Lld_worst = 0;
        for (int ij = 0; ij < reqrdnum; ij++) {
            if (abs(Lld[ij]) > Lld_worst) {
                Lld_worst = abs(Lld[ij]);
            }
        }
        dbeta_max = abs(dbeta[0]);
        for (int ij = 1; ij < totalnum; ij++) {
            if (abs(dbeta[ij]) > dbeta_max) {
                dbeta_max = abs(dbeta[ij]);
            }
        }
        if (Lld_worst < deriv_epsilon) {  //  ends if the derivatives are low enough
            iter_stop = 1;
            convgd = TRUE;
        }
        Ll_comp[1] = Ll[0];
        if (step_max < epsilon) {  //  if the maximum change is too low, then it ends
            iter_stop = 1;
        }
        if (dbeta_max < epsilon) {  //  if the maximum change is too low, then it ends
            iter_stop = 1;
        }
        if (model_bool["single"]) {
            iter_stop = 1;
        }
    }
    if (Lld_worst < deriv_epsilon) {  //  ends if the derivatives are low enough
        iter_stop = 1;
        convgd = TRUE;
    }
    //  -----------------------------------------------
    //  Performing Full Calculation to get full second derivative matrix
    //  -----------------------------------------------
    bool true_gradient = model_bool["gradient"];
    model_bool["gradient"] = false;
    Cox_Term_Risk_Calc(modelform, tform, term_n, totalnum, fir, dfc, term_tot, T0, Td0, Tdd0, Te, R, Rd, Rdd, Dose, nonDose, beta_0, df0, thres_step_max, step_max, TTerm, nonDose_LIN, nonDose_PLIN, nonDose_LOGLIN, RdR, RddR, nthreads, KeepConstant, verbose, model_bool, gmix_theta, gmix_term);
    LinkCovertRP(model_bool, reqrdnum, R, Rd, Rdd, RdR, RddR, P, Pd, Pdd, Pnot, PdP, PddP, PnotdP, PnotddP);
    //
    Calc_LogLik_Logist(model_bool, nthreads, totalnum, CountEvent, P, Pnot, Pd, Pdd, PdP, PnotdP, PddP, PnotddP, Ll, Lld, Lldd, KeepConstant);
    //
//    dev_temp.col(0) = CountEvent.col(0).array() * ( CountEvent.col(0).array() * P.col(0).array().pow(-1).array() ).array().log().array();
//    dev_temp.col(0) = (dev_temp.col(0).array().isFinite()).select(dev_temp.col(0), 0);
//    dev_temp.col(1) = CountEvent.col(1).array() - CountEvent.col(0).array();
//    dev_temp.col(2) = dev_temp.col(1).array() * ( dev_temp.col(1).array()* CountEvent.col(1).array().pow(-1).array() * Pnot.col(0).array().pow(-1).array() ).array().log().array();
//    dev_temp.col(2) = (dev_temp.col(2).array().isFinite()).select(dev_temp.col(2), 0);
//    dev = (dev_temp.col(0).array() + dev_temp.col(2).array()).array().sum();
    dev = -2*Ll[0];
    //
    model_bool["gradient"] = true_gradient;
    //
    if ((Ll_abs_best > 0) || (Ll_abs_best < Ll[ind0])) {
        Ll_abs_best = Ll[ind0];
        beta_abs_best = beta_c;
    }
    //
    List para_list = List::create(_["term_n"] = term_n, _["tforms"] = tform);  //  stores the term information
    List control_list = List::create(_["Iteration"] = iteration, _["Maximum Step"] = dbeta_max, _["Derivative Limiting"] = Lld_worst);  //  stores the total number of iterations used
    //
    NumericVector Lldd_vec(reqrdnum * reqrdnum);
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
    }
    Lldd_vec.attr("dim") = Dimension(reqrdnum, reqrdnum);
    const Map<MatrixXd> Lldd_mat(as<Map<MatrixXd> >(Lldd_vec));
    //
    MatrixXd cov;
    VectorXd stdev = VectorXd::Zero(totalnum);
    if (model_bool["observed_info"]) {
        cov = - 1 * Lldd_mat.inverse().matrix();  //  uses inverse information matrix to calculate the standard deviation
        for (int ij = 0; ij < totalnum; ij++) {
            if (KeepConstant[ij] == 0) {
                int pij_ind = ij - sum(head(KeepConstant, ij));
                stdev(ij) = sqrt(cov(pij_ind, pij_ind));
            }
        }
    } else {
        vector<double> InMa(pow(reqrdnum, 2), 0.0);
        Expected_Inform_Matrix_Logist(nthreads, totalnum, CountEvent, PdP, PnotdP, InMa, KeepConstant);
        NumericVector InMa_vec(reqrdnum * reqrdnum);  //  simplfied information matrix
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
            InMa_vec[ij * reqrdnum + jk] = InMa[ij*reqrdnum+jk];
            InMa_vec[jk * reqrdnum + ij] = InMa[ij * reqrdnum + jk];
        }
        InMa_vec.attr("dim") = Dimension(reqrdnum, reqrdnum);
        const Map<MatrixXd> InMa_mat(as<Map<MatrixXd> >(InMa_vec));    //
        cov = InMa_mat.inverse().matrix();  //  uses inverse information matrix to calculate the standard deviation
        for (int ij = 0; ij < totalnum; ij++) {
            if (KeepConstant[ij] == 0) {
                int pij_ind = ij - sum(head(KeepConstant, ij));
                stdev(ij) = sqrt(cov(pij_ind, pij_ind));
            }
        }
    }
    //
    res_list = List::create(_["LogLik"] = wrap(Ll[0]), _["Deviation"] = wrap(dev), _["First_Der"] = wrap(Lld), _["Second_Der"] = Lldd_vec, _["beta_0"] = wrap(beta_0), _["Standard_Deviation"] = wrap(stdev), _["Covariance"] = wrap(cov), _["AIC"] = 2*(totalnum-accumulate(KeepConstant.begin(), KeepConstant.end(), 0.0))-2*Ll[0], _["BIC"] = (totalnum-accumulate(KeepConstant.begin(), KeepConstant.end(), 0.0))*log(mat_row)-2*Ll[0], _["Parameter_Lists"] = para_list, _["Control_List"] = control_list, _["Converged"] = convgd, _["Status"] = "PASSED");
    //  returns a list of results
    return res_list;
}

