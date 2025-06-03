#include <RcppEigen.h>
#ifdef _OPENMP
#include <omp.h>
#endif
#include "Main_Functions.h"
#include "Omnibus_Pieces.h"
#include "Calc_Repeated.h"
#include "Subterms_Risk.h"
#include "Step_Calc.h"
#include "Colossus_types.h"
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <chrono>
#include <random>
#include <ctime>
#include <Eigen/Core>


// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::plugins(openmp)]]
using namespace std;
using namespace Rcpp;
using namespace Eigen;
using namespace std::chrono;

using Eigen::Map;
using Eigen::MatrixXd;
using Eigen::SparseMatrix;
using Eigen::VectorXd;
using Rcpp::as;

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
// [[Rcpp::export]]
bool Check_Risk(IntegerVector term_n, StringVector tform, NumericVector a_n, NumericMatrix& x_all, IntegerVector dfc, int fir, string modelform, int verbose, IntegerVector KeepConstant, int term_tot, int nthreads, const double gmix_theta, const IntegerVector gmix_term) {
    //
    List temp_list = List::create(_["Status"] = "FAILED");  // used as a dummy return value for code checking
    if (verbose >= 3) {
        Rcout << "C++ Note: START_RISK_CHECK" << endl;
    }
    //
    const Map<MatrixXd> df0(as<Map<MatrixXd> >(x_all));
    //
    int totalnum = term_n.size();
    //
    Rcout.precision(7);  // forces higher precision numbers printed to terminal
    // ---------------------------------------------
    // To Start, needs to seperate the derivative terms
    // ---------------------------------------------
    //
    Map<VectorXd> beta_0(as<Map<VectorXd> >(a_n));
    const int mat_row = df0.rows();
    MatrixXd T0 = MatrixXd::Zero(mat_row, totalnum);  // preallocates matrix for Term column
    //
    MatrixXd Te = MatrixXd::Zero(mat_row, 1);  // preallocates matrix for column terms used for temporary storage
    MatrixXd R = MatrixXd::Zero(mat_row, 1);  // preallocates matrix for Risks
    //
    MatrixXd Dose = MatrixXd::Constant(mat_row, term_tot, 0.0);  // matrix of the total dose term values
    MatrixXd nonDose = MatrixXd::Constant(mat_row, term_tot, 1.0);  // matrix of the total non-dose term values
    MatrixXd nonDose_LIN = MatrixXd::Constant(mat_row, term_tot, 0.0);  // matrix of Linear subterm values
    MatrixXd nonDose_PLIN = MatrixXd::Constant(mat_row, term_tot, 1.0);  // matrix of Loglinear subterm values
    MatrixXd nonDose_LOGLIN = MatrixXd::Constant(mat_row, term_tot, 1.0);  // matrix of Product linear subterm values
    MatrixXd TTerm = MatrixXd::Zero(mat_row, term_tot);  // matrix of term values
    //
    // Calculates the subterm and term values
    Make_subterms_Single(totalnum, term_n, tform, dfc, fir, T0, Dose, nonDose, TTerm, nonDose_LIN, nonDose_PLIN, nonDose_LOGLIN, beta_0, df0, nthreads, KeepConstant);
    // ---------------------------------------------------------
    // Prints off a series of calculations to check at what point values are changing
    // ---------------------------------------------------------
    // Calculates the risk for each row
    Make_Risks_Single(modelform, tform, term_n, totalnum, fir, T0, Te, R, Dose, nonDose, TTerm, nonDose_LIN, nonDose_PLIN, nonDose_LOGLIN, nthreads, KeepConstant, gmix_theta, gmix_term);
    //
    // Removes infinite values
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
// [[Rcpp::export]]
List LogLik_Cox_PH_Omnibus(IntegerVector term_n, StringVector tform, NumericMatrix& a_ns, NumericMatrix& x_all, IntegerVector dfc, int fir, string modelform, double lr, List optim_para, NumericVector maxiters, int guesses, int halfmax, double epsilon, double abs_max, double dose_abs_max, double deriv_epsilon, NumericMatrix df_groups, NumericVector tu, int double_step, int verbose, IntegerVector KeepConstant, int term_tot, string ties_method, int nthreads, NumericVector& Strata_vals, const VectorXd& cens_weight, List model_bool, const double gmix_theta, const IntegerVector gmix_term, const MatrixXd Lin_Sys, const VectorXd Lin_Res) {
    //
    List temp_list = List::create(_["Status"] = "TEMP");  // used as a dummy return value for code checking
    // Time durations are measured from this point on in microseconds
//    time_point<system_clock> start_point, end_point;
//    start_point = system_clock::now();
//    auto start = time_point_cast<microseconds>(start_point).time_since_epoch().count();
//    end_point = system_clock::now();
//    auto ending = time_point_cast<microseconds>(end_point).time_since_epoch().count();  // the time duration is tracked
    //
    // df0: covariate data
    // ntime: number of event times for Cox PH
    // totalnum: number of terms used
    //
    // ------------------------------------------------------------------------- // initialize
    const Map<MatrixXd> df0(as<Map<MatrixXd> >(x_all));
    const int mat_row = df0.rows();
    int ntime = tu.size();
    int totalnum;
    int reqrdnum;
    // ------------------------------------------------------------------------- // initialize
    if (!model_bool["null"]) {
        totalnum = term_n.size();
        reqrdnum = totalnum - sum(KeepConstant);
    } else {
        totalnum = 1;
        reqrdnum = 1;
    }
    //
    // cout.precision: controls the number of significant digits printed
    // nthreads: number of threads used for parallel operations
    //
    Rcout.precision(7);  // forces higher precision numbers printed to terminal
    //
    // Lld_worst: stores the highest magnitude log-likelihood derivative
    //
    //
    double Lld_worst = 0.0;  // stores derivative value used to determine if every parameter is near convergence
    //
    // ---------------------------------------------
    // To Start, needs to seperate the derivative terms
    // ---------------------------------------------
    //
    // ------------------------------------------------------------------------- // initialize
    NumericVector a_n = a_ns.row(0);
    Map<VectorXd> beta_0(as<Map<VectorXd> >(a_n));
    MatrixXd T0;
    MatrixXd Td0;
    MatrixXd Tdd0;
    //
    MatrixXd Te;
    MatrixXd R;
    ColXd Rd;
    ColXd Rdd;
    //
    MatrixXd Dose;
    MatrixXd nonDose;
    MatrixXd nonDose_LIN;
    MatrixXd nonDose_PLIN;
    MatrixXd nonDose_LOGLIN;
    MatrixXd TTerm;
    double dint = 0.0;  // the amount of change used to calculate derivatives in threshold paramters
    double dslp = 0.0;
    ColXd RdR;
    ColXd RddR;
    // ------------------------------------------------------------------------- // initialize
    if (!model_bool["null"]) {
        // ---------------------------------------------
        // To Start, needs to seperate the derivative terms
        // ---------------------------------------------
        //
        Cox_Refresh_R_TERM(totalnum, reqrdnum, term_tot, dint, dslp, dose_abs_max, abs_max, df0, T0, Td0, Tdd0, Te, R, Rd, Rdd, Dose, nonDose, TTerm, nonDose_LIN, nonDose_PLIN, nonDose_LOGLIN, RdR, RddR, model_bool);
    } else {
        R = MatrixXd::Constant(mat_row, 1, 1.0);
    }
    // ------------------------------------------------------------------------- // initialize
    IntegerMatrix RiskFail;
    vector<vector<int> > RiskPairs(ntime);
    vector<vector<vector<int> > > RiskPairs_Strata(ntime, vector<vector<int>>(Strata_vals.size()));
    const Map<MatrixXd> df_m(as<Map<MatrixXd> >(df_groups));
    //
    // ------------------------------------------------------------------------- // initialize
    if (model_bool["strata"]) {
        RiskFail = IntegerMatrix(ntime, 2*Strata_vals.size());  // vector giving the event rows
        // Creates matrices used to identify the event risk groups
        if (model_bool["cr"]) {
            Make_Groups_Strata_CR(ntime, df_m, RiskFail, RiskPairs_Strata, tu, nthreads, Strata_vals, cens_weight);
        } else {
            Make_Groups_Strata(ntime, df_m, RiskFail, RiskPairs_Strata, tu, nthreads, Strata_vals);
        }
    } else {
        RiskFail = IntegerMatrix(ntime, 2);  // vector giving the event rows
        // Creates matrices used to identify the event risk groups
        if (model_bool["cr"]) {
            Make_Groups_CR(ntime, df_m, RiskFail, RiskPairs, tu, cens_weight, nthreads);
        } else {
            Make_Groups(ntime, df_m, RiskFail, RiskPairs, tu, nthreads);
        }
    }
    // ------------------------------------------------------------------------- // initialize
    MatrixXd Rls1;
    MatrixXd Lls1;
    MatrixXd Rls2;
    MatrixXd Rls3;
    MatrixXd Lls2;
    MatrixXd Lls3;
    vector<double> Ll(reqrdnum, 0.0);  // log-likelihood values
    vector<double> Lld(reqrdnum, 0.0);  // log-likelihood derivative values
    vector<double> Lldd(pow(reqrdnum, 2), 0.0);  // the second derivative matrix has room for every combination, but only the lower triangle is calculated initially
    // ------------------------------------------------------------------------- // initialize
    if (model_bool["null"]) {
        if (model_bool["strata"]) {
            Rls1 = MatrixXd::Zero(ntime, Strata_vals.size());  // precomputes a series of sums used frequently in the log-liklihood calculations
            Lls1 = MatrixXd::Zero(ntime, Strata_vals.size());
            Calculate_Null_Sides_Strata(RiskFail, RiskPairs_Strata, ntime, R, Rls1, Lls1, Strata_vals, nthreads);
            Calc_Null_LogLik_Strata(nthreads, RiskFail, RiskPairs_Strata, ntime, R, Rls1, Lls1, Strata_vals, Ll, ties_method);
        } else {
            Rls1 = MatrixXd::Zero(ntime, 1);  // precomputes a series of sums used frequently in the log-liklihood calculations
            Lls1 = MatrixXd::Zero(ntime, 1);
            // the log-likelihood is calculated in parallel over the risk groups
            //
            Calculate_Null_Sides(RiskFail, RiskPairs, ntime, R, Rls1, Lls1, nthreads);
            Calc_Null_LogLik(nthreads, RiskFail, RiskPairs, ntime, R, Rls1, Lls1, Ll, ties_method);
        }
        //
        List res_list = List::create(_["LogLik"] = wrap(Ll[0]), _["AIC"]=-2*Ll[0], _["BIC"]=-2*Ll[0], _["Status"] = "PASSED");
        // returns a list of results
        return res_list;
    }
    Cox_Refresh_R_SIDES(reqrdnum, ntime, Rls1, Rls2, Rls3, Lls1, Lls2, Lls3, Strata_vals, model_bool);
    // the log-likelihood is calculated in parallel over the risk groups
    vector <double> Ll_comp(2, Ll[0]);  // vector to compare values
    double abs_max0 = abs_max;
    double dose_abs_max0 = dose_abs_max;
    //
    vector<double> dbeta(totalnum, 0.0);
    NumericVector m_g_store(reqrdnum);
    NumericVector v_beta_store(reqrdnum);
    //
    // --------------------------
    // always starts from initial guess
    // --------------------------
    vector<double> beta_c(totalnum, 0.0);
    vector<double> beta_a(totalnum, 0.0);
    vector<double> beta_best(totalnum, 0.0);
    vector<double> beta_p(totalnum, 0.0);
    VectorXd::Map(&beta_p[0], beta_0.size()) = beta_0;  // stores previous parameters
    VectorXd::Map(&beta_c[0], beta_0.size()) = beta_0;  // stores current parameters
    VectorXd::Map(&beta_a[0], beta_0.size()) = beta_0;  // stores a refrence value for parameters
    VectorXd::Map(&beta_best[0], beta_0.size()) = beta_0;  // stores the best parameters
    double halves = 0;  // number of half-steps taken
    int ind0 = fir;  // used for validations
    int iteration = 0;  // iteration number
    int maxiter = 0;
    //
    bool convgd = FALSE;
    int iter_stop  = 0;  // tracks if the iterations should be stopped for convergence
    int iter_check = 0;  // signal to check for convergence
    //
    NumericMatrix beta_fin(a_ns.rows(), a_ns.cols());
    NumericVector LL_fin(a_ns.rows());
    //
    double Ll_abs_best = 10;
    vector<double> beta_abs_best(totalnum, 0.0);
    int guess_abs_best = - 1;
    double Ll_iter_best = 10;
    ///
    // Variables that are used for the risk check function shared across cox, poisson, and log bound functions
    double dev = 0.0;
    MatrixXd dev_temp = MatrixXd::Zero(1, 1);
    double Lstar = 0.0;
    MatrixXd PyrC = MatrixXd::Zero(1, 1);
    //
//    end_point = system_clock::now();
//    ending = time_point_cast<microseconds>(end_point).time_since_epoch().count();
//    Rcout << "Prep step: " << ending - start << endl;
//    start_point = system_clock::now();
//    start = time_point_cast<microseconds>(start_point).time_since_epoch().count();
    //
    for (int guess = 0; guess <guesses; guess++) {
        Cox_Refresh_R_TERM(totalnum, reqrdnum, term_tot, dint, dslp, dose_abs_max, abs_max, df0, T0, Td0, Tdd0, Te, R, Rd, Rdd, Dose, nonDose, TTerm, nonDose_LIN, nonDose_PLIN, nonDose_LOGLIN, RdR, RddR, model_bool);
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
        abs_max = abs_max0;
        dose_abs_max = dose_abs_max0;
        iter_stop = 0;
        halves = 0;
        iteration = 0;
        halves = 0;  // number of half-steps taken
        ind0 = fir;  // used for validations
        iteration = 0;  // iteration number
        Ll_iter_best = 10;
        //
        convgd = FALSE;
        iter_stop  = 0;  // tracks if the iterations should be stopped for convergence
        iter_check = 0;  // signal to check for convergence
        //
        maxiter = maxiters[guess];
        a_n = a_ns.row(guess);
        for (int i = 0; i < beta_0.size(); i++) {
            beta_0[i] = a_n[i];
        }
        if (verbose >= 4) {
            Rcout << "C++ Note: starting guess " << guess << endl;
        }
        Cox_Term_Risk_Calc(modelform, tform, term_n, totalnum, fir, dfc, term_tot, T0, Td0, Tdd0, Te, R, Rd, Rdd, Dose, nonDose, beta_0, df0, dose_abs_max, abs_max, TTerm, nonDose_LIN, nonDose_PLIN, nonDose_LOGLIN, RdR, RddR, nthreads, KeepConstant, verbose, model_bool, gmix_theta, gmix_term);
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
        // -------------------------------------------------------------------------------------------
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
            // calculates the initial change in parameter
            if (model_bool["basic"]) {
                Calc_Change_Basic(double_step, nthreads, totalnum, lr, abs_max, Ll, Lld, Lldd, dbeta, KeepConstant);
            } else if (model_bool["gradient"]) {
                Calc_Change_Gradient(nthreads, model_bool, totalnum, optim_para, iteration, abs_max, Lld, m_g_store, v_beta_store, dbeta, KeepConstant);
            } else {
                if (model_bool["constraint"]) {
                    Calc_Change_Cons(Lin_Sys, Lin_Res, beta_0, nthreads, totalnum, dose_abs_max, lr, abs_max, Ll, Lld, Lldd, dbeta, tform, dose_abs_max, abs_max, KeepConstant);
                } else {
                    Calc_Change(double_step, nthreads, totalnum, dose_abs_max, lr, abs_max, Ll, Lld, Lldd, dbeta, tform, dose_abs_max, abs_max, KeepConstant);
                }
                Intercept_Bound(nthreads, totalnum, beta_0, dbeta, dfc, df0, KeepConstant, tform);
            }
            if ((Ll_iter_best > 0) || (Ll_iter_best < Ll[ind0])) {
                Ll_iter_best = Ll[ind0];
            }
            //
            if (model_bool["gradient"]){
                //
                for (int ij = 0; ij < totalnum; ij++) {
                    beta_0[ij] = beta_a[ij] + dbeta[ij];
                    beta_c[ij] = beta_0[ij];
                }
                //
                Cox_Term_Risk_Calc(modelform, tform, term_n, totalnum, fir, dfc, term_tot, T0, Td0, Tdd0, Te, R, Rd, Rdd, Dose, nonDose, beta_0, df0, dose_abs_max, abs_max, TTerm, nonDose_LIN, nonDose_PLIN, nonDose_LOGLIN, RdR, RddR, nthreads, KeepConstant, verbose, model_bool, gmix_theta, gmix_term);
                //
                Cox_Pois_Check_Continue(model_bool, beta_0, beta_best, beta_c, cens_weight, dbeta, dev, dev_temp, fir, halfmax, halves, ind0, iter_stop, KeepConstant, Ll, Ll_iter_best, Lld, Lldd, Lls1, Lls2, Lls3, Lstar, nthreads, ntime, PyrC, R, Rd, Rdd, RddR, RdR, reqrdnum, tform, RiskFail,  RiskPairs, RiskPairs_Strata, Rls1, Rls2, Rls3, Strata_vals, term_n, ties_method, totalnum, TTerm, verbose);
            } else {
                halves = 0;
                while ((Ll[ind0] <= Ll_iter_best) && (halves < halfmax)) {  // repeats until half-steps maxed or an improvement
    //                Cox_Refresh_R_TERM(totalnum, reqrdnum, term_tot, dint, dslp, dose_abs_max, abs_max, df0, T0, Td0, Tdd0, Te, R, Rd, Rdd, Dose, nonDose, TTerm, nonDose_LIN, nonDose_PLIN, nonDose_LOGLIN, RdR, RddR, model_bool);
                    for (int ij = 0; ij < totalnum; ij++) {
                        beta_0[ij] = beta_a[ij] + dbeta[ij];
                        beta_c[ij] = beta_0[ij];
                    }
                    // ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------//
                    // The same subterm, risk, sides, and log-likelihood calculations are performed every half-step and iteration
                    // ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------//
                    Cox_Term_Risk_Calc(modelform, tform, term_n, totalnum, fir, dfc, term_tot, T0, Td0, Tdd0, Te, R, Rd, Rdd, Dose, nonDose, beta_0, df0, dose_abs_max, abs_max, TTerm, nonDose_LIN, nonDose_PLIN, nonDose_LOGLIN, RdR, RddR, nthreads, KeepConstant, verbose, model_bool, gmix_theta, gmix_term);
                    //
                    Cox_Pois_Check_Continue(model_bool, beta_0, beta_best, beta_c, cens_weight, dbeta, dev, dev_temp, fir, halfmax, halves, ind0, iter_stop, KeepConstant, Ll, Ll_iter_best, Lld, Lldd, Lls1, Lls2, Lls3, Lstar, nthreads, ntime, PyrC, R, Rd, Rdd, RddR, RdR, reqrdnum, tform, RiskFail,  RiskPairs, RiskPairs_Strata, Rls1, Rls2, Rls3, Strata_vals, term_n, ties_method, totalnum, TTerm, verbose);
                }
                if (beta_best != beta_c) {  // if the risk matrices aren't the optimal values, then they must be recalculated
                    // If it goes through every half step without improvement, then the maximum change needs to be decreased
                    abs_max = abs_max*pow(0.5, halfmax);  // reduces the step sizes
                    dose_abs_max = dose_abs_max*pow(0.5, halfmax);
                    iter_check = 1;
                    //
                    beta_p = beta_best;  //
                    beta_a = beta_best;  //
                    beta_c = beta_best;  //
                    for (int ij = 0; ij < totalnum; ij++) {
                        beta_0[ij] = beta_best[ij];
                    }
                    Cox_Term_Risk_Calc(modelform, tform, term_n, totalnum, fir, dfc, term_tot, T0, Td0, Tdd0, Te, R, Rd, Rdd, Dose, nonDose, beta_0, df0, dose_abs_max, abs_max, TTerm, nonDose_LIN, nonDose_PLIN, nonDose_LOGLIN, RdR, RddR, nthreads, KeepConstant, verbose, model_bool, gmix_theta, gmix_term);
                }
            }
            if ((iteration % (reqrdnum)) || (iter_check == 1)) {  // checks every set number of iterations
                iter_check = 0;
                if (Lld_worst < deriv_epsilon) {  // ends if the derivatives are low enough
                    iter_stop = 1;
                }
                if (abs_max < epsilon/10) {  // if the maximum change is too low, then it ends
                    iter_stop = 1;
                }
            }
        }
        // -----------------------------------------------
        // Performing Full Calculation to get full second derivative matrix
        // -----------------------------------------------
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
    abs_max = abs_max0;
    dose_abs_max = dose_abs_max0;
    iter_stop = 0;
    halves = 0;
    iteration = 0;
    halves = 0;  // number of half-steps taken
    ind0 = fir;  // used for validations
    iteration = 0;  // iteration number
    //
    convgd = FALSE;
    iter_stop  = 0;  // tracks if the iterations should be stopped for convergence
    iter_check = 0;  // signal to check for convergence
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
    // Calculates the subterm and term values
    Cox_Term_Risk_Calc(modelform, tform, term_n, totalnum, fir, dfc, term_tot, T0, Td0, Tdd0, Te, R, Rd, Rdd, Dose, nonDose, beta_0, df0, dose_abs_max, abs_max, TTerm, nonDose_LIN, nonDose_PLIN, nonDose_LOGLIN, RdR, RddR, nthreads, KeepConstant, verbose, model_bool, gmix_theta, gmix_term);
    //
    // -------------------------------------------------------------------------------------------
    // Calculates the side sum terms used
    Cox_Side_LL_Calc(reqrdnum, ntime, tform, RiskFail,  RiskPairs, RiskPairs_Strata, totalnum, fir, R, Rd, Rdd, Rls1, Rls2, Rls3, Lls1, Lls2, Lls3, cens_weight, Strata_vals, beta_0, RdR, RddR, Ll, Lld, Lldd, nthreads, KeepConstant, ties_method, verbose, model_bool, iter_stop);
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
        // calculates the initial change in parameter
        if (model_bool["basic"]) {
            Calc_Change_Basic(double_step, nthreads, totalnum, lr, abs_max, Ll, Lld, Lldd, dbeta, KeepConstant);
        } else if (model_bool["gradient"]) {
                Calc_Change_Gradient(nthreads, model_bool, totalnum, optim_para, iteration, abs_max, Lld, m_g_store, v_beta_store, dbeta, KeepConstant);
        } else {
            if (model_bool["constraint"]) {
                Calc_Change_Cons(Lin_Sys, Lin_Res, beta_0, nthreads, totalnum, dose_abs_max, lr, abs_max, Ll, Lld, Lldd, dbeta, tform, dose_abs_max, abs_max, KeepConstant);
            } else {
                Calc_Change(double_step, nthreads, totalnum, dose_abs_max, lr, abs_max, Ll, Lld, Lldd, dbeta, tform, dose_abs_max, abs_max, KeepConstant);
            }
            Intercept_Bound(nthreads, totalnum, beta_0, dbeta, dfc, df0, KeepConstant, tform);
        }
        //
        if ((Ll_abs_best > 0) || (Ll_abs_best < Ll[ind0])) {
            Ll_abs_best = Ll[ind0];
            beta_abs_best = beta_c;
        }
        //
        if (model_bool["gradient"]){
            //
            for (int ij = 0; ij < totalnum; ij++) {
                beta_0[ij] = beta_a[ij] + dbeta[ij];
                beta_c[ij] = beta_0[ij];
            }
            //
            Cox_Term_Risk_Calc(modelform, tform, term_n, totalnum, fir, dfc, term_tot, T0, Td0, Tdd0, Te, R, Rd, Rdd, Dose, nonDose, beta_0, df0, dose_abs_max, abs_max, TTerm, nonDose_LIN, nonDose_PLIN, nonDose_LOGLIN, RdR, RddR, nthreads, KeepConstant, verbose, model_bool, gmix_theta, gmix_term);
            //
            Cox_Pois_Check_Continue(model_bool, beta_0, beta_best, beta_c, cens_weight, dbeta, dev, dev_temp, fir, halfmax, halves, ind0, iter_stop, KeepConstant, Ll, Ll_abs_best, Lld, Lldd, Lls1, Lls2, Lls3, Lstar, nthreads, ntime, PyrC, R, Rd, Rdd, RddR, RdR, reqrdnum, tform, RiskFail,  RiskPairs, RiskPairs_Strata, Rls1, Rls2, Rls3, Strata_vals, term_n, ties_method, totalnum, TTerm, verbose);
        } else {
            halves = 0;
            while ((Ll[ind0] <= Ll_abs_best) && (halves < halfmax)) {  // repeats until half-steps maxed or an improvement
                for (int ij = 0; ij < totalnum; ij++) {
                    beta_0[ij] = beta_a[ij] + dbeta[ij];
                    beta_c[ij] = beta_0[ij];
                }
                // ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------//
                // The same subterm, risk, sides, and log-likelihood calculations are performed every half-step and iteration
                // ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------//
                Cox_Term_Risk_Calc(modelform, tform, term_n, totalnum, fir, dfc, term_tot, T0, Td0, Tdd0, Te, R, Rd, Rdd, Dose, nonDose, beta_0, df0, dose_abs_max, abs_max, TTerm, nonDose_LIN, nonDose_PLIN, nonDose_LOGLIN, RdR, RddR, nthreads, KeepConstant, verbose, model_bool, gmix_theta, gmix_term);
                Cox_Pois_Check_Continue(model_bool, beta_0, beta_best, beta_c, cens_weight, dbeta, dev, dev_temp, fir, halfmax, halves, ind0, iter_stop, KeepConstant, Ll, Ll_abs_best, Lld, Lldd, Lls1, Lls2, Lls3, Lstar, nthreads, ntime, PyrC, R, Rd, Rdd, RddR, RdR, reqrdnum, tform, RiskFail,  RiskPairs, RiskPairs_Strata, Rls1, Rls2, Rls3, Strata_vals, term_n, ties_method, totalnum, TTerm, verbose);
            }
            if (beta_best != beta_c) {  // if the risk matrices aren't the optimal values, then they must be recalculated
                // If it goes through every half step without improvement, then the maximum change needs to be decreased
                abs_max = abs_max*pow(0.5, halfmax);  // reduces the step sizes
                dose_abs_max = dose_abs_max*pow(0.5, halfmax);
                iter_check = 1;
                beta_p = beta_best;  //
                beta_a = beta_best;  //
                beta_c = beta_best;  //
                for (int ij = 0; ij < totalnum; ij++) {
                    beta_0[ij] = beta_best[ij];
                }
                Cox_Term_Risk_Calc(modelform, tform, term_n, totalnum, fir, dfc, term_tot, T0, Td0, Tdd0, Te, R, Rd, Rdd, Dose, nonDose, beta_0, df0, dose_abs_max, abs_max, TTerm, nonDose_LIN, nonDose_PLIN, nonDose_LOGLIN, RdR, RddR, nthreads, KeepConstant, verbose, model_bool, gmix_theta, gmix_term);
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
        if (iteration > reqrdnum) {  // doesn't check the first several iterations for convergence
            if ((iteration % (reqrdnum)) || (iter_check == 1)) {  // checks every set number of iterations
                iter_check = 0;
                if (Lld_worst < deriv_epsilon) {  // ends if the derivatives are low enough
                    iter_stop = 1;
                    convgd = TRUE;
                }
                Ll_comp[1] = Ll[0];
                if (abs_max < epsilon/10) {  // if the maximum change is too low, then it ends
                    iter_stop = 1;
                }
            }
        }
    }
    if (Lld_worst < deriv_epsilon) {  // ends if the derivatives are low enough
        iter_stop = 1;
        convgd = TRUE;
    }
    // -----------------------------------------------
    // Performing Full Calculation to get full second derivative matrix
    // -----------------------------------------------
    fill(Ll.begin(), Ll.end(), 0.0);
    if (!model_bool["single"]) {
        fill(Lld.begin(), Lld.end(), 0.0);
        fill(Lldd.begin(), Lldd.end(), 0.0);
    }
    Cox_Side_LL_Calc(reqrdnum, ntime, tform, RiskFail,  RiskPairs, RiskPairs_Strata, totalnum, fir, R, Rd, Rdd, Rls1, Rls2, Rls3, Lls1, Lls2, Lls3, cens_weight, Strata_vals, beta_0, RdR, RddR, Ll, Lld, Lldd, nthreads, KeepConstant, ties_method, verbose, model_bool, iter_stop);
    //
    if ((Ll_abs_best > 0) || (Ll_abs_best < Ll[ind0])) {
        Ll_abs_best = Ll[ind0];
        beta_abs_best = beta_c;
    }
    //
    List res_list;
    //
    if (model_bool["single"]) {
        res_list = List::create(_["LogLik"] = wrap(Ll[0]), _["beta_0"] = wrap(beta_0), _["AIC"] = 2*(totalnum-accumulate(KeepConstant.begin(), KeepConstant.end(), 0.0))-2*Ll[0], _["BIC"] = (totalnum-accumulate(KeepConstant.begin(), KeepConstant.end(), 0.0))*log(mat_row)-2*Ll[0], _["Status"] = "PASSED");
        // returns a list of results
        return res_list;
    }
    List para_list;
    if (!model_bool["basic"]) {
        para_list = List::create(_["term_n"] = term_n, _["tforms"] = tform);  // stores the term information
    }
    List control_list = List::create(_["Iteration"] = iteration, _["Maximum Step"]= abs_max, _["Derivative Limiting"] = Lld_worst);  // stores the total number of iterations used
    //
    if (model_bool["gradient"]) {
        res_list = List::create(_["LogLik"] = wrap(Ll[0]), _["First_Der"] = wrap(Lld), _["beta_0"] = wrap(beta_0), _["AIC"] = 2*(totalnum-accumulate(KeepConstant.begin(), KeepConstant.end(), 0.0))-2*Ll[0], _["BIC"] = (totalnum-accumulate(KeepConstant.begin(), KeepConstant.end(), 0.0))*log(mat_row)-2*Ll[0], _["Parameter_Lists"] = para_list, _["Control_List"] = control_list, _["Converged"] = convgd, _["Status"] = "PASSED");
        return res_list;
    }
    //
    NumericVector Lldd_vec(reqrdnum * reqrdnum);  // simplfied information matrix
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
    if (model_bool["oberved_info"]){
        cov = - 1 * Lldd_mat.inverse().matrix();  // uses inverse information matrix to calculate the standard deviation
        for (int ij = 0; ij < totalnum; ij++) {
            if (KeepConstant[ij] == 0) {
                int pij_ind = ij - sum(head(KeepConstant, ij));
                stdev(ij) = sqrt(cov(pij_ind, pij_ind));
            }
        }
    } else {
        vector<double> InMa(pow(reqrdnum, 2), 0.0);
        if (model_bool["strata"]){
            if (model_bool["cr"]){
                Expected_Inform_Matrix_Cox_Strata_CR(nthreads, RiskFail, RiskPairs_Strata, totalnum, ntime, R, Rd, RdR, cens_weight, InMa, Strata_vals, KeepConstant);
            } else {
                Expected_Inform_Matrix_Cox_Strata(nthreads, RiskFail, RiskPairs_Strata, totalnum, ntime, R, Rd, RdR, InMa, Strata_vals, KeepConstant);
            }
        } else {
            if (model_bool["cr"]){
                Expected_Inform_Matrix_Cox_CR(nthreads, RiskFail, RiskPairs, totalnum, ntime, R, Rd, RdR, cens_weight, InMa, KeepConstant);
            } else {
                Expected_Inform_Matrix_Cox(nthreads, RiskFail, RiskPairs, totalnum, ntime, R, Rd, RdR, InMa, KeepConstant);
            }
        }
        NumericVector InMa_vec(reqrdnum * reqrdnum);  // simplfied information matrix
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
        cov = InMa_mat.inverse().matrix();  // uses inverse information matrix to calculate the standard deviation
        for (int ij = 0; ij < totalnum; ij++) {
            if (KeepConstant[ij] == 0) {
                int pij_ind = ij - sum(head(KeepConstant, ij));
                stdev(ij) = sqrt(cov(pij_ind, pij_ind));
            }
        }
    }
//    start_point = system_clock::now();
//    start = time_point_cast<microseconds>(start_point).time_since_epoch().count();
    //
    if (model_bool["basic"]) {
        res_list = List::create(_["LogLik"] = wrap(Ll[0]), _["First_Der"] = wrap(Lld), _["Second_Der"] = Lldd_vec, _["beta_0"] = wrap(beta_0), _["Standard_Deviation"] = wrap(stdev), _["Covariance"] = wrap(cov), _["AIC"] = 2*(totalnum-accumulate(KeepConstant.begin(), KeepConstant.end(), 0.0))-2*Ll[0], _["BIC"] = (totalnum-accumulate(KeepConstant.begin(), KeepConstant.end(), 0.0))*log(df0.rows())-2*Ll[0], _["Control_List"] = control_list, _["Converged"] = convgd, _["Status"] = "PASSED");
    } else {
        res_list = List::create(_["LogLik"] = wrap(Ll[0]), _["First_Der"] = wrap(Lld), _["Second_Der"] = Lldd_vec, _["beta_0"] = wrap(beta_0), _["Standard_Deviation"] = wrap(stdev), _["Covariance"] = wrap(cov), _["AIC"] = 2*(totalnum-accumulate(KeepConstant.begin(), KeepConstant.end(), 0.0))-2*Ll[0], _["BIC"] = (totalnum-accumulate(KeepConstant.begin(), KeepConstant.end(), 0.0))*log(df0.rows())-2*Ll[0], _["Parameter_Lists"] = para_list, _["Control_List"] = control_list, _["Converged"] = convgd, _["Status"] = "PASSED");
    }
    // returns a list of results
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
// [[Rcpp::export]]
List LogLik_Pois_Omnibus(const MatrixXd& PyrC, IntegerVector term_n, StringVector tform, NumericMatrix& a_ns, NumericMatrix& x_all, IntegerVector dfc, int fir, string modelform, double lr, List optim_para, NumericVector maxiters, int guesses, int halfmax, double epsilon, double abs_max, double dose_abs_max, double deriv_epsilon, int double_step, int verbose, IntegerVector KeepConstant, int term_tot, int nthreads, const MatrixXd& dfs, List model_bool, const double gmix_theta, const IntegerVector gmix_term, const MatrixXd Lin_Sys, const VectorXd Lin_Res) {
    //
    List temp_list = List::create(_["Status"] = "FAILED");  // used as a dummy return value for code checking
    //
    // Time durations are measured from this point on in microseconds
    //
    // df0: covariate data
    // ntime: number of event times for Cox PH
    // totalnum: number of terms used
    //
    // ------------------------------------------------------------------------- // initialize
    const Map<MatrixXd> df0(as<Map<MatrixXd> >(x_all));
    const int mat_row = df0.rows();
    //
    int totalnum = term_n.size();
    int reqrdnum = totalnum - sum(KeepConstant);
    // cout.precision: controls the number of significant digits printed
    // nthreads: number of threads used for parallel operations
    //
    Rcout.precision(7);  // forces higher precision numbers printed to terminal
    // int nthreads = Eigen::nbThreads() - 1;  // stores how many threads are allocated
    //
    // Lld_worst: stores the highest magnitude log-likelihood derivative
    //
    //
    double Lld_worst = 0.0;  // stores derivative value used to determine if every parameter is near convergence
    //
    // ---------------------------------------------
    // To Start, needs to seperate the derivative terms
    // ---------------------------------------------
    //
    // ------------------------------------------------------------------------- // initialize
    NumericVector a_n = a_ns.row(0);
    Map<VectorXd> beta_0(as<Map<VectorXd> >(a_n));
    MatrixXd T0;
    MatrixXd Td0;
    MatrixXd Tdd0;
    //
    MatrixXd Te;
    MatrixXd R;
    ColXd Rd;
    ColXd Rdd;
    //
    MatrixXd Dose;
    MatrixXd nonDose;
    MatrixXd nonDose_LIN;
    MatrixXd nonDose_PLIN;
    MatrixXd nonDose_LOGLIN;
    MatrixXd TTerm;
    double dint = 0.0;  // the amount of change used to calculate derivatives in threshold paramters
    double dslp = 0.0;
    ColXd RdR;
    ColXd RddR;
    // ------------------------------------------------------------------------- // initialize
    // ---------------------------------------------
    // To Start, needs to seperate the derivative terms
    // ---------------------------------------------
    //
    T0 = MatrixXd::Zero(mat_row, totalnum);  // preallocates matrix for Term column
    Te = MatrixXd::Zero(mat_row, 1);  // preallocates matrix for column terms used for temporary storage
    R = MatrixXd::Zero(mat_row, 1);  // preallocates matrix for Risks
    //
    Dose = MatrixXd::Constant(mat_row, term_tot, 0.0);  // matrix of the total dose term values
    nonDose = MatrixXd::Constant(mat_row, term_tot, 1.0);  // matrix of the total non-dose term values
    nonDose_LIN = MatrixXd::Constant(mat_row, term_tot, 0.0);  // matrix of Linear subterm values
    nonDose_PLIN = MatrixXd::Constant(mat_row, term_tot, 1.0);  // matrix of Loglinear subterm values
    nonDose_LOGLIN = MatrixXd::Constant(mat_row, term_tot, 1.0);  // matrix of Product linear subterm values
    TTerm = MatrixXd::Zero(mat_row, term_tot);  // matrix of term values
    if (model_bool["single"]) {
    } else {
        Td0 = MatrixXd::Zero(mat_row, reqrdnum);  // preallocates matrix for Term derivative columns
        Tdd0 = MatrixXd::Zero(mat_row, reqrdnum*(reqrdnum + 1)/2);  // preallocates matrix for Term second derivative columns
        //
        Rd = MatrixXd::Zero(mat_row, reqrdnum);  // preallocates matrix for Risk derivatives
        Rdd = MatrixXd::Zero(mat_row, reqrdnum*(reqrdnum + 1)/2);  // preallocates matrix for Risk second derivatives
        //
        dint = dose_abs_max;  // the amount of change used to calculate derivatives in threshold paramters
        dslp = abs_max;
        RdR = MatrixXd::Zero(mat_row, reqrdnum);  // preallocates matrix for Risk to derivative ratios
        RddR = MatrixXd::Zero(mat_row, reqrdnum*(reqrdnum + 1)/2);  // preallocates matrix for Risk to second derivative ratios
    }
    VectorXd s_weights;
    if (model_bool["strata"]) {
        s_weights = VectorXd::Zero(mat_row);
        Gen_Strat_Weight(modelform, dfs, PyrC, s_weights, nthreads, tform, term_n, term_tot);
    }
    // ------------------------------------------------------------------------- // initialize
    vector<double> Ll(reqrdnum, 0.0);  // log-likelihood values
    vector<double> Lld(reqrdnum, 0.0);  // log-likelihood derivative values
    vector<double> Lldd(pow(reqrdnum, 2), 0.0);  // the second derivative matrix has room for every combination, but only the lower triangle is calculated initially
    MatrixXd dev_temp = MatrixXd::Zero(PyrC.rows(), 2);
    double dev = 0.0;
    // the log-likelihood is calculated in parallel over the risk groups
    vector <double> Ll_comp(2, Ll[0]);  // vector to compare values
    double abs_max0 = abs_max;
    double dose_abs_max0 = dose_abs_max;
    //
    vector<double> dbeta(totalnum, 0.0);
    NumericVector m_g_store(reqrdnum);
    NumericVector v_beta_store(reqrdnum);
    //
    // --------------------------
    // always starts from initial guess
    // --------------------------
    vector<double> beta_c(totalnum, 0.0);
    vector<double> beta_a(totalnum, 0.0);
    vector<double> beta_best(totalnum, 0.0);
    vector<double> beta_p(totalnum, 0.0);
    VectorXd::Map(&beta_p[0], beta_0.size()) = beta_0;  // stores previous parameters
    VectorXd::Map(&beta_c[0], beta_0.size()) = beta_0;  // stores current parameters
    VectorXd::Map(&beta_a[0], beta_0.size()) = beta_0;  // stores a refrence value for parameters
    VectorXd::Map(&beta_best[0], beta_0.size()) = beta_0;  // stores the best parameters
    double halves = 0;  // number of half-steps taken
    int ind0 = fir;  // used for validations
    int iteration = 0;  // iteration number
    int maxiter = 0;
    //
    bool convgd = FALSE;
    int iter_stop  = 0;  // tracks if the iterations should be stopped for convergence
    int iter_check = 0;  // signal to check for convergence
    //
    NumericMatrix beta_fin(a_ns.rows(), a_ns.cols());
    NumericVector LL_fin(a_ns.rows());
    //
    double Ll_abs_best = 10;
    vector<double> beta_abs_best(totalnum, 0.0);
    int guess_abs_best = - 1;
    double Ll_iter_best = 10;
    // Variables that are used for the risk check function shared across cox, poisson, and log bound functions
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
        abs_max = abs_max0;
        dose_abs_max = dose_abs_max0;
        iter_check = 0;
        iter_stop = 0;
        halves = 0;
        iteration = 0;
        halves = 0;  // number of half-steps taken
        ind0 = fir;  // used for validations
        iteration = 0;  // iteration number
        Ll_iter_best = 10;
        //
        convgd = FALSE;
        iter_stop  = 0;  // tracks if the iterations should be stopped for convergence
        iter_check = 0;  // signal to check for convergence
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
        // -------------------------------------------------------------------------------------------
        //
        // Calculates log-likelihood
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
            // calculates the initial change in parameter
            if (model_bool["gradient"]) {
                Calc_Change_Gradient(nthreads, model_bool, totalnum, optim_para, iteration, abs_max, Lld, m_g_store, v_beta_store, dbeta, KeepConstant);
            } else if (model_bool["constraint"]) {
                Calc_Change_Cons(Lin_Sys, Lin_Res, beta_0, nthreads, totalnum, dose_abs_max, lr, abs_max, Ll, Lld, Lldd, dbeta, tform, dose_abs_max, abs_max, KeepConstant);
            } else {
                Calc_Change(double_step, nthreads, totalnum, dose_abs_max, lr, abs_max, Ll, Lld, Lldd, dbeta, tform, dose_abs_max, abs_max, KeepConstant);
            }
            Intercept_Bound(nthreads, totalnum, beta_0, dbeta, dfc, df0, KeepConstant, tform);
            if ((Ll_iter_best > 0) || (Ll_iter_best < Ll[ind0])) {
                Ll_iter_best = Ll[ind0];
            }
            if (model_bool["gradient"]){
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
                while ((Ll[ind0] <= Ll_iter_best) && (halves < halfmax)) {  // repeats until half-steps maxed or an improvement
                    for (int ij = 0; ij < totalnum; ij++) {
                        beta_0[ij] = beta_a[ij] + dbeta[ij];
                        beta_c[ij] = beta_0[ij];
                    }
                    // ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------//
                    // The same subterm, risk, sides, and log-likelihood calculations are performed every half-step and iteration
                    // ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------//
                    Pois_Term_Risk_Calc(modelform, tform, term_n, totalnum, fir, dfc, term_tot, T0, Td0, Tdd0, Te, R, Rd, Rdd, Dose, nonDose, beta_0, df0, dint, dslp, TTerm, nonDose_LIN, nonDose_PLIN, nonDose_LOGLIN, RdR, RddR, s_weights, nthreads, KeepConstant, verbose, model_bool, gmix_theta, gmix_term);
                    Cox_Pois_Check_Continue(model_bool, beta_0, beta_best, beta_c, cens_weight, dbeta, dev, dev_temp, fir, halfmax, halves, ind0, iter_stop, KeepConstant, Ll, Ll_iter_best, Lld, Lldd, Lls1, Lls2, Lls3, Lstar, nthreads, ntime, PyrC, R, Rd, Rdd, RddR, RdR, reqrdnum, tform, RiskFail,  RiskPairs, RiskPairs_Strata, Rls1, Rls2, Rls3, Strata_vals, term_n, ties_method, totalnum, TTerm, verbose);
                }
                if (beta_best != beta_c) {  // if the risk matrices aren't the optimal values, then they must be recalculated
                    // If it goes through every half step without improvement, then the maximum change needs to be decreased
                    abs_max = abs_max*pow(0.5, halfmax);  // reduces the step sizes
                    dose_abs_max = dose_abs_max*pow(0.5, halfmax);
                    iter_check = 1;
                    beta_p = beta_best;  //
                    beta_a = beta_best;  //
                    beta_c = beta_best;  //
                    for (int ij = 0; ij < totalnum; ij++) {
                        beta_0[ij] = beta_best[ij];
                    }
                    Pois_Term_Risk_Calc(modelform, tform, term_n, totalnum, fir, dfc, term_tot, T0, Td0, Tdd0, Te, R, Rd, Rdd, Dose, nonDose, beta_0, df0, dint, dslp, TTerm, nonDose_LIN, nonDose_PLIN, nonDose_LOGLIN, RdR, RddR, s_weights, nthreads, KeepConstant, verbose, model_bool, gmix_theta, gmix_term);
                    // Calculates log-likelihood
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
            if (iteration > reqrdnum) {  // doesn't check the first several iterations for convergence
                if ((iteration % (reqrdnum)) || (iter_check == 1)) {  // checks every set number of iterations
                    iter_check = 0;
                    if (Lld_worst < deriv_epsilon) {  // ends if the derivatives are low enough
                        iter_stop = 1;
                    }
                    Ll_comp[1] = Ll[0];
                    if (abs_max < epsilon/10) {  // if the maximum change is too low, then it ends
                        iter_stop = 1;
                    }
                }
            }
            if (model_bool["single"]) {
                iter_stop = 1;
            } else {}
        }
        // -----------------------------------------------
        // Performing Full Calculation to get full second derivative matrix
        // -----------------------------------------------
        fill(Ll.begin(), Ll.end(), 0.0);
        if (!model_bool["single"]) {
            fill(Lld.begin(), Lld.end(), 0.0);
            fill(Lldd.begin(), Lldd.end(), 0.0);
        }
        // Calculates log-likelihood
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
    abs_max = abs_max0;
    dose_abs_max = dose_abs_max0;
    iter_check = 0;
    iter_stop = 0;
    halves = 0;
    iteration = 0;
    halves = 0;  // number of half-steps taken
    ind0 = fir;  // used for validations
    iteration = 0;  // iteration number
    //
    convgd = FALSE;
    iter_stop  = 0;  // tracks if the iterations should be stopped for convergence
    iter_check = 0;  // signal to check for convergence
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
    // Calculates the subterm and term values
    Pois_Term_Risk_Calc(modelform, tform, term_n, totalnum, fir, dfc, term_tot, T0, Td0, Tdd0, Te, R, Rd, Rdd, Dose, nonDose, beta_0, df0, dint, dslp, TTerm, nonDose_LIN, nonDose_PLIN, nonDose_LOGLIN, RdR, RddR, s_weights, nthreads, KeepConstant, verbose, model_bool, gmix_theta, gmix_term);
    // Calculates log-likelihood
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
        // calculates the initial change in parameter
        if (model_bool["gradient"]) {
            Calc_Change_Gradient(nthreads, model_bool, totalnum, optim_para, iteration, abs_max, Lld, m_g_store, v_beta_store, dbeta, KeepConstant);
        } else if (model_bool["constraint"]) {
            Calc_Change_Cons(Lin_Sys, Lin_Res, beta_0, nthreads, totalnum, dose_abs_max, lr, abs_max, Ll, Lld, Lldd, dbeta, tform, dose_abs_max, abs_max, KeepConstant);
        } else {
            Calc_Change(double_step, nthreads, totalnum, dose_abs_max, lr, abs_max, Ll, Lld, Lldd, dbeta, tform, dose_abs_max, abs_max, KeepConstant);
        }
        Intercept_Bound(nthreads, totalnum, beta_0, dbeta, dfc, df0, KeepConstant, tform);
        //
        if ((Ll_abs_best > 0) || (Ll_abs_best < Ll[ind0])) {
            Ll_abs_best = Ll[ind0];
            beta_abs_best = beta_c;
        }
        //
        if (model_bool["gradient"]){
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
            while ((Ll[ind0] <= Ll_abs_best) && (halves < halfmax)) {  // repeats until half-steps maxed or an improvement
                for (int ij = 0; ij < totalnum; ij++) {
                    beta_0[ij] = beta_a[ij] + dbeta[ij];
                    beta_c[ij] = beta_0[ij];
                }
                // ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------//
                // The same subterm, risk, sides, and log-likelihood calculations are performed every half-step and iteration
                // ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------//
                Pois_Term_Risk_Calc(modelform, tform, term_n, totalnum, fir, dfc, term_tot, T0, Td0, Tdd0, Te, R, Rd, Rdd, Dose, nonDose, beta_0, df0, dint, dslp, TTerm, nonDose_LIN, nonDose_PLIN, nonDose_LOGLIN, RdR, RddR, s_weights, nthreads, KeepConstant, verbose, model_bool, gmix_theta, gmix_term);
                Cox_Pois_Check_Continue(model_bool, beta_0, beta_best, beta_c, cens_weight, dbeta, dev, dev_temp, fir, halfmax, halves, ind0, iter_stop, KeepConstant, Ll, Ll_abs_best, Lld, Lldd, Lls1, Lls2, Lls3, Lstar, nthreads, ntime, PyrC, R, Rd, Rdd, RddR, RdR, reqrdnum, tform, RiskFail,  RiskPairs, RiskPairs_Strata, Rls1, Rls2, Rls3, Strata_vals, term_n, ties_method, totalnum, TTerm, verbose);
            }
            if (beta_best != beta_c) {  // if the risk matrices aren't the optimal values, then they must be recalculated
                // If it goes through every half step without improvement, then the maximum change needs to be decreased
                abs_max = abs_max*pow(0.5, halfmax);  // reduces the step sizes
                dose_abs_max = dose_abs_max*pow(0.5, halfmax);
                iter_check = 1;
                //
                beta_p = beta_best;  //
                beta_a = beta_best;  //
                beta_c = beta_best;  //
                for (int ij = 0; ij < totalnum; ij++) {
                    beta_0[ij] = beta_best[ij];
                }
                Pois_Term_Risk_Calc(modelform, tform, term_n, totalnum, fir, dfc, term_tot, T0, Td0, Tdd0, Te, R, Rd, Rdd, Dose, nonDose, beta_0, df0, dint, dslp, TTerm, nonDose_LIN, nonDose_PLIN, nonDose_LOGLIN, RdR, RddR, s_weights, nthreads, KeepConstant, verbose, model_bool, gmix_theta, gmix_term);
                // Calculates log-likelihood
                Pois_Dev_LL_Calc(reqrdnum, totalnum, fir, R, Rd, Rdd, beta_0, RdR, RddR, Ll, Lld, Lldd, PyrC, dev_temp, nthreads, KeepConstant, verbose, model_bool, iter_stop, dev);
            }
        }
        Lld_worst = 0;
        for (int ij = 0; ij < reqrdnum; ij++) {
            if (abs(Lld[ij]) > Lld_worst) {
                Lld_worst = abs(Lld[ij]);
            }
        }
        if (iteration > reqrdnum) {  // doesn't check the first several iterations for convergence
            if ((iteration % (reqrdnum)) || (iter_check == 1)) {  // checks every set number of iterations
                iter_check = 0;
                if (Lld_worst < deriv_epsilon) {  // ends if the derivatives are low enough
                    iter_stop = 1;
                    convgd = TRUE;
                }
                Ll_comp[1] = Ll[0];
                if (abs_max < epsilon/10) {  // if the maximum change is too low, then it ends
                    iter_stop = 1;
                }
            }
        }
        if (model_bool["single"]) {
            iter_stop = 1;
        }
    }
    if (Lld_worst < deriv_epsilon) {  // ends if the derivatives are low enough
        iter_stop = 1;
        convgd = TRUE;
    }
    // -----------------------------------------------
    // Performing Full Calculation to get full second derivative matrix
    // -----------------------------------------------
    Pois_Dev_LL_Calc(reqrdnum, totalnum, fir, R, Rd, Rdd, beta_0, RdR, RddR, Ll, Lld, Lldd, PyrC, dev_temp, nthreads, KeepConstant, verbose, model_bool, iter_stop, dev);
    //
    if ((Ll_abs_best > 0) || (Ll_abs_best < Ll[ind0])) {
        Ll_abs_best = Ll[ind0];
        beta_abs_best = beta_c;
    }
    //
    List res_list;
    //
    if (model_bool["single"]) {
        res_list = List::create(_["LogLik"] = wrap(Ll[0]), _["beta_0"] = wrap(beta_0), _["AIC"] = 2*(totalnum-accumulate(KeepConstant.begin(), KeepConstant.end(), 0.0))+dev, _["BIC"] = (totalnum-accumulate(KeepConstant.begin(), KeepConstant.end(), 0.0))*log(mat_row)-2*Ll[0], _["Status"] = "PASSED");
        // returns a list of results
        return res_list;
    }
    List para_list = List::create(_["term_n"] = term_n, _["tforms"] = tform);  // stores the term information
    List control_list = List::create(_["Iteration"] = iteration, _["Maximum Step"] = abs_max, _["Derivative Limiting"] = Lld_worst);  // stores the total number of iterations used
    //
    if (model_bool["gradient"]) {
        res_list = List::create(_["LogLik"] = wrap(Ll[0]), _["First_Der"] = wrap(Lld), _["beta_0"] = wrap(beta_0), _["AIC"] = 2*(totalnum-accumulate(KeepConstant.begin(), KeepConstant.end(), 0.0))-2*Ll[0], _["BIC"] = (totalnum-accumulate(KeepConstant.begin(), KeepConstant.end(), 0.0))*log(mat_row)-2*Ll[0], _["Parameter_Lists"] = para_list, _["Control_List"] = control_list, _["Converged"] = convgd, _["Status"] = "PASSED");
        return res_list;
    }
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
    if (model_bool["oberved_info"]){
        cov = - 1 * Lldd_mat.inverse().matrix();  // uses inverse information matrix to calculate the standard deviation
        for (int ij = 0; ij < totalnum; ij++) {
            if (KeepConstant[ij] == 0) {
                int pij_ind = ij - sum(head(KeepConstant, ij));
                stdev(ij) = sqrt(cov(pij_ind, pij_ind));
            }
        }
    } else {
        vector<double> InMa(pow(reqrdnum, 2), 0.0);
        Expected_Inform_Matrix_Poisson(nthreads, totalnum, PyrC, R, Rd, RdR, InMa, KeepConstant);
        NumericVector InMa_vec(reqrdnum * reqrdnum);  // simplfied information matrix
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
        cov = InMa_mat.inverse().matrix();  // uses inverse information matrix to calculate the standard deviation
        for (int ij = 0; ij < totalnum; ij++) {
            if (KeepConstant[ij] == 0) {
                int pij_ind = ij - sum(head(KeepConstant, ij));
                stdev(ij) = sqrt(cov(pij_ind, pij_ind));
            }
        }
    }
    //
    res_list = List::create(_["LogLik"] = wrap(Ll[0]), _["First_Der"] = wrap(Lld), _["Second_Der"] = Lldd_vec, _["beta_0"] = wrap(beta_0), _["Standard_Deviation"] = wrap(stdev), _["Covariance"] = wrap(cov), _["AIC"] = 2*(totalnum-accumulate(KeepConstant.begin(), KeepConstant.end(), 0.0))+dev, _["BIC"] = (totalnum-accumulate(KeepConstant.begin(), KeepConstant.end(), 0.0))*log(mat_row)-2*Ll[0], _["Deviation"] = dev, _["Parameter_Lists"] = para_list, _["Control_List"] = control_list, _["Converged"] = convgd, _["Status"] = "PASSED");
    // returns a list of results
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
// [[Rcpp::export]]
List LogLik_Cox_PH_Omnibus_Log_Bound(IntegerVector term_n, StringVector tform, NumericVector a_n, NumericMatrix& x_all, IntegerVector dfc, int fir, string modelform, double lr, NumericVector maxiters, int guesses, int halfmax, double epsilon, double abs_max, double dose_abs_max, double deriv_epsilon, NumericMatrix df_groups, NumericVector tu, int verbose, IntegerVector KeepConstant, int term_tot, string ties_method, int nthreads, NumericVector& Strata_vals, const VectorXd& cens_weight, List model_bool, const double gmix_theta, const IntegerVector gmix_term, const MatrixXd Lin_Sys, const VectorXd Lin_Res, double qchi, int para_number, int maxstep, double mult) {
    //
    List temp_list = List::create(_["Status"] = "TEMP");  // used as a dummy return value for code checking
    // Time durations are measured from this point on in microseconds
    //
    // df0: covariate data
    // ntime: number of event times for Cox PH
    // totalnum: number of terms used
    //
    // ------------------------------------------------------------------------- // initialize
    const Map<MatrixXd> df0(as<Map<MatrixXd> >(x_all));
    const int mat_row = df0.rows();
    int ntime = tu.size();
    int totalnum = term_n.size();
    int reqrdnum = totalnum - sum(KeepConstant);
    // ------------------------------------------------------------------------- // initialize
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
    // cout.precision: controls the number of significant digits printed
    // nthreads: number of threads used for parallel operations
    //
    Rcout.precision(7);  // forces higher precision numbers printed to terminal
    //
    // Lld_worst: stores the highest magnitude log-likelihood derivative
    // ---------------------------------------------
    // To Start, needs to seperate the derivative terms
    // ---------------------------------------------
    //
    // ------------------------------------------------------------------------- // initialize
    Map<VectorXd> beta_0(as<Map<VectorXd> >(a_n));
    MatrixXd T0;
    MatrixXd Td0;
    MatrixXd Tdd0;
    //
    MatrixXd Te;
    MatrixXd R;
    ColXd Rd;
    ColXd Rdd;
    //
    MatrixXd Dose;
    MatrixXd nonDose;
    MatrixXd nonDose_LIN;
    MatrixXd nonDose_PLIN;
    MatrixXd nonDose_LOGLIN;
    MatrixXd TTerm;
    double dint = 0.0;  // the amount of change used to calculate derivatives in threshold paramters
    double dslp = 0.0;
    ColXd RdR;
    ColXd RddR;
    // ------------------------------------------------------------------------- // initialize
    // ---------------------------------------------
    // To Start, needs to seperate the derivative terms
    // ---------------------------------------------
    //
    Cox_Refresh_R_TERM(totalnum, reqrdnum, term_tot, dint, dslp, dose_abs_max, abs_max, df0, T0, Td0, Tdd0, Te, R, Rd, Rdd, Dose, nonDose, TTerm, nonDose_LIN, nonDose_PLIN, nonDose_LOGLIN, RdR, RddR, model_bool);
    // ------------------------------------------------------------------------- // initialize
    IntegerMatrix RiskFail;
    vector<vector<int> > RiskPairs(ntime);
    vector<vector<vector<int> > > RiskPairs_Strata(ntime, vector<vector<int>>(Strata_vals.size()));
    const Map<MatrixXd> df_m(as<Map<MatrixXd> >(df_groups));
    // ------------------------------------------------------------------------- // initialize
    if (model_bool["strata"]) {
        RiskFail = IntegerMatrix(ntime, 2*Strata_vals.size());  // vector giving the event rows
        //
        // Creates matrices used to identify the event risk groups
        if (model_bool["cr"]) {
            Make_Groups_Strata_CR(ntime, df_m, RiskFail, RiskPairs_Strata, tu, nthreads, Strata_vals, cens_weight);
        } else {
            Make_Groups_Strata(ntime, df_m, RiskFail, RiskPairs_Strata, tu, nthreads, Strata_vals);
        }
    } else {
        RiskFail = IntegerMatrix(ntime, 2);  // vector giving the event rows
        //
        // Creates matrices used to identify the event risk groups
        if (model_bool["cr"]) {
            Make_Groups_CR(ntime, df_m, RiskFail, RiskPairs, tu, cens_weight, nthreads);
        } else {
            Make_Groups(ntime, df_m, RiskFail, RiskPairs, tu, nthreads);
        }
    }
    // ------------------------------------------------------------------------- // initialize
    MatrixXd Rls1;
    MatrixXd Lls1;
    MatrixXd Rls2;
    MatrixXd Rls3;
    MatrixXd Lls2;
    MatrixXd Lls3;
    vector<double> Ll(reqrdnum, 0.0);  // log-likelihood values
    vector<double> Lld(reqrdnum, 0.0);  // log-likelihood derivative values
    vector<double> Lldd(pow(reqrdnum, 2), 0.0);  // the second derivative matrix has room for every combination, but only the lower triangle is calculated initially
    // ------------------------------------------------------------------------- // initialize
    Cox_Refresh_R_SIDES(reqrdnum, ntime, Rls1, Rls2, Rls3, Lls1, Lls2, Lls3, Strata_vals, model_bool);
    // the log-likelihood is calculated in parallel over the risk groups
    vector <double> Ll_comp(2, Ll[0]);  // vector to compare values
    double abs_max0 = abs_max;
    double dose_abs_max0 = dose_abs_max;
    //
    vector<double> dbeta(totalnum, 0.0);
    //
    // --------------------------
    // always starts from initial guess
    // --------------------------
    vector<double> beta_c(totalnum, 0.0);
    vector<double> beta_a(totalnum, 0.0);
    vector<double> beta_best(totalnum, 0.0);
    vector<double> beta_peak(totalnum, 0.0);
    vector<double> beta_p(totalnum, 0.0);
    VectorXd::Map(&beta_p[0], beta_0.size()) = beta_0;  // stores previous parameters
    VectorXd::Map(&beta_c[0], beta_0.size()) = beta_0;  // stores current parameters
    VectorXd::Map(&beta_a[0], beta_0.size()) = beta_0;  // stores a refrence value for parameters
    VectorXd::Map(&beta_best[0], beta_0.size()) = beta_0;  // stores the best parameters
    VectorXd::Map(&beta_peak[0], beta_0.size()) = beta_0;  // stores the best parameters
    int iter_stop  = 0;  // tracks if the iterations should be stopped for convergence
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
    abs_max = abs_max0;
    dose_abs_max = dose_abs_max0;
    //
    for (int i = 0; i < beta_0.size(); i++) {
        beta_0[i] = a_n[i];
    }
    Cox_Term_Risk_Calc(modelform, tform, term_n, totalnum, fir, dfc, term_tot, T0, Td0, Tdd0, Te, R, Rd, Rdd, Dose, nonDose, beta_0, df0, dose_abs_max, abs_max, TTerm, nonDose_LIN, nonDose_PLIN, nonDose_LOGLIN, RdR, RddR, nthreads, KeepConstant, verbose, model_bool, gmix_theta, gmix_term);
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
    IntegerVector KeepConstant_trouble (totalnum);
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
    // variables added for log loop code
    ///
    VectorXd s_weights(1);
    int bound_val = 1;
    //
    while ((step < maxstep) && (iter_continue)) {
        step++;
        trouble = false;
        half_check = 0;
        Log_Bound(deriv_max, Lldd_mat, Lld_vec, Lstar, qchi, Ll[0], para_number, nthreads, totalnum, reqrdnum, KeepConstant, term_tot, step, dbeta, beta_0, upper, trouble, verbose, mult);
        if (trouble) {
            Calc_Change_trouble(para_number, nthreads, totalnum, dose_abs_max, lr, abs_max, Ll, Lld, Lldd, dbeta, tform, dose_abs_max, abs_max, KeepConstant_trouble);
        }
        beta_p = beta_c;  //
        beta_a = beta_c;  //
        Cox_Refresh_R_SIDES(reqrdnum, ntime, Rls1, Rls2, Rls3, Lls1, Lls2, Lls3, Strata_vals, model_bool);
        for (int ijk = 0; ijk < totalnum; ijk++) {
            if (KeepConstant[ijk] == 0) {
                if ((tform[ijk] == "lin_quad_int") || (tform[ijk] == "lin_exp_int") || (tform[ijk] == "step_int") || (tform[ijk] == "lin_int")) {  // the threshold values use different maximum deviation values
                    if (abs(dbeta[ijk]) > dose_abs_max) {
                        dbeta[ijk] = dose_abs_max * sign(dbeta[ijk]);
                    }
                }else{
                    if (abs(dbeta[ijk]) > abs_max) {
                        dbeta[ijk] = abs_max * sign(dbeta[ijk]);
                    }
                }
            } else {
                dbeta[ijk] = 0;
            }
        }
        max_change = abs(dbeta[0]);
        for (int ij = 0; ij < totalnum; ij++) {
            if (ij == para_number) {
                // prevent the parameter estimate from crossing the optimum
                // issue is beta_0[para_number] <= beta_best[para_number]
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
        // ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------//
        // The same subterm, risk, sides, and log-likelihood calculations are performed every half-step and iteration
        // ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------//
        Cox_Term_Risk_Calc(modelform, tform, term_n, totalnum, fir, dfc, term_tot, T0, Td0, Tdd0, Te, R, Rd, Rdd, Dose, nonDose, beta_0, df0, dose_abs_max, abs_max, TTerm, nonDose_LIN, nonDose_PLIN, nonDose_LOGLIN, RdR, RddR, nthreads, KeepConstant, verbose, model_bool, gmix_theta, gmix_term);
        bound_val = 1;
        Cox_Pois_Log_Loop(abs_max, model_bool, beta_0, beta_a, beta_c, bound_val, dbeta, df0, dfc, dint, Dose, dose_abs_max, dslp, fir, gmix_term, gmix_theta, half_check, halfmax, KeepConstant, limit_hit, lr, modelform, nonDose, nonDose_LIN, nonDose_LOGLIN, nonDose_PLIN, nthreads, R, Rd, Rdd, RddR, RdR, s_weights, T0, Td0, Tdd0, Te, term_n, term_tot, tform, totalnum, TTerm, verbose);
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
        Map<MatrixXd> Lldd_mat(as<Map<MatrixXd> >(Lldd_vec));
        Map<VectorXd> Lld_vec(as<Map<VectorXd> >(Lld_vecc));
        limits[1] = beta_0[para_number];
        ll_final[1] = Ll[0];
        if ((max_change < epsilon) && (deriv_max < deriv_epsilon)) {
            iter_continue = false;
            limit_converged[1] = TRUE;
        }
    }
    //
    // if (verbose >= 4) {
    //     Rcout << "C++ Note: STARTING Lower Bound" << endl;
    // }
    beta_p = beta_peak;  //
    beta_a = beta_peak;  //
    beta_c = beta_peak;  //
    for (int ij = 0; ij < totalnum; ij++) {
        beta_0[ij] = beta_a[ij];
    }
    Cox_Refresh_R_SIDES(reqrdnum, ntime, Rls1, Rls2, Rls3, Lls1, Lls2, Lls3, Strata_vals, model_bool);
    //
    Cox_Term_Risk_Calc(modelform, tform, term_n, totalnum, fir, dfc, term_tot, T0, Td0, Tdd0, Te, R, Rd, Rdd, Dose, nonDose, beta_0, df0, dose_abs_max, abs_max, TTerm, nonDose_LIN, nonDose_PLIN, nonDose_LOGLIN, RdR, RddR, nthreads, KeepConstant, verbose, model_bool, gmix_theta, gmix_term);
    //
    // -------------------------------------------------------------------------------------------
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
    abs_max = abs_max0;
    dose_abs_max = dose_abs_max0;
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
            Calc_Change_trouble(para_number, nthreads, totalnum, dose_abs_max, lr, abs_max, Ll, Lld, Lldd, dbeta, tform, dose_abs_max, abs_max, KeepConstant_trouble);
        }
        //
        beta_p = beta_c;  //
        beta_a = beta_c;  //
        Cox_Refresh_R_SIDES(reqrdnum, ntime, Rls1, Rls2, Rls3, Lls1, Lls2, Lls3, Strata_vals, model_bool);
        //
        for (int ijk = 0; ijk < totalnum; ijk++) {
            if (KeepConstant[ijk] == 0) {
                if ((tform[ijk] == "lin_quad_int") || (tform[ijk] == "lin_exp_int") || (tform[ijk] == "step_int") || (tform[ijk] == "lin_int")) {  // the threshold values use different maximum deviation values
                    if (abs(dbeta[ijk]) > dose_abs_max) {
                        dbeta[ijk] = dose_abs_max * sign(dbeta[ijk]);
                    }
                }else{
                    if (abs(dbeta[ijk]) > abs_max) {
                        dbeta[ijk] = abs_max * sign(dbeta[ijk]);
                    }
                }
            } else {
                dbeta[ijk] = 0;
            }
        }
        max_change = abs(dbeta[0]);
        for (int ij = 0; ij < totalnum; ij++) {
            if (ij == para_number) {
                // prevent the parameter estimate from crossing the optimum
                // issue is beta_0[para_number] <= beta_best[para_number]
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
        // ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------//
        // The same subterm, risk, sides, and log-likelihood calculations are performed every half-step and iteration
        // ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------//
        Cox_Term_Risk_Calc(modelform, tform, term_n, totalnum, fir, dfc, term_tot, T0, Td0, Tdd0, Te, R, Rd, Rdd, Dose, nonDose, beta_0, df0, dose_abs_max, abs_max, TTerm, nonDose_LIN, nonDose_PLIN, nonDose_LOGLIN, RdR, RddR, nthreads, KeepConstant, verbose, model_bool, gmix_theta, gmix_term);
        bound_val = 0;
        Cox_Pois_Log_Loop(abs_max, model_bool, beta_0, beta_a, beta_c, bound_val, dbeta, df0, dfc, dint, Dose, dose_abs_max, dslp, fir, gmix_term, gmix_theta, half_check, halfmax, KeepConstant, limit_hit, lr, modelform, nonDose, nonDose_LIN, nonDose_LOGLIN, nonDose_PLIN, nthreads, R, Rd, Rdd, RddR, RdR, s_weights, T0, Td0, Tdd0, Te, term_n, term_tot, tform, totalnum, TTerm, verbose);
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
        Map<MatrixXd> Lldd_mat(as<Map<MatrixXd> >(Lldd_vec));
        Map<VectorXd> Lld_vec(as<Map<VectorXd> >(Lld_vecc));
        limits[0] = beta_0[para_number];
        ll_final[0] = Ll[0];
        if ((max_change < epsilon) && (deriv_max < deriv_epsilon)) {
            iter_continue = false;
            limit_converged[0] = TRUE;
        }
    }
    //
    res_list = List::create(_["Parameter_Limits"] = wrap(limits), _["Negative_Risk_Limit_Hit"] = wrap(limit_hit), _["Likelihood_Boundary"] = wrap(ll_final), _["Likelihood_Goal"] = wrap(Lstar), _["Limit_Converged"] = wrap(limit_converged), _["Status"] = "PASSED");
    // returns a list of results
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
// [[Rcpp::export]]
List LogLik_Cox_PH_Omnibus_Log_Bound_Search(IntegerVector term_n, StringVector tform, NumericVector a_n, NumericMatrix& x_all, IntegerVector dfc, int fir, string modelform, double lr, NumericVector maxiters, int guesses, int halfmax, double epsilon, double abs_max, double dose_abs_max, double deriv_epsilon, NumericMatrix df_groups, NumericVector tu, int verbose, IntegerVector KeepConstant, int term_tot, string ties_method, int nthreads, NumericVector& Strata_vals, const VectorXd& cens_weight, List model_bool, const double gmix_theta, const IntegerVector gmix_term, const MatrixXd Lin_Sys, const VectorXd Lin_Res, double qchi, int para_number, int maxstep, double mult) {
    //
    List temp_list = List::create(_["Status"] = "TEMP");  // used as a dummy return value for code checking
    // Time durations are measured from this point on in microseconds
    //
    // df0: covariate data
    // ntime: number of event times for Cox PH
    // totalnum: number of terms used
    //
    // ------------------------------------------------------------------------- // initialize
    const Map<MatrixXd> df0(as<Map<MatrixXd> >(x_all));
    int ntime = tu.size();
    int totalnum = term_n.size();
    int reqrdnum = totalnum - sum(KeepConstant);
    // ------------------------------------------------------------------------- // initialize
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
    if (model_bool["constraint"]) {
        if (verbose >= 1) {
            Rcout << "linear constataints are not compatable with Case-Control model calculation" << endl;
        }
        temp_list = List::create(_["Status"] = "FAILED_WITH_BAD_MODEL_CONSTRAINT", _["LogLik"] = R_NaN);
        return temp_list;
    }
    //
    // cout.precision: controls the number of significant digits printed
    // nthreads: number of threads used for parallel operations
    //
    Rcout.precision(7);  // forces higher precision numbers printed to terminal
    // ------------------------------------------------------------------------- // initialize
    Map<VectorXd> beta_0(as<Map<VectorXd> >(a_n));
    VectorXd beta_max = beta_0;
    MatrixXd T0;
    MatrixXd Td0;
    MatrixXd Tdd0;
    //
    MatrixXd Te;
    MatrixXd R;
    ColXd Rd;
    ColXd Rdd;
    //
    MatrixXd Dose;
    MatrixXd nonDose;
    MatrixXd nonDose_LIN;
    MatrixXd nonDose_PLIN;
    MatrixXd nonDose_LOGLIN;
    MatrixXd TTerm;
    double dint = 0.0;  // the amount of change used to calculate derivatives in threshold paramters
    double dslp = 0.0;
    ColXd RdR;
    ColXd RddR;
    // ------------------------------------------------------------------------- // initialize
    // ---------------------------------------------
    // To Start, needs to seperate the derivative terms
    // ---------------------------------------------
    //
    Cox_Refresh_R_TERM(totalnum, reqrdnum, term_tot, dint, dslp, dose_abs_max, abs_max, df0, T0, Td0, Tdd0, Te, R, Rd, Rdd, Dose, nonDose, TTerm, nonDose_LIN, nonDose_PLIN, nonDose_LOGLIN, RdR, RddR, model_bool);
    // ------------------------------------------------------------------------- // initialize
    IntegerMatrix RiskFail;
    vector<vector<int> > RiskPairs(ntime);
    vector<vector<vector<int> > > RiskPairs_Strata(ntime, vector<vector<int>>(Strata_vals.size()));
    const Map<MatrixXd> df_m(as<Map<MatrixXd> >(df_groups));
    // ------------------------------------------------------------------------- // initialize
    if (model_bool["strata"]) {
        RiskFail = IntegerMatrix(ntime, 2*Strata_vals.size());  // vector giving the event rows
        // Creates matrices used to identify the event risk groups
        if (model_bool["cr"]) {
            Make_Groups_Strata_CR(ntime, df_m, RiskFail, RiskPairs_Strata, tu, nthreads, Strata_vals, cens_weight);
        } else {
            Make_Groups_Strata(ntime, df_m, RiskFail, RiskPairs_Strata, tu, nthreads, Strata_vals);
        }
    } else {
        RiskFail = IntegerMatrix(ntime, 2);  // vector giving the event rows
        // Creates matrices used to identify the event risk groups
        if (model_bool["cr"]) {
            Make_Groups_CR(ntime, df_m, RiskFail, RiskPairs, tu, cens_weight, nthreads);
        } else {
            Make_Groups(ntime, df_m, RiskFail, RiskPairs, tu, nthreads);
        }
    }
    // ------------------------------------------------------------------------- // initialize
    MatrixXd Rls1;
    MatrixXd Lls1;
    MatrixXd Rls2;
    MatrixXd Rls3;
    MatrixXd Lls2;
    MatrixXd Lls3;
    vector<double> Ll(reqrdnum, 0.0);  // log-likelihood values
    vector<double> Lld(reqrdnum, 0.0);  // log-likelihood derivative values
    vector<double> Lldd(pow(reqrdnum, 2), 0.0);  // the second derivative matrix has room for every combination, but only the lower triangle is calculated initially
    // ------------------------------------------------------------------------- // initialize
    Cox_Refresh_R_SIDES(reqrdnum, ntime, Rls1, Rls2, Rls3, Lls1, Lls2, Lls3, Strata_vals, model_bool);
    // the log-likelihood is calculated in parallel over the risk groups
    vector <double> Ll_comp(2, Ll[0]);  // vector to compare values
    double abs_max0 = abs_max;
    double dose_abs_max0 = dose_abs_max;
    //
    vector<double> dbeta(totalnum, 0.0);
    vector<double> dbeta_start(totalnum, 0.0);
    //
    // --------------------------
    // always starts from initial guess
    // --------------------------
    vector<double> beta_c(totalnum, 0.0);
    vector<double> beta_a(totalnum, 0.0);
    vector<double> beta_best(totalnum, 0.0);
    vector<double> beta_p(totalnum, 0.0);
    vector<double> beta_peak(totalnum, 0.0);
    VectorXd::Map(&beta_p[0], beta_0.size()) = beta_0;  // stores previous parameters
    VectorXd::Map(&beta_c[0], beta_0.size()) = beta_0;  // stores current parameters
    VectorXd::Map(&beta_a[0], beta_0.size()) = beta_0;  // stores a refrence value for parameters
    VectorXd::Map(&beta_best[0], beta_0.size()) = beta_0;  // stores the best parameters
    VectorXd::Map(&beta_peak[0], beta_0.size()) = beta_0;  // stores the peak parameters
    int iter_stop  = 0;  // tracks if the iterations should be stopped for convergence
    vector<double> beta_abs_best(totalnum, 0.0);
    //
    Cox_Refresh_R_SIDES(reqrdnum, ntime, Rls1, Rls2, Rls3, Lls1, Lls2, Lls3, Strata_vals, model_bool);
    fill(Ll.begin(), Ll.end(), 0.0);
    fill(Lld.begin(), Lld.end(), 0.0);
    fill(Lldd.begin(), Lldd.end(), 0.0);
    beta_p = beta_peak;  //
    beta_a = beta_peak;  //
    beta_c = beta_peak;  //
    abs_max = abs_max0;
    dose_abs_max = dose_abs_max0;
    //
    for (int i = 0; i < beta_0.size(); i++) {
        beta_0[i] = a_n[i];
    }
    Cox_Term_Risk_Calc(modelform, tform, term_n, totalnum, fir, dfc, term_tot, T0, Td0, Tdd0, Te, R, Rd, Rdd, Dose, nonDose, beta_0, df0, dose_abs_max, abs_max, TTerm, nonDose_LIN, nonDose_PLIN, nonDose_LOGLIN, RdR, RddR, nthreads, KeepConstant, verbose, model_bool, gmix_theta, gmix_term);
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
    IntegerVector KeepConstant_trouble (totalnum);
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
    // //
    if (verbose >= 4) {
         Rcout << "C++ Note: STARTING Upper Bound" << endl;
    }
    upper = true;
    // Now define the list of points to check
    trouble = false;
    Log_Bound(deriv_max, Lldd_mat, Lld_vec, Lstar, qchi, Lmax, para_number, nthreads, totalnum, reqrdnum, KeepConstant, term_tot, 0, dbeta_start, beta_0, upper, trouble, verbose, mult);
    NumericMatrix a_ns(guesses, totalnum);
    // now the dbeta holds the range to check
    // note that currently the guesses are not bounded
    for (int i = 0; i < guesses; i++) {
        for (int j = 0; j < totalnum; j++) {
            // use dbeta to assign a_n guess i, parameter j
            // assign such that i=guesses - 1 gives mult*dbeta
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
    // now we have the points to test
    double halves = 0;  // number of half-steps taken
    int ind0 = fir;  // used for validations
    int iteration = 0;  // iteration number
    int maxiter = 0;
    //
    iter_stop  = 0;  // tracks if the iterations should be stopped for convergence
    int iter_check = 0;  // signal to check for convergence
    //
    NumericMatrix beta_fin(a_ns.rows(), a_ns.cols());
    NumericVector LL_fin(a_ns.rows());
    //
    double Ll_iter_best = 10;
    // Variables that are used for the risk check function shared across cox, poisson, and log bound functions
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
        abs_max = abs_max0;
        dose_abs_max = dose_abs_max0;
        iter_stop = 0;
        halves = 0;
        iteration = 0;
        halves = 0;  // number of half-steps taken
        ind0 = fir;  // used for validations
        iteration = 0;  // iteration number
        Ll_iter_best = 10;
        //
        iter_stop  = 0;  // tracks if the iterations should be stopped for convergence
        iter_check = 0;  // signal to check for convergence
        //
        maxiter = maxiters[0];
        a_n = a_ns.row(guess);
        for (int i = 0; i < beta_0.size(); i++) {
            beta_0[i] = a_n[i];
        }
        Cox_Term_Risk_Calc(modelform, tform, term_n, totalnum, fir, dfc, term_tot, T0, Td0, Tdd0, Te, R, Rd, Rdd, Dose, nonDose, beta_0, df0, dose_abs_max, abs_max, TTerm, nonDose_LIN, nonDose_PLIN, nonDose_LOGLIN, RdR, RddR, nthreads, KeepConstant, verbose, model_bool, gmix_theta, gmix_term);
        if ((R.minCoeff() <= 0) || (R.hasNaN())) {
            iter_stop = 1;
            Ll[0] = 404;  // note to consider this point too far
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
            // calculates the initial change in parameter
            Calc_Change_trouble(para_number, nthreads, totalnum, dose_abs_max, lr, abs_max, Ll, Lld, Lldd, dbeta, tform, dose_abs_max, abs_max, KeepConstant_trouble);
            Intercept_Bound(nthreads, totalnum, beta_0, dbeta, dfc, df0, KeepConstant, tform);
            //
            if ((Ll_iter_best > 0) || (Ll_iter_best < Ll[ind0])) {
                Ll_iter_best = Ll[ind0];
                beta_abs_best = beta_c;
            }
            //
            halves = 0;
            while ((Ll[ind0] <= Ll_iter_best) && (halves < halfmax)) {  // repeats until half-steps maxed or an improvement
                for (int ij = 0; ij < totalnum; ij++) {
                    beta_0[ij] = beta_a[ij] + dbeta[ij];
                    beta_c[ij] = beta_0[ij];
                }
                // ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------//
                // The same subterm, risk, sides, and log-likelihood calculations are performed every half-step and iteration
                // ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------//
                Cox_Term_Risk_Calc(modelform, tform, term_n, totalnum, fir, dfc, term_tot, T0, Td0, Tdd0, Te, R, Rd, Rdd, Dose, nonDose, beta_0, df0, dose_abs_max, abs_max, TTerm, nonDose_LIN, nonDose_PLIN, nonDose_LOGLIN, RdR, RddR, nthreads, KeepConstant, verbose, model_bool, gmix_theta, gmix_term);
                Cox_Pois_Check_Continue(model_bool, beta_0, beta_best, beta_c, cens_weight, dbeta, dev, dev_temp, fir, halfmax, halves, ind0, iter_stop, KeepConstant, Ll, Ll_iter_best, Lld, Lldd, Lls1, Lls2, Lls3, Lstar, nthreads, ntime, PyrC, R, Rd, Rdd, RddR, RdR, reqrdnum, tform, RiskFail,  RiskPairs, RiskPairs_Strata, Rls1, Rls2, Rls3, Strata_vals, term_n, ties_method, totalnum, TTerm, verbose);
            }
            if (beta_best != beta_c) {  // if the risk matrices aren't the optimal values, then they must be recalculated
                // If it goes through every half step without improvement, then the maximum change needs to be decreased
                abs_max = abs_max*pow(0.5, halfmax);  // reduces the step sizes
                dose_abs_max = dose_abs_max*pow(0.5, halfmax);
                iter_check = 1;
                //
                beta_p = beta_best;  //
                beta_a = beta_best;  //
                beta_c = beta_best;  //
                for (int ij = 0; ij < totalnum; ij++) {
                    beta_0[ij] = beta_best[ij];
                }
                Cox_Term_Risk_Calc(modelform, tform, term_n, totalnum, fir, dfc, term_tot, T0, Td0, Tdd0, Te, R, Rd, Rdd, Dose, nonDose, beta_0, df0, dose_abs_max, abs_max, TTerm, nonDose_LIN, nonDose_PLIN, nonDose_LOGLIN, RdR, RddR, nthreads, KeepConstant, verbose, model_bool, gmix_theta, gmix_term);
            }
            if ((iteration % (reqrdnum)) || (iter_check == 1)) {  // checks every set number of iterations
                iter_check = 0;
                if (abs_max < epsilon/10) {  // if the maximum change is too low, then it ends
                    iter_stop = 1;
                }
            }
        }
        // -----------------------------------------------
        // Performing Full Calculation to get full second derivative matrix
        // -----------------------------------------------
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
    abs_max = abs_max0;
    dose_abs_max = dose_abs_max0;
    iter_stop = 0;
    halves = 0;
    iteration = 0;
    halves = 0;  // number of half-steps taken
    ind0 = fir;  // used for validations
    iteration = 0;  // iteration number
    //
    iter_stop  = 0;  // tracks if the iterations should be stopped for convergence
    iter_check = 0;  // signal to check for convergence
    //
    // next we need to figure out what point to start at
    int best_guess = - 1;
    for (int guess = guesses - 1; guess> - 1; guess--) {
        // we want to find the closest point at which the loglikelihood is below lstar
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
        beta_temp = wrap(beta_fin.row(best_guess));  // the first point was closest, no lower bound
    } else if (best_guess == - 1) {
        beta_temp = wrap(beta_fin.row(guesses - 1));  // the last point was closest, no upper bound
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
    // Calculates the subterm and term values
    Cox_Term_Risk_Calc(modelform, tform, term_n, totalnum, fir, dfc, term_tot, T0, Td0, Tdd0, Te, R, Rd, Rdd, Dose, nonDose, beta_0, df0, dose_abs_max, abs_max, TTerm, nonDose_LIN, nonDose_PLIN, nonDose_LOGLIN, RdR, RddR, nthreads, KeepConstant, verbose, model_bool, gmix_theta, gmix_term);
    //
    // -------------------------------------------------------------------------------------------
    // Calculates the side sum terms used
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
    // variables added for log loop code
    ///
    VectorXd s_weights(1);
    int bound_val = 1;
    //
    while ((step < maxstep) && (iter_continue)) {
        step++;
        trouble = false;
        Log_Bound(deriv_max, Lldd_mat, Lld_vec, Lstar, qchi, Ll[0], para_number, nthreads, totalnum, reqrdnum, KeepConstant, term_tot, step, dbeta, beta_0, upper, trouble, verbose, mult);
        if (trouble) {
            Calc_Change_trouble(para_number, nthreads, totalnum, dose_abs_max, lr, abs_max, Ll, Lld, Lldd, dbeta, tform, dose_abs_max, abs_max, KeepConstant_trouble);
        }
        //
        beta_p = beta_c;  //
        beta_a = beta_c;  //
        Cox_Refresh_R_SIDES(reqrdnum, ntime, Rls1, Rls2, Rls3, Lls1, Lls2, Lls3, Strata_vals, model_bool);
        //
        for (int ijk = 0; ijk < totalnum; ijk++) {
            if (KeepConstant[ijk] == 0) {
                //
                if ((tform[ijk] == "lin_quad_int") || (tform[ijk] == "lin_exp_int") || (tform[ijk] == "step_int") || (tform[ijk] == "lin_int")) {  // the threshold values use different maximum deviation values
                    if (abs(dbeta[ijk]) > dose_abs_max) {
                        dbeta[ijk] = dose_abs_max * sign(dbeta[ijk]);
                    }
                }else{
                    if (abs(dbeta[ijk]) > abs_max) {
                        dbeta[ijk] = abs_max * sign(dbeta[ijk]);
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
                // prevent the parameter estimate from crossing the optimum
                // issue is beta_0[para_number] <= beta_best[para_number]
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
        // ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------//
        // The same subterm, risk, sides, and log-likelihood calculations are performed every half-step and iteration
        // ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------//
        Cox_Term_Risk_Calc(modelform, tform, term_n, totalnum, fir, dfc, term_tot, T0, Td0, Tdd0, Te, R, Rd, Rdd, Dose, nonDose, beta_0, df0, dose_abs_max, abs_max, TTerm, nonDose_LIN, nonDose_PLIN, nonDose_LOGLIN, RdR, RddR, nthreads, KeepConstant, verbose, model_bool, gmix_theta, gmix_term);
        bound_val = 1;
        Cox_Pois_Log_Loop(abs_max, model_bool, beta_0, beta_a, beta_c, bound_val, dbeta, df0, dfc, dint, Dose, dose_abs_max, dslp, fir, gmix_term, gmix_theta, half_check, halfmax, KeepConstant, limit_hit, lr, modelform, nonDose, nonDose_LIN, nonDose_LOGLIN, nonDose_PLIN, nthreads, R, Rd, Rdd, RddR, RdR, s_weights, T0, Td0, Tdd0, Te, term_n, term_tot, tform, totalnum, TTerm, verbose);
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
        Map<MatrixXd> Lldd_mat(as<Map<MatrixXd> >(Lldd_vec));
        Map<VectorXd> Lld_vec(as<Map<VectorXd> >(Lld_vecc));
        limits[1] = beta_0[para_number];
        ll_final[1] = Ll[0];
        if ((max_change < epsilon) && (deriv_max < deriv_epsilon)) {
           iter_continue = false;
           limit_converged[1] = TRUE;
        }
    }
    // Now refresh matrices back to the maximum point
    // if (verbose >= 4) {
    //     Rcout << "C++ Note: STARTING Lower Bound" << endl;
    // }
    beta_p = beta_peak;  //
    beta_a = beta_peak;  //
    beta_c = beta_peak;  //
    abs_max = abs_max0;
    dose_abs_max = dose_abs_max0;
    for (int ij = 0; ij < totalnum; ij++) {
        beta_0[ij] = beta_a[ij];
    }
    Cox_Refresh_R_SIDES(reqrdnum, ntime, Rls1, Rls2, Rls3, Lls1, Lls2, Lls3, Strata_vals, model_bool);
    Cox_Term_Risk_Calc(modelform, tform, term_n, totalnum, fir, dfc, term_tot, T0, Td0, Tdd0, Te, R, Rd, Rdd, Dose, nonDose, beta_0, df0, dose_abs_max, abs_max, TTerm, nonDose_LIN, nonDose_PLIN, nonDose_LOGLIN, RdR, RddR, nthreads, KeepConstant, verbose, model_bool, gmix_theta, gmix_term);
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
    // now the dbeta holds the range to check
    // note that currently the guesses are not bounded
    for (int i = 0; i < guesses; i++) {
        for (int j = 0; j < totalnum; j++) {
            // use dbeta to assign a_n guess i, parameter j
            // assign such that i=guesses - 1 gives mult*dbeta
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
    // now we have the points to test
    halves = 0;  // number of half-steps taken
    ind0 = fir;  // used for validations
    iteration = 0;  // iteration number
    maxiter = 0;
    //
    iter_stop  = 0;  // tracks if the iterations should be stopped for convergence
    iter_check = 0;  // signal to check for convergence
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
        abs_max = abs_max0;
        dose_abs_max = dose_abs_max0;
        iter_stop = 0;
        halves = 0;
        iteration = 0;
        Ll_iter_best = 10;
        halves = 0;  // number of half-steps taken
        ind0 = fir;  // used for validations
        iteration = 0;  // iteration number
        iter_stop  = 0;  // tracks if the iterations should be stopped for convergence
        iter_check = 0;  // signal to check for convergence
        maxiter = maxiters[0];
        a_n = a_ns.row(guess);
        for (int i = 0; i < beta_0.size(); i++) {
            beta_0[i] = a_n[i];
        }
        Cox_Term_Risk_Calc(modelform, tform, term_n, totalnum, fir, dfc, term_tot, T0, Td0, Tdd0, Te, R, Rd, Rdd, Dose, nonDose, beta_0, df0, dose_abs_max, abs_max, TTerm, nonDose_LIN, nonDose_PLIN, nonDose_LOGLIN, RdR, RddR, nthreads, KeepConstant, verbose, model_bool, gmix_theta, gmix_term);
        if ((R.minCoeff() <= 0) || (R.hasNaN())) {
            iter_stop = 1;
            Ll[0] = 404;  // note to consider this point too far
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
            // calculates the initial change in parameter
            Calc_Change_trouble(para_number, nthreads, totalnum, dose_abs_max, lr, abs_max, Ll, Lld, Lldd, dbeta, tform, dose_abs_max, abs_max, KeepConstant_trouble);
            Intercept_Bound(nthreads, totalnum, beta_0, dbeta, dfc, df0, KeepConstant, tform);
            if ((Ll_iter_best > 0) || (Ll_iter_best < Ll[ind0])) {
                Ll_iter_best = Ll[ind0];
                beta_abs_best = beta_c;
            }
            halves = 0;
            while ((Ll[ind0] <= Ll_iter_best) && (halves < halfmax)) {  // repeats until half-steps maxed or an improvement
                for (int ij = 0; ij < totalnum; ij++) {
                    beta_0[ij] = beta_a[ij] + dbeta[ij];
                    beta_c[ij] = beta_0[ij];
                }
                // ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------//
                // The same subterm, risk, sides, and log-likelihood calculations are performed every half-step and iteration
                // ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------//
                Cox_Term_Risk_Calc(modelform, tform, term_n, totalnum, fir, dfc, term_tot, T0, Td0, Tdd0, Te, R, Rd, Rdd, Dose, nonDose, beta_0, df0, dose_abs_max, abs_max, TTerm, nonDose_LIN, nonDose_PLIN, nonDose_LOGLIN, RdR, RddR, nthreads, KeepConstant, verbose, model_bool, gmix_theta, gmix_term);
                Cox_Pois_Check_Continue(model_bool, beta_0, beta_best, beta_c, cens_weight, dbeta, dev, dev_temp, fir, halfmax, halves, ind0, iter_stop, KeepConstant, Ll, Ll_iter_best, Lld, Lldd, Lls1, Lls2, Lls3, Lstar, nthreads, ntime, PyrC, R, Rd, Rdd, RddR, RdR, reqrdnum, tform, RiskFail,  RiskPairs, RiskPairs_Strata, Rls1, Rls2, Rls3, Strata_vals, term_n, ties_method, totalnum, TTerm, verbose);
            }
            if (beta_best != beta_c) {  // if the risk matrices aren't the optimal values, then they must be recalculated
                // If it goes through every half step without improvement, then the maximum change needs to be decreased
                abs_max = abs_max*pow(0.5, halfmax);  // reduces the step sizes
                dose_abs_max = dose_abs_max*pow(0.5, halfmax);
                iter_check = 1;
                //
                beta_p = beta_best;  //
                beta_a = beta_best;  //
                beta_c = beta_best;  //
                for (int ij = 0; ij < totalnum; ij++) {
                    beta_0[ij] = beta_best[ij];
                }
                Cox_Term_Risk_Calc(modelform, tform, term_n, totalnum, fir, dfc, term_tot, T0, Td0, Tdd0, Te, R, Rd, Rdd, Dose, nonDose, beta_0, df0, dose_abs_max, abs_max, TTerm, nonDose_LIN, nonDose_PLIN, nonDose_LOGLIN, RdR, RddR, nthreads, KeepConstant, verbose, model_bool, gmix_theta, gmix_term);
            }
            if ((iteration % (reqrdnum)) || (iter_check == 1)) {  // checks every set number of iterations
                iter_check = 0;
                if (abs_max < epsilon/10) {  // if the maximum change is too low, then it ends
                    iter_stop = 1;
                }
            }
        }
        // -----------------------------------------------
        // Performing Full Calculation to get full second derivative matrix
        // -----------------------------------------------
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
    abs_max = abs_max0;
    dose_abs_max = dose_abs_max0;
    iter_stop = 0;
    halves = 0;
    iteration = 0;
    halves = 0;  // number of half-steps taken
    ind0 = fir;  // used for validations
    iteration = 0;  // iteration number
    iter_stop  = 0;  // tracks if the iterations should be stopped for convergence
    iter_check = 0;  // signal to check for convergence
    // next we need to figure out what point to start at
    best_guess = - 1;
    for (int guess = guesses - 1; guess> - 1; guess--) {
        // we want to find the closest point at which the loglikelihood is below lstar
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
        beta_temp = wrap(beta_fin.row(best_guess));  // the first point was closest, no lower bound
    } else if (best_guess == - 1) {
        beta_temp = wrap(beta_fin.row(guesses - 1));  // the last point was closest, no upper bound
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
    // Calculates the subterm and term values
    Cox_Term_Risk_Calc(modelform, tform, term_n, totalnum, fir, dfc, term_tot, T0, Td0, Tdd0, Te, R, Rd, Rdd, Dose, nonDose, beta_0, df0, dose_abs_max, abs_max, TTerm, nonDose_LIN, nonDose_PLIN, nonDose_LOGLIN, RdR, RddR, nthreads, KeepConstant, verbose, model_bool, gmix_theta, gmix_term);
    // Calculates the side sum terms used
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
        step++;
        trouble = false;
        Log_Bound(deriv_max, Lldd_mat, Lld_vec, Lstar, qchi, Ll[0], para_number, nthreads, totalnum, reqrdnum, KeepConstant, term_tot, step, dbeta, beta_0, upper, trouble, verbose, mult);
        if (trouble) {
            Calc_Change_trouble(para_number, nthreads, totalnum, dose_abs_max, lr, abs_max, Ll, Lld, Lldd, dbeta, tform, dose_abs_max, abs_max, KeepConstant_trouble);
        }
        beta_p = beta_c;  //
        beta_a = beta_c;  //
        Cox_Refresh_R_SIDES(reqrdnum, ntime, Rls1, Rls2, Rls3, Lls1, Lls2, Lls3, Strata_vals, model_bool);
        //
        for (int ijk = 0; ijk < totalnum; ijk++) {
            if (KeepConstant[ijk] == 0) {
                //
                if ((tform[ijk] == "lin_quad_int") || (tform[ijk] == "lin_exp_int") || (tform[ijk] == "step_int") || (tform[ijk] == "lin_int")) {  // the threshold values use different maximum deviation values
                    if (abs(dbeta[ijk]) > dose_abs_max) {
                        dbeta[ijk] = dose_abs_max * sign(dbeta[ijk]);
                    }
                }else{
                    if (abs(dbeta[ijk]) > abs_max) {
                        dbeta[ijk] = abs_max * sign(dbeta[ijk]);
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
                // prevent the parameter estimate from crossing the optimum
                // issue is beta_0[para_number] <= beta_best[para_number]
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
        // ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------//
        // The same subterm, risk, sides, and log-likelihood calculations are performed every half-step and iteration
        // ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------//
        Cox_Term_Risk_Calc(modelform, tform, term_n, totalnum, fir, dfc, term_tot, T0, Td0, Tdd0, Te, R, Rd, Rdd, Dose, nonDose, beta_0, df0, dose_abs_max, abs_max, TTerm, nonDose_LIN, nonDose_PLIN, nonDose_LOGLIN, RdR, RddR, nthreads, KeepConstant, verbose, model_bool, gmix_theta, gmix_term);
        bound_val = 0;
        Cox_Pois_Log_Loop(abs_max, model_bool, beta_0, beta_a, beta_c, bound_val, dbeta, df0, dfc, dint, Dose, dose_abs_max, dslp, fir, gmix_term, gmix_theta, half_check, halfmax, KeepConstant, limit_hit, lr, modelform, nonDose, nonDose_LIN, nonDose_LOGLIN, nonDose_PLIN, nthreads, R, Rd, Rdd, RddR, RdR, s_weights, T0, Td0, Tdd0, Te, term_n, term_tot, tform, totalnum, TTerm, verbose);
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
        Map<MatrixXd> Lldd_mat(as<Map<MatrixXd> >(Lldd_vec));
        Map<VectorXd> Lld_vec(as<Map<VectorXd> >(Lld_vecc));
        limits[0] = beta_0[para_number];
        ll_final[0] = Ll[0];
        if ((max_change < epsilon) && (deriv_max < deriv_epsilon)) {
            iter_continue = false;
            limit_converged[0] = TRUE;
        }
    }
    //
    res_list = List::create(_["Parameter_Limits"] = wrap(limits), _["Negative_Risk_Limit_Hit"] = wrap(limit_hit), _["Likelihood_Boundary"] = wrap(ll_final), _["Likelihood_Goal"] = wrap(Lstar), _["Limit_Converged"] = wrap(limit_converged), _["Status"] = "PASSED");
    // returns a list of results
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
// [[Rcpp::export]]
List LogLik_Poisson_Omnibus_Log_Bound(const MatrixXd& PyrC, const MatrixXd& dfs, IntegerVector term_n, StringVector tform, NumericVector a_n, NumericMatrix& x_all, IntegerVector dfc, int fir, string modelform, double lr, NumericVector maxiters, int guesses, int halfmax, double epsilon, double abs_max, double dose_abs_max, double deriv_epsilon, int verbose, IntegerVector KeepConstant, int term_tot, int nthreads, List model_bool, const double gmix_theta, const IntegerVector gmix_term, const MatrixXd Lin_Sys, const VectorXd Lin_Res, double qchi, int para_number, int maxstep, double mult) {
    //
    List temp_list = List::create(_["Status"] = "TEMP");  // used as a dummy return value for code checking
    //
    // Time durations are measured from this point on in microseconds
    //
    // df0: covariate data
    // totalnum: number of terms used
    //
    // ------------------------------------------------------------------------- // initialize
    const Map<MatrixXd> df0(as<Map<MatrixXd> >(x_all));
    const int mat_row = df0.rows();
    int totalnum = term_n.size();
    int reqrdnum = totalnum - sum(KeepConstant);
    // ------------------------------------------------------------------------- // initialize
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
    // cout.precision: controls the number of significant digits printed
    // nthreads: number of threads used for parallel operations
    //
    Rcout.precision(7);  // forces higher precision numbers printed to terminal
    // ------------------------------------------------------------------------- // initialize
    Map<VectorXd> beta_0(as<Map<VectorXd> >(a_n));
    MatrixXd T0;
    MatrixXd Td0;
    MatrixXd Tdd0;
    MatrixXd Te;
    MatrixXd R;
    ColXd Rd;
    ColXd Rdd;
    MatrixXd Dose;
    MatrixXd nonDose;
    MatrixXd nonDose_LIN;
    MatrixXd nonDose_PLIN;
    MatrixXd nonDose_LOGLIN;
    MatrixXd TTerm;
    double dint = 0.0;  // the amount of change used to calculate derivatives in threshold paramters
    double dslp = 0.0;
    ColXd RdR;
    ColXd RddR;
    // ------------------------------------------------------------------------- // initialize
    // ---------------------------------------------
    // To Start, needs to seperate the derivative terms
    // ---------------------------------------------
    //
    T0 = MatrixXd::Zero(mat_row, totalnum);  // preallocates matrix for Term column
    Te = MatrixXd::Zero(mat_row, 1);  // preallocates matrix for column terms used for temporary storage
    R = MatrixXd::Zero(mat_row, 1);  // preallocates matrix for Risks
    //
    Dose = MatrixXd::Constant(mat_row, term_tot, 0.0);  // matrix of the total dose term values
    nonDose = MatrixXd::Constant(mat_row, term_tot, 1.0);  // matrix of the total non-dose term values
    nonDose_LIN = MatrixXd::Constant(mat_row, term_tot, 0.0);  // matrix of Linear subterm values
    nonDose_PLIN = MatrixXd::Constant(mat_row, term_tot, 1.0);  // matrix of Loglinear subterm values
    nonDose_LOGLIN = MatrixXd::Constant(mat_row, term_tot, 1.0);  // matrix of Product linear subterm values
    TTerm = MatrixXd::Zero(mat_row, term_tot);  // matrix of term values
    if (model_bool["single"]) {
    } else {
        Td0 = MatrixXd::Zero(mat_row, reqrdnum);  // preallocates matrix for Term derivative columns
        Tdd0 = MatrixXd::Zero(mat_row, reqrdnum*(reqrdnum + 1)/2);  // preallocates matrix for Term second derivative columns
        //
        Rd = MatrixXd::Zero(mat_row, reqrdnum);  // preallocates matrix for Risk derivatives
        Rdd = MatrixXd::Zero(mat_row, reqrdnum*(reqrdnum + 1)/2);  // preallocates matrix for Risk second derivatives
        //
        dint = dose_abs_max;  // the amount of change used to calculate derivatives in threshold paramters
        dslp = abs_max;
        RdR = MatrixXd::Zero(mat_row, reqrdnum);  // preallocates matrix for Risk to derivative ratios
        RddR = MatrixXd::Zero(mat_row, reqrdnum*(reqrdnum + 1)/2);  // preallocates matrix for Risk to second derivative ratios
    }
    VectorXd s_weights;
    if (model_bool["strata"]) {
        s_weights = VectorXd::Zero(mat_row);
        Gen_Strat_Weight(modelform, dfs, PyrC, s_weights, nthreads, tform, term_n, term_tot);
    }
    // ------------------------------------------------------------------------- // initialize
    vector<double> Ll(reqrdnum, 0.0);  // log-likelihood values
    vector<double> Lld(reqrdnum, 0.0);  // log-likelihood derivative values
    vector<double> Lldd(pow(reqrdnum, 2), 0.0);  // the second derivative matrix has room for every combination, but only the lower triangle is calculated initially
    MatrixXd dev_temp = MatrixXd::Zero(PyrC.rows(), 2);
    double dev = 0;
    // the log-likelihood is calculated in parallel over the risk groups
    vector <double> Ll_comp(2, Ll[0]);  // vector to compare values
    double abs_max0 = abs_max;
    double dose_abs_max0 = dose_abs_max;
    //
    vector<double> dbeta(totalnum, 0.0);
    //
    // --------------------------
    // always starts from initial guess
    // --------------------------
    vector<double> beta_c(totalnum, 0.0);
    vector<double> beta_a(totalnum, 0.0);
    vector<double> beta_best(totalnum, 0.0);
    vector<double> beta_peak(totalnum, 0.0);
    vector<double> beta_p(totalnum, 0.0);
    VectorXd::Map(&beta_p[0], beta_0.size()) = beta_0;  // stores previous parameters
    VectorXd::Map(&beta_c[0], beta_0.size()) = beta_0;  // stores current parameters
    VectorXd::Map(&beta_a[0], beta_0.size()) = beta_0;  // stores a refrence value for parameters
    VectorXd::Map(&beta_best[0], beta_0.size()) = beta_0;  // stores the best parameters
    VectorXd::Map(&beta_peak[0], beta_0.size()) = beta_0;  // stores the best parameters
    int iter_stop  = 0;  // tracks if the iterations should be stopped for convergence
    vector<double> beta_abs_best(totalnum, 0.0);
    //
    fill(Ll.begin(), Ll.end(), 0.0);
    fill(Lld.begin(), Lld.end(), 0.0);
    fill(Lldd.begin(), Lldd.end(), 0.0);
    beta_p = beta_peak;  //
    beta_a = beta_peak;  //
    beta_c = beta_peak;  //
    abs_max = abs_max0;
    dose_abs_max = dose_abs_max0;
    //
    for (int i = 0; i < beta_0.size(); i++) {
        beta_0[i] = a_n[i];
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
    IntegerVector KeepConstant_trouble (totalnum);
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
    // variables added for log loop code
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
            Calc_Change_trouble(para_number, nthreads, totalnum, dose_abs_max, lr, abs_max, Ll, Lld, Lldd, dbeta, tform, dose_abs_max, abs_max, KeepConstant_trouble);
        }
        //
        beta_p = beta_c;  //
        beta_a = beta_c;  //

        //
        for (int ijk = 0; ijk < totalnum; ijk++) {
            if (KeepConstant[ijk] == 0) {
                //
                if ((tform[ijk] == "lin_quad_int") || (tform[ijk] == "lin_exp_int") || (tform[ijk] == "step_int") || (tform[ijk] == "lin_int")) {  // the threshold values use different maximum deviation values
                    if (abs(dbeta[ijk]) > dose_abs_max) {
                        dbeta[ijk] = dose_abs_max * sign(dbeta[ijk]);
                    }
                }else{
                    if (abs(dbeta[ijk]) > abs_max) {
                        dbeta[ijk] = abs_max * sign(dbeta[ijk]);
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
                // prevent the parameter estimate from crossing the optimum
                // issue is beta_0[para_number] <= beta_best[para_number]
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
        // ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------//
        // The same subterm, risk, sides, and log-likelihood calculations are performed every half-step and iteration
        // ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------//

        Pois_Term_Risk_Calc(modelform, tform, term_n, totalnum, fir, dfc, term_tot, T0, Td0, Tdd0, Te, R, Rd, Rdd, Dose, nonDose, beta_0, df0, dint, dslp, TTerm, nonDose_LIN, nonDose_PLIN, nonDose_LOGLIN, RdR, RddR, s_weights, nthreads, KeepConstant, verbose, model_bool, gmix_theta, gmix_term);
        bound_val = 1;
        Cox_Pois_Log_Loop(abs_max, model_bool, beta_0, beta_a, beta_c, bound_val, dbeta, df0, dfc, dint, Dose, dose_abs_max, dslp, fir, gmix_term, gmix_theta, half_check, halfmax, KeepConstant, limit_hit, lr, modelform, nonDose, nonDose_LIN, nonDose_LOGLIN, nonDose_PLIN, nthreads, R, Rd, Rdd, RddR, RdR, s_weights, T0, Td0, Tdd0, Te, term_n, term_tot, tform, totalnum, TTerm, verbose);
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
        Map<MatrixXd> Lldd_mat(as<Map<MatrixXd> >(Lldd_vec));
        Map<VectorXd> Lld_vec(as<Map<VectorXd> >(Lld_vecc));
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
    abs_max = abs_max0;
    dose_abs_max = dose_abs_max0;
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
            Calc_Change_trouble(para_number, nthreads, totalnum, dose_abs_max, lr, abs_max, Ll, Lld, Lldd, dbeta, tform, dose_abs_max, abs_max, KeepConstant_trouble);
        }
        beta_p = beta_c;  //
        beta_a = beta_c;  //
        for (int ijk = 0; ijk < totalnum; ijk++) {
            if (KeepConstant[ijk] == 0) {
                //
                if ((tform[ijk] == "lin_quad_int") || (tform[ijk] == "lin_exp_int") || (tform[ijk] == "step_int") || (tform[ijk] == "lin_int")) {  // the threshold values use different maximum deviation values
                    if (abs(dbeta[ijk]) > dose_abs_max) {
                        dbeta[ijk] = dose_abs_max * sign(dbeta[ijk]);
                    }
                }else{
                    if (abs(dbeta[ijk]) > abs_max) {
                        dbeta[ijk] = abs_max * sign(dbeta[ijk]);
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
                // prevent the parameter estimate from crossing the optimum 
                // issue is beta_0[para_number] <= beta_best[para_number]
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
        // ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------//
        // The same subterm, risk, sides, and log-likelihood calculations are performed every half-step and iteration
        // ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------//
        Pois_Term_Risk_Calc(modelform, tform, term_n, totalnum, fir, dfc, term_tot, T0, Td0, Tdd0, Te, R, Rd, Rdd, Dose, nonDose, beta_0, df0, dint, dslp, TTerm, nonDose_LIN, nonDose_PLIN, nonDose_LOGLIN, RdR, RddR, s_weights, nthreads, KeepConstant, verbose, model_bool, gmix_theta, gmix_term);
        bound_val = 0;
        Cox_Pois_Log_Loop(abs_max, model_bool, beta_0, beta_a, beta_c, bound_val, dbeta, df0, dfc, dint, Dose, dose_abs_max, dslp, fir, gmix_term, gmix_theta, half_check, halfmax, KeepConstant, limit_hit, lr, modelform, nonDose, nonDose_LIN, nonDose_LOGLIN, nonDose_PLIN, nthreads, R, Rd, Rdd, RddR, RdR, s_weights, T0, Td0, Tdd0, Te, term_n, term_tot, tform, totalnum, TTerm, verbose);
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
        Map<MatrixXd> Lldd_mat(as<Map<MatrixXd> >(Lldd_vec));
        Map<VectorXd> Lld_vec(as<Map<VectorXd> >(Lld_vecc));
        limits[0] = beta_0[para_number];
        ll_final[0] = Ll[0];
        if ((max_change < epsilon) && (deriv_max < deriv_epsilon)) {
            iter_continue = false;
            limit_converged[0] = TRUE;
        }
    }
    res_list = List::create(_["Parameter_Limits"] = wrap(limits), _["Negative_Risk_Limit_Hit"] = wrap(limit_hit), _["Likelihood_Boundary"] = wrap(ll_final), _["Likelihood_Goal"] = wrap(Lstar), _["Limit_Converged"] = wrap(limit_converged), _["Status"] = "PASSED");
    // returns a list of results
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
// [[Rcpp::export]]
List LogLik_Poisson_Omnibus_Log_Bound_Search(const MatrixXd& PyrC, const MatrixXd& dfs, IntegerVector term_n, StringVector tform, NumericVector a_n, NumericMatrix& x_all, IntegerVector dfc, int fir, string modelform, double lr, NumericVector maxiters, int guesses, int halfmax, double epsilon, double abs_max, double dose_abs_max, double deriv_epsilon, int verbose, IntegerVector KeepConstant, int term_tot, int nthreads, List model_bool, const double gmix_theta, const IntegerVector gmix_term, const MatrixXd Lin_Sys, const VectorXd Lin_Res, double qchi, int para_number, int maxstep, double mult) {
    //
    List temp_list = List::create(_["Status"] = "TEMP");  // used as a dummy return value for code checking
    //
    // Time durations are measured from this point on in microseconds
    //
    // df0: covariate data
    // ntime: number of event times for Cox PH
    // totalnum: number of terms used
    //
    // ------------------------------------------------------------------------- // initialize
    const Map<MatrixXd> df0(as<Map<MatrixXd> >(x_all));
    const int mat_row = df0.rows();
    int totalnum = term_n.size();
    int reqrdnum = totalnum - sum(KeepConstant);
    // ------------------------------------------------------------------------- // initialize
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
    // cout.precision: controls the number of significant digits printed
    // nthreads: number of threads used for parallel operations
    //
    Rcout.precision(7);  // forces higher precision numbers printed to terminal
    // ------------------------------------------------------------------------- // initialize
    Map<VectorXd> beta_0(as<Map<VectorXd> >(a_n));
    VectorXd beta_max = beta_0;
    MatrixXd T0;
    MatrixXd Td0;
    MatrixXd Tdd0;
    //
    MatrixXd Te;
    MatrixXd R;
    ColXd Rd;
    ColXd Rdd;
    //
    MatrixXd Dose;
    MatrixXd nonDose;
    MatrixXd nonDose_LIN;
    MatrixXd nonDose_PLIN;
    MatrixXd nonDose_LOGLIN;
    MatrixXd TTerm;
    double dint = 0.0;  // the amount of change used to calculate derivatives in threshold paramters
    double dslp = 0.0;
    ColXd RdR;
    ColXd RddR;
    // ------------------------------------------------------------------------- // initialize
    // ---------------------------------------------
    // To Start, needs to seperate the derivative terms
    // ---------------------------------------------
    //
    T0 = MatrixXd::Zero(mat_row, totalnum);  // preallocates matrix for Term column
    Te = MatrixXd::Zero(mat_row, 1);  // preallocates matrix for column terms used for temporary storage
    R = MatrixXd::Zero(mat_row, 1);  // preallocates matrix for Risks
    //
    Dose = MatrixXd::Constant(mat_row, term_tot, 0.0);  // matrix of the total dose term values
    nonDose = MatrixXd::Constant(mat_row, term_tot, 1.0);  // matrix of the total non-dose term values
    nonDose_LIN = MatrixXd::Constant(mat_row, term_tot, 0.0);  // matrix of Linear subterm values
    nonDose_PLIN = MatrixXd::Constant(mat_row, term_tot, 1.0);  // matrix of Loglinear subterm values
    nonDose_LOGLIN = MatrixXd::Constant(mat_row, term_tot, 1.0);  // matrix of Product linear subterm values
    TTerm = MatrixXd::Zero(mat_row, term_tot);  // matrix of term values
    if (model_bool["single"]) {
    } else {
        Td0 = MatrixXd::Zero(mat_row, reqrdnum);  // preallocates matrix for Term derivative columns
        Tdd0 = MatrixXd::Zero(mat_row, reqrdnum*(reqrdnum + 1)/2);  // preallocates matrix for Term second derivative columns
        Rd = MatrixXd::Zero(mat_row, reqrdnum);  // preallocates matrix for Risk derivatives
        Rdd = MatrixXd::Zero(mat_row, reqrdnum*(reqrdnum + 1)/2);  // preallocates matrix for Risk second derivatives
        dint = dose_abs_max;  // the amount of change used to calculate derivatives in threshold paramters
        dslp = abs_max;
        RdR = MatrixXd::Zero(mat_row, reqrdnum);  // preallocates matrix for Risk to derivative ratios
        RddR = MatrixXd::Zero(mat_row, reqrdnum*(reqrdnum + 1)/2);  // preallocates matrix for Risk to second derivative ratios
    }
    VectorXd s_weights;
    if (model_bool["strata"]) {
        s_weights = VectorXd::Zero(mat_row);
        Gen_Strat_Weight(modelform, dfs, PyrC, s_weights, nthreads, tform, term_n, term_tot);
    }
    // ------------------------------------------------------------------------- // initialize
    MatrixXd Rls1;
    MatrixXd Lls1;
    MatrixXd Rls2;
    MatrixXd Rls3;
    MatrixXd Lls2;
    MatrixXd Lls3;
    vector<double> Ll(reqrdnum, 0.0);  // log-likelihood values
    vector<double> Lld(reqrdnum, 0.0);  // log-likelihood derivative values
    vector<double> Lldd(pow(reqrdnum, 2), 0.0);  // the second derivative matrix has room for every combination, but only the lower triangle is calculated initially
    MatrixXd dev_temp = MatrixXd::Zero(PyrC.rows(), 2);
    double dev = 0;
    // ------------------------------------------------------------------------- // initialize
    // the log-likelihood is calculated in parallel over the risk groups
    vector <double> Ll_comp(2, Ll[0]);  // vector to compare values
    double abs_max0 = abs_max;
    double dose_abs_max0 = dose_abs_max;
    //
    vector<double> dbeta(totalnum, 0.0);
    vector<double> dbeta_start(totalnum, 0.0);
    //
    // --------------------------
    // always starts from initial guess
    // --------------------------
    vector<double> beta_c(totalnum, 0.0);
    vector<double> beta_a(totalnum, 0.0);
    vector<double> beta_best(totalnum, 0.0);
    vector<double> beta_p(totalnum, 0.0);
    vector<double> beta_peak(totalnum, 0.0);
    VectorXd::Map(&beta_p[0], beta_0.size()) = beta_0;  // stores previous parameters
    VectorXd::Map(&beta_c[0], beta_0.size()) = beta_0;  // stores current parameters
    VectorXd::Map(&beta_a[0], beta_0.size()) = beta_0;  // stores a refrence value for parameters
    VectorXd::Map(&beta_best[0], beta_0.size()) = beta_0;  // stores the best parameters
    VectorXd::Map(&beta_peak[0], beta_0.size()) = beta_0;  // stores the peak parameters
    int iter_stop  = 0;  // tracks if the iterations should be stopped for convergence
    vector<double> beta_abs_best(totalnum, 0.0);
    fill(Ll.begin(), Ll.end(), 0.0);
    fill(Lld.begin(), Lld.end(), 0.0);
    fill(Lldd.begin(), Lldd.end(), 0.0);
    beta_p = beta_peak;  //
    beta_a = beta_peak;  //
    beta_c = beta_peak;  //
    abs_max = abs_max0;
    dose_abs_max = dose_abs_max0;
    for (int i = 0; i < beta_0.size(); i++) {
        beta_0[i] = a_n[i];
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
    IntegerVector KeepConstant_trouble (totalnum);
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
    // if (verbose >= 4) {
    //     Rcout << "C++ Note: STARTING Upper Bound" << endl;
    // }
    upper = true;
    // Now define the list of points to check
    trouble = false;
    Log_Bound(deriv_max, Lldd_mat, Lld_vec, Lstar, qchi, Lmax, para_number, nthreads, totalnum, reqrdnum, KeepConstant, term_tot, 0, dbeta_start, beta_0, upper, trouble, verbose, mult);
    NumericMatrix a_ns(guesses, totalnum);
    // now the dbeta holds the range to check
    // note that currently the guesses are not bounded
    for (int i = 0; i < guesses; i++) {
        for (int j = 0; j < totalnum; j++) {
            // use dbeta to assign a_n guess i, parameter j
            // assign such that i=guesses - 1 gives mult*dbeta
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
    // now we have the points to test
    double halves = 0;  // number of half-steps taken
    int ind0 = fir;  // used for validations
    int iteration = 0;  // iteration number
    int maxiter = 0;
    iter_stop  = 0;  // tracks if the iterations should be stopped for convergence
    int iter_check = 0;  // signal to check for convergence
    NumericMatrix beta_fin(a_ns.rows(), a_ns.cols());
    NumericVector LL_fin(a_ns.rows());
    double Ll_iter_best = 10;
    // Variables that are used for the risk check function shared across cox, poisson, and log bound functions
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
        abs_max = abs_max0;
        dose_abs_max = dose_abs_max0;
        iter_stop = 0;
        halves = 0;
        iteration = 0;
        halves = 0;  // number of half-steps taken
        ind0 = fir;  // used for validations
        iteration = 0;  // iteration number
        Ll_iter_best = 10;
        iter_stop  = 0;  // tracks if the iterations should be stopped for convergence
        iter_check = 0;  // signal to check for convergence
        maxiter = maxiters[0];
        a_n = a_ns.row(guess);
        for (int i = 0; i < beta_0.size(); i++) {
            beta_0[i] = a_n[i];
        }
        Pois_Term_Risk_Calc(modelform, tform, term_n, totalnum, fir, dfc, term_tot, T0, Td0, Tdd0, Te, R, Rd, Rdd, Dose, nonDose, beta_0, df0, dint, dslp, TTerm, nonDose_LIN, nonDose_PLIN, nonDose_LOGLIN, RdR, RddR, s_weights, nthreads, KeepConstant, verbose, model_bool, gmix_theta, gmix_term);
        if ((R.minCoeff() <= 0) || (R.hasNaN())) {
            iter_stop = 1;
            Ll[0] = 404;  // note to consider this point too far
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
            // calculates the initial change in parameter
            Calc_Change_trouble(para_number, nthreads, totalnum, dose_abs_max, lr, abs_max, Ll, Lld, Lldd, dbeta, tform, dose_abs_max, abs_max, KeepConstant_trouble);
            Intercept_Bound(nthreads, totalnum, beta_0, dbeta, dfc, df0, KeepConstant, tform);
            //
            if ((Ll_iter_best > 0) || (Ll_iter_best < Ll[ind0])) {
                Ll_iter_best = Ll[ind0];
                beta_abs_best = beta_c;
            }
            halves = 0;
            while ((Ll[ind0] <= Ll_iter_best) && (halves < halfmax)) {  // repeats until half-steps maxed or an improvement
                for (int ij = 0; ij < totalnum; ij++) {
                    beta_0[ij] = beta_a[ij] + dbeta[ij];
                    beta_c[ij] = beta_0[ij];
                }
                // ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------//
                // The same subterm, risk, sides, and log-likelihood calculations are performed every half-step and iteration
                // ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------//
                Pois_Term_Risk_Calc(modelform, tform, term_n, totalnum, fir, dfc, term_tot, T0, Td0, Tdd0, Te, R, Rd, Rdd, Dose, nonDose, beta_0, df0, dint, dslp, TTerm, nonDose_LIN, nonDose_PLIN, nonDose_LOGLIN, RdR, RddR, s_weights, nthreads, KeepConstant, verbose, model_bool, gmix_theta, gmix_term);
                Cox_Pois_Check_Continue(model_bool, beta_0, beta_best, beta_c, cens_weight, dbeta, dev, dev_temp, fir, halfmax, halves, ind0, iter_stop, KeepConstant, Ll, Ll_iter_best, Lld, Lldd, Lls1, Lls2, Lls3, Lstar, nthreads, ntime, PyrC, R, Rd, Rdd, RddR, RdR, reqrdnum, tform, RiskFail,  RiskPairs, RiskPairs_Strata, Rls1, Rls2, Rls3, Strata_vals, term_n, ties_method, totalnum, TTerm, verbose);
            }
            if (beta_best != beta_c) {  // if the risk matrices aren't the optimal values, then they must be recalculated
                // If it goes through every half step without improvement, then the maximum change needs to be decreased
                abs_max = abs_max*pow(0.5, halfmax);  // reduces the step sizes
                dose_abs_max = dose_abs_max*pow(0.5, halfmax);
                iter_check = 1;
                //
                beta_p = beta_best;  //
                beta_a = beta_best;  //
                beta_c = beta_best;  //
                for (int ij = 0; ij < totalnum; ij++) {
                    beta_0[ij] = beta_best[ij];
                }
                Pois_Term_Risk_Calc(modelform, tform, term_n, totalnum, fir, dfc, term_tot, T0, Td0, Tdd0, Te, R, Rd, Rdd, Dose, nonDose, beta_0, df0, dint, dslp, TTerm, nonDose_LIN, nonDose_PLIN, nonDose_LOGLIN, RdR, RddR, s_weights, nthreads, KeepConstant, verbose, model_bool, gmix_theta, gmix_term);
            }
            if ((iteration % (reqrdnum)) || (iter_check == 1)) {  // checks every set number of iterations
                iter_check = 0;
                if (abs_max < epsilon/10) {  // if the maximum change is too low, then it ends
                    iter_stop = 1;
                }
            }
        }
        // -----------------------------------------------
        // Performing Full Calculation to get full second derivative matrix
        // -----------------------------------------------
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
    abs_max = abs_max0;
    dose_abs_max = dose_abs_max0;
    iter_stop = 0;
    halves = 0;
    iteration = 0;
    halves = 0;  // number of half-steps taken
    ind0 = fir;  // used for validations
    iteration = 0;  // iteration number
    //
    iter_stop  = 0;  // tracks if the iterations should be stopped for convergence
    iter_check = 0;  // signal to check for convergence
    //
    // next we need to figure out what point to start at
    int best_guess = - 1;
    for (int guess = guesses - 1; guess> - 1; guess--) {
        // we want to find the closest point at which the loglikelihood is below lstar
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
        beta_temp = wrap(beta_fin.row(best_guess));  // the first point was closest, no lower bound
    } else if (best_guess == - 1) {
        beta_temp = wrap(beta_fin.row(guesses - 1));  // the last point was closest, no upper bound
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
    // Calculates the subterm and term values
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
    // variables added for log loop code
    ///
    int bound_val = 1;
    //
    while ((step < maxstep) && (iter_continue)) {
        step++;
        trouble = false;
        Log_Bound(deriv_max, Lldd_mat, Lld_vec, Lstar, qchi, Ll[0], para_number, nthreads, totalnum, reqrdnum, KeepConstant, term_tot, step, dbeta, beta_0, upper, trouble, verbose, mult);
        if (trouble) {
            Calc_Change_trouble(para_number, nthreads, totalnum, dose_abs_max, lr, abs_max, Ll, Lld, Lldd, dbeta, tform, dose_abs_max, abs_max, KeepConstant_trouble);
        }
        beta_p = beta_c;  //
        beta_a = beta_c;  //
        for (int ijk = 0; ijk < totalnum; ijk++) {
            if (KeepConstant[ijk] == 0) {
                //
                if ((tform[ijk] == "lin_quad_int") || (tform[ijk] == "lin_exp_int") || (tform[ijk] == "step_int") || (tform[ijk] == "lin_int")) {  // the threshold values use different maximum deviation values
                    if (abs(dbeta[ijk]) > dose_abs_max) {
                        dbeta[ijk] = dose_abs_max * sign(dbeta[ijk]);
                    }
                }else{
                    if (abs(dbeta[ijk]) > abs_max) {
                        dbeta[ijk] = abs_max * sign(dbeta[ijk]);
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
               // prevent the parameter estimate from crossing the optimum
               // issue is beta_0[para_number] <= beta_best[para_number]
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
        // ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------//
        // The same subterm, risk, sides, and log-likelihood calculations are performed every half-step and iteration
        // ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------//
        Pois_Term_Risk_Calc(modelform, tform, term_n, totalnum, fir, dfc, term_tot, T0, Td0, Tdd0, Te, R, Rd, Rdd, Dose, nonDose, beta_0, df0, dint, dslp, TTerm, nonDose_LIN, nonDose_PLIN, nonDose_LOGLIN, RdR, RddR, s_weights, nthreads, KeepConstant, verbose, model_bool, gmix_theta, gmix_term);
        bound_val = 1;
        Cox_Pois_Log_Loop(abs_max, model_bool, beta_0, beta_a, beta_c, bound_val, dbeta, df0, dfc, dint, Dose, dose_abs_max, dslp, fir, gmix_term, gmix_theta, half_check, halfmax, KeepConstant, limit_hit, lr, modelform, nonDose, nonDose_LIN, nonDose_LOGLIN, nonDose_PLIN, nthreads, R, Rd, Rdd, RddR, RdR, s_weights, T0, Td0, Tdd0, Te, term_n, term_tot, tform, totalnum, TTerm, verbose);
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
        Map<MatrixXd> Lldd_mat(as<Map<MatrixXd> >(Lldd_vec));
        Map<VectorXd> Lld_vec(as<Map<VectorXd> >(Lld_vecc));
        limits[1] = beta_0[para_number];
        ll_final[1] = Ll[0];
        if ((max_change < epsilon) && (deriv_max < deriv_epsilon)) {
           iter_continue = false;
           limit_converged[1] = TRUE;
        }
    }
    // Now refresh matrices back to the maximum point
    // if (verbose >= 4) {
    //     Rcout << "C++ Note: STARTING Lower Bound" << endl;
    // }
    beta_p = beta_peak;  //
    beta_a = beta_peak;  //
    beta_c = beta_peak;  //
    abs_max = abs_max0;
    dose_abs_max = dose_abs_max0;
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
    // now the dbeta holds the range to check
    // note that currently the guesses are not bounded
    for (int i = 0; i < guesses; i++) {
        for (int j = 0; j < totalnum; j++) {
            // use dbeta to assign a_n guess i, parameter j
            // assign such that i=guesses - 1 gives mult*dbeta
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
    // now we have the points to test
    halves = 0;  // number of half-steps taken
    ind0 = fir;  // used for validations
    iteration = 0;  // iteration number
    maxiter = 0;
    iter_stop  = 0;  // tracks if the iterations should be stopped for convergence
    iter_check = 0;  // signal to check for convergence
    Ll_iter_best = 10;
    for (int guess = 0; guess <guesses; guess++) {
        fill(Ll.begin(), Ll.end(), 0.0);
        fill(Lld.begin(), Lld.end(), 0.0);
        fill(Lldd.begin(), Lldd.end(), 0.0);
        beta_p = beta_peak;  //
        beta_a = beta_peak;  //
        beta_c = beta_peak;  //
        abs_max = abs_max0;
        dose_abs_max = dose_abs_max0;
        iter_stop = 0;
        halves = 0;
        iteration = 0;
        halves = 0;  // number of half-steps taken
        ind0 = fir;  // used for validations
        iteration = 0;  // iteration number
        Ll_iter_best = 10;
        iter_stop  = 0;  // tracks if the iterations should be stopped for convergence
        iter_check = 0;  // signal to check for convergence
        maxiter = maxiters[0];
        a_n = a_ns.row(guess);
        for (int i = 0; i < beta_0.size(); i++) {
            beta_0[i] = a_n[i];
        }
        Pois_Term_Risk_Calc(modelform, tform, term_n, totalnum, fir, dfc, term_tot, T0, Td0, Tdd0, Te, R, Rd, Rdd, Dose, nonDose, beta_0, df0, dint, dslp, TTerm, nonDose_LIN, nonDose_PLIN, nonDose_LOGLIN, RdR, RddR, s_weights, nthreads, KeepConstant, verbose, model_bool, gmix_theta, gmix_term);
        if ((R.minCoeff() <= 0) || (R.hasNaN())) {
            iter_stop = 1;
            Ll[0] = 404;  // note to consider this point too far
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
            // calculates the initial change in parameter
            Calc_Change_trouble(para_number, nthreads, totalnum, dose_abs_max, lr, abs_max, Ll, Lld, Lldd, dbeta, tform, dose_abs_max, abs_max, KeepConstant_trouble);
            Intercept_Bound(nthreads, totalnum, beta_0, dbeta, dfc, df0, KeepConstant, tform);
            if ((Ll_iter_best > 0) || (Ll_iter_best < Ll[ind0])) {
                Ll_iter_best = Ll[ind0];
                beta_abs_best = beta_c;
            }
            halves = 0;
            while ((Ll[ind0] <= Ll_iter_best) && (halves < halfmax)) {  // repeats until half-steps maxed or an improvement
                for (int ij = 0; ij < totalnum; ij++) {
                    beta_0[ij] = beta_a[ij] + dbeta[ij];
                    beta_c[ij] = beta_0[ij];
                }
                // ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------//
                // The same subterm, risk, sides, and log-likelihood calculations are performed every half-step and iteration
                // ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------//
                Pois_Term_Risk_Calc(modelform, tform, term_n, totalnum, fir, dfc, term_tot, T0, Td0, Tdd0, Te, R, Rd, Rdd, Dose, nonDose, beta_0, df0, dint, dslp, TTerm, nonDose_LIN, nonDose_PLIN, nonDose_LOGLIN, RdR, RddR, s_weights, nthreads, KeepConstant, verbose, model_bool, gmix_theta, gmix_term);
                Cox_Pois_Check_Continue(model_bool, beta_0, beta_best, beta_c, cens_weight, dbeta, dev, dev_temp, fir, halfmax, halves, ind0, iter_stop, KeepConstant, Ll, Ll_iter_best, Lld, Lldd, Lls1, Lls2, Lls3, Lstar, nthreads, ntime, PyrC, R, Rd, Rdd, RddR, RdR, reqrdnum, tform, RiskFail,  RiskPairs, RiskPairs_Strata, Rls1, Rls2, Rls3, Strata_vals, term_n, ties_method, totalnum, TTerm, verbose);
            }
            if (beta_best != beta_c) {  // if the risk matrices aren't the optimal values, then they must be recalculated
                // If it goes through every half step without improvement, then the maximum change needs to be decreased
                abs_max = abs_max*pow(0.5, halfmax);  // reduces the step sizes
                dose_abs_max = dose_abs_max*pow(0.5, halfmax);
                iter_check = 1;
                beta_p = beta_best;  //
                beta_a = beta_best;  //
                beta_c = beta_best;  //
                for (int ij = 0; ij < totalnum; ij++) {
                    beta_0[ij] = beta_best[ij];
                }
                Pois_Term_Risk_Calc(modelform, tform, term_n, totalnum, fir, dfc, term_tot, T0, Td0, Tdd0, Te, R, Rd, Rdd, Dose, nonDose, beta_0, df0, dint, dslp, TTerm, nonDose_LIN, nonDose_PLIN, nonDose_LOGLIN, RdR, RddR, s_weights, nthreads, KeepConstant, verbose, model_bool, gmix_theta, gmix_term);
            }
            if ((iteration % (reqrdnum)) || (iter_check == 1)) {  // checks every set number of iterations
                iter_check = 0;
                if (abs_max < epsilon/10) {  // if the maximum change is too low, then it ends
                    iter_stop = 1;
                }
            }
        }
        // -----------------------------------------------
        // Performing Full Calculation to get full second derivative matrix
        // -----------------------------------------------
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
    abs_max = abs_max0;
    dose_abs_max = dose_abs_max0;
    iter_stop = 0;
    halves = 0;
    iteration = 0;
    halves = 0;  // number of half-steps taken
    ind0 = fir;  // used for validations
    iteration = 0;  // iteration number
    iter_stop  = 0;  // tracks if the iterations should be stopped for convergence
    iter_check = 0;  // signal to check for convergence
    // next we need to figure out what point to start at
    best_guess = - 1;
    for (int guess = guesses - 1; guess> - 1; guess--) {
        // we want to find the closest point at which the loglikelihood is below lstar
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
        beta_temp = wrap(beta_fin.row(best_guess));  // the first point was closest, no lower bound
    } else if (best_guess == - 1) {
        beta_temp = wrap(beta_fin.row(guesses - 1));  // the last point was closest, no upper bound
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
    // Calculates the subterm and term values
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
            Calc_Change_trouble(para_number, nthreads, totalnum, dose_abs_max, lr, abs_max, Ll, Lld, Lldd, dbeta, tform, dose_abs_max, abs_max, KeepConstant_trouble);
        }
        //
        beta_p = beta_c;  //
        beta_a = beta_c;  //
        //
        for (int ijk = 0; ijk < totalnum; ijk++) {
            if (KeepConstant[ijk] == 0) {
                //
                if ((tform[ijk] == "lin_quad_int") || (tform[ijk] == "lin_exp_int") || (tform[ijk] == "step_int") || (tform[ijk] == "lin_int")) {  // the threshold values use different maximum deviation values
                    if (abs(dbeta[ijk]) > dose_abs_max) {
                        dbeta[ijk] = dose_abs_max * sign(dbeta[ijk]);
                    }
                }else{
                    if (abs(dbeta[ijk]) > abs_max) {
                        dbeta[ijk] = abs_max * sign(dbeta[ijk]);
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
                // prevent the parameter estimate from crossing the optimum
                // issue is beta_0[para_number] <= beta_best[para_number]
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
        // ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------//
        // The same subterm, risk, sides, and log-likelihood calculations are performed every half-step and iteration
        // ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------//
        Pois_Term_Risk_Calc(modelform, tform, term_n, totalnum, fir, dfc, term_tot, T0, Td0, Tdd0, Te, R, Rd, Rdd, Dose, nonDose, beta_0, df0, dint, dslp, TTerm, nonDose_LIN, nonDose_PLIN, nonDose_LOGLIN, RdR, RddR, s_weights, nthreads, KeepConstant, verbose, model_bool, gmix_theta, gmix_term);
        bound_val = 0;
        Cox_Pois_Log_Loop(abs_max, model_bool, beta_0, beta_a, beta_c, bound_val, dbeta, df0, dfc, dint, Dose, dose_abs_max, dslp, fir, gmix_term, gmix_theta, half_check, halfmax, KeepConstant, limit_hit, lr, modelform, nonDose, nonDose_LIN, nonDose_LOGLIN, nonDose_PLIN, nthreads, R, Rd, Rdd, RddR, RdR, s_weights, T0, Td0, Tdd0, Te, term_n, term_tot, tform, totalnum, TTerm, verbose);
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
        Map<MatrixXd> Lldd_mat(as<Map<MatrixXd> >(Lldd_vec));
        Map<VectorXd> Lld_vec(as<Map<VectorXd> >(Lld_vecc));
        limits[0] = beta_0[para_number];
        ll_final[0] = Ll[0];
        if ((max_change < epsilon) && (deriv_max < deriv_epsilon)) {
           iter_continue = false;
           limit_converged[0] = TRUE;
        }
    }
    //
    res_list = List::create(_["Parameter_Limits"] = wrap(limits), _["Negative_Risk_Limit_Hit"] = wrap(limit_hit), _["Likelihood_Boundary"] = wrap(ll_final), _["Likelihood_Goal"] = wrap(Lstar), _["Limit_Converged"] = wrap(limit_converged), _["Status"] = "PASSED");
    // returns a list of results
    return res_list;
}

//' Primary Cox PH regression with multiple distributed dose columns and optional combinations of null, stratification, competing risks, multiplicative log-linear model, and no derivative calculation.
//'
//' \code{LogLik_Cox_PH_Multidose_Omnibus_Serial} Performs the calls to calculation functions, Structures the Cox PH regression, With verbose option prints out time stamps and intermediate sums of terms and derivatives
//'
//' @inheritParams CPP_template
//'
//' @return List of final results: Log-likelihood of optimum, standard error, and convergence for each realization
//' @noRd
//'
// [[Rcpp::export]]
List LogLik_Cox_PH_Multidose_Omnibus_Serial(IntegerVector term_n, StringVector tform, NumericVector a_n, NumericMatrix& x_all, NumericMatrix& dose_all, IntegerMatrix dose_cols, IntegerVector dose_index, IntegerVector dfc, int fir, string modelform, double lr, List optim_para, int maxiter, int halfmax, double epsilon, double abs_max, double dose_abs_max, double deriv_epsilon, NumericMatrix df_groups, NumericVector tu, int double_step, int verbose, IntegerVector KeepConstant, int term_tot, string ties_method, int nthreads, NumericVector& Strata_vals, const VectorXd& cens_weight, List model_bool, const double gmix_theta, const IntegerVector gmix_term, const MatrixXd Lin_Sys, const VectorXd Lin_Res) {
    //
    List temp_list = List::create(_["Status"] = "TEMP");  // used as a dummy return value for code checking
    // Time durations are measured from this point on in microseconds
    time_point<system_clock> start_point, end_point;
    start_point = system_clock::now();
    auto start = time_point_cast<microseconds>(start_point).time_since_epoch().count();
    end_point = system_clock::now();
    auto ending = time_point_cast<microseconds>(end_point).time_since_epoch().count();  // the time duration is tracked
    //
    // df0: covariate data
    // ntime: number of event times for Cox PH
    // totalnum: number of terms used
    //
    // ------------------------------------------------------------------------- // initialize
    Map<MatrixXd> df0(as<Map<MatrixXd> >(x_all));
    const Map<MatrixXd> df1(as<Map<MatrixXd> >(dose_all));
    int ntime = tu.size();
    int totalnum = term_n.size();
    int reqrdnum = totalnum - sum(KeepConstant);
    // ------------------------------------------------------------------------- // initialize
    if (model_bool["null"]) {
        if (verbose >= 1) {
            Rcout << "null model is not compatable with multi-realization method" << endl;
        }
        temp_list = List::create(_["Status"] = "FAILED_WITH_BAD_MODEL_NULL", _["LogLik"] = R_NaN);
        return temp_list;
    }
    if (model_bool["single"]) {
        if (verbose >= 1) {
            Rcout << "non-derivative model calculation is not compatable with multi-realization method" << endl;
        }
        temp_list = List::create(_["Status"] = "FAILED_WITH_BAD_MODEL_SINGLE", _["LogLik"] = R_NaN);
        return temp_list;
    }
    // if (model_bool["gradient"]) {
    //     if (verbose >= 1) {
    //         Rcout << "gradient descent model calculation is not CURRENTLY compatable with multi-realization method" << endl;
    //     }
    //     temp_list = List::create(_["Status"] = "FAILED_WITH_BAD_MODEL_GRADIENT", _["LogLik"] = R_NaN);
    //     return temp_list;
    // }
    //
    // cout.precision: controls the number of significant digits printed
    // nthreads: number of threads used for parallel operations
    //
    Rcout.precision(7);  // forces higher precision numbers printed to terminal
    //
    // Lld_worst: stores the highest magnitude log-likelihood derivative
    double Lld_worst = 0.0;  // stores derivative value used to determine if every parameter is near convergence
    //
    // ---------------------------------------------
    // To Start, needs to seperate the derivative terms
    // ---------------------------------------------
    // ------------------------------------------------------------------------- // initialize
    Map<VectorXd> beta_0(as<Map<VectorXd> >(a_n));
    MatrixXd T0;
    MatrixXd Td0;
    MatrixXd Tdd0;
    MatrixXd Te;
    MatrixXd R;
    ColXd Rd;
    ColXd Rdd;
    MatrixXd Dose;
    MatrixXd nonDose;
    MatrixXd nonDose_LIN;
    MatrixXd nonDose_PLIN;
    MatrixXd nonDose_LOGLIN;
    MatrixXd TTerm;
    double dint = 0.0;  // the amount of change used to calculate derivatives in threshold paramters
    double dslp = 0.0;
    ColXd RdR;
    ColXd RddR;
    // ------------------------------------------------------------------------- // initialize
    // ---------------------------------------------
    // To Start, needs to seperate the derivative terms
    // ---------------------------------------------
    // Manually define the full memory reserve
    const int mat_row = df0.rows();
    T0 = MatrixXd::Zero(mat_row, totalnum);  // preallocates matrix for Term column
    if (model_bool["basic"]) {
        R = MatrixXd::Zero(mat_row, 1);  // preallocates matrix for Risks
        Rd = MatrixXd::Zero(mat_row, reqrdnum);  // preallocates matrix for Risk derivatives
        Rdd = MatrixXd::Zero(mat_row, reqrdnum*(reqrdnum + 1)/2);  // preallocates matrix for Risk second derivatives
        RdR = MatrixXd::Zero(mat_row, reqrdnum);  // preallocates matrix for Risk to derivative ratios
        TTerm = MatrixXd::Zero(mat_row, 1);  // matrix of term values
    } else if (model_bool["linear_err"]) {
        R = MatrixXd::Zero(mat_row, 1);  // preallocates matrix for Risks
        Rd = MatrixXd::Zero(mat_row, reqrdnum);  // preallocates matrix for Risk derivatives
        Rdd = MatrixXd::Zero(mat_row, reqrdnum*(reqrdnum + 1)/2);  // preallocates matrix for Risk second derivatives
        nonDose_PLIN = MatrixXd::Constant(mat_row, 1, 1.0);  // matrix of Loglinear subterm values
        nonDose_LOGLIN = MatrixXd::Constant(mat_row, 1, 1.0);  // matrix of Product linear subterm values
        TTerm = MatrixXd::Zero(mat_row, 1);  // matrix of term values
        RdR = MatrixXd::Zero(mat_row, reqrdnum);  // preallocates matrix for Risk to derivative ratios
        RddR = MatrixXd::Zero(mat_row, reqrdnum*(reqrdnum + 1)/2);  // preallocates matrix for Risk to second derivative ratios
    } else {
        Td0 = MatrixXd::Zero(mat_row, reqrdnum);  // preallocates matrix for Term derivative columns
        Tdd0 = MatrixXd::Zero(mat_row, reqrdnum*(reqrdnum + 1)/2);  // preallocates matrix for Term second derivative columns
        Te = MatrixXd::Zero(mat_row, 1);  // preallocates matrix for column terms used for temporary storage
        R = MatrixXd::Zero(mat_row, 1);  // preallocates matrix for Risks
        Rd = MatrixXd::Zero(mat_row, reqrdnum);  // preallocates matrix for Risk derivatives
        Rdd = MatrixXd::Zero(mat_row, reqrdnum*(reqrdnum + 1)/2);  // preallocates matrix for Risk second derivatives
        Dose = MatrixXd::Constant(mat_row, term_tot, 0.0);  // matrix of the total dose term values
        nonDose = MatrixXd::Constant(mat_row, term_tot, 1.0);  // matrix of the total non-dose term values
        nonDose_LIN = MatrixXd::Constant(mat_row, term_tot, 0.0);  // matrix of Linear subterm values
        nonDose_PLIN = MatrixXd::Constant(mat_row, term_tot, 1.0);  // matrix of Loglinear subterm values
        nonDose_LOGLIN = MatrixXd::Constant(mat_row, term_tot, 1.0);  // matrix of Product linear subterm values
        TTerm = MatrixXd::Zero(mat_row, term_tot);  // matrix of term values
        dint = dose_abs_max;  // the amount of change used to calculate derivatives in threshold paramters
        dslp = abs_max;
        RdR = MatrixXd::Zero(mat_row, reqrdnum);  // preallocates matrix for Risk to derivative ratios
        RddR = MatrixXd::Zero(mat_row, reqrdnum*(reqrdnum + 1)/2);  // preallocates matrix for Risk to second derivative ratios
    }
    //
    // Cox_Refresh_R_TERM(totalnum, reqrdnum, term_tot, dint, dslp, dose_abs_max, abs_max, df0, T0, Td0, Tdd0, Te, R, Rd, Rdd, Dose, nonDose, TTerm, nonDose_LIN, nonDose_PLIN, nonDose_LOGLIN, RdR, RddR, model_bool);
    IntegerMatrix RiskFail;
    vector<vector<int> > RiskPairs(ntime);
    vector<vector<vector<int> > > RiskPairs_Strata(ntime, vector<vector<int>>(Strata_vals.size()));
    const Map<MatrixXd> df_m(as<Map<MatrixXd> >(df_groups));
    if (model_bool["strata"]) {
        RiskFail = IntegerMatrix(ntime, 2*Strata_vals.size());  // vector giving the event rows
        // Creates matrices used to identify the event risk groups
        if (model_bool["cr"]) {
            Make_Groups_Strata_CR(ntime, df_m, RiskFail, RiskPairs_Strata, tu, nthreads, Strata_vals, cens_weight);
        } else {
            Make_Groups_Strata(ntime, df_m, RiskFail, RiskPairs_Strata, tu, nthreads, Strata_vals);
        }
    } else {
        RiskFail = IntegerMatrix(ntime, 2);  // vector giving the event rows
        // Creates matrices used to identify the event risk groups
        if (model_bool["cr"]) {
            Make_Groups_CR(ntime, df_m, RiskFail, RiskPairs, tu, cens_weight, nthreads);
        } else {
            Make_Groups(ntime, df_m, RiskFail, RiskPairs, tu, nthreads);
        }
    }
    // ------------------------------------------------------------------------- // initialize
    MatrixXd Rls1;
    MatrixXd Lls1;
    MatrixXd Rls2;
    MatrixXd Rls3;
    MatrixXd Lls2;
    MatrixXd Lls3;
    vector<double> Ll(reqrdnum, 0.0);  // log-likelihood values
    vector<double> Lld(reqrdnum, 0.0);  // log-likelihood derivative values
    vector<double> Lldd(pow(reqrdnum, 2), 0.0);  // the second derivative matrix has room for every combination, but only the lower triangle is calculated initially
    // ------------------------------------------------------------------------- // initialize
    // Manually reserve the amount of memory needed
    if (model_bool["strata"]) {
        Rls1 = MatrixXd::Zero(ntime, Strata_vals.size());  // precomputes a series of sums used frequently in the log-liklihood calculations
        Lls1 = MatrixXd::Zero(ntime, Strata_vals.size());  // the log-likelihood calculation has a Right and Left sum used
        Rls2 = MatrixXd::Zero(ntime, reqrdnum*Strata_vals.size());  // many are repeated due to the same risk groups and derivatives being used at mulitple points
        Lls2 = MatrixXd::Zero(ntime, reqrdnum*Strata_vals.size());
        Rls3 = MatrixXd::Zero(ntime, reqrdnum*(reqrdnum + 1)/2*Strata_vals.size());  // sum and its derivatives are precomputed
        Lls3 = MatrixXd::Zero(ntime, reqrdnum*(reqrdnum + 1)/2*Strata_vals.size());
    } else {
        Rls1 = MatrixXd::Zero(ntime, 1);  // precomputes a series of sums used frequently in the log-liklihood calculations
        Lls1 = MatrixXd::Zero(ntime, 1);  // the log-likelihood calculation has a Right and Left sum used
        Rls2 = MatrixXd::Zero(ntime, reqrdnum);  // many are repeated due to the same risk groups and derivatives being used at mulitple points
        Lls2 = MatrixXd::Zero(ntime, reqrdnum);
        Rls3 = MatrixXd::Zero(ntime, reqrdnum*(reqrdnum + 1)/2);  // sum and its derivatives are precomputed
        Lls3 = MatrixXd::Zero(ntime, reqrdnum*(reqrdnum + 1)/2);
    }
    // Cox_Refresh_R_SIDES(reqrdnum, ntime, Rls1, Rls2, Rls3, Lls1, Lls2, Lls3, Strata_vals, model_bool);
    // the log-likelihood is calculated in parallel over the risk groups
    vector <double> Ll_comp(2, Ll[0]);  // vector to compare values
    double abs_max0 = abs_max;
    double dose_abs_max0 = dose_abs_max;
    vector<double> dbeta(totalnum, 0.0);
    NumericVector m_g_store(reqrdnum);
    NumericVector v_beta_store(reqrdnum);
    // --------------------------
    // always starts from initial guess
    // --------------------------
    vector<double> beta_c(totalnum, 0.0);
    vector<double> beta_a(totalnum, 0.0);
    vector<double> beta_best(totalnum, 0.0);
    vector<double> beta_p(totalnum, 0.0);
    VectorXd::Map(&beta_p[0], beta_0.size()) = beta_0;  // stores previous parameters
    VectorXd::Map(&beta_c[0], beta_0.size()) = beta_0;  // stores current parameters
    VectorXd::Map(&beta_a[0], beta_0.size()) = beta_0;  // stores a refrence value for parameters
    VectorXd::Map(&beta_best[0], beta_0.size()) = beta_0;  // stores the best parameters
    double halves = 0;  // number of half-steps taken
    int ind0 = fir;  // used for validations
    int iteration = 0;  // iteration number
    int iter_stop  = 0;  // tracks if the iterations should be stopped for convergence
    int iter_check = 0;  // signal to check for convergence
    int guesses = dose_cols.cols();
    NumericMatrix beta_fin(dose_cols.cols(), totalnum);
    NumericVector LL_fin(dose_cols.cols());
    NumericVector AIC_fin(dose_cols.cols());
    LogicalVector conv_fin(dose_cols.cols());
    NumericMatrix std_fin(dose_cols.cols(), totalnum);
    IntegerVector dfc_0(dfc.length());
    for (int i = 0;  i<dfc.length(); i++) {
        dfc_0[i] = dfc[i];
    }
    //
    double Ll_abs_best = 10;
    vector<double> beta_abs_best(totalnum, 0.0);
    // Variables that are used for the risk check function shared across cox, poisson, and log bound functions
    double dev = 0.0;
    MatrixXd dev_temp = MatrixXd::Zero(1, 1);
    double Lstar = 0.0;
    MatrixXd PyrC = MatrixXd::Zero(1, 1);
    bool convgd = FALSE;
    //
//    end_point = system_clock::now();
//    ending = time_point_cast<microseconds>(end_point).time_since_epoch().count();
//    Rcout << "Prep step: " << ending - start << endl;
    List out_list;
    for (int guess = 0; guess <guesses; guess++) {
//        start_point = system_clock::now();
//        start = time_point_cast<microseconds>(start_point).time_since_epoch().count();
        Cox_Refresh_R_SIDES(reqrdnum, ntime, Rls1, Rls2, Rls3, Lls1, Lls2, Lls3, Strata_vals, model_bool);
        fill(Ll.begin(), Ll.end(), 0.0);
        fill(Lld.begin(), Lld.end(), 0.0);
        if (!model_bool["gradient"]){
            fill(Lldd.begin(), Lldd.end(), 0.0);
        }
        if (model_bool["gradient"]) {
            m_g_store.fill(0);
            v_beta_store.fill(0);
        }
        beta_p = beta_best;  //
        beta_a = beta_best;  //
        beta_c = beta_best;  //
        abs_max = abs_max0;
        dose_abs_max = dose_abs_max0;
        iter_stop = 0;
        halves = 0;
        convgd = FALSE;
        iteration = 0;
        halves = 0;  // number of half-steps taken
        ind0 = fir;  // used for validations
        iteration = 0;  // iteration number
        Ll_abs_best = 10;
        iter_stop  = 0;  // tracks if the iterations should be stopped for convergence
        iter_check = 0;  // signal to check for convergence
        for (int i = 0; i < beta_0.size(); i++) {
            beta_0[i] = beta_best[i];
        }
        for (int i = 0; i<dose_index.size(); i++) {
            for (int j = 0; j < totalnum; j++) {
                if (dfc_0[j] == dose_index[i]) {
                    df0.col(dfc[j]-1) = df1.col(dose_cols(i, guess)-1);
//                    dfc[j] = dose_cols(i, guess);
                }
            }
        }
//        end_point = system_clock::now();
//        ending = time_point_cast<microseconds>(end_point).time_since_epoch().count();
//        Rcout << "Step Prep " << guess << " step: " << ending - start << endl;
//        start_point = system_clock::now();
//        start = time_point_cast<microseconds>(start_point).time_since_epoch().count();
        Cox_Term_Risk_Calc(modelform, tform, term_n, totalnum, fir, dfc, term_tot, T0, Td0, Tdd0, Te, R, Rd, Rdd, Dose, nonDose, beta_0, df0, dose_abs_max, abs_max, TTerm, nonDose_LIN, nonDose_PLIN, nonDose_LOGLIN, RdR, RddR, nthreads, KeepConstant, verbose, model_bool, gmix_theta, gmix_term);
        int check_loop = 0;
        while (((R.minCoeff() <= 0) || (R.hasNaN())) && (check_loop < 3)) { // If a linear column changes, then a previous optimum can give a negative risk
            check_loop = check_loop + 1; // the parameters in a negative term are reduced to zero, which should generally scale the term value closer to zero as well
            for (int ijk = 0; ijk < totalnum; ijk++) {
                int tij = term_n[ijk];
                if (TTerm.col(tij).minCoeff()<=0) {
                    beta_0[ijk] = beta_0[ijk] / 2.0;
                } else if (isinf(TTerm.col(tij).maxCoeff())) {
                    beta_0[ijk] = beta_0[ijk] / 2.0;
                } else if (isnan(TTerm.col(tij).minCoeff())) {
                    beta_0[ijk] = beta_0[ijk] / 2.0;
                }
            }
            Cox_Term_Risk_Calc(modelform, tform, term_n, totalnum, fir, dfc, term_tot, T0, Td0, Tdd0, Te, R, Rd, Rdd, Dose, nonDose, beta_0, df0, dose_abs_max, abs_max, TTerm, nonDose_LIN, nonDose_PLIN, nonDose_LOGLIN, RdR, RddR, nthreads, KeepConstant, verbose, model_bool, gmix_theta, gmix_term);
        }
        if ((R.minCoeff() <= 0) || (R.hasNaN())) {
            if (verbose >= 1) {
                Rcout << "C++ Error: A non-positive risk was detected: " << R.minCoeff() << endl;
                Rcout << "C++ Warning: final failing values ";
                for (int ijk = 0; ijk < totalnum; ijk++) {
                    Rcout << beta_0[ijk] << " ";
                }
                Rcout << " " << endl;
            }
            temp_list = List::create(_["beta_0"] = wrap(beta_0), _["Deviation"] = R_NaN, _["Status"] = "FAILED_WITH_NEGATIVE_RISK", _["LogLik"] = R_NaN);
            return temp_list;
        }
        Cox_Side_LL_Calc(reqrdnum, ntime, tform, RiskFail,  RiskPairs, RiskPairs_Strata, totalnum, fir, R, Rd, Rdd, Rls1, Rls2, Rls3, Lls1, Lls2, Lls3, cens_weight, Strata_vals, beta_0, RdR, RddR, Ll, Lld, Lldd, nthreads, KeepConstant, ties_method, verbose, model_bool, iter_stop);
        //
        for (int i = 0; i < beta_0.size(); i++) {
            beta_c[i] = beta_0[i];
        }
//        end_point = system_clock::now();
//        ending = time_point_cast<microseconds>(end_point).time_since_epoch().count();
//        Rcout << "First Guess " << guess << " step: " << ending - start << endl;
//        start_point = system_clock::now();
//        start = time_point_cast<microseconds>(start_point).time_since_epoch().count();
        // out_list = Cox_Full_Run(reqrdnum, ntime, tform, RiskFail,  totalnum, fir, R, Rd, Rdd, Rls1, Rls2, Rls3, Lls1, Lls2, Lls3, cens_weight, Strata_vals, beta_0, RdR, RddR, Ll, Lld, Lldd, nthreads, KeepConstant, ties_method, verbose, model_bool, iter_stop, term_tot, dint, dslp, dose_abs_max, abs_max, df0, T0, Td0, Tdd0, Te, Dose, nonDose, TTerm, nonDose_LIN, nonDose_PLIN, nonDose_LOGLIN, modelform, gmix_theta, gmix_term, convgd, lr, optim_para, maxiter, double_step, Lin_Sys, Lin_Res, term_n, dfc, halfmax, epsilon, deriv_epsilon);
        while ((iteration < maxiter) && (iter_stop == 0)) {
            iteration++;
            beta_p = beta_c;  //
            beta_a = beta_c;  //
            beta_best = beta_c;  //
            // calculates the initial change in parameter
            if (model_bool["basic"]) {
                Calc_Change_Basic(double_step, nthreads, totalnum, lr, abs_max, Ll, Lld, Lldd, dbeta, KeepConstant);
            } else if (model_bool["gradient"]) {
                Calc_Change_Gradient(nthreads, model_bool, totalnum, optim_para, iteration, abs_max, Lld, m_g_store, v_beta_store, dbeta, KeepConstant);
            } else {
                if (model_bool["constraint"]) {
                    Calc_Change_Cons(Lin_Sys, Lin_Res, beta_0, nthreads, totalnum, dose_abs_max, lr, abs_max, Ll, Lld, Lldd, dbeta, tform, dose_abs_max, abs_max, KeepConstant);
                } else {
                    Calc_Change(double_step, nthreads, totalnum, dose_abs_max, lr, abs_max, Ll, Lld, Lldd, dbeta, tform, dose_abs_max, abs_max, KeepConstant);
                }
                Intercept_Bound(nthreads, totalnum, beta_0, dbeta, dfc, df0, KeepConstant, tform);
            }
            //
            if ((Ll_abs_best > 0) || (Ll_abs_best < Ll[ind0])) {
                Ll_abs_best = Ll[ind0];
                beta_abs_best = beta_c;
            }
            //
            if (model_bool["gradient"]){
                //
                for (int ij = 0; ij < totalnum; ij++) {
                    beta_0[ij] = beta_a[ij] + dbeta[ij];
                    beta_c[ij] = beta_0[ij];
                }
                Cox_Term_Risk_Calc(modelform, tform, term_n, totalnum, fir, dfc, term_tot, T0, Td0, Tdd0, Te, R, Rd, Rdd, Dose, nonDose, beta_0, df0, dose_abs_max, abs_max, TTerm, nonDose_LIN, nonDose_PLIN, nonDose_LOGLIN, RdR, RddR, nthreads, KeepConstant, verbose, model_bool, gmix_theta, gmix_term);
                Cox_Pois_Check_Continue(model_bool, beta_0, beta_best, beta_c, cens_weight, dbeta, dev, dev_temp, fir, halfmax, halves, ind0, iter_stop, KeepConstant, Ll, Ll_abs_best, Lld, Lldd, Lls1, Lls2, Lls3, Lstar, nthreads, ntime, PyrC, R, Rd, Rdd, RddR, RdR, reqrdnum, tform, RiskFail,  RiskPairs, RiskPairs_Strata, Rls1, Rls2, Rls3, Strata_vals, term_n, ties_method, totalnum, TTerm, verbose);
            } else {
                halves = 0;
                while ((Ll[ind0] <= Ll_abs_best) && (halves < halfmax)) {  // repeats until half-steps maxed or an improvement
                    for (int ij = 0; ij < totalnum; ij++) {
                        beta_0[ij] = beta_a[ij] + dbeta[ij];
                        beta_c[ij] = beta_0[ij];
                    }
                    // ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------//
                    // The same subterm, risk, sides, and log-likelihood calculations are performed every half-step and iteration
                    // ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------//
                    Cox_Term_Risk_Calc(modelform, tform, term_n, totalnum, fir, dfc, term_tot, T0, Td0, Tdd0, Te, R, Rd, Rdd, Dose, nonDose, beta_0, df0, dose_abs_max, abs_max, TTerm, nonDose_LIN, nonDose_PLIN, nonDose_LOGLIN, RdR, RddR, nthreads, KeepConstant, verbose, model_bool, gmix_theta, gmix_term);
                    Cox_Pois_Check_Continue(model_bool, beta_0, beta_best, beta_c, cens_weight, dbeta, dev, dev_temp, fir, halfmax, halves, ind0, iter_stop, KeepConstant, Ll, Ll_abs_best, Lld, Lldd, Lls1, Lls2, Lls3, Lstar, nthreads, ntime, PyrC, R, Rd, Rdd, RddR, RdR, reqrdnum, tform, RiskFail,  RiskPairs, RiskPairs_Strata, Rls1, Rls2, Rls3, Strata_vals, term_n, ties_method, totalnum, TTerm, verbose);
                }
                if (beta_best != beta_c) {  // if the risk matrices aren't the optimal values, then they must be recalculated
                    // If it goes through every half step without improvement, then the maximum change needs to be decreased
                    abs_max = abs_max*pow(0.5, halfmax);  // reduces the step sizes
                    dose_abs_max = dose_abs_max*pow(0.5, halfmax);
                    iter_check = 1;
                    beta_p = beta_best;  //
                    beta_a = beta_best;  //
                    beta_c = beta_best;  //
                    for (int ij = 0; ij < totalnum; ij++) {
                        beta_0[ij] = beta_best[ij];
                    }
                    Cox_Term_Risk_Calc(modelform, tform, term_n, totalnum, fir, dfc, term_tot, T0, Td0, Tdd0, Te, R, Rd, Rdd, Dose, nonDose, beta_0, df0, dose_abs_max, abs_max, TTerm, nonDose_LIN, nonDose_PLIN, nonDose_LOGLIN, RdR, RddR, nthreads, KeepConstant, verbose, model_bool, gmix_theta, gmix_term);
                }
            }
            Lld_worst = 0;
            for (int ij = 0; ij < reqrdnum; ij++) {
                if (abs(Lld[ij]) > Lld_worst) {
                    Lld_worst = abs(Lld[ij]);
                }
            }
            if ((iteration % (reqrdnum)) || (iter_check == 1)) {  // checks every set number of iterations
                iter_check = 0;
                if (Lld_worst < deriv_epsilon) {  // ends if the derivatives are low enough
                    iter_stop = 1;
                }
                if (abs_max < epsilon/10) {  // if the maximum change is too low, then it ends
                    iter_stop = 1;
                }
            }
        }
        if (Lld_worst < deriv_epsilon) {  // ends if the derivatives are low enough
            iter_stop = 1;
            convgd = TRUE;
        }
        conv_fin[guess] = convgd;
//        end_point = system_clock::now();
//        ending = time_point_cast<microseconds>(end_point).time_since_epoch().count();
//        Rcout << "Iteration " << guess << " step: " << ending - start << endl;
//        start_point = system_clock::now();
//        start = time_point_cast<microseconds>(start_point).time_since_epoch().count();
        // -----------------------------------------------
        // Performing Full Calculation to get full second derivative matrix
        // -----------------------------------------------
        fill(Ll.begin(), Ll.end(), 0.0);
        fill(Lld.begin(), Lld.end(), 0.0);
        fill(Lldd.begin(), Lldd.end(), 0.0);
        //
        if (model_bool["basic"]) {
            Make_subterms_Basic(totalnum, dfc, T0, beta_0, df0, nthreads);
            Make_Risks_Basic(totalnum, T0, R, Rd, Rdd, RdR, nthreads, df0, dfc, KeepConstant);
            RdR = (RdR.array().isFinite()).select(RdR, 0);
        } else if (model_bool["linear_err"]) {
            Make_subterms_Linear_ERR(totalnum, tform, dfc, nonDose_PLIN, nonDose_LOGLIN, beta_0, df0, nthreads, KeepConstant);
            Make_Risks_Linear_ERR(tform, dfc, df0, totalnum, R, Rd, Rdd, nonDose_PLIN, nonDose_LOGLIN, RdR, RddR, nthreads, KeepConstant);
            RdR = (RdR.array().isFinite()).select(RdR, 0);
            RddR = (RddR.array().isFinite()).select(RddR, 0);
            TTerm = R.col(0).array();
        } else {
            Make_subterms(totalnum, term_n, tform, dfc, fir, T0, Td0, Tdd0, Dose, nonDose, TTerm, nonDose_LIN, nonDose_PLIN, nonDose_LOGLIN, beta_0, df0, dint, dslp, nthreads, KeepConstant);
            Make_Risks(modelform, tform, term_n, totalnum, fir, T0, Td0, Tdd0, Te, R, Rd, Rdd, Dose, nonDose, TTerm, nonDose_LIN, nonDose_PLIN, nonDose_LOGLIN, RdR, RddR, nthreads, KeepConstant, gmix_theta, gmix_term);
            RdR = (RdR.array().isFinite()).select(RdR, 0);
            RddR = (RddR.array().isFinite()).select(RddR, 0);
        }
        //
        if (model_bool["strata"]) {
            if (model_bool["cr"]) {
                // strata_CR or strata_CR_single
                Calculate_Sides_Strata_CR(model_bool, RiskFail, RiskPairs_Strata, totalnum, ntime, R, Rd, Rdd, Rls1, Rls2, Rls3, Lls1, Lls2, Lls3, cens_weight, nthreads, Strata_vals, KeepConstant);  // strata_cr
            } else {
                Calculate_Sides_Strata(model_bool, RiskFail, RiskPairs_Strata, totalnum, ntime, R, Rd, Rdd, Rls1, Rls2, Rls3, Lls1, Lls2, Lls3, nthreads, Strata_vals, KeepConstant);
            }
        } else if (model_bool["cr"]) {
            Calculate_Sides_CR(model_bool, RiskFail, RiskPairs, totalnum, ntime, R, Rd, Rdd, Rls1, Rls2, Rls3, Lls1, Lls2, Lls3, cens_weight, nthreads, KeepConstant);
        } else {
            Calculate_Sides(model_bool, RiskFail, RiskPairs, totalnum, ntime, R, Rd, Rdd, Rls1, Rls2, Rls3, Lls1, Lls2, Lls3, nthreads, KeepConstant);
        }
        if (model_bool["strata"]) {
            if (model_bool["basic"]) {
                Calc_LogLik_Strata_Basic(model_bool, nthreads, RiskFail, RiskPairs_Strata, totalnum, ntime, R, Rd, Rdd, RdR, Rls1, Rls2, Rls3, Lls1, Lls2, Lls3, Ll, Lld, Lldd, ties_method, Strata_vals, KeepConstant);
            } else if (model_bool["linear_err"]) {
                Calc_LogLik_Strata_Linear_ERR(model_bool, tform, nthreads, RiskFail, RiskPairs_Strata, totalnum, ntime, R, Rd, Rdd, RdR, RddR, Rls1, Rls2, Rls3, Lls1, Lls2, Lls3, Ll, Lld, Lldd, ties_method, Strata_vals, KeepConstant);
            } else {
                Calc_LogLik_Strata(model_bool, nthreads, RiskFail, RiskPairs_Strata, totalnum, ntime, R, Rd, Rdd, RdR, RddR, Rls1, Rls2, Rls3, Lls1, Lls2, Lls3, Ll, Lld, Lldd, ties_method, Strata_vals, KeepConstant);
            }
        } else {
            if (model_bool["basic"]) {
                Calc_LogLik_Basic(model_bool, nthreads, RiskFail, RiskPairs, totalnum, ntime, R, Rd, Rdd, RdR, Rls1, Rls2, Rls3, Lls1, Lls2, Lls3, Ll, Lld, Lldd, ties_method, KeepConstant);
            } else if (model_bool["linear_err"]) {
                Calc_LogLik_Linear_ERR(model_bool, tform, nthreads, RiskFail, RiskPairs, totalnum, ntime, R, Rd, Rdd, RdR, RddR, Rls1, Rls2, Rls3, Lls1, Lls2, Lls3, Ll, Lld, Lldd, ties_method, KeepConstant);
            } else {
                Calc_LogLik(model_bool, nthreads, RiskFail, RiskPairs, totalnum, ntime, R, Rd, Rdd, RdR, RddR, Rls1, Rls2, Rls3, Lls1, Lls2, Lls3, Ll, Lld, Lldd, ties_method, KeepConstant);
            }
        }
//        end_point = system_clock::now();
//        ending = time_point_cast<microseconds>(end_point).time_since_epoch().count();
//        Rcout << "Last Guess " << guess << " step: " << ending - start << endl;
//        start_point = system_clock::now();
//        start = time_point_cast<microseconds>(start_point).time_since_epoch().count();
        // Cox_Term_Risk_Calc(modelform, tform, term_n, totalnum, fir, dfc, term_tot, T0, Td0, Tdd0, Te, R, Rd, Rdd, Dose, nonDose, beta_0, df0, dose_abs_max, abs_max, TTerm, nonDose_LIN, nonDose_PLIN, nonDose_LOGLIN, RdR, RddR, nthreads, KeepConstant, verbose, model_bool, gmix_theta, gmix_term);
        // Cox_Side_LL_Calc(reqrdnum, ntime, tform, RiskFail,  totalnum, fir, R, Rd, Rdd, Rls1, Rls2, Rls3, Lls1, Lls2, Lls3, cens_weight, Strata_vals, beta_0, RdR, RddR, Ll, Lld, Lldd, nthreads, KeepConstant, ties_method, verbose, model_bool, iter_stop);
        a_n = beta_0;
        beta_fin(guess, _) = a_n;
        LL_fin[guess] = Ll[0];
        AIC_fin[guess] = 2*(totalnum-accumulate(KeepConstant.begin(), KeepConstant.end(), 0.0))-2*Ll[0];
        MatrixXd cov;
        NumericVector stdev(totalnum);
        if (model_bool["oberved_info"]){
            NumericVector Lldd_vec(reqrdnum * reqrdnum);  // simplfied information matrix
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
            cov = - 1 * Lldd_mat.inverse().matrix();  // uses inverse information matrix to calculate the standard deviation
            for (int ij = 0; ij < totalnum; ij++) {
                if (KeepConstant[ij] == 0) {
                    int pij_ind = ij - sum(head(KeepConstant, ij));
                    stdev(ij) = sqrt(cov(pij_ind, pij_ind));
                }
            }
        } else {
            //
            vector<double> InMa(pow(reqrdnum, 2), 0.0);
            if (model_bool["strata"]){
                if (model_bool["cr"]){
                    Expected_Inform_Matrix_Cox_Strata_CR(nthreads, RiskFail, RiskPairs_Strata, totalnum, ntime, R, Rd, RdR, cens_weight, InMa, Strata_vals, KeepConstant);
                } else {
                    Expected_Inform_Matrix_Cox_Strata(nthreads, RiskFail, RiskPairs_Strata, totalnum, ntime, R, Rd, RdR, InMa, Strata_vals, KeepConstant);
                }
            } else {
                if (model_bool["cr"]){
                    Expected_Inform_Matrix_Cox_CR(nthreads, RiskFail, RiskPairs, totalnum, ntime, R, Rd, RdR, cens_weight, InMa, KeepConstant);
                } else {
                    Expected_Inform_Matrix_Cox(nthreads, RiskFail, RiskPairs, totalnum, ntime, R, Rd, RdR, InMa, KeepConstant);
                }
            }
            NumericVector InMa_vec(reqrdnum * reqrdnum);  // simplfied information matrix
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
            cov = InMa_mat.inverse().matrix();  // uses inverse information matrix to calculate the standard deviation
            for (int ij = 0; ij < totalnum; ij++) {
                if (KeepConstant[ij] == 0) {
                    int pij_ind = ij - sum(head(KeepConstant, ij));
                    stdev(ij) = sqrt(cov(pij_ind, pij_ind));
                }
            }
        }
//        end_point = system_clock::now();
//        ending = time_point_cast<microseconds>(end_point).time_since_epoch().count();
//        Rcout << "Deviation " << guess << " step: " << ending - start << endl;
//        start_point = system_clock::now();
//        start = time_point_cast<microseconds>(start_point).time_since_epoch().count();
        std_fin(guess, _) = stdev;
//        end_point = system_clock::now();
//        ending = time_point_cast<microseconds>(end_point).time_since_epoch().count();
//        Rcout << "Assigned " << guess << " step: " << ending - start << endl;
    }
//    start_point = system_clock::now();
//    start = time_point_cast<microseconds>(start_point).time_since_epoch().count();
    List para_list;
    if (!model_bool["basic"]) {
        para_list = List::create(_["term_n"] = term_n, _["tforms"] = tform);  // stores the term information
    }
    List res_list;  // = List::create(_["LogLik"] = wrap(LL_fin), _["parameters"] = wrap(beta_fin), _["error"] = wrap(std_fin));
    if (model_bool["basic"]) {
        res_list = List::create(_["LogLik"] = wrap(LL_fin), _["AIC"] = wrap(AIC_fin), _["Parameters"] = wrap(beta_fin), _["Standard_Error"] = wrap(std_fin), _["Convergance"] = wrap(conv_fin), _["Status"] = "PASSED");
    } else {
        res_list = List::create(_["LogLik"] = wrap(LL_fin), _["AIC"] = wrap(AIC_fin), _["Parameters"] = wrap(beta_fin), _["Standard_Error"] = wrap(std_fin), _["Parameter_Lists"] = para_list, _["Convergance"] = wrap(conv_fin), _["Status"] = "PASSED");
    }
//    end_point = system_clock::now();
//    ending = time_point_cast<microseconds>(end_point).time_since_epoch().count();
//
//    Rcout << "Wrapping step: " << ending - start << endl;
    // returns a list of results
    return res_list;
}

//' Primary Cox PH regression with multiple distributed dose columns and optional combinations of null, stratification, competing risks, multiplicative log-linear model, and no derivative calculation.
//'
//' \code{LogLik_Cox_PH_Multidose_Omnibus_Integrated} Performs the calls to calculation functions, Structures the Cox PH regression, With verbose option prints out time stamps and intermediate sums of terms and derivatives
//'
//' @inheritParams CPP_template
//'
//' @return List of final results: standard cox outputs for the integrated solution
//' @noRd
//'
// [[Rcpp::export]]
List LogLik_Cox_PH_Multidose_Omnibus_Integrated(IntegerVector term_n, StringVector tform, NumericVector a_n, NumericMatrix& x_all, NumericMatrix& dose_all, IntegerMatrix dose_cols, IntegerVector dose_index, IntegerVector dfc, int fir, string modelform, double lr, List optim_para, int maxiter, int halfmax, double epsilon, double abs_max, double dose_abs_max, double deriv_epsilon, NumericMatrix df_groups, NumericVector tu, int double_step, int verbose, IntegerVector KeepConstant, int term_tot, string ties_method, int nthreads, NumericVector& Strata_vals, const VectorXd& cens_weight, List model_bool, const double gmix_theta, const IntegerVector gmix_term, const MatrixXd Lin_Sys, const VectorXd Lin_Res) {
    //
    List temp_list = List::create(_["Status"] = "TEMP");  // used as a dummy return value for code checking
    // Time durations are measured from this point on in microseconds
    //
    // df0: covariate data
    // ntime: number of event times for Cox PH
    // totalnum: number of terms used
    //
    // ------------------------------------------------------------------------- // initialize
    Map<MatrixXd> df0(as<Map<MatrixXd> >(x_all));
    const Map<MatrixXd> df1(as<Map<MatrixXd> >(dose_all));
    const int mat_row = df0.rows();
    int ntime = tu.size();
    int totalnum = term_n.size();
    int reqrdnum = totalnum - sum(KeepConstant);
    // ------------------------------------------------------------------------- // initialize
    if (model_bool["null"]) {
        if (verbose >= 1) {
            Rcout << "null model is not compatable with multi-realization method" << endl;
        }
        temp_list = List::create(_["Status"] = "FAILED_WITH_BAD_MODEL_NULL", _["LogLik"] = R_NaN);
        return temp_list;
    }
    if (model_bool["single"]) {
        if (verbose >= 1) {
            Rcout << "non-derivative model calculation is not compatable with multi-realization method" << endl;
        }
        temp_list = List::create(_["Status"] = "FAILED_WITH_BAD_MODEL_SINGLE", _["LogLik"] = R_NaN);
        return temp_list;
    }
    // if (model_bool["gradient"]) {
    //     if (verbose >= 1) {
    //         Rcout << "gradient descent model calculation is not CURRENTLY compatable with multi-realization method" << endl;
    //     }
    //     temp_list = List::create(_["Status"] = "FAILED_WITH_BAD_MODEL_GRADIENT", _["LogLik"] = R_NaN);
    //     return temp_list;
    // }
    //
    // cout.precision: controls the number of significant digits printed
    // nthreads: number of threads used for parallel operations
    //
    Rcout.precision(7);  // forces higher precision numbers printed to terminal
    //
    // Lld_worst: stores the highest magnitude log-likelihood derivative
    double Lld_worst = 0.0;  // stores derivative value used to determine if every parameter is near convergence
    //
    // ---------------------------------------------
    // To Start, needs to seperate the derivative terms
    // ---------------------------------------------
    // ------------------------------------------------------------------------- // initialize
    Map<VectorXd> beta_0(as<Map<VectorXd> >(a_n));
    MatrixXd T0;
    MatrixXd Td0;
    MatrixXd Tdd0;
    MatrixXd Te;
    MatrixXd R;
    ColXd Rd;
    ColXd Rdd;
    MatrixXd Dose;
    MatrixXd nonDose;
    MatrixXd nonDose_LIN;
    MatrixXd nonDose_PLIN;
    MatrixXd nonDose_LOGLIN;
    MatrixXd TTerm;
    double dint = 0.0;  // the amount of change used to calculate derivatives in threshold paramters
    double dslp = 0.0;
    ColXd RdR;
    ColXd RddR;
    // ------------------------------------------------------------------------- // initialize
    // ---------------------------------------------
    // To Start, needs to seperate the derivative terms
    // ---------------------------------------------
    //
    Cox_Refresh_R_TERM(totalnum, reqrdnum, term_tot, dint, dslp, dose_abs_max, abs_max, df0, T0, Td0, Tdd0, Te, R, Rd, Rdd, Dose, nonDose, TTerm, nonDose_LIN, nonDose_PLIN, nonDose_LOGLIN, RdR, RddR, model_bool);
    IntegerMatrix RiskFail;
    vector<vector<int> > RiskPairs(ntime);
    vector<vector<vector<int> > > RiskPairs_Strata(ntime, vector<vector<int>>(Strata_vals.size()));
    const Map<MatrixXd> df_m(as<Map<MatrixXd> >(df_groups));
    if (model_bool["strata"]) {
        RiskFail = IntegerMatrix(ntime, 2*Strata_vals.size());  // vector giving the event rows
        // Creates matrices used to identify the event risk groups
        if (model_bool["cr"]) {
            Make_Groups_Strata_CR(ntime, df_m, RiskFail, RiskPairs_Strata, tu, nthreads, Strata_vals, cens_weight);
        } else {
            Make_Groups_Strata(ntime, df_m, RiskFail, RiskPairs_Strata, tu, nthreads, Strata_vals);
        }
    } else {
        RiskFail = IntegerMatrix(ntime, 2);  // vector giving the event rows
        // Creates matrices used to identify the event risk groups
        if (model_bool["cr"]) {
            Make_Groups_CR(ntime, df_m, RiskFail, RiskPairs, tu, cens_weight, nthreads);
        } else {
            Make_Groups(ntime, df_m, RiskFail, RiskPairs, tu, nthreads);
        }
    }
    // ------------------------------------------------------------------------- // initialize
    MatrixXd Rls1;
    MatrixXd Lls1;
    MatrixXd Rls2;
    MatrixXd Rls3;
    MatrixXd Lls2;
    MatrixXd Lls3;
    vector<double> Ll(reqrdnum, 0.0);  // log-likelihood values
    vector<double> Lld(reqrdnum, 0.0);  // log-likelihood derivative values
    vector<double> Lldd(pow(reqrdnum, 2), 0.0);  // the second derivative matrix has room for every combination, but only the lower triangle is calculated initially
    //
    vector<double> Ll_Total(reqrdnum, 0.0);  // log-likelihood values
    vector<double> Lld_Total(reqrdnum, 0.0);  // log-likelihood derivative values
    vector<double> Lldd_Total(pow(reqrdnum, 2), 0.0);  // the second derivative matrix has room for every combination, but only the lower triangle is calculated initially
    // ------------------------------------------------------------------------- // initialize
    Cox_Refresh_R_SIDES(reqrdnum, ntime, Rls1, Rls2, Rls3, Lls1, Lls2, Lls3, Strata_vals, model_bool);
    // the log-likelihood is calculated in parallel over the risk groups
    vector <double> Ll_comp(2, Ll[0]);  // vector to compare values
    double abs_max0 = abs_max;
    double dose_abs_max0 = dose_abs_max;
    vector<double> dbeta(totalnum, 0.0);
    NumericVector m_g_store(reqrdnum);
    NumericVector v_beta_store(reqrdnum);
    // --------------------------
    // always starts from initial guess
    // --------------------------
    vector<double> beta_c(totalnum, 0.0);
    vector<double> beta_a(totalnum, 0.0);
    vector<double> beta_best(totalnum, 0.0);
    vector<double> beta_p(totalnum, 0.0);
    VectorXd::Map(&beta_p[0], beta_0.size()) = beta_0;  // stores previous parameters
    VectorXd::Map(&beta_c[0], beta_0.size()) = beta_0;  // stores current parameters
    VectorXd::Map(&beta_a[0], beta_0.size()) = beta_0;  // stores a refrence value for parameters
    VectorXd::Map(&beta_best[0], beta_0.size()) = beta_0;  // stores the best parameters
    double halves = 0;  // number of half-steps taken
    int ind0 = fir;  // used for validations
    int iteration = 0;  // iteration number
    int iter_stop  = 0;  // tracks if the iterations should be stopped for convergence
    int iter_check = 0;  // signal to check for convergence
    int guesses = dose_cols.cols();
//    NumericMatrix beta_fin(dose_cols.cols(), totalnum);
//    NumericVector LL_fin(dose_cols.cols());
//    LogicalVector conv_fin(dose_cols.cols());
//    NumericMatrix std_fin(dose_cols.cols(), totalnum);
    IntegerVector dfc_0(dfc.length());
    for (int i = 0;  i<dfc.length(); i++) {
        dfc_0[i] = dfc[i];
    }
    //
    double Ll_abs_best = 10;
    vector<double> beta_abs_best(totalnum, 0.0);
    // Variables that are used for the risk check function shared across cox, poisson, and log bound functions
    // double dev = 0.0;
    MatrixXd dev_temp = MatrixXd::Zero(1, 1);
    // double Lstar = 0.0;
    MatrixXd PyrC = MatrixXd::Zero(1, 1);
    bool convgd = FALSE;
    //
    Cox_Refresh_R_SIDES(reqrdnum, ntime, Rls1, Rls2, Rls3, Lls1, Lls2, Lls3, Strata_vals, model_bool);
    fill(Ll.begin(), Ll.end(), 0.0);
    fill(Lld.begin(), Lld.end(), 0.0);
    if (!model_bool["gradient"]){
        fill(Lldd.begin(), Lldd.end(), 0.0);
    }
    if (model_bool["gradient"]) {
        m_g_store.fill(0);
        v_beta_store.fill(0);
    }
    beta_p = beta_best;  //
    beta_a = beta_best;  //
    beta_c = beta_best;  //
    abs_max = abs_max0;
    dose_abs_max = dose_abs_max0;
    iter_stop = 0;
    halves = 0;
    convgd = FALSE;
    iteration = 0;
    halves = 0;  // number of half-steps taken
    ind0 = fir;  // used for validations
    iteration = 0;  // iteration number
    Ll_abs_best = 10;
    iter_stop  = 0;  // tracks if the iterations should be stopped for convergence
    iter_check = 0;  // signal to check for convergence
    for (int i = 0; i < beta_0.size(); i++) {
        beta_0[i] = beta_best[i];
    }
    for (int guess = 0; guess <guesses; guess++) { // First loop to identify a possible beta value
        for (int i = 0; i<dose_index.size(); i++) {
            for (int j = 0; j < totalnum; j++) {
                if (dfc_0[j] == dose_index[i]) {
//                    dfc[j] = dose_cols(i, guess);
                    df0.col(dfc[j]-1) = df1.col(dose_cols(i, guess)-1);
                }
            }
        }
        Cox_Term_Risk_Calc(modelform, tform, term_n, totalnum, fir, dfc, term_tot, T0, Td0, Tdd0, Te, R, Rd, Rdd, Dose, nonDose, beta_0, df0, dose_abs_max, abs_max, TTerm, nonDose_LIN, nonDose_PLIN, nonDose_LOGLIN, RdR, RddR, nthreads, KeepConstant, verbose, model_bool, gmix_theta, gmix_term);
        int check_loop = 0;
        while (((R.minCoeff() <= 0) || (R.hasNaN())) && (check_loop < 3)) { // If a linear column changes, then a previous optimum can give a negative risk
            check_loop = check_loop + 1; // the parameters in a negative term are reduced to zero, which should generally scale the term value closer to zero as well
            for (int ijk = 0; ijk < totalnum; ijk++) {
                int tij = term_n[ijk];
                if (TTerm.col(tij).minCoeff()<=0) {
                    beta_0[ijk] = beta_0[ijk] / 2.0;
                } else if (isinf(TTerm.col(tij).maxCoeff())) {
                    beta_0[ijk] = beta_0[ijk] / 2.0;
                } else if (isnan(TTerm.col(tij).minCoeff())) {
                    beta_0[ijk] = beta_0[ijk] / 2.0;
                }
            }
            Cox_Term_Risk_Calc(modelform, tform, term_n, totalnum, fir, dfc, term_tot, T0, Td0, Tdd0, Te, R, Rd, Rdd, Dose, nonDose, beta_0, df0, dose_abs_max, abs_max, TTerm, nonDose_LIN, nonDose_PLIN, nonDose_LOGLIN, RdR, RddR, nthreads, KeepConstant, verbose, model_bool, gmix_theta, gmix_term);
        }
        if ((R.minCoeff() <= 0) || (R.hasNaN())) {
            if (verbose >= 1) {
                Rcout << "C++ Error: A non-positive risk was detected: " << R.minCoeff() << endl;
                Rcout << "C++ Warning: final failing values ";
                for (int ijk = 0; ijk < totalnum; ijk++) {
                    Rcout << beta_0[ijk] << " ";
                }
                Rcout << " " << endl;
            }
            temp_list = List::create(_["beta_0"] = wrap(beta_0), _["Deviation"] = R_NaN, _["Status"] = "FAILED_WITH_NEGATIVE_RISK", _["LogLik"] = R_NaN);
            return temp_list;
        }
    }
    fill(Ll_Total.begin(), Ll_Total.end(), 0.0);
    fill(Lld_Total.begin(), Lld_Total.end(), 0.0);
    if (!model_bool["gradient"]){
        fill(Lldd_Total.begin(), Lldd_Total.end(), 0.0);
    }
    for (int guess = 0; guess <guesses; guess++) { // second loop to run the likelihood calculations
        for (int i = 0; i<dose_index.size(); i++) {
            for (int j = 0; j < totalnum; j++) {
                if (dfc_0[j] == dose_index[i]) {
//                    dfc[j] = dose_cols(i, guess);
                    df0.col(dfc[j]-1) = df1.col(dose_cols(i, guess)-1);
                }
            }
        }
        Cox_Refresh_R_SIDES(reqrdnum, ntime, Rls1, Rls2, Rls3, Lls1, Lls2, Lls3, Strata_vals, model_bool);
        Cox_Term_Risk_Calc(modelform, tform, term_n, totalnum, fir, dfc, term_tot, T0, Td0, Tdd0, Te, R, Rd, Rdd, Dose, nonDose, beta_0, df0, dose_abs_max, abs_max, TTerm, nonDose_LIN, nonDose_PLIN, nonDose_LOGLIN, RdR, RddR, nthreads, KeepConstant, verbose, model_bool, gmix_theta, gmix_term);
        fill(Ll.begin(), Ll.end(), 0.0);
        fill(Lld.begin(), Lld.end(), 0.0);
        if (!model_bool["gradient"]){
            fill(Lldd.begin(), Lldd.end(), 0.0);
        }
        Cox_Side_LL_Calc(reqrdnum, ntime, tform, RiskFail,  RiskPairs, RiskPairs_Strata, totalnum, fir, R, Rd, Rdd, Rls1, Rls2, Rls3, Lls1, Lls2, Lls3, cens_weight, Strata_vals, beta_0, RdR, RddR, Ll, Lld, Lldd, nthreads, KeepConstant, ties_method, verbose, model_bool, iter_stop);
        //
        for (int i=0; i < reqrdnum; i ++){
            Ll_Total[i] += Ll[i] / guesses;
            Lld_Total[i] += Lld[i] / guesses;
        }
        for (int i=0; i < pow(reqrdnum,2); i ++){
            Lldd_Total[i] += Lldd[i] / guesses;
        }
    }
    //
    Print_LL(reqrdnum, totalnum, beta_0, Ll_Total, Lld_Total, Lldd_Total, verbose, model_bool);
    //
    for (int i = 0; i < beta_0.size(); i++) {
        beta_c[i] = beta_0[i];
    }
    while ((iteration < maxiter) && (iter_stop == 0)) {
        iteration++;
        beta_p = beta_c;  //
        beta_a = beta_c;  //
        beta_best = beta_c;  //
        // calculates the initial change in parameter
        if (model_bool["basic"]) {
            Calc_Change_Basic(double_step, nthreads, totalnum, lr, abs_max, Ll_Total, Lld_Total, Lldd_Total, dbeta, KeepConstant);
        } else if (model_bool["gradient"]) {
            Calc_Change_Gradient(nthreads, model_bool, totalnum, optim_para, iteration, abs_max, Lld_Total, m_g_store, v_beta_store, dbeta, KeepConstant);
        } else {
            if (model_bool["constraint"]) {
                Calc_Change_Cons(Lin_Sys, Lin_Res, beta_0, nthreads, totalnum, dose_abs_max, lr, abs_max, Ll_Total, Lld_Total, Lldd_Total, dbeta, tform, dose_abs_max, abs_max, KeepConstant);
            } else {
                Calc_Change(double_step, nthreads, totalnum, dose_abs_max, lr, abs_max, Ll_Total, Lld_Total, Lldd_Total, dbeta, tform, dose_abs_max, abs_max, KeepConstant);
            }
            Intercept_Bound(nthreads, totalnum, beta_0, dbeta, dfc, df0, KeepConstant, tform);
        }
        //
        if ((Ll_abs_best > 0) || (Ll_abs_best < Ll_Total[ind0])) {
            Ll_abs_best = Ll_Total[ind0];
            beta_abs_best = beta_c;
        }
        // Run the integrated feasible value check
        bool feasible_pass = TRUE;
        if (model_bool["gradient"]){
            for (int ij = 0; ij < totalnum; ij++) {
                beta_0[ij] = beta_a[ij] + dbeta[ij];
                beta_c[ij] = beta_0[ij];
            }
            for (int guess = 0; guess <guesses; guess++) {
                for (int i = 0; i<dose_index.size(); i++) {
                    for (int j = 0; j < totalnum; j++) {
                        if (dfc_0[j] == dose_index[i]) {
//                            dfc[j] = dose_cols(i, guess);
                            df0.col(dfc[j]-1) = df1.col(dose_cols(i, guess)-1);
                        }
                    }
                }
                Cox_Term_Risk_Calc(modelform, tform, term_n, totalnum, fir, dfc, term_tot, T0, Td0, Tdd0, Te, R, Rd, Rdd, Dose, nonDose, beta_0, df0, dose_abs_max, abs_max, TTerm, nonDose_LIN, nonDose_PLIN, nonDose_LOGLIN, RdR, RddR, nthreads, KeepConstant, verbose, model_bool, gmix_theta, gmix_term);
                if ((R.minCoeff() <= 0) || (R.hasNaN())) {
                    #ifdef _OPENMP
                    #pragma omp parallel for schedule(dynamic) num_threads(nthreads)
                    #endif
                    for (int ijk = 0; ijk < totalnum; ijk++) {
                        int tij = term_n[ijk];
                        if (TTerm.col(tij).minCoeff()<=0) {
                            dbeta[ijk] = dbeta[ijk] / 2.0;
                        } else if (isinf(TTerm.col(tij).maxCoeff())) {
                            dbeta[ijk] = dbeta[ijk] / 2.0;
                        } else if (isnan(TTerm.col(tij).minCoeff())) {
                            dbeta[ijk] = dbeta[ijk] / 2.0;
                        }
                    }
                    halves+=0.2;
                    feasible_pass = FALSE;
                    guess = guesses;
                }
            }
            fill(Ll_Total.begin(), Ll_Total.end(), 0.0);
            fill(Lld_Total.begin(), Lld_Total.end(), 0.0);
            if (!model_bool["gradient"]){
                fill(Lldd_Total.begin(), Lldd_Total.end(), 0.0);
            }
            if (feasible_pass){
                halves++;
                for (int ij = 0; ij < totalnum; ij++) {
                    beta_0[ij] = beta_a[ij] + dbeta[ij];
                    beta_c[ij] = beta_0[ij];
                }
                // the point was feasible for every realization
                // gather the new combined score
                for (int guess = 0; guess <guesses; guess++) {
                    for (int i = 0; i<dose_index.size(); i++) {
                        for (int j = 0; j < totalnum; j++) {
                            if (dfc_0[j] == dose_index[i]) {
//                                dfc[j] = dose_cols(i, guess);
                                df0.col(dfc[j]-1) = df1.col(dose_cols(i, guess)-1);
                            }
                        }
                    }
                    Cox_Term_Risk_Calc(modelform, tform, term_n, totalnum, fir, dfc, term_tot, T0, Td0, Tdd0, Te, R, Rd, Rdd, Dose, nonDose, beta_0, df0, dose_abs_max, abs_max, TTerm, nonDose_LIN, nonDose_PLIN, nonDose_LOGLIN, RdR, RddR, nthreads, KeepConstant, verbose, model_bool, gmix_theta, gmix_term);
                    fill(Ll.begin(), Ll.end(), 0.0);
                    fill(Lld.begin(), Lld.end(), 0.0);
                    if (!model_bool["gradient"]){
                        fill(Lldd.begin(), Lldd.end(), 0.0);
                    }
                    Cox_Refresh_R_SIDES(reqrdnum, ntime, Rls1, Rls2, Rls3, Lls1, Lls2, Lls3, Strata_vals, model_bool);
                    Cox_Side_LL_Calc(reqrdnum, ntime, tform, RiskFail,  RiskPairs, RiskPairs_Strata, totalnum, fir, R, Rd, Rdd, Rls1, Rls2, Rls3, Lls1, Lls2, Lls3, cens_weight, Strata_vals, beta_0, RdR, RddR, Ll, Lld, Lldd, nthreads, KeepConstant, ties_method, verbose, model_bool, iter_stop);
                    for (int i=0; i < reqrdnum; i ++){
                        Ll_Total[i] += Ll[i] / guesses;
                        Lld_Total[i] += Lld[i] / guesses;
                    }
                    for (int i=0; i < pow(reqrdnum,2); i ++){
                        Lldd_Total[i] += Lldd[i] / guesses;
                    }
                }
                Print_LL(reqrdnum, totalnum, beta_0, Ll_Total, Lld_Total, Lldd_Total, verbose, model_bool);
            }
                // check if any have issues
        } else {
            halves = 0;
            while ((Ll_Total[ind0] <= Ll_abs_best) && (halves < halfmax)) {  // repeats until half-steps maxed or an improvement
                for (int ij = 0; ij < totalnum; ij++) {
                    beta_0[ij] = beta_a[ij] + dbeta[ij];
                    beta_c[ij] = beta_0[ij];
                }
                for (int guess = 0; guess <guesses; guess++) {
                    for (int i = 0; i<dose_index.size(); i++) {
                        for (int j = 0; j < totalnum; j++) {
                            if (dfc_0[j] == dose_index[i]) {
//                                dfc[j] = dose_cols(i, guess);
                                df0.col(dfc[j]-1) = df1.col(dose_cols(i, guess)-1);
                            }
                        }
                    }
                    Cox_Term_Risk_Calc(modelform, tform, term_n, totalnum, fir, dfc, term_tot, T0, Td0, Tdd0, Te, R, Rd, Rdd, Dose, nonDose, beta_0, df0, dose_abs_max, abs_max, TTerm, nonDose_LIN, nonDose_PLIN, nonDose_LOGLIN, RdR, RddR, nthreads, KeepConstant, verbose, model_bool, gmix_theta, gmix_term);
                    if ((R.minCoeff() <= 0) || (R.hasNaN())) {
                        #ifdef _OPENMP
                        #pragma omp parallel for schedule(dynamic) num_threads(nthreads)
                        #endif
                        for (int ijk = 0; ijk < totalnum; ijk++) {
                            int tij = term_n[ijk];
                            if (TTerm.col(tij).minCoeff()<=0) {
                                dbeta[ijk] = dbeta[ijk] / 2.0;
                            } else if (isinf(TTerm.col(tij).maxCoeff())) {
                                dbeta[ijk] = dbeta[ijk] / 2.0;
                            } else if (isnan(TTerm.col(tij).minCoeff())) {
                                dbeta[ijk] = dbeta[ijk] / 2.0;
                            }
                        }
                        halves+=0.2;
                        feasible_pass = FALSE;
                        guess = guesses;
                    }
                }
                if (feasible_pass){
                    halves++;
                    fill(Ll_Total.begin(), Ll_Total.end(), 0.0);
                    fill(Lld_Total.begin(), Lld_Total.end(), 0.0);
                    if (!model_bool["gradient"]){
                        fill(Lldd_Total.begin(), Lldd_Total.end(), 0.0);
                    }
                    // the point was feasible for every realization
                    // gather the new combined score
                    for (int guess = 0; guess <guesses; guess++) {
                        for (int i = 0; i<dose_index.size(); i++) {
                            for (int j = 0; j < totalnum; j++) {
                                if (dfc_0[j] == dose_index[i]) {
//                                    dfc[j] = dose_cols(i, guess);
                                    df0.col(dfc[j]-1) = df1.col(dose_cols(i, guess)-1);
                                }
                            }
                        }
                        Cox_Term_Risk_Calc(modelform, tform, term_n, totalnum, fir, dfc, term_tot, T0, Td0, Tdd0, Te, R, Rd, Rdd, Dose, nonDose, beta_0, df0, dose_abs_max, abs_max, TTerm, nonDose_LIN, nonDose_PLIN, nonDose_LOGLIN, RdR, RddR, nthreads, KeepConstant, verbose, model_bool, gmix_theta, gmix_term);
                        fill(Ll.begin(), Ll.end(), 0.0);
                        fill(Lld.begin(), Lld.end(), 0.0);
                        if (!model_bool["gradient"]){
                            fill(Lldd.begin(), Lldd.end(), 0.0);
                        }
                        Cox_Refresh_R_SIDES(reqrdnum, ntime, Rls1, Rls2, Rls3, Lls1, Lls2, Lls3, Strata_vals, model_bool);
                        Cox_Side_LL_Calc(reqrdnum, ntime, tform, RiskFail,  RiskPairs, RiskPairs_Strata, totalnum, fir, R, Rd, Rdd, Rls1, Rls2, Rls3, Lls1, Lls2, Lls3, cens_weight, Strata_vals, beta_0, RdR, RddR, Ll, Lld, Lldd, nthreads, KeepConstant, ties_method, verbose, model_bool, iter_stop);
                        for (int i=0; i < reqrdnum; i ++){
                            Ll_Total[i] += Ll[i] / guesses;
                            Lld_Total[i] += Lld[i] / guesses;
                        }
                        for (int i=0; i < pow(reqrdnum,2); i ++){
                            Lldd_Total[i] += Lldd[i] / guesses;
                        }
                    }
                    Print_LL(reqrdnum, totalnum, beta_0, Ll_Total, Lld_Total, Lldd_Total, verbose, model_bool);
                    if (Ll_Total[ind0] <= Ll_abs_best) {  // if a better point wasn't found, takes a half-step
                        #ifdef _OPENMP
                        #pragma omp parallel for schedule(dynamic) num_threads(nthreads)
                        #endif
                        for (int ijk = 0; ijk < totalnum; ijk++) {
                            dbeta[ijk] = dbeta[ijk] * 0.5;  //
                        }
                    } else{  // if improved, updates the best vector
                        #ifdef _OPENMP
                        #pragma omp parallel for schedule(dynamic) num_threads(nthreads)
                        #endif
                        for (int ijk = 0; ijk < totalnum; ijk++) {
                            beta_best[ijk] = beta_c[ijk];
                        }
                    }
                }
            }
            if (beta_best != beta_c) {  // if the risk matrices aren't the optimal values, then they must be recalculated
                // If it goes through every half step without improvement, then the maximum change needs to be decreased
                abs_max = abs_max*pow(0.5, halfmax);  // reduces the step sizes
                dose_abs_max = dose_abs_max*pow(0.5, halfmax);
                iter_check = 1;
                beta_p = beta_best;  //
                beta_a = beta_best;  //
                beta_c = beta_best;  //
                for (int ij = 0; ij < totalnum; ij++) {
                    beta_0[ij] = beta_best[ij];
                }
            }
        }
        Lld_worst = 0;
        for (int ij = 0; ij < reqrdnum; ij++) {
            if (abs(Lld_Total[ij]) > Lld_worst) {
                Lld_worst = abs(Lld_Total[ij]);
            }
        }
        if ((iteration % (reqrdnum)) || (iter_check == 1)) {  // checks every set number of iterations
            iter_check = 0;
            if (Lld_worst < deriv_epsilon) {  // ends if the derivatives are low enough
                iter_stop = 1;
            }
            if (abs_max < epsilon/10) {  // if the maximum change is too low, then it ends
                iter_stop = 1;
            }
        }
    }
    if (Lld_worst < deriv_epsilon) {  // ends if the derivatives are low enough
        iter_stop = 1;
        convgd = TRUE;
    }
    fill(Ll_Total.begin(), Ll_Total.end(), 0.0);
    fill(Lld_Total.begin(), Lld_Total.end(), 0.0);
    if (!model_bool["gradient"]){
        fill(Lldd_Total.begin(), Lldd_Total.end(), 0.0);
    }
    for (int guess = 0; guess <guesses; guess++) {
        for (int i = 0; i<dose_index.size(); i++) {
            for (int j = 0; j < totalnum; j++) {
                if (dfc_0[j] == dose_index[i]) {
//                    dfc[j] = dose_cols(i, guess);
                    df0.col(dfc[j]-1) = df1.col(dose_cols(i, guess)-1);
                }
            }
        }
        fill(Ll.begin(), Ll.end(), 0.0);
        fill(Lld.begin(), Lld.end(), 0.0);
        if (!model_bool["gradient"]){
            fill(Lldd.begin(), Lldd.end(), 0.0);
        }
        Cox_Refresh_R_SIDES(reqrdnum, ntime, Rls1, Rls2, Rls3, Lls1, Lls2, Lls3, Strata_vals, model_bool);
        Cox_Term_Risk_Calc(modelform, tform, term_n, totalnum, fir, dfc, term_tot, T0, Td0, Tdd0, Te, R, Rd, Rdd, Dose, nonDose, beta_0, df0, dose_abs_max, abs_max, TTerm, nonDose_LIN, nonDose_PLIN, nonDose_LOGLIN, RdR, RddR, nthreads, KeepConstant, verbose, model_bool, gmix_theta, gmix_term);
        Cox_Side_LL_Calc(reqrdnum, ntime, tform, RiskFail,  RiskPairs, RiskPairs_Strata, totalnum, fir, R, Rd, Rdd, Rls1, Rls2, Rls3, Lls1, Lls2, Lls3, cens_weight, Strata_vals, beta_0, RdR, RddR, Ll, Lld, Lldd, nthreads, KeepConstant, ties_method, verbose, model_bool, iter_stop);
        for (int i=0; i < reqrdnum; i ++){
            Ll_Total[i] += Ll[i] / guesses;
            Lld_Total[i] += Lld[i] / guesses;
        }
        for (int i=0; i < pow(reqrdnum,2); i ++){
            Lldd_Total[i] += Lldd[i] / guesses;
        }
    }
    Print_LL(reqrdnum, totalnum, beta_0, Ll_Total, Lld_Total, Lldd_Total, verbose, model_bool);
    if ((Ll_abs_best > 0) || (Ll_abs_best < Ll[ind0])) {
        Ll_abs_best = Ll[ind0];
        beta_abs_best = beta_c;
    }
    //
    List res_list;
    //
    if (model_bool["single"]) {
        res_list = List::create(_["LogLik"] = wrap(Ll_Total[0]), _["beta_0"] = wrap(beta_0), _["AIC"] = 2*(totalnum-accumulate(KeepConstant.begin(), KeepConstant.end(), 0.0))-2*Ll_Total[0], _["BIC"] = (totalnum-accumulate(KeepConstant.begin(), KeepConstant.end(), 0.0))*log(mat_row)-2*Ll_Total[0], _["Status"] = "PASSED");
        // returns a list of results
        return res_list;
    }
    List para_list;
    if (!model_bool["basic"]) {
        para_list = List::create(_["term_n"] = term_n, _["tforms"] = tform);  // stores the term information
    }
    List control_list = List::create(_["Iteration"] = iteration, _["Maximum Step"]= abs_max, _["Derivative Limiting"] = Lld_worst);  // stores the total number of iterations used
    //
    if (model_bool["gradient"]) {
        res_list = List::create(_["LogLik"] = wrap(Ll_Total[0]), _["First_Der"] = wrap(Lld_Total), _["beta_0"] = wrap(beta_0), _["AIC"] = 2*(totalnum-accumulate(KeepConstant.begin(), KeepConstant.end(), 0.0))-2*Ll_Total[0], _["BIC"] = (totalnum-accumulate(KeepConstant.begin(), KeepConstant.end(), 0.0))*log(mat_row)-2*Ll_Total[0], _["Parameter_Lists"] = para_list, _["Control_List"] = control_list, _["Converged"] = convgd, _["Status"] = "PASSED");
        return res_list;
    }
    //
    NumericVector Lldd_vec(reqrdnum * reqrdnum);  // simplfied information matrix
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
        Lldd_vec[ij * reqrdnum + jk] = Lldd_Total[ij*reqrdnum+jk];
        Lldd_vec[jk * reqrdnum + ij] = Lldd_vec[ij * reqrdnum + jk];
    }
    Lldd_vec.attr("dim") = Dimension(reqrdnum, reqrdnum);
    MatrixXd cov;
    VectorXd stdev = VectorXd::Zero(totalnum);
    if (model_bool["oberved_info"]){
        const Map<MatrixXd> Lldd_mat(as<Map<MatrixXd> >(Lldd_vec));
        //
        cov = - 1 * Lldd_mat.inverse().matrix();  // uses inverse information matrix to calculate the standard deviation
        for (int ij = 0; ij < totalnum; ij++) {
            if (KeepConstant[ij] == 0) {
                int pij_ind = ij - sum(head(KeepConstant, ij));
                stdev(ij) = sqrt(cov(pij_ind, pij_ind));
            }
        }
    } else {
        //
        vector<double> InMa(pow(reqrdnum, 2), 0.0);
        if (model_bool["strata"]){
            if (model_bool["cr"]){
                Expected_Inform_Matrix_Cox_Strata_CR(nthreads, RiskFail, RiskPairs_Strata, totalnum, ntime, R, Rd, RdR, cens_weight, InMa, Strata_vals, KeepConstant);
            } else {
                Expected_Inform_Matrix_Cox_Strata(nthreads, RiskFail, RiskPairs_Strata, totalnum, ntime, R, Rd, RdR, InMa, Strata_vals, KeepConstant);
            }
        } else {
            if (model_bool["cr"]){
                Expected_Inform_Matrix_Cox_CR(nthreads, RiskFail, RiskPairs, totalnum, ntime, R, Rd, RdR, cens_weight, InMa, KeepConstant);
            } else {
                Expected_Inform_Matrix_Cox(nthreads, RiskFail, RiskPairs, totalnum, ntime, R, Rd, RdR, InMa, KeepConstant);
            }
        }
        NumericVector InMa_vec(reqrdnum * reqrdnum);  // simplfied information matrix
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
        cov = InMa_mat.inverse().matrix();  // uses inverse information matrix to calculate the standard deviation
        for (int ij = 0; ij < totalnum; ij++) {
            if (KeepConstant[ij] == 0) {
                int pij_ind = ij - sum(head(KeepConstant, ij));
                stdev(ij) = sqrt(cov(pij_ind, pij_ind));
            }
        }
    }
    //
    if (model_bool["basic"]) {
        res_list = List::create(_["LogLik"] = wrap(Ll_Total[0]), _["First_Der"] = wrap(Lld_Total), _["Second_Der"] = Lldd_vec, _["beta_0"] = wrap(beta_0), _["Standard_Deviation"] = wrap(stdev), _["Covariance"] = wrap(cov), _["AIC"] = 2*(totalnum-accumulate(KeepConstant.begin(), KeepConstant.end(), 0.0))-2*Ll_Total[0], _["BIC"] = (totalnum-accumulate(KeepConstant.begin(), KeepConstant.end(), 0.0))*log(mat_row)-2*Ll_Total[0], _["Control_List"] = control_list, _["Converged"] = convgd, _["Status"] = "PASSED");
    } else {
        res_list = List::create(_["LogLik"] = wrap(Ll_Total[0]), _["First_Der"] = wrap(Lld_Total), _["Second_Der"] = Lldd_vec, _["beta_0"] = wrap(beta_0), _["Standard_Deviation"] = wrap(stdev), _["Covariance"] = wrap(cov), _["AIC"] = 2*(totalnum-accumulate(KeepConstant.begin(), KeepConstant.end(), 0.0))-2*Ll_Total[0], _["BIC"] = (totalnum-accumulate(KeepConstant.begin(), KeepConstant.end(), 0.0))*log(mat_row)-2*Ll_Total[0], _["Parameter_Lists"] = para_list, _["Control_List"] = control_list, _["Converged"] = convgd, _["Status"] = "PASSED");
    }
    // returns a list of results
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
// [[Rcpp::export]]
List LogLik_Cox_PH_Omnibus_Log_Bound_CurveSearch(IntegerVector term_n, StringVector tform, NumericVector a_n, NumericMatrix& x_all, IntegerVector dfc, int fir, string modelform, double lr, List optim_para, int maxiter, int double_step, int halfmax, double epsilon, double abs_max, double dose_abs_max, double deriv_epsilon, NumericMatrix df_groups, NumericVector tu, int verbose, IntegerVector KeepConstant, int term_tot, string ties_method, int nthreads, NumericVector& Strata_vals, const VectorXd& cens_weight, List model_bool, const double gmix_theta, const IntegerVector gmix_term, const MatrixXd Lin_Sys, const VectorXd Lin_Res, double qchi, int para_number, int maxstep, double step_size) {
    //
    List temp_list = List::create(_["Status"] = "TEMP");  // used as a dummy return value for code checking
    // Time durations are measured from this point on in microseconds
    //
    // df0: covariate data
    // ntime: number of event times for Cox PH
    // totalnum: number of terms used
    //
    // ------------------------------------------------------------------------- // initialize
    const Map<MatrixXd> df0(as<Map<MatrixXd> >(x_all));
    int ntime = tu.size();
    int totalnum = term_n.size();
    int reqrdnum = totalnum - sum(KeepConstant);
    // ------------------------------------------------------------------------- // initialize
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
//    if (model_bool["single"]) {
//        if (verbose >= 1) {
//            Rcout << "non-derivative model calculation is not compatable with log-based bound calculation" << endl;
//        }
//        temp_list = List::create(_["Status"] = "FAILED_WITH_BAD_MODEL_SINGLE", _["LogLik"] = R_NaN);
//        return temp_list;
//    }
//    if (model_bool["gradient"]) {
//        if (verbose >= 1) {
//            Rcout << "gradient descent model calculation is not compatable with log-based bound calculation" << endl;
//        }
//        temp_list = List::create(_["Status"] = "FAILED_WITH_BAD_MODEL_GRADIENT", _["LogLik"] = R_NaN);
//        return temp_list;
//    }
    //
    // cout.precision: controls the number of significant digits printed
    // nthreads: number of threads used for parallel operations
    //
    Rcout.precision(7);  // forces higher precision numbers printed to terminal
    //
    // Lld_worst: stores the highest magnitude log-likelihood derivative
    // ---------------------------------------------
    // To Start, needs to seperate the derivative terms
    // ---------------------------------------------
    //
    // ------------------------------------------------------------------------- // initialize
    Map<VectorXd> beta_0(as<Map<VectorXd> >(a_n));
    MatrixXd T0;
    MatrixXd Td0;
    MatrixXd Tdd0;
    //
    MatrixXd Te;
    MatrixXd R;
    ColXd Rd;
    ColXd Rdd;
    //
    MatrixXd Dose;
    MatrixXd nonDose;
    MatrixXd nonDose_LIN;
    MatrixXd nonDose_PLIN;
    MatrixXd nonDose_LOGLIN;
    MatrixXd TTerm;
    double dint = 0.0;  // the amount of change used to calculate derivatives in threshold paramters
    double dslp = 0.0;
    ColXd RdR;
    ColXd RddR;
    // ------------------------------------------------------------------------- // initialize
    // ---------------------------------------------
    // To Start, needs to seperate the derivative terms
    // ---------------------------------------------
    //
    Cox_Refresh_R_TERM(totalnum, reqrdnum, term_tot, dint, dslp, dose_abs_max, abs_max, df0, T0, Td0, Tdd0, Te, R, Rd, Rdd, Dose, nonDose, TTerm, nonDose_LIN, nonDose_PLIN, nonDose_LOGLIN, RdR, RddR, model_bool);
    // ------------------------------------------------------------------------- // initialize
    IntegerMatrix RiskFail;
    vector<vector<int> > RiskPairs(ntime);
    vector<vector<vector<int> > > RiskPairs_Strata(ntime, vector<vector<int>>(Strata_vals.size()));
    const Map<MatrixXd> df_m(as<Map<MatrixXd> >(df_groups));
    // ------------------------------------------------------------------------- // initialize
    if (model_bool["strata"]) {
        RiskFail = IntegerMatrix(ntime, 2*Strata_vals.size());  // vector giving the event rows
        //
        // Creates matrices used to identify the event risk groups
        if (model_bool["cr"]) {
            Make_Groups_Strata_CR(ntime, df_m, RiskFail, RiskPairs_Strata, tu, nthreads, Strata_vals, cens_weight);
        } else {
            Make_Groups_Strata(ntime, df_m, RiskFail, RiskPairs_Strata, tu, nthreads, Strata_vals);
        }
    } else {
        RiskFail = IntegerMatrix(ntime, 2);  // vector giving the event rows
        //
        // Creates matrices used to identify the event risk groups
        if (model_bool["cr"]) {
            Make_Groups_CR(ntime, df_m, RiskFail, RiskPairs, tu, cens_weight, nthreads);
        } else {
            Make_Groups(ntime, df_m, RiskFail, RiskPairs, tu, nthreads);
        }
    }
    // ------------------------------------------------------------------------- // initialize
    MatrixXd Rls1;
    MatrixXd Lls1;
    MatrixXd Rls2;
    MatrixXd Rls3;
    MatrixXd Lls2;
    MatrixXd Lls3;
    vector<double> Ll(reqrdnum, 0.0);  // log-likelihood values
    vector<double> Lld(reqrdnum, 0.0);  // log-likelihood derivative values
    vector<double> Lldd(pow(reqrdnum, 2), 0.0);  // the second derivative matrix has room for every combination, but only the lower triangle is calculated initially
    // ------------------------------------------------------------------------- // initialize
    Cox_Refresh_R_SIDES(reqrdnum, ntime, Rls1, Rls2, Rls3, Lls1, Lls2, Lls3, Strata_vals, model_bool);
    // the log-likelihood is calculated in parallel over the risk groups
    vector <double> Ll_comp(2, Ll[0]);  // vector to compare values
    double abs_max0 = abs_max;
    double dose_abs_max0 = dose_abs_max;
    //
    vector<double> dbeta(totalnum, 0.0);
    //
    // --------------------------
    // always starts from initial guess
    // --------------------------
    vector<double> beta_c(totalnum, 0.0);
    vector<double> beta_a(totalnum, 0.0);
    vector<double> beta_best(totalnum, 0.0);
    vector<double> beta_peak(totalnum, 0.0);
    vector<double> beta_p(totalnum, 0.0);
    VectorXd::Map(&beta_p[0], beta_0.size()) = beta_0;  // stores previous parameters
    VectorXd::Map(&beta_c[0], beta_0.size()) = beta_0;  // stores current parameters
    VectorXd::Map(&beta_a[0], beta_0.size()) = beta_0;  // stores a refrence value for parameters
    VectorXd::Map(&beta_best[0], beta_0.size()) = beta_0;  // stores the best parameters
    VectorXd::Map(&beta_peak[0], beta_0.size()) = beta_0;  // stores the best parameters
    int iter_stop  = 0;  // tracks if the iterations should be stopped for convergence
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
    abs_max = abs_max0;
    dose_abs_max = dose_abs_max0;
    //
    for (int i = 0; i < beta_0.size(); i++) {
        beta_0[i] = a_n[i];
    }
    Cox_Term_Risk_Calc(modelform, tform, term_n, totalnum, fir, dfc, term_tot, T0, Td0, Tdd0, Te, R, Rd, Rdd, Dose, nonDose, beta_0, df0, dose_abs_max, abs_max, TTerm, nonDose_LIN, nonDose_PLIN, nonDose_LOGLIN, RdR, RddR, nthreads, KeepConstant, verbose, model_bool, gmix_theta, gmix_term);
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
    double Lpeak = Ll[0];
    // bool upper = true;
    // int half_check = 0;
    // bool trouble = false;
    IntegerVector KeepConstant_trouble (totalnum);
    for (int i = 0;  i < totalnum; i++) {
        KeepConstant_trouble[i] = KeepConstant[i];
    }
    KeepConstant_trouble[para_number] = 1;
    //
    vector<double> limits(2, 0.0);
    vector<bool>   limit_hit(2, FALSE);
    vector<bool>   limit_converged(2, FALSE);
    vector<double> ll_final(2, 0.0);
    vector<double> width_final(2,0.0);
    vector<int>    step_final(2,0.0);
    List res_list;
    //
    if (verbose >= 4) {
        Rcout << "C++ Note: STARTING Upper Bound" << endl;
    }
    // upper = true;
    // int step = -1;
    // bool iter_continue = true;
    // double max_change = 100;
    // double deriv_max = 100;
    bool convgd = false;
    ///
    // variables added for log loop code
    ///
    VectorXd s_weights(1);
    // int bound_val = 1;
    // We need the values reserved for the upper, middle, lower estimates and scores
    vector<double> beta_L(totalnum, 0.0);
    vector<double> beta_M(totalnum, 0.0);
    vector<double> beta_H(totalnum, 0.0);
    double L_L= 0.0;
    double L_M= 0.0;
    double L_H= 0.0;
    NumericVector reg_beta(totalnum);
    // First we need to establish the first interval estimates
    for (int ij = 0; ij < totalnum; ij++) {
        beta_L[ij] = beta_peak[ij];
        beta_H[ij] = beta_peak[ij];
    }
    L_L = Lpeak;
    bool loop_check = true;
    double temp_step = step_size;
    NumericVector temp_L(1);
    List reg_out;
    while ((loop_check) && (temp_step > 1e-3)){
        // assign new high point
        beta_H[para_number] = beta_L[para_number] + temp_step;
        for (int ij = 0; ij < totalnum; ij++) {
            beta_0[ij] = beta_H[ij];
        }
        temp_step = temp_step * 0.5;
        //
        reg_out = Cox_Full_Run(reqrdnum, ntime, tform, RiskFail,  RiskPairs, RiskPairs_Strata, totalnum, fir, R, Rd, Rdd, Rls1, Rls2, Rls3, Lls1, Lls2, Lls3, cens_weight, Strata_vals, beta_0, RdR, RddR, Ll, Lld, Lldd, nthreads, KeepConstant, ties_method, verbose, model_bool, iter_stop, term_tot, dint, dslp, dose_abs_max, abs_max, df0, T0, Td0, Tdd0, Te, Dose, nonDose, TTerm, nonDose_LIN, nonDose_PLIN, nonDose_LOGLIN, modelform, gmix_theta, gmix_term, convgd, lr, optim_para, maxiter, double_step, Lin_Sys, Lin_Res, term_n, dfc, halfmax, epsilon, deriv_epsilon);
        temp_L[0] = reg_out["LogLik"];
        if (!is_nan(temp_L)[0]){
        // if (reg_status == "PASSED"){
            loop_check = false;
        }
    }
    if (loop_check){
        limit_hit[1] = true;
        limits[1] = 0;
        ll_final[1] = 0;
        limit_converged[1] = false;
    } else {
        // Now we can run the actual algorithm
        reg_beta = reg_out["beta_0"];
        for (int ij = 0; ij < totalnum; ij++) {
            beta_H[ij] = reg_beta[ij];
        }
        L_H = reg_out["LogLik"];
        // now set the mid point value
        for (int ij = 0; ij < totalnum; ij++) {
            beta_M[ij] = (beta_H[ij] + beta_L[ij])/2;
            beta_0[ij] = beta_M[ij];
        }
        reg_out = Cox_Full_Run(reqrdnum, ntime, tform, RiskFail,  RiskPairs, RiskPairs_Strata, totalnum, fir, R, Rd, Rdd, Rls1, Rls2, Rls3, Lls1, Lls2, Lls3, cens_weight, Strata_vals, beta_0, RdR, RddR, Ll, Lld, Lldd, nthreads, KeepConstant, ties_method, verbose, model_bool, iter_stop, term_tot, dint, dslp, dose_abs_max, abs_max, df0, T0, Td0, Tdd0, Te, Dose, nonDose, TTerm, nonDose_LIN, nonDose_PLIN, nonDose_LOGLIN, modelform, gmix_theta, gmix_term, convgd, lr, optim_para, maxiter, double_step, Lin_Sys, Lin_Res, term_n, dfc, halfmax, epsilon, deriv_epsilon);
        reg_beta = reg_out["beta_0"];
        for (int ij = 0; ij < totalnum; ij++) {
            beta_M[ij] = reg_beta[ij];
        }
        L_M = reg_out["LogLik"];
        //
        int step = 0;
        // now run the bisection until stopping point
        // while ((step < step_limit) & (abs(beta_low[para_num] - beta_high[para_num]) > control$epsilon) & (!Limit_Hit[2])) {
        while ((step < maxstep) && (abs(beta_L[para_number] - beta_H[para_number]) > epsilon) && (! limit_hit[1])){
            step = step + 1;
            if (L_L < Lstar) {
                throw invalid_argument("The lower estimate is too high?");
            } else if (L_H < Lstar){
                // midpoint is in between the two
                if (L_M < Lstar){
                    // the mid point is past the optimum
                    for (int ij = 0; ij < totalnum; ij++) {
                        beta_H[ij] = beta_M[ij];
                    }
                    L_H = L_M;
                } else {
                    // the mid point is before the optimum
                    for (int ij = 0; ij < totalnum; ij++) {
                        beta_L[ij] = beta_M[ij];
                    }
                    L_L = L_M;
                }
            } else if (L_M < Lstar){
                for (int ij = 0; ij < totalnum; ij++) {
                    beta_H[ij] = beta_M[ij];
                }
                L_H = L_M;
            } else {
                // the upper estimate needs to be shifted up
                for (int ij = 0; ij < totalnum; ij++) {
                    beta_L[ij] = beta_H[ij];
                }
                L_L = L_H;
                // check new high point
                loop_check = true;
                temp_step = step_size;
                while ((loop_check) && (temp_step > 1e-3)){
                    // assign new high point
                    beta_H[para_number] = beta_L[para_number] + temp_step;
                    for (int ij = 0; ij < totalnum; ij++) {
                        beta_0[ij] = beta_H[ij];
                    }
                    temp_step = temp_step * 0.5;
                    //
                    reg_out = Cox_Full_Run(reqrdnum, ntime, tform, RiskFail,  RiskPairs, RiskPairs_Strata, totalnum, fir, R, Rd, Rdd, Rls1, Rls2, Rls3, Lls1, Lls2, Lls3, cens_weight, Strata_vals, beta_0, RdR, RddR, Ll, Lld, Lldd, nthreads, KeepConstant, ties_method, verbose, model_bool, iter_stop, term_tot, dint, dslp, dose_abs_max, abs_max, df0, T0, Td0, Tdd0, Te, Dose, nonDose, TTerm, nonDose_LIN, nonDose_PLIN, nonDose_LOGLIN, modelform, gmix_theta, gmix_term, convgd, lr, optim_para, maxiter, double_step, Lin_Sys, Lin_Res, term_n, dfc, halfmax, epsilon, deriv_epsilon);
                    temp_L[0] = reg_out["LogLik"];
                    if (!is_nan(temp_L)[0]){
                        loop_check = false;
                    }
                }
                if (loop_check){
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
            if (!limit_hit[1]){
                // now set the mid point value
                for (int ij = 0; ij < totalnum; ij++) {
                    beta_M[ij] = (beta_H[ij] + beta_L[ij])/2;
                    beta_0[ij] = beta_M[ij];
                }
                reg_out = Cox_Full_Run(reqrdnum, ntime, tform, RiskFail,  RiskPairs, RiskPairs_Strata, totalnum, fir, R, Rd, Rdd, Rls1, Rls2, Rls3, Lls1, Lls2, Lls3, cens_weight, Strata_vals, beta_0, RdR, RddR, Ll, Lld, Lldd, nthreads, KeepConstant, ties_method, verbose, model_bool, iter_stop, term_tot, dint, dslp, dose_abs_max, abs_max, df0, T0, Td0, Tdd0, Te, Dose, nonDose, TTerm, nonDose_LIN, nonDose_PLIN, nonDose_LOGLIN, modelform, gmix_theta, gmix_term, convgd, lr, optim_para, maxiter, double_step, Lin_Sys, Lin_Res, term_n, dfc, halfmax, epsilon, deriv_epsilon);
                reg_beta = reg_out["beta_0"];
                for (int ij = 0; ij < totalnum; ij++) {
                    beta_M[ij] = reg_beta[ij];
                }
                L_M = reg_out["LogLik"];
            }
        }
        if ((abs(beta_L[para_number] - beta_H[para_number]) < epsilon) && (! limit_hit[1])){
            limit_converged[1] = true;
        }
        limits[1] = beta_M[para_number];
        ll_final[1] = L_M;
        width_final[1] = abs(beta_L[para_number] - beta_H[para_number]);
        step_final[1] = step;
    }
    // upper limit found, now solve lower limit
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
    while ((loop_check) && (temp_step > 1e-3)){
        // assign new high point
        beta_L[para_number] = beta_H[para_number] - temp_step;
        for (int ij = 0; ij < totalnum; ij++) {
            beta_0[ij] = beta_L[ij];
        }
        temp_step = temp_step * 0.5;
        //
        reg_out = Cox_Full_Run(reqrdnum, ntime, tform, RiskFail,  RiskPairs, RiskPairs_Strata, totalnum, fir, R, Rd, Rdd, Rls1, Rls2, Rls3, Lls1, Lls2, Lls3, cens_weight, Strata_vals, beta_0, RdR, RddR, Ll, Lld, Lldd, nthreads, KeepConstant, ties_method, verbose, model_bool, iter_stop, term_tot, dint, dslp, dose_abs_max, abs_max, df0, T0, Td0, Tdd0, Te, Dose, nonDose, TTerm, nonDose_LIN, nonDose_PLIN, nonDose_LOGLIN, modelform, gmix_theta, gmix_term, convgd, lr, optim_para, maxiter, double_step, Lin_Sys, Lin_Res, term_n, dfc, halfmax, epsilon, deriv_epsilon);
        temp_L[0] = reg_out["LogLik"];
        if (!is_nan(temp_L)[0]){
            loop_check = false;
        }
    }
    if (loop_check){
        limit_hit[0] = true;
        limits[0] = 0;
        ll_final[0] = 0;
        limit_converged[0] = false;
    } else {
        // Now we can run the actual algorithm
        reg_beta = reg_out["beta_0"];
        for (int ij = 0; ij < totalnum; ij++) {
            beta_L[ij] = reg_beta[ij];
        }
        L_L = reg_out["LogLik"];
        // now set the mid point value
        for (int ij = 0; ij < totalnum; ij++) {
            beta_M[ij] = (beta_H[ij] + beta_L[ij])/2;
            beta_0[ij] = beta_M[ij];
        }
        reg_out = Cox_Full_Run(reqrdnum, ntime, tform, RiskFail,  RiskPairs, RiskPairs_Strata, totalnum, fir, R, Rd, Rdd, Rls1, Rls2, Rls3, Lls1, Lls2, Lls3, cens_weight, Strata_vals, beta_0, RdR, RddR, Ll, Lld, Lldd, nthreads, KeepConstant, ties_method, verbose, model_bool, iter_stop, term_tot, dint, dslp, dose_abs_max, abs_max, df0, T0, Td0, Tdd0, Te, Dose, nonDose, TTerm, nonDose_LIN, nonDose_PLIN, nonDose_LOGLIN, modelform, gmix_theta, gmix_term, convgd, lr, optim_para, maxiter, double_step, Lin_Sys, Lin_Res, term_n, dfc, halfmax, epsilon, deriv_epsilon);
        reg_beta = reg_out["beta_0"];
        for (int ij = 0; ij < totalnum; ij++) {
            beta_M[ij] = reg_beta[ij];
        }
        L_M = reg_out["LogLik"];
        //
        int step = 0;
        // now run the bisection until stopping point
        // while ((step < step_limit) & (abs(beta_low[para_num] - beta_high[para_num]) > control$epsilon) & (!Limit_Hit[2])) {
        while ((step < maxstep) && (abs(beta_L[para_number] - beta_H[para_number]) > epsilon) && (! limit_hit[0])){
            step = step + 1;
            if (L_H < Lstar) {
                throw invalid_argument("The upper estimate is too high?");
            } else if (L_L < Lstar){
                // midpoint is in between the two
                if (L_M > Lstar){
                    // the mid point is past the optimum
                    for (int ij = 0; ij < totalnum; ij++) {
                        beta_H[ij] = beta_M[ij];
                    }
                    L_H = L_M;
                } else {
                    // the mid point is before the optimum
                    for (int ij = 0; ij < totalnum; ij++) {
                        beta_L[ij] = beta_M[ij];
                    }
                    L_L = L_M;
                }
            } else if (L_M < Lstar){
                for (int ij = 0; ij < totalnum; ij++) {
                    beta_L[ij] = beta_M[ij];
                }
                L_L = L_M;
            } else {
                // the lower estimate needs to be shifted down
                for (int ij = 0; ij < totalnum; ij++) {
                    beta_H[ij] = beta_L[ij];
                }
                L_H = L_L;
                // check new high point
                loop_check = true;
                temp_step = step_size;
                while ((loop_check) && (temp_step > 1e-3)){
                    // assign new high point
                    beta_L[para_number] = beta_H[para_number] - temp_step;
                    for (int ij = 0; ij < totalnum; ij++) {
                        beta_0[ij] = beta_L[ij];
                    }
                    temp_step = temp_step * 0.5;
                    //
                    reg_out = Cox_Full_Run(reqrdnum, ntime, tform, RiskFail,  RiskPairs, RiskPairs_Strata, totalnum, fir, R, Rd, Rdd, Rls1, Rls2, Rls3, Lls1, Lls2, Lls3, cens_weight, Strata_vals, beta_0, RdR, RddR, Ll, Lld, Lldd, nthreads, KeepConstant, ties_method, verbose, model_bool, iter_stop, term_tot, dint, dslp, dose_abs_max, abs_max, df0, T0, Td0, Tdd0, Te, Dose, nonDose, TTerm, nonDose_LIN, nonDose_PLIN, nonDose_LOGLIN, modelform, gmix_theta, gmix_term, convgd, lr, optim_para, maxiter, double_step, Lin_Sys, Lin_Res, term_n, dfc, halfmax, epsilon, deriv_epsilon);
                    temp_L[0] = reg_out["LogLik"];
                    if (!is_nan(temp_L)[0]){
                        loop_check = false;
                    }
                }
                if (loop_check){
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
            if (!limit_hit[0]){
                // now set the mid point value
                for (int ij = 0; ij < totalnum; ij++) {
                    beta_M[ij] = (beta_H[ij] + beta_L[ij])/2;
                    beta_0[ij] = beta_M[ij];
                }
                reg_out = Cox_Full_Run(reqrdnum, ntime, tform, RiskFail,  RiskPairs, RiskPairs_Strata, totalnum, fir, R, Rd, Rdd, Rls1, Rls2, Rls3, Lls1, Lls2, Lls3, cens_weight, Strata_vals, beta_0, RdR, RddR, Ll, Lld, Lldd, nthreads, KeepConstant, ties_method, verbose, model_bool, iter_stop, term_tot, dint, dslp, dose_abs_max, abs_max, df0, T0, Td0, Tdd0, Te, Dose, nonDose, TTerm, nonDose_LIN, nonDose_PLIN, nonDose_LOGLIN, modelform, gmix_theta, gmix_term, convgd, lr, optim_para, maxiter, double_step, Lin_Sys, Lin_Res, term_n, dfc, halfmax, epsilon, deriv_epsilon);
                reg_beta = reg_out["beta_0"];
                for (int ij = 0; ij < totalnum; ij++) {
                    beta_M[ij] = reg_beta[ij];
                }
                L_M = reg_out["LogLik"];
            }
        }
        if ((abs(beta_L[para_number] - beta_H[para_number]) < epsilon) && (! limit_hit[0])){
            limit_converged[0] = true;
        }
        limits[0] = beta_M[para_number];
        ll_final[0] = L_M;
        width_final[0] = abs(beta_L[para_number] - beta_H[para_number]);
        step_final[0] = step;
    }
    res_list = List::create(_["Parameter_Limits"] = wrap(limits), _["Negative_Risk_Limit_Hit"] = wrap(limit_hit), _["Likelihood_Boundary"] = wrap(ll_final), _["Likelihood_Goal"] = wrap(Lstar), _["Limit_Converged"] = wrap(limit_converged), _["Final_Window_Width"] = wrap(width_final), _["Final_Step"] = wrap(step_final), _["Status"] = "PASSED");
    // returns a list of results
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
// [[Rcpp::export]]
List LogLik_Poisson_Omnibus_Log_Bound_CurveSearch(const MatrixXd& PyrC, const MatrixXd& dfs, IntegerVector term_n, StringVector tform, NumericVector a_n, NumericMatrix& x_all, IntegerVector dfc, int fir, string modelform, double lr, List optim_para, int maxiter, int doublestep, int halfmax, double epsilon, double abs_max, double dose_abs_max, double deriv_epsilon, int verbose, IntegerVector KeepConstant, int term_tot, int nthreads, List model_bool, const double gmix_theta, const IntegerVector gmix_term, const MatrixXd Lin_Sys, const VectorXd Lin_Res, double qchi, int para_number, int maxstep, double step_size) {
    //
    List temp_list = List::create(_["Status"] = "TEMP");  // used as a dummy return value for code checking
    // Time durations are measured from this point on in microseconds
    //
    // df0: covariate data
    // ntime: number of event times for Cox PH
    // totalnum: number of terms used
    //
    // ------------------------------------------------------------------------- // initialize
    const Map<MatrixXd> df0(as<Map<MatrixXd> >(x_all));
    const int mat_row = df0.rows();
    // int ntime = tu.size();
    int totalnum = term_n.size();
    int reqrdnum = totalnum - sum(KeepConstant);
    // ------------------------------------------------------------------------- // initialize
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
    // cout.precision: controls the number of significant digits printed
    // nthreads: number of threads used for parallel operations
    //
    Rcout.precision(7);  // forces higher precision numbers printed to terminal
    //
    // Lld_worst: stores the highest magnitude log-likelihood derivative
    // ---------------------------------------------
    // To Start, needs to seperate the derivative terms
    // ---------------------------------------------
    //
    // ------------------------------------------------------------------------- // initialize
    Map<VectorXd> beta_0(as<Map<VectorXd> >(a_n));
    VectorXd beta_max = beta_0;
    MatrixXd T0;
    MatrixXd Td0;
    MatrixXd Tdd0;
    //
    MatrixXd Te;
    MatrixXd R;
    ColXd Rd;
    ColXd Rdd;
    //
    MatrixXd Dose;
    MatrixXd nonDose;
    MatrixXd nonDose_LIN;
    MatrixXd nonDose_PLIN;
    MatrixXd nonDose_LOGLIN;
    MatrixXd TTerm;
    double dint = 0.0;  // the amount of change used to calculate derivatives in threshold paramters
    double dslp = 0.0;
    ColXd RdR;
    ColXd RddR;
    // ------------------------------------------------------------------------- // initialize
    // ---------------------------------------------
    // To Start, needs to seperate the derivative terms
    // ---------------------------------------------
    //
    T0 = MatrixXd::Zero(mat_row, totalnum);  // preallocates matrix for Term column
    Te = MatrixXd::Zero(mat_row, 1);  // preallocates matrix for column terms used for temporary storage
    R = MatrixXd::Zero(mat_row, 1);  // preallocates matrix for Risks
    //
    Dose = MatrixXd::Constant(mat_row, term_tot, 0.0);  // matrix of the total dose term values
    nonDose = MatrixXd::Constant(mat_row, term_tot, 1.0);  // matrix of the total non-dose term values
    nonDose_LIN = MatrixXd::Constant(mat_row, term_tot, 0.0);  // matrix of Linear subterm values
    nonDose_PLIN = MatrixXd::Constant(mat_row, term_tot, 1.0);  // matrix of Loglinear subterm values
    nonDose_LOGLIN = MatrixXd::Constant(mat_row, term_tot, 1.0);  // matrix of Product linear subterm values
    TTerm = MatrixXd::Zero(mat_row, term_tot);  // matrix of term values
    if (model_bool["single"]) {
    } else {
        Td0 = MatrixXd::Zero(mat_row, reqrdnum);  // preallocates matrix for Term derivative columns
        Tdd0 = MatrixXd::Zero(mat_row, reqrdnum*(reqrdnum + 1)/2);  // preallocates matrix for Term second derivative columns
        Rd = MatrixXd::Zero(mat_row, reqrdnum);  // preallocates matrix for Risk derivatives
        Rdd = MatrixXd::Zero(mat_row, reqrdnum*(reqrdnum + 1)/2);  // preallocates matrix for Risk second derivatives
        dint = dose_abs_max;  // the amount of change used to calculate derivatives in threshold paramters
        dslp = abs_max;
        RdR = MatrixXd::Zero(mat_row, reqrdnum);  // preallocates matrix for Risk to derivative ratios
        RddR = MatrixXd::Zero(mat_row, reqrdnum*(reqrdnum + 1)/2);  // preallocates matrix for Risk to second derivative ratios
    }
    VectorXd s_weights;
    if (model_bool["strata"]) {
        s_weights = VectorXd::Zero(mat_row);
        Gen_Strat_Weight(modelform, dfs, PyrC, s_weights, nthreads, tform, term_n, term_tot);
    }
    // ------------------------------------------------------------------------- // initialize
    MatrixXd Rls1;
    MatrixXd Lls1;
    MatrixXd Rls2;
    MatrixXd Rls3;
    MatrixXd Lls2;
    MatrixXd Lls3;
    vector<double> Ll(reqrdnum, 0.0);  // log-likelihood values
    vector<double> Lld(reqrdnum, 0.0);  // log-likelihood derivative values
    vector<double> Lldd(pow(reqrdnum, 2), 0.0);  // the second derivative matrix has room for every combination, but only the lower triangle is calculated initially
    MatrixXd dev_temp = MatrixXd::Zero(PyrC.rows(), 2);
    double dev = 0;
    // ------------------------------------------------------------------------- // initialize
    // the log-likelihood is calculated in parallel over the risk groups
    vector <double> Ll_comp(2, Ll[0]);  // vector to compare values
    double abs_max0 = abs_max;
    double dose_abs_max0 = dose_abs_max;
    //
    vector<double> dbeta(totalnum, 0.0);
    vector<double> dbeta_start(totalnum, 0.0);
    //
    // --------------------------
    // always starts from initial guess
    // --------------------------
    vector<double> beta_c(totalnum, 0.0);
    vector<double> beta_a(totalnum, 0.0);
    vector<double> beta_best(totalnum, 0.0);
    vector<double> beta_p(totalnum, 0.0);
    vector<double> beta_peak(totalnum, 0.0);
    VectorXd::Map(&beta_p[0], beta_0.size()) = beta_0;  // stores previous parameters
    VectorXd::Map(&beta_c[0], beta_0.size()) = beta_0;  // stores current parameters
    VectorXd::Map(&beta_a[0], beta_0.size()) = beta_0;  // stores a refrence value for parameters
    VectorXd::Map(&beta_best[0], beta_0.size()) = beta_0;  // stores the best parameters
    VectorXd::Map(&beta_peak[0], beta_0.size()) = beta_0;  // stores the peak parameters
    int iter_stop  = 0;  // tracks if the iterations should be stopped for convergence
    vector<double> beta_abs_best(totalnum, 0.0);
    fill(Ll.begin(), Ll.end(), 0.0);
    fill(Lld.begin(), Lld.end(), 0.0);
    fill(Lldd.begin(), Lldd.end(), 0.0);
    beta_p = beta_peak;  //
    beta_a = beta_peak;  //
    beta_c = beta_peak;  //
    abs_max = abs_max0;
    dose_abs_max = dose_abs_max0;
    for (int i = 0; i < beta_0.size(); i++) {
        beta_0[i] = a_n[i];
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
    double Lpeak = Ll[0];
    // bool upper = true;
    // int half_check = 0;
    // bool trouble = false;
    IntegerVector KeepConstant_trouble (totalnum);
    for (int i = 0;  i < totalnum; i++) {
        KeepConstant_trouble[i] = KeepConstant[i];
    }
    KeepConstant_trouble[para_number] = 1;
    //
    vector<double> limits(2, 0.0);
    vector<bool>   limit_hit(2, FALSE);
    vector<bool>   limit_converged(2, FALSE);
    vector<double> ll_final(2, 0.0);
    vector<double> width_final(2,0.0);
    vector<int>    step_final(2,0.0);
    List res_list;
    //
    if (verbose >= 4) {
        Rcout << "C++ Note: STARTING Upper Bound" << endl;
    }
    // upper = true;
    // int step = -1;
    // bool iter_continue = true;
    // double max_change = 100;
    // double deriv_max = 100;
    bool convgd = false;
    int double_step = 1;
    ///
    // variables added for log loop code
    ///
    // VectorXd s_weights(1);
    // int bound_val = 1;
    // We need the values reserved for the upper, middle, lower estimates and scores
    vector<double> beta_L(totalnum, 0.0);
    vector<double> beta_M(totalnum, 0.0);
    vector<double> beta_H(totalnum, 0.0);
    double L_L= 0.0;
    double L_M= 0.0;
    double L_H= 0.0;
    NumericVector reg_beta(totalnum);
    // First we need to establish the first interval estimates
    for (int ij = 0; ij < totalnum; ij++) {
        beta_L[ij] = beta_peak[ij];
        beta_H[ij] = beta_peak[ij];
    }
    L_L = Lpeak;
    bool loop_check = true;
    double temp_step = step_size;
    NumericVector temp_L(1);
    List reg_out;
    while ((loop_check) && (temp_step > 1e-3)){
        // assign new high point
        beta_H[para_number] = beta_L[para_number] + temp_step;
        for (int ij = 0; ij < totalnum; ij++) {
            beta_0[ij] = beta_H[ij];
        }
        temp_step = temp_step * 0.5;
        //
        reg_out = Pois_Full_Run(PyrC, reqrdnum, tform, totalnum, fir, R, Rd, Rdd, s_weights, beta_0, RdR, RddR, Ll, Lld, Lldd, nthreads, KeepConstant, verbose, model_bool, iter_stop, term_tot, dint, dslp, dose_abs_max, abs_max, df0, T0, Td0, Tdd0, Te, Dose, nonDose, TTerm, nonDose_LIN, nonDose_PLIN, nonDose_LOGLIN, modelform, gmix_theta, gmix_term, convgd, lr, optim_para, maxiter, double_step, Lin_Sys, Lin_Res, term_n, dfc, halfmax, epsilon, deriv_epsilon);
        temp_L[0] = reg_out["LogLik"];
        if (!is_nan(temp_L)[0]){
        // if (reg_status == "PASSED"){
            loop_check = false;
        }
    }
    if (loop_check){
        limit_hit[1] = true;
        limits[1] = 0;
        ll_final[1] = 0;
        limit_converged[1] = false;
    } else {
        // Now we can run the actual algorithm
        reg_beta = reg_out["beta_0"];
        for (int ij = 0; ij < totalnum; ij++) {
            beta_H[ij] = reg_beta[ij];
        }
        L_H = reg_out["LogLik"];
        // now set the mid point value
        for (int ij = 0; ij < totalnum; ij++) {
            beta_M[ij] = (beta_H[ij] + beta_L[ij])/2;
            beta_0[ij] = beta_M[ij];
        }
        reg_out = Pois_Full_Run(PyrC, reqrdnum, tform, totalnum, fir, R, Rd, Rdd, s_weights, beta_0, RdR, RddR, Ll, Lld, Lldd, nthreads, KeepConstant, verbose, model_bool, iter_stop, term_tot, dint, dslp, dose_abs_max, abs_max, df0, T0, Td0, Tdd0, Te, Dose, nonDose, TTerm, nonDose_LIN, nonDose_PLIN, nonDose_LOGLIN, modelform, gmix_theta, gmix_term, convgd, lr, optim_para, maxiter, double_step, Lin_Sys, Lin_Res, term_n, dfc, halfmax, epsilon, deriv_epsilon);
        reg_beta = reg_out["beta_0"];
        for (int ij = 0; ij < totalnum; ij++) {
            beta_M[ij] = reg_beta[ij];
        }
        L_M = reg_out["LogLik"];
        //
        int step = 0;
        // now run the bisection until stopping point
        // while ((step < step_limit) & (abs(beta_low[para_num] - beta_high[para_num]) > control$epsilon) & (!Limit_Hit[2])) {
        while ((step < maxstep) && (abs(beta_L[para_number] - beta_H[para_number]) > epsilon) && (! limit_hit[1])){
            step = step + 1;
            if (L_L < Lstar) {
                throw invalid_argument("The lower estimate is too high?");
            } else if (L_H < Lstar){
                // midpoint is in between the two
                if (L_M < Lstar){
                    // the mid point is past the optimum
                    for (int ij = 0; ij < totalnum; ij++) {
                        beta_H[ij] = beta_M[ij];
                    }
                    L_H = L_M;
                } else {
                    // the mid point is before the optimum
                    for (int ij = 0; ij < totalnum; ij++) {
                        beta_L[ij] = beta_M[ij];
                    }
                    L_L = L_M;
                }
            } else if (L_M < Lstar){
                for (int ij = 0; ij < totalnum; ij++) {
                    beta_H[ij] = beta_M[ij];
                }
                L_H = L_M;
            } else {
                // the upper estimate needs to be shifted up
                for (int ij = 0; ij < totalnum; ij++) {
                    beta_L[ij] = beta_H[ij];
                }
                L_L = L_H;
                // check new high point
                loop_check = true;
                temp_step = step_size;
                while ((loop_check) && (temp_step > 1e-3)){
                    // assign new high point
                    beta_H[para_number] = beta_L[para_number] + temp_step;
                    for (int ij = 0; ij < totalnum; ij++) {
                        beta_0[ij] = beta_H[ij];
                    }
                    temp_step = temp_step * 0.5;
                    //
                    reg_out = Pois_Full_Run(PyrC, reqrdnum, tform, totalnum, fir, R, Rd, Rdd, s_weights, beta_0, RdR, RddR, Ll, Lld, Lldd, nthreads, KeepConstant, verbose, model_bool, iter_stop, term_tot, dint, dslp, dose_abs_max, abs_max, df0, T0, Td0, Tdd0, Te, Dose, nonDose, TTerm, nonDose_LIN, nonDose_PLIN, nonDose_LOGLIN, modelform, gmix_theta, gmix_term, convgd, lr, optim_para, maxiter, double_step, Lin_Sys, Lin_Res, term_n, dfc, halfmax, epsilon, deriv_epsilon);
                    temp_L[0] = reg_out["LogLik"];
                    if (!is_nan(temp_L)[0]){
                        loop_check = false;
                    }
                }
                if (loop_check){
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
            if (!limit_hit[1]){
                // now set the mid point value
                for (int ij = 0; ij < totalnum; ij++) {
                    beta_M[ij] = (beta_H[ij] + beta_L[ij])/2;
                    beta_0[ij] = beta_M[ij];
                }
                reg_out = Pois_Full_Run(PyrC, reqrdnum, tform, totalnum, fir, R, Rd, Rdd, s_weights, beta_0, RdR, RddR, Ll, Lld, Lldd, nthreads, KeepConstant, verbose, model_bool, iter_stop, term_tot, dint, dslp, dose_abs_max, abs_max, df0, T0, Td0, Tdd0, Te, Dose, nonDose, TTerm, nonDose_LIN, nonDose_PLIN, nonDose_LOGLIN, modelform, gmix_theta, gmix_term, convgd, lr, optim_para, maxiter, double_step, Lin_Sys, Lin_Res, term_n, dfc, halfmax, epsilon, deriv_epsilon);
                reg_beta = reg_out["beta_0"];
                for (int ij = 0; ij < totalnum; ij++) {
                    beta_M[ij] = reg_beta[ij];
                }
                L_M = reg_out["LogLik"];
            }
        }
        if ((abs(beta_L[para_number] - beta_H[para_number]) < epsilon) && (! limit_hit[1])){
            limit_converged[1] = true;
        }
        limits[1] = beta_M[para_number];
        ll_final[1] = L_M;
        width_final[1] = abs(beta_L[para_number] - beta_H[para_number]);
        step_final[1] = step;
    }
    // upper limit found, now solve lower limit
    for (int ij = 0; ij < totalnum; ij++) {
        beta_L[ij] = beta_peak[ij];
        beta_H[ij] = beta_peak[ij];
    }
    L_H = Lpeak;
    loop_check = true;
    temp_step = step_size;
    while ((loop_check) && (temp_step > 1e-3)){
        // assign new high point
        beta_L[para_number] = beta_H[para_number] - temp_step;
        for (int ij = 0; ij < totalnum; ij++) {
            beta_0[ij] = beta_L[ij];
        }
        temp_step = temp_step * 0.5;
        //
        reg_out = Pois_Full_Run(PyrC, reqrdnum, tform, totalnum, fir, R, Rd, Rdd, s_weights, beta_0, RdR, RddR, Ll, Lld, Lldd, nthreads, KeepConstant, verbose, model_bool, iter_stop, term_tot, dint, dslp, dose_abs_max, abs_max, df0, T0, Td0, Tdd0, Te, Dose, nonDose, TTerm, nonDose_LIN, nonDose_PLIN, nonDose_LOGLIN, modelform, gmix_theta, gmix_term, convgd, lr, optim_para, maxiter, double_step, Lin_Sys, Lin_Res, term_n, dfc, halfmax, epsilon, deriv_epsilon);
        temp_L[0] = reg_out["LogLik"];
        if (!is_nan(temp_L)[0]){
            loop_check = false;
        }
    }
    if (loop_check){
        limit_hit[0] = true;
        limits[0] = 0;
        ll_final[0] = 0;
        limit_converged[0] = false;
    } else {
        // Now we can run the actual algorithm
        reg_beta = reg_out["beta_0"];
        for (int ij = 0; ij < totalnum; ij++) {
            beta_L[ij] = reg_beta[ij];
        }
        L_L = reg_out["LogLik"];
        // now set the mid point value
        for (int ij = 0; ij < totalnum; ij++) {
            beta_M[ij] = (beta_H[ij] + beta_L[ij])/2;
            beta_0[ij] = beta_M[ij];
        }
        reg_out = Pois_Full_Run(PyrC, reqrdnum, tform, totalnum, fir, R, Rd, Rdd, s_weights, beta_0, RdR, RddR, Ll, Lld, Lldd, nthreads, KeepConstant, verbose, model_bool, iter_stop, term_tot, dint, dslp, dose_abs_max, abs_max, df0, T0, Td0, Tdd0, Te, Dose, nonDose, TTerm, nonDose_LIN, nonDose_PLIN, nonDose_LOGLIN, modelform, gmix_theta, gmix_term, convgd, lr, optim_para, maxiter, double_step, Lin_Sys, Lin_Res, term_n, dfc, halfmax, epsilon, deriv_epsilon);
        reg_beta = reg_out["beta_0"];
        for (int ij = 0; ij < totalnum; ij++) {
            beta_M[ij] = reg_beta[ij];
        }
        L_M = reg_out["LogLik"];
        //
        int step = 0;
        // now run the bisection until stopping point
        // while ((step < step_limit) & (abs(beta_low[para_num] - beta_high[para_num]) > control$epsilon) & (!Limit_Hit[2])) {
        while ((step < maxstep) && (abs(beta_L[para_number] - beta_H[para_number]) > epsilon) && (! limit_hit[0])){
            step = step + 1;
            if (L_H < Lstar) {
                throw invalid_argument("The upper estimate is too high?");
            } else if (L_L < Lstar){
                // midpoint is in between the two
                if (L_M > Lstar){
                    // the mid point is past the optimum
                    for (int ij = 0; ij < totalnum; ij++) {
                        beta_H[ij] = beta_M[ij];
                    }
                    L_H = L_M;
                } else {
                    // the mid point is before the optimum
                    for (int ij = 0; ij < totalnum; ij++) {
                        beta_L[ij] = beta_M[ij];
                    }
                    L_L = L_M;
                }
            } else if (L_M < Lstar){
                for (int ij = 0; ij < totalnum; ij++) {
                    beta_L[ij] = beta_M[ij];
                }
                L_L = L_M;
            } else {
                // the lower estimate needs to be shifted down
                for (int ij = 0; ij < totalnum; ij++) {
                    beta_H[ij] = beta_L[ij];
                }
                L_H = L_L;
                // check new high point
                loop_check = true;
                temp_step = step_size;
                while ((loop_check) && (temp_step > 1e-3)){
                    // assign new high point
                    beta_L[para_number] = beta_H[para_number] - temp_step;
                    for (int ij = 0; ij < totalnum; ij++) {
                        beta_0[ij] = beta_L[ij];
                    }
                    temp_step = temp_step * 0.5;
                    //
                    reg_out = Pois_Full_Run(PyrC, reqrdnum, tform, totalnum, fir, R, Rd, Rdd, s_weights, beta_0, RdR, RddR, Ll, Lld, Lldd, nthreads, KeepConstant, verbose, model_bool, iter_stop, term_tot, dint, dslp, dose_abs_max, abs_max, df0, T0, Td0, Tdd0, Te, Dose, nonDose, TTerm, nonDose_LIN, nonDose_PLIN, nonDose_LOGLIN, modelform, gmix_theta, gmix_term, convgd, lr, optim_para, maxiter, double_step, Lin_Sys, Lin_Res, term_n, dfc, halfmax, epsilon, deriv_epsilon);
                    temp_L[0] = reg_out["LogLik"];
                    if (!is_nan(temp_L)[0]){
                        loop_check = false;
                    }
                }
                if (loop_check){
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
            if (!limit_hit[0]){
                // now set the mid point value
                for (int ij = 0; ij < totalnum; ij++) {
                    beta_M[ij] = (beta_H[ij] + beta_L[ij])/2;
                    beta_0[ij] = beta_M[ij];
                }
                reg_out = Pois_Full_Run(PyrC, reqrdnum, tform, totalnum, fir, R, Rd, Rdd, s_weights, beta_0, RdR, RddR, Ll, Lld, Lldd, nthreads, KeepConstant, verbose, model_bool, iter_stop, term_tot, dint, dslp, dose_abs_max, abs_max, df0, T0, Td0, Tdd0, Te, Dose, nonDose, TTerm, nonDose_LIN, nonDose_PLIN, nonDose_LOGLIN, modelform, gmix_theta, gmix_term, convgd, lr, optim_para, maxiter, double_step, Lin_Sys, Lin_Res, term_n, dfc, halfmax, epsilon, deriv_epsilon);
                reg_beta = reg_out["beta_0"];
                for (int ij = 0; ij < totalnum; ij++) {
                    beta_M[ij] = reg_beta[ij];
                }
                L_M = reg_out["LogLik"];
            }
        }
        if ((abs(beta_L[para_number] - beta_H[para_number]) < epsilon) && (! limit_hit[0])){
            limit_converged[0] = true;
        }
        limits[0] = beta_M[para_number];
        ll_final[0] = L_M;
        width_final[0] = abs(beta_L[para_number] - beta_H[para_number]);
        step_final[0] = step;
    }
    res_list = List::create(_["Parameter_Limits"] = wrap(limits), _["Negative_Risk_Limit_Hit"] = wrap(limit_hit), _["Likelihood_Boundary"] = wrap(ll_final), _["Likelihood_Goal"] = wrap(Lstar), _["Limit_Converged"] = wrap(limit_converged), _["Final_Window_Width"] = wrap(width_final), _["Final_Step"] = wrap(step_final), _["Status"] = "PASSED");
    // returns a list of results
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
// [[Rcpp::export]]
List LogLik_CaseCon_Omnibus(IntegerVector term_n, StringVector tform, NumericMatrix& a_ns, NumericMatrix& x_all, IntegerVector dfc, int fir, string modelform, double lr, List optim_para, NumericVector maxiters, int guesses, int halfmax, double epsilon, double abs_max, double dose_abs_max, double deriv_epsilon, NumericMatrix df_groups, NumericVector tu, int double_step, int verbose, IntegerVector KeepConstant, int term_tot, string ties_method, int nthreads, NumericVector& Strata_vals, List model_bool, const double gmix_theta, const IntegerVector gmix_term, const MatrixXd Lin_Sys, const VectorXd Lin_Res) {
    //
    List temp_list = List::create(_["Status"] = "TEMP");  // used as a dummy return value for code checking
    // Time durations are measured from this point on in microseconds
//    time_point<system_clock> start_point, end_point;
//    start_point = system_clock::now();
//    auto start = time_point_cast<microseconds>(start_point).time_since_epoch().count();
//    end_point = system_clock::now();
//    auto ending = time_point_cast<microseconds>(end_point).time_since_epoch().count();  // the time duration is tracked
    if (model_bool["constraint"]) {
        if (verbose >= 1) {
            Rcout << "linear constataints are currently not compatable with Case-Control model calculation" << endl;
        }
        temp_list = List::create(_["Status"] = "FAILED_WITH_BAD_MODEL_CONSTRAINT", _["LogLik"] = R_NaN);
        return temp_list;
    }
    //
    // df0: covariate data
    // ntime: number of event times for Cox PH
    // totalnum: number of terms used
    //
    // ------------------------------------------------------------------------- // initialize
    const Map<MatrixXd> df0(as<Map<MatrixXd> >(x_all));
    const int mat_row = df0.rows();
    int ntime = tu.size();
    int totalnum;
    int reqrdnum;
    // ------------------------------------------------------------------------- // initialize
    if (!model_bool["null"]) {
        totalnum = term_n.size();
        reqrdnum = totalnum - sum(KeepConstant);
    } else {
        totalnum = 1;
        reqrdnum = 1;
    }
    //
    // cout.precision: controls the number of significant digits printed
    // nthreads: number of threads used for parallel operations
    //
    Rcout.precision(7);  // forces higher precision numbers printed to terminal
    //
    // Lld_worst: stores the highest magnitude log-likelihood derivative
    //
    //
    double Lld_worst = 0.0;  // stores derivative value used to determine if every parameter is near convergence
    //
    // ---------------------------------------------
    // To Start, needs to seperate the derivative terms
    // ---------------------------------------------
    //
    // ------------------------------------------------------------------------- // initialize
    NumericVector a_n = a_ns.row(0);
    Map<VectorXd> beta_0(as<Map<VectorXd> >(a_n));
    MatrixXd T0;
    MatrixXd Td0;
    MatrixXd Tdd0;
    //
    MatrixXd Te;
    MatrixXd R;
    ColXd Rd;
    ColXd Rdd;
    //
    MatrixXd Dose;
    MatrixXd nonDose;
    MatrixXd nonDose_LIN;
    MatrixXd nonDose_PLIN;
    MatrixXd nonDose_LOGLIN;
    MatrixXd TTerm;
    double dint = 0.0;  // the amount of change used to calculate derivatives in threshold paramters
    double dslp = 0.0;
    ColXd RdR;
    ColXd RddR;
    // ------------------------------------------------------------------------- // initialize
    if (!model_bool["null"]) {
        // ---------------------------------------------
        // To Start, needs to seperate the derivative terms
        // ---------------------------------------------
        //
        Cox_Refresh_R_TERM(totalnum, reqrdnum, term_tot, dint, dslp, dose_abs_max, abs_max, df0, T0, Td0, Tdd0, Te, R, Rd, Rdd, Dose, nonDose, TTerm, nonDose_LIN, nonDose_PLIN, nonDose_LOGLIN, RdR, RddR, model_bool);
    } else {
        R = MatrixXd::Constant(mat_row, 1, 1.0);
    }
    // ------------------------------------------------------------------------- // initialize
    int group_num = ntime*Strata_vals.size();
    IntegerMatrix RiskFail(group_num, 2);
    vector<vector<int> > RiskPairs(group_num);
    const Map<MatrixXd> df_m(as<Map<MatrixXd> >(df_groups));
    //
    vector<vector<double> > Recur_Base(group_num);
    vector<vector<vector<double> > > Recur_First(group_num, vector<vector<double>>(reqrdnum));
    vector<vector<vector<double> > > Recur_Second(group_num, vector<vector<double>>(reqrdnum*(reqrdnum + 1)/2));
    vector<double> strata_odds(group_num, 0.0);
    vector<double> strata_def(group_num, 0.0);
    vector<int> strata_cond(group_num, 1);
    // ------------------------------------------------------------------------- // initialize
    if (model_bool["time_risk"]){
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
    for (int i=0; i<group_num; i++){
        strata_def[i] = strata_odds[i];
    }
    int reqrdcond = group_num - std::reduce(strata_cond.begin(), strata_cond.end());
//    return temp_list;
    // ------------------------------------------------------------------------- // initialize
    vector<double> Ll(reqrdnum, 0.0);  // log-likelihood values
    vector<double> Lld(reqrdnum, 0.0);  // log-likelihood derivative values
    vector<double> Lldd(pow(reqrdnum, 2), 0.0);  // the second derivative matrix has room for every combination, but only the lower triangle is calculated initially
    //
    vector<double> LldOdds(reqrdcond, 0.0); // storing the first derivative with respect to each strata odds ratio
    vector<double> LlddOdds(reqrdcond, 0.0); // storing the second derivative with respect to each strata odds ratio
    vector<double> LlddOddsBeta(reqrdcond*reqrdnum, 0.0); // storing the second derivative with respect to each strata odds ratio and model parameter
    double dev = 0.0; // deviance needs to be calculated seperately, because the conditional and unconditional portions have different formula?
    if (model_bool["null"]){
        Calculate_Recursive(model_bool, group_num, RiskFail, RiskPairs, totalnum, ntime, R, Rd, Rdd, Recur_Base, Recur_First, Recur_Second, nthreads, KeepConstant);
        Calc_Recur_LogLik(model_bool, group_num, RiskFail, RiskPairs, totalnum, ntime, R, Rd, Rdd, RdR, RddR, dev, Ll, Lld, Lldd, Recur_Base, Recur_First, Recur_Second, strata_odds, nthreads, KeepConstant, strata_cond, LldOdds, LlddOdds, LlddOddsBeta);
        //
        List res_list = List::create(_["LogLik"] = wrap(Ll[0]), _["Deviance"] = wrap(dev), _["Status"] = "PASSED", _["StrataOdds"]=wrap(strata_odds), _["FreeParameters"]=wrap(reqrdnum), _["FreeSets"]=wrap(reqrdcond));
        // returns a list of results
        return res_list;
    }
    // ------------------------------------------------------------------------- // initialize
    // the log-likelihood is calculated in parallel over the risk groups
    vector <double> Ll_comp(2, Ll[0]);  // vector to compare values
    double abs_max0 = abs_max;
    double dose_abs_max0 = dose_abs_max;
    //
    vector<double> dbeta(totalnum, 0.0);
    vector<double> dstrata(group_num, 0.0);
    NumericVector m_g_store(reqrdnum+reqrdcond);
    NumericVector v_beta_store(reqrdnum+reqrdcond);
    //
    // --------------------------
    // always starts from initial guess
    // --------------------------
    vector<double> beta_c(totalnum, 0.0);
    vector<double> beta_a(totalnum, 0.0);
    vector<double> beta_best(totalnum, 0.0);
    vector<double> beta_p(totalnum, 0.0);
    VectorXd::Map(&beta_p[0], beta_0.size()) = beta_0;  // stores previous parameters
    VectorXd::Map(&beta_c[0], beta_0.size()) = beta_0;  // stores current parameters
    VectorXd::Map(&beta_a[0], beta_0.size()) = beta_0;  // stores a refrence value for parameters
    VectorXd::Map(&beta_best[0], beta_0.size()) = beta_0;  // stores the best parameters
    //
    vector<double> strata_c(group_num, 0.0);
    vector<double> strata_a(group_num, 0.0);
    vector<double> strata_best(group_num, 0.0);
    vector<double> strata_p(group_num, 0.0);
    for (int i=0; i<group_num; i++){
        strata_p[i] = strata_odds[i];
        strata_c[i] = strata_odds[i];
        strata_a[i] = strata_odds[i];
        strata_best[i] = strata_odds[i];
    }
    //
    double halves = 0;  // number of half-steps taken
    int ind0 = fir;  // used for validations
    int iteration = 0;  // iteration number
    int maxiter = 0;
    //
    bool convgd = FALSE;
    int iter_stop  = 0;  // tracks if the iterations should be stopped for convergence
    int iter_check = 0;  // signal to check for convergence
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
    // Variables that are used for the risk check function shared across cox, poisson, and log bound functions
    MatrixXd dev_temp = MatrixXd::Zero(1, 1);
    double Lstar = 0.0;
    MatrixXd PyrC = MatrixXd::Zero(1, 1);
    //
    //
    for (int guess = 0; guess <guesses; guess++) {
        Cox_Refresh_R_TERM(totalnum, reqrdnum, term_tot, dint, dslp, dose_abs_max, abs_max, df0, T0, Td0, Tdd0, Te, R, Rd, Rdd, Dose, nonDose, TTerm, nonDose_LIN, nonDose_PLIN, nonDose_LOGLIN, RdR, RddR, model_bool);
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
        abs_max = abs_max0;
        dose_abs_max = dose_abs_max0;
        iter_stop = 0;
        halves = 0;
        iteration = 0;
        halves = 0;  // number of half-steps taken
        ind0 = fir;  // used for validations
        iteration = 0;  // iteration number
        Ll_iter_best = 10;
        //
        convgd = FALSE;
        iter_stop  = 0;  // tracks if the iterations should be stopped for convergence
        iter_check = 0;  // signal to check for convergence
        //
        maxiter = maxiters[guess];
        a_n = a_ns.row(guess);
        for (int i = 0; i < beta_0.size(); i++) {
            beta_0[i] = a_n[i];
        }
        for (int i=0; i<group_num; i++){
            strata_odds[i] = strata_def[i];
        }
        if (verbose >= 4) {
            Rcout << "C++ Note: starting guess " << guess << endl;
        }
        Cox_Term_Risk_Calc(modelform, tform, term_n, totalnum, fir, dfc, term_tot, T0, Td0, Tdd0, Te, R, Rd, Rdd, Dose, nonDose, beta_0, df0, dose_abs_max, abs_max, TTerm, nonDose_LIN, nonDose_PLIN, nonDose_LOGLIN, RdR, RddR, nthreads, KeepConstant, verbose, model_bool, gmix_theta, gmix_term);
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
            temp_list = List::create(_["beta_0"] = wrap(beta_0), _["StrataOdds"]=wrap(strata_odds), _["Deviation"] = R_NaN, _["Status"] = "FAILED_WITH_NEGATIVE_RISK", _["LogLik"] = R_NaN);
            return temp_list;
        }
        Calculate_Recursive(model_bool, group_num, RiskFail, RiskPairs, totalnum, ntime, R, Rd, Rdd, Recur_Base, Recur_First, Recur_Second , nthreads, KeepConstant);
        Calc_Recur_LogLik(model_bool, group_num, RiskFail, RiskPairs, totalnum, ntime, R, Rd, Rdd, RdR, RddR, dev, Ll, Lld, Lldd, Recur_Base, Recur_First, Recur_Second, strata_odds, nthreads, KeepConstant, strata_cond, LldOdds, LlddOdds, LlddOddsBeta);
        Print_LL(reqrdnum, totalnum, beta_0, Ll, Lld, Lldd, verbose, model_bool);
        Print_LL_Background(reqrdnum, totalnum, group_num, reqrdcond, strata_odds, LldOdds, LlddOdds, LlddOddsBeta, verbose, model_bool);
        // NOW WE RUN THE ITERATIONS
        for (int i = 0; i < beta_0.size(); i++) {
            beta_c[i] = beta_0[i];
        }
        for (int i = 0; i < strata_odds.size(); i++) {
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
                Calc_Change_Background_Gradient(nthreads, model_bool, totalnum, group_num, optim_para, iteration, abs_max, Lld, m_g_store, v_beta_store, dbeta, KeepConstant, strata_cond, LldOdds, dstrata);
            } else {
                Calc_Change_Background(double_step, nthreads, totalnum, group_num, dose_abs_max, lr, abs_max, Ll, Lld, Lldd, dbeta, tform, dose_abs_max, abs_max, KeepConstant, strata_cond, LldOdds, LlddOdds, LlddOddsBeta, dstrata);
                Intercept_Bound(nthreads, totalnum, beta_0, dbeta, dfc, df0, KeepConstant, tform);
            }
            if ((Ll_iter_best > 0) || (Ll_iter_best < Ll[ind0])) {
                Ll_iter_best = Ll[ind0];
            }
            //
            if (model_bool["gradient"]){
                //
                for (int ij = 0; ij < totalnum; ij++) {
                    beta_0[ij] = beta_a[ij] + dbeta[ij];
                    beta_c[ij] = beta_0[ij];
                }
                for (int ij = 0; ij < group_num; ij++) {
                    strata_odds[ij] = strata_a[ij] + dstrata[ij];
                    strata_c[ij] = strata_odds[ij];
                }
                Cox_Term_Risk_Calc(modelform, tform, term_n, totalnum, fir, dfc, term_tot, T0, Td0, Tdd0, Te, R, Rd, Rdd, Dose, nonDose, beta_0, df0, dose_abs_max, abs_max, TTerm, nonDose_LIN, nonDose_PLIN, nonDose_LOGLIN, RdR, RddR, nthreads, KeepConstant, verbose, model_bool, gmix_theta, gmix_term);
                //
                if ((R.minCoeff() <= 0) || (R.hasNaN())) {
                    #ifdef _OPENMP
                    #pragma omp parallel for schedule(dynamic) num_threads(nthreads)
                    #endif
                    for (int ijk = 0; ijk < totalnum; ijk++) {
                        int tij = term_n[ijk];
                        if (TTerm.col(tij).minCoeff()<=0) {
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
                    if (Ll[ind0] <= Ll_abs_best) {  // if a better point wasn't found, takes a half-step
                        #ifdef _OPENMP
                        #pragma omp parallel for schedule(dynamic) num_threads(nthreads)
                        #endif
                        for (int ijk = 0; ijk < totalnum; ijk++) {
                            dbeta[ijk] = dbeta[ijk] * 0.5;  //
                        }
                        #ifdef _OPENMP
                        #pragma omp parallel for schedule(dynamic) num_threads(nthreads)
                        #endif
                        for (int ijk = 0; ijk < group_num; ijk++) {
                            dstrata[ijk] = dstrata[ijk] * 0.5;  //
                        }
                    } else{  // if improved, updates the best vector
                        #ifdef _OPENMP
                        #pragma omp parallel for schedule(dynamic) num_threads(nthreads)
                        #endif
                        for (int ijk = 0; ijk < totalnum; ijk++) {
                            beta_best[ijk] = beta_c[ijk];
                        }
                        #ifdef _OPENMP
                        #pragma omp parallel for schedule(dynamic) num_threads(nthreads)
                        #endif
                        for (int ijk = 0; ijk < group_num; ijk++) {
                            strata_best[ijk] = strata_odds[ijk];  //
                        }
                    }
                    #ifdef _OPENMP
                    #pragma omp parallel for schedule(dynamic) num_threads(nthreads)
                    #endif
                    for (int ijk = 0; ijk < totalnum; ijk++) {  // totalnum*(totalnum + 1)/2
                        beta_0[ijk] = beta_c[ijk];
                    }
                    #ifdef _OPENMP
                    #pragma omp parallel for schedule(dynamic) num_threads(nthreads)
                    #endif
                    for (int ijk = 0; ijk < group_num; ijk++) {
                        strata_odds[ijk] = strata_c[ijk];  //
                    }
                }
            } else {
                halves = 0;
                while ((Ll[ind0] <= Ll_iter_best) && (halves < halfmax)) {  // repeats until half-steps maxed or an improvement
                    for (int ij = 0; ij < totalnum; ij++) {
                        beta_0[ij] = beta_a[ij] + dbeta[ij];
                        beta_c[ij] = beta_0[ij];
                    }
                    for (int ij = 0; ij < group_num; ij++) {
                        strata_odds[ij] = strata_a[ij] + dstrata[ij];
                        strata_c[ij] = strata_odds[ij];
                    }
                    // ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------//
                    // The same subterm, risk, sides, and log-likelihood calculations are performed every half-step and iteration
                    // ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------//
                    Cox_Term_Risk_Calc(modelform, tform, term_n, totalnum, fir, dfc, term_tot, T0, Td0, Tdd0, Te, R, Rd, Rdd, Dose, nonDose, beta_0, df0, dose_abs_max, abs_max, TTerm, nonDose_LIN, nonDose_PLIN, nonDose_LOGLIN, RdR, RddR, nthreads, KeepConstant, verbose, model_bool, gmix_theta, gmix_term);
                    //
                    if ((R.minCoeff() <= 0) || (R.hasNaN())) {
                        #ifdef _OPENMP
                        #pragma omp parallel for schedule(dynamic) num_threads(nthreads)
                        #endif
                        for (int ijk = 0; ijk < totalnum; ijk++) {
                            int tij = term_n[ijk];
                            if (TTerm.col(tij).minCoeff()<=0) {
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
                        if (Ll[ind0] <= Ll_abs_best) {  // if a better point wasn't found, takes a half-step
                            #ifdef _OPENMP
                            #pragma omp parallel for schedule(dynamic) num_threads(nthreads)
                            #endif
                            for (int ijk = 0; ijk < totalnum; ijk++) {
                                dbeta[ijk] = dbeta[ijk] * 0.5;  //
                            }
                            #ifdef _OPENMP
                            #pragma omp parallel for schedule(dynamic) num_threads(nthreads)
                            #endif
                            for (int ijk = 0; ijk < group_num; ijk++) {
                                dstrata[ijk] = dstrata[ijk] * 0.5;  //
                            }
                        } else{  // if improved, updates the best vector
                            #ifdef _OPENMP
                            #pragma omp parallel for schedule(dynamic) num_threads(nthreads)
                            #endif
                            for (int ijk = 0; ijk < totalnum; ijk++) {
                                beta_best[ijk] = beta_c[ijk];
                            }
                            #ifdef _OPENMP
                            #pragma omp parallel for schedule(dynamic) num_threads(nthreads)
                            #endif
                            for (int ijk = 0; ijk < group_num; ijk++) {
                                strata_best[ijk] = strata_odds[ijk];  //
                            }
                        }
                        #ifdef _OPENMP
                        #pragma omp parallel for schedule(dynamic) num_threads(nthreads)
                        #endif
                        for (int ijk = 0; ijk < totalnum; ijk++) {  // totalnum*(totalnum + 1)/2
                            beta_0[ijk] = beta_c[ijk];
                        }
                        #ifdef _OPENMP
                        #pragma omp parallel for schedule(dynamic) num_threads(nthreads)
                        #endif
                        for (int ijk = 0; ijk < group_num; ijk++) {
                            strata_odds[ijk] = strata_c[ijk];  //
                        }
                    }
                }
                if (beta_best != beta_c) {  // if the risk matrices aren't the optimal values, then they must be recalculated
                    // If it goes through every half step without improvement, then the maximum change needs to be decreased
                    abs_max = abs_max*pow(0.5, halfmax);  // reduces the step sizes
                    dose_abs_max = dose_abs_max*pow(0.5, halfmax);
                    iter_check = 1;
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
                    Cox_Term_Risk_Calc(modelform, tform, term_n, totalnum, fir, dfc, term_tot, T0, Td0, Tdd0, Te, R, Rd, Rdd, Dose, nonDose, beta_0, df0, dose_abs_max, abs_max, TTerm, nonDose_LIN, nonDose_PLIN, nonDose_LOGLIN, RdR, RddR, nthreads, KeepConstant, verbose, model_bool, gmix_theta, gmix_term);
                }
            }
            for (int ij = 0; ij < reqrdnum; ij++) {
                if (abs(Lld[ij]) > Lld_worst) {
                    Lld_worst = abs(Lld[ij]);
                }
            }
            if ((iteration % (reqrdnum)) || (iter_check == 1)) {  // checks every set number of iterations
                iter_check = 0;
                if (Lld_worst < deriv_epsilon) {  // ends if the derivatives are low enough
                    iter_stop = 1;
                }
                if (abs_max < epsilon/10) {  // if the maximum change is too low, then it ends
                    iter_stop = 1;
                }
            }
        }
        // -----------------------------------------------
        // Performing Full Calculation to get full second derivative matrix
        // -----------------------------------------------
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
        for (int i=0; i<group_num; i++){
            strata_fin(guess,i) = strata_odds[i];
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
    // ----------------------------------------------------------------------------------- //
    //             NOW WE NEED TO FIND THE BEST, ASSIGN, AND CONTINUE REGRESSION
    // ----------------------------------------------------------------------------------- //
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
    abs_max = abs_max0;
    dose_abs_max = dose_abs_max0;
    iter_stop = 0;
    halves = 0;
    iteration = 0;
    halves = 0;  // number of half-steps taken
    ind0 = fir;  // used for validations
    iteration = 0;  // iteration number
    //
    convgd = FALSE;
    iter_stop  = 0;  // tracks if the iterations should be stopped for convergence
    iter_check = 0;  // signal to check for convergence
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
    // Calculates the subterm and term values
    Cox_Term_Risk_Calc(modelform, tform, term_n, totalnum, fir, dfc, term_tot, T0, Td0, Tdd0, Te, R, Rd, Rdd, Dose, nonDose, beta_0, df0, dose_abs_max, abs_max, TTerm, nonDose_LIN, nonDose_PLIN, nonDose_LOGLIN, RdR, RddR, nthreads, KeepConstant, verbose, model_bool, gmix_theta, gmix_term);
    //
    // -------------------------------------------------------------------------------------------
    Calculate_Recursive(model_bool, group_num, RiskFail, RiskPairs, totalnum, ntime, R, Rd, Rdd, Recur_Base, Recur_First, Recur_Second , nthreads, KeepConstant);
    Calc_Recur_LogLik(model_bool, group_num, RiskFail, RiskPairs, totalnum, ntime, R, Rd, Rdd, RdR, RddR, dev, Ll, Lld, Lldd, Recur_Base, Recur_First, Recur_Second, strata_odds, nthreads, KeepConstant, strata_cond, LldOdds, LlddOdds, LlddOddsBeta);
    Print_LL(reqrdnum, totalnum, beta_0, Ll, Lld, Lldd, verbose, model_bool);
    Print_LL_Background(reqrdnum, totalnum, group_num, reqrdcond, strata_odds, LldOdds, LlddOdds, LlddOddsBeta, verbose, model_bool);
    //
    while ((iteration < maxiter) && (iter_stop == 0)) {
//        Rcout << " " << endl;
//        Rcout << "Strata val:";
//        for (int ij = 0; ij < group_num; ij++) {
//             Rcout << " " << strata_odds[ij];
//        }
//        Rcout << " " << endl;
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
            Calc_Change_Background_Gradient(nthreads, model_bool, totalnum, group_num, optim_para, iteration, abs_max, Lld, m_g_store, v_beta_store, dbeta, KeepConstant, strata_cond, LldOdds, dstrata);
        } else {
            Calc_Change_Background(double_step, nthreads, totalnum, group_num, dose_abs_max, lr, abs_max, Ll, Lld, Lldd, dbeta, tform, dose_abs_max, abs_max, KeepConstant, strata_cond, LldOdds, LlddOdds, LlddOddsBeta, dstrata);
            Intercept_Bound(nthreads, totalnum, beta_0, dbeta, dfc, df0, KeepConstant, tform);
        }
//        Rcout << "Beta change:";
//        for (int ij = 0; ij < totalnum; ij++) {
//            Rcout << " " << dbeta[ij];
//        }
//        Rcout << " " << endl;
//        Rcout << "Strata change:";
//        for (int ij = 0; ij < group_num; ij++) {
//             Rcout << " " << dstrata[ij];
//        }
//        Rcout << " " << endl;
        if ((Ll_iter_best > 0) || (Ll_iter_best < Ll[ind0])) {
            Ll_iter_best = Ll[ind0];
        }
        //
        if (model_bool["gradient"]){
            //
            for (int ij = 0; ij < totalnum; ij++) {
                beta_0[ij] = beta_a[ij] + dbeta[ij];
                beta_c[ij] = beta_0[ij];
            }
            for (int ij = 0; ij < group_num; ij++) {
                strata_odds[ij] = strata_a[ij] + dstrata[ij];
                strata_c[ij] = strata_odds[ij];
            }
            Cox_Term_Risk_Calc(modelform, tform, term_n, totalnum, fir, dfc, term_tot, T0, Td0, Tdd0, Te, R, Rd, Rdd, Dose, nonDose, beta_0, df0, dose_abs_max, abs_max, TTerm, nonDose_LIN, nonDose_PLIN, nonDose_LOGLIN, RdR, RddR, nthreads, KeepConstant, verbose, model_bool, gmix_theta, gmix_term);
            //
            if ((R.minCoeff() <= 0) || (R.hasNaN())) {
                #ifdef _OPENMP
                #pragma omp parallel for schedule(dynamic) num_threads(nthreads)
                #endif
                for (int ijk = 0; ijk < totalnum; ijk++) {
                    int tij = term_n[ijk];
                    if (TTerm.col(tij).minCoeff()<=0) {
                        dbeta[ijk] = dbeta[ijk] / 2.0;
                    } else if (isinf(TTerm.col(tij).maxCoeff())) {
                        dbeta[ijk] = dbeta[ijk] / 2.0;
                    } else if (isnan(TTerm.col(tij).minCoeff())) {
                        dbeta[ijk] = dbeta[ijk] / 2.0;
                    }
                }
            } else {
                Calculate_Recursive(model_bool, group_num, RiskFail, RiskPairs, totalnum, ntime, R, Rd, Rdd, Recur_Base, Recur_First, Recur_Second , nthreads, KeepConstant);
                Calc_Recur_LogLik(model_bool, group_num, RiskFail, RiskPairs, totalnum, ntime, R, Rd, Rdd, RdR, RddR, dev, Ll, Lld, Lldd, Recur_Base, Recur_First, Recur_Second, strata_odds, nthreads, KeepConstant, strata_cond, LldOdds, LlddOdds, LlddOddsBeta);
                Print_LL(reqrdnum, totalnum, beta_0, Ll, Lld, Lldd, verbose, model_bool);
                Print_LL_Background(reqrdnum, totalnum, group_num, reqrdcond, strata_odds, LldOdds, LlddOdds, LlddOddsBeta, verbose, model_bool);
                //
                if (Ll[ind0] <= Ll_iter_best) {  // if a better point wasn't found, takes a half-step
                    #ifdef _OPENMP
                    #pragma omp parallel for schedule(dynamic) num_threads(nthreads)
                    #endif
                    for (int ijk = 0; ijk < totalnum; ijk++) {
                        dbeta[ijk] = dbeta[ijk] * 0.5;  //
                    }
                    #ifdef _OPENMP
                    #pragma omp parallel for schedule(dynamic) num_threads(nthreads)
                    #endif
                    for (int ijk = 0; ijk < group_num; ijk++) {
                        dstrata[ijk] = dstrata[ijk] * 0.5;  //
                    }
                } else{  // if improved, updates the best vector
                    #ifdef _OPENMP
                    #pragma omp parallel for schedule(dynamic) num_threads(nthreads)
                    #endif
                    for (int ijk = 0; ijk < totalnum; ijk++) {
                        beta_best[ijk] = beta_c[ijk];
                    }
                    #ifdef _OPENMP
                    #pragma omp parallel for schedule(dynamic) num_threads(nthreads)
                    #endif
                    for (int ijk = 0; ijk < group_num; ijk++) {
                        strata_best[ijk] = strata_odds[ijk];  //
                    }
                }
                #ifdef _OPENMP
                #pragma omp parallel for schedule(dynamic) num_threads(nthreads)
                #endif
                for (int ijk = 0; ijk < totalnum; ijk++) {  // totalnum*(totalnum + 1)/2
                    beta_0[ijk] = beta_c[ijk];
                }
                #ifdef _OPENMP
                #pragma omp parallel for schedule(dynamic) num_threads(nthreads)
                #endif
                for (int ijk = 0; ijk < group_num; ijk++) {
                    strata_odds[ijk] = strata_c[ijk];  //
                }
            }
        } else {
            halves = 0;
            while ((Ll[ind0] <= Ll_iter_best) && (halves < halfmax)) {  // repeats until half-steps maxed or an improvement
                for (int ij = 0; ij < totalnum; ij++) {
                    beta_0[ij] = beta_a[ij] + dbeta[ij];
                    beta_c[ij] = beta_0[ij];
                }
                for (int ij = 0; ij < group_num; ij++) {
                    strata_odds[ij] = strata_a[ij] + dstrata[ij];
                    strata_c[ij] = strata_odds[ij];
                }
//                Rcout << "Beta val:";
//                for (int ij = 0; ij < totalnum; ij++) {
//                    Rcout << " " << beta_0[ij];
//                }
//                Rcout << " " << endl;
//                Rcout << "Strata val:";
//                for (int ij = 0; ij < group_num; ij++) {
//                     Rcout << " " << strata_odds[ij];
//                }
//                Rcout << " " << endl;
                // ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------//
                // The same subterm, risk, sides, and log-likelihood calculations are performed every half-step and iteration
                // ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------//
                Cox_Term_Risk_Calc(modelform, tform, term_n, totalnum, fir, dfc, term_tot, T0, Td0, Tdd0, Te, R, Rd, Rdd, Dose, nonDose, beta_0, df0, dose_abs_max, abs_max, TTerm, nonDose_LIN, nonDose_PLIN, nonDose_LOGLIN, RdR, RddR, nthreads, KeepConstant, verbose, model_bool, gmix_theta, gmix_term);
                //
                if ((R.minCoeff() <= 0) || (R.hasNaN())) {
                    #ifdef _OPENMP
                    #pragma omp parallel for schedule(dynamic) num_threads(nthreads)
                    #endif
                    for (int ijk = 0; ijk < totalnum; ijk++) {
                        int tij = term_n[ijk];
                        if (TTerm.col(tij).minCoeff()<=0) {
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
                    if (Ll[ind0] <= Ll_iter_best) {  // if a better point wasn't found, takes a half-step
                        #ifdef _OPENMP
                        #pragma omp parallel for schedule(dynamic) num_threads(nthreads)
                        #endif
                        for (int ijk = 0; ijk < totalnum; ijk++) {
                            dbeta[ijk] = dbeta[ijk] * 0.5;  //
                        }
                        #ifdef _OPENMP
                        #pragma omp parallel for schedule(dynamic) num_threads(nthreads)
                        #endif
                        for (int ijk = 0; ijk < group_num; ijk++) {
                            dstrata[ijk] = dstrata[ijk] * 0.5;  //
                        }
                    } else{  // if improved, updates the best vector
                        #ifdef _OPENMP
                        #pragma omp parallel for schedule(dynamic) num_threads(nthreads)
                        #endif
                        for (int ijk = 0; ijk < totalnum; ijk++) {
                            beta_best[ijk] = beta_c[ijk];
                        }
                        #ifdef _OPENMP
                        #pragma omp parallel for schedule(dynamic) num_threads(nthreads)
                        #endif
                        for (int ijk = 0; ijk < group_num; ijk++) {
                            strata_best[ijk] = strata_odds[ijk];  //
                        }
                    }
                    #ifdef _OPENMP
                    #pragma omp parallel for schedule(dynamic) num_threads(nthreads)
                    #endif
                    for (int ijk = 0; ijk < totalnum; ijk++) {  // totalnum*(totalnum + 1)/2
                        beta_0[ijk] = beta_c[ijk];
                    }
                    #ifdef _OPENMP
                    #pragma omp parallel for schedule(dynamic) num_threads(nthreads)
                    #endif
                    for (int ijk = 0; ijk < group_num; ijk++) {
                        strata_odds[ijk] = strata_c[ijk];  //
                    }
                }
            }
            if (beta_best != beta_c) {  // if the risk matrices aren't the optimal values, then they must be recalculated
                // If it goes through every half step without improvement, then the maximum change needs to be decreased
                abs_max = abs_max*pow(0.5, halfmax);  // reduces the step sizes
                dose_abs_max = dose_abs_max*pow(0.5, halfmax);
                iter_check = 1;
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
                Cox_Term_Risk_Calc(modelform, tform, term_n, totalnum, fir, dfc, term_tot, T0, Td0, Tdd0, Te, R, Rd, Rdd, Dose, nonDose, beta_0, df0, dose_abs_max, abs_max, TTerm, nonDose_LIN, nonDose_PLIN, nonDose_LOGLIN, RdR, RddR, nthreads, KeepConstant, verbose, model_bool, gmix_theta, gmix_term);
            }
        }
        for (int ij = 0; ij < reqrdnum; ij++) {
            if (abs(Lld[ij]) > Lld_worst) {
                Lld_worst = abs(Lld[ij]);
            }
        }
        if ((iteration % (reqrdnum)) || (iter_check == 1)) {  // checks every set number of iterations
            iter_check = 0;
            if (Lld_worst < deriv_epsilon) {  // ends if the derivatives are low enough
                iter_stop = 1;
            }
            if (abs_max < epsilon/10) {  // if the maximum change is too low, then it ends
                iter_stop = 1;
            }
        }
    }
    // ----------------------------------------------------------------------------------- //
    //               NOW WE WRAP UP
    // ----------------------------------------------------------------------------------- //
    for (int ij = 0; ij < totalnum; ij++) {
        beta_0[ij] = beta_best[ij];
    }
    for (int ij = 0; ij < group_num; ij++) {
        strata_odds[ij] = strata_best[ij];
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
    Calculate_Recursive(model_bool, group_num, RiskFail, RiskPairs, totalnum, ntime, R, Rd, Rdd, Recur_Base, Recur_First, Recur_Second , nthreads, KeepConstant);
    Calc_Recur_LogLik(model_bool, group_num, RiskFail, RiskPairs, totalnum, ntime, R, Rd, Rdd, RdR, RddR, dev, Ll, Lld, Lldd, Recur_Base, Recur_First, Recur_Second, strata_odds, nthreads, KeepConstant, strata_cond, LldOdds, LlddOdds, LlddOddsBeta);
    Print_LL(reqrdnum, totalnum, beta_0, Ll, Lld, Lldd, verbose, model_bool);
    Print_LL_Background(reqrdnum, totalnum, group_num, reqrdcond, strata_odds, LldOdds, LlddOdds, LlddOddsBeta, verbose, model_bool);
    //
    //
    List res_list;
    //
    if (model_bool["single"]) {
        res_list = List::create(_["LogLik"] = wrap(Ll[0]), _["Deviance"] = wrap(dev), _["beta_0"] = wrap(beta_0), _["StrataOdds"]=wrap(strata_odds), _["FreeParameters"]=wrap(reqrdnum), _["FreeSets"]=wrap(reqrdcond), _["Status"] = "PASSED");
        // returns a list of results
        return res_list;
    }
    List para_list;
    if (!model_bool["basic"]) {
        para_list = List::create(_["term_n"] = term_n, _["tforms"] = tform);  // stores the term information
    }
    List control_list = List::create(_["Iteration"] = iteration, _["Maximum Step"]= abs_max, _["Derivative Limiting"] = Lld_worst);  // stores the total number of iterations used
    //
    if (model_bool["gradient"]) {
        res_list = List::create(_["LogLik"] = wrap(Ll[0]), _["Deviance"] = wrap(dev), _["First_Der"] = wrap(Lld), _["beta_0"] = wrap(beta_0), _["StrataOdds"]=wrap(strata_odds), _["Parameter_Lists"] = para_list, _["Control_List"] = control_list, _["Converged"] = convgd, _["FreeParameters"]=wrap(reqrdnum), _["FreeSets"]=wrap(reqrdcond), _["Status"] = "PASSED");
        return res_list;
    }
    //
    NumericVector Lldd_vec(reqrdnum * reqrdnum);  // simplfied information matrix
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
    if (model_bool["oberved_info"]){
        cov = - 1 * Lldd_mat.inverse().matrix();  // uses inverse information matrix to calculate the standard deviation
        for (int ij = 0; ij < totalnum; ij++) {
            if (KeepConstant[ij] == 0) {
                int pij_ind = ij - sum(head(KeepConstant, ij));
                stdev(ij) = sqrt(cov(pij_ind, pij_ind));
            }
        }
    } else { }
    //
    res_list = List::create(_["LogLik"] = wrap(Ll[0]), _["Deviance"] = wrap(dev), _["First_Der"] = wrap(Lld), _["Second_Der"] = Lldd_vec, _["beta_0"] = wrap(beta_0), _["StrataOdds"]=wrap(strata_odds), _["Standard_Deviation"] = wrap(stdev), _["Covariance"] = wrap(cov), _["Parameter_Lists"] = para_list, _["Control_List"] = control_list, _["Converged"] = convgd, _["FreeParameters"]=wrap(reqrdnum), _["FreeSets"]=wrap(reqrdcond), _["Status"] = "PASSED");
    // returns a list of results
    return res_list;
}
