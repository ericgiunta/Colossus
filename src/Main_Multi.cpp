#include <RcppEigen.h>
#ifdef _OPENMP
#include <omp.h>
#endif
#include "Main_Multi.h"
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
