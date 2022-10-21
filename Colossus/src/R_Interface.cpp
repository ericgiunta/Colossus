#include <RcppEigen.h>
#include <RcppParallel.h>
#include <omp.h>
#include "R_Interface.h"
#include "Main_Functions.h"
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
// [[Rcpp::depends(RcppParallel)]]
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
void visit_lambda(const Mat& m, const Func& f)
{
    lambda_as_visitor_wrapper<Func> visitor(f);
    m.visit(visitor);
}


//' Interface between R code and the Cox PH regression
//' \code{cox_ph_transition} Called directly from R, Defines the control variables and calls the regression function
//' @param Term_n Term numbers
//' @param tform subterm types
//' @param a_n starting values
//' @param dfc covariate column numbers
//' @param x_all covariate matrix
//' @param fir first term number
//' @param der_iden subterm number for derivative tests
//' @param modelform model string
//' @param Control control list
//' @param df_groups time and event matrix
//' @param tu event times
//' @param KeepConstant vector of parameters to keep constant
//' @param term_tot total number of terms
//'
//' @return LogLik_Cox_PH output : Log-likelihood of optimum, first derivative of log-likelihood, second derivative matrix, parameter list, standard deviation estimate, AIC, model information
// [[Rcpp::export]]
List cox_ph_transition(IntegerVector Term_n, StringVector tform, NumericVector a_n,IntegerVector dfc,NumericMatrix x_all, int fir, int der_iden,string modelform, List Control, NumericMatrix df_groups, NumericVector tu, IntegerVector KeepConstant, int term_tot){
    bool change_all = Control["change_all"];
    bool verbose = Control["verbose"];
    bool debugging = FALSE;
    double lr = Control["lr"];
    int maxiter = Control["maxiter"];
    int halfmax = Control["halfmax"];
    double epsilon = Control["epsilon"];
    double dbeta_cap = Control["dbeta_max"];
    double abs_max = Control["abs_max"];
    double dose_abs_max = Control["dose_abs_max"];
    double deriv_epsilon =Control["deriv_epsilon"];
    string ties_method =Control["ties"];
    //
    // Performs regression
    //----------------------------------------------------------------------------------------------------------------//
    List res = LogLik_Cox_PH(Term_n, tform, a_n, x_all, dfc,fir, der_iden,modelform, lr, maxiter, halfmax, epsilon, dbeta_cap, abs_max,dose_abs_max, deriv_epsilon, df_groups, tu, change_all,verbose, debugging, KeepConstant, term_tot, ties_method);
    //----------------------------------------------------------------------------------------------------------------//
    return res;
}

//' Interface between R code and the Cox PH plotting
//' \code{cox_ph_plot} Called directly from R, Defines the control variables and calls the correct plotting function
//' @param Term_n Term numbers
//' @param tform subterm types
//' @param a_n starting values
//' @param dfc covariate column numbers
//' @param x_all covariate matrix
//' @param fir first term number
//' @param der_iden subterm number for derivative tests
//' @param modelform model string
//' @param Control control list
//' @param df_groups time and event matrix
//' @param tu event times
//' @param KeepConstant vector of parameters to keep constant
//' @param term_tot total number of terms
//' @param Plot_Type string specifying which plot type
//' @param uniq_v total number of unique covariate values
//'
//' @return Cox_PH_PLOT_SURV : ( baseline harzard, risk for each row) or Cox_PH_PLOT_RISK output : (covariate values, risks for each row)
// [[Rcpp::export]]
List cox_ph_plot(IntegerVector Term_n, StringVector tform, NumericVector a_n,IntegerVector dfc,NumericMatrix x_all, int fir, int der_iden,string modelform, List Control, NumericMatrix df_groups, NumericVector tu, IntegerVector KeepConstant, int term_tot, vector<string> Plot_Type ,int uniq_v){
    bool verbose = Control["verbose"];
    bool debugging = FALSE;
    double abs_max = Control["abs_max"];
    double dose_abs_max = Control["dose_abs_max"];
    string ties_method =Control["ties"];
    List res;
    // there are two types of plots that can be generated, survival curve and risk by covariate value
    //----------------------------------------------------------------------------------------------------------------//
    if (Plot_Type[0]=="SURV"){
        res = Cox_PH_PLOT_SURV(Term_n, tform, a_n, x_all, dfc,fir, der_iden,modelform, abs_max,dose_abs_max, df_groups, tu,verbose, debugging, KeepConstant, term_tot);
    }else if (Plot_Type[0]=="RISK"){
        res = Cox_PH_PLOT_RISK(Term_n, tform, a_n, x_all, dfc,fir, der_iden,modelform, abs_max,dose_abs_max, df_groups, tu,verbose, debugging, KeepConstant, term_tot, uniq_v);
    } else {
        throw invalid_argument("Invalid plot type");
    }
    return res;
}

//' Interface between R code and the schoenfeld residual calculation
//' \code{cox_ph_schoenfeld_transition} Called directly from R, Defines the control variables and calls the calculation function
//' @param Term_n Term numbers
//' @param tform subterm types
//' @param a_n starting values
//' @param dfc covariate column numbers
//' @param x_all covariate matrix
//' @param fir first term number
//' @param der_iden subterm number for derivative tests
//' @param modelform model string
//' @param Control control list
//' @param df_groups time and event matrix
//' @param tu event times
//' @param KeepConstant vector of parameters to keep constant
//' @param term_tot total number of terms
//'
//' @return Schoenfeld_Cox_PH output: scaled schoenfeld residuals
// [[Rcpp::export]]
NumericMatrix cox_ph_schoenfeld_transition(IntegerVector Term_n, StringVector tform, NumericVector a_n,IntegerVector dfc,NumericMatrix x_all, int fir, int der_iden,string modelform, List Control, NumericMatrix df_groups, NumericVector tu, IntegerVector KeepConstant, int term_tot){
    //
    bool verbose = Control["verbose"];
    bool debugging = FALSE;
    double abs_max = Control["abs_max"];
    double dose_abs_max = Control["dose_abs_max"];
    string ties_method =Control["ties"];
    //
    // performs schoenfeld residual calculation
    NumericMatrix res = Schoenfeld_Cox_PH(Term_n, tform, a_n, x_all, dfc,fir, der_iden,modelform, abs_max,dose_abs_max, df_groups, tu,verbose, debugging, KeepConstant, term_tot, ties_method);
    //----------------------------------------------------------------------------------------------------------------//
    return res;
}

//' Interface between R code and the poisson regression
//' \code{poisson_transition} Called directly from R, Defines the control variables and calls the regression function
//' @param dfe Matrix with person-year/event count information
//' @param Term_n Term numbers
//' @param tform subterm types
//' @param a_n starting values
//' @param dfc covariate column numbers
//' @param x_all covariate matrix
//' @param fir first term number
//' @param der_iden subterm number for derivative tests
//' @param modelform model string
//' @param Control control list
//' @param KeepConstant vector of parameters to keep constant
//' @param term_tot total number of terms
//'
//' @return LogLik_Poisson output : Log-likelihood of optimum, first derivative of log-likelihood, second derivative matrix, parameter list, standard deviation estimate, AIC, deviance, model information
// [[Rcpp::export]]
List poisson_transition(NumericMatrix dfe, IntegerVector Term_n, StringVector tform, NumericVector a_n,IntegerVector dfc,NumericMatrix x_all, int fir, int der_iden,string modelform, List Control, IntegerVector KeepConstant, int term_tot){
    //----------------------------------------------------------------------------------------------------------------//
    const Map<MatrixXd> PyrC(as<Map<MatrixXd> >(dfe));
    //
    bool change_all = Control["change_all"];
    double lr = Control["lr"];
    int maxiter = Control["maxiter"];
    int halfmax = Control["halfmax"];
    double epsilon = Control["epsilon"];
    double dbeta_cap = Control["dbeta_max"];
    double abs_max = Control["abs_max"];
    double dose_abs_max = Control["dose_abs_max"];
    double deriv_epsilon =Control["deriv_epsilon"];
    bool verbose = Control["verbose"];
    bool debugging = FALSE;
    // calculates the poisson regression
    //----------------------------------------------------------------------------------------------------------------//
    List res = LogLik_Poisson(PyrC,Term_n, tform, a_n, x_all, dfc,fir, der_iden,modelform, lr, maxiter, halfmax, epsilon, dbeta_cap, abs_max,dose_abs_max, deriv_epsilon, change_all,verbose, debugging, KeepConstant, term_tot);
    //----------------------------------------------------------------------------------------------------------------//
    return res;
}

//' Interface between R code and the Cox PH stress tester
//' \code{Stress_Test} Called directly from R, Defines the verbosity and tie method variables, Calls the calculation function
//' @param Term_n Term numbers
//' @param tform subterm types
//' @param a_n starting values
//' @param dfc covariate column numbers
//' @param x_all covariate matrix
//' @param fir first term number
//' @param der_iden subterm number for derivative tests
//' @param modelform model string
//' @param Control control list
//' @param df_groups time and event matrix
//' @param tu event times
//' @param KeepConstant vector of parameters to keep constant
//' @param term_tot total number of terms
//' @param test_point string vector of functions to test further
//'
//' @return NULL
// [[Rcpp::export]]
void Stress_Test(IntegerVector Term_n, StringVector tform, NumericVector a_n,IntegerVector dfc,NumericMatrix x_all, int fir, int der_iden,string modelform, List Control, NumericMatrix df_groups, NumericVector tu, IntegerVector KeepConstant, int term_tot, StringVector test_point){
    bool change_all = Control["change_all"];
    bool verbose = FALSE;
    bool debugging = FALSE;
    double lr = Control["lr"];
    int maxiter = Control["maxiter"];
    int halfmax = Control["halfmax"];
    double epsilon = Control["epsilon"];
    double dbeta_cap = Control["dbeta_max"];
    double abs_max = Control["abs_max"];
    double dose_abs_max = Control["dose_abs_max"];
    double deriv_epsilon =Control["deriv_epsilon"];
    string ties_method =Control["ties"];

    Stress_Run(Term_n, tform, a_n, x_all, dfc,fir, der_iden,modelform, lr, maxiter, halfmax, epsilon, dbeta_cap, abs_max,dose_abs_max, deriv_epsilon, df_groups, tu, change_all,verbose, debugging, KeepConstant, term_tot,test_point, ties_method );
    return;
}

//' Interface between R code and the Cox PH null model function
//' \code{Stress_Test} Called directly from R, Defines the control variables and calls the calculation function
//' @param ntime number of unique event times used
//' @param Control control list
//' @param df_groups time and event matrix
//' @param tu event times
//'
//' @return LogLike_Cox_PH_null output : Log-likelihood of optimum, AIC
// [[Rcpp::export]]
List cox_ph_null(int ntime, List Control, NumericMatrix df_groups, NumericVector tu){
    //----------------------------------------------------------------------------------------------------------------//
    //
    // Calculates null model
    //
    // Converts from Rcpp types to efficient Eigen types
    bool verbose = Control["verbose"];
    string ties_method =Control["ties"];
    //
    //----------------------------------------------------------------------------------------------------------------//
    List res = LogLik_Cox_PH_null(ntime, df_groups, tu, verbose, ties_method);
    //----------------------------------------------------------------------------------------------------------------//
    return res;
}

//' Interface between R code and the reference risk calculation
//' \code{cox_ph_risk_sub} Called directly from R, Defines the control variables and calls the calculation function
//' @param Term_n Term numbers
//' @param tform subterm types
//' @param a_n starting values
//' @param dfc covariate column numbers
//' @param x_all covariate matrix
//' @param fir first term number
//' @param modelform model string
//' @param Control control list
//' @param term_tot total number of terms
//'
//' @return RISK_SUBSET output: Risk at the reference
// [[Rcpp::export]]
NumericVector cox_ph_risk_sub(IntegerVector Term_n, StringVector tform, NumericVector a_n,IntegerVector dfc,NumericMatrix x_all, int fir,string modelform, List Control, int term_tot){
    //
    bool verbose = Control["verbose"];
    bool debugging = FALSE;
    NumericVector res;
    //----------------------------------------------------------------------------------------------------------------//
    // calculates risk for a reference vector
    res = RISK_SUBSET(Term_n, tform, a_n, x_all, dfc,fir,modelform,verbose, debugging, term_tot);
    return res;
}
