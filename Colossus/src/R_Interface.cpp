#include <RcppEigen.h>
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
    int double_step = Control["double_step"];
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
    int nthreads = Control["Ncores"];
    //
    // Performs regression
    //----------------------------------------------------------------------------------------------------------------//
    List res = LogLik_Cox_PH(Term_n, tform, a_n, x_all, dfc,fir, der_iden,modelform, lr, maxiter, halfmax, epsilon, dbeta_cap, abs_max,dose_abs_max, deriv_epsilon, df_groups, tu, double_step, change_all,verbose, debugging, KeepConstant, term_tot, ties_method, nthreads);
    //----------------------------------------------------------------------------------------------------------------//
    return res;
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
//' @param cens_vec censoring weight list
//' @param KeepConstant vector of parameters to keep constant
//' @param term_tot total number of terms
//'
//' @return LogLik_Cox_PH output : Log-likelihood of optimum, first derivative of log-likelihood, second derivative matrix, parameter list, standard deviation estimate, AIC, model information
// [[Rcpp::export]]
List cox_ph_transition_CR(IntegerVector Term_n, StringVector tform, NumericVector a_n,IntegerVector dfc,NumericMatrix x_all, int fir, int der_iden,string modelform, List Control, NumericMatrix df_groups, NumericVector tu, NumericVector cens_vec, IntegerVector KeepConstant, int term_tot){
    bool change_all = Control["change_all"];
    int double_step = Control["double_step"];
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
    int nthreads = Control["Ncores"];
	double cens_thres = Control["cens_thres"];
    //
	const Map<VectorXd> cens_weight(as<Map<VectorXd> >(cens_vec));
	//
    // Performs regression
    //----------------------------------------------------------------------------------------------------------------//
    List res = LogLik_Cox_PH_CR(Term_n, tform, a_n, x_all, dfc,fir, der_iden,modelform, lr, maxiter, halfmax, epsilon, dbeta_cap, abs_max,dose_abs_max, deriv_epsilon, df_groups, tu, cens_weight, cens_thres, double_step, change_all,verbose, debugging, KeepConstant, term_tot, ties_method, nthreads);
    //----------------------------------------------------------------------------------------------------------------//
    return res;
}

//' Interface between R code and the Cox PH calculation
//' \code{cox_ph_transition_single} Called directly from R, Defines the control variables and calls the function which only calculates the log-likelihood
//' @param Term_n Term numbers
//' @param tform subterm types
//' @param a_n starting values
//' @param dfc covariate column numbers
//' @param x_all covariate matrix
//' @param fir first term number
//' @param modelform model string
//' @param Control control list
//' @param df_groups time and event matrix
//' @param tu event times
//' @param KeepConstant vector of parameters to keep constant
//' @param term_tot total number of terms
//'
//' @return LogLik_Cox_PH output : Log-likelihood of optimum, first derivative of log-likelihood, second derivative matrix, parameter list, standard deviation estimate, AIC, model information
// [[Rcpp::export]]
List cox_ph_transition_single(IntegerVector Term_n, StringVector tform, NumericVector a_n,IntegerVector dfc,NumericMatrix x_all, int fir,string modelform, List Control, NumericMatrix df_groups, NumericVector tu, IntegerVector KeepConstant, int term_tot){
    bool change_all = Control["change_all"];
    bool verbose = Control["verbose"];
    bool debugging = FALSE;
    int maxiter = Control["maxiter"];
    int halfmax = Control["halfmax"];
    string ties_method =Control["ties"];
    int nthreads = Control["Ncores"];
    //
    // Performs regression
    //----------------------------------------------------------------------------------------------------------------//
    List res = LogLik_Cox_PH_Single(Term_n, tform, a_n, x_all, dfc,fir,modelform, df_groups, tu,verbose, debugging, KeepConstant, term_tot, ties_method, nthreads);
    //----------------------------------------------------------------------------------------------------------------//
    return res;
}

//' Interface between R code and the Cox PH regression for basic model
//' \code{cox_ph_transition_basic} Called directly from R, Defines the control variables and calls the regression function
//' @param a_n starting values
//' @param dfc covariate column numbers
//' @param x_all covariate matrix
//' @param der_iden subterm number for derivative tests
//' @param Control control list
//' @param df_groups time and event matrix
//' @param tu event times
//' @param KeepConstant vector of parameters to keep constant
//'
//' @return LogLik_Cox_PH output : Log-likelihood of optimum, first derivative of log-likelihood, second derivative matrix, parameter list, standard deviation estimate, AIC, model information
// [[Rcpp::export]]
List cox_ph_transition_basic( NumericVector a_n,IntegerVector dfc,NumericMatrix x_all, int der_iden, List Control, NumericMatrix df_groups, NumericVector tu, IntegerVector KeepConstant){
    bool change_all = Control["change_all"];
    int double_step = Control["double_step"];
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
    int nthreads = Control["Ncores"];
    // Performs regression
    //----------------------------------------------------------------------------------------------------------------//
    List res = LogLik_Cox_PH_basic(a_n, x_all, dfc, der_iden, lr, maxiter, halfmax, epsilon, dbeta_cap, abs_max,dose_abs_max, deriv_epsilon, df_groups, tu, double_step, change_all,verbose, debugging, KeepConstant, ties_method, nthreads);
    //----------------------------------------------------------------------------------------------------------------//
    return res;
}

//' Interface between R code and the Cox PH regression with STRATA
//' \code{cox_ph_STRATA} Called directly from R, Defines the control variables and calls the regression function
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
//' @param STRATA_vals vector of strata identifier values
//'
//' @return LogLik_Cox_PH output : Log-likelihood of optimum, first derivative of log-likelihood, second derivative matrix, parameter list, standard deviation estimate, AIC, model information
// [[Rcpp::export]]
List cox_ph_STRATA(IntegerVector Term_n, StringVector tform, NumericVector a_n,IntegerVector dfc,NumericMatrix x_all, int fir, int der_iden,string modelform, List Control, NumericMatrix df_groups, NumericVector tu, IntegerVector KeepConstant, int term_tot, IntegerVector STRATA_vals){
    bool change_all = Control["change_all"];
    int double_step = Control["double_step"];
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
    int nthreads = Control["Ncores"];
    // Performs regression
    //----------------------------------------------------------------------------------------------------------------//
    List res = LogLik_Cox_PH_STRATA(Term_n, tform, a_n, x_all, dfc,fir, der_iden,modelform, lr, maxiter, halfmax, epsilon, dbeta_cap, abs_max,dose_abs_max, deriv_epsilon, df_groups, tu, double_step, change_all,verbose, debugging, KeepConstant, term_tot, ties_method, STRATA_vals, nthreads);
    //----------------------------------------------------------------------------------------------------------------//
    return res;
}

//' Interface between R code and the Cox PH plotting
//' \code{cox_ph_plot} Called directly from R, Defines the control variables and calls the correct plotting function
//' @param Term_n Term numbers
//' @param tform subterm types
//' @param a_n optimal values
//' @param a_er optimal value standard error
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
//' @return Cox_PH_PLOT_SURV : ( baseline hazard, risk for each row) or Cox_PH_PLOT_RISK output : (covariate values, risks for each row)
// [[Rcpp::export]]
List cox_ph_plot(IntegerVector Term_n, StringVector tform, NumericVector a_n, NumericVector a_er,IntegerVector dfc,NumericMatrix x_all, int fir, int der_iden,string modelform, List Control, NumericMatrix df_groups, NumericVector tu, IntegerVector KeepConstant, int term_tot, vector<string> Plot_Type ,int uniq_v){
    bool verbose = Control["verbose"];
    bool debugging = FALSE;
    double abs_max = Control["abs_max"];
    double dose_abs_max = Control["dose_abs_max"];
    string ties_method =Control["ties"];
    List res;
    int nthreads = Control["Ncores"];
    // there are two types of plots that can be generated, survival curve and risk by covariate value
    //----------------------------------------------------------------------------------------------------------------//
    if (Plot_Type[0]=="SURV"){
        res = Cox_PH_PLOT_SURV(Term_n, tform, a_n, a_er, x_all, dfc,fir, der_iden,modelform, abs_max,dose_abs_max, df_groups, tu,verbose, debugging, KeepConstant, term_tot, nthreads);
    }else if (Plot_Type[0]=="RISK"){
        res = Cox_PH_PLOT_RISK(Term_n, tform, a_n, x_all, dfc,fir, der_iden,modelform, abs_max,dose_abs_max, df_groups, tu,verbose, debugging, KeepConstant, term_tot, uniq_v, nthreads);
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
List cox_ph_schoenfeld_transition(IntegerVector Term_n, StringVector tform, NumericVector a_n,IntegerVector dfc,NumericMatrix x_all, int fir, int der_iden,string modelform, List Control, NumericMatrix df_groups, NumericVector tu, IntegerVector KeepConstant, int term_tot){
    //
    bool verbose = Control["verbose"];
    bool debugging = FALSE;
    double abs_max = Control["abs_max"];
    double dose_abs_max = Control["dose_abs_max"];
    string ties_method =Control["ties"];
    //
    int nthreads = Control["Ncores"];
    // performs schoenfeld residual calculation
    List res = Schoenfeld_Cox_PH(Term_n, tform, a_n, x_all, dfc,fir, der_iden,modelform, abs_max,dose_abs_max, df_groups, tu,verbose, debugging, KeepConstant, term_tot, ties_method, nthreads);
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
    int double_step = Control["double_step"];
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
    int nthreads = Control["Ncores"];
    //----------------------------------------------------------------------------------------------------------------//
    List res = LogLik_Poisson(PyrC,Term_n, tform, a_n, x_all, dfc,fir, der_iden,modelform, lr, maxiter, halfmax, epsilon, dbeta_cap, abs_max,dose_abs_max, deriv_epsilon, double_step, change_all,verbose, debugging, KeepConstant, term_tot, nthreads);
    //----------------------------------------------------------------------------------------------------------------//
    return res;
}

//' Interface between R code and the poisson calculation
//' \code{poisson_transition_single} Called directly from R, Defines the control variables and calls the calculation function
//' @param dfe Matrix with person-year/event count information
//' @param Term_n Term numbers
//' @param tform subterm types
//' @param a_n starting values
//' @param dfc covariate column numbers
//' @param x_all covariate matrix
//' @param fir first term number
//' @param modelform model string
//' @param Control control list
//' @param KeepConstant vector of parameters to keep constant
//' @param term_tot total number of terms
//'
//' @return LogLik_Poisson output : Log-likelihood of optimum, first derivative of log-likelihood, second derivative matrix, parameter list, standard deviation estimate, AIC, deviance, model information
// [[Rcpp::export]]
List poisson_transition_single(NumericMatrix dfe, IntegerVector Term_n, StringVector tform, NumericVector a_n,IntegerVector dfc,NumericMatrix x_all, int fir,string modelform, List Control, IntegerVector KeepConstant, int term_tot){
    //----------------------------------------------------------------------------------------------------------------//
    const Map<MatrixXd> PyrC(as<Map<MatrixXd> >(dfe));
    //
    bool verbose = Control["verbose"];
    bool debugging = FALSE;
    // calculates the poisson regression
    int nthreads = Control["Ncores"];
    //----------------------------------------------------------------------------------------------------------------//
    List res = LogLik_Poisson_Single(PyrC,Term_n, tform, a_n, x_all, dfc,fir,modelform,verbose, debugging, KeepConstant, term_tot, nthreads);
    //----------------------------------------------------------------------------------------------------------------//
    return res;
}

//' Interface between R code and the poisson regression with strata effect
//' \code{poisson_strata_transition} Called directly from R, Defines the control variables and calls the regression function with strata effect
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
//' @param STRATA_vals vector of strata identifier values
//'
//' @return LogLik_Poisson output : Log-likelihood of optimum, first derivative of log-likelihood, second derivative matrix, parameter list, standard deviation estimate, AIC, deviance, model information
// [[Rcpp::export]]
List poisson_strata_transition(NumericMatrix dfe, IntegerVector Term_n, StringVector tform, NumericVector a_n,IntegerVector dfc,NumericMatrix x_all, int fir, int der_iden,string modelform, List Control, IntegerVector KeepConstant, int term_tot, IntegerVector STRATA_vals){
    //----------------------------------------------------------------------------------------------------------------//
    const Map<MatrixXd> PyrC(as<Map<MatrixXd> >(dfe));
    //
    bool change_all = Control["change_all"];
    bool keep_strata = Control["keep_strata"];
    int double_step = Control["double_step"];
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
    int nthreads = Control["Ncores"];
    //----------------------------------------------------------------------------------------------------------------//
    List res = LogLik_Poisson_STRATA(PyrC,Term_n, tform, a_n, x_all, dfc,fir, der_iden,modelform, lr, maxiter, halfmax, epsilon, dbeta_cap, abs_max,dose_abs_max, deriv_epsilon, double_step, change_all,verbose, debugging, KeepConstant, term_tot, STRATA_vals,keep_strata, nthreads);
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
    int double_step = Control["double_step"];
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
    int nthreads = Control["Ncores"];
    Stress_Run(Term_n, tform, a_n, x_all, dfc,fir, der_iden,modelform, lr, maxiter, halfmax, epsilon, dbeta_cap, abs_max,dose_abs_max, deriv_epsilon, df_groups, tu, double_step, change_all,verbose, debugging, KeepConstant, term_tot,test_point, ties_method , nthreads);
    return;
}

//' Interface between R code and the Cox PH null model function
//' \code{Stress_Test} Called directly from R, Defines the control variables and calls the calculation function
//' @param Control control list
//' @param df_groups time and event matrix
//' @param tu event times
//'
//' @return LogLike_Cox_PH_null output : Log-likelihood of optimum, AIC
// [[Rcpp::export]]
List cox_ph_null( List Control, NumericMatrix df_groups, NumericVector tu){
    //----------------------------------------------------------------------------------------------------------------//
    //
    // Calculates null model
    //
    // Converts from Rcpp types to efficient Eigen types
    bool verbose = Control["verbose"];
    string ties_method =Control["ties"];
    //
    int nthreads = Control["Ncores"];
    //----------------------------------------------------------------------------------------------------------------//
    List res = LogLik_Cox_PH_null( df_groups, tu, verbose, ties_method, nthreads);
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
    int nthreads = Control["Ncores"];
    //----------------------------------------------------------------------------------------------------------------//
    // calculates risk for a reference vector
    res = RISK_SUBSET(Term_n, tform, a_n, x_all, dfc,fir,modelform,verbose, debugging, term_tot, nthreads);
    return res;
}


//' Generates csv file with time-dependent columns
//' \code{Write_Time_Dep} Called directly from R, Defines a new matrix which interpolates time-dependent values on a grid
//' @param df0_Times Matrix with (starting time, ending time)
//' @param df0_dep matrix with pairs of (covariate at start, covariate at end) for each time-dependent covariate
//' @param df0_const matrix with values that are held constant
//' @param df0_event matrix with event status, zero up to the last entry for each original row
//' @param dt spacing in time
//' @param filename file to save the data to
//' @param tform vector with types of time dependent variables
//' @param tu list of event times
//' @param iscox boolean of cox formatting is used
//'
//' @return saves a dataframe to be used with time-dependent covariate analysis
// [[Rcpp::export]]
void Write_Time_Dep(const NumericMatrix df0_Times, const NumericMatrix df0_dep, const NumericMatrix df0_const, const NumericVector df0_event,double dt, string filename, StringVector tform, NumericVector tu, bool iscox){
    const Map<MatrixXd> df_Times(as<Map<MatrixXd> >(df0_Times));
    const Map<MatrixXd> df_dep(as<Map<MatrixXd> >(df0_dep));
    const Map<MatrixXd> df_const(as<Map<MatrixXd> >(df0_const));
    Rcout.precision(10); //forces higher precision numbers printed to terminal
    int nthreads = Eigen::nbThreads()-1; //stores how many threads are allocated
    if (df_dep.cols() % 2 !=0 ){
        Rcout << "Dependent columns not even" << endl;
        return;
    }
    //double epsilon=1e-2 * dt;
    int tot_covs = ceil(2 + df_dep.cols()/2 + df_const.cols() + 1);
    int max_rows = 0;
    if (iscox){
        max_rows = tu.size();
        dt = tu[1] - tu[0];
        for (int i =0; i<tu.size()-1;i++){
            if (dt > (tu[i+1] - tu[i])){
                dt = tu[i+1] - tu[i];
            }
        }
    } else {
        max_rows = ceil( ((df_Times.col(1).array() - df_Times.col(0).array()).array().abs().maxCoeff()) / dt);
    }
    // Rcout << tu << " " << dt << endl;
    int True_Rows=0;
    VectorXd row_store = VectorXd::Zero(tot_covs);
    MatrixXd new_row_store = MatrixXd::Zero(max_rows, tot_covs);
    //
    const static IOFormat CSVFormat(FullPrecision, DontAlignCols, ", ", "\n");
    ofstream file(filename);
    //
    int serial_0 = 0;
    int serial_1 = 0;
    // Rcout << "step 0" << endl;
    for (int i_row=0; i_row<df_Times.rows(); i_row++){
        new_row_store = MatrixXd::Zero(max_rows, tot_covs);
        if (iscox){
            // Rcout << "step 1 " << i_row << endl;
            True_Rows=0;
            serial_0 = 0;
            serial_1 = 0;
            for (int i=0;i<tu.size(); i++){
                if (df_Times.coeff(i_row,1) >= tu[i]){
                    serial_1 = i;
                }
                if (df_Times.coeff(i_row,0) > tu[i]){
                    serial_0 = i+1;
                }
            }
            True_Rows = serial_1 - serial_0 + 1;
//            Rcout << tu[serial_0] << " " << df_Times.coeff(i_row,0)  << " " << df_Times.coeff(i_row,1) << " " << tu[serial_1] << endl;
            // Rcout << "step 2 " << i_row << endl;
            #pragma omp parallel for schedule(dynamic) num_threads(nthreads)
            for (int i_inner=serial_0;i_inner<serial_1+1;i_inner++){
                VectorXd dep_temp = VectorXd::Zero(df_dep.cols()/2);
                double t0 = tu[i_inner]- dt/2;
                double t1 = tu[i_inner];
                double ratio = 0;
                if ((df_Times.coeff(i_row,1) - df_Times.coeff(i_row,0)) > 0){
                    ratio = (t1 - df_Times.coeff(i_row,0))/(df_Times.coeff(i_row,1) - df_Times.coeff(i_row,0));
                }
                string func_id = "";
                string delim = "?";
                size_t pos = 0;
                string token = "";
                double temp_tok =0;
                int gather_val=0;
                char tok_char = 'a';
                // step string is either g, l, a ,b for >=, <=, >, <
                for (int i=0; i<dep_temp.size();i++){
                    func_id = as<std::string>(tform[i]);
                    if (func_id=="lin"){
                        dep_temp[i] = ratio * df_dep.coeff(i_row,2*i+1) + (1 - ratio) * df_dep.coeff(i_row,2*i);
                    } else {
                        pos = func_id.find(delim);
                        token = func_id.substr(0, pos);
                        if (token=="step"){
                            func_id.erase(0, pos + delim.length());
                            gather_val=0;
                            while ((pos = func_id.find(delim)) != std::string::npos) {
                                token = func_id.substr(0, pos);
                                //
                                tok_char = token[token.length()-1];
//                                Rcout << tok_char << " " << gather_val << " " << stod(token) << " " << t1 << endl;
                                if (tok_char == 'g'){
                                    token.pop_back();
                                    temp_tok = stod(token);
                                    if (t1>temp_tok){
                                        gather_val = gather_val + 1;
                                    }
                                    if (t1==temp_tok){
                                        gather_val = gather_val + 1;
                                    }
                                } else if (tok_char == 'l'){
                                    token.pop_back();
                                    temp_tok = stod(token);
                                    if (t1<temp_tok){
                                        gather_val = gather_val + 1;
                                    }
                                    if (t1==temp_tok){
                                        gather_val = gather_val + 1;
                                    }
                                } else if (tok_char == 'a'){
                                    token.pop_back();
                                    temp_tok = stod(token);
                                    if (t1>temp_tok){
                                        gather_val = gather_val + 1;
                                    }
                                } else if (tok_char == 'b'){
                                    token.pop_back();
                                    temp_tok = stod(token);
                                    if (t1<temp_tok){
                                        gather_val = gather_val + 1;
                                    }
                                } else {
                                    ;
                                }
//                                Rcout << tok_char << " " << gather_val << " " << stod(token) << " " << t1 << endl;
                                //
                                func_id.erase(0, pos + delim.length());
                            }
                            dep_temp[i] = gather_val;
                        } else {
                            Rcout << func_id << " _:_ " << token << endl;
                            throw invalid_argument( "time dependent identifier is bad" );
                        }
                    }
                }
                int event0 = 0;
                if (i_inner==True_Rows-1){
//                    t1 = df_Times.coeff(i_row,1);
                    event0 = df0_event[i_row];
                }
                new_row_store.row(i_inner) << t0, t1, dep_temp.transpose(), df_const.row(i_row), event0;
            }
        } else {
            True_Rows = ceil( (df_Times.coeff(i_row,1) - df_Times.coeff(i_row,0))/dt);
            #pragma omp parallel for schedule(dynamic) num_threads(nthreads)
            for (int i_inner=0;i_inner<True_Rows;i_inner++){
                VectorXd dep_temp = VectorXd::Zero(df_dep.cols()/2);
                double ratio = (i_inner+0.5)/True_Rows;
                double t0 = df_Times.coeff(i_row,0) + i_inner * dt;
                double t1 = t0 + dt;
                string func_id = "";
                string delim = "?";
                size_t pos = 0;
                string token = "";
                double temp_tok =0;
                int gather_val=0;
                char tok_char = 'a';
                for (int i=0; i<dep_temp.size();i++){
                    func_id = as<std::string>(tform[i]);
                    if (func_id=="lin"){
                        dep_temp[i] = ratio * df_dep.coeff(i_row,2*i+1) + (1 - ratio) * df_dep.coeff(i_row,2*i);
                    } else {
                        pos = func_id.find(delim);
                        token = func_id.substr(0, pos);
                        if (token=="step"){
                            func_id.erase(0, pos + delim.length());
                            gather_val=0;
                            while ((pos = func_id.find(delim)) != std::string::npos) {
                                token = func_id.substr(0, pos);
                                //
                                tok_char = token[token.length()-1];
                                if (tok_char == 'u'){
                                    token.pop_back();
                                    temp_tok = stod(token);
                                    if (t0>=temp_tok){
                                        gather_val = gather_val + 1;
                                    }
                                } else if (tok_char == 'l'){
                                    token.pop_back();
                                    temp_tok = stod(token);
                                    if (t1<=temp_tok){
                                        gather_val = gather_val + 1;
                                    }
                                } else {
                                    ;
                                }
                                //
                                func_id.erase(0, pos + delim.length());
                            }
                            dep_temp[i] = gather_val;
                        }
                    }
                }
                // Rcout << "step 3 " << i_row << endl;
                int event0 = 0;
                if (i_inner==True_Rows-1){
                    t1 = df_Times.coeff(i_row,1);
                    event0 = df0_event[i_row];
                }
                new_row_store.row(i_inner) << t0 + (t1-t0)*1e-1, t1, dep_temp, df_const.row(i_row), event0;
            }
        }
        // Rcout << "step 4: " << True_Rows << " " << tot_covs << endl;
        if (file.is_open()){
            file << new_row_store.block(0,0,True_Rows,tot_covs).format(CSVFormat);
            file << "\n";
        }
        // Rcout << "step 5" << endl;
    }
    if (file.is_open()){
        file.close();
    }
}

//' Generates factored columns in parallel
//' \code{Gen_Fac_Par} Called directly from R, returns a matrix with factored columns
//' @param df0 Matrix with columns to factor, assumed to be numeric
//' @param vals list of values for each column, single continuous list
//' @param cols list of column identifiers, single continuous list
//' @param nthreads number of threads to use
//'
//' @return saves a dataframe to be used with time-dependent covariate analysis
// [[Rcpp::export]]
NumericMatrix Gen_Fac_Par(const NumericMatrix df0, const NumericVector vals, const NumericVector cols, const int nthreads){
    const Map<MatrixXd> df(as<Map<MatrixXd> >(df0));
    MatrixXd Mat_Fac = MatrixXd::Zero(df.rows(), vals.size());
    //
    #pragma omp parallel for schedule(dynamic) num_threads(1)
    for (int ijk=0;ijk<vals.size();ijk++){
        double col_c = cols[ijk];
        double val_c = vals[ijk];
        VectorXi select_ind_all = ((df.col(col_c).array() == val_c)).cast<int>(); //indices at risk
        //
        //
        int th = 1;
        visit_lambda(select_ind_all,
            [&Mat_Fac, ijk, th](double v, int i, int j) {
                if (v==th)
                    Mat_Fac(i,ijk)=1;
            });
        //
    }
    //
    return (wrap(Mat_Fac));
}





