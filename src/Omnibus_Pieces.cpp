//  Copyright 2022 - 2025, Eric Giunta and the project collaborators, Please see main R package for license and usage details

#include <RcppEigen.h>

#include "Omnibus_Pieces.h"
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
#include <functional>
#include <algorithm>

#include "Step_Calc.h"
#include "Step_Grad.h"
#include "Step_Newton.h"
#include "Calc_Repeated.h"
#include "Subterms_Risk.h"
#include "Colossus_types.h"


//  [[Rcpp::depends(RcppEigen)]]
//  [[Rcpp::plugins(openmp)]]

using std::endl;
using std::string;
using std::vector;
using std::transform;
using std::plus;
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

template <typename T> int sign(T val) {
    return (T(0) < val) - (val < T(0));
}


//' Utility function to refresh risk and subterm matrices for Cox Omnibus function
//'
//' \code{Cox_Refresh_R_TERM} Called to update matrices
//' @inheritParams CPP_template
//'
//' @return Updates matrices in place
//' @noRd
//'
void Cox_Refresh_R_TERM(const int& totalnum, const int& reqrdnum, const int& term_tot, double& dint, double& dslp, double& thres_step_max, double& step_max, const Ref<const MatrixXd>& df0, MatrixXd& T0, MatrixXd& Td0, MatrixXd& Tdd0, MatrixXd& Te, MatrixXd& R, MatrixXd& Rd, MatrixXd& Rdd, MatrixXd& Dose, MatrixXd& nonDose, MatrixXd& TTerm, MatrixXd& nonDose_LIN, MatrixXd& nonDose_PLIN, MatrixXd& nonDose_LOGLIN, MatrixXd& RdR, MatrixXd& RddR, List& model_bool) {
    const int mat_row = df0.rows();
    T0 = MatrixXd::Zero(mat_row, totalnum);  //  preallocates matrix for Term column
    R = MatrixXd::Zero(mat_row, 1);  //  preallocates matrix for Risks
    if (model_bool["basic"]) {
        Rd = MatrixXd::Zero(mat_row, reqrdnum);  //  preallocates matrix for Risk derivatives
        Rdd = MatrixXd::Zero(mat_row, reqrdnum*(reqrdnum + 1)/2);  //  preallocates matrix for Risk second derivatives
        RdR = MatrixXd::Zero(mat_row, reqrdnum);  //  preallocates matrix for Risk to derivative ratios
    } else if (model_bool["linear_err"]) {
        Rd = MatrixXd::Zero(mat_row, reqrdnum);  //  preallocates matrix for Risk derivatives
        Rdd = MatrixXd::Zero(mat_row, reqrdnum*(reqrdnum + 1)/2);  //  preallocates matrix for Risk second derivatives
        nonDose_PLIN = MatrixXd::Constant(mat_row, 1, 1.0);  //  matrix of Loglinear subterm values
        nonDose_LOGLIN = MatrixXd::Constant(mat_row, 1, 1.0);  //  matrix of Product linear subterm values
        TTerm = MatrixXd::Zero(mat_row, 1);  //  matrix of term values
        RdR = MatrixXd::Zero(mat_row, reqrdnum);  //  preallocates matrix for Risk to derivative ratios
        RddR = MatrixXd::Zero(mat_row, reqrdnum*(reqrdnum + 1)/2);  //  preallocates matrix for Risk to second derivative ratios
    } else if (model_bool["single"]) {
        Te = MatrixXd::Zero(mat_row, 1);  //  preallocates matrix for column terms used for temporary storage
        Dose = MatrixXd::Constant(mat_row, term_tot, 0.0);  //  matrix of the total dose term values
        nonDose = MatrixXd::Constant(mat_row, term_tot, 1.0);  //  matrix of the total non-dose term values
        nonDose_LIN = MatrixXd::Constant(mat_row, term_tot, 0.0);  //  matrix of Linear subterm values
        nonDose_PLIN = MatrixXd::Constant(mat_row, term_tot, 1.0);  //  matrix of Loglinear subterm values
        nonDose_LOGLIN = MatrixXd::Constant(mat_row, term_tot, 1.0);  //  matrix of Product linear subterm values
        TTerm = MatrixXd::Zero(mat_row, term_tot);  //  matrix of term values
    } else {
        Td0 = MatrixXd::Zero(mat_row, reqrdnum);  //  preallocates matrix for Term derivative columns
        Tdd0 = MatrixXd::Zero(mat_row, reqrdnum*(reqrdnum + 1)/2);  //  preallocates matrix for Term second derivative columns
        R = MatrixXd::Zero(mat_row, 1);  //  preallocates matrix for Risks
        Rd = MatrixXd::Zero(mat_row, reqrdnum);  //  preallocates matrix for Risk derivatives
        Rdd = MatrixXd::Zero(mat_row, reqrdnum*(reqrdnum + 1)/2);  //  preallocates matrix for Risk second derivatives
        Dose = MatrixXd::Constant(mat_row, term_tot, 0.0);  //  matrix of the total dose term values
        nonDose = MatrixXd::Constant(mat_row, term_tot, 1.0);  //  matrix of the total non-dose term values
        nonDose_LIN = MatrixXd::Constant(mat_row, term_tot, 0.0);  //  matrix of Linear subterm values
        nonDose_PLIN = MatrixXd::Constant(mat_row, term_tot, 1.0);  //  matrix of Loglinear subterm values
        nonDose_LOGLIN = MatrixXd::Constant(mat_row, term_tot, 1.0);  //  matrix of Product linear subterm values
        TTerm = MatrixXd::Zero(mat_row, term_tot);  //  matrix of term values
        dint = thres_step_max;  //  the amount of change used to calculate derivatives in threshold paramters
        dslp = step_max;
        RdR = MatrixXd::Zero(mat_row, reqrdnum);  //  preallocates matrix for Risk to derivative ratios
        RddR = MatrixXd::Zero(mat_row, reqrdnum*(reqrdnum + 1)/2);  //  preallocates matrix for Risk to second derivative ratios
    }
    return;
}

//' Utility function to refresh side matrices for Cox Omnibus
//'
//' \code{Cox_Refresh_R_SIDES} Called to fresh repeated sum calculation matrices
//' @inheritParams CPP_template
//'
//' @return Updates matrices in place: risk storage matrices
//' @noRd
//'
void Cox_Refresh_R_SIDES(const int& reqrdnum, const int& ntime, MatrixXd& Rls1, MatrixXd& Rls2, MatrixXd& Rls3, MatrixXd& Lls1, MatrixXd& Lls2, MatrixXd& Lls3, NumericVector& Strata_vals, List& model_bool) {
    if (model_bool["strata"]) {
        Rls1 = MatrixXd::Zero(ntime, Strata_vals.size());  //  precomputes a series of sums used frequently in the log-liklihood calculations
        Lls1 = MatrixXd::Zero(ntime, Strata_vals.size());  //  the log-likelihood calculation has a Right and Left sum used
        if (!model_bool["single"]) {
            Rls2 = MatrixXd::Zero(ntime, reqrdnum*Strata_vals.size());  //  many are repeated due to the same risk groups and derivatives being used at mulitple points
            Lls2 = MatrixXd::Zero(ntime, reqrdnum*Strata_vals.size());
            Rls3 = MatrixXd::Zero(ntime, reqrdnum*(reqrdnum + 1)/2*Strata_vals.size());  //  sum and its derivatives are precomputed
            Lls3 = MatrixXd::Zero(ntime, reqrdnum*(reqrdnum + 1)/2*Strata_vals.size());
        }
    } else {
        Rls1 = MatrixXd::Zero(ntime, 1);  //  precomputes a series of sums used frequently in the log-liklihood calculations
        Lls1 = MatrixXd::Zero(ntime, 1);  //  the log-likelihood calculation has a Right and Left sum used
        if (!model_bool["single"]) {
            Rls2 = MatrixXd::Zero(ntime, reqrdnum);  //  many are repeated due to the same risk groups and derivatives being used at mulitple points
            Lls2 = MatrixXd::Zero(ntime, reqrdnum);
            Rls3 = MatrixXd::Zero(ntime, reqrdnum*(reqrdnum + 1)/2);  //  sum and its derivatives are precomputed
            Lls3 = MatrixXd::Zero(ntime, reqrdnum*(reqrdnum + 1)/2);
        }
    }
    return;
}



//' Utility function to perform calculation of terms and risks for Cox Omnibus
//'
//' \code{Cox_Term_Risk_Calc} Called to perform repeated term and risk calculations
//' @inheritParams CPP_template
//'
//' @return Updates matrices in place: risk storage matrices
//' @noRd
//'
void Cox_Term_Risk_Calc(string modelform, const StringVector& tform, const IntegerVector& term_n, const int& totalnum, const int& fir, const IntegerVector& dfc, int term_tot, MatrixXd& T0, MatrixXd& Td0, MatrixXd& Tdd0, MatrixXd& Te, MatrixXd& R, MatrixXd& Rd, MatrixXd& Rdd, MatrixXd& Dose, MatrixXd& nonDose, VectorXd beta_0, const Ref<const MatrixXd>& df0, const double& dint, const double& dslp, MatrixXd& TTerm, MatrixXd& nonDose_LIN, MatrixXd& nonDose_PLIN, MatrixXd& nonDose_LOGLIN, MatrixXd& RdR, MatrixXd& RddR, const int& nthreads, const IntegerVector& KeepConstant, int verbose, List& model_bool, const double gmix_theta, const IntegerVector& gmix_term) {
    if (model_bool["basic"]) {
        //  Calculates the subterm and term values
        Make_subterms_Basic(totalnum, dfc, T0, beta_0, df0, nthreads);
        //  Calculates the risk for each row
        Make_Risks_Basic(totalnum, T0, R, Rd, Rdd, RdR, nthreads, df0, dfc, KeepConstant);
        //  Removes infinite values
        RdR = (RdR.array().isFinite()).select(RdR, 0);
        //
    } else if (model_bool["linear_err"]) {
        Make_subterms_Linear_ERR(totalnum, tform, dfc, nonDose_PLIN, nonDose_LOGLIN, beta_0, df0, nthreads, KeepConstant);
        //
        Make_Risks_Linear_ERR(tform, dfc, df0, totalnum, R, Rd, Rdd, nonDose_PLIN, nonDose_LOGLIN, RdR, RddR, nthreads, KeepConstant);
        RdR = (RdR.array().isFinite()).select(RdR, 0);
        RddR = (RddR.array().isFinite()).select(RddR, 0);
        //
        TTerm = R.col(0).array();
    } else if (model_bool["single"]) {
        //  Calculates the subterm and term values
        Make_subterms_Single(totalnum, term_n, tform, dfc, fir, T0, Dose, nonDose, TTerm, nonDose_LIN, nonDose_PLIN, nonDose_LOGLIN, beta_0, df0, nthreads, KeepConstant);
        //
        //  Calculates the risk for each row
        Make_Risks_Single(modelform, tform, term_n, totalnum, fir, T0, Te, R, Dose, nonDose, TTerm, nonDose_LIN, nonDose_PLIN, nonDose_LOGLIN, nthreads, KeepConstant, gmix_theta, gmix_term);
        //
    } else if (model_bool["gradient"]) {
        //
        Make_subterms_Gradient(totalnum, term_n, tform, dfc, fir, T0, Td0, Dose, nonDose, TTerm, nonDose_LIN, nonDose_PLIN, nonDose_LOGLIN, beta_0, df0, dint, dslp, nthreads, KeepConstant);
        Make_Risks_Gradient(modelform, tform, term_n, totalnum, fir, T0, Td0, Te, R, Rd, Dose, nonDose, TTerm, nonDose_LIN, nonDose_PLIN, nonDose_LOGLIN, RdR, nthreads, KeepConstant, gmix_theta, gmix_term);
        //
        RdR = (RdR.array().isFinite()).select(RdR, 0);
        //
    } else {
        //
        //  Calculates the subterm and term values
        //
        Make_subterms(totalnum, term_n, tform, dfc, fir, T0, Td0, Tdd0, Dose, nonDose, TTerm, nonDose_LIN, nonDose_PLIN, nonDose_LOGLIN, beta_0, df0, dint, dslp, nthreads, KeepConstant);
        //  ---------------------------------------------------------
        //  Prints off a series of calculations to check at what point values are changing
        //  ---------------------------------------------------------
        //  Calculates the risk for each row
        Make_Risks(modelform, tform, term_n, totalnum, fir, T0, Td0, Tdd0, Te, R, Rd, Rdd, Dose, nonDose, TTerm, nonDose_LIN, nonDose_PLIN, nonDose_LOGLIN, RdR, RddR, nthreads, KeepConstant, gmix_theta, gmix_term);
        //
        //  Removes infinite values
        RdR = (RdR.array().isFinite()).select(RdR, 0);
        RddR = (RddR.array().isFinite()).select(RddR, 0);
    }
    return;
}

//' Utility function to perform calculation of Repeated Calculations and Log-Likelihood for Cox Omnibus
//'
//' \code{Cox_Side_LL_Calc} Called to perform repeated term and risk calculations
//' @inheritParams CPP_template
//'
//' @return Updates matrices in place: risk storage matrices
//' @noRd
//'
void Cox_Side_LL_Calc(const int& reqrdnum, const int& ntime, const StringVector& tform, const IntegerMatrix& RiskFail, const vector<vector<int> >& RiskPairs, const vector<vector<vector<int> > >& RiskPairs_Strata, const int& totalnum, const int& fir, MatrixXd& R, MatrixXd& Rd, MatrixXd& Rdd, MatrixXd& Rls1, MatrixXd& Rls2, MatrixXd& Rls3, MatrixXd& Lls1, MatrixXd& Lls2, MatrixXd& Lls3, const VectorXd& cens_weight, NumericVector& Strata_vals, VectorXd beta_0, MatrixXd& RdR, MatrixXd& RddR, vector<double>& Ll, vector<double>& Lld, vector<double>& Lldd, const int& nthreads, const IntegerVector& KeepConstant, string ties_method, int verbose, List& model_bool, int iter_stop) {
    //  Calculates the side sum terms used
    if (model_bool["outcome_prob"]) {
        Calculate_Sides_PO(model_bool, RiskFail, RiskPairs, totalnum, ntime, R, Rd, Rdd, Rls1, Rls2, Rls3, Lls1, Lls2, Lls3, cens_weight, nthreads, KeepConstant);
    } else if (model_bool["strata"]) {
        if (model_bool["cr"]) {
            Calculate_Sides_Strata_CR(model_bool, RiskFail, RiskPairs_Strata, totalnum, ntime, R, Rd, Rdd, Rls1, Rls2, Rls3, Lls1, Lls2, Lls3, cens_weight, nthreads, Strata_vals, KeepConstant);
        } else {
            Calculate_Sides_Strata(model_bool, RiskFail, RiskPairs_Strata, totalnum, ntime, R, Rd, Rdd, Rls1, Rls2, Rls3, Lls1, Lls2, Lls3, nthreads, Strata_vals, KeepConstant);
        }
    } else if (model_bool["cr"]) {
        Calculate_Sides_CR(model_bool, RiskFail, RiskPairs, totalnum, ntime, R, Rd, Rdd, Rls1, Rls2, Rls3, Lls1, Lls2, Lls3, cens_weight, nthreads, KeepConstant);
    } else {
        Calculate_Sides(model_bool, RiskFail, RiskPairs, totalnum, ntime, R, Rd, Rdd, Rls1, Rls2, Rls3, Lls1, Lls2, Lls3, nthreads, KeepConstant);
    }
    //  Calculates log-likelihood
    fill(Ll.begin(), Ll.end(), 0.0);
    if (!model_bool["single"]) {
        fill(Lld.begin(), Lld.end(), 0.0);
        if (!model_bool["gradient"]) {
            fill(Lldd.begin(), Lldd.end(), 0.0);
        }
    }
    if (model_bool["outcome_prob"]) {
        Calc_LogLik_PO(model_bool, nthreads, RiskFail, RiskPairs, totalnum, ntime, R, Rd, Rdd, RdR, RddR, Rls1, Rls2, Rls3, Lls1, Lls2, Lls3, cens_weight, Ll, Lld, Lldd, ties_method, KeepConstant);
    } else if (model_bool["strata"]) {
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
    //
    if (model_bool["single"]) {
        iter_stop = 1;
    }
}

//' Utility function to print likelihood and derivatives
//'
//' \code{Print_LL} Called to print likelihood and derivatives
//' @inheritParams CPP_template
//'
//' @return Noting
//' @noRd
//'
void Print_LL(const int& reqrdnum, const int& totalnum, VectorXd beta_0, vector<double>& Ll, vector<double>& Lld, vector<double>& Lldd, int verbose, List& model_bool) {
    if (verbose >= 4) {
        Rcout << "C++ Note: df101 ";  //  prints the log-likelihoods
        for (int ij = 0; ij < reqrdnum; ij++) {
            Rcout << Ll[ij] << " ";
        }
        Rcout << " " << endl;
        if (!model_bool["single"]) {
            Rcout << "C++ Note: df102 ";  //  prints the first derivatives
            for (int ij = 0; ij < reqrdnum; ij++) {
                Rcout << Lld[ij] << " ";
            }
            Rcout << " " << endl;
            if (!model_bool["gradient"]) {
               Rcout << "C++ Note: df103 ";  //  prints the second derivatives
               for (int ij = 0; ij < reqrdnum; ij++) {
                   Rcout << Lldd[ij*reqrdnum+ij] << " ";
               }
               Rcout << " " << endl;
               Rcout << "C++ Note: ALL df103 ";  //  prints the second derivatives
               for (int ijk = 0; ijk < reqrdnum*reqrdnum; ijk++) {
                   Rcout << Lldd[ijk] << " ";
               }
                Rcout << " " << endl;
            }
        }
        if (!model_bool["null"]) {
            Rcout << "C++ Note: df104 ";  //  prints parameter values
            for (int ij = 0; ij < totalnum; ij++) {
                Rcout << beta_0[ij] << " ";
            }
            Rcout << " " << endl;
        }
    }
}

//' Utility function to print likelihood and derivatives
//'
//' \code{Print_LL_Background} Called to print likelihood and derivatives for background model
//' @inheritParams CPP_template
//'
//' @return Noting
//' @noRd
//'
void Print_LL_Background(const int& reqrdnum, const int& totalnum, const int& group_num, const int& reqrdcond, vector<double> strata_odds, vector<double>& LldOdds, vector<double>& LlddOdds, vector<double>& LlddOddsBeta, int verbose, List& model_bool) {
    if (verbose >= 4) {
        if (!model_bool["single"]) {
            Rcout << "C++ Note: df105 ";  //  prints the first derivatives
            for (int ij = 0; ij < reqrdcond; ij++) {
                Rcout << LldOdds[ij] << " ";
            }
            Rcout << " " << endl;
            if (!model_bool["gradient"]) {
               Rcout << "C++ Note: df106 ";  //  prints the second derivatives
               for (int ij = 0; ij < reqrdcond; ij++) {
                   Rcout << LlddOdds[ij] << " ";
               }
               Rcout << " " << endl;
               Rcout << "C++ Note: df107 ";  //  prints the second derivatives
               for (int ijk = 0; ijk < reqrdnum*reqrdcond; ijk++) {
                   Rcout << LlddOddsBeta[ijk] << " ";
               }
               Rcout << " " << endl;
            }
        }
        if (!model_bool["null"]) {
            Rcout << "C++ Note: df108 ";  //  prints parameter values
            for (int ij = 0; ij < group_num; ij++) {
                Rcout << strata_odds[ij] << " ";
            }
            Rcout << " " << endl;
        }
    }
}

//' Utility function to perform calculation of terms and risks for Poisson Omnibus
//'
//' \code{Pois_Term_Risk_Calc} Called to perform repeated term and risk calculations
//' @inheritParams CPP_template
//'
//' @return Updates matrices in place: risk storage matrices
//' @noRd
//'
void Pois_Term_Risk_Calc(string modelform, const StringVector& tform, const IntegerVector& term_n, const int& totalnum, const int& fir, const IntegerVector& dfc, int term_tot, MatrixXd& T0, MatrixXd& Td0, MatrixXd& Tdd0, MatrixXd& Te, MatrixXd& R, MatrixXd& Rd, MatrixXd& Rdd, MatrixXd& Dose, MatrixXd& nonDose, VectorXd beta_0, const Ref<const MatrixXd>& df0, const double& dint, const double& dslp, MatrixXd& TTerm, MatrixXd& nonDose_LIN, MatrixXd& nonDose_PLIN, MatrixXd& nonDose_LOGLIN, MatrixXd& RdR, MatrixXd& RddR, const VectorXd& s_weights, const int& nthreads, const IntegerVector& KeepConstant, int verbose, List& model_bool, const double gmix_theta, const IntegerVector& gmix_term) {
    int reqrdnum = totalnum - sum(KeepConstant);
    if (model_bool["single"]) {
        //  Calculates the subterm and term values
        Make_subterms_Single(totalnum, term_n, tform, dfc, fir, T0, Dose, nonDose, TTerm, nonDose_LIN, nonDose_PLIN, nonDose_LOGLIN, beta_0, df0, nthreads, KeepConstant);
        //  Calculates the risk for each row
        if (model_bool["strata"]) {
            Make_Risks_Weighted_Single(modelform, tform, term_n, totalnum, fir, s_weights, T0, Te, R, Dose, nonDose, TTerm, nonDose_LIN, nonDose_PLIN, nonDose_LOGLIN, nthreads, KeepConstant, gmix_theta, gmix_term);
        } else {
            Make_Risks_Single(modelform, tform, term_n, totalnum, fir, T0, Te, R, Dose, nonDose, TTerm, nonDose_LIN, nonDose_PLIN, nonDose_LOGLIN, nthreads, KeepConstant, gmix_theta, gmix_term);
        }
        //
    } else if (model_bool["gradient"]) {
        Make_subterms_Gradient(totalnum, term_n, tform, dfc, fir, T0, Td0, Dose, nonDose, TTerm, nonDose_LIN, nonDose_PLIN, nonDose_LOGLIN, beta_0, df0, dint, dslp, nthreads, KeepConstant);
        //
        if (model_bool["strata"]) {
            Make_Risks_Weighted_Gradient(modelform, tform, term_n, totalnum, fir, s_weights, T0, Td0, Te, R, Rd, Dose, nonDose, TTerm, nonDose_LIN, nonDose_PLIN, nonDose_LOGLIN, RdR, nthreads, KeepConstant, gmix_theta, gmix_term);
        } else {
            Make_Risks_Gradient(modelform, tform, term_n, totalnum, fir, T0, Td0, Te, R, Rd, Dose, nonDose, TTerm, nonDose_LIN, nonDose_PLIN, nonDose_LOGLIN, RdR, nthreads, KeepConstant, gmix_theta, gmix_term);
        }
        //
        RdR = (RdR.array().isFinite()).select(RdR, 0);
    } else {
        //
        //  Calculates the subterm and term values
        //
        Make_subterms(totalnum, term_n, tform, dfc, fir, T0, Td0, Tdd0, Dose, nonDose, TTerm, nonDose_LIN, nonDose_PLIN, nonDose_LOGLIN, beta_0, df0, dint, dslp, nthreads, KeepConstant);
        //  Calculates the risk for each row
        if (model_bool["strata"]) {
            Make_Risks_Weighted(modelform, tform, term_n, totalnum, fir, s_weights, T0, Td0, Tdd0, Te, R, Rd, Rdd, Dose, nonDose, TTerm, nonDose_LIN, nonDose_PLIN, nonDose_LOGLIN, RdR, RddR, nthreads, KeepConstant, gmix_theta, gmix_term);
        } else {
            Make_Risks(modelform, tform, term_n, totalnum, fir, T0, Td0, Tdd0, Te, R, Rd, Rdd, Dose, nonDose, TTerm, nonDose_LIN, nonDose_PLIN, nonDose_LOGLIN, RdR, RddR, nthreads, KeepConstant, gmix_theta, gmix_term);
        }
        //
        //  Removes infinite values
        RdR = (RdR.array().isFinite()).select(RdR, 0);
        RddR = (RddR.array().isFinite()).select(RddR, 0);
        //
        //
        if (R.minCoeff() <= 0) {
            if (verbose >= 4) { Rcout << "C++ Warning: risk mininum " << R.minCoeff() << " " << endl; }
        }
    }
    return;
}


//' Utility function to perform calculation of Log-Likelihood and Deviation for Poisson Omnibus
//'
//' \code{Pois_Dev_LL_Calc} Called to perform repeated term and risk calculations
//' @inheritParams CPP_template
//' @param dev_temp temporary storage for deviation calculation
//' @param dev model deviation
//'
//' @return Updates matrices in place: risk storage matrices
//' @noRd
//'
void Pois_Dev_LL_Calc(const int& reqrdnum, const int& totalnum, const int& fir, MatrixXd& R, MatrixXd& Rd, MatrixXd& Rdd, VectorXd beta_0, MatrixXd& RdR, MatrixXd& RddR, vector<double>& Ll, vector<double>& Lld, vector<double>& Lldd, const Ref<const MatrixXd>& PyrC, MatrixXd& dev_temp, const int& nthreads, const IntegerVector& KeepConstant, int verbose, List& model_bool, int iter_stop, double& dev) {
    fill(Ll.begin(), Ll.end(), 0.0);
    if (!model_bool["single"]) {
        fill(Lld.begin(), Lld.end(), 0.0);
        fill(Lldd.begin(), Lldd.end(), 0.0);
    }
    Poisson_LogLik(model_bool, nthreads, totalnum, PyrC, R, Rd, Rdd, RdR, RddR, Ll, Lld, Lldd, KeepConstant);
    //
    dev_temp.col(0) = PyrC.col(0).array() * R.col(0).array();
    dev_temp.col(0) = PyrC.col(1).array() * dev_temp.col(0).array().pow(- 1).array();
    dev_temp.col(0) = dev_temp.col(0).array().log().array();
    dev_temp.col(0) = PyrC.col(1).array() * dev_temp.col(0).array();
    dev_temp.col(1) = PyrC.col(1).array() - PyrC.col(0).array() * R.col(0).array();
    //
    dev_temp.col(0) = (dev_temp.col(0).array().isFinite()).select(dev_temp.col(0), 0);
    //
    dev_temp.col(0) = dev_temp.col(0).array() - dev_temp.col(1).array();
    dev_temp.col(0) = (2 * dev_temp.col(0).array()).array();
    dev_temp = (dev_temp.array().isFinite()).select(dev_temp, 0);
    dev_temp.col(0) = (R.col(0).array() < 0).select(0, dev_temp.col(0));
    dev = dev_temp.col(0).sum();  //  deviation calculation is split into steps
    //
    if (model_bool["single"]) {
        iter_stop = 1;
    }
    return;
}



//' Utility function to check if risk is valid, and if so continue
//'
//' \code{Cox_Pois_Check_Continue} Called to perform repeated risk checks
//' @inheritParams CPP_template
//'
//' @return Updates matrices in place: risk, scores, etc storage matrices
//' @noRd
//'
void Cox_Pois_Check_Continue(List& model_bool, VectorXd beta_0, vector<double>& beta_best, vector<double>& beta_c, const VectorXd& cens_weight, vector<double>& dbeta, double& dev, MatrixXd& dev_temp, const int fir, const int halfmax, double& halves, int& ind0, int& iter_stop, const IntegerVector& KeepConstant, vector<double>& Ll, double& Ll_abs_best, vector<double>& Lld, vector<double>& Lldd, MatrixXd& Lls1, MatrixXd& Lls2, MatrixXd& Lls3, const double& Lstar, const int& nthreads, const int& ntime, const Ref<const MatrixXd>& PyrC, MatrixXd& R, MatrixXd& Rd, MatrixXd& Rdd, MatrixXd& RddR, MatrixXd& RdR, const int& reqrdnum, const StringVector& tform, const IntegerMatrix& RiskFail, const vector<vector<int> >& RiskPairs, const vector<vector<vector<int> > >& RiskPairs_Strata, MatrixXd& Rls1, MatrixXd& Rls2, MatrixXd& Rls3, NumericVector& Strata_vals, const IntegerVector& term_n, const string ties_method, const int totalnum, MatrixXd& TTerm, const int verbose) {
    if ((R.minCoeff() <= 0) || (R.hasNaN())) {
        #ifdef _OPENMP
        #pragma omp parallel for schedule(dynamic) num_threads(nthreads)
        #endif
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
        if (model_bool["cox"]) {
            Cox_Side_LL_Calc(reqrdnum, ntime, tform, RiskFail, RiskPairs, RiskPairs_Strata, totalnum, fir, R, Rd, Rdd, Rls1, Rls2, Rls3, Lls1, Lls2, Lls3, cens_weight, Strata_vals, beta_0, RdR, RddR, Ll, Lld, Lldd, nthreads, KeepConstant, ties_method, verbose, model_bool, iter_stop);
        } else {
            Pois_Dev_LL_Calc(reqrdnum, totalnum, fir, R, Rd, Rdd, beta_0, RdR, RddR, Ll, Lld, Lldd, PyrC, dev_temp, nthreads, KeepConstant, verbose, model_bool, iter_stop, dev);
        }
        //
        if (model_bool["log_bound"]) {
            if (Ll[0] > Lstar) {
                //  If it has gone beyond Lstar, then this isn't the point
                iter_stop = 1;
                halves = halfmax;
            } else {
                if (Ll[ind0] <= Ll_abs_best) {  //  if a better point wasn't found, takes a half-step
                    for (int ijk = 0; ijk < totalnum; ijk++) {
                        dbeta[ijk] = dbeta[ijk] * 0.5;  //
                    }
                } else {  //  if improved, updates the best vector
                    for (int ijk = 0; ijk < totalnum; ijk++) {
                        beta_best[ijk] = beta_c[ijk];
                    }
                }
            }
        } else {
            if (Ll[ind0] <= Ll_abs_best) {  //  if a better point wasn't found, takes a half-step
                for (int ijk = 0; ijk < totalnum; ijk++) {
                    dbeta[ijk] = dbeta[ijk] * 0.5;
                }
            } else {  //  if improved, updates the best vector
                for (int ijk = 0; ijk < totalnum; ijk++) {
                    beta_best[ijk] = beta_c[ijk];
                }
            }
        }
        for (int ijk = 0; ijk < totalnum; ijk++) {
            beta_0[ijk] = beta_c[ijk];
        }
    }
    return;
}

//' Utility function to check if risk is valid, and if so continue
//'
//' \code{Cox_Pois_Check_Continue} Called to perform repeated risk checks
//' @inheritParams CPP_template
//'
//' @return Updates matrices in place: risk, scores, etc storage matrices
//' @noRd
//'
void Cox_Pois_Log_Loop(double& step_max, List& model_bool, VectorXd beta_0, vector<double>& beta_a, vector<double>& beta_c, int& bound_val, vector<double>& dbeta, const Ref<const MatrixXd>& df0, IntegerVector& dfc, double& dint, MatrixXd& Dose, double& thres_step_max, double& dslp, const int fir, const IntegerVector& gmix_term, const double& gmix_theta, int& half_check, const int halfmax, const IntegerVector& KeepConstant, vector<bool>& limit_hit, double& lr, string& modelform, MatrixXd& nonDose, MatrixXd& nonDose_LIN, MatrixXd& nonDose_LOGLIN, MatrixXd& nonDose_PLIN, const int& nthreads, MatrixXd& R, MatrixXd& Rd, MatrixXd& Rdd, MatrixXd& RddR, MatrixXd& RdR, VectorXd& s_weights, MatrixXd& T0, MatrixXd& Td0, MatrixXd& Tdd0, MatrixXd& Te, const IntegerVector& term_n, int& term_tot, StringVector& tform, const int totalnum, MatrixXd& TTerm, const int verbose) {
    while ((R.minCoeff() <= 0) || (R.hasNaN())) {
        half_check++;
        if (half_check > halfmax) {
            limit_hit[bound_val] = TRUE;
            break;
        } else {
            #ifdef _OPENMP
            #pragma omp parallel for schedule(dynamic) num_threads(nthreads)
            #endif
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
            for (int ij = 0; ij < totalnum; ij++) {
                beta_0[ij] = beta_a[ij] + lr*dbeta[ij];
                beta_c[ij] = beta_0[ij];
            }
            if (model_bool["cox"]) {
                Cox_Term_Risk_Calc(modelform, tform, term_n, totalnum, fir, dfc, term_tot, T0, Td0, Tdd0, Te, R, Rd, Rdd, Dose, nonDose, beta_0, df0, thres_step_max, step_max, TTerm, nonDose_LIN, nonDose_PLIN, nonDose_LOGLIN, RdR, RddR, nthreads, KeepConstant, verbose, model_bool, gmix_theta, gmix_term);
            } else {
                Pois_Term_Risk_Calc(modelform, tform, term_n, totalnum, fir, dfc, term_tot, T0, Td0, Tdd0, Te, R, Rd, Rdd, Dose, nonDose, beta_0, df0, dint, dslp, TTerm, nonDose_LIN, nonDose_PLIN, nonDose_LOGLIN, RdR, RddR, s_weights, nthreads, KeepConstant, verbose, model_bool, gmix_theta, gmix_term);
            }
        }
    }
    return;
}

//' Run a complete regression for a cox model
//'
//' \code{Cox_Full_Run} Called to perform one full regression
//' @inheritParams CPP_template
//'
//' @return Updates everything in place
//' @noRd
//'
List Cox_Full_Run(const int& reqrdnum, const int& ntime, const StringVector& tform, const IntegerMatrix& RiskFail, const vector<vector<int> >& RiskPairs, const vector<vector<vector<int> > >& RiskPairs_Strata, const int& totalnum, const int& fir, MatrixXd& R, MatrixXd& Rd, MatrixXd& Rdd, MatrixXd& Rls1, MatrixXd& Rls2, MatrixXd& Rls3, MatrixXd& Lls1, MatrixXd& Lls2, MatrixXd& Lls3, const VectorXd cens_weight, NumericVector& Strata_vals, VectorXd beta_0, MatrixXd& RdR, MatrixXd& RddR, vector<double>& Ll, vector<double>& Lld, vector<double>& Lldd, const int& nthreads, const IntegerVector& KeepConstant, string ties_method, int verbose, List& model_bool, int iter_stop, const int& term_tot, double& dint, double& dslp, double thres_step_max, double step_max, const Ref<const MatrixXd>& df0, MatrixXd& T0, MatrixXd& Td0, MatrixXd& Tdd0, MatrixXd& Te, MatrixXd& Dose, MatrixXd& nonDose, MatrixXd& TTerm, MatrixXd& nonDose_LIN, MatrixXd& nonDose_PLIN, MatrixXd& nonDose_LOGLIN, string modelform, const double gmix_theta, const IntegerVector& gmix_term, bool& convgd, double lr, List optim_para, int maxiter, const Ref<const MatrixXd>& Lin_Sys, const Ref<const VectorXd>& Lin_Res, const IntegerVector& term_n, const IntegerVector& dfc, const int halfmax, double epsilon, double deriv_epsilon) {
    //
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
    int iter_check = 0;  //  signal to check for convergence
    vector <double> Ll_comp(2, Ll[0]);  //  vector to compare values
    double Ll_abs_best = 10;
    vector<double> beta_abs_best(totalnum, 0.0);
    double Lld_worst = 0.0;
    //
    vector<double> dbeta(totalnum, 0.0);
    NumericVector m_g_store(reqrdnum);
    NumericVector v_beta_store(reqrdnum);
    //  Variables that are used for the risk check function shared across cox, poisson, and log bound functions
    double dev = 0.0;
    MatrixXd dev_temp = MatrixXd::Zero(1, 1);
    double Lstar = 0.0;
    MatrixXd PyrC = MatrixXd::Zero(1, 1);
    //  Calculates the subterm and term values
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
        return List::create(_["beta_0"] = wrap(beta_0), _["Deviation"] = R_NaN, _["Status"] = "FAILED_WITH_NEGATIVE_RISK", _["LogLik"] = R_NaN);
    }
    //
    //  -------------------------------------------------------------------------------------------
    //  Calculates the side sum terms used
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
        if (model_bool["basic"]) {
            Calc_Change_Basic(nthreads, totalnum, lr, step_max, Ll, Lld, Lldd, dbeta, KeepConstant);
        } else if (model_bool["gradient"]) {
                Calc_Change_Gradient(nthreads, model_bool, totalnum, optim_para, iteration, step_max, Lld, m_g_store, v_beta_store, dbeta, KeepConstant);
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
            Cox_Pois_Check_Continue(model_bool, beta_0, beta_best, beta_c, cens_weight, dbeta, dev, dev_temp, fir, halfmax, halves, ind0, iter_stop, KeepConstant, Ll, Ll_abs_best, Lld, Lldd, Lls1, Lls2, Lls3, Lstar, nthreads, ntime, PyrC, R, Rd, Rdd, RddR, RdR, reqrdnum, tform, RiskFail, RiskPairs, RiskPairs_Strata, Rls1, Rls2, Rls3, Strata_vals, term_n, ties_method, totalnum, TTerm, verbose);
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
                Cox_Pois_Check_Continue(model_bool, beta_0, beta_best, beta_c, cens_weight, dbeta, dev, dev_temp, fir, halfmax, halves, ind0, iter_stop, KeepConstant, Ll, Ll_abs_best, Lld, Lldd, Lls1, Lls2, Lls3, Lstar, nthreads, ntime, PyrC, R, Rd, Rdd, RddR, RdR, reqrdnum, tform, RiskFail, RiskPairs, RiskPairs_Strata, Rls1, Rls2, Rls3, Strata_vals, term_n, ties_method, totalnum, TTerm, verbose);
            }
            if (beta_best != beta_c) {  //  if the risk matrices aren't the optimal values, then they must be recalculated
                //  If it goes through every half step without improvement, then the maximum change needs to be decreased
                step_max = step_max*pow(0.5, halfmax);  //  reduces the step sizes
                thres_step_max = thres_step_max*pow(0.5, halfmax);
                iter_check = 1;
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
        Cox_Side_LL_Calc(reqrdnum, ntime, tform, RiskFail, RiskPairs, RiskPairs_Strata, totalnum, fir, R, Rd, Rdd, Rls1, Rls2, Rls3, Lls1, Lls2, Lls3, cens_weight, Strata_vals, beta_0, RdR, RddR, Ll, Lld, Lldd, nthreads, KeepConstant, ties_method, verbose, model_bool, iter_stop);
        Lld_worst = 0;
        for (int ij = 0; ij < reqrdnum; ij++) {
            if (abs(Lld[ij]) > Lld_worst) {
                Lld_worst = abs(Lld[ij]);
            }
        }
        if (iteration > reqrdnum) {  //  doesn't check the first several iterations for convergence
            if ((iteration % (reqrdnum)) || (iter_check == 1)) {  //  checks every set number of iterations
                iter_check = 0;
                if (Lld_worst < deriv_epsilon) {  //  ends if the derivatives are low enough
                    iter_stop = 1;
                    convgd = TRUE;
                }
                Ll_comp[1] = Ll[0];
                if (step_max < epsilon/10) {  //  if the maximum change is too low, then it ends
                    iter_stop = 1;
                }
            }
        }
    }
    if (Lld_worst < deriv_epsilon) {  //  ends if the derivatives are low enough
        iter_stop = 1;
        convgd = TRUE;
    }
    List res_list = List::create(_["LogLik"] = wrap(Ll[0]), _["beta_0"] = wrap(beta_0), _["Converged"] = convgd, _["Status"] = "PASSED");
    //  returns a list of results
    return res_list;
}

//' Run a complete regression for a cox model
//'
//' \code{Pois_Full_Run} Called to perform one full regression
//' @inheritParams CPP_template
//'
//' @return Updates everything in place
//' @noRd
//'
//
List Pois_Full_Run(const Ref<const MatrixXd>& PyrC, const int& reqrdnum, const StringVector& tform, const int& totalnum, const int& fir, MatrixXd& R, MatrixXd& Rd, MatrixXd& Rdd, const VectorXd& s_weights, VectorXd beta_0, MatrixXd& RdR, MatrixXd& RddR, vector<double>& Ll, vector<double>& Lld, vector<double>& Lldd, const int& nthreads, const IntegerVector& KeepConstant, int verbose, List& model_bool, int iter_stop, const int& term_tot, double& dint, double& dslp, double thres_step_max, double step_max, const Ref<const MatrixXd>& df0, MatrixXd& T0, MatrixXd& Td0, MatrixXd& Tdd0, MatrixXd& Te, MatrixXd& Dose, MatrixXd& nonDose, MatrixXd& TTerm, MatrixXd& nonDose_LIN, MatrixXd& nonDose_PLIN, MatrixXd& nonDose_LOGLIN, string modelform, const double gmix_theta, const IntegerVector& gmix_term, bool& convgd, double lr, List optim_para, int maxiter, const Ref<const MatrixXd>& Lin_Sys, const Ref<const VectorXd>& Lin_Res, const IntegerVector& term_n, const IntegerVector& dfc, const int halfmax, double epsilon, double deriv_epsilon) {
    //
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
    int iter_check = 0;  //  signal to check for convergence
    vector <double> Ll_comp(2, Ll[0]);  //  vector to compare values
    double Ll_abs_best = 10;
    vector<double> beta_abs_best(totalnum, 0.0);
    double Lld_worst = 0.0;
    double dev = 0.0;
    MatrixXd dev_temp = MatrixXd::Zero(PyrC.rows(), 2);
    //
    vector<double> dbeta(totalnum, 0.0);
    NumericVector m_g_store(reqrdnum);
    NumericVector v_beta_store(reqrdnum);
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
    //  Calculates the subterm and term values
    Pois_Term_Risk_Calc(modelform, tform, term_n, totalnum, fir, dfc, term_tot, T0, Td0, Tdd0, Te, R, Rd, Rdd, Dose, nonDose, beta_0, df0, dint, dslp, TTerm, nonDose_LIN, nonDose_PLIN, nonDose_LOGLIN, RdR, RddR, s_weights, nthreads, KeepConstant, verbose, model_bool, gmix_theta, gmix_term);
    if ((R.minCoeff() <= 0) || (R.hasNaN())) {
        if (verbose >= 1) {
            Rcout << "C++ Error: A non-positive risk was detected: " << R.minCoeff() << endl;
            Rcout << "C++ Warning: final failing values ";
            for (int ijk = 0; ijk < totalnum; ijk++) {
                Rcout << beta_0[ijk] << " ";
            }
            Rcout << " " << endl;
        }
        return List::create(_["beta_0"] = wrap(beta_0), _["Deviation"] = R_NaN, _["Status"] = "FAILED_WITH_NEGATIVE_RISK", _["LogLik"] = R_NaN);
    }
    //
    //  -------------------------------------------------------------------------------------------
    //  Calculates the side sum terms used
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
        //  calculates the initial change in parameter
        if (model_bool["gradient"]) {
            Calc_Change_Gradient(nthreads, model_bool, totalnum, optim_para, iteration, step_max, Lld, m_g_store, v_beta_store, dbeta, KeepConstant);
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
            Pois_Term_Risk_Calc(modelform, tform, term_n, totalnum, fir, dfc, term_tot, T0, Td0, Tdd0, Te, R, Rd, Rdd, Dose, nonDose, beta_0, df0, dint, dslp, TTerm, nonDose_LIN, nonDose_PLIN, nonDose_LOGLIN, RdR, RddR, s_weights, nthreads, KeepConstant, verbose, model_bool, gmix_theta, gmix_term);
            //
            Cox_Pois_Check_Continue(model_bool, beta_0, beta_best, beta_c, cens_weight, dbeta, dev, dev_temp, fir, halfmax, halves, ind0, iter_stop, KeepConstant, Ll, Ll_abs_best, Lld, Lldd, Lls1, Lls2, Lls3, Lstar, nthreads, ntime, PyrC, R, Rd, Rdd, RddR, RdR, reqrdnum, tform, RiskFail, RiskPairs, RiskPairs_Strata, Rls1, Rls2, Rls3, Strata_vals, term_n, ties_method, totalnum, TTerm, verbose);
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
                Pois_Term_Risk_Calc(modelform, tform, term_n, totalnum, fir, dfc, term_tot, T0, Td0, Tdd0, Te, R, Rd, Rdd, Dose, nonDose, beta_0, df0, dint, dslp, TTerm, nonDose_LIN, nonDose_PLIN, nonDose_LOGLIN, RdR, RddR, s_weights, nthreads, KeepConstant, verbose, model_bool, gmix_theta, gmix_term);
                Cox_Pois_Check_Continue(model_bool, beta_0, beta_best, beta_c, cens_weight, dbeta, dev, dev_temp, fir, halfmax, halves, ind0, iter_stop, KeepConstant, Ll, Ll_abs_best, Lld, Lldd, Lls1, Lls2, Lls3, Lstar, nthreads, ntime, PyrC, R, Rd, Rdd, RddR, RdR, reqrdnum, tform, RiskFail, RiskPairs, RiskPairs_Strata, Rls1, Rls2, Rls3, Strata_vals, term_n, ties_method, totalnum, TTerm, verbose);
            }
            if (beta_best != beta_c) {  //  if the risk matrices aren't the optimal values, then they must be recalculated
                //  If it goes through every half step without improvement, then the maximum change needs to be decreased
                step_max = step_max*pow(0.5, halfmax);  //  reduces the step sizes
                thres_step_max = thres_step_max*pow(0.5, halfmax);
                iter_check = 1;
                beta_p = beta_best;  //
                beta_a = beta_best;  //
                beta_c = beta_best;  //
                for (int ij = 0; ij < totalnum; ij++) {
                    beta_0[ij] = beta_best[ij];
                }
                Pois_Term_Risk_Calc(modelform, tform, term_n, totalnum, fir, dfc, term_tot, T0, Td0, Tdd0, Te, R, Rd, Rdd, Dose, nonDose, beta_0, df0, dint, dslp, TTerm, nonDose_LIN, nonDose_PLIN, nonDose_LOGLIN, RdR, RddR, s_weights, nthreads, KeepConstant, verbose, model_bool, gmix_theta, gmix_term);
                //
            }
        }
        Pois_Dev_LL_Calc(reqrdnum, totalnum, fir, R, Rd, Rdd, beta_0, RdR, RddR, Ll, Lld, Lldd, PyrC, dev_temp, nthreads, KeepConstant, verbose, model_bool, iter_stop, dev);
        Lld_worst = 0;
        for (int ij = 0; ij < reqrdnum; ij++) {
            if (abs(Lld[ij]) > Lld_worst) {
                Lld_worst = abs(Lld[ij]);
            }
        }
        if (iteration > reqrdnum) {  //  doesn't check the first several iterations for convergence
            if ((iteration % (reqrdnum)) || (iter_check == 1)) {  //  checks every set number of iterations
                iter_check = 0;
                if (Lld_worst < deriv_epsilon) {  //  ends if the derivatives are low enough
                    iter_stop = 1;
                    convgd = TRUE;
                }
                Ll_comp[1] = Ll[0];
                if (step_max < epsilon/10) {  //  if the maximum change is too low, then it ends
                    iter_stop = 1;
                }
            }
        }
    }
    if (Lld_worst < deriv_epsilon) {  //  ends if the derivatives are low enough
        iter_stop = 1;
        convgd = TRUE;
    }
    List res_list = List::create(_["LogLik"] = wrap(Ll[0]), _["beta_0"] = wrap(beta_0), _["Converged"] = convgd, _["Status"] = "PASSED");
    //  returns a list of results
    return res_list;
}


//' Utility function to calculate Information Matrix, from Epicure manual
//'
//' \code{Expected_Inform_Matrix_Cox} Called to update information matrix
//' @inheritParams CPP_template
//'
//' @return Updates matrices in place: information matrix
//' @noRd
//'
void Expected_Inform_Matrix_Cox(const int& nthreads, const IntegerMatrix& RiskFail, const vector<vector<int> >& RiskPairs, const int& totalnum, const int& ntime, const MatrixXd& R, const MatrixXd& Rd, const MatrixXd& RdR, vector<double>& InMa, const IntegerVector& KeepConstant) {
    int reqrdnum = totalnum - sum(KeepConstant);
    #ifdef _OPENMP
    #pragma omp declare reduction(vec_double_plus : vector<double> : \
        transform(omp_out.begin(), omp_out.end(), omp_in.begin(), omp_out.begin(), plus<double>())) \
        initializer(omp_priv = omp_orig)
    #pragma omp parallel for schedule(dynamic) num_threads(nthreads) reduction(vec_double_plus:InMa) collapse(2)
    #endif
    for (int ijk = 0; ijk < reqrdnum*(reqrdnum + 1)/2; ijk++) {  //  performs log-likelihood calculations for every derivative combination and risk group
        for (int j = 0; j < ntime; j++) {
            //  I need to collect four values, sum of risk, sum of first derivative 1, sum of first derivative 2, and the combined products of risk and derivative ratios
            //  Start by initializing the four values and getting the derivative columns
            double r_sum = 0;
            double rd_sum0 = 0;
            double rd_sum1 = 0;
            double r_rd_prod = 0;
            int ij = 0;
            int jk = ijk;
            while (jk > ij) {
                ij++;
                jk -= ij;
            }
            //  Next I need to iterate through the rows at risk
            vector<int> InGroup = RiskPairs[j];
            //  now has the grouping pairs and how many events per group
            int dj = RiskFail(j, 1)-RiskFail(j, 0) + 1;
            for (vector<double>::size_type i = 0; i < InGroup.size() - 1; i = i+2) {
                //
                r_sum += R.block(InGroup[i] - 1, 0, InGroup[i + 1]-InGroup[i] + 1, 1).sum();
                //
                rd_sum0 += Rd.block(InGroup[i] - 1, ij, InGroup[i + 1]-InGroup[i] + 1, 1).sum();
                rd_sum1 += Rd.block(InGroup[i] - 1, jk, InGroup[i + 1]-InGroup[i] + 1, 1).sum();
                //
                r_rd_prod += (RdR.block(InGroup[i] - 1, ij, InGroup[i + 1]-InGroup[i] + 1, 1).array() * RdR.block(InGroup[i] - 1, jk, InGroup[i + 1]-InGroup[i] + 1, 1).array() * R.block(InGroup[i] - 1, 0, InGroup[i + 1]-InGroup[i] + 1, 1).array()).sum();
            }
            //  now they should be combined
            double Rs1 = r_rd_prod/r_sum;
            double Ls1 = rd_sum0*rd_sum1/pow(r_sum, 2);
            InMa[ij*reqrdnum+jk] += dj*(Rs1 - Ls1);
        }
    }
    #ifdef _OPENMP
    #pragma omp parallel for schedule(dynamic) num_threads(nthreads)
    #endif
    for (int ijk = 0; ijk < reqrdnum*(reqrdnum + 1)/2; ijk++) {  //  fills second-derivative matrix
        int ij = 0;
        int jk = ijk;
        while (jk > ij) {
            ij++;
            jk -= ij;
        }
        InMa[jk*reqrdnum+ij] = InMa[ij*reqrdnum+jk];
    }
    return;
}

//' Utility function to calculate Information Matrix with strata, from Epicure manual
//'
//' \code{Expected_Inform_Matrix_Cox_Strata} Called to update information matrix with strata
//' @inheritParams CPP_template
//'
//' @return Updates matrices in place: information matrix
//' @noRd
//'
void Expected_Inform_Matrix_Cox_Strata(const int& nthreads, const IntegerMatrix& RiskFail, const vector<vector<vector<int> > >& RiskPairs_Strata, const int& totalnum, const int& ntime, const MatrixXd& R, const MatrixXd& Rd, const MatrixXd& RdR, vector<double>& InMa, NumericVector& Strata_vals, const IntegerVector& KeepConstant) {
    int reqrdnum = totalnum - sum(KeepConstant);
    #ifdef _OPENMP
    #pragma omp declare reduction(vec_double_plus : vector<double> : \
        transform(omp_out.begin(), omp_out.end(), omp_in.begin(), omp_out.begin(), plus<double>())) \
        initializer(omp_priv = omp_orig)
    #pragma omp parallel for schedule(dynamic) num_threads(nthreads) reduction(vec_double_plus:InMa) collapse(3)
    #endif
    for (int ijk = 0; ijk < reqrdnum*(reqrdnum + 1)/2; ijk++) {  //  performs log-likelihood calculations for every derivative combination and risk group
        for (int j = 0; j < ntime; j++) {
            for (int s_ij = 0; s_ij < Strata_vals.size(); s_ij++) {
                //  I need to collect four values, sum of risk, sum of first derivative 1, sum of first derivative 2, and the combined products of risk and derivative ratios
                //  Start by initializing the four values and getting the derivative columns
                double r_sum = 0;
                double rd_sum0 = 0;
                double rd_sum1 = 0;
                double r_rd_prod = 0;
                int ij = 0;
                int jk = ijk;
                while (jk > ij) {
                    ij++;
                    jk -= ij;
                }
                //  Next I need to iterate through the rows at risk
                vector<int> InGroup = RiskPairs_Strata[j][s_ij];
                if (RiskFail(j, 2*s_ij + 1) >  - 1) {
                    //  now has the grouping pairs and how many events per group
                    int dj = RiskFail(j, 2*s_ij + 1)-RiskFail(j, 2*s_ij + 0) + 1;
                    for (vector<double>::size_type i = 0; i < InGroup.size() - 1; i = i+2) {
                        //
                        r_sum += R.block(InGroup[i] - 1, 0, InGroup[i + 1]-InGroup[i] + 1, 1).sum();
                        //
                        rd_sum0 += Rd.block(InGroup[i] - 1, ij, InGroup[i + 1]-InGroup[i] + 1, 1).sum();
                        rd_sum1 += Rd.block(InGroup[i] - 1, jk, InGroup[i + 1]-InGroup[i] + 1, 1).sum();
                        //
                        r_rd_prod += (RdR.block(InGroup[i] - 1, ij, InGroup[i + 1]-InGroup[i] + 1, 1).array() * RdR.block(InGroup[i] - 1, jk, InGroup[i + 1]-InGroup[i] + 1, 1).array() * R.block(InGroup[i] - 1, 0, InGroup[i + 1]-InGroup[i] + 1, 1).array()).sum();
                    }
                    //  now they should be combined
                    double Rs1 = r_rd_prod/r_sum;
                    double Ls1 = rd_sum0*rd_sum1/pow(r_sum, 2);
                    InMa[ij*reqrdnum+jk] += dj*(Rs1 - Ls1);
                }
            }
        }
    }
    #ifdef _OPENMP
    #pragma omp parallel for schedule(dynamic) num_threads(nthreads)
    #endif
    for (int ijk = 0; ijk < reqrdnum*(reqrdnum + 1)/2; ijk++) {  //  fills second-derivative matrix
        int ij = 0;
        int jk = ijk;
        while (jk > ij) {
            ij++;
            jk -= ij;
        }
        InMa[jk*reqrdnum+ij] = InMa[ij*reqrdnum+jk];
    }
    return;
}

//' Utility function to calculate Information Matrix with competing risks, adapted from Epicure manual
//'
//' \code{Expected_Inform_Matrix_Cox_CR} Called to update information matrix with competing risks
//' @inheritParams CPP_template
//'
//' @return Updates matrices in place: information matrix
//' @noRd
//'
void Expected_Inform_Matrix_Cox_CR(const int& nthreads, const IntegerMatrix& RiskFail, const vector<vector<int> >& RiskPairs, const int& totalnum, const int& ntime, const MatrixXd& R, const MatrixXd& Rd, const MatrixXd& RdR, const VectorXd& cens_weight, vector<double>& InMa, const IntegerVector& KeepConstant) {
    int reqrdnum = totalnum - sum(KeepConstant);
    #ifdef _OPENMP
    #pragma omp declare reduction(vec_double_plus : vector<double> : \
        transform(omp_out.begin(), omp_out.end(), omp_in.begin(), omp_out.begin(), plus<double>())) \
        initializer(omp_priv = omp_orig)
    #pragma omp parallel for schedule(dynamic) num_threads(nthreads) reduction(vec_double_plus:InMa) collapse(2)
    #endif
    for (int ijk = 0; ijk < reqrdnum*(reqrdnum + 1)/2; ijk++) {  //  performs log-likelihood calculations for every derivative combination and risk group
        for (int j = 0; j < ntime; j++) {
            //  I need to collect four values, sum of risk, sum of first derivative 1, sum of first derivative 2, and the combined products of risk and derivative ratios
            //  Start by initializing the four values and getting the derivative columns
            double r_sum = 0;
            double rd_sum0 = 0;
            double rd_sum1 = 0;
            double r_rd_prod = 0;
            int ij = 0;
            int jk = ijk;
            while (jk > ij) {
                ij++;
                jk -= ij;
            }
            //  Next I need to iterate through the rows at risk
            vector<int> InGroup = RiskPairs[j];
            //  now has the grouping pairs and how many events per group
            int dj = RiskFail(j, 1)-RiskFail(j, 0) + 1;
            double cens_0 = cens_weight[RiskFail(j, 0)];
            VectorXd weighting = VectorXd::Zero(InGroup[1]-InGroup[0] + 1);
            for (vector<double>::size_type i = 0; i < InGroup.size() - 1; i = i+2) {
                if (weighting.size() != InGroup[i + 1]-InGroup[i] + 1) {
                    weighting.resize(InGroup[i + 1]-InGroup[i] + 1);
                }
                weighting.head(InGroup[i + 1]-InGroup[i] + 1) << cens_weight.segment(InGroup[i] - 1, InGroup[i + 1]-InGroup[i] + 1);
                weighting = weighting / cens_0;
                weighting = (weighting.array() < 1).select(weighting, 1);
                //
                r_sum += (R.block(InGroup[i] - 1, 0, InGroup[i + 1]-InGroup[i] + 1, 1).array() * weighting.head(InGroup[i + 1]-InGroup[i] + 1).array()).sum();
                //
                rd_sum0 += (Rd.block(InGroup[i] - 1, ij, InGroup[i + 1]-InGroup[i] + 1, 1).array() * weighting.head(InGroup[i + 1]-InGroup[i] + 1).array()).sum();
                rd_sum1 += (Rd.block(InGroup[i] - 1, jk, InGroup[i + 1]-InGroup[i] + 1, 1).array() * weighting.head(InGroup[i + 1]-InGroup[i] + 1).array()).sum();
                //
                r_rd_prod += (RdR.block(InGroup[i] - 1, ij, InGroup[i + 1]-InGroup[i] + 1, 1).array() * RdR.block(InGroup[i] - 1, jk, InGroup[i + 1]-InGroup[i] + 1, 1).array() * R.block(InGroup[i] - 1, 0, InGroup[i + 1]-InGroup[i] + 1, 1).array() * weighting.head(InGroup[i + 1]-InGroup[i] + 1).array()).sum();
            }
            //  now they should be combined
            double Rs1 = r_rd_prod/r_sum;
            double Ls1 = rd_sum0*rd_sum1/pow(r_sum, 2);
            InMa[ij*reqrdnum+jk] += dj*(Rs1 - Ls1);
        }
    }
    #ifdef _OPENMP
    #pragma omp parallel for schedule(dynamic) num_threads(nthreads)
    #endif
    for (int ijk = 0; ijk < reqrdnum*(reqrdnum + 1)/2; ijk++) {  //  fills second-derivative matrix
        int ij = 0;
        int jk = ijk;
        while (jk > ij) {
            ij++;
            jk -= ij;
        }
        InMa[jk*reqrdnum+ij] = InMa[ij*reqrdnum+jk];
    }
    return;
}

//' Utility function to calculate Information Matrix with strata and competing risks, adapted from from Epicure manual
//'
//' \code{Expected_Inform_Matrix_Cox_Strata_CR} Called to update information matrix with strata and competing risks
//' @inheritParams CPP_template
//'
//' @return Updates matrices in place: information matrix
//' @noRd
//'
void Expected_Inform_Matrix_Cox_Strata_CR(const int& nthreads, const IntegerMatrix& RiskFail, const vector<vector<vector<int> > >& RiskPairs_Strata, const int& totalnum, const int& ntime, const MatrixXd& R, const MatrixXd& Rd, const MatrixXd& RdR, const VectorXd& cens_weight, vector<double>& InMa, NumericVector& Strata_vals, const IntegerVector& KeepConstant) {
    int reqrdnum = totalnum - sum(KeepConstant);
    #ifdef _OPENMP
    #pragma omp declare reduction(vec_double_plus : vector<double> : \
        transform(omp_out.begin(), omp_out.end(), omp_in.begin(), omp_out.begin(), plus<double>())) \
        initializer(omp_priv = omp_orig)
    #pragma omp parallel for schedule(dynamic) num_threads(nthreads) reduction(vec_double_plus:InMa) collapse(3)
    #endif
    for (int ijk = 0; ijk < reqrdnum*(reqrdnum + 1)/2; ijk++) {  //  performs log-likelihood calculations for every derivative combination and risk group
        for (int j = 0; j < ntime; j++) {
            for (int s_ij = 0; s_ij < Strata_vals.size(); s_ij++) {
                //  I need to collect four values, sum of risk, sum of first derivative 1, sum of first derivative 2, and the combined products of risk and derivative ratios
                //  Start by initializing the four values and getting the derivative columns
                double r_sum = 0;
                double rd_sum0 = 0;
                double rd_sum1 = 0;
                double r_rd_prod = 0;
                int ij = 0;
                int jk = ijk;
                while (jk > ij) {
                    ij++;
                    jk -= ij;
                }
                //  Next I need to iterate through the rows at risk
                vector<int> InGroup = RiskPairs_Strata[j][s_ij];
                if (RiskFail(j, 2*s_ij + 1) >  - 1) {
                    //  now has the grouping pairs and how many events per group
                    int dj = RiskFail(j, 2*s_ij + 1)-RiskFail(j, 2*s_ij + 0) + 1;
                    double cens_0 = cens_weight[RiskFail(j, 2*s_ij)];
                    VectorXd weighting = VectorXd::Zero(InGroup[1]-InGroup[0] + 1);
                    for (vector<double>::size_type i = 0; i < InGroup.size() - 1; i = i+2) {
                        if (weighting.size() < InGroup[i + 1]-InGroup[i] + 1) {
                            weighting.resize(InGroup[i + 1]-InGroup[i] + 1);
                        }
                        weighting.head(InGroup[i + 1]-InGroup[i] + 1) << cens_weight.block(InGroup[i] - 1, 0, InGroup[i + 1]-InGroup[i] + 1, 1);
                        weighting = weighting / cens_0;
                        weighting = (weighting.array() < 1).select(weighting, 1);
                        //
                        r_sum += (R.block(InGroup[i] - 1, 0, InGroup[i + 1]-InGroup[i] + 1, 1).array() * weighting.head(InGroup[i + 1]-InGroup[i] + 1).array()).sum();
                        //
                        rd_sum0 += (Rd.block(InGroup[i] - 1, ij, InGroup[i + 1]-InGroup[i] + 1, 1).array() * weighting.head(InGroup[i + 1]-InGroup[i] + 1).array()).sum();
                        rd_sum1 += (Rd.block(InGroup[i] - 1, jk, InGroup[i + 1]-InGroup[i] + 1, 1).array() * weighting.head(InGroup[i + 1]-InGroup[i] + 1).array()).sum();
                        //
                        r_rd_prod += (RdR.block(InGroup[i] - 1, ij, InGroup[i + 1]-InGroup[i] + 1, 1).array() * RdR.block(InGroup[i] - 1, jk, InGroup[i + 1]-InGroup[i] + 1, 1).array() * R.block(InGroup[i] - 1, 0, InGroup[i + 1]-InGroup[i] + 1, 1).array() * weighting.head(InGroup[i + 1]-InGroup[i] + 1).array()).sum();
                    }
                    //  now they should be combined
                    double Rs1 = r_rd_prod/r_sum;
                    double Ls1 = rd_sum0*rd_sum1/pow(r_sum, 2);
                    InMa[ij*reqrdnum+jk] += dj*(Rs1 - Ls1);
                }
            }
        }
    }
    #ifdef _OPENMP
    #pragma omp parallel for schedule(dynamic) num_threads(nthreads)
    #endif
    for (int ijk = 0; ijk < reqrdnum*(reqrdnum + 1)/2; ijk++) {  //  fills second-derivative matrix
        int ij = 0;
        int jk = ijk;
        while (jk > ij) {
            ij++;
            jk -= ij;
        }
        InMa[jk*reqrdnum+ij] = InMa[ij*reqrdnum+jk];
    }
    return;
}


//' Utility function to calculate poisson expected information matrix
//'
//' \code{Expected_Inform_Matrix_Poisson} Called to update information matrix
//' @inheritParams CPP_template
//'
//' @return Updates matrices in place: Log-likelihood vectors/matrix
//' @noRd
//'
void Expected_Inform_Matrix_Poisson(const int& nthreads, const int& totalnum, const Ref<const MatrixXd>& PyrC, const MatrixXd& R, const MatrixXd& Rd, const MatrixXd& RdR, vector<double>& InMa, const IntegerVector& KeepConstant) {
    int reqrdnum = totalnum - sum(KeepConstant);
    #ifdef _OPENMP
    #pragma omp parallel for schedule(dynamic) num_threads(nthreads)
    #endif
    for (int ijk = 0; ijk < reqrdnum*(reqrdnum + 1)/2; ijk++) {  //  totalnum*(totalnum + 1)/2
        int ij = 0;
        int jk = ijk;
        while (jk > ij) {
            ij++;
            jk -= ij;
        }
        InMa[ij*reqrdnum+jk] = (PyrC.col(0).array() * R.col(0).array() * RdR.col(ij).array() * RdR.col(jk).array()).array().sum();
    }
    #ifdef _OPENMP
    #pragma omp parallel for schedule(dynamic) num_threads(nthreads)
    #endif
    for (int ijk = 0; ijk < reqrdnum*(reqrdnum + 1)/2; ijk++) {  //  fills second-derivative matrix
        int ij = 0;
        int jk = ijk;
        while (jk > ij) {
            ij++;
            jk -= ij;
        }
        InMa[jk*reqrdnum+ij] = InMa[ij*reqrdnum+jk];
    }
    return;
}

//' Utility function to calculate logistic model expected information matrix
//'
//' \code{Expected_Inform_Matrix_Logist} Called to update information matrix
//' @inheritParams CPP_template
//'
//' @return Updates matrices in place: Log-likelihood vectors/matrix
//' @noRd
//'
void Expected_Inform_Matrix_Logist(const int& nthreads, const int& totalnum, const Ref<const MatrixXd>& CountEvent, const MatrixXd& PdP, const MatrixXd& PnotdP, vector<double>& InMa, const IntegerVector& KeepConstant) {
    int reqrdnum = totalnum - sum(KeepConstant);
    #ifdef _OPENMP
    #pragma omp parallel for schedule(dynamic) num_threads(nthreads)
    #endif
    for (int ijk = 0; ijk < reqrdnum*(reqrdnum + 1)/2; ijk++) {  //  totalnum*(totalnum + 1)/2
        int ij = 0;
        int jk = ijk;
        while (jk > ij) {
            ij++;
            jk -= ij;
        }
        InMa[ij*reqrdnum+jk] = (CountEvent.col(1).array() * PdP.col(ij).array() * PnotdP.col(jk).array()).array().sum();
    }
    #ifdef _OPENMP
    #pragma omp parallel for schedule(dynamic) num_threads(nthreads)
    #endif
    for (int ijk = 0; ijk < reqrdnum*(reqrdnum + 1)/2; ijk++) {  //  fills second-derivative matrix
        int ij = 0;
        int jk = ijk;
        while (jk > ij) {
            ij++;
            jk -= ij;
        }
        InMa[jk*reqrdnum+ij] = InMa[ij*reqrdnum+jk];
    }
    return;
}

//' Utility function to apply a link function to a Logistic Regression Risk model
//'
//' \code{LinkCovertRP} Called to update the probability matrices
//' @inheritParams CPP_template
//'
//' @return Updates matrices in place: p, 1-p, and the associated derivatives/ratios
//' @noRd
//'
void LinkCovertRP(List& model_bool, const int& reqrdnum, const MatrixXd& R, const MatrixXd& Rd, const MatrixXd& Rdd, const MatrixXd& RdR, const MatrixXd& RddR, MatrixXd& P, MatrixXd& Pd, MatrixXd& Pdd, MatrixXd& Pnot, MatrixXd& PdP, MatrixXd& PddP, MatrixXd& PnotdP, MatrixXd& PnotddP) {
    bool odds = model_bool["odds"];
    bool ident = model_bool["ident"];
    bool loglink = model_bool["loglink"];

    if (ident == true) {
        P.col(0) = R.col(0);
        Pnot.col(0) = 1 - P.col(0).array();
        MatrixXd Ftemp = MatrixXd::Zero(R.rows(), 1);
        Ftemp.col(0) = Pnot.col(0).array().pow(-1);  //  1+r
        for (int ij=0; ij< reqrdnum; ij++) {
            Pd.col(ij) = Rd.col(ij);
            PdP.col(ij) = RdR.col(ij).array();
            PnotdP.col(ij) = Rd.col(ij).array() * Ftemp.col(0).array();
        }
        for (int ijk=0; ijk< reqrdnum*(reqrdnum + 1)/2; ijk++) {
            Pdd.col(ijk) = Rdd.col(ijk);
            PddP.col(ijk) = RddR.col(ijk).array();
            PnotddP.col(ijk) = Rdd.col(ijk).array() * Ftemp.col(0).array();
        }
    }  else if (loglink == true) {
        P.col(0) = (-1*R.col(0).array()).exp();
        Pnot.col(0) = 1 - P.col(0).array();
        //
        MatrixXd Ftemp = MatrixXd::Zero(R.rows(), 2);
        Ftemp.col(0) = P.col(0).array().pow(-1);   //
        Ftemp.col(1) = Pnot.col(0).array().pow(-1);  //
        //
        for (int ij=0; ij< reqrdnum; ij++) {
            Pd.col(ij) = -1 * Rd.col(ij).array() * P.col(0).array();
            PdP.col(ij) = -1 * Rd.col(ij).array();
            PnotdP.col(ij) = Pd.col(ij).array() * Ftemp.col(1).array();
        }
        for (int ijk=0; ijk< reqrdnum*(reqrdnum + 1)/2; ijk++) {
            int ij = 0;
            int jk = ijk;
            while (jk > ij) {
                ij++;
                jk -= ij;
            }
            PddP.col(ijk) = PdP.col(ij).array()*PdP.col(jk).array() - Rdd.col(ijk).array();
            Pdd.col(ijk) = PddP.col(ijk).array() * P.col(0).array();
            PnotddP.col(ijk) = Pdd.col(ijk).array() * Ftemp.col(1).array();
        }
    } else if (odds == true) {
        MatrixXd Ftemp = MatrixXd::Zero(R.rows(), 4);
        Ftemp.col(0) = 1 + R.col(0).array();  //  1+r
        Ftemp.col(1) = Ftemp.col(0).array().pow(-1);  // (1+r)^-1
        Ftemp.col(2) = Ftemp.col(1).array().pow(2);  //  (1+r)^-2
        Ftemp.col(3) = Ftemp.col(1).array().pow(3);  //  (1+r)^-3
        //
        P.col(0) = R.col(0).array() * Ftemp.col(1).array();
        Pnot.col(0) = 1 - P.col(0).array();
        for (int ij=0; ij< reqrdnum; ij++) {
            Pd.col(ij) = Rd.col(ij).array() * Ftemp.col(2).array();
            PdP.col(ij) = RdR.col(ij).array() * Ftemp.col(1).array();
            PnotdP.col(ij) = Rd.col(ij).array() * Ftemp.col(1).array();
        }
        for (int ijk=0; ijk< reqrdnum*(reqrdnum + 1)/2; ijk++) {
            int ij = 0;
            int jk = ijk;
            while (jk > ij) {
                ij++;
                jk -= ij;
            }
            Pdd.col(ijk) = Rdd.col(ijk).array() * Ftemp.col(2).array() - 2 * Rd.col(ij).array() * Rd.col(jk).array() * Ftemp.col(3).array();
            PddP.col(ijk) = RddR.col(ijk).array() * Ftemp.col(1).array() - 2 * RdR.col(ij).array() * Rd.col(jk).array() * Ftemp.col(2).array();
            PnotddP.col(ijk) = Rdd.col(ijk).array() * Ftemp.col(1).array() - 2 * Rd.col(ij).array() * Rd.col(jk).array() * Ftemp.col(2).array();
        }
    }
    return;
}

//' Fills out recursive vectors for matched case-control logistic regression expected information matrix
//'
//' \code{Calc_Recursive_Exp} Called to update the recursive vectors, uses model_bool list to select which vectors to update.
//'
//' @return Updates matrices in place: risk storage matrices
//' @noRd
//'
void Calc_Recursive_Exp(List& model_bool, const int& group_num, const IntegerMatrix& RiskFail, const vector<vector<int> >& RiskPairs, const int& totalnum, const int& ntime, const MatrixXd& R, const MatrixXd& Rd, vector<vector<double> >& Recur_Base, vector<vector<vector<double> > >& Recur_First, vector<vector<vector<double> > >& Recur_Second, const int& nthreads, const IntegerVector& KeepConstant) {
    int reqrdnum = 1;
    reqrdnum = totalnum - sum(KeepConstant);
    double cond_thres = model_bool["cond_thres"];
    //  We will be overwritting the Recur_Second matrix with the expected information matrix data
    #ifdef _OPENMP
    #pragma omp parallel for schedule(dynamic) num_threads(nthreads) collapse(2)
    #endif
    for (int group_ij = 0; group_ij < group_num; group_ij++) {
        for (int der_ijk = 0; der_ijk < reqrdnum*(reqrdnum + 1)/2; der_ijk++) {
            //  get the derivative column numbers
            int der_ij = 0;
            int der_jk = der_ijk;
            while (der_jk > der_ij) {
                der_ij++;
                der_jk -= der_ij;
            }
            //  we start by getting a vector of risks
            vector<double> risk_list;
            vector<double> riskd0_list;
            vector<double> riskd1_list;
            vector<int> InGroup = RiskPairs[group_ij];
            //  now has the grouping pairs and number of events
            if (InGroup.size() > 0) {
                int dj = RiskFail(group_ij, 1)-RiskFail(group_ij, 0) + 1;
                if (dj <= cond_thres) {
                    for (vector<double>::size_type i = 0; i < InGroup.size() - 1; i = i+2) {
                        int i0 = InGroup[i] - 1;
                        int i1 = InGroup[i + 1] - 1;
                        for (int i_inter = i0; i_inter <= i1; i_inter ++) {
                            risk_list.push_back(R(i_inter, 0));
                            riskd0_list.push_back(Rd(i_inter, der_ij));
                            riskd1_list.push_back(Rd(i_inter, der_jk));
                        }
                    }
                    int risk_size = risk_list.size();
                    //  we need to start filling out the recusion
                    //  we have risk_size elements and are selecting dj items, filling object B(m,n)
                    //  we start with filling out B(1, 1) to Recur_Base[group_ij][0], up to the final entry for m=1, B(1, risk_size-dj+1) to Recur_Base[group_ij][risk_size-dj]
                    double r_sum = 0;
                    int nm_dif = static_cast<int>(risk_size - dj + 1);
                    for (int i = 0; i< nm_dif; i++) {
                        //  start by incrementing the sum
                        r_sum += riskd0_list[i] * riskd1_list[i] / risk_list[i];
                        Recur_Second[group_ij][der_ijk][i] = r_sum;
                    }
                    //  now we need to progress through the remaining entries
                    for (int i_index = 1; i_index < dj; i_index ++) {
                        //  our increment in m
                        for (int j_index = 0; j_index < nm_dif; j_index ++) {
                            //  our increment in n
                            int recur_index = (i_index)*(nm_dif) + j_index;  //  the index of the value we are trying to fill
                            int risk_index = j_index + i_index;  //  the index of the risk value at this n
                            int t0 = recur_index - 1;  //  index for B(m, n-1)
                            int t1 = recur_index - nm_dif;  //  index for B(m-1, n-1)
                            //  the filled value is either an edge case, dB(m,n) = rn*dB(m-1, n-1) + drn*B(m-1, n-1), or the full case
                            if (j_index == 0) {
                                //  edge case
                                Recur_Second[group_ij][der_ijk][recur_index] = riskd0_list[risk_index] * riskd1_list[risk_index] / risk_list[risk_index] * Recur_Base[group_ij][t1] + riskd0_list[risk_index] * Recur_First[group_ij][der_jk][t1] + riskd1_list[risk_index] * Recur_First[group_ij][der_ij][t1] + risk_list[risk_index] * Recur_Second[group_ij][der_ijk][t1];
                            } else {
                                Recur_Second[group_ij][der_ijk][recur_index] = Recur_Second[group_ij][der_ijk][t0] + riskd0_list[risk_index] * riskd1_list[risk_index] / risk_list[risk_index] * Recur_Base[group_ij][t1] + riskd0_list[risk_index] * Recur_First[group_ij][der_jk][t1] + riskd1_list[risk_index] * Recur_First[group_ij][der_ij][t1] + risk_list[risk_index] * Recur_Second[group_ij][der_ijk][t1];
                            }
                        }
                    }
                }
            }
        }
    }
    return;
}

//' Fills out the expected information matrix for a case-control regression
//'
//' \code{Expected_Inform_Matrix_CaseCon} Called to update the expected information matrix
//'
//' @return Updates matrices in place: risk storage matrices
//' @noRd
//'
void Expected_Inform_Matrix_CaseCon(List& model_bool, const int& group_num, const IntegerMatrix& RiskFail, const vector<vector<int> >& RiskPairs, const int& totalnum, const int& ntime, const MatrixXd& R, const MatrixXd& Rd, const MatrixXd& Rdd, const MatrixXd& RdR, const MatrixXd& RddR, vector<double>& InMa, vector<vector<double> >& Recur_Base, vector<vector<vector<double> > >& Recur_First, vector<vector<vector<double> > >& Recur_Second, vector<double>& strata_odds, const int& nthreads, const IntegerVector& KeepConstant, vector<int>& strata_cond) {
    int reqrdnum = totalnum - sum(KeepConstant);
//    int kept_strata = group_num - reduce(strata_cond.begin(), strata_cond.end());
//    int total_val = reqrdnum + kept_strata;
    #ifdef _OPENMP
    #pragma omp declare reduction(vec_double_plus : vector<double> : \
        transform(omp_out.begin(), omp_out.end(), omp_in.begin(), omp_out.begin(), plus<double>())) \
        initializer(omp_priv = omp_orig)
    #endif
    double cond_thres = model_bool["cond_thres"];
    //  now we can calculate the expected information matrix
    #ifdef _OPENMP
    #pragma omp parallel for schedule(dynamic) num_threads(nthreads) \
        reduction(vec_double_plus:InMa) collapse(2)
    #endif
    for (int group_ij = 0; group_ij < group_num; group_ij++) {
        for (int der_ijk = 0; der_ijk < reqrdnum*(reqrdnum + 1)/2; der_ijk++) {
            int der_ij = 0;
            int der_jk = der_ijk;
            while (der_jk > der_ij) {
                der_ij++;
                der_jk -= der_ij;
            }
            //
            int recur_index = static_cast<int>(Recur_Base[group_ij].size() - 1);
            int dj = RiskFail(group_ij, 1)-RiskFail(group_ij, 0) + 1;
            if (recur_index > -1) {
                //  calculates the right-hand side terms
                if (dj <= cond_thres) {
                    double b_0 = Recur_Base[group_ij][recur_index];
                    double b_1 = Recur_First[group_ij][der_ij][recur_index];
                    double b_2 = Recur_First[group_ij][der_jk][recur_index];
                    double b_3 = Recur_Second[group_ij][der_ijk][recur_index];
                    //
                    InMa[der_ij*reqrdnum+der_jk] += b_3 / b_0 - b_1 / b_0 * b_2 / b_0;
                } else {
                    vector<int> InGroup = RiskPairs[group_ij];
                    double Ld3 = 0.0;
                    double Rs3 = 0.0;
                    double Rs3l = 0.0;
                    double Rs3r = 0.0;
                    for (vector<double>::size_type i = 0; i < InGroup.size() - 1; i = i+2) {
                        MatrixXd temp1 = MatrixXd::Zero(InGroup[i + 1]-InGroup[i] + 1, 1);
                        MatrixXd temp2 = MatrixXd::Zero(InGroup[i + 1]-InGroup[i] + 1, 1);
                        //
                        temp1 = RdR.block(InGroup[i] - 1, der_ij, InGroup[i + 1]-InGroup[i] + 1, 1).array() * RdR.block(InGroup[i] - 1, der_jk, InGroup[i + 1]-InGroup[i] + 1, 1).array();
                        temp1 = RddR.block(InGroup[i] - 1, der_ijk, InGroup[i + 1]-InGroup[i] + 1, 1).array() - temp1.array();  //  r''/r - r'/r x r'/r
                        temp2 = exp(strata_odds[group_ij]) * R.block(InGroup[i] - 1, 0, InGroup[i + 1]-InGroup[i] + 1, 1).array();  //  the odds ratio
                        temp2 = temp2.array() * (1 + temp2.array()).array().pow(-1).array();  //  the approximated probability of a case
                        temp1 = temp1.array() * temp2.array();
                        Ld3 += (temp1.array().isFinite()).select(temp1, 0).sum();
                        //
                        Rs3r += (Rd.block(InGroup[i] - 1, der_ij, InGroup[i + 1]-InGroup[i] + 1, 1).array() * Rd.block(InGroup[i] - 1, der_jk, InGroup[i + 1]-InGroup[i] + 1, 1).array() * (1.0 + exp(strata_odds[group_ij]) * R.block(InGroup[i] - 1, 0, InGroup[i + 1]-InGroup[i] + 1, 1).array()).pow(-2).array() ).sum();
                        Rs3l += (Rdd.block(InGroup[i] - 1, der_ijk, InGroup[i + 1]-InGroup[i] + 1, 1).array() * (1.0 + exp(strata_odds[group_ij]) * R.block(InGroup[i] - 1, 0, InGroup[i + 1]-InGroup[i] + 1, 1).array()).pow(-1).array()).sum();
                    }
                    Rs3 = exp(strata_odds[group_ij]) * (Rs3l - exp(strata_odds[group_ij])*Rs3r);
                    InMa[der_ij*reqrdnum+der_jk] -= Ld3 - Rs3;
                }
            }
        }
    }
    #ifdef _OPENMP
    #pragma omp parallel for schedule(dynamic) num_threads(nthreads)
    #endif
    for (int ijk = 0; ijk < reqrdnum*(reqrdnum + 1)/2; ijk++) {  //  fills second-derivative matrix
        int ij = 0;
        int jk = ijk;
        while (jk > ij) {
            ij++;
            jk -= ij;
        }
        InMa[jk*reqrdnum+ij] = InMa[ij*reqrdnum+jk];
    }
    return;
}
