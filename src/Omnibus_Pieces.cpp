#include <RcppEigen.h>
#ifdef _OPENMP
#include <omp.h>
#endif
#include "Omnibus_Pieces.h"
#include "Calc_Repeated.h"
#include "Subterms_Risk.h"
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


//' Utility function to refresh risk and subterm matrices for Cox Omnibus function
//'
//' \code{Cox_Refresh_R_TERM} Called to update matrices
//' @inheritParams CPP_template
//'
//' @return Updates matrices in place
//' @noRd
//'
// [[Rcpp::export]]
void Cox_Refresh_R_TERM(const int& totalnum, const int& reqrdnum, const int& term_tot, double& dint, double& dslp, double& dose_abs_max, double& abs_max, const MatrixXd& df0, MatrixXd& T0, MatrixXd& Td0, MatrixXd& Tdd0, MatrixXd& Te, MatrixXd& R, MatrixXd& Rd, MatrixXd& Rdd, MatrixXd& Dose, MatrixXd& nonDose, MatrixXd& TTerm, MatrixXd& nonDose_LIN, MatrixXd& nonDose_PLIN, MatrixXd& nonDose_LOGLIN, MatrixXd& RdR, MatrixXd& RddR, List& model_bool) {
    T0 = MatrixXd::Zero(df0.rows(), totalnum);  // preallocates matrix for Term column
    if (model_bool["basic"]) {
        //
        R = MatrixXd::Zero(df0.rows(), 1);  // preallocates matrix for Risks
        Rd = MatrixXd::Zero(df0.rows(), reqrdnum);  // preallocates matrix for Risk derivatives
        Rdd = MatrixXd::Zero(df0.rows(), reqrdnum*(reqrdnum + 1)/2);  // preallocates matrix for Risk second derivatives
        RdR = MatrixXd::Zero(df0.rows(), reqrdnum);  // preallocates matrix for Risk to derivative ratios
        TTerm = MatrixXd::Zero(df0.rows(), 1);  // matrix of term values
    } else if (model_bool["linear_err"]) {
        R = MatrixXd::Zero(df0.rows(), 1);  // preallocates matrix for Risks
        Rd = MatrixXd::Zero(df0.rows(), reqrdnum);  // preallocates matrix for Risk derivatives
        Rdd = MatrixXd::Zero(df0.rows(), reqrdnum*(reqrdnum + 1)/2);  // preallocates matrix for Risk second derivatives
        //
        nonDose_PLIN = MatrixXd::Constant(df0.rows(), 1, 1.0);  // matrix of Loglinear subterm values
        nonDose_LOGLIN = MatrixXd::Constant(df0.rows(), 1, 1.0);  // matrix of Product linear subterm values
        TTerm = MatrixXd::Zero(df0.rows(), 1);  // matrix of term values
        RdR = MatrixXd::Zero(df0.rows(), reqrdnum);  // preallocates matrix for Risk to derivative ratios
        RddR = MatrixXd::Zero(df0.rows(), reqrdnum*(reqrdnum + 1)/2);  // preallocates matrix for Risk to second derivative ratios
    } else if (model_bool["single"]) {
        //
        Te = MatrixXd::Zero(df0.rows(), 1);  // preallocates matrix for column terms used for temporary storage
        R = MatrixXd::Zero(df0.rows(), 1);  // preallocates matrix for Risks
        //
        Dose = MatrixXd::Constant(df0.rows(), term_tot, 0.0);  // matrix of the total dose term values
        nonDose = MatrixXd::Constant(df0.rows(), term_tot, 1.0);  // matrix of the total non-dose term values
        nonDose_LIN = MatrixXd::Constant(df0.rows(), term_tot, 0.0);  // matrix of Linear subterm values
        nonDose_PLIN = MatrixXd::Constant(df0.rows(), term_tot, 1.0);  // matrix of Loglinear subterm values
        nonDose_LOGLIN = MatrixXd::Constant(df0.rows(), term_tot, 1.0);  // matrix of Product linear subterm values
        TTerm = MatrixXd::Zero(Dose.rows(), Dose.cols());  // matrix of term values
    } else if (model_bool["gradient"]){
        Td0 = MatrixXd::Zero(df0.rows(), reqrdnum);  // preallocates matrix for Term derivative columns
        //
        Te = MatrixXd::Zero(df0.rows(), 1);  // preallocates matrix for column terms used for temporary storage
        R = MatrixXd::Zero(df0.rows(), 1);  // preallocates matrix for Risks
        Rd = MatrixXd::Zero(df0.rows(), reqrdnum);  // preallocates matrix for Risk derivatives
        //
        nonDose = MatrixXd::Constant(df0.rows(), term_tot, 1.0);  // matrix of the total non-dose term values
        nonDose_LIN = MatrixXd::Constant(df0.rows(), term_tot, 0.0);  // matrix of Linear subterm values
        nonDose_PLIN = MatrixXd::Constant(df0.rows(), term_tot, 1.0);  // matrix of Loglinear subterm values
        nonDose_LOGLIN = MatrixXd::Constant(df0.rows(), term_tot, 1.0);  // matrix of Product linear subterm values
        TTerm = MatrixXd::Zero(df0.rows(), term_tot);  // matrix of term values
        RdR = MatrixXd::Zero(df0.rows(), reqrdnum);  // preallocates matrix for Risk to derivative ratios
    } else {
        Td0 = MatrixXd::Zero(df0.rows(), reqrdnum);  // preallocates matrix for Term derivative columns
        Tdd0 = MatrixXd::Zero(df0.rows(), reqrdnum*(reqrdnum + 1)/2);  // preallocates matrix for Term second derivative columns
        //
        Te = MatrixXd::Zero(df0.rows(), 1);  // preallocates matrix for column terms used for temporary storage
        R = MatrixXd::Zero(df0.rows(), 1);  // preallocates matrix for Risks
        Rd = MatrixXd::Zero(df0.rows(), reqrdnum);  // preallocates matrix for Risk derivatives
        Rdd = MatrixXd::Zero(df0.rows(), reqrdnum*(reqrdnum + 1)/2);  // preallocates matrix for Risk second derivatives
        //
        Dose = MatrixXd::Constant(df0.rows(), term_tot, 0.0);  // matrix of the total dose term values
        nonDose = MatrixXd::Constant(df0.rows(), term_tot, 1.0);  // matrix of the total non-dose term values
        nonDose_LIN = MatrixXd::Constant(df0.rows(), term_tot, 0.0);  // matrix of Linear subterm values
        nonDose_PLIN = MatrixXd::Constant(df0.rows(), term_tot, 1.0);  // matrix of Loglinear subterm values
        nonDose_LOGLIN = MatrixXd::Constant(df0.rows(), term_tot, 1.0);  // matrix of Product linear subterm values
        TTerm = MatrixXd::Zero(Dose.rows(), Dose.cols());  // matrix of term values
        dint = dose_abs_max;  // the amount of change used to calculate derivatives in threshold paramters
        dslp = abs_max;
        RdR = MatrixXd::Zero(df0.rows(), reqrdnum);  // preallocates matrix for Risk to derivative ratios
        RddR = MatrixXd::Zero(df0.rows(), reqrdnum*(reqrdnum + 1)/2);  // preallocates matrix for Risk to second derivative ratios
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
// [[Rcpp::export]]
void Cox_Refresh_R_SIDES(const int& reqrdnum, const int& ntime, MatrixXd& Rls1, MatrixXd& Rls2, MatrixXd& Rls3, MatrixXd& Lls1, MatrixXd& Lls2, MatrixXd& Lls3, NumericVector& Strata_vals, List& model_bool) {
    if (model_bool["strata"]) {
        Rls1 = MatrixXd::Zero(ntime, Strata_vals.size());  // precomputes a series of sums used frequently in the log-liklihood calculations
        Lls1 = MatrixXd::Zero(ntime, Strata_vals.size());  // the log-likelihood calculation has a Right and Left sum used
        if (!model_bool["single"]) {
            Rls2 = MatrixXd::Zero(ntime, reqrdnum*Strata_vals.size());  // many are repeated due to the same risk groups and derivatives being used at mulitple points
            Lls2 = MatrixXd::Zero(ntime, reqrdnum*Strata_vals.size());
            if (!model_bool["gradient"]){
                Rls3 = MatrixXd::Zero(ntime, reqrdnum*(reqrdnum + 1)/2*Strata_vals.size());  // sum and its derivatives are precomputed
                Lls3 = MatrixXd::Zero(ntime, reqrdnum*(reqrdnum + 1)/2*Strata_vals.size());
            }
        }
    } else {
        Rls1 = MatrixXd::Zero(ntime, 1);  // precomputes a series of sums used frequently in the log-liklihood calculations
        Lls1 = MatrixXd::Zero(ntime, 1);  // the log-likelihood calculation has a Right and Left sum used
        if (!model_bool["single"]) {
            Rls2 = MatrixXd::Zero(ntime, reqrdnum);  // many are repeated due to the same risk groups and derivatives being used at mulitple points
            Lls2 = MatrixXd::Zero(ntime, reqrdnum);
            if (!model_bool["gradient"]){
                Rls3 = MatrixXd::Zero(ntime, reqrdnum*(reqrdnum + 1)/2);  // sum and its derivatives are precomputed
                Lls3 = MatrixXd::Zero(ntime, reqrdnum*(reqrdnum + 1)/2);
            }
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
// [[Rcpp::export]]
void Cox_Term_Risk_Calc(string modelform, const StringVector& tform, const IntegerVector& term_n, const int& totalnum, const int& fir, const IntegerVector& dfc, int term_tot, MatrixXd& T0, MatrixXd& Td0, MatrixXd& Tdd0, MatrixXd& Te, MatrixXd& R, MatrixXd& Rd, MatrixXd& Rdd, MatrixXd& Dose, MatrixXd& nonDose, VectorXd beta_0, const  MatrixXd& df0, const double& dint, const double& dslp, MatrixXd& TTerm, MatrixXd& nonDose_LIN, MatrixXd& nonDose_PLIN, MatrixXd& nonDose_LOGLIN, MatrixXd& RdR, MatrixXd& RddR, const int& nthreads, bool debugging, const IntegerVector& KeepConstant, int verbose, List& model_bool, const double gmix_theta, const IntegerVector& gmix_term) {
    // time_point<system_clock> comp_point, end_point;
    // end_point = system_clock::now();
    // auto ending = time_point_cast<microseconds>(end_point).time_since_epoch().count();  // the time duration is tracked
    // comp_point = system_clock::now();
    // auto comp = time_point_cast<microseconds>(end_point).time_since_epoch().count();  // the time duration is tracked
    int reqrdnum = totalnum - sum(KeepConstant);
    if (model_bool["basic"]) {
        // Calculates the subterm and term values
        Make_subterms_Basic(totalnum, dfc, T0, beta_0, df0, nthreads, debugging);
        // ---------------------------------------------------------
        // Prints off a series of calculations to check at what point values are changing
        // ---------------------------------------------------------
        //
//        if (verbose >= 4) {
//            Rcout << "C++ Note: values checked ";
//            for (int ijk = 0; ijk < totalnum; ijk++) {
//                Rcout << beta_0[ijk] << " ";
//            }
//            Rcout << " " << endl;
//            Rcout << "C++ Note: sums checked ";
//            for (int ijk = 0; ijk < totalnum; ijk++) {
//                Rcout << T0.col(ijk).sum() << " ";
//            }
//            Rcout << " " << endl;
//        }
        //
        // Calculates the risk for each row
        Make_Risks_Basic(totalnum, T0, R, Rd, Rdd, RdR, nthreads, debugging, df0, dfc, KeepConstant);
        //
//        Rcout << R.block(3300900, 0, 1, 1 ) << " " << Rd.block(3300900, 7, 1, 1 ).sum() << " " << Rd.block(3300900, 7, 1, 1 ).sum() << " " << Rdd.block(3300900, 35, 1, 1 ).sum() << endl;
        // Removes infinite values
        RdR = (RdR.array().isFinite()).select(RdR, 0);
        //
//        Rcout << "C++ Note: values checked ";
//        for (int ijk = 0; ijk < totalnum; ijk++) {
//            Rcout << beta_0[ijk] << " ";
//        }
//        Rcout << " " << endl;
//        Rcout << "C++ Note: RdR checked ";
//        for (int ijk = 0; ijk < reqrdnum; ijk++) {
//            Rcout << RdR.col(ijk).sum() << " ";
//        }
//        Rcout << " " << endl;
        //
        if (R.minCoeff() <= 0) {
            if (verbose >= 4) {
                Rcout << "C++ Warning: risk mininum " << R.minCoeff() << " " << endl;
            }
        } else if (verbose >= 4) {
            Rcout << "C++ Note: risk checked ";
            for (int ijk = 0; ijk < 1; ijk++) {
                Rcout << R.col(0).sum() << " ";
            }
            Rcout << " " << endl;
            Rcout << "C++ Note: risk1 checked ";
            for (int ijk = 0; ijk < reqrdnum; ijk++) {
                Rcout << Rd.col(ijk).sum() << " ";
            }
            Rcout << " " << endl;
            Rcout << "C++ Note: rdr checked ";
            for (int ijk = 0; ijk < reqrdnum; ijk++) {
                Rcout << RdR.col(ijk).sum() << " ";
            }
            Rcout << " " << endl;
            Rcout << "C++ Note: risk2 checked ";
            for (int ijk = 0; ijk < reqrdnum; ijk++) {
                Rcout << Rdd.col(ijk*(ijk + 1)/2+ijk).sum() << " ";
            }
            Rcout << " " << endl;
            //
            Rcout << "C++ Note: ALL risk2 checked ";
            for (int ijk = 0; ijk < reqrdnum*(reqrdnum + 1)/2; ijk++) {
                Rcout << Rdd.col(ijk).sum() << " ";
            }
            Rcout << " " << endl;
            //
        }
    } else if (model_bool["linear_err"]) {
        Make_subterms_Linear_ERR(totalnum, tform, dfc, nonDose_PLIN, nonDose_LOGLIN, beta_0, df0, nthreads, debugging, KeepConstant);
        //
        Make_Risks_Linear_ERR(tform, dfc, df0, totalnum, R, Rd, Rdd, nonDose_PLIN, nonDose_LOGLIN, RdR, RddR, nthreads, debugging, KeepConstant);
        RdR = (RdR.array().isFinite()).select(RdR, 0);
        RddR = (RddR.array().isFinite()).select(RddR, 0);
        //
        TTerm = R.col(0).array();
//        if (R.minCoeff() <= 0) {
//            if (verbose >= 4) {
//                Rcout << "C++ Warning: risk mininum " << R.minCoeff() << " " << endl;
//            }
//        } else if (verbose >= 4) {
//            Rcout << "C++ Note: risk checked ";
//            for (int ijk = 0; ijk < 1; ijk++) {
//                Rcout << R.col(0).sum() << " ";
//            }
//            Rcout << " " << endl;
//            Rcout << "C++ Note: risk1 checked ";
//            for (int ijk = 0; ijk < reqrdnum; ijk++) {
//                Rcout << Rd.col(ijk).sum() << " ";
//            }
//            Rcout << " " << endl;
//            Rcout << "C++ Note: risk2 checked ";
//            for (int ijk = 0; ijk < reqrdnum; ijk++) {
//                Rcout << Rdd.col(ijk*(ijk + 1)/2+ijk).sum() << " ";
//            }
//            Rcout << " " << endl;
//            //
//            Rcout << "C++ Note: ALL risk2 checked ";
//            for (int ijk = 0; ijk < reqrdnum*(reqrdnum + 1)/2; ijk++) {
//                Rcout << Rdd.col(ijk).sum() << " ";
//            }
//            Rcout << " " << endl;
//            //
//        }
    } else if (model_bool["single"]) {
        // Calculates the subterm and term values
        Make_subterms_Single(totalnum, term_n, tform, dfc, fir, T0, Dose, nonDose, TTerm, nonDose_LIN, nonDose_PLIN, nonDose_LOGLIN, beta_0, df0, nthreads, debugging, KeepConstant);
        // ---------------------------------------------------------
        // Prints off a series of calculations to check at what point values are changing
        // ---------------------------------------------------------
        //
        // if (verbose >= 4) {
        //     Rcout << "C++ Note: values checked ";
        //     for (int ijk = 0; ijk < totalnum; ijk++) {
        //         Rcout << beta_0[ijk] << " ";
        //     }
        //     Rcout << " " << endl;
        //     Rcout << "C++ Note: sums checked ";
        //     for (int ijk = 0; ijk < totalnum; ijk++) {
        //         Rcout << T0.col(ijk).sum() << " ";
        //     }
        //     Rcout << " " << endl;
        //     Rcout << "C++ Note: dose checked ";
        //     for (int ijk = 0; ijk < term_tot; ijk++) {
        //         Rcout << Dose.col(ijk).array().sum() << " ";
        //     }
        //     Rcout << " " << endl;
        //     Rcout << "C++ Note: non-dose checked ";
        //     for (int ijk = 0; ijk < term_tot; ijk++) {
        //         Rcout << nonDose.col(ijk).array().sum() << " ";
        //     }
        //     Rcout << " " << endl;
        //     Rcout << "C++ Note: LIN_non-dose checked ";
        //     for (int ijk = 0; ijk < term_tot; ijk++) {
        //         Rcout << nonDose_LIN.col(ijk).array().sum() << " ";
        //     }
        //     Rcout << " " << endl;
        //     Rcout << "C++ Note: PLIN_non-dose checked ";
        //     for (int ijk = 0; ijk < term_tot; ijk++) {
        //         Rcout << nonDose_PLIN.col(ijk).array().sum() << " ";
        //     }
        //     Rcout << " " << endl;
        //     Rcout << "C++ Note: LOGLIN_non-dose checked ";
        //     for (int ijk = 0; ijk < term_tot; ijk++) {
        //         Rcout << nonDose_LOGLIN.col(ijk).array().sum() << " ";
        //     }
        //     Rcout << " " << endl;
        // }
        //
        //
        // Calculates the risk for each row
        Make_Risks_Single(modelform, tform, term_n, totalnum, fir, T0, Te, R, Dose, nonDose, TTerm, nonDose_LIN, nonDose_PLIN, nonDose_LOGLIN, nthreads, debugging, KeepConstant, gmix_theta, gmix_term);
        //
        // Removes infinite values
        //
        //
        // if (R.minCoeff() <= 0) {
        //     if (verbose >= 4) {
        //         Rcout << "C++ Warning: risk mininum " << R.minCoeff() << " " << endl;
        //     }
        // } else if (verbose >= 4) {
        //     Rcout << "C++ Note: risk checked ";
        //     for (int ijk = 0; ijk < 1; ijk++) {
        //         Rcout << R.col(0).sum() << " ";
        //     }
        //     Rcout << " " << endl;
        // }
    } else if (model_bool["gradient"]) {
        //
        Make_subterms_Gradient(totalnum, term_n, tform, dfc, fir, T0, Td0, nonDose, TTerm, nonDose_LIN, nonDose_PLIN, nonDose_LOGLIN, beta_0, df0, nthreads, debugging, KeepConstant);
        // if (verbose >= 4) {
        //    Rcout << "C++ Note: values checked ";
        //    for (int ijk = 0; ijk < totalnum; ijk++) {
        //        Rcout << beta_0[ijk] << " ";
        //    }
        //    Rcout << " " << endl;
        //    Rcout << "C++ Note: sums checked ";
        //    for (int ijk = 0; ijk < totalnum; ijk++) {
        //        Rcout << T0.col(ijk).sum() << " ";
        //    }
        //    Rcout << " " << endl;
        //    Rcout << "C++ Note: derivs checked ";
        //    for (int ijk = 0; ijk < reqrdnum; ijk++) {
        //        Rcout << Td0.col(ijk).sum() << " ";
        //    }
        //    Rcout << " " << endl;
        //    //
        //    Rcout << "C++ Note: non-dose checked ";
        //    for (int ijk = 0; ijk < term_tot; ijk++) {
        //        Rcout << nonDose.col(ijk).array().sum() << " ";
        //    }
        //    Rcout << " " << endl;
        //    Rcout << "C++ Note: LIN_non-dose checked ";
        //    for (int ijk = 0; ijk < term_tot; ijk++) {
        //        Rcout << nonDose_LIN.col(ijk).array().sum() << " ";
        //    }
        //    Rcout << " " << endl;
        //    Rcout << "C++ Note: PLIN_non-dose checked ";
        //    for (int ijk = 0; ijk < term_tot; ijk++) {
        //        Rcout << nonDose_PLIN.col(ijk).array().sum() << " ";
        //    }
        //    Rcout << " " << endl;
        //    Rcout << "C++ Note: LOGLIN_non-dose checked ";
        //    for (int ijk = 0; ijk < term_tot; ijk++) {
        //        Rcout << nonDose_LOGLIN.col(ijk).array().sum() << " ";
        //    }
        //    Rcout << " " << endl;
        // }
        Make_Risks_Gradient(modelform, tform, term_n, totalnum, fir, T0, Td0, Te, R, Rd, nonDose, TTerm, nonDose_LIN, nonDose_PLIN, nonDose_LOGLIN, RdR, nthreads, debugging, KeepConstant);
        // if (R.minCoeff() <= 0) {
        //    if (verbose >= 4) {
        //        Rcout << "C++ Warning: risk mininum " << R.minCoeff() << " " << endl;
        //    }
        // } else if (verbose >= 4) {
        //    Rcout << "C++ Note: risk checked ";
        //    for (int ijk = 0; ijk < 1; ijk++) {
        //        Rcout << R.col(0).sum() << " ";
        //    }
        //    Rcout << " " << endl;
        //    Rcout << "C++ Note: risk1 checked ";
        //    for (int ijk = 0; ijk < reqrdnum; ijk++) {
        //        Rcout << Rd.col(ijk).sum() << " ";
        //    }
        //    Rcout << " " << endl;
        // }
        //
        RdR = (RdR.array().isFinite()).select(RdR, 0);
        //
    } else {
        //
        // Calculates the subterm and term values
        //
        // comp_point = system_clock::now();
        // comp = time_point_cast<microseconds>(comp_point).time_since_epoch().count();
        Make_subterms(totalnum, term_n, tform, dfc, fir, T0, Td0, Tdd0, Dose, nonDose, TTerm, nonDose_LIN, nonDose_PLIN, nonDose_LOGLIN, beta_0, df0, dint, dslp, nthreads, debugging, KeepConstant);
        // end_point = system_clock::now();
        // ending = time_point_cast<microseconds>(end_point).time_since_epoch().count();
//        Rcout << "C++ Note: subterm " << (ending-comp) * 1e-6  <<endl;
        // ---------------------------------------------------------
        // Prints off a series of calculations to check at what point values are changing
        // ---------------------------------------------------------
        //
        //
//        if (verbose >= 4) {
//            Rcout << "C++ Note: values checked ";
//            for (int ijk = 0; ijk < totalnum; ijk++) {
//                Rcout << beta_0[ijk] << " ";
//            }
//            Rcout << " " << endl;
//            Rcout << "C++ Note: sums checked ";
//            for (int ijk = 0; ijk < totalnum; ijk++) {
//                Rcout << T0.col(ijk).sum() << " ";
//            }
//            Rcout << " " << endl;
//            Rcout << "C++ Note: derivs checked ";
//            for (int ijk = 0; ijk < reqrdnum; ijk++) {
//                Rcout << Td0.col(ijk).sum() << " ";
//            }
//            Rcout << " " << endl;
//            Rcout << "C++ Note: second derivs checked ";
//            for (int ijk = 0; ijk < reqrdnum; ijk++) {
//                Rcout << Tdd0.col(ijk*(ijk + 1)/2+ijk).sum() << " ";
//            }
//            Rcout << " " << endl;
//            //
//            Rcout << "C++ Note: ALL second derivs checked ";
//            for (int ijk = 0; ijk < reqrdnum*(reqrdnum + 1)/2; ijk++) {
//                Rcout << Tdd0.col(ijk).sum() << " ";
//            }
//            Rcout << " " << endl;
//            //
//            Rcout << "C++ Note: dose checked ";
//            for (int ijk = 0; ijk < term_tot; ijk++) {
//                Rcout << Dose.col(ijk).array().sum() << " ";
//            }
//            Rcout << " " << endl;
//            Rcout << "C++ Note: non-dose checked ";
//            for (int ijk = 0; ijk < term_tot; ijk++) {
//                Rcout << nonDose.col(ijk).array().sum() << " ";
//            }
//            Rcout << " " << endl;
//            Rcout << "C++ Note: LIN_non-dose checked ";
//            for (int ijk = 0; ijk < term_tot; ijk++) {
//                Rcout << nonDose_LIN.col(ijk).array().sum() << " ";
//            }
//            Rcout << " " << endl;
//            Rcout << "C++ Note: PLIN_non-dose checked ";
//            for (int ijk = 0; ijk < term_tot; ijk++) {
//                Rcout << nonDose_PLIN.col(ijk).array().sum() << " ";
//            }
//            Rcout << " " << endl;
//            Rcout << "C++ Note: LOGLIN_non-dose checked ";
//            for (int ijk = 0; ijk < term_tot; ijk++) {
//                Rcout << nonDose_LOGLIN.col(ijk).array().sum() << " ";
//            }
//            Rcout << " " << endl;
//        }
        //
        //
        // Calculates the risk for each row
        // comp_point = system_clock::now();
        // comp = time_point_cast<microseconds>(comp_point).time_since_epoch().count();
        Make_Risks(modelform, tform, term_n, totalnum, fir, T0, Td0, Tdd0, Te, R, Rd, Rdd, Dose, nonDose, TTerm, nonDose_LIN, nonDose_PLIN, nonDose_LOGLIN, RdR, RddR, nthreads, debugging, KeepConstant, gmix_theta, gmix_term);
        // end_point = system_clock::now();
        // ending = time_point_cast<microseconds>(end_point).time_since_epoch().count();
//        Rcout << "C++ Note: risk " << (ending-comp) * 1e-6  <<endl;
        //
        // Removes infinite values
        RdR = (RdR.array().isFinite()).select(RdR, 0);
        RddR = (RddR.array().isFinite()).select(RddR, 0);
//        Rcout << R.block(3300900, 0, 1, 1 ) << " " << Rd.block(3300900, 7, 1, 1 ).sum() << " " << Rd.block(3300900, 7, 1, 1 ).sum() << " " << Rdd.block(3300900, 35, 1, 1 ).sum() << endl;
//        Rcout << RdR.block(3300900, 7, 1, 1 ).sum() << " " << RdR.block(3300900, 7, 1, 1 ).sum() << " " << RddR.block(3300900, 35, 1, 1 ).sum() << " " << R.block(3300900, 0, 1, 1 ) << " " << Rd.block(3300900, 7, 1, 1 ).sum() << " " << Rd.block(3300900, 7, 1, 1 ).sum() << " " << Rdd.block(3300900, 35, 1, 1 ).sum() << endl;
        // 3,764,306
//        int start_point = 800000;
//        Rcout << start_point;
//        int skip = 10000;
//        for (int ijk=start_point; ijk+skip < 1000000; ijk = ijk + skip) {
//            Rcout << " " << RdR.block(ijk, 1, skip, 1).sum();
//        }
//        Rcout << " " << endl;
        //
//        Rcout << "C++ Note: values checked ";
//        for (int ijk = 0; ijk < totalnum; ijk++) {
//            Rcout << beta_0[ijk] << " ";
//        }
//        Rcout << " " << endl;
//        Rcout << "C++ Note: RdR checked ";
//        for (int ijk = 0; ijk < reqrdnum; ijk++) {
//            Rcout << RdR.col(ijk).sum() << " ";
//        }
//        Rcout << " " << endl;
//        Rcout << "C++ Note: RddR checked ";
//        for (int ijk = 0; ijk < reqrdnum*(reqrdnum + 1)/2; ijk++) {
//            Rcout << RddR.col(ijk).sum() << " ";
//        }
//        Rcout << " " << endl;
        //
        if (R.minCoeff() <= 0) {
            if (verbose >= 4) {
                Rcout << "C++ Warning: risk mininum " << R.minCoeff() << " " << endl;
            }
        } else if (verbose >= 4) {
            Rcout << "C++ Note: risk checked ";
            for (int ijk = 0; ijk < 1; ijk++) {
                Rcout << R.col(0).sum() << " ";
            }
            Rcout << " " << endl;
            Rcout << "C++ Note: risk1 checked ";
            for (int ijk = 0; ijk < reqrdnum; ijk++) {
                Rcout << Rd.col(ijk).sum() << " ";
            }
            Rcout << " " << endl;
            Rcout << "C++ Note: rdr checked ";
            for (int ijk = 0; ijk < reqrdnum; ijk++) {
                Rcout << RdR.col(ijk).sum() << " ";
            }
            Rcout << " " << endl;
            Rcout << "C++ Note: risk2 checked ";
            for (int ijk = 0; ijk < reqrdnum; ijk++) {
                Rcout << Rdd.col(ijk*(ijk + 1)/2+ijk).sum() << " ";
            }
            Rcout << " " << endl;
            //
            Rcout << "C++ Note: ALL risk2 checked ";
            for (int ijk = 0; ijk < reqrdnum*(reqrdnum + 1)/2; ijk++) {
                Rcout << Rdd.col(ijk).sum() << " ";
            }
            Rcout << " " << endl;
            //
        }
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
// [[Rcpp::export]]
void Cox_Side_LL_Calc(const int& reqrdnum, const int& ntime, const StringVector& tform, const IntegerMatrix& RiskFail, const StringMatrix&  RiskGroup_Strata, const vector<string>& RiskGroup, const int& totalnum, const int& fir, MatrixXd& R, MatrixXd& Rd, MatrixXd& Rdd, MatrixXd& Rls1, MatrixXd& Rls2, MatrixXd& Rls3, MatrixXd& Lls1, MatrixXd& Lls2, MatrixXd& Lls3, const VectorXd& cens_weight, NumericVector& Strata_vals, VectorXd beta_0, MatrixXd& RdR, MatrixXd& RddR, vector<double>& Ll, vector<double>& Lld, vector<double>& Lldd, const int& nthreads, bool debugging, const IntegerVector& KeepConstant, string ties_method, int verbose, List& model_bool, int iter_stop) {
    // Calculates the side sum terms used
    //
    // time_point<system_clock> comp_point, end_point;
    // end_point = system_clock::now();
    // auto ending = time_point_cast<microseconds>(end_point).time_since_epoch().count();  // the time duration is tracked
    // comp_point = system_clock::now();
    // auto comp = time_point_cast<microseconds>(end_point).time_since_epoch().count();  // the time duration is tracked
    // //
    // comp_point = system_clock::now();
    // comp = time_point_cast<microseconds>(comp_point).time_since_epoch().count();
    if (model_bool["strata"]) {
        if (model_bool["cr"]) {
            // strata_CR or strata_CR_single
            if (model_bool["single"]) {
                Calculate_Sides_Strata_Single_CR(RiskFail, RiskGroup_Strata, totalnum, ntime, R, Rls1, Lls1, cens_weight, nthreads, debugging, Strata_vals, KeepConstant);  // strata_CR_single
            } else {
                if (model_bool["gradient"]){
                    Calculate_Sides_Strata_CR_Gradient(RiskFail, RiskGroup_Strata, totalnum, ntime, R, Rd, Rls1, Rls2, Lls1, Lls2, cens_weight, nthreads, debugging, Strata_vals, KeepConstant);
                } else {
                    Calculate_Sides_Strata_CR(RiskFail, RiskGroup_Strata, totalnum, ntime, R, Rd, Rdd, Rls1, Rls2, Rls3, Lls1, Lls2, Lls3, cens_weight, nthreads, debugging, Strata_vals, KeepConstant);  // strata_cr
                }
            }
        } else if (model_bool["single"]) {
            Calculate_Sides_Strata_Single(RiskFail, RiskGroup_Strata, totalnum, ntime, R, Rls1, Lls1, nthreads, debugging, Strata_vals, KeepConstant);  // strata_single
        } else if (model_bool["gradient"]) {
            Calculate_Sides_Strata_Gradient(RiskFail, RiskGroup_Strata, totalnum, ntime, R, Rd, Rls1, Rls2, Lls1, Lls2, nthreads, debugging, Strata_vals, KeepConstant);
        } else {
            Calculate_Sides_Strata(RiskFail, RiskGroup_Strata, totalnum, ntime, R, Rd, Rdd, Rls1, Rls2, Rls3, Lls1, Lls2, Lls3, nthreads, debugging, Strata_vals, KeepConstant);
        }
    } else if (model_bool["cr"]) {
        if (model_bool["single"]) {
            Calculate_Sides_CR_SINGLE(RiskFail, RiskGroup, totalnum, ntime, R, Rls1, Lls1, cens_weight, nthreads, debugging, KeepConstant);
        } else if (model_bool["gradient"]) {
            Calculate_Sides_CR_Gradient(RiskFail, RiskGroup, totalnum, ntime, R, Rd, Rls1, Rls2, Lls1, Lls2, cens_weight, nthreads, debugging, KeepConstant);
        } else {
            Calculate_Sides_CR(RiskFail, RiskGroup, totalnum, ntime, R, Rd, Rdd, Rls1, Rls2, Rls3, Lls1, Lls2, Lls3, cens_weight, nthreads, debugging, KeepConstant);
        }
    } else if (model_bool["single"]) {
        Calculate_Sides_Single(RiskFail, RiskGroup, totalnum, ntime, R, Rls1, Lls1, nthreads, debugging);
    } else if (model_bool["gradient"]) {
        Calculate_Sides_Gradient(RiskFail, RiskGroup, totalnum, ntime, R, Rd, Rls1, Rls2, Lls1, Lls2, nthreads, debugging, KeepConstant);
    } else {
        Calculate_Sides(RiskFail, RiskGroup, totalnum, ntime, R, Rd, Rdd, Rls1, Rls2, Rls3, Lls1, Lls2, Lls3, nthreads, debugging, KeepConstant);
    }
//    end_point = system_clock::now();
//    ending = time_point_cast<microseconds>(end_point).time_since_epoch().count();
//    Rcout << "C++ Note: sides " << (ending-comp) * 1e-6  <<endl;
    //
//    if (model_bool["strata"]) {
//        if (verbose >= 4) {
//            Rcout << "C++ Note: riskr checked ";
//            for (int ijk = 0; ijk < Strata_vals.size(); ijk++) {
//                Rcout << Rls1.col(ijk).sum() << " ";
//            }
//            Rcout << " " << endl;
//            //
//            Rcout << "C++ Note: riskl checked ";
//            for (int ijk = 0; ijk < Strata_vals.size(); ijk++) {
//                Rcout << Lls1.col(ijk).sum() << " ";
//            }
//            Rcout << " " << endl;
//        }
//    } else {
//        if (verbose >= 4) {
//            Rcout << "C++ Note: riskr checked ";
//            for (int ijk = 0; ijk < 1; ijk++) {
//                Rcout << Rls1.col(0).sum() << " ";
//            }
//            Rcout << " " << endl;
//            if (!model_bool["single"]) {
//                Rcout << "C++ Note: risk1r checked ";
//                for (int ijk = 0; ijk < reqrdnum; ijk++) {
//                    Rcout << Rls2.col(ijk).sum() << " ";
//                }
//                Rcout << " " << endl;
//                if (!model_bool["gradient"]){
//                    Rcout << "C++ Note: risk2r checked ";
//                    for (int ijk = 0; ijk < reqrdnum; ijk++) {
//                        Rcout << Rls3.col(ijk*(ijk + 1)/2+ijk).sum() << " ";
//                    }
//                    Rcout << " " << endl;
//                    //
//                    Rcout << "C++ Note: ALL risk2r checked ";
//                    for (int ijk = 0; ijk < reqrdnum*(reqrdnum + 1)/2; ijk++) {
//                        Rcout << Rls3.col(ijk).sum() << " ";
//                    }
//                    Rcout << " " << endl;
//                }
//                //
//            }
//            //
//            Rcout << "C++ Note: riskl checked ";
//            for (int ijk = 0; ijk < 1; ijk++) {
//                Rcout << Lls1.col(0).sum() << " ";
//            }
//            Rcout << " " << endl;
//            if (!model_bool["single"]) {
//                Rcout << "C++ Note: risk1l checked ";
//                for (int ijk = 0; ijk < reqrdnum; ijk++) {
//                    Rcout << Lls2.col(ijk).sum() << " ";
//                }
//                Rcout << " " << endl;
//                if (!model_bool["gradient"]){
//                    Rcout << "C++ Note: risk2l checked ";
//                    for (int ijk = 0; ijk < reqrdnum; ijk++) {
//                        Rcout << Lls3.col(ijk*(ijk + 1)/2+ijk).sum() << " ";
//                    }
//                    Rcout << " " << endl;
//                    //
//                    Rcout << "C++ Note: ALL risk2l checked ";
//                    for (int ijk = 0; ijk < reqrdnum*(reqrdnum + 1)/2; ijk++) {
//                        Rcout << Lls3.col(ijk).sum() << " ";
//                    }
//                    Rcout << " " << endl;
//                }
//                //
//            }
//        }
//    }
    // Calculates log-likelihood
    fill(Ll.begin(), Ll.end(), 0.0);
    if (!model_bool["single"]) {
        fill(Lld.begin(), Lld.end(), 0.0);
        if (!model_bool["gradient"]){
            fill(Lldd.begin(), Lldd.end(), 0.0);
        }
    }
    // comp_point = system_clock::now();
    // comp = time_point_cast<microseconds>(comp_point).time_since_epoch().count();
    if (model_bool["strata"]) {
        if (model_bool["single"]) {
            Calc_LogLik_Strata_SINGLE(nthreads, RiskFail, RiskGroup_Strata, totalnum, ntime, R, Rls1, Lls1, Ll, debugging, ties_method, Strata_vals, KeepConstant);  // strata_single
        } else if (model_bool["basic"]) {
            Calc_LogLik_Strata_BASIC(nthreads, RiskFail, RiskGroup_Strata, totalnum, ntime, R, Rd, Rdd, RdR, Rls1, Rls2, Rls3, Lls1, Lls2, Lls3, Ll, Lld, Lldd, debugging, ties_method, Strata_vals, KeepConstant);
        } else if (model_bool["linear_err"]) {
            Calc_LogLik_Strata_Linear_ERR(tform, nthreads, RiskFail, RiskGroup_Strata, totalnum, ntime, R, Rd, Rdd, RdR, RddR, Rls1, Rls2, Rls3, Lls1, Lls2, Lls3, Ll, Lld, Lldd, debugging, ties_method, Strata_vals, KeepConstant);
        } else if (model_bool["gradient"]) {
            Calc_LogLik_Strata_Gradient(nthreads, RiskFail, RiskGroup_Strata, totalnum, ntime, R, Rd, RdR, Rls1, Rls2, Lls1, Lls2, Ll, Lld, debugging, ties_method, Strata_vals, KeepConstant);
        } else {
            Calc_LogLik_Strata(nthreads, RiskFail, RiskGroup_Strata, totalnum, ntime, R, Rd, Rdd, RdR, RddR, Rls1, Rls2, Rls3, Lls1, Lls2, Lls3, Ll, Lld, Lldd, debugging, ties_method, Strata_vals, KeepConstant);
        }
    } else {
        if (model_bool["single"]) {
            Calc_LogLik_Single(nthreads, RiskFail, RiskGroup, totalnum, ntime, R, Rls1, Lls1, Ll, debugging, ties_method);  // single
        } else if (model_bool["basic"]) {
            Calc_LogLik_Basic(nthreads, RiskFail, RiskGroup, totalnum, ntime, R, Rd, Rdd, RdR, Rls1, Rls2, Rls3, Lls1, Lls2, Lls3, Ll, Lld, Lldd, debugging, ties_method, KeepConstant);
        } else if (model_bool["linear_err"]) {
            Calc_LogLik_Linear_ERR(tform, nthreads, RiskFail, RiskGroup, totalnum, ntime, R, Rd, Rdd, RdR, RddR, Rls1, Rls2, Rls3, Lls1, Lls2, Lls3, Ll, Lld, Lldd, debugging, ties_method, KeepConstant);
        } else if (model_bool["gradient"]) {
            Calc_LogLik_Gradient(nthreads, RiskFail, RiskGroup, totalnum, ntime, R, Rd, RdR, Rls1, Rls2, Lls1, Lls2, Ll, Lld, debugging, ties_method, KeepConstant);
        } else {
            Calc_LogLik(nthreads, RiskFail, RiskGroup, totalnum, ntime, R, Rd, Rdd, RdR, RddR, Rls1, Rls2, Rls3, Lls1, Lls2, Lls3, Ll, Lld, Lldd, debugging, ties_method, KeepConstant);
        }
    }
    // end_point = system_clock::now();
    // ending = time_point_cast<microseconds>(end_point).time_since_epoch().count();
//    Rcout << "C++ Note: log-lik " << (ending-comp) * 1e-6  <<endl;
    //
    if (model_bool["single"]) {
        iter_stop = 1;
        // if (verbose >= 4) {
        //     Rcout << "C++ Note: df101 ";  // prints the log-likelihoods
        //     for (int ij = 0; ij < reqrdnum; ij++) {
        //         Rcout << Ll[ij] << " ";
        //     }
        //     Rcout << " " << endl;
        //     Rcout << "C++ Note: df104 ";  // prints parameter values
        //     for (int ij = 0; ij < totalnum; ij++) {
        //         Rcout << beta_0[ij] << " ";
        //     }
        //     Rcout << " " << endl;
        // }
    } else {
        if (verbose >= 4) {
            Rcout << "C++ Note: df101 ";  // prints the log-likelihoods
            for (int ij = 0; ij < reqrdnum; ij++) {
                Rcout << Ll[ij] << " ";
            }
            Rcout << " " << endl;
            Rcout << "C++ Note: df102 ";  // prints the first derivatives
            for (int ij = 0; ij < reqrdnum; ij++) {
                Rcout << Lld[ij] << " ";
            }
            Rcout << " " << endl;
            if (!model_bool["gradient"]){
               Rcout << "C++ Note: df103 ";  // prints the second derivatives
               for (int ij = 0; ij < reqrdnum; ij++) {
                   Rcout << Lldd[ij*reqrdnum+ij] << " ";
               }
               Rcout << " " << endl;
               Rcout << "C++ Note: ALL df103 ";  // prints the second derivatives
               for (int ijk = 0; ijk < reqrdnum*reqrdnum; ijk++) {
                   Rcout << Lldd[ijk] << " ";
               }
                Rcout << " " << endl;
            }
            Rcout << "C++ Note: df104 ";  // prints parameter values
            for (int ij = 0; ij < totalnum; ij++) {
                Rcout << beta_0[ij] << " ";
            }
            Rcout << " " << endl;
        }
    }
//    if (verbose >= 4){
//        Rcout << "C++ Note: df101 ";  // prints the log-likelihoods
//        for (int ij = 0; ij < reqrdnum; ij++) {
//            Rcout << Ll[ij] << " ";
//        }
//        Rcout << " " << endl;
//        Rcout << "C++ Note: df102 ";  // prints the first derivatives
//        for (int ij = 0; ij < reqrdnum; ij++) {
//            Rcout << Lld[ij] << " ";
//        }
//        Rcout << " " << endl;
//        Rcout << "C++ Note: df104 ";  // prints parameter values
//        for (int ij = 0; ij < totalnum; ij++) {
//            Rcout << beta_0[ij] << " ";
//        }
//        Rcout << " " << endl;
//    }
}

//' Utility function to perform calculation of terms and risks for Poisson Omnibus
//'
//' \code{Pois_Term_Risk_Calc} Called to perform repeated term and risk calculations
//' @inheritParams CPP_template
//'
//' @return Updates matrices in place: risk storage matrices
//' @noRd
//'
// [[Rcpp::export]]
void Pois_Term_Risk_Calc(string modelform, const StringVector& tform, const IntegerVector& term_n, const int& totalnum, const int& fir, const IntegerVector& dfc, int term_tot, MatrixXd& T0, MatrixXd& Td0, MatrixXd& Tdd0, MatrixXd& Te, MatrixXd& R, MatrixXd& Rd, MatrixXd& Rdd, MatrixXd& Dose, MatrixXd& nonDose, VectorXd beta_0, const  MatrixXd& df0, const double& dint, const double& dslp, MatrixXd& TTerm, MatrixXd& nonDose_LIN, MatrixXd& nonDose_PLIN, MatrixXd& nonDose_LOGLIN, MatrixXd& RdR, MatrixXd& RddR, const MatrixXd& s_weights, const int& nthreads, bool debugging, const IntegerVector& KeepConstant, int verbose, List& model_bool, const double gmix_theta, const IntegerVector& gmix_term) {
    int reqrdnum = totalnum - sum(KeepConstant);
    if (model_bool["single"]) {
        // Calculates the subterm and term values
        Make_subterms_Single(totalnum, term_n, tform, dfc, fir, T0, Dose, nonDose, TTerm, nonDose_LIN, nonDose_PLIN, nonDose_LOGLIN, beta_0, df0, nthreads, debugging, KeepConstant);
        // ---------------------------------------------------------
        // Prints off a series of calculations to check at what point values are changing
        // ---------------------------------------------------------
        //
        // if (verbose >= 4) {
        //     Rcout << "C++ Note: values checked ";
        //     for (int ijk = 0; ijk < totalnum; ijk++) {
        //         Rcout << beta_0[ijk] << " ";
        //     }
        //     Rcout << " " << endl;
        //     Rcout << "C++ Note: sums checked ";
        //     for (int ijk = 0; ijk < totalnum; ijk++) {
        //         Rcout << T0.col(ijk).sum() << " ";
        //     }
        //     Rcout << " " << endl;
        //     Rcout << "C++ Note: dose checked ";
        //     for (int ijk = 0; ijk < term_tot; ijk++) {
        //         Rcout << Dose.col(ijk).array().sum() << " ";
        //     }
        //     Rcout << " " << endl;
        //     Rcout << "C++ Note: non-dose checked ";
        //     for (int ijk = 0; ijk < term_tot; ijk++) {
        //         Rcout << nonDose.col(ijk).array().sum() << " ";
        //     }
        //     Rcout << " " << endl;
        //     Rcout << "C++ Note: LIN_non-dose checked ";
        //     for (int ijk = 0; ijk < term_tot; ijk++) {
        //         Rcout << nonDose_LIN.col(ijk).array().sum() << " ";
        //     }
        //     Rcout << " " << endl;
        //     Rcout << "C++ Note: PLIN_non-dose checked ";
        //     for (int ijk = 0; ijk < term_tot; ijk++) {
        //         Rcout << nonDose_PLIN.col(ijk).array().sum() << " ";
        //     }
        //     Rcout << " " << endl;
        //     Rcout << "C++ Note: LOGLIN_non-dose checked ";
        //     for (int ijk = 0; ijk < term_tot; ijk++) {
        //         Rcout << nonDose_LOGLIN.col(ijk).array().sum() << " ";
        //     }
        //     Rcout << " " << endl;
        // }
        //
        // Calculates the risk for each row
        if (model_bool["strata"]) {
            Make_Risks_Weighted_Single(modelform, tform, term_n, totalnum, fir, s_weights, T0, Te, R, Dose, nonDose, TTerm, nonDose_LIN, nonDose_PLIN, nonDose_LOGLIN, nthreads, debugging, KeepConstant, gmix_theta, gmix_term);
        } else {
            Make_Risks_Single(modelform, tform, term_n, totalnum, fir, T0, Te, R, Dose, nonDose, TTerm, nonDose_LIN, nonDose_PLIN, nonDose_LOGLIN, nthreads, debugging, KeepConstant, gmix_theta, gmix_term);
        }
        //
        // Removes infinite values
        //
        //
        // if (R.minCoeff() <= 0) {
        //     if (verbose >= 4) {
        //         Rcout << "C++ Warning: risk mininum " << R.minCoeff() << " " << endl;
        //     }
        // } else if (verbose >= 4) {
        //     Rcout << "C++ Note: risk checked ";
        //     for (int ijk = 0; ijk < 1; ijk++) {
        //         Rcout << R.col(0).sum() << " ";
        //     }
        //     Rcout << " " << endl;
        //     //
        // }
    } else if (model_bool["gradient"]) {
        Make_subterms_Gradient(totalnum, term_n, tform, dfc, fir, T0, Td0, nonDose, TTerm, nonDose_LIN, nonDose_PLIN, nonDose_LOGLIN, beta_0, df0, nthreads, debugging, KeepConstant);
        //
        if (model_bool["strata"]) {
            Make_Risks_Weighted_Gradient(modelform, tform, term_n, totalnum, fir, s_weights, T0, Td0, Te, R, Rd, nonDose, TTerm, nonDose_LIN, nonDose_PLIN, nonDose_LOGLIN, RdR, nthreads, debugging, KeepConstant);
        } else {
            Make_Risks_Gradient(modelform, tform, term_n, totalnum, fir, T0, Td0, Te, R, Rd, nonDose, TTerm, nonDose_LIN, nonDose_PLIN, nonDose_LOGLIN, RdR, nthreads, debugging, KeepConstant);
        }
        //
        RdR = (RdR.array().isFinite()).select(RdR, 0);
    } else {
        //
        // Calculates the subterm and term values
        //
        Make_subterms(totalnum, term_n, tform, dfc, fir, T0, Td0, Tdd0, Dose, nonDose, TTerm, nonDose_LIN, nonDose_PLIN, nonDose_LOGLIN, beta_0, df0, dint, dslp, nthreads, debugging, KeepConstant);
        // ---------------------------------------------------------
        // Prints off a series of calculations to check at what point values are changing
        // ---------------------------------------------------------
        //
        //
        // if (verbose >= 4) {
        //     Rcout << "C++ Note: values checked ";
        //     for (int ijk = 0; ijk < totalnum; ijk++) {
        //         Rcout << beta_0[ijk] << " ";
        //     }
        //     Rcout << " " << endl;
        //     Rcout << "C++ Note: sums checked ";
        //     for (int ijk = 0; ijk < totalnum; ijk++) {
        //         Rcout << T0.col(ijk).sum() << " ";
        //     }
        //     Rcout << " " << endl;
        //     Rcout << "C++ Note: derivs checked ";
        //     for (int ijk = 0; ijk < reqrdnum; ijk++) {
        //         Rcout << Td0.col(ijk).sum() << " ";
        //     }
        //     Rcout << " " << endl;
        //     Rcout << "C++ Note: second derivs checked ";
        //     for (int ijk = 0; ijk < reqrdnum; ijk++) {
        //         Rcout << Tdd0.col(ijk*(ijk + 1)/2+ijk).sum() << " ";
        //     }
        //     Rcout << " " << endl;
        //     //
        //     Rcout << "C++ Note: ALL second derivs checked ";
        //     for (int ijk = 0; ijk < reqrdnum*(reqrdnum + 1)/2; ijk++) {
        //         Rcout << Tdd0.col(ijk).sum() << " ";
        //     }
        //     Rcout << " " << endl;
        //     //
        //     Rcout << "C++ Note: dose checked ";
        //     for (int ijk = 0; ijk < term_tot; ijk++) {
        //         Rcout << Dose.col(ijk).array().sum() << " ";
        //     }
        //     Rcout << " " << endl;
        //     Rcout << "C++ Note: non-dose checked ";
        //     for (int ijk = 0; ijk < term_tot; ijk++) {
        //         Rcout << nonDose.col(ijk).array().sum() << " ";
        //     }
        //     Rcout << " " << endl;
        //     Rcout << "C++ Note: LIN_non-dose checked ";
        //     for (int ijk = 0; ijk < term_tot; ijk++) {
        //         Rcout << nonDose_LIN.col(ijk).array().sum() << " ";
        //     }
        //     Rcout << " " << endl;
        //     Rcout << "C++ Note: PLIN_non-dose checked ";
        //     for (int ijk = 0; ijk < term_tot; ijk++) {
        //         Rcout << nonDose_PLIN.col(ijk).array().sum() << " ";
        //     }
        //     Rcout << " " << endl;
        //     Rcout << "C++ Note: LOGLIN_non-dose checked ";
        //     for (int ijk = 0; ijk < term_tot; ijk++) {
        //         Rcout << nonDose_LOGLIN.col(ijk).array().sum() << " ";
        //     }
        //     Rcout << " " << endl;
        // }
        // Calculates the risk for each row
        if (model_bool["strata"]) {
            Make_Risks_Weighted(modelform, tform, term_n, totalnum, fir, s_weights, T0, Td0, Tdd0, Te, R, Rd, Rdd, Dose, nonDose, TTerm, nonDose_LIN, nonDose_PLIN, nonDose_LOGLIN, RdR, RddR, nthreads, debugging, KeepConstant, gmix_theta, gmix_term);
        } else {
            Make_Risks(modelform, tform, term_n, totalnum, fir, T0, Td0, Tdd0, Te, R, Rd, Rdd, Dose, nonDose, TTerm, nonDose_LIN, nonDose_PLIN, nonDose_LOGLIN, RdR, RddR, nthreads, debugging, KeepConstant, gmix_theta, gmix_term);
        }
        //
        // Removes infinite values
        RdR = (RdR.array().isFinite()).select(RdR, 0);
        RddR = (RddR.array().isFinite()).select(RddR, 0);
        //
        //
        if (R.minCoeff() <= 0) {
            if (verbose >= 4) {
                Rcout << "C++ Warning: risk mininum " << R.minCoeff() << " " << endl;
            }
        } else if (verbose >= 4) {
            Rcout << "C++ Note: risk checked ";
            for (int ijk = 0; ijk < 1; ijk++) {
                Rcout << R.col(0).sum() << " ";
            }
            Rcout << " " << endl;
            Rcout << "C++ Note: risk1 checked ";
            for (int ijk = 0; ijk < reqrdnum; ijk++) {
                Rcout << Rd.col(ijk).sum() << " ";
            }
            Rcout << " " << endl;
            Rcout << "C++ Note: rdr checked ";
            for (int ijk = 0; ijk < reqrdnum; ijk++) {
                Rcout << RdR.col(ijk).sum() << " ";
            }
            Rcout << " " << endl;
            Rcout << "C++ Note: risk2 checked ";
            for (int ijk = 0; ijk < reqrdnum; ijk++) {
                Rcout << Rdd.col(ijk*(ijk + 1)/2+ijk).sum() << " ";
            }
            Rcout << " " << endl;
            //
            Rcout << "C++ Note: ALL risk2 checked ";
            for (int ijk = 0; ijk < reqrdnum*(reqrdnum + 1)/2; ijk++) {
               Rcout << Rdd.col(ijk).sum() << " ";
            }
            Rcout << " " << endl;
            //
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
// [[Rcpp::export]]
void Pois_Dev_LL_Calc(const int& reqrdnum, const int& totalnum, const int& fir, MatrixXd& R, MatrixXd& Rd, MatrixXd& Rdd, VectorXd beta_0, MatrixXd& RdR, MatrixXd& RddR, vector<double>& Ll, vector<double>& Lld, vector<double>& Lldd, const MatrixXd& PyrC, MatrixXd& dev_temp, const int& nthreads, bool debugging, const IntegerVector& KeepConstant, int verbose, List& model_bool, int iter_stop, double& dev) {
    fill(Ll.begin(), Ll.end(), 0.0);
    if (!model_bool["single"]) {
        fill(Lld.begin(), Lld.end(), 0.0);
        if (!model_bool["gradient"]) {
            fill(Lldd.begin(), Lldd.end(), 0.0);
        }
    }
    if (model_bool["single"]) {
        Poisson_LogLik_Single(nthreads, totalnum, PyrC, R, Ll, debugging);
    } else if (model_bool["gradient"]) {
        Poisson_LogLik_Gradient(nthreads, totalnum, PyrC, R, Rd, RdR, Ll, Lld, debugging, KeepConstant);
    } else {
        Poisson_LogLik(nthreads, totalnum, PyrC, R, Rd, Rdd, RdR, RddR, Ll, Lld, Lldd, debugging, KeepConstant);
    }
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
    dev_temp.col(0) = (2 * dev_temp.col(0).array()).array();  //.sqrt();
    dev_temp = (dev_temp.array().isFinite()).select(dev_temp, 0);
    dev_temp.col(0) = (R.col(0).array()<0).select(0, dev_temp.col(0));
    dev = dev_temp.col(0).sum();  // deviation calculation is split into steps
    //
    if (model_bool["single"]) {
        iter_stop = 1;
        // if (verbose >= 4) {
        //     Rcout << "C++ Note: df101 ";  // prints the log-likelihoods
        //     for (int ij = 0; ij < reqrdnum; ij++) {
        //         Rcout << Ll[ij] << " ";
        //     }
        //     Rcout << " " << endl;
        //     Rcout << "C++ Note: df104 ";  // prints parameter values
        //     for (int ij = 0; ij < totalnum; ij++) {
        //         Rcout << beta_0[ij] << " ";
        //     }
        //     Rcout << " " << endl;
        //     Rcout << "C++ Note: Checking Deviance " << dev << endl;
        // }
    } else {
        // if (verbose >= 4) {
        //     Rcout << "C++ Note: df101 ";  // prints the log-likelihoods
        //     for (int ij = 0; ij < reqrdnum; ij++) {
        //         Rcout << Ll[ij] << " ";
        //     }
        //     Rcout << " " << endl;
        //     Rcout << "C++ Note: df102 ";  // prints the first derivatives
        //     for (int ij = 0; ij < reqrdnum; ij++) {
        //         Rcout << Lld[ij] << " ";
        //     }
        //     Rcout << " " << endl;
        //     Rcout << "C++ Note: df103 ";  // prints the second derivatives
        //     for (int ij = 0; ij < reqrdnum; ij++) {
        //         Rcout << Lldd[ij*reqrdnum+ij] << " ";
        //     }
        //     Rcout << " " << endl;
        //     Rcout << "C++ Note: df104 ";  // prints parameter values
        //     for (int ij = 0; ij < totalnum; ij++) {
        //         Rcout << beta_0[ij] << " ";
        //     }
        //     Rcout << " " << endl;
        //     Rcout << "C++ Note: Checking Deviance " << dev << endl;
        // }
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
// [[Rcpp::export]]
void Cox_Pois_Check_Continue(List& model_bool, VectorXd beta_0, vector<double>& beta_best, vector<double>& beta_c, const VectorXd& cens_weight, const bool change_all, vector<double>& dbeta, const bool debugging, double& dev, MatrixXd& dev_temp, const int fir, const int halfmax, double& halves, int& ind0, int& iter_stop, const IntegerVector& KeepConstant, vector<double>& Ll, double& Ll_abs_best, vector<double>& Lld, vector<double>& Lldd, MatrixXd& Lls1, MatrixXd& Lls2, MatrixXd& Lls3, const double& Lstar, const int& nthreads, const int& ntime, const MatrixXd& PyrC, MatrixXd& R, MatrixXd& Rd, MatrixXd& Rdd, MatrixXd& RddR, MatrixXd& RdR, const int& reqrdnum, const StringVector& tform, const IntegerMatrix& RiskFail, const vector<string>& RiskGroup, const StringMatrix& RiskGroup_Strata, MatrixXd& Rls1, MatrixXd& Rls2, MatrixXd& Rls3, NumericVector& Strata_vals, const IntegerVector& term_n, const string ties_method, const int totalnum, MatrixXd& TTerm, const int verbose) {
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
        // if (verbose >= 4) {
        //     Rcout << "C++ Warning: A non-positive risk was detected: " << R.minCoeff() << endl;
        // }
        halves+=0.2;
    } else {
        halves++;
        if (model_bool["cox"]) {
            Cox_Side_LL_Calc(reqrdnum, ntime, tform, RiskFail, RiskGroup_Strata, RiskGroup, totalnum, fir, R, Rd, Rdd, Rls1, Rls2, Rls3, Lls1, Lls2, Lls3, cens_weight, Strata_vals, beta_0, RdR, RddR, Ll, Lld, Lldd, nthreads, debugging, KeepConstant, ties_method, verbose, model_bool, iter_stop);
        } else {
            Pois_Dev_LL_Calc(reqrdnum, totalnum, fir, R, Rd, Rdd, beta_0, RdR, RddR, Ll, Lld, Lldd, PyrC, dev_temp, nthreads, debugging, KeepConstant, verbose, model_bool, iter_stop, dev);
        }
        //
        if (model_bool["log_bound"]) {
            if (Ll[0] > Lstar) {
                // If it has gone beyond Lstar, then this isn't the point
                iter_stop = 1;
                halves = halfmax;
            } else {
            //
                if (Ll[ind0] <= Ll_abs_best) {  // if a better point wasn't found, takes a half-step
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
        } else {
            if (change_all) {  // if every covariate is to be changed
                if (Ll[ind0] <= Ll_abs_best) {  // if a better point wasn't found, takes a half-step
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
            } else {  // for validation, the step is always carried over
                // used if a single parameter is being changed to trick program
                #ifdef _OPENMP
                #pragma omp parallel for schedule(dynamic) num_threads(nthreads)
                #endif
                for (int ijk = 0; ijk < totalnum; ijk++) {
                    beta_best[ijk] = beta_c[ijk];
                }
            }
        }
        #ifdef _OPENMP
        #pragma omp parallel for schedule(dynamic) num_threads(nthreads)
        #endif
        for (int ijk = 0; ijk < totalnum; ijk++) {  // totalnum*(totalnum + 1)/2
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
// [[Rcpp::export]]
void Cox_Pois_Log_Loop(double& abs_max, List& model_bool, VectorXd beta_0, vector<double>& beta_a, vector<double>& beta_c, int& bound_val, vector<double>& dbeta, const bool debugging, const MatrixXd& df0, IntegerVector& dfc, double& dint, MatrixXd& Dose, double& dose_abs_max, double& dslp, const int fir, const IntegerVector& gmix_term, const double& gmix_theta, int& half_check, const int halfmax, const IntegerVector& KeepConstant, vector<bool>& limit_hit, double& lr, string& modelform, MatrixXd& nonDose, MatrixXd& nonDose_LIN, MatrixXd& nonDose_LOGLIN, MatrixXd& nonDose_PLIN, const int& nthreads, MatrixXd& R, MatrixXd& Rd, MatrixXd& Rdd, MatrixXd& RddR, MatrixXd& RdR, VectorXd& s_weights, MatrixXd& T0, MatrixXd& Td0, MatrixXd& Tdd0, MatrixXd& Te, const IntegerVector& term_n, int& term_tot, StringVector& tform, const int totalnum, MatrixXd& TTerm, const int verbose) {
    while ((R.minCoeff() <= 0) || (R.hasNaN())) {
        half_check++;
        if (half_check>halfmax) {
            limit_hit[bound_val] = TRUE;
            break;
        } else {
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
            for (int ij = 0; ij < totalnum; ij++) {
                beta_0[ij] = beta_a[ij] + lr*dbeta[ij];
                beta_c[ij] = beta_0[ij];
            }
            if (model_bool["cox"]) {
                Cox_Term_Risk_Calc(modelform, tform, term_n, totalnum, fir, dfc, term_tot, T0, Td0, Tdd0, Te, R, Rd, Rdd, Dose, nonDose, beta_0, df0, dose_abs_max, abs_max, TTerm, nonDose_LIN, nonDose_PLIN, nonDose_LOGLIN, RdR, RddR, nthreads, debugging, KeepConstant, verbose, model_bool, gmix_theta, gmix_term);
            } else {
                Pois_Term_Risk_Calc(modelform, tform, term_n, totalnum, fir, dfc, term_tot, T0, Td0, Tdd0, Te, R, Rd, Rdd, Dose, nonDose, beta_0, df0, dint, dslp, TTerm, nonDose_LIN, nonDose_PLIN, nonDose_LOGLIN, RdR, RddR, s_weights, nthreads, debugging, KeepConstant, verbose, model_bool, gmix_theta, gmix_term);
            }
        }
    }
    return;
}


// //' Utility function to calculate Cox Log-Likelihood and derivatives
// //'
// //' \code{Simplified_Inform_Matrix} Called to update log-likelihoods, Uses list of event rows, risk matrices, and repeated sums, Sums the log-likelihood contribution from each event time
// //' @inheritParams CPP_template
// //'
// //' @return Updates matrices in place: Log-likelihood vectors/matrix
// //' @noRd
// //'
// // [[Rcpp::export]]
// void Simplified_Inform_Matrix(const int& nthreads, const IntegerMatrix& RiskFail, const vector<string>&  RiskGroup, const int& totalnum, const int& ntime, const MatrixXd& R, const MatrixXd& Rd, const MatrixXd& Rdd, const MatrixXd& RdR, const MatrixXd& RddR, const MatrixXd& Rls1, const MatrixXd& Rls2, const MatrixXd& Rls3, const MatrixXd& Lls1, const MatrixXd& Lls2, const MatrixXd& Lls3, vector<double>& InMa, bool debugging, string ties_method, const IntegerVector& KeepConstant) {
//     int reqrdnum = totalnum - sum(KeepConstant);
//     #ifdef _OPENMP
//     #pragma omp declare reduction(vec_double_plus : std::vector<double> : \
//         std::transform(omp_out.begin(), omp_out.end(), omp_in.begin(), omp_out.begin(), std::plus<double>())) \
//         initializer(omp_priv = omp_orig)
//     #pragma omp parallel for schedule(dynamic) num_threads(nthreads) reduction(vec_double_plus:InMa) collapse(2)
//     #endif
//     for (int ijk = 0; ijk < reqrdnum*(reqrdnum + 1)/2; ijk++) {  // performs log-likelihood calculations for every derivative combination and risk group
//         for (int j = 0; j < ntime; j++) {
//             // rcout << ijk << ", " << j << ", " << " start" << endl;
//             int ij = 0;
//             int jk = ijk;
//             while (jk > ij) {
//                 ij++;
//                 jk -= ij;
//             }
//             double Rs1 = Rls1(j, 0);
//             double Rs2 = Rls2(j, ij);
//             double Rs2t = Rls2(j, jk);
//             double Rs3 = Rls3(j, ijk);
//             //
//             int dj = RiskFail(j, 1)-RiskFail(j, 0) + 1;
//             MatrixXd Ld = MatrixXd::Zero(dj, 4);
//             // rcout << 0 << endl;
//             // ld << R.block(RiskFail(j, 0), 0, dj, 1), RdR.block(RiskFail(j, 0), ij, dj, 1), RdR.block(RiskFail(j, 0), jk, dj, 1), RddR.block(RiskFail(j, 0), ijk, dj, 1);  // rows with events
//             //
//             MatrixXd Ldm = MatrixXd::Zero(dj, 4);
//             Vector4d Ldcs;
//             if (ties_method == "efron") {
//                 Ldcs << Lls1(j, 0), Lls2(j, ij), Lls2(j, jk), Lls3(j, ijk);
//                 for (int i = 0; i < dj; i++) {  // adds in the efron approximation terms
//                     Ldm.row(i) = (-double(i) / double(dj)) *Ldcs.array();
//                 }
//             }
//             Ldm.col(0) = Ldm.col(0).array() + Rs1;
//             Ldm.col(1) = Ldm.col(1).array() + Rs2;
//             Ldm.col(2) = Ldm.col(2).array() + Rs2t;
//             Ldm.col(3) = Ldm.col(3).array() + Rs3;
//             // Calculates the left-hand side terms
//             //
//             double Ld1;
//             double Ld3;
//             //
//             MatrixXd temp1 = MatrixXd::Zero(dj, 1);
//             MatrixXd temp2 = MatrixXd::Zero(dj, 1);
//             // calculates the right-hand side terms
//             temp1 = Ldm.col(1).array() * (Ldm.col(0).array().pow(- 1).array());
//             temp2 = Ldm.col(2).array() * (Ldm.col(0).array().pow(- 1).array());
//             temp1 = temp1.array() * temp2.array();
//             Rs3 = (temp1.array().isFinite()).select(temp1, 0).sum();
//             //
//             vector<int> InGroup;
//             string Groupstr = RiskGroup[j];
//             stringstream ss(Groupstr);
//             //
//             //
//             //
//             for (int i; ss >> i;) {
//                 InGroup.push_back(i);
//                 if (ss.peek() == ',')
//                     ss.ignore();
//             }
//             Ld3 = 0;
//             Ld1 = 0;
//             // rcout << 1 << endl;
//             // now has the grouping pairs
//             for (vector<double>::size_type i = 0; i < InGroup.size() - 1; i = i + 2) {
//                 // rcout << 2 << ", " << i << endl;
//                 Ld3 += (RdR.block(InGroup[i] - 1, ij, InGroup[i + 1]-InGroup[i] + 1, 1).array() * RdR.block(InGroup[i] - 1, jk, InGroup[i + 1]-InGroup[i] + 1, 1).array()* R.block(InGroup[i] - 1, 0, InGroup[i + 1]-InGroup[i] + 1, 1).array()).sum();
//                 Ld1 += R.block(InGroup[i] - 1, 0, InGroup[i + 1]-InGroup[i] + 1, 1).sum();
//             }
//             Ld3 = Ld3 / Ld1;
//             // rcout << 3 << endl;
//             //
//             InMa[ij*reqrdnum+jk] += Ld3 - Rs3;  // sums the log-likelihood and derivatives
//             // rcout << 4 << endl;
//         }
//     }
//     #ifdef _OPENMP
//     #pragma omp parallel for schedule(dynamic) num_threads(nthreads)
//     #endif
//     for (int ijk = 0; ijk < reqrdnum*(reqrdnum + 1)/2; ijk++) {  // fills second-derivative matrix
//         int ij = 0;
//         int jk = ijk;
//         while (jk > ij) {
//             ij++;
//             jk -= ij;
//         }
//         InMa[jk*reqrdnum+ij] = InMa[ij*reqrdnum+jk];
//     }
//     return;
// }

// //' Utility function to calculate Cox Log-Likelihood and derivatives
// //'
// //' \code{Simplified_Inform_Matrix_Strata} Called to update log-likelihoods, Uses list of event rows, risk matrices, and repeated sums, Sums the log-likelihood contribution from each event time
// //' @inheritParams CPP_template
// //'
// //' @return Updates matrices in place: Log-likelihood vectors/matrix
// //' @noRd
// //'
// // [[Rcpp::export]]
// void Simplified_Inform_Matrix_Strata(const int& nthreads, const IntegerMatrix& RiskFail, const StringMatrix&  RiskGroup, const int& totalnum, const int& ntime, const MatrixXd& R, const MatrixXd& Rd, const MatrixXd& Rdd, const MatrixXd& RdR, const MatrixXd& RddR, const MatrixXd& Rls1, const MatrixXd& Rls2, const MatrixXd& Rls3, const MatrixXd& Lls1, const MatrixXd& Lls2, const MatrixXd& Lls3, vector<double>& InMa, bool debugging, string ties_method, NumericVector& Strata_vals, const IntegerVector& KeepConstant) {
//     int reqrdnum = totalnum - sum(KeepConstant);
//     #ifdef _OPENMP
//     #pragma omp declare reduction(vec_double_plus : std::vector<double> : \
//         std::transform(omp_out.begin(), omp_out.end(), omp_in.begin(), omp_out.begin(), std::plus<double>())) \
//         initializer(omp_priv = omp_orig)
//     #pragma omp parallel for schedule(dynamic) num_threads(nthreads) reduction(vec_double_plus:InMa) collapse(3)
//     #endif
//     for (int ijk = 0; ijk < reqrdnum*(reqrdnum + 1)/2; ijk++) {  // performs log-likelihood calculations for every derivative combination and risk group
//         for (int j = 0; j < ntime; j++) {
//             for (int s_ij = 0; s_ij < Strata_vals.size(); s_ij++) {
//                 int ij = 0;
//                 int jk = ijk;
//                 while (jk > ij) {
//                     ij++;
//                     jk -= ij;
//                 }
//                 double Rs1 = Rls1(j, 0);
//                 double Rs2 = Rls2(j, ij*Strata_vals.size() + s_ij);
//                 double Rs2t = Rls2(j, jk*Strata_vals.size() + s_ij);
//                 double Rs3 = Rls3(j, ijk*Strata_vals.size() + s_ij);
//                 //
//                 int dj = RiskFail(j, 2*s_ij + 1)-RiskFail(j, 2*s_ij + 0) + 1;
//                 if (RiskFail(j, 2*s_ij + 1)> - 1) {
//                     MatrixXd Ld = MatrixXd::Zero(dj, 4);
//                     // ld << R.block(RiskFail(j, 2*s_ij), 0, dj, 1), RdR.block(RiskFail(j, 2*s_ij), ij, dj, 1), RdR.block(RiskFail(j, 2*s_ij), jk, dj, 1), RddR.block(RiskFail(j, 2*s_ij), ijk, dj, 1);  // rows with events
//                     //
//                     MatrixXd Ldm = MatrixXd::Zero(dj, 4);
//                     Vector4d Ldcs;
//                     if (ties_method == "efron") {
//                         Ldcs << Lls1(j, s_ij), Lls2(j, ij*Strata_vals.size() + s_ij), Lls2(j, jk*Strata_vals.size() + s_ij), Lls3(j, ijk*Strata_vals.size() + s_ij);
//                         for (int i = 0; i < dj; i++) {  // adds in the efron approximation terms
//                             Ldm.row(i) = (-double(i) / double(dj)) *Ldcs.array();
//                         }
//                     }
//                     Ldm.col(0) = Ldm.col(0).array() + Rs1;
//                     Ldm.col(1) = Ldm.col(1).array() + Rs2;
//                     Ldm.col(2) = Ldm.col(2).array() + Rs2t;
//                     Ldm.col(3) = Ldm.col(3).array() + Rs3;
//                     // Calculates the left-hand side terms
//                     //
//                     double Ld1;
//                     double Ld3;
//                     //
//                     MatrixXd temp1 = MatrixXd::Zero(dj, 1);
//                     MatrixXd temp2 = MatrixXd::Zero(dj, 1);
//                     // calculates the right-hand side terms
//                     temp1 = Ldm.col(1).array() * (Ldm.col(0).array().pow(- 1).array());
//                     temp2 = Ldm.col(2).array() * (Ldm.col(0).array().pow(- 1).array());
//                     temp1 = temp1.array() * temp2.array();
//                     Rs3 = (temp1.array().isFinite()).select(temp1, 0).sum();
//                     //
//                     vector<int> InGroup;
//                     string Groupstr = as<std::string>(RiskGroup(j, s_ij));
//                     stringstream ss(Groupstr);
//                     for (int i; ss >> i;) {
//                         InGroup.push_back(i);
//                         if (ss.peek() == ',')
//                             ss.ignore();
//                     }
//                     Ld3 = 0;
//                     Ld1 = 0;
//                     // now has the grouping pairs
//                     for (vector<double>::size_type i = 0; i < InGroup.size() - 1; i = i + 2) {
//                         Ld3 += (RdR.block(InGroup[i] - 1, ij, InGroup[i + 1]-InGroup[i] + 1, 1).array() * RdR.block(InGroup[i] - 1, jk, InGroup[i + 1]-InGroup[i] + 1, 1).array()* R.block(InGroup[i] - 1, 0, InGroup[i + 1]-InGroup[i] + 1, 1).array()).sum();
//                         Ld1 += R.block(InGroup[i] - 1, 0, InGroup[i + 1]-InGroup[i] + 1, 1).sum();
//                     }
//                     Ld3 = Ld3 / Ld1;
//                     //
//                     InMa[ij*reqrdnum+jk] += Ld3 - Rs3;  // sums the log-likelihood and derivatives
//                 }
//             }
//         }
//     }
//     #ifdef _OPENMP
//     #pragma omp parallel for schedule(dynamic) num_threads(nthreads)
//     #endif
//     for (int ijk = 0; ijk < reqrdnum*(reqrdnum + 1)/2; ijk++) {  // fills second-derivative matrix
//         int ij = 0;
//         int jk = ijk;
//         while (jk > ij) {
//             ij++;
//             jk -= ij;
//         }
//         InMa[jk*reqrdnum+ij] = InMa[ij*reqrdnum+jk];
//     }
//     return;
// }
