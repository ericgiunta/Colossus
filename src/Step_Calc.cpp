//  Copyright 2022 - 2025, Eric Giunta and the project collaborators, Please see main R package for license and usage details

#include <RcppEigen.h>

#include "Step_Calc.h"
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
#include <set>

#include "Calc_Repeated.h"
#include "Subterms_Risk.h"
#include "Colossus_types.h"


//  [[Rcpp::depends(RcppEigen)]]
//  [[Rcpp::plugins(openmp)]]

using std::string;
using std::vector;
using std::advance;
using std::reduce;
using std::endl;
using std::isinf;
using std::isnan;
using std::set;

using Eigen::Map;
using Eigen::Ref;
using Eigen::MatrixXd;
using Eigen::SparseMatrix;
using Eigen::VectorXd;

using Rcpp::as;
using Rcpp::IntegerVector;
using Rcpp::NumericVector;
using Rcpp::NumericMatrix;
using Rcpp::StringVector;
using Rcpp::List;
using Rcpp::Rcout;
using Rcpp::_;
using Rcpp::Dimension;

template <typename T> int sign(T val) {
    return (T(0) < val) - (val < T(0));
}

//' Utility function to keep intercept parameters within the range of possible values
//'
//' \code{Intercept_Bound} Called to update the parameter list in the event that intercepts leave the bounds of possible values
//' @inheritParams CPP_template
//'
//' @return Updates vector in place: parameter vector
//' @noRd
//'
//
void Intercept_Bound(const int& nthreads, const int& totalnum, const VectorXd& beta_0, vector<double>& dbeta, const IntegerVector& dfc, const Ref<const MatrixXd>& df0, const IntegerVector& KeepConstant, const StringVector&  tform) {
    set<string> Dose_Iden;  //  list of dose subterms
    Dose_Iden.insert("lin_int");
    Dose_Iden.insert("step_int");
    Dose_Iden.insert("lin_quad_int");
    Dose_Iden.insert("lin_exp_int");
    #ifdef _OPENMP
    #pragma omp parallel for schedule(dynamic) num_threads(nthreads)
    #endif
    for (int ij = 0; ij < totalnum; ij++) {
        if ((Dose_Iden.find(as<string>(tform[ij])) != Dose_Iden.end()) && (KeepConstant[ij] == 0)) {
            int df0_c = dfc[ij] - 1;
            double pmin = (df0.col(df0_c)).array().minCoeff();
            double pmax = (df0.col(df0_c)).array().maxCoeff();
            double db_temp = beta_0[ij] + dbeta[ij];
            if (db_temp < pmin) {
                dbeta[ij] = pmin-beta_0[ij];
            } else if (db_temp > pmax) {
                dbeta[ij] = pmax-beta_0[ij];
            }
        }
    }
    return;
}

//' Utility function to calculate steps for a likelihood based bound
//'
//' \code{Log_Bound} Called to perform likelihood bound steps
//' @inheritParams CPP_template
//' @param Lstar likelihood goal
//' @param L0 current likelihood
//'
//' @return Updates matrices in place: risk storage matrices
//' @noRd
//'
//
void Log_Bound(double& deriv_max, const MatrixXd& Lldd_mat, const VectorXd& Lld_vec, const double& Lstar, const double& qchi, const double& L0, const int& para_number, const int& nthreads, const int& totalnum, const int& reqrdnum, IntegerVector KeepConstant, const int& term_tot, const int& step, vector<double>& dbeta, const VectorXd& beta_0, bool upper, bool& trouble, int verbose, double mult) {
    //  starts with solved likelihoods and derivatives
    //  store the second derivative as D0
    MatrixXd D0 = Lldd_mat;
    deriv_max = 100;
    if (step == 0) {
        //  initial step, calculate dom/dbet and h
        MatrixXd dOmdBeta = Lldd_mat.col(para_number).matrix();
        removeRow(D0, para_number);
        removeColumn(D0, para_number);
        removeRow(dOmdBeta, para_number);
        D0 = D0.inverse().matrix();
        dOmdBeta = - 1 * D0 * dOmdBeta;
        //
        MatrixXd dLdBdO = Lldd_mat.row(para_number).matrix();
        removeColumn(dLdBdO, para_number);
        double h = Lldd_mat(para_number, para_number) - (dLdBdO.matrix() * D0 * dLdBdO.matrix().transpose().matrix())(0, 0);
        h = mult * pow(qchi/(- 1*h), 0.5);
        if (upper) {
            h = abs(h)/2;
        } else {
            h = abs(h)/-2;
        }
        //  calculate first step
        int j = 0;
        for (int ij = 0; ij < totalnum; ij++) {
            if (KeepConstant[ij] == 0) {
                int pij_ind = ij - sum(head(KeepConstant, ij));
                if (pij_ind == para_number) {
                    dbeta[ij] = h;
                } else {
                    dbeta[ij] = h * dOmdBeta(j);
                    j = j + 1;
                }
            }
        }
    } else {
        MatrixXd G = MatrixXd::Zero(reqrdnum, reqrdnum);
        VectorXd v = VectorXd::Zero(reqrdnum);
        v[para_number] = L0 - Lstar;
        G.row(para_number) = Lld_vec;
        for (int j = 0;  j < reqrdnum; j++) {
            if (j != para_number) {
                G.row(j) = D0.row(j);
                v[j] = Lld_vec[j];
            }
        }
        //  At this point, we have the standard newton-raphson equation defined
        deriv_max = abs(v[0]);
        for (int ij = 0; ij < reqrdnum; ij++) {
            if (abs(v[ij]) > deriv_max) {
                deriv_max = abs(v[ij]);
            }
        }
        //
        if (abs(G.determinant()) < 1e-6) {
            //  The inverted matrix does not exist
            for (int ij = 0; ij < totalnum; ij++) {
                if (KeepConstant[ij] == 0) {
                    int pij_ind = ij - sum(head(KeepConstant, ij));
                    dbeta[ij] = -v[pij_ind]/G(pij_ind, pij_ind);
                }
            }
        } else {
            G = G.inverse().matrix();
            v = G.matrix() * v.matrix();
            VectorXd g1 = G.col(para_number);
            //  we now must solve for the roots
            double as2 = g1.matrix().transpose() * D0 * g1.matrix();
            double bs1 = 2*v.matrix().transpose() *D0 * g1.matrix() - 2;
            double cs0 = v.matrix().transpose() * D0 * v.matrix();
            //
            if (pow(bs1, 2)-4*as2*cs0 >= 0) {
                double s0 = pow(bs1, 2)-4*as2*cs0;
                double s1 = (-bs1 - pow(s0, 0.5))/(2*as2);
                s0 = (-bs1 + pow(s0, 0.5))/(2*as2);
                //  check which is closer
                double s00 = (v + s0*g1).matrix().transpose() * D0 * (v + s0*g1).matrix();
                double s11 = (v + s1*g1).matrix().transpose() * D0 * (v + s1*g1).matrix();
                //
                if (abs(s00) < abs(s11)) {
                    //  s1 is further away
                    for (int ij = 0; ij < totalnum; ij++) {
                        if (KeepConstant[ij] == 0) {
                            int pij_ind = ij - sum(head(KeepConstant, ij));
                            dbeta[ij] = -v[pij_ind] - g1[pij_ind]*s0;
                        }
                    }
                } else {
                    //  s0 is further away
                    for (int ij = 0; ij < totalnum; ij++) {
                        if (KeepConstant[ij] == 0) {
                            int pij_ind = ij - sum(head(KeepConstant, ij));
                            dbeta[ij] = -v[pij_ind] - g1[pij_ind]*s1;
                        }
                    }
                }
            } else {
                //  there are no real solutions, needs a more conservative step?
                //  currently will take a step to optimize with constant beta
                trouble = true;
                double s0 = -bs1/2/as2;
                for (int ij = 0; ij < totalnum; ij++) {
                    if (KeepConstant[ij] == 0) {
                        int pij_ind = ij - sum(head(KeepConstant, ij));
                        dbeta[ij] = -v[pij_ind] - g1[pij_ind]*s0;
                    }
                }
            }
        }
    }
    return;
}

//' Utility function to calculate the change to make each iteration
//'
//' \code{Calc_Change_trouble} Called to update the parameter changes, Uses log-likelihoods and control parameters, Applies newton steps and change limitations
//' @inheritParams CPP_template
//'
//' @return Updates matrices in place: parameter change matrix
//' @noRd
//'
//
void Calc_Change_trouble(const int& para_number, const int& nthreads, const int& totalnum, const double& thres_step_max, const double& lr, const double& step_max, const vector<double>& Ll, const vector<double>& Lld, const vector<double>& Lldd, vector<double>& dbeta, const StringVector&   tform, const double& dint, const double& dslp, IntegerVector KeepConstant_trouble) {
    int kept_covs = totalnum - sum(KeepConstant_trouble);
    NumericVector Lldd_vec(kept_covs * kept_covs);
    NumericVector Lld_vec(kept_covs);
    #ifdef _OPENMP
    #pragma omp parallel for schedule(dynamic) num_threads(nthreads)
    #endif
    for (int ijk = 0; ijk < kept_covs*(kept_covs + 1)/2; ijk++) {
        int ij = 0;
        int jk = ijk;
        while (jk > ij) {
            ij++;
            jk -= ij;
        }
        int ij0 = ij;
        int jk0 = jk;
        if (ij >= para_number) {
            ij0++;
        }
        if (jk >= para_number) {
            jk0++;
        }
        Lldd_vec[jk * kept_covs + ij] = Lldd[jk0 * (kept_covs + 1) + ij0];
        if (ij == jk) {
            Lld_vec[ij] = Lld[ij0];
        } else {
            Lldd_vec[ij * kept_covs + jk] = Lldd_vec[jk0 * (kept_covs + 1) + ij0];
        }
    }
    Lldd_vec.attr("dim") = Dimension(kept_covs, kept_covs);
    const Map<MatrixXd> Lldd_mat(as<Map<MatrixXd> >(Lldd_vec));
    const Map<VectorXd> Lld_mat(as<Map<VectorXd> >(Lld_vec));
    VectorXd Lldd_solve0 = Lldd_mat.colPivHouseholderQr().solve(- 1*Lld_mat);
    VectorXd Lldd_solve = VectorXd::Zero(totalnum);
    for (int ij = 0; ij < totalnum; ij++) {
        if (KeepConstant_trouble[ij] == 0) {
            int pij_ind = ij - sum(head(KeepConstant_trouble, ij));
            Lldd_solve(ij) = Lldd_solve0(pij_ind);
        }
    }
    //
    #ifdef _OPENMP
    #pragma omp parallel for schedule(dynamic) num_threads(nthreads)
    #endif
    for (int ijk = 0; ijk < totalnum; ijk++) {
        if (KeepConstant_trouble[ijk] == 0) {
            int pjk_ind = ijk - sum(head(KeepConstant_trouble, ijk));
            if (isnan(Lldd_solve(ijk))) {
                if (Lldd[pjk_ind*kept_covs+pjk_ind] != 0) {
                    dbeta[ijk] = -lr * Lld[pjk_ind] / Lldd[pjk_ind*kept_covs+pjk_ind];
                } else {
                    dbeta[ijk] = 0;
                }
            } else {
                dbeta[ijk] = lr * Lldd_solve(ijk);
            }
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
    return;
}

