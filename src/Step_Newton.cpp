//  Copyright 2022 - 2025, Eric Giunta and the project collaborators, Please see main R package for license and usage details

#include <RcppEigen.h>

#include "Step_Newton.h"
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

//' Utility function to calculate the change to make each iteration, applying linear constraints
//'
//' \code{Calc_Change_Cons} Called to update the parameter changes, Uses log-likelihoods and control parameters, Applies newton steps and change limitations with a system of constraints
//' @inheritParams CPP_template
//'
//' @return Updates matrices in place: parameter change matrix
//' @noRd
//'
//
void Calc_Change_Cons(const MatrixXd& Lin_Sys, const VectorXd& Lin_Res, const  VectorXd& beta_0, const int& nthreads, const int& totalnum, const double& thres_step_max, const double& lr, const double& step_max, const vector<double>& Ll, const vector<double>& Lld, const vector<double>& Lldd, vector<double>& dbeta, const StringVector&   tform, const double& dint, const double& dslp, IntegerVector KeepConstant) {
    //
    int kept_covs = totalnum - sum(KeepConstant);
    //
    VectorXd beta_1(kept_covs);
    for (int ij = 0; ij < totalnum; ij++) {
        if (KeepConstant[ij] == 0) {
            int pij_ind = ij - sum(head(KeepConstant, ij));
            beta_1(pij_ind) = beta_0(ij);
        }
    }
    VectorXd Lin_Dif = Lin_Sys * beta_1 - Lin_Res;
    //
    int total_covs = kept_covs + Lin_Sys.rows();
    //
    NumericVector Lldd_vec(total_covs*total_covs);
    NumericVector Lld_vec(total_covs);
    #ifdef _OPENMP
    #pragma omp parallel for schedule(dynamic) num_threads(nthreads)
    #endif
    for (int ijk = 0; ijk < total_covs*(total_covs + 1)/2; ijk++) {
        int ij = 0;
        int jk = ijk;
        while (jk > ij) {
            ij++;
            jk -= ij;
        }
        if (ij < kept_covs) {
            Lldd_vec[jk * total_covs + ij] = Lldd[jk * kept_covs + ij];
            if (ij == jk) {
                Lld_vec[ij] = Lld[ij];
            } else {
                Lldd_vec[ij * total_covs + jk] = Lldd_vec[jk * kept_covs + ij];
            }
        } else {
            if (jk < kept_covs) {
                Lldd_vec[jk * total_covs + ij] = Lin_Sys(ij-kept_covs, jk);
            } else {
                Lldd_vec[jk * total_covs + ij] = 0.0;
            }
            if (ij == jk) {
                Lld_vec[ij] = Lin_Dif(ij-kept_covs);
            } else {
                Lldd_vec[ij * total_covs + jk] = Lldd_vec[jk * total_covs + ij];
            }
        }
    }
    //
    //
    Lldd_vec.attr("dim") = Dimension(total_covs, total_covs);
    const Map<MatrixXd> Lldd_mat(as<Map<MatrixXd> >(Lldd_vec));
    const Map<VectorXd> Lld_mat(as<Map<VectorXd> >(Lld_vec));
    //
    //
    VectorXd Lldd_solve0 = Lldd_mat.colPivHouseholderQr().solve(- 1*Lld_mat);
    VectorXd Lldd_solve = VectorXd::Zero(totalnum);
    for (int ij = 0; ij < totalnum; ij++) {
        if (KeepConstant[ij] == 0) {
            int pij_ind = ij - sum(head(KeepConstant, ij));
            Lldd_solve(ij) = Lldd_solve0(pij_ind);
        }
    }
    //
    #ifdef _OPENMP
    #pragma omp parallel for schedule(dynamic) num_threads(nthreads)
    #endif
    for (int ijk = 0; ijk < totalnum; ijk++) {
        if (KeepConstant[ijk] == 0) {
            int pjk_ind = ijk - sum(head(KeepConstant, ijk));
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

//' Utility function to calculate the change to make each iteration
//'
//' \code{Calc_Change} Called to update the parameter changes, Uses log-likelihoods and control parameters, Applies newton steps and change limitations
//' @inheritParams CPP_template
//'
//' @return Updates matrices in place: parameter change matrix
//' @noRd
//'
//
void Calc_Change(const int& nthreads, const int& totalnum, const double& thres_step_max, const double& lr, const double& step_max, const vector<double>& Ll, const vector<double>& Lld, const vector<double>& Lldd, vector<double>& dbeta, const StringVector&   tform, const double& dint, const double& dslp, IntegerVector KeepConstant) {
    int kept_covs = totalnum - sum(KeepConstant);
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
        Lldd_vec[jk * kept_covs + ij] = Lldd[jk * kept_covs + ij];
        if (ij == jk) {
            Lld_vec[ij] = Lld[ij];
        } else {
            Lldd_vec[ij * kept_covs + jk] = Lldd_vec[jk * kept_covs + ij];
        }
    }
    //
    Lldd_vec.attr("dim") = Dimension(kept_covs, kept_covs);
    const Map<MatrixXd> Lldd_mat(as<Map<MatrixXd> >(Lldd_vec));
    const Map<VectorXd> Lld_mat(as<Map<VectorXd> >(Lld_vec));
    //
    VectorXd Lldd_solve0 = Lldd_mat.colPivHouseholderQr().solve(- 1*Lld_mat);
    double ll_change = (Lldd_solve0.array() * Lld_mat.array()).sum();  //  We want to make sure it is moving toward the a maximum
    for (int i = 0; i < kept_covs; i++) {  //  The predicted change should be positive, accounting for the second order taylor expansion
        for (int j = 0; j < kept_covs; j++) {
            if (i == j) {
                ll_change += 1/2 * Lldd_solve0(i) * Lldd_solve0(j) * Lldd_mat(i, j);
            } else {
                ll_change += Lldd_solve0(i) * Lldd_solve0(j) * Lldd_mat(i, j);
            }
        }
    }
    if (ll_change < 0) {
        Lldd_solve0 *= -1;  //  If it is moving in the wrong direction, turn it around
    }
    VectorXd Lldd_solve = VectorXd::Zero(totalnum);
    for (int ij = 0; ij < totalnum; ij++) {
        if (KeepConstant[ij] == 0) {
            int pij_ind = ij - sum(head(KeepConstant, ij));
            Lldd_solve(ij) = Lldd_solve0(pij_ind);
        }
    }
    //
    #ifdef _OPENMP
    #pragma omp parallel for schedule(dynamic) num_threads(nthreads)
    #endif
    for (int ijk = 0; ijk < totalnum; ijk++) {
        if (KeepConstant[ijk] == 0) {
            int pjk_ind = ijk - sum(head(KeepConstant, ijk));
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

//' Utility function to calculate the change to make each iteration, with basic model
//'
//' \code{Calc_Change_Basic} Called to update the parameter changes, Uses log-likelihoods and control parameters, Applies newton steps and change limitations
//' @inheritParams CPP_template
//'
//' @return Updates matrices in place: parameter change matrix
//' @noRd
//'
//
void Calc_Change_Basic(const int& nthreads, const int& totalnum, const double& lr, const double& step_max, const vector<double>& Ll, const vector<double>& Lld, const vector<double>& Lldd, vector<double>& dbeta, IntegerVector KeepConstant) {
    //
    int kept_covs = totalnum - sum(KeepConstant);
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
        Lldd_vec[jk * kept_covs + ij] = Lldd[jk * kept_covs + ij];
        if (ij == jk) {
            Lld_vec[ij] = Lld[ij];
        } else {
            Lldd_vec[ij * kept_covs + jk] = Lldd_vec[jk * kept_covs + ij];
        }
    }
    Lldd_vec.attr("dim") = Dimension(kept_covs, kept_covs);
    const Map<MatrixXd> Lldd_mat(as<Map<MatrixXd> >(Lldd_vec));
    const Map<VectorXd> Lld_mat(as<Map<VectorXd> >(Lld_vec));
    VectorXd Lldd_solve0 = Lldd_mat.colPivHouseholderQr().solve(- 1*Lld_mat);
    double ll_change = (Lldd_solve0.array() * Lld_mat.array()).sum();  //  We want to make sure it is moving toward the a maximum
    for (int i = 0; i < kept_covs; i++) {  //  The predicted change should be positive, accounting for the second order taylor expansion
        for (int j = 0; j < kept_covs; j++) {
            if (i == j) {
                ll_change += 1/2 * Lldd_solve0(i) * Lldd_solve0(j) * Lldd_mat(i, j);
            } else {
                ll_change += Lldd_solve0(i) * Lldd_solve0(j) * Lldd_mat(i, j);
            }
        }
    }
    if (ll_change < 0) {
        Lldd_solve0 *= -1;  //  If it is moving in the wrong direction, turn it around
    }
    VectorXd Lldd_solve = VectorXd::Zero(totalnum);
    for (int ij = 0; ij < totalnum; ij++) {
        if (KeepConstant[ij] == 0) {
            int pij_ind = ij - sum(head(KeepConstant, ij));
            Lldd_solve(ij) = Lldd_solve0(pij_ind);
        }
    }
    //
    #ifdef _OPENMP
    #pragma omp parallel for schedule(dynamic) num_threads(nthreads)
    #endif
    for (int ijk = 0; ijk < totalnum; ijk++) {
        if (KeepConstant[ijk] == 0) {
            //
            int pjk_ind = ijk - sum(head(KeepConstant, ijk));
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
            //
            if (abs(dbeta[ijk]) > step_max) {
                dbeta[ijk] = step_max * sign(dbeta[ijk]);
            }
        } else {
            dbeta[ijk] = 0;
        }
    }
    return;
}

//' Utility function to calculate the change to make each iteration, with basic model and constraints added
//'
//' \code{Calc_Change_Basic_Cons} Called to update the parameter changes, Uses log-likelihoods and control parameters, Applies newton steps and change limitations
//' @inheritParams CPP_template
//'
//' @return Updates matrices in place: parameter change matrix
//' @noRd
//'
//
void Calc_Change_Basic_Cons(const MatrixXd& Lin_Sys, const VectorXd& Lin_Res, const  VectorXd& beta_0, const int& nthreads, const int& totalnum, const double& lr, const double& step_max, const vector<double>& Ll, const vector<double>& Lld, const vector<double>& Lldd, vector<double>& dbeta, IntegerVector KeepConstant) {
    int kept_covs = totalnum - sum(KeepConstant);
    //
    VectorXd beta_1(kept_covs);
    for (int ij = 0; ij < totalnum; ij++) {
        if (KeepConstant[ij] == 0) {
            int pij_ind = ij - sum(head(KeepConstant, ij));
            beta_1(pij_ind) = beta_0(ij);
        }
    }
    VectorXd Lin_Dif = Lin_Sys * beta_1 - Lin_Res;
    //
    int total_covs = kept_covs + Lin_Sys.rows();
    //
    NumericVector Lldd_vec(total_covs*total_covs);
    NumericVector Lld_vec(total_covs);
    #ifdef _OPENMP
    #pragma omp parallel for schedule(dynamic) num_threads(nthreads)
    #endif
    for (int ijk = 0; ijk < total_covs*(total_covs + 1)/2; ijk++) {
        int ij = 0;
        int jk = ijk;
        while (jk > ij) {
            ij++;
            jk -= ij;
        }
        if (ij < kept_covs) {
            Lldd_vec[jk * total_covs + ij] = Lldd[jk * kept_covs + ij];
            if (ij == jk) {
                Lld_vec[ij] = Lld[ij];
            } else {
                Lldd_vec[ij * total_covs + jk] = Lldd_vec[jk * kept_covs + ij];
            }
        } else {
            if (jk < kept_covs) {
                Lldd_vec[jk * total_covs + ij] = Lin_Sys(ij-kept_covs, jk);
            } else {
                Lldd_vec[jk * total_covs + ij] = 0.0;
            }
            if (ij == jk) {
                Lld_vec[ij] = Lin_Dif(ij-kept_covs);
            } else {
                Lldd_vec[ij * total_covs + jk] = Lldd_vec[jk * total_covs + ij];
            }
        }
    }
    //
    //
    Lldd_vec.attr("dim") = Dimension(total_covs, total_covs);
    const Map<MatrixXd> Lldd_mat(as<Map<MatrixXd> >(Lldd_vec));
    const Map<VectorXd> Lld_mat(as<Map<VectorXd> >(Lld_vec));
    VectorXd Lldd_solve0 = Lldd_mat.colPivHouseholderQr().solve(- 1*Lld_mat);
    VectorXd Lldd_solve = VectorXd::Zero(totalnum);
    for (int ij = 0; ij < totalnum; ij++) {
        if (KeepConstant[ij] == 0) {
            int pij_ind = ij - sum(head(KeepConstant, ij));
            Lldd_solve(ij) = Lldd_solve0(pij_ind);
        }
    }
    //
    #ifdef _OPENMP
    #pragma omp parallel for schedule(dynamic) num_threads(nthreads)
    #endif
    for (int ijk = 0; ijk < totalnum; ijk++) {
        if (KeepConstant[ijk] == 0) {
            //
            int pjk_ind = ijk - sum(head(KeepConstant, ijk));
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
            //
            if (abs(dbeta[ijk]) > step_max) {
                dbeta[ijk] = step_max * sign(dbeta[ijk]);
            }
        } else {
            dbeta[ijk] = 0;
        }
    }
    return;
}

//' Utility function to calculate the change to make each iteration, applies to background terms as well
//'
//' \code{Calc_Change_Background} Called to update the parameter changes, Uses log-likelihoods and control parameters, Applies newton steps and change limitations
//' @inheritParams CPP_template
//'
//' @return Updates matrices in place: parameter change matrix
//' @noRd
//'
//
void Calc_Change_Background(const int& nthreads, const int& totalnum, const int& group_num, const double& thres_step_max, const double& lr, const double& step_max, const vector<double>& Ll, const vector<double>& Lld, const vector<double>& Lldd, vector<double>& dbeta, const StringVector& tform, const double& dint, const double& dslp, IntegerVector KeepConstant, vector<int>& strata_cond, vector<double>& LldOdds, vector<double>& LlddOdds, vector<double>& LlddOddsBeta, vector<double>& dstrata) {
    int kept_covs = totalnum - sum(KeepConstant);
    int kept_strata = group_num - reduce(strata_cond.begin(), strata_cond.end());
    int total_val = kept_covs + kept_strata;
    NumericVector Lldd_vec(total_val * total_val);
    NumericVector Lld_vec(total_val);
    #ifdef _OPENMP
    #pragma omp parallel for schedule(dynamic) num_threads(nthreads)
    #endif
    for (int ijk = 0; ijk < total_val*(total_val + 1)/2; ijk++) {
        int ij = 0;
        int jk = ijk;
        while (jk > ij) {
            ij++;
            jk -= ij;
        }
        if ((ij < kept_covs) && (jk < kept_covs)) {
            //  Both are within the model parameters
            Lldd_vec[jk * total_val + ij] = Lldd[jk * kept_covs + ij];
            if (ij == jk) {
                Lld_vec[ij] = Lld[ij];
            } else {
                Lldd_vec[ij * total_val + jk] = Lldd_vec[jk * total_val + ij];
            }
        } else if ((ij >= kept_covs) && (jk < kept_covs)) {
            //  ij is a strata parameter, jk is a model parameter
            int ij_strata = ij - kept_covs;
            Lldd_vec[jk * total_val + ij] = LlddOddsBeta[ij_strata*kept_covs + jk];
            Lldd_vec[ij * total_val + jk] = Lldd_vec[jk * total_val + ij];
        } else {
            //  Both are strata parameters
            if (ij == jk) {
                //  We only want diagonal terms
                int ij_strata = ij - kept_covs;
                Lld_vec[ij] = LldOdds[ij_strata];
                Lldd_vec[jk * total_val + ij] = LlddOdds[ij_strata];
            }
        }
    }
    //
    Lldd_vec.attr("dim") = Dimension(total_val, total_val);
    const Map<MatrixXd> Lldd_mat(as<Map<MatrixXd> >(Lldd_vec));
    const Map<VectorXd> Lld_mat(as<Map<VectorXd> >(Lld_vec));
    VectorXd Lldd_solve0 = Lldd_mat.colPivHouseholderQr().solve(- 1*Lld_mat);
    double ll_change = (Lldd_solve0.array() * Lld_mat.array()).sum();  //  We want to make sure it is moving toward the a maximum
    for (int i = 0; i < kept_covs; i++) {  //  The predicted change should be positive, accounting for the second order taylor expansion
        for (int j = 0; j < kept_covs; j++) {
            if (i == j) {
                ll_change += 1/2 * Lldd_solve0(i) * Lldd_solve0(j) * Lldd_mat(i, j);
            } else {
                ll_change += Lldd_solve0(i) * Lldd_solve0(j) * Lldd_mat(i, j);
            }
        }
    }
    if (ll_change < 0) {
        Lldd_solve0 *= -1;  //  If it is moving in the wrong direction, turn it around
    }
    VectorXd Lldd_beta_solve = VectorXd::Zero(totalnum);
    VectorXd Lldd_strata_solve = VectorXd::Zero(group_num);
    for (int ij = 0; ij < totalnum; ij++) {
        if (KeepConstant[ij] == 0) {
            int pij_ind = ij - sum(head(KeepConstant, ij));
            Lldd_beta_solve(ij) = Lldd_solve0(pij_ind);
        }
    }
    vector<int>::iterator it_end = strata_cond.begin();
    for (int ij = 0; ij < group_num; ij++) {
        if (strata_cond[ij] == 0) {
            it_end = strata_cond.begin();
            advance(it_end, ij);
            int pij_ind = ij - reduce(strata_cond.begin(), it_end) + kept_covs;
            Lldd_strata_solve(ij) = Lldd_solve0(pij_ind);
        }
    }
    //
    #ifdef _OPENMP
    #pragma omp parallel for schedule(dynamic) num_threads(nthreads)
    #endif
    for (int ijk = 0; ijk < totalnum; ijk++) {
        if (KeepConstant[ijk] == 0) {
            int pjk_ind = ijk - sum(head(KeepConstant, ijk));
            if (isnan(Lldd_beta_solve(ijk))) {
                if (Lldd[pjk_ind*kept_covs+pjk_ind] != 0) {
                    dbeta[ijk] = -lr * Lld[pjk_ind] / Lldd[pjk_ind*kept_covs+pjk_ind];
                } else {
                    dbeta[ijk] = 0;
                }
            } else {
                dbeta[ijk] = lr * Lldd_beta_solve(ijk);
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
    #ifdef _OPENMP
    #pragma omp parallel for schedule(dynamic) num_threads(nthreads)
    #endif
    for (int ijk = 0; ijk < group_num; ijk++) {
        if (strata_cond[ijk] == 0) {
            vector<int>::iterator it_end = strata_cond.begin();
            advance(it_end, ijk);
            int pjk_ind = ijk - reduce(strata_cond.begin(), it_end);
            if (isnan(Lldd_strata_solve(ijk))) {
                if (LlddOdds[pjk_ind*kept_strata+pjk_ind] != 0) {
                    dstrata[ijk] = -lr * LldOdds[pjk_ind] / LlddOdds[pjk_ind*kept_strata+pjk_ind];
                } else {
                    dstrata[ijk] = 0;
                }
            } else {
                dstrata[ijk] = lr * Lldd_strata_solve(ijk);
            }
            if (abs(dstrata[ijk]) > step_max) {
                dstrata[ijk] = step_max * sign(dstrata[ijk]);
            }
            //
        } else {
            dstrata[ijk] = 0;
        }
    }
    return;
}

