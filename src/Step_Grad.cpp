//  Copyright 2022 - 2025, Eric Giunta and the project collaborators, Please see main R package for license and usage details

#include <RcppEigen.h>

#include "Step_Grad.h"
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

//' Utility function to calculate the change to make each iteration with gradient step
//'
//' \code{Calc_Change_Gradient} Called to update the parameter changes, Uses log-likelihoods and control parameters, Applies gradient normalization and change limitations
//' @inheritParams CPP_template
//'
//' @return Updates matrices in place: parameter change matrix
//' @noRd
//'
//
void Calc_Change_Gradient(const int& nthreads, List& model_bool, const int& totalnum, List& optim_para, int& iteration, const double& step_max, const vector<double>& Lld, NumericVector& m_g_store, NumericVector& v_beta_store, vector<double>& dbeta, IntegerVector KeepConstant) {
    int kept_covs = totalnum - sum(KeepConstant);
    NumericVector Lld_vec(kept_covs);
    for (int ij = 0; ij < kept_covs; ij++) {
        Lld_vec[ij] = Lld[ij];
    }
    //
    //  Written for the sake of preparing what variables will be needed
    bool momentum_bool = model_bool["momentum"];  //  assumed I will define booleans to pick which one is used
    bool adadelta_bool = model_bool["adadelta"];
    bool adam_bool     = model_bool["adam"];
    //  grouped the necessary variables
    double lr = optim_para["lr"];
    double decay1 = optim_para["momentum_decay"];
    double decay2 = optim_para["learning_decay"];
    double epsilon_momentum = optim_para["epsilon_decay"];
    //  required vectors for storage
    if (momentum_bool) {
        //
        for (int ijk = 0; ijk < totalnum; ijk++) {
            if (KeepConstant[ijk] == 0) {
                int pjk_ind = ijk - sum(head(KeepConstant, ijk));
                dbeta[ijk] = decay1 * m_g_store[pjk_ind] + lr * Lld_vec[pjk_ind];
                m_g_store[pjk_ind] = dbeta[ijk];
                if (abs(dbeta[ijk]) > step_max) {
                    dbeta[ijk] = step_max * sign(dbeta[ijk]);
                }
            } else {
                dbeta[ijk] = 0;
            }
        }
    } else if (adadelta_bool) {
        //
        for (int ijk = 0; ijk < totalnum; ijk++) {
            if (KeepConstant[ijk] == 0) {
                int pjk_ind = ijk - sum(head(KeepConstant, ijk));
                m_g_store[pjk_ind] = decay1 * m_g_store[pjk_ind] + (1-decay1)*pow(Lld_vec[pjk_ind], 2);
                v_beta_store[pjk_ind] = decay1 * v_beta_store[pjk_ind] + (1-decay1)*pow(dbeta[ijk], 2);
                dbeta[ijk] = pow((v_beta_store[pjk_ind] + epsilon_momentum) / (m_g_store[pjk_ind] +epsilon_momentum), 0.5)* Lld_vec[pjk_ind];
                if (abs(dbeta[ijk]) > step_max) {
                    dbeta[ijk] = step_max * sign(dbeta[ijk]);
                }
            } else {
                dbeta[ijk] = 0;
            }
        }
    } else if (adam_bool) {
        //
        for (int ijk = 0; ijk < totalnum; ijk++) {
            if (KeepConstant[ijk] == 0) {
                int pjk_ind = ijk - sum(head(KeepConstant, ijk));
                m_g_store[pjk_ind] = decay1 * m_g_store[pjk_ind] + (1-decay1)*Lld_vec[pjk_ind];
                v_beta_store[pjk_ind] = decay2 * v_beta_store[pjk_ind] + (1-decay2)*pow(Lld_vec[pjk_ind], 2);
                //
                double m_t_bias = m_g_store[pjk_ind] / (1 - pow(decay1, iteration));
                double v_t_bias = v_beta_store[pjk_ind] / (1 - pow(decay2, iteration));
                //
                dbeta[ijk] = lr / (pow(v_t_bias, 0.5)+epsilon_momentum) * m_t_bias;
                if (abs(dbeta[ijk]) > step_max) {
                    dbeta[ijk] = step_max * sign(dbeta[ijk]);
                }
            } else {
                dbeta[ijk] = 0;
            }
        }
    } else {
        //
        #ifdef _OPENMP
        #pragma omp parallel for schedule(dynamic) num_threads(nthreads)
        #endif
        for (int ijk = 0; ijk < totalnum; ijk++) {
            if (KeepConstant[ijk] == 0) {
                int pjk_ind = ijk - sum(head(KeepConstant, ijk));
                dbeta[ijk] = lr * Lld_vec[pjk_ind];
                if (abs(dbeta[ijk]) > step_max) {
                    dbeta[ijk] = step_max * sign(dbeta[ijk]);
                }
            } else {
                dbeta[ijk] = 0;
            }
        }
    }
    return;
}

//' Utility function to calculate the change to make each iteration with gradient step following a linear constraint
//'
//' \code{Calc_Change_Gradient_Cons} Called to update the parameter changes, Uses log-likelihoods and control parameters, Applies gradient normalization and change limitations
//' @inheritParams CPP_template
//'
//' @return Updates matrices in place: parameter change matrix
//' @noRd
//'
//
void Calc_Change_Gradient_Cons(const MatrixXd& Lin_Sys, const VectorXd& Lin_Res, const int& nthreads, List& model_bool, const int& totalnum, List& optim_para, int& iteration, const double& step_max, const vector<double>& Ll, const vector<double>& Lld, NumericVector& m_g_store, NumericVector& v_beta_store, const VectorXd& beta_0, vector<double>& dbeta, IntegerVector KeepConstant) {
    int kept_covs = totalnum - sum(KeepConstant);
    // We need to get the penalized function
    double penalty_weight = optim_para["penalty_weight"];
    penalty_weight *= pow(1.1, iteration);  // scaling the penalty weight up with additional iterations
    VectorXd pred_sys = Lin_Sys * beta_0.matrix();
    VectorXd pred_delta = pred_sys.array() - Lin_Res.array();
//    double S = Ll[0] - penalty_weight * (pred_delta.array().pow(2).sum());
    VectorXd dbeta_vec = VectorXd::Zero(totalnum);
    NumericVector Lld_vec(kept_covs);
    double score_adjust = 0;
    for (int ij = 0; ij < kept_covs; ij++) {
        score_adjust = 2 * penalty_weight * (pred_delta.array() * Lin_Sys.col(ij).array()).sum();
        Lld_vec[ij] = Lld[ij] - score_adjust;
    }
    //
    //  Written for the sake of preparing what variables will be needed
    bool momentum_bool = model_bool["momentum"];  //  assumed I will define booleans to pick which one is used
    bool adadelta_bool = model_bool["adadelta"];
    bool adam_bool     = model_bool["adam"];
    //  grouped the necessary variables
    double lr = optim_para["lr"];
    double decay1 = optim_para["momentum_decay"];
    double decay2 = optim_para["learning_decay"];
    double epsilon_momentum = optim_para["epsilon_decay"];
    //  required vectors for storage
    if (momentum_bool) {
        //
        for (int ijk = 0; ijk < totalnum; ijk++) {
            if (KeepConstant[ijk] == 0) {
                int pjk_ind = ijk - sum(head(KeepConstant, ijk));
                dbeta[ijk] = decay1 * m_g_store[pjk_ind] + lr * Lld_vec[pjk_ind];
                m_g_store[pjk_ind] = dbeta[ijk];
                if (abs(dbeta[ijk]) > step_max) {
                    dbeta[ijk] = step_max * sign(dbeta[ijk]);
                }
            } else {
                dbeta[ijk] = 0;
            }
        }
    } else if (adadelta_bool) {
        //
        for (int ijk = 0; ijk < totalnum; ijk++) {
            if (KeepConstant[ijk] == 0) {
                int pjk_ind = ijk - sum(head(KeepConstant, ijk));
                m_g_store[pjk_ind] = decay1 * m_g_store[pjk_ind] + (1-decay1)*pow(Lld_vec[pjk_ind], 2);
                v_beta_store[pjk_ind] = decay1 * v_beta_store[pjk_ind] + (1-decay1)*pow(dbeta[ijk], 2);
                dbeta[ijk] = pow((v_beta_store[pjk_ind] + epsilon_momentum) / (m_g_store[pjk_ind] +epsilon_momentum), 0.5)* Lld_vec[pjk_ind];
                if (abs(dbeta[ijk]) > step_max) {
                    dbeta[ijk] = step_max * sign(dbeta[ijk]);
                }
            } else {
                dbeta[ijk] = 0;
            }
        }
    } else if (adam_bool) {
        //
        for (int ijk = 0; ijk < totalnum; ijk++) {
            if (KeepConstant[ijk] == 0) {
                int pjk_ind = ijk - sum(head(KeepConstant, ijk));
                m_g_store[pjk_ind] = decay1 * m_g_store[pjk_ind] + (1-decay1)*Lld_vec[pjk_ind];
                v_beta_store[pjk_ind] = decay2 * v_beta_store[pjk_ind] + (1-decay2)*pow(Lld_vec[pjk_ind], 2);
                //
                double m_t_bias = m_g_store[pjk_ind] / (1 - pow(decay1, iteration));
                double v_t_bias = v_beta_store[pjk_ind] / (1 - pow(decay2, iteration));
                //
                dbeta[ijk] = lr / (pow(v_t_bias, 0.5)+epsilon_momentum) * m_t_bias;
                if (abs(dbeta[ijk]) > step_max) {
                    dbeta[ijk] = step_max * sign(dbeta[ijk]);
                }
            } else {
                dbeta[ijk] = 0;
            }
        }
    } else {
        //
        #ifdef _OPENMP
        #pragma omp parallel for schedule(dynamic) num_threads(nthreads)
        #endif
        for (int ijk = 0; ijk < totalnum; ijk++) {
            if (KeepConstant[ijk] == 0) {
                int pjk_ind = ijk - sum(head(KeepConstant, ijk));
                dbeta[ijk] = lr * Lld_vec[pjk_ind];
                if (abs(dbeta[ijk]) > step_max) {
                    dbeta[ijk] = step_max * sign(dbeta[ijk]);
                }
            } else {
                dbeta[ijk] = 0;
            }
        }
    }
    for (int ijk = 0; ijk < totalnum; ijk++) {
        dbeta_vec[ijk] = dbeta[ijk];
    }
    pred_sys = Lin_Sys * (beta_0.array() + dbeta_vec.array()).matrix();
    pred_delta = Lin_Res.array() - pred_sys.array();
    VectorXd dbeta_adjust = Lin_Sys.fullPivHouseholderQr().solve(pred_delta);
    for (int ijk = 0; ijk < totalnum; ijk++) {
        dbeta[ijk] += dbeta_adjust[ijk];
        if (KeepConstant[ijk] == 0) {
            if (abs(dbeta[ijk]) > step_max) {
                dbeta[ijk] = step_max * sign(dbeta[ijk]);
            }
        } else {
            dbeta[ijk] = 0;
        }
    }
    return;
}

//' Utility function to calculate the change to make each iteration with gradient step and background terms
//'
//' \code{Calc_Change_Background_Gradient} Called to update the parameter changes, Uses log-likelihoods and control parameters, Applies gradient normalization and change limitations
//' @inheritParams CPP_template
//'
//' @return Updates matrices in place: parameter change matrix
//' @noRd
//'
//
void Calc_Change_Background_Gradient(const int& nthreads, List& model_bool, const int& totalnum, const int& group_num, List& optim_para, int& iteration, const double& step_max, const vector<double>& Lld, NumericVector& m_g_store, NumericVector& v_beta_store, vector<double>& dbeta, IntegerVector KeepConstant, vector<int>& strata_cond, vector<double>& LldOdds, vector<double>& dstrata) {
    int kept_covs = totalnum - sum(KeepConstant);
    int kept_strata = group_num - reduce(strata_cond.begin(), strata_cond.end());
    int total_val = kept_covs + kept_strata;
    NumericVector Lld_vec(total_val);
    for (int ij = 0; ij < kept_covs; ij++) {
        Lld_vec[ij] = Lld[ij];
    }
    for (int ij = 0; ij < kept_strata; ij++) {
        Lld_vec[ij+kept_covs] = LldOdds[ij];
    }
    //
    //  Written for the sake of preparing what variables will be needed
    bool momentum_bool = model_bool["momentum"];  //  assumed I will define booleans to pick which one is used
    bool adadelta_bool = model_bool["adadelta"];
    bool adam_bool     = model_bool["adam"];
    //  grouped the necessary variables
    double lr = optim_para["lr"];
    double decay1 = optim_para["momentum_decay"];
    double decay2 = optim_para["learning_decay"];
    double epsilon_momentum = optim_para["epsilon_decay"];
    //  required vectors for storage
    if (momentum_bool) {
        //
        for (int ijk = 0; ijk < totalnum; ijk++) {
            if (KeepConstant[ijk] == 0) {
                int pjk_ind = ijk - sum(head(KeepConstant, ijk));
                dbeta[ijk] = decay1 * m_g_store[pjk_ind] + lr * Lld_vec[pjk_ind];
                m_g_store[pjk_ind] = dbeta[ijk];
                if (abs(dbeta[ijk]) > step_max) {
                    dbeta[ijk] = step_max * sign(dbeta[ijk]);
                }
            } else {
                dbeta[ijk] = 0;
            }
        }
        for (int ijk = 0; ijk < group_num; ijk++) {
            if (strata_cond[ijk] == 0) {
                vector<int>::iterator it_end = strata_cond.begin();
                advance(it_end, ijk);
                int pjk_ind = ijk - reduce(strata_cond.begin(), it_end) + kept_covs;
                dstrata[ijk] = decay1 * m_g_store[pjk_ind] + lr * Lld_vec[pjk_ind];
                m_g_store[pjk_ind] = dstrata[ijk];
                if (abs(dstrata[ijk]) > step_max) {
                    dstrata[ijk] = step_max * sign(dstrata[ijk]);
                }
            } else {
                dstrata[ijk] = 0;
            }
        }
    } else if (adadelta_bool) {
        //
        for (int ijk = 0; ijk < totalnum; ijk++) {
            if (KeepConstant[ijk] == 0) {
                int pjk_ind = ijk - sum(head(KeepConstant, ijk));
                m_g_store[pjk_ind] = decay1 * m_g_store[pjk_ind] + (1-decay1)*pow(Lld_vec[pjk_ind], 2);
                v_beta_store[pjk_ind] = decay1 * v_beta_store[pjk_ind] + (1-decay1)*pow(dbeta[ijk], 2);
                dbeta[ijk] = pow((v_beta_store[pjk_ind] + epsilon_momentum) / (m_g_store[pjk_ind] + epsilon_momentum), 0.5)* Lld_vec[pjk_ind];
                if (abs(dbeta[ijk]) > step_max) {
                    dbeta[ijk] = step_max * sign(dbeta[ijk]);
                }
            } else {
                dbeta[ijk] = 0;
            }
        }
        for (int ijk = 0; ijk < group_num; ijk++) {
            if (strata_cond[ijk] == 0) {
                vector<int>::iterator it_end = strata_cond.begin();
                advance(it_end, ijk);
                int pjk_ind = ijk - reduce(strata_cond.begin(), it_end) + kept_covs;
                m_g_store[pjk_ind] = decay1 * m_g_store[pjk_ind] + (1-decay1)*pow(Lld_vec[pjk_ind], 2);
                v_beta_store[pjk_ind] = decay1 * v_beta_store[pjk_ind] + (1-decay1)*pow(dstrata[ijk], 2);
                dstrata[ijk] = pow((v_beta_store[pjk_ind] + epsilon_momentum) / (m_g_store[pjk_ind] + epsilon_momentum), 0.5)* Lld_vec[pjk_ind];
                if (abs(dstrata[ijk]) > step_max) {
                    dstrata[ijk] = step_max * sign(dstrata[ijk]);
                }
            } else {
                dstrata[ijk] = 0;
            }
        }
    } else if (adam_bool) {
        //
        for (int ijk = 0; ijk < totalnum; ijk++) {
            if (KeepConstant[ijk] == 0) {
                int pjk_ind = ijk - sum(head(KeepConstant, ijk));
                m_g_store[pjk_ind] = decay1 * m_g_store[pjk_ind] + (1-decay1)*Lld_vec[pjk_ind];
                v_beta_store[pjk_ind] = decay2 * v_beta_store[pjk_ind] + (1-decay2)*pow(Lld_vec[pjk_ind], 2);
                //
                double m_t_bias = m_g_store[pjk_ind] / (1 - pow(decay1, iteration));
                double v_t_bias = v_beta_store[pjk_ind] / (1 - pow(decay2, iteration));
                //
                dbeta[ijk] = lr / (pow(v_t_bias, 0.5)+epsilon_momentum) * m_t_bias;
                if (abs(dbeta[ijk]) > step_max) {
                    dbeta[ijk] = step_max * sign(dbeta[ijk]);
                }
            } else {
                dbeta[ijk] = 0;
            }
        }
        for (int ijk = 0; ijk < group_num; ijk++) {
            if (strata_cond[ijk] == 0) {
                vector<int>::iterator it_end = strata_cond.begin();
                advance(it_end, ijk);
                int pjk_ind = ijk - reduce(strata_cond.begin(), it_end) + kept_covs;
                m_g_store[pjk_ind] = decay1 * m_g_store[pjk_ind] + (1-decay1)*Lld_vec[pjk_ind];
                v_beta_store[pjk_ind] = decay2 * v_beta_store[pjk_ind] + (1-decay2)*pow(Lld_vec[pjk_ind], 2);
                //
                double m_t_bias = m_g_store[pjk_ind] / (1 - pow(decay1, iteration));
                double v_t_bias = v_beta_store[pjk_ind] / (1 - pow(decay2, iteration));
                //
                dstrata[ijk] = lr / (pow(v_t_bias, 0.5)+epsilon_momentum) * m_t_bias;
                if (abs(dstrata[ijk]) > step_max) {
                    dstrata[ijk] = step_max * sign(dstrata[ijk]);
                }
            } else {
                dstrata[ijk] = 0;
            }
        }
    } else {
        //
        for (int ijk = 0; ijk < totalnum; ijk++) {
            if (KeepConstant[ijk] == 0) {
                int pjk_ind = ijk - sum(head(KeepConstant, ijk));
                dbeta[ijk] = lr * Lld_vec[pjk_ind];
                if (abs(dbeta[ijk]) > step_max) {
                    dbeta[ijk] = step_max * sign(dbeta[ijk]);
                }
            } else {
                dbeta[ijk] = 0;
            }
        }
        for (int ijk = 0; ijk < group_num; ijk++) {
            if (strata_cond[ijk] == 0) {
                vector<int>::iterator it_end = strata_cond.begin();
                advance(it_end, ijk);
                int pjk_ind = ijk - reduce(strata_cond.begin(), it_end) + kept_covs;
                dstrata[ijk] = lr * Lld_vec[pjk_ind];
                if (abs(dstrata[ijk]) > step_max) {
                    dstrata[ijk] = step_max * sign(dstrata[ijk]);
                }
            } else {
                dstrata[ijk] = 0;
            }
        }
    }
    return;
}
