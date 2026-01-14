//  Copyright 2022 - 2025, Eric Giunta and the project collaborators, Please see main R package for license and usage details

#include <RcppEigen.h>

#include "Grouping.h"
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
#include <numeric>

#include "Colossus_types.h"


//  [[Rcpp::depends(RcppEigen)]]
//  [[Rcpp::plugins(openmp)]]

using std::endl;
using std::string;
using std::vector;
using std::transform;
using std::plus;
using std::advance;
using std::reduce;

using Eigen::Map;
using Eigen::Ref;
using Eigen::ArrayXd;
using Eigen::MatrixXd;
using Eigen::SparseMatrix;
using Eigen::VectorXd;
using Eigen::VectorXi;
using Eigen::Vector2d;
using Eigen::Vector4d;

using Rcpp::as;
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

//' Utility function to define risk groups
//'
//' \code{Make_Groups} Called to update lists of risk groups, Uses list of event times and row time/event information, Matrices store starting/stopping row indices for each group
//' @inheritParams CPP_template
//'
//' @return Updates matrices in place: Matrix of event rows for each event time, vectors of strings with rows at risk for each event time
//' @noRd
//'
void Make_Groups(const int& ntime, const Ref<const MatrixXd>& df_m, IntegerMatrix& RiskFail, vector<vector<int> >& RiskPairs, NumericVector& tu, const int& nthreads) {
    #ifdef _OPENMP
    #pragma omp parallel for schedule(dynamic) num_threads(nthreads) shared(df_m, RiskPairs, RiskFail)
    #endif
    for (int ijk = 0; ijk < ntime; ijk++) {
        double t0 = tu[ijk];
        VectorXi select_ind_all = ( ((df_m.col(0).array() < t0) && (df_m.col(1).array() >= t0)) || ( (df_m.col(0).array() == df_m.col(1).array()) &&  (df_m.col(0).array() == t0))).cast<int>();  //  indices at risk
        vector<int> indices_all;
        int th = 1;
        visit_lambda(select_ind_all,
            [&indices_all, th](double v, int i, int j) {
                if (v == th)
                    indices_all.push_back(i + 1);
            });
        vector<int> indices;  //  generates vector of (start, end) pairs for indices at risk
        for (auto it = begin(indices_all); it != end(indices_all); ++it) {
            if (indices.size() == 0) {
                indices.push_back(*it);
                indices.push_back(*it);
            } else if (indices[indices.size() - 1] + 1 < *it) {
                indices.push_back(*it);
                indices.push_back(*it);
            } else {
                indices[indices.size() - 1] = *it;
            }
        }
        RiskPairs[ijk] = indices;
        select_ind_all = ((df_m.col(2).array() == 1) && (df_m.col(1).array() == t0)).cast<int>();  //  indices with events
        indices_all.clear();
        visit_lambda(select_ind_all,
            [&indices_all, th](double v, int i, int j) {
                if (v == th)
                    indices_all.push_back(i + 1);
            });
        RiskFail(ijk, 0) = indices_all[0] - 1;  //  Due to the sorting method, there is a continuous block of event rows
        RiskFail(ijk, 1) = indices_all[indices_all.size() - 1] - 1;
    }
    return;
}

//' Utility function to define risk groups with competing risks
//'
//' \code{Make_Groups_CR} Called to update lists of risk groups, Uses list of event times and row time/event information, Matrices store starting/stopping row indices for each group, adds rows with event = 2 past the event time
//' @inheritParams CPP_template
//'
//' @return Updates matrices in place: Matrix of event rows for each event time, vectors of strings with rows at risk for each event time
//' @noRd
//'
void Make_Groups_CR(const int& ntime, const Ref<const MatrixXd>& df_m, IntegerMatrix& RiskFail, vector<vector<int> >& RiskPairs, NumericVector& tu, const VectorXd& cens_weight, const int& nthreads) {
//    vector<vector<int> > RiskPairs(ntime);
    #ifdef _OPENMP
    #pragma omp parallel for schedule(dynamic) num_threads(nthreads) shared(df_m, RiskPairs, RiskFail)
    #endif
    for (int ijk = 0; ijk < ntime; ijk++) {
        double t0 = tu[ijk];
        VectorXi select_ind_all = ((((df_m.col(0).array() < t0) || (df_m.col(0).array() == df_m.col(1).array())) && (df_m.col(1).array() >= t0)) || ((df_m.col(2).array() == 2) && (df_m.col(1).array() <= t0))).cast<int>();  //  indices at risk
        vector<int> indices_all;
        //
        int th = 1;
        //
        visit_lambda(select_ind_all,
            [&indices_all, th](double v, int i, int j) {
                if (v == th)
                    indices_all.push_back(i + 1);
            });
        //
        vector<int> indices;  //  generates vector of (start, end) pairs for indices at risk
        for (auto it = begin(indices_all); it != end(indices_all); ++it) {
            if (indices.size() == 0) {
                indices.push_back(*it);
                indices.push_back(*it);
            } else if (indices[indices.size() - 1] + 1 < *it) {
                indices.push_back(*it);
                indices.push_back(*it);
            } else {
                indices[indices.size() - 1] = *it;
            }
        }
        RiskPairs[ijk] = indices;
        //
        select_ind_all = ((df_m.col(2).array() == 1) && (df_m.col(1).array() == t0)).cast<int>();  //  indices with events
        indices_all.clear();
        visit_lambda(select_ind_all,
            [&indices_all, th](double v, int i, int j) {
                if (v == th)
                    indices_all.push_back(i + 1);
            });
        RiskFail(ijk, 0) = indices_all[0] - 1;  //  due to the sorting method, there is a continuous block of event rows
        RiskFail(ijk, 1) = indices_all[indices_all.size() - 1] - 1;
        //
    }
    return;
}

//' Utility function to define risk groups with Strata
//'
//' \code{Make_Groups_Strata} Called to update lists of risk groups, Uses list of event times and row time/event information, Matrices store starting/stopping row indices for each group
//' @inheritParams CPP_template
//'
//' @return Updates matrices in place: Matrix of event rows for each event time, vectors of strings with rows at risk for each event time
//' @noRd
//'
void Make_Groups_Strata(const int& ntime, const Ref<const MatrixXd>& df_m, IntegerMatrix& RiskFail, vector<vector<vector<int> > >& RiskPairs_Strata, NumericVector& tu, const int& nthreads, NumericVector& Strata_vals) {
    //
    #ifdef _OPENMP
    #pragma omp parallel for schedule(dynamic) num_threads(nthreads) collapse(2) shared(df_m, RiskPairs_Strata, RiskFail, Strata_vals)
    #endif
    for (int s_ij = 0; s_ij < Strata_vals.size(); s_ij++) {
        for (int ijk = 0; ijk < ntime; ijk++) {
            double t0 = tu[ijk];
            VectorXi select_ind_end = ((df_m.col(2).array() == 1) && (df_m.col(1).array() == t0) && (df_m.col(3).array() == Strata_vals[s_ij])).cast<int>();  //  indices with events
            vector<int> indices_end;
            //
            //
            int th = 1;
            visit_lambda(select_ind_end,
                [&indices_end, th](double v, int i, int j) {
                    if (v == th)
                        indices_end.push_back(i + 1);
                });
            //
            vector<int> indices;  //  generates vector of (start, end) pairs for indices at risk
            if (indices_end.size() > 0) {
                RiskFail(ijk, 2*s_ij + 0) = indices_end[0] - 1;  //  due to the sorting method, there is a continuous block of event rows
                RiskFail(ijk, 2*s_ij + 1) = indices_end[indices_end.size() - 1] - 1;
                //
                select_ind_end = (((df_m.col(0).array() < t0) || (df_m.col(0).array() == df_m.col(1).array())) && (df_m.col(1).array() >= t0) && (df_m.col(3).array() == Strata_vals[s_ij])).cast<int>();  //  indices at risk
                indices_end.clear();
                visit_lambda(select_ind_end,
                    [&indices_end, th](double v, int i, int j) {
                        if (v == th)
                            indices_end.push_back(i + 1);
                    });
                for (auto it = begin(indices_end); it != end(indices_end); ++it) {
                    if (indices.size() == 0) {
                        indices.push_back(*it);
                        indices.push_back(*it);
                    } else if (indices[indices.size() - 1] + 1 < *it) {
                        indices.push_back(*it);
                        indices.push_back(*it);
                    } else {
                        indices[indices.size() - 1] = *it;
                    }
                }
                //
                RiskPairs_Strata[ijk][s_ij] = indices;
            } else {
                RiskFail(ijk, 2*s_ij + 0) = - 1;
                RiskFail(ijk, 2*s_ij + 1) = - 1;
            }
        }
    }
    return;
}

//' Utility function to define strata groups
//'
//' \code{Make_Strata} Called to update lists of strata groups, Uses list of strata, Matrices store starting/stopping row indices for each group
//' @inheritParams CPP_template
//'
//' @return Updates matrices in place: Vectors of strings with rows at risk for each strata
//' @noRd
//'
void Make_Strata(NumericVector& Strata_vals, const Ref<const MatrixXd>& dfs, vector<vector<int> >& RiskPairs_Strata, const int& nthreads) {
    //
    #ifdef _OPENMP
    #pragma omp parallel for schedule(dynamic) num_threads(nthreads) shared(dfs, RiskPairs_Strata, Strata_vals)
    #endif
    for (int s_ij = 0; s_ij < Strata_vals.size(); s_ij++) {
        VectorXi select_ind_end = ((dfs.col(0).array() == Strata_vals[s_ij])).cast<int>();  //  indices with events
        vector<int> indices_end;
        int th = 1;
        visit_lambda(select_ind_end,
            [&indices_end, th](double v, int i, int j) {
                if (v == th)
                    indices_end.push_back(i + 1);
            });
        vector<int> indices;
        for (auto it = begin(indices_end); it != end(indices_end); ++it) {
            if (indices.size() == 0) {
                indices.push_back(*it);
                indices.push_back(*it);
            } else if (indices[indices.size() - 1] + 1 < *it) {
                indices.push_back(*it);
                indices.push_back(*it);
            } else {
                indices[indices.size() - 1] = *it;
            }
        }
        //
        RiskPairs_Strata[s_ij] = indices;
    }
    return;
}

//' Utility function to define risk groups with Strata and competing risks
//'
//' \code{Make_Groups_Strata_CR} Called to update lists of risk groups, Uses list of event times and row time/event information, Matrices store starting/stopping row indices for each group, adds competing risks
//' @inheritParams CPP_template
//'
//' @return Updates matrices in place: Matrix of event rows for each event time, vectors of strings with rows at risk for each event time
//' @noRd
//'
void Make_Groups_Strata_CR(const int& ntime, const Ref<const MatrixXd>& df_m, IntegerMatrix& RiskFail, vector<vector<vector<int> > >& RiskPairs_Strata, NumericVector& tu, const int& nthreads, NumericVector& Strata_vals, const VectorXd& cens_weight) {
    #ifdef _OPENMP
    #pragma omp parallel for schedule(dynamic) num_threads(nthreads) collapse(2) shared(df_m, RiskPairs_Strata, RiskFail, Strata_vals)
    #endif
    for (int s_ij = 0; s_ij < Strata_vals.size(); s_ij++) {
        for (int ijk = 0; ijk < ntime; ijk++) {
            double t0 = tu[ijk];
            VectorXi select_ind_end = ((df_m.col(2).array() == 1) && (df_m.col(1).array() == t0) && (df_m.col(3).array() == Strata_vals[s_ij])).cast<int>();  //  indices with events
            vector<int> indices_end;
            //
            int th = 1;
            visit_lambda(select_ind_end,
                [&indices_end, th](double v, int i, int j) {
                    if (v == th)
                        indices_end.push_back(i + 1);
                });
            //
            vector<int> indices;  //  generates vector of (start, end) pairs for indices at risk
            if (indices_end.size() > 0) {
                RiskFail(ijk, 2*s_ij + 0) = indices_end[0] - 1;  //  due to the sorting method, there is a continuous block of event rows
                RiskFail(ijk, 2*s_ij + 1) = indices_end[indices_end.size() - 1] - 1;
                //
                select_ind_end = (((((df_m.col(0).array() < t0) || (df_m.col(0).array() == df_m.col(1).array())) && (df_m.col(1).array() >= t0)) || ((df_m.col(2).array() == 2) && (df_m.col(1).array() <= t0))) && (df_m.col(3).array() == Strata_vals[s_ij])).cast<int>();  //  indices at risk
                indices_end.clear();
                visit_lambda(select_ind_end,
                    [&indices_end, th](double v, int i, int j) {
                        if (v == th)
                                indices_end.push_back(i + 1);
                    });
                //
                for (auto it = begin(indices_end); it != end(indices_end); ++it) {
                    if (indices.size() == 0) {
                        indices.push_back(*it);
                        indices.push_back(*it);
                    } else if (indices[indices.size() - 1] + 1 < *it) {
                        indices.push_back(*it);
                        indices.push_back(*it);
                    } else {
                        indices[indices.size() - 1] = *it;
                    }
                }
                //
                RiskPairs_Strata[ijk][s_ij] = indices;
            } else {
                RiskFail(ijk, 2*s_ij + 0) = - 1;
                RiskFail(ijk, 2*s_ij + 1) = - 1;
            }
        }
    }
    return;
}

//' Utility function to define matched risk groups
//'
//' \code{Make_Match} Called to update lists of risk groups, assumes the data is matched into one group and df_m only contains the event indicator
//' @inheritParams CPP_template
//'
//' @return Updates matrices in place: Matrix of event rows for each event time, vectors of strings with rows at risk for each event time, and the various recursive matrices initialized
//' @noRd
//'
void Make_Match(List& model_bool, const Ref<const MatrixXd>& df_m, IntegerMatrix& RiskFail, vector<vector<int> >& RiskPairs, vector<vector<double> >& Recur_Base, vector<vector<vector<double> > >& Recur_First, vector<vector<vector<double> > >& Recur_Second, vector<double>& strata_odds, vector<int>& strata_cond, const int& nthreads) {
    double cond_thres = model_bool["cond_thres"];
    vector<int> indices = {1, static_cast<int>(df_m.rows())};
    int nstar = static_cast<int>(df_m.rows());
    RiskPairs[0] = indices;
    //
    VectorXi select_ind_all = (df_m.col(0).array() == 1).cast<int>();  //  indices with events
    vector<int> indices_all;
    int th = 1;
    visit_lambda(select_ind_all,
        [&indices_all, th](double v, int i, int j) {
            if (v == th)
                indices_all.push_back(i + 1);
        });
    RiskFail(0, 0) = indices_all[0] - 1;  //  Due to the sorting method, there is a continuous block of event rows
    RiskFail(0, 1) = indices_all[indices_all.size() - 1] - 1;
    //
    int dj = RiskFail(0, 1) - RiskFail(0, 0) + 1;
    int m = static_cast<int>((nstar - dj + 1)*dj);
    vector<double> risk_initial(m, 0.0);
    Recur_Base[0] = risk_initial;
    if (dj > cond_thres) {
        if (nstar > dj) {
            strata_odds[0] = log(static_cast<double>(dj) / static_cast<double>(nstar - dj));
            strata_cond[0] = 0;
        }
    } else {
        if (!model_bool["single"]) {
            for (vector<vector<double> >::size_type i = 0; i< Recur_First[0].size(); i++) {
                Recur_First[0][i] = risk_initial;
            }
            for (vector<vector<double> >::size_type i = 0; i< Recur_Second[0].size(); i++) {
                Recur_Second[0][i] = risk_initial;
            }
        }
    }
    //
    return;
}

//' Utility function to define matched risk groups by strata
//'
//' \code{Make_Match_Strata} Called to update lists of risk groups, assumes the data is matched into strata and df_m only contains the strata value and then the event indicator
//' @inheritParams CPP_template
//'
//' @return Updates matrices in place: Matrix of event rows for each event time, vectors of strings with rows at risk for each event time, and the various recursive matrices initialized
//' @noRd
//'
void Make_Match_Strata(List& model_bool, const Ref<const MatrixXd>& df_m, IntegerMatrix& RiskFail, vector<vector<int> >& RiskPairs, vector<vector<double> >& Recur_Base, vector<vector<vector<double> > >& Recur_First, vector<vector<vector<double> > >& Recur_Second, vector<double>& strata_odds, vector<int>& strata_cond, const int& nthreads, NumericVector& Strata_vals) {
    double cond_thres = model_bool["cond_thres"];
    if (model_bool["single"]) {
        #ifdef _OPENMP
        #pragma omp parallel for schedule(dynamic) num_threads(nthreads) shared(df_m, RiskPairs, RiskFail, Strata_vals, Recur_Base, strata_odds)
        #endif
        for (int s_ij = 0; s_ij < Strata_vals.size(); s_ij++) {
            VectorXi select_ind_end = ((df_m.col(1).array() == 1) && (df_m.col(0).array() == Strata_vals[s_ij])).cast<int>();  //  indices with events
            vector<int> indices_end;
            //
            //
            int th = 1;
            visit_lambda(select_ind_end,
                [&indices_end, th](double v, int i, int j) {
                    if (v == th)
                        indices_end.push_back(i + 1);
                });
            //
            vector<int> indices;  //  generates vector of (start, end) pairs for indices at risk
            if (indices_end.size() > 0) {
                RiskFail(s_ij, 0) = indices_end[0] - 1;  //  due to the sorting method, there is a continuous block of event rows
                RiskFail(s_ij, 1) = indices_end[indices_end.size() - 1] - 1;
                //
                int dj = RiskFail(s_ij, 1) - RiskFail(s_ij, 0) + 1;
                select_ind_end = (df_m.col(0).array() == Strata_vals[s_ij]).cast<int>();  //  indices at risk
                indices_end.clear();
                visit_lambda(select_ind_end,
                    [&indices_end, th](double v, int i, int j) {
                        if (v == th)
                            indices_end.push_back(i + 1);
                    });
                for (auto it = begin(indices_end); it != end(indices_end); ++it) {
                    if (indices.size() == 0) {
                        indices.push_back(*it);
                        indices.push_back(*it);
                    } else if (indices[indices.size() - 1] + 1 < *it) {
                        indices.push_back(*it);
                        indices.push_back(*it);
                    } else {
                        indices[indices.size() - 1] = *it;
                    }
                }
                //
                RiskPairs[s_ij] = indices;
                int nstar = indices_end.size();
                int m = static_cast<int>((nstar - dj + 1)*dj);
                vector<double> risk_initial(m, 0.0);
                Recur_Base[s_ij] = risk_initial;
                if (dj > cond_thres) {
                    if (nstar > dj) {
                        strata_odds[s_ij] = log(static_cast<double>(dj) / static_cast<double>(nstar - dj));
                        strata_cond[s_ij] = 0;
                    }
                }
                //
            } else {
                RiskFail(s_ij, 0) = - 1;
                RiskFail(s_ij, 1) = - 1;
            }
        }
    } else {
        #ifdef _OPENMP
        #pragma omp parallel for schedule(dynamic) num_threads(nthreads) shared(df_m, RiskPairs, RiskFail, Strata_vals, Recur_Base, Recur_First, Recur_Second, strata_odds)
        #endif
        for (int s_ij = 0; s_ij < Strata_vals.size(); s_ij++) {
            VectorXi select_ind_end = ((df_m.col(1).array() == 1) && (df_m.col(0).array() == Strata_vals[s_ij])).cast<int>();  //  indices with events
            //
            vector<int> indices_end;
            int th = 1;
            visit_lambda(select_ind_end,
                [&indices_end, th](double v, int i, int j) {
                    if (v == th)
                        indices_end.push_back(i + 1);
                });
            //
            if (indices_end.size() > 0) {
                RiskFail(s_ij, 0) = indices_end[0] - 1;  //  due to the sorting method, there is a continuous block of event rows
                RiskFail(s_ij, 1) = indices_end[indices_end.size() - 1] - 1;
                //
                int dj = RiskFail(s_ij, 1) - RiskFail(s_ij, 0) + 1;
                select_ind_end = (df_m.col(0).array() == Strata_vals[s_ij]).cast<int>();  //  indices at risk
                indices_end.clear();
                visit_lambda(select_ind_end,
                    [&indices_end, th](double v, int i, int j) {
                        if (v == th)
                            indices_end.push_back(i + 1);
                    });
                vector<int> indices;  //  generates vector of (start, end) pairs for indices at risk
                for (auto it = begin(indices_end); it != end(indices_end); ++it) {
                    if (indices.size() == 0) {
                        indices.push_back(*it);
                        indices.push_back(*it);
                    } else if (indices[indices.size() - 1] + 1 < *it) {
                        indices.push_back(*it);
                        indices.push_back(*it);
                    } else {
                        indices[indices.size() - 1] = *it;
                    }
                }
                //
                RiskPairs[s_ij] = indices;
                int nstar = indices_end.size();
                int m = static_cast<int>((nstar - dj + 1)*dj);
                vector<double> risk_initial(m, 0.0);
                Recur_Base[s_ij] = risk_initial;
                if (dj > cond_thres) {
                    if (nstar > dj) {
                        strata_odds[s_ij] = log(static_cast<double>(dj) / static_cast<double>(nstar - dj));
                        strata_cond[s_ij] = 0;
                    }
                } else {
                    for (vector<vector<double> >::size_type i = 0; i< Recur_First[s_ij].size(); i++) {
                        Recur_First[s_ij][i] = risk_initial;
                    }
                    for (vector<vector<double> >::size_type i = 0; i< Recur_Second[s_ij].size(); i++) {
                        Recur_Second[s_ij][i] = risk_initial;
                    }
                }
                //
            } else {
                RiskFail(s_ij, 0) = - 1;
                RiskFail(s_ij, 1) = - 1;
            }
        }
    }
    return;
}

//' Utility function to define matched risk groups by time
//'
//' \code{Make_Match_Time} Called to update lists of risk groups, assumes the data is matched into groups by time at risk and df_m contains the interval times, and then event indicator
//' @inheritParams CPP_template
//'
//' @return Updates matrices in place: Matrix of event rows for each event time, vectors of strings with rows at risk for each event time, and the various recursive matrices initialized
//' @noRd
//'
void Make_Match_Time(List& model_bool, const int& ntime, const Ref<const MatrixXd>& df_m, IntegerMatrix& RiskFail, vector<vector<int> >& RiskPairs, vector<vector<double> >& Recur_Base, vector<vector<vector<double> > >& Recur_First, vector<vector<vector<double> > >& Recur_Second, vector<double>& strata_odds, vector<int>& strata_cond, const int& nthreads, NumericVector& tu) {
    double cond_thres = model_bool["cond_thres"];
    if (model_bool["single"]) {
        #ifdef _OPENMP
        #pragma omp parallel for schedule(dynamic) num_threads(nthreads) shared(df_m, RiskPairs, RiskFail, Recur_Base, strata_odds)
        #endif
        for (int ijk = 0; ijk < ntime; ijk++) {
            double t0 = tu[ijk];
            VectorXi select_ind_all = ( ((df_m.col(0).array() < t0) && (df_m.col(1).array() >= t0)) || ( (df_m.col(0).array() == df_m.col(1).array()) &&  (df_m.col(0).array() == t0))).cast<int>();  //  indices at risk
            vector<int> indices_all;
            int th = 1;
            visit_lambda(select_ind_all,
                [&indices_all, th](double v, int i, int j) {
                    if (v == th)
                        indices_all.push_back(i + 1);
                });
            vector<int> indices;  //  generates vector of (start, end) pairs for indices at risk
            for (auto it = begin(indices_all); it != end(indices_all); ++it) {
                if (indices.size() == 0) {
                    indices.push_back(*it);
                    indices.push_back(*it);
                } else if (indices[indices.size() - 1] + 1 < *it) {
                    indices.push_back(*it);
                    indices.push_back(*it);
                } else {
                    indices[indices.size() - 1] = *it;
                }
            }
            RiskPairs[ijk] = indices;
            int nstar = indices_all.size();
            select_ind_all = ((df_m.col(2).array() == 1) && (df_m.col(1).array() == t0)).cast<int>();  //  indices with events
            indices_all.clear();
            visit_lambda(select_ind_all,
                [&indices_all, th](double v, int i, int j) {
                    if (v == th)
                        indices_all.push_back(i + 1);
                });
            RiskFail(ijk, 0) = indices_all[0] - 1;  //  Due to the sorting method, there is a continuous block of event rows
            RiskFail(ijk, 1) = indices_all[indices_all.size() - 1] - 1;
            //
            int dj = RiskFail(ijk, 1) - RiskFail(ijk, 0) + 1;
            int m = static_cast<int>((nstar - dj + 1)*dj);
            vector<double> risk_initial(m, 0.0);
            Recur_Base[ijk] = risk_initial;
            if (dj > cond_thres) {
                if (nstar > dj) {
                    strata_odds[ijk] = log(static_cast<double>(dj) / static_cast<double>(nstar - dj));
                    strata_cond[ijk] = 0;
                }
            }
        }
    } else {
        #ifdef _OPENMP
        #pragma omp parallel for schedule(dynamic) num_threads(nthreads) shared(df_m, RiskPairs, RiskFail, Recur_Base, Recur_First, Recur_Second, strata_odds)
        #endif
        for (int ijk = 0; ijk < ntime; ijk++) {
            double t0 = tu[ijk];
            VectorXi select_ind_all = ( ((df_m.col(0).array() < t0) && (df_m.col(1).array() >= t0)) || ( (df_m.col(0).array() == df_m.col(1).array()) &&  (df_m.col(0).array() == t0))).cast<int>();  //  indices at risk
            vector<int> indices_all;
            int th = 1;
            visit_lambda(select_ind_all,
                [&indices_all, th](double v, int i, int j) {
                    if (v == th)
                        indices_all.push_back(i + 1);
                });
            vector<int> indices;  //  generates vector of (start, end) pairs for indices at risk
            for (auto it = begin(indices_all); it != end(indices_all); ++it) {
                if (indices.size() == 0) {
                    indices.push_back(*it);
                    indices.push_back(*it);
                } else if (indices[indices.size() - 1] + 1 < *it) {
                    indices.push_back(*it);
                    indices.push_back(*it);
                } else {
                    indices[indices.size() - 1] = *it;
                }
            }
            RiskPairs[ijk] = indices;
            int nstar = indices_all.size();
            select_ind_all = ((df_m.col(2).array() == 1) && (df_m.col(1).array() == t0)).cast<int>();  //  indices with events
            indices_all.clear();
            visit_lambda(select_ind_all,
                [&indices_all, th](double v, int i, int j) {
                    if (v == th)
                        indices_all.push_back(i + 1);
                });
            RiskFail(ijk, 0) = indices_all[0] - 1;  //  Due to the sorting method, there is a continuous block of event rows
            RiskFail(ijk, 1) = indices_all[indices_all.size() - 1] - 1;
            //
            int dj = RiskFail(ijk, 1) - RiskFail(ijk, 0) + 1;
            int m = static_cast<int>((nstar - dj + 1)*dj);
            vector<double> risk_initial(m, 0.0);
            Recur_Base[ijk] = risk_initial;
            if (dj > cond_thres) {
                if (nstar > dj) {
                    strata_odds[ijk] = log(static_cast<double>(dj) / static_cast<double>(nstar - dj));
                    strata_cond[ijk] = 0;
                }
            } else {
                for (vector<vector<double> >::size_type i = 0; i< Recur_First[ijk].size(); i++) {
                    Recur_First[ijk][i] = risk_initial;
                }
                for (vector<vector<double> >::size_type i = 0; i< Recur_Second[ijk].size(); i++) {
                    Recur_Second[ijk][i] = risk_initial;
                }
            }
        }
    }
    return;
}

//' Utility function to define matched risk groups by time and strata
//'
//' \code{Make_Match_Time_Strata} Called to update lists of risk groups, assumes the data is matched into groups by time at risk as well as strata and df_m contains the interval times, the strata value, and then event indicator
//' @inheritParams CPP_template
//'
//' @return Updates matrices in place: Matrix of event rows for each event time, vectors of strings with rows at risk for each event time, and the various recursive matrices initialized
//' @noRd
//'
void Make_Match_Time_Strata(List& model_bool, const int& ntime, const Ref<const MatrixXd>& df_m, IntegerMatrix& RiskFail, vector<vector<int> >& RiskPairs, vector<vector<double> >& Recur_Base, vector<vector<vector<double> > >& Recur_First, vector<vector<vector<double> > >& Recur_Second, vector<double>& strata_odds, vector<int>& strata_cond, const int& nthreads, NumericVector& tu, NumericVector& Strata_vals) {
    double cond_thres = model_bool["cond_thres"];
    if (model_bool["single"]) {
        #ifdef _OPENMP
        #pragma omp parallel for schedule(dynamic) num_threads(nthreads) collapse(2) shared(df_m, RiskPairs, RiskFail, Strata_vals, Recur_Base, strata_odds)
        #endif
        for (int s_ij = 0; s_ij < Strata_vals.size(); s_ij++) {
            for (int ijk = 0; ijk < ntime; ijk++) {
                double t0 = tu[ijk];
                VectorXi select_ind_end = ((df_m.col(3).array() == 1) && (df_m.col(1).array() == t0) && (df_m.col(2).array() == Strata_vals[s_ij])).cast<int>();  //  indices with events
                //
                vector<int> indices_end;
                int th = 1;
                visit_lambda(select_ind_end,
                    [&indices_end, th](double v, int i, int j) {
                        if (v == th)
                            indices_end.push_back(i + 1);
                    });
                //
                if (indices_end.size() > 0) {
                    RiskFail(s_ij*ntime+ijk, 0) = indices_end[0] - 1;  //  due to the sorting method, there is a continuous block of event rows
                    RiskFail(s_ij*ntime+ijk, 1) = indices_end[indices_end.size() - 1] - 1;
                    //
                    int dj = RiskFail(s_ij*ntime+ijk, 1) - RiskFail(s_ij*ntime+ijk, 0) + 1;
                    //
                    select_ind_end = (((df_m.col(0).array() < t0) || (df_m.col(0).array() == df_m.col(1).array())) && (df_m.col(1).array() >= t0) && (df_m.col(2).array() == Strata_vals[s_ij])).cast<int>();  //  indices at risk
                    indices_end.clear();
                    visit_lambda(select_ind_end,
                        [&indices_end, th](double v, int i, int j) {
                            if (v == th)
                                indices_end.push_back(i + 1);
                        });
                    vector<int> indices;  //  generates vector of (start, end) pairs for indices at risk
                    for (auto it = begin(indices_end); it != end(indices_end); ++it) {
                        if (indices.size() == 0) {
                            indices.push_back(*it);
                            indices.push_back(*it);
                        } else if (indices[indices.size() - 1] + 1 < *it) {
                            indices.push_back(*it);
                            indices.push_back(*it);
                        } else {
                            indices[indices.size() - 1] = *it;
                        }
                    }
                    //
                    RiskPairs[s_ij*ntime+ijk] = indices;
                    int nstar = indices_end.size();
                    int m = static_cast<int>((nstar - dj + 1)*dj);
                    vector<double> risk_initial(m, 0.0);
                    Recur_Base[s_ij*ntime+ijk] = risk_initial;
                    if (dj > cond_thres) {
                        if (nstar > dj) {
                            strata_odds[s_ij*ntime+ijk] = log(static_cast<double>(dj) / static_cast<double>(nstar - dj));
                            strata_cond[s_ij*ntime+ijk] = 0;
                        }
                    }
                } else {
                    RiskFail(s_ij*ntime+ijk, 0) = - 1;
                    RiskFail(s_ij*ntime+ijk, 1) = - 1;
                }
            }
        }
    } else {
        #ifdef _OPENMP
        #pragma omp parallel for schedule(dynamic) num_threads(nthreads) collapse(2) shared(df_m, RiskPairs, RiskFail, Strata_vals, Recur_Base, Recur_First, Recur_Second, strata_odds)
        #endif
        for (int s_ij = 0; s_ij < Strata_vals.size(); s_ij++) {
            for (int ijk = 0; ijk < ntime; ijk++) {
                double t0 = tu[ijk];
                VectorXi select_ind_end = ((df_m.col(3).array() == 1) && (df_m.col(1).array() == t0) && (df_m.col(2).array() == Strata_vals[s_ij])).cast<int>();  //  indices with events
                vector<int> indices_end;
                //
                //
                int th = 1;
                visit_lambda(select_ind_end,
                    [&indices_end, th](double v, int i, int j) {
                        if (v == th)
                            indices_end.push_back(i + 1);
                    });
                //
                vector<int> indices;  //  generates vector of (start, end) pairs for indices at risk
                if (indices_end.size() > 0) {
                    RiskFail(s_ij*ntime+ijk, 0) = indices_end[0] - 1;  //  due to the sorting method, there is a continuous block of event rows
                    RiskFail(s_ij*ntime+ijk, 1) = indices_end[indices_end.size() - 1] - 1;
                    //
                    int dj = RiskFail(s_ij*ntime+ijk, 1) - RiskFail(s_ij*ntime+ijk, 0) + 1;
                    //
                    select_ind_end = (((df_m.col(0).array() < t0) || (df_m.col(0).array() == df_m.col(1).array())) && (df_m.col(1).array() >= t0) && (df_m.col(2).array() == Strata_vals[s_ij])).cast<int>();  //  indices at risk
                    indices_end.clear();
                    visit_lambda(select_ind_end,
                        [&indices_end, th](double v, int i, int j) {
                            if (v == th)
                                indices_end.push_back(i + 1);
                        });
                    for (auto it = begin(indices_end); it != end(indices_end); ++it) {
                        if (indices.size() == 0) {
                            indices.push_back(*it);
                            indices.push_back(*it);
                        } else if (indices[indices.size() - 1] + 1 < *it) {
                            indices.push_back(*it);
                            indices.push_back(*it);
                        } else {
                            indices[indices.size() - 1] = *it;
                        }
                    }
                    //
                    RiskPairs[s_ij*ntime+ijk] = indices;
                    int nstar = indices_end.size();
                    int m = static_cast<int>((nstar - dj + 1)*dj);
                    vector<double> risk_initial(m, 0.0);
                    Recur_Base[s_ij*ntime+ijk] = risk_initial;
                    if (dj > cond_thres) {
                        if (nstar > dj) {
                            strata_odds[s_ij*ntime+ijk] = log(static_cast<double>(dj) / static_cast<double>(nstar - dj));
                            strata_cond[s_ij*ntime+ijk] = 0;
                        }
                    } else {
                        for (vector<vector<double> >::size_type i = 0; i< Recur_First[s_ij*ntime+ijk].size(); i++) {
                            Recur_First[s_ij*ntime+ijk][i] = risk_initial;
                        }
                        for (vector<vector<double> >::size_type i = 0; i< Recur_Second[s_ij*ntime+ijk].size(); i++) {
                            Recur_Second[s_ij*ntime+ijk][i] = risk_initial;
                        }
                    }
                } else {
                    RiskFail(s_ij*ntime+ijk, 0) = - 1;
                    RiskFail(s_ij*ntime+ijk, 1) = - 1;
                }
            }
        }
    }
    return;
}
