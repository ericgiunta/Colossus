//  Copyright 2022 - 2025, Eric Giunta and the project collaborators, Please see main R package for license and usage details

#include <RcppEigen.h>

#include "Calc_Repeated.h"
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

#include "Subterms_Risk.h"
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

void removeRow(MatrixXd& matrix, unsigned int rowToRemove) {
    unsigned int numRows = matrix.rows() - 1;
    unsigned int numCols = matrix.cols();

    if (rowToRemove < numRows)
        matrix.block(rowToRemove, 0, numRows-rowToRemove, numCols) = matrix.block(rowToRemove + 1, 0, numRows-rowToRemove, numCols);

    matrix.conservativeResize(numRows, numCols);
}

void removeColumn(MatrixXd& matrix, unsigned int colToRemove) {
    unsigned int numRows = matrix.rows();
    unsigned int numCols = matrix.cols() - 1;

    if (colToRemove < numCols)
        matrix.block(0, colToRemove, numRows, numCols-colToRemove) = matrix.block(0, colToRemove + 1, numRows, numCols-colToRemove);

    matrix.conservativeResize(numRows, numCols);
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

//' Utility function to calculate repeated values used in Cox Log-Likelihood calculation
//'
//' \code{Calculate_Sides} Called to update repeated sum calculations, Uses list of event rows and risk matrices, Performs calculation of sums of risk in each group
//' @inheritParams CPP_template
//'
//' @return Updates matrices in place: risk storage matrices
//' @noRd
//'
void Calculate_Sides(List& model_bool, const IntegerMatrix& RiskFail, const vector<vector<int> >& RiskPairs, const int& totalnum, const int& ntime, const MatrixXd& R, const MatrixXd& Rd, const MatrixXd& Rdd, MatrixXd& Rls1, MatrixXd& Rls2, MatrixXd& Rls3, MatrixXd& Lls1, MatrixXd& Lls2, MatrixXd& Lls3, const int& nthreads, const IntegerVector& KeepConstant) {
    int reqrdnum = totalnum - sum(KeepConstant);
    //
    #ifdef _OPENMP
    #pragma omp parallel for schedule(dynamic) num_threads(nthreads) shared(RiskPairs, RiskFail, R, Rls1, Lls1)
    #endif
    for (int j = 0; j < ntime; j++) {
        double Rs1 = 0;
        //
        //
        vector<int> InGroup = RiskPairs[j];
        //  now has the grouping pairs
        int dj = RiskFail(j, 1)-RiskFail(j, 0) + 1;
        for (vector<double>::size_type i = 0; i < InGroup.size() - 1; i = i+2) {
            Rs1 += R.block(InGroup[i] - 1, 0, InGroup[i + 1]-InGroup[i] + 1, 1).sum();
        }  //  precalculates the sums of risk groups
        //  only assigns values once
        Rls1(j, 0) = Rs1;
        Lls1(j, 0) = R.block(RiskFail(j, 0), 0, dj, 1).sum();
    }
    if (!model_bool["single"]) {
        //
        #ifdef _OPENMP
        #pragma omp parallel for schedule(dynamic) num_threads(nthreads) collapse(2) shared(RiskPairs, RiskFail, R, Rd, Rls1, Lls1, Rls2, Lls2)
        #endif
        for (int ij = 0; ij < reqrdnum; ij++) {
            for (int j = 0; j < ntime; j++) {
                double Rs2 = 0;
                //
                vector<int> InGroup = RiskPairs[j];
                //  now has the grouping pairs
                int dj = RiskFail(j, 1)-RiskFail(j, 0) + 1;
                for (vector<double>::size_type i = 0; i < InGroup.size() - 1; i = i+2) {
                    Rs2 += Rd.block(InGroup[i] - 1, ij, InGroup[i + 1]-InGroup[i] + 1, 1).sum();
                }  //  precalculates the sums of risk groups
                //  only assigns values once
                Rls2(j, ij) = Rs2;
                Lls2(j, ij) = Rd.block(RiskFail(j, 0), ij, dj, 1).sum();
            }
        }
        if (!model_bool["gradient"]) {
            //
            #ifdef _OPENMP
            #pragma omp parallel for schedule(dynamic) num_threads(nthreads) collapse(2) shared(RiskPairs, RiskFail, R, Rd, Rdd, Rls1, Lls1, Rls2, Lls2, Rls3, Lls3)
            #endif
            for (int ijk = 0; ijk < reqrdnum*(reqrdnum + 1)/2; ijk++) {
                for (int j = 0; j < ntime; j++) {
                    double Rs3 = 0;
                    //
                    vector<int> InGroup = RiskPairs[j];
                    //  now has the grouping pairs
                    int dj = RiskFail(j, 1)-RiskFail(j, 0) + 1;
                    for (vector<double>::size_type i = 0; i < InGroup.size() - 1; i = i+2) {
                        Rs3 += Rdd.block(InGroup[i] - 1, ijk, InGroup[i + 1]-InGroup[i] + 1, 1).sum();
                    }  //  precalculates the sums of risk groups
                    //  only assigns values once
                    Rls3(j, ijk) = Rs3;
                    Lls3(j, ijk) = Rdd.block(RiskFail(j, 0), ijk, dj, 1).sum();
                }
            }
        }
    }
    return;
}

//' Utility function to calculate repeated values used in Cox Log-Likelihood calculation
//'
//' \code{Calculate_Sides_PO} Called to update repeated sum calculations, Uses list of event rows and risk matrices, Performs calculation of sums of risk in each group
//' @inheritParams CPP_template
//'
//' @return Updates matrices in place: risk storage matrices
//' @noRd
//'
void Calculate_Sides_PO(List& model_bool, const IntegerMatrix& RiskFail, const vector<vector<int> >& RiskPairs, const int& totalnum, const int& ntime, const MatrixXd& R, const MatrixXd& Rd, const MatrixXd& Rdd, MatrixXd& Rls1, MatrixXd& Rls2, MatrixXd& Rls3, MatrixXd& Lls1, MatrixXd& Lls2, MatrixXd& Lls3, const VectorXd& cens_weight, const int& nthreads, const IntegerVector& KeepConstant) {
    int reqrdnum = totalnum - sum(KeepConstant);
    #ifdef _OPENMP
    #pragma omp parallel for schedule(dynamic) num_threads(nthreads) shared(RiskPairs, RiskFail, R, Rls1, Lls1, cens_weight)
    #endif
    for (int j = 0; j < ntime; j++) {
        double Rs1 = 0;
        //
        vector<int> InGroup = RiskPairs[j];
        //  now has the grouping pairs
        int dj = RiskFail(j, 1)-RiskFail(j, 0) + 1;
        for (vector<double>::size_type i = 0; i < InGroup.size() - 1; i = i+2) {
            //
            Rs1 += R.block(InGroup[i] - 1, 0, InGroup[i + 1]-InGroup[i] + 1, 1).array().sum();
        }  //  precalculates the sums of risk groups
        VectorXd weighting = VectorXd::Zero(dj);
        weighting.head(dj) << cens_weight.segment(RiskFail(j, 0), dj);
        //  only assigns values once
        Rls1(j, 0) = Rs1;
        Lls1(j, 0) = (R.block(RiskFail(j, 0), 0, dj, 1).array() * weighting.array()).sum();
    }
    //
    if (!model_bool["single"]) {
        #ifdef _OPENMP
        #pragma omp parallel for schedule(dynamic) num_threads(nthreads) collapse(2) shared(RiskPairs, RiskFail, R, Rd, Rls1, Lls1, Rls2, Lls2, cens_weight)
        #endif
        for (int ij = 0; ij < reqrdnum; ij++) {
            for (int j = 0; j < ntime; j++) {
                double Rs2 = 0;
                //
                vector<int> InGroup = RiskPairs[j];
                //  now has the grouping pairs
                int dj = RiskFail(j, 1)-RiskFail(j, 0) + 1;
                for (vector<double>::size_type i = 0; i < InGroup.size() - 1; i = i+2) {
                    //
                    Rs2 += Rd.block(InGroup[i] - 1, ij, InGroup[i + 1]-InGroup[i] + 1, 1).array().sum();
                }  //  precalculates the sums of risk groups
                VectorXd weighting = VectorXd::Zero(dj);
                weighting.head(dj) << cens_weight.segment(RiskFail(j, 0), dj);
                //  only assigns values once
                Rls2(j, ij) = Rs2;
                Lls2(j, ij) = (Rd.block(RiskFail(j, 0), ij, dj, 1).array() * weighting.array()).sum();
            }
        }
        if (!model_bool["gradient"]) {
            //
            #ifdef _OPENMP
            #pragma omp parallel for schedule(dynamic) num_threads(nthreads) collapse(2) shared(RiskPairs, RiskFail, R, Rd, Rdd, Rls1, Lls1, Rls2, Lls2, Rls3, Lls3, cens_weight)
            #endif
            for (int ijk = 0; ijk < reqrdnum*(reqrdnum + 1)/2; ijk++) {
                for (int j = 0; j < ntime; j++) {
                    double Rs3 = 0;
                    //
                    vector<int> InGroup = RiskPairs[j];
                    //  now has the grouping pairs
                    int dj = RiskFail(j, 1)-RiskFail(j, 0) + 1;
                    for (vector<double>::size_type i = 0; i < InGroup.size() - 1; i = i + 2) {
                        //
                        Rs3 += Rdd.block(InGroup[i] - 1, ijk, InGroup[i + 1]-InGroup[i] + 1, 1).array().sum();
                    }  //  precalculates the sums of risk groups
                    VectorXd weighting = VectorXd::Zero(dj);
                    weighting.head(dj) << cens_weight.segment(RiskFail(j, 0), dj);
                    //  only assigns values once
                    Rls3(j, ijk) = Rs3;
                    Lls3(j, ijk) = (Rdd.block(RiskFail(j, 0), ijk, dj, 1).array() * weighting.array()).sum();
                }
            }
        }
    }
    return;
}

//' Utility function to calculate repeated values used in Cox Log-Likelihood calculation
//'
//' \code{Calculate_Sides_CR} Called to update repeated sum calculations, Uses list of event rows and risk matrices, Performs calculation of sums of risk in each group
//' @inheritParams CPP_template
//'
//' @return Updates matrices in place: risk storage matrices
//' @noRd
//'
void Calculate_Sides_CR(List& model_bool, const IntegerMatrix& RiskFail, const vector<vector<int> >& RiskPairs, const int& totalnum, const int& ntime, const MatrixXd& R, const MatrixXd& Rd, const MatrixXd& Rdd, MatrixXd& Rls1, MatrixXd& Rls2, MatrixXd& Rls3, MatrixXd& Lls1, MatrixXd& Lls2, MatrixXd& Lls3, const VectorXd& cens_weight, const int& nthreads, const IntegerVector& KeepConstant) {
    int reqrdnum = totalnum - sum(KeepConstant);
    #ifdef _OPENMP
    #pragma omp parallel for schedule(dynamic) num_threads(nthreads) shared(RiskPairs, RiskFail, R, Rls1, Lls1, cens_weight)
    #endif
    for (int j = 0; j < ntime; j++) {
        double Rs1 = 0;
        //
        vector<int> InGroup = RiskPairs[j];
        //  now has the grouping pairs
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
            Rs1 += (R.block(InGroup[i] - 1, 0, InGroup[i + 1]-InGroup[i] + 1, 1).array() * weighting.head(InGroup[i + 1]-InGroup[i] + 1).array()).sum();
        }  //  precalculates the sums of risk groups
        //  only assigns values once
        Rls1(j, 0) = Rs1;
        Lls1(j, 0) = R.block(RiskFail(j, 0), 0, dj, 1).sum();
    }
    if (!model_bool["single"]) {
        //
        #ifdef _OPENMP
        #pragma omp parallel for schedule(dynamic) num_threads(nthreads) collapse(2) shared(RiskPairs, RiskFail, R, Rd, Rls1, Lls1, Rls2, Lls2, cens_weight)
        #endif
        for (int ij = 0; ij < reqrdnum; ij++) {
            for (int j = 0; j < ntime; j++) {
                double Rs2 = 0;
                //
                vector<int> InGroup = RiskPairs[j];
                //  now has the grouping pairs
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
                    Rs2 += (Rd.block(InGroup[i] - 1, ij, InGroup[i + 1]-InGroup[i] + 1, 1).array() * weighting.head(InGroup[i + 1]-InGroup[i] + 1).array()).sum();
                }  //  precalculates the sums of risk groups
                //  only assigns values once
                Rls2(j, ij) = Rs2;
                Lls2(j, ij) = Rd.block(RiskFail(j, 0), ij, dj, 1).sum();
            }
        }
        if (!model_bool["gradient"]) {
            //
            #ifdef _OPENMP
            #pragma omp parallel for schedule(dynamic) num_threads(nthreads) collapse(2) shared(RiskPairs, RiskFail, R, Rd, Rdd, Rls1, Lls1, Rls2, Lls2, Rls3, Lls3, cens_weight)
            #endif
            for (int ijk = 0; ijk < reqrdnum*(reqrdnum + 1)/2; ijk++) {
                for (int j = 0; j < ntime; j++) {
                    double Rs3 = 0;
                    //
                    vector<int> InGroup = RiskPairs[j];
                    //  now has the grouping pairs
                    int dj = RiskFail(j, 1)-RiskFail(j, 0) + 1;
                    double cens_0 = cens_weight[RiskFail(j, 0)];
                    VectorXd weighting = VectorXd::Zero(InGroup[1]-InGroup[0] + 1);
                    for (vector<double>::size_type i = 0; i < InGroup.size() - 1; i = i + 2) {
                        if (weighting.size() != InGroup[i + 1]-InGroup[i] + 1) {
                            weighting.resize(InGroup[i + 1]-InGroup[i] + 1);
                        }
                        weighting.head(InGroup[i + 1]-InGroup[i] + 1) << cens_weight.segment(InGroup[i] - 1, InGroup[i + 1]-InGroup[i] + 1);
                        weighting = weighting / cens_0;
                        weighting = (weighting.array() < 1).select(weighting, 1);
                        //
                        Rs3 += (Rdd.block(InGroup[i] - 1, ijk, InGroup[i + 1]-InGroup[i] + 1, 1).array() * weighting.head(InGroup[i + 1]-InGroup[i] + 1).array()).sum();
                    }  //  precalculates the sums of risk groups
                    //  only assigns values once
                    Rls3(j, ijk) = Rs3;
                    Lls3(j, ijk) = Rdd.block(RiskFail(j, 0), ijk, dj, 1).sum();
                }
            }
        }
    }
    return;
}

//' Utility function to calculate repeated values used in Cox Log-Likelihood calculation with Strata
//'
//' \code{Calculate_Sides_Strata} Called to update repeated sum calculations, Uses list of event rows and risk matrices, Performs calculation of sums of risk in each group
//' @inheritParams CPP_template
//'
//' @return Updates matrices in place: risk storage matrices
//' @noRd
//'
void Calculate_Sides_Strata(List& model_bool, const IntegerMatrix& RiskFail, const vector<vector<vector<int> > >& RiskPairs_Strata, const int& totalnum, const int& ntime, const MatrixXd& R, const MatrixXd& Rd, const MatrixXd& Rdd, MatrixXd& Rls1, MatrixXd& Rls2, MatrixXd& Rls3, MatrixXd& Lls1, MatrixXd& Lls2, MatrixXd& Lls3, const int& nthreads, NumericVector& Strata_vals, const IntegerVector& KeepConstant) {
    int reqrdnum = totalnum - sum(KeepConstant);
    #ifdef _OPENMP
    #pragma omp parallel for schedule(dynamic) num_threads(nthreads) collapse(2)
    #endif
    for (int j = 0; j < ntime; j++) {
        for (int s_ij = 0; s_ij < Strata_vals.size(); s_ij++) {
            double Rs1 = 0;
            //
            vector<int> InGroup = RiskPairs_Strata[j][s_ij];
            //  now has the grouping pairs
            if (RiskFail(j, 2*s_ij + 1)> - 1) {
                int dj = RiskFail(j, 2*s_ij + 1)-RiskFail(j, 2*s_ij + 0) + 1;
                for (vector<double>::size_type i = 0; i < InGroup.size() - 1; i = i + 2) {
                    //
                    Rs1 += R.block(InGroup[i] - 1, 0, InGroup[i + 1]-InGroup[i] + 1, 1).sum();
                }  //  precalculates the sums of risk groups
                //  only assigns values once
                Rls1(j, s_ij) = Rs1;
                Lls1(j, s_ij) = R.block(RiskFail(j, 2*s_ij), 0, dj, 1).sum();
            }
        }
    }
    //
    if (!model_bool["single"]) {
        #ifdef _OPENMP
        #pragma omp parallel for schedule(dynamic) num_threads(nthreads) collapse(3)
        #endif
        for (int ij = 0; ij < reqrdnum; ij++) {  //  totalnum*(totalnum + 1)/2
            for (int j = 0; j < ntime; j++) {
                for (int s_ij = 0; s_ij < Strata_vals.size(); s_ij++) {
                    double Rs2 = 0;
                    //
                    vector<int> InGroup = RiskPairs_Strata[j][s_ij];
                    //  now has the grouping pairs
                    if (RiskFail(j, 2*s_ij + 1)> - 1) {
                        int dj = RiskFail(j, 2*s_ij + 1)-RiskFail(j, 2*s_ij + 0) + 1;
                        for (vector<double>::size_type i = 0; i < InGroup.size() - 1; i = i + 2) {
                            Rs2 += Rd.block(InGroup[i] - 1, ij, InGroup[i + 1]-InGroup[i] + 1, 1).sum();
                        }  //  precalculates the sums of risk groups
                        //  only assigns values once
                        Rls2(j, ij*Strata_vals.size() + s_ij) = Rs2;
                        Lls2(j, ij*Strata_vals.size() + s_ij) = Rd.block(RiskFail(j, 2*s_ij), ij, dj, 1).sum();
                    }
                }
            }
        }
        if (!model_bool["gradient"]) {
            //
            #ifdef _OPENMP
            #pragma omp parallel for schedule(dynamic) num_threads(nthreads) collapse(3)
            #endif
            for (int ijk = 0; ijk < reqrdnum*(reqrdnum + 1)/2; ijk++) {  //  totalnum*(totalnum + 1)/2
                for (int j = 0; j < ntime; j++) {
                    for (int s_ij = 0; s_ij < Strata_vals.size(); s_ij++) {
                        int ij = 0;
                        int jk = ijk;
                        while (jk > ij) {
                            ij++;
                            jk -= ij;
                        }
                        double Rs3 = 0;
                        //
                        vector<int> InGroup = RiskPairs_Strata[j][s_ij];
                        //  now has the grouping pairs
                        if (RiskFail(j, 2*s_ij + 1) >  - 1) {
                            int dj = RiskFail(j, 2*s_ij + 1)-RiskFail(j, 2*s_ij + 0) + 1;
                            for (vector<double>::size_type i = 0; i < InGroup.size() - 1; i = i + 2) {
                                Rs3 += Rdd.block(InGroup[i] - 1, ijk, InGroup[i + 1]-InGroup[i] + 1, 1).sum();
                            }  //  precalculates the sums of risk groups
                            //  only assigns values once
                            Rls3(j, ijk*Strata_vals.size() + s_ij) = Rs3;
                            Lls3(j, ijk*Strata_vals.size() + s_ij) = Rdd.block(RiskFail(j, 2*s_ij), ijk, dj, 1).sum();
                        }
                    }
                }
            }
        }
    }
    return;
}

//' Utility function to calculate repeated values used in Cox Log-Likelihood calculation with Strata and competing risks
//'
//' \code{Calculate_Sides_Strata_CR} Called to update repeated sum calculations, Uses list of event rows and risk matrices, Performs calculation of sums of risk in each group and competing risks
//' @inheritParams CPP_template
//'
//' @return Updates matrices in place: risk storage matrices
//' @noRd
//'
void Calculate_Sides_Strata_CR(List& model_bool, const IntegerMatrix& RiskFail, const vector<vector<vector<int> > >& RiskPairs_Strata, const int& totalnum, const int& ntime, const MatrixXd& R, const MatrixXd& Rd, const MatrixXd& Rdd, MatrixXd& Rls1, MatrixXd& Rls2, MatrixXd& Rls3, MatrixXd& Lls1, MatrixXd& Lls2, MatrixXd& Lls3, const VectorXd& cens_weight, const int& nthreads, NumericVector& Strata_vals, const IntegerVector& KeepConstant) {
    int reqrdnum = totalnum - sum(KeepConstant);
    #ifdef _OPENMP
    #pragma omp parallel for schedule(dynamic) num_threads(nthreads) collapse(2)
    #endif
    for (int j = 0; j < ntime; j++) {
        for (int s_ij = 0; s_ij < Strata_vals.size(); s_ij++) {
            double Rs1 = 0;
            //
            //
            vector<int> InGroup = RiskPairs_Strata[j][s_ij];
            //  now has the grouping pairs
            if (RiskFail(j, 2*s_ij + 1)> - 1) {
                int dj = RiskFail(j, 2*s_ij + 1)-RiskFail(j, 2*s_ij + 0) + 1;
                double cens_0 = cens_weight[RiskFail(j, 2*s_ij)];
                VectorXd weighting = VectorXd::Zero(InGroup[1]-InGroup[0] + 1);
                for (vector<double>::size_type i = 0; i < InGroup.size() - 1; i = i + 2) {
                    if (weighting.size() < InGroup[i + 1]-InGroup[i] + 1) {
                        weighting.resize(InGroup[i + 1]-InGroup[i] + 1);
                    }
                    weighting.head(InGroup[i + 1]-InGroup[i] + 1) << cens_weight.block(InGroup[i] - 1, 0, InGroup[i + 1]-InGroup[i] + 1, 1);
                    weighting = weighting / cens_0;
                    weighting = (weighting.array() < 1).select(weighting, 1);
                    //
                    Rs1 += (R.block(InGroup[i] - 1, 0, InGroup[i + 1]-InGroup[i] + 1, 1).array() * weighting.head(InGroup[i + 1]-InGroup[i] + 1).array()).sum();
                }  //  precalculates the sums of risk groups
                //  only assigns values once
                Rls1(j, s_ij) = Rs1;
                Lls1(j, s_ij) = R.block(RiskFail(j, 2*s_ij), 0, dj, 1).sum();
            }
        }
    }
    //
    if (!model_bool["single"]) {
        #ifdef _OPENMP
        #pragma omp parallel for schedule(dynamic) num_threads(nthreads) collapse(3)
        #endif
        for (int ij = 0; ij < reqrdnum; ij++) {  //  totalnum*(totalnum + 1)/2
            for (int j = 0; j < ntime; j++) {
                for (int s_ij = 0; s_ij < Strata_vals.size(); s_ij++) {
                    double Rs2 = 0;
                    //
                    vector<int> InGroup = RiskPairs_Strata[j][s_ij];
                    //  now has the grouping pairs
                    if (RiskFail(j, 2*s_ij + 1)> - 1) {
                        int dj = RiskFail(j, 2*s_ij + 1)-RiskFail(j, 2*s_ij + 0) + 1;
                        double cens_0 = cens_weight[RiskFail(j, 2*s_ij)];
                        VectorXd weighting = VectorXd::Zero(InGroup[1]-InGroup[0] + 1);
                        for (vector<double>::size_type i = 0; i < InGroup.size() - 1; i = i + 2) {
                            if (weighting.size() < InGroup[i + 1]-InGroup[i] + 1) {
                                weighting.resize(InGroup[i + 1]-InGroup[i] + 1);
                            }
                            weighting.head(InGroup[i + 1]-InGroup[i] + 1) << cens_weight.block(InGroup[i] - 1, 0, InGroup[i + 1]-InGroup[i] + 1, 1);
                            weighting = weighting / cens_0;
                            weighting = (weighting.array() < 1).select(weighting, 1);
                            //
                            Rs2 += (Rd.block(InGroup[i] - 1, ij, InGroup[i + 1]-InGroup[i] + 1, 1).array() * weighting.head(InGroup[i + 1]-InGroup[i] + 1).array()).sum();
                        }  //  precalculates the sums of risk groups
                        //  only assigns values once
                        Rls2(j, ij*Strata_vals.size() + s_ij) = Rs2;
                        Lls2(j, ij*Strata_vals.size() + s_ij) = Rd.block(RiskFail(j, 2*s_ij), ij, dj, 1).sum();
                    }
                }
            }
        }
        if (!model_bool["gradient"]) {
            //
            #ifdef _OPENMP
            #pragma omp parallel for schedule(dynamic) num_threads(nthreads) collapse(3)
            #endif
            for (int ijk = 0; ijk < reqrdnum*(reqrdnum + 1)/2; ijk++) {  //  totalnum*(totalnum + 1)/2
                for (int j = 0; j < ntime; j++) {
                    for (int s_ij = 0; s_ij < Strata_vals.size(); s_ij++) {
                        int ij = 0;
                        int jk = ijk;
                        while (jk > ij) {
                            ij++;
                            jk -= ij;
                        }
                        double Rs3 = 0;
                        //
                        vector<int> InGroup = RiskPairs_Strata[j][s_ij];
                        //  now has the grouping pairs
                        if (RiskFail(j, 2*s_ij + 1)> - 1) {
                            int dj = RiskFail(j, 2*s_ij + 1)-RiskFail(j, 2*s_ij + 0) + 1;
                            double cens_0 = cens_weight[RiskFail(j, 2*s_ij)];
                            VectorXd weighting = VectorXd::Zero(InGroup[1]-InGroup[0] + 1);
                            for (vector<double>::size_type i = 0; i < InGroup.size() - 1; i = i + 2) {
                                if (weighting.size() < InGroup[i + 1]-InGroup[i] + 1) {
                                    weighting.resize(InGroup[i + 1]-InGroup[i] + 1);
                                }
                                weighting.head(InGroup[i + 1]-InGroup[i] + 1) << cens_weight.block(InGroup[i] - 1, 0, InGroup[i + 1]-InGroup[i] + 1, 1);
                                weighting = weighting / cens_0;
                                weighting = (weighting.array() < 1).select(weighting, 1);
                                //
                                Rs3 += (Rdd.block(InGroup[i] - 1, ijk, InGroup[i + 1]-InGroup[i] + 1, 1).array() * weighting.head(InGroup[i + 1]-InGroup[i] + 1).array()).sum();
                            }  //  precalculates the sums of risk groups
                            //  only assigns values once
                            Rls3(j, ijk*Strata_vals.size() + s_ij) = Rs3;
                            Lls3(j, ijk*Strata_vals.size() + s_ij) = Rdd.block(RiskFail(j, 2*s_ij), ijk, dj, 1).sum();
                        }
                    }
                }
            }
        }
    }
    return;
}

//' Utility function to calculate Cox Log-Likelihood and derivatives
//'
//' \code{Calc_LogLik} Called to update log-likelihoods, Uses list of event rows, risk matrices, and repeated sums, Sums the log-likelihood contribution from each event time
//' @inheritParams CPP_template
//'
//' @return Updates matrices in place: Log-likelihood vectors/matrix
//' @noRd
//'
void Calc_LogLik(List& model_bool, const int& nthreads, const IntegerMatrix& RiskFail, const vector<vector<int> >& RiskPairs, const int& totalnum, const int& ntime, const MatrixXd& R, const MatrixXd& Rd, const MatrixXd& Rdd, const MatrixXd& RdR, const MatrixXd& RddR, const MatrixXd& Rls1, const MatrixXd& Rls2, const MatrixXd& Rls3, const MatrixXd& Lls1, const MatrixXd& Lls2, const MatrixXd& Lls3, vector<double>& Ll, vector<double>& Lld, vector<double>& Lldd, string ties_method, const IntegerVector& KeepConstant) {
    int reqrdnum = totalnum - sum(KeepConstant);
    fill(Ll.begin(), Ll.end(), 0.0);
    if (!model_bool["single"]) {
        fill(Lld.begin(), Lld.end(), 0.0);
        if (!model_bool["gradient"]) {
            fill(Lldd.begin(), Lldd.end(), 0.0);
        }
    }
    #ifdef _OPENMP
    #pragma omp declare reduction(vec_double_plus : vector<double> : \
        transform(omp_out.begin(), omp_out.end(), omp_in.begin(), omp_out.begin(), plus<double>())) \
        initializer(omp_priv = omp_orig)
    #endif
    if (model_bool["single"]) {
        #ifdef _OPENMP
        #pragma omp parallel for schedule(dynamic) num_threads(nthreads) reduction(vec_double_plus:Ll)
        #endif
        for (int j = 0; j < ntime; j++) {
            double Rs1 = Rls1(j, 0);
            //
            int dj = RiskFail(j, 1)-RiskFail(j, 0) + 1;
            MatrixXd Ldm = MatrixXd::Zero(dj, 1);
            double Ldcs = 0.0;
            if (ties_method == "efron") {
                Ldcs = Lls1(j, 0);
                for (int i = 0; i < dj; i++) {  //  adds in the efron approximation terms
                    Ldm(i, 0) = (-static_cast<double>(i) / static_cast<double>(dj)) * Ldcs;
                }
            }
            Ldm.col(0) = Ldm.col(0).array() + Rs1;
            //  Calculates the left-hand side terms
            MatrixXd temp1 = MatrixXd::Zero(dj, 1);
            temp1 = R.block(RiskFail(j, 0), 0, dj, 1).array().log();
            double Ld1 = (temp1.array().isFinite()).select(temp1, 0).sum();
            //  calculates the right-hand side terms
            temp1 = Ldm.col(0).array().log();
            Rs1 = (temp1.array().isFinite()).select(temp1, 0).sum();
            //
            Ll[0] += Ld1 - Rs1;
        }
    } else if (model_bool["gradient"]) {
        #ifdef _OPENMP
        #pragma omp parallel for schedule(dynamic) num_threads(nthreads) reduction(vec_double_plus:Ll, Lld) collapse(2)
        #endif
        for (int ij = 0; ij < reqrdnum; ij++) {  //  performs log-likelihood calculations for every derivative combination and risk group
            for (int j = 0; j < ntime; j++) {
                double Rs1 = Rls1(j, 0);
                double Rs2 = Rls2(j, ij);
                //
                int dj = RiskFail(j, 1)-RiskFail(j, 0) + 1;
                MatrixXd Ldm = MatrixXd::Zero(dj, 2);
                Vector2d Ldcs;
                if (ties_method == "efron") {
                    Ldcs << Lls1(j, 0), Lls2(j, ij);
                    for (int i = 0; i < dj; i++) {  //  adds in the efron approximation terms
                        Ldm.row(i) = (-static_cast<double>(i) / static_cast<double>(dj)) *Ldcs.array();
                    }
                }
                Ldm.col(0) = Ldm.col(0).array() + Rs1;
                Ldm.col(1) = Ldm.col(1).array() + Rs2;
                //  Calculates the left-hand side terms
                //
                double Ld1 = 0.0;
                double Ld2 = 0.0;
                //
                MatrixXd temp1 = MatrixXd::Zero(dj, 1);
                temp1 = R.block(RiskFail(j, 0), 0, dj, 1).array().log();
                Ld1 = (temp1.array().isFinite()).select(temp1, 0).sum();
                temp1 = RdR.block(RiskFail(j, 0), ij, dj, 1).array();
                Ld2 = (temp1.array().isFinite()).select(temp1, 0).sum();
                //  calculates the right-hand side terms
                temp1 = Ldm.col(0).array().log();
                Rs1 = (temp1.array().isFinite()).select(temp1, 0).sum();
                temp1 = Ldm.col(1).array() * (Ldm.col(0).array().pow(- 1).array());
                Rs2 = (temp1.array().isFinite()).select(temp1, 0).sum();
                //
                Ll[ij] += Ld1 - Rs1;
                Lld[ij] += Ld2 - Rs2;
            }
        }
    } else {
        #ifdef _OPENMP
        #pragma omp parallel for schedule(dynamic) num_threads(nthreads) reduction(vec_double_plus:Ll, Lld, Lldd) collapse(2)
        #endif
        for (int ijk = 0; ijk < reqrdnum*(reqrdnum + 1)/2; ijk++) {  //  performs log-likelihood calculations for every derivative combination and risk group
            for (int j = 0; j < ntime; j++) {
                int ij = 0;
                int jk = ijk;
                while (jk > ij) {
                    ij++;
                    jk -= ij;
                }
                double Rs1 = Rls1(j, 0);
                double Rs2 = Rls2(j, ij);
                double Rs2t = Rls2(j, jk);
                double Rs3 = Rls3(j, ijk);
                //
                int dj = RiskFail(j, 1)-RiskFail(j, 0) + 1;
                MatrixXd Ldm = MatrixXd::Zero(dj, 4);
                Vector4d Ldcs;
                if (ties_method == "efron") {
                    Ldcs << Lls1(j, 0), Lls2(j, ij), Lls2(j, jk), Lls3(j, ijk);
                    for (int i = 0; i < dj; i++) {  //  adds in the efron approximation terms
                        Ldm.row(i) = (-static_cast<double>(i) / static_cast<double>(dj)) *Ldcs.array();
                    }
                }
                Ldm.col(0) = Ldm.col(0).array() + Rs1;
                Ldm.col(1) = Ldm.col(1).array() + Rs2;
                Ldm.col(2) = Ldm.col(2).array() + Rs2t;
                Ldm.col(3) = Ldm.col(3).array() + Rs3;
                //  Calculates the left-hand side terms
                //
                double Ld1 = 0.0;
                double Ld2 = 0.0;
                double Ld3 = 0.0;
                //
                MatrixXd temp1 = MatrixXd::Zero(dj, 1);
                MatrixXd temp2 = MatrixXd::Zero(dj, 1);
                if (ij == jk) {
                    temp1 = R.block(RiskFail(j, 0), 0, dj, 1).array().log();
                    Ld1 = (temp1.array().isFinite()).select(temp1, 0).sum();
                }
                temp1 = RdR.block(RiskFail(j, 0), ij, dj, 1).array();
                if (ij == jk) {
                    Ld2 = (temp1.array().isFinite()).select(temp1, 0).sum();
                }
                temp2 = RdR.block(RiskFail(j, 0), jk, dj, 1).array();
                temp1 = RddR.block(RiskFail(j, 0), ijk, dj, 1).array() - (temp1.array() * temp2.array());
                Ld3 = (temp1.array().isFinite()).select(temp1, 0).sum();
                //  calculates the right-hand side terms
                if (ij == jk) {
                    temp1 = Ldm.col(0).array().log();
                    Rs1 = (temp1.array().isFinite()).select(temp1, 0).sum();
                }
                temp1 = Ldm.col(1).array() * (Ldm.col(0).array().pow(- 1).array());
                temp2 = Ldm.col(2).array() * (Ldm.col(0).array().pow(- 1).array());
                if (ij == jk) {
                    Rs2 = (temp1.array().isFinite()).select(temp1, 0).sum();
                }
                temp1 = Ldm.col(3).array() * (Ldm.col(0).array().pow(- 1).array()) - temp1.array() * temp2.array();
                Rs3 = (temp1.array().isFinite()).select(temp1, 0).sum();
                //
                if (ij == jk) {
                    Ll[ij] += Ld1 - Rs1;
                    Lld[ij] += Ld2 - Rs2;
                }
                Lldd[ij*reqrdnum+jk] += Ld3 - Rs3;  //  sums the log-likelihood and derivatives
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
            Lldd[jk*reqrdnum+ij] = Lldd[ij*reqrdnum+jk];
        }
    }
    double LogLik = 0;
    for (int i = 0; i < reqrdnum; i++) {
        if (Ll[i] != 0) {
            LogLik = Ll[i];
            break;
        }
    }
    fill(Ll.begin(), Ll.end(), LogLik);
    return;
}

//' Utility function to calculate Cox Log-Likelihood and derivatives with outcome probability
//'
//' \code{Calc_LogLik} Called to update log-likelihoods, Uses list of event rows, risk matrices, and repeated sums, Sums the log-likelihood contribution from each event time
//' @inheritParams CPP_template
//'
//' @return Updates matrices in place: Log-likelihood vectors/matrix
//' @noRd
//'
void Calc_LogLik_PO(List& model_bool, const int& nthreads, const IntegerMatrix& RiskFail, const vector<vector<int> >& RiskPairs, const int& totalnum, const int& ntime, const MatrixXd& R, const MatrixXd& Rd, const MatrixXd& Rdd, const MatrixXd& RdR, const MatrixXd& RddR, const MatrixXd& Rls1, const MatrixXd& Rls2, const MatrixXd& Rls3, const MatrixXd& Lls1, const MatrixXd& Lls2, const MatrixXd& Lls3, const VectorXd& cens_weight, vector<double>& Ll, vector<double>& Lld, vector<double>& Lldd, string ties_method, const IntegerVector& KeepConstant) {
    int reqrdnum = totalnum - sum(KeepConstant);
    fill(Ll.begin(), Ll.end(), 0.0);
    if (!model_bool["single"]) {
        fill(Lld.begin(), Lld.end(), 0.0);
        if (!model_bool["gradient"]) {
            fill(Lldd.begin(), Lldd.end(), 0.0);
        }
    }
    #ifdef _OPENMP
    #pragma omp declare reduction(vec_double_plus : vector<double> : \
        transform(omp_out.begin(), omp_out.end(), omp_in.begin(), omp_out.begin(), plus<double>())) \
        initializer(omp_priv = omp_orig)
    #endif
    if (model_bool["single"]) {
        #ifdef _OPENMP
        #pragma omp parallel for schedule(dynamic) num_threads(nthreads) reduction(vec_double_plus:Ll)
        #endif
        for (int j = 0; j < ntime; j++) {
            double Rs1 = Rls1(j, 0);
            //
            int dj = RiskFail(j, 1)-RiskFail(j, 0) + 1;
            VectorXd weighting = VectorXd::Zero(dj);
            weighting.head(dj) << cens_weight.segment(RiskFail(j, 0), dj);
            MatrixXd Ldm = MatrixXd::Zero(dj, 1);
            double Ldcs;
            if (ties_method == "efron") {
                Ldcs = Lls1(j, 0);
                for (int i = 0; i < dj; i++) {  //  adds in the efron approximation terms
                    Ldm(i, 0) = (-static_cast<double>(i) / static_cast<double>(dj)) * Ldcs;
                }
            }
            Ldm.col(0) = Ldm.col(0).array() + Rs1;
            //  Calculates the left-hand side terms
            //
            double Ld1 = 0.0;
            //
            MatrixXd temp1 = MatrixXd::Zero(dj, 1);
            temp1 = R.block(RiskFail(j, 0), 0, dj, 1).array().log();
            Ld1 = ((temp1.array().isFinite()).select(temp1, 0).array() * weighting.array()).sum();
            //  calculates the right-hand side terms
            temp1 = Ldm.col(0).array().log();
            Rs1 = ((temp1.array().isFinite()).select(temp1, 0).array() * weighting.array()).sum();
            //
            Ll[0] += Ld1 - Rs1;
        }
    } else if (model_bool["gradient"]) {
        #ifdef _OPENMP
        #pragma omp parallel for schedule(dynamic) num_threads(nthreads) reduction(vec_double_plus:Ll, Lld) collapse(2)
        #endif
        for (int ij = 0; ij < reqrdnum; ij++) {  //  performs log-likelihood calculations for every derivative combination and risk group
            for (int j = 0; j < ntime; j++) {
                double Rs1 = Rls1(j, 0);
                double Rs2 = Rls2(j, ij);
                //
                int dj = RiskFail(j, 1)-RiskFail(j, 0) + 1;
                VectorXd weighting = VectorXd::Zero(dj);
                weighting.head(dj) << cens_weight.segment(RiskFail(j, 0), dj);
                MatrixXd Ldm = MatrixXd::Zero(dj, 2);
                Vector2d Ldcs;
                if (ties_method == "efron") {
                    Ldcs << Lls1(j, 0), Lls2(j, ij);
                    for (int i = 0; i < dj; i++) {  //  adds in the efron approximation terms
                        Ldm.row(i) = (-static_cast<double>(i) / static_cast<double>(dj)) *Ldcs.array();
                    }
                }
                Ldm.col(0) = Ldm.col(0).array() + Rs1;
                Ldm.col(1) = Ldm.col(1).array() + Rs2;
                //  Calculates the left-hand side terms
                //
                double Ld1 = 0.0;
                double Ld2 = 0.0;
                //
                MatrixXd temp1 = MatrixXd::Zero(dj, 1);
                temp1 = R.block(RiskFail(j, 0), 0, dj, 1).array().log();
                Ld1 = ((temp1.array().isFinite()).select(temp1, 0).array() * weighting.array()).sum();
                temp1 = RdR.block(RiskFail(j, 0), ij, dj, 1).array();
                Ld2 = ((temp1.array().isFinite()).select(temp1, 0).array() * weighting.array()).sum();
                //  calculates the right-hand side terms
                temp1 = Ldm.col(0).array().log();
                Rs1 = ((temp1.array().isFinite()).select(temp1, 0).array() * weighting.array()).sum();
                temp1 = Ldm.col(1).array() * (Ldm.col(0).array().pow(- 1).array());
                Rs2 = ((temp1.array().isFinite()).select(temp1, 0).array() * weighting.array()).sum();
                //
                Ll[ij] += Ld1 - Rs1;
                Lld[ij] += Ld2 - Rs2;
            }
        }
    } else {
        #ifdef _OPENMP
        #pragma omp parallel for schedule(dynamic) num_threads(nthreads) reduction(vec_double_plus:Ll, Lld, Lldd) collapse(2)
        #endif
        for (int ijk = 0; ijk < reqrdnum*(reqrdnum + 1)/2; ijk++) {  //  performs log-likelihood calculations for every derivative combination and risk group
            for (int j = 0; j < ntime; j++) {
                int ij = 0;
                int jk = ijk;
                while (jk > ij) {
                    ij++;
                    jk -= ij;
                }
                double Rs1 = Rls1(j, 0);
                double Rs2 = Rls2(j, ij);
                double Rs2t = Rls2(j, jk);
                double Rs3 = Rls3(j, ijk);
                //
                int dj = RiskFail(j, 1)-RiskFail(j, 0) + 1;
                VectorXd weighting = VectorXd::Zero(dj);
                weighting.head(dj) << cens_weight.segment(RiskFail(j, 0), dj);
                MatrixXd Ldm = MatrixXd::Zero(dj, 4);
                Vector4d Ldcs;
                if (ties_method == "efron") {
                    Ldcs << Lls1(j, 0), Lls2(j, ij), Lls2(j, jk), Lls3(j, ijk);
                    for (int i = 0; i < dj; i++) {  //  adds in the efron approximation terms
                        Ldm.row(i) = (-static_cast<double>(i) / static_cast<double>(dj)) *Ldcs.array();
                    }
                }
                Ldm.col(0) = Ldm.col(0).array() + Rs1;
                Ldm.col(1) = Ldm.col(1).array() + Rs2;
                Ldm.col(2) = Ldm.col(2).array() + Rs2t;
                Ldm.col(3) = Ldm.col(3).array() + Rs3;
                //  Calculates the left-hand side terms
                //
                double Ld1 = 0.0;
                double Ld2 = 0.0;
                double Ld3 = 0.0;
                //
                MatrixXd temp1 = MatrixXd::Zero(dj, 1);
                MatrixXd temp2 = MatrixXd::Zero(dj, 1);
                if (ij == jk) {
                    temp1 = R.block(RiskFail(j, 0), 0, dj, 1).array().log();
                    Ld1 = ((temp1.array().isFinite()).select(temp1, 0).array() * weighting.array()).sum();
                }
                temp1 = RdR.block(RiskFail(j, 0), ij, dj, 1).array();
                if (ij == jk) {
                    Ld2 = ((temp1.array().isFinite()).select(temp1, 0).array() * weighting.array()).sum();
                }
                temp2 = RdR.block(RiskFail(j, 0), jk, dj, 1).array();
                temp1 = RddR.block(RiskFail(j, 0), ijk, dj, 1).array() - (temp1.array() * temp2.array());
                Ld3 = ((temp1.array().isFinite()).select(temp1, 0).array() * weighting.array()).sum();
                //  calculates the right-hand side terms
                if (ij == jk) {
                    temp1 = Ldm.col(0).array().log();
                    Rs1 = ((temp1.array().isFinite()).select(temp1, 0).array() * weighting.array()).sum();
                }
                temp1 = Ldm.col(1).array() * (Ldm.col(0).array().pow(- 1).array());
                temp2 = Ldm.col(2).array() * (Ldm.col(0).array().pow(- 1).array());
                if (ij == jk) {
                    Rs2 = ((temp1.array().isFinite()).select(temp1, 0).array() * weighting.array()).sum();
                }
                temp1 = Ldm.col(3).array() * (Ldm.col(0).array().pow(- 1).array()) - temp1.array() * temp2.array();
                Rs3 = ((temp1.array().isFinite()).select(temp1, 0).array() * weighting.array()).sum();
                //
                if (ij == jk) {
                    Ll[ij] += Ld1 - Rs1;
                    Lld[ij] += Ld2 - Rs2;
                }
                Lldd[ij*reqrdnum+jk] += Ld3 - Rs3;  //  sums the log-likelihood and derivatives
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
            Lldd[jk*reqrdnum+ij] = Lldd[ij*reqrdnum+jk];
        }
    }
    double LogLik = 0;
    for (int i = 0; i < reqrdnum; i++) {
        if (Ll[i] != 0) {
            LogLik = Ll[i];
            break;
        }
    }
    fill(Ll.begin(), Ll.end(), LogLik);
    return;
}

//' Utility function to calculate Cox Log-Likelihood and derivatives, basic model
//'
//' \code{Calc_LogLik_Basic} Basic model, Called to update log-likelihoods, Uses list of event rows, risk matrices, and repeated sums, Sums the log-likelihood contribution from each event time
//' @inheritParams CPP_template
//'
//' @return Updates matrices in place: Log-likelihood vectors/matrix
//' @noRd
//'
void Calc_LogLik_Basic(List& model_bool, const int& nthreads, const IntegerMatrix& RiskFail, const vector<vector<int> >& RiskPairs, const int& totalnum, const int& ntime, const MatrixXd& R, const MatrixXd& Rd, const MatrixXd& Rdd, const MatrixXd& RdR, const MatrixXd& Rls1, const MatrixXd& Rls2, const MatrixXd& Rls3, const MatrixXd& Lls1, const MatrixXd& Lls2, const MatrixXd& Lls3, vector<double>& Ll, vector<double>& Lld, vector<double>& Lldd, string ties_method, const IntegerVector& KeepConstant) {
    int reqrdnum = totalnum - sum(KeepConstant);
    fill(Ll.begin(), Ll.end(), 0.0);
    if (!model_bool["single"]) {
        fill(Lld.begin(), Lld.end(), 0.0);
        if (!model_bool["gradient"]) {
            fill(Lldd.begin(), Lldd.end(), 0.0);
        }
    }
    #ifdef _OPENMP
    #pragma omp declare reduction(vec_double_plus : vector<double> : \
        transform(omp_out.begin(), omp_out.end(), omp_in.begin(), omp_out.begin(), plus<double>())) \
        initializer(omp_priv = omp_orig)
    #endif
    if (model_bool["single"]) {
        #ifdef _OPENMP
        #pragma omp parallel for schedule(dynamic) num_threads(nthreads) reduction(vec_double_plus:Ll)
        #endif
        for (int j = 0; j < ntime; j++) {
            double Rs1 = Rls1(j, 0);
            //
            int dj = RiskFail(j, 1)-RiskFail(j, 0) + 1;
            MatrixXd Ldm = MatrixXd::Zero(dj, 1);
            double Ldcs;
            if (ties_method == "efron") {
                Ldcs = Lls1(j, 0);
                for (int i = 0; i < dj; i++) {  //  adds in the efron approximation terms
                    Ldm(i, 0) = (-static_cast<double>(i) / static_cast<double>(dj)) * Ldcs;
                }
            }
            Ldm.col(0) = Ldm.col(0).array() + Rs1;

            //  Calculates the left-hand side terms
            //
            double Ld1 = 0.0;

            //
            MatrixXd temp1 = MatrixXd::Zero(dj, 1);
            temp1 = R.block(RiskFail(j, 0), 0, dj, 1).array().log();
            Ld1 = (temp1.array().isFinite()).select(temp1, 0).sum();

            //  calculates the right-hand side terms
            temp1 = Ldm.col(0).array().log();
            Rs1 = (temp1.array().isFinite()).select(temp1, 0).sum();

            Ll[0] += Ld1 - Rs1;
        }
    } else if (model_bool["gradient"]) {
        #ifdef _OPENMP
        #pragma omp parallel for schedule(dynamic) num_threads(nthreads) reduction(vec_double_plus:Ll, Lld) collapse(2)
        #endif
        for (int ij = 0; ij < reqrdnum; ij++) {  //  performs log-likelihood calculations for every derivative combination and risk group
            for (int j = 0; j < ntime; j++) {
                double Rs1 = Rls1(j, 0);
                double Rs2 = Rls2(j, ij);
                //
                int dj = RiskFail(j, 1)-RiskFail(j, 0) + 1;
                MatrixXd Ldm = MatrixXd::Zero(dj, 2);
                Vector2d Ldcs;
                if (ties_method == "efron") {
                    Ldcs << Lls1(j, 0), Lls2(j, ij);
                    for (int i = 0; i < dj; i++) {  //  adds in the efron approximation terms
                        Ldm.row(i) = (-static_cast<double>(i) / static_cast<double>(dj)) *Ldcs.array();
                    }
                }
                Ldm.col(0) = Ldm.col(0).array() + Rs1;
                Ldm.col(1) = Ldm.col(1).array() + Rs2;
                //  Calculates the left-hand side terms
                //
                double Ld1 = 0.0;
                double Ld2 = 0.0;
                //
                MatrixXd temp1 = MatrixXd::Zero(dj, 1);
                temp1 = R.block(RiskFail(j, 0), 0, dj, 1).array().log();
                Ld1 = (temp1.array().isFinite()).select(temp1, 0).sum();
                temp1 = RdR.block(RiskFail(j, 0), ij, dj, 1).array();
                Ld2 = (temp1.array().isFinite()).select(temp1, 0).sum();
                //  calculates the right-hand side terms
                temp1 = Ldm.col(0).array().log();
                Rs1 = (temp1.array().isFinite()).select(temp1, 0).sum();
                temp1 = Ldm.col(1).array() * (Ldm.col(0).array().pow(- 1).array());
                Rs2 = (temp1.array().isFinite()).select(temp1, 0).sum();
                Ll[ij] += Ld1 - Rs1;
                Lld[ij] += Ld2 - Rs2;
            }
        }
    } else {
        #ifdef _OPENMP
        #pragma omp parallel for schedule(dynamic) num_threads(nthreads) reduction(vec_double_plus:Ll, Lld, Lldd) collapse(2)
        #endif
        for (int ijk = 0; ijk < reqrdnum*(reqrdnum + 1)/2; ijk++) {  //  performs log-likelihood calculations for every derivative combination and risk group
            for (int j = 0; j < ntime; j++) {
                int ij = 0;
                int jk = ijk;
                while (jk > ij) {
                    ij++;
                    jk -= ij;
                }
                double Rs1 = Rls1(j, 0);
                double Rs2 = Rls2(j, ij);
                double Rs2t = Rls2(j, jk);
                double Rs3 = Rls3(j, ijk);
                //
                int dj = RiskFail(j, 1)-RiskFail(j, 0) + 1;
                MatrixXd Ldm = MatrixXd::Zero(dj, 4);
                Vector4d Ldcs;
                if (ties_method == "efron") {
                    Ldcs << Lls1(j, 0), Lls2(j, ij), Lls2(j, jk), Lls3(j, ijk);
                    for (int i = 0; i < dj; i++) {  //  adds in the efron approximation terms
                        Ldm.row(i) = (-static_cast<double>(i) / static_cast<double>(dj)) *Ldcs.array();
                    }
                }
                Ldm.col(0) = Ldm.col(0).array() + Rs1;
                Ldm.col(1) = Ldm.col(1).array() + Rs2;
                Ldm.col(2) = Ldm.col(2).array() + Rs2t;
                Ldm.col(3) = Ldm.col(3).array() + Rs3;
                //  Calculates the left-hand side terms
                //
                double Ld1 = 0.0;
                double Ld2 = 0.0;
                //
                MatrixXd temp1 = MatrixXd::Zero(dj, 1);
                MatrixXd temp2 = MatrixXd::Zero(dj, 1);
                if (ij == jk) {
                    temp1 = R.block(RiskFail(j, 0), 0, dj, 1).array().log();
                    Ld1 = (temp1.array().isFinite()).select(temp1, 0).sum();
                    temp1 = RdR.block(RiskFail(j, 0), ij, dj, 1).array();
                    Ld2 = (temp1.array().isFinite()).select(temp1, 0).sum();
                }
                //  calculates the right-hand side terms
                if (ij == jk) {
                    temp1 = Ldm.col(0).array().log();
                    Rs1 = (temp1.array().isFinite()).select(temp1, 0).sum();
                }
                temp1 = Ldm.col(1).array() * (Ldm.col(0).array().pow(- 1).array());
                temp2 = Ldm.col(2).array() * (Ldm.col(0).array().pow(- 1).array());
                if (ij == jk) {
                    Rs2 = (temp1.array().isFinite()).select(temp1, 0).sum();
                }
                temp1 = Ldm.col(3).array() * (Ldm.col(0).array().pow(- 1).array()) - temp1.array() * temp2.array();
                Rs3 = (temp1.array().isFinite()).select(temp1, 0).sum();
                //
                if (ij == jk) {
                    Ll[ij] += Ld1 - Rs1;
                    Lld[ij] += Ld2 - Rs2;
                }
                Lldd[ij*reqrdnum+jk] += 0 - Rs3;  //  sums the log-likelihood and derivatives
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
            Lldd[jk*reqrdnum+ij] = Lldd[ij*reqrdnum+jk];
        }
    }
    double LogLik = 0;
    for (int i = 0; i < reqrdnum; i++) {
        if (Ll[i] != 0) {
            LogLik = Ll[i];
            break;
        }
    }
    fill(Ll.begin(), Ll.end(), LogLik);
    return;
}

//' Utility function to calculate Cox Log-Likelihood and derivatives with linear ERR simplification
//'
//' \code{Calc_LogLik_Linear_ERR} Called to update log-likelihoods, Uses list of event rows, risk matrices, and repeated sums, Sums the log-likelihood contribution from each event time
//' @inheritParams CPP_template
//'
//' @return Updates matrices in place: Log-likelihood vectors/matrix
//' @noRd
//'
void Calc_LogLik_Linear_ERR(List& model_bool, const StringVector& tform, const int& nthreads, const IntegerMatrix& RiskFail, const vector<vector<int> >& RiskPairs, const int& totalnum, const int& ntime, const MatrixXd& R, const MatrixXd& Rd, const MatrixXd& Rdd, const MatrixXd& RdR, const MatrixXd& RddR, const MatrixXd& Rls1, const MatrixXd& Rls2, const MatrixXd& Rls3, const MatrixXd& Lls1, const MatrixXd& Lls2, const MatrixXd& Lls3, vector<double>& Ll, vector<double>& Lld, vector<double>& Lldd, string ties_method, const IntegerVector& KeepConstant) {
    int reqrdnum = totalnum - sum(KeepConstant);
    fill(Ll.begin(), Ll.end(), 0.0);
    if (!model_bool["single"]) {
        fill(Lld.begin(), Lld.end(), 0.0);
        if (!model_bool["gradient"]) {
            fill(Lldd.begin(), Lldd.end(), 0.0);
        }
    }
    #ifdef _OPENMP
    #pragma omp declare reduction(vec_double_plus : vector<double> : \
        transform(omp_out.begin(), omp_out.end(), omp_in.begin(), omp_out.begin(), plus<double>())) \
        initializer(omp_priv = omp_orig)
    #endif
    if (model_bool["single"]) {
        #ifdef _OPENMP
        #pragma omp parallel for schedule(dynamic) num_threads(nthreads) reduction(vec_double_plus:Ll)
        #endif
        for (int j = 0; j < ntime; j++) {
            double Rs1 = Rls1(j, 0);
            //
            int dj = RiskFail(j, 1)-RiskFail(j, 0) + 1;
            MatrixXd Ldm = MatrixXd::Zero(dj, 1);
            double Ldcs;
            if (ties_method == "efron") {
                Ldcs = Lls1(j, 0);
                for (int i = 0; i < dj; i++) {  //  adds in the efron approximation terms
                    Ldm(i, 0) = (-static_cast<double>(i) / static_cast<double>(dj)) * Ldcs;
                }
            }
            Ldm.col(0) = Ldm.col(0).array() + Rs1;
            //  Calculates the left-hand side terms
            //
            double Ld1 = 0;
            //
            MatrixXd temp1 = MatrixXd::Zero(dj, 1);
            temp1 = R.block(RiskFail(j, 0), 0, dj, 1).array().log();
            Ld1 = (temp1.array().isFinite()).select(temp1, 0).sum();
            //  calculates the right-hand side terms
            temp1 = Ldm.col(0).array().log();
            Rs1 = (temp1.array().isFinite()).select(temp1, 0).sum();
            Ll[0] += Ld1 - Rs1;
        }
    } else if (model_bool["gradient"]) {
        #ifdef _OPENMP
        #pragma omp parallel for schedule(dynamic) num_threads(nthreads) reduction(vec_double_plus:Ll, Lld) collapse(2)
        #endif
        for (int t_ij = 0; t_ij < reqrdnum; t_ij++) {  //  performs log-likelihood calculations for every derivative combination and risk group
            for (int j = 0; j < ntime; j++) {
                if (KeepConstant[t_ij] == 0) {
                    int ij = t_ij - sum(head(KeepConstant, t_ij));
                    //
                    double Rs1 = Rls1(j, 0);
                    double Rs2 = Rls2(j, ij);
                    //
                    int dj = RiskFail(j, 1)-RiskFail(j, 0) + 1;
                    MatrixXd Ldm = MatrixXd::Zero(dj, 2);
                    Vector2d Ldcs;
                    if (ties_method == "efron") {
                        Ldcs << Lls1(j, 0), Lls2(j, ij);
                        for (int i = 0; i < dj; i++) {  //  adds in the efron approximation terms
                            Ldm.row(i) = (-static_cast<double>(i) / static_cast<double>(dj)) *Ldcs.array();
                        }
                    }
                    Ldm.col(0) = Ldm.col(0).array() + Rs1;
                    Ldm.col(1) = Ldm.col(1).array() + Rs2;
                    //  Calculates the left-hand side terms
                    //
                    double Ld1 = 0;
                    double Ld2 = 0;
                    //
                    MatrixXd temp1 = MatrixXd::Zero(dj, 1);
                    temp1 = R.block(RiskFail(j, 0), 0, dj, 1).array().log();
                    Ld1 = (temp1.array().isFinite()).select(temp1, 0).sum();
                    temp1 = RdR.block(RiskFail(j, 0), ij, dj, 1).array();
                    Ld2 = (temp1.array().isFinite()).select(temp1, 0).sum();
                    //  calculates the right-hand side terms
                    temp1 = Ldm.col(0).array().log();
                    Rs1 = (temp1.array().isFinite()).select(temp1, 0).sum();
                    temp1 = Ldm.col(1).array() * (Ldm.col(0).array().pow(- 1).array());
                    Rs2 = (temp1.array().isFinite()).select(temp1, 0).sum();
                    Ll[ij] += Ld1 - Rs1;
                    Lld[ij] += Ld2 - Rs2;
                }
            }
        }
    } else {
        #ifdef _OPENMP
        #pragma omp parallel for schedule(dynamic) num_threads(nthreads) reduction(vec_double_plus:Ll, Lld, Lldd) collapse(2)
        #endif
        for (int t_ijk = 0; t_ijk < reqrdnum*(reqrdnum + 1)/2; t_ijk++) {  //  performs log-likelihood calculations for every derivative combination and risk group
            for (int j = 0; j < ntime; j++) {
                int t_ij = 0;
                int t_jk = t_ijk;
                while (t_jk > t_ij) {
                    t_ij++;
                    t_jk -= t_ij;
                }
                if (KeepConstant[t_ij] + KeepConstant[t_jk] == 0) {
                    int ij = t_ij - sum(head(KeepConstant, t_ij));
                    int jk = t_jk - sum(head(KeepConstant, t_jk));
                    int ijk = ij*(ij + 1)/2 + jk;
                    //
                    double Rs1 = Rls1(j, 0);
                    double Rs2 = Rls2(j, ij);
                    double Rs2t = Rls2(j, jk);
                    double Rs3 = Rls3(j, ijk);
                    //
                    int dj = RiskFail(j, 1)-RiskFail(j, 0) + 1;
                    MatrixXd Ldm = MatrixXd::Zero(dj, 4);
                    Vector4d Ldcs;
                    if (ties_method == "efron") {
                        Ldcs << Lls1(j, 0), Lls2(j, ij), Lls2(j, jk), Lls3(j, ijk);
                        for (int i = 0; i < dj; i++) {  //  adds in the efron approximation terms
                            Ldm.row(i) = (-static_cast<double>(i) / static_cast<double>(dj)) *Ldcs.array();
                        }
                    }
                    Ldm.col(0) = Ldm.col(0).array() + Rs1;
                    Ldm.col(1) = Ldm.col(1).array() + Rs2;
                    Ldm.col(2) = Ldm.col(2).array() + Rs2t;
                    Ldm.col(3) = Ldm.col(3).array() + Rs3;
                    //  Calculates the left-hand side terms
                    //
                    double Ld1 = 0;
                    double Ld2 = 0;
                    double Ld3 = 0;
                    //
                    MatrixXd temp1 = MatrixXd::Zero(dj, 1);
                    MatrixXd temp2 = MatrixXd::Zero(dj, 1);
                    if (ij == jk) {
                        temp1 = R.block(RiskFail(j, 0), 0, dj, 1).array().log();
                        Ld1 = (temp1.array().isFinite()).select(temp1, 0).sum();
                    }
                    temp1 = RdR.block(RiskFail(j, 0), ij, dj, 1).array();
                    temp2 = RdR.block(RiskFail(j, 0), jk, dj, 1).array();
                    if (ij == jk) {
                        Ld2 = (temp1.array().isFinite()).select(temp1, 0).sum();
                    }
                    if (t_ij == t_jk) {
                        if (tform[t_ij] != "loglin") {
                            temp1 = RddR.block(RiskFail(j, 0), ijk, dj, 1).array() - (temp1.array() * temp2.array());
                            Ld3 = (temp1.array().isFinite()).select(temp1, 0).sum();
                        }
                    }
                    //  calculates the right-hand side terms
                    if (ij == jk) {
                        temp1 = Ldm.col(0).array().log();
                        Rs1 = (temp1.array().isFinite()).select(temp1, 0).sum();
                    }
                    temp1 = Ldm.col(1).array() * (Ldm.col(0).array().pow(- 1).array());
                    temp2 = Ldm.col(2).array() * (Ldm.col(0).array().pow(- 1).array());
                    if (ij == jk) {
                        Rs2 = (temp1.array().isFinite()).select(temp1, 0).sum();
                    }
                    temp1 = Ldm.col(3).array() * (Ldm.col(0).array().pow(- 1).array()) - temp1.array() * temp2.array();
                    Rs3 = (temp1.array().isFinite()).select(temp1, 0).sum();
                    //
                    if (ij == jk) {
                        Ll[ij] += Ld1 - Rs1;
                        Lld[ij] += Ld2 - Rs2;
                    }
                    Lldd[ij*reqrdnum+jk] += Ld3 - Rs3;  //  sums the log-likelihood and derivatives
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
            Lldd[jk*reqrdnum+ij] = Lldd[ij*reqrdnum+jk];
        }
    }
    double LogLik = 0;
    for (int i = 0; i < reqrdnum; i++) {
        if (Ll[i] != 0) {
            LogLik = Ll[i];
            break;
        }
    }
    fill(Ll.begin(), Ll.end(), LogLik);
    return;
}

//' Utility function to calculate Cox Log-Likelihood and derivatives with Strata and Linear ERR simplification
//'
//' \code{Calc_LogLik_Strata_Linear_ERR} Called to update log-likelihoods, Uses list of event rows, risk matrices, and repeated sums, Sums the log-likelihood contribution from each event time
//' @inheritParams CPP_template
//'
//' @return Updates matrices in place: Log-likelihood vectors/matrix
//' @noRd
//'
void Calc_LogLik_Strata_Linear_ERR(List& model_bool, const StringVector& tform, const int& nthreads, const IntegerMatrix& RiskFail, const vector<vector<vector<int> > >& RiskPairs_Strata, const int& totalnum, const int& ntime, const MatrixXd& R, const MatrixXd& Rd, const MatrixXd& Rdd, const MatrixXd& RdR, const MatrixXd& RddR, const MatrixXd& Rls1, const MatrixXd& Rls2, const MatrixXd& Rls3, const MatrixXd& Lls1, const MatrixXd& Lls2, const MatrixXd& Lls3, vector<double>& Ll, vector<double>& Lld, vector<double>& Lldd, string ties_method, NumericVector& Strata_vals, const IntegerVector& KeepConstant) {
    int reqrdnum = totalnum - sum(KeepConstant);
    fill(Ll.begin(), Ll.end(), 0.0);
    if (!model_bool["single"]) {
        fill(Lld.begin(), Lld.end(), 0.0);
        if (!model_bool["gradient"]) {
            fill(Lldd.begin(), Lldd.end(), 0.0);
        }
    }
    #ifdef _OPENMP
    #pragma omp declare reduction(vec_double_plus : vector<double> : \
        transform(omp_out.begin(), omp_out.end(), omp_in.begin(), omp_out.begin(), plus<double>())) \
        initializer(omp_priv = omp_orig)
    #endif
    if (model_bool["single"]) {
        #ifdef _OPENMP
        #pragma omp parallel for schedule(dynamic) num_threads(nthreads) reduction(vec_double_plus:Ll) collapse(2)
        #endif
        for (int j = 0; j < ntime; j++) {
            for (int s_ij = 0; s_ij < Strata_vals.size(); s_ij++) {
                //
                double Rs1 = Rls1(j, s_ij);
                //
                int dj = RiskFail(j, 2*s_ij + 1)-RiskFail(j, 2*s_ij + 0) + 1;
                if (RiskFail(j, 2*s_ij + 1)> - 1) {
                    //
                    MatrixXd Ldm = MatrixXd::Zero(dj, 1);
                    double Ldcs;
                    if (ties_method == "efron") {
                        Ldcs = Lls1(j, s_ij);
                        for (int i = 0; i < dj; i++) {  //  adds in the efron approximation terms
                            Ldm(i, 0) = (-static_cast<double>(i) / static_cast<double>(dj)) * Ldcs;
                        }
                    }
                    Ldm.col(0) = Ldm.col(0).array() + Rs1;
                    //  Calculates the left-hand side terms
                    //
                    double Ld1 = 0.0;
                    //
                    MatrixXd temp1 = MatrixXd::Zero(dj, 1);
                    temp1 = R.block(RiskFail(j, 2*s_ij), 0, dj, 1).array().log();
                    Ld1 = (temp1.array().isFinite()).select(temp1, 0).sum();
                    temp1 = Ldm.col(0).array().log();
                    Rs1 = (temp1.array().isFinite()).select(temp1, 0).sum();
                    Ll[0] += Ld1 - Rs1;
                }
            }
        }
    } else if (model_bool["gradient"]) {
        #ifdef _OPENMP
        #pragma omp parallel for schedule(dynamic) num_threads(nthreads) reduction(vec_double_plus:Ll, Lld) collapse(3)
        #endif
        for (int t_ij = 0; t_ij < totalnum; t_ij++) {  //  performs log-likelihood calculations for every derivative combination and risk group
            for (int j = 0; j < ntime; j++) {
                for (int s_ij = 0; s_ij < Strata_vals.size(); s_ij++) {
                    if (KeepConstant[t_ij] == 0) {
                        int ij = t_ij - sum(head(KeepConstant, t_ij));
                        //
                        double Rs1 = Rls1(j, s_ij);
                        double Rs2 = Rls2(j, ij*Strata_vals.size() + s_ij);
                        //
                        int dj = RiskFail(j, 2*s_ij + 1)-RiskFail(j, 2*s_ij + 0) + 1;
                        if (RiskFail(j, 2*s_ij + 1)> - 1) {
                            //
                            MatrixXd Ldm = MatrixXd::Zero(dj, 2);
                            Vector2d Ldcs;
                            if (ties_method == "efron") {
                                Ldcs << Lls1(j, s_ij), Lls2(j, ij*Strata_vals.size() + s_ij);
                                for (int i = 0; i < dj; i++) {  //  adds in the efron approximation terms
                                    Ldm.row(i) = (-static_cast<double>(i) / static_cast<double>(dj)) *Ldcs.array();
                                }
                            }
                            Ldm.col(0) = Ldm.col(0).array() + Rs1;
                            Ldm.col(1) = Ldm.col(1).array() + Rs2;
                            //  Calculates the left-hand side terms
                            //
                            double Ld1 = 0.0;
                            double Ld2 = 0.0;
                            //
                            MatrixXd temp1 = MatrixXd::Zero(dj, 1);
                            temp1 = R.block(RiskFail(j, 2*s_ij), 0, dj, 1).array().log();
                            Ld1 = (temp1.array().isFinite()).select(temp1, 0).sum();
                            temp1 = RdR.block(RiskFail(j, 2*s_ij), ij, dj, 1).array();
                            Ld2 = (temp1.array().isFinite()).select(temp1, 0).sum();
                            temp1 = Ldm.col(0).array().log();
                            Rs1 = (temp1.array().isFinite()).select(temp1, 0).sum();
                            temp1 = Ldm.col(1).array() * (Ldm.col(0).array().pow(- 1).array());
                            Rs2 = (temp1.array().isFinite()).select(temp1, 0).sum();
                            Ll[ij] += Ld1 - Rs1;
                            Lld[ij] += Ld2 - Rs2;
                        }
                    }
                }
            }
        }
    } else {
        #ifdef _OPENMP
        #pragma omp parallel for schedule(dynamic) num_threads(nthreads) reduction(vec_double_plus:Ll, Lld, Lldd) collapse(3)
        #endif
        for (int t_ijk = 0; t_ijk < totalnum*(totalnum + 1)/2; t_ijk++) {  //  performs log-likelihood calculations for every derivative combination and risk group
            for (int j = 0; j < ntime; j++) {
                for (int s_ij = 0; s_ij < Strata_vals.size(); s_ij++) {
                    int t_ij = 0;
                    int t_jk = t_ijk;
                    while (t_jk > t_ij) {
                        t_ij++;
                        t_jk -= t_ij;
                    }
                    if (KeepConstant[t_ij] + KeepConstant[t_jk] == 0) {
                        int ij = t_ij - sum(head(KeepConstant, t_ij));
                        int jk = t_jk - sum(head(KeepConstant, t_jk));
                        int ijk = ij*(ij + 1)/2 + jk;
                        //
                        double Rs1 = Rls1(j, s_ij);
                        double Rs2 = Rls2(j, ij*Strata_vals.size() + s_ij);
                        double Rs2t = Rls2(j, jk*Strata_vals.size() + s_ij);
                        double Rs3 = Rls3(j, ijk*Strata_vals.size() + s_ij);
                        //
                        int dj = RiskFail(j, 2*s_ij + 1)-RiskFail(j, 2*s_ij + 0) + 1;
                        if (RiskFail(j, 2*s_ij + 1)> - 1) {
                            //
                            MatrixXd Ldm = MatrixXd::Zero(dj, 4);
                            Vector4d Ldcs;
                            if (ties_method == "efron") {
                                Ldcs << Lls1(j, s_ij), Lls2(j, ij*Strata_vals.size() + s_ij), Lls2(j, jk*Strata_vals.size() + s_ij), Lls3(j, ijk*Strata_vals.size() + s_ij);
                                for (int i = 0; i < dj; i++) {  //  adds in the efron approximation terms
                                    Ldm.row(i) = (-static_cast<double>(i) / static_cast<double>(dj)) *Ldcs.array();
                                }
                            }
                            Ldm.col(0) = Ldm.col(0).array() + Rs1;
                            Ldm.col(1) = Ldm.col(1).array() + Rs2;
                            Ldm.col(2) = Ldm.col(2).array() + Rs2t;
                            Ldm.col(3) = Ldm.col(3).array() + Rs3;
                            //  Calculates the left-hand side terms
                            //
                            double Ld1 = 0.0;
                            double Ld2 = 0.0;
                            double Ld3 = 0.0;
                            //
                            MatrixXd temp1 = MatrixXd::Zero(dj, 1);
                            MatrixXd temp2 = MatrixXd::Zero(dj, 1);
                            if (ij == jk) {
                                temp1 = R.block(RiskFail(j, 2*s_ij), 0, dj, 1).array().log();
                                Ld1 = (temp1.array().isFinite()).select(temp1, 0).sum();
                            }
                            temp1 = RdR.block(RiskFail(j, 2*s_ij), ij, dj, 1).array();
                            temp2 = RdR.block(RiskFail(j, 2*s_ij), jk, dj, 1).array();
                            if (ij == jk) {
                                Ld2 = (temp1.array().isFinite()).select(temp1, 0).sum();
                            }
                            if (t_ij == t_jk) {
                                if (tform[t_ij] != "loglin") {
                                    temp1 = RddR.block(RiskFail(j, 2*s_ij), ijk, dj, 1).array() - (temp1.array() * temp2.array());
                                    Ld3 = (temp1.array().isFinite()).select(temp1, 0).sum();
                                }
                            }
                            //  calculates the right-hand side terms
                            if (ij == jk) {
                                temp1 = Ldm.col(0).array().log();
                                Rs1 = (temp1.array().isFinite()).select(temp1, 0).sum();
                            }
                            temp1 = Ldm.col(1).array() * (Ldm.col(0).array().pow(- 1).array());
                            temp2 = Ldm.col(2).array() * (Ldm.col(0).array().pow(- 1).array());
                            if (ij == jk) {
                                Rs2 = (temp1.array().isFinite()).select(temp1, 0).sum();
                            }
                            temp1 = Ldm.col(3).array() * (Ldm.col(0).array().pow(- 1).array()) - temp1.array() * temp2.array();
                            Rs3 = (temp1.array().isFinite()).select(temp1, 0).sum();
                            //
                            if (ij == jk) {
                                Ll[ij] += Ld1 - Rs1;
                                Lld[ij] += Ld2 - Rs2;
                            }
                            Lldd[ij*reqrdnum+jk] += Ld3 - Rs3;  //  sums the log-likelihood and derivatives
                        }
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
            Lldd[jk*reqrdnum+ij] = Lldd[ij*reqrdnum+jk];
        }
    }
    double LogLik = 0;
    for (int i = 0; i < reqrdnum; i++) {
        if (Ll[i] != 0) {
            LogLik = Ll[i];
            break;
        }
    }
    fill(Ll.begin(), Ll.end(), LogLik);
    return;
}

//' Utility function to calculate Cox Log-Likelihood and derivatives with Strata
//'
//' \code{Calc_LogLik_Strata} Called to update log-likelihoods, Uses list of event rows, risk matrices, and repeated sums, Sums the log-likelihood contribution from each event time
//' @inheritParams CPP_template
//'
//' @return Updates matrices in place: Log-likelihood vectors/matrix
//' @noRd
//'
void Calc_LogLik_Strata(List& model_bool, const int& nthreads, const IntegerMatrix& RiskFail, const vector<vector<vector<int> > >& RiskPairs_Strata, const int& totalnum, const int& ntime, const MatrixXd& R, const MatrixXd& Rd, const MatrixXd& Rdd, const MatrixXd& RdR, const MatrixXd& RddR, const MatrixXd& Rls1, const MatrixXd& Rls2, const MatrixXd& Rls3, const MatrixXd& Lls1, const MatrixXd& Lls2, const MatrixXd& Lls3, vector<double>& Ll, vector<double>& Lld, vector<double>& Lldd, string ties_method, NumericVector& Strata_vals, const IntegerVector& KeepConstant) {
    int reqrdnum = totalnum - sum(KeepConstant);
    fill(Ll.begin(), Ll.end(), 0.0);
    if (!model_bool["single"]) {
        fill(Lld.begin(), Lld.end(), 0.0);
        if (!model_bool["gradient"]) {
            fill(Lldd.begin(), Lldd.end(), 0.0);
        }
    }
    #ifdef _OPENMP
    #pragma omp declare reduction(vec_double_plus : vector<double> : \
        transform(omp_out.begin(), omp_out.end(), omp_in.begin(), omp_out.begin(), plus<double>())) \
        initializer(omp_priv = omp_orig)
    #endif
    if (model_bool["single"]) {
        #ifdef _OPENMP
        #pragma omp parallel for schedule(dynamic) num_threads(nthreads) reduction(vec_double_plus:Ll) collapse(2)
        #endif
        for (int j = 0; j < ntime; j++) {
            for (int s_ij = 0; s_ij < Strata_vals.size(); s_ij++) {
                double Rs1 = Rls1(j, s_ij);
                //
                int dj = RiskFail(j, 2*s_ij + 1)-RiskFail(j, 2*s_ij + 0) + 1;
                if (RiskFail(j, 2*s_ij + 1)> - 1) {
                    //
                    MatrixXd Ldm = MatrixXd::Zero(dj, 1);
                    double Ldcs = 0.0;
                    if (ties_method == "efron") {
                        Ldcs = Lls1(j, s_ij);
                        for (int i = 0; i < dj; i++) {  //  adds in the efron approximation terms
                            Ldm(i, 0) = (-static_cast<double>(i) / static_cast<double>(dj)) * Ldcs;
                        }
                    }
                    Ldm.col(0) = Ldm.col(0).array() + Rs1;
                    //  Calculates the left-hand side terms
                    //
                    double Ld1 = 0.0;
                    //
                    MatrixXd temp1 = MatrixXd::Zero(dj, 1);
                    temp1 = R.block(RiskFail(j, 2*s_ij), 0, dj, 1).array().log();
                    Ld1 = (temp1.array().isFinite()).select(temp1, 0).sum();
                    //
                    temp1 = Ldm.col(0).array().log();
                    Rs1 = (temp1.array().isFinite()).select(temp1, 0).sum();
                    Ll[0] += Ld1 - Rs1;
                }
            }
        }
    } else if (model_bool["gradient"]) {
        #ifdef _OPENMP
        #pragma omp parallel for schedule(dynamic) num_threads(nthreads) reduction(vec_double_plus:Ll, Lld) collapse(3)
        #endif
        for (int ij = 0; ij < reqrdnum; ij++) {  //  performs log-likelihood calculations for every derivative combination and risk group
            for (int j = 0; j < ntime; j++) {
                for (int s_ij = 0; s_ij < Strata_vals.size(); s_ij++) {
                    double Rs1 = Rls1(j, s_ij);
                    double Rs2 = Rls2(j, ij*Strata_vals.size() + s_ij);
                    //
                    int dj = RiskFail(j, 2*s_ij + 1)-RiskFail(j, 2*s_ij + 0) + 1;
                    if (RiskFail(j, 2*s_ij + 1)> - 1) {
                        //
                        MatrixXd Ldm = MatrixXd::Zero(dj, 2);
                        Vector2d Ldcs;
                        if (ties_method == "efron") {
                            Ldcs << Lls1(j, s_ij), Lls2(j, ij*Strata_vals.size() + s_ij);
                            for (int i = 0; i < dj; i++) {  //  adds in the efron approximation terms
                                Ldm.row(i) = (-static_cast<double>(i) / static_cast<double>(dj)) *Ldcs.array();
                            }
                        }
                        Ldm.col(0) = Ldm.col(0).array() + Rs1;
                        Ldm.col(1) = Ldm.col(1).array() + Rs2;
                        //  Calculates the left-hand side terms
                        //
                        double Ld1 = 0.0;
                        double Ld2 = 0.0;
                        //
                        MatrixXd temp1 = MatrixXd::Zero(dj, 1);
                        temp1 = R.block(RiskFail(j, 2*s_ij), 0, dj, 1).array().log();
                        Ld1 = (temp1.array().isFinite()).select(temp1, 0).sum();
                        temp1 = RdR.block(RiskFail(j, 2*s_ij), ij, dj, 1).array();
                        Ld2 = (temp1.array().isFinite()).select(temp1, 0).sum();
                        //  calculates the right-hand side terms
                        temp1 = Ldm.col(0).array().log();
                        Rs1 = (temp1.array().isFinite()).select(temp1, 0).sum();
                        temp1 = Ldm.col(1).array() * (Ldm.col(0).array().pow(- 1).array());
                        Rs2 = (temp1.array().isFinite()).select(temp1, 0).sum();
                        //
                        Ll[ij] += Ld1 - Rs1;
                        Lld[ij] += Ld2 - Rs2;
                    }
                }
            }
        }
    } else {
        #ifdef _OPENMP
        #pragma omp parallel for schedule(dynamic) num_threads(nthreads) reduction(vec_double_plus:Ll, Lld, Lldd) collapse(3)
        #endif
        for (int ijk = 0; ijk < reqrdnum*(reqrdnum + 1)/2; ijk++) {  //  performs log-likelihood calculations for every derivative combination and risk group
            for (int j = 0; j < ntime; j++) {
                for (int s_ij = 0; s_ij < Strata_vals.size(); s_ij++) {
                    int ij = 0;
                    int jk = ijk;
                    while (jk > ij) {
                        ij++;
                        jk -= ij;
                    }
                    double Rs1 = Rls1(j, s_ij);
                    double Rs2 = Rls2(j, ij*Strata_vals.size() + s_ij);
                    double Rs2t = Rls2(j, jk*Strata_vals.size() + s_ij);
                    double Rs3 = Rls3(j, ijk*Strata_vals.size() + s_ij);
                    //
                    int dj = RiskFail(j, 2*s_ij + 1)-RiskFail(j, 2*s_ij + 0) + 1;
                    if (RiskFail(j, 2*s_ij + 1)> - 1) {
                        //
                        MatrixXd Ldm = MatrixXd::Zero(dj, 4);
                        Vector4d Ldcs;
                        if (ties_method == "efron") {
                            Ldcs << Lls1(j, s_ij), Lls2(j, ij*Strata_vals.size() + s_ij), Lls2(j, jk*Strata_vals.size() + s_ij), Lls3(j, ijk*Strata_vals.size() + s_ij);
                            for (int i = 0; i < dj; i++) {  //  adds in the efron approximation terms
                                Ldm.row(i) = (-static_cast<double>(i) / static_cast<double>(dj)) *Ldcs.array();
                            }
                        }
                        Ldm.col(0) = Ldm.col(0).array() + Rs1;
                        Ldm.col(1) = Ldm.col(1).array() + Rs2;
                        Ldm.col(2) = Ldm.col(2).array() + Rs2t;
                        Ldm.col(3) = Ldm.col(3).array() + Rs3;
                        //  Calculates the left-hand side terms
                        //
                        double Ld1 = 0.0;
                        double Ld2 = 0.0;
                        double Ld3 = 0.0;
                        //
                        MatrixXd temp1 = MatrixXd::Zero(dj, 1);
                        MatrixXd temp2 = MatrixXd::Zero(dj, 1);
                        if (ij == jk) {
                            temp1 = R.block(RiskFail(j, 2*s_ij), 0, dj, 1).array().log();
                            Ld1 = (temp1.array().isFinite()).select(temp1, 0).sum();
                        }
                        temp1 = RdR.block(RiskFail(j, 2*s_ij), ij, dj, 1).array();
                        temp2 = RdR.block(RiskFail(j, 2*s_ij), jk, dj, 1).array();
                        if (ij == jk) {
                            Ld2 = (temp1.array().isFinite()).select(temp1, 0).sum();
                        }
                        temp1 = RddR.block(RiskFail(j, 2*s_ij), ijk, dj, 1).array() - (temp1.array() * temp2.array());
                        Ld3 = (temp1.array().isFinite()).select(temp1, 0).sum();
                        //  calculates the right-hand side terms
                        if (ij == jk) {
                            temp1 = Ldm.col(0).array().log();
                            Rs1 = (temp1.array().isFinite()).select(temp1, 0).sum();
                        }
                        temp1 = Ldm.col(1).array() * (Ldm.col(0).array().pow(- 1).array());
                        temp2 = Ldm.col(2).array() * (Ldm.col(0).array().pow(- 1).array());
                        if (ij == jk) {
                            Rs2 = (temp1.array().isFinite()).select(temp1, 0).sum();
                        }
                        temp1 = Ldm.col(3).array() * (Ldm.col(0).array().pow(- 1).array()) - temp1.array() * temp2.array();
                        Rs3 = (temp1.array().isFinite()).select(temp1, 0).sum();
                        //
                        if (ij == jk) {
                            Ll[ij] += Ld1 - Rs1;
                            Lld[ij] += Ld2 - Rs2;
                        }
                        Lldd[ij*reqrdnum+jk] += Ld3 - Rs3;  //  sums the log-likelihood and derivatives
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
            Lldd[jk*reqrdnum+ij] = Lldd[ij*reqrdnum+jk];
        }
    }
    double LogLik = 0;
    for (int i = 0; i < reqrdnum; i++) {
        if (Ll[i] != 0) {
            LogLik = Ll[i];
            break;
        }
    }
    fill(Ll.begin(), Ll.end(), LogLik);
    return;
}

//' Utility function to calculate Cox Log-Likelihood and derivatives with Strata, basic model
//'
//' \code{Calc_LogLik_Strata_BASIC} Called to update log-likelihoods, Uses list of event rows, risk matrices, and repeated sums, Sums the log-likelihood contribution from each event time, basic model
//' @inheritParams CPP_template
//'
//' @return Updates matrices in place: Log-likelihood vectors/matrix
//' @noRd
//'
void Calc_LogLik_Strata_Basic(List& model_bool, const int& nthreads, const IntegerMatrix& RiskFail, const vector<vector<vector<int> > >& RiskPairs_Strata, const int& totalnum, const int& ntime, const MatrixXd& R, const MatrixXd& Rd, const MatrixXd& Rdd, const MatrixXd& RdR, const MatrixXd& Rls1, const MatrixXd& Rls2, const MatrixXd& Rls3, const MatrixXd& Lls1, const MatrixXd& Lls2, const MatrixXd& Lls3, vector<double>& Ll, vector<double>& Lld, vector<double>& Lldd, string ties_method, NumericVector& Strata_vals, const IntegerVector& KeepConstant) {
    int reqrdnum = totalnum - sum(KeepConstant);
    fill(Ll.begin(), Ll.end(), 0.0);
    if (!model_bool["single"]) {
        fill(Lld.begin(), Lld.end(), 0.0);
        if (!model_bool["gradient"]) {
            fill(Lldd.begin(), Lldd.end(), 0.0);
        }
    }
    #ifdef _OPENMP
    #pragma omp declare reduction(vec_double_plus : vector<double> : \
        transform(omp_out.begin(), omp_out.end(), omp_in.begin(), omp_out.begin(), plus<double>())) \
        initializer(omp_priv = omp_orig)
    #endif
    if (model_bool["single"]) {
        #ifdef _OPENMP
        #pragma omp parallel for schedule(dynamic) num_threads(nthreads) reduction(vec_double_plus:Ll) collapse(2)
        #endif
        for (int j = 0; j < ntime; j++) {
            for (int s_ij = 0; s_ij < Strata_vals.size(); s_ij++) {
                if (RiskFail(j, 2*s_ij + 1)> - 1) {
                    double Rs1 = Rls1(j, s_ij);
                    //
                    int dj = RiskFail(j, 2*s_ij + 1)-RiskFail(j, 2*s_ij + 0) + 1;
                    //
                    MatrixXd Ldm = MatrixXd::Zero(dj, 1);
                    double Ldcs;
                    if (ties_method == "efron") {
                        Ldcs = Lls1(j, s_ij);
                        for (int i = 0; i < dj; i++) {  //  adds in the efron approximation terms
                            Ldm(i, 0) = (-static_cast<double>(i) / static_cast<double>(dj)) * Ldcs;
                        }
                    }
                    Ldm.col(0) = Ldm.col(0).array() + Rs1;
                    //  Calculates the left-hand side terms
                    //
                    double Ld1 = 0.0;
                    //
                    MatrixXd temp1 = MatrixXd::Zero(dj, 1);
                    temp1 = R.block(RiskFail(j, 2*s_ij), 0, dj, 1).array().log();
                    Ld1 = (temp1.array().isFinite()).select(temp1, 0).sum();
                    //  calculates the right-hand side terms
                    temp1 = Ldm.col(0).array().log();
                    Rs1 = (temp1.array().isFinite()).select(temp1, 0).sum();
                    Ll[0] += Ld1 - Rs1;
                }
            }
        }
    } else if (model_bool["gradient"]) {
        #ifdef _OPENMP
        #pragma omp parallel for schedule(dynamic) num_threads(nthreads) reduction(vec_double_plus:Ll, Lld) collapse(3)
        #endif
        for (int ij = 0; ij < reqrdnum; ij++) {  //  performs log-likelihood calculations for every derivative combination and risk group
            for (int j = 0; j < ntime; j++) {
                for (int s_ij = 0; s_ij < Strata_vals.size(); s_ij++) {
                    if (RiskFail(j, 2*s_ij + 1)> - 1) {
                        double Rs1 = Rls1(j, s_ij);
                        double Rs2 = Rls2(j, ij*Strata_vals.size() + s_ij);
                        //
                        int dj = RiskFail(j, 2*s_ij + 1)-RiskFail(j, 2*s_ij + 0) + 1;
                        //
                        MatrixXd Ldm = MatrixXd::Zero(dj, 2);
                        Vector2d Ldcs;
                        if (ties_method == "efron") {
                            Ldcs << Lls1(j, s_ij), Lls2(j, ij*Strata_vals.size() + s_ij);
                            for (int i = 0; i < dj; i++) {  //  adds in the efron approximation terms
                                Ldm.row(i) = (-static_cast<double>(i) / static_cast<double>(dj)) *Ldcs.array();
                            }
                        }
                        Ldm.col(0) = Ldm.col(0).array() + Rs1;
                        Ldm.col(1) = Ldm.col(1).array() + Rs2;
                        //  Calculates the left-hand side terms
                        //
                        double Ld1 = 0.0;
                        double Ld2 = 0.0;
                        //
                        MatrixXd temp1 = MatrixXd::Zero(dj, 1);
                        temp1 = R.block(RiskFail(j, 2*s_ij), 0, dj, 1).array().log();
                        Ld1 = (temp1.array().isFinite()).select(temp1, 0).sum();
                        temp1 = RdR.block(RiskFail(j, 2*s_ij), ij, dj, 1).array();
                        Ld2 = (temp1.array().isFinite()).select(temp1, 0).sum();
                        //  calculates the right-hand side terms
                        temp1 = Ldm.col(0).array().log();
                        Rs1 = (temp1.array().isFinite()).select(temp1, 0).sum();
                        temp1 = Ldm.col(1).array() * (Ldm.col(0).array().pow(- 1).array());
                        Rs2 = (temp1.array().isFinite()).select(temp1, 0).sum();
                        Ll[ij] += Ld1 - Rs1;
                        Lld[ij] += Ld2 - Rs2;
                    }
                }
            }
        }
    } else {
        #ifdef _OPENMP
        #pragma omp parallel for schedule(dynamic) num_threads(nthreads) reduction(vec_double_plus:Ll, Lld, Lldd) collapse(3)
        #endif
        for (int ijk = 0; ijk < reqrdnum*(reqrdnum + 1)/2; ijk++) {  //  performs log-likelihood calculations for every derivative combination and risk group
            for (int j = 0; j < ntime; j++) {
                for (int s_ij = 0; s_ij < Strata_vals.size(); s_ij++) {
                    int ij = 0;
                    int jk = ijk;
                    while (jk > ij) {
                        ij++;
                        jk -= ij;
                    }
                    if (RiskFail(j, 2*s_ij + 1)> - 1) {
                        double Rs1 = Rls1(j, s_ij);
                        double Rs2 = Rls2(j, ij*Strata_vals.size() + s_ij);
                        double Rs2t = Rls2(j, jk*Strata_vals.size() + s_ij);
                        double Rs3 = Rls3(j, ijk*Strata_vals.size() + s_ij);
                        //
                        int dj = RiskFail(j, 2*s_ij + 1)-RiskFail(j, 2*s_ij + 0) + 1;
                        //
                        MatrixXd Ldm = MatrixXd::Zero(dj, 4);
                        Vector4d Ldcs;
                        if (ties_method == "efron") {
                            Ldcs << Lls1(j, s_ij), Lls2(j, ij*Strata_vals.size() + s_ij), Lls2(j, jk*Strata_vals.size() + s_ij), Lls3(j, ijk*Strata_vals.size() + s_ij);
                            for (int i = 0; i < dj; i++) {  //  adds in the efron approximation terms
                                Ldm.row(i) = (-static_cast<double>(i) / static_cast<double>(dj)) *Ldcs.array();
                            }
                        }
                        Ldm.col(0) = Ldm.col(0).array() + Rs1;
                        Ldm.col(1) = Ldm.col(1).array() + Rs2;
                        Ldm.col(2) = Ldm.col(2).array() + Rs2t;
                        Ldm.col(3) = Ldm.col(3).array() + Rs3;
                        //  Calculates the left-hand side terms
                        //
                        double Ld1 = 0.0;
                        double Ld2 = 0.0;
                        //
                        MatrixXd temp1 = MatrixXd::Zero(dj, 1);
                        MatrixXd temp2 = MatrixXd::Zero(dj, 1);
                        if (ij == jk) {
                            temp1 = R.block(RiskFail(j, 2*s_ij), 0, dj, 1).array().log();
                            Ld1 = (temp1.array().isFinite()).select(temp1, 0).sum();
                            temp1 = RdR.block(RiskFail(j, 2*s_ij), ij, dj, 1).array();
                            Ld2 = (temp1.array().isFinite()).select(temp1, 0).sum();
                        }
                        //  calculates the right-hand side terms
                        if (ij == jk) {
                            temp1 = Ldm.col(0).array().log();
                            Rs1 = (temp1.array().isFinite()).select(temp1, 0).sum();
                        }
                        temp1 = Ldm.col(1).array() * (Ldm.col(0).array().pow(- 1).array());
                        temp2 = Ldm.col(2).array() * (Ldm.col(0).array().pow(- 1).array());
                        if (ij == jk) {
                            Rs2 = (temp1.array().isFinite()).select(temp1, 0).sum();
                        }
                        temp1 = Ldm.col(3).array() * (Ldm.col(0).array().pow(- 1).array()) - temp1.array() * temp2.array();
                        Rs3 = (temp1.array().isFinite()).select(temp1, 0).sum();
                        //
                        if (ij == jk) {
                            Ll[ij] += Ld1 - Rs1;
                            Lld[ij] += Ld2 - Rs2;
                        }
                        Lldd[ij*reqrdnum+jk] += 0 - Rs3;  //  sums the log-likelihood and derivatives
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
            Lldd[jk*reqrdnum+ij] = Lldd[ij*reqrdnum+jk];
        }
    }
    double LogLik = 0;
    for (int i = 0; i < reqrdnum; i++) {
        if (Ll[i] != 0) {
            LogLik = Ll[i];
            break;
        }
    }
    fill(Ll.begin(), Ll.end(), LogLik);
    return;
}

//' Utility function to calculate poisson log-likelihood and derivatives
//'
//' \code{Poisson_LogLik} Called to update log-likelihoods, Uses list risk matrices and person-years, Sums the log-likelihood contribution from each row
//' @inheritParams CPP_template
//'
//' @return Updates matrices in place: Log-likelihood vectors/matrix
//' @noRd
//'
void Poisson_LogLik(List& model_bool, const int& nthreads, const int& totalnum, const Ref<const MatrixXd>& PyrC, const MatrixXd& R, const MatrixXd& Rd, const MatrixXd& Rdd, const MatrixXd& RdR, const MatrixXd& RddR, vector<double>& Ll, vector<double>& Lld, vector<double>& Lldd, const IntegerVector& KeepConstant) {
    int reqrdnum = totalnum - sum(KeepConstant);
    fill(Ll.begin(), Ll.end(), 0.0);
    if (!model_bool["single"]) {
        fill(Lld.begin(), Lld.end(), 0.0);
        if (!model_bool["gradient"]) {
            fill(Lldd.begin(), Lldd.end(), 0.0);
        }
    }
    MatrixXd temp(Rd.rows(), Rd.cols());
    temp = (PyrC.col(1).array() * (PyrC.col(0).array() * R.col(0).array()).array().log()).array() - (PyrC.col(0).array() * R.col(0).array());
    fill(Ll.begin(), Ll.end(), (temp.array().isFinite()).select(temp, 0).sum());
    if (!model_bool["single"]) {
        VectorXd CoL = VectorXd::Zero(Rd.rows());
        CoL = PyrC.col(1).array() * R.col(0).array().pow(- 1).array();
        if (model_bool["gradient"]) {
            #ifdef _OPENMP
            #pragma omp parallel for schedule(dynamic) num_threads(nthreads)
            #endif
            for (int ij = 0; ij < reqrdnum; ij++) {  //
                VectorXd temp(Rd.rows(), 1);
                temp = Rd.col(ij).array() * (CoL.array() - PyrC.col(0).array());
                Lld[ij] = (temp.array().isFinite()).select(temp, 0).sum();
            }
        } else {
            #ifdef _OPENMP
            #pragma omp parallel for schedule(dynamic) num_threads(nthreads)
            #endif
            for (int ijk = 0; ijk < reqrdnum*(reqrdnum + 1)/2; ijk++) {  //
                int ij = 0;
                int jk = ijk;
                while (jk > ij) {
                    ij++;
                    jk -= ij;
                }
                VectorXd temp(Rdd.rows(), 1);
                temp = Rdd.col(ijk).array() * (CoL.array() - PyrC.col(0).array()) - PyrC.col(1).array() * RdR.col(ij).array() * RdR.col(jk).array();
                Lldd[ij*reqrdnum+jk] = (temp.array().isFinite()).select(temp, 0).sum();
                if (ij != jk) {
                    Lldd[jk*reqrdnum+ij] = (temp.array().isFinite()).select(temp, 0).sum();
                } else {
                    temp = Rd.col(ij).array() * (CoL.array() - PyrC.col(0).array());
                    Lld[ij] = (temp.array().isFinite()).select(temp, 0).sum();
                }
            }
        }
    }
    return;
}

//' Utility function to perform null model equivalent of Calculate_Sides
//'
//' \code{Calculate_Null_Sides} Called to update repeated sum calculations, Uses list of event rows, Performs calculation of counts in each group
//' @inheritParams CPP_template
//'
//' @return Updates matrices in place: risk storage matrices
//' @noRd
//'
void Calculate_Null_Sides(const IntegerMatrix& RiskFail, const vector<vector<int> >& RiskPairs, const int& ntime, const MatrixXd& R, MatrixXd& Rls1, MatrixXd& Lls1, const int& nthreads) {
    #ifdef _OPENMP
    #pragma omp parallel for schedule(dynamic) num_threads(nthreads)
    #endif
    for (int j = 0; j < ntime; j++) {
        double Rs1 = 0;
        //
        vector<int> InGroup = RiskPairs[j];
        //  now has the grouping pairs
        int dj = RiskFail(j, 1)-RiskFail(j, 0) + 1;
        for (vector<double>::size_type i = 0; i < InGroup.size() - 1; i = i + 2) {
            Rs1 += R.block(InGroup[i] - 1, 0, InGroup[i + 1]-InGroup[i] + 1, 1).sum();
        }  //  precalculates the sums of risk groups
        //  only assigns values once
        Rls1(j, 0) = Rs1;
        Lls1(j, 0) = R.block(RiskFail(j, 0), 0, dj, 1).sum();
    }
    return;
}


//' Utility function to perform null model equivalent of Calc_LogLik
//'
//' \code{Calc_Null_LogLik} Called to update log-likelihoods, Uses list of event rows and repeated sums, Sums the log-likelihood contribution from each event time
//' @inheritParams CPP_template
//'
//' @return Updates matrices in place: Log-likelihood vectors/matrix
//' @noRd
//'
void Calc_Null_LogLik(const int& nthreads, const IntegerMatrix& RiskFail, const vector<vector<int> >& RiskPairs, const int& ntime, const MatrixXd& R, const MatrixXd& Rls1, const MatrixXd& Lls1, vector<double>& Ll, string ties_method) {
    #ifdef _OPENMP
    #pragma omp declare reduction(vec_double_plus : vector<double> : \
        transform(omp_out.begin(), omp_out.end(), omp_in.begin(), omp_out.begin(), plus<double>())) \
        initializer(omp_priv = omp_orig)
    #pragma omp parallel for schedule(dynamic) num_threads(nthreads) reduction(vec_double_plus:Ll)
    #endif
    for (int j = 0; j < ntime; j++) {
        double Rs1 = Rls1(j, 0);
        int dj = RiskFail(j, 1)-RiskFail(j, 0) + 1;
        //
        MatrixXd Ldm = MatrixXd::Zero(dj, 1);
        Vector4d Ldcs;
        if (ties_method == "efron") {
            Ldcs << Lls1(j, 0);
            for (int i = 0; i < dj; i++) {  //  adds in the efron approximation terms
                Ldm.row(i) = (-static_cast<double>(i) / static_cast<double>(dj)) *Ldcs.array();
            }
        }
        Ldm.col(0) = Ldm.col(0).array() + Rs1;
        //  Calculates the left-hand side terms
        MatrixXd temp1 = R.block(RiskFail(j, 0), 0, dj, 1).array().log();
        double Ld1 = (temp1.array().isFinite()).select(temp1, 0).sum();
        //  calculates the right-hand side terms
        temp1 = Ldm.col(0).array().log();
        Rs1 = (temp1.array().isFinite()).select(temp1, 0).sum();
        //
        Ll[0] += Ld1 - Rs1;
    }
    return;
}

//' Utility function to perform null model equivalent of Calculate_Sides with strata
//'
//' \code{Calculate_Null_Sides_Strata} Called to update repeated sum calculations, Uses list of event rows, Performs calculation of counts in each group
//' @inheritParams CPP_template
//'
//' @return Updates matrices in place: risk storage matrices
//' @noRd
//'
void Calculate_Null_Sides_Strata(const IntegerMatrix& RiskFail, const vector<vector<vector<int> > >& RiskPairs_Strata, const int& ntime, const MatrixXd& R, MatrixXd& Rls1, MatrixXd& Lls1, NumericVector& Strata_vals, const int& nthreads) {
    #ifdef _OPENMP
    #pragma omp parallel for schedule(dynamic) num_threads(nthreads) collapse(2)
    #endif
    for (int j = 0; j < ntime; j++) {
        for (int s_ij = 0; s_ij < Strata_vals.size(); s_ij++) {
            double Rs1 = 0;
            //
            //
            vector<int> InGroup = RiskPairs_Strata[j][s_ij];
            //  now has the grouping pairs
            if (RiskFail(j, 2*s_ij + 1)> - 1) {
                int dj = RiskFail(j, 2*s_ij + 1)-RiskFail(j, 2*s_ij + 0) + 1;
                for (vector<double>::size_type i = 0; i < InGroup.size() - 1; i = i + 2) {
                    //
                    Rs1 += InGroup[i + 1]-InGroup[i] + 1;
                }  //  precalculates the sums of risk groups
                //  only assigns values once
                Rls1(j, s_ij) = Rs1;
                Lls1(j, s_ij) = dj;
            }
        }
    }
    return;
}




//' Utility function to perform null model equivalent of Calc_LogLik
//'
//' \code{Calc_Null_LogLik_Strata} Called to update log-likelihoods, Uses list of event rows and repeated sums, Sums the log-likelihood contribution from each event time
//' @inheritParams CPP_template
//'
//' @return Updates matrices in place: Log-likelihood vectors/matrix
//' @noRd
//'
void Calc_Null_LogLik_Strata(const int& nthreads, const IntegerMatrix& RiskFail, const vector<vector<vector<int> > >& RiskPairs_Strata, const int& ntime, const MatrixXd& R, const MatrixXd& Rls1, const MatrixXd& Lls1, NumericVector& Strata_vals, vector<double>& Ll, string ties_method) {
    #ifdef _OPENMP
    #pragma omp declare reduction(vec_double_plus : vector<double> : \
        transform(omp_out.begin(), omp_out.end(), omp_in.begin(), omp_out.begin(), plus<double>())) \
        initializer(omp_priv = omp_orig)
    #pragma omp parallel for schedule(dynamic) num_threads(nthreads) reduction(vec_double_plus:Ll) collapse(2)
    #endif
    for (int s_ij = 0; s_ij < Strata_vals.size(); s_ij++) {
        for (int j = 0; j < ntime; j++) {
            double Rs1 = Rls1(j, s_ij);
            int dj = RiskFail(j, 2*s_ij + 1)-RiskFail(j, 2*s_ij + 0) + 1;
            if (RiskFail(j, 2*s_ij + 1)> - 1) {
                //
                //
                MatrixXd Ldm = MatrixXd::Zero(dj, 1);
                double Ldcs = 0.0;
                if (ties_method == "efron") {
                    Ldcs = Lls1(j, s_ij);
                    for (int i = 0; i < dj; i++) {  //  adds in the efron approximation terms
                        Ldm(i, 0) = (-static_cast<double>(i) / static_cast<double>(dj)) *Ldcs;
                    }
                }
                Ldm.col(0) = Ldm.col(0).array() + Rs1;
                //  Calculates the left-hand side terms
                //  calculates the right-hand side terms
                MatrixXd temp1 = Ldm.col(0).array().log();
                Rs1 = (temp1.array().isFinite()).select(temp1, 0).sum();
                //
                Ll[0] += 0.0 - Rs1;
            }
        }
    }
    return;
}

//' Fills out recursive vectors for matched case-control logistic regression
//'
//' \code{Calculate_Recursive} Called to update the recursive vectors, uses model_bool list to select which vectors to update.
//'
//' @return Updates matrices in place: risk storage matrices
//' @noRd
//'
void Calculate_Recursive(List& model_bool, const int& group_num, const IntegerMatrix& RiskFail, const vector<vector<int> >& RiskPairs, const int& totalnum, const int& ntime, const MatrixXd& R, const MatrixXd& Rd, const MatrixXd& Rdd, vector<vector<double> >& Recur_Base, vector<vector<vector<double> > >& Recur_First, vector<vector<vector<double> > >& Recur_Second, const int& nthreads, const IntegerVector& KeepConstant) {
    int reqrdnum = 1;
    if (!model_bool["null"]) {
        reqrdnum = totalnum - sum(KeepConstant);
    }
    double cond_thres = model_bool["cond_thres"];
    //  We need B vectors
    //  start with the basic B matrix vector
    #ifdef _OPENMP
    #pragma omp parallel for schedule(dynamic) num_threads(nthreads)
    #endif
    for (int group_ij = 0; group_ij < group_num; group_ij++) {
        //  we start by getting a vector of risks
        vector<double> risk_list;
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
                    r_sum += risk_list[i];
                    Recur_Base[group_ij][i] = r_sum;
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
                        //  the filled value is either an edge case, B(m,n) = rn*B(m-1, n-1), or the full case
                        if (j_index == 0) {
                            //  edge case
                            Recur_Base[group_ij][recur_index] = risk_list[risk_index] * Recur_Base[group_ij][t1];
                        } else {
                            Recur_Base[group_ij][recur_index] = Recur_Base[group_ij][t0] + risk_list[risk_index] * Recur_Base[group_ij][t1];
                        }
                    }
                }
            }
        }
    }
    if (!model_bool["single"]) {
        //  next we want the first derivative B matrix vector
        #ifdef _OPENMP
        #pragma omp parallel for schedule(dynamic) num_threads(nthreads) collapse(2)
        #endif
        for (int group_ij = 0; group_ij < group_num; group_ij++) {
            for (int der_ij = 0; der_ij < reqrdnum; der_ij++) {
                //  we start by getting a vector of risks
                vector<double> risk_list;
                vector<double> riskd_list;
                vector<int> InGroup = RiskPairs[group_ij];
                if (InGroup.size() > 0) {
                    int dj = RiskFail(group_ij, 1)-RiskFail(group_ij, 0) + 1;
                    if (dj <= cond_thres) {
                        for (vector<double>::size_type i = 0; i < InGroup.size() - 1; i = i+2) {
                            int i0 = InGroup[i] - 1;
                            int i1 = InGroup[i + 1] - 1;
                            for (int i_inter = i0; i_inter <= i1; i_inter ++) {
                                risk_list.push_back(R(i_inter, 0));
                                riskd_list.push_back(Rd(i_inter, der_ij));
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
                            r_sum += riskd_list[i];
                            Recur_First[group_ij][der_ij][i] = r_sum;
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
                                    Recur_First[group_ij][der_ij][recur_index] = risk_list[risk_index] * Recur_First[group_ij][der_ij][t1] + riskd_list[risk_index] * Recur_Base[group_ij][t1];
                                } else {
                                    Recur_First[group_ij][der_ij][recur_index] = Recur_First[group_ij][der_ij][t0] + risk_list[risk_index] * Recur_First[group_ij][der_ij][t1] + riskd_list[risk_index] * Recur_Base[group_ij][t1];
                                }
                            }
                        }
                    }
                }
            }
        }
        if (!model_bool["gradient"]) {
            //  finally we want the second derivative B matrix vector
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
                    vector<double> riskdd_list;
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
                                    riskdd_list.push_back(Rdd(i_inter, der_ijk));
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
                                r_sum += riskdd_list[i];
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
                                        Recur_Second[group_ij][der_ijk][recur_index] = risk_list[risk_index] * Recur_Second[group_ij][der_ijk][t1] + riskd0_list[risk_index] * Recur_First[group_ij][der_jk][t1] + riskd1_list[risk_index] * Recur_First[group_ij][der_ij][t1] + riskdd_list[risk_index] * Recur_Base[group_ij][t1];
                                    } else {
                                        Recur_Second[group_ij][der_ijk][recur_index] = Recur_Second[group_ij][der_ijk][t0] + risk_list[risk_index] * Recur_Second[group_ij][der_ijk][t1] + riskd0_list[risk_index] * Recur_First[group_ij][der_jk][t1] + riskd1_list[risk_index] * Recur_First[group_ij][der_ij][t1] + riskdd_list[risk_index] * Recur_Base[group_ij][t1];
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }
    return;
}

//' Fills out recursive vectors for matched case-control logistic regression
//'
//' \code{Calc_Recur_LogLik} Called to update the recursive vectors, uses model_bool list to select which vectors to update.
//'
//' @return Updates matrices in place: risk storage matrices
//' @noRd
//'
void Calc_Recur_LogLik(List& model_bool, const int& group_num, const IntegerMatrix& RiskFail, const vector<vector<int> >& RiskPairs, const int& totalnum, const int& ntime, const MatrixXd& R, const MatrixXd& Rd, const MatrixXd& Rdd, const MatrixXd& RdR, const MatrixXd& RddR, double& dev, vector<double>& Ll, vector<double>& Lld, vector<double>& Lldd, vector<vector<double> >& Recur_Base, vector<vector<vector<double> > >& Recur_First, vector<vector<vector<double> > >& Recur_Second, vector<double>& strata_odds, const int& nthreads, const IntegerVector& KeepConstant, vector<int>& strata_cond, vector<double>& LldOdds, vector<double>& LlddOdds, vector<double>& LlddOddsBeta) {
    int reqrdnum = 1;
    if (!model_bool["null"]) {
        reqrdnum = totalnum - sum(KeepConstant);
    }
    fill(Ll.begin(), Ll.end(), 0.0);
    if (!model_bool["single"]) {
        fill(Lld.begin(), Lld.end(), 0.0);
        fill(LldOdds.begin(), LldOdds.end(), 0.0);
        if (!model_bool["gradient"]) {
            fill(Lldd.begin(), Lldd.end(), 0.0);
            fill(LlddOdds.begin(), LlddOdds.end(), 0.0);
            fill(LlddOddsBeta.begin(), LlddOddsBeta.end(), 0.0);
        }
    }
    #ifdef _OPENMP
    #pragma omp declare reduction(vec_double_plus : vector<double> : \
        transform(omp_out.begin(), omp_out.end(), omp_in.begin(), omp_out.begin(), plus<double>())) \
        initializer(omp_priv = omp_orig)
    #endif
    double cond_thres = model_bool["cond_thres"];
    //  we need to get the repeated values for unconditional likelihood calculation
    if (model_bool["single"]) {
        //  now we can calculate the loglikelihoods
        #ifdef _OPENMP
        #pragma omp parallel for schedule(dynamic) num_threads(nthreads) \
            reduction(vec_double_plus:Ll)  reduction(+:dev)
        #endif
        for (int group_ij = 0; group_ij < group_num; group_ij++) {
            //
            int recur_index = static_cast<int>(Recur_Base[group_ij].size() - 1);
            int dj = RiskFail(group_ij, 1)-RiskFail(group_ij, 0) + 1;
            if (recur_index > -1) {
                //
                double Ld1 = 0.0;
                //
                MatrixXd temp1 = MatrixXd::Zero(dj, 1);
                temp1 = R.block(RiskFail(group_ij, 0), 0, dj, 1).array().log();
                Ld1 = (temp1.array().isFinite()).select(temp1, 0).sum();
                double Rs1 = 0.0;
                if (dj <= cond_thres) {
                    //  calculates the right-hand side terms
                    double b_0 = Recur_Base[group_ij][recur_index];
                    //
                    Rs1 = log(b_0);
                    //
                } else {
                    Ld1 = strata_odds[group_ij]*dj + Ld1;
                    vector<int> InGroup = RiskPairs[group_ij];
                    for (vector<double>::size_type i = 0; i < InGroup.size() - 1; i = i+2) {
                        Rs1 += (1.0 + exp(strata_odds[group_ij]) * R.block(InGroup[i] - 1, 0, InGroup[i + 1]-InGroup[i] + 1, 1).array()).log().sum();
                    }
                }
                dev += -2*(Ld1 - Rs1);
                Ll[0] += Ld1 - Rs1;
            }
        }
    } else if (model_bool["gradient"]) {
        //  now we can calculate the loglikelihoods first derivative
        #ifdef _OPENMP
        #pragma omp parallel for schedule(dynamic) num_threads(nthreads) \
            reduction(vec_double_plus:Ll, Lld)  reduction(+:dev) collapse(2)
        #endif
        for (int group_ij = 0; group_ij < group_num; group_ij++) {
            for (int der_ij = 0; der_ij < reqrdnum; der_ij++) {
                //
                int recur_index = static_cast<int>(Recur_Base[group_ij].size() - 1);
                int dj = RiskFail(group_ij, 1)-RiskFail(group_ij, 0) + 1;
                if (recur_index > -1) {
                    //
                    double Ld1 = 0.0;
                    double Ld2 = 0.0;
                    //
                    MatrixXd temp1 = MatrixXd::Zero(dj, 1);
                    MatrixXd temp2 = MatrixXd::Zero(dj, 1);
                    temp1 = R.block(RiskFail(group_ij, 0), 0, dj, 1).array().log();
                    Ld1 = (temp1.array().isFinite()).select(temp1, 0).sum();
                    temp1 = RdR.block(RiskFail(group_ij, 0), der_ij, dj, 1).array();
                    Ld2 = (temp1.array().isFinite()).select(temp1, 0).sum();
                    //  calculates the right-hand side terms
                    double Rs1 = 0.0;
                    double Rs2 = 0.0;
                    if (dj <= cond_thres) {
                        double b_0 = Recur_Base[group_ij][recur_index];
                        double b_1 = Recur_First[group_ij][der_ij][recur_index];
                        //
                        Rs1 = log(b_0);
                        Rs2 = b_1 / b_0;
                    } else {
                        vector<int> InGroup = RiskPairs[group_ij];
                        Ld1 = strata_odds[group_ij]*dj + Ld1;
                        for (vector<double>::size_type i = 0; i < InGroup.size() - 1; i = i+2) {
                            Rs1 += (1.0 + exp(strata_odds[group_ij]) * R.block(InGroup[i] - 1, 0, InGroup[i + 1]-InGroup[i] + 1, 1).array()).log().sum();
                            Rs2 += (Rd.block(InGroup[i] - 1, der_ij, InGroup[i + 1]-InGroup[i] + 1, 1).array() * (1.0 + exp(strata_odds[group_ij]) * R.block(InGroup[i] - 1, 0, InGroup[i + 1]-InGroup[i] + 1, 1).array()).pow(-1).array()).sum();
                        }
                        Rs2 *= exp(strata_odds[group_ij]);
                    }
                    if (der_ij == 0) {
                        dev += -2*(Ld1 - Rs1);
                    }
                    //
                    Ll[der_ij] += Ld1 - Rs1;
                    Lld[der_ij] += Ld2 - Rs2;
                }
            }
        }
        #ifdef _OPENMP
        #pragma omp parallel for schedule(dynamic) num_threads(nthreads) \
            reduction(vec_double_plus:LldOdds)
        #endif
        for (int group_ij = 0; group_ij < group_num; group_ij++) {
            if (strata_cond[group_ij] == 0) {
                vector<int>::iterator it_end = strata_cond.begin();
                advance(it_end, group_ij);
                int group_jk = group_ij - reduce(strata_cond.begin(), it_end);
                //
                int dj = RiskFail(group_ij, 1)-RiskFail(group_ij, 0) + 1;
                vector<int> InGroup = RiskPairs[group_ij];
                double Rs1 = 0.0;
                int num_row = 0;
                for (vector<double>::size_type i = 0; i < InGroup.size() - 1; i = i+2) {
                    Rs1 += (1.0 + exp(strata_odds[group_ij]) * R.block(InGroup[i] - 1, 0, InGroup[i + 1]-InGroup[i] + 1, 1).array()).pow(-1).array().sum();
                    num_row += InGroup[i + 1]-InGroup[i] + 1;
                }
                double Ls1 = dj - num_row;
                LldOdds[group_jk] += Ls1 + Rs1;
            }
        }
    } else {
        //  now we can calculate the loglikelihoods second derivatives
        #ifdef _OPENMP
        #pragma omp parallel for schedule(dynamic) num_threads(nthreads) \
            reduction(vec_double_plus:Ll, Lld, Lldd) reduction(+:dev) collapse(2)
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
                    MatrixXd Ld = MatrixXd::Zero(dj, 4);
                    Ld << R.block(RiskFail(group_ij, 0), 0, dj, 1), RdR.block(RiskFail(group_ij, 0), der_ij, dj, 1), RdR.block(RiskFail(group_ij, 0), der_jk, dj, 1), RddR.block(RiskFail(group_ij, 0), der_ijk, dj, 1);  //  rows with events
                    //
                    double Ld1 = 0.0;
                    double Ld2 = 0.0;
                    double Ld3 = 0.0;
                    //
                    MatrixXd temp1 = MatrixXd::Zero(dj, 1);
                    MatrixXd temp2 = MatrixXd::Zero(dj, 1);
                    if (der_ij == der_jk) {
                        temp1 = R.block(RiskFail(group_ij, 0), 0, dj, 1).array().log();
                        Ld1 = (temp1.array().isFinite()).select(temp1, 0).sum();
                    }
                    temp1 = RdR.block(RiskFail(group_ij, 0), der_ij, dj, 1).array();
                    if (der_ij == der_jk) {
                        Ld2 = (temp1.array().isFinite()).select(temp1, 0).sum();
                    }
                    temp2 = RdR.block(RiskFail(group_ij, 0), der_jk, dj, 1).array();
                    temp1 = RddR.block(RiskFail(group_ij, 0), der_ijk, dj, 1).array() - (temp1.array() * temp2.array());
                    Ld3 = (temp1.array().isFinite()).select(temp1, 0).sum();
                    //  calculates the right-hand side terms
                    double Rs1 = 0.0;
                    double Rs2 = 0.0;
                    double Rs3 = 0.0;
                    if (dj <= cond_thres) {
                        double b_0 = Recur_Base[group_ij][recur_index];
                        double b_1 = Recur_First[group_ij][der_ij][recur_index];
                        double b_2 = Recur_First[group_ij][der_jk][recur_index];
                        double b_3 = Recur_Second[group_ij][der_ijk][recur_index];
                        //
                        if (der_ij == der_jk) {
                            Rs1 = log(b_0);
                            Rs2 = b_1 / b_0;
                        }
                        Rs3 = b_3 / b_0 - b_1 / b_0 * b_2 / b_0;
                    } else {
                        vector<int> InGroup = RiskPairs[group_ij];
                        Ld1 = strata_odds[group_ij]*dj + Ld1;
                        if (der_ij == der_jk) {
                            for (vector<double>::size_type i = 0; i < InGroup.size() - 1; i = i+2) {
                                Rs1 += (1.0 + exp(strata_odds[group_ij]) * R.block(InGroup[i] - 1, 0, InGroup[i + 1]-InGroup[i] + 1, 1).array()).log().sum();
                                Rs2 += (Rd.block(InGroup[i] - 1, der_ij, InGroup[i + 1]-InGroup[i] + 1, 1).array() * (1.0 + exp(strata_odds[group_ij]) * R.block(InGroup[i] - 1, 0, InGroup[i + 1]-InGroup[i] + 1, 1).array()).pow(-1).array()).sum();
                            }
                            Rs2 *= exp(strata_odds[group_ij]);
                        }
                        double Rs3l = 0.0;
                        double Rs3r = 0.0;
                        for (vector<double>::size_type i = 0; i < InGroup.size() - 1; i = i+2) {
                            Rs3r += (Rd.block(InGroup[i] - 1, der_ij, InGroup[i + 1]-InGroup[i] + 1, 1).array() * Rd.block(InGroup[i] - 1, der_jk, InGroup[i + 1]-InGroup[i] + 1, 1).array() * (1.0 + exp(strata_odds[group_ij]) * R.block(InGroup[i] - 1, 0, InGroup[i + 1]-InGroup[i] + 1, 1).array()).pow(-2).array() ).sum();
                            Rs3l += (Rdd.block(InGroup[i] - 1, der_ijk, InGroup[i + 1]-InGroup[i] + 1, 1).array() * (1.0 + exp(strata_odds[group_ij]) * R.block(InGroup[i] - 1, 0, InGroup[i + 1]-InGroup[i] + 1, 1).array()).pow(-1).array()).sum();
                        }
                        Rs3 = exp(strata_odds[group_ij]) * (Rs3l - exp(strata_odds[group_ij])*Rs3r);
                    }
                    //
                    if (der_ij == der_jk) {
                        Ll[der_ij] += Ld1 - Rs1;
                        Lld[der_ij] += Ld2 - Rs2;
                        if (der_ij == 0) {
                            dev += -2*(Ld1 - Rs1);
                        }
                    }
                    Lldd[der_ij*reqrdnum+der_jk] += Ld3 - Rs3;  //  sums the log-likelihood and derivatives
                }
            }
        }
        #ifdef _OPENMP
        #pragma omp parallel for schedule(dynamic) num_threads(nthreads) \
            reduction(vec_double_plus:LldOdds, LlddOdds)
        #endif
        for (int group_ij = 0; group_ij < group_num; group_ij++) {
            if (strata_cond[group_ij] == 0) {
                vector<int>::iterator it_end = strata_cond.begin();
                advance(it_end, group_ij);
                int group_jk = group_ij - reduce(strata_cond.begin(), it_end);
                //
                int dj = RiskFail(group_ij, 1)-RiskFail(group_ij, 0) + 1;
                int num_row = 0;
                vector<int> InGroup = RiskPairs[group_ij];
                double Rs1 = 0.0;
                double Rs2 = 0.0;
                for (vector<double>::size_type i = 0; i < InGroup.size() - 1; i = i+2) {
                    Rs1 += (1.0 + exp(strata_odds[group_ij]) * R.block(InGroup[i] - 1, 0, InGroup[i + 1]-InGroup[i] + 1, 1).array()).pow(-1).array().sum();
                    Rs2 += (R.block(InGroup[i] - 1, 0, InGroup[i + 1]-InGroup[i] + 1, 1).array() * (1.0 + exp(strata_odds[group_ij]) * R.block(InGroup[i] - 1, 0, InGroup[i + 1]-InGroup[i] + 1, 1).array()).pow(-2).array()).sum();
                    num_row += InGroup[i + 1]-InGroup[i] + 1;
                }
                double Ls1 = dj-num_row;
                LldOdds[group_jk] += Ls1 + Rs1;
                LlddOdds[group_jk] -= exp(strata_odds[group_ij])*Rs2;
            }
        }
        #ifdef _OPENMP
        #pragma omp parallel for schedule(dynamic) num_threads(nthreads) \
            reduction(vec_double_plus:LlddOddsBeta) collapse(2)
        #endif
        for (int group_ij = 0; group_ij < group_num; group_ij++) {
            for (int der_ij = 0; der_ij < reqrdnum; der_ij++) {
                if (strata_cond[group_ij] == 0) {
                    vector<int>::iterator it_end = strata_cond.begin();
                    advance(it_end, group_ij);
                    int group_jk = group_ij - reduce(strata_cond.begin(), it_end);
                    //
                    vector<int> InGroup = RiskPairs[group_ij];
                    double Rs2 = 0.0;
                    for (vector<double>::size_type i = 0; i < InGroup.size() - 1; i = i+2) {
                        Rs2 += (Rd.block(InGroup[i] - 1, der_ij, InGroup[i + 1]-InGroup[i] + 1, 1).array() * (1.0 + exp(strata_odds[group_ij]) * R.block(InGroup[i] - 1, 0, InGroup[i + 1]-InGroup[i] + 1, 1).array()).pow(-2).array()).sum();
                    }
                    LlddOddsBeta[group_jk*reqrdnum + der_ij] -= exp(strata_odds[group_ij])*Rs2;
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
            Lldd[jk*reqrdnum+ij] = Lldd[ij*reqrdnum+jk];
        }
    }
    double LogLik = 0;
    for (int i = 0; i < reqrdnum; i++) {
        if (Ll[i] != 0) {
            LogLik = Ll[i];
            break;
        }
    }
    fill(Ll.begin(), Ll.end(), LogLik);
    return;
}

//' Utility function to calculate Logistic Log-Likelihood and derivatives
//'
//' \code{Calc_LogLik_Logist} Called to update log-likelihoods, Uses probability matrices Sums the log-likelihood contribution from each row
//' @inheritParams CPP_template
//'
//' @return Updates matrices in place: Log-likelihood vectors/matrix
//' @noRd
//'
void Calc_LogLik_Logist(List& model_bool, const int& nthreads, const int& totalnum, const Ref<const MatrixXd>& CountEvent, const MatrixXd& P, const MatrixXd& Pnot, const MatrixXd& Pd, const MatrixXd& Pdd, const MatrixXd& PdP, const MatrixXd& PnotdP, const MatrixXd& PddP, const MatrixXd& PnotddP, vector<double>& Ll, vector<double>& Lld, vector<double>& Lldd, const IntegerVector& KeepConstant) {
    int reqrdnum = totalnum - sum(KeepConstant);
    fill(Ll.begin(), Ll.end(), 0.0);
    if (!model_bool["single"]) {
        fill(Lld.begin(), Lld.end(), 0.0);
        if (!model_bool["gradient"]) {
            fill(Lldd.begin(), Lldd.end(), 0.0);
        }
    }
    MatrixXd temp(Pd.rows(), 2);
    temp.col(0) = CountEvent.col(1).array() - CountEvent.col(0).array();  // N_i - y_i
    temp.col(1) = CountEvent.col(0).array() * P.col(0).array().log().array() + temp.col(0).array() * Pnot.col(0).array().log().array();
    fill(Ll.begin(), Ll.end(), (temp.col(1).array().isFinite()).select(temp.col(1), 0).sum());
    if (!model_bool["single"]) {
        // first derivatives
        #ifdef _OPENMP
        #pragma omp parallel for schedule(dynamic) num_threads(nthreads)
        #endif
        for (int ij = 0; ij < reqrdnum; ij++) {  //  performs log-likelihood calculations for every first derivative
            double dl = (CountEvent.col(0).array() * PdP.col(ij).array() - temp.col(0).array() * PnotdP.col(ij).array()).array().sum();
            Lld[ij] = dl;
        }
        if (!model_bool["gradient"]) {
            // second derivatives
            #ifdef _OPENMP
            #pragma omp parallel for schedule(dynamic) num_threads(nthreads)
            #endif
            for (int ijk = 0; ijk < reqrdnum*(reqrdnum + 1)/2; ijk++) {  //  performs log-likelihood calculations for every derivative combination and risk group
                int ij = 0;
                int jk = ijk;
                while (jk > ij) {
                    ij++;
                    jk -= ij;
                }
                double ddl = (CountEvent.col(0).array() * (PddP.col(ijk).array() - PdP.col(ij).array() * PdP.col(jk).array()).array() - temp.col(0).array() * (PnotddP.col(ijk).array() + PnotdP.col(ij).array() * PnotdP.col(jk).array()).array() ).array().sum();
                Lldd[ij*reqrdnum+jk] = ddl;
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
                Lldd[jk*reqrdnum+ij] = Lldd[ij*reqrdnum+jk];
            }
        }
    }
    return;
}
