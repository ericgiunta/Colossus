#include <RcppEigen.h>
#include <RcppParallel.h>
#include <omp.h>
#include "Main_Functions.h"
#include "Calc_Repeated.h"
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

//' Primary Cox PH regression
//' \code{LogLik_Cox_PH} Performs the calls to calculation functions, Structures the Cox PH regression, With verbose option prints out time stamps and intermediate sums of terms and derivatives
//'
//' @param     Term_n    Term numbers
//' @param     tform    subterm types
//' @param     a_n    starting values
//' @param     x_all    covariate matrix
//' @param     dfc    covariate column numbers
//' @param     fir    first term number
//' @param     der_iden    subterm number for derivative tests
//' @param     modelform    model string
//' @param     lr    learning rate for newton step toward 0 derivative
//' @param     maxiter    maximum number of iterations
//' @param     halfmax    maximum number of half steps
//' @param     epsilon    minimum acceptable maximum parameter change
//' @param     dbeta_cap    learning rate for newton step toward 0 log-likelihood
//' @param     abs_max    Maximum allowed parameter change
//' @param     dose_abs_max    Maximum allowed threshold parameter change
//' @param     deriv_epsilon    threshold for near-zero derivative
//' @param     df_groups    matrix with time and event information
//' @param     tu    event times
//' @param     change_all    boolean if every parameter is being updated
//' @param     verbose    verbosity boolean
//' @param     debugging    debugging boolean
//' @param     KeepConstant    vector identifying constant parameters
//' @param     term_tot    total number of terms
//' @param     ties_method    ties method
//'
//' @return List of results: Log-likelihood of optimum, first derivative of log-likelihood, second derivative matrix, parameter list, standard deviation estimate, AIC, model information
// [[Rcpp::export]]
List LogLik_Cox_PH( IntegerVector Term_n, StringVector tform, NumericVector a_n,NumericMatrix x_all,IntegerVector dfc,int fir, int der_iden,string modelform, double lr, int maxiter, int halfmax, double epsilon, double dbeta_cap, double abs_max,double dose_abs_max, double deriv_epsilon, NumericMatrix df_groups, NumericVector tu ,bool change_all, bool verbose, bool debugging, IntegerVector KeepConstant, int term_tot, string ties_method){
    ;
    //
    List temp_list = List::create(_["Status"]="FAILED"); //used as a dummy return value for code checking
    if (verbose){
        Rcout << "START_NEW" << endl;
    }
    time_point<system_clock> start_point, end_point;
    start_point = system_clock::now();
    auto start = time_point_cast<microseconds>(start_point).time_since_epoch().count();
    end_point = system_clock::now();
    auto ending = time_point_cast<microseconds>(end_point).time_since_epoch().count(); //The time duration is tracked
    //
    auto gibtime = system_clock::to_time_t(system_clock::now());
    if (verbose){
        Rcout << ctime(&gibtime) << endl;
    }
    //
    // Time durations are measured from this point on in microseconds
    //
    // df0: covariate data
    // ntime: number of event times for Cox PH
    // totalnum: number of terms used
    //
    const Map<MatrixXd> df0(as<Map<MatrixXd> >(x_all));
    int ntime = tu.size();
    //
    int totalnum = Term_n.size();
    //
    if (verbose){
        Rcout << "Term checked ";
        for (int ij=0;ij<totalnum;ij++){
            Rcout << Term_n[ij] << " ";
        }
        Rcout << " " << endl;
    }
    //
    // cout.precision: controls the number of significant digits printed
    // nthreads: number of threads used for parallel operations
    //
    Rcout.precision(7); //forces higher precision numbers printed to terminal
    int nthreads = Eigen::nbThreads()-1; //stores how many threads are allocated
    //
    //
    // Lld_worst: stores the highest magnitude log-likelihood derivative
    // totem: number of rows needed
    //
    //
    double Lld_worst = 0.0; //stores derivative value used to determine if every parameter is near convergence
    //
    if (verbose){
        end_point = system_clock::now();
        ending = time_point_cast<microseconds>(end_point).time_since_epoch().count();
        Rcout <<"df99,"<<(ending-start)<<",Starting"<<endl;
        gibtime = system_clock::to_time_t(system_clock::now());
        Rcout << ctime(&gibtime) << endl;
    }
    // ---------------------------------------------
    // To Start, needs to seperate the derivative terms
    // ---------------------------------------------
    //
    Map<VectorXd> beta_0(as<Map<VectorXd> >(a_n));
    MatrixXd T0 = MatrixXd::Zero(df0.rows(), totalnum); //preallocates matrix for Term column
    MatrixXd Td0 = MatrixXd::Zero(df0.rows(), totalnum); //preallocates matrix for Term derivative columns
    MatrixXd Tdd0 = MatrixXd::Zero(df0.rows(), totalnum*(totalnum+1)/2); //preallocates matrix for Term second derivative columns
    //
    MatrixXd Te = MatrixXd::Zero(df0.rows(), 1); //preallocates matrix for column terms used for temporary storage
    MatrixXd R = MatrixXd::Zero(df0.rows(), 1); //preallocates matrix for Risks
    MatrixXd Rd = MatrixXd::Zero(df0.rows(), totalnum); //preallocates matrix for Risk derivatives
    MatrixXd Rdd = MatrixXd::Zero(df0.rows(), totalnum*(totalnum+1)/2); //preallocates matrix for Risk second derivatives
    //
    MatrixXd Dose = MatrixXd::Constant(df0.rows(),term_tot,0.0); //Matrix of the total dose term values
    MatrixXd nonDose = MatrixXd::Constant(df0.rows(),term_tot,1.0); //Matrix of the total non-dose term values
    MatrixXd nonDose_LIN = MatrixXd::Constant(df0.rows(),term_tot,0.0); //matrix of Linear subterm values
    MatrixXd nonDose_PLIN = MatrixXd::Constant(df0.rows(),term_tot,1.0); //matrix of Loglinear subterm values
    MatrixXd nonDose_LOGLIN = MatrixXd::Constant(df0.rows(),term_tot,1.0); //matrix of Product linear subterm values
    MatrixXd TTerm=MatrixXd::Zero(Dose.rows(),Dose.cols()); //matrix of term values
    double dint = dose_abs_max; //The amount of change used to calculate derivatives in threshold paramters
    //
    if (verbose){
        Rcout << "starting subterms " << term_tot << endl;
    }
    // Calculates the subterm and term values
    Make_Subterms( totalnum, Term_n, tform, dfc, fir, T0, Td0, Tdd0, Dose, nonDose, TTerm,  nonDose_LIN, nonDose_PLIN, nonDose_LOGLIN ,beta_0, df0,dint,nthreads, debugging);
    // ---------------------------------------------------------
    // Prints off a series of calculations to check at what point values are changing
    // ---------------------------------------------------------
    if (verbose){
        Rcout << "values checked ";
        for (int ijk=0;ijk<totalnum;ijk++){
            Rcout << beta_0[ijk] << " ";
        }
        Rcout << " " << endl;
        Rcout << "sums checked ";
        for (int ijk=0;ijk<totalnum;ijk++){
            Rcout << T0.col(ijk).sum() << " ";
        }
        Rcout << " " << endl;
        Rcout << "derivs checked ";
        for (int ijk=0;ijk<totalnum;ijk++){
            Rcout << Td0.col(ijk).sum() << " ";
        }
        Rcout << " " << endl;
        Rcout << "second derivs checked ";
        for (int ijk=0;ijk<totalnum;ijk++){
            Rcout << Tdd0.col(ijk*(ijk+1)/2+ijk).sum() << " ";
        }
        Rcout << " " << endl;
        Rcout << "dose checked ";
        for (int ijk=0;ijk<term_tot;ijk++){
            Rcout << Dose.col(ijk).array().sum() << " ";
        }
        Rcout << " " << endl;
        Rcout << "non-dose checked ";
        for (int ijk=0;ijk<term_tot;ijk++){
            Rcout << nonDose.col(ijk).array().sum() << " ";
        }
        Rcout << " " << endl;
        Rcout << "LIN_non-dose checked ";
        for (int ijk=0;ijk<term_tot;ijk++){
            Rcout << nonDose_LIN.col(ijk).array().sum() << " ";
        }
        Rcout << " " << endl;
        Rcout << "PLIN_non-dose checked ";
        for (int ijk=0;ijk<term_tot;ijk++){
            Rcout << nonDose_PLIN.col(ijk).array().sum() << " ";
        }
        Rcout << " " << endl;
        Rcout << "LOGLIN_non-dose checked ";
        for (int ijk=0;ijk<term_tot;ijk++){
            Rcout << nonDose_LOGLIN.col(ijk).array().sum() << " ";
        }
        Rcout << " " << endl;
    }
    //
    MatrixXd RdR = MatrixXd::Zero(df0.rows(), totalnum); //preallocates matrix for Risk to derivative ratios
    MatrixXd RddR = MatrixXd::Zero(df0.rows(), totalnum*(totalnum+1)/2); //preallocates matrix for Risk to second derivative ratios
    //
    if (verbose){
        end_point = system_clock::now();
        ending = time_point_cast<microseconds>(end_point).time_since_epoch().count();
        Rcout <<"df99,"<<(ending-start)<<",Prep_Terms"<<endl;
        gibtime = system_clock::to_time_t(system_clock::now());
        Rcout << ctime(&gibtime) << endl;
    }
    //
    // Calculates the risk for each row
    Make_Risks(modelform, tform, Term_n, totalnum, fir, T0, Td0, Tdd0, Te, R, Rd, Rdd, Dose, nonDose, TTerm, nonDose_LIN, nonDose_PLIN, nonDose_LOGLIN, RdR, RddR, nthreads, debugging);
    //
    // Removes infinite values
    RdR = (RdR.array().isFinite()).select(RdR,0);
    RddR = (RddR.array().isFinite()).select(RddR,0);
    //
    if (R.minCoeff()<=0){
        //
        // The risk cannot be negative as the code currently is
        //
        for (int ijk=0;ijk<totalnum;ijk++){
            if (T0.col(ijk).minCoeff()<=0){
                Rcout << ijk << " had a non-positive term" << endl;
            }
        }
        for (int ijk=0;ijk<term_tot;ijk++){
            if (nonDose_LIN.col(ijk).minCoeff()<=0){
                Rcout << ijk << " had a non-positive Lin term" << endl;
            }
            if (nonDose_LOGLIN.col(ijk).minCoeff()<=0){
                Rcout << ijk << " had a non-positive loglin term" << endl;
            }
            if (nonDose_PLIN.col(ijk).minCoeff()<=0){
                Rcout << ijk << " had a non-positive plin term" << endl;
            }
            if (Dose.col(ijk).minCoeff()<=0){
                Rcout << ijk << " had a non-positive dose term" << endl;
            }
            if (nonDose.col(ijk).minCoeff()<=0){
                Rcout << ijk << " had a non-positive nondose term" << endl;
            }
            if (TTerm.col(ijk).minCoeff()<=0){
                Rcout << ijk << " had a non-positive total term" << endl;
            }
        }
        Rcout << R.sum() << endl;
        Rcout << "A non-positive risk was detected: " << R.minCoeff() << endl;
        return temp_list;
    }
    //
    if (verbose){
        Rcout << "risk checked ";
        for (int ijk=0;ijk<totalnum;ijk++){
            Rcout << R.col(0).sum() << " ";
        }
        Rcout << " " << endl;
        Rcout << "risk1 checked ";
        for (int ijk=0;ijk<totalnum;ijk++){
            Rcout << Rd.col(ijk).sum() << " ";
        }
        Rcout << " " << endl;
        Rcout << "risk2 checked ";
        for (int ijk=0;ijk<totalnum;ijk++){
            Rcout << Rdd.col(ijk*(ijk+1)/2+ijk).sum() << " ";
        }
        Rcout << " " << endl;
        //
        end_point = system_clock::now();
        ending = time_point_cast<microseconds>(end_point).time_since_epoch().count();
        Rcout <<"df100 "<<(ending-start)<<" "<<0<<" "<<0<<" "<<-1<<",Prep_R"<<endl;
        gibtime = system_clock::to_time_t(system_clock::now());
        Rcout << ctime(&gibtime) << endl;
    }
    //
    // -------------------------------------------------------------------------------------------
    //
    vector<string>  RiskGroup(ntime); //vector of strings detailing the rows
    IntegerMatrix RiskFail(ntime,2); //vector giving the event rows
    //
    const Map<MatrixXd> df_m(as<Map<MatrixXd> >(df_groups));
    // Creates matrices used to identify the event risk groups
    Make_Groups( ntime, df_m, RiskFail, RiskGroup, tu, nthreads, debugging);
    //
    if (verbose){
        end_point = system_clock::now();
        ending = time_point_cast<microseconds>(end_point).time_since_epoch().count();
        Rcout <<"df100 "<<(ending-start)<<" "<<0<<" "<<0<<" "<<-1<<",Prep_List"<<endl;
        gibtime = system_clock::to_time_t(system_clock::now());
        Rcout << ctime(&gibtime) << endl;
    }
    //
    // --------------------------
    // now a vector exists with row locations
    // --------------------------
    MatrixXd Rls1 =MatrixXd::Zero(ntime, 1); //precomputes a series of sums used frequently in the log-liklihood calculations
    MatrixXd Rls2 =MatrixXd::Zero(ntime, totalnum); //Many are repeated due to the same risk groups and derivatives being used at mulitple points
    MatrixXd Rls3 =MatrixXd::Zero(ntime, totalnum*(totalnum+1)/2); //Sum and its derivatives are precomputed
    MatrixXd Lls1 =MatrixXd::Zero(ntime, 1); //The log-likelihood calculation has a Right and Left sum used
    MatrixXd Lls2 =MatrixXd::Zero(ntime, totalnum);
    MatrixXd Lls3 =MatrixXd::Zero(ntime, totalnum*(totalnum+1)/2);
    //The log-likelihood is calculated in parallel over the risk groups
    vector<double> Ll(totalnum,0.0); //Log-likelihood values
    vector<double> Lld(totalnum,0.0); //Log-likelihood derivative values
    vector<double> Lldd(pow(totalnum,2),0.0);//The second derivative matrix has room for every combination, but only the lower triangle is calculated initially
    //
    // Calculates the side sum terms used
    Calculate_Sides( RiskFail, RiskGroup, totalnum, ntime, R, Rd, Rdd, Rls1, Rls2, Rls3, Lls1, Lls2, Lls3,nthreads, debugging);
    //
    //
    if (verbose){
        Rcout << "riskr checked ";
        for (int ijk=0;ijk<totalnum;ijk++){
            Rcout << Rls1.col(0).sum() << " ";
        }
        Rcout << " " << endl;
        Rcout << "risk1r checked ";
        for (int ijk=0;ijk<totalnum;ijk++){
            Rcout << Rls2.col(ijk).sum() << " ";
        }
        Rcout << " " << endl;
        Rcout << "risk2r checked ";
        for (int ijk=0;ijk<totalnum;ijk++){
            Rcout << Rls3.col(ijk*(ijk+1)/2+ijk).sum() << " ";
        }
        Rcout << " " << endl;
        //
        Rcout << "riskl checked ";
        for (int ijk=0;ijk<totalnum;ijk++){
            Rcout << Lls1.col(0).sum() << " ";
        }
        Rcout << " " << endl;
        Rcout << "risk1l checked ";
        for (int ijk=0;ijk<totalnum;ijk++){
            Rcout << Lls2.col(ijk).sum() << " ";
        }
        Rcout << " " << endl;
        Rcout << "risk2l checked ";
        for (int ijk=0;ijk<totalnum;ijk++){
            Rcout << Lls3.col(ijk*(ijk+1)/2+ijk).sum() << " ";
        }
        Rcout << " " << endl;
    }
    //
    //
    // Calculates log-likelihood
    Calc_LogLik( nthreads, RiskFail, RiskGroup, totalnum, ntime, R, Rd, Rdd,RdR,RddR, Rls1, Rls2, Rls3, Lls1, Lls2, Lls3, Ll, Lld, Lldd, debugging, ties_method);
    //
    vector <double> Ll_comp(2,Ll[0]); //vector to compare values
    //
    if (verbose){
        end_point = system_clock::now();
        ending = time_point_cast<microseconds>(end_point).time_since_epoch().count();
        Rcout <<"df100 "<<(ending-start)<<" "<<0<<" "<<0<<" "<<0<<",Calc"<<endl;//prints the time
        gibtime = system_clock::to_time_t(system_clock::now());
        Rcout << ctime(&gibtime) << endl;
        Rcout << "df101 ";//prints the log-likelihoods
        for (int ij=0;ij<totalnum;ij++){
            Rcout << Ll[ij] << " ";
        }
        Rcout << " " << endl;
        Rcout << "df102 ";//prints the first derivatives
        for (int ij=0;ij<totalnum;ij++){
            Rcout << Lld[ij] << " ";
        }
        Rcout << " " << endl;
        Rcout << "df103 ";//prints the second derivatives
        for (int ij=0;ij<totalnum;ij++){
            Rcout << Lldd[ij*totalnum+ij] << " ";
        }
        for (int ij=0;ij<totalnum;ij++){//locates highest magnitude derivative
            if (abs(Lld[ij]) > Lld_worst){
                Lld_worst = abs(Lld[ij]);
            }
        }
        Rcout << " " << endl;
        Rcout << "df104 ";//prints parameter values
        for (int ij=0;ij<totalnum;ij++){
            Rcout << beta_0[ij] << " ";
        }
        Rcout << " " << endl;
        Rcout << "df105 ";
        for (int ij=0;ij<totalnum;ij++){//prints the newton step value for zero derivative
            Rcout << Lld[ij]/Lldd[ij*totalnum+ij] << " ";
        }
        Rcout << " " << endl;
        Rcout << "df106 ";
        for (int ij=0;ij<totalnum;ij++){//prints the newton step value for zero log-likelihood
            Rcout << Ll[ij]/Lld[ij] << " ";
        }
        Rcout << " " << endl;
        Rcout << "df107 " << abs_max << " " << dose_abs_max << " " << Ll_comp[0] << " " << Ll_comp[1] << endl;//prints several convergence terms
    }
    //
    vector<double> dbeta(totalnum,0.0);
    //
    // --------------------------
    // always starts from intial guess
    // --------------------------
    vector<double> beta_p(totalnum,0.0);
    vector<double> beta_c(totalnum,0.0);
    vector<double> beta_a(totalnum,0.0);
    vector<double> beta_best(totalnum,0.0);
    VectorXd::Map(&beta_p[0], beta_0.size()) = beta_0;// stores previous parameters
    VectorXd::Map(&beta_c[0], beta_0.size()) = beta_0;// stores current parameters
    VectorXd::Map(&beta_a[0], beta_0.size()) = beta_0;// stores a refrence value for parameters
    VectorXd::Map(&beta_best[0], beta_0.size()) = beta_0;// stores the best parameters
    double Ll_best = 0.0; //a comparison log-likelihood
    int halves = 0; //number of half-steps taken
    int ind0 = fir; //used for validations
    int i = ind0;
    int iteration=0; //iteration number
    //
    int iter_stop =0;
    //
    while ((iteration < maxiter)&&(iter_stop==0)){
        iteration++;
        beta_p = beta_c;//
        beta_a = beta_c;//
        beta_best = beta_c;//
        //
        // Calcualtes the initial change in parameter
        Calc_Change( nthreads, totalnum, fir, der_iden, dbeta_cap, dose_abs_max, lr, abs_max, Ll, Lld, Lldd, dbeta, change_all, tform, dint, KeepConstant, debugging);
        if (verbose){
            Rcout << "Starting Halves"<<endl;//prints the final changes for validation
        }
        //
        Ll_best = Ll[ind0];
        i = ind0;
        //
        halves=0;
        while ((Ll[ind0] <= Ll_best)&&(halves<halfmax)){ //repeats until half-steps maxed or an improvement
            beta_p = beta_c;//
            beta_a = beta_c;//
            beta_best = beta_c;//
            halves++;
            //Refreshes the matrices used
            Dose = MatrixXd::Zero(df0.rows(),term_tot);
            nonDose = MatrixXd::Constant(df0.rows(),term_tot,0.0);
            nonDose_LIN = MatrixXd::Constant(df0.rows(),term_tot,0.0);
            nonDose_PLIN = MatrixXd::Constant(df0.rows(),term_tot,1.0);
            nonDose_LOGLIN = MatrixXd::Constant(df0.rows(),term_tot,1.0);
            T0 = MatrixXd::Zero(df0.rows(), totalnum); //preallocates matrix for Term column
            Td0 = MatrixXd::Zero(df0.rows(), totalnum); //preallocates matrix for Term derivative columns
            Tdd0 = MatrixXd::Zero(df0.rows(), totalnum*(totalnum+1)/2); //preallocates matrix for Term second derivative columns
            for (int ij=0;ij<totalnum;ij++){
                beta_0[ij] = beta_a[ij] + dbeta[ij];
                beta_c[ij] = beta_0[ij];
            }
            // ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------//
            // The same subterm, risk, sides, and log-likelihood calculations are performed every half-step and iteration
            // ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------//
            Make_Subterms( totalnum, Term_n, tform, dfc, fir, T0, Td0, Tdd0, Dose,  nonDose, TTerm,  nonDose_LIN, nonDose_PLIN, nonDose_LOGLIN ,beta_0, df0,dint,nthreads, debugging);
            if (verbose){
                Rcout << "values checked ";
                for (int ijk=0;ijk<totalnum;ijk++){
                    Rcout << beta_c[ijk] << " ";
                }
                Rcout << " " << endl;
                Rcout << "sums checked ";
                for (int ijk=0;ijk<totalnum;ijk++){
                    Rcout << T0.col(ijk).sum() << " ";
                }
                Rcout << " " << endl;
                Rcout << "derivs checked ";
                for (int ijk=0;ijk<totalnum;ijk++){
                    Rcout << Td0.col(ijk).sum() << " ";
                }
                Rcout << " " << endl;
                Rcout << "second derivs checked ";
                for (int ijk=0;ijk<totalnum;ijk++){
                    Rcout << Tdd0.col(ijk*(ijk+1)/2+ijk).sum() << " ";
                }
                Rcout << " " << endl;
                Rcout << "dose checked ";
                for (int ijk=0;ijk<term_tot;ijk++){
                    Rcout << Dose.col(ijk).array().sum() << " ";
                }
                Rcout << " " << endl;
                Rcout << "non-dose checked ";
                for (int ijk=0;ijk<term_tot;ijk++){
                    Rcout << nonDose.col(ijk).array().sum() << " ";
                }
                Rcout << " " << endl;
                Rcout << "LIN_non-dose checked ";
                for (int ijk=0;ijk<term_tot;ijk++){
                    Rcout << nonDose_LIN.col(ijk).array().sum() << " ";
                }
                Rcout << " " << endl;
                Rcout << "PLIN_non-dose checked ";
                for (int ijk=0;ijk<term_tot;ijk++){
                    Rcout << nonDose_PLIN.col(ijk).array().sum() << " ";
                }
                Rcout << " " << endl;
                Rcout << "LOGLIN_non-dose checked ";
                for (int ijk=0;ijk<term_tot;ijk++){
                    Rcout << nonDose_LOGLIN.col(ijk).array().sum() << " ";
                }
                Rcout << " " << endl;
            }
            //
            RdR = MatrixXd::Zero(df0.rows(), totalnum);
            RddR = MatrixXd::Zero(df0.rows(), totalnum*(totalnum+1)/2);
            //
            //
            Make_Risks(modelform, tform, Term_n, totalnum, fir, T0, Td0, Tdd0, Te, R, Rd, Rdd, Dose, nonDose, TTerm, nonDose_LIN, nonDose_PLIN, nonDose_LOGLIN, RdR, RddR, nthreads, debugging);
            //
            RdR = (RdR.array().isFinite()).select(RdR,0);
            RddR = (RddR.array().isFinite()).select(RddR,0);
            //
            if (verbose){
                Rcout << "risk checked ";
                for (int ijk=0;ijk<totalnum;ijk++){
                    Rcout << R.col(0).sum() << " ";
                }
                Rcout << " " << endl;
                Rcout << "risk1 checked ";
                for (int ijk=0;ijk<totalnum;ijk++){
                    Rcout << Rd.col(ijk).sum() << " ";
                }
                Rcout << " " << endl;
                Rcout << "risk2 checked ";
                for (int ijk=0;ijk<totalnum;ijk++){
                    Rcout << Rdd.col(ijk*(ijk+1)/2+ijk).sum() << " ";
                }
                Rcout << " " << endl;
                end_point = system_clock::now();
                ending = time_point_cast<microseconds>(end_point).time_since_epoch().count();
                Rcout <<"df100 "<<(ending-start)<<" "<<halves<<" "<<iteration<<" "<<ind0<<",Update_R"<<endl;
                //
                gibtime = system_clock::to_time_t(system_clock::now());
                Rcout << ctime(&gibtime) << endl;
            }
            fill(Ll.begin(), Ll.end(), 0.0);
            fill(Lld.begin(), Lld.end(), 0.0);
            fill(Lldd.begin(), Lldd.end(), 0.0);
            Rls1 =MatrixXd::Zero(ntime, 1); //precomputes a series of sums used frequently in the log-liklihood calculations
            Rls2 =MatrixXd::Zero(ntime, totalnum); //Many are repeated due to the same risk groups and derivatives being used at mulitple points
            Rls3 =MatrixXd::Zero(ntime, totalnum*(totalnum+1)/2);
            Lls1 =MatrixXd::Zero(ntime, 1);
            Lls2 =MatrixXd::Zero(ntime, totalnum);
            Lls3 =MatrixXd::Zero(ntime, totalnum*(totalnum+1)/2);
            Calculate_Sides( RiskFail, RiskGroup, totalnum, ntime, R, Rd, Rdd, Rls1, Rls2, Rls3, Lls1, Lls2, Lls3,nthreads, debugging);
            //
            if (verbose){
                Rcout << "riskr checked ";
                for (int ijk=0;ijk<totalnum;ijk++){
                    Rcout << Rls1.col(0).sum() << " ";
                }
                Rcout << " " << endl;
                Rcout << "risk1r checked ";
                for (int ijk=0;ijk<totalnum;ijk++){
                    Rcout << Rls2.col(ijk).sum() << " ";
                }
                Rcout << " " << endl;
                Rcout << "risk2r checked ";
                for (int ijk=0;ijk<totalnum;ijk++){
                    Rcout << Rls3.col(ijk*(ijk+1)/2+ijk).sum() << " ";
                }
                Rcout << " " << endl;
                //
                //
                Rcout << "riskl checked ";
                for (int ijk=0;ijk<totalnum;ijk++){
                    Rcout << Lls1.col(0).sum() << " ";
                }
                Rcout << " " << endl;
                Rcout << "risk1l checked ";
                for (int ijk=0;ijk<totalnum;ijk++){
                    Rcout << Lls2.col(ijk).sum() << " ";
                }
                Rcout << " " << endl;
                Rcout << "risk2l checked ";
                for (int ijk=0;ijk<totalnum;ijk++){
                    Rcout << Lls3.col(ijk*(ijk+1)/2+ijk).sum() << " ";
                }
                Rcout << " " << endl;
            }
            //
            Calc_LogLik( nthreads, RiskFail, RiskGroup, totalnum, ntime, R, Rd, Rdd,RdR,RddR, Rls1, Rls2, Rls3, Lls1, Lls2, Lls3, Ll, Lld, Lldd, debugging, ties_method);
            
            if (change_all){ //If every covariate is to be changed
                if (Ll[ind0] <= Ll_best){//takes a half-step if needed
                    #pragma omp parallel for num_threads(nthreads)
                    for (int ijk=0;ijk<totalnum;ijk++){
                        dbeta[ijk] = dbeta[ijk] / 2.0;
                    }
                } else{//If improved, updates the best vector
                    #pragma omp parallel for num_threads(nthreads)
                    for (int ijk=0;ijk<totalnum;ijk++){
                        beta_best[ijk] = beta_c[ijk];
                    }
                }
            } else {//For validation, the step is always carried over
                //used if a single parameter is being changed
                #pragma omp parallel for num_threads(nthreads)
                for (int ijk=0;ijk<totalnum;ijk++){
                    beta_best[ijk] = beta_c[ijk];
                }
            }
            if (verbose){
                end_point = system_clock::now();
                ending = time_point_cast<microseconds>(end_point).time_since_epoch().count();
                Rcout <<"df100 "<<(ending-start)<<" "<<halves<<" "<<iteration<<" "<<ind0<<",Update_calc"<<endl;
                gibtime = system_clock::to_time_t(system_clock::now());
                Rcout << ctime(&gibtime) << endl;
                Rcout << "df101 ";
                for (int ij=0;ij<totalnum;ij++){
                    Rcout << Ll[ij] << " ";
                }
                Rcout << " " << endl;
                Rcout << "df102 ";
                for (int ij=0;ij<totalnum;ij++){
                    Rcout << Lld[ij] << " ";
                }
                Rcout << " " << endl;
                Rcout << "df103 ";
                for (int ij=0;ij<totalnum;ij++){
                    Rcout << Lldd[ij*totalnum+ij] << " ";
                }
                for (int ij=0;ij<totalnum;ij++){
                    if (abs(Lld[ij]) > Lld_worst){
                        Lld_worst = abs(Lld[ij]);
                    }
                }
                Rcout << " " << endl;
                Rcout << "df104 ";
                for (int ij=0;ij<totalnum;ij++){
                    Rcout << beta_c[ij] << " ";
                }
                Rcout << " " << endl;
                Rcout << "df105 ";
                for (int ij=0;ij<totalnum;ij++){
                    Rcout << Lld[ij]/Lldd[ij*totalnum+ij] << " ";
                }
                Rcout << " " << endl;
                Rcout << "df106 ";
                for (int ij=0;ij<totalnum;ij++){
                    Rcout << Ll[ij]/Lld[ij] << " ";
                }
                Rcout << " " << endl;
                Rcout << "df107 " << abs_max << " " << dose_abs_max << " " << Ll_comp[0] << " " << Ll_comp[1] << endl;
            }
            #pragma omp parallel for num_threads(nthreads)
            for (int ijk=0;ijk<totalnum;ijk++){//totalnum*(totalnum+1)/2
                beta_0[ijk] = beta_c[ijk];
            }
        }
        if (beta_best!=beta_c){//if the risk matrices aren't the optimal values, then they must be recalculated
            if (verbose){
                Rcout << "Changing back to best"<<endl;
            }
            beta_p = beta_c;//
            beta_a = beta_c;//
            beta_best = beta_c;//
            Dose = MatrixXd::Zero(df0.rows(),term_tot);
            nonDose = MatrixXd::Constant(df0.rows(),term_tot,0.0);
            nonDose_LIN = MatrixXd::Constant(df0.rows(),term_tot,0.0);
            nonDose_PLIN = MatrixXd::Constant(df0.rows(),term_tot,1.0);
            nonDose_LOGLIN = MatrixXd::Constant(df0.rows(),term_tot,1.0);
            for (int ij=0;ij<totalnum;ij++){
                beta_0[ij] = beta_a[ij] + dbeta[ij];
                beta_c[ij] = beta_0[ij];
            }
            Make_Subterms( totalnum, Term_n, tform, dfc, fir, T0, Td0, Tdd0, Dose, nonDose, TTerm,  nonDose_LIN, nonDose_PLIN, nonDose_LOGLIN,beta_0, df0,dint,nthreads, debugging);;
            if (verbose){
                Rcout << "values checked ";
                for (int ijk=0;ijk<totalnum;ijk++){
                    Rcout << beta_c[ijk] << " ";
                }
                Rcout << " " << endl;
                //
                //
                Rcout << "sums checked ";
                for (int ijk=0;ijk<totalnum;ijk++){
                    Rcout << T0.col(ijk).sum() << " ";
                }
                Rcout << " " << endl;
                Rcout << "derivs checked ";
                for (int ijk=0;ijk<totalnum;ijk++){
                    Rcout << Td0.col(ijk).sum() << " ";
                }
                Rcout << " " << endl;
                Rcout << "second derivs checked ";
                for (int ijk=0;ijk<totalnum;ijk++){
                    Rcout << Tdd0.col(ijk*(ijk+1)/2+ijk).sum() << " ";
                }
                Rcout << " " << endl;
                Rcout << "dose checked ";
                for (int ijk=0;ijk<1;ijk++){
                    Rcout << Dose.array().sum() << " ";
                }
                Rcout << " " << endl;
            }
            //
            RdR = MatrixXd::Zero(df0.rows(), totalnum);
            RddR = MatrixXd::Zero(df0.rows(), totalnum*(totalnum+1)/2);
            //
            //
            Make_Risks(modelform, tform, Term_n, totalnum, fir, T0, Td0, Tdd0, Te, R, Rd, Rdd, Dose, nonDose, TTerm, nonDose_LIN, nonDose_PLIN, nonDose_LOGLIN, RdR, RddR, nthreads, debugging);
            R = (R.array().isFinite()).select(R,0);
            Rd = (Rd.array().isFinite()).select(Rd,0);
            Rdd = (Rdd.array().isFinite()).select(Rdd,0);
            RdR = (RdR.array().isFinite()).select(RdR,0);
            RddR = (RddR.array().isFinite()).select(RddR,0);
            //
            temp_list = List::create(_["betas"]=wrap(beta_0),_["Status"]="FAILED"); //used as a dummy return value for code checking
            if (R.minCoeff()<=0){
                for (int ijk=0;ijk<totalnum;ijk++){
                    if (T0.col(ijk).minCoeff()<=0){
                        Rcout << ijk << " had a non-positive term" << endl;
                    }
                }
                for (int ijk=0;ijk<term_tot;ijk++){
                    if (nonDose_LIN.col(ijk).minCoeff()<=0){
                        Rcout << ijk << " had a non-positive Lin term" << endl;
                    }
                    if (nonDose_LOGLIN.col(ijk).minCoeff()<=0){
                        Rcout << ijk << " had a non-positive loglin term" << endl;
                    }
                    if (nonDose_PLIN.col(ijk).minCoeff()<=0){
                        Rcout << ijk << " had a non-positive plin term" << endl;
                    }
                    if (Dose.col(ijk).minCoeff()<=0){
                        Rcout << ijk << " had a non-positive dose term" << endl;
                    }
                    if (nonDose.col(ijk).minCoeff()<=0){
                        Rcout << ijk << " had a non-positive nondose term" << endl;
                    }
                    if (TTerm.col(ijk).minCoeff()<=0){
                        Rcout << ijk << " had a non-positive total term" << endl;
                    }
                }
                Rcout << R.sum() << endl;
                Rcout << "A non-positive risk was detected: " << R.minCoeff() << endl;
                return temp_list;
            }
            if (verbose){
                Rcout << "risk checked ";
                for (int ijk=0;ijk<totalnum;ijk++){
                    Rcout << R.col(0).sum() << " ";
                }
                Rcout << " " << endl;
                Rcout << "risk1 checked ";
                for (int ijk=0;ijk<totalnum;ijk++){
                    Rcout << Rd.col(ijk).sum() << " ";
                }
                Rcout << " " << endl;
                Rcout << "risk2 checked ";
                for (int ijk=0;ijk<totalnum;ijk++){
                    Rcout << Rdd.col(ijk*(ijk+1)/2+ijk).sum() << " ";
                }
                Rcout << " " << endl;
                end_point = system_clock::now();
                ending = time_point_cast<microseconds>(end_point).time_since_epoch().count();
                Rcout <<"df100 "<<(ending-start)<<" "<<halves<<" "<<iteration<<" "<<ind0<<",Update_R"<<endl;
                //
                gibtime = system_clock::to_time_t(system_clock::now());
                Rcout << ctime(&gibtime) << endl;
            }
            fill(Ll.begin(), Ll.end(), 0.0);
            fill(Lld.begin(), Lld.end(), 0.0);
            fill(Lldd.begin(), Lldd.end(), 0.0);
            Rls1 =MatrixXd::Zero(ntime, 1); //precomputes a series of sums used frequently in the log-liklihood calculations
            Rls2 =MatrixXd::Zero(ntime, totalnum); //Many are repeated due to the same risk groups and derivatives being used at mulitple points
            Rls3 =MatrixXd::Zero(ntime, totalnum*(totalnum+1)/2);
            Lls1 =MatrixXd::Zero(ntime, 1);
            Lls2 =MatrixXd::Zero(ntime, totalnum);
            Lls3 =MatrixXd::Zero(ntime, totalnum*(totalnum+1)/2);
            Calculate_Sides( RiskFail, RiskGroup, totalnum, ntime, R, Rd, Rdd, Rls1, Rls2, Rls3, Lls1, Lls2, Lls3,nthreads, debugging);
            //
            if (verbose){
                Rcout << "riskr checked ";
                for (int ijk=0;ijk<totalnum;ijk++){
                    Rcout << Rls1.col(0).sum() << " ";
                }
                Rcout << " " << endl;
                Rcout << "risk1r checked ";
                for (int ijk=0;ijk<totalnum;ijk++){
                    Rcout << Rls2.col(ijk).sum() << " ";
                }
                Rcout << " " << endl;
                Rcout << "risk2r checked ";
                for (int ijk=0;ijk<totalnum;ijk++){
                    Rcout << Rls3.col(ijk*(ijk+1)/2+ijk).sum() << " ";
                }
                Rcout << " " << endl;
                //
                //
                Rcout << "riskl checked ";
                for (int ijk=0;ijk<totalnum;ijk++){
                    Rcout << Lls1.col(0).sum() << " ";
                }
                Rcout << " " << endl;
                Rcout << "risk1l checked ";
                for (int ijk=0;ijk<totalnum;ijk++){
                    Rcout << Lls2.col(ijk).sum() << " ";
                }
                Rcout << " " << endl;
                Rcout << "risk2l checked ";
                for (int ijk=0;ijk<totalnum;ijk++){
                    Rcout << Lls3.col(ijk*(ijk+1)/2+ijk).sum() << " ";
                }
                Rcout << " " << endl;
            }
            //
            Calc_LogLik( nthreads, RiskFail, RiskGroup, totalnum, ntime, R, Rd, Rdd,RdR,RddR, Rls1, Rls2, Rls3, Lls1, Lls2, Lls3, Ll, Lld, Lldd, debugging, ties_method);
        }
        for (int ij=0;ij<totalnum;ij++){
            if (abs(Lld[ij]) > Lld_worst){
                Lld_worst = abs(Lld[ij]);
            }
        }
        if (iteration > totalnum){//Doesn't check the first several iterations for convergence
            if (iteration % (2*totalnum)){//Checks every set number of iterations
                if (Lld_worst < deriv_epsilon){//ends if the derivatives are low enough
                    iter_stop = 1;
                }
                Ll_comp[1]=Ll[0];
                if (abs(Ll_comp[1]-Ll_comp[0])/abs(Ll_comp[1])<.01){//if the change in log-likelihood isn't high enough, the maximum step size is reduced
                    abs_max = abs_max*0.1; // reduces the step sizes
                    dose_abs_max = dose_abs_max*0.5;
                }
                if (abs_max < epsilon/10){//if the maximum change is too low, then it ends
                    iter_stop = 1;
                }
            }
        }
        if (verbose){
            end_point = system_clock::now();
            ending = time_point_cast<microseconds>(end_point).time_since_epoch().count();
            Rcout <<"df100 "<<(ending-start)<<" "<<halves<<" "<<iteration<<" "<<ind0<<",Recalc"<<endl;
            gibtime = system_clock::to_time_t(system_clock::now());
            Rcout << ctime(&gibtime) << endl;
            Rcout << "df101 ";
            for (int ij=0;ij<totalnum;ij++){
                Rcout << Ll[ij] << " ";
            }
            Rcout << " " << endl;
            Rcout << "df102 ";
            for (int ij=0;ij<totalnum;ij++){
                Rcout << Lld[ij] << " ";
            }
            Rcout << " " << endl;
            Rcout << "df103 ";
            for (int ij=0;ij<totalnum;ij++){
                Rcout << Lldd[ij*totalnum+ij] << " ";
            }
            Rcout << " " << endl;
            Rcout << "df104 ";
            for (int ij=0;ij<totalnum;ij++){
                Rcout << beta_c[ij] << " ";
            }
            Rcout << " " << endl;
            Rcout << "df105 ";
            for (int ij=0;ij<totalnum;ij++){//prints the newton step value for zero derivative
                Rcout << Lld[ij]/Lldd[ij*totalnum+ij] << " ";
            }
            Rcout << " " << endl;
            Rcout << "df106 ";
            for (int ij=0;ij<totalnum;ij++){//prints the newton step value for zero log-likelihood
                Rcout << Ll[ij]/Lld[ij] << " ";
            }
            Rcout << " " << endl;
            Rcout << "df107 " << abs_max << " " << dose_abs_max << " " << Ll_comp[0] << " " << Ll_comp[1] << endl;//prints several convergence terms
            Rcout << "Finshed iteration" << endl;
        }
    }
    // -----------------------------------------------
    // Performing Full Calculation to get full second derivative matrix
    // -----------------------------------------------
    fill(Ll.begin(), Ll.end(), 0.0);
    fill(Lld.begin(), Lld.end(), 0.0);
    fill(Lldd.begin(), Lldd.end(), 0.0);
    Calculate_Sides( RiskFail, RiskGroup, totalnum, ntime, R, Rd, Rdd, Rls1, Rls2, Rls3, Lls1, Lls2, Lls3,nthreads, debugging);
    //
    if (verbose){
        Rcout << "Wrapping up" << endl;
    }
    //
    Calc_LogLik( nthreads, RiskFail, RiskGroup, totalnum, ntime, R, Rd, Rdd,RdR,RddR, Rls1, Rls2, Rls3, Lls1, Lls2, Lls3, Ll, Lld, Lldd, debugging, ties_method);
    //
    List para_list = List::create(_["Term_n"]=Term_n,_["tforms"]=tform); //stores the term information
    List control_list = List::create(_["Iteration"]=iteration); //stores the total number of iterations used
    NumericVector Lldd_vec = wrap(Lldd);//
    Lldd_vec.attr("dim") = Dimension(totalnum, totalnum);
    //
    const Map<MatrixXd> Lldd_mat(as<Map<MatrixXd> >(Lldd_vec));
    //
    MatrixXd Lldd_inv = -1 * Lldd_mat.inverse().matrix(); //uses inverse information matrix to calculate the standard deviation
    //
    List res_list = List::create(_["LogLik"]=wrap(Ll),_["First_Der"]=wrap(Lld),_["Second_Der"]=Lldd_vec,_["beta_0"]=wrap(beta_0) ,_["Standard_Deviation"]=wrap(Lldd_inv.diagonal().cwiseSqrt()) ,_["AIC"]=2*totalnum-2*Ll[fir],_["Parameter_Lists"]=para_list,_["Control_List"]=control_list);
    // returns a list of results
    return res_list;
}


//' Primary Cox PH baseline hazard function
//' \code{Cox_PH_PLOT_SURV} Performs the calls to calculation functions, Uses calculated risks and risk groups to approximate the baseline, With verbose option prints out time stamps and intermediate sums of terms and derivatives
//'
//' @param     Term_n    Term numbers
//' @param     tform    subterm types
//' @param     a_n    starting values
//' @param     x_all    covariate matrix
//' @param     dfc    covariate column numbers
//' @param     fir    first term number
//' @param     der_iden    subterm number for derivative tests
//' @param     modelform    model string
//' @param     abs_max    Maximum allowed parameter change
//' @param     dose_abs_max    Maximum allowed threshold parameter change
//' @param     df_groups    matrix with time and event information
//' @param     tu    event times
//' @param     verbose    verbosity boolean
//' @param     debugging    debugging boolean
//' @param     KeepConstant    vector identifying constant parameters
//' @param     term_tot    total number of terms
//'
//' @return List of results: baseline harzard, risk for each row
// [[Rcpp::export]]
List Cox_PH_PLOT_SURV(IntegerVector Term_n, StringVector tform, NumericVector a_n,NumericMatrix x_all,IntegerVector dfc,int fir, int der_iden,string modelform, double abs_max,double dose_abs_max, NumericMatrix df_groups, NumericVector tu , bool verbose, bool debugging, IntegerVector KeepConstant, int term_tot){
    //
    // Calculates the baseline hazard
    //
    using namespace std::chrono;
    if (verbose){
        Rcout << "START_NEW" << endl;
    }
    time_point<system_clock> start_point, end_point;
    start_point = system_clock::now();
    auto start = time_point_cast<microseconds>(start_point).time_since_epoch().count();
    end_point = system_clock::now();
    auto ending = time_point_cast<microseconds>(end_point).time_since_epoch().count(); //The time duration is tracked
    //
    auto gibtime = system_clock::to_time_t(system_clock::now());
    if (verbose){
        Rcout << ctime(&gibtime) << endl;
    }
    //
    // Time durations are measured from this point on in microseconds
    //
    // df0: covariate data
    // ntime: number of event times for Cox PH
    // totalnum: number of terms used
    //
    //
    const Map<MatrixXd> df0(as<Map<MatrixXd> >(x_all));
    int ntime = tu.size();
    //
    int totalnum = Term_n.size();
    //
    //
    // cout.precision: controls the number of significant digits printed
    // nthreads: number of threads used for parallel operations
    //
    Rcout.precision(10); //forces higher precision numbers printed to terminal
    int nthreads = Eigen::nbThreads()-1; //stores how many threads are allocated
    //
    //
    //
    // Lld_worst: stores the highest magnitude log-likelihood derivative
    // totem: number of rows needed
    //
    //


    //
    //
    if (verbose){
        end_point = system_clock::now();
        ending = time_point_cast<microseconds>(end_point).time_since_epoch().count();
        Rcout <<"df99,"<<(ending-start)<<",Starting"<<endl;
        gibtime = system_clock::to_time_t(system_clock::now());
        Rcout << ctime(&gibtime) << endl;
    }
    // ---------------------------------------------
    // To Start, needs to seperate the derivative terms
    // ---------------------------------------------
    //
    Map<VectorXd> beta_0(as<Map<VectorXd> >(a_n));
    MatrixXd T0 = MatrixXd::Zero(df0.rows(), totalnum); //preallocates matrix for Term column
    MatrixXd Td0 = MatrixXd::Zero(df0.rows(), totalnum); //preallocates matrix for Term derivative columns
    MatrixXd Tdd0 = MatrixXd::Zero(df0.rows(), totalnum*(totalnum+1)/2); //preallocates matrix for Term second derivative columns
    //
    MatrixXd Te = MatrixXd::Zero(df0.rows(), 1); //preallocates matrix for column terms used for temporary storage
    MatrixXd R = MatrixXd::Zero(df0.rows(), 1); //preallocates matrix for Risks
    MatrixXd Rd = MatrixXd::Zero(df0.rows(), totalnum); //preallocates matrix for Risk derivatives
    MatrixXd Rdd = MatrixXd::Zero(df0.rows(), totalnum*(totalnum+1)/2); //preallocates matrix for Risk second derivatives
    //
    MatrixXd Dose = MatrixXd::Zero(df0.rows(),term_tot); //Matrix of the total dose term values
    MatrixXd nonDose = MatrixXd::Constant(df0.rows(),term_tot,1.0); //Matrix of the total non-dose term values
    MatrixXd nonDose_LIN = MatrixXd::Constant(df0.rows(),term_tot,0.0); //matrix of Linear subterm values
    MatrixXd nonDose_PLIN = MatrixXd::Constant(df0.rows(),term_tot,1.0); //matrix of Loglinear subterm values
    MatrixXd nonDose_LOGLIN = MatrixXd::Constant(df0.rows(),term_tot,1.0); //matrix of Product linear subterm values
    MatrixXd TTerm=MatrixXd::Zero(Dose.rows(),Dose.cols()); //matrix of term values
    double dint = dose_abs_max; //The amount of change used to calculate derivatives in threshold paramters    //
    // totalnum,& Term_n,  tform, dfc,& fir,& T0,& Td0,& Tdd0,& Dose,& nonDose,& beta_0,& df0, dint, nthreads,  debugging
    if (verbose){
        Rcout << "starting subterms " << term_tot << endl;
    }
    Make_Subterms( totalnum, Term_n, tform, dfc, fir, T0, Td0, Tdd0, Dose, nonDose, TTerm,  nonDose_LIN, nonDose_PLIN, nonDose_LOGLIN ,beta_0, df0,dint,nthreads, debugging);
    // ---------------------------------------------------------
    // Prints off a series of calculations to check at what point values are changing
    // ---------------------------------------------------------    
    if (verbose){
        Rcout << "values checked ";
        for (int ijk=0;ijk<totalnum;ijk++){
            Rcout << beta_0[ijk] << " ";
        }
        Rcout << " " << endl;
        Rcout << "sums checked ";
        for (int ijk=0;ijk<totalnum;ijk++){
            Rcout << T0.col(ijk).sum() << " ";
        }
        Rcout << " " << endl;
        Rcout << "derivs checked ";
        for (int ijk=0;ijk<totalnum;ijk++){
            Rcout << Td0.col(ijk).sum() << " ";
        }
        Rcout << " " << endl;
        Rcout << "second derivs checked ";
        for (int ijk=0;ijk<totalnum;ijk++){
            Rcout << Tdd0.col(ijk*(ijk+1)/2+ijk).sum() << " ";
        }
        Rcout << " " << endl;
        Rcout << "dose checked ";
        for (int ijk=0;ijk<term_tot;ijk++){
            Rcout << Dose.col(ijk).array().sum() << " ";
        }
        Rcout << " " << endl;
        Rcout << "non-dose checked ";
        for (int ijk=0;ijk<term_tot;ijk++){
            Rcout << nonDose.col(ijk).array().sum() << " ";
        }
        Rcout << " " << endl;
        Rcout << "LIN_non-dose checked ";
        for (int ijk=0;ijk<term_tot;ijk++){
            Rcout << nonDose_LIN.col(ijk).array().sum() << " ";
        }
        Rcout << " " << endl;
        Rcout << "PLIN_non-dose checked ";
        for (int ijk=0;ijk<term_tot;ijk++){
            Rcout << nonDose_PLIN.col(ijk).array().sum() << " ";
        }
        Rcout << " " << endl;
        Rcout << "LOGLIN_non-dose checked ";
        for (int ijk=0;ijk<term_tot;ijk++){
            Rcout << nonDose_LOGLIN.col(ijk).array().sum() << " ";
        }
        Rcout << " " << endl;
    }
    //
    MatrixXd RdR = MatrixXd::Zero(df0.rows(), totalnum); //preallocates matrix for Risk to derivative ratios
    MatrixXd RddR = MatrixXd::Zero(df0.rows(), totalnum*(totalnum+1)/2); //preallocates matrix for Risk to second derivative ratios
    //
    if (verbose){
        end_point = system_clock::now();
        ending = time_point_cast<microseconds>(end_point).time_since_epoch().count();
        Rcout <<"df99,"<<(ending-start)<<",Prep_Terms"<<endl;
        gibtime = system_clock::to_time_t(system_clock::now());
        Rcout << ctime(&gibtime) << endl;
    }
    //
    // Calculates the risk for each row
    Make_Risks(modelform, tform, Term_n, totalnum, fir, T0, Td0, Tdd0, Te, R, Rd, Rdd, Dose, nonDose, TTerm, nonDose_LIN, nonDose_PLIN, nonDose_LOGLIN, RdR, RddR, nthreads, debugging);
    //
    
    RdR = (RdR.array().isFinite()).select(RdR,0);
    RddR = (RddR.array().isFinite()).select(RddR,0);
    if (verbose){
        end_point = system_clock::now();
        ending = time_point_cast<microseconds>(end_point).time_since_epoch().count();
        Rcout <<"df100 "<<(ending-start)<<" "<<0<<" "<<0<<" "<<-1<<",Prep_R"<<endl;
        gibtime = system_clock::to_time_t(system_clock::now());
        Rcout << ctime(&gibtime) << endl;
    }
    //
    // -------------------------------------------------------------------------------------------
    //
    vector<double> baseline(ntime,0.0);
    //
    // Iterates through the risk groups and approximates the baseline
    //
    const Map<MatrixXd> df_m(as<Map<MatrixXd> >(df_groups));
    #pragma omp parallel for schedule(dynamic) num_threads(nthreads)
    for (int ijk=0;ijk<ntime;ijk++){
        double t0 = tu[ijk];
        VectorXi select_ind_all = ((df_m.col(0).array() <= t0)&&(df_m.col(1).array()>=t0)).cast<int>(); //indices at risk
        vector<int> indices_all;
        VectorXi select_ind_end = ((df_m.col(2).array() == 1)&&(df_m.col(1).array()==t0)).cast<int>(); //indices with events
        vector<int> indices_end;
        //
        //
        int th = 1;
        visit_lambda(select_ind_all,
            [&indices_all, th](double v, int i, int j) {
                if (v==th)
                    indices_all.push_back(i+1);
            });
        visit_lambda(select_ind_end,
            [&indices_end, th](double v, int i, int j) {
                if (v==th)
                    indices_end.push_back(i+1);
            });
        //
        vector<int> indices; //generates vector of (start,end) pairs for indices at risk
        for (auto it = begin (indices_all); it != end (indices_all); ++it) {
            if (indices.size()==0){
                indices.push_back(*it);
                indices.push_back(*it);
            } else if (indices[indices.size()-1]+1<*it){
                indices.push_back(*it);
                indices.push_back(*it);
            } else {
                indices[indices.size()-1] = *it;
            }
        }
        int dj = indices_end[indices_end.size()-1] - indices_end[0] + 1;// number of events
        double Rs1 = 0; //total risk
        for (vector<double>::size_type i = 0; i < indices.size()-1; i=i+2){
            Rs1 += R.block(indices[i]-1,0,indices[i+1]-indices[i]+1,1).sum();
        }
        baseline[ijk] = dj / Rs1; //approximates the baseline hazard
        //
    }
    //
    // returns the baseline approximates and the risk information
    List res_list = List::create(_["baseline"]=wrap(baseline), _["Risks"]=wrap(R));
    //
    return res_list;
}

//' Primary Cox PH risk plotting function
//' \code{Cox_PH_PLOT_RISK} Performs the calls to calculation functions, Uses formula and a generated list of covariate values to calculate risk over a grid, With verbose option prints out time stamps and intermediate sums of terms and derivatives
//'
//' @param     Term_n    Term numbers
//' @param     tform    subterm types
//' @param     a_n    starting values
//' @param     x_all    covariate matrix
//' @param     dfc    covariate column numbers
//' @param     fir    first term number
//' @param     der_iden    subterm number for derivative tests
//' @param     modelform    model string
//' @param     abs_max    Maximum allowed parameter change
//' @param     dose_abs_max    Maximum allowed threshold parameter change
//' @param     df_groups    matrix with time and event information
//' @param     tu    event times
//' @param     verbose    verbosity boolean
//' @param     debugging    debugging boolean
//' @param     KeepConstant    vector identifying constant parameters
//' @param     term_tot    total number of terms
//' @param     uniq_v    number of unqiue covariate values
//'
//' @return List of results: covariate values, risks for each row
// [[Rcpp::export]]
List Cox_PH_PLOT_RISK(IntegerVector Term_n, StringVector tform, NumericVector a_n,NumericMatrix x_all,IntegerVector dfc,int fir, int der_iden,string modelform, double abs_max,double dose_abs_max, NumericMatrix df_groups, NumericVector tu, bool verbose, bool debugging, IntegerVector KeepConstant, int term_tot, int uniq_v){
    ;
    //
    // Plots the risk over a series of covariate values
    //
    using namespace std::chrono;
    if (verbose){
        Rcout << "START_NEW" << endl;
    }
    time_point<system_clock> start_point, end_point;
    start_point = system_clock::now();
    auto start = time_point_cast<microseconds>(start_point).time_since_epoch().count();
    end_point = system_clock::now();
    auto ending = time_point_cast<microseconds>(end_point).time_since_epoch().count(); //The time duration is tracked
    //
    auto gibtime = system_clock::to_time_t(system_clock::now());
    if (verbose){
        Rcout << ctime(&gibtime) << endl;
    }
    //
    // Time durations are measured from this point on in microseconds
    //
    //
    // df1: covariate data
    // totalnum: number of terms used
    const Map<MatrixXd> df1(as<Map<MatrixXd> >(x_all));
    //
    int totalnum = Term_n.size();
    //
    if (verbose){
        Rcout << "Term checked ";
        for (int ij=0;ij<totalnum;ij++){
            Rcout << Term_n[ij] << " ";
        }
        Rcout << " " << endl;
    }
    //
    //
    // cout.precision: controls the number of significant digits printed
    // nthreads: number of threads used for parallel operations
    //
    Rcout.precision(7); //forces higher precision numbers printed to terminal
    int nthreads = Eigen::nbThreads()-1; //stores how many threads are allocated
    //
    //
    //
    if (verbose){
        end_point = system_clock::now();
        ending = time_point_cast<microseconds>(end_point).time_since_epoch().count();
        Rcout <<"df99,"<<(ending-start)<<",Starting"<<endl;
        gibtime = system_clock::to_time_t(system_clock::now());
        Rcout << ctime(&gibtime) << endl;
    }
    // ---------------------------------------------
    // To Start, needs to seperate the derivative terms
    // ---------------------------------------------
    //
    double dint = dose_abs_max; //The amount of change used to calculate derivatives in threshold paramters
    //
    // totalnum,& Term_n,  tform, dfc,& fir,& T0,& Td0,& Tdd0,& Dose,& nonDose,& beta_0,& df0, dint, nthreads,  debugging
    if (verbose){
        Rcout << "starting subterms " << term_tot << endl;
    }
    // ---------------------------------------------
    // To Start, needs to seperate the derivative terms
    // ---------------------------------------------
    //
    float dx = 0;
    if (der_iden >=0){
        ;
    } else {
        throw invalid_argument( "Incorrect parameter to plot by" );
    }
    vector<float> vv; //stores the covariate values
    if (uniq_v > 10){
        vv.resize(100); //continuous covariates use 100 steps
    } else{
        vv.resize(uniq_v); //factor covariates use the number of factors
    }
    MatrixXd df0 = MatrixXd::Zero(vv.size(), df1.cols()); // stores memory for the derivative term parameters and columns
    //
    Map<VectorXd> beta_0(as<Map<VectorXd> >(a_n));
    MatrixXd T0 = MatrixXd::Zero(df0.rows(), totalnum); //preallocates matrix for Term column
    MatrixXd Td0 = MatrixXd::Zero(df0.rows(), totalnum); //preallocates matrix for Term derivative columns
    MatrixXd Tdd0 = MatrixXd::Zero(df0.rows(), totalnum*(totalnum+1)/2); //preallocates matrix for Term second derivative columns
    //
    MatrixXd Te = MatrixXd::Zero(df0.rows(), 1); //preallocates matrix for column terms used for temporary storage
    MatrixXd R = MatrixXd::Zero(df0.rows(), 1); //preallocates matrix for Risks
    MatrixXd Rd = MatrixXd::Zero(df0.rows(), totalnum); //preallocates matrix for Risk derivatives
    MatrixXd Rdd = MatrixXd::Zero(df0.rows(), totalnum*(totalnum+1)/2); //preallocates matrix for Risk second derivatives
    //
    MatrixXd Dose = MatrixXd::Zero(df0.rows(),term_tot); //Matrix of the total dose term values
    MatrixXd nonDose = MatrixXd::Constant(df0.rows(),term_tot,1.0); //Matrix of the total non-dose term values
    MatrixXd nonDose_LIN = MatrixXd::Constant(df0.rows(),term_tot,0.0); //matrix of Linear subterm values
    MatrixXd nonDose_PLIN = MatrixXd::Constant(df0.rows(),term_tot,1.0); //matrix of Loglinear subterm values
    MatrixXd nonDose_LOGLIN = MatrixXd::Constant(df0.rows(),term_tot,1.0); //matrix of Product linear subterm values
    MatrixXd TTerm=MatrixXd::Zero(Dose.rows(),Dose.cols()); //matrix of term values
    //
    int ijk=der_iden;
    dx = (df1.col(ijk).maxCoeff() - df1.col(ijk).minCoeff())/(vv.size()-1);//varies from max to minimum
    vv[0] = df1.col(ijk).minCoeff();
    generate(vv.begin(), vv.end(), [n = 0, &dx]() mutable { return n++ * dx; });
    #pragma omp parallel for schedule(dynamic) num_threads(nthreads)
    for (vector<float>::size_type ij=0;ij<vv.size();ij++){
        df0(ij,der_iden)=vv[ij]; //fills the column with varying values
    }
    //
    // Calculates terms
    Make_Subterms( totalnum, Term_n, tform, dfc, fir, T0, Td0, Tdd0, Dose, nonDose, TTerm,  nonDose_LIN, nonDose_PLIN, nonDose_LOGLIN ,beta_0, df0,dint,nthreads, debugging);
    MatrixXd RdR = MatrixXd::Zero(df0.rows(), totalnum); //preallocates matrix for Risk to derivative ratios
    MatrixXd RddR = MatrixXd::Zero(df0.rows(), totalnum*(totalnum+1)/2);
    //
    // Calculates risk
    Make_Risks(modelform, tform, Term_n, totalnum, fir, T0, Td0, Tdd0, Te, R, Rd, Rdd, Dose, nonDose, TTerm, nonDose_LIN, nonDose_PLIN, nonDose_LOGLIN, RdR, RddR, nthreads, debugging);
    //
    // Returns the values used and risks calculated
    List res_list = List::create(_["x"]=wrap(df1.col(der_iden)), _["y"]=wrap(R.col(0)));//returns list of covariate values and risk
    return res_list;
}

//' Primary Cox PH schoenfeld residual function
//' \code{Schoenfeld_Cox_PH} Performs the calls to calculation functions, Uses calculated risks and risk groups to calculate the residuals, With verbose option prints out time stamps and intermediate sums of terms and derivatives
//'
//' @param     Term_n    Term numbers
//' @param     tform    subterm types
//' @param     a_n    starting values
//' @param     x_all    covariate matrix
//' @param     dfc    covariate column numbers
//' @param     fir    first term number
//' @param     der_iden    subterm number for derivative tests
//' @param     modelform    model string
//' @param     abs_max    Maximum allowed parameter change
//' @param     dose_abs_max    Maximum allowed threshold parameter change
//' @param     df_groups    matrix with time and event information
//' @param     tu    event times
//' @param     verbose    verbosity boolean
//' @param     debugging    debugging boolean
//' @param     KeepConstant    vector identifying constant parameters
//' @param     term_tot    total number of terms
//' @param     ties_method    ties method
//'
//' @return List of results: scaled schoenfeld residuals
// [[Rcpp::export]]
NumericMatrix Schoenfeld_Cox_PH( IntegerVector Term_n, StringVector tform, NumericVector a_n,NumericMatrix x_all,IntegerVector dfc,int fir, int der_iden,string modelform, double abs_max,double dose_abs_max, NumericMatrix df_groups, NumericVector tu , bool verbose, bool debugging, IntegerVector KeepConstant, int term_tot, string ties_method){
    ;
    //
    // Calculates the schoenfeld residuals
    //
    List temp_list = List::create(_["Status"]="FAILED"); //used as a dummy return value for code checking
    if (verbose){
        Rcout << "START_NEW" << endl;
    }
    //
    // Time durations are measured from this point on in microseconds
    //
    time_point<system_clock> start_point, end_point;
    start_point = system_clock::now();
    auto start = time_point_cast<microseconds>(start_point).time_since_epoch().count();
    end_point = system_clock::now();
    auto ending = time_point_cast<microseconds>(end_point).time_since_epoch().count(); //The time duration is tracked
    //
    auto gibtime = system_clock::to_time_t(system_clock::now());
    if (verbose){
        Rcout << ctime(&gibtime) << endl;
    }
    //
    //
    // df0: covariate data
    // ntime: number of event times for Cox PH
    // totalnum: number of terms used
    //
    const Map<MatrixXd> df0(as<Map<MatrixXd> >(x_all));
    int ntime = tu.size();
    //
    int totalnum = Term_n.size();
    //
    if (verbose){
        Rcout << "Term checked ";
        for (int ij=0;ij<totalnum;ij++){
            Rcout << Term_n[ij] << " ";
        }
        Rcout << " " << endl;
    }
    //
    //
    // cout.precision: controls the number of significant digits printed
    // nthreads: number of threads used for parallel operations
    //
    Rcout.precision(7); //forces higher precision numbers printed to terminal
    int nthreads = Eigen::nbThreads()-1; //stores how many threads are allocated
    //
    //
    //
    // Lld_worst: stores the highest magnitude log-likelihood derivative
    // totem: number of rows needed
    //
    //
    if (verbose){
        end_point = system_clock::now();
        ending = time_point_cast<microseconds>(end_point).time_since_epoch().count();
        Rcout <<"df99,"<<(ending-start)<<",Starting"<<endl;
        gibtime = system_clock::to_time_t(system_clock::now());
        Rcout << ctime(&gibtime) << endl;
    }
    // ---------------------------------------------
    // To Start, needs to seperate the derivative terms
    // ---------------------------------------------
    //
    Map<VectorXd> beta_0(as<Map<VectorXd> >(a_n));
    MatrixXd T0 = MatrixXd::Zero(df0.rows(), totalnum); //preallocates matrix for Term column
    MatrixXd Td0 = MatrixXd::Zero(df0.rows(), totalnum); //preallocates matrix for Term derivative columns
    MatrixXd Tdd0 = MatrixXd::Zero(df0.rows(), totalnum*(totalnum+1)/2); //preallocates matrix for Term second derivative columns
    //
    MatrixXd Te = MatrixXd::Zero(df0.rows(), 1); //preallocates matrix for column terms used for temporary storage
    MatrixXd R = MatrixXd::Zero(df0.rows(), 1); //preallocates matrix for Risks
    MatrixXd Rd = MatrixXd::Zero(df0.rows(), totalnum); //preallocates matrix for Risk derivatives
    MatrixXd Rdd = MatrixXd::Zero(df0.rows(), totalnum*(totalnum+1)/2); //preallocates matrix for Risk second derivatives
    //
    MatrixXd Dose = MatrixXd::Zero(df0.rows(),term_tot); //Matrix of the total dose term values
    MatrixXd nonDose = MatrixXd::Constant(df0.rows(),term_tot,1.0); //Matrix of the total non-dose term values
    MatrixXd nonDose_LIN = MatrixXd::Constant(df0.rows(),term_tot,0.0); //matrix of Linear subterm values
    MatrixXd nonDose_PLIN = MatrixXd::Constant(df0.rows(),term_tot,1.0); //matrix of Loglinear subterm values
    MatrixXd nonDose_LOGLIN = MatrixXd::Constant(df0.rows(),term_tot,1.0); //matrix of Product linear subterm values
    MatrixXd TTerm=MatrixXd::Zero(Dose.rows(),Dose.cols()); //matrix of term values
    double dint = dose_abs_max; //The amount of change used to calculate derivatives in threshold paramters    //
    // Calculates the subterm and term values
    if (verbose){
        Rcout << "starting subterms " << term_tot << endl;
    }
    Make_Subterms( totalnum, Term_n, tform, dfc, fir, T0, Td0, Tdd0, Dose, nonDose, TTerm,  nonDose_LIN, nonDose_PLIN, nonDose_LOGLIN ,beta_0, df0,dint,nthreads, debugging);
    // ---------------------------------------------------------
    // Prints off a series of calculations to check at what point values are changing
    // ---------------------------------------------------------
    //
    MatrixXd RdR = MatrixXd::Zero(df0.rows(), totalnum); //preallocates matrix for Risk to derivative ratios
    MatrixXd RddR = MatrixXd::Zero(df0.rows(), totalnum*(totalnum+1)/2);
    //
    Make_Risks(modelform, tform, Term_n, totalnum, fir, T0, Td0, Tdd0, Te, R, Rd, Rdd, Dose, nonDose, TTerm, nonDose_LIN, nonDose_PLIN, nonDose_LOGLIN, RdR, RddR, nthreads, debugging);
    //
    vector<string>  RiskGroup(ntime); //vector of strings detailing the rows
    IntegerMatrix RiskFail(ntime,2); //vector giving the event rows
    //
    const Map<MatrixXd> df_m(as<Map<MatrixXd> >(df_groups));
    // Creates matrices used to identify the event risk groups
    Make_Groups( ntime, df_m, RiskFail, RiskGroup, tu, nthreads, debugging);
    //
    //
    // --------------------------
    // now a vector exists with row locations
    // --------------------------
    MatrixXd Rls1 =MatrixXd::Zero(ntime, 1); //precomputes a series of sums used frequently in the log-liklihood calculations
    MatrixXd Rls2 =MatrixXd::Zero(ntime, totalnum); //Many are repeated due to the same risk groups and derivatives being used at mulitple points
    MatrixXd Rls3 =MatrixXd::Zero(ntime, totalnum*(totalnum+1)/2); //Sum and its derivatives are precomputed
    MatrixXd Lls1 =MatrixXd::Zero(ntime, 1); //The log-likelihood calculation has a Right and Left sum used
    MatrixXd Lls2 =MatrixXd::Zero(ntime, totalnum);
    MatrixXd Lls3 =MatrixXd::Zero(ntime, totalnum*(totalnum+1)/2);
    //The log-likelihood is calculated in parallel over the risk groups
    vector<double> Ll(totalnum,0.0); //Log-likelihood values
    vector<double> Lld(totalnum,0.0); //Log-likelihood derivative values
    vector<double> Lldd(pow(totalnum,2),0.0);//The second derivative matrix has room for every combination, but only the lower triangle is calculated initially
    //
    // Calculates the side sum terms used
    Calculate_Sides( RiskFail, RiskGroup, totalnum, ntime, R, Rd, Rdd, Rls1, Rls2, Rls3, Lls1, Lls2, Lls3,nthreads, debugging);
    //
    // Calculates log-likelihood
    Calc_LogLik( nthreads, RiskFail, RiskGroup, totalnum, ntime, R, Rd, Rdd,RdR,RddR, Rls1, Rls2, Rls3, Lls1, Lls2, Lls3, Ll, Lld, Lldd, debugging, ties_method);
    //
    NumericVector Lldd_vec = wrap(Lldd);//
    Lldd_vec.attr("dim") = Dimension(totalnum, totalnum);
    //
    const Map<MatrixXd> Lldd_mat(as<Map<MatrixXd> >(Lldd_vec));
    //
    MatrixXd Lldd_inv = -1 * Lldd_mat.inverse().matrix(); //uses inverse information matrix to calculate the standard deviation
    //
    NumericVector stdev = wrap(Lldd_inv.diagonal().cwiseSqrt()); //vector of standard deviations
    // --------------------------
    // now a vector exists with row locations
    // --------------------------
    MatrixXd residuals = MatrixXd::Zero(ntime,totalnum);
    #pragma omp parallel for schedule(dynamic) num_threads(nthreads) collapse(2)
    for (int ijk=0;ijk<totalnum;ijk++){//totalnum*(totalnum+1)/2
        for (int j=0;j<ntime;j++){
            //
            vector<int> InGroup;
            string Groupstr = RiskGroup[j];
            stringstream ss(Groupstr);
            //
            for (int i; ss >> i;) {
                InGroup.push_back(i);    
                if (ss.peek() == ',')
                    ss.ignore();
            }
            double t_sum =0;
            double x_expect =0;
            double Vscale=0;
            //
            // calculates the total term value
            //
            for (vector<double>::size_type i = 0; i < InGroup.size()-1; i=i+2){
                t_sum += R.block(InGroup[i]-1,0,InGroup[i+1]-InGroup[i]+1,1).sum();
                x_expect +=  (df0.block(InGroup[i]-1,ijk,InGroup[i+1]-InGroup[i]+1,1).array() * R.block(InGroup[i]-1,0,InGroup[i+1]-InGroup[i]+1,1).array()).sum();
            }
            x_expect = x_expect / t_sum; //calculates the averaged covariate value
            int dj = RiskFail(j,1)-RiskFail(j,0)+1;
            double x_risks = df0.block(RiskFail(j,0),ijk,dj,1).sum()/dj; //calculate the average covariate value with events
            //
            Vscale = dj/stdev(ijk); //the residual is scaled by the number of events divided by a measure of the variation
            //
            residuals(j,ijk) = (x_risks - x_expect)*Vscale;
        }
    }
    // returns residuals
    return wrap(residuals);
}

//' Primary poisson regression function
//' \code{LogLik_Poisson} Performs the calls to calculation functions, Structures the poisson regression, With verbose option prints out time stamps and intermediate sums of terms and derivatives
//'
//' @param     PyrC    person-year matrix
//' @param     Term_n    Term numbers
//' @param     tform    subterm types
//' @param     a_n    starting values
//' @param     x_all    covariate matrix
//' @param     dfc    covariate column numbers
//' @param     fir    first term number
//' @param     der_iden    subterm number for derivative tests
//' @param     modelform    model string
//' @param     lr    learning rate for newton step toward 0 derivative
//' @param     maxiter    maximum number of iterations
//' @param     halfmax    maximum number of half steps
//' @param     epsilon    minimum acceptable maximum parameter change
//' @param     dbeta_cap    learning rate for newton step toward 0 log-likelihood
//' @param     abs_max    Maximum allowed parameter change
//' @param     dose_abs_max    Maximum allowed threshold parameter change
//' @param     deriv_epsilon    threshold for near-zero derivative
//' @param     change_all    boolean if every parameter is being updated
//' @param     verbose    verbosity boolean
//' @param     debugging    debugging boolean
//' @param     KeepConstant    vector identifying constant parameters
//' @param     term_tot    total number of terms
//'
//' @return List of results: Log-likelihood of optimum, first derivative of log-likelihood, second derivative matrix, parameter list, standard deviation estimate, AIC, deviance, model information
// [[Rcpp::export]]
List LogLik_Poisson( MatrixXd PyrC, IntegerVector Term_n, StringVector tform, NumericVector a_n,NumericMatrix x_all,IntegerVector dfc,int fir, int der_iden,string modelform, double lr, int maxiter, int halfmax, double epsilon, double dbeta_cap, double abs_max,double dose_abs_max, double deriv_epsilon,bool change_all, bool verbose, bool debugging, IntegerVector KeepConstant, int term_tot){
    ;
    //
    List temp_list = List::create(_["Status"]="FAILED"); //used as a dummy return value for code checking
    if (verbose){
        Rcout << "START_NEW" << endl;
    }
    time_point<system_clock> start_point, end_point;
    start_point = system_clock::now();
    auto start = time_point_cast<microseconds>(start_point).time_since_epoch().count();
    end_point = system_clock::now();
    auto ending = time_point_cast<microseconds>(end_point).time_since_epoch().count(); //The time duration is tracked
    //
    //
    // Time durations are measured from this point on in microseconds
    //
    auto gibtime = system_clock::to_time_t(system_clock::now());
    if (verbose){
        Rcout << ctime(&gibtime) << endl;
    }
    //
    //
    // df0: covariate data
    // totalnum: number of terms used
    //
    const Map<MatrixXd> df0(as<Map<MatrixXd> >(x_all));
    //
    int totalnum = Term_n.size();
    //
    if (verbose){
        Rcout << "Term checked ";
        for (int ij=0;ij<totalnum;ij++){
            Rcout << Term_n[ij] << " ";
        }
        Rcout << " " << endl;
    }
    //
    //
    // cout.precision: controls the number of significant digits printed
    // nthreads: number of threads used for parallel operations
    //
    Rcout.precision(7); //forces higher precision numbers printed to terminal
    int nthreads = Eigen::nbThreads()-1; //stores how many threads are allocated
    //
    //
    // Lld_worst: stores the highest magnitude log-likelihood derivative
    // totem: number of rows needed
    //
    double Lld_worst = 0.0; //stores derivative value used to determine if every parameter is near convergence    //
    if (verbose){
        end_point = system_clock::now();
        ending = time_point_cast<microseconds>(end_point).time_since_epoch().count();
        Rcout <<"df99,"<<(ending-start)<<",Starting"<<endl;
        gibtime = system_clock::to_time_t(system_clock::now());
        Rcout << ctime(&gibtime) << endl;
    }
    // ---------------------------------------------
    // To Start, needs to seperate the derivative terms
    // ---------------------------------------------
    //
    Map<VectorXd> beta_0(as<Map<VectorXd> >(a_n));
    MatrixXd T0 = MatrixXd::Zero(df0.rows(), totalnum); //preallocates matrix for Term column
    MatrixXd Td0 = MatrixXd::Zero(df0.rows(), totalnum); //preallocates matrix for Term derivative columns
    MatrixXd Tdd0 = MatrixXd::Zero(df0.rows(), totalnum*(totalnum+1)/2); //preallocates matrix for Term second derivative columns
    //
    MatrixXd Te = MatrixXd::Zero(df0.rows(), 1); //preallocates matrix for column terms used for temporary storage
    MatrixXd R = MatrixXd::Zero(df0.rows(), 1); //preallocates matrix for Risks
    MatrixXd Rd = MatrixXd::Zero(df0.rows(), totalnum); //preallocates matrix for Risk derivatives
    MatrixXd Rdd = MatrixXd::Zero(df0.rows(), totalnum*(totalnum+1)/2); //preallocates matrix for Risk second derivatives
    //
    MatrixXd Dose = MatrixXd::Zero(df0.rows(),term_tot); //Matrix of the total dose term values
    MatrixXd nonDose = MatrixXd::Constant(df0.rows(),term_tot,1.0); //Matrix of the total non-dose term values
    MatrixXd nonDose_LIN = MatrixXd::Constant(df0.rows(),term_tot,0.0); //matrix of Linear subterm values
    MatrixXd nonDose_PLIN = MatrixXd::Constant(df0.rows(),term_tot,1.0); //matrix of Loglinear subterm values
    MatrixXd nonDose_LOGLIN = MatrixXd::Constant(df0.rows(),term_tot,1.0); //matrix of Product linear subterm values
    MatrixXd TTerm=MatrixXd::Zero(Dose.rows(),Dose.cols()); //matrix of term values
    double dint = dose_abs_max; //The amount of change used to calculate derivatives in threshold paramters    //
    if (verbose){
        Rcout << "starting subterms " << term_tot << endl;
    }
    // Calculates the subterm and term values
    Make_Subterms( totalnum, Term_n, tform, dfc, fir, T0, Td0, Tdd0, Dose, nonDose, TTerm,  nonDose_LIN, nonDose_PLIN, nonDose_LOGLIN ,beta_0, df0,dint,nthreads, debugging);
    // ---------------------------------------------------------
    // Prints off a series of calculations to check at what point values are changing
    // ---------------------------------------------------------
    if (verbose){
        Rcout << "values checked ";
        for (int ijk=0;ijk<totalnum;ijk++){
            Rcout << beta_0[ijk] << " ";
        }
        Rcout << " " << endl;
        Rcout << "sums checked ";
        for (int ijk=0;ijk<totalnum;ijk++){
            Rcout << T0.col(ijk).sum() << " ";
        }
        Rcout << " " << endl;
        Rcout << "derivs checked ";
        for (int ijk=0;ijk<totalnum;ijk++){
            Rcout << Td0.col(ijk).sum() << " ";
        }
        Rcout << " " << endl;
        Rcout << "second derivs checked ";
        for (int ijk=0;ijk<totalnum;ijk++){
            Rcout << Tdd0.col(ijk*(ijk+1)/2+ijk).sum() << " ";
        }
        Rcout << " " << endl;
        Rcout << "dose checked ";
        for (int ijk=0;ijk<term_tot;ijk++){
            Rcout << Dose.col(ijk).array().sum() << " ";
        }
        Rcout << " " << endl;
        Rcout << "non-dose checked ";
        for (int ijk=0;ijk<term_tot;ijk++){
            Rcout << nonDose.col(ijk).array().sum() << " ";
        }
        Rcout << " " << endl;
        Rcout << "LIN_non-dose checked ";
        for (int ijk=0;ijk<term_tot;ijk++){
            Rcout << nonDose_LIN.col(ijk).array().sum() << " ";
        }
        Rcout << " " << endl;
        Rcout << "PLIN_non-dose checked ";
        for (int ijk=0;ijk<term_tot;ijk++){
            Rcout << nonDose_PLIN.col(ijk).array().sum() << " ";
        }
        Rcout << " " << endl;
        Rcout << "LOGLIN_non-dose checked ";
        for (int ijk=0;ijk<term_tot;ijk++){
            Rcout << nonDose_LOGLIN.col(ijk).array().sum() << " ";
        }
        Rcout << " " << endl;
    }
    //
    MatrixXd RdR = MatrixXd::Zero(df0.rows(), totalnum); //preallocates matrix for Risk to derivative ratios
    MatrixXd RddR = MatrixXd::Zero(df0.rows(), totalnum*(totalnum+1)/2);
    //
    if (verbose){
        end_point = system_clock::now();
        ending = time_point_cast<microseconds>(end_point).time_since_epoch().count();
        Rcout <<"df99,"<<(ending-start)<<",Prep_Terms"<<endl;
        gibtime = system_clock::to_time_t(system_clock::now());
        Rcout << ctime(&gibtime) << endl;
    }
    //
    // Calculates the risk for each row
    Make_Risks(modelform, tform, Term_n, totalnum, fir, T0, Td0, Tdd0, Te, R, Rd, Rdd, Dose, nonDose, TTerm, nonDose_LIN, nonDose_PLIN, nonDose_LOGLIN, RdR, RddR, nthreads, debugging);
    //
    
    RdR = (RdR.array().isFinite()).select(RdR,0);
    RddR = (RddR.array().isFinite()).select(RddR,0);
    //
    if (R.minCoeff()<=0){
        for (int ijk=0;ijk<totalnum;ijk++){
            if (T0.col(ijk).minCoeff()<=0){
                Rcout << ijk << " had a non-positive term" << endl;
            }
        }
        for (int ijk=0;ijk<term_tot;ijk++){
            if (nonDose_LIN.col(ijk).minCoeff()<=0){
                Rcout << ijk << " had a non-positive Lin term" << endl;
            }
            if (nonDose_LOGLIN.col(ijk).minCoeff()<=0){
                Rcout << ijk << " had a non-positive loglin term" << endl;
            }
            if (nonDose_PLIN.col(ijk).minCoeff()<=0){
                Rcout << ijk << " had a non-positive plin term" << endl;
            }
            if (Dose.col(ijk).minCoeff()<=0){
                Rcout << ijk << " had a non-positive dose term" << endl;
            }
            if (nonDose.col(ijk).minCoeff()<=0){
                Rcout << ijk << " had a non-positive nondose term" << endl;
            }
            if (TTerm.col(ijk).minCoeff()<=0){
                Rcout << ijk << " had a non-positive total term" << endl;
            }
        }
        Rcout << R.sum() << endl;
        Rcout << "A non-positive risk was detected: " << R.minCoeff() << endl;
        return temp_list;
    }
    //
    if (verbose){
        Rcout << "risk checked ";
        for (int ijk=0;ijk<totalnum;ijk++){
            Rcout << R.col(0).sum() << " ";
        }
        Rcout << " " << endl;
        Rcout << "risk1 checked ";
        for (int ijk=0;ijk<totalnum;ijk++){
            Rcout << Rd.col(ijk).sum() << " ";
        }
        Rcout << " " << endl;
        Rcout << "risk2 checked ";
        for (int ijk=0;ijk<totalnum;ijk++){
            Rcout << Rdd.col(ijk*(ijk+1)/2+ijk).sum() << " ";
        }
        Rcout << " " << endl;
        //
        end_point = system_clock::now();
        ending = time_point_cast<microseconds>(end_point).time_since_epoch().count();
        Rcout <<"df100 "<<(ending-start)<<" "<<0<<" "<<0<<" "<<-1<<",Prep_R"<<endl;
        gibtime = system_clock::to_time_t(system_clock::now());
        Rcout << ctime(&gibtime) << endl;
    }
    //
    // -------------------------------------------------------------------------------------------
    //
    vector<double> Ll(totalnum,0.0);
    vector<double> Lld(totalnum,0.0);
    vector<double> Lldd(pow(totalnum,2),0.0);
    //
    //
    Poisson_LogLik( nthreads, totalnum, PyrC, R, Rd, Rdd, RdR, RddR, Ll, Lld, Lldd, debugging);
    MatrixXd dev_temp = MatrixXd::Zero(PyrC.rows(),2);
    dev_temp.col(0) = PyrC.col(0).array() * R.col(0).array();
    dev_temp.col(0) = PyrC.col(1).array() * dev_temp.col(0).array().pow(-1).array();
    dev_temp = (dev_temp.array().isFinite()).select(dev_temp,0);
    dev_temp.col(0) = dev_temp.col(0).array().log().array();
    dev_temp = (dev_temp.array().isFinite()).select(dev_temp,0);
    dev_temp.col(0) = PyrC.col(1).array() * dev_temp.col(0).array();
    dev_temp.col(1) = PyrC.col(1).array() - PyrC.col(0).array() * R.col(0).array();
    dev_temp = (dev_temp.array().isFinite()).select(dev_temp,0);
    double dev = 2*(dev_temp.col(0).sum() - dev_temp.col(1).sum()); //deviation calculation is split into steps
    //
    vector <double> Ll_comp(2,Ll[0]); //vector to compare values
    //
    if (verbose){
        end_point = system_clock::now();
        ending = time_point_cast<microseconds>(end_point).time_since_epoch().count();
        Rcout <<"df100 "<<(ending-start)<<" "<<0<<" "<<0<<" "<<0<<",Calc"<<endl;//prints the time
        gibtime = system_clock::to_time_t(system_clock::now());
        Rcout << ctime(&gibtime) << endl;
        Rcout << "df101 ";//prints the log-likelihoods
        for (int ij=0;ij<totalnum;ij++){
            Rcout << Ll[ij] << " ";
        }
        Rcout << " " << endl;
        Rcout << "df102 ";//prints the first derivatives
        for (int ij=0;ij<totalnum;ij++){
            Rcout << Lld[ij] << " ";
        }
        Rcout << " " << endl;
        Rcout << "df103 ";//prints the second derivatives
        for (int ij=0;ij<totalnum;ij++){
            Rcout << Lldd[ij*totalnum+ij] << " ";
        }
        for (int ij=0;ij<totalnum;ij++){//locates highest magnitude derivative
            if (abs(Lld[ij]) > Lld_worst){
                Lld_worst = abs(Lld[ij]);
            }
        }
        Rcout << " " << endl;
        Rcout << "df104 ";//prints parameter values
        for (int ij=0;ij<totalnum;ij++){
            Rcout << beta_0[ij] << " ";
        }
        Rcout << " " << endl;
        Rcout << "Checking Deviance " << dev << endl;
        Rcout << "df105 ";
        for (int ij=0;ij<totalnum;ij++){//prints the newton step value for zero derivative
            Rcout << Lld[ij]/Lldd[ij*totalnum+ij] << " ";
        }
        Rcout << " " << endl;
        Rcout << "df106 ";
        for (int ij=0;ij<totalnum;ij++){//prints the newton step value for zero log-likelihood
            Rcout << Ll[ij]/Lld[ij] << " ";
        }
        Rcout << " " << endl;
        Rcout << "df107 " << abs_max << " " << dose_abs_max << " " << Ll_comp[0] << " " << Ll_comp[1] << endl;//prints several convergence terms
    }
    //
    vector<double> dbeta(totalnum,0.0);
    //
    // --------------------------
    // always starts from intial guess
    // --------------------------
    vector<double> beta_p(totalnum,0.0);
    vector<double> beta_c(totalnum,0.0);
    vector<double> beta_a(totalnum,0.0);
    vector<double> beta_best(totalnum,0.0);
    VectorXd::Map(&beta_p[0], beta_0.size()) = beta_0;// stores previous parameters
    VectorXd::Map(&beta_c[0], beta_0.size()) = beta_0;// stores current parameters
    VectorXd::Map(&beta_a[0], beta_0.size()) = beta_0;// stores a refrence value for parameters
    VectorXd::Map(&beta_best[0], beta_0.size()) = beta_0;// stores the best parameters
    double Ll_best = 0.0; //a comparison log-likelihood
    int halves = 0; //number of half-steps taken
    int ind0 = fir; //used for validations
    int iteration=0; //iteration number
    //
    //
    while (iteration < maxiter){
        iteration++;
        beta_p = beta_c;//
        beta_a = beta_c;//
        beta_best = beta_c;//
        //
        // Calcualtes the initial change in parameter
        Calc_Change( nthreads, totalnum, fir, der_iden, dbeta_cap, dose_abs_max, lr, abs_max, Ll, Lld, Lldd, dbeta, change_all, tform, dint, KeepConstant, debugging);
        if (verbose){
            Rcout << "Starting Halves"<<endl;//prints the final changes for validation
        }
        //
        Ll_best = Ll[ind0];
        //
        halves=0;
        while ((Ll[ind0] <= Ll_best)&&(halves<halfmax)){ //repeats until half-steps maxed or an improvement
            beta_p = beta_c;//
            beta_a = beta_c;//
            beta_best = beta_c;//
            halves++;
            Dose = MatrixXd::Zero(df0.rows(),term_tot);
            nonDose = MatrixXd::Constant(df0.rows(),term_tot,0.0);
            nonDose_LIN = MatrixXd::Constant(df0.rows(),term_tot,0.0);
            nonDose_PLIN = MatrixXd::Constant(df0.rows(),term_tot,1.0);
            nonDose_LOGLIN = MatrixXd::Constant(df0.rows(),term_tot,1.0);
            T0 = MatrixXd::Zero(df0.rows(), totalnum); //preallocates matrix for Term column
            Td0 = MatrixXd::Zero(df0.rows(), totalnum); //preallocates matrix for Term derivative columns
            Tdd0 = MatrixXd::Zero(df0.rows(), totalnum*(totalnum+1)/2); //preallocates matrix for Term second derivative columns
            for (int ij=0;ij<totalnum;ij++){
                beta_0[ij] = beta_a[ij] + dbeta[ij];
                beta_c[ij] = beta_0[ij];
            }
            // ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------//
            // The same subterm, risk, sides, and log-likelihood calculations are performed every half-step and iteration
            // ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------//
            Make_Subterms( totalnum, Term_n, tform, dfc, fir, T0, Td0, Tdd0, Dose,  nonDose, TTerm,  nonDose_LIN, nonDose_PLIN, nonDose_LOGLIN ,beta_0, df0,dint,nthreads, debugging);
            if (verbose){
                Rcout << "values checked ";
                for (int ijk=0;ijk<totalnum;ijk++){
                    Rcout << beta_c[ijk] << " ";
                }
                Rcout << " " << endl;
                Rcout << "sums checked ";
                for (int ijk=0;ijk<totalnum;ijk++){
                    Rcout << T0.col(ijk).sum() << " ";
                }
                Rcout << " " << endl;
                Rcout << "derivs checked ";
                for (int ijk=0;ijk<totalnum;ijk++){
                    Rcout << Td0.col(ijk).sum() << " ";
                }
                Rcout << " " << endl;
                Rcout << "second derivs checked ";
                for (int ijk=0;ijk<totalnum;ijk++){
                    Rcout << Tdd0.col(ijk*(ijk+1)/2+ijk).sum() << " ";
                }
                Rcout << " " << endl;
                Rcout << "dose checked ";
                for (int ijk=0;ijk<term_tot;ijk++){
                    Rcout << Dose.col(ijk).array().sum() << " ";
                }
                Rcout << " " << endl;
                Rcout << "non-dose checked ";
                for (int ijk=0;ijk<term_tot;ijk++){
                    Rcout << nonDose.col(ijk).array().sum() << " ";
                }
                Rcout << " " << endl;
                Rcout << "LIN_non-dose checked ";
                for (int ijk=0;ijk<term_tot;ijk++){
                    Rcout << nonDose_LIN.col(ijk).array().sum() << " ";
                }
                Rcout << " " << endl;
                Rcout << "PLIN_non-dose checked ";
                for (int ijk=0;ijk<term_tot;ijk++){
                    Rcout << nonDose_PLIN.col(ijk).array().sum() << " ";
                }
                Rcout << " " << endl;
                Rcout << "LOGLIN_non-dose checked ";
                for (int ijk=0;ijk<term_tot;ijk++){
                    Rcout << nonDose_LOGLIN.col(ijk).array().sum() << " ";
                }
                Rcout << " " << endl;
            }
            //
            RdR = MatrixXd::Zero(df0.rows(), totalnum);
            RddR = MatrixXd::Zero(df0.rows(), totalnum*(totalnum+1)/2);
            //
            //
            Make_Risks(modelform, tform, Term_n, totalnum, fir, T0, Td0, Tdd0, Te, R, Rd, Rdd, Dose, nonDose, TTerm, nonDose_LIN, nonDose_PLIN, nonDose_LOGLIN, RdR, RddR, nthreads, debugging);
            //
            RdR = (RdR.array().isFinite()).select(RdR,0);
            RddR = (RddR.array().isFinite()).select(RddR,0);
            //
            if (verbose){
                Rcout << "risk checked ";
                for (int ijk=0;ijk<totalnum;ijk++){
                    Rcout << R.col(0).sum() << " ";
                }
                Rcout << " " << endl;
                Rcout << "risk1 checked ";
                for (int ijk=0;ijk<totalnum;ijk++){
                    Rcout << Rd.col(ijk).sum() << " ";
                }
                Rcout << " " << endl;
                Rcout << "risk2 checked ";
                for (int ijk=0;ijk<totalnum;ijk++){
                    Rcout << Rdd.col(ijk*(ijk+1)/2+ijk).sum() << " ";
                }
                Rcout << " " << endl;
                end_point = system_clock::now();
                ending = time_point_cast<microseconds>(end_point).time_since_epoch().count();
                Rcout <<"df100 "<<(ending-start)<<" "<<halves<<" "<<iteration<<" "<<ind0<<",Update_R"<<endl;
                //
                gibtime = system_clock::to_time_t(system_clock::now());
                Rcout << ctime(&gibtime) << endl;
            }
            fill(Ll.begin(), Ll.end(), 0.0);
            fill(Lld.begin(), Lld.end(), 0.0);
            fill(Lldd.begin(), Lldd.end(), 0.0);
            Poisson_LogLik( nthreads, totalnum, PyrC, R, Rd, Rdd, RdR, RddR, Ll, Lld, Lldd, debugging);
            dev_temp.col(0) = PyrC.col(0).array() * R.col(0).array();
            dev_temp.col(0) = PyrC.col(1).array() * dev_temp.col(0).array().pow(-1).array();
            dev_temp = (dev_temp.array().isFinite()).select(dev_temp,0);
            dev_temp.col(0) = dev_temp.col(0).array().log().array();
            dev_temp = (dev_temp.array().isFinite()).select(dev_temp,0);
            dev_temp.col(0) = PyrC.col(1).array() * dev_temp.col(0).array();
            dev_temp.col(1) = PyrC.col(1).array() - PyrC.col(0).array() * R.col(0).array();
            dev_temp = (dev_temp.array().isFinite()).select(dev_temp,0);
            dev = 2*(dev_temp.col(0).sum() - dev_temp.col(1).sum());
            
            if (change_all){
                if (Ll[ind0] <= Ll_best){//takes a half-step if needed
                    #pragma omp parallel for num_threads(nthreads)
                    for (int ijk=0;ijk<totalnum;ijk++){
                        dbeta[ijk] = dbeta[ijk] / 2.0;
                    }
                } else{//If improved, updates the best vector
                    #pragma omp parallel for num_threads(nthreads)
                    for (int ijk=0;ijk<totalnum;ijk++){
                        beta_best[ijk] = beta_c[ijk];
                    }
                }
            } else {//For validation, the step is always carried over
                #pragma omp parallel for num_threads(nthreads)
                for (int ijk=0;ijk<totalnum;ijk++){
                    beta_best[ijk] = beta_c[ijk];
                }
            }
            if (verbose){
                end_point = system_clock::now();
                ending = time_point_cast<microseconds>(end_point).time_since_epoch().count();
                Rcout <<"df100 "<<(ending-start)<<" "<<halves<<" "<<iteration<<" "<<ind0<<",Update_calc"<<endl;
                gibtime = system_clock::to_time_t(system_clock::now());
                Rcout << ctime(&gibtime) << endl;
                Rcout << "df101 ";
                for (int ij=0;ij<totalnum;ij++){
                    Rcout << Ll[ij] << " ";
                }
                Rcout << " " << endl;
                Rcout << "df102 ";
                for (int ij=0;ij<totalnum;ij++){
                    Rcout << Lld[ij] << " ";
                }
                Rcout << " " << endl;
                Rcout << "df103 ";
                for (int ij=0;ij<totalnum;ij++){
                    Rcout << Lldd[ij*totalnum+ij] << " ";
                }
                for (int ij=0;ij<totalnum;ij++){
                    if (abs(Lld[ij]) > Lld_worst){
                        Lld_worst = abs(Lld[ij]);
                    }
                }
                Rcout << " " << endl;
                Rcout << "df104 ";
                for (int ij=0;ij<totalnum;ij++){
                    Rcout << beta_c[ij] << " ";
                }
                Rcout << " " << endl;
                Rcout << "Checking Deviance " << dev << endl;
                Rcout << "df105 ";
                for (int ij=0;ij<totalnum;ij++){
                    Rcout << Lld[ij]/Lldd[ij*totalnum+ij] << " ";
                }
                Rcout << " " << endl;
                Rcout << "df106 ";
                for (int ij=0;ij<totalnum;ij++){
                    Rcout << Ll[ij]/Lld[ij] << " ";
                }
                Rcout << " " << endl;
                Rcout << "df107 " << abs_max << " " << dose_abs_max << " " << Ll_comp[0] << " " << Ll_comp[1] << endl;
            }
            #pragma omp parallel for num_threads(nthreads)
            for (int ijk=0;ijk<totalnum;ijk++){//totalnum*(totalnum+1)/2
                beta_0[ijk] = beta_c[ijk];
            }
            beta_best = beta_c;
        }
        if (beta_best!=beta_c){//if the risk matrices aren't the optimal values, then they must be recalculated
            if (verbose){
                Rcout << "Changing back to best"<<endl;
            }
            beta_p = beta_c;//
            beta_a = beta_c;//
            beta_best = beta_c;//
            Dose = MatrixXd::Zero(df0.rows(),term_tot);
            nonDose = MatrixXd::Constant(df0.rows(),term_tot,0.0);
            nonDose_LIN = MatrixXd::Constant(df0.rows(),term_tot,0.0);
            nonDose_PLIN = MatrixXd::Constant(df0.rows(),term_tot,1.0);
            nonDose_LOGLIN = MatrixXd::Constant(df0.rows(),term_tot,1.0);
            for (int ij=0;ij<totalnum;ij++){
                beta_0[ij] = beta_a[ij] + dbeta[ij];
                beta_c[ij] = beta_0[ij];
            }
            Make_Subterms( totalnum, Term_n, tform, dfc, fir, T0, Td0, Tdd0, Dose, nonDose, TTerm,  nonDose_LIN, nonDose_PLIN, nonDose_LOGLIN,beta_0, df0,dint,nthreads, debugging);;
            if (verbose){
                Rcout << "values checked ";
                for (int ijk=0;ijk<totalnum;ijk++){
                    Rcout << beta_c[ijk] << " ";
                }
                Rcout << " " << endl;
                //
                //
                Rcout << "sums checked ";
                for (int ijk=0;ijk<totalnum;ijk++){
                    Rcout << T0.col(ijk).sum() << " ";
                }
                Rcout << " " << endl;
                Rcout << "derivs checked ";
                for (int ijk=0;ijk<totalnum;ijk++){
                    Rcout << Td0.col(ijk).sum() << " ";
                }
                Rcout << " " << endl;
                Rcout << "second derivs checked ";
                for (int ijk=0;ijk<totalnum;ijk++){
                    Rcout << Tdd0.col(ijk*(ijk+1)/2+ijk).sum() << " ";
                }
                Rcout << " " << endl;
                Rcout << "dose checked ";
                for (int ijk=0;ijk<1;ijk++){
                    Rcout << Dose.array().sum() << " ";
                }
                Rcout << " " << endl;
            }
            //
            RdR = MatrixXd::Zero(df0.rows(), totalnum);
            RddR = MatrixXd::Zero(df0.rows(), totalnum*(totalnum+1)/2);
            //
            //
            Make_Risks(modelform, tform, Term_n, totalnum, fir, T0, Td0, Tdd0, Te, R, Rd, Rdd, Dose, nonDose, TTerm, nonDose_LIN, nonDose_PLIN, nonDose_LOGLIN, RdR, RddR, nthreads, debugging);
            R = (R.array().isFinite()).select(R,0);
            Rd = (Rd.array().isFinite()).select(Rd,0);
            Rdd = (Rdd.array().isFinite()).select(Rdd,0);
            RdR = (RdR.array().isFinite()).select(RdR,0);
            RddR = (RddR.array().isFinite()).select(RddR,0);
            //
            temp_list = List::create(_["betas"]=wrap(beta_0),_["Status"]="FAILED"); //used as a dummy return value for code checking
            if (R.minCoeff()<=0){
                for (int ijk=0;ijk<totalnum;ijk++){
                    if (T0.col(ijk).minCoeff()<=0){
                        Rcout << ijk << " had a non-positive term" << endl;
                    }
                }
                for (int ijk=0;ijk<term_tot;ijk++){
                    if (nonDose_LIN.col(ijk).minCoeff()<=0){
                        Rcout << ijk << " had a non-positive Lin term" << endl;
                    }
                    if (nonDose_LOGLIN.col(ijk).minCoeff()<=0){
                        Rcout << ijk << " had a non-positive loglin term" << endl;
                    }
                    if (nonDose_PLIN.col(ijk).minCoeff()<=0){
                        Rcout << ijk << " had a non-positive plin term" << endl;
                    }
                    if (Dose.col(ijk).minCoeff()<=0){
                        Rcout << ijk << " had a non-positive dose term" << endl;
                    }
                    if (nonDose.col(ijk).minCoeff()<=0){
                        Rcout << ijk << " had a non-positive nondose term" << endl;
                    }
                    if (TTerm.col(ijk).minCoeff()<=0){
                        Rcout << ijk << " had a non-positive total term" << endl;
                    }
                }
                Rcout << R.sum() << endl;
                Rcout << "A non-positive risk was detected: " << R.minCoeff() << endl;
                return temp_list;
            }
            if (verbose){
                Rcout << "risk checked ";
                for (int ijk=0;ijk<totalnum;ijk++){
                    Rcout << R.col(0).sum() << " ";
                }
                Rcout << " " << endl;
                Rcout << "risk1 checked ";
                for (int ijk=0;ijk<totalnum;ijk++){
                    Rcout << Rd.col(ijk).sum() << " ";
                }
                Rcout << " " << endl;
                Rcout << "risk2 checked ";
                for (int ijk=0;ijk<totalnum;ijk++){
                    Rcout << Rdd.col(ijk*(ijk+1)/2+ijk).sum() << " ";
                }
                Rcout << " " << endl;
                end_point = system_clock::now();
                ending = time_point_cast<microseconds>(end_point).time_since_epoch().count();
                Rcout <<"df100 "<<(ending-start)<<" "<<halves<<" "<<iteration<<" "<<ind0<<",Update_R"<<endl;
                //
                gibtime = system_clock::to_time_t(system_clock::now());
                Rcout << ctime(&gibtime) << endl;
            }
            fill(Ll.begin(), Ll.end(), 0.0);
            fill(Lld.begin(), Lld.end(), 0.0);
            fill(Lldd.begin(), Lldd.end(), 0.0);
            Poisson_LogLik( nthreads, totalnum, PyrC, R, Rd, Rdd, RdR, RddR, Ll, Lld, Lldd, debugging);
            dev_temp.col(0) = PyrC.col(0).array() * R.col(0).array();
            dev_temp.col(0) = PyrC.col(1).array() * dev_temp.col(0).array().pow(-1).array();
            dev_temp = (dev_temp.array().isFinite()).select(dev_temp,0);
            dev_temp.col(0) = dev_temp.col(0).array().log().array();
            dev_temp = (dev_temp.array().isFinite()).select(dev_temp,0);
            dev_temp.col(0) = PyrC.col(1).array() * dev_temp.col(0).array();
            dev_temp.col(1) = PyrC.col(1).array() - PyrC.col(0).array() * R.col(0).array();
            dev_temp = (dev_temp.array().isFinite()).select(dev_temp,0);
            dev = 2*(dev_temp.col(0).sum() - dev_temp.col(1).sum());
            
        }
        for (int ij=0;ij<totalnum;ij++){            if (abs(Lld[ij]) > Lld_worst){
                Lld_worst = abs(Lld[ij]);
            }
        }
        if (iteration > totalnum){//Sets the minimum number of iterations
            if (iteration % (2*totalnum)){//Checks every set number of iterations
                if (Lld_worst < deriv_epsilon){//ends if the derivatives are low enough
                    iteration = maxiter;
                }
                Ll_comp[1]=Ll[0];
                if (abs(Ll_comp[1]-Ll_comp[0])/abs(Ll_comp[1])<.01){//if the change in log-likelihood isn't high enough, the maximum step size if reduced
                    abs_max = abs_max*0.1;
                    dose_abs_max = dose_abs_max*0.5;
                }
                if (abs_max < epsilon/10){//if the maximum change is too low, then it ends
                    iteration = maxiter;
                }
            }
        }
        if (verbose){
            end_point = system_clock::now();
            ending = time_point_cast<microseconds>(end_point).time_since_epoch().count();
            Rcout <<"df100 "<<(ending-start)<<" "<<halves<<" "<<iteration<<" "<<ind0<<",Recalc"<<endl;
            gibtime = system_clock::to_time_t(system_clock::now());
            Rcout << ctime(&gibtime) << endl;
            Rcout << "df101 ";
            for (int ij=0;ij<totalnum;ij++){
                Rcout << Ll[ij] << " ";
            }
            Rcout << " " << endl;
            Rcout << "df102 ";
            for (int ij=0;ij<totalnum;ij++){
                Rcout << Lld[ij] << " ";
            }
            Rcout << " " << endl;
            Rcout << "df103 ";
            for (int ij=0;ij<totalnum;ij++){
                Rcout << Lldd[ij*totalnum+ij] << " ";
            }
            Rcout << " " << endl;
            Rcout << "df104 ";
            for (int ij=0;ij<totalnum;ij++){
                Rcout << beta_c[ij] << " ";
            }
            Rcout << " " << endl;
            Rcout << "Checking Deviance " << dev << endl;
            Rcout << "Finshed iteration" << endl;
            Rcout << "df105 ";
            for (int ij=0;ij<totalnum;ij++){
                Rcout << Lld[ij]/Lldd[ij*totalnum+ij] << " ";
            }
            Rcout << " " << endl;
            Rcout << "df106 ";
            for (int ij=0;ij<totalnum;ij++){
                Rcout << Ll[ij]/Lld[ij] << " ";
            }
            Rcout << " " << endl;
            Rcout << "df107 " << abs_max << " " << dose_abs_max << " " << Ll_comp[0] << " " << Ll_comp[1] << endl;
        }
    }
    // -----------------------------------------------
    // Performing Full Calculation to get full second derivative matrix
    // -----------------------------------------------
    fill(Ll.begin(), Ll.end(), 0.0);
    fill(Lld.begin(), Lld.end(), 0.0);
    fill(Lldd.begin(), Lldd.end(), 0.0);
    Poisson_LogLik( nthreads, totalnum, PyrC, R, Rd, Rdd, RdR, RddR, Ll, Lld, Lldd, debugging);
    dev_temp.col(0) = PyrC.col(0).array() * R.col(0).array();
    dev_temp.col(0) = PyrC.col(1).array() * dev_temp.col(0).array().pow(-1).array();
    dev_temp = (dev_temp.array().isFinite()).select(dev_temp,0);
    dev_temp.col(0) = dev_temp.col(0).array().log().array();
    dev_temp = (dev_temp.array().isFinite()).select(dev_temp,0);
    dev_temp.col(0) = PyrC.col(1).array() * dev_temp.col(0).array();
    dev_temp.col(1) = PyrC.col(1).array() - PyrC.col(0).array() * R.col(0).array();
    dev_temp = (dev_temp.array().isFinite()).select(dev_temp,0);
    dev = 2*(dev_temp.col(0).sum() - dev_temp.col(1).sum());
    if (verbose){
        end_point = system_clock::now();
        ending = time_point_cast<microseconds>(end_point).time_since_epoch().count();
        Rcout <<"df100 "<<(ending-start)<<" "<<halves<<" "<<iteration<<" "<<ind0<<",Recalc"<<endl;
        gibtime = system_clock::to_time_t(system_clock::now());
        Rcout << ctime(&gibtime) << endl;
        Rcout << "df101 ";
        for (int ij=0;ij<totalnum;ij++){
            Rcout << Ll[ij] << " ";
        }
        Rcout << " " << endl;
        Rcout << "df102 ";
        for (int ij=0;ij<totalnum;ij++){
            Rcout << Lld[ij] << " ";
        }
        Rcout << " " << endl;
        Rcout << "df103 ";
        for (int ij=0;ij<totalnum;ij++){
            Rcout << Lldd[ij*totalnum+ij] << " ";
        }
        Rcout << " " << endl;
        Rcout << "df104 ";
        for (int ij=0;ij<totalnum;ij++){
            Rcout << beta_c[ij] << " ";
        }
        Rcout << " " << endl;
        Rcout << "df105 ";
        for (int ij=0;ij<totalnum;ij++){//prints the newton step value for zero derivative
            Rcout << Lld[ij]/Lldd[ij*totalnum+ij] << " ";
        }
        Rcout << " " << endl;
        Rcout << "df106 ";
        for (int ij=0;ij<totalnum;ij++){//prints the newton step value for zero log-likelihood
            Rcout << Ll[ij]/Lld[ij] << " ";
        }
        Rcout << " " << endl;
        Rcout << "df107 " << abs_max << " " << dose_abs_max << " " << Ll_comp[0] << " " << Ll_comp[1] << endl;//prints several convergence terms
        Rcout << "Checking Deviance " << dev << endl;
        Rcout << "Finshed iteration" << endl;
    }
    // Changes the parameter back into the original form
    List para_list = List::create(_["Term_n"]=Term_n,_["tforms"]=tform);
    List control_list = List::create(_["Iteration"]=iteration);
    NumericVector Lldd_vec = wrap(Lldd);//creates list of dose parameters
    Lldd_vec.attr("dim") = Dimension(totalnum, totalnum);
    //
    const Map<MatrixXd> Lldd_mat(as<Map<MatrixXd> >(Lldd_vec));
    //
    MatrixXd Lldd_inv = -1 * Lldd_mat.inverse().matrix(); //uses inverse information matrix to calculate the standard deviation
    //
    List res_list = List::create(_["LogLik"]=wrap(Ll),_["First_Der"]=wrap(Lld),_["Second_Der"]=Lldd_vec,_["beta_0"]=wrap(beta_0) ,_["Standard_Deviation"]=wrap(Lldd_inv.diagonal().cwiseSqrt()) ,_["AIC"]=2*totalnum-2*Ll[fir],_["Deviation"]=dev,_["Parameter_Lists"]=para_list,_["Control_List"]=control_list);
    // returns a list of results
    return res_list;
}

//' Primary Cox PH stress test function
//' \code{Stress_Run} Performs the calls to calculation functions, Structures running one iteration of the cox PH regression, With verbose option prints out time stamps and intermediate sums of terms and derivatives, Prints out further timestamps within calculation calls
//'
//' @param     Term_n    Term numbers
//' @param     tform    subterm types
//' @param     a_n    starting values
//' @param     x_all    covariate matrix
//' @param     dfc    covariate column numbers
//' @param     fir    first term number
//' @param     der_iden    subterm number for derivative tests
//' @param     modelform    model string
//' @param     lr    learning rate for newton step toward 0 derivative
//' @param     maxiter    maximum number of iterations
//' @param     halfmax    maximum number of half steps
//' @param     epsilon    minimum acceptable maximum parameter change
//' @param     dbeta_cap    learning rate for newton step toward 0 log-likelihood
//' @param     abs_max    Maximum allowed parameter change
//' @param     dose_abs_max    Maximum allowed threshold parameter change
//' @param     deriv_epsilon    threshold for near-zero derivative
//' @param     df_groups    matrix with time and event information
//' @param     tu    event times
//' @param     change_all    boolean if every parameter is being updated
//' @param     verbose    verbosity boolean
//' @param     debugging    debugging boolean
//' @param     KeepConstant    vector identifying constant parameters
//' @param     term_tot    total number of terms
//' @param     debug_checks    string vector of functions to test
//' @param     ties_method    ties method
//'
//' @return NULL
// [[Rcpp::export]]
void Stress_Run( IntegerVector Term_n, StringVector tform, NumericVector a_n,NumericMatrix x_all,IntegerVector dfc,int fir, int der_iden,string modelform, double lr, int maxiter, int halfmax, double epsilon, double dbeta_cap, double abs_max,double dose_abs_max, double deriv_epsilon, NumericMatrix df_groups, NumericVector tu ,bool change_all, bool verbose, bool debugging, IntegerVector KeepConstant, int term_tot, StringVector debug_checks, string ties_method){
;
    //
    // Runs through a single calculation with some functions printing additional information printed
    //
    if (verbose){
        Rcout << "START_NEW" << endl;
    }
    //
    // The labels for the calculations which can have additional information printed
    //
    StringVector aval_str = StringVector::create("MakeMatrix","MakeRisk","MakeGroup","CalcSide","CalcLL","CalcChange","UpdateRisk","IterMakeRisk","IterCalcSide","IterCalcLL");
    LogicalVector Debug_It = in(aval_str,debug_checks);
    if (verbose){
        Rcout << Debug_It << endl; //Booleans for which functions to check further
    }
    //
    // Time durations are measured from this point on in microseconds
    //
    time_point<system_clock> start_point, end_point;
    start_point = system_clock::now();
    auto start = time_point_cast<microseconds>(start_point).time_since_epoch().count();
    end_point = system_clock::now();
    auto ending = time_point_cast<microseconds>(end_point).time_since_epoch().count(); //The time duration is tracked
    //
    auto gibtime = system_clock::to_time_t(system_clock::now());
    if (verbose){
        Rcout << ctime(&gibtime) << endl;
    }
    //
    // df0: covariate data
    // ntime: number of event times for Cox PH
    // totalnum: number of terms used
    //
    const Map<MatrixXd> df0(as<Map<MatrixXd> >(x_all));
    int ntime = tu.size();
    //
    int totalnum = Term_n.size();
    //
    if (verbose){
        Rcout << "Term checked ";
        for (int ij=0;ij<totalnum;ij++){
            Rcout << Term_n[ij] << " ";
        }
        Rcout << " " << endl;
    }
    //
    // cout.precision: controls the number of significant digits printed
    // nthreads: number of threads used for parallel operations
    //
    Rcout.precision(10); //forces higher precision numbers printed to terminal
    int nthreads = Eigen::nbThreads()-1; //stores how many threads are allocated    //
    // Lld_worst: stores the highest magnitude log-likelihood derivative
    // totem: number of rows needed
    //
    double Lld_worst = 0.0; //stores derivative value used to determine if every parameter is near convergence    //
    //
    if (verbose){
        end_point = system_clock::now();
        ending = time_point_cast<microseconds>(end_point).time_since_epoch().count();
        Rcout <<"df99,"<<(ending-start)<<",Starting"<<endl;
        gibtime = system_clock::to_time_t(system_clock::now());
        Rcout << ctime(&gibtime) << endl;
    }
    // ---------------------------------------------
    // To Start, needs to seperate the derivative terms
    // ---------------------------------------------
    //
    Map<VectorXd> beta_0(as<Map<VectorXd> >(a_n));
    MatrixXd T0 = MatrixXd::Zero(df0.rows(), totalnum); //preallocates matrix for Term column
    MatrixXd Td0 = MatrixXd::Zero(df0.rows(), totalnum); //preallocates matrix for Term derivative columns
    MatrixXd Tdd0 = MatrixXd::Zero(df0.rows(), totalnum*(totalnum+1)/2); //preallocates matrix for Term second derivative columns
    //
    MatrixXd Te = MatrixXd::Zero(df0.rows(), 1); //preallocates matrix for column terms used for temporary storage
    MatrixXd R = MatrixXd::Zero(df0.rows(), 1); //preallocates matrix for Risks
    MatrixXd Rd = MatrixXd::Zero(df0.rows(), totalnum); //preallocates matrix for Risk derivatives
    MatrixXd Rdd = MatrixXd::Zero(df0.rows(), totalnum*(totalnum+1)/2); //preallocates matrix for Risk second derivatives
    //
    MatrixXd Dose = MatrixXd::Zero(df0.rows(),term_tot); //Matrix of the total dose term values
    MatrixXd nonDose = MatrixXd::Constant(df0.rows(),term_tot,1.0); //Matrix of the total non-dose term values
    MatrixXd nonDose_LIN = MatrixXd::Constant(df0.rows(),term_tot,0.0); //matrix of Linear subterm values
    MatrixXd nonDose_PLIN = MatrixXd::Constant(df0.rows(),term_tot,1.0); //matrix of Loglinear subterm values
    MatrixXd nonDose_LOGLIN = MatrixXd::Constant(df0.rows(),term_tot,1.0); //matrix of Product linear subterm values
    MatrixXd TTerm=MatrixXd::Zero(Dose.rows(),Dose.cols()); //matrix of term values
    double dint = dose_abs_max; //The amount of change used to calculate derivatives in threshold paramters    //
    //
    // Calculates the subterm and term values
    if (verbose){
        Rcout << "starting subterms " << term_tot << endl;
    }
    if (Debug_It[0]){
            Make_Subterms( totalnum, Term_n, tform, dfc, fir, T0, Td0, Tdd0, Dose, nonDose, TTerm,  nonDose_LIN, nonDose_PLIN, nonDose_LOGLIN ,beta_0, df0,dint,nthreads, TRUE);
    } else {
            Make_Subterms( totalnum, Term_n, tform, dfc, fir, T0, Td0, Tdd0, Dose, nonDose, TTerm,  nonDose_LIN, nonDose_PLIN, nonDose_LOGLIN ,beta_0, df0,dint,nthreads, FALSE);
    }
    // ---------------------------------------------------------
    // Prints off a series of calculations to check at what point values are changing
    // ---------------------------------------------------------    
    if (verbose){
        Rcout << "values checked ";
        for (int ijk=0;ijk<totalnum;ijk++){
            Rcout << beta_0[ijk] << " ";
        }
        Rcout << " " << endl;
        Rcout << "sums checked ";
        for (int ijk=0;ijk<totalnum;ijk++){
            Rcout << T0.col(ijk).sum() << " ";
        }
        Rcout << " " << endl;
        Rcout << "derivs checked ";
        for (int ijk=0;ijk<totalnum;ijk++){
            Rcout << Td0.col(ijk).sum() << " ";
        }
        Rcout << " " << endl;
        Rcout << "second derivs checked ";
        for (int ijk=0;ijk<totalnum;ijk++){
            Rcout << Tdd0.col(ijk*(ijk+1)/2+ijk).sum() << " ";
        }
        Rcout << " " << endl;
        Rcout << "dose checked ";
        for (int ijk=0;ijk<term_tot;ijk++){
            Rcout << Dose.col(ijk).array().sum() << " ";
        }
        Rcout << " " << endl;
        Rcout << "non-dose checked ";
        for (int ijk=0;ijk<term_tot;ijk++){
            Rcout << nonDose.col(ijk).array().sum() << " ";
        }
        Rcout << " " << endl;
        Rcout << "LIN_non-dose checked ";
        for (int ijk=0;ijk<term_tot;ijk++){
            Rcout << nonDose_LIN.col(ijk).array().sum() << " ";
        }
        Rcout << " " << endl;
        Rcout << "PLIN_non-dose checked ";
        for (int ijk=0;ijk<term_tot;ijk++){
            Rcout << nonDose_PLIN.col(ijk).array().sum() << " ";
        }
        Rcout << " " << endl;
        Rcout << "LOGLIN_non-dose checked ";
        for (int ijk=0;ijk<term_tot;ijk++){
            Rcout << nonDose_LOGLIN.col(ijk).array().sum() << " ";
        }
        Rcout << " " << endl;
    }
    //
    List temp_list = List::create(_["betas"]=wrap(beta_0),_["Status"]="FAILED"); //used as a dummy return value for code checking
    MatrixXd RdR = MatrixXd::Zero(df0.rows(), totalnum); //preallocates matrix for Risk to derivative ratios
    MatrixXd RddR = MatrixXd::Zero(df0.rows(), totalnum*(totalnum+1)/2); //preallocates matrix for Risk to second derivative ratios
    //
    if (verbose){
        end_point = system_clock::now();
        ending = time_point_cast<microseconds>(end_point).time_since_epoch().count();
        Rcout <<"df99,"<<(ending-start)<<",Prep_Terms"<<endl;
        gibtime = system_clock::to_time_t(system_clock::now());
        Rcout << ctime(&gibtime) << endl;
    }
    //
    //
    // Calculates the risk for each row
    if (Debug_It[1]){
        Make_Risks(modelform, tform, Term_n, totalnum, fir, T0, Td0, Tdd0, Te, R, Rd, Rdd, Dose, nonDose, TTerm, nonDose_LIN, nonDose_PLIN, nonDose_LOGLIN, RdR, RddR, nthreads, TRUE);
    } else {
        Make_Risks(modelform, tform, Term_n, totalnum, fir, T0, Td0, Tdd0, Te, R, Rd, Rdd, Dose, nonDose, TTerm, nonDose_LIN, nonDose_PLIN, nonDose_LOGLIN, RdR, RddR, nthreads, FALSE);
    }
    //
    //
    // Removes infinite values
    RdR = (RdR.array().isFinite()).select(RdR,0);
    RddR = (RddR.array().isFinite()).select(RddR,0);
    //
    temp_list = List::create(_["betas"]=wrap(beta_0),_["Status"]="FAILED"); //used as a dummy return value for code checking
    //
    if (verbose){
        Rcout << "risk checked ";
        for (int ijk=0;ijk<totalnum;ijk++){
            Rcout << R.col(0).sum() << " ";
        }
        Rcout << " " << endl;
    //    return temp_list;
        Rcout << "risk1 checked ";
        for (int ijk=0;ijk<totalnum;ijk++){
            Rcout << Rd.col(ijk).sum() << " ";
        }
        Rcout << " " << endl;
        Rcout << "risk2 checked ";
        for (int ijk=0;ijk<totalnum;ijk++){
            Rcout << Rdd.col(ijk*(ijk+1)/2+ijk).sum() << " ";
        }
        Rcout << " " << endl;
        //
        end_point = system_clock::now();
        ending = time_point_cast<microseconds>(end_point).time_since_epoch().count();
        Rcout <<"df100 "<<(ending-start)<<" "<<0<<" "<<0<<" "<<-1<<",Prep_R"<<endl;
        gibtime = system_clock::to_time_t(system_clock::now());
        Rcout << ctime(&gibtime) << endl;
    }
    //
    // -------------------------------------------------------------------------------------------
    //
    vector<string>  RiskGroup(ntime); //vector of strings detailing the rows
    IntegerMatrix RiskFail(ntime,2); //vector giving the event rows
    //
    const Map<MatrixXd> df_m(as<Map<MatrixXd> >(df_groups));
    //
    //
    // Creates matrices used to identify the event risk groups
    if (Debug_It[2]){
        Make_Groups( ntime, df_m, RiskFail, RiskGroup, tu, nthreads, TRUE);
    } else {
        Make_Groups( ntime, df_m, RiskFail, RiskGroup, tu, nthreads, FALSE);
    }
    //
    if (verbose){
        end_point = system_clock::now();
        ending = time_point_cast<microseconds>(end_point).time_since_epoch().count();
        Rcout <<"df100 "<<(ending-start)<<" "<<0<<" "<<0<<" "<<-1<<",Prep_List"<<endl;
        gibtime = system_clock::to_time_t(system_clock::now());
        Rcout << ctime(&gibtime) << endl;
    }
    //
    // --------------------------
    // now a vector exists with row locations
    // --------------------------
    MatrixXd Rls1 =MatrixXd::Zero(ntime, 1); //precomputes a series of sums used frequently in the log-liklihood calculations
    MatrixXd Rls2 =MatrixXd::Zero(ntime, totalnum); //Many are repeated due to the same risk groups and derivatives being used at mulitple points
    MatrixXd Rls3 =MatrixXd::Zero(ntime, totalnum*(totalnum+1)/2); //Sum and its derivatives are precomputed
    MatrixXd Lls1 =MatrixXd::Zero(ntime, 1); //The log-likelihood calculation has a Right and Left sum used
    MatrixXd Lls2 =MatrixXd::Zero(ntime, totalnum);
    MatrixXd Lls3 =MatrixXd::Zero(ntime, totalnum*(totalnum+1)/2);
    //The log-likelihood is calculated in parallel over the risk groups
    vector<double> Ll(totalnum,0.0); //Log-likelihood values
    vector<double> Lld(totalnum,0.0); //Log-likelihood derivative values
    vector<double> Lldd(pow(totalnum,2),0.0);//The second derivative matrix has room for every combination, but only the lower triangle is calculated initially
    //
    //
    // Calculates the side sum terms used
    if (Debug_It[3]){
        Calculate_Sides( RiskFail, RiskGroup, totalnum, ntime, R, Rd, Rdd, Rls1, Rls2, Rls3, Lls1, Lls2, Lls3,nthreads, TRUE);
    } else {
        Calculate_Sides( RiskFail, RiskGroup, totalnum, ntime, R, Rd, Rdd, Rls1, Rls2, Rls3, Lls1, Lls2, Lls3,nthreads, FALSE);
    }
    //
    //
    if (verbose){
        Rcout << "riskr checked ";
        for (int ijk=0;ijk<totalnum;ijk++){
            Rcout << Rls1.col(0).sum() << " ";
        }
        Rcout << " " << endl;
        Rcout << "risk1r checked ";
        for (int ijk=0;ijk<totalnum;ijk++){
            Rcout << Rls2.col(ijk).sum() << " ";
        }
        Rcout << " " << endl;
        Rcout << "risk2r checked ";
        for (int ijk=0;ijk<totalnum;ijk++){
            Rcout << Rls3.col(ijk*(ijk+1)/2+ijk).sum() << " ";
        }
        Rcout << " " << endl;
        //
        Rcout << "riskl checked ";
        for (int ijk=0;ijk<totalnum;ijk++){
            Rcout << Lls1.col(0).sum() << " ";
        }
        Rcout << " " << endl;
        Rcout << "risk1l checked ";
        for (int ijk=0;ijk<totalnum;ijk++){
            Rcout << Lls2.col(ijk).sum() << " ";
        }
        Rcout << " " << endl;
        Rcout << "risk2l checked ";
        for (int ijk=0;ijk<totalnum;ijk++){
            Rcout << Lls3.col(ijk*(ijk+1)/2+ijk).sum() << " ";
        }
        Rcout << " " << endl;
    }
    //
    //
    // Calculates log-likelihood
    if (Debug_It[4]){
        Calc_LogLik( nthreads, RiskFail, RiskGroup, totalnum, ntime, R, Rd, Rdd,RdR,RddR, Rls1, Rls2, Rls3, Lls1, Lls2, Lls3, Ll, Lld, Lldd, TRUE, ties_method);
    } else {
        Calc_LogLik( nthreads, RiskFail, RiskGroup, totalnum, ntime, R, Rd, Rdd,RdR,RddR, Rls1, Rls2, Rls3, Lls1, Lls2, Lls3, Ll, Lld, Lldd, FALSE, ties_method);
    }
    //
    vector <double> Ll_comp(2,Ll[0]); //vector to compare values
    //
    if (verbose){
        end_point = system_clock::now();
        ending = time_point_cast<microseconds>(end_point).time_since_epoch().count();
        Rcout <<"df100 "<<(ending-start)<<" "<<0<<" "<<0<<" "<<0<<",Calc"<<endl;//prints the time
        gibtime = system_clock::to_time_t(system_clock::now());
        Rcout << ctime(&gibtime) << endl;
        Rcout << "df101 ";//prints the log-likelihoods
        for (int ij=0;ij<totalnum;ij++){
            Rcout << Ll[ij] << " ";
        }
        Rcout << " " << endl;
        Rcout << "df102 ";//prints the first derivatives
        for (int ij=0;ij<totalnum;ij++){
            Rcout << Lld[ij] << " ";
        }
        Rcout << " " << endl;
        Rcout << "df103 ";//prints the second derivatives
        for (int ij=0;ij<totalnum;ij++){
            Rcout << Lldd[ij*totalnum+ij] << " ";
        }
        for (int ij=0;ij<totalnum;ij++){//locates highest magnitude derivative
            if (abs(Lld[ij]) > Lld_worst){
                Lld_worst = abs(Lld[ij]);
            }
        }
        Rcout << " " << endl;
        Rcout << "df104 ";//prints parameter values
        for (int ij=0;ij<totalnum;ij++){
            Rcout << beta_0[ij] << " ";
        }
        Rcout << " " << endl;
        Rcout << "df105 ";
        for (int ij=0;ij<totalnum;ij++){//prints the newton step value for zero derivative
            Rcout << Lld[ij]/Lldd[ij*totalnum+ij] << " ";
        }
        Rcout << " " << endl;
        Rcout << "df106 ";
        for (int ij=0;ij<totalnum;ij++){//prints the newton step value for zero log-likelihood
            Rcout << Ll[ij]/Lld[ij] << " ";
        }
        Rcout << " " << endl;
        Rcout << "df107 " << abs_max << " " << dose_abs_max << " " << Ll_comp[0] << " " << Ll_comp[1] << endl;//prints several convergence terms
    }
    //
    vector<double> dbeta(totalnum,0.0);
    //
    // --------------------------
    // always starts from intial guess
    // --------------------------
    vector<double> beta_p(totalnum,0.0);
    vector<double> beta_c(totalnum,0.0);
    vector<double> beta_a(totalnum,0.0);
    vector<double> beta_best(totalnum,0.0);
    VectorXd::Map(&beta_p[0], beta_0.size()) = beta_0;// stores previous parameters
    VectorXd::Map(&beta_c[0], beta_0.size()) = beta_0;// stores current parameters
    VectorXd::Map(&beta_a[0], beta_0.size()) = beta_0;// stores a refrence value for parameters
    VectorXd::Map(&beta_best[0], beta_0.size()) = beta_0;// stores the best parameters
    double Ll_best = 0.0; //a comparison log-likelihood
    int halves = 0; //number of half-steps taken
    int ind0 = fir; //used for validations
    int iteration=0; //iteration number
    //
    //
    while (iteration < maxiter){
        iteration++;
        beta_p = beta_c;//
        beta_a = beta_c;//
        beta_best = beta_c;//
        //
        //
        // Calcualtes the initial change in parameter
        if (Debug_It[5]){
            Calc_Change( nthreads, totalnum, fir, der_iden, dbeta_cap, dose_abs_max, lr, abs_max, Ll, Lld, Lldd, dbeta, change_all, tform, dint, KeepConstant, TRUE);
        } else {
            Calc_Change( nthreads, totalnum, fir, der_iden, dbeta_cap, dose_abs_max, lr, abs_max, Ll, Lld, Lldd, dbeta, change_all, tform, dint, KeepConstant, FALSE);
        }
        if (verbose){
            Rcout << "Starting Halves"<<endl;//prints the final changes for validation
        }
        //
        Ll_best = Ll[ind0];
        //
        halves=0;
        while ((Ll[ind0] <= Ll_best)&&(halves<halfmax)){ //repeats until half-steps maxed or an improvement
            beta_p = beta_c;//
            beta_a = beta_c;//
            beta_best = beta_c;//
            halves++;
            Dose = MatrixXd::Zero(df0.rows(),term_tot); //Refreshes the matrices used
            nonDose = MatrixXd::Constant(df0.rows(),term_tot,0.0);
            nonDose_LIN = MatrixXd::Constant(df0.rows(),term_tot,0.0);
            nonDose_PLIN = MatrixXd::Constant(df0.rows(),term_tot,1.0);
            nonDose_LOGLIN = MatrixXd::Constant(df0.rows(),term_tot,1.0);
            T0 = MatrixXd::Zero(df0.rows(), totalnum); 
            Td0 = MatrixXd::Zero(df0.rows(), totalnum); 
            Tdd0 = MatrixXd::Zero(df0.rows(), totalnum*(totalnum+1)/2); 
            for (int ij=0;ij<totalnum;ij++){
                beta_0[ij] = beta_a[ij] + dbeta[ij];//applies the parameter changes
                beta_c[ij] = beta_0[ij];
            }
            // ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------//
            // The same subterm, risk, sides, and log-likelihood calculations are performed every half-step and iteration
            // ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------//
            if (Debug_It[4]){
                Make_Subterms( totalnum, Term_n, tform, dfc, fir, T0, Td0, Tdd0, Dose,  nonDose, TTerm,  nonDose_LIN, nonDose_PLIN, nonDose_LOGLIN ,beta_0, df0,dint,nthreads, TRUE);
            } else {
                Make_Subterms( totalnum, Term_n, tform, dfc, fir, T0, Td0, Tdd0, Dose,  nonDose, TTerm,  nonDose_LIN, nonDose_PLIN, nonDose_LOGLIN ,beta_0, df0,dint,nthreads, FALSE);
            }
            if (verbose){
                Rcout << "values checked ";
                for (int ijk=0;ijk<totalnum;ijk++){
                    Rcout << beta_c[ijk] << " ";
                }
                Rcout << " " << endl;
                Rcout << "sums checked ";
                for (int ijk=0;ijk<totalnum;ijk++){
                    Rcout << T0.col(ijk).sum() << " ";
                }
                Rcout << " " << endl;
                Rcout << "derivs checked ";
                for (int ijk=0;ijk<totalnum;ijk++){
                    Rcout << Td0.col(ijk).sum() << " ";
                }
                Rcout << " " << endl;
                Rcout << "second derivs checked ";
                for (int ijk=0;ijk<totalnum;ijk++){
                    Rcout << Tdd0.col(ijk).sum() << " ";
                }
                Rcout << " " << endl;
                Rcout << "dose checked ";
                for (int ijk=0;ijk<term_tot;ijk++){
                    Rcout << Dose.col(ijk).array().sum() << " ";
                }
                Rcout << " " << endl;
                Rcout << "non-dose checked ";
                for (int ijk=0;ijk<term_tot;ijk++){
                    Rcout << nonDose.col(ijk).array().sum() << " ";
                }
                Rcout << " " << endl;
                Rcout << "LIN_non-dose checked ";
                for (int ijk=0;ijk<term_tot;ijk++){
                    Rcout << nonDose_LIN.col(ijk).array().sum() << " ";
                }
                Rcout << " " << endl;
                Rcout << "PLIN_non-dose checked ";
                for (int ijk=0;ijk<term_tot;ijk++){
                    Rcout << nonDose_PLIN.col(ijk).array().sum() << " ";
                }
                Rcout << " " << endl;
                Rcout << "LOGLIN_non-dose checked ";
                for (int ijk=0;ijk<term_tot;ijk++){
                    Rcout << nonDose_LOGLIN.col(ijk).array().sum() << " ";
                }
                Rcout << " " << endl;
            }
            //
            RdR = MatrixXd::Zero(df0.rows(), totalnum);
            RddR = MatrixXd::Zero(df0.rows(), totalnum*(totalnum+1)/2);
            //
            //
            if (Debug_It[6]){
                Make_Risks(modelform, tform, Term_n, totalnum, fir, T0, Td0, Tdd0, Te, R, Rd, Rdd, Dose, nonDose, TTerm, nonDose_LIN, nonDose_PLIN, nonDose_LOGLIN, RdR, RddR, nthreads, TRUE);
            } else {
                Make_Risks(modelform, tform, Term_n, totalnum, fir, T0, Td0, Tdd0, Te, R, Rd, Rdd, Dose, nonDose, TTerm, nonDose_LIN, nonDose_PLIN, nonDose_LOGLIN, RdR, RddR, nthreads, FALSE);
            }
            //
            RdR = (RdR.array().isFinite()).select(RdR,0);
            RddR = (RddR.array().isFinite()).select(RddR,0);
            //
            temp_list = List::create(_["betas"]=wrap(beta_0),_["Status"]="FAILED"); //used as a dummy return value for code checking
            if (verbose){
                Rcout << "risk checked ";
                for (int ijk=0;ijk<totalnum;ijk++){
                    Rcout << R.col(0).sum() << " ";
                }
                Rcout << " " << endl;
                Rcout << "risk1 checked ";
                for (int ijk=0;ijk<totalnum;ijk++){
                    Rcout << Rd.col(ijk).sum() << " ";
                }
                Rcout << " " << endl;
                Rcout << "risk2 checked ";
                for (int ijk=0;ijk<totalnum;ijk++){
                    Rcout << Rdd.col(ijk*(ijk+1)/2+ijk).sum() << " ";
                }
                Rcout << " " << endl;
                end_point = system_clock::now();
                ending = time_point_cast<microseconds>(end_point).time_since_epoch().count();
                Rcout <<"df100 "<<(ending-start)<<" "<<halves<<" "<<iteration<<" "<<ind0<<",Update_R"<<endl;
                //
                gibtime = system_clock::to_time_t(system_clock::now());
                Rcout << ctime(&gibtime) << endl;
            }
            fill(Ll.begin(), Ll.end(), 0.0);
            fill(Lld.begin(), Lld.end(), 0.0);
            fill(Lldd.begin(), Lldd.end(), 0.0);
            Rls1 =MatrixXd::Zero(ntime, 1); //precomputes a series of sums used frequently in the log-liklihood calculations
            Rls2 =MatrixXd::Zero(ntime, totalnum); //Many are repeated due to the same risk groups and derivatives being used at mulitple points
            Rls3 =MatrixXd::Zero(ntime, totalnum*(totalnum+1)/2);
            Lls1 =MatrixXd::Zero(ntime, 1);
            Lls2 =MatrixXd::Zero(ntime, totalnum);
            Lls3 =MatrixXd::Zero(ntime, totalnum*(totalnum+1)/2);
            if (Debug_It[7]){
                Calculate_Sides( RiskFail, RiskGroup, totalnum, ntime, R, Rd, Rdd, Rls1, Rls2, Rls3, Lls1, Lls2, Lls3,nthreads, TRUE);
            } else {
                Calculate_Sides( RiskFail, RiskGroup, totalnum, ntime, R, Rd, Rdd, Rls1, Rls2, Rls3, Lls1, Lls2, Lls3,nthreads, FALSE);
            }
            //
            if (verbose){
                Rcout << "riskr checked ";
                for (int ijk=0;ijk<totalnum;ijk++){
                    Rcout << Rls1.col(0).sum() << " ";
                }
                Rcout << " " << endl;
                Rcout << "risk1r checked ";
                for (int ijk=0;ijk<totalnum;ijk++){
                    Rcout << Rls2.col(ijk).sum() << " ";
                }
                Rcout << " " << endl;
                Rcout << "risk2r checked ";
                for (int ijk=0;ijk<totalnum;ijk++){
                    Rcout << Rls3.col(ijk*(ijk+1)/2+ijk).sum() << " ";
                }
                Rcout << " " << endl;
                //
                //
                Rcout << "riskl checked ";
                for (int ijk=0;ijk<totalnum;ijk++){
                    Rcout << Lls1.col(0).sum() << " ";
                }
                Rcout << " " << endl;
                Rcout << "risk1l checked ";
                for (int ijk=0;ijk<totalnum;ijk++){
                    Rcout << Lls2.col(ijk).sum() << " ";
                }
                Rcout << " " << endl;
                Rcout << "risk2l checked ";
                for (int ijk=0;ijk<totalnum;ijk++){
                    Rcout << Lls3.col(ijk*(ijk+1)/2+ijk).sum() << " ";
                }
                Rcout << " " << endl;
            }
            //
            if (Debug_It[8]){
                Calc_LogLik( nthreads, RiskFail, RiskGroup, totalnum, ntime, R, Rd, Rdd,RdR,RddR, Rls1, Rls2, Rls3, Lls1, Lls2, Lls3, Ll, Lld, Lldd, TRUE, ties_method);
            } else {
                Calc_LogLik( nthreads, RiskFail, RiskGroup, totalnum, ntime, R, Rd, Rdd,RdR,RddR, Rls1, Rls2, Rls3, Lls1, Lls2, Lls3, Ll, Lld, Lldd, FALSE, ties_method);
            }
            
            if (change_all){
                if (Ll[ind0] <= Ll_best){//takes a half-step if needed
                    #pragma omp parallel for num_threads(nthreads)
                    for (int ijk=0;ijk<totalnum;ijk++){
                        dbeta[ijk] = dbeta[ijk] / 2.0;
                    }
                } else{//If improved, updates the best vector
                    #pragma omp parallel for num_threads(nthreads)
                    for (int ijk=0;ijk<totalnum;ijk++){
                        beta_best[ijk] = beta_c[ijk];
                    }
                }
            } else {//For validation, the step is always carried over
                #pragma omp parallel for num_threads(nthreads)
                for (int ijk=0;ijk<totalnum;ijk++){
                    beta_best[ijk] = beta_c[ijk];
                }
            }
            if (verbose){
                end_point = system_clock::now();
                ending = time_point_cast<microseconds>(end_point).time_since_epoch().count();
                Rcout <<"df100 "<<(ending-start)<<" "<<halves<<" "<<iteration<<" "<<ind0<<",Update_calc"<<endl;
                gibtime = system_clock::to_time_t(system_clock::now());
                Rcout << ctime(&gibtime) << endl;
                Rcout << "df101 ";
                for (int ij=0;ij<totalnum;ij++){
                    Rcout << Ll[ij] << " ";
                }
                Rcout << " " << endl;
                Rcout << "df102 ";
                for (int ij=0;ij<totalnum;ij++){
                    Rcout << Lld[ij] << " ";
                }
                Rcout << " " << endl;
                Rcout << "df103 ";
                for (int ij=0;ij<totalnum;ij++){
                    Rcout << Lldd[ij*totalnum+ij] << " ";
                }
                for (int ij=0;ij<totalnum;ij++){
                    if (abs(Lld[ij]) > Lld_worst){
                        Lld_worst = abs(Lld[ij]);
                    }
                }
                Rcout << " " << endl;
                Rcout << "df104 ";
                for (int ij=0;ij<totalnum;ij++){
                    Rcout << beta_c[ij] << " ";
                }
                Rcout << " " << endl;
                Rcout << "df105 ";
                for (int ij=0;ij<totalnum;ij++){
                    Rcout << Lld[ij]/Lldd[ij*totalnum+ij] << " ";
                }
                Rcout << " " << endl;
                Rcout << "df106 ";
                for (int ij=0;ij<totalnum;ij++){
                    Rcout << Ll[ij]/Lld[ij] << " ";
                }
                Rcout << " " << endl;
                Rcout << "df107 " << abs_max << " " << dose_abs_max << " " << Ll_comp[0] << " " << Ll_comp[1] << endl;
            }
            #pragma omp parallel for num_threads(nthreads)
            for (int ijk=0;ijk<totalnum;ijk++){//totalnum*(totalnum+1)/2
                beta_0[ijk] = beta_c[ijk];
            }
            beta_best = beta_c;
        }
        if (verbose){
            end_point = system_clock::now();
            ending = time_point_cast<microseconds>(end_point).time_since_epoch().count();
            Rcout <<"df100 "<<(ending-start)<<" "<<halves<<" "<<iteration<<" "<<ind0<<",Recalc"<<endl;
            gibtime = system_clock::to_time_t(system_clock::now());
            Rcout << ctime(&gibtime) << endl;
            Rcout << "df101 ";
            for (int ij=0;ij<totalnum;ij++){
                Rcout << Ll[ij] << " ";
            }
            Rcout << " " << endl;
            Rcout << "df102 ";
            for (int ij=0;ij<totalnum;ij++){
                Rcout << Lld[ij] << " ";
            }
            Rcout << " " << endl;
            Rcout << "df103 ";
            for (int ij=0;ij<totalnum;ij++){
                Rcout << Lldd[ij*totalnum+ij] << " ";
            }
            Rcout << " " << endl;
            Rcout << "df104 ";
            for (int ij=0;ij<totalnum;ij++){
                Rcout << beta_c[ij] << " ";
            }
            Rcout << " " << endl;
            Rcout << "df105 ";
            for (int ij=0;ij<totalnum;ij++){
                Rcout << Lld[ij]/Lldd[ij*totalnum+ij] << " ";
            }
            Rcout << " " << endl;
            Rcout << "df106 ";
            for (int ij=0;ij<totalnum;ij++){
                Rcout << Ll[ij]/Lld[ij] << " ";
            }
            Rcout << " " << endl;
            Rcout << "df107 " << abs_max << " " << dose_abs_max << " " << Ll_comp[0] << " " << Ll_comp[1] << endl;
            Rcout << "Finshed iteration" << endl;
        }
    }
    //
    // Does not return anything
    //
    return;
}

//' Primary Cox PH null model function
//' \code{LogLik_Cox_PH_null} Performs the calls to calculation functions, Structures the Cox PH null model, With verbose option prints out time stamps and intermediate sums of terms and derivatives
//'
//' @param     ntime    number of event times
//' @param     df_groups    status/time matrix
//' @param     tu    vector of event times
//' @param     verbose    verbose boolean
//' @param     ties_method    tied event method
//'
//' @return List of results: Log-likelihood of optimum, AIC
// [[Rcpp::export]]
List LogLik_Cox_PH_null( int ntime, NumericMatrix df_groups, NumericVector tu, bool verbose, string ties_method){
    ;
    //
    // null model value calculation
    // --------------------------------------------------------------------------- //
    // The same steps are taken for the non-null case, with the exception of derivative calculations
    // --------------------------------------------------------------------------- //
    //
    if (verbose){
        Rcout << "START_NEW" << endl;
    }
    int totalnum=1;
    time_point<system_clock> start_point, end_point;
    start_point = system_clock::now();
    auto start = time_point_cast<microseconds>(start_point).time_since_epoch().count();
    end_point = system_clock::now();
    auto ending = time_point_cast<microseconds>(end_point).time_since_epoch().count(); //The time duration is tracked
    //
    auto gibtime = system_clock::to_time_t(system_clock::now());
    if (verbose){
        Rcout << ctime(&gibtime) << endl;
    }
    //
    Rcout.precision(10); //forces higher precision numbers printed to terminal
    int nthreads = Eigen::nbThreads()-1; //stores how many threads are allocated
    //
    // ---------------------------------------------
    // To Start, needs to seperate the derivative terms
    // ---------------------------------------------
    //
    const Map<MatrixXd> df_m(as<Map<MatrixXd> >(df_groups));
    //
    MatrixXd R = MatrixXd::Zero(df_m.rows(), 1); //preallocates matrix for Risks
    if (verbose){
        end_point = system_clock::now();
        ending = time_point_cast<microseconds>(end_point).time_since_epoch().count();
        Rcout <<"df99,"<<(ending-start)<<",Prep_Terms"<<endl;
        gibtime = system_clock::to_time_t(system_clock::now());
        Rcout << ctime(&gibtime) << endl;
    }
    //
    R = R.array() + 1.0;
    //
    if (verbose){
        Rcout << "risk checked ";
        for (int ijk=0;ijk<totalnum;ijk++){
            Rcout << R.col(0).sum() << " ";
        }
        Rcout << " " << endl;
        //
        end_point = system_clock::now();
        ending = time_point_cast<microseconds>(end_point).time_since_epoch().count();
        Rcout <<"df100 "<<(ending-start)<<" "<<0<<" "<<0<<" "<<-1<<",Prep_R"<<endl;
        gibtime = system_clock::to_time_t(system_clock::now());
        Rcout << ctime(&gibtime) << endl;
    }
    //
    // -------------------------------------------------------------------------------------------
    //
    vector<string>  RiskGroup(ntime); //vector of strings detailing the rows
    IntegerMatrix RiskFail(ntime,2); //vector giving the event rows
    //
    Make_Groups( ntime, df_m, RiskFail, RiskGroup, tu, nthreads, FALSE);
    //
    if (verbose){
        end_point = system_clock::now();
        ending = time_point_cast<microseconds>(end_point).time_since_epoch().count();
        Rcout <<"df100 "<<(ending-start)<<" "<<0<<" "<<0<<" "<<-1<<",Prep_List"<<endl;
        gibtime = system_clock::to_time_t(system_clock::now());
        Rcout << ctime(&gibtime) << endl;
    }
    //
    // --------------------------
    // now a vector exists with row locations
    // --------------------------
    MatrixXd Rls1 =MatrixXd::Zero(ntime, 1); //precomputes a series of sums used frequently in the log-liklihood calculations
    MatrixXd Lls1 =MatrixXd::Zero(ntime, 1);
    //The log-likelihood is calculated in parallel over the risk groups
    vector<double> Ll(1,0.0);
    //
    Calculate_Null_Sides( RiskFail, RiskGroup, ntime, R, Rls1, Lls1,nthreads);
    //
    //
    if (verbose){
        Rcout << "riskr checked ";
        for (int ijk=0;ijk<totalnum;ijk++){
            Rcout << Rls1.col(0).sum() << " ";
        }
        Rcout << " " << endl;
        //
        Rcout << "riskl checked ";
        for (int ijk=0;ijk<totalnum;ijk++){
            Rcout << Lls1.col(0).sum() << " ";
        }
        Rcout << " " << endl;
    }
    //
    //
    Calc_Null_LogLik( nthreads, RiskFail, RiskGroup, ntime, R, Rls1, Lls1, Ll, ties_method);
    //
    List res_list = List::create(_["LogLik"]=wrap(Ll),_["AIC"]=-2*Ll[0]);
    // returns a list of results
    return res_list;
}

//' Primary reference vector risk function
//' \code{RISK_SUBSET} Performs the calls to calculation functions, Calculates risk for a set of reference vectors, With verbose option prints out time stamps and intermediate sums of terms and derivatives
//'
//' @param     Term_n    Term numbers
//' @param     tform    subterm types
//' @param     a_n    starting values
//' @param     x_all    covariate matrix
//' @param     dfc    covariate column numbers
//' @param     fir    first term number
//' @param     modelform    model string
//' @param     verbose    verbosity boolean
//' @param     debugging    debugging boolean
//' @param     term_tot    total number of terms
//'
//' @return List of results: Risk at the reference
// [[Rcpp::export]]
NumericVector RISK_SUBSET(IntegerVector Term_n, StringVector tform, NumericVector a_n,NumericMatrix x_all,IntegerVector dfc,int fir,string modelform, bool verbose, bool debugging, int term_tot){
    ;
    //
    // Calculates the terms and risks for a reference vector
    //
    using namespace std::chrono;
    if (verbose){
        Rcout << "START_NEW" << endl;
    }
    time_point<system_clock> start_point, end_point;
    start_point = system_clock::now();
    auto start = time_point_cast<microseconds>(start_point).time_since_epoch().count();
    end_point = system_clock::now();
    auto ending = time_point_cast<microseconds>(end_point).time_since_epoch().count(); //The time duration is tracked
    //
    auto gibtime = system_clock::to_time_t(system_clock::now());
    if (verbose){
        Rcout << ctime(&gibtime) << endl;
    }
    // df0: covariate data
    const Map<MatrixXd> df0(as<Map<MatrixXd> >(x_all));
    //
    int totalnum = Term_n.size();
    //
    if (verbose){
        Rcout << "Term checked ";
        for (int ij=0;ij<totalnum;ij++){
            Rcout << Term_n[ij] << " ";
        }
        Rcout << " " << endl;
    }
    //
    //
    // cout.precision: controls the number of significant digits printed
    // nthreads: number of threads used for parallel operations
    //
    Rcout.precision(7); //forces higher precision numbers printed to terminal
    int nthreads = Eigen::nbThreads()-1; //stores how many threads are allocated
    //
    //
    if (verbose){
        end_point = system_clock::now();
        ending = time_point_cast<microseconds>(end_point).time_since_epoch().count();
        Rcout <<"df99,"<<(ending-start)<<",Starting"<<endl;
        gibtime = system_clock::to_time_t(system_clock::now());
        Rcout << ctime(&gibtime) << endl;
    }
    // ---------------------------------------------
    // To Start, needs to seperate the derivative terms
    // ---------------------------------------------
    //    int dint = 1;
    //
    // totalnum,& Term_n,  tform, dfc,& fir,& T0,& Td0,& Tdd0,& Dose,& nonDose,& beta_0,& df0, dint, nthreads,  debugging
    if (verbose){
        Rcout << "starting subterms " << term_tot << endl;
    }
    // ---------------------------------------------
    // To Start, needs to seperate the derivative terms
    // ---------------------------------------------
    //
    //
    Map<VectorXd> beta_0(as<Map<VectorXd> >(a_n));
    MatrixXd T0 = MatrixXd::Zero(df0.rows(), totalnum); //preallocates matrix for Term column
    MatrixXd Td0 = MatrixXd::Zero(df0.rows(), totalnum); //preallocates matrix for Term derivative columns
    MatrixXd Tdd0 = MatrixXd::Zero(df0.rows(), totalnum*(totalnum+1)/2); //preallocates matrix for Term second derivative columns
    //
    MatrixXd Te = MatrixXd::Zero(df0.rows(), 1); //preallocates matrix for column terms used for temporary storage
    MatrixXd R = MatrixXd::Zero(df0.rows(), 1); //preallocates matrix for Risks
    MatrixXd Rd = MatrixXd::Zero(df0.rows(), totalnum); //preallocates matrix for Risk derivatives
    MatrixXd Rdd = MatrixXd::Zero(df0.rows(), totalnum*(totalnum+1)/2); //preallocates matrix for Risk second derivatives
    //
    MatrixXd Dose = MatrixXd::Zero(df0.rows(),term_tot); //Matrix of the total dose term values
    MatrixXd nonDose = MatrixXd::Constant(df0.rows(),term_tot,1.0); //Matrix of the total non-dose term values
    MatrixXd nonDose_LIN = MatrixXd::Constant(df0.rows(),term_tot,0.0); //matrix of Linear subterm values
    MatrixXd nonDose_PLIN = MatrixXd::Constant(df0.rows(),term_tot,1.0); //matrix of Loglinear subterm values
    MatrixXd nonDose_LOGLIN = MatrixXd::Constant(df0.rows(),term_tot,1.0); //matrix of Product linear subterm values
    MatrixXd TTerm=MatrixXd::Zero(Dose.rows(),Dose.cols()); //matrix of term values
    //
    // Calculates terms
    //
    Make_Subterms( totalnum, Term_n, tform, dfc, fir, T0, Td0, Tdd0, Dose, nonDose, TTerm,  nonDose_LIN, nonDose_PLIN, nonDose_LOGLIN ,beta_0, df0,1.0,nthreads, debugging);
    MatrixXd RdR = MatrixXd::Zero(df0.rows(), totalnum); //preallocates matrix for Risk to derivative ratios
    MatrixXd RddR = MatrixXd::Zero(df0.rows(), totalnum*(totalnum+1)/2);
    //
    // Calculates risk
    //
    Make_Risks(modelform, tform, Term_n, totalnum, fir, T0, Td0, Tdd0, Te, R, Rd, Rdd, Dose, nonDose, TTerm, nonDose_LIN, nonDose_PLIN, nonDose_LOGLIN, RdR, RddR, nthreads, debugging);
    if (verbose){
        Rcout << "values checked ";
        for (int ijk=0;ijk<totalnum;ijk++){
            Rcout << beta_0[ijk] << " ";
        }
        Rcout << " " << endl;
        Rcout << "sums checked ";
        for (int ijk=0;ijk<totalnum;ijk++){
            Rcout << T0.col(ijk).sum() << " ";
        }
        Rcout << " " << endl;
        Rcout << "derivs checked ";
        for (int ijk=0;ijk<totalnum;ijk++){
            Rcout << Td0.col(ijk).sum() << " ";
        }
        Rcout << " " << endl;
        Rcout << "second derivs checked ";
        for (int ijk=0;ijk<totalnum;ijk++){
            Rcout << Tdd0.col(ijk*(ijk+1)/2+ijk).sum() << " ";
        }
        Rcout << " " << endl;
        Rcout << "dose checked ";
        for (int ijk=0;ijk<term_tot;ijk++){
            Rcout << Dose.col(ijk).array().sum() << " ";
        }
        Rcout << " " << endl;
        Rcout << "non-dose checked ";
        for (int ijk=0;ijk<term_tot;ijk++){
            Rcout << nonDose.col(ijk).array().sum() << " ";
        }
        Rcout << " " << endl;
        Rcout << "LIN_non-dose checked ";
        for (int ijk=0;ijk<term_tot;ijk++){
            Rcout << nonDose_LIN.col(ijk).array().sum() << " ";
        }
        Rcout << " " << endl;
        Rcout << "PLIN_non-dose checked ";
        for (int ijk=0;ijk<term_tot;ijk++){
            Rcout << nonDose_PLIN.col(ijk).array().sum() << " ";
        }
        Rcout << " " << endl;
        Rcout << "LOGLIN_non-dose checked ";
        for (int ijk=0;ijk<term_tot;ijk++){
            Rcout << nonDose_LOGLIN.col(ijk).array().sum() << " ";
        }
        Rcout << " " << endl;
    }
    return wrap(R.col(0));
}


