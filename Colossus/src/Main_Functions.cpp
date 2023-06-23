#include <RcppEigen.h>
#include <omp.h>
#include "Main_Functions.h"
#include "Calc_Repeated.h"
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
//' @param     double_step controls the step calculation, 0 for independent changes, 1 for solving b=Ax with complete matrices
//' @param     change_all    boolean if every parameter is being updated
//' @param     verbose    verbosity boolean
//' @param     debugging    debugging boolean
//' @param     KeepConstant    vector identifying constant parameters
//' @param     term_tot    total number of terms
//' @param     ties_method    ties method
//' @param     nthreads number of threads to use
//'
//' @return List of results: Log-likelihood of optimum, first derivative of log-likelihood, second derivative matrix, parameter list, standard deviation estimate, AIC, model information
// [[Rcpp::export]]
List LogLik_Cox_PH( IntegerVector Term_n, StringVector tform, NumericVector a_n,NumericMatrix x_all,IntegerVector dfc,int fir, int der_iden,string modelform, double lr, int maxiter, int halfmax, double epsilon, double dbeta_cap, double abs_max,double dose_abs_max, double deriv_epsilon, NumericMatrix df_groups, NumericVector tu, int double_step ,bool change_all, bool verbose, bool debugging, IntegerVector KeepConstant, int term_tot, string ties_method, int nthreads){
    ;
    //
    List temp_list = List::create(_["Status"]="FAILED"); //used as a dummy return value for code checking
    if (verbose){
        Rcout << "START_LOGLIK_COX" << endl;
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
    int reqrdnum = totalnum - sum(KeepConstant);
    //
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
    // int nthreads = Eigen::nbThreads()-1; //stores how many threads are allocated
    //
    // Lld_worst: stores the highest magnitude log-likelihood derivative
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
    MatrixXd Td0 = MatrixXd::Zero(df0.rows(), reqrdnum); //preallocates matrix for Term derivative columns
    MatrixXd Tdd0 = MatrixXd::Zero(df0.rows(), reqrdnum*(reqrdnum+1)/2); //preallocates matrix for Term second derivative columns
    //
    MatrixXd Te = MatrixXd::Zero(df0.rows(), 1); //preallocates matrix for column terms used for temporary storage
    MatrixXd R = MatrixXd::Zero(df0.rows(), 1); //preallocates matrix for Risks
    ColXd Rd = MatrixXd::Zero(df0.rows(), reqrdnum); //preallocates matrix for Risk derivatives
    ColXd Rdd = MatrixXd::Zero(df0.rows(), reqrdnum*(reqrdnum+1)/2); //preallocates matrix for Risk second derivatives
    //
    MatrixXd Dose = MatrixXd::Constant(df0.rows(),term_tot,0.0); //Matrix of the total dose term values
    MatrixXd nonDose = MatrixXd::Constant(df0.rows(),term_tot,1.0); //Matrix of the total non-dose term values
    MatrixXd nonDose_LIN = MatrixXd::Constant(df0.rows(),term_tot,0.0); //matrix of Linear subterm values
    MatrixXd nonDose_PLIN = MatrixXd::Constant(df0.rows(),term_tot,1.0); //matrix of Loglinear subterm values
    MatrixXd nonDose_LOGLIN = MatrixXd::Constant(df0.rows(),term_tot,1.0); //matrix of Product linear subterm values
    MatrixXd TTerm = MatrixXd::Zero(Dose.rows(),Dose.cols()); //matrix of term values
    double dint = dose_abs_max; //The amount of change used to calculate derivatives in threshold paramters
    double dslp = abs_max;
    //
    //
    if (verbose){
        Rcout << "starting subterms " << term_tot << endl;
    }
    //
    // Calculates the subterm and term values
    Make_subterms( totalnum, Term_n, tform, dfc, fir, T0, Td0, Tdd0, Dose, nonDose, TTerm,  nonDose_LIN, nonDose_PLIN, nonDose_LOGLIN ,beta_0, df0,dint,dslp,nthreads, debugging,KeepConstant);
    // ---------------------------------------------------------
    // Prints off a series of calculations to check at what point values are changing
    // ---------------------------------------------------------
    //
    //
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
        for (int ijk=0;ijk<reqrdnum;ijk++){
            Rcout << Td0.col(ijk).sum() << " ";
        }
        Rcout << " " << endl;
        Rcout << "second derivs checked ";
        for (int ijk=0;ijk<reqrdnum;ijk++){
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
    ColXd RdR = MatrixXd::Zero(df0.rows(), reqrdnum); //preallocates matrix for Risk to derivative ratios
    ColXd RddR = MatrixXd::Zero(df0.rows(), reqrdnum*(reqrdnum+1)/2); //preallocates matrix for Risk to second derivative ratios
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
    Make_Risks(modelform, tform, Term_n, totalnum, fir, T0, Td0, Tdd0, Te, R, Rd, Rdd, Dose, nonDose, TTerm, nonDose_LIN, nonDose_PLIN, nonDose_LOGLIN, RdR, RddR, nthreads, debugging,KeepConstant);
    //
    // Removes infinite values
    RdR = (RdR.array().isFinite()).select(RdR,0);
    RddR = (RddR.array().isFinite()).select(RddR,0);
    //
    if (R.minCoeff()<0){
        Rcout << "A non-positive risk was detected: " << R.minCoeff() << endl;
        Rcout << "final failing values ";
        for (int ijk=0;ijk<totalnum;ijk++){
            Rcout << beta_0[ijk] << " ";
        }
        Rcout << " " << endl;
        //
        Rcout << "final failing terms ";
        for (int ijk=0;ijk<totalnum;ijk++){
            Rcout << tform[ijk] << " ";
        }
        Rcout << " " << endl;
        temp_list = List::create(_["beta_0"]=wrap(beta_0) ,_["Deviation"]=R_NaN,_["Status"]="FAILED",_["LogLik"]=R_NaN);
        return temp_list;
    }
    //
    if (verbose){
        Rcout << "risk checked ";
        for (int ijk=0;ijk<1;ijk++){
            Rcout << R.col(0).sum() << " ";
        }
        Rcout << " " << endl;
        Rcout << "risk1 checked ";
        for (int ijk=0;ijk<reqrdnum;ijk++){
            Rcout << Rd.col(ijk).sum() << " ";
        }
        Rcout << " " << endl;
        Rcout << "risk2 checked ";
        for (int ijk=0;ijk<reqrdnum;ijk++){
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
    //
    MatrixXd Rls1 =MatrixXd::Zero(ntime, 1); //precomputes a series of sums used frequently in the log-liklihood calculations
    MatrixXd Rls2 =MatrixXd::Zero(ntime, reqrdnum); //Many are repeated due to the same risk groups and derivatives being used at mulitple points
    MatrixXd Rls3 =MatrixXd::Zero(ntime, reqrdnum*(reqrdnum+1)/2); //Sum and its derivatives are precomputed
    MatrixXd Lls1 =MatrixXd::Zero(ntime, 1); //The log-likelihood calculation has a Right and Left sum used
    MatrixXd Lls2 =MatrixXd::Zero(ntime, reqrdnum);
    MatrixXd Lls3 =MatrixXd::Zero(ntime, reqrdnum*(reqrdnum+1)/2);
    //The log-likelihood is calculated in parallel over the risk groups
    vector<double> Ll(reqrdnum,0.0); //Log-likelihood values
    vector<double> Lld(reqrdnum,0.0); //Log-likelihood derivative values
    vector<double> Lldd(pow(reqrdnum,2),0.0);//The second derivative matrix has room for every combination, but only the lower triangle is calculated initially
    //
    if (verbose){
        Rcout << "Made Risk Side Lists" << endl;
    }
    // Calculates the side sum terms used
    Calculate_Sides( RiskFail, RiskGroup, totalnum, ntime, R, Rd, Rdd, Rls1, Rls2, Rls3, Lls1, Lls2, Lls3,nthreads, debugging,KeepConstant);
    //
    if (verbose){
        end_point = system_clock::now();
        ending = time_point_cast<microseconds>(end_point).time_since_epoch().count();
        Rcout <<"df100 "<<(ending-start)<<" "<<0<<" "<<0<<" "<<-1<<",Prep_Sides"<<endl;
        gibtime = system_clock::to_time_t(system_clock::now());
        Rcout << ctime(&gibtime) << endl;
    }
    //
    if (verbose){
        Rcout << "riskr checked ";
        for (int ijk=0;ijk<1;ijk++){
            Rcout << Rls1.col(0).sum() << " ";
        }
        Rcout << " " << endl;
        Rcout << "risk1r checked ";
        for (int ijk=0;ijk<reqrdnum;ijk++){
            Rcout << Rls2.col(ijk).sum() << " ";
        }
        Rcout << " " << endl;
        Rcout << "risk2r checked ";
        for (int ijk=0;ijk<reqrdnum;ijk++){
            Rcout << Rls3.col(ijk*(ijk+1)/2+ijk).sum() << " ";
        }
        Rcout << " " << endl;
        //
        Rcout << "riskl checked ";
        for (int ijk=0;ijk<1;ijk++){
            Rcout << Lls1.col(0).sum() << " ";
        }
        Rcout << " " << endl;
        Rcout << "risk1l checked ";
        for (int ijk=0;ijk<reqrdnum;ijk++){
            Rcout << Lls2.col(ijk).sum() << " ";
        }
        Rcout << " " << endl;
        Rcout << "risk2l checked ";
        for (int ijk=0;ijk<reqrdnum;ijk++){
            Rcout << Lls3.col(ijk*(ijk+1)/2+ijk).sum() << " ";
        }
        Rcout << " " << endl;
    }
    //
    //
    // Calculates log-likelihood
    Calc_LogLik( nthreads, RiskFail, RiskGroup, totalnum, ntime, R, Rd, Rdd,RdR,RddR, Rls1, Rls2, Rls3, Lls1, Lls2, Lls3, Ll, Lld, Lldd, debugging, ties_method,KeepConstant);
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
        for (int ij=0;ij<reqrdnum;ij++){
            Rcout << Ll[ij] << " ";
        }
        Rcout << " " << endl;
        Rcout << "df102 ";//prints the first derivatives
        for (int ij=0;ij<reqrdnum;ij++){
            Rcout << Lld[ij] << " ";
        }
        Rcout << " " << endl;
        Rcout << "df103 ";//prints the second derivatives
        for (int ij=0;ij<reqrdnum;ij++){
            Rcout << Lldd[ij*reqrdnum+ij] << " ";
        }
        Lld_worst=0;
        for (int ij=0;ij<reqrdnum;ij++){//locates highest magnitude derivative
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
        for (int ij=0;ij<reqrdnum;ij++){//prints the newton step value for zero derivative
            Rcout << Lld[ij]/Lldd[ij*reqrdnum+ij] << " ";
        }
        Rcout << " " << endl;
        Rcout << "df106 ";
        for (int ij=0;ij<reqrdnum;ij++){//prints the newton step value for zero log-likelihood
            Rcout << Ll[ij]/Lld[ij] << " ";
        }
        Rcout << " " << endl;
        Rcout << "df107 " << double_step << " " << abs_max << " " << dose_abs_max << " " << Ll_comp[0] << " " << Ll_comp[1] << endl;//prints several convergence terms
    }
    //
    vector<double> dbeta(totalnum,0.0);
    //
    // --------------------------
    // always starts from initial guess
    // --------------------------
    vector<double> beta_c(totalnum,0.0);
    vector<double> beta_a(totalnum,0.0);
    vector<double> beta_best(totalnum,0.0);
    vector<double> beta_p(totalnum,0.0);
    VectorXd::Map(&beta_p[0], beta_0.size()) = beta_0;// stores previous parameters
    VectorXd::Map(&beta_c[0], beta_0.size()) = beta_0;// stores current parameters
    VectorXd::Map(&beta_a[0], beta_0.size()) = beta_0;// stores a refrence value for parameters
    VectorXd::Map(&beta_best[0], beta_0.size()) = beta_0;// stores the best parameters
    double halves = 0; //number of half-steps taken
    int ind0 = fir; //used for validations
    int iteration=0; //iteration number
    //
    bool convgd = FALSE;
    int iter_stop =0; //tracks if the iterations should be stopped for convergence
    int iter_check=0; //signal to check for convergence
    //
    double Ll_abs_best = 10;
    vector<double> beta_abs_best(totalnum,0.0);
    while ((iteration < maxiter)&&(iter_stop==0)){
        iteration++;
        if (verbose){
            Rcout << "Starting Iteration " << iteration <<endl;//prints the final changes for validation
        }
        beta_p = beta_c;//
        beta_a = beta_c;//
        beta_best = beta_c;//
        //
        // Calcualtes the initial change in parameter
        Calc_Change( double_step, nthreads, totalnum, fir, der_iden, dbeta_cap, dose_abs_max, lr, abs_max, Ll, Lld, Lldd, dbeta, change_all, tform, dint,dslp, KeepConstant, debugging);
        Intercept_Bound(nthreads, totalnum, beta_0, dbeta, dfc, df0, KeepConstant, debugging, tform);
        if (verbose){
            Rcout << "Starting Halves"<<endl;//prints the final changes for validation
        }
        //
        //
        if ((Ll_abs_best>0)||(Ll_abs_best < Ll[ind0])){
            Ll_abs_best = Ll[ind0];
            beta_abs_best = beta_c;
        }
        //
        if (verbose){
            Rcout << "df501 " << Ll_abs_best << endl;
            Rcout << "df504 ";//prints parameter values
            for (int ij=0;ij<totalnum;ij++){
                Rcout << beta_abs_best[ij] << " ";
            }
            Rcout << " " << endl;
        }
        halves=0;
        while ((Ll[ind0] <= Ll_abs_best)&&(halves<halfmax)){ //repeats until half-steps maxed or an improvement
            //Refreshes the matrices used
            Dose = MatrixXd::Zero(df0.rows(),term_tot);
            nonDose = MatrixXd::Constant(df0.rows(),term_tot,0.0);
            nonDose_LIN = MatrixXd::Constant(df0.rows(),term_tot,0.0);
            nonDose_PLIN = MatrixXd::Constant(df0.rows(),term_tot,1.0);
            nonDose_LOGLIN = MatrixXd::Constant(df0.rows(),term_tot,1.0);
            T0 = MatrixXd::Zero(df0.rows(), totalnum); //preallocates matrix for Term column
            Td0 = MatrixXd::Zero(df0.rows(), reqrdnum); //preallocates matrix for Term derivative columns
            Tdd0 = MatrixXd::Zero(df0.rows(), reqrdnum*(reqrdnum+1)/2); //preallocates matrix for Term second derivative columns
            for (int ij=0;ij<totalnum;ij++){
                beta_0[ij] = beta_a[ij] + dbeta[ij];
                beta_c[ij] = beta_0[ij];
            }
            // ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------//
            // The same subterm, risk, sides, and log-likelihood calculations are performed every half-step and iteration
            // ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------//
            Make_subterms( totalnum, Term_n, tform, dfc, fir, T0, Td0, Tdd0, Dose,  nonDose, TTerm,  nonDose_LIN, nonDose_PLIN, nonDose_LOGLIN ,beta_0, df0,dint,dslp,nthreads, debugging,KeepConstant);
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
                for (int ijk=0;ijk<reqrdnum;ijk++){
                    Rcout << Td0.col(ijk).sum() << " ";
                }
                Rcout << " " << endl;
                Rcout << "second derivs checked ";
                for (int ijk=0;ijk<reqrdnum;ijk++){
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
            RdR = MatrixXd::Zero(df0.rows(), reqrdnum);
            RddR = MatrixXd::Zero(df0.rows(), reqrdnum*(reqrdnum+1)/2);
            //
            //
            Make_Risks(modelform, tform, Term_n, totalnum, fir, T0, Td0, Tdd0, Te, R, Rd, Rdd, Dose, nonDose, TTerm, nonDose_LIN, nonDose_PLIN, nonDose_LOGLIN, RdR, RddR, nthreads, debugging,KeepConstant);
            //
            if (R.minCoeff()<0){
                #pragma omp parallel for num_threads(nthreads)
                for (int ijk=0;ijk<totalnum;ijk++){
                    int tij = Term_n[ijk];
                    if (TTerm.col(tij).minCoeff()<=0){
                        dbeta[ijk] = dbeta[ijk] / 2.0;
                    } else if (isinf(TTerm.col(tij).maxCoeff())){
                        dbeta[ijk] = dbeta[ijk] / 2.0;
                    }
                }
                if (verbose){
                    Rcout << "A non-positive risk was detected: " << R.minCoeff() << endl;
                }
                halves+=0.2;
            } else {
                halves++;
                RdR = (RdR.array().isFinite()).select(RdR,0);
                RddR = (RddR.array().isFinite()).select(RddR,0);
                //
                if (verbose){
                    Rcout << "risk checked ";
                    for (int ijk=0;ijk<1;ijk++){
                        Rcout << R.col(0).sum() << " ";
                    }
                    Rcout << " " << endl;
                    Rcout << "risk1 checked ";
                    for (int ijk=0;ijk<reqrdnum;ijk++){
                        Rcout << Rd.col(ijk).sum() << " ";
                    }
                    Rcout << " " << endl;
                    Rcout << "risk2 checked ";
                    for (int ijk=0;ijk<reqrdnum;ijk++){
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
                Rls2 =MatrixXd::Zero(ntime, reqrdnum); //Many are repeated due to the same risk groups and derivatives being used at mulitple points
                Rls3 =MatrixXd::Zero(ntime, reqrdnum*(reqrdnum+1)/2);
                Lls1 =MatrixXd::Zero(ntime, 1);
                Lls2 =MatrixXd::Zero(ntime, reqrdnum);
                Lls3 =MatrixXd::Zero(ntime, reqrdnum*(reqrdnum+1)/2);
                Calculate_Sides( RiskFail, RiskGroup, totalnum, ntime, R, Rd, Rdd, Rls1, Rls2, Rls3, Lls1, Lls2, Lls3,nthreads, debugging,KeepConstant);
                //
                
                if (verbose){
                    end_point = system_clock::now();
                    ending = time_point_cast<microseconds>(end_point).time_since_epoch().count();
                    Rcout <<"df100 "<<(ending-start)<<" "<<halves<<" "<<iteration<<" "<<ind0<<",Update_Sides"<<endl;
                    gibtime = system_clock::to_time_t(system_clock::now());
                    Rcout << ctime(&gibtime) << endl;
                    Rcout << "riskr checked ";
                    for (int ijk=0;ijk<1;ijk++){
                        Rcout << Rls1.col(0).sum() << " ";
                    }
                    Rcout << " " << endl;
                    Rcout << "risk1r checked ";
                    for (int ijk=0;ijk<reqrdnum;ijk++){
                        Rcout << Rls2.col(ijk).sum() << " ";
                    }
                    Rcout << " " << endl;
                    Rcout << "risk2r checked ";
                    for (int ijk=0;ijk<reqrdnum;ijk++){
                        Rcout << Rls3.col(ijk*(ijk+1)/2+ijk).sum() << " ";
                    }
                    Rcout << " " << endl;
                    //
                    //
                    Rcout << "riskl checked ";
                    for (int ijk=0;ijk<1;ijk++){
                        Rcout << Lls1.col(0).sum() << " ";
                    }
                    Rcout << " " << endl;
                    Rcout << "risk1l checked ";
                    for (int ijk=0;ijk<reqrdnum;ijk++){
                        Rcout << Lls2.col(ijk).sum() << " ";
                    }
                    Rcout << " " << endl;
                    Rcout << "risk2l checked ";
                    for (int ijk=0;ijk<reqrdnum;ijk++){
                        Rcout << Lls3.col(ijk*(ijk+1)/2+ijk).sum() << " ";
                    }
                    Rcout << " " << endl;
                }
                //
                Calc_LogLik( nthreads, RiskFail, RiskGroup, totalnum, ntime, R, Rd, Rdd,RdR,RddR, Rls1, Rls2, Rls3, Lls1, Lls2, Lls3, Ll, Lld, Lldd, debugging, ties_method,KeepConstant);
                if (verbose){
                    end_point = system_clock::now();
                    ending = time_point_cast<microseconds>(end_point).time_since_epoch().count();
                    Rcout <<"df100 "<<(ending-start)<<" "<<halves<<" "<<iteration<<" "<<ind0<<", Half"<<endl;
                    gibtime = system_clock::to_time_t(system_clock::now());
                    Rcout << ctime(&gibtime) << endl;
                    Rcout << "df101 ";
                    for (int ij=0;ij<reqrdnum;ij++){
                        Rcout << Ll[ij] << " ";
                    }
                    Rcout << " " << endl;
                    Rcout << "df102 ";
                    for (int ij=0;ij<reqrdnum;ij++){
                        Rcout << Lld[ij] << " ";
                    }
                    Rcout << " " << endl;
                    Rcout << "df103 ";
                    for (int ij=0;ij<reqrdnum;ij++){
                        Rcout << Lldd[ij*reqrdnum+ij] << " ";
                    }
                    Rcout << " " << endl;
                    Rcout << "df104 ";
                    for (int ij=0;ij<totalnum;ij++){
                        Rcout << beta_c[ij] << " ";
                    }
                    Rcout << " " << endl;
                    Rcout << "df105 ";
                    for (int ij=0;ij<reqrdnum;ij++){//prints the newton step value for zero derivative
                        Rcout << Lld[ij]/Lldd[ij*reqrdnum+ij] << " ";
                    }
                    Rcout << " " << endl;
                    Rcout << "df106 ";
                    for (int ij=0;ij<reqrdnum;ij++){//prints the newton step value for zero log-likelihood
                        Rcout << Ll[ij]/Lld[ij] << " ";
                    }
                    Rcout << " " << endl;
                    Rcout << "df107 " << double_step << " " << abs_max << " " << dose_abs_max << " " << Ll_comp[0] << " " << Ll_comp[1] << endl;//prints several convergence terms
                }
                //
                if (change_all){ //If every covariate is to be changed
                    if (Ll[ind0] <= Ll_abs_best){//takes a half-step if needed
                        #pragma omp parallel for num_threads(nthreads)
                        for (int ijk=0;ijk<totalnum;ijk++){
                            dbeta[ijk] = dbeta[ijk] * 0.5; //
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
                }
                #pragma omp parallel for num_threads(nthreads)
                for (int ijk=0;ijk<totalnum;ijk++){//totalnum*(totalnum+1)/2
                    beta_0[ijk] = beta_c[ijk];
                }
            }
        }
        if (beta_best!=beta_c){//if the risk matrices aren't the optimal values, then they must be recalculated
            if (verbose){
                Rcout << "Changing back to best"<<endl;
            }
            // If it goes through every half step without improvement, then the maximum change needs to be decreased
            abs_max = abs_max*pow(0.5,halfmax); // reduces the step sizes
            dose_abs_max = dose_abs_max*pow(0.5,halfmax);
            iter_check = 1;
            //
            beta_p = beta_best;//
            beta_a = beta_best;//
            beta_c = beta_best;//
            Dose = MatrixXd::Zero(df0.rows(),term_tot);
            nonDose = MatrixXd::Constant(df0.rows(),term_tot,0.0);
            nonDose_LIN = MatrixXd::Constant(df0.rows(),term_tot,0.0);
            nonDose_PLIN = MatrixXd::Constant(df0.rows(),term_tot,1.0);
            nonDose_LOGLIN = MatrixXd::Constant(df0.rows(),term_tot,1.0);
            for (int ij=0;ij<totalnum;ij++){
                beta_0[ij] = beta_best[ij];
            }
            Make_subterms( totalnum, Term_n, tform, dfc, fir, T0, Td0, Tdd0, Dose, nonDose, TTerm,  nonDose_LIN, nonDose_PLIN, nonDose_LOGLIN,beta_0, df0,dint,dslp,nthreads, debugging, KeepConstant);
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
                for (int ijk=0;ijk<reqrdnum;ijk++){
                    Rcout << Td0.col(ijk).sum() << " ";
                }
                Rcout << " " << endl;
                Rcout << "second derivs checked ";
                for (int ijk=0;ijk<reqrdnum;ijk++){
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
            RdR = MatrixXd::Zero(df0.rows(), reqrdnum);
            RddR = MatrixXd::Zero(df0.rows(), reqrdnum*(reqrdnum+1)/2);
            //
            //
            Make_Risks(modelform, tform, Term_n, totalnum, fir, T0, Td0, Tdd0, Te, R, Rd, Rdd, Dose, nonDose, TTerm, nonDose_LIN, nonDose_PLIN, nonDose_LOGLIN, RdR, RddR, nthreads, debugging,KeepConstant);
            R = (R.array().isFinite()).select(R,0);
            Rd = (Rd.array().isFinite()).select(Rd,0);
            Rdd = (Rdd.array().isFinite()).select(Rdd,0);
            RdR = (RdR.array().isFinite()).select(RdR,0);
            RddR = (RddR.array().isFinite()).select(RddR,0);
            //
            temp_list = List::create(_["betas"]=wrap(beta_0),_["Status"]="FAILED"); //used as a dummy return value for code checking
            if (verbose){
                Rcout << "risk checked ";
                for (int ijk=0;ijk<1;ijk++){
                    Rcout << R.col(0).sum() << " ";
                }
                Rcout << " " << endl;
                Rcout << "risk1 checked ";
                for (int ijk=0;ijk<reqrdnum;ijk++){
                    Rcout << Rd.col(ijk).sum() << " ";
                }
                Rcout << " " << endl;
                Rcout << "risk2 checked ";
                for (int ijk=0;ijk<reqrdnum;ijk++){
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
            Rls2 =MatrixXd::Zero(ntime, reqrdnum); //Many are repeated due to the same risk groups and derivatives being used at mulitple points
            Rls3 =MatrixXd::Zero(ntime, reqrdnum*(reqrdnum+1)/2);
            Lls1 =MatrixXd::Zero(ntime, 1);
            Lls2 =MatrixXd::Zero(ntime, reqrdnum);
            Lls3 =MatrixXd::Zero(ntime, reqrdnum*(reqrdnum+1)/2);
            Calculate_Sides( RiskFail, RiskGroup, totalnum, ntime, R, Rd, Rdd, Rls1, Rls2, Rls3, Lls1, Lls2, Lls3,nthreads, debugging,KeepConstant);
            //
            if (verbose){
                end_point = system_clock::now();
                ending = time_point_cast<microseconds>(end_point).time_since_epoch().count();
                Rcout <<"df100 "<<(ending-start)<<" "<<halves<<" "<<iteration<<" "<<ind0<<",Update_Sides"<<endl;
                gibtime = system_clock::to_time_t(system_clock::now());
                Rcout << ctime(&gibtime) << endl;
                Rcout << "riskr checked ";
                for (int ijk=0;ijk<1;ijk++){
                    Rcout << Rls1.col(0).sum() << " ";
                }
                Rcout << " " << endl;
                Rcout << "risk1r checked ";
                for (int ijk=0;ijk<reqrdnum;ijk++){
                    Rcout << Rls2.col(ijk).sum() << " ";
                }
                Rcout << " " << endl;
                Rcout << "risk2r checked ";
                for (int ijk=0;ijk<reqrdnum;ijk++){
                    Rcout << Rls3.col(ijk*(ijk+1)/2+ijk).sum() << " ";
                }
                Rcout << " " << endl;
                //
                //
                Rcout << "riskl checked ";
                for (int ijk=0;ijk<1;ijk++){
                    Rcout << Lls1.col(0).sum() << " ";
                }
                Rcout << " " << endl;
                Rcout << "risk1l checked ";
                for (int ijk=0;ijk<reqrdnum;ijk++){
                    Rcout << Lls2.col(ijk).sum() << " ";
                }
                Rcout << " " << endl;
                Rcout << "risk2l checked ";
                for (int ijk=0;ijk<reqrdnum;ijk++){
                    Rcout << Lls3.col(ijk*(ijk+1)/2+ijk).sum() << " ";
                }
                Rcout << " " << endl;
            }
            //
            Calc_LogLik( nthreads, RiskFail, RiskGroup, totalnum, ntime, R, Rd, Rdd,RdR,RddR, Rls1, Rls2, Rls3, Lls1, Lls2, Lls3, Ll, Lld, Lldd, debugging, ties_method,KeepConstant);
        }
        Lld_worst=0;
        for (int ij=0;ij<reqrdnum;ij++){
            if (abs(Lld[ij]) > Lld_worst){
                Lld_worst = abs(Lld[ij]);
            }
        }
        if (iteration > reqrdnum){//Doesn't check the first several iterations for convergence
            if ((iteration % (reqrdnum))||(iter_check==1)){//Checks every set number of iterations
                iter_check=0;
                if (Lld_worst < deriv_epsilon){//ends if the derivatives are low enough
                    iter_stop = 1;
                    convgd = TRUE;
                }
                Ll_comp[1]=Ll[0];
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
            for (int ij=0;ij<reqrdnum;ij++){
                Rcout << Ll[ij] << " ";
            }
            Rcout << " " << endl;
            Rcout << "df102 ";
            for (int ij=0;ij<reqrdnum;ij++){
                Rcout << Lld[ij] << " ";
            }
            Rcout << " " << endl;
            Rcout << "df103 ";
            for (int ij=0;ij<reqrdnum;ij++){
                Rcout << Lldd[ij*reqrdnum+ij] << " ";
            }
            Rcout << " " << endl;
            Rcout << "df104 ";
            for (int ij=0;ij<totalnum;ij++){
                Rcout << beta_c[ij] << " ";
            }
            Rcout << " " << endl;
            Rcout << "df105 ";
            for (int ij=0;ij<reqrdnum;ij++){//prints the newton step value for zero derivative
                Rcout << Lld[ij]/Lldd[ij*reqrdnum+ij] << " ";
            }
            Rcout << " " << endl;
            Rcout << "df106 ";
            for (int ij=0;ij<reqrdnum;ij++){//prints the newton step value for zero log-likelihood
                Rcout << Ll[ij]/Lld[ij] << " ";
            }
            Rcout << " " << endl;
            Rcout << "df107 " << double_step << " " << abs_max << " " << dose_abs_max << " " << Ll_comp[0] << " " << Ll_comp[1] << endl;//prints several convergence terms
            Rcout << "Finshed iteration" << endl;
        }
    }
    // -----------------------------------------------
    // Performing Full Calculation to get full second derivative matrix
    // -----------------------------------------------
    fill(Ll.begin(), Ll.end(), 0.0);
    fill(Lld.begin(), Lld.end(), 0.0);
    fill(Lldd.begin(), Lldd.end(), 0.0);
    Calculate_Sides( RiskFail, RiskGroup, totalnum, ntime, R, Rd, Rdd, Rls1, Rls2, Rls3, Lls1, Lls2, Lls3,nthreads, debugging,KeepConstant);
    //
    if (verbose){
        end_point = system_clock::now();
        ending = time_point_cast<microseconds>(end_point).time_since_epoch().count();
        Rcout <<"df100 "<<(ending-start)<<" "<<halves<<" "<<iteration<<" "<<ind0<<",Update_Sides"<<endl;
        gibtime = system_clock::to_time_t(system_clock::now());
        Rcout << ctime(&gibtime) << endl;
        Rcout << "Wrapping up" << endl;
    }
    //
    Calc_LogLik( nthreads, RiskFail, RiskGroup, totalnum, ntime, R, Rd, Rdd,RdR,RddR, Rls1, Rls2, Rls3, Lls1, Lls2, Lls3, Ll, Lld, Lldd, debugging, ties_method,KeepConstant);
    //
    if ((Ll_abs_best>0)||(Ll_abs_best < Ll[ind0])){
        Ll_abs_best = Ll[ind0];
        beta_abs_best = beta_c;
    }
    //
    if (verbose){
        Rcout << "df501 " << Ll_abs_best << endl;
        Rcout << "df504 ";//prints parameter values
        for (int ij=0;ij<totalnum;ij++){
            Rcout << beta_abs_best[ij] << " ";
        }
        Rcout << " " << endl;
    }
    //
    List para_list = List::create(_["Term_n"]=Term_n,_["tforms"]=tform); //stores the term information
    List control_list = List::create(_["Iteration"]=iteration, _["Maximum Step"]=abs_max, _["Derivative Limiting"]=Lld_worst); //stores the total number of iterations used
    //
    int kept_covs = totalnum - sum(KeepConstant); //does not base the standard deviation off of constant parameters
    NumericVector Lldd_vec(kept_covs * kept_covs);
    #pragma omp parallel for schedule(dynamic) num_threads(nthreads)
    for (int ijk=0;ijk<totalnum*(totalnum+1)/2;ijk++){
        int ij = 0;
        int jk = ijk;
        int pij_ind=-100;
        int pjk_ind=-100;
        while (jk>ij){
            ij++;
            jk-=ij;
        }
        if (KeepConstant[ij]==0){
            pij_ind = ij - sum(head(KeepConstant,ij));
            if (KeepConstant[jk]==0){
                pjk_ind = jk - sum(head(KeepConstant,jk));
                Lldd_vec[pij_ind * kept_covs + pjk_ind]=Lldd[pij_ind * kept_covs + pjk_ind];
            }
        }
    }
    for (int ijk=0;ijk<kept_covs*(kept_covs+1)/2;ijk++){
        int ij = 0;
        int jk = ijk;
        while (jk>ij){
            ij++;
            jk-=ij;
        }
        Lldd_vec[jk * kept_covs + ij]=Lldd_vec[ij * kept_covs + jk];
    }
    Lldd_vec.attr("dim") = Dimension(kept_covs, kept_covs);
    const Map<MatrixXd> Lldd_mat(as<Map<MatrixXd> >(Lldd_vec));
    //
    MatrixXd Lldd_inv = -1 * Lldd_mat.inverse().matrix(); //uses inverse information matrix to calculate the standard deviation
    VectorXd stdev = VectorXd::Zero(totalnum);
    for (int ij=0;ij<totalnum;ij++){
        if (KeepConstant[ij]==0){
            int pij_ind = ij - sum(head(KeepConstant,ij));
            stdev(ij) = sqrt(Lldd_inv(pij_ind,pij_ind));
        }
    }
    //
    //
    List res_list = List::create(_["LogLik"]=wrap(Ll[0]),_["First_Der"]=wrap(Lld),_["Second_Der"]=Lldd_vec,_["beta_0"]=wrap(beta_0) ,_["Standard_Deviation"]=wrap(stdev) ,_["AIC"]=2*(totalnum-accumulate(KeepConstant.begin(),KeepConstant.end(), 0.0))-2*Ll[0],_["Parameter_Lists"]=para_list,_["Control_List"]=control_list,_["Convgerged"]=convgd);
    // returns a list of results
    return res_list;
}

//' Primary Cox PH regression with basic model
//' \code{LogLik_Cox_PH_basic} Performs the calls to calculation functions, Structures the Cox PH regression, With verbose option prints out time stamps and intermediate sums of terms and derivatives
//'
//' @param     a_n    starting values
//' @param     x_all    covariate matrix
//' @param     dfc    covariate column numbers
//' @param     der_iden    subterm number for derivative tests
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
//' @param     double_step controls the step calculation, 0 for independent changes, 1 for solving b=Ax with complete matrices
//' @param     change_all    boolean if every parameter is being updated
//' @param     verbose    verbosity boolean
//' @param     debugging    debugging boolean
//' @param     KeepConstant    vector identifying constant parameters
//' @param     ties_method    ties method
//' @param     nthreads number of threads to use
//'
//' @return List of results: Log-likelihood of optimum, first derivative of log-likelihood, second derivative matrix, parameter list, standard deviation estimate, AIC, model information
// [[Rcpp::export]]
List LogLik_Cox_PH_basic( NumericVector a_n,NumericMatrix x_all,IntegerVector dfc, int der_iden, double lr, int maxiter, int halfmax, double epsilon, double dbeta_cap, double abs_max,double dose_abs_max, double deriv_epsilon, NumericMatrix df_groups, NumericVector tu, int double_step ,bool change_all, bool verbose, bool debugging, IntegerVector KeepConstant, string ties_method, int nthreads){
    ;
    //
    List temp_list = List::create(_["Status"]="FAILED"); //used as a dummy return value for code checking
    if (verbose){
        Rcout << "START_LOGLIK_COX_BASIC" << endl;
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
    int totalnum = a_n.size();
    int reqrdnum = totalnum - sum(KeepConstant);
    //
    // cout.precision: controls the number of significant digits printed
    // nthreads: number of threads used for parallel operations
    //
    Rcout.precision(7); //forces higher precision numbers printed to terminal
    // int nthreads = Eigen::nbThreads()-1; //stores how many threads are allocated
    //
    int First_NonConst = -1;
    for (int i=0;i<totalnum;i++){
        if (KeepConstant[i]==0){
            First_NonConst=i;
            break;
        }
    }
    //
    //
    //
    // Lld_worst: stores the highest magnitude log-likelihood derivative
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
    //
    MatrixXd R = MatrixXd::Zero(df0.rows(), 1); //preallocates matrix for Risks
    ColXd Rd = MatrixXd::Zero(df0.rows(), reqrdnum); //preallocates matrix for Risk derivatives
    ColXd Rdd = MatrixXd::Zero(df0.rows(), reqrdnum*(reqrdnum+1)/2); //preallocates matrix for Risk second derivatives
    //
    //
    // Calculates the subterm and term values
    Make_subterms_Basic( totalnum,  dfc,  T0 ,beta_0, df0,nthreads, debugging);
    // ---------------------------------------------------------
    // Prints off a series of calculations to check at what point values are changing
    // ---------------------------------------------------------
    //
    //
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
    }
    //
    ColXd RdR = MatrixXd::Zero(df0.rows(), reqrdnum); //preallocates matrix for Risk to derivative ratios
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
    Make_Risks_Basic(totalnum, T0, R, Rd, Rdd, RdR, nthreads, debugging,df0,dfc,KeepConstant);
    //
    // Removes infinite values
    RdR = (RdR.array().isFinite()).select(RdR,0);
    //
    if (R.minCoeff()<0){
        Rcout << "A non-positive risk was detected: " << R.minCoeff() << endl;
        Rcout << "final failing values ";
        for (int ijk=0;ijk<totalnum;ijk++){
            Rcout << beta_0[ijk] << " ";
        }
        Rcout << " " << endl;
        //
        Rcout << " " << endl;
        temp_list = List::create(_["beta_0"]=wrap(beta_0) ,_["Deviation"]=R_NaN,_["Status"]="FAILED",_["LogLik"]=R_NaN);
        return temp_list;
    }
    //
    if (verbose){
        Rcout << "risk checked ";
        for (int ijk=0;ijk<1;ijk++){
            Rcout << R.col(0).sum() << " ";
        }
        Rcout << " " << endl;
        Rcout << "risk1 checked ";
        for (int ijk=0;ijk<reqrdnum;ijk++){
            Rcout << Rd.col(ijk).sum() << " ";
        }
        Rcout << " " << endl;
        Rcout << "risk2 checked ";
        for (int ijk=0;ijk<reqrdnum;ijk++){
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
    //
    MatrixXd Rls1 =MatrixXd::Zero(ntime, 1); //precomputes a series of sums used frequently in the log-liklihood calculations
    MatrixXd Rls2 =MatrixXd::Zero(ntime, reqrdnum); //Many are repeated due to the same risk groups and derivatives being used at mulitple points
    MatrixXd Rls3 =MatrixXd::Zero(ntime, reqrdnum*(reqrdnum+1)/2); //Sum and its derivatives are precomputed
    MatrixXd Lls1 =MatrixXd::Zero(ntime, 1); //The log-likelihood calculation has a Right and Left sum used
    MatrixXd Lls2 =MatrixXd::Zero(ntime, reqrdnum);
    MatrixXd Lls3 =MatrixXd::Zero(ntime, reqrdnum*(reqrdnum+1)/2);
    //The log-likelihood is calculated in parallel over the risk groups
    vector<double> Ll(reqrdnum,0.0); //Log-likelihood values
    vector<double> Lld(reqrdnum,0.0); //Log-likelihood derivative values
    vector<double> Lldd(pow(reqrdnum,2),0.0);//The second derivative matrix has room for every combination, but only the lower triangle is calculated initially
    //
    if (verbose){
        Rcout << "Made Risk Side Lists" << endl;
    }
    // Calculates the side sum terms used
    Calculate_Sides( RiskFail, RiskGroup, totalnum, ntime, R, Rd, Rdd, Rls1, Rls2, Rls3, Lls1, Lls2, Lls3,nthreads, debugging,KeepConstant);
    //
    //
    if (verbose){
        end_point = system_clock::now();
        ending = time_point_cast<microseconds>(end_point).time_since_epoch().count();
        Rcout <<"df100 "<<(ending-start)<<" "<<0<<" "<<0<<" "<<-1<<",Prep_Sides"<<endl;
        gibtime = system_clock::to_time_t(system_clock::now());
        Rcout << ctime(&gibtime) << endl;
        Rcout << "riskr checked ";
        for (int ijk=0;ijk<1;ijk++){
            Rcout << Rls1.col(0).sum() << " ";
        }
        Rcout << " " << endl;
        Rcout << "risk1r checked ";
        for (int ijk=0;ijk<reqrdnum;ijk++){
            Rcout << Rls2.col(ijk).sum() << " ";
        }
        Rcout << " " << endl;
        Rcout << "risk2r checked ";
        for (int ijk=0;ijk<reqrdnum;ijk++){
            Rcout << Rls3.col(ijk*(ijk+1)/2+ijk).sum() << " ";
        }
        Rcout << " " << endl;
        //
        Rcout << "riskl checked ";
        for (int ijk=0;ijk<1;ijk++){
            Rcout << Lls1.col(0).sum() << " ";
        }
        Rcout << " " << endl;
        Rcout << "risk1l checked ";
        for (int ijk=0;ijk<reqrdnum;ijk++){
            Rcout << Lls2.col(ijk).sum() << " ";
        }
        Rcout << " " << endl;
        Rcout << "risk2l checked ";
        for (int ijk=0;ijk<reqrdnum;ijk++){
            Rcout << Lls3.col(ijk*(ijk+1)/2+ijk).sum() << " ";
        }
        Rcout << " " << endl;
    }
    //
    //
    // Calculates log-likelihood
    Calc_LogLik_Basic( nthreads, RiskFail, RiskGroup, totalnum, ntime, R, Rd, Rdd,RdR, Rls1, Rls2, Rls3, Lls1, Lls2, Lls3, Ll, Lld, Lldd, debugging, ties_method,KeepConstant);
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
        for (int ij=0;ij<reqrdnum;ij++){
            Rcout << Ll[ij] << " ";
        }
        Rcout << " " << endl;
        Rcout << "df102 ";//prints the first derivatives
        for (int ij=0;ij<reqrdnum;ij++){
            Rcout << Lld[ij] << " ";
        }
        Rcout << " " << endl;
        Rcout << "df103 ";//prints the second derivatives
        for (int ij=0;ij<reqrdnum;ij++){
            Rcout << Lldd[ij*reqrdnum+ij] << " ";
        }
        Lld_worst=0;
        for (int ij=0;ij<reqrdnum;ij++){//locates highest magnitude derivative
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
        for (int ij=0;ij<reqrdnum;ij++){//prints the newton step value for zero derivative
            Rcout << Lld[ij]/Lldd[ij*reqrdnum+ij] << " ";
        }
        Rcout << " " << endl;
        Rcout << "df106 ";
        for (int ij=0;ij<reqrdnum;ij++){//prints the newton step value for zero log-likelihood
            Rcout << Ll[ij]/Lld[ij] << " ";
        }
        Rcout << " " << endl;
        Rcout << "df107 " << double_step << " " << abs_max << " " << dose_abs_max << " " << Ll_comp[0] << " " << Ll_comp[1] << endl;//prints several convergence terms
    }
    //
    vector<double> dbeta(totalnum,0.0);
    //
    // --------------------------
    // always starts from initial guess
    // --------------------------
    vector<double> beta_c(totalnum,0.0);
    vector<double> beta_a(totalnum,0.0);
    vector<double> beta_best(totalnum,0.0);
    vector<double> beta_p(totalnum,0.0);
    VectorXd::Map(&beta_p[0], beta_0.size()) = beta_0;// stores previous parameters
    VectorXd::Map(&beta_c[0], beta_0.size()) = beta_0;// stores current parameters
    VectorXd::Map(&beta_a[0], beta_0.size()) = beta_0;// stores a refrence value for parameters
    VectorXd::Map(&beta_best[0], beta_0.size()) = beta_0;// stores the best parameters
    double halves = 0; //number of half-steps taken
    int ind0 = 0; //used for validations
    int iteration=0; //iteration number
    //
    int iter_stop =0; //tracks if the iterations should be stopped for convergence
    int iter_check=0; //signal to check for convergence
    bool convgd = FALSE;
    //
    double Ll_abs_best = 10;
    vector<double> beta_abs_best(totalnum,0.0);
    while ((iteration < maxiter)&&(iter_stop==0)){
        iteration++;
        beta_p = beta_c;//
        beta_a = beta_c;//
        beta_best = beta_c;//
        //
        // Calcualtes the initial change in parameter
        Calc_Change_Basic( double_step, nthreads, totalnum, der_iden, dbeta_cap, lr, abs_max, Ll, Lld, Lldd, dbeta, change_all, KeepConstant, debugging);
        if (verbose){
            Rcout << "Starting Halves"<<endl;//prints the final changes for validation
        }
        //
        //
        if ((Ll_abs_best>0)||(Ll_abs_best < Ll[ind0])){
            Ll_abs_best = Ll[ind0];
            beta_abs_best = beta_c;
        }
        //
        if (verbose){
            Rcout << "df501 " << Ll_abs_best << endl;
            Rcout << "df504 ";//prints parameter values
            for (int ij=0;ij<totalnum;ij++){
                Rcout << beta_abs_best[ij] << " ";
            }
            Rcout << " " << endl;
        }
        halves=0;
        while ((Ll[ind0] <= Ll_abs_best)&&(halves<halfmax)){ //repeats until half-steps maxed or an improvement
            //Refreshes the matrices used
            T0 = MatrixXd::Zero(df0.rows(), totalnum); //preallocates matrix for Term column
            for (int ij=0;ij<totalnum;ij++){
                beta_0[ij] = beta_a[ij] + dbeta[ij];
                beta_c[ij] = beta_0[ij];
            }
            // ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------//
            // The same subterm, risk, sides, and log-likelihood calculations are performed every half-step and iteration
            // ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------//
            Make_subterms_Basic( totalnum,  dfc,  T0 ,beta_0, df0,nthreads, debugging);
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
            }
            //
            RdR = MatrixXd::Zero(df0.rows(), reqrdnum);
            //
            //
            Make_Risks_Basic(totalnum, T0, R, Rd, Rdd, RdR, nthreads, debugging,df0,dfc,KeepConstant);
            //
            halves++;
            RdR = (RdR.array().isFinite()).select(RdR,0);
            //
            if (verbose){
                Rcout << "risk checked ";
                for (int ijk=0;ijk<1;ijk++){
                    Rcout << R.col(0).sum() << " ";
                }
                Rcout << " " << endl;
                Rcout << "risk1 checked ";
                for (int ijk=0;ijk<reqrdnum;ijk++){
                    Rcout << Rd.col(ijk).sum() << " ";
                }
                Rcout << " " << endl;
                Rcout << "risk2 checked ";
                for (int ijk=0;ijk<reqrdnum;ijk++){
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
            Rls2 =MatrixXd::Zero(ntime, reqrdnum); //Many are repeated due to the same risk groups and derivatives being used at mulitple points
            Rls3 =MatrixXd::Zero(ntime, reqrdnum*(reqrdnum+1)/2);
            Lls1 =MatrixXd::Zero(ntime, 1);
            Lls2 =MatrixXd::Zero(ntime, reqrdnum);
            Lls3 =MatrixXd::Zero(ntime, reqrdnum*(reqrdnum+1)/2);
            Calculate_Sides( RiskFail, RiskGroup, totalnum, ntime, R, Rd, Rdd, Rls1, Rls2, Rls3, Lls1, Lls2, Lls3,nthreads, debugging,KeepConstant);
            //
            if (verbose){
                end_point = system_clock::now();
                ending = time_point_cast<microseconds>(end_point).time_since_epoch().count();
                Rcout <<"df100 "<<(ending-start)<<" "<<halves<<" "<<iteration<<" "<<ind0<<",Update_Sides"<<endl;
                //
                gibtime = system_clock::to_time_t(system_clock::now());
                Rcout << ctime(&gibtime) << endl;
                Rcout << "riskr checked ";
                for (int ijk=0;ijk<1;ijk++){
                    Rcout << Rls1.col(0).sum() << " ";
                }
                Rcout << " " << endl;
                Rcout << "risk1r checked ";
                for (int ijk=0;ijk<reqrdnum;ijk++){
                    Rcout << Rls2.col(ijk).sum() << " ";
                }
                Rcout << " " << endl;
                Rcout << "risk2r checked ";
                for (int ijk=0;ijk<reqrdnum;ijk++){
                    Rcout << Rls3.col(ijk*(ijk+1)/2+ijk).sum() << " ";
                }
                Rcout << " " << endl;
                //
                //
                Rcout << "riskl checked ";
                for (int ijk=0;ijk<1;ijk++){
                    Rcout << Lls1.col(0).sum() << " ";
                }
                Rcout << " " << endl;
                Rcout << "risk1l checked ";
                for (int ijk=0;ijk<reqrdnum;ijk++){
                    Rcout << Lls2.col(ijk).sum() << " ";
                }
                Rcout << " " << endl;
                Rcout << "risk2l checked ";
                for (int ijk=0;ijk<reqrdnum;ijk++){
                    Rcout << Lls3.col(ijk*(ijk+1)/2+ijk).sum() << " ";
                }
                Rcout << " " << endl;
            }
            //
            Calc_LogLik_Basic( nthreads, RiskFail, RiskGroup, totalnum, ntime, R, Rd, Rdd,RdR, Rls1, Rls2, Rls3, Lls1, Lls2, Lls3, Ll, Lld, Lldd, debugging, ties_method,KeepConstant);
            if (verbose){
                end_point = system_clock::now();
                ending = time_point_cast<microseconds>(end_point).time_since_epoch().count();
                Rcout <<"df100 "<<(ending-start)<<" "<<0<<" "<<0<<" "<<0<<",Half Calc"<<endl;//prints the time
                gibtime = system_clock::to_time_t(system_clock::now());
                Rcout << ctime(&gibtime) << endl;
                Rcout << "df101 ";//prints the log-likelihoods
                for (int ij=0;ij<reqrdnum;ij++){
                    Rcout << Ll[ij] << " ";
                }
                Rcout << " " << endl;
                Rcout << "df102 ";//prints the first derivatives
                for (int ij=0;ij<reqrdnum;ij++){
                    Rcout << Lld[ij] << " ";
                }
                Rcout << " " << endl;
                Rcout << "df103 ";//prints the second derivatives
                for (int ij=0;ij<reqrdnum;ij++){
                    Rcout << Lldd[ij*reqrdnum+ij] << " ";
                }
                Lld_worst=0;
                for (int ij=0;ij<reqrdnum;ij++){//locates highest magnitude derivative
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
                for (int ij=0;ij<reqrdnum;ij++){//prints the newton step value for zero derivative
                    Rcout << Lld[ij]/Lldd[ij*reqrdnum+ij] << " ";
                }
                Rcout << " " << endl;
                Rcout << "df106 ";
                for (int ij=0;ij<reqrdnum;ij++){//prints the newton step value for zero log-likelihood
                    Rcout << Ll[ij]/Lld[ij] << " ";
                }
                Rcout << " " << endl;
                Rcout << "df107 " << double_step << " " << abs_max << " " << dose_abs_max << " " << Ll_comp[0] << " " << Ll_comp[1] << endl;//prints several convergence terms
            }
            if (change_all){ //If every covariate is to be changed
                if (Ll[ind0] <= Ll_abs_best){//takes a half-step if needed
                    #pragma omp parallel for num_threads(nthreads)
                    for (int ijk=0;ijk<totalnum;ijk++){
                        dbeta[ijk] = dbeta[ijk] * 0.5; //
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
            // If it goes through every half step without improvement, then the maximum change needs to be decreased
            abs_max = abs_max*pow(0.5,halfmax); // reduces the step sizes
            dose_abs_max = dose_abs_max*pow(0.5,halfmax);
            iter_check = 1;
            //
            beta_p = beta_best;//
            beta_a = beta_best;//
            beta_c = beta_best;//
            for (int ij=0;ij<totalnum;ij++){
                beta_0[ij] = beta_best[ij];
            }
            Make_subterms_Basic( totalnum,  dfc,  T0 ,beta_0, df0,nthreads, debugging);
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
            }
            //
            RdR = MatrixXd::Zero(df0.rows(), reqrdnum);
            //
            //
            Make_Risks_Basic(totalnum, T0, R, Rd, Rdd, RdR, nthreads, debugging,df0,dfc,KeepConstant);
            R = (R.array().isFinite()).select(R,0);
            Rd = (Rd.array().isFinite()).select(Rd,0);
            Rdd = (Rdd.array().isFinite()).select(Rdd,0);
            RdR = (RdR.array().isFinite()).select(RdR,0);
            //
            temp_list = List::create(_["betas"]=wrap(beta_0),_["Status"]="FAILED"); //used as a dummy return value for code checking
            if (verbose){
                Rcout << "risk checked ";
                for (int ijk=0;ijk<1;ijk++){
                    Rcout << R.col(0).sum() << " ";
                }
                Rcout << " " << endl;
                Rcout << "risk1 checked ";
                for (int ijk=0;ijk<reqrdnum;ijk++){
                    Rcout << Rd.col(ijk).sum() << " ";
                }
                Rcout << " " << endl;
                Rcout << "risk2 checked ";
                for (int ijk=0;ijk<reqrdnum;ijk++){
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
            Rls2 =MatrixXd::Zero(ntime, reqrdnum); //Many are repeated due to the same risk groups and derivatives being used at mulitple points
            Rls3 =MatrixXd::Zero(ntime, reqrdnum*(reqrdnum+1)/2);
            Lls1 =MatrixXd::Zero(ntime, 1);
            Lls2 =MatrixXd::Zero(ntime, reqrdnum);
            Lls3 =MatrixXd::Zero(ntime, reqrdnum*(reqrdnum+1)/2);
            Calculate_Sides( RiskFail, RiskGroup, totalnum, ntime, R, Rd, Rdd, Rls1, Rls2, Rls3, Lls1, Lls2, Lls3,nthreads, debugging,KeepConstant);
            //
            if (verbose){
                end_point = system_clock::now();
                ending = time_point_cast<microseconds>(end_point).time_since_epoch().count();
                Rcout <<"df100 "<<(ending-start)<<" "<<halves<<" "<<iteration<<" "<<ind0<<",Update_Sides"<<endl;
                //
                gibtime = system_clock::to_time_t(system_clock::now());
                Rcout << ctime(&gibtime) << endl;
                Rcout << "riskr checked ";
                for (int ijk=0;ijk<1;ijk++){
                    Rcout << Rls1.col(0).sum() << " ";
                }
                Rcout << " " << endl;
                Rcout << "risk1r checked ";
                for (int ijk=0;ijk<reqrdnum;ijk++){
                    Rcout << Rls2.col(ijk).sum() << " ";
                }
                Rcout << " " << endl;
                Rcout << "risk2r checked ";
                for (int ijk=0;ijk<reqrdnum;ijk++){
                    Rcout << Rls3.col(ijk*(ijk+1)/2+ijk).sum() << " ";
                }
                Rcout << " " << endl;
                //
                //
                Rcout << "riskl checked ";
                for (int ijk=0;ijk<1;ijk++){
                    Rcout << Lls1.col(0).sum() << " ";
                }
                Rcout << " " << endl;
                Rcout << "risk1l checked ";
                for (int ijk=0;ijk<reqrdnum;ijk++){
                    Rcout << Lls2.col(ijk).sum() << " ";
                }
                Rcout << " " << endl;
                Rcout << "risk2l checked ";
                for (int ijk=0;ijk<reqrdnum;ijk++){
                    Rcout << Lls3.col(ijk*(ijk+1)/2+ijk).sum() << " ";
                }
                Rcout << " " << endl;
            }
            //
            Calc_LogLik_Basic( nthreads, RiskFail, RiskGroup, totalnum, ntime, R, Rd, Rdd,RdR, Rls1, Rls2, Rls3, Lls1, Lls2, Lls3, Ll, Lld, Lldd, debugging, ties_method,KeepConstant);
        }
        Lld_worst=0;
        for (int ij=0;ij<reqrdnum;ij++){
            if (abs(Lld[ij]) > Lld_worst){
                Lld_worst = abs(Lld[ij]);
            }
        }
        if (iteration > reqrdnum){//Doesn't check the first several iterations for convergence
            if ((iteration % (reqrdnum))||(iter_check==1)){//Checks every set number of iterations
                iter_check=0;
                if (Lld_worst < deriv_epsilon){//ends if the derivatives are low enough
                    iter_stop = 1;
                    convgd = TRUE;
                }
                Ll_comp[1]=Ll[0];
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
            for (int ij=0;ij<reqrdnum;ij++){
                Rcout << Ll[ij] << " ";
            }
            Rcout << " " << endl;
            Rcout << "df102 ";
            for (int ij=0;ij<reqrdnum;ij++){
                Rcout << Lld[ij] << " ";
            }
            Rcout << " " << endl;
            Rcout << "df103 ";
            for (int ij=0;ij<reqrdnum;ij++){
                Rcout << Lldd[ij*reqrdnum+ij] << " ";
            }
            Rcout << " " << endl;
            Rcout << "df104 ";
            for (int ij=0;ij<totalnum;ij++){
                Rcout << beta_c[ij] << " ";
            }
            Rcout << " " << endl;
            Rcout << "df105 ";
            for (int ij=0;ij<reqrdnum;ij++){//prints the newton step value for zero derivative
                Rcout << Lld[ij]/Lldd[ij*reqrdnum+ij] << " ";
            }
            Rcout << " " << endl;
            Rcout << "df106 ";
            for (int ij=0;ij<reqrdnum;ij++){//prints the newton step value for zero log-likelihood
                Rcout << Ll[ij]/Lld[ij] << " ";
            }
            Rcout << " " << endl;
            Rcout << "df107 " << double_step << " " << abs_max << " " << dose_abs_max << " " << Ll_comp[0] << " " << Ll_comp[1] << endl;//prints several convergence terms
            Rcout << "Finshed iteration" << endl;
        }
    }
    // -----------------------------------------------
    // Performing Full Calculation to get full second derivative matrix
    // -----------------------------------------------
    fill(Ll.begin(), Ll.end(), 0.0);
    fill(Lld.begin(), Lld.end(), 0.0);
    fill(Lldd.begin(), Lldd.end(), 0.0);
    Calculate_Sides( RiskFail, RiskGroup, totalnum, ntime, R, Rd, Rdd, Rls1, Rls2, Rls3, Lls1, Lls2, Lls3,nthreads, debugging,KeepConstant);
    //
    if (verbose){
        end_point = system_clock::now();
        ending = time_point_cast<microseconds>(end_point).time_since_epoch().count();
        Rcout <<"df100 "<<(ending-start)<<" "<<halves<<" "<<iteration<<" "<<ind0<<",Update_Sides"<<endl;
        //
        gibtime = system_clock::to_time_t(system_clock::now());
        Rcout << ctime(&gibtime) << endl;
        Rcout << "Wrapping up" << endl;
    }
    //
    Calc_LogLik_Basic( nthreads, RiskFail, RiskGroup, totalnum, ntime, R, Rd, Rdd,RdR, Rls1, Rls2, Rls3, Lls1, Lls2, Lls3, Ll, Lld, Lldd, debugging, ties_method,KeepConstant);
    //
    List control_list = List::create(_["Iteration"]=iteration, _["Maximum Step"]=abs_max, _["Derivative Limiting"]=Lld_worst);//stores the total number of iterations used
    //
    int kept_covs = totalnum - sum(KeepConstant); //does not base the standard deviation off of constant parameters
    NumericVector Lldd_vec(kept_covs * kept_covs);
    #pragma omp parallel for schedule(dynamic) num_threads(nthreads)
    for (int ijk=0;ijk<totalnum*(totalnum+1)/2;ijk++){
        int ij = 0;
        int jk = ijk;
        int pij_ind=-100;
        int pjk_ind=-100;
        while (jk>ij){
            ij++;
            jk-=ij;
        }
        if (KeepConstant[ij]==0){
            pij_ind = ij - sum(head(KeepConstant,ij));
            if (KeepConstant[jk]==0){
                pjk_ind = jk - sum(head(KeepConstant,jk));
                Lldd_vec[pij_ind * kept_covs + pjk_ind]=Lldd[ij*totalnum+jk];
            }
        }
    }
    for (int ijk=0;ijk<kept_covs*(kept_covs+1)/2;ijk++){
        int ij = 0;
        int jk = ijk;
        while (jk>ij){
            ij++;
            jk-=ij;
        }
        Lldd_vec[jk * kept_covs + ij]=Lldd_vec[ij * kept_covs + jk];
    }
    Lldd_vec.attr("dim") = Dimension(kept_covs, kept_covs);
    const Map<MatrixXd> Lldd_mat(as<Map<MatrixXd> >(Lldd_vec));
    //
    MatrixXd Lldd_inv = -1 * Lldd_mat.inverse().matrix(); //uses inverse information matrix to calculate the standard deviation
    VectorXd stdev = VectorXd::Zero(totalnum);
    for (int ij=0;ij<totalnum;ij++){
        if (KeepConstant[ij]==0){
            int pij_ind = ij - sum(head(KeepConstant,ij));
            stdev(ij) = sqrt(Lldd_inv(pij_ind,pij_ind));
        }
    }
    //
    List res_list = List::create(_["LogLik"]=wrap(Ll[0]),_["First_Der"]=wrap(Lld),_["Second_Der"]=Lldd_vec,_["beta_0"]=wrap(beta_0) ,_["Standard_Deviation"]=wrap(stdev) ,_["AIC"]=2*(totalnum-accumulate(KeepConstant.begin(),KeepConstant.end(), 0.0))-2*Ll[0],_["Control_List"]=control_list,_["Converged"]=convgd);
    // returns a list of results
    return res_list;
}

//' Primary Cox PH regression with STRATA effect
//' \code{LogLik_Cox_PH_STRATA} Performs the calls to calculation functions, Structures the Cox PH regression, With verbose option prints out time stamps and intermediate sums of terms and derivatives
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
//' @param     double_step controls the step calculation, 0 for independent changes, 1 for solving b=Ax with complete matrices
//' @param     change_all    boolean if every parameter is being updated
//' @param     verbose    verbosity boolean
//' @param     debugging    debugging boolean
//' @param     KeepConstant    vector identifying constant parameters
//' @param     term_tot    total number of terms
//' @param     ties_method    ties method
//' @param     STRATA_vals vector of strata identifier values
//' @param     nthreads number of threads to use
//'
//' @return List of results: Log-likelihood of optimum, first derivative of log-likelihood, second derivative matrix, parameter list, standard deviation estimate, AIC, model information
// [[Rcpp::export]]
List LogLik_Cox_PH_STRATA( IntegerVector Term_n, StringVector tform, NumericVector a_n,NumericMatrix x_all,IntegerVector dfc,int fir, int der_iden,string modelform, double lr, int maxiter, int halfmax, double epsilon, double dbeta_cap, double abs_max,double dose_abs_max, double deriv_epsilon, NumericMatrix df_groups, NumericVector tu, int double_step ,bool change_all, bool verbose, bool debugging, IntegerVector KeepConstant, int term_tot, string ties_method, NumericVector& STRATA_vals, int nthreads){
    ;
    //
    List temp_list = List::create(_["Status"]="FAILED"); //used as a dummy return value for code checking
    if (verbose){
        Rcout << "START_COX_STRATA" << endl;
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
    int reqrdnum = totalnum - sum(KeepConstant);
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
    // int nthreads = Eigen::nbThreads()-1; //stores how many threads are allocated
    //
    //
    // Lld_worst: stores the highest magnitude log-likelihood derivative
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
    MatrixXd Td0 = MatrixXd::Zero(df0.rows(), reqrdnum); //preallocates matrix for Term derivative columns
    MatrixXd Tdd0 = MatrixXd::Zero(df0.rows(), reqrdnum*(reqrdnum+1)/2); //preallocates matrix for Term second derivative columns
    //
    MatrixXd Te = MatrixXd::Zero(df0.rows(), 1); //preallocates matrix for column terms used for temporary storage
    MatrixXd R = MatrixXd::Zero(df0.rows(), 1); //preallocates matrix for Risks
    ColXd Rd = MatrixXd::Zero(df0.rows(), reqrdnum); //preallocates matrix for Risk derivatives
    ColXd Rdd = MatrixXd::Zero(df0.rows(), reqrdnum*(reqrdnum+1)/2); //preallocates matrix for Risk second derivatives
    //
    MatrixXd Dose = MatrixXd::Constant(df0.rows(),term_tot,0.0); //Matrix of the total dose term values
    MatrixXd nonDose = MatrixXd::Constant(df0.rows(),term_tot,1.0); //Matrix of the total non-dose term values
    MatrixXd nonDose_LIN = MatrixXd::Constant(df0.rows(),term_tot,0.0); //matrix of Linear subterm values
    MatrixXd nonDose_PLIN = MatrixXd::Constant(df0.rows(),term_tot,1.0); //matrix of Loglinear subterm values
    MatrixXd nonDose_LOGLIN = MatrixXd::Constant(df0.rows(),term_tot,1.0); //matrix of Product linear subterm values
    MatrixXd TTerm = MatrixXd::Zero(Dose.rows(),Dose.cols()); //matrix of term values
    double dint = dose_abs_max; //The amount of change used to calculate derivatives in threshold paramters
    double dslp = abs_max;
    //
    if (verbose){
        Rcout << "starting subterms " << term_tot << endl;
    }
    // Calculates the subterm and term values
    Make_subterms( totalnum, Term_n, tform, dfc, fir, T0, Td0, Tdd0, Dose, nonDose, TTerm,  nonDose_LIN, nonDose_PLIN, nonDose_LOGLIN ,beta_0, df0,dint,dslp,nthreads, debugging,KeepConstant);
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
        for (int ijk=0;ijk<reqrdnum;ijk++){
            Rcout << Td0.col(ijk).sum() << " ";
        }
        Rcout << " " << endl;
        Rcout << "second derivs checked ";
        for (int ijk=0;ijk<reqrdnum;ijk++){
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
    ColXd RdR = MatrixXd::Zero(df0.rows(), reqrdnum); //preallocates matrix for Risk to derivative ratios
    ColXd RddR = MatrixXd::Zero(df0.rows(), reqrdnum*(reqrdnum+1)/2); //preallocates matrix for Risk to second derivative ratios
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
    Make_Risks(modelform, tform, Term_n, totalnum, fir, T0, Td0, Tdd0, Te, R, Rd, Rdd, Dose, nonDose, TTerm, nonDose_LIN, nonDose_PLIN, nonDose_LOGLIN, RdR, RddR, nthreads, debugging,KeepConstant);
    //
    // Removes infinite values
    RdR = (RdR.array().isFinite()).select(RdR,0);
    RddR = (RddR.array().isFinite()).select(RddR,0);
    //
    if (R.minCoeff()<0){
        Rcout << "A non-positive risk was detected: " << R.minCoeff() << endl;
        Rcout << "final failing values ";
        for (int ijk=0;ijk<totalnum;ijk++){
            Rcout << beta_0[ijk] << " ";
        }
        Rcout << " " << endl;
        //
        Rcout << "final failing terms ";
        for (int ijk=0;ijk<totalnum;ijk++){
            Rcout << tform[ijk] << " ";
        }
        Rcout << " " << endl;
        temp_list = List::create(_["beta_0"]=wrap(beta_0) ,_["Deviation"]=R_NaN,_["Status"]="FAILED",_["LogLik"]=R_NaN);
        return temp_list;
    }
    //
    if (verbose){
        Rcout << "risk checked ";
        for (int ijk=0;ijk<1;ijk++){
            Rcout << R.col(0).sum() << " ";
        }
        Rcout << " " << endl;
        Rcout << "risk1 checked ";
        for (int ijk=0;ijk<reqrdnum;ijk++){
            Rcout << Rd.col(ijk).sum() << " ";
        }
        Rcout << " " << endl;
        Rcout << "risk2 checked ";
        for (int ijk=0;ijk<reqrdnum;ijk++){
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
    StringMatrix RiskGroup(ntime,STRATA_vals.size()); //vector of strings detailing the rows
    IntegerMatrix RiskFail(ntime,2*STRATA_vals.size()); //vector giving the event rows
    //
    if (verbose){
        Rcout << "Grouping Start" << endl;
    }
    const Map<MatrixXd> df_m(as<Map<MatrixXd> >(df_groups));
    // Creates matrices used to identify the event risk groups
    Make_Groups_STRATA( ntime, df_m, RiskFail, RiskGroup, tu, nthreads, debugging,STRATA_vals);
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
    MatrixXd Rls1 =MatrixXd::Zero(ntime, STRATA_vals.size()); //precomputes a series of sums used frequently in the log-liklihood calculations
    MatrixXd Rls2 =MatrixXd::Zero(ntime, reqrdnum*STRATA_vals.size()); //Many are repeated due to the same risk groups and derivatives being used at mulitple points
    MatrixXd Rls3 =MatrixXd::Zero(ntime, reqrdnum*(reqrdnum+1)/2*STRATA_vals.size()); //Sum and its derivatives are precomputed
    MatrixXd Lls1 =MatrixXd::Zero(ntime, STRATA_vals.size()); //The log-likelihood calculation has a Right and Left sum used
    MatrixXd Lls2 =MatrixXd::Zero(ntime, reqrdnum*STRATA_vals.size());
    MatrixXd Lls3 =MatrixXd::Zero(ntime, reqrdnum*(reqrdnum+1)/2*STRATA_vals.size());
    //The log-likelihood is calculated in parallel over the risk groups
    vector<double> Ll( reqrdnum,0.0); //Log-likelihood values
    vector<double> Lld(reqrdnum,0.0); //Log-likelihood derivative values
    vector<double> Lldd(pow(reqrdnum,2),0.0);//The second derivative matrix has room for every combination, but only the lower triangle is calculated initially
    //
    // Calculates the side sum terms used
    Calculate_Sides_STRATA( RiskFail, RiskGroup, totalnum, ntime, R, Rd, Rdd, Rls1, Rls2, Rls3, Lls1, Lls2, Lls3,nthreads, debugging, STRATA_vals,KeepConstant);
    //
    // Calculates log-likelihood
    Calc_LogLik_STRATA( nthreads, RiskFail, RiskGroup, totalnum, ntime, R, Rd, Rdd,RdR,RddR, Rls1, Rls2, Rls3, Lls1, Lls2, Lls3, Ll, Lld, Lldd, debugging, ties_method, STRATA_vals,KeepConstant);
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
        for (int ij=0;ij<reqrdnum;ij++){
            Rcout << Ll[ij] << " ";
        }
        Rcout << " " << endl;
        Rcout << "df102 ";//prints the first derivatives
        for (int ij=0;ij<reqrdnum;ij++){
            Rcout << Lld[ij] << " ";
        }
        Rcout << " " << endl;
        Rcout << "df103 ";//prints the second derivatives
        for (int ij=0;ij<reqrdnum;ij++){
            Rcout << Lldd[ij*reqrdnum+ij] << " ";
        }
        Lld_worst=0;
        for (int ij=0;ij<reqrdnum;ij++){//locates highest magnitude derivative
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
        for (int ij=0;ij<reqrdnum;ij++){//prints the newton step value for zero derivative
            Rcout << Lld[ij]/Lldd[ij*reqrdnum+ij] << " ";
        }
        Rcout << " " << endl;
        Rcout << "df106 ";
        for (int ij=0;ij<reqrdnum;ij++){//prints the newton step value for zero log-likelihood
            Rcout << Ll[ij]/Lld[ij] << " ";
        }
        Rcout << " " << endl;
        Rcout << "df107 " << double_step << " " << abs_max << " " << dose_abs_max << " " << Ll_comp[0] << " " << Ll_comp[1] << endl;//prints several convergence terms
    }
    //
    vector<double> dbeta(totalnum,0.0);
    //
    // --------------------------
    // always starts from initial guess
    // --------------------------
    vector<double> beta_c(totalnum,0.0);
    vector<double> beta_a(totalnum,0.0);
    vector<double> beta_best(totalnum,0.0);
    vector<double> beta_p(totalnum,0.0);
    //
    //
    //int //risk_check_iter=0;
    VectorXd::Map(&beta_p[0], beta_0.size()) = beta_0;// stores previous parameters
    VectorXd::Map(&beta_c[0], beta_0.size()) = beta_0;// stores current parameters
    VectorXd::Map(&beta_a[0], beta_0.size()) = beta_0;// stores a refrence value for parameters
    VectorXd::Map(&beta_best[0], beta_0.size()) = beta_0;// stores the best parameters
    double halves = 0; //number of half-steps taken
    int ind0 = fir; //used for validations
    //int //i = ind0;
    int iteration=0; //iteration number
    //
    int iter_stop =0; //tracks if the iterations should be stopped for convergence
    int iter_check=0; //signal to check for convergence
    bool convgd = FALSE;
    //
    double Ll_abs_best = 10;
    vector<double> beta_abs_best(totalnum,0.0);
    while ((iteration < maxiter)&&(iter_stop==0)){
        iteration++;
        beta_p = beta_c;//
        beta_a = beta_c;//
        beta_best = beta_c;//
        Calc_Change( double_step, nthreads, totalnum, fir, der_iden, dbeta_cap, dose_abs_max, lr, abs_max, Ll, Lld, Lldd, dbeta, change_all, tform, dint,dslp, KeepConstant, debugging);
        Intercept_Bound(nthreads, totalnum, beta_0, dbeta, dfc, df0, KeepConstant, debugging, tform);
        if (verbose){
            Rcout << "Starting Halves"<<endl;//prints the final changes for validation
        }
        //
        //
        if ((Ll_abs_best>0)||(Ll_abs_best < Ll[ind0])){
            Ll_abs_best = Ll[ind0];
            beta_abs_best = beta_c;
        }
        //
        if (verbose){
            Rcout << "df501 " << Ll_abs_best << endl;
            Rcout << "df504 ";//prints parameter values
            for (int ij=0;ij<totalnum;ij++){
                Rcout << beta_abs_best[ij] << " ";
            }
            Rcout << " " << endl;
        }
        halves=0;
        while ((Ll[ind0] <= Ll_abs_best)&&(halves<halfmax)){ //repeats until half-steps maxed or an improvement
            //Refreshes the matrices used
            Dose = MatrixXd::Zero(df0.rows(),term_tot);
            nonDose = MatrixXd::Constant(df0.rows(),term_tot,0.0);
            nonDose_LIN = MatrixXd::Constant(df0.rows(),term_tot,0.0);
            nonDose_PLIN = MatrixXd::Constant(df0.rows(),term_tot,1.0);
            nonDose_LOGLIN = MatrixXd::Constant(df0.rows(),term_tot,1.0);
            T0 = MatrixXd::Zero(df0.rows(), totalnum); //preallocates matrix for Term column
            Td0 = MatrixXd::Zero(df0.rows(), reqrdnum); //preallocates matrix for Term derivative columns
            Tdd0 = MatrixXd::Zero(df0.rows(), reqrdnum*(reqrdnum+1)/2); //preallocates matrix for Term second derivative columns
            for (int ij=0;ij<totalnum;ij++){
                beta_0[ij] = beta_a[ij] + dbeta[ij];
                beta_c[ij] = beta_0[ij];
            }
            // ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------//
            // The same subterm, risk, sides, and log-likelihood calculations are performed every half-step and iteration
            // ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------//
            Make_subterms( totalnum, Term_n, tform, dfc, fir, T0, Td0, Tdd0, Dose,  nonDose, TTerm,  nonDose_LIN, nonDose_PLIN, nonDose_LOGLIN ,beta_0, df0,dint,dslp,nthreads, debugging,KeepConstant);
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
                for (int ijk=0;ijk<reqrdnum;ijk++){
                    Rcout << Td0.col(ijk).sum() << " ";
                }
                Rcout << " " << endl;
                Rcout << "second derivs checked ";
                for (int ijk=0;ijk<reqrdnum;ijk++){
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
            RdR = MatrixXd::Zero(df0.rows(), reqrdnum);
            RddR = MatrixXd::Zero(df0.rows(), reqrdnum*(reqrdnum+1)/2);
            //
            //
            Make_Risks(modelform, tform, Term_n, totalnum, fir, T0, Td0, Tdd0, Te, R, Rd, Rdd, Dose, nonDose, TTerm, nonDose_LIN, nonDose_PLIN, nonDose_LOGLIN, RdR, RddR, nthreads, debugging,KeepConstant);
            //
            if (R.minCoeff()<0){
                #pragma omp parallel for num_threads(nthreads)
                for (int ijk=0;ijk<totalnum;ijk++){
                    int tij = Term_n[ijk];
                    if (TTerm.col(tij).minCoeff()<=0){
                        dbeta[ijk] = dbeta[ijk] / 2.0;
                    } else if (isinf(TTerm.col(tij).maxCoeff())){
                        dbeta[ijk] = dbeta[ijk] / 2.0;
                    }
                }
                if (verbose){
                    Rcout << "A non-positive risk was detected: " << R.minCoeff() << endl;
                }
                halves+=0.2;
            } else {
                halves++;
                RdR = (RdR.array().isFinite()).select(RdR,0);
                RddR = (RddR.array().isFinite()).select(RddR,0);
                //
                if (verbose){
                    Rcout << "risk checked ";
                    for (int ijk=0;ijk<1;ijk++){
                        Rcout << R.col(0).sum() << " ";
                    }
                    Rcout << " " << endl;
                    Rcout << "risk1 checked ";
                    for (int ijk=0;ijk<reqrdnum;ijk++){
                        Rcout << Rd.col(ijk).sum() << " ";
                    }
                    Rcout << " " << endl;
                    Rcout << "risk2 checked ";
                    for (int ijk=0;ijk<reqrdnum;ijk++){
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
                fill(Ll.begin(),   Ll.end(), 0.0);
                fill(Lld.begin(),  Lld.end(), 0.0);
                fill(Lldd.begin(), Lldd.end(), 0.0);
                Rls1 =MatrixXd::Zero(ntime, STRATA_vals.size()); //precomputes a series of sums used frequently in the log-liklihood calculations
                Rls2 =MatrixXd::Zero(ntime, reqrdnum*STRATA_vals.size()); //Many are repeated due to the same risk groups and derivatives being used at mulitple points
                Rls3 =MatrixXd::Zero(ntime, reqrdnum*(reqrdnum+1)/2*STRATA_vals.size());
                Lls1 =MatrixXd::Zero(ntime, STRATA_vals.size());
                Lls2 =MatrixXd::Zero(ntime, reqrdnum*STRATA_vals.size());
                Lls3 =MatrixXd::Zero(ntime, reqrdnum*(reqrdnum+1)/2*STRATA_vals.size());
                //
                Calculate_Sides_STRATA( RiskFail, RiskGroup, totalnum, ntime, R, Rd, Rdd, Rls1, Rls2, Rls3, Lls1, Lls2, Lls3,nthreads, debugging, STRATA_vals,KeepConstant);
                //
                Calc_LogLik_STRATA( nthreads, RiskFail, RiskGroup, totalnum, ntime, R, Rd, Rdd,RdR,RddR, Rls1, Rls2, Rls3, Lls1, Lls2, Lls3, Ll, Lld, Lldd, debugging, ties_method, STRATA_vals,KeepConstant);
                if (verbose){
                    end_point = system_clock::now();
                    ending = time_point_cast<microseconds>(end_point).time_since_epoch().count();
                    Rcout <<"df100 "<<(ending-start)<<" "<<0<<" "<<0<<" "<<0<<",Half Calc"<<endl;//prints the time
                    gibtime = system_clock::to_time_t(system_clock::now());
                    Rcout << ctime(&gibtime) << endl;
                    Rcout << "df101 ";//prints the log-likelihoods
                    for (int ij=0;ij<reqrdnum;ij++){
                        Rcout << Ll[ij] << " ";
                    }
                    Rcout << " " << endl;
                    Rcout << "df102 ";//prints the first derivatives
                    for (int ij=0;ij<reqrdnum;ij++){
                        Rcout << Lld[ij] << " ";
                    }
                    Rcout << " " << endl;
                    Rcout << "df103 ";//prints the second derivatives
                    for (int ij=0;ij<reqrdnum;ij++){
                        Rcout << Lldd[ij*reqrdnum+ij] << " ";
                    }
                    Lld_worst=0;
                    for (int ij=0;ij<reqrdnum;ij++){//locates highest magnitude derivative
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
                    for (int ij=0;ij<reqrdnum;ij++){//prints the newton step value for zero derivative
                        Rcout << Lld[ij]/Lldd[ij*reqrdnum+ij] << " ";
                    }
                    Rcout << " " << endl;
                    Rcout << "df106 ";
                    for (int ij=0;ij<reqrdnum;ij++){//prints the newton step value for zero log-likelihood
                        Rcout << Ll[ij]/Lld[ij] << " ";
                    }
                    Rcout << " " << endl;
                    Rcout << "df107 " << double_step << " " << abs_max << " " << dose_abs_max << " " << Ll_comp[0] << " " << Ll_comp[1] << endl;//prints several convergence terms
                }
                if (change_all){ //If every covariate is to be changed
                    if (Ll[ind0] <= Ll_abs_best){//takes a half-step if needed
                        #pragma omp parallel for num_threads(nthreads)
                        for (int ijk=0;ijk<totalnum;ijk++){
                            dbeta[ijk] = dbeta[ijk] * 0.5; //
                        }
                    } else{//If improved, updates the best vector
                        if (verbose){
                            Rcout << "Optimal" << endl;
                        }
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
                }
                #pragma omp parallel for num_threads(nthreads)
                for (int ijk=0;ijk<totalnum;ijk++){//totalnum*(totalnum+1)/2
                    beta_0[ijk] = beta_c[ijk];
                }
            }
        }
        if (beta_best!=beta_c){//if the risk matrices aren't the optimal values, then they must be recalculated
            if (verbose){
                Rcout << "Changing back to best"<<endl;
            }
            // If it goes through every half step without improvement, then the maximum change needs to be decreased
            abs_max = abs_max*pow(0.5,halfmax); // reduces the step sizes
            dose_abs_max = dose_abs_max*pow(0.5,halfmax);
            iter_check = 1;
            //
            beta_p = beta_best;//
            beta_a = beta_best;//
            beta_c = beta_best;//
            Dose = MatrixXd::Zero(df0.rows(),term_tot);
            nonDose = MatrixXd::Constant(df0.rows(),term_tot,0.0);
            nonDose_LIN = MatrixXd::Constant(df0.rows(),term_tot,0.0);
            nonDose_PLIN = MatrixXd::Constant(df0.rows(),term_tot,1.0);
            nonDose_LOGLIN = MatrixXd::Constant(df0.rows(),term_tot,1.0);
            for (int ij=0;ij<totalnum;ij++){
                beta_0[ij] = beta_best[ij];
            }
            Make_subterms( totalnum, Term_n, tform, dfc, fir, T0, Td0, Tdd0, Dose, nonDose, TTerm,  nonDose_LIN, nonDose_PLIN, nonDose_LOGLIN,beta_0, df0,dint,dslp,nthreads, debugging,KeepConstant);
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
                for (int ijk=0;ijk<reqrdnum;ijk++){
                    Rcout << Td0.col(ijk).sum() << " ";
                }
                Rcout << " " << endl;
                Rcout << "second derivs checked ";
                for (int ijk=0;ijk<reqrdnum;ijk++){
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
            RdR = MatrixXd::Zero(df0.rows(), reqrdnum);
            RddR = MatrixXd::Zero(df0.rows(), reqrdnum*(reqrdnum+1)/2);
            //
            //
            Make_Risks(modelform, tform, Term_n, totalnum, fir, T0, Td0, Tdd0, Te, R, Rd, Rdd, Dose, nonDose, TTerm, nonDose_LIN, nonDose_PLIN, nonDose_LOGLIN, RdR, RddR, nthreads, debugging,KeepConstant);
            R = (R.array().isFinite()).select(R,0);
            Rd = (Rd.array().isFinite()).select(Rd,0);
            Rdd = (Rdd.array().isFinite()).select(Rdd,0);
            RdR = (RdR.array().isFinite()).select(RdR,0);
            RddR = (RddR.array().isFinite()).select(RddR,0);
            //
            temp_list = List::create(_["betas"]=wrap(beta_0),_["Status"]="FAILED"); //used as a dummy return value for code checking
            //risk_check_iter=0;
            if (verbose){
                Rcout << "risk checked ";
                for (int ijk=0;ijk<1;ijk++){
                    Rcout << R.col(0).sum() << " ";
                }
                Rcout << " " << endl;
                Rcout << "risk1 checked ";
                for (int ijk=0;ijk<reqrdnum;ijk++){
                    Rcout << Rd.col(ijk).sum() << " ";
                }
                Rcout << " " << endl;
                Rcout << "risk2 checked ";
                for (int ijk=0;ijk<reqrdnum;ijk++){
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
            Rls1 =MatrixXd::Zero(ntime, STRATA_vals.size()); //precomputes a series of sums used frequently in the log-liklihood calculations
            Rls2 =MatrixXd::Zero(ntime, reqrdnum*STRATA_vals.size()); //Many are repeated due to the same risk groups and derivatives being used at mulitple points
            Rls3 =MatrixXd::Zero(ntime, reqrdnum*(reqrdnum+1)/2*STRATA_vals.size());
            Lls1 =MatrixXd::Zero(ntime, STRATA_vals.size());
            Lls2 =MatrixXd::Zero(ntime, reqrdnum*STRATA_vals.size());
            Lls3 =MatrixXd::Zero(ntime, reqrdnum*(reqrdnum+1)/2*STRATA_vals.size());
            Calculate_Sides_STRATA( RiskFail, RiskGroup, totalnum, ntime, R, Rd, Rdd, Rls1, Rls2, Rls3, Lls1, Lls2, Lls3,nthreads, debugging, STRATA_vals,KeepConstant);
            //
            //
            Calc_LogLik_STRATA( nthreads, RiskFail, RiskGroup, totalnum, ntime, R, Rd, Rdd,RdR,RddR, Rls1, Rls2, Rls3, Lls1, Lls2, Lls3, Ll, Lld, Lldd, debugging, ties_method, STRATA_vals,KeepConstant);
        }
        Lld_worst=0;
        for (int ij=0;ij<reqrdnum;ij++){
            if (abs(Lld[ij]) > Lld_worst){
                Lld_worst = abs(Lld[ij]);
            }
        }
        if (iteration > reqrdnum){//Doesn't check the first several iterations for convergence
            if ((iteration % (reqrdnum))||(iter_check==1)){//Checks every set number of iterations
                iter_check=0;
                if (Lld_worst < deriv_epsilon){//ends if the derivatives are low enough
                    iter_stop = 1;
                    convgd = TRUE;
                }
                Ll_comp[1]=Ll[0];
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
            for (int ij=0;ij<reqrdnum;ij++){
                Rcout << Ll[ij] << " ";
            }
            Rcout << " " << endl;
            Rcout << "df102 ";
            for (int ij=0;ij<reqrdnum;ij++){
                Rcout << Lld[ij] << " ";
            }
            Rcout << " " << endl;
            Rcout << "df103 ";
            for (int ij=0;ij<reqrdnum;ij++){
                Rcout << Lldd[ij*reqrdnum+ij] << " ";
            }
            Rcout << " " << endl;
            Rcout << "df104 ";
            for (int ij=0;ij<totalnum;ij++){
                Rcout << beta_c[ij] << " ";
            }
            Rcout << " " << endl;
            Rcout << "df105 ";
            for (int ij=0;ij<reqrdnum;ij++){//prints the newton step value for zero derivative
                Rcout << Lld[ij]/Lldd[ij*reqrdnum+ij] << " ";
            }
            Rcout << " " << endl;
            Rcout << "df106 ";
            for (int ij=0;ij<reqrdnum;ij++){//prints the newton step value for zero log-likelihood
                Rcout << Ll[ij]/Lld[ij] << " ";
            }
            Rcout << " " << endl;
            Rcout << "df107 " << double_step << " " << abs_max << " " << dose_abs_max << " " << Ll_comp[0] << " " << Ll_comp[1] << endl;//prints several convergence terms
            Rcout << "Finshed iteration" << endl;
        }
    }
    // -----------------------------------------------
    // Performing Full Calculation to get full second derivative matrix
    // -----------------------------------------------
    fill(Ll.begin(), Ll.end(), 0.0);
    fill(Lld.begin(), Lld.end(), 0.0);
    fill(Lldd.begin(), Lldd.end(), 0.0);
    Calculate_Sides_STRATA( RiskFail, RiskGroup, totalnum, ntime, R, Rd, Rdd, Rls1, Rls2, Rls3, Lls1, Lls2, Lls3,nthreads, debugging, STRATA_vals,KeepConstant);
    if (verbose){
        Rcout << "Wrapping up" << endl;
    }
    //
    Calc_LogLik_STRATA( nthreads, RiskFail, RiskGroup, totalnum, ntime, R, Rd, Rdd,RdR,RddR, Rls1, Rls2, Rls3, Lls1, Lls2, Lls3, Ll, Lld, Lldd, debugging, ties_method, STRATA_vals,KeepConstant);
    //
    if ((Ll_abs_best>0)||(Ll_abs_best < Ll[ind0])){
        Ll_abs_best = Ll[ind0];
        beta_abs_best = beta_c;
    }
    //
    if (verbose){
        Rcout << "df501 " << Ll_abs_best << endl;
        Rcout << "df504 ";//prints parameter values
        for (int ij=0;ij<totalnum;ij++){
            Rcout << beta_abs_best[ij] << " ";
        }
        Rcout << " " << endl;
    }
    //
    List para_list = List::create(_["Term_n"]=Term_n,_["tforms"]=tform); //stores the term information
    List control_list = List::create(_["Iteration"]=iteration, _["Maximum Step"]=abs_max, _["Derivative Limiting"]=Lld_worst); //stores the total number of iterations used
    //
    int kept_covs = totalnum - sum(KeepConstant);
    NumericVector Lldd_vec(kept_covs * kept_covs);
    #pragma omp parallel for schedule(dynamic) num_threads(nthreads)
    for (int ijk=0;ijk<totalnum*(totalnum+1)/2;ijk++){
        int ij = 0;
        int jk = ijk;
        int pij_ind=-100;
        int pjk_ind=-100;
        while (jk>ij){
            ij++;
            jk-=ij;
        }
        if (KeepConstant[ij]==0){
            pij_ind = ij - sum(head(KeepConstant,ij));
            if (KeepConstant[jk]==0){
                pjk_ind = jk - sum(head(KeepConstant,jk));
                Lldd_vec[pij_ind * kept_covs + pjk_ind]=Lldd[ij*totalnum+jk];
            }
        }
    }
    for (int ijk=0;ijk<kept_covs*(kept_covs+1)/2;ijk++){
        int ij = 0;
        int jk = ijk;
        while (jk>ij){
            ij++;
            jk-=ij;
        }
        Lldd_vec[jk * kept_covs + ij]=Lldd_vec[ij * kept_covs + jk];
    }
    Lldd_vec.attr("dim") = Dimension(kept_covs, kept_covs);
    const Map<MatrixXd> Lldd_mat(as<Map<MatrixXd> >(Lldd_vec));
    //
    MatrixXd Lldd_inv = -1 * Lldd_mat.inverse().matrix(); //uses inverse information matrix to calculate the standard deviation
    VectorXd stdev = VectorXd::Zero(totalnum);
    for (int ij=0;ij<totalnum;ij++){
        if (KeepConstant[ij]==0){
            int pij_ind = ij - sum(head(KeepConstant,ij));
            stdev(ij) = sqrt(Lldd_inv(pij_ind,pij_ind));
        }
    }
    //
    List res_list = List::create(_["LogLik"]=wrap(Ll[0]),_["First_Der"]=wrap(Lld),_["Second_Der"]=Lldd_vec,_["beta_0"]=wrap(beta_0) ,_["Standard_Deviation"]=wrap(stdev) ,_["AIC"]=2*(totalnum-accumulate(KeepConstant.begin(),KeepConstant.end(), 0.0))-2*Ll[0],_["Parameter_Lists"]=para_list,_["Control_List"]=control_list,_["Converged"]=convgd);
    // returns a list of results
    return res_list;
}


//' Primary Cox PH baseline hazard function
//' \code{Cox_PH_PLOT_SURV} Performs the calls to calculation functions, Uses calculated risks and risk groups to approximate the baseline, With verbose option prints out time stamps and intermediate sums of terms and derivatives
//'
//' @param     Term_n    Term numbers
//' @param     tform    subterm types
//' @param     a_n    Optimal values
//' @param     a_er    Optimal value standard error
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
//' @param     nthreads number of threads to use
//'
//' @return List of results: baseline hazard, risk for each row
// [[Rcpp::export]]
List Cox_PH_PLOT_SURV(IntegerVector Term_n, StringVector tform, NumericVector a_n, NumericVector a_er,NumericMatrix x_all,IntegerVector dfc,int fir, int der_iden,string modelform, double abs_max,double dose_abs_max, NumericMatrix df_groups, NumericVector tu , bool verbose, bool debugging, IntegerVector KeepConstant, int term_tot, int nthreads){
    //
    // Calculates the baseline hazard
    //
    using namespace std::chrono;
    if (verbose){
        Rcout << "START_COX_SURV" << endl;
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
    int reqrdnum = totalnum - sum(KeepConstant);
    //
    //
    // cout.precision: controls the number of significant digits printed
    // nthreads: number of threads used for parallel operations
    //
    Rcout.precision(10); //forces higher precision numbers printed to terminal
    // int nthreads = Eigen::nbThreads()-1; //stores how many threads are allocated
    //
    //
    //
    // Lld_worst: stores the highest magnitude log-likelihood derivative
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
    MatrixXd Td0 = MatrixXd::Zero(df0.rows(), reqrdnum); //preallocates matrix for Term derivative columns
    MatrixXd Tdd0 = MatrixXd::Zero(df0.rows(), reqrdnum*(reqrdnum+1)/2); //preallocates matrix for Term second derivative columns
    //
    MatrixXd Te = MatrixXd::Zero(df0.rows(), 1); //preallocates matrix for column terms used for temporary storage
    MatrixXd R = MatrixXd::Zero(df0.rows(), 1); //preallocates matrix for Risks
    ColXd Rd = MatrixXd::Zero(df0.rows(), reqrdnum); //preallocates matrix for Risk derivatives
    ColXd Rdd = MatrixXd::Zero(df0.rows(), reqrdnum*(reqrdnum+1)/2); //preallocates matrix for Risk second derivatives
    //
    MatrixXd Dose = MatrixXd::Zero(df0.rows(),term_tot); //Matrix of the total dose term values
    MatrixXd nonDose = MatrixXd::Constant(df0.rows(),term_tot,1.0); //Matrix of the total non-dose term values
    MatrixXd nonDose_LIN = MatrixXd::Constant(df0.rows(),term_tot,0.0); //matrix of Linear subterm values
    MatrixXd nonDose_PLIN = MatrixXd::Constant(df0.rows(),term_tot,1.0); //matrix of Loglinear subterm values
    MatrixXd nonDose_LOGLIN = MatrixXd::Constant(df0.rows(),term_tot,1.0); //matrix of Product linear subterm values
    MatrixXd TTerm = MatrixXd::Zero(Dose.rows(),Dose.cols()); //matrix of term values
    double dint = dose_abs_max; //The amount of change used to calculate derivatives in threshold paramters    //
    double dslp = abs_max;
    // totalnum,& Term_n,  tform, dfc,& fir,& T0,& Td0,& Tdd0,& Dose,& nonDose,& beta_0,& df0, dint, nthreads,  debugging
    if (verbose){
        Rcout << "starting subterms " << term_tot << endl;
    }
    Make_subterms( totalnum, Term_n, tform, dfc, fir, T0, Td0, Tdd0, Dose, nonDose, TTerm,  nonDose_LIN, nonDose_PLIN, nonDose_LOGLIN ,beta_0, df0,dint,dslp,nthreads, debugging,KeepConstant);
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
        for (int ijk=0;ijk<reqrdnum;ijk++){
            Rcout << Td0.col(ijk).sum() << " ";
        }
        Rcout << " " << endl;
        Rcout << "second derivs checked ";
        for (int ijk=0;ijk<reqrdnum;ijk++){
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
    ColXd RdR = MatrixXd::Zero(df0.rows(), reqrdnum); //preallocates matrix for Risk to derivative ratios
    ColXd RddR = MatrixXd::Zero(df0.rows(), reqrdnum*(reqrdnum+1)/2); //preallocates matrix for Risk to second derivative ratios
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
    Make_Risks(modelform, tform, Term_n, totalnum, fir, T0, Td0, Tdd0, Te, R, Rd, Rdd, Dose, nonDose, TTerm, nonDose_LIN, nonDose_PLIN, nonDose_LOGLIN, RdR, RddR, nthreads, debugging,KeepConstant);
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
    vector<double> hazard_error(ntime,0.0);
    #pragma omp parallel for schedule(dynamic) num_threads(nthreads)
    for (int ijk=0;ijk<reqrdnum;ijk++){
        Rd.col(ijk) = Rd.col(ijk).array().pow(2).array() * pow(a_er[ijk],2);
    }
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
        double Rds1 = 0; //total weighted risk derivative squared
        for (vector<double>::size_type i = 0; i < indices.size()-1; i=i+2){
            Rs1 += R.block(indices[i]-1,0,indices[i+1]-indices[i]+1,1).sum();
            Rds1 += Rd.block(indices[i]-1,0,indices[i+1]-indices[i]+1,reqrdnum).sum();
        }
        baseline[ijk] = dj / Rs1; //approximates the baseline hazard
        hazard_error[ijk] = dj / pow(Rs1,2);
        //
    }
    //
    if (verbose){
        end_point = system_clock::now();
        ending = time_point_cast<microseconds>(end_point).time_since_epoch().count();
        Rcout <<"df100 "<<(ending-start)<<" " << R.rows() <<" Finishing Baseline"<<endl;
        gibtime = system_clock::to_time_t(system_clock::now());
        Rcout << ctime(&gibtime) << endl;
    }
    NumericVector w_base = wrap(baseline);
    NumericVector w_base_er = wrap(hazard_error);
    NumericVector w_R = wrap(R.col(0));
    // returns the baseline approximates and the risk information
    List res_list = List::create(_["baseline"]=w_base, _["standard_error"]=w_base_er, _["Risks"]=w_R);
    //
    if (verbose){
        Rcout << "returning" << endl;
    }
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
//' @param     uniq_v    number of unique covariate values
//' @param     nthreads number of threads to use
//'
//' @return List of results: covariate values, risks for each row
// [[Rcpp::export]]
List Cox_PH_PLOT_RISK(IntegerVector Term_n, StringVector tform, NumericVector a_n,NumericMatrix x_all,IntegerVector dfc,int fir, int der_iden,string modelform, double abs_max,double dose_abs_max, NumericMatrix df_groups, NumericVector tu, bool verbose, bool debugging, IntegerVector KeepConstant, int term_tot, int uniq_v, int nthreads){
    ;
    //
    // Plots the risk over a series of covariate values
    //
    using namespace std::chrono;
    if (verbose){
        Rcout << "START_COX_RISK" << endl;
    }
    //
    List res_temp = List::create(_["eh"]=wrap(Term_n));
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
    int reqrdnum = totalnum - sum(KeepConstant);
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
    // int nthreads = Eigen::nbThreads()-1; //stores how many threads are allocated
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
    df0 = df0.array();
    //
    Map<VectorXd> beta_0(as<Map<VectorXd> >(a_n));
    MatrixXd T0 = MatrixXd::Zero(df0.rows(), totalnum); //preallocates matrix for Term column
    MatrixXd Td0 = MatrixXd::Zero(df0.rows(), reqrdnum); //preallocates matrix for Term derivative columns
    MatrixXd Tdd0 = MatrixXd::Zero(df0.rows(), reqrdnum*(reqrdnum+1)/2); //preallocates matrix for Term second derivative columns
    //
    MatrixXd Te = MatrixXd::Zero(df0.rows(), 1); //preallocates matrix for column terms used for temporary storage
    MatrixXd R = MatrixXd::Zero(df0.rows(), 1); //preallocates matrix for Risks
    ColXd Rd = MatrixXd::Zero(df0.rows(), reqrdnum); //preallocates matrix for Risk derivatives
    ColXd Rdd = MatrixXd::Zero(df0.rows(), reqrdnum*(reqrdnum+1)/2); //preallocates matrix for Risk second derivatives
    //
    MatrixXd Dose = MatrixXd::Constant(df0.rows(),term_tot,0.0); //Matrix of the total dose term values
    MatrixXd nonDose = MatrixXd::Constant(df0.rows(),term_tot,1.0); //Matrix of the total non-dose term values
    MatrixXd nonDose_LIN = MatrixXd::Constant(df0.rows(),term_tot,0.0); //matrix of Linear subterm values
    MatrixXd nonDose_PLIN = MatrixXd::Constant(df0.rows(),term_tot,1.0); //matrix of Loglinear subterm values
    MatrixXd nonDose_LOGLIN = MatrixXd::Constant(df0.rows(),term_tot,1.0); //matrix of Product linear subterm values
    MatrixXd TTerm = MatrixXd::Zero(Dose.rows(),Dose.cols()); //matrix of term values
    double dint = dose_abs_max; //The amount of change used to calculate derivatives in threshold paramters
    double dslp = abs_max;
    //
    if (verbose){
        Rcout << "Matrices made "<< endl;
    }
    int ijk= dfc[der_iden]-1;
    dx = (df1.col(ijk).maxCoeff() - df1.col(ijk).minCoeff())/(vv.size()-1);//varies from max to minimum
    vv[0] = df1.col(ijk).minCoeff();
    generate(vv.begin(), vv.end(), [n = 0, &dx]() mutable { return n++ * dx; });
    #pragma omp parallel for schedule(dynamic) num_threads(nthreads)
    for (vector<float>::size_type ij=0;ij<vv.size();ij++){
        df0(ij,ijk)=vv[ij]; //fills the column with varying values
    }
    //
    if (verbose){
        Rcout << "Making subterm "<< endl;
    }
    // Calculates terms
    //
    //
    Make_subterms( totalnum, Term_n, tform, dfc, fir, T0, Td0, Tdd0, Dose, nonDose, TTerm,  nonDose_LIN, nonDose_PLIN, nonDose_LOGLIN ,beta_0, df0,dint,dslp,nthreads, debugging,KeepConstant);
    ColXd RdR = MatrixXd::Zero(df0.rows(), reqrdnum); //preallocates matrix for Risk to derivative ratios
    ColXd RddR = MatrixXd::Zero(df0.rows(), reqrdnum*(reqrdnum+1)/2);
    //
    if (verbose){
        Rcout << "Making Risk"<< endl;
    }
    // Calculates risk
    Make_Risks(modelform, tform, Term_n, totalnum, fir, T0, Td0, Tdd0, Te, R, Rd, Rdd, Dose, nonDose, TTerm, nonDose_LIN, nonDose_PLIN, nonDose_LOGLIN, RdR, RddR, nthreads, debugging,KeepConstant);
    //
    // Returns the values used and risks calculated
    List res_list = List::create(_["x"]=wrap(df0.col(ijk)), _["y"]=wrap(R.col(0)));//returns list of covariate values and risk
    if (verbose){
        Rcout << "Returning"<< endl;
    }
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
//' @param     nthreads number of threads to use
//'
//' @return List of results: scaled schoenfeld residuals
// [[Rcpp::export]]
List Schoenfeld_Cox_PH( IntegerVector Term_n, StringVector tform, NumericVector a_n,NumericMatrix x_all,IntegerVector dfc,int fir, int der_iden,string modelform, double abs_max,double dose_abs_max, NumericMatrix df_groups, NumericVector tu , bool verbose, bool debugging, IntegerVector KeepConstant, int term_tot, string ties_method, int nthreads){
    ;
    //
    // Calculates the schoenfeld residuals
    //
    List temp_list = List::create(_["Status"]="FAILED"); //used as a dummy return value for code checking
    if (verbose){
        Rcout << "START_SCHOENFELD" << endl;
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
    int reqrdnum = totalnum - sum(KeepConstant);
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
    // int nthreads = Eigen::nbThreads()-1; //stores how many threads are allocated
    //
    //
    //
    // Lld_worst: stores the highest magnitude log-likelihood derivative
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
    MatrixXd Td0 = MatrixXd::Zero(df0.rows(), reqrdnum); //preallocates matrix for Term derivative columns
    MatrixXd Tdd0 = MatrixXd::Zero(df0.rows(), reqrdnum*(reqrdnum+1)/2); //preallocates matrix for Term second derivative columns
    //
    MatrixXd Te = MatrixXd::Zero(df0.rows(), 1); //preallocates matrix for column terms used for temporary storage
    MatrixXd R = MatrixXd::Zero(df0.rows(), 1); //preallocates matrix for Risks
    ColXd Rd = MatrixXd::Zero(df0.rows(), reqrdnum); //preallocates matrix for Risk derivatives
    ColXd Rdd = MatrixXd::Zero(df0.rows(), reqrdnum*(reqrdnum+1)/2); //preallocates matrix for Risk second derivatives
    //
    MatrixXd Dose = MatrixXd::Zero(df0.rows(),term_tot); //Matrix of the total dose term values
    MatrixXd nonDose = MatrixXd::Constant(df0.rows(),term_tot,1.0); //Matrix of the total non-dose term values
    MatrixXd nonDose_LIN = MatrixXd::Constant(df0.rows(),term_tot,0.0); //matrix of Linear subterm values
    MatrixXd nonDose_PLIN = MatrixXd::Constant(df0.rows(),term_tot,1.0); //matrix of Loglinear subterm values
    MatrixXd nonDose_LOGLIN = MatrixXd::Constant(df0.rows(),term_tot,1.0); //matrix of Product linear subterm values
    MatrixXd TTerm = MatrixXd::Zero(Dose.rows(),Dose.cols()); //matrix of term values
    double dint = dose_abs_max; //The amount of change used to calculate derivatives in threshold paramters    //
    double dslp = abs_max;
    // Calculates the subterm and term values
    if (verbose){
        Rcout << "starting subterms " << term_tot << endl;
    }
    Make_subterms( totalnum, Term_n, tform, dfc, fir, T0, Td0, Tdd0, Dose, nonDose, TTerm,  nonDose_LIN, nonDose_PLIN, nonDose_LOGLIN ,beta_0, df0,dint,dslp,nthreads, debugging,KeepConstant);
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
        for (int ijk=0;ijk<reqrdnum;ijk++){
            Rcout << Td0.col(ijk).sum() << " ";
        }
        Rcout << " " << endl;
        Rcout << "second derivs checked ";
        for (int ijk=0;ijk<reqrdnum;ijk++){
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
    // ---------------------------------------------------------
    // Prints off a series of calculations to check at what point values are changing
    // ---------------------------------------------------------
    //
    ColXd RdR = MatrixXd::Zero(df0.rows(), reqrdnum); //preallocates matrix for Risk to derivative ratios
    ColXd RddR = MatrixXd::Zero(df0.rows(), reqrdnum*(reqrdnum+1)/2);
    //
    Make_Risks(modelform, tform, Term_n, totalnum, fir, T0, Td0, Tdd0, Te, R, Rd, Rdd, Dose, nonDose, TTerm, nonDose_LIN, nonDose_PLIN, nonDose_LOGLIN, RdR, RddR, nthreads, debugging,KeepConstant);
    //
    if (verbose){
        Rcout << "risk checked ";
        for (int ijk=0;ijk<1;ijk++){
            Rcout << R.col(0).sum() << " ";
        }
        Rcout << " " << endl;
        Rcout << "risk1 checked ";
        for (int ijk=0;ijk<reqrdnum;ijk++){
            Rcout << Rd.col(ijk).sum() << " ";
        }
        Rcout << " " << endl;
        Rcout << "risk2 checked ";
        for (int ijk=0;ijk<reqrdnum;ijk++){
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
    MatrixXd Rls2 =MatrixXd::Zero(ntime, reqrdnum); //Many are repeated due to the same risk groups and derivatives being used at mulitple points
    MatrixXd Rls3 =MatrixXd::Zero(ntime, reqrdnum*(reqrdnum+1)/2); //Sum and its derivatives are precomputed
    MatrixXd Lls1 =MatrixXd::Zero(ntime, 1); //The log-likelihood calculation has a Right and Left sum used
    MatrixXd Lls2 =MatrixXd::Zero(ntime, reqrdnum);
    MatrixXd Lls3 =MatrixXd::Zero(ntime, reqrdnum*(reqrdnum+1)/2);
    //The log-likelihood is calculated in parallel over the risk groups
    vector<double> Ll(reqrdnum,0.0); //Log-likelihood values
    vector<double> Lld(reqrdnum,0.0); //Log-likelihood derivative values
    vector<double> Lldd(pow(reqrdnum,2),0.0);//The second derivative matrix has room for every combination, but only the lower triangle is calculated initially
    //
    // Calculates the side sum terms used
    Calculate_Sides( RiskFail, RiskGroup, totalnum, ntime, R, Rd, Rdd, Rls1, Rls2, Rls3, Lls1, Lls2, Lls3,nthreads, debugging,KeepConstant);
    //
    if (verbose){
        Rcout << "riskr checked ";
        for (int ijk=0;ijk<1;ijk++){
            Rcout << Rls1.col(0).sum() << " ";
        }
        Rcout << " " << endl;
        Rcout << "risk1r checked ";
        for (int ijk=0;ijk<reqrdnum;ijk++){
            Rcout << Rls2.col(ijk).sum() << " ";
        }
        Rcout << " " << endl;
        Rcout << "risk2r checked ";
        for (int ijk=0;ijk<reqrdnum;ijk++){
            Rcout << Rls3.col(ijk*(ijk+1)/2+ijk).sum() << " ";
        }
        Rcout << " " << endl;
        //
        Rcout << "riskl checked ";
        for (int ijk=0;ijk<1;ijk++){
            Rcout << Lls1.col(0).sum() << " ";
        }
        Rcout << " " << endl;
        Rcout << "risk1l checked ";
        for (int ijk=0;ijk<reqrdnum;ijk++){
            Rcout << Lls2.col(ijk).sum() << " ";
        }
        Rcout << " " << endl;
        Rcout << "risk2l checked ";
        for (int ijk=0;ijk<reqrdnum;ijk++){
            Rcout << Lls3.col(ijk*(ijk+1)/2+ijk).sum() << " ";
        }
        Rcout << " " << endl;
    }
    // Calculates log-likelihood
    Calc_LogLik( nthreads, RiskFail, RiskGroup, totalnum, ntime, R, Rd, Rdd,RdR,RddR, Rls1, Rls2, Rls3, Lls1, Lls2, Lls3, Ll, Lld, Lldd, debugging, ties_method,KeepConstant);
    //
    if (verbose){
        end_point = system_clock::now();
        ending = time_point_cast<microseconds>(end_point).time_since_epoch().count();
        Rcout <<"df100 "<<(ending-start)<<" "<<0<<" "<<0<<" "<<0<<",Calc"<<endl;//prints the time
        gibtime = system_clock::to_time_t(system_clock::now());
        Rcout << ctime(&gibtime) << endl;
        Rcout << "df101 ";//prints the log-likelihoods
        for (int ij=0;ij<reqrdnum;ij++){
            Rcout << Ll[ij] << " ";
        }
        Rcout << " " << endl;
        Rcout << "df102 ";//prints the first derivatives
        for (int ij=0;ij<reqrdnum;ij++){
            Rcout << Lld[ij] << " ";
        }
        Rcout << " " << endl;
        Rcout << "df103 ";//prints the second derivatives
        for (int ij=0;ij<reqrdnum;ij++){
            Rcout << Lldd[ij*reqrdnum+ij] << " ";
        }
        Rcout << " " << endl;
        Rcout << "df104 ";//prints parameter values
        for (int ij=0;ij<totalnum;ij++){
            Rcout << beta_0[ij] << " ";
        }
        Rcout << " " << endl;
    }
    //
    NumericVector Lldd_vec = wrap(Lldd);//
    Lldd_vec.attr("dim") = Dimension(reqrdnum, reqrdnum);
    //
    const Map<MatrixXd> Lldd_mat(as<Map<MatrixXd> >(Lldd_vec));
    //
    MatrixXd Lldd_inv = -1 * Lldd_mat.inverse().matrix(); //uses inverse information matrix to calculate the standard deviation
    //
    NumericVector stdev = wrap(Lldd_inv.diagonal().cwiseSqrt()); //vector of standard deviations
    // --------------------------
    // now a vector exists with row locations
    // --------------------------
    if (verbose){
        Rcout << "starting plot data " << endl;
    }
    MatrixXd residuals = MatrixXd::Zero(ntime,reqrdnum);
    MatrixXd res_scale = MatrixXd::Zero(ntime,reqrdnum);
    VectorXd res_df = VectorXd::Zero(ntime);
    //
    VectorXd req_beta = VectorXd::Zero(reqrdnum);
    for (int i=0; i<totalnum; i++){
        if (KeepConstant[i]==0){
            int j = i - sum(head(KeepConstant,i));
            req_beta[j] = beta_0[i];
        }
    }
    #pragma omp parallel for schedule(dynamic) num_threads(nthreads) collapse(2)
    for (int ijk=0;ijk<totalnum;ijk++){//totalnum*(totalnum+1)/2
        for (int j=0;j<ntime;j++){
            //
            if (KeepConstant[ijk]==0){
                //
                int ij = ijk - sum(head(KeepConstant,ijk));
                int df0_c = dfc[ijk]-1;
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
                //
                //
                // calculates the total term value
                //
                for (vector<double>::size_type i = 0; i < InGroup.size()-1; i=i+2){
                    t_sum += R.block(InGroup[i]-1,0,InGroup[i+1]-InGroup[i]+1,1).sum();
                    x_expect +=  (df0.block(InGroup[i]-1,df0_c,InGroup[i+1]-InGroup[i]+1,1).array() * R.block(InGroup[i]-1,0,InGroup[i+1]-InGroup[i]+1,1).array()).sum();
                }
                int dj = RiskFail(j,1)-RiskFail(j,0)+1;
                double x_risks = df0.block(RiskFail(j,0),df0_c,dj,1).sum()/dj; //calculate the average covariate value with events
                x_expect = x_expect / t_sum / dj; //calculates the averaged covariate value
                //
                residuals(j,ij) = (x_risks - x_expect);
                if (ij==0){
                    res_df(j) = dj;
                }
            }
        }
    }
    //
    res_scale = ((residuals * Lldd_inv) * ntime).array() + req_beta.transpose().replicate(residuals.rows(),1).array();
    //
    List res_list = List::create(_["residuals"]=wrap(residuals), _["scaled"]=wrap(res_scale), _["df"]=wrap(res_df));
    // returns residuals
    return res_list;
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
//' @param     double_step controls the step calculation, 0 for independent changes, 1 for solving b=Ax with complete matrices
//' @param     change_all    boolean if every parameter is being updated
//' @param     verbose    verbosity boolean
//' @param     debugging    debugging boolean
//' @param     KeepConstant    vector identifying constant parameters
//' @param     term_tot    total number of terms
//' @param     nthreads number of threads to use
//'
//' @return List of results: Log-likelihood of optimum, first derivative of log-likelihood, second derivative matrix, parameter list, standard deviation estimate, AIC, deviance, model information
// [[Rcpp::export]]
List LogLik_Poisson( MatrixXd PyrC, IntegerVector Term_n, StringVector tform, NumericVector a_n,NumericMatrix x_all,IntegerVector dfc,int fir, int der_iden,string modelform, double lr, int maxiter, int halfmax, double epsilon, double dbeta_cap, double abs_max,double dose_abs_max, double deriv_epsilon, int double_step,bool change_all, bool verbose, bool debugging, IntegerVector KeepConstant, int term_tot, int nthreads){
    ;
    //
    List temp_list = List::create(_["Status"]="FAILED"); //used as a dummy return value for code checking
    if (verbose){
        Rcout << "START_POISSON" << endl;
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
    int reqrdnum = totalnum - sum(KeepConstant);
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
    // int nthreads = Eigen::nbThreads()-1; //stores how many threads are allocated
    //
    //
    // Lld_worst: stores the highest magnitude log-likelihood derivative
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
    MatrixXd Td0 = MatrixXd::Zero(df0.rows(), reqrdnum); //preallocates matrix for Term derivative columns
    MatrixXd Tdd0 = MatrixXd::Zero(df0.rows(), reqrdnum*(reqrdnum+1)/2); //preallocates matrix for Term second derivative columns
    //
    MatrixXd Te = MatrixXd::Zero(df0.rows(), 1); //preallocates matrix for column terms used for temporary storage
    MatrixXd R = MatrixXd::Zero(df0.rows(), 1); //preallocates matrix for Risks
    ColXd Rd = MatrixXd::Zero(df0.rows(), reqrdnum); //preallocates matrix for Risk derivatives
    ColXd Rdd = MatrixXd::Zero(df0.rows(), reqrdnum*(reqrdnum+1)/2); //preallocates matrix for Risk second derivatives
    //
    MatrixXd Dose = MatrixXd::Constant(df0.rows(),term_tot,0.0);
    MatrixXd nonDose = MatrixXd::Constant(df0.rows(),term_tot,1.0); //Matrix of the total non-dose term values
    MatrixXd nonDose_LIN = MatrixXd::Constant(df0.rows(),term_tot,0.0);
    MatrixXd nonDose_PLIN = MatrixXd::Constant(df0.rows(),term_tot,1.0); //matrix of Loglinear subterm values
    MatrixXd nonDose_LOGLIN = MatrixXd::Constant(df0.rows(),term_tot,1.0); //matrix of Product linear subterm values
    MatrixXd TTerm = MatrixXd::Zero(Dose.rows(),Dose.cols()); //matrix of term values
    double dint = dose_abs_max; //The amount of change used to calculate derivatives in threshold paramters    //
    double dslp = abs_max;
    if (verbose){
        Rcout << "starting subterms " << term_tot << endl;
    }
    // Calculates the subterm and term values
    Make_subterms( totalnum, Term_n, tform, dfc, fir, T0, Td0, Tdd0, Dose, nonDose, TTerm,  nonDose_LIN, nonDose_PLIN, nonDose_LOGLIN ,beta_0, df0,dint,dslp,nthreads, debugging,KeepConstant);
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
        for (int ijk=0;ijk<reqrdnum;ijk++){
            Rcout << Td0.col(ijk).sum() << " ";
        }
        Rcout << " " << endl;
        Rcout << "second derivs checked ";
        for (int ijk=0;ijk<reqrdnum;ijk++){
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
    ColXd RdR = MatrixXd::Zero(df0.rows(), reqrdnum); //preallocates matrix for Risk to derivative ratios
    ColXd RddR = MatrixXd::Zero(df0.rows(), reqrdnum*(reqrdnum+1)/2);
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
    Make_Risks(modelform, tform, Term_n, totalnum, fir, T0, Td0, Tdd0, Te, R, Rd, Rdd, Dose, nonDose, TTerm, nonDose_LIN, nonDose_PLIN, nonDose_LOGLIN, RdR, RddR, nthreads, debugging,KeepConstant);
    //
    
    RdR = (RdR.array().isFinite()).select(RdR,0);
    RddR = (RddR.array().isFinite()).select(RddR,0);
    //
    //
    if (verbose){
        Rcout << "risk checked ";
        for (int ijk=0;ijk<1;ijk++){
            Rcout << R.col(0).sum() << " ";
        }
        Rcout << " " << endl;
        Rcout << "risk1 checked ";
        for (int ijk=0;ijk<reqrdnum;ijk++){
            Rcout << Rd.col(ijk).sum() << " ";
        }
        Rcout << " " << endl;
        Rcout << "risk2 checked ";
        for (int ijk=0;ijk<reqrdnum;ijk++){
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
    vector<double> Ll(reqrdnum,0.0);
    vector<double> Lld(reqrdnum,0.0);
    vector<double> Lldd(pow(reqrdnum,2),0.0);
    //
    //
    Poisson_LogLik( nthreads, totalnum, PyrC, R, Rd, Rdd, RdR, RddR, Ll, Lld, Lldd, debugging,KeepConstant);
    //
    MatrixXd dev_temp = MatrixXd::Zero(PyrC.rows(),2);
    dev_temp.col(0) = PyrC.col(0).array() * R.col(0).array();
    dev_temp.col(0) = PyrC.col(1).array() * dev_temp.col(0).array().pow(-1).array();
    dev_temp.col(0) = dev_temp.col(0).array().log().array();
    dev_temp.col(0) = PyrC.col(1).array() * dev_temp.col(0).array();
    dev_temp.col(1) = PyrC.col(1).array() - PyrC.col(0).array() * R.col(0).array();
    //
    dev_temp.col(0) = (dev_temp.col(0).array().isFinite()).select(dev_temp.col(0),0);
    //
    dev_temp.col(0) = dev_temp.col(0).array() - dev_temp.col(1).array();
    dev_temp.col(0) = (2 * dev_temp.col(0).array()).array();//.sqrt();
    dev_temp = (dev_temp.array().isFinite()).select(dev_temp,0);
    dev_temp.col(0) = (R.col(0).array()<0).select(0,dev_temp.col(0));
    double dev = dev_temp.col(0).sum(); //deviation calculation is split into steps
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
        for (int ij=0;ij<reqrdnum;ij++){
            Rcout << Ll[ij] << " ";
        }
        Rcout << " " << endl;
        Rcout << "df102 ";//prints the first derivatives
        for (int ij=0;ij<reqrdnum;ij++){
            Rcout << Lld[ij] << " ";
        }
        Rcout << " " << endl;
        Rcout << "df103 ";//prints the second derivatives
        for (int ij=0;ij<reqrdnum;ij++){
            Rcout << Lldd[ij*reqrdnum+ij] << " ";
        }
        Lld_worst=0;
        for (int ij=0;ij<reqrdnum;ij++){//locates highest magnitude derivative
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
        for (int ij=0;ij<reqrdnum;ij++){//prints the newton step value for zero derivative
            Rcout << Lld[ij]/Lldd[ij*reqrdnum+ij] << " ";
        }
        Rcout << " " << endl;
        Rcout << "df106 ";
        for (int ij=0;ij<reqrdnum;ij++){//prints the newton step value for zero log-likelihood
            Rcout << Ll[ij]/Lld[ij] << " ";
        }
        Rcout << " " << endl;
        Rcout << "df107 " << double_step << " " << abs_max << " " << dose_abs_max << " " << Ll_comp[0] << " " << Ll_comp[1] << endl;//prints several convergence terms
    }
    if (R.minCoeff()<0){
        Rcout << "A non-positive risk was detected: " << R.minCoeff() << endl;
        Rcout << "final failing values ";
        for (int ijk=0;ijk<totalnum;ijk++){
            Rcout << beta_0[ijk] << " ";
        }
        Rcout << " " << endl;
        //
        Rcout << "final failing terms ";
        for (int ijk=0;ijk<totalnum;ijk++){
            Rcout << tform[ijk] << " ";
        }
        Rcout << " " << endl;
        temp_list = List::create(_["LogLik"]=R_NaN,_["beta_0"]=wrap(beta_0) ,_["Deviation"]=R_NaN,_["Status"]="FAILED");
        return temp_list;
    }
    //
    vector<double> dbeta(totalnum,0.0);
    //
    // --------------------------
    // always starts from initial guess
    // --------------------------
    vector<double> beta_c(totalnum,0.0);
    vector<double> beta_a(totalnum,0.0);
    vector<double> beta_best(totalnum,0.0);
    vector<double> beta_p(totalnum,0.0);
    //int //risk_check_iter=0;
    VectorXd::Map(&beta_p[0], beta_0.size()) = beta_0;// stores previous parameters
    VectorXd::Map(&beta_c[0], beta_0.size()) = beta_0;// stores current parameters
    VectorXd::Map(&beta_a[0], beta_0.size()) = beta_0;// stores a refrence value for parameters
    VectorXd::Map(&beta_best[0], beta_0.size()) = beta_0;// stores the best parameters
    double halves = 0; //number of half-steps taken
    int ind0 = fir; //used for validations
    int iteration=0; //iteration number
    //
    int iter_stop =0; //tracks if the iterations should be stopped for convergence
    int iter_check=0; //signal to check for convergence
    bool convgd = FALSE;
    //
    if (sum(KeepConstant)==totalnum){
        iter_stop = 1;
        if (verbose){
            Rcout << "No parameters to change" << endl;
        }
    }
    //
    double Ll_abs_best = 10;
    vector<double> beta_abs_best(totalnum,0.0);
    while ((iteration < maxiter)&&(iter_stop==0)){
        iteration++;
        beta_p = beta_c;//
        beta_a = beta_c;//
        beta_best = beta_c;//
        //
        // Calcualtes the initial change in parameter
        Calc_Change( double_step, nthreads, totalnum, fir, der_iden, dbeta_cap, dose_abs_max, lr, abs_max, Ll, Lld, Lldd, dbeta, change_all, tform, dint,dslp, KeepConstant, debugging);
        Intercept_Bound(nthreads, totalnum, beta_0, dbeta, dfc, df0, KeepConstant, debugging, tform);
        if (verbose){
            Rcout << "Starting Halves"<<endl;//prints the final changes for validation
        }
        //
        //
        if ((Ll_abs_best>0)||(Ll_abs_best < Ll[ind0])){
            Ll_abs_best = Ll[ind0];
            beta_abs_best = beta_c;
        }
        //
        if (verbose){
            Rcout << "df501 " << Ll_abs_best << endl;
            Rcout << "df504 ";//prints parameter values
            for (int ij=0;ij<totalnum;ij++){
                Rcout << beta_abs_best[ij] << " ";
            }
            Rcout << " " << endl;
        }
        halves=0;
        while ((Ll[ind0] <= Ll_abs_best)&&(halves<halfmax)){ //repeats until half-steps maxed or an improvement
            Dose = MatrixXd::Zero(df0.rows(),term_tot);
            nonDose = MatrixXd::Constant(df0.rows(),term_tot,0.0);
            nonDose_LIN = MatrixXd::Constant(df0.rows(),term_tot,0.0);
            nonDose_PLIN = MatrixXd::Constant(df0.rows(),term_tot,1.0);
            nonDose_LOGLIN = MatrixXd::Constant(df0.rows(),term_tot,1.0);
            T0 = MatrixXd::Zero(df0.rows(), totalnum); //preallocates matrix for Term column
            Td0 = MatrixXd::Zero(df0.rows(), reqrdnum); //preallocates matrix for Term derivative columns
            Tdd0 = MatrixXd::Zero(df0.rows(), reqrdnum*(reqrdnum+1)/2); //preallocates matrix for Term second derivative columns
            for (int ij=0;ij<totalnum;ij++){
                beta_0[ij] = beta_a[ij] + dbeta[ij];
                beta_c[ij] = beta_0[ij];
            }
            // ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------//
            // The same subterm, risk, sides, and log-likelihood calculations are performed every half-step and iteration
            // ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------//
            Make_subterms( totalnum, Term_n, tform, dfc, fir, T0, Td0, Tdd0, Dose,  nonDose, TTerm,  nonDose_LIN, nonDose_PLIN, nonDose_LOGLIN ,beta_0, df0,dint,dslp,nthreads, debugging,KeepConstant);
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
                for (int ijk=0;ijk<reqrdnum;ijk++){
                    Rcout << Td0.col(ijk).sum() << " ";
                }
                Rcout << " " << endl;
                Rcout << "second derivs checked ";
                for (int ijk=0;ijk<reqrdnum;ijk++){
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
            RdR = MatrixXd::Zero(df0.rows(), reqrdnum);
            RddR = MatrixXd::Zero(df0.rows(), reqrdnum*(reqrdnum+1)/2);
            //
            //
            Make_Risks(modelform, tform, Term_n, totalnum, fir, T0, Td0, Tdd0, Te, R, Rd, Rdd, Dose, nonDose, TTerm, nonDose_LIN, nonDose_PLIN, nonDose_LOGLIN, RdR, RddR, nthreads, debugging,KeepConstant);
            //
            RdR = (RdR.array().isFinite()).select(RdR,0);
            RddR = (RddR.array().isFinite()).select(RddR,0);
            //
            if ((R.minCoeff()<0)&&(TRUE)){
                #pragma omp parallel for num_threads(nthreads)
                for (int ijk=0;ijk<totalnum;ijk++){
                    int tij = Term_n[ijk];
                    if (TTerm.col(tij).minCoeff()<=0){
                        dbeta[ijk] = dbeta[ijk] / 2.0;
                    } else if (isinf(TTerm.col(tij).maxCoeff())){
                        dbeta[ijk] = dbeta[ijk] / 2.0;
                    }
                }
                if (verbose){
                    Rcout << "A non-positive risk was detected: " << R.minCoeff() << endl;
                }
                halves+=0.2;
            } else {
                halves++;
                if (verbose){
                    Rcout << "risk checked ";
                    for (int ijk=0;ijk<1;ijk++){
                        Rcout << R.col(0).sum() << " ";
                    }
                    Rcout << " " << endl;
                    Rcout << "risk1 checked ";
                    for (int ijk=0;ijk<reqrdnum;ijk++){
                        Rcout << Rd.col(ijk).sum() << " ";
                    }
                    Rcout << " " << endl;
                    Rcout << "risk2 checked ";
                    for (int ijk=0;ijk<reqrdnum;ijk++){
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
                Poisson_LogLik( nthreads, totalnum, PyrC, R, Rd, Rdd, RdR, RddR, Ll, Lld, Lldd, debugging,KeepConstant);
                dev_temp.col(0) = PyrC.col(0).array() * R.col(0).array();
                dev_temp.col(0) = PyrC.col(1).array() * dev_temp.col(0).array().pow(-1).array();
                dev_temp.col(0) = dev_temp.col(0).array().log().array();
                dev_temp.col(0) = PyrC.col(1).array() * dev_temp.col(0).array();
                dev_temp.col(1) = PyrC.col(1).array() - PyrC.col(0).array() * R.col(0).array();
                //
                dev_temp.col(0) = (dev_temp.col(0).array().isFinite()).select(dev_temp.col(0),0);
                //
                dev_temp.col(0) = dev_temp.col(0).array() - dev_temp.col(1).array();
                dev_temp.col(0) = (2 * dev_temp.col(0).array()).array();//.sqrt();
                dev_temp = (dev_temp.array().isFinite()).select(dev_temp,0);
                dev_temp.col(0) = (R.col(0).array()<0).select(0,dev_temp.col(0));
                dev = dev_temp.col(0).sum(); //deviation calculation is split into steps
                if (verbose){
                    end_point = system_clock::now();
                    ending = time_point_cast<microseconds>(end_point).time_since_epoch().count();
                    Rcout <<"df100 "<<(ending-start)<<" "<<0<<" "<<0<<" "<<0<<",Half Calc"<<endl;//prints the time
                    gibtime = system_clock::to_time_t(system_clock::now());
                    Rcout << ctime(&gibtime) << endl;
                    Rcout << "df101 ";//prints the log-likelihoods
                    for (int ij=0;ij<reqrdnum;ij++){
                        Rcout << Ll[ij] << " ";
                    }
                    Rcout << " " << endl;
                    Rcout << "df102 ";//prints the first derivatives
                    for (int ij=0;ij<reqrdnum;ij++){
                        Rcout << Lld[ij] << " ";
                    }
                    Rcout << " " << endl;
                    Rcout << "df103 ";//prints the second derivatives
                    for (int ij=0;ij<reqrdnum;ij++){
                        Rcout << Lldd[ij*reqrdnum+ij] << " ";
                    }
                    Lld_worst=0;
                    for (int ij=0;ij<reqrdnum;ij++){//locates highest magnitude derivative
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
                    for (int ij=0;ij<reqrdnum;ij++){//prints the newton step value for zero derivative
                        Rcout << Lld[ij]/Lldd[ij*reqrdnum+ij] << " ";
                    }
                    Rcout << " " << endl;
                    Rcout << "df106 ";
                    for (int ij=0;ij<reqrdnum;ij++){//prints the newton step value for zero log-likelihood
                        Rcout << Ll[ij]/Lld[ij] << " ";
                    }
                    Rcout << " " << endl;
                    Rcout << "df107 " << double_step << " " << abs_max << " " << dose_abs_max << " " << Ll_comp[0] << " " << Ll_comp[1] << endl;//prints several convergence terms
                }
                if (change_all){
                    if (Ll[ind0] <= Ll_abs_best){//takes a half-step if needed
                        #pragma omp parallel for num_threads(nthreads)
                        for (int ijk=0;ijk<totalnum;ijk++){
                            dbeta[ijk] = dbeta[ijk] * 0.5; //
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
                }
                #pragma omp parallel for num_threads(nthreads)
                for (int ijk=0;ijk<totalnum;ijk++){//totalnum*(totalnum+1)/2
                    beta_0[ijk] = beta_c[ijk];
                }
            }
        }
        if (beta_best!=beta_c){//if the risk matrices aren't the optimal values, then they must be recalculated
            if (verbose){
                Rcout << "Changing back to best"<<endl;
            }
            // If it goes through every half step without improvement, then the maximum change needs to be decreased
            abs_max = abs_max*pow(0.5,halfmax);
            dose_abs_max = dose_abs_max*pow(0.5,halfmax);
            iter_check = 1;
            //
            beta_p = beta_best;//
            beta_a = beta_best;//
            beta_c = beta_best;//
            Dose = MatrixXd::Zero(df0.rows(),term_tot);
            nonDose = MatrixXd::Constant(df0.rows(),term_tot,0.0);
            nonDose_LIN = MatrixXd::Constant(df0.rows(),term_tot,0.0);
            nonDose_PLIN = MatrixXd::Constant(df0.rows(),term_tot,1.0);
            nonDose_LOGLIN = MatrixXd::Constant(df0.rows(),term_tot,1.0);
            for (int ij=0;ij<totalnum;ij++){
                beta_0[ij] = beta_best[ij];
            }
            Make_subterms( totalnum, Term_n, tform, dfc, fir, T0, Td0, Tdd0, Dose, nonDose, TTerm,  nonDose_LIN, nonDose_PLIN, nonDose_LOGLIN,beta_0, df0,dint,dslp,nthreads, debugging,KeepConstant);
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
                for (int ijk=0;ijk<reqrdnum;ijk++){
                    Rcout << Td0.col(ijk).sum() << " ";
                }
                Rcout << " " << endl;
                Rcout << "second derivs checked ";
                for (int ijk=0;ijk<reqrdnum;ijk++){
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
            RdR = MatrixXd::Zero(df0.rows(), reqrdnum);
            RddR = MatrixXd::Zero(df0.rows(), reqrdnum*(reqrdnum+1)/2);
            //
            //
            Make_Risks(modelform, tform, Term_n, totalnum, fir, T0, Td0, Tdd0, Te, R, Rd, Rdd, Dose, nonDose, TTerm, nonDose_LIN, nonDose_PLIN, nonDose_LOGLIN, RdR, RddR, nthreads, debugging,KeepConstant);
            R = (R.array().isFinite()).select(R,0);
            Rd = (Rd.array().isFinite()).select(Rd,0);
            Rdd = (Rdd.array().isFinite()).select(Rdd,0);
            RdR = (RdR.array().isFinite()).select(RdR,0);
            RddR = (RddR.array().isFinite()).select(RddR,0);
            //
            temp_list = List::create(_["betas"]=wrap(beta_0),_["Status"]="FAILED"); //used as a dummy return value for code checking
            if (verbose){
                Rcout << "risk checked ";
                for (int ijk=0;ijk<1;ijk++){
                    Rcout << R.col(0).sum() << " ";
                }
                Rcout << " " << endl;
                Rcout << "risk1 checked ";
                for (int ijk=0;ijk<reqrdnum;ijk++){
                    Rcout << Rd.col(ijk).sum() << " ";
                }
                Rcout << " " << endl;
                Rcout << "risk2 checked ";
                for (int ijk=0;ijk<reqrdnum;ijk++){
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
            Poisson_LogLik( nthreads, totalnum, PyrC, R, Rd, Rdd, RdR, RddR, Ll, Lld, Lldd, debugging,KeepConstant);
            dev_temp.col(0) = PyrC.col(0).array() * R.col(0).array();
            dev_temp.col(0) = PyrC.col(1).array() * dev_temp.col(0).array().pow(-1).array();
            dev_temp.col(0) = dev_temp.col(0).array().log().array();
            dev_temp.col(0) = PyrC.col(1).array() * dev_temp.col(0).array();
            dev_temp.col(1) = PyrC.col(1).array() - PyrC.col(0).array() * R.col(0).array();
            //
            dev_temp.col(0) = (dev_temp.col(0).array().isFinite()).select(dev_temp.col(0),0);
            //
            dev_temp.col(0) = dev_temp.col(0).array() - dev_temp.col(1).array();
            dev_temp.col(0) = (2 * dev_temp.col(0).array()).array();//.sqrt();
            dev_temp = (dev_temp.array().isFinite()).select(dev_temp,0);
            dev_temp.col(0) = (R.col(0).array()<0).select(0,dev_temp.col(0));
            dev = dev_temp.col(0).sum(); //deviation calculation is split into steps
            
        }
        Lld_worst=0;
        for (int ij=0;ij<reqrdnum;ij++){
            if (abs(Lld[ij]) > Lld_worst){
                Lld_worst = abs(Lld[ij]);
            }
        }
        if (iteration > reqrdnum){//Sets the minimum number of iterations
            if ((iteration % (reqrdnum))||(iter_check==1)){//Checks every set number of iterations
                iter_check=0;
                if (Lld_worst < deriv_epsilon){//ends if the derivatives are low enough
                    iter_stop = 1;
                    convgd = TRUE;
                }
                Ll_comp[1]=Ll[0];
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
            for (int ij=0;ij<reqrdnum;ij++){
                Rcout << Ll[ij] << " ";
            }
            Rcout << " " << endl;
            Rcout << "df102 ";
            for (int ij=0;ij<reqrdnum;ij++){
                Rcout << Lld[ij] << " ";
            }
            Rcout << " " << endl;
            Rcout << "df103 ";
            for (int ij=0;ij<reqrdnum;ij++){
                Rcout << Lldd[ij*reqrdnum+ij] << " ";
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
            for (int ij=0;ij<reqrdnum;ij++){
                Rcout << Lld[ij]/Lldd[ij*reqrdnum+ij] << " ";
            }
            Rcout << " " << endl;
            Rcout << "df106 ";
            for (int ij=0;ij<reqrdnum;ij++){
                Rcout << Ll[ij]/Lld[ij] << " ";
            }
            Rcout << " " << endl;
            Rcout << "df107 " << double_step << " " << abs_max << " " << dose_abs_max << " " << Ll_comp[0] << " " << Ll_comp[1] << endl;
        }
    }
    // -----------------------------------------------
    // Performing Full Calculation to get full second derivative matrix
    // -----------------------------------------------
    fill(Ll.begin(), Ll.end(), 0.0);
    fill(Lld.begin(), Lld.end(), 0.0);
    fill(Lldd.begin(), Lldd.end(), 0.0);
    Poisson_LogLik( nthreads, totalnum, PyrC, R, Rd, Rdd, RdR, RddR, Ll, Lld, Lldd, debugging,KeepConstant);
    dev_temp.col(0) = PyrC.col(0).array() * R.col(0).array();
    dev_temp.col(0) = PyrC.col(1).array() * dev_temp.col(0).array().pow(-1).array();
    dev_temp.col(0) = dev_temp.col(0).array().log().array();
    dev_temp.col(0) = PyrC.col(1).array() * dev_temp.col(0).array();
    dev_temp.col(1) = PyrC.col(1).array() - PyrC.col(0).array() * R.col(0).array();
    //
    dev_temp.col(0) = (dev_temp.col(0).array().isFinite()).select(dev_temp.col(0),0);
    //
    dev_temp.col(0) = dev_temp.col(0).array() - dev_temp.col(1).array();
    dev_temp.col(0) = (2 * dev_temp.col(0).array()).array();//.sqrt();
    dev_temp = (dev_temp.array().isFinite()).select(dev_temp,0);
    dev_temp.col(0) = (R.col(0).array()<0).select(0,dev_temp.col(0));
    dev = dev_temp.col(0).sum(); //deviation calculation is split into steps
    if (verbose){
        end_point = system_clock::now();
        ending = time_point_cast<microseconds>(end_point).time_since_epoch().count();
        Rcout <<"df100 "<<(ending-start)<<" "<<halves<<" "<<iteration<<" "<<ind0<<",Recalc"<<endl;
        gibtime = system_clock::to_time_t(system_clock::now());
        Rcout << ctime(&gibtime) << endl;
        Rcout << "df101 ";
        for (int ij=0;ij<reqrdnum;ij++){
            Rcout << Ll[ij] << " ";
        }
        Rcout << " " << endl;
        Rcout << "df102 ";
        for (int ij=0;ij<reqrdnum;ij++){
            Rcout << Lld[ij] << " ";
        }
        Rcout << " " << endl;
        Rcout << "df103 ";
        for (int ij=0;ij<reqrdnum;ij++){
            Rcout << Lldd[ij*reqrdnum+ij] << " ";
        }
        Rcout << " " << endl;
        Rcout << "df104 ";
        for (int ij=0;ij<totalnum;ij++){
            Rcout << beta_c[ij] << " ";
        }
        Rcout << " " << endl;
        Rcout << "df105 ";
        for (int ij=0;ij<reqrdnum;ij++){//prints the newton step value for zero derivative
            Rcout << Lld[ij]/Lldd[ij*reqrdnum+ij] << " ";
        }
        Rcout << " " << endl;
        Rcout << "df106 ";
        for (int ij=0;ij<reqrdnum;ij++){//prints the newton step value for zero log-likelihood
            Rcout << Ll[ij]/Lld[ij] << " ";
        }
        Rcout << " " << endl;
        Rcout << "df107 " << double_step << " " << abs_max << " " << dose_abs_max << " " << Ll_comp[0] << " " << Ll_comp[1] << endl;//prints several convergence terms
        Rcout << "Checking Deviance " << dev << endl;
        Rcout << "Finshed iteration" << endl;
    }
    //
    if ((Ll_abs_best>0)||(Ll_abs_best < Ll[ind0])){
        Ll_abs_best = Ll[ind0];
        beta_abs_best = beta_c;
    }
    //
    if (verbose){
        Rcout << "df501 " << Ll_abs_best << endl;
        Rcout << "df504 ";//prints parameter values
        for (int ij=0;ij<totalnum;ij++){
            Rcout << beta_abs_best[ij] << " ";
        }
        Rcout << " " << endl;
    }
    // Changes the parameter back into the original form
    List para_list = List::create(_["Term_n"]=Term_n,_["tforms"]=tform);
    List control_list = List::create(_["Iteration"]=iteration, _["Maximum Step"]=abs_max, _["Derivative Limiting"]=Lld_worst);
    //
    int kept_covs = totalnum - sum(KeepConstant);
    NumericVector Lldd_vec(kept_covs * kept_covs);
    #pragma omp parallel for schedule(dynamic) num_threads(nthreads)
    for (int ijk=0;ijk<totalnum*(totalnum+1)/2;ijk++){
        int ij = 0;
        int jk = ijk;
        int pij_ind=-100;
        int pjk_ind=-100;
        while (jk>ij){
            ij++;
            jk-=ij;
        }
        if (KeepConstant[ij]==0){
            pij_ind = ij - sum(head(KeepConstant,ij));
            if (KeepConstant[jk]==0){
                pjk_ind = jk - sum(head(KeepConstant,jk));
                Lldd_vec[pij_ind * kept_covs + pjk_ind]=Lldd[ij*totalnum+jk];
            }
        }
    }
    for (int ijk=0;ijk<kept_covs*(kept_covs+1)/2;ijk++){
        int ij = 0;
        int jk = ijk;
        while (jk>ij){
            ij++;
            jk-=ij;
        }
        Lldd_vec[jk * kept_covs + ij]=Lldd_vec[ij * kept_covs + jk];
    }
    Lldd_vec.attr("dim") = Dimension(kept_covs, kept_covs);
    const Map<MatrixXd> Lldd_mat(as<Map<MatrixXd> >(Lldd_vec));
    //
    MatrixXd Lldd_inv = -1 * Lldd_mat.inverse().matrix(); //uses inverse information matrix to calculate the standard deviation
    VectorXd stdev = VectorXd::Zero(totalnum);
    for (int ij=0;ij<totalnum;ij++){
        if (KeepConstant[ij]==0){
            int pij_ind = ij - sum(head(KeepConstant,ij));
            stdev(ij) = sqrt(Lldd_inv(pij_ind,pij_ind));
        }
    }
    //
    List res_list = List::create(_["LogLik"]=wrap(Ll[0]),_["First_Der"]=wrap(Lld),_["Second_Der"]=Lldd_vec,_["beta_0"]=wrap(beta_0) ,_["Standard_Deviation"]=wrap(stdev) ,_["AIC"]=2*(totalnum-accumulate(KeepConstant.begin(),KeepConstant.end(), 0.0))+dev,_["Deviation"]=dev,_["Parameter_Lists"]=para_list,_["Control_List"]=control_list,_["Converged"]=convgd);
    // returns a list of results
    return res_list;
}

//' Primary poisson regression function with strata effect
//' \code{LogLik_Poisson_STRATA} Performs the calls to calculation functions, Structures the poisson regression, With verbose option prints out time stamps and intermediate sums of terms and derivatives, all with strata effects
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
//' @param     double_step controls the step calculation, 0 for independent changes, 1 for solving b=Ax with complete matrices
//' @param     change_all    boolean if every parameter is being updated
//' @param     verbose    verbosity boolean
//' @param     debugging    debugging boolean
//' @param     KeepConstant    vector identifying constant parameters
//' @param     term_tot    total number of terms
//' @param     dfs matrix of stratification variables
//' @param     keep_strata boolean to return the strata parameter values
//' @param     nthreads number of threads to use
//'
//' @return List of results: Log-likelihood of optimum, first derivative of log-likelihood, second derivative matrix, parameter list, standard deviation estimate, AIC, deviance, model information
// [[Rcpp::export]]
List LogLik_Poisson_STRATA( MatrixXd PyrC, IntegerVector Term_n, StringVector tform, NumericVector a_n,NumericMatrix x_all,IntegerVector dfc,int fir, int der_iden,string modelform, double lr, int maxiter, int halfmax, double epsilon, double dbeta_cap, double abs_max,double dose_abs_max, double deriv_epsilon, int double_step,bool change_all, bool verbose, bool debugging, IntegerVector KeepConstant, int term_tot, const MatrixXd& dfs, bool keep_strata, int nthreads){
    ;
    //
    List temp_list = List::create(_["Status"]="FAILED"); //used as a dummy return value for code checking
    if (verbose){
        Rcout << "START_POISSON_STRATA" << endl;
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
    int reqrdnum = totalnum - sum(KeepConstant);
    //
    if (verbose){
        Rcout << "Term checked ";
        for (int ij=0;ij<totalnum;ij++){
            Rcout << Term_n[ij] << " ";
        }
        Rcout << " " << endl;
        Rcout << "Constant checked ";
        for (int ij=0;ij<totalnum;ij++){
            Rcout << KeepConstant[ij] << " ";
        }
        Rcout << " " << endl;
    }
    //
    //
    // cout.precision: controls the number of significant digits printed
    // nthreads: number of threads used for parallel operations
    //
    Rcout.precision(7); //forces higher precision numbers printed to terminal
    // int nthreads = Eigen::nbThreads()-1; //stores how many threads are allocated
    //
    //
    // Lld_worst: stores the highest magnitude log-likelihood derivative
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
    MatrixXd Td0 = MatrixXd::Zero(df0.rows(), reqrdnum); //preallocates matrix for Term derivative columns
    MatrixXd Tdd0 = MatrixXd::Zero(df0.rows(), reqrdnum*(reqrdnum+1)/2); //preallocates matrix for Term second derivative columns
    //
    MatrixXd Te = MatrixXd::Zero(df0.rows(), 1); //preallocates matrix for column terms used for temporary storage
    MatrixXd R = MatrixXd::Zero(df0.rows(), 1); //preallocates matrix for Risks
    ColXd Rd = MatrixXd::Zero(df0.rows(), reqrdnum); //preallocates matrix for Risk derivatives
    ColXd Rdd = MatrixXd::Zero(df0.rows(), reqrdnum*(reqrdnum+1)/2); //preallocates matrix for Risk second derivatives
    //
    MatrixXd Dose = MatrixXd::Constant(df0.rows(),term_tot,0.0);
    MatrixXd nonDose = MatrixXd::Constant(df0.rows(),term_tot,1.0); //Matrix of the total non-dose term values
    MatrixXd nonDose_LIN = MatrixXd::Constant(df0.rows(),term_tot,0.0);
    MatrixXd nonDose_PLIN = MatrixXd::Constant(df0.rows(),term_tot,1.0); //matrix of Loglinear subterm values
    MatrixXd nonDose_LOGLIN = MatrixXd::Constant(df0.rows(),term_tot,1.0); //matrix of Product linear subterm values
    MatrixXd TTerm = MatrixXd::Zero(Dose.rows(),Dose.cols()); //matrix of term values
    double dint = dose_abs_max; //The amount of change used to calculate derivatives in threshold paramters    //
    double dslp = abs_max;
    if (verbose){
        Rcout << "starting subterms " << term_tot << endl;
    }
    // Calculates the subterm and term values
    Make_subterms( totalnum, Term_n, tform, dfc, fir, T0, Td0, Tdd0, Dose, nonDose, TTerm,  nonDose_LIN, nonDose_PLIN, nonDose_LOGLIN ,beta_0, df0,dint,dslp,nthreads, debugging,KeepConstant);
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
        for (int ijk=0;ijk<reqrdnum;ijk++){
            Rcout << Td0.col(ijk).sum() << " ";
        }
        Rcout << " " << endl;
        Rcout << "second derivs checked ";
        for (int ijk=0;ijk<reqrdnum;ijk++){
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
    ColXd RdR = MatrixXd::Zero(df0.rows(), reqrdnum); //preallocates matrix for Risk to derivative ratios
    ColXd RddR = MatrixXd::Zero(df0.rows(), reqrdnum*(reqrdnum+1)/2);
    //
    if (verbose){
        end_point = system_clock::now();
        ending = time_point_cast<microseconds>(end_point).time_since_epoch().count();
        Rcout <<"df99,"<<(ending-start)<<",Prep_Terms"<<endl;
        gibtime = system_clock::to_time_t(system_clock::now());
        Rcout << ctime(&gibtime) << endl;
    }
    // generate weight
    VectorXd s_weights = VectorXd::Zero(df0.rows());
    Gen_Strat_Weight(modelform, dfs, PyrC, s_weights, nthreads, tform, Term_n, term_tot);
    //
    // Calculates the risk for each row
    Make_Risks_Weighted(modelform, tform, Term_n, totalnum, fir, s_weights, T0, Td0, Tdd0, Te, R, Rd, Rdd, Dose, nonDose, TTerm, nonDose_LIN, nonDose_PLIN, nonDose_LOGLIN, RdR, RddR, nthreads, debugging,KeepConstant);
    //
    
    RdR = (RdR.array().isFinite()).select(RdR,0);
    RddR = (RddR.array().isFinite()).select(RddR,0);
    //
    //
    if (verbose){
        Rcout << "risk checked ";
        for (int ijk=0;ijk<1;ijk++){
            Rcout << R.col(0).sum() << " ";
        }
        Rcout << " " << endl;
        Rcout << "risk1 checked ";
        for (int ijk=0;ijk<reqrdnum;ijk++){
            Rcout << Rd.col(ijk).sum() << " ";
        }
        Rcout << " " << endl;
        Rcout << "risk2 checked ";
        for (int ijk=0;ijk<reqrdnum;ijk++){
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
    vector<double> Ll(reqrdnum,0.0);
    vector<double> Lld(reqrdnum,0.0);
    vector<double> Lldd(pow(reqrdnum,2),0.0);
    //
    //
    Poisson_LogLik( nthreads, totalnum, PyrC, R, Rd, Rdd, RdR, RddR, Ll, Lld, Lldd, debugging,KeepConstant);
    //
    MatrixXd dev_temp = MatrixXd::Zero(PyrC.rows(),2);
    dev_temp.col(0) = PyrC.col(0).array() * R.col(0).array();
    dev_temp.col(0) = PyrC.col(1).array() * dev_temp.col(0).array().pow(-1).array();
    dev_temp.col(0) = dev_temp.col(0).array().log().array();
    dev_temp.col(0) = PyrC.col(1).array() * dev_temp.col(0).array();
    dev_temp.col(1) = PyrC.col(1).array() - PyrC.col(0).array() * R.col(0).array();
    //
    dev_temp.col(0) = (dev_temp.col(0).array().isFinite()).select(dev_temp.col(0),0);
    //
    dev_temp.col(0) = dev_temp.col(0).array() - dev_temp.col(1).array();
    dev_temp.col(0) = (2 * dev_temp.col(0).array()).array();//.sqrt();
    dev_temp = (dev_temp.array().isFinite()).select(dev_temp,0);
    dev_temp.col(0) = (R.col(0).array()<0).select(0,dev_temp.col(0));
    double dev = dev_temp.col(0).sum(); //deviation calculation is split into steps
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
        for (int ij=0;ij<reqrdnum;ij++){
            Rcout << Ll[ij] << " ";
        }
        Rcout << " " << endl;
        Rcout << "df102 ";//prints the first derivatives
        for (int ij=0;ij<reqrdnum;ij++){
            Rcout << Lld[ij] << " ";
        }
        Rcout << " " << endl;
        Rcout << "df103 ";//prints the second derivatives
        for (int ij=0;ij<reqrdnum;ij++){
            Rcout << Lldd[ij*reqrdnum+ij] << " ";
        }
        Lld_worst=0;
        for (int ij=0;ij<reqrdnum;ij++){//locates highest magnitude derivative
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
        Rcout << "Checking Best Deviance " << dev << endl;
        Rcout << "df105 ";
        for (int ij=0;ij<reqrdnum;ij++){//prints the newton step value for zero derivative
            Rcout << Lld[ij]/Lldd[ij*reqrdnum+ij] << " ";
        }
        Rcout << " " << endl;
        Rcout << "df106 ";
        for (int ij=0;ij<reqrdnum;ij++){//prints the newton step value for zero log-likelihood
            Rcout << Ll[ij]/Lld[ij] << " ";
        }
        Rcout << " " << endl;
        Rcout << "df107 " << double_step << " " << abs_max << " " << dose_abs_max << " " << Ll_comp[0] << " " << Ll_comp[1] << endl;//prints several convergence terms
    }
    if (R.minCoeff()<0){
        Rcout << "A non-positive risk was detected: " << R.minCoeff() << endl;
        Rcout << "final failing values ";
        for (int ijk=0;ijk<totalnum;ijk++){
            Rcout << beta_0[ijk] << " ";
        }
        Rcout << " " << endl;
        //
        Rcout << "final failing terms ";
        for (int ijk=0;ijk<totalnum;ijk++){
            Rcout << tform[ijk] << " ";
        }
        Rcout << " " << endl;
        //
        //
        temp_list = List::create(_["LogLik"]=R_NaN,_["First_Der"]=wrap(Lld),_["beta_0"]=wrap(beta_0) ,_["Deviation"]=R_NaN,_["Status"]="FAILED");
        return temp_list;
    }
    //
    vector<double> dbeta(totalnum,0.0);
    //
    // --------------------------
    // always starts from initial guess
    // --------------------------
    vector<double> beta_c(totalnum,0.0);
    vector<double> beta_a(totalnum,0.0);
    vector<double> beta_best(totalnum,0.0);
    vector<double> beta_p(totalnum,0.0);
    //int //risk_check_iter=0;
    VectorXd::Map(&beta_p[0], beta_0.size()) = beta_0;// stores previous parameters
    VectorXd::Map(&beta_c[0], beta_0.size()) = beta_0;// stores current parameters
    VectorXd::Map(&beta_a[0], beta_0.size()) = beta_0;// stores a refrence value for parameters
    VectorXd::Map(&beta_best[0], beta_0.size()) = beta_0;// stores the best parameters
    double halves = 0; //number of half-steps taken
    int ind0 = fir; //used for validations
    int iteration=0; //iteration number
    //
    int iter_stop =0; //tracks if the iterations should be stopped for convergence
    int iter_check=0; //signal to check for convergence
    bool convgd = FALSE;
    //
    if (sum(KeepConstant)==totalnum){
        iter_stop = 1;
        if (verbose){
            Rcout << "No parameters to change" << endl;
        }
    }
    //
    double Ll_abs_best = 10;
    vector<double> beta_abs_best(totalnum,0.0);
    while ((iteration < maxiter)&&(iter_stop==0)){
        iteration++;
        beta_p = beta_c;//
        beta_a = beta_c;//
        beta_best = beta_c;//
        //
        // Calcualtes the initial change in parameter
        Calc_Change( double_step, nthreads, totalnum, fir, der_iden, dbeta_cap, dose_abs_max, lr, abs_max, Ll, Lld, Lldd, dbeta, change_all, tform, dint,dslp, KeepConstant, debugging);
        Intercept_Bound(nthreads, totalnum, beta_0, dbeta, dfc, df0, KeepConstant, debugging, tform);
        if (verbose){
            Rcout << "Starting Halves"<<endl;//prints the final changes for validation
//            for (int ijk=0;ijk<totalnum;ijk++){
//                Rcout << dbeta[ijk] << " ";
//            }
//            Rcout << " " << endl;
        }
        //
        //
        if ((Ll_abs_best>0)||(Ll_abs_best < Ll[ind0])){
            Ll_abs_best = Ll[ind0];
            beta_abs_best = beta_c;
        }
        //
        if (verbose){
            Rcout << "df501 " << Ll_abs_best << endl;
            Rcout << "df504 ";//prints parameter values
            for (int ij=0;ij<totalnum;ij++){
                Rcout << beta_abs_best[ij] << " ";
            }
            Rcout << " " << endl;
        }
        halves=0;
        while ((Ll[ind0] <= Ll_abs_best)&&(halves<halfmax)){ //repeats until half-steps maxed or an improvement
            Dose = MatrixXd::Zero(df0.rows(),term_tot);
            nonDose = MatrixXd::Constant(df0.rows(),term_tot,0.0);
            nonDose_LIN = MatrixXd::Constant(df0.rows(),term_tot,0.0);
            nonDose_PLIN = MatrixXd::Constant(df0.rows(),term_tot,1.0);
            nonDose_LOGLIN = MatrixXd::Constant(df0.rows(),term_tot,1.0);
            T0 = MatrixXd::Zero(df0.rows(), totalnum); //preallocates matrix for Term column
            Td0 = MatrixXd::Zero(df0.rows(), reqrdnum); //preallocates matrix for Term derivative columns
            Tdd0 = MatrixXd::Zero(df0.rows(), reqrdnum*(reqrdnum+1)/2); //preallocates matrix for Term second derivative columns
            if (verbose){
                Rcout << "Starting Half Step "<<halves<<endl;//prints the final changes for validation
//                for (int ijk=0;ijk<totalnum;ijk++){
//                    Rcout << dbeta[ijk] << " ";
//                }
//                Rcout << " " << endl;
            }
            for (int ij=0;ij<totalnum;ij++){
                beta_0[ij] = beta_a[ij] + dbeta[ij];
                beta_c[ij] = beta_0[ij];
            }
            // ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------//
            // The same subterm, risk, sides, and log-likelihood calculations are performed every half-step and iteration
            // ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------//
            Make_subterms( totalnum, Term_n, tform, dfc, fir, T0, Td0, Tdd0, Dose,  nonDose, TTerm,  nonDose_LIN, nonDose_PLIN, nonDose_LOGLIN ,beta_0, df0,dint,dslp,nthreads, debugging,KeepConstant);
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
                for (int ijk=0;ijk<reqrdnum;ijk++){
                    Rcout << Td0.col(ijk).sum() << " ";
                }
                Rcout << " " << endl;
                Rcout << "second derivs checked ";
                for (int ijk=0;ijk<reqrdnum;ijk++){
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
            RdR = MatrixXd::Zero(df0.rows(), reqrdnum);
            RddR = MatrixXd::Zero(df0.rows(), reqrdnum*(reqrdnum+1)/2);
            //
            //
            Make_Risks_Weighted(modelform, tform, Term_n, totalnum, fir, s_weights, T0, Td0, Tdd0, Te, R, Rd, Rdd, Dose, nonDose, TTerm, nonDose_LIN, nonDose_PLIN, nonDose_LOGLIN, RdR, RddR, nthreads, debugging,KeepConstant);
            //
            RdR = (RdR.array().isFinite()).select(RdR,0);
            RddR = (RddR.array().isFinite()).select(RddR,0);
            //
            if ((R.minCoeff()<0)&&(TRUE)){
                #pragma omp parallel for num_threads(nthreads)
                for (int ijk=0;ijk<totalnum;ijk++){
                    int tij = Term_n[ijk];
                    if (TTerm.col(tij).minCoeff()<=0){
                        dbeta[ijk] = dbeta[ijk] / 2.0;
                    } else if (isinf(TTerm.col(tij).maxCoeff())){
                        dbeta[ijk] = dbeta[ijk] / 2.0;
                    }
                }
                if (verbose){
                    Rcout << "A non-positive risk was detected: " << R.minCoeff() << endl;
                }
                halves+=0.2;
            } else {
                halves++;
                if (verbose){
                    Rcout << "risk checked ";
                    for (int ijk=0;ijk<1;ijk++){
                        Rcout << R.col(0).sum() << " ";
                    }
                    Rcout << " " << endl;
                    Rcout << "risk1 checked ";
                    for (int ijk=0;ijk<reqrdnum;ijk++){
                        Rcout << Rd.col(ijk).sum() << " ";
                    }
                    Rcout << " " << endl;
                    Rcout << "risk2 checked ";
                    for (int ijk=0;ijk<reqrdnum;ijk++){
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
                Poisson_LogLik( nthreads, totalnum, PyrC, R, Rd, Rdd, RdR, RddR, Ll, Lld, Lldd, debugging,KeepConstant);
                dev_temp.col(0) = PyrC.col(0).array() * R.col(0).array();
                dev_temp.col(0) = PyrC.col(1).array() * dev_temp.col(0).array().pow(-1).array();
                dev_temp.col(0) = dev_temp.col(0).array().log().array();
                dev_temp.col(0) = PyrC.col(1).array() * dev_temp.col(0).array();
                dev_temp.col(1) = PyrC.col(1).array() - PyrC.col(0).array() * R.col(0).array();
                //
                dev_temp.col(0) = (dev_temp.col(0).array().isFinite()).select(dev_temp.col(0),0);
                //
                dev_temp.col(0) = dev_temp.col(0).array() - dev_temp.col(1).array();
                dev_temp.col(0) = (2 * dev_temp.col(0).array()).array();//.sqrt();
                dev_temp = (dev_temp.array().isFinite()).select(dev_temp,0);
                dev_temp.col(0) = (R.col(0).array()<0).select(0,dev_temp.col(0));
                dev = dev_temp.col(0).sum(); //deviation calculation is split into steps
                //
                if (verbose){
                    end_point = system_clock::now();
                    ending = time_point_cast<microseconds>(end_point).time_since_epoch().count();
                    Rcout <<"df100 "<<(ending-start)<<" "<<0<<" "<<0<<" "<<0<<",Half Calc"<<endl;//prints the time
                    gibtime = system_clock::to_time_t(system_clock::now());
                    Rcout << ctime(&gibtime) << endl;
                    Rcout << "df101 ";//prints the log-likelihoods
                    for (int ij=0;ij<reqrdnum;ij++){
                        Rcout << Ll[ij] << " ";
                    }
                    Rcout << " " << endl;
                    Rcout << "df102 ";//prints the first derivatives
                    for (int ij=0;ij<reqrdnum;ij++){
                        Rcout << Lld[ij] << " ";
                    }
                    Rcout << " " << endl;
                    Rcout << "df103 ";//prints the second derivatives
                    for (int ij=0;ij<reqrdnum;ij++){
                        Rcout << Lldd[ij*reqrdnum+ij] << " ";
                    }
                    Lld_worst=0;
                    for (int ij=0;ij<reqrdnum;ij++){//locates highest magnitude derivative
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
                    for (int ij=0;ij<reqrdnum;ij++){//prints the newton step value for zero derivative
                        Rcout << Lld[ij]/Lldd[ij*reqrdnum+ij] << " ";
                    }
                    Rcout << " " << endl;
                    Rcout << "df106 ";
                    for (int ij=0;ij<reqrdnum;ij++){//prints the newton step value for zero log-likelihood
                        Rcout << Ll[ij]/Lld[ij] << " ";
                    }
                    Rcout << " " << endl;
                    Rcout << "df107 " << double_step << " " << abs_max << " " << dose_abs_max << " " << Ll_comp[0] << " " << Ll_comp[1] << endl;//prints several convergence terms
                }
                if (change_all){
                    if (Ll[ind0] <= Ll_abs_best){//takes a half-step if needed
                        #pragma omp parallel for num_threads(nthreads)
                        for (int ijk=0;ijk<totalnum;ijk++){
                            dbeta[ijk] = dbeta[ijk] * 0.5;//
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
                }
                #pragma omp parallel for num_threads(nthreads)
                for (int ijk=0;ijk<totalnum;ijk++){//totalnum*(totalnum+1)/2
                    beta_0[ijk] = beta_c[ijk];
                }
            }
        }
        if (beta_best!=beta_c){//if the risk matrices aren't the optimal values, then they must be recalculated
            if (verbose){
                Rcout << "Changing back to best"<<endl;
            }
            // If it goes through every half step without improvement, then the maximum change needs to be decreased
            abs_max = abs_max*pow(0.5,halfmax);
            dose_abs_max = dose_abs_max*pow(0.5,halfmax);
            iter_check = 1;
            //
            beta_p = beta_best;//
            beta_a = beta_best;//
            beta_c = beta_best;//
            Dose = MatrixXd::Zero(df0.rows(),term_tot);
            nonDose = MatrixXd::Constant(df0.rows(),term_tot,0.0);
            nonDose_LIN = MatrixXd::Constant(df0.rows(),term_tot,0.0);
            nonDose_PLIN = MatrixXd::Constant(df0.rows(),term_tot,1.0);
            nonDose_LOGLIN = MatrixXd::Constant(df0.rows(),term_tot,1.0);
            for (int ij=0;ij<totalnum;ij++){
                beta_0[ij] = beta_best[ij];
            }
            Make_subterms( totalnum, Term_n, tform, dfc, fir, T0, Td0, Tdd0, Dose, nonDose, TTerm,  nonDose_LIN, nonDose_PLIN, nonDose_LOGLIN,beta_0, df0,dint,dslp,nthreads, debugging,KeepConstant);
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
                for (int ijk=0;ijk<reqrdnum;ijk++){
                    Rcout << Td0.col(ijk).sum() << " ";
                }
                Rcout << " " << endl;
                Rcout << "second derivs checked ";
                for (int ijk=0;ijk<reqrdnum;ijk++){
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
            RdR = MatrixXd::Zero(df0.rows(), reqrdnum);
            RddR = MatrixXd::Zero(df0.rows(), reqrdnum*(reqrdnum+1)/2);
            //
            //
            Make_Risks_Weighted(modelform, tform, Term_n, totalnum, fir, s_weights, T0, Td0, Tdd0, Te, R, Rd, Rdd, Dose, nonDose, TTerm, nonDose_LIN, nonDose_PLIN, nonDose_LOGLIN, RdR, RddR, nthreads, debugging,KeepConstant);
            R = (R.array().isFinite()).select(R,0);
            Rd = (Rd.array().isFinite()).select(Rd,0);
            Rdd = (Rdd.array().isFinite()).select(Rdd,0);
            RdR = (RdR.array().isFinite()).select(RdR,0);
            RddR = (RddR.array().isFinite()).select(RddR,0);
            //
            temp_list = List::create(_["betas"]=wrap(beta_0),_["Status"]="FAILED"); //used as a dummy return value for code checking
            //risk_check_iter=0;
            if (verbose){
                Rcout << "risk checked ";
                for (int ijk=0;ijk<1;ijk++){
                    Rcout << R.col(0).sum() << " ";
                }
                Rcout << " " << endl;
                Rcout << "risk1 checked ";
                for (int ijk=0;ijk<reqrdnum;ijk++){
                    Rcout << Rd.col(ijk).sum() << " ";
                }
                Rcout << " " << endl;
                Rcout << "risk2 checked ";
                for (int ijk=0;ijk<reqrdnum;ijk++){
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
            Poisson_LogLik( nthreads, totalnum, PyrC, R, Rd, Rdd, RdR, RddR, Ll, Lld, Lldd, debugging,KeepConstant);
            dev_temp.col(0) = PyrC.col(0).array() * R.col(0).array();
            dev_temp.col(0) = PyrC.col(1).array() * dev_temp.col(0).array().pow(-1).array();
            dev_temp.col(0) = dev_temp.col(0).array().log().array();
            dev_temp.col(0) = PyrC.col(1).array() * dev_temp.col(0).array();
            dev_temp.col(1) = PyrC.col(1).array() - PyrC.col(0).array() * R.col(0).array();
            //
            dev_temp.col(0) = (dev_temp.col(0).array().isFinite()).select(dev_temp.col(0),0);
            //
            dev_temp.col(0) = dev_temp.col(0).array() - dev_temp.col(1).array();
            dev_temp.col(0) = (2 * dev_temp.col(0).array()).array();//.sqrt();
            dev_temp = (dev_temp.array().isFinite()).select(dev_temp,0);
            dev_temp.col(0) = (R.col(0).array()<0).select(0,dev_temp.col(0));
            dev = dev_temp.col(0).sum(); //deviation calculation is split into steps
            
        }
        Lld_worst=0;
        for (int ij=0;ij<reqrdnum;ij++){
            if (abs(Lld[ij]) > Lld_worst){
                Lld_worst = abs(Lld[ij]);
            }
        }
        if (iteration > reqrdnum){//Sets the minimum number of iterations
            if ((iteration % (reqrdnum))||(iter_check==1)){//Checks every set number of iterations
                iter_check=0;
                if (Lld_worst < deriv_epsilon){//ends if the derivatives are low enough
                    iter_stop = 1;
                    convgd = TRUE;
                }
                Ll_comp[1]=Ll[0];
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
            for (int ij=0;ij<reqrdnum;ij++){
                Rcout << Ll[ij] << " ";
            }
            Rcout << " " << endl;
            Rcout << "df102 ";
            for (int ij=0;ij<reqrdnum;ij++){
                Rcout << Lld[ij] << " ";
            }
            Rcout << " " << endl;
            Rcout << "df103 ";
            for (int ij=0;ij<reqrdnum;ij++){
                Rcout << Lldd[ij*reqrdnum+ij] << " ";
            }
            Rcout << " " << endl;
            Rcout << "df104 ";
            for (int ij=0;ij<totalnum;ij++){
                Rcout << beta_c[ij] << " ";
            }
            Rcout << " " << endl;
            Rcout << "Checking Best Deviance " << dev << endl;
            Rcout << "Finshed iteration" << endl;
            Rcout << "df105 ";
            for (int ij=0;ij<reqrdnum;ij++){
                Rcout << Lld[ij]/Lldd[ij*reqrdnum+ij] << " ";
            }
            Rcout << " " << endl;
            Rcout << "df106 ";
            for (int ij=0;ij<reqrdnum;ij++){
                Rcout << Ll[ij]/Lld[ij] << " ";
            }
            Rcout << " " << endl;
            Rcout << "df107 " << double_step << " " << abs_max << " " << dose_abs_max << " " << Ll_comp[0] << " " << Ll_comp[1] << endl;
        }
    }
    // -----------------------------------------------
    // Performing Full Calculation to get full second derivative matrix
    // -----------------------------------------------
    fill(Ll.begin(), Ll.end(), 0.0);
    fill(Lld.begin(), Lld.end(), 0.0);
    fill(Lldd.begin(), Lldd.end(), 0.0);
    Poisson_LogLik( nthreads, totalnum, PyrC, R, Rd, Rdd, RdR, RddR, Ll, Lld, Lldd, debugging,KeepConstant);
    dev_temp.col(0) = PyrC.col(0).array() * R.col(0).array();
    dev_temp.col(0) = PyrC.col(1).array() * dev_temp.col(0).array().pow(-1).array();
    dev_temp.col(0) = dev_temp.col(0).array().log().array();
    dev_temp.col(0) = PyrC.col(1).array() * dev_temp.col(0).array();
    dev_temp.col(1) = PyrC.col(1).array() - PyrC.col(0).array() * R.col(0).array();
    //
    dev_temp.col(0) = (dev_temp.col(0).array().isFinite()).select(dev_temp.col(0),0);
    //
    dev_temp.col(0) = dev_temp.col(0).array() - dev_temp.col(1).array();
    dev_temp.col(0) = (2 * dev_temp.col(0).array()).array();//.sqrt();
    dev_temp = (dev_temp.array().isFinite()).select(dev_temp,0);
    dev_temp.col(0) = (R.col(0).array()<0).select(0,dev_temp.col(0));
    dev = dev_temp.col(0).sum(); //deviation calculation is split into steps
    if (verbose){
        end_point = system_clock::now();
        ending = time_point_cast<microseconds>(end_point).time_since_epoch().count();
        Rcout <<"df100 "<<(ending-start)<<" "<<halves<<" "<<iteration<<" "<<ind0<<",Recalc"<<endl;
        gibtime = system_clock::to_time_t(system_clock::now());
        Rcout << ctime(&gibtime) << endl;
        Rcout << "df101 ";
        for (int ij=0;ij<reqrdnum;ij++){
            Rcout << Ll[ij] << " ";
        }
        Rcout << " " << endl;
        Rcout << "df102 ";
        for (int ij=0;ij<reqrdnum;ij++){
            Rcout << Lld[ij] << " ";
        }
        Rcout << " " << endl;
        Rcout << "df103 ";
        for (int ij=0;ij<reqrdnum;ij++){
            Rcout << Lldd[ij*reqrdnum+ij] << " ";
        }
        Rcout << " " << endl;
        Rcout << "df104 ";
        for (int ij=0;ij<totalnum;ij++){
            Rcout << beta_c[ij] << " ";
        }
        Rcout << " " << endl;
        Rcout << "df105 ";
        for (int ij=0;ij<reqrdnum;ij++){//prints the newton step value for zero derivative
            Rcout << Lld[ij]/Lldd[ij*reqrdnum+ij] << " ";
        }
        Rcout << " " << endl;
        Rcout << "df106 ";
        for (int ij=0;ij<reqrdnum;ij++){//prints the newton step value for zero log-likelihood
            Rcout << Ll[ij]/Lld[ij] << " ";
        }
        Rcout << " " << endl;
        Rcout << "df107 " << double_step << " " << abs_max << " " << dose_abs_max << " " << Ll_comp[0] << " " << Ll_comp[1] << endl;//prints several convergence terms
        Rcout << "Checking Best Deviance " << dev << endl;
        Rcout << "Finshed iteration" << endl;
    }
    //
    if ((Ll_abs_best>0)||(Ll_abs_best < Ll[ind0])){
        Ll_abs_best = Ll[ind0];
        beta_abs_best = beta_c;
    }
    //
    if (verbose){
        Rcout << "df501 " << Ll_abs_best << endl;
        Rcout << "df504 ";//prints parameter values
        for (int ij=0;ij<totalnum;ij++){
            Rcout << beta_abs_best[ij] << " ";
        }
        Rcout << " " << endl;
    }
    // Account for the strata parameters
    // Changes the parameter back into the original form
    List para_list = List::create(_["Term_n"]=Term_n,_["tforms"]=tform);
    List control_list = List::create(_["Iteration"]=iteration, _["Maximum Step"]=abs_max, _["Derivative Limiting"]=Lld_worst);
    //
    int kept_covs = totalnum - sum(KeepConstant);
    NumericVector Lldd_vec(kept_covs * kept_covs);
    #pragma omp parallel for schedule(dynamic) num_threads(nthreads)
    for (int ijk=0;ijk<totalnum*(totalnum+1)/2;ijk++){
        int ij = 0;
        int jk = ijk;
        int pij_ind=-100;
        int pjk_ind=-100;
        while (jk>ij){
            ij++;
            jk-=ij;
        }
        if (KeepConstant[ij]==0){
            pij_ind = ij - sum(head(KeepConstant,ij));
            if (KeepConstant[jk]==0){
                pjk_ind = jk - sum(head(KeepConstant,jk));
                Lldd_vec[pij_ind * kept_covs + pjk_ind]=Lldd[ij*totalnum+jk];
            }
        }
    }
    for (int ijk=0;ijk<kept_covs*(kept_covs+1)/2;ijk++){
        int ij = 0;
        int jk = ijk;
        while (jk>ij){
            ij++;
            jk-=ij;
        }
        Lldd_vec[jk * kept_covs + ij]=Lldd_vec[ij * kept_covs + jk];
    }
    Lldd_vec.attr("dim") = Dimension(kept_covs, kept_covs);
    const Map<MatrixXd> Lldd_mat(as<Map<MatrixXd> >(Lldd_vec));
    //
    MatrixXd Lldd_inv = -1 * Lldd_mat.inverse().matrix(); //uses inverse information matrix to calculate the standard deviation
    VectorXd stdev = VectorXd::Zero(totalnum);
    for (int ij=0;ij<totalnum;ij++){
        if (KeepConstant[ij]==0){
            int pij_ind = ij - sum(head(KeepConstant,ij));
            stdev(ij) = sqrt(Lldd_inv(pij_ind,pij_ind));
        }
    }
    //
    List res_list = List::create(_["LogLik"]=wrap(Ll[0]),_["First_Der"]=wrap(Lld),_["Second_Der"]=Lldd_vec,_["beta_0"]=wrap(beta_0) ,_["Standard_Deviation"]=wrap(stdev) ,_["AIC"]=2*(totalnum-accumulate(KeepConstant.begin(),KeepConstant.end(), 0.0))+dev,_["Deviation"]=dev,_["Parameter_Lists"]=para_list,_["Control_List"]=control_list,_["Converged"]=convgd);
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
//' @param     double_step controls the step calculation, 0 for independent changes, 1 for solving b=Ax with complete matrices
//' @param     change_all    boolean if every parameter is being updated
//' @param     verbose    verbosity boolean
//' @param     debugging    debugging boolean
//' @param     KeepConstant    vector identifying constant parameters
//' @param     term_tot    total number of terms
//' @param     debug_checks    string vector of functions to test
//' @param     ties_method    ties method
//' @param     nthreads number of threads to use
//'
//' @return NULL
// [[Rcpp::export]]
void Stress_Run( IntegerVector Term_n, StringVector tform, NumericVector a_n,NumericMatrix x_all,IntegerVector dfc,int fir, int der_iden,string modelform, double lr, int maxiter, int halfmax, double epsilon, double dbeta_cap, double abs_max,double dose_abs_max, double deriv_epsilon, NumericMatrix df_groups, NumericVector tu, int double_step ,bool change_all, bool verbose, bool debugging, IntegerVector KeepConstant, int term_tot, StringVector debug_checks, string ties_method, int nthreads){
;
    //
    // Runs through a single calculation with some functions printing additional information printed
    //
    if (verbose){
        Rcout << "START_STRESS" << endl;
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
    int reqrdnum = totalnum - sum(KeepConstant);
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
    // int nthreads = Eigen::nbThreads()-1; //stores how many threads are allocated    //
    // Lld_worst: stores the highest magnitude log-likelihood derivative
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
    MatrixXd Td0 = MatrixXd::Zero(df0.rows(), reqrdnum); //preallocates matrix for Term derivative columns
    MatrixXd Tdd0 = MatrixXd::Zero(df0.rows(), reqrdnum*(reqrdnum+1)/2); //preallocates matrix for Term second derivative columns
    //
    MatrixXd Te = MatrixXd::Zero(df0.rows(), 1); //preallocates matrix for column terms used for temporary storage
    MatrixXd R = MatrixXd::Zero(df0.rows(), 1); //preallocates matrix for Risks
    ColXd Rd = MatrixXd::Zero(df0.rows(), reqrdnum); //preallocates matrix for Risk derivatives
    ColXd Rdd = MatrixXd::Zero(df0.rows(), reqrdnum*(reqrdnum+1)/2); //preallocates matrix for Risk second derivatives
    //
    MatrixXd Dose = MatrixXd::Zero(df0.rows(),term_tot); //Matrix of the total dose term values
    MatrixXd nonDose = MatrixXd::Constant(df0.rows(),term_tot,1.0); //Matrix of the total non-dose term values
    MatrixXd nonDose_LIN = MatrixXd::Constant(df0.rows(),term_tot,0.0); //matrix of Linear subterm values
    MatrixXd nonDose_PLIN = MatrixXd::Constant(df0.rows(),term_tot,1.0); //matrix of Loglinear subterm values
    MatrixXd nonDose_LOGLIN = MatrixXd::Constant(df0.rows(),term_tot,1.0); //matrix of Product linear subterm values
    MatrixXd TTerm = MatrixXd::Zero(Dose.rows(),Dose.cols()); //matrix of term values
    double dint = dose_abs_max; //The amount of change used to calculate derivatives in threshold paramters    //
    double dslp = abs_max;
    //
    // Calculates the subterm and term values
    if (verbose){
        Rcout << "starting subterms " << term_tot << endl;
    }
    if (Debug_It[0]){
            Make_subterms( totalnum, Term_n, tform, dfc, fir, T0, Td0, Tdd0, Dose, nonDose, TTerm,  nonDose_LIN, nonDose_PLIN, nonDose_LOGLIN ,beta_0, df0,dint,dslp,nthreads, TRUE,KeepConstant);
    } else {
            Make_subterms( totalnum, Term_n, tform, dfc, fir, T0, Td0, Tdd0, Dose, nonDose, TTerm,  nonDose_LIN, nonDose_PLIN, nonDose_LOGLIN ,beta_0, df0,dint,dslp,nthreads, FALSE,KeepConstant);
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
        for (int ijk=0;ijk<reqrdnum;ijk++){
            Rcout << Td0.col(ijk).sum() << " ";
        }
        Rcout << " " << endl;
        Rcout << "second derivs checked ";
        for (int ijk=0;ijk<reqrdnum;ijk++){
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
    ColXd RdR = MatrixXd::Zero(df0.rows(), reqrdnum); //preallocates matrix for Risk to derivative ratios
    ColXd RddR = MatrixXd::Zero(df0.rows(), reqrdnum*(reqrdnum+1)/2); //preallocates matrix for Risk to second derivative ratios
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
        Make_Risks(modelform, tform, Term_n, totalnum, fir, T0, Td0, Tdd0, Te, R, Rd, Rdd, Dose, nonDose, TTerm, nonDose_LIN, nonDose_PLIN, nonDose_LOGLIN, RdR, RddR, nthreads, TRUE,KeepConstant);
    } else {
        Make_Risks(modelform, tform, Term_n, totalnum, fir, T0, Td0, Tdd0, Te, R, Rd, Rdd, Dose, nonDose, TTerm, nonDose_LIN, nonDose_PLIN, nonDose_LOGLIN, RdR, RddR, nthreads, FALSE,KeepConstant);
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
        for (int ijk=0;ijk<1;ijk++){
            Rcout << R.col(0).sum() << " ";
        }
        Rcout << " " << endl;
    //    return temp_list;
        Rcout << "risk1 checked ";
        for (int ijk=0;ijk<reqrdnum;ijk++){
            Rcout << Rd.col(ijk).sum() << " ";
        }
        Rcout << " " << endl;
        Rcout << "risk2 checked ";
        for (int ijk=0;ijk<reqrdnum;ijk++){
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
    MatrixXd Rls2 =MatrixXd::Zero(ntime, reqrdnum); //Many are repeated due to the same risk groups and derivatives being used at mulitple points
    MatrixXd Rls3 =MatrixXd::Zero(ntime, reqrdnum*(reqrdnum+1)/2); //Sum and its derivatives are precomputed
    MatrixXd Lls1 =MatrixXd::Zero(ntime, 1); //The log-likelihood calculation has a Right and Left sum used
    MatrixXd Lls2 =MatrixXd::Zero(ntime, reqrdnum);
    MatrixXd Lls3 =MatrixXd::Zero(ntime, reqrdnum*(reqrdnum+1)/2);
    //The log-likelihood is calculated in parallel over the risk groups
    vector<double> Ll(reqrdnum,0.0); //Log-likelihood values
    vector<double> Lld(reqrdnum,0.0); //Log-likelihood derivative values
    vector<double> Lldd(pow(reqrdnum,2),0.0);//The second derivative matrix has room for every combination, but only the lower triangle is calculated initially
    //
    //
    // Calculates the side sum terms used
    if (Debug_It[3]){
        Calculate_Sides( RiskFail, RiskGroup, totalnum, ntime, R, Rd, Rdd, Rls1, Rls2, Rls3, Lls1, Lls2, Lls3,nthreads, TRUE,KeepConstant);
    } else {
        Calculate_Sides( RiskFail, RiskGroup, totalnum, ntime, R, Rd, Rdd, Rls1, Rls2, Rls3, Lls1, Lls2, Lls3,nthreads, FALSE,KeepConstant);
    }
    //
    //
    if (verbose){
        Rcout << "riskr checked ";
        for (int ijk=0;ijk<reqrdnum;ijk++){
            Rcout << Rls1.col(0).sum() << " ";
        }
        Rcout << " " << endl;
        Rcout << "risk1r checked ";
        for (int ijk=0;ijk<reqrdnum;ijk++){
            Rcout << Rls2.col(ijk).sum() << " ";
        }
        Rcout << " " << endl;
        Rcout << "risk2r checked ";
        for (int ijk=0;ijk<reqrdnum;ijk++){
            Rcout << Rls3.col(ijk*(ijk+1)/2+ijk).sum() << " ";
        }
        Rcout << " " << endl;
        //
        Rcout << "riskl checked ";
        for (int ijk=0;ijk<reqrdnum;ijk++){
            Rcout << Lls1.col(0).sum() << " ";
        }
        Rcout << " " << endl;
        Rcout << "risk1l checked ";
        for (int ijk=0;ijk<reqrdnum;ijk++){
            Rcout << Lls2.col(ijk).sum() << " ";
        }
        Rcout << " " << endl;
        Rcout << "risk2l checked ";
        for (int ijk=0;ijk<reqrdnum;ijk++){
            Rcout << Lls3.col(ijk*(ijk+1)/2+ijk).sum() << " ";
        }
        Rcout << " " << endl;
    }
    //
    //
    // Calculates log-likelihood
    if (Debug_It[4]){
        Calc_LogLik( nthreads, RiskFail, RiskGroup, totalnum, ntime, R, Rd, Rdd,RdR,RddR, Rls1, Rls2, Rls3, Lls1, Lls2, Lls3, Ll, Lld, Lldd, TRUE, ties_method,KeepConstant);
    } else {
        Calc_LogLik( nthreads, RiskFail, RiskGroup, totalnum, ntime, R, Rd, Rdd,RdR,RddR, Rls1, Rls2, Rls3, Lls1, Lls2, Lls3, Ll, Lld, Lldd, FALSE, ties_method,KeepConstant);
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
        for (int ij=0;ij<reqrdnum;ij++){
            Rcout << Ll[ij] << " ";
        }
        Rcout << " " << endl;
        Rcout << "df102 ";//prints the first derivatives
        for (int ij=0;ij<reqrdnum;ij++){
            Rcout << Lld[ij] << " ";
        }
        Rcout << " " << endl;
        Rcout << "df103 ";//prints the second derivatives
        for (int ij=0;ij<reqrdnum;ij++){
            Rcout << Lldd[ij*reqrdnum+ij] << " ";
        }
        Lld_worst=0;
        for (int ij=0;ij<reqrdnum;ij++){//locates highest magnitude derivative
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
        for (int ij=0;ij<reqrdnum;ij++){//prints the newton step value for zero derivative
            Rcout << Lld[ij]/Lldd[ij*reqrdnum+ij] << " ";
        }
        Rcout << " " << endl;
        Rcout << "df106 ";
        for (int ij=0;ij<reqrdnum;ij++){//prints the newton step value for zero log-likelihood
            Rcout << Ll[ij]/Lld[ij] << " ";
        }
        Rcout << " " << endl;
        Rcout << "df107 " << double_step << " " << abs_max << " " << dose_abs_max << " " << Ll_comp[0] << " " << Ll_comp[1] << endl;//prints several convergence terms
    }
    //
    vector<double> dbeta(totalnum,0.0);
    //
    // --------------------------
    // always starts from initial guess
    // --------------------------
    vector<double> beta_p(totalnum,0.0);
    vector<double> beta_c(totalnum,0.0);
    vector<double> beta_a(totalnum,0.0);
    vector<double> beta_best(totalnum,0.0);
    VectorXd::Map(&beta_p[0], beta_0.size()) = beta_0;// stores previous parameters
    VectorXd::Map(&beta_c[0], beta_0.size()) = beta_0;// stores current parameters
    VectorXd::Map(&beta_a[0], beta_0.size()) = beta_0;// stores a refrence value for parameters
    VectorXd::Map(&beta_best[0], beta_0.size()) = beta_0;// stores the best parameters
    double halves = 0; //number of half-steps taken
    int ind0 = fir; //used for validations
    int iteration=0; //iteration number
    //
    if (sum(KeepConstant)==totalnum){
        iteration=maxiter;
        if (verbose){
            Rcout << "No parameters to change" << endl;
        }
    }
    double Ll_abs_best = 10;
    vector<double> beta_abs_best(totalnum,0.0);
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
            Calc_Change( double_step, nthreads, totalnum, fir, der_iden, dbeta_cap, dose_abs_max, lr, abs_max, Ll, Lld, Lldd, dbeta, change_all, tform, dint,dslp, KeepConstant, TRUE);
        } else {
            Calc_Change( double_step, nthreads, totalnum, fir, der_iden, dbeta_cap, dose_abs_max, lr, abs_max, Ll, Lld, Lldd, dbeta, change_all, tform, dint,dslp, KeepConstant, FALSE);
        }
        Intercept_Bound(nthreads, totalnum, beta_0, dbeta, dfc, df0, KeepConstant, debugging, tform);
        if (verbose){
            Rcout << "Starting Halves"<<endl;//prints the final changes for validation
        }
        //
        //
        if ((Ll_abs_best>0)||(Ll_abs_best < Ll[ind0])){
            Ll_abs_best = Ll[ind0];
            beta_abs_best = beta_c;
        }
        //
        if (verbose){
            Rcout << "df501 " << Ll_abs_best << endl;
            Rcout << "df504 ";//prints parameter values
            for (int ij=0;ij<totalnum;ij++){
                Rcout << beta_abs_best[ij] << " ";
            }
            Rcout << " " << endl;
        }
        halves=0;
        while ((Ll[ind0] <= Ll_abs_best)&&(halves<halfmax)){ //repeats until half-steps maxed or an improvement
            halves++;
            Dose = MatrixXd::Zero(df0.rows(),term_tot); //Refreshes the matrices used
            nonDose = MatrixXd::Constant(df0.rows(),term_tot,0.0);
            nonDose_LIN = MatrixXd::Constant(df0.rows(),term_tot,0.0);
            nonDose_PLIN = MatrixXd::Constant(df0.rows(),term_tot,1.0);
            nonDose_LOGLIN = MatrixXd::Constant(df0.rows(),term_tot,1.0);
            T0 = MatrixXd::Zero(df0.rows(), totalnum); 
            Td0 = MatrixXd::Zero(df0.rows(), reqrdnum); 
            Tdd0 = MatrixXd::Zero(df0.rows(), reqrdnum*(reqrdnum+1)/2); 
            for (int ij=0;ij<totalnum;ij++){
                beta_0[ij] = beta_a[ij] + dbeta[ij];//applies the parameter changes
                beta_c[ij] = beta_0[ij];
            }
            // ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------//
            // The same subterm, risk, sides, and log-likelihood calculations are performed every half-step and iteration
            // ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------//
            if (Debug_It[4]){
                Make_subterms( totalnum, Term_n, tform, dfc, fir, T0, Td0, Tdd0, Dose,  nonDose, TTerm,  nonDose_LIN, nonDose_PLIN, nonDose_LOGLIN ,beta_0, df0,dint,dslp,nthreads, TRUE,KeepConstant);
            } else {
                Make_subterms( totalnum, Term_n, tform, dfc, fir, T0, Td0, Tdd0, Dose,  nonDose, TTerm,  nonDose_LIN, nonDose_PLIN, nonDose_LOGLIN ,beta_0, df0,dint,dslp,nthreads, FALSE,KeepConstant);
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
                for (int ijk=0;ijk<reqrdnum;ijk++){
                    Rcout << Td0.col(ijk).sum() << " ";
                }
                Rcout << " " << endl;
                Rcout << "second derivs checked ";
                for (int ijk=0;ijk<reqrdnum;ijk++){
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
            RdR = MatrixXd::Zero(df0.rows(), reqrdnum);
            RddR = MatrixXd::Zero(df0.rows(), reqrdnum*(reqrdnum+1)/2);
            //
            //
            if (Debug_It[6]){
                Make_Risks(modelform, tform, Term_n, totalnum, fir, T0, Td0, Tdd0, Te, R, Rd, Rdd, Dose, nonDose, TTerm, nonDose_LIN, nonDose_PLIN, nonDose_LOGLIN, RdR, RddR, nthreads, TRUE,KeepConstant);
            } else {
                Make_Risks(modelform, tform, Term_n, totalnum, fir, T0, Td0, Tdd0, Te, R, Rd, Rdd, Dose, nonDose, TTerm, nonDose_LIN, nonDose_PLIN, nonDose_LOGLIN, RdR, RddR, nthreads, FALSE,KeepConstant);
            }
            //
            RdR = (RdR.array().isFinite()).select(RdR,0);
            RddR = (RddR.array().isFinite()).select(RddR,0);
            //
            temp_list = List::create(_["betas"]=wrap(beta_0),_["Status"]="FAILED"); //used as a dummy return value for code checking
            if (verbose){
                Rcout << "risk checked ";
                for (int ijk=0;ijk<1;ijk++){
                    Rcout << R.col(0).sum() << " ";
                }
                Rcout << " " << endl;
                Rcout << "risk1 checked ";
                for (int ijk=0;ijk<reqrdnum;ijk++){
                    Rcout << Rd.col(ijk).sum() << " ";
                }
                Rcout << " " << endl;
                Rcout << "risk2 checked ";
                for (int ijk=0;ijk<reqrdnum;ijk++){
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
            Rls2 =MatrixXd::Zero(ntime, reqrdnum); //Many are repeated due to the same risk groups and derivatives being used at mulitple points
            Rls3 =MatrixXd::Zero(ntime, reqrdnum*(reqrdnum+1)/2);
            Lls1 =MatrixXd::Zero(ntime, 1);
            Lls2 =MatrixXd::Zero(ntime, reqrdnum);
            Lls3 =MatrixXd::Zero(ntime, reqrdnum*(reqrdnum+1)/2);
            if (Debug_It[7]){
                Calculate_Sides( RiskFail, RiskGroup, totalnum, ntime, R, Rd, Rdd, Rls1, Rls2, Rls3, Lls1, Lls2, Lls3,nthreads, TRUE,KeepConstant);
            } else {
                Calculate_Sides( RiskFail, RiskGroup, totalnum, ntime, R, Rd, Rdd, Rls1, Rls2, Rls3, Lls1, Lls2, Lls3,nthreads, FALSE,KeepConstant);
            }
            //
            if (verbose){
                Rcout << "riskr checked ";
                for (int ijk=0;ijk<1;ijk++){
                    Rcout << Rls1.col(0).sum() << " ";
                }
                Rcout << " " << endl;
                Rcout << "risk1r checked ";
                for (int ijk=0;ijk<reqrdnum;ijk++){
                    Rcout << Rls2.col(ijk).sum() << " ";
                }
                Rcout << " " << endl;
                Rcout << "risk2r checked ";
                for (int ijk=0;ijk<reqrdnum;ijk++){
                    Rcout << Rls3.col(ijk*(ijk+1)/2+ijk).sum() << " ";
                }
                Rcout << " " << endl;
                //
                //
                Rcout << "riskl checked ";
                for (int ijk=0;ijk<1;ijk++){
                    Rcout << Lls1.col(0).sum() << " ";
                }
                Rcout << " " << endl;
                Rcout << "risk1l checked ";
                for (int ijk=0;ijk<reqrdnum;ijk++){
                    Rcout << Lls2.col(ijk).sum() << " ";
                }
                Rcout << " " << endl;
                Rcout << "risk2l checked ";
                for (int ijk=0;ijk<reqrdnum;ijk++){
                    Rcout << Lls3.col(ijk*(ijk+1)/2+ijk).sum() << " ";
                }
                Rcout << " " << endl;
            }
            //
            if (Debug_It[8]){
                Calc_LogLik( nthreads, RiskFail, RiskGroup, totalnum, ntime, R, Rd, Rdd,RdR,RddR, Rls1, Rls2, Rls3, Lls1, Lls2, Lls3, Ll, Lld, Lldd, TRUE, ties_method,KeepConstant);
            } else {
                Calc_LogLik( nthreads, RiskFail, RiskGroup, totalnum, ntime, R, Rd, Rdd,RdR,RddR, Rls1, Rls2, Rls3, Lls1, Lls2, Lls3, Ll, Lld, Lldd, FALSE, ties_method,KeepConstant);
            }
            if (verbose){
                end_point = system_clock::now();
                ending = time_point_cast<microseconds>(end_point).time_since_epoch().count();
                Rcout <<"df100 "<<(ending-start)<<" "<<0<<" "<<0<<" "<<0<<",Half Calc"<<endl;//prints the time
                gibtime = system_clock::to_time_t(system_clock::now());
                Rcout << ctime(&gibtime) << endl;
                Rcout << "df101 ";//prints the log-likelihoods
                for (int ij=0;ij<reqrdnum;ij++){
                    Rcout << Ll[ij] << " ";
                }
                Rcout << " " << endl;
                Rcout << "df102 ";//prints the first derivatives
                for (int ij=0;ij<reqrdnum;ij++){
                    Rcout << Lld[ij] << " ";
                }
                Rcout << " " << endl;
                Rcout << "df103 ";//prints the second derivatives
                for (int ij=0;ij<reqrdnum;ij++){
                    Rcout << Lldd[ij*reqrdnum+ij] << " ";
                }
                Lld_worst=0;
                for (int ij=0;ij<reqrdnum;ij++){//locates highest magnitude derivative
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
                for (int ij=0;ij<reqrdnum;ij++){//prints the newton step value for zero derivative
                    Rcout << Lld[ij]/Lldd[ij*reqrdnum+ij] << " ";
                }
                Rcout << " " << endl;
                Rcout << "df106 ";
                for (int ij=0;ij<reqrdnum;ij++){//prints the newton step value for zero log-likelihood
                    Rcout << Ll[ij]/Lld[ij] << " ";
                }
                Rcout << " " << endl;
                Rcout << "df107 " << double_step << " " << abs_max << " " << dose_abs_max << " " << Ll_comp[0] << " " << Ll_comp[1] << endl;//prints several convergence terms
            }
            if (change_all){
                if (Ll[ind0] <= Ll_abs_best){//takes a half-step if needed
                    #pragma omp parallel for num_threads(nthreads)
                    for (int ijk=0;ijk<totalnum;ijk++){
                        dbeta[ijk] = dbeta[ijk] * 0.5; //
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
            for (int ij=0;ij<reqrdnum;ij++){
                Rcout << Ll[ij] << " ";
            }
            Rcout << " " << endl;
            Rcout << "df102 ";
            for (int ij=0;ij<reqrdnum;ij++){
                Rcout << Lld[ij] << " ";
            }
            Rcout << " " << endl;
            Rcout << "df103 ";
            for (int ij=0;ij<reqrdnum;ij++){
                Rcout << Lldd[ij*reqrdnum+ij] << " ";
            }
            Rcout << " " << endl;
            Rcout << "df104 ";
            for (int ij=0;ij<totalnum;ij++){
                Rcout << beta_c[ij] << " ";
            }
            Rcout << " " << endl;
            Rcout << "df105 ";
            for (int ij=0;ij<reqrdnum;ij++){
                Rcout << Lld[ij]/Lldd[ij*reqrdnum+ij] << " ";
            }
            Rcout << " " << endl;
            Rcout << "df106 ";
            for (int ij=0;ij<reqrdnum;ij++){
                Rcout << Ll[ij]/Lld[ij] << " ";
            }
            Rcout << " " << endl;
            Rcout << "df107 " << double_step << " " << abs_max << " " << dose_abs_max << " " << Ll_comp[0] << " " << Ll_comp[1] << endl;
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
//' @param     df_groups    status/time matrix
//' @param     tu    vector of event times
//' @param     verbose    verbose boolean
//' @param     ties_method    tied event method
//' @param     nthreads number of threads to use
//'
//' @return List of results: Log-likelihood of optimum, AIC
// [[Rcpp::export]]
List LogLik_Cox_PH_null( NumericMatrix df_groups, NumericVector tu, bool verbose, string ties_method, int nthreads){
    ;
    //
    // null model value calculation
    // --------------------------------------------------------------------------- //
    // The same steps are taken for the non-null case, with the exception of derivative calculations
    // --------------------------------------------------------------------------- //
    //
    if (verbose){
        Rcout << "START_COX_NULL" << endl;
    }
    int totalnum=1;
    int reqrdnum = 1;
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
    // int nthreads = Eigen::nbThreads()-1; //stores how many threads are allocated
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
    int ntime = tu.size();
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
        for (int ijk=0;ijk<reqrdnum;ijk++){
            Rcout << Rls1.col(0).sum() << " ";
        }
        Rcout << " " << endl;
        //
        Rcout << "riskl checked ";
        for (int ijk=0;ijk<reqrdnum;ijk++){
            Rcout << Lls1.col(0).sum() << " ";
        }
        Rcout << " " << endl;
    }
    //
    //
    Calc_Null_LogLik( nthreads, RiskFail, RiskGroup, ntime, R, Rls1, Lls1, Ll, ties_method);
    //
    List res_list = List::create(_["LogLik"]=wrap(Ll[0]),_["AIC"]=-2*Ll[0]);
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
//' @param     nthreads number of threads to use
//'
//' @return List of results: Risk at the reference
// [[Rcpp::export]]
NumericVector RISK_SUBSET(IntegerVector Term_n, StringVector tform, NumericVector a_n,NumericMatrix x_all,IntegerVector dfc,int fir,string modelform, bool verbose, bool debugging, int term_tot, int nthreads){
    ;
    //
    // Calculates the terms and risks for a reference vector
    //
    using namespace std::chrono;
    if (verbose){
        Rcout << "START_RISK" << endl;
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
    IntegerVector KeepConstant(totalnum);
    int reqrdnum = totalnum - sum(KeepConstant);
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
    // int nthreads = Eigen::nbThreads()-1; //stores how many threads are allocated
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
    MatrixXd Td0 = MatrixXd::Zero(df0.rows(), reqrdnum); //preallocates matrix for Term derivative columns
    MatrixXd Tdd0 = MatrixXd::Zero(df0.rows(), reqrdnum*(reqrdnum+1)/2); //preallocates matrix for Term second derivative columns
    //
    MatrixXd Te = MatrixXd::Zero(df0.rows(), 1); //preallocates matrix for column terms used for temporary storage
    MatrixXd R = MatrixXd::Zero(df0.rows(), 1); //preallocates matrix for Risks
    ColXd Rd = MatrixXd::Zero(df0.rows(), reqrdnum); //preallocates matrix for Risk derivatives
    ColXd Rdd = MatrixXd::Zero(df0.rows(), reqrdnum*(reqrdnum+1)/2); //preallocates matrix for Risk second derivatives
    //
    MatrixXd Dose = MatrixXd::Zero(df0.rows(),term_tot); //Matrix of the total dose term values
    MatrixXd nonDose = MatrixXd::Constant(df0.rows(),term_tot,1.0); //Matrix of the total non-dose term values
    MatrixXd nonDose_LIN = MatrixXd::Constant(df0.rows(),term_tot,0.0); //matrix of Linear subterm values
    MatrixXd nonDose_PLIN = MatrixXd::Constant(df0.rows(),term_tot,1.0); //matrix of Loglinear subterm values
    MatrixXd nonDose_LOGLIN = MatrixXd::Constant(df0.rows(),term_tot,1.0); //matrix of Product linear subterm values
    MatrixXd TTerm = MatrixXd::Zero(Dose.rows(),Dose.cols()); //matrix of term values
    //
    // Calculates terms
    //
    if (verbose){
        Rcout << "Starting Terms" << endl;
    }
    Make_subterms( totalnum, Term_n, tform, dfc, fir, T0, Td0, Tdd0, Dose, nonDose, TTerm,  nonDose_LIN, nonDose_PLIN, nonDose_LOGLIN ,beta_0, df0,1.0,1.0,nthreads, debugging,KeepConstant);
    ColXd RdR = MatrixXd::Zero(df0.rows(), reqrdnum); //preallocates matrix for Risk to derivative ratios
    ColXd RddR = MatrixXd::Zero(df0.rows(), reqrdnum*(reqrdnum+1)/2);
    //
    // Calculates risk
    //
    if (verbose){
        Rcout << "Starting Risks" << endl;
    }
    Make_Risks(modelform, tform, Term_n, totalnum, fir, T0, Td0, Tdd0, Te, R, Rd, Rdd, Dose, nonDose, TTerm, nonDose_LIN, nonDose_PLIN, nonDose_LOGLIN, RdR, RddR, nthreads, debugging,KeepConstant);
    //
    //
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
        for (int ijk=0;ijk<reqrdnum;ijk++){
            Rcout << Td0.col(ijk).sum() << " ";
        }
        Rcout << " " << endl;
        Rcout << "second derivs checked ";
        for (int ijk=0;ijk<reqrdnum;ijk++){
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

//' Primary Cox PH calculation
//' \code{LogLik_Cox_PH_Single} Calculates just the log-likelihood
//'
//' @param     Term_n    Term numbers
//' @param     tform    subterm types
//' @param     a_n    starting values
//' @param     x_all    covariate matrix
//' @param     dfc    covariate column numbers
//' @param     fir    first term number
//' @param     modelform    model string
//' @param     df_groups    matrix with time and event information
//' @param     tu    event times
//' @param     verbose    verbosity boolean
//' @param     debugging    debugging boolean
//' @param     KeepConstant    vector identifying constant parameters
//' @param     term_tot    total number of terms
//' @param     ties_method    ties method
//' @param     nthreads number of threads to use
//'
//' @return List of results: Log-likelihood, parameter list, AIC, model information
// [[Rcpp::export]]
List LogLik_Cox_PH_Single( IntegerVector Term_n, StringVector tform, NumericVector a_n,NumericMatrix x_all,IntegerVector dfc,int fir,string modelform, NumericMatrix df_groups, NumericVector tu, bool verbose, bool debugging, IntegerVector KeepConstant, int term_tot, string ties_method, int nthreads){
    ;
    //
    List temp_list = List::create(_["Status"]="FAILED"); //used as a dummy return value for code checking
    if (verbose){
        Rcout << "START_COX_SINGLE" << endl;
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
    const Map<MatrixXd> df0(as<Map<MatrixXd> >(x_all));
    int ntime = tu.size();
    //
    int totalnum = Term_n.size();
    int reqrdnum = totalnum - sum(KeepConstant);
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
    Rcout.precision(7); //forces higher precision numbers printed to terminal
    // int nthreads = Eigen::nbThreads()-1; //stores how many threads are allocated
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
    //
    MatrixXd Te = MatrixXd::Zero(df0.rows(), 1); //preallocates matrix for column terms used for temporary storage
    MatrixXd R = MatrixXd::Zero(df0.rows(), 1); //preallocates matrix for Risks
    //
    MatrixXd Dose = MatrixXd::Constant(df0.rows(),term_tot,0.0); //Matrix of the total dose term values
    MatrixXd nonDose = MatrixXd::Constant(df0.rows(),term_tot,1.0); //Matrix of the total non-dose term values
    MatrixXd nonDose_LIN = MatrixXd::Constant(df0.rows(),term_tot,0.0); //matrix of Linear subterm values
    MatrixXd nonDose_PLIN = MatrixXd::Constant(df0.rows(),term_tot,1.0); //matrix of Loglinear subterm values
    MatrixXd nonDose_LOGLIN = MatrixXd::Constant(df0.rows(),term_tot,1.0); //matrix of Product linear subterm values
    MatrixXd TTerm = MatrixXd::Zero(Dose.rows(),Dose.cols()); //matrix of term values
    //
    if (verbose){
        Rcout << "starting subterms " << term_tot << endl;
    }
    // Calculates the subterm and term values
    Make_subterms_Single( totalnum, Term_n, tform, dfc, fir, T0, Dose, nonDose, TTerm,  nonDose_LIN, nonDose_PLIN, nonDose_LOGLIN ,beta_0, df0,nthreads, debugging,KeepConstant);
    // ---------------------------------------------------------
    // Prints off a series of calculations to check at what point values are changing
    // ---------------------------------------------------------
    //
    //
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
    Make_Risks_Single(modelform, tform, Term_n, totalnum, fir, T0, Te, R, Dose, nonDose, TTerm, nonDose_LIN, nonDose_PLIN, nonDose_LOGLIN, nthreads, debugging,KeepConstant);
    //
    // Removes infinite values
    //
    if (R.minCoeff()<0){
        Rcout << "A non-positive risk was detected: " << R.minCoeff() << endl;
        Rcout << "final failing values ";
        for (int ijk=0;ijk<totalnum;ijk++){
            Rcout << beta_0[ijk] << " ";
        }
        Rcout << " " << endl;
        //
        Rcout << "final failing terms ";
        for (int ijk=0;ijk<totalnum;ijk++){
            Rcout << tform[ijk] << " ";
        }
        Rcout << " " << endl;
        temp_list = List::create(_["beta_0"]=wrap(beta_0) ,_["Deviation"]=R_NaN,_["Status"]="FAILED",_["LogLik"]=R_NaN);
        return temp_list;
    }
    //
    if (verbose){
        Rcout << "risk checked ";
        for (int ijk=0;ijk<1;ijk++){
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
    //
    MatrixXd Rls1 =MatrixXd::Zero(ntime, 1); //precomputes a series of sums used frequently in the log-liklihood calculations
    MatrixXd Lls1 =MatrixXd::Zero(ntime, 1); //The log-likelihood calculation has a Right and Left sum used
    //The log-likelihood is calculated in parallel over the risk groups
    vector<double> Ll(1,0.0); //Log-likelihood values
    //
    if (verbose){
        Rcout << "Made Risk Side Lists" << endl;
    }
    // Calculates the side sum terms used
    Calculate_Sides_Single( RiskFail, RiskGroup, totalnum, ntime, R, Rls1, Lls1,nthreads, debugging);
    //
    if (verbose){
        end_point = system_clock::now();
        ending = time_point_cast<microseconds>(end_point).time_since_epoch().count();
        Rcout <<"df100 "<<(ending-start)<<" "<<0<<" "<<0<<" "<<-1<<",Prep_Sides"<<endl;
        gibtime = system_clock::to_time_t(system_clock::now());
        Rcout << ctime(&gibtime) << endl;
    }
    //
    if (verbose){
        Rcout << "riskr checked ";
        for (int ijk=0;ijk<reqrdnum;ijk++){
            Rcout << Rls1.col(0).sum() << " ";
        }
        Rcout << " " << endl;
        //
        Rcout << "riskl checked ";
        for (int ijk=0;ijk<reqrdnum;ijk++){
            Rcout << Lls1.col(0).sum() << " ";
        }
        Rcout << " " << endl;
    }
    //
    //
    // Calculates log-likelihood
    Calc_LogLik_Single( nthreads, RiskFail, RiskGroup, totalnum, ntime, R, Rls1, Lls1, Ll, debugging, ties_method);
    //
    //
    if (verbose){
        end_point = system_clock::now();
        ending = time_point_cast<microseconds>(end_point).time_since_epoch().count();
        Rcout <<"df100 "<<(ending-start)<<" "<<0<<" "<<0<<" "<<0<<",Calc"<<endl;//prints the time
        gibtime = system_clock::to_time_t(system_clock::now());
        Rcout << ctime(&gibtime) << endl;
        Rcout << "df101 ";//prints the log-likelihoods
        for (int ij=0;ij<reqrdnum;ij++){
            Rcout << Ll[ij] << " ";
        }
        Rcout << " " << endl;
        Rcout << "df104 ";//prints parameter values
        for (int ij=0;ij<reqrdnum;ij++){
            Rcout << beta_0[ij] << " ";
        }
        Rcout << " " << endl;
    }
    //
    //
    List res_list = List::create(_["LogLik"]=wrap(Ll[0]),_["beta_0"]=wrap(beta_0) ,_["AIC"]=2*(totalnum-2*Ll[0]));
    // returns a list of results
    return res_list;
}

//' Primary poisson calculation function
//' \code{LogLik_Poisson_Single} Performs the calls to calculation functions, Structures the poisson calculation
//'
//' @param     PyrC    person-year matrix
//' @param     Term_n    Term numbers
//' @param     tform    subterm types
//' @param     a_n    starting values
//' @param     x_all    covariate matrix
//' @param     dfc    covariate column numbers
//' @param     fir    first term number
//' @param     modelform    model string
//' @param     verbose    verbosity boolean
//' @param     debugging    debugging boolean
//' @param     KeepConstant    vector identifying constant parameters
//' @param     term_tot    total number of terms
//' @param     nthreads number of threads to use
//'
//' @return List of results: Log-likelihood of optimum, first derivative of log-likelihood, second derivative matrix, parameter list, standard deviation estimate, AIC, deviance, model information
// [[Rcpp::export]]
List LogLik_Poisson_Single( MatrixXd PyrC, IntegerVector Term_n, StringVector tform, NumericVector a_n,NumericMatrix x_all,IntegerVector dfc,int fir,string modelform, bool verbose, bool debugging, IntegerVector KeepConstant, int term_tot, int nthreads){
    ;
    //
    List temp_list = List::create(_["Status"]="FAILED"); //used as a dummy return value for code checking
    if (verbose){
        Rcout << "START_POISSON_SINGLE" << endl;
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
    int reqrdnum = totalnum - sum(KeepConstant);
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
    // int nthreads = Eigen::nbThreads()-1; //stores how many threads are allocated
    //
    //
    // Lld_worst: stores the highest magnitude log-likelihood derivative
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
    //
    MatrixXd Te = MatrixXd::Zero(df0.rows(), 1); //preallocates matrix for column terms used for temporary storage
    MatrixXd R = MatrixXd::Zero(df0.rows(), 1); //preallocates matrix for Risks
    //
    MatrixXd Dose = MatrixXd::Constant(df0.rows(),term_tot,0.0);
    MatrixXd nonDose = MatrixXd::Constant(df0.rows(),term_tot,1.0); //Matrix of the total non-dose term values
    MatrixXd nonDose_LIN = MatrixXd::Constant(df0.rows(),term_tot,0.0);
    MatrixXd nonDose_PLIN = MatrixXd::Constant(df0.rows(),term_tot,1.0); //matrix of Loglinear subterm values
    MatrixXd nonDose_LOGLIN = MatrixXd::Constant(df0.rows(),term_tot,1.0); //matrix of Product linear subterm values
    MatrixXd TTerm = MatrixXd::Zero(Dose.rows(),Dose.cols()); //matrix of term values
    if (verbose){
        Rcout << "starting subterms " << term_tot << endl;
    }
    // Calculates the subterm and term values
    Make_subterms_Single( totalnum, Term_n, tform, dfc, fir, T0, Dose, nonDose, TTerm,  nonDose_LIN, nonDose_PLIN, nonDose_LOGLIN ,beta_0, df0,nthreads, debugging,KeepConstant);
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
    ColXd RdR = MatrixXd::Zero(df0.rows(), reqrdnum); //preallocates matrix for Risk to derivative ratios
    ColXd RddR = MatrixXd::Zero(df0.rows(), reqrdnum*(reqrdnum+1)/2);
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
    Make_Risks_Single(modelform, tform, Term_n, totalnum, fir, T0, Te, R, Dose, nonDose, TTerm, nonDose_LIN, nonDose_PLIN, nonDose_LOGLIN, nthreads, debugging,KeepConstant);
    //
    //
    //
    if (verbose){
        Rcout << "risk checked ";
        for (int ijk=0;ijk<1;ijk++){
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
    vector<double> Ll(reqrdnum,0.0);
    //
    //
    Poisson_LogLik_Single( nthreads, totalnum, PyrC, R, Ll, debugging);
    //
    MatrixXd dev_temp = MatrixXd::Zero(PyrC.rows(),2);
    dev_temp.col(0) = PyrC.col(0).array() * R.col(0).array();
    dev_temp.col(0) = PyrC.col(1).array() * dev_temp.col(0).array().pow(-1).array();
    dev_temp.col(0) = dev_temp.col(0).array().log().array();
    dev_temp.col(0) = PyrC.col(1).array() * dev_temp.col(0).array();
    dev_temp.col(1) = PyrC.col(1).array() - PyrC.col(0).array() * R.col(0).array();
    //
    dev_temp.col(0) = (dev_temp.col(0).array().isFinite()).select(dev_temp.col(0),0);
    //
    dev_temp.col(0) = dev_temp.col(0).array() - dev_temp.col(1).array();
    dev_temp.col(0) = (2 * dev_temp.col(0).array()).array();//.sqrt();
    dev_temp = (dev_temp.array().isFinite()).select(dev_temp,0);
    dev_temp.col(0) = (R.col(0).array()<0).select(0,dev_temp.col(0));
    double dev = dev_temp.col(0).sum(); //deviation calculation is split into steps
    //
    //
    if (verbose){
        end_point = system_clock::now();
        ending = time_point_cast<microseconds>(end_point).time_since_epoch().count();
        Rcout <<"df100 "<<(ending-start)<<" "<<0<<" "<<0<<" "<<0<<",Calc"<<endl;//prints the time
        gibtime = system_clock::to_time_t(system_clock::now());
        Rcout << ctime(&gibtime) << endl;
        Rcout << "df101 ";//prints the log-likelihoods
        for (int ij=0;ij<reqrdnum;ij++){
            Rcout << Ll[ij] << " ";
        }
        Rcout << " " << endl;
        Rcout << "df104 ";//prints parameter values
        for (int ij=0;ij<totalnum;ij++){
            Rcout << beta_0[ij] << " ";
        }
        Rcout << " " << endl;
    }
    if (R.minCoeff()<0){
        Rcout << "A non-positive risk was detected: " << R.minCoeff() << endl;
        Rcout << "final failing values ";
        for (int ijk=0;ijk<totalnum;ijk++){
            Rcout << beta_0[ijk] << " ";
        }
        Rcout << " " << endl;
        //
        Rcout << "final failing terms ";
        for (int ijk=0;ijk<totalnum;ijk++){
            Rcout << tform[ijk] << " ";
        }
        Rcout << " " << endl;
        temp_list = List::create(_["LogLik"]=R_NaN,_["beta_0"]=wrap(beta_0) ,_["Deviation"]=R_NaN,_["Status"]="FAILED");
        return temp_list;
    }
    //
    List res_list = List::create(_["LogLik"]=wrap(Ll[0]),_["beta_0"]=wrap(beta_0) ,_["AIC"]=2*(totalnum)+dev,_["Deviation"]=dev);
    // returns a list of results
    return res_list;
}



//' Primary Cox PH regression with competing risks, event=2
//' \code{LogLik_Cox_PH_CR} Performs the calls to calculation functions, Structures the Cox PH regression with competing events weighted by censoring rates, With verbose option prints out time stamps and intermediate sums of terms and derivatives
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
//' @param     cens_weight censoring weight list
//' @param     cens_thres threshold to add competing event to risk group
//' @param     double_step controls the step calculation, 0 for independent changes, 1 for solving b=Ax with complete matrices
//' @param     change_all    boolean if every parameter is being updated
//' @param     verbose    verbosity boolean
//' @param     debugging    debugging boolean
//' @param     KeepConstant    vector identifying constant parameters
//' @param     term_tot    total number of terms
//' @param     ties_method    ties method
//' @param     nthreads number of threads to use
//'
//' @return List of results: Log-likelihood of optimum, first derivative of log-likelihood, second derivative matrix, parameter list, standard deviation estimate, AIC, model information
// [[Rcpp::export]]
List LogLik_Cox_PH_CR( IntegerVector Term_n, StringVector tform, NumericVector a_n,NumericMatrix x_all,IntegerVector dfc,int fir, int der_iden,string modelform, double lr, int maxiter, int halfmax, double epsilon, double dbeta_cap, double abs_max,double dose_abs_max, double deriv_epsilon, NumericMatrix df_groups, NumericVector tu, const VectorXd cens_weight, const double cens_thres, int double_step ,bool change_all, bool verbose, bool debugging, IntegerVector KeepConstant, int term_tot, string ties_method, int nthreads){
    ;
    //
    List temp_list = List::create(_["Status"]="FAILED"); //used as a dummy return value for code checking
    if (verbose){
        Rcout << "START_NEW COMPETING RISKS" << endl;
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
    int reqrdnum = totalnum - sum(KeepConstant);
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
    // int nthreads = Eigen::nbThreads()-1; //stores how many threads are allocated
    //
    //
    // Lld_worst: stores the highest magnitude log-likelihood derivative
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
    MatrixXd Td0 = MatrixXd::Zero(df0.rows(), reqrdnum); //preallocates matrix for Term derivative columns
    MatrixXd Tdd0 = MatrixXd::Zero(df0.rows(), reqrdnum*(reqrdnum+1)/2); //preallocates matrix for Term second derivative columns
    //
    MatrixXd Te = MatrixXd::Zero(df0.rows(), 1); //preallocates matrix for column terms used for temporary storage
    MatrixXd R = MatrixXd::Zero(df0.rows(), 1); //preallocates matrix for Risks
    ColXd Rd = MatrixXd::Zero(df0.rows(), reqrdnum); //preallocates matrix for Risk derivatives
    ColXd Rdd = MatrixXd::Zero(df0.rows(), reqrdnum*(reqrdnum+1)/2); //preallocates matrix for Risk second derivatives
    //
    MatrixXd Dose = MatrixXd::Constant(df0.rows(),term_tot,0.0); //Matrix of the total dose term values
    MatrixXd nonDose = MatrixXd::Constant(df0.rows(),term_tot,1.0); //Matrix of the total non-dose term values
    MatrixXd nonDose_LIN = MatrixXd::Constant(df0.rows(),term_tot,0.0); //matrix of Linear subterm values
    MatrixXd nonDose_PLIN = MatrixXd::Constant(df0.rows(),term_tot,1.0); //matrix of Loglinear subterm values
    MatrixXd nonDose_LOGLIN = MatrixXd::Constant(df0.rows(),term_tot,1.0); //matrix of Product linear subterm values
    MatrixXd TTerm = MatrixXd::Zero(Dose.rows(),Dose.cols()); //matrix of term values
    double dint = dose_abs_max; //The amount of change used to calculate derivatives in threshold paramters
    double dslp = abs_max;
    //
    if (verbose){
        Rcout << "starting subterms " << term_tot << endl;
    }
    //
    // Calculates the subterm and term values
    Make_subterms( totalnum, Term_n, tform, dfc, fir, T0, Td0, Tdd0, Dose, nonDose, TTerm,  nonDose_LIN, nonDose_PLIN, nonDose_LOGLIN ,beta_0, df0,dint,dslp,nthreads, debugging,KeepConstant);
    // ---------------------------------------------------------
    // Prints off a series of calculations to check at what point values are changing
    // ---------------------------------------------------------
    //
    //
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
        for (int ijk=0;ijk<reqrdnum;ijk++){
            Rcout << Td0.col(ijk).sum() << " ";
        }
        Rcout << " " << endl;
        Rcout << "second derivs checked ";
        for (int ijk=0;ijk<reqrdnum;ijk++){
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
    ColXd RdR = MatrixXd::Zero(df0.rows(), reqrdnum); //preallocates matrix for Risk to derivative ratios
    ColXd RddR = MatrixXd::Zero(df0.rows(), reqrdnum*(reqrdnum+1)/2); //preallocates matrix for Risk to second derivative ratios
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
    Make_Risks(modelform, tform, Term_n, totalnum, fir, T0, Td0, Tdd0, Te, R, Rd, Rdd, Dose, nonDose, TTerm, nonDose_LIN, nonDose_PLIN, nonDose_LOGLIN, RdR, RddR, nthreads, debugging,KeepConstant);
    //
    // Removes infinite values
    RdR = (RdR.array().isFinite()).select(RdR,0);
    RddR = (RddR.array().isFinite()).select(RddR,0);
    //
    if (R.minCoeff()<0){
        Rcout << "A non-positive risk was detected: " << R.minCoeff() << endl;
        Rcout << "final failing values ";
        for (int ijk=0;ijk<totalnum;ijk++){
            Rcout << beta_0[ijk] << " ";
        }
        Rcout << " " << endl;
        //
        Rcout << "final failing terms ";
        for (int ijk=0;ijk<totalnum;ijk++){
            Rcout << tform[ijk] << " ";
        }
        Rcout << " " << endl;
        temp_list = List::create(_["beta_0"]=wrap(beta_0) ,_["Deviation"]=R_NaN,_["Status"]="FAILED",_["LogLik"]=R_NaN);
        return temp_list;
    }
    //
    if (verbose){
        Rcout << "risk checked ";
        for (int ijk=0;ijk<1;ijk++){
            Rcout << R.col(0).sum() << " ";
        }
        Rcout << " " << endl;
        Rcout << "risk1 checked ";
        for (int ijk=0;ijk<reqrdnum;ijk++){
            Rcout << Rd.col(ijk).sum() << " ";
        }
        Rcout << " " << endl;
        Rcout << "risk2 checked ";
        for (int ijk=0;ijk<reqrdnum;ijk++){
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
    Make_Groups_CR( ntime, df_m, RiskFail, RiskGroup, tu,cens_weight,cens_thres, nthreads, debugging);
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
    //
    MatrixXd Rls1 =MatrixXd::Zero(ntime, 1); //precomputes a series of sums used frequently in the log-liklihood calculations
    MatrixXd Rls2 =MatrixXd::Zero(ntime, reqrdnum); //Many are repeated due to the same risk groups and derivatives being used at mulitple points
    MatrixXd Rls3 =MatrixXd::Zero(ntime, reqrdnum*(reqrdnum+1)/2); //Sum and its derivatives are precomputed
    MatrixXd Lls1 =MatrixXd::Zero(ntime, 1); //The log-likelihood calculation has a Right and Left sum used
    MatrixXd Lls2 =MatrixXd::Zero(ntime, reqrdnum);
    MatrixXd Lls3 =MatrixXd::Zero(ntime, reqrdnum*(reqrdnum+1)/2);
    //The log-likelihood is calculated in parallel over the risk groups
    vector<double> Ll(reqrdnum,0.0); //Log-likelihood values
    vector<double> Lld(reqrdnum,0.0); //Log-likelihood derivative values
    vector<double> Lldd(pow(reqrdnum,2),0.0);//The second derivative matrix has room for every combination, but only the lower triangle is calculated initially
    //
    if (verbose){
        Rcout << "Made Risk Side Lists" << endl;
    }
    // Calculates the side sum terms used
    Calculate_Sides_CR( RiskFail, RiskGroup, totalnum, ntime, R, Rd, Rdd, Rls1, Rls2, Rls3, Lls1, Lls2, Lls3,cens_weight,nthreads, debugging,KeepConstant);
    //
    if (verbose){
        end_point = system_clock::now();
        ending = time_point_cast<microseconds>(end_point).time_since_epoch().count();
        Rcout <<"df100 "<<(ending-start)<<" "<<0<<" "<<0<<" "<<-1<<",Prep_Sides"<<endl;
        gibtime = system_clock::to_time_t(system_clock::now());
        Rcout << ctime(&gibtime) << endl;
    }
    //
    if (verbose){
        Rcout << "riskr checked ";
        for (int ijk=0;ijk<reqrdnum;ijk++){
            Rcout << Rls1.col(0).sum() << " ";
        }
        Rcout << " " << endl;
        Rcout << "risk1r checked ";
        for (int ijk=0;ijk<reqrdnum;ijk++){
            Rcout << Rls2.col(ijk).sum() << " ";
        }
        Rcout << " " << endl;
        Rcout << "risk2r checked ";
        for (int ijk=0;ijk<reqrdnum;ijk++){
            Rcout << Rls3.col(ijk*(ijk+1)/2+ijk).sum() << " ";
        }
        Rcout << " " << endl;
        //
        Rcout << "riskl checked ";
        for (int ijk=0;ijk<reqrdnum;ijk++){
            Rcout << Lls1.col(0).sum() << " ";
        }
        Rcout << " " << endl;
        Rcout << "risk1l checked ";
        for (int ijk=0;ijk<reqrdnum;ijk++){
            Rcout << Lls2.col(ijk).sum() << " ";
        }
        Rcout << " " << endl;
        Rcout << "risk2l checked ";
        for (int ijk=0;ijk<reqrdnum;ijk++){
            Rcout << Lls3.col(ijk*(ijk+1)/2+ijk).sum() << " ";
        }
        Rcout << " " << endl;
    }
    //
    //
    // Calculates log-likelihood
    Calc_LogLik( nthreads, RiskFail, RiskGroup, totalnum, ntime, R, Rd, Rdd,RdR,RddR, Rls1, Rls2, Rls3, Lls1, Lls2, Lls3, Ll, Lld, Lldd, debugging, ties_method,KeepConstant);
    //
    vector <double> Ll_comp(2,Ll[0]); //vector to compare values//
    if (verbose){
        end_point = system_clock::now();
        ending = time_point_cast<microseconds>(end_point).time_since_epoch().count();
        Rcout <<"df100 "<<(ending-start)<<" "<<0<<" "<<0<<" "<<0<<",Calc"<<endl;//prints the time
        gibtime = system_clock::to_time_t(system_clock::now());
        Rcout << ctime(&gibtime) << endl;
        Rcout << "df101 ";//prints the log-likelihoods
        for (int ij=0;ij<reqrdnum;ij++){
            Rcout << Ll[ij] << " ";
        }
        Rcout << " " << endl;
        Rcout << "df102 ";//prints the first derivatives
        for (int ij=0;ij<reqrdnum;ij++){
            Rcout << Lld[ij] << " ";
        }
        Rcout << " " << endl;
        Rcout << "df103 ";//prints the second derivatives
        for (int ij=0;ij<reqrdnum;ij++){
            Rcout << Lldd[ij*reqrdnum+ij] << " ";
        }
        Lld_worst=0;
        for (int ij=0;ij<reqrdnum;ij++){//locates highest magnitude derivative
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
        for (int ij=0;ij<reqrdnum;ij++){//prints the newton step value for zero derivative
            Rcout << Lld[ij]/Lldd[ij*reqrdnum+ij] << " ";
        }
        Rcout << " " << endl;
        Rcout << "df106 ";
        for (int ij=0;ij<reqrdnum;ij++){//prints the newton step value for zero log-likelihood
            Rcout << Ll[ij]/Lld[ij] << " ";
        }
        Rcout << " " << endl;
        Rcout << "df107 " << double_step << " " << abs_max << " " << dose_abs_max << " " << Ll_comp[0] << " " << Ll_comp[1] << endl;//prints several convergence terms
    }
    //
    vector<double> dbeta(totalnum,0.0);
    //
    // --------------------------
    // always starts from initial guess
    // --------------------------
    vector<double> beta_c(totalnum,0.0);
    vector<double> beta_a(totalnum,0.0);
    vector<double> beta_best(totalnum,0.0);
    vector<double> beta_p(totalnum,0.0);
    //int //risk_check_iter=0;
    VectorXd::Map(&beta_p[0], beta_0.size()) = beta_0;// stores previous parameters
    VectorXd::Map(&beta_c[0], beta_0.size()) = beta_0;// stores current parameters
    VectorXd::Map(&beta_a[0], beta_0.size()) = beta_0;// stores a refrence value for parameters
    VectorXd::Map(&beta_best[0], beta_0.size()) = beta_0;// stores the best parameters
    double halves = 0; //number of half-steps taken
    int ind0 = fir; //used for validations
    //int //i = ind0;
    int iteration=0; //iteration number
    //
    bool convgd = FALSE;
    int iter_stop =0; //tracks if the iterations should be stopped for convergence
    int iter_check=0; //signal to check for convergence
    //
    double Ll_abs_best = 10;
    vector<double> beta_abs_best(totalnum,0.0);
    while ((iteration < maxiter)&&(iter_stop==0)){
        iteration++;
        beta_p = beta_c;//
        beta_a = beta_c;//
        beta_best = beta_c;//
        //
        // Calculates the initial change in parameter
        Calc_Change( double_step, nthreads, totalnum, fir, der_iden, dbeta_cap, dose_abs_max, lr, abs_max, Ll, Lld, Lldd, dbeta, change_all, tform, dint,dslp, KeepConstant, debugging);
        Intercept_Bound(nthreads, totalnum, beta_0, dbeta, dfc, df0, KeepConstant, debugging, tform);
        if (verbose){
            Rcout << "Starting Halves"<<endl;//prints the final changes for validation
        }
        //
        //
        halves=0;
        //
        if ((Ll_abs_best>0)||(Ll_abs_best < Ll[ind0])){
            Ll_abs_best = Ll[ind0];
            beta_abs_best = beta_c;
        }
        //
        if (verbose){
            Rcout << "df501 " << Ll_abs_best << endl;
            Rcout << "df504 ";//prints parameter values
            for (int ij=0;ij<totalnum;ij++){
                Rcout << beta_abs_best[ij] << " ";
            }
            Rcout << " " << endl;
        }
        //
        while ((Ll[ind0] <= Ll_abs_best)&&(halves<halfmax)){ //repeats until half-steps maxed or an improvement
            //Refreshes the matrices used
            Dose = MatrixXd::Zero(df0.rows(),term_tot);
            nonDose = MatrixXd::Constant(df0.rows(),term_tot,0.0);
            nonDose_LIN = MatrixXd::Constant(df0.rows(),term_tot,0.0);
            nonDose_PLIN = MatrixXd::Constant(df0.rows(),term_tot,1.0);
            nonDose_LOGLIN = MatrixXd::Constant(df0.rows(),term_tot,1.0);
            T0 = MatrixXd::Zero(df0.rows(), totalnum); //preallocates matrix for Term column
            Td0 = MatrixXd::Zero(df0.rows(), reqrdnum); //preallocates matrix for Term derivative columns
            Tdd0 = MatrixXd::Zero(df0.rows(), reqrdnum*(reqrdnum+1)/2); //preallocates matrix for Term second derivative columns
            for (int ij=0;ij<totalnum;ij++){
                beta_0[ij] = beta_a[ij] + dbeta[ij];
                beta_c[ij] = beta_0[ij];
            }
            // ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------//
            // The same subterm, risk, sides, and log-likelihood calculations are performed every half-step and iteration
            // ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------//
            Make_subterms( totalnum, Term_n, tform, dfc, fir, T0, Td0, Tdd0, Dose,  nonDose, TTerm,  nonDose_LIN, nonDose_PLIN, nonDose_LOGLIN ,beta_0, df0,dint,dslp,nthreads, debugging,KeepConstant);
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
                for (int ijk=0;ijk<reqrdnum;ijk++){
                    Rcout << Td0.col(ijk).sum() << " ";
                }
                Rcout << " " << endl;
                Rcout << "second derivs checked ";
                for (int ijk=0;ijk<reqrdnum;ijk++){
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
            RdR = MatrixXd::Zero(df0.rows(), reqrdnum);
            RddR = MatrixXd::Zero(df0.rows(), reqrdnum*(reqrdnum+1)/2);
            //
            //
            Make_Risks(modelform, tform, Term_n, totalnum, fir, T0, Td0, Tdd0, Te, R, Rd, Rdd, Dose, nonDose, TTerm, nonDose_LIN, nonDose_PLIN, nonDose_LOGLIN, RdR, RddR, nthreads, debugging,KeepConstant);
            //
            if (R.minCoeff()<0){
                #pragma omp parallel for num_threads(nthreads)
                for (int ijk=0;ijk<totalnum;ijk++){
                    int tij = Term_n[ijk];
                    if (TTerm.col(tij).minCoeff()<=0){
                        dbeta[ijk] = dbeta[ijk] / 2.0;
                    } else if (isinf(TTerm.col(tij).maxCoeff())){
                        dbeta[ijk] = dbeta[ijk] / 2.0;
                    }
                }
                if (verbose){
                    Rcout << "A non-positive risk was detected: " << R.minCoeff() << endl;
                }
                halves+=0.2;
            } else {
                halves++;
                RdR = (RdR.array().isFinite()).select(RdR,0);
                RddR = (RddR.array().isFinite()).select(RddR,0);
                //
                if (verbose){
                    Rcout << "risk checked ";
                    for (int ijk=0;ijk<1;ijk++){
                        Rcout << R.col(0).sum() << " ";
                    }
                    Rcout << " " << endl;
                    Rcout << "risk1 checked ";
                    for (int ijk=0;ijk<reqrdnum;ijk++){
                        Rcout << Rd.col(ijk).sum() << " ";
                    }
                    Rcout << " " << endl;
                    Rcout << "risk2 checked ";
                    for (int ijk=0;ijk<reqrdnum;ijk++){
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
                Rls2 =MatrixXd::Zero(ntime, reqrdnum); //Many are repeated due to the same risk groups and derivatives being used at mulitple points
                Rls3 =MatrixXd::Zero(ntime, reqrdnum*(reqrdnum+1)/2);
                Lls1 =MatrixXd::Zero(ntime, 1);
                Lls2 =MatrixXd::Zero(ntime, reqrdnum);
                Lls3 =MatrixXd::Zero(ntime, reqrdnum*(reqrdnum+1)/2);
                Calculate_Sides_CR( RiskFail, RiskGroup, totalnum, ntime, R, Rd, Rdd, Rls1, Rls2, Rls3, Lls1, Lls2, Lls3,cens_weight,nthreads, debugging,KeepConstant);
                //
                
                if (verbose){
                    end_point = system_clock::now();
                    ending = time_point_cast<microseconds>(end_point).time_since_epoch().count();
                    Rcout <<"df100 "<<(ending-start)<<" "<<halves<<" "<<iteration<<" "<<ind0<<",Update_Sides"<<endl;
                    gibtime = system_clock::to_time_t(system_clock::now());
                    Rcout << ctime(&gibtime) << endl;
                    Rcout << "riskr checked ";
                    for (int ijk=0;ijk<reqrdnum;ijk++){
                        Rcout << Rls1.col(0).sum() << " ";
                    }
                    Rcout << " " << endl;
                    Rcout << "risk1r checked ";
                    for (int ijk=0;ijk<reqrdnum;ijk++){
                        Rcout << Rls2.col(ijk).sum() << " ";
                    }
                    Rcout << " " << endl;
                    Rcout << "risk2r checked ";
                    for (int ijk=0;ijk<reqrdnum;ijk++){
                        Rcout << Rls3.col(ijk*(ijk+1)/2+ijk).sum() << " ";
                    }
                    Rcout << " " << endl;
                    //
                    //
                    Rcout << "riskl checked ";
                    for (int ijk=0;ijk<reqrdnum;ijk++){
                        Rcout << Lls1.col(0).sum() << " ";
                    }
                    Rcout << " " << endl;
                    Rcout << "risk1l checked ";
                    for (int ijk=0;ijk<reqrdnum;ijk++){
                        Rcout << Lls2.col(ijk).sum() << " ";
                    }
                    Rcout << " " << endl;
                    Rcout << "risk2l checked ";
                    for (int ijk=0;ijk<reqrdnum;ijk++){
                        Rcout << Lls3.col(ijk*(ijk+1)/2+ijk).sum() << " ";
                    }
                    Rcout << " " << endl;
                }
                //
                Calc_LogLik( nthreads, RiskFail, RiskGroup, totalnum, ntime, R, Rd, Rdd,RdR,RddR, Rls1, Rls2, Rls3, Lls1, Lls2, Lls3, Ll, Lld, Lldd, debugging, ties_method,KeepConstant);
                if (verbose){
                    end_point = system_clock::now();
                    ending = time_point_cast<microseconds>(end_point).time_since_epoch().count();
                    Rcout <<"df100 "<<(ending-start)<<" "<<0<<" "<<0<<" "<<0<<",Half Calc"<<endl;//prints the time
                    gibtime = system_clock::to_time_t(system_clock::now());
                    Rcout << ctime(&gibtime) << endl;
                    Rcout << "df101 ";//prints the log-likelihoods
                    for (int ij=0;ij<reqrdnum;ij++){
                        Rcout << Ll[ij] << " ";
                    }
                    Rcout << " " << endl;
                    Rcout << "df102 ";//prints the first derivatives
                    for (int ij=0;ij<reqrdnum;ij++){
                        Rcout << Lld[ij] << " ";
                    }
                    Rcout << " " << endl;
                    Rcout << "df103 ";//prints the second derivatives
                    for (int ij=0;ij<reqrdnum;ij++){
                        Rcout << Lldd[ij*reqrdnum+ij] << " ";
                    }
                    Lld_worst=0;
                    for (int ij=0;ij<reqrdnum;ij++){//locates highest magnitude derivative
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
                    for (int ij=0;ij<reqrdnum;ij++){//prints the newton step value for zero derivative
                        Rcout << Lld[ij]/Lldd[ij*reqrdnum+ij] << " ";
                    }
                    Rcout << " " << endl;
                    Rcout << "df106 ";
                    for (int ij=0;ij<reqrdnum;ij++){//prints the newton step value for zero log-likelihood
                        Rcout << Ll[ij]/Lld[ij] << " ";
                    }
                    Rcout << " " << endl;
                    Rcout << "df107 " << double_step << " " << abs_max << " " << dose_abs_max << " " << Ll_comp[0] << " " << Ll_comp[1] << endl;//prints several convergence terms
                }
                if (change_all){ //If every covariate is to be changed
                    if (Ll[ind0] <= Ll_abs_best){//takes a half-step if needed
                        #pragma omp parallel for num_threads(nthreads)
                        for (int ijk=0;ijk<totalnum;ijk++){
                            dbeta[ijk] = dbeta[ijk] * 0.5; //
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
                }
                #pragma omp parallel for num_threads(nthreads)
                for (int ijk=0;ijk<totalnum;ijk++){//totalnum*(totalnum+1)/2
                    beta_0[ijk] = beta_c[ijk];
                }
            }
        }
        if (beta_best!=beta_c){//if the risk matrices aren't the optimal values, then they must be recalculated
            if (verbose){
                Rcout << "Changing back to best"<<endl;
            }
            // If it goes through every half step without improvement, then the maximum change needs to be decreased
            abs_max = abs_max*pow(0.5,halfmax); // reduces the step sizes
            dose_abs_max = dose_abs_max*pow(0.5,halfmax);
            iter_check = 1;
            //
            beta_p = beta_best;//
            beta_a = beta_best;//
            beta_c = beta_best;//
            Dose = MatrixXd::Zero(df0.rows(),term_tot);
            nonDose = MatrixXd::Constant(df0.rows(),term_tot,0.0);
            nonDose_LIN = MatrixXd::Constant(df0.rows(),term_tot,0.0);
            nonDose_PLIN = MatrixXd::Constant(df0.rows(),term_tot,1.0);
            nonDose_LOGLIN = MatrixXd::Constant(df0.rows(),term_tot,1.0);
            for (int ij=0;ij<totalnum;ij++){
                beta_0[ij] = beta_best[ij];
            }
            Make_subterms( totalnum, Term_n, tform, dfc, fir, T0, Td0, Tdd0, Dose, nonDose, TTerm,  nonDose_LIN, nonDose_PLIN, nonDose_LOGLIN,beta_0, df0,dint,dslp,nthreads, debugging, KeepConstant);
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
                for (int ijk=0;ijk<reqrdnum;ijk++){
                    Rcout << Td0.col(ijk).sum() << " ";
                }
                Rcout << " " << endl;
                Rcout << "second derivs checked ";
                for (int ijk=0;ijk<reqrdnum;ijk++){
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
            RdR = MatrixXd::Zero(df0.rows(), reqrdnum);
            RddR = MatrixXd::Zero(df0.rows(), reqrdnum*(reqrdnum+1)/2);
            //
            //
            Make_Risks(modelform, tform, Term_n, totalnum, fir, T0, Td0, Tdd0, Te, R, Rd, Rdd, Dose, nonDose, TTerm, nonDose_LIN, nonDose_PLIN, nonDose_LOGLIN, RdR, RddR, nthreads, debugging,KeepConstant);
            R = (R.array().isFinite()).select(R,0);
            Rd = (Rd.array().isFinite()).select(Rd,0);
            Rdd = (Rdd.array().isFinite()).select(Rdd,0);
            RdR = (RdR.array().isFinite()).select(RdR,0);
            RddR = (RddR.array().isFinite()).select(RddR,0);
            //
            temp_list = List::create(_["betas"]=wrap(beta_0),_["Status"]="FAILED"); //used as a dummy return value for code checking
            if (verbose){
                Rcout << "risk checked ";
                for (int ijk=0;ijk<1;ijk++){
                    Rcout << R.col(0).sum() << " ";
                }
                Rcout << " " << endl;
                Rcout << "risk1 checked ";
                for (int ijk=0;ijk<reqrdnum;ijk++){
                    Rcout << Rd.col(ijk).sum() << " ";
                }
                Rcout << " " << endl;
                Rcout << "risk2 checked ";
                for (int ijk=0;ijk<reqrdnum;ijk++){
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
            Rls2 =MatrixXd::Zero(ntime, reqrdnum); //Many are repeated due to the same risk groups and derivatives being used at mulitple points
            Rls3 =MatrixXd::Zero(ntime, reqrdnum*(reqrdnum+1)/2);
            Lls1 =MatrixXd::Zero(ntime, 1);
            Lls2 =MatrixXd::Zero(ntime, reqrdnum);
            Lls3 =MatrixXd::Zero(ntime, reqrdnum*(reqrdnum+1)/2);
            Calculate_Sides_CR( RiskFail, RiskGroup, totalnum, ntime, R, Rd, Rdd, Rls1, Rls2, Rls3, Lls1, Lls2, Lls3,cens_weight,nthreads, debugging,KeepConstant);
            //
            if (verbose){
                end_point = system_clock::now();
                ending = time_point_cast<microseconds>(end_point).time_since_epoch().count();
                Rcout <<"df100 "<<(ending-start)<<" "<<halves<<" "<<iteration<<" "<<ind0<<",Update_Sides"<<endl;
                gibtime = system_clock::to_time_t(system_clock::now());
                Rcout << ctime(&gibtime) << endl;
                Rcout << "riskr checked ";
                for (int ijk=0;ijk<1;ijk++){
                    Rcout << Rls1.col(0).sum() << " ";
                }
                Rcout << " " << endl;
                Rcout << "risk1r checked ";
                for (int ijk=0;ijk<reqrdnum;ijk++){
                    Rcout << Rls2.col(ijk).sum() << " ";
                }
                Rcout << " " << endl;
                Rcout << "risk2r checked ";
                for (int ijk=0;ijk<reqrdnum;ijk++){
                    Rcout << Rls3.col(ijk*(ijk+1)/2+ijk).sum() << " ";
                }
                Rcout << " " << endl;
                //
                //
                Rcout << "riskl checked ";
                for (int ijk=0;ijk<1;ijk++){
                    Rcout << Lls1.col(0).sum() << " ";
                }
                Rcout << " " << endl;
                Rcout << "risk1l checked ";
                for (int ijk=0;ijk<reqrdnum;ijk++){
                    Rcout << Lls2.col(ijk).sum() << " ";
                }
                Rcout << " " << endl;
                Rcout << "risk2l checked ";
                for (int ijk=0;ijk<reqrdnum;ijk++){
                    Rcout << Lls3.col(ijk*(ijk+1)/2+ijk).sum() << " ";
                }
                Rcout << " " << endl;
            }
            //
            Calc_LogLik( nthreads, RiskFail, RiskGroup, totalnum, ntime, R, Rd, Rdd,RdR,RddR, Rls1, Rls2, Rls3, Lls1, Lls2, Lls3, Ll, Lld, Lldd, debugging, ties_method,KeepConstant);
        }
        Lld_worst=0;
        for (int ij=0;ij<reqrdnum;ij++){
            if (abs(Lld[ij]) > Lld_worst){
                Lld_worst = abs(Lld[ij]);
            }
        }
        if (iteration > reqrdnum){//Doesn't check the first several iterations for convergence
            if ((iteration % (reqrdnum))||(iter_check==1)){//Checks every set number of iterations
                iter_check=0;
                if (Lld_worst < deriv_epsilon){//ends if the derivatives are low enough
                    iter_stop = 1;
                    convgd = TRUE;
                }
                Ll_comp[1]=Ll[0];
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
            for (int ij=0;ij<reqrdnum;ij++){
                Rcout << Ll[ij] << " ";
            }
            Rcout << " " << endl;
            Rcout << "df102 ";
            for (int ij=0;ij<reqrdnum;ij++){
                Rcout << Lld[ij] << " ";
            }
            Rcout << " " << endl;
            Rcout << "df103 ";
            for (int ij=0;ij<reqrdnum;ij++){
                Rcout << Lldd[ij*reqrdnum+ij] << " ";
            }
            Rcout << " " << endl;
            Rcout << "df104 ";
            for (int ij=0;ij<totalnum;ij++){
                Rcout << beta_c[ij] << " ";
            }
            Rcout << " " << endl;
            Rcout << "df105 ";
            for (int ij=0;ij<reqrdnum;ij++){//prints the newton step value for zero derivative
                Rcout << Lld[ij]/Lldd[ij*reqrdnum+ij] << " ";
            }
            Rcout << " " << endl;
            Rcout << "df106 ";
            for (int ij=0;ij<reqrdnum;ij++){//prints the newton step value for zero log-likelihood
                Rcout << Ll[ij]/Lld[ij] << " ";
            }
            Rcout << " " << endl;
            Rcout << "df107 " << double_step << " " << abs_max << " " << dose_abs_max << " " << Ll_comp[0] << " " << Ll_comp[1] << endl;//prints several convergence terms
            Rcout << "Finshed iteration" << endl;
        }
    }
    // -----------------------------------------------
    // Performing Full Calculation to get full second derivative matrix
    // -----------------------------------------------
    fill(Ll.begin(), Ll.end(), 0.0);
    fill(Lld.begin(), Lld.end(), 0.0);
    fill(Lldd.begin(), Lldd.end(), 0.0);
    Calculate_Sides_CR( RiskFail, RiskGroup, totalnum, ntime, R, Rd, Rdd, Rls1, Rls2, Rls3, Lls1, Lls2, Lls3,cens_weight,nthreads, debugging,KeepConstant);
    //
    if (verbose){
        end_point = system_clock::now();
        ending = time_point_cast<microseconds>(end_point).time_since_epoch().count();
        Rcout <<"df100 "<<(ending-start)<<" "<<halves<<" "<<iteration<<" "<<ind0<<",Update_Sides"<<endl;
        gibtime = system_clock::to_time_t(system_clock::now());
        Rcout << ctime(&gibtime) << endl;
        Rcout << "Wrapping up" << endl;
    }
    //
    Calc_LogLik( nthreads, RiskFail, RiskGroup, totalnum, ntime, R, Rd, Rdd,RdR,RddR, Rls1, Rls2, Rls3, Lls1, Lls2, Lls3, Ll, Lld, Lldd, debugging, ties_method,KeepConstant);
    //
    if ((Ll_abs_best>0)||(Ll_abs_best < Ll[ind0])){
        Ll_abs_best = Ll[ind0];
        beta_abs_best = beta_c;
    }
    //
    if (verbose){
        Rcout << "df501 " << Ll_abs_best << endl;
        Rcout << "df504 ";//prints parameter values
        for (int ij=0;ij<totalnum;ij++){
            Rcout << beta_abs_best[ij] << " ";
        }
        Rcout << " " << endl;
    }
    //
    List para_list = List::create(_["Term_n"]=Term_n,_["tforms"]=tform); //stores the term information
    List control_list = List::create(_["Iteration"]=iteration, _["Maximum Step"]=abs_max, _["Derivative Limiting"]=Lld_worst); //stores the total number of iterations used
    //
    int kept_covs = totalnum - sum(KeepConstant);
    NumericVector Lldd_vec(kept_covs * kept_covs);
    #pragma omp parallel for schedule(dynamic) num_threads(nthreads)
    for (int ijk=0;ijk<totalnum*(totalnum+1)/2;ijk++){
        int ij = 0;
        int jk = ijk;
        int pij_ind=-100;
        int pjk_ind=-100;
        while (jk>ij){
            ij++;
            jk-=ij;
        }
        if (KeepConstant[ij]==0){
            pij_ind = ij - sum(head(KeepConstant,ij));
            if (KeepConstant[jk]==0){
                pjk_ind = jk - sum(head(KeepConstant,jk));
                Lldd_vec[pij_ind * kept_covs + pjk_ind]=Lldd[ij*totalnum+jk];
            }
        }
    }
    for (int ijk=0;ijk<kept_covs*(kept_covs+1)/2;ijk++){
        int ij = 0;
        int jk = ijk;
        while (jk>ij){
            ij++;
            jk-=ij;
        }
        Lldd_vec[jk * kept_covs + ij]=Lldd_vec[ij * kept_covs + jk];
    }
    Lldd_vec.attr("dim") = Dimension(kept_covs, kept_covs);
    const Map<MatrixXd> Lldd_mat(as<Map<MatrixXd> >(Lldd_vec));
    //
    MatrixXd Lldd_inv = -1 * Lldd_mat.inverse().matrix(); //uses inverse information matrix to calculate the standard deviation
    VectorXd stdev = VectorXd::Zero(totalnum);
    for (int ij=0;ij<totalnum;ij++){
        if (KeepConstant[ij]==0){
            int pij_ind = ij - sum(head(KeepConstant,ij));
            stdev(ij) = sqrt(Lldd_inv(pij_ind,pij_ind));
        }
    }
    //
    //
    List res_list = List::create(_["LogLik"]=wrap(Ll[0]),_["First_Der"]=wrap(Lld),_["Second_Der"]=Lldd_vec,_["beta_0"]=wrap(beta_0) ,_["Standard_Deviation"]=wrap(stdev) ,_["AIC"]=2*(totalnum-accumulate(KeepConstant.begin(),KeepConstant.end(), 0.0))-2*Ll[0],_["Parameter_Lists"]=para_list,_["Control_List"]=control_list,_["Convgerged"]=convgd);
    // returns a list of results
    return res_list;
}




//' Primary Cox PH regression with multiple starting points
//' \code{LogLik_Cox_PH_Guess} Performs the calls to calculation functions, Structures the Cox PH regression, With verbose option prints out time stamps and intermediate sums of terms and derivatives
//'
//' @param     Term_n    Term numbers
//' @param     tform    subterm types
//' @param     a_ns    matrix starting values
//' @param     x_all    covariate matrix
//' @param     dfc    covariate column numbers
//' @param     fir    first term number
//' @param     der_iden    subterm number for derivative tests
//' @param     modelform    model string
//' @param     lr    learning rate for newton step toward 0 derivative
//' @param     maxiters    list of maximum number of iterations
//' @param     guesses    number of initial guesses
//' @param     halfmax    maximum number of half steps
//' @param     epsilon    minimum acceptable maximum parameter change
//' @param     dbeta_cap    learning rate for newton step toward 0 log-likelihood
//' @param     abs_max    Maximum allowed parameter change
//' @param     dose_abs_max    Maximum allowed threshold parameter change
//' @param     deriv_epsilon    threshold for near-zero derivative
//' @param     df_groups    matrix with time and event information
//' @param     tu    event times
//' @param     double_step controls the step calculation, 0 for independent changes, 1 for solving b=Ax with complete matrices
//' @param     change_all    boolean if every parameter is being updated
//' @param     verbose    verbosity boolean
//' @param     debugging    debugging boolean
//' @param     KeepConstant    vector identifying constant parameters
//' @param     term_tot    total number of terms
//' @param     ties_method    ties method
//' @param     nthreads number of threads to use
//'
//' @return List of final results: Log-likelihood of optimum, first derivative of log-likelihood, second derivative matrix, parameter list, standard deviation estimate, AIC, model information
// [[Rcpp::export]]
List LogLik_Cox_PH_Guess( IntegerVector Term_n, StringVector tform, NumericMatrix a_ns,NumericMatrix x_all,IntegerVector dfc,int fir, int der_iden,string modelform, double lr, NumericVector maxiters, int guesses, int halfmax, double epsilon, double dbeta_cap, double abs_max,double dose_abs_max, double deriv_epsilon, NumericMatrix df_groups, NumericVector tu, int double_step ,bool change_all, bool verbose, bool debugging, IntegerVector KeepConstant, int term_tot, string ties_method, int nthreads){
    ;
    //
    List temp_list = List::create(_["Status"]="FAILED"); //used as a dummy return value for code checking
    if (verbose){
        Rcout << "START_COX_GUESS" << endl;
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
    int reqrdnum = totalnum - sum(KeepConstant);
    //
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
    // int nthreads = Eigen::nbThreads()-1; //stores how many threads are allocated
    //
    // Lld_worst: stores the highest magnitude log-likelihood derivative
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
    NumericVector a_n = a_ns.row(0);
    Map<VectorXd> beta_0(as<Map<VectorXd> >(a_n));
    int maxiter = maxiters[0];
    MatrixXd T0 = MatrixXd::Zero(df0.rows(), totalnum); //preallocates matrix for Term column
    MatrixXd Td0 = MatrixXd::Zero(df0.rows(), reqrdnum); //preallocates matrix for Term derivative columns
    MatrixXd Tdd0 = MatrixXd::Zero(df0.rows(), reqrdnum*(reqrdnum+1)/2); //preallocates matrix for Term second derivative columns
    //
    MatrixXd Te = MatrixXd::Zero(df0.rows(), 1); //preallocates matrix for column terms used for temporary storage
    MatrixXd R = MatrixXd::Zero(df0.rows(), 1); //preallocates matrix for Risks
    ColXd Rd = MatrixXd::Zero(df0.rows(), reqrdnum); //preallocates matrix for Risk derivatives
    ColXd Rdd = MatrixXd::Zero(df0.rows(), reqrdnum*(reqrdnum+1)/2); //preallocates matrix for Risk second derivatives
    //
    MatrixXd Dose = MatrixXd::Constant(df0.rows(),term_tot,0.0); //Matrix of the total dose term values
    MatrixXd nonDose = MatrixXd::Constant(df0.rows(),term_tot,1.0); //Matrix of the total non-dose term values
    MatrixXd nonDose_LIN = MatrixXd::Constant(df0.rows(),term_tot,0.0); //matrix of Linear subterm values
    MatrixXd nonDose_PLIN = MatrixXd::Constant(df0.rows(),term_tot,1.0); //matrix of Loglinear subterm values
    MatrixXd nonDose_LOGLIN = MatrixXd::Constant(df0.rows(),term_tot,1.0); //matrix of Product linear subterm values
    MatrixXd TTerm = MatrixXd::Zero(Dose.rows(),Dose.cols()); //matrix of term values
    double dint = dose_abs_max; //The amount of change used to calculate derivatives in threshold paramters
    double dslp = abs_max;
    ColXd RdR = MatrixXd::Zero(df0.rows(), reqrdnum); //preallocates matrix for Risk to derivative ratios
    ColXd RddR = MatrixXd::Zero(df0.rows(), reqrdnum*(reqrdnum+1)/2); //preallocates matrix for Risk to second derivative ratios
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
    //
    MatrixXd Rls1 =MatrixXd::Zero(ntime, 1); //precomputes a series of sums used frequently in the log-liklihood calculations
    MatrixXd Rls2 =MatrixXd::Zero(ntime, reqrdnum); //Many are repeated due to the same risk groups and derivatives being used at mulitple points
    MatrixXd Rls3 =MatrixXd::Zero(ntime, reqrdnum*(reqrdnum+1)/2); //Sum and its derivatives are precomputed
    MatrixXd Lls1 =MatrixXd::Zero(ntime, 1); //The log-likelihood calculation has a Right and Left sum used
    MatrixXd Lls2 =MatrixXd::Zero(ntime, reqrdnum);
    MatrixXd Lls3 =MatrixXd::Zero(ntime, reqrdnum*(reqrdnum+1)/2);
    //The log-likelihood is calculated in parallel over the risk groups
    vector<double> Ll(reqrdnum,0.0); //Log-likelihood values
    vector<double> Lld(reqrdnum,0.0); //Log-likelihood derivative values
    vector<double> Lldd(pow(reqrdnum,2),0.0);//The second derivative matrix has room for every combination, but only the lower triangle is calculated initially
    //
    vector <double> Ll_comp(2,Ll[0]); //vector to compare values
    double abs_max0 = abs_max;
    double dose_abs_max0 = dose_abs_max;
    //
    vector<double> dbeta(totalnum,0.0);
    //
    // --------------------------
    // always starts from initial guess
    // --------------------------
    vector<double> beta_c(totalnum,0.0);
    vector<double> beta_a(totalnum,0.0);
    vector<double> beta_best(totalnum,0.0);
    vector<double> beta_p(totalnum,0.0);
    VectorXd::Map(&beta_p[0], beta_0.size()) = beta_0;// stores previous parameters
    VectorXd::Map(&beta_c[0], beta_0.size()) = beta_0;// stores current parameters
    VectorXd::Map(&beta_a[0], beta_0.size()) = beta_0;// stores a refrence value for parameters
    VectorXd::Map(&beta_best[0], beta_0.size()) = beta_0;// stores the best parameters
    double halves = 0; //number of half-steps taken
    int ind0 = fir; //used for validations
    int iteration=0; //iteration number
    //
    bool convgd = FALSE;
    int iter_stop =0; //tracks if the iterations should be stopped for convergence
    int iter_check=0; //signal to check for convergence
    //
    NumericMatrix beta_fin(a_ns.rows(), a_ns.cols());
    NumericVector LL_fin(a_ns.rows());
    //
        //
    double Ll_abs_best = 10;
    vector<double> beta_abs_best(totalnum,0.0);
    int guess_abs_best =-1;
    for (int guess=0; guess<=guesses;guess++){
        //
        Dose = MatrixXd::Zero(df0.rows(),term_tot);
        nonDose = MatrixXd::Constant(df0.rows(),term_tot,0.0);
        nonDose_LIN = MatrixXd::Constant(df0.rows(),term_tot,0.0);
        nonDose_PLIN = MatrixXd::Constant(df0.rows(),term_tot,1.0);
        nonDose_LOGLIN = MatrixXd::Constant(df0.rows(),term_tot,1.0);
        T0 = MatrixXd::Zero(df0.rows(), totalnum); //preallocates matrix for Term column
        Td0 = MatrixXd::Zero(df0.rows(), reqrdnum); //preallocates matrix for Term derivative columns
        Tdd0 = MatrixXd::Zero(df0.rows(), reqrdnum*(reqrdnum+1)/2); //preallocates matrix for Term second derivative columns
        TTerm = MatrixXd::Zero(Dose.rows(),Dose.cols());
        RdR = MatrixXd::Zero(df0.rows(), reqrdnum);
        RddR = MatrixXd::Zero(df0.rows(), reqrdnum*(reqrdnum+1)/2);
        fill(Ll.begin(), Ll.end(), 0.0);
        fill(Lld.begin(), Lld.end(), 0.0);
        fill(Lldd.begin(), Lldd.end(), 0.0);
        Rls1 =MatrixXd::Zero(ntime, 1); //precomputes a series of sums used frequently in the log-liklihood calculations
        Rls2 =MatrixXd::Zero(ntime, reqrdnum); //Many are repeated due to the same risk groups and derivatives being used at mulitple points
        Rls3 =MatrixXd::Zero(ntime, reqrdnum*(reqrdnum+1)/2);
        Lls1 =MatrixXd::Zero(ntime, 1);
        Lls2 =MatrixXd::Zero(ntime, reqrdnum);
        Lls3 =MatrixXd::Zero(ntime, reqrdnum*(reqrdnum+1)/2);
        beta_p = beta_best;//
        beta_a = beta_best;//
        beta_c = beta_best;//
        abs_max = abs_max0;
        dose_abs_max = dose_abs_max0;
        iter_check=0;
        iter_stop = 0;
        halves=0;
        iteration=0;
        halves = 0; //number of half-steps taken
        ind0 = fir; //used for validations
        iteration=0; //iteration number
        //
        convgd = FALSE;
        iter_stop =0; //tracks if the iterations should be stopped for convergence
        iter_check=0; //signal to check for convergence
        //
        maxiter = maxiters[guess];
        a_n = a_ns.row(guess);
        for (int i=0;i<beta_0.size();i++){
            beta_0[i] = a_n[i];
        }
        if (verbose){
            Rcout << "starting subterms " << term_tot << " in guess " << guess << endl;
        }
        //
        // Calculates the subterm and term values
        Make_subterms( totalnum, Term_n, tform, dfc, fir, T0, Td0, Tdd0, Dose, nonDose, TTerm,  nonDose_LIN, nonDose_PLIN, nonDose_LOGLIN ,beta_0, df0,dint,dslp,nthreads, debugging,KeepConstant);
        // ---------------------------------------------------------
        // Prints off a series of calculations to check at what point values are changing
        // ---------------------------------------------------------
        //
        //
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
            for (int ijk=0;ijk<reqrdnum;ijk++){
                Rcout << Td0.col(ijk).sum() << " ";
            }
            Rcout << " " << endl;
            Rcout << "second derivs checked ";
            for (int ijk=0;ijk<reqrdnum;ijk++){
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
        Make_Risks(modelform, tform, Term_n, totalnum, fir, T0, Td0, Tdd0, Te, R, Rd, Rdd, Dose, nonDose, TTerm, nonDose_LIN, nonDose_PLIN, nonDose_LOGLIN, RdR, RddR, nthreads, debugging,KeepConstant);
        //
        // Removes infinite values
        RdR = (RdR.array().isFinite()).select(RdR,0);
        RddR = (RddR.array().isFinite()).select(RddR,0);
        //
        if (verbose){
            Rcout << "risk checked ";
            for (int ijk=0;ijk<1;ijk++){
                Rcout << R.col(0).sum() << " ";
            }
            Rcout << " " << endl;
            Rcout << "risk1 checked ";
            for (int ijk=0;ijk<reqrdnum;ijk++){
                Rcout << Rd.col(ijk).sum() << " ";
            }
            Rcout << " " << endl;
            Rcout << "risk2 checked ";
            for (int ijk=0;ijk<reqrdnum;ijk++){
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
        if (verbose){
            Rcout << "Made Risk Side Lists" << endl;
        }
        // Calculates the side sum terms used
        Calculate_Sides( RiskFail, RiskGroup, totalnum, ntime, R, Rd, Rdd, Rls1, Rls2, Rls3, Lls1, Lls2, Lls3,nthreads, debugging,KeepConstant);
        //
        if (verbose){
            end_point = system_clock::now();
            ending = time_point_cast<microseconds>(end_point).time_since_epoch().count();
            Rcout <<"df100 "<<(ending-start)<<" "<<0<<" "<<0<<" "<<-1<<",Prep_Sides"<<endl;
            gibtime = system_clock::to_time_t(system_clock::now());
            Rcout << ctime(&gibtime) << endl;
        }
        //
        if (verbose){
            Rcout << "riskr checked ";
            for (int ijk=0;ijk<reqrdnum;ijk++){
                Rcout << Rls1.col(0).sum() << " ";
            }
            Rcout << " " << endl;
            Rcout << "risk1r checked ";
            for (int ijk=0;ijk<reqrdnum;ijk++){
                Rcout << Rls2.col(ijk).sum() << " ";
            }
            Rcout << " " << endl;
            Rcout << "risk2r checked ";
            for (int ijk=0;ijk<reqrdnum;ijk++){
                Rcout << Rls3.col(ijk*(ijk+1)/2+ijk).sum() << " ";
            }
            Rcout << " " << endl;
            //
            Rcout << "riskl checked ";
            for (int ijk=0;ijk<reqrdnum;ijk++){
                Rcout << Lls1.col(0).sum() << " ";
            }
            Rcout << " " << endl;
            Rcout << "risk1l checked ";
            for (int ijk=0;ijk<reqrdnum;ijk++){
                Rcout << Lls2.col(ijk).sum() << " ";
            }
            Rcout << " " << endl;
            Rcout << "risk2l checked ";
            for (int ijk=0;ijk<reqrdnum;ijk++){
                Rcout << Lls3.col(ijk*(ijk+1)/2+ijk).sum() << " ";
            }
            Rcout << " " << endl;
        }
        // Calculates log-likelihood
        Calc_LogLik( nthreads, RiskFail, RiskGroup, totalnum, ntime, R, Rd, Rdd,RdR,RddR, Rls1, Rls2, Rls3, Lls1, Lls2, Lls3, Ll, Lld, Lldd, debugging, ties_method,KeepConstant);
        //
        if (verbose){
            end_point = system_clock::now();
            ending = time_point_cast<microseconds>(end_point).time_since_epoch().count();
            Rcout <<"df100 "<<(ending-start)<<" "<<0<<" "<<0<<" "<<0<<",Calc"<<endl;//prints the time
            gibtime = system_clock::to_time_t(system_clock::now());
            Rcout << ctime(&gibtime) << endl;
            Rcout << "df101 ";//prints the log-likelihoods
            for (int ij=0;ij<reqrdnum;ij++){
                Rcout << Ll[ij] << " ";
            }
            Rcout << " " << endl;
            Rcout << "df102 ";//prints the first derivatives
            for (int ij=0;ij<reqrdnum;ij++){
                Rcout << Lld[ij] << " ";
            }
            Rcout << " " << endl;
            Rcout << "df103 ";//prints the second derivatives
            for (int ij=0;ij<reqrdnum;ij++){
                Rcout << Lldd[ij*reqrdnum+ij] << " ";
            }
            Lld_worst=0;
            for (int ij=0;ij<reqrdnum;ij++){//locates highest magnitude derivative
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
            for (int ij=0;ij<reqrdnum;ij++){//prints the newton step value for zero derivative
                Rcout << Lld[ij]/Lldd[ij*reqrdnum+ij] << " ";
            }
            Rcout << " " << endl;
            Rcout << "df106 ";
            for (int ij=0;ij<reqrdnum;ij++){//prints the newton step value for zero log-likelihood
                Rcout << Ll[ij]/Lld[ij] << " ";
            }
            Rcout << " " << endl;
            Rcout << "df107 " << double_step << " " << abs_max << " " << dose_abs_max << " " << Ll_comp[0] << " " << Ll_comp[1] << endl;//prints several convergence terms
        }
        //
        for (int i=0;i<beta_0.size();i++){
            beta_c[i] = beta_0[i];
        }
        while ((iteration < maxiter)&&(iter_stop==0)){
            iteration++;
            beta_p = beta_c;//
            beta_a = beta_c;//
            beta_best = beta_c;//
            //
            // Calcualtes the initial change in parameter
            Calc_Change( double_step, nthreads, totalnum, fir, der_iden, dbeta_cap, dose_abs_max, lr, abs_max, Ll, Lld, Lldd, dbeta, change_all, tform, dint,dslp, KeepConstant, debugging);
            Intercept_Bound(nthreads, totalnum, beta_0, dbeta, dfc, df0, KeepConstant, debugging, tform);
            if (verbose){
                Rcout << "Starting Halves"<<endl;//prints the final changes for validation
            }
            //
            //
            if ((Ll_abs_best>0)||(Ll_abs_best < Ll[ind0])){
                Ll_abs_best = Ll[ind0];
                beta_abs_best = beta_c;
                guess_abs_best=guess;
            }
            //
            if (verbose){
                Rcout << "df501 " << Ll_abs_best << endl;
                Rcout << "df504 ";//prints parameter values
                for (int ij=0;ij<totalnum;ij++){
                    Rcout << beta_abs_best[ij] << " ";
                }
                Rcout << " " << endl;
            }
            halves=0;
            while ((Ll[ind0] <= Ll_abs_best)&&(halves<halfmax)){ //repeats until half-steps maxed or an improvement
                //Refreshes the matrices used
                Dose = MatrixXd::Zero(df0.rows(),term_tot);
                nonDose = MatrixXd::Constant(df0.rows(),term_tot,0.0);
                nonDose_LIN = MatrixXd::Constant(df0.rows(),term_tot,0.0);
                nonDose_PLIN = MatrixXd::Constant(df0.rows(),term_tot,1.0);
                nonDose_LOGLIN = MatrixXd::Constant(df0.rows(),term_tot,1.0);
                T0 = MatrixXd::Zero(df0.rows(), totalnum); //preallocates matrix for Term column
                Td0 = MatrixXd::Zero(df0.rows(), reqrdnum); //preallocates matrix for Term derivative columns
                Tdd0 = MatrixXd::Zero(df0.rows(), reqrdnum*(reqrdnum+1)/2); //preallocates matrix for Term second derivative columns
                for (int ij=0;ij<totalnum;ij++){
                    beta_0[ij] = beta_a[ij] + dbeta[ij];
                    beta_c[ij] = beta_0[ij];
                }
                // ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------//
                // The same subterm, risk, sides, and log-likelihood calculations are performed every half-step and iteration
                // ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------//
                Make_subterms( totalnum, Term_n, tform, dfc, fir, T0, Td0, Tdd0, Dose,  nonDose, TTerm,  nonDose_LIN, nonDose_PLIN, nonDose_LOGLIN ,beta_0, df0,dint,dslp,nthreads, debugging,KeepConstant);
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
                    for (int ijk=0;ijk<reqrdnum;ijk++){
                        Rcout << Td0.col(ijk).sum() << " ";
                    }
                    Rcout << " " << endl;
                    Rcout << "second derivs checked ";
                    for (int ijk=0;ijk<reqrdnum;ijk++){
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
                RdR = MatrixXd::Zero(df0.rows(), reqrdnum);
                RddR = MatrixXd::Zero(df0.rows(), reqrdnum*(reqrdnum+1)/2);
                //
                //
                Make_Risks(modelform, tform, Term_n, totalnum, fir, T0, Td0, Tdd0, Te, R, Rd, Rdd, Dose, nonDose, TTerm, nonDose_LIN, nonDose_PLIN, nonDose_LOGLIN, RdR, RddR, nthreads, debugging,KeepConstant);
                //
                if (R.minCoeff()<0){
                    #pragma omp parallel for num_threads(nthreads)
                    for (int ijk=0;ijk<totalnum;ijk++){
                        int tij = Term_n[ijk];
                        if (TTerm.col(tij).minCoeff()<=0){
                            dbeta[ijk] = dbeta[ijk] / 2.0;
                        } else if (isinf(TTerm.col(tij).maxCoeff())){
                            dbeta[ijk] = dbeta[ijk] / 2.0;
                        }
                    }
                    if (verbose){
                        Rcout << "A non-positive risk was detected: " << R.minCoeff() << endl;
                    }
                    halves+=0.2;
                } else {
                    halves++;
                    RdR = (RdR.array().isFinite()).select(RdR,0);
                    RddR = (RddR.array().isFinite()).select(RddR,0);
                    //
                    if (verbose){
                        Rcout << "risk checked ";
                        for (int ijk=0;ijk<1;ijk++){
                            Rcout << R.col(0).sum() << " ";
                        }
                        Rcout << " " << endl;
                        Rcout << "risk1 checked ";
                        for (int ijk=0;ijk<reqrdnum;ijk++){
                            Rcout << Rd.col(ijk).sum() << " ";
                        }
                        Rcout << " " << endl;
                        Rcout << "risk2 checked ";
                        for (int ijk=0;ijk<reqrdnum;ijk++){
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
                    Rls2 =MatrixXd::Zero(ntime, reqrdnum); //Many are repeated due to the same risk groups and derivatives being used at mulitple points
                    Rls3 =MatrixXd::Zero(ntime, reqrdnum*(reqrdnum+1)/2);
                    Lls1 =MatrixXd::Zero(ntime, 1);
                    Lls2 =MatrixXd::Zero(ntime, reqrdnum);
                    Lls3 =MatrixXd::Zero(ntime, reqrdnum*(reqrdnum+1)/2);
                    Calculate_Sides( RiskFail, RiskGroup, totalnum, ntime, R, Rd, Rdd, Rls1, Rls2, Rls3, Lls1, Lls2, Lls3,nthreads, debugging,KeepConstant);
                    //
                    
                    if (verbose){
                        end_point = system_clock::now();
                        ending = time_point_cast<microseconds>(end_point).time_since_epoch().count();
                        Rcout <<"df100 "<<(ending-start)<<" "<<halves<<" "<<iteration<<" "<<ind0<<",Update_Sides"<<endl;
                        gibtime = system_clock::to_time_t(system_clock::now());
                        Rcout << ctime(&gibtime) << endl;
                        Rcout << "riskr checked ";
                        for (int ijk=0;ijk<reqrdnum;ijk++){
                            Rcout << Rls1.col(0).sum() << " ";
                        }
                        Rcout << " " << endl;
                        Rcout << "risk1r checked ";
                        for (int ijk=0;ijk<reqrdnum;ijk++){
                            Rcout << Rls2.col(ijk).sum() << " ";
                        }
                        Rcout << " " << endl;
                        Rcout << "risk2r checked ";
                        for (int ijk=0;ijk<reqrdnum;ijk++){
                            Rcout << Rls3.col(ijk*(ijk+1)/2+ijk).sum() << " ";
                        }
                        Rcout << " " << endl;
                        //
                        //
                        Rcout << "riskl checked ";
                        for (int ijk=0;ijk<reqrdnum;ijk++){
                            Rcout << Lls1.col(0).sum() << " ";
                        }
                        Rcout << " " << endl;
                        Rcout << "risk1l checked ";
                        for (int ijk=0;ijk<reqrdnum;ijk++){
                            Rcout << Lls2.col(ijk).sum() << " ";
                        }
                        Rcout << " " << endl;
                        Rcout << "risk2l checked ";
                        for (int ijk=0;ijk<reqrdnum;ijk++){
                            Rcout << Lls3.col(ijk*(ijk+1)/2+ijk).sum() << " ";
                        }
                        Rcout << " " << endl;
                    }
                    //
                    Calc_LogLik( nthreads, RiskFail, RiskGroup, totalnum, ntime, R, Rd, Rdd,RdR,RddR, Rls1, Rls2, Rls3, Lls1, Lls2, Lls3, Ll, Lld, Lldd, debugging, ties_method,KeepConstant);
                    if (verbose){
                        end_point = system_clock::now();
                        ending = time_point_cast<microseconds>(end_point).time_since_epoch().count();
                        Rcout <<"df100 "<<(ending-start)<<" "<<halves<<" "<<iteration<<" "<<ind0<<", Half"<<endl;
                        gibtime = system_clock::to_time_t(system_clock::now());
                        Rcout << ctime(&gibtime) << endl;
                        Rcout << "df101 ";
                        for (int ij=0;ij<reqrdnum;ij++){
                            Rcout << Ll[ij] << " ";
                        }
                        Rcout << " " << endl;
                        Rcout << "df102 ";
                        for (int ij=0;ij<reqrdnum;ij++){
                            Rcout << Lld[ij] << " ";
                        }
                        Rcout << " " << endl;
                        Rcout << "df103 ";
                        for (int ij=0;ij<reqrdnum;ij++){
                            Rcout << Lldd[ij*reqrdnum+ij] << " ";
                        }
                        Rcout << " " << endl;
                        Rcout << "df104 ";
                        for (int ij=0;ij<totalnum;ij++){
                            Rcout << beta_c[ij] << " ";
                        }
                        Rcout << " " << endl;
                        Rcout << "df105 ";
                        for (int ij=0;ij<reqrdnum;ij++){//prints the newton step value for zero derivative
                            Rcout << Lld[ij]/Lldd[ij*reqrdnum+ij] << " ";
                        }
                        Rcout << " " << endl;
                        Rcout << "df106 ";
                        for (int ij=0;ij<reqrdnum;ij++){//prints the newton step value for zero log-likelihood
                            Rcout << Ll[ij]/Lld[ij] << " ";
                        }
                        Rcout << " " << endl;
                        Rcout << "df107 " << double_step << " " << abs_max << " " << dose_abs_max << " " << Ll_comp[0] << " " << Ll_comp[1] << endl;//prints several convergence terms
                    }
                    //
                    if (change_all){ //If every covariate is to be changed
                        if (Ll[ind0] <= Ll_abs_best){//takes a half-step if needed
                            #pragma omp parallel for num_threads(nthreads)
                            for (int ijk=0;ijk<totalnum;ijk++){
                                dbeta[ijk] = dbeta[ijk] * 0.5; //
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
                    }
                    #pragma omp parallel for num_threads(nthreads)
                    for (int ijk=0;ijk<totalnum;ijk++){//totalnum*(totalnum+1)/2
                        beta_0[ijk] = beta_c[ijk];
                    }
                }
            }
            if (beta_best!=beta_c){//if the risk matrices aren't the optimal values, then they must be recalculated
                if (verbose){
                    Rcout << "Changing back to best"<<endl;
                }
                // If it goes through every half step without improvement, then the maximum change needs to be decreased
                abs_max = abs_max*pow(0.5,halfmax); // reduces the step sizes
                dose_abs_max = dose_abs_max*pow(0.5,halfmax);
                iter_check = 1;
                //
                beta_p = beta_best;//
                beta_a = beta_best;//
                beta_c = beta_best;//
                Dose = MatrixXd::Zero(df0.rows(),term_tot);
                nonDose = MatrixXd::Constant(df0.rows(),term_tot,0.0);
                nonDose_LIN = MatrixXd::Constant(df0.rows(),term_tot,0.0);
                nonDose_PLIN = MatrixXd::Constant(df0.rows(),term_tot,1.0);
                nonDose_LOGLIN = MatrixXd::Constant(df0.rows(),term_tot,1.0);
                for (int ij=0;ij<totalnum;ij++){
                    beta_0[ij] = beta_best[ij];
                }
                Make_subterms( totalnum, Term_n, tform, dfc, fir, T0, Td0, Tdd0, Dose, nonDose, TTerm,  nonDose_LIN, nonDose_PLIN, nonDose_LOGLIN,beta_0, df0,dint,dslp,nthreads, debugging, KeepConstant);
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
                    for (int ijk=0;ijk<reqrdnum;ijk++){
                        Rcout << Td0.col(ijk).sum() << " ";
                    }
                    Rcout << " " << endl;
                    Rcout << "second derivs checked ";
                    for (int ijk=0;ijk<reqrdnum;ijk++){
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
                RdR = MatrixXd::Zero(df0.rows(), reqrdnum);
                RddR = MatrixXd::Zero(df0.rows(), reqrdnum*(reqrdnum+1)/2);
                //
                //
                Make_Risks(modelform, tform, Term_n, totalnum, fir, T0, Td0, Tdd0, Te, R, Rd, Rdd, Dose, nonDose, TTerm, nonDose_LIN, nonDose_PLIN, nonDose_LOGLIN, RdR, RddR, nthreads, debugging,KeepConstant);
                R = (R.array().isFinite()).select(R,0);
                Rd = (Rd.array().isFinite()).select(Rd,0);
                Rdd = (Rdd.array().isFinite()).select(Rdd,0);
                RdR = (RdR.array().isFinite()).select(RdR,0);
                RddR = (RddR.array().isFinite()).select(RddR,0);
                //
                temp_list = List::create(_["betas"]=wrap(beta_0),_["Status"]="FAILED"); //used as a dummy return value for code checking
                if (verbose){
                    Rcout << "risk checked ";
                    for (int ijk=0;ijk<1;ijk++){
                        Rcout << R.col(0).sum() << " ";
                    }
                    Rcout << " " << endl;
                    Rcout << "risk1 checked ";
                    for (int ijk=0;ijk<reqrdnum;ijk++){
                        Rcout << Rd.col(ijk).sum() << " ";
                    }
                    Rcout << " " << endl;
                    Rcout << "risk2 checked ";
                    for (int ijk=0;ijk<reqrdnum;ijk++){
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
                Rls2 =MatrixXd::Zero(ntime, reqrdnum); //Many are repeated due to the same risk groups and derivatives being used at mulitple points
                Rls3 =MatrixXd::Zero(ntime, reqrdnum*(reqrdnum+1)/2);
                Lls1 =MatrixXd::Zero(ntime, 1);
                Lls2 =MatrixXd::Zero(ntime, reqrdnum);
                Lls3 =MatrixXd::Zero(ntime, reqrdnum*(reqrdnum+1)/2);
                Calculate_Sides( RiskFail, RiskGroup, totalnum, ntime, R, Rd, Rdd, Rls1, Rls2, Rls3, Lls1, Lls2, Lls3,nthreads, debugging,KeepConstant);
                //
                if (verbose){
                    end_point = system_clock::now();
                    ending = time_point_cast<microseconds>(end_point).time_since_epoch().count();
                    Rcout <<"df100 "<<(ending-start)<<" "<<halves<<" "<<iteration<<" "<<ind0<<",Update_Sides"<<endl;
                    gibtime = system_clock::to_time_t(system_clock::now());
                    Rcout << ctime(&gibtime) << endl;
                    Rcout << "riskr checked ";
                    for (int ijk=0;ijk<reqrdnum;ijk++){
                        Rcout << Rls1.col(0).sum() << " ";
                    }
                    Rcout << " " << endl;
                    Rcout << "risk1r checked ";
                    for (int ijk=0;ijk<reqrdnum;ijk++){
                        Rcout << Rls2.col(ijk).sum() << " ";
                    }
                    Rcout << " " << endl;
                    Rcout << "risk2r checked ";
                    for (int ijk=0;ijk<reqrdnum;ijk++){
                        Rcout << Rls3.col(ijk*(ijk+1)/2+ijk).sum() << " ";
                    }
                    Rcout << " " << endl;
                    //
                    //
                    Rcout << "riskl checked ";
                    for (int ijk=0;ijk<reqrdnum;ijk++){
                        Rcout << Lls1.col(0).sum() << " ";
                    }
                    Rcout << " " << endl;
                    Rcout << "risk1l checked ";
                    for (int ijk=0;ijk<reqrdnum;ijk++){
                        Rcout << Lls2.col(ijk).sum() << " ";
                    }
                    Rcout << " " << endl;
                    Rcout << "risk2l checked ";
                    for (int ijk=0;ijk<reqrdnum;ijk++){
                        Rcout << Lls3.col(ijk*(ijk+1)/2+ijk).sum() << " ";
                    }
                    Rcout << " " << endl;
                }
                //
                Calc_LogLik( nthreads, RiskFail, RiskGroup, totalnum, ntime, R, Rd, Rdd,RdR,RddR, Rls1, Rls2, Rls3, Lls1, Lls2, Lls3, Ll, Lld, Lldd, debugging, ties_method,KeepConstant);
            }
            Lld_worst=0;
            for (int ij=0;ij<reqrdnum;ij++){
                if (abs(Lld[ij]) > Lld_worst){
                    Lld_worst = abs(Lld[ij]);
                }
            }
            if (iteration > reqrdnum){//Doesn't check the first several iterations for convergence
                if ((iteration % (reqrdnum))||(iter_check==1)){//Checks every set number of iterations
                    iter_check=0;
                    if (Lld_worst < deriv_epsilon){//ends if the derivatives are low enough
                        iter_stop = 1;
                        convgd = TRUE;
                    }
                    Ll_comp[1]=Ll[0];
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
                for (int ij=0;ij<reqrdnum;ij++){
                    Rcout << Ll[ij] << " ";
                }
                Rcout << " " << endl;
                Rcout << "df102 ";
                for (int ij=0;ij<reqrdnum;ij++){
                    Rcout << Lld[ij] << " ";
                }
                Rcout << " " << endl;
                Rcout << "df103 ";
                for (int ij=0;ij<reqrdnum;ij++){
                    Rcout << Lldd[ij*reqrdnum+ij] << " ";
                }
                Rcout << " " << endl;
                Rcout << "df104 ";
                for (int ij=0;ij<totalnum;ij++){
                    Rcout << beta_c[ij] << " ";
                }
                Rcout << " " << endl;
                Rcout << "df105 ";
                for (int ij=0;ij<reqrdnum;ij++){//prints the newton step value for zero derivative
                    Rcout << Lld[ij]/Lldd[ij*reqrdnum+ij] << " ";
                }
                Rcout << " " << endl;
                Rcout << "df106 ";
                for (int ij=0;ij<reqrdnum;ij++){//prints the newton step value for zero log-likelihood
                    Rcout << Ll[ij]/Lld[ij] << " ";
                }
                Rcout << " " << endl;
                Rcout << "df107 " << double_step << " " << abs_max << " " << dose_abs_max << " " << Ll_comp[0] << " " << Ll_comp[1] << endl;//prints several convergence terms
                Rcout << "Finshed iteration" << endl;
            }
        }
        // -----------------------------------------------
        // Performing Full Calculation to get full second derivative matrix
        // -----------------------------------------------
        fill(Ll.begin(), Ll.end(), 0.0);
        fill(Lld.begin(), Lld.end(), 0.0);
        fill(Lldd.begin(), Lldd.end(), 0.0);
        Calculate_Sides( RiskFail, RiskGroup, totalnum, ntime, R, Rd, Rdd, Rls1, Rls2, Rls3, Lls1, Lls2, Lls3,nthreads, debugging,KeepConstant);
        //
        if (verbose){
            end_point = system_clock::now();
            ending = time_point_cast<microseconds>(end_point).time_since_epoch().count();
            Rcout <<"df100 "<<(ending-start)<<" "<<halves<<" "<<iteration<<" "<<ind0<<",Update_Sides"<<endl;
            gibtime = system_clock::to_time_t(system_clock::now());
            Rcout << ctime(&gibtime) << endl;
            Rcout << "Wrapping up" << endl;
        }
        //
        Calc_LogLik( nthreads, RiskFail, RiskGroup, totalnum, ntime, R, Rd, Rdd,RdR,RddR, Rls1, Rls2, Rls3, Lls1, Lls2, Lls3, Ll, Lld, Lldd, debugging, ties_method,KeepConstant);
        //
        a_n = beta_0;
        //
        beta_fin(guess, _) = a_n;
        LL_fin[guess] = Ll[0];
        if ((Ll_abs_best>0)||(Ll_abs_best < Ll[ind0])){
            Ll_abs_best = Ll[ind0];
            beta_abs_best = beta_c;
            guess_abs_best=guess;
        }
        //
        if (verbose){
            Rcout << "df501 " << Ll_abs_best << endl;
            Rcout << "df504 ";//prints parameter values
            for (int ij=0;ij<totalnum;ij++){
                Rcout << beta_abs_best[ij] << " ";
            }
            Rcout << " " << endl;
        }
    }
    //
    Dose = MatrixXd::Zero(df0.rows(),term_tot);
    nonDose = MatrixXd::Constant(df0.rows(),term_tot,0.0);
    nonDose_LIN = MatrixXd::Constant(df0.rows(),term_tot,0.0);
    nonDose_PLIN = MatrixXd::Constant(df0.rows(),term_tot,1.0);
    nonDose_LOGLIN = MatrixXd::Constant(df0.rows(),term_tot,1.0);
    T0 = MatrixXd::Zero(df0.rows(), totalnum); //preallocates matrix for Term column
    Td0 = MatrixXd::Zero(df0.rows(), reqrdnum); //preallocates matrix for Term derivative columns
    Tdd0 = MatrixXd::Zero(df0.rows(), reqrdnum*(reqrdnum+1)/2); //preallocates matrix for Term second derivative columns
    TTerm = MatrixXd::Zero(Dose.rows(),Dose.cols());
    RdR = MatrixXd::Zero(df0.rows(), reqrdnum);
    RddR = MatrixXd::Zero(df0.rows(), reqrdnum*(reqrdnum+1)/2);
    fill(Ll.begin(), Ll.end(), 0.0);
    fill(Lld.begin(), Lld.end(), 0.0);
    fill(Lldd.begin(), Lldd.end(), 0.0);
    Rls1 =MatrixXd::Zero(ntime, 1); //precomputes a series of sums used frequently in the log-liklihood calculations
    Rls2 =MatrixXd::Zero(ntime, reqrdnum); //Many are repeated due to the same risk groups and derivatives being used at mulitple points
    Rls3 =MatrixXd::Zero(ntime, reqrdnum*(reqrdnum+1)/2);
    Lls1 =MatrixXd::Zero(ntime, 1);
    Lls2 =MatrixXd::Zero(ntime, reqrdnum);
    Lls3 =MatrixXd::Zero(ntime, reqrdnum*(reqrdnum+1)/2);
    beta_p = beta_best;//
    beta_a = beta_best;//
    beta_c = beta_best;//
    abs_max = abs_max0;
    dose_abs_max = dose_abs_max0;
    iter_check=0;
    iter_stop = 0;
    halves=0;
    iteration=0;
    halves = 0; //number of half-steps taken
    ind0 = fir; //used for validations
    iteration=0; //iteration number
    //
    convgd = FALSE;
    iter_stop =0; //tracks if the iterations should be stopped for convergence
    iter_check=0; //signal to check for convergence
    //
    int guess_max=guess_abs_best;
    if (verbose){
        Rcout << "Guess number, parameter values, Log-Likelihood" << endl;
        NumericVector beta_temp;
        for (int i=0; i<=guesses; i++){
            beta_temp = wrap(beta_fin.row(i));
            if (i==guess_max){
                Rcout << i << ", " << beta_temp << ", " << LL_fin[i] << "<-- Best Guess" << endl;
            } else {
                Rcout << i << ", " << beta_temp << ", " << LL_fin[i] << endl;
            }
        }
    }
    //
    maxiter = maxiters[guesses+1];
    a_n = beta_abs_best;
    for (int i=0;i<beta_0.size();i++){
        beta_0[i] = a_n[i];
    }
    for (int i=0;i<beta_0.size();i++){
        beta_c[i] = beta_0[i];
    }
    //
    if (verbose){
        Rcout << "starting subterms " << term_tot << " in best guess " << guess_max << endl;
    }
    //
    // Calculates the subterm and term values
    Make_subterms( totalnum, Term_n, tform, dfc, fir, T0, Td0, Tdd0, Dose, nonDose, TTerm,  nonDose_LIN, nonDose_PLIN, nonDose_LOGLIN ,beta_0, df0,dint,dslp,nthreads, debugging,KeepConstant);
    // ---------------------------------------------------------
    // Prints off a series of calculations to check at what point values are changing
    // ---------------------------------------------------------
    //
    //
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
        for (int ijk=0;ijk<reqrdnum;ijk++){
            Rcout << Td0.col(ijk).sum() << " ";
        }
        Rcout << " " << endl;
        Rcout << "second derivs checked ";
        for (int ijk=0;ijk<reqrdnum;ijk++){
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
    Make_Risks(modelform, tform, Term_n, totalnum, fir, T0, Td0, Tdd0, Te, R, Rd, Rdd, Dose, nonDose, TTerm, nonDose_LIN, nonDose_PLIN, nonDose_LOGLIN, RdR, RddR, nthreads, debugging,KeepConstant);
    //
    // Removes infinite values
    RdR = (RdR.array().isFinite()).select(RdR,0);
    RddR = (RddR.array().isFinite()).select(RddR,0);
    //
    //
    if (verbose){
        Rcout << "risk checked ";
        for (int ijk=0;ijk<1;ijk++){
            Rcout << R.col(0).sum() << " ";
        }
        Rcout << " " << endl;
        Rcout << "risk1 checked ";
        for (int ijk=0;ijk<reqrdnum;ijk++){
            Rcout << Rd.col(ijk).sum() << " ";
        }
        Rcout << " " << endl;
        Rcout << "risk2 checked ";
        for (int ijk=0;ijk<reqrdnum;ijk++){
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
    if (verbose){
        Rcout << "Made Risk Side Lists" << endl;
    }
    // Calculates the side sum terms used
    Calculate_Sides( RiskFail, RiskGroup, totalnum, ntime, R, Rd, Rdd, Rls1, Rls2, Rls3, Lls1, Lls2, Lls3,nthreads, debugging,KeepConstant);
    //
    if (verbose){
        end_point = system_clock::now();
        ending = time_point_cast<microseconds>(end_point).time_since_epoch().count();
        Rcout <<"df100 "<<(ending-start)<<" "<<0<<" "<<0<<" "<<-1<<",Prep_Sides"<<endl;
        gibtime = system_clock::to_time_t(system_clock::now());
        Rcout << ctime(&gibtime) << endl;
    }
    //
    if (verbose){
        Rcout << "riskr checked ";
        for (int ijk=0;ijk<reqrdnum;ijk++){
            Rcout << Rls1.col(0).sum() << " ";
        }
        Rcout << " " << endl;
        Rcout << "risk1r checked ";
        for (int ijk=0;ijk<reqrdnum;ijk++){
            Rcout << Rls2.col(ijk).sum() << " ";
        }
        Rcout << " " << endl;
        Rcout << "risk2r checked ";
        for (int ijk=0;ijk<reqrdnum;ijk++){
            Rcout << Rls3.col(ijk*(ijk+1)/2+ijk).sum() << " ";
        }
        Rcout << " " << endl;
        //
        Rcout << "riskl checked ";
        for (int ijk=0;ijk<reqrdnum;ijk++){
            Rcout << Lls1.col(0).sum() << " ";
        }
        Rcout << " " << endl;
        Rcout << "risk1l checked ";
        for (int ijk=0;ijk<reqrdnum;ijk++){
            Rcout << Lls2.col(ijk).sum() << " ";
        }
        Rcout << " " << endl;
        Rcout << "risk2l checked ";
        for (int ijk=0;ijk<reqrdnum;ijk++){
            Rcout << Lls3.col(ijk*(ijk+1)/2+ijk).sum() << " ";
        }
        Rcout << " " << endl;
    }
    // Calculates log-likelihood
    Calc_LogLik( nthreads, RiskFail, RiskGroup, totalnum, ntime, R, Rd, Rdd,RdR,RddR, Rls1, Rls2, Rls3, Lls1, Lls2, Lls3, Ll, Lld, Lldd, debugging, ties_method,KeepConstant);
    //
    //
    if (verbose){
        end_point = system_clock::now();
        ending = time_point_cast<microseconds>(end_point).time_since_epoch().count();
        Rcout <<"df100 "<<(ending-start)<<" "<<0<<" "<<0<<" "<<0<<",Calc"<<endl;//prints the time
        gibtime = system_clock::to_time_t(system_clock::now());
        Rcout << ctime(&gibtime) << endl;
        Rcout << "df101 ";//prints the log-likelihoods
        for (int ij=0;ij<reqrdnum;ij++){
            Rcout << Ll[ij] << " ";
        }
        Rcout << " " << endl;
        Rcout << "df102 ";//prints the first derivatives
        for (int ij=0;ij<reqrdnum;ij++){
            Rcout << Lld[ij] << " ";
        }
        Rcout << " " << endl;
        Rcout << "df103 ";//prints the second derivatives
        for (int ij=0;ij<reqrdnum;ij++){
            Rcout << Lldd[ij*reqrdnum+ij] << " ";
        }
        Lld_worst=0;
        for (int ij=0;ij<reqrdnum;ij++){//locates highest magnitude derivative
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
        for (int ij=0;ij<reqrdnum;ij++){//prints the newton step value for zero derivative
            Rcout << Lld[ij]/Lldd[ij*reqrdnum+ij] << " ";
        }
        Rcout << " " << endl;
        Rcout << "df106 ";
        for (int ij=0;ij<reqrdnum;ij++){//prints the newton step value for zero log-likelihood
            Rcout << Ll[ij]/Lld[ij] << " ";
        }
        Rcout << " " << endl;
        Rcout << "df107 " << double_step << " " << abs_max << " " << dose_abs_max << " " << Ll_comp[0] << " " << Ll_comp[1] << endl;//prints several convergence terms
    }
    //
    while ((iteration < maxiter)&&(iter_stop==0)){
        iteration++;
        beta_p = beta_c;//
        beta_a = beta_c;//
        beta_best = beta_c;//
        //
        // Calculates the initial change in parameter
        Calc_Change( double_step, nthreads, totalnum, fir, der_iden, dbeta_cap, dose_abs_max, lr, abs_max, Ll, Lld, Lldd, dbeta, change_all, tform, dint,dslp, KeepConstant, debugging);
        Intercept_Bound(nthreads, totalnum, beta_0, dbeta, dfc, df0, KeepConstant, debugging, tform);
        if (verbose){
            Rcout << "Starting Halves"<<endl;//prints the final changes for validation
        }
        //
        //
        if ((Ll_abs_best>0)||(Ll_abs_best < Ll[ind0])){
            Ll_abs_best = Ll[ind0];
            beta_abs_best = beta_c;
        }
        //
        if (verbose){
            Rcout << "df501 " << Ll_abs_best << endl;
            Rcout << "df504 ";//prints parameter values
            for (int ij=0;ij<totalnum;ij++){
                Rcout << beta_abs_best[ij] << " ";
            }
            Rcout << " " << endl;
        }
        halves=0;
        while ((Ll[ind0] <= Ll_abs_best)&&(halves<halfmax)){ //repeats until half-steps maxed or an improvement
            //Refreshes the matrices used
            Dose = MatrixXd::Zero(df0.rows(),term_tot);
            nonDose = MatrixXd::Constant(df0.rows(),term_tot,0.0);
            nonDose_LIN = MatrixXd::Constant(df0.rows(),term_tot,0.0);
            nonDose_PLIN = MatrixXd::Constant(df0.rows(),term_tot,1.0);
            nonDose_LOGLIN = MatrixXd::Constant(df0.rows(),term_tot,1.0);
            T0 = MatrixXd::Zero(df0.rows(), totalnum); //preallocates matrix for Term column
            Td0 = MatrixXd::Zero(df0.rows(), reqrdnum); //preallocates matrix for Term derivative columns
            Tdd0 = MatrixXd::Zero(df0.rows(), reqrdnum*(reqrdnum+1)/2); //preallocates matrix for Term second derivative columns
            for (int ij=0;ij<totalnum;ij++){
                beta_0[ij] = beta_a[ij] + dbeta[ij];
                beta_c[ij] = beta_0[ij];
            }
            // ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------//
            // The same subterm, risk, sides, and log-likelihood calculations are performed every half-step and iteration
            // ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------//
            Make_subterms( totalnum, Term_n, tform, dfc, fir, T0, Td0, Tdd0, Dose,  nonDose, TTerm,  nonDose_LIN, nonDose_PLIN, nonDose_LOGLIN ,beta_0, df0,dint,dslp,nthreads, debugging,KeepConstant);
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
                for (int ijk=0;ijk<reqrdnum;ijk++){
                    Rcout << Td0.col(ijk).sum() << " ";
                }
                Rcout << " " << endl;
                Rcout << "second derivs checked ";
                for (int ijk=0;ijk<reqrdnum;ijk++){
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
            RdR = MatrixXd::Zero(df0.rows(), reqrdnum);
            RddR = MatrixXd::Zero(df0.rows(), reqrdnum*(reqrdnum+1)/2);
            //
            //
            Make_Risks(modelform, tform, Term_n, totalnum, fir, T0, Td0, Tdd0, Te, R, Rd, Rdd, Dose, nonDose, TTerm, nonDose_LIN, nonDose_PLIN, nonDose_LOGLIN, RdR, RddR, nthreads, debugging,KeepConstant);
            //
            if (R.minCoeff()<0){
                #pragma omp parallel for num_threads(nthreads)
                for (int ijk=0;ijk<totalnum;ijk++){
                    int tij = Term_n[ijk];
                    if (TTerm.col(tij).minCoeff()<=0){
                        dbeta[ijk] = dbeta[ijk] / 2.0;
                    } else if (isinf(TTerm.col(tij).maxCoeff())){
                        dbeta[ijk] = dbeta[ijk] / 2.0;
                    }
                }
                if (verbose){
                    Rcout << "A non-positive risk was detected: " << R.minCoeff() << endl;
                }
                halves+=0.2;
            } else {
                halves++;
                RdR = (RdR.array().isFinite()).select(RdR,0);
                RddR = (RddR.array().isFinite()).select(RddR,0);
                //
                if (verbose){
                    Rcout << "risk checked ";
                    for (int ijk=0;ijk<1;ijk++){
                        Rcout << R.col(0).sum() << " ";
                    }
                    Rcout << " " << endl;
                    Rcout << "risk1 checked ";
                    for (int ijk=0;ijk<reqrdnum;ijk++){
                        Rcout << Rd.col(ijk).sum() << " ";
                    }
                    Rcout << " " << endl;
                    Rcout << "risk2 checked ";
                    for (int ijk=0;ijk<reqrdnum;ijk++){
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
                Rls2 =MatrixXd::Zero(ntime, reqrdnum); //Many are repeated due to the same risk groups and derivatives being used at mulitple points
                Rls3 =MatrixXd::Zero(ntime, reqrdnum*(reqrdnum+1)/2);
                Lls1 =MatrixXd::Zero(ntime, 1);
                Lls2 =MatrixXd::Zero(ntime, reqrdnum);
                Lls3 =MatrixXd::Zero(ntime, reqrdnum*(reqrdnum+1)/2);
                Calculate_Sides( RiskFail, RiskGroup, totalnum, ntime, R, Rd, Rdd, Rls1, Rls2, Rls3, Lls1, Lls2, Lls3,nthreads, debugging,KeepConstant);
                //
                
                if (verbose){
                    end_point = system_clock::now();
                    ending = time_point_cast<microseconds>(end_point).time_since_epoch().count();
                    Rcout <<"df100 "<<(ending-start)<<" "<<halves<<" "<<iteration<<" "<<ind0<<",Update_Sides"<<endl;
                    gibtime = system_clock::to_time_t(system_clock::now());
                    Rcout << ctime(&gibtime) << endl;
                    Rcout << "riskr checked ";
                    for (int ijk=0;ijk<reqrdnum;ijk++){
                        Rcout << Rls1.col(0).sum() << " ";
                    }
                    Rcout << " " << endl;
                    Rcout << "risk1r checked ";
                    for (int ijk=0;ijk<reqrdnum;ijk++){
                        Rcout << Rls2.col(ijk).sum() << " ";
                    }
                    Rcout << " " << endl;
                    Rcout << "risk2r checked ";
                    for (int ijk=0;ijk<reqrdnum;ijk++){
                        Rcout << Rls3.col(ijk*(ijk+1)/2+ijk).sum() << " ";
                    }
                    Rcout << " " << endl;
                    //
                    //
                    Rcout << "riskl checked ";
                    for (int ijk=0;ijk<reqrdnum;ijk++){
                        Rcout << Lls1.col(0).sum() << " ";
                    }
                    Rcout << " " << endl;
                    Rcout << "risk1l checked ";
                    for (int ijk=0;ijk<reqrdnum;ijk++){
                        Rcout << Lls2.col(ijk).sum() << " ";
                    }
                    Rcout << " " << endl;
                    Rcout << "risk2l checked ";
                    for (int ijk=0;ijk<reqrdnum;ijk++){
                        Rcout << Lls3.col(ijk*(ijk+1)/2+ijk).sum() << " ";
                    }
                    Rcout << " " << endl;
                }
                //
                Calc_LogLik( nthreads, RiskFail, RiskGroup, totalnum, ntime, R, Rd, Rdd,RdR,RddR, Rls1, Rls2, Rls3, Lls1, Lls2, Lls3, Ll, Lld, Lldd, debugging, ties_method,KeepConstant);
                if (verbose){
                    end_point = system_clock::now();
                    ending = time_point_cast<microseconds>(end_point).time_since_epoch().count();
                    Rcout <<"df100 "<<(ending-start)<<" "<<halves<<" "<<iteration<<" "<<ind0<<", Half"<<endl;
                    gibtime = system_clock::to_time_t(system_clock::now());
                    Rcout << ctime(&gibtime) << endl;
                    Rcout << "df101 ";
                    for (int ij=0;ij<reqrdnum;ij++){
                        Rcout << Ll[ij] << " ";
                    }
                    Rcout << " " << endl;
                    Rcout << "df102 ";
                    for (int ij=0;ij<reqrdnum;ij++){
                        Rcout << Lld[ij] << " ";
                    }
                    Rcout << " " << endl;
                    Rcout << "df103 ";
                    for (int ij=0;ij<reqrdnum;ij++){
                        Rcout << Lldd[ij*reqrdnum+ij] << " ";
                    }
                    Rcout << " " << endl;
                    Rcout << "df104 ";
                    for (int ij=0;ij<totalnum;ij++){
                        Rcout << beta_c[ij] << " ";
                    }
                    Rcout << " " << endl;
                    Rcout << "df105 ";
                    for (int ij=0;ij<reqrdnum;ij++){//prints the newton step value for zero derivative
                        Rcout << Lld[ij]/Lldd[ij*reqrdnum+ij] << " ";
                    }
                    Rcout << " " << endl;
                    Rcout << "df106 ";
                    for (int ij=0;ij<reqrdnum;ij++){//prints the newton step value for zero log-likelihood
                        Rcout << Ll[ij]/Lld[ij] << " ";
                    }
                    Rcout << " " << endl;
                    Rcout << "df107 " << double_step << " " << abs_max << " " << dose_abs_max << " " << Ll_comp[0] << " " << Ll_comp[1] << endl;//prints several convergence terms
                }
                //
                if (change_all){ //If every covariate is to be changed
                    if (Ll[ind0] <= Ll_abs_best){//takes a half-step if needed
                        #pragma omp parallel for num_threads(nthreads)
                        for (int ijk=0;ijk<totalnum;ijk++){
                            dbeta[ijk] = dbeta[ijk] * 0.5; //
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
                }
                #pragma omp parallel for num_threads(nthreads)
                for (int ijk=0;ijk<totalnum;ijk++){//totalnum*(totalnum+1)/2
                    beta_0[ijk] = beta_c[ijk];
                }
            }
        }
        if (beta_best!=beta_c){//if the risk matrices aren't the optimal values, then they must be recalculated
            if (verbose){
                Rcout << "Changing back to best"<<endl;
            }
            // If it goes through every half step without improvement, then the maximum change needs to be decreased
            abs_max = abs_max*pow(0.5,halfmax); // reduces the step sizes
            dose_abs_max = dose_abs_max*pow(0.5,halfmax);
            iter_check = 1;
            //
            beta_p = beta_best;//
            beta_a = beta_best;//
            beta_c = beta_best;//
            Dose = MatrixXd::Zero(df0.rows(),term_tot);
            nonDose = MatrixXd::Constant(df0.rows(),term_tot,0.0);
            nonDose_LIN = MatrixXd::Constant(df0.rows(),term_tot,0.0);
            nonDose_PLIN = MatrixXd::Constant(df0.rows(),term_tot,1.0);
            nonDose_LOGLIN = MatrixXd::Constant(df0.rows(),term_tot,1.0);
            for (int ij=0;ij<totalnum;ij++){
                beta_0[ij] = beta_best[ij];
            }
            Make_subterms( totalnum, Term_n, tform, dfc, fir, T0, Td0, Tdd0, Dose, nonDose, TTerm,  nonDose_LIN, nonDose_PLIN, nonDose_LOGLIN,beta_0, df0,dint,dslp,nthreads, debugging, KeepConstant);
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
                for (int ijk=0;ijk<reqrdnum;ijk++){
                    Rcout << Td0.col(ijk).sum() << " ";
                }
                Rcout << " " << endl;
                Rcout << "second derivs checked ";
                for (int ijk=0;ijk<reqrdnum;ijk++){
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
            RdR = MatrixXd::Zero(df0.rows(), reqrdnum);
            RddR = MatrixXd::Zero(df0.rows(), reqrdnum*(reqrdnum+1)/2);
            //
            //
            Make_Risks(modelform, tform, Term_n, totalnum, fir, T0, Td0, Tdd0, Te, R, Rd, Rdd, Dose, nonDose, TTerm, nonDose_LIN, nonDose_PLIN, nonDose_LOGLIN, RdR, RddR, nthreads, debugging,KeepConstant);
            R = (R.array().isFinite()).select(R,0);
            Rd = (Rd.array().isFinite()).select(Rd,0);
            Rdd = (Rdd.array().isFinite()).select(Rdd,0);
            RdR = (RdR.array().isFinite()).select(RdR,0);
            RddR = (RddR.array().isFinite()).select(RddR,0);
            //
            temp_list = List::create(_["betas"]=wrap(beta_0),_["Status"]="FAILED"); //used as a dummy return value for code checking
            if (verbose){
                Rcout << "risk checked ";
                for (int ijk=0;ijk<1;ijk++){
                    Rcout << R.col(0).sum() << " ";
                }
                Rcout << " " << endl;
                Rcout << "risk1 checked ";
                for (int ijk=0;ijk<reqrdnum;ijk++){
                    Rcout << Rd.col(ijk).sum() << " ";
                }
                Rcout << " " << endl;
                Rcout << "risk2 checked ";
                for (int ijk=0;ijk<reqrdnum;ijk++){
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
            Rls2 =MatrixXd::Zero(ntime, reqrdnum); //Many are repeated due to the same risk groups and derivatives being used at mulitple points
            Rls3 =MatrixXd::Zero(ntime, reqrdnum*(reqrdnum+1)/2);
            Lls1 =MatrixXd::Zero(ntime, 1);
            Lls2 =MatrixXd::Zero(ntime, reqrdnum);
            Lls3 =MatrixXd::Zero(ntime, reqrdnum*(reqrdnum+1)/2);
            Calculate_Sides( RiskFail, RiskGroup, totalnum, ntime, R, Rd, Rdd, Rls1, Rls2, Rls3, Lls1, Lls2, Lls3,nthreads, debugging,KeepConstant);
            //
            if (verbose){
                end_point = system_clock::now();
                ending = time_point_cast<microseconds>(end_point).time_since_epoch().count();
                Rcout <<"df100 "<<(ending-start)<<" "<<halves<<" "<<iteration<<" "<<ind0<<",Update_Sides"<<endl;
                gibtime = system_clock::to_time_t(system_clock::now());
                Rcout << ctime(&gibtime) << endl;
                Rcout << "riskr checked ";
                for (int ijk=0;ijk<1;ijk++){
                    Rcout << Rls1.col(0).sum() << " ";
                }
                Rcout << " " << endl;
                Rcout << "risk1r checked ";
                for (int ijk=0;ijk<reqrdnum;ijk++){
                    Rcout << Rls2.col(ijk).sum() << " ";
                }
                Rcout << " " << endl;
                Rcout << "risk2r checked ";
                for (int ijk=0;ijk<reqrdnum;ijk++){
                    Rcout << Rls3.col(ijk*(ijk+1)/2+ijk).sum() << " ";
                }
                Rcout << " " << endl;
                //
                //
                Rcout << "riskl checked ";
                for (int ijk=0;ijk<1;ijk++){
                    Rcout << Lls1.col(0).sum() << " ";
                }
                Rcout << " " << endl;
                Rcout << "risk1l checked ";
                for (int ijk=0;ijk<reqrdnum;ijk++){
                    Rcout << Lls2.col(ijk).sum() << " ";
                }
                Rcout << " " << endl;
                Rcout << "risk2l checked ";
                for (int ijk=0;ijk<reqrdnum;ijk++){
                    Rcout << Lls3.col(ijk*(ijk+1)/2+ijk).sum() << " ";
                }
                Rcout << " " << endl;
            }
            //
            Calc_LogLik( nthreads, RiskFail, RiskGroup, totalnum, ntime, R, Rd, Rdd,RdR,RddR, Rls1, Rls2, Rls3, Lls1, Lls2, Lls3, Ll, Lld, Lldd, debugging, ties_method,KeepConstant);
        }
        Lld_worst=0;
        for (int ij=0;ij<reqrdnum;ij++){
            if (abs(Lld[ij]) > Lld_worst){
                Lld_worst = abs(Lld[ij]);
            }
        }
        if (iteration > reqrdnum){//Doesn't check the first several iterations for convergence
            if ((iteration % (reqrdnum))||(iter_check==1)){//Checks every set number of iterations
                iter_check=0;
                if (Lld_worst < deriv_epsilon){//ends if the derivatives are low enough
                    iter_stop = 1;
                    convgd = TRUE;
                }
                Ll_comp[1]=Ll[0];
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
            for (int ij=0;ij<reqrdnum;ij++){
                Rcout << Ll[ij] << " ";
            }
            Rcout << " " << endl;
            Rcout << "df102 ";
            for (int ij=0;ij<reqrdnum;ij++){
                Rcout << Lld[ij] << " ";
            }
            Rcout << " " << endl;
            Rcout << "df103 ";
            for (int ij=0;ij<reqrdnum;ij++){
                Rcout << Lldd[ij*reqrdnum+ij] << " ";
            }
            Rcout << " " << endl;
            Rcout << "df104 ";
            for (int ij=0;ij<totalnum;ij++){
                Rcout << beta_c[ij] << " ";
            }
            Rcout << " " << endl;
            Rcout << "df105 ";
            for (int ij=0;ij<reqrdnum;ij++){//prints the newton step value for zero derivative
                Rcout << Lld[ij]/Lldd[ij*reqrdnum+ij] << " ";
            }
            Rcout << " " << endl;
            Rcout << "df106 ";
            for (int ij=0;ij<reqrdnum;ij++){//prints the newton step value for zero log-likelihood
                Rcout << Ll[ij]/Lld[ij] << " ";
            }
            Rcout << " " << endl;
            Rcout << "df107 " << double_step << " " << abs_max << " " << dose_abs_max << " " << Ll_comp[0] << " " << Ll_comp[1] << endl;//prints several convergence terms
            Rcout << "Finshed iteration" << endl;
        }
    }
    // -----------------------------------------------
    // Performing Full Calculation to get full second derivative matrix
    // -----------------------------------------------
    fill(Ll.begin(), Ll.end(), 0.0);
    fill(Lld.begin(), Lld.end(), 0.0);
    fill(Lldd.begin(), Lldd.end(), 0.0);
    Calculate_Sides( RiskFail, RiskGroup, totalnum, ntime, R, Rd, Rdd, Rls1, Rls2, Rls3, Lls1, Lls2, Lls3,nthreads, debugging,KeepConstant);
    //
    if (verbose){
        end_point = system_clock::now();
        ending = time_point_cast<microseconds>(end_point).time_since_epoch().count();
        Rcout <<"df100 "<<(ending-start)<<" "<<halves<<" "<<iteration<<" "<<ind0<<",Update_Sides"<<endl;
        gibtime = system_clock::to_time_t(system_clock::now());
        Rcout << ctime(&gibtime) << endl;
        Rcout << "Wrapping up" << endl;
    }
    //
    Calc_LogLik( nthreads, RiskFail, RiskGroup, totalnum, ntime, R, Rd, Rdd,RdR,RddR, Rls1, Rls2, Rls3, Lls1, Lls2, Lls3, Ll, Lld, Lldd, debugging, ties_method,KeepConstant);
    //
    if ((Ll_abs_best>0)||(Ll_abs_best < Ll[ind0])){
        Ll_abs_best = Ll[ind0];
        beta_abs_best = beta_c;
    }
    //
    if (verbose){
        Rcout << "df501 " << Ll_abs_best << endl;
        Rcout << "df504 ";//prints parameter values
        for (int ij=0;ij<totalnum;ij++){
            Rcout << beta_abs_best[ij] << " ";
        }
        Rcout << " " << endl;
    }
    //
    List para_list = List::create(_["Term_n"]=Term_n,_["tforms"]=tform); //stores the term information
    List control_list = List::create(_["Iteration"]=iteration, _["Maximum Step"]=abs_max, _["Derivative Limiting"]=Lld_worst); //stores the total number of iterations used
    //
    int kept_covs = totalnum - sum(KeepConstant); //does not base the standard deviation off of constant parameters
    NumericVector Lldd_vec(kept_covs * kept_covs);
    #pragma omp parallel for schedule(dynamic) num_threads(nthreads)
    for (int ijk=0;ijk<totalnum*(totalnum+1)/2;ijk++){
        int ij = 0;
        int jk = ijk;
        int pij_ind=-100;
        int pjk_ind=-100;
        while (jk>ij){
            ij++;
            jk-=ij;
        }
        if (KeepConstant[ij]==0){
            pij_ind = ij - sum(head(KeepConstant,ij));
            if (KeepConstant[jk]==0){
                pjk_ind = jk - sum(head(KeepConstant,jk));
                Lldd_vec[pij_ind * kept_covs + pjk_ind]=Lldd[ij*totalnum+jk];
            }
        }
    }
    for (int ijk=0;ijk<kept_covs*(kept_covs+1)/2;ijk++){
        int ij = 0;
        int jk = ijk;
        while (jk>ij){
            ij++;
            jk-=ij;
        }
        Lldd_vec[jk * kept_covs + ij]=Lldd_vec[ij * kept_covs + jk];
    }
    Lldd_vec.attr("dim") = Dimension(kept_covs, kept_covs);
    const Map<MatrixXd> Lldd_mat(as<Map<MatrixXd> >(Lldd_vec));
    //
    MatrixXd Lldd_inv = -1 * Lldd_mat.inverse().matrix(); //uses inverse information matrix to calculate the standard deviation
    VectorXd stdev = VectorXd::Zero(totalnum);
    for (int ij=0;ij<totalnum;ij++){
        if (KeepConstant[ij]==0){
            int pij_ind = ij - sum(head(KeepConstant,ij));
            stdev(ij) = sqrt(Lldd_inv(pij_ind,pij_ind));
        }
    }
    //
    //
    List res_list = List::create(_["LogLik"]=wrap(Ll[0]),_["First_Der"]=wrap(Lld),_["Second_Der"]=Lldd_vec,_["beta_0"]=wrap(beta_0) ,_["Standard_Deviation"]=wrap(stdev) ,_["AIC"]=2*(totalnum-accumulate(KeepConstant.begin(),KeepConstant.end(), 0.0))-2*Ll[0],_["Parameter_Lists"]=para_list,_["Control_List"]=control_list,_["Convgerged"]=convgd);
    // returns a list of results
    return res_list;
}

//' Primary Cox PH regression with multiple starting points and strata effect
//' \code{LogLik_Cox_PH_Guess_Strata} Performs the calls to calculation functions, Structures the Cox PH regression, With verbose option prints out time stamps and intermediate sums of terms and derivatives, with strata effect
//'
//' @param     Term_n    Term numbers
//' @param     tform    subterm types
//' @param     a_ns    matrix starting values
//' @param     x_all    covariate matrix
//' @param     dfc    covariate column numbers
//' @param     fir    first term number
//' @param     der_iden    subterm number for derivative tests
//' @param     modelform    model string
//' @param     lr    learning rate for newton step toward 0 derivative
//' @param     maxiters    list of maximum number of iterations
//' @param     guesses    number of initial guesses
//' @param     halfmax    maximum number of half steps
//' @param     epsilon    minimum acceptable maximum parameter change
//' @param     dbeta_cap    learning rate for newton step toward 0 log-likelihood
//' @param     abs_max    Maximum allowed parameter change
//' @param     dose_abs_max    Maximum allowed threshold parameter change
//' @param     deriv_epsilon    threshold for near-zero derivative
//' @param     df_groups    matrix with time and event information
//' @param     tu    event times
//' @param     double_step controls the step calculation, 0 for independent changes, 1 for solving b=Ax with complete matrices
//' @param     change_all    boolean if every parameter is being updated
//' @param     verbose    verbosity boolean
//' @param     debugging    debugging boolean
//' @param     KeepConstant    vector identifying constant parameters
//' @param     term_tot    total number of terms
//' @param     ties_method    ties method
//' @param     STRATA_vals vector of strata identifier values
//' @param     nthreads number of threads to use
//'
//' @return List of final results: Log-likelihood of optimum, first derivative of log-likelihood, second derivative matrix, parameter list, standard deviation estimate, AIC, model information
// [[Rcpp::export]]
List LogLik_Cox_PH_Guess_Strata( IntegerVector Term_n, StringVector tform, NumericMatrix a_ns,NumericMatrix x_all,IntegerVector dfc,int fir, int der_iden,string modelform, double lr, NumericVector maxiters, int guesses, int halfmax, double epsilon, double dbeta_cap, double abs_max,double dose_abs_max, double deriv_epsilon, NumericMatrix df_groups, NumericVector tu, int double_step ,bool change_all, bool verbose, bool debugging, IntegerVector KeepConstant, int term_tot, string ties_method, NumericVector& STRATA_vals, int nthreads){
    ;
    //
    List temp_list = List::create(_["Status"]="FAILED"); //used as a dummy return value for code checking
    if (verbose){
        Rcout << "START_COX_GUESS_STRATA" << endl;
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
    int reqrdnum = totalnum - sum(KeepConstant);
    //
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
    // int nthreads = Eigen::nbThreads()-1; //stores how many threads are allocated
    //
    // Lld_worst: stores the highest magnitude log-likelihood derivative
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
    NumericVector a_n = a_ns.row(0);
    Map<VectorXd> beta_0(as<Map<VectorXd> >(a_n));
    int maxiter = maxiters[0];
    MatrixXd T0 = MatrixXd::Zero(df0.rows(), totalnum); //preallocates matrix for Term column
    MatrixXd Td0 = MatrixXd::Zero(df0.rows(), reqrdnum); //preallocates matrix for Term derivative columns
    MatrixXd Tdd0 = MatrixXd::Zero(df0.rows(), reqrdnum*(reqrdnum+1)/2); //preallocates matrix for Term second derivative columns
    //
    MatrixXd Te = MatrixXd::Zero(df0.rows(), 1); //preallocates matrix for column terms used for temporary storage
    MatrixXd R = MatrixXd::Zero(df0.rows(), 1); //preallocates matrix for Risks
    ColXd Rd = MatrixXd::Zero(df0.rows(), reqrdnum); //preallocates matrix for Risk derivatives
    ColXd Rdd = MatrixXd::Zero(df0.rows(), reqrdnum*(reqrdnum+1)/2); //preallocates matrix for Risk second derivatives
    //
    MatrixXd Dose = MatrixXd::Constant(df0.rows(),term_tot,0.0); //Matrix of the total dose term values
    MatrixXd nonDose = MatrixXd::Constant(df0.rows(),term_tot,1.0); //Matrix of the total non-dose term values
    MatrixXd nonDose_LIN = MatrixXd::Constant(df0.rows(),term_tot,0.0); //matrix of Linear subterm values
    MatrixXd nonDose_PLIN = MatrixXd::Constant(df0.rows(),term_tot,1.0); //matrix of Loglinear subterm values
    MatrixXd nonDose_LOGLIN = MatrixXd::Constant(df0.rows(),term_tot,1.0); //matrix of Product linear subterm values
    MatrixXd TTerm = MatrixXd::Zero(Dose.rows(),Dose.cols()); //matrix of term values
    double dint = dose_abs_max; //The amount of change used to calculate derivatives in threshold paramters
    double dslp = abs_max;
    ColXd RdR = MatrixXd::Zero(df0.rows(), reqrdnum); //preallocates matrix for Risk to derivative ratios
    ColXd RddR = MatrixXd::Zero(df0.rows(), reqrdnum*(reqrdnum+1)/2); //preallocates matrix for Risk to second derivative ratios
    //
    StringMatrix RiskGroup(ntime,STRATA_vals.size()); //vector of strings detailing the rows
    IntegerMatrix RiskFail(ntime,2*STRATA_vals.size()); //vector giving the event rows
    //
    const Map<MatrixXd> df_m(as<Map<MatrixXd> >(df_groups));
    // Creates matrices used to identify the event risk groups
    Make_Groups_STRATA( ntime, df_m, RiskFail, RiskGroup, tu, nthreads, debugging,STRATA_vals);
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
    //
    MatrixXd Rls1 =MatrixXd::Zero(ntime, STRATA_vals.size()); //precomputes a series of sums used frequently in the log-liklihood calculations
    MatrixXd Rls2 =MatrixXd::Zero(ntime, reqrdnum*STRATA_vals.size()); //Many are repeated due to the same risk groups and derivatives being used at mulitple points
    MatrixXd Rls3 =MatrixXd::Zero(ntime, reqrdnum*(reqrdnum+1)/2*STRATA_vals.size()); //Sum and its derivatives are precomputed
    MatrixXd Lls1 =MatrixXd::Zero(ntime, STRATA_vals.size()); //The log-likelihood calculation has a Right and Left sum used
    MatrixXd Lls2 =MatrixXd::Zero(ntime, reqrdnum*STRATA_vals.size());
    MatrixXd Lls3 =MatrixXd::Zero(ntime, reqrdnum*(reqrdnum+1)/2*STRATA_vals.size());
    //The log-likelihood is calculated in parallel over the risk groups
    vector<double> Ll(reqrdnum,0.0); //Log-likelihood values
    vector<double> Lld(reqrdnum,0.0); //Log-likelihood derivative values
    vector<double> Lldd(pow(reqrdnum,2),0.0);//The second derivative matrix has room for every combination, but only the lower triangle is calculated initially
    //
    vector <double> Ll_comp(2,Ll[0]); //vector to compare values
    double abs_max0 = abs_max;
    double dose_abs_max0 = dose_abs_max;
    //
    vector<double> dbeta(totalnum,0.0);
    //
    // --------------------------
    // always starts from initial guess
    // --------------------------
    vector<double> beta_c(totalnum,0.0);
    vector<double> beta_a(totalnum,0.0);
    vector<double> beta_best(totalnum,0.0);
    vector<double> beta_p(totalnum,0.0);
    VectorXd::Map(&beta_p[0], beta_0.size()) = beta_0;// stores previous parameters
    VectorXd::Map(&beta_c[0], beta_0.size()) = beta_0;// stores current parameters
    VectorXd::Map(&beta_a[0], beta_0.size()) = beta_0;// stores a refrence value for parameters
    VectorXd::Map(&beta_best[0], beta_0.size()) = beta_0;// stores the best parameters
    double halves = 0; //number of half-steps taken
    int ind0 = fir; //used for validations
    int iteration=0; //iteration number
    //
    bool convgd = FALSE;
    int iter_stop =0; //tracks if the iterations should be stopped for convergence
    int iter_check=0; //signal to check for convergence
    //
    NumericMatrix beta_fin(a_ns.rows(), a_ns.cols());
    NumericVector LL_fin(a_ns.rows());
    //
        //
    double Ll_abs_best = 10;
    vector<double> beta_abs_best(totalnum,0.0);
    int guess_abs_best =-1;
    for (int guess=0; guess<=guesses;guess++){
        //
        Dose = MatrixXd::Zero(df0.rows(),term_tot);
        nonDose = MatrixXd::Constant(df0.rows(),term_tot,0.0);
        nonDose_LIN = MatrixXd::Constant(df0.rows(),term_tot,0.0);
        nonDose_PLIN = MatrixXd::Constant(df0.rows(),term_tot,1.0);
        nonDose_LOGLIN = MatrixXd::Constant(df0.rows(),term_tot,1.0);
        T0 = MatrixXd::Zero(df0.rows(), totalnum); //preallocates matrix for Term column
        Td0 = MatrixXd::Zero(df0.rows(), reqrdnum); //preallocates matrix for Term derivative columns
        Tdd0 = MatrixXd::Zero(df0.rows(), reqrdnum*(reqrdnum+1)/2); //preallocates matrix for Term second derivative columns
        TTerm = MatrixXd::Zero(Dose.rows(),Dose.cols());
        RdR = MatrixXd::Zero(df0.rows(), reqrdnum);
        RddR = MatrixXd::Zero(df0.rows(), reqrdnum*(reqrdnum+1)/2);
        fill(Ll.begin(), Ll.end(), 0.0);
        fill(Lld.begin(), Lld.end(), 0.0);
        fill(Lldd.begin(), Lldd.end(), 0.0);
        Rls1 =MatrixXd::Zero(ntime, STRATA_vals.size()); //precomputes a series of sums used frequently in the log-liklihood calculations
        Rls2 =MatrixXd::Zero(ntime, reqrdnum*STRATA_vals.size()); //Many are repeated due to the same risk groups and derivatives being used at mulitple points
        Rls3 =MatrixXd::Zero(ntime, reqrdnum*(reqrdnum+1)/2*STRATA_vals.size()); //Sum and its derivatives are precomputed
        Lls1 =MatrixXd::Zero(ntime, STRATA_vals.size()); //The log-likelihood calculation has a Right and Left sum used
        Lls2 =MatrixXd::Zero(ntime, reqrdnum*STRATA_vals.size());
        Lls3 =MatrixXd::Zero(ntime, reqrdnum*(reqrdnum+1)/2*STRATA_vals.size());
        beta_p = beta_best;//
        beta_a = beta_best;//
        beta_c = beta_best;//
        abs_max = abs_max0;
        dose_abs_max = dose_abs_max0;
        iter_check=0;
        iter_stop = 0;
        halves=0;
        iteration=0;
        halves = 0; //number of half-steps taken
        ind0 = fir; //used for validations
        iteration=0; //iteration number
        //
        convgd = FALSE;
        iter_stop =0; //tracks if the iterations should be stopped for convergence
        iter_check=0; //signal to check for convergence
        //
        maxiter = maxiters[guess];
        a_n = a_ns.row(guess);
        for (int i=0;i<beta_0.size();i++){
            beta_0[i] = a_n[i];
        }
        if (verbose){
            Rcout << "starting subterms " << term_tot << " in guess " << guess << endl;
        }
        //
        // Calculates the subterm and term values
        Make_subterms( totalnum, Term_n, tform, dfc, fir, T0, Td0, Tdd0, Dose, nonDose, TTerm,  nonDose_LIN, nonDose_PLIN, nonDose_LOGLIN ,beta_0, df0,dint,dslp,nthreads, debugging,KeepConstant);
        // ---------------------------------------------------------
        // Prints off a series of calculations to check at what point values are changing
        // ---------------------------------------------------------
        //
        //
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
            for (int ijk=0;ijk<reqrdnum;ijk++){
                Rcout << Td0.col(ijk).sum() << " ";
            }
            Rcout << " " << endl;
            Rcout << "second derivs checked ";
            for (int ijk=0;ijk<reqrdnum;ijk++){
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
        Make_Risks(modelform, tform, Term_n, totalnum, fir, T0, Td0, Tdd0, Te, R, Rd, Rdd, Dose, nonDose, TTerm, nonDose_LIN, nonDose_PLIN, nonDose_LOGLIN, RdR, RddR, nthreads, debugging,KeepConstant);
        //
        // Removes infinite values
        RdR = (RdR.array().isFinite()).select(RdR,0);
        RddR = (RddR.array().isFinite()).select(RddR,0);
        //
        if (verbose){
            Rcout << "risk checked ";
            for (int ijk=0;ijk<1;ijk++){
                Rcout << R.col(0).sum() << " ";
            }
            Rcout << " " << endl;
            Rcout << "risk1 checked ";
            for (int ijk=0;ijk<reqrdnum;ijk++){
                Rcout << Rd.col(ijk).sum() << " ";
            }
            Rcout << " " << endl;
            Rcout << "risk2 checked ";
            for (int ijk=0;ijk<reqrdnum;ijk++){
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
        if (verbose){
            Rcout << "Made Risk Side Lists" << endl;
        }
        // Calculates the side sum terms used
        Calculate_Sides_STRATA( RiskFail, RiskGroup, totalnum, ntime, R, Rd, Rdd, Rls1, Rls2, Rls3, Lls1, Lls2, Lls3,nthreads, debugging, STRATA_vals,KeepConstant);
        //
        if (verbose){
            end_point = system_clock::now();
            ending = time_point_cast<microseconds>(end_point).time_since_epoch().count();
            Rcout <<"df100 "<<(ending-start)<<" "<<0<<" "<<0<<" "<<-1<<",Prep_Sides"<<endl;
            gibtime = system_clock::to_time_t(system_clock::now());
            Rcout << ctime(&gibtime) << endl;
        }
        // Calculates log-likelihood
        Calc_LogLik_STRATA( nthreads, RiskFail, RiskGroup, totalnum, ntime, R, Rd, Rdd,RdR,RddR, Rls1, Rls2, Rls3, Lls1, Lls2, Lls3, Ll, Lld, Lldd, debugging, ties_method, STRATA_vals,KeepConstant);
        //
        if (verbose){
            end_point = system_clock::now();
            ending = time_point_cast<microseconds>(end_point).time_since_epoch().count();
            Rcout <<"df100 "<<(ending-start)<<" "<<0<<" "<<0<<" "<<0<<",Calc"<<endl;//prints the time
            gibtime = system_clock::to_time_t(system_clock::now());
            Rcout << ctime(&gibtime) << endl;
            Rcout << "df101 ";//prints the log-likelihoods
            for (int ij=0;ij<reqrdnum;ij++){
                Rcout << Ll[ij] << " ";
            }
            Rcout << " " << endl;
            Rcout << "df102 ";//prints the first derivatives
            for (int ij=0;ij<reqrdnum;ij++){
                Rcout << Lld[ij] << " ";
            }
            Rcout << " " << endl;
            Rcout << "df103 ";//prints the second derivatives
            for (int ij=0;ij<reqrdnum;ij++){
                Rcout << Lldd[ij*reqrdnum+ij] << " ";
            }
            Lld_worst=0;
            for (int ij=0;ij<reqrdnum;ij++){//locates highest magnitude derivative
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
            for (int ij=0;ij<reqrdnum;ij++){//prints the newton step value for zero derivative
                Rcout << Lld[ij]/Lldd[ij*reqrdnum+ij] << " ";
            }
            Rcout << " " << endl;
            Rcout << "df106 ";
            for (int ij=0;ij<reqrdnum;ij++){//prints the newton step value for zero log-likelihood
                Rcout << Ll[ij]/Lld[ij] << " ";
            }
            Rcout << " " << endl;
            Rcout << "df107 " << double_step << " " << abs_max << " " << dose_abs_max << " " << Ll_comp[0] << " " << Ll_comp[1] << endl;//prints several convergence terms
        }
        //
        for (int i=0;i<beta_0.size();i++){
            beta_c[i] = beta_0[i];
        }
        while ((iteration < maxiter)&&(iter_stop==0)){
            iteration++;
            beta_p = beta_c;//
            beta_a = beta_c;//
            beta_best = beta_c;//
            //
            // Calcualtes the initial change in parameter
            Calc_Change( double_step, nthreads, totalnum, fir, der_iden, dbeta_cap, dose_abs_max, lr, abs_max, Ll, Lld, Lldd, dbeta, change_all, tform, dint,dslp, KeepConstant, debugging);
            Intercept_Bound(nthreads, totalnum, beta_0, dbeta, dfc, df0, KeepConstant, debugging, tform);
            if (verbose){
                Rcout << "Starting Halves"<<endl;//prints the final changes for validation
            }
            //
            //
            if ((Ll_abs_best>0)||(Ll_abs_best < Ll[ind0])){
                Ll_abs_best = Ll[ind0];
                beta_abs_best = beta_c;
                guess_abs_best=guess;
            }
            //
            if (verbose){
                Rcout << "df501 " << Ll_abs_best << endl;
                Rcout << "df504 ";//prints parameter values
                for (int ij=0;ij<totalnum;ij++){
                    Rcout << beta_abs_best[ij] << " ";
                }
                Rcout << " " << endl;
            }
            halves=0;
            while ((Ll[ind0] <= Ll_abs_best)&&(halves<halfmax)){ //repeats until half-steps maxed or an improvement
                //Refreshes the matrices used
                Dose = MatrixXd::Zero(df0.rows(),term_tot);
                nonDose = MatrixXd::Constant(df0.rows(),term_tot,0.0);
                nonDose_LIN = MatrixXd::Constant(df0.rows(),term_tot,0.0);
                nonDose_PLIN = MatrixXd::Constant(df0.rows(),term_tot,1.0);
                nonDose_LOGLIN = MatrixXd::Constant(df0.rows(),term_tot,1.0);
                T0 = MatrixXd::Zero(df0.rows(), totalnum); //preallocates matrix for Term column
                Td0 = MatrixXd::Zero(df0.rows(), reqrdnum); //preallocates matrix for Term derivative columns
                Tdd0 = MatrixXd::Zero(df0.rows(), reqrdnum*(reqrdnum+1)/2); //preallocates matrix for Term second derivative columns
                for (int ij=0;ij<totalnum;ij++){
                    beta_0[ij] = beta_a[ij] + dbeta[ij];
                    beta_c[ij] = beta_0[ij];
                }
                // ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------//
                // The same subterm, risk, sides, and log-likelihood calculations are performed every half-step and iteration
                // ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------//
                Make_subterms( totalnum, Term_n, tform, dfc, fir, T0, Td0, Tdd0, Dose,  nonDose, TTerm,  nonDose_LIN, nonDose_PLIN, nonDose_LOGLIN ,beta_0, df0,dint,dslp,nthreads, debugging,KeepConstant);
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
                    for (int ijk=0;ijk<reqrdnum;ijk++){
                        Rcout << Td0.col(ijk).sum() << " ";
                    }
                    Rcout << " " << endl;
                    Rcout << "second derivs checked ";
                    for (int ijk=0;ijk<reqrdnum;ijk++){
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
                RdR = MatrixXd::Zero(df0.rows(), reqrdnum);
                RddR = MatrixXd::Zero(df0.rows(), reqrdnum*(reqrdnum+1)/2);
                //
                //
                Make_Risks(modelform, tform, Term_n, totalnum, fir, T0, Td0, Tdd0, Te, R, Rd, Rdd, Dose, nonDose, TTerm, nonDose_LIN, nonDose_PLIN, nonDose_LOGLIN, RdR, RddR, nthreads, debugging,KeepConstant);
                //
                if (R.minCoeff()<0){
                    #pragma omp parallel for num_threads(nthreads)
                    for (int ijk=0;ijk<totalnum;ijk++){
                        int tij = Term_n[ijk];
                        if (TTerm.col(tij).minCoeff()<=0){
                            dbeta[ijk] = dbeta[ijk] / 2.0;
                        } else if (isinf(TTerm.col(tij).maxCoeff())){
                            dbeta[ijk] = dbeta[ijk] / 2.0;
                        }
                    }
                    if (verbose){
                        Rcout << "A non-positive risk was detected: " << R.minCoeff() << endl;
                    }
                    halves+=0.2;
                } else {
                    halves++;
                    RdR = (RdR.array().isFinite()).select(RdR,0);
                    RddR = (RddR.array().isFinite()).select(RddR,0);
                    //
                    if (verbose){
                        Rcout << "risk checked ";
                        for (int ijk=0;ijk<1;ijk++){
                            Rcout << R.col(0).sum() << " ";
                        }
                        Rcout << " " << endl;
                        Rcout << "risk1 checked ";
                        for (int ijk=0;ijk<reqrdnum;ijk++){
                            Rcout << Rd.col(ijk).sum() << " ";
                        }
                        Rcout << " " << endl;
                        Rcout << "risk2 checked ";
                        for (int ijk=0;ijk<reqrdnum;ijk++){
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
                    Rls1 =MatrixXd::Zero(ntime, STRATA_vals.size()); //precomputes a series of sums used frequently in the log-liklihood calculations
                    Rls2 =MatrixXd::Zero(ntime, reqrdnum*STRATA_vals.size()); //Many are repeated due to the same risk groups and derivatives being used at mulitple points
                    Rls3 =MatrixXd::Zero(ntime, reqrdnum*(reqrdnum+1)/2*STRATA_vals.size()); //Sum and its derivatives are precomputed
                    Lls1 =MatrixXd::Zero(ntime, STRATA_vals.size()); //The log-likelihood calculation has a Right and Left sum used
                    Lls2 =MatrixXd::Zero(ntime, reqrdnum*STRATA_vals.size());
                    Lls3 =MatrixXd::Zero(ntime, reqrdnum*(reqrdnum+1)/2*STRATA_vals.size());
                    Calculate_Sides_STRATA( RiskFail, RiskGroup, totalnum, ntime, R, Rd, Rdd, Rls1, Rls2, Rls3, Lls1, Lls2, Lls3,nthreads, debugging, STRATA_vals,KeepConstant);
                    //
                    
                    if (verbose){
                        end_point = system_clock::now();
                        ending = time_point_cast<microseconds>(end_point).time_since_epoch().count();
                        Rcout <<"df100 "<<(ending-start)<<" "<<halves<<" "<<iteration<<" "<<ind0<<",Update_Sides"<<endl;
                        gibtime = system_clock::to_time_t(system_clock::now());
                        Rcout << ctime(&gibtime) << endl;
                    }
                    //
                    Calc_LogLik_STRATA( nthreads, RiskFail, RiskGroup, totalnum, ntime, R, Rd, Rdd,RdR,RddR, Rls1, Rls2, Rls3, Lls1, Lls2, Lls3, Ll, Lld, Lldd, debugging, ties_method, STRATA_vals,KeepConstant);
                    if (verbose){
                        end_point = system_clock::now();
                        ending = time_point_cast<microseconds>(end_point).time_since_epoch().count();
                        Rcout <<"df100 "<<(ending-start)<<" "<<halves<<" "<<iteration<<" "<<ind0<<", Half"<<endl;
                        gibtime = system_clock::to_time_t(system_clock::now());
                        Rcout << ctime(&gibtime) << endl;
                        Rcout << "df101 ";
                        for (int ij=0;ij<reqrdnum;ij++){
                            Rcout << Ll[ij] << " ";
                        }
                        Rcout << " " << endl;
                        Rcout << "df102 ";
                        for (int ij=0;ij<reqrdnum;ij++){
                            Rcout << Lld[ij] << " ";
                        }
                        Rcout << " " << endl;
                        Rcout << "df103 ";
                        for (int ij=0;ij<reqrdnum;ij++){
                            Rcout << Lldd[ij*reqrdnum+ij] << " ";
                        }
                        Rcout << " " << endl;
                        Rcout << "df104 ";
                        for (int ij=0;ij<totalnum;ij++){
                            Rcout << beta_c[ij] << " ";
                        }
                        Rcout << " " << endl;
                        Rcout << "df105 ";
                        for (int ij=0;ij<reqrdnum;ij++){//prints the newton step value for zero derivative
                            Rcout << Lld[ij]/Lldd[ij*reqrdnum+ij] << " ";
                        }
                        Rcout << " " << endl;
                        Rcout << "df106 ";
                        for (int ij=0;ij<reqrdnum;ij++){//prints the newton step value for zero log-likelihood
                            Rcout << Ll[ij]/Lld[ij] << " ";
                        }
                        Rcout << " " << endl;
                        Rcout << "df107 " << double_step << " " << abs_max << " " << dose_abs_max << " " << Ll_comp[0] << " " << Ll_comp[1] << endl;//prints several convergence terms
                    }
                    //
                    if (change_all){ //If every covariate is to be changed
                        if (Ll[ind0] <= Ll_abs_best){//takes a half-step if needed
                            #pragma omp parallel for num_threads(nthreads)
                            for (int ijk=0;ijk<totalnum;ijk++){
                                dbeta[ijk] = dbeta[ijk] * 0.5; //
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
                    }
                    #pragma omp parallel for num_threads(nthreads)
                    for (int ijk=0;ijk<totalnum;ijk++){//totalnum*(totalnum+1)/2
                        beta_0[ijk] = beta_c[ijk];
                    }
                }
            }
            if (beta_best!=beta_c){//if the risk matrices aren't the optimal values, then they must be recalculated
                if (verbose){
                    Rcout << "Changing back to best"<<endl;
                }
                // If it goes through every half step without improvement, then the maximum change needs to be decreased
                abs_max = abs_max*pow(0.5,halfmax); // reduces the step sizes
                dose_abs_max = dose_abs_max*pow(0.5,halfmax);
                iter_check = 1;
                //
                beta_p = beta_best;//
                beta_a = beta_best;//
                beta_c = beta_best;//
                Dose = MatrixXd::Zero(df0.rows(),term_tot);
                nonDose = MatrixXd::Constant(df0.rows(),term_tot,0.0);
                nonDose_LIN = MatrixXd::Constant(df0.rows(),term_tot,0.0);
                nonDose_PLIN = MatrixXd::Constant(df0.rows(),term_tot,1.0);
                nonDose_LOGLIN = MatrixXd::Constant(df0.rows(),term_tot,1.0);
                for (int ij=0;ij<totalnum;ij++){
                    beta_0[ij] = beta_best[ij];
                }
                Make_subterms( totalnum, Term_n, tform, dfc, fir, T0, Td0, Tdd0, Dose, nonDose, TTerm,  nonDose_LIN, nonDose_PLIN, nonDose_LOGLIN,beta_0, df0,dint,dslp,nthreads, debugging, KeepConstant);
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
                    for (int ijk=0;ijk<reqrdnum;ijk++){
                        Rcout << Td0.col(ijk).sum() << " ";
                    }
                    Rcout << " " << endl;
                    Rcout << "second derivs checked ";
                    for (int ijk=0;ijk<reqrdnum;ijk++){
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
                RdR = MatrixXd::Zero(df0.rows(), reqrdnum);
                RddR = MatrixXd::Zero(df0.rows(), reqrdnum*(reqrdnum+1)/2);
                //
                //
                Make_Risks(modelform, tform, Term_n, totalnum, fir, T0, Td0, Tdd0, Te, R, Rd, Rdd, Dose, nonDose, TTerm, nonDose_LIN, nonDose_PLIN, nonDose_LOGLIN, RdR, RddR, nthreads, debugging,KeepConstant);
                R = (R.array().isFinite()).select(R,0);
                Rd = (Rd.array().isFinite()).select(Rd,0);
                Rdd = (Rdd.array().isFinite()).select(Rdd,0);
                RdR = (RdR.array().isFinite()).select(RdR,0);
                RddR = (RddR.array().isFinite()).select(RddR,0);
                //
                temp_list = List::create(_["betas"]=wrap(beta_0),_["Status"]="FAILED"); //used as a dummy return value for code checking
                if (verbose){
                    Rcout << "risk checked ";
                    for (int ijk=0;ijk<1;ijk++){
                        Rcout << R.col(0).sum() << " ";
                    }
                    Rcout << " " << endl;
                    Rcout << "risk1 checked ";
                    for (int ijk=0;ijk<reqrdnum;ijk++){
                        Rcout << Rd.col(ijk).sum() << " ";
                    }
                    Rcout << " " << endl;
                    Rcout << "risk2 checked ";
                    for (int ijk=0;ijk<reqrdnum;ijk++){
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
                Rls1 =MatrixXd::Zero(ntime, STRATA_vals.size()); //precomputes a series of sums used frequently in the log-liklihood calculations
                Rls2 =MatrixXd::Zero(ntime, reqrdnum*STRATA_vals.size()); //Many are repeated due to the same risk groups and derivatives being used at mulitple points
                Rls3 =MatrixXd::Zero(ntime, reqrdnum*(reqrdnum+1)/2*STRATA_vals.size()); //Sum and its derivatives are precomputed
                Lls1 =MatrixXd::Zero(ntime, STRATA_vals.size()); //The log-likelihood calculation has a Right and Left sum used
                Lls2 =MatrixXd::Zero(ntime, reqrdnum*STRATA_vals.size());
                Lls3 =MatrixXd::Zero(ntime, reqrdnum*(reqrdnum+1)/2*STRATA_vals.size());
                Calculate_Sides_STRATA( RiskFail, RiskGroup, totalnum, ntime, R, Rd, Rdd, Rls1, Rls2, Rls3, Lls1, Lls2, Lls3,nthreads, debugging, STRATA_vals,KeepConstant);
                //
                if (verbose){
                    end_point = system_clock::now();
                    ending = time_point_cast<microseconds>(end_point).time_since_epoch().count();
                    Rcout <<"df100 "<<(ending-start)<<" "<<halves<<" "<<iteration<<" "<<ind0<<",Update_Sides"<<endl;
                    gibtime = system_clock::to_time_t(system_clock::now());
                    Rcout << ctime(&gibtime) << endl;
                }
                //
                Calc_LogLik_STRATA( nthreads, RiskFail, RiskGroup, totalnum, ntime, R, Rd, Rdd,RdR,RddR, Rls1, Rls2, Rls3, Lls1, Lls2, Lls3, Ll, Lld, Lldd, debugging, ties_method, STRATA_vals,KeepConstant);
            }
            Lld_worst=0;
            for (int ij=0;ij<reqrdnum;ij++){
                if (abs(Lld[ij]) > Lld_worst){
                    Lld_worst = abs(Lld[ij]);
                }
            }
            if (iteration > reqrdnum){//Doesn't check the first several iterations for convergence
                if ((iteration % (reqrdnum))||(iter_check==1)){//Checks every set number of iterations
                    iter_check=0;
                    if (Lld_worst < deriv_epsilon){//ends if the derivatives are low enough
                        iter_stop = 1;
                        convgd = TRUE;
                    }
                    Ll_comp[1]=Ll[0];
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
                for (int ij=0;ij<reqrdnum;ij++){
                    Rcout << Ll[ij] << " ";
                }
                Rcout << " " << endl;
                Rcout << "df102 ";
                for (int ij=0;ij<reqrdnum;ij++){
                    Rcout << Lld[ij] << " ";
                }
                Rcout << " " << endl;
                Rcout << "df103 ";
                for (int ij=0;ij<reqrdnum;ij++){
                    Rcout << Lldd[ij*reqrdnum+ij] << " ";
                }
                Rcout << " " << endl;
                Rcout << "df104 ";
                for (int ij=0;ij<totalnum;ij++){
                    Rcout << beta_c[ij] << " ";
                }
                Rcout << " " << endl;
                Rcout << "df105 ";
                for (int ij=0;ij<reqrdnum;ij++){//prints the newton step value for zero derivative
                    Rcout << Lld[ij]/Lldd[ij*reqrdnum+ij] << " ";
                }
                Rcout << " " << endl;
                Rcout << "df106 ";
                for (int ij=0;ij<reqrdnum;ij++){//prints the newton step value for zero log-likelihood
                    Rcout << Ll[ij]/Lld[ij] << " ";
                }
                Rcout << " " << endl;
                Rcout << "df107 " << double_step << " " << abs_max << " " << dose_abs_max << " " << Ll_comp[0] << " " << Ll_comp[1] << endl;//prints several convergence terms
                Rcout << "Finshed iteration" << endl;
            }
        }
        // -----------------------------------------------
        // Performing Full Calculation to get full second derivative matrix
        // -----------------------------------------------
        fill(Ll.begin(), Ll.end(), 0.0);
        fill(Lld.begin(), Lld.end(), 0.0);
        fill(Lldd.begin(), Lldd.end(), 0.0);
        Calculate_Sides_STRATA( RiskFail, RiskGroup, totalnum, ntime, R, Rd, Rdd, Rls1, Rls2, Rls3, Lls1, Lls2, Lls3,nthreads, debugging, STRATA_vals,KeepConstant);
        //
        if (verbose){
            end_point = system_clock::now();
            ending = time_point_cast<microseconds>(end_point).time_since_epoch().count();
            Rcout <<"df100 "<<(ending-start)<<" "<<halves<<" "<<iteration<<" "<<ind0<<",Update_Sides"<<endl;
            gibtime = system_clock::to_time_t(system_clock::now());
            Rcout << ctime(&gibtime) << endl;
            Rcout << "Wrapping up" << endl;
        }
        //
        Calc_LogLik_STRATA( nthreads, RiskFail, RiskGroup, totalnum, ntime, R, Rd, Rdd,RdR,RddR, Rls1, Rls2, Rls3, Lls1, Lls2, Lls3, Ll, Lld, Lldd, debugging, ties_method, STRATA_vals,KeepConstant);
        //
        a_n = beta_0;
        //
        beta_fin(guess, _) = a_n;
        LL_fin[guess] = Ll[0];
        if ((Ll_abs_best>0)||(Ll_abs_best < Ll[ind0])){
            Ll_abs_best = Ll[ind0];
            beta_abs_best = beta_c;
            guess_abs_best=guess;
        }
        //
        if (verbose){
            Rcout << "df501 " << Ll_abs_best << endl;
            Rcout << "df504 ";//prints parameter values
            for (int ij=0;ij<totalnum;ij++){
                Rcout << beta_abs_best[ij] << " ";
            }
            Rcout << " " << endl;
        }
    }
    //
    Dose = MatrixXd::Zero(df0.rows(),term_tot);
    nonDose = MatrixXd::Constant(df0.rows(),term_tot,0.0);
    nonDose_LIN = MatrixXd::Constant(df0.rows(),term_tot,0.0);
    nonDose_PLIN = MatrixXd::Constant(df0.rows(),term_tot,1.0);
    nonDose_LOGLIN = MatrixXd::Constant(df0.rows(),term_tot,1.0);
    T0 = MatrixXd::Zero(df0.rows(), totalnum); //preallocates matrix for Term column
    Td0 = MatrixXd::Zero(df0.rows(), reqrdnum); //preallocates matrix for Term derivative columns
    Tdd0 = MatrixXd::Zero(df0.rows(), reqrdnum*(reqrdnum+1)/2); //preallocates matrix for Term second derivative columns
    TTerm = MatrixXd::Zero(Dose.rows(),Dose.cols());
    RdR = MatrixXd::Zero(df0.rows(), reqrdnum);
    RddR = MatrixXd::Zero(df0.rows(), reqrdnum*(reqrdnum+1)/2);
    fill(Ll.begin(), Ll.end(), 0.0);
    fill(Lld.begin(), Lld.end(), 0.0);
    fill(Lldd.begin(), Lldd.end(), 0.0);
    Rls1 =MatrixXd::Zero(ntime, STRATA_vals.size()); //precomputes a series of sums used frequently in the log-liklihood calculations
    Rls2 =MatrixXd::Zero(ntime, reqrdnum*STRATA_vals.size()); //Many are repeated due to the same risk groups and derivatives being used at mulitple points
    Rls3 =MatrixXd::Zero(ntime, reqrdnum*(reqrdnum+1)/2*STRATA_vals.size()); //Sum and its derivatives are precomputed
    Lls1 =MatrixXd::Zero(ntime, STRATA_vals.size()); //The log-likelihood calculation has a Right and Left sum used
    Lls2 =MatrixXd::Zero(ntime, reqrdnum*STRATA_vals.size());
    Lls3 =MatrixXd::Zero(ntime, reqrdnum*(reqrdnum+1)/2*STRATA_vals.size());
    beta_p = beta_best;//
    beta_a = beta_best;//
    beta_c = beta_best;//
    abs_max = abs_max0;
    dose_abs_max = dose_abs_max0;
    iter_check=0;
    iter_stop = 0;
    halves=0;
    iteration=0;
    halves = 0; //number of half-steps taken
    ind0 = fir; //used for validations
    iteration=0; //iteration number
    //
    convgd = FALSE;
    iter_stop =0; //tracks if the iterations should be stopped for convergence
    iter_check=0; //signal to check for convergence
    //
    int guess_max=guess_abs_best;
    if (verbose){
        Rcout << "Guess number, parameter values, Log-Likelihood" << endl;
        NumericVector beta_temp;
        for (int i=0; i<=guesses; i++){
            beta_temp = wrap(beta_fin.row(i));
            if (i==guess_max){
                Rcout << i << ", " << beta_temp << ", " << LL_fin[i] << "<-- Best Guess" << endl;
            } else {
                Rcout << i << ", " << beta_temp << ", " << LL_fin[i] << endl;
            }
        }
    }
    //
    maxiter = maxiters[guesses+1];
    a_n = beta_abs_best;
    for (int i=0;i<beta_0.size();i++){
        beta_0[i] = a_n[i];
    }
    for (int i=0;i<beta_0.size();i++){
        beta_c[i] = beta_0[i];
    }
    //
    if (verbose){
        Rcout << "starting subterms " << term_tot << " in best guess " << guess_max << endl;
    }
    //
    // Calculates the subterm and term values
    Make_subterms( totalnum, Term_n, tform, dfc, fir, T0, Td0, Tdd0, Dose, nonDose, TTerm,  nonDose_LIN, nonDose_PLIN, nonDose_LOGLIN ,beta_0, df0,dint,dslp,nthreads, debugging,KeepConstant);
    // ---------------------------------------------------------
    // Prints off a series of calculations to check at what point values are changing
    // ---------------------------------------------------------
    //
    //
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
        for (int ijk=0;ijk<reqrdnum;ijk++){
            Rcout << Td0.col(ijk).sum() << " ";
        }
        Rcout << " " << endl;
        Rcout << "second derivs checked ";
        for (int ijk=0;ijk<reqrdnum;ijk++){
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
    Make_Risks(modelform, tform, Term_n, totalnum, fir, T0, Td0, Tdd0, Te, R, Rd, Rdd, Dose, nonDose, TTerm, nonDose_LIN, nonDose_PLIN, nonDose_LOGLIN, RdR, RddR, nthreads, debugging,KeepConstant);
    //
    // Removes infinite values
    RdR = (RdR.array().isFinite()).select(RdR,0);
    RddR = (RddR.array().isFinite()).select(RddR,0);
    //
    //
    if (verbose){
        Rcout << "risk checked ";
        for (int ijk=0;ijk<1;ijk++){
            Rcout << R.col(0).sum() << " ";
        }
        Rcout << " " << endl;
        Rcout << "risk1 checked ";
        for (int ijk=0;ijk<reqrdnum;ijk++){
            Rcout << Rd.col(ijk).sum() << " ";
        }
        Rcout << " " << endl;
        Rcout << "risk2 checked ";
        for (int ijk=0;ijk<reqrdnum;ijk++){
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
    if (verbose){
        Rcout << "Made Risk Side Lists" << endl;
    }
    // Calculates the side sum terms used
    Calculate_Sides_STRATA( RiskFail, RiskGroup, totalnum, ntime, R, Rd, Rdd, Rls1, Rls2, Rls3, Lls1, Lls2, Lls3,nthreads, debugging, STRATA_vals,KeepConstant);
    //
    if (verbose){
        end_point = system_clock::now();
        ending = time_point_cast<microseconds>(end_point).time_since_epoch().count();
        Rcout <<"df100 "<<(ending-start)<<" "<<0<<" "<<0<<" "<<-1<<",Prep_Sides"<<endl;
        gibtime = system_clock::to_time_t(system_clock::now());
        Rcout << ctime(&gibtime) << endl;
    }
    //
    // Calculates log-likelihood
    Calc_LogLik_STRATA( nthreads, RiskFail, RiskGroup, totalnum, ntime, R, Rd, Rdd,RdR,RddR, Rls1, Rls2, Rls3, Lls1, Lls2, Lls3, Ll, Lld, Lldd, debugging, ties_method, STRATA_vals,KeepConstant);
    //
    //
    if (verbose){
        end_point = system_clock::now();
        ending = time_point_cast<microseconds>(end_point).time_since_epoch().count();
        Rcout <<"df100 "<<(ending-start)<<" "<<0<<" "<<0<<" "<<0<<",Calc"<<endl;//prints the time
        gibtime = system_clock::to_time_t(system_clock::now());
        Rcout << ctime(&gibtime) << endl;
        Rcout << "df101 ";//prints the log-likelihoods
        for (int ij=0;ij<reqrdnum;ij++){
            Rcout << Ll[ij] << " ";
        }
        Rcout << " " << endl;
        Rcout << "df102 ";//prints the first derivatives
        for (int ij=0;ij<reqrdnum;ij++){
            Rcout << Lld[ij] << " ";
        }
        Rcout << " " << endl;
        Rcout << "df103 ";//prints the second derivatives
        for (int ij=0;ij<reqrdnum;ij++){
            Rcout << Lldd[ij*reqrdnum+ij] << " ";
        }
        Lld_worst=0;
        for (int ij=0;ij<reqrdnum;ij++){//locates highest magnitude derivative
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
        for (int ij=0;ij<reqrdnum;ij++){//prints the newton step value for zero derivative
            Rcout << Lld[ij]/Lldd[ij*reqrdnum+ij] << " ";
        }
        Rcout << " " << endl;
        Rcout << "df106 ";
        for (int ij=0;ij<reqrdnum;ij++){//prints the newton step value for zero log-likelihood
            Rcout << Ll[ij]/Lld[ij] << " ";
        }
        Rcout << " " << endl;
        Rcout << "df107 " << double_step << " " << abs_max << " " << dose_abs_max << " " << Ll_comp[0] << " " << Ll_comp[1] << endl;//prints several convergence terms
    }
    //
    while ((iteration < maxiter)&&(iter_stop==0)){
        iteration++;
        beta_p = beta_c;//
        beta_a = beta_c;//
        beta_best = beta_c;//
        //
        // Calculates the initial change in parameter
        Calc_Change( double_step, nthreads, totalnum, fir, der_iden, dbeta_cap, dose_abs_max, lr, abs_max, Ll, Lld, Lldd, dbeta, change_all, tform, dint,dslp, KeepConstant, debugging);
        Intercept_Bound(nthreads, totalnum, beta_0, dbeta, dfc, df0, KeepConstant, debugging, tform);
        if (verbose){
            Rcout << "Starting Halves"<<endl;//prints the final changes for validation
        }
        //
        //
        if ((Ll_abs_best>0)||(Ll_abs_best < Ll[ind0])){
            Ll_abs_best = Ll[ind0];
            beta_abs_best = beta_c;
        }
        //
        if (verbose){
            Rcout << "df501 " << Ll_abs_best << endl;
            Rcout << "df504 ";//prints parameter values
            for (int ij=0;ij<totalnum;ij++){
                Rcout << beta_abs_best[ij] << " ";
            }
            Rcout << " " << endl;
        }
        halves=0;
        while ((Ll[ind0] <= Ll_abs_best)&&(halves<halfmax)){ //repeats until half-steps maxed or an improvement
            //Refreshes the matrices used
            Dose = MatrixXd::Zero(df0.rows(),term_tot);
            nonDose = MatrixXd::Constant(df0.rows(),term_tot,0.0);
            nonDose_LIN = MatrixXd::Constant(df0.rows(),term_tot,0.0);
            nonDose_PLIN = MatrixXd::Constant(df0.rows(),term_tot,1.0);
            nonDose_LOGLIN = MatrixXd::Constant(df0.rows(),term_tot,1.0);
            T0 = MatrixXd::Zero(df0.rows(), totalnum); //preallocates matrix for Term column
            Td0 = MatrixXd::Zero(df0.rows(), reqrdnum); //preallocates matrix for Term derivative columns
            Tdd0 = MatrixXd::Zero(df0.rows(), reqrdnum*(reqrdnum+1)/2); //preallocates matrix for Term second derivative columns
            for (int ij=0;ij<totalnum;ij++){
                beta_0[ij] = beta_a[ij] + dbeta[ij];
                beta_c[ij] = beta_0[ij];
            }
            // ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------//
            // The same subterm, risk, sides, and log-likelihood calculations are performed every half-step and iteration
            // ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------//
            Make_subterms( totalnum, Term_n, tform, dfc, fir, T0, Td0, Tdd0, Dose,  nonDose, TTerm,  nonDose_LIN, nonDose_PLIN, nonDose_LOGLIN ,beta_0, df0,dint,dslp,nthreads, debugging,KeepConstant);
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
                for (int ijk=0;ijk<reqrdnum;ijk++){
                    Rcout << Td0.col(ijk).sum() << " ";
                }
                Rcout << " " << endl;
                Rcout << "second derivs checked ";
                for (int ijk=0;ijk<reqrdnum;ijk++){
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
            RdR = MatrixXd::Zero(df0.rows(), reqrdnum);
            RddR = MatrixXd::Zero(df0.rows(), reqrdnum*(reqrdnum+1)/2);
            //
            //
            Make_Risks(modelform, tform, Term_n, totalnum, fir, T0, Td0, Tdd0, Te, R, Rd, Rdd, Dose, nonDose, TTerm, nonDose_LIN, nonDose_PLIN, nonDose_LOGLIN, RdR, RddR, nthreads, debugging,KeepConstant);
            //
            if (R.minCoeff()<0){
                #pragma omp parallel for num_threads(nthreads)
                for (int ijk=0;ijk<totalnum;ijk++){
                    int tij = Term_n[ijk];
                    if (TTerm.col(tij).minCoeff()<=0){
                        dbeta[ijk] = dbeta[ijk] / 2.0;
                    } else if (isinf(TTerm.col(tij).maxCoeff())){
                        dbeta[ijk] = dbeta[ijk] / 2.0;
                    }
                }
                if (verbose){
                    Rcout << "A non-positive risk was detected: " << R.minCoeff() << endl;
                }
                halves+=0.2;
            } else {
                halves++;
                RdR = (RdR.array().isFinite()).select(RdR,0);
                RddR = (RddR.array().isFinite()).select(RddR,0);
                //
                if (verbose){
                    Rcout << "risk checked ";
                    for (int ijk=0;ijk<1;ijk++){
                        Rcout << R.col(0).sum() << " ";
                    }
                    Rcout << " " << endl;
                    Rcout << "risk1 checked ";
                    for (int ijk=0;ijk<reqrdnum;ijk++){
                        Rcout << Rd.col(ijk).sum() << " ";
                    }
                    Rcout << " " << endl;
                    Rcout << "risk2 checked ";
                    for (int ijk=0;ijk<reqrdnum;ijk++){
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
                Rls1 =MatrixXd::Zero(ntime, STRATA_vals.size()); //precomputes a series of sums used frequently in the log-liklihood calculations
                Rls2 =MatrixXd::Zero(ntime, reqrdnum*STRATA_vals.size()); //Many are repeated due to the same risk groups and derivatives being used at mulitple points
                Rls3 =MatrixXd::Zero(ntime, reqrdnum*(reqrdnum+1)/2*STRATA_vals.size()); //Sum and its derivatives are precomputed
                Lls1 =MatrixXd::Zero(ntime, STRATA_vals.size()); //The log-likelihood calculation has a Right and Left sum used
                Lls2 =MatrixXd::Zero(ntime, reqrdnum*STRATA_vals.size());
                Lls3 =MatrixXd::Zero(ntime, reqrdnum*(reqrdnum+1)/2*STRATA_vals.size());
                Calculate_Sides_STRATA( RiskFail, RiskGroup, totalnum, ntime, R, Rd, Rdd, Rls1, Rls2, Rls3, Lls1, Lls2, Lls3,nthreads, debugging, STRATA_vals,KeepConstant);
                //
                
                if (verbose){
                    end_point = system_clock::now();
                    ending = time_point_cast<microseconds>(end_point).time_since_epoch().count();
                    Rcout <<"df100 "<<(ending-start)<<" "<<halves<<" "<<iteration<<" "<<ind0<<",Update_Sides"<<endl;
                    gibtime = system_clock::to_time_t(system_clock::now());
                    Rcout << ctime(&gibtime) << endl;
                }
                //
                Calc_LogLik_STRATA( nthreads, RiskFail, RiskGroup, totalnum, ntime, R, Rd, Rdd,RdR,RddR, Rls1, Rls2, Rls3, Lls1, Lls2, Lls3, Ll, Lld, Lldd, debugging, ties_method, STRATA_vals,KeepConstant);
                if (verbose){
                    end_point = system_clock::now();
                    ending = time_point_cast<microseconds>(end_point).time_since_epoch().count();
                    Rcout <<"df100 "<<(ending-start)<<" "<<halves<<" "<<iteration<<" "<<ind0<<", Half"<<endl;
                    gibtime = system_clock::to_time_t(system_clock::now());
                    Rcout << ctime(&gibtime) << endl;
                    Rcout << "df101 ";
                    for (int ij=0;ij<reqrdnum;ij++){
                        Rcout << Ll[ij] << " ";
                    }
                    Rcout << " " << endl;
                    Rcout << "df102 ";
                    for (int ij=0;ij<reqrdnum;ij++){
                        Rcout << Lld[ij] << " ";
                    }
                    Rcout << " " << endl;
                    Rcout << "df103 ";
                    for (int ij=0;ij<reqrdnum;ij++){
                        Rcout << Lldd[ij*reqrdnum+ij] << " ";
                    }
                    Rcout << " " << endl;
                    Rcout << "df104 ";
                    for (int ij=0;ij<totalnum;ij++){
                        Rcout << beta_c[ij] << " ";
                    }
                    Rcout << " " << endl;
                    Rcout << "df105 ";
                    for (int ij=0;ij<reqrdnum;ij++){//prints the newton step value for zero derivative
                        Rcout << Lld[ij]/Lldd[ij*reqrdnum+ij] << " ";
                    }
                    Rcout << " " << endl;
                    Rcout << "df106 ";
                    for (int ij=0;ij<reqrdnum;ij++){//prints the newton step value for zero log-likelihood
                        Rcout << Ll[ij]/Lld[ij] << " ";
                    }
                    Rcout << " " << endl;
                    Rcout << "df107 " << double_step << " " << abs_max << " " << dose_abs_max << " " << Ll_comp[0] << " " << Ll_comp[1] << endl;//prints several convergence terms
                }
                //
                if (change_all){ //If every covariate is to be changed
                    if (Ll[ind0] <= Ll_abs_best){//takes a half-step if needed
                        #pragma omp parallel for num_threads(nthreads)
                        for (int ijk=0;ijk<totalnum;ijk++){
                            dbeta[ijk] = dbeta[ijk] * 0.5; //
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
                }
                #pragma omp parallel for num_threads(nthreads)
                for (int ijk=0;ijk<totalnum;ijk++){//totalnum*(totalnum+1)/2
                    beta_0[ijk] = beta_c[ijk];
                }
            }
        }
        if (beta_best!=beta_c){//if the risk matrices aren't the optimal values, then they must be recalculated
            if (verbose){
                Rcout << "Changing back to best"<<endl;
            }
            // If it goes through every half step without improvement, then the maximum change needs to be decreased
            abs_max = abs_max*pow(0.5,halfmax); // reduces the step sizes
            dose_abs_max = dose_abs_max*pow(0.5,halfmax);
            iter_check = 1;
            //
            beta_p = beta_best;//
            beta_a = beta_best;//
            beta_c = beta_best;//
            Dose = MatrixXd::Zero(df0.rows(),term_tot);
            nonDose = MatrixXd::Constant(df0.rows(),term_tot,0.0);
            nonDose_LIN = MatrixXd::Constant(df0.rows(),term_tot,0.0);
            nonDose_PLIN = MatrixXd::Constant(df0.rows(),term_tot,1.0);
            nonDose_LOGLIN = MatrixXd::Constant(df0.rows(),term_tot,1.0);
            for (int ij=0;ij<totalnum;ij++){
                beta_0[ij] = beta_best[ij];
            }
            Make_subterms( totalnum, Term_n, tform, dfc, fir, T0, Td0, Tdd0, Dose, nonDose, TTerm,  nonDose_LIN, nonDose_PLIN, nonDose_LOGLIN,beta_0, df0,dint,dslp,nthreads, debugging, KeepConstant);
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
                for (int ijk=0;ijk<reqrdnum;ijk++){
                    Rcout << Td0.col(ijk).sum() << " ";
                }
                Rcout << " " << endl;
                Rcout << "second derivs checked ";
                for (int ijk=0;ijk<reqrdnum;ijk++){
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
            RdR = MatrixXd::Zero(df0.rows(), reqrdnum);
            RddR = MatrixXd::Zero(df0.rows(), reqrdnum*(reqrdnum+1)/2);
            //
            //
            Make_Risks(modelform, tform, Term_n, totalnum, fir, T0, Td0, Tdd0, Te, R, Rd, Rdd, Dose, nonDose, TTerm, nonDose_LIN, nonDose_PLIN, nonDose_LOGLIN, RdR, RddR, nthreads, debugging,KeepConstant);
            R = (R.array().isFinite()).select(R,0);
            Rd = (Rd.array().isFinite()).select(Rd,0);
            Rdd = (Rdd.array().isFinite()).select(Rdd,0);
            RdR = (RdR.array().isFinite()).select(RdR,0);
            RddR = (RddR.array().isFinite()).select(RddR,0);
            //
            temp_list = List::create(_["betas"]=wrap(beta_0),_["Status"]="FAILED"); //used as a dummy return value for code checking
            if (verbose){
                Rcout << "risk checked ";
                for (int ijk=0;ijk<1;ijk++){
                    Rcout << R.col(0).sum() << " ";
                }
                Rcout << " " << endl;
                Rcout << "risk1 checked ";
                for (int ijk=0;ijk<reqrdnum;ijk++){
                    Rcout << Rd.col(ijk).sum() << " ";
                }
                Rcout << " " << endl;
                Rcout << "risk2 checked ";
                for (int ijk=0;ijk<reqrdnum;ijk++){
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
            Rls1 =MatrixXd::Zero(ntime, STRATA_vals.size()); //precomputes a series of sums used frequently in the log-liklihood calculations
            Rls2 =MatrixXd::Zero(ntime, reqrdnum*STRATA_vals.size()); //Many are repeated due to the same risk groups and derivatives being used at mulitple points
            Rls3 =MatrixXd::Zero(ntime, reqrdnum*(reqrdnum+1)/2*STRATA_vals.size()); //Sum and its derivatives are precomputed
            Lls1 =MatrixXd::Zero(ntime, STRATA_vals.size()); //The log-likelihood calculation has a Right and Left sum used
            Lls2 =MatrixXd::Zero(ntime, reqrdnum*STRATA_vals.size());
            Lls3 =MatrixXd::Zero(ntime, reqrdnum*(reqrdnum+1)/2*STRATA_vals.size());
            Calculate_Sides_STRATA( RiskFail, RiskGroup, totalnum, ntime, R, Rd, Rdd, Rls1, Rls2, Rls3, Lls1, Lls2, Lls3,nthreads, debugging, STRATA_vals,KeepConstant);
            //
            if (verbose){
                end_point = system_clock::now();
                ending = time_point_cast<microseconds>(end_point).time_since_epoch().count();
                Rcout <<"df100 "<<(ending-start)<<" "<<halves<<" "<<iteration<<" "<<ind0<<",Update_Sides"<<endl;
                gibtime = system_clock::to_time_t(system_clock::now());
                Rcout << ctime(&gibtime) << endl;
            }
            //
            Calc_LogLik_STRATA( nthreads, RiskFail, RiskGroup, totalnum, ntime, R, Rd, Rdd,RdR,RddR, Rls1, Rls2, Rls3, Lls1, Lls2, Lls3, Ll, Lld, Lldd, debugging, ties_method, STRATA_vals,KeepConstant);
        }
        Lld_worst=0;
        for (int ij=0;ij<reqrdnum;ij++){
            if (abs(Lld[ij]) > Lld_worst){
                Lld_worst = abs(Lld[ij]);
            }
        }
        if (iteration > reqrdnum){//Doesn't check the first several iterations for convergence
            if ((iteration % (reqrdnum))||(iter_check==1)){//Checks every set number of iterations
                iter_check=0;
                if (Lld_worst < deriv_epsilon){//ends if the derivatives are low enough
                    iter_stop = 1;
                    convgd = TRUE;
                }
                Ll_comp[1]=Ll[0];
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
            for (int ij=0;ij<reqrdnum;ij++){
                Rcout << Ll[ij] << " ";
            }
            Rcout << " " << endl;
            Rcout << "df102 ";
            for (int ij=0;ij<reqrdnum;ij++){
                Rcout << Lld[ij] << " ";
            }
            Rcout << " " << endl;
            Rcout << "df103 ";
            for (int ij=0;ij<reqrdnum;ij++){
                Rcout << Lldd[ij*reqrdnum+ij] << " ";
            }
            Rcout << " " << endl;
            Rcout << "df104 ";
            for (int ij=0;ij<totalnum;ij++){
                Rcout << beta_c[ij] << " ";
            }
            Rcout << " " << endl;
            Rcout << "df105 ";
            for (int ij=0;ij<reqrdnum;ij++){//prints the newton step value for zero derivative
                Rcout << Lld[ij]/Lldd[ij*reqrdnum+ij] << " ";
            }
            Rcout << " " << endl;
            Rcout << "df106 ";
            for (int ij=0;ij<reqrdnum;ij++){//prints the newton step value for zero log-likelihood
                Rcout << Ll[ij]/Lld[ij] << " ";
            }
            Rcout << " " << endl;
            Rcout << "df107 " << double_step << " " << abs_max << " " << dose_abs_max << " " << Ll_comp[0] << " " << Ll_comp[1] << endl;//prints several convergence terms
            Rcout << "Finshed iteration" << endl;
        }
    }
    // -----------------------------------------------
    // Performing Full Calculation to get full second derivative matrix
    // -----------------------------------------------
    fill(Ll.begin(), Ll.end(), 0.0);
    fill(Lld.begin(), Lld.end(), 0.0);
    fill(Lldd.begin(), Lldd.end(), 0.0);
    Calculate_Sides_STRATA( RiskFail, RiskGroup, totalnum, ntime, R, Rd, Rdd, Rls1, Rls2, Rls3, Lls1, Lls2, Lls3,nthreads, debugging, STRATA_vals,KeepConstant);
    //
    if (verbose){
        end_point = system_clock::now();
        ending = time_point_cast<microseconds>(end_point).time_since_epoch().count();
        Rcout <<"df100 "<<(ending-start)<<" "<<halves<<" "<<iteration<<" "<<ind0<<",Update_Sides"<<endl;
        gibtime = system_clock::to_time_t(system_clock::now());
        Rcout << ctime(&gibtime) << endl;
        Rcout << "Wrapping up" << endl;
    }
    //
    Calc_LogLik_STRATA( nthreads, RiskFail, RiskGroup, totalnum, ntime, R, Rd, Rdd,RdR,RddR, Rls1, Rls2, Rls3, Lls1, Lls2, Lls3, Ll, Lld, Lldd, debugging, ties_method, STRATA_vals,KeepConstant);
    //
    if ((Ll_abs_best>0)||(Ll_abs_best < Ll[ind0])){
        Ll_abs_best = Ll[ind0];
        beta_abs_best = beta_c;
    }
    //
    if (verbose){
        Rcout << "df501 " << Ll_abs_best << endl;
        Rcout << "df504 ";//prints parameter values
        for (int ij=0;ij<totalnum;ij++){
            Rcout << beta_abs_best[ij] << " ";
        }
        Rcout << " " << endl;
    }
    //
    List para_list = List::create(_["Term_n"]=Term_n,_["tforms"]=tform); //stores the term information
    List control_list = List::create(_["Iteration"]=iteration, _["Maximum Step"]=abs_max, _["Derivative Limiting"]=Lld_worst); //stores the total number of iterations used
    //
    int kept_covs = totalnum - sum(KeepConstant); //does not base the standard deviation off of constant parameters
    NumericVector Lldd_vec(kept_covs * kept_covs);
    #pragma omp parallel for schedule(dynamic) num_threads(nthreads)
    for (int ijk=0;ijk<totalnum*(totalnum+1)/2;ijk++){
        int ij = 0;
        int jk = ijk;
        int pij_ind=-100;
        int pjk_ind=-100;
        while (jk>ij){
            ij++;
            jk-=ij;
        }
        if (KeepConstant[ij]==0){
            pij_ind = ij - sum(head(KeepConstant,ij));
            if (KeepConstant[jk]==0){
                pjk_ind = jk - sum(head(KeepConstant,jk));
                Lldd_vec[pij_ind * kept_covs + pjk_ind]=Lldd[ij*totalnum+jk];
            }
        }
    }
    for (int ijk=0;ijk<kept_covs*(kept_covs+1)/2;ijk++){
        int ij = 0;
        int jk = ijk;
        while (jk>ij){
            ij++;
            jk-=ij;
        }
        Lldd_vec[jk * kept_covs + ij]=Lldd_vec[ij * kept_covs + jk];
    }
    Lldd_vec.attr("dim") = Dimension(kept_covs, kept_covs);
    const Map<MatrixXd> Lldd_mat(as<Map<MatrixXd> >(Lldd_vec));
    //
    MatrixXd Lldd_inv = -1 * Lldd_mat.inverse().matrix(); //uses inverse information matrix to calculate the standard deviation
    VectorXd stdev = VectorXd::Zero(totalnum);
    for (int ij=0;ij<totalnum;ij++){
        if (KeepConstant[ij]==0){
            int pij_ind = ij - sum(head(KeepConstant,ij));
            stdev(ij) = sqrt(Lldd_inv(pij_ind,pij_ind));
        }
    }
    //
    //
    List res_list = List::create(_["LogLik"]=wrap(Ll[0]),_["First_Der"]=wrap(Lld),_["Second_Der"]=Lldd_vec,_["beta_0"]=wrap(beta_0) ,_["Standard_Deviation"]=wrap(stdev) ,_["AIC"]=2*(totalnum-accumulate(KeepConstant.begin(),KeepConstant.end(), 0.0))-2*Ll[0],_["Parameter_Lists"]=para_list,_["Control_List"]=control_list,_["Convgerged"]=convgd);
    // returns a list of results
    return res_list;
}


//' checks if the model is viable
//' \code{Check_Risk} Calculates risks and checks for negative values
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
//' @param     KeepConstant    vector identifying constant parameters
//' @param     term_tot    total number of terms
//' @param     nthreads number of threads to use
//'
//' @return True for viable point, False for negative error
// [[Rcpp::export]]
bool Check_Risk( IntegerVector Term_n, StringVector tform, NumericVector a_n,NumericMatrix x_all,IntegerVector dfc,int fir,string modelform, bool verbose, bool debugging, IntegerVector KeepConstant, int term_tot, int nthreads){
    ;
    //
    List temp_list = List::create(_["Status"]="FAILED"); //used as a dummy return value for code checking
    if (verbose){
        Rcout << "START_RISK_CHECK" << endl;
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
    Rcout.precision(7); //forces higher precision numbers printed to terminal
    // int nthreads = Eigen::nbThreads()-1; //stores how many threads are allocated
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
    //
    MatrixXd Te = MatrixXd::Zero(df0.rows(), 1); //preallocates matrix for column terms used for temporary storage
    MatrixXd R = MatrixXd::Zero(df0.rows(), 1); //preallocates matrix for Risks
    //
    MatrixXd Dose = MatrixXd::Constant(df0.rows(),term_tot,0.0); //Matrix of the total dose term values
    MatrixXd nonDose = MatrixXd::Constant(df0.rows(),term_tot,1.0); //Matrix of the total non-dose term values
    MatrixXd nonDose_LIN = MatrixXd::Constant(df0.rows(),term_tot,0.0); //matrix of Linear subterm values
    MatrixXd nonDose_PLIN = MatrixXd::Constant(df0.rows(),term_tot,1.0); //matrix of Loglinear subterm values
    MatrixXd nonDose_LOGLIN = MatrixXd::Constant(df0.rows(),term_tot,1.0); //matrix of Product linear subterm values
    MatrixXd TTerm = MatrixXd::Zero(Dose.rows(),Dose.cols()); //matrix of term values
    //
    if (verbose){
        Rcout << "starting subterms " << term_tot << endl;
    }
    // Calculates the subterm and term values
    Make_subterms_Single( totalnum, Term_n, tform, dfc, fir, T0, Dose, nonDose, TTerm,  nonDose_LIN, nonDose_PLIN, nonDose_LOGLIN ,beta_0, df0,nthreads, debugging,KeepConstant);
    // ---------------------------------------------------------
    // Prints off a series of calculations to check at what point values are changing
    // ---------------------------------------------------------
    //
    //
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
    Make_Risks_Single(modelform, tform, Term_n, totalnum, fir, T0, Te, R, Dose, nonDose, TTerm, nonDose_LIN, nonDose_PLIN, nonDose_LOGLIN, nthreads, debugging,KeepConstant);
    //
    // Removes infinite values
    //
    if (R.minCoeff()<0){
        Rcout << "A non-positive risk was detected: " << R.minCoeff() << endl;
        Rcout << "final failing values ";
        for (int ijk=0;ijk<totalnum;ijk++){
            Rcout << beta_0[ijk] << " ";
        }
        Rcout << " " << endl;
        //
        Rcout << "final failing terms ";
        for (int ijk=0;ijk<totalnum;ijk++){
            Rcout << tform[ijk] << " ";
        }
        Rcout << " " << endl;
        return FALSE;
    }
    return TRUE;
}
