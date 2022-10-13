#include <RcppEigen.h>
#include <RcppParallel.h>
#include <omp.h>
#include "PEANUT_MODEL.h"
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



// [[Rcpp::export]]
void Stress_Test(IntegerVector Term_n, StringVector tform, NumericVector a_n,IntegerVector dfc,NumericMatrix x_all, int fir, int der_iden,string modelform, List Control, NumericMatrix df_groups, NumericVector tu, IntegerVector KeepConstant, int term_tot, StringVector test_point){
    //----------------------------------------------------------------------------------------------------------------//
    //
    // Stress_Test runs through the standard calculations with a list of calculations to turn DEBUG on for
    //
    // change_all: list of 0/1s that control if a parameter will be kept constant
    // verbose: controls printing intermediate values and time-stamps
    // debugbing: controls printing additional time duration information
    // lr: learning rate for the second-derivative to first derivative ratio step
    // maxiter: maxmimum iterations used
    // halfmax: maximum number of half-steps used
    // epsilon: lowest acceptable maximum parameter change
    // dbeta_cap: learning rate for the first-derivative to log-likelihood ratio step
    // abs_max: starting upper limit on parameter change
    // dose_abs_max: starting upper limit on threshold parameter change
    // deriv_epsilon: threshold for lowest absolute value of derivative
    // ties_method: specifies breslow or efron tie method
    //
    //
    bool change_all = Control["change_all"];
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

    Stress_Run(Term_n, tform, a_n, x_all, dfc,fir, der_iden,modelform, lr, maxiter, halfmax, epsilon, dbeta_cap, abs_max,dose_abs_max, deriv_epsilon, df_groups, tu, change_all,verbose, debugging, KeepConstant, term_tot,test_point, ties_method );
    return;
}


void Stress_Run( IntegerVector Term_n, StringVector tform, NumericVector a_n,NumericMatrix x_all,IntegerVector dfc,int fir, int der_iden,string modelform, double lr, int maxiter, int halfmax, double epsilon, double dbeta_cap, double abs_max,double dose_abs_max, double deriv_epsilon, NumericMatrix df_groups, NumericVector tu ,bool change_all, bool verbose, bool debugging, IntegerVector KeepConstant, int term_tot, StringVector debug_checks, string ties_method){
srand (time(NULL));
    //
    // Runs through a single calculation with some functions printing additional information printed
    //
    if (verbose){
        cout << "START_NEW" << endl;
    }
    //
    // The labels for the calculations which can have additional information printed
    //
    StringVector aval_str = StringVector::create("MakeMatrix","MakeRisk","MakeGroup","CalcSide","CalcLL","CalcChange","UpdateRisk","IterMakeRisk","IterCalcSide","IterCalcLL");
    LogicalVector Debug_It = in(aval_str,debug_checks);
    cout << Debug_It << endl; //Booleans for which functions to check further
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
        cout << ctime(&gibtime) << endl;
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
    cout << "Term checked ";
    for (int ij=0;ij<totalnum;ij++){
        cout << Term_n[ij] << " ";
    }
    cout << " " << endl;
    //
    // cout.precision: controls the number of significant digits printed
    // nthreads: number of threads used for parallel operations
    //
    cout.precision(10); //forces higher precision numbers printed to terminal
    int nthreads = Eigen::nbThreads()-1; //stores how many threads are allocated

    //
    // Lld_worst: stores the highest magnitude log-likelihood derivative
    // totem: number of rows needed
    //
    double Lld_worst = 0.0; //stores derivative value used to determine if every parameter is near convergence
    double totem = df0.rows();//precalculates how many rows are needed
    //
    //
    if (verbose){
        end_point = system_clock::now();
        ending = time_point_cast<microseconds>(end_point).time_since_epoch().count();
        cout<<"df99,"<<(ending-start)<<",Starting"<<endl;
        gibtime = system_clock::to_time_t(system_clock::now());
        cout << ctime(&gibtime) << endl;
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
    double dint = dose_abs_max; //The amount of change used to calculate derivatives in threshold paramters
    int total_dose=0; //used later on for a section summing the dose terms
    //
    //
    // Calculates the subterm and term values
    cout << "starting subterms " << term_tot << endl;
    if (Debug_It[0]){
            Make_Subterms( totalnum, Term_n, tform, dfc, fir, T0, Td0, Tdd0, Dose, nonDose, TTerm,  nonDose_LIN, nonDose_PLIN, nonDose_LOGLIN ,beta_0, df0,dint,nthreads, TRUE);
    } else {
            Make_Subterms( totalnum, Term_n, tform, dfc, fir, T0, Td0, Tdd0, Dose, nonDose, TTerm,  nonDose_LIN, nonDose_PLIN, nonDose_LOGLIN ,beta_0, df0,dint,nthreads, FALSE);
    }
    // ---------------------------------------------------------
    // Prints off a series of calculations to check at what point values are changing
    // ---------------------------------------------------------
    int row_check = 10;
    if (verbose){
        cout << "values checked ";
        for (int ijk=0;ijk<totalnum;ijk++){
            cout << beta_0[ijk] << " ";
        }
        cout << " " << endl;
        cout << "sums checked ";
        for (int ijk=0;ijk<totalnum;ijk++){
            cout << T0.col(ijk).sum() << " ";
        }
        cout << " " << endl;
        cout << "derivs checked ";
        for (int ijk=0;ijk<totalnum;ijk++){
            cout << Td0.col(ijk).sum() << " ";
        }
        cout << " " << endl;
        cout << "second derivs checked ";
        for (int ijk=0;ijk<totalnum;ijk++){
            cout << Tdd0.col(ijk*(ijk+1)/2+ijk).sum() << " ";
        }
        cout << " " << endl;
        cout << "dose checked ";
        for (int ijk=0;ijk<term_tot;ijk++){
            cout << Dose.col(ijk).array().sum() << " ";
        }
        cout << " " << endl;
        cout << "non-dose checked ";
        for (int ijk=0;ijk<term_tot;ijk++){
            cout << nonDose.col(ijk).array().sum() << " ";
        }
        cout << " " << endl;
        cout << "LIN_non-dose checked ";
        for (int ijk=0;ijk<term_tot;ijk++){
            cout << nonDose_LIN.col(ijk).array().sum() << " ";
        }
        cout << " " << endl;
        cout << "PLIN_non-dose checked ";
        for (int ijk=0;ijk<term_tot;ijk++){
            cout << nonDose_PLIN.col(ijk).array().sum() << " ";
        }
        cout << " " << endl;
        cout << "LOGLIN_non-dose checked ";
        for (int ijk=0;ijk<term_tot;ijk++){
            cout << nonDose_LOGLIN.col(ijk).array().sum() << " ";
        }
        cout << " " << endl;
    }
    //
    List temp_list = List::create(_["betas"]=wrap(beta_0),_["Status"]="FAILED"); //used as a dummy return value for code checking
//    return temp_list;
    MatrixXd RdR = MatrixXd::Zero(df0.rows(), totalnum); //preallocates matrix for Risk to derivative ratios
    MatrixXd RddR = MatrixXd::Zero(df0.rows(), totalnum*(totalnum+1)/2); //preallocates matrix for Risk to second derivative ratios
    //
    if (verbose){
        end_point = system_clock::now();
        ending = time_point_cast<microseconds>(end_point).time_since_epoch().count();
        cout<<"df99,"<<(ending-start)<<",Prep_Terms"<<endl;
        gibtime = system_clock::to_time_t(system_clock::now());
        cout << ctime(&gibtime) << endl;
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
//    return temp_list;
    if (verbose){
        cout << "risk checked ";
        for (int ijk=0;ijk<totalnum;ijk++){
            cout << R.col(0).sum() << " ";
        }
        cout << " " << endl;
    //    return temp_list;
        cout << "risk1 checked ";
        for (int ijk=0;ijk<totalnum;ijk++){
            cout << Rd.col(ijk).sum() << " ";
        }
        cout << " " << endl;
        cout << "risk2 checked ";
        for (int ijk=0;ijk<totalnum;ijk++){
            cout << Rdd.col(ijk*(ijk+1)/2+ijk).sum() << " ";
        }
        cout << " " << endl;
        //
        end_point = system_clock::now();
        ending = time_point_cast<microseconds>(end_point).time_since_epoch().count();
        cout<<"df100 "<<(ending-start)<<" "<<0<<" "<<0<<" "<<-1<<",Prep_R"<<endl;
        gibtime = system_clock::to_time_t(system_clock::now());
        cout << ctime(&gibtime) << endl;
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
        cout<<"df100 "<<(ending-start)<<" "<<0<<" "<<0<<" "<<-1<<",Prep_List"<<endl;
        gibtime = system_clock::to_time_t(system_clock::now());
        cout << ctime(&gibtime) << endl;
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
        cout << "riskr checked ";
        for (int ijk=0;ijk<totalnum;ijk++){
            cout << Rls1.col(0).sum() << " ";
        }
        cout << " " << endl;
        cout << "risk1r checked ";
        for (int ijk=0;ijk<totalnum;ijk++){
            cout << Rls2.col(ijk).sum() << " ";
        }
        cout << " " << endl;
        cout << "risk2r checked ";
        for (int ijk=0;ijk<totalnum;ijk++){
            cout << Rls3.col(ijk*(ijk+1)/2+ijk).sum() << " ";
        }
        cout << " " << endl;
        //
        cout << "riskl checked ";
        for (int ijk=0;ijk<totalnum;ijk++){
            cout << Lls1.col(0).sum() << " ";
        }
        cout << " " << endl;
        cout << "risk1l checked ";
        for (int ijk=0;ijk<totalnum;ijk++){
            cout << Lls2.col(ijk).sum() << " ";
        }
        cout << " " << endl;
        cout << "risk2l checked ";
        for (int ijk=0;ijk<totalnum;ijk++){
            cout << Lls3.col(ijk*(ijk+1)/2+ijk).sum() << " ";
        }
        cout << " " << endl;
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
        cout<<"df100 "<<(ending-start)<<" "<<0<<" "<<0<<" "<<0<<",Calc"<<endl;//prints the time
        gibtime = system_clock::to_time_t(system_clock::now());
        cout << ctime(&gibtime) << endl;
        cout << "df101 ";//prints the log-likelihoods
        for (int ij=0;ij<totalnum;ij++){
            cout << Ll[ij] << " ";
        }
        cout << " " << endl;
        cout << "df102 ";//prints the first derivatives
        for (int ij=0;ij<totalnum;ij++){
            cout << Lld[ij] << " ";
        }
        cout << " " << endl;
        cout << "df103 ";//prints the second derivatives
        for (int ij=0;ij<totalnum;ij++){
            cout << Lldd[ij*totalnum+ij] << " ";
        }
        for (int ij=0;ij<totalnum;ij++){//locates highest magnitude derivative
            if (abs(Lld[ij]) > Lld_worst){
                Lld_worst = abs(Lld[ij]);
            }
        }
        cout << " " << endl;
        cout << "df104 ";//prints parameter values
        for (int ij=0;ij<totalnum;ij++){
            cout << beta_0[ij] << " ";
        }
        cout << " " << endl;
        cout << "df105 ";
        for (int ij=0;ij<totalnum;ij++){//prints the newton step value for zero derivative
            cout << Lld[ij]/Lldd[ij*totalnum+ij] << " ";
        }
        cout << " " << endl;
        cout << "df106 ";
        for (int ij=0;ij<totalnum;ij++){//prints the newton step value for zero log-likelihood
            cout << Ll[ij]/Lld[ij] << " ";
        }
        cout << " " << endl;
        cout << "df107 " << abs_max << " " << dose_abs_max << " " << Ll_comp[0] << " " << Ll_comp[1] << endl;//prints several convergence terms
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
            cout << "Starting Halves"<<endl;//prints the final changes for validation
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
                cout << "values checked ";
                for (int ijk=0;ijk<totalnum;ijk++){
                    cout << beta_c[ijk] << " ";
                }
                cout << " " << endl;
                cout << "sums checked ";
                for (int ijk=0;ijk<totalnum;ijk++){
                    cout << T0.col(ijk).sum() << " ";
                }
                cout << " " << endl;
                cout << "derivs checked ";
                for (int ijk=0;ijk<totalnum;ijk++){
                    cout << Td0.col(ijk).sum() << " ";
                }
                cout << " " << endl;
                cout << "second derivs checked ";
                for (int ijk=0;ijk<totalnum;ijk++){
                    cout << Tdd0.col(ijk).sum() << " ";
                }
                cout << " " << endl;
                cout << "dose checked ";
                for (int ijk=0;ijk<term_tot;ijk++){
                    cout << Dose.col(ijk).array().sum() << " ";
                }
                cout << " " << endl;
                cout << "non-dose checked ";
                for (int ijk=0;ijk<term_tot;ijk++){
                    cout << nonDose.col(ijk).array().sum() << " ";
                }
                cout << " " << endl;
                cout << "LIN_non-dose checked ";
                for (int ijk=0;ijk<term_tot;ijk++){
                    cout << nonDose_LIN.col(ijk).array().sum() << " ";
                }
                cout << " " << endl;
                cout << "PLIN_non-dose checked ";
                for (int ijk=0;ijk<term_tot;ijk++){
                    cout << nonDose_PLIN.col(ijk).array().sum() << " ";
                }
                cout << " " << endl;
                cout << "LOGLIN_non-dose checked ";
                for (int ijk=0;ijk<term_tot;ijk++){
                    cout << nonDose_LOGLIN.col(ijk).array().sum() << " ";
                }
                cout << " " << endl;
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
                cout << "risk checked ";
                for (int ijk=0;ijk<totalnum;ijk++){
                    cout << R.col(0).sum() << " ";
                }
                cout << " " << endl;
                cout << "risk1 checked ";
                for (int ijk=0;ijk<totalnum;ijk++){
                    cout << Rd.col(ijk).sum() << " ";
                }
                cout << " " << endl;
                cout << "risk2 checked ";
                for (int ijk=0;ijk<totalnum;ijk++){
                    cout << Rdd.col(ijk*(ijk+1)/2+ijk).sum() << " ";
                }
                cout << " " << endl;
                end_point = system_clock::now();
                ending = time_point_cast<microseconds>(end_point).time_since_epoch().count();
                cout<<"df100 "<<(ending-start)<<" "<<halves<<" "<<iteration<<" "<<ind0<<",Update_R"<<endl;
                //
                gibtime = system_clock::to_time_t(system_clock::now());
                cout << ctime(&gibtime) << endl;
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
                cout << "riskr checked ";
                for (int ijk=0;ijk<totalnum;ijk++){
                    cout << Rls1.col(0).sum() << " ";
                }
                cout << " " << endl;
                cout << "risk1r checked ";
                for (int ijk=0;ijk<totalnum;ijk++){
                    cout << Rls2.col(ijk).sum() << " ";
                }
                cout << " " << endl;
                cout << "risk2r checked ";
                for (int ijk=0;ijk<totalnum;ijk++){
                    cout << Rls3.col(ijk*(ijk+1)/2+ijk).sum() << " ";
                }
                cout << " " << endl;
                //
                //
                cout << "riskl checked ";
                for (int ijk=0;ijk<totalnum;ijk++){
                    cout << Lls1.col(0).sum() << " ";
                }
                cout << " " << endl;
                cout << "risk1l checked ";
                for (int ijk=0;ijk<totalnum;ijk++){
                    cout << Lls2.col(ijk).sum() << " ";
                }
                cout << " " << endl;
                cout << "risk2l checked ";
                for (int ijk=0;ijk<totalnum;ijk++){
                    cout << Lls3.col(ijk*(ijk+1)/2+ijk).sum() << " ";
                }
                cout << " " << endl;
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
                cout<<"df100 "<<(ending-start)<<" "<<halves<<" "<<iteration<<" "<<ind0<<",Update_calc"<<endl;
                gibtime = system_clock::to_time_t(system_clock::now());
                cout << ctime(&gibtime) << endl;
                cout << "df101 ";
                for (int ij=0;ij<totalnum;ij++){
                    cout << Ll[ij] << " ";
                }
                cout << " " << endl;
                cout << "df102 ";
                for (int ij=0;ij<totalnum;ij++){
                    cout << Lld[ij] << " ";
                }
                cout << " " << endl;
                cout << "df103 ";
                for (int ij=0;ij<totalnum;ij++){
                    cout << Lldd[ij*totalnum+ij] << " ";
                }
                for (int ij=0;ij<totalnum;ij++){
                    if (abs(Lld[ij]) > Lld_worst){
                        Lld_worst = abs(Lld[ij]);
                    }
                }
                cout << " " << endl;
                cout << "df104 ";
                for (int ij=0;ij<totalnum;ij++){
                    cout << beta_c[ij] << " ";
                }
                cout << " " << endl;
                cout << "df105 ";
                for (int ij=0;ij<totalnum;ij++){
                    cout << Lld[ij]/Lldd[ij*totalnum+ij] << " ";
                }
                cout << " " << endl;
                cout << "df106 ";
                for (int ij=0;ij<totalnum;ij++){
                    cout << Ll[ij]/Lld[ij] << " ";
                }
                cout << " " << endl;
                cout << "df107 " << abs_max << " " << dose_abs_max << " " << Ll_comp[0] << " " << Ll_comp[1] << endl;
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
            cout<<"df100 "<<(ending-start)<<" "<<halves<<" "<<iteration<<" "<<ind0<<",Recalc"<<endl;
            gibtime = system_clock::to_time_t(system_clock::now());
            cout << ctime(&gibtime) << endl;
            cout << "df101 ";
            for (int ij=0;ij<totalnum;ij++){
                cout << Ll[ij] << " ";
            }
            cout << " " << endl;
            cout << "df102 ";
            for (int ij=0;ij<totalnum;ij++){
                cout << Lld[ij] << " ";
            }
            cout << " " << endl;
            cout << "df103 ";
            for (int ij=0;ij<totalnum;ij++){
                cout << Lldd[ij*totalnum+ij] << " ";
            }
            cout << " " << endl;
            cout << "df104 ";
            for (int ij=0;ij<totalnum;ij++){
                cout << beta_c[ij] << " ";
            }
            cout << " " << endl;
            cout << "Finshed iteration" << endl;
            cout << "df105 ";
            for (int ij=0;ij<totalnum;ij++){
                cout << Lld[ij]/Lldd[ij*totalnum+ij] << " ";
            }
            cout << " " << endl;
            cout << "df106 ";
            for (int ij=0;ij<totalnum;ij++){
                cout << Ll[ij]/Lld[ij] << " ";
            }
            cout << " " << endl;
            cout << "df107 " << abs_max << " " << dose_abs_max << " " << Ll_comp[0] << " " << Ll_comp[1] << endl;
        }
    }
    //
    // Does not return anything
    //
    return;
}

// [[Rcpp::export]]
List peanut_transition(IntegerVector Term_n, StringVector tform, NumericVector a_n,IntegerVector dfc,NumericMatrix x_all, int fir, int der_iden,string modelform, List Control, NumericMatrix df_groups, NumericVector tu, IntegerVector KeepConstant, int term_tot){
    // Term_n,STerm_n,tform,a_n,dfc,x_all, fir, modelform,length(tu), control,as.matrix(df[,..ce]),tu,keep_constant,term_tot
    //----------------------------------------------------------------------------------------------------------------//
    // change_all: list of 0/1s that control if a parameter will be kept constant
    // verbose: controls printing intermediate values and time-stamps
    // debugbing: controls printing additional time duration information
    // lr: learning rate for the second-derivative to first derivative ratio step
    // maxiter: maxmimum iterations used
    // halfmax: maximum number of half-steps used
    // epsilon: lowest acceptable maximum parameter change
    // dbeta_cap: learning rate for the first-derivative to log-likelihood ratio step
    // abs_max: starting upper limit on parameter change
    // dose_abs_max: starting upper limit on threshold parameter change
    // deriv_epsilon: threshold for lowest absolute value of derivative
    // ties_method: specifies breslow or efron tie method
    //
    bool change_all = Control["change_all"];
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
    // Performs regression
    //----------------------------------------------------------------------------------------------------------------//
    List res = LogLik_PEANUT(Term_n, tform, a_n, x_all, dfc,fir, der_iden,modelform, lr, maxiter, halfmax, epsilon, dbeta_cap, abs_max,dose_abs_max, deriv_epsilon, df_groups, tu, change_all,verbose, debugging, KeepConstant, term_tot, ties_method);
    //----------------------------------------------------------------------------------------------------------------//
    return res;
}

/*
// --------------------------------------------------------------------------------------------------------------- //
// bound calculations are not currently supported, but the code is kept for later reincorporation
// --------------------------------------------------------------------------------------------------------------- //
// [[Rcpp::export]]
List peanut_bounds_transition(double q1, NumericVector a_lin,NumericVector a_loglin,NumericVector a_plin, NumericMatrix x_lin, NumericMatrix x_loglin, NumericMatrix x_plin, NumericMatrix x_dose,int fir,string modelform,int ntime, NumericVector include_bool, List Control, List Dose_paras, NumericMatrix df_groups, NumericVector tu){
    //----------------------------------------------------------------------------------------------------------------//
    Map<VectorXd> beta_lin(as<Map<VectorXd> >(a_lin));
    Map<VectorXd> beta_loglin(as<Map<VectorXd> >(a_loglin));
    Map<VectorXd> beta_plin(as<Map<VectorXd> >(a_plin));
    const Map<MatrixXd> df_lin(as<Map<MatrixXd> >(x_lin));
    const Map<MatrixXd> df_loglin(as<Map<MatrixXd> >(x_loglin));
    const Map<MatrixXd> df_plin(as<Map<MatrixXd> >(x_plin));
    const Map<MatrixXd> df_dose(as<Map<MatrixXd> >(x_dose));
    // Converts from Rcpp types to efficient Eigen types
    double lr = Control["lr"];
    bool verbose = Control["verbose"]
    int maxiter = Control["maxiter"];
	int halfmax = Control["halfmax"];
	double epsilon = Control["epsilon"];
	double dbeta_cap = Control["dbeta_max"];
	double abs_max = Control["abs_max"];
	double dose_abs_max = Control["dose_abs_max"];
	double deriv_epsilon =Control["deriv_epsilon"];
	List beta_loglin_slope = Dose_paras["beta_loglin_slope"];
    List beta_loglin_top  = Dose_paras["beta_loglin_top"];
    List beta_lin_slope  = Dose_paras["beta_lin_slope"];
    List beta_lin_int  = Dose_paras["beta_lin_int"];
    List beta_quad  = Dose_paras["beta_quad"];
    List beta_step_slope  = Dose_paras["beta_step_slope"];
    List beta_step_int  = Dose_paras["beta_step_int"];
    //
    //----------------------------------------------------------------------------------------------------------------//
    List res = LogLik_Bounds(q1, beta_lin,beta_loglin,beta_plin,df_lin,df_loglin,df_plin, df_dose,fir,modelform,ntime,include_bool, lr, maxiter, halfmax, epsilon, dbeta_cap, abs_max, dose_abs_max, deriv_epsilon,beta_loglin_slope, beta_loglin_top , beta_lin_slope , beta_lin_int , beta_quad , beta_step_slope , beta_step_int, df_groups,tu,verbose);
    //----------------------------------------------------------------------------------------------------------------//
    return res;
}
*/

// [[Rcpp::export]]
List peanut_plot(IntegerVector Term_n, StringVector tform, NumericVector a_n,IntegerVector dfc,NumericMatrix x_all, int fir, int der_iden,string modelform, List Control, NumericMatrix df_groups, NumericVector tu, IntegerVector KeepConstant, int term_tot, vector<string> Plot_Type ,int uniq_v){
    //----------------------------------------------------------------------------------------------------------------//
    // change_all: list of 0/1s that control if a parameter will be kept constant
    // verbose: controls printing intermediate values and time-stamps
    // debugbing: controls printing additional time duration information
    // lr: learning rate for the second-derivative to first derivative ratio step
    // maxiter: maxmimum iterations used
    // halfmax: maximum number of half-steps used
    // epsilon: lowest acceptable maximum parameter change
    // dbeta_cap: learning rate for the first-derivative to log-likelihood ratio step
    // abs_max: starting upper limit on parameter change
    // dose_abs_max: starting upper limit on threshold parameter change
    // deriv_epsilon: threshold for lowest absolute value of derivative
    // ties_method: specifies breslow or efron tie method
    double lr = Control["lr"];
    bool verbose = Control["verbose"];
    bool debugging = FALSE;
    int maxiter = Control["maxiter"];
	int halfmax = Control["halfmax"];
	double epsilon = Control["epsilon"];
	double dbeta_cap = Control["dbeta_max"];
	double abs_max = Control["abs_max"];
	double dose_abs_max = Control["dose_abs_max"];
	double deriv_epsilon =Control["deriv_epsilon"];
	string ties_method =Control["ties"];
    List res;
    // there are two types of plots that can be generated, survival curve and risk by covariate value
    //----------------------------------------------------------------------------------------------------------------//
    if (Plot_Type[0]=="SURV"){
        res = PEANUT_PLOT_SURV(Term_n, tform, a_n, x_all, dfc,fir, der_iden,modelform, abs_max,dose_abs_max, df_groups, tu,verbose, debugging, KeepConstant, term_tot);
    }else if (Plot_Type[0]=="RISK"){
        res = PEANUT_PLOT_RISK(Term_n, tform, a_n, x_all, dfc,fir, der_iden,modelform, abs_max,dose_abs_max, df_groups, tu,verbose, debugging, KeepConstant, term_tot, uniq_v);
    } else {
        throw invalid_argument("Invalid plot type");
    }
    return res;
}

// [[Rcpp::export]]
NumericMatrix peanut_schoenfeld_transition(IntegerVector Term_n, StringVector tform, NumericVector a_n,IntegerVector dfc,NumericMatrix x_all, int fir, int der_iden,string modelform, List Control, NumericMatrix df_groups, NumericVector tu, IntegerVector KeepConstant, int term_tot){
        //----------------------------------------------------------------------------------------------------------------//
    // Converts from Rcpp types to efficient Eigen types
    //
    //----------------------------------------------------------------------------------------------------------------//
    //
    // change_all: list of 0/1s that control if a parameter will be kept constant
    // verbose: controls printing intermediate values and time-stamps
    // debugbing: controls printing additional time duration information
    // abs_max: starting upper limit on parameter change
    // dose_abs_max: starting upper limit on threshold parameter change
    //
    bool change_all = Control["change_all"];
    bool verbose = Control["verbose"];
    bool debugging = FALSE;
	double abs_max = Control["abs_max"];
	double dose_abs_max = Control["dose_abs_max"];
	string ties_method =Control["ties"];
	//
	// performs schoenfeld residual calculation
    NumericMatrix res = Schoenfeld_PEANUT(Term_n, tform, a_n, x_all, dfc,fir, der_iden,modelform, abs_max,dose_abs_max, df_groups, tu,verbose, debugging, KeepConstant, term_tot, ties_method);
    //----------------------------------------------------------------------------------------------------------------//
    return res;
}

// [[Rcpp::export]]
List amfit_transition(NumericMatrix dfe, IntegerVector Term_n, StringVector tform, NumericVector a_n,IntegerVector dfc,NumericMatrix x_all, int fir, int der_iden,string modelform, List Control, IntegerVector KeepConstant, int term_tot){
    //----------------------------------------------------------------------------------------------------------------//
    const Map<MatrixXd> PyrC(as<Map<MatrixXd> >(dfe));
    //
    //
    // change_all: list of 0/1s that control if a parameter will be kept constant
    // verbose: controls printing intermediate values and time-stamps
    // debugbing: controls printing additional time duration information
    // lr: learning rate for the second-derivative to first derivative ratio step
    // maxiter: maxmimum iterations used
    // halfmax: maximum number of half-steps used
    // epsilon: lowest acceptable maximum parameter change
    // dbeta_cap: learning rate for the first-derivative to log-likelihood ratio step
    // abs_max: starting upper limit on parameter change
    // dose_abs_max: starting upper limit on threshold parameter change
    // deriv_epsilon: threshold for lowest absolute value of derivative
    //
    bool change_all = Control["change_all"];
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
    //----------------------------------------------------------------------------------------------------------------//
    List res = LogLik_AMFIT(PyrC,Term_n, tform, a_n, x_all, dfc,fir, der_iden,modelform, lr, maxiter, halfmax, epsilon, dbeta_cap, abs_max,dose_abs_max, deriv_epsilon, change_all,verbose, debugging, KeepConstant, term_tot);
    //----------------------------------------------------------------------------------------------------------------//
    return res;
}

List LogLik_PEANUT( IntegerVector Term_n, StringVector tform, NumericVector a_n,NumericMatrix x_all,IntegerVector dfc,int fir, int der_iden,string modelform, double lr, int maxiter, int halfmax, double epsilon, double dbeta_cap, double abs_max,double dose_abs_max, double deriv_epsilon, NumericMatrix df_groups, NumericVector tu ,bool change_all, bool verbose, bool debugging, IntegerVector KeepConstant, int term_tot, string ties_method){
    srand (time(NULL));
    //
    List temp_list = List::create(_["Status"]="FAILED"); //used as a dummy return value for code checking
    if (verbose){
        cout << "START_NEW" << endl;
    }
    time_point<system_clock> start_point, end_point;
    start_point = system_clock::now();
    auto start = time_point_cast<microseconds>(start_point).time_since_epoch().count();
    end_point = system_clock::now();
    auto ending = time_point_cast<microseconds>(end_point).time_since_epoch().count(); //The time duration is tracked
    //
    auto gibtime = system_clock::to_time_t(system_clock::now());
    if (verbose){
        cout << ctime(&gibtime) << endl;
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
    cout << "Term checked ";
    for (int ij=0;ij<totalnum;ij++){
        cout << Term_n[ij] << " ";
    }
    cout << " " << endl;
    //
    // cout.precision: controls the number of significant digits printed
    // nthreads: number of threads used for parallel operations
    //
    cout.precision(7); //forces higher precision numbers printed to terminal
    int nthreads = Eigen::nbThreads()-1; //stores how many threads are allocated
    //
    //
    // Lld_worst: stores the highest magnitude log-likelihood derivative
    // totem: number of rows needed
    //
    //
    double Lld_worst = 0.0; //stores derivative value used to determine if every parameter is near convergence
    double totem = df0.rows();//precalculates how many rows are needed
    //
    if (verbose){
        end_point = system_clock::now();
        ending = time_point_cast<microseconds>(end_point).time_since_epoch().count();
        cout<<"df99,"<<(ending-start)<<",Starting"<<endl;
        gibtime = system_clock::to_time_t(system_clock::now());
        cout << ctime(&gibtime) << endl;
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
    double dint = dose_abs_max; //The amount of change used to calculate derivatives in threshold paramters
    int total_dose=0; //used later on for a section summing the dose terms
    //
    if (verbose){
        cout << "starting subterms " << term_tot << endl;
    }
    // Calculates the subterm and term values
    Make_Subterms( totalnum, Term_n, tform, dfc, fir, T0, Td0, Tdd0, Dose, nonDose, TTerm,  nonDose_LIN, nonDose_PLIN, nonDose_LOGLIN ,beta_0, df0,dint,nthreads, debugging);
    // ---------------------------------------------------------
    // Prints off a series of calculations to check at what point values are changing
    // ---------------------------------------------------------
    int row_check = 10;
    if (verbose){
        cout << "values checked ";
        for (int ijk=0;ijk<totalnum;ijk++){
            cout << beta_0[ijk] << " ";
        }
        cout << " " << endl;
        cout << "sums checked ";
        for (int ijk=0;ijk<totalnum;ijk++){
            cout << T0.col(ijk).sum() << " ";
        }
        cout << " " << endl;
        cout << "derivs checked ";
        for (int ijk=0;ijk<totalnum;ijk++){
            cout << Td0.col(ijk).sum() << " ";
        }
        cout << " " << endl;
        cout << "second derivs checked ";
        for (int ijk=0;ijk<totalnum;ijk++){
            cout << Tdd0.col(ijk*(ijk+1)/2+ijk).sum() << " ";
        }
        cout << " " << endl;
        cout << "dose checked ";
        for (int ijk=0;ijk<term_tot;ijk++){
            cout << Dose.col(ijk).array().sum() << " ";
        }
        cout << " " << endl;
        cout << "non-dose checked ";
        for (int ijk=0;ijk<term_tot;ijk++){
            cout << nonDose.col(ijk).array().sum() << " ";
        }
        cout << " " << endl;
        cout << "LIN_non-dose checked ";
        for (int ijk=0;ijk<term_tot;ijk++){
            cout << nonDose_LIN.col(ijk).array().sum() << " ";
        }
        cout << " " << endl;
        cout << "PLIN_non-dose checked ";
        for (int ijk=0;ijk<term_tot;ijk++){
            cout << nonDose_PLIN.col(ijk).array().sum() << " ";
        }
        cout << " " << endl;
        cout << "LOGLIN_non-dose checked ";
        for (int ijk=0;ijk<term_tot;ijk++){
            cout << nonDose_LOGLIN.col(ijk).array().sum() << " ";
        }
        cout << " " << endl;
    }
    //
    MatrixXd RdR = MatrixXd::Zero(df0.rows(), totalnum); //preallocates matrix for Risk to derivative ratios
    MatrixXd RddR = MatrixXd::Zero(df0.rows(), totalnum*(totalnum+1)/2); //preallocates matrix for Risk to second derivative ratios
    //
    if (verbose){
        end_point = system_clock::now();
        ending = time_point_cast<microseconds>(end_point).time_since_epoch().count();
        cout<<"df99,"<<(ending-start)<<",Prep_Terms"<<endl;
        gibtime = system_clock::to_time_t(system_clock::now());
        cout << ctime(&gibtime) << endl;
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
                cout << ijk << " had a non-positive term" << endl;
            }
        }
        for (int ijk=0;ijk<term_tot;ijk++){
            if (nonDose_LIN.col(ijk).minCoeff()<=0){
                cout << ijk << " had a non-positive Lin term" << endl;
            }
            if (nonDose_LOGLIN.col(ijk).minCoeff()<=0){
                cout << ijk << " had a non-positive loglin term" << endl;
            }
            if (nonDose_PLIN.col(ijk).minCoeff()<=0){
                cout << ijk << " had a non-positive plin term" << endl;
            }
            if (Dose.col(ijk).minCoeff()<=0){
                cout << ijk << " had a non-positive dose term" << endl;
            }
            if (nonDose.col(ijk).minCoeff()<=0){
                cout << ijk << " had a non-positive nondose term" << endl;
            }
        }
        cout << R.sum() << endl;
        cout << "A non-positive risk was detected: " << R.minCoeff() << endl;
        return temp_list;
    }
    //
    if (verbose){
        cout << "risk checked ";
        for (int ijk=0;ijk<totalnum;ijk++){
            cout << R.col(0).sum() << " ";
        }
        cout << " " << endl;
        cout << "risk1 checked ";
        for (int ijk=0;ijk<totalnum;ijk++){
            cout << Rd.col(ijk).sum() << " ";
        }
        cout << " " << endl;
        cout << "risk2 checked ";
        for (int ijk=0;ijk<totalnum;ijk++){
            cout << Rdd.col(ijk*(ijk+1)/2+ijk).sum() << " ";
        }
        cout << " " << endl;
        //
        end_point = system_clock::now();
        ending = time_point_cast<microseconds>(end_point).time_since_epoch().count();
        cout<<"df100 "<<(ending-start)<<" "<<0<<" "<<0<<" "<<-1<<",Prep_R"<<endl;
        gibtime = system_clock::to_time_t(system_clock::now());
        cout << ctime(&gibtime) << endl;
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
        cout<<"df100 "<<(ending-start)<<" "<<0<<" "<<0<<" "<<-1<<",Prep_List"<<endl;
        gibtime = system_clock::to_time_t(system_clock::now());
        cout << ctime(&gibtime) << endl;
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
        cout << "riskr checked ";
        for (int ijk=0;ijk<totalnum;ijk++){
            cout << Rls1.col(0).sum() << " ";
        }
        cout << " " << endl;
        cout << "risk1r checked ";
        for (int ijk=0;ijk<totalnum;ijk++){
            cout << Rls2.col(ijk).sum() << " ";
        }
        cout << " " << endl;
        cout << "risk2r checked ";
        for (int ijk=0;ijk<totalnum;ijk++){
            cout << Rls3.col(ijk*(ijk+1)/2+ijk).sum() << " ";
        }
        cout << " " << endl;
        //
        cout << "riskl checked ";
        for (int ijk=0;ijk<totalnum;ijk++){
            cout << Lls1.col(0).sum() << " ";
        }
        cout << " " << endl;
        cout << "risk1l checked ";
        for (int ijk=0;ijk<totalnum;ijk++){
            cout << Lls2.col(ijk).sum() << " ";
        }
        cout << " " << endl;
        cout << "risk2l checked ";
        for (int ijk=0;ijk<totalnum;ijk++){
            cout << Lls3.col(ijk*(ijk+1)/2+ijk).sum() << " ";
        }
        cout << " " << endl;
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
        cout<<"df100 "<<(ending-start)<<" "<<0<<" "<<0<<" "<<0<<",Calc"<<endl;//prints the time
        gibtime = system_clock::to_time_t(system_clock::now());
        cout << ctime(&gibtime) << endl;
        cout << "df101 ";//prints the log-likelihoods
        for (int ij=0;ij<totalnum;ij++){
            cout << Ll[ij] << " ";
        }
        cout << " " << endl;
        cout << "df102 ";//prints the first derivatives
        for (int ij=0;ij<totalnum;ij++){
            cout << Lld[ij] << " ";
        }
        cout << " " << endl;
        cout << "df103 ";//prints the second derivatives
        for (int ij=0;ij<totalnum;ij++){
            cout << Lldd[ij*totalnum+ij] << " ";
        }
        for (int ij=0;ij<totalnum;ij++){//locates highest magnitude derivative
            if (abs(Lld[ij]) > Lld_worst){
                Lld_worst = abs(Lld[ij]);
            }
        }
        cout << " " << endl;
        cout << "df104 ";//prints parameter values
        for (int ij=0;ij<totalnum;ij++){
            cout << beta_0[ij] << " ";
        }
        cout << " " << endl;
        cout << "df105 ";
        for (int ij=0;ij<totalnum;ij++){//prints the newton step value for zero derivative
            cout << Lld[ij]/Lldd[ij*totalnum+ij] << " ";
        }
        cout << " " << endl;
        cout << "df106 ";
        for (int ij=0;ij<totalnum;ij++){//prints the newton step value for zero log-likelihood
            cout << Ll[ij]/Lld[ij] << " ";
        }
        cout << " " << endl;
        cout << "df107 " << abs_max << " " << dose_abs_max << " " << Ll_comp[0] << " " << Ll_comp[1] << endl;//prints several convergence terms
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
            cout << "Starting Halves"<<endl;//prints the final changes for validation
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
                cout << "values checked ";
                for (int ijk=0;ijk<totalnum;ijk++){
                    cout << beta_c[ijk] << " ";
                }
                cout << " " << endl;
                cout << "sums checked ";
                for (int ijk=0;ijk<totalnum;ijk++){
                    cout << T0.col(ijk).sum() << " ";
                }
                cout << " " << endl;
                cout << "derivs checked ";
                for (int ijk=0;ijk<totalnum;ijk++){
                    cout << Td0.col(ijk).sum() << " ";
                }
                cout << " " << endl;
                cout << "second derivs checked ";
                for (int ijk=0;ijk<totalnum;ijk++){
                    cout << Tdd0.col(ijk*(ijk+1)/2+ijk).sum() << " ";
                }
                cout << " " << endl;
                cout << "dose checked ";
                for (int ijk=0;ijk<term_tot;ijk++){
                    cout << Dose.col(ijk).array().sum() << " ";
                }
                cout << " " << endl;
                cout << "non-dose checked ";
                for (int ijk=0;ijk<term_tot;ijk++){
                    cout << nonDose.col(ijk).array().sum() << " ";
                }
                cout << " " << endl;
                cout << "LIN_non-dose checked ";
                for (int ijk=0;ijk<term_tot;ijk++){
                    cout << nonDose_LIN.col(ijk).array().sum() << " ";
                }
                cout << " " << endl;
                cout << "PLIN_non-dose checked ";
                for (int ijk=0;ijk<term_tot;ijk++){
                    cout << nonDose_PLIN.col(ijk).array().sum() << " ";
                }
                cout << " " << endl;
                cout << "LOGLIN_non-dose checked ";
                for (int ijk=0;ijk<term_tot;ijk++){
                    cout << nonDose_LOGLIN.col(ijk).array().sum() << " ";
                }
                cout << " " << endl;
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
                cout << "risk checked ";
                for (int ijk=0;ijk<totalnum;ijk++){
                    cout << R.col(0).sum() << " ";
                }
                cout << " " << endl;
                cout << "risk1 checked ";
                for (int ijk=0;ijk<totalnum;ijk++){
                    cout << Rd.col(ijk).sum() << " ";
                }
                cout << " " << endl;
                cout << "risk2 checked ";
                for (int ijk=0;ijk<totalnum;ijk++){
                    cout << Rdd.col(ijk*(ijk+1)/2+ijk).sum() << " ";
                }
                cout << " " << endl;
                end_point = system_clock::now();
                ending = time_point_cast<microseconds>(end_point).time_since_epoch().count();
                cout<<"df100 "<<(ending-start)<<" "<<halves<<" "<<iteration<<" "<<ind0<<",Update_R"<<endl;
                //
                gibtime = system_clock::to_time_t(system_clock::now());
                cout << ctime(&gibtime) << endl;
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
                cout << "riskr checked ";
                for (int ijk=0;ijk<totalnum;ijk++){
                    cout << Rls1.col(0).sum() << " ";
                }
                cout << " " << endl;
                cout << "risk1r checked ";
                for (int ijk=0;ijk<totalnum;ijk++){
                    cout << Rls2.col(ijk).sum() << " ";
                }
                cout << " " << endl;
                cout << "risk2r checked ";
                for (int ijk=0;ijk<totalnum;ijk++){
                    cout << Rls3.col(ijk*(ijk+1)/2+ijk).sum() << " ";
                }
                cout << " " << endl;
                //
                //
                cout << "riskl checked ";
                for (int ijk=0;ijk<totalnum;ijk++){
                    cout << Lls1.col(0).sum() << " ";
                }
                cout << " " << endl;
                cout << "risk1l checked ";
                for (int ijk=0;ijk<totalnum;ijk++){
                    cout << Lls2.col(ijk).sum() << " ";
                }
                cout << " " << endl;
                cout << "risk2l checked ";
                for (int ijk=0;ijk<totalnum;ijk++){
                    cout << Lls3.col(ijk*(ijk+1)/2+ijk).sum() << " ";
                }
                cout << " " << endl;
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
                cout<<"df100 "<<(ending-start)<<" "<<halves<<" "<<iteration<<" "<<ind0<<",Update_calc"<<endl;
                gibtime = system_clock::to_time_t(system_clock::now());
                cout << ctime(&gibtime) << endl;
                cout << "df101 ";
                for (int ij=0;ij<totalnum;ij++){
                    cout << Ll[ij] << " ";
                }
                cout << " " << endl;
                cout << "df102 ";
                for (int ij=0;ij<totalnum;ij++){
                    cout << Lld[ij] << " ";
                }
                cout << " " << endl;
                cout << "df103 ";
                for (int ij=0;ij<totalnum;ij++){
                    cout << Lldd[ij*totalnum+ij] << " ";
                }
                for (int ij=0;ij<totalnum;ij++){
                    if (abs(Lld[ij]) > Lld_worst){
                        Lld_worst = abs(Lld[ij]);
                    }
                }
                cout << " " << endl;
                cout << "df104 ";
                for (int ij=0;ij<totalnum;ij++){
                    cout << beta_c[ij] << " ";
                }
                cout << " " << endl;
                cout << "df105 ";
                for (int ij=0;ij<totalnum;ij++){
                    cout << Lld[ij]/Lldd[ij*totalnum+ij] << " ";
                }
                cout << " " << endl;
                cout << "df106 ";
                for (int ij=0;ij<totalnum;ij++){
                    cout << Ll[ij]/Lld[ij] << " ";
                }
                cout << " " << endl;
                cout << "df107 " << abs_max << " " << dose_abs_max << " " << Ll_comp[0] << " " << Ll_comp[1] << endl;
            }
            #pragma omp parallel for num_threads(nthreads)
            for (int ijk=0;ijk<totalnum;ijk++){//totalnum*(totalnum+1)/2
                beta_0[ijk] = beta_c[ijk];
            }
        }
        if (beta_best!=beta_c){//if the risk matrices aren't the optimal values, then they must be recalculated
            if (verbose){
                cout << "Changing back to best"<<endl;
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
                cout << "values checked ";
                for (int ijk=0;ijk<totalnum;ijk++){
                    cout << beta_c[ijk] << " ";
                }
                cout << " " << endl;
                //
                //
                cout << "sums checked ";
                for (int ijk=0;ijk<totalnum;ijk++){
                    cout << T0.col(ijk).sum() << " ";
                }
                cout << " " << endl;
                cout << "derivs checked ";
                for (int ijk=0;ijk<totalnum;ijk++){
                    cout << Td0.col(ijk).sum() << " ";
                }
                cout << " " << endl;
                cout << "second derivs checked ";
                for (int ijk=0;ijk<totalnum;ijk++){
                    cout << Tdd0.col(ijk*(ijk+1)/2+ijk).sum() << " ";
                }
                cout << " " << endl;
                cout << "dose checked ";
                for (int ijk=0;ijk<1;ijk++){
                    cout << Dose.array().sum() << " ";
                }
                cout << " " << endl;
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
                        cout << ijk << " had a non-positive term" << endl;
                    }
                }
                for (int ijk=0;ijk<term_tot;ijk++){
                    if (nonDose_LIN.col(ijk).minCoeff()<=0){
                        cout << ijk << " had a non-positive Lin term" << endl;
                    }
                    if (nonDose_LOGLIN.col(ijk).minCoeff()<=0){
                        cout << ijk << " had a non-positive loglin term" << endl;
                    }
                    if (nonDose_PLIN.col(ijk).minCoeff()<=0){
                        cout << ijk << " had a non-positive plin term" << endl;
                    }
                    if (Dose.col(ijk).minCoeff()<=0){
                        cout << ijk << " had a non-positive dose term" << endl;
                    }
                    if (nonDose.col(ijk).minCoeff()<=0){
                        cout << ijk << " had a non-positive nondose term" << endl;
                    }
                }
                cout << R.sum() << endl;
                cout << "A non-positive risk was detected: " << R.minCoeff() << endl;
                return temp_list;
            }
            if (verbose){
                cout << "risk checked ";
                for (int ijk=0;ijk<totalnum;ijk++){
                    cout << R.col(0).sum() << " ";
                }
                cout << " " << endl;
                cout << "risk1 checked ";
                for (int ijk=0;ijk<totalnum;ijk++){
                    cout << Rd.col(ijk).sum() << " ";
                }
                cout << " " << endl;
                cout << "risk2 checked ";
                for (int ijk=0;ijk<totalnum;ijk++){
                    cout << Rdd.col(ijk*(ijk+1)/2+ijk).sum() << " ";
                }
                cout << " " << endl;
                end_point = system_clock::now();
                ending = time_point_cast<microseconds>(end_point).time_since_epoch().count();
                cout<<"df100 "<<(ending-start)<<" "<<halves<<" "<<iteration<<" "<<ind0<<",Update_R"<<endl;
                //
                gibtime = system_clock::to_time_t(system_clock::now());
                cout << ctime(&gibtime) << endl;
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
                cout << "riskr checked ";
                for (int ijk=0;ijk<totalnum;ijk++){
                    cout << Rls1.col(0).sum() << " ";
                }
                cout << " " << endl;
                cout << "risk1r checked ";
                for (int ijk=0;ijk<totalnum;ijk++){
                    cout << Rls2.col(ijk).sum() << " ";
                }
                cout << " " << endl;
                cout << "risk2r checked ";
                for (int ijk=0;ijk<totalnum;ijk++){
                    cout << Rls3.col(ijk*(ijk+1)/2+ijk).sum() << " ";
                }
                cout << " " << endl;
                //
                //
                cout << "riskl checked ";
                for (int ijk=0;ijk<totalnum;ijk++){
                    cout << Lls1.col(0).sum() << " ";
                }
                cout << " " << endl;
                cout << "risk1l checked ";
                for (int ijk=0;ijk<totalnum;ijk++){
                    cout << Lls2.col(ijk).sum() << " ";
                }
                cout << " " << endl;
                cout << "risk2l checked ";
                for (int ijk=0;ijk<totalnum;ijk++){
                    cout << Lls3.col(ijk*(ijk+1)/2+ijk).sum() << " ";
                }
                cout << " " << endl;
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
            cout<<"df100 "<<(ending-start)<<" "<<halves<<" "<<iteration<<" "<<ind0<<",Recalc"<<endl;
            gibtime = system_clock::to_time_t(system_clock::now());
            cout << ctime(&gibtime) << endl;
            cout << "df101 ";
            for (int ij=0;ij<totalnum;ij++){
                cout << Ll[ij] << " ";
            }
            cout << " " << endl;
            cout << "df102 ";
            for (int ij=0;ij<totalnum;ij++){
                cout << Lld[ij] << " ";
            }
            cout << " " << endl;
            cout << "df103 ";
            for (int ij=0;ij<totalnum;ij++){
                cout << Lldd[ij*totalnum+ij] << " ";
            }
            cout << " " << endl;
            cout << "df104 ";
            for (int ij=0;ij<totalnum;ij++){
                cout << beta_c[ij] << " ";
            }
            cout << " " << endl;
            cout << "Finshed iteration" << endl;
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
        cout << "Wrapping up" << endl;
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


List PEANUT_PLOT_SURV(IntegerVector Term_n, StringVector tform, NumericVector a_n,NumericMatrix x_all,IntegerVector dfc,int fir, int der_iden,string modelform, double abs_max,double dose_abs_max, NumericMatrix df_groups, NumericVector tu , bool verbose, bool debugging, IntegerVector KeepConstant, int term_tot){
    srand (time(NULL));
    //
    // Calculates the baseline hazard
    //
    using namespace std::chrono;
    if (verbose){
        cout << "START_NEW" << endl;
    }
    time_point<system_clock> start_point, end_point;
    start_point = system_clock::now();
    auto start = time_point_cast<microseconds>(start_point).time_since_epoch().count();
    end_point = system_clock::now();
    auto ending = time_point_cast<microseconds>(end_point).time_since_epoch().count(); //The time duration is tracked
    //
    auto gibtime = system_clock::to_time_t(system_clock::now());
    if (verbose){
        cout << ctime(&gibtime) << endl;
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
    cout.precision(10); //forces higher precision numbers printed to terminal
    int nthreads = Eigen::nbThreads()-1; //stores how many threads are allocated
    //
    //
    //
    // Lld_worst: stores the highest magnitude log-likelihood derivative
    // totem: number of rows needed
    //
    //
    double Lld_worst = 0.0; //stores derivative value used to determine if every parameter is near convergence
    double totem = df0.rows();//precalculates how many rows are needed
    //
    //
    if (verbose){
        end_point = system_clock::now();
        ending = time_point_cast<microseconds>(end_point).time_since_epoch().count();
        cout<<"df99,"<<(ending-start)<<",Starting"<<endl;
        gibtime = system_clock::to_time_t(system_clock::now());
        cout << ctime(&gibtime) << endl;
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
    double dint = dose_abs_max; //The amount of change used to calculate derivatives in threshold paramters
    int total_dose=0; //used later on for a section summing the dose terms
    //
    // totalnum,& Term_n,  tform, dfc,& fir,& T0,& Td0,& Tdd0,& Dose,& nonDose,& beta_0,& df0, dint, nthreads,  debugging
    cout << "starting subterms " << term_tot << endl;
    Make_Subterms( totalnum, Term_n, tform, dfc, fir, T0, Td0, Tdd0, Dose, nonDose, TTerm,  nonDose_LIN, nonDose_PLIN, nonDose_LOGLIN ,beta_0, df0,dint,nthreads, debugging);
    // ---------------------------------------------------------
    // Prints off a series of calculations to check at what point values are changing
    // ---------------------------------------------------------
    int row_check = 10;
    if (verbose){
        cout << "values checked ";
        for (int ijk=0;ijk<totalnum;ijk++){
            cout << beta_0[ijk] << " ";
        }
        cout << " " << endl;
        cout << "sums checked ";
        for (int ijk=0;ijk<totalnum;ijk++){
            cout << T0.col(ijk).sum() << " ";
        }
        cout << " " << endl;
        cout << "derivs checked ";
        for (int ijk=0;ijk<totalnum;ijk++){
            cout << Td0.col(ijk).sum() << " ";
        }
        cout << " " << endl;
        cout << "second derivs checked ";
        for (int ijk=0;ijk<totalnum;ijk++){
            cout << Tdd0.col(ijk*(ijk+1)/2+ijk).sum() << " ";
        }
        cout << " " << endl;
        cout << "dose checked ";
        for (int ijk=0;ijk<term_tot;ijk++){
            cout << Dose.col(ijk).array().sum() << " ";
        }
        cout << " " << endl;
        cout << "non-dose checked ";
        for (int ijk=0;ijk<term_tot;ijk++){
            cout << nonDose.col(ijk).array().sum() << " ";
        }
        cout << " " << endl;
        cout << "LIN_non-dose checked ";
        for (int ijk=0;ijk<term_tot;ijk++){
            cout << nonDose_LIN.col(ijk).array().sum() << " ";
        }
        cout << " " << endl;
        cout << "PLIN_non-dose checked ";
        for (int ijk=0;ijk<term_tot;ijk++){
            cout << nonDose_PLIN.col(ijk).array().sum() << " ";
        }
        cout << " " << endl;
        cout << "LOGLIN_non-dose checked ";
        for (int ijk=0;ijk<term_tot;ijk++){
            cout << nonDose_LOGLIN.col(ijk).array().sum() << " ";
        }
        cout << " " << endl;
    }
    //
    MatrixXd RdR = MatrixXd::Zero(df0.rows(), totalnum); //preallocates matrix for Risk to derivative ratios
    MatrixXd RddR = MatrixXd::Zero(df0.rows(), totalnum*(totalnum+1)/2); //preallocates matrix for Risk to second derivative ratios
    //
    if (verbose){
        end_point = system_clock::now();
        ending = time_point_cast<microseconds>(end_point).time_since_epoch().count();
        cout<<"df99,"<<(ending-start)<<",Prep_Terms"<<endl;
        gibtime = system_clock::to_time_t(system_clock::now());
        cout << ctime(&gibtime) << endl;
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
        cout<<"df100 "<<(ending-start)<<" "<<0<<" "<<0<<" "<<-1<<",Prep_R"<<endl;
        gibtime = system_clock::to_time_t(system_clock::now());
        cout << ctime(&gibtime) << endl;
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
        for (int i = 0; i < indices.size()-1; i=i+2){
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

// [[Rcpp::export]]
NumericVector peanut_risk_sub(IntegerVector Term_n, StringVector tform, NumericVector a_n,IntegerVector dfc,NumericMatrix x_all, int fir,string modelform, List Control, int term_tot){
    //----------------------------------------------------------------------------------------------------------------//
    // verbose: controls printing intermediate values and time-stamps
    // debugbing: controls printing additional time duration information
    // lr: learning rate for the second-derivative to first derivative ratio step
    // maxiter: maxmimum iterations used
    // halfmax: maximum number of half-steps used
    // epsilon: lowest acceptable maximum parameter change
    // dbeta_cap: learning rate for the first-derivative to log-likelihood ratio step
    // abs_max: starting upper limit on parameter change
    // dose_abs_max: starting upper limit on threshold parameter change
    // deriv_epsilon: threshold for lowest absolute value of derivative
    // ties_method: specifies breslow or efron tie method
    //
    double lr = Control["lr"];
    bool verbose = Control["verbose"];
    bool debugging = FALSE;
    int maxiter = Control["maxiter"];
	int halfmax = Control["halfmax"];
	double epsilon = Control["epsilon"];
	double dbeta_cap = Control["dbeta_max"];
	double abs_max = Control["abs_max"];
	double dose_abs_max = Control["dose_abs_max"];
	double deriv_epsilon =Control["deriv_epsilon"];
    NumericVector res;
    //----------------------------------------------------------------------------------------------------------------//
    // calculates risk for a reference vector
    res = RISK_SUBSET(Term_n, tform, a_n, x_all, dfc,fir,modelform,verbose, debugging, term_tot);
    return res;
}

NumericVector RISK_SUBSET(IntegerVector Term_n, StringVector tform, NumericVector a_n,NumericMatrix x_all,IntegerVector dfc,int fir,string modelform, bool verbose, bool debugging, int term_tot){
    srand (time(NULL));
    //
    // Calculates the terms and risks for a reference vector
    //
    using namespace std::chrono;
    if (verbose){
        cout << "START_NEW" << endl;
    }
    time_point<system_clock> start_point, end_point;
    start_point = system_clock::now();
    auto start = time_point_cast<microseconds>(start_point).time_since_epoch().count();
    end_point = system_clock::now();
    auto ending = time_point_cast<microseconds>(end_point).time_since_epoch().count(); //The time duration is tracked
    //
    auto gibtime = system_clock::to_time_t(system_clock::now());
    if (verbose){
        cout << ctime(&gibtime) << endl;
    }
    // df0: covariate data
    const Map<MatrixXd> df0(as<Map<MatrixXd> >(x_all));
    //
    int totalnum = Term_n.size();
    //
    cout << "Term checked ";
    for (int ij=0;ij<totalnum;ij++){
        cout << Term_n[ij] << " ";
    }
    cout << " " << endl;
    //
    //
    // cout.precision: controls the number of significant digits printed
    // nthreads: number of threads used for parallel operations
    //
    cout.precision(7); //forces higher precision numbers printed to terminal
    int nthreads = Eigen::nbThreads()-1; //stores how many threads are allocated
    //
    //
    if (verbose){
        end_point = system_clock::now();
        ending = time_point_cast<microseconds>(end_point).time_since_epoch().count();
        cout<<"df99,"<<(ending-start)<<",Starting"<<endl;
        gibtime = system_clock::to_time_t(system_clock::now());
        cout << ctime(&gibtime) << endl;
    }
    // ---------------------------------------------
    // To Start, needs to seperate the derivative terms
    // ---------------------------------------------
    //
    int total_dose=0; //used later on for a section summing the dose terms
    int dint = 1;
    //
    // totalnum,& Term_n,  tform, dfc,& fir,& T0,& Td0,& Tdd0,& Dose,& nonDose,& beta_0,& df0, dint, nthreads,  debugging
    cout << "starting subterms " << term_tot << endl;
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
    Make_Subterms( totalnum, Term_n, tform, dfc, fir, T0, Td0, Tdd0, Dose, nonDose, TTerm,  nonDose_LIN, nonDose_PLIN, nonDose_LOGLIN ,beta_0, df0,dint,nthreads, debugging);
    MatrixXd RdR = MatrixXd::Zero(df0.rows(), totalnum); //preallocates matrix for Risk to derivative ratios
    MatrixXd RddR = MatrixXd::Zero(df0.rows(), totalnum*(totalnum+1)/2);
    //
    // Calculates risk
    //
    Make_Risks(modelform, tform, Term_n, totalnum, fir, T0, Td0, Tdd0, Te, R, Rd, Rdd, Dose, nonDose, TTerm, nonDose_LIN, nonDose_PLIN, nonDose_LOGLIN, RdR, RddR, nthreads, debugging);
    if (verbose){
        cout << "values checked ";
        for (int ijk=0;ijk<totalnum;ijk++){
            cout << beta_0[ijk] << " ";
        }
        cout << " " << endl;
        cout << "sums checked ";
        for (int ijk=0;ijk<totalnum;ijk++){
            cout << T0.col(ijk).sum() << " ";
        }
        cout << " " << endl;
        cout << "derivs checked ";
        for (int ijk=0;ijk<totalnum;ijk++){
            cout << Td0.col(ijk).sum() << " ";
        }
        cout << " " << endl;
        cout << "second derivs checked ";
        for (int ijk=0;ijk<totalnum;ijk++){
            cout << Tdd0.col(ijk*(ijk+1)/2+ijk).sum() << " ";
        }
        cout << " " << endl;
        cout << "dose checked ";
        for (int ijk=0;ijk<term_tot;ijk++){
            cout << Dose.col(ijk).array().sum() << " ";
        }
        cout << " " << endl;
        cout << "non-dose checked ";
        for (int ijk=0;ijk<term_tot;ijk++){
            cout << nonDose.col(ijk).array().sum() << " ";
        }
        cout << " " << endl;
        cout << "LIN_non-dose checked ";
        for (int ijk=0;ijk<term_tot;ijk++){
            cout << nonDose_LIN.col(ijk).array().sum() << " ";
        }
        cout << " " << endl;
        cout << "PLIN_non-dose checked ";
        for (int ijk=0;ijk<term_tot;ijk++){
            cout << nonDose_PLIN.col(ijk).array().sum() << " ";
        }
        cout << " " << endl;
        cout << "LOGLIN_non-dose checked ";
        for (int ijk=0;ijk<term_tot;ijk++){
            cout << nonDose_LOGLIN.col(ijk).array().sum() << " ";
        }
        cout << " " << endl;
    }
    return wrap(R.col(0));
}


List PEANUT_PLOT_RISK(IntegerVector Term_n, StringVector tform, NumericVector a_n,NumericMatrix x_all,IntegerVector dfc,int fir, int der_iden,string modelform, double abs_max,double dose_abs_max, NumericMatrix df_groups, NumericVector tu, bool verbose, bool debugging, IntegerVector KeepConstant, int term_tot, int uniq_v){
    srand (time(NULL));
    //
    // Plots the risk over a series of covariate values
    //
    using namespace std::chrono;
    if (verbose){
        cout << "START_NEW" << endl;
    }
    time_point<system_clock> start_point, end_point;
    start_point = system_clock::now();
    auto start = time_point_cast<microseconds>(start_point).time_since_epoch().count();
    end_point = system_clock::now();
    auto ending = time_point_cast<microseconds>(end_point).time_since_epoch().count(); //The time duration is tracked
    //
    auto gibtime = system_clock::to_time_t(system_clock::now());
    if (verbose){
        cout << ctime(&gibtime) << endl;
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
    cout << "Term checked ";
    for (int ij=0;ij<totalnum;ij++){
        cout << Term_n[ij] << " ";
    }
    cout << " " << endl;
    //
    //
    // cout.precision: controls the number of significant digits printed
    // nthreads: number of threads used for parallel operations
    //
    cout.precision(7); //forces higher precision numbers printed to terminal
    int nthreads = Eigen::nbThreads()-1; //stores how many threads are allocated
    //
    //
    //
    if (verbose){
        end_point = system_clock::now();
        ending = time_point_cast<microseconds>(end_point).time_since_epoch().count();
        cout<<"df99,"<<(ending-start)<<",Starting"<<endl;
        gibtime = system_clock::to_time_t(system_clock::now());
        cout << ctime(&gibtime) << endl;
    }
    // ---------------------------------------------
    // To Start, needs to seperate the derivative terms
    // ---------------------------------------------
    //
    double dint = dose_abs_max; //The amount of change used to calculate derivatives in threshold paramters
    int total_dose=0; //used later on for a section summing the dose terms
    int ntime = tu.size();
    //
    // totalnum,& Term_n,  tform, dfc,& fir,& T0,& Td0,& Tdd0,& Dose,& nonDose,& beta_0,& df0, dint, nthreads,  debugging
    cout << "starting subterms " << term_tot << endl;
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
    for (int ij=0;ij<vv.size();ij++){
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




void removeRow(Eigen::MatrixXd& matrix, unsigned int rowToRemove){
    //
    //Used to resize with removed rows
    unsigned int numRows = matrix.rows()-1;
    unsigned int numCols = matrix.cols();

    if( rowToRemove < numRows )
        matrix.block(rowToRemove,0,numRows-rowToRemove,numCols) = matrix.block(rowToRemove+1,0,numRows-rowToRemove,numCols);

    matrix.conservativeResize(numRows,numCols);
}

void removeColumn(Eigen::MatrixXd& matrix, unsigned int colToRemove){
    //
    //Used to resize with removed columns
    unsigned int numRows = matrix.rows();
    unsigned int numCols = matrix.cols()-1;

    if( colToRemove < numCols )
        matrix.block(0,colToRemove,numRows,numCols-colToRemove) = matrix.block(0,colToRemove+1,numRows,numCols-colToRemove);

    matrix.conservativeResize(numRows,numCols);
}


//
// -------------------------------------------------------------------------------------------------------------------------------------------- //
// Bound function is not used by kept for reference later on
// -------------------------------------------------------------------------------------------------------------------------------------------- //
//



/*
List LogLik_Bounds( double q1, VectorXd beta_linT,VectorXd beta_loglinT,VectorXd beta_plinT,MatrixXd df_lin,MatrixXd df_loglin,MatrixXd df_plin, MatrixXd df_dose,int fir,string modelform,int ntime, NumericVector include_bool, double lr, int maxiter, int halfmax, double epsilon, double dbeta_cap, double abs_max, double dose_abs_max, double deriv_epsilon,List beta_loglin_slopes, List beta_loglin_tops , List beta_lin_slopes , List beta_lin_ints , List beta_quads , List beta_step_slopes , List beta_step_ints, NumericMatrix df_groups, NumericVector tu ,bool verbose){
    srand (time(NULL));
    //
    using namespace std::chrono;
    cout << "START_NEW" << endl;
    time_point<system_clock> start_point, end_point;
    start_point = system_clock::now();
    auto start = time_point_cast<microseconds>(start_point).time_since_epoch().count();
    end_point = system_clock::now();
    auto ending = time_point_cast<microseconds>(end_point).time_since_epoch().count(); //The time duration is tracked
    //
    auto gibtime = system_clock::to_time_t(system_clock::now());
    cout << ctime(&gibtime) << endl;
    //
    vector<double> cumulative_dose_num(df_dose.cols(),0);
    int dose_num_tot=0;
    int dose_term_tot=0;
    int loglin_size=0;
    int lin_size=0;
    int quad_size=0;
    int step_size=0;
    vector<int> dose_breaks(df_dose.cols(),0);
    vector<vector<double>> beta_loglin_slopes_CPP;
    vector<vector<double>> beta_loglin_tops_CPP;
    vector<vector<double>> beta_lin_slopes_CPP;
    vector<vector<double>> beta_lin_ints_CPP;
    vector<vector<double>> beta_quads_CPP;
    vector<vector<double>> beta_step_slopes_CPP;
    vector<vector<double>> beta_step_ints_CPP;
    //
    for (int ijk=0;ijk<df_dose.cols();ijk++){
        cumulative_dose_num[ijk] = dose_num_tot;
        NumericVector beta_loglin_slope = beta_loglin_slopes[ijk];
        NumericVector beta_lin_slope = beta_lin_slopes[ijk];
        NumericVector beta_quad = beta_quads[ijk];
        NumericVector beta_step_slope = beta_step_slopes[ijk];
        beta_loglin_slopes_CPP.push_back(as<vector<double> >(beta_loglin_slopes[ijk]));
        beta_loglin_tops_CPP.push_back(as<vector<double> >(beta_loglin_tops[ijk]));
        beta_lin_slopes_CPP.push_back(as<vector<double> >(beta_lin_slopes[ijk]));
        beta_lin_ints_CPP.push_back(as<vector<double> >(beta_lin_ints[ijk]));
        beta_quads_CPP.push_back(as<vector<double> >(beta_quads[ijk]));
        beta_step_slopes_CPP.push_back(as<vector<double> >(beta_step_slopes[ijk]));
        beta_step_ints_CPP.push_back(as<vector<double> >(beta_step_ints[ijk]));
        //
        if ((beta_loglin_slope.size()==1)&&(beta_loglin_slope[0]==0.0)){
            ;
        } else {
            if ((beta_loglin_slope.size()==1)&&(beta_loglin_slope[0]==1)){
                dose_num_tot += beta_loglin_slope.size();
                dose_breaks[ijk] += beta_loglin_slope.size();
            } else {
                dose_num_tot += beta_loglin_slope.size()*2;
                dose_breaks[ijk] += beta_loglin_slope.size();
            }
        }
        if ((beta_lin_slope.size()==1)&&(beta_lin_slope[0]==0.0)){
            ;
        } else {
            dose_num_tot += beta_lin_slope.size()*2;
            dose_breaks[ijk] += beta_lin_slope.size();
        }
        if ((beta_quad.size()==1)&&(beta_quad[0]==0.0)){
            ;
        } else {
            dose_num_tot += beta_quad.size();
            dose_breaks[ijk] += beta_quad.size();
        }
        if ((beta_step_slope.size()==1)&&(beta_step_slope[0]==0.0)){
            ;
        } else {
            dose_num_tot += beta_step_slope.size()*2;
            dose_breaks[ijk] += beta_step_slope.size();
        }
        dose_term_tot += dose_breaks[ijk];
        //
    }
    //
    int totalnum = dose_num_tot;
    //
    cout.precision(10); //forces higher precision numbers printed to terminal
    int nthreads = Eigen::nbThreads(); //stores how many threads are allocated
    //
    VectorXd beta_lin;
    VectorXd beta_loglin; //The vectors of parameter used
    VectorXd beta_plin;
        //
    if (include_bool[0]==1){
        beta_lin = beta_linT.tail(beta_linT.size()-1);
    }
    if (include_bool[1]==1){
        beta_loglin = beta_loglinT.tail(beta_loglinT.size()-1); //creates the used vectors
    }
    if (include_bool[2]==1){
        beta_plin = beta_plinT.tail(beta_plinT.size()-1);
    }
    //
    if (include_bool[0]==1){
        totalnum = totalnum + beta_lin.size();
    }
    if (include_bool[1]==1){
        totalnum = totalnum + beta_loglin.size(); //determines how many parameters are needed
    }
    if (include_bool[2]==1){
        totalnum = totalnum + beta_plin.size();
    }
    //
    VectorXd res(totalnum); //preallocates a vector of final parameters
    //
    double Lld_worst = 0.0;
    vector <string> tform(totalnum);
    double totem = df_loglin.rows();//precalculates how many rows
    //
    //
    end_point = system_clock::now();
    ending = time_point_cast<microseconds>(end_point).time_since_epoch().count();
    cout<<"df99,"<<(ending-start)<<",Starting"<<endl;
    gibtime = system_clock::to_time_t(system_clock::now());
    cout << ctime(&gibtime) << endl;
    // ---------------------------------------------
    // To Start, needs to seperate the derivative terms
    // ---------------------------------------------
    VectorXd beta_0(totalnum);
    MatrixXd df0 = MatrixXd::Zero(df_lin.rows(), totalnum); // stores memory for the derivative term parameters and columns
    MatrixXd T0 = MatrixXd::Zero(df0.rows(), totalnum); //preallocates matrix for Derivative column terms
    MatrixXd Td0 = MatrixXd::Zero(df0.rows(), totalnum); //preallocates matrix for Derivative column terms
    MatrixXd Tdd0 = MatrixXd::Zero(df0.rows(), totalnum); //preallocates matrix for Derivative column terms
    //
    MatrixXd Te = MatrixXd::Zero(df0.rows(), 1); //preallocates matrix for non-Derivative column terms
    MatrixXd R = MatrixXd::Zero(df0.rows(), 1); //preallocates matrix for Risks and derivatives
    MatrixXd Rd = MatrixXd::Zero(df0.rows(), totalnum); //preallocates matrix for Risks and derivatives
    MatrixXd Rdd = MatrixXd::Zero(df0.rows(), totalnum*(totalnum+1)/2);
    //
    MatrixXd De = MatrixXd::Zero(df_dose.rows(),dose_num_tot);
    MatrixXd Dde = MatrixXd::Zero(df_dose.rows(),dose_num_tot);
    MatrixXd Ddde = MatrixXd::Zero(df_dose.rows(),dose_num_tot*(dose_num_tot+1)/2);
    double dint = dose_abs_max;
    int total_dose=0;
    //
    //
    int dub_off=0;
    VectorXd Dose = VectorXd::Zero(df_dose.rows()); //Matrix of the total dose term values
    #pragma omp declare reduction (eig_plus: VectorXd: omp_out=omp_out+omp_in) initializer(omp_priv=VectorXd::Zero(omp_orig.size()))
    #pragma omp parallel for schedule(dynamic) num_threads(nthreads) reduction(eig_plus:Dose)
    for (int ij=0;ij<(totalnum-dose_num_tot+dose_term_tot);ij++){
        if (ij < dose_term_tot){
            int ind0 = ij;
            int ijk=0;
            while (ind0>dose_breaks[ijk]){
                ind0=ind0 - dose_breaks[ijk];
                ijk++;
            }
            //
            NumericVector beta_loglin_slope;
            NumericVector beta_loglin_top;
            NumericVector beta_lin_slope;
            NumericVector beta_lin_int;
            NumericVector beta_quad;
            NumericVector beta_step_slope;
            NumericVector beta_step_int;
            //
            beta_loglin_slope = beta_loglin_slopes[ijk];
            beta_loglin_top = beta_loglin_tops[ijk];
            beta_lin_slope = beta_lin_slopes[ijk];
            beta_lin_int = beta_lin_ints[ijk];
            beta_quad = beta_quads[ijk];
            beta_step_slope = beta_step_slopes[ijk];
            beta_step_int = beta_step_ints[ijk];
            //
            if ((beta_loglin_slope.size()==1)&&(beta_loglin_slope[0]==0.0)){
                ;
            } else {
                loglin_size = beta_loglin_slope.size();
            }
            if ((beta_lin_slope.size()==1)&&(beta_lin_slope[0]==0.0)){
                ;
            } else {
                lin_size = beta_lin_slope.size();
            }
            if ((beta_quad.size()==1)&&(beta_quad[0]==0.0)){
                ;
            } else {
                quad_size = beta_quad.size();
            }
            if ((beta_step_slope.size()==1)&&(beta_step_slope[0]==0.0)){
                ;
            } else {
                step_size = beta_step_slope.size();
            }
            //
            dub_off=0;
            if ((beta_loglin_slope[0]==1)&&(loglin_size==1)){
                dub_off=1;
            }
            if (ind0 < loglin_size){
                ArrayXd temp = (beta_loglin_top[ind0] * df_dose.col(ijk)).array().exp();
                ArrayXd temp1 = beta_loglin_slope[ind0] * temp;
                //
                if ((beta_loglin_slope[ind0]==1)&&(loglin_size==1)){
                    int ind1 = cumulative_dose_num[ijk]+ind0;
                    //
                    beta_0[ind1] = beta_loglin_top[ind0];
                    tform[ind1] = "loglin_top";
                    df0.col(ind1) = df_dose.col(ijk);
                    //
                    De.col(ind1) = temp1;
                    Dose = Dose.array() + temp1.array();
                    Dde.col(ind1) = temp1.array() * df_dose.col(ijk).array();
                    Ddde.col(ind1 * (ind1+1)/2) = temp1.array() * df_dose.col(ijk).array().square().array();
                    //
                } else {
                    int ind1 = cumulative_dose_num[ijk]+2*ind0;
                    //
                    beta_0[ind1] = beta_loglin_slope[ind0];
                    beta_0[ind1 + 1] = beta_loglin_top[ind0];
                    tform[ind1] = "loglin_slope";
                    tform[ind1 + 1] = "loglin_top";
                    df0.col(ind1) = df_dose.col(ijk);
                    df0.col(ind1 + 1) = df_dose.col(ijk);
                    //
                    De.col(ind1) = temp1;
                    De.col(ind1 + 1) = temp1;
                    Dose = Dose.array() + temp1.array();
                    Dde.col(ind1) = temp.array();
                    Dde.col(ind1 + 1) = temp1.array() * df_dose.col(ijk).array();
                    Ddde.col((ind1 + 1) * (ind1 + 2)/2 + ind1) = temp.array() * df_dose.col(ijk).array();
                    Ddde.col((ind1 + 1) * (ind1 + 2)/2 + ind1 + 1) = temp1.array() * df_dose.col(ijk).array().square().array();
                }
                
            } else if (ind0 < loglin_size + lin_size){
                int jk = ind0 - loglin_size;
                ArrayXd temp = (df_dose.col(ijk).array() - beta_lin_int[jk]);
                ArrayXd temp0 = (df_dose.col(ijk).array() - beta_lin_int[jk]+dint);
                ArrayXd temp1 = (df_dose.col(ijk).array() - beta_lin_int[jk]-dint);
                //
                int ind1 = cumulative_dose_num[ijk]+2*loglin_size - dub_off + 2*jk;
                //
                beta_0[ind1] = beta_lin_slope[jk];
                beta_0[(ind1 + 1)] = beta_lin_int[jk];
                tform[ind1] = "lin_slope";
                tform[(ind1 + 1)] = "lin_int";
                df0.col(ind1) = df_dose.col(ijk);
                df0.col((ind1 + 1)) = df_dose.col(ijk);
                //
                temp = (temp.array() < 0).select(0.0, temp);
                temp0 = (temp0.array() < 0).select(0.0, temp0);
                temp1 = (temp1.array() < 0).select(0.0, temp1);
                //
                De.col(ind1) = beta_lin_slope[jk] * temp.array();
                De.col((ind1 + 1)) = beta_lin_slope[jk] * temp.array();
                Dose = Dose.array() + De.col(ind1).array();
                Dde.col(ind1) = temp.array();
                Dde.col((ind1 + 1)) = beta_lin_slope[jk] * (temp1.array() - temp0.array()) / 2.0/dint;
                //
                Ddde.col((ind1 + 1) * ((ind1 + 1)+1)/2 + ind1) = (temp1.array() - temp0.array()) / 2.0/dint;
                Ddde.col((ind1 + 1) * ((ind1 + 1)+1)/2 + (ind1 + 1)) = beta_lin_slope[jk] * (temp1.array()-2*temp.array()+temp0.array()) / pow(dint,2);
            } else if (ind0 < loglin_size + lin_size + quad_size){
                int jk = ind0 - loglin_size - lin_size;
                ArrayXd temp = df_dose.col(ijk).array().square();
                int ind1 = cumulative_dose_num[ijk]+2*loglin_size - dub_off  + 2*lin_size+jk;
                //
                beta_0[ind1] = beta_quad[jk];
                tform[ind1] = "quad_slope";
                df0.col(ind1) = df_dose.col(ijk);
                //
                De.col(ind1) = beta_quad[jk] * temp.array();
                Dde.col(ind1) = temp.array();
                Dose = Dose.array() + De.col(ind1).array();
            } else {
                int jk = ind0 - loglin_size - lin_size - quad_size;
                ArrayXd temp = (df_dose.col(ijk).array() - beta_step_int[jk]);
                ArrayXd temp0 = (df_dose.col(ijk).array() - beta_step_int[jk]+dint);
                ArrayXd temp1 = (df_dose.col(ijk).array() - beta_step_int[jk]-dint);
                //
                int ind1 = cumulative_dose_num[ijk]+2*loglin_size - dub_off  + 2*lin_size + quad_size + 2*jk;
    //                int ind2 = ind1 + 1;
                //
                beta_0[ind1] = beta_step_slope[jk];
                beta_0[(ind1 + 1)] = beta_step_int[jk];
                tform[ind1] = "step_slope";
                tform[(ind1 + 1)] = "step_int";
                df0.col(ind1) = df_dose.col(ijk);
                df0.col((ind1 + 1)) = df_dose.col(ijk);
                //
                temp = (temp.array() < 0).select(0.0, MatrixXd::Zero(temp.rows(),temp.cols()).array()+1.0);
                temp0 = (temp0.array() < 0).select(0.0, MatrixXd::Zero(temp.rows(),temp0.cols()).array()+1.0);
                temp1 = (temp1.array() < 0).select(0.0, MatrixXd::Zero(temp.rows(),temp1.cols()).array()+1.0);
                //
                De.col(ind1) = beta_step_slope[jk] * temp.array();
                De.col((ind1 + 1)) = beta_step_slope[jk] * temp.array();
                Dde.col(ind1) = temp.array();
                Dde.col((ind1 + 1)) = beta_step_slope[jk] * (temp1.array() - temp0.array()) / 2.0/dint;
                //
                Ddde.col((ind1 + 1) * ((ind1 + 1)+1)/2 + ind1) = (temp1.array() - temp0.array()) / 2.0/dint;
                Ddde.col((ind1 + 1) * ((ind1 + 1)+1)/2 + (ind1 + 1)) = beta_step_slope[jk] * (temp1.array()-2*temp.array()+temp0.array()) / pow(dint,2);
                Dose = Dose.array() + De.col(ind1).array();
            }
            //
        } else {
            int ind0 = ij-dose_term_tot;
            int ind1 = ind0 + dose_num_tot;
            if (include_bool[0]==1){
                if (ind0 < beta_lin.size()){
                    // one exists and is one
                    beta_0[ind1] = beta_lin[ind0];
                    df0.col(ind1) = df_lin.col(ind0);
                    tform[ind1] = "lin";
                    //
                } else {
                    //one exists and its not one
                    ind0 = ind0 - beta_lin.size();
                    if (include_bool[1]==1){
                        if (ind0 < beta_loglin.size()){
                            //one and two exists and is two
                            beta_0[ind1] = beta_loglin[ind0];
                            df0.col(ind1) = df_loglin.col(ind0);
                            tform[ind1] = "loglin";
                            //
                        } else{
                            //one exists, two does, must be three
                            if (include_bool[2]!=1){
                                throw invalid_argument( "Are all three used? 0" );
                            }
                            ind0 = ind0 - beta_loglin.size();
                            beta_0[ind1] = beta_plin[ind0];
                            df0.col(ind1) = df_plin.col(ind0);
                            tform[ind1] = "plin";
                            //
                        }
                    } else{
                        //one exists, and two doesn't exist, must be three
                        if (include_bool[2]!=1){
                            throw invalid_argument( "Are all first and third used?" );
                        }
                        beta_0[ind1] = beta_plin[ind0];
                        df0.col(ind1) = df_plin.col(ind0);
                        tform[ind1] = "plin";
                        //
                    }
                }
            }else{
                //one doesn't exist
                if (include_bool[1]==1){
                    if (ind0 < beta_loglin.size()){
                        //one doesn't exist and two exists and is two
                        beta_0[ind1] = beta_loglin[ind0];
                        df0.col(ind1) = df_loglin.col(ind0);
                        tform[ind1] = "loglin";
                        //
                    } else{
                        //one doesn't exist, two does, must be three
                        if (include_bool[2]!=1){
                            throw invalid_argument( "Are all three used? 1" );
                        }
                        ind0 = ind0 - beta_loglin.size();
                        beta_0[ind1] = beta_plin[ind0];
                        df0.col(ind1) = df_plin.col(ind0);
                        tform[ind1] = "plin";
                        //
                    }
                } else{
                    //one doesn't exist, and two doesn't exist, must be three
                    if (include_bool[2]!=1){
                        throw invalid_argument( "Are all first and third used?" );
                    }
                    beta_0[ind1] = beta_plin[ind0];
                    df0.col(ind1) = df_plin.col(ind0);
                    tform[ind1] = "plin";
                    //
                }
            }
            T0.col(ind1) = (df0.col(ind1).array() * beta_0[ind1]).matrix();
            if (tform[ind1]=="lin") {
                Td0.col(ind1) = df0.col(ind1);
            } else if (tform[ind1]=="loglin") {
                T0.col(ind1) = T0.col(ind1).array().exp();
                Td0.col(ind1) = df0.col(ind1).array() * T0.col(ind1).array();
                Tdd0.col(ind1) = df0.col(ind1).array() * Td0.col(ind1).array();
            } else if (tform[ind1]=="plin") {
                T0.col(ind1) = 1 + T0.col(ind1).array();
                Td0.col(ind1) = df0.col(ind1);
            } else {
                cout << tform[ind1] << " is invalid" << endl;
                throw invalid_argument( "Invalid term type" );
            }
        }
    }
    //
//    const SparseMatrix df_s0 = df0.sparseView();
    //
    //
    #pragma omp parallel for schedule(dynamic) num_threads(nthreads)
    for (int ijk=0;ijk<dose_num_tot;ijk++){
        Tdd0.col(ijk) = Ddde.col(ijk*(ijk+1)/2+ijk);
        Td0.col(ijk) = Dde.col(ijk);
        T0.col(ijk) = De.col(ijk);
    }
    //
    cout << "values checked ";
    for (int ijk=0;ijk<totalnum;ijk++){
        cout << beta_0[ijk] << " ";
    }
    cout << " " << endl;
    cout << "sums checked ";
    for (int ijk=0;ijk<totalnum;ijk++){
        cout << T0.col(ijk).sum() << " ";
    }
    cout << " " << endl;
    cout << "derivs checked ";
    for (int ijk=0;ijk<totalnum;ijk++){
        cout << Td0.col(ijk).sum() << " ";
    }
    cout << " " << endl;
    //
    MatrixXd RdR = MatrixXd::Zero(df0.rows(), totalnum); //preallocates matrix for Risk to derivative ratios
    MatrixXd RddR = MatrixXd::Zero(df0.rows(), totalnum*(totalnum+1)/2);
    //
    end_point = system_clock::now();
    ending = time_point_cast<microseconds>(end_point).time_since_epoch().count();
    cout<<"df99,"<<(ending-start)<<",Prep_Terms"<<endl;
    gibtime = system_clock::to_time_t(system_clock::now());
    cout << ctime(&gibtime) << endl;
    //
    Make_Risks(modelform, dose_num_tot, totalnum, fir, T0, Td0, Tdd0, Te, R, Rd, Rdd, De, Dde, Ddde, Dose, RdR, RddR, nthreads);
    R = (R.array().isFinite()).select(R,0);
    Rd = (Rd.array().isFinite()).select(Rd,0);
    Rdd = (Rdd.array().isFinite()).select(Rdd,0);
    RdR = (RdR.array().isFinite()).select(RdR,0);
    RddR = (RddR.array().isFinite()).select(RddR,0);
    //
    cout << "risk checked ";
    for (int ijk=0;ijk<totalnum;ijk++){
        cout << R.col(0).sum() << " ";
    }
    cout << " " << endl;
    cout << "risk1 checked ";
    for (int ijk=0;ijk<totalnum;ijk++){
        cout << Rd.col(ijk).sum() << " ";
    }
    cout << " " << endl;
    cout << "risk2 checked ";
    for (int ijk=0;ijk<totalnum;ijk++){
        cout << Rdd.col(ijk*(ijk+1)/2+ijk).sum() << " ";
    }
    cout << " " << endl;
    //
    end_point = system_clock::now();
    ending = time_point_cast<microseconds>(end_point).time_since_epoch().count();
    cout<<"df100 "<<(ending-start)<<" "<<0<<" "<<0<<" "<<-1<<",Prep_R"<<endl;
    gibtime = system_clock::to_time_t(system_clock::now());
    cout << ctime(&gibtime) << endl;
    //
    // -------------------------------------------------------------------------------------------
    //
//    cout << 0 << endl;
    vector<string>  RiskGroup(ntime); //vector of strings detailing the rows
    IntegerMatrix RiskFail(ntime,2); //vector giving the event rows
    //
    const Map<MatrixXd> df_m(as<Map<MatrixXd> >(df_groups));
    #pragma omp parallel for schedule(dynamic) num_threads(nthreads)
    for (int ijk=0;ijk<ntime;ijk++){
        double t0 = tu[ijk];
        VectorXi select_ind_all = ((df_m.col(0).array() <= t0)&&(df_m.col(1).array()>=t0)).cast<int>();
        vector<int> indices_all;
        VectorXi select_ind_end = ((df_m.col(2).array() == 1)&&(df_m.col(1).array()==t0)).cast<int>();
        vector<int> indices_end;
        //
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
        vector<int> indices;
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
        RiskFail(ijk,0)=indices_end[0]-1;
        RiskFail(ijk,1)=indices_end[indices_end.size()-1]-1;
        //
        ostringstream oss;
        copy(indices.begin(), indices.end(),
            std::ostream_iterator<int>(oss, ","));
        RiskGroup[ijk] = oss.str();
    }
    //
    end_point = system_clock::now();
    ending = time_point_cast<microseconds>(end_point).time_since_epoch().count();
    cout<<"df100 "<<(ending-start)<<" "<<0<<" "<<0<<" "<<-1<<",Prep_List"<<endl;
    gibtime = system_clock::to_time_t(system_clock::now());
    cout << ctime(&gibtime) << endl;
    //
    // --------------------------
    // now a vector exists with row locations
    // --------------------------
    MatrixXd Rls1 =MatrixXd::Zero(ntime, 1);
    MatrixXd Rls2 =MatrixXd::Zero(ntime, totalnum);
    MatrixXd Rls3 =MatrixXd::Zero(ntime, totalnum*(totalnum+1)/2);
    MatrixXd Lls1 =MatrixXd::Zero(ntime, 1);
    MatrixXd Lls2 =MatrixXd::Zero(ntime, totalnum);
    MatrixXd Lls3 =MatrixXd::Zero(ntime, totalnum*(totalnum+1)/2);
    //The log-likelihood is calculated in parallel over the risk groups
    vector<double> Ll(totalnum,0.0);
    vector<double> Lld(totalnum,0.0);
    vector<double> Lldd(pow(totalnum,2),0.0);
    Calculate_Sides( RiskFail, RiskGroup, totalnum, ntime, R, Rd, Rdd, Rls1, Rls2, Rls3, Lls1, Lls2, Lls3,nthreads, debugging);
    //
    cout << "riskr checked ";
    for (int ijk=0;ijk<totalnum;ijk++){
        cout << Rls1.col(0).sum() << " ";
    }
    cout << " " << endl;
    cout << "risk1r checked ";
    for (int ijk=0;ijk<totalnum;ijk++){
        cout << Rls2.col(ijk).sum() << " ";
    }
    cout << " " << endl;
    cout << "risk2r checked ";
    for (int ijk=0;ijk<totalnum;ijk++){
        cout << Rls3.col(ijk*(ijk+1)/2+ijk).sum() << " ";
    }
    cout << " " << endl;
    //
    //
    cout << "riskl checked ";
    for (int ijk=0;ijk<totalnum;ijk++){
        cout << Lls1.col(0).sum() << " ";
    }
    cout << " " << endl;
    cout << "risk1l checked ";
    for (int ijk=0;ijk<totalnum;ijk++){
        cout << Lls2.col(ijk).sum() << " ";
    }
    cout << " " << endl;
    cout << "risk2l checked ";
    for (int ijk=0;ijk<totalnum;ijk++){
        cout << Lls3.col(ijk*(ijk+1)/2+ijk).sum() << " ";
    }
    cout << " " << endl;
    //
    Calc_LogLik( nthreads, RiskFail, RiskGroup, totalnum, ntime, R, Rd, Rdd,RdR,RddR, Rls1, Rls2, Rls3, Lls1, Lls2, Lls3, Ll, Lld, Lldd, debugging);
    //
    vector <double> Ll_comp(2,Ll[0]);
    //
    end_point = system_clock::now();
    ending = time_point_cast<microseconds>(end_point).time_since_epoch().count();
    cout<<"df100 "<<(ending-start)<<" "<<0<<" "<<0<<" "<<0<<",Calc"<<endl;
    gibtime = system_clock::to_time_t(system_clock::now());
    cout << ctime(&gibtime) << endl;
    cout << "df101 ";
    for (int ij=0;ij<totalnum;ij++){
        cout << Ll[ij] << " ";
    }
    cout << " " << endl;
    cout << "df102 ";
    for (int ij=0;ij<totalnum;ij++){
        cout << Lld[ij] << " ";
    }
    cout << " " << endl;
    cout << "df103 ";
    for (int ij=0;ij<totalnum;ij++){
        cout << Lldd[ij*totalnum+ij] << " ";
    }
    for (int ij=0;ij<totalnum;ij++){
        if (abs(Lld[ij]) > Lld_worst){
            Lld_worst = abs(Lld[ij]);
        }
    }
    cout << " " << endl;
    cout << "df104 ";
    for (int ij=0;ij<totalnum;ij++){
        cout << beta_0[ij] << " ";
    }
    cout << " " << endl;
    cout << "df105 ";
    for (int ij=0;ij<totalnum;ij++){
        cout << Lld[ij]/Lldd[ij*totalnum+ij] << " ";
    }
    cout << " " << endl;
    cout << "df106 ";
    for (int ij=0;ij<totalnum;ij++){
        cout << Ll[ij]/Lld[ij] << " ";
    }
    cout << " " << endl;
    cout << "df107 " << abs_max << " " << dose_abs_max << " " << Ll_comp[0] << " " << Ll_comp[1] << endl;
    //-------------------------------------------------------------------------------------------------------------------------
    NumericVector Lldd_vec = wrap(Lldd);
    Lldd_vec.attr("dim") = Dimension(totalnum, totalnum);
    //
    const Map<MatrixXd> Lldd_mat(as<Map<MatrixXd> >(Lldd_vec));
    //
    VectorXd Lld_mat = VectorXd::Map(Lld.data(), Lld.size());
//    const Map<VectorXd> Lld_mat(as<Map<VectorXd> >(Lld));
    end_point = system_clock::now();
    ending = time_point_cast<microseconds>(end_point).time_since_epoch().count();
    cout<<"df99 "<<(ending-start)<<",Start_h+"<<endl;
    gibtime = system_clock::to_time_t(system_clock::now());
    cout << ctime(&gibtime) << endl;
    double lstar = Ll[0] - q1 / 2;
    MatrixXd D0 = MatrixXd::Zero(Lldd_mat.rows(), Lldd_mat.cols());
    D0 << Lldd_mat;
    MatrixXd D1 = MatrixXd::Zero(Lldd_mat.rows(), Lldd_mat.cols());
    D1 << Lldd_mat;
    MatrixXd Theta = MatrixXd::Zero(totalnum, totalnum);
    MatrixXd Theta_0 = MatrixXd::Zero(totalnum, totalnum);
    Theta << beta_0.replicate(totalnum,1).array();
    Theta_0 << beta_0.replicate(totalnum,1).array();
    MatrixXd Dbeta = MatrixXd::Zero(totalnum,totalnum);
    VectorXd v_step(totalnum);
    vector<double> quad_coefs(3,0);
    #pragma omp parallel for schedule(dynamic) num_threads(nthreads)
    for (int ijk=0;ijk<totalnum;ijk++){
        cout << "__________________" << endl;
        MatrixXd dL2dOm2  = MatrixXd::Zero(totalnum,totalnum);
        VectorXd dL2dOmdBet(totalnum-1);
        VectorXd dOmdBet(totalnum-1);
        double h=0;
        //
        dL2dOm2 << Lldd_mat;
        dL2dOmdBet << Lld_mat.head(ijk), Lld_mat.tail(totalnum-1-ijk);
        removeRow(dL2dOm2,ijk);
        removeColumn(dL2dOm2,ijk);
        //
        dOmdBet = -1 * dL2dOm2.inverse().matrix() * dL2dOmdBet;
        //
        //
        if (ijk==0){
            cout << h << endl;
            h = Lld_mat.coeff(ijk,ijk) - dL2dOmdBet.transpose().matrix() * dL2dOm2.inverse().matrix() * dL2dOmdBet.matrix();
            cout << h << endl;
            h = abs(sqrt( q1 / h)) / 2;
            //
            //
            Dbeta(ijk,ijk) = h;
            cout << h << endl;
        } else {
            h = Lld_mat.coeff(ijk,ijk) - dL2dOmdBet.transpose().matrix() * dL2dOm2.inverse().matrix() * dL2dOmdBet.matrix();
            h = abs(sqrt( q1 / h)) / 2;
            //
            //
            Dbeta(ijk,ijk) = h;
        }
        if (ijk==0){
            Dbeta.block(1,ijk,totalnum-1,1) =  h * dOmdBet.col(0);
        } else if (ijk==totalnum-1){
            Dbeta.block(0,ijk,totalnum-1,1) =  h * dOmdBet.col(0);
            //
        } else {
            ;
            //
            Dbeta.block(0,ijk,ijk,1) = h * dOmdBet.block(0,0,ijk,1);
            Dbeta.block(ijk+1,ijk,totalnum-ijk-1,1) =  h * dOmdBet.block(ijk,0,totalnum-ijk-2,1);
        }
    }
    //
    MatrixXd bound_results = MatrixXd::Zero(totalnum,2);
    vector<double> beta_p(totalnum,0.0);
    vector<double> beta_c(totalnum,0.0);
    vector<double> beta_a(totalnum,0.0);
    vector<double> beta_best(totalnum,0.0);
    int iteration = 0;
    //--------------------------------------------------------------------------------------------------------------------------
    end_point = system_clock::now();
    ending = time_point_cast<microseconds>(end_point).time_since_epoch().count();
    cout<<"df99 "<<(ending-start)<<" "<<",Start_iter+"<<endl;
    gibtime = system_clock::to_time_t(system_clock::now());
    cout << ctime(&gibtime) << endl;
    //
    //
    double goal_running=0;
//    for (int ijk_ind=0;ijk_ind<totalnum;ijk_ind++){
    for (int ijk_ind=0;ijk_ind<1;ijk_ind++){
        iteration = 0;
        cout << "------------------------------------------------------------------------------------" << endl;
        cout << "---------------------------------+"<<ijk_ind<<"+---------------------------------------------" << endl;
        cout << "------------------------------------------------------------------------------------" << endl;
        goal_running=0;
        for (int i=0;i<totalnum;i++){
            if (i!=ijk_ind){
                goal_running+=pow(Lld[i],2);
            }
        }
        goal_running+=pow(Ll[0]-lstar,2);
        goal_running = sqrt(goal_running);
        fill(Ll_comp.begin(), Ll_comp.end(),goal_running);
        end_point = system_clock::now();
        ending = time_point_cast<microseconds>(end_point).time_since_epoch().count();
        cout<<"df100 "<<(ending-start)<<" "<<",Start_para "<<ijk_ind<<endl;
        gibtime = system_clock::to_time_t(system_clock::now());
        cout << ctime(&gibtime) << endl;
        cout << "df104 "<< Theta.col(ijk_ind).transpose() << endl;
        while (iteration < maxiter){
            iteration++;
            end_point = system_clock::now();
            ending = time_point_cast<microseconds>(end_point).time_since_epoch().count();
            cout<<"df100 "<<(ending-start)<<" "<<",Start_iter "<<iteration<<endl;
            gibtime = system_clock::to_time_t(system_clock::now());
            cout << ctime(&gibtime) << endl;
            VectorXd::Map(&beta_p[0], beta_p.size()) = Theta.col(ijk_ind);//wrap(beta_0);
            VectorXd::Map(&beta_c[0], beta_p.size()) = Theta.col(ijk_ind);//beta_c = wrap(beta_0);
            VectorXd::Map(&beta_a[0], beta_p.size()) = Theta.col(ijk_ind);//beta_a = wrap(beta_0);
            VectorXd::Map(&beta_best[0], beta_p.size()) = Theta.col(ijk_ind);//beta_best = wrap(beta_0);
            //
            //
            cout << "df111, " << Dbeta.col(ijk_ind).transpose() << endl;
            #pragma omp parallel for schedule(dynamic) num_threads(nthreads)
            for (int ijk=0;ijk<totalnum;ijk++){
                double lim_0 = abs((lstar - Ll[0]) / Lld[ijk]);
                double lim_1 = abs((Lld[ijk]) / Lldd[ijk+totalnum*ijk]);
                if (abs(Dbeta(ijk,ijk_ind)) > lim_0){
                    Dbeta(ijk, ijk_ind) = abs(lim_0) * sign(Dbeta(ijk, ijk_ind));
                }
                if (ijk!=ijk_ind){
                    if (abs(Dbeta(ijk,ijk_ind)) > lim_1){
                        Dbeta(ijk, ijk_ind) = abs(lim_1) * sign(Dbeta(ijk, ijk_ind));
                    }
                }
                if ((tform[ijk]=="step_int")||(tform[ijk]=="lin_int")){
                    if (abs(Dbeta(ijk, ijk_ind))>dose_abs_max){
                        Dbeta(ijk, ijk_ind) = dose_abs_max * sign(Dbeta(ijk, ijk_ind));
                    }
                }else{
                    if (abs(Dbeta(ijk, ijk_ind))>abs_max){
                        Dbeta(ijk, ijk_ind) = abs_max * sign(Dbeta(ijk, ijk_ind));
                    }
                }
//                if (abs(Dbeta(ijk, ijk_ind)) > abs(Theta_0(ijk,ijk_ind))){
//                    Dbeta(ijk, ijk_ind) = abs(Theta_0(ijk,ijk_ind)) * sign(Dbeta(ijk, ijk_ind));
//                }
            }
            cout << "df111, " << Dbeta.col(ijk_ind).transpose() << endl;
            if (true){
                throw invalid_argument( "stopping" );
            }
            //
            Dose = VectorXd::Zero(df_dose.rows());
            #pragma omp parallel for schedule(dynamic) num_threads(nthreads) reduction(eig_plus:Dose)
            for (int ijk=0;ijk<totalnum;ijk++){
                //
                if (tform[ijk]=="lin"){
                    beta_c[ijk] = beta_a[ijk] + Dbeta.coeff(ijk,ijk_ind);
                    T0.col(ijk) = T0.col(ijk).array() * (beta_c[ijk] / beta_p[ijk]);
                } else if (tform[ijk]=="plin"){
                    beta_c[ijk] = beta_a[ijk] + Dbeta.coeff(ijk,ijk_ind);
                    T0.col(ijk) = T0.col(ijk).array() * (1 + beta_c[ijk] * df0.col(ijk).array()) / (1 + beta_p[ijk] * df0.col(ijk).array());
                } else if (tform[ijk]=="loglin") {
                    beta_c[ijk] = beta_a[ijk] + Dbeta.coeff(ijk,ijk_ind);
                    T0.col(ijk) = T0.col(ijk).array() * ((beta_c[ijk] - beta_p[ijk]) * df0.col(ijk)).array().exp().array();
                    Td0.col(ijk) = Td0.col(ijk).array() * ((beta_c[ijk] - beta_p[ijk]) * df0.col(ijk)).array().exp().array();
                    Tdd0.col(ijk) = Tdd0.col(ijk).array() * ((beta_c[ijk] - beta_p[ijk]) * df0.col(ijk)).array().exp().array();
                } else if (tform[ijk]=="loglin_slope"){
                    beta_c[ijk] = beta_a[ijk] + Dbeta.coeff(ijk,ijk_ind);
                    beta_c[ijk+1] = beta_a[ijk+1] + Dbeta.coeff(ijk+1,ijk_ind);
                    double ach = beta_c[ijk]/beta_p[ijk];
                    MatrixXd bch = ((beta_c[ijk+1] - beta_p[ijk+1]) * df0.col(ijk)).array().exp().array();
                    //
                    De.col(ijk) = ach * bch.array() * De.col(ijk).array();
                    De.col(ijk+1) = ach * bch.array() * De.col(ijk+1).array();
                    Dde.col(ijk) = bch.array() * Dde.col(ijk).array();
                    Dde.col(ijk+1) = ach * bch.array() * Dde.col(ijk+1).array();
                    Ddde.col((ijk+1)*(ijk+2)/2+ijk) = bch.array() * Ddde.col((ijk+1)*(ijk+2)/2+ijk).array();
                    Ddde.col((ijk+1)*(ijk+2)/2+ijk+1) = ach * bch.array() * Ddde.col((ijk+1)*(ijk+2)/2+ijk+1).array();
                    Dose = Dose.array() + De.col(ijk).array();
                    //
                } else if (tform[ijk]=="loglin_top"){
                    if (ijk==0){
                        beta_c[ijk] = beta_a[ijk] + Dbeta.coeff(ijk,ijk_ind);
                        MatrixXd bch = ((beta_c[ijk] - beta_p[ijk]) * df0.col(ijk)).array().exp().array();
                        //
                        De.col(ijk) = De.col(ijk).array() * bch.array();
                        Dde.col(ijk) = Dde.col(ijk).array() * bch.array();
                        Ddde.col((ijk)*(ijk+1)/2+ijk) = Ddde.col((ijk)*(ijk+1)/2+ijk).array() * bch.array();
                        Dose = Dose.array() + De.col(ijk).array();
                        //
                    } else if (tform[ijk-1]!="loglin_slope"){
                        beta_c[ijk] = beta_a[ijk] + Dbeta.coeff(ijk,ijk_ind);
                        MatrixXd bch = ((beta_c[ijk] - beta_p[ijk]) * df0.col(ijk)).array().exp().array();
                        De.col(ijk) = De.col(ijk).array() * bch.array();
                        Dde.col(ijk) = Dde.col(ijk).array() * bch.array();
                        Ddde.col((ijk)*(ijk+1)/2+ijk) = Ddde.col((ijk)*(ijk+1)/2+ijk).array() * bch.array();
                        Dose = Dose.array() + De.col(ijk).array();
                    } else {
                        ;
                    }
                } else if (tform[ijk]=="lin_slope"){
                    beta_c[ijk] = beta_a[ijk] + Dbeta.coeff(ijk,ijk_ind);
                    beta_c[ijk+1] = beta_a[ijk+1] + Dbeta.coeff(ijk+1,ijk_ind);
                    //
                    ArrayXd temp = (df0.col(ijk).array() - beta_c[ijk+1]);
                    ArrayXd temp0 = (df0.col(ijk).array() - beta_c[ijk+1]-dint);
                    ArrayXd temp1 = (df0.col(ijk).array() - beta_c[ijk+1]+dint);
                    //
                    temp = (temp.array() < 0).select(0.0, temp);
                    temp0 = (temp0.array() < 0).select(0.0, temp0);
                    temp1 = (temp1.array() < 0).select(0.0, temp1);
                    //
                    De.col(ijk) = beta_c[ijk] * temp.array();
                    De.col((ijk + 1)) = beta_c[ijk] * temp.array();
                    //
                    Dde.col(ijk) = temp.array();
                    Dde.col((ijk + 1)) = beta_c[ijk] * (temp1.array()-temp0.array()) / 2.0/dint;
                    //
                    Ddde.col((ijk + 1) * ((ijk + 1)+1)/2 + ijk) = (temp1.array()-temp0.array()) / 2.0/dint;
                    Ddde.col((ijk + 1) * ((ijk + 1)+1)/2 + (ijk + 1)) = beta_c[ijk] * (temp1.array()-2*temp.array()+temp0.array()) / pow(dint,2);
                    Dose = Dose.array() + De.col(ijk).array();
                    //
                } else if (tform[ijk]=="quad_slope"){
                    beta_c[ijk] = beta_a[ijk] + Dbeta.coeff(ijk,ijk_ind);
                    De.col(ijk) = (beta_c[ijk]/beta_p[ijk]) * De.col(ijk).array();
                    Dose = Dose.array() + De.col(ijk).array();
                    //
                } else if (tform[ijk]=="step_slope"){
                    beta_c[ijk] = beta_a[ijk] + Dbeta.coeff(ijk,ijk_ind);
                    beta_c[ijk+1] = beta_a[ijk+1] + Dbeta.coeff(ijk+1,ijk_ind);
                    
                    ArrayXd temp = (df0.col(ijk).array() - beta_0[ijk+1]);
                    ArrayXd temp0 = (df0.col(ijk).array() - beta_0[ijk+1]-dint);
                    ArrayXd temp1 = (df0.col(ijk).array() - beta_0[ijk+1]+dint);
                    //
                    temp = (temp.array() < 0).select(0.0, MatrixXd::Zero(temp.rows(),temp.cols()).array()+1.0);
                    temp0 = (temp0.array() < 0).select(0.0, MatrixXd::Zero(temp.rows(),temp0.cols()).array()+1.0);
                    temp1 = (temp1.array() < 0).select(0.0, MatrixXd::Zero(temp.rows(),temp1.cols()).array()+1.0);
                    //
                    De.col(ijk) = beta_c[ijk] * temp.array();
                    De.col((ijk + 1)) = beta_c[ijk] * temp.array();
                    //
                    Dde.col(ijk) = temp.array();
                    Dde.col((ijk + 1)) = beta_c[ijk] * (temp1.array()-temp0.array()) / 2.0/dint;
                    //
                    Ddde.col((ijk + 1) * ((ijk + 1)+1)/2 + ijk) = (temp1.array()-temp0.array()) / 2.0/dint;
                    Ddde.col((ijk + 1) * ((ijk + 1)+1)/2 + (ijk + 1)) = beta_c[ijk] * (temp1.array()-2*temp.array()+temp0.array()) / pow(dint,2);
                    Dose = Dose.array() + De.col(ijk).array();
                    //
                } else {
                    ;
                }
                Theta(ijk, ijk_ind) = beta_c[ijk];
                //
            }
            //
            //
            //
            Make_Risks(modelform, dose_num_tot, totalnum, fir, T0, Td0, Tdd0, Te, R, Rd, Rdd, De, Dde, Ddde, Dose, RdR, RddR, nthreads);
            R = (R.array().isFinite()).select(R,0);
            Rd = (Rd.array().isFinite()).select(Rd,0);
            Rdd = (Rdd.array().isFinite()).select(Rdd,0);
            RdR = (RdR.array().isFinite()).select(RdR,0);
            RddR = (RddR.array().isFinite()).select(RddR,0);
            end_point = system_clock::now();
            ending = time_point_cast<microseconds>(end_point).time_since_epoch().count();
            fill(Ll.begin(), Ll.end(), 0.0);
            fill(Lld.begin(), Lld.end(), 0.0);
            fill(Lldd.begin(), Lldd.end(), 0.0);
            Calculate_Sides( RiskFail, RiskGroup, totalnum, ntime, R, Rd, Rdd, Rls1, Rls2, Rls3, Lls1, Lls2, Lls3,nthreads, debugging);
            //
            Calc_LogLik( nthreads, RiskFail, RiskGroup, totalnum, ntime, R, Rd, Rdd,RdR,RddR, Rls1, Rls2, Rls3, Lls1, Lls2, Lls3, Ll, Lld, Lldd, debugging);
            //
            goal_running=0;
            for (int i=0;i<totalnum;i++){
                if (i!=ijk_ind){
                    goal_running+=pow(Lld[i],2);
                }
            }
            goal_running+=pow(Ll[0]-lstar,2);
            goal_running = sqrt(goal_running);
//            Ll_comp[0] = goal_running;
            //
            //
            end_point = system_clock::now();
            ending = time_point_cast<microseconds>(end_point).time_since_epoch().count();
            cout<<"df100 "<<(ending-start)<<",step_calc+"<<endl;
            gibtime = system_clock::to_time_t(system_clock::now());
            cout << ctime(&gibtime) << endl;
            //
            cout << "df101 ";
            for (int ij=0;ij<totalnum;ij++){
                cout << Ll[ij] << " ";
            }
            cout << " " << endl;
            cout << "df102 ";
            for (int ij=0;ij<totalnum;ij++){
                cout << Lld[ij] << " ";
            }
            cout << " " << endl;
            cout << "df103 ";
            for (int ij=0;ij<totalnum;ij++){
                cout << Lldd[ij*totalnum+ij] << " ";
            }
            cout << " " << endl;
            cout << "df104 ";
            for (int ij=0;ij<totalnum;ij++){
                cout << beta_c[ij] << " ";
            }
            cout << " " << endl;
            cout << "df105 ";
            for (int ij=0;ij<totalnum;ij++){
                cout << Lld[ij]/Lldd[ij*totalnum+ij] << " ";
            }
            cout << " " << endl;
            cout << "df106 ";
            for (int ij=0;ij<totalnum;ij++){
                cout << Ll[ij]/Lld[ij] << " ";
            }
            cout << " " << endl;
            cout << "df107 " << abs_max << " " << Ll_comp[0] << " " << goal_running << endl;
            //
            if ((goal_running>Ll_comp[0])&(iteration>1)&FALSE){
                Dbeta.col(ijk_ind) = Dbeta.col(ijk_ind) * 0.5;
            } else {
                Ll_comp[0] = goal_running;
                Lldd_vec = wrap(Lldd);
                Lldd_vec.attr("dim") = Dimension(totalnum, totalnum);
                //
                Map<MatrixXd> Lldd_mat(as<Map<MatrixXd> >(Lldd_vec));
                Lld_mat = VectorXd::Map(Lld.data(), Lld.size());
                //
                D1 << Lldd_mat;
                //
                Lldd_mat.row(ijk_ind) = Lld_mat.matrix().row(0);
                Lld_mat[ijk_ind] = Ll[0] - lstar;
                v_step << Lldd_mat.inverse().matrix() * Lld_mat;
                //
                //
                quad_coefs[0] = Lldd_mat.col(ijk_ind).transpose().matrix() * D1.matrix() * Lldd_mat.col(ijk_ind).matrix();
                quad_coefs[1] = 2 * (v_step.transpose() * D1 * Lldd_mat.col(ijk_ind) -  1);
                quad_coefs[2] = v_step.transpose() * D1 * v_step;
                //
                double smallest_factor =1;
                if (abs(quad_coefs[0]) < abs(quad_coefs[1])){
                    if (abs(quad_coefs[0]) < abs(quad_coefs[2])){
                        smallest_factor = abs(quad_coefs[0]);
                    } else {
                        smallest_factor = abs(quad_coefs[2]);
                    }
                } else {
                    if (abs(quad_coefs[1]) < abs(quad_coefs[2])){
                        smallest_factor = abs(quad_coefs[1]);
                    } else {
                        smallest_factor = abs(quad_coefs[2]);
                    }
                }
                //
                quad_coefs[0] = quad_coefs[0] / smallest_factor;
                quad_coefs[1] = quad_coefs[1] / smallest_factor;
                quad_coefs[2] = quad_coefs[2] / smallest_factor;
                //
                double temp1 = pow(quad_coefs[1],2) - 4*quad_coefs[0]*quad_coefs[1];
                double temp2 = -quad_coefs[1]/2/quad_coefs[0];
    //            cout << "df110, " << quad_coefs[0] << ", " <<quad_coefs[1] << ", " <<quad_coefs[2] << ", " <<temp1 << ", " <<temp2 <<endl;
                vector<double> s_res(2,0);
                if (abs(quad_coefs[0])<1e-10){
                    s_res[0] = -1*quad_coefs[2] / quad_coefs[1];
                    Dbeta.col(ijk_ind) = -1*v_step.matrix() - s_res[0] * Lldd_mat.col(ijk_ind);
                } else if (temp1 > 0){
                    s_res[0] = temp2 + sqrt(temp1)/2/quad_coefs[0];
                    s_res[1] = temp2 - sqrt(temp1)/2/quad_coefs[0];
                    temp1 = (v_step + s_res[0] * Lldd_mat.col(ijk_ind)).transpose() * D0 * (v_step + s_res[0] * Lldd_mat.col(ijk_ind));
                    temp2 = (v_step + s_res[1] * Lldd_mat.col(ijk_ind)).transpose() * D0 * (v_step + s_res[1] * Lldd_mat.col(ijk_ind));
                    if (temp1<temp2){
                        Dbeta.col(ijk_ind) = -1*v_step.matrix() - s_res[0] * Lldd_mat.col(ijk_ind);
                    } else {
                        Dbeta.col(ijk_ind) = -1*v_step.matrix() - s_res[1] * Lldd_mat.col(ijk_ind);
                    }
                } else {
                    Dbeta.col(ijk_ind) = -1*v_step;
                }
            }
            cout << "df111, " << Dbeta.col(ijk_ind).transpose() << endl;
            //
            if (iteration > 5){
                if (iteration % (3)){
                    if (Dbeta.col(ijk_ind).array().abs().maxCoeff() < 1e-10){
                        iteration = maxiter;
                    }
//                    if (abs(Ll_comp[1]-Ll_comp[0])<10){
//                        abs_max = abs_max*0.1;
//                    }
                    Ll_comp[1] = Ll_comp[0];
                    if (abs_max < epsilon/10){
                        iteration = maxiter;
                    }
                }
            }
        }
        //
        #pragma omp parallel for schedule(dynamic) num_threads(nthreads)
        for (int ij=0;ij<totalnum;ij++){
            beta_0[ij] = Theta_0.coeff(ij,ijk_ind);
        }
        Dose = VectorXd::Zero(df_dose.rows());
        Update_Risk( totalnum, dose_num_tot, beta_0, df0, De, Dde, Ddde, T0, Td0, Tdd0, Dose, tform, nthreads, dint, debugging);
    for (int ijk=0;ijk<totalnum;ijk++){
        bound_results(ijk,0) = Theta.coeff(ijk,ijk);
    }
    end_point = system_clock::now();
    ending = time_point_cast<microseconds>(end_point).time_since_epoch().count();
    cout<<"df100 "<<(ending-start)<<" "<<",Change_Back"<<endl;
    gibtime = system_clock::to_time_t(system_clock::now());
    cout << ctime(&gibtime) << endl;
    //
    //-----------------------------------------------------------------------------------------------------------------------------------------
    //
    if ((modelform=="A")||(modelform=="PA")||(modelform=="PAE")){ //same process used for all of the additive type models
        Te = T0.array().block(0,dose_num_tot,T0.rows(),T0.cols()-dose_num_tot).rowwise().sum().array() + Dose.array();
        // computes intial risk and derivatives
        if (modelform=="A"){
            R << Te.array();
            Rd << Dde.array(), Td0.array().block(0,dose_num_tot,T0.rows(),T0.cols()-dose_num_tot);
            #pragma omp parallel for schedule(dynamic) num_threads(nthreads)
            for (int ijk=0;ijk<totalnum*(totalnum+1)/2;ijk++){
                int ij = 0;
                int jk = ijk;
                while (jk>ij){
                    ij++;
                    jk-=ij;
                }
                if ((ij<dose_num_tot)&&(jk<dose_num_tot)){
                    Rdd.col(ijk) = Ddde.col(ijk);
                } else if (ij==jk) {
                    Rdd.col(ijk) = Tdd0.col(ij);
                }
            }
        } else if ((modelform=="PAE")||(modelform=="PA")){
            if (fir>=dose_num_tot){
                Te = Te.array() - T0.col(fir).array();
            } else {
                Te = Te.array() - Dose.array();
            }
            if (modelform=="PAE"){
                Te = Te.array() + 1;
            }
            if (fir>=dose_num_tot){
                R << T0.col(fir).array() * Te.array();
                Rd << Td0.array() * T0.col(fir).rowwise().replicate(totalnum).array();//, Td0.col(0).array() * Te.array(), Td0.col(1).array() * Te.array();
                Rd.col(fir) = Td0.col(fir).array() * Te.array();
            } else {
                R << Dose.array() * Te.array();
                Rd << Td0.array() * Dose.rowwise().replicate(totalnum).array();
                #pragma omp parallel for schedule(dynamic) num_threads(nthreads)
                for (int ij=0;ij<dose_num_tot;ij++){
                    Rd.col(ij) = Dde.col(ij).array() * Te.array();
                }
            }
            #pragma omp parallel for schedule(dynamic) num_threads(nthreads)
            for (int ijk=0;ijk<totalnum*(totalnum+1)/2;ijk++){
                int ij = 0;
                int jk = ijk;
                while (jk>ij){
                    ij++;
                    jk-=ij;
                }
                if (ij==jk){
                    if (fir>=dose_num_tot){
                        if (ij==fir){
                            Rdd.col(ijk) = Tdd0.col(ij).array() * Te.col(0).array();
                        } else {
                            if (ij<dose_num_tot){
                                Rdd.col(ijk) = Ddde.col(ijk).array() * T0.col(fir).array();
                            } else {
                                Rdd.col(ijk) = Tdd0.col(ij).array() * T0.col(fir).array();
                            }
                        }
                    } else {
                        if (ij<dose_num_tot){
                            Rdd.col(ijk) = Ddde.col(ijk).array() * Dose.array().pow(-1).array() * R.col(0).array();
                        } else {
                            Rdd.col(ijk) = Tdd0.col(ij).array() * De.col(0).array();
                        }
                    }
                } else {
                    if (fir!=0){
                        if ((ij==fir)||(jk==fir)){
                            if (ij<dose_num_tot){
                                Rdd.col(ijk) = Dde.col(ij).array() * Td0.col(jk).array();
                            } else {
                                Rdd.col(ijk) = Td0.col(ij).array() * Td0.col(jk).array();
                            }
                        } else if ((ij<dose_num_tot)&&(jk<dose_num_tot)){
                            Rdd.col(ijk) = Ddde.col(ijk).array() * Dde.col(jk).array();
                        }
                    }
                }
            }
        }
        RdR << R.rowwise().replicate(totalnum).array().pow(-1).array() * Rd.array();
        RddR << R.rowwise().replicate(totalnum*(totalnum+1)/2).array().pow(-1).array() * Rdd.array();
    }else if (modelform=="M"){
        Te = Te.array() * 0 + 1; //verifies the initial term product is 1
        //
        Te = T0.array().block(0,dose_num_tot,T0.rows(),T0.cols()-dose_num_tot).rowwise().prod().array() * Dose.array();
        // computes intial risk and derivatives
        R << Te.array();
        Rd = T0.array().pow(-1).array() * Te.rowwise().replicate(totalnum).array();
        Rd = Rd.array() * Td0.array();
        #pragma omp parallel for schedule(dynamic) num_threads(nthreads)
        for (int ijk=0;ijk<dose_num_tot;ijk++){
            Rd.col(ijk) = Rd.col(ijk).array() * T0.array().col(ijk).array() * Dose.array().pow(-1).array();
        }
        R = (R.array().isFinite()).select(R,0);
        Rd = (Rd.array().isFinite()).select(Rd,0);
        #pragma omp parallel for schedule(dynamic) num_threads(nthreads)
        for (int ijk=0;ijk<totalnum*(totalnum+1)/2;ijk++){
            int ij = 0;
            int jk = ijk;
            while (jk>ij){
                ij++;
                jk-=ij;
            }
            if (ij==jk){
                if (ij<dose_num_tot){
                    Rdd.col(ijk) = Ddde.col(ijk).array() * Dose.array().pow(-1).array() * R.col(0).array();
                } else {
                    Rdd.col(ijk) = Tdd0.col(jk).array() * T0.col(jk).array().pow(-1).array() * R.col(0).array();
                }
            } else {
                if ((ij<dose_num_tot)&&(jk<dose_num_tot)){
                    Rdd.col(ijk) = Ddde.col(ijk).array() * Dose.array().pow(-1).array() * R.col(0).array();
                } else if ((ij<dose_num_tot)||(jk<dose_num_tot)){
                    Rdd.col(ijk) = Dde.col(jk).array() * Dose.array().pow(-1).array() * Rd.col(ij).array();
                } else{
                    Rdd.col(ijk) = Td0.col(jk).array() * T0.col(jk).array().pow(-1).array() * Rd.col(ij).array();
                }
            }
        }
        RdR << R.rowwise().replicate(totalnum).array().pow(-1).array() * Rd.array();
        RddR << R.rowwise().replicate(totalnum*(totalnum+1)/2).array().pow(-1).array() * Rdd.array();
    } else if (modelform=="GM"){
        //currently isn't implemented, it can be calculated but not optimized the same way
        throw invalid_argument( "GM isn't implemented" );
    } else {
        throw invalid_argument( "Model isn't implemented" );
    }
    R = (R.array().isFinite()).select(R,0);
    Rd = (Rd.array().isFinite()).select(Rd,0);
    Rdd = (Rdd.array().isFinite()).select(Rdd,0);
    end_point = system_clock::now();
    ending = time_point_cast<microseconds>(end_point).time_since_epoch().count();
    cout<<"df100 "<<(ending-start)<<",Update_R"<<endl;
    gibtime = system_clock::to_time_t(system_clock::now());
    cout << ctime(&gibtime) << endl;
    fill(Ll.begin(), Ll.end(), 0.0);
    fill(Lld.begin(), Lld.end(), 0.0);
    fill(Lldd.begin(), Lldd.end(), 0.0);
    Calculate_Sides( RiskFail, RiskGroup, totalnum, ntime, R, Rd, Rdd, Rls1, Rls2, Rls3, Lls1, Lls2, Lls3,nthreads, debugging);
    //
    Calc_LogLik( nthreads, RiskFail, RiskGroup, totalnum, ntime, R, Rd, Rdd,RdR,RddR, Rls1, Rls2, Rls3, Lls1, Lls2, Lls3, Ll, Lld, Lldd, debugging);
    //-------------------------------------------------------------------------------------------------------------------------
    Lldd_vec = wrap(Lldd);
    Lldd_vec.attr("dim") = Dimension(totalnum, totalnum);
    //
    Lld_mat = VectorXd::Map(Lld.data(), Lld.size());
//    const Map<VectorXd> Lld_mat(as<Map<VectorXd> >(Lld));
    end_point = system_clock::now();
    ending = time_point_cast<microseconds>(end_point).time_since_epoch().count();
    cout<<"df100 "<<(ending-start)<<" "<<",Start_h-"<<endl;
    gibtime = system_clock::to_time_t(system_clock::now());
    cout << ctime(&gibtime) << endl;
    D1 << Lldd_mat;
    Theta << Theta_0;
    Dbeta = MatrixXd::Zero(totalnum,totalnum);
    #pragma omp parallel for schedule(dynamic) num_threads(nthreads)
    for (int ijk=0;ijk<totalnum;ijk++){
        MatrixXd dL2dOm2  = MatrixXd::Zero(totalnum,totalnum);
        VectorXd dL2dOmdBet(totalnum-1);
        VectorXd dOmdBet(totalnum-1);
        double h=0;
        //
        dL2dOm2 << Lldd_mat;
        dL2dOmdBet << Lld_mat.head(ijk), Lld_mat.tail(totalnum-1-ijk);
        removeRow(dL2dOm2,ijk);
        removeColumn(dL2dOm2,ijk);
        //
        dOmdBet = -1 * dL2dOm2.inverse().matrix() * dL2dOmdBet;
        //
        //
        h = Lld_mat.coeff(ijk,ijk) - dL2dOmdBet.transpose().matrix() * dL2dOm2.inverse().matrix() * dL2dOmdBet.matrix();
        h = -1 * abs(sqrt( q1 / h)) / 2;
        //
        //
//        Theta(ijk,ijk) = Theta(ijk,ijk) + h;
        Dbeta(ijk,ijk) = h;
        
        if (ijk==0){
            Dbeta.block(1,ijk,totalnum-1,1) =  h * dOmdBet.col(0);
        } else if (ijk==totalnum-1){
            Dbeta.block(0,ijk,totalnum-1,1) =  h * dOmdBet.col(0);
            //
        } else {
            ;
            //
            Dbeta.block(0,ijk,ijk,1) = h * dOmdBet.block(0,0,ijk,1);
            Dbeta.block(ijk+1,ijk,totalnum-ijk-1,1) =  h * dOmdBet.block(ijk,0,totalnum-ijk-2,1);
        }
    }
    iteration = 0;
    end_point = system_clock::now();
    ending = time_point_cast<microseconds>(end_point).time_since_epoch().count();
    cout<<"df100 "<<(ending-start)<<" "<<",Start_iter-"<<endl;
    gibtime = system_clock::to_time_t(system_clock::now());
    cout << ctime(&gibtime) << endl;
    fill(beta_p.begin(), beta_p.end(), 0.0);
    fill(beta_c.begin(), beta_c.end(), 0.0);
    fill(beta_a.begin(), beta_a.end(), 0.0);
    fill(beta_best.begin(), beta_best.end(), 0.0);
    //--------------------------------------------------------------------------------------------------------------------------
//    for (int ijk_ind=0;ijk_ind<totalnum;ijk_ind++){
    for (int ijk_ind=0;ijk_ind<1;ijk_ind++){
        iteration = 0;
        cout << "------------------------------------------------------------------------------------" << endl;
        cout << "--------------------------------_-"<<ijk_ind<<"-_--------------------------------------------" << endl;
        cout << "------------------------------------------------------------------------------------" << endl;
        goal_running=0;
        for (int i=0;i<totalnum;i++){
            if (i!=ijk_ind){
                goal_running+=pow(Lld[i],2);
            }
        }
        goal_running+=pow(Ll[0]-lstar,2);
        goal_running = sqrt(goal_running);
        fill(Ll_comp.begin(), Ll_comp.end(),goal_running);
        end_point = system_clock::now();
        ending = time_point_cast<microseconds>(end_point).time_since_epoch().count();
        cout<<"df100 "<<(ending-start)<<" "<<",Start_para "<<ijk_ind<<endl;
        gibtime = system_clock::to_time_t(system_clock::now());
        cout << ctime(&gibtime) << endl;
        cout << "df104 "<< Theta.col(ijk_ind).transpose() << endl;
        while (iteration < maxiter){
            iteration++;
            end_point = system_clock::now();
            ending = time_point_cast<microseconds>(end_point).time_since_epoch().count();
            cout<<"df100 "<<(ending-start)<<" "<<",Start_iter "<<iteration<<endl;
            gibtime = system_clock::to_time_t(system_clock::now());
            cout << ctime(&gibtime) << endl;
            VectorXd::Map(&beta_p[0], beta_p.size()) = Theta.col(ijk_ind);//wrap(beta_0);
            VectorXd::Map(&beta_c[0], beta_p.size()) = Theta.col(ijk_ind);//beta_c = wrap(beta_0);
            VectorXd::Map(&beta_a[0], beta_p.size()) = Theta.col(ijk_ind);//beta_a = wrap(beta_0);
            VectorXd::Map(&beta_best[0], beta_p.size()) = Theta.col(ijk_ind);//beta_best = wrap(beta_0);
            //
            //
            #pragma omp parallel for schedule(dynamic) num_threads(nthreads)
            for (int ijk=0;ijk<totalnum;ijk++){
                double lim_0 = abs((lstar - Ll[0]) / Lld[ijk]);
                double lim_1 = abs((Lld[ijk]) / Lldd[ijk+totalnum*ijk]);
                if (abs(Dbeta(ijk,ijk_ind)) > lim_0){
                    Dbeta(ijk, ijk_ind) = abs(lim_0) * sign(Dbeta(ijk, ijk_ind));
                }
                if (ijk!=ijk_ind){
                    if (abs(Dbeta(ijk,ijk_ind)) > lim_1){
                        Dbeta(ijk, ijk_ind) = abs(lim_1) * sign(Dbeta(ijk, ijk_ind));
                    }
                }
                if ((tform[ijk]=="step_int")||(tform[ijk]=="lin_int")){
                    if (abs(Dbeta(ijk, ijk_ind))>dose_abs_max){
                        Dbeta(ijk, ijk_ind) = dose_abs_max * sign(Dbeta(ijk, ijk_ind));
                    }
                }else{
                    if (abs(Dbeta(ijk, ijk_ind))>abs_max){
                        Dbeta(ijk, ijk_ind) = abs_max * sign(Dbeta(ijk, ijk_ind));
                    }
                }
//                if (abs(Dbeta(ijk, ijk_ind)) > abs(Theta_0(ijk,ijk_ind))){
//                    Dbeta(ijk, ijk_ind) = abs(Theta_0(ijk,ijk_ind)) * sign(Dbeta(ijk, ijk_ind));
//                }
            }
            cout << "df111, " << Dbeta.col(ijk_ind).transpose() << endl;
            //
            Dose = VectorXd::Zero(df_dose.rows());
            #pragma omp parallel for schedule(dynamic) num_threads(nthreads) reduction(eig_plus:Dose)
            for (int ijk=0;ijk<totalnum;ijk++){
                //
                if (tform[ijk]=="lin"){
                    beta_c[ijk] = beta_a[ijk] + Dbeta.coeff(ijk,ijk_ind);
                    T0.col(ijk) = T0.col(ijk).array() * (beta_c[ijk] / beta_p[ijk]);
                } else if (tform[ijk]=="plin"){
                    beta_c[ijk] = beta_a[ijk] + Dbeta.coeff(ijk,ijk_ind);
                    T0.col(ijk) = T0.col(ijk).array() * (1 + beta_c[ijk] * df0.col(ijk).array()) / (1 + beta_p[ijk] * df0.col(ijk).array());
                } else if (tform[ijk]=="loglin") {
                    beta_c[ijk] = beta_a[ijk] + Dbeta.coeff(ijk,ijk_ind);
                    T0.col(ijk) = T0.col(ijk).array() * ((beta_c[ijk] - beta_p[ijk]) * df0.col(ijk)).array().exp().array();
                    Td0.col(ijk) = Td0.col(ijk).array() * ((beta_c[ijk] - beta_p[ijk]) * df0.col(ijk)).array().exp().array();
                    Tdd0.col(ijk) = Tdd0.col(ijk).array() * ((beta_c[ijk] - beta_p[ijk]) * df0.col(ijk)).array().exp().array();
                } else if (tform[ijk]=="loglin_slope"){
                    beta_c[ijk] = beta_a[ijk] + Dbeta.coeff(ijk,ijk_ind);
                    beta_c[ijk+1] = beta_a[ijk+1] + Dbeta.coeff(ijk+1,ijk_ind);
                    double ach = beta_c[ijk]/beta_p[ijk];
                    MatrixXd bch = ((beta_c[ijk+1] - beta_p[ijk+1]) * df0.col(ijk)).array().exp().array();
                    //
                    De.col(ijk) = ach * bch.array() * De.col(ijk).array();
                    De.col(ijk+1) = ach * bch.array() * De.col(ijk+1).array();
                    Dde.col(ijk) = bch.array() * Dde.col(ijk).array();
                    Dde.col(ijk+1) = ach * bch.array() * Dde.col(ijk+1).array();
                    Ddde.col((ijk+1)*(ijk+2)/2+ijk) = bch.array() * Ddde.col((ijk+1)*(ijk+2)/2+ijk).array();
                    Ddde.col((ijk+1)*(ijk+2)/2+ijk+1) = ach * bch.array() * Ddde.col((ijk+1)*(ijk+2)/2+ijk+1).array();
                    Dose = Dose.array() + De.col(ijk).array();
                    //
                } else if (tform[ijk]=="loglin_top"){
                    if (ijk==0){
                        beta_c[ijk] = beta_a[ijk] + Dbeta.coeff(ijk,ijk_ind);
                        MatrixXd bch = ((beta_c[ijk] - beta_p[ijk]) * df0.col(ijk)).array().exp().array();
                        //
                        De.col(ijk) = De.col(ijk).array() * bch.array();
                        Dde.col(ijk) = Dde.col(ijk).array() * bch.array();
                        Ddde.col((ijk)*(ijk+1)/2+ijk) = Ddde.col((ijk)*(ijk+1)/2+ijk).array() * bch.array();
                        Dose = Dose.array() + De.col(ijk).array();
                        //
                    } else if (tform[ijk-1]!="loglin_slope"){
                        beta_c[ijk] = beta_a[ijk] + Dbeta.coeff(ijk,ijk_ind);
                        MatrixXd bch = ((beta_c[ijk] - beta_p[ijk]) * df0.col(ijk)).array().exp().array();                        
                        De.col(ijk) = De.col(ijk).array() * bch.array();
                        Dde.col(ijk) = Dde.col(ijk).array() * bch.array();
                        Ddde.col((ijk)*(ijk+1)/2+ijk) = Ddde.col((ijk)*(ijk+1)/2+ijk).array() * bch.array();
                        Dose = Dose.array() + De.col(ijk).array();
                    } else {
                        ;
                    }
                } else if (tform[ijk]=="lin_slope"){
                    beta_c[ijk] = beta_a[ijk] + Dbeta.coeff(ijk,ijk_ind);
                    beta_c[ijk+1] = beta_a[ijk+1] + Dbeta.coeff(ijk+1,ijk_ind);
                    //
                    ArrayXd temp = (df0.col(ijk).array() - beta_c[ijk+1]);
                    ArrayXd temp0 = (df0.col(ijk).array() - beta_c[ijk+1]-dint);
                    ArrayXd temp1 = (df0.col(ijk).array() - beta_c[ijk+1]+dint);
                    //
                    temp = (temp.array() < 0).select(0.0, temp);
                    temp0 = (temp0.array() < 0).select(0.0, temp0);
                    temp1 = (temp1.array() < 0).select(0.0, temp1);
                    //
                    De.col(ijk) = beta_c[ijk] * temp.array();
                    De.col((ijk + 1)) = beta_c[ijk] * temp.array();
                    //
                    Dde.col(ijk) = temp.array();
                    Dde.col((ijk + 1)) = beta_c[ijk] * (temp1.array()-temp0.array()) / 2.0/dint;
                    //
                    Ddde.col((ijk + 1) * ((ijk + 1)+1)/2 + ijk) = (temp1.array()-temp0.array()) / 2.0/dint;
                    Ddde.col((ijk + 1) * ((ijk + 1)+1)/2 + (ijk + 1)) = beta_c[ijk] * (temp1.array()-2*temp.array()+temp0.array()) / pow(dint,2);
                    Dose = Dose.array() + De.col(ijk).array();
                    //
                } else if (tform[ijk]=="quad_slope"){
                    beta_c[ijk] = beta_a[ijk] + Dbeta.coeff(ijk,ijk_ind);
                    De.col(ijk) = (beta_c[ijk]/beta_p[ijk]) * De.col(ijk).array();
                    Dose = Dose.array() + De.col(ijk).array();
                    //
                } else if (tform[ijk]=="step_slope"){
                    beta_c[ijk] = beta_a[ijk] + Dbeta.coeff(ijk,ijk_ind);
                    beta_c[ijk+1] = beta_a[ijk+1] + Dbeta.coeff(ijk+1,ijk_ind);
                    
                    ArrayXd temp = (df0.col(ijk).array() - beta_0[ijk+1]);
                    ArrayXd temp0 = (df0.col(ijk).array() - beta_0[ijk+1]-dint);
                    ArrayXd temp1 = (df0.col(ijk).array() - beta_0[ijk+1]+dint);
                    //
                    temp = (temp.array() < 0).select(0.0, MatrixXd::Zero(temp.rows(),temp.cols()).array()+1.0);
                    temp0 = (temp0.array() < 0).select(0.0, MatrixXd::Zero(temp.rows(),temp0.cols()).array()+1.0);
                    temp1 = (temp1.array() < 0).select(0.0, MatrixXd::Zero(temp.rows(),temp1.cols()).array()+1.0);
                    //
                    De.col(ijk) = beta_c[ijk] * temp.array();
                    De.col((ijk + 1)) = beta_c[ijk] * temp.array();
                    //
                    Dde.col(ijk) = temp.array();
                    Dde.col((ijk + 1)) = beta_c[ijk] * (temp1.array()-temp0.array()) / 2.0/dint;
                    //
                    Ddde.col((ijk + 1) * ((ijk + 1)+1)/2 + ijk) = (temp1.array()-temp0.array()) / 2.0/dint;
                    Ddde.col((ijk + 1) * ((ijk + 1)+1)/2 + (ijk + 1)) = beta_c[ijk] * (temp1.array()-2*temp.array()+temp0.array()) / pow(dint,2);
                    Dose = Dose.array() + De.col(ijk).array();
                    //
                } else {
                    ;
                }
                Theta(ijk, ijk_ind) = beta_c[ijk];
                //
            }
            //
            //
            //
            if ((modelform=="A")||(modelform=="PA")||(modelform=="PAE")){ //same process used for all of the additive type models
                Te = T0.array().block(0,dose_num_tot,T0.rows(),T0.cols()-dose_num_tot).rowwise().sum().array() + Dose.array();
                // computes intial risk and derivatives
                if (modelform=="A"){
                    R << Te.array();
                    Rd << Dde.array(), Td0.array().block(0,dose_num_tot,T0.rows(),T0.cols()-dose_num_tot);
                    #pragma omp parallel for schedule(dynamic) num_threads(nthreads)
                    for (int ijk=0;ijk<totalnum*(totalnum+1)/2;ijk++){
                        int ij = 0;
                        int jk = ijk;
                        while (jk>ij){
                            ij++;
                            jk-=ij;
                        }
                        if ((ij<dose_num_tot)&&(jk<dose_num_tot)){
                            Rdd.col(ijk) = Ddde.col(ijk);
                        } else if (ij==jk) {
                            Rdd.col(ijk) = Tdd0.col(ij);
                        }
                    }
                } else if ((modelform=="PAE")||(modelform=="PA")){
                    if (fir!=0){
                        Te = Te.array() - T0.col(fir).array();
                    } else {
                        Te = Te.array() - Dose.array();
                    }
                    if (modelform=="PAE"){
                        Te = Te.array() + 1;
                    }
                    if (fir!=0){
                        R << T0.col(fir).array() * Te.array();
                        Rd << Td0.array() * T0.col(fir).array();//, Td0.col(0).array() * Te.array(), Td0.col(1).array() * Te.array();
                        Rd.col(fir) = Td0.col(fir).array() * Te.array();
                    } else {
                        R << Dose.array() * Te.array();
                        Rd << Td0.array() * De.array();
                        #pragma omp parallel for schedule(dynamic) num_threads(nthreads)
                        for (int ij=0;ij<dose_num_tot;ij++){
                            Rd.col(ij) = Dde.col(ij).array() * Te.array();
                        }
                    }
                    #pragma omp parallel for schedule(dynamic) num_threads(nthreads)
                    for (int ijk=0;ijk<totalnum*(totalnum+1)/2;ijk++){
                        int ij = 0;
                        int jk = ijk;
                        while (jk>ij){
                            ij++;
                            jk-=ij;
                        }
                        if (ij==jk){
                            if (fir!=0){
                                if (ij==fir){
                                    Rdd.col(ijk) = Tdd0.col(ij).array() * Te.col(0).array();
                                } else {
                                    if (ij<dose_num_tot){
                                        Rdd.col(ijk) = Ddde.col(ijk).array() * T0.col(fir).array();
                                    } else {
                                        Rdd.col(ijk) = Tdd0.col(ij).array() * T0.col(fir).array();
                                    }
                                }
                            } else {
                                if (ij<dose_num_tot){
                                    Rdd.col(ijk) = Ddde.col(ijk).array() * Dose.array().pow(-1).array() * R.col(0).array();
                                } else {
                                    Rdd.col(ijk) = Tdd0.col(ij).array() * De.col(0).array();
                                }
                            }
                        } else {
                            if (fir!=0){
                                if ((ij==fir)||(jk==fir)){
                                    if (ij<dose_num_tot){
                                        Rdd.col(ijk) = Dde.col(ij).array() * Td0.col(jk).array();
                                    } else {
                                        Rdd.col(ijk) = Td0.col(ij).array() * Td0.col(jk).array();
                                    }
                                } else if ((ij<dose_num_tot)&&(jk<dose_num_tot)){
                                    Rdd.col(ijk) = Ddde.col(ijk).array() * Dde.col(jk).array();
                                }
                            }
                        }
                    }
                }
                RdR << R.rowwise().replicate(totalnum).array().pow(-1).array() * Rd.array();
                RddR << R.rowwise().replicate(totalnum*(totalnum+1)/2).array().pow(-1).array() * Rdd.array();
            }else if (modelform=="M"){
                Te = Te.array() * 0 + 1; //verifies the initial term product is 1
                //
                Te = T0.array().block(0,dose_num_tot,T0.rows(),T0.cols()-dose_num_tot).rowwise().prod().array() * Dose.array();
                // computes intial risk and derivatives
                R << Te.array();
                Rd = T0.array().pow(-1).array() * Te.rowwise().replicate(totalnum).array();
                Rd = Rd.array() * Td0.array();
                #pragma omp parallel for schedule(dynamic) num_threads(nthreads)
                for (int ijk=0;ijk<dose_num_tot;ijk++){
                    Rd.col(ijk) = Rd.col(ijk).array() * T0.array().col(ijk).array() * Dose.array().pow(-1).array();
                }
                R = (R.array().isFinite()).select(R,0);
                Rd = (Rd.array().isFinite()).select(Rd,0);
                #pragma omp parallel for schedule(dynamic) num_threads(nthreads)
                for (int ijk=0;ijk<totalnum*(totalnum+1)/2;ijk++){
                    int ij = 0;
                    int jk = ijk;
                    while (jk>ij){
                        ij++;
                        jk-=ij;
                    }
                    if (ij==jk){
                        if (ij<dose_num_tot){
                            Rdd.col(ijk) = Ddde.col(ijk).array() * De.col(jk).array().pow(-1).array() * R.col(0).array();
                        } else {
                            Rdd.col(ijk) = Tdd0.col(jk).array() * Td0.col(jk).array().pow(-1).array() * Rd.col(ij).array();
                        }
                    } else {
                        if ((ij<dose_num_tot)&&(jk<dose_num_tot)){
                            Rdd.col(ijk) = Ddde.col(ijk).array() * De.col(ij).array().pow(-1).array() * R.col(0).array();
                        } else if ((ij<dose_num_tot)||(jk<dose_num_tot)){
                            Rdd.col(ijk) = Dde.col(jk).array() * Dose.array().pow(-1).array() * Rd.col(ij).array();
                        } else{
                            Rdd.col(ijk) = Td0.col(jk).array() * T0.col(jk).array().pow(-1).array() * Rd.col(ij).array();
                        }
                    }
                }
                RdR << R.rowwise().replicate(totalnum).array().pow(-1).array() * Rd.array();
                RddR << R.rowwise().replicate(totalnum*(totalnum+1)/2).array().pow(-1).array() * Rdd.array();
            } else if (modelform=="GM"){
                //currently isn't implemented, it can be calculated but not optimized the same way
                throw invalid_argument( "GM isn't implemented" );
            } else {
                throw invalid_argument( "Model isn't implemented" );
            }
            R = (R.array().isFinite()).select(R,0);
            Rd = (Rd.array().isFinite()).select(Rd,0);
            Rdd = (Rdd.array().isFinite()).select(Rdd,0);
            end_point = system_clock::now();
            ending = time_point_cast<microseconds>(end_point).time_since_epoch().count();
            fill(Ll.begin(), Ll.end(), 0.0);
            fill(Lld.begin(), Lld.end(), 0.0);
            fill(Lldd.begin(), Lldd.end(), 0.0);
            Calculate_Sides( RiskFail, RiskGroup, totalnum, ntime, R, Rd, Rdd, Rls1, Rls2, Rls3, Lls1, Lls2, Lls3,nthreads, debugging);
            //
            Calc_LogLik( nthreads, RiskFail, RiskGroup, totalnum, ntime, R, Rd, Rdd,RdR,RddR, Rls1, Rls2, Rls3, Lls1, Lls2, Lls3, Ll, Lld, Lldd, debugging);
            //
            goal_running=0;
            for (int i=0;i<totalnum;i++){
                if (i!=ijk_ind){
                    goal_running+=pow(Lld[i],2);
                }
            }
            goal_running+=pow(Ll[0]-lstar,2);
            goal_running = sqrt(goal_running);
            //
            end_point = system_clock::now();
            ending = time_point_cast<microseconds>(end_point).time_since_epoch().count();
            cout<<"df100 "<<(ending-start)<<",step_calc+"<<endl;
            gibtime = system_clock::to_time_t(system_clock::now());
            cout << ctime(&gibtime) << endl;
            //
            cout << "df101 ";
            for (int ij=0;ij<totalnum;ij++){
                cout << Ll[ij] << " ";
            }
            cout << " " << endl;
            cout << "df102 ";
            for (int ij=0;ij<totalnum;ij++){
                cout << Lld[ij] << " ";
            }
            cout << " " << endl;
            cout << "df103 ";
            for (int ij=0;ij<totalnum;ij++){
                cout << Lldd[ij*totalnum+ij] << " ";
            }
            cout << " " << endl;
            cout << "df104 ";
            for (int ij=0;ij<totalnum;ij++){
                cout << beta_c[ij] << " ";
            }
            cout << " " << endl;
            cout << "df105 ";
            for (int ij=0;ij<totalnum;ij++){
                cout << Lld[ij]/Lldd[ij*totalnum+ij] << " ";
            }
            cout << " " << endl;
            cout << "df106 ";
            for (int ij=0;ij<totalnum;ij++){
                cout << Ll[ij]/Lld[ij] << " ";
            }
            cout << " " << endl;
            cout << "df107 " << abs_max << " " << Ll_comp[0] << " " << goal_running << endl;
            //
            if ((goal_running>Ll_comp[0])&(iteration>1)&FALSE){
                Dbeta.col(ijk_ind) = Dbeta.col(ijk_ind) * 0.5;
            } else {
                Ll_comp[0] = goal_running;
                Lldd_vec = wrap(Lldd);
                Lldd_vec.attr("dim") = Dimension(totalnum, totalnum);
                //
                Map<MatrixXd> Lldd_mat(as<Map<MatrixXd> >(Lldd_vec));
                Lld_mat = VectorXd::Map(Lld.data(), Lld.size());
                //
                D1 << Lldd_mat;
                //
                Lldd_mat.row(ijk_ind) = Lld_mat.matrix().row(0);
                Lld_mat[ijk_ind] = Ll[0] - lstar;
                v_step << Lldd_mat.inverse().matrix() * Lld_mat;
                //
                //
                quad_coefs[0] = Lldd_mat.col(ijk_ind).transpose().matrix() * D1.matrix() * Lldd_mat.col(ijk_ind).matrix();
                quad_coefs[1] = 2 * (v_step.transpose() * D1 * Lldd_mat.col(ijk_ind) -  1);
                quad_coefs[2] = v_step.transpose() * D1 * v_step;
                //
                double smallest_factor =1;
                if (abs(quad_coefs[0]) < abs(quad_coefs[1])){
                    if (abs(quad_coefs[0]) < abs(quad_coefs[2])){
                        smallest_factor = abs(quad_coefs[0]);
                    } else {
                        smallest_factor = abs(quad_coefs[2]);
                    }
                } else {
                    if (abs(quad_coefs[1]) < abs(quad_coefs[2])){
                        smallest_factor = abs(quad_coefs[1]);
                    } else {
                        smallest_factor = abs(quad_coefs[2]);
                    }
                }
                //
                quad_coefs[0] = quad_coefs[0] / smallest_factor;
                quad_coefs[1] = quad_coefs[1] / smallest_factor;
                quad_coefs[2] = quad_coefs[2] / smallest_factor;
                //
                double temp1 = pow(quad_coefs[1],2) - 4*quad_coefs[0]*quad_coefs[1];
                double temp2 = -quad_coefs[1]/2/quad_coefs[0];
    //            cout << "df110, " << quad_coefs[0] << ", " <<quad_coefs[1] << ", " <<quad_coefs[2] << ", " <<temp1 << ", " <<temp2 <<endl;
                vector<double> s_res(2,0);
                if (abs(quad_coefs[0])<1e-10){
                    s_res[0] = -1*quad_coefs[2] / quad_coefs[1];
                    Dbeta.col(ijk_ind) = -1*v_step.matrix() - s_res[0] * Lldd_mat.col(ijk_ind);
                } else if (temp1 > 0){
                    s_res[0] = temp2 + sqrt(temp1)/2/quad_coefs[0];
                    s_res[1] = temp2 - sqrt(temp1)/2/quad_coefs[0];
                    temp1 = (v_step + s_res[0] * Lldd_mat.col(ijk_ind)).transpose() * D0 * (v_step + s_res[0] * Lldd_mat.col(ijk_ind));
                    temp2 = (v_step + s_res[1] * Lldd_mat.col(ijk_ind)).transpose() * D0 * (v_step + s_res[1] * Lldd_mat.col(ijk_ind));
                    if (temp1<temp2){
                        Dbeta.col(ijk_ind) = -1*v_step.matrix() - s_res[0] * Lldd_mat.col(ijk_ind);
                    } else {
                        Dbeta.col(ijk_ind) = -1*v_step.matrix() - s_res[1] * Lldd_mat.col(ijk_ind);
                    }
                } else {
                    Dbeta.col(ijk_ind) = -1*v_step;
                }
            }
            cout << "df111, " << Dbeta.col(ijk_ind).transpose() << endl;
            //
            if (iteration > 5){
                if (iteration % (3)){
                    if (Dbeta.col(ijk_ind).array().abs().maxCoeff() < 1e-10){
                        iteration = maxiter;
                    }
//                    if (abs(Ll_comp[1]-Ll_comp[0])<10){
//                        abs_max = abs_max*0.1;
//                    }
                    Ll_comp[1] = Ll_comp[0];
                    if (abs_max < epsilon/10){
                        iteration = maxiter;
                    }
                }
            }
        }
        //
        #pragma omp parallel for schedule(dynamic) num_threads(nthreads)
        for (int ij=0;ij<totalnum;ij++){
            beta_0[ij] = Theta_0.coeff(ij,ijk_ind);
        }
        Dose = VectorXd::Zero(df_dose.rows());
        Update_Risk( totalnum, dose_num_tot, beta_0, df0, De, Dde, Ddde, T0, Td0, Tdd0, Dose, tform, nthreads, dint, debugging);
    }
    for (int ijk=0;ijk<totalnum;ijk++){
        bound_results(ijk,1) = Theta.coeff(ijk,ijk);
    }
    List res_list = List::create(_["beta_0"]=wrap(beta_0),_["wald_bounds"]=wrap(bound_results));
    //
    return res_list;
}

*/

NumericMatrix Schoenfeld_PEANUT(IntegerVector Term_n, StringVector tform, NumericVector a_n,NumericMatrix x_all,IntegerVector dfc,int fir, int der_iden,string modelform, double abs_max,double dose_abs_max, NumericMatrix df_groups, NumericVector tu , bool verbose, bool debugging, IntegerVector KeepConstant, int term_tot, string ties_method){
    srand (time(NULL));
    //
    // Calculates the schoenfeld residuals
    //
    List temp_list = List::create(_["Status"]="FAILED"); //used as a dummy return value for code checking
    if (verbose){
        cout << "START_NEW" << endl;
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
        cout << ctime(&gibtime) << endl;
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
    cout << "Term checked ";
    for (int ij=0;ij<totalnum;ij++){
        cout << Term_n[ij] << " ";
    }
    cout << " " << endl;
    //
    //
    // cout.precision: controls the number of significant digits printed
    // nthreads: number of threads used for parallel operations
    //
    cout.precision(7); //forces higher precision numbers printed to terminal
    int nthreads = Eigen::nbThreads()-1; //stores how many threads are allocated
    //
    //
    //
    // Lld_worst: stores the highest magnitude log-likelihood derivative
    // totem: number of rows needed
    //
    double Lld_worst = 0.0; //stores derivative value used to determine if every parameter is near convergence
    double totem = df0.rows();//precalculates how many rows are needed
    //
    if (verbose){
        end_point = system_clock::now();
        ending = time_point_cast<microseconds>(end_point).time_since_epoch().count();
        cout<<"df99,"<<(ending-start)<<",Starting"<<endl;
        gibtime = system_clock::to_time_t(system_clock::now());
        cout << ctime(&gibtime) << endl;
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
    double dint = dose_abs_max; //The amount of change used to calculate derivatives in threshold paramters
    int total_dose=0; //used later on for a section summing the dose terms
    //
    // Calculates the subterm and term values
    cout << "starting subterms " << term_tot << endl;
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
            for (int i = 0; i < InGroup.size()-1; i=i+2){
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

List LogLik_AMFIT(MatrixXd PyrC, IntegerVector Term_n, StringVector tform, NumericVector a_n,NumericMatrix x_all,IntegerVector dfc,int fir, int der_iden,string modelform, double lr, int maxiter, int halfmax, double epsilon, double dbeta_cap, double abs_max,double dose_abs_max, double deriv_epsilon,bool change_all, bool verbose, bool debugging, IntegerVector KeepConstant, int term_tot){
    srand (time(NULL));
    //
    List temp_list = List::create(_["Status"]="FAILED"); //used as a dummy return value for code checking
    if (verbose){
        cout << "START_NEW" << endl;
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
        cout << ctime(&gibtime) << endl;
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
    cout << "Term checked ";
    for (int ij=0;ij<totalnum;ij++){
        cout << Term_n[ij] << " ";
    }
    cout << " " << endl;
    //
    //
    // cout.precision: controls the number of significant digits printed
    // nthreads: number of threads used for parallel operations
    //
    cout.precision(7); //forces higher precision numbers printed to terminal
    int nthreads = Eigen::nbThreads()-1; //stores how many threads are allocated
    //
    //
    // Lld_worst: stores the highest magnitude log-likelihood derivative
    // totem: number of rows needed
    //
    double Lld_worst = 0.0; //stores derivative value used to determine if every parameter is near convergence
    double totem = df0.rows();//precalculates how many rows are needed
    //
    if (verbose){
        end_point = system_clock::now();
        ending = time_point_cast<microseconds>(end_point).time_since_epoch().count();
        cout<<"df99,"<<(ending-start)<<",Starting"<<endl;
        gibtime = system_clock::to_time_t(system_clock::now());
        cout << ctime(&gibtime) << endl;
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
    double dint = dose_abs_max; //The amount of change used to calculate derivatives in threshold paramters
    int total_dose=0; //used later on for a section summing the dose terms
    //
    cout << "starting subterms " << term_tot << endl;
    // Calculates the subterm and term values
    Make_Subterms( totalnum, Term_n, tform, dfc, fir, T0, Td0, Tdd0, Dose, nonDose, TTerm,  nonDose_LIN, nonDose_PLIN, nonDose_LOGLIN ,beta_0, df0,dint,nthreads, debugging);
    // ---------------------------------------------------------
    // Prints off a series of calculations to check at what point values are changing
    // ---------------------------------------------------------
    int row_check = 10;
    if (verbose){
        cout << "values checked ";
        for (int ijk=0;ijk<totalnum;ijk++){
            cout << beta_0[ijk] << " ";
        }
        cout << " " << endl;
        cout << "sums checked ";
        for (int ijk=0;ijk<totalnum;ijk++){
            cout << T0.col(ijk).sum() << " ";
        }
        cout << " " << endl;
        cout << "derivs checked ";
        for (int ijk=0;ijk<totalnum;ijk++){
            cout << Td0.col(ijk).sum() << " ";
        }
        cout << " " << endl;
        cout << "second derivs checked ";
        for (int ijk=0;ijk<totalnum;ijk++){
            cout << Tdd0.col(ijk*(ijk+1)/2+ijk).sum() << " ";
        }
        cout << " " << endl;
        cout << "dose checked ";
        for (int ijk=0;ijk<term_tot;ijk++){
            cout << Dose.col(ijk).array().sum() << " ";
        }
        cout << " " << endl;
        cout << "non-dose checked ";
        for (int ijk=0;ijk<term_tot;ijk++){
            cout << nonDose.col(ijk).array().sum() << " ";
        }
        cout << " " << endl;
        cout << "LIN_non-dose checked ";
        for (int ijk=0;ijk<term_tot;ijk++){
            cout << nonDose_LIN.col(ijk).array().sum() << " ";
        }
        cout << " " << endl;
        cout << "PLIN_non-dose checked ";
        for (int ijk=0;ijk<term_tot;ijk++){
            cout << nonDose_PLIN.col(ijk).array().sum() << " ";
        }
        cout << " " << endl;
        cout << "LOGLIN_non-dose checked ";
        for (int ijk=0;ijk<term_tot;ijk++){
            cout << nonDose_LOGLIN.col(ijk).array().sum() << " ";
        }
        cout << " " << endl;
    }
    //
    MatrixXd RdR = MatrixXd::Zero(df0.rows(), totalnum); //preallocates matrix for Risk to derivative ratios
    MatrixXd RddR = MatrixXd::Zero(df0.rows(), totalnum*(totalnum+1)/2);
    //
    if (verbose){
        end_point = system_clock::now();
        ending = time_point_cast<microseconds>(end_point).time_since_epoch().count();
        cout<<"df99,"<<(ending-start)<<",Prep_Terms"<<endl;
        gibtime = system_clock::to_time_t(system_clock::now());
        cout << ctime(&gibtime) << endl;
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
                cout << ijk << " had a non-positive term" << endl;
            }
        }
        for (int ijk=0;ijk<term_tot;ijk++){
            if (nonDose_LIN.col(ijk).minCoeff()<=0){
                cout << ijk << " had a non-positive Lin term" << endl;
            }
            if (nonDose_LOGLIN.col(ijk).minCoeff()<=0){
                cout << ijk << " had a non-positive loglin term" << endl;
            }
            if (nonDose_PLIN.col(ijk).minCoeff()<=0){
                cout << ijk << " had a non-positive plin term" << endl;
            }
            if (Dose.col(ijk).minCoeff()<=0){
                cout << ijk << " had a non-positive dose term" << endl;
            }
            if (nonDose.col(ijk).minCoeff()<=0){
                cout << ijk << " had a non-positive nondose term" << endl;
            }
        }
        cout << R.sum() << endl;
        cout << "A non-positive risk was detected: " << R.minCoeff() << endl;
        return temp_list;
    }
    //
    if (verbose){
        cout << "risk checked ";
        for (int ijk=0;ijk<totalnum;ijk++){
            cout << R.col(0).sum() << " ";
        }
        cout << " " << endl;
        cout << "risk1 checked ";
        for (int ijk=0;ijk<totalnum;ijk++){
            cout << Rd.col(ijk).sum() << " ";
        }
        cout << " " << endl;
        cout << "risk2 checked ";
        for (int ijk=0;ijk<totalnum;ijk++){
            cout << Rdd.col(ijk*(ijk+1)/2+ijk).sum() << " ";
        }
        cout << " " << endl;
        //
        end_point = system_clock::now();
        ending = time_point_cast<microseconds>(end_point).time_since_epoch().count();
        cout<<"df100 "<<(ending-start)<<" "<<0<<" "<<0<<" "<<-1<<",Prep_R"<<endl;
        gibtime = system_clock::to_time_t(system_clock::now());
        cout << ctime(&gibtime) << endl;
    }
    //
    // -------------------------------------------------------------------------------------------
    //
    vector<double> Ll(totalnum,0.0);
    vector<double> Lld(totalnum,0.0);
    vector<double> Lldd(pow(totalnum,2),0.0);
    //
    //
    Amfit_LogLik( nthreads, totalnum, PyrC, R, Rd, Rdd, RdR, RddR, Ll, Lld, Lldd, debugging);
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
        cout<<"df100 "<<(ending-start)<<" "<<0<<" "<<0<<" "<<0<<",Calc"<<endl;//prints the time
        gibtime = system_clock::to_time_t(system_clock::now());
        cout << ctime(&gibtime) << endl;
        cout << "df101 ";//prints the log-likelihoods
        for (int ij=0;ij<totalnum;ij++){
            cout << Ll[ij] << " ";
        }
        cout << " " << endl;
        cout << "df102 ";//prints the first derivatives
        for (int ij=0;ij<totalnum;ij++){
            cout << Lld[ij] << " ";
        }
        cout << " " << endl;
        cout << "df103 ";//prints the second derivatives
        for (int ij=0;ij<totalnum;ij++){
            cout << Lldd[ij*totalnum+ij] << " ";
        }
        for (int ij=0;ij<totalnum;ij++){//locates highest magnitude derivative
            if (abs(Lld[ij]) > Lld_worst){
                Lld_worst = abs(Lld[ij]);
            }
        }
        cout << " " << endl;
        cout << "df104 ";//prints parameter values
        for (int ij=0;ij<totalnum;ij++){
            cout << beta_0[ij] << " ";
        }
        cout << " " << endl;
        cout << "Checking Deviance " << dev << endl;
        cout << "df105 ";
        for (int ij=0;ij<totalnum;ij++){//prints the newton step value for zero derivative
            cout << Lld[ij]/Lldd[ij*totalnum+ij] << " ";
        }
        cout << " " << endl;
        cout << "df106 ";
        for (int ij=0;ij<totalnum;ij++){//prints the newton step value for zero log-likelihood
            cout << Ll[ij]/Lld[ij] << " ";
        }
        cout << " " << endl;
        cout << "df107 " << abs_max << " " << dose_abs_max << " " << Ll_comp[0] << " " << Ll_comp[1] << endl;//prints several convergence terms
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
            cout << "Starting Halves"<<endl;//prints the final changes for validation
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
                cout << "values checked ";
                for (int ijk=0;ijk<totalnum;ijk++){
                    cout << beta_c[ijk] << " ";
                }
                cout << " " << endl;
                cout << "sums checked ";
                for (int ijk=0;ijk<totalnum;ijk++){
                    cout << T0.col(ijk).sum() << " ";
                }
                cout << " " << endl;
                cout << "derivs checked ";
                for (int ijk=0;ijk<totalnum;ijk++){
                    cout << Td0.col(ijk).sum() << " ";
                }
                cout << " " << endl;
                cout << "second derivs checked ";
                for (int ijk=0;ijk<totalnum;ijk++){
                    cout << Tdd0.col(ijk*(ijk+1)/2+ijk).sum() << " ";
                }
                cout << " " << endl;
                cout << "dose checked ";
                for (int ijk=0;ijk<term_tot;ijk++){
                    cout << Dose.col(ijk).array().sum() << " ";
                }
                cout << " " << endl;
                cout << "non-dose checked ";
                for (int ijk=0;ijk<term_tot;ijk++){
                    cout << nonDose.col(ijk).array().sum() << " ";
                }
                cout << " " << endl;
                cout << "LIN_non-dose checked ";
                for (int ijk=0;ijk<term_tot;ijk++){
                    cout << nonDose_LIN.col(ijk).array().sum() << " ";
                }
                cout << " " << endl;
                cout << "PLIN_non-dose checked ";
                for (int ijk=0;ijk<term_tot;ijk++){
                    cout << nonDose_PLIN.col(ijk).array().sum() << " ";
                }
                cout << " " << endl;
                cout << "LOGLIN_non-dose checked ";
                for (int ijk=0;ijk<term_tot;ijk++){
                    cout << nonDose_LOGLIN.col(ijk).array().sum() << " ";
                }
                cout << " " << endl;
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
                cout << "risk checked ";
                for (int ijk=0;ijk<totalnum;ijk++){
                    cout << R.col(0).sum() << " ";
                }
                cout << " " << endl;
                cout << "risk1 checked ";
                for (int ijk=0;ijk<totalnum;ijk++){
                    cout << Rd.col(ijk).sum() << " ";
                }
                cout << " " << endl;
                cout << "risk2 checked ";
                for (int ijk=0;ijk<totalnum;ijk++){
                    cout << Rdd.col(ijk*(ijk+1)/2+ijk).sum() << " ";
                }
                cout << " " << endl;
                end_point = system_clock::now();
                ending = time_point_cast<microseconds>(end_point).time_since_epoch().count();
                cout<<"df100 "<<(ending-start)<<" "<<halves<<" "<<iteration<<" "<<ind0<<",Update_R"<<endl;
                //
                gibtime = system_clock::to_time_t(system_clock::now());
                cout << ctime(&gibtime) << endl;
            }
            fill(Ll.begin(), Ll.end(), 0.0);
            fill(Lld.begin(), Lld.end(), 0.0);
            fill(Lldd.begin(), Lldd.end(), 0.0);
            Amfit_LogLik( nthreads, totalnum, PyrC, R, Rd, Rdd, RdR, RddR, Ll, Lld, Lldd, debugging);
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
                cout<<"df100 "<<(ending-start)<<" "<<halves<<" "<<iteration<<" "<<ind0<<",Update_calc"<<endl;
                gibtime = system_clock::to_time_t(system_clock::now());
                cout << ctime(&gibtime) << endl;
                cout << "df101 ";
                for (int ij=0;ij<totalnum;ij++){
                    cout << Ll[ij] << " ";
                }
                cout << " " << endl;
                cout << "df102 ";
                for (int ij=0;ij<totalnum;ij++){
                    cout << Lld[ij] << " ";
                }
                cout << " " << endl;
                cout << "df103 ";
                for (int ij=0;ij<totalnum;ij++){
                    cout << Lldd[ij*totalnum+ij] << " ";
                }
                for (int ij=0;ij<totalnum;ij++){
                    if (abs(Lld[ij]) > Lld_worst){
                        Lld_worst = abs(Lld[ij]);
                    }
                }
                cout << " " << endl;
                cout << "df104 ";
                for (int ij=0;ij<totalnum;ij++){
                    cout << beta_c[ij] << " ";
                }
                cout << " " << endl;
                cout << "Checking Deviance " << dev << endl;
                cout << "df105 ";
                for (int ij=0;ij<totalnum;ij++){
                    cout << Lld[ij]/Lldd[ij*totalnum+ij] << " ";
                }
                cout << " " << endl;
                cout << "df106 ";
                for (int ij=0;ij<totalnum;ij++){
                    cout << Ll[ij]/Lld[ij] << " ";
                }
                cout << " " << endl;
                cout << "df107 " << abs_max << " " << dose_abs_max << " " << Ll_comp[0] << " " << Ll_comp[1] << endl;
            }
            #pragma omp parallel for num_threads(nthreads)
            for (int ijk=0;ijk<totalnum;ijk++){//totalnum*(totalnum+1)/2
                beta_0[ijk] = beta_c[ijk];
            }
            beta_best = beta_c;
        }
        if (beta_best!=beta_c){//if the risk matrices aren't the optimal values, then they must be recalculated
            if (verbose){
                cout << "Changing back to best"<<endl;
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
                cout << "values checked ";
                for (int ijk=0;ijk<totalnum;ijk++){
                    cout << beta_c[ijk] << " ";
                }
                cout << " " << endl;
                //
                //
                cout << "sums checked ";
                for (int ijk=0;ijk<totalnum;ijk++){
                    cout << T0.col(ijk).sum() << " ";
                }
                cout << " " << endl;
                cout << "derivs checked ";
                for (int ijk=0;ijk<totalnum;ijk++){
                    cout << Td0.col(ijk).sum() << " ";
                }
                cout << " " << endl;
                cout << "second derivs checked ";
                for (int ijk=0;ijk<totalnum;ijk++){
                    cout << Tdd0.col(ijk*(ijk+1)/2+ijk).sum() << " ";
                }
                cout << " " << endl;
                cout << "dose checked ";
                for (int ijk=0;ijk<1;ijk++){
                    cout << Dose.array().sum() << " ";
                }
                cout << " " << endl;
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
                        cout << ijk << " had a non-positive term" << endl;
                    }
                }
                for (int ijk=0;ijk<term_tot;ijk++){
                    if (nonDose_LIN.col(ijk).minCoeff()<=0){
                        cout << ijk << " had a non-positive Lin term" << endl;
                    }
                    if (nonDose_LOGLIN.col(ijk).minCoeff()<=0){
                        cout << ijk << " had a non-positive loglin term" << endl;
                    }
                    if (nonDose_PLIN.col(ijk).minCoeff()<=0){
                        cout << ijk << " had a non-positive plin term" << endl;
                    }
                    if (Dose.col(ijk).minCoeff()<=0){
                        cout << ijk << " had a non-positive dose term" << endl;
                    }
                    if (nonDose.col(ijk).minCoeff()<=0){
                        cout << ijk << " had a non-positive nondose term" << endl;
                    }
                }
                cout << R.sum() << endl;
                cout << "A non-positive risk was detected: " << R.minCoeff() << endl;
                return temp_list;
            }
            if (verbose){
                cout << "risk checked ";
                for (int ijk=0;ijk<totalnum;ijk++){
                    cout << R.col(0).sum() << " ";
                }
                cout << " " << endl;
                cout << "risk1 checked ";
                for (int ijk=0;ijk<totalnum;ijk++){
                    cout << Rd.col(ijk).sum() << " ";
                }
                cout << " " << endl;
                cout << "risk2 checked ";
                for (int ijk=0;ijk<totalnum;ijk++){
                    cout << Rdd.col(ijk*(ijk+1)/2+ijk).sum() << " ";
                }
                cout << " " << endl;
                end_point = system_clock::now();
                ending = time_point_cast<microseconds>(end_point).time_since_epoch().count();
                cout<<"df100 "<<(ending-start)<<" "<<halves<<" "<<iteration<<" "<<ind0<<",Update_R"<<endl;
                //
                gibtime = system_clock::to_time_t(system_clock::now());
                cout << ctime(&gibtime) << endl;
            }
            fill(Ll.begin(), Ll.end(), 0.0);
            fill(Lld.begin(), Lld.end(), 0.0);
            fill(Lldd.begin(), Lldd.end(), 0.0);
            Amfit_LogLik( nthreads, totalnum, PyrC, R, Rd, Rdd, RdR, RddR, Ll, Lld, Lldd, debugging);
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
        for (int ij=0;ij<totalnum;ij++){

            if (abs(Lld[ij]) > Lld_worst){
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
            cout<<"df100 "<<(ending-start)<<" "<<halves<<" "<<iteration<<" "<<ind0<<",Recalc"<<endl;
            gibtime = system_clock::to_time_t(system_clock::now());
            cout << ctime(&gibtime) << endl;
            cout << "df101 ";
            for (int ij=0;ij<totalnum;ij++){
                cout << Ll[ij] << " ";
            }
            cout << " " << endl;
            cout << "df102 ";
            for (int ij=0;ij<totalnum;ij++){
                cout << Lld[ij] << " ";
            }
            cout << " " << endl;
            cout << "df103 ";
            for (int ij=0;ij<totalnum;ij++){
                cout << Lldd[ij*totalnum+ij] << " ";
            }
            cout << " " << endl;
            cout << "df104 ";
            for (int ij=0;ij<totalnum;ij++){
                cout << beta_c[ij] << " ";
            }
            cout << " " << endl;
            cout << "Checking Deviance " << dev << endl;
            cout << "Finshed iteration" << endl;
            cout << "df105 ";
            for (int ij=0;ij<totalnum;ij++){
                cout << Lld[ij]/Lldd[ij*totalnum+ij] << " ";
            }
            cout << " " << endl;
            cout << "df106 ";
            for (int ij=0;ij<totalnum;ij++){
                cout << Ll[ij]/Lld[ij] << " ";
            }
            cout << " " << endl;
            cout << "df107 " << abs_max << " " << dose_abs_max << " " << Ll_comp[0] << " " << Ll_comp[1] << endl;
        }
    }
    // -----------------------------------------------
    // Performing Full Calculation to get full second derivative matrix
    // -----------------------------------------------
    fill(Ll.begin(), Ll.end(), 0.0);
    fill(Lld.begin(), Lld.end(), 0.0);
    fill(Lldd.begin(), Lldd.end(), 0.0);
    Amfit_LogLik( nthreads, totalnum, PyrC, R, Rd, Rdd, RdR, RddR, Ll, Lld, Lldd, debugging);
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
        cout<<"df100 "<<(ending-start)<<" "<<halves<<" "<<iteration<<" "<<ind0<<",Recalc"<<endl;
        gibtime = system_clock::to_time_t(system_clock::now());
        cout << ctime(&gibtime) << endl;
        cout << "df101 ";
        for (int ij=0;ij<totalnum;ij++){
            cout << Ll[ij] << " ";
        }
        cout << " " << endl;
        cout << "df102 ";
        for (int ij=0;ij<totalnum;ij++){
            cout << Lld[ij] << " ";
        }
        cout << " " << endl;
        cout << "df103 ";
        for (int ij=0;ij<totalnum;ij++){
            cout << Lldd[ij*totalnum+ij] << " ";
        }
        cout << " " << endl;
        cout << "df104 ";
        for (int ij=0;ij<totalnum;ij++){
            cout << beta_c[ij] << " ";
        }
        cout << " " << endl;
        cout << "Checking Deviance " << dev << endl;
        cout << "Finshed iteration" << endl;
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



void Make_Subterms(const int& totalnum, const IntegerVector& Term_n,const StringVector&  tform, const IntegerVector& dfc, const int& fir, MatrixXd& T0, MatrixXd& Td0, MatrixXd& Tdd0, MatrixXd& Dose, MatrixXd& nonDose,  MatrixXd& TTerm, MatrixXd& nonDose_LIN, MatrixXd& nonDose_PLIN, MatrixXd& nonDose_LOGLIN,const  VectorXd& beta_0,const  MatrixXd& df0,const double& dint, const int& nthreads, bool debugging){
    //
    //Make_Subterms( totalnum, dose_num_tot, dose_term_tot, dose_breaks, beta_loglin_slopes_CPP, beta_loglin_tops_CPP, beta_lin_slopes_CPP, beta_lin_ints_CPP, beta_quads_CPP, beta_step_slopes_CPP, beta_step_ints_CPP, beta_lin, beta_loglin, beta_plin, df_lin, df_loglin, df_plin, df_dose, De, Dde, Ddde, T0, Td0, Tdd0, Dose,cumulative_dose_num,beta_0, df0,dint,nthreads, tform,include_bool, debugging);
    //
    // Calculates the sub term values
    //
    #pragma omp declare reduction (eig_plus: MatrixXd: omp_out=omp_out.array() + omp_in.array()) initializer(omp_priv=MatrixXd::Constant(omp_orig.rows(),omp_orig.cols(),0.0))
    #pragma omp declare reduction (eig_mult: MatrixXd: omp_out=omp_out.array() * omp_in.array()) initializer(omp_priv=MatrixXd::Constant(omp_orig.rows(),omp_orig.cols(),1.0))
    #pragma omp parallel for schedule(dynamic) num_threads(nthreads) reduction(eig_plus:Dose,nonDose_LIN,nonDose_PLIN) reduction(eig_mult:nonDose_LOGLIN)
    for (int ij=0;ij<(totalnum);ij++){
        int df0_c = dfc[ij]-1;
        int tn = Term_n[ij];
        if (as< string>(tform[ij])=="loglin_slope"){
            ArrayXd temp = (beta_0[ij+1] * df0.col(df0_c)).array().exp();
            ArrayXd temp1 = beta_0[ij] * temp;
            //
            //
            T0.col(ij) = temp1;
            T0.col(ij+1) = temp1;
            Dose.col(tn) = Dose.col(tn).array() + T0.col(ij).array();
            
        } else if (as< string>(tform[ij])=="loglin_top"){
            if (ij==0){
                ArrayXd temp = (beta_0[ij] * df0.col(df0_c)).array().exp();
                T0.col(ij) = temp;
                Dose.col(tn) = Dose.col(tn).array() + T0.col(ij).array();

            } else if (tform[ij-1]!="loglin_slope"){
                ArrayXd temp = (beta_0[ij] * df0.col(df0_c)).array().exp();
                T0.col(ij) = temp;
                Dose.col(tn) = Dose.col(tn).array() + T0.col(ij).array();
                //

            } else {
                ;
            }
        } else if (as< string>(tform[ij])=="lin_slope"){
            ArrayXd temp = (df0.col(df0_c).array() - beta_0[ij+1]);
            ArrayXd temp0 = (df0.col(df0_c).array() - beta_0[ij+1]+dint);
            ArrayXd temp1 = (df0.col(df0_c).array() - beta_0[ij+1]-dint);
            //
            temp = (temp.array() < 0).select(0, temp);
            temp0 = (temp0.array() < 0).select(0, temp0);
            temp1 = (temp1.array() < 0).select(0, temp1);
            //
            T0.col(ij) = beta_0[ij] * temp.array();
            T0.col(ij+1) = beta_0[ij] * temp.array();
            //
            Dose.col(tn) = Dose.col(tn).array() + T0.col(ij).array();

        } else if (as< string>(tform[ij])=="quad_slope"){
            ArrayXd temp = df0.col(df0_c).array().square();
            //
            T0.col(ij) = beta_0[ij] * temp.array();
            Dose.col(tn) = Dose.col(tn).array() + T0.col(ij).array();
        } else if (as< string>(tform[ij])=="step_slope"){
            ArrayXd temp = (df0.col(df0_c).array() - beta_0[ij+1]);
            ArrayXd temp0 = (df0.col(df0_c).array() - beta_0[ij+1]+dint);
            ArrayXd temp1 = (df0.col(df0_c).array() - beta_0[ij+1]-dint);
            //
            temp = (temp.array() < 0).select(0.0, MatrixXd::Zero(temp.rows(),temp.cols()).array()+1.0);
            temp0 = (temp0.array() < 0).select(0.0, MatrixXd::Zero(temp.rows(),temp0.cols()).array()+1.0);
            temp1 = (temp1.array() < 0).select(0.0, MatrixXd::Zero(temp.rows(),temp1.cols()).array()+1.0);
            //
            T0.col(ij) = beta_0[ij] * temp.array();
            T0.col(ij+1) = beta_0[ij] * temp.array();
            Dose.col(tn) = Dose.col(tn).array() + T0.col(ij).array();

        } else if (as< string>(tform[ij])=="lin") {
            T0.col(ij) = (df0.col(df0_c).array() * beta_0[ij]).matrix();
            nonDose_LIN.col(tn) = nonDose_LIN.col(tn).array() + T0.col(ij).array();

        } else if (as< string>(tform[ij])=="loglin") {
            T0.col(ij) = (df0.col(df0_c).array() * beta_0[ij]).matrix();
            T0.col(ij) = T0.col(ij).array().exp();;
            nonDose_LOGLIN.col(tn) = nonDose_LOGLIN.col(tn).array() * T0.col(ij).array();

        } else if (as< string>(tform[ij])=="plin") {
            T0.col(ij) = (df0.col(df0_c).array() * beta_0[ij]).matrix();
            T0.col(ij) = 1 + T0.col(ij).array();
            nonDose_PLIN.col(tn) = nonDose_PLIN.col(tn).array() + T0.col(ij).array();

        } else {
            ;
        }
    }
    //
    // Calculates the terms and derivatives
    //
    Dose = (Dose.array() != 0).select(Dose,1.0);
    nonDose_LIN = (nonDose_LIN.array() != 0).select(nonDose_LIN,1.0);//replaces missing data with 1
    nonDose_LOGLIN = (nonDose_LOGLIN.array() != 0).select(nonDose_LOGLIN,1.0);
    nonDose_PLIN = (nonDose_PLIN.array() != 0).select(nonDose_PLIN,1.0);
    for (int ijk=0; ijk<nonDose.cols();ijk++){ //combines non-dose terms into a single term
        nonDose.col(ijk) = nonDose_LIN.col(ijk).array()  * nonDose_PLIN.col(ijk).array()  * nonDose_LOGLIN.col(ijk).array() ;
    }
    TTerm << Dose.array() * nonDose.array();
    #pragma omp parallel for schedule(dynamic) num_threads(nthreads)
    for (int ij=0;ij<(totalnum);ij++){
        int df0_c = dfc[ij]-1;
        int tn = Term_n[ij];
        if (as< string>(tform[ij])=="loglin_slope"){
            ArrayXd temp = (beta_0[ij+1] * df0.col(df0_c)).array().exp();
            ArrayXd temp1 = beta_0[ij] * temp;
            //
            //
            T0.col(ij) = Dose.col(tn);
            T0.col(ij+1) = Dose.col(tn);
            Td0.col(ij) = temp.array();
            Td0.col(ij+1) = temp1.array() * df0.col(df0_c).array();
            Tdd0.col((ij+1)*(ij+2)/2+ij) = temp.array() * df0.col(df0_c).array();
            Tdd0.col((ij+1)*(ij+2)/2+ij+1) = temp1.array() * df0.col(df0_c).array().square().array();
            
        } else if (as< string>(tform[ij])=="loglin_top"){
            if (ij==0){
                ArrayXd temp = (beta_0[ij] * df0.col(df0_c)).array().exp();
                T0.col(ij) = Dose.col(tn);
                Td0.col(ij) = temp.array() * df0.col(df0_c).array();
                Tdd0.col(ij * (ij+1)/2 + ij) = temp.array() * df0.col(df0_c).array().square().array();

            } else if (tform[ij-1]!="loglin_slope"){
                ArrayXd temp = (beta_0[ij] * df0.col(df0_c)).array().exp();
                T0.col(ij) = Dose.col(tn);
                Td0.col(ij) = temp.array() * df0.col(df0_c).array();
                Tdd0.col(ij * (ij+1)/2 + ij) = temp.array() * df0.col(df0_c).array().square().array();
                //

            } else {
                ;
            }
        } else if (as< string>(tform[ij])=="lin_slope"){
            ArrayXd temp = (df0.col(df0_c).array() - beta_0[ij+1]);
            ArrayXd temp0 = (df0.col(df0_c).array() - beta_0[ij+1]+dint);
            ArrayXd temp1 = (df0.col(df0_c).array() - beta_0[ij+1]-dint);
            //
            temp = (temp.array() < 0).select(0, temp);
            temp0 = (temp0.array() < 0).select(0, temp0);
            temp1 = (temp1.array() < 0).select(0, temp1);
            //
            T0.col(ij) = Dose.col(tn);
            T0.col(ij+1) = Dose.col(tn);
            Td0.col(ij) = temp.array();
            Td0.col(ij+1) = beta_0[ij] * (temp1.array()-temp0.array()) / 2/dint;
            //
            Tdd0.col((ij+1)*(ij+2)/2+ij) = (temp1.array()-temp0.array()) / 2/dint;
            Tdd0.col((ij+1)*(ij+2)/2+ij+1) = beta_0[ij] * (temp1.array()-2*temp.array()+temp0.array()) / pow(dint,2);

        } else if (as< string>(tform[ij])=="quad_slope"){
            ArrayXd temp = df0.col(df0_c).array().square();
            //
            T0.col(ij) = Dose.col(tn);
            Td0.col(ij) = temp.array();
        } else if (as< string>(tform[ij])=="step_slope"){
            ArrayXd temp = (df0.col(df0_c).array() - beta_0[ij+1]);
            ArrayXd temp0 = (df0.col(df0_c).array() - beta_0[ij+1]+dint);
            ArrayXd temp1 = (df0.col(df0_c).array() - beta_0[ij+1]-dint);
            //
            temp = (temp.array() < 0).select(0.0, MatrixXd::Zero(temp.rows(),temp.cols()).array()+1.0);
            temp0 = (temp0.array() < 0).select(0.0, MatrixXd::Zero(temp.rows(),temp0.cols()).array()+1.0);
            temp1 = (temp1.array() < 0).select(0.0, MatrixXd::Zero(temp.rows(),temp1.cols()).array()+1.0);
            //
            T0.col(ij) = Dose.col(tn);
            T0.col(ij+1) = Dose.col(tn);
            Td0.col(ij) = temp.array();
            Td0.col(ij+1) = beta_0[ij] * (temp1.array()-temp0.array()) / 2/dint;
            //
            Tdd0.col((ij+1)*(ij+2)/2+ij) = (temp1.array()-temp0.array()) / 2/dint;
            Tdd0.col((ij+1)*(ij+2)/2+ij+1) = beta_0[ij] * (temp1.array()-2*temp.array()+temp0.array()) / pow(dint,2);

        } else if (as< string>(tform[ij])=="lin") {
            T0.col(ij) = nonDose_LIN.col(tn);
            Td0.col(ij) = df0.col(df0_c);

        } else if (as< string>(tform[ij])=="loglin") {
            T0.col(ij) = nonDose_LOGLIN.col(tn);
            Td0.col(ij) = df0.col(df0_c).array() * T0.col(ij).array();
            Tdd0.col((ij)*(ij+1)/2+ij) = df0.col(df0_c).array() * Td0.col(ij).array();
        } else if (as< string>(tform[ij])=="plin") {
            T0.col(ij) = nonDose_PLIN.col(tn);
            Td0.col(ij) = df0.col(df0_c);
        } else {
            ;
        }
    }
    //
    // Adds in possible log-linear subterm second derivatives between DIFFERENT covariates
    //
    #pragma omp parallel for schedule(dynamic) num_threads(nthreads)
    for (int ijk=0;ijk<totalnum*(totalnum+1)/2;ijk++){
        int ij = 0;
        int jk = ijk;
        while (jk>ij){
            ij++;
            jk-=ij;
        }
        int tij = Term_n[ij];
        int tjk = Term_n[jk];
        int df0_ij = dfc[ij]-1;
        int df0_jk = dfc[jk]-1;
        if (tij==tjk){
            if (as< string>(tform[ij])=="loglin") {
                if (ij==jk){
                    T0.col(ij) = nonDose_LOGLIN.col(tij);
                    Td0.col(ij) = df0.col(df0_ij).array() * nonDose_LOGLIN.col(tij).array();
                    Tdd0.col((ij)*(ij+1)/2+ij) = df0.col(df0_ij).array().pow(2).array() * nonDose_LOGLIN.col(tij).array();
                } else if (as< string>(tform[jk])=="loglin") {
                    Tdd0.col((ij)*(ij+1)/2+jk) = df0.col(df0_ij).array() * df0.col(df0_jk).array() * nonDose_LOGLIN.col(tij).array();
                }
            }
        }
    }
    return;
}


void Make_Risks(string modelform, const StringVector& tform, const IntegerVector& Term_n, const int& totalnum, const int& fir, const MatrixXd& T0, const MatrixXd& Td0, const MatrixXd& Tdd0, MatrixXd& Te, MatrixXd& R, MatrixXd& Rd, MatrixXd& Rdd, MatrixXd& Dose, MatrixXd& nonDose,  MatrixXd& TTerm,  MatrixXd& nonDose_LIN, MatrixXd& nonDose_PLIN, MatrixXd& nonDose_LOGLIN, MatrixXd& RdR, MatrixXd& RddR, const int& nthreads, bool debugging){
    //
    //Make_Risks(modelform, tform, Term_n, totalnum, fir, T0, Td0, Tdd0, Te, R, Rd, Rdd, Dose, nonDose, RdR, RddR, nthreads, debugging);
    //
    set<string> Dose_Iden; //List of dose subterms
    Dose_Iden.insert("loglin_top");
    Dose_Iden.insert("loglin_slope");
    Dose_Iden.insert("lin_slope");
    Dose_Iden.insert( "lin_int");
    Dose_Iden.insert("quad_slope");
    Dose_Iden.insert("step_slope");
    Dose_Iden.insert("step_int");
    //
    if ((modelform=="A")||(modelform=="PA")||(modelform=="PAE")){ //same process used for all of the additive type models
        Te = TTerm.array().rowwise().sum().array();
        // computes intial risk and derivatives
        if (modelform=="A"){
            R << Te.array();
            #pragma omp parallel for schedule(dynamic) num_threads(nthreads)
            for (int ijk=0;ijk<totalnum*(totalnum+1)/2;ijk++){
                int ij = 0;
                int jk = ijk;
                while (jk>ij){
                    ij++;
                    jk-=ij;
                }
                int tij = Term_n[ij];
                int tjk = Term_n[jk];
                if (ij==jk){
                    if (Dose_Iden.find(as< string>(tform[ij])) != Dose_Iden.end()){
                        Rd.col(ij) =   TTerm.col(tij).array() * Dose.col(tij).array().pow(-1).array() *   Td0.col(ij).array();
                        Rdd.col(ijk) = TTerm.col(tij).array() * Dose.col(tij).array().pow(-1).array() *   Tdd0.col(ijk).array();
                    } else if (tform[ij]=="lin") {
                        Rd.col(ij) =   TTerm.col(tij).array() * nonDose_LIN.col(tij).array().pow(-1).array() *   Td0.col(ij).array();
                    } else if (tform[ij]=="plin") {
                        Rd.col(ij) =   TTerm.col(tij).array() * nonDose_PLIN.col(tij).array().pow(-1).array() *   Td0.col(ij).array();
                    } else if (tform[ij]=="loglin") {
                        Rd.col(ij) =   TTerm.col(tij).array() * nonDose_LOGLIN.col(tij).array().pow(-1).array() * Td0.col(ij).array();
                        Rdd.col(ijk) = TTerm.col(tij).array() * nonDose_LOGLIN.col(tij).array().pow(-1).array() * Tdd0.col(ijk).array();
                    }
                } else if (tij==tjk){
                    if (Dose_Iden.find(as< string>(tform[ij])) != Dose_Iden.end()){
                        if (Dose_Iden.find(as< string>(tform[jk])) != Dose_Iden.end()){
                            Rdd.col(ijk) = TTerm.col(tij).array() * Dose.col(tij).array().pow(-1).array() * Tdd0.col(ijk).array();
                        } else if (tform[jk]=="lin") {
                            Rdd.col(ijk) = TTerm.col(tij).array() * Dose.col(tij).array().pow(-1).array() * Td0.col(ij).array() * nonDose_LIN.col(tij).array().pow(-1).array()    * Td0.col(jk).array();
                        } else if (tform[jk]=="plin") {
                            Rdd.col(ijk) = TTerm.col(tij).array() * Dose.col(tij).array().pow(-1).array() * Td0.col(ij).array() * nonDose_PLIN.col(tij).array().pow(-1).array()   * Td0.col(jk).array();
                        } else if (tform[jk]=="loglin") {
                            Rdd.col(ijk) = TTerm.col(tij).array() * Dose.col(tij).array().pow(-1).array() * Td0.col(ij).array() * nonDose_LOGLIN.col(tij).array().pow(-1).array() * Td0.col(jk).array();
                        }
                    } else if (Dose_Iden.find(as< string>(tform[jk])) != Dose_Iden.end()){
                        if (tform[ij]=="lin") {
                            Rdd.col(ijk) = TTerm.col(tij).array() * Dose.col(tij).array().pow(-1).array() * Td0.col(jk).array() * nonDose_LIN.col(tij).array().pow(-1).array()    * Td0.col(ij).array();
                        } else if (tform[ij]=="plin") {
                            Rdd.col(ijk) = TTerm.col(tij).array() * Dose.col(tij).array().pow(-1).array() * Td0.col(jk).array() * nonDose_PLIN.col(tij).array().pow(-1).array()   * Td0.col(ij).array();
                        } else if (tform[ij]=="loglin") {
                            Rdd.col(ijk) = TTerm.col(tij).array() * Dose.col(tij).array().pow(-1).array() * Td0.col(jk).array() * nonDose_LOGLIN.col(tij).array().pow(-1).array() * Td0.col(ij).array();
                        }
                    } else if (tform[ij]=="loglin") {
                        if( tform[jk]=="lin") {
                            Rdd.col(ijk) = TTerm.col(tij).array() * nonDose_LOGLIN.col(tij).array().pow(-1).array() * Td0.col(ij).array() * nonDose_LIN.col(tij).array().pow(-1).array() * Td0.col(jk).array();
                        } else if (tform[jk]=="plin") {
                            Rdd.col(ijk) = TTerm.col(tij).array() * nonDose_LOGLIN.col(tij).array().pow(-1).array() * Td0.col(ij).array() * nonDose_PLIN.col(tij).array().pow(-1).array() * Td0.col(jk).array();
                        } else if (tform[jk]=="loglin") {
                            Rdd.col(ijk) = TTerm.col(tij).array() * nonDose_LOGLIN.col(tij).array().pow(-1).array() * Tdd0.col(ijk).array();
                        }
                    } else if (tform[jk]=="loglin") {
                        if( tform[ij]=="lin") {
                            Rdd.col(ijk) = TTerm.col(tij).array() * nonDose_LOGLIN.col(tjk).array().pow(-1).array() * Td0.col(jk).array() * nonDose_LIN.col(tij).array().pow(-1).array() * Td0.col(ij).array();
                        } else if (tform[ij]=="plin") {
                            Rdd.col(ijk) = TTerm.col(tij).array() * nonDose_LOGLIN.col(tjk).array().pow(-1).array() * Td0.col(jk).array() * nonDose_PLIN.col(tij).array().pow(-1).array() * Td0.col(ij).array();
                        }
                    } else if (tform[ij]=="lin") {
                        if( tform[jk]=="lin") {
                            ;
                        } else if (tform[jk]=="plin") {
                            Rdd.col(ijk) = TTerm.col(tij).array() * nonDose_LIN.col(tij).array().pow(-1).array() * Td0.col(ij).array() * nonDose_PLIN.col(tij).array().pow(-1).array() * Td0.col(jk).array();
                        }
                    } else if (tform[jk]=="lin") {
                        if (tform[ij]=="plin") {
                            Rdd.col(ijk) = TTerm.col(tij).array() * nonDose_LIN.col(tjk).array().pow(-1).array() * Td0.col(jk).array() * nonDose_PLIN.col(tij).array().pow(-1).array() * Td0.col(ij).array();
                        }
                    } else {
                        ;
                    }
                }
            }
        } else if ((modelform=="PAE")||(modelform=="PA")){
            Te = Te.array() - TTerm.col(fir).array();
            if (modelform=="PAE"){
                Te = Te.array() + 1;
            }
            R << TTerm.col(fir).array() * Te.array();
            #pragma omp parallel for schedule(dynamic) num_threads(nthreads)
            for (int ijk=0;ijk<totalnum*(totalnum+1)/2;ijk++){
                int ij = 0;
                int jk = ijk;
                while (jk>ij){
                    ij++;
                    jk-=ij;
                }
                int tij = Term_n[ij];
                int tjk = Term_n[jk];
                if (ij==jk){
                    if (tij==fir){
                        if (Dose_Iden.find(as< string>(tform[ij])) != Dose_Iden.end()){
                            Rd.col(ij) =   R.col(0).array() * Dose.col(tij).array().pow(-1).array() * Td0.col(ij).array();
                            Rdd.col(ij) =  R.col(0).array() * Dose.col(tij).array().pow(-1).array() * Tdd0.col(ijk).array();
                        } else if (tform[ij]=="lin") {
                            Rd.col(ij) =   R.col(0).array() * nonDose_LIN.col(tij).array().pow(-1).array() * Td0.col(ij).array();
                            Rdd.col(ij) =  R.col(0).array() * nonDose_LIN.col(tij).array().pow(-1).array() * Tdd0.col(ijk).array();
                        } else if (tform[ij]=="plin") {
                            Rd.col(ij) =   R.col(0).array() * nonDose_PLIN.col(tij).array().pow(-1).array() * Td0.col(ij).array();
                            Rdd.col(ij) =  R.col(0).array() * nonDose_PLIN.col(tij).array().pow(-1).array() * Tdd0.col(ijk).array();
                        } else if (tform[ij]=="loglin") {
                            Rd.col(ij) =   R.col(0).array() * nonDose_LOGLIN.col(tij).array().pow(-1).array() * Td0.col(ij).array();
                            Rdd.col(ij) =  R.col(0).array() * nonDose_LOGLIN.col(tij).array().pow(-1).array() * Tdd0.col(ijk).array();
                        }
                    } else {
                        if (Dose_Iden.find(as< string>(tform[ij])) != Dose_Iden.end()){
                            Rd.col(ij) =   TTerm.col(fir).array() * TTerm.col(tij).array() * Dose.col(tij).array().pow(-1).array() * Td0.col(ij).array();
                            Rdd.col(ijk) = TTerm.col(fir).array() * TTerm.col(tij).array() * Dose.col(tij).array().pow(-1).array() * Tdd0.col(ijk).array();
                        } else if (tform[ij]=="lin") {
                            Rd.col(ij) =   TTerm.col(fir).array() * TTerm.col(tij).array() * nonDose_LIN.col(tij).array().pow(-1).array() * Td0.col(ij).array();
                            Rdd.col(ijk) = TTerm.col(fir).array() * TTerm.col(tij).array() * nonDose_LIN.col(tij).array().pow(-1).array() * Tdd0.col(ijk).array();
                        } else if (tform[ij]=="plin") {
                            Rd.col(ij) =   TTerm.col(fir).array() * TTerm.col(tij).array() * nonDose_PLIN.col(tij).array().pow(-1).array() * Td0.col(ij).array();
                            Rdd.col(ijk) = TTerm.col(fir).array() * TTerm.col(tij).array() * nonDose_PLIN.col(tij).array().pow(-1).array() * Tdd0.col(ijk).array();
                        } else if (tform[ij]=="loglin") {
                            Rd.col(ij) =   TTerm.col(fir).array() * TTerm.col(tij).array() * nonDose_LOGLIN.col(tij).array().pow(-1).array() * Td0.col(ij).array();
                            Rdd.col(ijk) = TTerm.col(fir).array() * TTerm.col(tij).array() * nonDose_LOGLIN.col(tij).array().pow(-1).array() * Tdd0.col(ijk).array();
                        }
                    }
                } else {
                    if (tij==tjk){
                        if (tij==fir){
                            if (Dose_Iden.find(as< string>(tform[ij])) != Dose_Iden.end()){
                                if (Dose_Iden.find(as< string>(tform[jk])) != Dose_Iden.end()){
                                    Rdd.col(ijk) = R.col(0).array() * Dose.col(tij).array().pow(-1).array() * Tdd0.col(ijk).array();
                                } else if (tform[jk]=="lin") {
                                    Rdd.col(ijk) = R.col(0).array() * Dose.col(tij).array().pow(-1).array() * Td0.col(ij).array() * nonDose_LIN.col(tij).array().pow(-1).array()    * Td0.col(jk).array();
                                } else if (tform[jk]=="plin") {
                                    Rdd.col(ijk) = R.col(0).array() * Dose.col(tij).array().pow(-1).array() * Td0.col(ij).array() * nonDose_PLIN.col(tij).array().pow(-1).array()   * Td0.col(jk).array();
                                } else if (tform[jk]=="loglin") {
                                    Rdd.col(ijk) = R.col(0).array() * Dose.col(tij).array().pow(-1).array() * Td0.col(ij).array() * nonDose_LOGLIN.col(tij).array().pow(-1).array() * Td0.col(jk).array();
                                }
                            } else if (Dose_Iden.find(as< string>(tform[jk])) != Dose_Iden.end()){
                                if (tform[ij]=="lin") {
                                    Rdd.col(ijk) = R.col(0).array() * Dose.col(tij).array().pow(-1).array() * Td0.col(jk).array() * nonDose_LIN.col(tij).array().pow(-1).array()    * Td0.col(ij).array();
                                } else if (tform[ij]=="plin") {
                                    Rdd.col(ijk) = R.col(0).array() * Dose.col(tij).array().pow(-1).array() * Td0.col(jk).array() * nonDose_PLIN.col(tij).array().pow(-1).array()   * Td0.col(ij).array();
                                } else if (tform[ij]=="loglin") {
                                    Rdd.col(ijk) = R.col(0).array() * Dose.col(tij).array().pow(-1).array() * Td0.col(jk).array() * nonDose_LOGLIN.col(tij).array().pow(-1).array() * Td0.col(ij).array();
                                }
                            } else if (tform[ij]=="loglin") {
                                if( tform[jk]=="lin") {
                                    Rdd.col(ijk) = R.col(0).array() * nonDose_LOGLIN.col(tij).array().pow(-1).array() * Td0.col(ij).array() * nonDose_LIN.col(tij).array().pow(-1).array() * Td0.col(jk).array();
                                } else if (tform[jk]=="plin") {
                                    Rdd.col(ijk) = R.col(0).array() * nonDose_LOGLIN.col(tij).array().pow(-1).array() * Td0.col(ij).array() * nonDose_PLIN.col(tij).array().pow(-1).array() * Td0.col(jk).array();
                                } else if (tform[jk]=="loglin") {
                                    Rdd.col(ijk) = R.col(0).array() * nonDose_LOGLIN.col(tij).array().pow(-1).array() * Td0.col(ij).array() * nonDose_LOGLIN.col(tjk).array().pow(-1).array() * Td0.col(jk).array();
                                }
                            } else if (tform[jk]=="loglin") {
                                if( tform[ij]=="lin") {
                                    Rdd.col(ijk) = R.col(0).array() * nonDose_LOGLIN.col(tjk).array().pow(-1).array() * Td0.col(jk).array() * nonDose_LIN.col(tij).array().pow(-1).array() * Td0.col(ij).array();
                                } else if (tform[ij]=="plin") {
                                    Rdd.col(ijk) = R.col(0).array() * nonDose_LOGLIN.col(tjk).array().pow(-1).array() * Td0.col(jk).array() * nonDose_PLIN.col(tij).array().pow(-1).array() * Td0.col(ij).array();
                                }
                            } else if (tform[ij]=="lin") {
                                if( tform[jk]=="lin") {
                                    Rdd.col(ijk) = R.col(0).array() * nonDose_LIN.col(tij).array().pow(-1).array() * Td0.col(ij).array() * nonDose_LIN.col(tij).array().pow(-1).array() * Td0.col(jk).array();
                                } else if (tform[jk]=="plin") {
                                    Rdd.col(ijk) = R.col(0).array() * nonDose_LIN.col(tij).array().pow(-1).array() * Td0.col(ij).array() * nonDose_PLIN.col(tij).array().pow(-1).array() * Td0.col(jk).array();
                                }
                            } else if (tform[jk]=="lin") {
                                if (tform[ij]=="plin") {
                                    Rdd.col(ijk) = R.col(0).array() * nonDose_LIN.col(tjk).array().pow(-1).array() * Td0.col(jk).array() * nonDose_PLIN.col(tij).array().pow(-1).array() * Td0.col(ij).array();
                                }
                            } else {
                                Rdd.col(ijk) = R.col(0).array() * nonDose_PLIN.col(tjk).array().pow(-1).array() * Td0.col(jk).array() * nonDose_PLIN.col(tij).array().pow(-1).array() * Td0.col(ij).array();
                            }
                        } else {
                            if (Dose_Iden.find(as< string>(tform[ij])) != Dose_Iden.end()){
                                if (Dose_Iden.find(as< string>(tform[jk])) != Dose_Iden.end()){
                                    Rdd.col(ijk) = TTerm.col(fir).array() * TTerm.col(tij).array() * Dose.col(tij).array().pow(-1).array() * Tdd0.col(ijk).array();
                                } else if (tform[jk]=="lin") {
                                    Rdd.col(ijk) = TTerm.col(fir).array() * TTerm.col(tij).array() * Dose.col(tij).array().pow(-1).array() * Td0.col(ij).array() * nonDose_LIN.col(tij).array().pow(-1).array()    * Td0.col(jk).array();
                                } else if (tform[jk]=="plin") {
                                    Rdd.col(ijk) = TTerm.col(fir).array() * TTerm.col(tij).array() * Dose.col(tij).array().pow(-1).array() * Td0.col(ij).array() * nonDose_PLIN.col(tij).array().pow(-1).array()   * Td0.col(jk).array();
                                } else if (tform[jk]=="loglin") {
                                    Rdd.col(ijk) = TTerm.col(fir).array() * TTerm.col(tij).array() * Dose.col(tij).array().pow(-1).array() * Td0.col(ij).array() * nonDose_LOGLIN.col(tij).array().pow(-1).array() * Td0.col(jk).array();
                                }
                            } else if (Dose_Iden.find(as< string>(tform[jk])) != Dose_Iden.end()){
                                if (tform[ij]=="lin") {
                                    Rdd.col(ijk) = TTerm.col(fir).array() * TTerm.col(tij).array() * Dose.col(tij).array().pow(-1).array() * Td0.col(jk).array() * nonDose_LIN.col(tij).array().pow(-1).array()    * Td0.col(ij).array();
                                } else if (tform[ij]=="plin") {
                                    Rdd.col(ijk) = TTerm.col(fir).array() * TTerm.col(tij).array() * Dose.col(tij).array().pow(-1).array() * Td0.col(jk).array() * nonDose_PLIN.col(tij).array().pow(-1).array()   * Td0.col(ij).array();
                                } else if (tform[ij]=="loglin") {
                                    Rdd.col(ijk) = TTerm.col(fir).array() * TTerm.col(tij).array() * Dose.col(tij).array().pow(-1).array() * Td0.col(jk).array() * nonDose_LOGLIN.col(tij).array().pow(-1).array() * Td0.col(ij).array();
                                }
                            } else if (tform[ij]=="loglin") {
                                if( tform[jk]=="lin") {
                                    Rdd.col(ijk) = TTerm.col(fir).array() * TTerm.col(tij).array() * nonDose_LOGLIN.col(tij).array().pow(-1).array() * Td0.col(ij).array() * nonDose_LIN.col(tij).array().pow(-1).array() * Td0.col(jk).array();
                                } else if (tform[jk]=="plin") {
                                    Rdd.col(ijk) = TTerm.col(fir).array() * TTerm.col(tij).array() * nonDose_LOGLIN.col(tij).array().pow(-1).array() * Td0.col(ij).array() * nonDose_PLIN.col(tij).array().pow(-1).array() * Td0.col(jk).array();
                                } else if (tform[jk]=="loglin") {
                                    Rdd.col(ijk) = TTerm.col(fir).array() * TTerm.col(tij).array() * nonDose_LOGLIN.col(tij).array().pow(-1).array() * Td0.col(ij).array() * nonDose_LOGLIN.col(tjk).array().pow(-1).array() * Td0.col(jk).array();
                                }
                            } else if (tform[jk]=="loglin") {
                                if( tform[ij]=="lin") {
                                    Rdd.col(ijk) = TTerm.col(fir).array() * TTerm.col(tij).array() * nonDose_LOGLIN.col(tjk).array().pow(-1).array() * Td0.col(jk).array() * nonDose_LIN.col(tij).array().pow(-1).array() * Td0.col(ij).array();
                                } else if (tform[ij]=="plin") {
                                    Rdd.col(ijk) = TTerm.col(fir).array() * TTerm.col(tij).array() * nonDose_LOGLIN.col(tjk).array().pow(-1).array() * Td0.col(jk).array() * nonDose_PLIN.col(tij).array().pow(-1).array() * Td0.col(ij).array();
                                }
                            } else if (tform[ij]=="lin") {
                                if( tform[jk]=="lin") {
                                    Rdd.col(ijk) = TTerm.col(fir).array() * TTerm.col(tij).array() * nonDose_LIN.col(tij).array().pow(-1).array() * Td0.col(ij).array() * nonDose_LIN.col(tij).array().pow(-1).array() * Td0.col(jk).array();
                                } else if (tform[jk]=="plin") {
                                    Rdd.col(ijk) = TTerm.col(fir).array() * TTerm.col(tij).array() * nonDose_LIN.col(tij).array().pow(-1).array() * Td0.col(ij).array() * nonDose_PLIN.col(tij).array().pow(-1).array() * Td0.col(jk).array();
                                }
                            } else if (tform[jk]=="lin") {
                                if (tform[ij]=="plin") {
                                    Rdd.col(ijk) = TTerm.col(fir).array() * TTerm.col(tij).array() * nonDose_LIN.col(tjk).array().pow(-1).array() * Td0.col(jk).array() * nonDose_PLIN.col(tij).array().pow(-1).array() * Td0.col(ij).array();
                                }
                            } else {
                                Rdd.col(ijk) = TTerm.col(fir).array() * TTerm.col(tij).array() * nonDose_PLIN.col(tjk).array().pow(-1).array() * Td0.col(jk).array() * nonDose_PLIN.col(tij).array().pow(-1).array() * Td0.col(ij).array();
                            }
                        }  
                    } else if ((tij==fir)||(tjk==fir)){
                        if (Dose_Iden.find(as< string>(tform[ij])) != Dose_Iden.end()){
                            if (Dose_Iden.find(as< string>(tform[jk])) != Dose_Iden.end()){
                                Rdd.col(ijk) = TTerm.col(tjk).array() * TTerm.col(tij).array() * Dose.col(tij).array().pow(-1).array() * Tdd0.col(ijk).array();
                            } else if (tform[jk]=="lin") {
                                Rdd.col(ijk) = TTerm.col(tjk).array() * TTerm.col(tij).array() * Dose.col(tij).array().pow(-1).array() * Td0.col(ij).array() * nonDose_LIN.col(tij).array().pow(-1).array()    * Td0.col(jk).array();
                            } else if (tform[jk]=="plin") {
                                Rdd.col(ijk) = TTerm.col(tjk).array() * TTerm.col(tij).array() * Dose.col(tij).array().pow(-1).array() * Td0.col(ij).array() * nonDose_PLIN.col(tij).array().pow(-1).array()   * Td0.col(jk).array();
                            } else if (tform[jk]=="loglin") {
                                Rdd.col(ijk) = TTerm.col(tjk).array() * TTerm.col(tij).array() * Dose.col(tij).array().pow(-1).array() * Td0.col(ij).array() * nonDose_LOGLIN.col(tij).array().pow(-1).array() * Td0.col(jk).array();
                            }
                        } else if (Dose_Iden.find(as< string>(tform[jk])) != Dose_Iden.end()){
                            if (tform[ij]=="lin") {
                                Rdd.col(ijk) = TTerm.col(tjk).array() * TTerm.col(tij).array() * Dose.col(tij).array().pow(-1).array() * Td0.col(jk).array() * nonDose_LIN.col(tij).array().pow(-1).array()    * Td0.col(ij).array();
                            } else if (tform[ij]=="plin") {
                                Rdd.col(ijk) = TTerm.col(tjk).array() * TTerm.col(tij).array() * Dose.col(tij).array().pow(-1).array() * Td0.col(jk).array() * nonDose_PLIN.col(tij).array().pow(-1).array()   * Td0.col(ij).array();
                            } else if (tform[ij]=="loglin") {
                                Rdd.col(ijk) = TTerm.col(tjk).array() * TTerm.col(tij).array() * Dose.col(tij).array().pow(-1).array() * Td0.col(jk).array() * nonDose_LOGLIN.col(tij).array().pow(-1).array() * Td0.col(ij).array();
                            }
                        } else if (tform[ij]=="loglin") {
                            if( tform[jk]=="lin") {
                                Rdd.col(ijk) = TTerm.col(tjk).array() * TTerm.col(tij).array() * nonDose_LOGLIN.col(tij).array().pow(-1).array() * Td0.col(ij).array() * nonDose_LIN.col(tij).array().pow(-1).array() * Td0.col(jk).array();
                            } else if (tform[jk]=="plin") {
                                Rdd.col(ijk) = TTerm.col(tjk).array() * TTerm.col(tij).array() * nonDose_LOGLIN.col(tij).array().pow(-1).array() * Td0.col(ij).array() * nonDose_PLIN.col(tij).array().pow(-1).array() * Td0.col(jk).array();
                            } else if (tform[jk]=="loglin") {
                                Rdd.col(ijk) = TTerm.col(tjk).array() * TTerm.col(tij).array() * nonDose_LOGLIN.col(tij).array().pow(-1).array() * Td0.col(ij).array() * nonDose_LOGLIN.col(tjk).array().pow(-1).array() * Td0.col(jk).array();
                            }
                        } else if (tform[jk]=="loglin") {
                            if( tform[ij]=="lin") {
                                Rdd.col(ijk) = TTerm.col(tjk).array() * TTerm.col(tij).array() * nonDose_LOGLIN.col(tjk).array().pow(-1).array() * Td0.col(jk).array() * nonDose_LIN.col(tij).array().pow(-1).array() * Td0.col(ij).array();
                            } else if (tform[ij]=="plin") {
                                Rdd.col(ijk) = TTerm.col(tjk).array() * TTerm.col(tij).array() * nonDose_LOGLIN.col(tjk).array().pow(-1).array() * Td0.col(jk).array() * nonDose_PLIN.col(tij).array().pow(-1).array() * Td0.col(ij).array();
                            }
                        } else if (tform[ij]=="lin") {
                            if( tform[jk]=="lin") {
                                Rdd.col(ijk) = TTerm.col(tjk).array() * TTerm.col(tij).array() * nonDose_LIN.col(tij).array().pow(-1).array() * Td0.col(ij).array() * nonDose_LIN.col(tij).array().pow(-1).array() * Td0.col(jk).array();
                            } else if (tform[jk]=="plin") {
                                Rdd.col(ijk) = TTerm.col(tjk).array() * TTerm.col(tij).array() * nonDose_LIN.col(tij).array().pow(-1).array() * Td0.col(ij).array() * nonDose_PLIN.col(tij).array().pow(-1).array() * Td0.col(jk).array();
                            }
                        } else if (tform[jk]=="lin") {
                            if (tform[ij]=="plin") {
                                Rdd.col(ijk) = TTerm.col(tjk).array() * TTerm.col(tij).array() * nonDose_LIN.col(tjk).array().pow(-1).array() * Td0.col(jk).array() * nonDose_PLIN.col(tij).array().pow(-1).array() * Td0.col(ij).array();
                            }
                        } else {
                            Rdd.col(ijk) = TTerm.col(tjk).array() * TTerm.col(tij).array() * nonDose_PLIN.col(tjk).array().pow(-1).array() * Td0.col(jk).array() * nonDose_PLIN.col(tij).array().pow(-1).array() * Td0.col(ij).array();
                        }
                    }
                }
            }
        }
    }else if (modelform=="M"){
        //
        MatrixXd TTerm_p(TTerm.rows(),TTerm.cols());
        TTerm_p << TTerm.array() + 1.0;
        TTerm_p.col(fir) = TTerm_p.col(fir).array() -1.0;
        Te = TTerm_p.array().rowwise().prod().array();
        R << Te.array();
        //
        //
        Rd = Td0.array();
        
        #pragma omp parallel for schedule(dynamic) num_threads(nthreads)
        for (int ijk=0;ijk<totalnum;ijk++){
            int tij = Term_n[ijk];
            if (Dose_Iden.find(as< string>(tform[ijk])) != Dose_Iden.end()){
                Rd.col(ijk) = R.col(0).array() * TTerm_p.col(tij).array().pow(-1).array() * Td0.array().col(ijk).array() * TTerm.col(tij).array() * Dose.col(tij).array().pow(-1).array();
            } else if (tform[ijk]=="lin") {
                Rd.col(ijk) = R.col(0).array() * TTerm_p.col(tij).array().pow(-1).array() * Td0.array().col(ijk).array() * TTerm.col(tij).array() * nonDose_LIN.col(tij).array().pow(-1).array();
            } else if (tform[ijk]=="plin") {
                Rd.col(ijk) = R.col(0).array() * TTerm_p.col(tij).array().pow(-1).array() * Td0.array().col(ijk).array() * TTerm.col(tij).array() * nonDose_PLIN.col(tij).array().pow(-1).array();
            } else if (tform[ijk]=="loglin") {
                Rd.col(ijk) = R.col(0).array() * TTerm_p.col(tij).array().pow(-1).array() * Td0.array().col(ijk).array() * TTerm.col(tij).array() * nonDose_LOGLIN.col(tij).array().pow(-1).array() ;
            }
        }
        R = (R.array().isFinite()).select(R,0);
        Rd = (Rd.array().isFinite()).select(Rd,0);
        #pragma omp parallel for schedule(dynamic) num_threads(nthreads)
        for (int ijk=0;ijk<totalnum*(totalnum+1)/2;ijk++){
            int ij = 0;
            int jk = ijk;
            while (jk>ij){
                ij++;
                jk-=ij;
            }
            int tij = Term_n[ij];
            int tjk = Term_n[jk];
            if (ij==jk){
                if (Dose_Iden.find(as< string>(tform[ij])) != Dose_Iden.end()){
                    Rdd.col(ijk) = R.col(0).array() * TTerm.col(tij).array() * TTerm_p.col(tij).array().pow(-1).array() * Tdd0.array().col(ijk).array() * Dose.col(tij).array().pow(-1).array();
                } else if (tform[ij]=="lin") {
                    Rdd.col(ijk) = R.col(0).array() * TTerm.col(tij).array() * TTerm_p.col(tij).array().pow(-1).array() * Tdd0.array().col(ijk).array() * nonDose_LIN.col(tij).array().pow(-1).array();
                } else if (tform[ij]=="plin") {
                    Rdd.col(ijk) = R.col(0).array() * TTerm.col(tij).array() * TTerm_p.col(tij).array().pow(-1).array() * Tdd0.array().col(ijk).array() * nonDose_PLIN.col(tij).array().pow(-1).array();
                } else if (tform[ij]=="loglin") {
                    Rdd.col(ijk) = R.col(0).array() * TTerm.col(tij).array() * TTerm_p.col(tij).array().pow(-1).array() * Tdd0.array().col(ijk).array() * nonDose_LOGLIN.col(tij).array().pow(-1).array() ;
                }
            } else {
                Rdd.col(ijk) = R.col(0).array() * TTerm.col(tij).array() * TTerm_p.col(tij).array().pow(-1).array() * TTerm.col(tjk).array() * TTerm_p.col(tjk).array().pow(-1).array();
                if (Dose_Iden.find(as< string>(tform[ij])) != Dose_Iden.end()){
                    Rdd.col(ijk) = Rdd.col(ijk).array() * Td0.array().col(ij).array() * Dose.col(tij).array().pow(-1).array();
                } else if (tform[ij]=="lin") {
                    Rdd.col(ijk) = Rdd.col(ijk).array() * Td0.array().col(ij).array() * nonDose_LIN.col(tij).array().pow(-1).array();
                } else if (tform[ij]=="plin") {
                    Rdd.col(ijk) = Rdd.col(ijk).array() * Td0.array().col(ij).array() * nonDose_PLIN.col(tij).array().pow(-1).array();
                } else if (tform[ij]=="loglin") {
                    Rdd.col(ijk) = Rdd.col(ijk).array() * Td0.array().col(ij).array() * nonDose_LOGLIN.col(tij).array().pow(-1).array() ;
                }
                //
                if (Dose_Iden.find(as< string>(tform[jk])) != Dose_Iden.end()){
                    Rdd.col(ijk) = Rdd.col(ijk).array() * Td0.array().col(jk).array() * Dose.col(tjk).array().pow(-1).array();
                } else if (tform[jk]=="lin") {
                    Rdd.col(ijk) = Rdd.col(ijk).array() * Td0.array().col(jk).array() * nonDose_LIN.col(tjk).array().pow(-1).array();
                } else if (tform[jk]=="plin") {
                    Rdd.col(ijk) = Rdd.col(ijk).array() * Td0.array().col(jk).array() * nonDose_PLIN.col(tjk).array().pow(-1).array();
                } else if (tform[jk]=="loglin") {
                    Rdd.col(ijk) = Rdd.col(ijk).array() * Td0.array().col(jk).array() * nonDose_LOGLIN.col(tjk).array().pow(-1).array() ;
                }
                //
            }
        }
    } else if (modelform=="GM"){
        //currently isn't implemented, it can be calculated but not optimized the same way
        throw invalid_argument( "GM isn't implemented" );
    } else {
        throw invalid_argument( "Model isn't implemented" );
    }
    //
    //
    R = (R.array().isFinite()).select(R,0);
    Rd = (Rd.array().isFinite()).select(Rd,0);
    Rdd = (Rdd.array().isFinite()).select(Rdd,0);
    //
    for (int ijk=0;ijk<(totalnum*(totalnum+1)/2);ijk++){//Calculates ratios
        int ij = 0;
        int jk = ijk;
        while (jk>ij){
            ij++;
            jk-=ij;
        }
        if (ij==jk){
            RdR.col(ij)=R.col(0).array().pow(-1).array() * Rd.col(jk).array();
        }
        RddR.col(ijk)=R.col(0).array().pow(-1).array() * Rdd.col(ijk).array();
    }
    return;
}

void Make_Groups(const int& ntime, const MatrixXd& df_m, IntegerMatrix& RiskFail, vector<string>&  RiskGroup,  NumericVector& tu, const int& nthreads, bool debugging ){
    //
    //Make_Subterms( ntime, df_m, RiskFail, RiskGroup, tu, nthreads, debugging)
    //
    if (debugging){
        cout << "Starting Debug" << endl;
        vector<double> time_ref(6,0.0);
        vector<double> time_refs(6,0.0);
        cout << time_ref.size() << " " << time_refs.size() << endl;
        #pragma omp declare reduction(vec_double_plus : std::vector<double> : \
            std::transform(omp_out.begin(), omp_out.end(), omp_in.begin(), omp_out.begin(), std::plus<double>())) \
            initializer(omp_priv = omp_orig)
        #pragma omp parallel for schedule(dynamic) num_threads(nthreads) reduction(vec_double_plus:time_ref,time_refs) 
        for (int ijk=0;ijk<ntime;ijk++){
            //
            time_point<system_clock> start_point, end_point;
            start_point = system_clock::now();
            auto start = time_point_cast<microseconds>(start_point).time_since_epoch().count();
            end_point = system_clock::now();
            auto ending = time_point_cast<microseconds>(end_point).time_since_epoch().count(); //The time duration is tracked
            //
            double t0 = tu[ijk];
            VectorXi select_ind_all = ((df_m.col(0).array() <= t0)&&(df_m.col(1).array()>=t0)).cast<int>(); //indices at risk
            vector<int> indices_all;
            //
            end_point = system_clock::now();
            ending = time_point_cast<microseconds>(end_point).time_since_epoch().count();
            time_ref[0]+=ending-start;
            time_refs[0]+=pow(ending-start,2);
            start=ending;
            //
            VectorXi select_ind_end = ((df_m.col(2).array() == 1)&&((df_m.col(1).array()==t0))||(df_m.col(0).array()==t0)).cast<int>(); //indices with events
            vector<int> indices_end;
            //
            end_point = system_clock::now();
            ending = time_point_cast<microseconds>(end_point).time_since_epoch().count();
            time_ref[1]+=ending-start;
            time_refs[1]+=pow(ending-start,2);
            start=ending;
            //
            //
            //
            int th = 1;
            visit_lambda(select_ind_all,
                [&indices_all, th](double v, int i, int j) {
                    if (v==th)
                        indices_all.push_back(i+1);
                });
            //
            end_point = system_clock::now();
            ending = time_point_cast<microseconds>(end_point).time_since_epoch().count();
            time_ref[2]+=ending-start;
            time_refs[2]+=pow(ending-start,2);
            start=ending;
            //
            visit_lambda(select_ind_end,
                [&indices_end, th](double v, int i, int j) {
                    if (v==th)
                        indices_end.push_back(i+1);
                });
            //
            end_point = system_clock::now();
            ending = time_point_cast<microseconds>(end_point).time_since_epoch().count();
            time_ref[3]+=ending-start;
            time_refs[3]+=pow(ending-start,2);
            start=ending;
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
            //
            end_point = system_clock::now();
            ending = time_point_cast<microseconds>(end_point).time_since_epoch().count();
            time_ref[4]+=ending-start;
            time_refs[4]+=pow(ending-start,2);
            start=ending;
            //
            RiskFail(ijk,0)=indices_end[0]-1;//Due to the sorting method, there is a continuous block of event rows
            RiskFail(ijk,1)=indices_end[indices_end.size()-1]-1;
            //
            ostringstream oss;
            copy(indices.begin(), indices.end(),
                std::ostream_iterator<int>(oss, ","));
            RiskGroup[ijk] = oss.str();//stores risk groups in string
            //
            end_point = system_clock::now();
            ending = time_point_cast<microseconds>(end_point).time_since_epoch().count();
            time_ref[5]+=ending-start;
            time_refs[5]+=pow(ending-start,2);
            start=ending;
        }
        cout << "df204 ";
        for (int ijk=0;ijk<time_ref.size();ijk++){
            cout << (time_ref[ijk]/ntime)*1e-6 << " ";
        }
        cout << " " << endl;
        cout << "df205 ";
        for (int ijk=0;ijk<time_ref.size();ijk++){
            cout <<  sqrt(time_refs[ijk]/ntime - pow(time_ref[ijk]/ntime,2))*1e-6 << " ";
        }
        cout << " " << endl;
    } else {
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
            RiskFail(ijk,0)=indices_end[0]-1;//Due to the sorting method, there is a continuous block of event rows
            RiskFail(ijk,1)=indices_end[indices_end.size()-1]-1;
            //
            ostringstream oss;
            copy(indices.begin(), indices.end(),
                std::ostream_iterator<int>(oss, ","));
            RiskGroup[ijk] = oss.str();//stores risk groups in string
        }
    }
    return;
}

void Calculate_Sides(const IntegerMatrix& RiskFail, const vector<string>&  RiskGroup, const int& totalnum, const int& ntime, const MatrixXd& R, const MatrixXd& Rd, const MatrixXd& Rdd, MatrixXd& Rls1, MatrixXd& Rls2, MatrixXd& Rls3, MatrixXd& Lls1, MatrixXd& Lls2, MatrixXd& Lls3,const int& nthreads, bool debugging){
    //
    //Calculate_Sides( RiskFail, RiskGroup, totalnum, ntime, R, Rd, Rdd, Rls1, Rls2, Rls3, Lls1, Lls2, Lls3,nthreads, debugging);
    //
    if (debugging){
        cout << "Starting Debug" << endl;
        vector<double> time_ref(3,0.0);
        vector<double> time_refs(3,0.0);
        vector<int> time_count(3,0);
        cout << time_ref.size() << " " << time_refs.size() << endl;
        #pragma omp declare reduction(vec_double_plus : std::vector<double> : \
            std::transform(omp_out.begin(), omp_out.end(), omp_in.begin(), omp_out.begin(), std::plus<double>())) \
            initializer(omp_priv = omp_orig)
        #pragma omp declare reduction(vec_int_plus : std::vector<int> : \
            std::transform(omp_out.begin(), omp_out.end(), omp_in.begin(), omp_out.begin(), std::plus<int>())) \
            initializer(omp_priv = omp_orig)
        //
        #pragma omp parallel for schedule(dynamic) num_threads(nthreads) collapse(2) reduction(vec_double_plus:time_ref,time_refs) reduction(vec_int_plus:time_count)
        for (int ijk=0;ijk<totalnum*(totalnum+1)/2;ijk++){//totalnum*(totalnum+1)/2
            for (int j=0;j<ntime;j++){
                time_point<system_clock> start_point, end_point;
                start_point = system_clock::now();
                auto start = time_point_cast<microseconds>(start_point).time_since_epoch().count();
                end_point = system_clock::now();
                auto ending = time_point_cast<microseconds>(end_point).time_since_epoch().count(); //The time duration is tracked
                int ij = 0;
                int jk = ijk;
                while (jk>ij){
                    ij++;
                    jk-=ij;
                }
                double Rs1 = 0;
                double Rs2 = 0;
                double Rs2t = 0;
                double Rs3 = 0;
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
                //
                end_point = system_clock::now();
                ending = time_point_cast<microseconds>(end_point).time_since_epoch().count();
                time_ref[0]+=ending-start;
                time_refs[0]+=pow(ending-start,2);
                time_count[0]+=1;
                start=ending;
                //
                //Now has the grouping pairs
                int dj = RiskFail(j,1)-RiskFail(j,0)+1;
                for (int i = 0; i < InGroup.size()-1; i=i+2){
                    Rs1 += R.block(InGroup[i]-1,0,InGroup[i+1]-InGroup[i]+1,1).sum();
                    Rs2 += Rd.block(InGroup[i]-1,ij,InGroup[i+1]-InGroup[i]+1,1).sum();
                    Rs2t += Rd.block(InGroup[i]-1,jk,InGroup[i+1]-InGroup[i]+1,1).sum();
                    Rs3 += Rdd.block(InGroup[i]-1,ijk,InGroup[i+1]-InGroup[i]+1,1).sum();
                } //precalculates the sums of risk groups
                MatrixXd Ld = MatrixXd::Zero(dj,4);
                Ld << R.block(RiskFail(j,0),0,dj,1), Rd.block(RiskFail(j,0),ij,dj,1), Rd.block(RiskFail(j,0),jk,dj,1) ,Rdd.block(RiskFail(j,0),ijk,dj,1);//sum of risks in group
                //
                end_point = system_clock::now();
                ending = time_point_cast<microseconds>(end_point).time_since_epoch().count();
                time_ref[1]+=ending-start;
                time_refs[1]+=pow(ending-start,2);
                time_count[1]+=1;
                start=ending;
                //
                // only assigns values once
                if (ij==jk){
                    if (ij==0){
                        Rls1(j,0) = Rs1;
                        Lls1(j,0) = Ld.col(0).sum();
                    }
                    Rls2(j,ij) = Rs2;
                    Lls2(j,ij) = Ld.col(1).sum();
                }
                Rls3(j,ijk) = Rs3;
                Lls3(j,ijk) = Ld.col(3).sum();
                //
                end_point = system_clock::now();
                ending = time_point_cast<microseconds>(end_point).time_since_epoch().count();
                time_ref[2]+=ending-start;
                time_refs[2]+=pow(ending-start,2);
                time_count[2]+=1;
                start=ending;
                //
            }
        }
        cout << "df206 ";
        for (int ijk=0;ijk<time_ref.size();ijk++){
            cout << (time_ref[ijk]/time_count[ijk])*1e-6 << " ";
        }
        cout << " " << endl;
        cout << "df207 ";
        for (int ijk=0;ijk<time_ref.size();ijk++){
            cout <<  sqrt(time_refs[ijk]/time_count[ijk] - pow(time_ref[ijk]/time_count[ijk],2))*1e-6 << " ";
        }
        cout << " " << endl;
    } else {
        #pragma omp parallel for schedule(dynamic) num_threads(nthreads) collapse(2)
        for (int ijk=0;ijk<totalnum*(totalnum+1)/2;ijk++){//totalnum*(totalnum+1)/2
            for (int j=0;j<ntime;j++){
                int ij = 0;
                int jk = ijk;
                while (jk>ij){
                    ij++;
                    jk-=ij;
                }
                double Rs1 = 0;
                double Rs2 = 0;
                double Rs2t = 0;
                double Rs3 = 0;
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
                //Now has the grouping pairs
                int dj = RiskFail(j,1)-RiskFail(j,0)+1;
                for (int i = 0; i < InGroup.size()-1; i=i+2){
                    Rs1 += R.block(InGroup[i]-1,0,InGroup[i+1]-InGroup[i]+1,1).sum();
                    Rs2 += Rd.block(InGroup[i]-1,ij,InGroup[i+1]-InGroup[i]+1,1).sum();
                    Rs2t += Rd.block(InGroup[i]-1,jk,InGroup[i+1]-InGroup[i]+1,1).sum();
                    Rs3 += Rdd.block(InGroup[i]-1,ijk,InGroup[i+1]-InGroup[i]+1,1).sum();
                } //precalculates the sums of risk groups
                MatrixXd Ld = MatrixXd::Zero(dj,4);
                Ld << R.block(RiskFail(j,0),0,dj,1), Rd.block(RiskFail(j,0),ij,dj,1), Rd.block(RiskFail(j,0),jk,dj,1) ,Rdd.block(RiskFail(j,0),ijk,dj,1);//sum of risks in group
                // only assigns values once
                if (ij==jk){
                    if (ij==0){
                        Rls1(j,0) = Rs1;
                        Lls1(j,0) = Ld.col(0).sum();
                    }
                    Rls2(j,ij) = Rs2;
                    Lls2(j,ij) = Ld.col(1).sum();
                }
                Rls3(j,ijk) = Rs3;
                Lls3(j,ijk) = Ld.col(3).sum();
            }
        }
    }
    return;
}


void Calc_LogLik(const int& nthreads,const IntegerMatrix& RiskFail, const vector<string>&  RiskGroup, const int& totalnum, const int& ntime, const MatrixXd& R, const MatrixXd& Rd, const MatrixXd& Rdd, const MatrixXd& RdR, const MatrixXd& RddR,const MatrixXd& Rls1,const MatrixXd& Rls2,const MatrixXd& Rls3,const MatrixXd& Lls1,const MatrixXd& Lls2,const MatrixXd& Lls3, vector<double>& Ll, vector<double>& Lld, vector<double>& Lldd, bool debugging,string ties_method){
    //
    //Calc_LogLik( nthreads, RiskFail, RiskGroup, totalnum, ntime, R, Rd, Rdd,RdR,RddR, Rls1, Rls2, Rls3, Lls1, Lls2, Lls3, Ll, Lld, Lldd, debugging);
    //
    if (debugging){
        //
        cout << "Starting Debug" << endl;
        vector<double> time_ref(4,0.0);
        vector<double> time_refs(4,0.0);
        vector<int> time_count(4,0);
        cout << time_ref.size() << " " << time_refs.size() << endl;
        #pragma omp declare reduction(vec_int_plus : std::vector<int> : \
            std::transform(omp_out.begin(), omp_out.end(), omp_in.begin(), omp_out.begin(), std::plus<int>())) \
            initializer(omp_priv = omp_orig)
        //
        #pragma omp declare reduction(vec_double_plus : std::vector<double> : \
            std::transform(omp_out.begin(), omp_out.end(), omp_in.begin(), omp_out.begin(), std::plus<double>())) \
            initializer(omp_priv = omp_orig)
        #pragma omp parallel for schedule(dynamic) num_threads(nthreads) reduction(vec_double_plus:Ll,Lld,Lldd, time_ref,time_refs) reduction(vec_int_plus:time_count) collapse(2)
        for (int ijk=0;ijk<totalnum*(totalnum+1)/2;ijk++){//performs log-likelihood calculations for every derivative combination and risk group
            for (int j=0;j<ntime;j++){
                //
                time_point<system_clock> start_point, end_point;
                start_point = system_clock::now();
                auto start = time_point_cast<microseconds>(start_point).time_since_epoch().count();
                end_point = system_clock::now();
                auto ending = time_point_cast<microseconds>(end_point).time_since_epoch().count(); //The time duration is tracked
                //
                int ij = 0;
                int jk = ijk;
                while (jk>ij){ //Splits into two indices
                    ij++;
                    jk-=ij;
                }
                double Rs1 = Rls1(j,0);
                double Rs2 = Rls2(j,ij);
                double Rs2t = Rls2(j,jk);
                double Rs3 = Rls3(j,ijk);
                //
                int dj = RiskFail(j,1)-RiskFail(j,0)+1;
                MatrixXd Ld = MatrixXd::Zero(dj,4);
                Ld << R.block(RiskFail(j,0),0,dj,1), RdR.block(RiskFail(j,0),ij,dj,1), RdR.block(RiskFail(j,0),jk,dj,1) ,RddR.block(RiskFail(j,0),ijk,dj,1);//rows with events
                //
                MatrixXd Ldm = MatrixXd::Zero(dj,4);
                Vector4d Ldcs;
                if (ties_method=="efron"){
                    Ldcs << Lls1(j,0), Lls2(j,ij), Lls2(j,jk), Lls3(j,ijk);
                    for (int i = 0; i < dj; i++){ //adds in the efron approximation terms
                        Ldm.row(i) = (-double(i) / double(dj)) *Ldcs.array();
                    }
                }
                Ldm.col(0) = Ldm.col(0).array() + Rs1;
                Ldm.col(1) = Ldm.col(1).array() + Rs2;
                Ldm.col(2) = Ldm.col(2).array() + Rs2t;
                Ldm.col(3) = Ldm.col(3).array() + Rs3;
                //
                end_point = system_clock::now();
                ending = time_point_cast<microseconds>(end_point).time_since_epoch().count();
                time_ref[0]+=ending-start;
                time_refs[0]+=pow(ending-start,2);
                time_count[0]+=1;
                start=ending;
                //
                // Calculates the left-hand side terms
                MatrixXd temp1 = MatrixXd::Zero(Ld.rows(),1);
                MatrixXd temp2 = MatrixXd::Zero(Ld.rows(),1);
                temp1 = Ld.col(0).array().log();
                double Ld1 =  (temp1.array().isFinite()).select(temp1,0).sum();
                temp1 = Ld.col(1).array();
                temp2 = Ld.col(2).array();
                double Ld2 = (temp1.array().isFinite()).select(temp1,0).sum();
                temp1 = Ld.col(3).array() - (temp1.array() * temp2.array());
                double Ld3 = (temp1.array().isFinite()).select(temp1,0).sum();
                //
                end_point = system_clock::now();
                ending = time_point_cast<microseconds>(end_point).time_since_epoch().count();
                time_ref[1]+=ending-start;
                time_refs[1]+=pow(ending-start,2);
                time_count[1]+=1;
                start=ending;
                //
                // calculates the right-hand side terms
                temp1 = Ldm.col(0).array().log();
                Rs1 =  (temp1.array().isFinite()).select(temp1,0).sum();
                temp1 = Ldm.col(1).array() * (Ldm.col(0).array().pow(-1).array());
                temp2 = Ldm.col(2).array() * (Ldm.col(0).array().pow(-1).array());
                Rs2 = (temp1.array().isFinite()).select(temp1,0).sum();
                temp1 = Ldm.col(3).array() * (Ldm.col(0).array().pow(-1).array()) - temp1.array() * temp2.array();
                Rs3 = (temp1.array().isFinite()).select(temp1,0).sum();
                //
                //
                end_point = system_clock::now();
                ending = time_point_cast<microseconds>(end_point).time_since_epoch().count();
                time_ref[2]+=ending-start;
                time_refs[2]+=pow(ending-start,2);
                time_count[2]+=1;
                start=ending;
                //
                if (ij==jk){
                    Ll[ij] += Ld1 - Rs1;
                    Lld[ij] += Ld2 - Rs2;
                }
                Lldd[ij*totalnum+jk] += Ld3 - Rs3; //sums the log-likelihood and derivatives
                //
                end_point = system_clock::now();
                ending = time_point_cast<microseconds>(end_point).time_since_epoch().count();
                time_ref[3]+=ending-start;
                time_refs[3]+=pow(ending-start,2);
                time_count[3]+=1;
                start=ending;
                //
            }
        }
        cout << "df208 ";
        for (int ijk=0;ijk<time_ref.size();ijk++){
            cout << (time_ref[ijk]/time_count[ijk])*1e-6 << " ";
        }
        cout << " " << endl;
        cout << "df209 ";
        for (int ijk=0;ijk<time_ref.size();ijk++){
            cout <<  sqrt(time_refs[ijk]/time_count[ijk] - pow(time_ref[ijk]/time_count[ijk],2))*1e-6 << " ";
        }
        cout << " " << endl;
    } else {
        #pragma omp declare reduction(vec_double_plus : std::vector<double> : \
            std::transform(omp_out.begin(), omp_out.end(), omp_in.begin(), omp_out.begin(), std::plus<double>())) \
            initializer(omp_priv = omp_orig)
        #pragma omp parallel for schedule(dynamic) num_threads(nthreads) reduction(vec_double_plus:Ll,Lld,Lldd) collapse(2)
        for (int ijk=0;ijk<totalnum*(totalnum+1)/2;ijk++){//performs log-likelihood calculations for every derivative combination and risk group
            for (int j=0;j<ntime;j++){
                int ij = 0;
                int jk = ijk;
                while (jk>ij){
                    ij++;
                    jk-=ij;
                }
                double Rs1 = Rls1(j,0);
                double Rs2 = Rls2(j,ij);
                double Rs2t = Rls2(j,jk);
                double Rs3 = Rls3(j,ijk);
                //
                int dj = RiskFail(j,1)-RiskFail(j,0)+1;
                MatrixXd Ld = MatrixXd::Zero(dj,4);
                Ld << R.block(RiskFail(j,0),0,dj,1), RdR.block(RiskFail(j,0),ij,dj,1), RdR.block(RiskFail(j,0),jk,dj,1) ,RddR.block(RiskFail(j,0),ijk,dj,1);//rows with events
                //
                MatrixXd Ldm = MatrixXd::Zero(dj,4);
                Vector4d Ldcs;
                if (ties_method=="efron"){
                    Ldcs << Lls1(j,0), Lls2(j,ij), Lls2(j,jk), Lls3(j,ijk);
                    for (int i = 0; i < dj; i++){ //adds in the efron approximation terms
                        Ldm.row(i) = (-double(i) / double(dj)) *Ldcs.array();
                    }
                }
                Ldm.col(0) = Ldm.col(0).array() + Rs1;
                Ldm.col(1) = Ldm.col(1).array() + Rs2;
                Ldm.col(2) = Ldm.col(2).array() + Rs2t;
                Ldm.col(3) = Ldm.col(3).array() + Rs3;
                // Calculates the left-hand side terms
                MatrixXd temp1 = MatrixXd::Zero(Ld.rows(),1);
                MatrixXd temp2 = MatrixXd::Zero(Ld.rows(),1);
                temp1 = Ld.col(0).array().log();
                double Ld1 =  (temp1.array().isFinite()).select(temp1,0).sum();
                temp1 = Ld.col(1).array();
                temp2 = Ld.col(2).array();
                double Ld2 = (temp1.array().isFinite()).select(temp1,0).sum();
                temp1 = Ld.col(3).array() - (temp1.array() * temp2.array());
                double Ld3 = (temp1.array().isFinite()).select(temp1,0).sum();
                // calculates the right-hand side terms
                temp1 = Ldm.col(0).array().log();
                Rs1 =  (temp1.array().isFinite()).select(temp1,0).sum();
                temp1 = Ldm.col(1).array() * (Ldm.col(0).array().pow(-1).array());
                temp2 = Ldm.col(2).array() * (Ldm.col(0).array().pow(-1).array());
                Rs2 = (temp1.array().isFinite()).select(temp1,0).sum();
                temp1 = Ldm.col(3).array() * (Ldm.col(0).array().pow(-1).array()) - temp1.array() * temp2.array();
                Rs3 = (temp1.array().isFinite()).select(temp1,0).sum();
                //
                if (ij==jk){
                    Ll[ij] += Ld1 - Rs1;
                    Lld[ij] += Ld2 - Rs2;
                }
                Lldd[ij*totalnum+jk] += Ld3 - Rs3; //sums the log-likelihood and derivatives
            }
        }
        
    }
    #pragma omp parallel for num_threads(nthreads)
    for (int ijk=0;ijk<totalnum*(totalnum+1)/2;ijk++){//fills second-derivative matrix
        int ij = 0;
        int jk = ijk;
        while (jk>ij){
            ij++;
            jk-=ij;
        }
        Lldd[jk*totalnum+ij] = Lldd[ij*totalnum+jk];
    }
    return;
}

void Amfit_LogLik(const int& nthreads, const int& totalnum, const MatrixXd& PyrC, const MatrixXd& R, const MatrixXd& Rd, const MatrixXd& Rdd, const MatrixXd& RdR, const MatrixXd& RddR, vector<double>& Ll, vector<double>& Lld, vector<double>& Lldd, bool debugging){
    //
    // Amfit_LogLik( nthreads, totalnum, PyrC, R, Rd, Rdd, RdR, RddR, Ll, Lld, Lldd, debugging)
    //
    if (debugging){
        MatrixXd temp(Rd.rows(),Rd.cols());
        VectorXd CoL=VectorXd::Zero(Rd.rows());
        
        temp = (PyrC.col(1).array() * (PyrC.col(0).array() * R.col(0).array()).array().log()).array() - (PyrC.col(0).array() * R.col(0).array());
        fill(Ll.begin(), Ll.end(), (temp.array().isFinite()).select(temp,0).sum());
        
        CoL = PyrC.col(1).array() * R.col(0).array().pow(-1).array();

        #pragma omp parallel for schedule(dynamic) num_threads(nthreads)
        for (int ijk=0;ijk<totalnum*(totalnum+1)/2;ijk++){//totalnum*(totalnum+1)/2
            int ij = 0;
            int jk = ijk;
            while (jk>ij){
                ij++;
                jk-=ij;
            }
            VectorXd temp(Rdd.rows(),1);
            temp = Rdd.col(ijk).array() * ( CoL.array() - PyrC.col(0).array()) - PyrC.col(1).array() * RdR.col(ij).array() * RdR.col(jk).array();
            Lldd[ij*totalnum+jk] = (temp.array().isFinite()).select(temp,0).sum();
            if (ij!=jk){
                Lldd[jk*totalnum+ij] = (temp.array().isFinite()).select(temp,0).sum();
            } else{
                temp = Rd.col(ij).array() * ( CoL.array() - PyrC.col(0).array());
                Lld[ij] = (temp.array().isFinite()).select(temp,0).sum();
            }
        }
    } else {
        MatrixXd temp(Rd.rows(),Rd.cols());
        VectorXd CoL=VectorXd::Zero(Rd.rows());
        
        temp = (PyrC.col(1).array() * (PyrC.col(0).array() * R.col(0).array()).array().log()).array() - (PyrC.col(0).array() * R.col(0).array());
        fill(Ll.begin(), Ll.end(), (temp.array().isFinite()).select(temp,0).sum());
        
        CoL = PyrC.col(1).array() * R.col(0).array().pow(-1).array();

        #pragma omp parallel for schedule(dynamic) num_threads(nthreads)
        for (int ijk=0;ijk<totalnum*(totalnum+1)/2;ijk++){//totalnum*(totalnum+1)/2
            int ij = 0;
            int jk = ijk;
            while (jk>ij){
                ij++;
                jk-=ij;
            }
            VectorXd temp(Rdd.rows(),1);
            temp = Rdd.col(ijk).array() * ( CoL.array() - PyrC.col(0).array()) - PyrC.col(1).array() * RdR.col(ij).array() * RdR.col(jk).array();
            Lldd[ij*totalnum+jk] = (temp.array().isFinite()).select(temp,0).sum();
            if (ij!=jk){
                Lldd[jk*totalnum+ij] = (temp.array().isFinite()).select(temp,0).sum();
            } else{
                temp = Rd.col(ij).array() * ( CoL.array() - PyrC.col(0).array());
                Lld[ij] = (temp.array().isFinite()).select(temp,0).sum();
            }
        }
    }
    return;
}

void Calc_Change(const int& nthreads, const int& totalnum,const int& fir, const int& der_iden, const double& dbeta_cap, const double& dose_abs_max, const double& lr, const double& abs_max, const vector<double>& Ll, const vector<double>& Lld, const vector<double>& Lldd, vector<double>& dbeta, const bool change_all,const StringVector&   tform, const double& dint, IntegerVector KeepConstant, bool debugging){
    //
    //Calc_Change( nthreads, totalnum, fir, der_iden, dbeta_cap, dose_abs_max, lr, abs_max, Ll, Lld, Lldd, dbeta, change_all, tform, dint, KeepConstant, debugging);
    //
    if (debugging){
        #pragma omp parallel for schedule(dynamic) num_threads(nthreads)
        for (int ijk=0;ijk<totalnum;ijk++){
            if (change_all){
                if (KeepConstant[ijk]==0){
                    dbeta[ijk] = -lr * Lld[ijk] / Lldd[ijk*totalnum+ijk];
                    //
                    double dbeta_max = abs(Ll[ijk]/Lld[ijk] * dbeta_cap);//uses newtonian step for zero log-likelihood as a limit
                    if (abs(dbeta[ijk])>dbeta_max){
                        dbeta[ijk] = dbeta_max * sign(dbeta[ijk]);
                    }
                    if ((tform[ijk]=="step_int")||(tform[ijk]=="lin_int")){ //the threshold values use different maximum deviation values
                        if (abs(dbeta[ijk])>dose_abs_max){
                            dbeta[ijk] = dose_abs_max * sign(dbeta[ijk]);
                        }
                    }else{
                        if (abs(dbeta[ijk])>abs_max){
                            dbeta[ijk] = abs_max * sign(dbeta[ijk]);
                        }
                    }
                } else {
                    dbeta[ijk]=0;
                }
            }else{
                if (ijk!=der_iden){//Validation requires controlled changes
                    dbeta[ijk] = 0.0;
                } else {
                    if ((tform[ijk]=="step_int")||(tform[ijk]=="lin_int")){
                        dbeta[ijk] = dint;
                    } else {
                        dbeta[ijk] = 0.001;
                    }
                }
            }
        }
    } else {
        #pragma omp parallel for schedule(dynamic) num_threads(nthreads)
        for (int ijk=0;ijk<totalnum;ijk++){
            if (change_all){
                if (KeepConstant[ijk]==0){
                    dbeta[ijk] = -lr * Lld[ijk] / Lldd[ijk*totalnum+ijk];
                    //
                    double dbeta_max = abs(Ll[ijk]/Lld[ijk] * dbeta_cap);//uses newtonian step for zero log-likelihood as a limit
                    if (abs(dbeta[ijk])>dbeta_max){
                        dbeta[ijk] = dbeta_max * sign(dbeta[ijk]);
                    }
                    if ((tform[ijk]=="step_int")||(tform[ijk]=="lin_int")){ //the threshold values use different maximum deviation values
                        if (abs(dbeta[ijk])>dose_abs_max){
                            dbeta[ijk] = dose_abs_max * sign(dbeta[ijk]);
                        }
                    }else{
                        if (abs(dbeta[ijk])>abs_max){
                            dbeta[ijk] = abs_max * sign(dbeta[ijk]);
                        }
                    }
                } else {
                    dbeta[ijk]=0;
                }
            }else{
                if (ijk!=der_iden){//Validation requires controlled changes
                    dbeta[ijk] = 0.0;
                } else {
                    if ((tform[ijk]=="step_int")||(tform[ijk]=="lin_int")){
                        dbeta[ijk] = dint;
                    } else {
                        dbeta[ijk] = 0.01;
                    }
                }
            }
        }
    }
    return;
}


///
/// ------------------------------------------------------------------------------------- ///
///

// [[Rcpp::export]]
List peanut_null(int ntime, List Control, NumericMatrix df_groups, NumericVector tu){
    //----------------------------------------------------------------------------------------------------------------//
    //
    // Calculates null model
    //
    // Converts from Rcpp types to efficient Eigen types
    bool verbose = Control["verbose"];
    string ties_method =Control["ties"];
    //
    //----------------------------------------------------------------------------------------------------------------//
    List res = LogLik_PEANUT_null(ntime, df_groups, tu, verbose, ties_method);
    //----------------------------------------------------------------------------------------------------------------//
    return res;
}

List LogLik_PEANUT_null( int ntime, NumericMatrix df_groups, NumericVector tu, bool verbose, string ties_method){
    srand (time(NULL));
    //
    // null model value calculation
    // --------------------------------------------------------------------------- //
    // The same steps are taken for the non-null case, with the exception of derivative calculations
    // --------------------------------------------------------------------------- //
    //
    if (verbose){
        cout << "START_NEW" << endl;
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
        cout << ctime(&gibtime) << endl;
    }
    //
    cout.precision(10); //forces higher precision numbers printed to terminal
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
        cout<<"df99,"<<(ending-start)<<",Prep_Terms"<<endl;
        gibtime = system_clock::to_time_t(system_clock::now());
        cout << ctime(&gibtime) << endl;
    }
    //
    R = R.array() + 1.0;
    //
    if (verbose){
        cout << "risk checked ";
        for (int ijk=0;ijk<totalnum;ijk++){
            cout << R.col(0).sum() << " ";
        }
        cout << " " << endl;
        //
        end_point = system_clock::now();
        ending = time_point_cast<microseconds>(end_point).time_since_epoch().count();
        cout<<"df100 "<<(ending-start)<<" "<<0<<" "<<0<<" "<<-1<<",Prep_R"<<endl;
        gibtime = system_clock::to_time_t(system_clock::now());
        cout << ctime(&gibtime) << endl;
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
        cout<<"df100 "<<(ending-start)<<" "<<0<<" "<<0<<" "<<-1<<",Prep_List"<<endl;
        gibtime = system_clock::to_time_t(system_clock::now());
        cout << ctime(&gibtime) << endl;
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
        cout << "riskr checked ";
        for (int ijk=0;ijk<totalnum;ijk++){
            cout << Rls1.col(0).sum() << " ";
        }
        cout << " " << endl;
        //
        cout << "riskl checked ";
        for (int ijk=0;ijk<totalnum;ijk++){
            cout << Lls1.col(0).sum() << " ";
        }
        cout << " " << endl;
    }
    //
    //
    Calc_Null_LogLik( nthreads, RiskFail, RiskGroup, ntime, R, Rls1, Lls1, Ll, ties_method);
    //
    List res_list = List::create(_["LogLik"]=wrap(Ll),_["AIC"]=-2*Ll[0]);
    // returns a list of results
    return res_list;
}

void Calculate_Null_Sides(const IntegerMatrix& RiskFail, const vector<string>&  RiskGroup, const int& ntime, const MatrixXd& R, MatrixXd& Rls1, MatrixXd& Lls1,const int& nthreads){
    //
    //Calculate_Sides( RiskFail, RiskGroup, totalnum, ntime, R, Rd, Rdd, Rls1, Rls2, Rls3, Lls1, Lls2, Lls3,nthreads, debugging);
    //
    #pragma omp parallel for schedule(dynamic) num_threads(nthreads)
    for (int j=0;j<ntime;j++){
        double Rs1 = 0;
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
        //Now has the grouping pairs
        int dj = RiskFail(j,1)-RiskFail(j,0)+1;
        for (int i = 0; i < InGroup.size()-1; i=i+2){
            Rs1 += R.block(InGroup[i]-1,0,InGroup[i+1]-InGroup[i]+1,1).sum();
        } //precalculates the sums of risk groups
        MatrixXd Ld = MatrixXd::Zero(dj,1);
        Ld << R.block(RiskFail(j,0),0,dj,1);//sum of risks in group
        // only assigns values once
        Rls1(j,0) = Rs1;
        Lls1(j,0) = Ld.col(0).sum();
    }
    return;
}


void Calc_Null_LogLik(const int& nthreads,const IntegerMatrix& RiskFail, const vector<string>&  RiskGroup, const int& ntime, const MatrixXd& R, const MatrixXd& Rls1,const MatrixXd& Lls1, vector<double>& Ll, string ties_method){
    //
    //Calc_LogLik( nthreads, RiskFail, RiskGroup, totalnum, ntime, R, Rd, Rdd,RdR,RddR, Rls1, Rls2, Rls3, Lls1, Lls2, Lls3, Ll, Lld, Lldd, debugging);
    //
    #pragma omp declare reduction(vec_double_plus : std::vector<double> : \
        std::transform(omp_out.begin(), omp_out.end(), omp_in.begin(), omp_out.begin(), std::plus<double>())) \
        initializer(omp_priv = omp_orig)
    #pragma omp parallel for schedule(dynamic) num_threads(nthreads) reduction(vec_double_plus:Ll)
    for (int j=0;j<ntime;j++){
        double Rs1 = Rls1(j,0);
        //
        int dj = RiskFail(j,1)-RiskFail(j,0)+1;
        MatrixXd Ld = MatrixXd::Zero(dj,1);
        Ld << R.block(RiskFail(j,0),0,dj,1);//rows with events
        //
        MatrixXd Ldm = MatrixXd::Zero(dj,1);
        Vector4d Ldcs;
        if (ties_method=="efron"){
            Ldcs << Lls1(j,0);
            for (int i = 0; i < dj; i++){ //adds in the efron approximation terms
                Ldm.row(i) = (-double(i) / double(dj)) *Ldcs.array();
            }
        }
        Ldm.col(0) = Ldm.col(0).array() + Rs1;
        // Calculates the left-hand side terms
        MatrixXd temp1 = MatrixXd::Zero(Ld.rows(),1);
        temp1 = Ld.col(0).array().log();
        double Ld1 =  (temp1.array().isFinite()).select(temp1,0).sum();
        // calculates the right-hand side terms
        temp1 = Ldm.col(0).array().log();
        Rs1 =  (temp1.array().isFinite()).select(temp1,0).sum();
        //
        Ll[0] += Ld1 - Rs1;
    }
    return;
}

