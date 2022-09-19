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
List peanut_transition(IntegerVector Term_n, StringVector tform, NumericVector a_n,IntegerVector dfc,NumericMatrix x_all, int fir, int der_iden,string modelform, List Control, NumericMatrix df_groups, NumericVector tu, IntegerVector KeepConstant, int term_tot){
    // Term_n,STerm_n,tform,a_n,dfc,x_all, fir, modelform,length(tu), control,as.matrix(df[,..ce]),tu,keep_constant,term_tot
    //----------------------------------------------------------------------------------------------------------------//
    //
//    const Map<MatrixXd> test_df(as<Map<MatrixXd> >(x_lin));
//    const SparseMatrix<double> test_sp = test_df.sparseView();
    //
    // Converts from Rcpp types to efficient Eigen types
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
    //
//    List res = List::create(_["a"]=Term_n);
    //----------------------------------------------------------------------------------------------------------------//
    List res = LogLik_PEANUT(Term_n, tform, a_n, x_all, dfc,fir, der_iden,modelform, lr, maxiter, halfmax, epsilon, dbeta_cap, abs_max,dose_abs_max, deriv_epsilon, df_groups, tu, change_all,verbose, debugging, KeepConstant, term_tot);
    //----------------------------------------------------------------------------------------------------------------//
    cout << "it got out more" << endl;
    return res;
}


List LogLik_PEANUT( IntegerVector Term_n, StringVector tform, NumericVector a_n,NumericMatrix x_all,IntegerVector dfc,int fir, int der_iden,string modelform, double lr, int maxiter, int halfmax, double epsilon, double dbeta_cap, double abs_max,double dose_abs_max, double deriv_epsilon, NumericMatrix df_groups, NumericVector tu ,bool change_all, bool verbose, bool debugging, IntegerVector KeepConstant, int term_tot){
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
    cout.precision(7); //forces higher precision numbers printed to terminal
    int nthreads = Eigen::nbThreads()-1; //stores how many threads are allocated
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
    MatrixXd Te = MatrixXd::Zero(df0.rows(), 1); //preallocates matrix for non-Derivative column terms
    MatrixXd R = MatrixXd::Zero(df0.rows(), 1); //preallocates matrix for Risks
    MatrixXd Rd = MatrixXd::Zero(df0.rows(), totalnum); //preallocates matrix for Risk derivatives
    MatrixXd Rdd = MatrixXd::Zero(df0.rows(), totalnum*(totalnum+1)/2); //preallocates matrix for Risk second derivatives
    //
    MatrixXd Dose = MatrixXd::Zero(df0.rows(),term_tot); //Matrix of the total dose term values
    MatrixXd nonDose = MatrixXd::Constant(df0.rows(),term_tot,1.0); //Matrix of the total dose term values
    MatrixXd nonDose_LIN = MatrixXd::Constant(df0.rows(),term_tot,0.0);
    MatrixXd nonDose_PLIN = MatrixXd::Constant(df0.rows(),term_tot,1.0);
    MatrixXd nonDose_LOGLIN = MatrixXd::Constant(df0.rows(),term_tot,1.0);
    double dint = dose_abs_max; //The amount of change used to calculate derivatives in threshold paramters
    int total_dose=0; //used later on for a section summing the dose terms
    //
    // totalnum,& Term_n,  tform, dfc,& fir,& T0,& Td0,& Tdd0,& Dose,& nonDose,& beta_0,& df0, dint, nthreads,  debugging
    cout << "starting subterms " << term_tot << endl;
    Make_Subterms( totalnum, Term_n, tform, dfc, fir, T0, Td0, Tdd0, Dose, nonDose, nonDose_LIN, nonDose_PLIN, nonDose_LOGLIN ,beta_0, df0,dint,nthreads, debugging);
    // ---------------------------------------------------------
    // Prints off a series of calaculations to check at what point values are changing
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
    Make_Risks(modelform, tform, Term_n, totalnum, fir, T0, Td0, Tdd0, Te, R, Rd, Rdd, Dose, nonDose, nonDose_LIN, nonDose_PLIN, nonDose_LOGLIN, RdR, RddR, nthreads, debugging);
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
        cout << "A non-positive risk was detected: " << R.minCoeff() << endl;
        return temp_list;
    }
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
//    const static IOFormat CSVFormat(FullPrecision, DontAlignCols, ", ", "\n");
//    ofstream file("Error_Rows.csv");
    Make_Groups( ntime, df_m, RiskFail, RiskGroup, tu, nthreads, debugging);
//    file.close();
//    return temp_list;
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
    MatrixXd Rls3 =MatrixXd::Zero(ntime, totalnum*(totalnum+1)/2);
    MatrixXd Lls1 =MatrixXd::Zero(ntime, 1);
    MatrixXd Lls2 =MatrixXd::Zero(ntime, totalnum);
    MatrixXd Lls3 =MatrixXd::Zero(ntime, totalnum*(totalnum+1)/2);
    //The log-likelihood is calculated in parallel over the risk groups
    vector<double> Ll(totalnum,0.0);
    vector<double> Lld(totalnum,0.0);
    vector<double> Lldd(pow(totalnum,2),0.0);//The second derivative matrix has room for every combination, but only the lower triangle is calculated
    //
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
    Calc_LogLik( nthreads, RiskFail, RiskGroup, totalnum, ntime, R, Rd, Rdd,RdR,RddR, Rls1, Rls2, Rls3, Lls1, Lls2, Lls3, Ll, Lld, Lldd, debugging);
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
    }
//    cout << "df105 ";
//    for (int ij=0;ij<totalnum;ij++){//prints the newton step value for zero derivative
//        cout << Lld[ij]/Lldd[ij*totalnum+ij] << " ";
//    }
//    cout << " " << endl;
//    cout << "df106 ";
//    for (int ij=0;ij<totalnum;ij++){//prints the newton step value for zero log-likelihood
//        cout << Ll[ij]/Lld[ij] << " ";
//    }
//    cout << " " << endl;
//    cout << "df107 " << abs_max << " " << dose_abs_max << " " << Ll_comp[0] << " " << Ll_comp[1] << endl;//prints several convergence terms
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
            //Make_Subterms( totalnum, Term_n, tform, dfc, fir, T0, Td0, Tdd0, Dose, nonDose,beta_0, df0,dint,nthreads, debugging);
            Make_Subterms( totalnum, Term_n, tform, dfc, fir, T0, Td0, Tdd0, Dose,  nonDose, nonDose_LIN, nonDose_PLIN, nonDose_LOGLIN ,beta_0, df0,dint,nthreads, debugging);
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
            }
            //
            RdR = MatrixXd::Zero(df0.rows(), totalnum);
            RddR = MatrixXd::Zero(df0.rows(), totalnum*(totalnum+1)/2);
            //
            //
            Make_Risks(modelform, tform, Term_n, totalnum, fir, T0, Td0, Tdd0, Te, R, Rd, Rdd, Dose, nonDose, nonDose_LIN, nonDose_PLIN, nonDose_LOGLIN, RdR, RddR, nthreads, debugging);
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
            Calc_LogLik( nthreads, RiskFail, RiskGroup, totalnum, ntime, R, Rd, Rdd,RdR,RddR, Rls1, Rls2, Rls3, Lls1, Lls2, Lls3, Ll, Lld, Lldd, debugging);
            
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
            }
    //            cout << "df105 ";
    //            for (int ij=0;ij<totalnum;ij++){
    //                cout << Lld[ij]/Lldd[ij*totalnum+ij] << " ";
    //            }
    //            cout << " " << endl;
    //            cout << "df106 ";
    //            for (int ij=0;ij<totalnum;ij++){
    //                cout << Ll[ij]/Lld[ij] << " ";
    //            }
    //            cout << " " << endl;
    //            cout << "df107 " << abs_max << " " << dose_abs_max << " " << Ll_comp[0] << " " << Ll_comp[1] << endl;
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
            Make_Subterms( totalnum, Term_n, tform, dfc, fir, T0, Td0, Tdd0, Dose, nonDose, nonDose_LIN, nonDose_PLIN, nonDose_LOGLIN,beta_0, df0,dint,nthreads, debugging);;
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
            Make_Risks(modelform, tform, Term_n, totalnum, fir, T0, Td0, Tdd0, Te, R, Rd, Rdd, Dose, nonDose, nonDose_LIN, nonDose_PLIN, nonDose_LOGLIN, RdR, RddR, nthreads, debugging);
            R = (R.array().isFinite()).select(R,0);
            Rd = (Rd.array().isFinite()).select(Rd,0);
            Rdd = (Rdd.array().isFinite()).select(Rdd,0);
            RdR = (RdR.array().isFinite()).select(RdR,0);
            RddR = (RddR.array().isFinite()).select(RddR,0);
            //
            temp_list = List::create(_["betas"]=wrap(beta_0),_["Status"]="FAILED"); //used as a dummy return value for code checking
            if (R.minCoeff()<=0){
                cout << "A non-positive risk was detected" << endl;
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
            Calc_LogLik( nthreads, RiskFail, RiskGroup, totalnum, ntime, R, Rd, Rdd,RdR,RddR, Rls1, Rls2, Rls3, Lls1, Lls2, Lls3, Ll, Lld, Lldd, debugging);
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
    Calc_LogLik( nthreads, RiskFail, RiskGroup, totalnum, ntime, R, Rd, Rdd,RdR,RddR, Rls1, Rls2, Rls3, Lls1, Lls2, Lls3, Ll, Lld, Lldd, debugging);
    // Changes the parameter back into the original form
//    #pragma omp parallel for schedule(dynamic) num_threads(nthreads)
//    for (int ij=0;ij<totalnum;ij++){
//        int ind0 = ij;
//        if (ind0 < dose_num_tot){
//            ;
//        } else {
//            ind0 = ind0-dose_num_tot;
//            if (include_bool[0]==1){
//                if (ind0 < beta_lin.size()){
//                    beta_lin[ind0] = beta_0[ij];
//                    //
//                } else {
//                    //one exists and its not one
//                    ind0 = ind0 - beta_lin.size();
//                    if (include_bool[1]==1){
//                        if (ind0 < beta_loglin.size()){
//                            //one and two exists and is two
//                            beta_loglin[ind0] = beta_0[ij];
//                            //
//                        } else{
//                            //one exists, two does, must be three
//                            if (include_bool[2]!=1){
//                                throw invalid_argument( "Are all three used? 0" );
//                            }
//                            ind0 = ind0 - beta_loglin.size();
//                            beta_plin[ind0] = beta_0[ij];
//                            //
//                        }
//                    } else{
//                        //one exists, and two doesn't exist, must be three
//                        if (include_bool[2]!=1){
//                            throw invalid_argument( "Are all first and third used?" );
//                        }
//                        beta_plin[ind0] = beta_0[ij];
//                        //
//                    }
//                }
//            }else{
//                //one doesn't exist
//                if (include_bool[1]==1){
//                    if (ind0 < beta_loglin.size()){
//                        //one doesn't exist and two exists and is two
//                        beta_loglin[ind0] = beta_0[ij];
//                        //
//                    } else{
//                        //one doesn't exist, two does, must be three
//                        if (include_bool[2]!=1){
//                            throw invalid_argument( "Are all three used? 1" );
//                        }
//                        ind0 = ind0 - beta_loglin.size();
//                        beta_plin[ind0] = beta_0[ij];
//                        //
//                    }
//                } else{
//                    //one doesn't exist, and two doesn't exist, must be three
//                    if (include_bool[2]!=1){
//                        throw invalid_argument( "Are all first and third used?" );
//                    }
//                    //
//                }
//            }
//        }
//    }
//    //
//    if (verbose){
//        cout << "Reassign" << endl;
//    }
//    //
//    // --------------------------------
//    // Return the results
//    // --------------------------------
//    if (include_bool[3]==1){
//        vector<double> beta_loglin_slope;
//        vector<double> beta_loglin_top;
//        vector<double> beta_lin_slope;
//        vector<double> beta_lin_int;
//        vector<double> beta_quad;
//        vector<double> beta_step_slope;
//        vector<double> beta_step_int;
//        int loglin_size=0;
//        int lin_size=0;
//        int quad_size=0;
//        int step_size=0;
//        int dub_off=0;
//        for (int ijk=0;ijk<df_dose.cols();ijk++){
//            beta_loglin_slope = beta_loglin_slopes_CPP[ijk];
//            beta_loglin_top = beta_loglin_tops_CPP[ijk];
//            beta_lin_slope = beta_lin_slopes_CPP[ijk];
//            beta_lin_int = beta_lin_ints_CPP[ijk];
//            beta_quad = beta_quads_CPP[ijk];
//            beta_step_slope = beta_step_slopes_CPP[ijk];
//            beta_step_int = beta_step_ints_CPP[ijk];
//            //
//            if ((beta_loglin_slope.size()==1)&&(beta_loglin_slope[0]==0.0)){
//                ;
//            } else {
//                loglin_size = beta_loglin_slope.size();
//            }
//            if ((beta_lin_slope.size()==1)&&(beta_lin_slope[0]==0.0)){
//                ;
//            } else {
//                lin_size = beta_lin_slope.size();
//            }
//            if ((beta_quad.size()==1)&&(beta_quad[0]==0.0)){
//                ;
//            } else {
//                quad_size = beta_quad.size();
//            }
//            if ((beta_step_slope.size()==1)&&(beta_step_slope[0]==0.0)){
//                ;
//            } else {
//                step_size = beta_step_slope.size();
//            }
//            //
//            total_dose = loglin_size + lin_size + quad_size + step_size;//beta_loglin_slope.size() + beta_lin_slope.size() + beta_quad.size() + beta_step_int.size();
//            dub_off=0;
//            if ((beta_loglin_slope[0]==1)&&(loglin_size==1)){
//                dub_off=1;
//            }
//            #pragma omp parallel for schedule(dynamic) num_threads(nthreads)
//            for (int ij=0; ij<total_dose;ij++){
//                if (ij < loglin_size){
//                    if ((beta_loglin_slope[ij]==1)&&(loglin_size==1)){
//                        int ind0 = cumulative_dose_num[ijk]+ij;
//                        //
//                        beta_loglin_top[ij,0]=beta_0[ind0];
//                    } else {
//                        int ind0 = cumulative_dose_num[ijk]+2*ij;
//                        int ind1 = ind0 + 1; 
//                        //
//                        beta_loglin_slope[ij,0]=beta_0[ind0];
//                        beta_loglin_top[ij,0]=beta_0[ind1];
//                    }
//                    
//                } else if (ij < loglin_size + lin_size){
//                    int jk = ij - loglin_size;
//                    int ind0 = cumulative_dose_num[ijk]+2*loglin_size - dub_off  + 2*jk;
//                    int ind1 = ind0 + 1; 
//                    //
//                    beta_lin_slope[ij,0]=beta_0[ind0];
//                    beta_lin_int[ij,0]=beta_0[ind1];
//                } else if (ij < loglin_size + lin_size + quad_size){
//                    int jk = ij - loglin_size - lin_size;
//                    int ind0 = cumulative_dose_num[ijk]+2*loglin_size - dub_off  + 2*lin_size+jk;
//                    //
//                    beta_quad[ij,0]=beta_0[ind0];
//                } else {
//                    int jk = ij - loglin_size - lin_size - quad_size;
//                    int ind0 = cumulative_dose_num[ijk]+2*loglin_size - dub_off  + 2*lin_size + quad_size + 2*jk;
//                    int ind1 = ind0 + 1;
//                    //
//                    beta_step_slope[ij,0]=beta_0[ind0];
//                    beta_step_int[ij,0]=beta_0[ind1];
//                }
//            }
//            //
//            beta_loglin_slopes[ijk] = wrap(beta_loglin_slope);
//            beta_loglin_tops[ijk] = wrap(beta_loglin_top);
//            beta_lin_slopes[ijk] = wrap(beta_lin_slope);
//            beta_lin_ints[ijk] = wrap(beta_lin_int);
//            beta_quads[ijk] = wrap(beta_quad);
//            beta_step_slopes[ijk] = wrap(beta_step_slope);
//            beta_step_ints[ijk] = wrap(beta_step_int);
//            //
//        }
//    }
    List para_list = List::create(_["Term_n"]=Term_n,_["tforms"]=tform);
    List control_list = List::create(_["Iteration"]=iteration);
    NumericVector Lldd_vec = wrap(Lldd);//creates list of dose parameters
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
/*
*/



void removeRow(Eigen::MatrixXd& matrix, unsigned int rowToRemove){
    unsigned int numRows = matrix.rows()-1;
    unsigned int numCols = matrix.cols();

    if( rowToRemove < numRows )
        matrix.block(rowToRemove,0,numRows-rowToRemove,numCols) = matrix.block(rowToRemove+1,0,numRows-rowToRemove,numCols);

    matrix.conservativeResize(numRows,numCols);
}

void removeColumn(Eigen::MatrixXd& matrix, unsigned int colToRemove){
    unsigned int numRows = matrix.rows();
    unsigned int numCols = matrix.cols()-1;

    if( colToRemove < numCols )
        matrix.block(0,colToRemove,numRows,numCols-colToRemove) = matrix.block(0,colToRemove+1,numRows,numCols-colToRemove);

    matrix.conservativeResize(numRows,numCols);
}


void Make_Subterms(const int& totalnum, const IntegerVector& Term_n,const StringVector&  tform, const IntegerVector& dfc, const int& fir, MatrixXd& T0, MatrixXd& Td0, MatrixXd& Tdd0, MatrixXd& Dose, MatrixXd& nonDose,  MatrixXd& nonDose_LIN, MatrixXd& nonDose_PLIN, MatrixXd& nonDose_LOGLIN,const  VectorXd& beta_0,const  MatrixXd& df0,const double& dint, const int& nthreads, bool debugging){
    //
    //Make_Subterms( totalnum, dose_num_tot, dose_term_tot, dose_breaks, beta_loglin_slopes_CPP, beta_loglin_tops_CPP, beta_lin_slopes_CPP, beta_lin_ints_CPP, beta_quads_CPP, beta_step_slopes_CPP, beta_step_ints_CPP, beta_lin, beta_loglin, beta_plin, df_lin, df_loglin, df_plin, df_dose, De, Dde, Ddde, T0, Td0, Tdd0, Dose,cumulative_dose_num,beta_0, df0,dint,nthreads, tform,include_bool, debugging);
    //
    cout << "real start" << endl;
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
            Td0.col(ij) = temp.array();
            Td0.col(ij+1) = temp1.array() * df0.col(df0_c).array();
            Tdd0.col((ij+1)*(ij+2)/2+ij) = temp.array() * df0.col(df0_c).array();
            Tdd0.col((ij+1)*(ij+2)/2+ij+1) = temp1.array() * df0.col(df0_c).array().square().array();
            Dose.col(tn) = Dose.col(tn).array() + T0.col(ij).array();
            
        } else if (as< string>(tform[ij])=="loglin_top"){
            if (ij==0){
                ArrayXd temp = (beta_0[ij] * df0.col(df0_c)).array().exp();
                T0.col(ij) = temp;
                Td0.col(ij) = temp.array() * df0.col(df0_c).array();
                Tdd0.col(ij * (ij+1)/2 + ij) = temp.array() * df0.col(df0_c).array().square().array();
                Dose.col(tn) = Dose.col(tn).array() + T0.col(ij).array();

            } else if (tform[ij-1]!="loglin_slope"){
                ArrayXd temp = (beta_0[ij] * df0.col(df0_c)).array().exp();
                T0.col(ij) = temp;
                Td0.col(ij) = temp.array() * df0.col(df0_c).array();
                Tdd0.col(ij * (ij+1)/2 + ij) = temp.array() * df0.col(df0_c).array().square().array();
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
            Td0.col(ij) = temp.array();
            Td0.col(ij+1) = beta_0[ij] * (temp1.array()-temp0.array()) / 2/dint;
            //
            Tdd0.col((ij+1)*(ij+2)/2+ij) = (temp1.array()-temp0.array()) / 2/dint;
            Tdd0.col((ij+1)*(ij+2)/2+ij+1) = beta_0[ij] * (temp1.array()-2*temp.array()+temp0.array()) / pow(dint,2);
            Dose.col(tn) = Dose.col(tn).array() + T0.col(ij).array();

        } else if (as< string>(tform[ij])=="quad_slope"){
            ArrayXd temp = df0.col(df0_c).array().square();
            //
            T0.col(ij) = beta_0[ij] * temp.array();
            Td0.col(ij) = temp.array();
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
            Td0.col(ij) = temp.array();
            Td0.col(ij+1) = beta_0[ij] * (temp1.array()-temp0.array()) / 2/dint;
            //
            Tdd0.col((ij+1)*(ij+2)/2+ij) = (temp1.array()-temp0.array()) / 2/dint;
            Tdd0.col((ij+1)*(ij+2)/2+ij+1) = beta_0[ij] * (temp1.array()-2*temp.array()+temp0.array()) / pow(dint,2);
            Dose.col(tn) = Dose.col(tn).array() + T0.col(ij).array();

        } else if (as< string>(tform[ij])=="lin") {
            T0.col(ij) = (df0.col(df0_c).array() * beta_0[ij]).matrix();
            Td0.col(ij) = df0.col(df0_c);
            nonDose_LIN.col(tn) = nonDose_LIN.col(tn).array() + T0.col(ij).array();

        } else if (as< string>(tform[ij])=="loglin") {
            T0.col(ij) = (df0.col(df0_c).array() * beta_0[ij]).matrix();
            T0.col(ij) = T0.col(ij).array().exp();
            Td0.col(ij) = df0.col(df0_c).array() * T0.col(ij).array();
            Tdd0.col((ij)*(ij+1)/2+ij) = df0.col(df0_c).array() * Td0.col(ij).array();
            nonDose_LOGLIN.col(tn) = nonDose_LOGLIN.col(tn).array() * T0.col(ij).array();

        } else if (as< string>(tform[ij])=="plin") {
            T0.col(ij) = (df0.col(df0_c).array() * beta_0[ij]).matrix();
            T0.col(ij) = 1 + T0.col(ij).array();
            Td0.col(ij) = df0.col(df0_c);
            nonDose_PLIN.col(tn) = nonDose_PLIN.col(tn).array() + T0.col(ij).array();

        } else {
            ;
        }
    }
    nonDose_LIN = (nonDose_LIN.array() != 0).select(nonDose_LIN,1.0);
    nonDose_LOGLIN = (nonDose_LOGLIN.array() != 0).select(nonDose_LOGLIN,1.0);
    nonDose_PLIN = (nonDose_PLIN.array() != 0).select(nonDose_PLIN,1.0);
//    cout << "non-dose checked ";
//    for (int ijk=0;ijk<term_tot;ijk++){
//        cout << nonDose.col(ijk).array().sum() << " ";
//    }
//    cout << " " << endl;
//    cout << "LIN_non-dose checked ";
//    for (int ijk=0;ijk<term_tot;ijk++){
//        cout << nonDose_LIN.col(ijk).array().sum() << " ";
//    }
//    cout << " " << endl;
//    cout << "PLIN_non-dose checked ";
//    for (int ijk=0;ijk<term_tot;ijk++){
//        cout << nonDose_PLIN.col(ijk).array().sum() << " ";
//    }
//    cout << " " << endl;
//    cout << "LOGLIN_non-dose checked ";
//    for (int ijk=0;ijk<term_tot;ijk++){
//        cout << nonDose_LOGLIN.col(ijk).array().sum() << " ";
//    }
//    cout << " " << endl;
    for (int ijk=0; ijk<nonDose.cols();ijk++){
        nonDose.col(ijk) = nonDose_LIN.col(ijk).array()  * nonDose_PLIN.col(ijk).array()  * nonDose_LOGLIN.col(ijk).array() ;
    }
//    cout << "non-dose checked ";
//    for (int ijk=0;ijk<term_tot;ijk++){
//        cout << nonDose.col(ijk).array().sum() << " ";
//    }
//    cout << " " << endl;
//    cout << "LIN_non-dose checked ";
//    for (int ijk=0;ijk<term_tot;ijk++){
//        cout << nonDose_LIN.col(ijk).array().sum() << " ";
//    }
//    cout << " " << endl;
//    cout << "PLIN_non-dose checked ";
//    for (int ijk=0;ijk<term_tot;ijk++){
//        cout << nonDose_PLIN.col(ijk).array().sum() << " ";
//    }
//    cout << " " << endl;
//    cout << "LOGLIN_non-dose checked ";
//    for (int ijk=0;ijk<term_tot;ijk++){
//        cout << nonDose_LOGLIN.col(ijk).array().sum() << " ";
//    }
//    cout << " " << endl;
    return;
}


void Make_Risks(string modelform, const StringVector& tform, const IntegerVector& Term_n, const int& totalnum, const int& fir, const MatrixXd& T0, const MatrixXd& Td0, const MatrixXd& Tdd0, MatrixXd& Te, MatrixXd& R, MatrixXd& Rd, MatrixXd& Rdd, MatrixXd& Dose, MatrixXd& nonDose,  MatrixXd& nonDose_LIN, MatrixXd& nonDose_PLIN, MatrixXd& nonDose_LOGLIN, MatrixXd& RdR, MatrixXd& RddR, const int& nthreads, bool debugging){
    //
    //Make_Risks(modelform, tform, Term_n, totalnum, fir, T0, Td0, Tdd0, Te, R, Rd, Rdd, Dose, nonDose, RdR, RddR, nthreads, debugging);
    //
    Dose = (Dose.array() > 0).select(Dose,1.0);
    MatrixXd TTerm=MatrixXd::Zero(Dose.rows(),Dose.cols());
    TTerm << Dose.array() * nonDose.array();
    set<string> Dose_Iden;
    Dose_Iden.insert("loglin_top");
    Dose_Iden.insert("loglin_slope");
    Dose_Iden.insert("lin_slope");
    Dose_Iden.insert( "lin_int");
    Dose_Iden.insert("quad_slope");
    Dose_Iden.insert("step_slope");
    Dose_Iden.insert("step_int");
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
                        Rd.col(ij) =   TTerm.col(tij).array() * Dose.col(tij).array().pow(-1).array() * Td0.col(ij).array();
                        Rdd.col(ijk) = TTerm.col(tij).array() * Dose.col(tij).array().pow(-1).array() * Tdd0.col(ijk).array();
                    } else if (tform[ij]=="lin") {
                        Rd.col(ij) =   TTerm.col(tij).array() * nonDose_LIN.col(tij).array().pow(-1).array() * Td0.col(ij).array();
                        Rdd.col(ijk) = TTerm.col(tij).array() * nonDose_LIN.col(tij).array().pow(-1).array() * Tdd0.col(ijk).array();
                    } else if (tform[ij]=="plin") {
                        Rd.col(ij) =   TTerm.col(tij).array() * nonDose_PLIN.col(tij).array().pow(-1).array() * Td0.col(ij).array();
                        Rdd.col(ijk) = TTerm.col(tij).array() * nonDose_PLIN.col(tij).array().pow(-1).array() * Tdd0.col(ijk).array();
                    } else if (tform[ij]=="loglin") {
                        Rd.col(ij) =   TTerm.col(tij).array() * nonDose_LOGLIN.col(tij).array().pow(-1).array() * Td0.col(ij).array();
                        Rdd.col(ijk) = TTerm.col(tij).array() * nonDose_LOGLIN.col(tij).array().pow(-1).array() * Tdd0.col(ijk).array();
                    }
                } else {
                    if (tij==tjk){
                        if (Dose_Iden.find(as< string>(tform[ij])) != Dose_Iden.end()){
                            if (Dose_Iden.find(as< string>(tform[jk])) != Dose_Iden.end()){
                                ;
                            } else if (tform[jk]=="lin") {
                                Rdd.col(ijk) = TTerm.col(tij).array() * Dose.col(tij).array().pow(-1).array() * Td0.col(ij).array() * nonDose_LIN.col(tij).array().pow(-1).array() * Td0.col(jk).array();
                            } else if (tform[jk]=="plin") {
                                Rdd.col(ijk) = TTerm.col(tij).array() * Dose.col(tij).array().pow(-1).array() * Td0.col(ij).array() * nonDose_PLIN.col(tij).array().pow(-1).array() * Td0.col(jk).array();
                            } else if (tform[jk]=="loglin") {
                                Rdd.col(ijk) = TTerm.col(tij).array() * Dose.col(tij).array().pow(-1).array() * Td0.col(ij).array() * nonDose_LOGLIN.col(tij).array().pow(-1).array() * Td0.col(jk).array();
                            }
                        } else {
                            if (Dose_Iden.find(as< string>(tform[jk])) != Dose_Iden.end()){
                                if (tform[ij]=="lin") {
                                    Rdd.col(ijk) = TTerm.col(tij).array() * Dose.col(tij).array().pow(-1).array() * Td0.col(jk).array() * nonDose_LIN.col(tij).array().pow(-1).array() * Td0.col(ij).array();
                                } else if (tform[ij]=="plin") {
                                    Rdd.col(ijk) = TTerm.col(tij).array() * Dose.col(tij).array().pow(-1).array() * Td0.col(jk).array() * nonDose_PLIN.col(tij).array().pow(-1).array() * Td0.col(ij).array();
                                } else if (tform[ij]=="loglin") {
                                    Rdd.col(ijk) = TTerm.col(tij).array() * Dose.col(tij).array().pow(-1).array() * Td0.col(jk).array() * nonDose_LOGLIN.col(tij).array().pow(-1).array() * Td0.col(ij).array();
                                }
                            } else {
                                if (tform[ij]=="loglin") {
                                    if( tform[jk]=="lin") {
                                        Rdd.col(ijk) = TTerm.col(tij).array() * nonDose_LOGLIN.col(tij).array().pow(-1).array() * Td0.col(ij).array() * nonDose_LIN.col(tij).array().pow(-1).array() * Td0.col(jk).array();
                                    } else if (tform[jk]=="plin") {
                                        Rdd.col(ijk) = TTerm.col(tij).array() * nonDose_LOGLIN.col(tij).array().pow(-1).array() * Td0.col(ij).array() * nonDose_PLIN.col(tij).array().pow(-1).array() * Td0.col(jk).array();
                                    } else if (tform[jk]=="loglin") {
                                        Rdd.col(ijk) = TTerm.col(tij).array() * nonDose_LOGLIN.col(tij).array().pow(-1).array() * Td0.col(ij).array() * nonDose_LOGLIN.col(tij).array().pow(-1).array() * Td0.col(jk).array();
                                    }
                                } else if (tform[jk]=="loglin") {
                                    if( tform[ij]=="lin") {
                                        Rdd.col(ijk) = TTerm.col(tij).array() * nonDose_LOGLIN.col(tij).array().pow(-1).array() * Td0.col(jk).array() * nonDose_LIN.col(tij).array().pow(-1).array() * Td0.col(ij).array();
                                    } else if (tform[ij]=="plin") {
                                        Rdd.col(ijk) = TTerm.col(tij).array() * nonDose_LOGLIN.col(tij).array().pow(-1).array() * Td0.col(jk).array() * nonDose_PLIN.col(tij).array().pow(-1).array() * Td0.col(ij).array();
                                    }
                                }
                            }
                        }
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
                            Rd.col(ij) =  R.col(0).array() * Dose.col(tij).array().pow(-1).array() * Td0.col(ij).array();
                            Rdd.col(ij) =  R.col(0).array() * Dose.col(tij).array().pow(-1).array() * Tdd0.col(ijk).array();
                        } else if (tform[ij]=="lin") {
                            Rd.col(ij) =  R.col(0).array() * nonDose_LIN.col(tij).array().pow(-1).array() * Td0.col(ij).array();
                            Rdd.col(ij) =  R.col(0).array() * nonDose_LIN.col(tij).array().pow(-1).array() * Tdd0.col(ijk).array();
                        } else if (tform[ij]=="plin") {
                            Rd.col(ij) =  R.col(0).array() * nonDose_PLIN.col(tij).array().pow(-1).array() * Td0.col(ij).array();
                            Rdd.col(ij) =  R.col(0).array() * nonDose_PLIN.col(tij).array().pow(-1).array() * Tdd0.col(ijk).array();
                        } else if (tform[ij]=="loglin") {
                            Rd.col(ij) =  R.col(0).array() * nonDose_LOGLIN.col(tij).array().pow(-1).array() * Td0.col(ij).array();
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
                                    ;
                                } else if (tform[jk]=="lin") {
                                    Rdd.col(ijk) =   R.col(0).array() * Dose.col(tij).array().pow(-1).array() * Td0.col(ij).array() * nonDose_LIN.col(tij).array().pow(-1).array() * Td0.col(jk).array();
                                } else if (tform[jk]=="plin") {
                                    Rdd.col(ijk) =   R.col(0).array() * Dose.col(tij).array().pow(-1).array() * Td0.col(ij).array() * nonDose_PLIN.col(tij).array().pow(-1).array() * Td0.col(jk).array();
                                } else if (tform[jk]=="loglin") {
                                    Rdd.col(ijk) =   R.col(0).array() * Dose.col(tij).array().pow(-1).array() * Td0.col(ij).array() * nonDose_LOGLIN.col(tij).array().pow(-1).array() * Td0.col(jk).array();
                                }
                            } else {
                                if (Dose_Iden.find(as< string>(tform[jk])) != Dose_Iden.end()){
                                    if (tform[ij]=="lin") {
                                        Rdd.col(ij) =   R.col(0).array() * Dose.col(tjk).array().pow(-1).array() * Td0.col(jk).array() * nonDose_LIN.col(tij).array().pow(-1).array() * Td0.col(ij).array();
                                    } else if (tform[ij]=="plin") {
                                        Rdd.col(ij) =   R.col(0).array() * Dose.col(tjk).array().pow(-1).array() * Td0.col(jk).array() * nonDose_PLIN.col(tij).array().pow(-1).array() * Td0.col(ij).array();
                                    } else if (tform[ij]=="loglin") {
                                        Rdd.col(ij) =   R.col(0).array() * Dose.col(tjk).array().pow(-1).array() * Td0.col(jk).array() * nonDose_LOGLIN.col(tij).array().pow(-1).array() * Td0.col(ij).array();
                                    }
                                } else {
                                    if (tform[ij]=="loglin") {
                                        if( tform[jk]=="lin") {
                                            Rdd.col(ijk) = R.col(0).array() * nonDose_LOGLIN.col(tij).array().pow(-1).array() * Td0.col(ij).array() * nonDose_LIN.col(tij).array().pow(-1).array() * Td0.col(jk).array();
                                        } else if (tform[jk]=="plin") {
                                            Rdd.col(ijk) = R.col(0).array() * nonDose_LOGLIN.col(tij).array().pow(-1).array() * Td0.col(ij).array() * nonDose_PLIN.col(tij).array().pow(-1).array() * Td0.col(jk).array();
                                        } else if (tform[jk]=="loglin") {
                                            Rdd.col(ijk) = R.col(0).array() * nonDose_LOGLIN.col(tij).array().pow(-1).array() * Td0.col(ij).array() * nonDose_LOGLIN.col(tij).array().pow(-1).array() * Td0.col(jk).array();
                                        }
                                    } else if (tform[jk]=="loglin") {
                                        if( tform[ij]=="lin") {
                                            Rdd.col(ijk) = R.col(0).array() * nonDose_LOGLIN.col(tij).array().pow(-1).array() * Td0.col(jk).array() * nonDose_LIN.col(tij).array().pow(-1).array() * Td0.col(ij).array();
                                        } else if (tform[ij]=="plin") {
                                            Rdd.col(ijk) = R.col(0).array() * nonDose_LOGLIN.col(tij).array().pow(-1).array() * Td0.col(jk).array() * nonDose_PLIN.col(tij).array().pow(-1).array() * Td0.col(ij).array();
                                        }
                                    }
                                }
                            }
                        } else {
                            if (Dose_Iden.find(as< string>(tform[ij])) != Dose_Iden.end()){
                                if (Dose_Iden.find(as< string>(tform[jk])) != Dose_Iden.end()){
                                    ;
                                } else if (tform[jk]=="lin") {
                                    Rdd.col(ijk) =   TTerm.col(fir).array() * Dose.col(tij).array().pow(-1).array() * Td0.col(ij).array() * nonDose_LIN.col(tij).array().pow(-1).array() * Td0.col(jk).array();
                                } else if (tform[jk]=="plin") {
                                    Rdd.col(ijk) =   TTerm.col(fir).array() * Dose.col(tij).array().pow(-1).array() * Td0.col(ij).array() * nonDose_PLIN.col(tij).array().pow(-1).array() * Td0.col(jk).array();
                                } else if (tform[jk]=="loglin") {
                                    Rdd.col(ijk) =   TTerm.col(fir).array() * Dose.col(tij).array().pow(-1).array() * Td0.col(ij).array() * nonDose_LOGLIN.col(tij).array().pow(-1).array() * Td0.col(jk).array();
                                }
                            } else {
                                if (Dose_Iden.find(as< string>(tform[jk])) != Dose_Iden.end()){
                                    if (tform[ij]=="lin") {
                                        Rdd.col(ij) =   TTerm.col(fir).array() * Dose.col(tjk).array().pow(-1).array() * Td0.col(jk).array() * nonDose_LIN.col(tij).array().pow(-1).array() * Td0.col(ij).array();
                                    } else if (tform[ij]=="plin") {
                                        Rdd.col(ij) =   TTerm.col(fir).array() * Dose.col(tjk).array().pow(-1).array() * Td0.col(jk).array() * nonDose_PLIN.col(tij).array().pow(-1).array() * Td0.col(ij).array();
                                    } else if (tform[ij]=="loglin") {
                                        Rdd.col(ij) =   TTerm.col(fir).array() * Dose.col(tjk).array().pow(-1).array() * Td0.col(jk).array() * nonDose_LOGLIN.col(tij).array().pow(-1).array() * Td0.col(ij).array();
                                    }
                                } else {
                                    if (tform[ij]=="loglin") {
                                        if( tform[jk]=="lin") {
                                            Rdd.col(ijk) = TTerm.col(fir).array() * nonDose_LOGLIN.col(tij).array().pow(-1).array() * Td0.col(ij).array() * nonDose_LIN.col(tij).array().pow(-1).array() * Td0.col(jk).array();
                                        } else if (tform[jk]=="plin") {
                                            Rdd.col(ijk) = TTerm.col(fir).array() * nonDose_LOGLIN.col(tij).array().pow(-1).array() * Td0.col(ij).array() * nonDose_PLIN.col(tij).array().pow(-1).array() * Td0.col(jk).array();
                                        } else if (tform[jk]=="loglin") {
                                            Rdd.col(ijk) = TTerm.col(fir).array() * nonDose_LOGLIN.col(tij).array().pow(-1).array() * Td0.col(ij).array() * nonDose_LOGLIN.col(tij).array().pow(-1).array() * Td0.col(jk).array();
                                        }
                                    } else if (tform[jk]=="loglin") {
                                        if( tform[ij]=="lin") {
                                            Rdd.col(ijk) = TTerm.col(fir).array() * nonDose_LOGLIN.col(tij).array().pow(-1).array() * Td0.col(jk).array() * nonDose_LIN.col(tij).array().pow(-1).array() * Td0.col(ij).array();
                                        } else if (tform[ij]=="plin") {
                                            Rdd.col(ijk) = TTerm.col(fir).array() * nonDose_LOGLIN.col(tij).array().pow(-1).array() * Td0.col(jk).array() * nonDose_PLIN.col(tij).array().pow(-1).array() * Td0.col(ij).array();
                                        }
                                    }
                                }
                            }
                        }
                    } else if ((tij==fir)||(tjk==fir)){
                        if (Dose_Iden.find(as< string>(tform[ij])) != Dose_Iden.end()){
                            if (Dose_Iden.find(as< string>(tform[jk])) != Dose_Iden.end()){
                                Rdd.col(ijk) = TTerm.col(tij).array() * TTerm.col(tjk).array() * Dose.col(tij).array().pow(-1).array() * Td0.col(ij).array() * Dose.col(tjk).array().pow(-1).array() * Td0.col(jk).array();
                            } else if (tform[jk]=="lin") {
                                Rdd.col(ijk) = TTerm.col(tij).array() * TTerm.col(tjk).array() * Dose.col(tij).array().pow(-1).array() * Td0.col(ij).array() * Dose.col(tjk).array().pow(-1).array() * Td0.col(jk).array();
                            } else if (tform[jk]=="plin") {
                                Rdd.col(ijk) = TTerm.col(tij).array() * TTerm.col(tjk).array() * Dose.col(tij).array().pow(-1).array() * Td0.col(ij).array() * Dose.col(tjk).array().pow(-1).array() * Td0.col(jk).array();
                            } else if (tform[jk]=="loglin") {
                                Rdd.col(ijk) = TTerm.col(tij).array() * TTerm.col(tjk).array() * Dose.col(tij).array().pow(-1).array() * Td0.col(ij).array() * Dose.col(tjk).array().pow(-1).array() * Td0.col(jk).array();
                            }
                        } else {
                            if (Dose_Iden.find(as< string>(tform[jk])) != Dose_Iden.end()){
                                if (tform[ij]=="lin") {
                                    Rdd.col(ij) =   TTerm.col(tij).array() * TTerm.col(tjk).array() * Dose.col(tjk).array().pow(-1).array() * Td0.col(jk).array() * nonDose_LIN.col(tij).array().pow(-1).array() * Td0.col(ij).array();
                                } else if (tform[ij]=="plin") {
                                    Rdd.col(ij) =   TTerm.col(tij).array() * TTerm.col(tjk).array() * Dose.col(tjk).array().pow(-1).array() * Td0.col(jk).array() * nonDose_PLIN.col(tij).array().pow(-1).array() * Td0.col(ij).array();
                                } else if (tform[ij]=="loglin") {
                                    Rdd.col(ij) =   TTerm.col(tij).array() * TTerm.col(tjk).array() * Dose.col(tjk).array().pow(-1).array() * Td0.col(jk).array() * nonDose_LOGLIN.col(tij).array().pow(-1).array() * Td0.col(ij).array();
                                }
                            } else {
                                if (tform[ij]=="loglin") {
                                    if( tform[jk]=="lin") {
                                        Rdd.col(ijk) = TTerm.col(tij).array() * TTerm.col(tjk).array() * nonDose_LOGLIN.col(tij).array().pow(-1).array() * Td0.col(ij).array() * nonDose_LIN.col(tij).array().pow(-1).array() * Td0.col(jk).array();
                                    } else if (tform[jk]=="plin") {
                                        Rdd.col(ijk) = TTerm.col(tij).array() * TTerm.col(tjk).array() * nonDose_LOGLIN.col(tij).array().pow(-1).array() * Td0.col(ij).array() * nonDose_PLIN.col(tij).array().pow(-1).array() * Td0.col(jk).array();
                                    } else if (tform[jk]=="loglin") {
                                        Rdd.col(ijk) = TTerm.col(fir).array() * nonDose_LOGLIN.col(tij).array().pow(-1).array() * Td0.col(ij).array() * nonDose_LOGLIN.col(tij).array().pow(-1).array() * Td0.col(jk).array();
                                    }
                                } else if (tform[jk]=="loglin") {
                                    if( tform[ij]=="lin") {
                                        Rdd.col(ijk) = TTerm.col(tij).array() * TTerm.col(tjk).array() * nonDose_LOGLIN.col(tij).array().pow(-1).array() * Td0.col(jk).array() * nonDose_LIN.col(tij).array().pow(-1).array() * Td0.col(ij).array();
                                    } else if (tform[ij]=="plin") {
                                        Rdd.col(ijk) = TTerm.col(tij).array() * TTerm.col(tjk).array() * nonDose_LOGLIN.col(tij).array().pow(-1).array() * Td0.col(jk).array() * nonDose_PLIN.col(tij).array().pow(-1).array() * Td0.col(ij).array();
                                    }
                                }
                            }
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
        Rd = T0.array().pow(-1).array() * Td0.array() * Te.rowwise().replicate(totalnum).array();
        //
        //
        Rd = Td0.array();
        
        #pragma omp parallel for schedule(dynamic) num_threads(nthreads)
        for (int ijk=0;ijk<totalnum;ijk++){
            int tij = Term_n[ijk];
            if (Dose_Iden.find(as< string>(tform[ijk])) != Dose_Iden.end()){
                Rd.col(ijk) = R.col(0).array() * TTerm.col(tij).array() * TTerm_p.col(tij).array().pow(-1).array() * Td0.array().col(ijk).array() * Dose.col(tij).array().pow(-1).array();
            } else if (tform[ijk]=="lin") {
                Rd.col(ijk) = R.col(0).array() * TTerm.col(tij).array() * TTerm_p.col(tij).array().pow(-1).array() * Td0.array().col(ijk).array() * nonDose_LIN.col(tij).array().pow(-1).array();
            } else if (tform[ijk]=="plin") {
                Rd.col(ijk) = R.col(0).array() * TTerm.col(tij).array() * TTerm_p.col(tij).array().pow(-1).array() * Td0.array().col(ijk).array() * nonDose_PLIN.col(tij).array().pow(-1).array();
            } else if (tform[ijk]=="loglin") {
                Rd.col(ijk) = R.col(0).array() * TTerm.col(tij).array() * TTerm_p.col(tij).array().pow(-1).array() * Td0.array().col(ijk).array() * nonDose_LOGLIN.col(tij).array().pow(-1).array();
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
                    Rdd.col(ijk) = R.col(0).array() * TTerm.col(tij).array() * TTerm_p.col(tij).array().pow(-1).array() * Tdd0.array().col(ijk).array() * nonDose_LOGLIN.col(tij).array().pow(-1).array();
                }
            } else {
                if (tij!=tjk){
                    if (Dose_Iden.find(as< string>(tform[ij])) != Dose_Iden.end()){
                        Rdd.col(ijk) = Rd.col(jk).array() * TTerm.col(tij).array() * TTerm_p.col(tij).array().pow(-1).array() * Td0.array().col(ij).array() * Dose.col(tij).array().pow(-1).array();
                    } else if (tform[ij]=="lin") {
                        Rdd.col(ijk) = Rd.col(jk).array() * TTerm.col(tij).array() * TTerm_p.col(tij).array().pow(-1).array() * Td0.array().col(ij).array() * nonDose_LIN.col(tij).array().pow(-1).array();
                    } else if (tform[ij]=="plin") {
                        Rdd.col(ijk) = Rd.col(jk).array() * TTerm.col(tij).array() * TTerm_p.col(tij).array().pow(-1).array() * Td0.array().col(ij).array() * nonDose_PLIN.col(tij).array().pow(-1).array();
                    } else if (tform[ij]=="loglin") {
                        Rdd.col(ijk) = Rd.col(jk).array() * TTerm.col(tij).array() * TTerm_p.col(tij).array().pow(-1).array() * Tdd0.array().col(ij).array() * nonDose_LOGLIN.col(tij).array().pow(-1).array();
                    }
                } else {
                    if (Dose_Iden.find(as< string>(tform[ij])) != Dose_Iden.end()){
                        Rdd.col(ijk) = Rd.col(jk).array() * TTerm.col(tij).array() * TTerm_p.col(tij).array().pow(-1).array() * Td0.array().col(ij).array() * Dose.col(tij).array().pow(-1).array();
                    } else if ((tform[ij]=="lin")&&(tform[jk]=="loglin")) {
                        Rdd.col(ijk) = Rd.col(jk).array() * TTerm.col(tij).array() * TTerm_p.col(tij).array().pow(-1).array() * Td0.array().col(ij).array() * nonDose_LIN.col(tij).array().pow(-1).array();
                    } else if ((tform[ij]=="plin")&&(tform[jk]=="loglin")) {
                        Rdd.col(ijk) = Rd.col(jk).array() * TTerm.col(tij).array() * TTerm_p.col(tij).array().pow(-1).array() * Td0.array().col(ij).array() * nonDose_PLIN.col(tij).array().pow(-1).array();
                    } else if (tform[ij]=="loglin") {
                        Rdd.col(ijk) = Rd.col(jk).array() * TTerm.col(tij).array() * TTerm_p.col(tij).array().pow(-1).array() * Tdd0.array().col(ij).array() * nonDose_LOGLIN.col(tij).array().pow(-1).array();
                    }
                }
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
    for (int ijk=0;ijk<(totalnum*(totalnum+1)/2);ijk++){
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
//    cout << "entered" << endl;
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
    //            cout << ijk << " " << j << endl;
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
    //            cout << RiskFail(j,1) << " " << RiskFail(j,0) << " " << InGroup.size() << endl;
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
    //            cout << ijk << " " << j << endl;
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
    //            cout << RiskFail(j,1) << " " << RiskFail(j,0) << " " << InGroup.size() << endl;
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


void Calc_LogLik(const int& nthreads,const IntegerMatrix& RiskFail, const vector<string>&  RiskGroup, const int& totalnum, const int& ntime, const MatrixXd& R, const MatrixXd& Rd, const MatrixXd& Rdd, const MatrixXd& RdR, const MatrixXd& RddR,const MatrixXd& Rls1,const MatrixXd& Rls2,const MatrixXd& Rls3,const MatrixXd& Lls1,const MatrixXd& Lls2,const MatrixXd& Lls3, vector<double>& Ll, vector<double>& Lld, vector<double>& Lldd, bool debugging){
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
                Ldcs << Lls1(j,0), Lls2(j,ij), Lls2(j,jk), Lls3(j,ijk);
                for (int i = 0; i < dj; i++){ //adds in the efron approximation terms
                    Ldm.row(i) = (-double(i) / double(dj)) *Ldcs.array();
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
    //            if (ijk==0){
    //                LL1[j+1]=Ld1 - Rs1;
    //                LL2[j+1]=Ld2 - Rs2;
    //                LL3[j+1]=Ld3 - Rs3;
    //            }
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
                Ldcs << Lls1(j,0), Lls2(j,ij), Lls2(j,jk), Lls3(j,ijk);
                for (int i = 0; i < dj; i++){ //adds in the efron approximation terms
                    Ldm.row(i) = (-double(i) / double(dj)) *Ldcs.array();
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
    //            if (ijk==0){
    //                LL1[j+1]=Ld1 - Rs1;
    //                LL2[j+1]=Ld2 - Rs2;
    //                LL3[j+1]=Ld3 - Rs3;
    //            }
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

void Calc_Change(const int& nthreads, const int& totalnum,const int& fir, const int& der_iden, const double& dbeta_cap, const double& dose_abs_max, const double& lr, const double& abs_max, const vector<double>& Ll, const vector<double>& Lld, const vector<double>& Lldd, vector<double>& dbeta, const bool change_all,const StringVector&   tform, const double& dint, IntegerVector KeepConstant, bool debugging){
    //
    //Calc_Change( nthreads, totalnum, fir, dbeta_cap, dose_abs_max, lr, abs_max, Ll, Lld, Lldd, dbeta, change_all, tform, dint, KeepConstant, debugging);
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
    //                dbeta[ijk] = 0.0;
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
    //                dbeta[ijk] = 0.0;
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

