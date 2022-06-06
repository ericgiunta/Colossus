#include <RcppEigen.h>
//#include <RcppParallel.h>
#include <omp.h>
#include "AMFIT_MODEL.h"
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include<chrono>
#include<random>


// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::plugins(openmp)]]
using namespace std;
using namespace Rcpp;
using namespace Eigen;
//using namespace RcppParallel;

using Eigen::Map;
using Eigen::MatrixXd;
using Eigen::VectorXd;
//using Eigen::DenseBase::Random;
using Rcpp::as;


template <typename T> int sign(T val) {
    return (T(0) < val) - (val < T(0));
}

// [[Rcpp::export]]
NumericVector amfit_transition(NumericVector a_lin,NumericVector a_loglin,NumericVector a_plin, NumericMatrix x_lin, NumericMatrix x_loglin, NumericMatrix x_plin,int fir,string modelform,NumericMatrix dfe, NumericVector include_bool){
    //----------------------------------------------------------------------------------------------------------------//
    Map<VectorXd> beta_lin(as<Map<VectorXd> >(a_lin));
    Map<VectorXd> beta_loglin(as<Map<VectorXd> >(a_loglin));
    Map<VectorXd> beta_plin(as<Map<VectorXd> >(a_plin));
    const Map<MatrixXd> df_lin(as<Map<MatrixXd> >(x_lin));
    const Map<MatrixXd> df_loglin(as<Map<MatrixXd> >(x_loglin));
    const Map<MatrixXd> df_plin(as<Map<MatrixXd> >(x_plin));
    const Map<MatrixXd> PyrC(as<Map<MatrixXd> >(dfe));
    string IterStyle = "constant";
    // Converts from Rcpp types to efficient Eigen types
    //----------------------------------------------------------------------------------------------------------------//
    VectorXd res = LogLik_AMFIT(PyrC,beta_lin,beta_loglin,beta_plin,df_lin,df_loglin,df_plin,fir,modelform,IterStyle,include_bool);
    //----------------------------------------------------------------------------------------------------------------//
    return wrap(res);
}

VectorXd LogLik_AMFIT( MatrixXd PyrC, VectorXd beta_linT,VectorXd beta_loglinT,VectorXd beta_plinT,MatrixXd df_lin,MatrixXd df_loglin,MatrixXd df_plin,int fir,string modelform, string IterStyle, NumericVector include_bool){
//    cout << "start Risk" << endl;
    VectorXd beta_0(1);
    VectorXd df0(df_lin.rows(), 1); // stores memory for the derivative term parameter and column
    string tform;
    int totalnum = 0;
    int ind0 = fir; //starts by always using the derivative of the "first" term
    double lr = 0.95; //learning rate
    double Ll = 0.0;
    double Ll_old = 0.0;
    double Lld = 0.0;
    double Lldd = 0.0;
    int maxiter = 500; //maximum iterations
    int halfmax = 10; //maximum halfsteps/derivative steps per iteration
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
    vector <double> Ll_check(2); //preallocates a vector of log-likelihood progress comparisons
    double epsilon = 1e-5;
    double totem = PyrC.rows(); //precalculates how many rows
    //
    double dbeta_max = 0.1; //maximum parameter change allowed
    double dbeta_kill = 1e6; //maximum estimated parameter before skipping
    double deriv_epsilon = 1e-3;//derivative threshold for iteration stopping
    //
    srand (time(NULL));
    //
    using namespace std::chrono;
    cout << "START_NEW" << endl;
    time_point<high_resolution_clock> start_point, end_point;
    start_point = high_resolution_clock::now();
    auto start = time_point_cast<microseconds>(start_point).time_since_epoch().count();
    end_point = high_resolution_clock::now();
    auto end = time_point_cast<microseconds>(end_point).time_since_epoch().count(); //The time duration is tracked
    // ---------------------------------------------
    // To Start, needs to seperate the derivative term
    // ---------------------------------------------
    if (include_bool[0]==1){
        if (ind0 < beta_lin.size()){
            // one exists and is one
            beta_0 = beta_lin.segment(ind0,1);
            df0 = df_lin.col(ind0);
            tform = "lin";
            //
        } else {
            //one exists and its not one
            ind0 = ind0 - beta_lin.size();
            if (include_bool[1]==1){
                if (ind0 < beta_loglin.size()){
                    //one and two exists and is two
                    beta_0 = beta_loglin.segment(ind0,1);
                    df0 = df_loglin.col(ind0);
                    tform = "loglin";
                    //
                } else{
                    //one exists, two does, must be three
                    if (include_bool[2]!=1){
                        throw invalid_argument( "Are all three used?" );
                    }
                    ind0 = ind0 - beta_loglin.size();
                    beta_0 = beta_plin.segment(ind0,1);
                    df0 = df_plin.col(ind0);
                    tform = "plin";
                    //
                }
            } else{
                //one exists, and two doesn't exist, must be three
                if (include_bool[2]!=1){
                    throw invalid_argument( "Are all first and third used?" );
                }
                beta_0 = beta_plin.segment(ind0,1);
                df0 = df_plin.col(ind0);
                tform = "plin";
                //
            }
        }
    }else{
        //one doesn't exist
        if (include_bool[1]==1){
            if (ind0 < beta_loglin.size()){
                //one doesn't exist and two exists and is two
                beta_0 = beta_loglin.segment(ind0,1);
//                cout << ind0 << ", " << beta_0 << ", " << beta_loglin.transpose() << endl;
                df0 = df_loglin.col(ind0);
                tform = "loglin";
                //
            } else{
                //one doesn't exist, two does, must be three
                if (include_bool[2]!=1){
                    throw invalid_argument( "Are all three used?" );
                }
                ind0 = ind0 - beta_loglin.size();
                beta_0 = beta_plin.segment(ind0,1);
                df0 = df_plin.col(ind0);
                tform = "plin";
                //
            }
        } else{
            //one doesn't exist, and two doesn't exist, must be three
            if (include_bool[2]!=1){
                throw invalid_argument( "Are all first and third used?" );
            }
            beta_0 = beta_plin.segment(ind0,1);
            df0 = df_plin.col(ind0);
            tform = "plin";
            //
        }
    }
    //
    end_point = high_resolution_clock::now();
    end = time_point_cast<microseconds>(end_point).time_since_epoch().count();
    // ----------------------------------------------//
    // The terminal output is in the form:
    // df101 to identify the row
    // the elapsed time since starting
    // the number of half-steps/derivative steps taken
    // the current log-likelihood
    // the current log-likelihood first derivative
    // the current log-likelihood second derivative
    // the current best log-likelihood
    // the parameter value change in that half-step/derivative step
    // the current variance calculated
    //
    // zeros denote places without calculations
    // ----------------------------------------------//
    cout<<"df101,"<<(end-start)<<",0,0,0,0,0,0"<<endl;
    //
    MatrixXd T0 = MatrixXd::Zero(df0.rows(), 1); //preallocates matrix for Derivative column term
    MatrixXd Te = MatrixXd::Zero(df0.rows(), 1); //preallocates matrix for non-Derivative column terms
    MatrixXd Td0 = MatrixXd::Zero(df0.rows(), 2); //preallocates matrix for Derivative column term derivatives
    MatrixXd R = MatrixXd::Zero(df0.rows(), 3); //preallocates matrix for Risks and derivatives
    if ((modelform=="A")||(modelform=="PA")||(modelform=="PAE")){ //same process used for all of the additive type models
        Te = Te.array() * 0; //checks that the inital sum of terms are zero
        //
        T0 = df0 * beta_0; //calculates the derivative column values
        if (tform=="lin") {
            Td0 << df0, 0*df0;
        } else if (tform=="loglin") {
            T0 = T0.array().exp();
            Td0 << df0.array() * T0.array(),df0.array().pow(2) * T0.array();
        } else {
            T0 = 1 + T0.array();
            Td0 << df0, 0*df0;
        }
        // combines every used term values
        if (include_bool[0]==1){
            Te = Te.array() + (df_lin.array().rowwise() * beta_lin.transpose().array()).matrix().rowwise().sum().array();
        }
        if (include_bool[1]==1){
            Te = Te.array() + (df_loglin.array().rowwise() * beta_loglin.transpose().array()).array().exp().rowwise().sum().array();
        }
        if (include_bool[2]==1){
            Te = Te.array() + (1 + df_plin.array().rowwise() * beta_plin.transpose().array()).matrix().rowwise().sum().array();
        }
        // computes intial risk and derivatives
        if (modelform=="A"){
            R << Te.array(),Td0.col(0),Td0.col(1);
        } else if (modelform=="PA"){
            Te = Te.array() - T0.array();
            R << T0.array() * Te.array(), Td0.col(0).array() * Te.array(), Td0.col(1).array() * Te.array();
        } else if (modelform=="PAE"){
            Te = Te.array() + 1 - T0.array();
            R << T0.array() * Te.array(),  Td0.col(0).array() * Te.array(), Td0.col(1).array() * Te.array();
        }
    }else if (modelform=="M"){
        Te = Te.array() * 0 + 1; //verifies the initial term product is 1
        //
        T0 = df0 * beta_0; //calculates the derivative column values
        if (tform=="lin") {
            Td0 << df0, 0*df0;
        } else if (tform=="loglin") {
            T0 = T0.array().exp();
            Td0 << df0.array() * T0.array(),df0.array().pow(2) * T0.array();
        } else {
            T0 = 1 + T0.array();
            Td0 << df0, 0*df0;
        }
        // combines every used term values
        if (include_bool[0]==1){
            Te = Te.array() * (df_lin.array().rowwise() * beta_lin.transpose().array()).matrix().rowwise().prod().array();
        }
        if (include_bool[1]==1){
            Te = Te.array() * (df_loglin.array().rowwise() * beta_loglin.transpose().array()).array().rowwise().sum().array().exp().array();
        }
        if (include_bool[2]==1){
            Te = Te.array() * (1 + df_plin.array().rowwise() * beta_plin.transpose().array()).matrix().rowwise().prod().array();
        }
        // computes intial risk and derivatives
        Te = Te.array() * T0.array().pow(-1);
        R << T0.array() * Te.array(),Td0.col(0).array() * Te.array(), Td0.col(1).array() * Te.array();
    } else if (modelform=="GM"){
        //currently isn't implemented, it can be calculated but not optimized the same way
        throw invalid_argument( "GM isn't implemented" );
    } else {
        throw invalid_argument( "Model isn't implemented" );
    }
    R =(R.array().isFinite()).select(R,0);
    //
    Ll = 0.0;
    Lld = 0.0;
    Lldd = 0.0;
    //
    end_point = high_resolution_clock::now();
    end = time_point_cast<microseconds>(end_point).time_since_epoch().count();
    cout<<"df101,"<<(end-start)<<",0,0,0,0,0,0"<<endl;
    //
    // Calculates the log-likelihood and derivatives
    //
    VectorXd temp(R.rows(),R.cols());
    temp = (PyrC.col(1).array() * (PyrC.col(0).array() * R.col(0).array()).array().log()).array() - (PyrC.col(0).array() * R.col(0).array());
    Ll =  (temp.array().isFinite()).select(temp,0).sum();
    temp = (PyrC.col(1).array() * R.col(1).array() / R.col(0).array()).array() - (PyrC.col(0).array() * R.col(1).array());
    Lld = (temp.array().isFinite()).select(temp,0).sum();
    temp = (PyrC.col(1).array() * (R.col(0).pow(-1).array() * R.col(2).array() - (R.col(1).array() * R.col(0).pow(-1).array()).pow(2) )).array() - (PyrC.col(0).array() * R.col(2).array());
    Lldd = (temp.array().isFinite()).select(temp,0).sum();

    // Calculates the variance
    double dev = pow((PyrC.col(1).array() - PyrC.col(0).array() * R.col(0).array()).pow(2).sum(),.5)/totem;
    // --------------------------
    // stores the value to compare to in the iteration step
    // --------------------------
    Ll_check[0] = Ll;
    //
    end_point = high_resolution_clock::now();
    end = time_point_cast<microseconds>(end_point).time_since_epoch().count();
    cout<<"df101,"<<(end-start)<<",0,"<< Ll<<","<<Lld<<","<<Lldd<<","<<Ll<<",0,"<<dev<<endl;
    // --------------------------
    // uses the estimate of next step and direction increasing log-likelihood
    // --------------------------
    double dbeta=0;
    if ((Lld > 0)&&(Lldd<0)){
        dbeta = -lr * Lld / Lldd;// decreases left, zeros right, ratio left
    } else if ((Lld<0)&&(Lldd>0)){
        dbeta = lr * Lld / Lldd;// decreases right, zeros right, ratio left
    } else if ((Lld > 0)&&(Lldd>0)){
        dbeta = lr * Lld / Lldd;// decreases left, zeros left, ratio right
    } else if ((Lld<0)&&(Lldd<0)){
        dbeta = -lr * Lld / Lldd;// decreases right, zeros left, ratio right
    } else {
        dbeta=0.0;
    }
    //
    if (abs(dbeta)>dbeta_max){
        dbeta = dbeta_max * sign(dbeta);
    }
    //
    // --------------------------
    // always starts from intial guess
    // --------------------------
    double beta_p = beta_0[0];
    double beta_c = beta_p;
    double beta_a = beta_0[0];
    double beta_best = beta_p;
    double Ll_best = Ll;
    int halves = 0;
    int i = fir;
    double beta_b=dbeta;
    double beta_prev=dbeta;
    while (((halves==0)||(abs(beta_b)>epsilon))&&(halves<halfmax)&&(dbeta!=0.0)&&(abs(Lld)>deriv_epsilon)){
        halves++;
        beta_c = beta_a + dbeta;
        if (modelform=="A"){//The risks can be updated instead of recalculated for the models/terms allowed
            if (tform == "lin"){
                R.col(0) = R.col(0).array() + (beta_c - beta_p) * df0.col(0).array();
            } else if (tform == "plin"){
                R.col(0) = R.col(0).array() + (beta_c - beta_p) * df0.col(0).array();
            } else{
                R.col(0) = R.col(0).array() + (beta_c * df0.col(0)).array().exp() - (beta_p * df0.col(0)).array().exp();
                R.col(1) = R.col(1).array() * (beta_c * df0.col(0)).array().exp();
                R.col(2) = R.col(2).array() * (beta_c * df0.col(0)).array().exp();
            }
        } else if (modelform=="PA"){
            if (tform == "lin"){
                if (i==fir){
                    R.col(0) = R.col(0) * beta_c / beta_p;
                    T0 = T0 * beta_c / beta_p;
                } else {
                    R.col(0) = R.col(0).array() + (beta_c - beta_p) * df0.col(0).array() * T0.array();
                    Te = Te.array() + (beta_c - beta_p) * df0.col(0).array();
                }
            } else if (tform == "plin"){
                if (i==fir){
                    R.col(0) = (R.col(0).array() - Te.array()) * beta_c / beta_p + Te.array();
                    T0 = (T0.array() - 1) * beta_c / beta_p + 1;
                } else {
                    R.col(0) = R.col(0).array() + (beta_c - beta_p) * df0.col(0).array() * T0.array();
                    Te = Te.array() + (beta_c - beta_p) * df0.col(0).array();
                }
            } else{
                if (i==fir){
                    R.col(0) = R.col(0).array() * ((beta_c - beta_p) * df0.col(0)).array().exp();
                    R.col(1) = R.col(1).array() * ((beta_c - beta_p) * df0.col(0)).array().exp();
                    R.col(2) = R.col(2).array() * ((beta_c - beta_p) * df0.col(0)).array().exp();
                    T0 = T0.array() * ((beta_c - beta_p) * df0.col(0)).array().exp();
                } else {
                    R.col(0) = R.col(0).array() + ((beta_c * df0.col(0)).array().exp() - (beta_p * df0.col(0)).array().exp()).array() * T0.array();
                    R.col(1) = R.col(1).array() * ((beta_c - beta_p) * df0.col(0)).array().exp();
                    R.col(2) = R.col(2).array() * ((beta_c - beta_p) * df0.col(0)).array().exp();
                    Te = Te.array() + ((beta_c * df0.col(0)).array().exp() - (beta_p * df0.col(0)).array().exp()).array();
                }
            }
        } else if (modelform=="PAE"){
            if (tform == "lin"){
                if (i==fir){
                    R.col(0) = R.col(0) * beta_c / beta_p;
                    T0 = T0 * beta_c / beta_p;
                } else {
                    R.col(0) = R.col(0).array() + (beta_c - beta_p) * df0.col(0).array() * T0.array();
                    Te = Te.array() + (beta_c - beta_p) * df0.col(0).array();
                }
            } else if (tform == "plin"){
                if (i==fir){
                    R.col(0) = (R.col(0).array() - Te.array()) * beta_c / beta_p + Te.array();
                    T0 = (T0.array() - 1) * beta_c / beta_p + 1;
                } else {
                    R.col(0) = R.col(0).array() + (beta_c - beta_p) * df0.col(0).array() * T0.array();
                    Te = Te.array() + (beta_c - beta_p) * df0.col(0).array();
                }
            } else{
                if (i==fir){
                    R.col(0) = R.col(0).array() * ((beta_c - beta_p) * df0.col(0)).array().exp();
                    R.col(1) = R.col(1).array() * ((beta_c - beta_p) * df0.col(0)).array().exp();
                    R.col(2) = R.col(2).array() * ((beta_c - beta_p) * df0.col(0)).array().exp();
                    T0 = T0.array() * ((beta_c - beta_p) * df0.col(0)).array().exp();
                } else {
                    R.col(0) = R.col(0).array() + ((beta_c * df0.col(0)).array().exp() - (beta_p * df0.col(0)).array().exp()) * T0.array();
                    R.col(1) = R.col(1).array() * ((beta_c - beta_p) * df0.col(0)).array().exp();
                    R.col(2) = R.col(2).array() * ((beta_c - beta_p) * df0.col(0)).array().exp();
                    Te = Te.array() + ((beta_c * df0.col(0)).array().exp() - (beta_p * df0.col(0)).array().exp()).array();
                }
            }
        }else if (modelform=="M"){
            if (tform == "lin"){
                R.col(0) = R.col(0) * beta_c / beta_p;
                R.col(1) = R.col(0) / beta_c;
            } else if (tform == "plin"){
                R.col(0) = R.col(0).array() * (1 + beta_c * df0.col(0).array()) / (1 + beta_p * df0.col(0).array());
                R.col(1) = R.col(1).array() * (df0.col(0).array()) / (1 + beta_c * df0.col(0).array());
            } else{
                R.col(0) = R.col(0).array() * ((beta_c - beta_p) * df0.col(0)).array().exp();
                R.col(1) = R.col(1).array() * ((beta_c - beta_p) * df0.col(0)).array().exp();
                R.col(2) = R.col(2).array() * ((beta_c - beta_p) * df0.col(0)).array().exp();
            }
        }else if (modelform=="GM"){
            throw invalid_argument( "GM isn't implemented" );
        } else {
            throw invalid_argument( "Model isn't implemented" );
        }
        Ll = 0.0;
        Lld = 0.0;
        Lldd = 0.0;
        //
        temp = (PyrC.col(1).array() * (PyrC.col(0).array() * R.col(0).array()).array().log()).array() - (PyrC.col(0).array() * R.col(0).array());
        Ll =  (temp.array().isFinite()).select(temp,0).sum();
        temp = (PyrC.col(1).array() * R.col(1).array() / R.col(0).array()).array() - (PyrC.col(0).array() * R.col(1).array());
        Lld = (temp.array().isFinite()).select(temp,0).sum();
        temp = (PyrC.col(1).array() * (R.col(0).pow(-1).array() * R.col(2).array() - (R.col(1).array() * R.col(0).pow(-1).array()).pow(2) )).array() - (PyrC.col(0).array() * R.col(2).array());
        Lldd = (temp.array().isFinite()).select(temp,0).sum();
        //
        beta_b = 0.0;
        //
        // --------------------------
        // updates the change in initial guess
        // --------------------------
        if ((abs(lr*Lld/Lldd)<dbeta_kill)&&(halves<halfmax)&&(Lldd!=0.0)&&(abs(Lld)>deriv_epsilon)){
            if (Ll > Ll_best){ //if it improves, it trys the derivative again
                if ((Lld > 0)&&(Lldd<0)){
                    beta_b = - lr * Lld / Lldd;
                } else if ((Lld<0)&&(Lldd>0)){
                    beta_b = lr * Lld / Lldd;
                } else if ((Lld > 0)&&(Lldd>0)){
                    beta_b = lr * Lld / Lldd;
                } else if ((Lld<0)&&(Lldd<0)){
                    beta_b = - lr * Lld / Lldd;
                } else {
                    beta_b= 0.0;
                }
            } else{ //if it doesn't improve, it takes a half step
                if (dbeta!=0){
                    beta_b=-0.5 * beta_prev;
                }else{
                    beta_b=0.0;
                }
            }
            beta_b = beta_b * pow(.9,halves);
            if (abs(beta_b)>dbeta_max){
                beta_b = dbeta_max * sign(beta_b);
            }
            beta_prev = beta_b;
            dbeta = dbeta + beta_b;
        }
        // --------------------------
        // updates best guess
        // --------------------------
        if (Ll > Ll_best){
            Ll_best = Ll;
            beta_best = beta_c;
        }
        dev = pow((PyrC.col(1).array() - PyrC.col(0).array() * R.col(0).array()).pow(2).sum(),.5)/totem;
        end_point = high_resolution_clock::now();
        end = time_point_cast<microseconds>(end_point).time_since_epoch().count();
        cout<<"df101,"<<(end-start)<<","<<halves<<","<< Ll<<","<<Lld<<","<<Lldd<<","<<Ll_best<<","<<beta_b<<","<<dev<<endl;
        beta_p = beta_c;
    }
    halves = 0;
    //
    beta_0[0] = beta_best;
    beta_c = beta_best;
    // --------------------------
    // sets risks to best guess values
    // --------------------------
    if (modelform=="A"){
        if (tform == "lin"){
            R.col(0) = R.col(0).array() + (beta_c - beta_p) * df0.col(0).array();
        } else if (tform == "plin"){
            R.col(0) = R.col(0).array() + (beta_c - beta_p) * df0.col(0).array();
        } else{
            R.col(0) = R.col(0).array() + (beta_c * df0.col(0)).array().exp() - (beta_p * df0.col(0)).array().exp();
            R.col(1) = R.col(1).array() * ((beta_c - beta_p) * df0.col(0)).array().exp();
            R.col(2) = R.col(2).array() * ((beta_c - beta_p) * df0.col(0)).array().exp();
        }
    } else if (modelform=="PA"){
        if (tform == "lin"){
            if (i==fir){
                R.col(0) = R.col(0) * beta_c / beta_p;
                T0 = T0 * beta_c / beta_p;
            } else {
                R.col(0) = R.col(0).array() + (beta_c - beta_p) * df0.col(0).array() * T0.array();
                Te = Te.array() + (beta_c - beta_p) * df0.col(0).array();
            }
        } else if (tform == "plin"){
            if (i==fir){
                R.col(0) = (R.col(0).array() - Te.array()) * beta_c / beta_p + Te.array();
                T0 = (T0.array() - 1) * beta_c / beta_p + 1;
            } else {
                R.col(0) = R.col(0).array() + (beta_c - beta_p) * df0.col(0).array() * T0.array();
                Te = Te.array() + (beta_c - beta_p) * df0.col(0).array();
            }
        } else{
            if (i==fir){
                R.col(0) = R.col(0).array() * ((beta_c - beta_p) * df0.col(0)).array().exp();
                R.col(1) = R.col(1).array() * ((beta_c - beta_p) * df0.col(0)).array().exp();
                R.col(2) = R.col(2).array() * ((beta_c - beta_p) * df0.col(0)).array().exp();
                T0 = T0.array() * ((beta_c - beta_p) * df0.col(0)).array().exp();
            } else {
                R.col(0) = R.col(0).array() + ((beta_c * df0.col(0)).array().exp() - (beta_p * df0.col(0)).array().exp()) * T0.array();
                R.col(1) = R.col(1).array() * ((beta_c - beta_p) * df0.col(0)).array().exp();
                R.col(2) = R.col(2).array() * ((beta_c - beta_p) * df0.col(0)).array().exp();
                Te = Te.array() + ((beta_c * df0.col(0)).array().exp() - (beta_p * df0.col(0)).array().exp()).array();
            }
        }
    } else if (modelform=="PAE"){
        if (tform == "lin"){
            if (i==fir){
                R.col(0) = R.col(0) * beta_c / beta_p;
                T0 = T0 * beta_c / beta_p;
            } else {
                R.col(0) = R.col(0).array() + (beta_c - beta_p) * df0.col(0).array() * T0.array();
                Te = Te.array() + (beta_c - beta_p) * df0.col(0).array();
            }
        } else if (tform == "plin"){
            if (i==fir){
                R.col(0) = (R.col(0).array() - Te.array()) * beta_c / beta_p + Te.array();
                T0 = (T0.array() - 1) * beta_c / beta_p + 1;
            } else {
                R.col(0) = R.col(0).array() + (beta_c - beta_p) * df0.col(0).array() * T0.array();
                Te = Te.array() + (beta_c - beta_p) * df0.col(0).array();
            }
        } else{
            if (i==fir){
                R.col(0) = R.col(0).array() * ((beta_c - beta_p) * df0.col(0)).array().exp();
                R.col(1) = R.col(1).array() * ((beta_c - beta_p) * df0.col(0)).array().exp();
                R.col(2) = R.col(2).array() * ((beta_c - beta_p) * df0.col(0)).array().exp();
                T0 = T0.array() * ((beta_c - beta_p) * df0.col(0)).array().exp();
            } else {
                R.col(0) = R.col(0).array() + ((beta_c * df0.col(0)).array().exp() - (beta_p * df0.col(0)).array().exp()) * T0.array();
                R.col(1) = R.col(1).array() * ((beta_c - beta_p) * df0.col(0)).array().exp();
                R.col(2) = R.col(2).array() * ((beta_c - beta_p) * df0.col(0)).array().exp();
                Te = Te.array() + ((beta_c * df0.col(0)).exp() - (beta_p * df0.col(0)).exp()).array();
            }
        }
    }else if (modelform=="M"){
        if (tform == "lin"){
            R.col(0) = R.col(0) * beta_c / beta_p;
            R.col(1) = R.col(0) / beta_c;
        } else if (tform == "plin"){
            R.col(0) = R.col(0).array() * (1 + beta_c * df0.col(0).array()) / (1 + beta_p * df0.col(0).array());
            R.col(1) = R.col(1).array() * (df0.col(0).array()) / (1 + beta_c * df0.col(0).array());
        } else{
            R.col(0) = R.col(0).array() * ((beta_c - beta_p) * df0.col(0)).array().exp();
            R.col(1) = df0.col(0).array() * R.col(0).array();
            R.col(2) = df0.col(0).array() * R.col(1).array();
        }
    }else if (modelform=="GM"){
        throw invalid_argument( "GM isn't implemented" );
    } else {
        throw invalid_argument( "GM isn't implemented" );
    }
    // --------------------------
    // updates the list of parameters
    // --------------------------
    if (tform == "lin"){
        beta_lin[ind0] = beta_0[0];
    } else if (tform == "loglin"){
        beta_loglin[ind0] = beta_0[0];;
    } else {
        beta_plin[ind0] = beta_0[0];
    }
    Ll_check[1] = Ll; //stores the next log-likelihood to compare against
    end_point = high_resolution_clock::now();
    end = time_point_cast<microseconds>(end_point).time_since_epoch().count();
    dev = pow((PyrC.col(1).array() - PyrC.col(0).array() * R.col(0).array()).pow(2).sum(),.5)/totem;
    cout<<"df101,"<<(end-start)<<","<<halves<<","<< Ll<<","<<Lld<<","<<Lldd<<","<<Ll_best<<",0,"<<dev<<endl;
    int iteration=0;
    vector <double> Improve(totalnum,abs(Ll_best));
    vector<double> c_sum(totalnum);
    double lower_bound = 0;
    double upper_bound = 1;
    double rand_i;
    // --------------------------
    // repeats atleast 2 x the number of total parameters, then every 2 x number of parameters checks for convergence
    // --------------------------
    while (iteration < maxiter){
        // --------------------------
        // Uses last improvements to sample the next term used
        // As it converges, the probabilities converge to a uniform distribution
        // --------------------------
        partial_sum(Improve.begin(), Improve.end(), c_sum.begin());
        upper_bound = c_sum[c_sum.size()-1];
        rand_i = lower_bound + (double)(rand()) / ((double)(RAND_MAX/(upper_bound - lower_bound)));
        for (int j=0;j<c_sum.size();j++){
            if (c_sum[j]>rand_i){
                ind0=j;
                break;
            }
        }
        //
        i = ind0;
        iteration++;
        //
        // stores the new derivative column and parameter
        //
        if (include_bool[0]==1){
            if (ind0 < beta_lin.size()){
                // one exists and is one
                beta_0 = beta_lin.segment(ind0+1,1);
                df0 = df_lin.col(ind0);
                tform = "lin";
                //
            } else {
                //one exists and its not one
                ind0 = ind0 - beta_lin.size();
                if (include_bool[1]==1){
                    if (ind0 < beta_loglin.size()){
                        //one and two exists and is two
                        beta_0 = beta_loglin.segment(ind0+1,1);
                        df0 = df_loglin.col(ind0);
                        tform = "loglin";
                        //
                    } else{
                        //one exists, two does, must be three
                        if (include_bool[2]!=1){
                            throw invalid_argument( "Are all three used?" );
                        }
                        ind0 = ind0 - beta_loglin.size();
                        beta_0 = beta_plin.segment(ind0+1,1);
                        df0 = df_plin.col(ind0);
                        tform = "plin";
                        //
                    }
                } else{
                    //one exists, and two doesn't exist, must be three
                    if (include_bool[2]!=1){
                        throw invalid_argument( "Are all first and third used?" );
                    }
                    beta_0 = beta_plin.segment(ind0+1,1);
                    df0 = df_plin.col(ind0);
                    tform = "plin";
                    //
                }
            }
        }else{
            //one doesn't exist
            if (include_bool[1]==1){
                if (ind0 < beta_loglin.size()){
                    //one doesn't exist and two exists and is two
                    beta_0 = beta_loglin.segment(ind0+1,1);
                    df0 = df_loglin.col(ind0);
                    tform = "loglin";
                    //
                } else{
                    //one doesn't exist, two does, must be three
                    if (include_bool[2]!=1){
                        throw invalid_argument( "Are all three used?" );
                    }
                    ind0 = ind0 - beta_loglin.size();
                    beta_0 = beta_plin.segment(ind0+1,1);
                    df0 = df_plin.col(ind0);
                    tform = "plin";
                    //
                }
            } else{
                //one doesn't exist, and two doesn't exist, must be three
                if (include_bool[2]!=1){
                    throw invalid_argument( "Are all first and third used?" );
                }
                beta_0 = beta_plin.segment(ind0+1,1);
                df0 = df_plin.col(ind0);
                tform = "plin";
                //
            }
        }
        end_point = high_resolution_clock::now();
        end = time_point_cast<microseconds>(end_point).time_since_epoch().count();
        cout<<"df101,"<<(end-start)<<","<<halves<<","<< Ll<<","<<Lld<<","<<Lldd<<","<<Ll_best<<",0"<<endl;
        //
        beta_c = beta_0[0];
        beta_c = beta_p;
        // The risks don't change, but the derivatives are updated for the new term
        if (modelform=="A"){
            if (tform == "lin"){
                R.col(1) = df0.col(0);
                R.col(2) = 0.0 * df0.col(0);
            } else if (tform == "plin"){
                R.col(1) = df0.col(0);
                R.col(2) = 0.0 * df0.col(0);
            } else{
                R.col(1) = df0.col(0).array() * R.col(0).array();
                R.col(2) = df0.col(0).array() * R.col(1).array();
            }
        } else if (modelform=="PA"){
            if (tform == "lin"){
                if (i==fir){
                    R.col(1) = df0.col(0).array() * Te.array();
                    R.col(2) = 0.0 * df0.col(0);
                } else {
                    R.col(1) = df0.col(0).array() * T0.array();
                    R.col(2) = 0.0 * df0.col(0).array();
                }
            } else if (tform == "plin"){
                if (i==fir){
                    R.col(1) = df0.col(0).array() * Te.array();
                    R.col(2) = 0.0 * df0.col(0);
                } else {
                    R.col(1) = df0.col(0).array() * T0.array();
                    R.col(2) = 0.0 * df0.col(0);
                }
            } else{
                if (i==fir){
                    R.col(1) = df0.col(0).array() * R.col(0).array();
                    R.col(2) = df0.col(0).array() * R.col(1).array();
                } else {
                    R.col(1) =T0.array() * df0.col(0).array() * ((beta_c) * df0.col(0)).array().exp();
                    R.col(2) = df0.col(0).array() * R.col(1).array();
                }
            }
        } else if (modelform=="PAE"){
            if (tform == "lin"){
                if (i==fir){
                    R.col(1) = df0.col(0).array() * Te.array();
                    R.col(2) = 0.0 * df0.col(0);
                } else {
                    R.col(1) = df0.col(0).array() * T0.array();
                    R.col(2) = 0.0 * df0.col(0);
                }
            } else if (tform == "plin"){
                if (i==fir){
                    R.col(1) = df0.col(0).array() * Te.array();
                    R.col(2) = 0.0 * df0.col(0);
                } else {
                    R.col(1) = df0.col(0).array() * T0.array();
                    R.col(2) = 0.0 * df0.col(0);
                }
            } else{
                if (i==fir){
                    R.col(1) = df0.col(0).array() * R.col(0).array();
                    R.col(2) = df0.col(0).array() * R.col(1).array();
                } else {
                    R.col(1) = df0.col(0).array() * ((beta_c) * df0.col(0)).array().exp() * T0.array();
                    R.col(2) = df0.col(0).array() * R.col(1).array();
                }
            }
        }else if (modelform=="M"){
            if (tform == "lin"){
                R.col(1) = R.col(0) / beta_c;
                R.col(2) = 0.0 * df0.col(0);
            } else if (tform == "plin"){
                R.col(1) = R.col(1).array() * (df0.col(0).array()) / (1 + beta_c * df0.col(0).array());
                R.col(2) = 0.0 * df0.col(0);
            } else{
                R.col(1) = df0.col(0).array() * R.col(0).array();
                R.col(2) = df0.col(0).array() * R.col(1).array();
            }
        } else {
            throw invalid_argument( "GM isn't implemented" );
        }
        end_point = high_resolution_clock::now();
        end = time_point_cast<microseconds>(end_point).time_since_epoch().count();
        cout<<"df101,"<<(end-start)<<","<<halves<<","<< Ll<<","<<Lld<<","<<Lldd<<","<<Ll_best<<",0"<<endl;
        Ll = 0.0;
        Lld = 0.0;
        Lldd = 0.0;
        // Calculates log-likelihoods
        temp = (PyrC.col(1).array() * (PyrC.col(0).array() * R.col(0).array()).array().log()).array() - (PyrC.col(0).array() * R.col(0).array());
        Ll =  (temp.array().isFinite()).select(temp,0).sum();
        temp = (PyrC.col(1).array() * R.col(1).array() / R.col(0).array()).array() - (PyrC.col(0).array() * R.col(1).array());
        Lld = (temp.array().isFinite()).select(temp,0).sum();
        temp = (PyrC.col(1).array() * (R.col(0).pow(-1).array() * R.col(2).array() - (R.col(1).array() * R.col(0).pow(-1).array()).pow(2) )).array() - (PyrC.col(0).array() * R.col(2).array());
        Lldd = (temp.array().isFinite()).select(temp,0).sum();
        //
        dev = pow((PyrC.col(1).array() - PyrC.col(0).array() * R.col(0).array()).pow(2).sum(),.5)/totem;
        end_point = high_resolution_clock::now();
        end = time_point_cast<microseconds>(end_point).time_since_epoch().count();
        halves = 0;
        cout<<"df101,"<<(end-start)<<","<<halves<<","<< Ll<<","<<Lld<<","<<Lldd<<","<<Ll_best<<",0,"<<dev<<endl;
        if ((Lld > 0)&&(Lldd<0)){
            dbeta = -lr * Lld / Lldd;
        } else if ((Lld<0)&&(Lldd>0)){
            dbeta = lr * Lld / Lldd;
        } else if ((Lld > 0)&&(Lldd>0)){
            dbeta = lr * Lld / Lldd;
        } else if ((Lld<0)&&(Lldd<0)){
            dbeta = -lr * Lld / Lldd;
        } else {
            dbeta=0.0;
        }
        //
        if (abs(dbeta)>dbeta_max){
            dbeta = dbeta_max * sign(dbeta);
        }
        beta_b=dbeta;
        beta_prev=dbeta;
        //
        beta_p = beta_0[0];
        beta_c = beta_p;
        beta_best = beta_p;
        beta_a = beta_0[0];
        Ll_best = Ll;
        //
        //Stores the score before running, minus a value to prevent 0 probability
        //
        Improve[i] = Ll_best-.001;
        //
        halves = 0;
        while (((halves==0)||(abs(beta_b)>epsilon))&&(halves<halfmax)&&(dbeta!=0.0)&&(abs(Lld)>deriv_epsilon)){
            halves++;
            beta_c = beta_a + dbeta;
            if (modelform=="A"){
                if (tform == "lin"){
                    R.col(0) = R.col(0).array() + (beta_c - beta_p) * df0.col(0).array();
                } else if (tform == "plin"){
                    R.col(0) = R.col(0).array() + (beta_c - beta_p) * df0.col(0).array();
                } else{
                    R.col(0) = R.col(0).array() + (beta_c * df0.col(0)).array().exp() - (beta_p * df0.col(0)).array().exp();
                    R.col(1) = R.col(1).array() * ((beta_c - beta_p) * df0.col(0)).array().exp();
                    R.col(2) = R.col(2).array() * ((beta_c - beta_p) * df0.col(0)).array().exp();
                }
            } else if (modelform=="PA"){
                if (tform == "lin"){
                    if (i==fir){
                        R.col(0) = R.col(0) * beta_c / beta_p;
                        T0 = T0 * beta_c / beta_p;
                    } else {
                        R.col(0) = R.col(0).array() + (beta_c - beta_p) * df0.col(0).array() * T0.array();
                        Te = Te.array() + (beta_c - beta_p) * df0.col(0).array();
                    }
                } else if (tform == "plin"){
                    if (i==fir){
                        R.col(0) = (R.col(0).array() - Te.array()) * beta_c / beta_p + Te.array();
                        T0 = (T0.array() - 1) * beta_c / beta_p + 1;
                    } else {
                        R.col(0) = R.col(0).array() + (beta_c - beta_p) * df0.col(0).array() * T0.array();
                        Te = Te.array() + (beta_c - beta_p) * df0.col(0).array();
                    }
                } else{
                    if (i==fir){
                        R.col(0) = R.col(0).array() * ((beta_c - beta_p) * df0.col(0)).array().exp();
                        R.col(1) = R.col(1).array() * ((beta_c - beta_p) * df0.col(0)).array().exp();
                        R.col(2) = R.col(2).array() * ((beta_c - beta_p) * df0.col(0)).array().exp();
                        T0 = T0.array() * ((beta_c - beta_p) * df0.col(0)).array().exp();
                    } else {
                        R.col(0) = R.col(0).array() + ((beta_c * df0.col(0)).array().exp() - (beta_p * df0.col(0)).array().exp()).array() * T0.array();
                        R.col(1) = R.col(1).array() * ((beta_c - beta_p) * df0.col(0)).array().exp();
                        R.col(2) = R.col(2).array() * ((beta_c - beta_p) * df0.col(0)).array().exp();
                        Te = Te.array() + ((beta_c * df0.col(0)).exp() - (beta_p * df0.col(0)).exp()).array();
                    }
                }
            } else if (modelform=="PAE"){
                if (tform == "lin"){
                    if (i==fir){
                        R.col(0) = R.col(0) * beta_c / beta_p;
                        T0 = T0 * beta_c / beta_p;
                    } else {
                        R.col(0) = R.col(0).array() + (beta_c - beta_p) * df0.col(0).array() * T0.array();
                        Te = Te.array() + (beta_c - beta_p) * df0.col(0).array();
                    }
                } else if (tform == "plin"){
                    if (i==fir){
                        R.col(0) = (R.col(0).array() - Te.array()) * beta_c / beta_p + Te.array();
                        T0 = (T0.array() - 1) * beta_c / beta_p + 1;
                    } else {
                        R.col(0) = R.col(0).array() + (beta_c - beta_p) * df0.col(0).array() * T0.array();
                        Te = Te.array() + (beta_c - beta_p) * df0.col(0).array();
                    }
                } else{
                    if (i==fir){
                        R.col(0) = R.col(0).array() * ((beta_c - beta_p) * df0.col(0)).array().exp();
                        R.col(1) = R.col(1).array() * ((beta_c - beta_p) * df0.col(0)).array().exp();
                        R.col(2) = R.col(2).array() * ((beta_c - beta_p) * df0.col(0)).array().exp();
                        T0 = T0.array() * ((beta_c - beta_p) * df0.col(0)).array().exp();
                    } else {
                        R.col(0) = R.col(0).array() + ((beta_c * df0.col(0)).array().exp() - (beta_p * df0.col(0)).array().exp()).array() * T0.array();
                        R.col(1) = R.col(1).array() * ((beta_c - beta_p) * df0.col(0)).array().exp();
                        R.col(2) = R.col(2).array() * ((beta_c - beta_p) * df0.col(0)).array().exp();
                        Te = Te.array() + ((beta_c * df0.col(0)).array().exp() - (beta_p * df0.col(0)).array().exp()).array();
                    }
                }
            }else if (modelform=="M"){
                if (tform == "lin"){
                    R.col(0) = R.col(0) * beta_c / beta_p;
                    R.col(1) = R.col(0) / beta_c;
                } else if (tform == "plin"){
                    R.col(0) = R.col(0).array() * (1 + beta_c * df0.col(0).array()) / (1 + beta_p * df0.col(0).array());
                    R.col(1) = R.col(1).array() * (df0.col(0).array()) / (1 + beta_c * df0.col(0).array());
                } else{
                    R.col(0) = R.col(0).array() * ((beta_c - beta_p) * df0.col(0)).array().exp();
                    R.col(1) = R.col(1).array() * ((beta_c - beta_p) * df0.col(0)).array().exp();
                    R.col(2) = R.col(2).array() * ((beta_c - beta_p) * df0.col(0)).array().exp();
                }
            }else if (modelform=="GM"){
                throw invalid_argument( "GM isn't implemented" );
            } else {
                throw invalid_argument( "GM isn't implemented" );
            }
            Ll = 0.0;
            Lld = 0.0;
            Lldd = 0.0;
            //
            temp = (PyrC.col(1).array() * (PyrC.col(0).array() * R.col(0).array()).array().log()).array() - (PyrC.col(0).array() * R.col(0).array());
            Ll =  (temp.array().isFinite()).select(temp,0).sum();
            temp = (PyrC.col(1).array() * R.col(1).array() / R.col(0).array()).array() - (PyrC.col(0).array() * R.col(1).array());
            Lld = (temp.array().isFinite()).select(temp,0).sum();
            temp = (PyrC.col(1).array() * (R.col(0).pow(-1).array() * R.col(2).array() - (R.col(1).array() * R.col(0).pow(-1).array()).pow(2) )).array() - (PyrC.col(0).array() * R.col(2).array());
            Lldd = (temp.array().isFinite()).select(temp,0).sum();
            //
            beta_b = 0.0;
        //
        // --------------------------
        // updates the change in initial guess
        // --------------------------
        if ((abs(lr*Lld/Lldd)<dbeta_kill)&&(halves<halfmax)&&(Lldd!=0.0)&&(abs(Lld)>deriv_epsilon)){
            if (Ll > Ll_best){ //if it improves, it trys the derivative again
                if ((Lld > 0)&&(Lldd<0)){
                    beta_b = - lr * Lld / Lldd;
                } else if ((Lld<0)&&(Lldd>0)){
                    beta_b = lr * Lld / Lldd;
                } else if ((Lld > 0)&&(Lldd>0)){
                    beta_b = lr * Lld / Lldd;
                } else if ((Lld<0)&&(Lldd<0)){
                    beta_b = - lr * Lld / Lldd;
                } else {
                    beta_b= 0.0;
                }
            } else{ //if it doesn't improve, it takes a half step
                if (dbeta!=0){
                    beta_b=-0.5 * beta_prev;
                }else{
                    beta_b=0.0;
                }
            }
            beta_b = beta_b * pow(.9,halves);
            if (abs(beta_b)>dbeta_max){
                beta_b = dbeta_max * sign(beta_b);
            }
            beta_prev = beta_b;
            dbeta = dbeta + beta_b;
        }
        // --------------------------
        // updates best guess
        // --------------------------
        if (Ll > Ll_best){
            Ll_best = Ll;
            beta_best = beta_c;
        }
            dev = pow((PyrC.col(1).array() - PyrC.col(0).array() * R.col(0).array()).pow(2).sum(),.5)/totem;
            end_point = high_resolution_clock::now();
            end = time_point_cast<microseconds>(end_point).time_since_epoch().count();
            cout<<"df101,"<<(end-start)<<","<<halves<<","<< Ll<<","<<Lld<<","<<Lldd<<","<<Ll_best<<","<<beta_b<<","<<dev<<endl;
            beta_p = beta_c;
        }
        halves = 0;
        //
        beta_0[0] = beta_best;
        beta_c = beta_best;
        //
        //Stores how much improvement it made, uses to randomly sample better covariates, always positive non-zero
        //
        Improve[i] = Ll_best - Improve[i];
        //
        if (modelform=="A"){
            if (tform == "lin"){
                R.col(0) = R.col(0).array() + (beta_c - beta_p) * df0.col(0).array();
            } else if (tform == "plin"){
                R.col(0) = R.col(0).array() + (beta_c - beta_p) * df0.col(0).array();
            } else{
                R.col(0) = R.col(0).array() + (beta_c * df0.col(0)).array().exp() - (beta_p * df0.col(0)).array().exp();
                R.col(1) = R.col(1).array() * ((beta_c - beta_p) * df0.col(0)).array().exp();
                R.col(2) = R.col(2).array() * ((beta_c - beta_p) * df0.col(0)).array().exp();
            }
        } else if (modelform=="PA"){
            if (tform == "lin"){
                if (i==fir){
                    R.col(0) = R.col(0) * beta_c / beta_p;
                    T0 = T0 * beta_c / beta_p;
                } else {
                    R.col(0) = R.col(0).array() + (beta_c - beta_p) * df0.col(0).array() * T0.array();
                    Te = Te.array() + (beta_c - beta_p) * df0.col(0).array();
                }
            } else if (tform == "plin"){
                if (i==fir){
                    R.col(0) = (R.col(0) - Te) * beta_c / beta_p + Te;
                    T0 = (T0.array() - 1) * beta_c / beta_p + 1;
                } else {
                    R.col(0) = R.col(0).array() + (beta_c - beta_p) * df0.col(0).array() * T0.array();
                    Te = Te.array() + (beta_c - beta_p) * df0.col(0).array();
                }
            } else{
                if (i==fir){
                    R.col(0) = R.col(0).array() * ((beta_c - beta_p) * df0.col(0)).array().exp();
                    R.col(1) = R.col(1).array() * ((beta_c - beta_p) * df0.col(0)).array().exp();
                    R.col(2) = R.col(2).array() * ((beta_c - beta_p) * df0.col(0)).array().exp();
                    T0 = T0 * ((beta_c - beta_p) * df0.col(0)).exp();
                } else {
                    R.col(0) = R.col(0).array() + ((beta_c * df0.col(0)).array().exp() - (beta_p * df0.col(0)).array().exp()).array() * T0.array();
                    R.col(1) = R.col(1).array() * ((beta_c - beta_p) * df0.col(0)).array().exp();
                    R.col(2) = R.col(2).array() * ((beta_c - beta_p) * df0.col(0)).array().exp();
                    Te = Te.array() + ((beta_c * df0.col(0)).array().exp() - (beta_p * df0.col(0)).array().exp()).array();
                }
            }
        } else if (modelform=="PAE"){
            if (tform == "lin"){
                if (i==fir){
                    R.col(0) = R.col(0) * beta_c / beta_p;
                    T0 = T0 * beta_c / beta_p;
                } else {
                    R.col(0) = R.col(0).array() + (beta_c - beta_p) * df0.col(0).array() * T0.array();
                    Te = Te.array() + (beta_c - beta_p) * df0.col(0).array();
                }
            } else if (tform == "plin"){
                if (i==fir){
                    R.col(0) = (R.col(0).array() - Te.array()) * beta_c / beta_p + Te.array();
                    T0 = (T0.array() - 1) * beta_c / beta_p + 1;
                } else {
                    R.col(0) = R.col(0).array() + (beta_c - beta_p) * df0.col(0).array() * T0.array();
                    Te = Te.array() + (beta_c - beta_p) * df0.col(0).array();
                }
            } else{
                if (i==fir){
                    R.col(0) = R.col(0).array() * ((beta_c - beta_p) * df0.col(0)).array().exp();
                    R.col(1) = R.col(1).array() * ((beta_c - beta_p) * df0.col(0)).array().exp();
                    R.col(2) = R.col(2).array() * ((beta_c - beta_p) * df0.col(0)).array().exp();
                    T0 = T0.array() * ((beta_c - beta_p) * df0.col(0)).array().exp();
                } else {
                    R.col(0) = R.col(0).array() + ((beta_c * df0.col(0)).array().exp() - (beta_p * df0.col(0)).array().exp()).array() * T0.array();
                    R.col(1) = R.col(1).array() * ((beta_c - beta_p) * df0.col(0)).array().exp();
                    R.col(2) = R.col(2).array() * ((beta_c - beta_p) * df0.col(0)).array().exp();
                    Te = Te.array() + ((beta_c * df0.col(0)).exp() - (beta_p * df0.col(0)).exp()).array();
                }
            }
        }else if (modelform=="M"){
            if (tform == "lin"){
                R.col(0) = R.col(0) * beta_c / beta_p;
                R.col(1) = R.col(0) / beta_c;
            } else if (tform == "plin"){
                R.col(0) = R.col(0).array() * (1 + beta_c * df0.col(0).array()) / (1 + beta_p * df0.col(0).array());
                R.col(1) = R.col(1).array() * (df0.col(0).array()) / (1 + beta_c * df0.col(0).array());
            } else{
                R.col(0) = R.col(0).array() * ((beta_c - beta_p) * df0.col(0)).array().exp();
                R.col(1) = R.col(1).array() * ((beta_c - beta_p) * df0.col(0)).array().exp();
                R.col(2) = R.col(2).array() * ((beta_c - beta_p) * df0.col(0)).array().exp();
            }
        }else if (modelform=="GM"){
            throw invalid_argument( "GM isn't implemented" );
        } else {
            throw invalid_argument( "Model isn't implemented" );
        }
        //
        if (tform == "lin"){
            beta_lin[ind0] = beta_0[0];
        } else if (tform == "loglin"){
            beta_loglin[ind0] = beta_0[0];
        } else {
            beta_plin[ind0] = beta_0[0];
        }
        dev = pow((PyrC.col(1).array() - PyrC.col(0).array() * R.col(0).array()).pow(2).sum(),.5)/totem;
        end_point = high_resolution_clock::now();
        end = time_point_cast<microseconds>(end_point).time_since_epoch().count();
        cout<<"df101,"<<(end-start)<<","<<halves<<","<< Ll<<","<<Lld<<","<<Lldd<<","<<Ll_best<<",0,"<<dev<<endl;
        if ((iteration > totalnum * 2)&&(iteration % (totalnum*2) == 0)){ //checks if it should check for convergence
            Ll_check[0] = Ll_check[1];
            Ll_check[1] = Ll; //the values to compare are the current against the first/last time it checked
            //
            if (abs(Ll_check[0] - Ll_check[1]) / abs(Ll_check[0]) < epsilon){
                iteration = maxiter; //if it meets the threshold it just uses the iteration to exit
            }
        }
    }
    // --------------------------
    // returns the final set of parameters
    // --------------------------
    if (include_bool[0]==1){
        if (include_bool[1]==1){
            if (include_bool[2]==1){
                res << beta_lin, beta_loglin, beta_plin;
            } else{
                res << beta_lin, beta_loglin;
            }
        } else{
            if (include_bool[2]==1){
                res << beta_lin, beta_plin;
            } else{
                res << beta_lin;
            }
        }
    } else{
        if (include_bool[1]==1){
            if (include_bool[2]==1){
                res << beta_loglin, beta_plin;
            } else{
                res << beta_loglin;
            }
        } else{
            if (include_bool[2]==1){
                res << beta_plin;
            } else{
                throw invalid_argument( "There must be covariates" );
            }
        }
    }
    return res;
}

MatrixXd Residual_AMFIT( MatrixXd PyrC, VectorXd beta_linT,VectorXd beta_loglinT,VectorXd beta_plinT,MatrixXd df_lin,MatrixXd df_loglin,MatrixXd df_plin,string modelform, NumericVector include_bool){
    VectorXd beta_0(1); //Functionally the same as the log-likelihood code, just finds the residuals instead
    VectorXd df0(df_lin.rows(), 1);
    string tform;
    int totalnum = 0;
    int ind0 = 0;
    //
    cout.precision(10);
    int nthreads = Eigen::nbThreads();
    //
    VectorXd beta_lin;
    VectorXd beta_loglin;
    VectorXd beta_plin;
    //
    if (include_bool[0]==1){
        beta_lin = beta_linT.tail(beta_linT.size()-1);
    }
    if (include_bool[1]==1){
        beta_loglin = beta_loglinT.tail(beta_loglinT.size()-1);
    }
    if (include_bool[2]==1){
        beta_plin = beta_plinT.tail(beta_plinT.size()-1);
    }
    //
    if (include_bool[0]==1){
        totalnum = totalnum + beta_lin.size();
    }
    if (include_bool[1]==1){
        totalnum = totalnum + beta_loglin.size();
    }
    if (include_bool[2]==1){
        totalnum = totalnum + beta_plin.size();
    }
    //
    double totem = PyrC.rows();
    //
    //
    //
    srand (time(NULL));
    //
    using namespace std::chrono;
    cout << "START_NEW" << endl;
    time_point<high_resolution_clock> start_point, end_point;
    start_point = high_resolution_clock::now();
    auto start = time_point_cast<microseconds>(start_point).time_since_epoch().count();
    end_point = high_resolution_clock::now();
    auto end = time_point_cast<microseconds>(end_point).time_since_epoch().count();
    // ---------------------------------------------
    // To Start, needs to seperate the derivative term
    // ---------------------------------------------
    if (include_bool[0]==1){
        if (ind0 < beta_lin.size()){
            // one exists and is one
            beta_0 = beta_lin.segment(ind0,1);
            df0 = df_lin.col(ind0);
            tform = "lin";
            //
            VectorXd temp1(beta_lin.size()-1);
            MatrixXd temp2(df_lin.rows(), df_lin.cols()-1);
            temp1 << beta_lin.head(ind0), beta_lin.tail(beta_lin.size()-1-ind0);
            temp2 << df_lin.leftCols(ind0), df_lin.rightCols(df_lin.cols()-1-ind0);
            beta_lin.resize(beta_lin.size()-1);
            df_lin.resize(df_lin.rows(), df_lin.cols()-1);
            beta_lin=temp1;
            df_lin=temp2;
        } else {
            //one exists and its not one
            ind0 = ind0 - beta_lin.size();
            if (include_bool[1]==1){
                if (ind0 < beta_loglin.size()){
                    //one and two exists and is two
                    beta_0 = beta_loglin.segment(ind0,1);
                    df0 = df_loglin.col(ind0);
                    tform = "loglin";
                    //
                    VectorXd temp1(beta_loglin.size()-1);
                    MatrixXd temp2(df_loglin.rows(), df_loglin.cols()-1);
                    temp1 << beta_loglin.head(ind0), beta_loglin.tail(beta_loglin.size()-1-ind0);
                    temp2 << df_loglin.leftCols(ind0), df_loglin.rightCols(df_loglin.cols()-1-ind0);;
                    beta_loglin.resize(beta_loglin.size()-1);
                    df_loglin.resize(df_loglin.rows(), df_loglin.cols()-1);
                    beta_loglin=temp1;
                    df_loglin=temp2;
                } else{
                    //one exists, two does, must be three
                    if (include_bool[2]!=1){
                        throw invalid_argument( "Are all three used?" );
                    }
                    ind0 = ind0 - beta_loglin.size();
                    beta_0 = beta_plin.segment(ind0,1);
                    df0 = df_plin.col(ind0);
                    tform = "plin";
                    //
                    VectorXd temp1(beta_plin.size()-1);
                    MatrixXd temp2(df_plin.rows(), df_plin.cols()-1);
                    temp1 << beta_plin.head(ind0), beta_plin.tail(beta_plin.size()-1-ind0);
                    temp2 << df_plin.leftCols(ind0), df_plin.rightCols(df_plin.cols()-1-ind0);
                    beta_plin.resize(beta_plin.size()-1);
                    df_plin.resize(df_plin.rows(), df_plin.cols()-1);
                    beta_plin=temp1;
                    df_plin=temp2;
                }
            } else{
                //one exists, and two doesn't exist, must be three
                if (include_bool[2]!=1){
                    throw invalid_argument( "Are all first and third used?" );
                }
                beta_0 = beta_plin.segment(ind0,1);
                df0 = df_plin.col(ind0);
                tform = "plin";
                //
                VectorXd temp1(beta_plin.size()-1);
                MatrixXd temp2(df_plin.rows(), df_plin.cols()-1);
                temp1 << beta_plin.head(ind0), beta_plin.tail(beta_plin.size()-1-ind0);
                temp2 << df_plin.leftCols(ind0), df_plin.rightCols(df_plin.cols()-1-ind0);
                beta_plin.resize(beta_plin.size()-1);
                df_plin.resize(df_plin.rows(), df_plin.cols()-1);
                beta_plin=temp1;
                df_plin=temp2;
            }
        }
    }else{
        //one doesn't exist
        if (include_bool[1]==1){
            if (ind0 < beta_loglin.size()){
                //one doesn't exist and two exists and is two
                beta_0 = beta_loglin.segment(ind0,1);
//                cout << ind0 << ", " << beta_0 << ", " << beta_loglin.transpose() << endl;
                df0 = df_loglin.col(ind0);
                tform = "loglin";
                //
                VectorXd temp1(beta_loglin.size()-1);
                MatrixXd temp2(df_loglin.rows(), df_loglin.cols()-1);
                temp1 << beta_loglin.head(ind0), beta_loglin.tail(beta_loglin.size()-1-ind0);
                temp2 << df_loglin.leftCols(ind0), df_loglin.rightCols(df_loglin.cols()-1-ind0);;
                beta_loglin.resize(beta_loglin.size()-1);
                df_loglin.resize(df_loglin.rows(), df_loglin.cols()-1);
                beta_loglin=temp1;
                df_loglin=temp2;
            } else{
                //one doesn't exist, two does, must be three
                if (include_bool[2]!=1){
                    throw invalid_argument( "Are all three used?" );
                }
                ind0 = ind0 - beta_loglin.size();
                beta_0 = beta_plin.segment(ind0,1);
                df0 = df_plin.col(ind0);
                tform = "plin";
                //
                VectorXd temp1(beta_plin.size()-1);
                MatrixXd temp2(df_plin.rows(), df_plin.cols()-1);
                temp1 << beta_plin.head(ind0), beta_plin.tail(beta_plin.size()-1-ind0);
                temp2 << df_plin.leftCols(ind0), df_plin.rightCols(df_plin.cols()-1-ind0);
                beta_plin.resize(beta_plin.size()-1);
                df_plin.resize(df_plin.rows(), df_plin.cols()-1);
                beta_plin=temp1;
                df_plin=temp2;
            }
        } else{
            //one doesn't exist, and two doesn't exist, must be three
            if (include_bool[2]!=1){
                throw invalid_argument( "Are all first and third used?" );
            }
            beta_0 = beta_plin.segment(ind0,1);
            df0 = df_plin.col(ind0);
            tform = "plin";
            //
            VectorXd temp1(beta_plin.size()-1);
            MatrixXd temp2(df_plin.rows(), df_plin.cols()-1);
            temp1 << beta_plin.head(ind0), beta_plin.tail(beta_plin.size()-1-ind0);
            temp2 << df_plin.leftCols(ind0), df_plin.rightCols(df_plin.cols()-1-ind0);
            beta_plin.resize(beta_plin.size()-1);
            df_plin.resize(df_plin.rows(), df_plin.cols()-1);
            beta_plin=temp1;
            df_plin=temp2;
        }
    }
    //
    end_point = high_resolution_clock::now();
    end = time_point_cast<microseconds>(end_point).time_since_epoch().count();
    cout<<"df101,"<<(end-start)<<",0,0,0,0,0,0"<<endl;
    //
    MatrixXd T0 = MatrixXd::Zero(df0.rows(), 1);
    MatrixXd Te = MatrixXd::Zero(df0.rows(), 1);
    MatrixXd R = MatrixXd::Zero(df0.rows(), 1);
    if ((modelform=="A")||(modelform=="PA")||(modelform=="PAE")){
        Te = Te.array() * 0;
        //
        T0 = df0 * beta_0;
        if (tform=="lin") {
            T0 = T0;
        } else if (tform=="loglin") {
            T0 = T0.array().exp();;
        } else {
            T0 = 1 + T0.array();
        }
        //
        if (include_bool[0]==1){
            Te = Te.array() + (df_lin.array().rowwise() * beta_lin.transpose().array()).matrix().rowwise().sum().array();
        }
        if (include_bool[1]==1){
            Te = Te.array() + (df_loglin.array().rowwise() * beta_loglin.transpose().array()).array().exp().rowwise().sum().array();
        }
        if (include_bool[2]==1){
            Te = Te.array() + (1 + df_plin.array().rowwise() * beta_plin.transpose().array()).matrix().rowwise().sum().array();
        }
        if (modelform=="A"){
            R << Te;
        } else if (modelform=="PA"){
            R << T0.array() * Te.array();
        } else if (modelform=="PAE"){
            Te = Te.array() + 1;
            R << T0.array() * Te.array();
        }
    }else if (modelform=="M"){
        Te = Te.array() * 0 + 1;
        //
        T0 = df0 * beta_0;
        if (tform=="lin") {
            T0 = T0;
        } else if (tform=="loglin") {
            T0 = T0.array().exp();
        } else {
            T0 = 1 + T0.array();
        }
        //
        if (include_bool[0]==1){
            Te = Te.array() * (df_lin.array().rowwise() * beta_lin.transpose().array()).matrix().rowwise().prod().array();
        }
        if (include_bool[1]==1){
            Te = Te.array() * (df_loglin.array().rowwise() * beta_loglin.transpose().array()).array().rowwise().sum().array().exp().array();
        }
        if (include_bool[2]==1){
            Te = Te.array() * (1 + df_plin.array().rowwise() * beta_plin.transpose().array()).matrix().rowwise().prod().array();
        }
        R << T0.array() * Te.array();
    } else if (modelform=="GM"){
        throw invalid_argument( "GM isn't implemented" );
    }
    R =(R.array().isFinite()).select(R,0);
    //
    ArrayXd Residuals(R.rows(),2);
    //
    ArrayXd guess = PyrC.col(0).array() * R.col(0).array();
    ArrayXd dif = PyrC.col(1).array() - guess.array();
    // Finds the deviance and pearson residuals respectively
    Residuals.col(0) = pow(2,.5) * (dif.abs() * dif.pow(-1)) * (PyrC.col(1).array() * (PyrC.col(1).array() * guess.array().pow(-1)).array().log() - dif.array()).array().pow(.5);
    Residuals.col(1) = dif.array() * guess.array().pow(-.5);
    return Residuals;
}
