#include <RcppEigen.h>
//#include <RcppParallel.h>
#include <omp.h>
#include "AMFIT_MODEL.h"
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <chrono>
#include <random>
#include<ctime>


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
List amfit_transition(NumericVector a_lin,NumericVector a_loglin,NumericVector a_plin, NumericMatrix x_lin, NumericMatrix x_loglin, NumericMatrix x_plin,int fir,string modelform,NumericMatrix dfe, NumericVector include_bool, double lr, int maxiter, int halfmax, double epsilon, double dbeta_max, double deriv_epsilon, int batch_size){
    //----------------------------------------------------------------------------------------------------------------//
    Map<VectorXd> beta_lin(as<Map<VectorXd> >(a_lin));
    Map<VectorXd> beta_loglin(as<Map<VectorXd> >(a_loglin));
    Map<VectorXd> beta_plin(as<Map<VectorXd> >(a_plin));
    const Map<MatrixXd> df_lin(as<Map<MatrixXd> >(x_lin));
    const Map<MatrixXd> df_loglin(as<Map<MatrixXd> >(x_loglin));
    const Map<MatrixXd> df_plin(as<Map<MatrixXd> >(x_plin));
    const Map<MatrixXd> PyrC(as<Map<MatrixXd> >(dfe));
    // Converts from Rcpp types to efficient Eigen types
    //----------------------------------------------------------------------------------------------------------------//
    List res = LogLik_AMFIT(PyrC,beta_lin,beta_loglin,beta_plin,df_lin,df_loglin,df_plin,fir,modelform, include_bool, lr, maxiter, halfmax, epsilon, dbeta_max, deriv_epsilon, batch_size);
    //----------------------------------------------------------------------------------------------------------------//
    return wrap(res);
}

List LogLik_AMFIT( MatrixXd PyrC, VectorXd beta_linT,VectorXd beta_loglinT,VectorXd beta_plinT,MatrixXd df_lin,MatrixXd df_loglin,MatrixXd df_plin,int fir,string modelform, NumericVector include_bool, double lr, int maxiter, int halfmax, double epsilon, double dbeta_max, double deriv_epsilon, int batch_size){
    //    cout << "start Risk" << endl;
    //
    srand (time(NULL));
    //

    using namespace std::chrono;
    cout << "START_NEW" << endl;
    time_point<system_clock> start_point, end_point;
    start_point = system_clock::now();
    auto start = time_point_cast<microseconds>(start_point).time_since_epoch().count();
    end_point = system_clock::now();
    auto end = time_point_cast<microseconds>(end_point).time_since_epoch().count(); //The time duration is tracked
    //
    auto gibtime = system_clock::to_time_t(system_clock::now());
    cout << ctime(&gibtime) << endl;
    //
    int totalnum = 0;
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
    double totem = PyrC.rows(); //precalculates how many rows
    //
    vector <int> batch_cols;
    int i_temp=0;
    int j_temp=0;
    batch_cols.push_back(0);
    while (i_temp<totem-batch_size){
        if (j_temp==batch_size){
            batch_cols.push_back(i_temp);
            j_temp=0;
        }
        i_temp++;
        j_temp++;
    }
    if (totem-i_temp>batch_size/2){
        batch_cols.push_back(i_temp);
    }
    batch_cols.push_back(totem-1);
    //
    //
    end_point = system_clock::now();
    end = time_point_cast<microseconds>(end_point).time_since_epoch().count();
    cout<<"df99,"<<(end-start)<<",Batches"<<endl;
    gibtime = system_clock::to_time_t(system_clock::now());
    cout << ctime(&gibtime) << endl;
    // ---------------------------------------------
    // To Start, needs to seperate the derivative term
    // ---------------------------------------------
    VectorXd beta_0(totalnum);
    MatrixXd df0 = MatrixXd::Zero(df_lin.rows(), totalnum); // stores memory for the derivative term parameters and columns
    MatrixXd T0 = MatrixXd::Zero(df0.rows(), totalnum); //preallocates matrix for Derivative column terms
    MatrixXd Td0 = MatrixXd::Zero(df0.rows(), totalnum); //preallocates matrix for Derivative column terms
    MatrixXd Tdd0 = MatrixXd::Zero(df0.rows(), totalnum); //preallocates matrix for Derivative column terms
    #pragma omp parallel for num_threads(nthreads)
    for (int ij=0;ij<totalnum;ij++){
        int ind0 = ij;
        if (include_bool[0]==1){
            if (ind0 < beta_lin.size()){
                // one exists and is one
                beta_0[ij] = beta_lin[ind0];
                df0.col(ij) = df_lin.col(ind0);
                tform[ij] = "lin";
                //
            } else {
                //one exists and its not one
                ind0 = ind0 - beta_lin.size();
                if (include_bool[1]==1){
                    if (ind0 < beta_loglin.size()){
                        //one and two exists and is two
                        beta_0[ij] = beta_loglin[ind0];
                        df0.col(ij) = df_loglin.col(ind0);
                        tform[ij] = "loglin";
                        //
                    } else{
                        //one exists, two does, must be three
                        if (include_bool[2]!=1){
                            throw invalid_argument( "Are all three used?" );
                        }
                        ind0 = ind0 - beta_loglin.size();
                        beta_0[ij] = beta_plin[ind0];
                        df0.col(ij) = df_plin.col(ind0);
                        tform[ij] = "plin";
                        //
                    }
                } else{
                    //one exists, and two doesn't exist, must be three
                    if (include_bool[2]!=1){
                        throw invalid_argument( "Are all first and third used?" );
                    }
                    beta_0[ij] = beta_plin[ind0];
                    df0.col(ij) = df_plin.col(ind0);
                    tform[ij] = "plin";
                    //
                }
            }
        }else{
            //one doesn't exist
            if (include_bool[1]==1){
                if (ind0 < beta_loglin.size()){
                    //one doesn't exist and two exists and is two
                    beta_0[ij] = beta_loglin[ind0];
    //                cout << ind0 << ", " << beta_0 << ", " << beta_loglin.transpose() << endl;
                    df0.col(ij) = df_loglin.col(ind0);
                    tform[ij] = "loglin";
                    //
                } else{
                    //one doesn't exist, two does, must be three
                    if (include_bool[2]!=1){
                        throw invalid_argument( "Are all three used?" );
                    }
                    ind0 = ind0 - beta_loglin.size();
                    beta_0[ij] = beta_plin[ind0];
                    df0.col(ij) = df_plin.col(ind0);
                    tform[ij] = "plin";
                    //
                }
            } else{
                //one doesn't exist, and two doesn't exist, must be three
                if (include_bool[2]!=1){
                    throw invalid_argument( "Are all first and third used?" );
                }
                beta_0[ij] = beta_plin[ind0];
                df0.col(ij) = df_plin.col(ind0);
                tform[ij] = "plin";
                //
            }
        }
//        cout << df0.col(ij).array().transpose() << endl;
        T0.col(ij) = (df0.col(ij).array() * beta_0[ij]).matrix();
        if (tform[ij]=="lin") {
            Td0.col(ij) = df0.col(ij);
        } else if (tform[ij]=="loglin") {
            T0.col(ij) = T0.col(ij).array().exp();
            Td0.col(ij) = df0.col(ij).array() * T0.col(ij).array();
            Tdd0.col(ij) = df0.col(ij).array() * Td0.col(ij).array();
        } else {
            T0.col(ij) = 1 + T0.col(ij).array();
            Td0.col(ij) = df0.col(ij);
        }
    }
    //
    end_point = system_clock::now();
    end = time_point_cast<microseconds>(end_point).time_since_epoch().count();
    cout<<"df99,"<<(end-start)<<",Terms"<<endl;
    gibtime = system_clock::to_time_t(system_clock::now());
    cout << ctime(&gibtime) << endl;
    //
    MatrixXd Te = MatrixXd::Zero(df0.rows(), 1); //preallocates matrix for non-Derivative column terms
    MatrixXd R = MatrixXd::Zero(df0.rows(), 1); //preallocates matrix for Risks and derivatives
    MatrixXd Rd = MatrixXd::Zero(df0.rows(), totalnum); //preallocates matrix for Risks and derivatives
    MatrixXd Rdd = MatrixXd::Zero(df0.rows(), pow(totalnum,2));
    if ((modelform=="A")||(modelform=="PA")||(modelform=="PAE")){ //same process used for all of the additive type models
        Te = T0.array().rowwise().sum();
        // computes intial risk and derivatives
        if (modelform=="A"){
            R << Te.array();
            Rd << Td0.array();
            #pragma omp parallel for num_threads(nthreads)
            for (int ij=0;ij<totalnum;ij++){
                Rdd.col(ij*totalnum+ij) = Tdd0.col(ij);
            }
        } else if ((modelform=="PAE")||(modelform=="PA")){
            Te = Te.array() - T0.col(fir).array();
            if (modelform=="PAE"){
                Te = Te.array() + 1;
            }
            R << T0.col(fir).array() * Te.array();
            Rd << Td0.array() * T0.col(fir).array();//, Td0.col(0).array() * Te.array(), Td0.col(1).array() * Te.array();
            Rd.col(fir) = Td0.col(fir).array() * Te.array();
            #pragma omp parallel for schedule(dynamic) num_threads(nthreads) collapse(2)
            for (int bl=0;bl<batch_cols.size()-1;bl++){
                for (int ijk=0;ijk<totalnum*(totalnum+1)/2;ijk++){
                    int ij = 0;
                    int jk = ijk;
                    while (jk>ij){
                        ij++;
                        jk-=ij;
                    }
                    if (ij==jk){
                        if (ij==fir){
                            Rdd.block(batch_cols[bl],ij*totalnum+jk,batch_cols[bl+1]-batch_cols[bl],1) = Tdd0.block(batch_cols[bl],ij,batch_cols[bl+1]-batch_cols[bl],1).array() * Te.block(batch_cols[bl],1,batch_cols[bl+1]-batch_cols[bl],1).array();
                        } else {
                            Rdd.block(batch_cols[bl],ij*totalnum+jk,batch_cols[bl+1]-batch_cols[bl],1) = Tdd0.block(batch_cols[bl],ij,batch_cols[bl+1]-batch_cols[bl],1).array() * T0.block(batch_cols[bl],fir,batch_cols[bl+1]-batch_cols[bl],1).array();
                        }
                    } else {
                        if ((ij==fir)||(jk==fir)){
                            Rdd.block(batch_cols[bl],ij*totalnum+jk,batch_cols[bl+1]-batch_cols[bl],1) = Td0.block(batch_cols[bl],ij,batch_cols[bl+1]-batch_cols[bl],1).array() * Td0.block(batch_cols[bl],jk,batch_cols[bl+1]-batch_cols[bl],1).array();
                            Rdd.block(batch_cols[bl],ij*totalnum+jk,batch_cols[bl+1]-batch_cols[bl],1) = Td0.block(batch_cols[bl],ij,batch_cols[bl+1]-batch_cols[bl],1).array() * Td0.block(batch_cols[bl],jk,batch_cols[bl+1]-batch_cols[bl],1).array();
                        }
                    }
                }
            }
        }
    }else if (modelform=="M"){
        Te = Te.array() * 0 + 1; //verifies the initial term product is 1
        //
        Te = T0.array().rowwise().prod();
        // computes intial risk and derivatives
        R << Te.array();
        Rd = T0.array().pow(-1).array() * Te.colwise().replicate(totalnum).array();
        Rd = Rd.array() * Td0.array();
        #pragma omp parallel for schedule(dynamic) num_threads(nthreads) collapse(2)
        for (int bl=0;bl<batch_cols.size()-1;bl++){
            for (int ijk=0;ijk<totalnum*(totalnum+1)/2;ijk++){
                int ij = 0;
                int jk = ijk;
                while (jk>ij){
                    ij++;
                    jk-=ij;
                }
                if (ij==jk){
                    Rdd.block(batch_cols[bl],ij*totalnum+jk,batch_cols[bl+1]-batch_cols[bl],1) = Tdd0.block(batch_cols[bl],jk,batch_cols[bl+1]-batch_cols[bl],1).array() * Td0.block(batch_cols[bl],jk,batch_cols[bl+1]-batch_cols[bl],1).array().pow(-1).array() * Rd.block(batch_cols[bl],ij,batch_cols[bl+1]-batch_cols[bl],1).array();
                } else {
                    Rdd.block(batch_cols[bl],ij*totalnum+jk,batch_cols[bl+1]-batch_cols[bl],1) = Td0.block(batch_cols[bl],jk,batch_cols[bl+1]-batch_cols[bl],1).array() * T0.block(batch_cols[bl],jk,batch_cols[bl+1]-batch_cols[bl],1).array().pow(-1).array() * Rd.block(batch_cols[bl],ij,batch_cols[bl+1]-batch_cols[bl],1).array();
                    Rdd.block(batch_cols[bl],jk*totalnum+ij,batch_cols[bl+1]-batch_cols[bl],1) = Td0.block(batch_cols[bl],ij,batch_cols[bl+1]-batch_cols[bl],1).array() * T0.block(batch_cols[bl],ij,batch_cols[bl+1]-batch_cols[bl],1).array().pow(-1).array() * Rd.block(batch_cols[bl],jk,batch_cols[bl+1]-batch_cols[bl],1).array();
                }
            }
        }
    } else if (modelform=="GM"){
        //currently isn't implemented, it can be calculated but not optimized the same way
        throw invalid_argument( "GM isn't implemented" );
    } else {
        throw invalid_argument( "Model isn't implemented" );
    }
    R = (R.array().isFinite()).select(R,0);
    Rd = (Rd.array().isFinite()).select(Rd,0);
    Rdd = (Rdd.array().isFinite()).select(Rdd,0);
    //
    //
    end_point = system_clock::now();
    end = time_point_cast<microseconds>(end_point).time_since_epoch().count();
    cout<<"df100 "<<(end-start)<<" "<<0<<" "<<0<<" "<<-1<<",Prep_R"<<endl;
    gibtime = system_clock::to_time_t(system_clock::now());
    cout << ctime(&gibtime) << endl;
    vector<double> Ll(totalnum,1.0);
    vector<double> Lld(totalnum,1.0);
    vector<double> Lldd(pow(totalnum,2),1.0);
    //
    MatrixXd temp(Rd.rows(),Rd.cols());
    VectorXd tempd(1,totalnum);
    
    temp = (PyrC.col(1).array() * (PyrC.col(0).array() * R.col(0).array()).array().log()).array() - (PyrC.col(0).array() * R.col(0).array());
    fill(Ll.begin(), Ll.end(), (temp.array().isFinite()).select(temp,0).sum());

    temp = (PyrC.col(1).replicate(1,totalnum).array() * Rd.array() * R.replicate(1,totalnum).array().pow(-1).array()).array() - (PyrC.col(0).replicate(1,totalnum).array() * Rd.array());
    tempd << (temp.array().isFinite()).select(temp,0).colwise().sum();
    VectorXd::Map(&Lld[0], totalnum) = tempd;
    #pragma omp parallel for schedule(dynamic) num_threads(nthreads)
    for (int ijk=0;ijk<totalnum*(totalnum+1)/2;ijk++){//totalnum*(totalnum+1)/2
        int ij = 0;
        int jk = ijk;
        while (jk>ij){
            ij++;
            jk-=ij;
        }
        VectorXd temp(Rdd.rows(),Rdd.cols());
        temp = (PyrC.col(1).array() * (R.col(0).pow(-1).array() * Rdd.col(ij*totalnum+jk).array() - (Rd.col(ij).array() * Rd.col(jk).array() * R.col(0).pow(-1).array()).pow(2) )).array() - (PyrC.col(0).array() * Rdd.col(ij*totalnum+jk).array());
        Lldd[ij*totalnum+jk] = (temp.array().isFinite()).select(temp,0).sum();
        if (ij!=jk){
            Lldd[jk*totalnum+ij] = (temp.array().isFinite()).select(temp,0).sum();
        }
    }
    double dev = pow((PyrC.col(1).array() - PyrC.col(0).array() * R.col(0).array()).pow(2).sum(),.5)/totem;
    // --------------------------
    // stores the value to compare to in the iteration step
    // --------------------------
    end_point = system_clock::now();
    end = time_point_cast<microseconds>(end_point).time_since_epoch().count();
    cout<<"df100 "<<(end-start)<<" "<<0<<" "<<0<<" "<<0<<",Calc"<<endl;
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
    double dbeta=0;
    double beta_p = 0.0;
    double beta_c = 0.0;
    double beta_a = 0.0;
    double beta_best = 0.0;
    double Ll_best = 0.0;
    int halves = 0;
    int ind0 = fir;
    int i = ind0;
    double beta_b=dbeta;
    double beta_prev=dbeta;
    int iteration=0;
    //
    while (iteration < maxiter){
        iteration++;
        //
        ind0 = rand() % totalnum;
        if ((Lld[ind0] > 0)&&(Lldd[ind0*totalnum+ind0]<0)){
            dbeta = -lr * Lld[ind0] / Lldd[ind0*totalnum+ind0];// decreases left, zeros right, ratio left
        } else if ((Lld[ind0]<0)&&(Lldd[ind0*totalnum+ind0]>0)){
            dbeta = -lr * Lld[ind0] / Lldd[ind0*totalnum+ind0];// decreases right, zeros right, ratio left
        } else if ((Lld[ind0] > 0)&&(Lldd[ind0*totalnum+ind0]>0)){
            dbeta = -lr * Lld[ind0] / Lldd[ind0*totalnum+ind0];// decreases left, zeros left, ratio right
        } else if ((Lld[ind0]<0)&&(Lldd[ind0*totalnum+ind0]<0)){
            dbeta = -lr * Lld[ind0] / Lldd[ind0*totalnum+ind0];// decreases right, zeros left, ratio right
        } else {
            dbeta=0.0;
        }
        //
        if (abs(dbeta)>dbeta_max){
            dbeta = dbeta_max * sign(dbeta);
        }
        //
        beta_p = beta_0[ind0];
        beta_c = beta_p;
        beta_a = beta_0[ind0];
        beta_best = beta_p;
        Ll_best = Ll[ind0];
        // int halves = 0;
        i = ind0;
        beta_b=dbeta;
        beta_prev=dbeta;
//        cout << beta_best << ", " << Ll[ind0] << endl;
        //
        halves=0;
        while ((Ll[ind0] <= Ll_best)&&(halves<halfmax)&&(abs(dbeta) > epsilon)){
            halves++;
            beta_c = beta_a + dbeta;
            if (tform[ind0]=="lin"){
                T0.col(ind0) = T0.col(ind0).array() * (beta_c / beta_p);
            } else if (tform[ind0]=="plin"){
                T0.col(ind0) = T0.col(ind0).array() * (1 + beta_c * df0.col(ind0).array()) / (1 + beta_p * df0.col(ind0).array());
            } else {
                T0.col(ind0) = T0.col(ind0).array() * ((beta_c - beta_p) * df0.col(ind0)).array().exp().array();
                Td0.col(ind0) = Td0.col(ind0).array() * ((beta_c - beta_p) * df0.col(ind0)).array().exp().array();
                Tdd0.col(ind0) = Tdd0.col(ind0).array() * ((beta_c - beta_p) * df0.col(ind0)).array().exp().array();
            }
            if (modelform=="A"){//The risks can be updated instead of recalculated for the models/terms allowed
                if (tform[ind0] == "lin"){
                    R.col(0) = R.col(0).array() + (beta_c - beta_p) * df0.col(ind0).array();
                } else if (tform[ind0] == "plin"){
                    R.col(0) = R.col(0).array() + (beta_c - beta_p) * df0.col(ind0).array();
                } else{
                    R.col(0) = R.col(0).array() + (beta_c * df0.col(ind0)).array().exp() - (beta_p * df0.col(ind0)).array().exp();
                    Rd.col(ind0) = Rd.col(ind0).array() * ((beta_c - beta_p) * df0.col(ind0)).array().exp().array();
                    Rdd.col(ind0*totalnum+ind0) = Rdd.col(ind0*totalnum+ind0).array() * ((beta_c - beta_p) * df0.col(ind0)).array().exp().array();
                }
            } else if ((modelform=="PA")||(modelform=="PAE")){
                if (tform[ind0] == "lin"){
                    if (ind0==fir){
                        R.col(0) = R.col(0) * beta_c / beta_p;
                    } else {
                        R.col(0) = R.col(0).array() + (beta_c - beta_p) * df0.col(ind0).array() * T0.col(fir).array();
                        Te.col(0) = Te.col(0).array() + (beta_c - beta_p) * df0.col(ind0).array();
                    }
                } else if (tform[ind0] == "plin"){
                    if (ind0==fir){
                        R.col(0) = (R.col(0).array() - Te.array()) * beta_c / beta_p + Te.col(0).array();
                    } else {
                        R.col(0) = R.col(0).array() + (beta_c - beta_p) * df0.col(ind0).array() * T0.col(fir).array();
                        Te.col(0) = Te.col(0).array() + (beta_c - beta_p) * df0.col(ind0).array();
                    }
                } else{
                    if (ind0==fir){
                        R.col(0) = R.col(0).array() * ((beta_c - beta_p) * df0.col(ind0)).array().exp();
                    } else {
                        R.col(0) = R.col(0).array() + ((beta_c * df0.col(ind0)).array().exp() - (beta_p * df0.col(ind0)).array().exp()).array() * T0.col(fir).array();
                        Te.col(0) = Te.col(0).array() + ((beta_c * df0.col(ind0)).array().exp() - (beta_p * df0.col(ind0)).array().exp()).array();
                    }
                }
                #pragma omp parallel for schedule(dynamic) num_threads(nthreads) collapse(2)
                for (int bl=0;bl<batch_cols.size()-1;bl++){
                    for (int ijk=0;ijk<totalnum*(totalnum+1)/2;ijk++){
                        int ij = 0;
                        int jk = ijk;
                        while (jk>ij){
                            ij++;
                            jk-=ij;
                        }
                        if (ij==ind0){
                            if (ij==jk){
                                if (tform[ij]=="loglin"){
                                    Rd.block(batch_cols[bl],ij,batch_cols[bl+1]-batch_cols[bl],1) = Rd.block(batch_cols[bl],ij,batch_cols[bl+1]-batch_cols[bl],1).array() * ((beta_c - beta_p) * df0.block(batch_cols[bl],ind0,batch_cols[bl+1]-batch_cols[bl],1)).array().exp();
                                    Rdd.block(batch_cols[bl],ij*totalnum+jk,batch_cols[bl+1]-batch_cols[bl],1) = Rdd.block(batch_cols[bl],ij,batch_cols[bl+1]-batch_cols[bl],1).array() * ((beta_c - beta_p) * df0.block(batch_cols[bl],ind0,batch_cols[bl+1]-batch_cols[bl],1)).array().exp();
                                } else {
                                    ;
                                }
                            } else{
                                if ((ij!=fir)&&(jk!=fir)){
                                    ;
                                } else if (tform[ind0]=="loglin"){
                                    Rdd.block(batch_cols[bl],ij*totalnum+jk,batch_cols[bl+1]-batch_cols[bl],1) = Rdd.block(batch_cols[bl],ij*totalnum+jk,batch_cols[bl+1]-batch_cols[bl],1).array() * ((beta_c - beta_p) * df0.block(batch_cols[bl],ind0,batch_cols[bl+1]-batch_cols[bl],1)).array().exp();
                                    Rdd.block(batch_cols[bl],jk*totalnum+ij,batch_cols[bl+1]-batch_cols[bl],1) = Rdd.block(batch_cols[bl],jk*totalnum+ij,batch_cols[bl+1]-batch_cols[bl],1).array() * ((beta_c - beta_p) * df0.block(batch_cols[bl],ind0,batch_cols[bl+1]-batch_cols[bl],1)).array().exp();
                                }
                            }
                        } else if (jk==ind0){
                            if ((ij!=fir)&&(jk!=fir)){
                                ;
                            } else if (tform[jk]=="loglin"){
                                Rdd.block(batch_cols[bl],jk*totalnum+ij,batch_cols[bl+1]-batch_cols[bl],1) = Rdd.block(batch_cols[bl],jk*totalnum+ij,batch_cols[bl+1]-batch_cols[bl],1).array() * ((beta_c - beta_p) * df0.block(batch_cols[bl],ind0,batch_cols[bl+1]-batch_cols[bl],1)).array().exp();
                                Rdd.block(batch_cols[bl],ij*totalnum+jk,batch_cols[bl+1]-batch_cols[bl],1) = Rdd.block(batch_cols[bl],ij*totalnum+jk,batch_cols[bl+1]-batch_cols[bl],1).array() * ((beta_c - beta_p) * df0.block(batch_cols[bl],ind0,batch_cols[bl+1]-batch_cols[bl],1)).array().exp();
                            }
                        }
                    }
                }
            }else if (modelform=="M"){
                Te = Te.array() * 0 + 1; //verifies the initial term product is 1
                //
                Te = T0.array().rowwise().prod();
                // computes intial risk and derivatives
                R << Te.array();
                Rd = T0.array().pow(-1).array() * Te.colwise().replicate(totalnum).array();
                Rd = Rd.array() * Td0.array();
                #pragma omp parallel for schedule(dynamic) num_threads(nthreads) collapse(2)
                for (int bl=0;bl<batch_cols.size()-1;bl++){
                    for (int ijk=0;ijk<totalnum*(totalnum+1)/2;ijk++){
                        int ij = 0;
                        int jk = ijk;
                        while (jk>ij){
                            ij++;
                            jk-=ij;
                        }
                        if (ij==jk){
                            Rdd.block(batch_cols[bl],ij*totalnum+jk,batch_cols[bl+1]-batch_cols[bl],1) = Tdd0.block(batch_cols[bl],jk,batch_cols[bl+1]-batch_cols[bl],1).array() * Td0.block(batch_cols[bl],jk,batch_cols[bl+1]-batch_cols[bl],1).array().pow(-1).array() * Rd.block(batch_cols[bl],ij,batch_cols[bl+1]-batch_cols[bl],1).array();
                        } else {
                            Rdd.block(batch_cols[bl],ij*totalnum+jk,batch_cols[bl+1]-batch_cols[bl],1) = Td0.block(batch_cols[bl],jk,batch_cols[bl+1]-batch_cols[bl],1).array() * T0.block(batch_cols[bl],jk,batch_cols[bl+1]-batch_cols[bl],1).array().pow(-1).array() * Rd.block(batch_cols[bl],ij,batch_cols[bl+1]-batch_cols[bl],1).array();
                            Rdd.block(batch_cols[bl],jk*totalnum+ij,batch_cols[bl+1]-batch_cols[bl],1) = Td0.block(batch_cols[bl],ij,batch_cols[bl+1]-batch_cols[bl],1).array() * T0.block(batch_cols[bl],ij,batch_cols[bl+1]-batch_cols[bl],1).array().pow(-1).array() * Rd.block(batch_cols[bl],jk,batch_cols[bl+1]-batch_cols[bl],1).array();
                        }
                    }
                }
            }else if (modelform=="GM"){
                throw invalid_argument( "GM isn't implemented" );
            } else {
                throw invalid_argument( "Model isn't implemented" );
            }
            R = (R.array().isFinite()).select(R,0);
            Rd = (Rd.array().isFinite()).select(Rd,0);
            Rdd = (Rdd.array().isFinite()).select(Rdd,0);
            end_point = system_clock::now();
            end = time_point_cast<microseconds>(end_point).time_since_epoch().count();
            cout<<"df100 "<<(end-start)<<" "<<halves<<" "<<iteration<<" "<<ind0<<",Update_R"<<endl;
            gibtime = system_clock::to_time_t(system_clock::now());
            cout << ctime(&gibtime) << endl;
            fill(Ll.begin(), Ll.end(), 0.0);
            fill(Lld.begin(), Lld.end(), 0.0);
            fill(Lldd.begin(), Lldd.end(), 0.0);
            
            
            temp = (PyrC.col(1).array() * (PyrC.col(0).array() * R.col(0).array()).array().log()).array() - (PyrC.col(0).array() * R.col(0).array());
            fill(Ll.begin(), Ll.end(), (temp.array().isFinite()).select(temp,0).sum());

            temp = (PyrC.col(1).replicate(1,totalnum).array() * Rd.array() * R.replicate(1,totalnum).array().pow(-1).array()).array() - (PyrC.col(0).replicate(1,totalnum).array() * Rd.array());
            tempd << (temp.array().isFinite()).select(temp,0).colwise().sum();
            VectorXd::Map(&Lld[0], totalnum) = tempd;

            #pragma omp parallel for schedule(dynamic) num_threads(nthreads)
            for (int ijk=0;ijk<totalnum*(totalnum+1)/2;ijk++){//totalnum*(totalnum+1)/2
                int ij = 0;
                int jk = ijk;
                while (jk>ij){
                    ij++;
                    jk-=ij;
                }
                VectorXd temp(Rdd.rows(),Rdd.cols());
                temp = (PyrC.col(1).array() * (R.col(0).pow(-1).array() * Rdd.col(ij*totalnum+jk).array() - (Rd.col(ij).array() * Rd.col(jk).array() * R.col(0).pow(-1).array()).pow(2) )).array() - (PyrC.col(0).array() * Rdd.col(ij*totalnum+jk).array());
                Lldd[ij*totalnum+jk] = (temp.array().isFinite()).select(temp,0).sum();
                if (ij!=jk){
                    Lldd[jk*totalnum+ij] = (temp.array().isFinite()).select(temp,0).sum();
                }
            }
            if (Ll[ind0] <= Ll_best){
                dbeta = dbeta / 2.0;
            } else{
                beta_best = beta_c;
            }
            end_point = system_clock::now();
            end = time_point_cast<microseconds>(end_point).time_since_epoch().count();
            beta_p = beta_c;
            cout<<"df100 "<<(end-start)<<" "<<halves<<" "<<iteration<<" "<<ind0<<",Update_calc"<<endl;
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
            beta_0[ind0] = beta_c;
            cout << "df104 ";
            for (int ij=0;ij<totalnum;ij++){
                cout << beta_0[ij] << " ";
            }
            cout << " " << endl;
            beta_0[ind0] = beta_best;
//            cout << beta_best << ", " << Ll[ind0] << endl;
        }
        if (beta_best!=beta_p){
            beta_c = beta_best;
//            cout << beta_c << ":" << beta_p << endl;
            if (tform[ind0]=="lin"){
                T0.col(ind0) = T0.col(ind0).array() * (beta_c / beta_p);
            } else if (tform[ind0]=="plin"){
                T0.col(ind0) = T0.col(ind0).array() * (1 + beta_c * df0.col(ind0).array()) / (1 + beta_p * df0.col(ind0).array());
            } else {
                T0.col(ind0) = T0.col(ind0).array() * ((beta_c - beta_p) * df0.col(ind0)).array().exp().array();
                Td0.col(ind0) = Td0.col(ind0).array() * ((beta_c - beta_p) * df0.col(ind0)).array().exp().array();
                Tdd0.col(ind0) = Tdd0.col(ind0).array() * ((beta_c - beta_p) * df0.col(ind0)).array().exp().array();
            }
            if (modelform=="A"){//The risks can be updated instead of recalculated for the models/terms allowed
                if (tform[ind0] == "lin"){
                    R.col(0) = R.col(0).array() + (beta_c - beta_p) * df0.col(ind0).array();
                } else if (tform[ind0] == "plin"){
                    R.col(0) = R.col(0).array() + (beta_c - beta_p) * df0.col(ind0).array();
                } else{
                    R.col(0) = R.col(0).array() + (beta_c * df0.col(ind0)).array().exp() - (beta_p * df0.col(ind0)).array().exp();
                    Rd.col(ind0) = Rd.col(ind0).array() * ((beta_c - beta_p) * df0.col(ind0)).array().exp().array();
                    Rdd.col(ind0*totalnum+ind0) = Rdd.col(ind0*totalnum+ind0).array() * ((beta_c - beta_p) * df0.col(ind0)).array().exp().array();
                }
            } else if ((modelform=="PA")||(modelform=="PAE")){
                if (tform[ind0] == "lin"){
                    if (ind0==fir){
                        R.col(0) = R.col(0) * beta_c / beta_p;
                    } else {
                        R.col(0) = R.col(0).array() + (beta_c - beta_p) * df0.col(ind0).array() * T0.col(fir).array();
                        Te.col(0) = Te.col(0).array() + (beta_c - beta_p) * df0.col(ind0).array();
                    }
                } else if (tform[ind0] == "plin"){
                    if (ind0==fir){
                        R.col(0) = (R.col(0).array() - Te.array()) * beta_c / beta_p + Te.col(0).array();
                    } else {
                        R.col(0) = R.col(0).array() + (beta_c - beta_p) * df0.col(ind0).array() * T0.col(fir).array();
                        Te.col(0) = Te.col(0).array() + (beta_c - beta_p) * df0.col(ind0).array();
                    }
                } else{
                    if (ind0==fir){
                        R.col(0) = R.col(0).array() * ((beta_c - beta_p) * df0.col(ind0)).array().exp();
                    } else {
                        R.col(0) = R.col(0).array() + ((beta_c * df0.col(ind0)).array().exp() - (beta_p * df0.col(ind0)).array().exp()).array() * T0.col(fir).array();
                        Te.col(0) = Te.col(0).array() + ((beta_c * df0.col(ind0)).array().exp() - (beta_p * df0.col(ind0)).array().exp()).array();
                    }
                }
                #pragma omp parallel for schedule(dynamic) num_threads(nthreads) collapse(2)
                for (int bl=0;bl<batch_cols.size()-1;bl++){
                    for (int ijk=0;ijk<totalnum*(totalnum+1)/2;ijk++){
                        int ij = 0;
                        int jk = ijk;
                        while (jk>ij){
                            ij++;
                            jk-=ij;
                        }
                        if (ij==ind0){
                            if (ij==jk){
                                if (tform[ij]=="loglin"){
                                    Rd.block(batch_cols[bl],ij,batch_cols[bl+1]-batch_cols[bl],1) = Rd.block(batch_cols[bl],ij,batch_cols[bl+1]-batch_cols[bl],1).array() * ((beta_c - beta_p) * df0.block(batch_cols[bl],ind0,batch_cols[bl+1]-batch_cols[bl],1)).array().exp();
                                    Rdd.block(batch_cols[bl],ij*totalnum+jk,batch_cols[bl+1]-batch_cols[bl],1) = Rdd.block(batch_cols[bl],ij,batch_cols[bl+1]-batch_cols[bl],1).array() * ((beta_c - beta_p) * df0.block(batch_cols[bl],ind0,batch_cols[bl+1]-batch_cols[bl],1)).array().exp();
                                } else {
                                    ;
                                }
                            } else{
                                if ((ij!=fir)&&(jk!=fir)){
                                    ;
                                } else if (tform[ind0]=="loglin"){
                                    Rdd.block(batch_cols[bl],ij*totalnum+jk,batch_cols[bl+1]-batch_cols[bl],1) = Rdd.block(batch_cols[bl],ij*totalnum+jk,batch_cols[bl+1]-batch_cols[bl],1).array() * ((beta_c - beta_p) * df0.block(batch_cols[bl],ind0,batch_cols[bl+1]-batch_cols[bl],1)).array().exp();
                                    Rdd.block(batch_cols[bl],jk*totalnum+ij,batch_cols[bl+1]-batch_cols[bl],1) = Rdd.block(batch_cols[bl],jk*totalnum+ij,batch_cols[bl+1]-batch_cols[bl],1).array() * ((beta_c - beta_p) * df0.block(batch_cols[bl],ind0,batch_cols[bl+1]-batch_cols[bl],1)).array().exp();
                                }
                            }
                        } else if (jk==ind0){
                            if ((ij!=fir)&&(jk!=fir)){
                                ;
                            } else if (tform[jk]=="loglin"){
                                Rdd.block(batch_cols[bl],jk*totalnum+ij,batch_cols[bl+1]-batch_cols[bl],1) = Rdd.block(batch_cols[bl],jk*totalnum+ij,batch_cols[bl+1]-batch_cols[bl],1).array() * ((beta_c - beta_p) * df0.block(batch_cols[bl],ind0,batch_cols[bl+1]-batch_cols[bl],1)).array().exp();
                                Rdd.block(batch_cols[bl],ij*totalnum+jk,batch_cols[bl+1]-batch_cols[bl],1) = Rdd.block(batch_cols[bl],ij*totalnum+jk,batch_cols[bl+1]-batch_cols[bl],1).array() * ((beta_c - beta_p) * df0.block(batch_cols[bl],ind0,batch_cols[bl+1]-batch_cols[bl],1)).array().exp();
                            }
                        }
                    }
                }
            }else if (modelform=="M"){
                Te = Te.array() * 0 + 1; //verifies the initial term product is 1
                //
                Te = T0.array().rowwise().prod();
                // computes intial risk and derivatives
                R << Te.array();
                Rd = T0.array().pow(-1).array() * Te.colwise().replicate(totalnum).array();
                Rd = Rd.array() * Td0.array();
                #pragma omp parallel for schedule(dynamic) num_threads(nthreads) collapse(2)
                for (int bl=0;bl<batch_cols.size()-1;bl++){
                    for (int ijk=0;ijk<totalnum*(totalnum+1)/2;ijk++){
                        int ij = 0;
                        int jk = ijk;
                        while (jk>ij){
                            ij++;
                            jk-=ij;
                        }
                        if (ij==jk){
                            Rdd.block(batch_cols[bl],ij*totalnum+jk,batch_cols[bl+1]-batch_cols[bl],1) = Tdd0.block(batch_cols[bl],jk,batch_cols[bl+1]-batch_cols[bl],1).array() * Td0.block(batch_cols[bl],jk,batch_cols[bl+1]-batch_cols[bl],1).array().pow(-1).array() * Rd.block(batch_cols[bl],ij,batch_cols[bl+1]-batch_cols[bl],1).array();
                        } else {
                            Rdd.block(batch_cols[bl],ij*totalnum+jk,batch_cols[bl+1]-batch_cols[bl],1) = Td0.block(batch_cols[bl],jk,batch_cols[bl+1]-batch_cols[bl],1).array() * T0.block(batch_cols[bl],jk,batch_cols[bl+1]-batch_cols[bl],1).array().pow(-1).array() * Rd.block(batch_cols[bl],ij,batch_cols[bl+1]-batch_cols[bl],1).array();
                            Rdd.block(batch_cols[bl],jk*totalnum+ij,batch_cols[bl+1]-batch_cols[bl],1) = Td0.block(batch_cols[bl],ij,batch_cols[bl+1]-batch_cols[bl],1).array() * T0.block(batch_cols[bl],ij,batch_cols[bl+1]-batch_cols[bl],1).array().pow(-1).array() * Rd.block(batch_cols[bl],jk,batch_cols[bl+1]-batch_cols[bl],1).array();
                        }
                    }
                }
            }else if (modelform=="GM"){
                throw invalid_argument( "GM isn't implemented" );
            } else {
                throw invalid_argument( "Model isn't implemented" );
            }
            R = (R.array().isFinite()).select(R,0);
            Rd = (Rd.array().isFinite()).select(Rd,0);
            Rdd = (Rdd.array().isFinite()).select(Rdd,0);
            fill(Ll.begin(), Ll.end(), 0.0);
            fill(Lld.begin(), Lld.end(), 0.0);
            fill(Lldd.begin(), Lldd.end(), 0.0);
            end_point = system_clock::now();
            end = time_point_cast<microseconds>(end_point).time_since_epoch().count();
            cout<<"df100 "<<(end-start)<<" "<<halves<<" "<<iteration<<" "<<ind0<<",Revert"<<endl;
            gibtime = system_clock::to_time_t(system_clock::now());
            cout << ctime(&gibtime) << endl;
            
            temp = (PyrC.col(1).array() * (PyrC.col(0).array() * R.col(0).array()).array().log()).array() - (PyrC.col(0).array() * R.col(0).array());
            fill(Ll.begin(), Ll.end(), (temp.array().isFinite()).select(temp,0).sum());

            temp = (PyrC.col(1).replicate(1,totalnum).array() * Rd.array() * R.replicate(1,totalnum).array().pow(-1).array()).array() - (PyrC.col(0).replicate(1,totalnum).array() * Rd.array());
            tempd << (temp.array().isFinite()).select(temp,0).colwise().sum();
            VectorXd::Map(&Lld[0], totalnum) = tempd;

            #pragma omp parallel for schedule(dynamic) num_threads(nthreads)
            for (int ijk=0;ijk<totalnum*(totalnum+1)/2;ijk++){//totalnum*(totalnum+1)/2
                int ij = 0;
                int jk = ijk;
                while (jk>ij){
                    ij++;
                    jk-=ij;
                }
                VectorXd temp(Rdd.rows(),Rdd.cols());
                temp = (PyrC.col(1).array() * (R.col(0).pow(-1).array() * Rdd.col(ij*totalnum+jk).array() - (Rd.col(ij).array() * Rd.col(jk).array() * R.col(0).pow(-1).array()).pow(2) )).array() - (PyrC.col(0).array() * Rdd.col(ij*totalnum+jk).array());
                Lldd[ij*totalnum+jk] = (temp.array().isFinite()).select(temp,0).sum();
                if (ij!=jk){
                    Lldd[jk*totalnum+ij] = (temp.array().isFinite()).select(temp,0).sum();
                }
            }
            beta_0[ind0] = beta_best;
        }
        for (int ij=0;ij<totalnum;ij++){
            if (abs(Lld[ij]) > Lld_worst){
                Lld_worst = abs(Lld[ij]);
            }
        }
        beta_0[ind0] = beta_best;
        if (iteration > totalnum * 3){
            if (iteration % (2 * totalnum)){
                if (Lld_worst < deriv_epsilon){
                    iteration = maxiter;
                }
            }
        }
        end_point = system_clock::now();
        end = time_point_cast<microseconds>(end_point).time_since_epoch().count();
        cout<<"df100 "<<(end-start)<<" "<<halves<<" "<<iteration<<" "<<ind0<<",Recalc"<<endl;
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
            cout << beta_0[ij] << " ";
        }
        cout << " " << endl;
    }
//    NumericVector 
    NumericVector Lldd_vec = wrap(Lldd);
    Lldd_vec.attr("dim") = Dimension(totalnum, totalnum);
    //
    List res_list = List::create(_["LogLik"]=wrap(Ll),_["First_Der"]=wrap(Lld),_["Second_Der"]=Lldd_vec,_["beta_0"]=wrap(beta_0),_["AIC"]=2*totalnum-2*Ll[fir]);
    //
    return res_list;
}

//MatrixXd Residual_AMFIT( MatrixXd PyrC, VectorXd beta_linT,VectorXd beta_loglinT,VectorXd beta_plinT,MatrixXd df_lin,MatrixXd df_loglin,MatrixXd df_plin,string modelform, NumericVector include_bool){
//    using namespace std::chrono;
//    cout << "START_NEW" << endl;
//    time_point<system_clock> start_point, end_point;
//    start_point = system_clock::now();
//    auto start = time_point_cast<microseconds>(start_point).time_since_epoch().count();
//    end_point = system_clock::now();
//    auto end = time_point_cast<microseconds>(end_point).time_since_epoch().count(); //The time duration is tracked
//    //
//    auto gibtime = system_clock::to_time_t(system_clock::now());
//    cout << ctime(&gibtime) << endl;
//    //
//    int totalnum = 0;
//    //
//    cout.precision(10); //forces higher precision numbers printed to terminal
//    int nthreads = Eigen::nbThreads(); //stores how many threads are allocated
//    //
//    VectorXd beta_lin;
//    VectorXd beta_loglin; //The vectors of parameter used
//    VectorXd beta_plin;
//    //
//    if (include_bool[0]==1){
//        beta_lin = beta_linT.tail(beta_linT.size()-1);
//    }
//    if (include_bool[1]==1){
//        beta_loglin = beta_loglinT.tail(beta_loglinT.size()-1); //creates the used vectors
//    }
//    if (include_bool[2]==1){
//        beta_plin = beta_plinT.tail(beta_plinT.size()-1);
//    }
//    //
//    if (include_bool[0]==1){
//        totalnum = totalnum + beta_lin.size();
//    }
//    if (include_bool[1]==1){
//        totalnum = totalnum + beta_loglin.size(); //determines how many parameters are needed
//    }
//    if (include_bool[2]==1){
//        totalnum = totalnum + beta_plin.size();
//    }
//    //
//    VectorXd res(totalnum); //preallocates a vector of final parameters
//    //
//    double Lld_worst = 0.0;
//    vector <string> tform(totalnum);
//    double totem = PyrC.rows(); //precalculates how many rows
//    //
//    vector <int> batch_cols;
//    int i_temp=0;
//    int j_temp=0;
//    batch_cols.push_back(0);
//    while (i_temp<totem-batch_size){
//        if (j_temp==batch_size){
//            batch_cols.push_back(i_temp);
//            j_temp=0;
//        }
//        i_temp++;
//        j_temp++;
//    }
//    if (totem-i_temp>batch_size/2){
//        batch_cols.push_back(i_temp);
//    }
//    batch_cols.push_back(totem-1);
//    //
//    //
//    end_point = system_clock::now();
//    end = time_point_cast<microseconds>(end_point).time_since_epoch().count();
//    cout<<"df99,"<<(end-start)<<",Batches"<<endl;
//    gibtime = system_clock::to_time_t(system_clock::now());
//    cout << ctime(&gibtime) << endl;
//    // ---------------------------------------------
//    // To Start, needs to seperate the derivative term
//    // ---------------------------------------------
//    VectorXd beta_0(totalnum);
//    MatrixXd df0 = MatrixXd::Zero(df_lin.rows(), totalnum); // stores memory for the derivative term parameters and columns
//    MatrixXd T0 = MatrixXd::Zero(df0.rows(), totalnum); //preallocates matrix for Derivative column terms
//    MatrixXd Td0 = MatrixXd::Zero(df0.rows(), totalnum); //preallocates matrix for Derivative column terms
//    MatrixXd Tdd0 = MatrixXd::Zero(df0.rows(), totalnum); //preallocates matrix for Derivative column terms
//    #pragma omp parallel for num_threads(nthreads)
//    for (int ij=0;ij<totalnum;ij++){
//        int ind0 = ij;
//        if (include_bool[0]==1){
//            if (ind0 < beta_lin.size()){
//                // one exists and is one
//                beta_0[ij] = beta_lin[ind0];
//                df0.col(ij) = df_lin.col(ind0);
//                tform[ij] = "lin";
//                //
//            } else {
//                //one exists and its not one
//                ind0 = ind0 - beta_lin.size();
//                if (include_bool[1]==1){
//                    if (ind0 < beta_loglin.size()){
//                        //one and two exists and is two
//                        beta_0[ij] = beta_loglin[ind0];
//                        df0.col(ij) = df_loglin.col(ind0);
//                        tform[ij] = "loglin";
//                        //
//                    } else{
//                        //one exists, two does, must be three
//                        if (include_bool[2]!=1){
//                            throw invalid_argument( "Are all three used?" );
//                        }
//                        ind0 = ind0 - beta_loglin.size();
//                        beta_0[ij] = beta_plin[ind0];
//                        df0.col(ij) = df_plin.col(ind0);
//                        tform[ij] = "plin";
//                        //
//                    }
//                } else{
//                    //one exists, and two doesn't exist, must be three
//                    if (include_bool[2]!=1){
//                        throw invalid_argument( "Are all first and third used?" );
//                    }
//                    beta_0[ij] = beta_plin[ind0];
//                    df0.col(ij) = df_plin.col(ind0);
//                    tform[ij] = "plin";
//                    //
//                }
//            }
//        }else{
//            //one doesn't exist
//            if (include_bool[1]==1){
//                if (ind0 < beta_loglin.size()){
//                    //one doesn't exist and two exists and is two
//                    beta_0[ij] = beta_loglin[ind0];
//    //                cout << ind0 << ", " << beta_0 << ", " << beta_loglin.transpose() << endl;
//                    df0.col(ij) = df_loglin.col(ind0);
//                    tform[ij] = "loglin";
//                    //
//                } else{
//                    //one doesn't exist, two does, must be three
//                    if (include_bool[2]!=1){
//                        throw invalid_argument( "Are all three used?" );
//                    }
//                    ind0 = ind0 - beta_loglin.size();
//                    beta_0[ij] = beta_plin[ind0];
//                    df0.col(ij) = df_plin.col(ind0);
//                    tform[ij] = "plin";
//                    //
//                }
//            } else{
//                //one doesn't exist, and two doesn't exist, must be three
//                if (include_bool[2]!=1){
//                    throw invalid_argument( "Are all first and third used?" );
//                }
//                beta_0[ij] = beta_plin[ind0];
//                df0.col(ij) = df_plin.col(ind0);
//                tform[ij] = "plin";
//                //
//            }
//        }
//        cout << df0.col(ij).array().transpose() << endl;
//        T0.col(ij) = (df0.col(ij).array() * beta_0[ij]).matrix();
//        if (tform[ij]=="lin") {
//            Td0.col(ij) = df0.col(ij);
//        } else if (tform[ij]=="loglin") {
//            T0.col(ij) = T0.col(ij).array().exp();
//            Td0.col(ij) = df0.col(ij).array() * T0.col(ij).array();
//            Tdd0.col(ij) = df0.col(ij).array() * Td0.col(ij).array();
//        } else {
//            T0.col(ij) = 1 + T0.col(ij).array();
//            Td0.col(ij) = df0.col(ij);
//        }
//    }
//    //
//    end_point = system_clock::now();
//    end = time_point_cast<microseconds>(end_point).time_since_epoch().count();
//    cout<<"df99,"<<(end-start)<<",Terms"<<endl;
//    gibtime = system_clock::to_time_t(system_clock::now());
//    cout << ctime(&gibtime) << endl;
//    //
//    MatrixXd Te = MatrixXd::Zero(df0.rows(), 1); //preallocates matrix for non-Derivative column terms
//    MatrixXd R = MatrixXd::Zero(df0.rows(), 1); //preallocates matrix for Risks and derivatives
//    MatrixXd Rd = MatrixXd::Zero(df0.rows(), totalnum); //preallocates matrix for Risks and derivatives
//    MatrixXd Rdd = MatrixXd::Zero(df0.rows(), pow(totalnum,2));
//    if ((modelform=="A")||(modelform=="PA")||(modelform=="PAE")){ //same process used for all of the additive type models
//        Te = T0.array().rowwise().sum();
//        // computes intial risk and derivatives
//        if (modelform=="A"){
//            R << Te.array();
//            Rd << Td0.array();
//            #pragma omp parallel for num_threads(nthreads)
//            for (int ij=0;ij<totalnum;ij++){
//                Rdd.col(ij*totalnum+ij) = Tdd0.col(ij);
//            }
//        } else if ((modelform=="PAE")||(modelform=="PA")){
//            Te = Te.array() - T0.col(fir).array();
//            if (modelform=="PAE"){
//                Te = Te.array() + 1;
//            }
//            R << T0.col(fir).array() * Te.array();
//            Rd << Td0.array() * T0.col(fir).array();//, Td0.col(0).array() * Te.array(), Td0.col(1).array() * Te.array();
//            Rd.col(fir) = Td0.col(fir).array() * Te.array();
//            #pragma omp parallel for schedule(dynamic) num_threads(nthreads) collapse(2)
//            for (int bl=0;bl<batch_cols.size()-1;bl++){
//                for (int ijk=0;ijk<totalnum*(totalnum+1)/2;ijk++){
//                    int ij = 0;
//                    int jk = ijk;
//                    while (jk>ij){
//                        ij++;
//                        jk-=ij;
//                    }
//                    if (ij==jk){
//                        if (ij==fir){
//                            Rdd.block(batch_cols[bl],ij*totalnum+jk,batch_cols[bl+1]-batch_cols[bl],1) = Tdd0.block(batch_cols[bl],ij,batch_cols[bl+1]-batch_cols[bl],1).array() * Te.block(batch_cols[bl],1,batch_cols[bl+1]-batch_cols[bl],1).array();
//                        } else {
//                            Rdd.block(batch_cols[bl],ij*totalnum+jk,batch_cols[bl+1]-batch_cols[bl],1) = Tdd0.block(batch_cols[bl],ij,batch_cols[bl+1]-batch_cols[bl],1).array() * T0.block(batch_cols[bl],fir,batch_cols[bl+1]-batch_cols[bl],1).array();
//                        }
//                    } else {
//                        if ((ij==fir)||(jk==fir)){
//                            Rdd.block(batch_cols[bl],ij*totalnum+jk,batch_cols[bl+1]-batch_cols[bl],1) = Td0.block(batch_cols[bl],ij,batch_cols[bl+1]-batch_cols[bl],1).array() * Td0.block(batch_cols[bl],jk,batch_cols[bl+1]-batch_cols[bl],1).array();
//                            Rdd.block(batch_cols[bl],ij*totalnum+jk,batch_cols[bl+1]-batch_cols[bl],1) = Td0.block(batch_cols[bl],ij,batch_cols[bl+1]-batch_cols[bl],1).array() * Td0.block(batch_cols[bl],jk,batch_cols[bl+1]-batch_cols[bl],1).array();
//                        }
//                    }
//                }
//            }
//        }
//    }else if (modelform=="M"){
//        Te = Te.array() * 0 + 1; //verifies the initial term product is 1
//        //
//        Te = T0.array().rowwise().prod();
//        // computes intial risk and derivatives
//        R << Te.array();
//        Rd = T0.array().pow(-1).array() * Te.colwise().replicate(totalnum).array();
//        Rd = Rd.array() * Td0.array();
//        #pragma omp parallel for schedule(dynamic) num_threads(nthreads) collapse(2)
//        for (int bl=0;bl<batch_cols.size()-1;bl++){
//            for (int ijk=0;ijk<totalnum*(totalnum+1)/2;ijk++){
//                int ij = 0;
//                int jk = ijk;
//                while (jk>ij){
//                    ij++;
//                    jk-=ij;
//                }
//                if (ij==jk){
//                    Rdd.block(batch_cols[bl],ij*totalnum+jk,batch_cols[bl+1]-batch_cols[bl],1) = Tdd0.block(batch_cols[bl],jk,batch_cols[bl+1]-batch_cols[bl],1).array() * Td0.block(batch_cols[bl],jk,batch_cols[bl+1]-batch_cols[bl],1).array().pow(-1).array() * Rd.block(batch_cols[bl],ij,batch_cols[bl+1]-batch_cols[bl],1).array();
//                } else {
//                    Rdd.block(batch_cols[bl],ij*totalnum+jk,batch_cols[bl+1]-batch_cols[bl],1) = Td0.block(batch_cols[bl],jk,batch_cols[bl+1]-batch_cols[bl],1).array() * T0.block(batch_cols[bl],jk,batch_cols[bl+1]-batch_cols[bl],1).array().pow(-1).array() * Rd.block(batch_cols[bl],ij,batch_cols[bl+1]-batch_cols[bl],1).array();
//                    Rdd.block(batch_cols[bl],jk*totalnum+ij,batch_cols[bl+1]-batch_cols[bl],1) = Td0.block(batch_cols[bl],ij,batch_cols[bl+1]-batch_cols[bl],1).array() * T0.block(batch_cols[bl],ij,batch_cols[bl+1]-batch_cols[bl],1).array().pow(-1).array() * Rd.block(batch_cols[bl],jk,batch_cols[bl+1]-batch_cols[bl],1).array();
//                }
//            }
//        }
//    } else if (modelform=="GM"){
//        //currently isn't implemented, it can be calculated but not optimized the same way
//        throw invalid_argument( "GM isn't implemented" );
//    } else {
//        throw invalid_argument( "Model isn't implemented" );
//    }
//    R = (R.array().isFinite()).select(R,0);
//    Rd = (Rd.array().isFinite()).select(Rd,0);
//    Rdd = (Rdd.array().isFinite()).select(Rdd,0);
//    //
//    ArrayXd Residuals(R.rows(),2);
//    //
//    ArrayXd guess = PyrC.col(0).array() * R.col(0).array();
//    ArrayXd dif = PyrC.col(1).array() - guess.array();
//    // Finds the deviance and pearson residuals respectively
//    Residuals.col(0) = pow(2,.5) * (dif.abs() * dif.pow(-1)) * (PyrC.col(1).array() * (PyrC.col(1).array() * guess.array().pow(-1)).array().log() - dif.array()).array().pow(.5);
//    Residuals.col(1) = dif.array() * guess.array().pow(-.5);
//    return Residuals;
//}
