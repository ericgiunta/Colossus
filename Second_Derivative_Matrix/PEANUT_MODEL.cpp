#include <RcppEigen.h>
#include <omp.h>
#include "PEANUT_MODEL.h"
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

using Eigen::Map;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using Rcpp::as;


template <typename T> int sign(T val) {
    return (T(0) < val) - (val < T(0));
}

// [[Rcpp::export]]
List peanut_transition(NumericVector a_lin,NumericVector a_loglin,NumericVector a_plin, NumericVector a_dose, NumericMatrix x_lin, NumericMatrix x_loglin, NumericMatrix x_plin, NumericMatrix x_dose,int fir,string modelform, string doseform, StringVector dose_terms,int ntime, NumericVector include_bool, double lr, int maxiter, int halfmax, double epsilon, double dbeta_cap, double deriv_epsilon, int batch_size){
    //----------------------------------------------------------------------------------------------------------------//
    //
//    cout << ntime << endl;
    Map<VectorXd> beta_lin(as<Map<VectorXd> >(a_lin));
    Map<VectorXd> beta_loglin(as<Map<VectorXd> >(a_loglin));
    Map<VectorXd> beta_plin(as<Map<VectorXd> >(a_plin));
    Map<VectorXd> beta_dose(as<Map<VectorXd> >(a_dose));
    const Map<MatrixXd> df_lin(as<Map<MatrixXd> >(x_lin));
    const Map<MatrixXd> df_loglin(as<Map<MatrixXd> >(x_loglin));
    const Map<MatrixXd> df_plin(as<Map<MatrixXd> >(x_plin));
    const Map<MatrixXd> df_dose(as<Map<MatrixXd> >(x_dose));
    // Converts from Rcpp types to efficient Eigen types
    //----------------------------------------------------------------------------------------------------------------//
    List res = LogLik_PEANUT(beta_lin,beta_loglin,beta_plin, beta_dose,df_lin,df_loglin,df_plin, df_dose,fir,modelform, doseform, dose_terms,ntime,include_bool, lr, maxiter, halfmax, epsilon, dbeta_cap, deriv_epsilon, batch_size);
    //----------------------------------------------------------------------------------------------------------------//
    return res;
}

List LogLik_PEANUT( VectorXd beta_linT,VectorXd beta_loglinT,VectorXd beta_plinT,VectorXd beta_dose,MatrixXd df_lin,MatrixXd df_loglin,MatrixXd df_plin, MatrixXd df_dose,int fir,string modelform, string doseform, StringVector dose_terms,int ntime,NumericVector include_bool, double lr, int maxiter, int halfmax, double epsilon, double dbeta_cap, double deriv_epsilon, int batch_size){
    //    cout << "start Risk" << endl;
    //
    srand (time(NULL));
    //
//    cout << ntime << endl;
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
    int totalnum = beta_dose.size();
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
    //for (int ij=0; ij < batch_cols.size(); ij++){
    //    cout << batch_cols[ij] << endl;
    //}
    end_point = system_clock::now();
    end = time_point_cast<microseconds>(end_point).time_since_epoch().count();
    cout<<"df99,"<<(end-start)<<",Batches"<<endl;
    gibtime = system_clock::to_time_t(system_clock::now());
    cout << ctime(&gibtime) << endl;
    // ---------------------------------------------
    // To Start, needs to seperate the derivative terms
    // ---------------------------------------------
    //
    //
    VectorXd beta_0(totalnum);
    MatrixXd df0 = MatrixXd::Zero(df_lin.rows(), totalnum); // stores memory for the derivative term parameters and columns
    MatrixXd T0 = MatrixXd::Zero(df0.rows(), totalnum); //preallocates matrix for Derivative column terms
    MatrixXd Td0 = MatrixXd::Zero(df0.rows(), totalnum); //preallocates matrix for Derivative column terms
    MatrixXd Tdd0 = MatrixXd::Zero(df0.rows(), totalnum); //preallocates matrix for Derivative column terms
    #pragma omp parallel for num_threads(nthreads)
    for (int ij=0;ij<totalnum;ij++){
        int ind0 = ij;
        if (ind0 < beta_dose.size()){
            beta_0[ij] = beta_dose[ind0];
            df0.col(ij) = df_dose.col(ind0);
            tform[ij] = dose_terms[ij];
        } else {
            ind0 = ind0 - beta_dose.size();
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
        }
        //        cout << df0.col(ij).array().transpose() << endl;
        T0.col(ij) = (df0.col(ij).array() * beta_0[ij]).matrix();
        if (tform[ij]=="lin") {
            Td0.col(ij) = df0.col(ij);
        } else if (tform[ij]=="loglin") {
            T0.col(ij) = T0.col(ij).array().exp();
            Td0.col(ij) = df0.col(ij).array() * T0.col(ij).array();
            Tdd0.col(ij) = df0.col(ij).array() * Td0.col(ij).array();
        } else if (tform[ij]=="plin") {
            T0.col(ij) = 1 + T0.col(ij).array();
            Td0.col(ij) = df0.col(ij);
        } else {
            cout << tform[ij] << " is invalid" << endl;
            throw invalid_argument( "Invalid term type" );
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
    MatrixXd De = MatrixXd::Zero(df0.rows(), 1); //preallocates matrix for non-Derivative column terms
    MatrixXd Dde = MatrixXd::Zero(df0.rows(), beta_dose.size()); //preallocates matrix for non-Derivative column terms
    MatrixXd Ddde = MatrixXd::Zero(df0.rows(), pow(beta_dose.size(),2)); //preallocates matrix for non-Derivative column terms
    MatrixXd R = MatrixXd::Zero(df0.rows(), 1); //preallocates matrix for Risks and derivatives
    MatrixXd Rd = MatrixXd::Zero(df0.rows(), totalnum); //preallocates matrix for Risks and derivatives
    MatrixXd Rdd = MatrixXd::Zero(df0.rows(), pow(totalnum,2));
    if (doseform=="A"){
        De = T0.array().block(0,0,T0.rows(),beta_dose.size()).rowwise().sum();
        Dde = Td0.array().block(0,0,T0.rows(),beta_dose.size());
        #pragma omp parallel for num_threads(nthreads)
        for (int ij=0;ij<beta_dose.size();ij++){
            Ddde.col(ij*beta_dose.size()+ij) = Tdd0.col(ij);
        }
    } else {
        throw invalid_argument( "That Dose model isn't implemented" );
    }
    if ((modelform=="A")||(modelform=="PA")||(modelform=="PAE")){ //same process used for all of the additive type models
        Te = T0.array().block(0,beta_dose.size(),T0.rows(),T0.cols()-beta_dose.size()).rowwise().sum() + De.array();
        // computes intial risk and derivatives
        if (modelform=="A"){
            R << Te.array();
            Rd << Dde.array(), Td0.array().block(0,beta_dose.size(),T0.rows(),T0.cols()-beta_dose.size());
            #pragma omp parallel for num_threads(nthreads)
            for (int ij=0;ij<totalnum;ij++){
                if (ij<beta_dose.size()){
                    Rdd.col(ij*totalnum+ij) = Ddde.col(ij);
                } else {
                    Rdd.col(ij*totalnum+ij) = Tdd0.col(ij);
                }
            }
        } else if ((modelform=="PAE")||(modelform=="PA")){
            if (fir!=0){
                Te = Te.array() - T0.col(fir).array();
            } else {
                Te = Te.array() - De.array();
            }
            if (modelform=="PAE"){
                Te = Te.array() + 1;
            }
            if (fir!=0){
                R << T0.col(fir).array() * Te.array();
                Rd << Td0.array() * T0.col(fir).array();//, Td0.col(0).array() * Te.array(), Td0.col(1).array() * Te.array();
                Rd.col(fir) = Td0.col(fir).array() * Te.array();
            } else {
                R << De.array() * Te.array();
                Rd << Td0.array() * De.array();
                #pragma omp parallel for schedule(dynamic) num_threads(nthreads)
                for (int ij=0;ij<beta_dose.size();ij++){
                    Rd.col(ij) = Dde.col(ij).array() * Te.array();
                }
            }
            #pragma omp parallel for schedule(dynamic) num_threads(nthreads) collapse(2)
            for (int bl=0;bl<batch_cols.size()-1;bl++){
                for (int ijk=0;ijk<totalnum;ijk++){
                    int ij = ijk;
                    int jk = ijk;
                    while (jk>ij){
                        ij++;
                        jk-=ij;
                    }
                    if (ij==jk){
                        if (fir!=0){
                            if (ij==fir){
                                Rdd.block(batch_cols[bl],ij*totalnum+jk,batch_cols[bl+1]-batch_cols[bl],1) = Tdd0.block(batch_cols[bl],ij,batch_cols[bl+1]-batch_cols[bl],1).array() * Te.block(batch_cols[bl],0,batch_cols[bl+1]-batch_cols[bl],1).array();
                            } else {
                                if (ij<beta_dose.size()){
                                    Rdd.block(batch_cols[bl],ij*totalnum+jk,batch_cols[bl+1]-batch_cols[bl],1) = Ddde.block(batch_cols[bl],ij,batch_cols[bl+1]-batch_cols[bl],1).array() * T0.block(batch_cols[bl],fir,batch_cols[bl+1]-batch_cols[bl],1).array();
                                } else {
                                    Rdd.block(batch_cols[bl],ij*totalnum+jk,batch_cols[bl+1]-batch_cols[bl],1) = Tdd0.block(batch_cols[bl],ij,batch_cols[bl+1]-batch_cols[bl],1).array() * T0.block(batch_cols[bl],fir,batch_cols[bl+1]-batch_cols[bl],1).array();
                                }
                            }
                        } else {
                            if (ij<beta_dose.size()){
                                Rdd.block(batch_cols[bl],ij*totalnum+jk,batch_cols[bl+1]-batch_cols[bl],1) = Ddde.block(batch_cols[bl],ij,batch_cols[bl+1]-batch_cols[bl],1).array() * Te.block(batch_cols[bl],0,batch_cols[bl+1]-batch_cols[bl],1).array();
                            } else {
                                Rdd.block(batch_cols[bl],ij*totalnum+jk,batch_cols[bl+1]-batch_cols[bl],1) = Tdd0.block(batch_cols[bl],ij,batch_cols[bl+1]-batch_cols[bl],1).array() * De.block(batch_cols[bl],0,batch_cols[bl+1]-batch_cols[bl],1).array();
                            }
                        }
                    } else {
                        if (fir!=0){
                            if ((ij==fir)||(jk==fir)){
                                if (ij<beta_dose.size()){
                                    Rdd.block(batch_cols[bl],ij*totalnum+jk,batch_cols[bl+1]-batch_cols[bl],1) = Dde.block(batch_cols[bl],ij,batch_cols[bl+1]-batch_cols[bl],1).array() * Td0.block(batch_cols[bl],jk,batch_cols[bl+1]-batch_cols[bl],1).array();
                                    Rdd.block(batch_cols[bl],jk*totalnum+ij,batch_cols[bl+1]-batch_cols[bl],1) = Dde.block(batch_cols[bl],ij,batch_cols[bl+1]-batch_cols[bl],1).array() * Td0.block(batch_cols[bl],jk,batch_cols[bl+1]-batch_cols[bl],1).array();
                                } else {
                                    Rdd.block(batch_cols[bl],ij*totalnum+jk,batch_cols[bl+1]-batch_cols[bl],1) = Td0.block(batch_cols[bl],ij,batch_cols[bl+1]-batch_cols[bl],1).array() * Td0.block(batch_cols[bl],jk,batch_cols[bl+1]-batch_cols[bl],1).array();
                                    Rdd.block(batch_cols[bl],jk*totalnum+ij,batch_cols[bl+1]-batch_cols[bl],1) = Td0.block(batch_cols[bl],ij,batch_cols[bl+1]-batch_cols[bl],1).array() * Td0.block(batch_cols[bl],jk,batch_cols[bl+1]-batch_cols[bl],1).array();
                                }
                            } else if ((ij<beta_dose.size())&&(jk<beta_dose.size())){
                                Rdd.block(batch_cols[bl],ij*totalnum+jk,batch_cols[bl+1]-batch_cols[bl],1) = Dde.block(batch_cols[bl],ij,batch_cols[bl+1]-batch_cols[bl],1).array() * Dde.block(batch_cols[bl],jk,batch_cols[bl+1]-batch_cols[bl],1).array();
                                Rdd.block(batch_cols[bl],jk*totalnum+ij,batch_cols[bl+1]-batch_cols[bl],1) = Dde.block(batch_cols[bl],ij,batch_cols[bl+1]-batch_cols[bl],1).array() * Dde.block(batch_cols[bl],jk,batch_cols[bl+1]-batch_cols[bl],1).array();
                            }
                        }
                    }
                }
            }
        }
    }else if (modelform=="M"){
        Te = Te.array() * 0 + 1; //verifies the initial term product is 1
        //
        Te = T0.array().block(0,beta_dose.size(),T0.rows(),T0.cols()-beta_dose.size()).rowwise().prod() * De.array();
        // computes intial risk and derivatives
        R << Te.array();
        #pragma omp parallel for schedule(dynamic) num_threads(nthreads)
        for (int ijk=0;ijk<totalnum;ijk++){
            if (ijk<beta_dose.size()){
                Rd.col(ijk) = De.col(ijk).array().pow(-1).array()* Te.array();
                Rd.col(ijk) = Rd.col(ijk).array() * Dde.col(ijk).array();
            }
        }
        Rd = T0.array().pow(-1).array() * Te.colwise().replicate(totalnum).array();
        Rd = Rd.array() * Td0.array();
        #pragma omp parallel for schedule(dynamic) num_threads(nthreads) collapse(2)
        for (int bl=0;bl<batch_cols.size()-1;bl++){
            for (int ijk=0;ijk<totalnum;ijk++){
                int ij = ijk;
                int jk = ijk;
                while (jk>ij){
                    ij++;
                    jk-=ij;
                }
                if (ij==jk){
                    if (ij<beta_dose.size()){
                        Rdd.block(batch_cols[bl],ij*totalnum+jk,batch_cols[bl+1]-batch_cols[bl],1) = Ddde.block(batch_cols[bl],jk,batch_cols[bl+1]-batch_cols[bl],1).array() * Dde.block(batch_cols[bl],jk,batch_cols[bl+1]-batch_cols[bl],1).array().pow(-1).array() * Rd.block(batch_cols[bl],ij,batch_cols[bl+1]-batch_cols[bl],1).array();
                    } else {
                        Rdd.block(batch_cols[bl],ij*totalnum+jk,batch_cols[bl+1]-batch_cols[bl],1) = Tdd0.block(batch_cols[bl],jk,batch_cols[bl+1]-batch_cols[bl],1).array() * Td0.block(batch_cols[bl],jk,batch_cols[bl+1]-batch_cols[bl],1).array().pow(-1).array() * Rd.block(batch_cols[bl],ij,batch_cols[bl+1]-batch_cols[bl],1).array();
                    }
                } else {
                    if ((ij<beta_dose.size())||(ij<beta_dose.size())){
                        Rdd.block(batch_cols[bl],ij*totalnum+jk,batch_cols[bl+1]-batch_cols[bl],1) = Dde.block(batch_cols[bl],ij,batch_cols[bl+1]-batch_cols[bl],1).array() * De.block(batch_cols[bl],ij,batch_cols[bl+1]-batch_cols[bl],1).array().pow(-1).array() * Rd.block(batch_cols[bl],ij,batch_cols[bl+1]-batch_cols[bl],1).array();
                        Rdd.block(batch_cols[bl],jk*totalnum+ij,batch_cols[bl+1]-batch_cols[bl],1) = Dde.block(batch_cols[bl],ij,batch_cols[bl+1]-batch_cols[bl],1).array() * De.block(batch_cols[bl],ij,batch_cols[bl+1]-batch_cols[bl],1).array().pow(-1).array() * Rd.block(batch_cols[bl],jk,batch_cols[bl+1]-batch_cols[bl],1).array();
                    } else{
                        Rdd.block(batch_cols[bl],ij*totalnum+jk,batch_cols[bl+1]-batch_cols[bl],1) = Td0.block(batch_cols[bl],jk,batch_cols[bl+1]-batch_cols[bl],1).array() * T0.block(batch_cols[bl],jk,batch_cols[bl+1]-batch_cols[bl],1).array().pow(-1).array() * Rd.block(batch_cols[bl],ij,batch_cols[bl+1]-batch_cols[bl],1).array();
                        Rdd.block(batch_cols[bl],jk*totalnum+ij,batch_cols[bl+1]-batch_cols[bl],1) = Td0.block(batch_cols[bl],ij,batch_cols[bl+1]-batch_cols[bl],1).array() * T0.block(batch_cols[bl],ij,batch_cols[bl+1]-batch_cols[bl],1).array().pow(-1).array() * Rd.block(batch_cols[bl],jk,batch_cols[bl+1]-batch_cols[bl],1).array();
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
    R = (R.array().isFinite()).select(R,0);
    Rd = (Rd.array().isFinite()).select(Rd,0);
    Rdd = (Rdd.array().isFinite()).select(Rdd,0);
    //
    // -------------------------------------------------------------------------------------------
    //
    //
    string line;
    ifstream infile("test.txt"); //The file of risk group rows
    vector<string>  RiskGroup(ntime); //vector of strings detailing the rows
    IntegerMatrix RiskFail(ntime,2); //vector giving the event rows
    int j_iter = 0;
    //
    end_point = system_clock::now();
    end = time_point_cast<microseconds>(end_point).time_since_epoch().count();
    cout<<"df100 "<<(end-start)<<" "<<0<<" "<<0<<" "<<-1<<",Prep_R"<<endl;
    gibtime = system_clock::to_time_t(system_clock::now());
    cout << ctime(&gibtime) << endl;
    // --------------------------
    // A file is created previously that stores what rows belong to which risk sections
    // --------------------------
    while (getline(infile, line)){
        istringstream iss(line);
        vector <int> row;
        string lineStream;
        while (getline(iss, lineStream, ',')){  
            row.push_back(stoi(lineStream));
        }
        vector<int> rows = vector<int>(row.end() - 2, row.end());
//        cout << rows.size() << ", " << rows[0] << ", " << rows[1] << ", " << row[row.size()-1] << endl;
        // --------------------------
        // needs the rows of every event with the same end time
        // --------------------------
        RiskFail(j_iter,0)=rows[0]-1;
        RiskFail(j_iter,1)=rows[1]-1;
        ostringstream group_str;
        copy(row.begin(), row.end()-2, ostream_iterator<int>(group_str, ","));
        RiskGroup[j_iter] = group_str.str();
//        cout << group_str.str() << endl;
//        cout << j_iter << endl;
        j_iter++;
    }
    //
    // --------------------------
    // now a vector exists with row locations
    // --------------------------
    //The log-likelihood is calculated in parallel over the risk groups
    vector<double> Ll(totalnum,0.0);
    vector<double> Lld(totalnum,0.0);
    vector<double> Lldd(pow(totalnum,2),0.0);
    #pragma omp declare reduction(vec_double_plus : std::vector<double> : \
        std::transform(omp_out.begin(), omp_out.end(), omp_in.begin(), omp_out.begin(), std::plus<double>())) \
        initializer(omp_priv = omp_orig)
    #pragma omp parallel for schedule(dynamic) num_threads(nthreads) reduction(vec_double_plus:Ll,Lld,Lldd) collapse(2)
    for (int ijk=0;ijk<totalnum;ijk++){//totalnum*(totalnum+1)/2
        for (int j=0;j<ntime;j++){
//            cout << "_____________________________" << endl;
            int ij = ijk;
            int jk = ijk;
            while (jk>ij){
                ij++;
                jk-=ij;
            }
            double Rs1 = 0;
            double Rs2 = 0;
            double Rs2t = 0;
            double Rs3 = 0;
//            cout << 1 << endl;
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
//            cout << 1 << endl;
            //Now has the grouping pairs
            int dj = RiskFail(j,1)-RiskFail(j,0)+1;
//            cout << 2 << endl;
//            cout << Groupstr << endl;
//            cout << InGroup.size() << endl;
//            int at_risk = 0;
            for (int i = 0; i < InGroup.size()-1; i=i+2){
//                cout << InGroup[i] << "," << InGroup[i+1] << "," << ij << "," << jk << endl;
                Rs1 += R.block(InGroup[i]-1,0,InGroup[i+1]-InGroup[i]+1,1).sum();
                Rs2 += Rd.block(InGroup[i]-1,ij,InGroup[i+1]-InGroup[i]+1,1).sum();
                Rs2t += Rd.block(InGroup[i]-1,jk,InGroup[i+1]-InGroup[i]+1,1).sum();
                Rs3 += Rdd.block(InGroup[i]-1,ij*totalnum+jk,InGroup[i+1]-InGroup[i]+1,1).sum();
//                at_risk += InGroup[i+1]-InGroup[i]+1;
            }
//            cout << 1 << endl;
            //
            MatrixXd Ld = MatrixXd::Zero(dj,4);
            Ld << R.block(RiskFail(j,0),0,dj,1), Rd.block(RiskFail(j,0),ij,dj,1), Rd.block(RiskFail(j,0),jk,dj,1) ,Rdd.block(RiskFail(j,0),ij*totalnum+jk,dj,1);//sum of risks in group
            //
            MatrixXd Ldm = MatrixXd::Zero(dj,4);
            for (int i = 0; i < dj; i++){ //adds in the efron approximation terms
                Ldm.row(i) = (-double(i) / double(dj)) * Ld.colwise().sum().array();
            }
//            cout << 1 << endl;
            Ldm.col(0) = Ldm.col(0).array() + Rs1;
            Ldm.col(1) = Ldm.col(1).array() + Rs2;
            Ldm.col(2) = Ldm.col(2).array() + Rs2t;
            Ldm.col(3) = Ldm.col(3).array() + Rs3;
            //
//            cout << 1 << endl;
            MatrixXd temp1 = MatrixXd::Zero(Ld.rows(),1);
            MatrixXd temp2 = MatrixXd::Zero(Ld.rows(),1);
            temp1 = Ld.col(0).array().log();
            double Ld1 =  (temp1.array().isFinite()).select(temp1,0).sum();
            temp1 = Ld.col(1).array() * (Ld.col(0).array().pow(-1).array());
            temp2 = Ld.col(2).array() * (Ld.col(0).array().pow(-1).array());
            double Ld2 = (temp1.array().isFinite()).select(temp1,0).sum();
            temp1 = Ld.col(3).array() * (Ld.col(0).array().pow(-1).array()) - temp1.array() * temp2.array();
            double Ld3 = (temp1.array().isFinite()).select(temp1,0).sum();
            //
//            cout << 1 << endl;
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
//            cout << 1 << endl;
//            cout <<ij << "," << jk << "," << j << "," << Ld1 - Rs1 <<","<< Ld2 - Rs2 <<","<< Ld3 - Rs3 << endl;
        }
    }
//    cout << 1 << endl;
    #pragma omp parallel for num_threads(nthreads)
    for (int ijk=0;ijk<totalnum;ijk++){//totalnum*(totalnum+1)/2
        int ij = ijk;
        int jk = ijk;
        while (jk>ij){
            ij++;
            jk-=ij;
        }
        Lldd[jk*totalnum+ij] = Lldd[ij*totalnum+jk];
    }
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
    VectorXd::Map(&beta_p[0], beta_0.size()) = beta_0;//wrap(beta_0);
    VectorXd::Map(&beta_c[0], beta_0.size()) = beta_0;//beta_c = wrap(beta_0);
    VectorXd::Map(&beta_a[0], beta_0.size()) = beta_0;//beta_a = wrap(beta_0);
    VectorXd::Map(&beta_best[0], beta_0.size()) = beta_0;//beta_best = wrap(beta_0);
    double Ll_best = 0.0;
    int halves = 0;
    int ind0 = fir;
    int i = ind0;
    int iteration=0;
    //
    while (iteration < maxiter){
        iteration++;
        VectorXd::Map(&beta_p[0], beta_0.size()) = beta_0;//wrap(beta_0);
        VectorXd::Map(&beta_c[0], beta_0.size()) = beta_0;//beta_c = wrap(beta_0);
        VectorXd::Map(&beta_a[0], beta_0.size()) = beta_0;//beta_a = wrap(beta_0);
        VectorXd::Map(&beta_best[0], beta_0.size()) = beta_0;//beta_best = wrap(beta_0);
        //
        #pragma omp parallel for schedule(dynamic) num_threads(nthreads)
        for (int ijk=0;ijk<totalnum;ijk++){
            if ((Lld[ijk] > 0)&&(Lldd[ijk*totalnum+ijk]<0)){
                dbeta[ijk] = -lr * Lld[ijk] / Lldd[ijk*totalnum+ijk];// decreases left, zeros right, ratio left
            } else if ((Lld[ijk]<0)&&(Lldd[ijk*totalnum+ijk]>0)){
                dbeta[ijk] = -lr * Lld[ijk] / Lldd[ijk*totalnum+ijk];// decreases right, zeros right, ratio left
            } else if ((Lld[ijk] > 0)&&(Lldd[ijk*totalnum+ijk]>0)){
                dbeta[ijk] = -lr * Lld[ijk] / Lldd[ijk*totalnum+ijk];// decreases left, zeros left, ratio right
            } else if ((Lld[ijk]<0)&&(Lldd[ijk*totalnum+ijk]<0)){
                dbeta[ijk] = -lr * Lld[ijk] / Lldd[ijk*totalnum+ijk];// decreases right, zeros left, ratio right
            } else {
                dbeta[ijk]=0.0;
            }
            //
            double dbeta_max = abs(Ll[ijk]/Lld[ijk] * dbeta_cap);
            if (abs(dbeta[ijk])>dbeta_max){
                dbeta[ijk] = dbeta_max * sign(dbeta[ijk]);
            }
            if (ijk!=fir){
                dbeta[ijk] = 0.0;
            }
//            dbeta[ijk] = 0.0;
        }
        //
        VectorXd::Map(&beta_p[0], beta_0.size()) = beta_0;//wrap(beta_0);
        VectorXd::Map(&beta_c[0], beta_0.size()) = beta_0;//beta_c = wrap(beta_0);
        VectorXd::Map(&beta_a[0], beta_0.size()) = beta_0;//beta_a = wrap(beta_0);
        VectorXd::Map(&beta_best[0], beta_0.size()) = beta_0;//beta_best = wrap(beta_0);
        Ll_best = Ll[ind0];
        // int halves = 0;
        i = ind0;
//        cout << beta_best << ", " << Ll[ind0] << endl;
        //
        halves=0;
        while ((Ll[ind0] <= Ll_best)&&(halves<halfmax)&&(abs(dbeta[ind0]) > epsilon)){
            halves++;
            #pragma omp parallel for schedule(dynamic) num_threads(nthreads)
            for (int ijk=0;ijk<totalnum;ijk++){
                beta_c[ijk] = beta_a[ijk] + dbeta[ijk];
                if (tform[ijk]=="lin"){
                    T0.col(ijk) = T0.col(ijk).array() * (beta_c[ijk] / beta_p[ijk]);
                } else if (tform[ijk]=="plin"){
                    T0.col(ijk) = T0.col(ijk).array() * (1 + beta_c[ijk] * df0.col(ijk).array()) / (1 + beta_p[ijk] * df0.col(ijk).array());
                } else {
                    T0.col(ijk) = T0.col(ijk).array() * ((beta_c[ijk] - beta_p[ijk]) * df0.col(ijk)).array().exp().array();
                    Td0.col(ijk) = Td0.col(ijk).array() * ((beta_c[ijk] - beta_p[ijk]) * df0.col(ijk)).array().exp().array();
                    Tdd0.col(ijk) = Tdd0.col(ijk).array() * ((beta_c[ijk] - beta_p[ijk]) * df0.col(ijk)).array().exp().array();
                }
            }
            if (doseform=="A"){
                De = T0.array().block(0,0,T0.rows(),beta_dose.size()).rowwise().sum();
                Dde = Td0.array().block(0,0,T0.rows(),beta_dose.size());
                #pragma omp parallel for num_threads(nthreads)
                for (int ij=0;ij<beta_dose.size();ij++){
                    Ddde.col(ij*totalnum+ij) = Tdd0.col(ij);
                }
            }
            if ((modelform=="A")||(modelform=="PA")||(modelform=="PAE")){ //same process used for all of the additive type models
                Te = T0.array().block(0,beta_dose.size(),T0.rows(),T0.cols()-beta_dose.size()).rowwise().sum() + De.array();
                // computes intial risk and derivatives
                if (modelform=="A"){
                    R << Te.array();
                    Rd << Dde.array(), Td0.array().block(0,beta_dose.size(),T0.rows(),T0.cols()-beta_dose.size());
                    #pragma omp parallel for num_threads(nthreads)
                    for (int ij=0;ij<totalnum;ij++){
                        if (ij<beta_dose.size()){
                            Rdd.col(ij*totalnum+ij) = Ddde.col(ij);
                        } else {
                            Rdd.col(ij*totalnum+ij) = Tdd0.col(ij);
                        }
                    }
                } else if ((modelform=="PAE")||(modelform=="PA")){
                    if (fir!=0){
                        Te = Te.array() - T0.col(fir).array();
                    } else {
                        Te = Te.array() - De.array();
                    }
                    if (modelform=="PAE"){
                        Te = Te.array() + 1;
                    }
                    if (fir!=0){
                        R << T0.col(fir).array() * Te.array();
                        Rd << Td0.array() * T0.col(fir).array();//, Td0.col(0).array() * Te.array(), Td0.col(1).array() * Te.array();
                        Rd.col(fir) = Td0.col(fir).array() * Te.array();
                    } else {
                        R << De.array() * Te.array();
                        Rd << Td0.array() * De.array();
                        #pragma omp parallel for schedule(dynamic) num_threads(nthreads)
                        for (int ij=0;ij<beta_dose.size();ij++){
                            Rd.col(ij) = Dde.col(ij).array() * Te.array();
                        }
                    }
                    #pragma omp parallel for schedule(dynamic) num_threads(nthreads) collapse(2)
                    for (int bl=0;bl<batch_cols.size()-1;bl++){

                        for (int ijk=0;ijk<totalnum;ijk++){
                            int ij = ijk;
                            int jk = ijk;
                            while (jk>ij){
                                ij++;
                                jk-=ij;
                            }
                            if (ij==jk){
                                if (fir!=0){
                                    if (ij==fir){
                                        Rdd.block(batch_cols[bl],ij*totalnum+jk,batch_cols[bl+1]-batch_cols[bl],1) = Tdd0.block(batch_cols[bl],ij,batch_cols[bl+1]-batch_cols[bl],1).array() * Te.block(batch_cols[bl],0,batch_cols[bl+1]-batch_cols[bl],1).array();
                                    } else {
                                        if (ij<beta_dose.size()){
                                            Rdd.block(batch_cols[bl],ij*totalnum+jk,batch_cols[bl+1]-batch_cols[bl],1) = Ddde.block(batch_cols[bl],ij,batch_cols[bl+1]-batch_cols[bl],1).array() * T0.block(batch_cols[bl],fir,batch_cols[bl+1]-batch_cols[bl],1).array();
                                        } else {
                                            Rdd.block(batch_cols[bl],ij*totalnum+jk,batch_cols[bl+1]-batch_cols[bl],1) = Tdd0.block(batch_cols[bl],ij,batch_cols[bl+1]-batch_cols[bl],1).array() * T0.block(batch_cols[bl],fir,batch_cols[bl+1]-batch_cols[bl],1).array();
                                        }
                                    }
                                } else {
                                    if (ij<beta_dose.size()){
                                        Rdd.block(batch_cols[bl],ij*totalnum+jk,batch_cols[bl+1]-batch_cols[bl],1) = Ddde.block(batch_cols[bl],ij,batch_cols[bl+1]-batch_cols[bl],1).array() * Te.block(batch_cols[bl],0,batch_cols[bl+1]-batch_cols[bl],1).array();
                                    } else {
                                        Rdd.block(batch_cols[bl],ij*totalnum+jk,batch_cols[bl+1]-batch_cols[bl],1) = Tdd0.block(batch_cols[bl],ij,batch_cols[bl+1]-batch_cols[bl],1).array() * De.block(batch_cols[bl],0,batch_cols[bl+1]-batch_cols[bl],1).array();
                                    }
                                }
                            } else {
                                if (fir!=0){
                                    if ((ij==fir)||(jk==fir)){
                                        if (ij<beta_dose.size()){
                                            Rdd.block(batch_cols[bl],ij*totalnum+jk,batch_cols[bl+1]-batch_cols[bl],1) = Dde.block(batch_cols[bl],ij,batch_cols[bl+1]-batch_cols[bl],1).array() * Td0.block(batch_cols[bl],jk,batch_cols[bl+1]-batch_cols[bl],1).array();
                                            Rdd.block(batch_cols[bl],jk*totalnum+ij,batch_cols[bl+1]-batch_cols[bl],1) = Dde.block(batch_cols[bl],ij,batch_cols[bl+1]-batch_cols[bl],1).array() * Td0.block(batch_cols[bl],jk,batch_cols[bl+1]-batch_cols[bl],1).array();
                                        } else {
                                            Rdd.block(batch_cols[bl],ij*totalnum+jk,batch_cols[bl+1]-batch_cols[bl],1) = Td0.block(batch_cols[bl],ij,batch_cols[bl+1]-batch_cols[bl],1).array() * Td0.block(batch_cols[bl],jk,batch_cols[bl+1]-batch_cols[bl],1).array();
                                            Rdd.block(batch_cols[bl],jk*totalnum+ij,batch_cols[bl+1]-batch_cols[bl],1) = Td0.block(batch_cols[bl],ij,batch_cols[bl+1]-batch_cols[bl],1).array() * Td0.block(batch_cols[bl],jk,batch_cols[bl+1]-batch_cols[bl],1).array();
                                        }
                                    } else if ((ij<beta_dose.size())&&(jk<beta_dose.size())){
                                        Rdd.block(batch_cols[bl],ij*totalnum+jk,batch_cols[bl+1]-batch_cols[bl],1) = Dde.block(batch_cols[bl],ij,batch_cols[bl+1]-batch_cols[bl],1).array() * Dde.block(batch_cols[bl],jk,batch_cols[bl+1]-batch_cols[bl],1).array();
                                        Rdd.block(batch_cols[bl],jk*totalnum+ij,batch_cols[bl+1]-batch_cols[bl],1) = Dde.block(batch_cols[bl],ij,batch_cols[bl+1]-batch_cols[bl],1).array() * Dde.block(batch_cols[bl],jk,batch_cols[bl+1]-batch_cols[bl],1).array();
                                    }
                                }
                            }
                        }
                    }
                }
            }else if (modelform=="M"){
                Te = Te.array() * 0 + 1; //verifies the initial term product is 1
                //
                Te = T0.array().block(0,beta_dose.size(),T0.rows(),T0.cols()-beta_dose.size()).rowwise().sum() * De.array();
                // computes intial risk and derivatives
                R << Te.array();
                #pragma omp parallel for schedule(dynamic) num_threads(nthreads)
                for (int ijk=0;ijk<totalnum;ijk++){
                    if (ijk<beta_dose.size()){
                        Rd.col(ijk) = De.col(ijk).array().pow(-1).array()* Te.array();
                        Rd.col(ijk) = Rd.col(ijk).array() * Dde.col(ijk).array();
                    }
                }
                Rd = T0.array().pow(-1).array() * Te.colwise().replicate(totalnum).array();
                Rd = Rd.array() * Td0.array();
                #pragma omp parallel for schedule(dynamic) num_threads(nthreads) collapse(2)
                for (int bl=0;bl<batch_cols.size()-1;bl++){
                    for (int ijk=0;ijk<totalnum;ijk++){
                        int ij = ijk;
                        int jk = ijk;
                        while (jk>ij){
                            ij++;
                            jk-=ij;
                        }
                        if (ij==jk){
                            if (ij<beta_dose.size()){
                                Rdd.block(batch_cols[bl],ij*totalnum+jk,batch_cols[bl+1]-batch_cols[bl],1) = Ddde.block(batch_cols[bl],jk,batch_cols[bl+1]-batch_cols[bl],1).array() * Dde.block(batch_cols[bl],jk,batch_cols[bl+1]-batch_cols[bl],1).array().pow(-1).array() * Rd.block(batch_cols[bl],ij,batch_cols[bl+1]-batch_cols[bl],1).array();
                            } else {
                                Rdd.block(batch_cols[bl],ij*totalnum+jk,batch_cols[bl+1]-batch_cols[bl],1) = Tdd0.block(batch_cols[bl],jk,batch_cols[bl+1]-batch_cols[bl],1).array() * Td0.block(batch_cols[bl],jk,batch_cols[bl+1]-batch_cols[bl],1).array().pow(-1).array() * Rd.block(batch_cols[bl],ij,batch_cols[bl+1]-batch_cols[bl],1).array();
                            }
                        } else {
                            if ((ij<beta_dose.size())||(ij<beta_dose.size())){
                                Rdd.block(batch_cols[bl],ij*totalnum+jk,batch_cols[bl+1]-batch_cols[bl],1) = Dde.block(batch_cols[bl],ij,batch_cols[bl+1]-batch_cols[bl],1).array() * De.block(batch_cols[bl],ij,batch_cols[bl+1]-batch_cols[bl],1).array().pow(-1).array() * Rd.block(batch_cols[bl],ij,batch_cols[bl+1]-batch_cols[bl],1).array();
                                Rdd.block(batch_cols[bl],jk*totalnum+ij,batch_cols[bl+1]-batch_cols[bl],1) = Dde.block(batch_cols[bl],ij,batch_cols[bl+1]-batch_cols[bl],1).array() * De.block(batch_cols[bl],ij,batch_cols[bl+1]-batch_cols[bl],1).array().pow(-1).array() * Rd.block(batch_cols[bl],jk,batch_cols[bl+1]-batch_cols[bl],1).array();
                            } else{
                                Rdd.block(batch_cols[bl],ij*totalnum+jk,batch_cols[bl+1]-batch_cols[bl],1) = Td0.block(batch_cols[bl],jk,batch_cols[bl+1]-batch_cols[bl],1).array() * T0.block(batch_cols[bl],jk,batch_cols[bl+1]-batch_cols[bl],1).array().pow(-1).array() * Rd.block(batch_cols[bl],ij,batch_cols[bl+1]-batch_cols[bl],1).array();
                                Rdd.block(batch_cols[bl],jk*totalnum+ij,batch_cols[bl+1]-batch_cols[bl],1) = Td0.block(batch_cols[bl],ij,batch_cols[bl+1]-batch_cols[bl],1).array() * T0.block(batch_cols[bl],ij,batch_cols[bl+1]-batch_cols[bl],1).array().pow(-1).array() * Rd.block(batch_cols[bl],jk,batch_cols[bl+1]-batch_cols[bl],1).array();
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
            #pragma omp parallel for schedule(dynamic) num_threads(nthreads) reduction(vec_double_plus:Ll,Lld,Lldd) collapse(2)
            for (int ijk=0;ijk<totalnum;ijk++){
                for (int j=0;j<ntime;j++){
                    int ij = ijk;
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
//                    int at_risk = 0;
                    for (int i = 0; i < InGroup.size()-1; i=i+2){
                        Rs1 += R.block(InGroup[i]-1,0,InGroup[i+1]-InGroup[i]+1,1).sum();
                        Rs2 += Rd.block(InGroup[i]-1,ij,InGroup[i+1]-InGroup[i]+1,1).sum();
                        Rs2t += Rd.block(InGroup[i]-1,jk,InGroup[i+1]-InGroup[i]+1,1).sum();
                        Rs3 += Rdd.block(InGroup[i]-1,ij*totalnum+jk,InGroup[i+1]-InGroup[i]+1,1).sum();
//                        at_risk += InGroup[i+1]-InGroup[i]+1;
                    }
                    //
                    MatrixXd Ld = MatrixXd::Zero(dj,4);
                    Ld << R.block(RiskFail(j,0),0,dj,1), Rd.block(RiskFail(j,0),ij,dj,1), Rd.block(RiskFail(j,0),jk,dj,1) ,Rdd.block(RiskFail(j,0),ij*totalnum+jk,dj,1);//sum of risks in group
                    //
                    MatrixXd Ldm = MatrixXd::Zero(dj,4);
                    for (int i = 0; i < dj; i++){ //adds in the efron approximation terms
                        Ldm.row(i) = (-double(i) / double(dj)) * Ld.colwise().sum().array();
                    }
                    Ldm.col(0) = Ldm.col(0).array() + Rs1;
                    Ldm.col(1) = Ldm.col(1).array() + Rs2;
                    Ldm.col(2) = Ldm.col(2).array() + Rs2t;
                    Ldm.col(3) = Ldm.col(3).array() + Rs3;
                    //
                    MatrixXd temp1 = MatrixXd::Zero(Ld.rows(),1);
                    MatrixXd temp2 = MatrixXd::Zero(Ld.rows(),1);
                    temp1 = Ld.col(0).array().log();
                    double Ld1 =  (temp1.array().isFinite()).select(temp1,0).sum();
                    temp1 = Ld.col(1).array() * (Ld.col(0).array().pow(-1).array());
                    temp2 = Ld.col(2).array() * (Ld.col(0).array().pow(-1).array());
                    double Ld2 = (temp1.array().isFinite()).select(temp1,0).sum();
                    temp1 = Ld.col(3).array() * (Ld.col(0).array().pow(-1).array()) - temp1.array() * temp2.array();
                    double Ld3 = (temp1.array().isFinite()).select(temp1,0).sum();
                    //
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
        //            cout <<ij << "," << jk << "," << j << "," << Ld1 - Rs1 <<","<< Ld2 - Rs2 <<","<< Ld3 - Rs3 << endl;
                }
            }
            #pragma omp parallel for num_threads(nthreads)
            for (int ijk=0;ijk<totalnum;ijk++){//totalnum*(totalnum+1)/2
                int ij = ijk;
                int jk = ijk;
                while (jk>ij){
                    ij++;
                    jk-=ij;
                }
                Lldd[jk*totalnum+ij] = Lldd[ij*totalnum+jk];
            }
            if (Ll[ind0] <= Ll_best){
                #pragma omp parallel for num_threads(nthreads)
                for (int ijk=0;ijk<totalnum;ijk++){//totalnum*(totalnum+1)/2
                    dbeta[ijk] = dbeta[ijk] / 2.0;
                }
            } else{
                #pragma omp parallel for num_threads(nthreads)
                for (int ijk=0;ijk<totalnum;ijk++){//totalnum*(totalnum+1)/2
                    beta_best[ijk] = beta_c[ijk];
                }
            }
            end_point = system_clock::now();
            end = time_point_cast<microseconds>(end_point).time_since_epoch().count();
            #pragma omp parallel for num_threads(nthreads)
            for (int ijk=0;ijk<totalnum;ijk++){//totalnum*(totalnum+1)/2
                beta_p[ijk] = beta_c[ijk];
            }
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
            #pragma omp parallel for num_threads(nthreads)
            for (int ijk=0;ijk<totalnum;ijk++){//totalnum*(totalnum+1)/2
                beta_0[ijk] = beta_c[ijk];
            }
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
            #pragma omp parallel for num_threads(nthreads)
            for (int ijk=0;ijk<totalnum;ijk++){//totalnum*(totalnum+1)/2
                beta_0[ijk] = beta_best[ijk];
            }
//            cout << beta_best << ", " << Ll[ind0] << endl;
        }
        if (beta_best!=beta_p){
//            beta_c = beta_best;
//            cout << beta_c << ":" << beta_p << endl;
            #pragma omp parallel for schedule(dynamic) num_threads(nthreads)
            for (int ijk=0;ijk<totalnum;ijk++){
                beta_c[ijk] = beta_best[ijk];
                if (tform[ijk]=="lin"){
                    T0.col(ijk) = T0.col(ijk).array() * (beta_c[ijk] / beta_p[ijk]);
                } else if (tform[ijk]=="plin"){
                    T0.col(ijk) = T0.col(ijk).array() * (1 + beta_c[ijk] * df0.col(ijk).array()) / (1 + beta_p[ijk] * df0.col(ijk).array());
                } else {
                    T0.col(ijk) = T0.col(ijk).array() * ((beta_c[ijk] - beta_p[ijk]) * df0.col(ijk)).array().exp().array();
                    Td0.col(ijk) = Td0.col(ijk).array() * ((beta_c[ijk] - beta_p[ijk]) * df0.col(ijk)).array().exp().array();
                    Tdd0.col(ijk) = Tdd0.col(ijk).array() * ((beta_c[ijk] - beta_p[ijk]) * df0.col(ijk)).array().exp().array();
                }
            }
            if (doseform=="A"){
                De = T0.array().block(0,0,T0.rows(),beta_dose.size()).rowwise().sum();
                Dde = Td0.array().block(0,0,T0.rows(),beta_dose.size());
                #pragma omp parallel for num_threads(nthreads)
                for (int ij=0;ij<beta_dose.size();ij++){
                    Ddde.col(ij*totalnum+ij) = Tdd0.col(ij);
                }
            }
            if ((modelform=="A")||(modelform=="PA")||(modelform=="PAE")){ //same process used for all of the additive type models
                Te = T0.array().block(0,beta_dose.size(),T0.rows(),T0.cols()-beta_dose.size()).rowwise().sum() + De.array();
                // computes intial risk and derivatives
                if (modelform=="A"){
                    R << Te.array();
                    Rd << Dde.array(), Td0.array().block(0,beta_dose.size(),T0.rows(),T0.cols()-beta_dose.size());
                    #pragma omp parallel for num_threads(nthreads)
                    for (int ij=0;ij<totalnum;ij++){
                        if (ij<beta_dose.size()){
                            Rdd.col(ij*totalnum+ij) = Ddde.col(ij);
                        } else {
                            Rdd.col(ij*totalnum+ij) = Tdd0.col(ij);
                        }
                    }
                } else if ((modelform=="PAE")||(modelform=="PA")){
                    if (fir!=0){
                        Te = Te.array() - T0.col(fir).array();
                    } else {
                        Te = Te.array() - De.array();
                    }
                    if (modelform=="PAE"){
                        Te = Te.array() + 1;
                    }
                    if (fir!=0){
                        R << T0.col(fir).array() * Te.array();
                        Rd << Td0.array() * T0.col(fir).array();//, Td0.col(0).array() * Te.array(), Td0.col(1).array() * Te.array();
                        Rd.col(fir) = Td0.col(fir).array() * Te.array();
                    } else {
                        R << De.array() * Te.array();
                        Rd << Td0.array() * De.array();
                        #pragma omp parallel for schedule(dynamic) num_threads(nthreads)
                        for (int ij=0;ij<beta_dose.size();ij++){
                            Rd.col(ij) = Dde.col(ij).array() * Te.array();
                        }
                    }
                    #pragma omp parallel for schedule(dynamic) num_threads(nthreads) collapse(2)
                    for (int bl=0;bl<batch_cols.size()-1;bl++){
                        for (int ijk=0;ijk<totalnum;ijk++){
                            int ij = ijk;
                            int jk = ijk;
                            while (jk>ij){
                                ij++;
                                jk-=ij;
                            }
                            if (ij==jk){
                                if (fir!=0){
                                    if (ij==fir){
                                        Rdd.block(batch_cols[bl],ij*totalnum+jk,batch_cols[bl+1]-batch_cols[bl],1) = Tdd0.block(batch_cols[bl],ij,batch_cols[bl+1]-batch_cols[bl],1).array() * Te.block(batch_cols[bl],0,batch_cols[bl+1]-batch_cols[bl],1).array();
                                    } else {
                                        if (ij<beta_dose.size()){
                                            Rdd.block(batch_cols[bl],ij*totalnum+jk,batch_cols[bl+1]-batch_cols[bl],1) = Ddde.block(batch_cols[bl],ij,batch_cols[bl+1]-batch_cols[bl],1).array() * T0.block(batch_cols[bl],fir,batch_cols[bl+1]-batch_cols[bl],1).array();
                                        } else {
                                            Rdd.block(batch_cols[bl],ij*totalnum+jk,batch_cols[bl+1]-batch_cols[bl],1) = Tdd0.block(batch_cols[bl],ij,batch_cols[bl+1]-batch_cols[bl],1).array() * T0.block(batch_cols[bl],fir,batch_cols[bl+1]-batch_cols[bl],1).array();
                                        }
                                    }
                                } else {
                                    if (ij<beta_dose.size()){
                                        Rdd.block(batch_cols[bl],ij*totalnum+jk,batch_cols[bl+1]-batch_cols[bl],1) = Ddde.block(batch_cols[bl],ij,batch_cols[bl+1]-batch_cols[bl],1).array() * Te.block(batch_cols[bl],0,batch_cols[bl+1]-batch_cols[bl],1).array();
                                    } else {
                                        Rdd.block(batch_cols[bl],ij*totalnum+jk,batch_cols[bl+1]-batch_cols[bl],1) = Tdd0.block(batch_cols[bl],ij,batch_cols[bl+1]-batch_cols[bl],1).array() * De.block(batch_cols[bl],0,batch_cols[bl+1]-batch_cols[bl],1).array();
                                    }
                                }
                            } else {
                                if (fir!=0){
                                    if ((ij==fir)||(jk==fir)){
                                        if (ij<beta_dose.size()){
                                            Rdd.block(batch_cols[bl],ij*totalnum+jk,batch_cols[bl+1]-batch_cols[bl],1) = Dde.block(batch_cols[bl],ij,batch_cols[bl+1]-batch_cols[bl],1).array() * Td0.block(batch_cols[bl],jk,batch_cols[bl+1]-batch_cols[bl],1).array();
                                            Rdd.block(batch_cols[bl],jk*totalnum+ij,batch_cols[bl+1]-batch_cols[bl],1) = Dde.block(batch_cols[bl],ij,batch_cols[bl+1]-batch_cols[bl],1).array() * Td0.block(batch_cols[bl],jk,batch_cols[bl+1]-batch_cols[bl],1).array();
                                        } else {
                                            Rdd.block(batch_cols[bl],ij*totalnum+jk,batch_cols[bl+1]-batch_cols[bl],1) = Td0.block(batch_cols[bl],ij,batch_cols[bl+1]-batch_cols[bl],1).array() * Td0.block(batch_cols[bl],jk,batch_cols[bl+1]-batch_cols[bl],1).array();
                                            Rdd.block(batch_cols[bl],jk*totalnum+ij,batch_cols[bl+1]-batch_cols[bl],1) = Td0.block(batch_cols[bl],ij,batch_cols[bl+1]-batch_cols[bl],1).array() * Td0.block(batch_cols[bl],jk,batch_cols[bl+1]-batch_cols[bl],1).array();
                                        }
                                    } else if ((ij<beta_dose.size())&&(jk<beta_dose.size())){
                                        Rdd.block(batch_cols[bl],ij*totalnum+jk,batch_cols[bl+1]-batch_cols[bl],1) = Dde.block(batch_cols[bl],ij,batch_cols[bl+1]-batch_cols[bl],1).array() * Dde.block(batch_cols[bl],jk,batch_cols[bl+1]-batch_cols[bl],1).array();
                                        Rdd.block(batch_cols[bl],jk*totalnum+ij,batch_cols[bl+1]-batch_cols[bl],1) = Dde.block(batch_cols[bl],ij,batch_cols[bl+1]-batch_cols[bl],1).array() * Dde.block(batch_cols[bl],jk,batch_cols[bl+1]-batch_cols[bl],1).array();
                                    }
                                }
                            }
                        }
                    }
                }
            }else if (modelform=="M"){
                Te = Te.array() * 0 + 1; //verifies the initial term product is 1
                //
                Te = T0.array().block(0,beta_dose.size(),T0.rows(),T0.cols()-beta_dose.size()).rowwise().sum() * De.array();
                // computes intial risk and derivatives
                R << Te.array();
                #pragma omp parallel for schedule(dynamic) num_threads(nthreads)
                for (int ijk=0;ijk<totalnum;ijk++){
                    if (ijk<beta_dose.size()){
                        Rd.col(ijk) = De.col(ijk).array().pow(-1).array()* Te.array();
                        Rd.col(ijk) = Rd.col(ijk).array() * Dde.col(ijk).array();
                    }
                }
                Rd = T0.array().pow(-1).array() * Te.colwise().replicate(totalnum).array();
                Rd = Rd.array() * Td0.array();
                #pragma omp parallel for schedule(dynamic) num_threads(nthreads) collapse(2)
                for (int bl=0;bl<batch_cols.size()-1;bl++){
                    for (int ijk=0;ijk<totalnum;ijk++){
                        int ij = ijk;
                        int jk = ijk;
                        while (jk>ij){
                            ij++;
                            jk-=ij;
                        }
                        if (ij==jk){
                            if (ij<beta_dose.size()){
                                Rdd.block(batch_cols[bl],ij*totalnum+jk,batch_cols[bl+1]-batch_cols[bl],1) = Ddde.block(batch_cols[bl],jk,batch_cols[bl+1]-batch_cols[bl],1).array() * Dde.block(batch_cols[bl],jk,batch_cols[bl+1]-batch_cols[bl],1).array().pow(-1).array() * Rd.block(batch_cols[bl],ij,batch_cols[bl+1]-batch_cols[bl],1).array();
                            } else {
                                Rdd.block(batch_cols[bl],ij*totalnum+jk,batch_cols[bl+1]-batch_cols[bl],1) = Tdd0.block(batch_cols[bl],jk,batch_cols[bl+1]-batch_cols[bl],1).array() * Td0.block(batch_cols[bl],jk,batch_cols[bl+1]-batch_cols[bl],1).array().pow(-1).array() * Rd.block(batch_cols[bl],ij,batch_cols[bl+1]-batch_cols[bl],1).array();
                            }
                        } else {
                            if ((ij<beta_dose.size())||(ij<beta_dose.size())){
                                Rdd.block(batch_cols[bl],ij*totalnum+jk,batch_cols[bl+1]-batch_cols[bl],1) = Dde.block(batch_cols[bl],ij,batch_cols[bl+1]-batch_cols[bl],1).array() * De.block(batch_cols[bl],ij,batch_cols[bl+1]-batch_cols[bl],1).array().pow(-1).array() * Rd.block(batch_cols[bl],ij,batch_cols[bl+1]-batch_cols[bl],1).array();
                                Rdd.block(batch_cols[bl],jk*totalnum+ij,batch_cols[bl+1]-batch_cols[bl],1) = Dde.block(batch_cols[bl],ij,batch_cols[bl+1]-batch_cols[bl],1).array() * De.block(batch_cols[bl],ij,batch_cols[bl+1]-batch_cols[bl],1).array().pow(-1).array() * Rd.block(batch_cols[bl],jk,batch_cols[bl+1]-batch_cols[bl],1).array();
                            } else{
                                Rdd.block(batch_cols[bl],ij*totalnum+jk,batch_cols[bl+1]-batch_cols[bl],1) = Td0.block(batch_cols[bl],jk,batch_cols[bl+1]-batch_cols[bl],1).array() * T0.block(batch_cols[bl],jk,batch_cols[bl+1]-batch_cols[bl],1).array().pow(-1).array() * Rd.block(batch_cols[bl],ij,batch_cols[bl+1]-batch_cols[bl],1).array();
                                Rdd.block(batch_cols[bl],jk*totalnum+ij,batch_cols[bl+1]-batch_cols[bl],1) = Td0.block(batch_cols[bl],ij,batch_cols[bl+1]-batch_cols[bl],1).array() * T0.block(batch_cols[bl],ij,batch_cols[bl+1]-batch_cols[bl],1).array().pow(-1).array() * Rd.block(batch_cols[bl],jk,batch_cols[bl+1]-batch_cols[bl],1).array();
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
            #pragma omp parallel for schedule(dynamic) num_threads(nthreads) reduction(vec_double_plus:Ll,Lld,Lldd) collapse(2)
            for (int ijk=0;ijk<totalnum;ijk++){
                for (int j=0;j<ntime;j++){
                    int ij = ijk;
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
//                    int at_risk = 0;
                    for (int i = 0; i < InGroup.size()-1; i=i+2){
                        Rs1 += R.block(InGroup[i]-1,0,InGroup[i+1]-InGroup[i]+1,1).sum();
                        Rs2 += Rd.block(InGroup[i]-1,ij,InGroup[i+1]-InGroup[i]+1,1).sum();
                        Rs2t += Rd.block(InGroup[i]-1,jk,InGroup[i+1]-InGroup[i]+1,1).sum();
                        Rs3 += Rdd.block(InGroup[i]-1,ij*totalnum+jk,InGroup[i+1]-InGroup[i]+1,1).sum();
//                        at_risk += InGroup[i+1]-InGroup[i]+1;
                    }
                    //
                    MatrixXd Ld = MatrixXd::Zero(dj,4);
                    Ld << R.block(RiskFail(j,0),0,dj,1), Rd.block(RiskFail(j,0),ij,dj,1), Rd.block(RiskFail(j,0),jk,dj,1) ,Rdd.block(RiskFail(j,0),ij*totalnum+jk,dj,1);//sum of risks in group
                    //
                    MatrixXd Ldm = MatrixXd::Zero(dj,4);
                    for (int i = 0; i < dj; i++){ //adds in the efron approximation terms
                        Ldm.row(i) = (-double(i) / double(dj)) * Ld.colwise().sum().array();
                    }
                    Ldm.col(0) = Ldm.col(0).array() + Rs1;
                    Ldm.col(1) = Ldm.col(1).array() + Rs2;
                    Ldm.col(2) = Ldm.col(2).array() + Rs2t;
                    Ldm.col(3) = Ldm.col(3).array() + Rs3;
                    //
                    MatrixXd temp1 = MatrixXd::Zero(Ld.rows(),1);
                    MatrixXd temp2 = MatrixXd::Zero(Ld.rows(),1);
                    temp1 = Ld.col(0).array().log();
                    double Ld1 =  (temp1.array().isFinite()).select(temp1,0).sum();
                    temp1 = Ld.col(1).array() * (Ld.col(0).array().pow(-1).array());
                    temp2 = Ld.col(2).array() * (Ld.col(0).array().pow(-1).array());
                    double Ld2 = (temp1.array().isFinite()).select(temp1,0).sum();
                    temp1 = Ld.col(3).array() * (Ld.col(0).array().pow(-1).array()) - temp1.array() * temp2.array();
                    double Ld3 = (temp1.array().isFinite()).select(temp1,0).sum();
                    //
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
        //            cout <<ij << "," << jk << "," << j << "," << Ld1 - Rs1 <<","<< Ld2 - Rs2 <<","<< Ld3 - Rs3 << endl;
                }
            }
            #pragma omp parallel for num_threads(nthreads)
            for (int ijk=0;ijk<totalnum;ijk++){//totalnum*(totalnum+1)/2
                int ij = ijk;
                int jk = ijk;
                while (jk>ij){
                    ij++;
                    jk-=ij;
                }
                Lldd[jk*totalnum+ij] = Lldd[ij*totalnum+jk];
            }
            #pragma omp parallel for num_threads(nthreads)
            for (int ijk=0;ijk<totalnum;ijk++){//totalnum*(totalnum+1)/2
                beta_0[ijk] = beta_best[ijk];
            }
        }
        for (int ij=0;ij<totalnum;ij++){
            if (abs(Lld[ij]) > Lld_worst){
                Lld_worst = abs(Lld[ij]);
            }
        }
        #pragma omp parallel for num_threads(nthreads)
        for (int ijk=0;ijk<totalnum;ijk++){//totalnum*(totalnum+1)/2
            beta_0[ijk] = beta_best[ijk];
        }
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
    }
//    NumericVector 
    NumericVector Lldd_vec = wrap(Lldd);
    Lldd_vec.attr("dim") = Dimension(totalnum, totalnum);
    //
    List res_list = List::create(_["LogLik"]=wrap(Ll),_["First_Der"]=wrap(Lld),_["Second_Der"]=Lldd_vec,_["beta_0"]=wrap(beta_0),_["AIC"]=2*totalnum-2*Ll[fir]);
    //
    return res_list;
}
