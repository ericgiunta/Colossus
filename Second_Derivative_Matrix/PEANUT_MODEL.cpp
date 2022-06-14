#include <RcppEigen.h>
#include <omp.h>
#include "PEANUT_MODEL.h"
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

using Eigen::Map;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using Rcpp::as;


template <typename T> int sign(T val) {
    return (T(0) < val) - (val < T(0));
}

// [[Rcpp::export]]
List peanut_transition(NumericVector a_lin,NumericVector a_loglin,NumericVector a_plin, NumericMatrix x_lin, NumericMatrix x_loglin, NumericMatrix x_plin,int fir,string modelform,int ntime, NumericVector include_bool, double lr, int maxiter, int halfmax, double epsilon, double dbeta_max, double deriv_epsilon, int batch_size){
    //----------------------------------------------------------------------------------------------------------------//
    //
    Map<VectorXd> beta_lin(as<Map<VectorXd> >(a_lin));
    Map<VectorXd> beta_loglin(as<Map<VectorXd> >(a_loglin));
    Map<VectorXd> beta_plin(as<Map<VectorXd> >(a_plin));
    const Map<MatrixXd> df_lin(as<Map<MatrixXd> >(x_lin));
    const Map<MatrixXd> df_loglin(as<Map<MatrixXd> >(x_loglin));
    const Map<MatrixXd> df_plin(as<Map<MatrixXd> >(x_plin));
    // Converts from Rcpp types to efficient Eigen types
    //----------------------------------------------------------------------------------------------------------------//
    List res = LogLik_PEANUT(beta_lin,beta_loglin,beta_plin,df_lin,df_loglin,df_plin,fir,modelform,ntime,include_bool, lr, maxiter, halfmax, epsilon, dbeta_max, deriv_epsilon, batch_size);
    //----------------------------------------------------------------------------------------------------------------//
    return res;
}

List LogLik_PEANUT( VectorXd beta_linT,VectorXd beta_loglinT,VectorXd beta_plinT,MatrixXd df_lin,MatrixXd df_loglin,MatrixXd df_plin,int fir,string modelform,int ntime,NumericVector include_bool, double lr, int maxiter, int halfmax, double epsilon, double dbeta_max, double deriv_epsilon, int batch_size){
    //    cout << "start Risk" << endl;
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
    double totem = df_loglin.rows();//precalculates how many rows
    //
    vector <int> batch_cols;
    int i_temp=0;
    int j_temp=0;
    batch_cols.push_back(0);
    while (i_temp<totem-batch_size){
        if (j_temp==batch_size){
            batch_cols.push_back(i_temp);
        }
        i_temp++;
        j_temp++;
    }
    if (totem-i_temp>batch_size/2){
        batch_cols.push_back(i_temp);
    }
    batch_cols.push_back(totem-1);
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
    // To Start, needs to seperate the derivative terms
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
    end_point = high_resolution_clock::now();
    end = time_point_cast<microseconds>(end_point).time_since_epoch().count();
    cout <<"df100 "<<(end-start)<<" "<<0<<" "<<0<<" "<<-1<<",Prep_T"<<endl;
//    cout << T0 << endl;
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
            #pragma omp parallel for num_threads(nthreads) collapse(2)
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
        #pragma omp parallel for num_threads(nthreads) collapse(2)
        for (int bl=0;bl<batch_cols.size()-1;bl++){
            for (int ijk=0;ijk<totalnum*(totalnum+1)/2;ijk++){
                int ij = 0;
                int jk = ijk;
                while (jk>ij){
                    ij++;
                    jk-=ij;
                }
                Rdd.block(batch_cols[bl],ij*totalnum+jk,batch_cols[bl+1]-batch_cols[bl],1) = Tdd0.block(batch_cols[bl],jk,batch_cols[bl+1]-batch_cols[bl],1).array() * T0.block(batch_cols[bl],jk,batch_cols[bl+1]-batch_cols[bl],1).array().pow(-1).array() * Rd.block(batch_cols[bl],ij,batch_cols[bl+1]-batch_cols[bl],1).array();
                Rdd.block(batch_cols[bl],jk*totalnum+ij,batch_cols[bl+1]-batch_cols[bl],1) = Tdd0.block(batch_cols[bl],ij,batch_cols[bl+1]-batch_cols[bl],1).array() * T0.block(batch_cols[bl],ij,batch_cols[bl+1]-batch_cols[bl],1).array().pow(-1).array() * Rd.block(batch_cols[bl],jk,batch_cols[bl+1]-batch_cols[bl],1).array();
            }
        }
    } else if (modelform=="GM"){
        //currently isn't implemented, it can be calculated but not optimized the same way
        throw invalid_argument( "GM isn't implemented" );
    } else {
        throw invalid_argument( "Model isn't implemented" );
    }
    //
    // -------------------------------------------------------------------------------------------
    //
    string line;
    ifstream infile("test.txt"); //The file of risk group rows
    vector<string>  RiskGroup(ntime); //vector of strings detailing the rows
    IntegerMatrix RiskFail(ntime,2); //vector giving the event rows
    int j_iter = 0;
    //
    end_point = high_resolution_clock::now();
    end = time_point_cast<microseconds>(end_point).time_since_epoch().count();
    cout<<"df100 "<<(end-start)<<" "<<0<<" "<<0<<" "<<-1<<",Prep_R"<<endl;
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
    #pragma omp parallel for num_threads(nthreads) reduction(vec_double_plus:Ll,Lld,Lldd) collapse(2)
    for (int ijk=0;ijk<totalnum*(totalnum+1)/2;ijk++){
        for (int j=0;j<ntime;j++){
            int ij = 0;
            int jk = ijk;
            while (jk>ij){
                ij++;
                jk-=ij;
            }
            double Rs1 = 0;
            double Rs2 = 0;
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
            int at_risk = 0;
            for (int i = 0; i < InGroup.size()-1; i=i+2){
                Rs1 += R.block(InGroup[i]-1,0,InGroup[i+1]-InGroup[i]+1,1).sum();
                Rs2 += Rd.block(InGroup[i]-1,ij,InGroup[i+1]-InGroup[i]+1,1).sum();
                Rs3 += Rdd.block(InGroup[i]-1,ij*totalnum+jk,InGroup[i+1]-InGroup[i]+1,1).sum();
                at_risk += InGroup[i+1]-InGroup[i]+1;
            }
            //
            MatrixXd Ld = MatrixXd::Zero(dj,3);
            Ld << R.block(RiskFail(j,0),0,dj,1), Rd.block(RiskFail(j,0),ij,dj,1) ,Rdd.block(RiskFail(j,0),ij*totalnum+jk,dj,1);//sum of risks in group
            //
            MatrixXd Ldm = MatrixXd::Zero(dj,3);
            for (int i = 0; i < dj; i++){ //adds in the efron approximation terms
                Ldm.row(i) = (-double(i) / double(dj)) * Ld.colwise().sum().array();
            }
            Ldm.col(0) = Ldm.col(0).array() + Rs1;
            Ldm.col(1) = Ldm.col(1).array() + Rs2;
            Ldm.col(2) = Ldm.col(2).array() + Rs3;
            //
            MatrixXd temp1 = MatrixXd::Zero(Ld.rows(),1);
            temp1 = Ld.col(0).array().log();
            double Ld1 =  (temp1.array().isFinite()).select(temp1,0).sum();
            temp1 = Ld.col(1).array() * (Ld.col(0).array().pow(-1));
            double Ld2 = (temp1.array().isFinite()).select(temp1,0).sum();
            temp1 = Ld.col(2).array() * (Ld.col(0).array().pow(-1)) - temp1.array().pow(2);
            double Ld3 = (temp1.array().isFinite()).select(temp1,0).sum();
            //
            temp1 = Ldm.col(0).array().log();
            Rs1 =  (temp1.array().isFinite()).select(temp1,0).sum();
            temp1 = Ldm.col(1).array() * (Ldm.col(0).array().pow(-1));
            Rs2 = (temp1.array().isFinite()).select(temp1,0).sum();
            temp1 = Ldm.col(2).array() * (Ldm.col(0).array().pow(-1)) - temp1.array().pow(2);
            Rs3 = (temp1.array().isFinite()).select(temp1,0).sum();
            //
            if (ij==jk){
                Ll[ij] += Ld1 - Rs1;
                Lld[ij] += Ld2 - Rs2;
            }
            Lldd[ij*totalnum+jk] += Ld3 - Rs3; //sums the log-likelihood and derivatives
            Lldd[jk*totalnum+ij] += Ld3 - Rs3; //sums the log-likelihood and derivatives
        }
    }
    end_point = high_resolution_clock::now();
    end = time_point_cast<microseconds>(end_point).time_since_epoch().count();
    cout<<"df100 "<<(end-start)<<" "<<0<<" "<<0<<" "<<0<<",Calc"<<endl;
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
    //
    double dbeta=0;
    //
    // --------------------------
    // always starts from intial guess
    // --------------------------
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
            #pragma omp parallel for num_threads(nthreads)
            for (int bl=0;bl<batch_cols.size()-1;bl++){
                if (modelform=="A"){//The risks can be updated instead of recalculated for the models/terms allowed
                    if (tform[ind0] == "lin"){
                        R.block(batch_cols[bl],0,batch_cols[bl+1]-batch_cols[bl],1) = R.block(batch_cols[bl],0,batch_cols[bl+1]-batch_cols[bl],1).array() + (beta_c - beta_p) * df0.col(ind0).array();
                    } else if (tform[ind0] == "plin"){
                        R.block(batch_cols[bl],0,batch_cols[bl+1]-batch_cols[bl],1) = R.block(batch_cols[bl],0,batch_cols[bl+1]-batch_cols[bl],1).array() + (beta_c - beta_p) * df0.col(ind0).array();
                    } else{
                        R.block(batch_cols[bl],0,batch_cols[bl+1]-batch_cols[bl],1) = R.block(batch_cols[bl],0,batch_cols[bl+1]-batch_cols[bl],1).array() + (beta_c * df0.col(ind0)).array().exp() - (beta_p * df0.col(ind0)).array().exp();
                    }
                } else if (modelform=="PA"){
                    if (tform[ind0] == "lin"){
                        if (i==fir){
                            R.block(batch_cols[bl],0,batch_cols[bl+1]-batch_cols[bl],1) = R.block(batch_cols[bl],0,batch_cols[bl+1]-batch_cols[bl],1) * beta_c / beta_p;
                            T0.col(ind0) = T0.col(ind0) * beta_c / beta_p;
                        } else {
                            R.block(batch_cols[bl],0,batch_cols[bl+1]-batch_cols[bl],1) = R.block(batch_cols[bl],0,batch_cols[bl+1]-batch_cols[bl],1).array() + (beta_c - beta_p) * df0.col(ind0).array() * T0.col(ind0).array();
                            Te = Te.array() + (beta_c - beta_p) * df0.col(ind0).array();
                        }
                    } else if (tform[ind0] == "plin"){
                        if (i==fir){
                            R.block(batch_cols[bl],0,batch_cols[bl+1]-batch_cols[bl],1) = (R.block(batch_cols[bl],0,batch_cols[bl+1]-batch_cols[bl],1).array() - Te.array()) * beta_c / beta_p + Te.array();
                            T0.col(ind0) = (T0.col(ind0).array() - 1) * beta_c / beta_p + 1;
                        } else {
                            R.block(batch_cols[bl],0,batch_cols[bl+1]-batch_cols[bl],1) = R.block(batch_cols[bl],0,batch_cols[bl+1]-batch_cols[bl],1).array() + (beta_c - beta_p) * df0.col(ind0).array() * T0.col(ind0).array();
                            Te = Te.array() + (beta_c - beta_p) * df0.col(ind0).array();
                        }
                    } else{
                        if (i==fir){
                            R.block(batch_cols[bl],0,batch_cols[bl+1]-batch_cols[bl],1) = R.block(batch_cols[bl],0,batch_cols[bl+1]-batch_cols[bl],1).array() * ((beta_c - beta_p) * df0.col(ind0)).array().exp();
                        } else {
                            R.block(batch_cols[bl],0,batch_cols[bl+1]-batch_cols[bl],1) = R.block(batch_cols[bl],0,batch_cols[bl+1]-batch_cols[bl],1).array() + ((beta_c * df0.col(ind0)).array().exp() - (beta_p * df0.col(ind0)).array().exp()).array() * T0.col(ind0).array();
                            Te = Te.array() + ((beta_c * df0.col(ind0)).array().exp() - (beta_p * df0.col(ind0)).array().exp()).array();
                        }
                    }
                } else if (modelform=="PAE"){
                    if (tform[ind0] == "lin"){
                        if (i==fir){
                            R.block(batch_cols[bl],0,batch_cols[bl+1]-batch_cols[bl],1) = R.block(batch_cols[bl],0,batch_cols[bl+1]-batch_cols[bl],1) * beta_c / beta_p;
                            T0.col(ind0) = T0.col(ind0) * beta_c / beta_p;
                        } else {
                            R.block(batch_cols[bl],0,batch_cols[bl+1]-batch_cols[bl],1) = R.block(batch_cols[bl],0,batch_cols[bl+1]-batch_cols[bl],1).array() + (beta_c - beta_p) * df0.col(ind0).array() * T0.col(ind0).array();
                            Te = Te.array() + (beta_c - beta_p) * df0.col(ind0).array();
                        }
                    } else if (tform[ind0] == "plin"){
                        if (i==fir){
                            R.block(batch_cols[bl],0,batch_cols[bl+1]-batch_cols[bl],1) = (R.block(batch_cols[bl],0,batch_cols[bl+1]-batch_cols[bl],1).array() - Te.array()) * beta_c / beta_p + Te.array();
                            T0.col(ind0) = (T0.col(ind0).array() - 1) * beta_c / beta_p + 1;
                        } else {
                            R.block(batch_cols[bl],0,batch_cols[bl+1]-batch_cols[bl],1) = R.block(batch_cols[bl],0,batch_cols[bl+1]-batch_cols[bl],1).array() + (beta_c - beta_p) * df0.col(ind0).array() * T0.col(ind0).array();
                            Te = Te.array() + (beta_c - beta_p) * df0.col(ind0).array();
                        }
                    } else{
                        if (i==fir){
                            R.block(batch_cols[bl],0,batch_cols[bl+1]-batch_cols[bl],1) = R.block(batch_cols[bl],0,batch_cols[bl+1]-batch_cols[bl],1).array() * ((beta_c - beta_p) * df0.col(ind0)).array().exp();
                        } else {
                            R.block(batch_cols[bl],0,batch_cols[bl+1]-batch_cols[bl],1) = R.block(batch_cols[bl],0,batch_cols[bl+1]-batch_cols[bl],1).array() + ((beta_c * df0.col(ind0)).array().exp() - (beta_p * df0.col(ind0)).array().exp()) * T0.array();
                            Te = Te.array() + ((beta_c * df0.col(ind0)).array().exp() - (beta_p * df0.col(ind0)).array().exp()).array();
                        }
                    }
                }else if (modelform=="M"){
                    if (tform[ind0] == "lin"){
                        R.block(batch_cols[bl],0,batch_cols[bl+1]-batch_cols[bl],1) = R.block(batch_cols[bl],0,batch_cols[bl+1]-batch_cols[bl],1) * beta_c / beta_p;
                    } else if (tform[ind0] == "plin"){
                        R.block(batch_cols[bl],0,batch_cols[bl+1]-batch_cols[bl],1) = R.block(batch_cols[bl],0,batch_cols[bl+1]-batch_cols[bl],1).array() * (1 + beta_c * df0.col(ind0).array()) / (1 + beta_p * df0.col(ind0).array());
                    } else{
                        R.block(batch_cols[bl],0,batch_cols[bl+1]-batch_cols[bl],1) = R.block(batch_cols[bl],0,batch_cols[bl+1]-batch_cols[bl],1).array() * ((beta_c - beta_p) * df0.block(batch_cols[bl],0,batch_cols[bl+1]-batch_cols[bl],1)).array().exp();
                    }
                }else if (modelform=="GM"){
                    throw invalid_argument( "GM isn't implemented" );
                } else {
                    throw invalid_argument( "Model isn't implemented" );
                }
            }
            #pragma omp parallel for num_threads(nthreads) collapse(2)
            for (int bl=0;bl<batch_cols.size()-1;bl++){
                for (int ijk=0;ijk<totalnum*(totalnum+1)/2;ijk++){
                    int ij = 0;
                    int jk = ijk;
                    while (jk>ij){
                        ij++;
                        jk-=ij;
                    }
                    if (modelform=="A"){
                        if (ij==jk){
                            if (tform[ij] == "lin"){
                                ;
                            } else if (tform[ij] == "plin"){
                                ;
                            } else{
                                Rd.block(batch_cols[bl],ij,batch_cols[bl+1]-batch_cols[bl],1) = df0.block(batch_cols[bl],ij,batch_cols[bl+1]-batch_cols[bl],1).array() * R.block(batch_cols[bl],0,batch_cols[bl+1]-batch_cols[bl],1).array();
                                Rdd.block(batch_cols[bl],ij*totalnum+jk,batch_cols[bl+1]-batch_cols[bl],1) = df0.block(batch_cols[bl],ij,batch_cols[bl+1]-batch_cols[bl],1).array() * Rd.block(batch_cols[bl],ij,batch_cols[bl+1]-batch_cols[bl],1).array();
                            }
                        }
                    } else if (modelform=="PA"){
                        if (tform[ij] == "lin"){
                            if (ij==jk){
                                if (ij==fir){
                                    Rd.block(batch_cols[bl],ij,batch_cols[bl+1]-batch_cols[bl],1) = df0.block(batch_cols[bl],ij,batch_cols[bl+1]-batch_cols[bl],1).array() * Te.array();
                                } else {
                                    Rd.block(batch_cols[bl],ij,batch_cols[bl+1]-batch_cols[bl],1) = df0.block(batch_cols[bl],ij,batch_cols[bl+1]-batch_cols[bl],1).array() * T0.col(fir).array();
                                }
                            }
                        } else if (tform[ij] == "plin"){
                            if (ij==jk){
                                if (ij==fir){
                                    Rd.block(batch_cols[bl],ij,batch_cols[bl+1]-batch_cols[bl],1) = df0.block(batch_cols[bl],ij,batch_cols[bl+1]-batch_cols[bl],1).array() * Te.array();
                                } else {
                                    Rd.block(batch_cols[bl],ij,batch_cols[bl+1]-batch_cols[bl],1) = df0.block(batch_cols[bl],ij,batch_cols[bl+1]-batch_cols[bl],1).array() * T0.col(fir).array();
                                }
                            }
                        } else if (tform[jk] == "loglin"){
                                if (ij==fir){
                                    Rd.block(batch_cols[bl],ij,batch_cols[bl+1]-batch_cols[bl],1) = df0.block(batch_cols[bl],ij,batch_cols[bl+1]-batch_cols[bl],1).array() * R.block(batch_cols[bl],0,batch_cols[bl+1]-batch_cols[bl],1).array();
                                } else {
                                    Rd.block(batch_cols[bl],ij,batch_cols[bl+1]-batch_cols[bl],1) =T0.col(fir).array() * df0.block(batch_cols[bl],ij,batch_cols[bl+1]-batch_cols[bl],1).array() * ((beta_0[ij]) * df0.block(batch_cols[bl],ij,batch_cols[bl+1]-batch_cols[bl],1)).array().exp();
                                }
                                Rdd.block(batch_cols[bl],ij*totalnum+jk,batch_cols[bl+1]-batch_cols[bl],1) = df0.block(batch_cols[bl],jk,batch_cols[bl+1]-batch_cols[bl],1).array() * Rd.block(batch_cols[bl],ij,batch_cols[bl+1]-batch_cols[bl],1).array();
                                if (ij!=jk){
                                    Rdd.block(batch_cols[bl],jk*totalnum+ij,batch_cols[bl+1]-batch_cols[bl],1) = df0.block(batch_cols[bl],ij,batch_cols[bl+1]-batch_cols[bl],1).array() * Rd.block(batch_cols[bl],jk,batch_cols[bl+1]-batch_cols[bl],1).array();
                                }
                            }
                    } else if (modelform=="PAE"){
                        if (tform[ij] == "lin"){
                            if (ij==jk){
                                if (ij==fir){
                                    Rd.block(batch_cols[bl],ij,batch_cols[bl+1]-batch_cols[bl],1) = df0.block(batch_cols[bl],ij,batch_cols[bl+1]-batch_cols[bl],1).array() * Te.array();
                                } else {
                                    Rd.block(batch_cols[bl],ij,batch_cols[bl+1]-batch_cols[bl],1) = df0.block(batch_cols[bl],ij,batch_cols[bl+1]-batch_cols[bl],1).array() * T0.col(fir).array();
                                }
                            }
                        } else if (tform[ij] == "plin"){
                            if (ij==jk){
                                if (ij==fir){
                                    Rd.block(batch_cols[bl],ij,batch_cols[bl+1]-batch_cols[bl],1) = df0.block(batch_cols[bl],ij,batch_cols[bl+1]-batch_cols[bl],1).array() * Te.array();
                                } else {
                                    Rd.block(batch_cols[bl],ij,batch_cols[bl+1]-batch_cols[bl],1) = df0.block(batch_cols[bl],ij,batch_cols[bl+1]-batch_cols[bl],1).array() * T0.block(batch_cols[bl],ij,batch_cols[bl+1]-batch_cols[bl],1).array();
                                }
                            }
                        } else{
                            if (ij==fir){
                                Rd.block(batch_cols[bl],ij,batch_cols[bl+1]-batch_cols[bl],1) = df0.block(batch_cols[bl],ij,batch_cols[bl+1]-batch_cols[bl],1).array() * R.block(batch_cols[bl],0,batch_cols[bl+1]-batch_cols[bl],1).array();
                                Rdd.block(batch_cols[bl],ij,batch_cols[bl+1]-batch_cols[bl],1) = df0.block(batch_cols[bl],ij,batch_cols[bl+1]-batch_cols[bl],1).array() * Rd.block(batch_cols[bl],ij,batch_cols[bl+1]-batch_cols[bl],1).array();
                            } else if (tform[jk] == "loglin"){
                                    if (ij==fir){
                                        Rd.block(batch_cols[bl],ij,batch_cols[bl+1]-batch_cols[bl],1) = df0.block(batch_cols[bl],ij,batch_cols[bl+1]-batch_cols[bl],1).array() * R.block(batch_cols[bl],0,batch_cols[bl+1]-batch_cols[bl],1).array();
                                    } else {
                                        Rd.block(batch_cols[bl],ij,batch_cols[bl+1]-batch_cols[bl],1) =T0.col(fir).array() * df0.block(batch_cols[bl],ij,batch_cols[bl+1]-batch_cols[bl],1).array() * ((beta_0[ij]) * df0.block(batch_cols[bl],ij,batch_cols[bl+1]-batch_cols[bl],1)).array().exp();
                                    }
                                    Rdd.block(batch_cols[bl],ij*totalnum+jk,batch_cols[bl+1]-batch_cols[bl],1) = df0.block(batch_cols[bl],jk,batch_cols[bl+1]-batch_cols[bl],1).array() * Rd.block(batch_cols[bl],ij,batch_cols[bl+1]-batch_cols[bl],1).array();
                                    if (ij!=jk){
                                        Rdd.block(batch_cols[bl],jk*totalnum+ij,batch_cols[bl+1]-batch_cols[bl],1) = df0.block(batch_cols[bl],ij,batch_cols[bl+1]-batch_cols[bl],1).array() * Rd.block(batch_cols[bl],jk,batch_cols[bl+1]-batch_cols[bl],1).array();
                                    }
                                }
                            }
                    }else if (modelform=="M"){
                        if (tform[ij] == "lin"){
                            if (ij==jk){
                                Rd.block(batch_cols[bl],ij,batch_cols[bl+1]-batch_cols[bl],1) = R.block(batch_cols[bl],0,batch_cols[bl+1]-batch_cols[bl],1) / beta_0[ij];
                            }
                        } else if (tform[ij] == "plin"){
                            if (ij==jk){
                                Rd.block(batch_cols[bl],ij,batch_cols[bl+1]-batch_cols[bl],1) = R.block(batch_cols[bl],0,batch_cols[bl+1]-batch_cols[bl],1).array() * (df0.block(batch_cols[bl],ij,batch_cols[bl+1]-batch_cols[bl],1).array()) / (1 + beta_0[ij] * df0.block(batch_cols[bl],ij,batch_cols[bl+1]-batch_cols[bl],1).array());
                            }
                        } else{
                            if (ij==jk){
                                Rd.block(batch_cols[bl],ij,batch_cols[bl+1]-batch_cols[bl],1) = df0.block(batch_cols[bl],ij,batch_cols[bl+1]-batch_cols[bl],1).array() * R.block(batch_cols[bl],0,batch_cols[bl+1]-batch_cols[bl],1).array();
                            }
                            if (tform[jk] == "loglin"){
                                Rdd.block(batch_cols[bl],ij*totalnum+jk,batch_cols[bl+1]-batch_cols[bl],1) = df0.block(batch_cols[bl],ij,batch_cols[bl+1]-batch_cols[bl],1).array() * Rd.block(batch_cols[bl],jk,batch_cols[bl+1]-batch_cols[bl],1).array();
                                if (ij!=jk){
                                    Rdd.block(batch_cols[bl],jk*totalnum+ij,batch_cols[bl+1]-batch_cols[bl],1) = df0.block(batch_cols[bl],jk,batch_cols[bl+1]-batch_cols[bl],1).array() * Rd.block(batch_cols[bl],ij,batch_cols[bl+1]-batch_cols[bl],1).array();
                                }
                            }
                        }
                    } else {
                        throw invalid_argument( "GM isn't implemented" );
                    }
                }
            }
            end_point = high_resolution_clock::now();
            end = time_point_cast<microseconds>(end_point).time_since_epoch().count();
            cout<<"df100 "<<(end-start)<<" "<<halves<<" "<<iteration<<" "<<ind0<<",Update_R"<<endl;
            fill(Ll.begin(), Ll.end(), 0.0);
            fill(Lld.begin(), Lld.end(), 0.0);
            fill(Lldd.begin(), Lldd.end(), 0.0);
            #pragma omp parallel for num_threads(nthreads) reduction(vec_double_plus:Ll,Lld,Lldd) collapse(2)
            for (int j=0;j<ntime;j++){
                for (int ijk=0;ijk<totalnum*(totalnum+1)/2;ijk++){
                    int ij = 0;
                    int jk = ijk;
                    while (jk>ij){
                        ij++;
                        jk-=ij;
                    }
                    double Rs1 = 0;
                    double Rs2 = 0;
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
                    int at_risk = 0;
                    for (int i = 0; i < InGroup.size()-1; i=i+2){
                        Rs1 += R.block(InGroup[i]-1,0,InGroup[i+1]-InGroup[i]+1,1).sum();
                        Rs2 += Rd.block(InGroup[i]-1,ij,InGroup[i+1]-InGroup[i]+1,1).sum();
                        Rs3 += Rdd.block(InGroup[i]-1,ij*totalnum+jk,InGroup[i+1]-InGroup[i]+1,1).sum();
                        at_risk += InGroup[i+1]-InGroup[i]+1;
                    }
                    //
                    MatrixXd Ld = MatrixXd::Zero(dj,3);
                    Ld << R.block(RiskFail(j,0),0,dj,1), Rd.block(RiskFail(j,0),ij,dj,1) ,Rdd.block(RiskFail(j,0),ij*totalnum+jk,dj,1);//sum of risks in group
                    //
                    MatrixXd Ldm = MatrixXd::Zero(dj,3);
                    for (int i = 0; i < dj; i++){ //adds in the efron approximation terms
                        Ldm.row(i) = (-double(i) / double(dj)) * Ld.colwise().sum().array();
                    }
                    Ldm.col(0) = Ldm.col(0).array() + Rs1;
                    Ldm.col(1) = Ldm.col(1).array() + Rs2;
                    Ldm.col(2) = Ldm.col(2).array() + Rs3;
                    //
                    MatrixXd temp1 = MatrixXd::Zero(Ld.rows(),1);
                    temp1 = Ld.col(0).array().log();
                    double Ld1 =  (temp1.array().isFinite()).select(temp1,0).sum();
                    temp1 = Ld.col(1).array() * (Ld.col(0).array().pow(-1));
                    double Ld2 = (temp1.array().isFinite()).select(temp1,0).sum();
                    temp1 = Ld.col(2).array() * (Ld.col(0).array().pow(-1)) - temp1.array().pow(2);
                    double Ld3 = (temp1.array().isFinite()).select(temp1,0).sum();
                    //
                    temp1 = Ldm.col(0).array().log();
                    Rs1 =  (temp1.array().isFinite()).select(temp1,0).sum();
                    temp1 = Ldm.col(1).array() * (Ldm.col(0).array().pow(-1));
                    Rs2 = (temp1.array().isFinite()).select(temp1,0).sum();
                    temp1 = Ldm.col(2).array() * (Ldm.col(0).array().pow(-1)) - temp1.array().pow(2);
                    Rs3 = (temp1.array().isFinite()).select(temp1,0).sum();
                    //
                    if (ij==jk){
                        Ll[ij] += Ld1 - Rs1;
                        Lld[ij] += Ld2 - Rs2;
                    }
                    Lldd[ij*totalnum+jk] += Ld3 - Rs3; //sums the log-likelihood and derivatives
                    Lldd[jk*totalnum+ij] += Ld3 - Rs3; //sums the log-likelihood and derivatives
                }
            }
            if (Ll[ind0] <= Ll_best){
                dbeta = dbeta / 2.0;
            } else{
                beta_best = beta_c;
            }
            end_point = high_resolution_clock::now();
            end = time_point_cast<microseconds>(end_point).time_since_epoch().count();
            beta_p = beta_c;
            cout<<"df100 "<<(end-start)<<" "<<halves<<" "<<iteration<<" "<<ind0<<",Update_calc"<<endl;
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
            #pragma omp parallel for num_threads(nthreads)
            for (int bl=0;bl<batch_cols.size()-1;bl++){
                if (modelform=="A"){//The risks can be updated instead of recalculated for the models/terms allowed
                    if (tform[ind0] == "lin"){
                        R.block(batch_cols[bl],0,batch_cols[bl+1]-batch_cols[bl],1) = R.block(batch_cols[bl],0,batch_cols[bl+1]-batch_cols[bl],1).array() + (beta_c - beta_p) * df0.col(ind0).array();
                    } else if (tform[ind0] == "plin"){
                        R.block(batch_cols[bl],0,batch_cols[bl+1]-batch_cols[bl],1) = R.block(batch_cols[bl],0,batch_cols[bl+1]-batch_cols[bl],1).array() + (beta_c - beta_p) * df0.col(ind0).array();
                    } else{
                        R.block(batch_cols[bl],0,batch_cols[bl+1]-batch_cols[bl],1) = R.block(batch_cols[bl],0,batch_cols[bl+1]-batch_cols[bl],1).array() + (beta_c * df0.col(ind0)).array().exp() - (beta_p * df0.col(ind0)).array().exp();
                    }
                } else if (modelform=="PA"){
                    if (tform[ind0] == "lin"){
                        if (i==fir){
                            R.block(batch_cols[bl],0,batch_cols[bl+1]-batch_cols[bl],1) = R.block(batch_cols[bl],0,batch_cols[bl+1]-batch_cols[bl],1) * beta_c / beta_p;
                            T0.col(ind0) = T0.col(ind0) * beta_c / beta_p;
                        } else {
                            R.block(batch_cols[bl],0,batch_cols[bl+1]-batch_cols[bl],1) = R.block(batch_cols[bl],0,batch_cols[bl+1]-batch_cols[bl],1).array() + (beta_c - beta_p) * df0.col(ind0).array() * T0.col(ind0).array();
                            Te = Te.array() + (beta_c - beta_p) * df0.col(ind0).array();
                        }
                    } else if (tform[ind0] == "plin"){
                        if (i==fir){
                            R.block(batch_cols[bl],0,batch_cols[bl+1]-batch_cols[bl],1) = (R.block(batch_cols[bl],0,batch_cols[bl+1]-batch_cols[bl],1).array() - Te.array()) * beta_c / beta_p + Te.array();
                            T0.col(ind0) = (T0.col(ind0).array() - 1) * beta_c / beta_p + 1;
                        } else {
                            R.block(batch_cols[bl],0,batch_cols[bl+1]-batch_cols[bl],1) = R.block(batch_cols[bl],0,batch_cols[bl+1]-batch_cols[bl],1).array() + (beta_c - beta_p) * df0.col(ind0).array() * T0.col(ind0).array();
                            Te = Te.array() + (beta_c - beta_p) * df0.col(ind0).array();
                        }
                    } else{
                        if (i==fir){
                            R.block(batch_cols[bl],0,batch_cols[bl+1]-batch_cols[bl],1) = R.block(batch_cols[bl],0,batch_cols[bl+1]-batch_cols[bl],1).array() * ((beta_c - beta_p) * df0.col(ind0)).array().exp();
                        } else {
                            R.block(batch_cols[bl],0,batch_cols[bl+1]-batch_cols[bl],1) = R.block(batch_cols[bl],0,batch_cols[bl+1]-batch_cols[bl],1).array() + ((beta_c * df0.col(ind0)).array().exp() - (beta_p * df0.col(ind0)).array().exp()).array() * T0.col(ind0).array();
                            Te = Te.array() + ((beta_c * df0.col(ind0)).array().exp() - (beta_p * df0.col(ind0)).array().exp()).array();
                        }
                    }
                } else if (modelform=="PAE"){
                    if (tform[ind0] == "lin"){
                        if (i==fir){
                            R.block(batch_cols[bl],0,batch_cols[bl+1]-batch_cols[bl],1) = R.block(batch_cols[bl],0,batch_cols[bl+1]-batch_cols[bl],1) * beta_c / beta_p;
                            T0.col(ind0) = T0.col(ind0) * beta_c / beta_p;
                        } else {
                            R.block(batch_cols[bl],0,batch_cols[bl+1]-batch_cols[bl],1) = R.block(batch_cols[bl],0,batch_cols[bl+1]-batch_cols[bl],1).array() + (beta_c - beta_p) * df0.col(ind0).array() * T0.col(ind0).array();
                            Te = Te.array() + (beta_c - beta_p) * df0.col(ind0).array();
                        }
                    } else if (tform[ind0] == "plin"){
                        if (i==fir){
                            R.block(batch_cols[bl],0,batch_cols[bl+1]-batch_cols[bl],1) = (R.block(batch_cols[bl],0,batch_cols[bl+1]-batch_cols[bl],1).array() - Te.array()) * beta_c / beta_p + Te.array();
                            T0.col(ind0) = (T0.col(ind0).array() - 1) * beta_c / beta_p + 1;
                        } else {
                            R.block(batch_cols[bl],0,batch_cols[bl+1]-batch_cols[bl],1) = R.block(batch_cols[bl],0,batch_cols[bl+1]-batch_cols[bl],1).array() + (beta_c - beta_p) * df0.col(ind0).array() * T0.col(ind0).array();
                            Te = Te.array() + (beta_c - beta_p) * df0.col(ind0).array();
                        }
                    } else{
                        if (i==fir){
                            R.block(batch_cols[bl],0,batch_cols[bl+1]-batch_cols[bl],1) = R.block(batch_cols[bl],0,batch_cols[bl+1]-batch_cols[bl],1).array() * ((beta_c - beta_p) * df0.col(ind0)).array().exp();
                        } else {
                            R.block(batch_cols[bl],0,batch_cols[bl+1]-batch_cols[bl],1) = R.block(batch_cols[bl],0,batch_cols[bl+1]-batch_cols[bl],1).array() + ((beta_c * df0.col(ind0)).array().exp() - (beta_p * df0.col(ind0)).array().exp()) * T0.array();
                            Te = Te.array() + ((beta_c * df0.col(ind0)).array().exp() - (beta_p * df0.col(ind0)).array().exp()).array();
                        }
                    }
                }else if (modelform=="M"){
                    if (tform[ind0] == "lin"){
                        R.block(batch_cols[bl],0,batch_cols[bl+1]-batch_cols[bl],1) = R.block(batch_cols[bl],0,batch_cols[bl+1]-batch_cols[bl],1) * beta_c / beta_p;
                    } else if (tform[ind0] == "plin"){
                        R.block(batch_cols[bl],0,batch_cols[bl+1]-batch_cols[bl],1) = R.block(batch_cols[bl],0,batch_cols[bl+1]-batch_cols[bl],1).array() * (1 + beta_c * df0.col(ind0).array()) / (1 + beta_p * df0.col(ind0).array());
                    } else{
                        R.block(batch_cols[bl],0,batch_cols[bl+1]-batch_cols[bl],1) = R.block(batch_cols[bl],0,batch_cols[bl+1]-batch_cols[bl],1).array() * ((beta_c - beta_p) * df0.block(batch_cols[bl],0,batch_cols[bl+1]-batch_cols[bl],1)).array().exp();
                    }
                }else if (modelform=="GM"){
                    throw invalid_argument( "GM isn't implemented" );
                } else {
                    throw invalid_argument( "Model isn't implemented" );
                }
            }
            #pragma omp parallel for num_threads(nthreads) collapse(2)
            for (int bl=0;bl<batch_cols.size()-1;bl++){
                for (int ijk=0;ijk<totalnum*(totalnum+1)/2;ijk++){
                    int ij = 0;
                    int jk = ijk;
                    while (jk>ij){
                        ij++;
                        jk-=ij;
                    }
                    if (modelform=="A"){
                        if (ij==jk){
                            if (tform[ij] == "lin"){
                                ;
                            } else if (tform[ij] == "plin"){
                                ;
                            } else{
                                Rd.block(batch_cols[bl],ij,batch_cols[bl+1]-batch_cols[bl],1) = df0.block(batch_cols[bl],ij,batch_cols[bl+1]-batch_cols[bl],1).array() * R.block(batch_cols[bl],0,batch_cols[bl+1]-batch_cols[bl],1).array();
                                Rdd.block(batch_cols[bl],ij*totalnum+jk,batch_cols[bl+1]-batch_cols[bl],1) = df0.block(batch_cols[bl],ij,batch_cols[bl+1]-batch_cols[bl],1).array() * Rd.block(batch_cols[bl],ij,batch_cols[bl+1]-batch_cols[bl],1).array();
                            }
                        }
                    } else if (modelform=="PA"){
                        if (tform[ij] == "lin"){
                            if (ij==jk){
                                if (ij==fir){
                                    Rd.block(batch_cols[bl],ij,batch_cols[bl+1]-batch_cols[bl],1) = df0.block(batch_cols[bl],ij,batch_cols[bl+1]-batch_cols[bl],1).array() * Te.array();
                                } else {
                                    Rd.block(batch_cols[bl],ij,batch_cols[bl+1]-batch_cols[bl],1) = df0.block(batch_cols[bl],ij,batch_cols[bl+1]-batch_cols[bl],1).array() * T0.col(fir).array();
                                }
                            }
                        } else if (tform[ij] == "plin"){
                            if (ij==jk){
                                if (ij==fir){
                                    Rd.block(batch_cols[bl],ij,batch_cols[bl+1]-batch_cols[bl],1) = df0.block(batch_cols[bl],ij,batch_cols[bl+1]-batch_cols[bl],1).array() * Te.array();
                                } else {
                                    Rd.block(batch_cols[bl],ij,batch_cols[bl+1]-batch_cols[bl],1) = df0.block(batch_cols[bl],ij,batch_cols[bl+1]-batch_cols[bl],1).array() * T0.col(fir).array();
                                }
                            }
                        } else if (tform[jk] == "loglin"){
                                if (ij==fir){
                                    Rd.block(batch_cols[bl],ij,batch_cols[bl+1]-batch_cols[bl],1) = df0.block(batch_cols[bl],ij,batch_cols[bl+1]-batch_cols[bl],1).array() * R.block(batch_cols[bl],0,batch_cols[bl+1]-batch_cols[bl],1).array();
                                } else {
                                    Rd.block(batch_cols[bl],ij,batch_cols[bl+1]-batch_cols[bl],1) =T0.col(fir).array() * df0.block(batch_cols[bl],ij,batch_cols[bl+1]-batch_cols[bl],1).array() * ((beta_0[ij]) * df0.block(batch_cols[bl],ij,batch_cols[bl+1]-batch_cols[bl],1)).array().exp();
                                }
                                Rdd.block(batch_cols[bl],ij*totalnum+jk,batch_cols[bl+1]-batch_cols[bl],1) = df0.block(batch_cols[bl],jk,batch_cols[bl+1]-batch_cols[bl],1).array() * Rd.block(batch_cols[bl],ij,batch_cols[bl+1]-batch_cols[bl],1).array();
                                if (ij!=jk){
                                    Rdd.block(batch_cols[bl],jk*totalnum+ij,batch_cols[bl+1]-batch_cols[bl],1) = df0.block(batch_cols[bl],ij,batch_cols[bl+1]-batch_cols[bl],1).array() * Rd.block(batch_cols[bl],jk,batch_cols[bl+1]-batch_cols[bl],1).array();
                                }
                            }
                    } else if (modelform=="PAE"){
                        if (tform[ij] == "lin"){
                            if (ij==jk){
                                if (ij==fir){
                                    Rd.block(batch_cols[bl],ij,batch_cols[bl+1]-batch_cols[bl],1) = df0.block(batch_cols[bl],ij,batch_cols[bl+1]-batch_cols[bl],1).array() * Te.array();
                                } else {
                                    Rd.block(batch_cols[bl],ij,batch_cols[bl+1]-batch_cols[bl],1) = df0.block(batch_cols[bl],ij,batch_cols[bl+1]-batch_cols[bl],1).array() * T0.col(fir).array();
                                }
                            }
                        } else if (tform[ij] == "plin"){
                            if (ij==jk){
                                if (ij==fir){
                                    Rd.block(batch_cols[bl],ij,batch_cols[bl+1]-batch_cols[bl],1) = df0.block(batch_cols[bl],ij,batch_cols[bl+1]-batch_cols[bl],1).array() * Te.array();
                                } else {
                                    Rd.block(batch_cols[bl],ij,batch_cols[bl+1]-batch_cols[bl],1) = df0.block(batch_cols[bl],ij,batch_cols[bl+1]-batch_cols[bl],1).array() * T0.block(batch_cols[bl],ij,batch_cols[bl+1]-batch_cols[bl],1).array();
                                }
                            }
                        } else{
                            if (ij==fir){
                                Rd.block(batch_cols[bl],ij,batch_cols[bl+1]-batch_cols[bl],1) = df0.block(batch_cols[bl],ij,batch_cols[bl+1]-batch_cols[bl],1).array() * R.block(batch_cols[bl],0,batch_cols[bl+1]-batch_cols[bl],1).array();
                                Rdd.block(batch_cols[bl],ij,batch_cols[bl+1]-batch_cols[bl],1) = df0.block(batch_cols[bl],ij,batch_cols[bl+1]-batch_cols[bl],1).array() * Rd.block(batch_cols[bl],ij,batch_cols[bl+1]-batch_cols[bl],1).array();
                            } else if (tform[jk] == "loglin"){
                                    if (ij==fir){
                                        Rd.block(batch_cols[bl],ij,batch_cols[bl+1]-batch_cols[bl],1) = df0.block(batch_cols[bl],ij,batch_cols[bl+1]-batch_cols[bl],1).array() * R.block(batch_cols[bl],0,batch_cols[bl+1]-batch_cols[bl],1).array();
                                    } else {
                                        Rd.block(batch_cols[bl],ij,batch_cols[bl+1]-batch_cols[bl],1) =T0.col(fir).array() * df0.block(batch_cols[bl],ij,batch_cols[bl+1]-batch_cols[bl],1).array() * ((beta_0[ij]) * df0.block(batch_cols[bl],ij,batch_cols[bl+1]-batch_cols[bl],1)).array().exp();
                                    }
                                    Rdd.block(batch_cols[bl],ij*totalnum+jk,batch_cols[bl+1]-batch_cols[bl],1) = df0.block(batch_cols[bl],jk,batch_cols[bl+1]-batch_cols[bl],1).array() * Rd.block(batch_cols[bl],ij,batch_cols[bl+1]-batch_cols[bl],1).array();
                                    if (ij!=jk){
                                        Rdd.block(batch_cols[bl],jk*totalnum+ij,batch_cols[bl+1]-batch_cols[bl],1) = df0.block(batch_cols[bl],ij,batch_cols[bl+1]-batch_cols[bl],1).array() * Rd.block(batch_cols[bl],jk,batch_cols[bl+1]-batch_cols[bl],1).array();
                                    }
                                }
                            }
                    }else if (modelform=="M"){
                        if (tform[ij] == "lin"){
                            if (ij==jk){
                                Rd.block(batch_cols[bl],ij,batch_cols[bl+1]-batch_cols[bl],1) = R.block(batch_cols[bl],0,batch_cols[bl+1]-batch_cols[bl],1) / beta_0[ij];
                            }
                        } else if (tform[ij] == "plin"){
                            if (ij==jk){
                                Rd.block(batch_cols[bl],ij,batch_cols[bl+1]-batch_cols[bl],1) = R.block(batch_cols[bl],0,batch_cols[bl+1]-batch_cols[bl],1).array() * (df0.block(batch_cols[bl],ij,batch_cols[bl+1]-batch_cols[bl],1).array()) / (1 + beta_0[ij] * df0.block(batch_cols[bl],ij,batch_cols[bl+1]-batch_cols[bl],1).array());
                            }
                        } else{
                            if (ij==jk){
                                Rd.block(batch_cols[bl],ij,batch_cols[bl+1]-batch_cols[bl],1) = df0.block(batch_cols[bl],ij,batch_cols[bl+1]-batch_cols[bl],1).array() * R.block(batch_cols[bl],0,batch_cols[bl+1]-batch_cols[bl],1).array();
                            }
                            if (tform[jk] == "loglin"){
                                Rdd.block(batch_cols[bl],ij*totalnum+jk,batch_cols[bl+1]-batch_cols[bl],1) = df0.block(batch_cols[bl],ij,batch_cols[bl+1]-batch_cols[bl],1).array() * Rd.block(batch_cols[bl],jk,batch_cols[bl+1]-batch_cols[bl],1).array();
                                if (ij!=jk){
                                    Rdd.block(batch_cols[bl],jk*totalnum+ij,batch_cols[bl+1]-batch_cols[bl],1) = df0.block(batch_cols[bl],jk,batch_cols[bl+1]-batch_cols[bl],1).array() * Rd.block(batch_cols[bl],ij,batch_cols[bl+1]-batch_cols[bl],1).array();
                                }
                            }
                        }
                    } else {
                        throw invalid_argument( "GM isn't implemented" );
                    }
                }
            }
            fill(Ll.begin(), Ll.end(), 0.0);
            fill(Lld.begin(), Lld.end(), 0.0);
            fill(Lldd.begin(), Lldd.end(), 0.0);
            end_point = high_resolution_clock::now();
            end = time_point_cast<microseconds>(end_point).time_since_epoch().count();
            cout<<"df100 "<<(end-start)<<" "<<halves<<" "<<iteration<<" "<<ind0<<",Revert"<<endl;
            #pragma omp parallel for num_threads(nthreads) reduction(vec_double_plus:Ll,Lld,Lldd) collapse(2)
            for (int j=0;j<ntime;j++){
                for (int ijk=0;ijk<totalnum*(totalnum+1)/2;ijk++){
                    int ij = 0;
                    int jk = ijk;
                    while (jk>ij){
                        ij++;
                        jk-=ij;
                    }
                   
                    double Rs1 = 0;
                    double Rs2 = 0;
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
                    int at_risk = 0;
                    for (int i = 0; i < InGroup.size()-1; i=i+2){
                        Rs1 += R.block(InGroup[i]-1,0,InGroup[i+1]-InGroup[i]+1,1).sum();
                        Rs2 += Rd.block(InGroup[i]-1,ij,InGroup[i+1]-InGroup[i]+1,1).sum();
                        Rs3 += Rdd.block(InGroup[i]-1,ij*totalnum+jk,InGroup[i+1]-InGroup[i]+1,1).sum();
                        at_risk += InGroup[i+1]-InGroup[i]+1;
                    }
                    //
                    MatrixXd Ld = MatrixXd::Zero(dj,3);
                    Ld << R.block(RiskFail(j,0),0,dj,1), Rd.block(RiskFail(j,0),ij,dj,1) ,Rdd.block(RiskFail(j,0),ij*totalnum+jk,dj,1);//sum of risks in group
                    //
                    MatrixXd Ldm = MatrixXd::Zero(dj,3);
                    for (int i = 0; i < dj; i++){ //adds in the efron approximation terms
                        Ldm.row(i) = (-double(i) / double(dj)) * Ld.colwise().sum().array();
                    }
                    Ldm.col(0) = Ldm.col(0).array() + Rs1;
                    Ldm.col(1) = Ldm.col(1).array() + Rs2;
                    Ldm.col(2) = Ldm.col(2).array() + Rs3;
                    //
                    MatrixXd temp1 = MatrixXd::Zero(Ld.rows(),1);
                    temp1 = Ld.col(0).array().log();
                    double Ld1 =  (temp1.array().isFinite()).select(temp1,0).sum();
                    temp1 = Ld.col(1).array() * (Ld.col(0).array().pow(-1));
                    double Ld2 = (temp1.array().isFinite()).select(temp1,0).sum();
                    temp1 = Ld.col(2).array() * (Ld.col(0).array().pow(-1)) - temp1.array().pow(2);
                    double Ld3 = (temp1.array().isFinite()).select(temp1,0).sum();
                    //
                    temp1 = Ldm.col(0).array().log();
                    Rs1 =  (temp1.array().isFinite()).select(temp1,0).sum();
                    temp1 = Ldm.col(1).array() * (Ldm.col(0).array().pow(-1));
                    Rs2 = (temp1.array().isFinite()).select(temp1,0).sum();
                    temp1 = Ldm.col(2).array() * (Ldm.col(0).array().pow(-1)) - temp1.array().pow(2);
                    Rs3 = (temp1.array().isFinite()).select(temp1,0).sum();
                    //
                    if (ij==jk){
                        Ll[ij] += Ld1 - Rs1;
                        Lld[ij] += Ld2 - Rs2;
                    }
                    Lldd[ij*totalnum+jk] += Ld3 - Rs3; //sums the log-likelihood and derivatives
                    Lldd[jk*totalnum+ij] += Ld3 - Rs3; //sums the log-likelihood and derivatives
                }
            }
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
        end_point = high_resolution_clock::now();
        end = time_point_cast<microseconds>(end_point).time_since_epoch().count();
        cout<<"df100 "<<(end-start)<<" "<<halves<<" "<<iteration<<" "<<ind0<<",Recalc"<<endl;
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
    NumericVector Lldd_vec = wrap(Lldd);
    Lldd_vec.attr("dim") = Dimension(totalnum, totalnum);
    //
    List res_list = List::create(_["LogLik"]=wrap(Ll),_["First_Der"]=wrap(Lld),_["Second_Der"]=Lldd_vec,_["beta_0"]=wrap(beta_0),_["AIC"]=2*totalnum-2*Ll[fir]);
    //
    return res_list;
}
