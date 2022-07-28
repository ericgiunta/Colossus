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
void Write_Ind_File(NumericMatrix df, NumericVector tu){
    //----------------------------------------------------------------------------------------------------------------//
    const Map<MatrixXd> df_m(as<Map<MatrixXd> >(df));
    vector <string> ind_list(tu.size(),"");
    cout.precision(10); //forces higher precision numbers printed to terminal
    int nthreads = Eigen::nbThreads(); //stores how many threads are allocated
    #pragma omp parallel for schedule(dynamic) num_threads(nthreads)
    for (int ijk=0;ijk<tu.size();ijk++){
        double t0 = tu[ijk];
        VectorXi select_ind_all = ((df_m.col(0).array() < t0)&&(df_m.col(1).array()>=t0)).cast<int>();
        vector<int> indices_all;
        VectorXi select_ind_end = ((df_m.col(2).array() == 1)&&(df_m.col(1).array()==t0)).cast<int>();
        vector<int> indices_end;
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
        indices.push_back(indices_end[0]);
        indices.push_back(indices_end[indices_end.size()-1]);
        //
        ostringstream oss;
        copy(indices.begin(), indices.end(),
            std::ostream_iterator<int>(oss, ","));
        ind_list[ijk] = oss.str();
    }
    ofstream file_out;
    file_out.open("test.txt");
    for (int ijk=0;ijk<ind_list.size();ijk++){
        file_out << ind_list[ijk]<<endl;
    }
    file_out.close();
    //----------------------------------------------------------------------------------------------------------------//
    return;
}

// [[Rcpp::export]]
List peanut_transition(NumericVector a_lin,NumericVector a_loglin,NumericVector a_plin, NumericMatrix x_lin, NumericMatrix x_loglin, NumericMatrix x_plin, NumericMatrix x_dose,int fir,string modelform,int ntime, NumericVector include_bool, List Control, List Dose_paras){
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
    int maxiter = Control["maxiter"];
	int halfmax = Control["halfmax"];
	double epsilon = Control["epsilon"];
	double dbeta_cap = Control["dbeta_max"];
	double abs_max = Control["abs_max"];
	double deriv_epsilon =Control["deriv_epsilon"];
	int batch_size =Control["batch_size"];
	List beta_loglin_slope = Dose_paras["beta_loglin_slope"];
    List beta_loglin_top  = Dose_paras["beta_loglin_top"];
    List beta_lin_slope  = Dose_paras["beta_lin_slope"];
    List beta_lin_int  = Dose_paras["beta_lin_int"];
    List beta_quad  = Dose_paras["beta_quad"];
    List beta_step_slope  = Dose_paras["beta_step_slope"];
    List beta_step_int  = Dose_paras["beta_step_int"];
    //
    //----------------------------------------------------------------------------------------------------------------//
    List res = LogLik_PEANUT(beta_lin,beta_loglin,beta_plin,df_lin,df_loglin,df_plin, df_dose,fir,modelform,ntime,include_bool, lr, maxiter, halfmax, epsilon, dbeta_cap, abs_max, deriv_epsilon, batch_size,beta_loglin_slope, beta_loglin_top , beta_lin_slope , beta_lin_int , beta_quad , beta_step_slope , beta_step_int );
    //----------------------------------------------------------------------------------------------------------------//
    return res;
}

//List peanut_plot(NumericVector a_lin,NumericVector a_loglin,NumericVector a_plin, NumericVector a_dose, NumericMatrix x_lin, NumericMatrix x_loglin, NumericMatrix x_plin, NumericMatrix x_dose,int fir,string modelform, string doseform, StringVector dose_terms,int ntime, NumericVector include_bool, int batch_size, NumericVector beta_loglin_slope, NumericVector beta_loglin_top , NumericVector beta_lin_slope , NumericVector beta_lin_int , NumericVector beta_quad , NumericVector beta_step_slope , NumericVector beta_step_int ){
//    //----------------------------------------------------------------------------------------------------------------//
//    Map<VectorXd> beta_lin(as<Map<VectorXd> >(a_lin));
//    Map<VectorXd> beta_loglin(as<Map<VectorXd> >(a_loglin));
//    Map<VectorXd> beta_plin(as<Map<VectorXd> >(a_plin));
//    Map<VectorXd> beta_dose(as<Map<VectorXd> >(a_dose));
//    const Map<MatrixXd> df_lin(as<Map<MatrixXd> >(x_lin));
//    const Map<MatrixXd> df_loglin(as<Map<MatrixXd> >(x_loglin));
//    const Map<MatrixXd> df_plin(as<Map<MatrixXd> >(x_plin));
//    const Map<MatrixXd> df_dose(as<Map<MatrixXd> >(x_dose));
//    // Converts from Rcpp types to efficient Eigen types
//    //----------------------------------------------------------------------------------------------------------------//
//    List res = PEANUT_PLOT(beta_lin,beta_loglin,beta_plin, beta_dose,df_lin,df_loglin,df_plin, df_dose,fir,modelform, doseform, dose_terms,ntime,include_bool, batch_size);
//    //----------------------------------------------------------------------------------------------------------------//
//    return res;
//}

List LogLik_PEANUT( VectorXd beta_linT,VectorXd beta_loglinT,VectorXd beta_plinT,MatrixXd df_lin,MatrixXd df_loglin,MatrixXd df_plin, MatrixXd df_dose,int fir,string modelform,int ntime, NumericVector include_bool, double lr, int maxiter, int halfmax, double epsilon, double dbeta_cap, double abs_max, double deriv_epsilon, int batch_size,List beta_loglin_slopes, List beta_loglin_tops , List beta_lin_slopes , List beta_lin_ints , List beta_quads , List beta_step_slopes , List beta_step_ints ){
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
    vector<double> cumulative_dose_num(df_dose.cols(),0);
    int dose_num_tot=0;
    int loglin_size=0;
    int lin_size=0;
    int quad_size=0;
    int step_size=0;
    for (int ijk=0;ijk<df_dose.cols();ijk++){
        cumulative_dose_num[ijk] = dose_num_tot;
        NumericVector beta_loglin_slope = beta_loglin_slopes[ijk];
        NumericVector beta_lin_slope = beta_lin_slopes[ijk];
        NumericVector beta_quad = beta_quads[ijk];
        NumericVector beta_step_slope = beta_step_slopes[ijk];
        //
        if ((beta_loglin_slope.size()==1)&&(beta_loglin_slope[0]==0.0){
            ;
        } else {
            dose_num_tot += beta_loglin_slope.size();
            loglin_size = beta_loglin_slope.size();
        }
        if ((beta_lin_slope.size()==1)&&(beta_lin_slope[0]==0.0){
            ;
        } else {
            dose_num_tot += beta_lin_slope.size();
            lin_size = beta_lin_slope.size();
        }
        if ((beta_quad.size()==1)&&(beta_quad[0]==0.0){
            ;
        } else {
            dose_num_tot += beta_quad.size();
            quad_size = beta_quad.size();
        }
        if ((beta_step_slope.size()==1)&&(beta_step_slope[0]==0.0){
            ;
        } else {
            dose_num_tot += beta_step_slope.size();
            step_size = beta_step_int.size();
        }
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
    end_point = system_clock::now();
    end = time_point_cast<microseconds>(end_point).time_since_epoch().count();
    cout<<"df99,"<<(end-start)<<",Batches"<<endl;
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
    #pragma omp parallel for schedule(dynamic) num_threads(nthreads)
    for (int ij=0;ij<totalnum;ij++){
        int ind0 = ij;
        if (ind0 < dose_num_tot){
            ;
        } else {
            ind0 = ind0 - dose_num_tot;
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
            } else if (tform[ij]=="plin") {
                T0.col(ij) = 1 + T0.col(ij).array();
                Td0.col(ij) = df0.col(ij);
            } else {
                cout << tform[ij] << " is invalid" << endl;
                throw invalid_argument( "Invalid term type" );
            }
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
    NumericVector beta_loglin_slope;
    NumericVector beta_loglin_top;
    NumericVector beta_lin_slope;
    NumericVector beta_lin_int;
    NumericVector beta_quad;
    NumericVector beta_step_slope;
    NumericVector beta_step_int;
    //
    ArrayXd De(df_dose.rows(),dose_num_tot);
    ArrayXd Dde(df_dose.rows(),dose_num_tot);
    ArrayXd Ddde(df_dose.rows(),dose_num_tot*dose_num_tot);
    double dint = 0.1;
    int total_dose=0;
    cout << 1 << endl;
    for (int ijk=0;ijk<df_dose.cols();ijk++){
        beta_loglin_slope = beta_loglin_slopes[ijk];
        beta_loglin_top = beta_loglin_tops[ijk];
        beta_lin_slope = beta_lin_slopes[ijk];
        beta_lin_int = beta_lin_ints[ijk];
        beta_quad = beta_quads[ijk];
        beta_step_slope = beta_step_slopes[ijk];
        beta_step_int = beta_step_ints[ijk];
        //
        total_dose = loglin_size + lin_size + quad_size + step_size;//beta_loglin_slope.size() + beta_lin_slope.size() + beta_quad.size() + beta_step_int.size();
        cout << 1 << endl;
//        #pragma omp parallel for schedule(dynamic) num_threads(nthreads)
        int ind0;
        int ind1;
        int jk;
        ArrayXd temp;
        ArrayXd temp0;
        ArrayXd temp1;
        ArrayXd temp2;
        for (int ij=0; ij<total_dose;ij++){
            cout << ij << endl;
            if (ij < loglin_size){
                temp = (beta_loglin_top[ij] * df_dose.col(ijk)).array().exp();
                temp1 = beta_loglin_slope[ij] * temp;
                //
                ind0 = cumulative_dose_num[ijk]+2*ij;
                ind1 = ind0 + 1; 
                //
                beta_0[ind0] = beta_loglin_slope[ij];
                beta_0[ind1] = beta_loglin_top[ij];
                tform[ind0] = "loglin_slope";
                tform[ind1] = "loglin_top";
                df0.col(ind0) = df_dose.col(ijk);
                df0.col(ind1) = df_dose.col(ijk);
                //
                De.col(ind0) = temp1;
                De.col(ind1) = temp1;
                Dde.col(ind0) = temp.array();
                Dde.col(ind1) = temp1.array() * df_dose.col(ijk).array();
                Ddde.col(ind0 * dose_num_tot + ind1) = temp.array() * df_dose.col(ijk).array();
                Ddde.col(ind1 * dose_num_tot + ind1) = temp1.array() * df_dose.col(ijk).array().pow(2).array();
                
            } else if (ij < loglin_size + lin_size){
                jk = ij - loglin_size;
                temp = (df_dose.col(ijk).array() - beta_lin_int[jk]);
                temp0 = (df_dose.col(ijk).array() - beta_lin_int[jk]-dint);
                temp1 = (df_dose.col(ijk).array() - beta_lin_int[jk]+dint);
                temp2 = (temp1.array() < 0).select(0, temp1) - (temp0.array() < 0).select(0, temp0);
                //
                ind0 = cumulative_dose_num[ijk]+2*loglin_size + 2*jk;
                ind1 = ind0 + 1; 
                //
                beta_0[ind0] = beta_lin_slope[ij];
                beta_0[ind1] = beta_lin_int[ij];
                tform[ind0] = "lin_slope";
                tform[ind1] = "lin_ind";
                df0.col(ind0) = df_dose.col(ijk);
                df0.col(ind1) = df_dose.col(ijk);
                //
                De.col(ind0) = beta_lin_slope[jk] * (temp.array() < 0).select(0.0, temp);
                De.col(ind1) = beta_lin_slope[jk] * (temp.array() < 0).select(0.0, temp);
                Dde.col(ind0) = (temp.array() < 0).select(0.0, temp);
                Dde.col(ind1) = beta_lin_slope[jk] * (temp2) / 2/dint;
                //
                Ddde.col(ind0 * dose_num_tot + ind1) = (temp2) / 2/dint;
                Ddde.col(ind1 * dose_num_tot + ind1) = beta_lin_slope[jk] * (temp2) / pow(dint,2);
            } else if (ij < loglin_size + lin_size + quad_size){
                jk = ij - loglin_size - lin_size;
                temp = df_dose.col(ijk).array().pow(2);
                ind0 = cumulative_dose_num[ijk]+2*loglin_size + 2*lin_size+jk;
                //
                beta_0[ind0] = beta_quad[ij];
                tform[ind0] = "quad_slope";
                df0.col(ind0) = df_dose.col(ijk);
                //
                De.col(ind0) = beta_quad[jk] * temp.array();
                Dde.col(ind0) = temp.array();
            } else {
                jk = ij - loglin_size - lin_size - quad_size;
                temp = (df_dose.col(ijk).array() - beta_step_int[jk]);
                temp0 = (df_dose.col(ijk).array() - beta_step_int[jk]-dint);
                temp1 = (df_dose.col(ijk).array() - beta_step_int[jk]+dint);
                temp2 = (temp1.array() < 0).select(0.0, MatrixXd::Zero(temp.rows(),temp.cols()).array()+1.0) - (temp0.array() < 0).select(0.0, MatrixXd::Zero(temp.rows(),temp.cols()).array()+1.0);
                //
                ind0 = cumulative_dose_num[ijk]+2*loglin_size + 2*lin_size + quad_size + 2*jk;
                ind1 = ind0 + 1;
                //
                beta_0[ind0] = beta_step_slope[ij];
                beta_0[ind1] = beta_step_int[ij];
                tform[ind0] = "step_slope";
                tform[ind1] = "step_ind";
                df0.col(ind0) = df_dose.col(ijk);
                df0.col(ind1) = df_dose.col(ijk);
                //
                De.col(ind0) = beta_step_slope[jk] * (temp.array() < 0).select(0.0, MatrixXd::Zero(temp.rows(),temp.cols()).array()+1.0);
                De.col(ind1) = beta_step_slope[jk] * (temp.array() < 0).select(0.0, MatrixXd::Zero(temp.rows(),temp.cols()).array()+1.0);
                Dde.col(ind0) = (temp.array() < 0).select(0.0, MatrixXd::Zero(temp.rows(),temp.cols()).array()+1.0);
                Dde.col(ind1) = beta_step_slope[jk] * (temp2) / 2/dint;
                //
                Ddde.col(ind0 * dose_num_tot + ind1) = (temp2) / 2/dint;
                Ddde.col(ind1 * dose_num_tot + ind1) = beta_step_slope[jk] * (temp2) / pow(dint,2);
            }
        }
    }
    cout << 1 << endl;
    ArrayXd Dose(df_dose.rows(),1);
    Dose <<  De.rowwise().sum();
    //
    if ((modelform=="A")||(modelform=="PA")||(modelform=="PAE")){ //same process used for all of the additive type models
        Te = T0.array().block(0,dose_num_tot,T0.rows(),T0.cols()-dose_num_tot).rowwise().sum() + Dose.array();
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
                    Rdd.col(ij*totalnum+ij) = Ddde.col(ij*dose_num_tot+jk);
                } else if (ij==jk) {
                    Rdd.col(ij*totalnum+ij) = Tdd0.col(ij);
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
                        if (fir!=0){
                            if (ij==fir){
                                Rdd.block(batch_cols[bl],ij*totalnum+jk,batch_cols[bl+1]-batch_cols[bl],1) = Tdd0.block(batch_cols[bl],ij,batch_cols[bl+1]-batch_cols[bl],1).array() * Te.block(batch_cols[bl],0,batch_cols[bl+1]-batch_cols[bl],1).array();
                            } else {
                                if (ij<dose_num_tot){
                                    Rdd.block(batch_cols[bl],ij*totalnum+jk,batch_cols[bl+1]-batch_cols[bl],1) = Ddde.block(batch_cols[bl],ij*dose_num_tot+jk,batch_cols[bl+1]-batch_cols[bl],1).array() * T0.block(batch_cols[bl],fir,batch_cols[bl+1]-batch_cols[bl],1).array();
                                } else {
                                    Rdd.block(batch_cols[bl],ij*totalnum+jk,batch_cols[bl+1]-batch_cols[bl],1) = Tdd0.block(batch_cols[bl],ij,batch_cols[bl+1]-batch_cols[bl],1).array() * T0.block(batch_cols[bl],fir,batch_cols[bl+1]-batch_cols[bl],1).array();
                                }
                            }
                        } else {
                            if (ij<dose_num_tot){
                                Rdd.block(batch_cols[bl],ij*totalnum+jk,batch_cols[bl+1]-batch_cols[bl],1) = Ddde.block(batch_cols[bl],ij*dose_num_tot+jk,batch_cols[bl+1]-batch_cols[bl],1).array() * Te.block(batch_cols[bl],0,batch_cols[bl+1]-batch_cols[bl],1).array();
                            } else {
                                Rdd.block(batch_cols[bl],ij*totalnum+jk,batch_cols[bl+1]-batch_cols[bl],1) = Tdd0.block(batch_cols[bl],ij,batch_cols[bl+1]-batch_cols[bl],1).array() * De.block(batch_cols[bl],0,batch_cols[bl+1]-batch_cols[bl],1).array();
                            }
                        }
                    } else {
                        if (fir!=0){
                            if ((ij==fir)||(jk==fir)){
                                if (ij<dose_num_tot){
                                    Rdd.block(batch_cols[bl],ij*totalnum+jk,batch_cols[bl+1]-batch_cols[bl],1) = Dde.block(batch_cols[bl],ij,batch_cols[bl+1]-batch_cols[bl],1).array() * Td0.block(batch_cols[bl],jk,batch_cols[bl+1]-batch_cols[bl],1).array();
                                    Rdd.block(batch_cols[bl],jk*totalnum+ij,batch_cols[bl+1]-batch_cols[bl],1) = Dde.block(batch_cols[bl],ij,batch_cols[bl+1]-batch_cols[bl],1).array() * Td0.block(batch_cols[bl],jk,batch_cols[bl+1]-batch_cols[bl],1).array();
                                } else {
                                    Rdd.block(batch_cols[bl],ij*totalnum+jk,batch_cols[bl+1]-batch_cols[bl],1) = Td0.block(batch_cols[bl],ij,batch_cols[bl+1]-batch_cols[bl],1).array() * Td0.block(batch_cols[bl],jk,batch_cols[bl+1]-batch_cols[bl],1).array();
                                    Rdd.block(batch_cols[bl],jk*totalnum+ij,batch_cols[bl+1]-batch_cols[bl],1) = Td0.block(batch_cols[bl],ij,batch_cols[bl+1]-batch_cols[bl],1).array() * Td0.block(batch_cols[bl],jk,batch_cols[bl+1]-batch_cols[bl],1).array();
                                }
                            } else if ((ij<dose_num_tot)&&(jk<dose_num_tot)){
                                Rdd.block(batch_cols[bl],ij*totalnum+jk,batch_cols[bl+1]-batch_cols[bl],1) = Ddde.block(batch_cols[bl],ij*dose_num_tot+jk,batch_cols[bl+1]-batch_cols[bl],1).array() * Dde.block(batch_cols[bl],jk,batch_cols[bl+1]-batch_cols[bl],1).array();
                                Rdd.block(batch_cols[bl],jk*totalnum+ij,batch_cols[bl+1]-batch_cols[bl],1) = Ddde.block(batch_cols[bl],ij*dose_num_tot+jk,batch_cols[bl+1]-batch_cols[bl],1).array() * Dde.block(batch_cols[bl],jk,batch_cols[bl+1]-batch_cols[bl],1).array();
                            }
                        }
                    }
                }
            }
        }
    }else if (modelform=="M"){
        Te = Te.array() * 0 + 1; //verifies the initial term product is 1
        //
        cout << 1 << endl;
        Te = T0.array().block(0,dose_num_tot,T0.rows(),T0.cols()-dose_num_tot).rowwise().prod() * Dose.array();
        // computes intial risk and derivatives
        R << Te.array();
        Rd = T0.array().pow(-1).array() * Te.colwise().replicate(totalnum).array();
        Rd = Rd.array() * Td0.array();
        #pragma omp parallel for schedule(dynamic) num_threads(nthreads)
        for (int ijk=0;ijk<totalnum;ijk++){
            if (ijk<dose_num_tot){
                Rd.col(ijk) = De.col(ijk).array().pow(-1).array()* Te.array();
                Rd.col(ijk) = Rd.col(ijk).array() * Dde.col(ijk).array();
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
                if (ij==jk){
                    if (ij<dose_num_tot){
                        Rdd.block(batch_cols[bl],ij*totalnum+jk,batch_cols[bl+1]-batch_cols[bl],1) = Ddde.block(batch_cols[bl],ij*dose_num_tot+jk,batch_cols[bl+1]-batch_cols[bl],1).array() * De.block(batch_cols[bl],jk,batch_cols[bl+1]-batch_cols[bl],1).array().pow(-1).array() * R.block(batch_cols[bl],1,batch_cols[bl+1]-batch_cols[bl],1).array();
                    } else {
                        Rdd.block(batch_cols[bl],ij*totalnum+jk,batch_cols[bl+1]-batch_cols[bl],1) = Tdd0.block(batch_cols[bl],jk,batch_cols[bl+1]-batch_cols[bl],1).array() * Td0.block(batch_cols[bl],jk,batch_cols[bl+1]-batch_cols[bl],1).array().pow(-1).array() * Rd.block(batch_cols[bl],ij,batch_cols[bl+1]-batch_cols[bl],1).array();
                    }
                } else {
                    if ((ij<dose_num_tot)&&(jk<dose_num_tot)){
                        Rdd.block(batch_cols[bl],ij*totalnum+jk,batch_cols[bl+1]-batch_cols[bl],1) = Ddde.block(batch_cols[bl],ij*dose_num_tot+jk,batch_cols[bl+1]-batch_cols[bl],1).array() * De.block(batch_cols[bl],ij,batch_cols[bl+1]-batch_cols[bl],1).array().pow(-1).array() * R.block(batch_cols[bl],1,batch_cols[bl+1]-batch_cols[bl],1).array();
                        Rdd.block(batch_cols[bl],jk*totalnum+ij,batch_cols[bl+1]-batch_cols[bl],1) = Ddde.block(batch_cols[bl],ij*dose_num_tot+jk,batch_cols[bl+1]-batch_cols[bl],1).array() * De.block(batch_cols[bl],ij,batch_cols[bl+1]-batch_cols[bl],1).array().pow(-1).array() * R.block(batch_cols[bl],1,batch_cols[bl+1]-batch_cols[bl],1).array();
                    } else if ((ij<dose_num_tot)||(jk<dose_num_tot)){
                        Rdd.block(batch_cols[bl],ij*totalnum+jk,batch_cols[bl+1]-batch_cols[bl],1) = Dde.block(batch_cols[bl],jk,batch_cols[bl+1]-batch_cols[bl],1).array() * De.block(batch_cols[bl],jk,batch_cols[bl+1]-batch_cols[bl],1).array().pow(-1).array() * Rd.block(batch_cols[bl],ij,batch_cols[bl+1]-batch_cols[bl],1).array();
                        Rdd.block(batch_cols[bl],jk*totalnum+ij,batch_cols[bl+1]-batch_cols[bl],1) = Dde.block(batch_cols[bl],jk,batch_cols[bl+1]-batch_cols[bl],1).array() * De.block(batch_cols[bl],jk,batch_cols[bl+1]-batch_cols[bl],1).array().pow(-1).array() * Rd.block(batch_cols[bl],ij,batch_cols[bl+1]-batch_cols[bl],1).array();
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
    cout << 1 << endl;
    R = (R.array().isFinite()).select(R,0);
    Rd = (Rd.array().isFinite()).select(Rd,0);
    Rdd = (Rdd.array().isFinite()).select(Rdd,0);
    //
    // -------------------------------------------------------------------------------------------
    string line;
    ifstream infile("test.txt"); //The file of risk group rows
    vector<string>  RiskGroup(ntime); //vector of strings detailing the rows
    IntegerMatrix RiskFail(ntime,2); //vector giving the event rows
    int j_iter = 0;
    cout << 1 << endl;
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
    #pragma omp parallel for schedule(dynamic) num_threads(nthreads) reduction(vec_double_plus:Ll,Lld,Lldd) collapse(2)
    for (int ijk=0;ijk<totalnum;ijk++){//totalnum*(totalnum+1)/2
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
            for (int i = 0; i < InGroup.size()-1; i=i+2){
                Rs1 += R.block(InGroup[i]-1,0,InGroup[i+1]-InGroup[i]+1,1).sum();
                Rs2 += Rd.block(InGroup[i]-1,ij,InGroup[i+1]-InGroup[i]+1,1).sum();
                Rs2t += Rd.block(InGroup[i]-1,jk,InGroup[i+1]-InGroup[i]+1,1).sum();
                Rs3 += Rdd.block(InGroup[i]-1,ij*totalnum+jk,InGroup[i+1]-InGroup[i]+1,1).sum();
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
        }
    }
    #pragma omp parallel for num_threads(nthreads)
    for (int ijk=0;ijk<totalnum*(totalnum+1)/2;ijk++){//totalnum*(totalnum+1)/2
        int ij = 0;
        int jk = ijk;
        while (jk>ij){
            ij++;
            jk-=ij;
        }
        Lldd[jk*totalnum+ij] = Lldd[ij*totalnum+jk];
    }
    //
    vector <double> Ll_comp(2,Ll[0]);
    //
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
            if (abs(dbeta[ijk])>abs_max){
                dbeta[ijk] = abs_max * sign(dbeta[ijk]);
            }
        }
        //
        VectorXd::Map(&beta_p[0], beta_0.size()) = beta_0;//wrap(beta_0);
        VectorXd::Map(&beta_c[0], beta_0.size()) = beta_0;//beta_c = wrap(beta_0);
        VectorXd::Map(&beta_a[0], beta_0.size()) = beta_0;//beta_a = wrap(beta_0);
        VectorXd::Map(&beta_best[0], beta_0.size()) = beta_0;//beta_best = wrap(beta_0);
        Ll_best = Ll[ind0];
        i = ind0;
        //
        halves=0;
        while ((Ll[ind0] <= Ll_best)&&(halves<halfmax)&&(abs(dbeta[ind0]) > epsilon)){
            halves++;
            #pragma omp parallel for schedule(dynamic) num_threads(nthreads)
            for (int ijk=0;ijk<totalnum;ijk++){
                //
                if (tform[ijk]=="lin"){
                    beta_c[ijk] = beta_a[ijk] + dbeta[ijk];
                    T0.col(ijk) = T0.col(ijk).array() * (beta_c[ijk] / beta_p[ijk]);
                } else if (tform[ijk]=="plin"){
                    beta_c[ijk] = beta_a[ijk] + dbeta[ijk];
                    T0.col(ijk) = T0.col(ijk).array() * (1 + beta_c[ijk] * df0.col(ijk).array()) / (1 + beta_p[ijk] * df0.col(ijk).array());
                } else if (tform[ijk]=="loglin") {
                    beta_c[ijk] = beta_a[ijk] + dbeta[ijk];
                    T0.col(ijk) = T0.col(ijk).array() * ((beta_c[ijk] - beta_p[ijk]) * df0.col(ijk)).array().exp().array();
                    Td0.col(ijk) = Td0.col(ijk).array() * ((beta_c[ijk] - beta_p[ijk]) * df0.col(ijk)).array().exp().array();
                    Tdd0.col(ijk) = Tdd0.col(ijk).array() * ((beta_c[ijk] - beta_p[ijk]) * df0.col(ijk)).array().exp().array();
                } else if (tform[ijk]=="loglin_slope"){
                    beta_c[ijk] = beta_a[ijk] + dbeta[ijk];
                    beta_c[ijk+1] = beta_a[ijk+1] + dbeta[ijk+1];
                    double ach = beta_c[ijk]/beta_p[ijk];
                    MatrixXd bch = ((dbeta[ijk+1]) * df0.col(ijk)).array().exp().array();
                    //
                    De.col(ijk) = ach * bch.array() * De.col(ijk).array();
                    Dde.col(ijk) = bch.array() * Dde.col(ijk).array();
                    Dde.col(ijk+1) = ach * bch.array() * Dde.col(ijk+1).array();
                    Ddde.col(ijk*dose_num_tot+ijk+1) = bch.array() * Ddde.col(ijk*dose_num_tot+ijk+1).array();
                    Ddde.col(ijk*dose_num_tot+ijk+1+dose_num_tot) = ach * bch.array() * Ddde.col(ijk*dose_num_tot+ijk+1+dose_num_tot).array();
                    //
                } else if (tform[ijk]=="lin_slope"){
                    beta_c[ijk] = beta_a[ijk] + dbeta[ijk];
                    beta_c[ijk+1] = beta_a[ijk+1] + dbeta[ijk+1];
                    //
                    ArrayXd temp = (df0.col(ijk).array() - beta_c[ijk+1]);
                    ArrayXd temp0 = (df0.col(ijk).array() - beta_c[ijk+1]-dint);
                    ArrayXd temp1 = (df0.col(ijk).array() - beta_c[ijk+1]+dint);
                    ArrayXd temp2 = (temp1.array() < 0).select(0, temp1) - (temp0.array() < 0).select(0, temp0);
                    //
                    De.col(ijk) = beta_c[ijk] * (temp.array() < 0).select(0.0, temp);
                    De.col((ijk+1)) = beta_c[ijk] * (temp.array() < 0).select(0.0, temp);
                    Dde.col(ijk) = (temp.array() < 0).select(0.0, temp);
                    Dde.col((ijk+1)) = beta_c[ijk] * (temp2) / 2/dint;
                    //
                    Ddde.col(ijk * dose_num_tot + (ijk+1)) = (temp2) / 2/dint;
                    Ddde.col((ijk+1) * dose_num_tot + (ijk+1)) = beta_c[ijk] * (temp2) / pow(dint,2);
                    //
                } else if (tform[ijk]=="quad_slope"){
                    beta_c[ijk] = beta_a[ijk] + dbeta[ijk];
                    De.col(ijk) = beta_c[ijk] * df0.col(ijk);
                    //
                } else if (tform[ijk]=="step_slope"){
                    beta_c[ijk] = beta_a[ijk] + dbeta[ijk];
                    beta_c[ijk+1] = beta_a[ijk+1] + dbeta[ijk+1];
                    
                    ArrayXd temp = (df0.col(ijk).array() - beta_0[ijk+1]);
                    ArrayXd temp0 = (df0.col(ijk).array() - beta_0[ijk+1]-dint);
                    ArrayXd temp1 = (df0.col(ijk).array() - beta_0[ijk+1]+dint);
                    ArrayXd temp2 = (temp1.array() < 0).select(0.0, MatrixXd::Zero(temp.rows(),temp.cols()).array()+1.0) - (temp0.array() < 0).select(0.0, MatrixXd::Zero(temp.rows(),temp.cols()).array()+1.0);
                    //
                    De.col(ijk) = beta_0[ijk] * (temp.array() < 0).select(0.0, MatrixXd::Zero(temp.rows(),temp.cols()).array()+1.0);
                    De.col(ijk+1) = beta_0[ijk] * (temp.array() < 0).select(0.0, MatrixXd::Zero(temp.rows(),temp.cols()).array()+1.0);
                    Dde.col(ijk) = (temp.array() < 0).select(0.0, MatrixXd::Zero(temp.rows(),temp.cols()).array()+1.0);
                    Dde.col(ijk+1) = beta_0[ijk] * (temp2) / 2/dint;
                    //
                    Ddde.col(ijk * dose_num_tot + ijk+1) = (temp2) / 2/dint;
                    Ddde.col((ijk+1) * dose_num_tot + ijk+1) = beta_0[ijk] * (temp2) / pow(dint,2);
                    //
                } else {
                    ;
                }
                //
            }
            //
            Dose << De.rowwise().sum();
            //
            if ((modelform=="A")||(modelform=="PA")||(modelform=="PAE")){ //same process used for all of the additive type models
                Te = T0.array().block(0,dose_num_tot,T0.rows(),T0.cols()-dose_num_tot).rowwise().sum() + Dose.array();
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
                            Rdd.col(ij*totalnum+ij) = Ddde.col(ij*dose_num_tot+jk);
                        } else if (ij==jk) {
                            Rdd.col(ij*totalnum+ij) = Tdd0.col(ij);
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
                                if (fir!=0){
                                    if (ij==fir){
                                        Rdd.block(batch_cols[bl],ij*totalnum+jk,batch_cols[bl+1]-batch_cols[bl],1) = Tdd0.block(batch_cols[bl],ij,batch_cols[bl+1]-batch_cols[bl],1).array() * Te.block(batch_cols[bl],0,batch_cols[bl+1]-batch_cols[bl],1).array();
                                    } else {
                                        if (ij<dose_num_tot){
                                            Rdd.block(batch_cols[bl],ij*totalnum+jk,batch_cols[bl+1]-batch_cols[bl],1) = Ddde.block(batch_cols[bl],ij*dose_num_tot+jk,batch_cols[bl+1]-batch_cols[bl],1).array() * T0.block(batch_cols[bl],fir,batch_cols[bl+1]-batch_cols[bl],1).array();
                                        } else {
                                            Rdd.block(batch_cols[bl],ij*totalnum+jk,batch_cols[bl+1]-batch_cols[bl],1) = Tdd0.block(batch_cols[bl],ij,batch_cols[bl+1]-batch_cols[bl],1).array() * T0.block(batch_cols[bl],fir,batch_cols[bl+1]-batch_cols[bl],1).array();
                                        }
                                    }
                                } else {
                                    if (ij<dose_num_tot){
                                        Rdd.block(batch_cols[bl],ij*totalnum+jk,batch_cols[bl+1]-batch_cols[bl],1) = Ddde.block(batch_cols[bl],ij*dose_num_tot+jk,batch_cols[bl+1]-batch_cols[bl],1).array() * Te.block(batch_cols[bl],0,batch_cols[bl+1]-batch_cols[bl],1).array();
                                    } else {
                                        Rdd.block(batch_cols[bl],ij*totalnum+jk,batch_cols[bl+1]-batch_cols[bl],1) = Tdd0.block(batch_cols[bl],ij,batch_cols[bl+1]-batch_cols[bl],1).array() * De.block(batch_cols[bl],0,batch_cols[bl+1]-batch_cols[bl],1).array();
                                    }
                                }
                            } else {
                                if (fir!=0){
                                    if ((ij==fir)||(jk==fir)){
                                        if (ij<dose_num_tot){
                                            Rdd.block(batch_cols[bl],ij*totalnum+jk,batch_cols[bl+1]-batch_cols[bl],1) = Dde.block(batch_cols[bl],ij,batch_cols[bl+1]-batch_cols[bl],1).array() * Td0.block(batch_cols[bl],jk,batch_cols[bl+1]-batch_cols[bl],1).array();
                                            Rdd.block(batch_cols[bl],jk*totalnum+ij,batch_cols[bl+1]-batch_cols[bl],1) = Dde.block(batch_cols[bl],ij,batch_cols[bl+1]-batch_cols[bl],1).array() * Td0.block(batch_cols[bl],jk,batch_cols[bl+1]-batch_cols[bl],1).array();
                                        } else {
                                            Rdd.block(batch_cols[bl],ij*totalnum+jk,batch_cols[bl+1]-batch_cols[bl],1) = Td0.block(batch_cols[bl],ij,batch_cols[bl+1]-batch_cols[bl],1).array() * Td0.block(batch_cols[bl],jk,batch_cols[bl+1]-batch_cols[bl],1).array();
                                            Rdd.block(batch_cols[bl],jk*totalnum+ij,batch_cols[bl+1]-batch_cols[bl],1) = Td0.block(batch_cols[bl],ij,batch_cols[bl+1]-batch_cols[bl],1).array() * Td0.block(batch_cols[bl],jk,batch_cols[bl+1]-batch_cols[bl],1).array();
                                        }
                                    } else if ((ij<dose_num_tot)&&(jk<dose_num_tot)){
                                        Rdd.block(batch_cols[bl],ij*totalnum+jk,batch_cols[bl+1]-batch_cols[bl],1) = Ddde.block(batch_cols[bl],ij*dose_num_tot+jk,batch_cols[bl+1]-batch_cols[bl],1).array() * Dde.block(batch_cols[bl],jk,batch_cols[bl+1]-batch_cols[bl],1).array();
                                        Rdd.block(batch_cols[bl],jk*totalnum+ij,batch_cols[bl+1]-batch_cols[bl],1) = Ddde.block(batch_cols[bl],ij*dose_num_tot+jk,batch_cols[bl+1]-batch_cols[bl],1).array() * Dde.block(batch_cols[bl],jk,batch_cols[bl+1]-batch_cols[bl],1).array();
                                    }
                                }
                            }
                        }
                    }
                }
            }else if (modelform=="M"){
                Te = Te.array() * 0 + 1; //verifies the initial term product is 1
                //
                Te = T0.array().block(0,dose_num_tot,T0.rows(),T0.cols()-dose_num_tot).rowwise().prod() * Dose.array();
                // computes intial risk and derivatives
                R << Te.array();
                Rd = T0.array().pow(-1).array() * Te.colwise().replicate(totalnum).array();
                Rd = Rd.array() * Td0.array();
                #pragma omp parallel for schedule(dynamic) num_threads(nthreads)
                for (int ijk=0;ijk<totalnum;ijk++){
                    if (ijk<dose_num_tot){
                        Rd.col(ijk) = De.col(ijk).array().pow(-1).array()* Te.array();
                        Rd.col(ijk) = Rd.col(ijk).array() * Dde.col(ijk).array();
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
                        if (ij==jk){
                            if (ij<dose_num_tot){
                                Rdd.block(batch_cols[bl],ij*totalnum+jk,batch_cols[bl+1]-batch_cols[bl],1) = Ddde.block(batch_cols[bl],ij*dose_num_tot+jk,batch_cols[bl+1]-batch_cols[bl],1).array() * De.block(batch_cols[bl],jk,batch_cols[bl+1]-batch_cols[bl],1).array().pow(-1).array() * R.block(batch_cols[bl],1,batch_cols[bl+1]-batch_cols[bl],1).array();
                            } else {
                                Rdd.block(batch_cols[bl],ij*totalnum+jk,batch_cols[bl+1]-batch_cols[bl],1) = Tdd0.block(batch_cols[bl],jk,batch_cols[bl+1]-batch_cols[bl],1).array() * Td0.block(batch_cols[bl],jk,batch_cols[bl+1]-batch_cols[bl],1).array().pow(-1).array() * Rd.block(batch_cols[bl],ij,batch_cols[bl+1]-batch_cols[bl],1).array();
                            }
                        } else {
                            if ((ij<dose_num_tot)&&(jk<dose_num_tot)){
                                Rdd.block(batch_cols[bl],ij*totalnum+jk,batch_cols[bl+1]-batch_cols[bl],1) = Ddde.block(batch_cols[bl],ij*dose_num_tot+jk,batch_cols[bl+1]-batch_cols[bl],1).array() * De.block(batch_cols[bl],ij,batch_cols[bl+1]-batch_cols[bl],1).array().pow(-1).array() * R.block(batch_cols[bl],1,batch_cols[bl+1]-batch_cols[bl],1).array();
                                Rdd.block(batch_cols[bl],jk*totalnum+ij,batch_cols[bl+1]-batch_cols[bl],1) = Ddde.block(batch_cols[bl],ij*dose_num_tot+jk,batch_cols[bl+1]-batch_cols[bl],1).array() * De.block(batch_cols[bl],ij,batch_cols[bl+1]-batch_cols[bl],1).array().pow(-1).array() * R.block(batch_cols[bl],1,batch_cols[bl+1]-batch_cols[bl],1).array();
                            } else if ((ij<dose_num_tot)||(jk<dose_num_tot)){
                                Rdd.block(batch_cols[bl],ij*totalnum+jk,batch_cols[bl+1]-batch_cols[bl],1) = Dde.block(batch_cols[bl],jk,batch_cols[bl+1]-batch_cols[bl],1).array() * De.block(batch_cols[bl],jk,batch_cols[bl+1]-batch_cols[bl],1).array().pow(-1).array() * Rd.block(batch_cols[bl],ij,batch_cols[bl+1]-batch_cols[bl],1).array();
                                Rdd.block(batch_cols[bl],jk*totalnum+ij,batch_cols[bl+1]-batch_cols[bl],1) = Dde.block(batch_cols[bl],jk,batch_cols[bl+1]-batch_cols[bl],1).array() * De.block(batch_cols[bl],jk,batch_cols[bl+1]-batch_cols[bl],1).array().pow(-1).array() * Rd.block(batch_cols[bl],ij,batch_cols[bl+1]-batch_cols[bl],1).array();
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
                        Rs3 += Rdd.block(InGroup[i]-1,ij*totalnum+jk,InGroup[i+1]-InGroup[i]+1,1).sum();
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
                }
            }
            #pragma omp parallel for num_threads(nthreads)
            for (int ijk=0;ijk<totalnum*(totalnum+1)/2;ijk++){//totalnum*(totalnum+1)/2
                int ij = 0;
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
        }
        if (beta_best!=beta_p){
            #pragma omp parallel for schedule(dynamic) num_threads(nthreads)
            for (int ijk=0;ijk<totalnum;ijk++){
                //
                if (tform[ijk]=="lin"){
                    beta_c[ijk] = beta_a[ijk] + dbeta[ijk];
                    T0.col(ijk) = T0.col(ijk).array() * (beta_c[ijk] / beta_p[ijk]);
                } else if (tform[ijk]=="plin"){
                    beta_c[ijk] = beta_a[ijk] + dbeta[ijk];
                    T0.col(ijk) = T0.col(ijk).array() * (1 + beta_c[ijk] * df0.col(ijk).array()) / (1 + beta_p[ijk] * df0.col(ijk).array());
                } else if (tform[ijk]=="loglin") {
                    beta_c[ijk] = beta_a[ijk] + dbeta[ijk];
                    T0.col(ijk) = T0.col(ijk).array() * ((beta_c[ijk] - beta_p[ijk]) * df0.col(ijk)).array().exp().array();
                    Td0.col(ijk) = Td0.col(ijk).array() * ((beta_c[ijk] - beta_p[ijk]) * df0.col(ijk)).array().exp().array();
                    Tdd0.col(ijk) = Tdd0.col(ijk).array() * ((beta_c[ijk] - beta_p[ijk]) * df0.col(ijk)).array().exp().array();
                } else if (tform[ijk]=="loglin_slope"){
                    beta_c[ijk] = beta_a[ijk] + dbeta[ijk];
                    beta_c[ijk+1] = beta_a[ijk+1] + dbeta[ijk+1];
                    double ach = beta_c[ijk]/beta_p[ijk];
                    MatrixXd bch = ((dbeta[ijk+1]) * df0.col(ijk)).array().exp().array();
                    //
                    De.col(ijk) = ach * bch.array() * De.col(ijk).array();
                    Dde.col(ijk) = bch.array() * Dde.col(ijk).array();
                    Dde.col(ijk+1) = ach * bch.array() * Dde.col(ijk+1).array();
                    Ddde.col(ijk*dose_num_tot+ijk+1) = bch.array() * Ddde.col(ijk*dose_num_tot+ijk+1).array();
                    Ddde.col(ijk*dose_num_tot+ijk+1+dose_num_tot) = ach * bch.array() * Ddde.col(ijk*dose_num_tot+ijk+1+dose_num_tot).array();
                    //
                } else if (tform[ijk]=="lin_slope"){
                    beta_c[ijk] = beta_a[ijk] + dbeta[ijk];
                    beta_c[ijk+1] = beta_a[ijk+1] + dbeta[ijk+1];
                    //
                    ArrayXd temp = (df0.col(ijk).array() - beta_c[ijk+1]);
                    ArrayXd temp0 = (df0.col(ijk).array() - beta_c[ijk+1]-dint);
                    ArrayXd temp1 = (df0.col(ijk).array() - beta_c[ijk+1]+dint);
                    ArrayXd temp2 = (temp1.array() < 0).select(0, temp1) - (temp0.array() < 0).select(0, temp0);
                    //
                    De.col(ijk) = beta_c[ijk] * (temp.array() < 0).select(0.0, temp);
                    De.col((ijk+1)) = beta_c[ijk] * (temp.array() < 0).select(0.0, temp);
                    Dde.col(ijk) = (temp.array() < 0).select(0.0, temp);
                    Dde.col((ijk+1)) = beta_c[ijk] * (temp2) / 2/dint;
                    //
                    Ddde.col(ijk * dose_num_tot + (ijk+1)) = (temp2) / 2/dint;
                    Ddde.col((ijk+1) * dose_num_tot + (ijk+1)) = beta_c[ijk] * (temp2) / pow(dint,2);
                    //
                } else if (tform[ijk]=="quad_slope"){
                    beta_c[ijk] = beta_a[ijk] + dbeta[ijk];
                    De.col(ijk) = beta_c[ijk] * df0.col(ijk);
                    //
                } else if (tform[ijk]=="step_slope"){
                    beta_c[ijk] = beta_a[ijk] + dbeta[ijk];
                    beta_c[ijk+1] = beta_a[ijk+1] + dbeta[ijk+1];
                    
                    ArrayXd temp = (df0.col(ijk).array() - beta_0[ijk+1]);
                    ArrayXd temp0 = (df0.col(ijk).array() - beta_0[ijk+1]-dint);
                    ArrayXd temp1 = (df0.col(ijk).array() - beta_0[ijk+1]+dint);
                    ArrayXd temp2 = (temp1.array() < 0).select(0.0, MatrixXd::Zero(temp.rows(),temp.cols()).array()+1.0) - (temp0.array() < 0).select(0.0, MatrixXd::Zero(temp.rows(),temp.cols()).array()+1.0);
                    //
                    De.col(ijk) = beta_0[ijk] * (temp.array() < 0).select(0.0, MatrixXd::Zero(temp.rows(),temp.cols()).array()+1.0);
                    De.col(ijk+1) = beta_0[ijk] * (temp.array() < 0).select(0.0, MatrixXd::Zero(temp.rows(),temp.cols()).array()+1.0);
                    Dde.col(ijk) = (temp.array() < 0).select(0.0, MatrixXd::Zero(temp.rows(),temp.cols()).array()+1.0);
                    Dde.col(ijk+1) = beta_0[ijk] * (temp2) / 2/dint;
                    //
                    Ddde.col(ijk * dose_num_tot + ijk+1) = (temp2) / 2/dint;
                    Ddde.col((ijk+1) * dose_num_tot + ijk+1) = beta_0[ijk] * (temp2) / pow(dint,2);
                    //
                } else {
                    ;
                }
                //
            }
            //
            Dose << De.rowwise().sum();
            //
            if ((modelform=="A")||(modelform=="PA")||(modelform=="PAE")){ //same process used for all of the additive type models
                Te = T0.array().block(0,dose_num_tot,T0.rows(),T0.cols()-dose_num_tot).rowwise().sum() + Dose.array();
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
                            Rdd.col(ij*totalnum+ij) = Ddde.col(ij*dose_num_tot+jk);
                        } else if (ij==jk) {
                            Rdd.col(ij*totalnum+ij) = Tdd0.col(ij);
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
                                if (fir!=0){
                                    if (ij==fir){
                                        Rdd.block(batch_cols[bl],ij*totalnum+jk,batch_cols[bl+1]-batch_cols[bl],1) = Tdd0.block(batch_cols[bl],ij,batch_cols[bl+1]-batch_cols[bl],1).array() * Te.block(batch_cols[bl],0,batch_cols[bl+1]-batch_cols[bl],1).array();
                                    } else {
                                        if (ij<dose_num_tot){
                                            Rdd.block(batch_cols[bl],ij*totalnum+jk,batch_cols[bl+1]-batch_cols[bl],1) = Ddde.block(batch_cols[bl],ij*dose_num_tot+jk,batch_cols[bl+1]-batch_cols[bl],1).array() * T0.block(batch_cols[bl],fir,batch_cols[bl+1]-batch_cols[bl],1).array();
                                        } else {
                                            Rdd.block(batch_cols[bl],ij*totalnum+jk,batch_cols[bl+1]-batch_cols[bl],1) = Tdd0.block(batch_cols[bl],ij,batch_cols[bl+1]-batch_cols[bl],1).array() * T0.block(batch_cols[bl],fir,batch_cols[bl+1]-batch_cols[bl],1).array();
                                        }
                                    }
                                } else {
                                    if (ij<dose_num_tot){
                                        Rdd.block(batch_cols[bl],ij*totalnum+jk,batch_cols[bl+1]-batch_cols[bl],1) = Ddde.block(batch_cols[bl],ij*dose_num_tot+jk,batch_cols[bl+1]-batch_cols[bl],1).array() * Te.block(batch_cols[bl],0,batch_cols[bl+1]-batch_cols[bl],1).array();
                                    } else {
                                        Rdd.block(batch_cols[bl],ij*totalnum+jk,batch_cols[bl+1]-batch_cols[bl],1) = Tdd0.block(batch_cols[bl],ij,batch_cols[bl+1]-batch_cols[bl],1).array() * De.block(batch_cols[bl],0,batch_cols[bl+1]-batch_cols[bl],1).array();
                                    }
                                }
                            } else {
                                if (fir!=0){
                                    if ((ij==fir)||(jk==fir)){
                                        if (ij<dose_num_tot){
                                            Rdd.block(batch_cols[bl],ij*totalnum+jk,batch_cols[bl+1]-batch_cols[bl],1) = Dde.block(batch_cols[bl],ij,batch_cols[bl+1]-batch_cols[bl],1).array() * Td0.block(batch_cols[bl],jk,batch_cols[bl+1]-batch_cols[bl],1).array();
                                            Rdd.block(batch_cols[bl],jk*totalnum+ij,batch_cols[bl+1]-batch_cols[bl],1) = Dde.block(batch_cols[bl],ij,batch_cols[bl+1]-batch_cols[bl],1).array() * Td0.block(batch_cols[bl],jk,batch_cols[bl+1]-batch_cols[bl],1).array();
                                        } else {
                                            Rdd.block(batch_cols[bl],ij*totalnum+jk,batch_cols[bl+1]-batch_cols[bl],1) = Td0.block(batch_cols[bl],ij,batch_cols[bl+1]-batch_cols[bl],1).array() * Td0.block(batch_cols[bl],jk,batch_cols[bl+1]-batch_cols[bl],1).array();
                                            Rdd.block(batch_cols[bl],jk*totalnum+ij,batch_cols[bl+1]-batch_cols[bl],1) = Td0.block(batch_cols[bl],ij,batch_cols[bl+1]-batch_cols[bl],1).array() * Td0.block(batch_cols[bl],jk,batch_cols[bl+1]-batch_cols[bl],1).array();
                                        }
                                    } else if ((ij<dose_num_tot)&&(jk<dose_num_tot)){
                                        Rdd.block(batch_cols[bl],ij*totalnum+jk,batch_cols[bl+1]-batch_cols[bl],1) = Ddde.block(batch_cols[bl],ij*dose_num_tot+jk,batch_cols[bl+1]-batch_cols[bl],1).array() * Dde.block(batch_cols[bl],jk,batch_cols[bl+1]-batch_cols[bl],1).array();
                                        Rdd.block(batch_cols[bl],jk*totalnum+ij,batch_cols[bl+1]-batch_cols[bl],1) = Ddde.block(batch_cols[bl],ij*dose_num_tot+jk,batch_cols[bl+1]-batch_cols[bl],1).array() * Dde.block(batch_cols[bl],jk,batch_cols[bl+1]-batch_cols[bl],1).array();
                                    }
                                }
                            }
                        }
                    }
                }
            }else if (modelform=="M"){
                Te = Te.array() * 0 + 1; //verifies the initial term product is 1
                //
                Te = T0.array().block(0,dose_num_tot,T0.rows(),T0.cols()-dose_num_tot).rowwise().prod() * Dose.array();
                // computes intial risk and derivatives
                R << Te.array();
                Rd = T0.array().pow(-1).array() * Te.colwise().replicate(totalnum).array();
                Rd = Rd.array() * Td0.array();
                #pragma omp parallel for schedule(dynamic) num_threads(nthreads)
                for (int ijk=0;ijk<totalnum;ijk++){
                    if (ijk<dose_num_tot){
                        Rd.col(ijk) = De.col(ijk).array().pow(-1).array()* Te.array();
                        Rd.col(ijk) = Rd.col(ijk).array() * Dde.col(ijk).array();
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
                        if (ij==jk){
                            if (ij<dose_num_tot){
                                Rdd.block(batch_cols[bl],ij*totalnum+jk,batch_cols[bl+1]-batch_cols[bl],1) = Ddde.block(batch_cols[bl],ij*dose_num_tot+jk,batch_cols[bl+1]-batch_cols[bl],1).array() * De.block(batch_cols[bl],jk,batch_cols[bl+1]-batch_cols[bl],1).array().pow(-1).array() * R.block(batch_cols[bl],1,batch_cols[bl+1]-batch_cols[bl],1).array();
                            } else {
                                Rdd.block(batch_cols[bl],ij*totalnum+jk,batch_cols[bl+1]-batch_cols[bl],1) = Tdd0.block(batch_cols[bl],jk,batch_cols[bl+1]-batch_cols[bl],1).array() * Td0.block(batch_cols[bl],jk,batch_cols[bl+1]-batch_cols[bl],1).array().pow(-1).array() * Rd.block(batch_cols[bl],ij,batch_cols[bl+1]-batch_cols[bl],1).array();
                            }
                        } else {
                            if ((ij<dose_num_tot)&&(jk<dose_num_tot)){
                                Rdd.block(batch_cols[bl],ij*totalnum+jk,batch_cols[bl+1]-batch_cols[bl],1) = Ddde.block(batch_cols[bl],ij*dose_num_tot+jk,batch_cols[bl+1]-batch_cols[bl],1).array() * De.block(batch_cols[bl],ij,batch_cols[bl+1]-batch_cols[bl],1).array().pow(-1).array() * R.block(batch_cols[bl],1,batch_cols[bl+1]-batch_cols[bl],1).array();
                                Rdd.block(batch_cols[bl],jk*totalnum+ij,batch_cols[bl+1]-batch_cols[bl],1) = Ddde.block(batch_cols[bl],ij*dose_num_tot+jk,batch_cols[bl+1]-batch_cols[bl],1).array() * De.block(batch_cols[bl],ij,batch_cols[bl+1]-batch_cols[bl],1).array().pow(-1).array() * R.block(batch_cols[bl],1,batch_cols[bl+1]-batch_cols[bl],1).array();
                            } else if ((ij<dose_num_tot)||(jk<dose_num_tot)){
                                Rdd.block(batch_cols[bl],ij*totalnum+jk,batch_cols[bl+1]-batch_cols[bl],1) = Dde.block(batch_cols[bl],jk,batch_cols[bl+1]-batch_cols[bl],1).array() * De.block(batch_cols[bl],jk,batch_cols[bl+1]-batch_cols[bl],1).array().pow(-1).array() * Rd.block(batch_cols[bl],ij,batch_cols[bl+1]-batch_cols[bl],1).array();
                                Rdd.block(batch_cols[bl],jk*totalnum+ij,batch_cols[bl+1]-batch_cols[bl],1) = Dde.block(batch_cols[bl],jk,batch_cols[bl+1]-batch_cols[bl],1).array() * De.block(batch_cols[bl],jk,batch_cols[bl+1]-batch_cols[bl],1).array().pow(-1).array() * Rd.block(batch_cols[bl],ij,batch_cols[bl+1]-batch_cols[bl],1).array();
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
                        Rs3 += Rdd.block(InGroup[i]-1,ij*totalnum+jk,InGroup[i+1]-InGroup[i]+1,1).sum();
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
                }
            }
            #pragma omp parallel for num_threads(nthreads)
            for (int ijk=0;ijk<totalnum*(totalnum+1)/2;ijk++){//totalnum*(totalnum+1)/2
                int ij = 0;
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
        if (iteration > totalnum*2){
            if (iteration % (totalnum)){
                if (Lld_worst < deriv_epsilon){
                    iteration = maxiter;
                }
                Ll_comp[1]=Ll[0];
                if (abs(Ll_comp[1]-Ll_comp[0])<1){
                    abs_max = abs_max*0.1;
                }
                if (abs_max < epsilon/10){
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
    const Map<MatrixXd> Lldd_mat(as<Map<MatrixXd> >(Lldd_vec));
    MatrixXd Lldd_inv = Lldd_mat.inverse().matrix();
    //
    vector <double> Std_Er(totalnum,0);
    
    for (int ijk=0;ijk<totalnum;ijk++){
        Std_Er[ijk] = -1*Lldd_inv.coeff(ijk,ijk);
    }
    //
    List res_list = List::create(_["LogLik"]=wrap(Ll),_["First_Der"]=wrap(Lld),_["Second_Der"]=Lldd_vec,_["beta_0"]=wrap(beta_0),_["AIC"]=2*totalnum-2*Ll[fir], _["Standard_Error"] = wrap(Std_Er), _["Information_Matrix"] = wrap(Lldd_inv));
    //
    return res_list;
}

//List PEANUT_PLOT( VectorXd beta_linT,VectorXd beta_loglinT,VectorXd beta_plinT, VectorXd beta_dose,MatrixXd df_lin,MatrixXd df_loglin,MatrixXd df_plin, MatrixXd df_dose,int fir,string modelform, string doseform, StringVector dose_terms,int ntime, NumericVector include_bool, int batch_size){
//    //    cout << "start Risk" << endl;
//    //
//    srand (time(NULL));
//    //
//    cout << ntime << endl;
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
//    int totalnum = dose_num_tot;
//    //
//    cout.precision(10); //forces higher precision numbers printed to terminal
//    int nthreads = Eigen::nbThreads(); //stores how many threads are allocated
//    //
//    VectorXd beta_lin;
//    VectorXd beta_loglin; //The vectors of parameter used
//    VectorXd beta_plin;
//        //
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
//    double totem = df_loglin.rows();//precalculates how many rows
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
//    //for (int ij=0; ij < batch_cols.size(); ij++){
//    //    cout << batch_cols[ij] << endl;
//    //}
//    end_point = system_clock::now();
//    end = time_point_cast<microseconds>(end_point).time_since_epoch().count();
//    cout<<"df99,"<<(end-start)<<",Batches"<<endl;
//    gibtime = system_clock::to_time_t(system_clock::now());
//    cout << ctime(&gibtime) << endl;
//    // ---------------------------------------------
//    // To Start, needs to seperate the derivative terms
//    // ---------------------------------------------
//    //
//    //
//    VectorXd beta_0(totalnum);
//    MatrixXd df0 = MatrixXd::Zero(df_lin.rows(), totalnum); // stores memory for the derivative term parameters and columns
//    MatrixXd T0 = MatrixXd::Zero(df0.rows(), totalnum); //preallocates matrix for Derivative column terms
//    MatrixXd Td0 = MatrixXd::Zero(df0.rows(), totalnum); //preallocates matrix for Derivative column terms
//    MatrixXd Tdd0 = MatrixXd::Zero(df0.rows(), totalnum); //preallocates matrix for Derivative column terms
//    #pragma omp parallel for num_threads(nthreads)
//    for (int ij=0;ij<totalnum;ij++){
//        int ind0 = ij;
//        if (ind0 < dose_num_tot){
//            beta_0[ij] = beta_dose[ind0];
//            df0.col(ij) = df_dose.col(ind0);
//            tform[ij] = dose_terms[ij];
//        } else {
//            ind0 = ind0 - dose_num_tot;
//            if (include_bool[0]==1){
//                if (ind0 < beta_lin.size()){
//                    // one exists and is one
//                    beta_0[ij] = beta_lin[ind0];
//                    df0.col(ij) = df_lin.col(ind0);
//                    tform[ij] = "lin";
//                    //
//                } else {
//                    //one exists and its not one
//                    ind0 = ind0 - beta_lin.size();
//                    if (include_bool[1]==1){
//                        if (ind0 < beta_loglin.size()){
//                            //one and two exists and is two
//                            beta_0[ij] = beta_loglin[ind0];
//                            df0.col(ij) = df_loglin.col(ind0);
//                            tform[ij] = "loglin";
//                            //
//                        } else{
//                            //one exists, two does, must be three
//                            if (include_bool[2]!=1){
//                                throw invalid_argument( "Are all three used?" );
//                            }
//                            ind0 = ind0 - beta_loglin.size();
//                            beta_0[ij] = beta_plin[ind0];
//                            df0.col(ij) = df_plin.col(ind0);
//                            tform[ij] = "plin";
//                            //
//                        }
//                    } else{
//                        //one exists, and two doesn't exist, must be three
//                        if (include_bool[2]!=1){
//                            throw invalid_argument( "Are all first and third used?" );
//                        }
//                        beta_0[ij] = beta_plin[ind0];
//                        df0.col(ij) = df_plin.col(ind0);
//                        tform[ij] = "plin";
//                        //
//                    }
//                }
//            }else{
//                //one doesn't exist
//                if (include_bool[1]==1){
//                    if (ind0 < beta_loglin.size()){
//                        //one doesn't exist and two exists and is two
//                        beta_0[ij] = beta_loglin[ind0];
//        //                cout << ind0 << ", " << beta_0 << ", " << beta_loglin.transpose() << endl;
//                        df0.col(ij) = df_loglin.col(ind0);
//                        tform[ij] = "loglin";
//                        //
//                    } else{
//                        //one doesn't exist, two does, must be three
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
//                    //one doesn't exist, and two doesn't exist, must be three
//                    if (include_bool[2]!=1){
//                        throw invalid_argument( "Are all first and third used?" );
//                    }
//                    beta_0[ij] = beta_plin[ind0];
//                    df0.col(ij) = df_plin.col(ind0);
//                    tform[ij] = "plin";
//                    //
//                }
//            }
//        }
//        //        cout << df0.col(ij).array().transpose() << endl;
//        T0.col(ij) = (df0.col(ij).array() * beta_0[ij]).matrix();
//        if (tform[ij]=="lin") {
//            Td0.col(ij) = df0.col(ij);
//        } else if (tform[ij]=="loglin") {
//
//            T0.col(ij) = T0.col(ij).array().exp();
//            Td0.col(ij) = df0.col(ij).array() * T0.col(ij).array();
//            Tdd0.col(ij) = df0.col(ij).array() * Td0.col(ij).array();
//        } else if (tform[ij]=="plin") {
//            T0.col(ij) = 1 + T0.col(ij).array();
//            Td0.col(ij) = df0.col(ij);
//        } else {
//            cout << tform[ij] << " is invalid" << endl;
//            throw invalid_argument( "Invalid term type" );
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
//    MatrixXd De = MatrixXd::Zero(df0.rows(), 1); //preallocates matrix for non-Derivative column terms
//    MatrixXd Dde = MatrixXd::Zero(df0.rows(), dose_num_tot); //preallocates matrix for non-Derivative column terms
//    MatrixXd Ddde = MatrixXd::Zero(df0.rows(), pow(dose_num_tot,2)); //preallocates matrix for non-Derivative column terms
//    MatrixXd R = MatrixXd::Zero(df0.rows(), 1); //preallocates matrix for Risks and derivatives
//    MatrixXd Rd = MatrixXd::Zero(df0.rows(), totalnum); //preallocates matrix for Risks and derivatives
//    MatrixXd Rdd = MatrixXd::Zero(df0.rows(), pow(totalnum,2));
//    if (doseform=="A"){
//        De = T0.array().block(0,0,T0.rows(),dose_num_tot).rowwise().sum();
//        Dde = Td0.array().block(0,0,T0.rows(),dose_num_tot);
//        #pragma omp parallel for num_threads(nthreads)
//        for (int ij=0;ij<dose_num_tot;ij++){
//            Ddde.col(ij*dose_num_tot+ij) = Tdd0.col(ij);
//        }
//    } else {
//        throw invalid_argument( "That Dose model isn't implemented" );
//    }
//    if ((modelform=="A")||(modelform=="PA")||(modelform=="PAE")){ //same process used for all of the additive type models
//        Te = T0.array().block(0,dose_num_tot,T0.rows(),T0.cols()-dose_num_tot).rowwise().sum() + De.array();
//        // computes intial risk and derivatives
//        if (modelform=="A"){
//            R << Te.array();
//            Rd << Dde.array(), Td0.array().block(0,dose_num_tot,T0.rows(),T0.cols()-dose_num_tot);
//            #pragma omp parallel for num_threads(nthreads)
//            for (int ij=0;ij<totalnum;ij++){
//                if (ij<dose_num_tot){
//                    Rdd.col(ij*totalnum+ij) = Ddde.col(ij*dose_num_tot + ij);
//                } else {
//                    Rdd.col(ij*totalnum+ij) = Tdd0.col(ij);
//                }
//            }
//        } else if ((modelform=="PAE")||(modelform=="PA")){
//            if (fir!=0){
//                Te = Te.array() - T0.col(fir).array();
//            } else {
//                Te = Te.array() - De.array();
//            }
//            if (modelform=="PAE"){
//                Te = Te.array() + 1;
//            }
//            if (fir!=0){
//                R << T0.col(fir).array() * Te.array();
//                Rd << Td0.array() * T0.col(fir).array();//, Td0.col(0).array() * Te.array(), Td0.col(1).array() * Te.array();
//                Rd.col(fir) = Td0.col(fir).array() * Te.array();
//            } else {
//                R << De.array() * Te.array();
//                Rd << Td0.array() * De.array();
//                #pragma omp parallel for schedule(dynamic) num_threads(nthreads)
//                for (int ij=0;ij<dose_num_tot;ij++){
//                    Rd.col(ij) = Dde.col(ij).array() * Te.array();
//                }
//            }
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
//                        if (fir!=0){
//                            if (ij==fir){
//                                Rdd.block(batch_cols[bl],ij*totalnum+jk,batch_cols[bl+1]-batch_cols[bl],1) = Tdd0.block(batch_cols[bl],ij,batch_cols[bl+1]-batch_cols[bl],1).array() * Te.block(batch_cols[bl],0,batch_cols[bl+1]-batch_cols[bl],1).array();
//                            } else {
//                                if (ij<dose_num_tot){
//                                    Rdd.block(batch_cols[bl],ij*totalnum+jk,batch_cols[bl+1]-batch_cols[bl],1) = Ddde.block(batch_cols[bl],ij,batch_cols[bl+1]-batch_cols[bl],1).array() * T0.block(batch_cols[bl],fir,batch_cols[bl+1]-batch_cols[bl],1).array();
//                                } else {
//                                    Rdd.block(batch_cols[bl],ij*totalnum+jk,batch_cols[bl+1]-batch_cols[bl],1) = Tdd0.block(batch_cols[bl],ij,batch_cols[bl+1]-batch_cols[bl],1).array() * T0.block(batch_cols[bl],fir,batch_cols[bl+1]-batch_cols[bl],1).array();
//                                }
//                            }
//                        } else {
//                            if (ij<dose_num_tot){
//                                Rdd.block(batch_cols[bl],ij*totalnum+jk,batch_cols[bl+1]-batch_cols[bl],1) = Ddde.block(batch_cols[bl],ij,batch_cols[bl+1]-batch_cols[bl],1).array() * Te.block(batch_cols[bl],0,batch_cols[bl+1]-batch_cols[bl],1).array();
//                            } else {
//                                Rdd.block(batch_cols[bl],ij*totalnum+jk,batch_cols[bl+1]-batch_cols[bl],1) = Tdd0.block(batch_cols[bl],ij,batch_cols[bl+1]-batch_cols[bl],1).array() * De.block(batch_cols[bl],0,batch_cols[bl+1]-batch_cols[bl],1).array();
//                            }
//                        }
//                    } else {
//                        if (fir!=0){
//                            if ((ij==fir)||(jk==fir)){
//                                if (ij<dose_num_tot){
//                                    Rdd.block(batch_cols[bl],ij*totalnum+jk,batch_cols[bl+1]-batch_cols[bl],1) = Dde.block(batch_cols[bl],ij,batch_cols[bl+1]-batch_cols[bl],1).array() * Td0.block(batch_cols[bl],jk,batch_cols[bl+1]-batch_cols[bl],1).array();
//                                    Rdd.block(batch_cols[bl],jk*totalnum+ij,batch_cols[bl+1]-batch_cols[bl],1) = Dde.block(batch_cols[bl],ij,batch_cols[bl+1]-batch_cols[bl],1).array() * Td0.block(batch_cols[bl],jk,batch_cols[bl+1]-batch_cols[bl],1).array();
//                                } else {
//                                    Rdd.block(batch_cols[bl],ij*totalnum+jk,batch_cols[bl+1]-batch_cols[bl],1) = Td0.block(batch_cols[bl],ij,batch_cols[bl+1]-batch_cols[bl],1).array() * Td0.block(batch_cols[bl],jk,batch_cols[bl+1]-batch_cols[bl],1).array();
//                                    Rdd.block(batch_cols[bl],jk*totalnum+ij,batch_cols[bl+1]-batch_cols[bl],1) = Td0.block(batch_cols[bl],ij,batch_cols[bl+1]-batch_cols[bl],1).array() * Td0.block(batch_cols[bl],jk,batch_cols[bl+1]-batch_cols[bl],1).array();
//                                }
//                            } else if ((ij<dose_num_tot)&&(jk<dose_num_tot)){
//                                Rdd.block(batch_cols[bl],ij*totalnum+jk,batch_cols[bl+1]-batch_cols[bl],1) = Dde.block(batch_cols[bl],ij,batch_cols[bl+1]-batch_cols[bl],1).array() * Dde.block(batch_cols[bl],jk,batch_cols[bl+1]-batch_cols[bl],1).array();
//                                Rdd.block(batch_cols[bl],jk*totalnum+ij,batch_cols[bl+1]-batch_cols[bl],1) = Dde.block(batch_cols[bl],ij,batch_cols[bl+1]-batch_cols[bl],1).array() * Dde.block(batch_cols[bl],jk,batch_cols[bl+1]-batch_cols[bl],1).array();
//                            }
//                        }
//                    }
//                }
//            }
//        }
//    }else if (modelform=="M"){
//        Te = Te.array() * 0 + 1; //verifies the initial term product is 1
//        //
//        Te = T0.array().block(0,dose_num_tot,T0.rows(),T0.cols()-dose_num_tot).rowwise().prod() * De.array();
//
//        // computes intial risk and derivatives
//        R << Te.array();
//        #pragma omp parallel for schedule(dynamic) num_threads(nthreads)
//        for (int ijk=0;ijk<totalnum;ijk++){
//            if (ijk<dose_num_tot){
//                Rd.col(ijk) = De.col(ijk).array().pow(-1).array()* Te.array();
//                Rd.col(ijk) = Rd.col(ijk).array() * Dde.col(ijk).array();
//            }
//        }
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
//                    if (ij<dose_num_tot){
//                        Rdd.block(batch_cols[bl],ij*totalnum+jk,batch_cols[bl+1]-batch_cols[bl],1) = Ddde.block(batch_cols[bl],jk,batch_cols[bl+1]-batch_cols[bl],1).array() * Dde.block(batch_cols[bl],jk,batch_cols[bl+1]-batch_cols[bl],1).array().pow(-1).array() * Rd.block(batch_cols[bl],ij,batch_cols[bl+1]-batch_cols[bl],1).array();
//                    } else {
//                        Rdd.block(batch_cols[bl],ij*totalnum+jk,batch_cols[bl+1]-batch_cols[bl],1) = Tdd0.block(batch_cols[bl],jk,batch_cols[bl+1]-batch_cols[bl],1).array() * Td0.block(batch_cols[bl],jk,batch_cols[bl+1]-batch_cols[bl],1).array().pow(-1).array() * Rd.block(batch_cols[bl],ij,batch_cols[bl+1]-batch_cols[bl],1).array();
//                    }
//                } else {
//                    if ((ij<dose_num_tot)||(ij<dose_num_tot)){
//                        Rdd.block(batch_cols[bl],ij*totalnum+jk,batch_cols[bl+1]-batch_cols[bl],1) = Dde.block(batch_cols[bl],ij,batch_cols[bl+1]-batch_cols[bl],1).array() * De.block(batch_cols[bl],ij,batch_cols[bl+1]-batch_cols[bl],1).array().pow(-1).array() * Rd.block(batch_cols[bl],ij,batch_cols[bl+1]-batch_cols[bl],1).array();
//                        Rdd.block(batch_cols[bl],jk*totalnum+ij,batch_cols[bl+1]-batch_cols[bl],1) = Dde.block(batch_cols[bl],ij,batch_cols[bl+1]-batch_cols[bl],1).array() * De.block(batch_cols[bl],ij,batch_cols[bl+1]-batch_cols[bl],1).array().pow(-1).array() * Rd.block(batch_cols[bl],jk,batch_cols[bl+1]-batch_cols[bl],1).array();
//                    } else{
//                        Rdd.block(batch_cols[bl],ij*totalnum+jk,batch_cols[bl+1]-batch_cols[bl],1) = Td0.block(batch_cols[bl],jk,batch_cols[bl+1]-batch_cols[bl],1).array() * T0.block(batch_cols[bl],jk,batch_cols[bl+1]-batch_cols[bl],1).array().pow(-1).array() * Rd.block(batch_cols[bl],ij,batch_cols[bl+1]-batch_cols[bl],1).array();
//                        Rdd.block(batch_cols[bl],jk*totalnum+ij,batch_cols[bl+1]-batch_cols[bl],1) = Td0.block(batch_cols[bl],ij,batch_cols[bl+1]-batch_cols[bl],1).array() * T0.block(batch_cols[bl],ij,batch_cols[bl+1]-batch_cols[bl],1).array().pow(-1).array() * Rd.block(batch_cols[bl],jk,batch_cols[bl+1]-batch_cols[bl],1).array();
//                    }
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
//    // -------------------------------------------------------------------------------------------
//    string line;
//    ifstream infile("test.txt"); //The file of risk group rows
//    vector<string>  RiskGroup(ntime); //vector of strings detailing the rows
//    IntegerMatrix RiskFail(ntime,2); //vector giving the event rows
//    int j_iter = 0;
//    //
//    end_point = system_clock::now();
//    end = time_point_cast<microseconds>(end_point).time_since_epoch().count();
//    cout<<"df100 "<<(end-start)<<" "<<0<<" "<<0<<" "<<-1<<",Prep_R"<<endl;
//    gibtime = system_clock::to_time_t(system_clock::now());
//    cout << ctime(&gibtime) << endl;
//    // --------------------------
//    // A file is created previously that stores what rows belong to which risk sections
//    // --------------------------
//    while (getline(infile, line)){
//        istringstream iss(line);
//        vector <int> row;
//        string lineStream;
//        while (getline(iss, lineStream, ',')){  
//            row.push_back(stoi(lineStream));
//        }
//        vector<int> rows = vector<int>(row.end() - 2, row.end());
//        cout << rows.size() << ", " << rows[0] << ", " << rows[1] << ", " << row[row.size()-1] << endl;
//        // --------------------------
//        // needs the rows of every event with the same end time
//        // --------------------------
//        RiskFail(j_iter,0)=rows[0]-1;
//        RiskFail(j_iter,1)=rows[1]-1;
//        ostringstream group_str;
//        copy(row.begin(), row.end()-2, ostream_iterator<int>(group_str, ","));
//        RiskGroup[j_iter] = group_str.str();
//        cout << group_str.str() << endl;
//        cout << j_iter << endl;
//        j_iter++;
//    }
//    //
//    // --------------------------
//    // now a vector exists with row locations
//    // --------------------------
//    //The log-likelihood is calculated in parallel over the risk groups
//    vector<double> baseline(ntime,0.0);
//    #pragma omp parallel for schedule(dynamic) num_threads(nthreads)
//    for (int j=0;j<ntime;j++){
//        double Rs1 = 0;
//            cout << 1 << endl;
//        //
//        vector<int> InGroup;
//        string Groupstr = RiskGroup[j];
//        stringstream ss(Groupstr);
//        //
//        for (int i; ss >> i;) {
//            InGroup.push_back(i);    
//            if (ss.peek() == ',')
//                ss.ignore();
//        }
//            cout << 1 << endl;
//        //Now has the grouping pairs
//        int dj = RiskFail(j,1)-RiskFail(j,0)+1;
//            cout << 2 << endl;
//            cout << Groupstr << endl;
//            cout << InGroup.size() << endl;
//            int at_risk = 0;
//        for (int i = 0; i < InGroup.size()-1; i=i+2){
//                cout << InGroup[i] << "," << InGroup[i+1] << "," << ij << "," << jk << endl;
//            Rs1 += R.block(InGroup[i]-1,0,InGroup[i+1]-InGroup[i]+1,1).sum();
//                at_risk += InGroup[i+1]-InGroup[i]+1;
//        }
//        baseline[j] = dj / Rs1;
//    }
//    List res_list = List::create(_["baseline"]=wrap(baseline));
//    //
//    return res_list;
//}
