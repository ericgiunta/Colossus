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
void Write_Ind_File(NumericMatrix df, NumericVector tu){
    //----------------------------------------------------------------------------------------------------------------//
    const Map<MatrixXd> df_m(as<Map<MatrixXd> >(df));
    vector <string> ind_list(tu.size(),"");
    cout.precision(10); //forces higher precision numbers printed to terminal
    int nthreads = Eigen::nbThreads()-1; //stores how many threads are allocated
    cout << "Start Write" << endl;
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
    cout << "End Write" << endl;
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
List Martingale_Residuals(StringVector IDS, NumericMatrix dfr, NumericVector event_times){
    const Map<MatrixXd> df_data(as<Map<MatrixXd> >(dfr));
//    cout << df_data.rows() << " " << df_data.cols() << endl;   
    cout.precision(10); //forces higher precision numbers printed to terminal
    int nthreads = Eigen::nbThreads()-1; //stores how many threads are allocated
    vector<float> vv(40); //stores the covariate values
    double dx = (max(event_times)-min(event_times))/vv.size();//varies from max to minimum
    vv[0] = min(event_times);
    generate(vv.begin(), vv.end(), [n = 0, &dx]() mutable { return n++ * dx; });
    int tot_size = vv.size()*IDS.size();
    VectorXd times = VectorXd::Zero(tot_size);
    VectorXd MRes = VectorXd::Zero(tot_size);
//    cout << MRes.rows() << endl;
    #pragma omp parallel for schedule(dynamic) num_threads(nthreads) collapse(2)
    for (int ij=0;ij<(vv.size());ij++){
        for (int jk=0;jk<(IDS.size());jk++){
            int ind = jk*vv.size() + ij;
//            cout << ind << endl;
            times(ind) = event_times[ij];
//            cout << ind << endl;
            if ((df_data.col(0).array() < vv[ij]).count()>0){
//                cout << ind << endl;
                auto mask = (df_data.col(0).array() < vv[ij]);
//                cout << ind << endl;
                MRes(ind) = mask.select(df_data.col(1), 0).sum();
            }
        }
    }
    List temp_list = List::create(_["times"]=wrap(times),_["Martingale_Error"]=wrap(MRes));
    return temp_list;
}


// [[Rcpp::export]]
List peanut_transition(NumericVector a_lin,NumericVector a_loglin,NumericVector a_plin, NumericMatrix x_lin, NumericMatrix x_loglin, NumericMatrix x_plin, NumericMatrix x_dose,int fir,string modelform,int ntime, NumericVector include_bool, List Control, List Dose_paras, NumericMatrix df_groups, NumericVector tu){
    //----------------------------------------------------------------------------------------------------------------//
    Map<VectorXd> beta_lin(as<Map<VectorXd> >(a_lin));
    Map<VectorXd> beta_loglin(as<Map<VectorXd> >(a_loglin));
    Map<VectorXd> beta_plin(as<Map<VectorXd> >(a_plin));
    const Map<MatrixXd> df_lin(as<Map<MatrixXd> >(x_lin));
    const Map<MatrixXd> df_loglin(as<Map<MatrixXd> >(x_loglin));
    const Map<MatrixXd> df_plin(as<Map<MatrixXd> >(x_plin));
    const Map<MatrixXd> df_dose(as<Map<MatrixXd> >(x_dose));
    //
//    const Map<MatrixXd> test_df(as<Map<MatrixXd> >(x_lin));
//    const SparseMatrix<double> test_sp = test_df.sparseView();
    //
    // Converts from Rcpp types to efficient Eigen types
    bool change_all = Control["change_all"];
    double lr = Control["lr"];
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
    List res = LogLik_PEANUT(beta_lin,beta_loglin,beta_plin,df_lin,df_loglin,df_plin, df_dose,fir,modelform,ntime,include_bool, lr, maxiter, halfmax, epsilon, dbeta_cap, abs_max,dose_abs_max, deriv_epsilon,beta_loglin_slope, beta_loglin_top , beta_lin_slope , beta_lin_int , beta_quad , beta_step_slope , beta_step_int, df_groups, tu, change_all );
    //----------------------------------------------------------------------------------------------------------------//
    return res;
}

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
    List res = LogLik_Bounds(q1, beta_lin,beta_loglin,beta_plin,df_lin,df_loglin,df_plin, df_dose,fir,modelform,ntime,include_bool, lr, maxiter, halfmax, epsilon, dbeta_cap, abs_max, dose_abs_max, deriv_epsilon,beta_loglin_slope, beta_loglin_top , beta_lin_slope , beta_lin_int , beta_quad , beta_step_slope , beta_step_int, df_groups,tu);
    //----------------------------------------------------------------------------------------------------------------//
    return res;
}

// [[Rcpp::export]]
List peanut_plot(NumericVector a_lin,NumericVector a_loglin,NumericVector a_plin, NumericMatrix x_lin, NumericMatrix x_loglin, NumericMatrix x_plin, NumericMatrix x_dose,int fir,string modelform,int ntime, NumericVector include_bool, List Control, List Dose_paras, NumericMatrix df_groups, NumericVector tu , vector<string> Plot_Type, int uniq_v){
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
	double dose_abs_max = Control["dose_abs_max"];
	double deriv_epsilon =Control["deriv_epsilon"];
	List beta_loglin_slope = Dose_paras["beta_loglin_slope"];
    List beta_loglin_top  = Dose_paras["beta_loglin_top"];
    List beta_lin_slope  = Dose_paras["beta_lin_slope"];
    List beta_lin_int  = Dose_paras["beta_lin_int"];
    List beta_quad  = Dose_paras["beta_quad"];
    List beta_step_slope  = Dose_paras["beta_step_slope"];
    List beta_step_int  = Dose_paras["beta_step_int"];
    // Converts from Rcpp types to efficient Eigen types
    List res;
    // there are two types of plots that can be generated
    //----------------------------------------------------------------------------------------------------------------//
    if (Plot_Type[0]=="SURV"){
        res = PEANUT_PLOT_SURV(beta_lin,beta_loglin,beta_plin,df_lin,df_loglin,df_plin, df_dose,fir,modelform,ntime,include_bool,beta_loglin_slope, beta_loglin_top , beta_lin_slope , beta_lin_int , beta_quad , beta_step_slope , beta_step_int, df_groups,tu,dose_abs_max);
    }else if (Plot_Type[0]=="RISK"){
        res = PEANUT_PLOT_RISK(beta_lin,beta_loglin,beta_plin,df_lin,df_loglin,df_plin, df_dose,fir,modelform,ntime,include_bool,beta_loglin_slope, beta_loglin_top , beta_lin_slope , beta_lin_int , beta_quad , beta_step_slope , beta_step_int, df_groups,tu, uniq_v,dose_abs_max);
    } else {
        throw invalid_argument("Invalid plot type");
    }
    return res;
}

// [[Rcpp::export]]
NumericMatrix peanut_schoenfeld_transition( NumericVector a_lin,NumericVector a_loglin,NumericVector a_plin, NumericMatrix x_lin, NumericMatrix x_loglin, NumericMatrix x_plin, NumericMatrix x_dose,int fir,string modelform,int ntime, NumericVector include_bool, List Control, List Dose_paras, NumericMatrix df_groups, NumericVector tu){
        //----------------------------------------------------------------------------------------------------------------//
    Map<VectorXd> beta_lin(as<Map<VectorXd> >(a_lin));
    Map<VectorXd> beta_loglin(as<Map<VectorXd> >(a_loglin));
    Map<VectorXd> beta_plin(as<Map<VectorXd> >(a_plin));
    const Map<MatrixXd> df_lin(as<Map<MatrixXd> >(x_lin));
    const Map<MatrixXd> df_loglin(as<Map<MatrixXd> >(x_loglin));
    const Map<MatrixXd> df_plin(as<Map<MatrixXd> >(x_plin));
    const Map<MatrixXd> df_dose(as<Map<MatrixXd> >(x_dose));
	double dose_abs_max = Control["dose_abs_max"];
    // Converts from Rcpp types to efficient Eigen types
	List beta_loglin_slope = Dose_paras["beta_loglin_slope"];
    List beta_loglin_top  = Dose_paras["beta_loglin_top"];
    List beta_lin_slope  = Dose_paras["beta_lin_slope"];
    List beta_lin_int  = Dose_paras["beta_lin_int"];
    List beta_quad  = Dose_paras["beta_quad"];
    List beta_step_slope  = Dose_paras["beta_step_slope"];
    List beta_step_int  = Dose_paras["beta_step_int"];
    //
    //----------------------------------------------------------------------------------------------------------------//
    NumericMatrix res = Schoenfeld_PEANUT(beta_lin,beta_loglin,beta_plin,df_lin,df_loglin,df_plin, df_dose,fir,modelform,ntime,include_bool,beta_loglin_slope, beta_loglin_top , beta_lin_slope , beta_lin_int , beta_quad , beta_step_slope , beta_step_int, df_groups,tu, dose_abs_max);
    //----------------------------------------------------------------------------------------------------------------//
    return res;
}


List LogLik_PEANUT( VectorXd beta_linT,VectorXd beta_loglinT,VectorXd beta_plinT,MatrixXd df_lin,MatrixXd df_loglin,MatrixXd df_plin, MatrixXd df_dose,int fir,string modelform,int ntime, NumericVector include_bool, double lr, int maxiter, int halfmax, double epsilon, double dbeta_cap, double abs_max,double dose_abs_max, double deriv_epsilon,List beta_loglin_slopes, List beta_loglin_tops , List beta_lin_slopes , List beta_lin_ints , List beta_quads , List beta_step_slopes , List beta_step_ints, NumericMatrix df_groups, NumericVector tu, bool change_all){
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
    vector<int> dose_breaks(df_dose.cols(),0);
    int dose_num_tot=0;
    int dose_term_tot=0;
    // Converts from a list of Rcpp vectors to a std::vector of std::vectors
    // There were possible issues with Rcpp vectors being created in OMP sections
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
        // Establishes the number of dose terms needed
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
        // Gathers the dose terms per dose column
        dose_term_tot += dose_breaks[ijk];
        //
    }
    //
    int totalnum = dose_num_tot;
    //
    cout.precision(10); //forces higher precision numbers printed to terminal
    int nthreads = Eigen::nbThreads()-1; //stores how many threads are allocated
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
    //
    double Lld_worst = 0.0; //stores derivative value used to determine if every parameter is near convergence
    vector <string> tform(totalnum);// list of term types
    double totem = df_loglin.rows();//precalculates how many rows are needed
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
    //
    VectorXd beta_0(totalnum);
    MatrixXd df0 = MatrixXd::Zero(df_lin.rows(), totalnum); // stores memory for the parameter columns
    MatrixXd T0 = MatrixXd::Zero(df0.rows(), totalnum); //preallocates matrix for Term column
    MatrixXd Td0 = MatrixXd::Zero(df0.rows(), totalnum); //preallocates matrix for Term derivative columns
    MatrixXd Tdd0 = MatrixXd::Zero(df0.rows(), totalnum); //preallocates matrix for Term second derivative columns
    //
    MatrixXd Te = MatrixXd::Zero(df0.rows(), 1); //preallocates matrix for non-Derivative column terms
    MatrixXd R = MatrixXd::Zero(df0.rows(), 1); //preallocates matrix for Risks
    MatrixXd Rd = MatrixXd::Zero(df0.rows(), totalnum); //preallocates matrix for Risk derivatives
    MatrixXd Rdd = MatrixXd::Zero(df0.rows(), totalnum*(totalnum+1)/2); //preallocates matrix for Risk second derivatives
    //
    MatrixXd De = MatrixXd::Zero(df_dose.rows(),dose_num_tot); //matrix of dose term values
    MatrixXd Dde = MatrixXd::Zero(df_dose.rows(),dose_num_tot); //matrix of dose term derivatives
    MatrixXd Ddde = MatrixXd::Zero(df_dose.rows(),dose_num_tot*(dose_num_tot+1)/2); //matrix of dose term second derivatives
    double dint = dose_abs_max; //The amount of change used to calculate derivatives in threshold paramters
    int total_dose=0; //used later on for a section summing the dose terms
    //
    VectorXd Dose = VectorXd::Zero(df_dose.rows()); //Matrix of the total dose term values
    #pragma omp declare reduction (eig_plus: VectorXd: omp_out=omp_out+omp_in) initializer(omp_priv=VectorXd::Zero(omp_orig.size()))
    #pragma omp parallel for schedule(dynamic) num_threads(nthreads) reduction(eig_plus:Dose)
    for (int ij=0;ij<(totalnum-dose_num_tot+dose_term_tot);ij++){
        if (ij < dose_term_tot){//Dose terms are indexed originally by complete term, i.e. step function not the slope and threshold seperately
            int ind0 = ij;
            int ijk=0;
            while (ind0>dose_breaks[ijk]){
                ind0=ind0 - dose_breaks[ijk];
                ijk++;
            } //determines which dose column is used
            //
            int loglin_size=0;
            int lin_size=0;
            int quad_size=0;
            int step_size=0;
            vector<double> beta_loglin_slope = beta_loglin_slopes_CPP[ijk];
            vector<double> beta_loglin_top = beta_loglin_tops_CPP[ijk];
            vector<double> beta_lin_slope = beta_lin_slopes_CPP[ijk];
            vector<double> beta_lin_int = beta_lin_ints_CPP[ijk];
            vector<double> beta_quad = beta_quads_CPP[ijk];
            vector<double> beta_step_slope = beta_step_slopes_CPP[ijk];
            vector<double> beta_step_int = beta_step_ints_CPP[ijk];
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
            int dub_off=0;
            if ((beta_loglin_slope[0]==1)&&(loglin_size==1)){
                dub_off=1; //if the exponential term has a slope of 1, it isn't changed. This corrects the number of parameters
            }
            if (ind0 < loglin_size){
                ArrayXd temp = (beta_loglin_top[ind0] * df_dose.col(ijk)).array().exp();
                ArrayXd temp1 = beta_loglin_slope[ind0] * temp;
                // in both cases the total exponential term is calculated
                //
                if ((beta_loglin_slope[ind0]==1)&&(loglin_size==1)){
                    int ind1 = cumulative_dose_num[ijk]+ind0; //skips slope term
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
                    int ind1 = cumulative_dose_num[ijk]+2*ind0; //does not skip dose terms
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
                int jk = ind0 - loglin_size; // The linear threshold model has no closed form equation for derivatives
                ArrayXd temp = (df_dose.col(ijk).array() - beta_lin_int[jk]);
                ArrayXd temp0 = (df_dose.col(ijk).array() - beta_lin_int[jk]-dint);//derivatives are approximated
                ArrayXd temp1 = (df_dose.col(ijk).array() - beta_lin_int[jk]+dint);
                //
                int ind1 = cumulative_dose_num[ijk]+2*loglin_size - dub_off + 2*jk;
//                int ind2 = (ind1 + 1); 
                //
                beta_0[ind1] = beta_lin_slope[jk];
                beta_0[(ind1 + 1)] = beta_lin_int[jk];
                tform[ind1] = "lin_slope";
                tform[(ind1 + 1)] = "lin_ind";
                df0.col(ind1) = df_dose.col(ijk);
                df0.col((ind1 + 1)) = df_dose.col(ijk);
                //
                temp = (temp.array() < 0).select(0.0, temp);
                temp0 = (temp0.array() < 0).select(0.0, temp0);
                temp1 = (temp1.array() < 0).select(0.0, temp1); //applies thresholds
                //
                De.col(ind1) = beta_lin_slope[jk] * temp.array();
                De.col((ind1 + 1)) = beta_lin_slope[jk] * temp.array();
                Dose = Dose.array() + De.col(ind1).array();
                Dde.col(ind1) = temp.array();
                Dde.col((ind1 + 1)) = beta_lin_slope[jk] * (temp1.array() -temp0.array()) / 2.0/dint;
                //
                Ddde.col((ind1 + 1) * ((ind1 + 1)+1)/2 + ind1) = (temp1.array() -temp0.array()) / 2.0/dint;
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
                ArrayXd temp = (df_dose.col(ijk).array() - beta_step_int[jk]);// Step function also has approximated derivatives
                ArrayXd temp0 = (df_dose.col(ijk).array() - beta_step_int[jk]-dint);
                ArrayXd temp1 = (df_dose.col(ijk).array() - beta_step_int[jk]+dint);
                //
                int ind1 = cumulative_dose_num[ijk]+2*loglin_size - dub_off  + 2*lin_size + quad_size + 2*jk;
//                int ind2 = ind1 + 1;
                //
                beta_0[ind1] = beta_step_slope[jk];
                beta_0[(ind1 + 1)] = beta_step_int[jk];
                tform[ind1] = "step_slope";
                tform[(ind1 + 1)] = "step_ind";
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
        } else {//past this point, non-dose terms are calcualted
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
            T0.col(ind1) = (df0.col(ind1).array() * beta_0[ind1]).matrix();//a similar form is used for every form
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
//    T0.block(0,0,T0.rows(),dose_num_tot) = De.block(0,0,T0.rows(),dose_num_tot);
//    Td0.block(0,0,T0.rows(),dose_num_tot) = Dde.block(0,0,T0.rows(),dose_num_tot);
    //
    #pragma omp parallel for schedule(dynamic) num_threads(nthreads)
    for (int ijk=0;ijk<dose_num_tot;ijk++){//inserts the dose terms into the total term arrays
        Tdd0.col(ijk) = Ddde.col(ijk*(ijk+1)/2+ijk);
        Td0.col(ijk) = Dde.col(ijk);
        T0.col(ijk) = De.col(ijk);
    }
    // ---------------------------------------------------------
    // Prints off a series of calaculations to check at what point values are changing
    // ---------------------------------------------------------
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
//    cout << "doses checked ";
//    for (int ijk=0;ijk<totalnum;ijk++){
//        cout << df0.col(ijk).sum() << " ";
//    }
//    cout << " " << endl;
    cout << "derivs checked ";
    for (int ijk=0;ijk<totalnum;ijk++){
        cout << Td0.col(ijk).sum() << " ";
    }
    cout << " " << endl;
//    cout << "forms checked ";
//    for (int ijk=0;ijk<totalnum;ijk++){
//        cout << tform[ijk] << " ";
//    }
//    cout << " " << endl;
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
                            Rdd.block(0,ijk,df0.rows(),1) = Tdd0.block(0,ij,df0.rows(),1).array() * Te.block(0,0,df0.rows(),1).array();
                        } else {
                            Rdd.block(0,ijk,df0.rows(),1) = Tdd0.block(0,ij,df0.rows(),1).array() * T0.block(0,fir,df0.rows(),1).array();
                        }
                    } else {
                        if (ij<dose_num_tot){
                            Rdd.block(0,ijk,df0.rows(),1) = Ddde.block(0,ijk,df0.rows(),1).array() * Te.block(0,0,df0.rows(),1).array();
                        } else {
                            Rdd.block(0,ijk,df0.rows(),1) = Tdd0.block(0,ij,df0.rows(),1).array() * Dose.block(0,0,df0.rows(),1).array();
                        }
                    }
                } else {
                    if (fir>=dose_num_tot){
                        if ((ij==fir)||(jk==fir)){
                            Rdd.block(0,ijk,df0.rows(),1) = Td0.block(0,ij,df0.rows(),1).array() * Td0.block(0,jk,df0.rows(),1).array();
                        } else if ((ij<dose_num_tot)&&(jk<dose_num_tot)){
                            Rdd.block(0,ijk,df0.rows(),1) = Ddde.block(0,ijk,df0.rows(),1).array() * T0.block(0,fir,df0.rows(),1).array();
                        }
                    } else {
                        if ((ij<dose_num_tot)&&(jk<dose_num_tot)){
                            Rdd.block(0,ijk,df0.rows(),1) = Te.block(0,0,df0.rows(),1).array() * Ddde.block(0,ijk,df0.rows(),1).array();
                        } else if ((ij>=dose_num_tot)&&(jk<dose_num_tot)){
                            Rdd.block(0,ijk,df0.rows(),1) = Td0.block(0,ij,df0.rows(),1).array() * Td0.block(0,jk,df0.rows(),1).array();
                        }
                    }
                }
            }
        }
        RdR << R.rowwise().replicate(totalnum).array().pow(-1).array() * Rd.array();
        RddR << R.rowwise().replicate(totalnum*(totalnum+1)/2).array().pow(-1).array() * Rdd.array();
    }else if (modelform=="M"){
        Te = Te.array() * 0 + 1; //verifies the initial term product is 1
//        cout << "Te checked " << Te.array().sum() << " " << T0.array().block(0,dose_num_tot,T0.rows(),T0.cols()-dose_num_tot).rowwise().prod().sum() << " " << Dose.array().sum() << " " << T0.array().rowwise().prod().sum() << endl;
        //
        Te = T0.array().block(0,dose_num_tot,T0.rows(),T0.cols()-dose_num_tot).rowwise().prod() * Dose.array();
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
                    Rdd.block(0,ijk,df0.rows(),1) = Ddde.block(0,ijk,df0.rows(),1).array() * Te.block(0,0,df0.rows(),1).array();
                } else {
                    Rdd.block(0,ijk,df0.rows(),1) = Tdd0.block(0,jk,df0.rows(),1).array() * T0.block(0,jk,df0.rows(),1).array().pow(-1).array() * R.block(0,0,df0.rows(),1).array();
                }
            } else {
                if ((ij<dose_num_tot)&&(jk<dose_num_tot)){
                    Rdd.block(0,ijk,df0.rows(),1) = Ddde.block(0,ijk,df0.rows(),1).array() * Te.block(0,0,df0.rows(),1).array();
                } else if ((ij<dose_num_tot)||(jk<dose_num_tot)){
                    Rdd.block(0,ijk,df0.rows(),1) = Dde.block(0,jk,df0.rows(),1).array() * De.block(0,jk,df0.rows(),1).array().pow(-1).array() * Rd.block(0,ij,df0.rows(),1).array();
                } else{
                    Rdd.block(0,ijk,df0.rows(),1) = Td0.block(0,jk,df0.rows(),1).array() * T0.block(0,jk,df0.rows(),1).array().pow(-1).array() * Rd.block(0,ij,df0.rows(),1).array();
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
    //
    List temp_list = List::create(_["betas"]=wrap(beta_0),_["Status"]="FAILED"); //used as a dummy return value for code checking
    if (R.minCoeff()<0){
        cout << "A negative risk was detected" << endl;
        return temp_list;
    }
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
//    cout << "risk2 checked ";
//    for (int ijk=0;ijk<totalnum;ijk++){
//        cout << Rdd.col(ijk*(ijk+1)/2+ijk).sum() << " ";
//    }
//    cout << " " << endl;
    //
    end_point = system_clock::now();
    ending = time_point_cast<microseconds>(end_point).time_since_epoch().count();
    cout<<"df100 "<<(ending-start)<<" "<<0<<" "<<0<<" "<<-1<<",Prep_R"<<endl;
    gibtime = system_clock::to_time_t(system_clock::now());
    cout << ctime(&gibtime) << endl;
    //
    // -------------------------------------------------------------------------------------------
    //
    vector<string>  RiskGroup(ntime); //vector of strings detailing the rows
    IntegerMatrix RiskFail(ntime,2); //vector giving the event rows
    //
    const Map<MatrixXd> df_m(as<Map<MatrixXd> >(df_groups));
    #pragma omp parallel for schedule(dynamic) num_threads(nthreads)
    for (int ijk=0;ijk<ntime;ijk++){
        double t0 = tu[ijk];
        VectorXi select_ind_all = ((df_m.col(0).array() < t0)&&(df_m.col(1).array()>=t0)).cast<int>(); //indices at risk
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
    //
    //
//    cout << "riskr checked ";
//    for (int ijk=0;ijk<totalnum;ijk++){
//        cout << Rls1.col(0).sum() << " ";
//    }
//    cout << " " << endl;
//    cout << "risk1r checked ";
//    for (int ijk=0;ijk<totalnum;ijk++){
//        cout << Rls2.col(ijk).sum() << " ";
//    }
//    cout << " " << endl;
//    cout << "risk2r checked ";
//    for (int ijk=0;ijk<totalnum;ijk++){
//        cout << Rls3.col(ijk*(ijk+1)/2+ijk).sum() << " ";
//    }
//    cout << " " << endl;
//    //
//    //
//    cout << "riskl checked ";
//    for (int ijk=0;ijk<totalnum;ijk++){
//        cout << Lls1.col(0).sum() << " ";
//    }
//    cout << " " << endl;
//    cout << "risk1l checked ";
//    for (int ijk=0;ijk<totalnum;ijk++){
//        cout << Lls2.col(ijk).sum() << " ";
//    }
//    cout << " " << endl;
//    cout << "risk2l checked ";
//    for (int ijk=0;ijk<totalnum;ijk++){
//        cout << Lls3.col(ijk*(ijk+1)/2+ijk).sum() << " ";
//    }
//    cout << " " << endl;
    //
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
            temp1 = Ld.col(3).array() - temp1.array() * temp2.array();
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
    //
    vector <double> Ll_comp(2,Ll[0]); //vector to compare values
    //
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
        VectorXd::Map(&beta_p[0], beta_0.size()) = beta_0;//
        VectorXd::Map(&beta_c[0], beta_0.size()) = beta_0;//
        VectorXd::Map(&beta_a[0], beta_0.size()) = beta_0;//
        VectorXd::Map(&beta_best[0], beta_0.size()) = beta_0;//
        //
        #pragma omp parallel for schedule(dynamic) num_threads(nthreads)
        for (int ijk=0;ijk<totalnum;ijk++){
            if (change_all){
                dbeta[ijk] = -lr * Lld[ijk] / Lldd[ijk*totalnum+ijk];
                //
                double dbeta_max = abs(Ll[ijk]/Lld[ijk] * dbeta_cap);//uses newtonian step for zero log-likelihood as a limit
                if (abs(dbeta[ijk])>dbeta_max){
                    dbeta[ijk] = dbeta_max * sign(dbeta[ijk]);
                }
                if ((tform[ijk]=="step_ind")||(tform[ijk]=="lin_ind")){ //the threshold values use different maximum deviation values
                    if (abs(dbeta[ijk])>dose_abs_max){
                        dbeta[ijk] = dose_abs_max * sign(dbeta[ijk]);
                    }
                }else{
                    if (abs(dbeta[ijk])>abs_max){
                        dbeta[ijk] = abs_max * sign(dbeta[ijk]);
                    }
                }
            }else{
                if (ijk!=fir){//Validation requires controlled changes
                    dbeta[ijk] = 0.0;
                } else {
                    if ((tform[ijk]=="step_ind")||(tform[ijk]=="lin_ind")){
                        dbeta[ijk] = dose_abs_max/5.0;
                    } else {
                        dbeta[ijk] = 0.01;
                    }
                }
            }
        }
        cout << "testing dbeta: ";//prints the final changes for validation
        for (int ijk=0;ijk<totalnum;ijk++){
            cout << dbeta[ijk] << " ";
        }
        cout << " " << endl;
        //
        Ll_best = Ll[ind0];
        i = ind0;
        //
        halves=0;
        while ((Ll[ind0] <= Ll_best)&&(halves<halfmax)){ //repeats until half-steps maxed or an improvement
            halves++;
            if (change_all){ //calculates the terms for the new parameters
                Dose = VectorXd::Zero(df_dose.rows());
                #pragma omp parallel for schedule(dynamic) num_threads(nthreads) reduction(eig_plus:Dose)
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
                    } else if (tform[ijk]=="loglin_slope"){//the existence of this term form implies that the slope exists
                        beta_c[ijk] = beta_a[ijk] + dbeta[ijk];
                        beta_c[ijk+1] = beta_a[ijk+1] + dbeta[ijk+1];//slope and top are changed at the same time
                        double ach = beta_c[ijk]/beta_p[ijk];//change due to slope
                        MatrixXd bch = ((beta_c[ijk+1] - beta_p[ijk+1]) * df0.col(ijk)).array().exp().array();//change due to exponential
                        //
                        De.col(ijk) =   ach * bch.array() * De.col(ijk).array();
                        De.col(ijk+1) = ach * bch.array() * De.col(ijk+1).array();
                        Dde.col(ijk) = bch.array() * Dde.col(ijk).array();
                        Dde.col(ijk+1) = ach * bch.array() * Dde.col(ijk+1).array();
                        Ddde.col((ijk+1)*(ijk+2)/2+ijk) = bch.array() * Ddde.col((ijk+1)*(ijk+2)/2+ijk).array();
                        Ddde.col((ijk+1)*(ijk+2)/2+ijk+1) = ach * bch.array() * Ddde.col((ijk+1)*(ijk+2)/2+ijk+1).array();
                        Dose = Dose.array() + De.col(ijk).array();
                        //
                    } else if (tform[ijk]=="loglin_top"){//The change is only made if a slope value isn't changed
                        if (ijk==0){
                            beta_c[ijk] = beta_a[ijk] + dbeta[ijk];
                            MatrixXd bch = ((beta_c[ijk] - beta_p[ijk]) * df0.col(ijk)).array().exp().array();
                            //
                            De.col(ijk) = De.col(ijk).array() * bch.array();
                            Dde.col(ijk) = Dde.col(ijk).array() * bch.array();
                            Ddde.col((ijk)*(ijk+1)/2+ijk) = Ddde.col((ijk)*(ijk+1)/2+ijk).array() * bch.array();
                            Dose = Dose.array() + De.col(ijk).array();
                            //
                        } else if (tform[ijk-1]!="loglin_slope"){
                            beta_c[ijk] = beta_a[ijk] + dbeta[ijk];
                            MatrixXd bch = ((beta_c[ijk] - beta_p[ijk]) * df0.col(ijk)).array().exp().array();
                            //
                            De.col(ijk) = De.col(ijk).array() * bch.array();
                            Dde.col(ijk) = Dde.col(ijk).array() * bch.array();
                            Ddde.col((ijk)*(ijk+1)/2+ijk) = Ddde.col((ijk)*(ijk+1)/2+ijk).array() * bch.array();
                            Dose = Dose.array() + De.col(ijk).array();
                        } else {
                            ;
                        }
                    } else if (tform[ijk]=="lin_slope"){//threshold terms are recalculted
                        beta_c[ijk] = beta_a[ijk] + dbeta[ijk];
                        beta_c[ijk+1] = beta_a[ijk+1] + dbeta[ijk+1];
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
                        beta_c[ijk] = beta_a[ijk] + dbeta[ijk];
                        De.col(ijk) = (beta_c[ijk]/beta_p[ijk]) * De.col(ijk).array();
                        Dose = Dose.array() + De.col(ijk).array();
                        //
                    } else if (tform[ijk]=="step_slope") {
                        beta_c[ijk] = beta_a[ijk] + dbeta[ijk];
                        beta_c[ijk+1] = beta_a[ijk+1] + dbeta[ijk+1];
                        
                        ArrayXd temp = (df0.col(ijk).array() - beta_c[ijk+1]);
                        ArrayXd temp0 = (df0.col(ijk).array() - beta_c[ijk+1]-dint);
                        ArrayXd temp1 = (df0.col(ijk).array() - beta_c[ijk+1]+dint);
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
                        cout << ijk << " not changed, with form " << tform[ijk] << endl;//prints which aren't changed, should only be the threshold parameters
                    }
                }
            } else {
                #pragma omp parallel for schedule(dynamic) num_threads(nthreads) reduction(eig_plus:Dose)
                for (int ijk=fir;ijk<fir+1;ijk++){//If only one is changed, then the identification is done differently
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
                        double ach = beta_c[ijk]/beta_p[ijk];
                        //
                        De.col(ijk) =   ach * De.col(ijk).array();
                        De.col(ijk+1) = ach * De.col(ijk+1).array();
                        Dde.col(ijk+1) = ach * Dde.col(ijk+1).array();
                        Ddde.col((ijk+1)*(ijk+2)/2+ijk+1) = ach * Ddde.col((ijk+1)*(ijk+2)/2+ijk+1).array();
                        Dose = Dose.array() + De.col(ijk).array();
                        //
                    } else if (tform[ijk]=="loglin_top"){
                        if (ijk==0){
                            beta_c[ijk] = beta_a[ijk] + dbeta[ijk];
                            MatrixXd bch = ((beta_c[ijk] - beta_p[ijk]) * df0.col(ijk)).array().exp().array();
                            //
                            De.col(ijk) = De.col(ijk).array() * bch.array();
                            Dde.col(ijk) = Dde.col(ijk).array() * bch.array();
                            Ddde.col((ijk)*(ijk+1)/2+ijk) = Ddde.col((ijk)*(ijk+1)/2+ijk).array() * bch.array();
                            Dose = Dose.array() + De.col(ijk).array();
                            //
                        } else if (tform[ijk-1]!="loglin_slope"){
                            beta_c[ijk] = beta_a[ijk] + dbeta[ijk];
                            MatrixXd bch = ((beta_c[ijk] - beta_p[ijk]) * df0.col(ijk)).array().exp().array();
                            //
                            De.col(ijk) = De.col(ijk).array() * bch.array();
                            Dde.col(ijk) = Dde.col(ijk).array() * bch.array();
                            Ddde.col((ijk)*(ijk+1)/2+ijk) = Ddde.col((ijk)*(ijk+1)/2+ijk).array() * bch.array();
                            Dose = Dose.array() + De.col(ijk).array();
                        } else {
                            ;
                        }
                    } else if (tform[ijk]=="lin_slope"){
                        beta_c[ijk] = beta_a[ijk] + dbeta[ijk];
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
                    } else if (tform[ijk]=="lin_ind"){
                        beta_c[ijk] = beta_a[ijk] + dbeta[ijk];
                        //
                        ArrayXd temp = (df0.col(ijk-1).array() - beta_c[ijk]);
                        ArrayXd temp0 = (df0.col(ijk-1).array() - beta_c[ijk]-dint);
                        ArrayXd temp1 = (df0.col(ijk-1).array() - beta_c[ijk]+dint);
                        //
                        temp = (temp.array() < 0).select(0.0, temp);
                        temp0 = (temp0.array() < 0).select(0.0, temp0);
                        temp1 = (temp1.array() < 0).select(0.0, temp1);
                        //
                        De.col(ijk-1) = beta_c[ijk-1] * temp.array();
                        De.col((ijk)) = beta_c[ijk-1] * temp.array();
                        //
                        Dde.col(ijk-1) = temp.array();
                        Dde.col((ijk)) = beta_c[ijk] * (temp1.array()-temp0.array()) / 2.0/dint;
                        //
                        Ddde.col((ijk) * ((ijk)+1)/2 + ijk-1) = (temp1.array()-temp0.array()) / 2.0/dint;
                        Ddde.col((ijk) * ((ijk)+1)/2 + (ijk)) = beta_c[ijk-1] * (temp1.array()-2*temp.array()+temp0.array()) / pow(dint,2);
                        Dose = Dose.array() + De.col(ijk).array();
                        //
                    } else if (tform[ijk]=="quad_slope"){
                        beta_c[ijk] = beta_a[ijk] + dbeta[ijk];
                        De.col(ijk) = beta_c[ijk] * df0.col(ijk).array().square().array();
                        Dose = Dose.array() + De.col(ijk).array();
                        //
                    } else if (tform[ijk]=="step_slope") {
                        beta_c[ijk] = beta_a[ijk] + dbeta[ijk];
                        
                        ArrayXd temp = (df0.col(ijk).array() - beta_c[ijk+1]);
                        ArrayXd temp0 = (df0.col(ijk).array() - beta_c[ijk+1]-dint);
                        ArrayXd temp1 = (df0.col(ijk).array() - beta_c[ijk+1]+dint);
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
                    } else if (tform[ijk]=="step_ind") {
                        beta_c[ijk] = beta_a[ijk] + dbeta[ijk];
                        
                        ArrayXd temp = (df0.col(ijk-1).array() - beta_c[ijk]);
                        ArrayXd temp0 = (df0.col(ijk-1).array() - beta_c[ijk]-dint);
                        ArrayXd temp1 = (df0.col(ijk-1).array() - beta_c[ijk]+dint);
                        //
                        temp = (temp.array() < 0).select(0.0, MatrixXd::Zero(temp.rows(),temp.cols()).array()+1.0);
                        temp0 = (temp0.array() < 0).select(0.0, MatrixXd::Zero(temp.rows(),temp0.cols()).array()+1.0);
                        temp1 = (temp1.array() < 0).select(0.0, MatrixXd::Zero(temp.rows(),temp1.cols()).array()+1.0);
                        //
                        De.col(ijk-1) = beta_c[ijk-1] * temp.array();
                        De.col((ijk)) = beta_c[ijk-1] * temp.array();
                        //
                        Dde.col(ijk-1) = temp.array();
                        Dde.col((ijk)) = beta_c[ijk-1] * (temp1.array()-temp0.array()) / 2.0/dint;
                        //
                        Ddde.col((ijk) * ((ijk)+1)/2 + ijk-1) = (temp1.array()-temp0.array()) / 2.0/dint;
                        Ddde.col((ijk) * ((ijk)+1)/2 + (ijk)) = beta_c[ijk-1] * (temp1.array()-2*temp.array()+temp0.array()) / pow(dint,2);
                        Dose = Dose.array() + De.col(ijk).array();
                        //
                    } else {
                        cout << ijk << " not changed, with form " << tform[ijk] << endl;
                    }
                }
            }
            #pragma omp parallel for schedule(dynamic) num_threads(nthreads)
            for (int ijk=0;ijk<dose_num_tot;ijk++){//merges dose and all term matrices
                Tdd0.col(ijk) = Ddde.col(ijk*(ijk+1)/2+ijk);
                Td0.col(ijk) = Dde.col(ijk);
                T0.col(ijk) = De.col(ijk);
            }
//            cout << "values checked ";
//            for (int ijk=0;ijk<totalnum;ijk++){
//                cout << beta_0[ijk] << " ";
//            }
//            cout << " " << endl;
//            //
//            cout << "values_a checked ";
//            for (int ijk=0;ijk<totalnum;ijk++){
//                cout << beta_a[ijk] << " ";
//            }
//            cout << " " << endl;
//            cout << "values_p checked ";
//            for (int ijk=0;ijk<totalnum;ijk++){
//                cout << beta_p[ijk] << " ";
//            }
//            cout << " " << endl;
//            cout << "values_c checked ";
//            for (int ijk=0;ijk<totalnum;ijk++){
//                cout << beta_c[ijk] << " ";
//            }
//            cout << " " << endl;
//            //
//            cout << "sums checked ";
//            for (int ijk=0;ijk<totalnum;ijk++){
//                cout << T0.col(ijk).sum() << " ";
//            }
//            cout << " " << endl;
//            cout << "doses checked ";
//            for (int ijk=0;ijk<totalnum;ijk++){
//                cout << df0.col(ijk).sum() << " ";
//            }
//            cout << " " << endl;
//            cout << "derivs checked ";
//            for (int ijk=0;ijk<totalnum;ijk++){
//                cout << Td0.col(ijk).sum() << " ";
//            }
//            cout << " " << endl;
            //
//            Dose << De.rowwise().sum();
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
                                    Rdd.block(0,ijk,df0.rows(),1) = Tdd0.block(0,ij,df0.rows(),1).array() * Te.block(0,0,df0.rows(),1).array();
                                } else {
                                    Rdd.block(0,ijk,df0.rows(),1) = Tdd0.block(0,ij,df0.rows(),1).array() * T0.block(0,fir,df0.rows(),1).array();
                                }
                            } else {
                                if (ij<dose_num_tot){
                                    Rdd.block(0,ijk,df0.rows(),1) = Ddde.block(0,ijk,df0.rows(),1).array() * Te.block(0,0,df0.rows(),1).array();
                                } else {
                                    Rdd.block(0,ijk,df0.rows(),1) = Tdd0.block(0,ij,df0.rows(),1).array() * Dose.block(0,0,df0.rows(),1).array();
                                }
                            }
                        } else {
                            if (fir>=dose_num_tot){
                                if ((ij==fir)||(jk==fir)){
                                    Rdd.block(0,ijk,df0.rows(),1) = Td0.block(0,ij,df0.rows(),1).array() * Td0.block(0,jk,df0.rows(),1).array();
                                } else if ((ij<dose_num_tot)&&(jk<dose_num_tot)){
                                    Rdd.block(0,ijk,df0.rows(),1) = Ddde.block(0,ijk,df0.rows(),1).array() * T0.block(0,fir,df0.rows(),1).array();
                                }
                            } else {
                                if ((ij<dose_num_tot)&&(jk<dose_num_tot)){
                                    Rdd.block(0,ijk,df0.rows(),1) = Te.block(0,0,df0.rows(),1).array() * Ddde.block(0,ijk,df0.rows(),1).array();
                                } else if ((ij>=dose_num_tot)&&(jk<dose_num_tot)){
                                    Rdd.block(0,ijk,df0.rows(),1) = Td0.block(0,ij,df0.rows(),1).array() * Td0.block(0,jk,df0.rows(),1).array();
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
//                cout << "Te checked " << Te.array().sum() << " " << T0.array().block(0,dose_num_tot,T0.rows(),T0.cols()-dose_num_tot).rowwise().prod().sum() << " " << Dose.array().sum() << " " << T0.array().rowwise().prod().sum() << endl;
                Te = T0.array().block(0,dose_num_tot,T0.rows(),T0.cols()-dose_num_tot).rowwise().prod() * Dose.array();
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
                            Rdd.block(0,ijk,df0.rows(),1) = Ddde.block(0,ijk,df0.rows(),1).array() * Te.block(0,0,df0.rows(),1).array();
                        } else {
                            Rdd.block(0,ijk,df0.rows(),1) = Tdd0.block(0,jk,df0.rows(),1).array() * T0.block(0,jk,df0.rows(),1).array().pow(-1).array() * R.block(0,0,df0.rows(),1).array();
                        }
                    } else {
                        if ((ij<dose_num_tot)&&(jk<dose_num_tot)){
                            Rdd.block(0,ijk,df0.rows(),1) = Ddde.block(0,ijk,df0.rows(),1).array() * Te.block(0,0,df0.rows(),1).array();
                        } else if ((ij<dose_num_tot)||(jk<dose_num_tot)){
                            Rdd.block(0,ijk,df0.rows(),1) = Dde.block(0,jk,df0.rows(),1).array() * De.block(0,jk,df0.rows(),1).array().pow(-1).array() * Rd.block(0,ij,df0.rows(),1).array();
                        } else{
                            Rdd.block(0,ijk,df0.rows(),1) = Td0.block(0,jk,df0.rows(),1).array() * T0.block(0,jk,df0.rows(),1).array().pow(-1).array() * Rd.block(0,ij,df0.rows(),1).array();
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
            //
            temp_list = List::create(_["betas"]=wrap(beta_0),_["Status"]="FAILED"); //used as a dummy return value for code checking
            if (R.minCoeff()<0){
                cout << "A negative risk was detected" << endl;
                return temp_list;
            }
//            cout << "risk checked ";
//            for (int ijk=0;ijk<totalnum;ijk++){
//                cout << R.col(0).sum() << " ";
//            }
//            cout << " " << endl;
//            cout << "risk1 checked ";
//            for (int ijk=0;ijk<totalnum;ijk++){
//                cout << Rd.col(ijk).sum() << " ";
//            }
//            cout << " " << endl;
//            cout << "risk2 checked ";
//            for (int ijk=0;ijk<totalnum;ijk++){
//                cout << Rdd.col(ijk*(ijk+1)/2+ijk).sum() << " ";
//            }
//            cout << " " << endl;
            end_point = system_clock::now();
            ending = time_point_cast<microseconds>(end_point).time_since_epoch().count();
            cout<<"df100 "<<(ending-start)<<" "<<halves<<" "<<iteration<<" "<<ind0<<",Update_R"<<endl;
            //
            gibtime = system_clock::to_time_t(system_clock::now());
            cout << ctime(&gibtime) << endl;
            fill(Ll.begin(), Ll.end(), 0.0);
            fill(Lld.begin(), Lld.end(), 0.0);
            fill(Lldd.begin(), Lldd.end(), 0.0);
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
                    }
                    MatrixXd Ld = MatrixXd::Zero(dj,4);
                    Ld << R.block(RiskFail(j,0),0,dj,1), Rd.block(RiskFail(j,0),ij,dj,1), Rd.block(RiskFail(j,0),jk,dj,1) ,Rdd.block(RiskFail(j,0),ijk,dj,1);//sum of risks in group
                    //
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
            //
            #pragma omp declare reduction(vec_double_plus : std::vector<double> : \
                std::transform(omp_out.begin(), omp_out.end(), omp_in.begin(), omp_out.begin(), std::plus<double>())) \
                initializer(omp_priv = omp_orig)
            #pragma omp parallel for schedule(dynamic) num_threads(nthreads) reduction(vec_double_plus:Ll,Lld,Lldd) collapse(2)
            for (int ijk=0;ijk<totalnum*(totalnum+1)/2;ijk++){//totalnum*(totalnum+1)/2
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
                    Ld << R.block(RiskFail(j,0),0,dj,1), RdR.block(RiskFail(j,0),ij,dj,1), RdR.block(RiskFail(j,0),jk,dj,1) ,RddR.block(RiskFail(j,0),ijk,dj,1);//sum of risks in group
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
                    MatrixXd temp1 = MatrixXd::Zero(Ld.rows(),1);
                    MatrixXd temp2 = MatrixXd::Zero(Ld.rows(),1);
                    temp1 = Ld.col(0).array().log();
                    double Ld1 =  (temp1.array().isFinite()).select(temp1,0).sum();
                    temp1 = Ld.col(1).array();
                    temp2 = Ld.col(2).array();
                    double Ld2 = (temp1.array().isFinite()).select(temp1,0).sum();
                    temp1 = Ld.col(3).array() - temp1.array() * temp2.array();
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
            #pragma omp parallel for num_threads(nthreads)
            for (int ijk=0;ijk<totalnum;ijk++){//updates the previous point
                beta_p[ijk] = beta_c[ijk];
            }
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
            #pragma omp parallel for num_threads(nthreads)
            for (int ijk=0;ijk<totalnum;ijk++){//totalnum*(totalnum+1)/2
                beta_0[ijk] = beta_best[ijk];
            }
        }
        if (beta_best!=beta_p){//if the risk matrices aren't the optimal values, then they must be recalculated
            Dose = VectorXd::Zero(df_dose.rows());
            if (change_all){ //calculates the terms for the new parameters
                #pragma omp parallel for schedule(dynamic) num_threads(nthreads) reduction(eig_plus:Dose)
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
                    } else if (tform[ijk]=="loglin_slope"){//the existence of this term form implies that the slope exists
                        beta_c[ijk] = beta_a[ijk] + dbeta[ijk];
                        beta_c[ijk+1] = beta_a[ijk+1] + dbeta[ijk+1];//slope and top are changed at the same time
                        double ach = beta_c[ijk]/beta_p[ijk];//change due to slope
                        MatrixXd bch = ((beta_c[ijk+1] - beta_p[ijk+1]) * df0.col(ijk)).array().exp().array();//change due to exponential
                        //
                        De.col(ijk) =   ach * bch.array() * De.col(ijk).array();
                        De.col(ijk+1) = ach * bch.array() * De.col(ijk+1).array();
                        Dde.col(ijk) = bch.array() * Dde.col(ijk).array();
                        Dde.col(ijk+1) = ach * bch.array() * Dde.col(ijk+1).array();
                        Ddde.col((ijk+1)*(ijk+2)/2+ijk) = bch.array() * Ddde.col((ijk+1)*(ijk+2)/2+ijk).array();
                        Ddde.col((ijk+1)*(ijk+2)/2+ijk+1) = ach * bch.array() * Ddde.col((ijk+1)*(ijk+2)/2+ijk+1).array();
                        Dose = Dose.array() + De.col(ijk).array();
                        //
                    } else if (tform[ijk]=="loglin_top"){//The change is only made if a slope value isn't changed
                        if (ijk==0){
                            beta_c[ijk] = beta_a[ijk] + dbeta[ijk];
                            MatrixXd bch = ((beta_c[ijk] - beta_p[ijk]) * df0.col(ijk)).array().exp().array();
                            //
                            De.col(ijk) = De.col(ijk).array() * bch.array();
                            Dde.col(ijk) = Dde.col(ijk).array() * bch.array();
                            Ddde.col((ijk)*(ijk+1)/2+ijk) = Ddde.col((ijk)*(ijk+1)/2+ijk).array() * bch.array();
                            Dose = Dose.array() + De.col(ijk).array();
                            //
                        } else if (tform[ijk-1]!="loglin_slope"){
                            beta_c[ijk] = beta_a[ijk] + dbeta[ijk];
                            MatrixXd bch = ((beta_c[ijk] - beta_p[ijk]) * df0.col(ijk)).array().exp().array();
                            //
                            De.col(ijk) = De.col(ijk).array() * bch.array();
                            Dde.col(ijk) = Dde.col(ijk).array() * bch.array();
                            Ddde.col((ijk)*(ijk+1)/2+ijk) = Ddde.col((ijk)*(ijk+1)/2+ijk).array() * bch.array();
                            Dose = Dose.array() + De.col(ijk).array();
                        } else {
                            ;
                        }
                    } else if (tform[ijk]=="lin_slope"){//threshold terms are recalculted
                        beta_c[ijk] = beta_a[ijk] + dbeta[ijk];
                        beta_c[ijk+1] = beta_a[ijk+1] + dbeta[ijk+1];
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
                        beta_c[ijk] = beta_a[ijk] + dbeta[ijk];
                        De.col(ijk) = (beta_c[ijk]/beta_p[ijk]) * De.col(ijk).array();
                        Dose = Dose.array() + De.col(ijk).array();
                        //
                    } else if (tform[ijk]=="step_slope") {
                        beta_c[ijk] = beta_a[ijk] + dbeta[ijk];
                        beta_c[ijk+1] = beta_a[ijk+1] + dbeta[ijk+1];
                        
                        ArrayXd temp = (df0.col(ijk).array() - beta_c[ijk+1]);
                        ArrayXd temp0 = (df0.col(ijk).array() - beta_c[ijk+1]-dint);
                        ArrayXd temp1 = (df0.col(ijk).array() - beta_c[ijk+1]+dint);
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
                        cout << ijk << " not changed, with form " << tform[ijk] << endl;//prints which aren't changed, should only be the threshold parameters
                    }
                }
            } else {
                #pragma omp parallel for schedule(dynamic) num_threads(nthreads) reduction(eig_plus:Dose)
                for (int ijk=fir;ijk<fir+1;ijk++){//If only one is changed, then the identification is done differently

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
                        double ach = beta_c[ijk]/beta_p[ijk];
                        //
                        De.col(ijk) =   ach * De.col(ijk).array();
                        De.col(ijk+1) = ach * De.col(ijk+1).array();
                        Dde.col(ijk+1) = ach * Dde.col(ijk+1).array();
                        Ddde.col((ijk+1)*(ijk+2)/2+ijk+1) = ach * Ddde.col((ijk+1)*(ijk+2)/2+ijk+1).array();
                        Dose = Dose.array() + De.col(ijk).array();
                        //
                    } else if (tform[ijk]=="loglin_top"){
                        if (ijk==0){
                            beta_c[ijk] = beta_a[ijk] + dbeta[ijk];
                            MatrixXd bch = ((beta_c[ijk] - beta_p[ijk]) * df0.col(ijk)).array().exp().array();
                            //
                            De.col(ijk) = De.col(ijk).array() * bch.array();
                            Dde.col(ijk) = Dde.col(ijk).array() * bch.array();
                            Ddde.col((ijk)*(ijk+1)/2+ijk) = Ddde.col((ijk)*(ijk+1)/2+ijk).array() * bch.array();
                            Dose = Dose.array() + De.col(ijk).array();
                            //
                        } else if (tform[ijk-1]!="loglin_slope"){
                            beta_c[ijk] = beta_a[ijk] + dbeta[ijk];
                            MatrixXd bch = ((beta_c[ijk] - beta_p[ijk]) * df0.col(ijk)).array().exp().array();
                            //
                            De.col(ijk) = De.col(ijk).array() * bch.array();
                            Dde.col(ijk) = Dde.col(ijk).array() * bch.array();
                            Ddde.col((ijk)*(ijk+1)/2+ijk) = Ddde.col((ijk)*(ijk+1)/2+ijk).array() * bch.array();
                            Dose = Dose.array() + De.col(ijk).array();
                        } else {
                            ;
                        }
                    } else if (tform[ijk]=="lin_slope"){
                        beta_c[ijk] = beta_a[ijk] + dbeta[ijk];
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
                    } else if (tform[ijk]=="lin_ind"){
                        beta_c[ijk] = beta_a[ijk] + dbeta[ijk];
                        //
                        ArrayXd temp = (df0.col(ijk-1).array() - beta_c[ijk]);
                        ArrayXd temp0 = (df0.col(ijk-1).array() - beta_c[ijk]-dint);
                        ArrayXd temp1 = (df0.col(ijk-1).array() - beta_c[ijk]+dint);
                        //
                        temp = (temp.array() < 0).select(0.0, temp);
                        temp0 = (temp0.array() < 0).select(0.0, temp0);
                        temp1 = (temp1.array() < 0).select(0.0, temp1);
                        //
                        De.col(ijk-1) = beta_c[ijk-1] * temp.array();
                        De.col((ijk)) = beta_c[ijk-1] * temp.array();
                        //
                        Dde.col(ijk-1) = temp.array();
                        Dde.col((ijk)) = beta_c[ijk] * (temp1.array()-temp0.array()) / 2.0/dint;
                        //
                        Ddde.col((ijk) * ((ijk)+1)/2 + ijk-1) = (temp1.array()-temp0.array()) / 2.0/dint;
                        Ddde.col((ijk) * ((ijk)+1)/2 + (ijk)) = beta_c[ijk-1] * (temp1.array()-2*temp.array()+temp0.array()) / pow(dint,2);
                        Dose = Dose.array() + De.col(ijk).array();
                        //
                    } else if (tform[ijk]=="quad_slope"){
                        beta_c[ijk] = beta_a[ijk] + dbeta[ijk];
                        De.col(ijk) = beta_c[ijk] * df0.col(ijk).array().square().array();
                        Dose = Dose.array() + De.col(ijk).array();
                        //
                    } else if (tform[ijk]=="step_slope") {
                        beta_c[ijk] = beta_a[ijk] + dbeta[ijk];
                        
                        ArrayXd temp = (df0.col(ijk).array() - beta_c[ijk+1]);
                        ArrayXd temp0 = (df0.col(ijk).array() - beta_c[ijk+1]-dint);
                        ArrayXd temp1 = (df0.col(ijk).array() - beta_c[ijk+1]+dint);
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
                    } else if (tform[ijk]=="step_ind") {
                        beta_c[ijk] = beta_a[ijk] + dbeta[ijk];
                        
                        ArrayXd temp = (df0.col(ijk-1).array() - beta_c[ijk]);
                        ArrayXd temp0 = (df0.col(ijk-1).array() - beta_c[ijk]-dint);
                        ArrayXd temp1 = (df0.col(ijk-1).array() - beta_c[ijk]+dint);
                        //
                        temp = (temp.array() < 0).select(0.0, MatrixXd::Zero(temp.rows(),temp.cols()).array()+1.0);
                        temp0 = (temp0.array() < 0).select(0.0, MatrixXd::Zero(temp.rows(),temp0.cols()).array()+1.0);
                        temp1 = (temp1.array() < 0).select(0.0, MatrixXd::Zero(temp.rows(),temp1.cols()).array()+1.0);
                        //
                        De.col(ijk-1) = beta_c[ijk-1] * temp.array();
                        De.col((ijk)) = beta_c[ijk-1] * temp.array();
                        //
                        Dde.col(ijk-1) = temp.array();
                        Dde.col((ijk)) = beta_c[ijk-1] * (temp1.array()-temp0.array()) / 2.0/dint;
                        //
                        Ddde.col((ijk) * ((ijk)+1)/2 + ijk-1) = (temp1.array()-temp0.array()) / 2.0/dint;
                        Ddde.col((ijk) * ((ijk)+1)/2 + (ijk)) = beta_c[ijk-1] * (temp1.array()-2*temp.array()+temp0.array()) / pow(dint,2);
                        Dose = Dose.array() + De.col(ijk).array();
                        //
                    } else {
                        cout << ijk << " not changed, with form " << tform[ijk] << endl;
                    }
                }
            }
            #pragma omp parallel for schedule(dynamic) num_threads(nthreads)
            for (int ijk=0;ijk<dose_num_tot;ijk++){
                Tdd0.col(ijk) = Ddde.col(ijk*(ijk+1)/2+ijk);
                Td0.col(ijk) = Dde.col(ijk);
                T0.col(ijk) = De.col(ijk);
            }
//            cout << "values checked ";
//            for (int ijk=0;ijk<totalnum;ijk++){
//                cout << beta_0[ijk] << " ";
//            }
//            cout << " " << endl;
//            cout << "sums checked ";
//            for (int ijk=0;ijk<totalnum;ijk++){
//                cout << T0.col(ijk).sum() << " ";
//            }
//            cout << " " << endl;
//            cout << "doses checked ";
//            for (int ijk=0;ijk<totalnum;ijk++){
//                cout << df0.col(ijk).sum() << " ";
//            }
//            cout << " " << endl;
//            cout << "derivs checked ";
//            for (int ijk=0;ijk<totalnum;ijk++){
//                cout << Td0.col(ijk).sum() << " ";
//            }
//            cout << " " << endl;
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
                                    Rdd.block(0,ijk,df0.rows(),1) = Tdd0.block(0,ij,df0.rows(),1).array() * Te.block(0,0,df0.rows(),1).array();
                                } else {
                                    Rdd.block(0,ijk,df0.rows(),1) = Tdd0.block(0,ij,df0.rows(),1).array() * T0.block(0,fir,df0.rows(),1).array();
                                }
                            } else {
                                if (ij<dose_num_tot){
                                    Rdd.block(0,ijk,df0.rows(),1) = Ddde.block(0,ijk,df0.rows(),1).array() * Te.block(0,0,df0.rows(),1).array();
                                } else {
                                    Rdd.block(0,ijk,df0.rows(),1) = Tdd0.block(0,ij,df0.rows(),1).array() * Dose.block(0,0,df0.rows(),1).array();
                                }
                            }
                        } else {
                            if (fir>=dose_num_tot){
                                if ((ij==fir)||(jk==fir)){
                                    Rdd.block(0,ijk,df0.rows(),1) = Td0.block(0,ij,df0.rows(),1).array() * Td0.block(0,jk,df0.rows(),1).array();
                                } else if ((ij<dose_num_tot)&&(jk<dose_num_tot)){
                                    Rdd.block(0,ijk,df0.rows(),1) = Ddde.block(0,ijk,df0.rows(),1).array() * T0.block(0,fir,df0.rows(),1).array();
                                }
                            } else {
                                if ((ij<dose_num_tot)&&(jk<dose_num_tot)){
                                    Rdd.block(0,ijk,df0.rows(),1) = Te.block(0,0,df0.rows(),1).array() * Ddde.block(0,ijk,df0.rows(),1).array();
                                } else if ((ij>=dose_num_tot)&&(jk<dose_num_tot)){
                                    Rdd.block(0,ijk,df0.rows(),1) = Td0.block(0,ij,df0.rows(),1).array() * Td0.block(0,jk,df0.rows(),1).array();
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
                Te = T0.array().block(0,dose_num_tot,T0.rows(),T0.cols()-dose_num_tot).rowwise().prod() * Dose.array();
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
                            Rdd.block(0,ijk,df0.rows(),1) = Ddde.block(0,ijk,df0.rows(),1).array() * Te.block(0,0,df0.rows(),1).array();
                        } else {
                            Rdd.block(0,ijk,df0.rows(),1) = Tdd0.block(0,jk,df0.rows(),1).array() * T0.block(0,jk,df0.rows(),1).array().pow(-1).array() * R.block(0,0,df0.rows(),1).array();
                        }
                    } else {
                        if ((ij<dose_num_tot)&&(jk<dose_num_tot)){
                            Rdd.block(0,ijk,df0.rows(),1) = Ddde.block(0,ijk,df0.rows(),1).array() * Te.block(0,0,df0.rows(),1).array();
                        } else if ((ij<dose_num_tot)||(jk<dose_num_tot)){
                            Rdd.block(0,ijk,df0.rows(),1) = Dde.block(0,jk,df0.rows(),1).array() * De.block(0,jk,df0.rows(),1).array().pow(-1).array() * Rd.block(0,ij,df0.rows(),1).array();
                        } else{
                            Rdd.block(0,ijk,df0.rows(),1) = Td0.block(0,jk,df0.rows(),1).array() * T0.block(0,jk,df0.rows(),1).array().pow(-1).array() * Rd.block(0,ij,df0.rows(),1).array();
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
            fill(Ll.begin(), Ll.end(), 0.0);
            fill(Lld.begin(), Lld.end(), 0.0);
            fill(Lldd.begin(), Lldd.end(), 0.0);
            end_point = system_clock::now();
            ending = time_point_cast<microseconds>(end_point).time_since_epoch().count();
            cout<<"df100 "<<(ending-start)<<" "<<halves<<" "<<iteration<<" "<<ind0<<",Revert"<<endl;
            gibtime = system_clock::to_time_t(system_clock::now());
            cout << ctime(&gibtime) << endl;
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
                    }
                    MatrixXd Ld = MatrixXd::Zero(dj,4);
                    Ld << R.block(RiskFail(j,0),0,dj,1), Rd.block(RiskFail(j,0),ij,dj,1), Rd.block(RiskFail(j,0),jk,dj,1) ,Rdd.block(RiskFail(j,0),ijk,dj,1);//sum of risks in group
                    //
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
            //
            #pragma omp declare reduction(vec_double_plus : std::vector<double> : \
                std::transform(omp_out.begin(), omp_out.end(), omp_in.begin(), omp_out.begin(), std::plus<double>())) \
                initializer(omp_priv = omp_orig)
            #pragma omp parallel for schedule(dynamic) num_threads(nthreads) reduction(vec_double_plus:Ll,Lld,Lldd) collapse(2)
            for (int ijk=0;ijk<totalnum*(totalnum+1)/2;ijk++){//totalnum*(totalnum+1)/2
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
                    Ld << R.block(RiskFail(j,0),0,dj,1), RdR.block(RiskFail(j,0),ij,dj,1), RdR.block(RiskFail(j,0),jk,dj,1) ,RddR.block(RiskFail(j,0),ijk,dj,1);//sum of risks in group
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
                    MatrixXd temp1 = MatrixXd::Zero(Ld.rows(),1);
                    MatrixXd temp2 = MatrixXd::Zero(Ld.rows(),1);
                    temp1 = Ld.col(0).array().log();
                    double Ld1 =  (temp1.array().isFinite()).select(temp1,0).sum();
                    temp1 = Ld.col(1).array();
                    temp2 = Ld.col(2).array();
                    double Ld2 = (temp1.array().isFinite()).select(temp1,0).sum();
                    temp1 = Ld.col(3).array() - temp1.array() * temp2.array();
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
        }
        for (int ij=0;ij<totalnum;ij++){
            if (abs(Lld[ij]) > Lld_worst){
                Lld_worst = abs(Lld[ij]);
            }
        }
        if (iteration > totalnum){//Sets the minimum number of iterations
            if (iteration % (totalnum)){//Checks every set number of iterations
                if (Lld_worst < deriv_epsilon){//ends if the derivatives are low enough
                    iteration = maxiter;
                }
                Ll_comp[1]=Ll[0];
                if (abs(Ll_comp[1]-Ll_comp[0])<10){//if the change in log-likelihood isn't high enough, the maximum step size if reduced
                    abs_max = abs_max*0.1;
                    dose_abs_max = dose_abs_max*0.5;
                }
                if (abs_max < epsilon/10){//if the maximum change is too low, then it ends
                    iteration = maxiter;
                }
            }
        }
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
    }
    // Changes the parameter back into the original form
    #pragma omp parallel for schedule(dynamic) num_threads(nthreads)
    for (int ij=0;ij<totalnum;ij++){
        int ind0 = ij;
        if (ind0 < dose_num_tot){
            ;
        } else {
            ind0 = ind0-dose_num_tot;
            if (include_bool[0]==1){
                if (ind0 < beta_lin.size()){
                    beta_lin[ind0] = beta_0[ij];
                    //
                } else {
                    //one exists and its not one
                    ind0 = ind0 - beta_lin.size();
                    if (include_bool[1]==1){
                        if (ind0 < beta_loglin.size()){
                            //one and two exists and is two
                            beta_loglin[ind0] = beta_0[ij];
                            //
                        } else{
                            //one exists, two does, must be three
                            if (include_bool[2]!=1){
                                throw invalid_argument( "Are all three used? 0" );
                            }
                            ind0 = ind0 - beta_loglin.size();
                            beta_plin[ind0] = beta_0[ij];
                            //
                        }
                    } else{
                        //one exists, and two doesn't exist, must be three
                        if (include_bool[2]!=1){
                            throw invalid_argument( "Are all first and third used?" );
                        }
                        beta_plin[ind0] = beta_0[ij];
                        //
                    }
                }
            }else{
                //one doesn't exist
                if (include_bool[1]==1){
                    if (ind0 < beta_loglin.size()){
                        //one doesn't exist and two exists and is two
                        beta_loglin[ind0] = beta_0[ij];
                        //
                    } else{
                        //one doesn't exist, two does, must be three
                        if (include_bool[2]!=1){
                            throw invalid_argument( "Are all three used? 1" );
                        }
                        ind0 = ind0 - beta_loglin.size();
                        beta_plin[ind0] = beta_0[ij];
                        //
                    }
                } else{
                    //one doesn't exist, and two doesn't exist, must be three
                    if (include_bool[2]!=1){
                        throw invalid_argument( "Are all first and third used?" );
                    }
                    //
                }
            }
        }
    }
//    cout << "Making the list" << endl;
    vector<double> beta_loglin_slope;
    vector<double> beta_loglin_top;
    vector<double> beta_lin_slope;
    vector<double> beta_lin_int;
    vector<double> beta_quad;
    vector<double> beta_step_slope;
    vector<double> beta_step_int;
    int loglin_size=0;
    int lin_size=0;
    int quad_size=0;
    int step_size=0;
    int dub_off=0;
    for (int ijk=0;ijk<df_dose.cols();ijk++){
        beta_loglin_slope = beta_loglin_slopes_CPP[ijk];
        beta_loglin_top = beta_loglin_tops_CPP[ijk];
        beta_lin_slope = beta_lin_slopes_CPP[ijk];
        beta_lin_int = beta_lin_ints_CPP[ijk];
        beta_quad = beta_quads_CPP[ijk];
        beta_step_slope = beta_step_slopes_CPP[ijk];
        beta_step_int = beta_step_ints_CPP[ijk];
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
        total_dose = loglin_size + lin_size + quad_size + step_size;//beta_loglin_slope.size() + beta_lin_slope.size() + beta_quad.size() + beta_step_int.size();
        dub_off=0;
        if ((beta_loglin_slope[0]==1)&&(loglin_size==1)){
            dub_off=1;
        }
        #pragma omp parallel for schedule(dynamic) num_threads(nthreads)
        for (int ij=0; ij<total_dose;ij++){
            if (ij < loglin_size){
                if ((beta_loglin_slope[ij]==1)&&(loglin_size==1)){
                    int ind0 = cumulative_dose_num[ijk]+ij;
                    //
                    beta_loglin_top[ij,0]=beta_0[ind0];
                } else {
                    int ind0 = cumulative_dose_num[ijk]+2*ij;
                    int ind1 = ind0 + 1; 
                    //
                    beta_loglin_slope[ij,0]=beta_0[ind0];
                    beta_loglin_top[ij,0]=beta_0[ind1];
                }
                
            } else if (ij < loglin_size + lin_size){
                int jk = ij - loglin_size;
                int ind0 = cumulative_dose_num[ijk]+2*loglin_size - dub_off  + 2*jk;
                int ind1 = ind0 + 1; 
                //
                beta_lin_slope[ij,0]=beta_0[ind0];
                beta_lin_int[ij,0]=beta_0[ind1];
            } else if (ij < loglin_size + lin_size + quad_size){
                int jk = ij - loglin_size - lin_size;
                int ind0 = cumulative_dose_num[ijk]+2*loglin_size - dub_off  + 2*lin_size+jk;
                //
                beta_quad[ij,0]=beta_0[ind0];
            } else {
                int jk = ij - loglin_size - lin_size - quad_size;
                int ind0 = cumulative_dose_num[ijk]+2*loglin_size - dub_off  + 2*lin_size + quad_size + 2*jk;
                int ind1 = ind0 + 1;
                //
                beta_step_slope[ij,0]=beta_0[ind0];
                beta_step_int[ij,0]=beta_0[ind1];
            }
        }
        //
        beta_loglin_slopes[ijk] = wrap(beta_loglin_slope);
        beta_loglin_tops[ijk] = wrap(beta_loglin_top);
        beta_lin_slopes[ijk] = wrap(beta_lin_slope);
        beta_lin_ints[ijk] = wrap(beta_lin_int);
        beta_quads[ijk] = wrap(beta_quad);
        beta_step_slopes[ijk] = wrap(beta_step_slope);
        beta_step_ints[ijk] = wrap(beta_step_int);
        //
    }
    List para_list = List::create(_["beta_loglin_slopes"]=beta_loglin_slopes,_["beta_loglin_tops"]=beta_loglin_tops,_["beta_lin_slopes"]=beta_lin_slopes,_["beta_lin_ints"]=beta_lin_ints, _["beta_quads"]=beta_quads,_["beta_step_slopes"]=beta_step_slopes,_["beta_step_ints"]=beta_step_ints, _["beta_lin"]=beta_lin,_["beta_loglin"]=beta_loglin,_["beta_plin"]=beta_plin);
    NumericVector Lldd_vec = wrap(Lldd);//creates list of dose parameters
    Lldd_vec.attr("dim") = Dimension(totalnum, totalnum);
    //
    const Map<MatrixXd> Lldd_mat(as<Map<MatrixXd> >(Lldd_vec));
    //
    MatrixXd Lldd_inv = -1 * Lldd_mat.inverse().matrix(); //uses inverse information matrix to calculate the standard deviation
    //
    List res_list = List::create(_["LogLik"]=wrap(Ll),_["First_Der"]=wrap(Lld),_["Second_Der"]=Lldd_vec,_["beta_0"]=wrap(beta_0) ,_["Standard_Deviation"]=wrap(Lldd_inv.diagonal().cwiseSqrt()) ,_["AIC"]=2*totalnum-2*Ll[fir],_["Parameter_Lists"]=para_list);
    // returns a list of results
    return res_list;
}


List PEANUT_PLOT_SURV(VectorXd beta_linT,VectorXd beta_loglinT,VectorXd beta_plinT,MatrixXd df_lin,MatrixXd df_loglin,MatrixXd df_plin, MatrixXd df_dose,int fir,string modelform,int ntime, NumericVector include_bool,List beta_loglin_slopes, List beta_loglin_tops , List beta_lin_slopes , List beta_lin_ints , List beta_quads , List beta_step_slopes , List beta_step_ints, NumericMatrix df_groups, NumericVector tu ,double dose_abs_max){
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
    vector<int> dose_breaks(df_dose.cols(),0);
    int dose_num_tot=0;
    int dose_term_tot=0;
    int loglin_size=0;
    int lin_size=0;
    int quad_size=0;
    int step_size=0;
    // Converts from a list of Rcpp vectors to a std::vector of std::vectors
// There were possible issues with Rcpp vectors being created in OMP sections
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
        // Establishes the number of dose terms needed
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
        // Gathers the dose terms per dose column
        dose_term_tot += dose_breaks[ijk];
        //
    }
    //
    int totalnum = dose_num_tot;
    //
    cout.precision(10); //forces higher precision numbers printed to terminal
    int nthreads = Eigen::nbThreads()-1; //stores how many threads are allocated
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
    //
    double Lld_worst = 0.0; //stores derivative value used to determine if every parameter is near convergence
    vector <string> tform(totalnum);// list of term types
    double totem = df_loglin.rows();//precalculates how many rows are needed
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
    //
    VectorXd beta_0(totalnum);
    MatrixXd df0 = MatrixXd::Zero(df_lin.rows(), totalnum); // stores memory for the parameter columns
    MatrixXd T0 = MatrixXd::Zero(df0.rows(), totalnum); //preallocates matrix for Term column
    MatrixXd Td0 = MatrixXd::Zero(df0.rows(), totalnum); //preallocates matrix for Term derivative columns
    MatrixXd Tdd0 = MatrixXd::Zero(df0.rows(), totalnum); //preallocates matrix for Term second derivative columns
    //
    MatrixXd Te = MatrixXd::Zero(df0.rows(), 1); //preallocates matrix for non-Derivative column terms
    MatrixXd R = MatrixXd::Zero(df0.rows(), 1); //preallocates matrix for Risks
    MatrixXd Rd = MatrixXd::Zero(df0.rows(), totalnum); //preallocates matrix for Risk derivatives
    MatrixXd Rdd = MatrixXd::Zero(df0.rows(), totalnum*(totalnum+1)/2); //preallocates matrix for Risk second derivatives
    //
    MatrixXd De = MatrixXd::Zero(df_dose.rows(),dose_num_tot); //matrix of dose term values
    MatrixXd Dde = MatrixXd::Zero(df_dose.rows(),dose_num_tot); //matrix of dose term derivatives
    MatrixXd Ddde = MatrixXd::Zero(df_dose.rows(),dose_num_tot*(dose_num_tot+1)/2); //matrix of dose term second derivatives
    double dint = dose_abs_max; //The amount of change used to calculate derivatives in threshold paramters
    int total_dose=0; //used later on for a section summing the dose terms
    //
    ArrayXd Dose(df_dose.rows(),1); //Matrix of the total dose term values
    #pragma omp parallel for schedule(dynamic) num_threads(nthreads)
    for (int ij=0;ij<(totalnum-dose_num_tot+dose_term_tot);ij++){
        if (ij < dose_term_tot){//Dose terms are indexed originally by complete term, i.e. step function not the slope and threshold seperately
            int ind0 = ij;
            int ijk=0;
            while (ind0>dose_breaks[ijk]){
                ind0=ind0 - dose_breaks[ijk];
                ijk++;
            } //determines which dose column is used
            //
            int loglin_size=0;
            int lin_size=0;
            int quad_size=0;
            int step_size=0;
            vector<double> beta_loglin_slope = beta_loglin_slopes_CPP[ijk];
            vector<double> beta_loglin_top = beta_loglin_tops_CPP[ijk];
            vector<double> beta_lin_slope = beta_lin_slopes_CPP[ijk];
            vector<double> beta_lin_int = beta_lin_ints_CPP[ijk];
            vector<double> beta_quad = beta_quads_CPP[ijk];
            vector<double> beta_step_slope = beta_step_slopes_CPP[ijk];
            vector<double> beta_step_int = beta_step_ints_CPP[ijk];
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
            int dub_off=0;
            if ((beta_loglin_slope[0]==1)&&(loglin_size==1)){
                dub_off=1; //if the exponential term has a slope of 1, it isn't changed. This corrects the number of parameters
            }
            if (ind0 < loglin_size){
                ArrayXd temp = (beta_loglin_top[ind0] * df_dose.col(ijk)).array().exp();
                ArrayXd temp1 = beta_loglin_slope[ind0] * temp;
                // in both cases the total exponential term is calculated
                //
                if ((beta_loglin_slope[ind0]==1)&&(loglin_size==1)){
                    int ind1 = cumulative_dose_num[ijk]+ind0; //skips slope term
                    //
                    beta_0[ind1] = beta_loglin_top[ind0];
                    tform[ind1] = "loglin_top";
                    df0.col(ind1) = df_dose.col(ijk);
                    //
                    De.col(ind1) = temp1;
                    Dose = Dose + temp1.array();
                    Dde.col(ind1) = temp1.array() * df_dose.col(ijk).array();
                    Ddde.col(ind1 * (ind1+1)/2) = temp1.array() * df_dose.col(ijk).array().square().array();
                    //
                } else {
                    int ind1 = cumulative_dose_num[ijk]+2*ind0; //does not skip dose terms
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
                    Dose = Dose + temp1.array();
                    Dde.col(ind1) = temp.array();
                    Dde.col(ind1 + 1) = temp1.array() * df_dose.col(ijk).array();
                    Ddde.col((ind1 + 1) * (ind1 + 2)/2 + ind1) = temp.array() * df_dose.col(ijk).array();
                    Ddde.col((ind1 + 1) * (ind1 + 2)/2 + ind1 + 1) = temp1.array() * df_dose.col(ijk).array().square().array();
                }
                
            } else if (ind0 < loglin_size + lin_size){
                int jk = ind0 - loglin_size; // The linear threshold model has no closed form equation for derivatives
                ArrayXd temp = (df_dose.col(ijk).array() - beta_lin_int[jk]);
                ArrayXd temp0 = (df_dose.col(ijk).array() - beta_lin_int[jk]-dint);//derivatives are approximated
                ArrayXd temp1 = (df_dose.col(ijk).array() - beta_lin_int[jk]+dint);
                //
                int ind1 = cumulative_dose_num[ijk]+2*loglin_size - dub_off + 2*jk;
    //                int ind2 = (ind1 + 1); 
                //
                beta_0[ind1] = beta_lin_slope[jk];
                beta_0[(ind1 + 1)] = beta_lin_int[jk];
                tform[ind1] = "lin_slope";
                tform[(ind1 + 1)] = "lin_ind";
                df0.col(ind1) = df_dose.col(ijk);
                df0.col((ind1 + 1)) = df_dose.col(ijk);
                //
                temp = (temp.array() < 0).select(0.0, temp);
                temp0 = (temp0.array() < 0).select(0.0, temp0);
                temp1 = (temp1.array() < 0).select(0.0, temp1); //applies thresholds
                //
                De.col(ind1) = beta_lin_slope[jk] * temp.array();
                De.col((ind1 + 1)) = beta_lin_slope[jk] * temp.array();
                Dose = Dose + De.col(ind1).array();
                Dde.col(ind1) = temp.array();
                Dde.col((ind1 + 1)) = beta_lin_slope[jk] * (temp1.array() -temp0.array()) / 2.0/dint;
                //
                Ddde.col((ind1 + 1) * ((ind1 + 1)+1)/2 + ind1) = (temp1.array() -temp0.array()) / 2.0/dint;
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
                Dose = Dose + De.col(ind1).array();
            } else {
                int jk = ind0 - loglin_size - lin_size - quad_size;
                ArrayXd temp = (df_dose.col(ijk).array() - beta_step_int[jk]);// Step function also has approximated derivatives
                ArrayXd temp0 = (df_dose.col(ijk).array() - beta_step_int[jk]-dint);
                ArrayXd temp1 = (df_dose.col(ijk).array() - beta_step_int[jk]+dint);
                //
                int ind1 = cumulative_dose_num[ijk]+2*loglin_size - dub_off  + 2*lin_size + quad_size + 2*jk;
    //                int ind2 = ind1 + 1;
                //
                beta_0[ind1] = beta_step_slope[jk];
                beta_0[(ind1 + 1)] = beta_step_int[jk];
                tform[ind1] = "step_slope";
                tform[(ind1 + 1)] = "step_ind";
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
                Dose = Dose + De.col(ind1).array();
            }
            //
        } else {//past this point, non-dose terms are calcualted
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
            T0.col(ind1) = (df0.col(ind1).array() * beta_0[ind1]).matrix();//a similar form is used for every form
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
    T0.block(0,0,T0.rows(),dose_num_tot) = De.block(0,0,T0.rows(),dose_num_tot);
    Td0.block(0,0,T0.rows(),dose_num_tot) = Dde.block(0,0,T0.rows(),dose_num_tot);
    //
    #pragma omp parallel for schedule(dynamic) num_threads(nthreads)
    for (int ijk=0;ijk<dose_num_tot;ijk++){//inserts the dose terms into the total term arrays
        Tdd0.col(ijk) = Ddde.col(ijk*(ijk+1)/2+ijk);
        Td0.col(ijk) = Dde.col(ijk);
        T0.col(ijk) = De.col(ijk);
    }
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
                            Rdd.block(0,ijk,df0.rows(),1) = Tdd0.block(0,ij,df0.rows(),1).array() * Te.block(0,0,df0.rows(),1).array();
                        } else {
                            Rdd.block(0,ijk,df0.rows(),1) = Tdd0.block(0,ij,df0.rows(),1).array() * T0.block(0,fir,df0.rows(),1).array();
                        }
                    } else {
                        if (ij<dose_num_tot){
                            Rdd.block(0,ijk,df0.rows(),1) = Ddde.block(0,ijk,df0.rows(),1).array() * Te.block(0,0,df0.rows(),1).array();
                        } else {
                            Rdd.block(0,ijk,df0.rows(),1) = Tdd0.block(0,ij,df0.rows(),1).array() * Dose.block(0,0,df0.rows(),1).array();
                        }
                    }
                } else {
                    if (fir>=dose_num_tot){
                        if ((ij==fir)||(jk==fir)){
                            Rdd.block(0,ijk,df0.rows(),1) = Td0.block(0,ij,df0.rows(),1).array() * Td0.block(0,jk,df0.rows(),1).array();
                        } else if ((ij<dose_num_tot)&&(jk<dose_num_tot)){
                            Rdd.block(0,ijk,df0.rows(),1) = Ddde.block(0,ijk,df0.rows(),1).array() * T0.block(0,fir,df0.rows(),1).array();
                        }
                    } else {
                        if ((ij<dose_num_tot)&&(jk<dose_num_tot)){
                            Rdd.block(0,ijk,df0.rows(),1) = Te.block(0,0,df0.rows(),1).array() * Ddde.block(0,ijk,df0.rows(),1).array();
                        } else if ((ij>=dose_num_tot)&&(jk<dose_num_tot)){
                            Rdd.block(0,ijk,df0.rows(),1) = Td0.block(0,ij,df0.rows(),1).array() * Td0.block(0,jk,df0.rows(),1).array();
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
        Te = T0.array().block(0,dose_num_tot,T0.rows(),T0.cols()-dose_num_tot).rowwise().prod() * Dose.array();
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
                    Rdd.block(0,ijk,df0.rows(),1) = Ddde.block(0,ijk,df0.rows(),1).array() * Te.block(0,0,df0.rows(),1).array();
                } else {
                    Rdd.block(0,ijk,df0.rows(),1) = Tdd0.block(0,jk,df0.rows(),1).array() * T0.block(0,jk,df0.rows(),1).array().pow(-1).array() * R.block(0,0,df0.rows(),1).array();
                }
            } else {
                if ((ij<dose_num_tot)&&(jk<dose_num_tot)){
                    Rdd.block(0,ijk,df0.rows(),1) = Ddde.block(0,ijk,df0.rows(),1).array() * Te.block(0,0,df0.rows(),1).array();
                } else if ((ij<dose_num_tot)||(jk<dose_num_tot)){
                    Rdd.block(0,ijk,df0.rows(),1) = Dde.block(0,jk,df0.rows(),1).array() * De.block(0,jk,df0.rows(),1).array().pow(-1).array() * Rd.block(0,ij,df0.rows(),1).array();
                } else{
                    Rdd.block(0,ijk,df0.rows(),1) = Td0.block(0,jk,df0.rows(),1).array() * T0.block(0,jk,df0.rows(),1).array().pow(-1).array() * Rd.block(0,ij,df0.rows(),1).array();
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
    cout<<"df100 "<<(ending-start)<<" "<<0<<" "<<0<<" "<<-1<<",Prep_R"<<endl;
    gibtime = system_clock::to_time_t(system_clock::now());
    cout << ctime(&gibtime) << endl;
    //
    // -------------------------------------------------------------------------------------------
    //
    vector<double> baseline(ntime,0.0);
    //
    const Map<MatrixXd> df_m(as<Map<MatrixXd> >(df_groups));
    #pragma omp parallel for schedule(dynamic) num_threads(nthreads)
    for (int ijk=0;ijk<ntime;ijk++){
        double t0 = tu[ijk];
        VectorXi select_ind_all = ((df_m.col(0).array() < t0)&&(df_m.col(1).array()>=t0)).cast<int>(); //indices at risk
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
    List res_list = List::create(_["baseline"]=wrap(baseline), _["Risks"]=wrap(R));
    //
    return res_list;
}

List PEANUT_PLOT_RISK(VectorXd beta_linT,VectorXd beta_loglinT,VectorXd beta_plinT,MatrixXd df_lin,MatrixXd df_loglin,MatrixXd df_plin, MatrixXd df_dose,int fir,string modelform,int ntime, NumericVector include_bool,List beta_loglin_slopes, List beta_loglin_tops , List beta_lin_slopes , List beta_lin_ints , List beta_quads , List beta_step_slopes , List beta_step_ints, NumericMatrix df_groups, NumericVector tu,int uniq_v ,double dose_abs_max){
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
    vector<int> dose_breaks(df_dose.cols(),0);
    int dose_num_tot=0;
    int dose_term_tot=0;
    int loglin_size=0;
    int lin_size=0;
    int quad_size=0;
    int step_size=0;
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
        // Establishes the number of dose terms needed
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
        // Gathers the dose terms per dose column
        dose_term_tot += dose_breaks[ijk];
        //
    }
    //
    int totalnum = dose_num_tot;
    //
    cout.precision(10); //forces higher precision numbers printed to terminal
    int nthreads = Eigen::nbThreads()-1; //stores how many threads are allocated
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
    //
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
    double dint = dose_abs_max;
    int dub_off=0;
    int total_dose=0;
    double beta_0;
    string tform;
    //
    float dx = 0;
    if (fir >=0){
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
    MatrixXd df0 = MatrixXd::Zero(vv.size(), 1); // stores memory for the derivative term parameters and columns
    MatrixXd T0 = MatrixXd::Zero(vv.size(), totalnum); //preallocates matrix for Derivative column terms
    if (fir < df_dose.cols()){;
        int ijk=fir;
        dx = (df_dose.col(ijk).maxCoeff() - df_dose.col(ijk).minCoeff())/vv.size();//varies from max to minimum
        vv[0] = df_dose.col(ijk).minCoeff();
        generate(vv.begin(), vv.end(), [n = 0, &dx]() mutable { return n++ * dx; });
        #pragma omp parallel for schedule(dynamic) num_threads(nthreads)
        for (int ij=0;ij<vv.size();ij++){
            df0(ij,0)=vv[ij];
        }
        //
        vector<double> beta_loglin_slope = beta_loglin_slopes_CPP[ijk];
        vector<double> beta_loglin_top = beta_loglin_tops_CPP[ijk];
        vector<double> beta_lin_slope = beta_lin_slopes_CPP[ijk];
        vector<double> beta_lin_int = beta_lin_ints_CPP[ijk];
        vector<double> beta_quad = beta_quads_CPP[ijk];
        vector<double> beta_step_slope = beta_step_slopes_CPP[ijk];
        vector<double> beta_step_int = beta_step_ints_CPP[ijk];
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
        dub_off=0;
        if ((beta_loglin_slope[0]==1)&&(loglin_size==1)){
            dub_off=1; //if the exponential term has a slope of 1, it isn't changed. This corrects the number of parameters
        }
        //
        #pragma omp parallel for schedule(dynamic) num_threads(nthreads)
        for (int ind0=0;ind0<dose_breaks[ijk];ind0++){
            if (ind0 < loglin_size){
                ArrayXd temp = (beta_loglin_top[ind0] * df0.col(0)).array().exp();
                ArrayXd temp1 = beta_loglin_slope[ind0] * temp;
                //
                T0.col(ind0) = temp1.array();
                
            } else if (ind0 < loglin_size + lin_size){
                int jk = ind0 - loglin_size;
                ArrayXd temp = (df0.col(0).array() - beta_lin_int[jk]);
                //
                T0.col(ind0) = beta_lin_slope[jk] * (temp.array() < 0).select(0.0, temp);
            } else if (ind0 < loglin_size + lin_size + quad_size){
                int jk = ind0 - loglin_size - lin_size;
                ArrayXd temp = df0.col(0).array().square();
                //
                T0.col(ind0) = beta_quad[jk] * temp.array();
            } else {
                int jk = ind0 - loglin_size - lin_size - quad_size;
                ArrayXd temp = (df0.col(0).array() - beta_step_int[jk]);
                //
                T0.col(ind0) = beta_step_slope[jk] * (temp.array() < 0).select(0.0, MatrixXd::Zero(temp.rows(),temp.cols()).array()+1.0);
            }
        }
        //
    } else {
        int ind0 = fir-df_dose.cols();
        double vmax;
        double vmin;
        if (include_bool[0]==1){
            if (ind0 < beta_lin.size()){
                // one exists and is one
                beta_0 = beta_lin[ind0];
                vmax=df_lin.col(ind0).maxCoeff();
                vmin=df_lin.col(ind0).minCoeff();
                tform="lin";
                //
            } else {
                //one exists and its not one
                ind0 = ind0 - beta_lin.size();
                if (include_bool[1]==1){
                    if (ind0 < beta_loglin.size()){
                        //one and two exists and is two
                        beta_0 = beta_loglin[ind0];
                        vmax=df_loglin.col(ind0).maxCoeff();
                        vmin=df_loglin.col(ind0).minCoeff();
                        tform="loglin";
                        //
                    } else{
                        //one exists, two does, must be three
                        if (include_bool[2]!=1){
                            throw invalid_argument( "Are all three used? 0" );
                        }
                        ind0 = ind0 - beta_loglin.size();
                        beta_0 = beta_plin[ind0];
                        vmax=df_plin.col(ind0).maxCoeff();
                        vmin=df_plin.col(ind0).minCoeff();
                        tform="plin";
                        //
                    }
                } else{
                    //one exists, and two doesn't exist, must be three
                    if (include_bool[2]!=1){
                        throw invalid_argument( "Are all first and third used?" );
                    }
                    beta_0 = beta_plin[ind0];
                    vmax=df_plin.col(ind0).maxCoeff();
                    vmin=df_plin.col(ind0).minCoeff();
                    tform="plin";
                    //
                }
            }
        }else{
            //one doesn't exist
            if (include_bool[1]==1){
                if (ind0 < beta_loglin.size()){
                    //one doesn't exist and two exists and is two
                    beta_0 = beta_loglin[ind0];
                    vmax=df_loglin.col(ind0).maxCoeff();
                    vmin=df_loglin.col(ind0).minCoeff();
                    tform="loglin";
                    //
                } else{
                    //one doesn't exist, two does, must be three
                    if (include_bool[2]!=1){
                        throw invalid_argument( "Are all three used? 1" );
                    }
                    ind0 = ind0 - beta_loglin.size();
                    beta_0 = beta_plin[ind0];
                    vmax=df_plin.col(ind0).maxCoeff();
                    vmin=df_plin.col(ind0).minCoeff();
                    tform="plin";
                    //
                }
            } else{
                //one doesn't exist, and two doesn't exist, must be three
                if (include_bool[2]!=1){
                    throw invalid_argument( "Are all first and third used?" );
                }
                beta_0 = beta_plin[ind0];
                vmax=df_plin.col(ind0).maxCoeff();
                vmin=df_plin.col(ind0).minCoeff();
                tform="plin";

                //
            }
        }
        dx = (vmax - vmin)/vv.size();
        vv[0] = vmin;
        generate(vv.begin(), vv.end(), [n = 0, &dx]() mutable { return n++ * dx; });
        #pragma omp parallel for schedule(dynamic) num_threads(nthreads)
        for (int ij=0;ij<vv.size();ij++){
            df0(ij,0)=vv[ij];
        }
        //
        T0.col(fir) = (df0.array() * beta_0).matrix();

        if (tform=="lin") {
            ;
        } else if (tform=="loglin") {
            T0.col(fir) = T0.col(fir).array().exp();
        } else if (tform=="plin") {
            T0.col(fir) = 1 + T0.col(fir).array();
        } else {
            cout << tform << " is invalid" << endl;
            throw invalid_argument( "Invalid term type" );
        }
    }
    List res_list = List::create(_["x"]=wrap(df0), _["y"]=wrap(T0.rowwise().sum()));//returns list of covariate values and risk
    return res_list;
}




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


//
// --------------------------------------------------------------------------------------------------------------------------------------------
//




List LogLik_Bounds( double q1, VectorXd beta_linT,VectorXd beta_loglinT,VectorXd beta_plinT,MatrixXd df_lin,MatrixXd df_loglin,MatrixXd df_plin, MatrixXd df_dose,int fir,string modelform,int ntime, NumericVector include_bool, double lr, int maxiter, int halfmax, double epsilon, double dbeta_cap, double abs_max, double dose_abs_max, double deriv_epsilon,List beta_loglin_slopes, List beta_loglin_tops , List beta_lin_slopes , List beta_lin_ints , List beta_quads , List beta_step_slopes , List beta_step_ints, NumericMatrix df_groups, NumericVector tu ){
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
    ArrayXd Dose(df_dose.rows(),1);
    //
    int dub_off=0;
    #pragma omp parallel for schedule(dynamic) num_threads(nthreads)
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
                    Dose = Dose + temp1.array();
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
                    Dose = Dose + temp1.array();
                    Dde.col(ind1) = temp.array();
                    Dde.col(ind1 + 1) = temp1.array() * df_dose.col(ijk).array();
                    Ddde.col((ind1 + 1) * (ind1 + 2)/2 + ind1) = temp.array() * df_dose.col(ijk).array();
                    Ddde.col((ind1 + 1) * (ind1 + 2)/2 + ind1 + 1) = temp1.array() * df_dose.col(ijk).array().square().array();
                }
                
            } else if (ind0 < loglin_size + lin_size){
                int jk = ind0 - loglin_size;
                ArrayXd temp = (df_dose.col(ijk).array() - beta_lin_int[jk]);
                ArrayXd temp0 = (df_dose.col(ijk).array() - beta_lin_int[jk]-dint);
                ArrayXd temp1 = (df_dose.col(ijk).array() - beta_lin_int[jk]+dint);
                //
                int ind1 = cumulative_dose_num[ijk]+2*loglin_size - dub_off + 2*jk;
                //
                beta_0[ind1] = beta_lin_slope[jk];
                beta_0[(ind1 + 1)] = beta_lin_int[jk];
                tform[ind1] = "lin_slope";
                tform[(ind1 + 1)] = "lin_ind";
                df0.col(ind1) = df_dose.col(ijk);
                df0.col((ind1 + 1)) = df_dose.col(ijk);
                //
                temp = (temp.array() < 0).select(0.0, temp);
                temp0 = (temp0.array() < 0).select(0.0, temp0);
                temp1 = (temp1.array() < 0).select(0.0, temp1);
                //
                De.col(ind1) = beta_lin_slope[jk] * temp.array();
                De.col((ind1 + 1)) = beta_lin_slope[jk] * temp.array();
                Dose = Dose + De.col(ind1).array();
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
                Dose = Dose + De.col(ind1).array();
            } else {
                int jk = ind0 - loglin_size - lin_size - quad_size;
                ArrayXd temp = (df_dose.col(ijk).array() - beta_step_int[jk]);
                ArrayXd temp0 = (df_dose.col(ijk).array() - beta_step_int[jk]-dint);
                ArrayXd temp1 = (df_dose.col(ijk).array() - beta_step_int[jk]+dint);
                //
                int ind1 = cumulative_dose_num[ijk]+2*loglin_size - dub_off  + 2*lin_size + quad_size + 2*jk;
    //                int ind2 = ind1 + 1;
                //
                beta_0[ind1] = beta_step_slope[jk];
                beta_0[(ind1 + 1)] = beta_step_int[jk];
                tform[ind1] = "step_slope";
                tform[(ind1 + 1)] = "step_ind";
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
                Dose = Dose + De.col(ind1).array();
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
    T0.block(0,0,T0.rows(),dose_num_tot) = De.block(0,0,T0.rows(),dose_num_tot);
    Td0.block(0,0,T0.rows(),dose_num_tot) = Dde.block(0,0,T0.rows(),dose_num_tot);
    //
    #pragma omp parallel for schedule(dynamic) num_threads(nthreads)
    for (int ijk=0;ijk<dose_num_tot;ijk++){
        Tdd0.col(ijk) = Ddde.col(ijk*(ijk+1)/2+ijk);
        Td0.col(ijk) = Dde.col(ijk);
        T0.col(ijk) = De.col(ijk);
    }
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
                            Rdd.block(0,ijk,df0.rows(),1) = Tdd0.block(0,ij,df0.rows(),1).array() * Te.block(0,0,df0.rows(),1).array();
                        } else {
                            Rdd.block(0,ijk,df0.rows(),1) = Tdd0.block(0,ij,df0.rows(),1).array() * T0.block(0,fir,df0.rows(),1).array();
                        }
                    } else {
                        if (ij<dose_num_tot){
                            Rdd.block(0,ijk,df0.rows(),1) = Ddde.block(0,ijk,df0.rows(),1).array() * Te.block(0,0,df0.rows(),1).array();
                        } else {
                            Rdd.block(0,ijk,df0.rows(),1) = Tdd0.block(0,ij,df0.rows(),1).array() * Dose.block(0,0,df0.rows(),1).array();
                        }
                    }
                } else {
                    if (fir>=dose_num_tot){
                        if ((ij==fir)||(jk==fir)){
                            Rdd.block(0,ijk,df0.rows(),1) = Td0.block(0,ij,df0.rows(),1).array() * Td0.block(0,jk,df0.rows(),1).array();
                        } else if ((ij<dose_num_tot)&&(jk<dose_num_tot)){
                            Rdd.block(0,ijk,df0.rows(),1) = Ddde.block(0,ijk,df0.rows(),1).array() * T0.block(0,fir,df0.rows(),1).array();
                        }
                    } else {
                        if ((ij<dose_num_tot)&&(jk<dose_num_tot)){
                            Rdd.block(0,ijk,df0.rows(),1) = Te.block(0,0,df0.rows(),1).array() * Ddde.block(0,ijk,df0.rows(),1).array();
                        } else if ((ij>=dose_num_tot)&&(jk<dose_num_tot)){
                            Rdd.block(0,ijk,df0.rows(),1) = Td0.block(0,ij,df0.rows(),1).array() * Td0.block(0,jk,df0.rows(),1).array();
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
        Te = T0.array().block(0,dose_num_tot,T0.rows(),T0.cols()-dose_num_tot).rowwise().prod() * Dose.array();
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
    //        #pragma omp parallel for schedule(dynamic) num_threads(nthreads)
    //        for (int ijk=0;ijk<totalnum;ijk++){
    //            if (ijk<dose_num_tot){
    //                Rd.col(ijk) = De.col(ijk).array().pow(-1).array()* Te.array();
    //                Rd.col(ijk) = Rd.col(ijk).array() * Dde.col(ijk).array();
    //            }
    //        }
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
                    Rdd.block(0,ijk,df0.rows(),1) = Ddde.block(0,ijk,df0.rows(),1).array() * Te.block(0,0,df0.rows(),1).array();
                } else {
                    Rdd.block(0,ijk,df0.rows(),1) = Tdd0.block(0,jk,df0.rows(),1).array() * T0.block(0,jk,df0.rows(),1).array().pow(-1).array() * R.block(0,0,df0.rows(),1).array();
                }
            } else {
                if ((ij<dose_num_tot)&&(jk<dose_num_tot)){
                    Rdd.block(0,ijk,df0.rows(),1) = Ddde.block(0,ijk,df0.rows(),1).array() * Te.block(0,0,df0.rows(),1).array();
                } else if ((ij<dose_num_tot)||(jk<dose_num_tot)){
                    Rdd.block(0,ijk,df0.rows(),1) = Dde.block(0,jk,df0.rows(),1).array() * De.block(0,jk,df0.rows(),1).array().pow(-1).array() * Rd.block(0,ij,df0.rows(),1).array();
                } else{
                    Rdd.block(0,ijk,df0.rows(),1) = Td0.block(0,jk,df0.rows(),1).array() * T0.block(0,jk,df0.rows(),1).array().pow(-1).array() * Rd.block(0,ij,df0.rows(),1).array();
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
    cout<<"df100 "<<(ending-start)<<" "<<0<<" "<<0<<" "<<-1<<",Prep_R"<<endl;
    gibtime = system_clock::to_time_t(system_clock::now());
    cout << ctime(&gibtime) << endl;
    //
    // -------------------------------------------------------------------------------------------
    //
    vector<string>  RiskGroup(ntime); //vector of strings detailing the rows
    IntegerMatrix RiskFail(ntime,2); //vector giving the event rows
    //
    const Map<MatrixXd> df_m(as<Map<MatrixXd> >(df_groups));
    #pragma omp parallel for schedule(dynamic) num_threads(nthreads)
    for (int ijk=0;ijk<ntime;ijk++){
        double t0 = tu[ijk];
        VectorXi select_ind_all = ((df_m.col(0).array() < t0)&&(df_m.col(1).array()>=t0)).cast<int>();
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
            }
            MatrixXd Ld = MatrixXd::Zero(dj,4);
            Ld << R.block(RiskFail(j,0),0,dj,1), Rd.block(RiskFail(j,0),ij,dj,1), Rd.block(RiskFail(j,0),jk,dj,1) ,Rdd.block(RiskFail(j,0),ijk,dj,1);//sum of risks in group
            //
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
    //
    #pragma omp declare reduction(vec_double_plus : std::vector<double> : \
        std::transform(omp_out.begin(), omp_out.end(), omp_in.begin(), omp_out.begin(), std::plus<double>())) \
        initializer(omp_priv = omp_orig)
    #pragma omp parallel for schedule(dynamic) num_threads(nthreads) reduction(vec_double_plus:Ll,Lld,Lldd) collapse(2)
    for (int ijk=0;ijk<totalnum*(totalnum+1)/2;ijk++){//totalnum*(totalnum+1)/2
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
            Ld << R.block(RiskFail(j,0),0,dj,1), RdR.block(RiskFail(j,0),ij,dj,1), RdR.block(RiskFail(j,0),jk,dj,1) ,RddR.block(RiskFail(j,0),ijk,dj,1);//sum of risks in group
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
            MatrixXd temp1 = MatrixXd::Zero(Ld.rows(),1);
            MatrixXd temp2 = MatrixXd::Zero(Ld.rows(),1);
            temp1 = Ld.col(0).array().log();
            double Ld1 =  (temp1.array().isFinite()).select(temp1,0).sum();
            temp1 = Ld.col(1).array();
            temp2 = Ld.col(2).array();
            double Ld2 = (temp1.array().isFinite()).select(temp1,0).sum();
            temp1 = Ld.col(3).array() - temp1.array() * temp2.array();
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
        h = abs(sqrt( q1 / h)) / 2;
        //
        //
//        Theta(ijk,ijk) = Theta(ijk,ijk) + h;
        Dbeta(ijk,ijk) = h;
        
        if (ijk==0){
//            Theta.block(1,ijk,totalnum-1,1) = Theta.block(1,ijk,totalnum-1,1) + h * dOmdBet.block(0,ijk,totalnum-1,1);
            Dbeta.block(1,ijk,totalnum-1,1) =  h * dOmdBet.col(0);
        } else if (ijk==totalnum-1){
//            Theta.block(0,ijk,totalnum-1,1) = Theta.block(0,ijk,totalnum-1,1) + h * dOmdBet.block(0,ijk,totalnum-1,1);
            Dbeta.block(0,ijk,totalnum-1,1) =  h * dOmdBet.col(0);
            //
        } else {
            ;
//            Theta.block(0,ijk,ijk,1) = Theta.block(0,ijk,ijk,1) + h * dOmdBet.block(0,0,ijk,1);
//            Theta.block(ijk+1,ijk,totalnum-ijk-1,1) = Theta.block(1,ijk+1,totalnum-ijk-1,1) + h * dOmdBet.block(ijk,0,totalnum-ijk-2,1);
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
                if ((tform[ijk]=="step_ind")||(tform[ijk]=="lin_ind")){
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
            #pragma omp parallel for schedule(dynamic) num_threads(nthreads)
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
                    MatrixXd bch = ((beta_c[ijk] - beta_p[ijk]) * df0.col(ijk)).array().exp().array();
                    //
                    De.col(ijk) = ach * bch.array() * De.col(ijk).array();
                    De.col(ijk+1) = ach * bch.array() * De.col(ijk+1).array();
                    Dde.col(ijk) = bch.array() * Dde.col(ijk).array();
                    Dde.col(ijk+1) = ach * bch.array() * Dde.col(ijk+1).array();
                    Ddde.col((ijk+1)*(ijk+2)/2+ijk) = bch.array() * Ddde.col((ijk+1)*(ijk+2)/2+ijk).array();
                    Ddde.col((ijk+1)*(ijk+2)/2+ijk+1) = ach * bch.array() * Ddde.col((ijk+1)*(ijk+2)/2+ijk+1).array();
                    //
                } else if (tform[ijk]=="loglin_top"){
                    if (ijk==0){
                        beta_c[ijk] = beta_a[ijk] + Dbeta.coeff(ijk,ijk_ind);
                        //
                        De.col(ijk) = De.col(ijk).array() * ((beta_c[ijk] - beta_p[ijk]) * df0.col(ijk)).array().exp().array();
                        Dde.col(ijk) = Dde.col(ijk).array() * ((beta_c[ijk] - beta_p[ijk]) * df0.col(ijk)).array().exp().array();
                        Ddde.col((ijk)*(ijk+1)/2+ijk) = Ddde.col((ijk)*(ijk+1)/2+ijk).array() * ((beta_c[ijk] - beta_p[ijk]) * df0.col(ijk)).array().exp().array();
                        //
                    } else if (tform[ijk-1]!="loglin_slope"){
                        beta_c[ijk] = beta_a[ijk] + Dbeta.coeff(ijk,ijk_ind);
                        //
                        De.col(ijk) = De.col(ijk).array() * ((beta_c[ijk] - beta_p[ijk]) * df0.col(ijk)).array().exp().array();
                        Dde.col(ijk) = Dde.col(ijk).array() * ((beta_c[ijk] - beta_p[ijk]) * df0.col(ijk)).array().exp().array();
                        Ddde.col((ijk)*(ijk+1)/2+ijk) = Ddde.col((ijk)*(ijk+1)/2+ijk).array() * ((beta_c[ijk] - beta_p[ijk]) * df0.col(ijk)).array().exp().array();
                        //
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
                    ArrayXd temp2 = (temp1.array() < 0).select(0, temp1) - (temp0.array() < 0).select(0, temp0);
                    //
                    De.col(ijk) = beta_c[ijk] * (temp.array() < 0).select(0.0, temp);
                    De.col((ijk+1)) = beta_c[ijk] * (temp.array() < 0).select(0.0, temp);
                    Dde.col(ijk) = (temp.array() < 0).select(0.0, temp);
                    Dde.col((ijk+1)) = beta_c[ijk] * (temp2) / 2/dint;
                    //
                    Ddde.col((ijk+1)*(ijk+2)/2+ijk) = (temp2) / 2/dint;
                    Ddde.col((ijk+1)*(ijk+2)/2+ijk+1) = beta_c[ijk] * (temp2) / pow(dint,2);
                    //
                } else if (tform[ijk]=="quad_slope"){
                    beta_c[ijk] = beta_a[ijk] + Dbeta.coeff(ijk,ijk_ind);
                    De.col(ijk) = beta_c[ijk] * df0.col(ijk).array().square().array();
                    //
                } else if (tform[ijk]=="step_slope"){
                    beta_c[ijk] = beta_a[ijk] + Dbeta.coeff(ijk,ijk_ind);
                    beta_c[ijk+1] = beta_a[ijk+1] + Dbeta.coeff(ijk+1,ijk_ind);
                    
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
                    Ddde.col((ijk+1)*(ijk+2)/2+ijk) = (temp2) / 2/dint;
                    Ddde.col((ijk+1)*(ijk+2)/2+ijk+1) = beta_0[ijk] * (temp2) / pow(dint,2);
                    //
                } else {
                    ;
                }
                Theta(ijk, ijk_ind) = beta_c[ijk];
                //
            }
            //
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
                                    Rdd.block(0,ijk,df0.rows(),1) = Tdd0.block(0,ij,df0.rows(),1).array() * Te.block(0,0,df0.rows(),1).array();
                                } else {
                                    if (ij<dose_num_tot){
                                        Rdd.block(0,ijk,df0.rows(),1) = Ddde.block(0,ijk,df0.rows(),1).array() * T0.block(0,fir,df0.rows(),1).array();
                                    } else {
                                        Rdd.block(0,ijk,df0.rows(),1) = Tdd0.block(0,ij,df0.rows(),1).array() * T0.block(0,fir,df0.rows(),1).array();
                                    }
                                }
                            } else {
                                if (ij<dose_num_tot){
                                    Rdd.block(0,ijk,df0.rows(),1) = Ddde.block(0,ijk,df0.rows(),1).array() * Te.block(0,0,df0.rows(),1).array();
                                } else {
                                    Rdd.block(0,ijk,df0.rows(),1) = Tdd0.block(0,ij,df0.rows(),1).array() * De.block(0,0,df0.rows(),1).array();
                                }
                            }
                        } else {
                            if (fir!=0){
                                if ((ij==fir)||(jk==fir)){
                                    if (ij<dose_num_tot){
                                        Rdd.block(0,ijk,df0.rows(),1) = Dde.block(0,ij,df0.rows(),1).array() * Td0.block(0,jk,df0.rows(),1).array();
                                    } else {
                                        Rdd.block(0,ijk,df0.rows(),1) = Td0.block(0,ij,df0.rows(),1).array() * Td0.block(0,jk,df0.rows(),1).array();
                                    }
                                } else if ((ij<dose_num_tot)&&(jk<dose_num_tot)){
                                    Rdd.block(0,ijk,df0.rows(),1) = Ddde.block(0,ijk,df0.rows(),1).array() * Dde.block(0,jk,df0.rows(),1).array();
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
                Te = T0.array().block(0,dose_num_tot,T0.rows(),T0.cols()-dose_num_tot).rowwise().prod() * Dose.array();
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
                            Rdd.block(0,ijk,df0.rows(),1) = Ddde.block(0,ijk,df0.rows(),1).array() * De.block(0,jk,df0.rows(),1).array().pow(-1).array() * R.block(0,0,df0.rows(),1).array();
                        } else {
                            Rdd.block(0,ijk,df0.rows(),1) = Tdd0.block(0,jk,df0.rows(),1).array() * Td0.block(0,jk,df0.rows(),1).array().pow(-1).array() * Rd.block(0,ij,df0.rows(),1).array();
                        }
                    } else {
                        if ((ij<dose_num_tot)&&(jk<dose_num_tot)){
                            Rdd.block(0,ijk,df0.rows(),1) = Ddde.block(0,ijk,df0.rows(),1).array() * De.block(0,ij,df0.rows(),1).array().pow(-1).array() * R.block(0,0,df0.rows(),1).array();
                        } else if ((ij<dose_num_tot)||(jk<dose_num_tot)){
                            Rdd.block(0,ijk,df0.rows(),1) = Dde.block(0,jk,df0.rows(),1).array() * De.block(0,jk,df0.rows(),1).array().pow(-1).array() * Rd.block(0,ij,df0.rows(),1).array();
                        } else{
                            Rdd.block(0,ijk,df0.rows(),1) = Td0.block(0,jk,df0.rows(),1).array() * T0.block(0,jk,df0.rows(),1).array().pow(-1).array() * Rd.block(0,ij,df0.rows(),1).array();
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
                    }
                    MatrixXd Ld = MatrixXd::Zero(dj,4);
                    Ld << R.block(RiskFail(j,0),0,dj,1), Rd.block(RiskFail(j,0),ij,dj,1), Rd.block(RiskFail(j,0),jk,dj,1) ,Rdd.block(RiskFail(j,0),ijk,dj,1);//sum of risks in group
                    //
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
            //
            #pragma omp declare reduction(vec_double_plus : std::vector<double> : \
                std::transform(omp_out.begin(), omp_out.end(), omp_in.begin(), omp_out.begin(), std::plus<double>())) \
                initializer(omp_priv = omp_orig)
            #pragma omp parallel for schedule(dynamic) num_threads(nthreads) reduction(vec_double_plus:Ll,Lld,Lldd) collapse(2)
            for (int ijk=0;ijk<totalnum*(totalnum+1)/2;ijk++){//totalnum*(totalnum+1)/2
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
                    Ld << R.block(RiskFail(j,0),0,dj,1), RdR.block(RiskFail(j,0),ij,dj,1), RdR.block(RiskFail(j,0),jk,dj,1) ,RddR.block(RiskFail(j,0),ijk,dj,1);//sum of risks in group
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
                    MatrixXd temp1 = MatrixXd::Zero(Ld.rows(),1);
                    MatrixXd temp2 = MatrixXd::Zero(Ld.rows(),1);
                    temp1 = Ld.col(0).array().log();
                    double Ld1 =  (temp1.array().isFinite()).select(temp1,0).sum();
                    temp1 = Ld.col(1).array();
                    temp2 = Ld.col(2).array();
                    double Ld2 = (temp1.array().isFinite()).select(temp1,0).sum();
                    temp1 = Ld.col(3).array() - temp1.array() * temp2.array();
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
        #pragma omp parallel for schedule(dynamic) num_threads(nthreads)
        for (int ij=0;ij<totalnum;ij++){
            int ind0 = ij;
            if (ind0 < dose_num_tot){
                if (tform[ij]=="loglin_slope"){
                    ArrayXd temp = (beta_0[ij+1] * df0.col(ij)).array().exp();
                    ArrayXd temp1 = beta_0[ij] * temp;
                    //
                    //
                    De.col(ij) = temp1;
                    De.col(ij+1) = temp1;
                    Dde.col(ij) = temp.array();
                    Dde.col(ij+1) = temp1.array() * df0.col(ij).array();
                    Ddde.col((ij+1)*(ij+2)/2+ij) = temp.array() * df0.col(ij).array();
                    Ddde.col((ij+1)*(ij+2)/2+ij+1) = temp1.array() * df0.col(ij).array().square().array();
                    
                } else if (tform[ij]=="loglin_top"){
                    if (ij==0){
                        ArrayXd temp = (beta_0[ij] * df0.col(ij)).array().exp();
                        De.col(ij) = temp;
                        Dde.col(ij) = temp.array() * df0.col(ij).array();
                        Ddde.col(ij * (ij+1)/2.0 + ij) = temp.array() * df0.col(ij).array().square().array();
                    } else if (tform[ij-1]!="loglin_slope"){
                        ArrayXd temp = (beta_0[ij] * df0.col(ij)).array().exp();
                        De.col(ij) = temp;
                        Dde.col(ij) = temp.array() * df0.col(ij).array();
                        Ddde.col(ij * (ij+1)/2.0 + ij) = temp.array() * df0.col(ij).array().square().array();
                        //
                    } else {
                        ;
                    }
                } else if (tform[ij]=="lin_slope"){
                    ArrayXd temp = (df0.col(ij).array() - beta_0[ij]);
                    ArrayXd temp0 = (df0.col(ij).array() - beta_0[ij]-dint);
                    ArrayXd temp1 = (df0.col(ij).array() - beta_0[ij]+dint);
                    //
                    temp = (temp.array() < 0).select(0, temp);
                    temp0 = (temp0.array() < 0).select(0, temp0);
                    temp1 = (temp1.array() < 0).select(0, temp1);
                    //
                    De.col(ij) = beta_0[ij] * temp.array();
                    De.col(ij+1) = beta_0[ij] * temp.array();
                    Dde.col(ij) = temp.array();
                    Dde.col(ij+1) = beta_0[ij] * (temp1.array()-temp0.array()) / 2/dint;
                    //
                    Ddde.col((ij+1)*(ij+2)/2+ij) = (temp1.array()-temp0.array()) / 2/dint;
                    Ddde.col((ij+1)*(ij+2)/2+ij+1) = beta_0[ij] * (temp1.array()-2*temp.array()+temp0.array()) / pow(dint,2);
                } else if (tform[ij]=="quad_slope"){
                    ArrayXd temp = df0.col(ij).array().square();
                    //
                    De.col(ij) = beta_0[ij] * temp.array();
                    Dde.col(ij) = temp.array();
                } else if (tform[ij]=="step_slope"){
                    ArrayXd temp = (df0.col(ij).array() - beta_0[ij]);
                    ArrayXd temp0 = (df0.col(ij).array() - beta_0[ij]-dint);
                    ArrayXd temp1 = (df0.col(ij).array() - beta_0[ij]+dint);
                    //
                    temp = (temp.array() < 0).select(0.0, MatrixXd::Zero(temp.rows(),temp.cols()).array()+1.0);
                    temp0 = (temp0.array() < 0).select(0.0, MatrixXd::Zero(temp.rows(),temp0.cols()).array()+1.0);
                    temp1 = (temp1.array() < 0).select(0.0, MatrixXd::Zero(temp.rows(),temp1.cols()).array()+1.0);
                    //
                    De.col(ij) = beta_0[ij] * temp.array();
                    De.col(ij+1) = beta_0[ij] * temp.array();
                    Dde.col(ij) = temp.array();
                    Dde.col(ij+1) = beta_0[ij] * (temp1.array()-temp0.array()) / 2/dint;
                    //
                    Ddde.col((ij+1)*(ij+2)/2+ij) = (temp1.array()-temp0.array()) / 2/dint;
                    Ddde.col((ij+1)*(ij+2)/2+ij+1) = beta_0[ij] * (temp1.array()-2*temp.array()+temp0.array()) / pow(dint,2);
                } else {
                    ;
                }
            } else {
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
    }
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
                            Rdd.block(0,ijk,df0.rows(),1) = Tdd0.block(0,ij,df0.rows(),1).array() * Te.block(0,0,df0.rows(),1).array();
                        } else {
                            if (ij<dose_num_tot){
                                Rdd.block(0,ijk,df0.rows(),1) = Ddde.block(0,ijk,df0.rows(),1).array() * T0.block(0,fir,df0.rows(),1).array();
                            } else {
                                Rdd.block(0,ijk,df0.rows(),1) = Tdd0.block(0,ij,df0.rows(),1).array() * T0.block(0,fir,df0.rows(),1).array();
                            }
                        }
                    } else {
                        if (ij<dose_num_tot){
                            Rdd.block(0,ijk,df0.rows(),1) = Ddde.block(0,ijk,df0.rows(),1).array() * Te.block(0,0,df0.rows(),1).array();
                        } else {
                            Rdd.block(0,ijk,df0.rows(),1) = Tdd0.block(0,ij,df0.rows(),1).array() * De.block(0,0,df0.rows(),1).array();
                        }
                    }
                } else {
                    if (fir!=0){
                        if ((ij==fir)||(jk==fir)){
                            if (ij<dose_num_tot){
                                Rdd.block(0,ijk,df0.rows(),1) = Dde.block(0,ij,df0.rows(),1).array() * Td0.block(0,jk,df0.rows(),1).array();
                            } else {
                                Rdd.block(0,ijk,df0.rows(),1) = Td0.block(0,ij,df0.rows(),1).array() * Td0.block(0,jk,df0.rows(),1).array();
                            }
                        } else if ((ij<dose_num_tot)&&(jk<dose_num_tot)){
                            Rdd.block(0,ijk,df0.rows(),1) = Ddde.block(0,ijk,df0.rows(),1).array() * Dde.block(0,jk,df0.rows(),1).array();
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
        Te = T0.array().block(0,dose_num_tot,T0.rows(),T0.cols()-dose_num_tot).rowwise().prod() * Dose.array();
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
                    Rdd.block(0,ijk,df0.rows(),1) = Ddde.block(0,ijk,df0.rows(),1).array() * Te.block(0,0,df0.rows(),1).array();
                } else {
                    Rdd.block(0,ijk,df0.rows(),1) = Tdd0.block(0,jk,df0.rows(),1).array() * T0.block(0,jk,df0.rows(),1).array().pow(-1).array() * R.block(0,0,df0.rows(),1).array();
                }
            } else {
                if ((ij<dose_num_tot)&&(jk<dose_num_tot)){
                    Rdd.block(0,ijk,df0.rows(),1) = Ddde.block(0,ijk,df0.rows(),1).array() * Te.block(0,0,df0.rows(),1).array();
                } else if ((ij<dose_num_tot)||(jk<dose_num_tot)){
                    Rdd.block(0,ijk,df0.rows(),1) = Dde.block(0,jk,df0.rows(),1).array() * De.block(0,jk,df0.rows(),1).array().pow(-1).array() * Rd.block(0,ij,df0.rows(),1).array();
                } else{
                    Rdd.block(0,ijk,df0.rows(),1) = Td0.block(0,jk,df0.rows(),1).array() * T0.block(0,jk,df0.rows(),1).array().pow(-1).array() * Rd.block(0,ij,df0.rows(),1).array();
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
            }
            MatrixXd Ld = MatrixXd::Zero(dj,4);
            Ld << R.block(RiskFail(j,0),0,dj,1), Rd.block(RiskFail(j,0),ij,dj,1), Rd.block(RiskFail(j,0),jk,dj,1) ,Rdd.block(RiskFail(j,0),ijk,dj,1);//sum of risks in group
            //
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
    //
    #pragma omp parallel for schedule(dynamic) num_threads(nthreads) reduction(vec_double_plus:Ll,Lld,Lldd) collapse(2)
    for (int ijk=0;ijk<totalnum*(totalnum+1)/2;ijk++){//totalnum*(totalnum+1)/2
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
            Ld << R.block(RiskFail(j,0),0,dj,1), RdR.block(RiskFail(j,0),ij,dj,1), RdR.block(RiskFail(j,0),jk,dj,1) ,RddR.block(RiskFail(j,0),ijk,dj,1);//sum of risks in group
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
            MatrixXd temp1 = MatrixXd::Zero(Ld.rows(),1);
            MatrixXd temp2 = MatrixXd::Zero(Ld.rows(),1);
            temp1 = Ld.col(0).array().log();
            double Ld1 =  (temp1.array().isFinite()).select(temp1,0).sum();
            temp1 = Ld.col(1).array();
            temp2 = Ld.col(2).array();
            double Ld2 = (temp1.array().isFinite()).select(temp1,0).sum();
            temp1 = Ld.col(3).array() - temp1.array() * temp2.array();
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
                if ((tform[ijk]=="step_ind")||(tform[ijk]=="lin_ind")){
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
            #pragma omp parallel for schedule(dynamic) num_threads(nthreads)
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
                    MatrixXd bch = ((beta_c[ijk] - beta_p[ijk]) * df0.col(ijk)).array().exp().array();
                    //
                    De.col(ijk) = ach * bch.array() * De.col(ijk).array();
                    De.col(ijk+1) = ach * bch.array() * De.col(ijk+1).array();
                    Dde.col(ijk) = bch.array() * Dde.col(ijk).array();
                    Dde.col(ijk+1) = ach * bch.array() * Dde.col(ijk+1).array();
                    Ddde.col((ijk+1)*(ijk+2)/2+ijk) = bch.array() * Ddde.col((ijk+1)*(ijk+2)/2+ijk).array();
                    Ddde.col((ijk+1)*(ijk+2)/2+ijk+1) = ach * bch.array() * Ddde.col((ijk+1)*(ijk+2)/2+ijk+1).array();
                    //
                } else if (tform[ijk]=="loglin_top"){
                    if (ijk==0){
                        beta_c[ijk] = beta_a[ijk] + Dbeta.coeff(ijk,ijk_ind);
                        //
                        De.col(ijk) = De.col(ijk).array() * ((beta_c[ijk] - beta_p[ijk]) * df0.col(ijk)).array().exp().array();
                        Dde.col(ijk) = Dde.col(ijk).array() * ((beta_c[ijk] - beta_p[ijk]) * df0.col(ijk)).array().exp().array();
                        Ddde.col((ijk)*(ijk+1)/2+ijk) = Ddde.col((ijk)*(ijk+1)/2+ijk).array() * ((beta_c[ijk] - beta_p[ijk]) * df0.col(ijk)).array().exp().array();
                        //
                    } else if (tform[ijk-1]!="loglin_slope"){
                        beta_c[ijk] = beta_a[ijk] + Dbeta.coeff(ijk,ijk_ind);
                        //
                        De.col(ijk) = De.col(ijk).array() * ((beta_c[ijk] - beta_p[ijk]) * df0.col(ijk)).array().exp().array();
                        Dde.col(ijk) = Dde.col(ijk).array() * ((beta_c[ijk] - beta_p[ijk]) * df0.col(ijk)).array().exp().array();
                        Ddde.col((ijk)*(ijk+1)/2+ijk) = Ddde.col((ijk)*(ijk+1)/2+ijk).array() * ((beta_c[ijk] - beta_p[ijk]) * df0.col(ijk)).array().exp().array();
                        //
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
                    ArrayXd temp2 = (temp1.array() < 0).select(0, temp1) - (temp0.array() < 0).select(0, temp0);
                    //
                    De.col(ijk) = beta_c[ijk] * (temp.array() < 0).select(0.0, temp);
                    De.col((ijk+1)) = beta_c[ijk] * (temp.array() < 0).select(0.0, temp);
                    Dde.col(ijk) = (temp.array() < 0).select(0.0, temp);
                    Dde.col((ijk+1)) = beta_c[ijk] * (temp2) / 2/dint;
                    //
                    Ddde.col((ijk+1)*(ijk+2)/2+ijk) = (temp2) / 2/dint;
                    Ddde.col((ijk+1)*(ijk+2)/2+ijk+1) = beta_c[ijk] * (temp2) / pow(dint,2);
                    //
                } else if (tform[ijk]=="quad_slope"){
                    beta_c[ijk] = beta_a[ijk] + Dbeta.coeff(ijk,ijk_ind);
                    De.col(ijk) = beta_c[ijk] * df0.col(ijk).array().square().array();
                    //
                } else if (tform[ijk]=="step_slope"){
                    beta_c[ijk] = beta_a[ijk] + Dbeta.coeff(ijk,ijk_ind);
                    beta_c[ijk+1] = beta_a[ijk+1] + Dbeta.coeff(ijk+1,ijk_ind);
                    
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
                    Ddde.col((ijk+1)*(ijk+2)/2+ijk) = (temp2) / 2/dint;
                    Ddde.col((ijk+1)*(ijk+2)/2+ijk+1) = beta_0[ijk] * (temp2) / pow(dint,2);
                    //
                } else {
                    ;
                }
                Theta(ijk, ijk_ind) = beta_c[ijk];
                //
            }
            //
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
                                    Rdd.block(0,ijk,df0.rows(),1) = Tdd0.block(0,ij,df0.rows(),1).array() * Te.block(0,0,df0.rows(),1).array();
                                } else {
                                    if (ij<dose_num_tot){
                                        Rdd.block(0,ijk,df0.rows(),1) = Ddde.block(0,ijk,df0.rows(),1).array() * T0.block(0,fir,df0.rows(),1).array();
                                    } else {
                                        Rdd.block(0,ijk,df0.rows(),1) = Tdd0.block(0,ij,df0.rows(),1).array() * T0.block(0,fir,df0.rows(),1).array();
                                    }
                                }
                            } else {
                                if (ij<dose_num_tot){
                                    Rdd.block(0,ijk,df0.rows(),1) = Ddde.block(0,ijk,df0.rows(),1).array() * Te.block(0,0,df0.rows(),1).array();
                                } else {
                                    Rdd.block(0,ijk,df0.rows(),1) = Tdd0.block(0,ij,df0.rows(),1).array() * De.block(0,0,df0.rows(),1).array();
                                }
                            }
                        } else {
                            if (fir!=0){
                                if ((ij==fir)||(jk==fir)){
                                    if (ij<dose_num_tot){
                                        Rdd.block(0,ijk,df0.rows(),1) = Dde.block(0,ij,df0.rows(),1).array() * Td0.block(0,jk,df0.rows(),1).array();
                                    } else {
                                        Rdd.block(0,ijk,df0.rows(),1) = Td0.block(0,ij,df0.rows(),1).array() * Td0.block(0,jk,df0.rows(),1).array();
                                    }
                                } else if ((ij<dose_num_tot)&&(jk<dose_num_tot)){
                                    Rdd.block(0,ijk,df0.rows(),1) = Ddde.block(0,ijk,df0.rows(),1).array() * Dde.block(0,jk,df0.rows(),1).array();
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
                Te = T0.array().block(0,dose_num_tot,T0.rows(),T0.cols()-dose_num_tot).rowwise().prod() * Dose.array();
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
                            Rdd.block(0,ijk,df0.rows(),1) = Ddde.block(0,ijk,df0.rows(),1).array() * De.block(0,jk,df0.rows(),1).array().pow(-1).array() * R.block(0,0,df0.rows(),1).array();
                        } else {
                            Rdd.block(0,ijk,df0.rows(),1) = Tdd0.block(0,jk,df0.rows(),1).array() * Td0.block(0,jk,df0.rows(),1).array().pow(-1).array() * Rd.block(0,ij,df0.rows(),1).array();
                        }
                    } else {
                        if ((ij<dose_num_tot)&&(jk<dose_num_tot)){
                            Rdd.block(0,ijk,df0.rows(),1) = Ddde.block(0,ijk,df0.rows(),1).array() * De.block(0,ij,df0.rows(),1).array().pow(-1).array() * R.block(0,0,df0.rows(),1).array();
                        } else if ((ij<dose_num_tot)||(jk<dose_num_tot)){
                            Rdd.block(0,ijk,df0.rows(),1) = Dde.block(0,jk,df0.rows(),1).array() * De.block(0,jk,df0.rows(),1).array().pow(-1).array() * Rd.block(0,ij,df0.rows(),1).array();
                        } else{
                            Rdd.block(0,ijk,df0.rows(),1) = Td0.block(0,jk,df0.rows(),1).array() * T0.block(0,jk,df0.rows(),1).array().pow(-1).array() * Rd.block(0,ij,df0.rows(),1).array();
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
                    }
                    MatrixXd Ld = MatrixXd::Zero(dj,4);
                    Ld << R.block(RiskFail(j,0),0,dj,1), Rd.block(RiskFail(j,0),ij,dj,1), Rd.block(RiskFail(j,0),jk,dj,1) ,Rdd.block(RiskFail(j,0),ijk,dj,1);//sum of risks in group
                    //
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
            //
            #pragma omp declare reduction(vec_double_plus : std::vector<double> : \
                std::transform(omp_out.begin(), omp_out.end(), omp_in.begin(), omp_out.begin(), std::plus<double>())) \
                initializer(omp_priv = omp_orig)
            #pragma omp parallel for schedule(dynamic) num_threads(nthreads) reduction(vec_double_plus:Ll,Lld,Lldd) collapse(2)
            for (int ijk=0;ijk<totalnum*(totalnum+1)/2;ijk++){//totalnum*(totalnum+1)/2
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
                    Ld << R.block(RiskFail(j,0),0,dj,1), RdR.block(RiskFail(j,0),ij,dj,1), RdR.block(RiskFail(j,0),jk,dj,1) ,RddR.block(RiskFail(j,0),ijk,dj,1);//sum of risks in group
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
                    MatrixXd temp1 = MatrixXd::Zero(Ld.rows(),1);
                    MatrixXd temp2 = MatrixXd::Zero(Ld.rows(),1);
                    temp1 = Ld.col(0).array().log();
                    double Ld1 =  (temp1.array().isFinite()).select(temp1,0).sum();
                    temp1 = Ld.col(1).array();
                    temp2 = Ld.col(2).array();
                    double Ld2 = (temp1.array().isFinite()).select(temp1,0).sum();
                    temp1 = Ld.col(3).array() - temp1.array() * temp2.array();
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
        #pragma omp parallel for schedule(dynamic) num_threads(nthreads)
        for (int ij=0;ij<totalnum;ij++){
            int ind0 = ij;
            if (ind0 < dose_num_tot){
                if (tform[ij]=="loglin_slope"){
                    ArrayXd temp = (beta_0[ij+1] * df0.col(ij)).array().exp();
                    ArrayXd temp1 = beta_0[ij] * temp;
                    //
                    //
                    De.col(ij) = temp1;
                    De.col(ij+1) = temp1;
                    Dde.col(ij) = temp.array();
                    Dde.col(ij+1) = temp1.array() * df0.col(ij).array();
                    Ddde.col((ij+1)*(ij+2)/2+ij) = temp.array() * df0.col(ij).array();
                    Ddde.col((ij+1)*(ij+2)/2+ij+1) = temp1.array() * df0.col(ij).array().square().array();
                    
                } else if (tform[ij]=="loglin_top"){
                    if (ij==0){
                        ArrayXd temp = (beta_0[ij] * df0.col(ij)).array().exp();
                        De.col(ij) = temp;
                        Dde.col(ij) = temp.array() * df0.col(ij).array();
                        Ddde.col(ij * (ij+1)/2.0 + ij) = temp.array() * df0.col(ij).array().square().array();
                    } else if (tform[ij-1]!="loglin_slope"){
                        ArrayXd temp = (beta_0[ij] * df0.col(ij)).array().exp();
                        De.col(ij) = temp;
                        Dde.col(ij) = temp.array() * df0.col(ij).array();
                        Ddde.col(ij * (ij+1)/2.0 + ij) = temp.array() * df0.col(ij).array().square().array();
                        //
                    } else {
                        ;
                    }
                } else if (tform[ij]=="lin_slope"){
                    ArrayXd temp = (df0.col(ij).array() - beta_0[ij]);
                    ArrayXd temp0 = (df0.col(ij).array() - beta_0[ij]-dint);
                    ArrayXd temp1 = (df0.col(ij).array() - beta_0[ij]+dint);
                    //
                    temp = (temp.array() < 0).select(0, temp);
                    temp0 = (temp0.array() < 0).select(0, temp0);
                    temp1 = (temp1.array() < 0).select(0, temp1);
                    //
                    De.col(ij) = beta_0[ij] * temp.array();
                    De.col(ij+1) = beta_0[ij] * temp.array();
                    Dde.col(ij) = temp.array();
                    Dde.col(ij+1) = beta_0[ij] * (temp1.array()-temp0.array()) / 2/dint;
                    //
                    Ddde.col((ij+1)*(ij+2)/2+ij) = (temp1.array()-temp0.array()) / 2/dint;
                    Ddde.col((ij+1)*(ij+2)/2+ij+1) = beta_0[ij] * (temp1.array()-2*temp.array()+temp0.array()) / pow(dint,2);
                } else if (tform[ij]=="quad_slope"){
                    ArrayXd temp = df0.col(ij).array().square();
                    //
                    De.col(ij) = beta_0[ij] * temp.array();
                    Dde.col(ij) = temp.array();
                } else if (tform[ij]=="step_slope"){
                    ArrayXd temp = (df0.col(ij).array() - beta_0[ij]);
                    ArrayXd temp0 = (df0.col(ij).array() - beta_0[ij]-dint);
                    ArrayXd temp1 = (df0.col(ij).array() - beta_0[ij]+dint);
                    //
                    temp = (temp.array() < 0).select(0.0, MatrixXd::Zero(temp.rows(),temp.cols()).array()+1.0);
                    temp0 = (temp0.array() < 0).select(0.0, MatrixXd::Zero(temp.rows(),temp0.cols()).array()+1.0);
                    temp1 = (temp1.array() < 0).select(0.0, MatrixXd::Zero(temp.rows(),temp1.cols()).array()+1.0);
                    //
                    De.col(ij) = beta_0[ij] * temp.array();
                    De.col(ij+1) = beta_0[ij] * temp.array();
                    Dde.col(ij) = temp.array();
                    Dde.col(ij+1) = beta_0[ij] * (temp1.array()-temp0.array()) / 2/dint;
                    //
                    Ddde.col((ij+1)*(ij+2)/2+ij) = (temp1.array()-temp0.array()) / 2/dint;
                    Ddde.col((ij+1)*(ij+2)/2+ij+1) = beta_0[ij] * (temp1.array()-2*temp.array()+temp0.array()) / pow(dint,2);
                } else {
                    ;
                }
            } else {
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
    }
    for (int ijk=0;ijk<totalnum;ijk++){
        bound_results(ijk,1) = Theta.coeff(ijk,ijk);
    }
    List res_list = List::create(_["beta_0"]=wrap(beta_0),_["wald_bounds"]=wrap(bound_results));
    //
    return res_list;
}

NumericMatrix Schoenfeld_PEANUT( VectorXd beta_linT,VectorXd beta_loglinT,VectorXd beta_plinT,MatrixXd df_lin,MatrixXd df_loglin,MatrixXd df_plin, MatrixXd df_dose,int fir,string modelform,int ntime, NumericVector include_bool,List beta_loglin_slopes, List beta_loglin_tops , List beta_lin_slopes , List beta_lin_ints , List beta_quads , List beta_step_slopes , List beta_step_ints, NumericMatrix df_groups, NumericVector tu, double dose_abs_max){
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
        //as<std::vector<double> >(a);
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
    int nthreads = Eigen::nbThreads()-1; //stores how many threads are allocated
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
    //
    VectorXd beta_0(totalnum);
    MatrixXd df0 = MatrixXd::Zero(df_lin.rows(), totalnum); // stores memory for the derivative term parameters and columns
    MatrixXd T0 = MatrixXd::Zero(df0.rows(), totalnum); //preallocates matrix for Derivative column terms
//    MatrixXd Td0 = MatrixXd::Zero(df0.rows(), totalnum); //preallocates matrix for Derivative column terms
//    MatrixXd Tdd0 = MatrixXd::Zero(df0.rows(), totalnum); //preallocates matrix for Derivative column terms
    //
//    MatrixXd Te = MatrixXd::Zero(df0.rows(), 1); //preallocates matrix for non-Derivative column terms
//    MatrixXd R = MatrixXd::Zero(df0.rows(), 1); //preallocates matrix for Risks and derivatives
//    MatrixXd Rd = MatrixXd::Zero(df0.rows(), totalnum); //preallocates matrix for Risks and derivatives
//    MatrixXd Rdd = MatrixXd::Zero(df0.rows(), totalnum*(totalnum+1)/2);
    //
    MatrixXd De = MatrixXd::Zero(df_dose.rows(),dose_num_tot);
//    MatrixXd Dde = MatrixXd::Zero(df_dose.rows(),dose_num_tot);
    /*
    for (int ijk=0;ijk<totalnum*(totalnum+1)/2;ijk++){
        int ij = 0;
        int jk = ijk;
        while (jk>ij){
            ij++;
            jk-=ij;
        }
    */
//    MatrixXd Ddde = MatrixXd::Zero(df_dose.rows(),dose_num_tot*(dose_num_tot+1)/2);
    double dint = dose_abs_max;
    int total_dose=0;
    //
    ArrayXd Dose(df_dose.rows(),1);
    //
    #pragma omp parallel for schedule(dynamic) num_threads(nthreads)
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

            if (ind0 < loglin_size){
                ArrayXd temp = (beta_loglin_top[ind0] * df_dose.col(ijk)).array().exp();
                ArrayXd temp1 = beta_loglin_slope[ind0] * temp;
                //
                if ((beta_loglin_slope[ind0]==1)&&(loglin_size=1)){
                    int ind1 = cumulative_dose_num[ijk]+ind0;
                    //
                    beta_0[ind1] = beta_loglin_top[ind0];
                    tform[ind1] = "loglin_top";
                    df0.col(ind1) = df_dose.col(ijk);
                    //
                    De.col(ind1) = temp1;
                    Dose = Dose + temp1.array();
//                    Dde.col(ind1) = temp1.array() * df_dose.col(ijk).array();
//                    Ddde.col(ind1 * (ind1+1)/2) = temp1.array() * df_dose.col(ijk).array().square().array();
                    //
                } else {
                    int ind1 = cumulative_dose_num[ijk]+2*ind0;
//                    int ind2 = ind1 + 1; 
                    //
                    beta_0[ind1] = beta_loglin_slope[ind0];
                    beta_0[(ind1 + 1)] = beta_loglin_top[ind0];
                    tform[ind1] = "loglin_slope";
                    tform[(ind1 + 1)] = "loglin_top";
                    df0.col(ind1) = df_dose.col(ijk);
                    df0.col((ind1 + 1)) = df_dose.col(ijk);
                    //
                    De.col(ind1) = temp1;
                    De.col((ind1 + 1)) = temp1;
                    Dose = Dose + temp1.array();
//                    Dde.col(ind1) = temp.array();
//                    Dde.col((ind1 + 1)) = temp1.array() * df_dose.col(ijk).array();
//                    Ddde.col((ind1 + 1) * ((ind1 + 1)+1)/2 + ind1) = temp.array() * df_dose.col(ijk).array();
//                    Ddde.col((ind1 + 1) * ((ind1 + 1)+1)/2 + (ind1 + 1)) = temp1.array() * df_dose.col(ijk).array().square().array();
                }
                
            } else if (ind0 < loglin_size + lin_size){
                int jk = ind0 - loglin_size;
                ArrayXd temp = (df_dose.col(ijk).array() - beta_lin_int[jk]);
                ArrayXd temp0 = (df_dose.col(ijk).array() - beta_lin_int[jk]-dint);
                ArrayXd temp1 = (df_dose.col(ijk).array() - beta_lin_int[jk]+dint);
                ArrayXd temp2 = (temp1.array() < 0).select(0, temp1) - (temp0.array() < 0).select(0, temp0);
                //
                int ind1 = cumulative_dose_num[ijk]+2*loglin_size + 2*jk;
//                int ind2 = ind1 + 1; 
                //
                beta_0[ind1] = beta_lin_slope[jk];
                beta_0[(ind1 + 1)] = beta_lin_int[jk];
                tform[ind1] = "lin_slope";
                tform[(ind1 + 1)] = "lin_ind";
                df0.col(ind1) = df_dose.col(ijk);
                df0.col((ind1 + 1)) = df_dose.col(ijk);
                //
                De.col(ind1) = beta_lin_slope[jk] * (temp.array() < 0).select(0.0, temp);
                De.col((ind1 + 1)) = beta_lin_slope[jk] * (temp.array() < 0).select(0.0, temp);
                Dose = Dose + De.col(ind1).array();
//                Dde.col(ind1) = (temp.array() < 0).select(0.0, temp);
//                Dde.col((ind1 + 1)) = beta_lin_slope[jk] * (temp2) / 2/dint;
                //
//                Ddde.col((ind1 + 1) * ((ind1 + 1)+1)/2 + ind1) = (temp2) / 2/dint;
//                Ddde.col((ind1 + 1) * ((ind1 + 1)+1)/2 + (ind1 + 1)) = beta_lin_slope[jk] * (temp2) / pow(dint,2);
            } else if (ind0 < loglin_size + lin_size + quad_size){
                int jk = ind0 - loglin_size - lin_size;
                ArrayXd temp = df_dose.col(ijk).array().square();
                int ind1 = cumulative_dose_num[ijk]+2*loglin_size + 2*lin_size+jk;
                //
                beta_0[ind1] = beta_quad[jk];
                tform[ind1] = "quad_slope";
                df0.col(ind1) = df_dose.col(ijk);
                //
//                De.col(ind1) = beta_quad[jk] * temp.array();
//                Dde.col(ind1) = temp.array();
                Dose = Dose + De.col(ind1).array();
            } else {
                int jk = ind0 - loglin_size - lin_size - quad_size;
                ArrayXd temp = (df_dose.col(ijk).array() - beta_step_int[jk]);
                ArrayXd temp0 = (df_dose.col(ijk).array() - beta_step_int[jk]-dint);
                ArrayXd temp1 = (df_dose.col(ijk).array() - beta_step_int[jk]+dint);
                ArrayXd temp2 = (temp1.array() < 0).select(0.0, MatrixXd::Zero(temp.rows(),temp.cols()).array()+1.0) - (temp0.array() < 0).select(0.0, MatrixXd::Zero(temp.rows(),temp.cols()).array()+1.0);
                //
                int ind1 = cumulative_dose_num[ijk]+2*loglin_size + 2*lin_size + quad_size + 2*jk;
//                int ind2 = ind1 + 1;
                //
                beta_0[ind1] = beta_step_slope[jk];
                beta_0[(ind1 + 1)] = beta_step_int[jk];
                tform[ind1] = "step_slope";
                tform[(ind1 + 1)] = "step_ind";
                df0.col(ind1) = df_dose.col(ijk);
                df0.col((ind1 + 1)) = df_dose.col(ijk);
                //
                De.col(ind1) = beta_step_slope[jk] * (temp.array() < 0).select(0.0, MatrixXd::Zero(temp.rows(),temp.cols()).array()+1.0);
                De.col((ind1 + 1)) = beta_step_slope[jk] * (temp.array() < 0).select(0.0, MatrixXd::Zero(temp.rows(),temp.cols()).array()+1.0);
//                Dde.col(ind1) = (temp.array() < 0).select(0.0, MatrixXd::Zero(temp.rows(),temp.cols()).array()+1.0);
//                Dde.col((ind1 + 1)) = beta_step_slope[jk] * (temp2) / 2/dint;
                //
//                Ddde.col((ind1 + 1) * ((ind1 + 1)+1)/2 + ind1) = (temp2) / 2/dint;
//                Ddde.col((ind1 + 1) * ((ind1 + 1)+1)/2 + (ind1 + 1)) = beta_step_slope[jk] * (temp2) / pow(dint,2);
                Dose = Dose + De.col(ind1).array();
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
                ;
//                Td0.col(ind1) = df0.col(ind1);
            } else if (tform[ind1]=="loglin") {
                T0.col(ind1) = T0.col(ind1).array().exp();
//                Td0.col(ind1) = df0.col(ind1).array() * T0.col(ind1).array();
//                Tdd0.col(ind1) = df0.col(ind1).array() * Td0.col(ind1).array();
            } else if (tform[ind1]=="plin") {
                T0.col(ind1) = 1 + T0.col(ind1).array();
//                Td0.col(ind1) = df0.col(ind1);
            } else {
                cout << tform[ind1] << " is invalid" << endl;
                throw invalid_argument( "Invalid term type" );
            }
        }
    }
    //
    T0.block(0,0,T0.rows(),dose_num_tot) = De.block(0,0,T0.rows(),dose_num_tot);
    //
    // -------------------------------------------------------------------------------------------
    //
    vector<string>  RiskGroup(ntime); //vector of strings detailing the rows
    IntegerMatrix RiskFail(ntime,2); //vector giving the event rows
    //
    const Map<MatrixXd> df_m(as<Map<MatrixXd> >(df_groups));
    #pragma omp parallel for schedule(dynamic) num_threads(nthreads)
    for (int ijk=0;ijk<ntime;ijk++){
        double t0 = tu[ijk];
        VectorXi select_ind_all = ((df_m.col(0).array() < t0)&&(df_m.col(1).array()>=t0)).cast<int>();
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
    //
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
            double x2_expect=0;
            for (int i = 0; i < InGroup.size()-1; i=i+2){
                t_sum += T0.block(InGroup[i]-1,ijk,InGroup[i+1]-InGroup[i]+1,1).sum();
                x_expect += (df0.block(InGroup[i]-1,ijk,InGroup[i+1]-InGroup[i]+1,1).array() * T0.block(InGroup[i]-1,ijk,InGroup[i+1]-InGroup[i]+1,1).array()).sum();
                x2_expect += (df0.block(InGroup[i]-1,ijk,InGroup[i+1]-InGroup[i]+1,1).array().square().array() * T0.block(InGroup[i]-1,ijk,InGroup[i+1]-InGroup[i]+1,1).array()).sum();
            }
            x_expect = x_expect / t_sum;
            //Now has the grouping pairs
            int dj = RiskFail(j,1)-RiskFail(j,0)+1;
            double x_risks = df0.block(RiskFail(j,0),ijk,dj,1).sum()/dj;
            //
            x2_expect = x2_expect / t_sum - pow(x_expect,2);
            if (x2_expect==0){
                x2_expect=1;
            }
            residuals(j,ijk) = (x_risks - x_expect)/x2_expect;
        }
    }
    cout << residuals.block(1,1,5,5) << endl;
    return wrap(residuals);
}








