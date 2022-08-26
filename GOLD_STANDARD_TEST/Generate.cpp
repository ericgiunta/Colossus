#include <RcppEigen.h>
#include <RcppParallel.h>
#include <omp.h>
#include "Generate.h"
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
using namespace R;
using namespace Eigen;

using Eigen::Map;
using Eigen::MatrixXd;
using Eigen::SparseMatrix;
using Eigen::VectorXd;
using Rcpp::as;

// [[Rcpp::export]]
void RunGenerator(string fname, NumericVector Pbs, int nm_xtra_rows, NumericVector Pcs, StringVector RowNames,int ngoal){
    Map<VectorXd> Pc(as<Map<VectorXd> >(Pcs));
    Map<VectorXd> Pb(as<Map<VectorXd> >(Pbs));
    cout << Pc.transpose() << endl;
    cout << Pb.transpose() << endl;
    //
//    cout << R::ppois(70+1,Pb[0],true,false) << " " << R::ppois(70,Pb[0],true,false) << endl;
    //
    cout.precision(10); //forces higher precision numbers printed to terminal
    int nthreads = Eigen::nbThreads()-1; //stores how many threads are allocated
    //
    class People {
        public:
            int YOB;
            int age_entry;
            int sexm;
            int edu;
            int duration;
            double cumulative_dose=0.0;
            int lung=0;
            int base_lung=0;
            int dead=0;
            double risk=1.0;
            double base_risk=1.0;
            int ID;
            vector <double> lung_x;
        void set_person( int nm_xtra_rows,int index){
            unsigned seed1 = std::chrono::system_clock::now().time_since_epoch().count();
            mt19937_64 generator(seed1);
            uniform_int_distribution<int> distribution_30(0,30);
            uniform_int_distribution<int> distribution_20(0,20);
            uniform_int_distribution<int> distribution_5(0,5);
            uniform_int_distribution<int> distribution_1_int(0,1);
            uniform_real_distribution<double> distribution_1_real(0.0,1.0);
            //
            ID=index;
            YOB=distribution_30(generator) + 1920;
            age_entry=distribution_20(generator) + 30;
            sexm=distribution_1_int(generator);
            edu=distribution_1_int(generator);
            lung_x.resize(nm_xtra_rows);
            cumulative_dose = distribution_1_real(generator) * 1400.0 * ((double)age_entry-20.0)/(100.0-20.0);
            //
            auto gen = [&distribution_1_real, &generator](){
               return distribution_1_real(generator);
            };
            //
            generate(lung_x.begin(), lung_x.end(), gen);
        }
        void RollRisk(VectorXd cov_vals,VectorXd base_vals,VectorXd para_vals, int age_change){
            unsigned seed1 = std::chrono::system_clock::now().time_since_epoch().count();
            mt19937_64 generator(seed1);
            uniform_real_distribution<double> distribution_1_real(0.0,1.0);
            duration=age_change;
            base_risk=BaseHazard_CDF(age_entry,base_vals[0],age_change) * (1.0 + 0.1*(distribution_1_real(generator)-0.5));
            risk=HazardRatio(cov_vals,para_vals) * (1.0 + 0.2*(distribution_1_real(generator)-0.5));
            double risk_comp = distribution_1_real(generator);
            if (risk*base_risk>risk_comp){
                lung= 1;
                dead=1;
            }else{
                lung= 0;
            }
            if (base_risk>risk_comp){
                base_lung= 1;
                dead=1;
            }else{
                base_lung= 0;
            }
        }
        void UpdateSubject(){
            unsigned seed1 = std::chrono::system_clock::now().time_since_epoch().count();
            mt19937_64 generator(seed1);
            uniform_real_distribution<double> distribution_1_real(0.0,1.0);
            age_entry = age_entry + duration;
            cumulative_dose = cumulative_dose + distribution_1_real(generator)*20.0;
        }
        void RollDeath(double pd){
            unsigned seed1 = std::chrono::system_clock::now().time_since_epoch().count();
            mt19937_64 generator(seed1);
            uniform_real_distribution<double> distribution_1_real(0.0,1.0);
            if (distribution_1_real(generator)<pd){
                dead=1;
            }
        }
    };
    int all_col = Pc.size()+5;
    int indv_tot = 100;
    int indv_time = 50;
    int running_tot =0;
//    int ngoal = 1e5;
    const static IOFormat CSVFormat(FullPrecision, DontAlignCols, ", ", "\n");
    ofstream file(fname);
    for (int i=0;i<RowNames.size()-1;i++){
        file << RowNames[i] << ",";
    }
    file << RowNames[RowNames.size()-1] << endl;
    MatrixXd Subject_Record = MatrixXd::Zero(indv_tot*indv_time,all_col);
//    cout << Subject_Record.sum() << endl;
//    vector <int> end_locs(indv_tot,0);
    VectorXi keep_row;
    MatrixXd mat_sel;
    while (running_tot<ngoal){
//        cout << running_tot << endl;
//        fill(end_locs.begin(),end_locs.end(),0);
        #pragma omp parallel for schedule(dynamic) num_threads(nthreads)
        for (int j=0;j<indv_tot;j++){
            VectorXd cov_vec = VectorXd::Zero(Pc.size());
            VectorXd add_vec = VectorXd::Zero(all_col);
            People subject;
            subject.set_person(nm_xtra_rows,running_tot+j);
            unsigned seed1 = std::chrono::system_clock::now().time_since_epoch().count();
            mt19937_64 generator(seed1);
            uniform_real_distribution<double> distribution_1_real(0.0,1.0);
            int k=0;
            while (k<indv_time){
                //
                k++;
                cov_vec = VectorXd::Zero(Pc.size());
                add_vec = VectorXd::Zero(all_col);
                //
                cov_vec(0) = subject.YOB;
                cov_vec(1) = subject.age_entry;
                cov_vec(2) = subject.sexm;
                cov_vec(3) = subject.edu;
                cov_vec(4) = subject.duration;
                cov_vec(5) = subject.cumulative_dose;
                vector<double> xtra_lung = subject.lung_x;
                for (int i=0;i<nm_xtra_rows;i++){
                    cov_vec(6+i) = xtra_lung[i];
                }
                //
                subject.RollRisk(cov_vec,Pb,Pc,1);
                //
                add_vec.block(1,0,cov_vec.size(),1) << cov_vec;
                add_vec(0) = subject.ID+1;
                add_vec(7+nm_xtra_rows) = subject.base_risk;
                add_vec(7+nm_xtra_rows+1) = subject.risk;
                add_vec(7+nm_xtra_rows+2) = subject.lung;
                add_vec(7+nm_xtra_rows+3) = subject.base_lung;
                Subject_Record.row(j*indv_time+k) << add_vec.transpose();
//                cout << k << " " << subject.age_entry << endl;
                //
                if (subject.age_entry>100){
                    subject.RollDeath(0.2);
                }
                if (subject.dead==1){
                    break;
                }
                subject.UpdateSubject();
            }
//            end_locs[j]=k;
        }
        //
        keep_row = (Subject_Record.col(0).array() > 0.0).cast<int>();
        mat_sel.resize(keep_row.sum(), Subject_Record.cols());
        int cur_row = 0;
        for (int i=0; i<Subject_Record.rows(); ++i) {
            if (keep_row[i]) {       
               mat_sel.row(cur_row) = Subject_Record.row(i);
               cur_row++;
            }
        }
        //
        file << mat_sel.format(CSVFormat);
        file << " " << endl;
        running_tot+=mat_sel.rows();
        Subject_Record = MatrixXd::Zero(indv_tot*indv_time,all_col);
    }
    file.close();
}



double BaseHazard_CDF(int t, int a, int d){
    return R::ppois(t+1,a,true,false)-R::ppois(t,a,true,false);
}

double HazardRatio(VectorXd C, VectorXd Pc){
    return exp((C.array() * Pc.array()).sum());
}
