#include <RcppEigen.h>
#include <RcppParallel.h>
#include <omp.h>
#include "Calc_Repeated.h"
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

//' Utility function to take the norm of a vector
//' \code{vec_norm} Called within the code, takes a vector and length and finds the norm
//' @param     x    std::vector to take norm of, assumed doubles
//' @param     totalnum    vector length
//'
//' @return double of the norm
// [[Rcpp::export]]
double vec_norm(const vector<double>& x,int totalnum){
  double total = 0;
  for (int i = 0; i < totalnum; i++)
  {
    total += x[i] * x[i];
  }
  return sqrt(total);
}

//' Utility function to remove rows
//' \code{removeRow} Called within the code, removes and resizes the matrix
//' @param     matrix    Eigen matrix to change
//' @param     rowToRemove    row index to remove
//'
//' @return eigen matrix with the row removed
// [[Rcpp::export]]
void removeRow(Eigen::MatrixXd& matrix, unsigned int rowToRemove){
    //
    //Used to resize with removed rows
    unsigned int numRows = matrix.rows()-1;
    unsigned int numCols = matrix.cols();

    if( rowToRemove < numRows )
        matrix.block(rowToRemove,0,numRows-rowToRemove,numCols) = matrix.block(rowToRemove+1,0,numRows-rowToRemove,numCols);

    matrix.conservativeResize(numRows,numCols);
}

//' Utility function to remove rows
//' \code{removeColumn} Called within the code, removes and resizes the matrix
//' @param     matrix    Eigen matrix to change
//' @param     colToRemove    column index to remove
//'
//' @return eigen matrix with the column removed
// [[Rcpp::export]]
void removeColumn(Eigen::MatrixXd& matrix, unsigned int colToRemove){
    //
    //Used to resize with removed columns
    unsigned int numRows = matrix.rows();
    unsigned int numCols = matrix.cols()-1;

    if( colToRemove < numCols )
        matrix.block(0,colToRemove,numRows,numCols-colToRemove) = matrix.block(0,colToRemove+1,numRows,numCols-colToRemove);

    matrix.conservativeResize(numRows,numCols);
}


//' Utility function to calculate the term and subterm values
//' \code{Make_Subterms} Called to update term matrices, Uses lists of term numbers and types to apply formulas
//' @param     totalnum    Total number of terms
//' @param     Term_n    Term numbers
//' @param     tform    subterm types
//' @param     dfc    covariate column numbers
//' @param     fir    first term number
//' @param     T0    Term by subterm matrix
//' @param     Td0    Term by subterm derivative matrix
//' @param     Tdd0    Term by subterm second derivative matrix
//' @param     Dose    Dose term matrix
//' @param     nonDose    nonDose term matrix
//' @param     TTerm    Total term matrix
//' @param     nonDose_LIN    Linear term matrix
//' @param     nonDose_PLIN    Product linear term matrix
//' @param     nonDose_LOGLIN    Loglinear term matrix
//' @param     beta_0    parameter list
//' @param     df0    covariate matrix
//' @param     dint    value used for threshold derivative finite step
//' @param     nthreads    number of threads to use
//' @param     debugging    debugging boolean
//'
//' @return Updates matrices in place: Sub-term matrices, Term matrices
// [[Rcpp::export]]
void Make_Subterms(const int& totalnum, const IntegerVector& Term_n,const StringVector&  tform, const IntegerVector& dfc, const int& fir, MatrixXd& T0, MatrixXd& Td0, MatrixXd& Tdd0, MatrixXd& Dose, MatrixXd& nonDose,  MatrixXd& TTerm, MatrixXd& nonDose_LIN, MatrixXd& nonDose_PLIN, MatrixXd& nonDose_LOGLIN,const  VectorXd& beta_0,const  MatrixXd& df0,const double& dint, const int& nthreads, bool debugging){
    //
    //Make_Subterms( totalnum, dose_num_tot, dose_term_tot, dose_breaks, beta_loglin_slopes_CPP, beta_loglin_tops_CPP, beta_lin_slopes_CPP, beta_lin_ints_CPP, beta_quads_CPP, beta_step_slopes_CPP, beta_step_ints_CPP, beta_lin, beta_loglin, beta_plin, df_lin, df_loglin, df_plin, df_dose, De, Dde, Ddde, T0, Td0, Tdd0, Dose,cumulative_dose_num,beta_0, df0,dint,nthreads, tform,include_bool, debugging);
    //
    // Calculates the sub term values
    //
    vector<int> lin_count(nonDose.cols(),0);
    vector<int> dose_count(nonDose.cols(),0);
    #pragma omp declare reduction (eig_plus: MatrixXd: omp_out=omp_out.array() + omp_in.array()) initializer(omp_priv=MatrixXd::Constant(omp_orig.rows(),omp_orig.cols(),0.0))
    #pragma omp declare reduction (eig_mult: MatrixXd: omp_out=omp_out.array() * omp_in.array()) initializer(omp_priv=MatrixXd::Constant(omp_orig.rows(),omp_orig.cols(),1.0))
    #pragma omp declare reduction(vec_int_plus : std::vector<int> : \
            std::transform(omp_out.begin(), omp_out.end(), omp_in.begin(), omp_out.begin(), std::plus<int>())) \
            initializer(omp_priv = omp_orig)
    #pragma omp parallel for schedule(dynamic) num_threads(nthreads) reduction(eig_plus:Dose,nonDose_LIN,nonDose_PLIN) reduction(eig_mult:nonDose_LOGLIN) reduction(vec_int_plus:lin_count,dose_count)
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
            dose_count[tn]=dose_count[tn]+1;
            
        } else if (as< string>(tform[ij])=="loglin_top"){
            if (ij==0){
                ArrayXd temp = (beta_0[ij] * df0.col(df0_c)).array().exp();
                T0.col(ij) = temp;
                Dose.col(tn) = Dose.col(tn).array() + T0.col(ij).array();
                dose_count[tn]=dose_count[tn]+1;

            } else if (tform[ij-1]!="loglin_slope"){
                ArrayXd temp = (beta_0[ij] * df0.col(df0_c)).array().exp();
                T0.col(ij) = temp;
                Dose.col(tn) = Dose.col(tn).array() + T0.col(ij).array();
                dose_count[tn]=dose_count[tn]+1;
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
            dose_count[tn]=dose_count[tn]+1;

        } else if (as< string>(tform[ij])=="quad_slope"){
            ArrayXd temp = df0.col(df0_c).array().square();
            //
            T0.col(ij) = beta_0[ij] * temp.array();
            Dose.col(tn) = Dose.col(tn).array() + T0.col(ij).array();
            dose_count[tn]=dose_count[tn]+1;
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
            dose_count[tn]=dose_count[tn]+1;

        } else if (as< string>(tform[ij])=="lin") {
            T0.col(ij) = (df0.col(df0_c).array() * beta_0[ij]).matrix();
            nonDose_LIN.col(tn) = nonDose_LIN.col(tn).array() + T0.col(ij).array();
            lin_count[tn]=lin_count[tn]+1;

        } else if (as< string>(tform[ij])=="loglin") {
            T0.col(ij) = (df0.col(df0_c).array() * beta_0[ij]).matrix();
            T0.col(ij) = T0.col(ij).array().exp();;
            nonDose_LOGLIN.col(tn) = nonDose_LOGLIN.col(tn).array() * T0.col(ij).array();

        } else if (as< string>(tform[ij])=="plin") {
            T0.col(ij) = (df0.col(df0_c).array() * beta_0[ij]).matrix();
            T0.col(ij) = 1 + T0.col(ij).array();
            nonDose_PLIN.col(tn) = nonDose_PLIN.col(tn).array() + T0.col(ij).array();

        } else {
            throw invalid_argument( "incorrect subterm type" );
        }
    }
    //
    // Calculates the terms and derivatives
    //
    //
    //
//    nonDose_LOGLIN = (nonDose_LOGLIN.array() != 1).select(nonDose_LOGLIN,1.0);
//    nonDose_PLIN = (nonDose_PLIN.array() != 0).select(nonDose_PLIN,1.0);
    for (int ijk=0; ijk<nonDose.cols();ijk++){ //combines non-dose terms into a single term
//        Rcout << dose_count[ijk] << " " << lin_count[ijk] << endl;
        if (dose_count[ijk]==0){
            Dose.col(ijk) = Dose.col(ijk).array() * 0.0 + 1;
        }
        if (lin_count[ijk]==0){
            nonDose_LIN.col(ijk) = nonDose_LIN.col(ijk).array() * 0.0 + 1;//replaces missing data with 1
        }
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
            throw invalid_argument( "incorrect subterm type" );
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
                    Tdd0.col((ij)*(ij+1)/2+ij) = df0.col(df0_ij).array().pow(2).array() * nonDose_LOGLIN.col(tij).array();
                } else if (as< string>(tform[jk])=="loglin") {
                    Tdd0.col((ij)*(ij+1)/2+jk) = df0.col(df0_ij).array() * df0.col(df0_jk).array() * nonDose_LOGLIN.col(tij).array();
                }
            }
        }
    }
    return;
}

//' Utility function to calculate the term and subterm values with the basic model
//' \code{Make_Subterms_Basic} Called to update term matrices, Uses lists of term numbers and types to apply formulas
//' @param     totalnum    Total number of terms
//' @param     dfc    covariate column numbers
//' @param     T0    subterm values
//' @param     beta_0    parameter list
//' @param     df0    covariate matrix
//' @param     nthreads    number of threads to use
//' @param     debugging    debugging boolean
//'
//' @return Updates matrices in place: Sub-term matrices, Term matrices
// [[Rcpp::export]]
void Make_Subterms_Basic(const int& totalnum, const IntegerVector& dfc, MatrixXd& T0, const VectorXd& beta_0,const MatrixXd& df0, const int& nthreads, bool debugging){
    //
    //Make_Subterms( totalnum, dose_num_tot, dose_term_tot, dose_breaks, beta_loglin_slopes_CPP, beta_loglin_tops_CPP, beta_lin_slopes_CPP, beta_lin_ints_CPP, beta_quads_CPP, beta_step_slopes_CPP, beta_step_ints_CPP, beta_lin, beta_loglin, beta_plin, df_lin, df_loglin, df_plin, df_dose, De, Dde, Ddde, T0, Td0, Tdd0, Dose,cumulative_dose_num,beta_0, df0,dint,nthreads, tform,include_bool, debugging);
    //
    // Calculates the sub term values
    //

    #pragma omp parallel for schedule(dynamic) num_threads(nthreads)
    for (int ij=0;ij<(totalnum);ij++){
        int df0_c = dfc[ij]-1;
        T0.col(ij) = (df0.col(df0_c).array() * beta_0[ij]).matrix();
        T0.col(ij) = T0.col(ij).array().exp();
    }
    return;
}


//' Utility function to calculate the risk and risk ratios
//' \code{Make_Risks} Called to update risk matrices, Splits into cases based on model form, Uses lists of term numbers and types to apply different derivative formulas    
//' @param     modelform    Model string
//' @param     tform    subterm types
//' @param     Term_n    term numbers
//' @param     totalnum    total number of terms
//' @param     fir    first term number
//' @param     T0    Term by subterm matrix
//' @param     Td0    Term by subterm derivative matrix
//' @param     Tdd0    Term by subterm second derivative matrix
//' @param     Te    Temporary term storage matrix
//' @param     R    Risk matrix
//' @param     Rd    Risk first derivative matrix
//' @param     Rdd    Risk second derivative matrix
//' @param     Dose    Dose term matrix
//' @param     nonDose    nonDose term matrix
//' @param     TTerm    Total term matrix
//' @param     nonDose_LIN    Linear term matrix
//' @param     nonDose_PLIN    Product linear term matrix
//' @param     nonDose_LOGLIN    Loglinear term matrix
//' @param     RdR    Risk to first derivative ratio matrix
//' @param     RddR    Risk to second derivative ratio matrix
//' @param     nthreads    number of threads to use
//' @param     debugging    debugging boolean
//'
//' @return Updates matrices in place: Risk, Risk ratios
// [[Rcpp::export]]
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
        MatrixXd TTerm_p = MatrixXd::Zero(TTerm.rows(),TTerm.cols());
        TTerm_p << TTerm.array() + 1.0;
        TTerm_p.col(fir) = TTerm.col(fir).array();
        Te = TTerm_p.array().rowwise().prod().array();
        R << Te.array();
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

//' Utility function to calculate the risk and risk ratios for the basic model
//' \code{Make_Risks_Basic} Called to update risk matrices, Splits into cases based on model form, Uses lists of term numbers and types to apply different derivative formulas    
//' @param     totalnum    total number of terms
//' @param     T0    Term by subterm matrix
//' @param     R    Risk matrix
//' @param     Rd    Risk first derivative matrix
//' @param     Rdd    Risk second derivative matrix
//' @param     RdR    Risk to first derivative ratio matrix
//' @param     RddR    Risk to second derivative ratio matrix
//' @param     nthreads    number of threads to use
//' @param     debugging    debugging boolean
//' @param     df0    covariate matrix
//' @param     dfc    covariate column numbers
//'
//' @return Updates matrices in place: Risk, Risk ratios
// [[Rcpp::export]]
void Make_Risks_Basic(const int& totalnum, const MatrixXd& T0, MatrixXd& R, MatrixXd& Rd, MatrixXd& Rdd, MatrixXd& RdR, MatrixXd& RddR, const int& nthreads, bool debugging,const MatrixXd& df0, const IntegerVector& dfc){
    //
    //Make_Risks(modelform, tform, Term_n, totalnum, fir, T0, Td0, Tdd0, Te, R, Rd, Rdd, Dose, nonDose, RdR, RddR, nthreads, debugging);
    //
    //
    //
    
    #pragma omp parallel for schedule(dynamic) num_threads(nthreads)
    for (int ijk=0;ijk<totalnum;ijk++){
        int df0_c = dfc[ijk]-1;
        Rd.col(ijk) = R.col(0).array() * df0.col(df0_c).array() ;
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
        int df0_c = dfc[ij]-1;
        Rdd.col(ijk) = Rd.col(jk).array() * df0.col(df0_c).array();
    }
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
        int df0_ij = dfc[ij]-1;
        int df0_jk = dfc[jk]-1;
        if (ij==jk){
            RdR.col(ij)=df0.col(df0_ij).array();
        }
        RddR.col(ijk)=df0.col(df0_ij).array() * df0.col(df0_jk).array();
    }
    return;
}

//' Utility function to define risk groups
//' \code{Make_Groups} Called to update lists of risk groups, Uses list of event times and row time/event information, Matrices store starting/stopping row indices for each group    
//' @param     ntime    number of event times
//' @param     df_m    event/time matrix
//' @param     RiskFail    Matrix of event rows for each event time
//' @param     RiskGroup    vectors of strings with rows at risk for each event time
//' @param     tu    Event time vector
//' @param     nthreads    number of threads
//' @param     debugging    debugging boolean
//'
//' @return Updates matrices in place: Matrix of event rows for each event time, vectors of strings with rows at risk for each event time
// [[Rcpp::export]]
void Make_Groups(const int& ntime, const MatrixXd& df_m, IntegerMatrix& RiskFail, vector<string>&  RiskGroup,  NumericVector& tu, const int& nthreads, bool debugging ){
    //
    //Make_Subterms( ntime, df_m, RiskFail, RiskGroup, tu, nthreads, debugging)
    //
    if (debugging){
        Rcout << "Starting Debug" << endl;
        vector<double> time_ref(6,0.0);
        vector<double> time_refs(6,0.0);
        Rcout << time_ref.size() << " " << time_refs.size() << endl;
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
            VectorXi select_ind_all = (((df_m.col(0).array() < t0)||(df_m.col(0).array()==df_m.col(1).array()))&&(df_m.col(1).array()>=t0)).cast<int>(); //indices at risk
            vector<int> indices_all;
            //
            end_point = system_clock::now();
            ending = time_point_cast<microseconds>(end_point).time_since_epoch().count();
            time_ref[0]+=ending-start;
            time_refs[0]+=pow(ending-start,2);
            start=ending;
            //
            VectorXi select_ind_end = ((df_m.col(2).array() == 1)&&(((df_m.col(1).array()==t0))||(df_m.col(0).array()==t0))).cast<int>(); //indices with events
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
        Rcout << "df204 ";
        for (vector<double>::size_type ijk=0;ijk<time_ref.size();ijk++){ 
            Rcout << (time_ref[ijk]/ntime)*1e-6 << " ";
        }
        Rcout << " " << endl;
        Rcout << "df205 ";
        for (vector<double>::size_type ijk=0;ijk<time_ref.size();ijk++){ 
            Rcout <<  sqrt(time_refs[ijk]/ntime - pow(time_ref[ijk]/ntime,2))*1e-6 << " ";
        }
        Rcout << " " << endl;
    } else {
        #pragma omp parallel for schedule(dynamic) num_threads(nthreads)
        for (int ijk=0;ijk<ntime;ijk++){
            double t0 = tu[ijk];
            VectorXi select_ind_all = (((df_m.col(0).array() < t0)||(df_m.col(0).array()==df_m.col(1).array()))&&(df_m.col(1).array()>=t0)).cast<int>(); //indices at risk
            vector<int> indices_all;
            VectorXi select_ind_end = ((df_m.col(2).array() == 1)&&(df_m.col(1).array()==t0)).cast<int>(); //indices with events
            vector<int> indices_end;
            //
            // for (int i=0;i<df_m.rows();i++){
            //     Rcout << df_m.row(i) << " " << select_ind_all[i] << " " << select_ind_end[i] << " " << t0 << endl;
            // }
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
            // Rcout << oss.str() << ":" << RiskGroup[ijk] << endl;
        }
    }
    return;
}

//' Utility function to define risk groups with STRATA
//' \code{Make_Groups_STRATA} Called to update lists of risk groups, Uses list of event times and row time/event information, Matrices store starting/stopping row indices for each group    
//' @param     ntime    number of event times
//' @param     df_m    event/time matrix
//' @param     RiskFail    Matrix of event rows for each event time
//' @param     RiskGroup    vectors of strings with rows at risk for each event time
//' @param     tu    Event time vector
//' @param     nthreads    number of threads
//' @param     debugging    debugging boolean
//' @param     STRATA_vals vector of strata identifier values
//'
//' @return Updates matrices in place: Matrix of event rows for each event time, vectors of strings with rows at risk for each event time
// [[Rcpp::export]]
void Make_Groups_STRATA(const int& ntime, const MatrixXd& df_m, IntegerMatrix& RiskFail, StringMatrix&  RiskGroup,  NumericVector& tu, const int& nthreads, bool debugging, IntegerVector& STRATA_vals){
    //
    //Make_Subterms( ntime, df_m, RiskFail, RiskGroup, tu, nthreads, debugging)
    //
    //
    vector<vector<int>> safe_fail(ntime);
    vector<vector<string>> safe_group(ntime);
    for (int i=0;i<ntime;i++){
        safe_fail[i] = vector<int>(RiskFail.cols(),0);
        safe_group[i] = vector<string>(RiskGroup.cols(),"");
    }
    //
    #pragma omp parallel for schedule(dynamic) num_threads(nthreads) collapse(2)
    for (int s_ij=0;s_ij<STRATA_vals.size();s_ij++){
        for (int ijk=0;ijk<ntime;ijk++){
            double t0 = tu[ijk];
            VectorXi select_ind_all = (((df_m.col(0).array() < t0)||(df_m.col(0).array()==df_m.col(1).array()))&&(df_m.col(1).array()>=t0)&&(df_m.col(3).array()==STRATA_vals[s_ij])).cast<int>(); //indices at risk
            vector<int> indices_all;
            VectorXi select_ind_end = ((df_m.col(2).array() == 1)&&(df_m.col(1).array()==t0)&&(df_m.col(3).array()==STRATA_vals[s_ij])).cast<int>(); //indices with events
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
            if (indices_end.size()>0){
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
//                Rcout << s_ij << " " << ijk << endl;
                safe_fail[ijk][2*s_ij+0] = indices_end[0]-1;//Due to the sorting method, there is a continuous block of event rows
                safe_fail[ijk][2*s_ij+1] = indices_end[indices_end.size()-1]-1;
//                RiskFail(ijk,2*s_ij + 0)=indices_end[0]-1;//Due to the sorting method, there is a continuous block of event rows
//                RiskFail(ijk,2*s_ij + 1)=indices_end[indices_end.size()-1]-1;
                //
                ostringstream oss;
                copy(indices.begin(), indices.end(),
                    std::ostream_iterator<int>(oss, ","));
                safe_group[ijk][s_ij] = oss.str();//stores risk groups in string
//                RiskGroup(ijk,s_ij) = oss.str();//stores risk groups in string
            } else {
                safe_fail[ijk][2*s_ij+0] = -1;
                safe_fail[ijk][2*s_ij+1] = -1;
//                RiskFail(ijk,2*s_ij + 0)=-1;
//                RiskFail(ijk,2*s_ij + 1)=-1;
            }
        }
    }
//    Rcout << "Ending" << endl;
    for (int s_ij=0;s_ij<STRATA_vals.size();s_ij++){
        for (int ijk=0;ijk<ntime;ijk++){
            RiskFail(ijk,2*s_ij + 0)= safe_fail[ijk][2*s_ij+0];
            RiskFail(ijk,2*s_ij + 1)= safe_fail[ijk][2*s_ij+1];
            RiskGroup(ijk,s_ij) = safe_group[ijk][s_ij];
        }
    }
    return;
}


//' Utility function to calculate repeated values used in Cox Log-Likelihood calculation
//' \code{Calculate_Sides} Called to update repeated sum calculations, Uses list of event rows and risk matrices, Performs calculation of sums of risk in each group
//' @param     RiskFail    Matrix of event rows for each event time
//' @param     RiskGroup    vectors of strings with rows at risk for each event time
//' @param     totalnum    total number of parameters
//' @param     ntime    number of event times
//' @param     R    Risk matrix
//' @param     Rd    Risk derivative matrix
//' @param     Rdd    Risk second derivative matrix
//' @param     Rls1    First Risk sum storage
//' @param     Rls2    First Risk sum derivative storage
//' @param     Rls3    First Risk sum second derivative storage
//' @param     Lls1    Second Risk sum storage
//' @param     Lls2    Second Risk sum derivative storage
//' @param     Lls3    Second Risk sum second derivative storage
//' @param     nthreads    number of threads
//' @param     debugging    debugging boolean
//'
//' @return Updates matrices in place: risk storage matrices
// [[Rcpp::export]]
void Calculate_Sides(const IntegerMatrix& RiskFail, const vector<string>&  RiskGroup, const int& totalnum, const int& ntime, const MatrixXd& R, const MatrixXd& Rd, const MatrixXd& Rdd, MatrixXd& Rls1, MatrixXd& Rls2, MatrixXd& Rls3, MatrixXd& Lls1, MatrixXd& Lls2, MatrixXd& Lls3,const int& nthreads, bool debugging){
    //
    //Calculate_Sides( RiskFail, RiskGroup, totalnum, ntime, R, Rd, Rdd, Rls1, Rls2, Rls3, Lls1, Lls2, Lls3,nthreads, debugging);
    //
    if (debugging){
        Rcout << "Starting Debug" << endl;
        vector<double> time_ref(3,0.0);
        vector<double> time_refs(3,0.0);
        vector<int> time_count(3,0);
        Rcout << time_ref.size() << " " << time_refs.size() << endl;
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
                for (vector<double>::size_type i = 0; i < InGroup.size()-1; i=i+2){
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
        Rcout << "df206 ";
        for (vector<double>::size_type ijk=0;ijk<time_ref.size();ijk++){ 
            Rcout << (time_ref[ijk]/time_count[ijk])*1e-6 << " ";
        }
        Rcout << " " << endl;
        Rcout << "df207 ";
        for (vector<double>::size_type ijk=0;ijk<time_ref.size();ijk++){ 
            Rcout <<  sqrt(time_refs[ijk]/time_count[ijk] - pow(time_ref[ijk]/time_count[ijk],2))*1e-6 << " ";
        }
        Rcout << " " << endl;
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
                for (vector<double>::size_type i = 0; i < InGroup.size()-1; i=i+2){
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

//' Utility function to calculate repeated values used in Cox Log-Likelihood calculation with STRATA
//' \code{Calculate_Sides_STRATA} Called to update repeated sum calculations, Uses list of event rows and risk matrices, Performs calculation of sums of risk in each group
//' @param     RiskFail    Matrix of event rows for each event time
//' @param     RiskGroup    vectors of strings with rows at risk for each event time
//' @param     totalnum    total number of parameters
//' @param     ntime    number of event times
//' @param     R    Risk matrix
//' @param     Rd    Risk derivative matrix
//' @param     Rdd    Risk second derivative matrix
//' @param     Rls1    First Risk sum storage
//' @param     Rls2    First Risk sum derivative storage
//' @param     Rls3    First Risk sum second derivative storage
//' @param     Lls1    Second Risk sum storage
//' @param     Lls2    Second Risk sum derivative storage
//' @param     Lls3    Second Risk sum second derivative storage
//' @param     nthreads    number of threads
//' @param     debugging    debugging boolean
//' @param     STRATA_vals vector of strata identifier values
//'
//' @return Updates matrices in place: risk storage matrices
// [[Rcpp::export]]
void Calculate_Sides_STRATA(const IntegerMatrix& RiskFail, const StringMatrix&  RiskGroup, const int& totalnum, const int& ntime, const MatrixXd& R, const MatrixXd& Rd, const MatrixXd& Rdd, MatrixXd& Rls1, MatrixXd& Rls2, MatrixXd& Rls3, MatrixXd& Lls1, MatrixXd& Lls2, MatrixXd& Lls3,const int& nthreads, bool debugging, IntegerVector& STRATA_vals){
    //
    //Calculate_Sides( RiskFail, RiskGroup, totalnum, ntime, R, Rd, Rdd, Rls1, Rls2, Rls3, Lls1, Lls2, Lls3,nthreads, debugging);
    //
    #pragma omp parallel for schedule(dynamic) num_threads(1) collapse(3)
    for (int ijk=0;ijk<totalnum*(totalnum+1)/2;ijk++){//totalnum*(totalnum+1)/2
        for (int j=0;j<ntime;j++){
            for (int s_ij=0;s_ij<STRATA_vals.size();s_ij++){
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
                //Now has the grouping pairs
//                Rcout << ijk << " " << j << " " << s_ij << endl;
                if (RiskFail(j,2*s_ij + 1)>-1){
                    string Groupstr = as<std::string>(RiskGroup(j,s_ij));
                    stringstream ss(Groupstr);
                    //
                    for (int i; ss >> i;) {
                        InGroup.push_back(i);    
                        if (ss.peek() == ',')
                            ss.ignore();
                    }
                    int dj = RiskFail(j,2*s_ij + 1)-RiskFail(j,2*s_ij + 0)+1;
                    for (vector<double>::size_type i = 0; i < InGroup.size()-1; i=i+2){
                        Rs1 += R.block(  InGroup[i]-1, 0,  InGroup[i+1]-InGroup[i]+1,1).sum();
                        Rs2 += Rd.block( InGroup[i]-1, ij, InGroup[i+1]-InGroup[i]+1,1).sum();
                        Rs2t += Rd.block(InGroup[i]-1, jk, InGroup[i+1]-InGroup[i]+1,1).sum();
                        Rs3 += Rdd.block(InGroup[i]-1, ijk,InGroup[i+1]-InGroup[i]+1,1).sum();
                    } //precalculates the sums of risk groups
                    MatrixXd Ld = MatrixXd::Zero(dj,4);
                    Ld << R.block(RiskFail(j,2*s_ij),0,dj,1), Rd.block(RiskFail(j,2*s_ij),ij,dj,1), Rd.block(RiskFail(j,2*s_ij),jk,dj,1) ,Rdd.block(RiskFail(j,2*s_ij),ijk,dj,1);//sum of risks in group
                    // only assigns values once
                    if (ij==jk){
                        if (ij==0){
                            Rls1(j,s_ij) = Rs1;
                            Lls1(j,s_ij) = Ld.col(0).sum();
                        }
                        Rls2(j,ij*STRATA_vals.size() + s_ij) = Rs2;
                        Lls2(j,ij*STRATA_vals.size() + s_ij) = Ld.col(1).sum();
                    }
                    Rls3(j,ijk*STRATA_vals.size() + s_ij) = Rs3;
                    Lls3(j,ijk*STRATA_vals.size() + s_ij) = Ld.col(3).sum();
                }
//                Rcout << ijk << " " << j << " " << s_ij << endl;
            }
        }
    }
    return;
}

//' Utility function to calculate Cox Log-Likelihood and derivatives
//' \code{Calc_LogLik} Called to update log-likelihoods, Uses list of event rows, risk matrices, and repeated sums, Sums the log-likelihood contribution from each event time
//' @param     nthreads    number of threads
//' @param     RiskFail    Matrix of event rows for each event time
//' @param     RiskGroup    vectors of strings with rows at risk for each event time
//' @param     totalnum    total number of parameters
//' @param     ntime    number of event times
//' @param     R    Risk matrix
//' @param     Rd    Risk derivative matrix
//' @param     Rdd    Risk second derivative matrix
//' @param     RdR    Risk to first derivative ratio matrix
//' @param     RddR    Risk to second derivative ratio matrix
//' @param     Rls1    First Risk sum storage
//' @param     Rls2    First Risk sum derivative storage
//' @param     Rls3    First Risk sum second derivative storage
//' @param     Lls1    Second Risk sum storage
//' @param     Lls2    Second Risk sum derivative storage
//' @param     Lls3    Second Risk sum second derivative storage
//' @param     Ll    Log-likelihood vector
//' @param     Lld    Log-likelihood first derivative vector
//' @param     Lldd    Log-likelihood second derivative matrix
//' @param     debugging    debugging boolean
//' @param     ties_method    Ties method
//'
//' @return Updates matrices in place: Log-likelihood vectors/matrix
// [[Rcpp::export]]
void Calc_LogLik(const int& nthreads,const IntegerMatrix& RiskFail, const vector<string>&  RiskGroup, const int& totalnum, const int& ntime, const MatrixXd& R, const MatrixXd& Rd, const MatrixXd& Rdd, const MatrixXd& RdR, const MatrixXd& RddR,const MatrixXd& Rls1,const MatrixXd& Rls2,const MatrixXd& Rls3,const MatrixXd& Lls1,const MatrixXd& Lls2,const MatrixXd& Lls3, vector<double>& Ll, vector<double>& Lld, vector<double>& Lldd, bool debugging,string ties_method){
    //
    //Calc_LogLik( nthreads, RiskFail, RiskGroup, totalnum, ntime, R, Rd, Rdd,RdR,RddR, Rls1, Rls2, Rls3, Lls1, Lls2, Lls3, Ll, Lld, Lldd, debugging);
    //
    if (debugging){
        //
        Rcout << "Starting Debug" << endl;
        vector<double> time_ref(4,0.0);
        vector<double> time_refs(4,0.0);
        vector<int> time_count(4,0);
        Rcout << time_ref.size() << " " << time_refs.size() << endl;
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
        Rcout << "df208 ";
        for (vector<double>::size_type ijk=0;ijk<time_ref.size();ijk++){ 
            Rcout << (time_ref[ijk]/time_count[ijk])*1e-6 << " ";
        }
        Rcout << " " << endl;
        Rcout << "df209 ";
        for (vector<double>::size_type ijk=0;ijk<time_ref.size();ijk++){ 
            Rcout <<  sqrt(time_refs[ijk]/time_count[ijk] - pow(time_ref[ijk]/time_count[ijk],2))*1e-6 << " ";
        }
        Rcout << " " << endl;
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

//' Utility function to calculate Cox Log-Likelihood and derivatives with STRATA
//' \code{Calc_LogLik_STRATA} Called to update log-likelihoods, Uses list of event rows, risk matrices, and repeated sums, Sums the log-likelihood contribution from each event time
//' @param     nthreads    number of threads
//' @param     RiskFail    Matrix of event rows for each event time
//' @param     RiskGroup    vectors of strings with rows at risk for each event time
//' @param     totalnum    total number of parameters
//' @param     ntime    number of event times
//' @param     R    Risk matrix
//' @param     Rd    Risk derivative matrix
//' @param     Rdd    Risk second derivative matrix
//' @param     RdR    Risk to first derivative ratio matrix
//' @param     RddR    Risk to second derivative ratio matrix
//' @param     Rls1    First Risk sum storage
//' @param     Rls2    First Risk sum derivative storage
//' @param     Rls3    First Risk sum second derivative storage
//' @param     Lls1    Second Risk sum storage
//' @param     Lls2    Second Risk sum derivative storage
//' @param     Lls3    Second Risk sum second derivative storage
//' @param     Ll    Log-likelihood vector
//' @param     Lld    Log-likelihood first derivative vector
//' @param     Lldd    Log-likelihood second derivative matrix
//' @param     debugging    debugging boolean
//' @param     ties_method    Ties method
//' @param     STRATA_vals vector of strata identifier values
//'
//' @return Updates matrices in place: Log-likelihood vectors/matrix
// [[Rcpp::export]]
void Calc_LogLik_STRATA(const int& nthreads,const IntegerMatrix& RiskFail, const StringMatrix& RiskGroup, const int& totalnum, const int& ntime, const MatrixXd& R, const MatrixXd& Rd, const MatrixXd& Rdd, const MatrixXd& RdR, const MatrixXd& RddR,const MatrixXd& Rls1,const MatrixXd& Rls2,const MatrixXd& Rls3,const MatrixXd& Lls1,const MatrixXd& Lls2,const MatrixXd& Lls3, vector<double>& Ll, vector<double>& Lld, vector<double>& Lldd, bool debugging,string ties_method, IntegerVector& STRATA_vals){
    //
    //Calc_LogLik( nthreads, RiskFail, RiskGroup, totalnum, ntime, R, Rd, Rdd,RdR,RddR, Rls1, Rls2, Rls3, Lls1, Lls2, Lls3, Ll, Lld, Lldd, debugging);
    //
    #pragma omp declare reduction(vec_double_plus : std::vector<double> : \
        std::transform(omp_out.begin(), omp_out.end(), omp_in.begin(), omp_out.begin(), std::plus<double>())) \
        initializer(omp_priv = omp_orig)
    #pragma omp parallel for schedule(dynamic) num_threads(nthreads) reduction(vec_double_plus:Ll,Lld,Lldd) collapse(3)
    for (int ijk=0;ijk<totalnum*(totalnum+1)/2;ijk++){//performs log-likelihood calculations for every derivative combination and risk group
        for (int j=0;j<ntime;j++){
            for (int s_ij=0;s_ij<STRATA_vals.size();s_ij++){
                int ij = 0;
                int jk = ijk;
                while (jk>ij){
                    ij++;
                    jk-=ij;
                }
                double Rs1 = Rls1(j,s_ij);
                double Rs2 =  Rls2(j,ij*STRATA_vals.size() + s_ij);
                double Rs2t = Rls2(j,jk*STRATA_vals.size() + s_ij);
                double Rs3 = Rls3(j,ijk*STRATA_vals.size() + s_ij);
                //
                int dj = RiskFail(j,2*s_ij + 1)-RiskFail(j,2*s_ij + 0)+1;
                if (RiskFail(j,2*s_ij + 1)>-1){
                    MatrixXd Ld = MatrixXd::Zero(dj,4);
                    Ld << R.block(RiskFail(j,2*s_ij),0,dj,1), RdR.block(RiskFail(j,2*s_ij),ij,dj,1), RdR.block(RiskFail(j,2*s_ij),jk,dj,1) ,RddR.block(RiskFail(j,2*s_ij),ijk,dj,1);//rows with events
                    //
                    MatrixXd Ldm = MatrixXd::Zero(dj,4);
                    Vector4d Ldcs;
                    if (ties_method=="efron"){
                        Ldcs << Lls1(j,s_ij), Lls2(j,ij*STRATA_vals.size() + s_ij), Lls2(j,jk*STRATA_vals.size() + s_ij), Lls3(j,ijk*STRATA_vals.size() + s_ij);
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


//' Utility function to calculate poisson log-likelihood and derivatives
//' \code{Poisson_LogLik} Called to update log-likelihoods, Uses list risk matrices and person-years, Sums the log-likelihood contribution from each row
//' @param     nthreads    number of threads
//' @param     totalnum    total number of parameters
//' @param     PyrC    person-year matrix
//' @param     R    Risk matrix
//' @param     Rd    Risk derivative matrix
//' @param     Rdd    Risk second derivative matrix
//' @param     RdR    Risk to first derivative ratio matrix
//' @param     RddR    Risk to second derivative ratio matrix
//' @param     Ll    Log-likelihood vector
//' @param     Lld    Log-likelihood first derivative vector
//' @param     Lldd    Log-likelihood second derivative matrix
//' @param     debugging    debugging boolean
//'
//' @return Updates matrices in place: Log-likelihood vectors/matrix
// [[Rcpp::export]]
void Poisson_LogLik(const int& nthreads, const int& totalnum, const MatrixXd& PyrC, const MatrixXd& R, const MatrixXd& Rd, const MatrixXd& Rdd, const MatrixXd& RdR, const MatrixXd& RddR, vector<double>& Ll, vector<double>& Lld, vector<double>& Lldd, bool debugging){
    //
    // Poisson_LogLik( nthreads, totalnum, PyrC, R, Rd, Rdd, RdR, RddR, Ll, Lld, Lldd, debugging)
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


//' Utility function to calculate the change to make each iteration
//' \code{Calc_Change} Called to update the parameter changes, Uses log-likelihoods and control parameters, Applys newton steps and change limitations    
//' @param     double_step controls the step calculation, 0 for independent changes, 1 for solving b=Ax with complete matrices
//' @param     nthreads    number of threads
//' @param     totalnum    total number of parameter
//' @param     fir    first term number
//' @param     der_iden    subterm number for derivative tests
//' @param     dbeta_cap    learning rate for newton step toward 0 log-likelihood
//' @param     dose_abs_max    Maximum allowed threshold parameter change
//' @param     lr    learning rate fo newton step toward 0 derivative
//' @param     abs_max    Maximum allowed parameter change
//' @param     Ll    Log-Likelihood
//' @param     Lld    Log-Likelihood first derivative
//' @param     Lldd    Log-Likelihood second derivative
//' @param     dbeta    parameter change vector
//' @param     change_all    boolean to change every parameter
//' @param     tform    subterm type
//' @param     dint    value used for threshold derivative calculation
//' @param     KeepConstant    vector of parameters to keep constant
//' @param     debugging    debugging boolean
//'
//' @return Updates matrices in place: parameter change matrix
// [[Rcpp::export]]
void Calc_Change(const int& double_step, const int& nthreads, const int& totalnum,const int& fir, const int& der_iden, const double& dbeta_cap, const double& dose_abs_max, const double& lr, const double& abs_max, const vector<double>& Ll, const vector<double>& Lld, const vector<double>& Lldd, vector<double>& dbeta, const bool change_all,const StringVector&   tform, const double& dint, IntegerVector KeepConstant, bool debugging){
    //
    //Calc_Change( nthreads, totalnum, fir, der_iden, dbeta_cap, dose_abs_max, lr, abs_max, Ll, Lld, Lldd, dbeta, change_all, tform, dint, KeepConstant, debugging);
    //
    if (debugging){
        if (double_step==1){
            //
            NumericVector Lldd_vec = wrap(Lldd);//
            NumericVector Lld_vec  = wrap(Lld);//
            Lldd_vec.attr("dim") = Dimension(totalnum, totalnum);
            //
            const Map<MatrixXd> Lldd_mat(as<Map<MatrixXd> >(Lldd_vec));
            const Map<VectorXd> Lld_mat(as<Map<VectorXd> >(Lld_vec));
            //
            VectorXd Lldd_solve = Lldd_mat.colPivHouseholderQr().solve(-1*Lld_mat);
            Rcout << Lldd_solve.transpose() << endl;
            //
//            int kept_covs = totalnum - sum(KeepConstant);
//            NumericVector Lldd_vec(kept_covs * kept_cov);
//            NumericVector Lld_vec(kept_covs);
//            #pragma omp parallel for schedule(dynamic) num_threads(nthreads)
//            for (int ijk=0;ijk<totalnum*(totalnum+1)/2;ijk++){
//                int ij = 0;
//                int jk = ijk;
//                while (jk>ij){
//                    ij++;
//                    jk-=ij;
//                }
//                if (KeepConstant[ij]==0){
//                    int pij_ind = ij - sum(head(KeepConstant,ij));
//                    if (ij==jk){
//                        Lld_vec[pij_ind]=Lld[ij];
//                    }
//                    if (KeepConstant[jk]==0){
//                        int pjk_ind = jk - sum(head(KeepConstant,jk));
//                        Lldd_vec[pij_ind * kept_covs + pjk_ind]=Lldd_vec[ijk];
//                    }
//                }
//            }
//            for (int ijk=0;ijk<kept_covs*(kept_covs+1)/2;ijk++){
//                Lldd_vec[pij_ind * kept_covs + pjk_ind]=Lldd_vec[pjk_ind * kept_covs + pij_ind];
//            }
//            Lldd_vec.attr("dim") = Dimension(kept_covs, kept_covs);
//            const Map<MatrixXd> Lldd_mat(as<Map<MatrixXd> >(Lldd_vec));
//            const Map<VectorXd> Lld_mat(as<Map<VectorXd> >(Lld_vec));
//            VectorXd Lldd_solve0 = Lldd_mat.colPivHouseholderQr().solve(-1*Lld_mat);
//            VectorXd Lldd_solve = VectorXd::Zero(totalnum);
//            for (int ij=0;ij<totalnum;ij++){
//                int pij_ind = ij - sum(head(KeepConstant,ij));
//                Lldd_solve(ij) = Lldd_solve0(pij_ind)
//            }
            //
            #pragma omp parallel for schedule(dynamic) num_threads(nthreads)
            for (int ijk=0;ijk<totalnum;ijk++){
                if (change_all){
                    if (KeepConstant[ijk]==0){
                        dbeta[ijk] = Lldd_solve(ijk);//-lr * Lld[ijk] / Lldd[ijk*totalnum+ijk];
                        //
                        double dbeta_max;
                        if (Lld[ijk]!=0){
                            dbeta_max = abs(Ll[ijk]/Lld[ijk] * dbeta_cap);//uses newtonian step for zero log-likelihood as a limit
                        }else{
                            dbeta_max = 0;
                        }
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
                        if (Lldd[ijk*totalnum+ijk] != 0 ){
                            dbeta[ijk] = -lr * Lld[ijk] / Lldd[ijk*totalnum+ijk];
                        } else {
                            dbeta[ijk] = 0;
                        }
                        //
                        double dbeta_max;
                        if (Lld[ijk]!=0){
                            dbeta_max = abs(Ll[ijk]/Lld[ijk] * dbeta_cap);//uses newtonian step for zero log-likelihood as a limit
                        }else{
                            dbeta_max = 0;
                        }
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
        }
    } else {
        if (double_step==1){
//            
//            NumericVector Lldd_vec = wrap(Lldd);//
//            NumericVector Lld_vec  = wrap(Lld);//
//            Lldd_vec.attr("dim") = Dimension(totalnum, totalnum);
//            //
//            const Map<MatrixXd> Lldd_mat(as<Map<MatrixXd> >(Lldd_vec));
//            const Map<VectorXd> Lld_mat(as<Map<VectorXd> >(Lld_vec));
//            //
//            VectorXd Lldd_solve = Lldd_mat.colPivHouseholderQr().solve(-1*Lld_mat);
//            Rcout << Lldd_solve.transpose() << endl;
            //
            int kept_covs = totalnum - sum(KeepConstant);
            NumericVector Lldd_vec(kept_covs * kept_covs);
            NumericVector Lld_vec(kept_covs);
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
                if (KeepConstant[jk]==0){
                    pjk_ind = jk - sum(head(KeepConstant,jk));
                }
                if (KeepConstant[ij]==0){
                    pij_ind = ij - sum(head(KeepConstant,ij));
                    if (ij==jk){
                        Lld_vec[pij_ind]=Lld[ij];
                    }
                    if (KeepConstant[jk]==0){
                        pjk_ind = jk - sum(head(KeepConstant,jk));
                        Lldd_vec[pij_ind * kept_covs + pjk_ind]=Lldd[ij*totalnum+jk];
                    }
                }
//                Rcout << "(" << ij <<"," << jk << ") (" << pij_ind << "," << pjk_ind << ") (" << KeepConstant[ij] << "," << KeepConstant[jk] << ") " << Lldd[ij*totalnum+jk] << endl;
            }
            for (int ijk=0;ijk<kept_covs*(kept_covs+1)/2;ijk++){
                int ij = 0;
                int jk = ijk;
                while (jk>ij){
                    ij++;
                    jk-=ij;
                }
                Lldd_vec[ij * kept_covs + jk]=Lldd_vec[jk * kept_covs + ij];
            }
//            Rcout << Lld_vec << endl;
//            Rcout << Lldd_vec << endl;
//            for (int i=0;i<totalnum;i++){
//                Rcout << Lld[i] << " ";
//            }
//            Rcout << " " << endl;
//            for (int i=0;i<totalnum*totalnum;i++){
//                Rcout << Lldd[i] << " ";
//            }
//            Rcout << " " << endl;
            Lldd_vec.attr("dim") = Dimension(kept_covs, kept_covs);
            const Map<MatrixXd> Lldd_mat(as<Map<MatrixXd> >(Lldd_vec));
            const Map<VectorXd> Lld_mat(as<Map<VectorXd> >(Lld_vec));
            VectorXd Lldd_solve0 = Lldd_mat.colPivHouseholderQr().solve(-1*Lld_mat);
            VectorXd Lldd_solve = VectorXd::Zero(totalnum);
//            Rcout << Lldd_solve0.transpose() << endl;
            for (int ij=0;ij<totalnum;ij++){
                if (KeepConstant[ij]==0){
                    int pij_ind = ij - sum(head(KeepConstant,ij));
                    Lldd_solve(ij) = Lldd_solve0(pij_ind);
                }
            }
//            Rcout << Lldd_solve.transpose() << endl;
            //
            #pragma omp parallel for schedule(dynamic) num_threads(nthreads)
            for (int ijk=0;ijk<totalnum;ijk++){
                if (change_all){
                    if (KeepConstant[ijk]==0){
                        dbeta[ijk] = Lldd_solve(ijk);//-lr * Lld[ijk] / Lldd[ijk*totalnum+ijk];
                        //
//                        double dbeta_max;
//                        if (Lld[ijk]!=0){
//                            dbeta_max = abs(Ll[ijk]/Lld[ijk] * dbeta_cap);//uses newtonian step for zero log-likelihood as a limit
//                        }else{
//                            dbeta_max = 0;
//                        }
//                        if (abs(dbeta[ijk])>dbeta_max){
//                            dbeta[ijk] = dbeta_max * sign(dbeta[ijk]);
//                        }
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
                        if (Lldd[ijk*totalnum+ijk] != 0 ){
                            dbeta[ijk] = -lr * Lld[ijk] / Lldd[ijk*totalnum+ijk];
                        } else {
                            dbeta[ijk] = 0;
                        }
                        //
                        double dbeta_max;
                        if (Lld[ijk]!=0){
                            dbeta_max = abs(Ll[ijk]/Lld[ijk] * dbeta_cap);//uses newtonian step for zero log-likelihood as a limit
                        }else{
                            dbeta_max = 0;
                        }
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
        }
    }
    return;
}

//' Utility function to calculate the change to make each iteration, with basic model
//' \code{Calc_Change_Basic} Called to update the parameter changes, Uses log-likelihoods and control parameters, Applys newton steps and change limitations    
//' @param     double_step controls the step calculation, 0 for independent changes, 1 for solving b=Ax with complete matrices
//' @param     nthreads    number of threads
//' @param     totalnum    total number of parameter
//' @param     der_iden    subterm number for derivative tests
//' @param     dbeta_cap    learning rate for newton step toward 0 log-likelihood
//' @param     lr    learning rate fo newton step toward 0 derivative
//' @param     abs_max    Maximum allowed parameter change
//' @param     Ll    Log-Likelihood
//' @param     Lld    Log-Likelihood first derivative
//' @param     Lldd    Log-Likelihood second derivative
//' @param     dbeta    parameter change vector
//' @param     change_all    boolean to change every parameter
//' @param     KeepConstant    vector of parameters to keep constant
//' @param     debugging    debugging boolean
//'
//' @return Updates matrices in place: parameter change matrix
// [[Rcpp::export]]
void Calc_Change_Basic(const int& double_step, const int& nthreads, const int& totalnum, const int& der_iden, const double& dbeta_cap, const double& lr, const double& abs_max, const vector<double>& Ll, const vector<double>& Lld, const vector<double>& Lldd, vector<double>& dbeta, const bool change_all, IntegerVector KeepConstant, bool debugging){
    if (double_step==1){
        //
        int kept_covs = totalnum - sum(KeepConstant);
        NumericVector Lldd_vec(kept_covs * kept_covs);
        NumericVector Lld_vec(kept_covs);
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
            if (KeepConstant[jk]==0){
                pjk_ind = jk - sum(head(KeepConstant,jk));
            }
            if (KeepConstant[ij]==0){
                pij_ind = ij - sum(head(KeepConstant,ij));
                if (ij==jk){
                    Lld_vec[pij_ind]=Lld[ij];
                }
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
            Lldd_vec[ij * kept_covs + jk]=Lldd_vec[jk * kept_covs + ij];
        }
        Lldd_vec.attr("dim") = Dimension(kept_covs, kept_covs);
        const Map<MatrixXd> Lldd_mat(as<Map<MatrixXd> >(Lldd_vec));
        const Map<VectorXd> Lld_mat(as<Map<VectorXd> >(Lld_vec));
        VectorXd Lldd_solve0 = Lldd_mat.colPivHouseholderQr().solve(-1*Lld_mat);
        VectorXd Lldd_solve = VectorXd::Zero(totalnum);
        for (int ij=0;ij<totalnum;ij++){
            if (KeepConstant[ij]==0){
                int pij_ind = ij - sum(head(KeepConstant,ij));
                Lldd_solve(ij) = Lldd_solve0(pij_ind);
            }
        }
        //
        #pragma omp parallel for schedule(dynamic) num_threads(nthreads)
        for (int ijk=0;ijk<totalnum;ijk++){
            if (change_all){
                if (KeepConstant[ijk]==0){
                    dbeta[ijk] = Lldd_solve(ijk);//-lr * Lld[ijk] / Lldd[ijk*totalnum+ijk];
                    //
                    if (abs(dbeta[ijk])>abs_max){
                        dbeta[ijk] = abs_max * sign(dbeta[ijk]);
                    }
                } else {
                    dbeta[ijk]=0;
                }
            }else{
                if (ijk!=der_iden){//Validation requires controlled changes
                    dbeta[ijk] = 0.0;
                }
            }
        }
    } else {
        #pragma omp parallel for schedule(dynamic) num_threads(nthreads)
        for (int ijk=0;ijk<totalnum;ijk++){
            if (change_all){
                if (KeepConstant[ijk]==0){
                    if (Lldd[ijk*totalnum+ijk] != 0 ){
                        dbeta[ijk] = -lr * Lld[ijk] / Lldd[ijk*totalnum+ijk];
                    } else {
                        dbeta[ijk] = 0;
                    }
                    //
                    double dbeta_max;
                    if (Lld[ijk]!=0){
                        dbeta_max = abs(Ll[ijk]/Lld[ijk] * dbeta_cap);//uses newtonian step for zero log-likelihood as a limit
                    }else{
                        dbeta_max = 0;
                    }
                    if (abs(dbeta[ijk])>dbeta_max){
                        dbeta[ijk] = dbeta_max * sign(dbeta[ijk]);
                    }
                    if (abs(dbeta[ijk])>abs_max){
                        dbeta[ijk] = abs_max * sign(dbeta[ijk]);
                    }
                } else {
                    dbeta[ijk]=0;
                }
            }else{
                if (ijk!=der_iden){//Validation requires controlled changes
                    dbeta[ijk] = 0.0;
                }
            }
        }
    }
    return;
}

//' Utility function to perform null model equivalent of Calculate_Sides
//' \code{Calculate_Null_Sides} Called to update repeated sum calculations, Uses list of event rows, Performs calculation of counts in each group
//' @param     RiskFail    Matrix of event rows for each event time
//' @param     RiskGroup    vectors of strings with rows at risk for each event time
//' @param     ntime    number of event times
//' @param     R    Risk matrix
//' @param     Rls1    First Risk sum storage
//' @param     Lls1    Second Risk sum storage
//' @param     nthreads    number of threads
//'
//' @return Updates matrices in place: risk storage matrices
// [[Rcpp::export]]
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


//' Utility function to perform null model equivalent of Calc_LogLik
//' \code{Calc_Null_LogLik} Called to update log-likelihoods, Uses list of event rows and repeated sums, Sums the log-likelihood contribution from each event time
//' @param     nthreads    number of threads
//' @param     RiskFail    Matrix of event rows for each event time
//' @param     RiskGroup    vectors of strings with rows at risk for each event time
//' @param     ntime    number of event times
//' @param     R    Risk matrix
//' @param     Rls1    First Risk sum storage
//' @param     Lls1    Second Risk sum storage
//' @param     Ll    Log-likelihood vector
//' @param     ties_method    Ties method
//'
//' @return Updates matrices in place: Log-likelihood vectors/matrix
// [[Rcpp::export]]
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










