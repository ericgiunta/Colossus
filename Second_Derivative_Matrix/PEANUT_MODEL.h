using namespace std;
using namespace Rcpp;

using Eigen::Map;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using Rcpp::as;

//----------------------------------------------------------------------------//
//  The transition function is directly called from R
//  This function converts the R input into a matrix format that c++ is more efficient in
//  Returns the log-liklihood
List peanut_transition(NumericVector a_lin,NumericVector a_loglin,NumericVector a_plin,NumericVector a_dose, NumericMatrix x_lin, NumericMatrix x_loglin, NumericMatrix x_plin, NumericMatrix x_dose,int fir,string modelform, string doseform, StringVector dose_terms,int ntime, NumericVector include_bool, double lr, int maxiter, int halfmax, double epsilon, double dbeta_cap, double deriv_epsilon, int batch_size);

//  The LogLik function is directly called from transition
//  This function converts the R input into a matrix format that c++ is more efficient in
//  Currently return the final beta vector
List LogLik_PEANUT( VectorXd beta_linT,VectorXd beta_loglinT,VectorXd beta_plinT, VectorXd beta_dose,MatrixXd df_lin,MatrixXd df_loglin,MatrixXd df_plin, MatrixXd df_dose,int fir,string modelform, string doseform, StringVector dose_terms,int ntime, NumericVector include_bool, double lr, int maxiter, int halfmax, double epsilon, double dbeta_cap, double deriv_epsilon, int batch_size);

template <typename T> int sign(T val);
//----------------------------------------------------------------------------//
