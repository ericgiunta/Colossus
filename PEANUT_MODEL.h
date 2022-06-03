using namespace std;
using namespace Rcpp;

using Eigen::Map;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using Rcpp::as;

//----------------------------------------------------------------------------//
//  The transition function is directly called from R
//  This function converts the R input into a matrix format that c++ is more efficient in
//  Returns the optimal parameters
NumericVector peanut_transition(NumericVector a_lin,NumericVector a_loglin,NumericVector a_plin, NumericMatrix x_lin, NumericMatrix x_loglin, NumericMatrix x_plin,int fir,string modelform,int ntime, NumericVector include_bool);

//  The LogLik function is directly called from transition
//  Currently return the final beta vector
VectorXd LogLik_PEANUT( VectorXd beta_linT,VectorXd beta_loglinT,VectorXd beta_plinT,MatrixXd df_lin,MatrixXd df_loglin,MatrixXd df_plin,int fir,string modelform, string IterStyle,int ntime, NumericVector include_bool);

template <typename T> int sign(T val);
//----------------------------------------------------------------------------//
