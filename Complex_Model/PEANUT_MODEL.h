using namespace std;
using namespace Rcpp;

using Eigen::Map;
using Eigen::MatrixXd;
using Eigen::SparseMatrix;
using Eigen::VectorXd;
using Rcpp::as;

//----------------------------------------------------------------------------//
//  The transition function is directly called from R
//  This function converts the R input into a matrix format that c++ is more efficient in
//  Returns the log-liklihood

List Martingale_Residuals(StringVector IDS, NumericMatrix dfr, NumericVector event_times);

List peanut_transition(NumericVector a_lin,NumericVector a_loglin,NumericVector a_plin, NumericMatrix x_lin, NumericMatrix x_loglin, NumericMatrix x_plin, NumericMatrix x_dose,int fir,string modelform,int ntime, NumericVector include_bool, List Control, List Dose_paras, NumericMatrix df_groups, NumericVector tu);

List peanut_bounds_transition(double q1, NumericVector a_lin,NumericVector a_loglin,NumericVector a_plin, NumericMatrix x_lin, NumericMatrix x_loglin, NumericMatrix x_plin, NumericMatrix x_dose,int fir,string modelform,int ntime, NumericVector include_bool, List Control, List Dose_paras, NumericMatrix df_groups, NumericVector tu);

List peanut_plot(NumericVector a_lin,NumericVector a_loglin,NumericVector a_plin, NumericMatrix x_lin, NumericMatrix x_loglin, NumericMatrix x_plin, NumericMatrix x_dose,int fir,string modelform,int ntime, NumericVector include_bool, List Control, List Dose_paras, NumericMatrix df_groups, NumericVector tu ,int uniq_v);

NumericMatrix peanut_schoenfeld_transition(double q1, NumericVector a_lin,NumericVector a_loglin,NumericVector a_plin, NumericMatrix x_lin, NumericMatrix x_loglin, NumericMatrix x_plin, NumericMatrix x_dose,int fir,string modelform,int ntime, NumericVector include_bool, List Control, List Dose_paras, NumericMatrix df_groups, NumericVector tu);


//  The LogLik function is directly called from transition
//  This function converts the R input into a matrix format that c++ is more efficient in
//  Currently return the final beta vector
List LogLik_PEANUT( VectorXd beta_linT,VectorXd beta_loglinT,VectorXd beta_plinT,MatrixXd df_lin,MatrixXd df_loglin,MatrixXd df_plin, MatrixXd df_dose,int fir,string modelform,int ntime, NumericVector include_bool, double lr, int maxiter, int halfmax, double epsilon, double dbeta_cap, double abs_max,double dose_abs_max, double deriv_epsilon,List beta_loglin_slopes, List beta_loglin_tops , List beta_lin_slopes , List beta_lin_ints , List beta_quads , List beta_step_slopes , List beta_step_ints, NumericMatrix df_groups, NumericVector tu ,bool change_all);

List PEANUT_PLOT_SURV(VectorXd beta_linT,VectorXd beta_loglinT,VectorXd beta_plinT,MatrixXd df_lin,MatrixXd df_loglin,MatrixXd df_plin, MatrixXd df_dose,int fir,string modelform,int ntime, NumericVector include_bool,List beta_loglin_slopes, List beta_loglin_tops , List beta_lin_slopes , List beta_lin_ints , List beta_quads , List beta_step_slopes , List beta_step_ints, NumericMatrix df_groups, NumericVector tu ,double dose_abs_max);

List PEANUT_PLOT_RISK(VectorXd beta_linT,VectorXd beta_loglinT,VectorXd beta_plinT,MatrixXd df_lin,MatrixXd df_loglin,MatrixXd df_plin, MatrixXd df_dose,int fir,string modelform,int ntime, NumericVector include_bool,List beta_loglin_slopes, List beta_loglin_tops , List beta_lin_slopes , List beta_lin_ints , List beta_quads , List beta_step_slopes , List beta_step_ints, NumericMatrix df_groups, NumericVector tu ,int uniq_v,double dose_abs_max);

List LogLik_Bounds( double q1, VectorXd beta_linT,VectorXd beta_loglinT,VectorXd beta_plinT,MatrixXd df_lin,MatrixXd df_loglin,MatrixXd df_plin, MatrixXd df_dose,int fir,string modelform,int ntime, NumericVector include_bool, double lr, int maxiter, int halfmax, double epsilon, double dbeta_cap, double abs_max,double dose_abs_max, double deriv_epsilon,List beta_loglin_slopes, List beta_loglin_tops , List beta_lin_slopes , List beta_lin_ints , List beta_quads , List beta_step_slopes , List beta_step_ints, NumericMatrix df_groups, NumericVector tu );

NumericMatrix Schoenfeld_PEANUT( VectorXd beta_linT,VectorXd beta_loglinT,VectorXd beta_plinT,MatrixXd df_lin,MatrixXd df_loglin,MatrixXd df_plin, MatrixXd df_dose,int fir,string modelform,int ntime, NumericVector include_bool,List beta_loglin_slopes, List beta_loglin_tops , List beta_lin_slopes , List beta_lin_ints , List beta_quads , List beta_step_slopes , List beta_step_ints, NumericMatrix df_groups, NumericVector tu, double dose_abs_max );

void Write_Ind_File(NumericMatrix df, NumericVector tu);


template <typename T> int sign(T val);

void removeRow(Eigen::MatrixXd& matrix, unsigned int rowToRemove);
void removeColumn(Eigen::MatrixXd& matrix, unsigned int colToRemove);
//----------------------------------------------------------------------------//
