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


List peanut_transition(NumericVector a_lin,NumericVector a_loglin,NumericVector a_plin, NumericMatrix x_lin, NumericMatrix x_loglin, NumericMatrix x_plin, NumericMatrix x_dose,int fir,string modelform,int ntime, NumericVector include_bool, List Control, List Dose_paras, NumericMatrix df_groups, NumericVector tu, NumericVector KeepConstant);

//List peanut_bounds_transition(double q1, NumericVector a_lin,NumericVector a_loglin,NumericVector a_plin, NumericMatrix x_lin, NumericMatrix x_loglin, NumericMatrix x_plin, NumericMatrix x_dose,int fir,string modelform,int ntime, NumericVector include_bool, List Control, List Dose_paras, NumericMatrix df_groups, NumericVector tu);

List peanut_plot(NumericVector a_lin,NumericVector a_loglin,NumericVector a_plin, NumericMatrix x_lin, NumericMatrix x_loglin, NumericMatrix x_plin, NumericMatrix x_dose,int fir,string modelform,int ntime, NumericVector include_bool, List Control, List Dose_paras, NumericMatrix df_groups, NumericVector tu ,int uniq_v);

NumericMatrix peanut_schoenfeld_transition(double q1, NumericVector a_lin,NumericVector a_loglin,NumericVector a_plin, NumericMatrix x_lin, NumericMatrix x_loglin, NumericMatrix x_plin, NumericMatrix x_dose,int fir,string modelform,int ntime, NumericVector include_bool, List Control, List Dose_paras, NumericMatrix df_groups, NumericVector tu);

List amfit_transition(NumericMatrix dfe, NumericVector a_lin,NumericVector a_loglin,NumericVector a_plin, NumericMatrix x_lin, NumericMatrix x_loglin, NumericMatrix x_plin, NumericMatrix x_dose,int fir,string modelform, NumericVector include_bool, List Control, List Dose_paras, NumericVector KeepConstant);


//  The LogLik function is directly called from transition
//  This function converts the R input into a matrix format that c++ is more efficient in
//  Currently return the final beta vector
List LogLik_PEANUT( VectorXd beta_linT,VectorXd beta_loglinT,VectorXd beta_plinT,MatrixXd df_lin,MatrixXd df_loglin,MatrixXd df_plin, MatrixXd df_dose,int fir,string modelform,int ntime, NumericVector include_bool, double lr, int maxiter, int halfmax, double epsilon, double dbeta_cap, double abs_max,double dose_abs_max, double deriv_epsilon,List beta_loglin_slopes, List beta_loglin_tops , List beta_lin_slopes , List beta_lin_ints , List beta_quads , List beta_step_slopes , List beta_step_ints, NumericMatrix df_groups, NumericVector tu ,bool change_all, bool verbose, bool debugging, NumericVector KeepConstant);

List PEANUT_PLOT_SURV(VectorXd beta_linT,VectorXd beta_loglinT,VectorXd beta_plinT,MatrixXd df_lin,MatrixXd df_loglin,MatrixXd df_plin, MatrixXd df_dose,int fir,string modelform,int ntime, NumericVector include_bool,List beta_loglin_slopes, List beta_loglin_tops , List beta_lin_slopes , List beta_lin_ints , List beta_quads , List beta_step_slopes , List beta_step_ints, NumericMatrix df_groups, NumericVector tu ,double dose_abs_max, bool verbose, bool debugging);

List PEANUT_PLOT_RISK(VectorXd beta_linT,VectorXd beta_loglinT,VectorXd beta_plinT,MatrixXd df_lin,MatrixXd df_loglin,MatrixXd df_plin, MatrixXd df_dose,int fir,string modelform,int ntime, NumericVector include_bool,List beta_loglin_slopes, List beta_loglin_tops , List beta_lin_slopes , List beta_lin_ints , List beta_quads , List beta_step_slopes , List beta_step_ints, NumericMatrix df_groups, NumericVector tu ,int uniq_v,double dose_abs_max, bool verbose, bool debugging);

//List LogLik_Bounds( double q1, VectorXd beta_linT,VectorXd beta_loglinT,VectorXd beta_plinT,MatrixXd df_lin,MatrixXd df_loglin,MatrixXd df_plin, MatrixXd df_dose,int fir,string modelform,int ntime, NumericVector include_bool, double lr, int maxiter, int halfmax, double epsilon, double dbeta_cap, double abs_max,double dose_abs_max, double deriv_epsilon,List beta_loglin_slopes, List beta_loglin_tops , List beta_lin_slopes , List beta_lin_ints , List beta_quads , List beta_step_slopes , List beta_step_ints, NumericMatrix df_groups, NumericVector tu, bool verbose );

NumericMatrix Schoenfeld_PEANUT( VectorXd beta_linT,VectorXd beta_loglinT,VectorXd beta_plinT,MatrixXd df_lin,MatrixXd df_loglin,MatrixXd df_plin, MatrixXd df_dose,int fir,string modelform,int ntime, NumericVector include_bool,List beta_loglin_slopes, List beta_loglin_tops , List beta_lin_slopes , List beta_lin_ints , List beta_quads , List beta_step_slopes , List beta_step_ints, NumericMatrix df_groups, NumericVector tu, double dose_abs_max , bool verbose, bool debugging);

List LogLik_AMFIT( MatrixXd PyrC, VectorXd beta_linT,VectorXd beta_loglinT,VectorXd beta_plinT,MatrixXd df_lin,MatrixXd df_loglin,MatrixXd df_plin, MatrixXd df_dose,int fir,string modelform, NumericVector include_bool, double lr, int maxiter, int halfmax, double epsilon, double dbeta_cap, double abs_max, double dose_abs_max, double deriv_epsilon,List beta_loglin_slopes, List beta_loglin_tops , List beta_lin_slopes , List beta_lin_ints , List beta_quads , List beta_step_slopes , List beta_step_ints,bool change_all, bool verbose, bool debugging, NumericVector KeepConstant);

void Write_Ind_File(NumericMatrix df, NumericVector tu);


template <typename T> int sign(T val);

void removeRow(Eigen::MatrixXd& matrix, unsigned int rowToRemove);
void removeColumn(Eigen::MatrixXd& matrix, unsigned int colToRemove);
//----------------------------------------------------------------------------//
void Stress_Test(NumericVector a_lin,NumericVector a_loglin,NumericVector a_plin, NumericMatrix x_lin, NumericMatrix x_loglin, NumericMatrix x_plin, NumericMatrix x_dose,int fir,string modelform,int ntime, NumericVector include_bool, List Control, List Dose_paras, NumericMatrix df_groups, NumericVector tu, NumericVector KeepConstant, string test_point);

void Stress_Run( VectorXd beta_linT,VectorXd beta_loglinT,VectorXd beta_plinT,MatrixXd df_lin,MatrixXd df_loglin,MatrixXd df_plin, MatrixXd df_dose,int fir,string modelform,int ntime, NumericVector include_bool, double lr, int maxiter, int halfmax, double epsilon, double dbeta_cap, double abs_max,double dose_abs_max, double deriv_epsilon,List beta_loglin_slopes, List beta_loglin_tops , List beta_lin_slopes , List beta_lin_ints , List beta_quads , List beta_step_slopes , List beta_step_ints, NumericMatrix df_groups, NumericVector tu, bool change_all, bool verbose, bool debugging, NumericVector KeepConstant, StringVector debug_checks);
//----------------------------------------------------------------------------//
void Create_Matrices(const int& totalnum, const int& dose_num_tot, const int& dose_term_tot, vector<int> dose_breaks, vector<vector<double>> beta_loglin_slopes_CPP, vector<vector<double>> beta_loglin_tops_CPP, vector<vector<double>> beta_lin_slopes_CPP, vector<vector<double>> beta_lin_ints_CPP, vector<vector<double>> beta_quads_CPP, vector<vector<double>> beta_step_slopes_CPP, vector<vector<double>> beta_step_ints_CPP, VectorXd& beta_lin,VectorXd& beta_loglin,VectorXd& beta_plin,MatrixXd& df_lin,MatrixXd& df_loglin,MatrixXd& df_plin, MatrixXd& df_dose,  MatrixXd& De, MatrixXd& Dde, MatrixXd& Ddde, MatrixXd& T0, MatrixXd& Td0, MatrixXd& Tdd0, VectorXd& Dose,const vector<double>& cumulative_dose_num, VectorXd& beta_0, MatrixXd& df0, double& dint, int& nthreads, vector<string>&  tform, NumericVector include_bool, bool debugging);

void Make_Risks(string modelform, const int& dose_num_tot, const int& totalnum, const int& fir, const MatrixXd& T0, const MatrixXd& Td0, const MatrixXd& Tdd0, MatrixXd& Te, MatrixXd& R, MatrixXd& Rd, MatrixXd& Rdd, const MatrixXd& De, const MatrixXd& Dde, const MatrixXd& Ddde, const VectorXd& Dose, MatrixXd& RdR, MatrixXd& RddR, const int& nthreads, bool debugging);

void Make_Groups(const int& ntime, const MatrixXd& df_m, IntegerMatrix& RiskFail, vector<string>&  RiskGroup,  NumericVector& tu, const int& nthreads, bool debugging );

void Calculate_Sides(const IntegerMatrix& RiskFail, const vector<string>&  RiskGroup, const int& totalnum, const int& ntime, const MatrixXd& R, const MatrixXd& Rd, const MatrixXd& Rdd, MatrixXd& Rls1, MatrixXd& Rls2, MatrixXd& Rls3, MatrixXd& Lls1, MatrixXd& Lls2, MatrixXd& Lls3,const int& nthreads, bool debugging);

void Calc_LogLik(const int& nthreads,const IntegerMatrix& RiskFail, const vector<string>&  RiskGroup, const int& totalnum, const int& ntime, const MatrixXd& R, const MatrixXd& Rd, const MatrixXd& Rdd, const MatrixXd& RdR, const MatrixXd& RddR,const MatrixXd& Rls1,const MatrixXd& Rls2,const MatrixXd& Rls3,const MatrixXd& Lls1,const MatrixXd& Lls2,const MatrixXd& Lls3, vector<double>& Ll, vector<double>& Lld, vector<double>& Lldd, bool debugging);

void Amfit_LogLik(const int& nthreads, const int& totalnum, const MatrixXd& PyrC, const MatrixXd& R, const MatrixXd& Rd, const MatrixXd& Rdd, const MatrixXd& RdR, const MatrixXd& RddR, vector<double>& Ll, vector<double>& Lld, vector<double>& Lldd, bool debugging);

void Calc_Change(const int& nthreads, const int& totalnum,const int& fir, const double& dbeta_cap, const double& dose_abs_max, const double& lr, const double& abs_max, const vector<double>& Ll, const vector<double>& Lld, const vector<double>& Lldd, vector<double>& dbeta, const bool change_all,vector<string>&  tform, const double& dint, NumericVector KeepConstant, bool debugging);

void Update_Risk(const int& totalnum, const int& dose_num_tot,const VectorXd& beta_0,const MatrixXd& df0, MatrixXd& De, MatrixXd& Dde, MatrixXd& Ddde, MatrixXd& T0, MatrixXd& Td0, MatrixXd& Tdd0, VectorXd& Dose, vector<string>&  tform, const int& nthreads, const double& dint, bool debugging);

//----------------------------------------------------------------------------//
List peanut_null(int ntime, List Control, NumericMatrix df_groups, NumericVector tu);
List LogLik_PEANUT_null( int ntime, NumericMatrix df_groups, NumericVector tu, bool verbose);
void Calculate_Null_Sides(const IntegerMatrix& RiskFail, const vector<string>&  RiskGroup, const int& ntime, const MatrixXd& R, MatrixXd& Rls1, MatrixXd& Lls1,const int& nthreads);
void Calc_Null_LogLik(const int& nthreads,const IntegerMatrix& RiskFail, const vector<string>&  RiskGroup, const int& ntime, const MatrixXd& R, const MatrixXd& Rls1,const MatrixXd& Lls1, vector<double>& Ll);


