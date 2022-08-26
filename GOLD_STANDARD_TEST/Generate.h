using namespace std;
using namespace Rcpp;

using Eigen::Map;
using Eigen::MatrixXd;
using Eigen::SparseMatrix;
using Eigen::VectorXd;
using Rcpp::as;




void RunGenerator(string fname, NumericVector Pb, int nm_xtra_rows, NumericVector Pcs, StringVector RowNames,int ngoal);
double BaseHazard_CDF(int t, int a, int d);
double HazardRatio(VectorXd C, VectorXd Pc);
