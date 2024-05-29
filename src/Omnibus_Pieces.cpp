#include <RcppEigen.h>
#ifdef _OPENMP
#include <omp.h>
#endif
#include "Omnibus_Pieces.h"
#include "Calc_Repeated.h"
#include "Subterms_Risk.h"
#include "Colossus_types.h"
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


//' Utility function to refresh risk and subterm matrices for Cox Omnibus function
//'
//' \code{Cox_Refresh_R_TERM} Called to update matrices
//' @inheritParams CPP_template
//'
//' @return Updates matrices in place
//' @noRd
//'
// [[Rcpp::export]]
void Cox_Refresh_R_TERM(const int& totalnum, const int& reqrdnum, const int& term_tot, double& dint, double& dslp,double& dose_abs_max, double& abs_max, const MatrixXd& df0, MatrixXd& T0, MatrixXd& Td0, MatrixXd& Tdd0, MatrixXd& Te, MatrixXd& R, MatrixXd& Rd, MatrixXd& Rdd, MatrixXd& Dose, MatrixXd& nonDose,  MatrixXd& TTerm,  MatrixXd& nonDose_LIN, MatrixXd& nonDose_PLIN, MatrixXd& nonDose_LOGLIN, MatrixXd& RdR, MatrixXd& RddR, bool basic_bool, bool single_bool){
    T0 = MatrixXd::Zero(df0.rows(), totalnum); //preallocates matrix for Term column
	if (basic_bool){
		//
		R = MatrixXd::Zero(df0.rows(), 1); //preallocates matrix for Risks
		Rd = MatrixXd::Zero(df0.rows(), reqrdnum); //preallocates matrix for Risk derivatives
		Rdd = MatrixXd::Zero(df0.rows(), reqrdnum*(reqrdnum+1)/2); //preallocates matrix for Risk second derivatives
        RdR = MatrixXd::Zero(df0.rows(), reqrdnum); //preallocates matrix for Risk to derivative ratios
        TTerm = MatrixXd::Zero(df0.rows(),1); //matrix of term values
    } else if (single_bool){
		//
		Te = MatrixXd::Zero(df0.rows(), 1); //preallocates matrix for column terms used for temporary storage
		R = MatrixXd::Zero(df0.rows(), 1); //preallocates matrix for Risks
		//
		Dose = MatrixXd::Constant(df0.rows(),term_tot,0.0); //Matrix of the total dose term values
		nonDose = MatrixXd::Constant(df0.rows(),term_tot,1.0); //Matrix of the total non-dose term values
		nonDose_LIN = MatrixXd::Constant(df0.rows(),term_tot,0.0); //matrix of Linear subterm values
		nonDose_PLIN = MatrixXd::Constant(df0.rows(),term_tot,1.0); //matrix of Loglinear subterm values
		nonDose_LOGLIN = MatrixXd::Constant(df0.rows(),term_tot,1.0); //matrix of Product linear subterm values
		TTerm = MatrixXd::Zero(Dose.rows(),Dose.cols()); //matrix of term values
    } else {
		Td0 = MatrixXd::Zero(df0.rows(), reqrdnum); //preallocates matrix for Term derivative columns
		Tdd0 = MatrixXd::Zero(df0.rows(), reqrdnum*(reqrdnum+1)/2); //preallocates matrix for Term second derivative columns
		//
		Te = MatrixXd::Zero(df0.rows(), 1); //preallocates matrix for column terms used for temporary storage
		R = MatrixXd::Zero(df0.rows(), 1); //preallocates matrix for Risks
		Rd = MatrixXd::Zero(df0.rows(), reqrdnum); //preallocates matrix for Risk derivatives
		Rdd = MatrixXd::Zero(df0.rows(), reqrdnum*(reqrdnum+1)/2); //preallocates matrix for Risk second derivatives
		//
		Dose = MatrixXd::Constant(df0.rows(),term_tot,0.0); //Matrix of the total dose term values
		nonDose = MatrixXd::Constant(df0.rows(),term_tot,1.0); //Matrix of the total non-dose term values
		nonDose_LIN = MatrixXd::Constant(df0.rows(),term_tot,0.0); //matrix of Linear subterm values
		nonDose_PLIN = MatrixXd::Constant(df0.rows(),term_tot,1.0); //matrix of Loglinear subterm values
		nonDose_LOGLIN = MatrixXd::Constant(df0.rows(),term_tot,1.0); //matrix of Product linear subterm values
		TTerm = MatrixXd::Zero(Dose.rows(),Dose.cols()); //matrix of term values
		dint = dose_abs_max; //The amount of change used to calculate derivatives in threshold paramters
		dslp = abs_max;
        RdR = MatrixXd::Zero(df0.rows(), reqrdnum); //preallocates matrix for Risk to derivative ratios
		RddR = MatrixXd::Zero(df0.rows(), reqrdnum*(reqrdnum+1)/2); //preallocates matrix for Risk to second derivative ratios
    }
    return;
}

//' Utility function to refresh side matrices for Cox Omnibus
//'
//' \code{Cox_Refresh_R_SIDES} Called to fresh repeated sum calculation matrices
//' @inheritParams CPP_template
//'
//' @return Updates matrices in place: risk storage matrices
//' @noRd
//'
// [[Rcpp::export]]
void Cox_Refresh_R_SIDES( const int& reqrdnum, const int& ntime, MatrixXd& Rls1, MatrixXd& Rls2, MatrixXd& Rls3, MatrixXd& Lls1, MatrixXd& Lls2, MatrixXd& Lls3, NumericVector& STRATA_vals, bool strata_bool, bool single_bool){
    if (strata_bool){
		Rls1 =MatrixXd::Zero(ntime, STRATA_vals.size()); //precomputes a series of sums used frequently in the log-liklihood calculations
		Lls1 =MatrixXd::Zero(ntime, STRATA_vals.size()); //The log-likelihood calculation has a Right and Left sum used
		if (!single_bool){
			Rls2 =MatrixXd::Zero(ntime, reqrdnum*STRATA_vals.size()); //Many are repeated due to the same risk groups and derivatives being used at mulitple points
			Rls3 =MatrixXd::Zero(ntime, reqrdnum*(reqrdnum+1)/2*STRATA_vals.size()); //Sum and its derivatives are precomputed
			Lls2 =MatrixXd::Zero(ntime, reqrdnum*STRATA_vals.size());
			Lls3 =MatrixXd::Zero(ntime, reqrdnum*(reqrdnum+1)/2*STRATA_vals.size());
		}
	} else {
		Rls1 =MatrixXd::Zero(ntime, 1); //precomputes a series of sums used frequently in the log-liklihood calculations
		Lls1 =MatrixXd::Zero(ntime, 1); //The log-likelihood calculation has a Right and Left sum used
		if (!single_bool){
			Rls2 =MatrixXd::Zero(ntime, reqrdnum); //Many are repeated due to the same risk groups and derivatives being used at mulitple points
			Rls3 =MatrixXd::Zero(ntime, reqrdnum*(reqrdnum+1)/2); //Sum and its derivatives are precomputed
			Lls2 =MatrixXd::Zero(ntime, reqrdnum);
			Lls3 =MatrixXd::Zero(ntime, reqrdnum*(reqrdnum+1)/2);
		}
	}
	return;
}



//' Utility function to perform calculation of terms and risks for Cox Omnibus
//'
//' \code{Cox_Term_Risk_Calc} Called to perform repeated term and risk calculations
//' @inheritParams CPP_template
//'
//' @return Updates matrices in place: risk storage matrices
//' @noRd
//'
// [[Rcpp::export]]
void Cox_Term_Risk_Calc(string modelform, const StringVector& tform, const IntegerVector& Term_n, const int& totalnum, const int& fir, const IntegerVector& dfc, int term_tot, MatrixXd& T0, MatrixXd& Td0, MatrixXd& Tdd0, MatrixXd& Te, MatrixXd& R, MatrixXd& Rd, MatrixXd& Rdd, MatrixXd& Dose, MatrixXd& nonDose, VectorXd beta_0,const  MatrixXd& df0,const double& dint, const double& dslp,  MatrixXd& TTerm,  MatrixXd& nonDose_LIN, MatrixXd& nonDose_PLIN, MatrixXd& nonDose_LOGLIN, MatrixXd& RdR, MatrixXd& RddR, const int& nthreads, bool debugging, const IntegerVector& KeepConstant, bool verbose, bool basic_bool, bool single_bool, int start, const double gmix_theta, const IntegerVector& gmix_term){
    int reqrdnum = totalnum - sum(KeepConstant);
    if (basic_bool){
		// Calculates the subterm and term values
		Make_subterms_Basic( totalnum,  dfc,  T0 ,beta_0, df0,nthreads, debugging);
		// ---------------------------------------------------------
		// Prints off a series of calculations to check at what point values are changing
		// ---------------------------------------------------------
		//
		if (verbose){
			Rcout << "C++ Note: values checked ";
			for (int ijk=0;ijk<totalnum;ijk++){
				Rcout << beta_0[ijk] << " ";
			}
			Rcout << " " << endl;
			Rcout << "C++ Note: sums checked ";
			for (int ijk=0;ijk<totalnum;ijk++){
				Rcout << T0.col(ijk).sum() << " ";
			}
			Rcout << " " << endl;
		}
		//
		// Calculates the risk for each row
		Make_Risks_Basic(totalnum, T0, R, Rd, Rdd, RdR, nthreads, debugging,df0,dfc,KeepConstant);
		//
		// Removes infinite values
		RdR = (RdR.array().isFinite()).select(RdR,0);
		//
		//
		if (R.minCoeff()<=0){
		    ;
		} else if (verbose){
			Rcout << "C++ Note: risk checked ";
			for (int ijk=0;ijk<1;ijk++){
				Rcout << R.col(0).sum() << " ";
			}
			Rcout << " " << endl;
			Rcout << "C++ Note: risk1 checked ";
			for (int ijk=0;ijk<reqrdnum;ijk++){
				Rcout << Rd.col(ijk).sum() << " ";
			}
			Rcout << " " << endl;
			Rcout << "C++ Note: risk2 checked ";
			for (int ijk=0;ijk<reqrdnum;ijk++){
				Rcout << Rdd.col(ijk*(ijk+1)/2+ijk).sum() << " ";
			}
			Rcout << " " << endl;
			//
		}
	} else if (single_bool){
		// Calculates the subterm and term values
		Make_subterms_Single( totalnum, Term_n, tform, dfc, fir, T0, Dose, nonDose, TTerm,  nonDose_LIN, nonDose_PLIN, nonDose_LOGLIN ,beta_0, df0,nthreads, debugging,KeepConstant);
		// ---------------------------------------------------------
		// Prints off a series of calculations to check at what point values are changing
		// ---------------------------------------------------------
		//
		if (verbose){
			Rcout << "C++ Note: values checked ";
			for (int ijk=0;ijk<totalnum;ijk++){
				Rcout << beta_0[ijk] << " ";
			}
			Rcout << " " << endl;
			Rcout << "C++ Note: sums checked ";
			for (int ijk=0;ijk<totalnum;ijk++){
				Rcout << T0.col(ijk).sum() << " ";
			}
			Rcout << " " << endl;
			Rcout << "C++ Note: dose checked ";
			for (int ijk=0;ijk<term_tot;ijk++){
				Rcout << Dose.col(ijk).array().sum() << " ";
			}
			Rcout << " " << endl;
			Rcout << "C++ Note: non-dose checked ";
			for (int ijk=0;ijk<term_tot;ijk++){
				Rcout << nonDose.col(ijk).array().sum() << " ";
			}
			Rcout << " " << endl;
			Rcout << "C++ Note: LIN_non-dose checked ";
			for (int ijk=0;ijk<term_tot;ijk++){
				Rcout << nonDose_LIN.col(ijk).array().sum() << " ";
			}
			Rcout << " " << endl;
			Rcout << "C++ Note: PLIN_non-dose checked ";
			for (int ijk=0;ijk<term_tot;ijk++){
				Rcout << nonDose_PLIN.col(ijk).array().sum() << " ";
			}
			Rcout << " " << endl;
			Rcout << "C++ Note: LOGLIN_non-dose checked ";
			for (int ijk=0;ijk<term_tot;ijk++){
				Rcout << nonDose_LOGLIN.col(ijk).array().sum() << " ";
			}
			Rcout << " " << endl;
		}
		//
		//
		// Calculates the risk for each row
        Make_Risks_Single(modelform, tform, Term_n, totalnum, fir, T0, Te, R, Dose, nonDose, TTerm, nonDose_LIN, nonDose_PLIN, nonDose_LOGLIN, nthreads, debugging,KeepConstant,gmix_theta, gmix_term);
		//
		// Removes infinite values
		//
		//
		if (R.minCoeff()<=0){
		    ;
		} else if (verbose){
			Rcout << "C++ Note: risk checked ";
			for (int ijk=0;ijk<1;ijk++){
				Rcout << R.col(0).sum() << " ";
			}
			Rcout << " " << endl;
		}
	} else {
		//
		// Calculates the subterm and term values
		//
		Make_subterms( totalnum, Term_n, tform, dfc, fir, T0, Td0, Tdd0, Dose, nonDose, TTerm,  nonDose_LIN, nonDose_PLIN, nonDose_LOGLIN ,beta_0, df0,dint,dslp,nthreads, debugging,KeepConstant);
		// ---------------------------------------------------------
		// Prints off a series of calculations to check at what point values are changing
		// ---------------------------------------------------------
		//
		//
		if (verbose){
			Rcout << "C++ Note: values checked ";
			for (int ijk=0;ijk<totalnum;ijk++){
				Rcout << beta_0[ijk] << " ";
			}
			Rcout << " " << endl;
			Rcout << "C++ Note: sums checked ";
			for (int ijk=0;ijk<totalnum;ijk++){
				Rcout << T0.col(ijk).sum() << " ";
			}
			Rcout << " " << endl;
			Rcout << "C++ Note: derivs checked ";
			for (int ijk=0;ijk<reqrdnum;ijk++){
				Rcout << Td0.col(ijk).sum() << " ";
			}
			Rcout << " " << endl;
			Rcout << "C++ Note: second derivs checked ";
			for (int ijk=0;ijk<reqrdnum;ijk++){
				Rcout << Tdd0.col(ijk*(ijk+1)/2+ijk).sum() << " ";
			}
			Rcout << " " << endl;
			Rcout << "C++ Note: dose checked ";
			for (int ijk=0;ijk<term_tot;ijk++){
				Rcout << Dose.col(ijk).array().sum() << " ";
			}
			Rcout << " " << endl;
			Rcout << "C++ Note: non-dose checked ";
			for (int ijk=0;ijk<term_tot;ijk++){
				Rcout << nonDose.col(ijk).array().sum() << " ";
			}
			Rcout << " " << endl;
			Rcout << "C++ Note: LIN_non-dose checked ";
			for (int ijk=0;ijk<term_tot;ijk++){
				Rcout << nonDose_LIN.col(ijk).array().sum() << " ";
			}
			Rcout << " " << endl;
			Rcout << "C++ Note: PLIN_non-dose checked ";
			for (int ijk=0;ijk<term_tot;ijk++){
				Rcout << nonDose_PLIN.col(ijk).array().sum() << " ";
			}
			Rcout << " " << endl;
			Rcout << "C++ Note: LOGLIN_non-dose checked ";
			for (int ijk=0;ijk<term_tot;ijk++){
				Rcout << nonDose_LOGLIN.col(ijk).array().sum() << " ";
			}
			Rcout << " " << endl;
		}
		//
		//
		// Calculates the risk for each row
        Make_Risks(modelform, tform, Term_n, totalnum, fir, T0, Td0, Tdd0, Te, R, Rd, Rdd, Dose, nonDose, TTerm, nonDose_LIN, nonDose_PLIN, nonDose_LOGLIN, RdR, RddR, nthreads, debugging,KeepConstant,gmix_theta, gmix_term);
		//
		// Removes infinite values
		RdR = (RdR.array().isFinite()).select(RdR,0);
		RddR = (RddR.array().isFinite()).select(RddR,0);
		//
		//
		if (R.minCoeff()<=0){
		    ;
		} else if (verbose){
			Rcout << "C++ Note: risk checked ";
			for (int ijk=0;ijk<1;ijk++){
				Rcout << R.col(0).sum() << " ";
			}
			Rcout << " " << endl;
			Rcout << "C++ Note: risk1 checked ";
			for (int ijk=0;ijk<reqrdnum;ijk++){
				Rcout << Rd.col(ijk).sum() << " ";
			}
			Rcout << " " << endl;
			Rcout << "C++ Note: risk2 checked ";
			for (int ijk=0;ijk<reqrdnum;ijk++){
				Rcout << Rdd.col(ijk*(ijk+1)/2+ijk).sum() << " ";
			}
			Rcout << " " << endl;
			//
		}
	}
	return;
}

//' Utility function to perform calculation of Repeated Calculations and Log-Likelihood for Cox Omnibus
//'
//' \code{Cox_Side_LL_Calc} Called to perform repeated term and risk calculations
//' @inheritParams CPP_template
//'
//' @return Updates matrices in place: risk storage matrices
//' @noRd
//'
// [[Rcpp::export]]
void Cox_Side_LL_Calc(const int& reqrdnum, const int& ntime, const IntegerMatrix& RiskFail, const StringMatrix&  RiskGroup_Strata, const vector<string>& RiskGroup, const int& totalnum, const int& fir, MatrixXd& R, MatrixXd& Rd, MatrixXd& Rdd,  MatrixXd& Rls1, MatrixXd& Rls2, MatrixXd& Rls3, MatrixXd& Lls1, MatrixXd& Lls2, MatrixXd& Lls3, const VectorXd& cens_weight, NumericVector& STRATA_vals, VectorXd beta_0 , MatrixXd& RdR, MatrixXd& RddR, vector<double>& Ll, vector<double>& Lld, vector<double>& Lldd, const int& nthreads, bool debugging, const IntegerVector& KeepConstant,string ties_method, bool verbose,bool strata_bool, bool CR_bool, bool basic_bool, bool single_bool, int start, int iter_stop){
    // Calculates the side sum terms used
    if (strata_bool){
        if (CR_bool){
            ;//strata_CR or strata_CR_single
            if (single_bool){
                Calculate_Sides_STRATA_Single_CR( RiskFail, RiskGroup_Strata, totalnum, ntime, R, Rls1, Lls1,cens_weight, nthreads, debugging, STRATA_vals,KeepConstant);//strata_CR_single
            } else {
                Calculate_Sides_STRATA_CR( RiskFail, RiskGroup_Strata, totalnum, ntime, R, Rd, Rdd, Rls1, Rls2, Rls3, Lls1, Lls2, Lls3,cens_weight,nthreads, debugging, STRATA_vals,KeepConstant);//strata_cr
            }
        } else if (single_bool){
            Calculate_Sides_STRATA_Single( RiskFail, RiskGroup_Strata, totalnum, ntime, R, Rls1, Lls1, nthreads, debugging, STRATA_vals,KeepConstant);//strata_single
        } else {
            Calculate_Sides_STRATA( RiskFail, RiskGroup_Strata, totalnum, ntime, R, Rd, Rdd, Rls1, Rls2, Rls3, Lls1, Lls2, Lls3,nthreads, debugging, STRATA_vals,KeepConstant);
        }
    } else if (CR_bool) {
        if (single_bool){
            Calculate_Sides_CR_SINGLE( RiskFail, RiskGroup, totalnum, ntime, R, Rls1, Lls1,cens_weight,nthreads, debugging,KeepConstant);
        } else {
            Calculate_Sides_CR( RiskFail, RiskGroup, totalnum, ntime, R, Rd, Rdd, Rls1, Rls2, Rls3, Lls1, Lls2, Lls3,cens_weight,nthreads, debugging,KeepConstant);
        }
    } else if (single_bool) {
        Calculate_Sides_Single( RiskFail, RiskGroup, totalnum, ntime, R, Rls1, Lls1,nthreads, debugging);
    } else {
        Calculate_Sides( RiskFail, RiskGroup, totalnum, ntime, R, Rd, Rdd, Rls1, Rls2, Rls3, Lls1, Lls2, Lls3,nthreads, debugging,KeepConstant);
    }
    //
    if (strata_bool){
        if (verbose){
            Rcout << "C++ Note: riskr checked ";
            for (int ijk=0;ijk<STRATA_vals.size();ijk++){
                Rcout << Rls1.col(ijk).sum() << " ";
            }
            Rcout << " " << endl;
            //
            Rcout << "C++ Note: riskl checked ";
            for (int ijk=0;ijk<STRATA_vals.size();ijk++){
                Rcout << Lls1.col(ijk).sum() << " ";
            }
            Rcout << " " << endl;
        }
    } else {
        if (verbose){
            Rcout << "C++ Note: riskr checked ";
            for (int ijk=0;ijk<1;ijk++){
                Rcout << Rls1.col(0).sum() << " ";
            }
            Rcout << " " << endl;
            if (!single_bool){
                Rcout << "C++ Note: risk1r checked ";
                for (int ijk=0;ijk<reqrdnum;ijk++){
                    Rcout << Rls2.col(ijk).sum() << " ";
                }
                Rcout << " " << endl;
                Rcout << "C++ Note: risk2r checked ";
                for (int ijk=0;ijk<reqrdnum;ijk++){
                    Rcout << Rls3.col(ijk*(ijk+1)/2+ijk).sum() << " ";
                }
                Rcout << " " << endl;
            }
            //
            Rcout << "C++ Note: riskl checked ";
            for (int ijk=0;ijk<1;ijk++){
                Rcout << Lls1.col(0).sum() << " ";
            }
            Rcout << " " << endl;
            if (!single_bool){
                Rcout << "C++ Note: risk1l checked ";
                for (int ijk=0;ijk<reqrdnum;ijk++){
                    Rcout << Lls2.col(ijk).sum() << " ";
                }
                Rcout << " " << endl;
                Rcout << "C++ Note: risk2l checked ";
                for (int ijk=0;ijk<reqrdnum;ijk++){
                    Rcout << Lls3.col(ijk*(ijk+1)/2+ijk).sum() << " ";
                }
                Rcout << " " << endl;
            }
        }
    }
    // Calculates log-likelihood
    fill(Ll.begin(), Ll.end(), 0.0);
    if (!single_bool){
        fill(Lld.begin(), Lld.end(), 0.0);
        fill(Lldd.begin(), Lldd.end(), 0.0);
    }
    if (strata_bool){
        if (basic_bool){
            if (single_bool){
                Calc_LogLik_STRATA_BASIC_SINGLE( nthreads, RiskFail, RiskGroup_Strata, totalnum, ntime, R, Rls1, Lls1, Ll, debugging, ties_method, STRATA_vals,KeepConstant);
            } else {
                Calc_LogLik_STRATA_BASIC( nthreads, RiskFail, RiskGroup_Strata, totalnum, ntime, R, Rd, Rdd,RdR, Rls1, Rls2, Rls3, Lls1, Lls2, Lls3, Ll, Lld, Lldd, debugging, ties_method, STRATA_vals,KeepConstant);
            }
        } else if (single_bool){
            Calc_LogLik_STRATA_SINGLE( nthreads, RiskFail, RiskGroup_Strata, totalnum, ntime, R, Rls1, Lls1, Ll, debugging, ties_method, STRATA_vals,KeepConstant);//strata_single
        } else {
            Calc_LogLik_STRATA( nthreads, RiskFail, RiskGroup_Strata, totalnum, ntime, R, Rd, Rdd,RdR,RddR, Rls1, Rls2, Rls3, Lls1, Lls2, Lls3, Ll, Lld, Lldd, debugging, ties_method, STRATA_vals,KeepConstant);
        }
    } else {
        if (basic_bool){
            if (single_bool){
                Calc_LogLik_Basic_Single( nthreads, RiskFail,  RiskGroup, totalnum, ntime, R, Rls1, Lls1, Ll, debugging, ties_method, KeepConstant);//basic_single
            } else {
                Calc_LogLik_Basic( nthreads, RiskFail, RiskGroup, totalnum, ntime, R, Rd, Rdd,RdR, Rls1, Rls2, Rls3, Lls1, Lls2, Lls3, Ll, Lld, Lldd, debugging, ties_method,KeepConstant);
            }
        } else {
            if (single_bool){
                Calc_LogLik_Single( nthreads, RiskFail, RiskGroup, totalnum, ntime, R, Rls1, Lls1, Ll, debugging, ties_method);//single
            } else {
                Calc_LogLik( nthreads, RiskFail, RiskGroup, totalnum, ntime, R, Rd, Rdd,RdR,RddR, Rls1, Rls2, Rls3, Lls1, Lls2, Lls3, Ll, Lld, Lldd, debugging, ties_method,KeepConstant);
            }
        }
    }
    //
    if (single_bool){
        iter_stop=1;
        if (verbose){
            Rcout << "C++ Note: df101 ";//prints the log-likelihoods
            for (int ij=0;ij<reqrdnum;ij++){
                Rcout << Ll[ij] << " ";
            }
            Rcout << " " << endl;
            Rcout << "C++ Note: df104 ";//prints parameter values
            for (int ij=0;ij<totalnum;ij++){
                Rcout << beta_0[ij] << " ";
            }
            Rcout << " " << endl;
        }
    } else {
        if (verbose){
            Rcout << "C++ Note: df101 ";//prints the log-likelihoods
            for (int ij=0;ij<reqrdnum;ij++){
                Rcout << Ll[ij] << " ";
            }
            Rcout << " " << endl;
            Rcout << "C++ Note: df102 ";//prints the first derivatives
            for (int ij=0;ij<reqrdnum;ij++){
                Rcout << Lld[ij] << " ";
            }
            Rcout << " " << endl;
            Rcout << "C++ Note: df103 ";//prints the second derivatives
            for (int ij=0;ij<reqrdnum;ij++){
                Rcout << Lldd[ij*reqrdnum+ij] << " ";
            }
            Rcout << " " << endl;
            Rcout << "C++ Note: df104 ";//prints parameter values
            for (int ij=0;ij<totalnum;ij++){
                Rcout << beta_0[ij] << " ";
            }
            Rcout << " " << endl;
            Rcout << "C++ Note: df105 ";
            for (int ij=0;ij<reqrdnum;ij++){//prints the newton step value for zero derivative
                Rcout << Lld[ij]/Lldd[ij*reqrdnum+ij] << " ";
            }
            Rcout << " " << endl;
            Rcout << "C++ Note: df106 ";
            for (int ij=0;ij<reqrdnum;ij++){//prints the newton step value for zero log-likelihood
                Rcout << Ll[ij]/Lld[ij] << " ";
            }
            Rcout << " " << endl;
        }
    }
}

//' Utility function to perform calculation of terms and risks for Poisson Omnibus
//'
//' \code{Pois_Term_Risk_Calc} Called to perform repeated term and risk calculations
//' @inheritParams CPP_template
//'
//' @return Updates matrices in place: risk storage matrices
//' @noRd
//'
// [[Rcpp::export]]
void Pois_Term_Risk_Calc(string modelform, const StringVector& tform, const IntegerVector& Term_n, const int& totalnum, const int& fir, const IntegerVector& dfc, int term_tot, MatrixXd& T0, MatrixXd& Td0, MatrixXd& Tdd0, MatrixXd& Te, MatrixXd& R, MatrixXd& Rd, MatrixXd& Rdd, MatrixXd& Dose, MatrixXd& nonDose, VectorXd beta_0,const  MatrixXd& df0,const double& dint, const double& dslp,  MatrixXd& TTerm,  MatrixXd& nonDose_LIN, MatrixXd& nonDose_PLIN, MatrixXd& nonDose_LOGLIN, MatrixXd& RdR, MatrixXd& RddR, const MatrixXd& s_weights, const int& nthreads, bool debugging, const IntegerVector& KeepConstant, bool verbose, bool strata_bool, bool single_bool, int start, const double gmix_theta, const IntegerVector& gmix_term){
    int reqrdnum = totalnum - sum(KeepConstant);
    if (single_bool){
		// Calculates the subterm and term values
		Make_subterms_Single( totalnum, Term_n, tform, dfc, fir, T0, Dose, nonDose, TTerm,  nonDose_LIN, nonDose_PLIN, nonDose_LOGLIN ,beta_0, df0,nthreads, debugging,KeepConstant);
		// ---------------------------------------------------------
		// Prints off a series of calculations to check at what point values are changing
		// ---------------------------------------------------------
		//
		if (verbose){
			Rcout << "C++ Note: values checked ";
			for (int ijk=0;ijk<totalnum;ijk++){
				Rcout << beta_0[ijk] << " ";
			}
			Rcout << " " << endl;
			Rcout << "C++ Note: sums checked ";
			for (int ijk=0;ijk<totalnum;ijk++){
				Rcout << T0.col(ijk).sum() << " ";
			}
			Rcout << " " << endl;
			Rcout << "C++ Note: dose checked ";
			for (int ijk=0;ijk<term_tot;ijk++){
				Rcout << Dose.col(ijk).array().sum() << " ";
			}
			Rcout << " " << endl;
			Rcout << "C++ Note: non-dose checked ";
			for (int ijk=0;ijk<term_tot;ijk++){
				Rcout << nonDose.col(ijk).array().sum() << " ";
			}
			Rcout << " " << endl;
			Rcout << "C++ Note: LIN_non-dose checked ";
			for (int ijk=0;ijk<term_tot;ijk++){
				Rcout << nonDose_LIN.col(ijk).array().sum() << " ";
			}
			Rcout << " " << endl;
			Rcout << "C++ Note: PLIN_non-dose checked ";
			for (int ijk=0;ijk<term_tot;ijk++){
				Rcout << nonDose_PLIN.col(ijk).array().sum() << " ";
			}
			Rcout << " " << endl;
			Rcout << "C++ Note: LOGLIN_non-dose checked ";
			for (int ijk=0;ijk<term_tot;ijk++){
				Rcout << nonDose_LOGLIN.col(ijk).array().sum() << " ";
			}
			Rcout << " " << endl;
		}
		//
		// Calculates the risk for each row
		if (strata_bool){
            Make_Risks_Weighted_Single(modelform, tform, Term_n, totalnum, fir, s_weights, T0, Te, R, Dose, nonDose, TTerm, nonDose_LIN, nonDose_PLIN, nonDose_LOGLIN, nthreads, debugging,KeepConstant,gmix_theta, gmix_term);
        } else {
            Make_Risks_Single(modelform, tform, Term_n, totalnum, fir, T0, Te, R, Dose, nonDose, TTerm, nonDose_LIN, nonDose_PLIN, nonDose_LOGLIN, nthreads, debugging,KeepConstant,gmix_theta, gmix_term);
        }
		//
		// Removes infinite values
		//
		//
		if (R.minCoeff()<=0){
		    ;
		} else if (verbose){
			Rcout << "C++ Note: risk checked ";
			for (int ijk=0;ijk<1;ijk++){
				Rcout << R.col(0).sum() << " ";
			}
			Rcout << " " << endl;
			//
		}
	} else {
		//
		// Calculates the subterm and term values
		//
		Make_subterms( totalnum, Term_n, tform, dfc, fir, T0, Td0, Tdd0, Dose, nonDose, TTerm,  nonDose_LIN, nonDose_PLIN, nonDose_LOGLIN ,beta_0, df0,dint,dslp,nthreads, debugging,KeepConstant);
		// ---------------------------------------------------------
		// Prints off a series of calculations to check at what point values are changing
		// ---------------------------------------------------------
		//
		//
		if (verbose){
			Rcout << "C++ Note: values checked ";
			for (int ijk=0;ijk<totalnum;ijk++){
				Rcout << beta_0[ijk] << " ";
			}
			Rcout << " " << endl;
			Rcout << "C++ Note: sums checked ";
			for (int ijk=0;ijk<totalnum;ijk++){
				Rcout << T0.col(ijk).sum() << " ";
			}
			Rcout << " " << endl;
			Rcout << "C++ Note: derivs checked ";
			for (int ijk=0;ijk<reqrdnum;ijk++){
				Rcout << Td0.col(ijk).sum() << " ";
			}
			Rcout << " " << endl;
			Rcout << "C++ Note: second derivs checked ";
			for (int ijk=0;ijk<reqrdnum;ijk++){
				Rcout << Tdd0.col(ijk*(ijk+1)/2+ijk).sum() << " ";
			}
			Rcout << " " << endl;
			Rcout << "C++ Note: dose checked ";
			for (int ijk=0;ijk<term_tot;ijk++){
				Rcout << Dose.col(ijk).array().sum() << " ";
			}
			Rcout << " " << endl;
			Rcout << "C++ Note: non-dose checked ";
			for (int ijk=0;ijk<term_tot;ijk++){
				Rcout << nonDose.col(ijk).array().sum() << " ";
			}
			Rcout << " " << endl;
			Rcout << "C++ Note: LIN_non-dose checked ";
			for (int ijk=0;ijk<term_tot;ijk++){
				Rcout << nonDose_LIN.col(ijk).array().sum() << " ";
			}
			Rcout << " " << endl;
			Rcout << "C++ Note: PLIN_non-dose checked ";
			for (int ijk=0;ijk<term_tot;ijk++){
				Rcout << nonDose_PLIN.col(ijk).array().sum() << " ";
			}
			Rcout << " " << endl;
			Rcout << "C++ Note: LOGLIN_non-dose checked ";
			for (int ijk=0;ijk<term_tot;ijk++){
				Rcout << nonDose_LOGLIN.col(ijk).array().sum() << " ";
			}
			Rcout << " " << endl;
		}
		// Calculates the risk for each row
		if (strata_bool){
            Make_Risks_Weighted(modelform, tform, Term_n, totalnum, fir, s_weights, T0, Td0, Tdd0, Te, R, Rd, Rdd, Dose, nonDose, TTerm, nonDose_LIN, nonDose_PLIN, nonDose_LOGLIN, RdR, RddR, nthreads, debugging,KeepConstant,gmix_theta, gmix_term);
        } else {
            Make_Risks(modelform, tform, Term_n, totalnum, fir, T0, Td0, Tdd0, Te, R, Rd, Rdd, Dose, nonDose, TTerm, nonDose_LIN, nonDose_PLIN, nonDose_LOGLIN, RdR, RddR, nthreads, debugging,KeepConstant,gmix_theta, gmix_term);
        }
		//
		// Removes infinite values
		RdR = (RdR.array().isFinite()).select(RdR,0);
		RddR = (RddR.array().isFinite()).select(RddR,0);
		//
		//
		if (R.minCoeff()<=0){
		    ;
		} else if (verbose){
			Rcout << "C++ Note: risk checked ";
			for (int ijk=0;ijk<1;ijk++){
				Rcout << R.col(0).sum() << " ";
			}
			Rcout << " " << endl;
			Rcout << "C++ Note: risk1 checked ";
			for (int ijk=0;ijk<reqrdnum;ijk++){
				Rcout << Rd.col(ijk).sum() << " ";
			}
			Rcout << " " << endl;
			Rcout << "C++ Note: risk2 checked ";
			for (int ijk=0;ijk<reqrdnum;ijk++){
				Rcout << Rdd.col(ijk*(ijk+1)/2+ijk).sum() << " ";
			}
			Rcout << " " << endl;
			//
		}
	}
	return;
}


//' Utility function to perform calculation of Log-Likelihood and Deviation for Poisson Omnibus
//'
//' \code{Pois_Dev_LL_Calc} Called to perform repeated term and risk calculations
//' @inheritParams CPP_template
//' @param dev_temp temporary storage for deviation calculation
//' @param dev model deviation
//'
//' @return Updates matrices in place: risk storage matrices
//' @noRd
//'
// [[Rcpp::export]]
void Pois_Dev_LL_Calc(const int& reqrdnum, const int& totalnum, const int& fir, MatrixXd& R, MatrixXd& Rd, MatrixXd& Rdd, VectorXd beta_0 , MatrixXd& RdR, MatrixXd& RddR, vector<double>& Ll, vector<double>& Lld, vector<double>& Lldd, const MatrixXd& PyrC, MatrixXd& dev_temp, const int& nthreads, bool debugging, const IntegerVector& KeepConstant, bool verbose, bool single_bool, int start, int iter_stop, double& dev){
    fill(Ll.begin(), Ll.end(), 0.0);
    if (!single_bool){
        fill(Lld.begin(), Lld.end(), 0.0);
        fill(Lldd.begin(), Lldd.end(), 0.0);
    }
    if (single_bool){
        Poisson_LogLik_Single( nthreads, totalnum, PyrC, R, Ll, debugging);
    } else {
        Poisson_LogLik( nthreads, totalnum, PyrC, R, Rd, Rdd, RdR, RddR, Ll, Lld, Lldd, debugging,KeepConstant);
    }
    //
    dev_temp.col(0) = PyrC.col(0).array() * R.col(0).array();
    dev_temp.col(0) = PyrC.col(1).array() * dev_temp.col(0).array().pow(-1).array();
    dev_temp.col(0) = dev_temp.col(0).array().log().array();
    dev_temp.col(0) = PyrC.col(1).array() * dev_temp.col(0).array();
    dev_temp.col(1) = PyrC.col(1).array() - PyrC.col(0).array() * R.col(0).array();
    //
    dev_temp.col(0) = (dev_temp.col(0).array().isFinite()).select(dev_temp.col(0),0);
    //
    dev_temp.col(0) = dev_temp.col(0).array() - dev_temp.col(1).array();
    dev_temp.col(0) = (2 * dev_temp.col(0).array()).array();//.sqrt();
    dev_temp = (dev_temp.array().isFinite()).select(dev_temp,0);
    dev_temp.col(0) = (R.col(0).array()<0).select(0,dev_temp.col(0));
    dev = dev_temp.col(0).sum(); //deviation calculation is split into steps
    //
    if (single_bool){
        iter_stop=1;
        if (verbose){
            Rcout << "C++ Note: df101 ";//prints the log-likelihoods
            for (int ij=0;ij<reqrdnum;ij++){
                Rcout << Ll[ij] << " ";
            }
            Rcout << " " << endl;
            Rcout << "C++ Note: df104 ";//prints parameter values
            for (int ij=0;ij<totalnum;ij++){
                Rcout << beta_0[ij] << " ";
            }
            Rcout << " " << endl;
            Rcout << "C++ Note: Checking Deviance " << dev << endl;
        }
    } else {
        if (verbose){
            Rcout << "C++ Note: df101 ";//prints the log-likelihoods
            for (int ij=0;ij<reqrdnum;ij++){
                Rcout << Ll[ij] << " ";
            }
            Rcout << " " << endl;
            Rcout << "C++ Note: df102 ";//prints the first derivatives
            for (int ij=0;ij<reqrdnum;ij++){
                Rcout << Lld[ij] << " ";
            }
            Rcout << " " << endl;
            Rcout << "C++ Note: df103 ";//prints the second derivatives
            for (int ij=0;ij<reqrdnum;ij++){
                Rcout << Lldd[ij*reqrdnum+ij] << " ";
            }
            Rcout << " " << endl;
            Rcout << "C++ Note: df104 ";//prints parameter values
            for (int ij=0;ij<totalnum;ij++){
                Rcout << beta_0[ij] << " ";
            }
            Rcout << " " << endl;
            Rcout << "C++ Note: Checking Deviance " << dev << endl;
            Rcout << "C++ Note: df105 ";
            for (int ij=0;ij<reqrdnum;ij++){//prints the newton step value for zero derivative
                Rcout << Lld[ij]/Lldd[ij*reqrdnum+ij] << " ";
            }
            Rcout << " " << endl;
            Rcout << "C++ Note: df106 ";
            for (int ij=0;ij<reqrdnum;ij++){//prints the newton step value for zero log-likelihood
                Rcout << Ll[ij]/Lld[ij] << " ";
            }
            Rcout << " " << endl;
        }
    }
    return;
}

//' Utility function to calculate steps for a likelihood based bound
//'
//' \code{Log_Bound} Called to perform likelihood bound steps
//' @inheritParams CPP_template
//' @param Lstar likelihood goal
//' @param L0 current likelihood
//'
//' @return Updates matrices in place: risk storage matrices
//' @noRd
//'
// [[Rcpp::export]]
void Log_Bound(const MatrixXd& Lldd_mat, const VectorXd& Lld_vec, const double& Lstar, const double& qchi, const double& L0, const int& para_number, const int& nthreads, const int& totalnum, const int& reqrdnum, IntegerVector KeepConstant, const int& term_tot, const int& step, vector<double>& dbeta, const VectorXd& beta_0, bool upper, bool& trouble, bool verbose){
    // starts with solved likelihoods and derivatives
    // store the second derivative as D0
    MatrixXd D0 = Lldd_mat;
    if (step==0){
        if (verbose){
            Rcout << "C++ Note: df201 " << L0 << " " << Lstar << " " << endl;
            Rcout << "C++ Note: df204 ";//prints parameter values
            for (int ij=0;ij<totalnum;ij++){
                Rcout << beta_0[ij] << " ";
            }
            Rcout << " " << endl;
        }
        //Initial step, calculate dom/dbet and h
        VectorXd dOmdBeta = Lldd_mat.col(para_number);
        removeRow(D0, para_number);
        removeColumn(D0, para_number);
        D0 = D0.inverse().matrix();
        dOmdBeta = -1 * dOmdBeta * D0;
        //
        MatrixXd dLdBdO = Lldd_mat.row(para_number).matrix();
        removeColumn(dLdBdO, para_number);
        double h = Lldd_mat(para_number, para_number) - (dLdBdO.matrix() * D0 * dLdBdO.matrix().transpose().matrix())(0,0);
        h = pow(qchi/(-1*h),0.5);
        if (upper){
            h = h/2;
        } else {
            h = h/-2;
        }
        // calculate first step
        int j=0;
        for (int ij=0;ij<totalnum;ij++){
            if (KeepConstant[ij]==0){
                int pij_ind = ij - sum(head(KeepConstant,ij));
                if (pij_ind == para_number){
                    dbeta[ij] = h;
                } else {
                    dbeta[ij] = h * dOmdBeta(j);
                    j=j+1;
                }
            }
        }
//        Rcout << "finished step 0" << endl;
    } else {
        //Td0 = MatrixXd::Zero(df0.rows(), reqrdnum);
        MatrixXd G = MatrixXd::Zero(reqrdnum, reqrdnum);
        VectorXd v = VectorXd::Zero(reqrdnum);
        v[para_number] = L0 - Lstar;
        G.row(para_number) = Lld_vec;
        for (int j=0; j<reqrdnum;j++){
            if (j!=para_number){
                G.row(j) = D0.row(j);
                v[j] = Lld_vec[j];
            }
        }
        // At this point, we have the standard newton-raphson equation defined
        if (verbose){
            Rcout << "C++ Note: df201 " << L0 << " " << Lstar << " " << endl;
            Rcout << "C++ Note: df202 ";//prints the first derivatives
            for (int ij=0;ij<reqrdnum;ij++){
                Rcout << v[ij] << " ";
            }
            Rcout << " " << endl;
            Rcout << "C++ Note: df203 ";//prints the second derivatives
            for (int ij=0;ij<reqrdnum;ij++){
                Rcout << G(ij, ij) << " ";
            }
            Rcout << " " << endl;
            Rcout << "C++ Note: df204 ";//prints parameter values
            for (int ij=0;ij<totalnum;ij++){
                Rcout << beta_0[ij] << " ";
            }
            Rcout << " " << endl;
            Rcout << "C++ Note: Second Derivative Determinant: " << G.determinant() << endl;
        }
        //
        if (abs(G.determinant()) < 1e-6){
            // The inverted matrix does not exist
            for (int ij=0;ij<totalnum;ij++){
                if (KeepConstant[ij]==0){
                    int pij_ind = ij - sum(head(KeepConstant,ij));
                    dbeta[ij] = -v[pij_ind]/G(pij_ind,pij_ind);
                }
            }
        } else {
            G = G.inverse().matrix();
            v = G.matrix() * v.matrix();
    //        Rcout << "reached G" << endl;
            VectorXd g1 = G.col(para_number);
            // we now must solve for the roots
            double as2 = g1.matrix().transpose() * D0 * g1.matrix();
            double bs1 = 2*v.matrix().transpose() *D0 * g1.matrix() - 2;
            double cs0 = v.matrix().transpose() * D0 * v.matrix();
//            Rcout << as2 << "s^2 + " << bs1 << "s + " << cs0 << endl;
            //
            if (pow(bs1,2)-4*as2*cs0 >= 0){
                double s0 = pow(bs1,2)-4*as2*cs0;
                double s1 = (-bs1 - pow(s0, 0.5))/(2*as2);
                s0 = (-bs1 + pow(s0, 0.5))/(2*as2);
                // check which is closer
                double s00 = (v + s0*g1).matrix().transpose() * D0 * (v + s0*g1).matrix();
                double s11 = (v + s1*g1).matrix().transpose() * D0 * (v + s1*g1).matrix();
                //
                if (abs(s00)<abs(s11)){
                    //s1 is further away
                    for (int ij=0;ij<totalnum;ij++){
                        if (KeepConstant[ij]==0){
                            int pij_ind = ij - sum(head(KeepConstant,ij));
                            dbeta[ij] = -v[pij_ind] - g1[pij_ind]*s0;
                        }
                    }
                } else {
                    //s0 is further away
                    for (int ij=0;ij<totalnum;ij++){
                        if (KeepConstant[ij]==0){
                            int pij_ind = ij - sum(head(KeepConstant,ij));
                            dbeta[ij] = -v[pij_ind] - g1[pij_ind]*s1;
                        }
                    }
                }
            } else {
                // there are no real solutions, needs a more conservative step?
                // currently will just use the extrema
                double s0 = -bs1/2/as2;
                for (int ij=0;ij<totalnum;ij++){
                    if (KeepConstant[ij]==0){
                        int pij_ind = ij - sum(head(KeepConstant,ij));
                        dbeta[ij] = -v[pij_ind] - g1[pij_ind]*s0;
                    }
                }
            }
        }
        if (verbose){
            Rcout << "C++ Note: df205 ";//prints parameter values
            for (int ij=0;ij<totalnum;ij++){
                Rcout << dbeta[ij] << " ";
            }
            Rcout << " " << endl;
        }
//        Rcout << as2 << ", " << bs1 << ", " << cs0 << endl;
    }
    return;
}

//' Utility function to calculate the change to make each iteration
//'
//' \code{Calc_Change_trouble} Called to update the parameter changes, Uses log-likelihoods and control parameters, Applies newton steps and change limitations    
//' @inheritParams CPP_template
//'
//' @return Updates matrices in place: parameter change matrix
//' @noRd
//'
// [[Rcpp::export]]
void Calc_Change_trouble(const int& para_number, const int& nthreads, const int& totalnum, const int& der_iden, const double& dbeta_cap, const double& dose_abs_max, const double& lr, const double& abs_max, const vector<double>& Ll, const vector<double>& Lld, const vector<double>& Lldd, vector<double>& dbeta,const StringVector&   tform, const double& dint, const double& dslp, IntegerVector KeepConstant_trouble, bool debugging){
    int kept_covs = totalnum - sum(KeepConstant_trouble);
    NumericVector Lldd_vec(kept_covs * kept_covs);
    NumericVector Lld_vec(kept_covs);
    #ifdef _OPENMP
    #pragma omp parallel for schedule(dynamic) num_threads(nthreads)
    #endif
    for (int ijk=0;ijk<kept_covs*(kept_covs+1)/2;ijk++){
        int ij = 0;
        int jk = ijk;
        while (jk>ij){
            ij++;
            jk-=ij;
        }
        int ij0 = ij;
        int jk0 = jk;
        if (ij >= para_number){
            ij0++;
        }
        if (jk >= para_number){
            jk0++;
        }
        Lldd_vec[jk * kept_covs + ij]=Lldd[jk0 * (kept_covs+1) + ij0];
        if (ij==jk){
            Lld_vec[ij]=Lld[ij0];
        } else {
            Lldd_vec[ij * kept_covs + jk]=Lldd_vec[jk0 * (kept_covs+1) + ij0];
        }
    }
    //
    //
    Lldd_vec.attr("dim") = Dimension(kept_covs, kept_covs);
    const Map<MatrixXd> Lldd_mat(as<Map<MatrixXd> >(Lldd_vec));
    const Map<VectorXd> Lld_mat(as<Map<VectorXd> >(Lld_vec));
    VectorXd Lldd_solve0 = Lldd_mat.colPivHouseholderQr().solve(-1*Lld_mat);
    VectorXd Lldd_solve = VectorXd::Zero(totalnum);
    for (int ij=0;ij<totalnum;ij++){
        if (KeepConstant_trouble[ij]==0){
            int pij_ind = ij - sum(head(KeepConstant_trouble,ij));
            Lldd_solve(ij) = Lldd_solve0(pij_ind);
        }
    }
    //
    #ifdef _OPENMP
    #pragma omp parallel for schedule(dynamic) num_threads(nthreads)
    #endif
    for (int ijk=0;ijk<totalnum;ijk++){
        if (KeepConstant_trouble[ijk]==0){
            int pjk_ind = ijk - sum(head(KeepConstant_trouble,ijk));
            if (isnan(Lldd_solve(ijk))){
                if (Lldd[pjk_ind*kept_covs+pjk_ind] != 0 ){
                    dbeta[ijk] = -lr * Lld[pjk_ind] / Lldd[pjk_ind*kept_covs+pjk_ind];
                } else {
                    dbeta[ijk] = 0;
                }
            } else {
                dbeta[ijk] = lr * Lldd_solve(ijk);//-lr * Lld[ijk] / Lldd[ijk*totalnum+ijk];
            }
//            //
            if ((tform[ijk]=="lin_quad_int")||(tform[ijk]=="lin_exp_int")||(tform[ijk]=="step_int")||(tform[ijk]=="lin_int")){ //the threshold values use different maximum deviation values
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
    }
    return;
}


//' Utility function to calculate Cox Log-Likelihood and derivatives
//'
//' \code{Simplified_Inform_Matrix} Called to update log-likelihoods, Uses list of event rows, risk matrices, and repeated sums, Sums the log-likelihood contribution from each event time
//' @inheritParams CPP_template
//'
//' @return Updates matrices in place: Log-likelihood vectors/matrix
//' @noRd
//'
// [[Rcpp::export]]
void Simplified_Inform_Matrix(const int& nthreads,const IntegerMatrix& RiskFail, const vector<string>&  RiskGroup, const int& totalnum, const int& ntime, const MatrixXd& R, const MatrixXd& Rd, const MatrixXd& Rdd, const MatrixXd& RdR, const MatrixXd& RddR,const MatrixXd& Rls1,const MatrixXd& Rls2,const MatrixXd& Rls3,const MatrixXd& Lls1,const MatrixXd& Lls2,const MatrixXd& Lls3, vector<double>& InMa, bool debugging,string ties_method, const IntegerVector& KeepConstant){
    int reqrdnum = totalnum - sum(KeepConstant);
    #ifdef _OPENMP
    #pragma omp declare reduction(vec_double_plus : std::vector<double> : \
        std::transform(omp_out.begin(), omp_out.end(), omp_in.begin(), omp_out.begin(), std::plus<double>())) \
        initializer(omp_priv = omp_orig)
    #pragma omp parallel for schedule(dynamic) num_threads(nthreads) reduction(vec_double_plus:InMa) collapse(2)
    #endif
    for (int ijk=0;ijk<reqrdnum*(reqrdnum+1)/2;ijk++){//performs log-likelihood calculations for every derivative combination and risk group
        for (int j=0;j<ntime;j++){
			//Rcout << ijk << ", " << j << ", " << " start" << endl;
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
			//Rcout << 0 << endl;
            //Ld << R.block(RiskFail(j,0),0,dj,1), RdR.block(RiskFail(j,0),ij,dj,1), RdR.block(RiskFail(j,0),jk,dj,1) ,RddR.block(RiskFail(j,0),ijk,dj,1);//rows with events
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
            //
            double Ld1;
            double Ld3;
            //
            MatrixXd temp1 = MatrixXd::Zero(dj,1);
            MatrixXd temp2 = MatrixXd::Zero(dj,1);
            // calculates the right-hand side terms
            temp1 = Ldm.col(1).array() * (Ldm.col(0).array().pow(-1).array());
            temp2 = Ldm.col(2).array() * (Ldm.col(0).array().pow(-1).array());
            temp1 = temp1.array() * temp2.array();
            Rs3 = (temp1.array().isFinite()).select(temp1,0).sum();
            //
            vector<int> InGroup;
            string Groupstr = RiskGroup[j];
            stringstream ss(Groupstr);
            //
            //
            //
            for (int i; ss >> i;) {
                InGroup.push_back(i);    
                if (ss.peek() == ',')
                    ss.ignore();
            }
            Ld3 = 0;
            Ld1 = 0;
			//Rcout << 1 << endl;
            //Now has the grouping pairs
            for (vector<double>::size_type i = 0; i < InGroup.size()-1; i=i+2){
				//Rcout << 2 << ", " << i << endl;
                Ld3 += (RdR.block(InGroup[i]-1,ij,InGroup[i+1]-InGroup[i]+1,1).array() * RdR.block(InGroup[i]-1,jk,InGroup[i+1]-InGroup[i]+1,1).array()* R.block(InGroup[i]-1,0,InGroup[i+1]-InGroup[i]+1,1).array()).sum();
                Ld1 += R.block(InGroup[i]-1,0,InGroup[i+1]-InGroup[i]+1,1).sum();
            }
            Ld3 = Ld3 / Ld1;
			//Rcout << 3 << endl;
            //
            InMa[ij*reqrdnum+jk] += Ld3 - Rs3; //sums the log-likelihood and derivatives
			//Rcout << 4 << endl;
        }
    }
    #ifdef _OPENMP
    #pragma omp parallel for num_threads(nthreads)
    #endif
    for (int ijk=0;ijk<reqrdnum*(reqrdnum+1)/2;ijk++){//fills second-derivative matrix
        int ij = 0;
        int jk = ijk;
        while (jk>ij){
            ij++;
            jk-=ij;
        }
        InMa[jk*reqrdnum+ij] = InMa[ij*reqrdnum+jk];
    }
    return;
}

//' Utility function to calculate Cox Log-Likelihood and derivatives
//'
//' \code{Simplified_Inform_Matrix_STRATA} Called to update log-likelihoods, Uses list of event rows, risk matrices, and repeated sums, Sums the log-likelihood contribution from each event time
//' @inheritParams CPP_template
//'
//' @return Updates matrices in place: Log-likelihood vectors/matrix
//' @noRd
//'
// [[Rcpp::export]]
void Simplified_Inform_Matrix_STRATA(const int& nthreads,const IntegerMatrix& RiskFail, const StringMatrix&  RiskGroup, const int& totalnum, const int& ntime, const MatrixXd& R, const MatrixXd& Rd, const MatrixXd& Rdd, const MatrixXd& RdR, const MatrixXd& RddR,const MatrixXd& Rls1,const MatrixXd& Rls2,const MatrixXd& Rls3,const MatrixXd& Lls1,const MatrixXd& Lls2,const MatrixXd& Lls3, vector<double>& InMa, bool debugging,string ties_method, NumericVector& STRATA_vals, const IntegerVector& KeepConstant){
    int reqrdnum = totalnum - sum(KeepConstant);
    #ifdef _OPENMP
    #pragma omp declare reduction(vec_double_plus : std::vector<double> : \
        std::transform(omp_out.begin(), omp_out.end(), omp_in.begin(), omp_out.begin(), std::plus<double>())) \
        initializer(omp_priv = omp_orig)
    #pragma omp parallel for schedule(dynamic) num_threads(nthreads) reduction(vec_double_plus:InMa) collapse(3)
    #endif
    for (int ijk=0;ijk<reqrdnum*(reqrdnum+1)/2;ijk++){//performs log-likelihood calculations for every derivative combination and risk group
        for (int j=0;j<ntime;j++){
            for (int s_ij=0;s_ij<STRATA_vals.size();s_ij++){
                int ij = 0;
                int jk = ijk;
                while (jk>ij){
                    ij++;
                    jk-=ij;
                }
                double Rs1 = Rls1(j,0);
                double Rs2 =  Rls2(j,ij*STRATA_vals.size() + s_ij);
                double Rs2t = Rls2(j,jk*STRATA_vals.size() + s_ij);
                double Rs3 =  Rls3(j,ijk*STRATA_vals.size() + s_ij);
                //
                int dj = RiskFail(j,2*s_ij + 1)-RiskFail(j,2*s_ij + 0)+1;
                if (RiskFail(j,2*s_ij + 1)>-1){
                    MatrixXd Ld = MatrixXd::Zero(dj,4);
                    //Ld << R.block(RiskFail(j,2*s_ij),0,dj,1), RdR.block(RiskFail(j,2*s_ij),ij,dj,1), RdR.block(RiskFail(j,2*s_ij),jk,dj,1) ,RddR.block(RiskFail(j,2*s_ij),ijk,dj,1);//rows with events
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
                    //
                    double Ld1;
                    double Ld3;
                    //
                    MatrixXd temp1 = MatrixXd::Zero(dj,1);
                    MatrixXd temp2 = MatrixXd::Zero(dj,1);
                    // calculates the right-hand side terms
                    temp1 = Ldm.col(1).array() * (Ldm.col(0).array().pow(-1).array());
                    temp2 = Ldm.col(2).array() * (Ldm.col(0).array().pow(-1).array());
                    temp1 = temp1.array() * temp2.array();
                    Rs3 = (temp1.array().isFinite()).select(temp1,0).sum();
                    //
                    vector<int> InGroup;
                    string Groupstr = as<std::string>(RiskGroup(j,s_ij));
                    stringstream ss(Groupstr);
                    for (int i; ss >> i;) {
                        InGroup.push_back(i);    
                        if (ss.peek() == ',')
                            ss.ignore();
                    }
                    Ld3 = 0;
                    Ld1 = 0;
                    //Now has the grouping pairs
                    for (vector<double>::size_type i = 0; i < InGroup.size()-1; i=i+2){
                        Ld3 += (RdR.block(InGroup[i]-1,ij,InGroup[i+1]-InGroup[i]+1,1).array() * RdR.block(InGroup[i]-1,jk,InGroup[i+1]-InGroup[i]+1,1).array()* R.block(InGroup[i]-1,0,InGroup[i+1]-InGroup[i]+1,1).array()).sum();
                        Ld1 += R.block(InGroup[i]-1,0,InGroup[i+1]-InGroup[i]+1,1).sum();
                    }
                    Ld3 = Ld3 / Ld1;
                    //
                    InMa[ij*reqrdnum+jk] += Ld3 - Rs3; //sums the log-likelihood and derivatives
                }
            }
        }
    }
    #ifdef _OPENMP
    #pragma omp parallel for num_threads(nthreads)
    #endif
    for (int ijk=0;ijk<reqrdnum*(reqrdnum+1)/2;ijk++){//fills second-derivative matrix
        int ij = 0;
        int jk = ijk;
        while (jk>ij){
            ij++;
            jk-=ij;
        }
        InMa[jk*reqrdnum+ij] = InMa[ij*reqrdnum+jk];
    }
    return;
}
