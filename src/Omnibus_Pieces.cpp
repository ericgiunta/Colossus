#include <RcppEigen.h>
#ifdef _OPENMP
#include <omp.h>
#endif
#include "Omnibus_Pieces.h"
#include "Calc_Repeated.h"
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

