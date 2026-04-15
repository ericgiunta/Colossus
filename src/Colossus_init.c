#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME: 
   Check these declarations against the C/Fortran source code.
*/

/* .Call calls */
extern SEXP _Colossus_Assigned_Event_Poisson_transition(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _Colossus_caco_Omnibus_transition(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _Colossus_cox_ph_multidose_Omnibus_transition(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _Colossus_cox_ph_Omnibus_Bounds_transition(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _Colossus_cox_ph_Omnibus_CurveSearch_transition(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _Colossus_cox_ph_Omnibus_transition(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _Colossus_logist_Omnibus_transition(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _Colossus_OMP_Check(void);
extern SEXP _Colossus_Plot_Omnibus_transition(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _Colossus_pois_multidose_Omnibus_transition(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _Colossus_pois_Omnibus_Bounds_transition(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _Colossus_pois_Omnibus_CurveSearch_transition(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _Colossus_pois_Omnibus_transition(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _Colossus_pois_Residual_transition(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _Colossus_Write_Time_Dep(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);

static const R_CallMethodDef CallEntries[] = {
    {"_Colossus_Assigned_Event_Poisson_transition",     (DL_FUNC) &_Colossus_Assigned_Event_Poisson_transition,     14},
    {"_Colossus_caco_Omnibus_transition",               (DL_FUNC) &_Colossus_caco_Omnibus_transition,               16},
    {"_Colossus_cox_ph_multidose_Omnibus_transition",   (DL_FUNC) &_Colossus_cox_ph_multidose_Omnibus_transition,   20},
    {"_Colossus_cox_ph_Omnibus_Bounds_transition",      (DL_FUNC) &_Colossus_cox_ph_Omnibus_Bounds_transition,      17},
    {"_Colossus_cox_ph_Omnibus_CurveSearch_transition", (DL_FUNC) &_Colossus_cox_ph_Omnibus_CurveSearch_transition, 17},
    {"_Colossus_cox_ph_Omnibus_transition",             (DL_FUNC) &_Colossus_cox_ph_Omnibus_transition,             17},
    {"_Colossus_logist_Omnibus_transition",             (DL_FUNC) &_Colossus_logist_Omnibus_transition,             14},
    {"_Colossus_OMP_Check",                             (DL_FUNC) &_Colossus_OMP_Check,                              0},
    {"_Colossus_Plot_Omnibus_transition",               (DL_FUNC) &_Colossus_Plot_Omnibus_transition,               16},
    {"_Colossus_pois_multidose_Omnibus_transition",     (DL_FUNC) &_Colossus_pois_multidose_Omnibus_transition,     19},
    {"_Colossus_pois_Omnibus_Bounds_transition",        (DL_FUNC) &_Colossus_pois_Omnibus_Bounds_transition,        16},
    {"_Colossus_pois_Omnibus_CurveSearch_transition",   (DL_FUNC) &_Colossus_pois_Omnibus_CurveSearch_transition,   16},
    {"_Colossus_pois_Omnibus_transition",               (DL_FUNC) &_Colossus_pois_Omnibus_transition,               16},
    {"_Colossus_pois_Residual_transition",              (DL_FUNC) &_Colossus_pois_Residual_transition,              14},
    {"_Colossus_Write_Time_Dep",                        (DL_FUNC) &_Colossus_Write_Time_Dep,                        10},
    {NULL, NULL, 0}
};

void R_init_Colossus(DllInfo *dll)
{
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}

