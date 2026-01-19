# Package index

## All functions

- [`CaseControlRun()`](CaseControlRun.md) : Fully runs a case-control
  regression model, returning the model and results
- [`ColossusCoxSurv()`](ColossusCoxSurv.md) : Interprets basic cox
  survival formula RHS
- [`ColossusLogitSurv()`](ColossusLogitSurv.md) : Interprets basic
  logistic survival formula RHS with no grouping
- [`ColossusPoisSurv()`](ColossusPoisSurv.md) : Interprets basic poisson
  survival formula RHS
- [`CoxRun()`](CoxRun.md) : Fully runs a cox or fine-gray regression
  model, returning the model and results
- [`CoxRunMulti()`](CoxRunMulti.md) : Fully runs a cox or fine-gray
  regression model with multiple column realizations, returning the
  model and results
- [`Date_Shift()`](Date_Shift.md) : Automates creating a date difference
  column
- [`EventAssignment()`](EventAssignment.md) : Generic background/excess
  event calculation function
- [`EventAssignment(`*`<default>`*`)`](EventAssignment.default.md) :
  Predicts how many events are due to baseline vs excess
- [`EventAssignment(`*`<poisres>`*`)`](EventAssignment.poisres.md) :
  Predicts how many events are due to baseline vs excess for a completed
  poisson model
- [`EventAssignment(`*`<poisresbound>`*`)`](EventAssignment.poisresbound.md)
  : Predicts how many events are due to baseline vs excess for a
  completed poisson likelihood boundary regression
- [`Event_Count_Gen()`](Event_Count_Gen.md) : uses a table, list of
  categories, and list of event summaries to generate person-count
  tables
- [`Event_Time_Gen()`](Event_Time_Gen.md) : uses a table, list of
  categories, list of summaries, list of events, and person-year
  information to generate person-time tables
- [`Joint_Multiple_Events()`](Joint_Multiple_Events.md) : Automates
  creating data for a joint competing risks analysis
- [`LikelihoodBound()`](LikelihoodBound.md) : Generic likelihood
  boundary calculation function
- [`LikelihoodBound(`*`<coxres>`*`)`](LikelihoodBound.coxres.md) :
  Calculates the likelihood boundary for a completed cox model
- [`LikelihoodBound(`*`<default>`*`)`](LikelihoodBound.default.md) :
  Generic likelihood boundary calculation function, default option
- [`LikelihoodBound(`*`<poisres>`*`)`](LikelihoodBound.poisres.md) :
  Calculates the likelihood boundary for a completed Poisson model
- [`Likelihood_Ratio_Test()`](Likelihood_Ratio_Test.md) : Defines the
  likelihood ratio test
- [`Linked_Dose_Formula()`](Linked_Dose_Formula.md) : Calculates Full
  Parameter list for Special Dose Formula
- [`Linked_Lin_Exp_Para()`](Linked_Lin_Exp_Para.md) : Calculates The
  Additional Parameter For a linear-exponential formula with known
  maximum
- [`LogisticRun()`](LogisticRun.md) : Fully runs a logistic regression
  model, returning the model and results
- [`OMP_Check()`](OMP_Check.md) : Checks the OMP flag
- [`PoisRun()`](PoisRun.md) : Fully runs a poisson regression model,
  returning the model and results
- [`PoisRunJoint()`](PoisRunJoint.md) : Fully runs a joint poisson
  regression model, returning the model and results
- [`PoisRunMulti()`](PoisRunMulti.md) : Fully runs a poisson regression
  model with multiple column realizations, returning the model and
  results
- [`RelativeRisk()`](RelativeRisk.md) : Generic relative risk
  calculation function
- [`RelativeRisk(`*`<coxres>`*`)`](RelativeRisk.coxres.md) : Calculates
  hazard ratios for a reference vector
- [`RelativeRisk(`*`<default>`*`)`](RelativeRisk.default.md) : Generic
  relative risk calculation function, default option
- [`Replace_Missing()`](Replace_Missing.md) : Automatically assigns
  missing values in listed columns
- [`Residual()`](Residual.md) : Generic Residual calculation function
- [`Residual(`*`<default>`*`)`](Residual.default.md) : Generic Residual
  calculation function, default option
- [`Residual(`*`<poisres>`*`)`](Residual.poisres.md) : Calculates the
  Residuals for a completed poisson model
- [`System_Version()`](System_Version.md) : Checks OS, compilers, and
  OMP
- [`Time_Since()`](Time_Since.md) : Automates creating a date since a
  reference column
- [`factorize()`](factorize.md) : Splits a parameter into factors
- [`gen_time_dep()`](gen_time_dep.md) : Applies time dependence to
  parameters
- [`get_form()`](get_form.md) : Interprets a Colossus formula and makes
  necessary changes to data
- [`get_form_joint()`](get_form_joint.md) : Interprets a Poisson joint
  formula and makes necessary changes to data
- [`plot(`*`<coxres>`*`)`](plot.coxres.md) : Performs Cox Proportional
  Hazard model plots
- [`plotMartingale()`](plotMartingale.md) : Generic Martingale Residual
  Plotting function
- [`plotMartingale(`*`<coxres>`*`)`](plotMartingale.coxres.md) :
  Performs Cox Proportional Hazard model martingale residual plots
- [`plotMartingale(`*`<default>`*`)`](plotMartingale.default.md) :
  Generic Martingale Residual Plotting function, default option
- [`plotRisk()`](plotRisk.md) : Generic Risk Plotting function
- [`plotRisk(`*`<coxres>`*`)`](plotRisk.coxres.md) : Performs Cox
  Proportional Hazard model hazard ratio plots
- [`plotRisk(`*`<default>`*`)`](plotRisk.default.md) : Generic Risk
  Plotting function, default option
- [`plotSchoenfeld()`](plotSchoenfeld.md) : Generic Schoenfeld Residual
  Plotting function
- [`plotSchoenfeld(`*`<coxres>`*`)`](plotSchoenfeld.coxres.md) :
  Performs Cox Proportional Hazard model schoenfeld residual plots
- [`plotSchoenfeld(`*`<default>`*`)`](plotSchoenfeld.default.md) :
  Generic Schoenfeld Residual Plotting function, default option
- [`plotSurvival()`](plotSurvival.md) : Generic Survival Plotting
  function
- [`plotSurvival(`*`<coxres>`*`)`](plotSurvival.coxres.md) : Performs
  Cox Proportional Hazard model survival plots
- [`plotSurvival(`*`<default>`*`)`](plotSurvival.default.md) : Generic
  Survival Plotting function, default option
- [`print(`*`<caseconres>`*`)`](print.caseconres.md) : Prints a
  case-control regression output clearly
- [`print(`*`<coxres>`*`)`](print.coxres.md) : Prints a cox regression
  output clearly
- [`print(`*`<coxresbound>`*`)`](print.coxresbound.md) : Prints a cox
  likelihood boundary regression output clearly
- [`print(`*`<logitres>`*`)`](print.logitres.md) : Prints a logistic
  regression output clearly
- [`print(`*`<poisres>`*`)`](print.poisres.md) : Prints a poisson
  regression output clearly
- [`print(`*`<poisresbound>`*`)`](print.poisresbound.md) : Prints a
  poisson likelihood boundary regression output clearly
