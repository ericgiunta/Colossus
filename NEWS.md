# Colossus 0.9

* Added a `NEWS.md` file to track changes to the package.

# Colossus 1.0.0

* Initial submission
* C++ code modified to not apply OpenMP code if OpenMP isn't detected, to resolve MacOS installation failures

# Colossus 1.0.1

* Configuration improved to detect compiler information
* Linux configuration depends on the system wide c++ default and the compiler used to compile R
* OpenMP support is not used if the c++ default or the R compiler is clang

# Colossus 1.0.2

* utility checks updated to check that keep_constant is 0/1 values
* code will fail if keep_constant is not integer valued and 0/1, with explanation of why

# Colossus 1.0.3

* configuration script libraries moved from Suggested: to Imports:
* configuration script functions moved to non-exported functions in Colossus, circumvents note about imported libraries not being used and may be later used to provide the user with a function that informs them of why OpenMP is/isn't supported

# Colossus 1.0.4

* utility checks updated to check term numbers and subterm types in R side
* code will fail if term numbers are not integers, term numbers are missing, or subterm types are

# Colossus 1.0.5

* compilation flags changed to macros

# Colossus 1.1.0

* Default Makevars is now fully portable
* By default only windows uses OpenMP now
* The GitHub version will include instructions on activating configuration to use OpenMP on gcc Linux

# Colossus 1.1.1

* ggplot2 no longer required, now optional
* additional testing added for coverage
* default Makevars added via bash script

# Colossus 1.1.2

* Log-likelihood bound functionality added
* subject to usual convergence issues, manual search option

# Colossus 1.1.3

* Cox regression now removes rows that end before first event and start after last event
* Cox regression now sets constant rows to be constant, helps with aliasing
* Tests now use sink() to avoid printing as much excessive output to console. Tests now consolidated further.

# Colossus 1.1.4.1

* Cox plotting functions now return the tables used for plots (only last plot table returned)
* Plotting vignette updated to include more details and plots
* survival package listed as suggested for plotting vignette

# Colossus 1.1.4.2

* R errors and warnings sent to stop() and warning(). C++ errors and warnings still controlled by verbosity argument
* ggsave defaults to width/height = 7
* Updates started to fix possible OpenMP issues with fedora 36 running clang 18
* Unable to have debug printing as an option and cover the c++ files with testing, but still have test output be readable. debug output removed.

# Colossus 1.1.5

* Started adding simplifications to allow for faster iterations
* Added simplification for linear ERR model
* Started gradient descent code
* Started external rate comparison options
* Added person-count and person-time table generation code and vignette

# Colossus 1.1.5.5

* Added CoxCurveSolver function to solve likelihood boundaries via bisection method

# Colossus 1.1.6

* MacOS testing with OpenMP finished. MacOS use with OpenMP officially checked.
* By default the only systems that is forced to use single thread is linux using clang. This can be turned off by setting the "R_COLOSSUS_NOT_CRAN" environment variable.

# Colossus 1.1.7

* Gradient descent algorithms tested further and presented in vignette
* Multiple realization function tested further

# Colossus 1.1.8

* CurveSolve functions converted to c++ functions
* Testing scaled back to take up less time

# Colossus 1.1.9

* Cox based functions switched covariance matrix calculation from negative inverse of log-likelihood second derivative, to expected information matrix.

# Colossus 1.1.10

* Cox based functions updated to improve speed
* Additional CurveSolve output provided to give final window width and final step

# Colossus 1.2

* Both Cox and Poisson functions now all return the expected information matrix derived covariance

# Colossus 1.2.1

* Matched Case-Control base code and equation vignette added

# Colossus 1.2.2

* Additional code cleanup and minor utility function speed improvements

# Colossus 1.3.0

* Multiple realization code updated to improve speed
* lingering debugging variables removed: fir and der_iden

# Colossus 1.3.1

* Convergence check performed more often
* Checks the actual maximum step taken to compare to threshold
* Previous versions may have run more iterations than necessary to meet derivative and step size thresholds

# Colossus 1.4.0

* Switch to formula inputs and model classes
* Removed unused functions, to simplify documentation

# Colossus 1.4.1

* Formula input now allows more general applications of factor, ns(), bs(), I(var^n), interaction, etc.

# Colossus 1.4.2

* Option added to normalize covariates, either scaled by the mean or maximum values.
* Gradient descent option now calculates standard error and covariance through second derivative after regression.

# Colossus 1.4.3

* Logistic regression added
* Newton step calculation now checks the predicted change in score, and moves in opposite direction if the expected score is worse
