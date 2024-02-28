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
