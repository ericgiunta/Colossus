# Prepare your package for installation here.
# Use 'define()' to define configuration variables.
# Use 'configure_file()' to substitute configuration values.

if (Sys.info()['sysname']=="Linux"){
    define(PKG_CXXFLAGS = "PKG_CXXFLAGS=-fopenmp")
    define(PKG_LIBS = 'PKG_LIBS = `$(R_HOME)/bin/Rscript -e "Rcpp:::LdFlags()"` $(LAPACK_LIBS) $(BLAS_LIBS) $(FLIBS) -fopenmp')
    define(CPPFLAGS = "#CPPFLAGS += -Xclang -fopenmp")
    define(LDFLAGS = "#LDFLAGS += -lomp")
} else if (Sys.info()['sysname']=="Darwin"){
    define(PKG_CXXFLAGS = "#PKG_CXXFLAGS=-fopenmp")
    define(PKG_LIBS = '#PKG_LIBS = `$(R_HOME)/bin/Rscript -e "Rcpp:::LdFlags()"` $(LAPACK_LIBS) $(BLAS_LIBS) $(FLIBS) -fopenmp')
    define(CPPFLAGS = "CPPFLAGS += -Xclang -fopenmp")
    define(LDFLAGS = "LDFLAGS += -lomp")
} else {
    print(paste("OS",Sys.info()['sysname'],sep=" "))
    define(PKG_CXXFLAGS = "#PKG_CXXFLAGS=-fopenmp")
    define(PKG_LIBS = '#PKG_LIBS = `$(R_HOME)/bin/Rscript -e "Rcpp:::LdFlags()"` $(LAPACK_LIBS) $(BLAS_LIBS) $(FLIBS) -fopenmp')
    define(CPPFLAGS = "#CPPFLAGS += -Xclang -fopenmp")
    define(LDFLAGS = "#LDFLAGS += -lomp")
}

configure_file("src/Makevars.in")
