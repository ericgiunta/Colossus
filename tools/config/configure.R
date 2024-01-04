# Prepare your package for installation here.
# Use 'define()' to define configuration variables.
# Use 'configure_file()' to substitute configuration values.

get_os <- function(){
  sysinf <- Sys.info()
  if (!is.null(sysinf)){
    os <- sysinf['sysname']
    if (os == 'Darwin')
      os <- "osx"
  } else { ## mystery machine
    os <- .Platform$OS.type
    if (grepl("^darwin", R.version$os))
      os <- "osx"
    if (grepl("linux-gnu", R.version$os))
      os <- "linux"
  }
  tolower(os)
}

os <- get_os()

if (os=="linux"){
    define(PKG_CXXFLAGS = "PKG_CXXFLAGS=-fopenmp")
    define(PKG_LIBS = 'PKG_LIBS = `$(R_HOME)/bin/Rscript -e "Rcpp:::LdFlags()"` $(LAPACK_LIBS) $(BLAS_LIBS) $(FLIBS) -fopenmp')
    define(CPPFLAGS = "#CPPFLAGS += -Xclang -fopenmp")
    define(LDFLAGS = "#LDFLAGS += -lomp")
    configure_file("src/Makevars.in")
} else if (os=="osx"){
    define(PKG_CXXFLAGS = "#PKG_CXXFLAGS=-fopenmp")
    define(PKG_LIBS = '#PKG_LIBS = `$(R_HOME)/bin/Rscript -e "Rcpp:::LdFlags()"` $(LAPACK_LIBS) $(BLAS_LIBS) $(FLIBS) -fopenmp')
    define(CPPFLAGS = "CPPFLAGS += -Xclang -fopenmp")
    define(LDFLAGS = "LDFLAGS += -lomp")
    configure_file("src/Makevars.in")
} else {
    print(paste("OS",os,sep=" "))
    define(PKG_CXXFLAGS = "PKG_CXXFLAGS=-fopenmp")
    define(PKG_LIBS = 'PKG_LIBS = `$(R_HOME)/bin/Rscript -e "Rcpp:::LdFlags()"` $(LAPACK_LIBS) $(BLAS_LIBS) $(FLIBS) -fopenmp')
    configure_file("src/Makevars.win.in")
}
