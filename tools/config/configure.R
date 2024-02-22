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

gcc_version <- function() {
  out <- tryCatch(processx::run("c++", "-v", stderr_to_stdout = TRUE),
                  error = function(cnd) list(stdout = ""))
  out0 <- stringr::str_match(out$stdout, "gcc version")[1]
  if (!is.na(out0)){
  	out <- "gcc"
  } else {
    out0 <- stringr::str_match(out$stdout, "clang version")[1]
    if (!is.na(out0)){
      out <- "clang"
    } else {
      out <- out$stdout
    }
  }
  out
}

Rcomp_version <- function() {
  out <- callr::rcmd("config","CC")
  out0 <- stringr::str_match(out$stdout, "clang")[1]
  if (!is.na(out0)){
  	out <- "clang"
  } else {
    out0 <- stringr::str_match(out$stdout, "gcc")[1]
    if (!is.na(out0)){
      out <- "gcc"
    } else {
      out <- out$stdout
    }
  }
  out
}

Rcpp_version <- function() {
  out <- tryCatch(processx::run("head", "~/.R/Makevars", stderr_to_stdout = TRUE),
                  error = function(cnd) list(stdout = ""))
  out0 <- stringr::str_match(out$stdout, "clang")[1]
  if (!is.na(out0)){
  	out <- "clang"
  } else {
    out0 <- stringr::str_match(out$stdout, "gcc")[1]
    if (!is.na(out0)){
      out <- "gcc"
    } else {
      out <- out$stdout
    }
  }
  out
}

os <- get_os()
cpp_compiler <- gcc_version()
R_compiler <- Rcomp_version()
R_Make_Comp <- Rcpp_version()
print("-------------------------------CONFIG START------------------------")
print(cpp_compiler)
print(R_compiler)
print(R_Make_Comp)
print("-------------------------------CONFIG END------------------------")

if (os=="linux"){
	if (cpp_compiler=="gcc"){
	    if (R_compiler=="gcc"){
		    define(PKG_CXXFLAGS = "PKG_CXXFLAGS=-fopenmp")
		    define(PKG_LIBS = 'PKG_LIBS = `$(R_HOME)/bin/Rscript -e "Rcpp:::LdFlags()"` $(LAPACK_LIBS) $(BLAS_LIBS) $(FLIBS) -fopenmp')
		    define(CPPFLAGS = "#CPPFLAGS += -Xclang -fopenmp")
		    define(LDFLAGS = "#LDFLAGS += -lomp")
		    configure_file("src/Makevars.in")
	    } else {
	        define(PKG_CXXFLAGS = "#PKG_CXXFLAGS=-fopenmp")
		    define(PKG_LIBS = 'PKG_LIBS = `$(R_HOME)/bin/Rscript -e "Rcpp:::LdFlags()"` $(LAPACK_LIBS) $(BLAS_LIBS) $(FLIBS)')
		    define(CPPFLAGS = "#CPPFLAGS += -Xclang -fopenmp")
		    define(LDFLAGS = "#LDFLAGS += -lomp")
		    configure_file("src/Makevars.in")
	    }
	} else if (cpp_compiler=='clang'){
		define(PKG_CXXFLAGS = "#PKG_CXXFLAGS=-fopenmp")
		define(PKG_LIBS = 'PKG_LIBS = `$(R_HOME)/bin/Rscript -e "Rcpp:::LdFlags()"` $(LAPACK_LIBS) $(BLAS_LIBS) $(FLIBS)')
		define(CPPFLAGS = "#CPPFLAGS += -Xclang -fopenmp")
		define(LDFLAGS = "#LDFLAGS += -lomp")
		configure_file("src/Makevars.in")
	}
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
