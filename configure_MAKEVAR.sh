#!/bin/bash

#make the Makevars file
if [ ! -e "./src/Makevars" ]; then
touch ./src/Makevars
fi


case "$OSTYPE" in
  linux*)   echo 'PKG_CXXFLAGS=-fopenmp
PKG_LIBS = `$(R_HOME)/bin/Rscript -e "Rcpp:::LdFlags()"` $(LAPACK_LIBS) $(BLAS_LIBS) $(FLIBS) -fopenmp' > ./src/Makevars ;;
  darwin*)  echo 'CPPFLAGS += -Xclang -fopenmp
LDFLAGS += -lomp' > ./src/Makevars ;; 
  msys*)    echo "windows" ;;
  solaris*) echo "solaris" ;;
  bsd*)     echo "bsd" ;;
  *)        echo "unknown" ;;
esac

