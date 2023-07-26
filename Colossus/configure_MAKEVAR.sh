#!/bin/bash

#make the Makevars file
if [ ! -e "./src/Makevars" ]; then
touch ./src/Makevars
fi


case "$OSTYPE" in
  linux*)   echo 'PKG_CXXFLAGS=-fopenmp
PKG_LIBS = `$(R_HOME)/bin/Rscript -e "Rcpp:::LdFlags()"` $(LAPACK_LIBS) $(BLAS_LIBS) $(FLIBS) -fopenmp' > ./src/Makevars ;;
  darwin*)  echo 'PKG_CPPFLAGS='-Xclang -fopenmp'
LDFLAGS="-L/usr/local/opt/llvm/lib/c++ -Wl,-rpath,/usr/local/opt/llvm/lib/c++"
LDFLAGS += -L/opt/homebrew/opt/libomp/lib -lomp
CPPFLAGS += -I/opt/homebrew/opt/libomp/include -Xclang -fopenmp' > ./src/Makevars ;; 
  msys*)    echo "windows" ;;
  solaris*) echo "solaris" ;;
  bsd*)     echo "bsd" ;;
  *)        echo "unknown" ;;
esac

