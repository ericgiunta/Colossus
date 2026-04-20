.onAttach <- function(libname, pkgname) {
  # Checks if the system can use openmp
  syscheck <- Colossus::System_Version() # nocov
  # Starts by checking if openmp is linked with C++
  OpenMP <- syscheck[["OpenMP Enabled"]] # nocov
  if (!OpenMP) { # nocov
    Sys.setenv(ColossusOMP = "FALSE") # nocov
  } else { # nocov
    Sys.setenv(ColossusOMP = "TRUE") # nocov
    Sys.setenv(ColossusGCC = "TRUE") # nocov
    os <- syscheck[["Operating System"]] # nocov
    # Issues have been found on some older clang machines
    if (os == "linux") { # nocov
      cpp_compiler <- syscheck[["Default c++"]] # nocov
      if (cpp_compiler != "") { # nocov
        if (cpp_compiler == "gcc") { # nocov
          R_compiler <- syscheck[["R Compiler"]] # nocov
          if (R_compiler != "gcc") { # nocov
            Sys.setenv(ColossusGCC = "FALSE") # nocov
          } # nocov
        } else if (cpp_compiler == "clang") { # nocov
          Sys.setenv(ColossusGCC = "FALSE") # nocov
        } # nocov
      } else { # nocov
        R_compiler <- syscheck[["R Compiler"]] # nocov
        if (R_compiler != "gcc") { # nocov
          Sys.setenv(ColossusGCC = "FALSE") # nocov
        }
      }
    }
  }
}
