.onAttach <- function(libname, pkgname) {
  syscheck <- Colossus::System_Version()
  OpenMP <- syscheck[["OpenMP Enabled"]]
  if (!OpenMP) {
    Sys.setenv(ColossusOMP = "FALSE")
  } else {
    Sys.setenv(ColossusOMP = "TRUE")
    Sys.setenv(ColossusGCC = "TRUE")
    os <- syscheck[["Operating System"]]
    if (os == "linux") {
      cpp_compiler <- syscheck[["Default c++"]]
      if (cpp_compiler == "gcc") {
        R_compiler <- syscheck[["R Compiler"]]
        if (R_compiler != "gcc") { # nocov
          Sys.setenv(ColossusGCC = "FALSE")
        }
      } else if (cpp_compiler == "clang") { # nocov
        Sys.setenv(ColossusGCC = "FALSE")
      }
    }
  }
}
