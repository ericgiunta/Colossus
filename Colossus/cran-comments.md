── R CMD check results ─ Colossus 0.9.2 ─MacOS─
Duration: 3m 22s

❯ checking compilation flags in Makevars ... WARNING
Warning:   Variables overriding user/site settings:
    CPPFLAGS: -o /dev/null -Xclang -fopenmp

❯ checking for GNU extensions in Makefiles ... WARNING
Warning:   Found the following file(s) containing GNU extensions:
    src/Makevars
  Portable Makefiles do not use GNU extensions such as +=, :=, $(shell),
  $(wildcard), ifeq ... endif, .NOTPARALLEL See section ‘Writing portable
  packages’ in the ‘Writing R Extensions’ manual.

❯ checking top-level files ... NOTE
  Non-standard file/directory found at top level:
    ‘configure_MAKEVAR.sh’

0 errors ✔ | 2 warnings ✖ | 1 note ✖

── R CMD check results ─ Colossus 0.9.2 ─Windows─
Duration: 3m 32.3s

❯ checking compilation flags in Makevars ... WARNING
Warning:   Non-portable flags in variable 'PKG_CXXFLAGS':
    -fopenmp

❯ checking top-level files ... NOTE
  Non-standard file/directory found at top level:
    'configure_MAKEVAR.sh'

0 errors ✔ | 1 warning ✖ | 1 note ✖

 ── R CMD check results ─ Colossus 0.9.2 ─Linux─
Duration: 2m 29.6s

❯ checking compilation flags in Makevars ... WARNING
Warning:   Non-portable flags in variable 'PKG_CXXFLAGS':
    -fopenmp

❯ checking installed package size ... NOTE
    installed size is 62.0Mb
    sub-directories of 1Mb or more:
      libs  61.1Mb

❯ checking top-level files ... NOTE
  Non-standard file/directory found at top level:
    ‘configure_MAKEVAR.sh’

0 errors ✔ | 1 warning ✖ | 2 notes ✖
