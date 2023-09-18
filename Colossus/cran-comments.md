── R CMD check results ─ Colossus 0.9.2 ─MacOS─
 Duration: 4m 18.7s

❯ checking compilation flags in Makevars ... WARNING
Warning:   Variables overriding user/site settings:
    CPPFLAGS: -o /dev/null -Xclang -fopenmp

❯ checking for GNU extensions in Makefiles ... WARNING
Warning:   Found the following file(s) containing GNU extensions:
    src/Makevars
  Portable Makefiles do not use GNU extensions such as +=, :=, $(shell),
  $(wildcard), ifeq ... endif, .NOTPARALLEL See section ‘Writing portable
  packages’ in the ‘Writing R Extensions’ manual.

0 errors ✔ | 2 warnings ✖ | 0 notes ✔

── R CMD check results ─ Colossus 0.9.2 ─Windows─
Duration: 3m 17.6s

❯ checking compilation flags in Makevars ... WARNING
Warning:   Non-portable flags in variable 'PKG_CXXFLAGS':
    -fopenmp

0 errors ✔ | 1 warning ✖ | 0 notes ✔

 ── R CMD check results ─ Colossus 0.9.2 ─Linux─
Duration: 3m 15.2s

❯ checking compilation flags in Makevars ... WARNING
Warning:   Non-portable flags in variable 'PKG_CXXFLAGS':
    -fopenmp

❯ checking installed package size ... NOTE
    installed size is 62.0Mb
    sub-directories of 1Mb or more:
      libs  61.1Mb

0 errors ✔ | 1 warning ✖ | 1 note ✖
