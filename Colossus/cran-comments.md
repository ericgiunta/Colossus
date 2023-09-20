── R CMD check results ─ Colossus 0.9.3 ─MacOS─
 Duration: 4m 29.9s

❯ checking for GNU extensions in Makefiles ... WARNING
Warning:   Found the following file(s) containing GNU extensions:
    src/Makevars
  Portable Makefiles do not use GNU extensions such as +=, :=, $(shell),
  $(wildcard), ifeq ... endif, .NOTPARALLEL See section ‘Writing portable
  packages’ in the ‘Writing R Extensions’ manual.

0 errors ✔ | 1 warning ✖ | 0 notes ✔

── R CMD check results ─ Colossus 0.9.3 ─Windows─
 Duration: 3m 40.8s

❯ checking line endings in Makefiles ... WARNING
Warning:   Found the following Makefile(s) with CR or CRLF line endings:
    src/Makevars
  Some Unix 'make' programs require LF line endings.

0 errors ✔ | 1 warning ✖ | 0 notes ✔

 ── R CMD check results ─ Colossus 0.9.3 ─Linux─
 Duration: 3m 34.6s

❯ checking installed package size ... NOTE
    installed size is 62.0Mb
    sub-directories of 1Mb or more:
      libs  61.1Mb

0 errors ✔ | 0 warnings ✔ | 1 note ✖
