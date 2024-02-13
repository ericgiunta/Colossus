# Contributing to Colossus

First of all, thanks for considering contributing to Colossus.

## Code of conduct

Please note that this project is released with a Contributor Code of Conduct. By participating in this project you agree to abide by its terms.

## Contributing to Colossus

There are several ways you can contribute to this project:
    -Think Colossus is useful? Let others discover it, by telling them in person, via Twitter or a blog post.
    -Using Colossus for a paper you are writing? Consider citing it.
    -Using Colossus and got stuck? Browse the documentation to see if you can find a solution. Still stuck? Post your question as an issue on GitHub. While we cannot offer user support, we'll try to do our best to address it, as questions often lead to better documentation or the discovery of bugs.
    -Want to ask a question in private? Contact the package maintainer by [egiunta@ksu.edu].
    -Have an idea for a new Colossus feature? Take a look at the documentation and issue list to see if it isn't included or suggested yet. If not, suggest your idea as an issue on GitHub. While we can't promise to implement your idea, it helps to:
        -Explain in detail how it would work.
        -Keep the scope as narrow as possible.
    -See below if you want to contribute code for your idea as well.
    -Using Colossus and discovered a bug? Don't let others have the same experience and report it as an issue on GitHub so we can fix it. A good bug report makes it easier for us to do so, so please include:
        -Your operating system name and version (e.g. Mac OS 10.13.6).
        -Any details about your local setup that might be helpful in troubleshooting.
        -Detailed steps to reproduce the bug.
    -Noticed a typo on the website? Think a function could use a better example? Good documentation makes all the difference, so your help to improve it is very welcome!

We try to follow the GitHub flow for development.
    -Fork this repo and clone it to your computer. To learn more about this process, see this guide.
    -If you have forked and cloned the project before and it has been a while since you worked on it, pull changes from the original repo to your clone by using git pull upstream master.
    -Open the RStudio project file (.Rproj).
    -Make your changes:
        -Write your code.
        -Test your code (bonus points for adding unit tests).
        -Document your code (see function documentation above).
        -Check your code with devtools::check() and aim for 0 errors and warnings.
    -Commit and push your changes.
    -Every Github action is expected to pass
    -Submit a pull request.

## Roadmap
The development of Colossus is split into stages based on relative stages in the Million Person Study.
    -Stage 1: Validated Cox Proportional Hazards and Poisson regressions and plotting, for repeating previous study results and initial study
    -Stage 2: Validated Grey-Fine regression and additional testing capability
    -Stage 3: We expect additional capability to be added as Million Person Study needs evolve
Colossus is currently at stage 2. Modifications to underlying C++ code may occur, but user code with existing functions should not be impacted.
