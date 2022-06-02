# Contributing to the SCEC Broadband Platform (BBP) System

This document provides an overview on how to contribute to BBP, and will provide step-by-step instructions on a practical contribution workflow.

## Getting Started

* Make sure you have an active GitHub account
* Download and install git
* Read the git documentation
* Install the main branch of the SCECcode/bbp.git repository
* If you haven't worked with Git Forks before, make sure to read the documentation linked below.

## Submitting a Pull Request

Pull requests are great! Please submit them to us! Here's how:

1. Fork the repo. [Some helping info](https://help.github.com/en/github/collaborating-with-issues-and-pull-requests/working-with-forks)
2. Make a new branch. For features/additions base your new branch at `dev`.
3. Add a test! Only pull requests for documentation and refactoring do not require a test.
4. Make sure the tests pass. Run `./run_tests.sh` in the top-level directory of the repo.
5. Push your changes to your fork and submit a pull request. Make sure to set the branch to `ucvm:dev`.
6. Wait for our review. There may be some suggested changes or improvements. Changes can be made after
the pull request has been opening by simply adding more commits to your branch.

Pull requests can be changed after they are opened, so create a pull request as early as possible.
This allows us to provide feedback during development and to answer any questions.

Please make sure to set the correct branch for your pull request. Also, please do not include large files in your pull request.
If you feel that you need to add large files, let us know and we can figure something out.

## Submitting an Issue

Please open an issue if you want to ask a question about BBP.

* Please search through the past issues to see if your question or the bug has already been addressed
* Please apply the correct tag to your issue so others can search

If you want to submit a bug report, please provide the information below:
* BBP version, Python version, and Platform (Linux, Windows, Mac OSX, etc)
* How did you install BBP (Docker, from source...)
* Please provide a short, complete, and correct example that demonstrates the issue.
* If this broke in a recent update, please tell us when it used to work.

## Additional Resources
* [Working with Git Forks](https://help.github.com/en/github/collaborating-with-issues-and-pull-requests/working-with-forks)
* [Style Guide](http://google.github.io/styleguide/pyguide.html)
* [Docs or it doesn’t exist](https://lukeplant.me.uk/blog/posts/docs-or-it-doesnt-exist/)
* Performance Tips:
  * [Python](https://wiki.python.org/moin/PythonSpeed/PerformanceTips)
  * [NumPy and ctypes](https://scipy-cookbook.readthedocs.io/)
  * [SciPy](https://www.scipy.org/docs.html)
  * [NumPy Book](http://csc.ucdavis.edu/~chaos/courses/nlp/Software/NumPyBook.pdf)
