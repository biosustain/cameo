.. highlight:: shell

============
Contributing
============

Contributions are welcome, and they are greatly appreciated! Every little bit helps, and credit will always be given.

You can contribute in many ways:

Types of Contributions
----------------------

Report Bugs
~~~~~~~~~~~

Report bugs at https://github.com/biosustain/cameo/issues.

If you are reporting a bug, please include:

* Your operating system name and version.
* Any details about your local setup that might be helpful in troubleshooting.
* Detailed steps to reproduce the bug.

Fix Bugs
~~~~~~~~

Look through the GitHub issues for bugs. Anything tagged with "bug" and "help wanted" is open to whoever wants to
implement it.

Implement Features
~~~~~~~~~~~~~~~~~~

Look through the GitHub issues for features. Anything tagged with "enhancement" and "help wanted" is open to whoever
wants to implement it.

Write Documentation
~~~~~~~~~~~~~~~~~~~

cameo could always use more documentation, whether as part of the official cameo using data on iloop docs, in
docstrings, or even on the web in blog posts, articles, and such - all contributions are welcome!

Submit Feedback
~~~~~~~~~~~~~~~

The best way to send feedback is to file an issue at https://github.com/biosustain/cameo/issues.

If you are proposing a feature:

* Explain in detail how it would work.
* Keep the scope as narrow as possible, to make it easier to implement.
* Remember that this is a volunteer-driven project, and that contributions
  are welcome :)

Get Started!
------------

Want to contribute a new feature or improvement? Consider starting by raising an issue and assign it to yourself to
describe what you want to achieve. This way, we reduce the risk of duplicated efforts and you may also get
suggestions on how to best proceed, e.g. there may be half-finished work in some branch that you could start with.

Here's how to set up `cameo` for local development.

1. Fork the `cameo` repo on GitHub.
2. Clone your fork locally::

    git clone git@github.com:your_name_here/cameo.git

3. Install your local copy into a virtualenv. Assuming you have virtualenvwrapper installed, this is how you set up your fork for local development ::

    mkvirtualenv cameo
    cd cameo/
    python setup.py develop

4. Create a branch for local development::

    git checkout -b name-of-your-bugfix-or-feature

   Now you can make your changes locally.

5. When you are done making changes, check that your changes pass flake8 and the tests, including testing other Python versions with tox::

    tox

   or::

    tox -e py35

   To only test on python 3.5 (tox on all python versions takes quite long so may want to try one at a time to start with.
   To get flake8 and tox, just pip install them into your virtualenv.

6. Commit your changes and push your branch to GitHub::

    git add .
    git commit -m "Your detailed description of your changes."
    git push origin name-of-your-bugfix-or-feature

7. Submit a pull request through the GitHub website.

Pull Request Guidelines
-----------------------

Before you submit a pull request, check that it meets these guidelines:

1. The pull request should include tests. Except in very rare circumstances, code coverage must not decrease (as
   reported by codecov which runs automatically when you submit your pull request)
2. If the pull request adds functionality, the docs should be updated. Put
   your new functionality into a function with a docstring and consider creating a notebook that demonstrates the usage in `docs`
3. The pull request should work for Python 2.7, 3.4 and 3.5. Check
   https://travis-ci.org/biosustain/cameo/pull_requests
   and make sure that the tests pass for all supported Python versions.
4. Assign a reviewer to your pull request. If in doubt, assign Niko Sonnenschein. Your pull request must be
   approved by all reviewers before it can be merged.
