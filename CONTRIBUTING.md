Contributing
============

Bug reports, feature suggestions and other contributions are greatly
appreciated!  pysat and pysatCDAAC are community-driven projects that
welcome both feedback and contributions.

Short version
-------------

* Submit bug reports, feature requests, and questions at
`GitHub Issues <https://github.com/pysat/pysatCDAAC/issues>`_
* Make pull requests to the ``develop`` branch

More about Issues
-----------------

Bug reports, questions, and feature requests should all be made as GitHub
Issues.  Templates are provided for each type of issue, to help you include
all the necessary information.

Questions
^^^^^^^^^

Not sure how something works?  Ask away!  The more information you provide, the
easier the question will be to answer.  You can also interact with the pysat
developers on our `slack channel <https://pysat.slack.com>`_.  

Bug reports
^^^^^^^^^^^

When reporting a bug please include:

* Your operating system name and version
* Any details about your local setup that might be helpful in troubleshooting
* Detailed steps to reproduce the bug

Feature requests
^^^^^^^^^^^^^^^^

If you are proposing a new feature or a change in something that already exists:

* Explain in detail how it would work.
* Keep the scope as narrow as possible, to make it easier to implement.
* Remember that this is a volunteer-driven project, and that code contributions
  are welcome :)

More about Development
----------------------

To set up `pysatCDAAC` for local development:

1. Fork pysatCDAAC on
   `GitHub <https://github.com/pysat/pysatCDAAC/fork>`_.
2. Clone your fork locally::

    git clone git@github.com:your_name_here/pysatCDAAC.git

3. Create a branch for local development::

    git checkout -b name-of-your-bugfix-or-feature

4. Make your changes locally. Tests for new instruments are performed
   automatically.  Tests for custom functions should be added to the
   appropriately named file in ``pysatCDAAC/tests``.  For example,
   methods contained in ``pysatCDAAC/instruments/methods/usercode.py``
   should be named ``pysatCDAAC/tests/test_methods_usercode.py``.  If no test
   file exists, then you should create one.  This testing uses pytest, which
   will run tests on any python file in the test directory that starts with
   ``test``.  Test classes must begin with ``Test``, and test methods must also
   begin with ``test``.

5. When you're done making changes, run all the checks to ensure that nothing
   is broken on your local system::

    pytest -vs pysatCDAAC

6. Update/add documentation (in ``docs``).  Even if you don't think it's
   relevant, check to see if any existing examples have changed.

7. Add your name to the .zenodo.json file as an author

8. Commit your changes and push your branch to GitHub::

    git add .
    git commit -m "Brief description of your changes"
    git push origin name-of-your-bugfix-or-feature

9. Submit a pull request through the GitHub website. Pull requests should be
   made to the ``develop`` branch.

Pull Request Guidelines
^^^^^^^^^^^^^^^^^^^^^^^

If you need some code review or feedback while you're developing the code, just
make a pull request. Pull requests should be made to the ``develop`` branch.

For merging, you should:

1. Include an example for use
2. Add a note to ``CHANGELOG.md`` about the changes
3. Ensure that all checks passed (current checks include Travis-CI
   and Coveralls) [1]_

.. [1] If you don't have all the necessary Python versions available locally or
       have trouble building all the testing environments, you can rely on
       Travis to run the tests for each change you add in the pull request.
       Because testing here will delay tests by other developers, please ensure
       that the code passes all tests on your local system first.

Project Style Guidelines
^^^^^^^^^^^^^^^^^^^^^^^^

In general, pysat follows PEP8 and numpydoc guidelines.  Pytest runs the unit
and integration tests, flake8 checks for style, and sphinx-build performs
documentation tests.  However, there are certain additional style elements that
have been settled on to ensure the project maintains a consistent coding style.
These include:

* Line breaks should occur before a binary operator (ignoring flake8 W503)
* Combine long strings using `join`
* Preferably break long lines on open parentheses rather than using `\`
* Use no more than 80 characters per line
* Avoid using Instrument class key attribute names as unrelated variable names:
  `platform`, `name`, `tag`, and `inst_id`
* The pysat logger is imported into each sub-module and provides status updates
  at the info and warning levels (as appropriate)
* Several dependent packages have common nicknames, including:
  * `import datetime as dt`
  * `import numpy as np`
  * `import pandas as pds`
  * `import xarray as xr`
* All classes should have `__repr__` and `__str__` functions
* Docstrings use `Note` instead of `Notes`
* Try to avoid creating a try/except statement where except passes
* Use setup and teardown in test classes
* Use pytest parametrize in test classes when appropriate
* Provide testing class methods with informative failure statements and
  descriptive, one-line docstrings
* Block and inline comments should use proper English grammar and punctuation
  with the exception of single sentences in a block, which may then omit the
  final period
* When casting is necessary, use `np.int64` and `np.float64` to ensure operating 
   system agnosticism
