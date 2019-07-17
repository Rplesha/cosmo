Getting Started
===============
COSMO is intended to installed and it's monitors run either manually or via cronjob.
Use of COSMO outside of the scope of executing the monitors that are defined is not recommended, however, the data
models and their API can be used to access the monitoring data for use outside the scope of this project if needed.

Installing
----------
Before installing COSMO, be sure to create an environment with python 3.7+ at the minimum.
A good starting point would be::

    conda create -n cosmo_env python=3.7 stsci

For developers, also including ``coverage`` is also recommended (but not mandatory).

After an environment has been prepared, clone the repository::

    git clone https://github.com/spacetelescope/cosmo.git

Then install using pip::

    cd cosmo
    pip install .


Settings with a configuration file
----------------------------------
To manage configurations, COSMO uses a ``yaml`` configuration file.
Create a yaml configuration file with the following format:

.. code-block:: yaml

    # Settings for obtaining data from files
    filesystem:
      source: ''  # Path to the source files

    # Settings for obtaining data from sms files and setting up the sms database
    sms:
      source: ''  # This is the path where the sms files exist
      db_settings:  # These are database keyword arguments
        database: ''  # Path to sqlite database file
        pragmas:  # sqlite database connection configurations.
          journal_mode: 'wal'
          foreign_keys: 1
          ignore_check_constraints: 0
          synchronous: 0

    output: ''

For more information on sqlite pragma statements, see `this <https://www.sqlite.org/pragma.html>`_.

Once the file is ready, set it as an environment variable, ``COSMO_CONFIG``.

.. warning::

    Use proper precautions around your configuration file.
    It may or may not contain sensitive information, so please ensure that permissions on that file are restricted to
    the intended users.
    DON'T push it to GitHub!

Running tests
-------------
COSMO includes a suite of tests for the package.
For developers, it's a good idea to execute these tests whenever there are changes to the code or environment.

If you're in the project directory, you can execute the tests with::

    python -m pytest

For executing the tests with coverage (after ``coverage`` has been installed), use::

    coverage run -m pytest
