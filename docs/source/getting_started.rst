Getting Started
===============
COSMO is intended to be installed and its monitors run either manually or via cronjob.
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

Configuration
--------------
COSMO Settings with a configuration file
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
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

``monitorframe`` also requires a ``yaml`` configuration file with the following:

.. code-block:: yaml

    # Monitor data database
    data:
      db_settings:
        database: ''
        pragmas:
          journal_mode: 'wal'
          foreign_keys: 1
          ignore_check_constraints: 0
          synchronous: 0

    # Monitor status and results database
    results:
      db_settings:
        database: ''
        pragmas:
          journal_mode: 'wal'
          foreign_keys: 1
          ignore_check_constraints: 0
          synchronous: 0

This configuration file should be set to an environment variable called ``MONITOR_CONFIG``.

You can store these configurations in the same file and have both of the environment variables point to the same file.

.. warning::

    Use proper precautions around your configuration file.
    It may or may not contain sensitive information, so please ensure that permissions on that file are restricted to
    the intended users.
    DON'T push it to GitHub!

CRDS
^^^^
Some of the COSMO DataModels utilize data from reference files, and take advantage of ``crds`` to do so.
For configuration and setup instructions for using ``crds``, see
`the crds user manual <https://hst-crds.stsci.edu/static/users_guide/environment.html>`_.

At minimum, users will need access to a CRDS cache with the following reference file types:

- LAMPTAB
- WCPTAB

Since the COSMO monitors use data from reference files across time, it would be best to get all files of those types
available in the *active context*.

The easiest way to ensure that the local CRDS cache has everything required, users can use::

    crds sync --contexts hst-cos-operational --fetch-references

This command with download *all* COS reference files and mappings to the ``CRDS_CACHE`` (see the instructions mentioned
above).

.. warning::

    The command given above works well, but there's a caveat: it requires a large amount of available storage space at
    the cache location.

Running Tests
-------------
COSMO includes a suite of tests for the package.
For developers, it's a good idea to execute these tests whenever there are changes to the code or environment.

If you're in the project directory, you can execute the tests with::

    python -m pytest

For executing the tests with coverage (after ``coverage`` has been installed), use::

    coverage run -m pytest

Executing Monitors
------------------
Monitors can be executed by using the monitoring classes directly:

.. code-block:: python

    from cosmo.monitors import AcqImageMonitor

    monitor = AcqImageMonitor()

    # Run it
    monitor.monitor()

Or, they can be executed from the command line::

    (cosmoenv) mycomputer:~ user$ cosmo --monthly

For more command line options::

    (cosmoenv) mycomputer:~ user$ cosmo --help

