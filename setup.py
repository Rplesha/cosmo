from setuptools import setup, find_packages

setup(
    name = 'cos_monitoring',
    version = '0.0.1',
    description = 'Provide utilities and monotiring of cos data',
    author = 'Justin Ely',
    author_email = 'ely@stsci.edu',
    keywords = ['astronomy'],
    classifiers = ['Programming Language :: Python',
                   'Programming Language :: Python :: 3',
                   'Development Status :: 1 - Planning',
                   'Intended Audience :: Science/Research',
                   'Topic :: Scientific/Engineering :: Astronomy',
                   'Topic :: Scientific/Engineering :: Physics',
                   'Topic :: Software Development :: Libraries :: Python Modules'],
    packages = find_packages(),
    requires = ['numpy', 'scipy', 'astropy', 'matplotlib'],
    entry_points = {'console_scripts': ['clean_slate=cos_monitoring.database:clean_slate',
                                        'cm_ingest=cos_monitoring.database:do_all',
                                        'cm_monitors=cos_monitoring.database:run_all_monitors',
                                        'create_master_csv=scripts.create_master_csv:main',
                                        'find_new_cos_data=cos_monitoring.retrieval.find_new_cos_data:compare_tables',
					                    'cm_reports=cos_monitoring.database.report:query_all',
                                        'cm_delete=cos_monitoring.database.database:cm_delete',
                                        'cm_describe=cos_monitoring.database.database:cm_describe',
                                        'use_glue=cos_monitoring.database.glue_query:main',
                                        ],
    },
    install_requires = ['setuptools',
                        'numpy',
                        'astropy>=1.0.1',
                        'sqlalchemy>=1.0.12',
                        'pymysql',
                        'matplotlib',
                        'scipy',
                        'fitsio']
    )
