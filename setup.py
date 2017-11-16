#!/usr/bin/env python
"""
This is a setup script to install pyctscripts
"""

import setuptools

if __name__ == "__main__":
    setuptools.setup(
        name='scripts',
        version="0.1.1",
        description='Collection of programming scripts to automate variety of \
                     repetitive research tasks',
        author='Viswanath Pasumarthi',
        author_email='pasumart@buffalo.edu',
        url="https://github.com/vpasumarthi/scripts",
        packages=setuptools.find_packages(),
        install_requires=[
            'numpy>=1.7',
            'pyyaml'
        ],
        extras_require={
            'docs': [
                'sphinx==1.2.3',  # autodoc was broken in 1.3.1
                'sphinxcontrib-napoleon',
                'sphinx_rtd_theme',
                'numpydoc',
            ],
            'tests': [
                'pytest',
                'pytest-cov',
                'pytest-pep8',
                'tox',
            ],
        },

        tests_require=[
            'pytest',
            'pytest-cov',
            'pytest-pep8',
            'tox',
        ],

        classifiers=[
            'Development Status :: 4 - Beta',
            'Intended Audience :: Science/Research',
            'Programming Language :: Python :: 3',
        ],
        zip_safe=True)
