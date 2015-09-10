#!/usr/bin/env python
from setuptools import setup
from setuptools.command.test import test as TestCommand
import versioneer


# Inspired by the example at https://pytest.org/latest/goodpractises.html
class NoseCommand(TestCommand):
    def finalize_options(self):
        TestCommand.finalize_options(self)
        self.test_args = []
        self.test_suite = True

    def run_tests(self):
        # Run nose ensuring that argv simulates running nosetests directly
        import nose
        nose.run_exit(argv=['nosetests'])


desc = """
LPI.py: Lineage Probablity Index calculation
"""

setup_requires = [
    'nose',
    'coverage',
]

install_requires = [
    'six',
    'docopt',
]

test_requires = [
    'pep8',
    'pylint',
]

command_classes=versioneer.get_cmdclass()
command_classes['test'] =  NoseCommand

setup(
    name="lpi",
    packages=['lpi', 'lpi.tests', ],
    scripts=[
    ],
    version=versioneer.get_version(),
    cmdclass=command_classes,
    install_requires=install_requires,
    tests_require=test_requires,
    setup_requires=setup_requires,
    description=desc,
    author="Kevin Murray",
    author_email="spam@kdmurray.id.au",
    url="https://github.com/kdmurray91/lpipy",
    keywords=[
        "bioinformatics",
        "blast",
        "phylogenetics",
        "Next-gen Sequencing",
    ],
    classifiers=[
        "Programming Language :: Python",
        "Programming Language :: Python :: 2",
        "Programming Language :: Python :: 3",
        "Development Status :: 3 - Alpha",
        "Environment :: Console",
        "Intended Audience :: Developers",
        "Intended Audience :: Science/Research",
        "Operating System :: OS Independent",
        "Topic :: Software Development :: Libraries :: Python Modules",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
        ("License :: OSI Approved :: GNU Lesser General Public License v3 or "
         "later (LGPLv3+)"),
    ],
    test_suite="lpi.test",
)
