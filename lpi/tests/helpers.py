"""Misc helpers for LPI unit tests.
"""
from __future__ import print_function, absolute_import, division
import atexit
from pkg_resources import resource_filename, Requirement, ResolutionError
from os import path
import shutil

def get_data_file(filename):
    filepath = None
    try:
        filepath = resource_filename(Requirement.parse("lpipy"),
                                     "lpipy/lpi/tests/data/" + filename)
    except ResolutionError:
        pass
    if not filepath or not path.isfile(filepath):
        filepath = path.join(path.dirname(__file__), 'data', filename)
    return filepath


TEMP_DIRECTORIES = []

def get_temp_filename(filename, tempdir=None):
    if tempdir is None:
        tempdir = tempfile.mkdtemp(prefix='test-')
        TEMP_DIRECTORIES.append(tempdir)

    return path.join(tempdir, filename)


@atexit.register
def cleanup():
    global TEMP_DIRECTORIES
    while TEMP_DIRECTORIES:
        path = TEMP_DIRECTORIES.pop()
        shutil.rmtree(path, ignore_errors=True)


