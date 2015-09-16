import hashlib
import logging
import os
from os import path
import sys

import requests
from six.moves.configparser import ConfigParser


LOG = None
LOGLEVEL = logging.DEBUG


def get_config():
    '''Returns a dict of config options from the config file'''
    cfgparse = ConfigParser()
    cfgparse.read(['.lpipyrc', path.expanduser('~/.lpipyrc')])
    config = cfgparse.defaults()
    local = os.environ.get('VIRTUAL_ENV', path.expanduser('~/.local/'))
    datadir = config.get('path', path.join(local, 'var', 'lpi'))

    return {
        "datadir": datadir,
    }


def get_data_dir():
    '''
    Find (and create in requrired) a lpi-specific data directory under either
    the current virtualenv directory, falling back on '~/.local' if we're not
    in a virtualenv.
    '''
    datadir = path.expanduser(get_config()["datadir"])
    if not path.isdir(datadir):
        os.makedirs(datadir, 0o755)
    return datadir


def get_data_file(filename):
    '''Get a data file (or subdir) contained within the current lpi data dir'''
    datadir = path.join(get_data_dir(), filename)
    return datadir


def init_logger(name='lpi', quiet=False, level=LOGLEVEL):
    log = logging.getLogger(name)
    log.setLevel(level)

    if not quiet:
        stdout = logging.StreamHandler(stream=sys.stderr)
        stdout.setLevel(logging.DEBUG)
        log.addHandler(stdout)

    return log


def get_logger(quiet=False, level=LOGLEVEL):
    global LOG
    if not LOG:
        LOG = init_logger('lpi', quiet, level)
    return LOG


def download_file(url, filename, chunk_size=128 * 1024):
    log = get_logger()
    req = requests.get(url, stream=True)
    cnt = 0
    log.info("Downloading %s to %s", url, filename)
    with open(filename, 'wb') as fh:
        for chunk in req.iter_content(chunk_size):
            cnt += len(chunk)
            log.debug("{}MiB Downloaded".format(cnt / (1024 ** 2)))
            fh.write(chunk)
    log.debug("Downloaded {}".format(url))


def read_remote_file(url):
    req = requests.get(url)
    return req.text


def md5sum(filename):
    h = hashlib.md5()
    with open(filename, 'rb') as fh:
        while True:
            chunk = fh.read(1048576)
            h.update(chunk)
            if not chunk:
                break
    return h.hexdigest()
