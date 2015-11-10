"""
Parse blast tabular formats
"""

from __future__ import division, absolute_import, print_function
from collections import namedtuple
import pandas as pd

class SeqID(object):
    _id_types = []

    def __init__(self, seqid_str):
        # Remove whitespace and trailing |s
        seqid_str = seqid_str.strip().rstrip('|')
        # Split into a list of id types and values
        split = list(map(lambda x: x.strip(), seqid_str.split('|')))
        # Ensure that we have a set of pairs of id_type, value
        if len(split) % 2 != 0:
            raise ValueError("Odd number of 'id|value' pairs", seqid_str)
        # Form the list into a dict.
        pairs = {split[i]: split[i+1] for i in range(0, len(split), 2)}
        for id_type, value in pairs.items():
            # Parse integral IDs to ints.
            try:
                value = int(value)
            except ValueError:
                pass
            setattr(self, id_type, value)
            self._id_types.append(id_type)

    def __eq__(self, other):
        if self._id_types != other._id_types:
            return False
        for id_type in self._id_types:
            if getattr(self, id_type) != getattr(other, id_type):
                return False
        return True


BLAST_FIELDS = {
    'qseqid': str,
    'qgi': int,
    'qacc': float,
    'qaccver': float,
    'qlen': int,
    'sseqid': SeqID,
    'sallseqid': SeqID,
    'sgi': float,
    'sallgi': float,
    'sacc': float,
    'saccver': float,
    'sallacc': float,
    'slen': float,
    'qstart': int,
    'qend': int,
    'sstart': int,
    'send': int,
    'qseq': float,
    'sseq': float,
    'evalue': float,
    'bitscore': float,
    'score': float,
    'length': int,
    'pident': float,
    'nident': float,
    'mismatch': float,
    'positive': float,
    'gapopen': int,
    'gaps': float,
    'ppos': float,
    'frames': float,
    'qframe': float,
    'sframe': float,
    'btop': float,
    'staxids': float,
    'sscinames': float,
    'scomnames': float,
    'sblastnames': float,
    'sskingdoms': float,
    'stitle': float,
    'salltitles': float,
    'sstrand': float,
    'qcovs': float,
    'qcovhsp': float,
}


DEFAULT_BLAST_FIELDS = [
    'qseqid',
    'sseqid',
    'pident',
    'length',
    'mismatch',
    'gapopen',
    'qstart',
    'qend',
    'sstart',
    'send',
    'evalue',
    'bitscore',
]


class BlastFile(object):
    '''Returns a BlastRecord for each record in ``filename``'''

    fh = None

    def __init__(self, filename, fields=DEFAULT_BLAST_FIELDS):
        self.record = namedtuple("BlastRecord", fields)
        self.filename = filename
        self.fh = open(self.filename)
        self.fields = fields

    def __iter__(self):
        if self.fh:
            self.fh.seek(0)
        return self

    def __next__(self):
        line = self.fh.readline().strip()
        if not line:
            raise StopIteration()
        values = line.split('\t')
        # Convert each value given the blast field dictionary.
        for i, field in enumerate(self.fields):
            values[i] = BLAST_FIELDS[field](values[i])
        return self.record(*values)

    def __enter__(self):
        pass

    def __exit__(self):
        if self.fh:
            self.fh.close()


def parse_blast_queries(filename, fields=DEFAULT_BLAST_FIELDS, filter='',
                        chunksize=5000):
    '''Iterate over a blast table, parsing each query's hits into a separate
    pandas data frame. Assumes the table is sorted by query sequence ID.

    Field names must be specified if a format other than the standard '-outfmt
    6' tabular format is being parsed. An optional filter on table columns may
    be provided, e.g 'evalue < 1e-10'. See ``pandas.DataFrame.query()`` for
    more details on filtering.
    '''

    # Algorithm:
    #  - keep a dictonary of unfinished queries.
    #  - for each chunk, for all unfinished queries, check if they don't appear
    #    in the current chunk, if they don't, yield them.
    #  - otherwise, append or add the current queries to the dict of unfinished
    #    queries.
    #  - finally, return all query: dataframe pairs.

    lastgroups = {}
    for chunk in pd.read_table(filename, names=fields, chunksize=chunksize):
        if filter:
            chunk = chunk.query(filter)
        groups = {k:df for k, df in chunk.groupby(['qseqid',])}
        # To avoid changing dict len when popping in a for loop, we keep track
        # of things we have yeild-ed and therefore want to pop outside of the
        # loop (in the next one).
        to_pop = []
        for key in lastgroups:
            if key in groups:
                lastgroups[key] = pd.concat((lastgroups[key], df))
            else:
                to_pop.append(key)
                yield key, lastgroups[key]
        for key in to_pop:
            lastgroups.pop(key)
        for key, df in groups.items():
            lastgroups[key] = df
    for key, df in lastgroups.items():
        yield key, df
