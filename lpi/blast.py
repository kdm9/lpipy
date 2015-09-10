"""
Parse blast tabular formats
"""

from __future__ import division, absolute_import, print_function
from collections import namedtuple

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
