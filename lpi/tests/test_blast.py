"""
Test of lpi.blast module
"""

from __future__ import division, absolute_import, print_function
from nose.tools import assert_raises
from nose import tools as nt
import json
from lpi import blast
from . import helpers


def test_seqid_type():
    '''Test that the SeqID type behaves as expected.'''
    def _do_seqid(seqid, expected):
        res = blast.SeqID(seqid)
        for k, v in expected.items():
            assert(getattr(res, k) == v), (k, v)

    _do_seqid('gi|12345', {'gi': 12345})
    _do_seqid('gi|12345\n', {'gi': 12345})
    _do_seqid('gi|12345|', {'gi': 12345})

    _do_seqid('gi|12345|ref|WP_233421.2',
              {'gi': 12345, 'ref': 'WP_233421.2'})
    _do_seqid('gi|12345|ref|WP_233421.2\t',
              {'gi': 12345, 'ref': 'WP_233421.2'})

    with assert_raises(ValueError):
        _do_seqid('gi|12345|ref', {'gi': 12345})



def test_parse_blast_default():
    '''Test that the default tabular format columns are parsed correctly.'''
    def _do_bf(fname, expected):
        bf = blast.BlastFile(fname)
        for res, expt in zip(bf, expected):
            for k, v in expt.items():
                assert(getattr(res, k) == v), (k, v, getattr(res, k))

    _do_bf(helpers.get_data_file('default_blastx_table.tab'), [
           {'qseqid': 'seq1',
            'sseqid': blast.SeqID('gi|496001848|ref|WP_008726427.1|'),
            'pident': 48.57,
            'length': 70,
            'mismatch': 34,
            'gapopen': 1,
            'qstart': 25,
            'qend': 228,
            'sstart': 31,
            'send': 100,
            'evalue': 1e-10,
            'bitscore': 61.6,},
           {'qseqid': 'seq2',
            'sseqid': blast.SeqID('gi|489566164|ref|WP_003470669.1|'),
            'pident': 100.0,
            'length': 60,
            'mismatch': 0,
            'gapopen': 0,
            'qstart': 267,
            'qend': 88,
            'sstart': 216,
            'send': 275,
            'evalue': 2e-34,
            'bitscore': 125,},
           ])


def __do_test_parse_blast_queries(parser, json_sseqid_file):
    # ensure the order of the parsed queries is consistent
    parser = sorted((q, df) for q, df in parser)
    with open(helpers.get_data_file(json_sseqid_file)) as jsonfh:
        # Loads a dict of {'qseqid': ['sseqid1', 'sseqid2', ...]}
        expect = sorted((q, sseqs) for q, sseqs in json.load(jsonfh).items())

    for (query, df), (expt_query, expt_sseqids) in zip(parser, expect):
        assert query == expt_query
        assert all(df.qseqid == expt_query)
        assert all(df.sseqid == expt_sseqids)


def test_parse_blast_queries():
    '''Check that hits are grouped correctly into query-wise `pd.DataFrame`s.
    '''
    __do_test_parse_blast_queries(
        blast.parse_blast_groupby(helpers.get_data_file('50_blast_hits.tab')),
        helpers.get_data_file('50_blast_hits_sseqids.json')
    )


def test_parse_blast_queries_filtering():
    __do_test_parse_blast_queries(
        blast.parse_blast_groupby(helpers.get_data_file('50_blast_hits.tab'),
                                  filter='evalue < 1e-150'),
        helpers.get_data_file('50_blast_hits_sseqids_1e-150.json')
    )
