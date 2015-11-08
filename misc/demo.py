#!/usr/bin/env python
from __future__ import division, print_function
from sys import argv
from collections import defaultdict, namedtuple
from json import dumps
from lpi import BlastFile, NCBITaxonomyDB, ALL_RANKS

db = NCBITaxonomyDB()

Lineage = namedtuple("Lineage", ['query', 'gi', 'genus', 'phylum', 'bitscore'])

query = None
lineages = []
badgi = 0
for hit in BlastFile(argv[1]):
    qry = hit.qseqid
    if query is None or query.get('name', '') != qry:
        if query is not None:
            #max_bitscore = max(lineages, key=lambda x: x.bitscore).bitscore
            genus_table = dict()
            for lineage in lineages:
                score = lineage.bitscore
                genus_id = lineage.genus['taxid']
                if genus_id not in genus_table or \
                        score > genus_table[genus_id].bitscore:
                    genus_table[genus_id] = lineage
            phyla = {}
            for genus in genus_table.values():
                phylum = genus.phylum['name']
                if phylum not in phyla:
                    phyla[phylum] = {'num': 0, 'score': 0.0}
                phyla[phylum]['score'] += genus.bitscore
                phyla[phylum]['num'] += 1

            for phylum in phyla:
                score = phyla[phylum]['score']
                num = phyla[phylum]['num']
                query['taxa'][phylum] = score # / num
            query['top_taxon'] = max(query['taxa'])

            #print(dumps(query, sort_keys=True))
            print(query['top_taxon'])

        query = dict(name=qry, taxa=defaultdict(float), num_unassigned=0,
                     num_assigned=0)

    gi = hit.sseqid.gi
    try:
        tax = db.sequence_lineage(gi, ALL_RANKS)
        if len(tax) < 3:
            raise KeyError("Insufficent levels of tax tree found")
        query['num_assigned'] += 1
    except KeyError:
        query['num_unassigned'] += 1
        continue
    try:
        gen = tax['genus']
    except KeyError:
        gen = list(tax.values())[1]

    try:
        phy = tax['phylum']
    except KeyError:
        gen = list(tax.values())[-2]

    lineages.append(Lineage(qry, gi, gen, phy, hit.bitscore))
