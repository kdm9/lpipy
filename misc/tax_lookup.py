#!/usr/bin/env python3
from collections import OrderedDict
import sqlite3
import sys

CON = sqlite3.connect("lpi-new.db")

def get_merged_id(taxid):
    db = CON.cursor()
    qry = db.execute('SELECT newid FROM merges WHERE oldid = ?', (taxid, ))
    newid = qry.fetchone()
    if newid:
        return newid[0]
    return taxid

def get_taxid_parents(taxid):
    db = CON.cursor()
    qry = db.execute("""SELECT parents.parent, taxa.rank, taxa.name
                        FROM parents INNER JOIN taxa ON taxa.taxid = parents.parent
                        WHERE parents.taxid = {}
                        ORDER BY parents.level;
                     """.format(taxid))
    for parent, rank, name in qry:
        yield (parent, rank, name)

def get_gi_taxid(gi):
    db = CON.cursor()
    qry = db.execute("SELECT taxid FROM sequences WHERE gi = ?", (gi,))
    return qry.fetchone()[0]

DEFAULT_RANKS = set(['species', 'genus', 'family', 'order', 'class', 'phylum',
                     'kingdom', 'superkingdom'])

def taxtree(taxid, ranks=DEFAULT_RANKS):
    taxid = get_merged_id(taxid)
    namedheir = OrderedDict()
    for taxid, rank, name in get_taxid_parents(taxid):
        if rank in ranks:
            namedheir[rank] = {'taxid': taxid, 'name': name}
    return namedheir

for i in sys.argv[1:]:
    print("gi ", i)
    print("taxid ", get_gi_taxid(int(i)))
    print(taxtree(get_gi_taxid(int(i))))

CON.close()
