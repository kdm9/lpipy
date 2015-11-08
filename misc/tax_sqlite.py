#!/usr/bin/env python3

from collections import OrderedDict
import gzip
import sqlite3
import sys
import tarfile


def setup_db(dbfile="lpi.db"):
    conn = sqlite3.connect(dbfile)
    c = conn.cursor()
    c.execute("""CREATE TABLE merges (
                    oldid INTEGER PRIMARY KEY,
                    newid INTEGER
                    )""")
    c.execute("""CREATE TABLE sequences (
                    gi INTEGER PRIMARY KEY,
                    taxid INTEGER,
                    FOREIGN KEY (taxid) REFERENCES taxa(taxid)
                    )""")
    c.execute("""CREATE TABLE taxa (
                    taxid INTEGER PRIMARY KEY,
                    name TEXT,
                    rank TEXT
                    )""")
    c.execute("""CREATE TABLE parents (
                    taxid INTEGER,
                    parent INTEGER,
                    level INTEGER,
                    FOREIGN KEY (taxid) REFERENCES taxa(taxid),
                    FOREIGN KEY (parent) REFERENCES taxa(taxid)
                    )""")
    c.execute("CREATE INDEX parents_taxid ON parents (taxid)")
    conn.commit()
    return conn

CON = setup_db("lpi-new.db")
#CON = setup_db("lpi-dev.db")

def load_tax(tarball):
    tar = tarfile.open(tarball)
    db = CON.cursor()

    def dmpfile(filename):
        with tar.extractfile(filename) as fh:
            for line in fh:
                line = line.decode("utf-8")
                fields = [x.strip() for x in line.strip('\n|').split('\t|\t')]
                yield fields

    names = {}
    nodes = {}
    ranks = {}
    merges = {}

    for taxon in dmpfile('names.dmp'):
        taxid = int(taxon[0])
        name = taxon[1]
        nameclass = taxon[3]
        if nameclass == "scientific name":
            names[taxid] = name

    for taxon in dmpfile('nodes.dmp'):
        taxid = int(taxon[0])
        parent = int(taxon[1])
        rank = taxon[2]
        nodes[taxid] = parent
        ranks[taxid] = rank

    for ids in dmpfile("merged.dmp"):
        oldid, newid = map(int, ids)
        db.execute("INSERT INTO merges VALUES (?, ?)", (oldid, newid))

    print("Loaded data files")

    # The tree is stored a a dict of {taxid: [taxid, parent, ..., 1]}
    # for each parent in the tree to the root.
    # Start with the root node.
    tree = {1: [], }
    for taxon, parent in nodes.items():
        if parent in tree and taxon != 1:
            tree[taxon] = [taxon] + tree[parent]
        else:
            subtree = [taxon, parent]
            # A new list for just the new section of the tree
            newtree = [taxon, parent]
            thisid = parent
            verb = False
            while True:
                if thisid in tree:
                    # This will eventually be true, so long as everything
                    # is rooted at taxid 1
                    subtree.extend(tree[thisid])
                    break
                thisparent = nodes[thisid]
                subtree.append(thisparent)
                newtree.append(thisparent)

                thisid = thisparent
            # Add the new section of the tree
            for i in range(len(newtree)):
                thisid = subtree[i]
                heirarchy = subtree[i + 1:]
                tree[thisid] = heirarchy

    print("Created tree (with", len(tree), "records)")
    taxa = []
    count = 0
    for taxid, parents in tree.items():
        taxa.append({
                "_id": taxid,
                "parents": parents,
                "name": names[taxid],
                "rank": ranks[taxid],
            })
        parent_items = [(taxid, par, lvl) for lvl, par in  enumerate(parents)]
        db.executemany("INSERT INTO parents VALUES (?, ?, ?)", parent_items)
        db.execute("INSERT INTO taxa VALUES (?, ?, ?)",
                   (taxid, names[taxid], ranks[taxid]))
        if len(taxa) > 10000:
            count += 10
            print("\33[2KInserted ", count, "K records", sep="", end='\r')
            taxa = []

    print("Loaded tree into DB")
    print("Finished loading taxonomy")
    CON.commit()

load_tax("taxdump.tar.gz")

def load_gi_to_tax(file):
    db = CON.cursor()
    with gzip.open(file) as fh:
        pairs = []
        for i, line in enumerate(fh):
            gi, taxid = map(int, line.strip().split())
            pairs.append((gi, taxid))
            if i % 100000 == 0:
                db.executemany("INSERT INTO sequences VALUES (?, ?)", pairs)
                pairs[:] = []
                print("\33[2K\r", i/1000000, "M GIs", sep='', end='\r')
                sys.stdout.flush()
    print("\nFinished loading GIs")
    CON.commit()

#load_gi_to_tax('gi_taxid_nucl.dmp.gz')
load_gi_to_tax('gi_taxid_prot.dmp.gz')

CON.close()
