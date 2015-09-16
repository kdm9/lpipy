"""
NCBI Taxonomy database handling.
"""

from __future__ import absolute_import, division, print_function
from collections import OrderedDict
import gzip
import json
import os
from os import path
import sqlite3
import tarfile

import requests

from .utils import (
    download_file,
    get_data_dir,
    get_data_file,
    get_logger,
    md5sum,
    read_remote_file,
)


__all__ = ["NCBITaxonomyDB"]

DB_VERSION = 1
DEFAULT_RANKS = [
    'species',
    'genus',
    'family',
    'order',
    'class',
    'phylum',
    'kingdom',
    'superkingdom'
]

LOG = get_logger()

class NCBITaxonomyDB(object):
    '''
    Provides a local transparent connector to the NCBI taxonomy database.
    '''
    _db = None
    metadata = {
        "db_version": DB_VERSION,
        "db_complete": False,
        "taxdump_md5": '',
        "gi_taxid_nucl_md5": '',
        "gi_taxid_prot_md5": '',
    }


    def __init__(self, dbfile=None):

        self._dbfile = dbfile
        if not self._dbfile:
            self._dbfile = get_data_file('taxa.sqlite')

        self._connect()

    def _connect(self):
        if self._db:
            # No-op if we have a db handle
            return
        try:
            self._db = sqlite3.connect(self._dbfile)
            with open(self._dbfile + '.info') as fh:
                self.metadata = json.load(fh)
            if self.metadata["db_version"] != DB_VERSION:
                raise ValueError("Invalid DB version")
            if not self.metadata["db_complete"]:
                raise ValueError("Incomplete DB")
            LOG.debug('Taxonomy DB loaded successfully')
        except Exception:
            LOG.info('Taxonomy DB not valid, (re-)creating')
            self._wipe_db()
            self._create_db()

    def _write_metadata(self):
        with open(self._dbfile + '.info', 'w') as fh:
            json.dump(self.metadata, fh)

    def _wipe_db(self):
        LOG.debug('Wiping sqlite taxonomy database')
        try:
            os.remove(self._dbfile)
        except Exception:
            pass
        try:
            os.remove(self._dbfile + '.info')
        except Exception:
            pass

    def _create_db(self):
        """Creates the structure in the sqlite DB and opens a connection to the
        database."""
        LOG.debug('Creating sqlite taxonomy database')
        self._db = sqlite3.connect(self._dbfile)
        c = self._db.cursor()
        c.execute("""CREATE TABLE merges (
                        oldid INTEGER PRIMARY KEY,
                        newid INTEGER
                        );""")
        c.execute("""CREATE TABLE sequences (
                        gi INTEGER PRIMARY KEY,
                        taxid TEXT,
                        FOREIGN KEY (taxid) REFERENCES taxa(taxid)
                        );""")
        c.execute("""CREATE TABLE taxa (
                        taxid INTEGER PRIMARY KEY,
                        name TEXT,
                        rank TEXT
                        );""")
        c.execute("""CREATE TABLE parents (
                        taxid INTEGER,
                        parent INTEGER,
                        level INTEGER,
                        FOREIGN KEY (taxid) REFERENCES taxa(taxid),
                        FOREIGN KEY (parent) REFERENCES taxa(taxid)
                        );""")
        c.execute("CREATE INDEX parents_taxid ON parents (taxid)")
        self._db.commit()
        self._write_metadata()
        self._load_db()

    def _get_taxdump_file(self):
        '''Helper to download the NCBI taxdump tarball if we need it'''
        taxdump_url = 'http://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz'
        taxdump_file = get_data_file('taxdump.tar.gz')

        taxdump_remote_hash = read_remote_file(taxdump_url + '.md5').split()[0]

        if path.exists(taxdump_file):
            LOG.debug("Taxdump file found, calculating checksum")
            self.metadata['taxdump_md5'] = md5sum(taxdump_file)

        while taxdump_remote_hash != self.metadata.get('taxdump_md5', ''):
            LOG.info("Taxdump file out of date, downloading it")
            download_file(taxdump_url, taxdump_file)
            self.metadata['taxdump_md5'] = md5sum(taxdump_file)

        LOG.debug("We have a valid Taxdump file")

        self._write_metadata()

    def _get_gitax_file(self):
        '''Helper to download the NCBI gi_taxid mappings if we need to.'''
        url = 'http://ftp.ncbi.nlm.nih.gov/pub/taxonomy/gi_taxid_{sect}.dmp.gz'
        gi_tax_file = get_data_file('gi_taxid_{sect}.dmp.gz')

        for sect in ['prot', 'nucl']:
            sect_url = url.format(sect=sect)
            sect_remote_hash = read_remote_file(sect_url + '.md5').split()[0]
            sect_file = gi_tax_file.format(sect=sect)
            sect_hash_key = 'gi_taxid_{}_md5'.format(sect)
            sect_local_hash = self.metadata.get(sect_hash_key, '')

            if path.exists(sect_file):
                LOG.debug("GI to tax %s file found, calculating checksum",
                          sect)
                self.metadata[sect_hash_key] = md5sum(sect_file)

            while self.metadata[sect_hash_key] != sect_remote_hash:
                download_file(sect_url, sect_file)
                self.metadata[sect_hash_key] = md5sum(sect_file)
            LOG.debug("We have a valid gi_taxid_%s file", sect)
        self._write_metadata()

    def _load_db(self):
        # up to here, transfer it from the other files from ~/ws/lpi
        LOG.debug('Loading taxonomy database from NCBI dump files')
        self._get_taxdump_file()
        self._get_gitax_file()

        taxdump_file = get_data_file('taxdump.tar.gz')
        tar = tarfile.open(taxdump_file)
        cursor = self._db.cursor()

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
            cursor.execute("INSERT INTO merges VALUES (?, ?)", (oldid, newid))

        LOG.debug("Loaded data files")

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

        LOG.debug("Created taxa tree")

        for i, (taxid, parents) in enumerate(tree.items()):
            parent_items = [(taxid, par, i) for i, par in enumerate(parents)]
            cursor.executemany("INSERT INTO parents VALUES (?, ?, ?)",
                               parent_items)
            cursor.execute("INSERT INTO taxa VALUES (?, ?, ?)",
                           (taxid, names[taxid], ranks[taxid]))
            if i % 100000 == 0:
                LOG.debug("Inserted %0.2fM taxa", i / 1000000)

        LOG.debug("Created taxa database (with %d records)", len(tree))
        self._db.commit()

        LOG.debug("Creating sequence GI to taxid database")
        for sect in ['prot', 'nucl']:
            filename = get_data_file('gi_taxid_{}.dmp.gz').format(sect)
            LOG.debug("Loading 'gi_taxid_%s.dmp.gz'", sect)
            with gzip.open(filename) as fh:
                pairs = []
                for i, line in enumerate(fh):
                    gi, taxid = map(int, line.strip().split())
                    pairs.append((gi, taxid))
                    if i % 1000000 == 0:
                        cursor.executemany(
                            "INSERT INTO sequences VALUES (?, ?)",
                            pairs
                        )
                        pairs[:] = []
                        LOG.debug("%0.0fM %s GIs loaded", i/1000000, sect)
            self._db.commit()
        LOG.debug("Finished loading sequences ID to taxid database")
        LOG.info("Loaded taxonomy and sequence data into database")
        self.metadata["db_complete"] = True
        self._write_metadata()

    def get_merged_taxid(self, taxid):
        '''Convert an old taxid to it's up-to-date ID, or return the orginal
        taxid'''
        qry = self._db.execute('''SELECT newid FROM merges
                                 WHERE oldid = ?''', (taxid, ))
        newid = qry.fetchone()
        if newid:
            return newid
        return taxid

    def gi_to_taxid(self, gi):
        # Ensure GI is an integer
        gi = int(gi)

        c = self._db.cursor()
        qry = c.execute("SELECT taxid FROM sequences WHERE gi = ?", (gi,))
        res = qry.fetchone()
        if not res:
            raise ValueError("GI {} not found in database".format(gi))
        return res[0]

    def taxon_parents(self, taxid):
        c = self._db.cursor()
        qry = c.execute("""SELECT parents.parent, taxa.rank, taxa.name
                           FROM parents
                           INNER JOIN taxa ON taxa.taxid = parents.parent
                           WHERE parents.taxid = ?
                           ORDER BY parents.level;""",
                        (taxid, ))
        for parent, rank, name in qry:
            yield (parent, rank, name)

    def taxon_lineage(self, taxid, ranks=DEFAULT_RANKS):
        ranks = set(ranks)
        taxid = self.get_merged_taxid(taxid)
        namedheir = OrderedDict()
        for taxid, rank, name in self.taxon_parents(taxid):
            if rank in ranks:
                namedheir[rank] = {'taxid': taxid, 'name': name}
        return namedheir

    def sequence_lineage(self, gi, ranks=DEFAULT_RANKS):
        taxid = self.gi_to_taxid(gi)
        return self.taxon_lineage(taxid)

