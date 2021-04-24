from setup_indices import connect_SQL_db, extract_bed_columns
from config import config
from contextlib import contextmanager
import pandas as pd
import unittest
import mysql.connector
import requests
import random

config = config

class BasicTests(unittest.TestCase):

    ######################
    # setup and teardown #
    ######################

    # executed prior to each test
    def setUp(self):
        config.TESTING = True

    # executed after each test
    def tearDown(self):
        pass

    #############################
    # UCSC API Structure Tests #
    #############################
    # If any of these tests break try to check and see if UCSC has
    # modified API structure.
    def test_UCSC_API_genomes(self):
        # testing if genome directory exists
        r = requests.get(config.UCSC_API+'/list/ucscGenomes')
        self.assertEqual(r.status_code, requests.codes.ok)
        r = requests.get(config.UCSC_API+'/list/fakeURL')
        self.assertEqual(r.status_code, 400)

    def test_UCSC_API_track_structure(self):
        # testing if track structure holds
        genome = "rn6"
        r = requests.get('{}/list/tracks?genome={};trackLeavesOnly=1'.format(config.UCSC_API, "rn6"))
        self.assertEqual(r.status_code, requests.codes.ok)
        rn6_leaves = r.json()["rn6"]
        r = requests.get('{}/list/tracks?genome={}'.format(config.UCSC_API, "rn6"))
        rn6_full = r.json()["rn6"]
        self.assertTrue(len(rn6_leaves) > len(rn6_full))

    def test_UCSC_API_all_tracks(self):
        # testing that /genomes listed can link to /tracks
        r = requests.get(config.UCSC_API+'/list/ucscGenomes')
        self.assertEqual(r.status_code, requests.codes.ok)
        all_genomes = r.json()["ucscGenomes"]
        # Only chosing acouple to check because hg19 delay, normal delay
        for key in random.choices(list(all_genomes.keys()), k=3):
            track = requests.head('{}/list/tracks?genome={}'.format(config.UCSC_API, key), timeout=45.00)
            self.assertEqual(track.status_code, requests.codes.ok)

    #############################
    # UCSC MySQL Database Tests #
    #############################
    def test_false_SQL_connection(self):
        # test that false connection throws connection error
        with self.assertRaises(mysql.connector.Error):
            host = "fake"
            user = "password"
            with connect_SQL_db(host, user):
                pass

    def test_UCSC_SQL_connection(self):
        # test that current UCSC connection works
        host = config.UCSC_SQL_DB_HOST
        user = config.UCSC_SQL_DB_USER
        with connect_SQL_db(host, user) as db:
            mycursor = db.cursor()
            mycursor.execute("SHOW DATABASES")
            self.assertTrue(mycursor)

    def test_UCSC_SQL_databases(self):
        # test thatincorrect query throws
        host = config.UCSC_SQL_DB_HOST
        user = config.UCSC_SQL_DB_USER
        with connect_SQL_db(host, user) as db:
            mycursor = db.cursor()
            mycursor.execute("SELECT SCHEMA_NAME FROM INFORMATION_SCHEMA.SCHEMATA\
                              WHERE SCHEMA_NAME = 'hg19'")
        self.assertTrue(mycursor)

    ######################
    # Data Scapper Tests #
    ######################
    def test_extract_bed_columns(self):
        # Local test of extract_bed_columns
        file_type = "chain ailMel1"
        columns = ['bin', 'score', 'tName', 'tSize', 'tStart', 'tEnd', 'qName', 'qSize', 'qStrand', 'qStart', 'qEnd', 'id']
        expected_mapping = {"tName": "chrom", "tStart":"start", "tEnd":"end"}
        mapping, bs = extract_bed_columns(columns, file_type)
        self.assertEqual(expected_mapping, mapping)
        self.assertEqual(bs, 0)

        file_type = "netAlign panTro5 chainPanTro5"
        columns = ['bin', 'level', 'tName', 'tStart', 'tEnd', 'strand', 'qName', 'qStart', 'qEnd', 'chainId', 'ali', 'score', 'qOver', 'qFar', 'qDup', 'type', 'tN', 'qN', 'tR', 'qR', 'tNewR', 'qNewR', 'tOldR', 'qOldR', 'tTrf', 'qTrf']
        expected_mapping = {"tName":"chrom", "tStart":"start", "tEnd":"end"}
        mapping, bs = extract_bed_columns(columns, file_type)
        self.assertEqual(expected_mapping, mapping)
        self.assertEqual(bs, 0)

        file_type = "bed 4 +"
        columns = ['bin', 'chrom', 'chromStart', 'chromEnd', 'name', 'period', 'copyNum', 'consensusSize', 'perMatch', 'perIndel', 'score', 'A', 'C', 'G', 'T', 'entropy', 'sequence']
        expected_mapping = {"chrom": "chrom", "chromStart":"start", "chromEnd":"end"}
        mapping, bs = extract_bed_columns(columns, file_type)
        self.assertEqual(expected_mapping, mapping)
        self.assertEqual(bs, 0)

        file_type = "bed 4"
        columns = ['chrom', 'chromStart', 'chromEnd', 'name']
        expected_mapping = {"chrom":"chrom", "chromStart":"start", "chromEnd": "end"}
        mapping, bs = extract_bed_columns(columns, file_type)
        self.assertEqual(expected_mapping, mapping)
        self.assertEqual(bs, 0)

        file_type = "psl est"
        columns = ['bin', 'matches', 'misMatches', 'repMatches', 'nCount', 'qNumInsert', 'qBaseInsert', 'tNumInsert', 'tBaseInsert', 'strand', 'qName', 'qSize', 'qStart', 'qEnd', 'tName', 'tSize', 'tStart', 'tEnd', 'blockCount', 'blockSizes', 'qStarts', 'tStarts']
        expected_mapping = {"tName":"chrom", "tStart":"start", "tEnd":"end"}
        mapping, bs = extract_bed_columns(columns, file_type)
        self.assertEqual(expected_mapping, mapping)
        self.assertEqual(bs, 0)

        file_type = "genePred genscanPep"
        columns = ['bin', 'name', 'chrom', 'strand', 'txStart', 'txEnd', 'cdsStart', 'cdsEnd', 'exonCount', 'exonStarts', 'exonEnds']
        expected_mapping = {"chrom": "chrom", 'txStart':"start", 'txEnd':"end"}
        mapping, bs = extract_bed_columns(columns, file_type)
        self.assertEqual(expected_mapping, mapping)
        self.assertEqual(bs, 0)

        file_type = "rmsk"
        columns = ['bin', 'swScore', 'milliDiv', 'milliDel', 'milliIns', 'genoName', 'genoStart', 'genoEnd', 'genoLeft', 'strand', 'repName', 'repClass', 'repFamily', 'repStart', 'repEnd', 'repLeft', 'id']
        expected_mapping = {'genoName':"chrom", 'genoStart':"start", 'genoEnd':"end"}
        mapping, bs = extract_bed_columns(columns, file_type)
        self.assertEqual(expected_mapping, mapping)
        self.assertEqual(bs, 0)

        file_type = "GFF"
        columns = ["start", "end", "chrom"]
        expected_mapping = {"chrom": "chrom", "start": "start", "end": "end"}
        mapping, bs = extract_bed_columns(columns, file_type)
        self.assertEqual(expected_mapping, mapping)
        self.assertEqual(bs, -1)

        file_type = "MAF"
        columns = ["start", "end", "chrom"]
        expected_mapping = {"chrom": "chrom", "start": "start", "end": "end"}
        mapping, bs = extract_bed_columns(columns, file_type)
        self.assertEqual(expected_mapping, mapping)
        self.assertEqual(bs, 0)

    ####################
    # Config Variables #
    ####################

    def test_CONFIG_VAR_MAX_INTERVAL_PER_INDEX(self):
        """ Failure of this test implies that MAX_INTERVALS_PER_INDEX is negative"""
        self.assertGreater(config.MAX_INTERVALS_PER_INDEX, 0)
        # genome = "rn6"
        # for i in 
        # tracks_info = requests.get(url="{}/list/tracks?genome={};trackLeavesOnly=1".format(config.UCSC_API, genome)).json()[genome]
        # # unzip track details and names
        # for track, track_info in tracks_info.items():
        #     self.assertGreaterEqual(config.MAX_INTERVALS_PER_INDEX, track_info["itemCount"])
  
    ##################
    # Indexing Tests #
    ##################

    # def test_cluster_data(self):
    #     testing_size = config.MAX_INTERVALS_PER_INDEX/2
    #     files_info = {"bannana": {
    #                 "file_size": testing_size, 
    #                 "download_function": "ba", 
    #                 "download_params": [1,1,1]
    #                 },
    #              "mellon": {
    #                 "file_size": testing_size * 1.1, 
    #                 "download_function": "", 
    #                 "download_params": [2,2,2]
    #                 },

    #              "apple": {
    #                 "file_size": testing_size, 
    #                 "download_function": "", 
    #                 "download_params": [3,3,3]
    #             }}
    #     genome = "fruits"
    #     source = "kitchen"
    #     clusters = cluster_data(source, genome, files_info)
    #     expected_clusters = {source + "_" + genome+"_1": {
    #                             "full": True,
    #                             "files": {
    #                                 "bannana": {
    #                                     "file_size": testing_size, 
    #                                     "download_function": "ba", 
    #                                     "download_params": [1,1,1]
    #                                     },
    #                                 "mellon": {
    #                                     "file_size": testing_size * 1.1, 
    #                                     "download_function": "", 
    #                                     "download_params": [2,2,2]
    #                                 }},
    #                             "index_size": testing_size + 1.1*testing_size},
    #                          source + "_" + genome+"_2": {
    #                              "full": False,
    #                              "files": {
    #                                 "apple": {
    #                                     "file_size": testing_size, 
    #                                     "download_function": "", 
    #                                     "download_params": [3,3,3]
    #                                 }},
    #                             "index_size": testing_size}
    #                       }                     
    #     self.assertTrue(expected_clusters == clusters)


if __name__ == "__main__":
    unittest.main()
