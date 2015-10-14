import os
import seq_exp.seq_exp as seq_exp
import unittest
import tempfile
import ast

class EntrezTestCase(unittest.TestCase):
    def setUp(self):
        realapp, db = seq_exp.setup_api_and_db('sqlite:///:memory:')
        self.realapp = realapp
        self.db = db
        self.db_fd, realapp.config['DATABASE'] = tempfile.mkstemp()
        realapp.config['TESTING'] = True
        self.app = realapp.test_client()
        seq_exp.PROJECTS = {}

    def tearDown(self):
        self.db.close()
        os.close(self.db_fd)
        os.unlink(self.realapp.config['DATABASE'])

    def literal_eval(self, rv):
        resp_str = rv.data.decode("utf-8")
        return ast.literal_eval(resp_str)

    def test_fetch_four_human_dna(self):
        #kind of fragile since relies upon external web server
        rv = self.app.get('/entrez/nucleotide', data=dict(term='human', retmax='4'))
        resp = self.literal_eval(rv)
        self.assertEqual(4, resp['count'])

    def test_fetch_five_mouse_protein(self):
        #kind of fragile since relies upon external web server
        rv = self.app.get('/entrez/protein', data=dict(term='mouse', retmax='5'))
        resp = self.literal_eval(rv)
        self.assertEqual(5, resp['count'])
