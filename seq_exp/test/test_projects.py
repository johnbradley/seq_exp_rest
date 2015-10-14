import os
import seq_exp.seq_exp as seq_exp
import unittest
import tempfile
import ast

class ProjectsTestCase(unittest.TestCase):
    def setUp(self):
        realapp, db = seq_exp.setup_api_and_db('sqlite:///:memory:')
        self.realapp = realapp
        self.db = db
        self.db.create_tables()
        self.db_fd, self.realapp.config['DATABASE'] = tempfile.mkstemp()
        self.realapp.config['TESTING'] = True
        self.app = self.realapp.test_client()

    def tearDown(self):
        self.db.close()
        os.close(self.db_fd)
        os.unlink(self.realapp.config['DATABASE'])

    def literal_eval(self, rv):
        resp_str = rv.data.decode("utf-8")
        return ast.literal_eval(resp_str)

    def get_projects_data_and_count(self):
        rv = self.app.get('/projects')
        resp = self.literal_eval(rv)
        data = {}
        for project in resp['data']:
            data[project['id']] = project['title']
        return data, resp['count']

    def post_and_get_id(self, title):
        ret = self.app.post('/projects', data=dict(title=title))
        return self.literal_eval(ret)['id']

    def get_project_by_id(self, id):
        ret = self.app.get('/projects/' + str(id))
        if ret.status.startswith("404 "):
            return None
        return self.literal_eval(ret)

    def put_project_by_id(self, id, data):
        self.app.put('/projects/' + str(id), data=data)

    def delete_project_by_id(self, id):
        self.app.delete('/projects/' + str(id))

    def test_empty_projects(self):
        # Test GET /projects.
        data, count = self.get_projects_data_and_count()
        print(count)
        self.assertEqual(count, 0)
        self.assertEqual(data, {})

    def test_add_one(self):
        # Test POST /projects.
        id = self.post_and_get_id("Test 1")
        data, count = self.get_projects_data_and_count()
        self.assertEqual(count, 1)
        self.assertEqual(data[id], "Test 1")

    def test_add_two(self):
        # Test POST /projects twice.
        red_id = self.post_and_get_id("Red")
        green_id = self.post_and_get_id("Green")
        data, count = self.get_projects_data_and_count()
        self.assertEqual(count, 2)
        self.assertEqual(data[red_id], "Red")
        self.assertEqual(data[green_id], "Green")

    def test_add_two_remove_one(self):
        # Test DELETE /projects/<id>.
        red_id = self.post_and_get_id("Red")
        green_id = self.post_and_get_id("Green")
        self.delete_project_by_id(red_id)
        data, count = self.get_projects_data_and_count()
        self.assertEqual(count, 1)
        self.assertEqual(data[green_id], "Green")

    def test_add_then_get(self):
        # Test GET /projects/<id>.
        id = self.post_and_get_id("Test 1")
        ret = self.get_project_by_id(id)
        self.assertEqual(ret['id'], id)
        self.assertEqual(ret['title'], "Test 1")

    def test_add_then_put_then_get(self):
        # Test PUT /projects/<id>.
        id = self.post_and_get_id("Test 1")
        self.put_project_by_id(id, dict(title='Test 2'))
        ret = self.get_project_by_id(id)
        self.assertEqual(ret['id'], id)
        self.assertEqual(ret['title'], "Test 2")

    def test_get_invalid(self):
        ret = self.get_project_by_id(99)
        self.assertEqual(None, ret)

if __name__ == '__main__':
    unittest.main()
