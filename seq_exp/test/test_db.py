import unittest
from seq_exp.seq_exp import DB, ProjectRecord, ProjectRecordList

class ProjectRecordListTest(unittest.TestCase):
    def setUp(self):
        self.db = DB(url='sqlite:///:memory:')
        self.db.create_tables()
        self.project_list = ProjectRecordList(self.db)

    def tearDown(self):
        self.db.close()

    def assert_project_cnt(self, expected_cnt):
        cnt = len(self.project_list.get_all())
        self.assertEqual(expected_cnt, cnt)

    def save_project(self, title):
        proj = ProjectRecord()
        proj.title = title
        self.project_list.save(proj)
        return proj.id

    def test_create_one(self):
        self.assert_project_cnt(0)
        project_id = self.save_project("Some Title")
        self.assert_project_cnt(1)
        proj2 = self.project_list.get_one(project_id)
        self.assertEqual("Some Title", proj2.title)

    def test_create_two_and_delete_one(self):
        self.assert_project_cnt(0)
        project_id = self.save_project("First Title")
        project_id2 = self.save_project("Second Title")
        self.assert_project_cnt(2)
        self.project_list.delete_one(project_id)
        self.assert_project_cnt(1)

    def test_update(self):
        project_id = self.save_project("First Title")
        proj = self.project_list.get_one(project_id)
        self.assertEqual("First Title", proj.title)
        proj.title = "Second Title"
        proj.id = project_id
        self.project_list.save(proj)
        proj2 = self.project_list.get_one(project_id)
        self.assertEqual("Second Title", proj2.title)
