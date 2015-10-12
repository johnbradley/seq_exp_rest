from flask import Flask
from flask_restful import reqparse, abort, Api, Resource
from Bio import Entrez, SeqIO
import os

Entrez.email = os.environ['EMAIL']

app = Flask(__name__)
api = Api(app)

PROJECTS = {
}

parser = reqparse.RequestParser()
parser.add_argument('title')
parser.add_argument('term')
parser.add_argument('retmax')

# Todo
# shows a single todo item and lets you delete a todo item
class Project(Resource):
    def get(self, project_id):
        self.abort_if_project_doesnt_exist(project_id)
        return PROJECTS[project_id]

    def delete(self, project_id):
        self.abort_if_project_doesnt_exist(project_id)
        del PROJECTS[project_id]
        return '', 204

    def put(self, project_id):
        args = parser.parse_args()
        project = {'id': project_id, 'title': args['title']}
        PROJECTS[project_id] = project
        return project, 201

    def abort_if_project_doesnt_exist(self, project_id):
        if project_id not in PROJECTS:
            abort(404, message="Project {} doesn't exist".format(project_id))

# TodoList
# shows a list of all todos, and lets you POST to add new tasks
class ProjectList(Resource):
    def get(self):
        cnt = len(PROJECTS)
        return {'count':cnt, 'data':list(PROJECTS.values())}

    def post(self):
        args = parser.parse_args()
        project_id = len(PROJECTS.keys()) + 1
        project_id = '%i' % project_id
        PROJECTS[project_id] = {'id': project_id, 'title': args['title']}
        return PROJECTS[project_id], 201

class EntrezSummary(Resource):
    def get(self, db):
        return self.get_summary(db, parser.parse_args())

    def get_summary(self, db, args):
        term = str(args['term'])
        retmax = 10
        retmax_str = args['retmax']
        if retmax_str:
            retmax = int(retmax_str)
        handle = Entrez.esearch(db=db, term=term, retmax=retmax)
        record = Entrez.read(handle)
        total_cnt = record['Count']
        gi_str = ",".join(record["IdList"])
        count = len(record["IdList"])
        handle = Entrez.esummary(db=db, id=gi_str)
        record = Entrez.read(handle)
        return {'total_cnt': total_cnt, 'count': count, 'data': record}

class EntrezDetail(Resource):
    def get(self, db, id):
        handle = Entrez.esummary(db=db, id=id)
        record = Entrez.read(handle)
        return record[0]

##
## Actually setup the Api resource routing here
##
api.add_resource(ProjectList, '/projects')
api.add_resource(Project, '/projects/<project_id>')
api.add_resource(EntrezSummary, '/entrez/<db>')
api.add_resource(EntrezDetail, '/entrez/<db>/<id>')
#api.add_resource(EntrezProtein, '/entrez/protein')
#/entrez/nucleotide/summary?term=<search_term>&retmax=<retmax> - search for DNA info
#/project - list/create project
#/project/<project_id> - view/create specific project

if __name__ == '__main__': # pragma: no cover
    app.run(debug=True)
