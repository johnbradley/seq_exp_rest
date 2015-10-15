from flask import Flask, request
from flask_restful import abort, Api, Resource
from Bio import Entrez, SeqIO

from sqlalchemy.engine import create_engine
from sqlalchemy import schema, types, orm
import os

from seq_exp.db import DB, ProjectRecord, ProjectRecordList

Entrez.email = os.environ['EMAIL']


class Project(Resource):
    def __init__(self, db=None):
        self.project_list = ProjectRecordList(db)

    def get(self, project_id):
        proj = self.project_list.get_one(project_id);
        if proj:
            return proj.to_json()
        abort(404, message="Project {} doesn't exist".format(project_id))

    def delete(self, project_id):
        self.project_list.delete_one(project_id)
        return '', 204

    def put(self, project_id):
        proj = self.project_list.get_one(project_id);
        if proj:
            args = request.values
            proj.title = args['title']
            self.project_list.save(proj)
            return proj.to_json(), 201
        abort(404, message="Project {} doesn't exist".format(project_id))


class ProjectList(Resource):
    def __init__(self, db=None):
        self.project_list = ProjectRecordList(db)

    def get(self):
        data = [proj.to_json() for proj in self.project_list.get_all()]
        cnt = len(data)
        return {'items':data}, 201, {'Access-Control-Allow-Origin': '*'} 

    def post(self):
        args = request.values
        proj = ProjectRecord()
        proj.title = args['title']
        self.project_list.save(proj)
        return proj.to_json(), 201


class EntrezSummary(Resource):
    def __init__(self, db=None):
        self.db = db

    def get(self, db):
        return self.get_summary(db, request.values)

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
    def __init__(self, db=None):
        self.db = db

    def get(self, db, id):
        handle = Entrez.esummary(db=db, id=id)
        record = Entrez.read(handle)
        return record[0]


class EntrezDownload(Resource):
    def __init__(self, db=None):
        self.db = db

    def post(self):
        ids = request.form['ids']
        project_id = request.values['project_id']
        db = request.values['db']
        #print(ids, project_id,)
        return {'result':ids + db + project_id}




def setup_api_and_db(url):
    app = Flask(__name__, static_folder="website")
    api = Api(app)
    mydb = DB(url=url)
    api.add_resource(ProjectList, '/projects',
            resource_class_kwargs={ 'db': mydb })
    api.add_resource(Project, '/projects/<project_id>',
            resource_class_kwargs={ 'db': mydb })
    api.add_resource(EntrezSummary, '/entrez/<db>',
            resource_class_kwargs={ 'db': mydb })
    api.add_resource(EntrezDetail, '/entrez/<db>/<id>',
            resource_class_kwargs={ 'db': mydb })
    api.add_resource(EntrezDownload, '/fetch_request',
            resource_class_kwargs={ 'db': mydb })
    return app, mydb
