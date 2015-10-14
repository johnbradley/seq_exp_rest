from sqlalchemy.engine import create_engine
from sqlalchemy import schema, types, orm
import os

class DB(object):
    def __init__(self, url=None):
        if not url:
            url = os.environ['DBURL']
        self.setup_engine(url)
        self.setup_project()

    def setup_engine(self, url):
        self.engine = create_engine(url) #'sqlite:///:memory:', echo=True)
        self.metadata = schema.MetaData()
        self.metadata.bind = self.engine
        orm.clear_mappers()
        sm = orm.sessionmaker(bind=self.engine, autoflush=True, autocommit=False,
            expire_on_commit=True)
        self.session = orm.scoped_session(sm)

    def setup_project(self):
        self.project_table = schema.Table('project', self.metadata,
            schema.Column('id', types.Integer, schema.Sequence('project_id_seq'),
                    primary_key=True),
            schema.Column('title', types.Unicode()),
        )
        orm.mapper(ProjectRecord, self.project_table)

    def create_tables(self):
        self.metadata.create_all(checkfirst=True)

    def delete_all(self):
        self.session.query(ProjectRecord).delete()

    def close(self):
        self.session.close()
        self.engine.dispose()
        orm.clear_mappers()

class ProjectRecord(object):
    def to_json(self):
        return {'id':self.id, 'title' : self.title}

class ProjectRecordList(object):
    def __init__(self, db):
        self.db = db
        self.project_table = self.db.project_table

    def save(self, project_record):
        self.db.session.add(project_record)
        self.db.session.commit()
        self.db.session.flush()

    def get_all(self):
        query = self.db.session.query(ProjectRecord)
        query.order_by(self.project_table.c.title)
        return list(query)

    def get_one(self, project_id):
        query = self.db.session.query(ProjectRecord)
        for item in query.filter(self.project_table.c.id == project_id):
            print(str(item.id))
            return item
        return None

    def delete_one(self, project_id):
        item = self.get_one(project_id)
        if item:
            self.db.session.delete(item)
            self.db.session.commit()
