from datetime import datetime
from application import db


class Index(db.Model):
    id = db.Column(db.String(100), primary_key=True)
    timestamp = db.Column(db.DateTime, index=True, default=datetime.utcnow)
    source = db.Column(db.String(64), index=True)
    genome = db.Column(db.String(64), index=True)
    size = db.Column(db.Integer)
    full = db.Column(db.Boolean)

    def __repr__(self):
        return '<Index {} {}>'.format(self.id, self.full)


class Files(db.Model):
    tablename = db.Column(db.String(55))
    source = db.Column(db.String(55), index=True)
    genome = db.Column(db.String(64), index=True)
    size = db.Column(db.Integer)
    timestamp = db.Column(db.DateTime, index=True, default=datetime.utcnow)
    index_id = db.Column(db.String(100))
    id = db.Column(db.Integer, primary_key=True, autoincrement=True)

    def __repr__(self):
        return '<File {}>'.format(self.tablename)
