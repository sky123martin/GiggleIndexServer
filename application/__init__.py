from flask import Flask
from application.config import config
# Libraries for scheduling
from flask_apscheduler import APScheduler
import datetime
# Libraries for sql db
from flask_sqlalchemy import SQLAlchemy
from flask_migrate import Migrate

app = Flask(__name__)
app.config.from_object(config)
app.config['TESTING'] = False
app.config['SQLALCHEMY_DATABASE_URI'] = 'sqlite:////tmp/test.db'
db = SQLAlchemy(app)
migrate = Migrate(app, db)

from application import models, utility, routes

db.drop_all()
db.create_all()
if app.config['TESTING'] == False:
    utility.setup_indices()

# Scheduling logic for updating
# scheduler = APScheduler()
# scheduler.add_job(func=update_indices, args=[], trigger='interval', id='job', seconds=app.config['UPDATE_INDEX_INTERVAL'])
# scheduler.start()
# app.run(port = 8000)