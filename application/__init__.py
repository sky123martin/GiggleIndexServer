from flask import Flask
from application.config import config
# Libraries for scheduling
from flask_apscheduler import APScheduler
import datetime

app = Flask(__name__)
app.config.from_object(config)
app.config['TESTING'] = True

from application import utility
from application import routes

utility.setup_indices()

# Scheduling logic for updating
# scheduler = APScheduler()
# scheduler.add_job(func=update_indices, args=[], trigger='interval', id='job', seconds=app.config['UPDATE_INDEX_INTERVAL'])
# scheduler.start()
# app.run(port = 8000)