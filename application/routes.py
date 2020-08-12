from application import app
from flask import make_response

@app.route('/')
def home():
    return "Hello World!"

@app.route('/<page_name>')
def other_page(page_name):
    pass

@app.route('/sources')
def sources(page_name):
    pass

@app.route('/{}/genomes')
def genomes(source):
    pass

@app.route('/{}/{}/tracks')
def genomes(source):
    pass

@app.route('/{}/{}/indices')
def genomes(source):
    pass