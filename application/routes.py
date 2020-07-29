from application import app
from flask import make_response

@app.route('/')
def home():
    return "Hello World!"

@app.route('/<page_name>')
def other_page(page_name):
    response = make_response('The page named %s does not exist.' \
                             % page_name, 404)
    return response
