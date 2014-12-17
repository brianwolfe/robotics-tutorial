from flask import Flask, render_template, abort
from flask.ext.sqlalchemy import SQLAlchemy
import os
from os.path import join
from urllib.parse import urljoin

app = Flask(__name__)
app.config.from_object(os.environ['APP_SETTINGS'])
db = SQLAlchemy(app)


tutorial_path = 'tutorials'
tutorials = os.listdir(tutorial_path)

tutorial_pages = {}
for tutorial in tutorials:
    # cut off .html
    if tutorial[-5:] == '.html':
        print("loading tutorial: ", tutorial)
        tutorial_name = tutorial[:-5]
        tutorial_pages[tutorial_name] = \
                open(join(tutorial_path, tutorial)).read()

@app.route('/', methods=['GET', 'POST'])
def index():
    return render_template('index.html')

@app.route('/tutorial/<path:tutorial>', methods=['GET'])
def tutorial(tutorial):
    print("Rendering: ", tutorial)
    print("Pages: ", list(tutorial_pages.keys()))
    print(app.config['STATIC_URL'])
    if tutorial in tutorial_pages:
        jsname = app.config['STATIC_URL'] +  '/js/' + tutorial + '.js'
        return render_template('tutorial.html',
            body=tutorial_pages[tutorial],
            jsname=jsname)
    else:
        abort(404)

# Sending static assets from Flask directly.
# I'll update this to send from CDN URLs and
# S3 URLs eventually, but this works for now.
@app.route('/static/<path:filename>')
def get_static_assets(filename):
    return send_from_directory('/static/', filename)

if __name__ == '__main__':
    app.run()
