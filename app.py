from flask import Flask, render_template, abort, jsonify, request
from flask.ext.sqlalchemy import SQLAlchemy
import os
from os.path import join
from urllib.parse import urljoin
import numpy as np

app = Flask(__name__)
app.config.from_object(os.environ['APP_SETTINGS'])
db = SQLAlchemy(app)

tutorial_path = 'tutorials'
tutorials = []
tutorial_pages = {}

def load_tutorials():
    global tutorials
    global tutorial_pages
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
    if app.config['DEBUG']:
        load_tutorials()
    print("Rendering: ", tutorial)
    print("Pages: ", list(tutorial_pages.keys()))
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


ROBOT_RADIUS = 1

def movement_update(left_wheel, right_wheel, robot):
    mean = np.array(robot['mu'], dtype=np.float64)
    cov = np.array(robot['sigma'], dtype=np.float64)

    dd = (left_wheel + right_wheel) / 2;
    dt = (right_wheel - left_wheel) / ROBOT_RADIUS

    x1 = mean[0]
    y1 = mean[1]
    theta1 = mean[2]

    x2 = x1 + dd * np.cos(theta1)
    y2 = y1 + dd * np.sin(theta1)
    theta2 = theta1 + dt

    cov[0,0] += 1
    cov[1,1] += 1
    cov[2,2] += 0.1
    mean2 = np.array([x2, y2, theta2]).reshape((3, 1))
    print(mean)
    print(mean2)

    return {'mu': mean2.tolist(), 'sigma': cov.tolist()}

# TODO split the robots out into different modules
# staying here for now
@app.route('/robotupdate/<path:filtertype>', methods=['POST'])
def get_robot_update(filtertype='ukf'):
    data = request.get_json()
    if not data or 'leftwheel' not in data or \
            'rightwheel' not in data or 'robots' not in data:
        abort(404)

    leftwheel = float(data['leftwheel'])
    rightwheel = float(data['rightwheel'])
    print(leftwheel, rightwheel)
    new_positions = []
    for robot in data['robots']:
        new_positions.append(movement_update(
            leftwheel, rightwheel, robot))

    return jsonify({"robots": new_positions})


if __name__ == '__main__':
    load_tutorials()
    app.run()
