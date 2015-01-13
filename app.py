from flask import Flask, render_template, abort, jsonify, request
from flask.ext.sqlalchemy import SQLAlchemy
import os
import os.path
import sys
from urllib.parse import urljoin
from numpy import (array, concatenate, float64, cos, sin, exp, mean,
    diag)
from numpy.linalg import cholesky, inv

cmd_folder = os.path.dirname(os.path.abspath(__file__))
if cmd_folder not in sys.path:
     sys.path.insert(0, cmd_folder)

app = Flask(__name__)
app.config.from_object(os.environ['APP_SETTINGS'])
db = SQLAlchemy(app)

tutorial_path = os.path.join(cmd_folder, 'tutorials')
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
                    open(os.path.join(tutorial_path, tutorial)).read()

@app.route('/', methods=['GET', 'POST'])
def index():
    return render_template('index.html')

@app.route('/tutorial/<path:tutorial>', methods=['GET'])
def tutorial(tutorial):
    if app.config['DEBUG']:
        load_tutorials()
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
UKF_WEIGHT_ALPHA = 1.0
UKF_WEIGHT_BETA = 0.0
UKF_WEIGHT_KAPPA = 1.0

def ukf_measurement(x, y, Q, robot):
    # Previous mean
    mu = array(robot['mu']).reshape((-1,))
    # Previous covariance
    sigma = array(robot['sigma'])
    # New measurement
    z = array([x, y])

    # Compute the sigmapoints based on the previous state
    sigmapoints, wm, wc = ukf_get_sigmapoints(mu, sigma)

    # Make diagonal version of weights to make algebra easier
    wm_diag, wc_diag = diag(wm), diag(wc)

    # Expected measurements from sigmapoints
    Z_points = sigmapoints.dot(array([[1, 0], [0, 1], [0, 0]]))

    # Mean of expected measurements
    z_mean = wm.dot(Z_points)
    # Covariance of expected measurements
    S_cov = (Z_points - z_mean).T.dot(wm_diag).dot(Z_points - z_mean) + Q
    # Cross-covariance of position and expected measurements
    sigma_xz = (sigmapoints - mu).T.dot(wc_diag).dot(Z_points - z_mean)
    # Kalman gain
    K = sigma_xz.dot(inv(S_cov))

    # Updated mean
    mu_new = mu + K.dot(z_mean - z)
    # Updated covariance
    sigma_new = sigma - K.dot(S_cov).dot(K.T)
    # Compute the updated sigmapoints
    sigmapoints_new, wm, wc = ukf_get_sigmapoints(mu_new, sigma_new)

    robot = {
            'mu': mu_new.reshape((-1, 1)).tolist(),
            'sigma': sigma_new.tolist(),
            'previous_sigmapoints': sigmapoints.tolist(),
            'current_sigmapoints': sigmapoints_new.tolist()
            }

    return robot

@app.route('/measurementupdate/<path:filtertype>', methods=['POST'])
def measurement_update_request(filtertype='ukf'):
    data = request.get_json()

    if not data or \
            'robots' not in data or \
            'x' not in data or \
            'y' not in data or \
            'Q' not in data:
        abort(404)
    x = float(data['x'])
    y = float(data['y'])
    Q = array(data['Q'], dtype=float64)

    new_robots = []
    for robot in data['robots']:
        nr = ukf_measurement(x, y, Q, robot)
        new_robots.append(nr)

    return jsonify({'robots': new_robots})

def ukf_get_sigmapoints(mu, sigma):
    mu.reshape((-1,))

    sigma_offsets = cholesky(sigma)
    # Scale offsets?

    sigmapoints = concatenate([
        [mu],
        mu - sigma_offsets,
        mu + sigma_offsets
        ])

    alpha = UKF_WEIGHT_ALPHA
    beta = UKF_WEIGHT_BETA
    kappa = UKF_WEIGHT_KAPPA

    n = sigma.shape[0]
    ukf_lambda = alpha**2 * (n + kappa) - n

    w0m = ukf_lambda / (n + ukf_lambda)
    w0c = ukf_lambda / (n + ukf_lambda) + (1 - alpha**2 + beta)
    wi = 1 / (2 * (n + ukf_lambda))
    wm = array([w0m] + [wi]*(sigmapoints.shape[0] - 1))
    wc = array([w0c] + [wi]*(sigmapoints.shape[0] - 1))

    return sigmapoints, wm, wc

def ukf_movement(left_wheel, right_wheel, robot):
    mean = array(robot['mu'], dtype=float64)
    cov = array(robot['sigma'], dtype=float64)

    dd = (left_wheel + right_wheel) / 2;
    dt = (right_wheel - left_wheel) / ROBOT_RADIUS

    sigmapoints, wm, wc = ukf_get_sigmapoints(mean, cov)

    wm_diag, wc_diag = diag(wm), diag(wc)

    x1 = mean[0]
    y1 = mean[1]
    theta1 = mean[2]

    x2 = x1 + dd * cos(theta1)
    y2 = y1 + dd * sin(theta1)
    theta2 = theta1 + dt

    cov[0,0] += 1
    cov[1,1] += 1
    cov[2,2] += 0.1
    mean2 = array([x2, y2, theta2]).reshape((3, 1))
    print(mean)
    print(mean2)

    return {'mu': mean2.tolist(), 'sigma': cov.tolist()}

# TODO split the robots out into different modules
# staying here for now
@app.route('/movementupdate/<path:filtertype>', methods=['POST'])
def movement_update_request(filtertype='ukf'):
    data = request.get_json()

    if not data or 'leftwheel' not in data or \
            'rightwheel' not in data or 'robots' not in data:
        abort(404)

    leftwheel = float(data['leftwheel'])
    rightwheel = float(data['rightwheel'])
    new_positions = []
    for robot in data['robots']:
        new_positions.append(ukf_movement(
            leftwheel, rightwheel, robot))

    return jsonify({"robots": new_positions})

load_tutorials()

if __name__ == '__main__':
    app.run()
