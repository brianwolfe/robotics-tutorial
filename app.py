from flask import Flask, render_template, abort, jsonify, request
from flask.ext.sqlalchemy import SQLAlchemy
import os
import os.path
import sys
# from urllib.parse import urljoin
from numpy import (array, concatenate, cos, diag, exp, float64, sin,
                   zeros)
from numpy.linalg import cholesky, inv

from numpy.random import multivariate_normal

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
        jsname = app.config['STATIC_URL'] + '/js/' + tutorial + '.js'
        return render_template(
            'tutorial.html',
            body=tutorial_pages[tutorial],
            jsname=jsname)
    else:
        abort(404)


# Sending static assets from Flask directly.
# I'll update this to send from CDN URLs and
# S3 URLs eventually, but this works for now.
# @app.route('/static/<path:filename>')
# def get_static_assets(filename):
#     return send_from_directory('/static/', filename)


ROBOT_RADIUS = 1
UKF_WEIGHT_ALPHA = 1.0
UKF_WEIGHT_BETA = 2.0
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
    wc_diag = diag(wc)

    # Expected measurements from sigmapoints
    Z_points = sigmapoints.dot(array([[1, 0], [0, 1], [0, 0]]))

    # Mean of expected measurements
    z_mean = wm.dot(Z_points)

    # Covariance of expected measurements
    S_cov = (Z_points - z_mean).T.dot(wc_diag).dot(Z_points - z_mean) + Q
    # Cross-covariance of position and expected measurements
    sigma_xz = (sigmapoints - mu).T.dot(wc_diag).dot(Z_points - z_mean)
    # Kalman gain
    K = sigma_xz.dot(inv(S_cov))

    # Updated mean
    mu_new = mu + K.dot(z - z_mean)
    # Updated covariance
    sigma_new = sigma - sigma_xz.dot(K.T)
    # Compute the updated sigmapoints
    sigmapoints_new, wm, wc = ukf_get_sigmapoints(mu_new, sigma_new)

    robot = {
        'mu': mu_new.reshape((-1, 1)).tolist(),
        'sigma': sigma_new.tolist(),
        'type': 'ukf',
        'previous_sigmapoints': sigmapoints.tolist(),
        'current_sigmapoints': sigmapoints_new.tolist()
        }

    return robot

def histogram_measurement(x, y, Q, robot):
    if 'hist' not in robot or 'x' not in robot or 'y' not in robot or \
            'heading' not in robot:
        abort(404)

    hist = array(robot['hist'], ndmin=2)
    xbins = robot['x']
    ybins = robot['y']
    hbins = robot['h']

    if hist.shape != (len(xbins), len(ybins), len(hbins)):
        abort(404)

    Qinv = inv(Q)

    z = array((x, y))

    for i, xm in enumerate(xbins):
        for j, ym in enumerate(ybins):
            mu = array((xm, ym))
            hist[i, j, k] *= exp(-(z - mu).dot(Qinv).dot(z - mu)/2)
    hist /= sum(sum(sum(hist)))

    return {
        'type': 'hist',
        'x': xbins,
        'y': ybins,
        'heading': hbins,
        'hist': hist.tolist()
        }

def particle_measurement(x, y, Q, robot):
    if 'particles' not in robot:
        abort(404)

    hist = array(robot['hist'])

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
    Q = array(data['Q'], dtype=float64, ndmin=2)

    if Q.shape != (2, 2):
        abort(404)

    new_robots = []
    for robot in data['robots']:
        if robot['type'] == 'ukf':
            nr = ukf_measurement(x, y, Q, robot)
            new_robots.append(nr)
        elif robot['type'] == 'hist':
            nr = histogram_measurement(x, y, Q, robot)
            new_robots.append(nr)
        elif robot['type'] == 'particle':
            nr = particle_measurement(x, y, Q, robot)
            new_robots.append(nr)

    return jsonify({'robots': new_robots})


def ukf_get_sigmapoints(mu, sigma):

    alpha = UKF_WEIGHT_ALPHA
    beta = UKF_WEIGHT_BETA
    kappa = UKF_WEIGHT_KAPPA

    n = sigma.shape[0]
    ukf_lambda = alpha**2 * (n + kappa) - n

    mu = mu.reshape((-1,))
    sigma_offsets = cholesky((n + ukf_lambda) * sigma).T

    sigmapoints = concatenate([
        [mu],
        mu - sigma_offsets,
        mu + sigma_offsets
        ])

    w0m = ukf_lambda / (n + ukf_lambda)
    w0c = ukf_lambda / (n + ukf_lambda) + (1 - alpha**2 + beta)
    wi = 1 / (2 * (n + ukf_lambda))
    wm = array([w0m] + [wi]*(sigmapoints.shape[0] - 1))
    wc = array([w0c] + [wi]*(sigmapoints.shape[0] - 1))

    return sigmapoints, wm, wc


def ukf_movement(left_wheel, right_wheel, R, robot):
    mean = array(robot['mu'], dtype=float64)
    cov = array(robot['sigma'], dtype=float64)

    dd = (left_wheel + right_wheel) / 2
    dt = (right_wheel - left_wheel) / ROBOT_RADIUS

    xi, wm, wc = ukf_get_sigmapoints(mean, cov)

    wc_diag = diag(wc)

    x1 = xi[:, 0]
    y1 = xi[:, 1]
    theta1 = xi[:, 2]

    x2 = (x1 + dd * cos(theta1)).reshape((-1, 1))
    y2 = (y1 + dd * sin(theta1)).reshape((-1, 1))
    theta2 = (theta1 + dt).reshape((-1, 1))

    xi_new = concatenate([x2, y2, theta2], axis=1)

    mu_new = wm.dot(xi_new)
    sigma_new = (xi_new - mu_new).T.dot(wc_diag).dot(xi_new - mu_new) + R

    return {'mu': mu_new.reshape((-1, 1)).tolist(),
            'sigma': sigma_new.tolist(),
            'type': 'ukf'}


def histogram_movement(leftwheel, rightwheel, R, robot):
    if 'hist' not in robot or 'x' not in robot or 'y' not in robot or \
            'heading' not in robot:
        abort(404)

    xbins = robot['x']
    ybins = robot['y']
    hbins = robot['heading']

    hist = robot['hist']

    newhist = zeros(xbins.shape[0], ybins.shape[0], hbins.shape[0])

    dd = (leftwheel + rightwheel) / 2
    dh = (rightwheel - leftwheel) / ROBOT_RADIUS
    Rinv = inv(R)

    def p(x1, y1, h1, x2, y2, h2):
        mu = array([
            x1 + dd*cos(h1),
            y1 + dd*cos(h1),
            h1 + dh
            ])
        prob = exp(-mu.dot(Rinv).dot(mu)/2)
        return prob

    for i1, x1 in enumerate(xbins):
        for j1, y1 in enumerate(ybins):
            for k1, h1 in enumerate(hbins):
                for i2, x2 in enumerate(xbins):
                    for j2, y2 in enumerate(ybins):
                        for k2, h2 in enumerate(hbins):
                            newhist[i2, j2, k2] += hist[i1, j1, k1] *\
                                p(x1, y1, h1, x2, y2, h2)

    # Equally spaced, so all bins are equal
    newhist /= sum(sum(sum(newhist)))
    return {
        'type': 'hist',
        'x': xbins,
        'y': ybins,
        'heading': hbins,
        'hist': newhist.tolist()
        }

def particle_movement(leftwheel, rightwheel, R, robot):
    if not 'particles' in robot:
        abort(404)

    particles = array(robot['particles'])

    dd = (leftwheel + rightwheel) / 2
    dh = (rightwheel - leftwheel) / ROBOT_RADIUS

    z = zeros((len(particles[0]),))

    noises = multivariate_normal(z, R)

    new_particles = zeros(particles.shape)
    for i, (particle, noise) in enumerate(zip(particles, noises)):
        x = particle[0]
        y = particle[1]
        h = particle[2]
        new_particles[i] = array([
            x + dd*cos(h),
            y + dd*cos(h),
            h + dh
            ]) + noises

    return {
            'type': 'particle',
            'particles': new_particles.tolist()
            }

# TODO split the robots out into different modules
# staying here for now
@app.route('/movementupdate/<path:filtertype>', methods=['POST'])
def movement_update_request(filtertype='ukf'):
    data = request.get_json()

    if not data or 'leftwheel' not in data or \
            'rightwheel' not in data or 'robots' not in data or \
            'R' not in data:
        abort(404)

    leftwheel = float(data['leftwheel'])
    rightwheel = float(data['rightwheel'])
    R = array(data['R'], dtype=float64)

    new_positions = []
    for robot in data['robots']:
        if robot['type'] is not None:
            if robot['type'] == 'ukf':
                new_positions.append(ukf_movement(
                    leftwheel, rightwheel, R, robot))
            elif robot['type'] == 'histogram':
                new_positions.append(histogram_movement(
                    leftwheel, rightwheel, R, robot))
            elif robot['type'] == 'particle':
                new_positions.append(particle_movement(
                    leftwheel, rightwheel, R, robot))

    return jsonify({"robots": new_positions})

load_tutorials()

if __name__ == '__main__':
    app.run()
