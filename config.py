import os

class Config(object):
    DEBUG = False
    TESTING = False
    DEVELOPMENT = False
    CSRF_ENABLED = True
    SECRET_KEY = 'fakesecret'
    SQLALCHEMY_DATABASE_URI = os.environ['DATABASE_URL']

    # Sending static assets via flask. Yes, this is a bad idea,
    # but I don't have S3 set up correctly yet. Will serve from there
    # after I get things set up. This allows me to develop without an
    # internet connection.
    BOWER_PREFIX = '/static/bower_components'
    D3_URL = '/'.join([BOWER_PREFIX, 'd3', 'd3.js'])
    FOUNDATION_URL = '/'.join([BOWER_PREFIX, 'foundation', 'js', 'foundation.js'])
    JQUERY_URL = '/'.join([BOWER_PREFIX, 'jquery', 'dist', 'jquery.js'])
    HANDLEBARS_URL = '/'.join([BOWER_PREFIX, 'handlebars', 'handlebars.js'])
    MODERNIZR_URL = '/'.join([BOWER_PREFIX, 'modernizr', 'modernizr.js'])
    UNDERSCORE_URL = '/'.join([BOWER_PREFIX, 'underscore', 'underscore.js'])
    NVD3_URL = '/'.join([BOWER_PREFIX, 'nvd3', 'nv.d3.js'])
    MATHJAX_URL = '/'.join([BOWER_PREFIX, 'MathJax', 'MathJax.js'])

    STATIC_URL = '/'.join(['/static'])
    APP_URL = '/'.join([STATIC_URL, 'js', 'app.js'])

    APP_STYLESHEET_URL = '/'.join([STATIC_URL, 'stylesheets', 'app.css'])

class ProductionConfig(Config):
    DEBUG = False


class StagingConfig(Config):
    DEVELOPMENT = True


class DevelopmentConfig(Config):
    DEVELOPMENT = True
    DEBUG = True


class TestingConfig(Config):
    TESTING = True
