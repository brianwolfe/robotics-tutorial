from app import db
from sqlalchemy.dialects.postgresql import JSON

class SimulationResult(db.Model):
    __tablename__ = 'simresults'

    id = db.Column(db.Integer, primary_key=True)
    sim_results = db.Column(JSON)

    def __init__(self, sim_results):
        self.sim_results = sim_results

    def __repr__(self):
        return '<id {}>'.format(self.id)

