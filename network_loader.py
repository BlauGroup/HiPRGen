import sqlite3
import pickle


"""
class for dynamically loading a reaction network
"""

sql_get_reaction = """
    SELECT * FROM reactions WHERE reaction_id = ?;
"""

sql_get_trajectory = """
    SELECT * FROM trajectories;
"""

sql_get_initial_state = """
    SELECT * FROM initial_state;
"""

class NetworkLoader:

    def __init__(
            self,
            network_database,
            mol_entries):


        self.con = sqlite3.connect(network_database)
        self.mol_entries = mol_entries

        self.reactions = {}


    def index_to_reaction(self, reaction_index):

        """
        this method gets called a lot, so we cache the reactions to
        minimize database interaction
        """

        if reaction_index in self.reactions:
            return self.reactions[reaction_index]

        else:
            print("fetching data for reaction", reaction_index)
            cur = self.con.cursor()
            res = list(
                cur.execute(sql_get_reaction, (int(reaction_index),))
            )[0]
            reaction = {}
            reaction['number_of_reactants'] = res[1]
            reaction['number_of_products'] = res[2]
            reaction['reactants'] = res[3:5]
            reaction['products'] = res[5:7]
            reaction['rate'] = res[7]
            reaction['dG'] = res[8]
            self.reactions[reaction_index] = reaction
            return reaction

    def load_trajectories(self):

        cur = self.con.cursor()

        # trajectories[seed][step] = (reaction_id, time)
        trajectories = {}
        for row in cur.execute(sql_get_trajectory):
            seed = row[0]
            step = row[1]
            reaction_id = row[2]
            time = row[3]

            if seed not in trajectories:
                trajectories[seed] = {}

            trajectories[seed][step] = (reaction_id, time)

        self.trajectories = trajectories


    def load_initial_state(self):

        cur = self.con.cursor()
        initial_state = {}

        for row in cur.execute(sql_get_initial_state):
            initial_state[row[0]] = row[1]

        self.initial_state = initial_state
