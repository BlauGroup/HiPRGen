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

class NetworkLoader:

    def __init__(
            network_database,
            mol_entries_pickle):


        self.con = sqlite3.connect(network_database)

        with open(mol_entries_pickle, 'rb') as f:
            self.mol_entries = pickle.load(f)


        self.reactions = {}


    def index_to_reactions(self, reaction_index):

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



