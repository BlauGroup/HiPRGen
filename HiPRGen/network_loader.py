import sqlite3
import pickle
import numpy as np

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
            mol_entries_pickle,
            initial_state_database=None
    ):


        self.rn_con = sqlite3.connect(network_database)

        with open(mol_entries_pickle, 'rb') as f:
            self.mol_entries = pickle.load(f)

        cur = self.rn_con.cursor()
        metadata = list(cur.execute("SELECT * FROM metadata"))[0]
        self.number_of_species = metadata[0]
        self.number_of_reactions = metadata[1]


        if initial_state_database:
            self.initial_state_con = sqlite3.connect(initial_state_database)
            self.load_trajectories()
            self.load_initial_state()

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
            cur = self.rn_con.cursor()
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

        cur = self.initial_state_con.cursor()

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

        cur = self.initial_state_con.cursor()
        initial_state_dict = {}

        for row in cur.execute(sql_get_initial_state):
            initial_state_dict[row[0]] = row[1]

        initial_state_array = np.zeros(
            self.number_of_species,
            dtype=int
        )

        for i in range(self.number_of_species):
            initial_state_array[i] = initial_state_dict[i]


        self.initial_state_dict = initial_state_dict
        self.initial_state_array = initial_state_array
