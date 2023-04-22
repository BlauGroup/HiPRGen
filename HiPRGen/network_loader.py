import sqlite3
import pickle
import numpy as np

"""
class for dynamically loading a reaction network
"""

sql_get_reaction = """
    SELECT * FROM reactions WHERE reaction_id = ?;
"""

sql_get_reaction_range = """
    SELECT * FROM reactions WHERE ? <= reaction_id AND reaction_id < ?;
"""

sql_get_redox = """
    SELECT * FROM reactions WHERE is_redox = 1;
"""

def sql_get_coord(metal_id):
    return "SELECT * FROM reactions WHERE (number_of_reactants=2 AND number_of_products=1 AND (reactant_1={0} OR reactant_2={0})) ORDER BY dG DESC;".format(metal_id)

def sql_get_decoord(metal_id):
    return "SELECT * FROM reactions WHERE (number_of_reactants=1 AND number_of_products=2 AND (product_1={0} OR product_2={0})) ORDER BY dG DESC;".format(metal_id)


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

        self.reactions = {}
        self.trajectories = {}
        self.initial_state_dict = {}
        self.initial_state_array = {}

    def get_all_redox_reactions(self):
        redox_reactions = []
        cur = self.rn_con.cursor()
        for res in cur.execute(sql_get_redox):
            reaction = {}
            reaction['number_of_reactants'] = res[1]
            reaction['number_of_products'] = res[2]
            reaction['reactants'] = res[3:5]
            reaction['products'] = res[5:7]
            reaction['rate'] = res[7]
            reaction['dG'] = res[8]
            reaction['dG_barrier'] = res[9]
            redox_reactions.append(reaction)

        return redox_reactions


    def get_all_coordination_reactions(self, metal_id):
        coordination_reactions = []
        cur = self.rn_con.cursor()
        for res in cur.execute(sql_get_coord(metal_id)):
            reaction = {}
            reaction['number_of_reactants'] = res[1]
            reaction['number_of_products'] = res[2]
            reaction['reactants'] = res[3:5]
            reaction['products'] = res[5:7]
            reaction['rate'] = res[7]
            reaction['dG'] = res[8]
            reaction['dG_barrier'] = res[9]
            coordination_reactions.append(reaction)

        return coordination_reactions


    def get_all_decoordination_reactions(self, metal_id):
        decoordination_reactions = []
        cur = self.rn_con.cursor()
        for res in cur.execute(sql_get_decoord(metal_id)):
            reaction = {}
            reaction['number_of_reactants'] = res[1]
            reaction['number_of_products'] = res[2]
            reaction['reactants'] = res[3:5]
            reaction['products'] = res[5:7]
            reaction['rate'] = res[7]
            reaction['dG'] = res[8]
            reaction['dG_barrier'] = res[9]
            decoordination_reactions.append(reaction)

        return decoordination_reactions


    def get_reactions_in_range(self, lower_bound, upper_bound):
        """
        get range of reactions from database but don't cache them
        """
        cur = self.rn_con.cursor()
        for res in cur.execute(sql_get_reaction_range,
                               (lower_bound, upper_bound)):
            reaction = {}
            reaction['reaction_id'] = res[0]
            reaction['number_of_reactants'] = res[1]
            reaction['number_of_products'] = res[2]
            reaction['reactants'] = res[3:5]
            reaction['products'] = res[5:7]
            reaction['rate'] = res[7]
            reaction['dG'] = res[8]
            reaction['dG_barrier'] = res[9]
            yield reaction


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
                cur.execute(sql_get_reaction, (reaction_index,))
            )[0]
            reaction = {}
            reaction['number_of_reactants'] = res[1]
            reaction['number_of_products'] = res[2]
            reaction['reactants'] = res[3:5]
            reaction['products'] = res[5:7]
            reaction['rate'] = res[7]
            reaction['dG'] = res[8]
            reaction['dG_barrier'] = res[9]
            self.reactions[reaction_index] = reaction
            return reaction


    def load_initial_state_and_trajectories(self):

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

        if self.initial_state_dict == {}:
            self.initial_state_dict = initial_state_dict
        else:
            for i in range(self.number_of_species):
                if initial_state_dict[i] > self.initial_state_dict[i]:
                    self.initial_state_dict[i] = initial_state_dict[i]

        for row in cur.execute(sql_get_trajectory):
            seed = row[0]
            step = row[1]
            reaction_id = row[2]
            time = row[3]

            if seed not in self.trajectories:
                self.trajectories[seed] = {}
                self.initial_state_array[seed] = initial_state_array

            self.trajectories[seed][step] = (reaction_id, time)
                    


    def set_initial_state_db(self, initial_state_database):
        # NOTE: switching to a new initial state database and loading in trajectory
        # info from it will only work if the new database has different seeds!
        self.initial_state_con = sqlite3.connect(initial_state_database)
