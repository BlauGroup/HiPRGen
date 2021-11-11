from HiPRGen.report_generator import ReportGenerator
from HiPRGen.network_renderer import Renderer
from HiPRGen.network_loader import NetworkLoader
from HiPRGen.constants import ROOM_TEMP, KB
from HiPRGen.reaction_questions import marcus_barrier
from monty.serialization import dumpfn
import math
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from itertools import chain

def default_cost(free_energy):
    return math.exp(min(10.0, free_energy) / (ROOM_TEMP * KB)) + 1


def render_reactions_which_fired(network_loader, path):
    renderer = Renderer()
    reactions_which_fired = set()
    species_which_formed = set()
    for seed in network_loader.trajectories:
        for step in network_loader.trajectories[seed]:
            reaction_id = network_loader.trajectories[seed][step][0]
            reaction = network_loader.index_to_reaction(reaction_id)
            reactions_which_fired.add(reaction_id)

            for i in range(reaction['number_of_reactants']):
                reactant_id = reaction['reactants'][i]
                species_which_formed.add(reactant_id)

            for j in range(reaction['number_of_products']):
                product_id = reaction['products'][j]
                species_which_formed.add(product_id)


    for species_id in species_which_formed:
        renderer.new_node(species_id)

    for reaction_id in reactions_which_fired:
        reaction = network_loader.index_to_reaction(reaction_id)
        for i in range(reaction['number_of_reactants']):
            for j in range(reaction['number_of_products']):
                reactant_id = reaction['reactants'][i]
                product_id = reaction['products'][j]
                renderer.draw_edge(reactant_id, product_id)


    for species_id in species_which_formed:
        renderer.draw_node(species_id)

    renderer.render(path)


def render_top_trajectories(pathfinding, threshold, output_path):
    pass

def redox_report(
        network_loader,
        redox_report_path,
        params
):

    redox_reactions = network_loader.get_all_redox_reactions()

    for r in redox_reactions:
        marcus_barrier(r, network_loader.mol_entries, params)
        r['dG_barrier'] = r['marcus_barrier']

    redox_reactions = sorted(redox_reactions, key=lambda r: r['marcus_barrier'])

    report_generator = ReportGenerator(
        network_loader.mol_entries,
        redox_report_path,
        rebuild_mol_pictures=False)

    report_generator.emit_text("redox report")
    report_generator.emit_text("marcus barrier appears below the arrow")
    for reaction in redox_reactions:
        report_generator.emit_reaction(reaction)

    report_generator.finished()


def reaction_tally_report(
        network_loader,
        reaction_tally_report_path,
        cutoff=10):

    reaction_tally = {}
    species_set = set()
    for seed in network_loader.trajectories:
        for step in network_loader.trajectories[seed]:
            reaction_id = network_loader.trajectories[seed][step][0]

            reaction = network_loader.index_to_reaction(reaction_id)

            if reaction_id in reaction_tally:
                reaction_tally[reaction_id] += 1
            else:
                reaction_tally[reaction_id] = 1

            for i in range(reaction['number_of_reactants']):
                reactant_id = reaction['reactants'][i]
                species_set.add(reactant_id)

            for j in range(reaction['number_of_products']):
                product_id = reaction['products'][j]
                species_set.add(product_id)



    report_generator = ReportGenerator(
        network_loader.mol_entries,
        reaction_tally_report_path,
        rebuild_mol_pictures=False)

    report_generator.emit_text("reaction tally report")
    report_generator.emit_text(
        "total number of reactions: " +
        str(network_loader.number_of_reactions))
    report_generator.emit_text(
        "number of reactions which fired: " +
        str(len(reaction_tally.keys())))
    report_generator.emit_text(
        "total number of species: " +
        str(network_loader.number_of_species))
    report_generator.emit_text(
        "number of species observed: " +
        str(len(species_set)))


    for (reaction_index, number) in sorted(
            reaction_tally.items(), key=lambda pair: -pair[1]):
        if number > cutoff:
            report_generator.emit_text(str(number) + " occourances of:")
            report_generator.emit_reaction(
                network_loader.index_to_reaction(reaction_index),
                label=str(reaction_index)
            )
            report_generator.emit_newline()

    report_generator.finished()

def species_report(network_loader, species_report_path):
    report_generator = ReportGenerator(
        network_loader.mol_entries,
        species_report_path,
        rebuild_mol_pictures=False)

    report_generator.emit_text("species report")
    for i in range(network_loader.number_of_species):
        mol = network_loader.mol_entries[i]
        report_generator.emit_text(str(mol.entry_id))
        report_generator.emit_text(
            "formula: " + mol.formula)

        report_generator.emit_molecule(i)
        report_generator.emit_newline()

    report_generator.finished()


class Pathfinding:
    """
    Given a chemical system, we are interested in exploring the reaction
    pathways which produce species of interest. Reaction networks are a
    useful tool for approaching the problem because they can be simulated
    efficiently, even when the network has hundreds of millions of
    reactions. Unfortunately, since we collapse all spacial aspects of the
    system, identical molecules become indistinguishable (in reality,
    identical molecules can be distinguished by their locations in
    space). This creates the following problem: Suppose we are interested
    in the production of species G from A and have the following
    simulation trajectory:

      A -> Z + F
      F -> X
      A -> B + F
      A -> C + H
      A -> D + H
      C -> E
      D -> E
      E + F -> G

    It is impossible to decide between the two pathways

      A -> B + F    A -> C + H        C -> E        E + F -> G
      A -> B + F    A -> D + H        D -> E        E + F -> G

    If our model had a spacial aspect, we would be able to trace the
    specific E used back to either a C or a D. Fundamentally, this
    ambiguity is caused by sequence

      A -> C + H        C -> E         E -> D        D + H -> A

    which is called a deficiency loop. To avoid this problem, we extract a
    pathway which produces the target in the following way: Take the first
    reaction x which produces the target molecule. Then recursively, take
    the first reactions which produced the reactants of x. If a reactant
    is a starting molecule, then stop. Applying the procedure to the above
    sequence gives the pathway

    A -> C + H        C -> E         A -> Z + F        E + F -> G

    Intuitively, this procedure is producing pathways which don't take
    into account competition, since there is no guarantee that the first
    molecule which is produced isn't immediately consumed by some
    competing reaction. Since the focus of this software is monte carlo
    simulation with thermodynamic rates, this isn't a problem as all
    downhill reactions have the maximum possible rate constant, so
    reactions sans competition are also likely to occour.
    """

    def __init__(
            self,
            network_loader
    ):

        self.network_loader = network_loader
        self.pathways = {}


    def compute_pathway(
            self,
            species_id,
            trajectory):

        pathway = []
        target_produced = False

        for step in trajectory:
            reaction_id = trajectory[step][0]
            reaction = self.network_loader.index_to_reaction(reaction_id)
            if species_id in reaction['products']:
                target_produced = True
                break

        if target_produced:
            pathway.append(reaction_id)


            prefixes = []
            for i in range(reaction['number_of_reactants']):
                reactant_id = reaction['reactants'][i]
                if self.network_loader.initial_state_dict[reactant_id] == 0:
                    prefix = self.compute_pathway(reactant_id, trajectory)
                    prefixes.append(prefix)

                    # early return if one of the final reaction in one
                    # of the prefixes gives both of the reactants
                    prefix_final_reaction = self.network_loader.index_to_reaction(
                        prefix[-1])

                    if (sorted(reaction['reactants']) ==
                        sorted(prefix_final_reaction['products'])):
                        return prefix + pathway

            for prefix in prefixes:
                pathway = prefix + pathway


        return pathway



    def compute_pathways(self, species_id):

        if species_id not in self.pathways:
            reaction_pathway_list = []
            for seed in self.network_loader.trajectories:
                pathway = self.compute_pathway(
                    species_id,
                    self.network_loader.trajectories[seed]
                )

                if len(pathway) > 0:
                    reaction_pathway_list.append(pathway)

            self.pathways[species_id] = self.collect_duplicate_pathways(
                reaction_pathway_list)

        return self.pathways[species_id]


    def collect_duplicate_pathways(
        self, pathways
    ):
        pathway_dict = {}
        for pathway in pathways:
            key = frozenset(pathway)
            if key in pathway_dict:
                pathway_dict[key]["frequency"] += 1
            else:
                path_weight = self.compute_path_weight(pathway)
                pathway_dict[key] = {
                    "pathway": pathway,
                    "frequency": 1,
                    "weight": path_weight,
                }

        return pathway_dict


    def compute_path_weight(self, pathway):
        weight = 0.0
        for reaction_index in pathway:
            reaction = self.network_loader.index_to_reaction(reaction_index)
            weight += default_cost(reaction["dG"])
        return weight

def export_pathways_to_json(pathfinding, species_id, path):
    pathways = pathfinding.compute_pathways(species_id)
    reactions = {}
    for pathway in pathways:
        for reaction_id in pathway:
            db_reaction = pathfinding.network_loader.index_to_reaction(reaction_id)
            json_reactants = [ pathfinding.network_loader.mol_entries[i].entry_id
                               for i in db_reaction['reactants'] if i != -1]
            json_products = [ pathfinding.network_loader.mol_entries[i].entry_id
                               for i in db_reaction['products'] if i != -1]
            json_reaction = {
                'reactants' : json_reactants,
                'products' : json_products
                }
            reactions[reaction_id] = json_reaction

    dumpfn({
        'pathways' : list(pathways.values()),
        'reactions' : reactions}, path)



def generate_pathway_report(
        pathfinding,
        species_id,
        report_file_path,
        number_of_pathways=100,
        sort_by_frequency=False
):

    report_generator = ReportGenerator(
        pathfinding.network_loader.mol_entries,
        report_file_path,
        rebuild_mol_pictures=False)

    pathways = pathfinding.compute_pathways(species_id)


    report_generator.emit_text("pathway report for")
    report_generator.emit_molecule(species_id)

    if sort_by_frequency:
        report_generator.emit_text(
            "top " +
            str(number_of_pathways) +
            " pathways sorted by frequency")

    else:
        report_generator.emit_text(
            "top " +
            str(number_of_pathways) +
            " pathways sorted by cost")

    report_generator.emit_newline()
    report_generator.emit_initial_state(pathfinding.network_loader.initial_state_dict)
    report_generator.emit_newpage()

    if sort_by_frequency:
        def sort_function(item):
            return -item[1]["frequency"]

    else:
        def sort_function(item):
            return item[1]["weight"]


    count = 1
    for _, unique_pathway in sorted(pathways.items(), key=sort_function):

        frequency = unique_pathway["frequency"]
        weight = unique_pathway["weight"]

        report_generator.emit_text("pathway " + str(count))
        report_generator.emit_text("path weight: " + str(weight))
        report_generator.emit_text(str(frequency) + " occurrences:")

        for reaction_index in unique_pathway["pathway"]:
            reaction = pathfinding.network_loader.index_to_reaction(reaction_index)
            report_generator.emit_reaction(reaction, label=str(reaction_index))

        report_generator.emit_newpage()
        count += 1
        if count > number_of_pathways:
            break

    report_generator.finished()


class SimulationReplayer:
    """
    class for rerunning through all the simulations. This is
    relatively fast since we don't need to actually choose which
    reactions fire / update reaction propensities.
    """

    def __init__(self, network_loader):
        self.network_loader = network_loader


        self.compute_expected_final_state()
        self.compute_production_consumption_info()
        self.compute_sink_data()


    def compute_expected_final_state(self):
        self.expected_final_state = np.zeros(
            self.network_loader.number_of_species,
            dtype=int
        )

        for seed in self.network_loader.trajectories:
            state = np.copy(self.network_loader.initial_state_array)
            for step in self.network_loader.trajectories[seed]:
                reaction_index = self.network_loader.trajectories[seed][step][0]
                time = self.network_loader.trajectories[seed][step][1]
                reaction = self.network_loader.index_to_reaction(reaction_index)

                for i in range(reaction['number_of_reactants']):
                    reactant_index = reaction['reactants'][i]
                    state[reactant_index] -= 1

                for j in range(reaction['number_of_products']):
                    product_index = reaction['products'][j]
                    state[product_index] += 1

            self.expected_final_state += state

        self.expected_final_state = (
            self.expected_final_state / len(self.network_loader.trajectories))

    def compute_production_consumption_info(self):
        self.consuming_reactions = {}
        self.producing_reactions = {}

        for i in range(self.network_loader.number_of_species):
            self.consuming_reactions[i] = {}
            self.producing_reactions[i] = {}

        for seed in self.network_loader.trajectories:
            for step in self.network_loader.trajectories[seed]:
                reaction_index = self.network_loader.trajectories[seed][step][0]
                time = self.network_loader.trajectories[seed][step][1]
                reaction = self.network_loader.index_to_reaction(reaction_index)

                for i in range(reaction['number_of_reactants']):
                    reactant_index = reaction['reactants'][i]
                    if reaction_index not in self.consuming_reactions[reactant_index]:
                        self.consuming_reactions[reactant_index][reaction_index] = 1
                    else:
                        self.consuming_reactions[reactant_index][reaction_index] += 1


                for j in range(reaction['number_of_products']):
                    product_index = reaction['products'][j]
                    if reaction_index not in self.producing_reactions[product_index]:
                        self.producing_reactions[product_index][reaction_index] = 1
                    else:
                        self.producing_reactions[product_index][reaction_index] += 1

    def compute_state_time_series(self, seed):
        state_dimension_size = len(self.network_loader.initial_state_array)
        step_dimension_size = len(self.network_loader.trajectories[seed])
        time_series = np.zeros(
            (step_dimension_size, state_dimension_size),
            dtype=int)

        state = np.copy(self.network_loader.initial_state_array)
        for step in self.network_loader.trajectories[seed]:
            reaction_index = self.network_loader.trajectories[seed][step][0]
            time = self.network_loader.trajectories[seed][step][1]
            reaction = self.network_loader.index_to_reaction(reaction_index)

            for i in range(reaction['number_of_reactants']):
                reactant_index = reaction['reactants'][i]
                state[reactant_index] -= 1

            for j in range(reaction['number_of_products']):
                product_index = reaction['products'][j]
                state[product_index] += 1

            time_series[step] = state

        return time_series


    def time_series_graph(
            self,
            seeds,
            species_of_interest,
            path,
            colors = list(mcolors.TABLEAU_COLORS.values()),
            styles = ['solid', 'dotted', 'dashed', 'dashdot'],
            internal_index_labels=True
    ):


        max_trajectory_length = 0
        for seed in seeds:
            max_trajectory_length = max(
                max_trajectory_length,
                len(self.network_loader.trajectories[seed]))

        total_time_series = np.zeros(
            (max_trajectory_length, self.network_loader.number_of_species),
            dtype=int
            )

        for seed in seeds:
            total_time_series += pad_time_series(
                self.compute_state_time_series(seed),
                max_trajectory_length)

        total_time_series = total_time_series / len(seeds)

        background_species = set()
        for index in range(self.network_loader.number_of_species):
            for step in range(max_trajectory_length):
                if (total_time_series[step,index] > 0.1 and
                    index not in species_of_interest):
                    background_species.add(index)


        line_dict = {}
        i = 0
        for species_index in species_of_interest:
            r = i % len(colors)
            q = i // len(colors)
            line_dict[species_index] = (colors[r], styles[q])
            i += 1


        fig, (ax0, ax1, ax2) = plt.subplots(
            3, 1,
            figsize=(5,10),
            gridspec_kw={'height_ratios':[2,2,1]})

        y_max = 0
        for step in range(total_time_series.shape[0]):
            for species_index in species_of_interest:
                y_max = max(y_max, total_time_series[step,species_index])

        ax0.set_xlim([0,total_time_series.shape[0]])
        ax0.set_ylim([0,y_max+1])

        ax1.set_xlim([0,total_time_series.shape[0]])
        ax1.set_ylim([0,(y_max+1)/10])


        ticks = np.arange(0, total_time_series.shape[0])
        for i, species_index in enumerate(background_species):
            ax0.plot(ticks,
                     total_time_series[:, species_index],
                     color=mcolors.hsv_to_rgb((0,0,0.9))
                     )

            ax1.plot(ticks,
                     total_time_series[:, species_index],
                     color=mcolors.hsv_to_rgb((0,0,0.9))
                     )



        for species_index in species_of_interest:


            ax0.plot(ticks,
                     total_time_series[:, species_index],
                     color=line_dict[species_index][0],
                     linestyle=line_dict[species_index][1]
                     )

            ax1.plot(ticks,
                     total_time_series[:, species_index],
                     color=line_dict[species_index][0],
                     linestyle=line_dict[species_index][1]
                     )



        # creating a legend
        ax2.yaxis.set_visible(False)
        ax2.xaxis.set_visible(False)
        ax2.set_axis_off()
        ax2.set_xlim([0,1])
        ax2.set_ylim([0,1])


        i = 0
        for species_index in species_of_interest:
            if internal_index_labels:
                label = str(species_index)
            else:
                label = self.network_loader.mol_entries[species_index].entry_id

            r = i % 10
            q = i // 10
            pos_x = [q * 0.3 + 0, q * 0.3 + 0.05]
            pos_y = [r * 0.1 + 0.05, r * 0.1 + 0.05]
            text_pos_x = q * 0.3 + 0.06
            text_pos_y = r * 0.1 + 0.05

            ax2.plot(
                pos_x, pos_y,
                color=line_dict[species_index][0],
                linestyle=line_dict[species_index][1])

            ax2.text(text_pos_x, text_pos_y,
                     label,
                     fontsize=8,
                     horizontalalignment='left',
                     verticalalignment='center')
            i += 1

        fig.savefig(path, transparent=True)


    def compute_sink_data(self):
        max_ratio = 1e10
        self.sink_data = {}

        for i in range(self.network_loader.number_of_species):
            number_of_consuming_reactions = sum(
                self.consuming_reactions[i].values())
            number_of_distinct_consuming_reactions = len(
                self.consuming_reactions[i].keys())
            number_of_producing_reactions = sum(
                self.producing_reactions[i].values())
            number_of_distinct_producing_reactions = len(
                self.producing_reactions[i].keys())


            if number_of_consuming_reactions != 0:
                ratio = number_of_producing_reactions / number_of_consuming_reactions
            else:
                ratio = max_ratio

            expected_value = self.expected_final_state[i]

            self.sink_data[i] = {
                "species_index" : i,
                "number_of_consuming_reactions" : number_of_consuming_reactions,
                "number_of_distinct_consuming_reactions"
                : number_of_distinct_consuming_reactions,
                "number_of_producing_reactions" : number_of_producing_reactions,
                "number_of_distinct_producing_reactions"
                : number_of_distinct_producing_reactions,
                "ratio" : ratio,
                "expected_value" : expected_value
            }

        self.sinks = [
            i for i in range(self.network_loader.number_of_species)
            if self.sink_filter(i)]

    def sink_filter(self, species_index):
        sink_data = self.sink_data[species_index]

        number_of_consuming_reactions = sink_data["number_of_consuming_reactions"]
        number_of_distinct_consuming_reactions = sink_data[
            "number_of_distinct_consuming_reactions"]
        number_of_producing_reactions = sink_data["number_of_producing_reactions"]
        number_of_distinct_producing_reactions = sink_data[
            "number_of_distinct_producing_reactions"]
        ratio = sink_data["ratio"]
        expected_value = sink_data["expected_value"]

        mol = self.network_loader.mol_entries[species_index]
        if (number_of_consuming_reactions + number_of_producing_reactions > 0  and
            ratio > 1.5 and
            expected_value > 0.1 and
            mol.spin_multiplicity == 1):
            return True
        else:
            return False


def export_consumption_to_json(simulation_replayer, species_index, path):
    network_loader = simulation_replayer.network_loader
    producing_reactions = simulation_replayer.producing_reactions[species_index]
    consuming_reactions = simulation_replayer.consuming_reactions[species_index]
    reactions = {}
    for reaction_id in chain(producing_reactions, consuming_reactions) :
        db_reaction = network_loader.index_to_reaction(reaction_id)
        json_reactants = [ network_loader.mol_entries[i].entry_id
                           for i in db_reaction['reactants'] if i != -1]
        json_products = [ network_loader.mol_entries[i].entry_id
                           for i in db_reaction['products'] if i != -1]
        json_reaction = {
            'reactants' : json_reactants,
            'products' : json_products
            }
        reactions[reaction_id] = json_reaction

        dumpfn({
            'reactions' : reactions,
            'producing_reactions' : producing_reactions,
            'consuming_reactions' : consuming_reactions},
               path)

def pad_time_series(time_series, max_number_of_steps):
    num_steps = time_series.shape[0]
    state_size = time_series.shape[1]
    padded_time_series = np.zeros(
        (max_number_of_steps, state_size),
        dtype=int)

    for step in range(max_number_of_steps):
        if step < num_steps:
            padded_time_series[step] = time_series[step]
        else:
            padded_time_series[step] = time_series[num_steps-1]

    return padded_time_series



def export_sinks_to_json(simulation_replayer, path):
    sink_data = simulation_replayer.sink_data
    sink_data_json = {}
    for i in sink_data:
        mol = simulation_replayer.network_loader.mol_entries[i]
        sink_data_json[mol.entry_id] = sink_data[i]

    dumpfn(sink_data_json, path)

def export_species_report_to_json(network_loader, path):
    data = {}
    for i in range(network_loader.number_of_species):
        data[i] = network_loader.mol_entries[i].entry_id

    dumpfn(data, path)


def consumption_report(
        simulation_replayer,
        species_index,
        consumption_report_path
):

    sink_data = simulation_replayer.sink_data[species_index]

    producing_reactions = simulation_replayer.producing_reactions[species_index]
    consuming_reactions = simulation_replayer.consuming_reactions[species_index]

    report_generator = ReportGenerator(
        simulation_replayer.network_loader.mol_entries,
        consumption_report_path,
        rebuild_mol_pictures=False)


    report_generator.emit_text("P/C ratio: " +
                               str(sink_data["ratio"]))
    report_generator.emit_text("expected val: " +
                               str(sink_data["expected_value"]))
    report_generator.emit_newline()
    report_generator.emit_text("consuming reactions:")
    for (reaction_index, number) in sorted(
            consuming_reactions.items(),
            key=lambda item: -item[1]):

        reaction = simulation_replayer.network_loader.index_to_reaction(
            reaction_index)
        report_generator.emit_text(str(number) + " occourances:")
        report_generator.emit_reaction(reaction)


    report_generator.emit_text("producing reactions:")
    for (reaction_index, number) in sorted(
            producing_reactions.items(),
            key=lambda item: -item[1]):

        reaction = simulation_replayer.network_loader.index_to_reaction(
            reaction_index)

        report_generator.emit_text(str(number) + " occourances:")
        report_generator.emit_reaction(reaction)


    report_generator.finished()


def sink_report(
        simulation_replayer,
        sink_report_path
):


    report_generator = ReportGenerator(
        simulation_replayer.network_loader.mol_entries,
        sink_report_path,
        rebuild_mol_pictures=False)

    sinks_sorted = sorted(
        simulation_replayer.sinks,
        key = lambda i: -simulation_replayer.sink_data[i]["ratio"])

    for species_index in sinks_sorted:

        sink_data = simulation_replayer.sink_data[species_index]

        number_of_consuming_reactions = sink_data[
            "number_of_consuming_reactions"]
        number_of_distinct_consuming_reactions = sink_data[
            "number_of_distinct_consuming_reactions"]
        number_of_producing_reactions = sink_data[
            "number_of_producing_reactions"]
        number_of_distinct_producing_reactions = sink_data[
            "number_of_distinct_producing_reactions"]
        ratio = sink_data["ratio"]
        expected_value = sink_data["expected_value"]


        report_generator.emit_text("P/C ratio: " + str(ratio))
        report_generator.emit_text("expected val: " + str(expected_value))

        report_generator.emit_text(
            "produced: " + str(number_of_producing_reactions))

        report_generator.emit_text(
            str(number_of_distinct_producing_reactions) +
            " distinct producing reactions")

        report_generator.emit_text(
            "consumed: " + str(number_of_consuming_reactions))

        report_generator.emit_text(
            str(number_of_distinct_consuming_reactions) +
            " distinct consuming reactions")
        report_generator.emit_molecule(species_index)
        report_generator.emit_newline()

    report_generator.finished()
