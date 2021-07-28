from HiPRGen.report_generator import *
from HiPRGen.network_loader import *
from HiPRGen.constants import *
import math


def default_cost(free_energy):
    return math.exp(min(10.0, free_energy) / (ROOM_TEMP * KB)) + 1


def reaction_tally_report(
        network_loader,
        reaction_tally_report_path,
        cutoff=10):

    reaction_tally = {}
    for seed in network_loader.trajectories:
        for step in network_loader.trajectories[seed]:
            reaction_id = network_loader.trajectories[seed][step][0]

            if reaction_id in reaction_tally:
                reaction_tally[reaction_id] += 1
            else:
                reaction_tally[reaction_id] = 1


    report_generator = ReportGenerator(
        network_loader.mol_entries,
        reaction_tally_report_path,
        rebuild_mol_pictures=False)

    report_generator.emit_text("reaction tally report")

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
        report_generator.emit_molecule(i)
        report_generator.emit_newline()

    report_generator.finished()




class Pathfinding:

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


            for i in range(reaction['number_of_reactants']):
                reactant_id = reaction['reactants'][i]
                if self.network_loader.initial_state[reactant_id] == 0:
                    pathway = self.compute_pathway(reactant_id, trajectory) + pathway


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


    def generate_pathway_report(
            self,
            species_id,
            report_file_path,
            number_of_pathways=100,
            sort_by_frequency=True
    ):

        report_generator = ReportGenerator(
            self.network_loader.mol_entries,
            report_file_path,
            rebuild_mol_pictures=False)

        pathways = self.compute_pathways(species_id)


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
        report_generator.emit_initial_state(self.network_loader.initial_state)
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
                reaction = self.network_loader.index_to_reaction(reaction_index)
                report_generator.emit_reaction(reaction, label=str(reaction_index))

            report_generator.emit_newpage()
            count += 1
            if count > number_of_pathways:
                break

        report_generator.finished()

def sink_report(
        network_loader,
        sink_report_path
):


    consumed_dict = [[0,{}] for i in
                     range(network_loader.number_of_species)]

    produced_dict = [[0,{}] for i in
                     range(network_loader.number_of_species)]



    for seed in network_loader.trajectories:
        for step in network_loader.trajectories[seed]:
            reaction_index = network_loader.trajectories[seed][step][0]

            reaction = network_loader.index_to_reaction(reaction_index)

            for i in range(reaction['number_of_reactants']):
                reactant_index = reaction['reactants'][i]
                consumed_dict[reactant_index][0] += 1
                consumed_dict[reactant_index][1][reaction_index] = True

            for j in range(reaction['number_of_products']):
                product_index = reaction['products'][j]
                produced_dict[product_index][0] += 1
                produced_dict[product_index][1][reaction_index] = True

    max_ratio = 1e10
    ratio_dict = [max_ratio] * network_loader.number_of_species

    for i in range(network_loader.number_of_species):
        if consumed_dict[i][0] != 0:
            ratio_dict[i] = produced_dict[i][0] / consumed_dict[i][0]


    sink_data = sorted(
        list(enumerate(zip(consumed_dict, produced_dict, ratio_dict))),
        key=lambda item: -item[1][2])



    report_generator = ReportGenerator(
        network_loader.mol_entries,
        sink_report_path,
        rebuild_mol_pictures=False)

    for species_index, (c,p,r) in sink_data:
        if c[0] + p[0] > 0:
            report_generator.emit_text("ratio: " + str(r))
            report_generator.emit_text("produced: " + str(p[0]))
            report_generator.emit_text(str(len(p[1])) + " distinct producing reactions")
            report_generator.emit_text("consumed: " + str(c[0]))
            report_generator.emit_text(str(len(c[1])) + " distinct consuming reactions")
            report_generator.emit_molecule(species_index)
            report_generator.emit_newline()


    report_generator.finished()


def consumption_report(
        network_loader,
        species_index,
        consumption_report_path
):


    producing_reactions = {}
    consuming_reactions = {}

    for seed in network_loader.trajectories:
        for step in network_loader.trajectories[seed]:
            reaction_index = network_loader.trajectories[seed][step][0]

            reaction = network_loader.index_to_reaction(reaction_index)

            for i in range(reaction['number_of_reactants']):
                reactant_index = reaction['reactants'][i]
                if reactant_index == species_index:
                    if reaction_index not in consuming_reactions:
                        consuming_reactions[reaction_index] = 1
                    else:
                        consuming_reactions[reaction_index] += 1



            for j in range(reaction['number_of_products']):
                product_index = reaction['products'][j]
                if product_index == species_index:
                    if reaction_index not in producing_reactions:
                        producing_reactions[reaction_index] = 1
                    else:
                        producing_reactions[reaction_index] += 1


    report_generator = ReportGenerator(
        network_loader.mol_entries,
        consumption_report_path,
        rebuild_mol_pictures=False)


    report_generator.emit_text("consuming reactions:")
    for (reaction_index, number) in sorted(
            consuming_reactions.items(),
            key=lambda item: -item[1]):

        reaction = network_loader.index_to_reaction(reaction_index)
        report_generator.emit_text(str(number) + " occourances:")
        report_generator.emit_reaction(reaction)


    report_generator.emit_text("producing reactions:")
    for (reaction_index, number) in sorted(
            producing_reactions.items(),
            key=lambda item: -item[1]):

        reaction = network_loader.index_to_reaction(reaction_index)
        report_generator.emit_text(str(number) + " occourances:")
        report_generator.emit_reaction(reaction)


    report_generator.finished()
