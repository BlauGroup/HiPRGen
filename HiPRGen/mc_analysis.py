from HiPRGen.report_generator import *
from HiPRGen.network_loader import *
from HiPRGen.constants import *
import math



def default_cost(free_energy):
    return math.exp(min(10.0, free_energy) / (ROOM_TEMP * KB)) + 1


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
                report_generator.emit_reaction(reaction)

            report_generator.emit_newpage()
            count += 1
            if count > number_of_pathways:
                break

        report_generator.finished()
