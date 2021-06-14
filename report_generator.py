import networkx as nx
from copy import deepcopy
import os
from pathlib import Path

def visualize_molecule_entry(molecule_entry, path):
    """
    visualize a molecule using graphviz and
    output the resulting pdf to path
    """

    atom_colors = {"O": "red",
                   "H": "gray",
                   "C": "black",
                   "Li": "purple"}

    graph = deepcopy(molecule_entry.graph).to_undirected()

    nx.set_node_attributes(graph, "", "label")
    nx.set_node_attributes(graph, "filled", "style")
    nx.set_node_attributes(graph, "circle", "shape")
    nx.set_node_attributes(graph, "0.2", "width")
    nx.set_node_attributes(
        graph,
        dict(enumerate([atom_colors[a]
                        for a in molecule_entry.species])),
        "color",
    )

    charge = molecule_entry.charge
    agraph = nx.nx_agraph.to_agraph(graph)
    if charge != 0:
        agraph.add_node(
            "charge",
            label=str(charge),
            fontsize="25.0",
            shape="box",
            color="gray",
            style="dashed, rounded",
        )

    agraph.layout()
    agraph.draw(path.as_posix(), format="pdf")


def visualize_molecules(mol_entries, folder):

    folder.mkdir()
    for index, molecule_entry in enumerate(mol_entries):
        visualize_molecule_entry(
            molecule_entry,
            folder.joinpath(str(index) + ".pdf"))



class ReportGenerator:

    def __init__(
            self,
            mol_entries,
            report_file_path,
            rebuild_mol_pictures=True
    ):
        self.report_file_path = Path(report_file_path)
        self.mol_pictures_folder = self.report_file_path.parent.joinpath(
            'mol_pictures')

        if rebuild_mol_pictures:
            visualize_molecules(mol_entries, self.mol_pictures_folder)

        self.mol_entries = mol_entries
        self.f = self.report_file_path.open(mode='w')


        # write in header
        self.f.write("\\documentclass{article}\n")
        self.f.write("\\usepackage{graphicx}\n")
        self.f.write("\\usepackage[margin=1cm]{geometry}\n")
        self.f.write("\\usepackage{amsmath}\n")
        self.f.write("\\pagenumbering{gobble}\n")
        self.f.write("\\begin{document}\n")

    def finished(self):
        self.f.write("\\end{document}")
        self.f.close()

    def emit_molecule(self, species_index):
        self.f.write(str(species_index) + "\n")
        self.f.write(
            "\\raisebox{-.5\\height}{"
            + "\\includegraphics[scale=0.2]{"
            + './mol_pictures/'
            + str(species_index)
            + ".pdf}}\n"
        )

    def emit_newline(self):
        self.f.write(
            "\n\\vspace{1cm}\n")

    def emit_verbatim(self, s):
        self.f.write('\\begin{verbatim}\n')
        self.f.write(s)
        self.f.write('\n')
        self.f.write('\\end{verbatim}\n')


    def emit_reaction(self, reaction):
        reactants_filtered = [i for i in reaction['reactants']
                              if i != -1]

        products_filtered = [i for i in reaction['products']
                             if i != -1]
        self.f.write("$$\n")
        first = True

        for reactant_index in reactants_filtered:
            if first:
                first = False
            else:
                self.f.write("+\n")

            self.emit_molecule(reactant_index)

        if 'dG' in reaction:
            self.f.write(
                "\\xrightarrow{" + ("%.2f" % reaction["dG"]) + "}\n")
        else:
            self.f.write(
                "\\xrightarrow{}\n")

        first = True
        for product_index in products_filtered:
            if first:
                first = False
            else:
                self.f.write("+\n")

            self.emit_molecule(product_index)

        self.f.write("$$")
        self.f.write("\n\n\n")
