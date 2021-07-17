import networkx as nx
from copy import deepcopy
import os
from pathlib import Path

def visualize_molecule_entry(molecule_entry, path):
    """
    visualize a molecule using graphviz and
    output the resulting pdf to path
    """

    atom_colors = {"H": "gray",
                   "C": "black",
                   "Li": "#cc80ff",
                   "B": "#ffb5b5",
                   "N": "#3050f8",
                   "O": "#ff0d0d",
                   "F": "#90e050",
                   "Mg": "#8aff00",
                   "P": "#ff8000",
                   "S": "#ffff30",
                   "Cl": "#1ff01f"}

    graph = deepcopy(molecule_entry.graph)


    nx.set_node_attributes(graph, "filled", "style")
    nx.set_node_attributes(graph, "circle", "shape")
    nx.set_node_attributes(graph, "0.2", "width")
    nx.set_node_attributes(graph, "8.0", "fontsize")
    nx.set_node_attributes(graph, "white", "fontcolor")
    nx.set_node_attributes(graph, "true", "fixedsize")


    nx.set_node_attributes(
        graph,
        dict(enumerate([atom_colors[a]
                        for a in molecule_entry.species])),
        "color"
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
            mol_pictures_folder_name='mol_pictures',
            rebuild_mol_pictures=True
    ):
        self.report_file_path = Path(report_file_path)
        self.mol_pictures_folder_name = mol_pictures_folder_name
        self.mol_pictures_folder = self.report_file_path.parent.joinpath(
            mol_pictures_folder_name)


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
        self.f.write("\\setlength\\parindent{0pt}\n")

    def finished(self):
        self.f.write("\\end{document}")
        self.f.close()

    def emit_molecule(self, species_index, include_index=True):
        if include_index:
            self.f.write(str(species_index) + "\n")

        self.f.write(
            "\\raisebox{-.5\\height}{"
            + "\\includegraphics[scale=0.2]{"
            + self.mol_pictures_folder_name + '/'
            + str(species_index)
            + ".pdf}}\n"
        )

    def emit_newline(self):
        self.f.write(
            "\n\\vspace{1cm}\n")

    def emit_newpage(self):
        self.f.write("\\newpage\n\n\n")

    def emit_verbatim(self, s):
        self.f.write('\\begin{verbatim}\n')
        self.f.write(s)
        self.f.write('\n')
        self.f.write('\\end{verbatim}\n')

    def emit_text(self,s):
        self.f.write('\n\n' + s + '\n\n')

    def emit_initial_state(self, initial_state):
        self.emit_text("initial state:")
        for species_id in initial_state:
            num = initial_state[species_id]
            if num > 0:
                self.emit_text(str(num) + " molecules of")
                self.emit_molecule(species_id)
                self.emit_newline()


    def emit_reaction(self, reaction, label=None):
        reactants_filtered = [i for i in reaction['reactants']
                              if i != -1]

        products_filtered = [i for i in reaction['products']
                             if i != -1]

        self.f.write("$$\n")
        if label is not None:
            self.f.write(label + ":  \n")

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

    def emit_bond_breakage(self, reaction):
        if 'reactant_bonds_broken' in reaction:
            self.f.write("reactant bonds broken:")
            for bond in reaction['reactant_bonds_broken']:
                self.emit_verbatim(str(bond))

        if 'product_bonds_broken' in reaction:
            self.f.write("product bonds broken:")
            for bond in reaction['product_bonds_broken']:
                self.emit_verbatim(str(bond))

        self.f.write("\n\n\n")
