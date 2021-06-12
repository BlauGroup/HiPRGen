import networkx as nx
from copy import deepcopy

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
    agraph.draw(path, format="pdf")


def visualize_molecules(mol_entries, folder):

    os.mkdir(folder)
    for index, molecule_entry in enumerate(mol_entries):
        visualize_molecule_entry(
            molecule_entry,
            folder + "/" + str(index) + ".pdf")



class ReportGenerator:

    def __init__(
            self,
            mol_entries,
            report_file_name,    # use full file path
            mol_pictures_folder, # use full file path
            rebuild_mol_pictures=True
    ):

        if rebuild_mol_pictures:
            visualize_molecules(mol_entries, mol_pictures_folder)

        self.mol_entries = mol_entries
        self.f = open(report_file_name, 'w')
        self.mol_pictures_folder = mol_pictures_folder

        # write in header
        f.write("\\documentclass{article}\n")
        f.write("\\usepackage{graphicx}\n")
        f.write("\\usepackage[margin=1cm]{geometry}\n")
        f.write("\\usepackage{amsmath}\n")
        f.write("\\pagenumbering{gobble}\n")
        f.write("\\begin{document}\n")

    def finished(self):
        f.write("\\end{document}")
        f.close()

    def emit_molecule(self, species_index):
        self.f.write(str(species_index) + "\n")
        self.f.write(
            "\\raisebox{-.5\\height}{"
            + "\\includegraphics[scale=0.2]{"
            + self.mol_pictures_folder
            + str(species_index)
            + ".pdf}}\n"
        )


    def emit_reaction(self, reaction):
        self.f.write("$$\n")
        first = True

        for reactant_index in reaction["reactants"]:
            if first:
                first = False
            else:
                self.f.write("+\n")

            self.emit_molecule(reactant_index)

        self.f.write(
            "\\xrightarrow{" + ("%.2f" % reaction["dG"]) + "}\n")

        first = True
        for product_index in reaction["products"]:
            if first:
                first = False
            else:
                self.f.write("+\n")

            self.emit_molecule(product_index)

        self.f.write("$$")
        self.f.write("\n\n\n")
