from HiPRGen.network_loader import *
import cairo
import math
import random


class QuadTreeNode:
    """
    origin is at top left so to agree with
    the cairo canvas coordinates.

    Notice that this is a recursive initializer. It creates
    1 + 4 + ... + 4^(depth) = O(4^(depth + 1)) QuadTreeNodes,
    so don't go too deep!
    """
    def __init__(self, depth, x_min, x_max, y_min, y_max):

        self.x_min = x_min
        self.x_max = x_max
        self.y_min = y_min
        self.y_max = y_max

        # you either have quads or data
        # if you have quads, you are non terminal
        # if you have data you are terminal
        self.quads = None
        self.data = []
        self.branch(depth)


    def branch(self, depth):
        """
        break node into 4 nodes.
        """

        if depth > 0:
            self.data = None
            self.quads = [
                None, # top left
                None, # top right
                None, # bottom left
                None  # bottom right
            ]

            self.x_mid = (self.x_min + self.x_max) / 2
            self.y_mid = (self.y_min + self.y_max) / 2

            # top left
            self.quads[0] = QuadTreeNode(
                depth - 1,
                self.x_min,
                self.x_mid,
                self.y_min,
                self.y_mid)

            # top right
            self.quads[1] = QuadTreeNode(
                depth - 1,
                self.x_mid,
                self.x_max,
                self.y_min,
                self.y_mid)

            # bottom left
            self.quads[2] = QuadTreeNode(
                depth - 1,
                self.x_min,
                self.x_mid,
                self.y_mid,
                self.y_max)

            # bottom right
            self.quads[3] = QuadTreeNode(
                depth - 1,
                self.x_mid,
                self.x_max,
                self.y_mid,
                self.y_max)

    def find_neighborhood(self,x,y):
        """
        find all nodes adjacent to our point.
        doesn't return the node actually containing our point.
        """
        node = self.find_node(x,y)
        x_diff = node.x_max - node.x_min
        y_diff = node.y_max - node.y_min
        maybe_adjacent_nodes = [
            self.find_node(x + x_diff, y),
            self.find_node(x - x_diff, y),
            self.find_node(x, y + y_diff),
            self.find_node(x, y - y_diff),
            self.find_node(x + x_diff, y + y_diff),
            self.find_node(x - x_diff, y + y_diff),
            self.find_node(x + x_diff, y - y_diff),
            self.find_node(x - x_diff, y - y_diff)
        ]

        adjacent_nodes = [n for n in maybe_adjacent_nodes if n is not None]
        return adjacent_nodes

    def find_node(self, x, y):
        """
        find the terminal node so that
        x_min <= x < x_max
        y_min <= y < y_max
        return None if there is no node.
        Note: this gives the wrong answer if called from a terminal node.
        """
        if self.quads is not None:
            for quad in self.quads:
                if (quad.x_min <= x < quad.x_max and
                    quad.y_min <= y < quad.y_max):
                    return quad.find_node(x,y)

            return None

        else:
            return self

    def __str__(self):
        return (
            "x : [" + str(self.x_min) + ", " + str(self.x_max) + ")  " +
            "y : [" + str(self.y_min) + ", " + str(self.y_max) + ")"
        )

    def __repr__(self):
        return self.__str__()


class RepulsiveSampler:
    def __init__(self,
                 rejection_radius,
                 x_min,
                 x_max,
                 y_min,
                 y_max,
                 global_mask, # reject a sample if global mask returns false
                 quad_tree_depth=7,
                 seed=42,
                 ):

        self.quad_tree = QuadTreeNode(quad_tree_depth, x_min, x_max, y_min, y_max)
        self.rejection_radius = rejection_radius
        self.internal_sampler = random.Random(seed)
        self.global_mask = global_mask

    def sample(self):
        while (True):

            x = self.internal_sampler.uniform(
                self.quad_tree.x_min,
                self.quad_tree.x_max)

            y = self.internal_sampler.uniform(
                self.quad_tree.y_min,
                self.quad_tree.y_max)

            if not self.global_mask(x,y):
                continue

            node = self.quad_tree.find_node(x,y)
            neighborhood = self.quad_tree.find_neighborhood(x,y)
            neighborhood.append(node)

            too_close = False
            for adjacent_node in neighborhood:
                for point in adjacent_node.data:
                    if (point[0] - x)**2 + (point[1] - y)**2 < (self.rejection_radius **2):
                        too_close = True
                        break

                if too_close:
                    break

            if (not too_close):
                result = (x,y)
                print(result)
                node.data.append(result)
                return result


class NetworkRenderer:

    def __init__(
            self,
            network_loader,
            species_of_interest,
            reactions_of_interest,
            output_file,
            width=1024,
            height=1024,
            rejection_radius = 0.008,
            node_radius = 0.002,
            background_line_width=0.001,
            global_mask_radius=0.47,
            # size of reaction batches fetched from the database
            reaction_batch_size = 10000,
            # for the background, there is no need to render every reaction
            # instead we render reactions according to the reaction probability
            reaction_probability = 1,
            colors = [(x,x,x) for x in [0.3,0.4,0.5,0.6,0.7,0.8]],
    ):
        """
        species of interest is a dict mapping species ids to
        their canvas coordinates and styling.

        reactions of interest is a dict mapping reaction ids
        to their styling.
        """

        self.network_loader = network_loader
        self.species_of_interest = species_of_interest
        self.reactions_of_interest = reactions_of_interest
        self.output_file = output_file
        self.rejection_radius = rejection_radius
        self.width = width
        self.height = height
        self.reaction_probability = reaction_probability
        self.colors = colors
        self.repulsive_sampler = RepulsiveSampler(
            rejection_radius,
            0.0,
            1.0,
            0.0,
            1.0,
            lambda x, y: (
                True if (x - 0.5)**2 + (y - 0.5)**2 < global_mask_radius**2
                else False )
        )
        self.node_radius = node_radius
        self.background_line_width = background_line_width


        self.surface = cairo.ImageSurface(cairo.Format.ARGB32, width, height)
        self.context = cairo.Context(self.surface)
        self.context.scale(width, height)

        self.species_locations = {}

        if self.species_of_interest is not None:
            for species_id in self.species_of_interest:
                position = self.species_of_interest[species_id]
                self.species_locations[species_id] = position
                node = self.repulsive_sampler.quad_tree.find_node(*position)
                node.data.append(position)

        for i in range(self.network_loader.number_of_species):
            if i not in self.species_locations:
                self.species_locations[i] = self.repulsive_sampler.sample()


        self.reaction_batch_size = reaction_batch_size


    def render_reaction(self, reaction, color):
        context = self.context
        for i in range(reaction['number_of_reactants']):
            for j in range(reaction['number_of_products']):
                if i <= j:
                    reactant_index = reaction['reactants'][i]
                    product_index = reaction['products'][j]

                    context.set_source_rgb(*color)
                    context.move_to(*self.species_locations[reactant_index])
                    context.line_to(*self.species_locations[product_index])
                    context.stroke()



    def render(self):

        local_sampler = random.Random(42)
        context = self.context

        context.set_line_width(self.background_line_width)
        current_base_reaction = 0

        while (current_base_reaction < self.network_loader.number_of_reactions):

            reactions = self.network_loader.get_reactions_in_range(
                current_base_reaction,
                current_base_reaction + self.reaction_batch_size)


            for reaction in reactions:
                print(reaction['reaction_id'])
                if (local_sampler.random() < self.reaction_probability):
                    color = local_sampler.choice(self.colors)
                    self.render_reaction(reaction, color)

            current_base_reaction += self.reaction_batch_size


        if self.reactions_of_interest is not None:
            for reaction_id in self.reactions_of_interest:
                reaction = self.network_loader.index_to_reaction(reaction_id)
                self.render_reaction(reaction, self.reactions_of_interest[reaction_id])



        context.set_source_rgb(0,0,0)
        # plot species nodes
        for i in range(self.network_loader.number_of_species):
            context.arc(self.species_locations[i][0],
                        self.species_locations[i][1],
                        self.node_radius,
                        0,
                        2 * math.pi)

            context.fill()



        self.surface.write_to_png(self.output_file)

