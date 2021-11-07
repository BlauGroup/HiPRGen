from HiPRGen.network_loader import *
import numpy as np
import cairo
import math



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
        find all nodes adjacent to our point
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
        adjacent_nodes.append(node)
        return adjacent_nodes

    def find_node(self, x, y):
        """
        find the terminal node so that
        x_min <= x < x_max
        y_min <= y < y_max
        return None if there is no node

        note: does check the terminal node twice,
        but we would stack overflow before that ever posed
        a performance problem.
        """
        if self.quads is not None:
            for quad in self.quads:
                if (quad.x_min <= x < quad.x_max and
                    quad.y_min <= y < quad.y_max):
                    return quad.find_node(x,y)

            return None


        elif self.data is not None:
            if (self.x_min <= x < self.x_max and
                self.y_min <= y < self.y_max):

                return self
            else:
                breakpoint()
                return None


    def __str__(self):
        return (
            "x : [" + str(self.x_min) + ", " + str(self.x_max) + ")  " +
            "y : [" + str(self.y_min) + ", " + str(self.y_max) + ")"
        )

    def __repr__(self):
        return self.__str__()


class NetworkRenderer:

    def __init__(
            self,
            network_loader,
            species_of_interest,
            reactions_of_interest,
            output_file,
            width=1024,
            height=1024
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
        self.width = width
        self.height = height


        self.surface = cairo.ImageSurface(cairo.Format.ARGB32, width, height)
        self.context = cairo.Context(self.surface)
        self.context.scale(width, height)


        self.compute_species_locations()


    def compute_species_locations(self):
        rng = np.random.default_rng()
        self.species_locations = rng.uniform(
            0.1,
            0.9,
            self.network_loader.number_of_species * 2).reshape(
                self.network_loader.number_of_species, 2)


    def render(self):

        context = self.context

        context.set_source_rgb(0,0,0)
        for i in range(self.network_loader.number_of_species):
            context.arc(self.species_locations[i,0],
                        self.species_locations[i,1],
                        0.002,
                        0,
                        2 * math.pi)

            context.fill()

            context.set_line_width(0.01)
            context.move_to(0.5,0.5)
            context.line_to(0.99,0.01)

            context.stroke()

        self.surface.write_to_png(self.output_file)

