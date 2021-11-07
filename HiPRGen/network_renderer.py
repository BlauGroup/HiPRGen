from HiPRGen.network_loader import *
import numpy as np
import cairo
import math

class QuadTreeNode:
    """
    origin is at top left so to agree with
    the cairo canvas coordinates
    """
    def __init__(self, x_min, x_max, y_min, y_max):
        self.x_min = x_min
        self.x_max = x_max
        self.y_min = y_min
        self.y_max = y_max

        # top left
        self.quad_1 = []

        # top right
        self.quad_2 = []

        # bottom left
        self.quad_3 = []

        # bottom right
        self.quad_4 = []

    def branch(self):
        self.x_mid = (self.x_min + self.x_max) / 2
        self.x_mid = (self.y_min + self.y_max) / 2

        # top left
        self.quad_1 = QuadTreeNode(
            self.x_min,
            self.x_mid,
            self.y_min,
            self.y_mid)

        # top right
        self.quad_2 = QuadTreeNode(
            self.x_mid,
            self.x_max,
            self.y_min,
            self.y_mid)

        # bottom left
        self.quad_3 = QuadTreeNode(
            self.x_min,
            self.x_mid,
            self.y_mid,
            self.y_max)

        # bottom right
        self.quad_4 = QuadTreeNode(
            self.x_mid,
            self.x_max,
            self.y_mid,
            self.y_max)



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

