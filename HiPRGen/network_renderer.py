from HiPRGen.network_loader import *
import cairo

class NetworkRenderer:

    def __init__(
            self,
            network_loader,
            species_of_interest,
            reactions_of_interest,
            output_file,
            width=1024,
            height=1024):
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



    def render(self):

        context = self.context

        x, y, x1, y1 = 0.5, 0.5, 0.1, 0.1
        context.set_line_width(0.01)
        context.move_to(x, y)
        context.line_to(x1,y1)
        context.stroke()
        self.surface.write_to_png(self.output_file)

