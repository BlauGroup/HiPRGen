from HiPRGen.network_loader import NetworkLoader
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

    def insert(self, x, y, val):
        node = self.find_node(x,y)
        node.data.append(val)
        return val

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



class Renderer:

    def __init__(
            self,
            width=1024,
            height=1024,
            rejection_radius=0.01,
            global_mask_radius=0.47,
            colors = [(x,x,x) for x in [0.3,0.4,0.5,0.6,0.7,0.8]]
    ):

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

        self.local_sampler = random.Random(42)
        self.node_dict = {}

        self.width = width
        self.global_mask_radius = global_mask_radius
        self.height = height
        self.colors = colors

        self.surface = cairo.ImageSurface(cairo.Format.ARGB32, width, height)
        self.context = cairo.Context(self.surface)
        self.context.scale(width, height)

    def new_node(self, tag, point=None):
        # if point is None, a node position will be generated
        # note: if you provide a point, it will go exactly where you say, which
        # may be very close to other points. If tag already used, do nothing.
        if tag not in self.node_dict:

            if point is not None:
                self.node_dict[tag] = (
                    self.repulsive_sampler.quad_tree.insert(
                        point[0],
                        point[1],
                        point))

            else:
                self.node_dict[tag] = self.repulsive_sampler.sample()


    def new_node_boundary(self, tag, angle):
        point = (0.5 + self.global_mask_radius * math.cos(angle),
                 0.5 + self.global_mask_radius * math.sin(angle))

        self.new_node(tag, point=point)

    def draw_node(self, tag, color=(0,0,0), radius=0.0008):
        point = self.node_dict[tag]
        self.context.set_source_rgb(*color)
        self.context.arc(point[0], point[1], radius, 0, 2 * math.pi)
        self.context.fill()

    def draw_node_square(self, tag, color=(0,0,0), side=0.005):
        point = self.node_dict[tag]
        self.context.set_source_rgb(*color)
        self.context.rectangle(point[0] - side/2, point[1] - side/2, side, side)

    def draw_edge(self, tag1, tag2, color=None, width=0.001):
        if color is None:
            color = self.local_sampler.choice(self.colors)


        point1 = self.node_dict[tag1]
        point2 = self.node_dict[tag2]
        self.context.set_source_rgb(*color)
        self.context.set_line_width(width)
        self.context.move_to(*point1)
        self.context.line_to(*point2)
        self.context.stroke()

    def render(self, path):
        self.surface.write_to_png(path)
