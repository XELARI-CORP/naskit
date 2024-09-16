from .circular_graph import CircularGraph
from .svg import DrawSVG



class DrawNA(CircularGraph, DrawSVG):

    def draw(self):
        n = len(self)

        nb_coords, R = self.calculate_circular_coords(n)
        helix_radiuses = self.calculate_helix_radiuses(nb_coords, R)
        svg = self.make_svg(nb_coords, helix_radiuses, R)

        if not hasattr(self, 'svg'):
            self.__dict__['svg'] = svg
        
        return svg
    

    