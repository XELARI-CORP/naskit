import numpy as np
from typing import Optional



class PDBDraw:
    atom_colors_map = {"C":'#1D90DE', "H":"#BCC5E0", "S":"#E0D86B", "O":"#DE371D", "N":"#1DDE81", "P":"#DC8BE0"}
    default_color = "#111"
    atom_radius_map = {"C":0.67, "H":0.53, "S":0.88, "O":0.48, "N":0.56, "P":0.98}
    default_radius = 0.75


    def _draw_params(self, 
                     draw_hydrogens: bool = True, 
                     size_m: bool = None
                    ):
        if size_m is None:
            size_m = 100*np.exp(-(2.e-2*self.natoms)**.5) + 10
        
        if draw_hydrogens:
            c = self.coords
            color = [self.atom_colors_map.get(a.element, self.default_color) for a in self.atoms()]
            size = [self.atom_radius_map.get(a.element, self.default_radius)*size_m for a in self.atoms()]
            name = [f"{a.name} {a.moln} {a.mol_name}" for a in self.atoms()]
        else:
            c = np.array([a.coords for a in self.atoms() if a.element!='H'], dtype=np.float32)
            color = [self.atom_colors_map.get(a.element, self.default_color) for a in self.atoms() if a.element!='H']
            size = [self.atom_radius_map.get(a.element, self.default_radius)*size_m for a in self.atoms() if a.element!='H']
            name = [f"{a.name} {a.moln} {a.mol_name}" for a in self.atoms() if a.element!='H']

        return {'coords':c, 'colors':color, 'sizes':size, 'names':name}
        
        
    def draw(
            self,
            width: int = 600,
            height: int = 500,
            size_m: int = None,
            draw_hydrogens: bool = True,
            show_axis: bool = False,
            matplot: bool = False
            ):
        
        params = self._draw_params(draw_hydrogens, size_m)
        xs, ys, zs = params['coords'][:, 0], params['coords'][:, 1], params['coords'][:, 2]

        if matplot:
            import matplotlib.pyplot as plt
            fig = plt.figure(figsize=(width/100, height/100))
            ax = fig.add_subplot(projection='3d')
            if not show_axis:
                ax.set_axis_off()
                
            ax.scatter(xs, ys, zs, 
                       s=[8*s for s in params['sizes']], 
                       c=params['colors'],
                       alpha=1.
                      )
            ax.set_aspect('equal', 'box')

        else:
            import plotly.graph_objects as go
            scatter_go = go.Scatter3d(x=xs, y=ys, z=zs, 
                                      hoverinfo="text",
                                      hovertext=params['names'],
                                      mode='markers', 
                                      marker=dict(
                                          size=params['sizes'], 
                                          color=params['colors'], 
                                          opacity=1.
                                      )
                                     )
            fig = go.Figure(data=[scatter_go], layout=go.Layout(width=width, height=height))
            fig.update_scenes(dragmode="orbit", aspectmode="cube", aspectratio=dict(x=1, y=1, z=1))
            if not show_axis:
                fig.update_scenes(xaxis_visible=False, 
                                  yaxis_visible=False, 
                                  zaxis_visible=False)
            fig.show()



