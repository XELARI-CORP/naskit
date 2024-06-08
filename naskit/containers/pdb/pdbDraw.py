import numpy as np



class PDBDraw:
    atom_colors_map = {"C":'#1D90DE', "H":"#BCC5E0", "S":"#E0D86B", "O":"#DE371D", "N":"#1DDE81", "P":"#DC8BE0"}
    atom_radius_map = {"C":0.67, "H":0.53, "S":0.88, "O":0.48, "N":0.56, "P":0.98}
    
    def _marker_size_multiplier(self, natoms):
        return 100*np.exp(-(2.e-2*natoms)**.5) + 10
        
    def draw(
            self,
            width=800,
            height=600,
            size_m=None,
            show_axis=False
            ):
        
        import plotly.graph_objects as go
        
        if size_m is None:
            size_m = self._marker_size_multiplier(self.natoms)
            
        c = self.coords
        color = [self.atom_colors_map[a.element] for a in self.atoms()]
        size = [self.atom_radius_map[a.element]*size_m for a in self.atoms()]
        name = [f"{a.name} {a.moln} {a.mol_name}" for a in self.atoms()]

        scatter_go = go.Scatter3d(x=c[:, 0], y=c[:, 1], z=c[:, 2], 
                                  hoverinfo="text",
                                  hovertext=name,
                                  mode='markers', 
                                  marker=dict(size=size, color=color, opacity=1.)
                                 )
        fig = go.Figure(data=[scatter_go], layout=go.Layout(width=width, height=height))
        fig.update_scenes(dragmode="orbit")
        if not show_axis:
            fig.update_scenes(xaxis_visible=False, 
                              yaxis_visible=False, 
                              zaxis_visible=False)
        fig.show()