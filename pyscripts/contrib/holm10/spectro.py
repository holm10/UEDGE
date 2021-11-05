from matplotlib.pyplot import ion
ion()

def save_grid(grid):
    print('TBD: function to save processed grids - crm computationally intensive')

class Grid():
    def __init__(self):
        from uedge import com, bbb
        from numpy import zeros
        self.cells = []
        self.ue2arr = zeros((com.nx+1, com.ny+1))
        self.arr2ue = zeros((com.nx*com.ny,2))
        # TODO: Append cell objects from grid
        for ix in range(1,com.nx+1):
            for iy in range(1,com.ny+1):
                verts = []
                for iv in [1, 2, 4, 3]:
                    verts.append([com.rm[ix, iy, iv], com.zm[ix, iy, iv]])
                self.ue2arr[ix, iy] = len(self.cells)
                self.arr2ue[len(self.cells)] = [ix, iy]
                self.cells.append(Cell(verts, (ix,iy), te = bbb.te[ix,iy], 
                        ti = bbb.ti[ix, iy], ne = bbb.ne[ix,iy], ni = bbb.ni[ix,iy,0],
                        na = bbb.ng[ix,iy,0], nm = bbb.ng[ix,iy,1]))
        self.ue2arr = self.ue2arr.astype(int)
        self.arr2ue = self.arr2ue.astype(int)

    def plot_grid(self, ax=None):
        ''' Plots a polygrid of Patches '''
        from matplotlib.pyplot import subplots
        if ax is None:
            f, ax = subplots()
            ret = True
        else:
            ret = False

        for cell in self.cells:
            cell.plot_cell(ax)

        if ret is True:
            return f
    
    def get_cell_ue(self, ix, iy):
        return self.cells[ self.ue2arr[ix, iy]]
    
    def get_cell_arr(self, i):
        return self.cells[i]




class Cell():
    ''' Containter object for grid cell info '''
    
    # TODO: create function calculating and storing the spectra
    def __init__(self, vertices, indices, te=0, ti=0, ne=0, ni=0, na=0, nm=0, E=0.1, crmpath=None, crmfname=None):
        from shapely.geometry import Polygon
        self.vertices = vertices
        self.indices = indices
        self.polygon = Polygon(vertices)
        self.crmeval = False
        self.te = te/1.602e-19
        self.ti = ti/1.602e-19
        self.ne = ne*1e-6
        self.ni = ni*1e-6
        self.na = na*1e-6
        self.nm = nm*1e-6
        self.E = E

    def calc_emissivity(self, crm, nH=None):
        if nH is not None:
            crm.species['H(n=1)']['n'] = nH
        else:
            crm.species['H(n=1)']['n'] = self.na
        crm.species['H2(v=0)']['n'] = self.nm
        [self.H_emission, self.H2_emission] = crm.intensity(self.te,self.ne,self.E,self.ti,self.ni)
        self.H_emission[1] *= self.polygon.area
        self.H2_emission[1] *= self.polygon.area
        self.crmeval=True

    def plot_cell(self, ax):
        ax.plot(*self.polygon.exterior.xy, 'k-', linewidth=0.5) 
   



class Chord():
    ''' Creates a single chord object '''
    def __init__(self, points, width, crm=None):
        from shapely.geometry import Point
        self.points = (Point((points[0][0], points[0][1])), 
                        Point((points[1][0], points[1][1])))
        self.width = width
        self.crm = crm
        create_poly()

    def create_poly(self):
        from numpy import sqrt, cos, sin, arctan
        from shapely.geometry import Polygon

        L = sqrt((self.points[1].x - self.points[0].x)**2 + 
                        (self.points[1].y - self.points[0].y)**2)
        theta = arctan((self.points[1].x - self.points[0].x)/(self.points[1].y - 
                        self.points[0].y))
        
        # Elongate the LOS arbitrarily to make sure that it crosses the
        # wall boundary and calculate the new end coordinates, length
        # and end width of the LOS:
        r_elong = self.points[1].x - 1.0 * sin(theta)
        z_elong = self.points[1].y - 1.0 * cos(theta)
        L_elong = sqrt((r_elong - self.points[0].x)**2 + (z_elong - self.points[0].y)**2)
        w_elong = L_elong/L*self.width
        
        # Add the elongated line-of-sight to the list of LOS polygons: 
        self.poly = Polygon([(self.points[0].x, self.points[0].y), 
                                (r_elong - 0.5*w_elong*cos(theta), 
                                    z_elong + 0.5*w_elong*sin(theta)), 
                                (r_elong + 0.5*w_elong*cos(theta), 
                                    z_elong - 0.5*w_elong*sin(theta))])

    def plot_chord(self, ax, poly=True, line=True):
        if line is True:
            (p1, p2) = self.points
            ax.plot([p1.x, p2.x], [p1.y, p2.y], 'r-', linewidth=0.5)
        if poly is True:
             xs, ys = self.poly.exterior.xy    
             ax.fill(xs, ys, alpha=0.2, fc='r', ec='none')

    def plot_LOS_polygons(self, ax):
        for poly in self.LOS_polygons:
             xs, ys = poly.exterior.xy    
             ax.fill(xs, ys, fc='g', ec='none')
            
    def plot_LOS_spectra(self, ax=None, xlim=(500,8500), ylim=(None, None),yaxis='lin'):
        from matplotlib.pyplot import subplots
        from numpy import nonzero
        ret = False
        if ax is None:
            f, ax = subplots()
            ret = True

        if yaxis is 'lin':
            pl = ax.plot
        elif yaxis is 'log':
            pl = ax.semilogy
        
    
        for i in nonzero(self.LOS_H2_emission_integral)[0]:
            ene = self.LOS_H2_energies[i]
            pl([10*ene, 10*ene], [0, self.LOS_H2_emission_integral[i]], 'r-')

        for i in nonzero(self.LOS_H_emission_integral)[0]:
            ene = self.LOS_H_energies[i]
            pl([10*ene, 10*ene], [0, self.LOS_H_emission_integral[i]], 'b-')
        

        ax.set_xlabel(r'Wavelength [Å]')
        ax.set_ylabel('Intensity [counts]')
        ax.set_xlim(xlim)
        ax.set_ylim(ylim)

        if ret is True:
            return f

    def LOS_integral(self, grid, reevaluate=False, nH=None):
        from numpy import sqrt, sum, array
        from tqdm import tqdm
        self.LOS_polygons = []
        self.LOS_cells = [] 
        self.LOS_area = []
        self.LOS_H2_emission = []
        self.LOS_H_emission = []
        self.LOS_dL = []
        self.LOS_L = []
        
        
        # Loop through each inversion grid cell:
        for cell in tqdm(grid.cells):
            
            # Check if the grid cell is completely inside the line-of-
            # sight polygon:
            if self.poly.contains(cell.polygon):
                
                # Add the entire cell polygon to the list of grid cells
                # inside the LOS polygon:
                self.LOS_polygons.append(cell.polygon)
                self.LOS_cells.append(cell)
                
                # Calculate dL for the integration by postulating a
                # rectangular cell orthogonal to the line-of-sight with
                # area of the original cell and width determined by the
                # distance of the cell from the origin of the LOS and
                # the width of the LOS at its end:
                L_cell = sqrt((cell.polygon.centroid.x - self.points[0].x)**2 + 
                                (cell.polygon.centroid.y - self.points[0].y)**2)
                w_orthog = L_cell/sqrt((self.points[1].x - self.points[0].x)**2 + 
                                (self.points[1].y - self.points[0].y)**2)*self.width
                dL = cell.polygon.area/w_orthog
                
                # Store the emission in the grid cell multiplied by dL
                # in the list of emission inside the LOS polygon:
                # TODO: do CRM if it is not done
                if (cell.crmeval is False) or (reevaluate is True):
                    cell.calc_emissivity(self.crm, nH=nH)
                self.LOS_H_emission.append(cell.H_emission[1]*dL)
                self.LOS_H2_emission.append(cell.H2_emission[1]*dL)
                self.LOS_L.append(L_cell)
                self.LOS_dL.append(dL)
                
                self.LOS_area.append(cell.polygon.area)
                self.LOS_H_energies = cell.H_emission[0]
                self.LOS_H2_energies = cell.H2_emission[0]
            
            # Alternatively, check if part of the grid cell is inside
            # the line-of-sight polygon:
            elif self.poly.intersects(cell.polygon):
                
                # Determine the part of the grid cell polygon that
                # intersects with the line-of-sight polygon:
                grid_poly_clipped = self.poly.intersection(cell.polygon)
                
                # Add the intersecting part of the cell polygon to the
                # list of grid cells inside the LOS polygon:
                # TODO: How to do this using objects?
                self.LOS_polygons.append(grid_poly_clipped)
                self.LOS_cells.append(cell)
                
                # Calculate dL for the integration by postulating a
                # rectangular cell orthogonal to the line-of-sight with
                # area of the original cell and width determined by the
                # distance of the cell from the origin of the LOS and
                # the width of the LOS at its end: 
                L_cell = sqrt((grid_poly_clipped.centroid.x - self.points[0].x)**2 + 
                                (grid_poly_clipped.centroid.y - self.points[0].y)**2)
                w_orthog = L_cell/sqrt((self.points[1].x - self.points[0].x)**2 + 
                                (self.points[1].y - self.points[0].y)**2)*self.width
                dL = grid_poly_clipped.area/w_orthog
                
                # Store the emission in the grid cell multiplied by dL
                # in the list of emission inside the LOS polygon:
                if (cell.crmeval is False) or (reevaluate is True):
                    cell.calc_emissivity(self.crm, nH=nH)
                self.LOS_H_emission.append(cell.H_emission[1]*dL)
                self.LOS_H2_emission.append(cell.H2_emission[1]*dL)
                self.LOS_dL.append(dL)
                self.LOS_L.append(L_cell)
                
                self.LOS_area.append(cell.polygon.area)
                self.LOS_H_energies = cell.H_emission[0]
                self.LOS_H2_energies = cell.H2_emission[0]
        
        self.LOS_polygons = array(self.LOS_polygons)
        self.LOS_cells = array(self.LOS_cells)
        self.LOS_area = array(self.LOS_area)
        self.LOS_H_emission = array(self.LOS_H_emission)
        self.LOS_H2_emission = array(self.LOS_H2_emission)
        self.LOS_dL = array(self.LOS_dL)
        self.LOS_L = array(self.LOS_L)
        # Calculate the LOS-integrated emission as a sum of each
        # emission*dL element stored for the current LOS:   
        self.LOS_H_emission_integral = sum(self.LOS_H_emission, axis=0)
        self.LOS_H2_emission_integral = sum(self.LOS_H2_emission, axis=0)
        


class Spectrometer():

    def __init__(self, fname, displ=0, width = 0.017226946, norm_zmag=True, crmfname='input/psi22_spectro.dat', crmpath='/Users/holma2/Dropbox (Aalto)/Analysis/PSI22/'):
        ''' Creates polygons for the chords from a data file '''
        from CRUMPET import Crumpet
        self.crm = Crumpet(path=crmpath, fname=crmfname)

        self.read_chordfile(fname, norm_zmag=norm_zmag, displ=displ, width=width)

        # Loop through the lines-of-sight:
        for chord in self.chords:
            chord.create_poly()
        

    def read_chordfile(self, fname, norm_zmag=True, displ=0, width=0.005):
        from uedge import com
        datin = []
        with open(fname) as f:
            for l in f:
                datin.append(l)
        [x, y, name, date, time]=datin.pop(0).split()

        self.chords = []
        for i in range (2):
            datin.pop(0)
        for c in datin:
            chord = [float(x) for x in c.split()]
            # Displace the read chords by the magnetic axis
            chord[1] += (com.zmagx**norm_zmag + displ)
            chord[3] += (com.zmagx**norm_zmag + displ)
            self.chords.append(Chord(((chord[0], chord[1]), 
                                    (chord[2], chord[3])),
                                width, crm=self.crm))

    def plot_spectrometer(self, ax=None, poly=True, line=True):
        ''' Plots a polygrid of Patches '''
        from matplotlib.pyplot import subplots
        if ax is None:
            f, ax = subplots()
            ret = True
        else:
            ret = False

        for chord in self.chords:
            chord.plot_chord(ax, poly, line)

        if ret is True:
            return f

    def calculate_LOS_integral(self, grid):
        ''' Calculates the LOS integrals from a Grid object '''
        for chord in self.chords:
            chord.LOS_integral(grid)
 


