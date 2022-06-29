from matplotlib.pyplot import ion
ion()

def create_CRM(crmfname='input/psi22_spectro_vib.dat', crmpath='/Users/holma2/Dropbox (Aalto)/Analysis/PSI22/'):
        from CRUMPET import Crumpet
        return Crumpet(path=crmpath, fname=crmfname)
    
def save_grid(grid, savename):
    from pickle import dump
    with open(savename,'wb') as f:
        dump(grid,f)

def restore_grid(savename):
    from pickle import load
    with open(savename,'rb') as f:
        return load(f)

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
                        na = bbb.ng[ix,iy,0], nm = bbb.ng[ix,iy,1],
                        vol = com.vol[ix,iy]))
        self.ue2arr = self.ue2arr.astype(int)
        self.arr2ue = self.arr2ue.astype(int)

    def save(self, savename):
        from pickle import dump
        with open(savename,'wb') as f:
            dump(self,f)
    

    def plot_intensity(self, interval, ax=None,crm=None,zrange=(None,None), mol=True, cbar=True,zscale=1):
        from matplotlib.pyplot import subplots
        from numpy import array, nonzero, transpose, log10, floor
        from tqdm import tqdm
        from matplotlib.colors import Normalize,LogNorm
        from matplotlib.cm import ScalarMappable
        from matplotlib.pyplot import get_cmap,colorbar,figure

        if ax is None:
            f, ax = subplots()
            ret = True
        else:
            ret = False

        # Store all cell intensities in the requested interval to a list
        intensity = []
        for i in tqdm(range(len(self.cells))):
            cell = self.cells[i]
            if cell.crmeval is False:
                cell.calc_emissivity(crm)

            ci = [0,0]
            for i in nonzero(cell.H2_emission[1])[0]:
                ene = cell.H2_emission[0][i]
                if 10*ene > interval[0] and 10*ene < interval[1]:
                    ci[0] += cell.H2_emission[1][i]
            for i in nonzero(cell.H_emission[1])[0]:
                ene = cell.H_emission[0][i]
                if 10*ene > interval[0] and 10*ene < interval[1]:
                    ci[1] += cell.H_emission[1][i]
            intensity.append(ci)
        intensity = zscale*transpose(array(intensity))

        ind = 0
        if mol is False:
            ind = 1
        zmin, zmax = intensity[ind].min(), intensity[ind].max()
        if zrange[0] is not None:
            zmin=zrange[0]
        if zrange[1] is not None:
            zmax=zrange[1]
        Zcol=((log10(intensity[ind])-floor(log10(zmin)))/(floor(log10(zmax))-floor(log10(zmin))))
        cmap=get_cmap('magma')

        for i in range(len(self.cells)):
            cell = self.cells[i]
            xs, ys = cell.polygon.exterior.xy    
            col=cmap(Zcol[i])
            ax.fill(xs, ys, fc=col, ec='none')

        if cbar is True:
            norm = Normalize(vmin=floor(log10(zmin)),vmax=floor(log10(zmax)))
            norm = LogNorm(vmin=zmin,vmax=zmax)
            sm = ScalarMappable(cmap=cmap, norm=norm)
            sm.set_array([])
            cbar=colorbar(sm,ax=ax)
        ax.set_aspect('equal')
        return ax.get_figure()

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


        ax.set_aspect('equal')
        if ret is True:
            return f
    
    def get_cell_ue(self, ix, iy):
        return self.cells[ self.ue2arr[ix, iy]]
    
    def get_cell_arr(self, i):
        return self.cells[i]




    def crm_cell(self, i, cell, crm):
        import time
        # Do shit here
        cell.solve_crm(crm)
        return i, cell

    def store_result(self,result):
            # This is called whenever foo_pool(i) returns a result.
            # result_list is modified only by the main process, not the pool workers.
            self.para_result.append(result)

    def apply_grid_crm(self, crm):
        from multiprocessing import Pool, cpu_count
        self.para_result = []
        pool = Pool(cpu_count())
        total = len(self.cells)
        progress = len(range(1,101))
        for i,cell in enumerate(self.cells[:10]):
#            if (100*i/total>progress[0]):
#                print('{}% done...'.format(ptogress.pop(0)))
            pool.apply_async(self.crm_cell, args = (i, cell, crm), callback = self.store_result)
        pool.close()
        pool.join()
        for res in self.para_result:
            self.cells[res[0]] = res[1]
#        self.para_result.sort(key=lambda x: x[0])
#        self.cells = [x[1] for x in self.para_result] # Copy over ordered list w/ CR results








class Cell():
    ''' Containter object for grid cell info '''
    
    # TODO: create function calculating and storing the spectra
    def __init__(self, vertices, indices, te=0, ti=0, ne=0, ni=0, na=0, nm=0, E=0.1, vol=1, crmpath=None, crmfname=None):
        from shapely.geometry import Polygon
        self.vertices = vertices
        self.indices = indices
        self.polygon = Polygon(vertices)
        self.crmeval = False
        self.crmsol = False
        self.vol = vol
        self.te = te/1.602e-19
        self.ti = ti/1.602e-19
        self.ne = ne*1e-6
        self.ni = ni*1e-6
        self.na = na*1e-6
        self.nm = nm*1e-6
        self.E = E

    def app(self,i):
        self.res=i

    def calc_emissivity(self, crm, nH=None):
        # Calc and store emissivity in 1/s/cm^-3
        if self.crmeval is False:
            self.solve_crm(crm, nH)
        if nH is not None:
            crm.species['H(n=1)']['n'] = nH
        else:
            crm.species['H(n=1)']['n'] = self.na
        crm.species['H2(v=0)']['n'] = self.nm
        # TODO: Figure out a better way to store data to reduce size
        [self.H_emission, self.H2_emission] = crm.intensity(self.te, self.ne, 
                                self.E,self.ti,self.ni, n_ss=self.crm_steady_state)

    def solve_crm(self, crm, nH=None):
        if nH is not None:
            crm.species['H(n=1)']['n'] = nH
        else:
            crm.species['H(n=1)']['n'] = self.na
        crm.species['H2(v=0)']['n'] = self.nm
        self.crm_species = crm.slist
        self.crm_steady_state = crm.steady_state(self.te, self.ne, self.E, self.ti,self.ni)
        self.crmeval=True
            
    def plot_cell(self, ax=None):
        if ax is None:
            f, ax = subplots()
        ax.plot(*self.polygon.exterior.xy, 'k-', linewidth=0.5) 
   



class Chord():
    ''' Creates a single chord object '''
    def __init__(self, points, width, crm=None):
        from shapely.geometry import Point
        self.points = (Point((points[0][0], points[0][1])), 
                        Point((points[1][0], points[1][1])))
        self.width = width
        self.crm = crm
        self.create_poly()
        self.integrated = False

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

    def plot_chord(self, ax=None, poly=True, line=True, color='r', alpha=0.2):
        if ax is None:
            f, ax = subplots()
        if line is True:
            (p1, p2) = self.points
            ax.plot([p1.x, p2.x], [p1.y, p2.y], '-', linewidth=0.5,color=color)
        if poly is True:
             xs, ys = self.poly.exterior.xy    
             ax.fill(xs, ys, alpha=alpha, fc=color, ec='none')

    def plot_LOS_polygons(self, ax=None):
        if ax is None:
            f, ax = subplots()
        for poly in self.LOS_polygons:
             xs, ys = poly.exterior.xy    
             ax.fill(xs, ys, fc='g', ec='none')

    def calc_LOS_spectra(self,lower, upper):
        from numpy import nonzero
        ret = [0,0]
        for i in nonzero(self.LOS_H2_emission_integral)[0]:
            ene = self.LOS_H2_energies[i]
            if 10*ene > lower and 10*ene < upper:
                ret[0] += self.LOS_H2_emission_integral[i]
        for i in nonzero(self.LOS_H_emission_integral)[0]:
            ene = self.LOS_H_energies[i]
            if 10*ene > lower and 10*ene < upper:
                ret[1] += self.LOS_H_emission_integral[i]
        return ret     
            
    def plot_LOS_spectra(self, ax=None, xlim=(500,8500), ylim=(None, None),yaxis='lin',yscale=1):
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
            ene = 10*self.LOS_H2_energies[i]
            if (ene > xlim[0]) and (ene < xlim[1]):
                pl([ene, ene], [0, yscale*self.LOS_H2_emission_integral[i]], 'r-')

        for i in nonzero(self.LOS_H_emission_integral)[0]:
            ene = 10*self.LOS_H_energies[i]
            if (ene > xlim[0]) and (ene < xlim[1]):
                pl([ene, ene], [0, yscale*self.LOS_H_emission_integral[i]], 'b-')
        

        ax.set_xlabel(r'Wavelength [Å]')
        ax.set_ylabel(r'Intensity [ph/sr/s/$\rm cm^3$]')
        ax.set_xlim(xlim)
        ax.set_ylim(ylim)

        if ret is True:
            return f

    def PARA_LOS_integral(self, grid, reevaluate=False, nH=None):
        from multiprocessing import cpu_count
        from joblib import Parallel, delayed
        from tqdm import tqdm
        
        return Parallel(n_jobs=cpu_count())(delayed(
                    self.PARA_LOS_integral_cell)(grid.cells[i], reevaluate=reevaluate, 
                    nH=nH) for i in tqdm(range(len(grid.cells))))
        

    def PARA_LOS_integral_cell(self, cell, reevaluate=False, nH=None, vol=True):
        '''
        vol (bool) : normalize the emission to volume
        '''
        from numpy import sqrt, sum, array, pi
        from tqdm import tqdm
        self.LOS_polygons = []
        self.LOS_cells = [] 
        self.LOS_area = []
        self.LOS_H2_emission = []
        self.LOS_H_emission = []
        self.LOS_dL = []
        self.LOS_L = []
    
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
            self.LOS_H_emission.append(cell.H_emission[1]*dL/(pi*4))
            self.LOS_H2_emission.append(cell.H2_emission[1]*dL/(pi*4))
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
            self.LOS_H_emission.append(cell.H_emission[1]*dL/(pi*4))
            self.LOS_H2_emission.append(cell.H2_emission[1]*dL/(pi*4))
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
    
        return cell.indices


    def LOS_integral_value(self, values, grid):
        from numpy import sqrt, sum, array
        from tqdm import tqdm
        self.LOS_polygons = []
        self.LOS_cells = [] 
        self.LOS_area = []
        self.LOS_value = []
        self.LOS_dL = []
        self.LOS_L = []
        
        
        # Loop through each inversion grid cell:
        for i in tqdm(range(len(grid.cells))):
            cell = grid.cells[i]            

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
                self.LOS_value.append(values[i]*dL)
                self.LOS_L.append(L_cell)
                self.LOS_dL.append(dL)
                self.LOS_area.append(cell.polygon.area)
            
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
                self.LOS_value.append(values[i]*dL)
                self.LOS_dL.append(dL)
                self.LOS_L.append(L_cell)
                self.LOS_area.append(cell.polygon.area)
        
        self.LOS_polygons = array(self.LOS_polygons)
        self.LOS_cells = array(self.LOS_cells)
        self.LOS_area = array(self.LOS_area)
        self.LOS_value = array(self.LOS_value)
        self.LOS_dL = array(self.LOS_dL)
        self.LOS_L = array(self.LOS_L)
        # Calculate the LOS-integrated emission as a sum of each
        # emission*dL element stored for the current LOS:   
        self.LOS_value_integral = sum(self.LOS_value, axis=0)
        




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
        self.integrated = True
        


class Spectrometer():

    def __init__(self, fname, displ=0, width = 0.017226946, norm_zmag=True, crm=None):
        ''' Creates polygons for the chords from a data file '''
        from uedge import com
        from numpy import array
        self.crm=crm
        self.read_chordfile(fname, norm_zmag=norm_zmag, displ=displ, width=width)
        self.osep = array([[com.rm[-2, com.iysptrx,4], com.zm[-2,com.iysptrx,4]]]).transpose()
    

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
            chord[1] += (com.zmagx*norm_zmag + displ)
            chord[3] += (com.zmagx*norm_zmag + displ)
            self.chords.append(Chord(((chord[0], chord[1]), 
                                    (chord[2], chord[3])),
                                width, crm=self.crm))

    def plot_spectrometer(self, ax=None, **kwargs):
        ''' Plots a polygrid of Patches '''
        from matplotlib.pyplot import subplots
        if ax is None:
            f, ax = subplots()
            ret = True
        else:
            ret = False

        for chord in self.chords:
            chord.plot_chord(ax, **kwargs)

        if ret is True:
            return f

    def calculate_LOS_integral(self, grid, **kwargs):
        ''' Calculates the LOS integrals from a Grid object '''
        for chord in self.chords:
            if chord.integrated is False:
                chord.LOS_integral(grid, **kwargs)
    
    def calculate_LOS_integral_data(self, data, grid, **kwargs):
        ret = []
        for chord in self.chords:
            chord.LOS_integral_value(data, grid, **kwargs)
            ret.append(chord.LOS_value_integral)

        return ret

    def plot_chord_angle(self, grid, interval = (5980,6330), mol=True, ax=None, style='ko-',yscale=1,dtheta=0, **kwargs):
        from matplotlib.pyplot import subplots
        from numpy import array, mean, array#, arccos, dot, degrees
        from math import degrees, atan2
        from numpy.linalg import norm
        from shapely.geometry import Point
        if ax is None:
            f, ax = subplots()
        self.calculate_LOS_integral(grid)
        chord_ends = []
        detector = []
        angles = []
        intensities = []
        for chord in self.chords:
            detector.append(array(chord.points[0].xy))
            chord_ends.append(array(chord.points[1].xy))
            # Integrate chord intensity in interval
            intensities.append(chord.calc_LOS_spectra(interval[0], interval[1]))

        # Mean upper point
        detector = mean(array(detector), axis=0)
        for point in chord_ends:
#        for i in range(len(chord_ends)):
#            angles.append(degrees(atan2(chord_ends[i][1]-detector[i][1], chord_ends[i][0]-detector[i][0]) - 
#                    atan2(detector[i][1]-detector[i][1], detector[i][0]+1-detector[0][0])))
            angles.append(degrees(atan2(point[1]-detector[1], point[0]-detector[0]) - 
                    atan2(self.osep[1]-detector[1], self.osep[0]-detector[0])))
#            chord = point - detector
#            sep = self.osep - detector
#            angles.append(degrees(arccos(dot(chord[:,0], sep[:,0]) / (norm(chord) * norm(sep)))))
        # Calculate angle between line and sep
        ax.plot(array(angles)+dtheta, yscale*array(intensities)[:,0**mol],style,**kwargs)
        return ax.get_figure()
                
 

    def plot_chord_angle_data(self,data,  grid, ax=None, style='ko-',yscale=1, dtheta=0, **kwargs):
        from matplotlib.pyplot import subplots
        from numpy import array, mean, array#, arccos, dot, degrees
        from math import degrees, atan2
        from numpy.linalg import norm
        from shapely.geometry import Point
        if ax is None:
            f, ax = subplots()
        chord_ends = []
        detector = []
        angles = []
        for chord in self.chords:
            detector.append(array(chord.points[0].xy))
            chord_ends.append(array(chord.points[1].xy))

        # Mean upper point
        detector = mean(array(detector), axis=0)
#        for i in range(len(chord_ends)):
        for point in chord_ends:
#            angles.append(degrees(atan2(chord_ends[i][1]-detector[i][1], chord_ends[i][0]-1-detector[i][0]) - 
#                    atan2(detector[i][1]-detector[i][1], detector[i][0]+1-detector[0][0])))
            angles.append(degrees(atan2(point[1]-detector[1], point[0]-detector[0]) - 
                    atan2(self.osep[1]-detector[1], self.osep[0]-detector[0])))
#            chord = point - detector
#            sep = self.osep - detector
#            angles.append(degrees(arccos(dot(chord[:,0], sep[:,0]) / (norm(chord) * norm(sep)))))
        # Calculate angle between line and sep
        ax.plot(array(angles)+dtheta, yscale*array(data),style,**kwargs)
        return ax.get_figure()
                
 


