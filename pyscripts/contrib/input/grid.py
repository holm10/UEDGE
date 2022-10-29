# Grid package for UEDGE, grid.py
# Separated from 160299 master input by holm10
# Changelog
# 200213 - file created
from uedge import bbb,com,grd,flx


def snull(  gridpath='.', plate=None, nx_ileg=12, nx_oleg=10, nxcore_inside=14, nxcore_outside=14, nycore=15, nysol=16,
            gengrid=1, isnonog=1, psi0min1=.90, psi0min2=0.98, psi0sep=1.00005, psi0max=1.10,nxxpt=1,nxmod=3,
            alfxptu=0.75, alfcy=2.5, dmix0=0.5
         ):
    from numpy import array,pi
    '''
    Function for defining the grid setup in UEGDE
    '''
    
    # Set grid paths if V7.08 is being used: not available for earlier versions
    # For previous versions, place EFIT files in data folder
    try:
        com.aeqdskfname=gridpath+"/aeqdsk" # Efit equilibrium
        com.geqdskfname=gridpath+"/neqdsk" # Efit equilibrium
    except: 
        pass
    """=====================================================================================================
    GRID
    ====================================================================================================="""
    # Geometry definition
    #- - - - - - - - - - -
    bbb.ngrid=1             # Number of grids to allocate space for: odesetup.m L6488

    bbb.mhdgeo=1            # Grid geometry:
                                #=-2 mag mirror (FRC-annulus)   
                                #=-1 cartesian geometry
                                #=0: cylindrical geometry
                                #=1: toroidal MHD equilibrium
                                #=2: toroidal circular limiter

    com.geometry="snull"    # Magnetic configuration
                                #'snull': lower single null
                                #'uppersn': upper single null
                                #'dnbot': bottom half of double null
                                #'dnull': full double null
                                #'snowflake15': Ryutov's theta~15 deg
                                #'snowflake45': Ryutov's theta~45 deg
                                #'snowflake75': Ryutov's theta~75 deg
                                #'dnXtarget': dnbot with Xtarget

    com.isnonog=isnonog           # switch determining if 9-point differencing is on for non-orthogonal meshes


    # Generate grid
    bbb.gengrid=gengrid           #1= generates grid, 0=restores grid from gridue 


    # Magnetic fluxes
    #- - - - - - - - -
    flx.psi0min1=psi0min1        # Normalized flux value at innermost core flux surface
    flx.psi0min2=psi0min2        # Normalized flux value at innermost PF surface
    flx.psi0sep=psi0sep     # Normalized flux value at separatrix flux surface (just slightly outside)
    flx.psi0max=psi0max        # Normalized flux value at wall on outboard side of magnetic axis



    # Cells
    #- - - - 
    com.nxleg[0,]=[nx_ileg,nx_oleg]                 # Number of poloidal cells in 1st divertor leg [inside, outside]
    com.nxcore[0,]=[nxcore_inside,nxcore_outside]   # Number of poloidal cells along 1st core boundary [inside, outside]
    com.nycore[0]=nycore                            # Number of radial zones in 1st core of plasma
    com.nysol[0]=nysol                              # Number of radial zones in 1st SOL

    com.nxxpt=nxxpt             # Number of extra poloidal cells at X-point (per quadrant)
    grd.nxmod=nxmod             # Number of "upstream" poloidal cells (per quadrant) in the
                            #   original mesh that are modified by subroutine refinex

    bbb.isybdryog=1         #=1 sets fx0, fmx stencil to orthogonal values at iy=0 & iy=ny

    # Shaping
    #- - - - -
    grd.alfxptu=alfxptu        # Variable for extra X-point grid spacign above X-point:
                            #  frac=(i/(nxxpt+nxmod))**alfxptu 

    flx.alfcy=alfcy           # SOL flux contour distribution:
                            #  <1: uniform
                            #  >1: concentrated near sep

    grd.slpxt=1.2           # Slope enchantment factor for x(t) near core crown

    grd.kxmesh=1            # X-mesh definition model:
                                #=0: old model (manual def. of seed points)
                                #=1: linear*rational form for x(t)
                                #=2: linear*exponential form for x(t)
                                #=3: spline form for x(t)
                                #=4: exponential+spline form for x(t)

    flx.istchkon=1          # Switches for imposing limits on polar angle about (rmagx,zmagx)
    if flx.istchkon==1:     # Limts imposed
        flx.dtheta_exclude    = array([.75,.50])*pi
                            # angular width of region where SOL flux contours are excluded
                            #   [inboard flux contours, outboard flux contours]
        flx.dtheta_overlap_pf = array([.05,.01])*pi
                            # angular width over which p.f. flux contours can overlap
                            #   with flux contours in the adjacent region.
        flx.dtheta_overlap_sol= array([0.25,0.25])*pi
                            # angular width over which SOL flux contours can overlap
                            #   with flux contours in the adjacent region.
                            #   [inboard flux contours, outboard flux contours]
        # com.theta_split=np.pi/2   # Computed poloidal angle where in-/outboard mesh regions meet
        # flx.thetax=-1.88      # Computed poloidal angle of X-point relative to magnetic axis

    flx.altsearch=1         # Search path for PF surfaces:
                                #=0: search vertically up toward x-point
                                #=1: search vertically down from x-point
                                #=2: search diagonally down and in from x-point

    com.ismmon=3            # Mesh modification:
                                #=0: strictly orthogonal mesh and divertor plates
                                #=1: non-orthogonal mesh, compressed distrib'n on each surface
                                #=2: non-orthogonal mesh, standard distrib'n on all surfaces
                                #=3: combination of options 1 and 2 using weight factor wtmesh1
    if com.ismmon==3:       # Using weight factor
        grd.wtmesh1=0.0     # Weight factor; =1: ismmon=1, =0: ismmon2

    grd.istream=0           # Parameter dir fixed upstream reference surface
                                #=0: midplane+cut(ismmon=1) or top-of-mesh(ismmon=2)
                                #=1: user-defined upstream surface arrays

    grd.nsmooth=4           # Number of times to apply the smoothing algorithm to each
                            #   angle-like surface after non-orthogonal plate construction

    grd.dmix0=dmix0           # Normalized poloidal mixing length for combining mesh0 with mesh 12
                                #=0: abrupt  change from orthogonal mesh to mesh12 at upstream position
                                #=1: gradual change from orthogonal mesh to mesh12 between upstream an downstram pos    

    # Plates
    #- - - - 
    if plate==None:
        grd.iplate=0
    else:
        grd.iplate=1            # Divertor plate definition
                                    #=0: Orthogonal plates
                                    #=1: User-defined plates
        plate()  # Import plate geo


def uppersn(  gridpath='.', plate=None, nx_ileg=12, nx_oleg=10, nxcore_inside=14, nxcore_outside=14, nycore=15, nysol=16,
            gengrid=1, isnonog=1, psi0min1=.90, psi0min2=0.98, psi0sep=1.00005, psi0max=1.10,nxxpt=1,nxmod=3,
            alfxptu=0.75, alfcy=2.5, dmix0=0.5
         ):
    from numpy import array,pi
    '''
    Function for defining the grid setup in UEGDE
    '''
    
    # Set grid paths if V7.08 is being used: not available for earlier versions
    # For previous versions, place EFIT files in data folder
    try:
        com.aeqdskfname=gridpath+"/aeqdsk" # Efit equilibrium
        com.geqdskfname=gridpath+"/neqdsk" # Efit equilibrium
    except: 
        pass
    """=====================================================================================================
    GRID
    ====================================================================================================="""
    # Geometry definition
    #- - - - - - - - - - -
    bbb.ngrid=1             # Number of grids to allocate space for: odesetup.m L6488

    bbb.mhdgeo=1            # Grid geometry:
                                #=-2 mag mirror (FRC-annulus)   
                                #=-1 cartesian geometry
                                #=0: cylindrical geometry
                                #=1: toroidal MHD equilibrium
                                #=2: toroidal circular limiter

    com.geometry="uppersn"    # Magnetic configuration
                                #'snull': lower single null
                                #'uppersn': upper single null
                                #'dnbot': bottom half of double null
                                #'dnull': full double null
                                #'snowflake15': Ryutov's theta~15 deg
                                #'snowflake45': Ryutov's theta~45 deg
                                #'snowflake75': Ryutov's theta~75 deg
                                #'dnXtarget': dnbot with Xtarget

    com.isnonog=isnonog           # switch determining if 9-point differencing is on for non-orthogonal meshes


    # Generate grid
    bbb.gengrid=gengrid           #1= generates grid, 0=restores grid from gridue 


    # Magnetic fluxes
    #- - - - - - - - -
    flx.psi0min1=psi0min1        # Normalized flux value at innermost core flux surface
    flx.psi0min2=psi0min2        # Normalized flux value at innermost PF surface
    flx.psi0sep=psi0sep     # Normalized flux value at separatrix flux surface (just slightly outside)
    flx.psi0max=psi0max        # Normalized flux value at wall on outboard side of magnetic axis



    # Cells
    #- - - - 
    com.nxleg[0,]=[nx_ileg,nx_oleg]                 # Number of poloidal cells in 1st divertor leg [inside, outside]
    com.nxcore[0,]=[nxcore_inside,nxcore_outside]   # Number of poloidal cells along 1st core boundary [inside, outside]
    com.nycore[0]=nycore                            # Number of radial zones in 1st core of plasma
    com.nysol[0]=nysol                              # Number of radial zones in 1st SOL

    com.nxxpt=nxxpt             # Number of extra poloidal cells at X-point (per quadrant)
    grd.nxmod=nxmod             # Number of "upstream" poloidal cells (per quadrant) in the
                            #   original mesh that are modified by subroutine refinex

    bbb.isybdryog=1         #=1 sets fx0, fmx stencil to orthogonal values at iy=0 & iy=ny

    # Shaping
    #- - - - -
    grd.alfxptu=alfxptu        # Variable for extra X-point grid spacign above X-point:
                            #  frac=(i/(nxxpt+nxmod))**alfxptu 

    flx.alfcy=alfcy           # SOL flux contour distribution:
                            #  <1: uniform
                            #  >1: concentrated near sep

    grd.slpxt=1.2           # Slope enchantment factor for x(t) near core crown

    grd.kxmesh=1            # X-mesh definition model:
                                #=0: old model (manual def. of seed points)
                                #=1: linear*rational form for x(t)
                                #=2: linear*exponential form for x(t)
                                #=3: spline form for x(t)
                                #=4: exponential+spline form for x(t)

    flx.istchkon=1          # Switches for imposing limits on polar angle about (rmagx,zmagx)
    if flx.istchkon==1:     # Limts imposed
        flx.dtheta_exclude    = array([.75,.50])*pi
                            # angular width of region where SOL flux contours are excluded
                            #   [inboard flux contours, outboard flux contours]
        flx.dtheta_overlap_pf = array([.05,.01])*pi
                            # angular width over which p.f. flux contours can overlap
                            #   with flux contours in the adjacent region.
        flx.dtheta_overlap_sol= array([0.25,0.25])*pi
                            # angular width over which SOL flux contours can overlap
                            #   with flux contours in the adjacent region.
                            #   [inboard flux contours, outboard flux contours]
        # com.theta_split=np.pi/2   # Computed poloidal angle where in-/outboard mesh regions meet
        # flx.thetax=-1.88      # Computed poloidal angle of X-point relative to magnetic axis

    flx.altsearch=1         # Search path for PF surfaces:
                                #=0: search vertically up toward x-point
                                #=1: search vertically down from x-point
                                #=2: search diagonally down and in from x-point

    com.ismmon=3            # Mesh modification:
                                #=0: strictly orthogonal mesh and divertor plates
                                #=1: non-orthogonal mesh, compressed distrib'n on each surface
                                #=2: non-orthogonal mesh, standard distrib'n on all surfaces
                                #=3: combination of options 1 and 2 using weight factor wtmesh1
    if com.ismmon==3:       # Using weight factor
        grd.wtmesh1=0.0     # Weight factor; =1: ismmon=1, =0: ismmon2

    grd.istream=0           # Parameter dir fixed upstream reference surface
                                #=0: midplane+cut(ismmon=1) or top-of-mesh(ismmon=2)
                                #=1: user-defined upstream surface arrays

    grd.nsmooth=4           # Number of times to apply the smoothing algorithm to each
                            #   angle-like surface after non-orthogonal plate construction

    grd.dmix0=dmix0           # Normalized poloidal mixing length for combining mesh0 with mesh 12
                                #=0: abrupt  change from orthogonal mesh to mesh12 at upstream position
                                #=1: gradual change from orthogonal mesh to mesh12 between upstream an downstram pos    

    # Plates
    #- - - - 
    if plate==None:
        grd.iplate=0
    else:
        grd.iplate=1            # Divertor plate definition
                                    #=0: Orthogonal plates
                                    #=1: User-defined plates
        plate()  # Import plate geo


def slab_neumann(   radx=5e-2, rad0=0, radm=-1e-2, za0=0, zaxpt=2, zax=3, alfyt=-2, alfxt=4.0,
                    isadjalfxt=0, btfix=2, bpolfix=0.2, gengrid=1, nxxpt=0,nxmod=2,alfxptu=1, alfcy=1e-4,
                    nx_oleg=32, nxcore_outside=12, nycore=6, nysol=12, dmix0=0,dispx=0,dispy=0):
    '''
    Function for defining the grid setup in UEGDE
    '''
    """=====================================================================================================
    GRID
    ====================================================================================================="""
    # Geometry definition
    #- - - - - - - - - - -
    bbb.ngrid=1             # Number of grids to allocate space for: odesetup.m L6488

    bbb.mhdgeo=-1           # Grid geometry:
                                #=-2 mag mirror (FRC-annulus)   
                                #=-1 cartesian geometry
                                #=0: cylindrical geometry
                                #=1: toroidal MHD equilibrium
                                #=2: toroidal circular limiter


    bbb.isfixlb[0]=2    # =1 fixes values on left boundary (nib, upb, teb, tib, yylb0)
                        # =2 for symmetry point
    grd.radx=dispy+radx       # Maximum radius of cylinder or outer wall location for slab [m]
    grd.rad0=dispy+rad0       # Radial separatrix location for cylinder and slab [m]
    grd.radm=dispy+radm      # Minimum radius of cylinder or inner wall location for slab [m]
    grd.za0=dispx+za0           # Axial position of LHB
    grd.zaxpt=dispx+zaxpt         # Axial position of XPT
    grd.zax=dispx+zax           # Axial position of RHB
    grd.alfyt=alfyt        # Radial nonuniformity factor
    grd.alfxt=alfxt   # Axial nonuniformity factor
    grd.isadjalfxt=isadjalfxt    # Alter alfxt for smooth dx at XPT
    grd.btfix=btfix         # Total B-field for slab
    grd.bpolfix=bpolfix     # Poloidal B-field for slab

    com.geometry="snull"    # Magnetic configuration
                                #'snull': lower single null
                                #'uppersn': upper single null
                                #'dnbot': bottom half of double null
                                #'dnull': full double null
                                #'snowflake15': Ryutov's theta~15 deg
                                #'snowflake45': Ryutov's theta~45 deg
                                #'snowflake75': Ryutov's theta~75 deg
                                #'dnXtarget': dnbot with Xtarget

    com.isnonog=0           # switch determining if 9-point differencing is on for non-orthogonal meshes


    # Generate grid
    bbb.gengrid=gengrid           #1= generates grid, 0=restores grid from gridue 


    # Cells
    #- - - - 
    com.nxleg[0,]=[0,nx_oleg]                 # Number of poloidal cells in 1st divertor leg [inside, outside]
    com.nxcore[0,]=[0,nxcore_outside]   # Number of poloidal cells along 1st core boundary [inside, outside]
    com.nycore[0]=nycore                            # Number of radial zones in 1st core of plasma
    com.nysol[0]=nysol                              # Number of radial zones in 1st SOL

    com.nxxpt=nxxpt             # Number of extra poloidal cells at X-point (per quadrant)
    grd.nxmod=nxmod             # Number of "upstream" poloidal cells (per quadrant) in the
                            #   original mesh that are modified by subroutine refinex

    bbb.isybdryog=0         #=1 sets fx0, fmx stencil to orthogonal values at iy=0 & iy=ny

    # Shaping
    #- - - - -
    grd.alfxptu=alfxptu        # Variable for extra X-point grid spacign above X-point:
                            #  frac=(i/(nxxpt+nxmod))**alfxptu 

    flx.alfcy=alfcy           # SOL flux contour distribution:
                            #  <1: uniform
                            #  >1: concentrated near sep

    grd.slpxt=1           # Slope enchantment factor for x(t) near core crown

    grd.kxmesh=1            # X-mesh definition model:
                                #=0: old model (manual def. of seed points)
                                #=1: linear*rational form for x(t)
                                #=2: linear*exponential form for x(t)
                                #=3: spline form for x(t)
                                #=4: exponential+spline form for x(t)

    flx.istchkon=0

    flx.altsearch=0         # Search path for PF surfaces:
                                #=0: search vertically up toward x-point
                                #=1: search vertically down from x-point
                                #=2: search diagonally down and in from x-point

    com.ismmon=0            # Mesh modification:
                                #=0: strictly orthogonal mesh and divertor plates
                                #=1: non-orthogonal mesh, compressed distrib'n on each surface
                                #=2: non-orthogonal mesh, standard distrib'n on all surfaces
                                #=3: combination of options 1 and 2 using weight factor wtmesh1


    grd.istream=0           # Parameter dir fixed upstream reference surface
                                #=0: midplane+cut(ismmon=1) or top-of-mesh(ismmon=2)
                                #=1: user-defined upstream surface arrays

    grd.nsmooth=2           # Number of times to apply the smoothing algorithm to each
                            #   angle-like surface after non-orthogonal plate construction

    grd.dmix0=dmix0           # Normalized poloidal mixing length for combining mesh0 with mesh 12
                                #=0: abrupt  change from orthogonal mesh to mesh12 at upstream position
                                #=1: gradual change from orthogonal mesh to mesh12 between upstream an downstram pos    

    grd.iplate=0


def slab_neumann_SOLPS(   radx=8e-2, rad0=0, radm=-2e-2, za0=0, zaxpt=2, zax=3, alfyt=-2, alfxt=4.0,
                    isadjalfxt=0, btor=2, bpol=0.2, gengrid=1, nxxpt=0,nxmod=2,alfxptu=1, alfcy=1e-4,
                    nxleg=32, nxcore=12, nycore=6, nysol=12, dmix0=0,dispx=0,dispy=0):


    # Calculate tot b-field based on toroidal an poloidal components
    btot = (btor**2 + bpol**2)
    # Calculate correct slab length for mirroring in middle of upstream cell
    dc = ((zaxpt/(nxcore-0.5))*nxcore - zaxpt)
    # Add required extension
    zax += dc
    zaxpt += dc

    slab_neumann(radx=radx, rad0=rad0, radm=radm, za0=za0, zaxpt=zaxpt, zax=zax, alfyt=alfyt, alfxt=alfxt,
                    isadjalfxt=isadjalfxt, btfix=btot, bpolfix=bpol, gengrid=gengrid, nxxpt=nxxpt,nxmod=nxmod,alfxptu=alfxptu, alfcy=alfcy,
                    nx_oleg=nxleg, nxcore_outside=nxcore, nycore=nycore, nysol=nysol, dmix0=dmix0,dispx=dispx,dispy=dispy)



def mesh_seed_points():
    ''' Compute mesh seed points using subroutine exponseed
        
        Adapted for PyUEDGE based on Jaervinen's 174270 H-mode case
    '''
    grd.isupstreamx = 1 # Sets X-pt as limit for nonorthogonal mesh
    com.reset_core_og = 1 # Forces orthogonal diff above X-point
 #   grd.kxmesh=0 # Manual definition of seed points
#    grd.kxmesh = 0
    grd.fraclplt[1] = 0.6      #frac div leg with near-plate spacing
    grd.alfxdiv[1] = 0.21      #expon factor for cell-size expansion
    grd.alfxcore[1] = 0.4      #expon factor for cell-size expansion
    grd.shift_seed_leg[:2] = [0.2,0.25]  
    grd.shift_seed_core[:2] = [1.2,1.0]
    grd.nxlplt[:2] = [6,6]  #num poloidal cells in leg-plate region
    grd.nxlxpt[:2] = [4,4]    #num pol cells in leg-xpt region; note nxplt+nxlxpt=nxleg must be satisfied
    grd.fcorenunif=0.4 # Fraction of poloidal mesh with exponential scaling 
    
    bbb.issfon=0;bbb.ftol=1e20;bbb.exmain() # Generate all arrays
    bbb.issfon=1;bbb.ftol=1e-8 # Restore defaults
    grd.kxmesh=0
    grd.exponseed()
   

def plasma_box(w=1e-1,res=40,dispx=0,dispy=0):
    ''' Defines a square plasma box with uniform grid'''
 
    bbb.mhdgeo=-1		#set cartesian com.geometry
    # Set the com.geometry
    grd.radx = dispy+w		#outer "radial" wall
    grd.radm = dispy		#minimum "radial" position
    grd.rad0 = dispy+w/2
    grd.alfyt=1e-10		#radial nonuniformity factor < 0 => expanding
    grd.za0 = 0.		#poloidal symmetry plane location
    grd.zaxpt = dispx		#poloidal location of flx.x-point
    grd.zax = dispx+w		#poloidal location of divertor plate
    grd.alfxt=1e-10

    grd.btfix = 5.		#constant total B-field
    grd.bpolfix = 0.3     	#constant poloidal B-field
    bbb.ngrid = 1
    com.nxleg[0,0]=0
    com.nxcore[0,0]=0
    com.nxcore[0,1]=0
    com.nxleg[0,1]=res
    com.nycore[0]=0
    com.nysol=res


def plasma_tube(l=1,w=1e-2,nx=100,ny=10,dispx=0,dispy=0):
    ''' Defines a square plasma box with uniform grid'''
 
    bbb.mhdgeo=-1		#set cartesian com.geometry
    # Set the com.geometry
    grd.radx = dispy+w		#outer "radial" wall
    grd.radm = dispy		#minimum "radial" position
    grd.rad0 = dispy+w/2
    grd.alfyt=1e-10		#radial nonuniformity factor < 0 => expanding
    grd.za0 = 0.		#poloidal symmetry plane location
    grd.zaxpt = dispx		#poloidal location of flx.x-point
    grd.zax = dispx+l		#poloidal location of divertor plate
    grd.alfxt=3

    grd.btfix = 5.		#constant total B-field
    grd.bpolfix = 0.3     	#constant poloidal B-field
    bbb.ngrid = 1
    com.nxleg[0,0]=0
    com.nxcore[0,0]=0
    com.nxcore[0,1]=0
    com.nxleg[0,1]=nx
    com.nycore[0]=0
    com.nysol=ny



