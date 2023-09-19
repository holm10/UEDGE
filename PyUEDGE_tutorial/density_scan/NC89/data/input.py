# Import UEDGE and set paths
from uedge import *
from uedge.hdf5 import hdf5_save,hdf5_restore
from numpy import where,sqrt,linspace
from matplotlib.pyplot import figure
from uedge.contrib.holm10.plot import *



def restore_input():
    '''
    Reads a standard input deck
    Geometry: lower single null
    Machine: DIII-D
    Description: A simple setup on a coarse grid for tutorial purposes
    '''

    """====================================================================================================
    SET PROBLEM NAME    
    ===================================================================================================="""
    bbb.label[0] = "tutorial_nc89" # Define problem name

    """====================================================================================================
    SET REQUIRED DIRECTORY PATHS    
    ===================================================================================================="""


    # SOLVER AND GRID    
    solver()                            # Solver setup
    gridsetup(gridpath='../../grid')         # Grid definition: gridpath must be set if reading EFIT equilibriums

    # SPECIES
    plasma(aphpath='../../rates/')   # Plasma setup: aphpath must be set
    #molecules()                         # Molecular model setup
    allocate()                          # Allocate UEDGE arrays
    #carbon(apipath='../../rates/api')   # Carbon setup: apipath must be set

    # PHYSICS
    recycling_puff_pump()               # Recycling, puffing, and pumping 
    transport_fluxlim()                 # Defines the transport coefficients and sets flux limiting
    currpot()                           # Set currents and potentials: including drifts

    # EQUATIONS TO SOLVE
    equations()

    # RESTORE CASE
    restore(True)                       # Restores the save file 
    

    return 'Input deck successfully read.'
    


def solver():
    '''
    Solver setup for UEDGE case
    '''

    """=====================================================================================================
    SOLVER
    ====================================================================================================="""

    # Solver used
    #- - - - - - - 
    bbb.svrpkg="nksol"
    # Setup solver
    if bbb.svrpkg=="nksol":     # Nonlinear Krylov solver
        bbb.ftol=1.0e-8         # Stop tolerance
        bbb.itermx=30           # Maximum number of "nksol" iterations
        bbb.epscon1=.005        # Linear solve tolerance
        bbb.mfnksol=-3           #nksol method flag; 
                                    #=1 means dogleg strategy
                                    #=2 means linesearch with Arnoldi method,
                                    #=3 means linesearch with GMRES method.
                                    #negative mfnksol ignores global constaints

    elif bbb.svrpkg=="vodpk":   # Vodpk solver
        bbb.runtim=1.0e-07          # convergence for VODPK? (Time for first output)
                                        #=1, set uniform, fixed density, ncore
        bbb.delpy=1e-08             # Sets del: fractional change for finite differences derivative for "vodpk"
        bbb.jacflg=1                # Flag for computing Jacobian in "vodpk"
        bbb.jpre=1                  # Flag for using preconditioning in "vodpk"
        bbb.rtolv=1.e-4             # Relative tolerance vector usid in convert.m
        bbb.inopt=0                 # Resets options for vodpk and daspk solvers
    elif bbb.svrpkg=="newton": # Newton solver
        bbb.rwmin=1.0e-11           # convergence for newton solver
        bbb.nmaxnewt=20             # Maximum number of Newton iterations

    bbb.iscolnorm=3         # Column normalization setup for y's
                                # =0 for no implicit scaling (suscal=1)
                                # =1 for scaling by normalization constants
                                # =2 for scaling by max(abs(yl),floors)
                                # =3 combination of global scaling with nnorm,
                                    #    etc, followed by local scaling by each yl

    # Jacobian increments --
    bbb.xrinc=2             # Jacobian parameter in x-dir, right bound: cell indices
    bbb.xlinc=3             # Jacobian parameter in x-dir, left bound: cell indices
    bbb.yinc=2              # Jacobian parameter in y-dir: cell indices
    # Decomposition
    bbb.premeth="ilut"      # Preconditioning routine for linear iterations:
                                #"banded": use full banded jacobian as preconditiones
                                #"ilut": use ilut preconditioning
                                #"inel": use INEL ILU preconditioning

    bbb.scrit=1.e-3         # Upper limit for re-computing Jacobian and LU decomposition (when sum of residuals>scrit) [Maybe only Newton?]
    bbb.lfililut=200        # Fill-in parameter used in ilut: will allow up to lfililut additional nonzeros in each L/U row
    bbb.lenpfac=200         #  fudge factor to multiply neq by to get an estimate for the number of nonzeros in the preconditioner matrix.
    bbb.lenplufac=200       #  fudge factor to multiply neq by to get an estimate for the number of nonzeros in the factored preconditioner matrix.


    # Differencing methods
    #- - - - - - - - - - - -
    # Differencing method flags for schemes: x-dir scheme=mod(meth,10), y-dir=floor(meth/10)
    # Differencing methods: (Oderhs L100-L
    # 0: Central differencing, assuming no flows
    # 1: Upwind differencing, assuming dif=0
    # 2: Central differencing (harmonic average)
    # 3: Upwind differencing
    # 4: Simple hybrid scheme
    # 5: Fifth power scheme
    # 6: Logarithmic central differencing
    # 7: Regular upwind differencing: can be used fir methg=77 in nonorthogonall differencing
    # 8: Regular upwind differencing on an approximate staggered stencil (for velocities)
    bbb.methn=33    # Ion continuity equation differencing
    bbb.methu=33    # Ion momentum equation differencing
    bbb.methe=33    # Electron energy equation differencing
    bbb.methi=33    # Ion energy equation differencing
    bbb.methg=33    # Neutral gas continuity differencing
    bbb.methp=33    # Potential equation differencing

    # SOLVER OPTIONS AND COEFFICIENTS
    bbb.cngatol=1.  # Coefficient for absolute tolerance for LHS-ode - like routines
    bbb.rlx=0.4     # Fractional change allowed per iteration

    bbb.nurlx=1.e7  # Rate coefficient to relax boundary condition
def gridsetup(gridpath='.'):
    '''
    Function for defining the grid setup in UEGDE
    '''
    com.aeqdskfname=gridpath+"/aeqdsk" # Efit equilibrium
    com.geqdskfname=gridpath+"/neqdsk" # Efit equilibrium

    """=====================================================================================================
    GRID
    ====================================================================================================="""

    # Geometry definition
    #- - - - - - - - - - -
    bbb.ngrid=1             # Number of grids to allocate space for? odesetup.m L6488

    bbb.mhdgeo=1            # Grid geometry:
                                #=-2 mag mirror (FRC-annulus)   
                                #=-1 cartesian geometry
                                #=0: cylindrical geometry
                                #=1: toroidal MHD equilibrium
                                #=2: toroidal circular limiter

    if bbb.mhdgeo==-1:
        bbb.isfixlb[0]=2    # =1 fixes values on left boundary (nib, upb, teb, tib, yylb0)
                            # =2 for symmetry point
        grd.radx=5e-2       # Maximum radius of cylinder or outer wall location for slab [m]
        grd.rad0=0e-2       # Radial separatrix location for cylinder and slab [m]
        grd.radm=-1e-2      # Minimum radius of cylinder or inner wall location for slab [m]
        grd.za0=0           # Axial position of LHB
        grd.zaxpt=2         # Axial position of XPT
        grd.zax=3           # Axial position of RHB
        grd.alfyt=-2        # Radial nonuniformity factor
        grd.alfxt=4.0   # Axial nonuniformity factor
        grd.isadjalfxt=0    # Alter alfxt for smooth dx at XPT
        grd.btfix=2         # Total B-field for slab
        grd.bpolfix=0.2     # Poloidal B-field for slab
        


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



    # Magnetic fluxes
    #- - - - - - - - -
    flx.psi0min1=.98        # Normalized flux value at innermost core flux surface
    flx.psi0min2=.98        # Normalized flux value at innermost PF surface
    flx.psi0sep=1.00001     # Normalized flux value at separatrix flux surface (just slightly outside)
    flx.psi0max=1.07        # Normalized flux value at wall on outboard side of magnetic axis



    # Cells
    #- - - - 
    com.nxleg[0,]=[8,8]   # Number of poloidal cells in 1st divertor leg [inside, outside]
    com.nxcore[0,]=[8,8]  # Number of poloidal cells along 1st core boundary [inside, outside]
    com.nycore[0]=4        # Number of radial zones in 1st core of plasma
    com.nysol[0]=12        # Number of radial zones in 1st SOL

    com.nxxpt=0             # Number of extra poloidal cells at X-point (per quadrant)
    grd.nxmod=2             # Number of "upstream" poloidal cells (per quadrant) in the
                            #   original mesh that are modified by subroutine refinex

    bbb.isybdryog=0         #=1 sets fx0, fmx stencil to orthogonal values at iy=0 & iy=ny

    # Shaping
    #- - - - -
    grd.alfxptu=1        # Variable for extra X-point grid spacign above X-point:
                            #  frac=(i/(nxxpt+nxmod))**alfxptu 

    flx.alfcy=1e-4           # SOL flux contour distribution:
                            #  <1: uniform
                            #  >1: concentrated near sep

    grd.slpxt=1           # Slope enchantment factor for x(t) near core crown

    grd.kxmesh=1            # X-mesh definition model:
                                #=0: old model (manual def. of seed points)
                                #=1: linear*rational form for x(t)
                                #=2: linear*exponential form for x(t)
                                #=3: spline form for x(t)
                                #=4: exponential+spline form for x(t)

    flx.istchkon=0          # Switches for imposing limits on polar angle about (rmagx,zmagx)
    if flx.istchkon==1:     # Limts imposed
        flx.dtheta_exclude    = array([.75,.50])*np.pi
                            # angular width of region where SOL flux contours are excluded
                            #   [inboard flux contours, outboard flux contours]
        flx.dtheta_overlap_pf = array([.05,.01])*np.pi
                            # angular width over which p.f. flux contours can overlap
                            #   with flux contours in the adjacent region.
        flx.dtheta_overlap_sol= array([0.25,0.25])*np.pi
                            # angular width over which SOL flux contours can overlap
                            #   with flux contours in the adjacent region.
                            #   [inboard flux contours, outboard flux contours]
        # com.theta_split=np.pi/2   # Computed poloidal angle where in-/outboard mesh regions meet
        # flx.thetax=-1.88      # Computed poloidal angle of X-point relative to magnetic axis

    flx.altsearch=0         # Search path for PF surfaces:
                                #=0: search vertically up toward x-point
                                #=1: search vertically down from x-point
                                #=2: search diagonally down and in from x-point

    com.ismmon=0            # Mesh modification:
                                #=0: strictly orthogonal mesh and divertor plates
                                #=1: non-orthogonal mesh, compressed distrib'n on each surface
                                #=2: non-orthogonal mesh, standard distrib'n on all surfaces
                                #=3: combination of options 1 and 2 using weight factor wtmesh1
    if com.ismmon==3:       # Using weight factor
        grd.wtmesh1=0.0     # Weight factor; =1: ismmon=1, =0: ismmon2

    grd.istream=0           # Parameter dir fixed upstream reference surface
                                #=0: midplane+cut(ismmon=1) or top-of-mesh(ismmon=2)
                                #=1: user-defined upstream surface arrays

    grd.nsmooth=2           # Number of times to apply the smoothing algorithm to each
                            #   angle-like surface after non-orthogonal plate construction

    grd.dmix0=0           # Normalized poloidal mixing length for combining mesh0 with mesh 12
                                #=0: abrupt  change from orthogonal mesh to mesh12 at upstream position
                                #=1: gradual change from orthogonal mesh to mesh12 between upstream an downstram pos    

    # Plates
    #- - - - 
    # TODO: Specify plate file?
    grd.iplate=0            # Divertor plate definition
                                #=0: Orthogonal plates
                                #=1: User-defined plates
    if grd.iplate==1:       #User-defined plates
        import plate_d3d_10 as pl  # Import plate geo




    #fuzz=2                  # Number of decimals after point - not sure how to set in PyUEDGE

    # Generate grid
    bbb.gengrid=1           #1= generates grid, 0=restores grid from gridue 
def plasma(aphpath='.'):
    """
    Function setting up the UEDGE plasma: ions and atoms
    """


    """=====================================================================================================
    PLASMA SETUP
    ====================================================================================================="""
    # Initalize arrays
    com.nhsp=2          # N.o. hydrogenic species
    com.ngsp=1          # N.o. hydrogenic gaseous specie
    # Set path to rate data
    aph.aphdir=aphpath+'/aph' # Hydrogen rates


    """-----------------------------------------------------------------------------------------------------
    CHARGED SPECIES SETUP
    -----------------------------------------------------------------------------------------------------"""
    bbb.minu[0]=2     # H+ mass in AMU
    bbb.ziin[0]=1     # H+
    bbb.znuclin[0]=1  # H+ nuclear charge 

    # BOUNDARY CONDITIONS
    #-------------------------------------------------------------------------------------------------------

    # Core density BC
    #- - - - - - - - - - - - - -
    bbb.isnicore[0]=1    # Density BC:
                                #=0, set flux to curcore/sy locally in ix
                                #=1, set uniform, fixed density, ncore
                                #=2, set flux & ni over range
                                #=3, set icur=curcore-recycc*fngy, const ni
                                #=4, use impur. source terms (impur only)
                                #=5, set d(ni)/dy=-ni/lynicore          
    if bbb.isnicore[0]==1:# Set up uniform, fixed H+ core density
        bbb.ncore[0]=8.9e19   # H+ core density [m^-3]
        # TODO: Should not matter in these cases (plasma eq's turned off)
    # Core momentum BC
    #- - - - - - - - - - - - - -
    bbb.isupcore[0]=0 # Velocity BC:
                            #=0; Dirichlet BC: up=upcore
                            #=1; Neumann BC: d(up)/dy=0
                            #=2 sets d^2(up)/dy^2 = 0
                            #=3 sets poloidal velocity (uu) = 0
                            #=4 sets tor. ang mom flux = lzflux & n*uz/R=const
                            #=5 sets ave tor vel = utorave & n*uz/R=const

    # Core energy BC
    #- - - - - - - - - - - - - -
    bbb.iflcore=0   # Core power condition:
                        #=0; Temperature dirichlet BC: core Te,i=tcoree,i
                        #=1; Power dirichlet BC:  core power=pcoree,i
                        #=-1; Temperature Neumann BC:  core d(Te,i)/dy=0o

    if bbb.iflcore==0:  # Temperature Dirichlet setup
        bbb.tcoree = 100    # Core electron temperature
        bbb.tcorei = 100    # Core ion temperature 

    if bbb.iflcore==1:  # Power Dirichlet setup
        bbb.pcoree=2.5e5   # Electron power over core boundary
        bbb.pcorei=2.5e5   # Ion power over core boundary


    # PF wall density BC
    #- - - - - - - - - - - - - - - 
    bbb.isextrnpf=0 # Extrapolation BC:
                        #=0; no extrapolation for ni at PF wall
                        #=1; extrapolation BC for ni at PF wall
    bbb.lyni[0]=0.05    # Added outside of loop to get rid of warning on unused lyni!
    bbb.isnwconi[0]=1 # PF wall density BC:
                            #=0, old case; if ifluxni=0, dn/dy=0; if ifluxni=1, fniy=0 (default)
                            #=1, fixed density to nwalli(ix) array
                            #=2, extrapolation B.C.
                            #=3, approx gradient scale length

    if bbb.isnwconi[0]==1: # Fixed inner wall density
        bbb.nwalli[0]=1e18

    if bbb.isnwconi[0]==3:# Gradient scale length BC
        bbb.lyni[0]=0.05    # Fixed scale length at PF wall for BC
        bbb.nwimin[0]=1.e16   # Minimum density limit at PF boundary

    # Outer wall density BC
    #- - - - - - - - - - - - - - - - 
    bbb.isextrnw=0  # Extrapolation BC:
                        #=0; no extrapolation for ni at outer wall
                        #=1; extrapolation BC for ni at outer wall
    bbb.isnwcono[0]=1 # Outer wall density BC:
                            #=0, old case; if ifluxni=0, dn/dy=0; if ifluxni=1, fniy=0 (default)
                            #=1, fixed density to nwallo(ix) array
                            #=2, extrapolation B.C.
                            #=3, approx gradient scale length, limited by nwomin

    if bbb.isnwcono[0]==1: # Fixed inner wall density
        bbb.nwallo[0]=1e18

    if bbb.isnwcono[0]==3:# Gradient scale length BC
        bbb.lyni[1]=0.03    # Fixed scale length at outer wall for BC
        bbb.nwomin[0]=1.e16   # Minimum density limit at outer boundary



    # Wall BC extrapolation
    #- - - - - - - - - - - - - - - - 
    bbb.isextrtpf=0 # Extrapolation BC:
                        #=0; no extrapolation for Ti,Te at PF wall
                        #=1; extrapolation BC for Ti,Te at PF wall
    bbb.isextrtw=0  # Extrapolation BC:
                        #=0; no extrapolation for Ti,Te at outer wall
                        #=1; extrapolation BC for Ti,Te at outer wall


    # PF wall electron  energy BC
    #- - - - - - - - - - - - - - - - - - - 
    bbb.istepfc=3   # PF wall e- temperature BC:
                        #=0, zero energy flux
                        #=1, fixed temp to tedge or tewalli
                        #=2, extrapolation BC
                        #=3, Te scale length
                        #=4, feey = ~bceew*te*elec_flux
    if bbb.istepfc==3:  # Te scale length PF wll BC
        bbb.lyte[0]=0.03# Outer wall electron temperature scale length

    # PF wall ion  energy BC
    #- - - - - - - - - - - - - - - - -
    bbb.istipfc=3   # PF wall H+ temperature BC as above:
    if bbb.istipfc==3:  # Ti scale length PF wall BC
        bbb.lyti[0]=0.1# PF wall H+ temperature scale length


    # Outer wall electron  energy BC
    #- - - - - - - - - - - - - - - - - - - - - 
    bbb.istewc=3    # Outer wall e- temperature BC:
                        #=0, zero energy flux
                        #=1, fixed temp to tedge or tewallo
                        #=2, extrapolation BC
                        #=3, Te scale length
                        #=4, feey = ~bceew*te*elec_flux
    if bbb.istewc==3:   # Te scale length outer wll BC
        bbb.lyte[1]=0.1# Outer wall electron temperature scale length

    # Outer wall ion  energy BC
    #- - - - - - - - - - - - - - - - - - - -
    bbb.istiwc=3    # Outer wall H+ temperature BC as above
    if bbb.istiwc==3:   # Ti scale length outer wall BC
        bbb.lyti[1]=0.1# Outer wall H+ temperature scale length
        

    # Plate BC
    #- - - - - - - - - - 
    bbb.isupss=1    # Plate boundary condition
                        #=-1: dup/dx=0
                        #=0: up=cs
                        #=1: up>=1
    #Flux limits
    bbb.isplflxl=1  # Switch activating flux limits (flalfe/flalfi) at plates (=1)

    # Scale factors
    #- - - - - - - - - - -
    #TODO Figure out if these control anythong: oderhs.m L5536 seems to indicate these vars are only used if ng eq is on, and then only the first index is considered
    # All of the below seem to be considered only for first index (H0) - confirm and change?
    # Combining ion and CX-neutral energy equations (only for com.ngsp=2)??
    bbb.cngtgx=0.   # X&Y-flux coefficient for gaseous component i in ion energy equation
    bbb.cngtgy=0.   # Y-flux coefficient for gaseous component i in ion energy equation


    """-----------------------------------------------------------------------------------------------------
    ATOMIC SPECIS SETUP
    -----------------------------------------------------------------------------------------------------"""
    com.nhgsp=1     # Number of hydrogenic gas species

    bbb.ineudif=2   # Pressure driven neutral transport model
                        #=1 gas sub. neudif uses ng, tg for gas vel & fngx->fnix
                        #=2 gas sub. neudifgp uses pg for gas vel & fngx->fnix
                        #=3 gas sub. neudifl use log_ng, tg for gas vel
                        #  otherwise, old case has ug=ui (strong cx coupling)

    bbb.minu[1]=2    # H0 mass in AMU
    bbb.ziin[1]=0    # H0 inertial neutrals
    bbb.znuclin[1]=1 # H0 nuclear charge 
    bbb.istgcon[0]=0 # H0 temperature: tg=(1-istgcon)*rtg2ti*ti+istgcon*tgas*ev
    bbb.tgas[0]=0    # Neutral temperature scaling is istgcon>1
    if bbb.isupgon[0]==1:# If intertial neutral model

        bbb.cngmom[0]=0  # Momentum loss coefficient for diffusive H0 only
        bbb.cmwall[0]=0  # Momentum wall coefficient for neutral hydrogen only


    # CORE BC
    #- - - - - - - - - 
    # Density
    bbb.isngcore[0]=0    # Neutral gas core boundary conditions
                                #=0, set loc flux= -(1-albedoc)*ng*vtg/4
                                #=1, set uniform, fixed density, ngcore
                                #=2, set rad. grad. to sqrt(lam_i*lam_cx)
                                #=3, extrapolation, but limited
                                #=anything else, set zero deriv which was
                                # prev default inert hy
                                # anything else same as =0
    if bbb.isngcore[0]==1:# Uniform core boundary H0 density
        bbb.ngcore[0]=2.e13  # Core H0 density

    elif bbb.isngcore[0]==0: # Local core flux pumped
        bbb.albedoc[0]=0.5       # Core H0 albedo

    # Momentum
    bbb.isupcore[1]=0    # Velocity BC:
                                #=0; Dirichlet BC: up=upcore
                                #=1; Neumann BC: d(up)/dy=0
                                #=2 sets d^2(up)/dy^2 = 0
                                #=3 sets poloidal velocity (uu) = 0
                                #=4 sets tor. ang mom flux = lzflux & n*uz/R=const

    # Energy 
    bbb.istgcore[0] = 1  # H0 temperature core BC:
                                #=0; set core boundary temperature to ion temperature
                                #=1; set core boundary temperature to tgcore
    if bbb.istgcore[0]==1:   # Specify H0 core boundary temp
        bbb.tgcore[0]=100    # Core D0 boundary temperature


    # Wall boundary conditions
    #- - - - - - - - - - - - - - - - - 
    bbb.tgwall = 4.e-2      # TODO define and place?
    # DNESITY
    # PF wall
    bbb.isnwconi[1]=1    # PF wall density BC as above

    if bbb.isnwconi[1]==0:   # Old case BC
        bbb.ifluxni=1       # Switch for setting PF and outer wall fluxes to 0 (=1)

    if bbb.isnwconi[1]==1: # Fixed inner wall density
        bbb.nwalli[1]=1e18

    # Outer wall
    bbb.isnwcono[1]=1    # Outer wall density BC as above

    if bbb.isnwcono[1]==0:   # Old case BC
        bbb.ifluxni=1       # Switch for setting PF and outer wall fluxes to 0 (=1)

    elif bbb.isnwcono[1]==1: # Fixed inner wall density
        bbb.nwallo[1]=1e18

    elif bbb.isnwcono[1]==3:
        bbb.lyni[1]=3e-2

    # Background gas and normalization
    #- - - - - - - - - - - - - - - - - - - - - 
    bbb.ingb=2      # BG gas source scaling: BG gas source=nuiz*ngbackg*(0.9+0.1(ngbackg/ng)**ingb)
    bbb.ngbackg[0]=1.0e+14   # "soft" artificial floor for neutral densities


    # Scale factors
    #- - - - - - - - - - - 
    #TODO Figure out if these control anythong: oderhs.m L5536 seems to indicate these vars are only used if ng eq is on, and then only the first index is considered
    # All of the below seem to be considered only for first index (H0) - confirm and change?
    # Combining ion and CX-neutral energy equations (only for com.ngsp=2)??
    bbb.cngfx=1.    # Scale factor for flux from Grad_x temperature in neutral gas velocity eqn
    bbb.cngfy=1.    # Scale factor for flux from Grad_y temperature in neutral gas velocity eqn

    bbb.cngflox=1.  # Factor for x-flux from convection in gaseous continuity equation
    bbb.cngfloy=1.  # Factor for y-flux from convection in gaseous continuity equation


    bbb.xstscal=1 # Exponential scale-length with stretched coordinate decays from plates
    bbb.ngscal=0.1   # Ratio of initial gas density to ion density (restart=0)
    bbb.xgscal=1  # Exponential scale of initial gas (restart=0)



    # Inertial neutrals
    #bbb.cngmom=0    # Momentum loss coefficient for diffusive H0 # Defined above
    #bbb.cmwall=0    # Momentum wall coefficient for neutral hydrogen only # Defined above
    bbb.cfbgt=0 # B x Grad(T) coefficient
    bbb.kxn=0   # Poloidal CX-neutral heat diffusivity factor
    bbb.kyn=0   # Radial CX-neutral heat diffusivity factor


    """---------------------------------------------------------------------------------------------------- 
    RATE PARAMETERS
    -----------------------------------------------------------------------------------------------------"""

    # LOOK-UP TABLES FOR HYDROGENIC RATES
    com.istabon = 10    # Rate model/lookup tables:
                            #=0: simple analytic rates and constant energy loss per ionization
                            #=1: table look-up from ADPAK; rates.adpak
                            #=2: table look-up from STRAHL; rates.strahl
                            #=3: table look-up old DEGAS (created 84/09/14); eh.dat & atmc.dat
                            #=4,5,6: table look-up from new DEGAS (created 93/05/06): nwfits
                                #=4:  linear interpolation for  rsa  vs log(te) and log10(ne)
                                #=5:  spline fit for log10(rsa) vs log(te) and log10(ne) 
                                # temp=6:  Hindmarsh spline    log10(rsa) vs log(te) and log10(ne) 
                                #disabled=6_  spline fit for       rsa  vs log(te) and log10(ne)
                            #=7: Campbell's poly. fit for rsa, etc., vs log10(te) and log10(ne)
                            #=8: tab look-up from latest DEGAS (created 93/04/08); ehr1.dat
                            #=9: tab look-up; Stotler PPPL (~95/07/10) with log(Te)-sigv; ehr2.dat
                            #=10: tab look-up;Stotler PPPL (~95/07/10);log(Te)-log(sigv); ehr2.dat
                            #=11: same as istabon=10 with data for n=5-9 excited states; thin.dat
                            #=12: as istabon=11, ex. Lyman-alpha local absb (thick); thickLyA.dat
                            #=13: as istabon=11, ex. all Lyman lines loc absorbed; thickAllLy.dat
                            #=14: H.Scott data; add rtau=Ly-a opacity, lin. interp; ehrtau.dat
                            #=15: H.Scott data; add rtau=Ly-a opacity, log interp; ehrtau.dat
                            #=16: table look-up using 'b2frates_hyd' file data, log-log interp
    if com.istabon==16: # Table lookup used
        com.isrtndep=1      # Are table lookup parameters density dependent: check compabnility with imps!

    # CHARGE-EXCHANGE
    #- - - - - - - - - - -
    bbb.icnucx=0    # CX rate model:
                        #=0: variable nucx
                        #=1: constant nucx=cnucx
                        #=2: use sigcx; nucx~sqrt(Tg)
    if bbb.icnucx==1:   # Constant CX rates
        bbb.cnucx=0.        # Defined constant CX rate

    # IONIZATION
    #- - - - - - - - - - 
    bbb.icnuiz=0    # Ionization model:
                        #=0: variable nuiz
                        #=1: constant nuiz=cnuiz
                        #=2: Freezes (?)
    if bbb.icnuiz==1:   # Constant ionzation frequenzy
        bbb.cnuiz=5.e4      # Defined nuix

    # RECOMBINATION
    #- - - - - - - - - - -
    bbb.isrecmon = 1    # Switch for recombination: 1=on, 0=off
    if bbb.isrecmon==1: # If recombination on
        bbb.cfrecom=1       # Coefficient for recombination frequencies

    # DISSOCIATION
    #- - - - - - - -
    bbb.eion = 2.3            # Energy per atom from dissociation
    bbb.ediss = 2 * bbb.eion    # Electron dissociation loss: 
                                    #<2*eion, virtual energy source through dissociation
                                    #=2*eion, no radiation and power conserved
                                    #>2*eion, virtual energy sink through dissociation (radiation)
 
    """---------------------------------------------------------------------------------------------------- 
    VOLUMETRIC PLASMA SOURCES
    ----------------------------------------------------------------------------------------------------"""
    # Volumetric sources
    bbb.ivolcur=0.0     # Volume source [A] for EACH charge species
    bbb.zwni=1000.      # Width for volume source
    bbb.rwni=1000.      # Width for volume source
def molecules():
    '''
    Function setting up the UEDGE molecular model
    '''


    """====================================================================================================
    MOLECULAR HYDROGEN SETUP
    ===================================================================================================="""
    bbb.ishymol=1       # Includes molecules as 2nd gaseous species (index 1)

    com.nhgsp=com.nhgsp+1       # Allocate space for molecules in hygrogenic gas species array
    com.ngsp=com.ngsp+1     # Allocate space for hydrogen in gas species array
    igh2=com.ngsp-1         # Index for molecules in gaseous arrays
    bbb.istgcon[igh2] = -1.         # Don't reset tg[,,2] using istgcon



    """---------------------------------------------------------------------------------------------
     CORE BC
    ---------------------------------------------------------------------------------------------"""
    # Desntiy
    bbb.isngcore[igh2] = 0  # Hydrogen molecular core BC:
                                #=0, set loc flux= -(1-albedoc)*ng*vtg/4
                                #=1, set uniform, fixed density, ngcore
                                #=2, set rad. grad. to sqrt(lam_i*lam_cx)
                                #=3, extrapolation, but limited
                                #=anything else, set zero deriv which was
                                # prev default inert hy
                                # anything else same as =0
    if bbb.isngcore[igh2]==1:   # Uniform core H2 density
        bbb.ngcore[igh2] = 1.e12    # Core H2 density

    elif bbb.isngcore[igh2]==0: # Local core flux pumped
        bbb.albedoc[igh2]=0.5       # Core H0 albedo

    # Energy
    #- - - - - - - - 
    bbb.istgcore[igh2] = 2  # H2 temperature core BC
                                #=0; set core boundary temperature to ion temperature
                                #=1; set core boundary temperature to tgcore
                                #>1; set zero-temperature gradient over core boundary
    if bbb.istgcore[igh2]==1:# Specify H2 core boundary temp
        bbb.tgcore[igh2] = 100. # Core boundary temp




    """---------------------------------------------------------------------------------------------
    BACKGROUND AND NORMALIZATION
    ---------------------------------------------------------------------------------------------"""
    bbb.n0g[igh2]=1.e17     # Global hydrogen molecule normalization
    bbb.ngbackg[igh2] = 1.e14   # Set background H2 density

    """---------------------------------------------------------------------------------------------
    RATES
    ---------------------------------------------------------------------------------------------"""
    bbb.kelighg[igh2]=5e-16     # Elastic collision coefficient for gas i and hydrogen gas
    bbb.kelighi[igh2]=5e-16     # Elastic collision coeffisient for gas i and hydrogen ions
    
    """---------------------------------------------------------------------------------------------
    SCALE FACTORS
    ---------------------------------------------------------------------------------------------"""
    bbb.cngfx[igh2]=1       # Scale factor for flux from Grad(x)T_g in gas continuity equation
    bbb.cngfy[igh2]=1       # Scale factor for flux from Grad(y)T_g in gas continuity equation
    bbb.cfloxiplt=0         # Coefficient multiplying the neutral convected energy from plates
    # TODO does this belong under H2??
    bbb.cfcvtg=1e0         # Convective thermal transport scale factor
    bbb.cftgcond=1e0       # Conductive thermal transport scale factor
def carbon(apipath='.'):
    """=====================================================================================================
    CARBON SETUP
    ====================================================================================================="""
    # Turn impurity ions on
    bbb.isimpon=6       # Switch for activating impurities
                            #=0: no impurities
                            #=2: fixed-fraction model
                            #=5: Hirshman's reduced-ion model
                            #=6: force-balance model or nusp_imp > 0; see also isofric for full-Z drag term
                            #=7: for simultaneous fixed-fraction and multi-charge-state (isimpon=6) models
    # Set impurita rate path
    api.apidir=apipath+'api' # Impurity rates



    # Helper indices
    #- - - - - - - - 
    species,imps=0,6            # Impurity species and impurity charge states
    iicl=com.nhsp+species*imps      # Helper indices, lower imp index in ion arrays
    iicu=com.nhsp+imps*(species+1)      # Helper indices, upper imp index in ion arrays
    ingc=com.ngsp               # Helper index, index of neutral imps in gas array
    
    # Allocate arrays
    #- - - - - - - - -
    com.ngsp=com.ngsp+1         # Allocate arrays for impurity species in gaseous equatio   
    com.nzsp[species]=imps          # Number of impurity species for gas species species+1
                          # Determines nisp (nisp=nhsp+sum(nzsp)) and allocates arrays
    
    # Turn impurity gas on
    #- - - - - - - - - - - 
    bbb.istgcon[ingc]=0 # Impurity gas temperature: tg=(1-istgcon)*rtg2ti*ti+istgcon*tgas*ev

    # Impurity species parameters
    #- - - - - - - - - - - - - - -
    bbb.ziin[iicl:iicu]=range(1,imps+1)     # Impurity charge states
    bbb.minu[iicl:iicu]=12      # Atomic mass unit species mass
    bbb.znuclin[com.nhsp:com.nhsp+6]=6  # Nuclear charge of impurities


    allocate() # Reallocate arrays to accommodate for carbon

    """---------------------------------------------------------------------------------------------
    CARBON ION BC
    ---------------------------------------------------------------------------------------------"""

    # Carbon ion core BC
    #-----------------------------------------------------------------------------------------------
    bbb.isnicore[iicl:iicu]=3   # Core impurity ion BC model
                                    #=0: flux = curcore/sy locally in ix
                                    #=1: uniform, fixed density, ncore
                                    #=2: flux & ni over range
                                    #=3: icur=curcore-recycc*fngy, const ni
                                    #=4: impur. source terms (impur only)
                                    #=5: d(ni)/dy=-ni/lynicore at midp & ni constant poloidall
    if 1 in bbb.isnicore[iicl:iicu]:    # Constant imputiry dens on core vboundary
        bbb.ncore[iicl:iicu]=4.e15      # Core boundary density
    if 3 in bbb.isnicore[iicl:iicu]:    # Core BC set by fluxes
        bbb.curcore[iicl:iicu]=0        # Constant flux contribution over core BC

    bbb.isupcore[iicl:iicu]=1   # Velocity BC:
                                    #=0; Dirichlet BC: up=upcore
                                    #=1; Neumann BC: d(up)/dy=0
                                    #=2 sets d^2(up)/dy^2 = 0
                                    #=3 sets poloidal velocity (uu) = 0
                                    #=4 sets tor. ang mom flux = lzflux & n*uz/R=const
                                    #=5 sets ave tor vel = utorave & n*uz/R=const

    # Carbon ion wall BC
    #-----------------------------------------------------------------------------------------------
    # Outer wall
    #- - - - - - 
    bbb.isnwcono[iicl:iicu]=3# Outer wall BC:
                                #=0, old case; if ifluxni=0, dn/dy=0; if ifluxni=1, fniy=0 (default)
                                #=1, fixed density to nwallo(ix) array
                                #=2, extrapolation B.C.
                                #=3, approx gradient scale length, limited by nwomin
    if 3 in bbb.isnwcono[iicl:iicu]:    # Grad scale length BC: same as for plasma
        bbb.nwomin[iicl:iicu]=1.e7      # Minimum outer wall density
    
    # PFR wall
    #- - - - -
    bbb.isnwconi[iicl:iicu]=3# PFR wall BC: as above
    if 3 in bbb.isnwconi[iicl:iicu]:    # Grad scale length BC: same as for plasma
        bbb.nwimin[iicl:iicu]=1.e7      # Minimum PFR wall density

    # Carbon plate BC:s
    #-----------------------------------------------------------------------------------------------
    bbb.isbohmms=0      #0=single-species Bohm condition (H+)
                #1=multi-species Bohm condition (all ions)

    """---------------------------------------------------------------------------------------------
    NORMALIZATION AND BACKGROUND
    ---------------------------------------------------------------------------------------------"""
    bbb.n0[iicl:iicu]=1.e17 # Global impurity ion density normalization
    bbb.n0g[ingc]=1.e18     # Global impurity gas density normalization
    bbb.nzbackg=1.e10       # Background impurity ion density TODO fix indices?
    bbb.inzb=2          # Impurity floor scaling (nzbackg/ni)^inzb
    bbb.ngbackg[ingc]=1.e10 # Impurity gas background density 


    """---------------------------------------------------------------------------------------------
    CARBON RATES
    ---------------------------------------------------------------------------------------------"""
    # Setup rate model
    #- - - - - - - - -
    bbb.ismctab=2       # Define data to be used for multi-charge-state rates
                            #=1: tables originally generated by R. Campbell for D. Knoll,
                              # data file name is specified by inelmc=....
                              # corresponding rate evaluation routines are imprates and radimpmc.
                            #=2: tables generated by code from B. Braams,
                              # data file name is specified by mcfilename=...,
                              # corresponding rate evaluation routines are mcrates and radmc
    if bbb.ismctab==2:  # Braams tables
        com.mcfilename="C_rates.adas"   # Rate data to be used
        com.isrtndep=1          # Are table lookup parameters density dependent?
                            # Check compability with hydrogen if istabon=16
    
    # CX
    #- - -
    bbb.rcxighg[ingc]=0.0       # Turn off imp0 + H+ --> imp+ + H+
    
    com.iscxfit=2           # C-ion + H0 CX model:
                                #=0: analytic forms in Braams' rate package
                                #=1: polynomial fit to C.F. Maggi curves (1997)
                                #=2: same as =1, except Z=1 has lower rate from Pigarov

    # Scattering
    #- - - -  - -
    bbb.kelighi[ingc] = 5.e-16  # Elastic collision coefficient with H+
    bbb.kelighg[ingc] = 5.e-16  # Elastic collision coefficient with H0


    """---------------------------------------------------------------------------------------------
    CARBON SPUTTERING
    ---------------------------------------------------------------------------------------------"""
    
    # Chemical
    #- - - - - 
    bbb.isch_sput[ingc]=7   # Chemical sputtering model
                                #=0: Old
                                #=5: Roth, G-R
                                #=6: Haasz 97
                                #=7: Haasz 97 + Davis at low E

    bbb.fchemygwi=  1   # Factor multiplying chemical sputtering gas yield; PF wall
    bbb.fchemygwo=  1   # Factor multiplying chemical sputtering gas yield; Outer wall
    bbb.fchemylb=   1   # Factor multiplying chemical sputtering gas yield; Left plate
    bbb.fchemyrb=   1   # Factor multiplying chemical sputtering gas yield; Right plate




    # Physical
    #- - - - - - 
    bbb.isph_sput[ingc]=3   # Physical sputtering model
                                #=0: old fixed case
                                #=1: DIVIMP/JET physical sputtering fits
                                #=2: adds H+ chemical sputtering
                                #=3: adds H0 carbon sputtering
    bbb.crmb=bbb.minu[0]  # Mass of incident sputtering particles
    bbb.cizb=bbb.ziin[0]  # Max plasma charge state

    bbb.fphysylb=   1   # Factor multiplying physical sputtering gas yield; Left plate  
    bbb.fphysyrb=   1   # Factor multiplying physical sputtering gas yield; Right plate

    # Wall sputtering
    #- - - - - - 
    bbb.isi_sputw[ingc]=2   # Outer wall sputtering model
                                #=0: no ion sputtering
                                #=1: adds physical ion sputtering
                                #=2: adds chemical ion sputtering
    bbb.isi_sputpf[ingc]=2 # PF wall sputtering: as above
    
    bbb.t_wall=300      # Side wall temperatures
    bbb.t_plat=300      # Plate temperatures


def recycling_puff_pump():
    '''
    Function setting up recycling, puffing, and pumping in UEDGE
    '''

    """====================================================================================================
    WALL AND PLATE RECYCLING, PUFFING, PUMPING, ETC
    ===================================================================================================="""


    bbb.bcen=0          # Neutral energy transfer factor on plates
    bbb.bcenw=0         # Neutral energy transfer factor on walls

    """-----------------------------------------------------------------------------------------------------
    RECYCLING
    -----------------------------------------------------------------------------------------------------"""
    bbb.recycp[0]=0.98          # H0 plate recyc coeff
    bbb.recycw[0]=0.90          # H0 wall recyc coeff

    if bbb.ishymol==1:
        bbb.recycp[igh2]=1e-10      # H2 plate recyc coeff
        bbb.recycw[igh2]=1e-10      # H2 wall recyc coeff

    if bbb.isimpon>0:
        bbb.recycp[ingc]=1e-10        # C plate recyc coeff
        bbb.recycw[ingc]=1.e-10     # C wall recyc coeff (>0!)


    """-----------------------------------------------------------------------------------------------------
    ALBEDOS
    -----------------------------------------------------------------------------------------------------"""
    # Sets albedo-like pumping
    # Affects the one-sided maxwellian contribution to the flux
    # albedolb[gas spexcies, X-point]

    bbb.albedolb[0,:]=1       # Left plate H0 albedo
    bbb.albedorb[0,:]=1       # Right plate H0 albedo

    if bbb.ishymol==1:
        bbb.albedolb[igh2,:]=1   # Left plate H2 albedo
        bbb.albedorb[igh2,:]=1   # Right plate H2 albedo

    if bbb.isimpon>0:
        bbb.albedolb[ingc,:]=1.00   # Left plate C0 albedo
        bbb.albedorb[ingc,:]=1.00   # Right plate C0 albedo


    """-----------------------------------------------------------------------------------------------------
    HYDROGENIC STEADY-STATE PARTICLE BALANCE CALCULATION
    --------------------------------------------------------------------------------------------------------

    fnix(H0)=-0.01*ni(H0)*v_maxwellian(H0)*Area
    fngx(H0)=-0.01*ni(H0)*v_maxwellian(H0)*Area
    fngx(H2)=-0.5*( 0.01*fnix(H+) + fnix(H0) ) -> 0.01*fnix(H+)+fnix(H0) H-particles returned

    """






    """-----------------------------------------------------------------------------------------------------
    MOMENTUM RECYCLING
    -----------------------------------------------------------------------------------------------------"""
    # Momentum recycling for inertial H:
        # (-9.9,inf):   ydot=-nurlxu*( recycm*up(H+)+up(H0) ) / Normalization
        # (-10.1,-9.9]: Zero-Neumann BC
        # (-inf,10.1]:  Neutral thermal flux used
    bbb.recycm=-10     # Old comment: should be set to -0.9 TODO verify what should be used?




    """-----------------------------------------------------------------------------------------------------
    GAS PUFFS AND PUMPS
    -----------------------------------------------------------------------------------------------------"""
    bbb.nwsor=1     # Source regions at wall

    # Puffs: indices correspond to source # - limited by nwsor
    bbb.issorlb=    [1, 0, 0, 0, 0, 0, 0, 0, 0, 0]  # Source position origin flag; =0: right plate, =1: left plate
    # Outer wall
    bbb.igspsoro=   [1, 0, 0, 0, 0, 0, 0, 0, 0, 0]  # Puffed species index (PYTHON OR BASIS INDICES?)
    bbb.igaso=  [0, 0, 0, 0, 0, 0, 0, 0, 0, 0]  # Gas currents from outer wall [Amp]
    bbb.xgaso=  [0, 0, 0, 0, 0, 0, 0, 0, 0, 0]  # Location of wall source: origin set by issorlb
    bbb.wgaso=  [1e3, 0, 0, 0, 0, 0, 0, 0, 0, 0]  # Total Cosine widths of source
    bbb.albdso= [1.0, 1, 1, 1, 1, 1, 1, 1, 1, 1]  # Albedos at outer gas source locations (Pumps?)
    bbb.matwso= [1, 0, 0, 0, 0, 0, 0, 0, 0, 0]  # Material wall BC at outer gas puff on if =1 ???
    # Inner wall
    bbb.igspsori=   [1, 2, 0, 0, 0, 0, 0, 0, 0, 0]  # Puffed species index (PYTHON OR BASIS INDICES?)
    bbb.igasi=  [0, 0, 0, 0, 0, 0, 0, 0, 0, 0]  # Gas currents from inner wall [Amp]
    bbb.xgasi=  [0, 0, 0, 0, 0, 0, 0, 0, 0, 0]  # Location of wall source: origin set by issorlb
    bbb.wgasi=  [1e3, 1e3, 0, 0, 0, 0, 0, 0, 0, 0]  # Total Cosine widths of source
    bbb.albdsi= [1.0, 1.0, 1, 1, 1, 1, 1, 1, 1, 1]  # Albedos at inner gas source locations (Pumps?)
    bbb.matwsi= [1, 1, 0, 0, 0, 0, 0, 0, 0, 0]  # Material wall BC at inner gas puff on if =1 ???





    if bbb.isimpon>0:
        """---------------------------------------------------------------------------------------------
        10.6: Impurity ion wall puffing
        ---------------------------------------------------------------------------------------------"""
        api.nzsor=0                 # N.o. impurity wall sources
        if api.nzsor>0:
            api.iszsorlb=   [0, 0, 0, 0, 0, 0]      # Source location origin: 1=lb, 0=rb
            api.wimpi=  [1e3, 1e3, 1e3, 1e3, 1e3, 1e3]  # PF source zone width
            api.wimpo=  [1e3, 1e3, 1e3, 1e3, 1e3, 1e3]  # Outer wall source zone width
            api.impsori=    [0, 0, 0, 0, 0, 0]      # PF wall source strength [Amp]
            api.impsoro=    [0, 0, 0, 0, 0, 0]      #  wall source strength [Amp]


def transport_fluxlim():
    '''
    Function defining transport parameters and flux limiting
    '''

    """---------------------------------------------------------------------------------------------------- 
    TRANSPORT PARAMETERS
    -----------------------------------------------------------------------------------------------------"""

    #TODO place these
    bbb.sxgsol=1.0  # X-coord stretching factors in SOL & core for neutral diffusion
    bbb.sxgpr=1.0   # X-coord stretching factors in PFR for neutral diffusion

    # 8.1.1: Parallel transport
    #- - - - - - - - - -
    bbb.kxe=1.35      # Poloidal electron heat condictivity multiplication factor
    bbb.kxi=1.      # Poloidal ion heat conductivity multiplication factor

    # 8.1.4: Anomalous radial transport
    #- - - - - - - - - - - - - - - - - 

    bbb.isbohmcalc=1    # Bohm condition switch:
                            #=0: Uses the values set in dif_use, dif2_use,
                              # tra_use, kye_use, kyi_use, and difutm.
                            #=1: calc Bohm diff is facb...>0
                            #=2: Harmonic average of Bohm, difni, etc
                            #=3: D=difniv*(B0/B)**inbpdif, etc

    if bbb.isbohmcalc==0:   # Setting all diffusivities manually: Needs allocation and manual operation
                            # Not used now: easier to use isbohmcalc=3, allocate and populate arrs,
                            # and edit after an exmain
        bbb.fcdif=0.0       # Scaling factor for constant nomalous diffusivities (fcdif*kye+kye_use, etc)
                            # =0: use only *_use parameters
        bbb.dif_use[:,:,0]=1
        bbb.dif_use[:,:,1]=0

        bbb.tray_use=2.6
        bbb.ky_use=1
        bbb.kyi_use=1

        if bbb.isimpon>0:
            bbb.dif_use[:,:,iicl:iicu]=1

    if bbb.isbohmcalc==1:
        bbb.fcdif=1.0       # Scaling factor for constant nomalous diffusivities
        bbb.difni[0]=1
        bbb.kye=1
        bbb.kyi=1
        bbb.travis[0]=1
        bbb.travis[1]=0
        bbb.parvis[0]=1





    if bbb.isbohmcalc==2:   # Use harmonic ave of difni etc and bohm
        bbb.kye=2.6     # Radial electron heat diffusivity
        bbb.kyi=0.75        # Radial ion heat diffusivity

        # H+
        bbb.difni[0]=2.6          # Radial H+ density diffusivity cefficient
        bbb.travis[0]=2.6         # H+ perpendicular viscosity

        # H0
        bbb.difni[1]=0.5         # Radial H+ density diffusivity cefficient
        bbb.travis[1]=1          # H+ perpendicular viscosity

        # C
        if bbb.isimpon>0:
            bbb.difni[iicl:iicu]=2.6    # Radial H+ density diffusivity cefficient
            bbb.travis[0]=1       # H+ perpendicular viscosity

    if bbb.isbohmcalc==3:   # B-scaling
        bbb.inbtdif=0       # Exponential scaling of (B0/B)
        bbb.inbpdif=0       # Exponential scaling of (B0/B)

        bbb.kyev=1          # Radial electron heat diffusivity
        bbb.kyiv=0.75       # Radial ion heat diffusivity

        # H+
        bbb.difniv[:,0]=1           # Radial H+ diffusivity
        bbb.travisv[:,0]=2.6

        # H0
        bbb.difniv[:,1]=0            # Radial H0 diffusivity TODO Ceck if this one is indeed intended to be 0?
        bbb.travisv[:,1]=2.6

        # C
        if bbb.isimpon>0:
            for ind in range(iicl,iicu):
                bbb.difniv[:,ind]=1   # Radial C diffusivity
            bbb.travisv[:,iicl:iicu]=2.6    # Perpendicular viscosity

    """---------------------------------------------------------------------------------------------------- 
    FLUX LIMIT FACTORS
    ----------------------------------------------------------------------------------------------------""" 
    bbb.flalfv=1.0      # Parallel velocuty flux limit factor (linear)
    bbb.flalfe=0.21     # Parallel electron heat flux limit factor: 1/(1+abs(1/flalfe)**flgam)
    bbb.flalfi=0.21     # Parallel ion heat flux limit factor: 1/(1+abs(1/flalfi))
    bbb.flalfgx=1    # Poloidal gas diffusivity flux limit   
    bbb.flalfgy=1    # Radial gas diffusivity flux limit
    bbb.flalfgxy=1   # Nonorthogonal poloidal face gas flux limit
    bbb.flalftgx=1   # Poloidal gas temperature diffusivity flux limit
    bbb.flalftgy=1   # Radial gas temperature diffusivity flux limit
    bbb.flalftmx=1   # Poloidal molceular temperature diffusivity flux lim
    bbb.flalftmy=1   # Radial molceular temperature diffusivity flux lim
    bbb.lgmax=1e20       # Max gas scale length for calculating particle gaseous diffusivity coefficient
    bbb.lgtmax=1e20      # Max gas scale length for calculating thermal gas diffusivity coefficient
    bbb.lgvmax=1e20      # Max gas scale length for calculating viscous gas diffusivity coefficient
def currpot():
    '''  
    Function setting up the currents and potentials in UEDGE
    '''

    """-----------------------------------------------------------------------------------------------------
    CURRENT AND POTENTIAL PARAMETERS
    -----------------------------------------------------------------------------------------------------"""

    # SHEAT PARAMETERS
    #- - - - - - - - - - - - - 
    # TODO see if active
    bbb.newbcl=1        # Linear scaling factor for left plate sheath potential:
                            #=0: use bcei and bcee
                            #=1: new model
    bbb.newbcr=1        # Linear scaling factor for right plate sheath potential:
                            #=0: use bcei and bcee
                            #=1: new model
    bbb.bcei=2.5    # Ion sheath energy transmission factor
    bbb.bcee=4.0    # Electron sheath energy transmission factor


    # MAGNETIC FIELDS
    #- - - - - - - - - - - -
    bbb.b0=1    # Magnetic field scale factor (Buse=B/B0)
                    #>0: normal direction B field
                    #<0: rev B-field


    bbb.rsigpl=1.e-8    # Anomalous radial electrical conductivity
    bbb.cfjhf=1.        # Coefficient for convective current (fqp) heat flow
    bbb.cfjve=1.        # Coefficient for current contribution to electron velocity: vex=vix-cfjve*fqx
    bbb.cfjpy = 0       # Coefficient for diamagnetic drift in y-direction
    bbb.cfjp2 = 0       # Coefficient for diamagnetic drift in direction perpendicular to radial and parallel directions
    bbb.jhswitch=1      #Joule Heating term coefficient

    # DRIFTS
    #- - - - - - - - 
    # Diamagnetic
    bbb.isfdiax=1       # Diamagnetic drift sheat contribution coefficient
    bbb.cfqydbo=0       # Coefficient for diagmagetic current on core boundary (=1 forces j_r=0)
    bbb.cfydd=0.0       # Divergence-free diamagnetic drift coefficient in radial direction
    bbb.cf2dd=0.0       # Diamagnetic drift coefficient in direction perpendicular to parallel and radial direction
    bbb.cftdd=0.0       # Coefficient for diamagnetic contribution on toroiudal drift

    # ExB
    bbb.cfyef=1.0       # ExB drift coefficient in y-direction
    bbb.cftef=0.0       # ExB drift coefficient in toroidal direction
    bbb.cf2ef=1.0       # ExB drift in direction perpendicular to radial and parallel direction

    # Grad(B)
    bbb.cfybf=1.0       # Grad(B) radial drift coefficient 
    bbb.cf2bf=1.0       # Grad(B) drift in direction perpendicular no parallel and radial direcion
    bbb.cfqybf=1.0      # Grad(B) radial current coefficient 
    bbb.cfq2bf=1.0      # Grad(B) current in direction perpendicular no parallel and radial direcion
    bbb.cfqybbo=0       # Grad(B) current coefficient on core boundary
    bbb.cfniybbo=0.     # Grad(B) coefficient of velocity contribution to ion/neutral flux on core boundary
    bbb.cfeeybbo=0.     # Grad(B) coefficient of velocity contribution to electron flux on core boundary

    # Grad(P_i x B)
    bbb.cfniydbo=0.     # Grad(P_i x B) coefficient of velocity contribution to ion/neutral flux on core boundary
    bbb.cfeeydbo=0.     # Grad(P_i x B) coefficient of velocity contribution to electron flux on core boundary
    # BxGrad(T)
    bbb.cfeixdbo=0.     # Coeffcient for BxGrad[T] drift in plate BC
    bbb.cfeexdbo=0.     # Coefficient for diamagnetic drift in plate BC
    # Inertial correction
    bbb.cfqym=1.0       # Coefficient for spatial inertial radial current in Y-dir (fqy)



    # POTENTIAL
    #- - - - - - - - -
    bbb.isnewpot=1      # Potential model switch:
                            #=1: new potential model; radial curr dens from toroidal momentum balance
                            #=-2: Constant Phi on core boundary with total core current = icoreelec
    bbb.rnewpot=1.0     # Linear scaling factor for potential model: 1=new, 0=old

    if bbb.isnewpot==1: # New potential equation
        bbb.iphibcc=0       # Core potential boundary equation:
                                #=1: d^2(ey)/dy^2=0
                                #=2: te=constant & ey(ixmp,0)=eycore
                                #=3:  phi=constant & ey(ixmp,0)=eycore
                                #else:  dphi(ix,1)=dphi_iy1; isutcore controls ix=ixmp 
        if bbb.iphibcc not in [1,2,3]:  # No extrapolation BC for radial electric field on core boundary
            bbb.isutcore=2  # Model for determining Phi at midplane
                                #=0 toroidal momentum = lzcore on core
                                #=1 d<uz>/dy=0
                                #>1 d^2(Ey)/dy^2=0 at OMP

    bbb.iphibcwi=0      # Potential BC at PF wall:
                            #=0: d(ey)/dy=0
                            #=1: phi(ix,0) = phintewi*te(ix,0)/ev
                            #=3: d(phi)/dy/phi = 1/lyphi(1)
                            #=4: phi(ix,0)=phiwi(ix) in PF region
    bbb.iphibcwo=0      # Potential BC at outer wall: as above
def equations():
    '''
    Function specifying the equations to be solved by UEDGE
    '''

    """====================================================================================================
    DEFINE EQUATIONS TO BE SOLVED
    ===================================================================================================="""

    """----------------------------------------------------------------------------------------------------
    CONTINUITY EQUATION
    ----------------------------------------------------------------------------------------------------"""
    bbb.isnion[0]=      1   # H+
    bbb.isnion[1]=      1   # Inertial H0
    bbb.isngon[0]=      0   # =0 if isupgon[0]=1: Inertial H0 model!
    bbb.isngon[1]=      0   # Diffusive H2
    bbb.isngon[2]=      0   # Diffusive neutral C
    bbb.isnion[2:8]=    [0, 0, 0, 0, 0, 0]  # C ions


    # H+ continuity equation scalings
    #- - - - - - - - - - - - - - - - - -
    bbb.cnfx=1. # Coefficient for poloidal convection in ion continuity equation
    bbb.cnfy=1. # Coefficient for radial convection in ion continuity equation
    bbb.cnsor=1.    # Coefficient for particle source in ion continuity equation
    bbb.cnurn=1.    # Coefficient scaling nurlx (rate coefficient to relax BC) in ion continuity equation

    # H0 continuity equation scalings
    #- - - - - - - - - - - - - - - - - -
    bbb.cngsor=1.   # Coefficient for particle source in gaseous continuity equation
    bbb.cnurg=1.    # Coefficient scaling nurlx (rate coefficient to relax BC) in gaseous continuity equation


    """----------------------------------------------------------------------------------------------------
    MOMENTUM EQUATIONS
    ----------------------------------------------------------------------------------------------------"""
    bbb.isupon[0]=              1   # H+
    bbb.isupon[1]=              1   # H0
    bbb.isupgon[0]=             1   # H0: INERTIAL FLAG
    bbb.isupgon[1]=             0   # H2
    bbb.isupgon[2]=             0   # C
    bbb.isupon[2:8]=            [0, 0, 0, 0, 0, 0]  # C ions

    # Momentum equation options
    #- - - - - - - - - - - - - - - - - 
    bbb.cmfx=1. # Coefficient for poloidal convection in ion momentum equation
    bbb.cmfy=1. # Coefficient for radial convection in ion momentum equation
    bbb.cpgx=1. # Coefficient for pressure gradient in ion momentum equation
    bbb.cnuru=1.    # Coefficient scaling nurlx (rate coefficient to relax BC) in ion momentum equation

    # Check that the correct model is chosen
    #- - - - - - - - - - - - - - - - - 
    if bbb.isupgon[0]==1 and bbb.isngon[0]==1:
        print('WARNING! Both the switch for inertial atoms (insupgon[0]=1) and diffusive atoms (isngon[0]=1) are turned on!')
        print('UEDGE only allows one switch to be active as the options are conflicting. Aborting read...')
        return 'Failure to read'
        


    """----------------------------------------------------------------------------------------------------
    ENERGY EQUATION
    ----------------------------------------------------------------------------------------------------"""
    bbb.isteon=                 1   # Electrons
    bbb.istion=                 1   # H+
    bbb.istgon[1]=              0   # H2

    # Energy equation options
    #- - - - - - - - - - - - - 
    bbb.cfvisx=1.   # Coefficient for poloidal viscosity in ion energy equation
    bbb.cfvisy=1.   # Coefficient for radia viscosity in ion energy equation
    bbb.cvgp=1. # Coefficient for v.Grad(P) in ion energy equation
    bbb.cnure=1.    # Coefficient scaling nurlx (rate coefficient to relax BC) in electron energy equation
    bbb.cnuri=1.    # Coefficient scaling nurlx (rate coefficient to relax BC) in ion energy equation

    """----------------------------------------------------------------------------------------------------
    POTENTIAL EQUATION
    ----------------------------------------------------------------------------------------------------"""
    bbb.isphion=                1   # Potential equation
    bbb.isphiofft=              0      # Potential eq complement
def restore(restart):
    '''
    Function to restore UEDGE save file if argument is True
    '''

    if restart==True:
        bbb.restart=1       # Flag to restart from previous case (=1)
    if bbb.restart==1:      # Restart from previous case
        hdf5_restore("../solutions/"+bbb.label[0].decode('UTF-8')+".hdf5")    # Append hdf5 restart name and relative pat and relative pathh
        


