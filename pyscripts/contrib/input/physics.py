# Physics package for UEDGE, physics.py
# Separated from 160299 master input by holm10
# Changelog
# 200213 - file created
from uedge import bbb,com,grd,flx,aph


def powerdens_core(pcoree,pcorei,ncore,aphpath='.',lyni=[0.05,0.05],nwomin=1e15,lyte=[0.05,0.05],lyti=[0.05,0.05],nwimin=1e16,isplflxl=1,ngbackg=1e13,cngflox=1,cngfloy=1,
                    isngcore=1,albedoc=0.5,isupcore=1,xstscal=0.02,ngscal=1,xgscal=0.01,cfloxiplt=1

):
    """ Plasma setup with power and density BC at core
    plasma_Pcore_ncore(pcoree,pcorei,ncore,aphpath)

    Function setting up the UEDGE plasma with core-flux power
    boundary conditions for ions and electrons and constant 
    core boundary plasma density boundary condition.

    Keyword parameters:
    pcoree -  electron power crossing the core boundary
    pcorei -  ion power crossing the core boundary
    ncore -   plasma density at the core boundary
    aphpath - path to hydrogenic rate filestoms
    """
    # Set path to rate data
    aph.aphdir=aphpath # Hydrogen rates

    # Initalize arrays
    com.nhsp=2          # N.o. hydrogenic species
    com.ngsp=1          # N.o. hydrogenic gaseous species


    ''' Charged species setup '''
    bbb.minu[0]=2     # H+ mass in AMU
    bbb.ziin[0]=1     # H+
    bbb.znuclin[0]=1  # H+ nuclear charge 

    # Scale factors
    #- - - - - - - - - - -
    bbb.cngtgx=0.   # X-flux coefficient for gaseous component i in ion energy equation
    bbb.cngtgy=0.   # Y-flux coefficient for gaseous component i in ion energy equation

    corepower(pcoree,pcorei)    # Set power BC for ions at core bound
    coremomentum(isupcore)     # Set neumann momentum BC for plasma at core bound
    coredens(ncore)             # Set dirichlet BC for plasma at core bound
    
    wall_BC_scalelength(lyni=lyni,nwomin=nwomin,lyte=lyte,lyti=lyti,nwimin=nwimin,isplflxl=isplflxl)       # Set gradient scal-length wall BCs

    ''' Atom model setup '''
    inertial_atoms( ngbackg=ngbackg,cngflox=cngflox,cngfloy=cngfloy,isngcore=isngcore,albedoc=albedoc, 
                    isupcore=isupcore,xstscal=xstscal,ngscal=ngscal,xgscal=xgscal,cfloxiplt=cfloxiplt)    # Activate intertial atoms 

    ''' Rate model setup '''
    Stotler_loglog()    # Use Stotler log-log rate file

    ''' Allocate plasma arrays '''
    if bbb.pyrestart_file[0].decode('UTF-8')!='read':
        bbb.allocate()          



def Stotler_loglog(icnucx=0,cnucx=0,icnuiz=0,cnuiz=5e4,isrecmon=1,cfrecom=1,eion=5,ediss=10):
    """---------------------------------------------------------------------------------------------------- 
    RATE PARAMETERS
    -----------------------------------------------------------------------------------------------------"""

    # LOOK-UP TABLES FOR HYDROGENIC RATES
    com.istabon = 10    
    ''' Rate model/lookup tables:
        0: simple analytic rates and constant energy loss per ionization
        1: table look-up from ADPAK; rates.adpak
        2: table look-up from STRAHL; rates.strahl
        3: table look-up old DEGAS (created 84/09/14); eh.dat & atmc.dat
        4,5,6: table look-up from new DEGAS (created 93/05/06): nwfits
            4:  linear interpolation for  rsa  vs log(te) and log10(ne)
            5:  spline fit for log10(rsa) vs log(te) and log10(ne) 
            temp=6:  Hindmarsh spline    log10(rsa) vs log(te) and log10(ne) 
            disabled=6_  spline fit for       rsa  vs log(te) and log10(ne)
        7: Campbell's poly. fit for rsa, etc., vs log10(te) and log10(ne)
        8: tab look-up from latest DEGAS (created 93/04/08); ehr1.dat
        9: tab look-up; Stotler PPPL (~95/07/10) with log(Te)-sigv; ehr2.dat
        10: tab look-up;Stotler PPPL (~95/07/10);log(Te)-log(sigv); ehr2.dat
        11: same as istabon=10 with data for n=5-9 excited states; thin.dat
        12: as istabon=11, ex. Lyman-alpha local absb (thick); thickLyA.dat
        13: as istabon=11, ex. all Lyman lines loc absorbed; thickAllLy.dat
        14: H.Scott data; add rtau=Ly-a opacity, lin. interp; ehrtau.dat
        15: H.Scott data; add rtau=Ly-a opacity, log interp; ehrtau.dat
        16: table look-up using 'b2frates_hyd' file data, log-log interp
    '''

    if com.istabon==16: # Table lookup used
        com.isrtndep=1      # Are table lookup parameters density dependent: check compabnility with imps!

    # CHARGE-EXCHANGE
    #- - - - - - - - - - -
    bbb.icnucx=icnucx    # CX rate model:
                        #=0: variable nucx
                        #=1: constant nucx=cnucx
                        #=2: use sigcx; nucx~sqrt(Tg)
    if bbb.icnucx==1:   # Constant CX rates
        bbb.cnucx=cnucx        # Defined constant CX rate



    # IONIZATION
    #- - - - - - - - - - 
    bbb.icnuiz=icnuiz    # Ionization model:
                        #=0: variable nuiz
                        #=1: constant nuiz=cnuiz
                        #=2: Freezes
    if bbb.icnuiz==1:   # Constant ionzation frequenzy
        bbb.cnuiz=cnuiz      # Defined nuix

    # RECOMBINATION
    #- - - - - - - - - - -
    bbb.isrecmon = isrecmon    # Switch for recombination: 1=on, 0=off
    if bbb.isrecmon==1: # If recombination on
        bbb.cfrecom=cfrecom       # Coefficient for recombination frequencies

    # DISSOCIATION
    #- - - - - - - -
    bbb.eion = eion            # Energy per atom from dissociation
    bbb.ediss = ediss    # Electron dissociation loss: 
                                    #<2*eion, virtual energy source through dissociation
                                    #=2*eion, no radiation and power conserved
                                    #>2*eion, virtual energy sink through dissociation (radiation)







def inertial_atoms( istgcon=0,tgas=0,tgcore=100,istgcore=1,cfloxiplt=1, isngcore=1,albedoc=1,
                    cngfx=1,cngfy=1,cngflox=1,cngfloy=1,cngmom=0,cmwall=0,cfbgt=0,kxn=0,kyn=0,ngbackg=1e13,
                   isupcore=1,xstscal=0.02,ngscal=1,xgscal=0.01):
    """-----------------------------------------------------------------------------------------------------
    ATOMIC SPECIES SETUP
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
    bbb.istgcon[0]=istgcon # H0 temperature: tg=(1-istgcon)*rtg2ti*ti+istgcon*tgas*ev
    bbb.tgas[0]=tgas    # Neutral temperature scaling is istgcon>1



    # 5.2.1: CORE BC
    #- - - - - - - - - 
    # Density
    bbb.isnicore[1]=0    # Switch off neutral ion BC:s 
    bbb.isngcore[0]=isngcore    # Neutral gas core boundary conditions
                                #=0, set loc flux= -(1-albedoc)*ng*vtg/4
                                #=1, set uniform, fixed density, ngcore
                                #=2, set rad. grad. to sqrt(lam_i*lam_cx)
                                #=3, extrapolation, but limited
                                #=anything else, set zero deriv which was
                                # prev default inert hy
                                # anything else same as =0
    if bbb.isngcore[0]==0: 
        bbb.albedoc[0]=albedoc
    bbb.ngcore[0]=2.e13  # Core H0 density

    # Momentum
    bbb.isupcore[1]=isupcore    # Velocity BC:
                                #=0; Dirichlet BC: up=upcore
                                #=1; Neumann BC: d(up)/dy=0
                                #=2 sets d^2(up)/dy^2 = 0
                                #=3 sets poloidal velocity (uu) = 0
                                #=4 sets tor. ang mom flux = lzflux & n*uz/R=const

    # Energy 
    bbb.istgcore[0] = istgcore  # H0 temperature core BC:
                                #=0; set core boundary temperature to ion temperature
                                #=1; set core boundary temperature to tgcore
    bbb.tgcore[0]=tgcore    # Core D0 boundary temperature


    # Wall boundary conditions
    #- - - - - - - - - - - - - - - - - 
    bbb.tgwall[0] = 4.e-2      # Wall gas temperature when BC used, atoms
    # DNESITY
    # PF wall
    bbb.isnwconi[1]=0    # PF wall density BC as above
    bbb.ifluxni=1       # Switch for setting PF and outer wall fluxes to 0 (=1)

    # Outer wall
    bbb.isnwcono[1]=0    # Outer wall density BC as above
    bbb.ifluxni=1       # Switch for setting PF and outer wall fluxes to 0 (=1)

    # Plate BC
    #- - - - -
    bbb.cfloxiplt=cfloxiplt         # Coefficient considering the neutral convected energy from plates: 0=off, 1=on

    # Background gas and normalization
    #- - - - - - - - - - - - - - - - - - - - - 
    bbb.ingb=2      # BG gas source scaling: BG gas source=nuiz*ngbackg*(0.9+0.1(ngbackg/ng)**ingb)
    bbb.ngbackg[0]=ngbackg   # "soft" artificial floor for neutral densities


    # Scale factors
    #- - - - - - - - - - - 

    bbb.cngfx=cngfx    # Scale factor for flux from Grad_x temperature in neutral gas velocity eqn
    bbb.cngfy=cngfy    # Scale factor for flux from Grad_y temperature in neutral gas velocity eqn

    bbb.cngflox[0]=cngflox  # Factor for x-flux from convection in atom continuity equation
    bbb.cngfloy[0]=cngfloy  # Factor for y-flux from convection in atom continuity equation


    bbb.xstscal=.02 # Exponential scale-length with stretched coordinate decays from plates
    bbb.ngscal=1.   # Ratio of initial gas density to ion density (restart=0)
    bbb.xgscal=.01  # Exponential scale of initial gas (restart=0)



    # Inertial neutrals
    bbb.cngmom=cngmom    # Momentum loss coefficient for diffusive H0 # Defined above
    bbb.cmwall=cmwall    # Momentum wall coefficient for neutral hydrogen only # Defined above
    bbb.cfbgt=cfbgt # B x Grad(T) coefficient
    bbb.kxn=kxn   # Poloidal CX-neutral heat diffusivity factor
    bbb.kyn=kyn   # Radial CX-neutral heat diffusivity factor






def corepower(pcoree,pcorei):
    ''' Sets the core energy BC to e/i power over bound '''

    bbb.iflcore=1   # Core power condition:
                        #=0; Temperature dirichlet BC: core Te,i=tcoree,i
                        #=1; Power dirichlet BC:  core power=pcoree,i
                        #=-1; Temperature Neumann BC:  core d(Te,i)/dy=0 

    bbb.pcoree=pcoree   # Electron power over core boundary
    bbb.pcorei=pcorei   # Ion power over core boundary

def coretemp(tcoree,tcorei):
    ''' Sets the core energy BC to e/i temp at bound '''

    bbb.iflcore=0   # Core power condition:
                        #=0; Temperature dirichlet BC: core Te,i=tcoree,i
                        #=1; Power dirichlet BC:  core power=pcoree,i
                        #=-1; Temperature Neumann BC:  core d(Te,i)/dy=0 

    bbb.tcoree=tcoree   # Electron temperature at core boundary
    bbb.tcorei=tcorei   # Ion temperature at core boundary


def coredens(ncore):
    ''' Set density at core bound to ncore '''
    # Core density BC
    #- - - - - - - - - - - - - -
    bbb.isnicore[0]=1    # Density BC:
                                #=0, set flux to curcore/sy locally in ix
                                #=1, set uniform, fixed density, ncore
                                #=2, set flux & ni over range
                                #=3, set icur=curcore-recycc*fngy, const ni
                                #=4, use impur. source terms (impur only)
                                #=5, set d(ni)/dy=-ni/lynicore          
    bbb.ncore[0]=ncore   # H+ core density [m^-3]



def coremomentum(isupcore):
    ''' Set Neumann momentum BC for plasma at core bound '''
    # Core momentum BC
    #- - - - - - - - - - - - - -
    bbb.isupcore[0]=isupcore # Velocity BC:
                            #=0; Dirichlet BC: up=upcore
                            #=1; Neumann BC: d(up)/dy=0
                            #=2 sets d^2(up)/dy^2 = 0
                            #=3 sets poloidal velocity (uu) = 0
                            #=4 sets tor. ang mom flux = lzflux & n*uz/R=const
                            #=5 sets ave tor vel = utorave & n*uz/R=const

def wall_BC_scalelength(lyni=[0.05,0.05],nwomin=1e15,lyte=[0.05,0.05],lyti=[0.05,0.05],nwimin=1e16,isplflxl=1
):

    # PF wall density BC
    #- - - - - - - - - - - - - - - 
    bbb.isextrnpf=0 # Extrapolation BC:
                        #=0; no extrapolation for ni at PF wall
                        #=1; extrapolation BC for ni at PF wall

    bbb.isnwconi[0]=3 # PF wall density BC:
                            #=0, old case; if ifluxni=0, dn/dy=0; if ifluxni=1, fniy=0 (default)
                            #=1, fixed density to nwallo(ix) array
                            #=2, extrapolation B.C.
                            #=3, approx gradient scale length
    bbb.lyni[0]=lyni[0]    # Fixed scale length at PF wall for BC
    bbb.nwimin[0]=nwimin   # Minimum density limit at PF boundary

    # Outer wall density BC
    #- - - - - - - - - - - - - - - - 
    bbb.isextrnw=0  # Extrapolation BC:
                        #=0; no extrapolation for ni at outer wall
                        #=1; extrapolation BC for ni at outer wall
    bbb.isnwcono[0]=3 # Outer wall density BC:
                            #=0, old case; if ifluxni=0, dn/dy=0; if ifluxni=1, fniy=0 (default)
                            #=1, fixed density to nwallo(ix) array
                            #=2, extrapolation B.C.
                            #=3, approx gradient scale length, limited by nwomin
    bbb.lyni[1]=lyni[1]    # Fixed scale length at outer wall for BC
    bbb.nwomin[0]=nwomin   # Minimum density limit at outer boundary



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
    bbb.lyte[0]=lyte[0]# Outer wall electron temperature scale length

    # PF wall ion  energy BC
    #- - - - - - - - - - - - - - - - -
    bbb.istipfc=3   # PF wall H+ temperature BC as above:
    bbb.lyti[0]=lyti[0]# PF wall H+ temperature scale length


    # Outer wall electron  energy BC
    #- - - - - - - - - - - - - - - - - - - - - 
    bbb.istewc=3    # Outer wall e- temperature BC:
                        #=0, zero energy flux
                        #=1, fixed temp to tedge or tewallo
                        #=2, extrapolation BC
                        #=3, Te scale length
                        #=4, feey = ~bceew*te*elec_flux
    bbb.lyte[1]=lyte[1]# Outer wall electron temperature scale length

    # Outer wall ion  energy BC
    #- - - - - - - - - - - - - - - - - - - -
    bbb.istiwc=3    # Outer wall H+ temperature BC as above
    bbb.lyti[1]=lyti[1]# Outer wall H+ temperature scale length
        

    # Plate BC
    #- - - - - - - - - - 
    bbb.isupss=1    # Plate boundary condition
                        #=-1: dup/dx=0
                        #=0: up=cs
                        #=1: up>=1
    # Flux limits
    bbb.isplflxl=isplflxl  # Switch activating flux limits (flalfe/flalfi) at plates (=1)

def volsource(ivolcur,zwni=1000,rwni=1000):
    """---------------------------------------------------------------------------------------------------- 
    VOLUMETRIC PLASMA SOURCES
    ----------------------------------------------------------------------------------------------------"""
    # Volumetric sources
    bbb.ivolcur=ivolcur     # Volume source [A] for EACH charge species
    bbb.zwni=zwni      # Width for volume source
    bbb.rwni=rwni      # Width for volume source





