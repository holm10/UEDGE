# Molceule package for UEDGE, molecules.py
# Separated from 160299 master input by holm10
# Changelog
# 200213 - file created
from uedge import bbb,com,grd,flx


def activate_mol():
    ''' Turns on the molecular switches, returns the mol gas ind '''
    bbb.ishymol=1       # Includes molecules as 2nd gaseous species (index 1)

    com.nhgsp=com.nhgsp+1       # Allocate space for molecules in hygrogenic gas species array
    com.ngsp=com.ngsp+1     # Allocate space for hydrogen in gas species array
    return com.ngsp-1         # Index for molecules in gaseous arrays
    


def constantTm(Tm,n0g=1e17,ngbackg=1e11,kelighg=5e-16,kelighi=5e-16,cngfx=0,cngfy=0,cfcvtg=1,cftgcond=1,ismolcrm=1):

    """====================================================================================================
    MOLECULAR HYDROGEN SETUP FOR MOLECULES WITH SPATIALLY CONSTANT TEMPERAUR
    ===================================================================================================="""
    igh2=activate_mol()

    bbb.istgcon[igh2] = 1.         # Don't reset tg[,,2] using istgcon
    bbb.tgas[1]=Tm

    # Set parameters common to all mol models
    common_mol(igh2,n0g,ngbackg,kelighg,kelighi,cngfx,cngfy,cfcvtg,cftgcond,ismolcrm)


    """---------------------------------------------------------------------------------------------------- 
    ALLOCATE MOLECULAR ARRAYS
    ----------------------------------------------------------------------------------------------------"""
    if bbb.pyrestart_file[0].decode('UTF-8')!='read':
        bbb.allocate()          






def Emol_V707( n0g=1e17,ngbackg=1e11,kelighg=5e-16,kelighi=5e-16,cngfx=0,cngfy=0,cfcvtg=1,cftgcond=1, 
                    isngcore=1,albedoc=0.5,ngcore=1e12,istgcore=1,tgcore=100,tgwall=4e-2,ismolcrm=1):

    igh2=activate_mol()

    bbb.istgcon[igh2] = -1.         # Don't reset tg[,,2] using istgcon


    # Set the energy boundary conditions
    energy_bc_mol(igh2,isngcore,albedoc,ngcore,istgcore,tgcore,tgwall)

    # Set parameters common to all mol models
    common_mol(igh2,n0g,ngbackg,kelighg,kelighi,cngfx,cngfy,cfcvtg,cftgcond,ismolcrm)

    """---------------------------------------------------------------------------------------------------- 
    ALLOCATE MOLECULAR ARRAYS
    ----------------------------------------------------------------------------------------------------"""
    if bbb.pyrestart_file[0].decode('UTF-8')!='read':
        bbb.allocate()          




def Emol(  n0g=1e17,ngbackg=1e10,kelighg=5e-16,kelighi=5e-16,cngfx=1,cngfy=1,cfcvtg=1,cftgcond=1,
                isngcore=0,albedoc=0.5,ngcore=1e12,istgcore=2,tgcore=100,tgwall=4e-2,ismolcrm=1):


    igh2=activate_mol()

    bbb.istgcon[igh2] = -1.         # Don't reset tg[,,2] using istgcon

    # Set the energy boundary conditions
    energy_bc_mol(igh2,isngcore,albedoc,ngcore,istgcore,tgcore,tgwall)


    # Set parameters common to all mol models
    common_mol(igh2,n0g,ngbackg,kelighg,kelighi,cngfx,cngfy,cfcvtg,cftgcond,ismolcrm)

    """---------------------------------------------------------------------------------------------------- 
    ALLOCATE MOLECULAR ARRAYS
    ----------------------------------------------------------------------------------------------------"""
    if bbb.pyrestart_file[0].decode('UTF-8')!='read':
        bbb.allocate()          


def energy_bc_mol(igh2,isngcore,albedoc,ngcore,istgcore,tgcore,tgwall):

    ''' Core BC '''
    # Desntiy
    bbb.isngcore[igh2] = isngcore  # Hydrogen molecular core BC:
                                #=0, set loc flux= -(1-albedoc)*ng*vtg/4
                                #=1, set uniform, fixed density, ngcore
                                #=2, set rad. grad. to sqrt(lam_i*lam_cx)
                                #=3, extrapolation, but limited
                                #=anything else, set zero deriv which was
                                # prev default inert hy
                                # anything else same as =0
    if bbb.isngcore[igh2]==0:   # Local core flux pumped
        bbb.albedoc[igh2]=albedoc       # Core H2 albedo
    elif bbb.isngcore[igh2]==1: # Uniform core H2 density
        bbb.ngcore[igh2] = ngcore    # Core H2 density

    # 6.1.1: Energy
    #- - - - - - - - 
    bbb.istgcore[igh2] = istgcore  # H2 temperature core BC
                                #=0; set core boundary temperature to ion temperature
                                #=1; set core boundary temperature to tgcore
                                #>1; set zero-temperature gradient over the core boundary
    if bbb.istgcore[igh2]==1:# Specify H2 core boundary temp
        bbb.tgcore[igh2] = tgcore # Core boundary temp

 
    ''' Wall BC '''
    bbb.tgwall[1] = 4.e-2      # Wall gas temperature when BC used, molecules



def common_mol(igh2,n0g,ngbackg,kelighg,kelighi,cngfx,cngfy,cfcvtg,cftgcond,ismolcrm):
        
    """---------------------------------------------------------------------------------------------
    BACKGROUND AND NORMALIZATION
    ---------------------------------------------------------------------------------------------"""
    bbb.n0g[igh2]=n0g     # Global hydrogen molecule normalization
    bbb.ngbackg[igh2] =ngbackg   # Set background H2 density

    """---------------------------------------------------------------------------------------------
    RATES
    ---------------------------------------------------------------------------------------------"""
    bbb.kelighg[igh2]=kelighg     # Elastic collision coefficient for gas i and hydrogen gas
    bbb.kelighi[igh2]=kelighi     # Elastic collision coeffisient for gas i and hydrogen ions
    bbb.ismolcrm=ismolcrm         # Switch for molecular particle and energy rates
    
    """---------------------------------------------------------------------------------------------
    SCALE FACTORS
    ---------------------------------------------------------------------------------------------"""
    bbb.cngfx[igh2]=cngfx       # Scale factor for flux from Grad(x)T_g in gas continuity equation
    bbb.cngfy[igh2]=cngfy       # Scale factor for flux from Grad(y)T_g in gas continuity equation
    bbb.cfcvtg=cfcvtg            # Convective thermal transport scaling coefficient for gas temperature
    bbb.cftgcond=cftgcond          # Conductive thermal transport scaling coefficient for gas temperature




def volsource(curr=1e15,h=0.5,w=0.5,r=0.01):
    ''' Molecular volumetric source '''
    # Volumetric sources
    bbb.ivolcurg[1]=curr     # Volume source [A] for EACH charge species
    bbb.zwng=r      # Width for volume source
    bbb.rwng=r      # Width for volume source
    bbb.z0ng=w*grd.zax
    bbb.r0ng=h*grd.radx


    '''
    from numpy import array
    # Define a cylindrical source in middle of domain
    c=array([   [0,     0,      .25,    .25,    0,      0   ],
                [0,     .25,    .5,     .5,     .25,    0   ],
                [.25,   .5,     1,      1,      .5,     .25 ],
                [.25,   .5,     1,      1,      .5,     .25 ],
                [0,     .25,    .5,     .5,     .25,    0   ],
                [0,     0,      .25,    .25,    0,      0   ]   ]   )
    c=c/sum(sum(c)) # Normalize to one
    # Get central coordinates

    
    h=int(h*com.nx)+1 
    w=int(w*com.ny)+1
    r=int(len(c)/2)
    # Set source strength
    bbb.psgov_use[h-r:h+r,w-r:w+r,1]=c*curr
    '''


def molbox(T,ismolcrm=1,ngbackg=1e10,istgon=0):

    # Add molecules density & temp
    bbb.ismolcrm=ismolcrm

    igh2=activate_mol()
    bbb.isngon[igh2] = 1
    bbb.istgcon[igh2] = -1
    if bbb.pyrestart_file[0].decode('UTF-8')!='read':
        bbb.allocate()  
    bbb.tgas = 0.02
    bbb.tgwall = 0.04
    bbb.istgon=0
    bbb.istgon[1]=istgon
    bbb.tgs[:,:,igh2]=T*bbb.ev
    bbb.ngbackg[igh2] = ngbackg		#floor level where background neut source on


