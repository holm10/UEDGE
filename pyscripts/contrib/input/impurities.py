# Impurity package for UEDGE, impurities.py
# Separated from 160299 master input by holm10
# Changelog
# 200213 - file created
from uedge import bbb,com,grd,flx,api



def carbon_forcebalance(    apipath='.', isnicore=3, ncore=4e15, curcore=0, isupcore=1, isnwcono=3, nwomin=1e7, 
                            isnwconi=3, nwimin=1e7, tgwall=4e-2, isbohmms=0, n0=1e17, n0g=1e18, nzbackg=1e10, 
                            inzb=2,ngbackg=1e10, ismctab=2, mcfilename='C_rates.adas', isrtndep=1, rcxighg=0, 
                            iscxfit=2, kelighi=5e-16, kelighg=5e-16, isch_sput=7, fchemygwi=1, fchemygwo=1, 
                            fchemylb=1, fchemyrb=1, isph_sput=3, fphysylb=1, fphysyrb=1, isi_sputw=2, isi_sputpf=2, 
                            t_wall=300, t_plat=300,isngcore=0,ngcore=1e15
                       ):


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

    api.apidir=apipath # Impurity rates
    # IMPURITY SETUP
    #- - - - - - - -
    # Helper indices
    Ccs=6                   # Carbon charge states
    

    zind=com.ngsp-bbb.ishymol-1   # Carbon impurity species index
    ingc=com.ngsp
    iicl=com.nhsp           # Determine the lower ion index for carbon
    for z in com.nzsp:
        iicl+=z
    iicu=iicl+Ccs


    # Allocate space
    com.ngsp+=1             # 1 gaseous C species 
    com.nzsp[zind]=Ccs      # Number of impurity species for gas species species+1
                            # Determines nisp (nisp=nhsp+sum(nzsp)) and allocates arrays

    # Turn impurity gas on
    #- - - - - - - - - - - 
    #bbb.istgcon[com.ngsp]=0 # Impurity gas temperature: tg=(1-istgcon)*rtg2ti*ti+istgcon*tgas*ev

    # Impurity species parameters
    #- - - - - - - - - - - - - - -
    bbb.ziin[iicl:iicu]=range(1,Ccs+1)     # Impurity charge states
    bbb.minu[iicl:iicu]=12      # Atomic mass unit species mass
    bbb.znuclin[iicl:iicu]=6  # Nuclear charge of impurities

    # Core BC
    #-----------------------------------------------------------------------------------------------
    bbb.isnicore[iicl:iicu]=isnicore   # Core impurity ion BC model
                                    #=0: flux = curcore/sy locally in ix
                                    #=1: uniform, fixed density, ncore
                                    #=2: flux & ni over range
                                    #=3: icur=curcore-recycc*fngy, const ni
                                    #=4: impur. source terms (impur only)
                                    #=5: d(ni)/dy=-ni/lynicore at midp & ni constant poloidall
    if 1 in bbb.isnicore[iicl:iicu]:    # Constant imputiry dens on core vboundary
        bbb.ncore[iicl:iicu]=ncore      # Core boundary density
    if 3 in bbb.isnicore[iicl:iicu]:    # Core BC set by fluxes
        bbb.curcore[iicl:iicu]=curcore        # Constant flux contribution over core BC

    bbb.isupcore[iicl:iicu]=isupcore   # Velocity BC:
                                    #=0; Dirichlet BC: up=upcore
                                    #=1; Neumann BC: d(up)/dy=0
                                    #=2 sets d^2(up)/dy^2 = 0
                                    #=3 sets poloidal velocity (uu) = 0
                                    #=4 sets tor. ang mom flux = lzflux & n*uz/R=const
                                    #=5 sets ave tor vel = utorave & n*uz/R=const

    bbb.isngcore[ingc]=isngcore # Neutral gas density BC at core
				                # 0 - set loc flux= -(1-albedoc)*ng*vtg/4
				                #=1, set uniform, fixed density, ngcore
				                #=2, not available
				                #=3, extrapolation, but limited
				                #=anything else, set zero deriv which was
				                #prev default inert hy
				                # anything else same as =0

    if bbb.isngcore[ingc]==1:
        bbb.ngcore[ingc]=1e9    # Set Uniform density if BC requested


    # Carbon wall BC
    #-----------------------------------------------------------------------------------------------
    # Outer wall - ION
    #- - - - - - 


    bbb.isnwcono[iicl:iicu]=isnwcono# Outer wall BC:
                                #=0, old case; if ifluxni=0, dn/dy=0; if ifluxni=1, fniy=0 (default)
                                #=1, fixed density to nwallo(ix) array
                                #=2, extrapolation B.C.
                                #=3, approx gradient scale length, limited by nwomin
    if 3 in bbb.isnwcono[iicl:iicu]:    # Grad scale length BC: same as for plasma
        bbb.nwomin[iicl:iicu]=nwomin      # Minimum outer wall density
    
    # PFR wall - ION
    #- - - - -
    bbb.isnwconi[iicl:iicu]=isnwconi# PFR wall BC: as above
    if 3 in bbb.isnwconi[iicl:iicu]:    # Grad scale length BC: same as for plasma
        bbb.nwimin[iicl:iicu]=nwimin      # Minimum PFR wall density
    
    # Gas wall BC
    #- - - - - - 
    bbb.tgwall[ingc] = tgwall      # Wall gas temperature when BC used, carbon

    # Impurity plate BC:s
    #-----------------------------------------------------------------------------------------------
    bbb.isbohmms=isbohmms      #0=single-species Bohm condition (H+)
                #1=multi-species Bohm condition (all ions)

    # Normalization and background
    #-----------------------------------------------------------------------------------------------
    bbb.n0[iicl:iicu]=n0 # Global impurity ion density normalization
    bbb.n0g[ingc]=n0g     # Global impurity gas density normalization
    bbb.nzbackg=nzbackg       # Background impurity ion density
    bbb.inzb=inzb          # Impurity floor scaling (nzbackg/ni)^inzb
    bbb.ngbackg[ingc]=ngbackg # Impurity gas background density 




    # Impurity rates
    #-----------------------------------------------------------------------------------------------
    # Setup rate model
    #- - - - - - - - -
    bbb.ismctab=ismctab       # Define data to be used for multi-charge-state rates
                            #=1: tables originally generated by R. Campbell for D. Knoll,
                              # data file name is specified by inelmc=....
                              # corresponding rate evaluation routines are imprates and radimpmc.
                            #=2: tables generated by code from B. Braams,
                              # data file name is specified by mcfilename=...,
                              # corresponding rate evaluation routines are mcrates and radmc
    if bbb.ismctab==2:  # Braams tables
        com.mcfilename[0]=mcfilename   # Rate data to be used
        com.isrtndep=isrtndep          # Are table lookup parameters density dependent
                            # Check compability with hydrogen if istabon=16
    
    # CX
    #- - -
    bbb.rcxighg[ingc]=rcxighg       # Turn off imp0 + H+ --> imp+ + H+
    
    com.iscxfit=iscxfit           # C-ion + H0 CX model:
                                #=0: analytic forms in Braams' rate package
                                #=1: polynomial fit to C.F. Maggi curves (1997)
                                #=2: same as =1, except Z=1 has lower rate from Pigarov

    # Scattering
    #- - - -  - -
    bbb.kelighi[ingc] = kelighi  # Elastic collision coefficient with H+
    bbb.kelighg[ingc] = kelighg  # Elastic collision coefficient with H0


    # Impurity sputtering
    #-----------------------------------------------------------------------------------------------
    

    # Chemical
    #- - - - - 
    bbb.isch_sput[ingc]=isch_sput   # Chemical sputtering model
                                #=0: Old
                                #=5: Roth, G-R
                                #=6: Haasz 97
                                #=7: Haasz 97 + Davis at low E

    bbb.fchemygwi=  fchemygwi   # Factor multiplying chemical sputtering gas yield; PF wall
    bbb.fchemygwo=  fchemygwo   # Factor multiplying chemical sputtering gas yield; Outer wall
    bbb.fchemylb=   fchemylb   # Factor multiplying chemical sputtering gas yield; Left plate
    bbb.fchemyrb=   fchemyrb   # Factor multiplying chemical sputtering gas yield; Right plate




    # Physical
    #- - - - - - 
    bbb.isph_sput[ingc]=3   # Physical sputtering model
                                #=0: old fixed case
                                #=1: DIVIMP/JET physical sputtering fits
                                #=2: adds H+ chemical sputtering
                                #=3: adds H0 carbon sputtering
    bbb.crmb=bbb.minu[0]  # Mass of incident sputtering particles
    bbb.cizb=bbb.ziin[0]  # Max plasma charge state

    bbb.fphysylb=   fphysylb   # Factor multiplying physical sputtering gas yield; Left plate  
    bbb.fphysyrb=   fphysyrb   # Factor multiplying physical sputtering gas yield; Right plate

    # Wall sputtering
    #- - - - - - 
    bbb.isi_sputw[ingc]=isi_sputw   # Outer wall sputtering model
                                #=0: no ion sputtering
                                #=1: adds physical ion sputtering
                                #=2: adds chemical ion sputtering
    bbb.isi_sputpf[2]=isi_sputpf # PF wall sputtering: as above
    
    bbb.t_wall=t_wall      # Side wall temperatures
    bbb.t_plat=t_plat      # Plate temperatures


    """---------------------------------------------------------------------------------------------------- 
    ALLOCATE CARBON ARRAYS
    ----------------------------------------------------------------------------------------------------"""
    if bbb.pyrestart_file[0].decode('UTF-8')!='read':
        bbb.allocate()



def wall_source(nzsor=1):
    ''' Defines an impurity wall source '''
    # TODO: Needs fixing
    api.nzsor=nzsor         # number of wall source zones for impurity
    bbb.allocate()          # Allocate wall source arrays
    api.iszsorlb=0          # measure from left boundary
    api.wimpi=1000.         # width of impurity source zone in PF
    api.impsori[0]=0.           # PF wall source strength [A]
    api.wimpo=1000.         # width of impurity source zone on outer wall
    api.impsoro[0]=0.           # outer wall source strength [A]
    bbb.ivolcur=0.0         # volume source [A] for EACH charge species
    bbb.zwni=1000.          # width for volume source
    bbb.rwni=1000.          # width for volume source



