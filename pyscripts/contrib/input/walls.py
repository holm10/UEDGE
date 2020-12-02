# Wall model package for UEDGE, walls.py
# Separated from 160299 master input by holm10
# Changelog
# 200213 - file created
from uedge import bbb,com,grd,flx,api

def uniform_recycle_pump(recycH0=1e-10, recycH2=1, recycZ=1e-10, pumpH0=0.01, pumpH2=0.01, pumpZ=0,recycm=-0.9,nwsor=None):
    ''' Deck to set up recycling of gas species and spatially constant pumping of gaseous species
    H0H2_recycling_pump(recycH0,recycH2,pumpH0,pumpH2,pumpZ,recycm)

    Keyword parameters:
    recycH0[=1e-10]     Fraction of incident hydrogenic fluxes recycled as hydrogen atoms
    recycH2[=1]         Fraction of incident hydrogenic fluxes recycled as hydrogen molecules - only used if molecules are present
    recycZ[=1e-10]      Fraction of incident impurity fluxes recycled as impurity atoms
    pumpH0[=0.01]       Fraction of one-sided Maxwellian atom flux removed (albedo-like pumping)
    pumpH2[=0.01]       Fraction of one-sided Maxwellian molecule flux removed (albedo-like pumping) - only used if molecules are present
    pumpZ[=0.01]        Fraction of one-sided Maxwellian impurity flux removed (albedo-like pumping)
    recycm[=-0.9]       Momentum recycling of inertial atoms: default set by old comment

    This setup assumes the recycling and pumping coefficients to be the spatially constant and the 
    same for the targets and walls.
     recycH0+recycH2=1:  conserves particle balance
     recycH0+recycH2<1:  removes the the fraction 1-recycH0-recycH2 of the impiging fluxes
    '''
    """====================================================================================================
    WALL AND PLATE RECYCLING, PUFFING, PUMPING, ETC
    ===================================================================================================="""
    bbb.bcen=0          # Neutral energy transfer factor on plates
    bbb.bcenw=0         # Neutral energy transfer factor on walls

    """-----------------------------------------------------------------------------------------------------
    RECYCLING
    --------------------------------------------------------------------------------------------------------
    """
    bbb.recycp[0]=max(recycH0,1e-10)    # H0 plate recyc coeff
    bbb.recycw[0]=max(recycH0,1e-10)    # H0 wall recyc coeff

    # Only set molecular rates if molecules are present
    if bbb.ishymol:
        bbb.recycp[1]=max(recycH2,1e-10)    # H2 plate recyc coeff
        bbb.recycw[1]=max(recycH2,1e-10)    # H2 wall recyc coeff

    bbb.recycp[1+bbb.ishymol:]=max(recycZ,1e-10)     # C plate recyc coeff
    bbb.recycw[1+bbb.ishymol:]=max(recycZ,1e-10)     # C wall recyc coeff (>0!)





    """-----------------------------------------------------------------------------------------------------
    ALBEDOS - TARGET PUMPING
    -----------------------------------------------------------------------------------------------------"""
    # Sets albedo-like pumping
    # Affects the one-sided maxwellian contribution to the flux
    # albedolb[gas spexcies, X-point]

    bbb.albedolb[0,:]=1-pumpH0          # Left plate H0 albedo
    bbb.albedorb[0,:]=1-pumpH0          # Right plate H0 albedo

    # Only set molecular rates if molecules are present
    if bbb.ishymol:
        bbb.albedolb[1,:]=1-pumpH2          # Left plate H2 albedo
        bbb.albedorb[1,:]=1-pumpH2          # Right plate H2 albedo

    bbb.albedolb[1+bbb.ishymol:,:]=1-pumpZ           # Left plate Z albedo
    bbb.albedorb[1+bbb.ishymol:,:]=1-pumpZ           # Right plate Z albedo





    """-----------------------------------------------------------------------------------------------------
    MOMENTUM RECYCLING
    -----------------------------------------------------------------------------------------------------"""
    # Momentum recycling for inertial H:
        # (-9.9,inf):   ydot=-nurlxu*( recycm*up(H+)+up(H0) ) / Normalization
        # (-10.1,-9.9]: Zero-Neumann BC
        # (-inf,10.1]:  Neutral thermal flux used
    bbb.recycm=recycm     # Old comment: should be set to -0.9




    """-----------------------------------------------------------------------------------------------------
    GAS PUFFS AND PUMPS ON WALLS
    -----------------------------------------------------------------------------------------------------"""
    
    ''' ALL SPECIES INDICES BELOW USE BASIS INDEXING (FIRST INDEX=1) '''

    # Set the number of wall source regions 
    # Arrays below are only considered up to this BASIS index! [e.g. for nwsor=3 the 3 first indexes are considered)

    if nwsor is None:
        bbb.nwsor=1+bbb.ishymol     
        for ind in com.nzsp:
            if ind>0:
                bbb.nwsor+=1
    else:
        bbb.nwsor=nwsor
        
    pZ=1-pumpZ  # Impurity pumping rate helper


    # Puffs: indices correspond to source # - limited by nwsor
    bbb.issorlb=    [1, 1, 1, 1, 1, 1, 1, 1, 1, 1]  # Source position origin flag; =0: right plate, =1: left plate
    # Outer wall
    bbb.igspsoro=   [1, 2, 3, 4, 5, 6, 7, 8, 9, 10]  # Puffed species index
    bbb.igaso=  [0, 0, 0, 0, 0, 0, 0, 0, 0, 0]  # Gas currents from outer wall [Amp]
    bbb.xgaso=  [0, 0, 0, 0, 0, 0, 0, 0, 0, 0]  # Location of wall source: origin set by issorlb
    bbb.wgaso=  [1e3, 1e3, 1e3, 1e3, 1e3, 1e3, 1e3, 1e3, 1e3, 1e3]  # Total Cosine widths of source
    bbb.albdso= [1-pumpH0, (1-pumpH2)*bbb.ishymol+(1-bbb.ishymol)*pZ, pZ, pZ, pZ, pZ, pZ, pZ, pZ, pZ]  # Albedos at outer gas source locations
    bbb.matwso= [1, 1, 1, 1, 1, 1, 1, 1, 1, 1]  # Material wall BC at outer gas puff on if =1 - switch
    # Inner wall
    bbb.igspsori=   [1, 2, 3, 4, 5, 6, 7, 8, 9, 10]  # Puffed species index
    bbb.igasi=  [0, 0, 0, 0, 0, 0, 0, 0, 0, 0]  # Gas currents from inner wall [Amp]
    bbb.xgasi=  [0, 0, 0, 0, 0, 0, 0, 0, 0, 0]  # Location of wall source: origin set by issorlb
    bbb.wgasi=  [1e3, 1e3, 1e3, 1e3, 1e3, 1e3, 1e3, 1e3, 1e3, 1e3]  # Total Cosine widths of source
    bbb.albdsi= [1-pumpH0, (1-pumpH2)*bbb.ishymol+(1-bbb.ishymol)*pZ, pZ, pZ, pZ, pZ, pZ, pZ, pZ, pZ]  # Albedos at inner gas source locations
    bbb.matwsi= [1, 1, 1, 1, 1, 1, 1, 1, 1, 1]  # Material wall BC at inner gas puff on if =1 - switch




def target_wall_recycle_pump(   recycH0_wall=1e-10, recycH0_target=1e-10, recycH2_wall=1,recycH2_target=1, recycZ_target=1e-10, recycZ_wall=1e-10,
                                pumpH0_wall=0.01, pumpH0_target=0.01, pumpH2_wall=0.01, pumpH2_target=0.01, pumpZ_wall=0, pumpZ_target=0, recycm=-0.9, nwsor=None):

    ''' Deck to set up recycling of gas species: different rates at walls and targets
    H0H2_recycling_pump(recycH0,recycH2,pumpH0,pumpH2,pumpZ,recycm)

    Keyword parameters:
    recycH0[=1e-10]     Fraction of incident hydrogenic fluxes recycled as hydrogen atoms
    recycH2[=1]         Fraction of incident hydrogenic fluxes recycled as hydrogen molecules - only used if molecules are present
    recycZ[=1e-10]      Fraction of incident impurity fluxes recycled as impurity atoms
    pumpH0[=0.01]       Fraction of one-sided Maxwellian atom flux removed (albedo-like pumping)
    pumpH2[=0.01]       Fraction of one-sided Maxwellian molecule flux removed (albedo-like pumping) - only used if molecules are present
    pumpZ[=0.01]        Fraction of one-sided Maxwellian impurity flux removed (albedo-like pumping)
    recycm[=-0.9]       Momentum recycling of inertial atoms: default set by old comment

    This setup assumes the recycling and pumping coefficients to be the spatially constant and the 
    same for the targets and walls.
     recycH0+recycH2=1:  conserves particle balance
     recycH0+recycH2<1:  removes the the fraction 1-recycH0-recycH2 of the impiging fluxes
    '''
    """====================================================================================================
    WALL AND PLATE RECYCLING, PUFFING, PUMPING, ETC
    ===================================================================================================="""
    bbb.bcen=0          # Neutral energy transfer factor on plates
    bbb.bcenw=0         # Neutral energy transfer factor on walls

    """-----------------------------------------------------------------------------------------------------
    RECYCLING
    --------------------------------------------------------------------------------------------------------
    """
    bbb.recycp[0]=max(recycH0_target,1e-10)    # H0 plate recyc coeff
    bbb.recycw[0]=max(recycH0_wall,1e-10)    # H0 wall recyc coeff

    # Only set molecular rates if molecules are present
    if bbb.ishymol:
        bbb.recycp[1]=max(recycH2_target,1e-10)    # H2 plate recyc coeff
        bbb.recycw[1]=max(recycH2_wall,1e-10)    # H2 wall recyc coeff

    bbb.recycp[1+bbb.ishymol:]=max(recycZ_target,1e-10)     # C plate recyc coeff
    bbb.recycw[1+bbb.ishymol:]=max(recycZ_wall,1e-10)     # C wall recyc coeff (>0!)





    """-----------------------------------------------------------------------------------------------------
    ALBEDOS - TARGET PUMPING
    -----------------------------------------------------------------------------------------------------"""
    # Sets albedo-like pumping
    # Affects the one-sided maxwellian contribution to the flux
    # albedolb[gas spexcies, X-point]

    bbb.albedolb[0,:]=1-pumpH0_target          # Left plate H0 albedo
    bbb.albedorb[0,:]=1-pumpH0_target          # Right plate H0 albedo

    # Only set molecular rates if molecules are present
    if bbb.ishymol:
        bbb.albedolb[1,:]=1-pumpH2_target          # Left plate H2 albedo
        bbb.albedorb[1,:]=1-pumpH2_target          # Right plate H2 albedo

    bbb.albedolb[1+bbb.ishymol:,:]=1-pumpZ_target           # Left plate Z albedo
    bbb.albedorb[1+bbb.ishymol:,:]=1-pumpZ_target           # Right plate Z albedo





    """-----------------------------------------------------------------------------------------------------
    MOMENTUM RECYCLING
    -----------------------------------------------------------------------------------------------------"""
    # Momentum recycling for inertial H:
        # (-9.9,inf):   ydot=-nurlxu*( recycm*up(H+)+up(H0) ) / Normalization
        # (-10.1,-9.9]: Zero-Neumann BC
        # (-inf,10.1]:  Neutral thermal flux used
    bbb.recycm=recycm     # Old comment: should be set to -0.9




    """-----------------------------------------------------------------------------------------------------
    GAS PUFFS AND PUMPS ON WALLS
    -----------------------------------------------------------------------------------------------------"""
    
    ''' ALL SPECIES INDICES BELOW USE BASIS INDEXING (FIRST INDEX=1) '''

    # Set the number of wall source regions 
    # Arrays below are only considered up to this BASIS index! [e.g. for nwsor=3 the 3 first indexes are considered)
    if nwsor is None:
        bbb.nwsor=1+bbb.ishymol     
        for ind in com.nzsp:
            if ind>0:
                bbb.nwsor+=1
    else:
        bbb.nwsor=nwsor

    
    pZ=1-pumpZ_wall  # Impurity pumping rate helper


    # Puffs: indices correspond to source # - limited by nwsor
    bbb.issorlb=    [1, 1, 1, 1, 1, 1, 1, 1, 1, 1]  # Source position origin flag; =0: right plate, =1: left plate
    # Outer wall
    bbb.igspsoro=   [1, 2, 3, 4, 5, 6, 7, 8, 9, 10]  # Puffed species index
    bbb.igaso=  [0, 0, 0, 0, 0, 0, 0, 0, 0, 0]  # Gas currents from outer wall [Amp]
    bbb.xgaso=  [0, 0, 0, 0, 0, 0, 0, 0, 0, 0]  # Location of wall source: origin set by issorlb
    bbb.wgaso=  [1e3, 1e3, 1e3, 1e3, 1e3, 1e3, 1e3, 1e3, 1e3, 1e3]  # Total Cosine widths of source
    bbb.albdso= [1-pumpH0_wall, (1-pumpH2_wall)*bbb.ishymol+(1-bbb.ishymol)*pZ, pZ, pZ, pZ, pZ, pZ, pZ, pZ, pZ]  # Albedos at outer gas source locations
    bbb.matwso= [1, 1, 1, 1, 1, 1, 1, 1, 1, 1]  # Material wall BC at outer gas puff on if =1 - switch
    # Inner wall
    bbb.igspsori=   [1, 2, 3, 4, 5, 6, 7, 8, 9, 10]  # Puffed species index
    bbb.igasi=  [0, 0, 0, 0, 0, 0, 0, 0, 0, 0]  # Gas currents from inner wall [Amp]
    bbb.xgasi=  [0, 0, 0, 0, 0, 0, 0, 0, 0, 0]  # Location of wall source: origin set by issorlb
    bbb.wgasi=  [1e3, 1e3, 1e3, 1e3, 1e3, 1e3, 1e3, 1e3, 1e3, 1e3]  # Total Cosine widths of source
    bbb.albdsi= [1-pumpH0_wall, (1-pumpH2_wall)*bbb.ishymol+(1-bbb.ishymol)*pZ, pZ, pZ, pZ, pZ, pZ, pZ, pZ, pZ]  # Albedos at inner gas source locations
    bbb.matwsi= [1, 1, 1, 1, 1, 1, 1, 1, 1, 1]  # Material wall BC at inner gas puff on if =1 - switch







def impp_wallpuff(nzsor):
    ''' Set impurity wall puffing sources '''
    api.nzsor=nzsor                 # N.o. impurity wall sources
    if api.nzsor>0:
        api.iszsorlb=   [0, 0, 0, 0, 0, 0]      # Source location origin: 1=lb, 0=rb
        api.wimpi=  [1e3, 1e3, 1e3, 1e3, 1e3, 1e3]  # PF source zone width
        api.wimpo=  [1e3, 1e3, 1e3, 1e3, 1e3, 1e3]  # Outer wall source zone width
        api.impsori=    [0, 0, 0, 0, 0, 0]      # PF wall source strength [Amp]
        api.impsoro=    [0, 0, 0, 0, 0, 0]      #  wall source strength [Amp]


def box_walls(alba=1e-2,albm=0.99,recyca=1e-10,recycm=1):

    # Boundary conditions
    bbb.curcore[0] = 0.
    bbb.ncore[0]=3.0e19
    bbb.isnicore[0] = 1		#=1 gives bbb.ni=bbb.ncore on core bdry
    bbb.isngcore[0] = 0		#=0 gives api.zero flux
    bbb.isupcore[0:2] = 1	#=1 gives dup[,,1:2]/com.dy=0
    bbb.istgcore[1] = 2		#=2 gives dTg/com.dy=0 on core bdry
    bbb.iflcore = -1		#if=0, specify core temps if -1, dTe,i/com.dy=0 
    bbb.recycp[0] = recyca
    bbb.recycw[0] = recyca
    bbb.recycp[1] = recycm
    bbb.recycw[1] = recycm
    bbb.nwsor=2
    bbb.igspsoro[1]=2
    bbb.igspsori[1]=2
    bbb.albdsi[0]=alba
    bbb.albdso[0]=alba
    bbb.albdsi[1]=albm
    bbb.albdso[1]=albm
    bbb.albedorb[0] = alba  #pump 50% of atoms others reflected
    bbb.albedolb[0] = alba  #pump 50% of atoms others reflected
    bbb.albedorb[1] = albm  #pump 50% of atoms others reflected
    bbb.albedolb[1] = albm  #pump 50% of atoms others reflected
    bbb.istepfc=0
    bbb.istipfc=0	#bbb.priv. flux has api.zero temp. deriv.
    bbb.istewc=0
    bbb.istiwc=0	#wall has api.zero temp. deriv.
    bbb.istgpfc[ 1] = 2		#=3 gives dTg/com.dy=0 on PF wall
    bbb.istgwc[1] = 2		#=3 gives dTg/com.dy=0 on outer wall

def box_dirichlet(na=1e12,nm=1e12,psor=1e-2,recyca=1e-10,recycm=1e-10):

    # Boundary conditions
    bbb.curcore[0] = 0.
    bbb.ncore[0]=3.0e19
    bbb.isnicore[0] = 1		#=1 gives bbb.ni=bbb.ncore on core bdry
    bbb.isngcore[0] = 0		#=0 gives api.zero flux
    bbb.isupcore[0:2] = 1	#=1 gives dup[,,1:2]/com.dy=0
    bbb.istgcore[1] = 2		#=2 gives dTg/com.dy=0 on core bdry
    bbb.iflcore = -1		#if=0, specify core temps if -1, dTe,i/com.dy=0 
    bbb.recycp[0] = recyca
    bbb.recycw[0] = recyca
    bbb.recycp[1] = recycm
    bbb.recycw[1] = recycm
    bbb.nwsor=0
    bbb.istepfc=0
    bbb.istipfc=0	#bbb.priv. flux has api.zero temp. deriv.
    bbb.istewc=0
    bbb.istiwc=0	#wall has api.zero temp. deriv.
    bbb.istgpfc[ 1] = 2		#=3 gives dTg/com.dy=0 on PF wall
    bbb.istgwc[1] = 2		#=3 gives dTg/com.dy=0 on outer wall
    bbb.isfixlb=1
    bbb.isfixrb=1
    bbb.nib=na
    bbb.ywnio=1e4   
    bbb.ywnii=1e4   
    bbb.isnwconi=1
    bbb.isnwcono=1
    bbb.nwalli=na
    bbb.nwallo=na
    bbb.ngbackg[1]=nm



def tube_walls(alba=0.99,albm=0.99,albl=0.55):

    # Boundary conditions
    bbb.curcore[0] = 0.
    bbb.ncore[0]=3.0e19
    bbb.isnicore[0] = 1		#=1 gives bbb.ni=bbb.ncore on core bdry
    bbb.isngcore[0] = 0		#=0 gives api.zero flux
    bbb.isupcore[0:2] = 1	#=1 gives dup[,,1:2]/com.dy=0
    bbb.istgcore[1] = 2		#=2 gives dTg/com.dy=0 on core bdry
    bbb.iflcore = -1		#if=0, specify core temps if -1, dTe,i/com.dy=0 

    #bbb.isfixlb[0]=2 # Symmetry plane at LHS
    bbb.recycp[:2] = [1e-10,1]
    bbb.recycw[:2] = [1e-10,1]

    bbb.nwsor=2
    bbb.igspsoro[1]=2
    bbb.igspsori[1]=2
    bbb.albdsi[:2]=[alba,albm]
    bbb.albdso[:2]=[alba,albm]
    bbb.albedorb[0] = alba  #pump 50% of atoms others reflected
    bbb.albedolb[0] = albl  #pump 50% of atoms others reflected
    bbb.albedorb[1] = albm  #pump 50% of atoms others reflected
    bbb.albedolb[1] = albl  #pump 50% of atoms others reflected
    bbb.istepfc=0
    bbb.istipfc=0	#bbb.priv. flux has api.zero temp. deriv.
    bbb.istewc=0
    bbb.istiwc=0	#wall has api.zero temp. deriv.
    bbb.istgpfc[ 1] = 2		#=3 gives dTg/com.dy=0 on PF wall
    bbb.istgwc[1] = 2		#=3 gives dTg/com.dy=0 on outer wall



