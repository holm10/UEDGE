# Transport model package for UEDGE, transport.py
# Separated from 160299 master input by holm10
# Changelog
# 200213 - file created
from uedge import bbb,com,grd,flx


def const_diff( kye=2.6,kyi=0.75,
                difni=[2.6,0.5,2.6,2.6,2.6,2.6,2.6,2.6],travis=[2.6,1,1,1,1,1,1,1],
                flalfv=1.0, flalfe=0.21, flalfi=0.21, flalfgx=1.0, flalfgy=1.0, flalfgxy=1.0, 
                flalftgx=1.0, flalftgy=1.0, lgmax=0.5, lgtmax=0.1, lgvmax=0.1,
                kxe=1, kxi=1, parvis=1, fcdif=0,flalfvgx=1e20, flalfvgy=1e20, flalfvgxy=1e20
              ):
    ''' Sets the transport parameters and flux limits
        transport_fluxlim(*keys)
        
        Optional parameters
        kye (2.6)       -   Radial electron heat diffusivity
        kyi (0.75)      -   Radial ion heat diffusivity
        difni           -   Radial density diffusion coefficient for ion species ([2.6,0.5,2.6,2.6,2.6,2.6,2.6,2.6])
        travis          -   Perpendicular viscosity for ion species ([2.6,1,1,1,1,1,1,1])
        flalfv (1.0)    -   Parallel velocuty flux limit factor (linear)
        flalfe (0.21)   -   Parallel electron heat flux limit factor: 1/(1+abs(1/flalfe)**flgam)
        flalfi (0.21)   -   Parallel ion heat flux limit factor: 1/(1+abs(1/flalfi))
        flalfgx (1.0)   -   Poloidal gas diffusivity flux limit   
        flalfgy (1.0)   -   Radial gas diffusivity flux limit
        flalfgxy (1.0)  -   Nonorthogonal poloidal face gas flux limit
        flalftgx (1.0)  -   Poloidal gas temperature diffusivity flux limit
        flalftgy (1.0)  -   Radial gas temperature diffusivity flux limit
        lgmax (0.5)     -   Max gas scale length for calculating particle gaseous diffusivity coefficient
        lgtmax (0.1)    -   Max gas scale length for calculating thermal gas diffusivity coefficient
        lgvmax (0.1)    -   Max gas scale length for calculating viscous gas diffusivity coefficient
        kxe (1)         -   Poloidal electron heat condictivity multiplication factor
        kxi (1)         -   Poloidal ion heat conductivity multiplication factor
        parvis (1)      -   Paralell viscosity multiplication factor
        fcdif (0)       -   Scaling factor for constant nomalous diffusivities (fcdif*kye+kye_use, etc)
                                
                                

    '''

    """====================================================================================================
    TRANSPORT AND FLUX LIMITS
    ===================================================================================================="""


    # 8.1.1: Parallel transport
    #- - - - - - - - - -
    bbb.kxe=kxe      # Poloidal electron heat condictivity multiplication factor
    bbb.kxi=kxi      # Poloidal ion heat conductivity multiplication factor
    bbb.parvis[0]=parvis    # Paralell viscosity multiplication factor

    # 8.1.4: Anomalous radial transport
    #- - - - - - - - - - - - - - - - - 
    bbb.fcdif=fcdif       # Scaling factor for constant nomalous diffusivities (fcdif*kye+kye_use, etc)
                        # Not used for isbohmcalc==2

    bbb.isbohmcalc=2    # Bohm condition switch:
                            #=0: Uses the values set in dif_use, dif2_use,
                              # tra_use, kye_use, kyi_use, and difutm.
                            #=1: calc Bohm diff is facb...>0
                            #=2: Harmonic average of Bohm, difni, etc
                            #=3: D=difniv*(B0/B)**inbpdif, etc

    bbb.kye=kye     # Radial electron heat diffusivity
    bbb.kyi=kyi        # Radial ion heat diffusivity

    # H+
    bbb.difni[0]=difni[0]          # Radial H+ density diffusivity coefficient
    bbb.travis[0]=travis[0]         # H+ perpendicular viscosity

    # H0
    bbb.difni[1]=difni[1]         # Radial H+ density diffusivity coefficient
    bbb.travis[1]=travis[1]          # H+ perpendicular viscosity

    # C
    if bbb.isimpon>0:
        bbb.difni[com.nhsp+1:]=difni[com.nhsp+1:]    # Radial H+ density diffusivity cefficient
        bbb.travis[0]=travis[com.nhsp+1:]       # H+ perpendicular viscosity

    # Set fluxlimits
    fluxlim(    flalfv, flalfe, flalfi, flalfgx, flalfgy, flalfgxy, 
                flalftgx, flalftgy, flalfvgx, flalfvgy, flalfvgxy,
                lgmax, lgtmax, lgvmax)






def radial_diff(    kye=2.6,kyi=0.75,
                    difni=[2.6,0.5,2.6,2.6,2.6,2.6,2.6,2.6],travis=[2.6,1,1,1,1,1,1,1],
                    flalfv=1.0, flalfe=0.21, flalfi=0.21, flalfgx=1.0, flalfgy=1.0, flalfgxy=1.0, 
                    flalftgx=1.0, flalftgy=1.0, lgmax=0.5, lgtmax=0.1, lgvmax=0.1, 
                    kxe=1, kxi=1, parvis=1, fcdif=0,flalfvgx=1e20, flalfvgy=1e20, flalfvgxy=1e20

                      ):
    ''' Sets the transport parameters and flux limits
        transport_fluxlim(*keys)
        
        Optional parameters
        kye (2.6)       -   Radial electron heat diffusivity
                                -List of radial values optional 
        kyi (0.75)      -   Radial ion heat diffusivity
                                -List of radial values optional 
        difni           -   Radial density diffusion coefficient for ion species ([2.6,0.5,2.6,2.6,2.6,2.6,2.6,2.6])
                                -List of radial values otional (2D array)
        travis          -   Perpendicular viscosity for ion species ([2.6,1,1,1,1,1,1,1])
                                -List of radial values optional (2D array)
        flalfv (1.0)    -   Parallel velocuty flux limit factor (linear)
        flalfe (0.21)   -   Parallel electron heat flux limit factor: 1/(1+abs(1/flalfe)**flgam)
        flalfi (0.21)   -   Parallel ion heat flux limit factor: 1/(1+abs(1/flalfi))
        flalfgx (1.0)   -   Poloidal gas diffusivity flux limit   
        flalfgy (1.0)   -   Radial gas diffusivity flux limit
        flalfgxy (1.0)  -   Nonorthogonal poloidal face gas flux limit
        flalftgx (1.0)  -   Poloidal gas temperature diffusivity flux limit
        flalftgy (1.0)  -   Radial gas temperature diffusivity flux limit
        lgmax (0.5)     -   Max gas scale length for calculating particle gaseous diffusivity coefficient
        lgtmax (0.1)    -   Max gas scale length for calculating thermal gas diffusivity coefficient
        lgvmax (0.1)    -   Max gas scale length for calculating viscous gas diffusivity coefficient
        kxe (1)         -   Poloidal electron heat condictivity multiplication factor
        kxi (1)         -   Poloidal ion heat conductivity multiplication factor
        parvis (1)      -   Paralell viscosity multiplication factor
                                
                                

    '''



    """---------------------------------------------------------------------------------------------------- 
    TRANSPORT PARAMETERS
    -----------------------------------------------------------------------------------------------------"""


    # 8.1.1: Parallel transport
    #- - - - - - - - - -
    bbb.kxe=kxe      # Poloidal electron heat condictivity multiplication factor
    bbb.kxi=kxi      # Poloidal ion heat conductivity multiplication factor
    bbb.parvis[0]=parvis    # Paralell viscosity multiplication factor

    # 8.1.4: Anomalous radial transport
    #- - - - - - - - - - - - - - - - - 
    bbb.fcdif=fcdif       # Scaling factor for constant nomalous diffusivities (fcdif*kye+kye_use, etc)
                        # Not used for isbohmcalc==2

    bbb.isbohmcalc=3    # Bohm condition switch:
                            #=0: Uses the values set in dif_use, dif2_use,
                              # tra_use, kye_use, kyi_use, and difutm.
                            #=1: calc Bohm diff is facb...>0
                            #=2: Harmonic average of Bohm, difni, etc
                            #=3: D=difniv*(B0/B)**inbpdif, etc

    bbb.inbtdif=0       # Exponential scaling of (B0/B)
    bbb.inbpdif=0       # Exponential scaling of (B0/B)

    bbb.kyev=kye        # Radial electron heat diffusivity
    bbb.kyiv=kyi       # Radial ion heat diffusivity

    # H+
    bbb.difniv[:,0]=difni[0]           # Radial H+ diffusivity
    bbb.travisv[:,0]=travis[0]

    # H0
    bbb.difniv[:,1]=difni[1]            # Radial H0 diffusivity
    bbb.travisv[:,1]=travis[1]

    # C
    if bbb.isimpon>0:
        for ind in range(com.nhsp,bbb.difniv.shape[1]):
            bbb.difniv[:,ind]=difni[ind]   # Radial C diffusivity
            bbb.travisv[:,ind]=travis[ind]    # Perpendicular viscosity

    # Set fluxlimits
    fluxlim(    flalfv, flalfe, flalfi, flalfgx, flalfgy, flalfgxy, 
                flalftgx, flalftgy, flalfvgx, flalfvgy, flalfvgxy,
                lgmax, lgtmax, lgvmax)




def user_diff(    kye=2.6,kyi=0.75,
                    difni=[2.6,0.5,2.6,2.6,2.6,2.6,2.6,2.6],travis=[2.6,1,1,1,1,1,1,1],
                    flalfv=1.0, flalfe=0.21, flalfi=0.21, flalfgx=1.0, flalfgy=1.0, flalfgxy=1.0, 
                    flalftgx=1.0, flalftgy=1.0, lgmax=0.5, lgtmax=0.1, lgvmax=0.1, 
                    kxe=1, kxi=1, parvis=1, fcdif=0,flalfvgx=1e20, flalfvgy=1e20, flalfvgxy=1e20

                      ):
    ''' Sets the transport parameters and flux limits
        transport_fluxlim(*keys)
        
        Optional parameters
        kye (2.6)       -   Radial electron heat diffusivity
                                -List of radial values optional 
        kyi (0.75)      -   Radial ion heat diffusivity
                                -List of radial values optional 
        difni           -   Radial density diffusion coefficient for ion species ([2.6,0.5,2.6,2.6,2.6,2.6,2.6,2.6])
                                -List of radial values otional (2D array)
        travis          -   Perpendicular viscosity for ion species ([2.6,1,1,1,1,1,1,1])
                                -List of radial values optional (2D array)
        flalfv (1.0)    -   Parallel velocuty flux limit factor (linear)
        flalfe (0.21)   -   Parallel electron heat flux limit factor: 1/(1+abs(1/flalfe)**flgam)
        flalfi (0.21)   -   Parallel ion heat flux limit factor: 1/(1+abs(1/flalfi))
        flalfgx (1.0)   -   Poloidal gas diffusivity flux limit   
        flalfgy (1.0)   -   Radial gas diffusivity flux limit
        flalfgxy (1.0)  -   Nonorthogonal poloidal face gas flux limit
        flalftgx (1.0)  -   Poloidal gas temperature diffusivity flux limit
        flalftgy (1.0)  -   Radial gas temperature diffusivity flux limit
        lgmax (0.5)     -   Max gas scale length for calculating particle gaseous diffusivity coefficient
        lgtmax (0.1)    -   Max gas scale length for calculating thermal gas diffusivity coefficient
        lgvmax (0.1)    -   Max gas scale length for calculating viscous gas diffusivity coefficient
        kxe (1)         -   Poloidal electron heat condictivity multiplication factor
        kxi (1)         -   Poloidal ion heat conductivity multiplication factor
        parvis (1)      -   Paralell viscosity multiplication factor
                                
                                

    '''



    """---------------------------------------------------------------------------------------------------- 
    TRANSPORT PARAMETERS
    -----------------------------------------------------------------------------------------------------"""


    # 8.1.1: Parallel transport
    #- - - - - - - - - -
    bbb.kxe=kxe      # Poloidal electron heat condictivity multiplication factor
    bbb.kxi=kxi      # Poloidal ion heat conductivity multiplication factor
    bbb.parvis[0]=parvis    # Paralell viscosity multiplication factor

    # 8.1.4: Anomalous radial transport
    #- - - - - - - - - - - - - - - - - 
    bbb.fcdif=fcdif       # Scaling factor for constant nomalous diffusivities (fcdif*kye+kye_use, etc)
                        # Not used for isbohmcalc==2

    bbb.isbohmcalc=1    # Bohm condition switch:
                            #=0: Uses the values set in dif_use, dif2_use,
                              # tra_use, kye_use, kyi_use, and difutm.
                            #=1: calc Bohm diff is facb...>0
                            #=2: Harmonic average of Bohm, difni, etc
                            #=3: D=difniv*(B0/B)**inbpdif, etc


    bbb.kye=kye        # Radial electron heat diffusivity
    bbb.kyi=kyi       # Radial ion heat diffusivity

    # H+
    bbb.difni[0]=difni[0]           # Radial H+ diffusivity
    bbb.travis[0]=travis[0]

    # H0
    bbb.difni[1]=difni[1]            # Radial H0 diffusivity
    bbb.travis[1]=travis[1]

    # C
    if bbb.isimpon>0:
        for ind in range(com.nhsp,bbb.difniv.shape[1]):
            bbb.difni[ind]=difni[ind]   # Radial C diffusivity
            bbb.travis[ind]=travis[ind]    # Perpendicular viscosity

    # Set fluxlimits
    fluxlim(    flalfv, flalfe, flalfi, flalfgx, flalfgy, flalfgxy, 
                flalftgx, flalftgy, flalfvgx, flalfvgy, flalfvgxy,
                lgmax, lgtmax, lgvmax)







def fluxlim(flalfv, flalfe, flalfi, flalfgx, flalfgy, flalfgxy, 
            flalftgx, flalftgy, flalfvgx, flalfvgy, flalfvgxy,
            lgmax, lgtmax, lgvmax
           ):                      
    """---------------------------------------------------------------------------------------------------- 
    FLUX LIMIT FACTORS
    ----------------------------------------------------------------------------------------------------""" 
    bbb.flalfv=flalfv      # Parallel velocuty flux limit factor (linear)
    bbb.flalfe=flalfe     # Parallel electron heat flux limit factor: 1/(1+abs(1/flalfe)**flgam)
    bbb.flalfi=flalfi     # Parallel ion heat flux limit factor: 1/(1+abs(1/flalfi))
    bbb.flalfgx=flalfgx     # Poloidal gas diffusivity flux limit   
    bbb.flalfgy=flalfgy     # Radial gas diffusivity flux limit
    bbb.flalfgxy=flalfgxy    # Nonorthogonal poloidal face gas flux limit
    bbb.flalftgx=flalftgx    # Poloidal gas temperature diffusivity flux limit
    bbb.flalftgy=flalftgy    # Radial gas temperature diffusivity flux limit
    bbb.flalfvgx=flalfvgx   # Poloidal gas parallel viscosity flux limit 
    bbb.flalfvgy=flalfvgy   # Radial gas parallel viscosity flux limit
    bbb.flalfvgxy=flalfvgxy # Flux limit for neutral fmixy (co-perp to radial and parallel)
    bbb.lgmax=lgmax       # Max gas scale length for calculating particle gaseous diffusivity coefficient
    bbb.lgtmax=lgtmax      # Max gas scale length for calculating thermal gas diffusivity coefficient
    bbb.lgvmax=lgvmax      # Max gas scale length for calculating viscous gas diffusivity coefficient




def set_constant_div(D_perp=1,Chi_e=2):
    ''' Freezes the radial transport parameters and sets constant diffusivities and e- heat diff below X-point
        set_constant_div(*keys)
        
        Optional parameters
        D_perp (1)  -   Constant divertor diffusivities for all species
        Chi_e (2)   -   Constant electron heat diffusivities

    '''
    """====================================================================================================
    SET CONSTANT DIVERTOR DIFFUSIVITIES
    ===================================================================================================="""


    # Execute exmain without evaluating the convergence to generate and populate
    # the arrays used by UEDGE. 
    bbb.issfon=0
    bbb.ftol=1e20
    bbb.exmain()

    # Use *_use arrays for radial diffusivities
    bbb.isbohmcalc=0

    # Set diffusivities of the *_use arrays in the divertor SOL 
    # (common and private flux region) to 1.0 m2/s (D_perp)
    # and 2.0 m2/s (kye). Set the ion particle diffusivities 
    # in the divertor to 1 m2/s 
    # H+
    bbb.dif_use[:com.ixpt1[0]+1,:,0]=D_perp
    bbb.dif_use[com.ixpt2[0]+1:,:,0]=D_perp
    # Carbon
    bbb.dif_use[:com.ixpt1[0]+1,:,com.nhsp:]=D_perp
    bbb.dif_use[com.ixpt2[0]+1:,:,com.nhsp:]=D_perp

    # Set electron heat diffusivities in divertor to 2m2/s
    bbb.kye_use[:com.ixpt1[0]+1,:]=Chi_e
    bbb.kye_use[com.ixpt2[0]+1:,]=Chi_e

    # Switches to solve equations to convergence
    bbb.issfon=1
    bbb.ftol=1e-8
   







def dif_groth(    jtanh = 8.0,    jw = 1.5,   jexp = 18.0,    jg = 4.0,
            d1 = 0.15,      d2 = 0.75,  d3=0.25                     ):
    ''' Radial particle diffusivity profile
        dif(jtanh=8.0,jw=1.5,jexp=18.0,jg=4.0,d1=0.15,d2=0.75,d3=0.25)

    '''
    from uedge import com,tanh,exp
    from numpy import zeros

    ret=zeros((com.ny+2,))

    for i in range(com.ny+2):
        ret[i]=(d2 + (d1-d2) * (tanh( (i-jtanh)/jw ) + 1 ) / 2.0 + d3 * exp( ( i - jexp ) / jg ))
    
    return ret




def kye_groth(    jtanh = 12.0,   jw = 1.0,   chi_core = 1.6, chi_edge = 2.0  ):
    ''' Radial electron heat diffusivity profile
        kye(jtanh=12.0,jw=1.0,chi_core=1.6,chi_edge=2.0)
    '''
    from uedge import com,tanh,exp
    from numpy import zeros

    ret=zeros(com.ny+2)
    for i in range(com.ny + 2):
        ret[i]=( chi_core + 0.5 * ( tanh( ( i - jtanh ) / jw ) + 1 ) * (chi_edge-chi_core))

    return ret


def box_transport():

    bbb.flalfe=0.21
    bbb.flalfi=0.21
    bbb.flalfv=1.0
    bbb.flalfgx=1.
    bbb.flalfvgx=1.
    bbb.flalftgx=1.
    bbb.flalfvgy=1.
    bbb.flalftgy=1.
    bbb.isplflxlgx=0 # Turn off flux lim at plates -  causes issues when 
                     # atom temperatures decrease below a threshold temperature
    
        
