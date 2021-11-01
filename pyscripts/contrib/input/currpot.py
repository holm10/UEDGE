# Current and potential package for UEDGE, currpot.py
# Separated from 160299 master input by holm10
# Changelog
# 200213 - file created
# 200216 - Created common blocks for all b-directions
from uedge import bbb


def slabB(cfniybbo=0,cfniydbo=0,cfeeybbo=0,cfeeydbo=0,iphibcc=3,cfeexdbo=0):
    bbb.newbcl=0        # Linear scaling factor for left plate sheath potential:
                            #=0: use bcei and bcee
                            #=1: new model
    bbb.newbcr=0        # Linear scaling factor for right plate sheath potential:
                            #=0: use bcei and bcee
                            #=1: new model
    bbb.b0=1

    drift_const(    cfjhf=0,cfjve=0,jhswitch=0,isfdiax=0,cfqydbo=0,cftdd=0,cfyef=0,cftef=0,cf2ef=0,
                    cfybf=0, cf2bf=0, cfqybf=0, cfq2bf=0, cfniydbo=0, cfeeydbo=0, cfeixdbo=0, cfeexdbo=0,
                    cfniybbo=cfniybbo,cfeeybbo=cfeeybbo)
    potential(iphibcc=iphibcc)


def fwdB(cfniybbo=0,cfniydbo=1,cfeeybbo=0,cfeeydbo=1,iphibcc=3,cfeexdbo=1):
    bbb.newbcl=1        # Linear scaling factor for left plate sheath potential:
                            #=0: use bcei and bcee
                            #=1: new model
    bbb.newbcr=1        # Linear scaling factor for right plate sheath potential:
                            #=0: use bcei and bcee
                            #=1: new model
    bbb.b0=1

    drift_const(cfniybbo=cfniybbo,cfniydbo=cfniydbo,cfeeybbo=cfeeybbo,cfeeydbo=cfeeydbo,cfeexdbo=cfeexdbo)
    potential(iphibcc=iphibcc)

def revB(cfniybbo=0,cfniydbo=1,cfeeybbo=0,cfeeydbo=1,iphibcc=3,cfeexdbo=1):
    bbb.newbcl=1        # Linear scaling factor for left plate sheath potential:
                            #=0: use bcei and bcee
                            #=1: new model
    bbb.newbcr=1        # Linear scaling factor for right plate sheath potential:
                            #=0: use bcei and bcee
                            #=1: new model
    bbb.b0=-1


    drift_const(cfniybbo=cfniybbo,cfniydbo=cfniydbo,cfeeybbo=cfeeybbo,cfeeydbo=cfeeydbo,cfeexdbo=cfeexdbo)
    potential(iphibcc=iphibcc)

def noB(cfniybbo=0,cfniydbo=1,cfeeybbo=0,cfeeydbo=1,iphibcc=3, cfeexdbo=1):
    bbb.newbcl=1        # Linear scaling factor for left plate sheath potential:
                            #=0: use bcei and bcee
                            #=1: new model
    bbb.newbcr=1        # Linear scaling factor for right plate sheath potential:
                            #=0: use bcei and bcee
                            #=1: new model
    bbb.b0=100


    drift_const(cfniybbo=cfniybbo,cfniydbo=cfniydbo,cfeeybbo=cfeeybbo,cfeeydbo=cfeeydbo,cfeexdbo=cfeexdbo)
    potential(iphibcc=iphibcc)


def drift_const( rsigpl=1e-8,cfjhf=1,cfjve=1,jhswitch=1,isfdiax=1,cfqydbo=1,
                        cfydd=0,cf2dd=0,cftdd=1,cfyef=1,cftef=1,cf2ef=1,cfjpy=0,cfjp2=0,
                        cfybf=1, cf2bf=1, cfqybf=1, cfq2bf=1, cfqybbo=0, cfniybbo=0,cfeeybbo=0,
                        cfniydbo=1, cfeeydbo=1, cfeixdbo=1, cfeexdbo=1, cfqym=1
                      ):
    """-----------------------------------------------------------------------------------------------------
    CURRENT AND POTENTIAL PARAMETERS FOR UEDGE SIMULATIONS OF DIII-D SHOT 160299, FWD-B
    -----------------------------------------------------------------------------------------------------"""


    bbb.rsigpl=rsigpl    # Anomalous radial electrical conductivity
    bbb.cfjhf=cfjhf        # Coefficient for convective current (fqp) heat flow
    bbb.cfjve=cfjve        # Coefficient for current contribution to electron velocity: vex=vix-cfjve*fqx
    bbb.cfjpy = cfjpy       # Coefficient for diamagnetic drift in y-direction
    bbb.cfjp2 = cfjp2       # Coefficient for diamagnetic drift in direction perpendicular to radial and parallel directions
    bbb.jhswitch=jhswitch      #Joule Heating term coefficient

    # 5.3.3: DRIFTS
    #- - - - - - - - 
    # Diamagnetic
    bbb.isfdiax=isfdiax       # Diamagnetic drift sheat contribution coefficient
    bbb.cfqydbo=cfqydbo       # Coefficient for diagmagetic current on core boundary (=1 forces j_r=0)
    bbb.cfydd=cfydd       # Divergence-free diamagnetic drift coefficient in radial direction
    bbb.cf2dd=cf2dd       # Diamagnetic drift coefficient in direction perpendicular to parallel and radial direction
    bbb.cftdd=cftdd       # Coefficient for diamagnetic contribution on toroiudal drift

    # ExB
    bbb.cfyef=cfyef       # ExB drift coefficient in y-direction
    bbb.cftef=cftef       # ExB drift coefficient in toroidal direction
    bbb.cf2ef=cf2ef       # ExB drift in direction perpendicular to radial and parallel direction

    # Grad(B)
    bbb.cfybf=cfybf       # Grad(B) radial drift coefficient 
    bbb.cf2bf=cf2bf       # Grad(B) drift in direction perpendicular no parallel and radial direcion
    bbb.cfqybf=cfqybf      # Grad(B) radial current coefficient 
    bbb.cfq2bf=cfq2bf      # Grad(B) current in direction perpendicular no parallel and radial direcion
    bbb.cfqybbo=cfqybbo       # Grad(B) current coefficient on core boundary
    bbb.cfniybbo=cfniybbo     # Grad(B) coefficient of velocity contribution to ion/neutral flux on core boundary
    bbb.cfeeybbo=cfeeybbo     # Grad(B) coefficient of velocity contribution to electron flux on core boundary

    # Grad(P_i x B)
    bbb.cfniydbo=cfniydbo     # Grad(P_i x B) coefficient of velocity contribution to ion/neutral flux on core boundary
    bbb.cfeeydbo=cfeeydbo     # Grad(P_i x B) coefficient of velocity contribution to electron flux on core boundary
    # BxGrad(T)
    bbb.cfeixdbo=cfeixdbo     # Coeffcient for BxGrad[T] drift in plate BC
    bbb.cfeexdbo=cfeexdbo     # Coefficient for diamagnetic drift in plate BC
    # Inertial correction
    bbb.cfqym=cfqym       # Coefficient for spatial inertial radial current in Y-dir (fqy)



def potential(isnewpot=1,rnewpot=1,iphibcc=3,isutcore=2,iphibcwi=0,iphibcwo=0):
    # 5.3.4: POTENTIAL
    #- - - - - - - - -
    bbb.isnewpot=isnewpot      # Potential model switch:
                            #=1: new potential model; radial curr dens from toroidal momentum balance
                            #=-2: Constant Phi on core boundary with total core current = icoreelec
    bbb.rnewpot=rnewpot     # Linear scaling factor for potential model: 1=new, 0=old

    if bbb.isnewpot==1: # New potential equation
        bbb.iphibcc=iphibcc       # Core potential boundary equation:
                                #=1: d^2(ey)/dy^2=0
                                #=2: te=constant & ey(ixmp,0)=eycore
                                #=3:  phi=constant & ey(ixmp,0)=eycore
                                #else:  dphi(ix,1)=dphi_iy1; isutcore controls ix=ixmp 
        if bbb.iphibcc not in [1,2,3]:  # No extrapolation BC for radial electric field on core boundary
            bbb.isutcore=isutcore  # Model for determining Phi at midplane
                                #=0 toroidal momentum = lzcore on core
                                #=1 d<uz>/dy=0
                                #>1 d^2(Ey)/dy^2=0 at OMP

    bbb.iphibcwi=iphibcwi      # Potential BC at PF wall:
                            #=0: d(ey)/dy=0
                            #=1: phi(ix,0) = phintewi*te(ix,0)/ev
                            #=3: d(phi)/dy/phi = 1/lyphi(1)
                            #=4: phi(ix,0)=phiwi(ix) in PF region
    bbb.iphibcwo=iphibcwo      # Potential BC at outer wall: as above




