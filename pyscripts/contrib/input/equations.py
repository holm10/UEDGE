# Equation package for UEDGE, equations.py
# Separated from 160299 master input by holm10
# Changelog
# 200213 - file created
from uedge import bbb,com,grd,flx


def inertial_atoms( ni=1, na=1, nm=1,ngimp=1,niimp=1,
                    ui=1, ua=1, uiimp=1,
                    te=1, ti=1, tm=1, 
                    phi=1
                   ):
    """====================================================================================================
    EQUATIONS SOLVED
    ===================================================================================================="""

    """----------------------------------------------------------------------------------------------------
    CONTINUITY EQUATION
    ----------------------------------------------------------------------------------------------------"""
    bbb.isnion[0]=              ni   # H+
    bbb.isnion[1]=              na   # Inertial H0
    bbb.isngon[0]=              0   # =0 if isupgon[H0]=1 -> Inertial H0 model
    if bbb.ishymol==1:
        bbb.isngon[1]=          nm   # Diffusive H2

    bbb.isngon[com.nhgsp:]=     ngimp       # Diffusive neutral C
    bbb.isnion[com.nhgsp:]=     niimp  # Impurity ions

    # 10.1.1: H+ continuity equation scalings
    #- - - - - - - - - - - - - - - - - -
    bbb.cnfx=1. # Coefficient for poloidal convection in ion continuity equation
    bbb.cnfy=1. # Coefficient for radial convection in ion continuity equation
    bbb.cnsor=1.    # Coefficient for particle source in ion continuity equation
    bbb.cnurn=1.    # Coefficient scaling nurlx (rate coefficient to relax BC) in ion continuity equation

    # 10.1.1: H0 continuity equation scalings
    #- - - - - - - - - - - - - - - - - -
    bbb.cngsor=1.   # Coefficient for particle source in gaseous continuity equation
    bbb.cnurg=1.    # Coefficient scaling nurlx (rate coefficient to relax BC) in gaseous continuity equation



    """----------------------------------------------------------------------------------------------------
    MOMENTUM EQUATIONS
    ----------------------------------------------------------------------------------------------------"""
    bbb.isupon[0]=              1   # H+
    bbb.isupon[1]=              1   # H0
    bbb.isupgon[0]=             1   # =1 if isnfon[H0]=0 -> inertial H0 model
    if bbb.ishymol==1:
        bbb.isupgon[1]=         0   # H2: diffusive
    bbb.isupgon[com.nhgsp:]=    0   # C: diffusive
    bbb.isupon[com.nhgsp:]=     uiimp  # Impuritiy ions

    # 10.2.1: Momentum equation options
    #- - - - - - - - - - - - - - - - - 
    bbb.cmfx=1. # Coefficient for poloidal convection in ion momentum equation
    bbb.cmfy=1. # Coefficient for radial convection in ion momentum equation
    bbb.cpgx=1. # Coefficient for pressure gradient in ion momentum equation
    bbb.cnuru=1.    # Coefficient scaling nurlx (rate coefficient to relax BC) in ion momentum equation


    # Verify that the H0 model is set correctly
    if (bbb.isupgon[0]==1) and (bbb.isngon[0]==1):
        print('H0 set up incorrectly!')
        print('Inertial H0: isupgon[0]=1 & isngon=[0]=0 - recommended')
        print('Diffusive H0: isupgon[0]=0 & isngon[0]=1 - use with care, not maintained!')
        return('Input file not read!')


    """----------------------------------------------------------------------------------------------------
    ENERGY EQUATION
    ----------------------------------------------------------------------------------------------------"""
    bbb.isteon=                 te   # Electrons
    bbb.istion=                 ti  # H+
    bbb.istgon[1]=              tm  # H2

    # 10.3.1: Energy equation options
    #- - - - - - - - - - - - - 
    bbb.cfvisx=1.   # Coefficient for poloidal viscosity in ion energy equation
    bbb.cfvisy=1.   # Coefficient for radia viscosity in ion energy equation
    bbb.cvgp=1. # Coefficient for v.Grad(P) in ion energy equation
    bbb.cnure=1.    # Coefficient scaling nurlx (rate coefficient to relax BC) in electron energy equation
    bbb.cnuri=1.    # Coefficient scaling nurlx (rate coefficient to relax BC) in ion energy equation




    """----------------------------------------------------------------------------------------------------
    POTENTIAL EQUATION
    ----------------------------------------------------------------------------------------------------"""
    bbb.isphion=                phi   # Potential equation
    bbb.isphiofft=              0   # Potential eq complement





