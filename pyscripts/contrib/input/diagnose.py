# Diagnostic package for UEDGE, diagnostic.py
# Contains functions to diagnose the input deck
# Created from scratch by holm10
# Changelog
# 200213 - file created
from uedge import bbb,com,grd,flx

def diagnose():
    bbb.issfon=0;bbb.ftol=1e20;bbb.exmain();bbb.issfon=1;bbb.ftol=1e-8
    species()
    models()
    recycling()
    equations()
    

def species():
    ''' Script describing the plasma species present in the setup '''
    S=['X','He','Li','Be','B','C','N','O','F','Ne']
    S[0]='H'*(bbb.minu[0]==1)+'D'*(bbb.minu[0]==2)+'T'*(bbb.minu[0]==3)
    mu=bbb.mi[0]/bbb.minu[0] # Proton mass
    ng=list(dict.fromkeys(bbb.znucl)) # Gas species list

    print('')
    print('SPECIES NUMBER AND INDICES')
    print('++++++++++++++++++++++++++')
    print('')
    print('  NUMBER OF ION SPECIES: {}'.format(com.nisp))
    print('  Species index in zero-base (Python) indexing')
    for i in range(com.nisp):
        print('    {}: {}{}'.format(i,S[bbb.znucl[i]-1],''+str(int(bbb.zi[i-1])+1)*(bbb.zi[i]>1)+'+'*(bbb.zi[i]>0)))
    print('')
    print('  NUMBER OF GAS SPECIES: {}'.format(com.ngsp))
    print('  Species index in zero-base (Python) indexing')
    for i in range(com.ngsp):
        if i==1 and bbb.ishymol==1: # Molecules
            print('    {}: {}2'.format(i,S[0]))
        else:
            print('    {}: {}'.format(i,S[ng[i-bbb.ishymol*(i>0)]-1])) 

def models():

    S=['X','He','Li','Be','B','C','N','O','F','Ne']
    S[0]='H'*(bbb.minu[0]==1)+'D'*(bbb.minu[0]==2)+'T'*(bbb.minu[0]==3)
    ng=list(dict.fromkeys(bbb.znucl)) # Gas species list


    print('')
    print('PHYSICS MODELS')
    print('++++++++++++++')
    print('')
    print('  PLASMA MODEL')
    # Core BC: power/density
    # Bohm/momentum setup
    # Rate files
    # Ionization
    # CX
    # Recombination
    print('')
    print('  ATOM MODEL: {}'.format('INERTIAL'*(0 in bbb.zi)+'DIFFUSIVE'*(0 not in bbb.zi)))
    if bbb.isngon[0]==1:
        print('+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++')
        print('ERROR! Inertial atoms requested but diffusive model flag is on (bbb.isngon[0]==1)')
        print('+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++')
    if bbb.ishymol==1:
        print('')
        print('  MOLECULAR MODEL: {}'.format('CONSTANT-TEMPERATURE'*(bbb.istgcon[1]==1)+'ENERGY-EQUATION SOLVED'*(bbb.istgcon[1]==-1)))
        if bbb.istgcon[1]==1:
            print('    Min Tm={} eV'.format(bbb.tg[1:-1,1:-1,1].min()/bbb.ev))
            print('    Max Tm={} eV'.format(bbb.tg[1:-1,1:-1,1].max()/bbb.ev))
    if bbb.isimpon>0:
        print('')
        if bbb.isimpon==2: mod='FIXED FRACTION'
        elif bbb.isimpon==5: mod="HIRSCHMAN'S REDUCED-ION"
        elif bbb.isimpon==6: mod='FORCE-BALANCE'
        elif bbb.isimpon==7: mod='SIMULTAEOUS FIXED-FRACTION AND MULTI-CHARGE-STATE'
        else:
            print('+++++++++++++++++++++++++++++++++++++')
            print('ERROR! Invalid impurity model option!')
            print('+++++++++++++++++++++++++++++++++++++')
        print('  IMPURITY MODEL: {} MODEL'.format(mod))
        if bbb.isimpon in [2,7]:
            print('    Atomic impurity concentration: {}'.format(bbb.afrac.max()))

        for i in range(com.ngsp):
            if bbb.isch_sput[i]==0: 
                if bbb.chemsputi[i,:].max()==0 and bbb.chemsputo[i,:].max()==0:
                    continue
                else:
                    chsput='Old model, use chemsputi/o'
            elif bbb.isch_sput[i]==5: chsput='Roth, G-R'
            elif bbb.isch_sput[i]==6: chsput='Haasz 97'
            elif bbb.isch_sput[i]==7: chsput='Haasz 97+Davis (low E)'
            print('    {} chemical sputtering model: {} model'.format(S[ng[i-bbb.ishymol*(i>0)]-1],chsput))


        for i in range(com.nhgsp,com.ngsp):
            if bbb.isph_sput[i]==0: phsput='Old fixed case'
            elif bbb.isph_sput[i]==1: phsput='DIVIMP/JET physical sputtering fits'
            elif bbb.isph_sput[i]==2: phsput='DIVIMP/JET physical sputtering fits with hydrogen ion chemical sputtering'
            elif bbb.isph_sput[i]==3: phsput='DIVIMP/JET physical sputtering fits with hydrogen ion chemical and physical sputtering'
            print('    {} physical sputtering model: {}'.format(S[ng[i-bbb.ishymol*(i>0)]-1],phsput))

            if bbb.isi_sputw[i]>0 or bbb.isi_sputpf[i]>0: print('')

        for i in range(com.ngsp):
            if bbb.isi_sputw[i]>0:
                print('    OUTER WALL SPUTTRING OF {}'.format(S[ng[i-bbb.ishymol*(i>0)]-1]))
                if bbb.isi_sputw[i]==1: wsput='Physical ion sputtering'
                elif bbb.isi_sputw[i]==2: wsput='Chemical and physical ion sputtering'
                print('      Sputtering model: {}'.format(chsput))
            
            if bbb.isi_sputpf[i]>0:
                print('    PRIVATE-FLUX WALL SPUTTRING OF {}'.format(S[ng[i-bbb.ishymol*(i>0)]-1]))
                if bbb.isi_sputpf[i]==1: wsput='Physical ion sputtering'
                elif bbb.isi_sputpf[i]==2: wsput='Chemical and physical ion sputtering'
                print('      Sputtering model: {}'.format(chsput))

        #print('    

 #   bbb.rcxighg[ingc]=rcxighg       # Turn off imp0 + H+ --> imp+ + H+
    
 #   com.iscxfit=iscxfit           # C-ion + H0 CX model:
                                #=0: analytic forms in Braams' rate package
                                #=1: polynomial fit to C.F. Maggi curves (1997)
                                #=2: same as =1, except Z=1 has lower rate from Pigarov
        # TODO: Impurity puffs
        # TODO: CX
        # Bohm?

def recycling():
    S=['X','He','Li','Be','B','C','N','O','F','Ne']
    S[0]='H'*(bbb.minu[0]==1)+'D'*(bbb.minu[0]==2)+'T'*(bbb.minu[0]==3)
    ng=list(dict.fromkeys(bbb.znucl)) # Gas species list
    print('')
    print('RECYCLING MODEL')
    print('+++++++++++++++')

    print('')
    print('  PLATE RECYCLING')
    if bbb.ishymol==0:
        # No molecules, recycling is direct
        print('    Fraction of ions recycled as atoms: {0:.2f}'.format(bbb.recycp[0]))
    else:
        # Molecules, calculate recycling
        print('    Fraction of incident particles recycled: {0:.2f}'.format(bbb.recycp[1]))
        print('    Fraction of particles recycled as molecules: {0:.2f}'.format(mol_frac(bbb.recycp[0],bbb.recycp[1])))
    if bbb.isimpon>0:
        # Impurities are in the mix
        for i in range(com.nhgsp,com.ngsp):
            print('    Fraction of incident '+S[ng[i-bbb.ishymol]-1]+' ions recycled as '+S[ng[i-bbb.ishymol]-1]+' atoms: {0:.2f}'.format(bbb.recycp[i])) 
    print('')
    print('  WALL RECYCLING')
    if bbb.ishymol==0:
        # No molecules, recycling is direct
        print('    Fraction of ions recycled as atoms: {0:.2f}'.format(bbb.recycw[0]))
    else:
        # Molecules, calculate recycling
        print('    Fraction of incident particles recycled: {0:.2f}'.format(bbb.recycw[1]))
        print('    Fraction of particles recycled as molecules: {0:.2f}'.format(mol_frac(bbb.recycw[0],bbb.recycw[1])))
    if bbb.isimpon>0:
        # Impurities are in the mix
        for i in range(com.nhgsp,com.ngsp):
            print('    Fraction of incident '+S[ng[i-bbb.ishymol]-1]+' ions recycled as '+S[ng[i-bbb.ishymol]-1]+' atoms: {0:.2f}'.format(bbb.recycp[i])) 
    
    
    print('')
    print('ALBEDO-LIKE PUMPING OF ONE-SIDED MAXWELLIAN OF GAS')
    print('++++++++++++++++++++++++++++++++++++++++++++++++++')

    print('')
    print('  OUTER WALL')
    for i in range(bbb.nwsor): #Loop through all specified wall sources
        print('    SOURCE/SINK {}'.format(i+1))
        species=(S[0]+'2')*(bbb.igspsoro[i]==2)*(bbb.ishymol==1)+S[ng[bbb.igspsoro[i]-1]-1]*(bbb.igspsoro[i]!=2)
        print('      Source width, centered at {} m from the '.format(bbb.xgaso[i])+'left boundary'*(bbb.issorlb[i]==1)+'right boundary'*(bbb.issorlb[i]==0)+': {} m'.format(bbb.wgaso[i]))
        print('      Fraction of '+species+' one-sided Maxwellian removed: {0:.2f}'.format(1-bbb.albdso[i]))
        print('      '+species+' current: {0:.2f} A'.format(bbb.igaso[i]))


    print('')
    print('  INNER WALL')
    for i in range(bbb.nwsor): #Loop through all specified wall sources
        print('    SOURCE/SINK {}'.format(i+1))
        species=(S[0]+'2')*(bbb.igspsori[i]==2)*(bbb.ishymol==1)+S[ng[bbb.igspsori[i]-1]-1]*(bbb.igspsori[i]!=2)
        print('      Source width, centered at {} m from the '.format(bbb.xgasi[i])+'left boundary'*(bbb.issorlb[i]==1)+'right boundary'*(bbb.issorlb[i]==0)+': {} m'.format(bbb.wgasi[i]))
        print('      Fraction of '+species+' one-sided Maxwellian removed: {0:.2f}'.format(1-bbb.albdsi[i]))
        print('      '+species+' current: {0:.2f} A'.format(bbb.igasi[i]))

    print('')
    print('  LEFT PLATE')
    for i in range(com.nxpt):   # Loop over X-pts
        tab='    '
        if com.nxpt>1: 
            print('    X-POINT {}'.format(i+1))
            tab=tab+'  '
        for j in range(com.ngsp):
            species=(S[0]+'2')*(j==1)*(bbb.ishymol==1)+S[ng[j-1*(j>0)]-1]*(j!=1)
            print(tab+'Fraction of '+species+' one-sided Maxwellian removed: {0:.2f}'.format(1-bbb.albedolb[j,i]))

    print('')
    print('  RIGHT PLATE')
    for i in range(com.nxpt):   # Loop over X-pts
        tab='    '
        if com.nxpt>1: 
            print('    X-POINT {}'.format(i+1))
            tab=tab+'  '
        for j in range(com.ngsp):
            species=(S[0]+'2')*(j==1)*(bbb.ishymol==1)+S[ng[j-1*(j>0)]-1]*(j!=1)
            print(tab+'Fraction of '+species+' one-sided Maxwellian removed: {0:.2f}'.format(1-bbb.albedolb[j,i]))




def power_balance():
    print('TODO')
    # Core boundary
    # Walls
    # OT
    # IT
    # Rad
    # Other

def particle_balance():
    print('TODO')
    # Core boundary
    # Recycling
        # OT
        # IT
        # PFR wall
        # Outer wall

def currpot():
    print('TODO')
    # b0
    # Potential model


def equations():
    S=['X','He','Li','Be','B','C','N','O','F','Ne']
    S[0]='H'*(bbb.minu[0]==1)+'D'*(bbb.minu[0]==2)+'T'*(bbb.minu[0]==3)
    mu=bbb.mi[0]/bbb.minu[0] # Proton mass
    ng=list(dict.fromkeys(bbb.znucl)) # Gas species list

    print('')
    print('EQUATIONS EVOLVED')
    print('+++++++++++++++++')
    
    print('')
    print('  ENERGY EQUATIONS')
    
    print('    The ion energy is {}'.format('evolved'*(bbb.istion==1)+'frozen'*(bbb.istion==0)))
    print('    The electron energy is {}'.format('evolved'*(bbb.isteon==1)+'frozen'*(bbb.isteon==0)))

    for i in range(com.ngsp):
            if bbb.istgcon[i]==0:
                print('    {} atom temperature is taken to be {}*Ti'.format(S[ng[i-bbb.ishymol*(i>0)]-1],bbb.rtg2ti[i]))
            elif bbb.istgcon[i]==1:
                print('    {}'.format(S[ng[i-bbb.ishymol*(i>0)]-1])+' atom'*((bbb.ishymol!=1)+(i!=1))+'2'*(bbb.ishymol==1)*(i==1)+' temperature is taken to be {} eV'.format(bbb.tgas[i]))
            elif bbb.istgcon[i]==-1:
                if i==1 and bbb.ishymol==1:
                    print('    {}2 temperature is evolved'.format(S[ng[i-bbb.ishymol*(i>0)]-1]))
                else:
                    print('+++++++++++++++++++++++++++++++++++++++++')
                    print('ERROR! Invalid gas energy option istgcon!')
                    print('+++++++++++++++++++++++++++++++++++++++++')
                    
            else:
                print('    {} atom temperature is taken to be {}*Ti+{} eV'.format(S[ng[i-bbb.ishymol*(i>0)]-1],(1-bbb.istgcon[i])*bbb.rtg2ti[i],(bbb.istgcon[i]*bbb.tgas[i])))

    print('')
    print('  MOMENTUM EQUATIONS')
    for i in range(com.nisp):
        print('    {} momentum equation is {}'.format(S[bbb.znucl[i]-1]+str(int(bbb.zi[i-1])+1)*(bbb.zi[i]>1)+'+'*(bbb.zi[i]>0),'evolved'*(bbb.isupon[i]==1)+'frozen'*(bbb.isupon[i]==0)))

    print('')
    print('  CONTINUITY EQUATIONS')
    for i in range(com.nisp):
        print('    {} continuity equation is {}'.format(S[bbb.znucl[i]-1]+str(int(bbb.zi[i-1])+1)*(bbb.zi[i]>1)+'+'*(bbb.zi[i]>0),'evolved'*(bbb.isnion[i]==1)+'frozen'*(bbb.isnion[i]==0)))
    
    for i in range(com.ngsp):
        if i==0 and bbb.isupgon[0]==1:
            continue
        else:
            print('    {}'.format(S[ng[i-bbb.ishymol*(i>0)]-1])+' atom'*((bbb.ishymol!=1)+(i!=1))+'2'*(bbb.ishymol==1)*(i==1)+' continuity equation is {}'.format('evolved'*(bbb.isngon[i]==1)+'frozen'*(bbb.isngon[i]==0)))


    print('')
    print('  POTENTIAL EQUATION')
    print('    Potential equation is {}'.format('evolved'*bbb.isphion+'frozen'*(bbb.isphion==0)))

            
    




def mol_frac(Ra,Rm):
    ''' Calculates the molecular recycling fraction based on the recycling coefficients in UEDGE V7.08.04

        mol_frac(Ra,Rm)

        Variables:
        Ra:     Atomic recycling coefficient (recycp[0]/recycw[0])
        Rm:     Molecular recycling coefficient (recycp[1]/recycw[1])

        Returns:
        The fraction of returned flux (particles, not nuclei) that are molecules

        Note that these quantities only consider the incident ion fluxes, and does not consider
        albedo-like particle removal.

    '''
    return (Ra-Rm)/(Rm*(Ra-1)-2*Ra)

