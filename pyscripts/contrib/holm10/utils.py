# Utilities created to help emulate BASIS capabilities
# Created from scratch by holm10
#   191121 - Created 
from uedge import *

def readcase(casedir,dirpath='data'):
    '''
    Simple script reading the requested input file from folder casedir

    readcase(casedir,**keys)

    Required parameters:
    casedir     -   The name of the case folder containin dirpath/input.py,
                    or alternatively the path to this folder
    
    Keyword parameters:
    dirpath     -   Path to input.py within casedir, default is "data"
    '''
    from os import chdir,getcwd
    from os.path import exists, abspath, relpath, isabs
    from importlib import reload
    import sys

    # Check that path is set up correctly
    if exists(casedir) is False:
        print('Case folder {} not found!'.format(casedir))
        return
    elif exists('{}/{}'.format(casedir,dirpath)) is False:
        print('Input folder{}/{} not found!'.format(casedir,dirpath))
        return
    
    cwd=getcwd()                                        # Store original location
    chdir('{}/{}'.format(casedir,dirpath))  # Move to the requested location
    import input as i                              # This is ugly coding, but it is straightforward and works
    reload(i)                                           # Reload case: should not strictly be necessary
    i.restore_input()                                              # Read the new setup 

    aphabs=isabs(aph.aphdir[0].decode('utf-8'))
    apiabs=isabs(api.apidir[0].decode('utf-8'))


    # Fix the path dependencies to be absolute paths for reading case out of data
    if aphabs is False:
        newaphdir='{}/{}'.format(abspath(getcwd()),aph.aphdir[0].decode('utf-8'))
    if apiabs is False:
        newapidir='{}/{}'.format(abspath(getcwd()),api.apidir[0].decode('utf-8'))


    # Set grid path file in V7.08 and later
    try:
        newaeqdskfname='{}/{}'.format(abspath(getcwd()),com.aeqdskfname[0].decode('utf-8'))
        newgeqdskfname='{}/{}'.format(abspath(getcwd()),com.geqdskfname[0].decode('utf-8'))
    except:
        pass

    chdir(cwd)                                          # Go back to the original directory to enable recursive loads 

    # Shorten the concatenated path to fit into arrays
    if aphabs is False:
        aph.aphdir[0]=abspath(newaphdir)
    if apiabs is False:
        api.apidir[0]=abspath(newapidir)

    try:
        com.aeqdskfname[0]=abspath(newaeqdskfname)
        com.geqdskfname[0]=abspath(newgeqdskfname)
    except:
        pass



def write_csv(varlist,fname,header=None,probname='',descriptor=''):
    ''' Writes a list of 1D profiles to a csv file '''
    from numpy import savetxt,zeros 
    from time import ctime    

    # Set problem name if not given
    if probname is None:
        probname=fname


    # Find longest array to be written
    maxlen=0
    for col in varlist:
        maxlen=max(maxlen,len(col))
    
    # Create an array to accommodate all variables to be written
    vals=zeros((maxlen,len(varlist)))

    # Set up header
    headerp='Problem name: {}\n{}\n{}\n'.format(probname,ctime(),descriptor)    

    # Populate the array and header
    for ii in range(len(varlist)):
        vals[:len(varlist[ii]),ii]=varlist[ii]
        if header is not None:
            headerp+=(header[ii]+', '*(ii!=len(varlist)-1)).ljust(14)
    
    # Save to CSV file
    savetxt('{}.csv'.format(fname),vals,fmt='%.5E',header=headerp,comments='',delimiter=' , ')


def write_radial_csv(varlist, fname,probname=None,  location='mp',gc=False):
    ''' Writes a CSV file with radial data along columns 

    write_radial_csv(varlist, **keys)

    Required parameters:
    varlist     -   2D list containing the data to be written in the following form:
                    [   [title1,    data1],
                        [title2,    data2],
                        ...
                        [titleN,    dataN] ]
                    Here, title is a string and dataX is the data to be written.
                    The data can be a 1D or 2D array:
                        -1D arrays are directly written to the column
                        -The radial profile specified by location is written to the column

    Keyword parameters:
    probname    -   Problem name to be printed at first row, defaults to bbb.label
    fname       -   File name for printed file (script adds .csv)
    location    -   What location to pic from 2D arrays. The options are:
                        'mp' - midplane [default]
                        'ot' - outer target
                        'it' - inner target
                        int  - integer specifying index
    gc          -   Boolean whether to include guard cells or not
                        False - omits guard cells [default]
                        True  - includes guard cells

    ************************************
    WARNING! Fluxes are taken at the RHS boundary of the evaluated cell, e.g. out
    of the cell in front of the inner target for location='it'. To get the right
    fluxes, specify the flux as a 1D array at the appropriate cell index.
    ************************************
    '''
    from uedge import bbb,com
    from numpy import zeros,savetxt

    if probname==None:
        probname=bbb.label[0].decode('UTF-8')
    
    # Set the poloidal index
    if location=='mp':
        ix=bbb.ixmp
        descriptor='Midplane values'
    elif location=='ot':
        ix=com.nx
        descriptor='Outer target values'
    elif location=='it':
        ix=1
        descriptor='Inner target values'
    elif isinstance(location,int):
        ix=location
        descriptor='Row index {} values'.format(ix)
    else:
        print('Location "{}" not recognized!'.format(location))
        print('Use "mp", "ot", "it" or the row integer.')
        return

    # Initialize array for values
    vals=zeros((len(varlist),com.ny+2))
    header='Problem name: {}\n{}\n'.format(probname,descriptor)    
    

    # Loop throug list of titles and values
    for i in range(len(varlist)):
        if len(varlist[i][1].shape)>2:
            print('Variable {} has more than two dimensions!'.format(varlist[i][0]))
            print('Only supply 1D or 2D arrays. Aborting...')
            return
        # Append header
        header+=(varlist[i][0]+', '*(i!=len(varlist)-1)).ljust(14)
        # Store values
        if len(varlist[i][1].shape)==2:   # Store index from 2D profile
            vals[i,:]=varlist[i][1][ix,:]
        elif len(varlist[i][1].shape)==1: # 1D array, store full array
            vals[i,:]=varlist[i][1]

    # Drop guard cells unless they are requested
    if gc is False:
        vals=vals[:,1:-1]
    
    # Save to CSV file
    savetxt('{}.csv'.format(fname),vals.transpose(),fmt='%.5E',header=header,comments='',delimiter=' , ')


def default_csv_parameters():
    '''
    Default parameters to be written, based on MGroth file
    '''
    from uedge import bbb,com
    return  [   ['psin',    com.psinormc],
                ['dyc',     1/(com.gy+1e-30)],
                ['dyf',     1/(com.gyf+1e-30)],
                ['ne',      bbb.ne],
                ['ni',      bbb.ni[:,:,0]],
                ['te',      bbb.te],
                ['ti',      bbb.ti]    ]


def quickconverge():
    ''' Converges case w/o Jac eval '''
    
    bbb.issfon=0;bbb.ftol=1e20;bbb.exmain();bbb.issfon=1;bbb.ftol=1e-8



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

def calc_recyc(fmol,frem):
    ''' Calculates the atomic and molecular recycling coefficients to achieve fmom and frem in V7.08.04

        cacl_recyc(fmom,frem)
    
        Variables:
        fmom:   The fraction returned flux (particles, not nuclei) that are recycled as molecules
        frem:   The fraction incident ion flux that is removed
    
        Returns:
        Ra,Rm
        Ra:     Atomic recycling coefficient (recycp[0]/recycw[0])
        Rm:     Molecular recycling coefficient (recycp[1]/recycw[1])

    '''
    from numpy import sqrt    

    # Perturbate around node to avoid singularity
    fmol+=1e-12

    Ra=((fmol-1)*(frem-1))/(fmol+1)
    Rm=(2*fmol-2*fmol*frem)/(2*fmol+frem-frem*fmol)

    return Ra,Rm


def reset():
    ''' File that resets UEDGE to the default setup '''
    import h5py
    from numpy import array
    from uedge import com,aph,api,bbb,flx,grd,svr,wdf
    hf=h5py.File('/home/holm10/uedge/uescripts/reset.hdf5')
    # Set species indices
    com.nhsp=hf['com']['nhsp'][()]
    com.nhgsp=hf['com']['nhgsp'][()]
    com.ngsp=hf['com']['ngsp'][()]
    com.nzsp=hf['com']['nzsp'][()]
    com.nxleg=hf['com']['nxleg'][()]
    com.nxcore=hf['com']['nxcore'][()]
    com.nycore=hf['com']['nycore'][()]
    com.nysol=hf['com']['nysol'][()]
    bbb.mhdgeo=hf['bbb']['mhdgeo'][()]
    aph.mpe=hf['aph']['mpe'][()]
    aph.mpd=hf['aph']['mpd'][()]
    aph.mpr=hf['aph']['mpr'][()]
    bbb.igspsori=hf['bbb']['igspsori'][()]
    bbb.igspsoro=hf['bbb']['igspsoro'][()]
    bbb.isimpon=hf['bbb']['isimpon'][()]
    # Allocate the correct size of all arrays
    bbb.allocate()
    bbb.ishymol=0
    # Restore all values
    #bbb.issfon=0;bbb.ftol=1e20;exmain();bbb.issfon=1;bbb.ftol=1e-8
    bbb.gchange('Impurity_source',0)

    for p in hf.keys():
        for v in hf[p].keys():
            if hf[p][v].dtype in ['int32','int64','float32','float64']:
                if v=='del': continue
                if len(hf[p][v].shape)==0:
                    if array(hf[p][v]).max()==0:
                        exec('{}.{}=0'.format(p,v))
                    elif array(hf[p][v]).max()==1:
                        exec('{}.{}=1'.format(p,v))
                    elif array(hf[p][v]).max()==300:
                        exec('{}.{}=300'.format(p,v))
                    else:
                        a=hf[p][v][()]
                        exec('{}.{}=a'.format(p,v)) 
                    #exec('{}.{}=hf[{}][{}][()]'.format(p,v,p,v))
                else:
                    if array(hf[p][v]).max()==0:
                        exec('{}.{}=0'.format(p,v))
                    elif array(hf[p][v]).max()==1:
                        exec('{}.{}=1'.format(p,v))
                    elif array(hf[p][v]).max()==300:
                        exec('{}.{}=300'.format(p,v))
                    else:
                        a=hf[p][v][()]
                        exec('{}.{}=a'.format(p,v)) 
                    #exec('{}.{}=array(hf[{}][{}])'.format(p,v,p,v))

    # Restore all values
    bbb.issfon=0;bbb.ftol=1e20;exmain();bbb.issfon=1;bbb.ftol=1e-8


def div(varx,vary,ind=None):
    ''' Calculates the divergence of var in the numerical domain '''
    from numpy import zeros

    if len(varx.shape)!=len(vary.shape):
        print('Shape mismatch between varx and vary! Aborting...')
        return
    
    if len(varx.shape)==3:
        if ind is not None:
            varx=varx[:,:,ind]
            vary=vary[:,:,ind]
        else:
            print('Multi-dimensional array without specified 3rd dimension!')
            print('Calculating for 3rd dimension index 0')
            varx=varx[:,:,0]
            vary=vary[:,:,0]

    # TODO: FIX THIS
    (nx,ny)=bbb.te.shape
    ret=zeros(bbb.te.shape)
    for iy in range(1,ny):
        for ix in range(nx):
            ret[ix,iy]=-varx[ix,iy]+varx[bbb.ixm1[ix,iy],iy]+vary[ix,iy-1]-vary[ix,iy]
    return -ret



def create_EIRENE(path='.',subpath='data'):
    ''' Creates a database
        Parameters:
            savename        If set, saves dict as pickle named 'savename'
            sortlocation    Location for sorting by te: 'core' or 'mp'
            outpath         Path to location where pickle is saved
            path            Path to parent directory
            subpath         Path to input.py within child directories of path: 
                            path/*/supath/input.py
            commands        List of commands to be executed before restoring solution
            variables       List of all variable names, including package, to be stored
            ret             Boolean whether to return dict or not
    '''
    from os import getcwd,chdir,remove,walk
    from os.path import abspath  
    from uedge import bbb,com
    from importlib import reload
    from uedge.contrib.holm10.eireneinput import writeEIRENE



    def natsort(l): 
        from re import split
        convert = lambda text: int(text) if text.isdigit() else text.lower() 
        alphanum_key = lambda key: [ convert(c) for c in split('([0-9]+)', key) ] 
        return sorted(l, key = alphanum_key)

    chdir(path)                 # Go to the parent directory
    parent=getcwd()               # Get path to parent
    # Get list of subdirectories in path
    dirs=natsort(next(walk(path))[1])
    # Omit supporting file directories
    try:
        dirs.remove('grid')
    except:    
        pass
    try:
        dirs.remove('rates')
    except:
        pass
    try:
        dirs.remove('ignore')
    except:
        pass
    
    if len(dirs)==0:
        return 'No directories found! Aborting...'


    loop=0
    for child in dirs:      # Loop through all directories
        loop+=1


        # TODO: exmain every Nth step
        print('******************************')
        print('*** Directory: '+child+' ***')
        print('******************************')

        chdir('{}/{}'.format(child,subpath))
        
        import input as i
        reload(i)
        i.restore_input()


        # Read and repopulate all arrays
        bbb.issfon=0;bbb.ftol=1e20;bbb.exmain()
        writeEIRENE(verbatim=False)


        chdir(parent)
    

def interpolate(x, y, z):
    from scipy.interpolate import griddata
    return interp2d(x,y,z)
    

def increase_resolution(new_nxleg, new_nxcore, new_nysol, new_nycore):
    ''' WORK IN PROGRESS '''
    from uedge import bbb, com
    from copy import deepcopy    
    import matplotlib.pyplot as plt

    old_ftol = deepcopy(bbb.ftol)
    bbb.issfon=0; bbb.ftol=1e20
    exmain() # Generate jacobian
    bbb.issfon=1 # Turn on jacobian evaluation

    points_old = []
    for ix in range(com.nx+2):
        for iy in range(com.ny+2):
            points_old.append([com.rm[ix,iy,0], com.zm[ix,iy,0]])
 

    interp_vars = []
    restore_vars = [bbb.tes]
    for var in restore_vars:
        buff = []
        for ix in range(com.nx+2):
            for iy in range(com.ny+2):
                buff.append(var[ix,iy])


    com.nxleg[0] = new_nxleg
    com.nxcore[0] = new_nxcore
    com.nysol[0] = new_nysol
    com.nycore[0] = new_nycore

    # Turn off molecular temperature eq interpolation - not yet implemented
    if bbb.ishymol == 1:
        istgon_old = deepcopy(bbb.istgon)
        istgcon_old = deepcopy(bbb.istgcon)
        bbb.istgcon = 0
        bbb.istgon = 0

    exmain() # Interpolate locally

    # Restore molecular model 
    if bbb.ishymol:
        bbb.istgon = istgon_old
        bbb.istgcon = istgcon_old

    
        

    bbb.ftol = old_ftol
#    return interp_vars








