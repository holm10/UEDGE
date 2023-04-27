# Python routine for reading/writing HDF5-type grid files
# Routines are called from basis in favor of Fortran routines
from h5py import File

# TODO: Check that files exist, raise errors, flag for calling xerrab

def read_gridpars(fname=None):
    from uedge import bbb, com, grd
    dblxpts = ['dnull', 'snowflake15', 'snowflake45', 'snowflake75', 
        'dnXtarget', 'isoleg']
    geometry = com.geometry[0].decode('UTF-8').strip()
    if fname is None:
        fname = bbb.GridFileName[0].decode('UTF-8').strip()
    gridue = File('{}.hdf5'.format(fname), 'r')
    gue = gridue['gridue']
    com.nxm = gue['nxm'][()]
    com.nym = gue['nym'][()]
    if geometry in dblxpts:
        com.iysptrx1 = gue['iysptrx1'][()]
        com.iysptrx2 = gue['iysptrx2'][()]
        com.ixlb = gue['ixlb'][()]
        com.ixpt1 = gue['ixpt1'][()]
        com.ixmdp = gue['ixmdp'][()]
        com.ixpt2 = gue['ixpt2'][()]
        com.ixrb = gue['ixrb'][()]
        if geometry == 'dnXtarget':
            com.nxc = com.ixmdp[0]
    else:
        com.ixpt1[0] = gue['ixpt1'][()]
        com.ixpt2[0] = gue['ixpt2'][()]
        com.iysptrx1[0] = gue['iysptrx1'][()]
        com.iysptrx2[0] = com.iysptrx1[0]
        com.ixlb[0] = 0
        com.ixrb[0] = com.nxm
    try:
        com.simagxs = gue['simagxs'][()]
    except:
        pass
    try:
        com.sibdrys = gue['sibdrys'][()]
    except:
        pass
    gridue.close()
        
    

def read_gridue(fname=None):
    from uedge import bbb, com, grd
    from Forthon import gchange
    if fname is None:
        fname = bbb.GridFileName[0].decode('UTF-8').strip()

    if bbb.iprint != 0:
        print(' Reading grid data from {}.hdf5'.format(fname))
    read_gridpars(fname)
    gchange('RZ_grid_info')
    gridue = File('{}.hdf5'.format(fname), 'r')
    gue = gridue['gridue']
    com.rm = gue['rm'][()]
    com.zm = gue['zm'][()]
    com.psi = gue['psi'][()]
    com.br = gue['br'][()]
    com.bz = gue['bz'][()]
    com.bpol = gue['bpol'][()]
    com.bphi = gue['bphi'][()]
    com.b = gue['b'][()]
    com.runid = gue['runid'][()]
    try:
        com.nlim = gue['nlim'][()]
        gchange('Comflxgrd')
        com.xlim = gue['xlim'][()]
        com.ylim = gue['ylim'][()]
    except:
        pass
    try:
        grd.nplate1 = gue['nplate1'][()]
        gchange('Mmod')
        grd.rplate1 = gue['rplate1'][()]
        grd.zplate1 = gue['zplate1'][()]
    except:
        pass
    try:
        grd.nplate2 = gue['nplate2'][()]
        gchange('Mmod')
        grd.rplate2 = gue['rplate2'][()]
        grd.zplate2 = gue['zplate2'][()]
    except:
        pass
    gridue.close()
    if bbb.iprint != 0:
        print(' Grid data read successfully:')
        print('     file name:   {}.hdf5'.format(fname))
        print('     run-ID:      {}'.format(com.runid[0].decode('UTF-8')))

def write_gridue(fname=None, runid=None):
    ''' Writes HDF5 grid file with name GridFileName '''
    from uedge import bbb, com, grd
    if fname is None:
        fname = bbb.GridFileName[0].decode('UTF-8').strip()
    if runid is None:
        runid = com.runid[0].decode('UTF-8').strip()
    gridue = File('{}.hdf5'.format(fname), 'w')
    gridue.require_group('gridue')
    gue = gridue['gridue']
    gue.create_dataset('nxm', data=com.nxm)
    gue.create_dataset('nym', data=com.nym)
    gue.create_dataset('rm', data=com.rm)
    gue.create_dataset('zm', data=com.zm)
    gue.create_dataset('psi', data=com.psi)
    gue.create_dataset('br', data=com.br)
    gue.create_dataset('bz', data=com.bz)
    gue.create_dataset('bpol', data=com.bpol)
    gue.create_dataset('bphi', data=com.bphi)
    gue.create_dataset('b', data=com.b)
    gue.create_dataset('runid', data=runid)
    
    if com.geometry[0].decode('UTF-8').strip() == 'dnull':
        gue.create_dataset('ixpt1', data=com.ixpt1)
        gue.create_dataset('ixpt2', data=com.ixpt2)
        gue.create_dataset('iysptrx1', data=com.iysptrx1)
        gue.create_dataset('iysptrx2', data=com.iysptrx2)
        gue.create_dataset('ixlb', data=com.ixlb)
        gue.create_dataset('ixmdp', data=com.ixmdp)
        gue.create_dataset('ixrb', data=com.ixrb)
    else:
        gue.create_dataset('ixpt1', data=com.ixpt1[0])
        gue.create_dataset('ixpt2', data=com.ixpt2[0])
        gue.create_dataset('iysptrx1', data=com.iysptrx1[0])
        
    
    # Store extra data, such as limiter and plate data
    try:
        gue.create_dataset('simagxs', data=com.simagxs)
    except:
        pass
    try:
        gue.create_dataset('sibdrys', data=com.sibdrys)
    except:
        pass
    try:
        gue.create_dataset('nlim', data=com.nlim)
        gue.create_dataset('xlim', data=com.xlim)
        gue.create_dataset('ylim', data=com.ylim)
    except:
        pass
    try:
        gue.create_dataset('nplate1', data=grd.nplate1)
        gue.create_dataset('rplate1', data=grd.rplate1)
        gue.create_dataset('zplate1', data=grd.zplate1)
    except:
        pass
    try:
        gue.create_dataset('nplate2', data=grd.nplate1)
        gue.create_dataset('rplate2', data=grd.rplate1)
        gue.create_dataset('zplate2', data=grd.zplate1)
    except:
        pass
     

    gridue.close()
    if bbb.iprint != 0:
        print(' Wrote grid file successfully:')
        print('     file name:   {}.hdf5'.format(fname))
        print('     run-ID:      {}'.format(runid))
