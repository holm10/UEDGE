# Interpolates grid 

def interpolate(oldgrid, newgrid, savefile, geometry='snull'):
    from uedge import bbb, com, api, aph
    from uedge.hdf5 import hdf5_restore, hdf5_save
    from uedge.contrib.input.impurities import carbon_forcebalance
    from copy import deepcopy
    # set iprint
    if geometry in ['snull', 'uppersn']:
        bbb.mhdgeo=1
        com.geometry="snull"
        # Extract grid information
        with open(oldgrid) as f:
            [nxo, nyo, ixpt1o, ixpt2o, iysptrxo] = [int(x) for x \
                in f.readline().strip().split()]
        with open(newgrid) as f:
            [nxn, nyn, ixpt1n, ixpt2n, iysptrxn] = [int(x) for x \
                in f.readline().strip().split()]
        
        com.nxleg[0] = [ixpt1o, nxo - ixpt2o]
        ncore = ixpt2o - ixpt1o
        com.nxcore[0] = [round(ncore/2), ncore - round(ncore/2)]
        com.nycore[0] = iysptrxo
        com.nysol[0] = nyo - iysptrxo
        # Create grid
        # TODO: figure out how to do this
        gridpath = '/Users/holm10/Documents/fusion/uedge/runs/d3d/160299/2240/uegen/grid/'
        com.aeqdskfname=gridpath+"aeqdsk"
        com.geqdskfname=gridpath+"neqdsk" 
        com.nhsp=2
        carbon_forcebalance(apipath='/Users/holm10/Documents/fusion/uedge/runs/d3d/186841/rates/api')

        # Recreate old grid
        issfon_o = deepcopy(bbb.issfon)
        ftol_o = deepcopy(bbb.ftol)
        bbb.issfon=0
        bbb.ftol=1e20
        bbb.exmain()
    
        # Restore data to interpolate
        hdf5_restore(savefile)
        bbb.restart=1

        # Switch to new grid dimentions
        com.nxleg[0] = [ixpt1n, nxn - ixpt2n]
        ncore = ixpt2n - ixpt1n
        com.nxcore[0] = [round(ncore/2), ncore - round(ncore/2)]
        com.nycore[0] = iysptrxn
        com.nysol[0] = nyn - iysptrxn
        
        # Interpolate using UEDGE interpolation
        bbb.issfon=1
        bbb.exmain()
         
        bbb.issfon = issfon_o
        bbb.ftol = ftol_o 

        
        savefile = savefile.split('/')[-1]
        newname = savefile[:-5] + '_interp{}x{}'.format(nxn,nyn) + savefile[-5:]
        hdf5_save(newname)
        print('Wrote interpolated solution "{}"'.format(newname))
        
