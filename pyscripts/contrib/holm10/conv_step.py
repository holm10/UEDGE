# Holm10 Sep 10 2019, created from scratch


def conv_ncore_step(increment, name, var='ncore', ivar=0, stop=None, dtreal=1e-9, **kwargs):
    """ Simple function to incrementally change the density of a converged solution """
    from uedge import bbb
    from uedge.rundt import rundt
    from uedge.hdf5 import hdf5_save
    from copy import deepcopy
    

    if stop is None:
        try:
            stop = deepcopy(bbb.__getattribute__(var)[ivar])
        except:
            stop = deepcopy(bbb.__getattribute__(var))
        if increment < 0:
            stop /= 5
        else:
            stop *= 5


    # Setup dt-run
    bbb.t_stop=t_stop   # Set stop time
    bbb.ii1max=ii1max   # Set max outer loop iterations 

    while True:    
        bbb.dtreal=dtreal # Small step size to ensure convergence
        try:
            bbb.__getattribute__(var)[ivar] += increment
        except:
            _var = bbb.__getattribute__(var) 
            bbb.__setattr__(var, _var+increment)
       
#        bbb.ncore[0]+=increment # Increase step size
        bbb.label[0] = '{}_{}={:.3e}'.format(name, var, bbb.ncore[0])
        print('===================================')
        print('Solving for ncore[0]={:.2E}'.format(bbb.ncore[0]))
        print('===================================')
        '''
        bbb.exmain() # Check convergence at small dt
        # Check that we can get started
        if bbb.iterm!=1:    # The case did not converge
            bbb.isbcwdt=1       # Relax BC:s and try again
            bbb.exmain()            # Check convergence
            if bbb.iterm!=1:        # Case will not converge
                return "Case does not converge at dtreal=1e-9 w/ isbcwdt=1. Aborting..."

        # We have an initially converging case. Start dt-rund
        bbb.dt_tot=0        # Reset time
        if bbb.isbcwdt==1:  # We have relaxed BC:s - do a small step first
            # Advance to a micro-second
            bbb.t_stop=1e-5
            rundt(bbb.dtreal)
            # Set BCs and run to SS
            bbb.isbcwdt=0
            bbb.t_stop=t_stop
            rundt(bbb.dtreal)
        else:
            rundt(bbb.dtreal)
        '''
        bbb.dtreal = dtreal
        bbb.issfon=1
        bbb.ftol=1e-5
        bbb.exmain()
        rundt(dtreal=bbb.dtreal, **kwargs)
        # If run with BC relaxation, ensure convergence without relaxation
        
        # We should now have a steady-state solution: ensure this!
#        bbb.dtreal=1e20
#        bbb.itermx=30
#        bbb.icntnunk=0
#        bbb.issfon=0
#        bbb.ftol=1e-8
#        bbb.exmain()
        # Check if SS or not
        if bbb.iterm != 1:
            break
#        if bbb.iterm==1:   # SS solution
#            # Save to solutions with appropriate name
#            hdf5_save("../solutions/{}_{:.3E}_ss.hdf5".format(name,bbb.ncore[0]))
#        else:
#            hdf5_save("../solutions/{}_{:.2E}_failed.hdf5".format(name,bbb.ncore[0]))
#            print("Ramp did not converge for variable: {:.2E}. Aborting...".format(bbb.ncore[0]))
#            break
        try:
            currvar = deepcopy(bbb.__getattribute__(var)[ivar])
        except:
            currvar = deepcopy(bbb.__getattribute__(var))
        if increment > 0:
            if currvar > stop:
                break
        else:
            if currvar < stop:
                break






        

def conv_ivolcur_step(  d, name, t_stop=100,ii1max=100,ivolcurstop=100,
                        dtreal=1e-9,iisp=0,res=3):
    """ Simple function to incrementally change the density of a converged solution """
    from uedge import bbb
    from uedge.rundt import rundt
    from uedge.hdf5 import hdf5_save
    
    isbcwdt=bbb.isbcwdt

    # Check if we are increasing or decreasing!
    if d>0:
        increasing=True
    else:
        increasing=False



    # Setup dt-run
    bbb.t_stop=t_stop   # Set stop time
    bbb.ii1max=ii1max   # Set max outer loop iterations 

    while True:    
        bbb.dtreal=dtreal # Small step size to ensure convergence
        bbb.ivolcur[iisp]+=d # Increase step size
        print('===================================')
        print('Solving for ivolcur[0]={{:.{}E}}'.format(res).format(bbb.ivolcur[0]))
        print('===================================')
        bbb.exmain() # Check convergence at small dt
        '''
        # Check that we can get started
        if bbb.iterm!=1:    # The case did not converge
            bbb.isbcwdt=1       # Relax BC:s and try again
            bbb.exmain()            # Check convergence
            if bbb.iterm!=1:        # Case will not converge
                return "Case does not converge at dtreal=1e-9 w/ isbcwdt=1. Aborting..."

        # We have an initially converging case. Start dt-rund
        bbb.dt_tot=0        # Reset time
        if bbb.isbcwdt==1:  # We have relaxed BC:s - do a small step first
            # Advance to a micro-second
            bbb.t_stop=1e-5
            rundt(bbb.dtreal)
            # Set BCs and run to SS
            bbb.isbcwdt=0
            bbb.t_stop=t_stop
            rundt(bbb.dtreal)
        else:
            rundt(bbb.dtreal)
        '''
        rundt(dtreal=bbb.dtreal)
        # If run with BC relaxation, ensure convergence without relaxation
        if bbb.isbcwdt==1:
            bbb.isbcwdt=0
            rundt(dtreal=1e-6)

        
        # We should now have a steady-state solution: ensure this!
        bbb.dtreal=1e20
        bbb.itermx=30
        bbb.icntnunk=0
        bbb.ftol=1e-8
        bbb.exmain()
        # Check if SS or not
        if bbb.iterm==1:   # SS solution
            # Save to solutions with appropriate name
            hdf5_save("../solutions/{{}}_{{:.{}E}}_ss.hdf5".format(res).format(name,bbb.ivolcur[iisp]))
        else:
            hdf5_save("../solutions/{{}}_{{:.{}E}}_failed.hdf5".format(res).format(name,bbb.ivolcur[iisp]))
            print("Ramp did not converge for variable: {{:.{}E}}. Aborting...".format(res).format(bbb.ivolcur[iisp]))
            break
        if increasing:
            if bbb.ivolcur[iisp]>ivolcurstop:
                break
        else:
            if bbb.ivolcur[iisp]<ivolcurstop:
                break
        bbb.isbcwdt=isbcwdt






        

