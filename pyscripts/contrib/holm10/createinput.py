from uetools import Case



def createinput(nhgsp, nhsp, ngsp, nzsp, isimpon, gengrid, inputfunc, inputfile,
        outputfile='changed_vars.txt', casename='converted'):
    from uedge import bbb, com
    from uetools import watch
    from copy import deepcopy


    c=Case()

    _nhgsp = deepcopy(com.nhgsp)
    _nhsp = deepcopy(com.nhsp)
    _ngsp = deepcopy(com.ngsp)
    _nzsp = deepcopy(com.nzsp)

    com.nhgsp=nhgsp
    com.nhsp=nhsp
    com.ngsp=ngsp
    if isinstance(nzsp, int):
        com.nzsp[0]=nzsp
    elif isinstance(nzsp, int):
        com.nzsp.put(range(len(nzsp)), nzsp)
    else:
        print('nzsp type not recognized, aborting')
        return False
    bbb.isimpon=isimpon
    bbb.gengrid=gengrid
    bbb.allocate()
    
    com.nhgsp = _nhgsp
    com.ngsp = _ngsp
    com.nhsp = _nhsp
    com.nzsp = _nzsp
    
    # Create checkpoint
    checkpoint = watch.watchall()
    # Execute original input function
    inputfunc()
    # Retirieve changed variables
    changed = checkpoint.check()

    
    with open(inputfile,'r') as f:
        strinput = f.read()    

    inputvars = []
    for var in changed:
        if c.getpackage(var) in ['grd', 'flx']:
            pass
        else:
            if var in strinput:
                if var not in ['nis', 'ngs', 'tes', 'tis', 'phis', 
                        'ups', 'tgs', 'nysol', 'nycore', 'nxleg', 'nxcore',
                        'geqdskfname', 'aeqdskfname', 'aphdir', 'apidir']:
                    inputvars.append(var)

    c=Case()

    inputvars.append('isimpon')

    with open(outputfile, 'w') as f:
        f.write('casename: {}\n'.format(casename))
        f.write('\nspecies:\n')
        for var in ['nhgsp', 'nhsp', 'ngsp', 'nzsp']:
            f.write('    {}: {}\n'.format(var, (' '.join(str(c.get(var)).split()).replace(' ', ',').replace('[,','['))))
        f.write('\ngrid:\n')
        for var in ['mhdgeo', 'gengrid', 'geometry', 'isnonog']: 
            f.write('    {}: {}\n'.format(var, (' '.join(str(c.get(var)).split()).replace(' ', ',').replace('[,','['))))
        f.write('\nphysics:\n')
        for var in inputvars:
            if var not in ['nhgsp', 'nhsp', 'ngsp', 'nzsp','mhdgeo', 
                'gengrid', 'geometry', 'isnonog']: 
                f.write('    {}: {}\n'.format(var, (' '.join(str(c.get(var)).split()).replace(' ', ',').replace('[,','['))))
    return 

