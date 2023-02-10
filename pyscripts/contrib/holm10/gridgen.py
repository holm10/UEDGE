
''' PLOT GRID VISUALIZATION '''

def plot_flx(first=None, last=None, ax=None):
    ''' Plots flux surfaces as defined in grid setup
        Based on plotflx.bas by Rensink/Rognlien/Porter 
    '''
    from uedge import com, grd, flx
    from matplotlib.pyplot import subplots

    flx.flxrun()
    if ax is None:
        f, ax = subplots(figsize=(7,9))

    if first is None:
        first = 0
    if last is None:
        last = 2*(com.nycore[0]+com.nysol[0]+com.nyout[0])+2
    first = max(0, first)

    # Plot vessel if present in EQDSK files
    if com.nlim > 0:
        ax.plot(com.xlim, com.ylim, 'k-', linewidth=2)

    # Plot target plates, if specified
    try:
        ax.plot(grd.rplate1, grd.zplate1, 'r-')
    except:
        pass
    try:
        ax.plot(grd.rplate2, grd.zplate2, 'r-')
    except:
        pass

    # Plot the flux surfaces within the sepcified range
    for i in range(first,last+1):
        # Plot SOL flux surfaces for each half-mesh
        if ((i >= com.jmin[0]-1) and (i < com.jsptrx[0])) or \
            ((i >= com.jsptrx[1]) and (i <= com.jmax[1])):
            ax.plot(com.xcurve[:,i][abs(com.xcurve[:,i])>0], 
                    com.ycurve[:,i][abs(com.ycurve[:,i])>0], 'k-', linewidth=0.3)
        # Plot CORE/PFR flux surfaces for each half-mesh
        elif ((i >= com.jsptrx[0]) and (i <= com.jmax[0])) or \
                ((i >= com.jmin[1]-1) and (i <= com.jsptrx[1]+1)):
            ax.plot(com.xcurve[:,i][abs(com.xcurve[:,i])>0][:flx.ijumpf[i]], 
                    com.ycurve[:,i][abs(com.ycurve[:,i])>0][:flx.ijumpf[i]], 
                    'k-', linewidth=0.3)
            ax.plot(com.xcurve[:,i][abs(com.xcurve[:,i])>0][flx.ijumpf[i]:], 
                    com.ycurve[:,i][abs(com.ycurve[:,i])>0][flx.ijumpf[i]:], 
                    'k-', linewidth=0.3)
            
        

    ax.set_aspect('equal')
    ax.set_xlabel('Horizontal position [m]')
    ax.set_ylabel('Vertical position [m]')

def plot_poloidal_distribution(yaxis='log', ylim = [1e-5, 5]):
    ''' Plots poloidal settings of model 
        Plot before executing grdrun?
    '''
    from uedge import com, flx, grd
    from matplotlib.pyplot import figure
    from numpy import linspace
    
    f = figure(figsize=(14,9))
    gs = f.add_gridspec(2,2, height_ratios=[6,1], wspace=0, hspace=0.3, top=0.98, bottom=0.05)
    axul = f.add_subplot(gs[0,0])
    axur_buf = f.add_subplot(gs[0,1])
    axur = axur_buf.twinx()
    axb = f.add_subplot(gs[1,:])

#    f, ax = subplots(2, 2, figsize=(14,9), gridspec_kw={'height_ratios': [5, 1]})
    x = linspace(0,1,200)

    if yaxis == 'lin':
        pli = axul.plot
        plo = axur.plot
    elif yaxis == 'log':
        pli = axul.semilogy
        plo = axur.semilogy
    else:
        print('yaxis option not recognized, using linear')
        pli = axur.plot
        plo = axur.plot
        
    

    flx.flxrun()
    grd.grdrun()

    nxtot = sum(com.nxleg[0])+sum(grd.nxuse)
    #for ix in range(nxtot):
    #    ax[0].axvline(ix/nxtot, color='grey', linewidth=0.5)
    
    if grd.kxmesh == 0:
        print('Manual mesh seeding interface TBI')
        return None
    elif grd.kxmesh == 1:
        func = grd.xfcn
    elif grd.kxmesh == 2:
        func = grd.xfcn2
    elif grd.kxmesh == 3:
        func = grd.xfcn3
    elif grd.kxmesh == 4:
        func = grd.xfcn4
    else:
        print('grd.kxmesh option not recognized, aborting')
        return None
    

    if grd.kxmesh != 4:
        pli(x, [func(a) for a in x], 'k-')
        pli([a/nxtot for a in range(nxtot)], [func(a/nxtot) for a in range(nxtot)], 'r.')
        plo(x, [func(1)-func(a) for a in x], 'k-')
        plo([a/nxtot for a in range(nxtot)], [func(1)-func(a/nxtot) for a in range(nxtot)], 'r.')
        plot_distro(func, ax=axb)
    else:
        pli(x, [func(a, nxtot) for a in x], 'k-')
        pli([a/nxtot for a in range(nxtot)], [func(a/nxtot, nxtot) for a in range(nxtot)], 'r.')
        plo(x, [func(1, nxtot)-func(a, nxtot) for a in x], 'k-')
        plo([a/nxtot for a in range(nxtot)], [func(1, nxtot)-func(a/nxtot, nxtot) for a in range(nxtot)], 'r.')
        plot_distro(func, nxtot, ax=axb)
   
    for ax in [axul, axur]:
        ax.axvline(com.nxleg[0,0]/nxtot, color='r', linewidth=3)
        ax.axvline(1-com.nxleg[0,1]/nxtot, color='r', linewidth=3)
        ax.set_ylim(ylim)
    
    axul.set_xlabel('Normalized distance between targets, index space')
    axur_buf.set_xlabel('Normalized distance between targets, index space')
    axur_buf.set_yticks([])
    axul.set_ylabel('Distance from inner target [m]')
    axur.set_ylabel('Distance from outer target [m]')
    axul.set_xlim(0,0.5)
    axur.set_xlim(0.5-1e-6,1)

def plot_distro(func, *args, ax=None, fname='',
            outpath='/Users/holm10/Documents/fusion/analysis/gridcomp_22/figs'):
    from uedge import com, flx, grd
    from matplotlib.pyplot import subplots


    if ax is None:
        f, ax = subplots(1, 1, figsize=(10,3))
    ax2 = ax.twinx()    

    nxtot = sum(com.nxleg[0])+sum(grd.nxuse)
    for ix in range(nxtot):
        x = func(ix/nxtot,*args)
        ax.axvline(func(ix/nxtot,*args), color='k', linewidth=1)
    ax.axvline(func(com.nxleg[0,0]/nxtot, *args), color='r', linewidth=3)
    ax.axvline(func(1-com.nxleg[0,1]/nxtot, *args), color='r', linewidth=3)

    ax.set_xlabel('Poloidal distance along separatrix [m]')
    ax.set_ylabel('Inner target')
    ax2.set_ylabel('Outer target')    
    ax.set_yticks([])
    ax2.set_yticks([])
    ax.set_title('Poloidal cell distribution along flux-tube')

    ax.set_ylim(0,1)
    ax.set_xlim(0,func(1,*args))
    
    return ax.get_figure()

def plot_efit(aeqdskfname='aeqdsk', geqdskfname='neqdsk', ax=None, ncontour=80):
    from matplotlib.pyplot import subplots
    from copy import deepcopy
    from numpy import linspace
    from uedge import com, flx
    from os.path import exists
    from scipy.interpolate import interp2d

    # Backup original pointers
    oldaeqdskfname = deepcopy(com.aeqdskfname)
    oldgeqdskfname = deepcopy(com.geqdskfname)
    # Set new file paths
    com.aeqdskfname = aeqdskfname
    com.geqdskfname = geqdskfname


    # Check whether the aeqdsk file can be located: if not, do not execute aeqdsk()
    if exists(aeqdskfname):
        flx.aeqdsk()

    if exists(geqdskfname):
        flx.neqdsk()
    else:
        print('EFIT geqdsk file "{}" not found.\nAborting...'.format(geqdskfname))
        return

    if ax is None:
        f, ax = subplots(figsize=(7,9))

    # Reconstruct EFIT grid
    x = linspace(0, com.xdim, com.nxefit)+com.rgrid1
    y = linspace(0, com.zdim, com.nyefit)-(com.zdim*0.5-com.zmid)

    # Interpolate on EFIT grid to find X-points
    interp = interp2d(x, y, com.fold.transpose())

    ax.contour(x, y, com.fold.transpose(), ncontour, colors='grey', 
        linewidths=0.5, linestyles='solid') #, negative_linestyles='solid')

    # Check whether the upper X-point exists
    if (com.rseps2 >= x.min()) and (com.rseps2 <= x.max()):
        if (com.zseps2 >= y.min()) and (com.zseps2 <= y.max()):
            upperxpoint = interp(com.rseps2, com.zseps2)
            ax.contour(x, y, com.fold.transpose(), [upperxpoint], colors='k', 
                    linewidths=1, linestyles='solid') #, negative_linestyles='solid')
    
    # Check whether the lower X-point exists
    if (com.rseps >= x.min()) and (com.rseps <= x.max()):
        if (com.zseps >= y.min()) and (com.zseps <= y.max()):
            lowerxpoint = interp(com.rseps, com.zseps)
            ax.contour(x, y, com.fold.transpose(), [lowerxpoint], colors='k', 
                    linewidths=1, linestyles='solid') #, negative_linestyles='solid')
    

    try:
        ax.plot(com.xlim, com.ylim, 'k-', linewidth=2)
    except:
        pass
#    ax.set_xlim(com.xold.min(), com.xold.max())
#    ax.set_ylim(com.yold.min(), com.yold.max())
    ax.set_aspect('equal')
    ax.set_xlabel('Horizontal position [m]')
    ax.set_ylabel('Vertical position [m]')
    ax.set_title(com.runid[0].decode('UTF-8').strip())

    # Restore original pointers
    com.aeqdskfname = oldaeqdskfname
    com.geqdskfname = oldgeqdskfname

    return ax.get_figure()
