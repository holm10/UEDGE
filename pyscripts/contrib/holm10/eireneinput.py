def split(arr, N):
    return [arr[x:x+N] for x in range(0, len(arr), N)]



def fortranFormat(n,d=6):
    sign = '-'*(n<0)
    n=abs(n)
    a = '{{:.{}E}}'.format(d).format(float(n))
    if n==0:
        return '{{:.{}E}}'.format(d+1).format(float(n))
    e = a.find('E')
    eE = a[e:e+1]+a[e+1:e+2]*(abs(int(a[e+1:])*1+1)!=0)+'+'*(abs(int(a[e+1:])*1+1)==0)
    return sign+'0.{}{}{}{:02d}'.format(a[0],a[2:e],eE,abs(int(a[e+1:])*1+1))




def sfmt(n,d=6,nz=6):
    '''
    d (int) : decimal places shown
    nz (int) : non-zero decimal places shown
    '''
    sign = '-'*(n<0)
    n=abs(n)
    a = '{{:.{}E}}'.format(nz).format(float(n))
    if n==0:
        return '{{:.{}E}}'.format(d+1).format(float(n))
    e = a.find('E')
    eE = a[e:e+1]+a[e+1:e+2]*(abs(int(a[e+1:])*1)!=0)+'+'*(abs(int(a[e+1:])*1)==0)
    return sign+'{}.{}{}{}{:02d}'.format(a[0],a[2:2+nz],'0'*(d-nz+1),eE,abs(int(a[e+1:])*1))




def sonnet_gc(x,y,d=1e-7, decimals=6):
    from uedge import com
    from numpy import around

    x,y = x.copy(),y.copy()
    for ix in range(1, com.nx+1):    
        for p in [   (0,1,1,1,-d), (0,2,1,2,-d), (0,4,1,2,0), 
                     (0,3,1,1,0), (-1,1,-2,1,0), (-1,2,-2,2,0),
                     (-1,4,-2,2,d), (-1,3,-2,1,d)]:
            x[ix, p[0], p[1]] = x[ix, p[2], p[3]] 
            y[ix, p[0], p[1]] = y[ix, p[2], p[3]] + p[-1]


    for iy in range(1, com.ny+1):
        for p in [  (0,1,1,1,-d), (0,2,1,1,0), (0,4,1,3,0),
                    (0,3,1,3,-d), (-1,1,-2,1,0), (-1,2,-2,1,d),
                    (-1,4,-2,3,d), (-1,3,-2,3,0)
                    ]:

            x[p[0], iy, p[1]] = x[p[2], iy, p[3]] + p[-1]
            y[p[0], iy, p[1]] = y[p[2], iy, p[3]] 

    for p in [  (0,0,1,0,1,1), (0,0,2,1,0,1), (0,0,4,1,1,1), (0,0,3,0,1,1),
                (0,-1,1,0,-2,3), (0,-1,2,1,-2,1), (0,-1,4,1,-1,3), (0,-1,3,0,-2,1),
                (-1,0,1,-2,0,2), (-1,0,2,-1,1,2), (-1,0,4,-1,1,1), (-1,0,3,-1,1,1),
                (-1,-1,1,-2,-2,3), (-1,-1,2,-1,-2,4), (-1,1,4,-1,-2,4), (-1,-1,3,-2,-1,4)]:
        x[p[0], p[1], p[2]] = x[p[3], p[4], p[5]]

    for p in [  (0,0,1,1,0,1), (0,0,2,1,0,1), (0,0,4,1,1,1), (0,0,3,0,1,1),
                (0,-1,1,0,-2,3), (0,-1,2,1,-1,1), (0,-1,4,1,-1,3), (0,-1,3,0,-2,1),
                (-1,0,1,-2,0,2), (-1,0,2,-2,0,1), (-1,0,4,-1,1,1), (-1,0,3,-1,1,1),
                (-1,-1,1,-2,-2,3), (-1,-1,2,-1,-2,4), (-1,-1,4,-2,-1,4), (-1,-1,3,-2,-1,4)]:
        y[p[0], p[1], p[2]] = y[p[3], p[4], p[5]]

      
    return x, y
      


def write_cell(i, ix, iy, x, y, f):
    from uedge import com

#      x3 ----------- x2
#      |              |
#      |              |
#      x4 ----------- x1

    d=9
    f.write('\n   ELEMENT{:>4d} = ({:>4d},{:>4d}):'.format(i, ix, iy))
    f.write(' ({:>18},{:>18})    '.format(sfmt(x[3],d), sfmt(y[3],d)))
    f.write(' ({:>18},{:>18})\n'.format(sfmt(x[4],d), sfmt(y[4],d)))
    f.write('   FIELD RATIO = {:>18}'.format(sfmt(0,d)))
    f.write(' '*10+'({:>18},{:>18})\n'.format(sfmt(0.25*sum(x[1:]),d),sfmt(0.25*sum(y[1:]),d)))
    f.write(' '*30 +'({:>18},{:>18})'.format(sfmt(x[1],d), sfmt(y[1],d)))
    f.write('     ({:>18},{:>18})\n'.format(sfmt(x[2],d=d), sfmt(y[2],d=d)))
    f.write('--------------------------------')
    return



def getxy(rm, zm, flip=False, dispr=0, dispz=0):
    from uedge import com, bbb

    y, x = rm.copy(), zm.copy()
    
    if (bbb.mhdgeo==1) or (flip is True):
        x, y = y, x
        if flip is True:
            y = abs(y - y[1:-1,1:-1].max())
            
    x += dispr
    y += dispz

    return x, y




def writeSonnet_slab(fname='uedge.sonnet', runid=" UEDGE geometry data", flip=False, dgc=1e-7, d=6, sonnetGC=True, dispr=0, dispz=0, rm=None, zm=None):
      
    from uedge import bbb, com

    if rm is None:
        rm = com.rm
    if zm is None:
        zm = com.zm

    x, y = getxy(rm, zm, flip, dispr, dispz)

    with open(fname,'w') as f:
        f.write('{:<80}\nNO. CUT   0\n\n\n================================='.format(runid))

    
        i = 0
        for ix in range(com.nx+2):
            for iy in range(com.ny+2):
                i+=1
                write_cell(i, ix, iy, x[ix,iy], y[ix,iy], f)

    print(" *** Sonnet slab geometry file written for EIRENE ***")






def write30(rm=None, zm=None, fname='fort.30', runid='UEDGE geometry data', flip=False, dispr=0, dispz=0):
    from uedge import com, bbb
    from uedge.contrib.holm10.plot import grid,plot

    if rm is None:
        rm = com.rm
    if zm is None:
        zm = com.zm

    order = getxy(rm, zm, flip, dispr, dispz)

    plot(order[0][:,:,0],order[1][:,:,0],line='',marker='.')

    nncut=0+2*(bbb.mhdgeo==1)
    with open(fname, 'w') as f:
        f.write(runid+'\n\n')
        f.write('{:>21}{:>21}{:>21}\n'.format(com.nx, com.ny, nncut))
        if bbb.mhdgeo==1:
            f.write('{:>21}{:>21}{:>21}{:>21}'.format(com.ixpt1[0], com.ixpt2[0]+1, 0, com.iysptrx1[0]))
            f.write('{:>21}{:>21}{:>21}{:>21}'.format(com.ixpt2[0], com.ixpt1[0]+1, 0, com.iysptrx1[0]))
        f.write('\n\n')
        for ix in range(1,com.nx+1):
            for iy in range(1,com.ny+1):
                for o in order:
                    for n in [2, 4, 3, 1]:
                        f.write('{:>15}'.format(fortranFormat(o[ix,iy,n])))
                    f.write('\n')
    print(" *** geometry file written for EIRENE ***")
        


def write31(fname='fort.31'):
    from uedge import bbb, com


    def gfsub3(var,f):
        from uedge import com
        eps=1e-90
        lim = int(((com.nx+2)/5))*5 - 4
        for iy in range(com.ny+2):
            for ix in range(0,lim+1,5):
                for iii in range(1,6):
                    f.write('{:>16}'.format(fortranFormat(var[ix-1+iii,iy]+eps,d=7)))
                f.write('\n')
            if ( lim+4 < com.nx+2 ):
                for ix in range(lim+5, com.nx+3):
                    f.write('{:>16}'.format(fortranFormat(var[ix-1,iy]+eps,d=7)))
                f.write('\n')

    with open(fname, 'w') as f:
        for var in [bbb.ni, bbb.uu, bbb.vy, bbb.te, bbb.ti, bbb.pr, bbb.upi,
                    com.rr, bbb.fnix, bbb.fniy, bbb.feix, bbb.feiy, bbb.feex,
                    bbb.feey, com.vol, com.b]:
            if len(var.shape)==3: # Three-dimensional, test for charge
                for i in range(com.nisp):
                    if bbb.zi[i]>0:
                        gfsub3(var[:,:,i], f)
            else:
                gfsub3(var, f)
    print(" *** background plasma file written for EIRENE ***")



def write_neighbor(fname='triageom.neighbor'):
    from uedge import com
    from numpy import reshape, zeros

    N = com.nx*com.ny*2
    indar = reshape(list(range(1, (com.nx)*(com.ny)+1)), (com.nx, com.ny), 'F')

    neighb = zeros((com.nx*com.ny*2,3))
    neighs = zeros((com.nx*com.ny*2,3))
    neighr = zeros((com.nx*com.ny*2,3))
    

    NTRIA1=0
    for iy in range(com.ny):
        for ix in range(com.nx):
            # Shared side
            neighb[2*indar[ix,iy]-2+NTRIA1,0] = 2*indar[ix,iy] + NTRIA1
            neighs[2*indar[ix,iy]-2+NTRIA1,0] = 3
            neighr[2*indar[ix,iy]-2+NTRIA1,0] = 0

            neighb[2*indar[ix,iy]-1+NTRIA1,2] = 2*indar[ix,iy] -1 + NTRIA1
            neighs[2*indar[ix,iy]-1+NTRIA1,2] = 1
            neighr[2*indar[ix,iy]-1+NTRIA1,2] = 0
            # Top side
            neighb[2*indar[ix,iy]-2+NTRIA1,1] = (2*indar[ix,iy+(iy+1 != com.ny)] + NTRIA1)*(iy+1 != com.ny)
            neighs[2*indar[ix,iy]-2+NTRIA1,1] = 1*(iy+1 != com.ny)
            neighr[2*indar[ix,iy]-2+NTRIA1,1] = 2*(iy+1 == com.ny)
            # Left side
            neighb[2*indar[ix,iy]-2+NTRIA1,2] = (2*indar[ix-1,iy] + NTRIA1)*(ix != 0)
            neighs[2*indar[ix,iy]-2+NTRIA1,2] = 2*(ix != 0)
            neighr[2*indar[ix,iy]-2+NTRIA1,2] = 3*(ix==0)
            # bottom side
            neighb[2*indar[ix,iy]-1+NTRIA1,0] = (2*indar[ix,iy-1] - 1 + NTRIA1)*(iy != 0)
            neighs[2*indar[ix,iy]-1+NTRIA1,0] = 2*(iy != 0)
            neighr[2*indar[ix,iy]-1+NTRIA1,0] = 1*(iy == 0)
            # right side
            neighb[2*indar[ix,iy]-1+NTRIA1,1] = (2*indar[ix+(ix+1 != com.nx),iy] - 1 + NTRIA1)*(ix+1 != com.nx)
            neighs[2*indar[ix,iy]-1+NTRIA1,1] = 3*(ix+1 != com.nx)
            neighr[2*indar[ix,iy]-1+NTRIA1,1] = 4*(ix+1 == com.nx)


    with open(fname,'w') as f:
        f.write('{:>12}\n'.format(N))
    
        k = 0
        for j in range(com.ny):
            for i in range(com.nx):
                for l in range(2):
                    buff = [int(x) for x in [k+1,    
                                    neighb[k,0], neighs[k,0], neighr[k,0],
                                    neighb[k,1], neighs[k,1], neighr[k,1],
                                    neighb[k,2], neighs[k,2], neighr[k,2],
                                    i+1, j+1]]
                    f.write('{:>6d}  {:>6d}{:>6d}{:>6d}    {:>6d}{:>6d}{:>6d}    {:>6d}{:>6d}{:>6d}    {:>6d}{:>6d}\n'.format(*buff))
                    k += 1

    print(" *** Triageom neighbor file for slab case  written for EIRENE ***")





def write_elemente(fname='triageom.elemente'):
    from uedge import com
    from numpy import reshape

    indarr = reshape(list(range(1, (com.nx+1)*(com.ny+1)+1)), (com.nx+1, com.ny+1), 'F')

    with open(fname,'w') as f:
        f.write('{:>12}\n'.format(com.nx*com.ny*2))
    
        k=1
        for j in range(com.ny):
            for i in range(com.nx):
                f.write('{:>6d}  {:>6d}{:>6d}{:>6d}\n'.format(k, indarr[i,j],indarr[i+1,j+1],indarr[i,j+1]))
                f.write('{:>6d}  {:>6d}{:>6d}{:>6d}\n'.format(k+1, indarr[i,j],indarr[i+1,j],indarr[i+1,j+1]))
                k += 2

    print(" *** Triageom elemente file for slab case  written for EIRENE ***")

def write_nodes(fname='triageom.nodes',rm=None, zm=None, flip=False, dispr=0, dispz=0, d=6):
    from uedge import com
    from numpy import zeros, reshape, around

    if rm is None:
        rm = com.rm
    if zm is None:
        zm = com.zm

    x, y = getxy(rm, zm, flip, dispr, dispz)

    N = (com.nx+1)*(com.ny+1)

    xn = zeros((com.nx+1, com.ny+1))
    yn = zeros((com.nx+1, com.ny+1))
    
    xn[:-1,:] = x[1:-1,1:,1]
    xn[-1,:] = x[-1,1:,1]
    yn[:-1,:] = y[1:-1,1:,1]
    yn[-1,:] = y[-1,1:,1]

    xn = split(reshape(around(xn*100,6), (N,), 'F'), 4)
    yn = split(reshape(around(yn*100,6), (N,), 'F'), 4)
    



    with open(fname,'w') as f:
        f.write('{:>12}\n'.format(N))

        for ix in xn:
            f.write(('{{:>{}.{}E}}'.format(19-8+d,d)*len(ix)+'\n').format(*ix).replace('E','{}E'.format('0'*(8-d))))

        for iy in yn:
            f.write(('{{:>{}.{}E}}'.format(19-8+d,d)*len(iy)+'\n').format(*iy).replace('E','{}E'.format('0'*(8-d))))
    print(" *** Triageom node file for slab case  written for EIRENE ***")


def writeEIRENE(sonnet=True, dispr=0, dispz=0, flip=False):

    writeSonnet_slab(fname='fort.30', flip=flip, dispr=dispr, dispz=dispz)
    write31()
    write_nodes(flip=flip, dispr=dispr, dispz=dispz)
    write_elemente()
    write_neighbor() 
