from matplotlib.pyplot import ion
from uedge import bbb, com
ion()

def natsort(l): 
    from re import split
    convert = lambda text: int(text) if text.isdigit() else text.lower()
    alphanum_key = lambda key: [convert(c) for c in split('([0-9]+)', key)]
    return sorted(l, key=alphanum_key)

def nesepfrac():
    # Accounting for cell sizes
    sx = com.sx[bbb.ixmp, com.iysptrx:com.iysptrx+2]
    wt = [x/sum(sx) for x in sx]    

    return sum(wt * bbb.ne[bbb.ixmp, com.iysptrx:com.iysptrx+2])/\
        bbb.ne[bbb.ixmp,1]

    

def teomp():
    # Accounting for cell sizes
    sx = com.sx[bbb.ixmp, com.iysptrx:com.iysptrx+2]
    wt = [x/sum(sx) for x in sx]    

    return sum(wt * bbb.te[bbb.ixmp, com.iysptrx:com.iysptrx+2])/bbb.ev


def plot_nfit(ix=None, ishift=1):

    if ix is None:
        ix = com.ixpt2[0]

    n = bbb.ne
    return plot_fit(n, ix, ishift)


def plot_qfit(ix=None):

    if ix is None:
        ix = com.ixpt2[0]

    qpara =  bbb.feex + bbb.feix
    qpara /= abs(com.bpol/com.b)[:,:,0]
#    qpara /= abs(bbb.rbfbt)
    qpara /= com.sx

    return plot_fit(qpara, ix)

    
def plot_fit(y, ix=None, ishift=0):
    from matplotlib.pyplot import subplots
    from  numpy import log, exp, zeros, array, matmul, linspace, where
    from numpy.linalg import inv
    f, ax = subplots()
    
    if ix is None:
        ix = com.ixpt2[0]

    yfull = y[ix,2:-1]
    xfull = com.yyc[2:-1]
#    y = bbb.feex[ix, com.iysptrx+1:-1] + bbb.feix[ix, com.iysptrx+1:-1]
#    y /= abs(bbb.rbfbt[ix, com.iysptrx+1:-1])

    # Only analyze main SOL
    y = y[ix,com.iysptrx+1+ishift:-1]
    x = com.yyc[com.iysptrx+1+ishift:-1]
    # Start from peak, in case it has drifted outwards
    istart = where(y==y.max())[0][0]
    y = y[istart:]
    x = x[istart:]



    # Algorithm from Regressions et equations integrales, Jean Jaquelin, 2009
    S=[0]
    for i in range(1, len(y)):
        S.append(S[i-1]+0.5*((y[i]+y[i-1])*(x[i]-x[i-1])))

    mat1 = zeros((2,2))
    mat2 = zeros((2,1))
    mat3 = zeros((2,2))
    mat4 = zeros((2,1))
    mat3[0,0] = len(y)
    for i in range(len(y)):
        mat1[0,0] += (x[i]-x[0])**2
        mat1[0,1] += (x[i]-x[0])*S[i]
        mat1[1,0] += (x[i]-x[0])*S[i]
        mat1[1,1] += S[i]**2

        mat2[0] += (y[i]-y[0])*(x[i]-x[0])
        mat2[1] += (y[i]-y[0])*S[i]
    
    [_, c] = matmul(inv(mat1), mat2)
    for i in range(len(y)):
        mat3[0,1] += exp(c*x[i])
        mat3[1,0] += exp(c*x[i])
        mat3[1,1] += exp(2*c*x[i])

        mat4[0] += y[i]
        mat4[1] += y[i]*exp(c*x[i])

    [a, b] = matmul(inv(mat3), mat4)

    ax.plot(xfull, yfull, 'k.')
    
    xl = linspace(x[0],x[-1])
    ax.plot(xl, a+b*exp(c*xl), 'b-')
    ax.axvline(0, linewidth=0.5, color='k')

    return f


def lambdanfrac(ix=None, ishiftn=1):
    return calc_lambdan(ix, ishiftn)/calc_lambdaq(ix)

def calc_lambdan(ix=None, ishift=1):
    
    if ix is None:
        ix = com.ixpt2[0]

    n = bbb.ne

    return calc_lambda(n, ix, ishift)


def calc_lambdaq(ix=None):
    
    if ix is None:
        ix = com.ixpt2[0]

    qpara =  bbb.feex + bbb.feix
    qpara /= abs(com.bpol/com.b)[:,:,0]
#    qpara /= abs(bbb.rbfbt)
    qpara /= com.sx

    return calc_lambda(qpara, ix)

def calc_lambda(y, ix=None, ishift=0):
    from matplotlib.pyplot import subplots
    from  numpy import log, exp, zeros, array, matmul, linspace, where
    from numpy.linalg import inv

    if ix is None:
        ix = com.ixpt2[0]


    yfull = y[ix,2:-1]
    xfull = com.yyc[2:-1]
#    y = bbb.feex[ix, com.iysptrx+1:-1] + bbb.feix[ix, com.iysptrx+1:-1]
#    y /= abs(bbb.rbfbt[ix, com.iysptrx+1:-1])

    # Only analyze main SOL
    y = y[ix,com.iysptrx+1+ishift:-1]
    x = com.yyc[com.iysptrx+1+ishift:-1]
    # Start from peak, in case it has drifted outwards
    istart = where(y==y.max())[0][0]
    y = y[istart:]
    x = x[istart:]



    # Algorithm from Regressions et equations integrales, Jean Jaquelin, 2009
    S=[0]
    for i in range(1, len(y)):
        S.append(S[i-1]+0.5*((y[i]+y[i-1])*(x[i]-x[i-1])))

    mat1 = zeros((2,2))
    mat2 = zeros((2,1))
    for i in range(len(y)):
        mat1[0,0] += (x[i]-x[0])**2
        mat1[0,1] += (x[i]-x[0])*S[i]
        mat1[1,0] += (x[i]-x[0])*S[i]
        mat1[1,1] += S[i]**2

        mat2[0] += (y[i]-y[0])*(x[i]-x[0])
        mat2[1] += (y[i]-y[0])*S[i]

    [_, c] = matmul(inv(mat1), mat2)
    return -1/c[0]


def kappae(lnl=15, zeff=1.6):
    return 2.16*25000/(lnl*(1+0.27*zeff))



def solve_teosp(Psol, lambdaq, ne=None, alpha=0.3, kappa=None):
    from scipy.optimize import minimize
    
    def solve(Te, Psol, lambdaq, ne, alpha, kappa):
        from numpy import pi
        L = sum((((com.b/com.bpol)[:,:,0])/com.gx)[bbb.ixmp:,com.iysptrx+1])
        BpBpol = (com.b/com.bpol)[bbb.ixmp,com.iysptrx+1,0]
        R = com.rm[bbb.ixmp, com.iysptrx+1, 0]
        qpara = (Psol/(8*pi*R*lambdaq))*BpBpol
        if kappa is None:
            kappa = kappae()
        if ne is None:
            ne = bbb.ne[bbb.ixmp, com.iysptrx+1]
        qSp = (2*kappa*Te**(7/2))/(7*L)
        qfl = alpha*ne*1.602e-19*Te*(1.602e-19*Te/9.11e-31)**0.5
        qSpfl = 1/(1/qSp+1/qfl)
        return abs(qpara-qSpfl)
    
    return  minimize(solve, 100, args=(Psol, lambdaq, ne, alpha, kappa),
            method = 'Nelder-Mead')

def teosp_powerbalance(Psol, lambdaq, kappa=None):
    ''' Calculates Te from Psol [W] and lamda_q [m]'''
    from numpy import pi
    L = sum((((com.b/com.bpol)[:,:,0])/com.gx)[bbb.ixmp:,com.iysptrx+1])
    BpBpol = (com.b/com.bpol)[bbb.ixmp,com.iysptrx+1,0]
    R = com.rm[bbb.ixmp, com.iysptrx+1, 0]
    if kappa is None:
        kappa = kappae()
    print('qpara', (Psol/(8*pi*R*lambdaq))*BpBpol)
    print('qSp', 2*kappa/(7*L))
    return ((7*Psol*L*BpBpol)/(16*pi*R*lambdaq*kappa))**(2/7)


def restore():
    from copy import deepcopy
    issfon = deepcopy(bbb.issfon)
    ftol = deepcopy(bbb.ftol)    
    bbb.issfon=0
    bbb.ftol=1e20
    bbb.exmain()
    bbb.issfon = issfon
    bbb.ftol = ftol
    


def analyze_folder(path='../solutions', I=1.3, a=0.7, ncol=5, fontsize=10):
    from os import walk
    from random import shuffle
    from uedge.hdf5 import hdf5_restore
    from itertools import product
    from tqdm import tqdm
    from numpy import pi
    from input import restore_input
    bbb.iprint = 0
    restore_input() # Restores grid etc    
    restore()
    from matplotlib.pyplot import subplots
    nGW = 10**20 * I / pi / a**2
    f, ax = subplots(2,2, figsize=(12,8))

    colors = ['k', 'r', 'b', 'c', 'm', 'y', 'g']
    marker = ['o', 's', 'v', 'P', 'X', 'p', '*','^']

    markers = []
    for m in product(colors, marker):
        markers.append(m[0]+m[1])
    shuffle(markers)

    saves = []
    for (path, dirs, files) in walk(path):
        saves.extend(files)

    saves = natsort(saves)
    casename = []
    te = []
    nefrac = []
    lambdaq = []
    lambdafrac = []
    nenGW = []

    for save in tqdm(saves):
        if 'hdf5' in save.lower():
            hdf5_restore('{}/{}'.format(path,save))
            restore()
            casename.append(save.replace('.hdf5','').replace('_last_ii2',''))
            te.append(teomp())
            nefrac.append(nesepfrac())
            lambdaq.append(calc_lambdaq())
            lambdafrac.append(lambdanfrac())
            nenGW.append(bbb.ne[bbb.ixmp,1]/nGW)

    bbb.iprint = 1

    for i in range(len(casename)):
        ax[0,0].plot(nenGW[i], te[i], markers[i])
        ax[0,0].set_xlabel('ne/nGW')
        ax[0,0].set_ylabel('tesep')

        ax[0,1].plot(nenGW[i], nefrac[i], markers[i])
        ax[0,1].set_xlabel('ne/nGW')
        ax[0,1].set_ylabel('nesep/neped')

        ax[1,0].plot(nenGW[i], lambdaq[i], markers[i])
        ax[1,0].set_xlabel('ne/nGW')
        ax[1,0].set_ylabel('lamdaq')

        ax[1,1].plot(nenGW[i], lambdafrac[i], markers[i])
        ax[1,1].set_xlabel('ne/nGW')
        ax[1,1].set_ylabel('lambdan/lambdaq')
    
    f.subplots_adjust(bottom=0.3)
    ax[1,0].legend(casename, bbox_to_anchor=(2.4, -0.2), ncol=ncol, frameon=False, fontsize=fontsize)

    ax[0,0].axhline(100, color='k', linewidth=0.5)
    ax[0,1].axhline(0.3, color='k', linewidth=0.5)
    ax[0,1].axhline(0.5, color='k', linewidth=0.5)
    ax[1,0].axhline(0.002, color='k', linewidth=0.5)
    ax[1,1].axhline(10, color='k', linewidth=0.5)
    ax[1,1].axhline(3, color='k', linewidth=0.5)

    return f
