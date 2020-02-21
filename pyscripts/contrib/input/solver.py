# Solver package for UEDGE, solver.py
# Separated from 160299 master input by holm10
# Changelog
# 200213 - file created
from uedge import bbb,com,grd,flx

def nksol(   ftol=1e-8, itermx=30, mfnksol=3, mdif=0, incpset=5, premeth='ilut',lfililut=200,lenpfac=200,lenplufac=200,
                    xrinc=2, xlinc=3, yinc=2, scrit=1e-3, nurlx=1e7, rlx=0.9, iscolnorm=3
                 ):

    """=====================================================================================================
    NONLINEAR KRYLOV SOLVER SETUP
    ====================================================================================================="""

    # Solver used
    #- - - - - - - 
    bbb.svrpkg=="nksol"
    bbb.ftol=ftol         # Stop tolerance
    bbb.itermx=itermx           # Maximum number of "nksol" iterations
    bbb.mfnksol=mfnksol           #nksol method flag; 
                                #=1 means dogleg strategy
                                #=2 means linesearch with Arnoldi method,
                                #=3 means linesearch with GMRES method.
                                #negative mfnksol ignores global constaints
    bbb.mdif=mdif              # =0 Matrix vector multiplication J*v approximated by numerical finite difference RHS (def)
                            # =1 Matrix vector multiplication J*v is done direcctly with current J
    bbb.rlx=rlx             # Restricts relative change of temperature and density to rlx: increasing loosens convergence criteria
    bbb.incpset=incpset           # Nonlinear iterations before Jacobian is re-evaluated
    # Decomposition
    bbb.premeth=premeth      # Preconditioning routine for linear iterations:
                                #"banded": use full banded jacobian as preconditiones
                                #"ilut": use ilut preconditioning
                                #"inel": use INEL ILU preconditioning
    bbb.iscolnorm=iscolnorm         # Column normalization setup for y's
                                # =0 for no implicit scaling (suscal=1)
                                # =1 for scaling by normalization constants
                                # =2 for scaling by max(abs(yl),floors)
                                # =3 combination of global scaling with nnorm,
                                    #    etc, followed by local scaling by each yl
    if bbb.premeth=='ilut':
        bbb.lfililut=lfililut        # Fill-in parameter used in ilut: will allow up to lfililut additional nonzeros in each L/U row
    bbb.lenpfac=lenpfac         #  fudge factor to multiply neq by to get an estimate for the number of nonzeros in the preconditioner matrix.
    bbb.lenplufac=lenplufac       #  fudge factor to multiply neq by to get an estimate for the number of nonzeros in the factored preconditioner matrix.
        

    # Set common parameters
    common(xrinc, xlinc, yinc, scrit, nurlx)


    # Differencing methods
    differencing()







def vodpk(   runtim=1e-7, delpy=1e-8, jacflg=1,jpre=1,rtolv=1e-4,inopt=0,cngatol=1,
                    xrinc=2, xlinc=3, yinc=2, scrit=1e-3, nurlx=1e7
                ):

    """=====================================================================================================
    PRECONDITIONED KRYLOV SOLVER SETUP FOR ODES
    ====================================================================================================="""

    # Solver used
    #- - - - - - - 
    bbb.svrpkg=="vodpk"
    bbb.runtim=runtim          # convergence for VODPK (Time for first output)
                                    #=1, set uniform, fixed density, ncore
    bbb.delpy=delpy             # Sets del: fractional change for finite differences derivative for "vodpk"
    bbb.jacflg=facflg                # Flag for computing Jacobian in "vodpk"
    bbb.jpre=jpre                  # Flag for using preconditioning in "vodpk"
    bbb.rtolv=rtolv             # Relative tolerance vector usid in convert.m
    bbb.inopt=inopt                 # Resets options for vodpk and daspk solvers
    bbb.cngatol=cngatol  # Coefficient for absolute tolerance for LHS-ode - like routines

    
    # Set common parameters
    common(xrinc, xlinc, yinc, scrit, nurlx)

    # Differencing methods
    differencing()


   

def newton(rwmin=1e-11, nmaxnewt=20, xrinc=2, xlinc=3, yinc=2, scrit=1e-3, nurlx=1e7):

    """=====================================================================================================
    DIRECT NEWTON SOLVER SETUP FOR ODES
    ====================================================================================================="""

    # Solver used
    #- - - - - - - 
    bbb.svrpkg=="newton"
    bbb.rwmin=rwmin           # convergence for newton solver
    bbb.nmaxnewt=nmaxnewt             # Maximum number of Newton iterations

    # Set common parameters
    common(xrinc, xlinc, yinc, scrit, nurlx)

    # Differencing methods
    #- - - - - - - - - - - -
    differencing()


def differencing(methn=33,methu=83,methe=33,methi=33,methg=66,methp=33):
    ''' Sets up differencing methods
        differencing(*keys)
        
        Optional parameters
        methn (33)  -   Ion continuity equation differencing
        methu (83)  -   Ion momentum equation differencing
        methe (33)  -   Electron energy equation differencing
        methi (33)  -   Ion energy equation differencing
        methg (66)  -   Neutral gas continuity differencing
        methp (33)  -   Potential equation differencing
        


        Differencing method flags for schemes: x-dir scheme=mod(meth,10), y-dir=floor(meth/10)
        Differencing methods: (Oderhs L100)
        0: Central differencing, assuming no flows
        1: Upwind differencing, assuming dif=0
        2: Central differencing (harmonic average)
        3: Upwind differencing
        4: Simple hybrid scheme
        5: Fifth power scheme
        6: Logarithmic central differencing
        7: Regular upwind differencing: can be used fir methg=77 in nonorthogonall differencing
        8: Regular upwind differencing on an approximate staggered stencil (for velocities)
    '''
    bbb.methn=methn    # Ion continuity equation differencing
    bbb.methu=methu    # Ion momentum equation differencing
    bbb.methe=methe    # Electron energy equation differencing
    bbb.methi=methi    # Ion energy equation differencing
    bbb.methg=methg    # Neutral gas continuity differencing
    bbb.methp=methp    # Potential equation differencing


    

def common(xrinc, xlinc, yinc, scrit, nurlx):
 
    # JACOBIAN INCREMENTS
    bbb.xrinc=xrinc             # Jacobian parameter in x-dir, right bound: cell indices
    bbb.xlinc=xlinc             # Jacobian parameter in x-dir, left bound: cell indices
    bbb.yinc=yinc              # Jacobian parameter in y-dir: cell indices
    bbb.scrit=scrit         # Upper limit for re-computing Jacobian and LU decomposition (when sum of residuals>scrit)

    # SOLVER OPTIONS AND COEFFICIENTS
    bbb.nurlx=nurlx  # Rate coefficient to relax boundary condition





