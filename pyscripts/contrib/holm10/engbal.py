

def volavenv(ps):
    from uedge import bbb, com
    from numpy import zeros_like, log
    ps_ave = zeros_like(ps)
    fs0 = 1-4*bbb.fsprd
    for iy in range(1, com.ny+1):
        for ix in range(1, com.nx+1):
            ixm = bbb.ixm1[ix,iy]
            ixp = bbb.ixp1[ix,iy]
            iym = max(0, com.iy-1)
            iyp1 = min(com.ny+1, iy+1)
            if (ps[ix, iy]* ps[ixm, iy]* ps[ixp, iy]* \
                ps[ix, iym]* ps[ix, iyp] > 1e-200):
                signps = 1
                if (ps[ix, iy] < 0):
                    signps = -1
                ps_ave[ix, iy] = sign*exp(
                                fs0 * log(abs(ps[ix, iy])) +
                                fsprd * log(abs(ps[ixm,iy])) +
                                fsprd * log(abs(ps[ixp,iy])) +
                                fsprd * log(abs(ps[ix,iym])) +
                                fsprd * log(abs(ps[ix,iyp])) )
            
    return ps_ave


def engbal():
    from uedge import bbb, com, aph
    from numpy import zeros_like, cos, sum
    if bbb.ishosor == 1:
        volavenv(bbb.prad)
        for igsp in range(com.nhgsp+1, com.ngsp+1):
            jz = igsp - com.nhgsp
            for iimp in range(0, com.nzsp[jz]+1):
                volavenv(pradz[:,:,iimp,jz])

    fekinx = zeros_like(bbb.visy)
    fekiny = zeros_like(bbb.visy)
    fevisx = zeros_like(bbb.visy)
    fevisy = zeros_like(bbb.visy)
    for iy in range(0, com.ny+1):
        for ix in range(0, com.nx+1):
            ixm = bbb.ixm1[ix,iy]
            ixp = bbb.ixp1[ix,iy]
            for iid in range(0,com.nusp):
                thetaix = 0.5*(com.angfx[ixm, iy] + com.angfx[ix, iy])
                thetaixp = 0.5*(com.angfx[ix, iy] + com.angfx[ixp, iy])
                eta_dup2dy = 0.25*( bbb.visy[ix, iy+1, iid] * \
                    (bbb.upi[ixm, iy+1, iid] + bbb.upi[ix, iy+1, iid])**2 - \
                    bbb.visy[ix, iy, iid] * (bbb.upi[ixm, iy, iid] \
                    + bbb.upi[ix, iy, iid])**2)
                fekiny[ix, iy, iid] = (bbb.mi[iid]/32) * (bbb.upi[ixm, iy, iid] \
                    + bbb.upi[ix, iy, iid] + bbb.upi[ixm, iy+1, iid] \
                    + bbb.upi[ix, iy+1, iid] )**2*bbb.fniy[ix, iy, iid]
                fevisy[ix, iy, iid] = -bbb.cfvisy * 0.5 * com.sy[ix, iy] \
                    * com.gyf[ix, iy] * eta_dup2dy 
                fekinx[ix, iy, iid] = 0.5 * bbb.fnix[ix, iy, iid] * bbb.mi[iid]\
                    * bbb.upi[ix, iy, iid]**2 - bbb.upi[ix, iy, iid] \
                    * bbb.fmixy[ix, iy, iid]
                fevisx[ix, iy, iid] = -bbb.cfvisx * 0.25 * com.sx[ix, iy] * (\
                    bbb.visx[ix, iy, iid] * com.gx[ix, iy] * cos(thetaix) *\
                    (bbb.upi[ix, iy, iid]**2 - bbb.upi[ixm, iy, iid]**2) + \
                    bbb.visx[ixp, iy, iid] * com.gx[ixp, iy] * cos(thetaixp) *\
                    (bbb.upi[ixp, iy, iid]**2 - bbb.upi[ix, iy, iid]**2) )

    pmloss = zeros_like(bbb.ne)
    pmpot = zeros_like(bbb.ne)
    pmrada = zeros_like(bbb.ne)
    pmradm = zeros_like(bbb.ne)
    peirad = zeros_like(bbb.ne)
    pvmomcx = zeros_like(bbb.ne)
    pioniz = zeros_like(bbb.ne)
    precom = zeros_like(bbb.ne)
    pioniztot = zeros_like(bbb.ne)
    precomtot = zeros_like(bbb.ne)
    for ix in range(1, com.nx+1):
        for iy in range(1, com.ny+1):
            ixm = bbb.ixm1[ix, iy]
            ixp = bbb.ixp1[ix, iy]
            if bbb.ishymol == 1:
                pmloss[ix, iy] = (1 - bbb.ismolcrm) * bbb.cnsor * (bbb.ediss * \
                    bbb.ev * 0.5 * bbb.psordis[ix, iy, 1] + bbb.ceisor * \
                    bbb.eion * bbb.ev * bbb.psordis[ix, iy, 1] ) + bbb.ismolcrm * \
                    bbb.cnsor * (bbb.cmesori * bbb.emolia[ix, iy] + \
                    bbb.cmesore * bbb.edisse[ix, iy])
                pmpot[ix, iy] = bbb.ismolcrm * bbb.ng[ix, iy, 1] * com.vol[ix, iy] * \
                    aph.sv_crumpet(bbb.te[ix, iy], bbb.ne[ix, iy], 22)
                pmrada[ix, iy]  = bbb.ismolcrm * bbb.ng[ix, iy, 1] * com.vol[ix, iy] * \
                    aph.sv_crumpet(bbb.te[ix, iy], bbb.ne[ix, iy], 23)
                pmradm[ix, iy] = bbb.ismolcrm * bbb.ng[ix, iy, 1] * com.vol[ix, iy] * \
                    aph.sv_crumpet(bbb.te[ix, iy], bbb.ne[ix, iy], 24)

            pioniz[ix, iy] = bbb.cnsor * bbb.ebind * bbb.ev * bbb.psor[ix,iy,0]
            precom[ix, iy] = -bbb.cnsor * bbb.ebind * bbb.ev * bbb.psorrg[ix,iy, 0]
            pioniztot[ix, iy] = pioniz[ix,iy] + bbb.cnsor * bbb.erliz[ix, iy]
            precomtot[ix, iy] = precom[ix,iy] + bbb.cnsor * bbb.erlrc[ix, iy]
            peirad[ix, iy] = bbb.cnsor * bbb.pmloss[ix, iy] + pioniztot[ix, iy] + \
                precomtot[ix,iy]

#            peirad[ix, iy] = bbb.cnsor * (bbb.erliz[ix, iy] + bbb.erlrc[ix, iy] +\
 #               bbb.ebind * bbb.ev * bbb.psor[ix, iy, 0] - bbb.ebind * bbb.ev * \
  #              bbb.psorrg[ix, iy, 0] + pmloss[ix, iy])

            if bbb.isupgon[0] == 0:
                pvmomcx[ix, iy] =  bbb.cngmom[0] * bbb.up[ix, iy, 0] *\
                    com.sx[ix, iy] * com.rrv[ix, iy] * (bbb.ng[ixp, iy, 0] * \
                    bbb.tg[ixp, iy, 0] - bbb.ng[ix, iy, 0] * bbb.tg[ix, iy, 0])\
                    + bbb.cmwall[0] * 0.125*bbb.mi[0] * (bbb.up[ix, iy, 0] + \
                    bbb.up[ixm, iy, 0])**2 * bbb.ng[ix, iy, 0] * \
                    bbb.nucx[ix, iy, 0] * com.vol[ix, iy]

    if (abs(bbb.ckinfl-1) > 1e-10):
        for ix in range(1, com.nxpt):
            ixt = com.ixlb[ix]
            ixtp = ixt + 1 
            ixr = com.ixrb[ix]
            ixrm = ixr - 1
            for iy in range(0, com.ny+1):
                for iid in range(1, com.nfsp+1):
                    fevisx[ixt, iy, iid] = -bbb.ckinfl * 0.5 * com.sx[ixt, iy]\
                        * bbb.visx[ixtp, iy, iid] * com.gx[ixtp, iy] * \
                        (bbb.up[ixtp, iy, iid]**2 - bbb.up[ixt, iy, iid] **2)
                    fekinx[ixt, iy, iid] = 0.5 * bbb.mi[iid] * \
                        bbb.up[ixt, iy, iid]**2 * bbb.fnix[ixt, iy, iid]
                    fevisx[ixr, iy, iid] = -bbb.ckinfl * 0.5 * com.sx[ixr, iy]\
                        * bbb.visx[ixr, iy, iid] * com.gx[ixr, iy] *\
                        (bbb.up[ixr, iy, iid]**2 - up[ixrm, iy, iid]**2)
                    fekinx[ixr, iy, iid] = 0.5 * bbb.mi[iid] * \
                        bbb.up[ixr, iy, iid]**2 * bbb.fnix[ixr, iy, iid] 

    fetx = sum(fekinx, axis=2) + sum(fevisx, axis=2) + bbb.feix + bbb.feex
    fety = sum(fekiny, axis=2) + sum(fevisy, axis=2) + bbb.feiy + bbb.feey
    ptjdote = sum(bbb.wjdote)

    engerr = zeros_like(bbb.ne)
    divfekin = zeros_like(bbb.ne)
    divfevis = zeros_like(bbb.ne)
    divfe = zeros_like(bbb.ne)
    divfetot = zeros_like(bbb.ne)
    volloss = zeros_like(bbb.ne)
    pwrin = bbb.pcoree + bbb.pcorei # TODO: Expression to catch alternative 
                                    #       ways to define Pin!
    for ix in range(1, com.nx+1):
        for iy in range(1, com.ny+1):
            ixm = bbb.ixm1[ix, iy]
            divfekin[ix, iy] = sum(fekinx[ixm, iy] - fekinx[ix, iy] +\
                fekiny[ix, iy-1] - fekiny[ix, iy]) 
            divfevis[ix, iy] = sum(fevisx[ixm, iy] - fevisx[ix, iy] +\
                fevisy[ix, iy-1] - fevisy[ix, iy])
            divfe[ix, iy] = bbb.feex[ixm, iy] - bbb.feex[ix, iy] +\
                bbb.feey[ix, iy-1] - bbb.feey[ix, iy] + \
                sum(bbb.feix[ixm, iy] - bbb.feix[ix, iy] +\
                bbb.feiy[ix, iy-1] - bbb.feiy[ix, iy])
            divfetot[ix, iy] = divfekin[ix,iy] + divfevis[ix, iy] +\
                divfe[ix, iy] 
            volloss[ix, iy] = - peirad[ix, iy] - bbb.png2ni[ix, iy]
            if bbb.isimpon != 0:
                volloss[ix, iy] -= bbb.prad[ix, iy]*com.vol[ix, iy]

            engerr[ix, iy] = (divfetot[ix, iy] + volloss[ix, iy] ) / abs(pwrin)


    # TODO: add ionization and background sources

    return engerr, divfekin, divfevis, divfe, divfetot, pioniz, precom, pmloss
