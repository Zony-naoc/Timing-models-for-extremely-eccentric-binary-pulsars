#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
//  Copyright (C) 2006,2007,2008,2009, George Hobbs, Russell Edwards

/*
 *    This file is part of TEMPO2. 
 * 
 *    TEMPO2 is free software: you can redistribute it and/or modify 
 *    it under the terms of the GNU General Public License as published by 
 *    the Free Software Foundation, either version 3 of the License, or 
 *    (at your option) any later version. 
 *    TEMPO2 is distributed in the hope that it will be useful, 
 *    but WITHOUT ANY WARRANTY; without even the implied warranty of 
 *    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the 
 *    GNU General Public License for more details. 
 *    You should have received a copy of the GNU General Public License 
 *    along with TEMPO2.  If not, see <http://www.gnu.org/licenses/>. 
 */

/*
 *    If you use TEMPO2 then please acknowledge it by citing 
 *    Hobbs, Edwards & Manchester (2006) MNRAS, Vol 369, Issue 2, 
 *    pp. 655-672 (bibtex: 2006MNRAS.369..655H)
 *    or Edwards, Hobbs & Manchester (2006) MNRAS, VOl 372, Issue 4,
 *    pp. 1549-1574 (bibtex: 2006MNRAS.372.1549E) when discussing the
 *    timing model.
 */

#include <stdio.h>
#include <math.h>
#include "tempo2.h"
#include <stdlib.h>

/* Timing model      */
/* Based on bnrydd.f */

longdouble DDEmodel(pulsar *psr,int p,int ipos,int param) //DDE+
{
    longdouble an;
    longdouble pb,k;
    longdouble rad2deg = 180.0/M_PI;
    longdouble SUNMASS = 4.925490947e-6;
    longdouble m2,tt0,t0,x,ecc,er,xdot,edot,dr,dth,eth,am2,ct;
    longdouble pbdot,xpbdot,phase,u,du,gamma;
    longdouble orbits;
    int norbits;
    longdouble  cu,onemecu,cae,sae,ae,omega,omz,sw,cw,alpha,beta,bg,dre,drep,drepp,anhat,su;
    longdouble sqr1me2,cume,brace,si,dlogbr,ds,da,a0,b0,d2bar,torb;
    longdouble csigma,ce,cx,comega,cgamma,cm2,csi,cdth;
    
    /* DDE */
    longdouble zeta1,zeta2,zeta3;
    longdouble czeta1,czeta2,czeta3;
    longdouble xz,eccz;

    /* DDGR+ */
    longdouble xomdot,xk;
    

    const char *CVS_verNum = "$Id: b2338f13ba913de4009f1b09eaa6c3df4bbb5a38 $";

    if (displayCVSversion == 1) CVSdisplayVersion("DDmodel.C","DDmodel()",CVS_verNum);
    

    dr = 0.0; /* WHAT SHOULD THESE BE SET TO? */

    if (psr[p].param[param_dtheta].paramSet[0]==1) dth = getParameterValue(&psr[p],param_dtheta,0);
    else dth = 0.0; 

    if (psr[p].param[param_sini].paramSet[0]==1) si = getParameterValue(&psr[p],param_sini,0);
    else si = 0.0;

    if (si > 1.0)
    {
        displayMsg(1,"BIN1","SIN I > 1.0, setting to 1: should probably use DDS model","",psr[p].noWarnings);
        si = 1.0;
        psr[p].param[param_sini].val[0] = longdouble(1.0);
    }

    if (psr[p].param[param_m2].paramSet[0]==1) am2 = psr[p].param[param_m2].val[0];
    else am2 = 0.0;

    pb = psr[p].param[param_pb].val[0]*SECDAY;
    an = 2.0*M_PI/pb;
    k = psr[p].param[param_omdot].val[0]/(rad2deg*365.25*86400.0*an);
    xk = psr[p].param[param_xomdot].val[0]/(rad2deg*365.25*86400.0*an);

    m2 = am2*SUNMASS;
    t0 = psr[p].param[param_t0].val[0];
    ct = psr[p].obsn[ipos].bbat;    

    tt0 = (ct-t0)*SECDAY;

    if (psr[p].param[param_gamma].paramSet[0]==1)
        gamma = psr[p].param[param_gamma].val[0];
    else
        gamma = 0.0;
    if (psr[p].param[param_a0].paramSet[0]==1) a0 = 1e-6*psr[p].param[param_a0].val[0];
    else 
        a0 = 0.0; /* WHAT SHOULD THIS BE SET TO? */
    b0    = 0.0; /* WHAT SHOULD THIS BE SET TO? */

    // if (psr[p].param[param_om].paramSet[0]==1) omz = psr[p].param[param_om].val[0];
    // else omz = 0.0;

    if (psr[p].param[param_a1dot].paramSet[0]==1) xdot  = psr[p].param[param_a1dot].val[0];
    else xdot  = 0.0;

    if (psr[p].param[param_pbdot].paramSet[0] == 1) pbdot = psr[p].param[param_pbdot].val[0];
    else pbdot = 0.0;

    if (psr[p].param[param_edot].paramSet[0] == 1) edot = psr[p].param[param_edot].val[0];
    else edot = 0.0;

    if (psr[p].param[param_xpbdot].paramSet[0] == 1) xpbdot = psr[p].param[param_xpbdot].val[0];
    else xpbdot = 0.0;

    if (psr[p].param[param_zeta1].paramSet[0] == 1) zeta1 = psr[p].param[param_zeta1].val[0];
    else{
        printf("zeta1 no found!!");
        exit(0);
    }

    if (psr[p].param[param_zeta2].paramSet[0] == 1) zeta2 = psr[p].param[param_zeta2].val[0];
    else{
        printf("zeta2 no found!!");
        exit(0);
    }

    if (psr[p].param[param_zeta3].paramSet[0] == 1) zeta3 = psr[p].param[param_zeta3].val[0];
    else{
        printf("zeta3 no found!!");
        exit(0);
    }
    
    //zeta1=x*sw*(2-ecc);zeta2=x*cw*sqrt(1-ecc*ecc);zeta3=x*sw*(2*ecc-1)
    
    eccz = sqrt(1-zeta3*zeta3); 
    xz   = sqrt(powl(zeta1*eccz,2.0)+powl(zeta2/zeta3,2.0));
	omz  = atan2(zeta1*eccz,zeta2/zeta3);
    
    x = xz+xdot*tt0;
    ecc = eccz+edot*tt0;

    // x = psr[p].param[param_a1].val[0]+xdot*tt0;
    // ecc = psr[p].param[param_ecc].val[0]+edot*tt0;
    er = ecc*(1.0+dr);
    eth = ecc*(1.0+dth);

    if (ecc < 0.0 || ecc > 1.0)
    {
        ld_printf("DDmodel: problem with eccentricity = %Lg [%s]\n",psr[p].param[param_ecc].val[0],psr[p].name);
        exit(1);
    }

    orbits = tt0/pb - longdouble(0.5)*(pbdot+xpbdot)*(tt0/pb)*(tt0/pb);
    norbits = (int)orbits;
    if (orbits<0.0) norbits--;
    phase=2.0*M_PI*(orbits-norbits);
    /*  Compute eccentric anomaly u by iterating Kepler's equation. */
    u=M_PI;
    do {
        du=(phase-(u-ecc*sin(u)))/(1.0-ecc*cos(u));
        u=u+du;
    } while (fabs(du)>1.0e-12);

    /*  DD equations 17b, 17c, 29, and 46 through 52 */
    su=sin(u);
    cu=cos(u);
    onemecu=1.0-ecc*cu;
    cae=(cu-ecc)/onemecu;
    sae=sqrt(1.0-pow(ecc,2))*su/onemecu;
    ae=atan2(sae,cae);
    if(ae<0.0) ae=ae+2.0*M_PI;
    ae=2.0*M_PI*orbits + ae - phase;

    if (psr[p].param[param_xomdot].paramSet[0]==1) omega=omz + k*(ae-2.0*M_PI*orbits) + xk*2.0*M_PI*orbits;
    else omega=omz + k*ae;
    sw=sin(omega);
    cw=cos(omega);
    alpha=x*sw;
    beta=x*sqrt(1-pow(eth,2))*cw;
    bg=beta+gamma;
    dre=alpha*(cu-er) + bg*su;
    drep=-alpha*su + bg*cu;
    drepp=-alpha*cu - bg*su;
    anhat=an/onemecu;

    /* DD equations 26, 27, 57: */
    sqr1me2=sqrt(1-pow(ecc,2));
    cume=cu-ecc;
    brace=onemecu-si*(sw*cume+sqr1me2*cw*su);
    //  printf("GEORGE: si = %g, brace = %g\n",(double)si,(double)brace);
    dlogbr=log(brace);
    ds=-2*m2*dlogbr;
    da=a0*(sin(omega+ae) + ecc*sw) + b0*(cos(omega+ae) + ecc*cw);

    /*  Now compute d2bar, the orbital time correction in DD equation 42. */
    d2bar=dre*(1-anhat*drep+(pow(anhat,2))*(pow(drep,2) + 0.5*dre*drepp -
                0.5*ecc*su*dre*drep/onemecu)) + ds + da + 1.5*zeta1*(1-zeta3*zeta3);
    torb=-d2bar;

    if (param==-1) return torb;

    /*  Now we need the partial derivatives. Use DD equations 62a - 62k. */
    csigma=x*(-sw*su+sqr1me2*cw*cu)/onemecu;
    ce=su*csigma-x*sw-ecc*x*cw*su/sqr1me2;
    cx=sw*cume+sqr1me2*cw*su;
    comega=x*(cw*cume-sqr1me2*sw*su);
    cgamma=su;
    cdth=-ecc*ecc*x*cw*su/sqr1me2;
    cm2=-2*dlogbr;
    csi=2*m2*(sw*cume+sqr1me2*cw*su)/brace; 

    
    longdouble de2dzeta3;
    de2dzeta3=-zeta3/eccz;

    longdouble dx2dzeta1,dx2dzeta2,dx2dzeta3;
    dx2dzeta1=zeta1*eccz*eccz/xz;
    dx2dzeta2=zeta2/(xz*zeta3*zeta3);
    dx2dzeta3=-(zeta1*zeta1*zeta3+zeta2*zeta2/(zeta3*zeta3*zeta3))/xz;

    longdouble dw2dzeta1,dw2dzeta2,dw2dzeta3;
    dw2dzeta1=zeta2*eccz/(xz*xz*zeta3);
    dw2dzeta2=-eccz*zeta1/(xz*xz*zeta3);
    dw2dzeta3=zeta2/(xz*xz*zeta3*zeta3)*(eccz-zeta3*zeta3/eccz)*zeta1;

    czeta1=dx2dzeta1*cx+dw2dzeta1*comega;
    czeta2=dx2dzeta2*cx+dw2dzeta2*comega;
    czeta3=dx2dzeta3*cx+dw2dzeta3*comega+de2dzeta3*ce;

    if (param==param_pb)
        return -csigma*an*SECDAY*tt0/(pb*SECDAY); 
    else if (param==param_zeta1)
        return czeta1;
    else if (param==param_zeta2)
        return czeta2;
    else if (param==param_zeta3)
        return czeta3;
    else if (param==param_edot)
        return ce*tt0;
    else if (param==param_omdot){
        if (psr[p].param[param_xomdot].paramSet[0]==1) return (ae-2.0*M_PI*orbits)*comega/(an*360.0/(2.0*M_PI)*365.25*SECDAY);
        else return ae*comega/(an*360.0/(2.0*M_PI)*365.25*SECDAY);
    }
    else if (param==param_xomdot){       
        return 2.0*M_PI*orbits*comega/(an*360.0/(2.0*M_PI)*365.25*SECDAY);
    }
    else if (param==param_t0)
        return -csigma*an*SECDAY;
    else if (param==param_pbdot)
        return 0.5*tt0*(-csigma*an*SECDAY*tt0/(pb*SECDAY));
    else if (param==param_sini)
        return csi;
    else if (param==param_gamma)
        return cgamma;
    else if (param==param_dtheta)
        return cdth;
    else if (param==param_m2)
        return cm2*SUNMASS;
    else if (param==param_a1dot) /* Also known as xdot */
        return cx*tt0;

    return 0.0;
}


void updateDDE(pulsar *psr,double val,double err,int pos)
{
    if (pos==param_pb)
    {
        psr->param[param_pb].val[0] += val/SECDAY;
        psr->param[param_pb].err[0]  = err/SECDAY;
    }
    else if (pos==param_zeta1 || pos==param_zeta2 || pos==param_zeta3 || pos==param_t0 || pos==param_sini || pos==param_m2
            || pos == param_gamma || pos==param_edot)
    {
        psr->param[pos].val[0] += val;
        psr->param[pos].err[0]  = err;
    }
    // else if (pos==param_om)
    // {
    //     psr->param[pos].val[0] += val*180.0/M_PI;
    //     psr->param[pos].err[0]  = err*180.0/M_PI;
    // }
    else if (pos==param_pbdot)
    {
        psr->param[pos].val[0] += val;
        psr->param[pos].err[0]  = err;
    }
    else if (pos==param_a1dot)
    {
        psr->param[pos].val[0] += val;
        psr->param[pos].err[0]  = err;
    }
    else if (pos==param_dtheta)
    {
        psr->param[pos].val[0] += val;
        psr->param[pos].err[0]  = err;
    }
    else if (pos==param_omdot||pos==param_xomdot)
    {
        psr->param[pos].val[0] += val; /* *(SECDAY*365.25)*180.0/M_PI; */
        psr->param[pos].err[0]  = err; /* *(SECDAY*365.25)*180.0/M_PI; */
    }
}
