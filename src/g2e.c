/*
 * Copyright (C) 2013-  Qiming Sun <osirpt.sun@gmail.com>
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include "config.h"
#include "cint_bas.h"
#include "rys_roots.h"
#include "misc.h"
#include "g2e.h"

#define DEF_GXYZ(type, G, GX, GY, GZ) \
        type *GX = G; \
        type *GY = G + envs->g_size; \
        type *GZ = G + envs->g_size * 2

void CINTinit_int2e_EnvVars(CINTEnvVars *envs, FINT *ng, FINT *shls,
                            FINT *atm, FINT natm, FINT *bas, FINT nbas, double *env)
{
        envs->natm = natm;
        envs->nbas = nbas;
        envs->atm = atm;
        envs->bas = bas;
        envs->env = env;
        envs->shls = shls;

        const FINT i_sh = shls[0];
        const FINT j_sh = shls[1];
        const FINT k_sh = shls[2];
        const FINT l_sh = shls[3];
        envs->i_l = bas(ANG_OF, i_sh);
        envs->j_l = bas(ANG_OF, j_sh);
        envs->k_l = bas(ANG_OF, k_sh);
        envs->l_l = bas(ANG_OF, l_sh);
        envs->x_ctr[0] = bas(NCTR_OF, i_sh);
        envs->x_ctr[1] = bas(NCTR_OF, j_sh);
        envs->x_ctr[2] = bas(NCTR_OF, k_sh);
        envs->x_ctr[3] = bas(NCTR_OF, l_sh);
        envs->nfi = (envs->i_l+1)*(envs->i_l+2)/2;
        envs->nfj = (envs->j_l+1)*(envs->j_l+2)/2;
        envs->nfk = (envs->k_l+1)*(envs->k_l+2)/2;
        envs->nfl = (envs->l_l+1)*(envs->l_l+2)/2;
        envs->nf = envs->nfi * envs->nfk * envs->nfl * envs->nfj;

        envs->ri = env + atm(PTR_COORD, bas(ATOM_OF, i_sh));
        envs->rj = env + atm(PTR_COORD, bas(ATOM_OF, j_sh));
        envs->rk = env + atm(PTR_COORD, bas(ATOM_OF, k_sh));
        envs->rl = env + atm(PTR_COORD, bas(ATOM_OF, l_sh));

        envs->common_factor = (M_PI*M_PI*M_PI)*2/SQRTPI
                * CINTcommon_fac_sp(envs->i_l) * CINTcommon_fac_sp(envs->j_l)
                * CINTcommon_fac_sp(envs->k_l) * CINTcommon_fac_sp(envs->l_l);
        if (env[PTR_EXPCUTOFF] == 0) {
                envs->expcutoff = EXPCUTOFF;
        } else {
                // +1 to ensure accuracy. See comments in function CINT2e_loop_nopt
                envs->expcutoff = MAX(MIN_EXPCUTOFF, env[PTR_EXPCUTOFF]) + 1;
        }

        envs->gbits = ng[GSHIFT];
        envs->ncomp_e1 = ng[POS_E1];
        envs->ncomp_e2 = ng[POS_E2];
        envs->ncomp_tensor = ng[TENSOR];

        envs->li_ceil = envs->i_l + ng[IINC];
        envs->lj_ceil = envs->j_l + ng[JINC];
        envs->lk_ceil = envs->k_l + ng[KINC];
        envs->ll_ceil = envs->l_l + ng[LINC];
        int rys_order =(envs->li_ceil + envs->lj_ceil
                        + envs->lk_ceil + envs->ll_ceil)/2 + 1;
        int nrys_roots = rys_order;
#ifdef WITH_RANGE_COULOMB
        double omega = env[PTR_RANGE_OMEGA];
        if (omega < 0 && rys_order <= 3) {
                nrys_roots *= 2;
        }
#endif
        envs->rys_order = rys_order;
        envs->nrys_roots = nrys_roots;

        assert(i_sh < SHLS_MAX);
        assert(j_sh < SHLS_MAX);
        assert(k_sh < SHLS_MAX);
        assert(l_sh < SHLS_MAX);
        assert(envs->i_l < ANG_MAX);
        assert(envs->j_l < ANG_MAX);
        assert(envs->k_l < ANG_MAX);
        assert(envs->l_l < ANG_MAX);
        assert(bas(ATOM_OF,i_sh) >= 0);
        assert(bas(ATOM_OF,j_sh) >= 0);
        assert(bas(ATOM_OF,k_sh) >= 0);
        assert(bas(ATOM_OF,l_sh) >= 0);
        assert(bas(ATOM_OF,i_sh) < natm);
        assert(bas(ATOM_OF,j_sh) < natm);
        assert(bas(ATOM_OF,k_sh) < natm);
        assert(bas(ATOM_OF,l_sh) < natm);
        assert(rys_order < MXRYSROOTS);

        FINT dli, dlj, dlk, dll;
        FINT ibase = envs->li_ceil > envs->lj_ceil;
        FINT kbase = envs->lk_ceil > envs->ll_ceil;
        if (kbase) {
                dlk = envs->lk_ceil + envs->ll_ceil + 1;
                dll = envs->ll_ceil + 1;
        } else {
                dlk = envs->lk_ceil + 1;
                dll = envs->lk_ceil + envs->ll_ceil + 1;
        }

        if (ibase) {
                dli = envs->li_ceil + envs->lj_ceil + 1;
                dlj = envs->lj_ceil + 1;
        } else {
                dli = envs->li_ceil + 1;
                dlj = envs->li_ceil + envs->lj_ceil + 1;
        }
        envs->g_stride_i = nrys_roots;
        envs->g_stride_k = nrys_roots * dli;
        envs->g_stride_l = nrys_roots * dli * dlk;
        envs->g_stride_j = nrys_roots * dli * dlk * dll;
        envs->g_size     = nrys_roots * dli * dlk * dll * dlj;

        if (kbase) {
                envs->g2d_klmax = envs->g_stride_k;
                envs->rx_in_rklrx = envs->rk;
                envs->rkrl[0] = envs->rk[0] - envs->rl[0];
                envs->rkrl[1] = envs->rk[1] - envs->rl[1];
                envs->rkrl[2] = envs->rk[2] - envs->rl[2];
        } else {
                envs->g2d_klmax = envs->g_stride_l;
                envs->rx_in_rklrx = envs->rl;
                envs->rkrl[0] = envs->rl[0] - envs->rk[0];
                envs->rkrl[1] = envs->rl[1] - envs->rk[1];
                envs->rkrl[2] = envs->rl[2] - envs->rk[2];
        }

        if (ibase) {
                envs->g2d_ijmax = envs->g_stride_i;
                envs->rx_in_rijrx = envs->ri;
                envs->rirj[0] = envs->ri[0] - envs->rj[0];
                envs->rirj[1] = envs->ri[1] - envs->rj[1];
                envs->rirj[2] = envs->ri[2] - envs->rj[2];
        } else {
                envs->g2d_ijmax = envs->g_stride_j;
                envs->rx_in_rijrx = envs->rj;
                envs->rirj[0] = envs->rj[0] - envs->ri[0];
                envs->rirj[1] = envs->rj[1] - envs->ri[1];
                envs->rirj[2] = envs->rj[2] - envs->ri[2];
        }

        if (rys_order <= 2) {
                envs->f_g0_2d4d = &CINTg0_2e_2d4d_unrolled;
#ifdef WITH_RANGE_COULOMB
                if (rys_order != nrys_roots) {
                        envs->f_g0_2d4d = &CINTsrg0_2e_2d4d_unrolled;
                }
#endif
        } else if (kbase) {
                if (ibase) {
                        envs->f_g0_2d4d = &CINTg0_2e_ik2d4d;
                } else {
                        envs->f_g0_2d4d = &CINTg0_2e_kj2d4d;
                }
        } else {
                if (ibase) {
                        envs->f_g0_2d4d = &CINTg0_2e_il2d4d;
                } else {
                        envs->f_g0_2d4d = &CINTg0_2e_lj2d4d;
                }
        }
        envs->f_g0_2e = &CINTg0_2e;
}

void CINTg2e_index_xyz(FINT *idx, const CINTEnvVars *envs)
{
        const FINT i_l = envs->i_l;
        const FINT j_l = envs->j_l;
        const FINT k_l = envs->k_l;
        const FINT l_l = envs->l_l;
        const FINT nfi = envs->nfi;
        const FINT nfj = envs->nfj;
        const FINT nfk = envs->nfk;
        const FINT nfl = envs->nfl;
        const FINT di = envs->g_stride_i;
        const FINT dk = envs->g_stride_k;
        const FINT dl = envs->g_stride_l;
        const FINT dj = envs->g_stride_j;
        FINT i, j, k, l, n;
        FINT ofx, ofkx, oflx;
        FINT ofy, ofky, ofly;
        FINT ofz, ofkz, oflz;
        FINT i_nx[CART_MAX], i_ny[CART_MAX], i_nz[CART_MAX];
        FINT j_nx[CART_MAX], j_ny[CART_MAX], j_nz[CART_MAX];
        FINT k_nx[CART_MAX], k_ny[CART_MAX], k_nz[CART_MAX];
        FINT l_nx[CART_MAX], l_ny[CART_MAX], l_nz[CART_MAX];

        CINTcart_comp(i_nx, i_ny, i_nz, i_l);
        CINTcart_comp(j_nx, j_ny, j_nz, j_l);
        CINTcart_comp(k_nx, k_ny, k_nz, k_l);
        CINTcart_comp(l_nx, l_ny, l_nz, l_l);

        ofx = 0;
        ofy = envs->g_size;
        ofz = envs->g_size * 2;
        n = 0;
        for (j = 0; j < nfj; j++) {
                for (l = 0; l < nfl; l++) {
                        oflx = ofx + dj * j_nx[j] + dl * l_nx[l];
                        ofly = ofy + dj * j_ny[j] + dl * l_ny[l];
                        oflz = ofz + dj * j_nz[j] + dl * l_nz[l];
                        for (k = 0; k < nfk; k++) {
                                ofkx = oflx + dk * k_nx[k];
                                ofky = ofly + dk * k_ny[k];
                                ofkz = oflz + dk * k_nz[k];
                                switch (i_l) {
                                        case 0:
                                                idx[n+0] = ofkx;
                                                idx[n+1] = ofky;
                                                idx[n+2] = ofkz;
                                                n += 3;
                                                break;
                                        case 1:
                                                idx[n+0] = ofkx + di;
                                                idx[n+1] = ofky;
                                                idx[n+2] = ofkz;
                                                idx[n+3] = ofkx;
                                                idx[n+4] = ofky + di;
                                                idx[n+5] = ofkz;
                                                idx[n+6] = ofkx;
                                                idx[n+7] = ofky;
                                                idx[n+8] = ofkz + di;
                                                n += 9;
                                                break;
                                        case 2:
                                                idx[n+0 ] = ofkx + di*2;
                                                idx[n+1 ] = ofky;
                                                idx[n+2 ] = ofkz;
                                                idx[n+3 ] = ofkx + di;
                                                idx[n+4 ] = ofky + di;
                                                idx[n+5 ] = ofkz;
                                                idx[n+6 ] = ofkx + di;
                                                idx[n+7 ] = ofky;
                                                idx[n+8 ] = ofkz + di;
                                                idx[n+9 ] = ofkx;
                                                idx[n+10] = ofky + di*2;
                                                idx[n+11] = ofkz;
                                                idx[n+12] = ofkx;
                                                idx[n+13] = ofky + di;
                                                idx[n+14] = ofkz + di;
                                                idx[n+15] = ofkx;
                                                idx[n+16] = ofky;
                                                idx[n+17] = ofkz + di*2;
                                                n += 18;
                                                break;
                                        default:
                                                for (i = 0; i < nfi; i++) {
                                                        idx[n+0] = ofkx + di * i_nx[i]; //(:,ix,kx,lx,jx,1)
                                                        idx[n+1] = ofky + di * i_ny[i]; //(:,iy,ky,ly,jy,2)
                                                        idx[n+2] = ofkz + di * i_nz[i]; //(:,iz,kz,lz,jz,3)
                                                        n += 3;
                                                } // i
                                }
                        } // k
                } // l
        } // j
}


/*
 * g(nroots,0:nmax,0:mmax)
 */
void CINTg0_2e_2d(double *g, struct _BC *bc, CINTEnvVars *envs)
{
        const FINT nroots = envs->nrys_roots;
        const FINT nmax = envs->li_ceil + envs->lj_ceil;
        const FINT mmax = envs->lk_ceil + envs->ll_ceil;
        const FINT dm = envs->g2d_klmax;
        const FINT dn = envs->g2d_ijmax;
        FINT i, j, m, n, off;
        DEF_GXYZ(double, g, gx, gy, gz);
        const double *c00;
        const double *c0p;
        const double *b01 = bc->b01;
        const double *b00 = bc->b00;
        const double *b10 = bc->b10;
        double *p0x, *p0y, *p0z;
        const double *p1x, *p1y, *p1z, *p2x, *p2y, *p2z;

        for (i = 0; i < nroots; i++) {
                gx[i] = 1;
                gy[i] = 1;
                //gz[i] = w[i];
        }

        if (nmax > 0) {
                p0x = gx + dn;
                p0y = gy + dn;
                p0z = gz + dn;
                p1x = gx - dn;
                p1y = gy - dn;
                p1z = gz - dn;
                // gx(irys,0,1) = c00(irys) * gx(irys,0,0)
                for (c00 = bc->c00, i = 0; i < nroots; i++, c00+=3) {
                        p0x[i] = c00[0] * gx[i];
                        p0y[i] = c00[1] * gy[i];
                        p0z[i] = c00[2] * gz[i];
                }
                // gx(irys,0,n+1) = c00(irys)*gx(irys,0,n)
                // + n*b10(irys)*gx(irys,0,n-1)
                for (n = 1; n < nmax; n++) {
                        off = n * dn;
                        for (c00 = bc->c00, i = 0, j = off;
                             i < nroots; i++, j++, c00+=3) {
                                p0x[j] = c00[0] * gx[j] + n * b10[i] * p1x[j];
                                p0y[j] = c00[1] * gy[j] + n * b10[i] * p1y[j];
                                p0z[j] = c00[2] * gz[j] + n * b10[i] * p1z[j];
                        }
                }
        }

        if (mmax > 0) {
                p0x = gx + dm;
                p0y = gy + dm;
                p0z = gz + dm;
                p1x = gx - dm;
                p1y = gy - dm;
                p1z = gz - dm;
                // gx(irys,1,0) = c0p(irys) * gx(irys,0,0)
                for (c0p = bc->c0p, i = 0; i < nroots; i++, c0p+=3) {
                        p0x[i] = c0p[0] * gx[i];
                        p0y[i] = c0p[1] * gy[i];
                        p0z[i] = c0p[2] * gz[i];
                }
                // gx(irys,m+1,0) = c0p(irys)*gx(irys,m,0)
                // + m*b01(irys)*gx(irys,m-1,0)
                for (m = 1; m < mmax; m++) {
                        off = m * dm;
                        for (c0p = bc->c0p, i = 0, j = off;
                             i < nroots; i++, j++, c0p+=3) {
                                p0x[j] = c0p[0] * gx[j] + m * b01[i] * p1x[j];
                                p0y[j] = c0p[1] * gy[j] + m * b01[i] * p1y[j];
                                p0z[j] = c0p[2] * gz[j] + m * b01[i] * p1z[j];
                        }
                }
        }

        if (nmax > 0 && mmax > 0) {
                p0x = gx + dn;
                p0y = gy + dn;
                p0z = gz + dn;
                p1x = gx - dn;
                p1y = gy - dn;
                p1z = gz - dn;
                p2x = gx - dm;
                p2y = gy - dm;
                p2z = gz - dm;
                // gx(irys,1,1) = c0p(irys)*gx(irys,0,1)
                // + b00(irys)*gx(irys,0,0)
                for (c0p = bc->c0p, i = 0; i < nroots; i++, c0p+=3) {
                        p0x[i+dm] = c0p[0] * p0x[i] + b00[i] * gx[i];
                        p0y[i+dm] = c0p[1] * p0y[i] + b00[i] * gy[i];
                        p0z[i+dm] = c0p[2] * p0z[i] + b00[i] * gz[i];
                }

                // gx(irys,m+1,1) = c0p(irys)*gx(irys,m,1)
                // + m*b01(irys)*gx(irys,m-1,1)
                // + b00(irys)*gx(irys,m,0)
                for (m = 1; m < mmax; m++) {
                        off = m * dm + dn;
                        for (c0p = bc->c0p, i = 0, j = off;
                             i < nroots; i++, j++, c0p+=3) {
                                gx[j+dm] = c0p[0]*gx[j] + m*b01[i]*p2x[j] +b00[i]*p1x[j];
                                gy[j+dm] = c0p[1]*gy[j] + m*b01[i]*p2y[j] +b00[i]*p1y[j];
                                gz[j+dm] = c0p[2]*gz[j] + m*b01[i]*p2z[j] +b00[i]*p1z[j];
                        }
                }

                // gx(irys,m,n+1) = c00(irys)*gx(irys,m,n)
                // + n*b10(irys)*gx(irys,m,n-1)
                // + m*b00(irys)*gx(irys,m-1,n)
                for (m = 1; m <= mmax; m++) {
                        for (n = 1; n < nmax; n++) {
                                off = m * dm + n * dn;
                                for (c00 = bc->c00, i = 0, j = off;
                                     i < nroots; i++, j++, c00+=3) {
                                        p0x[j] = c00[0]*gx[j] +n*b10[i]*p1x[j] + m*b00[i]*p2x[j];
                                        p0y[j] = c00[1]*gy[j] +n*b10[i]*p1y[j] + m*b00[i]*p2y[j];
                                        p0z[j] = c00[2]*gz[j] +n*b10[i]*p1z[j] + m*b00[i]*p2z[j];
                                }
                        }
                }
        }
}


/*
 * g0[i,k,l,j] = < ik | lj > = ( i j | k l )
 */
/* 2d is based on l,j */
void CINTg0_lj2d_4d(double *restrict g, CINTEnvVars *envs)
{
        const FINT nmax = envs->li_ceil + envs->lj_ceil;
        const FINT mmax = envs->lk_ceil + envs->ll_ceil;
        const FINT li = envs->li_ceil;
        const FINT lk = envs->lk_ceil;
        //const FINT ll = envs->ll_ceil;
        const FINT lj = envs->lj_ceil;
        const FINT nroots = envs->nrys_roots;
        FINT i, j, k, l, ptr, n;
        const FINT di = envs->g_stride_i;
        const FINT dk = envs->g_stride_k;
        const FINT dl = envs->g_stride_l;
        const FINT dj = envs->g_stride_j;
        const double *rirj = envs->rirj;
        const double *rkrl = envs->rkrl;
        DEF_GXYZ(double, g, gx, gy, gz);
        const double *p1x, *p1y, *p1z, *p2x, *p2y, *p2z;
        double rx, ry, rz;

        // g(i,...,j) = rirj * g(i-1,...,j) +  g(i-1,...,j+1)
        rx = rirj[0];
        ry = rirj[1];
        rz = rirj[2];
        p1x = gx - di;
        p1y = gy - di;
        p1z = gz - di;
        p2x = gx - di + dj;
        p2y = gy - di + dj;
        p2z = gz - di + dj;
        for (i = 1; i <= li; i++) {
        for (j = 0; j <= nmax-i; j++) {
        for (l = 0; l <= mmax; l++) {
                ptr = j*dj + l*dl + i*di;
                for (n = ptr; n < ptr+nroots; n++) {
                        gx[n] = rx * p1x[n] + p2x[n];
                        gy[n] = ry * p1y[n] + p2y[n];
                        gz[n] = rz * p1z[n] + p2z[n];
                }
        } } }

        // g(...,k,l,..) = rkrl * g(...,k-1,l,..) + g(...,k-1,l+1,..)
        rx = rkrl[0];
        ry = rkrl[1];
        rz = rkrl[2];
        p1x = gx - dk;
        p1y = gy - dk;
        p1z = gz - dk;
        p2x = gx - dk + dl;
        p2y = gy - dk + dl;
        p2z = gz - dk + dl;
        for (j = 0; j <= lj; j++) {
        for (k = 1; k <= lk; k++) {
        for (l = 0; l <= mmax-k; l++) {
                ptr = j*dj + l*dl + k*dk;
                for (n = ptr; n < ptr+dk; n++) {
                        gx[n] = rx * p1x[n] + p2x[n];
                        gy[n] = ry * p1y[n] + p2y[n];
                        gz[n] = rz * p1z[n] + p2z[n];
                }
        } } }
}
/* 2d is based on k,j */
void CINTg0_kj2d_4d(double *restrict g, CINTEnvVars *envs)
{
        const FINT nmax = envs->li_ceil + envs->lj_ceil;
        const FINT mmax = envs->lk_ceil + envs->ll_ceil;
        const FINT li = envs->li_ceil;
        //const FINT lk = envs->lk_ceil;
        const FINT ll = envs->ll_ceil;
        const FINT lj = envs->lj_ceil;
        const FINT nroots = envs->nrys_roots;
        FINT i, j, k, l, ptr, n;
        const FINT di = envs->g_stride_i;
        const FINT dk = envs->g_stride_k;
        const FINT dl = envs->g_stride_l;
        const FINT dj = envs->g_stride_j;
        const double *rirj = envs->rirj;
        const double *rkrl = envs->rkrl;
        DEF_GXYZ(double, g, gx, gy, gz);
        const double *p1x, *p1y, *p1z, *p2x, *p2y, *p2z;
        double rx, ry, rz;

        // g(i,...,j) = rirj * g(i-1,...,j) +  g(i-1,...,j+1)
        rx = rirj[0];
        ry = rirj[1];
        rz = rirj[2];
        p1x = gx - di;
        p1y = gy - di;
        p1z = gz - di;
        p2x = gx - di + dj;
        p2y = gy - di + dj;
        p2z = gz - di + dj;
        for (i = 1; i <= li; i++) {
        for (j = 0; j <= nmax-i; j++) {
        for (k = 0; k <= mmax; k++) {
                ptr = j*dj + k*dk + i*di;
                for (n = ptr; n < ptr+nroots; n++) {
                        gx[n] = rx * p1x[n] + p2x[n];
                        gy[n] = ry * p1y[n] + p2y[n];
                        gz[n] = rz * p1z[n] + p2z[n];
                }
        } } }

        // g(...,k,l,..) = rkrl * g(...,k,l-1,..) + g(...,k+1,l-1,..)
        rx = rkrl[0];
        ry = rkrl[1];
        rz = rkrl[2];
        p1x = gx - dl;
        p1y = gy - dl;
        p1z = gz - dl;
        p2x = gx - dl + dk;
        p2y = gy - dl + dk;
        p2z = gz - dl + dk;
        for (j = 0; j <= lj; j++) {
        for (l = 1; l <= ll; l++) {
        for (k = 0; k <= mmax-l; k++) {
                ptr = j*dj + l*dl + k*dk;
                for (n = ptr; n < ptr+dk; n++) {
                        gx[n] = rx * p1x[n] + p2x[n];
                        gy[n] = ry * p1y[n] + p2y[n];
                        gz[n] = rz * p1z[n] + p2z[n];
                }
        } } }
}
/* 2d is based on i,l */
void CINTg0_il2d_4d(double *restrict g, CINTEnvVars *envs)
{
        const FINT nmax = envs->li_ceil + envs->lj_ceil;
        const FINT mmax = envs->lk_ceil + envs->ll_ceil;
        //const FINT li = envs->li_ceil;
        const FINT lk = envs->lk_ceil;
        const FINT ll = envs->ll_ceil;
        const FINT lj = envs->lj_ceil;
        const FINT nroots = envs->nrys_roots;
        FINT i, j, k, l, ptr, n;
        const FINT di = envs->g_stride_i;
        const FINT dk = envs->g_stride_k;
        const FINT dl = envs->g_stride_l;
        const FINT dj = envs->g_stride_j;
        const double *rirj = envs->rirj;
        const double *rkrl = envs->rkrl;
        DEF_GXYZ(double, g, gx, gy, gz);
        const double *p1x, *p1y, *p1z, *p2x, *p2y, *p2z;
        double rx, ry, rz;

        // g(...,k,l,..) = rkrl * g(...,k-1,l,..) + g(...,k-1,l+1,..)
        rx = rkrl[0];
        ry = rkrl[1];
        rz = rkrl[2];
        p1x = gx - dk;
        p1y = gy - dk;
        p1z = gz - dk;
        p2x = gx - dk + dl;
        p2y = gy - dk + dl;
        p2z = gz - dk + dl;
        for (k = 1; k <= lk; k++) {
        for (l = 0; l <= mmax-k; l++) {
        for (i = 0; i <= nmax; i++) {
                ptr = l*dl + k*dk + i*di;
                for (n = ptr; n < ptr+nroots; n++) {
                        gx[n] = rx * p1x[n] + p2x[n];
                        gy[n] = ry * p1y[n] + p2y[n];
                        gz[n] = rz * p1z[n] + p2z[n];
                }
        } } }

        // g(i,...,j) = rirj * g(i,...,j-1) +  g(i+1,...,j-1)
        rx = rirj[0];
        ry = rirj[1];
        rz = rirj[2];
        p1x = gx - dj;
        p1y = gy - dj;
        p1z = gz - dj;
        p2x = gx - dj + di;
        p2y = gy - dj + di;
        p2z = gz - dj + di;
        for (j = 1; j <= lj; j++) {
        for (l = 0; l <= ll; l++) {
        for (k = 0; k <= lk; k++) {
                ptr = j*dj + l*dl + k*dk;
                for (n = ptr; n < ptr+dk-di*j; n++) {
                        gx[n] = rx * p1x[n] + p2x[n];
                        gy[n] = ry * p1y[n] + p2y[n];
                        gz[n] = rz * p1z[n] + p2z[n];
                }
        } } }
}
/* 2d is based on i,k */
void CINTg0_ik2d_4d(double *restrict g, CINTEnvVars *envs)
{
        const FINT nmax = envs->li_ceil + envs->lj_ceil;
        const FINT mmax = envs->lk_ceil + envs->ll_ceil;
        //const FINT li = envs->li_ceil;
        const FINT lk = envs->lk_ceil;
        const FINT ll = envs->ll_ceil;
        const FINT lj = envs->lj_ceil;
        const FINT nroots = envs->nrys_roots;
        FINT i, j, k, l, ptr, n;
        const FINT di = envs->g_stride_i;
        const FINT dk = envs->g_stride_k;
        const FINT dl = envs->g_stride_l;
        const FINT dj = envs->g_stride_j;
        const double *rirj = envs->rirj;
        const double *rkrl = envs->rkrl;
        DEF_GXYZ(double, g, gx, gy, gz);
        const double *p1x, *p1y, *p1z, *p2x, *p2y, *p2z;
        double rx, ry, rz;

        // g(...,k,l,..) = rkrl * g(...,k,l-1,..) + g(...,k+1,l-1,..)
        rx = rkrl[0];
        ry = rkrl[1];
        rz = rkrl[2];
        p1x = gx - dl;
        p1y = gy - dl;
        p1z = gz - dl;
        p2x = gx - dl + dk;
        p2y = gy - dl + dk;
        p2z = gz - dl + dk;
        for (l = 1; l <= ll; l++) {
                // (:,i) is full, so loop:k and loop:n can be merged to
                // for(n = l*dl; n < ptr+dl-dk*l; n++)
                for (k = 0; k <= mmax-l; k++) {
                for (i = 0; i <= nmax; i++) {
                        ptr = l*dl + k*dk + i*di;
                        for (n = ptr; n < ptr+nroots; n++) {
                                gx[n] = rx * p1x[n] + p2x[n];
                                gy[n] = ry * p1y[n] + p2y[n];
                                gz[n] = rz * p1z[n] + p2z[n];
                        }
                } }
        }

        // g(i,...,j) = rirj * g(i,...,j-1) +  g(i+1,...,j-1)
        rx = rirj[0];
        ry = rirj[1];
        rz = rirj[2];
        p1x = gx - dj;
        p1y = gy - dj;
        p1z = gz - dj;
        p2x = gx - dj + di;
        p2y = gy - dj + di;
        p2z = gz - dj + di;
        for (j = 1; j <= lj; j++) {
        for (l = 0; l <= ll; l++) {
        for (k = 0; k <= lk; k++) {
                ptr = j*dj + l*dl + k*dk;
                for (n = ptr; n < ptr+dk-di*j; n++) {
                        gx[n] = rx * p1x[n] + p2x[n];
                        gy[n] = ry * p1y[n] + p2y[n];
                        gz[n] = rz * p1z[n] + p2z[n];
                }
        } } }
}

static inline void _g0_2d4d_0000(double *restrict g, struct _BC *bc, CINTEnvVars *envs)
{
        g[0] = 1;
        g[1] = 1;
        //g[2] = w[0];
}

static inline void _g0_2d4d_0001(double *restrict g, struct _BC *bc, CINTEnvVars *envs)
{
        double *cp = bc->c0p;
        g[0] = 1;
        g[1] = cp[0];
        g[2] = 1;
        g[3] = cp[1];
        //g[4] = w[0];
        g[5] = cp[2] * g[4];
}

static inline void _g0_2d4d_0002(double *restrict g, struct _BC *bc, CINTEnvVars *envs)
{
        double *cp = bc->c0p;
        double *b01 = bc->b01;
        g[0] = 1;
        g[1] = 1;
        g[2] = cp[0];
        g[3] = cp[3];
        g[4] = cp[0] * cp[0] + b01[0];
        g[5] = cp[3] * cp[3] + b01[1];
        g[6] = 1;
        g[7] = 1;
        g[8] = cp[1];
        g[9] = cp[4];
        g[10] = cp[1] * cp[1] + b01[0];
        g[11] = cp[4] * cp[4] + b01[1];
        //g[12] = w[0];
        //g[13] = w[1];
        g[14] = cp[2] * g[12];
        g[15] = cp[5] * g[13];
        g[16] = cp[2] * g[14] + b01[0] * g[12];
        g[17] = cp[5] * g[15] + b01[1] * g[13];
}

static inline void _g0_2d4d_0003(double *restrict g, struct _BC *bc, CINTEnvVars *envs)
{
        double *cp = bc->c0p;
        double *b01 = bc->b01;
        g[0] = 1;
        g[1] = 1;
        g[2] = cp[0];
        g[3] = cp[3];
        g[4] = cp[0] * cp[0] + b01[0];
        g[5] = cp[3] * cp[3] + b01[1];
        g[6] = cp[0] * (g[4] + 2 * b01[0]);
        g[7] = cp[3] * (g[5] + 2 * b01[1]);
        g[8] = 1;
        g[9] = 1;
        g[10] = cp[1];
        g[11] = cp[4];
        g[12] = cp[1] * cp[1] + b01[0];
        g[13] = cp[4] * cp[4] + b01[1];
        g[14] = cp[1] * (g[12] + 2 * b01[0]);
        g[15] = cp[4] * (g[13] + 2 * b01[1]);
        //g[16] = w[0];
        //g[17] = w[1];
        g[18] = cp[2] * g[16];
        g[19] = cp[5] * g[17];
        g[20] = cp[2] * g[18] + b01[0] * g[16];
        g[21] = cp[5] * g[19] + b01[1] * g[17];
        g[22] = cp[2] * g[20] + 2 * b01[0] * g[18];
        g[23] = cp[5] * g[21] + 2 * b01[1] * g[19];
}

static inline void _g0_2d4d_0010(double *restrict g, struct _BC *bc, CINTEnvVars *envs)
{
        double *cp = bc->c0p;
        g[0] = 1;
        g[1] = cp[0];
        g[2] = 1;
        g[3] = cp[1];
        //g[4] = w[0];
        g[5] = cp[2] * g[4];
}

static inline void _g0_2d4d_0011(double *restrict g, struct _BC *bc, CINTEnvVars *envs)
{
        double *cp = bc->c0p;
        double *b01 = bc->b01;
        double xkxl = envs->rkrl[0];
        double ykyl = envs->rkrl[1];
        double zkzl = envs->rkrl[2];
        g[0] = 1;
        g[1] = 1;
        g[4] = cp[0];
        g[5] = cp[3];
        g[6] = cp[0] * (xkxl + cp[0]) + b01[0];
        g[7] = cp[3] * (xkxl + cp[3]) + b01[1];
        g[2] = xkxl + cp[0];
        g[3] = xkxl + cp[3];
        g[12] = 1;
        g[13] = 1;
        g[16] = cp[1];
        g[17] = cp[4];
        g[18] = cp[1] * (ykyl + cp[1]) + b01[0];
        g[19] = cp[4] * (ykyl + cp[4]) + b01[1];
        g[14] = ykyl + cp[1];
        g[15] = ykyl + cp[4];
        //g[24] = w[0];
        //g[25] = w[1];
        g[28] = cp[2] * g[24];
        g[29] = cp[5] * g[25];
        g[30] = g[28] * (zkzl + cp[2]) + b01[0] * g[24];
        g[31] = g[29] * (zkzl + cp[5]) + b01[1] * g[25];
        g[26] = g[24] * (zkzl + cp[2]);
        g[27] = g[25] * (zkzl + cp[5]);
}

static inline void _g0_2d4d_0012(double *restrict g, struct _BC *bc, CINTEnvVars *envs)
{
        double *cp = bc->c0p;
        double *b01 = bc->b01;
        double xkxl = envs->rkrl[0];
        double ykyl = envs->rkrl[1];
        double zkzl = envs->rkrl[2];
        g[0] = 1;
        g[1] = 1;
        g[4] = cp[0];
        g[5] = cp[3];
        g[8] = cp[0] * cp[0] + b01[0];
        g[9] = cp[3] * cp[3] + b01[1];
        g[10] = g[8] * (xkxl + cp[0]) + cp[0] * 2 * b01[0];
        g[11] = g[9] * (xkxl + cp[3]) + cp[3] * 2 * b01[1];
        g[6] = cp[0] * (xkxl + cp[0]) + b01[0];
        g[7] = cp[3] * (xkxl + cp[3]) + b01[1];
        g[2] = xkxl + cp[0];
        g[3] = xkxl + cp[3];
        g[16] = 1;
        g[17] = 1;
        g[20] = cp[1];
        g[21] = cp[4];
        g[24] = cp[1] * cp[1] + b01[0];
        g[25] = cp[4] * cp[4] + b01[1];
        g[26] = g[24] * (ykyl + cp[1]) + cp[1] * 2 * b01[0];
        g[27] = g[25] * (ykyl + cp[4]) + cp[4] * 2 * b01[1];
        g[22] = cp[1] * (ykyl + cp[1]) + b01[0];
        g[23] = cp[4] * (ykyl + cp[4]) + b01[1];
        g[18] = ykyl + cp[1];
        g[19] = ykyl + cp[4];
        //g[32] = w[0];
        //g[33] = w[1];
        g[36] = cp[2] * g[32];
        g[37] = cp[5] * g[33];
        g[40] = cp[2] * g[36] + b01[0] * g[32];
        g[41] = cp[5] * g[37] + b01[1] * g[33];
        g[42] = g[40] * (zkzl + cp[2]) + 2 * b01[0] * g[36];
        g[43] = g[41] * (zkzl + cp[5]) + 2 * b01[1] * g[37];
        g[38] = g[36] * (zkzl + cp[2]) + b01[0] * g[32];
        g[39] = g[37] * (zkzl + cp[5]) + b01[1] * g[33];
        g[34] = g[32] * (zkzl + cp[2]);
        g[35] = g[33] * (zkzl + cp[5]);
}

static inline void _g0_2d4d_0020(double *restrict g, struct _BC *bc, CINTEnvVars *envs)
{
        double *cp = bc->c0p;
        double *b01 = bc->b01;
        g[0] = 1;
        g[1] = 1;
        g[2] = cp[0];
        g[3] = cp[3];
        g[4] = cp[0] * cp[0] + b01[0];
        g[5] = cp[3] * cp[3] + b01[1];
        g[6] = 1;
        g[7] = 1;
        g[8] = cp[1];
        g[9] = cp[4];
        g[10] = cp[1] * cp[1] + b01[0];
        g[11] = cp[4] * cp[4] + b01[1];
        //g[12] = w[0];
        //g[13] = w[1];
        g[14] = cp[2] * g[12];
        g[15] = cp[5] * g[13];
        g[16] = cp[2] * g[14] + b01[0] * g[12];
        g[17] = cp[5] * g[15] + b01[1] * g[13];
}

static inline void _g0_2d4d_0021(double *restrict g, struct _BC *bc, CINTEnvVars *envs)
{
        double *cp = bc->c0p;
        double *b01 = bc->b01;
        double xkxl = envs->rkrl[0];
        double ykyl = envs->rkrl[1];
        double zkzl = envs->rkrl[2];
        g[0] = 1;
        g[1] = 1;
        g[2] = cp[0];
        g[3] = cp[3];
        g[4] = cp[0] * cp[0] + b01[0];
        g[5] = cp[3] * cp[3] + b01[1];
        g[12] = g[4] * (xkxl + cp[0]) + cp[0] * 2 * b01[0];
        g[13] = g[5] * (xkxl + cp[3]) + cp[3] * 2 * b01[1];
        g[10] = cp[0] * (xkxl + cp[0]) + b01[0];
        g[11] = cp[3] * (xkxl + cp[3]) + b01[1];
        g[8] = xkxl + cp[0];
        g[9] = xkxl + cp[3];
        g[16] = 1;
        g[17] = 1;
        g[18] = cp[1];
        g[19] = cp[4];
        g[20] = cp[1] * cp[1] + b01[0];
        g[21] = cp[4] * cp[4] + b01[1];
        g[28] = g[20] * (ykyl + cp[1]) + cp[1] * 2 * b01[0];
        g[29] = g[21] * (ykyl + cp[4]) + cp[4] * 2 * b01[1];
        g[26] = cp[1] * (ykyl + cp[1]) + b01[0];
        g[27] = cp[4] * (ykyl + cp[4]) + b01[1];
        g[24] = ykyl + cp[1];
        g[25] = ykyl + cp[4];
        //g[32] = w[0];
        //g[33] = w[1];
        g[34] = cp[2] * g[32];
        g[35] = cp[5] * g[33];
        g[36] = cp[2] * g[34] + b01[0] * g[32];
        g[37] = cp[5] * g[35] + b01[1] * g[33];
        g[44] = g[36] * (zkzl + cp[2]) + 2 * b01[0] * g[34];
        g[45] = g[37] * (zkzl + cp[5]) + 2 * b01[1] * g[35];
        g[42] = g[34] * (zkzl + cp[2]) + b01[0] * g[32];
        g[43] = g[35] * (zkzl + cp[5]) + b01[1] * g[33];
        g[40] = g[32] * (zkzl + cp[2]);
        g[41] = g[33] * (zkzl + cp[5]);
}

static inline void _g0_2d4d_0030(double *restrict g, struct _BC *bc, CINTEnvVars *envs)
{
        double *cp = bc->c0p;
        double *b01 = bc->b01;
        g[0] = 1;
        g[1] = 1;
        g[2] = cp[0];
        g[3] = cp[3];
        g[4] = cp[0] * cp[0] + b01[0];
        g[5] = cp[3] * cp[3] + b01[1];
        g[6] = cp[0] * (g[4] + 2 * b01[0]);
        g[7] = cp[3] * (g[5] + 2 * b01[1]);
        g[8] = 1;
        g[9] = 1;
        g[10] = cp[1];
        g[11] = cp[4];
        g[12] = cp[1] * cp[1] + b01[0];
        g[13] = cp[4] * cp[4] + b01[1];
        g[14] = cp[1] * (g[12] + 2 * b01[0]);
        g[15] = cp[4] * (g[13] + 2 * b01[1]);
        //g[16] = w[0];
        //g[17] = w[1];
        g[18] = cp[2] * g[16];
        g[19] = cp[5] * g[17];
        g[20] = cp[2] * g[18] + b01[0] * g[16];
        g[21] = cp[5] * g[19] + b01[1] * g[17];
        g[22] = cp[2] * g[20] + 2 * b01[0] * g[18];
        g[23] = cp[5] * g[21] + 2 * b01[1] * g[19];
}

static inline void _g0_2d4d_0100(double *restrict g, struct _BC *bc, CINTEnvVars *envs)
{
        double *c0 = bc->c00;
        g[0] = 1;
        g[1] = c0[0];
        g[2] = 1;
        g[3] = c0[1];
        //g[4] = w[0];
        g[5] = c0[2] * g[4];
}

static inline void _g0_2d4d_0101(double *restrict g, struct _BC *bc, CINTEnvVars *envs)
{
        double *c0 = bc->c00;
        double *cp = bc->c0p;
        double *b00 = bc->b00;
        g[0] = 1;
        g[1] = 1;
        g[4] = c0[0];
        g[5] = c0[3];
        g[2] = cp[0];
        g[3] = cp[3];
        g[6] = cp[0] * c0[0] + b00[0];
        g[7] = cp[3] * c0[3] + b00[1];
        g[8] = 1;
        g[9] = 1;
        g[12] = c0[1];
        g[13] = c0[4];
        g[10] = cp[1];
        g[11] = cp[4];
        g[14] = cp[1] * c0[1] + b00[0];
        g[15] = cp[4] * c0[4] + b00[1];
        //g[16] = w[0];
        //g[17] = w[1];
        g[20] = c0[2] * g[16];
        g[21] = c0[5] * g[17];
        g[18] = cp[2] * g[16];
        g[19] = cp[5] * g[17];
        g[22] = cp[2] * g[20] + b00[0] * g[16];
        g[23] = cp[5] * g[21] + b00[1] * g[17];
}

static inline void _g0_2d4d_0102(double *restrict g, struct _BC *bc, CINTEnvVars *envs)
{
        double *c0 = bc->c00;
        double *cp = bc->c0p;
        double *b00 = bc->b00;
        double *b01 = bc->b01;
        g[0] = 1;
        g[1] = 1;
        g[6] = c0[0];
        g[7] = c0[3];
        g[2] = cp[0];
        g[3] = cp[3];
        g[4] = cp[0] * cp[0] + b01[0];
        g[5] = cp[3] * cp[3] + b01[1];
        g[8] = cp[0] * c0[0] + b00[0];
        g[9] = cp[3] * c0[3] + b00[1];
        g[10] = cp[0] * (g[8] + b00[0]) + b01[0] * c0[0];
        g[11] = cp[3] * (g[9] + b00[1]) + b01[1] * c0[3];
        g[12] = 1;
        g[13] = 1;
        g[18] = c0[1];
        g[19] = c0[4];
        g[14] = cp[1];
        g[15] = cp[4];
        g[16] = cp[1] * cp[1] + b01[0];
        g[17] = cp[4] * cp[4] + b01[1];
        g[20] = cp[1] * c0[1] + b00[0];
        g[21] = cp[4] * c0[4] + b00[1];
        g[22] = cp[1] * (g[20] + b00[0]) + b01[0] * c0[1];
        g[23] = cp[4] * (g[21] + b00[1]) + b01[1] * c0[4];
        //g[24] = w[0];
        //g[25] = w[1];
        g[30] = c0[2] * g[24];
        g[31] = c0[5] * g[25];
        g[26] = cp[2] * g[24];
        g[27] = cp[5] * g[25];
        g[28] = cp[2] * g[26] + b01[0] * g[24];
        g[29] = cp[5] * g[27] + b01[1] * g[25];
        g[32] = cp[2] * g[30] + b00[0] * g[24];
        g[33] = cp[5] * g[31] + b00[1] * g[25];
        g[34] = cp[2] * g[32] + b01[0] * g[30] + b00[0] * g[26];
        g[35] = cp[5] * g[33] + b01[1] * g[31] + b00[1] * g[27];
}

static inline void _g0_2d4d_0110(double *restrict g, struct _BC *bc, CINTEnvVars *envs)
{
        double *c0 = bc->c00;
        double *cp = bc->c0p;
        double *b00 = bc->b00;
        g[0] = 1;
        g[1] = 1;
        g[4] = c0[0];
        g[5] = c0[3];
        g[2] = cp[0];
        g[3] = cp[3];
        g[6] = cp[0] * c0[0] + b00[0];
        g[7] = cp[3] * c0[3] + b00[1];
        g[8] = 1;
        g[9] = 1;
        g[12] = c0[1];
        g[13] = c0[4];
        g[10] = cp[1];
        g[11] = cp[4];
        g[14] = cp[1] * c0[1] + b00[0];
        g[15] = cp[4] * c0[4] + b00[1];
        //g[16] = w[0];
        //g[17] = w[1];
        g[20] = c0[2] * g[16];
        g[21] = c0[5] * g[17];
        g[18] = cp[2] * g[16];
        g[19] = cp[5] * g[17];
        g[22] = cp[2] * g[20] + b00[0] * g[16];
        g[23] = cp[5] * g[21] + b00[1] * g[17];
}

static inline void _g0_2d4d_0111(double *restrict g, struct _BC *bc, CINTEnvVars *envs)
{
        double *c0 = bc->c00;
        double *cp = bc->c0p;
        double *b00 = bc->b00;
        double *b01 = bc->b01;
        double xkxl = envs->rkrl[0];
        double ykyl = envs->rkrl[1];
        double zkzl = envs->rkrl[2];
        g[0] = 1;
        g[1] = 1;
        g[12] = c0[0];
        g[13] = c0[3];
        g[4] = cp[0];
        g[5] = cp[3];
        g[16] = cp[0] * c0[0] + b00[0];
        g[17] = cp[3] * c0[3] + b00[1];
        g[6] = cp[0] * (xkxl + cp[0]) + b01[0];
        g[7] = cp[3] * (xkxl + cp[3]) + b01[1];
        g[18] = g[16] * (xkxl + cp[0]) + cp[0] * b00[0] + b01[0] * c0[0];
        g[19] = g[17] * (xkxl + cp[3]) + cp[3] * b00[1] + b01[1] * c0[3];
        g[2] = xkxl + cp[0];
        g[3] = xkxl + cp[3];
        g[14] = c0[0] * (xkxl + cp[0]) + b00[0];
        g[15] = c0[3] * (xkxl + cp[3]) + b00[1];
        g[24] = 1;
        g[25] = 1;
        g[36] = c0[1];
        g[37] = c0[4];
        g[28] = cp[1];
        g[29] = cp[4];
        g[40] = cp[1] * c0[1] + b00[0];
        g[41] = cp[4] * c0[4] + b00[1];
        g[30] = cp[1] * (ykyl + cp[1]) + b01[0];
        g[31] = cp[4] * (ykyl + cp[4]) + b01[1];
        g[42] = g[40] * (ykyl + cp[1]) + cp[1] * b00[0] + b01[0] * c0[1];
        g[43] = g[41] * (ykyl + cp[4]) + cp[4] * b00[1] + b01[1] * c0[4];
        g[26] = ykyl + cp[1];
        g[27] = ykyl + cp[4];
        g[38] = c0[1] * (ykyl + cp[1]) + b00[0];
        g[39] = c0[4] * (ykyl + cp[4]) + b00[1];
        //g[48] = w[0];
        //g[49] = w[1];
        g[60] = c0[2] * g[48];
        g[61] = c0[5] * g[49];
        g[52] = cp[2] * g[48];
        g[53] = cp[5] * g[49];
        g[64] = cp[2] * g[60] + b00[0] * g[48];
        g[65] = cp[5] * g[61] + b00[1] * g[49];
        g[54] = g[52] * (zkzl + cp[2]) + b01[0] * g[48];
        g[55] = g[53] * (zkzl + cp[5]) + b01[1] * g[49];
        g[66] = g[64] * (zkzl + cp[2]) + b01[0] * g[60] + b00[0] * g[52];
        g[67] = g[65] * (zkzl + cp[5]) + b01[1] * g[61] + b00[1] * g[53];
        g[50] = g[48] * (zkzl + cp[2]);
        g[51] = g[49] * (zkzl + cp[5]);
        g[62] = g[60] * (zkzl + cp[2]) + b00[0] * g[48];
        g[63] = g[61] * (zkzl + cp[5]) + b00[1] * g[49];
}

static inline void _g0_2d4d_0120(double *restrict g, struct _BC *bc, CINTEnvVars *envs)
{
        double *c0 = bc->c00;
        double *cp = bc->c0p;
        double *b00 = bc->b00;
        double *b01 = bc->b01;
        g[0] = 1;
        g[1] = 1;
        g[6] = c0[0];
        g[7] = c0[3];
        g[2] = cp[0];
        g[3] = cp[3];
        g[4] = cp[0] * cp[0] + b01[0];
        g[5] = cp[3] * cp[3] + b01[1];
        g[8] = cp[0] * c0[0] + b00[0];
        g[9] = cp[3] * c0[3] + b00[1];
        g[10] = cp[0] * (g[8] + b00[0]) + b01[0] * c0[0];
        g[11] = cp[3] * (g[9] + b00[1]) + b01[1] * c0[3];
        g[12] = 1;
        g[13] = 1;
        g[18] = c0[1];
        g[19] = c0[4];
        g[14] = cp[1];
        g[15] = cp[4];
        g[16] = cp[1] * cp[1] + b01[0];
        g[17] = cp[4] * cp[4] + b01[1];
        g[20] = cp[1] * c0[1] + b00[0];
        g[21] = cp[4] * c0[4] + b00[1];
        g[22] = cp[1] * (g[20] + b00[0]) + b01[0] * c0[1];
        g[23] = cp[4] * (g[21] + b00[1]) + b01[1] * c0[4];
        //g[24] = w[0];
        //g[25] = w[1];
        g[30] = c0[2] * g[24];
        g[31] = c0[5] * g[25];
        g[26] = cp[2] * g[24];
        g[27] = cp[5] * g[25];
        g[28] = cp[2] * g[26] + b01[0] * g[24];
        g[29] = cp[5] * g[27] + b01[1] * g[25];
        g[32] = cp[2] * g[30] + b00[0] * g[24];
        g[33] = cp[5] * g[31] + b00[1] * g[25];
        g[34] = cp[2] * g[32] + b01[0] * g[30] + b00[0] * g[26];
        g[35] = cp[5] * g[33] + b01[1] * g[31] + b00[1] * g[27];
}

static inline void _g0_2d4d_0200(double *restrict g, struct _BC *bc, CINTEnvVars *envs)
{
        double *c0 = bc->c00;
        double *b10 = bc->b10;
        g[0] = 1;
        g[1] = 1;
        g[2] = c0[0];
        g[3] = c0[3];
        g[4] = c0[0] * c0[0] + b10[0];
        g[5] = c0[3] * c0[3] + b10[1];
        g[6] = 1;
        g[7] = 1;
        g[8] = c0[1];
        g[9] = c0[4];
        g[10] = c0[1] * c0[1] + b10[0];
        g[11] = c0[4] * c0[4] + b10[1];
        //g[12] = w[0];
        //g[13] = w[1];
        g[14] = c0[2] * g[12];
        g[15] = c0[5] * g[13];
        g[16] = c0[2] * g[14] + b10[0] * g[12];
        g[17] = c0[5] * g[15] + b10[1] * g[13];
}

static inline void _g0_2d4d_0201(double *restrict g, struct _BC *bc, CINTEnvVars *envs)
{
        double *c0 = bc->c00;
        double *cp = bc->c0p;
        double *b00 = bc->b00;
        double *b10 = bc->b10;
        g[0] = 1;
        g[1] = 1;
        g[4] = c0[0];
        g[5] = c0[3];
        g[8] = c0[0] * c0[0] + b10[0];
        g[9] = c0[3] * c0[3] + b10[1];
        g[2] = cp[0];
        g[3] = cp[3];
        g[6] = cp[0] * c0[0] + b00[0];
        g[7] = cp[3] * c0[3] + b00[1];
        g[10] = c0[0] * (g[6] + b00[0]) + b10[0] * cp[0];
        g[11] = c0[3] * (g[7] + b00[1]) + b10[1] * cp[3];
        g[12] = 1;
        g[13] = 1;
        g[16] = c0[1];
        g[17] = c0[4];
        g[20] = c0[1] * c0[1] + b10[0];
        g[21] = c0[4] * c0[4] + b10[1];
        g[14] = cp[1];
        g[15] = cp[4];
        g[18] = cp[1] * c0[1] + b00[0];
        g[19] = cp[4] * c0[4] + b00[1];
        g[22] = c0[1] * (g[18] + b00[0]) + b10[0] * cp[1];
        g[23] = c0[4] * (g[19] + b00[1]) + b10[1] * cp[4];
        //g[24] = w[0];
        //g[25] = w[1];
        g[28] = c0[2] * g[24];
        g[29] = c0[5] * g[25];
        g[32] = c0[2] * g[28] + b10[0] * g[24];
        g[33] = c0[5] * g[29] + b10[1] * g[25];
        g[26] = cp[2] * g[24];
        g[27] = cp[5] * g[25];
        g[30] = cp[2] * g[28] + b00[0] * g[24];
        g[31] = cp[5] * g[29] + b00[1] * g[25];
        g[34] = c0[2] * g[30] + b10[0] * g[26] + b00[0] * g[28];
        g[35] = c0[5] * g[31] + b10[1] * g[27] + b00[1] * g[29];
}

static inline void _g0_2d4d_0210(double *restrict g, struct _BC *bc, CINTEnvVars *envs)
{
        double *c0 = bc->c00;
        double *cp = bc->c0p;
        double *b00 = bc->b00;
        double *b10 = bc->b10;
        g[0] = 1;
        g[1] = 1;
        g[4] = c0[0];
        g[5] = c0[3];
        g[8] = c0[0] * c0[0] + b10[0];
        g[9] = c0[3] * c0[3] + b10[1];
        g[2] = cp[0];
        g[3] = cp[3];
        g[6] = cp[0] * c0[0] + b00[0];
        g[7] = cp[3] * c0[3] + b00[1];
        g[10] = c0[0] * (g[6] + b00[0]) + b10[0] * cp[0];
        g[11] = c0[3] * (g[7] + b00[1]) + b10[1] * cp[3];
        g[12] = 1;
        g[13] = 1;
        g[16] = c0[1];
        g[17] = c0[4];
        g[20] = c0[1] * c0[1] + b10[0];
        g[21] = c0[4] * c0[4] + b10[1];
        g[14] = cp[1];
        g[15] = cp[4];
        g[18] = cp[1] * c0[1] + b00[0];
        g[19] = cp[4] * c0[4] + b00[1];
        g[22] = c0[1] * (g[18] + b00[0]) + b10[0] * cp[1];
        g[23] = c0[4] * (g[19] + b00[1]) + b10[1] * cp[4];
        //g[24] = w[0];
        //g[25] = w[1];
        g[28] = c0[2] * g[24];
        g[29] = c0[5] * g[25];
        g[32] = c0[2] * g[28] + b10[0] * g[24];
        g[33] = c0[5] * g[29] + b10[1] * g[25];
        g[26] = cp[2] * g[24];
        g[27] = cp[5] * g[25];
        g[30] = cp[2] * g[28] + b00[0] * g[24];
        g[31] = cp[5] * g[29] + b00[1] * g[25];
        g[34] = c0[2] * g[30] + b10[0] * g[26] + b00[0] * g[28];
        g[35] = c0[5] * g[31] + b10[1] * g[27] + b00[1] * g[29];
}

static inline void _g0_2d4d_0300(double *restrict g, struct _BC *bc, CINTEnvVars *envs)
{
        double *c0 = bc->c00;
        double *b10 = bc->b10;
        g[0] = 1;
        g[1] = 1;
        g[2] = c0[0];
        g[3] = c0[3];
        g[4] = c0[0] * c0[0] + b10[0];
        g[5] = c0[3] * c0[3] + b10[1];
        g[6] = c0[0] * (g[4] + 2 * b10[0]);
        g[7] = c0[3] * (g[5] + 2 * b10[1]);
        g[8] = 1;
        g[9] = 1;
        g[10] = c0[1];
        g[11] = c0[4];
        g[12] = c0[1] * c0[1] + b10[0];
        g[13] = c0[4] * c0[4] + b10[1];
        g[14] = c0[1] * (g[12] + 2 * b10[0]);
        g[15] = c0[4] * (g[13] + 2 * b10[1]);
        //g[16] = w[0];
        //g[17] = w[1];
        g[18] = c0[2] * g[16];
        g[19] = c0[5] * g[17];
        g[20] = c0[2] * g[18] + b10[0] * g[16];
        g[21] = c0[5] * g[19] + b10[1] * g[17];
        g[22] = c0[2] * g[20] + 2 * b10[0] * g[18];
        g[23] = c0[5] * g[21] + 2 * b10[1] * g[19];
}

static inline void _g0_2d4d_1000(double *restrict g, struct _BC *bc, CINTEnvVars *envs)
{
        double *c0 = bc->c00;
        g[0] = 1;
        g[1] = c0[0];
        g[2] = 1;
        g[3] = c0[1];
        //g[4] = w[0];
        g[5] = c0[2] * g[4];
}

static inline void _g0_2d4d_1001(double *restrict g, struct _BC *bc, CINTEnvVars *envs)
{
        double *c0 = bc->c00;
        double *cp = bc->c0p;
        double *b00 = bc->b00;
        g[0] = 1;
        g[1] = 1;
        g[2] = c0[0];
        g[3] = c0[3];
        g[4] = cp[0];
        g[5] = cp[3];
        g[6] = cp[0] * c0[0] + b00[0];
        g[7] = cp[3] * c0[3] + b00[1];
        g[8] = 1;
        g[9] = 1;
        g[10] = c0[1];
        g[11] = c0[4];
        g[12] = cp[1];
        g[13] = cp[4];
        g[14] = cp[1] * c0[1] + b00[0];
        g[15] = cp[4] * c0[4] + b00[1];
        //g[16] = w[0];
        //g[17] = w[1];
        g[18] = c0[2] * g[16];
        g[19] = c0[5] * g[17];
        g[20] = cp[2] * g[16];
        g[21] = cp[5] * g[17];
        g[22] = cp[2] * g[18] + b00[0] * g[16];
        g[23] = cp[5] * g[19] + b00[1] * g[17];
}

static inline void _g0_2d4d_1002(double *restrict g, struct _BC *bc, CINTEnvVars *envs)
{
        double *c0 = bc->c00;
        double *cp = bc->c0p;
        double *b00 = bc->b00;
        double *b01 = bc->b01;
        g[0] = 1;
        g[1] = 1;
        g[2] = c0[0];
        g[3] = c0[3];
        g[4] = cp[0];
        g[5] = cp[3];
        g[8] = cp[0] * cp[0] + b01[0];
        g[9] = cp[3] * cp[3] + b01[1];
        g[6] = cp[0] * c0[0] + b00[0];
        g[7] = cp[3] * c0[3] + b00[1];
        g[10] = cp[0] * (g[6] + b00[0]) + b01[0] * c0[0];
        g[11] = cp[3] * (g[7] + b00[1]) + b01[1] * c0[3];
        g[12] = 1;
        g[13] = 1;
        g[14] = c0[1];
        g[15] = c0[4];
        g[16] = cp[1];
        g[17] = cp[4];
        g[20] = cp[1] * cp[1] + b01[0];
        g[21] = cp[4] * cp[4] + b01[1];
        g[18] = cp[1] * c0[1] + b00[0];
        g[19] = cp[4] * c0[4] + b00[1];
        g[22] = cp[1] * (g[18] + b00[0]) + b01[0] * c0[1];
        g[23] = cp[4] * (g[19] + b00[1]) + b01[1] * c0[4];
        //g[24] = w[0];
        //g[25] = w[1];
        g[26] = c0[2] * g[24];
        g[27] = c0[5] * g[25];
        g[28] = cp[2] * g[24];
        g[29] = cp[5] * g[25];
        g[32] = cp[2] * g[28] + b01[0] * g[24];
        g[33] = cp[5] * g[29] + b01[1] * g[25];
        g[30] = cp[2] * g[26] + b00[0] * g[24];
        g[31] = cp[5] * g[27] + b00[1] * g[25];
        g[34] = cp[2] * g[30] + b01[0] * g[26] + b00[0] * g[28];
        g[35] = cp[5] * g[31] + b01[1] * g[27] + b00[1] * g[29];
}

static inline void _g0_2d4d_1010(double *restrict g, struct _BC *bc, CINTEnvVars *envs)
{
        double *c0 = bc->c00;
        double *cp = bc->c0p;
        double *b00 = bc->b00;
        g[0] = 1;
        g[1] = 1;
        g[2] = c0[0];
        g[3] = c0[3];
        g[4] = cp[0];
        g[5] = cp[3];
        g[6] = cp[0] * c0[0] + b00[0];
        g[7] = cp[3] * c0[3] + b00[1];
        g[8] = 1;
        g[9] = 1;
        g[10] = c0[1];
        g[11] = c0[4];
        g[12] = cp[1];
        g[13] = cp[4];
        g[14] = cp[1] * c0[1] + b00[0];
        g[15] = cp[4] * c0[4] + b00[1];
        //g[16] = w[0];
        //g[17] = w[1];
        g[18] = c0[2] * g[16];
        g[19] = c0[5] * g[17];
        g[20] = cp[2] * g[16];
        g[21] = cp[5] * g[17];
        g[22] = cp[2] * g[18] + b00[0] * g[16];
        g[23] = cp[5] * g[19] + b00[1] * g[17];
}

static inline void _g0_2d4d_1011(double *restrict g, struct _BC *bc, CINTEnvVars *envs)
{
        double *c0 = bc->c00;
        double *cp = bc->c0p;
        double *b00 = bc->b00;
        double *b01 = bc->b01;
        double xkxl = envs->rkrl[0];
        double ykyl = envs->rkrl[1];
        double zkzl = envs->rkrl[2];
        g[0] = 1;
        g[1] = 1;
        g[2] = c0[0];
        g[3] = c0[3];
        g[8] = cp[0];
        g[9] = cp[3];
        g[10] = cp[0] * c0[0] + b00[0];
        g[11] = cp[3] * c0[3] + b00[1];
        g[12] = cp[0] * (xkxl + cp[0]) + b01[0];
        g[13] = cp[3] * (xkxl + cp[3]) + b01[1];
        g[14] = g[10] * (xkxl + cp[0]) + cp[0] * b00[0] + b01[0] * c0[0];
        g[15] = g[11] * (xkxl + cp[3]) + cp[3] * b00[1] + b01[1] * c0[3];
        g[4] = xkxl + cp[0];
        g[5] = xkxl + cp[3];
        g[6] = c0[0] * (xkxl + cp[0]) + b00[0];
        g[7] = c0[3] * (xkxl + cp[3]) + b00[1];
        g[24] = 1;
        g[25] = 1;
        g[26] = c0[1];
        g[27] = c0[4];
        g[32] = cp[1];
        g[33] = cp[4];
        g[34] = cp[1] * c0[1] + b00[0];
        g[35] = cp[4] * c0[4] + b00[1];
        g[36] = cp[1] * (ykyl + cp[1]) + b01[0];
        g[37] = cp[4] * (ykyl + cp[4]) + b01[1];
        g[38] = g[34] * (ykyl + cp[1]) + cp[1] * b00[0] + b01[0] * c0[1];
        g[39] = g[35] * (ykyl + cp[4]) + cp[4] * b00[1] + b01[1] * c0[4];
        g[28] = ykyl + cp[1];
        g[29] = ykyl + cp[4];
        g[30] = c0[1] * (ykyl + cp[1]) + b00[0];
        g[31] = c0[4] * (ykyl + cp[4]) + b00[1];
        //g[48] = w[0];
        //g[49] = w[1];
        g[50] = c0[2] * g[48];
        g[51] = c0[5] * g[49];
        g[56] = cp[2] * g[48];
        g[57] = cp[5] * g[49];
        g[58] = cp[2] * g[50] + b00[0] * g[48];
        g[59] = cp[5] * g[51] + b00[1] * g[49];
        g[60] = g[56] * (zkzl + cp[2]) + b01[0] * g[48];
        g[61] = g[57] * (zkzl + cp[5]) + b01[1] * g[49];
        g[62] = g[58] * (zkzl + cp[2]) + b01[0] * g[50] + b00[0] * g[56];
        g[63] = g[59] * (zkzl + cp[5]) + b01[1] * g[51] + b00[1] * g[57];
        g[52] = g[48] * (zkzl + cp[2]);
        g[53] = g[49] * (zkzl + cp[5]);
        g[54] = g[50] * (zkzl + cp[2]) + b00[0] * g[48];
        g[55] = g[51] * (zkzl + cp[5]) + b00[1] * g[49];
}

static inline void _g0_2d4d_1020(double *restrict g, struct _BC *bc, CINTEnvVars *envs)
{
        double *c0 = bc->c00;
        double *cp = bc->c0p;
        double *b00 = bc->b00;
        double *b01 = bc->b01;
        g[0] = 1;
        g[1] = 1;
        g[2] = c0[0];
        g[3] = c0[3];
        g[4] = cp[0];
        g[5] = cp[3];
        g[8] = cp[0] * cp[0] + b01[0];
        g[9] = cp[3] * cp[3] + b01[1];
        g[6] = cp[0] * c0[0] + b00[0];
        g[7] = cp[3] * c0[3] + b00[1];
        g[10] = cp[0] * (g[6] + b00[0]) + b01[0] * c0[0];
        g[11] = cp[3] * (g[7] + b00[1]) + b01[1] * c0[3];
        g[12] = 1;
        g[13] = 1;
        g[14] = c0[1];
        g[15] = c0[4];
        g[16] = cp[1];
        g[17] = cp[4];
        g[20] = cp[1] * cp[1] + b01[0];
        g[21] = cp[4] * cp[4] + b01[1];
        g[18] = cp[1] * c0[1] + b00[0];
        g[19] = cp[4] * c0[4] + b00[1];
        g[22] = cp[1] * (g[18] + b00[0]) + b01[0] * c0[1];
        g[23] = cp[4] * (g[19] + b00[1]) + b01[1] * c0[4];
        //g[24] = w[0];
        //g[25] = w[1];
        g[26] = c0[2] * g[24];
        g[27] = c0[5] * g[25];
        g[28] = cp[2] * g[24];
        g[29] = cp[5] * g[25];
        g[32] = cp[2] * g[28] + b01[0] * g[24];
        g[33] = cp[5] * g[29] + b01[1] * g[25];
        g[30] = cp[2] * g[26] + b00[0] * g[24];
        g[31] = cp[5] * g[27] + b00[1] * g[25];
        g[34] = cp[2] * g[30] + b01[0] * g[26] + b00[0] * g[28];
        g[35] = cp[5] * g[31] + b01[1] * g[27] + b00[1] * g[29];
}

static inline void _g0_2d4d_1100(double *restrict g, struct _BC *bc, CINTEnvVars *envs)
{
        double *c0 = bc->c00;
        double *b10 = bc->b10;
        double xixj = envs->rirj[0];
        double yiyj = envs->rirj[1];
        double zizj = envs->rirj[2];
        g[0] = 1;
        g[1] = 1;
        g[4] = c0[0];
        g[5] = c0[3];
        g[6] = c0[0] * (xixj + c0[0]) + b10[0];
        g[7] = c0[3] * (xixj + c0[3]) + b10[1];
        g[2] = xixj + c0[0];
        g[3] = xixj + c0[3];
        g[12] = 1;
        g[13] = 1;
        g[16] = c0[1];
        g[17] = c0[4];
        g[18] = c0[1] * (yiyj + c0[1]) + b10[0];
        g[19] = c0[4] * (yiyj + c0[4]) + b10[1];
        g[14] = yiyj + c0[1];
        g[15] = yiyj + c0[4];
        //g[24] = w[0];
        //g[25] = w[1];
        g[28] = c0[2] * g[24];
        g[29] = c0[5] * g[25];
        g[30] = g[28] * (zizj + c0[2]) + b10[0] * g[24];
        g[31] = g[29] * (zizj + c0[5]) + b10[1] * g[25];
        g[26] = g[24] * (zizj + c0[2]);
        g[27] = g[25] * (zizj + c0[5]);
}

static inline void _g0_2d4d_1101(double *restrict g, struct _BC *bc, CINTEnvVars *envs)
{
        double *c0 = bc->c00;
        double *cp = bc->c0p;
        double *b00 = bc->b00;
        double *b10 = bc->b10;
        double xixj = envs->rirj[0];
        double yiyj = envs->rirj[1];
        double zizj = envs->rirj[2];
        g[0] = 1;
        g[1] = 1;
        g[8] = c0[0];
        g[9] = c0[3];
        g[4] = cp[0];
        g[5] = cp[3];
        g[12] = cp[0] * c0[0] + b00[0];
        g[13] = cp[3] * c0[3] + b00[1];
        g[10] = c0[0] * (xixj + c0[0]) + b10[0];
        g[11] = c0[3] * (xixj + c0[3]) + b10[1];
        g[2] = xixj + c0[0];
        g[3] = xixj + c0[3];
        g[14] = g[12] * (xixj + c0[0]) + c0[0] * b00[0] + b10[0] * cp[0];
        g[15] = g[13] * (xixj + c0[3]) + c0[3] * b00[1] + b10[1] * cp[3];
        g[6] = cp[0] * (xixj + c0[0]) + b00[0];
        g[7] = cp[3] * (xixj + c0[3]) + b00[1];
        g[24] = 1;
        g[25] = 1;
        g[32] = c0[1];
        g[33] = c0[4];
        g[28] = cp[1];
        g[29] = cp[4];
        g[36] = cp[1] * c0[1] + b00[0];
        g[37] = cp[4] * c0[4] + b00[1];
        g[34] = c0[1] * (yiyj + c0[1]) + b10[0];
        g[35] = c0[4] * (yiyj + c0[4]) + b10[1];
        g[26] = yiyj + c0[1];
        g[27] = yiyj + c0[4];
        g[38] = g[36] * (yiyj + c0[1]) + c0[1] * b00[0] + b10[0] * cp[1];
        g[39] = g[37] * (yiyj + c0[4]) + c0[4] * b00[1] + b10[1] * cp[4];
        g[30] = cp[1] * (yiyj + c0[1]) + b00[0];
        g[31] = cp[4] * (yiyj + c0[4]) + b00[1];
        //g[48] = w[0];
        //g[49] = w[1];
        g[56] = c0[2] * g[48];
        g[57] = c0[5] * g[49];
        g[52] = cp[2] * g[48];
        g[53] = cp[5] * g[49];
        g[60] = cp[2] * g[56] + b00[0] * g[48];
        g[61] = cp[5] * g[57] + b00[1] * g[49];
        g[58] = g[56] * (zizj + c0[2]) + b10[0] * g[48];
        g[59] = g[57] * (zizj + c0[5]) + b10[1] * g[49];
        g[50] = g[48] * (zizj + c0[2]);
        g[51] = g[49] * (zizj + c0[5]);
        g[62] = g[60] * (zizj + c0[2]) + b10[0] * g[52] + b00[0] * g[56];
        g[63] = g[61] * (zizj + c0[5]) + b10[1] * g[53] + b00[1] * g[57];
        g[54] = zizj * g[52] + cp[2] * g[56] + b00[0] * g[48];
        g[55] = zizj * g[53] + cp[5] * g[57] + b00[1] * g[49];
}

static inline void _g0_2d4d_1110(double *restrict g, struct _BC *bc, CINTEnvVars *envs)
{
        double *c0 = bc->c00;
        double *cp = bc->c0p;
        double *b00 = bc->b00;
        double *b10 = bc->b10;
        double xixj = envs->rirj[0];
        double yiyj = envs->rirj[1];
        double zizj = envs->rirj[2];
        g[0] = 1;
        g[1] = 1;
        g[8] = c0[0];
        g[9] = c0[3];
        g[4] = cp[0];
        g[5] = cp[3];
        g[12] = cp[0] * c0[0] + b00[0];
        g[13] = cp[3] * c0[3] + b00[1];
        g[10] = c0[0] * (xixj + c0[0]) + b10[0];
        g[11] = c0[3] * (xixj + c0[3]) + b10[1];
        g[2] = xixj + c0[0];
        g[3] = xixj + c0[3];
        g[14] = g[12] * (xixj + c0[0]) + c0[0] * b00[0] + b10[0] * cp[0];
        g[15] = g[13] * (xixj + c0[3]) + c0[3] * b00[1] + b10[1] * cp[3];
        g[6] = cp[0] * (xixj + c0[0]) + b00[0];
        g[7] = cp[3] * (xixj + c0[3]) + b00[1];
        g[24] = 1;
        g[25] = 1;
        g[32] = c0[1];
        g[33] = c0[4];
        g[28] = cp[1];
        g[29] = cp[4];
        g[36] = cp[1] * c0[1] + b00[0];
        g[37] = cp[4] * c0[4] + b00[1];
        g[34] = c0[1] * (yiyj + c0[1]) + b10[0];
        g[35] = c0[4] * (yiyj + c0[4]) + b10[1];
        g[26] = yiyj + c0[1];
        g[27] = yiyj + c0[4];
        g[38] = g[36] * (yiyj + c0[1]) + c0[1] * b00[0] + b10[0] * cp[1];
        g[39] = g[37] * (yiyj + c0[4]) + c0[4] * b00[1] + b10[1] * cp[4];
        g[30] = cp[1] * (yiyj + c0[1]) + b00[0];
        g[31] = cp[4] * (yiyj + c0[4]) + b00[1];
        //g[48] = w[0];
        //g[49] = w[1];
        g[56] = c0[2] * g[48];
        g[57] = c0[5] * g[49];
        g[52] = cp[2] * g[48];
        g[53] = cp[5] * g[49];
        g[60] = cp[2] * g[56] + b00[0] * g[48];
        g[61] = cp[5] * g[57] + b00[1] * g[49];
        g[58] = g[56] * (zizj + c0[2]) + b10[0] * g[48];
        g[59] = g[57] * (zizj + c0[5]) + b10[1] * g[49];
        g[50] = g[48] * (zizj + c0[2]);
        g[51] = g[49] * (zizj + c0[5]);
        g[62] = g[60] * (zizj + c0[2]) + b10[0] * g[52] + b00[0] * g[56];
        g[63] = g[61] * (zizj + c0[5]) + b10[1] * g[53] + b00[1] * g[57];
        g[54] = zizj * g[52] + cp[2] * g[56] + b00[0] * g[48];
        g[55] = zizj * g[53] + cp[5] * g[57] + b00[1] * g[49];
}

static inline void _g0_2d4d_1200(double *restrict g, struct _BC *bc, CINTEnvVars *envs)
{
        double *c0 = bc->c00;
        double *b10 = bc->b10;
        double xixj = envs->rirj[0];
        double yiyj = envs->rirj[1];
        double zizj = envs->rirj[2];
        g[0] = 1;
        g[1] = 1;
        g[4] = c0[0];
        g[5] = c0[3];
        g[8] = c0[0] * c0[0] + b10[0];
        g[9] = c0[3] * c0[3] + b10[1];
        g[10] = g[8] * (xixj + c0[0]) + c0[0] * 2 * b10[0];
        g[11] = g[9] * (xixj + c0[3]) + c0[3] * 2 * b10[1];
        g[6] = c0[0] * (xixj + c0[0]) + b10[0];
        g[7] = c0[3] * (xixj + c0[3]) + b10[1];
        g[2] = xixj + c0[0];
        g[3] = xixj + c0[3];
        g[16] = 1;
        g[17] = 1;
        g[20] = c0[1];
        g[21] = c0[4];
        g[24] = c0[1] * c0[1] + b10[0];
        g[25] = c0[4] * c0[4] + b10[1];
        g[26] = g[24] * (yiyj + c0[1]) + c0[1] * 2 * b10[0];
        g[27] = g[25] * (yiyj + c0[4]) + c0[4] * 2 * b10[1];
        g[22] = c0[1] * (yiyj + c0[1]) + b10[0];
        g[23] = c0[4] * (yiyj + c0[4]) + b10[1];
        g[18] = yiyj + c0[1];
        g[19] = yiyj + c0[4];
        //g[32] = w[0];
        //g[33] = w[1];
        g[36] = c0[2] * g[32];
        g[37] = c0[5] * g[33];
        g[40] = c0[2] * g[36] + b10[0] * g[32];
        g[41] = c0[5] * g[37] + b10[1] * g[33];
        g[42] = g[40] * (zizj + c0[2]) + 2 * b10[0] * g[36];
        g[43] = g[41] * (zizj + c0[5]) + 2 * b10[1] * g[37];
        g[38] = g[36] * (zizj + c0[2]) + b10[0] * g[32];
        g[39] = g[37] * (zizj + c0[5]) + b10[1] * g[33];
        g[34] = g[32] * (zizj + c0[2]);
        g[35] = g[33] * (zizj + c0[5]);
}

static inline void _g0_2d4d_2000(double *restrict g, struct _BC *bc, CINTEnvVars *envs)
{
        double *c0 = bc->c00;
        double *b10 = bc->b10;
        g[0] = 1;
        g[1] = 1;
        g[2] = c0[0];
        g[3] = c0[3];
        g[4] = c0[0] * c0[0] + b10[0];
        g[5] = c0[3] * c0[3] + b10[1];
        g[6] = 1;
        g[7] = 1;
        g[8] = c0[1];
        g[9] = c0[4];
        g[10] = c0[1] * c0[1] + b10[0];
        g[11] = c0[4] * c0[4] + b10[1];
        //g[12] = w[0];
        //g[13] = w[1];
        g[14] = c0[2] * g[12];
        g[15] = c0[5] * g[13];
        g[16] = c0[2] * g[14] + b10[0] * g[12];
        g[17] = c0[5] * g[15] + b10[1] * g[13];
}

static inline void _g0_2d4d_2001(double *restrict g, struct _BC *bc, CINTEnvVars *envs)
{
        double *c0 = bc->c00;
        double *cp = bc->c0p;
        double *b00 = bc->b00;
        double *b10 = bc->b10;
        g[0] = 1;
        g[1] = 1;
        g[2] = c0[0];
        g[3] = c0[3];
        g[4] = c0[0] * c0[0] + b10[0];
        g[5] = c0[3] * c0[3] + b10[1];
        g[6] = cp[0];
        g[7] = cp[3];
        g[8] = cp[0] * c0[0] + b00[0];
        g[9] = cp[3] * c0[3] + b00[1];
        g[10] = c0[0] * (g[8] + b00[0]) + b10[0] * cp[0];
        g[11] = c0[3] * (g[9] + b00[1]) + b10[1] * cp[3];
        g[12] = 1;
        g[13] = 1;
        g[14] = c0[1];
        g[15] = c0[4];
        g[16] = c0[1] * c0[1] + b10[0];
        g[17] = c0[4] * c0[4] + b10[1];
        g[18] = cp[1];
        g[19] = cp[4];
        g[20] = cp[1] * c0[1] + b00[0];
        g[21] = cp[4] * c0[4] + b00[1];
        g[22] = c0[1] * (g[20] + b00[0]) + b10[0] * cp[1];
        g[23] = c0[4] * (g[21] + b00[1]) + b10[1] * cp[4];
        //g[24] = w[0];
        //g[25] = w[1];
        g[26] = c0[2] * g[24];
        g[27] = c0[5] * g[25];
        g[28] = c0[2] * g[26] + b10[0] * g[24];
        g[29] = c0[5] * g[27] + b10[1] * g[25];
        g[30] = cp[2] * g[24];
        g[31] = cp[5] * g[25];
        g[32] = cp[2] * g[26] + b00[0] * g[24];
        g[33] = cp[5] * g[27] + b00[1] * g[25];
        g[34] = c0[2] * g[32] + b10[0] * g[30] + b00[0] * g[26];
        g[35] = c0[5] * g[33] + b10[1] * g[31] + b00[1] * g[27];
}

static inline void _g0_2d4d_2010(double *restrict g, struct _BC *bc, CINTEnvVars *envs)
{
        double *c0 = bc->c00;
        double *cp = bc->c0p;
        double *b00 = bc->b00;
        double *b10 = bc->b10;
        g[0] = 1;
        g[1] = 1;
        g[2] = c0[0];
        g[3] = c0[3];
        g[4] = c0[0] * c0[0] + b10[0];
        g[5] = c0[3] * c0[3] + b10[1];
        g[6] = cp[0];
        g[7] = cp[3];
        g[8] = cp[0] * c0[0] + b00[0];
        g[9] = cp[3] * c0[3] + b00[1];
        g[10] = c0[0] * (g[8] + b00[0]) + b10[0] * cp[0];
        g[11] = c0[3] * (g[9] + b00[1]) + b10[1] * cp[3];
        g[12] = 1;
        g[13] = 1;
        g[14] = c0[1];
        g[15] = c0[4];
        g[16] = c0[1] * c0[1] + b10[0];
        g[17] = c0[4] * c0[4] + b10[1];
        g[18] = cp[1];
        g[19] = cp[4];
        g[20] = cp[1] * c0[1] + b00[0];
        g[21] = cp[4] * c0[4] + b00[1];
        g[22] = c0[1] * (g[20] + b00[0]) + b10[0] * cp[1];
        g[23] = c0[4] * (g[21] + b00[1]) + b10[1] * cp[4];
        //g[24] = w[0];
        //g[25] = w[1];
        g[26] = c0[2] * g[24];
        g[27] = c0[5] * g[25];
        g[28] = c0[2] * g[26] + b10[0] * g[24];
        g[29] = c0[5] * g[27] + b10[1] * g[25];
        g[30] = cp[2] * g[24];
        g[31] = cp[5] * g[25];
        g[32] = cp[2] * g[26] + b00[0] * g[24];
        g[33] = cp[5] * g[27] + b00[1] * g[25];
        g[34] = c0[2] * g[32] + b10[0] * g[30] + b00[0] * g[26];
        g[35] = c0[5] * g[33] + b10[1] * g[31] + b00[1] * g[27];
}

static inline void _g0_2d4d_2100(double *restrict g, struct _BC *bc, CINTEnvVars *envs)
{
        double *c0 = bc->c00;
        double *b10 = bc->b10;
        double xixj = envs->rirj[0];
        double yiyj = envs->rirj[1];
        double zizj = envs->rirj[2];
        g[0] = 1;
        g[1] = 1;
        g[2] = c0[0];
        g[3] = c0[3];
        g[4] = c0[0] * c0[0] + b10[0];
        g[5] = c0[3] * c0[3] + b10[1];
        g[12] = g[4] * (xixj + c0[0]) + c0[0] * 2 * b10[0];
        g[13] = g[5] * (xixj + c0[3]) + c0[3] * 2 * b10[1];
        g[10] = c0[0] * (xixj + c0[0]) + b10[0];
        g[11] = c0[3] * (xixj + c0[3]) + b10[1];
        g[8] = xixj + c0[0];
        g[9] = xixj + c0[3];
        g[16] = 1;
        g[17] = 1;
        g[18] = c0[1];
        g[19] = c0[4];
        g[20] = c0[1] * c0[1] + b10[0];
        g[21] = c0[4] * c0[4] + b10[1];
        g[28] = g[20] * (yiyj + c0[1]) + c0[1] * 2 * b10[0];
        g[29] = g[21] * (yiyj + c0[4]) + c0[4] * 2 * b10[1];
        g[26] = c0[1] * (yiyj + c0[1]) + b10[0];
        g[27] = c0[4] * (yiyj + c0[4]) + b10[1];
        g[24] = yiyj + c0[1];
        g[25] = yiyj + c0[4];
        //g[32] = w[0];
        //g[33] = w[1];
        g[34] = c0[2] * g[32];
        g[35] = c0[5] * g[33];
        g[36] = c0[2] * g[34] + b10[0] * g[32];
        g[37] = c0[5] * g[35] + b10[1] * g[33];
        g[44] = g[36] * (zizj + c0[2]) + 2 * b10[0] * g[34];
        g[45] = g[37] * (zizj + c0[5]) + 2 * b10[1] * g[35];
        g[42] = g[34] * (zizj + c0[2]) + b10[0] * g[32];
        g[43] = g[35] * (zizj + c0[5]) + b10[1] * g[33];
        g[40] = g[32] * (zizj + c0[2]);
        g[41] = g[33] * (zizj + c0[5]);
}

static inline void _g0_2d4d_3000(double *restrict g, struct _BC *bc, CINTEnvVars *envs)
{
        double *c0 = bc->c00;
        double *b10 = bc->b10;
        g[0] = 1;
        g[1] = 1;
        g[2] = c0[0];
        g[3] = c0[3];
        g[4] = c0[0] * c0[0] + b10[0];
        g[5] = c0[3] * c0[3] + b10[1];
        g[6] = c0[0] * (g[4] + 2 * b10[0]);
        g[7] = c0[3] * (g[5] + 2 * b10[1]);
        g[8] = 1;
        g[9] = 1;
        g[10] = c0[1];
        g[11] = c0[4];
        g[12] = c0[1] * c0[1] + b10[0];
        g[13] = c0[4] * c0[4] + b10[1];
        g[14] = c0[1] * (g[12] + 2 * b10[0]);
        g[15] = c0[4] * (g[13] + 2 * b10[1]);
        //g[16] = w[0];
        //g[17] = w[1];
        g[18] = c0[2] * g[16];
        g[19] = c0[5] * g[17];
        g[20] = c0[2] * g[18] + b10[0] * g[16];
        g[21] = c0[5] * g[19] + b10[1] * g[17];
        g[22] = c0[2] * g[20] + 2 * b10[0] * g[18];
        g[23] = c0[5] * g[21] + 2 * b10[1] * g[19];
}

void CINTg0_2e_2d4d_unrolled(double *restrict g, struct _BC *bc, CINTEnvVars *envs)
{
        int type_ijkl = ((envs->li_ceil << 6) | (envs->lj_ceil << 4) |
                         (envs->lk_ceil << 2) | (envs->ll_ceil));
        switch (type_ijkl) {
                case 0b00000000: _g0_2d4d_0000(g, bc, envs); return;
                case 0b00000001: _g0_2d4d_0001(g, bc, envs); return;
                case 0b00000010: _g0_2d4d_0002(g, bc, envs); return;
                case 0b00000011: _g0_2d4d_0003(g, bc, envs); return;
                case 0b00000100: _g0_2d4d_0010(g, bc, envs); return;
                case 0b00000101: _g0_2d4d_0011(g, bc, envs); return;
                case 0b00000110: _g0_2d4d_0012(g, bc, envs); return;
                case 0b00001000: _g0_2d4d_0020(g, bc, envs); return;
                case 0b00001001: _g0_2d4d_0021(g, bc, envs); return;
                case 0b00001100: _g0_2d4d_0030(g, bc, envs); return;
                case 0b00010000: _g0_2d4d_0100(g, bc, envs); return;
                case 0b00010001: _g0_2d4d_0101(g, bc, envs); return;
                case 0b00010010: _g0_2d4d_0102(g, bc, envs); return;
                case 0b00010100: _g0_2d4d_0110(g, bc, envs); return;
                case 0b00010101: _g0_2d4d_0111(g, bc, envs); return;
                case 0b00011000: _g0_2d4d_0120(g, bc, envs); return;
                case 0b00100000: _g0_2d4d_0200(g, bc, envs); return;
                case 0b00100001: _g0_2d4d_0201(g, bc, envs); return;
                case 0b00100100: _g0_2d4d_0210(g, bc, envs); return;
                case 0b00110000: _g0_2d4d_0300(g, bc, envs); return;
                case 0b01000000: _g0_2d4d_1000(g, bc, envs); return;
                case 0b01000001: _g0_2d4d_1001(g, bc, envs); return;
                case 0b01000010: _g0_2d4d_1002(g, bc, envs); return;
                case 0b01000100: _g0_2d4d_1010(g, bc, envs); return;
                case 0b01000101: _g0_2d4d_1011(g, bc, envs); return;
                case 0b01001000: _g0_2d4d_1020(g, bc, envs); return;
                case 0b01010000: _g0_2d4d_1100(g, bc, envs); return;
                case 0b01010001: _g0_2d4d_1101(g, bc, envs); return;
                case 0b01010100: _g0_2d4d_1110(g, bc, envs); return;
                case 0b01100000: _g0_2d4d_1200(g, bc, envs); return;
                case 0b10000000: _g0_2d4d_2000(g, bc, envs); return;
                case 0b10000001: _g0_2d4d_2001(g, bc, envs); return;
                case 0b10000100: _g0_2d4d_2010(g, bc, envs); return;
                case 0b10010000: _g0_2d4d_2100(g, bc, envs); return;
                case 0b11000000: _g0_2d4d_3000(g, bc, envs); return;
        }
        fprintf(stderr, "Dimension error for CINTg0_2e_lj2d4d: iklj = %d %d %d %d",
               (int)envs->li_ceil, (int)envs->lk_ceil,
               (int)envs->ll_ceil, (int)envs->lj_ceil);
}

static inline void _srg0_2d4d_0000(double *restrict g, struct _BC *bc, CINTEnvVars *envs)
{
        g[0] = 1;
        g[1] = 1;
        g[2] = 1;
        g[3] = 1;
        //g[4] = w[0];
        //g[5] = w[1];
}

static inline void _srg0_2d4d_0001(double *restrict g, struct _BC *bc, CINTEnvVars *envs)
{
        double *cp = bc->c0p;
        g[0] = 1;
        g[1] = 1;
        g[2] = cp[0];
        g[3] = cp[3];
        g[4] = 1;
        g[5] = 1;
        g[6] = cp[1];
        g[7] = cp[4];
        //g[8] = w[0];
        //g[9] = w[1];
        g[10] = cp[2] * g[8];
        g[11] = cp[5] * g[9];
}

static inline void _srg0_2d4d_0002(double *restrict g, struct _BC *bc, CINTEnvVars *envs)
{
        double *cp = bc->c0p;
        double *b01 = bc->b01;
        g[0] = 1;
        g[1] = 1;
        g[2] = 1;
        g[3] = 1;
        g[4] = cp[0];
        g[5] = cp[3];
        g[6] = cp[6];
        g[7] = cp[9];
        g[8] = cp[0] * cp[0] + b01[0];
        g[9] = cp[3] * cp[3] + b01[1];
        g[10] = cp[6] * cp[6] + b01[2];
        g[11] = cp[9] * cp[9] + b01[3];
        g[12] = 1;
        g[13] = 1;
        g[14] = 1;
        g[15] = 1;
        g[16] = cp[1];
        g[17] = cp[4];
        g[18] = cp[7];
        g[19] = cp[10];
        g[20] = cp[1] * cp[1] + b01[0];
        g[21] = cp[4] * cp[4] + b01[1];
        g[22] = cp[7] * cp[7] + b01[2];
        g[23] = cp[10] * cp[10] + b01[3];
        //g[24] = w[0];
        //g[25] = w[1];
        //g[26] = w[2];
        //g[27] = w[3];
        g[28] = cp[2] * g[24];
        g[29] = cp[5] * g[25];
        g[30] = cp[8] * g[26];
        g[31] = cp[11] * g[27];
        g[32] = cp[2] * g[28] + b01[0] * g[24];
        g[33] = cp[5] * g[29] + b01[1] * g[25];
        g[34] = cp[8] * g[30] + b01[2] * g[26];
        g[35] = cp[11] * g[31] + b01[3] * g[27];
}

static inline void _srg0_2d4d_0003(double *restrict g, struct _BC *bc, CINTEnvVars *envs)
{
        double *cp = bc->c0p;
        double *b01 = bc->b01;
        g[0] = 1;
        g[1] = 1;
        g[2] = 1;
        g[3] = 1;
        g[4] = cp[0];
        g[5] = cp[3];
        g[6] = cp[6];
        g[7] = cp[9];
        g[8] = cp[0] * cp[0] + b01[0];
        g[9] = cp[3] * cp[3] + b01[1];
        g[10] = cp[6] * cp[6] + b01[2];
        g[11] = cp[9] * cp[9] + b01[3];
        g[12] = cp[0] * (g[8] + 2 * b01[0]);
        g[13] = cp[3] * (g[9] + 2 * b01[1]);
        g[14] = cp[6] * (g[10] + 2 * b01[2]);
        g[15] = cp[9] * (g[11] + 2 * b01[3]);
        g[16] = 1;
        g[17] = 1;
        g[18] = 1;
        g[19] = 1;
        g[20] = cp[1];
        g[21] = cp[4];
        g[22] = cp[7];
        g[23] = cp[10];
        g[24] = cp[1] * cp[1] + b01[0];
        g[25] = cp[4] * cp[4] + b01[1];
        g[26] = cp[7] * cp[7] + b01[2];
        g[27] = cp[10] * cp[10] + b01[3];
        g[28] = cp[1] * (g[24] + 2 * b01[0]);
        g[29] = cp[4] * (g[25] + 2 * b01[1]);
        g[30] = cp[7] * (g[26] + 2 * b01[2]);
        g[31] = cp[10] * (g[27] + 2 * b01[3]);
        //g[32] = w[0];
        //g[33] = w[1];
        //g[34] = w[2];
        //g[35] = w[3];
        g[36] = cp[2] * g[32];
        g[37] = cp[5] * g[33];
        g[38] = cp[8] * g[34];
        g[39] = cp[11] * g[35];
        g[40] = cp[2] * g[36] + b01[0] * g[32];
        g[41] = cp[5] * g[37] + b01[1] * g[33];
        g[42] = cp[8] * g[38] + b01[2] * g[34];
        g[43] = cp[11] * g[39] + b01[3] * g[35];
        g[44] = cp[2] * g[40] + 2 * b01[0] * g[36];
        g[45] = cp[5] * g[41] + 2 * b01[1] * g[37];
        g[46] = cp[8] * g[42] + 2 * b01[2] * g[38];
        g[47] = cp[11] * g[43] + 2 * b01[3] * g[39];
}

static inline void _srg0_2d4d_0010(double *restrict g, struct _BC *bc, CINTEnvVars *envs)
{
        double *cp = bc->c0p;
        g[0] = 1;
        g[1] = 1;
        g[2] = cp[0];
        g[3] = cp[3];
        g[4] = 1;
        g[5] = 1;
        g[6] = cp[1];
        g[7] = cp[4];
        //g[8] = w[0];
        //g[9] = w[1];
        g[10] = cp[2] * g[8];
        g[11] = cp[5] * g[9];
}

static inline void _srg0_2d4d_0011(double *restrict g, struct _BC *bc, CINTEnvVars *envs)
{
        double *cp = bc->c0p;
        double *b01 = bc->b01;
        double xkxl = envs->rkrl[0];
        double ykyl = envs->rkrl[1];
        double zkzl = envs->rkrl[2];
        g[0] = 1;
        g[1] = 1;
        g[2] = 1;
        g[3] = 1;
        g[8] = cp[0];
        g[9] = cp[3];
        g[10] = cp[6];
        g[11] = cp[9];
        g[12] = cp[0] * (xkxl + cp[0]) + b01[0];
        g[13] = cp[3] * (xkxl + cp[3]) + b01[1];
        g[14] = cp[6] * (xkxl + cp[6]) + b01[2];
        g[15] = cp[9] * (xkxl + cp[9]) + b01[3];
        g[4] = xkxl + cp[0];
        g[5] = xkxl + cp[3];
        g[6] = xkxl + cp[6];
        g[7] = xkxl + cp[9];
        g[24] = 1;
        g[25] = 1;
        g[26] = 1;
        g[27] = 1;
        g[32] = cp[1];
        g[33] = cp[4];
        g[34] = cp[7];
        g[35] = cp[10];
        g[36] = cp[1] * (ykyl + cp[1]) + b01[0];
        g[37] = cp[4] * (ykyl + cp[4]) + b01[1];
        g[38] = cp[7] * (ykyl + cp[7]) + b01[2];
        g[39] = cp[10] * (ykyl + cp[10]) + b01[3];
        g[28] = ykyl + cp[1];
        g[29] = ykyl + cp[4];
        g[30] = ykyl + cp[7];
        g[31] = ykyl + cp[10];
        //g[48] = w[0];
        //g[49] = w[1];
        //g[50] = w[2];
        //g[51] = w[3];
        g[56] = cp[2] * g[48];
        g[57] = cp[5] * g[49];
        g[58] = cp[8] * g[50];
        g[59] = cp[11] * g[51];
        g[60] = g[56] * (zkzl + cp[2]) + b01[0] * g[48];
        g[61] = g[57] * (zkzl + cp[5]) + b01[1] * g[49];
        g[62] = g[58] * (zkzl + cp[8]) + b01[2] * g[50];
        g[63] = g[59] * (zkzl + cp[11]) + b01[3] * g[51];
        g[52] = g[48] * (zkzl + cp[2]);
        g[53] = g[49] * (zkzl + cp[5]);
        g[54] = g[50] * (zkzl + cp[8]);
        g[55] = g[51] * (zkzl + cp[11]);
}

static inline void _srg0_2d4d_0012(double *restrict g, struct _BC *bc, CINTEnvVars *envs)
{
        double *cp = bc->c0p;
        double *b01 = bc->b01;
        double xkxl = envs->rkrl[0];
        double ykyl = envs->rkrl[1];
        double zkzl = envs->rkrl[2];
        g[0] = 1;
        g[1] = 1;
        g[2] = 1;
        g[3] = 1;
        g[8] = cp[0];
        g[9] = cp[3];
        g[10] = cp[6];
        g[11] = cp[9];
        g[16] = cp[0] * cp[0] + b01[0];
        g[17] = cp[3] * cp[3] + b01[1];
        g[18] = cp[6] * cp[6] + b01[2];
        g[19] = cp[9] * cp[9] + b01[3];
        g[20] = g[16] * (xkxl + cp[0]) + cp[0] * 2 * b01[0];
        g[21] = g[17] * (xkxl + cp[3]) + cp[3] * 2 * b01[1];
        g[22] = g[18] * (xkxl + cp[6]) + cp[6] * 2 * b01[2];
        g[23] = g[19] * (xkxl + cp[9]) + cp[9] * 2 * b01[3];
        g[12] = cp[0] * (xkxl + cp[0]) + b01[0];
        g[13] = cp[3] * (xkxl + cp[3]) + b01[1];
        g[14] = cp[6] * (xkxl + cp[6]) + b01[2];
        g[15] = cp[9] * (xkxl + cp[9]) + b01[3];
        g[4] = xkxl + cp[0];
        g[5] = xkxl + cp[3];
        g[6] = xkxl + cp[6];
        g[7] = xkxl + cp[9];
        g[32] = 1;
        g[33] = 1;
        g[34] = 1;
        g[35] = 1;
        g[40] = cp[1];
        g[41] = cp[4];
        g[42] = cp[7];
        g[43] = cp[10];
        g[48] = cp[1] * cp[1] + b01[0];
        g[49] = cp[4] * cp[4] + b01[1];
        g[50] = cp[7] * cp[7] + b01[2];
        g[51] = cp[10] * cp[10] + b01[3];
        g[52] = g[48] * (ykyl + cp[1]) + cp[1] * 2 * b01[0];
        g[53] = g[49] * (ykyl + cp[4]) + cp[4] * 2 * b01[1];
        g[54] = g[50] * (ykyl + cp[7]) + cp[7] * 2 * b01[2];
        g[55] = g[51] * (ykyl + cp[10]) + cp[10] * 2 * b01[3];
        g[44] = cp[1] * (ykyl + cp[1]) + b01[0];
        g[45] = cp[4] * (ykyl + cp[4]) + b01[1];
        g[46] = cp[7] * (ykyl + cp[7]) + b01[2];
        g[47] = cp[10] * (ykyl + cp[10]) + b01[3];
        g[36] = ykyl + cp[1];
        g[37] = ykyl + cp[4];
        g[38] = ykyl + cp[7];
        g[39] = ykyl + cp[10];
        //g[64] = w[0];
        //g[65] = w[1];
        //g[66] = w[2];
        //g[67] = w[3];
        g[72] = cp[2] * g[64];
        g[73] = cp[5] * g[65];
        g[74] = cp[8] * g[66];
        g[75] = cp[11] * g[67];
        g[80] = cp[2] * g[72] + b01[0] * g[64];
        g[81] = cp[5] * g[73] + b01[1] * g[65];
        g[82] = cp[8] * g[74] + b01[2] * g[66];
        g[83] = cp[11] * g[75] + b01[3] * g[67];
        g[84] = g[80] * (zkzl + cp[2]) + 2 * b01[0] * g[72];
        g[85] = g[81] * (zkzl + cp[5]) + 2 * b01[1] * g[73];
        g[86] = g[82] * (zkzl + cp[8]) + 2 * b01[2] * g[74];
        g[87] = g[83] * (zkzl + cp[11]) + 2 * b01[3] * g[75];
        g[76] = g[72] * (zkzl + cp[2]) + b01[0] * g[64];
        g[77] = g[73] * (zkzl + cp[5]) + b01[1] * g[65];
        g[78] = g[74] * (zkzl + cp[8]) + b01[2] * g[66];
        g[79] = g[75] * (zkzl + cp[11]) + b01[3] * g[67];
        g[68] = g[64] * (zkzl + cp[2]);
        g[69] = g[65] * (zkzl + cp[5]);
        g[70] = g[66] * (zkzl + cp[8]);
        g[71] = g[67] * (zkzl + cp[11]);
}

static inline void _srg0_2d4d_0020(double *restrict g, struct _BC *bc, CINTEnvVars *envs)
{
        double *cp = bc->c0p;
        double *b01 = bc->b01;
        g[0] = 1;
        g[1] = 1;
        g[2] = 1;
        g[3] = 1;
        g[4] = cp[0];
        g[5] = cp[3];
        g[6] = cp[6];
        g[7] = cp[9];
        g[8] = cp[0] * cp[0] + b01[0];
        g[9] = cp[3] * cp[3] + b01[1];
        g[10] = cp[6] * cp[6] + b01[2];
        g[11] = cp[9] * cp[9] + b01[3];
        g[12] = 1;
        g[13] = 1;
        g[14] = 1;
        g[15] = 1;
        g[16] = cp[1];
        g[17] = cp[4];
        g[18] = cp[7];
        g[19] = cp[10];
        g[20] = cp[1] * cp[1] + b01[0];
        g[21] = cp[4] * cp[4] + b01[1];
        g[22] = cp[7] * cp[7] + b01[2];
        g[23] = cp[10] * cp[10] + b01[3];
        //g[24] = w[0];
        //g[25] = w[1];
        //g[26] = w[2];
        //g[27] = w[3];
        g[28] = cp[2] * g[24];
        g[29] = cp[5] * g[25];
        g[30] = cp[8] * g[26];
        g[31] = cp[11] * g[27];
        g[32] = cp[2] * g[28] + b01[0] * g[24];
        g[33] = cp[5] * g[29] + b01[1] * g[25];
        g[34] = cp[8] * g[30] + b01[2] * g[26];
        g[35] = cp[11] * g[31] + b01[3] * g[27];
}

static inline void _srg0_2d4d_0021(double *restrict g, struct _BC *bc, CINTEnvVars *envs)
{
        double *cp = bc->c0p;
        double *b01 = bc->b01;
        double xkxl = envs->rkrl[0];
        double ykyl = envs->rkrl[1];
        double zkzl = envs->rkrl[2];
        g[0] = 1;
        g[1] = 1;
        g[2] = 1;
        g[3] = 1;
        g[4] = cp[0];
        g[5] = cp[3];
        g[6] = cp[6];
        g[7] = cp[9];
        g[8] = cp[0] * cp[0] + b01[0];
        g[9] = cp[3] * cp[3] + b01[1];
        g[10] = cp[6] * cp[6] + b01[2];
        g[11] = cp[9] * cp[9] + b01[3];
        g[24] = g[8] * (xkxl + cp[0]) + cp[0] * 2 * b01[0];
        g[25] = g[9] * (xkxl + cp[3]) + cp[3] * 2 * b01[1];
        g[26] = g[10] * (xkxl + cp[6]) + cp[6] * 2 * b01[2];
        g[27] = g[11] * (xkxl + cp[9]) + cp[9] * 2 * b01[3];
        g[20] = cp[0] * (xkxl + cp[0]) + b01[0];
        g[21] = cp[3] * (xkxl + cp[3]) + b01[1];
        g[22] = cp[6] * (xkxl + cp[6]) + b01[2];
        g[23] = cp[9] * (xkxl + cp[9]) + b01[3];
        g[16] = xkxl + cp[0];
        g[17] = xkxl + cp[3];
        g[18] = xkxl + cp[6];
        g[19] = xkxl + cp[9];
        g[32] = 1;
        g[33] = 1;
        g[34] = 1;
        g[35] = 1;
        g[36] = cp[1];
        g[37] = cp[4];
        g[38] = cp[7];
        g[39] = cp[10];
        g[40] = cp[1] * cp[1] + b01[0];
        g[41] = cp[4] * cp[4] + b01[1];
        g[42] = cp[7] * cp[7] + b01[2];
        g[43] = cp[10] * cp[10] + b01[3];
        g[56] = g[40] * (ykyl + cp[1]) + cp[1] * 2 * b01[0];
        g[57] = g[41] * (ykyl + cp[4]) + cp[4] * 2 * b01[1];
        g[58] = g[42] * (ykyl + cp[7]) + cp[7] * 2 * b01[2];
        g[59] = g[43] * (ykyl + cp[10]) + cp[10] * 2 * b01[3];
        g[52] = cp[1] * (ykyl + cp[1]) + b01[0];
        g[53] = cp[4] * (ykyl + cp[4]) + b01[1];
        g[54] = cp[7] * (ykyl + cp[7]) + b01[2];
        g[55] = cp[10] * (ykyl + cp[10]) + b01[3];
        g[48] = ykyl + cp[1];
        g[49] = ykyl + cp[4];
        g[50] = ykyl + cp[7];
        g[51] = ykyl + cp[10];
        //g[64] = w[0];
        //g[65] = w[1];
        //g[66] = w[2];
        //g[67] = w[3];
        g[68] = cp[2] * g[64];
        g[69] = cp[5] * g[65];
        g[70] = cp[8] * g[66];
        g[71] = cp[11] * g[67];
        g[72] = cp[2] * g[68] + b01[0] * g[64];
        g[73] = cp[5] * g[69] + b01[1] * g[65];
        g[74] = cp[8] * g[70] + b01[2] * g[66];
        g[75] = cp[11] * g[71] + b01[3] * g[67];
        g[88] = g[72] * (zkzl + cp[2]) + 2 * b01[0] * g[68];
        g[89] = g[73] * (zkzl + cp[5]) + 2 * b01[1] * g[69];
        g[90] = g[74] * (zkzl + cp[8]) + 2 * b01[2] * g[70];
        g[91] = g[75] * (zkzl + cp[11]) + 2 * b01[3] * g[71];
        g[84] = g[68] * (zkzl + cp[2]) + b01[0] * g[64];
        g[85] = g[69] * (zkzl + cp[5]) + b01[1] * g[65];
        g[86] = g[70] * (zkzl + cp[8]) + b01[2] * g[66];
        g[87] = g[71] * (zkzl + cp[11]) + b01[3] * g[67];
        g[80] = g[64] * (zkzl + cp[2]);
        g[81] = g[65] * (zkzl + cp[5]);
        g[82] = g[66] * (zkzl + cp[8]);
        g[83] = g[67] * (zkzl + cp[11]);
}

static inline void _srg0_2d4d_0030(double *restrict g, struct _BC *bc, CINTEnvVars *envs)
{
        double *cp = bc->c0p;
        double *b01 = bc->b01;
        g[0] = 1;
        g[1] = 1;
        g[2] = 1;
        g[3] = 1;
        g[4] = cp[0];
        g[5] = cp[3];
        g[6] = cp[6];
        g[7] = cp[9];
        g[8] = cp[0] * cp[0] + b01[0];
        g[9] = cp[3] * cp[3] + b01[1];
        g[10] = cp[6] * cp[6] + b01[2];
        g[11] = cp[9] * cp[9] + b01[3];
        g[12] = cp[0] * (g[8] + 2 * b01[0]);
        g[13] = cp[3] * (g[9] + 2 * b01[1]);
        g[14] = cp[6] * (g[10] + 2 * b01[2]);
        g[15] = cp[9] * (g[11] + 2 * b01[3]);
        g[16] = 1;
        g[17] = 1;
        g[18] = 1;
        g[19] = 1;
        g[20] = cp[1];
        g[21] = cp[4];
        g[22] = cp[7];
        g[23] = cp[10];
        g[24] = cp[1] * cp[1] + b01[0];
        g[25] = cp[4] * cp[4] + b01[1];
        g[26] = cp[7] * cp[7] + b01[2];
        g[27] = cp[10] * cp[10] + b01[3];
        g[28] = cp[1] * (g[24] + 2 * b01[0]);
        g[29] = cp[4] * (g[25] + 2 * b01[1]);
        g[30] = cp[7] * (g[26] + 2 * b01[2]);
        g[31] = cp[10] * (g[27] + 2 * b01[3]);
        //g[32] = w[0];
        //g[33] = w[1];
        //g[34] = w[2];
        //g[35] = w[3];
        g[36] = cp[2] * g[32];
        g[37] = cp[5] * g[33];
        g[38] = cp[8] * g[34];
        g[39] = cp[11] * g[35];
        g[40] = cp[2] * g[36] + b01[0] * g[32];
        g[41] = cp[5] * g[37] + b01[1] * g[33];
        g[42] = cp[8] * g[38] + b01[2] * g[34];
        g[43] = cp[11] * g[39] + b01[3] * g[35];
        g[44] = cp[2] * g[40] + 2 * b01[0] * g[36];
        g[45] = cp[5] * g[41] + 2 * b01[1] * g[37];
        g[46] = cp[8] * g[42] + 2 * b01[2] * g[38];
        g[47] = cp[11] * g[43] + 2 * b01[3] * g[39];
}

static inline void _srg0_2d4d_0100(double *restrict g, struct _BC *bc, CINTEnvVars *envs)
{
        double *c0 = bc->c00;
        g[0] = 1;
        g[1] = 1;
        g[2] = c0[0];
        g[3] = c0[3];
        g[4] = 1;
        g[5] = 1;
        g[6] = c0[1];
        g[7] = c0[4];
        //g[8] = w[0];
        //g[9] = w[1];
        g[10] = c0[2] * g[8];
        g[11] = c0[5] * g[9];
}

static inline void _srg0_2d4d_0101(double *restrict g, struct _BC *bc, CINTEnvVars *envs)
{
        double *c0 = bc->c00;
        double *cp = bc->c0p;
        double *b00 = bc->b00;
        g[0] = 1;
        g[1] = 1;
        g[2] = 1;
        g[3] = 1;
        g[8] = c0[0];
        g[9] = c0[3];
        g[10] = c0[6];
        g[11] = c0[9];
        g[4] = cp[0];
        g[5] = cp[3];
        g[6] = cp[6];
        g[7] = cp[9];
        g[12] = cp[0] * c0[0] + b00[0];
        g[13] = cp[3] * c0[3] + b00[1];
        g[14] = cp[6] * c0[6] + b00[2];
        g[15] = cp[9] * c0[9] + b00[3];
        g[16] = 1;
        g[17] = 1;
        g[18] = 1;
        g[19] = 1;
        g[24] = c0[1];
        g[25] = c0[4];
        g[26] = c0[7];
        g[27] = c0[10];
        g[20] = cp[1];
        g[21] = cp[4];
        g[22] = cp[7];
        g[23] = cp[10];
        g[28] = cp[1] * c0[1] + b00[0];
        g[29] = cp[4] * c0[4] + b00[1];
        g[30] = cp[7] * c0[7] + b00[2];
        g[31] = cp[10] * c0[10] + b00[3];
        //g[32] = w[0];
        //g[33] = w[1];
        //g[34] = w[2];
        //g[35] = w[3];
        g[40] = c0[2] * g[32];
        g[41] = c0[5] * g[33];
        g[42] = c0[8] * g[34];
        g[43] = c0[11] * g[35];
        g[36] = cp[2] * g[32];
        g[37] = cp[5] * g[33];
        g[38] = cp[8] * g[34];
        g[39] = cp[11] * g[35];
        g[44] = cp[2] * g[40] + b00[0] * g[32];
        g[45] = cp[5] * g[41] + b00[1] * g[33];
        g[46] = cp[8] * g[42] + b00[2] * g[34];
        g[47] = cp[11] * g[43] + b00[3] * g[35];
}

static inline void _srg0_2d4d_0102(double *restrict g, struct _BC *bc, CINTEnvVars *envs)
{
        double *c0 = bc->c00;
        double *cp = bc->c0p;
        double *b00 = bc->b00;
        double *b01 = bc->b01;
        g[0] = 1;
        g[1] = 1;
        g[2] = 1;
        g[3] = 1;
        g[12] = c0[0];
        g[13] = c0[3];
        g[14] = c0[6];
        g[15] = c0[9];
        g[4] = cp[0];
        g[5] = cp[3];
        g[6] = cp[6];
        g[7] = cp[9];
        g[8] = cp[0] * cp[0] + b01[0];
        g[9] = cp[3] * cp[3] + b01[1];
        g[10] = cp[6] * cp[6] + b01[2];
        g[11] = cp[9] * cp[9] + b01[3];
        g[16] = cp[0] * c0[0] + b00[0];
        g[17] = cp[3] * c0[3] + b00[1];
        g[18] = cp[6] * c0[6] + b00[2];
        g[19] = cp[9] * c0[9] + b00[3];
        g[20] = cp[0] * (g[16] + b00[0]) + b01[0] * c0[0];
        g[21] = cp[3] * (g[17] + b00[1]) + b01[1] * c0[3];
        g[22] = cp[6] * (g[18] + b00[2]) + b01[2] * c0[6];
        g[23] = cp[9] * (g[19] + b00[3]) + b01[3] * c0[9];
        g[24] = 1;
        g[25] = 1;
        g[26] = 1;
        g[27] = 1;
        g[36] = c0[1];
        g[37] = c0[4];
        g[38] = c0[7];
        g[39] = c0[10];
        g[28] = cp[1];
        g[29] = cp[4];
        g[30] = cp[7];
        g[31] = cp[10];
        g[32] = cp[1] * cp[1] + b01[0];
        g[33] = cp[4] * cp[4] + b01[1];
        g[34] = cp[7] * cp[7] + b01[2];
        g[35] = cp[10] * cp[10] + b01[3];
        g[40] = cp[1] * c0[1] + b00[0];
        g[41] = cp[4] * c0[4] + b00[1];
        g[42] = cp[7] * c0[7] + b00[2];
        g[43] = cp[10] * c0[10] + b00[3];
        g[44] = cp[1] * (g[40] + b00[0]) + b01[0] * c0[1];
        g[45] = cp[4] * (g[41] + b00[1]) + b01[1] * c0[4];
        g[46] = cp[7] * (g[42] + b00[2]) + b01[2] * c0[7];
        g[47] = cp[10] * (g[43] + b00[3]) + b01[3] * c0[10];
        //g[48] = w[0];
        //g[49] = w[1];
        //g[50] = w[2];
        //g[51] = w[3];
        g[60] = c0[2] * g[48];
        g[61] = c0[5] * g[49];
        g[62] = c0[8] * g[50];
        g[63] = c0[11] * g[51];
        g[52] = cp[2] * g[48];
        g[53] = cp[5] * g[49];
        g[54] = cp[8] * g[50];
        g[55] = cp[11] * g[51];
        g[56] = cp[2] * g[52] + b01[0] * g[48];
        g[57] = cp[5] * g[53] + b01[1] * g[49];
        g[58] = cp[8] * g[54] + b01[2] * g[50];
        g[59] = cp[11] * g[55] + b01[3] * g[51];
        g[64] = cp[2] * g[60] + b00[0] * g[48];
        g[65] = cp[5] * g[61] + b00[1] * g[49];
        g[66] = cp[8] * g[62] + b00[2] * g[50];
        g[67] = cp[11] * g[63] + b00[3] * g[51];
        g[68] = cp[2] * g[64] + b01[0] * g[60] + b00[0] * g[52];
        g[69] = cp[5] * g[65] + b01[1] * g[61] + b00[1] * g[53];
        g[70] = cp[8] * g[66] + b01[2] * g[62] + b00[2] * g[54];
        g[71] = cp[11] * g[67] + b01[3] * g[63] + b00[3] * g[55];
}

static inline void _srg0_2d4d_0110(double *restrict g, struct _BC *bc, CINTEnvVars *envs)
{
        double *c0 = bc->c00;
        double *cp = bc->c0p;
        double *b00 = bc->b00;
        g[0] = 1;
        g[1] = 1;
        g[2] = 1;
        g[3] = 1;
        g[8] = c0[0];
        g[9] = c0[3];
        g[10] = c0[6];
        g[11] = c0[9];
        g[4] = cp[0];
        g[5] = cp[3];
        g[6] = cp[6];
        g[7] = cp[9];
        g[12] = cp[0] * c0[0] + b00[0];
        g[13] = cp[3] * c0[3] + b00[1];
        g[14] = cp[6] * c0[6] + b00[2];
        g[15] = cp[9] * c0[9] + b00[3];
        g[16] = 1;
        g[17] = 1;
        g[18] = 1;
        g[19] = 1;
        g[24] = c0[1];
        g[25] = c0[4];
        g[26] = c0[7];
        g[27] = c0[10];
        g[20] = cp[1];
        g[21] = cp[4];
        g[22] = cp[7];
        g[23] = cp[10];
        g[28] = cp[1] * c0[1] + b00[0];
        g[29] = cp[4] * c0[4] + b00[1];
        g[30] = cp[7] * c0[7] + b00[2];
        g[31] = cp[10] * c0[10] + b00[3];
        //g[32] = w[0];
        //g[33] = w[1];
        //g[34] = w[2];
        //g[35] = w[3];
        g[40] = c0[2] * g[32];
        g[41] = c0[5] * g[33];
        g[42] = c0[8] * g[34];
        g[43] = c0[11] * g[35];
        g[36] = cp[2] * g[32];
        g[37] = cp[5] * g[33];
        g[38] = cp[8] * g[34];
        g[39] = cp[11] * g[35];
        g[44] = cp[2] * g[40] + b00[0] * g[32];
        g[45] = cp[5] * g[41] + b00[1] * g[33];
        g[46] = cp[8] * g[42] + b00[2] * g[34];
        g[47] = cp[11] * g[43] + b00[3] * g[35];
}

static inline void _srg0_2d4d_0111(double *restrict g, struct _BC *bc, CINTEnvVars *envs)
{
        double *c0 = bc->c00;
        double *cp = bc->c0p;
        double *b00 = bc->b00;
        double *b01 = bc->b01;
        double xkxl = envs->rkrl[0];
        double ykyl = envs->rkrl[1];
        double zkzl = envs->rkrl[2];
        g[0] = 1;
        g[1] = 1;
        g[2] = 1;
        g[3] = 1;
        g[24] = c0[0];
        g[25] = c0[3];
        g[26] = c0[6];
        g[27] = c0[9];
        g[8] = cp[0];
        g[9] = cp[3];
        g[10] = cp[6];
        g[11] = cp[9];
        g[32] = cp[0] * c0[0] + b00[0];
        g[33] = cp[3] * c0[3] + b00[1];
        g[34] = cp[6] * c0[6] + b00[2];
        g[35] = cp[9] * c0[9] + b00[3];
        g[12] = cp[0] * (xkxl + cp[0]) + b01[0];
        g[13] = cp[3] * (xkxl + cp[3]) + b01[1];
        g[14] = cp[6] * (xkxl + cp[6]) + b01[2];
        g[15] = cp[9] * (xkxl + cp[9]) + b01[3];
        g[36] = g[32] * (xkxl + cp[0]) + cp[0] * b00[0] + b01[0] * c0[0];
        g[37] = g[33] * (xkxl + cp[3]) + cp[3] * b00[1] + b01[1] * c0[3];
        g[38] = g[34] * (xkxl + cp[6]) + cp[6] * b00[2] + b01[2] * c0[6];
        g[39] = g[35] * (xkxl + cp[9]) + cp[9] * b00[3] + b01[3] * c0[9];
        g[4] = xkxl + cp[0];
        g[5] = xkxl + cp[3];
        g[6] = xkxl + cp[6];
        g[7] = xkxl + cp[9];
        g[28] = c0[0] * (xkxl + cp[0]) + b00[0];
        g[29] = c0[3] * (xkxl + cp[3]) + b00[1];
        g[30] = c0[6] * (xkxl + cp[6]) + b00[2];
        g[31] = c0[9] * (xkxl + cp[9]) + b00[3];
        g[48] = 1;
        g[49] = 1;
        g[50] = 1;
        g[51] = 1;
        g[72] = c0[1];
        g[73] = c0[4];
        g[74] = c0[7];
        g[75] = c0[10];
        g[56] = cp[1];
        g[57] = cp[4];
        g[58] = cp[7];
        g[59] = cp[10];
        g[80] = cp[1] * c0[1] + b00[0];
        g[81] = cp[4] * c0[4] + b00[1];
        g[82] = cp[7] * c0[7] + b00[2];
        g[83] = cp[10] * c0[10] + b00[3];
        g[60] = cp[1] * (ykyl + cp[1]) + b01[0];
        g[61] = cp[4] * (ykyl + cp[4]) + b01[1];
        g[62] = cp[7] * (ykyl + cp[7]) + b01[2];
        g[63] = cp[10] * (ykyl + cp[10]) + b01[3];
        g[84] = g[80] * (ykyl + cp[1]) + cp[1] * b00[0] + b01[0] * c0[1];
        g[85] = g[81] * (ykyl + cp[4]) + cp[4] * b00[1] + b01[1] * c0[4];
        g[86] = g[82] * (ykyl + cp[7]) + cp[7] * b00[2] + b01[2] * c0[7];
        g[87] = g[83] * (ykyl + cp[10]) + cp[10] * b00[3] + b01[3] * c0[10];
        g[52] = ykyl + cp[1];
        g[53] = ykyl + cp[4];
        g[54] = ykyl + cp[7];
        g[55] = ykyl + cp[10];
        g[76] = c0[1] * (ykyl + cp[1]) + b00[0];
        g[77] = c0[4] * (ykyl + cp[4]) + b00[1];
        g[78] = c0[7] * (ykyl + cp[7]) + b00[2];
        g[79] = c0[10] * (ykyl + cp[10]) + b00[3];
        //g[96] = w[0];
        //g[97] = w[1];
        //g[98] = w[2];
        //g[99] = w[3];
        g[120] = c0[2] * g[96];
        g[121] = c0[5] * g[97];
        g[122] = c0[8] * g[98];
        g[123] = c0[11] * g[99];
        g[104] = cp[2] * g[96];
        g[105] = cp[5] * g[97];
        g[106] = cp[8] * g[98];
        g[107] = cp[11] * g[99];
        g[128] = cp[2] * g[120] + b00[0] * g[96];
        g[129] = cp[5] * g[121] + b00[1] * g[97];
        g[130] = cp[8] * g[122] + b00[2] * g[98];
        g[131] = cp[11] * g[123] + b00[3] * g[99];
        g[108] = g[104] * (zkzl + cp[2]) + b01[0] * g[96];
        g[109] = g[105] * (zkzl + cp[5]) + b01[1] * g[97];
        g[110] = g[106] * (zkzl + cp[8]) + b01[2] * g[98];
        g[111] = g[107] * (zkzl + cp[11]) + b01[3] * g[99];
        g[132] = g[128] * (zkzl + cp[2]) + b01[0] * g[120] + b00[0] * g[104];
        g[133] = g[129] * (zkzl + cp[5]) + b01[1] * g[121] + b00[1] * g[105];
        g[134] = g[130] * (zkzl + cp[8]) + b01[2] * g[122] + b00[2] * g[106];
        g[135] = g[131] * (zkzl + cp[11]) + b01[3] * g[123] + b00[3] * g[107];
        g[100] = g[96] * (zkzl + cp[2]);
        g[101] = g[97] * (zkzl + cp[5]);
        g[102] = g[98] * (zkzl + cp[8]);
        g[103] = g[99] * (zkzl + cp[11]);
        g[124] = g[120] * (zkzl + cp[2]) + b00[0] * g[96];
        g[125] = g[121] * (zkzl + cp[5]) + b00[1] * g[97];
        g[126] = g[122] * (zkzl + cp[8]) + b00[2] * g[98];
        g[127] = g[123] * (zkzl + cp[11]) + b00[3] * g[99];
}

static inline void _srg0_2d4d_0120(double *restrict g, struct _BC *bc, CINTEnvVars *envs)
{
        double *c0 = bc->c00;
        double *cp = bc->c0p;
        double *b00 = bc->b00;
        double *b01 = bc->b01;
        g[0] = 1;
        g[1] = 1;
        g[2] = 1;
        g[3] = 1;
        g[12] = c0[0];
        g[13] = c0[3];
        g[14] = c0[6];
        g[15] = c0[9];
        g[4] = cp[0];
        g[5] = cp[3];
        g[6] = cp[6];
        g[7] = cp[9];
        g[8] = cp[0] * cp[0] + b01[0];
        g[9] = cp[3] * cp[3] + b01[1];
        g[10] = cp[6] * cp[6] + b01[2];
        g[11] = cp[9] * cp[9] + b01[3];
        g[16] = cp[0] * c0[0] + b00[0];
        g[17] = cp[3] * c0[3] + b00[1];
        g[18] = cp[6] * c0[6] + b00[2];
        g[19] = cp[9] * c0[9] + b00[3];
        g[20] = cp[0] * (g[16] + b00[0]) + b01[0] * c0[0];
        g[21] = cp[3] * (g[17] + b00[1]) + b01[1] * c0[3];
        g[22] = cp[6] * (g[18] + b00[2]) + b01[2] * c0[6];
        g[23] = cp[9] * (g[19] + b00[3]) + b01[3] * c0[9];
        g[24] = 1;
        g[25] = 1;
        g[26] = 1;
        g[27] = 1;
        g[36] = c0[1];
        g[37] = c0[4];
        g[38] = c0[7];
        g[39] = c0[10];
        g[28] = cp[1];
        g[29] = cp[4];
        g[30] = cp[7];
        g[31] = cp[10];
        g[32] = cp[1] * cp[1] + b01[0];
        g[33] = cp[4] * cp[4] + b01[1];
        g[34] = cp[7] * cp[7] + b01[2];
        g[35] = cp[10] * cp[10] + b01[3];
        g[40] = cp[1] * c0[1] + b00[0];
        g[41] = cp[4] * c0[4] + b00[1];
        g[42] = cp[7] * c0[7] + b00[2];
        g[43] = cp[10] * c0[10] + b00[3];
        g[44] = cp[1] * (g[40] + b00[0]) + b01[0] * c0[1];
        g[45] = cp[4] * (g[41] + b00[1]) + b01[1] * c0[4];
        g[46] = cp[7] * (g[42] + b00[2]) + b01[2] * c0[7];
        g[47] = cp[10] * (g[43] + b00[3]) + b01[3] * c0[10];
        //g[48] = w[0];
        //g[49] = w[1];
        //g[50] = w[2];
        //g[51] = w[3];
        g[60] = c0[2] * g[48];
        g[61] = c0[5] * g[49];
        g[62] = c0[8] * g[50];
        g[63] = c0[11] * g[51];
        g[52] = cp[2] * g[48];
        g[53] = cp[5] * g[49];
        g[54] = cp[8] * g[50];
        g[55] = cp[11] * g[51];
        g[56] = cp[2] * g[52] + b01[0] * g[48];
        g[57] = cp[5] * g[53] + b01[1] * g[49];
        g[58] = cp[8] * g[54] + b01[2] * g[50];
        g[59] = cp[11] * g[55] + b01[3] * g[51];
        g[64] = cp[2] * g[60] + b00[0] * g[48];
        g[65] = cp[5] * g[61] + b00[1] * g[49];
        g[66] = cp[8] * g[62] + b00[2] * g[50];
        g[67] = cp[11] * g[63] + b00[3] * g[51];
        g[68] = cp[2] * g[64] + b01[0] * g[60] + b00[0] * g[52];
        g[69] = cp[5] * g[65] + b01[1] * g[61] + b00[1] * g[53];
        g[70] = cp[8] * g[66] + b01[2] * g[62] + b00[2] * g[54];
        g[71] = cp[11] * g[67] + b01[3] * g[63] + b00[3] * g[55];
}

static inline void _srg0_2d4d_0200(double *restrict g, struct _BC *bc, CINTEnvVars *envs)
{
        double *c0 = bc->c00;
        double *b10 = bc->b10;
        g[0] = 1;
        g[1] = 1;
        g[2] = 1;
        g[3] = 1;
        g[4] = c0[0];
        g[5] = c0[3];
        g[6] = c0[6];
        g[7] = c0[9];
        g[8] = c0[0] * c0[0] + b10[0];
        g[9] = c0[3] * c0[3] + b10[1];
        g[10] = c0[6] * c0[6] + b10[2];
        g[11] = c0[9] * c0[9] + b10[3];
        g[12] = 1;
        g[13] = 1;
        g[14] = 1;
        g[15] = 1;
        g[16] = c0[1];
        g[17] = c0[4];
        g[18] = c0[7];
        g[19] = c0[10];
        g[20] = c0[1] * c0[1] + b10[0];
        g[21] = c0[4] * c0[4] + b10[1];
        g[22] = c0[7] * c0[7] + b10[2];
        g[23] = c0[10] * c0[10] + b10[3];
        //g[24] = w[0];
        //g[25] = w[1];
        //g[26] = w[2];
        //g[27] = w[3];
        g[28] = c0[2] * g[24];
        g[29] = c0[5] * g[25];
        g[30] = c0[8] * g[26];
        g[31] = c0[11] * g[27];
        g[32] = c0[2] * g[28] + b10[0] * g[24];
        g[33] = c0[5] * g[29] + b10[1] * g[25];
        g[34] = c0[8] * g[30] + b10[2] * g[26];
        g[35] = c0[11] * g[31] + b10[3] * g[27];
}

static inline void _srg0_2d4d_0201(double *restrict g, struct _BC *bc, CINTEnvVars *envs)
{
        double *c0 = bc->c00;
        double *cp = bc->c0p;
        double *b00 = bc->b00;
        double *b10 = bc->b10;
        g[0] = 1;
        g[1] = 1;
        g[2] = 1;
        g[3] = 1;
        g[8] = c0[0];
        g[9] = c0[3];
        g[10] = c0[6];
        g[11] = c0[9];
        g[16] = c0[0] * c0[0] + b10[0];
        g[17] = c0[3] * c0[3] + b10[1];
        g[18] = c0[6] * c0[6] + b10[2];
        g[19] = c0[9] * c0[9] + b10[3];
        g[4] = cp[0];
        g[5] = cp[3];
        g[6] = cp[6];
        g[7] = cp[9];
        g[12] = cp[0] * c0[0] + b00[0];
        g[13] = cp[3] * c0[3] + b00[1];
        g[14] = cp[6] * c0[6] + b00[2];
        g[15] = cp[9] * c0[9] + b00[3];
        g[20] = c0[0] * (g[12] + b00[0]) + b10[0] * cp[0];
        g[21] = c0[3] * (g[13] + b00[1]) + b10[1] * cp[3];
        g[22] = c0[6] * (g[14] + b00[2]) + b10[2] * cp[6];
        g[23] = c0[9] * (g[15] + b00[3]) + b10[3] * cp[9];
        g[24] = 1;
        g[25] = 1;
        g[26] = 1;
        g[27] = 1;
        g[32] = c0[1];
        g[33] = c0[4];
        g[34] = c0[7];
        g[35] = c0[10];
        g[40] = c0[1] * c0[1] + b10[0];
        g[41] = c0[4] * c0[4] + b10[1];
        g[42] = c0[7] * c0[7] + b10[2];
        g[43] = c0[10] * c0[10] + b10[3];
        g[28] = cp[1];
        g[29] = cp[4];
        g[30] = cp[7];
        g[31] = cp[10];
        g[36] = cp[1] * c0[1] + b00[0];
        g[37] = cp[4] * c0[4] + b00[1];
        g[38] = cp[7] * c0[7] + b00[2];
        g[39] = cp[10] * c0[10] + b00[3];
        g[44] = c0[1] * (g[36] + b00[0]) + b10[0] * cp[1];
        g[45] = c0[4] * (g[37] + b00[1]) + b10[1] * cp[4];
        g[46] = c0[7] * (g[38] + b00[2]) + b10[2] * cp[7];
        g[47] = c0[10] * (g[39] + b00[3]) + b10[3] * cp[10];
        //g[48] = w[0];
        //g[49] = w[1];
        //g[50] = w[2];
        //g[51] = w[3];
        g[56] = c0[2] * g[48];
        g[57] = c0[5] * g[49];
        g[58] = c0[8] * g[50];
        g[59] = c0[11] * g[51];
        g[64] = c0[2] * g[56] + b10[0] * g[48];
        g[65] = c0[5] * g[57] + b10[1] * g[49];
        g[66] = c0[8] * g[58] + b10[2] * g[50];
        g[67] = c0[11] * g[59] + b10[3] * g[51];
        g[52] = cp[2] * g[48];
        g[53] = cp[5] * g[49];
        g[54] = cp[8] * g[50];
        g[55] = cp[11] * g[51];
        g[60] = cp[2] * g[56] + b00[0] * g[48];
        g[61] = cp[5] * g[57] + b00[1] * g[49];
        g[62] = cp[8] * g[58] + b00[2] * g[50];
        g[63] = cp[11] * g[59] + b00[3] * g[51];
        g[68] = c0[2] * g[60] + b10[0] * g[52] + b00[0] * g[56];
        g[69] = c0[5] * g[61] + b10[1] * g[53] + b00[1] * g[57];
        g[70] = c0[8] * g[62] + b10[2] * g[54] + b00[2] * g[58];
        g[71] = c0[11] * g[63] + b10[3] * g[55] + b00[3] * g[59];
}

static inline void _srg0_2d4d_0210(double *restrict g, struct _BC *bc, CINTEnvVars *envs)
{
        double *c0 = bc->c00;
        double *cp = bc->c0p;
        double *b00 = bc->b00;
        double *b10 = bc->b10;
        g[0] = 1;
        g[1] = 1;
        g[2] = 1;
        g[3] = 1;
        g[8] = c0[0];
        g[9] = c0[3];
        g[10] = c0[6];
        g[11] = c0[9];
        g[16] = c0[0] * c0[0] + b10[0];
        g[17] = c0[3] * c0[3] + b10[1];
        g[18] = c0[6] * c0[6] + b10[2];
        g[19] = c0[9] * c0[9] + b10[3];
        g[4] = cp[0];
        g[5] = cp[3];
        g[6] = cp[6];
        g[7] = cp[9];
        g[12] = cp[0] * c0[0] + b00[0];
        g[13] = cp[3] * c0[3] + b00[1];
        g[14] = cp[6] * c0[6] + b00[2];
        g[15] = cp[9] * c0[9] + b00[3];
        g[20] = c0[0] * (g[12] + b00[0]) + b10[0] * cp[0];
        g[21] = c0[3] * (g[13] + b00[1]) + b10[1] * cp[3];
        g[22] = c0[6] * (g[14] + b00[2]) + b10[2] * cp[6];
        g[23] = c0[9] * (g[15] + b00[3]) + b10[3] * cp[9];
        g[24] = 1;
        g[25] = 1;
        g[26] = 1;
        g[27] = 1;
        g[32] = c0[1];
        g[33] = c0[4];
        g[34] = c0[7];
        g[35] = c0[10];
        g[40] = c0[1] * c0[1] + b10[0];
        g[41] = c0[4] * c0[4] + b10[1];
        g[42] = c0[7] * c0[7] + b10[2];
        g[43] = c0[10] * c0[10] + b10[3];
        g[28] = cp[1];
        g[29] = cp[4];
        g[30] = cp[7];
        g[31] = cp[10];
        g[36] = cp[1] * c0[1] + b00[0];
        g[37] = cp[4] * c0[4] + b00[1];
        g[38] = cp[7] * c0[7] + b00[2];
        g[39] = cp[10] * c0[10] + b00[3];
        g[44] = c0[1] * (g[36] + b00[0]) + b10[0] * cp[1];
        g[45] = c0[4] * (g[37] + b00[1]) + b10[1] * cp[4];
        g[46] = c0[7] * (g[38] + b00[2]) + b10[2] * cp[7];
        g[47] = c0[10] * (g[39] + b00[3]) + b10[3] * cp[10];
        //g[48] = w[0];
        //g[49] = w[1];
        //g[50] = w[2];
        //g[51] = w[3];
        g[56] = c0[2] * g[48];
        g[57] = c0[5] * g[49];
        g[58] = c0[8] * g[50];
        g[59] = c0[11] * g[51];
        g[64] = c0[2] * g[56] + b10[0] * g[48];
        g[65] = c0[5] * g[57] + b10[1] * g[49];
        g[66] = c0[8] * g[58] + b10[2] * g[50];
        g[67] = c0[11] * g[59] + b10[3] * g[51];
        g[52] = cp[2] * g[48];
        g[53] = cp[5] * g[49];
        g[54] = cp[8] * g[50];
        g[55] = cp[11] * g[51];
        g[60] = cp[2] * g[56] + b00[0] * g[48];
        g[61] = cp[5] * g[57] + b00[1] * g[49];
        g[62] = cp[8] * g[58] + b00[2] * g[50];
        g[63] = cp[11] * g[59] + b00[3] * g[51];
        g[68] = c0[2] * g[60] + b10[0] * g[52] + b00[0] * g[56];
        g[69] = c0[5] * g[61] + b10[1] * g[53] + b00[1] * g[57];
        g[70] = c0[8] * g[62] + b10[2] * g[54] + b00[2] * g[58];
        g[71] = c0[11] * g[63] + b10[3] * g[55] + b00[3] * g[59];
}

static inline void _srg0_2d4d_0300(double *restrict g, struct _BC *bc, CINTEnvVars *envs)
{
        double *c0 = bc->c00;
        double *b10 = bc->b10;
        g[0] = 1;
        g[1] = 1;
        g[2] = 1;
        g[3] = 1;
        g[4] = c0[0];
        g[5] = c0[3];
        g[6] = c0[6];
        g[7] = c0[9];
        g[8] = c0[0] * c0[0] + b10[0];
        g[9] = c0[3] * c0[3] + b10[1];
        g[10] = c0[6] * c0[6] + b10[2];
        g[11] = c0[9] * c0[9] + b10[3];
        g[12] = c0[0] * (g[8] + 2 * b10[0]);
        g[13] = c0[3] * (g[9] + 2 * b10[1]);
        g[14] = c0[6] * (g[10] + 2 * b10[2]);
        g[15] = c0[9] * (g[11] + 2 * b10[3]);
        g[16] = 1;
        g[17] = 1;
        g[18] = 1;
        g[19] = 1;
        g[20] = c0[1];
        g[21] = c0[4];
        g[22] = c0[7];
        g[23] = c0[10];
        g[24] = c0[1] * c0[1] + b10[0];
        g[25] = c0[4] * c0[4] + b10[1];
        g[26] = c0[7] * c0[7] + b10[2];
        g[27] = c0[10] * c0[10] + b10[3];
        g[28] = c0[1] * (g[24] + 2 * b10[0]);
        g[29] = c0[4] * (g[25] + 2 * b10[1]);
        g[30] = c0[7] * (g[26] + 2 * b10[2]);
        g[31] = c0[10] * (g[27] + 2 * b10[3]);
        //g[32] = w[0];
        //g[33] = w[1];
        //g[34] = w[2];
        //g[35] = w[3];
        g[36] = c0[2] * g[32];
        g[37] = c0[5] * g[33];
        g[38] = c0[8] * g[34];
        g[39] = c0[11] * g[35];
        g[40] = c0[2] * g[36] + b10[0] * g[32];
        g[41] = c0[5] * g[37] + b10[1] * g[33];
        g[42] = c0[8] * g[38] + b10[2] * g[34];
        g[43] = c0[11] * g[39] + b10[3] * g[35];
        g[44] = c0[2] * g[40] + 2 * b10[0] * g[36];
        g[45] = c0[5] * g[41] + 2 * b10[1] * g[37];
        g[46] = c0[8] * g[42] + 2 * b10[2] * g[38];
        g[47] = c0[11] * g[43] + 2 * b10[3] * g[39];
}

static inline void _srg0_2d4d_1000(double *restrict g, struct _BC *bc, CINTEnvVars *envs)
{
        double *c0 = bc->c00;
        g[0] = 1;
        g[1] = 1;
        g[2] = c0[0];
        g[3] = c0[3];
        g[4] = 1;
        g[5] = 1;
        g[6] = c0[1];
        g[7] = c0[4];
        //g[8] = w[0];
        //g[9] = w[1];
        g[10] = c0[2] * g[8];
        g[11] = c0[5] * g[9];
}

static inline void _srg0_2d4d_1001(double *restrict g, struct _BC *bc, CINTEnvVars *envs)
{
        double *c0 = bc->c00;
        double *cp = bc->c0p;
        double *b00 = bc->b00;
        g[0] = 1;
        g[1] = 1;
        g[2] = 1;
        g[3] = 1;
        g[4] = c0[0];
        g[5] = c0[3];
        g[6] = c0[6];
        g[7] = c0[9];
        g[8] = cp[0];
        g[9] = cp[3];
        g[10] = cp[6];
        g[11] = cp[9];
        g[12] = cp[0] * c0[0] + b00[0];
        g[13] = cp[3] * c0[3] + b00[1];
        g[14] = cp[6] * c0[6] + b00[2];
        g[15] = cp[9] * c0[9] + b00[3];
        g[16] = 1;
        g[17] = 1;
        g[18] = 1;
        g[19] = 1;
        g[20] = c0[1];
        g[21] = c0[4];
        g[22] = c0[7];
        g[23] = c0[10];
        g[24] = cp[1];
        g[25] = cp[4];
        g[26] = cp[7];
        g[27] = cp[10];
        g[28] = cp[1] * c0[1] + b00[0];
        g[29] = cp[4] * c0[4] + b00[1];
        g[30] = cp[7] * c0[7] + b00[2];
        g[31] = cp[10] * c0[10] + b00[3];
        //g[32] = w[0];
        //g[33] = w[1];
        //g[34] = w[2];
        //g[35] = w[3];
        g[36] = c0[2] * g[32];
        g[37] = c0[5] * g[33];
        g[38] = c0[8] * g[34];
        g[39] = c0[11] * g[35];
        g[40] = cp[2] * g[32];
        g[41] = cp[5] * g[33];
        g[42] = cp[8] * g[34];
        g[43] = cp[11] * g[35];
        g[44] = cp[2] * g[36] + b00[0] * g[32];
        g[45] = cp[5] * g[37] + b00[1] * g[33];
        g[46] = cp[8] * g[38] + b00[2] * g[34];
        g[47] = cp[11] * g[39] + b00[3] * g[35];
}

static inline void _srg0_2d4d_1002(double *restrict g, struct _BC *bc, CINTEnvVars *envs)
{
        double *c0 = bc->c00;
        double *cp = bc->c0p;
        double *b00 = bc->b00;
        double *b01 = bc->b01;
        g[0] = 1;
        g[1] = 1;
        g[2] = 1;
        g[3] = 1;
        g[4] = c0[0];
        g[5] = c0[3];
        g[6] = c0[6];
        g[7] = c0[9];
        g[8] = cp[0];
        g[9] = cp[3];
        g[10] = cp[6];
        g[11] = cp[9];
        g[16] = cp[0] * cp[0] + b01[0];
        g[17] = cp[3] * cp[3] + b01[1];
        g[18] = cp[6] * cp[6] + b01[2];
        g[19] = cp[9] * cp[9] + b01[3];
        g[12] = cp[0] * c0[0] + b00[0];
        g[13] = cp[3] * c0[3] + b00[1];
        g[14] = cp[6] * c0[6] + b00[2];
        g[15] = cp[9] * c0[9] + b00[3];
        g[20] = cp[0] * (g[12] + b00[0]) + b01[0] * c0[0];
        g[21] = cp[3] * (g[13] + b00[1]) + b01[1] * c0[3];
        g[22] = cp[6] * (g[14] + b00[2]) + b01[2] * c0[6];
        g[23] = cp[9] * (g[15] + b00[3]) + b01[3] * c0[9];
        g[24] = 1;
        g[25] = 1;
        g[26] = 1;
        g[27] = 1;
        g[28] = c0[1];
        g[29] = c0[4];
        g[30] = c0[7];
        g[31] = c0[10];
        g[32] = cp[1];
        g[33] = cp[4];
        g[34] = cp[7];
        g[35] = cp[10];
        g[40] = cp[1] * cp[1] + b01[0];
        g[41] = cp[4] * cp[4] + b01[1];
        g[42] = cp[7] * cp[7] + b01[2];
        g[43] = cp[10] * cp[10] + b01[3];
        g[36] = cp[1] * c0[1] + b00[0];
        g[37] = cp[4] * c0[4] + b00[1];
        g[38] = cp[7] * c0[7] + b00[2];
        g[39] = cp[10] * c0[10] + b00[3];
        g[44] = cp[1] * (g[36] + b00[0]) + b01[0] * c0[1];
        g[45] = cp[4] * (g[37] + b00[1]) + b01[1] * c0[4];
        g[46] = cp[7] * (g[38] + b00[2]) + b01[2] * c0[7];
        g[47] = cp[10] * (g[39] + b00[3]) + b01[3] * c0[10];
        //g[48] = w[0];
        //g[49] = w[1];
        //g[50] = w[2];
        //g[51] = w[3];
        g[52] = c0[2] * g[48];
        g[53] = c0[5] * g[49];
        g[54] = c0[8] * g[50];
        g[55] = c0[11] * g[51];
        g[56] = cp[2] * g[48];
        g[57] = cp[5] * g[49];
        g[58] = cp[8] * g[50];
        g[59] = cp[11] * g[51];
        g[64] = cp[2] * g[56] + b01[0] * g[48];
        g[65] = cp[5] * g[57] + b01[1] * g[49];
        g[66] = cp[8] * g[58] + b01[2] * g[50];
        g[67] = cp[11] * g[59] + b01[3] * g[51];
        g[60] = cp[2] * g[52] + b00[0] * g[48];
        g[61] = cp[5] * g[53] + b00[1] * g[49];
        g[62] = cp[8] * g[54] + b00[2] * g[50];
        g[63] = cp[11] * g[55] + b00[3] * g[51];
        g[68] = cp[2] * g[60] + b01[0] * g[52] + b00[0] * g[56];
        g[69] = cp[5] * g[61] + b01[1] * g[53] + b00[1] * g[57];
        g[70] = cp[8] * g[62] + b01[2] * g[54] + b00[2] * g[58];
        g[71] = cp[11] * g[63] + b01[3] * g[55] + b00[3] * g[59];
}

static inline void _srg0_2d4d_1010(double *restrict g, struct _BC *bc, CINTEnvVars *envs)
{
        double *c0 = bc->c00;
        double *cp = bc->c0p;
        double *b00 = bc->b00;
        g[0] = 1;
        g[1] = 1;
        g[2] = 1;
        g[3] = 1;
        g[4] = c0[0];
        g[5] = c0[3];
        g[6] = c0[6];
        g[7] = c0[9];
        g[8] = cp[0];
        g[9] = cp[3];
        g[10] = cp[6];
        g[11] = cp[9];
        g[12] = cp[0] * c0[0] + b00[0];
        g[13] = cp[3] * c0[3] + b00[1];
        g[14] = cp[6] * c0[6] + b00[2];
        g[15] = cp[9] * c0[9] + b00[3];
        g[16] = 1;
        g[17] = 1;
        g[18] = 1;
        g[19] = 1;
        g[20] = c0[1];
        g[21] = c0[4];
        g[22] = c0[7];
        g[23] = c0[10];
        g[24] = cp[1];
        g[25] = cp[4];
        g[26] = cp[7];
        g[27] = cp[10];
        g[28] = cp[1] * c0[1] + b00[0];
        g[29] = cp[4] * c0[4] + b00[1];
        g[30] = cp[7] * c0[7] + b00[2];
        g[31] = cp[10] * c0[10] + b00[3];
        //g[32] = w[0];
        //g[33] = w[1];
        //g[34] = w[2];
        //g[35] = w[3];
        g[36] = c0[2] * g[32];
        g[37] = c0[5] * g[33];
        g[38] = c0[8] * g[34];
        g[39] = c0[11] * g[35];
        g[40] = cp[2] * g[32];
        g[41] = cp[5] * g[33];
        g[42] = cp[8] * g[34];
        g[43] = cp[11] * g[35];
        g[44] = cp[2] * g[36] + b00[0] * g[32];
        g[45] = cp[5] * g[37] + b00[1] * g[33];
        g[46] = cp[8] * g[38] + b00[2] * g[34];
        g[47] = cp[11] * g[39] + b00[3] * g[35];
}

static inline void _srg0_2d4d_1011(double *restrict g, struct _BC *bc, CINTEnvVars *envs)
{
        double *c0 = bc->c00;
        double *cp = bc->c0p;
        double *b00 = bc->b00;
        double *b01 = bc->b01;
        double xkxl = envs->rkrl[0];
        double ykyl = envs->rkrl[1];
        double zkzl = envs->rkrl[2];
        g[0] = 1;
        g[1] = 1;
        g[2] = 1;
        g[3] = 1;
        g[4] = c0[0];
        g[5] = c0[3];
        g[6] = c0[6];
        g[7] = c0[9];
        g[16] = cp[0];
        g[17] = cp[3];
        g[18] = cp[6];
        g[19] = cp[9];
        g[20] = cp[0] * c0[0] + b00[0];
        g[21] = cp[3] * c0[3] + b00[1];
        g[22] = cp[6] * c0[6] + b00[2];
        g[23] = cp[9] * c0[9] + b00[3];
        g[24] = cp[0] * (xkxl + cp[0]) + b01[0];
        g[25] = cp[3] * (xkxl + cp[3]) + b01[1];
        g[26] = cp[6] * (xkxl + cp[6]) + b01[2];
        g[27] = cp[9] * (xkxl + cp[9]) + b01[3];
        g[28] = g[20] * (xkxl + cp[0]) + cp[0] * b00[0] + b01[0] * c0[0];
        g[29] = g[21] * (xkxl + cp[3]) + cp[3] * b00[1] + b01[1] * c0[3];
        g[30] = g[22] * (xkxl + cp[6]) + cp[6] * b00[2] + b01[2] * c0[6];
        g[31] = g[23] * (xkxl + cp[9]) + cp[9] * b00[3] + b01[3] * c0[9];
        g[8] = xkxl + cp[0];
        g[9] = xkxl + cp[3];
        g[10] = xkxl + cp[6];
        g[11] = xkxl + cp[9];
        g[12] = c0[0] * (xkxl + cp[0]) + b00[0];
        g[13] = c0[3] * (xkxl + cp[3]) + b00[1];
        g[14] = c0[6] * (xkxl + cp[6]) + b00[2];
        g[15] = c0[9] * (xkxl + cp[9]) + b00[3];
        g[48] = 1;
        g[49] = 1;
        g[50] = 1;
        g[51] = 1;
        g[52] = c0[1];
        g[53] = c0[4];
        g[54] = c0[7];
        g[55] = c0[10];
        g[64] = cp[1];
        g[65] = cp[4];
        g[66] = cp[7];
        g[67] = cp[10];
        g[68] = cp[1] * c0[1] + b00[0];
        g[69] = cp[4] * c0[4] + b00[1];
        g[70] = cp[7] * c0[7] + b00[2];
        g[71] = cp[10] * c0[10] + b00[3];
        g[72] = cp[1] * (ykyl + cp[1]) + b01[0];
        g[73] = cp[4] * (ykyl + cp[4]) + b01[1];
        g[74] = cp[7] * (ykyl + cp[7]) + b01[2];
        g[75] = cp[10] * (ykyl + cp[10]) + b01[3];
        g[76] = g[68] * (ykyl + cp[1]) + cp[1] * b00[0] + b01[0] * c0[1];
        g[77] = g[69] * (ykyl + cp[4]) + cp[4] * b00[1] + b01[1] * c0[4];
        g[78] = g[70] * (ykyl + cp[7]) + cp[7] * b00[2] + b01[2] * c0[7];
        g[79] = g[71] * (ykyl + cp[10]) + cp[10] * b00[3] + b01[3] * c0[10];
        g[56] = ykyl + cp[1];
        g[57] = ykyl + cp[4];
        g[58] = ykyl + cp[7];
        g[59] = ykyl + cp[10];
        g[60] = c0[1] * (ykyl + cp[1]) + b00[0];
        g[61] = c0[4] * (ykyl + cp[4]) + b00[1];
        g[62] = c0[7] * (ykyl + cp[7]) + b00[2];
        g[63] = c0[10] * (ykyl + cp[10]) + b00[3];
        //g[96] = w[0];
        //g[97] = w[1];
        //g[98] = w[2];
        //g[99] = w[3];
        g[100] = c0[2] * g[96];
        g[101] = c0[5] * g[97];
        g[102] = c0[8] * g[98];
        g[103] = c0[11] * g[99];
        g[112] = cp[2] * g[96];
        g[113] = cp[5] * g[97];
        g[114] = cp[8] * g[98];
        g[115] = cp[11] * g[99];
        g[116] = cp[2] * g[100] + b00[0] * g[96];
        g[117] = cp[5] * g[101] + b00[1] * g[97];
        g[118] = cp[8] * g[102] + b00[2] * g[98];
        g[119] = cp[11] * g[103] + b00[3] * g[99];
        g[120] = g[112] * (zkzl + cp[2]) + b01[0] * g[96];
        g[121] = g[113] * (zkzl + cp[5]) + b01[1] * g[97];
        g[122] = g[114] * (zkzl + cp[8]) + b01[2] * g[98];
        g[123] = g[115] * (zkzl + cp[11]) + b01[3] * g[99];
        g[124] = g[116] * (zkzl + cp[2]) + b01[0] * g[100] + b00[0] * g[112];
        g[125] = g[117] * (zkzl + cp[5]) + b01[1] * g[101] + b00[1] * g[113];
        g[126] = g[118] * (zkzl + cp[8]) + b01[2] * g[102] + b00[2] * g[114];
        g[127] = g[119] * (zkzl + cp[11]) + b01[3] * g[103] + b00[3] * g[115];
        g[104] = g[96] * (zkzl + cp[2]);
        g[105] = g[97] * (zkzl + cp[5]);
        g[106] = g[98] * (zkzl + cp[8]);
        g[107] = g[99] * (zkzl + cp[11]);
        g[108] = g[100] * (zkzl + cp[2]) + b00[0] * g[96];
        g[109] = g[101] * (zkzl + cp[5]) + b00[1] * g[97];
        g[110] = g[102] * (zkzl + cp[8]) + b00[2] * g[98];
        g[111] = g[103] * (zkzl + cp[11]) + b00[3] * g[99];
}

static inline void _srg0_2d4d_1020(double *restrict g, struct _BC *bc, CINTEnvVars *envs)
{
        double *c0 = bc->c00;
        double *cp = bc->c0p;
        double *b00 = bc->b00;
        double *b01 = bc->b01;
        g[0] = 1;
        g[1] = 1;
        g[2] = 1;
        g[3] = 1;
        g[4] = c0[0];
        g[5] = c0[3];
        g[6] = c0[6];
        g[7] = c0[9];
        g[8] = cp[0];
        g[9] = cp[3];
        g[10] = cp[6];
        g[11] = cp[9];
        g[16] = cp[0] * cp[0] + b01[0];
        g[17] = cp[3] * cp[3] + b01[1];
        g[18] = cp[6] * cp[6] + b01[2];
        g[19] = cp[9] * cp[9] + b01[3];
        g[12] = cp[0] * c0[0] + b00[0];
        g[13] = cp[3] * c0[3] + b00[1];
        g[14] = cp[6] * c0[6] + b00[2];
        g[15] = cp[9] * c0[9] + b00[3];
        g[20] = cp[0] * (g[12] + b00[0]) + b01[0] * c0[0];
        g[21] = cp[3] * (g[13] + b00[1]) + b01[1] * c0[3];
        g[22] = cp[6] * (g[14] + b00[2]) + b01[2] * c0[6];
        g[23] = cp[9] * (g[15] + b00[3]) + b01[3] * c0[9];
        g[24] = 1;
        g[25] = 1;
        g[26] = 1;
        g[27] = 1;
        g[28] = c0[1];
        g[29] = c0[4];
        g[30] = c0[7];
        g[31] = c0[10];
        g[32] = cp[1];
        g[33] = cp[4];
        g[34] = cp[7];
        g[35] = cp[10];
        g[40] = cp[1] * cp[1] + b01[0];
        g[41] = cp[4] * cp[4] + b01[1];
        g[42] = cp[7] * cp[7] + b01[2];
        g[43] = cp[10] * cp[10] + b01[3];
        g[36] = cp[1] * c0[1] + b00[0];
        g[37] = cp[4] * c0[4] + b00[1];
        g[38] = cp[7] * c0[7] + b00[2];
        g[39] = cp[10] * c0[10] + b00[3];
        g[44] = cp[1] * (g[36] + b00[0]) + b01[0] * c0[1];
        g[45] = cp[4] * (g[37] + b00[1]) + b01[1] * c0[4];
        g[46] = cp[7] * (g[38] + b00[2]) + b01[2] * c0[7];
        g[47] = cp[10] * (g[39] + b00[3]) + b01[3] * c0[10];
        //g[48] = w[0];
        //g[49] = w[1];
        //g[50] = w[2];
        //g[51] = w[3];
        g[52] = c0[2] * g[48];
        g[53] = c0[5] * g[49];
        g[54] = c0[8] * g[50];
        g[55] = c0[11] * g[51];
        g[56] = cp[2] * g[48];
        g[57] = cp[5] * g[49];
        g[58] = cp[8] * g[50];
        g[59] = cp[11] * g[51];
        g[64] = cp[2] * g[56] + b01[0] * g[48];
        g[65] = cp[5] * g[57] + b01[1] * g[49];
        g[66] = cp[8] * g[58] + b01[2] * g[50];
        g[67] = cp[11] * g[59] + b01[3] * g[51];
        g[60] = cp[2] * g[52] + b00[0] * g[48];
        g[61] = cp[5] * g[53] + b00[1] * g[49];
        g[62] = cp[8] * g[54] + b00[2] * g[50];
        g[63] = cp[11] * g[55] + b00[3] * g[51];
        g[68] = cp[2] * g[60] + b01[0] * g[52] + b00[0] * g[56];
        g[69] = cp[5] * g[61] + b01[1] * g[53] + b00[1] * g[57];
        g[70] = cp[8] * g[62] + b01[2] * g[54] + b00[2] * g[58];
        g[71] = cp[11] * g[63] + b01[3] * g[55] + b00[3] * g[59];
}

static inline void _srg0_2d4d_1100(double *restrict g, struct _BC *bc, CINTEnvVars *envs)
{
        double *c0 = bc->c00;
        double *b10 = bc->b10;
        double xixj = envs->rirj[0];
        double yiyj = envs->rirj[1];
        double zizj = envs->rirj[2];
        g[0] = 1;
        g[1] = 1;
        g[2] = 1;
        g[3] = 1;
        g[8] = c0[0];
        g[9] = c0[3];
        g[10] = c0[6];
        g[11] = c0[9];
        g[12] = c0[0] * (xixj + c0[0]) + b10[0];
        g[13] = c0[3] * (xixj + c0[3]) + b10[1];
        g[14] = c0[6] * (xixj + c0[6]) + b10[2];
        g[15] = c0[9] * (xixj + c0[9]) + b10[3];
        g[4] = xixj + c0[0];
        g[5] = xixj + c0[3];
        g[6] = xixj + c0[6];
        g[7] = xixj + c0[9];
        g[24] = 1;
        g[25] = 1;
        g[26] = 1;
        g[27] = 1;
        g[32] = c0[1];
        g[33] = c0[4];
        g[34] = c0[7];
        g[35] = c0[10];
        g[36] = c0[1] * (yiyj + c0[1]) + b10[0];
        g[37] = c0[4] * (yiyj + c0[4]) + b10[1];
        g[38] = c0[7] * (yiyj + c0[7]) + b10[2];
        g[39] = c0[10] * (yiyj + c0[10]) + b10[3];
        g[28] = yiyj + c0[1];
        g[29] = yiyj + c0[4];
        g[30] = yiyj + c0[7];
        g[31] = yiyj + c0[10];
        //g[48] = w[0];
        //g[49] = w[1];
        //g[50] = w[2];
        //g[51] = w[3];
        g[56] = c0[2] * g[48];
        g[57] = c0[5] * g[49];
        g[58] = c0[8] * g[50];
        g[59] = c0[11] * g[51];
        g[60] = g[56] * (zizj + c0[2]) + b10[0] * g[48];
        g[61] = g[57] * (zizj + c0[5]) + b10[1] * g[49];
        g[62] = g[58] * (zizj + c0[8]) + b10[2] * g[50];
        g[63] = g[59] * (zizj + c0[11]) + b10[3] * g[51];
        g[52] = g[48] * (zizj + c0[2]);
        g[53] = g[49] * (zizj + c0[5]);
        g[54] = g[50] * (zizj + c0[8]);
        g[55] = g[51] * (zizj + c0[11]);
}

static inline void _srg0_2d4d_1101(double *restrict g, struct _BC *bc, CINTEnvVars *envs)
{
        double *c0 = bc->c00;
        double *cp = bc->c0p;
        double *b00 = bc->b00;
        double *b10 = bc->b10;
        double xixj = envs->rirj[0];
        double yiyj = envs->rirj[1];
        double zizj = envs->rirj[2];
        g[0] = 1;
        g[1] = 1;
        g[2] = 1;
        g[3] = 1;
        g[16] = c0[0];
        g[17] = c0[3];
        g[18] = c0[6];
        g[19] = c0[9];
        g[8] = cp[0];
        g[9] = cp[3];
        g[10] = cp[6];
        g[11] = cp[9];
        g[24] = cp[0] * c0[0] + b00[0];
        g[25] = cp[3] * c0[3] + b00[1];
        g[26] = cp[6] * c0[6] + b00[2];
        g[27] = cp[9] * c0[9] + b00[3];
        g[20] = c0[0] * (xixj + c0[0]) + b10[0];
        g[21] = c0[3] * (xixj + c0[3]) + b10[1];
        g[22] = c0[6] * (xixj + c0[6]) + b10[2];
        g[23] = c0[9] * (xixj + c0[9]) + b10[3];
        g[4] = xixj + c0[0];
        g[5] = xixj + c0[3];
        g[6] = xixj + c0[6];
        g[7] = xixj + c0[9];
        g[28] = g[24] * (xixj + c0[0]) + c0[0] * b00[0] + b10[0] * cp[0];
        g[29] = g[25] * (xixj + c0[3]) + c0[3] * b00[1] + b10[1] * cp[3];
        g[30] = g[26] * (xixj + c0[6]) + c0[6] * b00[2] + b10[2] * cp[6];
        g[31] = g[27] * (xixj + c0[9]) + c0[9] * b00[3] + b10[3] * cp[9];
        g[12] = cp[0] * (xixj + c0[0]) + b00[0];
        g[13] = cp[3] * (xixj + c0[3]) + b00[1];
        g[14] = cp[6] * (xixj + c0[6]) + b00[2];
        g[15] = cp[9] * (xixj + c0[9]) + b00[3];
        g[48] = 1;
        g[49] = 1;
        g[50] = 1;
        g[51] = 1;
        g[64] = c0[1];
        g[65] = c0[4];
        g[66] = c0[7];
        g[67] = c0[10];
        g[56] = cp[1];
        g[57] = cp[4];
        g[58] = cp[7];
        g[59] = cp[10];
        g[72] = cp[1] * c0[1] + b00[0];
        g[73] = cp[4] * c0[4] + b00[1];
        g[74] = cp[7] * c0[7] + b00[2];
        g[75] = cp[10] * c0[10] + b00[3];
        g[68] = c0[1] * (yiyj + c0[1]) + b10[0];
        g[69] = c0[4] * (yiyj + c0[4]) + b10[1];
        g[70] = c0[7] * (yiyj + c0[7]) + b10[2];
        g[71] = c0[10] * (yiyj + c0[10]) + b10[3];
        g[52] = yiyj + c0[1];
        g[53] = yiyj + c0[4];
        g[54] = yiyj + c0[7];
        g[55] = yiyj + c0[10];
        g[76] = g[72] * (yiyj + c0[1]) + c0[1] * b00[0] + b10[0] * cp[1];
        g[77] = g[73] * (yiyj + c0[4]) + c0[4] * b00[1] + b10[1] * cp[4];
        g[78] = g[74] * (yiyj + c0[7]) + c0[7] * b00[2] + b10[2] * cp[7];
        g[79] = g[75] * (yiyj + c0[10]) + c0[10] * b00[3] + b10[3] * cp[10];
        g[60] = cp[1] * (yiyj + c0[1]) + b00[0];
        g[61] = cp[4] * (yiyj + c0[4]) + b00[1];
        g[62] = cp[7] * (yiyj + c0[7]) + b00[2];
        g[63] = cp[10] * (yiyj + c0[10]) + b00[3];
        //g[96] = w[0];
        //g[97] = w[1];
        //g[98] = w[2];
        //g[99] = w[3];
        g[112] = c0[2] * g[96];
        g[113] = c0[5] * g[97];
        g[114] = c0[8] * g[98];
        g[115] = c0[11] * g[99];
        g[104] = cp[2] * g[96];
        g[105] = cp[5] * g[97];
        g[106] = cp[8] * g[98];
        g[107] = cp[11] * g[99];
        g[120] = cp[2] * g[112] + b00[0] * g[96];
        g[121] = cp[5] * g[113] + b00[1] * g[97];
        g[122] = cp[8] * g[114] + b00[2] * g[98];
        g[123] = cp[11] * g[115] + b00[3] * g[99];
        g[116] = g[112] * (zizj + c0[2]) + b10[0] * g[96];
        g[117] = g[113] * (zizj + c0[5]) + b10[1] * g[97];
        g[118] = g[114] * (zizj + c0[8]) + b10[2] * g[98];
        g[119] = g[115] * (zizj + c0[11]) + b10[3] * g[99];
        g[100] = g[96] * (zizj + c0[2]);
        g[101] = g[97] * (zizj + c0[5]);
        g[102] = g[98] * (zizj + c0[8]);
        g[103] = g[99] * (zizj + c0[11]);
        g[124] = g[120] * (zizj + c0[2]) + b10[0] * g[104] + b00[0] * g[112];
        g[125] = g[121] * (zizj + c0[5]) + b10[1] * g[105] + b00[1] * g[113];
        g[126] = g[122] * (zizj + c0[8]) + b10[2] * g[106] + b00[2] * g[114];
        g[127] = g[123] * (zizj + c0[11]) + b10[3] * g[107] + b00[3] * g[115];
        g[108] = zizj * g[104] + cp[2] * g[112] + b00[0] * g[96];
        g[109] = zizj * g[105] + cp[5] * g[113] + b00[1] * g[97];
        g[110] = zizj * g[106] + cp[8] * g[114] + b00[2] * g[98];
        g[111] = zizj * g[107] + cp[11] * g[115] + b00[3] * g[99];
}

static inline void _srg0_2d4d_1110(double *restrict g, struct _BC *bc, CINTEnvVars *envs)
{
        double *c0 = bc->c00;
        double *cp = bc->c0p;
        double *b00 = bc->b00;
        double *b10 = bc->b10;
        double xixj = envs->rirj[0];
        double yiyj = envs->rirj[1];
        double zizj = envs->rirj[2];
        g[0] = 1;
        g[1] = 1;
        g[2] = 1;
        g[3] = 1;
        g[16] = c0[0];
        g[17] = c0[3];
        g[18] = c0[6];
        g[19] = c0[9];
        g[8] = cp[0];
        g[9] = cp[3];
        g[10] = cp[6];
        g[11] = cp[9];
        g[24] = cp[0] * c0[0] + b00[0];
        g[25] = cp[3] * c0[3] + b00[1];
        g[26] = cp[6] * c0[6] + b00[2];
        g[27] = cp[9] * c0[9] + b00[3];
        g[20] = c0[0] * (xixj + c0[0]) + b10[0];
        g[21] = c0[3] * (xixj + c0[3]) + b10[1];
        g[22] = c0[6] * (xixj + c0[6]) + b10[2];
        g[23] = c0[9] * (xixj + c0[9]) + b10[3];
        g[4] = xixj + c0[0];
        g[5] = xixj + c0[3];
        g[6] = xixj + c0[6];
        g[7] = xixj + c0[9];
        g[28] = g[24] * (xixj + c0[0]) + c0[0] * b00[0] + b10[0] * cp[0];
        g[29] = g[25] * (xixj + c0[3]) + c0[3] * b00[1] + b10[1] * cp[3];
        g[30] = g[26] * (xixj + c0[6]) + c0[6] * b00[2] + b10[2] * cp[6];
        g[31] = g[27] * (xixj + c0[9]) + c0[9] * b00[3] + b10[3] * cp[9];
        g[12] = cp[0] * (xixj + c0[0]) + b00[0];
        g[13] = cp[3] * (xixj + c0[3]) + b00[1];
        g[14] = cp[6] * (xixj + c0[6]) + b00[2];
        g[15] = cp[9] * (xixj + c0[9]) + b00[3];
        g[48] = 1;
        g[49] = 1;
        g[50] = 1;
        g[51] = 1;
        g[64] = c0[1];
        g[65] = c0[4];
        g[66] = c0[7];
        g[67] = c0[10];
        g[56] = cp[1];
        g[57] = cp[4];
        g[58] = cp[7];
        g[59] = cp[10];
        g[72] = cp[1] * c0[1] + b00[0];
        g[73] = cp[4] * c0[4] + b00[1];
        g[74] = cp[7] * c0[7] + b00[2];
        g[75] = cp[10] * c0[10] + b00[3];
        g[68] = c0[1] * (yiyj + c0[1]) + b10[0];
        g[69] = c0[4] * (yiyj + c0[4]) + b10[1];
        g[70] = c0[7] * (yiyj + c0[7]) + b10[2];
        g[71] = c0[10] * (yiyj + c0[10]) + b10[3];
        g[52] = yiyj + c0[1];
        g[53] = yiyj + c0[4];
        g[54] = yiyj + c0[7];
        g[55] = yiyj + c0[10];
        g[76] = g[72] * (yiyj + c0[1]) + c0[1] * b00[0] + b10[0] * cp[1];
        g[77] = g[73] * (yiyj + c0[4]) + c0[4] * b00[1] + b10[1] * cp[4];
        g[78] = g[74] * (yiyj + c0[7]) + c0[7] * b00[2] + b10[2] * cp[7];
        g[79] = g[75] * (yiyj + c0[10]) + c0[10] * b00[3] + b10[3] * cp[10];
        g[60] = cp[1] * (yiyj + c0[1]) + b00[0];
        g[61] = cp[4] * (yiyj + c0[4]) + b00[1];
        g[62] = cp[7] * (yiyj + c0[7]) + b00[2];
        g[63] = cp[10] * (yiyj + c0[10]) + b00[3];
        //g[96] = w[0];
        //g[97] = w[1];
        //g[98] = w[2];
        //g[99] = w[3];
        g[112] = c0[2] * g[96];
        g[113] = c0[5] * g[97];
        g[114] = c0[8] * g[98];
        g[115] = c0[11] * g[99];
        g[104] = cp[2] * g[96];
        g[105] = cp[5] * g[97];
        g[106] = cp[8] * g[98];
        g[107] = cp[11] * g[99];
        g[120] = cp[2] * g[112] + b00[0] * g[96];
        g[121] = cp[5] * g[113] + b00[1] * g[97];
        g[122] = cp[8] * g[114] + b00[2] * g[98];
        g[123] = cp[11] * g[115] + b00[3] * g[99];
        g[116] = g[112] * (zizj + c0[2]) + b10[0] * g[96];
        g[117] = g[113] * (zizj + c0[5]) + b10[1] * g[97];
        g[118] = g[114] * (zizj + c0[8]) + b10[2] * g[98];
        g[119] = g[115] * (zizj + c0[11]) + b10[3] * g[99];
        g[100] = g[96] * (zizj + c0[2]);
        g[101] = g[97] * (zizj + c0[5]);
        g[102] = g[98] * (zizj + c0[8]);
        g[103] = g[99] * (zizj + c0[11]);
        g[124] = g[120] * (zizj + c0[2]) + b10[0] * g[104] + b00[0] * g[112];
        g[125] = g[121] * (zizj + c0[5]) + b10[1] * g[105] + b00[1] * g[113];
        g[126] = g[122] * (zizj + c0[8]) + b10[2] * g[106] + b00[2] * g[114];
        g[127] = g[123] * (zizj + c0[11]) + b10[3] * g[107] + b00[3] * g[115];
        g[108] = zizj * g[104] + cp[2] * g[112] + b00[0] * g[96];
        g[109] = zizj * g[105] + cp[5] * g[113] + b00[1] * g[97];
        g[110] = zizj * g[106] + cp[8] * g[114] + b00[2] * g[98];
        g[111] = zizj * g[107] + cp[11] * g[115] + b00[3] * g[99];
}

static inline void _srg0_2d4d_1200(double *restrict g, struct _BC *bc, CINTEnvVars *envs)
{
        double *c0 = bc->c00;
        double *b10 = bc->b10;
        double xixj = envs->rirj[0];
        double yiyj = envs->rirj[1];
        double zizj = envs->rirj[2];
        g[0] = 1;
        g[1] = 1;
        g[2] = 1;
        g[3] = 1;
        g[8] = c0[0];
        g[9] = c0[3];
        g[10] = c0[6];
        g[11] = c0[9];
        g[16] = c0[0] * c0[0] + b10[0];
        g[17] = c0[3] * c0[3] + b10[1];
        g[18] = c0[6] * c0[6] + b10[2];
        g[19] = c0[9] * c0[9] + b10[3];
        g[20] = g[16] * (xixj + c0[0]) + c0[0] * 2 * b10[0];
        g[21] = g[17] * (xixj + c0[3]) + c0[3] * 2 * b10[1];
        g[22] = g[18] * (xixj + c0[6]) + c0[6] * 2 * b10[2];
        g[23] = g[19] * (xixj + c0[9]) + c0[9] * 2 * b10[3];
        g[12] = c0[0] * (xixj + c0[0]) + b10[0];
        g[13] = c0[3] * (xixj + c0[3]) + b10[1];
        g[14] = c0[6] * (xixj + c0[6]) + b10[2];
        g[15] = c0[9] * (xixj + c0[9]) + b10[3];
        g[4] = xixj + c0[0];
        g[5] = xixj + c0[3];
        g[6] = xixj + c0[6];
        g[7] = xixj + c0[9];
        g[32] = 1;
        g[33] = 1;
        g[34] = 1;
        g[35] = 1;
        g[40] = c0[1];
        g[41] = c0[4];
        g[42] = c0[7];
        g[43] = c0[10];
        g[48] = c0[1] * c0[1] + b10[0];
        g[49] = c0[4] * c0[4] + b10[1];
        g[50] = c0[7] * c0[7] + b10[2];
        g[51] = c0[10] * c0[10] + b10[3];
        g[52] = g[48] * (yiyj + c0[1]) + c0[1] * 2 * b10[0];
        g[53] = g[49] * (yiyj + c0[4]) + c0[4] * 2 * b10[1];
        g[54] = g[50] * (yiyj + c0[7]) + c0[7] * 2 * b10[2];
        g[55] = g[51] * (yiyj + c0[10]) + c0[10] * 2 * b10[3];
        g[44] = c0[1] * (yiyj + c0[1]) + b10[0];
        g[45] = c0[4] * (yiyj + c0[4]) + b10[1];
        g[46] = c0[7] * (yiyj + c0[7]) + b10[2];
        g[47] = c0[10] * (yiyj + c0[10]) + b10[3];
        g[36] = yiyj + c0[1];
        g[37] = yiyj + c0[4];
        g[38] = yiyj + c0[7];
        g[39] = yiyj + c0[10];
        //g[64] = w[0];
        //g[65] = w[1];
        //g[66] = w[2];
        //g[67] = w[3];
        g[72] = c0[2] * g[64];
        g[73] = c0[5] * g[65];
        g[74] = c0[8] * g[66];
        g[75] = c0[11] * g[67];
        g[80] = c0[2] * g[72] + b10[0] * g[64];
        g[81] = c0[5] * g[73] + b10[1] * g[65];
        g[82] = c0[8] * g[74] + b10[2] * g[66];
        g[83] = c0[11] * g[75] + b10[3] * g[67];
        g[84] = g[80] * (zizj + c0[2]) + 2 * b10[0] * g[72];
        g[85] = g[81] * (zizj + c0[5]) + 2 * b10[1] * g[73];
        g[86] = g[82] * (zizj + c0[8]) + 2 * b10[2] * g[74];
        g[87] = g[83] * (zizj + c0[11]) + 2 * b10[3] * g[75];
        g[76] = g[72] * (zizj + c0[2]) + b10[0] * g[64];
        g[77] = g[73] * (zizj + c0[5]) + b10[1] * g[65];
        g[78] = g[74] * (zizj + c0[8]) + b10[2] * g[66];
        g[79] = g[75] * (zizj + c0[11]) + b10[3] * g[67];
        g[68] = g[64] * (zizj + c0[2]);
        g[69] = g[65] * (zizj + c0[5]);
        g[70] = g[66] * (zizj + c0[8]);
        g[71] = g[67] * (zizj + c0[11]);
}

static inline void _srg0_2d4d_2000(double *restrict g, struct _BC *bc, CINTEnvVars *envs)
{
        double *c0 = bc->c00;
        double *b10 = bc->b10;
        g[0] = 1;
        g[1] = 1;
        g[2] = 1;
        g[3] = 1;
        g[4] = c0[0];
        g[5] = c0[3];
        g[6] = c0[6];
        g[7] = c0[9];
        g[8] = c0[0] * c0[0] + b10[0];
        g[9] = c0[3] * c0[3] + b10[1];
        g[10] = c0[6] * c0[6] + b10[2];
        g[11] = c0[9] * c0[9] + b10[3];
        g[12] = 1;
        g[13] = 1;
        g[14] = 1;
        g[15] = 1;
        g[16] = c0[1];
        g[17] = c0[4];
        g[18] = c0[7];
        g[19] = c0[10];
        g[20] = c0[1] * c0[1] + b10[0];
        g[21] = c0[4] * c0[4] + b10[1];
        g[22] = c0[7] * c0[7] + b10[2];
        g[23] = c0[10] * c0[10] + b10[3];
        //g[24] = w[0];
        //g[25] = w[1];
        //g[26] = w[2];
        //g[27] = w[3];
        g[28] = c0[2] * g[24];
        g[29] = c0[5] * g[25];
        g[30] = c0[8] * g[26];
        g[31] = c0[11] * g[27];
        g[32] = c0[2] * g[28] + b10[0] * g[24];
        g[33] = c0[5] * g[29] + b10[1] * g[25];
        g[34] = c0[8] * g[30] + b10[2] * g[26];
        g[35] = c0[11] * g[31] + b10[3] * g[27];
}

static inline void _srg0_2d4d_2001(double *restrict g, struct _BC *bc, CINTEnvVars *envs)
{
        double *c0 = bc->c00;
        double *cp = bc->c0p;
        double *b00 = bc->b00;
        double *b10 = bc->b10;
        g[0] = 1;
        g[1] = 1;
        g[2] = 1;
        g[3] = 1;
        g[4] = c0[0];
        g[5] = c0[3];
        g[6] = c0[6];
        g[7] = c0[9];
        g[8] = c0[0] * c0[0] + b10[0];
        g[9] = c0[3] * c0[3] + b10[1];
        g[10] = c0[6] * c0[6] + b10[2];
        g[11] = c0[9] * c0[9] + b10[3];
        g[12] = cp[0];
        g[13] = cp[3];
        g[14] = cp[6];
        g[15] = cp[9];
        g[16] = cp[0] * c0[0] + b00[0];
        g[17] = cp[3] * c0[3] + b00[1];
        g[18] = cp[6] * c0[6] + b00[2];
        g[19] = cp[9] * c0[9] + b00[3];
        g[20] = c0[0] * (g[16] + b00[0]) + b10[0] * cp[0];
        g[21] = c0[3] * (g[17] + b00[1]) + b10[1] * cp[3];
        g[22] = c0[6] * (g[18] + b00[2]) + b10[2] * cp[6];
        g[23] = c0[9] * (g[19] + b00[3]) + b10[3] * cp[9];
        g[24] = 1;
        g[25] = 1;
        g[26] = 1;
        g[27] = 1;
        g[28] = c0[1];
        g[29] = c0[4];
        g[30] = c0[7];
        g[31] = c0[10];
        g[32] = c0[1] * c0[1] + b10[0];
        g[33] = c0[4] * c0[4] + b10[1];
        g[34] = c0[7] * c0[7] + b10[2];
        g[35] = c0[10] * c0[10] + b10[3];
        g[36] = cp[1];
        g[37] = cp[4];
        g[38] = cp[7];
        g[39] = cp[10];
        g[40] = cp[1] * c0[1] + b00[0];
        g[41] = cp[4] * c0[4] + b00[1];
        g[42] = cp[7] * c0[7] + b00[2];
        g[43] = cp[10] * c0[10] + b00[3];
        g[44] = c0[1] * (g[40] + b00[0]) + b10[0] * cp[1];
        g[45] = c0[4] * (g[41] + b00[1]) + b10[1] * cp[4];
        g[46] = c0[7] * (g[42] + b00[2]) + b10[2] * cp[7];
        g[47] = c0[10] * (g[43] + b00[3]) + b10[3] * cp[10];
        //g[48] = w[0];
        //g[49] = w[1];
        //g[50] = w[2];
        //g[51] = w[3];
        g[52] = c0[2] * g[48];
        g[53] = c0[5] * g[49];
        g[54] = c0[8] * g[50];
        g[55] = c0[11] * g[51];
        g[56] = c0[2] * g[52] + b10[0] * g[48];
        g[57] = c0[5] * g[53] + b10[1] * g[49];
        g[58] = c0[8] * g[54] + b10[2] * g[50];
        g[59] = c0[11] * g[55] + b10[3] * g[51];
        g[60] = cp[2] * g[48];
        g[61] = cp[5] * g[49];
        g[62] = cp[8] * g[50];
        g[63] = cp[11] * g[51];
        g[64] = cp[2] * g[52] + b00[0] * g[48];
        g[65] = cp[5] * g[53] + b00[1] * g[49];
        g[66] = cp[8] * g[54] + b00[2] * g[50];
        g[67] = cp[11] * g[55] + b00[3] * g[51];
        g[68] = c0[2] * g[64] + b10[0] * g[60] + b00[0] * g[52];
        g[69] = c0[5] * g[65] + b10[1] * g[61] + b00[1] * g[53];
        g[70] = c0[8] * g[66] + b10[2] * g[62] + b00[2] * g[54];
        g[71] = c0[11] * g[67] + b10[3] * g[63] + b00[3] * g[55];
}

static inline void _srg0_2d4d_2010(double *restrict g, struct _BC *bc, CINTEnvVars *envs)
{
        double *c0 = bc->c00;
        double *cp = bc->c0p;
        double *b00 = bc->b00;
        double *b10 = bc->b10;
        g[0] = 1;
        g[1] = 1;
        g[2] = 1;
        g[3] = 1;
        g[4] = c0[0];
        g[5] = c0[3];
        g[6] = c0[6];
        g[7] = c0[9];
        g[8] = c0[0] * c0[0] + b10[0];
        g[9] = c0[3] * c0[3] + b10[1];
        g[10] = c0[6] * c0[6] + b10[2];
        g[11] = c0[9] * c0[9] + b10[3];
        g[12] = cp[0];
        g[13] = cp[3];
        g[14] = cp[6];
        g[15] = cp[9];
        g[16] = cp[0] * c0[0] + b00[0];
        g[17] = cp[3] * c0[3] + b00[1];
        g[18] = cp[6] * c0[6] + b00[2];
        g[19] = cp[9] * c0[9] + b00[3];
        g[20] = c0[0] * (g[16] + b00[0]) + b10[0] * cp[0];
        g[21] = c0[3] * (g[17] + b00[1]) + b10[1] * cp[3];
        g[22] = c0[6] * (g[18] + b00[2]) + b10[2] * cp[6];
        g[23] = c0[9] * (g[19] + b00[3]) + b10[3] * cp[9];
        g[24] = 1;
        g[25] = 1;
        g[26] = 1;
        g[27] = 1;
        g[28] = c0[1];
        g[29] = c0[4];
        g[30] = c0[7];
        g[31] = c0[10];
        g[32] = c0[1] * c0[1] + b10[0];
        g[33] = c0[4] * c0[4] + b10[1];
        g[34] = c0[7] * c0[7] + b10[2];
        g[35] = c0[10] * c0[10] + b10[3];
        g[36] = cp[1];
        g[37] = cp[4];
        g[38] = cp[7];
        g[39] = cp[10];
        g[40] = cp[1] * c0[1] + b00[0];
        g[41] = cp[4] * c0[4] + b00[1];
        g[42] = cp[7] * c0[7] + b00[2];
        g[43] = cp[10] * c0[10] + b00[3];
        g[44] = c0[1] * (g[40] + b00[0]) + b10[0] * cp[1];
        g[45] = c0[4] * (g[41] + b00[1]) + b10[1] * cp[4];
        g[46] = c0[7] * (g[42] + b00[2]) + b10[2] * cp[7];
        g[47] = c0[10] * (g[43] + b00[3]) + b10[3] * cp[10];
        //g[48] = w[0];
        //g[49] = w[1];
        //g[50] = w[2];
        //g[51] = w[3];
        g[52] = c0[2] * g[48];
        g[53] = c0[5] * g[49];
        g[54] = c0[8] * g[50];
        g[55] = c0[11] * g[51];
        g[56] = c0[2] * g[52] + b10[0] * g[48];
        g[57] = c0[5] * g[53] + b10[1] * g[49];
        g[58] = c0[8] * g[54] + b10[2] * g[50];
        g[59] = c0[11] * g[55] + b10[3] * g[51];
        g[60] = cp[2] * g[48];
        g[61] = cp[5] * g[49];
        g[62] = cp[8] * g[50];
        g[63] = cp[11] * g[51];
        g[64] = cp[2] * g[52] + b00[0] * g[48];
        g[65] = cp[5] * g[53] + b00[1] * g[49];
        g[66] = cp[8] * g[54] + b00[2] * g[50];
        g[67] = cp[11] * g[55] + b00[3] * g[51];
        g[68] = c0[2] * g[64] + b10[0] * g[60] + b00[0] * g[52];
        g[69] = c0[5] * g[65] + b10[1] * g[61] + b00[1] * g[53];
        g[70] = c0[8] * g[66] + b10[2] * g[62] + b00[2] * g[54];
        g[71] = c0[11] * g[67] + b10[3] * g[63] + b00[3] * g[55];
}

static inline void _srg0_2d4d_2100(double *restrict g, struct _BC *bc, CINTEnvVars *envs)
{
        double *c0 = bc->c00;
        double *b10 = bc->b10;
        double xixj = envs->rirj[0];
        double yiyj = envs->rirj[1];
        double zizj = envs->rirj[2];
        g[0] = 1;
        g[1] = 1;
        g[2] = 1;
        g[3] = 1;
        g[4] = c0[0];
        g[5] = c0[3];
        g[6] = c0[6];
        g[7] = c0[9];
        g[8] = c0[0] * c0[0] + b10[0];
        g[9] = c0[3] * c0[3] + b10[1];
        g[10] = c0[6] * c0[6] + b10[2];
        g[11] = c0[9] * c0[9] + b10[3];
        g[24] = g[8] * (xixj + c0[0]) + c0[0] * 2 * b10[0];
        g[25] = g[9] * (xixj + c0[3]) + c0[3] * 2 * b10[1];
        g[26] = g[10] * (xixj + c0[6]) + c0[6] * 2 * b10[2];
        g[27] = g[11] * (xixj + c0[9]) + c0[9] * 2 * b10[3];
        g[20] = c0[0] * (xixj + c0[0]) + b10[0];
        g[21] = c0[3] * (xixj + c0[3]) + b10[1];
        g[22] = c0[6] * (xixj + c0[6]) + b10[2];
        g[23] = c0[9] * (xixj + c0[9]) + b10[3];
        g[16] = xixj + c0[0];
        g[17] = xixj + c0[3];
        g[18] = xixj + c0[6];
        g[19] = xixj + c0[9];
        g[32] = 1;
        g[33] = 1;
        g[34] = 1;
        g[35] = 1;
        g[36] = c0[1];
        g[37] = c0[4];
        g[38] = c0[7];
        g[39] = c0[10];
        g[40] = c0[1] * c0[1] + b10[0];
        g[41] = c0[4] * c0[4] + b10[1];
        g[42] = c0[7] * c0[7] + b10[2];
        g[43] = c0[10] * c0[10] + b10[3];
        g[56] = g[40] * (yiyj + c0[1]) + c0[1] * 2 * b10[0];
        g[57] = g[41] * (yiyj + c0[4]) + c0[4] * 2 * b10[1];
        g[58] = g[42] * (yiyj + c0[7]) + c0[7] * 2 * b10[2];
        g[59] = g[43] * (yiyj + c0[10]) + c0[10] * 2 * b10[3];
        g[52] = c0[1] * (yiyj + c0[1]) + b10[0];
        g[53] = c0[4] * (yiyj + c0[4]) + b10[1];
        g[54] = c0[7] * (yiyj + c0[7]) + b10[2];
        g[55] = c0[10] * (yiyj + c0[10]) + b10[3];
        g[48] = yiyj + c0[1];
        g[49] = yiyj + c0[4];
        g[50] = yiyj + c0[7];
        g[51] = yiyj + c0[10];
        //g[64] = w[0];
        //g[65] = w[1];
        //g[66] = w[2];
        //g[67] = w[3];
        g[68] = c0[2] * g[64];
        g[69] = c0[5] * g[65];
        g[70] = c0[8] * g[66];
        g[71] = c0[11] * g[67];
        g[72] = c0[2] * g[68] + b10[0] * g[64];
        g[73] = c0[5] * g[69] + b10[1] * g[65];
        g[74] = c0[8] * g[70] + b10[2] * g[66];
        g[75] = c0[11] * g[71] + b10[3] * g[67];
        g[88] = g[72] * (zizj + c0[2]) + 2 * b10[0] * g[68];
        g[89] = g[73] * (zizj + c0[5]) + 2 * b10[1] * g[69];
        g[90] = g[74] * (zizj + c0[8]) + 2 * b10[2] * g[70];
        g[91] = g[75] * (zizj + c0[11]) + 2 * b10[3] * g[71];
        g[84] = g[68] * (zizj + c0[2]) + b10[0] * g[64];
        g[85] = g[69] * (zizj + c0[5]) + b10[1] * g[65];
        g[86] = g[70] * (zizj + c0[8]) + b10[2] * g[66];
        g[87] = g[71] * (zizj + c0[11]) + b10[3] * g[67];
        g[80] = g[64] * (zizj + c0[2]);
        g[81] = g[65] * (zizj + c0[5]);
        g[82] = g[66] * (zizj + c0[8]);
        g[83] = g[67] * (zizj + c0[11]);
}

static inline void _srg0_2d4d_3000(double *restrict g, struct _BC *bc, CINTEnvVars *envs)
{
        double *c0 = bc->c00;
        double *b10 = bc->b10;
        g[0] = 1;
        g[1] = 1;
        g[2] = 1;
        g[3] = 1;
        g[4] = c0[0];
        g[5] = c0[3];
        g[6] = c0[6];
        g[7] = c0[9];
        g[8] = c0[0] * c0[0] + b10[0];
        g[9] = c0[3] * c0[3] + b10[1];
        g[10] = c0[6] * c0[6] + b10[2];
        g[11] = c0[9] * c0[9] + b10[3];
        g[12] = c0[0] * (g[8] + 2 * b10[0]);
        g[13] = c0[3] * (g[9] + 2 * b10[1]);
        g[14] = c0[6] * (g[10] + 2 * b10[2]);
        g[15] = c0[9] * (g[11] + 2 * b10[3]);
        g[16] = 1;
        g[17] = 1;
        g[18] = 1;
        g[19] = 1;
        g[20] = c0[1];
        g[21] = c0[4];
        g[22] = c0[7];
        g[23] = c0[10];
        g[24] = c0[1] * c0[1] + b10[0];
        g[25] = c0[4] * c0[4] + b10[1];
        g[26] = c0[7] * c0[7] + b10[2];
        g[27] = c0[10] * c0[10] + b10[3];
        g[28] = c0[1] * (g[24] + 2 * b10[0]);
        g[29] = c0[4] * (g[25] + 2 * b10[1]);
        g[30] = c0[7] * (g[26] + 2 * b10[2]);
        g[31] = c0[10] * (g[27] + 2 * b10[3]);
        //g[32] = w[0];
        //g[33] = w[1];
        //g[34] = w[2];
        //g[35] = w[3];
        g[36] = c0[2] * g[32];
        g[37] = c0[5] * g[33];
        g[38] = c0[8] * g[34];
        g[39] = c0[11] * g[35];
        g[40] = c0[2] * g[36] + b10[0] * g[32];
        g[41] = c0[5] * g[37] + b10[1] * g[33];
        g[42] = c0[8] * g[38] + b10[2] * g[34];
        g[43] = c0[11] * g[39] + b10[3] * g[35];
        g[44] = c0[2] * g[40] + 2 * b10[0] * g[36];
        g[45] = c0[5] * g[41] + 2 * b10[1] * g[37];
        g[46] = c0[8] * g[42] + 2 * b10[2] * g[38];
        g[47] = c0[11] * g[43] + 2 * b10[3] * g[39];
}

void CINTsrg0_2e_2d4d_unrolled(double *restrict g, struct _BC *bc, CINTEnvVars *envs)
{
        int type_ijkl = ((envs->li_ceil << 6) | (envs->lj_ceil << 4) |
                         (envs->lk_ceil << 2) | (envs->ll_ceil));
        switch (type_ijkl) {
                case 0b00000000: _srg0_2d4d_0000(g, bc, envs); return;
                case 0b00000001: _srg0_2d4d_0001(g, bc, envs); return;
                case 0b00000010: _srg0_2d4d_0002(g, bc, envs); return;
                case 0b00000011: _srg0_2d4d_0003(g, bc, envs); return;
                case 0b00000100: _srg0_2d4d_0010(g, bc, envs); return;
                case 0b00000101: _srg0_2d4d_0011(g, bc, envs); return;
                case 0b00000110: _srg0_2d4d_0012(g, bc, envs); return;
                case 0b00001000: _srg0_2d4d_0020(g, bc, envs); return;
                case 0b00001001: _srg0_2d4d_0021(g, bc, envs); return;
                case 0b00001100: _srg0_2d4d_0030(g, bc, envs); return;
                case 0b00010000: _srg0_2d4d_0100(g, bc, envs); return;
                case 0b00010001: _srg0_2d4d_0101(g, bc, envs); return;
                case 0b00010010: _srg0_2d4d_0102(g, bc, envs); return;
                case 0b00010100: _srg0_2d4d_0110(g, bc, envs); return;
                case 0b00010101: _srg0_2d4d_0111(g, bc, envs); return;
                case 0b00011000: _srg0_2d4d_0120(g, bc, envs); return;
                case 0b00100000: _srg0_2d4d_0200(g, bc, envs); return;
                case 0b00100001: _srg0_2d4d_0201(g, bc, envs); return;
                case 0b00100100: _srg0_2d4d_0210(g, bc, envs); return;
                case 0b00110000: _srg0_2d4d_0300(g, bc, envs); return;
                case 0b01000000: _srg0_2d4d_1000(g, bc, envs); return;
                case 0b01000001: _srg0_2d4d_1001(g, bc, envs); return;
                case 0b01000010: _srg0_2d4d_1002(g, bc, envs); return;
                case 0b01000100: _srg0_2d4d_1010(g, bc, envs); return;
                case 0b01000101: _srg0_2d4d_1011(g, bc, envs); return;
                case 0b01001000: _srg0_2d4d_1020(g, bc, envs); return;
                case 0b01010000: _srg0_2d4d_1100(g, bc, envs); return;
                case 0b01010001: _srg0_2d4d_1101(g, bc, envs); return;
                case 0b01010100: _srg0_2d4d_1110(g, bc, envs); return;
                case 0b01100000: _srg0_2d4d_1200(g, bc, envs); return;
                case 0b10000000: _srg0_2d4d_2000(g, bc, envs); return;
                case 0b10000001: _srg0_2d4d_2001(g, bc, envs); return;
                case 0b10000100: _srg0_2d4d_2010(g, bc, envs); return;
                case 0b10010000: _srg0_2d4d_2100(g, bc, envs); return;
                case 0b11000000: _srg0_2d4d_3000(g, bc, envs); return;
        }
        fprintf(stderr, "Dimension error for CINTg0_2e_lj2d4d: iklj = %d %d %d %d",
               (int)envs->li_ceil, (int)envs->lk_ceil,
               (int)envs->ll_ceil, (int)envs->lj_ceil);
}

void CINTg0_2e_lj2d4d(double *restrict g, struct _BC *bc, CINTEnvVars *envs)
{
        CINTg0_2e_2d(g, bc, envs);
        CINTg0_lj2d_4d(g, envs);
}

void CINTg0_2e_kj2d4d(double *restrict g, struct _BC *bc, CINTEnvVars *envs)
{
        CINTg0_2e_2d(g, bc, envs);
        CINTg0_kj2d_4d(g, envs);
}
void CINTg0_2e_ik2d4d(double *restrict g, struct _BC *bc, CINTEnvVars *envs)
{
        CINTg0_2e_2d(g, bc, envs);
        CINTg0_ik2d_4d(g, envs);
}
void CINTg0_2e_il2d4d(double *restrict g, struct _BC *bc, CINTEnvVars *envs)
{
        CINTg0_2e_2d(g, bc, envs);
        CINTg0_il2d_4d(g, envs);
}

/*
 * g[i,k,l,j] = < ik | lj > = ( i j | k l )
 */
FINT CINTg0_2e(double *g, double *rij, double *rkl, double cutoff, CINTEnvVars *envs)
{
        FINT irys;
        FINT nroots = envs->nrys_roots;
        double aij = envs->ai[0] + envs->aj[0];
        double akl = envs->ak[0] + envs->al[0];
        double a0, a1, fac1, x;
        double u[MXRYSROOTS];
        double *w = g + envs->g_size * 2; // ~ gz
        double xij_kl = rij[0] - rkl[0];
        double yij_kl = rij[1] - rkl[1];
        double zij_kl = rij[2] - rkl[2];
        double rr = xij_kl * xij_kl + yij_kl * yij_kl + zij_kl * zij_kl;

        a1 = aij * akl;
        a0 = a1 / (aij + akl);

#ifdef WITH_RANGE_COULOMB
        const double omega = envs->env[PTR_RANGE_OMEGA];
        double theta = 0;
        if (omega == 0.) {
                x = a0 * rr;
                CINTrys_roots(nroots, x, u, w);
        } else if (omega < 0.) {
                // short-range part of range-separated Coulomb
                theta = omega * omega / (omega * omega + a0);
                x = a0 * rr;
                // very small erfc() leads to ~0 weights. They can cause
                // numerical issue in sr_rys_roots
                if (theta * x > cutoff || theta * x > EXPCUTOFF_SR) {
                        return 0;
                }
                int rorder = envs->rys_order;
                if (rorder == nroots) {
                        CINTsr_rys_roots(nroots, x, sqrt(theta), u, w);
                } else {
                        double theta1 = -sqrt(theta);
                        CINTrys_roots(rorder, x, u, w);
                        CINTrys_roots(rorder, theta*x, u+rorder, w+rorder);
                        if (envs->g_size == 2) {
                                fac1 = sqrt(a0 / (a1 * a1 * a1)) * envs->fac[0];
                                g[0] = 1;
                                g[1] = 1;
                                g[2] = 1;
                                g[3] = 1;
                                g[4] *= fac1;
                                g[5] *= fac1 * theta1;
                                return 1;
                        }
                        for (irys = rorder; irys < nroots; irys++) {
                                double ut = u[irys] * theta;
                                u[irys] = ut / (u[irys]+1.-ut);
                                w[irys] *= theta1;
                        }
                }
        } else { // omega > 0.
                // long-range part of range-separated Coulomb
                theta = omega * omega / (omega * omega + a0);
                a0 *= theta;
                x = a0 * rr;
                CINTrys_roots(nroots, x, u, w);
                /* u[:] = tau^2 / (1 - tau^2)
                 * omega^2u^2 = a0 * tau^2 / (theta^-1 - tau^2)
                 * transform u[:] to theta^-1 tau^2 / (theta^-1 - tau^2)
                 * so the rest code can be reused.
                 */
                for (irys = 0; irys < nroots; irys++) {
                        u[irys] /= u[irys] + 1 - u[irys] * theta;
                }
        }
#else
        x = a0 * rr;
        CINTrys_roots(nroots, x, u, w);
#endif
        fac1 = sqrt(a0 / (a1 * a1 * a1)) * envs->fac[0];
        if (envs->g_size == 1) {
                g[0] = 1;
                g[1] = 1;
                g[2] *= fac1;
                return 1;
        }

        double u2, tmp1, tmp2, tmp3, tmp4, tmp5;
        double rijrx[3];
        double rklrx[3];
        rijrx[0] = rij[0] - envs->rx_in_rijrx[0];
        rijrx[1] = rij[1] - envs->rx_in_rijrx[1];
        rijrx[2] = rij[2] - envs->rx_in_rijrx[2];
        rklrx[0] = rkl[0] - envs->rx_in_rklrx[0];
        rklrx[1] = rkl[1] - envs->rx_in_rklrx[1];
        rklrx[2] = rkl[2] - envs->rx_in_rklrx[2];
        struct _BC bc;
        double *c00 = bc.c00;
        double *c0p = bc.c0p;
        double *b00 = bc.b00;
        double *b10 = bc.b10;
        double *b01 = bc.b01;

        for (irys = 0; irys < nroots; irys++, c00+=3, c0p+=3) {
                /*
                 *u(irys) = t2/(1-t2)
                 *t2 = u(irys)/(1+u(irys))
                 *u2 = aij*akl/(aij+akl)*t2/(1-t2)
                 */
                u2 = a0 * u[irys];
                tmp4 = .5 / (u2 * (aij + akl) + a1);
                tmp5 = u2 * tmp4;
                tmp1 = 2. * tmp5;
                tmp2 = tmp1 * akl;
                tmp3 = tmp1 * aij;
                b00[irys] = tmp5;
                b10[irys] = tmp5 + tmp4 * akl;
                b01[irys] = tmp5 + tmp4 * aij;
                c00[0] = rijrx[0] - tmp2 * xij_kl;
                c00[1] = rijrx[1] - tmp2 * yij_kl;
                c00[2] = rijrx[2] - tmp2 * zij_kl;
                c0p[0] = rklrx[0] + tmp3 * xij_kl;
                c0p[1] = rklrx[1] + tmp3 * yij_kl;
                c0p[2] = rklrx[2] + tmp3 * zij_kl;
                w[irys] *= fac1;
        }

        (*envs->f_g0_2d4d)(g, &bc, envs);

        return 1;
}

/*
 * ( \nabla i j | kl )
 */
void CINTnabla1i_2e(double *f, const double *g,
                    const FINT li, const FINT lj, const FINT lk, const FINT ll,
                    const CINTEnvVars *envs)
{
        FINT i, j, k, l, n, ptr;
        const FINT di = envs->g_stride_i;
        const FINT dk = envs->g_stride_k;
        const FINT dl = envs->g_stride_l;
        const FINT dj = envs->g_stride_j;
        const FINT nroots = envs->nrys_roots;
        const double ai2 = -2 * envs->ai[0];
        DEF_GXYZ(const double, g, gx, gy, gz);
        DEF_GXYZ(double, f, fx, fy, fz);

        const double *p1x = gx - di;
        const double *p1y = gy - di;
        const double *p1z = gz - di;
        const double *p2x = gx + di;
        const double *p2y = gy + di;
        const double *p2z = gz + di;
        for (j = 0; j <= lj; j++)
        for (l = 0; l <= ll; l++)
        for (k = 0; k <= lk; k++) {
                ptr = dj * j + dl * l + dk * k;
                //f(...,0,...) = -2*ai*g(...,1,...)
                for (n = ptr; n < ptr+nroots; n++) {
                        fx[n] = ai2 * p2x[n];
                        fy[n] = ai2 * p2y[n];
                        fz[n] = ai2 * p2z[n];
                }
                ptr += di;
                //f(...,i,...) = i*g(...,i-1,...)-2*ai*g(...,i+1,...)
                for (i = 1; i <= li; i++) {
                        for (n = ptr; n < ptr+nroots; n++) {
                                fx[n] = i*p1x[n] + ai2*p2x[n];
                                fy[n] = i*p1y[n] + ai2*p2y[n];
                                fz[n] = i*p1z[n] + ai2*p2z[n];
                        }
                        ptr += di;
                }
        }
}


/*
 * ( i \nabla j | kl )
 */
void CINTnabla1j_2e(double *f, const double *g,
                    const FINT li, const FINT lj, const FINT lk, const FINT ll,
                    const CINTEnvVars *envs)
{
        FINT i, j, k, l, n, ptr;
        const FINT di = envs->g_stride_i;
        const FINT dk = envs->g_stride_k;
        const FINT dl = envs->g_stride_l;
        const FINT dj = envs->g_stride_j;
        const FINT nroots = envs->nrys_roots;
        const double aj2 = -2 * envs->aj[0];
        DEF_GXYZ(const double, g, gx, gy, gz);
        DEF_GXYZ(double, f, fx, fy, fz);

        const double *p1x = gx - dj;
        const double *p1y = gy - dj;
        const double *p1z = gz - dj;
        const double *p2x = gx + dj;
        const double *p2y = gy + dj;
        const double *p2z = gz + dj;
        //f(...,0,...) = -2*aj*g(...,1,...)
        for (l = 0; l <= ll; l++) {
        for (k = 0; k <= lk; k++) {
                ptr = dl * l + dk * k;
                for (i = 0; i <= li; i++) {
                        for (n = ptr; n < ptr+nroots; n++) {
                                fx[n] = aj2 * p2x[n];
                                fy[n] = aj2 * p2y[n];
                                fz[n] = aj2 * p2z[n];
                        }
                        ptr += di;
                }
        } }
        //f(...,j,...) = j*g(...,j-1,...)-2*aj*g(...,j+1,...)
        for (j = 1; j <= lj; j++) {
                for (l = 0; l <= ll; l++) {
                for (k = 0; k <= lk; k++) {
                        ptr = dj * j + dl * l + dk * k;
                        for (i = 0; i <= li; i++) {
                                for (n = ptr; n < ptr+nroots; n++) {
                                        fx[n] = j*p1x[n] + aj2*p2x[n];
                                        fy[n] = j*p1y[n] + aj2*p2y[n];
                                        fz[n] = j*p1z[n] + aj2*p2z[n];
                                }
                                ptr += di;
                        }
                } }
        }
}


/*
 * ( ij | \nabla k l )
 */
void CINTnabla1k_2e(double *f, const double *g,
                    const FINT li, const FINT lj, const FINT lk, const FINT ll,
                    const CINTEnvVars *envs)
{
        FINT i, j, k, l, n, ptr;
        const FINT di = envs->g_stride_i;
        const FINT dk = envs->g_stride_k;
        const FINT dl = envs->g_stride_l;
        const FINT dj = envs->g_stride_j;
        const FINT nroots = envs->nrys_roots;
        const double ak2 = -2 * envs->ak[0];
        DEF_GXYZ(const double, g, gx, gy, gz);
        DEF_GXYZ(double, f, fx, fy, fz);

        const double *p1x = gx - dk;
        const double *p1y = gy - dk;
        const double *p1z = gz - dk;
        const double *p2x = gx + dk;
        const double *p2y = gy + dk;
        const double *p2z = gz + dk;
        for (j = 0; j <= lj; j++)
        for (l = 0; l <= ll; l++) {
                ptr = dj * j + dl * l;
                //f(...,0,...) = -2*ak*g(...,1,...)
                for (i = 0; i <= li; i++) {
                        for (n = ptr; n < ptr+nroots; n++) {
                                fx[n] = ak2 * p2x[n];
                                fy[n] = ak2 * p2y[n];
                                fz[n] = ak2 * p2z[n];
                        }
                        ptr += di;
                }
                //f(...,k,...) = k*g(...,k-1,...)-2*ak*g(...,k+1,...)
                for (k = 1; k <= lk; k++) {
                        ptr = dj * j + dl * l + dk * k;
                        for (i = 0; i <= li; i++) {
                                for (n = ptr; n < ptr+nroots; n++) {
                                        fx[n] = k*p1x[n] + ak2*p2x[n];
                                        fy[n] = k*p1y[n] + ak2*p2y[n];
                                        fz[n] = k*p1z[n] + ak2*p2z[n];
                                }
                                ptr += di;
                        }
                }
        }
}


/*
 * ( ij | k \nabla l )
 */
void CINTnabla1l_2e(double *f, const double *g,
                    const FINT li, const FINT lj, const FINT lk, const FINT ll,
                    const CINTEnvVars *envs)
{
        FINT i, j, k, l, n, ptr;
        const FINT di = envs->g_stride_i;
        const FINT dk = envs->g_stride_k;
        const FINT dl = envs->g_stride_l;
        const FINT dj = envs->g_stride_j;
        const FINT nroots = envs->nrys_roots;
        const double al2 = -2 * envs->al[0];
        DEF_GXYZ(const double, g, gx, gy, gz);
        DEF_GXYZ(double, f, fx, fy, fz);

        const double *p1x = gx - dl;
        const double *p1y = gy - dl;
        const double *p1z = gz - dl;
        const double *p2x = gx + dl;
        const double *p2y = gy + dl;
        const double *p2z = gz + dl;
        for (j = 0; j <= lj; j++) {
                //f(...,0,...) = -2*al*g(...,1,...)
                for (k = 0; k <= lk; k++) {
                        ptr = dj * j + dk * k;
                        for (i = 0; i <= li; i++) {
                                for (n = ptr; n < ptr+nroots; n++) {
                                        fx[n] = al2 * p2x[n];
                                        fy[n] = al2 * p2y[n];
                                        fz[n] = al2 * p2z[n];
                                }
                                ptr += di;
                        }
                }
                //f(...,l,...) = l*g(...,l-1,...)-2*al*g(...,l+1,...)
                for (l = 1; l <= ll; l++) {
                        for (k = 0; k <= lk; k++) {
                                ptr = dj * j + dl * l + dk * k;
                                for (i = 0; i <= li; i++, ptr += di) {
                                for (n = ptr; n < ptr+nroots; n++) {
                                        fx[n] = l*p1x[n] + al2*p2x[n];
                                        fy[n] = l*p1y[n] + al2*p2y[n];
                                        fz[n] = l*p1z[n] + al2*p2z[n];
                                } }
                        }
                }
        }
}

/*
 * ( x^1 i j | kl )
 * ri is the shift from the center R_O to the center of |i>
 * r - R_O = (r-R_i) + ri, ri = R_i - R_O
 */
void CINTx1i_2e(double *f, const double *g, const double *ri,
                const FINT li, const FINT lj, const FINT lk, const FINT ll,
                const CINTEnvVars *envs)
{
        FINT i, j, k, l, n, ptr;
        const FINT di = envs->g_stride_i;
        const FINT dk = envs->g_stride_k;
        const FINT dl = envs->g_stride_l;
        const FINT dj = envs->g_stride_j;
        const FINT nroots = envs->nrys_roots;
        DEF_GXYZ(const double, g, gx, gy, gz);
        DEF_GXYZ(double, f, fx, fy, fz);

        const double *p1x = gx + di;
        const double *p1y = gy + di;
        const double *p1z = gz + di;
        for (j = 0; j <= lj; j++)
        for (l = 0; l <= ll; l++) {
        for (k = 0; k <= lk; k++) {
                //f(...,0:li,...) = g(...,1:li+1,...) + ri(1)*g(...,0:li,...)
                ptr = dj * j + dl * l + dk * k;
                for (i = 0; i <= li; i++) {
                        for (n = ptr; n < ptr+nroots; n++) {
                                fx[n] = p1x[n] + ri[0] * gx[n];
                                fy[n] = p1y[n] + ri[1] * gy[n];
                                fz[n] = p1z[n] + ri[2] * gz[n];
                        }
                        ptr += di;
                }
        } }
}


/*
 * ( i x^1 j | kl )
 */
void CINTx1j_2e(double *f, const double *g, const double *rj,
                const FINT li, const FINT lj, const FINT lk, const FINT ll,
                const CINTEnvVars *envs)
{
        FINT i, j, k, l, n, ptr;
        const FINT di = envs->g_stride_i;
        const FINT dk = envs->g_stride_k;
        const FINT dl = envs->g_stride_l;
        const FINT dj = envs->g_stride_j;
        const FINT nroots = envs->nrys_roots;
        DEF_GXYZ(const double, g, gx, gy, gz);
        DEF_GXYZ(double, f, fx, fy, fz);

        const double *p1x = gx + dj;
        const double *p1y = gy + dj;
        const double *p1z = gz + dj;
        for (j = 0; j <= lj; j++)
        for (l = 0; l <= ll; l++) {
        for (k = 0; k <= lk; k++) {
                // f(...,0:lj,...) = g(...,1:lj+1,...) + rj(1)*g(...,0:lj,...)
                ptr = dj * j + dl * l + dk * k;
                for (i = 0; i <= li; i++) {
                        for (n = ptr; n < ptr+nroots; n++) {
                                fx[n] = p1x[n] + rj[0] * gx[n];
                                fy[n] = p1y[n] + rj[1] * gy[n];
                                fz[n] = p1z[n] + rj[2] * gz[n];
                        }
                        ptr += di;
                }
        } }
}


/*
 * ( ij | x^1 k l )
 */
void CINTx1k_2e(double *f, const double *g, const double *rk,
                const FINT li, const FINT lj, const FINT lk, const FINT ll,
                const CINTEnvVars *envs)
{
        FINT i, j, k, l, n, ptr;
        const FINT di = envs->g_stride_i;
        const FINT dk = envs->g_stride_k;
        const FINT dl = envs->g_stride_l;
        const FINT dj = envs->g_stride_j;
        const FINT nroots = envs->nrys_roots;
        DEF_GXYZ(const double, g, gx, gy, gz);
        DEF_GXYZ(double, f, fx, fy, fz);

        const double *p1x = gx + dk;
        const double *p1y = gy + dk;
        const double *p1z = gz + dk;
        for (j = 0; j <= lj; j++)
        for (l = 0; l <= ll; l++) {
        for (k = 0; k <= lk; k++) {
                // f(...,0:lk,...) = g(...,1:lk+1,...) + rk(1)*g(...,0:lk,...)
                ptr = dj * j + dl * l + dk * k;
                for (i = 0; i <= li; i++) {
                        for (n = ptr; n < ptr+nroots; n++) {
                                fx[n] = p1x[n] + rk[0] * gx[n];
                                fy[n] = p1y[n] + rk[1] * gy[n];
                                fz[n] = p1z[n] + rk[2] * gz[n];
                        }
                        ptr += di;
                }
        } }
}


/*
 * ( i j | x^1 kl )
 */
void CINTx1l_2e(double *f, const double *g, const double *rl,
                const FINT li, const FINT lj, const FINT lk, const FINT ll,
                const CINTEnvVars *envs)
{
        FINT i, j, k, l, n, ptr;
        const FINT di = envs->g_stride_i;
        const FINT dk = envs->g_stride_k;
        const FINT dl = envs->g_stride_l;
        const FINT dj = envs->g_stride_j;
        const FINT nroots = envs->nrys_roots;
        DEF_GXYZ(const double, g, gx, gy, gz);
        DEF_GXYZ(double, f, fx, fy, fz);

        const double *p1x = gx + dl;
        const double *p1y = gy + dl;
        const double *p1z = gz + dl;
        for (j = 0; j <= lj; j++)
        for (l = 0; l <= ll; l++) {
        for (k = 0; k <= lk; k++) {
                // f(...,0:ll,...) = g(...,1:ll+1,...) + rl(1)*g(...,0:ll,...)
                ptr = dj * j + dl * l + dk * k;
                for (i = 0; i <= li; i++) {
                        for (n = ptr; n < ptr+nroots; n++) {
                                fx[n] = p1x[n] + rl[0] * gx[n];
                                fy[n] = p1y[n] + rl[1] * gy[n];
                                fz[n] = p1z[n] + rl[2] * gz[n];
                        }
                        ptr += di;
                }
        } }
}

