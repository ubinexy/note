
#include "udf.h"
#include "sg.h"

enum {
    RHO,
    ZETA,
    OMEGA,
    PHI,
    PSI,
    N_REQUIRED_UDS
};

DEFINE_UDS_UNSTEADY(time_rho, c, t, i, apu, su)
{
    real physical_dt, vol, rho0 = 998.2;
    vol = C_VOLUME(c, t);
    physical_dt = RP_Get_Real("physical-time-step");
    if(i == RHO) {
        *apu = - rho0 * vol / physical_dt;
        *su = rho0 * vol * C_STORAGE_R(c, t, SV_UDSI_M1(i))/ physical_dt;
    }
}

DEFINE_UDS_UNSTEADY(time_zeta, c, t, i, apu, su)
{
    real physical_dt, vol, rho0 = 998.2;
    vol = C_VOLUME(c, t);
    physical_dt = RP_Get_Real("physcial-time-step");

    if( i == ZETA ) {
        *apu = -rho0 * vol / physical_dt;
        *su = rho0 * vol * C_STORAGE_R(c, t, SV_UDSI_M1(i))/ physical_dt;
    }
}

DEFINE_UDS_UNSTEADY(time_omega, c, t, i, apu, su)
{
    real physical_dt, vol, rho0 = 998.2;
    vol = C_VOLUME(c, t);
    physical_dt = RP_Get_Real("physcial-time-step");
    
    if( i == OMEGA ) {
         *apu = - rho0 * vol / physical_dt; 
         *su = rho0 * vol *C_STORAGE_R(c, t, SV_UDSI_M1(i))/ physical_dt;
    }
}

DEFINE_UDS_UNSTEADY(time_phi, c, t, i, apu, su)
{
    if( i == PHI ) {
        *apu = 0.0;
        *su = 0.0;
    }
}

DEFINE_USD_UNSTEADY(time_psi, c, t, i, apu, su)
{
    if( i == PSI ) {
        *apu = 0.0;
        *su = 0.0;
    }  
}

//*********************************************************************************//

DEFINE_FLUX(flux_zeta, f, t, i)
{
    real k = 1483.0*1483.0;
    cell_t c0, c1;
    Thread *t1, *t0;
    real A[ND_ND], dRHO[ND_ND], dr0[ND_ND], dr1[ND_ND], es[ND_ND], A_by_es, ds, source;
    c0 = F_C0(f, t);
    t0 = THREAD_T0(t);

    if (i != ZETA) return 0.;

    if (BOUNDARY_FACE_THREAD_P(t)) {
        // \nabla \cdot (c^2 \nabla p)
        if(NULL != THREAD_STORAGE(t, SV_UDS_I(RHO)) ) {
            BOUNDARY_FACE_GEOMETRY(f, t, A, ds, es, A_by_es, dr0);
            source = k * (F_UDSI(f, t, RHO) - C_UDSI(c0, t0, RHO)) * A_by_es / ds;
            if( NULL != T_STORAGE_R_NV(t0, SV_UDSI_G(RHO)) ) {
                BOUNDARY_SECONDARY_GRADIENT_SOURCE(GSource, SV_UDSI_G(RHO), dRHO, es, A_by_es, k);
                source += GSource;
            }
        } else {
            source = 0.;
        }
        // \nabla \cdot (bar{v} zeta)   
        if(NULL != THREAD_STORAGE(t, SV_UDS_I(ZETA)) {
            zeta =F_UDSI(f, t, ZETA);
        } else if(NULL != T_STORAGE_R_NV(t0, SV_UDSI_G(ZETA))) {
            BOUNDARY_FACE_GEOMETRY(f, t, A, ds, es, A_by_es, dr0);
            zeta =C_UDSI(c0, t0, ZETA) + NV_DOT(C_UDSI_G(c0, t0, ZETA), dr0);
        } else {
            zeta =C_UDSI(c0, t0, ZETA);
        }
        NV_DS(V, =, F_U(f,t), F_V(c,t), 0, *, zeta);
        return rho0*NV_DOT(V, A) + source;
    } else {
        c1 = F_C1(f, t);
        t1 = THREAD_T1(t);
        // \nabla \cdot (c^2 \nabla p)
        INTERIOR_FACE_GEOMETRY(f, t, A, ds, es, A_by_es, dr0, dr1);
        source = k * (C_UDSI(c1, t1, rho) - C_UDSI(c0, t0, rho)) * A_by_es / ds;
        if( NULL != T_STORAGE_R_NV(t0, SV_UDSI_G(RHO)) ) {
            SECONDARY_GRADIENT_SOURCE(GSource, SV_UDSI_G(RHO), dRHO, es, A_by_es, k);
            source += GSource;
        }
        // \nabla \cdot (bar{v} zeta)  
        if(NULL != T_STORAGE_R_NV(t0, SV_UDSI_G(ZETA)) ) {
            //INTERIOR_FACE_GEOMETRY(f, t, A, ds, es, A_by_es, dr0, dr1);
            ZETA = C_UDSI(c0, t0, ZETA) + NV_DOT(C_UDSI_G(c0, t0, ZETA), dr0);
        } else {
            ZETA = (C_UDSI(c0, t0, ZETA) + C_UDSI(c1, t1, ZETA))/2.0;
        }
        NV_D(V, =, C_U(c0, t0), C_V(c0, t0), 0);
        NV_D(V,+=, C_U(c1, t1), C_V(c1, t1), 0);
        return rho0*NV_DOT(V, A)/2.0*zeta + source;
    }
}

DEFINE_FLUX(flux_omega, f, t, i)
{   
    cell_t c0, c1;
    Thread *t0, *t1;
    real A[ND_ND], dRHO[ND_ND], dr0[ND_ND], dr1[ND_ND], es[ND_ND], A_by_es, ds, source;

    c0 = F_C0(f, t);
    t0 = THREAD_T0(t);
    
    if (i != OMEGA) return 0.;

    if (BOUNDARY_FACE_THREAD_P(t)) {
        if(NULL != THREAD_STORAGE(t, SV_UDS_I(OMEGA)) ) {
            omega = F_UDSI(f, t, OMEGA);
        } else if ( NULL != T_STORAGE_R_NV(t0, SV_UDSI_G(OMEGA) ) {
            BOUNDARY_FACE_GEOMETRY(f, t, A, ds, es, A_by_es, dr0);
            omega = C_UDSI(c0, t0, OMEGA) + NV_DOT(C_UDSI_G(c0, t0, OMEGA), dr0);
        } else {
            omega = C_UDSI(c0, t0, OMEGA);
        }
        NV_D(V, =, F_U(f, t), F_V(f, t), 0, omega);
        NV_DS(VOMGA, =, C_UDMI(c, t, 1), C_UDMI(c, t, 2), 0, C_UDMI(c, t, 0));
        return rho0*NV_DOT(V, A) + rho0*NVDOT(VOMGA, A);
    } else {
        c1 = F_C1(f, t);
        t1 = THREAD_T1(t);
        if(NULL != T_STORAGE_R_NV(t0, SV_UDSI_G(OMEGA)) ) {
            INTERIOR_FACE_GEOMETRY(f, t, A, ds, es, A_by_es, dr0, dr1);
            omega = C_UDSI(c0, t0, OMEGA) + NV_DOT(C_UDSI_G(c0, t0, OMEGA), dr0);
        } else {
            omega = (C_UDSI(c0, t0, OMEGA) + C_UDSI(c1, t1, OMEGA))/2.0;
        }
        NV_D(V, =, C_U(c0, t0), C_V(c0, t0), 0);
        NV_D(V,+=, C_U(c1, t1), C_V(c1, t1), 0);
        NV_DS(VOMGA, =, F_UDMI(f, t, 1), C_UDMI(f, t, 2), 0, C_UDMI(f, t, 0));
        return rho0*NV_DOT(V, A)/2.0*omega + rho0*NVDOT(VOMGA, A);
    }
}


DEFINE_FLUX(flux_rho, f, t, i)
{
    cell_t c0, c1;
    Thread *t0, *t1;
    real A[ND_ND], dRHO[ND_ND], dr0[ND_ND], dr1[ND_ND], es[ND_ND], A_by_es, ds, source;

    c0 = F_C0(f, t);
    t0 = THREAD_T0(t);

    if (i != RHO) return 0.;

    if (BOUNDARY_FACE_THREAD_P(t)) {
        if(NULL != THREAD_STORAGE(t, SV_UDS_I(RHO)) ) {
            rho = F_UDSI(f, t);
        } else if ( NULL != T_STORAGE_R_NV(t0, SV_UDSI_G(RHO) ) {
            BOUNDARY_FACE_GEOMETRY(f, t, A, ds, es, A_by_es, dr0);
            rho = C_UDSI(c0, t0) + NV_DOT(C_UDSI_G(c0, t0, RHO), dr0);
        } else {
            rho = C_UDSI(c0, t0);
        }
        NV_D(V, =, F_U(f, t), F_V(f, t), 0, RHO);
        return NV_DOT(V, A);
    } else {
        c1 = F_C1(f, t);
        t1 = THREAD_T1(t);
        if(NULL != T_STORAGE_R_NV(t0, SV_UDSI_G(RHO))) {
            INTERIOR_FACE_GEOMETRY(f, t, A, ds, es, A_by_es, dr0, dr1);
            rho = C_UDSI(c0, t0, RHO) + NV_DOT(C_UDSI_G(c0, t0, RHO), dr0);
        } else {
            rho = (C_UDSI(c0, t0, RHO) + C_UDSI(c1, t1, RHO))/2.0;
        }
        NV_D(V, =, C_U(c0, t0), C_V(c0, t0), 0);
        NV_D(V,+=, C_U(c1, t1), C_V(c1, t1), 0);
        return NV_DOT(V, A)/2.0*rho;
    }
}

//**********************************************************************//
           
DEFINE_ADJUST(def_omega_bar, d)
{
    // define omega bar
    Thread *t, *t0;
    cell_t c, c0; 
    face_t f;
  
    thread_loop_c(t, d)
    {
        begin_c_loop(c, t)
        {
            C_UDMI(c, t, 0) = C_U_G(c, t)[1] - C_V_G(c, t)[0];
        }
        end_c_loop(c, t)
    }
  
    thread_loop_f(t, d)
    {
        c0 = F_C0(f, t);
        t0 = THREAD_T0(t);
        if(NULL != T_STORAGE_R_NV(t0, SV_U_G))
        {
            begin_f_loop(f, t)
            {
                F_UDMI(f, t, 0) = C_U_G(c, t)[1] - C_V_G(c, t)[0];  
            }
            end_f_loop(f, t)
        }
    }
    // define v^\prime
  
    thread_loop_c(t, d)
    {
      if(NULL != T_STORAGE_R_NV(t, SV_UDSI_G(PHI))) 
      {
        begin_c_loop(c, t)
        {
            C_CENTROID(x, c, t);
            C_UDMI(c, t, 1) = C_UDSI_G(c, t, PHI)[0] + C_UDSI_G(c, t, PSI)[1];
            C_UDMI(c, t, 2) = C_UDSI_G(c, t, PHI)[1] - C_UDSI_G(c, t, PSI)[0] - C_UDSI(c, t, PSI)/x[0];
        }
        end_c_loop(c, t)
      }
    }
    
    thread_loop_f(t, d)
    {
        c0 = F_C0(f, t);
        t0 = THREAD_T0(t);
        if(NULL != T_STORAGE_R_NV(t0, SV_UDSI_G(PHI))) 
        {
          begin_f_loop(f, t)
          {
              C_CENTROID(x, c0, t0);
              F_UDMI(f, t, 1) = C_UDSI_G(c0, t0, PHI)[0] + C_UDSI_G(c0, t0, PSI)[1];
              F_UDMI(f, t, 2) = C_UDSI_G(c0, t0, PHI)[1] - C_UDSI_G(c0, t0, PSI)[0] - C_UDSI(c0, t0, PSI)/x[0];
          }
          end_f_loop(f, t)
        }
    }
     
}
//**********************************************************************//
DEFINE_SOURCE(source_rho, c, t, i, dS, eqn)
{
    real rho0 = 998.2;
    return - rho0* C_UDSI(c, t, ZETA);
}

DEFINE_SOURCE(source_zeta, c, t, i, dS, eqn)
{
    //real rho0 = 998.2;
    // dS[eqn] = 0.0;
    
    return - 0.0;
}
           
DEFINE_SOURCE(source_omega, c, t, i, dS, eqn)
{
    real rho0 = 998.2;
    return rho0 * omega_bar * C_UDSI(c, t, ZETA);
}

DEFINE_SOURCE(source_phi, c, t, i, dS, eqn)
{
    return C_UDSI(c, t, ZETA);
}

DEFINE_SOURCE(source_psi, c, t, i, dS, eqn)
{
    return C_UDSI(c, t, OMEGA);
}
//*********************************************************************//
