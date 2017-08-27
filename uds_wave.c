#include "udf.h"
#include "sg.h"

enum
{
    RHO,
    ZETA,
    N_REQURIED_UDS
};

DEFINE_UDS_UNSTEADY(uds_time, c, t, i, apu, su)
{
    real physical_dt, vol, rho0 = 998.2, phi_old;
    physical_dt = RP_Get_Real("physical-time-step");
    vol = C_VOLUME(c,t);

    if(i == RHO) {
        *apu = -vol / physical_dt;
        *su = vol*C_STORAGE_R(c,t,SV_UDSI_M1(i))/physical_dt;
    } else if(i == ZETA){
        *apu = -rho0*vol / physical_dt;
        *su = rho0*vol*C_STORAGE_R(c,t,SV_UDSI_M1(i))/physical_dt;
    }
}

DEFINE_UDS_FLUX(uds_flux, f, t, i)
{
    real k = 1483.0*1483.0, source, Gsource;
    cell_t c0, c1;
    Thread *t1, *t0;
    real A[ND_ND], dRHO[ND_ND], dr0[ND_ND], dr1[ND_ND], es[ND_ND], A_by_es, ds;
    c0 = F_C0(f, t);
    t0 = THREAD_T0(t);
    if (i != ZETA) return 0.;

    if (BOUNDARY_FACE_THREAD_P(t))
    {
	if(NULL != THREAD_STORAGE(t, SV_UDS_I(RHO)) {
	    BOUNDARY_FACE_GEOMETRY(f, t, A, ds, es, A_by_es, dr0);
	    source = k * (F_UDSI(f, t, RHO) - C_UDSI(c0, t0, RHO)) * A_by_es / ds;
            if(NULL != T_STORAGE_R_NV(t0, SV_UDSI_G(RHO))) {
                BOUNDARY_SECONDARY_GRADIENT_SOURCE(Gsource, SV_UDSI_G(RHO), dRHO, es, A_by_es, k);
                source += Gsource;
	    }
        } else {
            source = 0.;
        }
    } else {
        c1 = F_C1(f, t);
        t1 = THREAD_T1(t);
        INTERIOR_FACE_GEOMETRY(f, t, A, ds, es, A_by_es, dr0, dr1);
        source = k * (C_UDSI(c1, t1, RHO) - C_UDSI(c0, t0, RHO)) * A_by_es / ds;
	if(NULL != T_STORAGE_R_NV(t0, SV_UDSI_G(RHO))) {
	    SECONDARY_GRADIENT_SOURCE(Gsource, SV_UDSI_G(RHO), dRHO, es, A_by_es, k);
	    source += Gsource;
	}
    }
    return source;
}

DEFINE_SOURCE(uds0, c, t, dS, eqn)
{
    real rho0 = 998.2;
    dS[eqn] = 0.;
    return - rho0 * C_UDSI(c, t, ZETA);
}

// boundary condition rho, \nabla \cdot V
DEFINE_PROFILE(myprof0, t, i) 
{
    face_t f;
    real freq = 2.0e+4, rho = 998.2, c = 1483.0;
    real omega = 2.0 * M_PI * freq;
    real A = 40e-6;
    real time = CURRENT_TIME;
    begin_f_loop(f, t) 
    {
	F_PROFILE(f, t, i) = A*omega*rho/c * (cos(omega*time) - 1.0);    
    }
    end_f_loop(f, t)
}
	   
DEFINE_PROFILE(myprof1, t, i)
{
    face_t f;
    real freq = 2.0e+4, rho = 998.2, c = 1483.0;
    real omega = 2.0 * M_PI * freq;
    real A = 40e-6;
    begin_f_loop(f, t)
    {
	F_PROFILE(f, t) = M_PI*pow(5e-3, 2)*A*omega * sin(omega*time);
    }
    end_f_loop(f, t)
}
