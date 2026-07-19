// rrtmg_lw_gpu_tables.h - Shared declarations for the RRTMG-LW GPU implmentation
#ifndef RRTMG_LW_GPU_TABLES_H
#define RRTMG_LW_GPU_TABLES_H

#define NBNDLW 16
#define NGPTLW 140
#define NTBL 10000
#define MAXLAY 128

// Per-band g-point counts and cumulative offsets
static const int LW_NG[16]     = {10,12,16,14,16,8,12,8,12,6,8,8,4,2,2,2};
static const int LW_NG_OFF[17] = {0,10,22,38,52,68,76,88,96,108,114,122,130,134,136,138,140};

// nspa/nspb for each band (control table indexing)
static const int LW_NSPA[16] = {1,1,9,9,9,1,9,1,9,1,1,9,9,1,9,9};
static const int LW_NSPB[16] = {1,1,5,5,5,0,1,1,1,1,1,0,0,1,0,0};

// ---- static coefficient tables (on device) --------------------------
struct LwTables {
    float heatfac, fluxfac, oneminus, abscld1, absliq0;

    // Global tables
    const int   *ngb;        // [140]  — band number per g-point
    const float *delwave;    // [16]   — spectral width per band
    const float *preflog;    // [59]   — reference pressure (log)
    const float *tref;       // [59]   — reference temperature
    const float *chi_mls;    // [7×59] — reference gas mixing ratios
    const float *totplnk;    // [181×16] — integrated Planck function
    const float *tau_tbl;    // [10001] — tau lookup table
    const float *exp_tbl;    // [10001] — exp(-tau) lookup table
    const float *tfn_tbl;    // [10001] — tau transition function table

    // Cloud optics tables
    const float *absice0;    // [2]
    const float *absice1;    // [2×5]
    const float *absice2;    // [43×16]
    const float *absice3;    // [46×16]
    const float *absliq1;    // [58×16]

    // Per-band k-tables
    const float *fracrefa[16], *fracrefb[16];
    const float *absa[16], *absb[16];
    const float *selfref[16], *forref[16];

    // Band-specific minor-gas tables
    const float *ka_mn2_1,   *kb_mn2_1;
    const float *ka_mn2o_3,  *kb_mn2o_3;
    const float *ka_mo3_5,   *ccl4_5;
    const float *ka_mco2_6,  *cfc11adj_6, *cfc12_6;
    const float *ka_mco2_7,  *kb_mco2_7;
    const float *ka_mco2_8,  *ka_mn2o_8, *ka_mo3_8, *kb_mco2_8, *kb_mn2o_8,
                *cfc12_8,    *cfc22adj_8;
    const float *ka_mn2o_9,  *kb_mn2o_9;
    const float *ka_mo2_11,  *kb_mo2_11;
    const float *ka_mco2_13, *ka_mco_13, *kb_mo3_13;
    const float *ka_mn2_15;
};

// Per-batch dimensions and option flags
struct LwDims {
    int nbatch, nlayers;
    int icld, inflglw, iceflglw, liqflglw;
};

// Device-side pointers to the packed column inputs
struct LwInputs {
    const float *play, *plev, *tlay, *tlev, *tsfc;
    const float *h2ovmr, *o3vmr, *co2vmr, *ch4vmr, *n2ovmr, *o2vmr;
    const float *cfc11vmr, *cfc12vmr, *cfc22vmr, *ccl4vmr;
    const float *emis;
    const float *cldfmcl, *ciwpmcl, *clwpmcl, *cswpmcl;
    const float *reicmcl, *relqmcl, *resnmcl;
    const float *tauaer;
    float *taucmcl;
};

// Per-batch scratch buffers (device memory)
struct LwScratch {    
    // ---- inatm outputs ----
    float *coldry, *wbrodl;     // (nb,L)
    float *wkl;                  // (nb,L,7)
    float *wx;                   // (nb,L,4)
    float *pwvcm;                // (nb)

    // ---- setcoef outputs ----
    int   *laytrop;              // (nb)
    int   *jp, *jt, *jt1;        // (nb,L)
    float *planklay;             // (nb,L,16)
    float *planklev;             // (nb,L+1,16)
    float *plankbnd;             // (nb,16)
    float *colh2o, *colco2, *colo3, *coln2o, *colco, *colch4, *colo2, *colbrd;
    float *fac00, *fac01, *fac10, *fac11;
    float *rat_h2oco2, *rat_h2oco2_1, *rat_h2oo3, *rat_h2oo3_1;
    float *rat_h2on2o, *rat_h2on2o_1, *rat_h2och4, *rat_h2och4_1;
    float *rat_n2oco2, *rat_n2oco2_1, *rat_o3co2, *rat_o3co2_1;
    float *selffac, *selffrac, *forfac, *forfrac;
    int   *indself, *indfor, *indminor;
    float *minorfrac, *scaleminor, *scaleminorn2;

    // ---- taumol outputs ----
    float *taug, *fracs;         // (140,nb,L)

    // ---- rtrnmc ----
    int   *icldlyr;              // (nb,L)
    float *urad_g, *drad_g;      // (140,nb,L+1)
    float *clrurad_g, *clrdrad_g;// (140,nb,L+1)
};

#endif // RRTMG_LW_GPU_TABLES_H