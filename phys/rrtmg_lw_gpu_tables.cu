// GPU table upload for RRTMG-LW.
// Called once from Fortran (wrf_gpu_lw_init_tables) after rrtmg_lw_ini has
// filled the rrlw_* module arrays.  Stores every pointer in a global
// LwTables struct that the pipeline kernels read.
#include "../inc/gpu_helpers.h"
#include "rrtmg_lw_gpu_tables.h"

static LwTables g_tab;
static bool     g_loaded = false;

const LwTables &lw_tables() {
    if (!g_loaded) { fprintf(stderr, "lw_tables: not loaded\n"); exit(1); }
    return g_tab;
}

#define UPLOAD(dst, src, n) do {                                         \
        size_t _n = (n);                                                  \
        float *_d;                                                        \
        CUDA_CHECK(cudaMalloc((void**)&_d, _n * sizeof(float)));          \
        CUDA_CHECK(cudaMemcpy(_d, (src), _n * sizeof(float),              \
                              cudaMemcpyHostToDevice));                   \
        (dst) = _d;                                                       \
    } while(0)

#define UPLOAD_INT(dst, src, n) do {                                      \
        size_t _n = (n);                                                  \
        int *_d;                                                          \
        CUDA_CHECK(cudaMalloc((void**)&_d, _n * sizeof(int)));            \
        CUDA_CHECK(cudaMemcpy(_d, (src), _n * sizeof(int),                \
                              cudaMemcpyHostToDevice));                   \
        (dst) = _d;                                                       \
    } while(0)

// ---- non-band-specific tables ------------------------------------------
extern "C" void wrf_gpu_lw_upload_globals(
    const float *heatfac, const float *fluxfac, const float *oneminus,
    const int   *ngb,
    const float *delwave,
    const float *preflog, const float *tref, const float *chi_mls,
    const float *totplnk,
    const float *tau_tbl, const float *exp_tbl, const float *tfn_tbl)
{
    g_tab.heatfac  = *heatfac;
    g_tab.fluxfac  = *fluxfac;
    g_tab.oneminus = *oneminus;

    UPLOAD_INT(g_tab.ngb, ngb, NGPTLW);
    UPLOAD(g_tab.delwave, delwave, NBNDLW);
    UPLOAD(g_tab.preflog, preflog, 59);
    UPLOAD(g_tab.tref,    tref, 59);
    UPLOAD(g_tab.chi_mls, chi_mls, 7 * 59);
    UPLOAD(g_tab.totplnk, totplnk, 181 * NBNDLW);
    UPLOAD(g_tab.tau_tbl, tau_tbl, NTBL + 1);
    UPLOAD(g_tab.exp_tbl, exp_tbl, NTBL + 1);
    UPLOAD(g_tab.tfn_tbl, tfn_tbl, NTBL + 1);
}

extern "C" void wrf_gpu_lw_upload_cloud(
    const float *abscld1,
    const float *absice0, const float *absice1,
    const float *absice2, const float *absice3,
    const float *absliq0, const float *absliq1)
{
    g_tab.abscld1 = *abscld1;
    UPLOAD(g_tab.absice0, absice0, 2);
    UPLOAD(g_tab.absice1, absice1, 2 * 5);
    UPLOAD(g_tab.absice2, absice2, 43 * NBNDLW);
    UPLOAD(g_tab.absice3, absice3, 46 * NBNDLW);
    g_tab.absliq0 = *absliq0;
    UPLOAD(g_tab.absliq1, absliq1, 58 * NBNDLW);
}

// ---- Per-band tables ----------------------------------------------------
// Called once per band from the Fortran init loop.
// Sets n_fb=0 for bands with no fracrefb, n_ab=0 for bands with no absb.
extern "C" void wrf_gpu_lw_upload_band(
    int band,
    const float *fracrefa, int n_fa,
    const float *fracrefb, int n_fb,
    const float *absa,     int n_aa,
    const float *absb,     int n_ab,
    const float *selfref,  int n_sr,
    const float *forref,   int n_fr)
{
    UPLOAD(g_tab.fracrefa[band], fracrefa, n_fa);
    if (n_fb) UPLOAD(g_tab.fracrefb[band], fracrefb, n_fb);
    else      g_tab.fracrefb[band] = nullptr;
    UPLOAD(g_tab.absa[band], absa, n_aa);
    if (n_ab) UPLOAD(g_tab.absb[band], absb, n_ab);
    else      g_tab.absb[band] = nullptr;
    UPLOAD(g_tab.selfref[band], selfref, n_sr);
    UPLOAD(g_tab.forref[band],  forref,  n_fr);
}

// ---- Band-specific minor-gas and cross-section tables ------------------

extern "C" void wrf_gpu_lw_upload_minor_b1(
    const float *ka_mn2_1, const float *kb_mn2_1, int ng) {
    UPLOAD(g_tab.ka_mn2_1,  ka_mn2_1, 19 * ng);
    UPLOAD(g_tab.kb_mn2_1,  kb_mn2_1, 19 * ng);
}

extern "C" void wrf_gpu_lw_upload_minor_b3(
    const float *ka, const float *kb, int ng) {
    UPLOAD(g_tab.ka_mn2o_3, ka, 9 * 19 * ng);
    UPLOAD(g_tab.kb_mn2o_3, kb, 5 * 19 * ng);
}

extern "C" void wrf_gpu_lw_upload_minor_b5(
    const float *ka_mo3, const float *ccl4, int ng) {
    UPLOAD(g_tab.ka_mo3_5, ka_mo3, 9 * 19 * ng);
    UPLOAD(g_tab.ccl4_5,   ccl4,   ng);
}

extern "C" void wrf_gpu_lw_upload_minor_b6(
    const float *ka_mco2, const float *cfc11, const float *cfc12, int ng) {
    UPLOAD(g_tab.ka_mco2_6,  ka_mco2, 19 * ng);
    UPLOAD(g_tab.cfc11adj_6, cfc11,   ng);
    UPLOAD(g_tab.cfc12_6,    cfc12,   ng);
}

extern "C" void wrf_gpu_lw_upload_minor_b7(
    const float *ka, const float *kb, int ng) {
    UPLOAD(g_tab.ka_mco2_7, ka, 9 * 19 * ng);
    UPLOAD(g_tab.kb_mco2_7, kb, 19 * ng);
}

extern "C" void wrf_gpu_lw_upload_minor_b8(
    const float *ka_co2, const float *ka_n2o, const float *ka_o3,
    const float *kb_co2, const float *kb_n2o,
    const float *cfc12, const float *cfc22, int ng) {
    UPLOAD(g_tab.ka_mco2_8, ka_co2, 19 * ng);
    UPLOAD(g_tab.ka_mn2o_8, ka_n2o, 19 * ng);
    UPLOAD(g_tab.ka_mo3_8,  ka_o3,  19 * ng);
    UPLOAD(g_tab.kb_mco2_8, kb_co2, 19 * ng);
    UPLOAD(g_tab.kb_mn2o_8, kb_n2o, 19 * ng);
    UPLOAD(g_tab.cfc12_8,   cfc12,  ng);
    UPLOAD(g_tab.cfc22adj_8,cfc22,  ng);
}

extern "C" void wrf_gpu_lw_upload_minor_b9(
    const float *ka, const float *kb, int ng) {
    UPLOAD(g_tab.ka_mn2o_9, ka, 9 * 19 * ng);
    UPLOAD(g_tab.kb_mn2o_9, kb, 19 * ng);
}

extern "C" void wrf_gpu_lw_upload_minor_b11(
    const float *ka, const float *kb, int ng) {
    UPLOAD(g_tab.ka_mo2_11, ka, 19 * ng);
    UPLOAD(g_tab.kb_mo2_11, kb, 19 * ng);
}

extern "C" void wrf_gpu_lw_upload_minor_b13(
    const float *ka_co2, const float *ka_co, const float *kb_o3, int ng) {
    UPLOAD(g_tab.ka_mco2_13, ka_co2, 9 * 19 * ng);
    UPLOAD(g_tab.ka_mco_13,  ka_co,  9 * 19 * ng);
    UPLOAD(g_tab.kb_mo3_13,  kb_o3,  19 * ng);
}

extern "C" void wrf_gpu_lw_upload_minor_b15(
    const float *ka, int ng) {
    UPLOAD(g_tab.ka_mn2_15, ka, 9 * 19 * ng);
}

// Mark tables as ready — called once after all uploads complete.
extern "C" void wrf_gpu_lw_tables_loaded() {
    g_loaded = true;
    printf("lw_tables: uploaded (heatfac=%g fluxfac=%g oneminus=%.7f)\n",
           g_tab.heatfac, g_tab.fluxfac, g_tab.oneminus);
}
