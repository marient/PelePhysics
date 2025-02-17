
#include <math.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>

#if defined(BL_FORT_USE_UPPERCASE)
#define CKINDX CKINDX
#define CKINIT CKINIT
#define CKFINALIZE CKFINALIZE
#define CKXNUM CKXNUM
#define CKSYME CKSYME
#define CKSYMS CKSYMS
#define CKRP CKRP
#define CKPX CKPX
#define CKPY CKPY
#define CKPC CKPC
#define CKRHOX CKRHOX
#define CKRHOY CKRHOY
#define CKRHOC CKRHOC
#define CKWT CKWT
#define CKAWT CKAWT
#define CKMMWY CKMMWY
#define CKMMWX CKMMWX
#define CKMMWC CKMMWC
#define CKYTX CKYTX
#define CKYTCP CKYTCP
#define CKYTCR CKYTCR
#define CKXTY CKXTY
#define CKXTCP CKXTCP
#define CKXTCR CKXTCR
#define CKCTX CKCTX
#define CKCTY CKCTY
#define CKCPOR CKCPOR
#define CKHORT CKHORT
#define CKSOR CKSOR
#define CKCVML CKCVML
#define CKCPML CKCPML
#define CKUML CKUML
#define CKHML CKHML
#define CKGML CKGML
#define CKAML CKAML
#define CKSML CKSML
#define CKCVMS CKCVMS
#define CKCPMS CKCPMS
#define CKUMS CKUMS
#define CKHMS CKHMS
#define CKGMS CKGMS
#define CKAMS CKAMS
#define CKSMS CKSMS
#define CKCPBL CKCPBL
#define CKCPBS CKCPBS
#define CKCVBL CKCVBL
#define CKCVBS CKCVBS
#define CKHBML CKHBML
#define CKHBMS CKHBMS
#define CKUBML CKUBML
#define CKUBMS CKUBMS
#define CKSBML CKSBML
#define CKSBMS CKSBMS
#define CKGBML CKGBML
#define CKGBMS CKGBMS
#define CKABML CKABML
#define CKABMS CKABMS
#define CKWC CKWC
#define CKWYP CKWYP
#define CKWXP CKWXP
#define CKWYR CKWYR
#define CKWXR CKWXR
#define CKQC CKQC
#define CKKFKR CKKFKR
#define CKQYP CKQYP
#define CKQXP CKQXP
#define CKQYR CKQYR
#define CKQXR CKQXR
#define CKNU CKNU
#define CKNCF CKNCF
#define CKABE CKABE
#define CKEQC CKEQC
#define CKEQYP CKEQYP
#define CKEQXP CKEQXP
#define CKEQYR CKEQYR
#define CKEQXR CKEQXR
#define DWDOT DWDOT
#define DWDOT_PRECOND DWDOT_PRECOND
#define SPARSITY_INFO SPARSITY_INFO
#define SPARSITY_INFO_PRECOND SPARSITY_INFO_PRECOND
#define SPARSITY_PREPROC SPARSITY_PREPROC
#define SPARSITY_PREPROC_PRECOND SPARSITY_PREPROC_PRECOND
#define VCKHMS VCKHMS
#define VCKPY VCKPY
#define VCKWYR VCKWYR
#define VCKYTX VCKYTX
#define GET_T_GIVEN_EY GET_T_GIVEN_EY
#define GET_T_GIVEN_HY GET_T_GIVEN_HY
#define GET_REACTION_MAP GET_REACTION_MAP
#define GET_CRITPARAMS GET_CRITPARAMS
#elif defined(BL_FORT_USE_LOWERCASE)
#define CKINDX ckindx
#define CKINIT ckinit
#define CKFINALIZE ckfinalize
#define CKXNUM ckxnum
#define CKSYME cksyme
#define CKSYMS cksyms
#define CKRP ckrp
#define CKPX ckpx
#define CKPY ckpy
#define CKPC ckpc
#define CKRHOX ckrhox
#define CKRHOY ckrhoy
#define CKRHOC ckrhoc
#define CKWT ckwt
#define CKAWT ckawt
#define CKMMWY ckmmwy
#define CKMMWX ckmmwx
#define CKMMWC ckmmwc
#define CKYTX ckytx
#define CKYTCP ckytcp
#define CKYTCR ckytcr
#define CKXTY ckxty
#define CKXTCP ckxtcp
#define CKXTCR ckxtcr
#define CKCTX ckctx
#define CKCTY ckcty
#define CKCPOR ckcpor
#define CKHORT ckhort
#define CKSOR cksor
#define CKCVML ckcvml
#define CKCPML ckcpml
#define CKUML ckuml
#define CKHML ckhml
#define CKGML ckgml
#define CKAML ckaml
#define CKSML cksml
#define CKCVMS ckcvms
#define CKCPMS ckcpms
#define CKUMS ckums
#define CKHMS ckhms
#define CKGMS ckgms
#define CKAMS ckams
#define CKSMS cksms
#define CKCPBL ckcpbl
#define CKCPBS ckcpbs
#define CKCVBL ckcvbl
#define CKCVBS ckcvbs
#define CKHBML ckhbml
#define CKHBMS ckhbms
#define CKUBML ckubml
#define CKUBMS ckubms
#define CKSBML cksbml
#define CKSBMS cksbms
#define CKGBML ckgbml
#define CKGBMS ckgbms
#define CKABML ckabml
#define CKABMS ckabms
#define CKWC ckwc
#define CKWYP ckwyp
#define CKWXP ckwxp
#define CKWYR ckwyr
#define CKWXR ckwxr
#define CKQC ckqc
#define CKKFKR ckkfkr
#define CKQYP ckqyp
#define CKQXP ckqxp
#define CKQYR ckqyr
#define CKQXR ckqxr
#define CKNU cknu
#define CKNCF ckncf
#define CKABE ckabe
#define CKEQC ckeqc
#define CKEQYP ckeqyp
#define CKEQXP ckeqxp
#define CKEQYR ckeqyr
#define CKEQXR ckeqxr
#define DWDOT dwdot
#define DWDOT_PRECOND dwdot_precond
#define SPARSITY_INFO sparsity_info
#define SPARSITY_INFO_PRECOND sparsity_info_precond
#define SPARSITY_PREPROC sparsity_preproc
#define SPARSITY_PREPROC_PRECOND sparsity_preproc_precond
#define VCKHMS vckhms
#define VCKPY vckpy
#define VCKWYR vckwyr
#define VCKYTX vckytx
#define GET_T_GIVEN_EY get_t_given_ey
#define GET_T_GIVEN_HY get_t_given_hy
#define GET_REACTION_MAP get_reaction_map
#define GET_CRITPARAMS get_critparams
#elif defined(BL_FORT_USE_UNDERSCORE)
#define CKINDX ckindx_
#define CKINIT ckinit_
#define CKFINALIZE ckfinalize_
#define CKXNUM ckxnum_
#define CKSYME cksyme_
#define CKSYMS cksyms_
#define CKRP ckrp_
#define CKPX ckpx_
#define CKPY ckpy_
#define CKPC ckpc_
#define CKRHOX ckrhox_
#define CKRHOY ckrhoy_
#define CKRHOC ckrhoc_
#define CKWT ckwt_
#define CKAWT ckawt_
#define CKMMWY ckmmwy_
#define CKMMWX ckmmwx_
#define CKMMWC ckmmwc_
#define CKYTX ckytx_
#define CKYTCP ckytcp_
#define CKYTCR ckytcr_
#define CKXTY ckxty_
#define CKXTCP ckxtcp_
#define CKXTCR ckxtcr_
#define CKCTX ckctx_
#define CKCTY ckcty_
#define CKCPOR ckcpor_
#define CKHORT ckhort_
#define CKSOR cksor_
#define CKCVML ckcvml_
#define CKCPML ckcpml_
#define CKUML ckuml_
#define CKHML ckhml_
#define CKGML ckgml_
#define CKAML ckaml_
#define CKSML cksml_
#define CKCVMS ckcvms_
#define CKCPMS ckcpms_
#define CKUMS ckums_
#define CKHMS ckhms_
#define CKGMS ckgms_
#define CKAMS ckams_
#define CKSMS cksms_
#define CKCPBL ckcpbl_
#define CKCPBS ckcpbs_
#define CKCVBL ckcvbl_
#define CKCVBS ckcvbs_
#define CKHBML ckhbml_
#define CKHBMS ckhbms_
#define CKUBML ckubml_
#define CKUBMS ckubms_
#define CKSBML cksbml_
#define CKSBMS cksbms_
#define CKGBML ckgbml_
#define CKGBMS ckgbms_
#define CKABML ckabml_
#define CKABMS ckabms_
#define CKWC ckwc_
#define CKWYP ckwyp_
#define CKWXP ckwxp_
#define CKWYR ckwyr_
#define CKWXR ckwxr_
#define CKQC ckqc_
#define CKKFKR ckkfkr_
#define CKQYP ckqyp_
#define CKQXP ckqxp_
#define CKQYR ckqyr_
#define CKQXR ckqxr_
#define CKNU cknu_
#define CKNCF ckncf_
#define CKABE ckabe_
#define CKEQC ckeqc_
#define CKEQYP ckeqyp_
#define CKEQXP ckeqxp_
#define CKEQYR ckeqyr_
#define CKEQXR ckeqxr_
#define DWDOT dwdot_
#define DWDOT_PRECOND dwdot_precond_
#define SPARSITY_INFO sparsity_info_
#define SPARSITY_INFO_PRECOND sparsity_info_precond_ 
#define SPARSITY_PREPROC sparsity_preproc_
#define SPARSITY_PREPROC_PRECOND sparsity_preproc_precond_
#define VCKHMS vckhms_
#define VCKPY vckpy_
#define VCKWYR vckwyr_
#define VCKYTX vckytx_
#define GET_T_GIVEN_EY get_t_given_ey_
#define GET_T_GIVEN_HY get_t_given_hy_
#define GET_REACTION_MAP get_reaction_map_
#define GET_CRITPARAMS get_critparams_
#endif

/*function declarations */
#if defined(BL_FORT_USE_UPPERCASE)
#define egtransetEPS EGTRANSETEPS
#elif defined(BL_FORT_USE_LOWERCASE)
#define egtransetEPS egtranseteps
#elif defined(BL_FORT_USE_UNDERSCORE)
#define egtransetEPS egtranseteps_
#endif
#if defined(BL_FORT_USE_UPPERCASE)
#define egtransetSIG EGTRANSETSIG
#elif defined(BL_FORT_USE_LOWERCASE)
#define egtransetSIG egtransetsig
#elif defined(BL_FORT_USE_UNDERSCORE)
#define egtransetSIG egtransetsig_
#endif
extern "C"
{
void egtransetEPS(double *  EPS);
void egtransetSIG(double* SIG);
void atomicWeight(double *  awt);
void molecularWeight(double *  wt);
void gibbs(double *  species, double *  tc);
void helmholtz(double *  species, double *  tc);
void speciesInternalEnergy(double *  species, double *  tc);
void speciesEnthalpy(double *  species, double *  tc);
void speciesEntropy(double *  species, double *  tc);
void cp_R(double *  species, double *  tc);
void cv_R(double *  species, double *  tc);
void equilibriumConstants(double *  kc, double *  g_RT, double T);
void productionRate(double *  wdot, double *  sc, double T);
void comp_k_f(double *  tc, double invT, double *  k_f);
void comp_Kc(double *  tc, double invT, double *  Kc);
void comp_qfqr(double *  q_f, double *  q_r, double *  sc, double *  tc, double invT);
void progressRate(double *  qdot, double *  speciesConc, double T);
void progressRateFR(double *  q_f, double *  q_r, double *  speciesConc, double T);
void CKINIT();
void CKFINALIZE();
void CKINDX(int * mm, int * kk, int * ii, int * nfit );
void CKXNUM(char * line, int * nexp, int * lout, int * nval, double *  rval, int * kerr, int lenline);
void CKSNUM(char * line, int * nexp, int * lout, char * kray, int * nn, int * knum, int * nval, double *  rval, int * kerr, int lenline, int lenkray);
void CKSYME(int * kname, int * lenkname);
void CKSYMS(int * kname, int * lenkname);
void CKRP(double *  ru, double *  ruc, double *  pa);
void CKPX(double *  rho, double *  T, double *  x, double *  P);
void CKPY(double *  rho, double *  T, double *  y, double *  P);
void CKPC(double *  rho, double *  T, double *  c, double *  P);
void CKRHOX(double *  P, double *  T, double *  x, double *  rho);
void CKRHOY(double *  P, double *  T, double *  y, double *  rho);
void CKRHOC(double *  P, double *  T, double *  c, double *  rho);
void CKWT(double *  wt);
void CKAWT(double *  awt);
void CKMMWY(double *  y, double *  wtm);
void CKMMWX(double *  x, double *  wtm);
void CKMMWC(double *  c, double *  wtm);
void CKYTX(double *  y, double *  x);
void CKYTCP(double *  P, double *  T, double *  y, double *  c);
void CKYTCR(double *  rho, double *  T, double *  y, double *  c);
void CKXTY(double *  x, double *  y);
void CKXTCP(double *  P, double *  T, double *  x, double *  c);
void CKXTCR(double *  rho, double *  T, double *  x, double *  c);
void CKCTX(double *  c, double *  x);
void CKCTY(double *  c, double *  y);
void CKCPOR(double *  T, double *  cpor);
void CKHORT(double *  T, double *  hort);
void CKSOR(double *  T, double *  sor);
void CKCVML(double *  T, double *  cvml);
void CKCPML(double *  T, double *  cvml);
void CKUML(double *  T, double *  uml);
void CKHML(double *  T, double *  uml);
void CKGML(double *  T, double *  gml);
void CKAML(double *  T, double *  aml);
void CKSML(double *  T, double *  sml);
void CKCVMS(double *  T, double *  cvms);
void CKCPMS(double *  T, double *  cvms);
void CKUMS(double *  T, double *  ums);
void CKHMS(double *  T, double *  ums);
void CKGMS(double *  T, double *  gms);
void CKAMS(double *  T, double *  ams);
void CKSMS(double *  T, double *  sms);
void CKCPBL(double *  T, double *  x, double *  cpbl);
void CKCPBS(double *  T, double *  y, double *  cpbs);
void CKCVBL(double *  T, double *  x, double *  cpbl);
void CKCVBS(double *  T, double *  y, double *  cpbs);
void CKHBML(double *  T, double *  x, double *  hbml);
void CKHBMS(double *  T, double *  y, double *  hbms);
void CKUBML(double *  T, double *  x, double *  ubml);
void CKUBMS(double *  T, double *  y, double *  ubms);
void CKSBML(double *  P, double *  T, double *  x, double *  sbml);
void CKSBMS(double *  P, double *  T, double *  y, double *  sbms);
void CKGBML(double *  P, double *  T, double *  x, double *  gbml);
void CKGBMS(double *  P, double *  T, double *  y, double *  gbms);
void CKABML(double *  P, double *  T, double *  x, double *  abml);
void CKABMS(double *  P, double *  T, double *  y, double *  abms);
void CKWC(double *  T, double *  C, double *  wdot);
void CKWYP(double *  P, double *  T, double *  y, double *  wdot);
void CKWXP(double *  P, double *  T, double *  x, double *  wdot);
void CKWYR(double *  rho, double *  T, double *  y, double *  wdot);
void CKWXR(double *  rho, double *  T, double *  x, double *  wdot);
void CKQC(double *  T, double *  C, double *  qdot);
void CKKFKR(double *  P, double *  T, double *  x, double *  q_f, double *  q_r);
void CKQYP(double *  P, double *  T, double *  y, double *  qdot);
void CKQXP(double *  P, double *  T, double *  x, double *  qdot);
void CKQYR(double *  rho, double *  T, double *  y, double *  qdot);
void CKQXR(double *  rho, double *  T, double *  x, double *  qdot);
void CKNU(int * kdim, int * nuki);
void CKNCF(int * mdim, int * ncf);
void CKABE(double *  a, double *  b, double *  e );
void CKEQC(double *  T, double *  C , double *  eqcon );
void CKEQYP(double *  P, double *  T, double *  y, double *  eqcon);
void CKEQXP(double *  P, double *  T, double *  x, double *  eqcon);
void CKEQYR(double *  rho, double *  T, double *  y, double *  eqcon);
void CKEQXR(double *  rho, double *  T, double *  x, double *  eqcon);
void DWDOT(double *  J, double *  sc, double *  T, int * consP);
void DWDOT_PRECOND(double *  J, double *  sc, double *  Tp, int * HP);
void SPARSITY_INFO(int * nJdata, int * consP, int NCELLS);
void SPARSITY_INFO_PRECOND(int * nJdata, int * consP);
void SPARSITY_PREPROC(int * rowVals, int * colPtrs, int * consP, int NCELLS);
void SPARSITY_PREPROC_PRECOND(int * rowVals, int * colPtrs, int * consP);
void aJacobian(double *  J, double *  sc, double T, int consP);
void aJacobian_precond(double *  J, double *  sc, double T, int HP);
void dcvpRdT(double *  species, double *  tc);
void GET_T_GIVEN_EY(double *  e, double *  y, double *  t, int *ierr);
void GET_T_GIVEN_HY(double *  h, double *  y, double *  t, int *ierr);
void GET_REACTION_MAP(int *  rmap);
/*vector version */
void vproductionRate(int npt, double *  wdot, double *  c, double *  T);
void VCKHMS(int *  np, double *  T, double *  ums);
void VCKPY(int *  np, double *  rho, double *  T, double *  y, double *  P);
void VCKWYR(int *  np, double *  rho, double *  T,
            double *  y,
            double *  wdot);
void VCKYTX(int *  np, double *  y, double *  x);
void vcomp_k_f(int npt, double *  k_f_s, double *  tc, double *  invT);
void vcomp_gibbs(int npt, double *  g_RT, double *  tc);
void vcomp_Kc(int npt, double *  Kc_s, double *  g_RT, double *  invT);
void GET_CRITPARAMS(double *  Tci, double *  ai, double *  bi, double *  acentric_i);
void vcomp_wdot(int npt, double *  wdot, double *  mixture, double *  sc,
                double *  k_f_s, double *  Kc_s,
                double *  tc, double *  invT, double *  T);

/* Inverse molecular weights */
static const double imw[7] = {
    1.0 / 16.043030,  /*CH4 */
    1.0 / 31.998800,  /*O2 */
    1.0 / 18.015340,  /*H2O */
    1.0 / 28.013400,  /*N2 */
    1.0 / 28.010550,  /*CO */
    1.0 / 44.009950,  /*CO2 */
    1.0 / 2.015940};  /*H2 */



static double fwd_A[4], fwd_beta[4], fwd_Ea[4];
static double low_A[4], low_beta[4], low_Ea[4];
static double rev_A[4], rev_beta[4], rev_Ea[4];
static double troe_a[4],troe_Ts[4], troe_Tss[4], troe_Tsss[4];
static double sri_a[4], sri_b[4], sri_c[4], sri_d[4], sri_e[4];
static double activation_units[4], prefactor_units[4], phase_units[4];
static int is_PD[4], troe_len[4], sri_len[4], nTB[4], *TBid[4];
static double *TB[4];

static double fwd_A_DEF[4], fwd_beta_DEF[4], fwd_Ea_DEF[4];
static double low_A_DEF[4], low_beta_DEF[4], low_Ea_DEF[4];
static double rev_A_DEF[4], rev_beta_DEF[4], rev_Ea_DEF[4];
static double troe_a_DEF[4],troe_Ts_DEF[4], troe_Tss_DEF[4], troe_Tsss_DEF[4];
static double sri_a_DEF[4], sri_b_DEF[4], sri_c_DEF[4], sri_d_DEF[4], sri_e_DEF[4];
static double activation_units_DEF[4], prefactor_units_DEF[4], phase_units_DEF[4];
static int is_PD_DEF[4], troe_len_DEF[4], sri_len_DEF[4], nTB_DEF[4], *TBid_DEF[4];
static double *TB_DEF[4];
static int rxn_map[4] = {0,1,2,3};

void GET_REACTION_MAP(int *rmap)
{
    for (int i=0; i<4; ++i) {
        rmap[i] = rxn_map[i];
    }
}


#include <ReactionData.H>
double* GetParamPtr(int                reaction_id,
                    REACTION_PARAMETER param_id,
                    int                species_id,
                    int                get_default)
{
  double* ret = 0;
  if (reaction_id<0 || reaction_id>=4) {
    printf("Bad reaction id = %d",reaction_id);
    abort();
  };
  int mrid = rxn_map[reaction_id];

  if (param_id == THIRD_BODY) {
    if (species_id<0 || species_id>=7) {
      printf("GetParamPtr: Bad species id = %d",species_id);
      abort();
    }
    if (get_default) {
      for (int i=0; i<nTB_DEF[mrid]; ++i) {
        if (species_id == TBid_DEF[mrid][i]) {
          ret = &(TB_DEF[mrid][i]);
        }
      }
    }
    else {
      for (int i=0; i<nTB[mrid]; ++i) {
        if (species_id == TBid[mrid][i]) {
          ret = &(TB[mrid][i]);
        }
      }
    }
    if (ret == 0) {
      printf("GetParamPtr: No TB for reaction id = %d",reaction_id);
      abort();
    }
  }
  else {
    if (     param_id == FWD_A)     {ret = (get_default ? &(fwd_A_DEF[mrid]) : &(fwd_A[mrid]));}
      else if (param_id == FWD_BETA)  {ret = (get_default ? &(fwd_beta_DEF[mrid]) : &(fwd_beta[mrid]));}
      else if (param_id == FWD_EA)    {ret = (get_default ? &(fwd_Ea_DEF[mrid]) : &(fwd_Ea[mrid]));}
      else if (param_id == LOW_A)     {ret = (get_default ? &(low_A_DEF[mrid]) : &(low_A[mrid]));}
      else if (param_id == LOW_BETA)  {ret = (get_default ? &(low_beta_DEF[mrid]) : &(low_beta[mrid]));}
      else if (param_id == LOW_EA)    {ret = (get_default ? &(low_Ea_DEF[mrid]) : &(low_Ea[mrid]));}
      else if (param_id == REV_A)     {ret = (get_default ? &(rev_A_DEF[mrid]) : &(rev_A[mrid]));}
      else if (param_id == REV_BETA)  {ret = (get_default ? &(rev_beta_DEF[mrid]) : &(rev_beta[mrid]));}
      else if (param_id == REV_EA)    {ret = (get_default ? &(rev_Ea_DEF[mrid]) : &(rev_Ea[mrid]));}
      else if (param_id == TROE_A)    {ret = (get_default ? &(troe_a_DEF[mrid]) : &(troe_a[mrid]));}
      else if (param_id == TROE_TS)   {ret = (get_default ? &(troe_Ts_DEF[mrid]) : &(troe_Ts[mrid]));}
      else if (param_id == TROE_TSS)  {ret = (get_default ? &(troe_Tss_DEF[mrid]) : &(troe_Tss[mrid]));}
      else if (param_id == TROE_TSSS) {ret = (get_default ? &(troe_Tsss_DEF[mrid]) : &(troe_Tsss[mrid]));}
      else if (param_id == SRI_A)     {ret = (get_default ? &(sri_a_DEF[mrid]) : &(sri_a[mrid]));}
      else if (param_id == SRI_B)     {ret = (get_default ? &(sri_b_DEF[mrid]) : &(sri_b[mrid]));}
      else if (param_id == SRI_C)     {ret = (get_default ? &(sri_c_DEF[mrid]) : &(sri_c[mrid]));}
      else if (param_id == SRI_D)     {ret = (get_default ? &(sri_d_DEF[mrid]) : &(sri_d[mrid]));}
      else if (param_id == SRI_E)     {ret = (get_default ? &(sri_e_DEF[mrid]) : &(sri_e[mrid]));}
    else {
      printf("GetParamPtr: Unknown parameter id");
      abort();
    }
  }
  return ret;
}

void ResetAllParametersToDefault()
{
    for (int i=0; i<4; i++) {
        if (nTB[i] != 0) {
            nTB[i] = 0;
            free(TB[i]);
            free(TBid[i]);
        }

        fwd_A[i]    = fwd_A_DEF[i];
        fwd_beta[i] = fwd_beta_DEF[i];
        fwd_Ea[i]   = fwd_Ea_DEF[i];

        low_A[i]    = low_A_DEF[i];
        low_beta[i] = low_beta_DEF[i];
        low_Ea[i]   = low_Ea_DEF[i];

        rev_A[i]    = rev_A_DEF[i];
        rev_beta[i] = rev_beta_DEF[i];
        rev_Ea[i]   = rev_Ea_DEF[i];

        troe_a[i]    = troe_a_DEF[i];
        troe_Ts[i]   = troe_Ts_DEF[i];
        troe_Tss[i]  = troe_Tss_DEF[i];
        troe_Tsss[i] = troe_Tsss_DEF[i];

        sri_a[i] = sri_a_DEF[i];
        sri_b[i] = sri_b_DEF[i];
        sri_c[i] = sri_c_DEF[i];
        sri_d[i] = sri_d_DEF[i];
        sri_e[i] = sri_e_DEF[i];

        is_PD[i]    = is_PD_DEF[i];
        troe_len[i] = troe_len_DEF[i];
        sri_len[i]  = sri_len_DEF[i];

        activation_units[i] = activation_units_DEF[i];
        prefactor_units[i]  = prefactor_units_DEF[i];
        phase_units[i]      = phase_units_DEF[i];

        nTB[i]  = nTB_DEF[i];
        if (nTB[i] != 0) {
           TB[i] = (double *) malloc(sizeof(double) * nTB[i]);
           TBid[i] = (int *) malloc(sizeof(int) * nTB[i]);
           for (int j=0; j<nTB[i]; j++) {
             TB[i][j] = TB_DEF[i][j];
             TBid[i][j] = TBid_DEF[i][j];
           }
        }
    }
}

void SetAllDefaults()
{
    for (int i=0; i<4; i++) {
        if (nTB_DEF[i] != 0) {
            nTB_DEF[i] = 0;
            free(TB_DEF[i]);
            free(TBid_DEF[i]);
        }

        fwd_A_DEF[i]    = fwd_A[i];
        fwd_beta_DEF[i] = fwd_beta[i];
        fwd_Ea_DEF[i]   = fwd_Ea[i];

        low_A_DEF[i]    = low_A[i];
        low_beta_DEF[i] = low_beta[i];
        low_Ea_DEF[i]   = low_Ea[i];

        rev_A_DEF[i]    = rev_A[i];
        rev_beta_DEF[i] = rev_beta[i];
        rev_Ea_DEF[i]   = rev_Ea[i];

        troe_a_DEF[i]    = troe_a[i];
        troe_Ts_DEF[i]   = troe_Ts[i];
        troe_Tss_DEF[i]  = troe_Tss[i];
        troe_Tsss_DEF[i] = troe_Tsss[i];

        sri_a_DEF[i] = sri_a[i];
        sri_b_DEF[i] = sri_b[i];
        sri_c_DEF[i] = sri_c[i];
        sri_d_DEF[i] = sri_d[i];
        sri_e_DEF[i] = sri_e[i];

        is_PD_DEF[i]    = is_PD[i];
        troe_len_DEF[i] = troe_len[i];
        sri_len_DEF[i]  = sri_len[i];

        activation_units_DEF[i] = activation_units[i];
        prefactor_units_DEF[i]  = prefactor_units[i];
        phase_units_DEF[i]      = phase_units[i];

        nTB_DEF[i]  = nTB[i];
        if (nTB_DEF[i] != 0) {
           TB_DEF[i] = (double *) malloc(sizeof(double) * nTB_DEF[i]);
           TBid_DEF[i] = (int *) malloc(sizeof(int) * nTB_DEF[i]);
           for (int j=0; j<nTB_DEF[i]; j++) {
             TB_DEF[i][j] = TB[i][j];
             TBid_DEF[i][j] = TBid[i][j];
           }
        }
    }
}

/* Finalizes parameter database */
void CKFINALIZE()
{
  for (int i=0; i<4; ++i) {
    free(TB[i]); TB[i] = 0; 
    free(TBid[i]); TBid[i] = 0;
    nTB[i] = 0;

    free(TB_DEF[i]); TB_DEF[i] = 0; 
    free(TBid_DEF[i]); TBid_DEF[i] = 0;
    nTB_DEF[i] = 0;
  }
}

/* Initializes parameter database */
void CKINIT()
{
    // (0):  2 CH4 + O2 => 2 CO + 4 H2
    fwd_A[0]     = 39100000000000;
    fwd_beta[0]  = 0;
    fwd_Ea[0]    = 30000;
    prefactor_units[0]  = 1.0000000000000002e-12;
    activation_units[0] = 0.50321666580471969;
    phase_units[0]      = 1e-18;
    is_PD[0] = 0;
    nTB[0] = 0;

    // (1):  CH4 + H2O <=> CO + 3 H2
    fwd_A[1]     = 300000000000;
    fwd_beta[1]  = 0;
    fwd_Ea[1]    = 30000;
    prefactor_units[1]  = 1.0000000000000002e-06;
    activation_units[1] = 0.50321666580471969;
    phase_units[1]      = 1e-12;
    is_PD[1] = 0;
    nTB[1] = 0;

    // (2):  2 H2 + O2 => 2 H2O
    fwd_A[2]     = 6.045e+18;
    fwd_beta[2]  = -1;
    fwd_Ea[2]    = 40000;
    prefactor_units[2]  = 1.0000000000000002e-12;
    activation_units[2] = 0.50321666580471969;
    phase_units[2]      = 1e-18;
    is_PD[2] = 0;
    nTB[2] = 0;

    // (3):  CO + H2O <=> CO2 + H2
    fwd_A[3]     = 2750000000000;
    fwd_beta[3]  = 0;
    fwd_Ea[3]    = 20000;
    prefactor_units[3]  = 1.0000000000000002e-06;
    activation_units[3] = 0.50321666580471969;
    phase_units[3]      = 1e-12;
    is_PD[3] = 0;
    nTB[3] = 0;

    SetAllDefaults();
}



/*A few mechanism parameters */
void CKINDX(int * mm, int * kk, int * ii, int * nfit)
{
    *mm = 4;
    *kk = 7;
    *ii = 4;
    *nfit = -1; /*Why do you need this anyway ?  */
}



/* ckxnum... for parsing strings  */
void CKXNUM(char * line, int * nexp, int * lout, int * nval, double *  rval, int * kerr, int lenline )
{
    int n,i; /*Loop Counters */
    char cstr[1000];
    char *saveptr;
    char *p; /*String Tokens */
    /* Strip Comments  */
    for (i=0; i<lenline; ++i) {
        if (line[i]=='!') {
            break;
        }
        cstr[i] = line[i];
    }
    cstr[i] = '\0';

    p = strtok_r(cstr," ", &saveptr);
    if (!p) {
        *nval = 0;
        *kerr = 1;
        return;
    }
    for (n=0; n<*nexp; ++n) {
        rval[n] = atof(p);
        p = strtok_r(NULL, " ", &saveptr);
        if (!p) break;
    }
    *nval = n+1;
    if (*nval < *nexp) *kerr = 1;
    return;
}


/* cksnum... for parsing strings  */
void CKSNUM(char * line, int * nexp, int * lout, char * kray, int * nn, int * knum, int * nval, double *  rval, int * kerr, int lenline, int lenkray)
{
    /*Not done yet ... */
}


/* Returns the char strings of element names */
void CKSYME(int * kname, int * plenkname )
{
    int i; /*Loop Counter */
    int lenkname = *plenkname;
    /*clear kname */
    for (i=0; i<lenkname*4; i++) {
        kname[i] = ' ';
    }

    /* C  */
    kname[ 0*lenkname + 0 ] = 'C';
    kname[ 0*lenkname + 1 ] = ' ';

    /* O  */
    kname[ 1*lenkname + 0 ] = 'O';
    kname[ 1*lenkname + 1 ] = ' ';

    /* H  */
    kname[ 2*lenkname + 0 ] = 'H';
    kname[ 2*lenkname + 1 ] = ' ';

    /* N  */
    kname[ 3*lenkname + 0 ] = 'N';
    kname[ 3*lenkname + 1 ] = ' ';

}


/* Returns the char strings of species names */
void CKSYMS(int * kname, int * plenkname )
{
    int i; /*Loop Counter */
    int lenkname = *plenkname;
    /*clear kname */
    for (i=0; i<lenkname*7; i++) {
        kname[i] = ' ';
    }

    /* CH4  */
    kname[ 0*lenkname + 0 ] = 'C';
    kname[ 0*lenkname + 1 ] = 'H';
    kname[ 0*lenkname + 2 ] = '4';
    kname[ 0*lenkname + 3 ] = ' ';

    /* O2  */
    kname[ 1*lenkname + 0 ] = 'O';
    kname[ 1*lenkname + 1 ] = '2';
    kname[ 1*lenkname + 2 ] = ' ';

    /* H2O  */
    kname[ 2*lenkname + 0 ] = 'H';
    kname[ 2*lenkname + 1 ] = '2';
    kname[ 2*lenkname + 2 ] = 'O';
    kname[ 2*lenkname + 3 ] = ' ';

    /* N2  */
    kname[ 3*lenkname + 0 ] = 'N';
    kname[ 3*lenkname + 1 ] = '2';
    kname[ 3*lenkname + 2 ] = ' ';

    /* CO  */
    kname[ 4*lenkname + 0 ] = 'C';
    kname[ 4*lenkname + 1 ] = 'O';
    kname[ 4*lenkname + 2 ] = ' ';

    /* CO2  */
    kname[ 5*lenkname + 0 ] = 'C';
    kname[ 5*lenkname + 1 ] = 'O';
    kname[ 5*lenkname + 2 ] = '2';
    kname[ 5*lenkname + 3 ] = ' ';

    /* H2  */
    kname[ 6*lenkname + 0 ] = 'H';
    kname[ 6*lenkname + 1 ] = '2';
    kname[ 6*lenkname + 2 ] = ' ';

}


/* Returns R, Rc, Patm */
void CKRP(double *  ru, double *  ruc, double *  pa)
{
     *ru  = 8.31451e+07; 
     *ruc = 1.98721558317399615845; 
     *pa  = 1.01325e+06; 
}


/*Compute P = rhoRT/W(x) */
void CKPX(double *  rho, double *  T, double *  x, double *  P)
{
    double XW = 0;/* To hold mean molecular wt */
    XW += x[0]*16.043030; /*CH4 */
    XW += x[1]*31.998800; /*O2 */
    XW += x[2]*18.015340; /*H2O */
    XW += x[3]*28.013400; /*N2 */
    XW += x[4]*28.010550; /*CO */
    XW += x[5]*44.009950; /*CO2 */
    XW += x[6]*2.015940; /*H2 */
    *P = *rho * 8.31451e+07 * (*T) / XW; /*P = rho*R*T/W */

    return;
}


/*Compute P = rhoRT/W(y) */
void CKPY(double *  rho, double *  T, double *  y,  double *  P)
{
    double YOW = 0;/* for computing mean MW */
    YOW += y[0]*imw[0]; /*CH4 */
    YOW += y[1]*imw[1]; /*O2 */
    YOW += y[2]*imw[2]; /*H2O */
    YOW += y[3]*imw[3]; /*N2 */
    YOW += y[4]*imw[4]; /*CO */
    YOW += y[5]*imw[5]; /*CO2 */
    YOW += y[6]*imw[6]; /*H2 */
    *P = *rho * 8.31451e+07 * (*T) * YOW; /*P = rho*R*T/W */

    return;
}


/*Compute P = rhoRT/W(y) */
void VCKPY(int *  np, double *  rho, double *  T, double *  y,  double *  P)
{
    double YOW[*np];
    for (int i=0; i<(*np); i++) {
        YOW[i] = 0.0;
    }

    for (int n=0; n<7; n++) {
        for (int i=0; i<(*np); i++) {
            YOW[i] += y[n*(*np)+i] * imw[n];
        }
    }

    for (int i=0; i<(*np); i++) {
        P[i] = rho[i] * 8.31451e+07 * T[i] * YOW[i]; /*P = rho*R*T/W */
    }

    return;
}


/*Compute P = rhoRT/W(c) */
void CKPC(double *  rho, double *  T, double *  c,  double *  P)
{
    int id; /*loop counter */
    /*See Eq 5 in CK Manual */
    double W = 0;
    double sumC = 0;
    W += c[0]*16.043030; /*CH4 */
    W += c[1]*31.998800; /*O2 */
    W += c[2]*18.015340; /*H2O */
    W += c[3]*28.013400; /*N2 */
    W += c[4]*28.010550; /*CO */
    W += c[5]*44.009950; /*CO2 */
    W += c[6]*2.015940; /*H2 */

    for (id = 0; id < 7; ++id) {
        sumC += c[id];
    }
    *P = *rho * 8.31451e+07 * (*T) * sumC / W; /*P = rho*R*T/W */

    return;
}


/*Compute rho = PW(x)/RT */
void CKRHOX(double *  P, double *  T, double *  x,  double *  rho)
{
    double XW = 0;/* To hold mean molecular wt */
    XW += x[0]*16.043030; /*CH4 */
    XW += x[1]*31.998800; /*O2 */
    XW += x[2]*18.015340; /*H2O */
    XW += x[3]*28.013400; /*N2 */
    XW += x[4]*28.010550; /*CO */
    XW += x[5]*44.009950; /*CO2 */
    XW += x[6]*2.015940; /*H2 */
    *rho = *P * XW / (8.31451e+07 * (*T)); /*rho = P*W/(R*T) */

    return;
}


/*Compute rho = P*W(y)/RT */
void CKRHOY(double *  P, double *  T, double *  y,  double *  rho)
{
    double YOW = 0;
    double tmp[7];

    for (int i = 0; i < 7; i++)
    {
        tmp[i] = y[i]*imw[i];
    }
    for (int i = 0; i < 7; i++)
    {
        YOW += tmp[i];
    }

    *rho = *P / (8.31451e+07 * (*T) * YOW);/*rho = P*W/(R*T) */
    return;
}


/*Compute rho = P*W(c)/(R*T) */
void CKRHOC(double *  P, double *  T, double *  c,  double *  rho)
{
    int id; /*loop counter */
    /*See Eq 5 in CK Manual */
    double W = 0;
    double sumC = 0;
    W += c[0]*16.043030; /*CH4 */
    W += c[1]*31.998800; /*O2 */
    W += c[2]*18.015340; /*H2O */
    W += c[3]*28.013400; /*N2 */
    W += c[4]*28.010550; /*CO */
    W += c[5]*44.009950; /*CO2 */
    W += c[6]*2.015940; /*H2 */

    for (id = 0; id < 7; ++id) {
        sumC += c[id];
    }
    *rho = *P * W / (sumC * (*T) * 8.31451e+07); /*rho = PW/(R*T) */

    return;
}


/*get molecular weight for all species */
void CKWT( double *  wt)
{
    molecularWeight(wt);
}


/*get atomic weight for all elements */
void CKAWT( double *  awt)
{
    atomicWeight(awt);
}


/*given y[species]: mass fractions */
/*returns mean molecular weight (gm/mole) */
void CKMMWY(double *  y,  double *  wtm)
{
    double YOW = 0;
    double tmp[7];

    for (int i = 0; i < 7; i++)
    {
        tmp[i] = y[i]*imw[i];
    }
    for (int i = 0; i < 7; i++)
    {
        YOW += tmp[i];
    }

    *wtm = 1.0 / YOW;
    return;
}


/*given x[species]: mole fractions */
/*returns mean molecular weight (gm/mole) */
void CKMMWX(double *  x,  double *  wtm)
{
    double XW = 0;/* see Eq 4 in CK Manual */
    XW += x[0]*16.043030; /*CH4 */
    XW += x[1]*31.998800; /*O2 */
    XW += x[2]*18.015340; /*H2O */
    XW += x[3]*28.013400; /*N2 */
    XW += x[4]*28.010550; /*CO */
    XW += x[5]*44.009950; /*CO2 */
    XW += x[6]*2.015940; /*H2 */
    *wtm = XW;

    return;
}


/*given c[species]: molar concentration */
/*returns mean molecular weight (gm/mole) */
void CKMMWC(double *  c,  double *  wtm)
{
    int id; /*loop counter */
    /*See Eq 5 in CK Manual */
    double W = 0;
    double sumC = 0;
    W += c[0]*16.043030; /*CH4 */
    W += c[1]*31.998800; /*O2 */
    W += c[2]*18.015340; /*H2O */
    W += c[3]*28.013400; /*N2 */
    W += c[4]*28.010550; /*CO */
    W += c[5]*44.009950; /*CO2 */
    W += c[6]*2.015940; /*H2 */

    for (id = 0; id < 7; ++id) {
        sumC += c[id];
    }
    /* CK provides no guard against divison by zero */
    *wtm = W/sumC;

    return;
}


/*convert y[species] (mass fracs) to x[species] (mole fracs) */
void CKYTX(double *  y,  double *  x)
{
    double YOW = 0;
    double tmp[7];

    for (int i = 0; i < 7; i++)
    {
        tmp[i] = y[i]*imw[i];
    }
    for (int i = 0; i < 7; i++)
    {
        YOW += tmp[i];
    }

    double YOWINV = 1.0/YOW;

    for (int i = 0; i < 7; i++)
    {
        x[i] = y[i]*imw[i]*YOWINV;
    }
    return;
}


/*convert y[npoints*species] (mass fracs) to x[npoints*species] (mole fracs) */
void VCKYTX(int *  np, double *  y,  double *  x)
{
    double YOW[*np];
    for (int i=0; i<(*np); i++) {
        YOW[i] = 0.0;
    }

    for (int n=0; n<7; n++) {
        for (int i=0; i<(*np); i++) {
            x[n*(*np)+i] = y[n*(*np)+i] * imw[n];
            YOW[i] += x[n*(*np)+i];
        }
    }

    for (int i=0; i<(*np); i++) {
        YOW[i] = 1.0/YOW[i];
    }

    for (int n=0; n<7; n++) {
        for (int i=0; i<(*np); i++) {
            x[n*(*np)+i] *=  YOW[i];
        }
    }
}


/*convert y[species] (mass fracs) to c[species] (molar conc) */
void CKYTCP(double *  P, double *  T, double *  y,  double *  c)
{
    double YOW = 0;
    double PWORT;

    /*Compute inverse of mean molecular wt first */
    for (int i = 0; i < 7; i++)
    {
        c[i] = y[i]*imw[i];
    }
    for (int i = 0; i < 7; i++)
    {
        YOW += c[i];
    }

    /*PW/RT (see Eq. 7) */
    PWORT = (*P)/(YOW * 8.31451e+07 * (*T)); 
    /*Now compute conversion */

    for (int i = 0; i < 7; i++)
    {
        c[i] = PWORT * y[i] * imw[i];
    }
    return;
}


/*convert y[species] (mass fracs) to c[species] (molar conc) */
void CKYTCR(double *  rho, double *  T, double *  y,  double *  c)
{
    for (int i = 0; i < 7; i++)
    {
        c[i] = (*rho)  * y[i] * imw[i];
    }
}


/*convert x[species] (mole fracs) to y[species] (mass fracs) */
void CKXTY(double *  x,  double *  y)
{
    double XW = 0; /*See Eq 4, 9 in CK Manual */
    /*Compute mean molecular wt first */
    XW += x[0]*16.043030; /*CH4 */
    XW += x[1]*31.998800; /*O2 */
    XW += x[2]*18.015340; /*H2O */
    XW += x[3]*28.013400; /*N2 */
    XW += x[4]*28.010550; /*CO */
    XW += x[5]*44.009950; /*CO2 */
    XW += x[6]*2.015940; /*H2 */
    /*Now compute conversion */
    double XWinv = 1.0/XW;
    y[0] = x[0]*16.043030*XWinv; 
    y[1] = x[1]*31.998800*XWinv; 
    y[2] = x[2]*18.015340*XWinv; 
    y[3] = x[3]*28.013400*XWinv; 
    y[4] = x[4]*28.010550*XWinv; 
    y[5] = x[5]*44.009950*XWinv; 
    y[6] = x[6]*2.015940*XWinv; 

    return;
}


/*convert x[species] (mole fracs) to c[species] (molar conc) */
void CKXTCP(double *  P, double *  T, double *  x,  double *  c)
{
    int id; /*loop counter */
    double PORT = (*P)/(8.31451e+07 * (*T)); /*P/RT */

    /*Compute conversion, see Eq 10 */
    for (id = 0; id < 7; ++id) {
        c[id] = x[id]*PORT;
    }

    return;
}


/*convert x[species] (mole fracs) to c[species] (molar conc) */
void CKXTCR(double *  rho, double *  T, double *  x, double *  c)
{
    int id; /*loop counter */
    double XW = 0; /*See Eq 4, 11 in CK Manual */
    double ROW; 
    /*Compute mean molecular wt first */
    XW += x[0]*16.043030; /*CH4 */
    XW += x[1]*31.998800; /*O2 */
    XW += x[2]*18.015340; /*H2O */
    XW += x[3]*28.013400; /*N2 */
    XW += x[4]*28.010550; /*CO */
    XW += x[5]*44.009950; /*CO2 */
    XW += x[6]*2.015940; /*H2 */
    ROW = (*rho) / XW;

    /*Compute conversion, see Eq 11 */
    for (id = 0; id < 7; ++id) {
        c[id] = x[id]*ROW;
    }

    return;
}


/*convert c[species] (molar conc) to x[species] (mole fracs) */
void CKCTX(double *  c, double *  x)
{
    int id; /*loop counter */
    double sumC = 0; 

    /*compute sum of c  */
    for (id = 0; id < 7; ++id) {
        sumC += c[id];
    }

    /* See Eq 13  */
    double sumCinv = 1.0/sumC;
    for (id = 0; id < 7; ++id) {
        x[id] = c[id]*sumCinv;
    }

    return;
}


/*convert c[species] (molar conc) to y[species] (mass fracs) */
void CKCTY(double *  c, double *  y)
{
    double CW = 0; /*See Eq 12 in CK Manual */
    /*compute denominator in eq 12 first */
    CW += c[0]*16.043030; /*CH4 */
    CW += c[1]*31.998800; /*O2 */
    CW += c[2]*18.015340; /*H2O */
    CW += c[3]*28.013400; /*N2 */
    CW += c[4]*28.010550; /*CO */
    CW += c[5]*44.009950; /*CO2 */
    CW += c[6]*2.015940; /*H2 */
    /*Now compute conversion */
    double CWinv = 1.0/CW;
    y[0] = c[0]*16.043030*CWinv; 
    y[1] = c[1]*31.998800*CWinv; 
    y[2] = c[2]*18.015340*CWinv; 
    y[3] = c[3]*28.013400*CWinv; 
    y[4] = c[4]*28.010550*CWinv; 
    y[5] = c[5]*44.009950*CWinv; 
    y[6] = c[6]*2.015940*CWinv; 

    return;
}


/*get Cp/R as a function of T  */
/*for all species (Eq 19) */
void CKCPOR(double *  T, double *  cpor)
{
    double tT = *T; /*temporary temperature */
    double tc[] = { 0, tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    cp_R(cpor, tc);
}


/*get H/RT as a function of T  */
/*for all species (Eq 20) */
void CKHORT(double *  T, double *  hort)
{
    double tT = *T; /*temporary temperature */
    double tc[] = { 0, tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    speciesEnthalpy(hort, tc);
}


/*get S/R as a function of T  */
/*for all species (Eq 21) */
void CKSOR(double *  T, double *  sor)
{
    double tT = *T; /*temporary temperature */
    double tc[] = { log(tT), tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    speciesEntropy(sor, tc);
}


/*get specific heat at constant volume as a function  */
/*of T for all species (molar units) */
void CKCVML(double *  T,  double *  cvml)
{
    int id; /*loop counter */
    double tT = *T; /*temporary temperature */
    double tc[] = { 0, tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    cv_R(cvml, tc);

    /*convert to chemkin units */
    for (id = 0; id < 7; ++id) {
        cvml[id] *= 8.31451e+07;
    }
}


/*get specific heat at constant pressure as a  */
/*function of T for all species (molar units) */
void CKCPML(double *  T,  double *  cpml)
{
    int id; /*loop counter */
    double tT = *T; /*temporary temperature */
    double tc[] = { 0, tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    cp_R(cpml, tc);

    /*convert to chemkin units */
    for (id = 0; id < 7; ++id) {
        cpml[id] *= 8.31451e+07;
    }
}


/*get internal energy as a function  */
/*of T for all species (molar units) */
void CKUML(double *  T,  double *  uml)
{
    int id; /*loop counter */
    double tT = *T; /*temporary temperature */
    double tc[] = { 0, tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    double RT = 8.31451e+07*tT; /*R*T */
    speciesInternalEnergy(uml, tc);

    /*convert to chemkin units */
    for (id = 0; id < 7; ++id) {
        uml[id] *= RT;
    }
}


/*get enthalpy as a function  */
/*of T for all species (molar units) */
void CKHML(double *  T,  double *  hml)
{
    int id; /*loop counter */
    double tT = *T; /*temporary temperature */
    double tc[] = { 0, tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    double RT = 8.31451e+07*tT; /*R*T */
    speciesEnthalpy(hml, tc);

    /*convert to chemkin units */
    for (id = 0; id < 7; ++id) {
        hml[id] *= RT;
    }
}


/*get standard-state Gibbs energy as a function  */
/*of T for all species (molar units) */
void CKGML(double *  T,  double *  gml)
{
    int id; /*loop counter */
    double tT = *T; /*temporary temperature */
    double tc[] = { log(tT), tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    double RT = 8.31451e+07*tT; /*R*T */
    gibbs(gml, tc);

    /*convert to chemkin units */
    for (id = 0; id < 7; ++id) {
        gml[id] *= RT;
    }
}


/*get standard-state Helmholtz free energy as a  */
/*function of T for all species (molar units) */
void CKAML(double *  T,  double *  aml)
{
    int id; /*loop counter */
    double tT = *T; /*temporary temperature */
    double tc[] = { log(tT), tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    double RT = 8.31451e+07*tT; /*R*T */
    helmholtz(aml, tc);

    /*convert to chemkin units */
    for (id = 0; id < 7; ++id) {
        aml[id] *= RT;
    }
}


/*Returns the standard-state entropies in molar units */
void CKSML(double *  T,  double *  sml)
{
    int id; /*loop counter */
    double tT = *T; /*temporary temperature */
    double tc[] = { log(tT), tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    speciesEntropy(sml, tc);

    /*convert to chemkin units */
    for (id = 0; id < 7; ++id) {
        sml[id] *= 8.31451e+07;
    }
}


/*Returns the specific heats at constant volume */
/*in mass units (Eq. 29) */
void CKCVMS(double *  T,  double *  cvms)
{
    double tT = *T; /*temporary temperature */
    double tc[] = { 0, tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    cv_R(cvms, tc);
    /*multiply by R/molecularweight */
    cvms[0] *= 5.182630712527496e+06; /*CH4 */
    cvms[1] *= 2.598381814318037e+06; /*O2 */
    cvms[2] *= 4.615239012974499e+06; /*H2O */
    cvms[3] *= 2.968047434442088e+06; /*N2 */
    cvms[4] *= 2.968349425484326e+06; /*CO */
    cvms[5] *= 1.889234139098090e+06; /*CO2 */
    cvms[6] *= 4.124383662212169e+07; /*H2 */
}


/*Returns the specific heats at constant pressure */
/*in mass units (Eq. 26) */
void CKCPMS(double *  T,  double *  cpms)
{
    double tT = *T; /*temporary temperature */
    double tc[] = { 0, tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    cp_R(cpms, tc);
    /*multiply by R/molecularweight */
    cpms[0] *= 5.182630712527496e+06; /*CH4 */
    cpms[1] *= 2.598381814318037e+06; /*O2 */
    cpms[2] *= 4.615239012974499e+06; /*H2O */
    cpms[3] *= 2.968047434442088e+06; /*N2 */
    cpms[4] *= 2.968349425484326e+06; /*CO */
    cpms[5] *= 1.889234139098090e+06; /*CO2 */
    cpms[6] *= 4.124383662212169e+07; /*H2 */
}


/*Returns internal energy in mass units (Eq 30.) */
void CKUMS(double *  T,  double *  ums)
{
    double tT = *T; /*temporary temperature */
    double tc[] = { 0, tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    double RT = 8.31451e+07*tT; /*R*T */
    speciesInternalEnergy(ums, tc);
    for (int i = 0; i < 7; i++)
    {
        ums[i] *= RT*imw[i];
    }
}


/*Returns enthalpy in mass units (Eq 27.) */
void CKHMS(double *  T,  double *  hms)
{
    double tT = *T; /*temporary temperature */
    double tc[] = { 0, tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    double RT = 8.31451e+07*tT; /*R*T */
    speciesEnthalpy(hms, tc);
    for (int i = 0; i < 7; i++)
    {
        hms[i] *= RT*imw[i];
    }
}


/*Returns enthalpy in mass units (Eq 27.) */
void VCKHMS(int *  np, double *  T,  double *  hms)
{
    double tc[5], h[7];

    for (int i=0; i<(*np); i++) {
        tc[0] = 0.0;
        tc[1] = T[i];
        tc[2] = T[i]*T[i];
        tc[3] = T[i]*T[i]*T[i];
        tc[4] = T[i]*T[i]*T[i]*T[i];

        speciesEnthalpy(h, tc);

        hms[0*(*np)+i] = h[0];
        hms[1*(*np)+i] = h[1];
        hms[2*(*np)+i] = h[2];
        hms[3*(*np)+i] = h[3];
        hms[4*(*np)+i] = h[4];
        hms[5*(*np)+i] = h[5];
        hms[6*(*np)+i] = h[6];
    }

    for (int n=0; n<7; n++) {
        for (int i=0; i<(*np); i++) {
            hms[n*(*np)+i] *= 8.31451e+07 * T[i] * imw[n];
        }
    }
}


/*Returns gibbs in mass units (Eq 31.) */
void CKGMS(double *  T,  double *  gms)
{
    double tT = *T; /*temporary temperature */
    double tc[] = { log(tT), tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    double RT = 8.31451e+07*tT; /*R*T */
    gibbs(gms, tc);
    for (int i = 0; i < 7; i++)
    {
        gms[i] *= RT*imw[i];
    }
}


/*Returns helmholtz in mass units (Eq 32.) */
void CKAMS(double *  T,  double *  ams)
{
    double tT = *T; /*temporary temperature */
    double tc[] = { log(tT), tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    double RT = 8.31451e+07*tT; /*R*T */
    helmholtz(ams, tc);
    for (int i = 0; i < 7; i++)
    {
        ams[i] *= RT*imw[i];
    }
}


/*Returns the entropies in mass units (Eq 28.) */
void CKSMS(double *  T,  double *  sms)
{
    double tT = *T; /*temporary temperature */
    double tc[] = { log(tT), tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    speciesEntropy(sms, tc);
    /*multiply by R/molecularweight */
    sms[0] *= 5.182630712527496e+06; /*CH4 */
    sms[1] *= 2.598381814318037e+06; /*O2 */
    sms[2] *= 4.615239012974499e+06; /*H2O */
    sms[3] *= 2.968047434442088e+06; /*N2 */
    sms[4] *= 2.968349425484326e+06; /*CO */
    sms[5] *= 1.889234139098090e+06; /*CO2 */
    sms[6] *= 4.124383662212169e+07; /*H2 */
}


/*Returns the mean specific heat at CP (Eq. 33) */
void CKCPBL(double *  T, double *  x,  double *  cpbl)
{
    int id; /*loop counter */
    double result = 0; 
    double tT = *T; /*temporary temperature */
    double tc[] = { 0, tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    double cpor[7]; /* temporary storage */
    cp_R(cpor, tc);

    /*perform dot product */
    for (id = 0; id < 7; ++id) {
        result += x[id]*cpor[id];
    }

    *cpbl = result * 8.31451e+07;
}


/*Returns the mean specific heat at CP (Eq. 34) */
void CKCPBS(double *  T, double *  y,  double *  cpbs)
{
    double result = 0; 
    double tT = *T; /*temporary temperature */
    double tc[] = { 0, tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    double cpor[7], tresult[7]; /* temporary storage */
    cp_R(cpor, tc);
    for (int i = 0; i < 7; i++)
    {
        tresult[i] = cpor[i]*y[i]*imw[i];

    }
    for (int i = 0; i < 7; i++)
    {
        result += tresult[i];
    }

    *cpbs = result * 8.31451e+07;
}


/*Returns the mean specific heat at CV (Eq. 35) */
void CKCVBL(double *  T, double *  x,  double *  cvbl)
{
    int id; /*loop counter */
    double result = 0; 
    double tT = *T; /*temporary temperature */
    double tc[] = { 0, tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    double cvor[7]; /* temporary storage */
    cv_R(cvor, tc);

    /*perform dot product */
    for (id = 0; id < 7; ++id) {
        result += x[id]*cvor[id];
    }

    *cvbl = result * 8.31451e+07;
}


/*Returns the mean specific heat at CV (Eq. 36) */
void CKCVBS(double *  T, double *  y,  double *  cvbs)
{
    double result = 0; 
    double tT = *T; /*temporary temperature */
    double tc[] = { 0, tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    double cvor[7]; /* temporary storage */
    cv_R(cvor, tc);
    /*multiply by y/molecularweight */
    result += cvor[0]*y[0]*imw[0]; /*CH4 */
    result += cvor[1]*y[1]*imw[1]; /*O2 */
    result += cvor[2]*y[2]*imw[2]; /*H2O */
    result += cvor[3]*y[3]*imw[3]; /*N2 */
    result += cvor[4]*y[4]*imw[4]; /*CO */
    result += cvor[5]*y[5]*imw[5]; /*CO2 */
    result += cvor[6]*y[6]*imw[6]; /*H2 */

    *cvbs = result * 8.31451e+07;
}


/*Returns the mean enthalpy of the mixture in molar units */
void CKHBML(double *  T, double *  x,  double *  hbml)
{
    int id; /*loop counter */
    double result = 0; 
    double tT = *T; /*temporary temperature */
    double tc[] = { 0, tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    double hml[7]; /* temporary storage */
    double RT = 8.31451e+07*tT; /*R*T */
    speciesEnthalpy(hml, tc);

    /*perform dot product */
    for (id = 0; id < 7; ++id) {
        result += x[id]*hml[id];
    }

    *hbml = result * RT;
}


/*Returns mean enthalpy of mixture in mass units */
void CKHBMS(double *  T, double *  y,  double *  hbms)
{
    double result = 0;
    double tT = *T; /*temporary temperature */
    double tc[] = { 0, tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    double hml[7], tmp[7]; /* temporary storage */
    double RT = 8.31451e+07*tT; /*R*T */
    speciesEnthalpy(hml, tc);
    int id;
    for (id = 0; id < 7; ++id) {
        tmp[id] = y[id]*hml[id]*imw[id];
    }
    for (id = 0; id < 7; ++id) {
        result += tmp[id];
    }

    *hbms = result * RT;
}


/*get mean internal energy in molar units */
void CKUBML(double *  T, double *  x,  double *  ubml)
{
    int id; /*loop counter */
    double result = 0; 
    double tT = *T; /*temporary temperature */
    double tc[] = { 0, tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    double uml[7]; /* temporary energy array */
    double RT = 8.31451e+07*tT; /*R*T */
    speciesInternalEnergy(uml, tc);

    /*perform dot product */
    for (id = 0; id < 7; ++id) {
        result += x[id]*uml[id];
    }

    *ubml = result * RT;
}


/*get mean internal energy in mass units */
void CKUBMS(double *  T, double *  y,  double *  ubms)
{
    double result = 0;
    double tT = *T; /*temporary temperature */
    double tc[] = { 0, tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    double ums[7]; /* temporary energy array */
    double RT = 8.31451e+07*tT; /*R*T */
    speciesInternalEnergy(ums, tc);
    /*perform dot product + scaling by wt */
    result += y[0]*ums[0]*imw[0]; /*CH4 */
    result += y[1]*ums[1]*imw[1]; /*O2 */
    result += y[2]*ums[2]*imw[2]; /*H2O */
    result += y[3]*ums[3]*imw[3]; /*N2 */
    result += y[4]*ums[4]*imw[4]; /*CO */
    result += y[5]*ums[5]*imw[5]; /*CO2 */
    result += y[6]*ums[6]*imw[6]; /*H2 */

    *ubms = result * RT;
}


/*get mixture entropy in molar units */
void CKSBML(double *  P, double *  T, double *  x,  double *  sbml)
{
    int id; /*loop counter */
    double result = 0; 
    /*Log of normalized pressure in cgs units dynes/cm^2 by Patm */
    double logPratio = log ( *P / 1013250.0 ); 
    double tT = *T; /*temporary temperature */
    double tc[] = { log(tT), tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    double sor[7]; /* temporary storage */
    speciesEntropy(sor, tc);

    /*Compute Eq 42 */
    for (id = 0; id < 7; ++id) {
        result += x[id]*(sor[id]-log((x[id]+1e-100))-logPratio);
    }

    *sbml = result * 8.31451e+07;
}


/*get mixture entropy in mass units */
void CKSBMS(double *  P, double *  T, double *  y,  double *  sbms)
{
    double result = 0; 
    /*Log of normalized pressure in cgs units dynes/cm^2 by Patm */
    double logPratio = log ( *P / 1013250.0 ); 
    double tT = *T; /*temporary temperature */
    double tc[] = { log(tT), tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    double sor[7]; /* temporary storage */
    double x[7]; /* need a ytx conversion */
    double YOW = 0; /*See Eq 4, 6 in CK Manual */
    /*Compute inverse of mean molecular wt first */
    YOW += y[0]*imw[0]; /*CH4 */
    YOW += y[1]*imw[1]; /*O2 */
    YOW += y[2]*imw[2]; /*H2O */
    YOW += y[3]*imw[3]; /*N2 */
    YOW += y[4]*imw[4]; /*CO */
    YOW += y[5]*imw[5]; /*CO2 */
    YOW += y[6]*imw[6]; /*H2 */
    /*Now compute y to x conversion */
    x[0] = y[0]/(16.043030*YOW); 
    x[1] = y[1]/(31.998800*YOW); 
    x[2] = y[2]/(18.015340*YOW); 
    x[3] = y[3]/(28.013400*YOW); 
    x[4] = y[4]/(28.010550*YOW); 
    x[5] = y[5]/(44.009950*YOW); 
    x[6] = y[6]/(2.015940*YOW); 
    speciesEntropy(sor, tc);
    /*Perform computation in Eq 42 and 43 */
    result += x[0]*(sor[0]-log((x[0]+1e-100))-logPratio);
    result += x[1]*(sor[1]-log((x[1]+1e-100))-logPratio);
    result += x[2]*(sor[2]-log((x[2]+1e-100))-logPratio);
    result += x[3]*(sor[3]-log((x[3]+1e-100))-logPratio);
    result += x[4]*(sor[4]-log((x[4]+1e-100))-logPratio);
    result += x[5]*(sor[5]-log((x[5]+1e-100))-logPratio);
    result += x[6]*(sor[6]-log((x[6]+1e-100))-logPratio);
    /*Scale by R/W */
    *sbms = result * 8.31451e+07 * YOW;
}


/*Returns mean gibbs free energy in molar units */
void CKGBML(double *  P, double *  T, double *  x,  double *  gbml)
{
    int id; /*loop counter */
    double result = 0; 
    /*Log of normalized pressure in cgs units dynes/cm^2 by Patm */
    double logPratio = log ( *P / 1013250.0 ); 
    double tT = *T; /*temporary temperature */
    double tc[] = { log(tT), tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    double RT = 8.31451e+07*tT; /*R*T */
    double gort[7]; /* temporary storage */
    /*Compute g/RT */
    gibbs(gort, tc);

    /*Compute Eq 44 */
    for (id = 0; id < 7; ++id) {
        result += x[id]*(gort[id]+log((x[id]+1e-100))+logPratio);
    }

    *gbml = result * RT;
}


/*Returns mixture gibbs free energy in mass units */
void CKGBMS(double *  P, double *  T, double *  y,  double *  gbms)
{
    double result = 0; 
    /*Log of normalized pressure in cgs units dynes/cm^2 by Patm */
    double logPratio = log ( *P / 1013250.0 ); 
    double tT = *T; /*temporary temperature */
    double tc[] = { log(tT), tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    double RT = 8.31451e+07*tT; /*R*T */
    double gort[7]; /* temporary storage */
    double x[7]; /* need a ytx conversion */
    double YOW = 0; /*To hold 1/molecularweight */
    /*Compute inverse of mean molecular wt first */
    YOW += y[0]*imw[0]; /*CH4 */
    YOW += y[1]*imw[1]; /*O2 */
    YOW += y[2]*imw[2]; /*H2O */
    YOW += y[3]*imw[3]; /*N2 */
    YOW += y[4]*imw[4]; /*CO */
    YOW += y[5]*imw[5]; /*CO2 */
    YOW += y[6]*imw[6]; /*H2 */
    /*Now compute y to x conversion */
    x[0] = y[0]/(16.043030*YOW); 
    x[1] = y[1]/(31.998800*YOW); 
    x[2] = y[2]/(18.015340*YOW); 
    x[3] = y[3]/(28.013400*YOW); 
    x[4] = y[4]/(28.010550*YOW); 
    x[5] = y[5]/(44.009950*YOW); 
    x[6] = y[6]/(2.015940*YOW); 
    gibbs(gort, tc);
    /*Perform computation in Eq 44 */
    result += x[0]*(gort[0]+log((x[0]+1e-100))+logPratio);
    result += x[1]*(gort[1]+log((x[1]+1e-100))+logPratio);
    result += x[2]*(gort[2]+log((x[2]+1e-100))+logPratio);
    result += x[3]*(gort[3]+log((x[3]+1e-100))+logPratio);
    result += x[4]*(gort[4]+log((x[4]+1e-100))+logPratio);
    result += x[5]*(gort[5]+log((x[5]+1e-100))+logPratio);
    result += x[6]*(gort[6]+log((x[6]+1e-100))+logPratio);
    /*Scale by RT/W */
    *gbms = result * RT * YOW;
}


/*Returns mean helmholtz free energy in molar units */
void CKABML(double *  P, double *  T, double *  x,  double *  abml)
{
    int id; /*loop counter */
    double result = 0; 
    /*Log of normalized pressure in cgs units dynes/cm^2 by Patm */
    double logPratio = log ( *P / 1013250.0 ); 
    double tT = *T; /*temporary temperature */
    double tc[] = { log(tT), tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    double RT = 8.31451e+07*tT; /*R*T */
    double aort[7]; /* temporary storage */
    /*Compute g/RT */
    helmholtz(aort, tc);

    /*Compute Eq 44 */
    for (id = 0; id < 7; ++id) {
        result += x[id]*(aort[id]+log((x[id]+1e-100))+logPratio);
    }

    *abml = result * RT;
}


/*Returns mixture helmholtz free energy in mass units */
void CKABMS(double *  P, double *  T, double *  y,  double *  abms)
{
    double result = 0; 
    /*Log of normalized pressure in cgs units dynes/cm^2 by Patm */
    double logPratio = log ( *P / 1013250.0 ); 
    double tT = *T; /*temporary temperature */
    double tc[] = { log(tT), tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    double RT = 8.31451e+07*tT; /*R*T */
    double aort[7]; /* temporary storage */
    double x[7]; /* need a ytx conversion */
    double YOW = 0; /*To hold 1/molecularweight */
    /*Compute inverse of mean molecular wt first */
    YOW += y[0]*imw[0]; /*CH4 */
    YOW += y[1]*imw[1]; /*O2 */
    YOW += y[2]*imw[2]; /*H2O */
    YOW += y[3]*imw[3]; /*N2 */
    YOW += y[4]*imw[4]; /*CO */
    YOW += y[5]*imw[5]; /*CO2 */
    YOW += y[6]*imw[6]; /*H2 */
    /*Now compute y to x conversion */
    x[0] = y[0]/(16.043030*YOW); 
    x[1] = y[1]/(31.998800*YOW); 
    x[2] = y[2]/(18.015340*YOW); 
    x[3] = y[3]/(28.013400*YOW); 
    x[4] = y[4]/(28.010550*YOW); 
    x[5] = y[5]/(44.009950*YOW); 
    x[6] = y[6]/(2.015940*YOW); 
    helmholtz(aort, tc);
    /*Perform computation in Eq 44 */
    result += x[0]*(aort[0]+log((x[0]+1e-100))+logPratio);
    result += x[1]*(aort[1]+log((x[1]+1e-100))+logPratio);
    result += x[2]*(aort[2]+log((x[2]+1e-100))+logPratio);
    result += x[3]*(aort[3]+log((x[3]+1e-100))+logPratio);
    result += x[4]*(aort[4]+log((x[4]+1e-100))+logPratio);
    result += x[5]*(aort[5]+log((x[5]+1e-100))+logPratio);
    result += x[6]*(aort[6]+log((x[6]+1e-100))+logPratio);
    /*Scale by RT/W */
    *abms = result * RT * YOW;
}


/*compute the production rate for each species */
void CKWC(double *  T, double *  C,  double *  wdot)
{
    int id; /*loop counter */

    /*convert to SI */
    for (id = 0; id < 7; ++id) {
        C[id] *= 1.0e6;
    }

    /*convert to chemkin units */
    productionRate(wdot, C, *T);

    /*convert to chemkin units */
    for (id = 0; id < 7; ++id) {
        C[id] *= 1.0e-6;
        wdot[id] *= 1.0e-6;
    }
}


/*Returns the molar production rate of species */
/*Given P, T, and mass fractions */
void CKWYP(double *  P, double *  T, double *  y,  double *  wdot)
{
    int id; /*loop counter */
    double c[7]; /*temporary storage */
    double YOW = 0; 
    double PWORT; 
    /*Compute inverse of mean molecular wt first */
    YOW += y[0]*imw[0]; /*CH4 */
    YOW += y[1]*imw[1]; /*O2 */
    YOW += y[2]*imw[2]; /*H2O */
    YOW += y[3]*imw[3]; /*N2 */
    YOW += y[4]*imw[4]; /*CO */
    YOW += y[5]*imw[5]; /*CO2 */
    YOW += y[6]*imw[6]; /*H2 */
    /*PW/RT (see Eq. 7) */
    PWORT = (*P)/(YOW * 8.31451e+07 * (*T)); 
    /*multiply by 1e6 so c goes to SI */
    PWORT *= 1e6; 
    /*Now compute conversion (and go to SI) */
    c[0] = PWORT * y[0]*imw[0]; 
    c[1] = PWORT * y[1]*imw[1]; 
    c[2] = PWORT * y[2]*imw[2]; 
    c[3] = PWORT * y[3]*imw[3]; 
    c[4] = PWORT * y[4]*imw[4]; 
    c[5] = PWORT * y[5]*imw[5]; 
    c[6] = PWORT * y[6]*imw[6]; 

    /*convert to chemkin units */
    productionRate(wdot, c, *T);

    /*convert to chemkin units */
    for (id = 0; id < 7; ++id) {
        wdot[id] *= 1.0e-6;
    }
}


/*Returns the molar production rate of species */
/*Given P, T, and mole fractions */
void CKWXP(double *  P, double *  T, double *  x,  double *  wdot)
{
    int id; /*loop counter */
    double c[7]; /*temporary storage */
    double PORT = 1e6 * (*P)/(8.31451e+07 * (*T)); /*1e6 * P/RT so c goes to SI units */

    /*Compute conversion, see Eq 10 */
    for (id = 0; id < 7; ++id) {
        c[id] = x[id]*PORT;
    }

    /*convert to chemkin units */
    productionRate(wdot, c, *T);

    /*convert to chemkin units */
    for (id = 0; id < 7; ++id) {
        wdot[id] *= 1.0e-6;
    }
}


/*Returns the molar production rate of species */
/*Given rho, T, and mass fractions */
void CKWYR(double *  rho, double *  T, double *  y,  double *  wdot)
{
    int id; /*loop counter */
    double c[7]; /*temporary storage */
    /*See Eq 8 with an extra 1e6 so c goes to SI */
    c[0] = 1e6 * (*rho) * y[0]*imw[0]; 
    c[1] = 1e6 * (*rho) * y[1]*imw[1]; 
    c[2] = 1e6 * (*rho) * y[2]*imw[2]; 
    c[3] = 1e6 * (*rho) * y[3]*imw[3]; 
    c[4] = 1e6 * (*rho) * y[4]*imw[4]; 
    c[5] = 1e6 * (*rho) * y[5]*imw[5]; 
    c[6] = 1e6 * (*rho) * y[6]*imw[6]; 

    /*call productionRate */
    productionRate(wdot, c, *T);

    /*convert to chemkin units */
    for (id = 0; id < 7; ++id) {
        wdot[id] *= 1.0e-6;
    }
}


/*Returns the molar production rate of species */
/*Given rho, T, and mass fractions */
void VCKWYR(int *  np, double *  rho, double *  T,
	    double *  y,
	    double *  wdot)
{
    double c[7*(*np)]; /*temporary storage */
    /*See Eq 8 with an extra 1e6 so c goes to SI */
    for (int n=0; n<7; n++) {
        for (int i=0; i<(*np); i++) {
            c[n*(*np)+i] = 1.0e6 * rho[i] * y[n*(*np)+i] * imw[n];
        }
    }

    /*call productionRate */
    vproductionRate(*np, wdot, c, T);

    /*convert to chemkin units */
    for (int i=0; i<7*(*np); i++) {
        wdot[i] *= 1.0e-6;
    }
}


/*Returns the molar production rate of species */
/*Given rho, T, and mole fractions */
void CKWXR(double *  rho, double *  T, double *  x,  double *  wdot)
{
    int id; /*loop counter */
    double c[7]; /*temporary storage */
    double XW = 0; /*See Eq 4, 11 in CK Manual */
    double ROW; 
    /*Compute mean molecular wt first */
    XW += x[0]*16.043030; /*CH4 */
    XW += x[1]*31.998800; /*O2 */
    XW += x[2]*18.015340; /*H2O */
    XW += x[3]*28.013400; /*N2 */
    XW += x[4]*28.010550; /*CO */
    XW += x[5]*44.009950; /*CO2 */
    XW += x[6]*2.015940; /*H2 */
    /*Extra 1e6 factor to take c to SI */
    ROW = 1e6*(*rho) / XW;

    /*Compute conversion, see Eq 11 */
    for (id = 0; id < 7; ++id) {
        c[id] = x[id]*ROW;
    }

    /*convert to chemkin units */
    productionRate(wdot, c, *T);

    /*convert to chemkin units */
    for (id = 0; id < 7; ++id) {
        wdot[id] *= 1.0e-6;
    }
}


/*Returns the rate of progress for each reaction */
void CKQC(double *  T, double *  C, double *  qdot)
{
    int id; /*loop counter */

    /*convert to SI */
    for (id = 0; id < 7; ++id) {
        C[id] *= 1.0e6;
    }

    /*convert to chemkin units */
    progressRate(qdot, C, *T);

    /*convert to chemkin units */
    for (id = 0; id < 7; ++id) {
        C[id] *= 1.0e-6;
    }

    for (id = 0; id < 4; ++id) {
        qdot[id] *= 1.0e-6;
    }
}


/*Returns the progress rates of each reactions */
/*Given P, T, and mole fractions */
void CKKFKR(double *  P, double *  T, double *  x, double *  q_f, double *  q_r)
{
    int id; /*loop counter */
    double c[7]; /*temporary storage */
    double PORT = 1e6 * (*P)/(8.31451e+07 * (*T)); /*1e6 * P/RT so c goes to SI units */

    /*Compute conversion, see Eq 10 */
    for (id = 0; id < 7; ++id) {
        c[id] = x[id]*PORT;
    }

    /*convert to chemkin units */
    progressRateFR(q_f, q_r, c, *T);

    /*convert to chemkin units */
    for (id = 0; id < 4; ++id) {
        q_f[id] *= 1.0e-6;
        q_r[id] *= 1.0e-6;
    }
}


/*Returns the progress rates of each reactions */
/*Given P, T, and mass fractions */
void CKQYP(double *  P, double *  T, double *  y, double *  qdot)
{
    int id; /*loop counter */
    double c[7]; /*temporary storage */
    double YOW = 0; 
    double PWORT; 
    /*Compute inverse of mean molecular wt first */
    YOW += y[0]*imw[0]; /*CH4 */
    YOW += y[1]*imw[1]; /*O2 */
    YOW += y[2]*imw[2]; /*H2O */
    YOW += y[3]*imw[3]; /*N2 */
    YOW += y[4]*imw[4]; /*CO */
    YOW += y[5]*imw[5]; /*CO2 */
    YOW += y[6]*imw[6]; /*H2 */
    /*PW/RT (see Eq. 7) */
    PWORT = (*P)/(YOW * 8.31451e+07 * (*T)); 
    /*multiply by 1e6 so c goes to SI */
    PWORT *= 1e6; 
    /*Now compute conversion (and go to SI) */
    c[0] = PWORT * y[0]*imw[0]; 
    c[1] = PWORT * y[1]*imw[1]; 
    c[2] = PWORT * y[2]*imw[2]; 
    c[3] = PWORT * y[3]*imw[3]; 
    c[4] = PWORT * y[4]*imw[4]; 
    c[5] = PWORT * y[5]*imw[5]; 
    c[6] = PWORT * y[6]*imw[6]; 

    /*convert to chemkin units */
    progressRate(qdot, c, *T);

    /*convert to chemkin units */
    for (id = 0; id < 4; ++id) {
        qdot[id] *= 1.0e-6;
    }
}


/*Returns the progress rates of each reactions */
/*Given P, T, and mole fractions */
void CKQXP(double *  P, double *  T, double *  x, double *  qdot)
{
    int id; /*loop counter */
    double c[7]; /*temporary storage */
    double PORT = 1e6 * (*P)/(8.31451e+07 * (*T)); /*1e6 * P/RT so c goes to SI units */

    /*Compute conversion, see Eq 10 */
    for (id = 0; id < 7; ++id) {
        c[id] = x[id]*PORT;
    }

    /*convert to chemkin units */
    progressRate(qdot, c, *T);

    /*convert to chemkin units */
    for (id = 0; id < 4; ++id) {
        qdot[id] *= 1.0e-6;
    }
}


/*Returns the progress rates of each reactions */
/*Given rho, T, and mass fractions */
void CKQYR(double *  rho, double *  T, double *  y, double *  qdot)
{
    int id; /*loop counter */
    double c[7]; /*temporary storage */
    /*See Eq 8 with an extra 1e6 so c goes to SI */
    c[0] = 1e6 * (*rho) * y[0]*imw[0]; 
    c[1] = 1e6 * (*rho) * y[1]*imw[1]; 
    c[2] = 1e6 * (*rho) * y[2]*imw[2]; 
    c[3] = 1e6 * (*rho) * y[3]*imw[3]; 
    c[4] = 1e6 * (*rho) * y[4]*imw[4]; 
    c[5] = 1e6 * (*rho) * y[5]*imw[5]; 
    c[6] = 1e6 * (*rho) * y[6]*imw[6]; 

    /*call progressRate */
    progressRate(qdot, c, *T);

    /*convert to chemkin units */
    for (id = 0; id < 4; ++id) {
        qdot[id] *= 1.0e-6;
    }
}


/*Returns the progress rates of each reactions */
/*Given rho, T, and mole fractions */
void CKQXR(double *  rho, double *  T, double *  x, double *  qdot)
{
    int id; /*loop counter */
    double c[7]; /*temporary storage */
    double XW = 0; /*See Eq 4, 11 in CK Manual */
    double ROW; 
    /*Compute mean molecular wt first */
    XW += x[0]*16.043030; /*CH4 */
    XW += x[1]*31.998800; /*O2 */
    XW += x[2]*18.015340; /*H2O */
    XW += x[3]*28.013400; /*N2 */
    XW += x[4]*28.010550; /*CO */
    XW += x[5]*44.009950; /*CO2 */
    XW += x[6]*2.015940; /*H2 */
    /*Extra 1e6 factor to take c to SI */
    ROW = 1e6*(*rho) / XW;

    /*Compute conversion, see Eq 11 */
    for (id = 0; id < 7; ++id) {
        c[id] = x[id]*ROW;
    }

    /*convert to chemkin units */
    progressRate(qdot, c, *T);

    /*convert to chemkin units */
    for (id = 0; id < 4; ++id) {
        qdot[id] *= 1.0e-6;
    }
}


/*Returns the stoichiometric coefficients */
/*of the reaction mechanism. (Eq 50) */
void CKNU(int * kdim,  int * nuki)
{
    int id; /*loop counter */
    int kd = (*kdim); 
    /*Zero nuki */
    for (id = 0; id < 7 * kd; ++ id) {
         nuki[id] = 0; 
    }

    /*reaction 1: 2 CH4 + O2 => 2 CO + 4 H2 */
    nuki[ 0 * kd + 0 ] += -2 ;
    nuki[ 1 * kd + 0 ] += -1 ;
    nuki[ 4 * kd + 0 ] += +2 ;
    nuki[ 6 * kd + 0 ] += +4 ;

    /*reaction 2: CH4 + H2O <=> CO + 3 H2 */
    nuki[ 0 * kd + 1 ] += -1 ;
    nuki[ 2 * kd + 1 ] += -1 ;
    nuki[ 4 * kd + 1 ] += +1 ;
    nuki[ 6 * kd + 1 ] += +3 ;

    /*reaction 3: 2 H2 + O2 => 2 H2O */
    nuki[ 6 * kd + 2 ] += -2 ;
    nuki[ 1 * kd + 2 ] += -1 ;
    nuki[ 2 * kd + 2 ] += +2 ;

    /*reaction 4: CO + H2O <=> CO2 + H2 */
    nuki[ 4 * kd + 3 ] += -1 ;
    nuki[ 2 * kd + 3 ] += -1 ;
    nuki[ 5 * kd + 3 ] += +1 ;
    nuki[ 6 * kd + 3 ] += +1 ;
}


/*Returns the elemental composition  */
/*of the speciesi (mdim is num of elements) */
void CKNCF(int * mdim,  int * ncf)
{
    int id; /*loop counter */
    int kd = (*mdim); 
    /*Zero ncf */
    for (id = 0; id < kd * 7; ++ id) {
         ncf[id] = 0; 
    }

    /*CH4 */
    ncf[ 0 * kd + 0 ] = 1; /*C */
    ncf[ 0 * kd + 2 ] = 4; /*H */

    /*O2 */
    ncf[ 1 * kd + 1 ] = 2; /*O */

    /*H2O */
    ncf[ 2 * kd + 2 ] = 2; /*H */
    ncf[ 2 * kd + 1 ] = 1; /*O */

    /*N2 */
    ncf[ 3 * kd + 3 ] = 2; /*N */

    /*CO */
    ncf[ 4 * kd + 0 ] = 1; /*C */
    ncf[ 4 * kd + 1 ] = 1; /*O */

    /*CO2 */
    ncf[ 5 * kd + 0 ] = 1; /*C */
    ncf[ 5 * kd + 1 ] = 2; /*O */

    /*H2 */
    ncf[ 6 * kd + 2 ] = 2; /*H */

}


/*Returns the arrehenius coefficients  */
/*for all reactions */
void CKABE( double *  a, double *  b, double *  e)
{
    for (int i=0; i<4; ++i) {
        a[i] = fwd_A[i];
        b[i] = fwd_beta[i];
        e[i] = fwd_Ea[i];
    }

    return;
}


/*Returns the equil constants for each reaction */
void CKEQC(double *  T, double *  C, double *  eqcon)
{
    double tT = *T; /*temporary temperature */
    double tc[] = { log(tT), tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    double gort[7]; /* temporary storage */

    /*compute the Gibbs free energy */
    gibbs(gort, tc);

    /*compute the equilibrium constants */
    equilibriumConstants(eqcon, gort, tT);

    /*reaction 1: 2 CH4 + O2 => 2 CO + 4 H2 */
    eqcon[0] *= 1e-18; 

    /*reaction 2: CH4 + H2O <=> CO + 3 H2 */
    eqcon[1] *= 1e-12; 

    /*reaction 3: 2 H2 + O2 => 2 H2O */
    eqcon[2] *= 1e+06; 

    /*reaction 4: CO + H2O <=> CO2 + H2 */
    /*eqcon[3] *= 1;  */
}


/*Returns the equil constants for each reaction */
/*Given P, T, and mass fractions */
void CKEQYP(double *  P, double *  T, double *  y, double *  eqcon)
{
    double tT = *T; /*temporary temperature */
    double tc[] = { log(tT), tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    double gort[7]; /* temporary storage */

    /*compute the Gibbs free energy */
    gibbs(gort, tc);

    /*compute the equilibrium constants */
    equilibriumConstants(eqcon, gort, tT);

    /*reaction 1: 2 CH4 + O2 => 2 CO + 4 H2 */
    eqcon[0] *= 1e-18; 

    /*reaction 2: CH4 + H2O <=> CO + 3 H2 */
    eqcon[1] *= 1e-12; 

    /*reaction 3: 2 H2 + O2 => 2 H2O */
    eqcon[2] *= 1e+06; 

    /*reaction 4: CO + H2O <=> CO2 + H2 */
    /*eqcon[3] *= 1;  */
}


/*Returns the equil constants for each reaction */
/*Given P, T, and mole fractions */
void CKEQXP(double *  P, double *  T, double *  x, double *  eqcon)
{
    double tT = *T; /*temporary temperature */
    double tc[] = { log(tT), tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    double gort[7]; /* temporary storage */

    /*compute the Gibbs free energy */
    gibbs(gort, tc);

    /*compute the equilibrium constants */
    equilibriumConstants(eqcon, gort, tT);

    /*reaction 1: 2 CH4 + O2 => 2 CO + 4 H2 */
    eqcon[0] *= 1e-18; 

    /*reaction 2: CH4 + H2O <=> CO + 3 H2 */
    eqcon[1] *= 1e-12; 

    /*reaction 3: 2 H2 + O2 => 2 H2O */
    eqcon[2] *= 1e+06; 

    /*reaction 4: CO + H2O <=> CO2 + H2 */
    /*eqcon[3] *= 1;  */
}


/*Returns the equil constants for each reaction */
/*Given rho, T, and mass fractions */
void CKEQYR(double *  rho, double *  T, double *  y, double *  eqcon)
{
    double tT = *T; /*temporary temperature */
    double tc[] = { log(tT), tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    double gort[7]; /* temporary storage */

    /*compute the Gibbs free energy */
    gibbs(gort, tc);

    /*compute the equilibrium constants */
    equilibriumConstants(eqcon, gort, tT);

    /*reaction 1: 2 CH4 + O2 => 2 CO + 4 H2 */
    eqcon[0] *= 1e-18; 

    /*reaction 2: CH4 + H2O <=> CO + 3 H2 */
    eqcon[1] *= 1e-12; 

    /*reaction 3: 2 H2 + O2 => 2 H2O */
    eqcon[2] *= 1e+06; 

    /*reaction 4: CO + H2O <=> CO2 + H2 */
    /*eqcon[3] *= 1;  */
}


/*Returns the equil constants for each reaction */
/*Given rho, T, and mole fractions */
void CKEQXR(double *  rho, double *  T, double *  x, double *  eqcon)
{
    double tT = *T; /*temporary temperature */
    double tc[] = { log(tT), tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    double gort[7]; /* temporary storage */

    /*compute the Gibbs free energy */
    gibbs(gort, tc);

    /*compute the equilibrium constants */
    equilibriumConstants(eqcon, gort, tT);

    /*reaction 1: 2 CH4 + O2 => 2 CO + 4 H2 */
    eqcon[0] *= 1e-18; 

    /*reaction 2: CH4 + H2O <=> CO + 3 H2 */
    eqcon[1] *= 1e-12; 

    /*reaction 3: 2 H2 + O2 => 2 H2O */
    eqcon[2] *= 1e+06; 

    /*reaction 4: CO + H2O <=> CO2 + H2 */
    /*eqcon[3] *= 1;  */
}

static double T_save = -1;
#ifdef _OPENMP
#pragma omp threadprivate(T_save)
#endif

static double k_f_save[4];
#ifdef _OPENMP
#pragma omp threadprivate(k_f_save)
#endif

static double Kc_save[4];
#ifdef _OPENMP
#pragma omp threadprivate(Kc_save)
#endif


/*compute the production rate for each species */
void productionRate(double *  wdot, double *  sc, double T)
{
    double tc[] = { log(T), T, T*T, T*T*T, T*T*T*T }; /*temperature cache */
    double invT = 1.0 / tc[1];

    if (T != T_save)
    {
        T_save = T;
        comp_k_f(tc,invT,k_f_save);
        comp_Kc(tc,invT,Kc_save);
    }

    double qdot, q_f[4], q_r[4];
    comp_qfqr(q_f, q_r, sc, tc, invT);

    for (int i = 0; i < 7; ++i) {
        wdot[i] = 0.0;
    }

    qdot = q_f[0]-q_r[0];
    wdot[0] -= 2 * qdot;
    wdot[1] -= qdot;
    wdot[4] += 2 * qdot;
    wdot[6] += 4 * qdot;

    qdot = q_f[1]-q_r[1];
    wdot[0] -= qdot;
    wdot[2] -= qdot;
    wdot[4] += qdot;
    wdot[6] += 3 * qdot;

    qdot = q_f[2]-q_r[2];
    wdot[1] -= qdot;
    wdot[2] += 2 * qdot;
    wdot[6] -= 2 * qdot;

    qdot = q_f[3]-q_r[3];
    wdot[2] -= qdot;
    wdot[4] -= qdot;
    wdot[5] += qdot;
    wdot[6] += qdot;

    return;
}

void comp_k_f(double *  tc, double invT, double *  k_f)
{
#ifdef __INTEL_COMPILER
    #pragma simd
#endif
    for (int i=0; i<4; ++i) {
        k_f[i] = prefactor_units[i] * fwd_A[i]
                    * exp(fwd_beta[i] * tc[0] - activation_units[i] * fwd_Ea[i] * invT);
    };
    return;
}

void comp_Kc(double *  tc, double invT, double *  Kc)
{
    /*compute the Gibbs free energy */
    double g_RT[7];
    gibbs(g_RT, tc);

    Kc[0] = 2*g_RT[0] + g_RT[1] - 2*g_RT[4] - 4*g_RT[6];
    Kc[1] = g_RT[0] + g_RT[2] - g_RT[4] - 3*g_RT[6];
    Kc[2] = g_RT[1] - 2*g_RT[2] + 2*g_RT[6];
    Kc[3] = g_RT[2] + g_RT[4] - g_RT[5] - g_RT[6];

#ifdef __INTEL_COMPILER
     #pragma simd
#endif
    for (int i=0; i<4; ++i) {
        Kc[i] = exp(Kc[i]);
    };

    /*reference concentration: P_atm / (RT) in inverse mol/m^3 */
    double refC = 101325 / 8.31451 * invT;
    double refCinv = 1 / refC;

    Kc[0] *= refC*refC*refC;
    Kc[1] *= refC*refC;
    Kc[2] *= refCinv;

    return;
}

void comp_qfqr(double *  qf, double *  qr, double *  sc, double *  tc, double invT)
{

    /*reaction 1: 2 CH4 + O2 => 2 CO + 4 H2 */
    qf[0] = sc[0]*sc[0]*sc[1];
    qr[0] = 0.0;

    /*reaction 2: CH4 + H2O <=> CO + 3 H2 */
    qf[1] = sc[0]*sc[2];
    qr[1] = sc[4]*sc[6]*sc[6]*sc[6];

    /*reaction 3: 2 H2 + O2 => 2 H2O */
    qf[2] = sc[1]*sc[6]*sc[6];
    qr[2] = 0.0;

    /*reaction 4: CO + H2O <=> CO2 + H2 */
    qf[3] = sc[2]*sc[4];
    qr[3] = sc[5]*sc[6];

    double T = tc[1];

    /*compute the mixture concentration */
    double mixture = 0.0;
    for (int i = 0; i < 7; ++i) {
        mixture += sc[i];
    }

    double Corr[4];
    for (int i = 0; i < 4; ++i) {
        Corr[i] = 1.0;
    }

    for (int i=0; i<4; i++)
    {
        qf[i] *= Corr[i] * k_f_save[i];
        qr[i] *= Corr[i] * k_f_save[i] / Kc_save[i];
    }

    return;
}


/*compute the production rate for each species */
void vproductionRate(int npt, double *  wdot, double *  sc, double *  T)
{
    double k_f_s[4*npt], Kc_s[4*npt], mixture[npt], g_RT[7*npt];
    double tc[5*npt], invT[npt];

#ifdef __INTEL_COMPILER
     #pragma simd
#endif
    for (int i=0; i<npt; i++) {
        tc[0*npt+i] = log(T[i]);
        tc[1*npt+i] = T[i];
        tc[2*npt+i] = T[i]*T[i];
        tc[3*npt+i] = T[i]*T[i]*T[i];
        tc[4*npt+i] = T[i]*T[i]*T[i]*T[i];
        invT[i] = 1.0 / T[i];
    }

    for (int i=0; i<npt; i++) {
        mixture[i] = 0.0;
    }

    for (int n=0; n<7; n++) {
        for (int i=0; i<npt; i++) {
            mixture[i] += sc[n*npt+i];
            wdot[n*npt+i] = 0.0;
        }
    }

    vcomp_k_f(npt, k_f_s, tc, invT);

    vcomp_gibbs(npt, g_RT, tc);

    vcomp_Kc(npt, Kc_s, g_RT, invT);

    vcomp_wdot(npt, wdot, mixture, sc, k_f_s, Kc_s, tc, invT, T);
}

void vcomp_k_f(int npt, double *  k_f_s, double *  tc, double *  invT)
{
#ifdef __INTEL_COMPILER
    #pragma simd
#endif
    for (int i=0; i<npt; i++) {
        k_f_s[0*npt+i] = prefactor_units[0] * fwd_A[0] * exp(fwd_beta[0] * tc[i] - activation_units[0] * fwd_Ea[0] * invT[i]);
        k_f_s[1*npt+i] = prefactor_units[1] * fwd_A[1] * exp(fwd_beta[1] * tc[i] - activation_units[1] * fwd_Ea[1] * invT[i]);
        k_f_s[2*npt+i] = prefactor_units[2] * fwd_A[2] * exp(fwd_beta[2] * tc[i] - activation_units[2] * fwd_Ea[2] * invT[i]);
        k_f_s[3*npt+i] = prefactor_units[3] * fwd_A[3] * exp(fwd_beta[3] * tc[i] - activation_units[3] * fwd_Ea[3] * invT[i]);
    }
}

void vcomp_gibbs(int npt, double *  g_RT, double *  tc)
{
    /*compute the Gibbs free energy */
    for (int i=0; i<npt; i++) {
        double tg[5], g[7];
        tg[0] = tc[0*npt+i];
        tg[1] = tc[1*npt+i];
        tg[2] = tc[2*npt+i];
        tg[3] = tc[3*npt+i];
        tg[4] = tc[4*npt+i];

        gibbs(g, tg);

        g_RT[0*npt+i] = g[0];
        g_RT[1*npt+i] = g[1];
        g_RT[2*npt+i] = g[2];
        g_RT[3*npt+i] = g[3];
        g_RT[4*npt+i] = g[4];
        g_RT[5*npt+i] = g[5];
        g_RT[6*npt+i] = g[6];
    }
}

void vcomp_Kc(int npt, double *  Kc_s, double *  g_RT, double *  invT)
{
#ifdef __INTEL_COMPILER
    #pragma simd
#endif
    for (int i=0; i<npt; i++) {
        /*reference concentration: P_atm / (RT) in inverse mol/m^3 */
        double refC = (101325. / 8.31451) * invT[i];
        double refCinv = 1.0 / refC;

        Kc_s[0*npt+i] = refC*refC*refC * exp((2 * g_RT[0*npt+i] + g_RT[1*npt+i]) - (2 * g_RT[4*npt+i] + 4 * g_RT[6*npt+i]));
        Kc_s[1*npt+i] = refC*refC * exp((g_RT[0*npt+i] + g_RT[2*npt+i]) - (g_RT[4*npt+i] + 3 * g_RT[6*npt+i]));
        Kc_s[2*npt+i] = refCinv * exp((g_RT[1*npt+i] + 2 * g_RT[6*npt+i]) - (2 * g_RT[2*npt+i]));
        Kc_s[3*npt+i] = exp((g_RT[2*npt+i] + g_RT[4*npt+i]) - (g_RT[5*npt+i] + g_RT[6*npt+i]));
    }
}

void vcomp_wdot(int npt, double *  wdot, double *  mixture, double *  sc,
		double *  k_f_s, double *  Kc_s,
		double *  tc, double *  invT, double *  T)
{
#ifdef __INTEL_COMPILER
    #pragma simd
#endif
    for (int i=0; i<npt; i++) {
        double qdot, q_f, q_r, phi_f, phi_r, k_f, k_r, Kc;

        /*reaction 1: 2 CH4 + O2 => 2 CO + 4 H2 */
        phi_f = sc[0*npt+i]*sc[0*npt+i]*sc[1*npt+i];
        k_f = k_f_s[0*npt+i];
        q_f = phi_f * k_f;
        q_r = 0.0;
        qdot = q_f - q_r;
        wdot[0*npt+i] -= 2 * qdot;
        wdot[1*npt+i] -= qdot;
        wdot[4*npt+i] += 2 * qdot;
        wdot[6*npt+i] += 4 * qdot;

        /*reaction 2: CH4 + H2O <=> CO + 3 H2 */
        phi_f = sc[0*npt+i]*sc[2*npt+i];
        k_f = k_f_s[1*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[4*npt+i]*sc[6*npt+i]*sc[6*npt+i]*sc[6*npt+i];
        Kc = Kc_s[1*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[0*npt+i] -= qdot;
        wdot[2*npt+i] -= qdot;
        wdot[4*npt+i] += qdot;
        wdot[6*npt+i] += 3 * qdot;

        /*reaction 3: 2 H2 + O2 => 2 H2O */
        phi_f = sc[1*npt+i]*sc[6*npt+i]*sc[6*npt+i];
        k_f = k_f_s[2*npt+i];
        q_f = phi_f * k_f;
        q_r = 0.0;
        qdot = q_f - q_r;
        wdot[1*npt+i] -= qdot;
        wdot[2*npt+i] += 2 * qdot;
        wdot[6*npt+i] -= 2 * qdot;

        /*reaction 4: CO + H2O <=> CO2 + H2 */
        phi_f = sc[2*npt+i]*sc[4*npt+i];
        k_f = k_f_s[3*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[5*npt+i]*sc[6*npt+i];
        Kc = Kc_s[3*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[2*npt+i] -= qdot;
        wdot[4*npt+i] -= qdot;
        wdot[5*npt+i] += qdot;
        wdot[6*npt+i] += qdot;
    }
}

/*compute an approx to the reaction Jacobian */
void DWDOT_PRECOND(double *  J, double *  sc, double *  Tp, int * HP)
{
    double c[7];

    for (int k=0; k<7; k++) {
        c[k] = 1.e6 * sc[k];
    }

    aJacobian_precond(J, c, *Tp, *HP);

    /* dwdot[k]/dT */
    /* dTdot/d[X] */
    for (int k=0; k<7; k++) {
        J[56+k] *= 1.e-6;
        J[k*8+7] *= 1.e6;
    }

    return;
}

/*compute the reaction Jacobian */
void DWDOT(double *  J, double *  sc, double *  Tp, int * consP)
{
    double c[7];

    for (int k=0; k<7; k++) {
        c[k] = 1.e6 * sc[k];
    }

    aJacobian(J, c, *Tp, *consP);

    /* dwdot[k]/dT */
    /* dTdot/d[X] */
    for (int k=0; k<7; k++) {
        J[56+k] *= 1.e-6;
        J[k*8+7] *= 1.e6;
    }

    return;
}

/*compute the sparsity pattern Jacobian */
void SPARSITY_INFO( int * nJdata, int * consP, int NCELLS)
{
    double c[7];
    double J[64];

    for (int k=0; k<7; k++) {
        c[k] = 1.0/ 7.000000 ;
    }

    aJacobian(J, c, 1500.0, *consP);

    int nJdata_tmp = 0;
    for (int k=0; k<8; k++) {
        for (int l=0; l<8; l++) {
            if(J[ 8 * k + l] != 0.0){
                nJdata_tmp = nJdata_tmp + 1;
            }
        }
    }

    *nJdata = NCELLS * nJdata_tmp;

    return;
}



/*compute the sparsity pattern of simplified Jacobian */
void SPARSITY_INFO_PRECOND( int * nJdata, int * consP)
{
    double c[7];
    double J[64];

    for (int k=0; k<7; k++) {
        c[k] = 1.0/ 7.000000 ;
    }

    aJacobian_precond(J, c, 1500.0, *consP);

    int nJdata_tmp = 0;
    for (int k=0; k<8; k++) {
        for (int l=0; l<8; l++) {
            if(k == l){
                nJdata_tmp = nJdata_tmp + 1;
            } else {
                if(J[ 8 * k + l] != 0.0){
                    nJdata_tmp = nJdata_tmp + 1;
                }
            }
        }
    }

    nJdata[0] = nJdata_tmp;

    return;
}


/*compute the sparsity pattern of the simplified precond Jacobian */
void SPARSITY_PREPROC_PRECOND(int * rowVals, int * colPtrs, int * consP)
{
    double c[7];
    double J[64];

    for (int k=0; k<7; k++) {
        c[k] = 1.0/ 7.000000 ;
    }

    aJacobian_precond(J, c, 1500.0, *consP);

    colPtrs[0] = 0;
    int nJdata_tmp = 0;
    for (int k=0; k<8; k++) {
        for (int l=0; l<8; l++) {
            if (k == l) {
                rowVals[nJdata_tmp] = l; 
                nJdata_tmp = nJdata_tmp + 1; 
            } else {
                if(J[8*k + l] != 0.0) {
                    rowVals[nJdata_tmp] = l; 
                    nJdata_tmp = nJdata_tmp + 1; 
                }
            }
        }
        colPtrs[k+1] = nJdata_tmp;
    }

    return;
}
/*compute the sparsity pattern of the Jacobian */
void SPARSITY_PREPROC(int *  rowVals, int *  colPtrs, int * consP, int NCELLS)
{
    double c[7];
    double J[64];
    int offset_row;
    int offset_col;

    for (int k=0; k<7; k++) {
        c[k] = 1.0/ 7.000000 ;
    }

    aJacobian(J, c, 1500.0, *consP);

    colPtrs[0] = 0;
    int nJdata_tmp = 0;
    for (int nc=0; nc<NCELLS; nc++) {
        offset_row = nc * 8;
        offset_col = nc * 8;
        for (int k=0; k<8; k++) {
            for (int l=0; l<8; l++) {
                if(J[8*k + l] != 0.0) {
                    rowVals[nJdata_tmp] = l + offset_row; 
                    nJdata_tmp = nJdata_tmp + 1; 
                }
            }
            colPtrs[offset_col + (k + 1)] = nJdata_tmp;
        }
    }

    return;
}

/*compute the reaction Jacobian */
void aJacobian(double *  J, double *  sc, double T, int consP)
{
    for (int i=0; i<64; i++) {
        J[i] = 0.0;
    }

    double wdot[7];
    for (int k=0; k<7; k++) {
        wdot[k] = 0.0;
    }

    double tc[] = { log(T), T, T*T, T*T*T, T*T*T*T }; /*temperature cache */
    double invT = 1.0 / tc[1];
    double invT2 = invT * invT;

    /*reference concentration: P_atm / (RT) in inverse mol/m^3 */
    double refC = 101325 / 8.31451 / T;
    double refCinv = 1.0 / refC;

    /*compute the mixture concentration */
    double mixture = 0.0;
    for (int k = 0; k < 7; ++k) {
        mixture += sc[k];
    }

    /*compute the Gibbs free energy */
    double g_RT[7];
    gibbs(g_RT, tc);

    /*compute the species enthalpy */
    double h_RT[7];
    speciesEnthalpy(h_RT, tc);

    double phi_f, k_f, k_r, phi_r, Kc, q, q_nocor, Corr, alpha;
    double dlnkfdT, dlnk0dT, dlnKcdT, dkrdT, dqdT;
    double dqdci, dcdc_fac, dqdc[7];
    double Pr, fPr, F, k_0, logPr;
    double logFcent, troe_c, troe_n, troePr_den, troePr, troe;
    double Fcent1, Fcent2, Fcent3, Fcent;
    double dlogFdc, dlogFdn, dlogFdcn_fac;
    double dlogPrdT, dlogfPrdT, dlogFdT, dlogFcentdT, dlogFdlogPr, dlnCorrdT;
    const double ln10 = log(10.0);
    const double log10e = 1.0/log(10.0);
    /*reaction 1: 2 CH4 + O2 => 2 CO + 4 H2 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[0]*sc[0]*sc[1];
    k_f = prefactor_units[0] * fwd_A[0]
                * exp(fwd_beta[0] * tc[0] - activation_units[0] * fwd_Ea[0] * invT);
    dlnkfdT = fwd_beta[0] * invT + activation_units[0] * fwd_Ea[0] * invT2;
    /* rate of progress */
    q = k_f*phi_f;
    dqdT = dlnkfdT*k_f*phi_f;
    /* update wdot */
    wdot[0] -= 2 * q; /* CH4 */
    wdot[1] -= q; /* O2 */
    wdot[4] += 2 * q; /* CO */
    wdot[6] += 4 * q; /* H2 */
    /* d()/d[CH4] */
    dqdci =  + k_f*2*sc[0]*sc[1];
    J[0] += -2 * dqdci;           /* dwdot[CH4]/d[CH4] */
    J[1] -= dqdci;                /* dwdot[O2]/d[CH4] */
    J[4] += 2 * dqdci;            /* dwdot[CO]/d[CH4] */
    J[6] += 4 * dqdci;            /* dwdot[H2]/d[CH4] */
    /* d()/d[O2] */
    dqdci =  + k_f*sc[0]*sc[0];
    J[8] += -2 * dqdci;           /* dwdot[CH4]/d[O2] */
    J[9] -= dqdci;                /* dwdot[O2]/d[O2] */
    J[12] += 2 * dqdci;           /* dwdot[CO]/d[O2] */
    J[14] += 4 * dqdci;           /* dwdot[H2]/d[O2] */
    /* d()/dT */
    J[56] += -2 * dqdT;           /* dwdot[CH4]/dT */
    J[57] -= dqdT;                /* dwdot[O2]/dT */
    J[60] += 2 * dqdT;            /* dwdot[CO]/dT */
    J[62] += 4 * dqdT;            /* dwdot[H2]/dT */

    /*reaction 2: CH4 + H2O <=> CO + 3 H2 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[0]*sc[2];
    k_f = prefactor_units[1] * fwd_A[1]
                * exp(fwd_beta[1] * tc[0] - activation_units[1] * fwd_Ea[1] * invT);
    dlnkfdT = fwd_beta[1] * invT + activation_units[1] * fwd_Ea[1] * invT2;
    /* reverse */
    phi_r = sc[4]*sc[6]*sc[6]*sc[6];
    Kc = refC*refC * exp(g_RT[0] + g_RT[2] - g_RT[4] - 3*g_RT[6]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[0] + h_RT[2]) + (h_RT[4] + 3*h_RT[6]) - 2);
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[0] -= q; /* CH4 */
    wdot[2] -= q; /* H2O */
    wdot[4] += q; /* CO */
    wdot[6] += 3 * q; /* H2 */
    /* d()/d[CH4] */
    dqdci =  + k_f*sc[2];
    J[0] -= dqdci;                /* dwdot[CH4]/d[CH4] */
    J[2] -= dqdci;                /* dwdot[H2O]/d[CH4] */
    J[4] += dqdci;                /* dwdot[CO]/d[CH4] */
    J[6] += 3 * dqdci;            /* dwdot[H2]/d[CH4] */
    /* d()/d[H2O] */
    dqdci =  + k_f*sc[0];
    J[16] -= dqdci;               /* dwdot[CH4]/d[H2O] */
    J[18] -= dqdci;               /* dwdot[H2O]/d[H2O] */
    J[20] += dqdci;               /* dwdot[CO]/d[H2O] */
    J[22] += 3 * dqdci;           /* dwdot[H2]/d[H2O] */
    /* d()/d[CO] */
    dqdci =  - k_r*sc[6]*sc[6]*sc[6];
    J[32] -= dqdci;               /* dwdot[CH4]/d[CO] */
    J[34] -= dqdci;               /* dwdot[H2O]/d[CO] */
    J[36] += dqdci;               /* dwdot[CO]/d[CO] */
    J[38] += 3 * dqdci;           /* dwdot[H2]/d[CO] */
    /* d()/d[H2] */
    dqdci =  - k_r*sc[4]*3*sc[6]*sc[6];
    J[48] -= dqdci;               /* dwdot[CH4]/d[H2] */
    J[50] -= dqdci;               /* dwdot[H2O]/d[H2] */
    J[52] += dqdci;               /* dwdot[CO]/d[H2] */
    J[54] += 3 * dqdci;           /* dwdot[H2]/d[H2] */
    /* d()/dT */
    J[56] -= dqdT;                /* dwdot[CH4]/dT */
    J[58] -= dqdT;                /* dwdot[H2O]/dT */
    J[60] += dqdT;                /* dwdot[CO]/dT */
    J[62] += 3 * dqdT;            /* dwdot[H2]/dT */

    /*reaction 3: 2 H2 + O2 => 2 H2O */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[1]*sc[6]*sc[6];
    k_f = prefactor_units[2] * fwd_A[2]
                * exp(fwd_beta[2] * tc[0] - activation_units[2] * fwd_Ea[2] * invT);
    dlnkfdT = fwd_beta[2] * invT + activation_units[2] * fwd_Ea[2] * invT2;
    /* rate of progress */
    q = k_f*phi_f;
    dqdT = dlnkfdT*k_f*phi_f;
    /* update wdot */
    wdot[1] -= q; /* O2 */
    wdot[2] += 2 * q; /* H2O */
    wdot[6] -= 2 * q; /* H2 */
    /* d()/d[O2] */
    dqdci =  + k_f*sc[6]*sc[6];
    J[9] -= dqdci;                /* dwdot[O2]/d[O2] */
    J[10] += 2 * dqdci;           /* dwdot[H2O]/d[O2] */
    J[14] += -2 * dqdci;          /* dwdot[H2]/d[O2] */
    /* d()/d[H2] */
    dqdci =  + k_f*sc[1]*2*sc[6];
    J[49] -= dqdci;               /* dwdot[O2]/d[H2] */
    J[50] += 2 * dqdci;           /* dwdot[H2O]/d[H2] */
    J[54] += -2 * dqdci;          /* dwdot[H2]/d[H2] */
    /* d()/dT */
    J[57] -= dqdT;                /* dwdot[O2]/dT */
    J[58] += 2 * dqdT;            /* dwdot[H2O]/dT */
    J[62] += -2 * dqdT;           /* dwdot[H2]/dT */

    /*reaction 4: CO + H2O <=> CO2 + H2 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[2]*sc[4];
    k_f = prefactor_units[3] * fwd_A[3]
                * exp(fwd_beta[3] * tc[0] - activation_units[3] * fwd_Ea[3] * invT);
    dlnkfdT = fwd_beta[3] * invT + activation_units[3] * fwd_Ea[3] * invT2;
    /* reverse */
    phi_r = sc[5]*sc[6];
    Kc = exp(g_RT[2] + g_RT[4] - g_RT[5] - g_RT[6]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[2] + h_RT[4]) + (h_RT[5] + h_RT[6]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[2] -= q; /* H2O */
    wdot[4] -= q; /* CO */
    wdot[5] += q; /* CO2 */
    wdot[6] += q; /* H2 */
    /* d()/d[H2O] */
    dqdci =  + k_f*sc[4];
    J[18] -= dqdci;               /* dwdot[H2O]/d[H2O] */
    J[20] -= dqdci;               /* dwdot[CO]/d[H2O] */
    J[21] += dqdci;               /* dwdot[CO2]/d[H2O] */
    J[22] += dqdci;               /* dwdot[H2]/d[H2O] */
    /* d()/d[CO] */
    dqdci =  + k_f*sc[2];
    J[34] -= dqdci;               /* dwdot[H2O]/d[CO] */
    J[36] -= dqdci;               /* dwdot[CO]/d[CO] */
    J[37] += dqdci;               /* dwdot[CO2]/d[CO] */
    J[38] += dqdci;               /* dwdot[H2]/d[CO] */
    /* d()/d[CO2] */
    dqdci =  - k_r*sc[6];
    J[42] -= dqdci;               /* dwdot[H2O]/d[CO2] */
    J[44] -= dqdci;               /* dwdot[CO]/d[CO2] */
    J[45] += dqdci;               /* dwdot[CO2]/d[CO2] */
    J[46] += dqdci;               /* dwdot[H2]/d[CO2] */
    /* d()/d[H2] */
    dqdci =  - k_r*sc[5];
    J[50] -= dqdci;               /* dwdot[H2O]/d[H2] */
    J[52] -= dqdci;               /* dwdot[CO]/d[H2] */
    J[53] += dqdci;               /* dwdot[CO2]/d[H2] */
    J[54] += dqdci;               /* dwdot[H2]/d[H2] */
    /* d()/dT */
    J[58] -= dqdT;                /* dwdot[H2O]/dT */
    J[60] -= dqdT;                /* dwdot[CO]/dT */
    J[61] += dqdT;                /* dwdot[CO2]/dT */
    J[62] += dqdT;                /* dwdot[H2]/dT */

    double c_R[7], dcRdT[7], e_RT[7];
    double * eh_RT;
    if (consP) {
        cp_R(c_R, tc);
        dcvpRdT(dcRdT, tc);
        eh_RT = &h_RT[0];
    }
    else {
        cv_R(c_R, tc);
        dcvpRdT(dcRdT, tc);
        speciesInternalEnergy(e_RT, tc);
        eh_RT = &e_RT[0];
    }

    double cmix = 0.0, ehmix = 0.0, dcmixdT=0.0, dehmixdT=0.0;
    for (int k = 0; k < 7; ++k) {
        cmix += c_R[k]*sc[k];
        dcmixdT += dcRdT[k]*sc[k];
        ehmix += eh_RT[k]*wdot[k];
        dehmixdT += invT*(c_R[k]-eh_RT[k])*wdot[k] + eh_RT[k]*J[56+k];
    }

    double cmixinv = 1.0/cmix;
    double tmp1 = ehmix*cmixinv;
    double tmp3 = cmixinv*T;
    double tmp2 = tmp1*tmp3;
    double dehmixdc;
    /* dTdot/d[X] */
    for (int k = 0; k < 7; ++k) {
        dehmixdc = 0.0;
        for (int m = 0; m < 7; ++m) {
            dehmixdc += eh_RT[m]*J[k*8+m];
        }
        J[k*8+7] = tmp2*c_R[k] - tmp3*dehmixdc;
    }
    /* dTdot/dT */
    J[63] = -tmp1 + tmp2*dcmixdT - tmp3*dehmixdT;
}

/*compute an approx to the reaction Jacobian */
void aJacobian_precond(double *  J, double *  sc, double T, int HP)
{
    for (int i=0; i<64; i++) {
        J[i] = 0.0;
    }

    double wdot[7];
    for (int k=0; k<7; k++) {
        wdot[k] = 0.0;
    }

    double tc[] = { log(T), T, T*T, T*T*T, T*T*T*T }; /*temperature cache */
    double invT = 1.0 / tc[1];
    double invT2 = invT * invT;

    /*reference concentration: P_atm / (RT) in inverse mol/m^3 */
    double refC = 101325 / 8.31451 / T;
    double refCinv = 1.0 / refC;

    /*compute the mixture concentration */
    double mixture = 0.0;
    for (int k = 0; k < 7; ++k) {
        mixture += sc[k];
    }

    /*compute the Gibbs free energy */
    double g_RT[7];
    gibbs(g_RT, tc);

    /*compute the species enthalpy */
    double h_RT[7];
    speciesEnthalpy(h_RT, tc);

    double phi_f, k_f, k_r, phi_r, Kc, q, q_nocor, Corr, alpha;
    double dlnkfdT, dlnk0dT, dlnKcdT, dkrdT, dqdT;
    double dqdci, dcdc_fac, dqdc[7];
    double Pr, fPr, F, k_0, logPr;
    double logFcent, troe_c, troe_n, troePr_den, troePr, troe;
    double Fcent1, Fcent2, Fcent3, Fcent;
    double dlogFdc, dlogFdn, dlogFdcn_fac;
    double dlogPrdT, dlogfPrdT, dlogFdT, dlogFcentdT, dlogFdlogPr, dlnCorrdT;
    const double ln10 = log(10.0);
    const double log10e = 1.0/log(10.0);
    /*reaction 1: 2 CH4 + O2 => 2 CO + 4 H2 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[0]*sc[0]*sc[1];
    k_f = prefactor_units[0] * fwd_A[0]
                * exp(fwd_beta[0] * tc[0] - activation_units[0] * fwd_Ea[0] * invT);
    dlnkfdT = fwd_beta[0] * invT + activation_units[0] * fwd_Ea[0] * invT2;
    /* rate of progress */
    q = k_f*phi_f;
    dqdT = dlnkfdT*k_f*phi_f;
    /* update wdot */
    wdot[0] -= 2 * q; /* CH4 */
    wdot[1] -= q; /* O2 */
    wdot[4] += 2 * q; /* CO */
    wdot[6] += 4 * q; /* H2 */
    /* d()/d[CH4] */
    dqdci =  + k_f*2*sc[0]*sc[1];
    J[0] += -2 * dqdci;           /* dwdot[CH4]/d[CH4] */
    J[1] -= dqdci;                /* dwdot[O2]/d[CH4] */
    J[4] += 2 * dqdci;            /* dwdot[CO]/d[CH4] */
    J[6] += 4 * dqdci;            /* dwdot[H2]/d[CH4] */
    /* d()/d[O2] */
    dqdci =  + k_f*sc[0]*sc[0];
    J[8] += -2 * dqdci;           /* dwdot[CH4]/d[O2] */
    J[9] -= dqdci;                /* dwdot[O2]/d[O2] */
    J[12] += 2 * dqdci;           /* dwdot[CO]/d[O2] */
    J[14] += 4 * dqdci;           /* dwdot[H2]/d[O2] */
    /* d()/dT */
    J[56] += -2 * dqdT;           /* dwdot[CH4]/dT */
    J[57] -= dqdT;                /* dwdot[O2]/dT */
    J[60] += 2 * dqdT;            /* dwdot[CO]/dT */
    J[62] += 4 * dqdT;            /* dwdot[H2]/dT */

    /*reaction 2: CH4 + H2O <=> CO + 3 H2 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[0]*sc[2];
    k_f = prefactor_units[1] * fwd_A[1]
                * exp(fwd_beta[1] * tc[0] - activation_units[1] * fwd_Ea[1] * invT);
    dlnkfdT = fwd_beta[1] * invT + activation_units[1] * fwd_Ea[1] * invT2;
    /* reverse */
    phi_r = sc[4]*sc[6]*sc[6]*sc[6];
    Kc = refC*refC * exp(g_RT[0] + g_RT[2] - g_RT[4] - 3*g_RT[6]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[0] + h_RT[2]) + (h_RT[4] + 3*h_RT[6]) - 2);
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[0] -= q; /* CH4 */
    wdot[2] -= q; /* H2O */
    wdot[4] += q; /* CO */
    wdot[6] += 3 * q; /* H2 */
    /* d()/d[CH4] */
    dqdci =  + k_f*sc[2];
    J[0] -= dqdci;                /* dwdot[CH4]/d[CH4] */
    J[2] -= dqdci;                /* dwdot[H2O]/d[CH4] */
    J[4] += dqdci;                /* dwdot[CO]/d[CH4] */
    J[6] += 3 * dqdci;            /* dwdot[H2]/d[CH4] */
    /* d()/d[H2O] */
    dqdci =  + k_f*sc[0];
    J[16] -= dqdci;               /* dwdot[CH4]/d[H2O] */
    J[18] -= dqdci;               /* dwdot[H2O]/d[H2O] */
    J[20] += dqdci;               /* dwdot[CO]/d[H2O] */
    J[22] += 3 * dqdci;           /* dwdot[H2]/d[H2O] */
    /* d()/d[CO] */
    dqdci =  - k_r*sc[6]*sc[6]*sc[6];
    J[32] -= dqdci;               /* dwdot[CH4]/d[CO] */
    J[34] -= dqdci;               /* dwdot[H2O]/d[CO] */
    J[36] += dqdci;               /* dwdot[CO]/d[CO] */
    J[38] += 3 * dqdci;           /* dwdot[H2]/d[CO] */
    /* d()/d[H2] */
    dqdci =  - k_r*sc[4]*3*sc[6]*sc[6];
    J[48] -= dqdci;               /* dwdot[CH4]/d[H2] */
    J[50] -= dqdci;               /* dwdot[H2O]/d[H2] */
    J[52] += dqdci;               /* dwdot[CO]/d[H2] */
    J[54] += 3 * dqdci;           /* dwdot[H2]/d[H2] */
    /* d()/dT */
    J[56] -= dqdT;                /* dwdot[CH4]/dT */
    J[58] -= dqdT;                /* dwdot[H2O]/dT */
    J[60] += dqdT;                /* dwdot[CO]/dT */
    J[62] += 3 * dqdT;            /* dwdot[H2]/dT */

    /*reaction 3: 2 H2 + O2 => 2 H2O */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[1]*sc[6]*sc[6];
    k_f = prefactor_units[2] * fwd_A[2]
                * exp(fwd_beta[2] * tc[0] - activation_units[2] * fwd_Ea[2] * invT);
    dlnkfdT = fwd_beta[2] * invT + activation_units[2] * fwd_Ea[2] * invT2;
    /* rate of progress */
    q = k_f*phi_f;
    dqdT = dlnkfdT*k_f*phi_f;
    /* update wdot */
    wdot[1] -= q; /* O2 */
    wdot[2] += 2 * q; /* H2O */
    wdot[6] -= 2 * q; /* H2 */
    /* d()/d[O2] */
    dqdci =  + k_f*sc[6]*sc[6];
    J[9] -= dqdci;                /* dwdot[O2]/d[O2] */
    J[10] += 2 * dqdci;           /* dwdot[H2O]/d[O2] */
    J[14] += -2 * dqdci;          /* dwdot[H2]/d[O2] */
    /* d()/d[H2] */
    dqdci =  + k_f*sc[1]*2*sc[6];
    J[49] -= dqdci;               /* dwdot[O2]/d[H2] */
    J[50] += 2 * dqdci;           /* dwdot[H2O]/d[H2] */
    J[54] += -2 * dqdci;          /* dwdot[H2]/d[H2] */
    /* d()/dT */
    J[57] -= dqdT;                /* dwdot[O2]/dT */
    J[58] += 2 * dqdT;            /* dwdot[H2O]/dT */
    J[62] += -2 * dqdT;           /* dwdot[H2]/dT */

    /*reaction 4: CO + H2O <=> CO2 + H2 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[2]*sc[4];
    k_f = prefactor_units[3] * fwd_A[3]
                * exp(fwd_beta[3] * tc[0] - activation_units[3] * fwd_Ea[3] * invT);
    dlnkfdT = fwd_beta[3] * invT + activation_units[3] * fwd_Ea[3] * invT2;
    /* reverse */
    phi_r = sc[5]*sc[6];
    Kc = exp(g_RT[2] + g_RT[4] - g_RT[5] - g_RT[6]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[2] + h_RT[4]) + (h_RT[5] + h_RT[6]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[2] -= q; /* H2O */
    wdot[4] -= q; /* CO */
    wdot[5] += q; /* CO2 */
    wdot[6] += q; /* H2 */
    /* d()/d[H2O] */
    dqdci =  + k_f*sc[4];
    J[18] -= dqdci;               /* dwdot[H2O]/d[H2O] */
    J[20] -= dqdci;               /* dwdot[CO]/d[H2O] */
    J[21] += dqdci;               /* dwdot[CO2]/d[H2O] */
    J[22] += dqdci;               /* dwdot[H2]/d[H2O] */
    /* d()/d[CO] */
    dqdci =  + k_f*sc[2];
    J[34] -= dqdci;               /* dwdot[H2O]/d[CO] */
    J[36] -= dqdci;               /* dwdot[CO]/d[CO] */
    J[37] += dqdci;               /* dwdot[CO2]/d[CO] */
    J[38] += dqdci;               /* dwdot[H2]/d[CO] */
    /* d()/d[CO2] */
    dqdci =  - k_r*sc[6];
    J[42] -= dqdci;               /* dwdot[H2O]/d[CO2] */
    J[44] -= dqdci;               /* dwdot[CO]/d[CO2] */
    J[45] += dqdci;               /* dwdot[CO2]/d[CO2] */
    J[46] += dqdci;               /* dwdot[H2]/d[CO2] */
    /* d()/d[H2] */
    dqdci =  - k_r*sc[5];
    J[50] -= dqdci;               /* dwdot[H2O]/d[H2] */
    J[52] -= dqdci;               /* dwdot[CO]/d[H2] */
    J[53] += dqdci;               /* dwdot[CO2]/d[H2] */
    J[54] += dqdci;               /* dwdot[H2]/d[H2] */
    /* d()/dT */
    J[58] -= dqdT;                /* dwdot[H2O]/dT */
    J[60] -= dqdT;                /* dwdot[CO]/dT */
    J[61] += dqdT;                /* dwdot[CO2]/dT */
    J[62] += dqdT;                /* dwdot[H2]/dT */

    double c_R[7], dcRdT[7], e_RT[7];
    double * eh_RT;
    if (HP) {
        cp_R(c_R, tc);
        dcvpRdT(dcRdT, tc);
        eh_RT = &h_RT[0];
    }
    else {
        cv_R(c_R, tc);
        dcvpRdT(dcRdT, tc);
        speciesInternalEnergy(e_RT, tc);
        eh_RT = &e_RT[0];
    }

    double cmix = 0.0, ehmix = 0.0, dcmixdT=0.0, dehmixdT=0.0;
    for (int k = 0; k < 7; ++k) {
        cmix += c_R[k]*sc[k];
        dcmixdT += dcRdT[k]*sc[k];
        ehmix += eh_RT[k]*wdot[k];
        dehmixdT += invT*(c_R[k]-eh_RT[k])*wdot[k] + eh_RT[k]*J[56+k];
    }

    double cmixinv = 1.0/cmix;
    double tmp1 = ehmix*cmixinv;
    double tmp3 = cmixinv*T;
    double tmp2 = tmp1*tmp3;
    double dehmixdc;
    /* dTdot/d[X] */
    for (int k = 0; k < 7; ++k) {
        dehmixdc = 0.0;
        for (int m = 0; m < 7; ++m) {
            dehmixdc += eh_RT[m]*J[k*8+m];
        }
        J[k*8+7] = tmp2*c_R[k] - tmp3*dehmixdc;
    }
    /* dTdot/dT */
    J[63] = -tmp1 + tmp2*dcmixdT - tmp3*dehmixdT;
}


/*compute d(Cp/R)/dT and d(Cv/R)/dT at the given temperature */
/*tc contains precomputed powers of T, tc[0] = log(T) */
void dcvpRdT(double *  species, double *  tc)
{

    /*temperature */
    double T = tc[1];

    /*species with midpoint at T=1000 kelvin */
    if (T < 1000) {
        /*species 0: CH4 */
        species[0] =
            -1.36709788e-02
            +9.83601198e-05 * tc[1]
            -1.45422908e-07 * tc[2]
            +6.66775824e-11 * tc[3];
        /*species 1: O2 */
        species[1] =
            -2.99673416e-03
            +1.96946040e-05 * tc[1]
            -2.90438853e-08 * tc[2]
            +1.29749135e-11 * tc[3];
        /*species 2: H2O */
        species[2] =
            -2.03643410e-03
            +1.30408042e-05 * tc[1]
            -1.64639119e-08 * tc[2]
            +7.08791268e-12 * tc[3];
        /*species 3: N2 */
        species[3] =
            +1.40824040e-03
            -7.92644400e-06 * tc[1]
            +1.69245450e-08 * tc[2]
            -9.77941600e-12 * tc[3];
        /*species 4: CO */
        species[4] =
            -6.10353680e-04
            +2.03362866e-06 * tc[1]
            +2.72101765e-09 * tc[2]
            -3.61769800e-12 * tc[3];
        /*species 5: CO2 */
        species[5] =
            +8.98459677e-03
            -1.42471254e-05 * tc[1]
            +7.37757066e-09 * tc[2]
            -5.74798192e-13 * tc[3];
        /*species 6: H2 */
        species[6] =
            +7.98052075e-03
            -3.89563020e-05 * tc[1]
            +6.04716282e-08 * tc[2]
            -2.95044704e-11 * tc[3];
    } else {
        /*species 0: CH4 */
        species[0] =
            +1.33909467e-02
            -1.14657162e-05 * tc[1]
            +3.66877605e-09 * tc[2]
            -4.07260920e-13 * tc[3];
        /*species 1: O2 */
        species[1] =
            +1.48308754e-03
            -1.51593334e-06 * tc[1]
            +6.28411665e-10 * tc[2]
            -8.66871176e-14 * tc[3];
        /*species 2: H2O */
        species[2] =
            +2.17691804e-03
            -3.28145036e-07 * tc[1]
            -2.91125961e-10 * tc[2]
            +6.72803968e-14 * tc[3];
        /*species 3: N2 */
        species[3] =
            +1.48797680e-03
            -1.13695200e-06 * tc[1]
            +3.02911140e-10 * tc[2]
            -2.70134040e-14 * tc[3];
        /*species 4: CO */
        species[4] =
            +2.06252743e-03
            -1.99765154e-06 * tc[1]
            +6.90159024e-10 * tc[2]
            -8.14590864e-14 * tc[3];
        /*species 5: CO2 */
        species[5] =
            +4.41437026e-03
            -4.42962808e-06 * tc[1]
            +1.57047056e-09 * tc[2]
            -1.88833666e-13 * tc[3];
        /*species 6: H2 */
        species[6] =
            -4.94024731e-05
            +9.98913556e-07 * tc[1]
            -5.38699182e-10 * tc[2]
            +8.01021504e-14 * tc[3];
    }
    return;
}


/*compute the progress rate for each reaction */
void progressRate(double *  qdot, double *  sc, double T)
{
    double tc[] = { log(T), T, T*T, T*T*T, T*T*T*T }; /*temperature cache */
    double invT = 1.0 / tc[1];

    if (T != T_save)
    {
        T_save = T;
        comp_k_f(tc,invT,k_f_save);
        comp_Kc(tc,invT,Kc_save);
    }

    double q_f[4], q_r[4];
    comp_qfqr(q_f, q_r, sc, tc, invT);

    for (int i = 0; i < 4; ++i) {
        qdot[i] = q_f[i] - q_r[i];
    }

    return;
}


/*compute the progress rate for each reaction */
void progressRateFR(double *  q_f, double *  q_r, double *  sc, double T)
{
    double tc[] = { log(T), T, T*T, T*T*T, T*T*T*T }; /*temperature cache */
    double invT = 1.0 / tc[1];

    if (T != T_save)
    {
        T_save = T;
        comp_k_f(tc,invT,k_f_save);
        comp_Kc(tc,invT,Kc_save);
    }

    comp_qfqr(q_f, q_r, sc, tc, invT);

    return;
}


/*compute the equilibrium constants for each reaction */
void equilibriumConstants(double *  kc, double *  g_RT, double T)
{
    /*reference concentration: P_atm / (RT) in inverse mol/m^3 */
    double refC = 101325 / 8.31451 / T;

    /*reaction 1: 2 CH4 + O2 => 2 CO + 4 H2 */
    kc[0] = refC*refC*refC * exp((2 * g_RT[0] + g_RT[1]) - (2 * g_RT[4] + 4 * g_RT[6]));

    /*reaction 2: CH4 + H2O <=> CO + 3 H2 */
    kc[1] = refC*refC * exp((g_RT[0] + g_RT[2]) - (g_RT[4] + 3 * g_RT[6]));

    /*reaction 3: 2 H2 + O2 => 2 H2O */
    kc[2] = 1.0 / (refC) * exp((2 * g_RT[6] + g_RT[1]) - (2 * g_RT[2]));

    /*reaction 4: CO + H2O <=> CO2 + H2 */
    kc[3] = exp((g_RT[4] + g_RT[2]) - (g_RT[5] + g_RT[6]));

    return;
}


/*compute the g/(RT) at the given temperature */
/*tc contains precomputed powers of T, tc[0] = log(T) */
void gibbs(double *  species, double *  tc)
{

    /*temperature */
    double T = tc[1];
    double invT = 1 / T;

    /*species with midpoint at T=1000 kelvin */
    if (T < 1000) {
        /*species 0: CH4 */
        species[0] =
            -1.024664760000000e+04 * invT
            +9.791179889999999e+00
            -5.149876130000000e+00 * tc[0]
            +6.835489400000000e-03 * tc[1]
            -8.196676650000000e-06 * tc[2]
            +4.039525216666667e-09 * tc[3]
            -8.334697800000000e-13 * tc[4];
        /*species 1: O2 */
        species[1] =
            -1.063943560000000e+03 * invT
            +1.247806300000001e-01
            -3.782456360000000e+00 * tc[0]
            +1.498367080000000e-03 * tc[1]
            -1.641217001666667e-06 * tc[2]
            +8.067745908333334e-10 * tc[3]
            -1.621864185000000e-13 * tc[4];
        /*species 2: H2O */
        species[2] =
            -3.029372670000000e+04 * invT
            +5.047672768000000e+00
            -4.198640560000000e+00 * tc[0]
            +1.018217050000000e-03 * tc[1]
            -1.086733685000000e-06 * tc[2]
            +4.573308850000000e-10 * tc[3]
            -8.859890850000000e-14 * tc[4];
        /*species 3: N2 */
        species[3] =
            -1.020899900000000e+03 * invT
            -6.516950000000001e-01
            -3.298677000000000e+00 * tc[0]
            -7.041202000000000e-04 * tc[1]
            +6.605369999999999e-07 * tc[2]
            -4.701262500000001e-10 * tc[3]
            +1.222427000000000e-13 * tc[4];
        /*species 4: CO */
        species[4] =
            -1.434408600000000e+04 * invT
            +7.112418999999992e-02
            -3.579533470000000e+00 * tc[0]
            +3.051768400000000e-04 * tc[1]
            -1.694690550000000e-07 * tc[2]
            -7.558382366666667e-11 * tc[3]
            +4.522122495000000e-14 * tc[4];
        /*species 5: CO2 */
        species[5] =
            -4.837196970000000e+04 * invT
            -7.544278700000000e+00
            -2.356773520000000e+00 * tc[0]
            -4.492298385000000e-03 * tc[1]
            +1.187260448333333e-06 * tc[2]
            -2.049325183333333e-10 * tc[3]
            +7.184977399999999e-15 * tc[4];
        /*species 6: H2 */
        species[6] =
            -9.179351730000000e+02 * invT
            +1.661320882000000e+00
            -2.344331120000000e+00 * tc[0]
            -3.990260375000000e-03 * tc[1]
            +3.246358500000000e-06 * tc[2]
            -1.679767450000000e-09 * tc[3]
            +3.688058805000000e-13 * tc[4];
    } else {
        /*species 0: CH4 */
        species[0] =
            -9.468344590000001e+03 * invT
            -1.836246650500000e+01
            -7.485149500000000e-02 * tc[0]
            -6.695473350000000e-03 * tc[1]
            +9.554763483333333e-07 * tc[2]
            -1.019104458333333e-10 * tc[3]
            +5.090761500000000e-15 * tc[4];
        /*species 1: O2 */
        species[1] =
            -1.088457720000000e+03 * invT
            -2.170693450000000e+00
            -3.282537840000000e+00 * tc[0]
            -7.415437700000000e-04 * tc[1]
            +1.263277781666667e-07 * tc[2]
            -1.745587958333333e-11 * tc[3]
            +1.083588970000000e-15 * tc[4];
        /*species 2: H2O */
        species[2] =
            -3.000429710000000e+04 * invT
            -1.932777610000000e+00
            -3.033992490000000e+00 * tc[0]
            -1.088459020000000e-03 * tc[1]
            +2.734541966666666e-08 * tc[2]
            +8.086832250000000e-12 * tc[3]
            -8.410049600000000e-16 * tc[4];
        /*species 3: N2 */
        species[3] =
            -9.227977000000000e+02 * invT
            -3.053888000000000e+00
            -2.926640000000000e+00 * tc[0]
            -7.439884000000000e-04 * tc[1]
            +9.474600000000001e-08 * tc[2]
            -8.414198333333333e-12 * tc[3]
            +3.376675500000000e-16 * tc[4];
        /*species 4: CO */
        species[4] =
            -1.415187240000000e+04 * invT
            -5.103502110000000e+00
            -2.715185610000000e+00 * tc[0]
            -1.031263715000000e-03 * tc[1]
            +1.664709618333334e-07 * tc[2]
            -1.917108400000000e-11 * tc[3]
            +1.018238580000000e-15 * tc[4];
        /*species 5: CO2 */
        species[5] =
            -4.875916600000000e+04 * invT
            +1.585822230000000e+00
            -3.857460290000000e+00 * tc[0]
            -2.207185130000000e-03 * tc[1]
            +3.691356733333334e-07 * tc[2]
            -4.362418233333334e-11 * tc[3]
            +2.360420820000000e-15 * tc[4];
        /*species 6: H2 */
        species[6] =
            -9.501589220000000e+02 * invT
            +6.542302510000000e+00
            -3.337279200000000e+00 * tc[0]
            +2.470123655000000e-05 * tc[1]
            -8.324279633333333e-08 * tc[2]
            +1.496386616666667e-11 * tc[3]
            -1.001276880000000e-15 * tc[4];
    }
    return;
}


/*compute the a/(RT) at the given temperature */
/*tc contains precomputed powers of T, tc[0] = log(T) */
void helmholtz(double *  species, double *  tc)
{

    /*temperature */
    double T = tc[1];
    double invT = 1 / T;

    /*species with midpoint at T=1000 kelvin */
    if (T < 1000) {
        /*species 0: CH4 */
        species[0] =
            -1.02466476e+04 * invT
            +8.79117989e+00
            -5.14987613e+00 * tc[0]
            +6.83548940e-03 * tc[1]
            -8.19667665e-06 * tc[2]
            +4.03952522e-09 * tc[3]
            -8.33469780e-13 * tc[4];
        /*species 1: O2 */
        species[1] =
            -1.06394356e+03 * invT
            -8.75219370e-01
            -3.78245636e+00 * tc[0]
            +1.49836708e-03 * tc[1]
            -1.64121700e-06 * tc[2]
            +8.06774591e-10 * tc[3]
            -1.62186418e-13 * tc[4];
        /*species 2: H2O */
        species[2] =
            -3.02937267e+04 * invT
            +4.04767277e+00
            -4.19864056e+00 * tc[0]
            +1.01821705e-03 * tc[1]
            -1.08673369e-06 * tc[2]
            +4.57330885e-10 * tc[3]
            -8.85989085e-14 * tc[4];
        /*species 3: N2 */
        species[3] =
            -1.02089990e+03 * invT
            -1.65169500e+00
            -3.29867700e+00 * tc[0]
            -7.04120200e-04 * tc[1]
            +6.60537000e-07 * tc[2]
            -4.70126250e-10 * tc[3]
            +1.22242700e-13 * tc[4];
        /*species 4: CO */
        species[4] =
            -1.43440860e+04 * invT
            -9.28875810e-01
            -3.57953347e+00 * tc[0]
            +3.05176840e-04 * tc[1]
            -1.69469055e-07 * tc[2]
            -7.55838237e-11 * tc[3]
            +4.52212249e-14 * tc[4];
        /*species 5: CO2 */
        species[5] =
            -4.83719697e+04 * invT
            -8.54427870e+00
            -2.35677352e+00 * tc[0]
            -4.49229839e-03 * tc[1]
            +1.18726045e-06 * tc[2]
            -2.04932518e-10 * tc[3]
            +7.18497740e-15 * tc[4];
        /*species 6: H2 */
        species[6] =
            -9.17935173e+02 * invT
            +6.61320882e-01
            -2.34433112e+00 * tc[0]
            -3.99026037e-03 * tc[1]
            +3.24635850e-06 * tc[2]
            -1.67976745e-09 * tc[3]
            +3.68805881e-13 * tc[4];
    } else {
        /*species 0: CH4 */
        species[0] =
            -9.46834459e+03 * invT
            -1.93624665e+01
            -7.48514950e-02 * tc[0]
            -6.69547335e-03 * tc[1]
            +9.55476348e-07 * tc[2]
            -1.01910446e-10 * tc[3]
            +5.09076150e-15 * tc[4];
        /*species 1: O2 */
        species[1] =
            -1.08845772e+03 * invT
            -3.17069345e+00
            -3.28253784e+00 * tc[0]
            -7.41543770e-04 * tc[1]
            +1.26327778e-07 * tc[2]
            -1.74558796e-11 * tc[3]
            +1.08358897e-15 * tc[4];
        /*species 2: H2O */
        species[2] =
            -3.00042971e+04 * invT
            -2.93277761e+00
            -3.03399249e+00 * tc[0]
            -1.08845902e-03 * tc[1]
            +2.73454197e-08 * tc[2]
            +8.08683225e-12 * tc[3]
            -8.41004960e-16 * tc[4];
        /*species 3: N2 */
        species[3] =
            -9.22797700e+02 * invT
            -4.05388800e+00
            -2.92664000e+00 * tc[0]
            -7.43988400e-04 * tc[1]
            +9.47460000e-08 * tc[2]
            -8.41419833e-12 * tc[3]
            +3.37667550e-16 * tc[4];
        /*species 4: CO */
        species[4] =
            -1.41518724e+04 * invT
            -6.10350211e+00
            -2.71518561e+00 * tc[0]
            -1.03126372e-03 * tc[1]
            +1.66470962e-07 * tc[2]
            -1.91710840e-11 * tc[3]
            +1.01823858e-15 * tc[4];
        /*species 5: CO2 */
        species[5] =
            -4.87591660e+04 * invT
            +5.85822230e-01
            -3.85746029e+00 * tc[0]
            -2.20718513e-03 * tc[1]
            +3.69135673e-07 * tc[2]
            -4.36241823e-11 * tc[3]
            +2.36042082e-15 * tc[4];
        /*species 6: H2 */
        species[6] =
            -9.50158922e+02 * invT
            +5.54230251e+00
            -3.33727920e+00 * tc[0]
            +2.47012365e-05 * tc[1]
            -8.32427963e-08 * tc[2]
            +1.49638662e-11 * tc[3]
            -1.00127688e-15 * tc[4];
    }
    return;
}


/*compute Cv/R at the given temperature */
/*tc contains precomputed powers of T, tc[0] = log(T) */
void cv_R(double *  species, double *  tc)
{

    /*temperature */
    double T = tc[1];

    /*species with midpoint at T=1000 kelvin */
    if (T < 1000) {
        /*species 0: CH4 */
        species[0] =
            +4.14987613e+00
            -1.36709788e-02 * tc[1]
            +4.91800599e-05 * tc[2]
            -4.84743026e-08 * tc[3]
            +1.66693956e-11 * tc[4];
        /*species 1: O2 */
        species[1] =
            +2.78245636e+00
            -2.99673416e-03 * tc[1]
            +9.84730201e-06 * tc[2]
            -9.68129509e-09 * tc[3]
            +3.24372837e-12 * tc[4];
        /*species 2: H2O */
        species[2] =
            +3.19864056e+00
            -2.03643410e-03 * tc[1]
            +6.52040211e-06 * tc[2]
            -5.48797062e-09 * tc[3]
            +1.77197817e-12 * tc[4];
        /*species 3: N2 */
        species[3] =
            +2.29867700e+00
            +1.40824040e-03 * tc[1]
            -3.96322200e-06 * tc[2]
            +5.64151500e-09 * tc[3]
            -2.44485400e-12 * tc[4];
        /*species 4: CO */
        species[4] =
            +2.57953347e+00
            -6.10353680e-04 * tc[1]
            +1.01681433e-06 * tc[2]
            +9.07005884e-10 * tc[3]
            -9.04424499e-13 * tc[4];
        /*species 5: CO2 */
        species[5] =
            +1.35677352e+00
            +8.98459677e-03 * tc[1]
            -7.12356269e-06 * tc[2]
            +2.45919022e-09 * tc[3]
            -1.43699548e-13 * tc[4];
        /*species 6: H2 */
        species[6] =
            +1.34433112e+00
            +7.98052075e-03 * tc[1]
            -1.94781510e-05 * tc[2]
            +2.01572094e-08 * tc[3]
            -7.37611761e-12 * tc[4];
    } else {
        /*species 0: CH4 */
        species[0] =
            -9.25148505e-01
            +1.33909467e-02 * tc[1]
            -5.73285809e-06 * tc[2]
            +1.22292535e-09 * tc[3]
            -1.01815230e-13 * tc[4];
        /*species 1: O2 */
        species[1] =
            +2.28253784e+00
            +1.48308754e-03 * tc[1]
            -7.57966669e-07 * tc[2]
            +2.09470555e-10 * tc[3]
            -2.16717794e-14 * tc[4];
        /*species 2: H2O */
        species[2] =
            +2.03399249e+00
            +2.17691804e-03 * tc[1]
            -1.64072518e-07 * tc[2]
            -9.70419870e-11 * tc[3]
            +1.68200992e-14 * tc[4];
        /*species 3: N2 */
        species[3] =
            +1.92664000e+00
            +1.48797680e-03 * tc[1]
            -5.68476000e-07 * tc[2]
            +1.00970380e-10 * tc[3]
            -6.75335100e-15 * tc[4];
        /*species 4: CO */
        species[4] =
            +1.71518561e+00
            +2.06252743e-03 * tc[1]
            -9.98825771e-07 * tc[2]
            +2.30053008e-10 * tc[3]
            -2.03647716e-14 * tc[4];
        /*species 5: CO2 */
        species[5] =
            +2.85746029e+00
            +4.41437026e-03 * tc[1]
            -2.21481404e-06 * tc[2]
            +5.23490188e-10 * tc[3]
            -4.72084164e-14 * tc[4];
        /*species 6: H2 */
        species[6] =
            +2.33727920e+00
            -4.94024731e-05 * tc[1]
            +4.99456778e-07 * tc[2]
            -1.79566394e-10 * tc[3]
            +2.00255376e-14 * tc[4];
    }
    return;
}


/*compute Cp/R at the given temperature */
/*tc contains precomputed powers of T, tc[0] = log(T) */
void cp_R(double *  species, double *  tc)
{

    /*temperature */
    double T = tc[1];

    /*species with midpoint at T=1000 kelvin */
    if (T < 1000) {
        /*species 0: CH4 */
        species[0] =
            +5.14987613e+00
            -1.36709788e-02 * tc[1]
            +4.91800599e-05 * tc[2]
            -4.84743026e-08 * tc[3]
            +1.66693956e-11 * tc[4];
        /*species 1: O2 */
        species[1] =
            +3.78245636e+00
            -2.99673416e-03 * tc[1]
            +9.84730201e-06 * tc[2]
            -9.68129509e-09 * tc[3]
            +3.24372837e-12 * tc[4];
        /*species 2: H2O */
        species[2] =
            +4.19864056e+00
            -2.03643410e-03 * tc[1]
            +6.52040211e-06 * tc[2]
            -5.48797062e-09 * tc[3]
            +1.77197817e-12 * tc[4];
        /*species 3: N2 */
        species[3] =
            +3.29867700e+00
            +1.40824040e-03 * tc[1]
            -3.96322200e-06 * tc[2]
            +5.64151500e-09 * tc[3]
            -2.44485400e-12 * tc[4];
        /*species 4: CO */
        species[4] =
            +3.57953347e+00
            -6.10353680e-04 * tc[1]
            +1.01681433e-06 * tc[2]
            +9.07005884e-10 * tc[3]
            -9.04424499e-13 * tc[4];
        /*species 5: CO2 */
        species[5] =
            +2.35677352e+00
            +8.98459677e-03 * tc[1]
            -7.12356269e-06 * tc[2]
            +2.45919022e-09 * tc[3]
            -1.43699548e-13 * tc[4];
        /*species 6: H2 */
        species[6] =
            +2.34433112e+00
            +7.98052075e-03 * tc[1]
            -1.94781510e-05 * tc[2]
            +2.01572094e-08 * tc[3]
            -7.37611761e-12 * tc[4];
    } else {
        /*species 0: CH4 */
        species[0] =
            +7.48514950e-02
            +1.33909467e-02 * tc[1]
            -5.73285809e-06 * tc[2]
            +1.22292535e-09 * tc[3]
            -1.01815230e-13 * tc[4];
        /*species 1: O2 */
        species[1] =
            +3.28253784e+00
            +1.48308754e-03 * tc[1]
            -7.57966669e-07 * tc[2]
            +2.09470555e-10 * tc[3]
            -2.16717794e-14 * tc[4];
        /*species 2: H2O */
        species[2] =
            +3.03399249e+00
            +2.17691804e-03 * tc[1]
            -1.64072518e-07 * tc[2]
            -9.70419870e-11 * tc[3]
            +1.68200992e-14 * tc[4];
        /*species 3: N2 */
        species[3] =
            +2.92664000e+00
            +1.48797680e-03 * tc[1]
            -5.68476000e-07 * tc[2]
            +1.00970380e-10 * tc[3]
            -6.75335100e-15 * tc[4];
        /*species 4: CO */
        species[4] =
            +2.71518561e+00
            +2.06252743e-03 * tc[1]
            -9.98825771e-07 * tc[2]
            +2.30053008e-10 * tc[3]
            -2.03647716e-14 * tc[4];
        /*species 5: CO2 */
        species[5] =
            +3.85746029e+00
            +4.41437026e-03 * tc[1]
            -2.21481404e-06 * tc[2]
            +5.23490188e-10 * tc[3]
            -4.72084164e-14 * tc[4];
        /*species 6: H2 */
        species[6] =
            +3.33727920e+00
            -4.94024731e-05 * tc[1]
            +4.99456778e-07 * tc[2]
            -1.79566394e-10 * tc[3]
            +2.00255376e-14 * tc[4];
    }
    return;
}


/*compute the e/(RT) at the given temperature */
/*tc contains precomputed powers of T, tc[0] = log(T) */
void speciesInternalEnergy(double *  species, double *  tc)
{

    /*temperature */
    double T = tc[1];
    double invT = 1 / T;

    /*species with midpoint at T=1000 kelvin */
    if (T < 1000) {
        /*species 0: CH4 */
        species[0] =
            +4.14987613e+00
            -6.83548940e-03 * tc[1]
            +1.63933533e-05 * tc[2]
            -1.21185757e-08 * tc[3]
            +3.33387912e-12 * tc[4]
            -1.02466476e+04 * invT;
        /*species 1: O2 */
        species[1] =
            +2.78245636e+00
            -1.49836708e-03 * tc[1]
            +3.28243400e-06 * tc[2]
            -2.42032377e-09 * tc[3]
            +6.48745674e-13 * tc[4]
            -1.06394356e+03 * invT;
        /*species 2: H2O */
        species[2] =
            +3.19864056e+00
            -1.01821705e-03 * tc[1]
            +2.17346737e-06 * tc[2]
            -1.37199266e-09 * tc[3]
            +3.54395634e-13 * tc[4]
            -3.02937267e+04 * invT;
        /*species 3: N2 */
        species[3] =
            +2.29867700e+00
            +7.04120200e-04 * tc[1]
            -1.32107400e-06 * tc[2]
            +1.41037875e-09 * tc[3]
            -4.88970800e-13 * tc[4]
            -1.02089990e+03 * invT;
        /*species 4: CO */
        species[4] =
            +2.57953347e+00
            -3.05176840e-04 * tc[1]
            +3.38938110e-07 * tc[2]
            +2.26751471e-10 * tc[3]
            -1.80884900e-13 * tc[4]
            -1.43440860e+04 * invT;
        /*species 5: CO2 */
        species[5] =
            +1.35677352e+00
            +4.49229839e-03 * tc[1]
            -2.37452090e-06 * tc[2]
            +6.14797555e-10 * tc[3]
            -2.87399096e-14 * tc[4]
            -4.83719697e+04 * invT;
        /*species 6: H2 */
        species[6] =
            +1.34433112e+00
            +3.99026037e-03 * tc[1]
            -6.49271700e-06 * tc[2]
            +5.03930235e-09 * tc[3]
            -1.47522352e-12 * tc[4]
            -9.17935173e+02 * invT;
    } else {
        /*species 0: CH4 */
        species[0] =
            -9.25148505e-01
            +6.69547335e-03 * tc[1]
            -1.91095270e-06 * tc[2]
            +3.05731338e-10 * tc[3]
            -2.03630460e-14 * tc[4]
            -9.46834459e+03 * invT;
        /*species 1: O2 */
        species[1] =
            +2.28253784e+00
            +7.41543770e-04 * tc[1]
            -2.52655556e-07 * tc[2]
            +5.23676387e-11 * tc[3]
            -4.33435588e-15 * tc[4]
            -1.08845772e+03 * invT;
        /*species 2: H2O */
        species[2] =
            +2.03399249e+00
            +1.08845902e-03 * tc[1]
            -5.46908393e-08 * tc[2]
            -2.42604967e-11 * tc[3]
            +3.36401984e-15 * tc[4]
            -3.00042971e+04 * invT;
        /*species 3: N2 */
        species[3] =
            +1.92664000e+00
            +7.43988400e-04 * tc[1]
            -1.89492000e-07 * tc[2]
            +2.52425950e-11 * tc[3]
            -1.35067020e-15 * tc[4]
            -9.22797700e+02 * invT;
        /*species 4: CO */
        species[4] =
            +1.71518561e+00
            +1.03126372e-03 * tc[1]
            -3.32941924e-07 * tc[2]
            +5.75132520e-11 * tc[3]
            -4.07295432e-15 * tc[4]
            -1.41518724e+04 * invT;
        /*species 5: CO2 */
        species[5] =
            +2.85746029e+00
            +2.20718513e-03 * tc[1]
            -7.38271347e-07 * tc[2]
            +1.30872547e-10 * tc[3]
            -9.44168328e-15 * tc[4]
            -4.87591660e+04 * invT;
        /*species 6: H2 */
        species[6] =
            +2.33727920e+00
            -2.47012365e-05 * tc[1]
            +1.66485593e-07 * tc[2]
            -4.48915985e-11 * tc[3]
            +4.00510752e-15 * tc[4]
            -9.50158922e+02 * invT;
    }
    return;
}


/*compute the h/(RT) at the given temperature (Eq 20) */
/*tc contains precomputed powers of T, tc[0] = log(T) */
void speciesEnthalpy(double *  species, double *  tc)
{

    /*temperature */
    double T = tc[1];
    double invT = 1 / T;

    /*species with midpoint at T=1000 kelvin */
    if (T < 1000) {
        /*species 0: CH4 */
        species[0] =
            +5.14987613e+00
            -6.83548940e-03 * tc[1]
            +1.63933533e-05 * tc[2]
            -1.21185757e-08 * tc[3]
            +3.33387912e-12 * tc[4]
            -1.02466476e+04 * invT;
        /*species 1: O2 */
        species[1] =
            +3.78245636e+00
            -1.49836708e-03 * tc[1]
            +3.28243400e-06 * tc[2]
            -2.42032377e-09 * tc[3]
            +6.48745674e-13 * tc[4]
            -1.06394356e+03 * invT;
        /*species 2: H2O */
        species[2] =
            +4.19864056e+00
            -1.01821705e-03 * tc[1]
            +2.17346737e-06 * tc[2]
            -1.37199266e-09 * tc[3]
            +3.54395634e-13 * tc[4]
            -3.02937267e+04 * invT;
        /*species 3: N2 */
        species[3] =
            +3.29867700e+00
            +7.04120200e-04 * tc[1]
            -1.32107400e-06 * tc[2]
            +1.41037875e-09 * tc[3]
            -4.88970800e-13 * tc[4]
            -1.02089990e+03 * invT;
        /*species 4: CO */
        species[4] =
            +3.57953347e+00
            -3.05176840e-04 * tc[1]
            +3.38938110e-07 * tc[2]
            +2.26751471e-10 * tc[3]
            -1.80884900e-13 * tc[4]
            -1.43440860e+04 * invT;
        /*species 5: CO2 */
        species[5] =
            +2.35677352e+00
            +4.49229839e-03 * tc[1]
            -2.37452090e-06 * tc[2]
            +6.14797555e-10 * tc[3]
            -2.87399096e-14 * tc[4]
            -4.83719697e+04 * invT;
        /*species 6: H2 */
        species[6] =
            +2.34433112e+00
            +3.99026037e-03 * tc[1]
            -6.49271700e-06 * tc[2]
            +5.03930235e-09 * tc[3]
            -1.47522352e-12 * tc[4]
            -9.17935173e+02 * invT;
    } else {
        /*species 0: CH4 */
        species[0] =
            +7.48514950e-02
            +6.69547335e-03 * tc[1]
            -1.91095270e-06 * tc[2]
            +3.05731338e-10 * tc[3]
            -2.03630460e-14 * tc[4]
            -9.46834459e+03 * invT;
        /*species 1: O2 */
        species[1] =
            +3.28253784e+00
            +7.41543770e-04 * tc[1]
            -2.52655556e-07 * tc[2]
            +5.23676387e-11 * tc[3]
            -4.33435588e-15 * tc[4]
            -1.08845772e+03 * invT;
        /*species 2: H2O */
        species[2] =
            +3.03399249e+00
            +1.08845902e-03 * tc[1]
            -5.46908393e-08 * tc[2]
            -2.42604967e-11 * tc[3]
            +3.36401984e-15 * tc[4]
            -3.00042971e+04 * invT;
        /*species 3: N2 */
        species[3] =
            +2.92664000e+00
            +7.43988400e-04 * tc[1]
            -1.89492000e-07 * tc[2]
            +2.52425950e-11 * tc[3]
            -1.35067020e-15 * tc[4]
            -9.22797700e+02 * invT;
        /*species 4: CO */
        species[4] =
            +2.71518561e+00
            +1.03126372e-03 * tc[1]
            -3.32941924e-07 * tc[2]
            +5.75132520e-11 * tc[3]
            -4.07295432e-15 * tc[4]
            -1.41518724e+04 * invT;
        /*species 5: CO2 */
        species[5] =
            +3.85746029e+00
            +2.20718513e-03 * tc[1]
            -7.38271347e-07 * tc[2]
            +1.30872547e-10 * tc[3]
            -9.44168328e-15 * tc[4]
            -4.87591660e+04 * invT;
        /*species 6: H2 */
        species[6] =
            +3.33727920e+00
            -2.47012365e-05 * tc[1]
            +1.66485593e-07 * tc[2]
            -4.48915985e-11 * tc[3]
            +4.00510752e-15 * tc[4]
            -9.50158922e+02 * invT;
    }
    return;
}


/*compute the S/R at the given temperature (Eq 21) */
/*tc contains precomputed powers of T, tc[0] = log(T) */
void speciesEntropy(double *  species, double *  tc)
{

    /*temperature */
    double T = tc[1];

    /*species with midpoint at T=1000 kelvin */
    if (T < 1000) {
        /*species 0: CH4 */
        species[0] =
            +5.14987613e+00 * tc[0]
            -1.36709788e-02 * tc[1]
            +2.45900299e-05 * tc[2]
            -1.61581009e-08 * tc[3]
            +4.16734890e-12 * tc[4]
            -4.64130376e+00 ;
        /*species 1: O2 */
        species[1] =
            +3.78245636e+00 * tc[0]
            -2.99673416e-03 * tc[1]
            +4.92365101e-06 * tc[2]
            -3.22709836e-09 * tc[3]
            +8.10932092e-13 * tc[4]
            +3.65767573e+00 ;
        /*species 2: H2O */
        species[2] =
            +4.19864056e+00 * tc[0]
            -2.03643410e-03 * tc[1]
            +3.26020105e-06 * tc[2]
            -1.82932354e-09 * tc[3]
            +4.42994543e-13 * tc[4]
            -8.49032208e-01 ;
        /*species 3: N2 */
        species[3] =
            +3.29867700e+00 * tc[0]
            +1.40824040e-03 * tc[1]
            -1.98161100e-06 * tc[2]
            +1.88050500e-09 * tc[3]
            -6.11213500e-13 * tc[4]
            +3.95037200e+00 ;
        /*species 4: CO */
        species[4] =
            +3.57953347e+00 * tc[0]
            -6.10353680e-04 * tc[1]
            +5.08407165e-07 * tc[2]
            +3.02335295e-10 * tc[3]
            -2.26106125e-13 * tc[4]
            +3.50840928e+00 ;
        /*species 5: CO2 */
        species[5] =
            +2.35677352e+00 * tc[0]
            +8.98459677e-03 * tc[1]
            -3.56178134e-06 * tc[2]
            +8.19730073e-10 * tc[3]
            -3.59248870e-14 * tc[4]
            +9.90105222e+00 ;
        /*species 6: H2 */
        species[6] =
            +2.34433112e+00 * tc[0]
            +7.98052075e-03 * tc[1]
            -9.73907550e-06 * tc[2]
            +6.71906980e-09 * tc[3]
            -1.84402940e-12 * tc[4]
            +6.83010238e-01 ;
    } else {
        /*species 0: CH4 */
        species[0] =
            +7.48514950e-02 * tc[0]
            +1.33909467e-02 * tc[1]
            -2.86642905e-06 * tc[2]
            +4.07641783e-10 * tc[3]
            -2.54538075e-14 * tc[4]
            +1.84373180e+01 ;
        /*species 1: O2 */
        species[1] =
            +3.28253784e+00 * tc[0]
            +1.48308754e-03 * tc[1]
            -3.78983334e-07 * tc[2]
            +6.98235183e-11 * tc[3]
            -5.41794485e-15 * tc[4]
            +5.45323129e+00 ;
        /*species 2: H2O */
        species[2] =
            +3.03399249e+00 * tc[0]
            +2.17691804e-03 * tc[1]
            -8.20362590e-08 * tc[2]
            -3.23473290e-11 * tc[3]
            +4.20502480e-15 * tc[4]
            +4.96677010e+00 ;
        /*species 3: N2 */
        species[3] =
            +2.92664000e+00 * tc[0]
            +1.48797680e-03 * tc[1]
            -2.84238000e-07 * tc[2]
            +3.36567933e-11 * tc[3]
            -1.68833775e-15 * tc[4]
            +5.98052800e+00 ;
        /*species 4: CO */
        species[4] =
            +2.71518561e+00 * tc[0]
            +2.06252743e-03 * tc[1]
            -4.99412886e-07 * tc[2]
            +7.66843360e-11 * tc[3]
            -5.09119290e-15 * tc[4]
            +7.81868772e+00 ;
        /*species 5: CO2 */
        species[5] =
            +3.85746029e+00 * tc[0]
            +4.41437026e-03 * tc[1]
            -1.10740702e-06 * tc[2]
            +1.74496729e-10 * tc[3]
            -1.18021041e-14 * tc[4]
            +2.27163806e+00 ;
        /*species 6: H2 */
        species[6] =
            +3.33727920e+00 * tc[0]
            -4.94024731e-05 * tc[1]
            +2.49728389e-07 * tc[2]
            -5.98554647e-11 * tc[3]
            +5.00638440e-15 * tc[4]
            -3.20502331e+00 ;
    }
    return;
}


/*save molecular weights into array */
void molecularWeight(double *  wt)
{
    wt[0] = 16.043030; /*CH4 */
    wt[1] = 31.998800; /*O2 */
    wt[2] = 18.015340; /*H2O */
    wt[3] = 28.013400; /*N2 */
    wt[4] = 28.010550; /*CO */
    wt[5] = 44.009950; /*CO2 */
    wt[6] = 2.015940; /*H2 */

    return;
}


/*save atomic weights into array */
void atomicWeight(double *  awt)
{
    awt[0] = 12.011150; /*C */
    awt[1] = 15.999400; /*O */
    awt[2] = 1.007970; /*H */
    awt[3] = 14.006700; /*N */

    return;
}
/* get temperature given internal energy in mass units and mass fracs */
void GET_T_GIVEN_EY(double *  e, double *  y, double *  t, int * ierr)
{
#ifdef CONVERGENCE
    const int maxiter = 5000;
    const double tol  = 1.e-12;
#else
    const int maxiter = 200;
    const double tol  = 1.e-6;
#endif
    double ein  = *e;
    double tmin = 90;/*max lower bound for thermo def */
    double tmax = 4000;/*min upper bound for thermo def */
    double e1,emin,emax,cv,t1,dt;
    int i;/* loop counter */
    CKUBMS(&tmin, y, &emin);
    CKUBMS(&tmax, y, &emax);
    if (ein < emin) {
        /*Linear Extrapolation below tmin */
        CKCVBS(&tmin, y, &cv);
        *t = tmin - (emin-ein)/cv;
        *ierr = 1;
        return;
    }
    if (ein > emax) {
        /*Linear Extrapolation above tmax */
        CKCVBS(&tmax, y, &cv);
        *t = tmax - (emax-ein)/cv;
        *ierr = 1;
        return;
    }
    t1 = *t;
    if (t1 < tmin || t1 > tmax) {
        t1 = tmin + (tmax-tmin)/(emax-emin)*(ein-emin);
    }
    for (i = 0; i < maxiter; ++i) {
        CKUBMS(&t1,y,&e1);
        CKCVBS(&t1,y,&cv);
        dt = (ein - e1) / cv;
        if (dt > 100.) { dt = 100.; }
        else if (dt < -100.) { dt = -100.; }
        else if (fabs(dt) < tol) break;
        else if (t1+dt == t1) break;
        t1 += dt;
    }
    *t = t1;
    *ierr = 0;
    return;
}

/* get temperature given enthalpy in mass units and mass fracs */
void GET_T_GIVEN_HY(double *  h, double *  y, double *  t, int * ierr)
{
#ifdef CONVERGENCE
    const int maxiter = 5000;
    const double tol  = 1.e-12;
#else
    const int maxiter = 200;
    const double tol  = 1.e-6;
#endif
    double hin  = *h;
    double tmin = 90;/*max lower bound for thermo def */
    double tmax = 4000;/*min upper bound for thermo def */
    double h1,hmin,hmax,cp,t1,dt;
    int i;/* loop counter */
    CKHBMS(&tmin, y, &hmin);
    CKHBMS(&tmax, y, &hmax);
    if (hin < hmin) {
        /*Linear Extrapolation below tmin */
        CKCPBS(&tmin, y, &cp);
        *t = tmin - (hmin-hin)/cp;
        *ierr = 1;
        return;
    }
    if (hin > hmax) {
        /*Linear Extrapolation above tmax */
        CKCPBS(&tmax, y, &cp);
        *t = tmax - (hmax-hin)/cp;
        *ierr = 1;
        return;
    }
    t1 = *t;
    if (t1 < tmin || t1 > tmax) {
        t1 = tmin + (tmax-tmin)/(hmax-hmin)*(hin-hmin);
    }
    for (i = 0; i < maxiter; ++i) {
        CKHBMS(&t1,y,&h1);
        CKCPBS(&t1,y,&cp);
        dt = (hin - h1) / cp;
        if (dt > 100.) { dt = 100.; }
        else if (dt < -100.) { dt = -100.; }
        else if (fabs(dt) < tol) break;
        else if (t1+dt == t1) break;
        t1 += dt;
    }
    *t = t1;
    *ierr = 0;
    return;
}


/*compute the critical parameters for each species */
void GET_CRITPARAMS(double *  Tci, double *  ai, double *  bi, double *  acentric_i)
{

    double   EPS[7];
    double   SIG[7];
    double    wt[7];
    double avogadro = 6.02214199e23;
    double boltzmann = 1.3806503e-16; //we work in CGS
    double Rcst = 83.144598; //in bar [CGS] !

    egtransetEPS(EPS);
    egtransetSIG(SIG);
    molecularWeight(wt);

    /*species 0: CH4 */
    /*Imported from NIST */
    Tci[0] = 190.560000 ; 
    ai[0] = 1e6 * 0.42748 * pow(Rcst,2.0) * pow(Tci[0],2.0) / (pow(16.043030,2.0) * 45.990000); 
    bi[0] = 0.08664 * Rcst * Tci[0] / (16.043030 * 45.990000); 
    acentric_i[0] = 0.011000 ;

    /*species 1: O2 */
    /*Imported from NIST */
    Tci[1] = 154.581000 ; 
    ai[1] = 1e6 * 0.42748 * pow(Rcst,2.0) * pow(Tci[1],2.0) / (pow(31.998800,2.0) * 50.430466); 
    bi[1] = 0.08664 * Rcst * Tci[1] / (31.998800 * 50.430466); 
    acentric_i[1] = 0.022200 ;

    /*species 2: H2O */
    /*Imported from NIST */
    Tci[2] = 647.096000 ; 
    ai[2] = 1e6 * 0.42748 * pow(Rcst,2.0) * pow(Tci[2],2.0) / (pow(18.015340,2.0) * 220.640000); 
    bi[2] = 0.08664 * Rcst * Tci[2] / (18.015340 * 220.640000); 
    acentric_i[2] = 0.344300 ;

    /*species 3: N2 */
    /*Imported from NIST */
    Tci[3] = 126.192000 ; 
    ai[3] = 1e6 * 0.42748 * pow(Rcst,2.0) * pow(Tci[3],2.0) / (pow(28.013400,2.0) * 33.958000); 
    bi[3] = 0.08664 * Rcst * Tci[3] / (28.013400 * 33.958000); 
    acentric_i[3] = 0.037200 ;

    /*species 4: CO */
    /*Imported from NIST */
    Tci[4] = 132.850000 ; 
    ai[4] = 1e6 * 0.42748 * pow(Rcst,2.0) * pow(Tci[4],2.0) / (pow(28.010000,2.0) * 34.940000); 
    bi[4] = 0.08664 * Rcst * Tci[4] / (28.010000 * 34.940000); 
    acentric_i[4] = 0.045000 ;

    /*species 5: CO2 */
    /*Imported from NIST */
    Tci[5] = 304.120000 ; 
    ai[5] = 1e6 * 0.42748 * pow(Rcst,2.0) * pow(Tci[5],2.0) / (pow(44.009950,2.0) * 73.740000); 
    bi[5] = 0.08664 * Rcst * Tci[5] / (44.009950 * 73.740000); 
    acentric_i[5] = 0.225000 ;

    /*species 6: H2 */
    /*Imported from NIST */
    Tci[6] = 33.145000 ; 
    ai[6] = 1e6 * 0.42748 * pow(Rcst,2.0) * pow(Tci[6],2.0) / (pow(2.015880,2.0) * 12.964000); 
    bi[6] = 0.08664 * Rcst * Tci[6] / (2.015880 * 12.964000); 
    acentric_i[6] = -0.219000 ;

    return;
}


#if defined(BL_FORT_USE_UPPERCASE)
#define egtransetLENIMC EGTRANSETLENIMC
#elif defined(BL_FORT_USE_LOWERCASE)
#define egtransetLENIMC egtransetlenimc
#elif defined(BL_FORT_USE_UNDERSCORE)
#define egtransetLENIMC egtransetlenimc_
#endif
void egtransetLENIMC(int* LENIMC ) {
    *LENIMC = 29;}


#if defined(BL_FORT_USE_UPPERCASE)
#define egtransetLENRMC EGTRANSETLENRMC
#elif defined(BL_FORT_USE_LOWERCASE)
#define egtransetLENRMC egtransetlenrmc
#elif defined(BL_FORT_USE_UNDERSCORE)
#define egtransetLENRMC egtransetlenrmc_
#endif
void egtransetLENRMC(int* LENRMC ) {
    *LENRMC = 1148;}


#if defined(BL_FORT_USE_UPPERCASE)
#define egtransetNO EGTRANSETNO
#elif defined(BL_FORT_USE_LOWERCASE)
#define egtransetNO egtransetno
#elif defined(BL_FORT_USE_UNDERSCORE)
#define egtransetNO egtransetno_
#endif
void egtransetNO(int* NO ) {
    *NO = 4;}


#if defined(BL_FORT_USE_UPPERCASE)
#define egtransetKK EGTRANSETKK
#elif defined(BL_FORT_USE_LOWERCASE)
#define egtransetKK egtransetkk
#elif defined(BL_FORT_USE_UNDERSCORE)
#define egtransetKK egtransetkk_
#endif
void egtransetKK(int* KK ) {
    *KK = 7;}


#if defined(BL_FORT_USE_UPPERCASE)
#define egtransetNLITE EGTRANSETNLITE
#elif defined(BL_FORT_USE_LOWERCASE)
#define egtransetNLITE egtransetnlite
#elif defined(BL_FORT_USE_UNDERSCORE)
#define egtransetNLITE egtransetnlite_
#endif
void egtransetNLITE(int* NLITE ) {
    *NLITE = 1;}


/*Patm in ergs/cm3 */
#if defined(BL_FORT_USE_UPPERCASE)
#define egtransetPATM EGTRANSETPATM
#elif defined(BL_FORT_USE_LOWERCASE)
#define egtransetPATM egtransetpatm
#elif defined(BL_FORT_USE_UNDERSCORE)
#define egtransetPATM egtransetpatm_
#endif
void egtransetPATM(double* PATM) {
    *PATM =   0.1013250000000000E+07;}


/*the molecular weights in g/mol */
#if defined(BL_FORT_USE_UPPERCASE)
#define egtransetWT EGTRANSETWT
#elif defined(BL_FORT_USE_LOWERCASE)
#define egtransetWT egtransetwt
#elif defined(BL_FORT_USE_UNDERSCORE)
#define egtransetWT egtransetwt_
#endif
void egtransetWT(double* WT ) {
    WT[0] = 1.60430300E+01;
    WT[1] = 3.19988000E+01;
    WT[2] = 1.80153400E+01;
    WT[3] = 2.80134000E+01;
    WT[4] = 2.80105500E+01;
    WT[5] = 4.40099500E+01;
    WT[6] = 2.01594000E+00;
}


/*the lennard-jones potential well depth eps/kb in K */
#if defined(BL_FORT_USE_UPPERCASE)
#define egtransetEPS EGTRANSETEPS
#elif defined(BL_FORT_USE_LOWERCASE)
#define egtransetEPS egtranseteps
#elif defined(BL_FORT_USE_UNDERSCORE)
#define egtransetEPS egtranseteps_
#endif
void egtransetEPS(double* EPS ) {
    EPS[2] = 5.72400000E+02;
    EPS[6] = 3.80000000E+01;
    EPS[0] = 1.41400000E+02;
    EPS[1] = 1.07400000E+02;
    EPS[4] = 9.81000000E+01;
    EPS[5] = 2.44000000E+02;
    EPS[3] = 9.75300000E+01;
}


/*the lennard-jones collision diameter in Angstroms */
#if defined(BL_FORT_USE_UPPERCASE)
#define egtransetSIG EGTRANSETSIG
#elif defined(BL_FORT_USE_LOWERCASE)
#define egtransetSIG egtransetsig
#elif defined(BL_FORT_USE_UNDERSCORE)
#define egtransetSIG egtransetsig_
#endif
void egtransetSIG(double* SIG ) {
    SIG[2] = 2.60500000E+00;
    SIG[6] = 2.92000000E+00;
    SIG[0] = 3.74600000E+00;
    SIG[1] = 3.45800000E+00;
    SIG[4] = 3.65000000E+00;
    SIG[5] = 3.76300000E+00;
    SIG[3] = 3.62100000E+00;
}


/*the dipole moment in Debye */
#if defined(BL_FORT_USE_UPPERCASE)
#define egtransetDIP EGTRANSETDIP
#elif defined(BL_FORT_USE_LOWERCASE)
#define egtransetDIP egtransetdip
#elif defined(BL_FORT_USE_UNDERSCORE)
#define egtransetDIP egtransetdip_
#endif
void egtransetDIP(double* DIP ) {
    DIP[2] = 1.84400000E+00;
    DIP[6] = 0.00000000E+00;
    DIP[0] = 0.00000000E+00;
    DIP[1] = 0.00000000E+00;
    DIP[4] = 0.00000000E+00;
    DIP[5] = 0.00000000E+00;
    DIP[3] = 0.00000000E+00;
}


/*the polarizability in cubic Angstroms */
#if defined(BL_FORT_USE_UPPERCASE)
#define egtransetPOL EGTRANSETPOL
#elif defined(BL_FORT_USE_LOWERCASE)
#define egtransetPOL egtransetpol
#elif defined(BL_FORT_USE_UNDERSCORE)
#define egtransetPOL egtransetpol_
#endif
void egtransetPOL(double* POL ) {
    POL[2] = 0.00000000E+00;
    POL[6] = 7.90000000E-01;
    POL[0] = 2.60000000E+00;
    POL[1] = 1.60000000E+00;
    POL[4] = 1.95000000E+00;
    POL[5] = 2.65000000E+00;
    POL[3] = 1.76000000E+00;
}


/*the rotational relaxation collision number at 298 K */
#if defined(BL_FORT_USE_UPPERCASE)
#define egtransetZROT EGTRANSETZROT
#elif defined(BL_FORT_USE_LOWERCASE)
#define egtransetZROT egtransetzrot
#elif defined(BL_FORT_USE_UNDERSCORE)
#define egtransetZROT egtransetzrot_
#endif
void egtransetZROT(double* ZROT ) {
    ZROT[2] = 4.00000000E+00;
    ZROT[6] = 2.80000000E+02;
    ZROT[0] = 1.30000000E+01;
    ZROT[1] = 3.80000000E+00;
    ZROT[4] = 1.80000000E+00;
    ZROT[5] = 2.10000000E+00;
    ZROT[3] = 4.00000000E+00;
}


/*0: monoatomic, 1: linear, 2: nonlinear */
#if defined(BL_FORT_USE_UPPERCASE)
#define egtransetNLIN EGTRANSETNLIN
#elif defined(BL_FORT_USE_LOWERCASE)
#define egtransetNLIN egtransetnlin
#elif defined(BL_FORT_USE_UNDERSCORE)
#define egtransetNLIN egtransetnlin_
#endif
void egtransetNLIN(int* NLIN) {
    NLIN[2] = 2;
    NLIN[6] = 1;
    NLIN[0] = 2;
    NLIN[1] = 1;
    NLIN[4] = 1;
    NLIN[5] = 1;
    NLIN[3] = 1;
}


/*Poly fits for the viscosities, dim NO*KK */
#if defined(BL_FORT_USE_UPPERCASE)
#define egtransetCOFETA EGTRANSETCOFETA
#elif defined(BL_FORT_USE_LOWERCASE)
#define egtransetCOFETA egtransetcofeta
#elif defined(BL_FORT_USE_UNDERSCORE)
#define egtransetCOFETA egtransetcofeta_
#endif
void egtransetCOFETA(double* COFETA) {
    COFETA[0] = -2.04642214E+01;
    COFETA[1] = 3.75363357E+00;
    COFETA[2] = -4.11756834E-01;
    COFETA[3] = 1.81766117E-02;
    COFETA[4] = -1.78915826E+01;
    COFETA[5] = 2.98311502E+00;
    COFETA[6] = -3.14105508E-01;
    COFETA[7] = 1.40500162E-02;
    COFETA[8] = -1.14441613E+01;
    COFETA[9] = -9.67162014E-01;
    COFETA[10] = 3.58651197E-01;
    COFETA[11] = -2.09789135E-02;
    COFETA[12] = -1.73976942E+01;
    COFETA[13] = 2.73482764E+00;
    COFETA[14] = -2.81916034E-01;
    COFETA[15] = 1.26588405E-02;
    COFETA[16] = -1.74469975E+01;
    COFETA[17] = 2.74728386E+00;
    COFETA[18] = -2.83509015E-01;
    COFETA[19] = 1.27267083E-02;
    COFETA[20] = -2.28110345E+01;
    COFETA[21] = 4.62954710E+00;
    COFETA[22] = -5.00689001E-01;
    COFETA[23] = 2.10012969E-02;
    COFETA[24] = -1.40419527E+01;
    COFETA[25] = 1.08789225E+00;
    COFETA[26] = -6.18592115E-02;
    COFETA[27] = 2.86838304E-03;
}


/*Poly fits for the conductivities, dim NO*KK */
#if defined(BL_FORT_USE_UPPERCASE)
#define egtransetCOFLAM EGTRANSETCOFLAM
#elif defined(BL_FORT_USE_LOWERCASE)
#define egtransetCOFLAM egtransetcoflam
#elif defined(BL_FORT_USE_UNDERSCORE)
#define egtransetCOFLAM egtransetcoflam_
#endif
void egtransetCOFLAM(double* COFLAM) {
    COFLAM[0] = 1.76957851E+01;
    COFLAM[1] = -6.71274475E+00;
    COFLAM[2] = 1.26299606E+00;
    COFLAM[3] = -6.62748072E-02;
    COFLAM[4] = 5.31578403E-01;
    COFLAM[5] = 1.87067453E+00;
    COFLAM[6] = -1.31586198E-01;
    COFLAM[7] = 5.22416151E-03;
    COFLAM[8] = 1.81612053E+01;
    COFLAM[9] = -6.74137053E+00;
    COFLAM[10] = 1.21372119E+00;
    COFLAM[11] = -6.11027962E-02;
    COFLAM[12] = 7.77703429E+00;
    COFLAM[13] = -1.30957605E+00;
    COFLAM[14] = 3.28842203E-01;
    COFLAM[15] = -1.69485484E-02;
    COFLAM[16] = 8.17515440E+00;
    COFLAM[17] = -1.53836440E+00;
    COFLAM[18] = 3.68036945E-01;
    COFLAM[19] = -1.90917513E-02;
    COFLAM[20] = -8.74831432E+00;
    COFLAM[21] = 4.79275291E+00;
    COFLAM[22] = -4.18685061E-01;
    COFLAM[23] = 1.35210242E-02;
    COFLAM[24] = 4.34729192E+00;
    COFLAM[25] = 1.55347646E+00;
    COFLAM[26] = -1.60615552E-01;
    COFLAM[27] = 9.89934485E-03;
}


/*Poly fits for the diffusion coefficients, dim NO*KK*KK */
#if defined(BL_FORT_USE_UPPERCASE)
#define egtransetCOFD EGTRANSETCOFD
#elif defined(BL_FORT_USE_LOWERCASE)
#define egtransetCOFD egtransetcofd
#elif defined(BL_FORT_USE_UNDERSCORE)
#define egtransetCOFD egtransetcofd_
#endif
void egtransetCOFD(double* COFD) {
    COFD[0] = -1.75618457E+01;
    COFD[1] = 4.30617914E+00;
    COFD[2] = -3.46490389E-01;
    COFD[3] = 1.51071405E-02;
    COFD[4] = -1.68102713E+01;
    COFD[5] = 4.01337907E+00;
    COFD[6] = -3.10488902E-01;
    COFD[7] = 1.36288975E-02;
    COFD[8] = -1.95657633E+01;
    COFD[9] = 4.78813636E+00;
    COFD[10] = -3.65976055E-01;
    COFD[11] = 1.42215137E-02;
    COFD[12] = -1.65706884E+01;
    COFD[13] = 3.92005093E+00;
    COFD[14] = -2.99040611E-01;
    COFD[15] = 1.31607610E-02;
    COFD[16] = -1.65940140E+01;
    COFD[17] = 3.92553905E+00;
    COFD[18] = -2.99706984E-01;
    COFD[19] = 1.31876655E-02;
    COFD[20] = -1.93937264E+01;
    COFD[21] = 4.87146645E+00;
    COFD[22] = -4.13323360E-01;
    COFD[23] = 1.77408400E-02;
    COFD[24] = -1.29544834E+01;
    COFD[25] = 2.96758239E+00;
    COFD[26] = -1.76586224E-01;
    COFD[27] = 7.90559536E-03;
    COFD[28] = -1.68102713E+01;
    COFD[29] = 4.01337907E+00;
    COFD[30] = -3.10488902E-01;
    COFD[31] = 1.36288975E-02;
    COFD[32] = -1.60936570E+01;
    COFD[33] = 3.70633871E+00;
    COFD[34] = -2.71897253E-01;
    COFD[35] = 1.20097588E-02;
    COFD[36] = -1.99035583E+01;
    COFD[37] = 5.01694644E+00;
    COFD[38] = -4.08963011E-01;
    COFD[39] = 1.66143416E-02;
    COFD[40] = -1.58214936E+01;
    COFD[41] = 3.60000113E+00;
    COFD[42] = -2.58255120E-01;
    COFD[43] = 1.14251480E-02;
    COFD[44] = -1.58458281E+01;
    COFD[45] = 3.60600362E+00;
    COFD[46] = -2.59019961E-01;
    COFD[47] = 1.14576923E-02;
    COFD[48] = -1.87634092E+01;
    COFD[49] = 4.61060397E+00;
    COFD[50] = -3.83564503E-01;
    COFD[51] = 1.66168246E-02;
    COFD[52] = -1.22181183E+01;
    COFD[53] = 2.70415313E+00;
    COFD[54] = -1.41236971E-01;
    COFD[55] = 6.32236816E-03;
    COFD[56] = -1.95657633E+01;
    COFD[57] = 4.78813636E+00;
    COFD[58] = -3.65976055E-01;
    COFD[59] = 1.42215137E-02;
    COFD[60] = -1.99035583E+01;
    COFD[61] = 5.01694644E+00;
    COFD[62] = -4.08963011E-01;
    COFD[63] = 1.66143416E-02;
    COFD[64] = -1.16123849E+01;
    COFD[65] = 8.27754782E-01;
    COFD[66] = 2.52262233E-01;
    COFD[67] = -1.62567414E-02;
    COFD[68] = -1.99472346E+01;
    COFD[69] = 5.05636623E+00;
    COFD[70] = -4.17733686E-01;
    COFD[71] = 1.71403506E-02;
    COFD[72] = -1.99647405E+01;
    COFD[73] = 5.05179386E+00;
    COFD[74] = -4.16351103E-01;
    COFD[75] = 1.70488551E-02;
    COFD[76] = -1.82187624E+01;
    COFD[77] = 3.93854160E+00;
    COFD[78] = -2.28424632E-01;
    COFD[79] = 7.18603342E-03;
    COFD[80] = -1.73864044E+01;
    COFD[81] = 4.71143088E+00;
    COFD[82] = -3.95288626E-01;
    COFD[83] = 1.70702272E-02;
    COFD[84] = -1.65706884E+01;
    COFD[85] = 3.92005093E+00;
    COFD[86] = -2.99040611E-01;
    COFD[87] = 1.31607610E-02;
    COFD[88] = -1.58214936E+01;
    COFD[89] = 3.60000113E+00;
    COFD[90] = -2.58255120E-01;
    COFD[91] = 1.14251480E-02;
    COFD[92] = -1.99472346E+01;
    COFD[93] = 5.05636623E+00;
    COFD[94] = -4.17733686E-01;
    COFD[95] = 1.71403506E-02;
    COFD[96] = -1.56019563E+01;
    COFD[97] = 3.51542686E+00;
    COFD[98] = -2.47677471E-01;
    COFD[99] = 1.09841319E-02;
    COFD[100] = -1.56221627E+01;
    COFD[101] = 3.51977302E+00;
    COFD[102] = -2.48210923E-01;
    COFD[103] = 1.10059241E-02;
    COFD[104] = -1.85068873E+01;
    COFD[105] = 4.52122572E+00;
    COFD[106] = -3.73088946E-01;
    COFD[107] = 1.62076520E-02;
    COFD[108] = -1.20381391E+01;
    COFD[109] = 2.61421687E+00;
    COFD[110] = -1.28887086E-01;
    COFD[111] = 5.75609167E-03;
    COFD[112] = -1.65940140E+01;
    COFD[113] = 3.92553905E+00;
    COFD[114] = -2.99706984E-01;
    COFD[115] = 1.31876655E-02;
    COFD[116] = -1.58458281E+01;
    COFD[117] = 3.60600362E+00;
    COFD[118] = -2.59019961E-01;
    COFD[119] = 1.14576923E-02;
    COFD[120] = -1.99647405E+01;
    COFD[121] = 5.05179386E+00;
    COFD[122] = -4.16351103E-01;
    COFD[123] = 1.70488551E-02;
    COFD[124] = -1.56221627E+01;
    COFD[125] = 3.51977302E+00;
    COFD[126] = -2.48210923E-01;
    COFD[127] = 1.10059241E-02;
    COFD[128] = -1.56423580E+01;
    COFD[129] = 3.52412711E+00;
    COFD[130] = -2.48745351E-01;
    COFD[131] = 1.10277551E-02;
    COFD[132] = -1.85324360E+01;
    COFD[133] = 4.52748688E+00;
    COFD[134] = -3.73847542E-01;
    COFD[135] = 1.62384117E-02;
    COFD[136] = -1.20607690E+01;
    COFD[137] = 2.61969379E+00;
    COFD[138] = -1.29638429E-01;
    COFD[139] = 5.79050588E-03;
    COFD[140] = -1.93937264E+01;
    COFD[141] = 4.87146645E+00;
    COFD[142] = -4.13323360E-01;
    COFD[143] = 1.77408400E-02;
    COFD[144] = -1.87634092E+01;
    COFD[145] = 4.61060397E+00;
    COFD[146] = -3.83564503E-01;
    COFD[147] = 1.66168246E-02;
    COFD[148] = -1.82187624E+01;
    COFD[149] = 3.93854160E+00;
    COFD[150] = -2.28424632E-01;
    COFD[151] = 7.18603342E-03;
    COFD[152] = -1.85068873E+01;
    COFD[153] = 4.52122572E+00;
    COFD[154] = -3.73088946E-01;
    COFD[155] = 1.62076520E-02;
    COFD[156] = -1.85324360E+01;
    COFD[157] = 4.52748688E+00;
    COFD[158] = -3.73847542E-01;
    COFD[159] = 1.62384117E-02;
    COFD[160] = -2.05810669E+01;
    COFD[161] = 5.07469434E+00;
    COFD[162] = -4.25340301E-01;
    COFD[163] = 1.76800795E-02;
    COFD[164] = -1.43978662E+01;
    COFD[165] = 3.49721576E+00;
    COFD[166] = -2.45465191E-01;
    COFD[167] = 1.08948372E-02;
    COFD[168] = -1.29544834E+01;
    COFD[169] = 2.96758239E+00;
    COFD[170] = -1.76586224E-01;
    COFD[171] = 7.90559536E-03;
    COFD[172] = -1.22181183E+01;
    COFD[173] = 2.70415313E+00;
    COFD[174] = -1.41236971E-01;
    COFD[175] = 6.32236816E-03;
    COFD[176] = -1.73864044E+01;
    COFD[177] = 4.71143088E+00;
    COFD[178] = -3.95288626E-01;
    COFD[179] = 1.70702272E-02;
    COFD[180] = -1.20381391E+01;
    COFD[181] = 2.61421687E+00;
    COFD[182] = -1.28887086E-01;
    COFD[183] = 5.75609167E-03;
    COFD[184] = -1.20607690E+01;
    COFD[185] = 2.61969379E+00;
    COFD[186] = -1.29638429E-01;
    COFD[187] = 5.79050588E-03;
    COFD[188] = -1.43978662E+01;
    COFD[189] = 3.49721576E+00;
    COFD[190] = -2.45465191E-01;
    COFD[191] = 1.08948372E-02;
    COFD[192] = -1.04285080E+01;
    COFD[193] = 2.23477534E+00;
    COFD[194] = -8.11809423E-02;
    COFD[195] = 3.77342041E-03;
}


/*List of specs with small weight, dim NLITE */
#if defined(BL_FORT_USE_UPPERCASE)
#define egtransetKTDIF EGTRANSETKTDIF
#elif defined(BL_FORT_USE_LOWERCASE)
#define egtransetKTDIF egtransetktdif
#elif defined(BL_FORT_USE_UNDERSCORE)
#define egtransetKTDIF egtransetktdif_
#endif
void egtransetKTDIF(int* KTDIF) {
    KTDIF[0] = 7;
}


/*Poly fits for thermal diff ratios, dim NO*NLITE*KK */
#if defined(BL_FORT_USE_UPPERCASE)
#define egtransetCOFTD EGTRANSETCOFTD
#elif defined(BL_FORT_USE_LOWERCASE)
#define egtransetCOFTD egtransetcoftd
#elif defined(BL_FORT_USE_UNDERSCORE)
#define egtransetCOFTD egtransetcoftd_
#endif
void egtransetCOFTD(double* COFTD) {
    COFTD[0] = 2.98973958E-01;
    COFTD[1] = 2.32230992E-04;
    COFTD[2] = -1.23675023E-07;
    COFTD[3] = 2.01029713E-11;
    COFTD[4] = 3.81864172E-01;
    COFTD[5] = 1.84117353E-04;
    COFTD[6] = -9.79617476E-08;
    COFTD[7] = 1.62542227E-11;
    COFTD[8] = 1.95123509E-03;
    COFTD[9] = 6.69470998E-04;
    COFTD[10] = -3.12148757E-07;
    COFTD[11] = 4.52949938E-11;
    COFTD[12] = 3.88181125E-01;
    COFTD[13] = 1.55380218E-04;
    COFTD[14] = -8.20880914E-08;
    COFTD[15] = 1.37104636E-11;
    COFTD[16] = 3.87405318E-01;
    COFTD[17] = 1.56883797E-04;
    COFTD[18] = -8.29309791E-08;
    COFTD[19] = 1.38460299E-11;
    COFTD[20] = 2.47129011E-01;
    COFTD[21] = 4.49395677E-04;
    COFTD[22] = -2.32030740E-07;
    COFTD[23] = 3.62578797E-11;
    COFTD[24] = 0.00000000E+00;
    COFTD[25] = 0.00000000E+00;
    COFTD[26] = 0.00000000E+00;
    COFTD[27] = 0.00000000E+00;
}

}

/* End of file  */




#if 0




\\
\\
\\  This is the mechanism file
\\
\\
ELEMENTS
C O H N
END
SPECIES
CH4 O2 H2O N2 CO CO2 H2
END
TRANS ALL
CH4                2   141.400     3.746     0.000     2.600    13.000          
O2                 1   107.400     3.458     0.000     1.600     3.800          
H2O                2   572.400     2.605     1.844     0.000     4.000          
N2                 1    97.530     3.621     0.000     1.760     4.000          
CO                 1    98.100     3.650     0.000     1.950     1.800          
CO2                1   244.000     3.763     0.000     2.650     2.100          
H2                 1    38.000     2.920     0.000     0.790   280.000          
END
REACTIONS
2CH4+O2=>2CO+4H2                3.91E13 0 30.0E3
  FORD /CH4 0.5/
  FORD /O2 1.25/
CH4+H2O=CO+3H2                  3.0E11 0 30.0E3
2H2+O2=>2H2O                    6.045E18 -1 40.0E3
  FORD /H2 0.25/
  FORD /O2 1.5/
CO+H2O=CO2+H2                   2.75E12 0 20.0E3
END

\\
\\
\\  This is the therm file
\\
\\
THERMO ALL
   300.000  1000.000  5000.000
O                 L 1/90O   1               G   200.000  3500.000  1000.000    1
 2.56942078E+00-8.59741137E-05 4.19484589E-08-1.00177799E-11 1.22833691E-15    2
 2.92175791E+04 4.78433864E+00 3.16826710E+00-3.27931884E-03 6.64306396E-06    3
-6.12806624E-09 2.11265971E-12 2.91222592E+04 2.05193346E+00                   4
O2                TPIS89O   2               G   200.000  3500.000  1000.000    1
 3.28253784E+00 1.48308754E-03-7.57966669E-07 2.09470555E-10-2.16717794E-14    2
-1.08845772E+03 5.45323129E+00 3.78245636E+00-2.99673416E-03 9.84730201E-06    3
-9.68129509E-09 3.24372837E-12-1.06394356E+03 3.65767573E+00                   4
H                 L 7/88H   1               G   200.000  3500.000  1000.000    1
 2.50000001E+00-2.30842973E-11 1.61561948E-14-4.73515235E-18 4.98197357E-22    2
 2.54736599E+04-4.46682914E-01 2.50000000E+00 7.05332819E-13-1.99591964E-15    3
 2.30081632E-18-9.27732332E-22 2.54736599E+04-4.46682853E-01                   4
H2                TPIS78H   2               G   200.000  3500.000  1000.000    1
 3.33727920E+00-4.94024731E-05 4.99456778E-07-1.79566394E-10 2.00255376E-14    2
-9.50158922E+02-3.20502331E+00 2.34433112E+00 7.98052075E-03-1.94781510E-05    3
 2.01572094E-08-7.37611761E-12-9.17935173E+02 6.83010238E-01                   4
OH                RUS 78O   1H   1          G   200.000  3500.000  1000.000    1
 3.09288767E+00 5.48429716E-04 1.26505228E-07-8.79461556E-11 1.17412376E-14    2
 3.85865700E+03 4.47669610E+00 3.99201543E+00-2.40131752E-03 4.61793841E-06    3
-3.88113333E-09 1.36411470E-12 3.61508056E+03-1.03925458E-01                   4
H2O               L 8/89H   2O   1          G   200.000  3500.000  1000.000    1
 3.03399249E+00 2.17691804E-03-1.64072518E-07-9.70419870E-11 1.68200992E-14    2
-3.00042971E+04 4.96677010E+00 4.19864056E+00-2.03643410E-03 6.52040211E-06    3
-5.48797062E-09 1.77197817E-12-3.02937267E+04-8.49032208E-01                   4
HO2               L 5/89H   1O   2          G   200.000  3500.000  1000.000    1
 4.01721090E+00 2.23982013E-03-6.33658150E-07 1.14246370E-10-1.07908535E-14    2
 1.11856713E+02 3.78510215E+00 4.30179801E+00-4.74912051E-03 2.11582891E-05    3
-2.42763894E-08 9.29225124E-12 2.94808040E+02 3.71666245E+00                   4
H2O2              L 7/88H   2O   2          G   200.000  3500.000  1000.000    1
 4.16500285E+00 4.90831694E-03-1.90139225E-06 3.71185986E-10-2.87908305E-14    2
-1.78617877E+04 2.91615662E+00 4.27611269E+00-5.42822417E-04 1.67335701E-05    3
-2.15770813E-08 8.62454363E-12-1.77025821E+04 3.43505074E+00                   4
C                 L11/88C   1               G   200.000  3500.000  1000.000    1
 2.49266888E+00 4.79889284E-05-7.24335020E-08 3.74291029E-11-4.87277893E-15    2
 8.54512953E+04 4.80150373E+00 2.55423955E+00-3.21537724E-04 7.33792245E-07    3
-7.32234889E-10 2.66521446E-13 8.54438832E+04 4.53130848E+00                   4
CH                TPIS79C   1H   1          G   200.000  3500.000  1000.000    1
 2.87846473E+00 9.70913681E-04 1.44445655E-07-1.30687849E-10 1.76079383E-14    2
 7.10124364E+04 5.48497999E+00 3.48981665E+00 3.23835541E-04-1.68899065E-06    3
 3.16217327E-09-1.40609067E-12 7.07972934E+04 2.08401108E+00                   4
CH2               L S/93C   1H   2          G   200.000  3500.000  1000.000    1
 2.87410113E+00 3.65639292E-03-1.40894597E-06 2.60179549E-10-1.87727567E-14    2
 4.62636040E+04 6.17119324E+00 3.76267867E+00 9.68872143E-04 2.79489841E-06    3
-3.85091153E-09 1.68741719E-12 4.60040401E+04 1.56253185E+00                   4
CH2S              L S/93C   1H   2          G   200.000  3500.000  1000.000    1
 2.29203842E+00 4.65588637E-03-2.01191947E-06 4.17906000E-10-3.39716365E-14    2
 5.09259997E+04 8.62650169E+00 4.19860411E+00-2.36661419E-03 8.23296220E-06    3
-6.68815981E-09 1.94314737E-12 5.04968163E+04-7.69118967E-01                   4
CH3               L11/89C   1H   3          G   200.000  3500.000  1000.000    1
 2.28571772E+00 7.23990037E-03-2.98714348E-06 5.95684644E-10-4.67154394E-14    2
 1.67755843E+04 8.48007179E+00 3.67359040E+00 2.01095175E-03 5.73021856E-06    3
-6.87117425E-09 2.54385734E-12 1.64449988E+04 1.60456433E+00                   4
CH4               L 8/88C   1H   4          G   200.000  3500.000  1000.000    1
 7.48514950E-02 1.33909467E-02-5.73285809E-06 1.22292535E-09-1.01815230E-13    2
-9.46834459E+03 1.84373180E+01 5.14987613E+00-1.36709788E-02 4.91800599E-05    3
-4.84743026E-08 1.66693956E-11-1.02466476E+04-4.64130376E+00                   4
CO                TPIS79C   1O   1          G   200.000  3500.000  1000.000    1
 2.71518561E+00 2.06252743E-03-9.98825771E-07 2.30053008E-10-2.03647716E-14    2
-1.41518724E+04 7.81868772E+00 3.57953347E+00-6.10353680E-04 1.01681433E-06    3
 9.07005884E-10-9.04424499E-13-1.43440860E+04 3.50840928E+00                   4
CO2               L 7/88C   1O   2          G   200.000  3500.000  1000.000    1
 3.85746029E+00 4.41437026E-03-2.21481404E-06 5.23490188E-10-4.72084164E-14    2
-4.87591660E+04 2.27163806E+00 2.35677352E+00 8.98459677E-03-7.12356269E-06    3
 2.45919022E-09-1.43699548E-13-4.83719697E+04 9.90105222E+00                   4
HCO               L12/89H   1C   1O   1     G   200.000  3500.000  1000.000    1
 2.77217438E+00 4.95695526E-03-2.48445613E-06 5.89161778E-10-5.33508711E-14    2
 4.01191815E+03 9.79834492E+00 4.22118584E+00-3.24392532E-03 1.37799446E-05    3
-1.33144093E-08 4.33768865E-12 3.83956496E+03 3.39437243E+00                   4
CH2O              L 8/88H   2C   1O   1     G   200.000  3500.000  1000.000    1
 1.76069008E+00 9.20000082E-03-4.42258813E-06 1.00641212E-09-8.83855640E-14    2
-1.39958323E+04 1.36563230E+01 4.79372315E+00-9.90833369E-03 3.73220008E-05    3
-3.79285261E-08 1.31772652E-11-1.43089567E+04 6.02812900E-01                   4
CH2OH             GUNL93C   1H   3O   1     G   200.000  3500.000  1000.000    1
 3.69266569E+00 8.64576797E-03-3.75101120E-06 7.87234636E-10-6.48554201E-14    2
-3.24250627E+03 5.81043215E+00 3.86388918E+00 5.59672304E-03 5.93271791E-06    3
-1.04532012E-08 4.36967278E-12-3.19391367E+03 5.47302243E+00                   4
CH3O              121686C   1H   3O   1     G   200.00   3000.00   1000.000    1
 0.03770799E+02 0.07871497E-01-0.02656384E-04 0.03944431E-08-0.02112616E-12    2
 0.12783252E+03 0.02929575E+02 0.02106204E+02 0.07216595E-01 0.05338472E-04    3
-0.07377636E-07 0.02075610E-10 0.09786011E+04 0.13152177E+02                   4
CH3OH             L 8/88C   1H   4O   1     G   200.000  3500.000  1000.000    1
 1.78970791E+00 1.40938292E-02-6.36500835E-06 1.38171085E-09-1.17060220E-13    2
-2.53748747E+04 1.45023623E+01 5.71539582E+00-1.52309129E-02 6.52441155E-05    3
-7.10806889E-08 2.61352698E-11-2.56427656E+04-1.50409823E+00                   4
C2H               L 1/91C   2H   1          G   200.000  3500.000  1000.000    1
 3.16780652E+00 4.75221902E-03-1.83787077E-06 3.04190252E-10-1.77232770E-14    2
 6.71210650E+04 6.63589475E+00 2.88965733E+00 1.34099611E-02-2.84769501E-05    3
 2.94791045E-08-1.09331511E-11 6.68393932E+04 6.22296438E+00                   4
C2H2              L 1/91C   2H   2          G   200.000  3500.000  1000.000    1
 4.14756964E+00 5.96166664E-03-2.37294852E-06 4.67412171E-10-3.61235213E-14    2
 2.59359992E+04-1.23028121E+00 8.08681094E-01 2.33615629E-02-3.55171815E-05    3
 2.80152437E-08-8.50072974E-12 2.64289807E+04 1.39397051E+01                   4
C2H3              L 2/92C   2H   3          G   200.000  3500.000  1000.000    1
 3.01672400E+00 1.03302292E-02-4.68082349E-06 1.01763288E-09-8.62607041E-14    2
 3.46128739E+04 7.78732378E+00 3.21246645E+00 1.51479162E-03 2.59209412E-05    3
-3.57657847E-08 1.47150873E-11 3.48598468E+04 8.51054025E+00                   4
C2H4              L 1/91C   2H   4          G   200.000  3500.000  1000.000    1
 2.03611116E+00 1.46454151E-02-6.71077915E-06 1.47222923E-09-1.25706061E-13    2
 4.93988614E+03 1.03053693E+01 3.95920148E+00-7.57052247E-03 5.70990292E-05    3
-6.91588753E-08 2.69884373E-11 5.08977593E+03 4.09733096E+00                   4
C2H5              L12/92C   2H   5          G   200.000  3500.000  1000.000    1
 1.95465642E+00 1.73972722E-02-7.98206668E-06 1.75217689E-09-1.49641576E-13    2
 1.28575200E+04 1.34624343E+01 4.30646568E+00-4.18658892E-03 4.97142807E-05    3
-5.99126606E-08 2.30509004E-11 1.28416265E+04 4.70720924E+00                   4
C2H6              L 8/88C   2H   6          G   200.000  3500.000  1000.000    1
 1.07188150E+00 2.16852677E-02-1.00256067E-05 2.21412001E-09-1.90002890E-13    2
-1.14263932E+04 1.51156107E+01 4.29142492E+00-5.50154270E-03 5.99438288E-05    3
-7.08466285E-08 2.68685771E-11-1.15222055E+04 2.66682316E+00                   4
CH2CO             L 5/90C   2H   2O   1     G   200.000  3500.000  1000.000    1
 4.51129732E+00 9.00359745E-03-4.16939635E-06 9.23345882E-10-7.94838201E-14    2
-7.55105311E+03 6.32247205E-01 2.13583630E+00 1.81188721E-02-1.73947474E-05    3
 9.34397568E-09-2.01457615E-12-7.04291804E+03 1.22156480E+01                   4
HCCO              SRIC91H   1C   2O   1     G   200.00   4000.00   1000.000    1
 0.56282058E+01 0.40853401E-02-0.15934547E-05 0.28626052E-09-0.19407832E-13    2
 0.19327215E+05-0.39302595E+01 0.22517214E+01 0.17655021E-01-0.23729101E-04    3
 0.17275759E-07-0.50664811E-11 0.20059449E+05 0.12490417E+02                   4
HCCOH              SRI91C   2O   1H   2     G   200.000  5000.000  1000.000    1
 0.59238291E+01 0.67923600E-02-0.25658564E-05 0.44987841E-09-0.29940101E-13    2
 0.72646260E+04-0.76017742E+01 0.12423733E+01 0.31072201E-01-0.50866864E-04    3
 0.43137131E-07-0.14014594E-10 0.80316143E+04 0.13874319E+02                   4
H2CN               41687H   2C   1N   1     G   200.00   4000.000  1000.000    1
 0.52097030E+01 0.29692911E-02-0.28555891E-06-0.16355500E-09 0.30432589E-13    2
 0.27677109E+05-0.44444780E+01 0.28516610E+01 0.56952331E-02 0.10711400E-05    3
-0.16226120E-08-0.23511081E-12 0.28637820E+05 0.89927511E+01                   4
HCN               GRI/98H   1C   1N   1     G   200.000  6000.000  1000.000    1
 0.38022392E+01 0.31464228E-02-0.10632185E-05 0.16619757E-09-0.97997570E-14    2
 0.14407292E+05 0.15754601E+01 0.22589886E+01 0.10051170E-01-0.13351763E-04    3
 0.10092349E-07-0.30089028E-11 0.14712633E+05 0.89164419E+01                   4
HNO               And93 H   1N   1O   1     G   200.000  6000.000  1000.000    1
 0.29792509E+01 0.34944059E-02-0.78549778E-06 0.57479594E-10-0.19335916E-15    2
 0.11750582E+05 0.86063728E+01 0.45334916E+01-0.56696171E-02 0.18473207E-04    3
-0.17137094E-07 0.55454573E-11 0.11548297E+05 0.17498417E+01                   4
N                 L 6/88N   1               G   200.000  6000.000  1000.000    1
 0.24159429E+01 0.17489065E-03-0.11902369E-06 0.30226245E-10-0.20360982E-14    2
 0.56133773E+05 0.46496096E+01 0.25000000E+01 0.00000000E+00 0.00000000E+00    3
 0.00000000E+00 0.00000000E+00 0.56104637E+05 0.41939087E+01                   4
NNH               T07/93N   2H   1          G   200.000  6000.000  1000.000    1
 0.37667544E+01 0.28915082E-02-0.10416620E-05 0.16842594E-09-0.10091896E-13    2
 0.28650697E+05 0.44705067E+01 0.43446927E+01-0.48497072E-02 0.20059459E-04    3
-0.21726464E-07 0.79469539E-11 0.28791973E+05 0.29779410E+01                   4
N2O               L 7/88N   2O   1          G   200.000  6000.000  1000.000    1
 0.48230729E+01 0.26270251E-02-0.95850874E-06 0.16000712E-09-0.97752303E-14    2
 0.80734048E+04-0.22017207E+01 0.22571502E+01 0.11304728E-01-0.13671319E-04    3
 0.96819806E-08-0.29307182E-11 0.87417744E+04 0.10757992E+02                   4
NH                And94 N   1H   1          G   200.000  6000.000  1000.000    1
 0.27836928E+01 0.13298430E-02-0.42478047E-06 0.78348501E-10-0.55044470E-14    2
 0.42120848E+05 0.57407799E+01 0.34929085E+01 0.31179198E-03-0.14890484E-05    3
 0.24816442E-08-0.10356967E-11 0.41880629E+05 0.18483278E+01                   4
NH2               And89 N   1H   2          G   200.000  6000.000  1000.000    1
 0.28347421E+01 0.32073082E-02-0.93390804E-06 0.13702953E-09-0.79206144E-14    2
 0.22171957E+05 0.65204163E+01 0.42040029E+01-0.21061385E-02 0.71068348E-05    3
-0.56115197E-08 0.16440717E-11 0.21885910E+05-0.14184248E+00                   4
NH3               J 6/77N   1H   3          G   200.000  6000.000  1000.000    1
 0.26344521E+01 0.56662560E-02-0.17278676E-05 0.23867161E-09-0.12578786E-13    2
-0.65446958E+04 0.65662928E+01 0.42860274E+01-0.46605230E-02 0.21718513E-04    3
-0.22808887E-07 0.82638046E-11-0.67417285E+04-0.62537277E+00                   4
NO                RUS 78N   1O   1          G   200.000  6000.000  1000.000    1
 0.32606056E+01 0.11911043E-02-0.42917048E-06 0.69457669E-10-0.40336099E-14    2
 0.99209746E+04 0.63693027E+01 0.42184763E+01-0.46389760E-02 0.11041022E-04    3
-0.93361354E-08 0.28035770E-11 0.98446230E+04 0.22808464E+01                   4
NO2               L 7/88N   1O   2          G   200.000  6000.000  1000.000    1
 0.48847542E+01 0.21723956E-02-0.82806906E-06 0.15747510E-09-0.10510895E-13    2
 0.23164983E+04-0.11741695E+00 0.39440312E+01-0.15854290E-02 0.16657812E-04    3
-0.20475426E-07 0.78350564E-11 0.28966179E+04 0.63119917E+01                   4
HCNO              BDEA94H   1N   1C   1O   1G   200.000  5000.000  1382.000    1
 6.59860456E+00 3.02778626E-03-1.07704346E-06 1.71666528E-10-1.01439391E-14    2
 1.79661339E+04-1.03306599E+01 2.64727989E+00 1.27505342E-02-1.04794236E-05    3
 4.41432836E-09-7.57521466E-13 1.92990252E+04 1.07332972E+01                   4
HOCN              BDEA94H   1N   1C   1O   1G   200.000  5000.000  1368.000    1
 5.89784885E+00 3.16789393E-03-1.11801064E-06 1.77243144E-10-1.04339177E-14    2
-3.70653331E+03-6.18167825E+00 3.78604952E+00 6.88667922E-03-3.21487864E-06    3
 5.17195767E-10 1.19360788E-14-2.82698400E+03 5.63292162E+00                   4
HNCO              BDEA94H   1N   1C   1O   1G   200.000  5000.000  1478.000    1
 6.22395134E+00 3.17864004E-03-1.09378755E-06 1.70735163E-10-9.95021955E-15    2
-1.66599344E+04-8.38224741E+00 3.63096317E+00 7.30282357E-03-2.28050003E-06    3
-6.61271298E-10 3.62235752E-13-1.55873636E+04 6.19457727E+00                   4
NCO               EA 93 N   1C   1O   1     G   200.000  6000.000  1000.000    1
 0.51521845E+01 0.23051761E-02-0.88033153E-06 0.14789098E-09-0.90977996E-14    2
 0.14004123E+05-0.25442660E+01 0.28269308E+01 0.88051688E-02-0.83866134E-05    3
 0.48016964E-08-0.13313595E-11 0.14682477E+05 0.95504646E+01                   4
CN                HBH92 C   1N   1          G   200.000  6000.000  1000.000    1
 0.37459805E+01 0.43450775E-04 0.29705984E-06-0.68651806E-10 0.44134173E-14    2
 0.51536188E+05 0.27867601E+01 0.36129351E+01-0.95551327E-03 0.21442977E-05    3
-0.31516323E-09-0.46430356E-12 0.51708340E+05 0.39804995E+01                   4
HCNN              SRI/94C   1N   2H   1     G   200.000  5000.000  1000.000    1
 0.58946362E+01 0.39895959E-02-0.15982380E-05 0.29249395E-09-0.20094686E-13    2
 0.53452941E+05-0.51030502E+01 0.25243194E+01 0.15960619E-01-0.18816354E-04    3
 0.12125540E-07-0.32357378E-11 0.54261984E+05 0.11675870E+02                   4
N2                121286N   2               G   200.000  5000.000  1000.000    1
 0.02926640E+02 0.14879768E-02-0.05684760E-05 0.10097038E-09-0.06753351E-13    2
-0.09227977E+04 0.05980528E+02 0.03298677E+02 0.14082404E-02-0.03963222E-04    3
 0.05641515E-07-0.02444854E-10-0.10208999E+04 0.03950372E+02                   4
AR                120186AR  1               G   200.000  5000.000  1000.000    1
 0.02500000E+02 0.00000000E+00 0.00000000E+00 0.00000000E+00 0.00000000E+00    2
-0.07453750E+04 0.04366000E+02 0.02500000E+02 0.00000000E+00 0.00000000E+00    3
 0.00000000E+00 0.00000000E+00-0.07453750E+04 0.04366000E+02                   4
C3H8              L 4/85C   3H   8          G   200.000  5000.000  1000.000    1
 0.75341368E+01 0.18872239E-01-0.62718491E-05 0.91475649E-09-0.47838069E-13    2
-0.16467516E+05-0.17892349E+02 0.93355381E+00 0.26424579E-01 0.61059727E-05    3
-0.21977499E-07 0.95149253E-11-0.13958520E+05 0.19201691E+02                   4
C3H7              L 9/84C   3H   7          G   200.000  5000.000  1000.000    1
 0.77026987E+01 0.16044203E-01-0.52833220E-05 0.76298590E-09-0.39392284E-13    2
 0.82984336E+04-0.15480180E+02 0.10515518E+01 0.25991980E-01 0.23800540E-05    3
-0.19609569E-07 0.93732470E-11 0.10631863E+05 0.21122559E+02                   4
CH3CHO            L 8/88C   2H   4O   1     G   200.000  6000.000  1000.000    1
 0.54041108E+01 0.11723059E-01-0.42263137E-05 0.68372451E-09-0.40984863E-13    2
-0.22593122E+05-0.34807917E+01 0.47294595E+01-0.31932858E-02 0.47534921E-04    3
-0.57458611E-07 0.21931112E-10-0.21572878E+05 0.41030159E+01                   4
CH2CHO            SAND86O   1H   3C   2     G   200.000  5000.000  1000.000    1
 0.05975670E+02 0.08130591E-01-0.02743624E-04 0.04070304E-08-0.02176017E-12    2
 0.04903218E+04-0.05045251E+02 0.03409062E+02 0.10738574E-01 0.01891492E-04    3
-0.07158583E-07 0.02867385E-10 0.15214766E+04 0.09558290E+02                   4
END





#endif
