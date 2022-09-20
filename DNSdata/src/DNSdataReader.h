/* DNSdataReader.h */
/* Wolfgang Tichy 8/2016 */

int DNS1_position_fileptr_after_str(FILE *in, char *str);
void DNS1_dataReader(CCTK_ARGUMENTS);
void DNS1_dataPars(CCTK_ARGUMENTS);
void DNS1_select_polytrope_n_kappa_k_of_hm1(double hm1,
                                           double *n, double *kappa, double *k);

void DNS1_set_TimeDeriv_in_inertFrame_assuming_HKV(CCTK_ARGUMENTS,
                                              CCTK_REAL *var, CCTK_REAL *dtvar,
                                              CCTK_REAL Omega);
