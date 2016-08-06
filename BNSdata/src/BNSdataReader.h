/* BNSdataReader.h */
/* Wolfgang Tichy 8/2016 */


void BNSdataReader(CCTK_ARGUMENTS);
void BNSdataPars(CCTK_ARGUMENTS);
void BNS_select_polytrope_n_kappa_k_of_hm1(double hm1,
                                           double *n, double *kappa, double *k);
