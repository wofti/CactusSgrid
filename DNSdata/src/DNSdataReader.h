/* BNSdataReader.h */
/* Wolfgang Tichy 8/2016 */


void DNSdataReader(CCTK_ARGUMENTS);
void DNSdataPars(CCTK_ARGUMENTS);
void DNS_select_polytrope_n_kappa_k_of_hm1(double hm1,
                                           double *n, double *kappa, double *k);
void DNS_set_TimeDeriv_in_inertFrame_assuming_HKV(CCTK_ARGUMENTS,
                                              CCTK_REAL *var, CCTK_REAL *dtvar,
                                              CCTK_REAL Omega);
int DNS_call_sgrid(const char *command);

/* functions and variables from sgrid */
int libsgrid_main(int argc, char **argv);
extern int SGRID_memory_persists;
int SGRID_grid_exists(void);
void SGRID_free_everything(void);

void SGRID_errorexits(char *file, int line, char *s, char *t);
#define SGRID_errorexits(s,t) SGRID_errorexits(__FILE__, __LINE__, (s), (t))

int SGRID_system2(char *s1, char *s2);
int SGRID_lock_curr_til_EOF(FILE *out);
int SGRID_construct_argv(char *str, char ***argv);

int SGRID_fgotonext(FILE *in, const char *label);
int SGRID_fgetparameter(FILE *in, const char *par, char *str);
int SGRID_extract_after_EQ(char *str);
int SGRID_extrstr_before_after_EQ(const char *str, char *before, char *after);
int SGRID_fscanline(FILE *in,char *str);
int SGRID_extrstr_before_after(const char *str, char *before, char *after, char z);
int SGRID_find_before_after(const char *str, char *before, char *after, const char *z);
int SGRID_pfind_before_after(const char *str,int p,char *before,char *after,const char *z);
int SGRID_sscan_word_at_p(const char *str, int p, char *word);
int SGRID_fscan_str_using_getc(FILE *in, char *str);
int SGRID_fscanf1(FILE *in, char *fmt, char *str);

void SGRID_EoS_T0_rho0_P_rhoE_from_hm1(double hm1,
                                       double *rho0, double *P, double *rhoE);
double SGRID_epsl_of_rho0_rhoE(double rho0, double rhoE);
