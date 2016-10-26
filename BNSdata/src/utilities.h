/* utilities.h */
/* Wolfgang Tichy, 2016 */

/* constants */
#ifdef PI
#undef PI
#endif
#define PI  3.14159265358979323846264338327950
#define PIh 1.57079632679489661923132169163975
#define PIq 0.785398163397448309615660845819876


/* snap effect for grid coordinates */
#define dequaleps 1e-10
#define dless(a,b) ((a)<(b)-dequaleps)
#define dequal(a,b) (!(dless(a,b)||dless(b,a)))
#define dgreater(a,b) ((a)>(b)+dequaleps)

/* approx. <= and >= with dequaleps tolerance */
#define dlesseq(a,b) ( (a)<(b)+dequaleps )
#define dgreatereq(a,b) ( (a)>(b)-dequaleps )

#define signum(v) ((v) > 0.0 ? (1.0) : ((v) < 0.0 ? (-1.0) : (0.0)))


/* utilities.c */
void errorexit(char *file, int line, char *s);
void errorexits(char *file, int line, char *s, char *t);
void errorexiti(char *file, int line, char *s, int i);
#define errorexit(s) errorexit(__FILE__, __LINE__, (s))
#define errorexits(s,t) errorexits(__FILE__, __LINE__, (s), (t))
#define errorexiti(s,i) errorexiti(__FILE__, __LINE__, (s), (i))

void yo(void);
void Yo(double x);
void prdivider(int n);
double min2(double x, double y);
double min3(double x, double y, double z);
double max2(double x, double y);
double max3(double x, double y, double z);
double min_in_1d_array(double *f, int n, int *imin);
double max_in_1d_array(double *f, int n, int *imax);
double min2_in_1d_array(double *f0, int n0, double *f1, int n1, 
                        int *ai, int *imin);
double max2_in_1d_array(double *f0, int n0, double *f1, int n1, 
                        int *ai, int *imax);
double min3_in_1d_array(double *f0, int n0, double *f1, int n1, double *f2, int n2,
                        int *ai, int *imin);
double max3_in_1d_array(double *f0, int n0, double *f1, int n1, double *f2, int n2,
                        int *ai, int *imax);
int copy_file_into_dir(char *fname, char *dir);
int system2(char *s1, char *s2);
int system3(char *s1, char *s2, char *s3);
int system_emu(const char *command);
int lock_curr_til_EOF(FILE *out);
int unlock_curr_til_EOF(FILE *out);
int construct_argv(char *str, char ***argv);
double *dmalloc(int n);
int *imalloc(int n);
char *cmalloc(int n);
void *pmalloc(int n);
