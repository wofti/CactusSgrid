/* utilities.h */
/* Wolfgang Tichy and Michal Pirog 8, 2022 */

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


int DNS2_run(const char *command);
int DNS2_lock_curr_til_EOF(FILE *out);
int DNS2_construct_argv(char *str, char ***argv);
double *DNS2_dmalloc(int n);
int *DNS2_imalloc(int n);
char *DNS2_cmalloc(int n);
void *DNS2_pmalloc(int n);
