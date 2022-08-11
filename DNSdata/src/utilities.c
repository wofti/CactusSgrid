/* utilities.c */
/* Wolfgang Tichy 2016*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>

/* for POSIX.1-2001 mkdir, opendir, fork, wait functions */
#include <unistd.h>     /* for fork */
#include <sys/stat.h>
#include <sys/types.h>  /* for pid_t */
#include <sys/wait.h>   /* for wait */
#include <dirent.h>     /* for opendir */

#include "utilities.h"


/* use opendir to scan through dir and remove the entire dir */
int th_remove_dir(char *which_dir)
{
  DIR           *d;
  struct dirent *dir;
  char file[1000];
  
  d = opendir(which_dir);
  if(d)
  {
    while ((dir = readdir(d)) != NULL)
    {
      /* exclude . and .. directories */
      if( strcmp( dir->d_name, "." ) == 0 || strcmp( dir->d_name, ".." ) == 0)
        continue;

      snprintf(file, 999, "./%s/%s", which_dir, dir->d_name);
      //printf("*"); //print * for every deleted file

      if(opendir(file)!=NULL) th_remove_dir(file);
      else
      {
        if(remove(file) != 0)
        {
          //printf("\n%s\n", file);
          //perror("Remove failed");
          closedir(d);
          return -2;
        }
      }
    } /* end of while loop */

    closedir(d);

    /* delete directory */
    if(remove(which_dir) != 0)
    {
      //printf("%s\n", which_dir);
      //perror("Remove failed");
      return -1;
    }
  }
  return 0;
}

/* ugh, but how universal are those built in functions? */
int th_system2(char *s1, char *s2) 
{
  return th_system3(s1, s2, "");
}
int th_system3(char *s1, char *s2, char *s3) 
{
  char command[10000];
  int status = 0;

  /* check for special cases where we can use c-functions */
  if( strcmp(s1,"mv")==0 || strcmp(s1,"mv -f")==0 ) /* use rename */
  {
    sprintf(command, "rename(\"%s\", \"%s\");", s2, s3);
    status = rename(s2, s3);
    printf("BNS_utilities_th_system3: ANSI C call: %s\n", command);
  }
  else if( strcmp(s1,"rm -rf")==0 ) /* use remove */
  {
    if(strlen(s2)>0)
    {
      status = th_remove_dir(s2);
      printf("BNS_utilities_th_system3: remove_dir(\"%s\");\n", s2);
    }
    if(strlen(s3)>0)
    {
      status = th_remove_dir(s3);
      printf("BNS_utilities_th_system3: remove_dir(\"%s\");\n", s3);
    }
  }
  else if( strcmp(s1,"mkdir")== 0 ) /* use POSIX.1-2001 mkdir function */
  {
    sprintf(command, "mkdir(\"%s\", S_IRWXU | S_IRWXG);", s2);
    status = mkdir(s2, S_IRWXU | S_IRWXG);
    printf("BNS_utilities_th_system3: POSIX.1-2001 call: %s\n", command);
  }
  else /* use system */
  { 
    sprintf(command, "%s %s %s", s1, s2, s3);
    status = system(command);
    printf("BNS_utilities_th_system3: System call: %s\n", command);
  }
  
  if(status!=0) printf("BNS_utilities_th_system3: -> WARNING: Return value = %d\n", status);
  return status;
}

/* construct an argv array from a string and return number of args */
/* NOTE: str is modified and used as mem for argv! */
int th_construct_argv(char *str, char ***argv)
{
  char *str1, *token, *saveptr;
  int count;

  *argv = NULL;
  for(count=0, str1=str; ; count++, str1=NULL)
  {
    *argv = (char **) realloc(*argv, sizeof(char *)*(count+1));
    token = strtok_r(str1, " ", &saveptr);
    //printf("token=%p:%s\n", token,token);
    (*argv)[count] = token;
    if(token == NULL) break;
  }
  //printf("saveptr=%p:%s\n", saveptr,saveptr);
  return count;
}

/* run a command, without a shell */
int th_system_emu(const char *command)
{
  char *com = strdup(command); /* duplicate since construct_argv modifies its args */
  int ret, status;

    printf("BNS_utilities_th_system_emu: th_system_emu: running command:\n%s\n", command);
    char **argv;
    construct_argv(com, &argv);
    ret = libsgrid_main(6, argv);

  status = ret;
  if(status!=0) printf("BNS_utilities_th_system_emu: -> WARNING: Return value = %d\n", status);
  free(com);
  return status;
}

/* Lock a file from current file position to the end. The lock will be
   released when the file is closed.
   fd is a file descriptor open for writing. */
int th_lock_curr_til_EOF(FILE *out)
{
  int fd = fileno(out); /* get file dscriptor */
  if(fd==-1) return fd; /* return -1 on error */
  return lockf(fd, F_LOCK, 0);
}


/* preliminary ... */
double *th_dmalloc(int n)
{
  double *p = (double *) malloc(sizeof(double) * n);
  
  if (!p) th_errorexiti("out of memory (%d double)", n);
  return p;
}

int *th_imalloc(int n)
{
  int *p = (int *) malloc(sizeof(int) * n);
  
  if (!p) th_errorexiti("out of memory (%d int)", n);
  return p;
}

char *th_cmalloc(int n)
{
  char *p = (char *) malloc(sizeof(char) * n);
  
  if (!p) th_errorexiti("out of memory (%d char)", n);
  return p;
}

void *th_pmalloc(int n)
{
  void *p = malloc(sizeof(void *) * n);
  
  if (!p) th_errorexiti("out of memory (%d void *)", n);
  return p;
}





/* the one function every program should have */
/* note that sgrid_main.h defines a macro so that the user does not have
   to specify __FILE__ and __LINE__ for location where the error occured
*/

#undef th_errorexits
#undef th_errorexiti

void th_errorexits(char *file, int line, char *s, char *t)
{
  fflush(stdout);
  fprintf(stderr, "Error: ");
  fprintf(stderr, s, t);
  fprintf(stderr, "  (%s, line %d)\n", file, line);
  fflush(stderr);
  //sgrid_MPI_Finalize();
  exit(1);
}

void th_errorexiti(char *file, int line, char *s, int i)
{
  fflush(stdout);
  fprintf(stderr, "Error: ");
  fprintf(stderr, s, i);
  fprintf(stderr, "  (%s, line %d)\n", file, line);
  fflush(stderr);
  //sgrid_MPI_Finalize();
  exit(1);
}

/* do not write functions beyond this line because the undef/define 
   method for the errorexit functions means that they should go last */
