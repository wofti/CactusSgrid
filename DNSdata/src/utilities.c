/* utilities.c */
/* Wolfgang Tichy and Michal Pirog 8, 2022 */

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
int DNS2_remove_dir(char *which_dir)
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

      if(opendir(file)!=NULL) DNS2_remove_dir(file);
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
//    printf("DNS2_remove_dir: done_06 \n");
  return 0;
}

/* construct an argv array from a string and return number of args */
/* NOTE: str is modified and used as mem for argv! */
int DNS2_construct_argv(char *str, char ***argv)
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
//  printf("DNS2_construct_argv: done_07 \n");
  return count;
}

/* run a command, without a shell */
int DNS2_run(const char *command)
{
  char *com = strdup(command); /* duplicate since construct_argv modifies its args */
  int ret, status;

  printf("DNS2_run: running command:\n%s\n", command);
  char **argv;
  DNS2_construct_argv(com, &argv);
  ret = libsgrid_main(6, argv);

  status = ret;
  if(status!=0) printf("DNS2_run: -> WARNING: Return value = %d\n", status);
  free(com);
    
//  printf("DNS2_run: done_08 \n");
  return status;
}

/* preliminary ... */
double *DNS2_dmalloc(int n)
{
  double *p = (double *) malloc(sizeof(double) * n);
  
  if (!p) DNS2_errorexiti("out of memory (%d double)", n);

//  printf("*DNS2_dmalloc: done_10a \n");
  return p;
}

int *DNS2_imalloc(int n)
{
  int *p = (int *) malloc(sizeof(int) * n);
  
  if (!p) DNS2_errorexiti("out of memory (%d int)", n);
  
//  printf("*DNS2_imalloc: done_10b \n");
  return p;
}

char *DNS2_cmalloc(int n)
{
  char *p = (char *) malloc(sizeof(char) * n);
  
  if (!p) DNS2_errorexiti("out of memory (%d char)", n);
  
//  printf("*DNS2_cmalloc: done_10c \n");
  return p;
}

void *DNS2_pmalloc(int n)
{
  void *p = malloc(sizeof(void *) * n);
  
  if (!p) DNS2_errorexiti("out of memory (%d void *)", n);
  
//  printf("*DNS2_pmalloc: done_10d \n");
  return p;
}

/* the one function every program should have */
/* note that sgrid_main.h defines a macro so that the user does not have
   to specify __FILE__ and __LINE__ for location where the error occured */

#undef DNS2_errorexits
#undef DNS2_errorexiti

void DNS2_errorexits(char *file, int line, char *s, char *t)
{
  fflush(stdout);
  fprintf(stderr, "Error: ");
  fprintf(stderr, s, t);
  fprintf(stderr, "  (%s, line %d)\n", file, line);
  fflush(stderr);
//  printf("DNS2_errorexits: done_11 \n");
  exit(1);
}

void DNS2_errorexiti(char *file, int line, char *s, int i)
{
  fflush(stdout);
  fprintf(stderr, "Error: ");
  fprintf(stderr, s, i);
  fprintf(stderr, "  (%s, line %d)\n", file, line);
  fflush(stderr);
//  printf("DNS2_errorexiti: done_12 \n");
  exit(1);
}

/* do not write functions beyond this line because the undef/define 
   method for the errorexit functions means that they should go last */
