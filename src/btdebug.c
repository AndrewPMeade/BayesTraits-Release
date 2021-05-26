#include "btdebug.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
//#include <stdarg.h>

#define MAX_INDENTATION 20
#define INDENT 4
#define MAX_SPACES (INDENT + MAX_DEPTH*INDENT)
#define MAX_DEPTH 15
#define MAX_DEBUG 15
#define MAX_FNAME 20

char debug_buffer[256];

typedef struct {
  double miliseconds;  // start time
  char fname[MAX_FNAME];
} stack_info;

typedef struct debug_info {
  char fname[MAX_FNAME];
  double total_time;
  unsigned int number_runs;
} debug_info;

static char spaces[MAX_SPACES];
static char* indentation;

static unsigned int depth = 0;  // [0,+inf]
static stack_info debugStack[MAX_DEPTH];

/* ********* debug info stack manipulation: start ************** */
int stack_isfull() {
  return (depth >= MAX_DEPTH);
}

stack_info* pushFrame(const char* f,double ms) {
  int l;
  stack_info* newframe = NULL;
  if (depth < MAX_DEPTH) {  // add function info to stack
    if (strlen(f) < MAX_FNAME) {
      l = strlen(f);
    } else {
      l = MAX_FNAME-1;
    }
    newframe = &debugStack[depth];
    debugStack[depth].miliseconds = ms;
    memcpy(debugStack[depth].fname, f, l);
    debugStack[depth].fname[l] = '\0';
    debugStack[depth].miliseconds = ms;
    indentation -= INDENT;
    depth++;
  }
  return newframe;
}

// returns the top of the stack (a valid index). Decrements stack size.
stack_info*  popFrame(const char* f) {
  // need to check if function nam eis valid

  if (depth <= 0)  // error
    return NULL;
  if (depth >= MAX_DEPTH) {
    depth--;
    return NULL;
  } else {
    indentation += INDENT;
    --depth;
    return &debugStack[depth];
  }
}
/* ********* debug info stack manipulation: end ************** */

void btdebug_init() {
  memset(spaces,' ',MAX_SPACES);
  indentation = (char*)&spaces[MAX_SPACES-1]; // last index
  *indentation = '\0';
  debug_buffer[0] = '\0';
  depth = 0;
}

//void btdebug_print(const char *fmt, ...)
//{
//    va_list args;
//
//    va_start(args, fmt);
//    printf("%s",indentation);
//    printf(fmt, args);
//    printf("\n");
//   va_end(args);
//

void btdebug_print(const char *msg)
{
    printf("%s",indentation);
    printf("%s\n",msg);
}

void btdebug_enter(const char* f) {
  stack_info* frame;
  double ms = btdebug_ms();
  printf("%s",indentation);
  printf("--> %s\n",f);
  frame = pushFrame(f,ms);
}

void btdebug_exit(const char* f) {
  stack_info* frame;
  double now;
  double elapsed = -1;
  now = btdebug_ms();
  frame = popFrame(f);
  if (strcmp(frame->fname,f) != 0) {
    printf("Error at debug_exit: wrong function name\n");
    exit(0);
  }
    
  if (frame) {
    elapsed = now - frame->miliseconds;
  }
  printf("%s",indentation);
  printf("<-- %s: elapsed %lf\n",f,elapsed);
}

double btdebug_ms() {
	//return (double)clock();
	return (double)clock()/(CLOCKS_PER_SEC/1000.0);
}
