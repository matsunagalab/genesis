
#include <stdio.h>
#include <setjmp.h>
#include <string.h>
#include "../config.h"

/*==========================================================================
 *  Library-mode error guard (setjmp/longjmp)
 *
 *  GENESIS' error_msg() aborts the whole process with exit(1). That is fine
 *  for the CLI tools and the MD engines, but it is fatal for the Python
 *  interface: a bad input file would tear down the host Python process
 *  instead of raising a catchable exception.
 *
 *  These helpers let a bind(C) wrapper run its body under a setjmp guard.
 *  While a guard is armed, error_msg() records the pending message/code
 *  (in messages_mod) and calls fi_error_signal(), which longjmps back to
 *  fi_error_guard_run_ctx() instead of exiting. Outside a guard (no armed
 *  frame) error_msg() keeps its original exit(1) behaviour, so the
 *  CLI/atdyn/spdyn binaries are completely unaffected.
 *
 *  The body receives an opaque context pointer so that it can be an ordinary
 *  (module) procedure. Do NOT pass c_funloc() of an internal procedure:
 *  gfortran implements that with a trampoline on the stack, which marks the
 *  shared library as requiring an executable stack, and dlopen() rejects
 *  such libraries on glibc >= 2.41.
 *
 *  Guards nest correctly: every frame keeps its own jmp_buf and the innermost
 *  armed one is the longjmp target. Not thread-safe (the interface as a
 *  whole is not).
 *
 *  These functions live here (rather than in a dedicated file) so that they
 *  are compiled into lib.a and therefore visible to every binary that links
 *  it, without touching the autotools-generated build files.
 *========================================================================*/

static jmp_buf *volatile fi_guard_cur = NULL;   /* innermost armed guard */

/* Is a guard currently armed on this call stack? (queried by error_msg) */
int fi_error_is_armed(void)
{
  return fi_guard_cur != NULL;
}

/* Abort the current guarded region. Does nothing (returns) when no guard is
 * armed, letting the caller fall through to its normal exit(1) path. */
void fi_error_signal(void)
{
  jmp_buf *target = fi_guard_cur;
  if (target == NULL)
    return;
  longjmp(*target, 1);
}

/* Run body(ctx) under a setjmp guard. Returns 0 when body() completes
 * normally and 1 when body() triggered fi_error_signal(). */
int fi_error_guard_run_ctx(void (*body)(void *), void *ctx)
{
  jmp_buf        buf;
  jmp_buf *const prev = fi_guard_cur;
  int            rc;

  fi_guard_cur = &buf;
  if (setjmp(buf) == 0) {
    body(ctx);
    rc = 0;
  } else {
    rc = 1;
  }
  fi_guard_cur = prev;
  return rc;
}

#define MaxFile  10

static FILE* Files [MaxFile] = { 0 };

void fd_open__(char* filename, const int* read_only, int* unit_no,
			  int fnlen)
{
  int i;
  for(i = fnlen-1; i >= 0; --i)
	if (filename[i] == ' ')
	  filename[i] = 0;
	else
	  break;
	  
  for(i = 0; i < MaxFile; ++i)
	if (Files[i] == 0)
	  {
		*unit_no = i;
		Files[i] = fopen(filename, (*read_only) == 1 ? "r" : "w");
		break;
	  }

  if (i == MaxFile)
	*unit_no = -1;
}

void fd_close__(int* unit_no)
{
  if (*unit_no < 0 || *unit_no >= MaxFile)
	return ;
  fclose(Files[*unit_no]);
  Files[*unit_no] = 0;
}

void fd_flen__(int* unit_no, long long* flen)
{
#ifdef KCOMP
  FILE* fp = Files[*unit_no];
  fseeko64(fp, 0, SEEK_END);
  *flen = ftello64(fp);
  fseeko64(fp, 0, SEEK_SET);
#else
  fprintf(stderr, "fd_flen__ doesn't work on the machine\n");
  *flen = -1;
#endif
}

void fd_read__(int* unit_no, long long *pos, unsigned char* b, long long *blen)
{
#ifdef KCOMP
  FILE* fp = Files[*unit_no];
  fseeko64(fp, (*pos)-1, SEEK_SET);
  fread(b, *blen, 1, fp);
#else
  fprintf(stderr, "fd_read__ doesn't work on the machine\n");
  *blen = -1;
#endif
}

void fd_write__(int* unit_no, unsigned char* b, long long *blen)
{
  FILE* fp = Files[*unit_no];
  fwrite(b, *blen, 1, fp);
}
