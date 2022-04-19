/*----------------------------------------------------------------------
  File    : seqwog.c
  Contents: find frequent sequences without gaps
  Author  : Christian Borgelt
  History : 2010.08.06 file created from relim.c
            2010.08.10 bug in sorting direction fixed (ascending)
            2010.08.22 adapted to modified modules tabread and tract
            2010.10.15 adapted to modified interface of module report
            2010.11.24 adapted to modified error reporting (tract)
            2010.12.11 adapted to a generic error reporting function
            2011.03.20 optional integer transaction weights added
            2011.03.24 bug in recurse() fixed (single item tails)
            2011.07.08 adapted to modified function tbg_recode()
            2011.08.28 output of sequence counters per size added
            2012.05.30 bug in closedness check w.r.t. head fixed
            2012.06.27 maximal sequence mining added (option -tm)
            2012.06.28 bug w.r.t. empty sequence fixed (closed/maximal)
            2013.04.03 adapted to type changes in module tract
            2013.10.08 generating simple sequence rules added
            2013.10.15 checks of result of isr_iset() added
            2013.10.18 optional pattern spectrum collection added
            2014.05.12 option -F# added (support border for filtering)
            2014.08.26 adapted to modified item set reporter interface
            2014.10.24 changed from LGPL license to MIT license
            2015.06.04 bug in reporting of long sequences fixed
            2015.06.05 bug in filtering closed/maximal sequences fixed
----------------------------------------------------------------------*/
#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <limits.h>
#include <float.h>
#include <math.h>
#include <time.h>
#include <assert.h>
#ifndef TA_READ
#define TA_READ
#endif
#include "tract.h"
#ifndef ISR_PATSPEC
#define ISR_PATSPEC
#endif
#ifndef PSP_REPORT
#define PSP_REPORT
#endif
#include "report.h"
#include "error.h"
#ifdef STORAGE
#include "storage.h"
#endif

/*----------------------------------------------------------------------
  Preprocessor Definitions
----------------------------------------------------------------------*/
#define PRGNAME     "seqwog"
#define DESCRIPTION "find frequent sequences without gaps"
#define VERSION     "version 3.16 (2016.10.15)        " \
                    "(c) 2010-2016   Christian Borgelt"

/* --- operation modes --- */
#define SEQ_ALLOCC   0x10       /* count all occurrences */

/* --- error codes --- */
/* error codes   0 to  -4 defined in tract.h */
#define E_STDIN      (-5)       /* double assignment of stdin */
#define E_OPTION     (-6)       /* unknown option */
#define E_OPTARG     (-7)       /* missing option argument */
#define E_ARGCNT     (-8)       /* too few/many arguments */
#define E_TARGET     (-9)       /* invalid target type */
#define E_SIZE      (-10)       /* invalid sequence length */
#define E_SUPPORT   (-11)       /* invalid minimum sequence support */
#define E_CONF      (-12)       /* invalid confidence */
/* error codes -15 to -25 defined in tract.h */

#ifndef QUIET                   /* if not quiet version, */
#define MSG         fprintf     /* print messages */
#define CLOCK(t)    ((t) = clock())
#else                           /* if quiet version, */
#define MSG(...)    ((void)0)   /* suppress messages */
#define CLOCK(t)    ((void)0)
#endif

#define SEC_SINCE(t)  ((double)(clock()-(t)) /(double)CLOCKS_PER_SEC)

/*----------------------------------------------------------------------
  Type Definitions
----------------------------------------------------------------------*/
typedef struct {                /* --- sequence occurrence --- */
  TID        tid;               /* transaction identifier */
  SUPP       supp;              /* support (number of occurrences) */
  ITEM       head;              /* item before current pattern */
  const ITEM *tail;             /* pattern occurrence tail */
} SEQOCC;                       /* (sequence occurrence) */

typedef struct {                /* --- recursion data --- */
  ITEMBASE   *base;             /* underlying item base */
  int        target;            /* target type (e.g. closed/maximal) */
  int        mode;              /* operation mode (e.g. pruning) */
  SUPP       smin;              /* minimum support of a sequence */
  SUPP       body;              /* minimal support of a rule */
  double     conf;              /* minimum confidence of a rule */
  TID        curr;              /* marker for support counting */
  ITEM       zmax;              /* maximum length of a sequence */
  TID        *marks;            /* transaction markers */
  SUPP       *buf;              /* buffer for maximality check */
  SUPP       *sall;             /* item support for all occurrences */
  ISREPORT   *report;           /* item set/sequence reporter */
} RECDATA;                      /* (recursion data) */

/*----------------------------------------------------------------------
  Constants
----------------------------------------------------------------------*/
#if !defined QUIET && defined SWOG_MAIN
/* --- error messages --- */
static const char *errmsgs[] = {
  /* E_NONE      0 */  "no error",
  /* E_NOMEM    -1 */  "not enough memory",
  /* E_FOPEN    -2 */  "cannot open file %s",
  /* E_FREAD    -3 */  "read error on file %s",
  /* E_FWRITE   -4 */  "write error on file %s",
  /* E_STDIN    -5 */  "double assignment of standard input",
  /* E_OPTION   -6 */  "unknown option -%c",
  /* E_OPTARG   -7 */  "missing option argument",
  /* E_ARGCNT   -8 */  "wrong number of arguments",
  /* E_TARGET   -9 */  "invalid target type '%c'",
  /* E_SIZE    -10 */  "invalid sequence length %"ITEM_FMT,
  /* E_SUPPORT -11 */  "invalid minimum support %g",
  /* E_CONF    -12 */  "invalid minimum confidence %g",
  /*    -13 to -14 */  NULL, NULL,
  /* E_NOITEMS -15 */  "no (frequent) items found",
  /*           -16 */  "unknown error"
};
#endif

/*----------------------------------------------------------------------
  Global Variables
----------------------------------------------------------------------*/
#ifdef SWOG_MAIN
#ifndef QUIET
static CCHAR    *prgname;       /* program name for error messages */
#endif
static TABREAD  *tread  = NULL; /* table/transaction reader */
static ITEMBASE *ibase  = NULL; /* item base */
static TABAG    *tabag  = NULL; /* transaction bag/multiset */
static ISREPORT *report = NULL; /* item set/sequence reporter */
static TABWRITE *twrite = NULL; /* table writer for pattern spectrum */
static double   *border = NULL; /* support border for filtering */
#endif

/*----------------------------------------------------------------------
  Auxiliary Functions for Debugging
----------------------------------------------------------------------*/
#ifndef NDEBUG

static void indent (int k)
{ while (--k >= 0) printf("   "); }

/*--------------------------------------------------------------------*/

static void show (const SEQOCC *occs, TID n, int ind)
{                               /* --- show an occurrence array */
  TID       i;                  /* loop variable */
  ITEM      h;                  /* head item */
  const int *s;                 /* to traverse the tail items */

  assert(occs && (n >= 0));     /* check the function arguments */
  for (i = 0; i < n; i++) {     /* traverse the array elements */
    indent(ind);                /* indent the output line */
    printf("%"TID_FMT "|", occs[i].tid);  /* print transaction id, */
    printf("%"SUPP_FMT"|", occs[i].supp); /* support and head item */
    if ((h = occs[i].head) < 0) printf("*/*:");
    else printf("%s/%"ITEM_FMT":", ib_name(ibase, h), h);
    for (s = occs[i].tail; *s >= 0; s++)
      printf(" %s/%"ITEM_FMT, ib_name(ibase, *s), *s);
    printf("\n");               /* traverse and print tail items */
  }                             /* finally terminate the output line */
  printf("%"TID_FMT" occurrence(s)\n", n);
}  /* show() */                 /* print the number of occurrences */

#endif
/*----------------------------------------------------------------------
  Sequence Mining without Gaps
----------------------------------------------------------------------*/

static int cmp (const void *p, const void *q)
{                               /* --- compare two sequence occs. */
  const ITEM *a = ((SEQOCC*)p)->tail;
  const ITEM *b = ((SEQOCC*)q)->tail;
  for ( ; 1; a++, b++) {        /* lexicographic comparison loop */
    if (*a < *b) return -1;     /* compare corresponding items */
    if (*a > *b) return +1;     /* and if one is greater, abort */
    if (*a <  0) return  0;     /* otherwise check for the sentinel */
  }                             /* and abort if it is reached */
}  /* cmp() */

/*--------------------------------------------------------------------*/
#if 1

static int maximal (const SEQOCC *occs, TID n, RECDATA *rd)
{                               /* --- check for a maximal sequence */
  ITEM k;                       /* buffer for an item */
  TID  i, m, t;                 /* loop variable, buffer */
  SUPP s, max = 0;              /* maximum support of a head item */

  assert(occs && (n > 0) && rd);/* check the function arguments */
  m = (rd->mode & SEQ_ALLOCC) ? 0 : ++rd->curr;
  for (i = 0; i < n; i++) {     /* traverse the sequence occurrences */
    if ((k = occs[i].head) < 0) continue;
    t = occs[i].tid;            /* get trans. info. for occurrence */
    if (rd->marks[t] != rd->curr) rd->buf[k] += occs[i].supp;
    rd->marks[t] = m;           /* update the head item support */
  }                             /* and the transaction marker */
  for (i = 0; i < n; i++) {     /* traverse the sequence occurrences */
    if ((k = occs[i].head) < 0) continue;
    if ((s = rd->buf[k]) > max) max = s;
    rd->buf[k] = 0;             /* find the maximum head item support */
  }                             /* and clear the support counters */
  return (max < rd->smin);      /* return whether sequence is maximal */
}  /* maximal() */

/*--------------------------------------------------------------------*/
#else

static int maximal (const SEQOCC *occs, TID n, RECDATA *rd)
{                               /* --- check for a maximal sequence */
  int  r;                       /* return value */
  ITEM k;                       /* buffer for an item */
  TID  i, m, t;                 /* loop variable, buffers */

  assert(occs && (n > 0) && rd);/* check the function arguments */
  m = (rd->mode & SEQ_ALLOCC) ? 0 : ++rd->curr;
  for (i = 0; i < n; i++) {     /* traverse the sequence occurrences */
    if ((k = occs[i].head) < 0) continue;
    t = occs[i].tid;            /* get trans. info. for occurrence */
    if ((rd->marks[t] != rd->curr)
    && ((rd->buf[k] += occs[i].supp) >= rd->smin)) break;
    rd->marks[t] = m;           /* update the head item support */
  }                             /* and the transaction marker */
  r = (i >= n);                 /* clear the support counters */
  for (i = (r) ? n : i+1; --i >= 0; )
    if ((k = occs[i].head) >= 0) rd->buf[k] = 0;
  return r;                     /* return whether sequence is maximal */
}  /* maximal() */

#endif
/*--------------------------------------------------------------------*/

static SUPP recurse (const SEQOCC *occs, TID n, ITEM len, RECDATA *rd)
{                               /* --- recursive part of the search */
  ITEM item;                    /* buffer for an item */
  TID  i, m, t;                 /* loop variables, buffers */
  SUPP s, supp, smax = 0;       /* support of current sequence */
  ITEM head;                    /* head of first occurrence */
  ITEM clsd;                    /* flag for closed w.r.t. head */

  assert(occs                   /* check the function arguments */
  &&    (n > 0) && (len > 0) && rd);
  for (i = 0; i < n; i++)       /* find next tail item (extension) */
    if (occs[i].tail[len] >= 0) break;
  while ((n -= i) > 0) {        /* traverse the pattern occurrences */
    occs += i;                  /* skip the already processed occs. */
    m = (rd->mode & SEQ_ALLOCC) ? 0 : ++rd->curr;
    rd->marks[occs[0].tid] = m; /* update the transaction marker */
    head = occs[0].head;        /* get the first head item */
    clsd = head < 0;            /* and check whether it exists */
    item = occs[0].tail[len];   /* get the next tail item (extension) */
    supp = occs[i = 0].supp;    /* traverse pattern occurrences */
    while ((++i < n)            /* with the same tail item */
    &&     (occs[i].tail[len] == item)) {
      clsd |= head ^ occs[i].head;
      t = occs[i].tid;          /* get trans. info. for occurrence */
      if (rd->marks[t] != rd->curr) supp += occs[i].supp;
      rd->marks[t] = m;         /* update the sequence support */
    }                           /* and the transaction marker */
    if ((supp < rd->smin)       /* if infrequent or not closed */
    ||  ((rd->target & (ISR_CLOSED|ISR_MAXIMAL)) && !clsd))
      continue;                 /* skip the extension item */
    if (supp > smax) smax = supp;
    s = (len >= rd->zmax) ? 0   /* recusively find freq. sequences */
      : recurse(occs, i, len+1, rd);
    if      (rd->target & ISR_MAXIMAL)
      m = (s < rd->smin) && maximal(occs, i, rd);
    else if (rd->target & ISR_CLOSED)
      m = (s < supp);           /* check whether to report */
    else m = -1;                /* the current sequence */
    if (!m) continue;           /* if not to report, continue loop */
    if (isr_iset(rd->report, occs[0].tail, len+1, supp, 0, 0) < 0)
      return -1;                /* report the current sequence */
  }
  return smax;                  /* return maximal extension support */
}  /* recurse() */

/*--------------------------------------------------------------------*/

static int rules (const SEQOCC *occs, TID n,
                  ITEM len, SUPP body, RECDATA *rd)
{                               /* --- recursive part for rules */
  ITEM   item;                  /* buffer for an item */
  TID    i, k;                  /* loop variables */
  SUPP   supp = 0;              /* support of current sequence */

  assert(occs                   /* check the function arguments */
  &&    (n > 0) && (len > 0) && rd);
  for (i = 0; i < n; i++)       /* find next tail item (extension) */
    if (occs[i].tail[len] >= 0) break;
  for (k = i; k < n; k += i) {  /* traverse the pattern occurrences */
    occs += i;                  /* skip the already processed occs. */
    item = occs[0].tail[len];   /* get the next tail item (extension) */
    supp = occs[i = 0].supp;    /* traverse pattern occurrences */
    while ((++i < n-k)          /* with the same tail item */
    &&     (occs[i].tail[len] == item))
      supp += occs[i].supp;     /* update the sequence support */
    if  (supp <  rd->smin)      /* if extension item is infrequent, */
      continue;                 /* skip the extension item */
    if ((supp >= rd->body)      /* if there could be larger rules */
    &&  (len  <  rd->zmax)      /* (enough support, below max. size) */
    &&  (rules(occs, i, len+1, supp, rd) < 0))
      return -1;                /* recurse to search for larger rules */
    if ((body >= rd->body)      /* if sufficient support and conf. */
    &&  ((double)supp >= rd->conf*(double)body)
    &&  (isr_seqrule(rd->report, occs[0].tail, len+1,
                     supp, body, rd->sall[item], 0) < 0))
      return -1;                /* report the current rule */
  }
  return 0;                     /* return 'ok' */
}  /* rules() */

/*--------------------------------------------------------------------*/

static int seqwog (TABAG *tabag, int target, SUPP smin, SUPP body,
                   double conf, int mode, ISREPORT *report)
{                               /* --- search for frequent sequences */
  int        e = 0;             /* error status */
  ITEM       i, k;              /* loop variable, number of items */
  TID        j, n, m, r;        /* loop variable, number of trans. */
  size_t     x;                 /* size of transaction counter array */
  ITEM       head;              /* head of first occurrence */
  ITEM       clsd;              /* flag for closed w.r.t. head */
  SUPP       supp, w, f;        /* weight/support buffer */
  SUPP       max;               /* maximum support of an item */
  TRACT      *t;                /* to traverse the transactions */
  TID        *c;                /* item occurrence counters */
  const ITEM *s, *d;            /* to traverse the items */
  SEQOCC     *occs;             /* to traverse the occurrence array */
  RECDATA    rd;                /* recursion data */

  assert(tabag && report);      /* check the function arguments */
  rd.target = target;           /* store target type and search mode */
  rd.mode   = mode;             /* check and adapt minimum support */
  rd.body   = (body > 0) ? body : 1;
  if ((tbg_wgt(tabag) < rd.body) && !(mode & SEQ_ALLOCC))
    return 0;                   /* check the total transaction weight */
  rd.smin   = (smin > 0) ? smin : 1;
  rd.conf   = conf;             /* note the minimum confidence */
  rd.base   = tbg_base(tabag);  /* initialize the recursion data */
  rd.report = report;           /* (note item base and reporter) */
  rd.zmax   = isr_zmax(report); /* get maximum seq. length to check */
  if ((rd.zmax < ITEM_MAX) && (target & (ISR_CLOSED|ISR_MAXIMAL)))
    rd.zmax++;                  /* one more for closed/maximal */
  k = tbg_itemcnt(tabag);       /* get the number of items */
  n = tbg_cnt(tabag);           /* and the number of transactions */
  if (n <= 0) {                 /* if there are no transactions */
    if (!(mode & SEQ_ALLOCC) && !(target & ISR_RULES))
      return isr_report(report);/* report the empty sequence */
  }                             /* and abort the search */
  x = ((size_t)k > (size_t)n) ? (size_t)k : (size_t)n;
  rd.marks = (TID*)malloc(x *sizeof(TID) +(size_t)(k+k) *sizeof(SUPP));
  if (!rd.marks) return -1;     /* allocate counter and marker array */
  memset(c = rd.marks, 0, (size_t)k *sizeof(TID));
  rd.buf  = (SUPP*)(c+x);       /* note buffers for closedness test */
  rd.sall = rd.buf +k;          /* and item support for all occs. */
  memset(rd.sall, 0, (size_t)k *sizeof(SUPP));
  for (j = 0; j < n; j++) {     /* traverse the transactions */
    t = tbg_tract(tabag, j); w = ta_wgt(t);
    for (s = ta_items(t); *s >= 0; s++) {
      c[*s]++;                  /* count the item occurrences */
      rd.sall[*s] += w;         /* compute the item support */
    }                           /* for all occurrences */
  }                             /* (also needed for rule heads) */
  for (m = 0, i = 0; i < k; i++)/* traverse the items and find the */
    if (c[i] > m) m = c[i];     /* maximum number of occurrences, */
  memset(c, 0, (size_t)n *sizeof(TID));    /* clear occ. counters */
  occs = (SEQOCC*)malloc((size_t)m *sizeof(SEQOCC));
  if (!occs) return -1;         /* create an occurrence array */
  rd.curr = 1; max = 0;         /* init. the transaction marker */
  while (--k >= 0) {            /* traverse items (sequence starts) */
    f = (mode & SEQ_ALLOCC) ? rd.sall[k] : ib_getfrq(rd.base, k);
    if (f < rd.smin) continue;  /* skip infrequent items */
    head = clsd = 0;            /* init. head and closed flag */
    m = (mode & SEQ_ALLOCC) ? 0 : ++rd.curr;
    for (supp = 0, j = r = 0; j < n; j++) {
      t = tbg_tract(tabag, j);  /* traverse the transactions */
      w = ta_wgt(t);            /* increment the current marker */
      for (s = d = ta_items(t); *s > TA_END; s++) {
        if (*s != k) continue;  /* traverse the item occurrences */
        if (rd.marks[j] != rd.curr) supp += w;
        rd.marks[j]  = m;       /* compute support and update marker */
        occs[r].tid  = j;       /* store the transaction identifier, */
        occs[r].supp = w;       /* the support and the head item */
        occs[r].head = (s > d) ? s[-1] : -1;
        if (r <= 0)  head = occs[r].head;
        else clsd |= head ^ occs[r].head;
        occs[r++].tail = s;     /* combine heads to check whether */
      }                         /* the item is closed w.r.t. head */
    }                           /* finally note the tail sequence */
    if (supp > max) max = supp; /* find maximal support of an item */
    if ((supp < rd.body)        /* if the item is infrequent */
    || ((target & ISR_CLOSED)   /* or not closed (if requested) */
    &&  (head >= 0) && !clsd))  /* (check for an existing head), */
      continue;                 /* skip the item */
    qsort(occs, (size_t)r, sizeof(SEQOCC), cmp);
    if (target & ISR_RULES) {   /* if to find sequence rules */
      e = rules(occs, r, 1, supp, &rd);
      if (e < 0) break; }       /* recursively find rules */
    else {                      /* if to find item sets */
      w = (rd.zmax <= 1) ? 0    /* frequent (closed/max.) sequences */
        : recurse(occs, r, 1, &rd);
      if      (target & ISR_MAXIMAL)
        m = (w < rd.smin) && maximal(occs, r, &rd);
      else if (target & ISR_CLOSED) /* check whether to report the */
        m = (w < supp);         /* the 1-item sequence, that is, */
      else m = -1;              /* whether it is closed or maximal */
      if (m && (isr_iset(report, &k, 1, supp, 0, 0) < 0)) {
        e = -1; break; }        /* report the item (1-sequence) */
    }
  }
  if (!(target & ISR_RULES)) {  /* if not to find sequence rules */
    if      (target & ISR_MAXIMAL) m = (max < rd.smin);
    else if (target & ISR_CLOSED)  m = (max < tbg_wgt(tabag));
    else                           m = !(mode & SEQ_ALLOCC);
    if (m) e = isr_report(report);
  }                             /* report the empty sequence */
  free(occs); free(rd.marks);   /* deallocate the arrays */
  return e;                     /* return the error status */
}  /* seqwog() */

/*----------------------------------------------------------------------
  Main Functions
----------------------------------------------------------------------*/
#ifdef SWOG_MAIN

static void help (void)
{                               /* --- print add. option information */
  #ifndef QUIET
  fprintf(stderr, "\n");        /* terminate startup message */
  printf("information output format characters (option -v#)\n");
  printf("  %%%%  a percent sign\n");
  printf("  %%i  number of items (sequence/rule length)\n");
  printf("  %%a  absolute sequence support\n");
  printf("  %%s  relative sequence support as a fraction\n");
  printf("  %%S  relative sequence support as a percentage\n");
  printf("with rules as the target (option -tr) also possible:\n");
  printf("  %%b  absolute body sequence support\n");
  printf("  %%x  relative body sequence support as a fraction\n");
  printf("  %%X  relative body sequence support as a percentage\n");
  printf("  %%h  absolute head item     support\n");
  printf("  %%y  relative head item     support as a fraction\n");
  printf("  %%Y  relative head item     support as a percentage\n");
  printf("  %%c  rule confidence as a fraction\n");
  printf("  %%C  rule confidence as a percentage\n");
  printf("  %%Q  total transaction weight (database size)\n");
  printf("All format characters can be preceded by the number\n");
  printf("of significant digits to be printed (at most 32 digits),\n");
  printf("even though this value is ignored for integer numbers.\n");
  #endif                        /* print help information */
  exit(0);                      /* abort the program */
}  /* help() */

/*--------------------------------------------------------------------*/

static ITEM getbdr (char *s, char **end, double **border)
{                               /* --- get the support border */
  ITEM   i, k;                  /* loop variables */
  double *b;                    /* support border */

  assert(s && end && border);   /* check the function arguments */
  for (i = k = 0; s[i]; i++)    /* traverse the string and */
    if (s[i] == ':') k++;       /* count the number separators */
  *border = b = (double*)malloc((size_t)++k *sizeof(double));
  if (!b) return -1;            /* allocate a support border */
  for (i = 0; i < k; i++) {     /* traverse the parameters */
    b[i] = strtod(s, end);      /* get the next parameter and */
    if (*end == s) break;       /* check for an empty parameter */
    s = *end; if (*s++ != ':') break;
  }                             /* check for a colon */
  if (++i < k)                  /* shrink support array if possible */
    *border = (double*)realloc(b, (size_t)i *sizeof(double));
  return i;                     /* return number of support values */
}  /* getbdr() */

/*--------------------------------------------------------------------*/

static int setbdr (ISREPORT *report, SUPP w, ITEM zmin,
                   double **border, ITEM n)
{                               /* --- set the support border */
  double s;                     /* to traverse the support values */

  assert(report                 /* check the function arguments */
  &&    (w > 0) && (zmin >= 0) && border && (*border || (n <= 0)));
  while (--n >= 0) {            /* traverse the support values */
    s = (*border)[n];           /* transform to absolute count */
    s = ceilsupp((s >= 0) ? s/100.0 *(double)w *(1-DBL_EPSILON) : -s);
    if (isr_setbdr(report, n+zmin, (RSUPP)s) < 0) return -1;
  }                             /* set support in item set reporter */
  if (*border) { free(*border); *border = NULL; }
  return 0;                     /* return 'ok' */
}  /* setbdr() */

/*--------------------------------------------------------------------*/

#ifndef NDEBUG                  /* if debug version */
  #undef  CLEANUP               /* clean up memory and close files */
  #define CLEANUP \
  if (twrite) twr_delete(twrite, 1); \
  if (report) isr_delete(report, 0); \
  if (tabag)  tbg_delete(tabag,  0); \
  if (tread)  trd_delete(tread,  1); \
  if (ibase)  ib_delete (ibase);     \
  if (border) free(border);
#endif

GENERROR(error, exit)           /* generic error reporting function */

/*--------------------------------------------------------------------*/

int main (int argc, char *argv[])
{                               /* --- main function */
  int     i, k = 0;             /* loop variables, counters */
  char    *s;                   /* to traverse the options */
  CCHAR   **optarg = NULL;      /* option argument */
  CCHAR   *fn_inp  = NULL;      /* name of input  file */
  CCHAR   *fn_out  = NULL;      /* name of output file */
  CCHAR   *fn_psp  = NULL;      /* name of pattern spectrum file */
  CCHAR   *recseps = NULL;      /* record  separators */
  CCHAR   *fldseps = NULL;      /* field   separators */
  CCHAR   *blanks  = NULL;      /* blank   characters */
  CCHAR   *comment = NULL;      /* comment characters */
  CCHAR   *hdr     = "";        /* record header  for output */
  CCHAR   *sep     = " ";       /* item separator for output */
  CCHAR   *imp     = " -> ";    /* implication sign for rules */
  CCHAR   *dflt    = " (%S)";   /* default format for check */
  CCHAR   *info    = dflt;      /* format for information output */
  int     target   = 's';       /* target type (closed/maximal) */
  ITEM    zmin     = 1;         /* minimum length of a sequence */
  ITEM    zmax     = ITEM_MAX;  /* maximum length of a sequence */
  double  supp     = 10;        /* minimum support    (in percent) */
  SUPP    smin     = 1;         /* minimum support of an item set */
  SUPP    body     = 10;        /* minimum support of a rule body */
  int     orig     = 0;         /* flag for rule support definition */
  double  conf     = 80;        /* minimum confidence (in percent) */
  int     mode     = 0;         /* operation mode (support counting) */
  int     mtar     = TA_DUPLICS;/* mode for transaction reading */
  int     scan     = 0;         /* flag for scanable item output */
  int     bdrcnt   = 0;         /* number of support values in border */
  int     stats    = 0;         /* flag for sequence statistics */
  int     *marks;               /* marker array for trimming */
  PATSPEC *psp;                 /* collected pattern spectrum */
  SUPP    frq;                  /* frequency of an item */
  ITEM    m, j;                 /* number of items, loop variable */
  TID     n;                    /* number of transactions */
  SUPP    w;                    /* total transaction weight */
  #ifndef QUIET                 /* if not quiet version */
  clock_t t;                    /* timer for measurements */

  prgname = argv[0];            /* get program name for error msgs. */

  /* --- print usage message --- */
  if (argc > 1) {               /* if arguments are given */
    fprintf(stderr, "%s - %s\n", argv[0], DESCRIPTION);
    fprintf(stderr, VERSION); } /* print a startup message */
  else {                        /* if no arguments given */
    printf("usage: %s [options] infile [outfile]\n", argv[0]);
    printf("%s\n", DESCRIPTION);
    printf("%s\n", VERSION);
    printf("-t#      target type                              "
                    "(default: %c)\n", target);
    printf("         (s: frequent, c: closed, m: maximal sequences, "
                     "r: rules)\n");
    printf("         (target type 'r' implies -a (all occurrences))\n");
    printf("-m#      minimum number of items per sequence     "
                    "(default: %"ITEM_FMT")\n", zmin);
    printf("-n#      maximum number of items per sequence     "
                    "(default: no limit)\n");
    printf("-s#      minimum support of a sequence            "
                    "(default: %g%%)\n", supp);
    printf("         (positive: percentage, "
                     "negative: absolute number)\n");
    printf("-o       use original rule support definition     "
                    "(body & head)\n");
    printf("-c#      minimum confidence of a     rule         "
                    "(default: %g%%)\n", conf);
    printf("-a       count all occurrences of a pattern       "
                    "(default: #sequences)\n");
    printf("-F#:#..  support border for filtering item sets   "
                    "(default: none)\n");
    printf("         (list of minimum support values, "
                    "one per item set size,\n");
    printf("         starting at the minimum size, "
                    "as given with option -m#)\n");
    printf("-P#      write pattern spectrum to a file\n");
    printf("-Z       print item set statistics "
                    "(number of item sets per size)\n");
    printf("-g       write output in scanable form "
                    "(quote certain characters)\n");
    printf("-h#      record header  for output                "
                    "(default: \"%s\")\n", hdr);
    printf("-k#      item separator for output                "
                    "(default: \"%s\")\n", sep);
    printf("-I#      implication sign for sequence rules      "
                    "(default: \"%s\")\n", imp);
    printf("-v#      output format for sequence information   "
                    "(default: \"%s\")\n", info);
    printf("-w       integer transaction weight in last field "
                    "(default: only items)\n");
    printf("-r#      record/transaction separators            "
                    "(default: \"\\n\")\n");
    printf("-f#      field /item        separators            "
                    "(default: \" \\t,\")\n");
    printf("-b#      blank   characters                       "
                    "(default: \" \\t\\r\")\n");
    printf("-C#      comment characters                       "
                    "(default: \"#\")\n");
    printf("-!       print additional option information\n");
    printf("infile   file to read transactions from           "
                    "[required]\n");
    printf("outfile  file to write frequent sequences to      "
                    "[optional]\n");
    return 0;                   /* print a usage message */
  }                             /* and abort the program */
  #endif  /* #ifndef QUIET */
  /* free option characters: cdeijlopquwxyz [A-Z]\[CFIPZ] */

  /* --- evaluate arguments --- */
  for (i = 1; i < argc; i++) {  /* traverse arguments */
    s = argv[i];                /* get option argument */
    if (optarg) { *optarg = s; optarg = NULL; continue; }
    if ((*s == '-') && *++s) {  /* -- if argument is an option */
      while (*s) {              /* traverse options */
        switch (*s++) {         /* evaluate switches */
          case '!': help();                          break;
          case 't': target = (*s) ? *s++ : 's';      break;
          case 'm': zmin   = (ITEM)strtol(s, &s, 0); break;
          case 'n': zmax   = (ITEM)strtol(s, &s, 0); break;
          case 's': supp   =       strtod(s, &s);    break;
          case 'o': orig   = 1;                      break;
          case 'c': conf   =       strtod(s, &s);    break;
          case 'a': mode  |= SEQ_ALLOCC;             break;
          case 'F': bdrcnt = getbdr(s, &s, &border); break;
          case 'P': optarg = &fn_psp;                break;
          case 'Z': stats  = 1;                      break;
          case 'g': scan   = 1;                      break;
          case 'h': optarg = &hdr;                   break;
          case 'k': optarg = &sep;                   break;
          case 'I': optarg = &imp;                   break;
          case 'v': optarg = &info;                  break;
          case 'w': mtar  |= TA_WEIGHT;              break;
          case 'r': optarg = &recseps;               break;
          case 'f': optarg = &fldseps;               break;
          case 'b': optarg = &blanks;                break;
          case 'C': optarg = &comment;               break;
          default : error(E_OPTION, *--s);           break;
        }                       /* set option variables */
        if (optarg && *s) { *optarg = s; optarg = NULL; break; }
      } }                       /* get option argument */
    else {                      /* -- if argument is no option */
      switch (k++) {            /* evaluate non-options */
        case  0: fn_inp = s;      break;
        case  1: fn_out = s;      break;
        default: error(E_ARGCNT); break;
      }                         /* note filenames */
    }
  }
  if (optarg)       error(E_OPTARG);     /* check option arguments */
  if (k      < 1)   error(E_ARGCNT);     /* and number of arguments */
  if (zmin   < 0)   error(E_SIZE, zmin); /* check the size limits */
  if (zmax   < 0)   error(E_SIZE, zmax); /* and the minimum support */
  if (supp   > 100) error(E_SUPPORT, supp);
  if (bdrcnt < 0)   error(E_NOMEM);
  if ((conf  < 0) || (conf > 100))
    error(E_CONF, conf);        /* check the minimum confidence */
  switch (target) {             /* check and translate target type */
    case 's': target = ISR_ALL;              break;
    case 'c': target = ISR_CLOSED;           break;
    case 'm': target = ISR_MAXIMAL;          break;
    case 'r': target = ISR_RULES;            break;
    default : error(E_TARGET, (char)target); break;
  }                             /* (get the target type code) */
  if (target < ISR_RULES)       /* remove rule specific settings */
    conf = 100;                 /* (to avoid problems with minsupp) */
  else                          /* for rules, however, all */
    mode |= SEQ_ALLOCC;         /* occurrences must be counted */
  if (info == dflt) {           /* if default info. format is used, */
    if (target != ISR_RULES)    /* set default according to target */
         info = (supp < 0) ? " (%a)"     : " (%S)";
    else info = (supp < 0) ? " (%b, %C)" : " (%X, %C)";
  }                             /* select absolute/relative support */
  MSG(stderr, "\n");            /* terminate the startup message */

  /* --- read transactions --- */
  ibase = ib_create(0, 0);      /* create an item base */
  if (!ibase) error(E_NOMEM);   /* to manage the items */
  tabag = tbg_create(ibase);    /* create a transaction bag */
  if (!tabag) error(E_NOMEM);   /* to store the transactions */
  tread = trd_create();         /* create a table reader and */
  if (!tread) error(E_NOMEM);   /* configure the characters */
  trd_allchs(tread, recseps, fldseps, blanks, "", comment);
  CLOCK(t);                     /* start timer, open input file */
  if (trd_open(tread, NULL, fn_inp) != 0)
    error(E_FOPEN, trd_name(tread));
  MSG(stderr, "reading %s ... ", trd_name(tread));
  k = tbg_read(tabag, tread, mtar);
  if (k < 0) error(-k, tbg_errmsg(tabag, NULL, 0));
  trd_delete(tread, 1);         /* read the transaction database, */
  tread = NULL;                 /* then delete the table reader */
  m = ib_cnt(ibase);            /* get the number of items, */
  n = tbg_cnt(tabag);           /* the number of transactions, */
  w = tbg_wgt(tabag);           /* the total transaction weight */
  MSG(stderr, "[%"ITEM_FMT" item(s), %"TID_FMT, m, n);
  if (w != (SUPP)n) { MSG(stderr, "/%"SUPP_FMT, w); }
  MSG(stderr, " transaction(s)] done [%.2fs].", SEC_SINCE(t));
  if ((m <= 0) || (n <= 0))     /* check for at least one item */
    error(E_NOITEMS);           /* and at least one transaction */
  MSG(stderr, "\n");            /* terminate the log message */
  conf /= 100.0;                /* scale the minimum confidence */
  supp = (supp >= 0) ? supp/100.0 *(double)w *(1-DBL_EPSILON) : -supp;
  body = (SUPP)ceilsupp(supp);  /* compute absolute support values */
  smin = (SUPP)ceilsupp((!(target & ISR_RULES) || orig)
       ? supp : ceilsupp(supp) *conf *(1-DBL_EPSILON));

  /* --- sort and recode items --- */
  CLOCK(t);                     /* start timer, print log message */
  MSG(stderr, "recoding items ... ");
  m = tbg_recode(tabag, 0, -1, -1, (mode & SEQ_ALLOCC) ? -2 : -1);
  if (m <  0) error(E_NOMEM);   /* recode items and transactions */
  if (m <= 0) error(E_NOITEMS); /* and check the number of items */
  MSG(stderr, "[%"ITEM_FMT" item(s)]", m);
  MSG(stderr, " done [%.2fs].\n", SEC_SINCE(t));

  /* --- trim and reduce transactions --- */
  CLOCK(t);                     /* start timer, print log message */
  MSG(stderr, "reducing and triming transactions ... ");
  marks = (int*)malloc((size_t)m *sizeof(int));
  if (!marks) error(E_NOMEM);   /* create a flag array for trimming */
  tbg_sort(tabag, 1, 0);        /* sort the trans. lexicographically */
  for (j = 0; j < m; j++) {     /* traverse the items */
    frq = (mode & SEQ_ALLOCC)   /* get the item frequency */
        ? ib_getxfq(ibase, j) : ib_getfrq(ibase, j);
    marks[j] = (frq >= (SUPP)supp);
  }                             /* set item markers for trimming */
  tbg_trim(tabag, zmin-1, marks, 0); /* trim the transaction database */
  n = tbg_reduce(tabag, 0);     /* reduce transactions to unique ones */
  free(marks);                  /* delete the flag array */
  MSG(stderr, "[%"TID_FMT, n);  /* print number of transactions */
  if (w != (SUPP)n) { MSG(stderr, "/%"SUPP_FMT, w); }
  MSG(stderr, " transaction(s)] done [%.2fs].\n", SEC_SINCE(t));

  /* --- find frequent sequences --- */
  report = isr_createx(ibase, tbg_max(tabag));
  if (!report) error(E_NOMEM);  /* create an item set reporter */
  isr_setsize(report, zmin, zmax);
  isr_setsupp(report, smin, SUPP_MAX);
  if (setbdr(report, w, zmin, &border, bdrcnt) != 0)
    error(E_NOMEM);             /* set the support border */
  if (fn_psp && (isr_addpsp(report, NULL) < 0))
    error(E_NOMEM);             /* set a pattern spectrum if req. */
  if (isr_setfmt(report, scan, hdr, sep, imp, info) != 0)
    error(E_NOMEM);             /* set the output format strings */
  k = isr_open(report, NULL, fn_out);
  if (k) error(k, isr_name(report)); /* open the item set file */
  if ((isr_settarg(report, target, ISR_NOFILTER|ISR_SEQUENCE, -1) < 0)
  ||  (isr_setup(report) < 0))  /* set target and reporting mode */
    error(E_NOMEM);             /* and set up the item set reporter */
  CLOCK(t);                     /* start timer, print log message */
  MSG(stderr, "writing %s ... ", isr_name(report));
  k = seqwog(tabag, target, smin, body, conf, mode, report);
  if (k < 0) error(E_NOMEM);    /* search for frequent sequences */
  MSG(stderr, "[%"ISCNT_FMT" %s(s)]", ISCOUT(isr_repcnt(report)),
              (target == ISR_RULES) ? "rule" : "sequence");
  MSG(stderr, " done [%.2fs].\n", SEC_SINCE(t));
  if (stats)                    /* print item set statistics */
    isr_prstats(report, stdout, 0);
  if (isr_close(report) != 0)   /* close the output file */
    error(E_FWRITE, isr_name(report));

  /* --- write pattern spectrum --- */
  if (fn_psp) {                 /* if to write a pattern spectrum */
    CLOCK(t);                   /* start timer, create table write */
    psp    = isr_getpsp(report);/* get the pattern spectrum */
    twrite = twr_create();      /* create a table writer and */
    if (!twrite) error(E_NOMEM);/* open the output file */
    if (twr_open(twrite, NULL, fn_psp) != 0)
      error(E_FOPEN,  twr_name(twrite));
    MSG(stderr, "writing %s ... ", twr_name(twrite));
    if (psp_report(psp, twrite, 1.0) != 0)
      error(E_FWRITE, twr_name(twrite));
    twr_delete(twrite, 1);      /* write the pattern spectrum */
    twrite = NULL;              /* and delete the table writer */
    MSG(stderr, "[%"SIZE_FMT" signature(s)]", psp_sigcnt(psp));
    MSG(stderr, " done [%.2fs].\n", SEC_SINCE(t));
  }                             /* write a log message */

  /* --- clean up --- */
  CLEANUP;                      /* clean up memory and close files */
  SHOWMEM;                      /* show (final) memory usage */
  return 0;                     /* return 'ok' */
}  /* main() */

#endif
