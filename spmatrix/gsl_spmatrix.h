#ifndef __GSL_SPMATRIX_H__
#define __GSL_SPMATRIX_H__

enum
{
  GSL_SPMATRIX_COO = 0, /* coordinate/triplet representation */
  GSL_SPMATRIX_CSC = 1, /* compressed sparse column */
  GSL_SPMATRIX_CSR = 2, /* compressed sparse row */
  GSL_SPMATRIX_TRIPLET = GSL_SPMATRIX_COO,
  GSL_SPMATRIX_CCS = GSL_SPMATRIX_CSC,
  GSL_SPMATRIX_CRS = GSL_SPMATRIX_CSR
};

/* memory pool for binary tree node allocation */
struct gsl_spmatrix_pool_node
{
  struct gsl_spmatrix_pool_node * next;
  void * block_ptr;          /* pointer to memory block, of size n*tree_node_size */
  unsigned char * free_slot; /* pointer to next available slot */
};

typedef struct gsl_spmatrix_pool_node gsl_spmatrix_pool;

#define GSL_SPMATRIX_ISCOO(m)         ((m)->sptype == GSL_SPMATRIX_COO)
#define GSL_SPMATRIX_ISCSC(m)         ((m)->sptype == GSL_SPMATRIX_CSC)
#define GSL_SPMATRIX_ISCSR(m)         ((m)->sptype == GSL_SPMATRIX_CSR)

#define GSL_SPMATRIX_ISTRIPLET(m)     GSL_SPMATRIX_ISCOO(m)
#define GSL_SPMATRIX_ISCCS(m)         GSL_SPMATRIX_ISCSC(m)
#define GSL_SPMATRIX_ISCRS(m)         GSL_SPMATRIX_ISCSR(m)

#define GSL_SPMATRIX_FLG_GROW         (1 << 0) /* allow size of matrix to grow as elements are added */
#define GSL_SPMATRIX_FLG_FIXED        (1 << 1) /* sparsity pattern is fixed */

/* compare matrix entries (ia,ja) and (ib,jb) - sort by rows first, then by columns */
#define GSL_SPMATRIX_COMPARE_ROWCOL(m,ia,ja,ib,jb)   ((ia) < (ib) ? -1 : ((ia) > (ib) ? 1 : ((ja) < (jb) ? -1 : ((ja) > (jb)))))

/* common/utility functions */

void gsl_spmatrix_cumsum(const size_t n, int * c);

#include <gsl/gsl_spmatrix_complex_long_double.h>
#include <gsl/gsl_spmatrix_complex_double.h>
#include <gsl/gsl_spmatrix_complex_float.h>

#include <gsl/gsl_spmatrix_long_double.h>
#include <gsl/gsl_spmatrix_double.h>
#include <gsl/gsl_spmatrix_float.h>

#include <gsl/gsl_spmatrix_ulong.h>
#include <gsl/gsl_spmatrix_long.h>

#include <gsl/gsl_spmatrix_uint.h>
#include <gsl/gsl_spmatrix_int.h>

#include <gsl/gsl_spmatrix_ushort.h>
#include <gsl/gsl_spmatrix_short.h>

#include <gsl/gsl_spmatrix_uchar.h>
#include <gsl/gsl_spmatrix_char.h>

#endif /* __GSL_SPMATRIX_H__ */
