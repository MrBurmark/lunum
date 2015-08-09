

#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <complex.h>
#include "lunum.h"

#define EXPR_UNM(T) {for(size_t i=0;i<N;++i)((T*)b)[i]=-((T*)a)[i];}
#define EXPR_BNOT(T) {for(size_t i=0;i<N;++i)((T*)b)[i]=~((T*)a)[i];}

#define EXPR_AN_ADD(T) {\
  if (array_first) {\
    for(size_t i=0;i<N;++i)((T*)c)[i]=((T*)a)[i]+((T*)b)[0];\
  } else {\
    for(size_t i=0;i<N;++i)((T*)c)[i]=((T*)b)[0]+((T*)a)[i];\
  }\
}
#define EXPR_AN_SUB(T) {\
  if (array_first) {\
    for(size_t i=0;i<N;++i)((T*)c)[i]=((T*)a)[i]-((T*)b)[0];\
  } else {\
    for(size_t i=0;i<N;++i)((T*)c)[i]=((T*)b)[0]-((T*)a)[i];\
  }\
}
#define EXPR_AN_MUL(T) {\
  if (array_first) {\
    for(size_t i=0;i<N;++i)((T*)c)[i]=((T*)a)[i]*((T*)b)[0];\
  } else {\
    for(size_t i=0;i<N;++i)((T*)c)[i]=((T*)b)[0]*((T*)a)[i];\
  }\
}
#define EXPR_AN_DIV(T) {\
  if (array_first) {\
    for(size_t i=0;i<N;++i)((T*)c)[i]=((T*)a)[i]/((T*)b)[0];\
  } else {\
    for(size_t i=0;i<N;++i)((T*)c)[i]=((T*)b)[0]/((T*)a)[i];\
  }\
}
#define EXPR_AN_IDIVF(T) {\
  if (array_first) {\
    for(size_t i=0;i<N;++i)((T*)c)[i]=(T)floor(((T*)a)[i]/((T*)b)[0]);\
  } else {\
    for(size_t i=0;i<N;++i)((T*)c)[i]=(T)floor(((T*)b)[0]/((T*)a)[i]);\
  }\
}
#define EXPR_AN_MOD(T) {\
  if (array_first) {\
    for(size_t i=0;i<N;++i)((T*)c)[i]=((T*)a)[i]%((T*)b)[0];\
  } else {\
    for(size_t i=0;i<N;++i)((T*)c)[i]=((T*)b)[0]%((T*)a)[i];\
  }\
}
#define EXPR_AN_MODF(T) {\
  if (array_first) {\
    for(size_t i=0;i<N;++i)((T*)c)[i]=((T*)a)[i]-(T)floor(((T*)a)[i]/((T*)b)[0])*((T*)b)[0];\
  } else {\
    for(size_t i=0;i<N;++i)((T*)c)[i]=((T*)b)[0]-(T)floor(((T*)b)[0]/((T*)a)[i])*((T*)a)[i];\
  }\
}
#define EXPR_AN_POW(T) {\
  if (array_first) {\
    for(size_t i=0;i<N;++i)((T*)c)[i]=(T)pow(((T*)a)[i],((T*)b)[0]);\
  } else {\
    for(size_t i=0;i<N;++i)((T*)c)[i]=(T)pow(((T*)b)[0],((T*)a)[i]);\
  }\
}
#define EXPR_AN_COW(T) {\
  if (array_first) {\
    for(size_t i=0;i<N;++i)((T*)c)[i]=(T)cpow(((T*)a)[i],((T*)b)[0]);\
  } else {\
    for(size_t i=0;i<N;++i)((T*)c)[i]=(T)cpow(((T*)b)[0],((T*)a)[i]);\
  }\
}
#define EXPR_AN_BAND(T) {\
  if (array_first) {\
    for(size_t i=0;i<N;++i)((T*)c)[i]=((T*)a)[i]&((T*)b)[0];\
  } else {\
    for(size_t i=0;i<N;++i)((T*)c)[i]=((T*)b)[0]&((T*)a)[i];\
  }\
}
#define EXPR_AN_BOR(T) {\
  if (array_first) {\
    for(size_t i=0;i<N;++i)((T*)c)[i]=((T*)a)[i]|((T*)b)[0];\
  } else {\
    for(size_t i=0;i<N;++i)((T*)c)[i]=((T*)b)[0]|((T*)a)[i];\
  }\
}
#define EXPR_AN_BXOR(T) {\
  if (array_first) {\
    for(size_t i=0;i<N;++i)((T*)c)[i]=((T*)a)[i]^((T*)b)[0];\
  } else {\
    for(size_t i=0;i<N;++i)((T*)c)[i]=((T*)b)[0]^((T*)a)[i];\
  }\
}
#define EXPR_AN_SHL(T) {\
  if (array_first) {\
    for(size_t i=0;i<N;++i)((T*)c)[i]=((T*)a)[i]<<((T*)b)[0];\
  } else {\
    for(size_t i=0;i<N;++i)((T*)c)[i]=((T*)b)[0]<<((T*)a)[i];\
  }\
}
#define EXPR_AN_SHR(T) {\
  if (array_first) {\
    for(size_t i=0;i<N;++i)((T*)c)[i]=((T*)a)[i]>>((T*)b)[0];\
  } else {\
    for(size_t i=0;i<N;++i)((T*)c)[i]=((T*)b)[0]>>((T*)a)[i];\
  }\
}

#define EXPR_AA_ADD(T) {for(size_t i=0;i<N;++i)((T*)c)[i]=((T*)a)[i]+((T*)b)[i];}
#define EXPR_AA_SUB(T) {for(size_t i=0;i<N;++i)((T*)c)[i]=((T*)a)[i]-((T*)b)[i];}
#define EXPR_AA_MUL(T) {for(size_t i=0;i<N;++i)((T*)c)[i]=((T*)a)[i]*((T*)b)[i];}
#define EXPR_AA_DIV(T) {for(size_t i=0;i<N;++i)((T*)c)[i]=((T*)a)[i]/((T*)b)[i];}
#define EXPR_AA_IDIVF(T) {for(size_t i=0;i<N;++i)((T*)c)[i]=(T)floor(((T*)a)[i]/((T*)b)[i]);}
#define EXPR_AA_MOD(T) {for(size_t i=0;i<N;++i)((T*)c)[i]=((T*)a)[i]%((T*)b)[i];}
#define EXPR_AA_MODF(T) {for(size_t i=0;i<N;++i)((T*)c)[i]=((T*)a)[i]-(T)floor(((T*)a)[i]/((T*)b)[i])*((T*)b)[i];}
#define EXPR_AA_POW(T) {for(size_t i=0;i<N;++i)((T*)c)[i]= pow(((T*)a)[i],((T*)b)[i]);}
#define EXPR_AA_COW(T) {for(size_t i=0;i<N;++i)((T*)c)[i]=cpow(((T*)a)[i],((T*)b)[i]);}
#define EXPR_AA_BAND(T) {for(size_t i=0;i<N;++i)((T*)c)[i]=((T*)a)[i]&((T*)b)[i];}
#define EXPR_AA_BOR(T) {for(size_t i=0;i<N;++i)((T*)c)[i]=((T*)a)[i]|((T*)b)[i];}
#define EXPR_AA_BXOR(T) {for(size_t i=0;i<N;++i)((T*)c)[i]=((T*)a)[i]^((T*)b)[i];}
#define EXPR_AA_SHL(T) {for(size_t i=0;i<N;++i)((T*)c)[i]=((T*)a)[i]<<((T*)b)[i];}
#define EXPR_AA_SHR(T) {for(size_t i=0;i<N;++i)((T*)c)[i]=((T*)a)[i]>>((T*)b)[i];}
#define EXPR_ERR(T) {luaL_error(L, "Invalid operation");}

#define ARRAY_OP(op) \
    switch (A->dtype) {\
      case ARRAY_TYPE_BOOL    : op(Bool   ) ; break;\
      case ARRAY_TYPE_CHAR    : op(char   ) ; break;\
      case ARRAY_TYPE_SHORT   : op(short  ) ; break;\
      case ARRAY_TYPE_INT     : op(int    ) ; break;\
      case ARRAY_TYPE_LONG    : op(long   ) ; break;\
      case ARRAY_TYPE_SIZE_T  : op(size_t ) ; break;\
      case ARRAY_TYPE_FLOAT   : op(float  ) ; break;\
      case ARRAY_TYPE_DOUBLE  : op(double ) ; break;\
      case ARRAY_TYPE_COMPLEX : op(Complex) ; break;\
    }

#define ARRAY_OPC(op, opc) \
    switch (A->dtype) {\
      case ARRAY_TYPE_BOOL    : op(Bool   ) ; break;\
      case ARRAY_TYPE_CHAR    : op(char   ) ; break;\
      case ARRAY_TYPE_SHORT   : op(short  ) ; break;\
      case ARRAY_TYPE_INT     : op(int    ) ; break;\
      case ARRAY_TYPE_LONG    : op(long   ) ; break;\
      case ARRAY_TYPE_SIZE_T  : op(size_t ) ; break;\
      case ARRAY_TYPE_FLOAT   : op(float  ) ; break;\
      case ARRAY_TYPE_DOUBLE  : op(double ) ; break;\
      case ARRAY_TYPE_COMPLEX : opc(Complex) ; break;\
    }

#define ARRAY_OPFC(op, opf, opc) \
    switch (A->dtype) {\
      case ARRAY_TYPE_BOOL    : op(Bool   ) ; break;\
      case ARRAY_TYPE_CHAR    : op(char   ) ; break;\
      case ARRAY_TYPE_SHORT   : op(short  ) ; break;\
      case ARRAY_TYPE_INT     : op(int    ) ; break;\
      case ARRAY_TYPE_LONG    : op(long   ) ; break;\
      case ARRAY_TYPE_SIZE_T  : op(size_t ) ; break;\
      case ARRAY_TYPE_FLOAT   : opf(float  ) ; break;\
      case ARRAY_TYPE_DOUBLE  : opf(double ) ; break;\
      case ARRAY_TYPE_COMPLEX : opc(Complex) ; break;\
    }

#define EXPR_ASSIGN0(T,val) {for(size_t i=0;i<N;++i)((T*)a)[i]=val;}
#define EXPR_ASSIGN1(T,val) {for(size_t i=0;i<N;++i)((T*)a)[i]=*((T*)val);}
#define EXPR_ASSIGN2(Ta,Tb) {for(size_t i=0;i<N;++i)((Ta*)a)[i]=((Tb*)b)[i];}

#define ARRAY_ASSIGN_OP1(sw, op, A1) \
    switch (sw) {\
      case ARRAY_TYPE_BOOL    : op(A1, Bool   ) ; break;\
      case ARRAY_TYPE_CHAR    : op(A1, char   ) ; break;\
      case ARRAY_TYPE_SHORT   : op(A1, short  ) ; break;\
      case ARRAY_TYPE_INT     : op(A1, int    ) ; break;\
      case ARRAY_TYPE_LONG    : op(A1, long   ) ; break;\
      case ARRAY_TYPE_SIZE_T  : op(A1, size_t ) ; break;\
      case ARRAY_TYPE_FLOAT   : op(A1, float  ) ; break;\
      case ARRAY_TYPE_DOUBLE  : op(A1, double ) ; break;\
      case ARRAY_TYPE_COMPLEX : op(A1, Complex) ; break;\
    }

#define ARRAY_ASSIGN_OP2(sw, op, A2) \
    switch (sw) {\
      case ARRAY_TYPE_BOOL    : op(Bool,    A2) ; break;\
      case ARRAY_TYPE_CHAR    : op(char,    A2) ; break;\
      case ARRAY_TYPE_SHORT   : op(short,   A2) ; break;\
      case ARRAY_TYPE_INT     : op(int,     A2) ; break;\
      case ARRAY_TYPE_LONG    : op(long,    A2) ; break;\
      case ARRAY_TYPE_SIZE_T  : op(size_t,  A2) ; break;\
      case ARRAY_TYPE_FLOAT   : op(float,   A2) ; break;\
      case ARRAY_TYPE_DOUBLE  : op(double,  A2) ; break;\
      case ARRAY_TYPE_COMPLEX : op(Complex, A2) ; break;\
    }



char *array_typename(ArrayType T)
{
  switch (T) {
  case ARRAY_TYPE_BOOL    : return "bool";
  case ARRAY_TYPE_CHAR    : return "char";
  case ARRAY_TYPE_SHORT   : return "short";
  case ARRAY_TYPE_INT     : return "int";
  case ARRAY_TYPE_LONG    : return "long";
  case ARRAY_TYPE_SIZE_T  : return "size_t";
  case ARRAY_TYPE_FLOAT   : return "float";
  case ARRAY_TYPE_DOUBLE  : return "double";
  case ARRAY_TYPE_COMPLEX : return "complex";
  }
  return NULL; // indicates invalid type
}

ArrayType array_typeflag(char c)
{
  switch (c) {
  case 'b': return ARRAY_TYPE_BOOL;
  case 'c': return ARRAY_TYPE_CHAR;
  case 's': return ARRAY_TYPE_SHORT;
  case 'i': return ARRAY_TYPE_INT;
  case 'l': return ARRAY_TYPE_LONG;
  case 't': return ARRAY_TYPE_SIZE_T;
  case 'f': return ARRAY_TYPE_FLOAT;
  case 'd': return ARRAY_TYPE_DOUBLE;
  case 'z': return ARRAY_TYPE_COMPLEX;
  }
  return -1; // indicates invalid type
}


Array array_new_zeros(size_t N, ArrayType T)
{
  Array A;

  A.data  = malloc(N*array_sizeof(T));
  A.owns  = 1;
  A.size  = N;
  A.dtype = T;
  A.shape = (size_t*) malloc(sizeof(size_t));
  A.ndims = 1;

  A.shape[0] = N;
  void *restrict a = A.data;

  ARRAY_ASSIGN_OP2(T, EXPR_ASSIGN0, 0);

  return A;
}

Array array_new_copy(const Array *B, ArrayType T)
{
  Array A = array_new_zeros(B->size, T);
  array_resize(&A, B->shape, B->ndims);
  array_assign_from_array(&A, B);
  return A;
}

void array_del(Array *A)
{
  if (A->data && A->owns) free(A->data);
  if (A->shape) free(A->shape);

  A->size = 0;
  A->data = NULL;
  A->shape = NULL;
}

int array_resize(Array *A, const size_t *N, int Nd)
{
  size_t ntot = 1;
  for (int d=0; d<Nd; ++d) ntot *= N[d];

  if (A->size != ntot) {
    return 1;
  }
  if (A->shape) free(A->shape);

  A->ndims = Nd;
  A->shape = (size_t*) malloc(Nd*sizeof(size_t));
  memcpy(A->shape, N, Nd*sizeof(size_t));

  return 0;
}


void array_unary_op(lua_State * L, const Array *A, Array *B, 
                    ArrayUnaryOperation op)
{
  const size_t          N = A->size;
  const void  *restrict a = A->data;
  void        *restrict b = B->data;

  switch (op) {
  case ARRAY_OP_UNM:
    ARRAY_OP(EXPR_UNM);
    break;
  case ARRAY_OP_BNOT:
    ARRAY_OPFC(EXPR_BNOT, EXPR_ERR, EXPR_ERR);
    break;
  }
}

void array_array_binary_op(lua_State * L, const Array *A, const Array *B,
                     Array *C, ArrayBinaryOperation op)
{
  const size_t          N = A->size;
  const void  *restrict a = A->data;
  const void  *restrict b = B->data;
  void        *restrict c = C->data;

  switch (op) {
    case ARRAY_OP_ADD:
      ARRAY_OP(EXPR_AA_ADD);
      break;
    case ARRAY_OP_SUB:
      ARRAY_OP(EXPR_AA_SUB);
      break;
    case ARRAY_OP_MUL:
      ARRAY_OP(EXPR_AA_MUL);
      break;
    case ARRAY_OP_DIV:
      ARRAY_OP(EXPR_AA_DIV);
      break;
    case ARRAY_OP_IDIV:
      ARRAY_OPFC(EXPR_AA_DIV, EXPR_AA_IDIVF, EXPR_ERR);
      break;
    case ARRAY_OP_MOD:
      ARRAY_OPFC(EXPR_AA_MOD, EXPR_AA_MODF, EXPR_ERR);
      break;
    case ARRAY_OP_POW:
      ARRAY_OPC(EXPR_AA_POW, EXPR_AA_COW);
      break;
    case ARRAY_OP_BAND:
      ARRAY_OPFC(EXPR_AA_BAND, EXPR_ERR, EXPR_ERR);
      break;
    case ARRAY_OP_BOR:
      ARRAY_OPFC(EXPR_AA_BOR, EXPR_ERR, EXPR_ERR);
      break;
    case ARRAY_OP_BXOR:
      ARRAY_OPFC(EXPR_AA_BXOR, EXPR_ERR, EXPR_ERR);
      break;
    case ARRAY_OP_SHL:
      ARRAY_OPFC(EXPR_AA_SHL, EXPR_ERR, EXPR_ERR);
      break;
    case ARRAY_OP_SHR:
      ARRAY_OPFC(EXPR_AA_SHR, EXPR_ERR, EXPR_ERR);
      break;
  }
}

void array_number_binary_op(lua_State * L, const Array *A, const void *B,
                     Array *C, ArrayBinaryOperation op, Bool array_first)
{
  const size_t          N = A->size;
  const void  *restrict a = A->data;
  const void  *restrict b = B;
  void        *restrict c = C->data;

  switch (op) {
    case ARRAY_OP_ADD:
      ARRAY_OP(EXPR_AN_ADD);
      break;
    case ARRAY_OP_SUB:
      ARRAY_OP(EXPR_AN_SUB);
      break;
    case ARRAY_OP_MUL:
      ARRAY_OP(EXPR_AN_MUL);
      break;
    case ARRAY_OP_DIV:
      ARRAY_OP(EXPR_AN_DIV);
      break;
    case ARRAY_OP_IDIV:
      ARRAY_OPFC(EXPR_AN_DIV, EXPR_AN_IDIVF, EXPR_ERR);
      break;
    case ARRAY_OP_MOD:
      ARRAY_OPFC(EXPR_AN_MOD, EXPR_AN_MODF, EXPR_ERR);
      break;
    case ARRAY_OP_POW:
      ARRAY_OPC(EXPR_AN_POW, EXPR_AN_COW);
      break;
    case ARRAY_OP_BAND:
      ARRAY_OPFC(EXPR_AN_BAND, EXPR_ERR, EXPR_ERR);
      break;
    case ARRAY_OP_BOR:
      ARRAY_OPFC(EXPR_AN_BOR, EXPR_ERR, EXPR_ERR);
      break;
    case ARRAY_OP_BXOR:
      ARRAY_OPFC(EXPR_AN_BXOR, EXPR_ERR, EXPR_ERR);
      break;
    case ARRAY_OP_SHL:
      ARRAY_OPFC(EXPR_AN_SHL, EXPR_ERR, EXPR_ERR);
      break;
    case ARRAY_OP_SHR:
      ARRAY_OPFC(EXPR_AN_SHR, EXPR_ERR, EXPR_ERR);
      break;
  }
}

int array_sizeof(ArrayType T)
{
  switch (T) {
    case ARRAY_TYPE_BOOL    : return sizeof(Bool);
    case ARRAY_TYPE_CHAR    : return sizeof(char);
    case ARRAY_TYPE_SHORT   : return sizeof(short);
    case ARRAY_TYPE_INT     : return sizeof(int);
    case ARRAY_TYPE_LONG    : return sizeof(long);
    case ARRAY_TYPE_SIZE_T  : return sizeof(size_t);
    case ARRAY_TYPE_FLOAT   : return sizeof(float);
    case ARRAY_TYPE_DOUBLE  : return sizeof(double);
    case ARRAY_TYPE_COMPLEX : return sizeof(Complex);
  }
  return 0;  // indicates invalid type
}

void array_assign_from_scalar(Array *A, const void *val)
{
  const size_t N = A->size;
  void        *restrict a = A->data;

  ARRAY_ASSIGN_OP2(A->dtype, EXPR_ASSIGN1, val);
}

void array_assign_from_array(Array *A, const Array *B)
{
  void          *restrict a = A->data;
  const void    *restrict b = B->data;
  const size_t   N = B->size;

  switch (A->dtype) {
    case ARRAY_TYPE_BOOL:
      ARRAY_ASSIGN_OP1(B->dtype, EXPR_ASSIGN2, Bool);
      break;
    case ARRAY_TYPE_CHAR:
      ARRAY_ASSIGN_OP1(B->dtype, EXPR_ASSIGN2, char);
      break;
    case ARRAY_TYPE_SHORT:
      ARRAY_ASSIGN_OP1(B->dtype, EXPR_ASSIGN2, short);
      break;
    case ARRAY_TYPE_INT:
      ARRAY_ASSIGN_OP1(B->dtype, EXPR_ASSIGN2, int);
      break;
    case ARRAY_TYPE_LONG:
      ARRAY_ASSIGN_OP1(B->dtype, EXPR_ASSIGN2, long);
      break;
    case ARRAY_TYPE_SIZE_T:
      ARRAY_ASSIGN_OP1(B->dtype, EXPR_ASSIGN2, size_t);
      break;
    case ARRAY_TYPE_FLOAT:
      ARRAY_ASSIGN_OP1(B->dtype, EXPR_ASSIGN2, float);
      break;
    case ARRAY_TYPE_DOUBLE:
      ARRAY_ASSIGN_OP1(B->dtype, EXPR_ASSIGN2, double);
      break;
    case ARRAY_TYPE_COMPLEX:
      ARRAY_ASSIGN_OP1(B->dtype, EXPR_ASSIGN2, Complex);
      break;
  }
}

#define POS_ROWMAJOR(pos) \
{\
  pos = indices[0];\
  for(int i = 1; i < ndims; ++i) {\
    pos = pos * shape[i] + indices[i];\
  }\
}

#define POS_COLMAJOR(pos) \
{\
  pos = indices[ndims-1];\
  for(int i = ndims-2; i >= 0; --i) {\
    pos = pos * shape[i] + indices[i];\
  }\
}

#define ADVANCE_INDICES() \
{\
  ++(indices[0]);\
  for(int i = 0; i < ndims-1; ++i) {\
    if (indices[i] >= shape[i]){\
      indices[i] = 0;\
      ++(indices[i+1]);\
    } else {\
      break;\
    }\
  }\
}

#define TRANSPOSE(T) \
{\
  while (indices[ndims-1] < shape[ndims-1])\
  {\
    POS_ROWMAJOR(rpos);\
    POS_COLMAJOR(cpos);\
    ((T*)b)[cpos]=((T*)a)[rpos];\
    ADVANCE_INDICES();\
  }\
}

void array_transpose(const Array *A, Array *B)
{
  const int           ndims = A->ndims;
  const void   *restrict a = A->data;
  void         *restrict b = B->data;

  if (ndims < 2) {
    /* copy linear array */
    memcpy(b, a, A->size * array_sizeof(A->dtype));
    return;
  }

  const size_t *restrict shape = A->shape;
  size_t       *restrict indices = (size_t *) calloc(ndims, sizeof(size_t));

  size_t rpos = 0;
  size_t cpos = 0;

  ARRAY_OP(TRANSPOSE);

  free(indices);
}

Array array_new_from_slice(const Array *B1,
				  size_t *restrict start, size_t *restrict stop, size_t *restrict skip, int Nd)
// -----------------------------------------------------------------------------
// Extracts a slice from B1, and returns it as the contiguous array 'B0'
// -----------------------------------------------------------------------------
// @start : starting indices into B1
// @stop  : upper bound on selection (non-inclusive)
// @skip  : distance between entries of B1 along each axis
// @Nd    : the number of axes in each array
// -----------------------------------------------------------------------------
{

  size_t *restrict J = (size_t*) malloc(Nd*sizeof(size_t)); // current indices into B1
  size_t *restrict N = (size_t*) malloc(Nd*sizeof(size_t)); // number of elements to select
  size_t *restrict S = (size_t*) malloc(Nd*sizeof(size_t)); // strides (in memory) along each axis

  size_t ntot = 1;

  for (int d=0; d<Nd; ++d) {
    J[d] = 0;
    N[d] = 1 + (stop[d] - start[d] - 1) / skip[d];
    ntot *= N[d];
  }

  S[Nd-1] = 1;
  for (int d=Nd-2; d>=0; --d) S[d] = S[d+1] * B1->shape[d+1];


  Array B0 = array_new_zeros(ntot, B1->dtype);
  array_resize(&B0, N, Nd);
  int sizeof_T = array_sizeof(B0.dtype);
  size_t m = 0; // indexes into B0, advanced uniformly


  char *restrict b0 = (char*) B0 .data;
  char *restrict b1 = (char*) B1->data;


  while (J[0] < N[0]) {

    size_t M = 0;
    for (int d=0; d<Nd; ++d) M += (J[d] * skip[d] + start[d]) * S[d];

    // ----- use the indices m,M -----
    memcpy(b0 + (m++)*sizeof_T, b1 + M*sizeof_T, sizeof_T);
    // -----                 -----

    ++J[Nd-1];
    for (int d=Nd-1; d!=0; --d) {
      if (J[d] == N[d]) {
	J[d] = 0;
	++J[d-1];
      }
    }
  }

  free(J);
  free(N);
  free(S);

  return B0;
}

Array array_new_from_mask(const Array *B1, Array *M)
// -----------------------------------------------------------------------------
// Extracts the indices of B1 for which M is true, and returns a 1d-array
// -----------------------------------------------------------------------------
// @M : Array of bool's, must have the same size as B1
// -----------------------------------------------------------------------------
{
  int sizeof_T = array_sizeof(B1->dtype);

  char *restrict b0 = (char*) malloc(sizeof_T);
  char *restrict b1 = (char*) B1->data;
  Bool *restrict Mp = (Bool*) M->data;

  size_t m = 0;

  for (size_t n=0; n<B1->size; ++n) {
    if (Mp[n]) {
      b0 = (char*) realloc(b0, (++m)*sizeof(double));
      memcpy(b0 + (m-1)*sizeof_T, b1 + n*sizeof_T, sizeof_T);
    }
  }

  Array B0 = array_new_zeros(m, B1->dtype);
  memcpy(B0.data, b0, m*sizeof_T);
  free(b0);

  return B0;
}
