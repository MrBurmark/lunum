

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
#define EXPR_ERR() {luaL_error(L, "Invalid operation");}

#define EXPR_ASSIGN0(T,val) {for(size_t i=0;i<N;++i)((T*)a)[i]=val;}
#define EXPR_ASSIGN1(T,val) {for(size_t i=0;i<N;++i)((T*)a)[i]=*((T*)val);}
#define EXPR_ASSIGN2(Ta,Tb) {for(size_t i=0;i<N;++i)((Ta*)a)[i]=((Tb*)b)[i];}



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

  switch (T) {
  case ARRAY_TYPE_BOOL    : EXPR_ASSIGN0(Bool   ,0) ; break;
  case ARRAY_TYPE_CHAR    : EXPR_ASSIGN0(char   ,0) ; break;
  case ARRAY_TYPE_SHORT   : EXPR_ASSIGN0(short  ,0) ; break;
  case ARRAY_TYPE_INT     : EXPR_ASSIGN0(int    ,0) ; break;
  case ARRAY_TYPE_LONG    : EXPR_ASSIGN0(long   ,0) ; break;
  case ARRAY_TYPE_SIZE_T  : EXPR_ASSIGN0(size_t ,0) ; break;
  case ARRAY_TYPE_FLOAT   : EXPR_ASSIGN0(float  ,0) ; break;
  case ARRAY_TYPE_DOUBLE  : EXPR_ASSIGN0(double ,0) ; break;
  case ARRAY_TYPE_COMPLEX : EXPR_ASSIGN0(Complex,0) ; break;
  }

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
    switch (A->dtype) {
    case ARRAY_TYPE_BOOL    : EXPR_UNM(Bool   ) ; break;
    case ARRAY_TYPE_CHAR    : EXPR_UNM(char   ) ; break;
    case ARRAY_TYPE_SHORT   : EXPR_UNM(short  ) ; break;
    case ARRAY_TYPE_INT     : EXPR_UNM(int    ) ; break;
    case ARRAY_TYPE_LONG    : EXPR_UNM(long   ) ; break;
    case ARRAY_TYPE_SIZE_T  : EXPR_UNM(size_t ) ; break;
    case ARRAY_TYPE_FLOAT   : EXPR_UNM(float  ) ; break;
    case ARRAY_TYPE_DOUBLE  : EXPR_UNM(double ) ; break;
    case ARRAY_TYPE_COMPLEX : EXPR_UNM(Complex) ; break;
    }
    break;
  case ARRAY_OP_BNOT:
    switch (A->dtype) {
    case ARRAY_TYPE_BOOL    : EXPR_BNOT(Bool   ) ; break;
    case ARRAY_TYPE_CHAR    : EXPR_BNOT(char   ) ; break;
    case ARRAY_TYPE_SHORT   : EXPR_BNOT(short  ) ; break;
    case ARRAY_TYPE_INT     : EXPR_BNOT(int    ) ; break;
    case ARRAY_TYPE_LONG    : EXPR_BNOT(long   ) ; break;
    case ARRAY_TYPE_SIZE_T  : EXPR_BNOT(size_t ) ; break;
    case ARRAY_TYPE_FLOAT   : EXPR_ERR(        ) ; break;
    case ARRAY_TYPE_DOUBLE  : EXPR_ERR(        ) ; break;
    case ARRAY_TYPE_COMPLEX : EXPR_ERR(        ) ; break;
    }
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
    switch (A->dtype) {
    case ARRAY_TYPE_BOOL    : EXPR_AA_ADD(Bool   ) ; break;
    case ARRAY_TYPE_CHAR    : EXPR_AA_ADD(char   ) ; break;
    case ARRAY_TYPE_SHORT   : EXPR_AA_ADD(short  ) ; break;
    case ARRAY_TYPE_INT     : EXPR_AA_ADD(int    ) ; break;
    case ARRAY_TYPE_LONG    : EXPR_AA_ADD(long   ) ; break;
    case ARRAY_TYPE_SIZE_T  : EXPR_AA_ADD(size_t ) ; break;
    case ARRAY_TYPE_FLOAT   : EXPR_AA_ADD(float  ) ; break;
    case ARRAY_TYPE_DOUBLE  : EXPR_AA_ADD(double ) ; break;
    case ARRAY_TYPE_COMPLEX : EXPR_AA_ADD(Complex) ; break;
    }
    break;
  case ARRAY_OP_SUB:
    switch (A->dtype) {
    case ARRAY_TYPE_BOOL    : EXPR_AA_SUB(Bool   ) ; break;
    case ARRAY_TYPE_CHAR    : EXPR_AA_SUB(char   ) ; break;
    case ARRAY_TYPE_SHORT   : EXPR_AA_SUB(short  ) ; break;
    case ARRAY_TYPE_INT     : EXPR_AA_SUB(int    ) ; break;
    case ARRAY_TYPE_LONG    : EXPR_AA_SUB(long   ) ; break;
    case ARRAY_TYPE_SIZE_T  : EXPR_AA_SUB(size_t ) ; break;
    case ARRAY_TYPE_FLOAT   : EXPR_AA_SUB(float  ) ; break;
    case ARRAY_TYPE_DOUBLE  : EXPR_AA_SUB(double ) ; break;
    case ARRAY_TYPE_COMPLEX : EXPR_AA_SUB(Complex) ; break;
    }
    break;
  case ARRAY_OP_MUL:
    switch (A->dtype) {
    case ARRAY_TYPE_BOOL    : EXPR_AA_MUL(Bool   ) ; break;
    case ARRAY_TYPE_CHAR    : EXPR_AA_MUL(char   ) ; break;
    case ARRAY_TYPE_SHORT   : EXPR_AA_MUL(short  ) ; break;
    case ARRAY_TYPE_INT     : EXPR_AA_MUL(int    ) ; break;
    case ARRAY_TYPE_LONG    : EXPR_AA_MUL(long   ) ; break;
    case ARRAY_TYPE_SIZE_T  : EXPR_AA_MUL(size_t ) ; break;
    case ARRAY_TYPE_FLOAT   : EXPR_AA_MUL(float  ) ; break;
    case ARRAY_TYPE_DOUBLE  : EXPR_AA_MUL(double ) ; break;
    case ARRAY_TYPE_COMPLEX : EXPR_AA_MUL(Complex) ; break;
    }
    break;
  case ARRAY_OP_DIV:
    switch (A->dtype) {
    case ARRAY_TYPE_BOOL    : EXPR_AA_DIV(Bool   ) ; break;
    case ARRAY_TYPE_CHAR    : EXPR_AA_DIV(char   ) ; break;
    case ARRAY_TYPE_SHORT   : EXPR_AA_DIV(short  ) ; break;
    case ARRAY_TYPE_INT     : EXPR_AA_DIV(int    ) ; break;
    case ARRAY_TYPE_LONG    : EXPR_AA_DIV(long   ) ; break;
    case ARRAY_TYPE_SIZE_T  : EXPR_AA_DIV(size_t ) ; break;
    case ARRAY_TYPE_FLOAT   : EXPR_AA_DIV(float  ) ; break;
    case ARRAY_TYPE_DOUBLE  : EXPR_AA_DIV(double ) ; break;
    case ARRAY_TYPE_COMPLEX : EXPR_AA_DIV(Complex) ; break;
    }
    break;
  case ARRAY_OP_IDIV:
    switch (A->dtype) {
    case ARRAY_TYPE_BOOL    : EXPR_AA_DIV(Bool    ) ; break;
    case ARRAY_TYPE_CHAR    : EXPR_AA_DIV(char    ) ; break;
    case ARRAY_TYPE_SHORT   : EXPR_AA_DIV(short   ) ; break;
    case ARRAY_TYPE_INT     : EXPR_AA_DIV(int     ) ; break;
    case ARRAY_TYPE_LONG    : EXPR_AA_DIV(long    ) ; break;
    case ARRAY_TYPE_SIZE_T  : EXPR_AA_DIV(size_t  ) ; break;
    case ARRAY_TYPE_FLOAT   : EXPR_AA_IDIVF(float ) ; break;
    case ARRAY_TYPE_DOUBLE  : EXPR_AA_IDIVF(double) ; break;
    case ARRAY_TYPE_COMPLEX : EXPR_ERR(        ) ; break;
    }
    break;
  case ARRAY_OP_MOD:
    switch (A->dtype) {
    case ARRAY_TYPE_BOOL    : EXPR_AA_MOD(Bool   ) ; break;
    case ARRAY_TYPE_CHAR    : EXPR_AA_MOD(char   ) ; break;
    case ARRAY_TYPE_SHORT   : EXPR_AA_MOD(short  ) ; break;
    case ARRAY_TYPE_INT     : EXPR_AA_MOD(int    ) ; break;
    case ARRAY_TYPE_LONG    : EXPR_AA_MOD(long   ) ; break;
    case ARRAY_TYPE_SIZE_T  : EXPR_AA_MOD(size_t ) ; break;
    case ARRAY_TYPE_FLOAT   : EXPR_AA_MODF(float ) ; break;
    case ARRAY_TYPE_DOUBLE  : EXPR_AA_MODF(double) ; break;
    case ARRAY_TYPE_COMPLEX : EXPR_ERR(       ) ; break;
    }
    break;
  case ARRAY_OP_POW:
    switch (A->dtype) {
    case ARRAY_TYPE_BOOL    : EXPR_AA_POW(Bool   ) ; break;
    case ARRAY_TYPE_CHAR    : EXPR_AA_POW(char   ) ; break;
    case ARRAY_TYPE_SHORT   : EXPR_AA_POW(short  ) ; break;
    case ARRAY_TYPE_INT     : EXPR_AA_POW(int    ) ; break;
    case ARRAY_TYPE_LONG    : EXPR_AA_POW(long   ) ; break;
    case ARRAY_TYPE_SIZE_T  : EXPR_AA_POW(size_t ) ; break;
    case ARRAY_TYPE_FLOAT   : EXPR_AA_POW(float  ) ; break;
    case ARRAY_TYPE_DOUBLE  : EXPR_AA_POW(double ) ; break;
    case ARRAY_TYPE_COMPLEX : EXPR_AA_COW(Complex) ; break;
    }
    break;
  case ARRAY_OP_BAND:
    switch (A->dtype) {
    case ARRAY_TYPE_BOOL    : EXPR_AA_BAND(Bool   ) ; break;
    case ARRAY_TYPE_CHAR    : EXPR_AA_BAND(char   ) ; break;
    case ARRAY_TYPE_SHORT   : EXPR_AA_BAND(short  ) ; break;
    case ARRAY_TYPE_INT     : EXPR_AA_BAND(int    ) ; break;
    case ARRAY_TYPE_LONG    : EXPR_AA_BAND(long   ) ; break;
    case ARRAY_TYPE_SIZE_T  : EXPR_AA_BAND(size_t ) ; break;
    case ARRAY_TYPE_FLOAT   : EXPR_ERR(        ) ; break;
    case ARRAY_TYPE_DOUBLE  : EXPR_ERR(        ) ; break;
    case ARRAY_TYPE_COMPLEX : EXPR_ERR(        ) ; break;
    }
    break;
  case ARRAY_OP_BOR:
    switch (A->dtype) {
    case ARRAY_TYPE_BOOL    : EXPR_AA_BOR(Bool    ) ; break;
    case ARRAY_TYPE_CHAR    : EXPR_AA_BOR(char    ) ; break;
    case ARRAY_TYPE_SHORT   : EXPR_AA_BOR(short   ) ; break;
    case ARRAY_TYPE_INT     : EXPR_AA_BOR(int     ) ; break;
    case ARRAY_TYPE_LONG    : EXPR_AA_BOR(long    ) ; break;
    case ARRAY_TYPE_SIZE_T  : EXPR_AA_BOR(size_t  ) ; break;
    case ARRAY_TYPE_FLOAT   : EXPR_ERR(        ) ; break;
    case ARRAY_TYPE_DOUBLE  : EXPR_ERR(        ) ; break;
    case ARRAY_TYPE_COMPLEX : EXPR_ERR(        ) ; break;
    }
    break;
  case ARRAY_OP_BXOR:
    switch (A->dtype) {
    case ARRAY_TYPE_BOOL    : EXPR_AA_BXOR(Bool   ) ; break;
    case ARRAY_TYPE_CHAR    : EXPR_AA_BXOR(char   ) ; break;
    case ARRAY_TYPE_SHORT   : EXPR_AA_BXOR(short  ) ; break;
    case ARRAY_TYPE_INT     : EXPR_AA_BXOR(int    ) ; break;
    case ARRAY_TYPE_LONG    : EXPR_AA_BXOR(long   ) ; break;
    case ARRAY_TYPE_SIZE_T  : EXPR_AA_BXOR(size_t ) ; break;
    case ARRAY_TYPE_FLOAT   : EXPR_ERR(        ) ; break;
    case ARRAY_TYPE_DOUBLE  : EXPR_ERR(        ) ; break;
    case ARRAY_TYPE_COMPLEX : EXPR_ERR(        ) ; break;
    }
    break;
  case ARRAY_OP_SHL:
    switch (A->dtype) {
    case ARRAY_TYPE_BOOL    : EXPR_AA_SHL(Bool    ) ; break;
    case ARRAY_TYPE_CHAR    : EXPR_AA_SHL(char    ) ; break;
    case ARRAY_TYPE_SHORT   : EXPR_AA_SHL(short   ) ; break;
    case ARRAY_TYPE_INT     : EXPR_AA_SHL(int     ) ; break;
    case ARRAY_TYPE_LONG    : EXPR_AA_SHL(long    ) ; break;
    case ARRAY_TYPE_SIZE_T  : EXPR_AA_SHL(size_t  ) ; break;
    case ARRAY_TYPE_FLOAT   : EXPR_ERR(        ) ; break;
    case ARRAY_TYPE_DOUBLE  : EXPR_ERR(        ) ; break;
    case ARRAY_TYPE_COMPLEX : EXPR_ERR(        ) ; break;
    }
    break;
  case ARRAY_OP_SHR:
    switch (A->dtype) {
    case ARRAY_TYPE_BOOL    : EXPR_AA_SHR(Bool    ) ; break;
    case ARRAY_TYPE_CHAR    : EXPR_AA_SHR(char    ) ; break;
    case ARRAY_TYPE_SHORT   : EXPR_AA_SHR(short   ) ; break;
    case ARRAY_TYPE_INT     : EXPR_AA_SHR(int     ) ; break;
    case ARRAY_TYPE_LONG    : EXPR_AA_SHR(long    ) ; break;
    case ARRAY_TYPE_SIZE_T  : EXPR_AA_SHR(size_t  ) ; break;
    case ARRAY_TYPE_FLOAT   : EXPR_ERR(        ) ; break;
    case ARRAY_TYPE_DOUBLE  : EXPR_ERR(        ) ; break;
    case ARRAY_TYPE_COMPLEX : EXPR_ERR(        ) ; break;
    }
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
    switch (A->dtype) {
    case ARRAY_TYPE_BOOL    : EXPR_AN_ADD(Bool   ) ; break;
    case ARRAY_TYPE_CHAR    : EXPR_AN_ADD(char   ) ; break;
    case ARRAY_TYPE_SHORT   : EXPR_AN_ADD(short  ) ; break;
    case ARRAY_TYPE_INT     : EXPR_AN_ADD(int    ) ; break;
    case ARRAY_TYPE_LONG    : EXPR_AN_ADD(long   ) ; break;
    case ARRAY_TYPE_SIZE_T  : EXPR_AN_ADD(size_t ) ; break;
    case ARRAY_TYPE_FLOAT   : EXPR_AN_ADD(float  ) ; break;
    case ARRAY_TYPE_DOUBLE  : EXPR_AN_ADD(double ) ; break;
    case ARRAY_TYPE_COMPLEX : EXPR_AN_ADD(Complex) ; break;
    }
    break;
  case ARRAY_OP_SUB:
    switch (A->dtype) {
    case ARRAY_TYPE_BOOL    : EXPR_AN_SUB(Bool   ) ; break;
    case ARRAY_TYPE_CHAR    : EXPR_AN_SUB(char   ) ; break;
    case ARRAY_TYPE_SHORT   : EXPR_AN_SUB(short  ) ; break;
    case ARRAY_TYPE_INT     : EXPR_AN_SUB(int    ) ; break;
    case ARRAY_TYPE_LONG    : EXPR_AN_SUB(long   ) ; break;
    case ARRAY_TYPE_SIZE_T  : EXPR_AN_SUB(size_t ) ; break;
    case ARRAY_TYPE_FLOAT   : EXPR_AN_SUB(float  ) ; break;
    case ARRAY_TYPE_DOUBLE  : EXPR_AN_SUB(double ) ; break;
    case ARRAY_TYPE_COMPLEX : EXPR_AN_SUB(Complex) ; break;
    }
    break;
  case ARRAY_OP_MUL:
    switch (A->dtype) {
    case ARRAY_TYPE_BOOL    : EXPR_AN_MUL(Bool   ) ; break;
    case ARRAY_TYPE_CHAR    : EXPR_AN_MUL(char   ) ; break;
    case ARRAY_TYPE_SHORT   : EXPR_AN_MUL(short  ) ; break;
    case ARRAY_TYPE_INT     : EXPR_AN_MUL(int    ) ; break;
    case ARRAY_TYPE_LONG    : EXPR_AN_MUL(long   ) ; break;
    case ARRAY_TYPE_SIZE_T  : EXPR_AN_MUL(size_t ) ; break;
    case ARRAY_TYPE_FLOAT   : EXPR_AN_MUL(float  ) ; break;
    case ARRAY_TYPE_DOUBLE  : EXPR_AN_MUL(double ) ; break;
    case ARRAY_TYPE_COMPLEX : EXPR_AN_MUL(Complex) ; break;
    }
    break;
  case ARRAY_OP_DIV:
    switch (A->dtype) {
    case ARRAY_TYPE_BOOL    : EXPR_AN_DIV(Bool   ) ; break;
    case ARRAY_TYPE_CHAR    : EXPR_AN_DIV(char   ) ; break;
    case ARRAY_TYPE_SHORT   : EXPR_AN_DIV(short  ) ; break;
    case ARRAY_TYPE_INT     : EXPR_AN_DIV(int    ) ; break;
    case ARRAY_TYPE_LONG    : EXPR_AN_DIV(long   ) ; break;
    case ARRAY_TYPE_SIZE_T  : EXPR_AN_DIV(size_t ) ; break;
    case ARRAY_TYPE_FLOAT   : EXPR_AN_DIV(float  ) ; break;
    case ARRAY_TYPE_DOUBLE  : EXPR_AN_DIV(double ) ; break;
    case ARRAY_TYPE_COMPLEX : EXPR_AN_DIV(Complex) ; break;
    }
    break;
  case ARRAY_OP_IDIV:
    switch (A->dtype) {
    case ARRAY_TYPE_BOOL    : EXPR_AN_DIV(Bool    ) ; break;
    case ARRAY_TYPE_CHAR    : EXPR_AN_DIV(char    ) ; break;
    case ARRAY_TYPE_SHORT   : EXPR_AN_DIV(short   ) ; break;
    case ARRAY_TYPE_INT     : EXPR_AN_DIV(int     ) ; break;
    case ARRAY_TYPE_LONG    : EXPR_AN_DIV(long    ) ; break;
    case ARRAY_TYPE_SIZE_T  : EXPR_AN_DIV(size_t  ) ; break;
    case ARRAY_TYPE_FLOAT   : EXPR_AN_IDIVF(float ) ; break;
    case ARRAY_TYPE_DOUBLE  : EXPR_AN_IDIVF(double) ; break;
    case ARRAY_TYPE_COMPLEX : EXPR_ERR(        ) ; break;
    }
    break;
  case ARRAY_OP_MOD:
    switch (A->dtype) {
    case ARRAY_TYPE_BOOL    : EXPR_AN_MOD(Bool   ) ; break;
    case ARRAY_TYPE_CHAR    : EXPR_AN_MOD(char   ) ; break;
    case ARRAY_TYPE_SHORT   : EXPR_AN_MOD(short  ) ; break;
    case ARRAY_TYPE_INT     : EXPR_AN_MOD(int    ) ; break;
    case ARRAY_TYPE_LONG    : EXPR_AN_MOD(long   ) ; break;
    case ARRAY_TYPE_SIZE_T  : EXPR_AN_MOD(size_t ) ; break;
    case ARRAY_TYPE_FLOAT   : EXPR_AN_MODF(float ) ; break;
    case ARRAY_TYPE_DOUBLE  : EXPR_AN_MODF(double) ; break;
    case ARRAY_TYPE_COMPLEX : EXPR_ERR(       ) ; break;
    }
    break;
  case ARRAY_OP_POW:
    switch (A->dtype) {
    case ARRAY_TYPE_BOOL    : EXPR_AN_POW(Bool   ) ; break;
    case ARRAY_TYPE_CHAR    : EXPR_AN_POW(char   ) ; break;
    case ARRAY_TYPE_SHORT   : EXPR_AN_POW(short  ) ; break;
    case ARRAY_TYPE_INT     : EXPR_AN_POW(int    ) ; break;
    case ARRAY_TYPE_LONG    : EXPR_AN_POW(long   ) ; break;
    case ARRAY_TYPE_SIZE_T  : EXPR_AN_POW(size_t ) ; break;
    case ARRAY_TYPE_FLOAT   : EXPR_AN_POW(float  ) ; break;
    case ARRAY_TYPE_DOUBLE  : EXPR_AN_POW(double ) ; break;
    case ARRAY_TYPE_COMPLEX : EXPR_AN_COW(Complex) ; break;
    }
    break;
  case ARRAY_OP_BAND:
    switch (A->dtype) {
    case ARRAY_TYPE_BOOL    : EXPR_AN_BAND(Bool   ) ; break;
    case ARRAY_TYPE_CHAR    : EXPR_AN_BAND(char   ) ; break;
    case ARRAY_TYPE_SHORT   : EXPR_AN_BAND(short  ) ; break;
    case ARRAY_TYPE_INT     : EXPR_AN_BAND(int    ) ; break;
    case ARRAY_TYPE_LONG    : EXPR_AN_BAND(long   ) ; break;
    case ARRAY_TYPE_SIZE_T  : EXPR_AN_BAND(size_t ) ; break;
    case ARRAY_TYPE_FLOAT   : EXPR_ERR(        ) ; break;
    case ARRAY_TYPE_DOUBLE  : EXPR_ERR(        ) ; break;
    case ARRAY_TYPE_COMPLEX : EXPR_ERR(        ) ; break;
    }
    break;
  case ARRAY_OP_BOR:
    switch (A->dtype) {
    case ARRAY_TYPE_BOOL    : EXPR_AN_BOR(Bool    ) ; break;
    case ARRAY_TYPE_CHAR    : EXPR_AN_BOR(char    ) ; break;
    case ARRAY_TYPE_SHORT   : EXPR_AN_BOR(short   ) ; break;
    case ARRAY_TYPE_INT     : EXPR_AN_BOR(int     ) ; break;
    case ARRAY_TYPE_LONG    : EXPR_AN_BOR(long    ) ; break;
    case ARRAY_TYPE_SIZE_T  : EXPR_AN_BOR(size_t  ) ; break;
    case ARRAY_TYPE_FLOAT   : EXPR_ERR(        ) ; break;
    case ARRAY_TYPE_DOUBLE  : EXPR_ERR(        ) ; break;
    case ARRAY_TYPE_COMPLEX : EXPR_ERR(        ) ; break;
    }
    break;
  case ARRAY_OP_BXOR:
    switch (A->dtype) {
    case ARRAY_TYPE_BOOL    : EXPR_AN_BXOR(Bool   ) ; break;
    case ARRAY_TYPE_CHAR    : EXPR_AN_BXOR(char   ) ; break;
    case ARRAY_TYPE_SHORT   : EXPR_AN_BXOR(short  ) ; break;
    case ARRAY_TYPE_INT     : EXPR_AN_BXOR(int    ) ; break;
    case ARRAY_TYPE_LONG    : EXPR_AN_BXOR(long   ) ; break;
    case ARRAY_TYPE_SIZE_T  : EXPR_AN_BXOR(size_t ) ; break;
    case ARRAY_TYPE_FLOAT   : EXPR_ERR(        ) ; break;
    case ARRAY_TYPE_DOUBLE  : EXPR_ERR(        ) ; break;
    case ARRAY_TYPE_COMPLEX : EXPR_ERR(        ) ; break;
    }
    break;
  case ARRAY_OP_SHL:
    switch (A->dtype) {
    case ARRAY_TYPE_BOOL    : EXPR_AN_SHL(Bool    ) ; break;
    case ARRAY_TYPE_CHAR    : EXPR_AN_SHL(char    ) ; break;
    case ARRAY_TYPE_SHORT   : EXPR_AN_SHL(short   ) ; break;
    case ARRAY_TYPE_INT     : EXPR_AN_SHL(int     ) ; break;
    case ARRAY_TYPE_LONG    : EXPR_AN_SHL(long    ) ; break;
    case ARRAY_TYPE_SIZE_T  : EXPR_AN_SHL(size_t  ) ; break;
    case ARRAY_TYPE_FLOAT   : EXPR_ERR(        ) ; break;
    case ARRAY_TYPE_DOUBLE  : EXPR_ERR(        ) ; break;
    case ARRAY_TYPE_COMPLEX : EXPR_ERR(        ) ; break;
    }
    break;
  case ARRAY_OP_SHR:
    switch (A->dtype) {
    case ARRAY_TYPE_BOOL    : EXPR_AN_SHR(Bool    ) ; break;
    case ARRAY_TYPE_CHAR    : EXPR_AN_SHR(char    ) ; break;
    case ARRAY_TYPE_SHORT   : EXPR_AN_SHR(short   ) ; break;
    case ARRAY_TYPE_INT     : EXPR_AN_SHR(int     ) ; break;
    case ARRAY_TYPE_LONG    : EXPR_AN_SHR(long    ) ; break;
    case ARRAY_TYPE_SIZE_T  : EXPR_AN_SHR(size_t  ) ; break;
    case ARRAY_TYPE_FLOAT   : EXPR_ERR(        ) ; break;
    case ARRAY_TYPE_DOUBLE  : EXPR_ERR(        ) ; break;
    case ARRAY_TYPE_COMPLEX : EXPR_ERR(        ) ; break;
    }
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
  return sizeof(int);
}

void array_assign_from_scalar(Array *A, const void *val)
{
  const size_t N = A->size;
  void        *restrict a = A->data;

  switch (A->dtype) {
  case ARRAY_TYPE_BOOL    : EXPR_ASSIGN1(Bool   , val) ; break;
  case ARRAY_TYPE_CHAR    : EXPR_ASSIGN1(char   , val) ; break;
  case ARRAY_TYPE_SHORT   : EXPR_ASSIGN1(short  , val) ; break;
  case ARRAY_TYPE_INT     : EXPR_ASSIGN1(int    , val) ; break;
  case ARRAY_TYPE_LONG    : EXPR_ASSIGN1(long   , val) ; break;
  case ARRAY_TYPE_SIZE_T  : EXPR_ASSIGN1(size_t , val) ; break;
  case ARRAY_TYPE_FLOAT   : EXPR_ASSIGN1(float  , val) ; break;
  case ARRAY_TYPE_DOUBLE  : EXPR_ASSIGN1(double , val) ; break;
  case ARRAY_TYPE_COMPLEX : EXPR_ASSIGN1(Complex, val) ; break;
  }
}

void array_assign_from_array(Array *A, const Array *B)
{
  void          *restrict a = A->data;
  const void    *restrict b = B->data;
  const size_t   N = B->size;

  switch (A->dtype) {
  case ARRAY_TYPE_BOOL:
    switch (B->dtype) {
      //                               (A->type, B->type)
    case ARRAY_TYPE_BOOL    : EXPR_ASSIGN2(Bool, Bool)    ; break;
    case ARRAY_TYPE_CHAR    : EXPR_ASSIGN2(Bool, char)    ; break;
    case ARRAY_TYPE_SHORT   : EXPR_ASSIGN2(Bool, short)   ; break;
    case ARRAY_TYPE_INT     : EXPR_ASSIGN2(Bool, int)     ; break;
    case ARRAY_TYPE_LONG    : EXPR_ASSIGN2(Bool, long)    ; break;
    case ARRAY_TYPE_SIZE_T  : EXPR_ASSIGN2(Bool, size_t)  ; break;
    case ARRAY_TYPE_FLOAT   : EXPR_ASSIGN2(Bool, float)   ; break;
    case ARRAY_TYPE_DOUBLE  : EXPR_ASSIGN2(Bool, double)  ; break;
    case ARRAY_TYPE_COMPLEX : EXPR_ASSIGN2(Bool, Complex) ; break;
    }
    break;

  case ARRAY_TYPE_CHAR:
    switch (B->dtype) {
    case ARRAY_TYPE_BOOL    : EXPR_ASSIGN2(char, Bool)    ; break;
    case ARRAY_TYPE_CHAR    : EXPR_ASSIGN2(char, char)    ; break;
    case ARRAY_TYPE_SHORT   : EXPR_ASSIGN2(char, short)   ; break;
    case ARRAY_TYPE_INT     : EXPR_ASSIGN2(char, int)     ; break;
    case ARRAY_TYPE_LONG    : EXPR_ASSIGN2(char, long)    ; break;
    case ARRAY_TYPE_SIZE_T  : EXPR_ASSIGN2(char, size_t)  ; break;
    case ARRAY_TYPE_FLOAT   : EXPR_ASSIGN2(char, float)   ; break;
    case ARRAY_TYPE_DOUBLE  : EXPR_ASSIGN2(char, double)  ; break;
    case ARRAY_TYPE_COMPLEX : EXPR_ASSIGN2(char, Complex) ; break;
    }
    break;

  case ARRAY_TYPE_SHORT:
    switch (B->dtype) {
    case ARRAY_TYPE_BOOL    : EXPR_ASSIGN2(short, Bool)    ; break;
    case ARRAY_TYPE_CHAR    : EXPR_ASSIGN2(short, char)    ; break;
    case ARRAY_TYPE_SHORT   : EXPR_ASSIGN2(short, short)   ; break;
    case ARRAY_TYPE_INT     : EXPR_ASSIGN2(short, int)     ; break;
    case ARRAY_TYPE_LONG    : EXPR_ASSIGN2(short, long)    ; break;
    case ARRAY_TYPE_SIZE_T  : EXPR_ASSIGN2(short, size_t)  ; break;
    case ARRAY_TYPE_FLOAT   : EXPR_ASSIGN2(short, float)   ; break;
    case ARRAY_TYPE_DOUBLE  : EXPR_ASSIGN2(short, double)  ; break;
    case ARRAY_TYPE_COMPLEX : EXPR_ASSIGN2(short, Complex) ; break;
    }
    break;

  case ARRAY_TYPE_INT:
    switch (B->dtype) {
    case ARRAY_TYPE_BOOL    : EXPR_ASSIGN2(int, Bool)    ; break;
    case ARRAY_TYPE_CHAR    : EXPR_ASSIGN2(int, char)    ; break;
    case ARRAY_TYPE_SHORT   : EXPR_ASSIGN2(int, short)   ; break;
    case ARRAY_TYPE_INT     : EXPR_ASSIGN2(int, int)     ; break;
    case ARRAY_TYPE_LONG    : EXPR_ASSIGN2(int, long)    ; break;
    case ARRAY_TYPE_SIZE_T  : EXPR_ASSIGN2(int, size_t)  ; break;
    case ARRAY_TYPE_FLOAT   : EXPR_ASSIGN2(int, float)   ; break;
    case ARRAY_TYPE_DOUBLE  : EXPR_ASSIGN2(int, double)  ; break;
    case ARRAY_TYPE_COMPLEX : EXPR_ASSIGN2(int, Complex) ; break;
    }
    break;

  case ARRAY_TYPE_LONG:
    switch (B->dtype) {
    case ARRAY_TYPE_BOOL    : EXPR_ASSIGN2(long, Bool)    ; break;
    case ARRAY_TYPE_CHAR    : EXPR_ASSIGN2(long, char)    ; break;
    case ARRAY_TYPE_SHORT   : EXPR_ASSIGN2(long, short)   ; break;
    case ARRAY_TYPE_INT     : EXPR_ASSIGN2(long, int)     ; break;
    case ARRAY_TYPE_LONG    : EXPR_ASSIGN2(long, long)    ; break;
    case ARRAY_TYPE_SIZE_T  : EXPR_ASSIGN2(long, size_t)  ; break;
    case ARRAY_TYPE_FLOAT   : EXPR_ASSIGN2(long, float)   ; break;
    case ARRAY_TYPE_DOUBLE  : EXPR_ASSIGN2(long, double)  ; break;
    case ARRAY_TYPE_COMPLEX : EXPR_ASSIGN2(long, Complex) ; break;
    }
    break;

  case ARRAY_TYPE_SIZE_T:
    switch (B->dtype) {
    case ARRAY_TYPE_BOOL    : EXPR_ASSIGN2(size_t, Bool)    ; break;
    case ARRAY_TYPE_CHAR    : EXPR_ASSIGN2(size_t, char)    ; break;
    case ARRAY_TYPE_SHORT   : EXPR_ASSIGN2(size_t, short)   ; break;
    case ARRAY_TYPE_INT     : EXPR_ASSIGN2(size_t, int)     ; break;
    case ARRAY_TYPE_LONG    : EXPR_ASSIGN2(size_t, long)    ; break;
    case ARRAY_TYPE_SIZE_T  : EXPR_ASSIGN2(size_t, size_t)  ; break;
    case ARRAY_TYPE_FLOAT   : EXPR_ASSIGN2(size_t, float)   ; break;
    case ARRAY_TYPE_DOUBLE  : EXPR_ASSIGN2(size_t, double)  ; break;
    case ARRAY_TYPE_COMPLEX : EXPR_ASSIGN2(size_t, Complex) ; break;
    }
    break;

  case ARRAY_TYPE_FLOAT:
    switch (B->dtype) {
    case ARRAY_TYPE_BOOL    : EXPR_ASSIGN2(float, Bool)    ; break;
    case ARRAY_TYPE_CHAR    : EXPR_ASSIGN2(float, char)    ; break;
    case ARRAY_TYPE_SHORT   : EXPR_ASSIGN2(float, short)   ; break;
    case ARRAY_TYPE_INT     : EXPR_ASSIGN2(float, int)     ; break;
    case ARRAY_TYPE_LONG    : EXPR_ASSIGN2(float, long)    ; break;
    case ARRAY_TYPE_SIZE_T  : EXPR_ASSIGN2(float, size_t)  ; break;
    case ARRAY_TYPE_FLOAT   : EXPR_ASSIGN2(float, float)   ; break;
    case ARRAY_TYPE_DOUBLE  : EXPR_ASSIGN2(float, double)  ; break;
    case ARRAY_TYPE_COMPLEX : EXPR_ASSIGN2(float, Complex) ; break;
    }
    break;

  case ARRAY_TYPE_DOUBLE:
    switch (B->dtype) {
    case ARRAY_TYPE_BOOL    : EXPR_ASSIGN2(double, Bool)    ; break;
    case ARRAY_TYPE_CHAR    : EXPR_ASSIGN2(double, char)    ; break;
    case ARRAY_TYPE_SHORT   : EXPR_ASSIGN2(double, short)   ; break;
    case ARRAY_TYPE_INT     : EXPR_ASSIGN2(double, int)     ; break;
    case ARRAY_TYPE_LONG    : EXPR_ASSIGN2(double, long)    ; break;
    case ARRAY_TYPE_SIZE_T  : EXPR_ASSIGN2(double, size_t)  ; break;
    case ARRAY_TYPE_FLOAT   : EXPR_ASSIGN2(double, float)   ; break;
    case ARRAY_TYPE_DOUBLE  : EXPR_ASSIGN2(double, double)  ; break;
    case ARRAY_TYPE_COMPLEX : EXPR_ASSIGN2(double, Complex) ; break;
    }
    break;

  case ARRAY_TYPE_COMPLEX:
    switch (B->dtype) {
    case ARRAY_TYPE_BOOL    : EXPR_ASSIGN2(Complex, Bool)    ; break;
    case ARRAY_TYPE_CHAR    : EXPR_ASSIGN2(Complex, char)    ; break;
    case ARRAY_TYPE_SHORT   : EXPR_ASSIGN2(Complex, short)   ; break;
    case ARRAY_TYPE_INT     : EXPR_ASSIGN2(Complex, int)     ; break;
    case ARRAY_TYPE_LONG    : EXPR_ASSIGN2(Complex, long)    ; break;
    case ARRAY_TYPE_SIZE_T  : EXPR_ASSIGN2(Complex, size_t)  ; break;
    case ARRAY_TYPE_FLOAT   : EXPR_ASSIGN2(Complex, float)   ; break;
    case ARRAY_TYPE_DOUBLE  : EXPR_ASSIGN2(Complex, double)  ; break;
    case ARRAY_TYPE_COMPLEX : EXPR_ASSIGN2(Complex, Complex) ; break;
    }
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
  for(int i = 0; i < ndims-1; ++i) {\
    ++(indices[i]);\
    if (indices[i] == shape[i]){\
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
  if (ndims < 2) return; /* nothing to do for linear array */

  const void   *restrict a = A->data;
  void         *restrict b = B->data;

  const size_t *restrict shape = A->shape;
  size_t       *restrict indices = (size_t *) calloc(ndims, sizeof(size_t));

  size_t rpos = 0;
  size_t cpos = 0;

  switch (A->dtype) {
    case ARRAY_TYPE_BOOL    : TRANSPOSE(Bool)    ; break;
    case ARRAY_TYPE_CHAR    : TRANSPOSE(char)    ; break;
    case ARRAY_TYPE_SHORT   : TRANSPOSE(short)   ; break;
    case ARRAY_TYPE_INT     : TRANSPOSE(int)     ; break;
    case ARRAY_TYPE_LONG    : TRANSPOSE(long)    ; break;
    case ARRAY_TYPE_SIZE_T  : TRANSPOSE(size_t)  ; break;
    case ARRAY_TYPE_FLOAT   : TRANSPOSE(float)   ; break;
    case ARRAY_TYPE_DOUBLE  : TRANSPOSE(double)  ; break;
    case ARRAY_TYPE_COMPLEX : TRANSPOSE(Complex) ; break;
  }
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
