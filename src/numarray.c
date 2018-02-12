

#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <complex.h>
#include <limits.h>
#include <stdint.h>
#include <float.h>

#include "lunum.h"

#define EXPR_UNM(Ta, Tb) {for(size_t i=0;i<N;++i)((Tb*)b)[i]=-((Ta*)a)[i];}
#define EXPR_BNOT(Ta, Tb) {for(size_t i=0;i<N;++i)((Tb*)b)[i]=~((Ta*)a)[i];}
#define EXPR_ERR_UNARY(Ta, Tb) {luaL_error(L, "Invalid operation"); return;}

#define EXPR_AN_ADD(Ta, Tb, Tc) {\
  if (array_first) {\
    for(size_t i=0;i<N;++i)((Tc*)c)[i]=((Ta*)a)[i]+((Tb*)b)[0];\
  } else {\
    for(size_t i=0;i<N;++i)((Tc*)c)[i]=((Tb*)b)[0]+((Ta*)a)[i];\
  }\
}
#define EXPR_AN_SUB(Ta, Tb, Tc) {\
  if (array_first) {\
    for(size_t i=0;i<N;++i)((Tc*)c)[i]=((Ta*)a)[i]-((Tb*)b)[0];\
  } else {\
    for(size_t i=0;i<N;++i)((Tc*)c)[i]=((Tb*)b)[0]-((Ta*)a)[i];\
  }\
}
#define EXPR_AN_MUL(Ta, Tb, Tc) {\
  if (array_first) {\
    for(size_t i=0;i<N;++i)((Tc*)c)[i]=((Ta*)a)[i]*((Tb*)b)[0];\
  } else {\
    for(size_t i=0;i<N;++i)((Tc*)c)[i]=((Tb*)b)[0]*((Ta*)a)[i];\
  }\
}
#define EXPR_AN_DIV(Ta, Tb, Tc) {\
  if (array_first) {\
    for(size_t i=0;i<N;++i)((Tc*)c)[i]=((Ta*)a)[i]/((Tb*)b)[0];\
  } else {\
    for(size_t i=0;i<N;++i)((Tc*)c)[i]=((Tb*)b)[0]/((Ta*)a)[i];\
  }\
}
#define EXPR_AN_FIDIVF(Ta, Tb, Tc) {\
  if (array_first) {\
    for(size_t i=0;i<N;++i)((Tc*)c)[i]=(Tc)floorf(((Ta*)a)[i]/((Tb*)b)[0]);\
  } else {\
    for(size_t i=0;i<N;++i)((Tc*)c)[i]=(Tc)floorf(((Tb*)b)[0]/((Ta*)a)[i]);\
  }\
}
#define EXPR_AN_FIDIV(Ta, Tb, Tc) {\
  if (array_first) {\
    for(size_t i=0;i<N;++i)((Tc*)c)[i]=(Tc)floor(((Ta*)a)[i]/((Tb*)b)[0]);\
  } else {\
    for(size_t i=0;i<N;++i)((Tc*)c)[i]=(Tc)floor(((Tb*)b)[0]/((Ta*)a)[i]);\
  }\
}
#define EXPR_AN_MOD(Ta, Tb, Tc) {\
  if (array_first) {\
    for(size_t i=0;i<N;++i)((Tc*)c)[i]=((Ta*)a)[i]-(Tc)(((Ta*)a)[i]/((Tb*)b)[0])*((Tb*)b)[0];\
  } else {\
    for(size_t i=0;i<N;++i)((Tc*)c)[i]=((Tb*)b)[0]-(Tc)(((Tb*)b)[0]/((Ta*)a)[i])*((Ta*)a)[i];\
  }\
}
#define EXPR_AN_FMODF(Ta, Tb, Tc) {\
  if (array_first) {\
    for(size_t i=0;i<N;++i)((Tc*)c)[i]=((Ta*)a)[i]-(Tc)floorf(((Ta*)a)[i]/((Tb*)b)[0])*((Tb*)b)[0];\
  } else {\
    for(size_t i=0;i<N;++i)((Tc*)c)[i]=((Tb*)b)[0]-(Tc)floorf(((Tb*)b)[0]/((Ta*)a)[i])*((Ta*)a)[i];\
  }\
}
#define EXPR_AN_FMOD(Ta, Tb, Tc) {\
  if (array_first) {\
    for(size_t i=0;i<N;++i)((Tc*)c)[i]=((Ta*)a)[i]-(Tc)floor(((Ta*)a)[i]/((Tb*)b)[0])*((Tb*)b)[0];\
  } else {\
    for(size_t i=0;i<N;++i)((Tc*)c)[i]=((Tb*)b)[0]-(Tc)floor(((Tb*)b)[0]/((Ta*)a)[i])*((Ta*)a)[i];\
  }\
}
#define EXPR_AN_POW(Ta, Tb, Tc) {\
  if (array_first) {\
    for(size_t i=0;i<N;++i)((Tc*)c)[i]=(Tc)pow(((Ta*)a)[i],((Tb*)b)[0]);\
  } else {\
    for(size_t i=0;i<N;++i)((Tc*)c)[i]=(Tc)pow(((Tb*)b)[0],((Ta*)a)[i]);\
  }\
}
#define EXPR_AN_FPOWF(Ta, Tb, Tc) {\
  if (array_first) {\
    for(size_t i=0;i<N;++i)((Tc*)c)[i]=(Tc)powf(((Ta*)a)[i],((Tb*)b)[0]);\
  } else {\
    for(size_t i=0;i<N;++i)((Tc*)c)[i]=(Tc)powf(((Tb*)b)[0],((Ta*)a)[i]);\
  }\
}
#define EXPR_AN_FPOW(Ta, Tb, Tc) {\
  if (array_first) {\
    for(size_t i=0;i<N;++i)((Tc*)c)[i]=(Tc)pow(((Ta*)a)[i],((Tb*)b)[0]);\
  } else {\
    for(size_t i=0;i<N;++i)((Tc*)c)[i]=(Tc)pow(((Tb*)b)[0],((Ta*)a)[i]);\
  }\
}
#define EXPR_AN_CPOW(Ta, Tb, Tc) {\
  if (array_first) {\
    for(size_t i=0;i<N;++i)((Tc*)c)[i]=(Tc)cpow(((Ta*)a)[i],((Tb*)b)[0]);\
  } else {\
    for(size_t i=0;i<N;++i)((Tc*)c)[i]=(Tc)cpow(((Tb*)b)[0],((Ta*)a)[i]);\
  }\
}
#define EXPR_AN_BAND(Ta, Tb, Tc) {\
  if (array_first) {\
    for(size_t i=0;i<N;++i)((Tc*)c)[i]=((Ta*)a)[i]&((Tb*)b)[0];\
  } else {\
    for(size_t i=0;i<N;++i)((Tc*)c)[i]=((Tb*)b)[0]&((Ta*)a)[i];\
  }\
}
#define EXPR_AN_BOR(Ta, Tb, Tc) {\
  if (array_first) {\
    for(size_t i=0;i<N;++i)((Tc*)c)[i]=((Ta*)a)[i]|((Tb*)b)[0];\
  } else {\
    for(size_t i=0;i<N;++i)((Tc*)c)[i]=((Tb*)b)[0]|((Ta*)a)[i];\
  }\
}
#define EXPR_AN_BXOR(Ta, Tb, Tc) {\
  if (array_first) {\
    for(size_t i=0;i<N;++i)((Tc*)c)[i]=((Ta*)a)[i]^((Tb*)b)[0];\
  } else {\
    for(size_t i=0;i<N;++i)((Tc*)c)[i]=((Tb*)b)[0]^((Ta*)a)[i];\
  }\
}
#define EXPR_AN_SHL(Ta, Tb, Tc) {\
  if (array_first) {\
    for(size_t i=0;i<N;++i)((Tc*)c)[i]=((Ta*)a)[i]<<((Tb*)b)[0];\
  } else {\
    for(size_t i=0;i<N;++i)((Tc*)c)[i]=((Tb*)b)[0]<<((Ta*)a)[i];\
  }\
}
#define EXPR_AN_SHR(Ta, Tb, Tc) {\
  if (array_first) {\
    for(size_t i=0;i<N;++i)((Tc*)c)[i]=((Ta*)a)[i]>>((Tb*)b)[0];\
  } else {\
    for(size_t i=0;i<N;++i)((Tc*)c)[i]=((Tb*)b)[0]>>((Ta*)a)[i];\
  }\
}

#define EXPR_AA_ADD(Ta, Tb, Tc) {for(size_t i=0;i<N;++i)((Tc*)c)[i]=((Ta*)a)[i]+((Tb*)b)[i];}
#define EXPR_AA_SUB(Ta, Tb, Tc) {for(size_t i=0;i<N;++i)((Tc*)c)[i]=((Ta*)a)[i]-((Tb*)b)[i];}
#define EXPR_AA_MUL(Ta, Tb, Tc) {for(size_t i=0;i<N;++i)((Tc*)c)[i]=((Ta*)a)[i]*((Tb*)b)[i];}
#define EXPR_AA_DIV(Ta, Tb, Tc) {for(size_t i=0;i<N;++i)((Tc*)c)[i]=((Ta*)a)[i]/((Tb*)b)[i];}
#define EXPR_AA_FIDIVF(Ta, Tb, Tc) {for(size_t i=0;i<N;++i)((Tc*)c)[i]=(Tc)floorf(((Ta*)a)[i]/((Tb*)b)[i]);}
#define EXPR_AA_FIDIV(Ta, Tb, Tc) {for(size_t i=0;i<N;++i)((Tc*)c)[i]=(Tc)floor(((Ta*)a)[i]/((Tb*)b)[i]);}
#define EXPR_AA_MOD(Ta, Tb, Tc) {for(size_t i=0;i<N;++i)((Tc*)c)[i]=((Ta*)a)[i]-(Tc)(((Ta*)a)[i]/((Tb*)b)[i])*((Tb*)b)[i];}
#define EXPR_AA_FMODF(Ta, Tb, Tc) {for(size_t i=0;i<N;++i)((Tc*)c)[i]=((Ta*)a)[i]-(Tc)floorf(((Ta*)a)[i]/((Tb*)b)[i])*((Tb*)b)[i];}
#define EXPR_AA_FMOD(Ta, Tb, Tc) {for(size_t i=0;i<N;++i)((Tc*)c)[i]=((Ta*)a)[i]-(Tc)floor(((Ta*)a)[i]/((Tb*)b)[i])*((Tb*)b)[i];}
#define EXPR_AA_POW(Ta, Tb, Tc) {for(size_t i=0;i<N;++i)((Tc*)c)[i]= pow(((Ta*)a)[i],((Tb*)b)[i]);}
#define EXPR_AA_FPOWF(Ta, Tb, Tc) {for(size_t i=0;i<N;++i)((Tc*)c)[i]=powf(((Ta*)a)[i],((Tb*)b)[i]);}
#define EXPR_AA_FPOW(Ta, Tb, Tc) {for(size_t i=0;i<N;++i)((Tc*)c)[i]= pow(((Ta*)a)[i],((Tb*)b)[i]);}
#define EXPR_AA_CPOW(Ta, Tb, Tc) {for(size_t i=0;i<N;++i)((Tc*)c)[i]=cpow(((Ta*)a)[i],((Tb*)b)[i]);}
#define EXPR_AA_BAND(Ta, Tb, Tc) {for(size_t i=0;i<N;++i)((Tc*)c)[i]=((Ta*)a)[i]&((Tb*)b)[i];}
#define EXPR_AA_BOR(Ta, Tb, Tc) {for(size_t i=0;i<N;++i)((Tc*)c)[i]=((Ta*)a)[i]|((Tb*)b)[i];}
#define EXPR_AA_BXOR(Ta, Tb, Tc) {for(size_t i=0;i<N;++i)((Tc*)c)[i]=((Ta*)a)[i]^((Tb*)b)[i];}
#define EXPR_AA_SHL(Ta, Tb, Tc) {for(size_t i=0;i<N;++i)((Tc*)c)[i]=((Ta*)a)[i]<<((Tb*)b)[i];}
#define EXPR_AA_SHR(Ta, Tb, Tc) {for(size_t i=0;i<N;++i)((Tc*)c)[i]=((Ta*)a)[i]>>((Tb*)b)[i];}
#define EXPR_ERR_BINARY(Ta, Tb, Tc) {luaL_error(L, "Invalid operation"); return;}

#define GET_TC(op, Ta, Tb) {\
  if (outA) {\
    op(Ta, Tb, Ta)\
  } else {\
    op(Ta, Tb, Tb)\
  }\
}

#define ARRAY_OP_BINARY_B(Ta, op, opf, opd, opc) \
    switch (TypeB) {\
      case ARRAY_TYPE_BOOL    : GET_TC(op,  Ta, Bool   ) ; break;\
      case ARRAY_TYPE_CHAR    : GET_TC(op,  Ta, char   ) ; break;\
      case ARRAY_TYPE_SHORT   : GET_TC(op,  Ta, short  ) ; break;\
      case ARRAY_TYPE_INT     : GET_TC(op,  Ta, int    ) ; break;\
      case ARRAY_TYPE_LONG    : GET_TC(op,  Ta, long   ) ; break;\
      case ARRAY_TYPE_SIZE_T  : GET_TC(op,  Ta, size_t ) ; break;\
      case ARRAY_TYPE_FLOAT   : GET_TC(opf, Ta, float  ) ; break;\
      case ARRAY_TYPE_DOUBLE  : GET_TC(opd, Ta, double ) ; break;\
      case ARRAY_TYPE_COMPLEX : GET_TC(opc, Ta, Complex) ; break;\
    }

#define ARRAY_OP_BINARY(op, opf, opd, opc) \
    switch (TypeA) {\
      case ARRAY_TYPE_BOOL    : ARRAY_OP_BINARY_B(Bool,    op,  opf, opd, opc) ; break;\
      case ARRAY_TYPE_CHAR    : ARRAY_OP_BINARY_B(char,    op,  opf, opd, opc) ; break;\
      case ARRAY_TYPE_SHORT   : ARRAY_OP_BINARY_B(short,   op,  opf, opd, opc) ; break;\
      case ARRAY_TYPE_INT     : ARRAY_OP_BINARY_B(int,     op,  opf, opd, opc) ; break;\
      case ARRAY_TYPE_LONG    : ARRAY_OP_BINARY_B(long,    op,  opf, opd, opc) ; break;\
      case ARRAY_TYPE_SIZE_T  : ARRAY_OP_BINARY_B(size_t,  op,  opf, opd, opc) ; break;\
      case ARRAY_TYPE_FLOAT   : ARRAY_OP_BINARY_B(float,   opf, opf, opd, opc) ; break;\
      case ARRAY_TYPE_DOUBLE  : ARRAY_OP_BINARY_B(double,  opd, opd, opd, opc) ; break;\
      case ARRAY_TYPE_COMPLEX : ARRAY_OP_BINARY_B(Complex, opc, opc, opc, opc) ; break;\
    }

#define ARRAY_OP_UNARY(op, opf, opd, opc) \
    switch (TypeA) {\
      case ARRAY_TYPE_BOOL    : op(Bool,     Bool   ) ; break;\
      case ARRAY_TYPE_CHAR    : op(char,     char   ) ; break;\
      case ARRAY_TYPE_SHORT   : op(short,    short  ) ; break;\
      case ARRAY_TYPE_INT     : op(int,      int    ) ; break;\
      case ARRAY_TYPE_LONG    : op(long,     long   ) ; break;\
      case ARRAY_TYPE_SIZE_T  : op(size_t,   size_t ) ; break;\
      case ARRAY_TYPE_FLOAT   : opf(float,   float  ) ; break;\
      case ARRAY_TYPE_DOUBLE  : opd(double,  double ) ; break;\
      case ARRAY_TYPE_COMPLEX : opc(Complex, Complex) ; break;\
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

double array_typemin(ArrayType T)
{
  switch (T) {
  case ARRAY_TYPE_BOOL    : return 0;
  case ARRAY_TYPE_CHAR    : return CHAR_MIN;
  case ARRAY_TYPE_SHORT   : return SHRT_MIN;
  case ARRAY_TYPE_INT     : return INT_MIN;
  case ARRAY_TYPE_LONG    : return LONG_MIN;
  case ARRAY_TYPE_SIZE_T  : return 0;
  case ARRAY_TYPE_FLOAT   : return -FLT_MAX;
  case ARRAY_TYPE_DOUBLE  : return -DBL_MAX;
  case ARRAY_TYPE_COMPLEX : return -DBL_MAX;
  }
  return -1; // indicates invalid type
}

double array_typemax(ArrayType T)
{
  switch (T) {
  case ARRAY_TYPE_BOOL    : return 1;
  case ARRAY_TYPE_CHAR    : return CHAR_MAX;
  case ARRAY_TYPE_SHORT   : return SHRT_MAX;
  case ARRAY_TYPE_INT     : return INT_MAX;
  case ARRAY_TYPE_LONG    : return LONG_MAX;
  case ARRAY_TYPE_SIZE_T  : return SIZE_MAX;
  case ARRAY_TYPE_FLOAT   : return FLT_MAX;
  case ARRAY_TYPE_DOUBLE  : return DBL_MAX;
  case ARRAY_TYPE_COMPLEX : return DBL_MAX;
  }
  return -1; // indicates invalid type
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
  array_resize_t(&A, B->shape, B->ndims);
  array_assign_from_array(&A, B);
  return A;
}

void array_del(Array *A)
{
  if (A->data) free(A->data);
  if (A->shape) free(A->shape);

  A->size = 0;
  A->data = NULL;
  A->shape = NULL;
}

int array_resize(Array *A, const int *N, int Nd)
{
  size_t ntot = 1;
  for (int d=0; d<Nd; ++d) ntot *= N[d];

  if (A->size != ntot) {
    return 1;
  }
  if (A->ndims != Nd) {
    if (A->shape) {free(A->shape); A->shape = NULL;}
    A->ndims = Nd;
  }
  if (A->shape == NULL) {
    A->shape = (size_t*) malloc(Nd*sizeof(size_t));
  }

  for (int d=0; d<Nd; ++d) A->shape[d] = N[d];

  return 0;
}

int array_resize_t(Array *A, const size_t *N, int Nd)
{
  size_t ntot = 1;
  for (int d=0; d<Nd; ++d) ntot *= N[d];

  if (A->size != ntot) {
    return 1;
  }
  if (A->ndims != Nd) {
    if (A->shape) {free(A->shape); A->shape = NULL;}
    A->ndims = Nd;
  }
  if (A->shape == NULL) {
    A->shape = (size_t*) malloc(Nd*sizeof(size_t));
  }

  memcpy(A->shape, N, Nd*sizeof(size_t));

  return 0;
}

/* Arrays A and B guaranteed to have same types */
void array_unary_op(lua_State * L, const Array *A, Array *B,
                    ArrayUnaryOperation op)
{
  const size_t          N = A->size;
  const void  *restrict a = A->data;
  void        *restrict b = B->data;

  const ArrayType TypeA = A->dtype;

  switch (op) {
  case ARRAY_OP_UNM:
    ARRAY_OP_UNARY(EXPR_UNM, EXPR_UNM, EXPR_UNM, EXPR_UNM);
    break;
  case ARRAY_OP_BNOT:
    ARRAY_OP_UNARY(EXPR_BNOT, EXPR_ERR_UNARY, EXPR_ERR_UNARY, EXPR_ERR_UNARY);
    break;
  }
}

/* Arrays A and B may have different types, C has higher of A and B types */
void array_array_binary_op(lua_State * L, const Array *A, const Array *B,
                     Array *C, ArrayBinaryOperation op)
{
  const size_t          N = A->size;
  const void  *restrict a = A->data;
  const void  *restrict b = B->data;
  void        *restrict c = C->data;

  const ArrayType TypeA = A->dtype;
  const ArrayType TypeB = B->dtype;
  const Bool      outA  = TypeA == C->dtype;

  switch (op) {
    case ARRAY_OP_ADD:
      ARRAY_OP_BINARY(EXPR_AA_ADD, EXPR_AA_ADD, EXPR_AA_ADD, EXPR_AA_ADD);
      break;
    case ARRAY_OP_SUB:
      ARRAY_OP_BINARY(EXPR_AA_SUB, EXPR_AA_SUB, EXPR_AA_SUB, EXPR_AA_SUB);
      break;
    case ARRAY_OP_MUL:
      ARRAY_OP_BINARY(EXPR_AA_MUL, EXPR_AA_MUL, EXPR_AA_MUL, EXPR_AA_MUL);
      break;
    case ARRAY_OP_DIV:
      ARRAY_OP_BINARY(EXPR_AA_DIV, EXPR_AA_DIV, EXPR_AA_DIV, EXPR_AA_DIV);
      break;
    case ARRAY_OP_IDIV:
      ARRAY_OP_BINARY(EXPR_AA_DIV, EXPR_AA_FIDIVF, EXPR_AA_FIDIV, EXPR_ERR_BINARY);
      break;
    case ARRAY_OP_MOD:
      ARRAY_OP_BINARY(EXPR_AA_MOD, EXPR_AA_FMODF, EXPR_AA_FMOD, EXPR_ERR_BINARY);
      break;
    case ARRAY_OP_POW:
      ARRAY_OP_BINARY(EXPR_AA_POW, EXPR_AA_FPOWF, EXPR_AA_FPOW, EXPR_AA_CPOW);
      break;
    case ARRAY_OP_BAND:
      ARRAY_OP_BINARY(EXPR_AA_BAND, EXPR_ERR_BINARY, EXPR_ERR_BINARY, EXPR_ERR_BINARY);
      break;
    case ARRAY_OP_BOR:
      ARRAY_OP_BINARY(EXPR_AA_BOR, EXPR_ERR_BINARY, EXPR_ERR_BINARY, EXPR_ERR_BINARY);
      break;
    case ARRAY_OP_BXOR:
      ARRAY_OP_BINARY(EXPR_AA_BXOR, EXPR_ERR_BINARY, EXPR_ERR_BINARY, EXPR_ERR_BINARY);
      break;
    case ARRAY_OP_SHL:
      ARRAY_OP_BINARY(EXPR_AA_SHL, EXPR_ERR_BINARY, EXPR_ERR_BINARY, EXPR_ERR_BINARY);
      break;
    case ARRAY_OP_SHR:
      ARRAY_OP_BINARY(EXPR_AA_SHR, EXPR_ERR_BINARY, EXPR_ERR_BINARY, EXPR_ERR_BINARY);
      break;
  }
}

/* Arrays A and B may have different types, C is the same type as B */
void array_number_binary_op(lua_State * L, const Array *A, const void *B,
                     Array *C, ArrayBinaryOperation op, Bool array_first)
{
  const size_t          N = A->size;
  const void  *restrict a = A->data;
  const void  *restrict b = B;
  void        *restrict c = C->data;

  const ArrayType TypeA = A->dtype;
  const ArrayType TypeB = C->dtype;
  const Bool      outA  = TypeA == C->dtype;

  switch (op) {
    case ARRAY_OP_ADD:
      ARRAY_OP_BINARY(EXPR_AN_ADD, EXPR_AN_ADD, EXPR_AN_ADD, EXPR_AN_ADD);
      break;
    case ARRAY_OP_SUB:
      ARRAY_OP_BINARY(EXPR_AN_SUB, EXPR_AN_SUB, EXPR_AN_SUB, EXPR_AN_SUB);
      break;
    case ARRAY_OP_MUL:
      ARRAY_OP_BINARY(EXPR_AN_MUL, EXPR_AN_MUL, EXPR_AN_MUL, EXPR_AN_MUL);
      break;
    case ARRAY_OP_DIV:
      ARRAY_OP_BINARY(EXPR_AN_DIV, EXPR_AN_DIV, EXPR_AN_DIV, EXPR_AN_DIV);
      break;
    case ARRAY_OP_IDIV:
      ARRAY_OP_BINARY(EXPR_AN_DIV, EXPR_AN_FIDIVF, EXPR_AN_FIDIV, EXPR_ERR_BINARY);
      break;
    case ARRAY_OP_MOD:
      ARRAY_OP_BINARY(EXPR_AN_MOD, EXPR_AN_FMODF, EXPR_AN_FMOD, EXPR_ERR_BINARY);
      break;
    case ARRAY_OP_POW:
      ARRAY_OP_BINARY(EXPR_AN_POW, EXPR_AN_FPOWF, EXPR_AN_FPOW, EXPR_AA_CPOW);
      break;
    case ARRAY_OP_BAND:
      ARRAY_OP_BINARY(EXPR_AN_BAND, EXPR_ERR_BINARY, EXPR_ERR_BINARY, EXPR_ERR_BINARY);
      break;
    case ARRAY_OP_BOR:
      ARRAY_OP_BINARY(EXPR_AN_BOR, EXPR_ERR_BINARY, EXPR_ERR_BINARY, EXPR_ERR_BINARY);
      break;
    case ARRAY_OP_BXOR:
      ARRAY_OP_BINARY(EXPR_AN_BXOR, EXPR_ERR_BINARY, EXPR_ERR_BINARY, EXPR_ERR_BINARY);
      break;
    case ARRAY_OP_SHL:
      ARRAY_OP_BINARY(EXPR_AN_SHL, EXPR_ERR_BINARY, EXPR_ERR_BINARY, EXPR_ERR_BINARY);
      break;
    case ARRAY_OP_SHR:
      ARRAY_OP_BINARY(EXPR_AN_SHR, EXPR_ERR_BINARY, EXPR_ERR_BINARY, EXPR_ERR_BINARY);
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

#define TRANSPOSE(Ta, Tb) \
{\
  rpos = 0;\
  while (indices[0] < shape[0])\
  {\
    POS_COLMAJOR(cpos, indices, shape, ndims);\
    ((Tb*)b)[cpos]=((Ta*)a)[rpos];\
    ADVANCE_INDICES(indices, shape, ndims);\
    ++rpos;\
  }\
}

void array_transpose(const Array *A, Array *B)
{
  const int        ndims = A->ndims;
  const void *restrict a = A->data;
  void       *restrict b = B->data;

  if (ndims < 2) {
    /* copy linear array */
    memcpy(b, a, A->size * array_sizeof(A->dtype));
    return;
  }

  const size_t *restrict   shape = A->shape;
  const ArrayType          TypeA = A->dtype;
  size_t       *restrict indices = (size_t *) calloc(ndims, sizeof(size_t));

  size_t rpos = 0;
  size_t cpos = 0;

  ARRAY_OP_UNARY(TRANSPOSE, TRANSPOSE, TRANSPOSE, TRANSPOSE);

  for (int d = 0; d < ndims; d++) {
    indices[d] = shape[ndims - d - 1];
  }
  array_resize_t(B, indices, ndims);

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
  array_resize_t(&B0, N, Nd);
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
