

#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <complex.h>
#include <limits.h>
#include <stdint.h>
#include <float.h>

#include "lunum.h"

#define EXPR_UN_MNS(Ta, Tb)  {const Ta*restrict _a=a;Tb*restrict _b=b;for(size_t i=0;i<N;++i)_b[i]=-_a[i];}
#define EXPR_UN_BNOT(Ta, Tb) {const Ta*restrict _a=a;Tb*restrict _b=b;for(size_t i=0;i<N;++i)_b[i]=~_a[i];}
#define EXPR_UN_ERR(Ta, Tb)  {luaL_error(L, "Invalid operation"); return;}

#define EXPR_BIN_ADD(Ta, a, ia, Tb, b, ib, Tc, c)    {Tc*restrict _c=c;const Ta*restrict _a=a;const Tb*restrict _b=b;for(size_t i=0;i<N;++i)_c[i]=(Tc)_a[ia]+(Tc)_b[ib];}
#define EXPR_BIN_SUB(Ta, a, ia, Tb, b, ib, Tc, c)    {Tc*restrict _c=c;const Ta*restrict _a=a;const Tb*restrict _b=b;for(size_t i=0;i<N;++i)_c[i]=(Tc)_a[ia]-(Tc)_b[ib];}
#define EXPR_BIN_MUL(Ta, a, ia, Tb, b, ib, Tc, c)    {Tc*restrict _c=c;const Ta*restrict _a=a;const Tb*restrict _b=b;for(size_t i=0;i<N;++i)_c[i]=(Tc)_a[ia]*(Tc)_b[ib];}
#define EXPR_BIN_DIV(Ta, a, ia, Tb, b, ib, Tc, c)    {Tc*restrict _c=c;const Ta*restrict _a=a;const Tb*restrict _b=b;for(size_t i=0;i<N;++i)_c[i]=(Tc)_a[ia]/(Tc)_b[ib];}
#define EXPR_BIN_FIDIVF(Ta, a, ia, Tb, b, ib, Tc, c) {Tc*restrict _c=c;const Ta*restrict _a=a;const Tb*restrict _b=b;for(size_t i=0;i<N;++i)_c[i]=(Tc)floorf((Tc)_a[ia]/(Tc)_b[ib]);}
#define EXPR_BIN_FIDIV(Ta, a, ia, Tb, b, ib, Tc, c)  {Tc*restrict _c=c;const Ta*restrict _a=a;const Tb*restrict _b=b;for(size_t i=0;i<N;++i)_c[i]=(Tc) floor((Tc)_a[ia]/(Tc)_b[ib]);}
#define EXPR_BIN_MOD(Ta, a, ia, Tb, b, ib, Tc, c)    {Tc*restrict _c=c;const Ta*restrict _a=a;const Tb*restrict _b=b;for(size_t i=0;i<N;++i)_c[i]=(Tc)_a[ia]-(Tc)      ((Tc)_a[ia]/(Tc)_b[ib])*(Tc)_b[ib];}
#define EXPR_BIN_FMODF(Ta, a, ia, Tb, b, ib, Tc, c)  {Tc*restrict _c=c;const Ta*restrict _a=a;const Tb*restrict _b=b;for(size_t i=0;i<N;++i)_c[i]=(Tc)_a[ia]-(Tc)floorf((Tc)_a[ia]/(Tc)_b[ib])*(Tc)_b[ib];}
#define EXPR_BIN_FMOD(Ta, a, ia, Tb, b, ib, Tc, c)   {Tc*restrict _c=c;const Ta*restrict _a=a;const Tb*restrict _b=b;for(size_t i=0;i<N;++i)_c[i]=(Tc)_a[ia]-(Tc) floor((Tc)_a[ia]/(Tc)_b[ib])*(Tc)_b[ib];}
#define EXPR_BIN_POW(Ta, a, ia, Tb, b, ib, Tc, c)    {Tc*restrict _c=c;const Ta*restrict _a=a;const Tb*restrict _b=b;for(size_t i=0;i<N;++i)_c[i]= pow((Tc)_a[ia],(Tc)_b[ib]);}
#define EXPR_BIN_FPOWF(Ta, a, ia, Tb, b, ib, Tc, c)  {Tc*restrict _c=c;const Ta*restrict _a=a;const Tb*restrict _b=b;for(size_t i=0;i<N;++i)_c[i]=powf((Tc)_a[ia],(Tc)_b[ib]);}
#define EXPR_BIN_FPOW(Ta, a, ia, Tb, b, ib, Tc, c)   {Tc*restrict _c=c;const Ta*restrict _a=a;const Tb*restrict _b=b;for(size_t i=0;i<N;++i)_c[i]= pow((Tc)_a[ia],(Tc)_b[ib]);}
#define EXPR_BIN_CPOW(Ta, a, ia, Tb, b, ib, Tc, c)   {Tc*restrict _c=c;const Ta*restrict _a=a;const Tb*restrict _b=b;for(size_t i=0;i<N;++i)_c[i]=cpow((Tc)_a[ia],(Tc)_b[ib]);}
#define EXPR_BIN_BAND(Ta, a, ia, Tb, b, ib, Tc, c)   {Tc*restrict _c=c;const Ta*restrict _a=a;const Tb*restrict _b=b;for(size_t i=0;i<N;++i)_c[i]=(Tc)_a[ia]&(Tc)_b[ib];}
#define EXPR_BIN_BOR(Ta, a, ia, Tb, b, ib, Tc, c)    {Tc*restrict _c=c;const Ta*restrict _a=a;const Tb*restrict _b=b;for(size_t i=0;i<N;++i)_c[i]=(Tc)_a[ia]|(Tc)_b[ib];}
#define EXPR_BIN_BXOR(Ta, a, ia, Tb, b, ib, Tc, c)   {Tc*restrict _c=c;const Ta*restrict _a=a;const Tb*restrict _b=b;for(size_t i=0;i<N;++i)_c[i]=(Tc)_a[ia]^(Tc)_b[ib];}
#define EXPR_BIN_SHL(Ta, a, ia, Tb, b, ib, Tc, c)    {Tc*restrict _c=c;const Ta*restrict _a=a;const Tb*restrict _b=b;for(size_t i=0;i<N;++i)_c[i]=(Tc)_a[ia]<<(Tc)_b[ib];}
#define EXPR_BIN_SHR(Ta, a, ia, Tb, b, ib, Tc, c)    {Tc*restrict _c=c;const Ta*restrict _a=a;const Tb*restrict _b=b;for(size_t i=0;i<N;++i)_c[i]=(Tc)_a[ia]>>(Tc)_b[ib];}
#define EXPR_BIN_ERR(Ta, a, ia, Tb, b, ib, Tc, c)    {luaL_error(L, "Invalid operation"); return;}

#define GET_ORDR(op, Ta, Tb) {\
  if (array_first) {\
    op(Ta, a, i, Tb, b, 0, Tb, c)\
  } else {\
    op(Tb, b, 0, Ta, a, i, Tb, c)\
  }\
}

#define GET_TYPEC(op, Ta, Tb) {\
  if (outA) {\
    op(Ta, a, i, Tb, b, i, Ta, c)\
  } else {\
    op(Ta, a, i, Tb, b, i, Tb, c)\
  }\
}

#define ARRAY_OP_BINARY_B(op_adpt, Ta, op, opf, opd, opc) \
    switch (TypeB) {\
      case ARRAY_TYPE_BOOL    : op_adpt(op,  Ta, Bool   ) ; break;\
      case ARRAY_TYPE_CHAR    : op_adpt(op,  Ta, char   ) ; break;\
      case ARRAY_TYPE_SHORT   : op_adpt(op,  Ta, short  ) ; break;\
      case ARRAY_TYPE_INT     : op_adpt(op,  Ta, int    ) ; break;\
      case ARRAY_TYPE_LONG    : op_adpt(op,  Ta, long   ) ; break;\
      case ARRAY_TYPE_SIZE_T  : op_adpt(op,  Ta, size_t ) ; break;\
      case ARRAY_TYPE_FLOAT   : op_adpt(opf, Ta, float  ) ; break;\
      case ARRAY_TYPE_DOUBLE  : op_adpt(opd, Ta, double ) ; break;\
      case ARRAY_TYPE_COMPLEX : op_adpt(opc, Ta, Complex) ; break;\
    }

#define ARRAY_OP_BINARY(op_adpt, op, opf, opd, opc) \
    switch (TypeA) {\
      case ARRAY_TYPE_BOOL    : ARRAY_OP_BINARY_B(op_adpt, Bool,    op,  opf, opd, opc) ; break;\
      case ARRAY_TYPE_CHAR    : ARRAY_OP_BINARY_B(op_adpt, char,    op,  opf, opd, opc) ; break;\
      case ARRAY_TYPE_SHORT   : ARRAY_OP_BINARY_B(op_adpt, short,   op,  opf, opd, opc) ; break;\
      case ARRAY_TYPE_INT     : ARRAY_OP_BINARY_B(op_adpt, int,     op,  opf, opd, opc) ; break;\
      case ARRAY_TYPE_LONG    : ARRAY_OP_BINARY_B(op_adpt, long,    op,  opf, opd, opc) ; break;\
      case ARRAY_TYPE_SIZE_T  : ARRAY_OP_BINARY_B(op_adpt, size_t,  op,  opf, opd, opc) ; break;\
      case ARRAY_TYPE_FLOAT   : ARRAY_OP_BINARY_B(op_adpt, float,   opf, opf, opd, opc) ; break;\
      case ARRAY_TYPE_DOUBLE  : ARRAY_OP_BINARY_B(op_adpt, double,  opd, opd, opd, opc) ; break;\
      case ARRAY_TYPE_COMPLEX : ARRAY_OP_BINARY_B(op_adpt, Complex, opc, opc, opc, opc) ; break;\
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


#define EXPR_ASSIGN0(Ta,val) {Ta*restrict _a=a;                         for(size_t i=0;i<N;++i)_a[i]=val;}
#define EXPR_ASSIGN1(Ta,val) {Ta*restrict _a=a;const Ta*restrict _b=val;for(size_t i=0;i<N;++i)_a[i]=_b[0];}
#define EXPR_ASSIGN2(Ta,Tb)  {Ta*restrict _a=a;const Tb*restrict _b=b;  for(size_t i=0;i<N;++i)_a[i]=_b[i];}

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
    ARRAY_OP_UNARY(EXPR_UN_MNS, EXPR_UN_MNS, EXPR_UN_MNS, EXPR_UN_MNS);
    break;
  case ARRAY_OP_BNOT:
    ARRAY_OP_UNARY(EXPR_UN_BNOT, EXPR_UN_ERR, EXPR_UN_ERR, EXPR_UN_ERR);
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
      ARRAY_OP_BINARY(GET_TYPEC, EXPR_BIN_ADD, EXPR_BIN_ADD, EXPR_BIN_ADD, EXPR_BIN_ADD); break;
    case ARRAY_OP_SUB:
      ARRAY_OP_BINARY(GET_TYPEC, EXPR_BIN_SUB, EXPR_BIN_SUB, EXPR_BIN_SUB, EXPR_BIN_SUB); break;
    case ARRAY_OP_MUL:
      ARRAY_OP_BINARY(GET_TYPEC, EXPR_BIN_MUL, EXPR_BIN_MUL, EXPR_BIN_MUL, EXPR_BIN_MUL); break;
    case ARRAY_OP_DIV:
      ARRAY_OP_BINARY(GET_TYPEC, EXPR_BIN_DIV, EXPR_BIN_DIV, EXPR_BIN_DIV, EXPR_BIN_DIV); break;
    case ARRAY_OP_IDIV:
      ARRAY_OP_BINARY(GET_TYPEC, EXPR_BIN_DIV, EXPR_BIN_FIDIVF, EXPR_BIN_FIDIV, EXPR_BIN_ERR); break;
    case ARRAY_OP_MOD:
      ARRAY_OP_BINARY(GET_TYPEC, EXPR_BIN_MOD, EXPR_BIN_FMODF, EXPR_BIN_FMOD, EXPR_BIN_ERR); break;
    case ARRAY_OP_POW:
      ARRAY_OP_BINARY(GET_TYPEC, EXPR_BIN_POW, EXPR_BIN_FPOWF, EXPR_BIN_FPOW, EXPR_BIN_CPOW); break;
    case ARRAY_OP_BAND:
      ARRAY_OP_BINARY(GET_TYPEC, EXPR_BIN_BAND, EXPR_BIN_ERR, EXPR_BIN_ERR, EXPR_BIN_ERR); break;
    case ARRAY_OP_BOR:
      ARRAY_OP_BINARY(GET_TYPEC, EXPR_BIN_BOR, EXPR_BIN_ERR, EXPR_BIN_ERR, EXPR_BIN_ERR); break;
    case ARRAY_OP_BXOR:
      ARRAY_OP_BINARY(GET_TYPEC, EXPR_BIN_BXOR, EXPR_BIN_ERR, EXPR_BIN_ERR, EXPR_BIN_ERR); break;
    case ARRAY_OP_SHL:
      ARRAY_OP_BINARY(GET_TYPEC, EXPR_BIN_SHL, EXPR_BIN_ERR, EXPR_BIN_ERR, EXPR_BIN_ERR); break;
    case ARRAY_OP_SHR:
      ARRAY_OP_BINARY(GET_TYPEC, EXPR_BIN_SHR, EXPR_BIN_ERR, EXPR_BIN_ERR, EXPR_BIN_ERR); break;
  }
}

/* Arrays A and B may have different types, C is the same type as B (B is upcast to the output type) */
void array_number_binary_op(lua_State * L, const Array *A, const void *B,
                     Array *C, ArrayBinaryOperation op, Bool array_first)
{
  const size_t          N = A->size;
  const void  *restrict a = A->data;
  const void  *restrict b = B;
  void        *restrict c = C->data;

  const ArrayType TypeA = A->dtype;
  const ArrayType TypeB = C->dtype;

  switch (op) {
    case ARRAY_OP_ADD:
      ARRAY_OP_BINARY(GET_ORDR, EXPR_BIN_ADD, EXPR_BIN_ADD, EXPR_BIN_ADD, EXPR_BIN_ADD); break;
    case ARRAY_OP_SUB:
      ARRAY_OP_BINARY(GET_ORDR, EXPR_BIN_SUB, EXPR_BIN_SUB, EXPR_BIN_SUB, EXPR_BIN_SUB); break;
    case ARRAY_OP_MUL:
      ARRAY_OP_BINARY(GET_ORDR, EXPR_BIN_MUL, EXPR_BIN_MUL, EXPR_BIN_MUL, EXPR_BIN_MUL); break;
    case ARRAY_OP_DIV:
      ARRAY_OP_BINARY(GET_ORDR, EXPR_BIN_DIV, EXPR_BIN_DIV, EXPR_BIN_DIV, EXPR_BIN_DIV); break;
    case ARRAY_OP_IDIV:
      ARRAY_OP_BINARY(GET_ORDR, EXPR_BIN_DIV, EXPR_BIN_FIDIVF, EXPR_BIN_FIDIV, EXPR_BIN_ERR); break;
    case ARRAY_OP_MOD:
      ARRAY_OP_BINARY(GET_ORDR, EXPR_BIN_MOD, EXPR_BIN_FMODF, EXPR_BIN_FMOD, EXPR_BIN_ERR); break;
    case ARRAY_OP_POW:
      ARRAY_OP_BINARY(GET_ORDR, EXPR_BIN_POW, EXPR_BIN_FPOWF, EXPR_BIN_FPOW, EXPR_BIN_CPOW); break;
    case ARRAY_OP_BAND:
      ARRAY_OP_BINARY(GET_ORDR, EXPR_BIN_BAND, EXPR_BIN_ERR, EXPR_BIN_ERR, EXPR_BIN_ERR); break;
    case ARRAY_OP_BOR:
      ARRAY_OP_BINARY(GET_ORDR, EXPR_BIN_BOR, EXPR_BIN_ERR, EXPR_BIN_ERR, EXPR_BIN_ERR); break;
    case ARRAY_OP_BXOR:
      ARRAY_OP_BINARY(GET_ORDR, EXPR_BIN_BXOR, EXPR_BIN_ERR, EXPR_BIN_ERR, EXPR_BIN_ERR); break;
    case ARRAY_OP_SHL:
      ARRAY_OP_BINARY(GET_ORDR, EXPR_BIN_SHL, EXPR_BIN_ERR, EXPR_BIN_ERR, EXPR_BIN_ERR); break;
    case ARRAY_OP_SHR:
      ARRAY_OP_BINARY(GET_ORDR, EXPR_BIN_SHR, EXPR_BIN_ERR, EXPR_BIN_ERR, EXPR_BIN_ERR); break;
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
  const Ta*restrict _a=a;Tb*restrict _b=b;\
  while (indices[0] < shape[0])\
  {\
    POS_COLMAJOR(cpos, indices, shape, ndims);\
    _b[cpos]=_a[rpos];\
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
