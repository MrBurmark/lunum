

#ifndef __Array_HEADER__
#define __Array_HEADER__

#include "lua.h"
#include "lualib.h"

typedef      unsigned char Bool;

typedef enum {
  ARRAY_TYPE_BOOL,
  ARRAY_TYPE_CHAR,
  ARRAY_TYPE_SHORT,
  ARRAY_TYPE_INT,
  ARRAY_TYPE_LONG,
  ARRAY_TYPE_SIZE_T,
  ARRAY_TYPE_FLOAT,
  ARRAY_TYPE_DOUBLE,
  ARRAY_TYPE_COMPLEX,
} ArrayType;

typedef enum {
  ARRAY_OP_ADD,
  ARRAY_OP_SUB,
  ARRAY_OP_MUL,
  ARRAY_OP_DIV,
  ARRAY_OP_IDIV,
  ARRAY_OP_MOD,
  ARRAY_OP_POW,
  ARRAY_OP_BAND,
  ARRAY_OP_BOR,
  ARRAY_OP_BXOR,
  ARRAY_OP_SHL,
  ARRAY_OP_SHR,
} ArrayBinaryOperation;

typedef enum {
  ARRAY_OP_UNM,
  ARRAY_OP_BNOT,
} ArrayUnaryOperation;

typedef struct {
  void *data;
  ArrayType dtype;
  Bool owns;
  int ndims;
  size_t size, *shape;
} Array;

char      *array_typename(ArrayType T);
double     array_typemin(ArrayType T);
double     array_typemax(ArrayType T);
ArrayType  array_typeflag(char c);
Array      array_new_zeros(size_t N, ArrayType T);
Array      array_new_copy(const Array *B, ArrayType T);
Array      array_new_from_slice(const Array *B1,
				   size_t *start, size_t *stop, size_t *skip, int Nd);
Array      array_new_from_mask(const Array *B1, Array *M);
void       array_del(Array *A);
void       array_assign_from_scalar(Array *A, const void *val);
void       array_assign_from_array(Array *A, const Array *B);
void       array_array_binary_op(lua_State * L, const Array *A,
			                const Array *B, Array *C, ArrayBinaryOperation op);
void       array_number_binary_op(lua_State * L, const Array *A, const void *B,
                     Array *C, ArrayBinaryOperation op, Bool array_first);
void       array_unary_op(lua_State * L, const Array *A,
                      Array *B, ArrayUnaryOperation op);
/* returns the sizeof( the array data type ) */
int        array_sizeof(ArrayType T);
int        array_resize(Array *A, const size_t *N, int Nd);
void       array_transpose(const Array *A, Array *B);
void       array_extract_slice(Array *B0, const Array *B1,
				   size_t *start, size_t *size, size_t *stride, int Nd);

#endif // __Array_HEADER__
