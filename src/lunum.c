
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "lunum.h"
#include "limits.h"



#define LUA_NEW_METAMETHOD(luastate, obj, funcname) {           \
    lua_pushfstring((luastate), "__%s", (#funcname));           \
    lua_pushcfunction((luastate), (luaC_##obj##__##funcname));  \
    lua_settable((luastate), -3);                               \
  }

#define LUA_NEW_MODULEMETHOD(luastate, obj, funcname) {         \
    lua_pushfstring((luastate), "%s", (#funcname));             \
    lua_pushcfunction((luastate), (luaC_##obj##_##funcname));   \
    lua_settable((luastate), -3);                               \
  }

#define LUA_NEW_MODULEDATA(luastate, obj, dataname) {   \
    lua_pushstring((luastate), (#dataname));            \
    lua_pushnumber((luastate), (obj));                  \
    lua_settable((luastate), -3);                       \
  }



static int luaC_lunum_array(lua_State *L);
static int luaC_lunum_zeros(lua_State *L);
static int luaC_lunum_range(lua_State *L);
static int luaC_lunum_resize(lua_State *L);


static int luaC_lunum_sin(lua_State *L);
static int luaC_lunum_cos(lua_State *L);
static int luaC_lunum_tan(lua_State *L);

static int luaC_lunum_sin(lua_State *L);
static int luaC_lunum_cos(lua_State *L);
static int luaC_lunum_tan(lua_State *L);

static int luaC_lunum_asin(lua_State *L);
static int luaC_lunum_acos(lua_State *L);
static int luaC_lunum_atan(lua_State *L);

static int luaC_lunum_sinh(lua_State *L);
static int luaC_lunum_cosh(lua_State *L);
static int luaC_lunum_tanh(lua_State *L);

static int luaC_lunum_asinh(lua_State *L);
static int luaC_lunum_acosh(lua_State *L);
static int luaC_lunum_atanh(lua_State *L);

static int luaC_lunum_exp(lua_State *L);
static int luaC_lunum_log(lua_State *L);
static int luaC_lunum_log10(lua_State *L);

static int luaC_lunum_transpose(lua_State *L);

static int luaC_lunum_conjugate(lua_State *L);
static int luaC_lunum_slice(lua_State *L);
static int luaC_lunum_loadtxt(lua_State *L);
static int luaC_lunum_fromfile(lua_State *L);



static int luaC_array__tostring(lua_State *L);
static int luaC_array__call(lua_State *L);
static int luaC_array__index(lua_State *L);
static int luaC_array__newindex(lua_State *L);
static int luaC_array__add(lua_State *L);
static int luaC_array__sub(lua_State *L);
static int luaC_array__mul(lua_State *L);
static int luaC_array__div(lua_State *L);
static int luaC_array__idiv(lua_State *L);
static int luaC_array__mod(lua_State *L);
static int luaC_array__pow(lua_State *L);
static int luaC_array__unm(lua_State *L);
static int luaC_array__bnot(lua_State *L);
static int luaC_array__band(lua_State *L);
static int luaC_array__bor(lua_State *L);
static int luaC_array__bxor(lua_State *L);
static int luaC_array__shl(lua_State *L);
static int luaC_array__shr(lua_State *L);
static int luaC_array__gc(lua_State *L);


static int luaC_complex__tostring(lua_State *L);
static int luaC_complex__add(lua_State *L);
static int luaC_complex__sub(lua_State *L);
static int luaC_complex__mul(lua_State *L);
static int luaC_complex__div(lua_State *L);
static int luaC_complex__pow(lua_State *L);
static int luaC_complex__unm(lua_State *L);
static int luaC_complex__eq(lua_State *L);
static int luaC_complex__lt(lua_State *L);
static int luaC_complex__le(lua_State *L);


static int _array_unary_op(lua_State *L, ArrayUnaryOperation op);

static int _array_binary_op(lua_State *L, ArrayBinaryOperation op);
static int _array_number_binary_op(lua_State *L, ArrayBinaryOperation op, Bool array_first);
static int _array_array_binary_op(lua_State *L, ArrayBinaryOperation op);

static int _complex_binary_op1(lua_State *L, ArrayBinaryOperation op);
static int _complex_binary_op2(lua_State *L, ArrayBinaryOperation op);

static void _unary_func(lua_State *L, double(*f)(double), Complex(*g)(Complex),
			int cast);
static void _push_value(lua_State *L, ArrayType T, void *v);
static size_t _get_index(lua_State *L, Array *A);


int luaopen_lunum(lua_State *L)
{
  lua_settop(L, 0); // start with an empty stack

  // Create the 'array' metatable
  // ---------------------------------------------------------------------------
  luaL_newmetatable(L, "array");
  LUA_NEW_METAMETHOD(L, array, tostring);
  LUA_NEW_METAMETHOD(L, array, call);
  LUA_NEW_METAMETHOD(L, array, index);
  LUA_NEW_METAMETHOD(L, array, newindex);
  LUA_NEW_METAMETHOD(L, array, add);
  LUA_NEW_METAMETHOD(L, array, sub);
  LUA_NEW_METAMETHOD(L, array, mul);
  LUA_NEW_METAMETHOD(L, array, div);
  LUA_NEW_METAMETHOD(L, array, idiv);
  LUA_NEW_METAMETHOD(L, array, mod);
  LUA_NEW_METAMETHOD(L, array, pow);
  LUA_NEW_METAMETHOD(L, array, band);
  LUA_NEW_METAMETHOD(L, array, bor);
  LUA_NEW_METAMETHOD(L, array, bxor);
  LUA_NEW_METAMETHOD(L, array, shl);
  LUA_NEW_METAMETHOD(L, array, shr);
  LUA_NEW_METAMETHOD(L, array, unm);
  LUA_NEW_METAMETHOD(L, array, bnot);
  LUA_NEW_METAMETHOD(L, array, gc);
  lua_pop(L, 1);

  // Create the 'complex' metatable
  // ---------------------------------------------------------------------------
  luaL_newmetatable(L, "complex");
  LUA_NEW_METAMETHOD(L, complex, tostring);
  LUA_NEW_METAMETHOD(L, complex, add);
  LUA_NEW_METAMETHOD(L, complex, sub);
  LUA_NEW_METAMETHOD(L, complex, mul);
  LUA_NEW_METAMETHOD(L, complex, div);
  LUA_NEW_METAMETHOD(L, complex, pow);
  LUA_NEW_METAMETHOD(L, complex, unm);
  LUA_NEW_METAMETHOD(L, complex, eq);
  LUA_NEW_METAMETHOD(L, complex, lt);
  LUA_NEW_METAMETHOD(L, complex, le);
  lua_pop(L, 1);


  // Create the 'lunum' table
  // ---------------------------------------------------------------------------
  lua_newtable(L);
  LUA_NEW_MODULEMETHOD(L, lunum, array);
  LUA_NEW_MODULEMETHOD(L, lunum, zeros);
  LUA_NEW_MODULEMETHOD(L, lunum, range);
  LUA_NEW_MODULEMETHOD(L, lunum, resize);

  LUA_NEW_MODULEMETHOD(L, lunum, sin);
  LUA_NEW_MODULEMETHOD(L, lunum, cos);
  LUA_NEW_MODULEMETHOD(L, lunum, tan);

  LUA_NEW_MODULEMETHOD(L, lunum, asin);
  LUA_NEW_MODULEMETHOD(L, lunum, acos);
  LUA_NEW_MODULEMETHOD(L, lunum, atan);

  LUA_NEW_MODULEMETHOD(L, lunum, sinh);
  LUA_NEW_MODULEMETHOD(L, lunum, cosh);
  LUA_NEW_MODULEMETHOD(L, lunum, tanh);

  LUA_NEW_MODULEMETHOD(L, lunum, asinh);
  LUA_NEW_MODULEMETHOD(L, lunum, acosh);
  LUA_NEW_MODULEMETHOD(L, lunum, atanh);

  LUA_NEW_MODULEMETHOD(L, lunum, exp);
  LUA_NEW_MODULEMETHOD(L, lunum, log);
  LUA_NEW_MODULEMETHOD(L, lunum, log10);

  LUA_NEW_MODULEMETHOD(L, lunum, transpose);

  LUA_NEW_MODULEMETHOD(L, lunum, conjugate);
  LUA_NEW_MODULEMETHOD(L, lunum, loadtxt);
  LUA_NEW_MODULEMETHOD(L, lunum, fromfile);

  LUA_NEW_MODULEMETHOD(L, lunum, slice);


  LUA_NEW_MODULEDATA(L, ARRAY_TYPE_BOOL   , bool);
  LUA_NEW_MODULEDATA(L, ARRAY_TYPE_CHAR   , char);
  LUA_NEW_MODULEDATA(L, ARRAY_TYPE_SHORT  , short);
  LUA_NEW_MODULEDATA(L, ARRAY_TYPE_INT    , int);
  LUA_NEW_MODULEDATA(L, ARRAY_TYPE_LONG   , long);
  LUA_NEW_MODULEDATA(L, ARRAY_TYPE_SIZE_T , size_t);
  LUA_NEW_MODULEDATA(L, ARRAY_TYPE_FLOAT  , float);
  LUA_NEW_MODULEDATA(L, ARRAY_TYPE_DOUBLE , double);
  LUA_NEW_MODULEDATA(L, ARRAY_TYPE_COMPLEX, complex);

  // Register the purely imaginary number 'I'
  lunum_pushcomplex(L, I);
  lua_setfield(L, 1, "I");

  lua_setglobal(L, "lunum");
#include "array_class.lc" // sets the lunum.__array_methods table

  lua_getglobal(L, "lunum");
  return 1;
}




// *****************************************************************************
// Implementation of lunum.array metatable
//
// *****************************************************************************
static int luaC_array__gc(lua_State *L)
{
  Array *A = lunum_checkarray1(L, 1);
  array_del(A);
  return 0;
}

static int luaC_array__tostring(lua_State *L)
{
  Array *A = lunum_checkarray1(L, 1);

  lua_pushstring(L, "  [ ");
  size_t nstr = 1;
  for (size_t n=0; n<A->size; ++n) {

    char s[64];

    switch (A->dtype) {
    case ARRAY_TYPE_BOOL    : sprintf(s, "%s" , ((Bool*)A->data)[n]?"true":"false"); break;
    case ARRAY_TYPE_CHAR    : sprintf(s, "%d" , ((char   *)A->data)[n]); break;
    case ARRAY_TYPE_SHORT   : sprintf(s, "%d" , ((short  *)A->data)[n]); break;
    case ARRAY_TYPE_INT     : sprintf(s, "%d" , ((int    *)A->data)[n]); break;
    case ARRAY_TYPE_LONG    : sprintf(s, "%ld", ((long   *)A->data)[n]); break;
    case ARRAY_TYPE_SIZE_T  : sprintf(s, "%lu", ((size_t *)A->data)[n]); break;
    case ARRAY_TYPE_FLOAT   : sprintf(s, "%g" , ((float  *)A->data)[n]); break;
    case ARRAY_TYPE_DOUBLE  : sprintf(s, "%g" , ((double *)A->data)[n]); break;
    case ARRAY_TYPE_COMPLEX : sprintf(s, "%g%s%gj",
                                      creal(((Complex*)A->data)[n]),
                                      cimag(((Complex*)A->data)[n]) >= 0.0 ? "+" : "-",
                                      fabs(cimag(((Complex*)A->data)[n]))); break;
    }

    if (n == A->size-1) {
      lua_pushfstring(L, "%s", s);
    }
    else {
      lua_pushfstring(L, "%s, ", s);
    }
    if ((n+1) % 10 == 0 && n != 0 && n != A->size-1) {
      lua_pushstring(L, "\n    "); ++nstr;
    }
  }
  lua_pushstring(L, " ]"); ++nstr;
  lua_concat(L, A->size + nstr);

  return 1;
}

static int luaC_array__call(lua_State *L)
{
  Array *A = lunum_checkarray1(L, 1);
  int nind = lua_gettop(L) - 1;

  if (nind != A->ndims) {
    luaL_error(L, "wrong number of indices (%d) for array of dimension %d",
               nind, A->ndims);
    return 0;
  }
  const int Nd = A->ndims;
  size_t *stride = (size_t*) malloc(A->ndims * sizeof(size_t));
  stride[Nd-1] = 1;

  for (int d=Nd-2; d>=0; --d) {
    stride[d] = stride[d+1] * A->shape[d+1];
  }

  size_t m = 0;
  for (int d=0; d<A->ndims; ++d) {
    lua_Integer i = lua_tointeger(L, d+2);
    m += i*stride[d];
  }

  _push_value(L, A->dtype, (char*)A->data + m*array_sizeof(A->dtype));
  free(stride);

  return 1;
}

static int luaC_array__index(lua_State *L)
{
  Array *A = lunum_checkarray1(L, 1);

  // Figure out what is the format of the input index. If it's a number or a
  // table of numbers, then pass it along to _get_index. If it's a table of
  // tables or numbers, then assume it's a slice. If it's an array of bools,
  // then use it as a mask.
  // ---------------------------------------------------------------------------

  if (lunum_hasmetatable(L, 2, "array")) {
    Array *M = lunum_checkarray1(L, 2);
    if (M->dtype != ARRAY_TYPE_BOOL) {
      luaL_error(L, "index array must be of type bool");
    }
    Array B = array_new_from_mask(A, M);
    lunum_pusharray1(L, &B);
    return 1;
  }
  else if (lua_type(L, 2) == LUA_TTABLE || lua_type(L, 2) == LUA_TSTRING) {

    lua_getglobal(L, "lunum");
    lua_getfield(L, -1, "__build_slice");
    lua_remove(L, -2);
    lua_pushvalue(L, 1);
    lua_pushvalue(L, 2);
    lua_call(L, 2, 1);

    return 1;
  }

  const size_t m = _get_index(L, A);
  _push_value(L, A->dtype, (char*)A->data + array_sizeof(A->dtype)*m);

  return 1;
}


static int luaC_array__newindex(lua_State *L)
{
  Array *A = lunum_checkarray1(L, 1);
  const size_t m = _get_index(L, A);

  const ArrayType T = A->dtype;

  void *val = lunum_tovalue(L, T);
  memcpy((char*)A->data + array_sizeof(T)*m, val, array_sizeof(T));
  free(val);

  return 0;
}

static int luaC_array__add(lua_State *L)  { return _array_binary_op(L, ARRAY_OP_ADD); }
static int luaC_array__sub(lua_State *L)  { return _array_binary_op(L, ARRAY_OP_SUB); }
static int luaC_array__mul(lua_State *L)  { return _array_binary_op(L, ARRAY_OP_MUL); }
static int luaC_array__div(lua_State *L)  { return _array_binary_op(L, ARRAY_OP_DIV); }
static int luaC_array__idiv(lua_State *L) { return _array_binary_op(L, ARRAY_OP_IDIV);}
static int luaC_array__mod(lua_State *L)  { return _array_binary_op(L, ARRAY_OP_MOD); }
static int luaC_array__pow(lua_State *L)  { return _array_binary_op(L, ARRAY_OP_POW); }
static int luaC_array__unm(lua_State *L)  { return _array_unary_op(L,  ARRAY_OP_UNM); }
static int luaC_array__bnot(lua_State *L) { return _array_unary_op(L,  ARRAY_OP_BNOT);}
static int luaC_array__band(lua_State *L) { return _array_binary_op(L, ARRAY_OP_BAND);}
static int luaC_array__bor(lua_State *L)  { return _array_binary_op(L, ARRAY_OP_BOR); }
static int luaC_array__bxor(lua_State *L) { return _array_binary_op(L, ARRAY_OP_BXOR);}
static int luaC_array__shl(lua_State *L)  { return _array_binary_op(L, ARRAY_OP_SHL); }
static int luaC_array__shr(lua_State *L)  { return _array_binary_op(L, ARRAY_OP_SHR); }


static int _array_unary_op(lua_State *L, ArrayUnaryOperation op)
{
  Array *A = lunum_checkarray1(L, 1);

  const size_t N = A->size;
  ArrayType T = A->dtype;

  Array B = array_new_zeros(N, T);
  array_resize(&B, A->shape, A->ndims);
  lunum_pusharray1(L, &B);

  array_unary_op(L, A, &B, op);

  return 1;
}

static int _array_binary_op(lua_State *L, ArrayBinaryOperation op)
{
  if (lua_istable(L, 1) && lua_istable(L, 2)) {
    /* both args are tables, upcast to arrays if not already */
    if (!lunum_hasmetatable(L, 1, "array")) {
      Array *B = lunum_checkarray1(L, 2);
      lunum_upcast(L, 1, B->dtype, B->size);
      lua_replace(L, 1);
      Array *A = lunum_checkarray1(L, 1);
      array_resize(A, B->shape, B->ndims);
    }
    if (!lunum_hasmetatable(L, 2, "array")) {
      Array *A = lunum_checkarray1(L, 1);
      lunum_upcast(L, 2, A->dtype, A->size);
      lua_replace(L, 2);
      Array *B = lunum_checkarray1(L, 2);
      array_resize(B, A->shape, A->ndims);
    }
    return _array_array_binary_op(L, op);
  } else {
    /* one arg is not a table(array) */
    return _array_number_binary_op(L, op, lunum_hasmetatable(L, 1, "array"));
  }
}

static int _array_number_binary_op(lua_State *L, ArrayBinaryOperation op, Bool array_first)
{
  const Array *A = array_first ? lunum_checkarray1(L, 1) : lunum_checkarray1(L, 2);
  Array A_;
  ArrayType T = A->dtype;

  int num_pos = array_first ? 2 : 1;
  union {
    Bool b;
    char c;
    short s;
    int i;
    long l;
    size_t t;
    float f;
    double d;
    Complex z;
    lua_Integer li;
    lua_Number ln;
  } num;

  if (lua_isboolean(L, num_pos)) {
    num.i = lua_toboolean(L, num_pos);
    /* number can't have a higher type, upgrade to type T */
    switch (T) {
      case ARRAY_TYPE_BOOL    : num.b = (Bool)num.i;    break;
      case ARRAY_TYPE_CHAR    : num.c = (char)num.i;    break;
      case ARRAY_TYPE_SHORT   : num.s = (short)num.i;   break;
      case ARRAY_TYPE_INT     : num.i = (int)num.i;     break;
      case ARRAY_TYPE_LONG    : num.l = (long)num.i;    break;
      case ARRAY_TYPE_SIZE_T  : num.t = (size_t)num.i;  break;
      case ARRAY_TYPE_FLOAT   : num.f = (float)num.i;   break;
      case ARRAY_TYPE_DOUBLE  : num.d = (double)num.i;  break;
      case ARRAY_TYPE_COMPLEX : num.z = (Complex)num.i; break;
    }
    A_ = *A;
  }
  else if (lua_isinteger(L, num_pos)) {
    num.li = lua_tointeger(L, num_pos);
    if (T >= ARRAY_TYPE_LONG) {
      /* A has higher type, upgrade to type T */
      switch (T) {
        case ARRAY_TYPE_BOOL    : num.b = (Bool)num.li;    break;
        case ARRAY_TYPE_CHAR    : num.c = (char)num.li;    break;
        case ARRAY_TYPE_SHORT   : num.s = (short)num.li;   break;
        case ARRAY_TYPE_INT     : num.i = (int)num.li;     break;
        case ARRAY_TYPE_LONG    : num.l = (long)num.li;    break;
        case ARRAY_TYPE_SIZE_T  : num.t = (size_t)num.li;  break;
        case ARRAY_TYPE_FLOAT   : num.f = (float)num.li;   break;
        case ARRAY_TYPE_DOUBLE  : num.d = (double)num.li;  break;
        case ARRAY_TYPE_COMPLEX : num.z = (Complex)num.li; break;
      }
      A_ = *A;
    } else {
      /* number has higher type, copy A to higher type */
      T = ARRAY_TYPE_LONG;
      A_ = array_new_copy(A, T);
    }
  }
  else if (lua_isnumber(L, num_pos)) {
    num.ln = lua_tonumber(L, num_pos);
    if (T >= ARRAY_TYPE_DOUBLE) {
      /* A has higher type, upgrade to type T */
      switch (T) {
        case ARRAY_TYPE_BOOL    : num.b = (Bool)num.ln;    break;
        case ARRAY_TYPE_CHAR    : num.c = (char)num.ln;    break;
        case ARRAY_TYPE_SHORT   : num.s = (short)num.ln;   break;
        case ARRAY_TYPE_INT     : num.i = (int)num.ln;     break;
        case ARRAY_TYPE_LONG    : num.l = (long)num.ln;    break;
        case ARRAY_TYPE_SIZE_T  : num.t = (size_t)num.ln;  break;
        case ARRAY_TYPE_FLOAT   : num.f = (float)num.ln;   break;
        case ARRAY_TYPE_DOUBLE  : num.d = (double)num.ln;  break;
        case ARRAY_TYPE_COMPLEX : num.z = (Complex)num.ln; break;
      }
      A_ = *A;
    } else {
      /* number has higher type, copy A to higher type */
      T = ARRAY_TYPE_DOUBLE;
      A_ = array_new_copy(A, T);
    }
  }
  else if (lunum_hasmetatable(L, num_pos, "complex")) {
    /* number complex, copy A to complex if not already complex */
    num.z = *((Complex*) lua_touserdata(L, num_pos));
    A_ = (T == ARRAY_TYPE_COMPLEX) ? *A : array_new_copy(A, T = ARRAY_TYPE_COMPLEX);
  } else {
    return luaL_error(L, "Invalid argument in Array binary op");
  }

  Array C = array_new_zeros(A->size, T);
  array_resize(&C, A->shape, A->ndims);
  lunum_pusharray1(L, &C);

  array_number_binary_op(L, &A_, (void *)&num, &C, op, array_first);

  if (A->dtype != T) array_del(&A_);

  return 1;
}

static int _array_array_binary_op(lua_State *L, ArrayBinaryOperation op)
{
  Array *A = lunum_checkarray1(L, 1);
  Array *B = lunum_checkarray1(L, 2);


  if (A->ndims != B->ndims) {
    luaL_error(L, "arrays have different dimensions");
  }
  for (int d=0; d<A->ndims; ++d) {
    if (A->shape[d] != B->shape[d]) {
      luaL_error(L, "arrays shapes do not agree");
    }
  }

  ArrayType T = (A->dtype >= B->dtype) ? A->dtype : B->dtype;

  Array A_ = (A->dtype == T) ? *A : array_new_copy(A, T);
  Array B_ = (B->dtype == T) ? *B : array_new_copy(B, T);

  Array C = array_new_zeros(A->size, T);
  array_resize(&C, A->shape, A->ndims);
  lunum_pusharray1(L, &C);

  array_array_binary_op(L, &A_, &B_, &C, op);

  if (A->dtype != T) array_del(&A_);
  if (B->dtype != T) array_del(&B_);

  return 1;
}



// *****************************************************************************
// Implementation of lunum.complex metatable
//
// *****************************************************************************
static int luaC_complex__tostring(lua_State *L)
{
  Complex z = *((Complex*) luaL_checkudata(L, 1, "complex"));

  lua_pushfstring(L, "%f%s%fj", creal(z), cimag(z)>=0.0?"+":"-", fabs(cimag(z)));
  return 1;
}



static int luaC_complex__add(lua_State *L) { return _complex_binary_op1(L, ARRAY_OP_ADD); }
static int luaC_complex__sub(lua_State *L) { return _complex_binary_op1(L, ARRAY_OP_SUB); }
static int luaC_complex__mul(lua_State *L) { return _complex_binary_op1(L, ARRAY_OP_MUL); }
static int luaC_complex__div(lua_State *L) { return _complex_binary_op1(L, ARRAY_OP_DIV); }
static int luaC_complex__pow(lua_State *L) { return _complex_binary_op1(L, ARRAY_OP_POW); }
static int luaC_complex__unm(lua_State *L) { 
  Complex v = *((Complex*) luaL_checkudata(L, 1, "complex"));
  Complex *z = (Complex*) lua_newuserdata(L, sizeof(Complex));
  luaL_setmetatable(L, "complex");
  *z = -v;
  return(1);
}


// -----------------------------------------------------------------------------
// Use a dictionary ordering on the complex numbers. Might not be useful too
// often, but it's better than having this behavior undefined.
// -----------------------------------------------------------------------------
#define LUA_COMPARISON(comp)				\
  {							\
    Complex z1 = lunum_checkcomplex(L, 1);		\
    Complex z2 = lunum_checkcomplex(L, 2);		\
							\
    if (creal(z1) != creal(z2)) {			\
      lua_pushboolean(L, creal(z1) comp creal(z2));	\
    }							\
    else {						\
      lua_pushboolean(L, cimag(z1) comp cimag(z2));	\
    }							\
    return 1;						\
  }							\


static int luaC_complex__lt(lua_State *L) LUA_COMPARISON(<);
static int luaC_complex__le(lua_State *L) LUA_COMPARISON(<=);
static int luaC_complex__eq(lua_State *L)
{
  Complex z1 = lunum_checkcomplex(L, 1);
  Complex z2 = lunum_checkcomplex(L, 2);
  lua_pushboolean(L, z1==z2);
  return 1;
}


static int _complex_binary_op1(lua_State *L, ArrayBinaryOperation op)
{
  if (lunum_hasmetatable(L, 1, "array") ||
      lunum_hasmetatable(L, 2, "array")) {
    return _array_binary_op(L, op);
  }
  if (!lunum_hasmetatable(L, 1, "complex")) {
    lunum_pushcomplex(L, lua_tonumber(L, 1));
    lua_replace(L, 1);
  }
  if (!lunum_hasmetatable(L, 2, "complex")) {
    lunum_pushcomplex(L, lua_tonumber(L, 2));
    lua_replace(L, 2);
  }
  return _complex_binary_op2(L, op);
}

static int _complex_binary_op2(lua_State *L, ArrayBinaryOperation op)
{
  Complex v = *((Complex*) luaL_checkudata(L, 1, "complex"));
  Complex w = *((Complex*) luaL_checkudata(L, 2, "complex"));

  Complex *z = (Complex*) lua_newuserdata(L, sizeof(Complex));
  luaL_setmetatable(L, "complex");

  switch (op) {
  case ARRAY_OP_ADD: *z = v + w; break;
  case ARRAY_OP_SUB: *z = v - w; break;
  case ARRAY_OP_MUL: *z = v * w; break;
  case ARRAY_OP_DIV: *z = v / w; break;
  case ARRAY_OP_POW: *z = cpow(v,w); break;
  case ARRAY_OP_IDIV:
  case ARRAY_OP_MOD:
  case ARRAY_OP_BAND:
  case ARRAY_OP_BOR:
  case ARRAY_OP_BXOR:
  case ARRAY_OP_SHL:
  case ARRAY_OP_SHR: luaL_error(L, "Invalid operation"); break;
  }

  return 1;
}








static int luaC_lunum_array(lua_State *L)
{
  if (lua_type(L, 2) == LUA_TSTRING) {
    const ArrayType T = array_typeflag(lua_tostring(L, 2)[0]);
    lunum_upcast(L, 1, T, 1);
  }
  else {
    const ArrayType T = (ArrayType) luaL_optinteger(L, 2, ARRAY_TYPE_DOUBLE);
    lunum_upcast(L, 1, T, 1);
  }
  return 1;
}

static int luaC_lunum_zeros(lua_State *L)
{
  if (lua_isnumber(L, 1)) {
    const lua_Integer N = luaL_checkinteger(L, 1);
    const ArrayType T = (ArrayType) luaL_optinteger(L, 2, ARRAY_TYPE_DOUBLE);
    Array A = array_new_zeros(N, T);
    lunum_pusharray1(L, &A);
    return 1;
  }
  else if (lua_istable(L, 1) || lunum_hasmetatable(L, 1, "array")) {

    size_t Nd_t;
    size_t *N = (size_t*) lunum_checkarray2(L, 1, ARRAY_TYPE_SIZE_T, &Nd_t);
    int Nd = (int)Nd_t;
    const ArrayType T = (ArrayType) luaL_optinteger(L, 2, ARRAY_TYPE_DOUBLE);

    size_t ntot = 1;
    for (int d=0; d<Nd; ++d) ntot *= N[d];
    Array A = array_new_zeros(ntot, T);

    array_resize(&A, N, Nd);
    lunum_pusharray1(L, &A);

    return 1;
  }
  else {
    luaL_error(L, "argument must be either number, table, or array");
    return 0;
  }
}

static int luaC_lunum_range(lua_State *L)
{
  const lua_Integer N = luaL_checkinteger(L, 1);
  if (N <= INT_MAX) {
    Array A = array_new_zeros(N, ARRAY_TYPE_INT);
    lunum_pusharray1(L, &A);
    for (size_t i=0; i<N; ++i) {
      ((int*)A.data)[i] = i;
    }
  } else {
    Array A = array_new_zeros(N, ARRAY_TYPE_LONG);
    lunum_pusharray1(L, &A);
    for (size_t i=0; i<N; ++i) {
      ((long*)A.data)[i] = i;
    }
  }
  
  return 1;
}

static int luaC_lunum_resize(lua_State *L)
{
  size_t Nd_t;
  Array *A = lunum_checkarray1(L, 1); // the array to resize
  size_t *N = (size_t*) lunum_checkarray2(L, 2, ARRAY_TYPE_SIZE_T, &Nd_t);
  int Nd = (int)Nd_t;

  size_t ntot = 1;
  for (int d=0; d<Nd; ++d) ntot *= N[d];

  if (A->size != ntot) {
    luaL_error(L, "new and old total sizes do not agree");
    return 0;
  }
  array_resize(A, N, Nd);

  return 0;
}

static int luaC_lunum_transpose(lua_State *L)
{
  const Array *A = lunum_checkarray1(L, 1); // the array to transpose

  Array B = array_new_zeros(A->size, A->dtype);
  /* copy transposed shape and size to B */
  B.ndims = A->ndims;
  B.shape = (size_t *)malloc(B.ndims * sizeof(size_t));
  for (int i = 0; i < A->ndims; i++)
    B.shape[i] = A->shape[A->ndims - i - 1];

  lunum_pusharray1(L, &B);

  array_transpose(A, &B);

  return 1;

}

static int luaC_lunum_slice(lua_State *L)
{

  // The first part of this function extracts a slice of the array 'A' according
  // to the convention start:stop:skip. The result is a contiguous array 'B'
  // having the same number of dimensions as 'A'.
  // ---------------------------------------------------------------------------
  size_t Nd0_t, Nd1_t, Nd2_t, Nd3_t;

  const Array *A = lunum_checkarray1(L, 1); // the array to resize
  size_t *start   = (size_t*) lunum_checkarray2(L, 2, ARRAY_TYPE_SIZE_T, &Nd0_t);
  size_t *stop    = (size_t*) lunum_checkarray2(L, 3, ARRAY_TYPE_SIZE_T, &Nd1_t);
  size_t *skip    = (size_t*) lunum_checkarray2(L, 4, ARRAY_TYPE_SIZE_T, &Nd2_t);
  size_t *squeeze = (size_t*) lunum_checkarray2(L, 5, ARRAY_TYPE_SIZE_T, &Nd3_t);

  int Nd0 = (int)Nd0_t, Nd1 = (int)Nd1_t, Nd2 = (int)Nd2_t, Nd3 = (int)Nd3_t;


  if (Nd0 != A->ndims || Nd1 != A->ndims || Nd2 != A->ndims || Nd3 != A->ndims) {
    luaL_error(L, "slice has wrong number of dimensions for array");
  }

  for (int d=0; d<A->ndims; ++d) {
    if (start[d] < 0 || stop[d] > A->shape[d]) {
      luaL_error(L, "slice not within array extent");
    }
  }
  Array B = array_new_from_slice(A, start, stop, skip, Nd0);


  // The rest of this function deals with squeezing out the size-1 dimensions of
  // 'B' which are marked by the 'squeeze' array.
  // ---------------------------------------------------------------------------
  size_t Nd_new = 0;
  for (int d=0; d<Nd0; ++d) Nd_new += !squeeze[d];

  // In case we're left with a 0-dimensional (scalar) slice
  if (Nd_new == 0) {
    _push_value(L, B.dtype, B.data);
    return 1;
  }
  // In case there are any dims to squeeze out
  else if (Nd_new != Nd0) {

    size_t *shape_new = (size_t*) malloc(Nd_new * sizeof(size_t));
    for (int d=0,e=0; d<Nd0; ++d) {
      if (B.shape[d] > 1 || !squeeze[d]) {
	shape_new[e] = B.shape[d];
	++e;
      }
    }
    array_resize(&B, shape_new, Nd_new);
    free(shape_new);
  }

  lunum_pusharray1(L, &B);
  return 1;
}


// Functions used by unary predicates
// -----------------------------------------------------------------------------
static double rconj(double x) { return x; }   // conj, but for real argument
// -----------------------------------------------------------------------------

static int luaC_lunum_sin(lua_State *L) { _unary_func(L, sin, csin, 1); return 1; }
static int luaC_lunum_cos(lua_State *L) { _unary_func(L, cos, ccos, 1); return 1; }
static int luaC_lunum_tan(lua_State *L) { _unary_func(L, tan, ctan, 1); return 1; }

static int luaC_lunum_asin(lua_State *L) { _unary_func(L, asin, casin, 1); return 1; }
static int luaC_lunum_acos(lua_State *L) { _unary_func(L, acos, cacos, 1); return 1; }
static int luaC_lunum_atan(lua_State *L) { _unary_func(L, atan, catan, 1); return 1; }

static int luaC_lunum_sinh(lua_State *L) { _unary_func(L, sinh, csinh, 1); return 1; }
static int luaC_lunum_cosh(lua_State *L) { _unary_func(L, cosh, ccosh, 1); return 1; }
static int luaC_lunum_tanh(lua_State *L) { _unary_func(L, tanh, ctanh, 1); return 1; }

static int luaC_lunum_asinh(lua_State *L) { _unary_func(L, asinh, casinh, 1); return 1; }
static int luaC_lunum_acosh(lua_State *L) { _unary_func(L, acosh, cacosh, 1); return 1; }
static int luaC_lunum_atanh(lua_State *L) { _unary_func(L, atanh, catanh, 1); return 1; }

static int luaC_lunum_exp(lua_State *L) { _unary_func(L, exp, cexp, 1); return 1; }
static int luaC_lunum_log(lua_State *L) { _unary_func(L, log, clog, 1); return 1; }
static int luaC_lunum_log10(lua_State *L) { _unary_func(L, log10, NULL, 1); return 1; }
static int luaC_lunum_conjugate(lua_State *L) { _unary_func(L, rconj, conj, 0); return 1; }


static int luaC_lunum_loadtxt(lua_State *L)
// -----------------------------------------------------------------------------
// Opens the text file 'fname' for reading, and parses the data
// line-by-line. It is assumed that the data is all floating point, and that
// only a space is used as a separator. If there are multiple columns then a 2d
// array is created. All rows must have the same number of entries, otherwise an
// error is generated.
// -----------------------------------------------------------------------------
{
  const char *fname = luaL_checkstring(L, 1);
  FILE *input = fopen(fname, "r");

  if (input == NULL) {
    luaL_error(L, "no such file %s", fname);
  }

  size_t nline = 0;
  size_t ncols = 0;
  size_t ntot = 0;
  double *data = NULL;

  char line[2048];

  while (fgets(line, sizeof(line), input)) {

    if (strlen(line) == 1) {
      continue;
    }

    size_t nvals = 0;
    double *vals = NULL;
    char *word = strtok(line, " \n");

    while (word) {
      vals = (double*) realloc(vals, ++nvals*sizeof(double));
      vals[nvals-1] = atof(word);
      word = strtok(NULL, " \n");
    }

    if (ncols == 0) ncols = nvals;
    if (ncols != nvals) {
      luaL_error(L, "wrong number of data on line %ld of %s", nline, fname);
    }

    data = (double*) realloc(data, (ntot+=nvals)*sizeof(double));
    memcpy(data+ntot-nvals, vals, nvals*sizeof(double));
    free(vals);

    ++nline;
  }
  fclose(input);

  lunum_pusharray2(L, data, ARRAY_TYPE_DOUBLE, ntot);
  Array *A = lunum_checkarray1(L, -1);

  size_t shape[2] = { nline, ncols };
  array_resize(A, shape, ncols == 1 ? 1 : 2);

  free(data);
  return 1;
}


static int luaC_lunum_fromfile(lua_State *L)
// -----------------------------------------------------------------------------
// Opens the binary file 'fname' for reading, and returns a 1d array from the
// data. The file size must be a multiple of the data type 'T'.
// -----------------------------------------------------------------------------
{
  const char *fname = luaL_checkstring(L, 1);
  const ArrayType T = luaL_optinteger(L, 2, ARRAY_TYPE_DOUBLE);
  const int sizeof_T = array_sizeof(T);

  FILE *input = fopen(fname, "rb");

  if (input == NULL) {
    luaL_error(L, "no such file %s", fname);
  }
  fseek(input, 0L, SEEK_END); const size_t sz = ftell(input);
  fseek(input, 0L, SEEK_SET);

  if (sz % sizeof_T != 0) {
    luaL_error(L, "file size must be a multiple of the data type size");
  }
  const size_t N = sz / sizeof_T;
  Array A = array_new_zeros(N, T);

  if (fread(A.data, N, sizeof_T, input) != sizeof_T) {
    fclose(input);
    return luaL_error(L, "Error while reading file %s", fname);
  }
  fclose(input);
  lunum_pusharray1(L, &A);

  return 1;
}


#define EXPR_EVALF(T,N,x) {for(size_t i=0;i<N;++i)((T*)(x))[i]=f(((T*)(x))[i]);}
#define EXPR_EVALG(T,N,x) {for(size_t i=0;i<N;++i)((T*)(x))[i]=g(((T*)(x))[i]);}

void _unary_func(lua_State *L, double(*f)(double), Complex(*g)(Complex), int cast)
{
  if (lua_isnumber(L, 1)) {
    const double x = lua_tonumber(L, 1);
    lua_pushnumber(L, f(x));
  }
  else if (lunum_hasmetatable(L, 1, "complex")) {

    if (g == NULL) {
      luaL_error(L, "complex operation not supported");
    }

    const Complex z = lunum_checkcomplex(L, 1);
    lunum_pushcomplex(L, g(z));
  }
  else if (lunum_hasmetatable(L, 1, "array")) {
    Array *A = (Array*) lunum_checkarray1(L, 1);

    if (cast == 0) {
      Array B = array_new_copy(A, A->dtype);

      switch (B.dtype) {
      case ARRAY_TYPE_BOOL    : EXPR_EVALF(Bool   , B.size, B.data); break;
      case ARRAY_TYPE_CHAR    : EXPR_EVALF(char   , B.size, B.data); break;
      case ARRAY_TYPE_SHORT   : EXPR_EVALF(short  , B.size, B.data); break;
      case ARRAY_TYPE_INT     : EXPR_EVALF(int    , B.size, B.data); break;
      case ARRAY_TYPE_LONG    : EXPR_EVALF(long   , B.size, B.data); break;
      case ARRAY_TYPE_SIZE_T  : EXPR_EVALF(size_t , B.size, B.data); break;
      case ARRAY_TYPE_FLOAT   : EXPR_EVALF(float  , B.size, B.data); break;
      case ARRAY_TYPE_DOUBLE  : EXPR_EVALF(double , B.size, B.data); break;
      case ARRAY_TYPE_COMPLEX : EXPR_EVALG(Complex, B.size, B.data); break;
      }

      lunum_pusharray1(L, &B);
    }
    else if (A->dtype <= ARRAY_TYPE_DOUBLE) {
      Array B = array_new_copy(A, ARRAY_TYPE_DOUBLE);
      double *b = (double*) B.data;
      for (size_t i=0; i<B.size; ++i) b[i] = f(b[i]);
      lunum_pusharray1(L, &B);
    }
    else if (A->dtype == ARRAY_TYPE_COMPLEX) {

      if (g == NULL) {
        luaL_error(L, "complex operation not supported");
      }

      Array B = array_new_copy(A, ARRAY_TYPE_COMPLEX);
      Complex *b = (Complex*) B.data;
      for (size_t i=0; i<B.size; ++i) b[i] = g(b[i]);
      lunum_pusharray1(L, &B);
    }
  }
}

size_t _get_index(lua_State *L, Array *A)
{
  lua_Integer m = 0;

  if (lua_isnumber(L, 2)) {
    m = luaL_checkinteger(L, 2);

    if (m >= A->size || m < 0) {
      luaL_error(L, "index %ld out of bounds on array of length %lu", m, A->size);
    }
  }
  else if (lua_istable(L, 2)) {
    size_t Nd_t;
    size_t *ind = (size_t*) lunum_checkarray2(L, 2, ARRAY_TYPE_SIZE_T, &Nd_t);
    int Nd = (int)Nd_t;

    if (A->ndims != Nd) {
      luaL_error(L, "wrong number of indices (%lu) on array of dimension %d",
                 Nd, A->ndims);
    }
    size_t *stride = (size_t*) malloc(A->ndims * sizeof(size_t));
    stride[Nd-1] = 1;

    for (int d=Nd-2; d>=0; --d) {
      stride[d] = stride[d+1] * A->shape[d+1];
    }

    for (int d=0; d<A->ndims; ++d) {
      if (ind[d] >= A->shape[d] || ind[d] < 0) {
        luaL_error(L, "array indexed out of bounds (%lu) on dimension %lu of size %lu",
                   ind[d], d, A->shape[d]);
      }
      m += ind[d]*stride[d];
    }
    free(stride);
  }
  return m;
}


void _push_value(lua_State *L, ArrayType T, void *v)
{
  switch (T) {
  case ARRAY_TYPE_BOOL    : lua_pushboolean(L,    *((Bool   *)v)); break;
  case ARRAY_TYPE_CHAR    : lua_pushnumber (L,    *((char   *)v)); break;
  case ARRAY_TYPE_SHORT   : lua_pushnumber (L,    *((short  *)v)); break;
  case ARRAY_TYPE_INT     : lua_pushnumber (L,    *((int    *)v)); break;
  case ARRAY_TYPE_LONG    : lua_pushnumber (L,    *((long   *)v)); break;
  case ARRAY_TYPE_SIZE_T  : lua_pushnumber (L,    *((size_t *)v)); break;
  case ARRAY_TYPE_FLOAT   : lua_pushnumber (L,    *((float  *)v)); break;
  case ARRAY_TYPE_DOUBLE  : lua_pushnumber (L,    *((double *)v)); break;
  case ARRAY_TYPE_COMPLEX : lunum_pushcomplex (L, *((Complex*)v)); break;
  }
}

