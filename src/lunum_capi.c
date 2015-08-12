

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "lunum.h"


void lunum_pusharray1(lua_State *L, Array *B)
{
  Array *A = (Array*) lua_newuserdata(L, sizeof(Array));
  *A = *B;

  luaL_setmetatable(L, "array");
}

void lunum_pusharray2(lua_State *L, void *data, ArrayType T, size_t N)
{
  Array A = array_new_zeros(N, T);
  memcpy(A.data, data, N*array_sizeof(T));
  lunum_pusharray1(L, &A);
}

Array *lunum_checkarray1(lua_State *L, int pos)
{
  return (Array*) luaL_checkudata(L, pos, "array");
}

void *lunum_checkarray2(lua_State *L, int pos, ArrayType T, size_t *N)
{
  if (lunum_upcast(L, pos, T, 1)) {
    lua_replace(L, pos);
  }
  Array *A = lunum_checkarray1(L, pos);
  if (N != NULL) *N = A->size;
  return A->data;
}


void lunum_astable(lua_State *L, int pos)
{
  Array *A = lunum_checkarray1(L, pos);
  const void *a = A->data;

  lua_createtable(L, A->size, 0);
  for (size_t i=0; i < A->size; ++i) {

    switch (A->dtype) {
      case ARRAY_TYPE_BOOL    : lua_pushboolean  (L, ((Bool    *)a)[i]); break;
      case ARRAY_TYPE_CHAR    : lua_pushinteger  (L, ((char    *)a)[i]); break;
      case ARRAY_TYPE_SHORT   : lua_pushinteger  (L, ((short   *)a)[i]); break;
      case ARRAY_TYPE_INT     : lua_pushinteger  (L, ((int     *)a)[i]); break;
      case ARRAY_TYPE_LONG    : lua_pushinteger  (L, ((long    *)a)[i]); break;
      case ARRAY_TYPE_SIZE_T  : lua_pushinteger  (L, ((size_t  *)a)[i]); break;
      case ARRAY_TYPE_FLOAT   : lua_pushnumber   (L, ((float   *)a)[i]); break;
      case ARRAY_TYPE_DOUBLE  : lua_pushnumber   (L, ((double  *)a)[i]); break;
      case ARRAY_TYPE_COMPLEX : lunum_pushcomplex(L, ((Complex *)a)[i]); break;
    }

    lua_seti(L, -2, i+1);
  }
}

void lunum_pushcomplex(lua_State *L, Complex z)
{
  Complex *w = (Complex*) lua_newuserdata(L, sizeof(Complex));
  luaL_setmetatable(L, "complex");
  *w = z;
}

Complex lunum_checkcomplex(lua_State *L, int n)
{
  Complex *w = (Complex*) luaL_checkudata(L, n, "complex");
  return *w;
}



int lunum_upcast(lua_State *L, int pos, ArrayType T, size_t N)
// -----------------------------------------------------------------------------
// If the object at position 'pos' is already an array of dtype 'T', then push
// nothing and return 0. If the dtype is not 'T', then return 1 and push a copy
// of that array with dtype 'T' onto the stack. If it is a table, then push an
// array of dtype 'T' having the length of the table. If it is a number or
// complex, then push an array of dtype double or complex respectively having
// length 'N'.
// -----------------------------------------------------------------------------
{
  if (array_typename(T) == NULL) {
    return luaL_error(L, "invalid array type");
  }

  // Deal with lunum.array
  // ---------------------------------------------------------------------------
  if (lunum_hasmetatable(L, pos, "array")) {

    Array *A = lunum_checkarray1(L, pos);

    if (A->dtype == T) {
      return 0;
    }

    else {

      Array A_ = array_new_copy(A, T);
      lunum_pusharray1(L, &A_);
      return 1;
    }
  }

  // Deal with Lua table
  // ---------------------------------------------------------------------------
  else if (lua_istable(L, pos)) {

    Array A = array_new_zeros(lua_rawlen(L, pos), T);

    for (size_t i=0; i<A.size; ++i) {

      lua_geti(L, pos, i+1);

      ArrayAllNum val;
      lunum_tovalue(L, T, &val);
      memcpy((char*)A.data + array_sizeof(T)*i, &val, array_sizeof(T));

      lua_pop(L, 1);
    }
    lunum_pusharray1(L, &A);

    return 1;
  }

  // Deal with Lua bool
  // ---------------------------------------------------------------------------
  else if (lua_isboolean(L, pos)) {
    const Bool x = lua_toboolean(L, pos);
    Array A = array_new_zeros(N, ARRAY_TYPE_BOOL);
    array_assign_from_scalar(&A, &x);
    lunum_pusharray1(L, &A);
    return 1;
  }

  // Deal with Lua numbers
  // ---------------------------------------------------------------------------
  else if (lua_isnumber(L, pos)) {
    const double x = lua_tonumber(L, pos);
    Array A = array_new_zeros(N, ARRAY_TYPE_DOUBLE);
    array_assign_from_scalar(&A, &x);
    Array B = array_new_copy(&A, T);
    array_del(&A);
    lunum_pusharray1(L, &B);
    return 1;
  }

  // Deal with lunum.complex
  // ---------------------------------------------------------------------------
  else if (lunum_hasmetatable(L, pos, "complex")) {

    const Complex z = *((Complex*) lua_touserdata(L, pos));
    Array A = array_new_zeros(N, ARRAY_TYPE_COMPLEX);
    array_assign_from_scalar(&A, &z);
    lunum_pusharray1(L, &A);
    return 1;
  }

  // Throw an error
  // ---------------------------------------------------------------------------
  else {
    return luaL_error(L, "cannot cast to array from object of dtype %s\n",
               lua_typename(L, lua_type(L, pos)));
  }
}


int lunum_hasmetatable(lua_State *L, int pos, const char *name)
{
  const int top = lua_gettop(L);
  if (lua_getmetatable(L, pos) == 0) return 0;
  luaL_getmetatable(L, name);
  const int eq = lua_rawequal(L, -2, -1);
  lua_settop(L, top);
  return eq;
}

#define ASSIGN_TO_VOID(T, n, x) \
    switch (T) {\
      case ARRAY_TYPE_BOOL    : *((Bool    *)n) = (Bool)   (x); break;\
      case ARRAY_TYPE_CHAR    : *((char    *)n) = (char)   (x); break;\
      case ARRAY_TYPE_SHORT   : *((short   *)n) = (short)  (x); break;\
      case ARRAY_TYPE_INT     : *((int     *)n) = (int)    (x); break;\
      case ARRAY_TYPE_LONG    : *((long    *)n) = (long)   (x); break;\
      case ARRAY_TYPE_SIZE_T  : *((size_t  *)n) = (size_t) (x); break;\
      case ARRAY_TYPE_FLOAT   : *((float   *)n) = (float)  (x); break;\
      case ARRAY_TYPE_DOUBLE  : *((double  *)n) = (double) (x); break;\
      case ARRAY_TYPE_COMPLEX : *((Complex *)n) = (Complex)(x); break;\
    }

void lunum_tovalue(lua_State *L, ArrayType T, void *num)
{
  int isnum;

  ArrayAllNum tmp_num;

  if (tmp_num.li = lua_tointegerx(L, -1, &isnum), isnum) {
    ASSIGN_TO_VOID(T, num, tmp_num.li);
  } else if (tmp_num.ln = lua_tonumberx(L, -1, &isnum), isnum) {
    ASSIGN_TO_VOID(T, num, tmp_num.ln);
  }
  else if (lua_isboolean(L, -1)) {
    tmp_num.b = lua_toboolean(L, -1);
    ASSIGN_TO_VOID(T, num, tmp_num.b);
  }
  else if (lunum_hasmetatable(L, -1, "complex")) {
    tmp_num.z = *((Complex*) lua_touserdata(L, -1));
    ASSIGN_TO_VOID(T, num, tmp_num.z);
  }
  else {
    luaL_error(L, "unkown data type");
  }
}

