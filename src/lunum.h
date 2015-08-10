

#ifndef __NumluaCapi_HEADER__
#define __NumluaCapi_HEADER__

#include "lua.h"
#include "lualib.h"
#include "lauxlib.h"

#include "numarray.h"

int luaopen_lunum(lua_State *L);

Array *lunum_checkarray1(lua_State *L, int pos);
/* returns a pointer to the data, and sets *N to the size */
void         *lunum_checkarray2(lua_State *L, int pos, ArrayType T, size_t *N);
void          lunum_pusharray1(lua_State *L, Array *B);
void          lunum_pusharray2(lua_State *L, void *data, ArrayType T, size_t N);
void          lunum_astable(lua_State *L, int pos);
int           lunum_upcast(lua_State *L, int pos, ArrayType T, size_t N);
int           lunum_hasmetatable(lua_State *L, int pos, const char *name);
void         *lunum_tovalue(lua_State *L, ArrayType T);

#ifndef LUNUM_API_NOCOMPLEX

#include <complex.h>
typedef double complex Complex;

Complex       lunum_checkcomplex(lua_State *L, int n);
void          lunum_pushcomplex(lua_State *L, Complex z);
#endif // LUNUM_API_COMPLEX

#ifdef LUNUM_PRIVATE_API
void         _lunum_register_array(lua_State *L, Array *B);
#endif // LUNUM_PRIVATE_API

#endif // __NumluaCapi_HEADER__
