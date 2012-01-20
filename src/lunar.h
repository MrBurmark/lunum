

#ifndef __NumluaCapi_HEADER__
#define __NumluaCapi_HEADER__

#include "numarray.h"
#include "lualib.h"

struct Array *lunar_checkarray1(lua_State *L, int pos);
void         *lunar_checkarray2(lua_State *L, int pos, enum ArrayType T, size_t *N);
void          lunar_pusharray1(lua_State *L, struct Array *A);
void          lunar_pusharray2(lua_State *L, void *data, enum ArrayType T, size_t N);
void          lunar_upcast(lua_State *L, int pos, enum ArrayType T, size_t N);

#endif // __NumluaCapi_HEADER__
