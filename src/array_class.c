

#define LUNUM_PRIVATE_API
#include <stdlib.h>
#include <string.h>
#include "lunum.h"


static int luaC_array_dtype(lua_State *L);
static int luaC_array_dtypemin(lua_State *L);
static int luaC_array_dtypemax(lua_State *L);
static int luaC_array_shape(lua_State *L);
static int luaC_array_size(lua_State *L);
static int luaC_array_astable(lua_State *L);
static int luaC_array_astype(lua_State *L);
static int luaC_array_tofile(lua_State *L);

void _lunum_register_array(lua_State *L, Array *B)
{
  lua_newtable(L);

  lua_pushcfunction(L, luaC_array_dtype);
  lua_setfield(L, -2, "dtype");

  lua_pushcfunction(L, luaC_array_dtypemin);
  lua_setfield(L, -2, "dtypemin");

  lua_pushcfunction(L, luaC_array_dtypemax);
  lua_setfield(L, -2, "dtypemax");

  lua_pushcfunction(L, luaC_array_shape);
  lua_setfield(L, -2, "shape");

  lua_pushcfunction(L, luaC_array_size);
  lua_setfield(L, -2, "size");

  lua_pushcfunction(L, luaC_array_astable);
  lua_setfield(L, -2, "astable");

  lua_pushcfunction(L, luaC_array_astype);
  lua_setfield(L, -2, "astype");

  lua_pushcfunction(L, luaC_array_tofile);
  lua_setfield(L, -2, "tofile");


  // Calls code written in Lua: lunum.__register_array(new_array) to add
  // additional methods to the new array.
  // ---------------------------------------------------------------------------
  lua_getglobal(L, "lunum");
  lua_getfield(L, -1, "__register_array");
  lua_pushvalue(L, -3);
  lua_call(L, 1, 0);
  lua_pop(L, 1);


  void *data = lua_newuserdata(L, B->size * array_sizeof(B->dtype));
  lua_setfield(L, -2, "__buffer");

  memcpy(data, B->data, B->size * array_sizeof(B->dtype));
  free(B->data);
  B->owns = 0;
  B->data = data;

  Array *A = (Array*) lua_newuserdata(L, sizeof(Array));
  lua_setfield(L, -2, "__cstruct");
  *A = *B;

  luaL_setmetatable(L, "array");
}




int luaC_array_dtype(lua_State *L)
// -----------------------------------------------------------------------------
// If there is no argument, return a string description of the data type. If
// the string 'enum' is given as the first argument, then return the enumated
// value of the Array's type.
// -----------------------------------------------------------------------------
{
  Array *A = lunum_checkarray1(L, 1);

  if (lua_isstring(L, 2)) {
    if (strcmp(lua_tostring(L, 2), "enum") == 0) {
      lua_pushinteger(L, A->dtype);
      return 1;
    }
  }

  lua_pushstring(L, array_typename(A->dtype));
  return 1;
}

int luaC_array_dtypemin(lua_State *L)
// -----------------------------------------------------------------------------
// If there is no argument, return a string description of the data type. If
// the string 'enum' is given as the first argument, then return the enumated
// value of the Array's type.
// -----------------------------------------------------------------------------
{
  Array *A = lunum_checkarray1(L, 1);

  lua_pushnumber(L, array_typemin(A->dtype));
  return 1;
}

int luaC_array_dtypemax(lua_State *L)
// -----------------------------------------------------------------------------
// If there is no argument, return a string description of the data type. If
// the string 'enum' is given as the first argument, then return the enumated
// value of the Array's type.
// -----------------------------------------------------------------------------
{
  Array *A = lunum_checkarray1(L, 1);

  lua_pushnumber(L, array_typemax(A->dtype));
  return 1;
}

int luaC_array_shape(lua_State *L)
// -----------------------------------------------------------------------------
// If there is no argument, return the shape as a table. If the string 'array'
// is given, return it as an array.
// -----------------------------------------------------------------------------
{
  Array *A = lunum_checkarray1(L, 1);
  lunum_pusharray2(L, A->shape, ARRAY_TYPE_SIZE_T, (size_t)A->ndims);

  if (lua_isstring(L, 2)) {
    if (strcmp(lua_tostring(L, 2), "array") == 0) {
      return 1;
    }
  }

  lunum_astable(L, 2);
  lua_replace(L, -2);
  return 1;
}

int luaC_array_size(lua_State *L)
{
  Array *A = lunum_checkarray1(L, 1);
  lua_pushinteger(L, A->size);
  return 1;
}

int luaC_array_astable(lua_State *L)
{
  lunum_astable(L, 1);
  return 1;
}

int luaC_array_astype(lua_State *L)
{
  Array *A = lunum_checkarray1(L, 1);
  ArrayType T;

  if (lua_type(L, 2) == LUA_TSTRING) {
    T = array_typeflag(lua_tostring(L, 2)[0]);
  }
  else {
    T = (ArrayType) luaL_checkinteger(L, 2);
  }

  Array B = array_new_copy(A, T);
  lunum_pusharray1(L, &B);
  return 1;
}



int luaC_array_tofile(lua_State *L)
// -----------------------------------------------------------------------------
// Writes the array 'A' as binary data to the file named 'fname'.
// -----------------------------------------------------------------------------
{
  Array *A = lunum_checkarray1(L, 1);
  const char *fname = luaL_checkstring(L, 2);
  FILE *output = fopen(fname, "wb");

  if (output == NULL) {
    luaL_error(L, "could not create file %s", fname);
  }
  fwrite(A->data, A->size, array_sizeof(A->dtype), output);
  fclose(output);

  return 0;
}
