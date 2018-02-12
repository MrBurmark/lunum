
require 'lunum'

-- test helper routines

-- checks if every value in a table is true
local function check_table_true(t)
   for k,v in pairs(t) do
      if (v ~= true) then
         print('value not true in tablecheck', k, v)
      end
      assert(v == true)
   end
end

local function raw_pairs(t)
   return next, t, nil
end

-- returns the higher typed of A, B
local function higher_typed(A, B)
   if (lunum[A:dtype()] >= lunum[B:dtype()]) then
      return A
   else
      return B
   end
end
-- returns the type max of the higher type of A, B
local function type_max(A, B)
   return higher_typed(A, B):dtypemax()
end
-- returns the type min of the higher type of A, B
local function type_min(A, B)
   return higher_typed(A, B):dtypemin()
end
-- returns true if lout is not representable in the higher typed of A, B
local function under_over_flow(A, B, lout)
   if type(A) == 'userdata' and getmetatable(A).__name == 'array' and type(B) == 'userdata' and getmetatable(B).__name == 'array' then
      return (lout > type_max(A, B)) or (lout < type_min(A, B))
   elseif type(A) == 'userdata' and getmetatable(A).__name == 'array' then
      return (lout > A:dtypemax()) or (lout < A:dtypemin())
   elseif type(B) == 'userdata' and getmetatable(B).__name == 'array' then
      return (lout > B:dtypemax()) or (lout < B:dtypemin())
   else
      return false
   end
end
-- returns true if both args are nan
local function both_nan(a, b)
   return (a ~= a) and (b ~= b)
end

local function int_typed(...)
   local ints = true
   for _,v in pairs({...}) do
      if type(v) == 'userdata' and getmetatable(v).__name == 'array' then
         ints = ints and lunum[v:dtype()] <= lunum.size_t
      elseif tonumber(v) then
         ints = ints and (math.floor(tonumber(v)) == tonumber(v))
      else
         ints = false
      end
   end
   return ints
end

-- tests for binary operations and unary operations
-- take a test name, index into the array being tested, Arrays, the array result at i and the lua result at i
-- doesn't print failure if overflow or underflow
local function test_binary_op(name, i, A, B, aout, lout)
   if int_typed(A, B, aout) then
      lout = lout >= 0 and math.floor(lout) or math.ceil(lout)
   end
   if (aout[i] ~= lout and not under_over_flow(A, B, lout) and not both_nan(aout[i], lout)) then
      if type(A) == 'userdata' and getmetatable(A).__name == 'array' and type(B) == 'userdata' and getmetatable(B).__name == 'array' then
         print('binary op failed', name, i, A:dtype(), B:dtype(), A[i], B[i], aout:dtype(), aout[i], '~=', lout)
      elseif type(A) == 'userdata' and getmetatable(A).__name == 'array' then
         print('binary op failed', name, i, A:dtype(), tonumber(B) and 'number' or 'unknown', A[i], B, aout:dtype(), aout[i], '~=', lout)
      elseif type(B) == 'userdata' and getmetatable(B).__name == 'array' then
         print('binary op failed', name, i, tonumber(A) and 'number' or 'unknown', B:dtype(), A, B[i], aout:dtype(), aout[i], '~=', lout)
      else
         print('binary op failed', name, i, tonumber(A) and 'number' or 'unknown', tonumber(B) and 'number' or 'unknown', A, B, aout:dtype(), aout[i], '~=', lout)
      end
   end
end
local function test_unary_op(name, i, A, aout, lout)
   if (aout[i] ~= lout and not both_nan(aout[i], lout)) then
      print('unary op failed', name, i, A:dtype(), A[i], aout:dtype(), aout[i], '~=', lout)
   end
end

-- tests that all the functions and constants are in the array metatable and the lunum table
local function test1()
   print('test1')
   local A = lunum.zeros(100)

   -- check array object table
   assert(type(A) == 'userdata')

   -- check array object metatable
   local expected_metafunctions = {__gc=false, __div=false, __index=false, __band=false, __add=false,
                                    __mod=false, __sub=false, __call=false, __unm=false, __pow=false,
                                    __mul=false, __shl=false, __bxor=false, __idiv=false, __bor=false,
                                    __bnot=false, __newindex=false, __tostring=false, __shr=false,
                                    dtypemin=false, dtypemax=false, dtype=false, shape=false, size=false,
                                    astable=false, astype=false, tofile=false, max=false, imag=false,
                                    eq=false, min=false, copy=false, real=false, indices=false, ne=false,
                                    resize=false, setasflat=false, reshape=false, gt=false, ge=false,
                                    le=false, conj=false, lt=false, __preserve=false}

   local expected_metastrings = {__name=false}

   for k,v in pairs(getmetatable(A)) do
      if (expected_metafunctions[k] == false) then
         assert(type(v) == 'function')
         expected_metafunctions[k] = true
      elseif (expected_metastrings[k] == false) then
         assert(type(v) == 'string')
         expected_metastrings[k] = true
      else
         print('Unknown pair in array metatable', k, v)
         assert(false)
      end
   end

   check_table_true(expected_metafunctions)
   check_table_true(expected_metastrings)


   -- check complex metatable
   local C = lunum.I;

   assert(type(C) == 'userdata')

   local expected_complex_metafunctions = {__div=false, __add=false, __sub=false, __unm=false, __pow=false,
                                    __mul=false, __tostring=false, __eq=false, __lt=false, __le=false,
                                    new=false, __preserve=false}

   local expected_complex_metastrings = {__name=false}

   for k,v in pairs(getmetatable(C)) do
      if (expected_complex_metafunctions[k] == false) then
         assert(type(v) == 'function')
         expected_complex_metafunctions[k] = true
      elseif (expected_complex_metastrings[k] == false) then
         assert(type(v) == 'string')
         expected_complex_metastrings[k] = true
      else
         print('Unknown pair in complex metatable', k, v)
         assert(false)
      end
   end

   check_table_true(expected_complex_metafunctions)
   check_table_true(expected_complex_metastrings)


   -- check lunum table
   local expected_functions = {tanh=false, atanh=false, __register_array_metafunctions=false,
                              sinh=false, asin=false, __build_slice=false, acosh=false, zeros=false,
                              conjugate=false, array=false, resize=false, sin=false, cos=false,
                              fromfile=false, exp=false, log=false, log10=false, range=false,
                              asinh=false, apply=false, slice=false, cosh=false, acos=false,
                              transpose=false, atan=false, tan=false, loadtxt=false, linear=false}

   local expected_numbers = {char=false, short=false, double=false, long=false, int=false,
                              complex=false, size_t=false, float=false, bool=false}

   local expected_userdata = {I=false}

   for k,v in pairs(lunum) do
      if (expected_functions[k] == false) then
         assert(type(v) == 'function')
         expected_functions[k] = true
      elseif (expected_numbers[k] == false) then
         assert(type(v) == 'number')
         expected_numbers[k] = true
      elseif (expected_userdata[k] == false) then
         assert(type(v) == 'userdata')
         expected_userdata[k] = true
      else
         print('Unknown pair in lunum', k, v)
         assert(false)
      end
   end

   check_table_true(expected_functions)
   check_table_true(expected_numbers)
   check_table_true(expected_userdata)

   assert(A:size() == 100)
end

-- checks if array creation from tables works, simple assignment works
local function test2()
   print('test2')
   local t = {0,1,2,3,4,5,6,7,8,9}
   local A = lunum.array(t)
   for k,v in pairs(t) do
      if (A[k-1] ~= t[k]) then
         print('creation from array failed at', k, A[k-1], "~=", t[k])
      end
      assert(A[k-1] == t[k])
   end
   A[4] = 4.5
   assert(A[4] == 4.5)
end

-- check basic array array operators
local function test3()
   print('test3')
   local A = lunum.array({0,1,2,3,4,5,6,7,8,9})
   local B = lunum.array({9,8,7,6,5,4,3,2,1,0})

   for i = 0,A:size()-1 do
      test_binary_op('+', i, A, B, (A + B), A[i] + B[i])
   end
   for i = 0,A:size()-1 do
      test_binary_op('-', i, A, B, (A - B), A[i] - B[i])
   end
   for i = 0,A:size()-1 do
      test_binary_op('*', i, A, B, (A * B), A[i] * B[i])
   end
   for i = 0,A:size()-1 do
      test_binary_op('/', i, A, B, (A / B), A[i] / B[i])
   end
   -- for i = 0,A:size()-1 do
   --    test_binary_op('//', i, A, B, (A // B), A[i] // B[i])
   -- end
   for i = 0,A:size()-1 do
      -- some results nan
      test_binary_op('%', i, A, B, (A % B), A[i] % B[i])
   end
   for i = 0,A:size()-1 do
      test_binary_op('^', i, A, B, (A ^ B), A[i] ^ B[i])
   end

   for i = 0,A:size()-1 do
      test_unary_op('-', i, A, (-A), -A[i])
   end

   A = A:astype(lunum.long)
   B = A:astype(lunum.long)

   -- for i = 0,A:size()-1 do
   --    test_binary_op('&', i, A, B, (A & B), A[i] & B[i])
   -- end
   -- for i = 0,A:size()-1 do
   --    test_binary_op('|', i, A, B, (A | B), A[i] | B[i])
   -- end
   -- for i = 0,A:size()-1 do
   --    test_binary_op('~', i, A, B, (A ~ B), A[i] ~ B[i])
   -- end
   -- for i = 0,A:size()-1 do
   --    test_binary_op('<<', i, A, B, (A << B), A[i] << B[i])
   -- end
   -- for i = 0,A:size()-1 do
   --    test_binary_op('>>', i, A, B, (A >> B), A[i] >> B[i])
   -- end

   -- for i = 0,A:size()-1 do
   --    test_unary_op('~', i, A, (~A), ~A[i])
   -- end
end

-- test operations on arrays of different types
local function test4()
   print('test4')
   -- no testing with booleans as lua doesn't allow boolean arithmetic
   -- tests use -5 to -1, 1 to 8 to avoid division by 0, loss of float precision
   local t = {
      lunum.array({1,2,3,4,5,6,7,8}, lunum.char),
      lunum.array({1,2,3,4,5,6,7,8}, lunum.short),
      lunum.array({1,2,3,4,5,6,7,8}, lunum.int),
      lunum.array({1,2,3,4,5,6,7,8}, lunum.long),
      lunum.array({1,2,3,4,5,6,7,8}, lunum.size_t)
   }

   -- integer only tests
   for i = 1,#t do
      for j = 1,#t do
         -- for k = 0,t[i]:size()-1 do
         --    test_binary_op('&', k, t[i], t[j], (t[i] & t[j]), t[i][k] & t[j][k])
         -- end
         -- for k = 0,t[i]:size()-1 do
         --    test_binary_op('|', k, t[i], t[j], (t[i] | t[j]), t[i][k] | t[j][k])
         -- end
         -- for k = 0,t[i]:size()-1 do
         --    test_binary_op('~', k, t[i], t[j], (t[i] ~ t[j]), t[i][k] ~ t[j][k])
         -- end
         -- for k = 0,t[i]:size()-1 do
         --    test_binary_op('<<', k, t[i], t[j], (t[i] << t[j]), t[i][k] << t[j][k])
         -- end
         -- for k = 0,t[i]:size()-1 do
         --    test_binary_op('>>', k, t[i], t[j], (t[i] >> t[j]), t[i][k] >> t[j][k])
         -- end
      end
   end

   t[#t+1] = lunum.array({1,2,3,4,5,6,7,8}, lunum.float)
   t[#t+1] = lunum.array({1,2,3,4,5,6,7,8}, lunum.double)

   -- real only tests
   for i = 1,#t do
      for j = 1,#t do
         -- for k = 0,t[i]:size()-1 do
         --    test_binary_op('//', k, t[i], t[j], (t[i] // t[j]), t[i][k] // t[j][k])
         -- end
         for k = 0,t[i]:size()-1 do
            test_binary_op('%', k, t[i], t[j], (t[i] % t[j]), t[i][k] % t[j][k])
         end
      end
   end

   local t = {
      lunum.array({-5,-4,-3,-2,-1,1,2,3,4,5,6,7,8}, lunum.char),
      lunum.array({-5,-4,-3,-2,-1,1,2,3,4,5,6,7,8}, lunum.short),
      lunum.array({-5,-4,-3,-2,-1,1,2,3,4,5,6,7,8}, lunum.int),
      lunum.array({-5,-4,-3,-2,-1,1,2,3,4,5,6,7,8}, lunum.long),
      lunum.array({ 1, 1, 1, 1, 1,1,2,3,4,5,6,7,8}, lunum.size_t),
      lunum.array({-5,-4,-3,-2,-1,1,2,3,4,5,6,7,8}, lunum.float),
      lunum.array({-5,-4,-3,-2,-1,1,2,3,4,5,6,7,8}, lunum.double),
      lunum.array({-5,-4,-3,-2,-1,1,2,3,4,5,6,7,8}, lunum.complex)
   }

   -- real and complex tests
   for i = 1,#t do
      for j = 1,#t do
         for k = 0,t[i]:size()-1 do
            test_binary_op('+', k, t[i], t[j], (t[i] + t[j]), t[i][k] + t[j][k])
         end
         for k = 0,t[i]:size()-1 do
            test_binary_op('-', k, t[i], t[j], (t[i] - t[j]), t[i][k] - t[j][k])
         end
         for k = 0,t[i]:size()-1 do
            test_binary_op('*', k, t[i], t[j], (t[i] * t[j]), t[i][k] * t[j][k])
         end
         for k = 0,t[i]:size()-1 do
            test_binary_op('/', k, t[i], t[j], (t[i] / t[j]), t[i][k] / t[j][k])
         end
         for k = 0,t[i]:size()-1 do
            test_binary_op('^', k, t[i], t[j], (t[i] ^ t[j]), t[i][k] ^ t[j][k])
         end
      end
   end
end

-- test operations with array and constant
local function test5()
   print('test5')
   -- no testing with booleans as lua doesn't allow boolean arithmetic
   -- tests use -5 to -1, 1 to 8 to avoid division by 0, loss of float precision
   local t = {
      lunum.array({1,2,3,4,5,6,7,8}, lunum.char),
      lunum.array({1,2,3,4,5,6,7,8}, lunum.short),
      lunum.array({1,2,3,4,5,6,7,8}, lunum.int),
      lunum.array({1,2,3,4,5,6,7,8}, lunum.long),
      lunum.array({1,2,3,4,5,6,7,8}, lunum.size_t)
   }

   -- integer only tests
   for i = 1,#t do
      for j = 1,#t do
         -- for k = 0,t[i]:size()-1 do
         --    test_binary_op('&', k, t[i], t[j][k], (t[i] & t[j][k]), t[i][k] & t[j][k])
         -- end
         -- for k = 0,t[i]:size()-1 do
         --    test_binary_op('|', k, t[i], t[j][k], (t[i] | t[j][k]), t[i][k] | t[j][k])
         -- end
         -- for k = 0,t[i]:size()-1 do
         --    test_binary_op('~', k, t[i], t[j][k], (t[i] ~ t[j][k]), t[i][k] ~ t[j][k])
         -- end
         -- for k = 0,t[i]:size()-1 do
         --    test_binary_op('<<', k, t[i], t[j][k], (t[i] << t[j][k]), t[i][k] << t[j][k])
         -- end
         -- for k = 0,t[i]:size()-1 do
         --    test_binary_op('>>', k, t[i], t[j][k], (t[i] >> t[j][k]), t[i][k] >> t[j][k])
         -- end
         -- for k = 0,t[i]:size()-1 do
         --    test_binary_op('&', k, t[i][k], t[j], (t[i][k] & t[j]), t[i][k] & t[j][k])
         -- end
         -- for k = 0,t[i]:size()-1 do
         --    test_binary_op('|', k, t[i][k], t[j], (t[i][k] | t[j]), t[i][k] | t[j][k])
         -- end
         -- for k = 0,t[i]:size()-1 do
         --    test_binary_op('~', k, t[i][k], t[j], (t[i][k] ~ t[j]), t[i][k] ~ t[j][k])
         -- end
         -- for k = 0,t[i]:size()-1 do
         --    test_binary_op('<<', k, t[i][k], t[j], (t[i][k] << t[j]), t[i][k] << t[j][k])
         -- end
         -- for k = 0,t[i]:size()-1 do
         --    test_binary_op('>>', k, t[i][k], t[j], (t[i][k] >> t[j]), t[i][k] >> t[j][k])
         -- end
      end
   end

   t[#t+1] = lunum.array({1,2,3,4,5,6,7,8}, lunum.float)
   t[#t+1] = lunum.array({1,2,3,4,5,6,7,8}, lunum.double)

   -- real only tests
   for i = 1,#t do
      for j = 1,#t do
         -- for k = 0,t[i]:size()-1 do
         --    test_binary_op('//', k, t[i], t[j][k], (t[i] // t[j][k]), t[i][k] // t[j][k])
         -- end
         for k = 0,t[i]:size()-1 do
            test_binary_op('%', k, t[i], t[j][k], (t[i] % t[j][k]), t[i][k] % t[j][k])
         end
         -- for k = 0,t[i]:size()-1 do
         --    test_binary_op('//', k, t[i][k], t[j], (t[i][k] // t[j]), t[i][k] // t[j][k])
         -- end
         for k = 0,t[i]:size()-1 do
            test_binary_op('%', k, t[i][k], t[j], (t[i][k] % t[j]), t[i][k] % t[j][k])
         end
      end
   end

   local t = {
      lunum.array({-5,-4,-3,-2,-1,1,2,3,4,5,6,7,8}, lunum.char),
      lunum.array({-5,-4,-3,-2,-1,1,2,3,4,5,6,7,8}, lunum.short),
      lunum.array({-5,-4,-3,-2,-1,1,2,3,4,5,6,7,8}, lunum.int),
      lunum.array({-5,-4,-3,-2,-1,1,2,3,4,5,6,7,8}, lunum.long),
      lunum.array({ 1, 1, 1, 1, 1,1,2,3,4,5,6,7,8}, lunum.size_t),
      lunum.array({-5,-4,-3,-2,-1,1,2,3,4,5,6,7,8}, lunum.float),
      lunum.array({-5,-4,-3,-2,-1,1,2,3,4,5,6,7,8}, lunum.double),
      lunum.array({-5,-4,-3,-2,-1,1,2,3,4,5,6,7,8}, lunum.complex)
   }

   -- real and complex tests
   for i = 1,#t do
      for j = 1,#t do
         for k = 0,t[i]:size()-1 do
            test_binary_op('+', k, t[i], t[j][k], (t[i] + t[j][k]), t[i][k] + t[j][k])
         end
         for k = 0,t[i]:size()-1 do
            test_binary_op('-', k, t[i], t[j][k], (t[i] - t[j][k]), t[i][k] - t[j][k])
         end
         for k = 0,t[i]:size()-1 do
            test_binary_op('*', k, t[i], t[j][k], (t[i] * t[j][k]), t[i][k] * t[j][k])
         end
         for k = 0,t[i]:size()-1 do
            test_binary_op('/', k, t[i], t[j][k], (t[i] / t[j][k]), t[i][k] / t[j][k])
         end
         for k = 0,t[i]:size()-1 do
            test_binary_op('^', k, t[i], t[j][k], (t[i] ^ t[j][k]), t[i][k] ^ t[j][k])
         end
         for k = 0,t[i]:size()-1 do
            test_binary_op('+', k, t[i][k], t[j], (t[i][k] + t[j]), t[i][k] + t[j][k])
         end
         for k = 0,t[i]:size()-1 do
            test_binary_op('-', k, t[i][k], t[j], (t[i][k] - t[j]), t[i][k] - t[j][k])
         end
         for k = 0,t[i]:size()-1 do
            test_binary_op('*', k, t[i][k], t[j], (t[i][k] * t[j]), t[i][k] * t[j][k])
         end
         for k = 0,t[i]:size()-1 do
            test_binary_op('/', k, t[i][k], t[j], (t[i][k] / t[j]), t[i][k] / t[j][k])
         end
         for k = 0,t[i]:size()-1 do
            test_binary_op('^', k, t[i][k], t[j], (t[i][k] ^ t[j]), t[i][k] ^ t[j][k])
         end
      end
   end
end

-- test complex assignment
local function test6()
   print('test6')
   local I = lunum.I
   local A = lunum.zeros(109, lunum.complex)
   A[4] = 1
   A[5] = I
   assert(A[4] == 1+0*I)
   assert(A[5] == I)
end

-- test other math functions
local function test7()
   print('test7')
   local I = lunum.I
   local A = lunum.array({1,2,3})

   for i = 0,A:size()-1 do
      test_unary_op('sin', i, A, lunum.sin(A), math.sin(A[i]))
   end
   for i = 0,A:size()-1 do
      test_unary_op('cos', i, A, lunum.cos(A), math.cos(A[i]))
   end
   for i = 0,A:size()-1 do
      test_unary_op('log', i, A, lunum.log(A), math.log(A[i]))
   end
   for i = 0,A:size()-1 do
      test_unary_op('atan', i, A, lunum.atan(A), math.atan(A[i]))
   end
end

-- test resizing
local function test8()
   print('test8')
   local A = lunum.array({1,2,3,4,5,6,7,8})
   local B = A:copy()
   lunum.resize(B, {4,2})
   local n = 0
   for I in B:indices('table') do
      assert(B[I] == A[n])
      n = n + 1
   end
   lunum.resize(B, {2,4})
   local i = 0
   for I in B:indices('table') do
      assert(B[I] == A[i])
      i = i + 1
   end
end

-- test indices iterators
local function test9()
   print('test9')
   local A = lunum.array({1,2,3,4,5,6,7,8,9,10,11,12})
   B = A:copy()
   lunum.resize(B, {2,1,3,2})

   local x = 0
   for i,j,k,l in B:indices() do
      assert(B(i,j,k,l) == A[x])
      x = x + 1
   end

   local i = 0
   for I in B:indices('table') do
      assert(B[I] == A[i])
      i = i + 1
   end
end

-- test apply, resize
local function test10()
   print('test10')

   local B = lunum.range(100)
   local C = lunum.apply(function(x,y,z) return x+y+z end, B, B, B)

   local B = lunum.zeros({10,10})
   local C = lunum.zeros({10,10}, lunum.complex)

   B[{1,1}] = 1
   C[{1,1}] = lunum.I

   local D = lunum.apply(function(x,y) return x+y end, B, C)
   assert(D:dtype() == 'complex')

   local C = lunum.zeros({4,4}, lunum.complex)
   assert(C:shape('array')[0] == 4)
   assert(C:shape('array')[1] == 4)

   lunum.resize(C, {4,4})
   assert(C:shape('array')[0] == 4)
   assert(C:shape('array')[1] == 4)
end

-- test bool arrays and logical operations
local function test11()
   print('test11')
   local C = lunum.array({0,1,true,false}, lunum.bool)
   assert(C[0] == false)
   assert(C[1] == true)
   assert(C[2] == true)
   assert(C[3] == false)

   local A = lunum.array({1,2,3})
   local B = lunum.array({3,2,1})

   local t = A:eq(B)
   assert(t[0] == false)
   assert(t[1] == true)
   assert(t[2] == false)

   local t = A:ne(B)
   assert(t[0] == true)
   assert(t[1] == false)
   assert(t[2] == true)

   local t = A:lt(B)
   assert(t[0] == true)
   assert(t[1] == false)
   assert(t[2] == false)

   local t = A:le(B)
   assert(t[0] == true)
   assert(t[1] == true)
   assert(t[2] == false)

   local t = A:gt(B)
   assert(t[0] == false)
   assert(t[1] == false)
   assert(t[2] == true)

   local t = A:ge(B)
   assert(t[0] == false)
   assert(t[1] == true)
   assert(t[2] == true)
end

-- test type flags
local function test12()
   print('test12')
   local I = lunum.I
   local C = lunum.array({1,2,3}, 'd')

   assert(C:dtype() == 'double')
   local z = C:astype('z')
   assert(z:dtype() == 'complex')
   assert(z[0] == 1 + 0*I)
   assert(z[1] == 2 + 0*I)
   assert(z[2] == 3 + 0*I)
end

-- test transpose
local function test13()
   print('test13')
   local A = lunum.range(12)
   local B = lunum.transpose(A)

   for i in A:indices() do
      assert(A[i] == B[i])
   end

   lunum.resize(A, {3,4})
   local B = lunum.transpose(A)
   for i,j in A:indices() do
      assert(A(i,j) == B(j,i))
   end

   lunum.resize(A, {2,3,2})
   local B = lunum.transpose(A)
   for i,j,k in A:indices() do
      assert(A(i,j,k) == B(k,j,i))
   end

   lunum.resize(A, {3,2,2})
   local B = lunum.transpose(A)
   for i,j,k in A:indices() do
      assert(A(i,j,k) == B(k,j,i))
   end

   lunum.resize(A, {3,1,2,2})
   local B = lunum.transpose(A)
   for i,j,k,l in A:indices() do
      assert(A(i,j,k,l) == B(l,k,j,i))
   end
end

-- test linear
local function test14()
   print('test14')
   local A = lunum.linear(-2, 2, 5, lunum.char)
   local B = lunum.linear(3, -4, 8)
   local C = lunum.linear(-math.pi/2, math.exp(1), 13)

   local n = -2
   for i in A:indices() do
      assert(A[i] == n)
      n = n + 1
   end

   local n = 3
   for i in B:indices() do
      assert(B[i] == n)
      n = n - 1
   end

   for i in C:indices() do
      local n = ((C:size() - 1 - i) * (-math.pi/2) + i * math.exp(1)) / 12
      assert(C[i] == n)
   end
end

test1()
test2()
test3()
test4()
test5()
test6()
test7()
test8()
test9()
test10()
test11()
test12()
test13()
test14()
