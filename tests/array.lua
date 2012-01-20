
require 'numlua'

local function test1()
   local A = numlua.zeros(100)

   for k,v in pairs(getmetatable(A)) do
      print (k,v)
   end

   for k,v in pairs(numlua) do
      print (k,v)
   end

   print(getmetatable(A).__len)
   print(getmetatable(A).__gc)
   print("length is", #A)
end

local function test2()
   local A = numlua.array({0,1,2,3,4,5,6,7,8,9})
   A[4] = 4.5
   print(A[4])
end

local function test3()
   local A = numlua.array({0,1,2,3,4,5,6,7,8,9})
   local B = numlua.array({9,8,7,6,5,4,3,2,1,0})
   print("(A + B)[3] = ", (A + B)[2])
   print("(A - B)[3] = ", (A - B)[2])
   print("(A * B)[3] = ", (A * B)[2])
   print("(A / B)[3] = ", (A / B)[2])
end

local function test4()
   local A = numlua.array({0,1,2,3,4,5,6,7,8,9}, numlua.float)
   local B = numlua.array({0,1,2,3,4,5,6,7,8,9}, numlua.double)
   print("(A + B)[3] = ", (A + B)[2])
end

local function test5()
   print("[char]    A*10 = ", numlua.array({0,1,2,3,4,5,6,7,8,9}, numlua.char) * 10)
   print("[float]   A    = ", numlua.array({0,1,2,3,4,5,6,7,8,9}, numlua.float))
   print("[complex] A    = ", numlua.array({0,1,2,3,4,5,6,7,8,9}, numlua.complex))
end

local function test6()
   local I = numlua.I
   local A = numlua.zeros(109, numlua.complex)
   A[4] = 1
   A[5] = I
   print(A[4], A[5])
   print(A)
end


test1()
test2()

test3()
test4()
test5()
test6()
