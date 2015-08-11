

LUA = $(LUA_HOME)/bin/lua

default : $(LUNUM_A) $(LUNUM_SO)

%.o : %.c
	$(CC) $(CFLAGS) -std=c99 -c $<

array_class.lc : array_class.lua
	$(LUA) bin2c.lua +$< > $@

lunum.o : lunum.c array_class.lc
	$(CC) $(CFLAGS) -std=c99 -c $< -I$(LUA_HOME)/include

lunum_capi.o : lunum_capi.c
	$(CC) $(CFLAGS) -std=c99 -c $< -I$(LUA_HOME)/include

$(LUNUM_A)  : numarray.o lunum.o lunum_capi.o
	$(AR) $@ $?

$(LUNUM_SO) : numarray.o lunum.o lunum_capi.o
	$(SO) -o $(LUNUM_SO) $^

clean :
	rm -f *.o $(LUNUM_A) $(LUNUM_SO) *.lc
