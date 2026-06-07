# Preliminary makefile

CFLAGS = -DWITH_SPLAYTREE -I../src -I../src/api -I../src/env -I../src/misc -I../src/mpl -I../nl-writer2/include -march=native -O2

OBJSET = \
mathprog.o \
api/mpl.o \
env/alloc.o \
env/dlsup.o \
env/env.o \
env/error.o \
env/stdc.o \
env/stream.o \
env/tls.o \
misc/avl.o \
misc/dmp.o \
misc/jd.o \
misc/splaytree.o \
misc/str2int.o \
misc/str2num.o \
misc/strspx.o \
misc/strtrim.o \
misc/rng.o \
misc/rng1.o \
mpl/mpl1.o \
mpl/mpl2.o \
mpl/mpl3.o \
mpl/mpl4.o \
mpl/mpl5.o \
mpl/mpl6.o \
mpl/mplsql.o \
../nl-writer2/src/nl-solver.o \
../nl-writer2/src/nl-utils.o \
../nl-writer2/src/nl-model-c.o \
../nl-writer2/src/nl-solver-c.o \
../nl-writer2/src/nl-writer2.o


.c.o:
	gcc $(CFLAGS) -c $*.c -o $*.o

.cc.o:
	gcc $(CFLAGS) -c $*.cc -o $*.o

all: mathprog

clean:
	  rm ../src/mathprog \
		    ../src/*.o \
		    ../src/api/*.o \
		    ../src/env/*.o \
		    ../src/misc/*.o \
		    ../src/mpl/*.o \
		    ../nl-writer2/src/*.o

mathprog: $(OBJSET)
	  gcc $(CFLAGS) $(OBJSET) -o mathprog -lglpk -lm -lz -lstdc++

check: mathprog
	  ./mathprog --version
	  ./mathprog -m ../examples/sudoku.mod

## eof ##
