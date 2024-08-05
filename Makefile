MAKE=make
CC=gcc
LD=ld
AR=ar

OPTFLAGS=-O3 -g
CFLAGS=$(OPTFLAGS) -Wall
LDFLAGS=-lpthread -lm -std=c++17 -D_GLIBCXX_PARALLEL -fopenmp
MY_LIBS=trifac.o areanorm.o sphfunc.o ellfit.o covsrt.o ludcmp.o lubksb.o mrqmin.o mrqcof.o\
        curv.o blmatrix.o conv.o gauss.o phasec.o matrix.o bright.o memory.o\
	dot_product.o

all: convexinv yorp bootstrap yorp_e minkowski standardtri shape2obj

libs: $(MY_LIBS)

convexinv: convexinv.c $(MY_LIBS)
	$(CC) $(CFLAGS) -o $@ $< $(MY_LIBS) $(LDFLAGS) -g

yorp: yorp.cpp
	g++ yorp.cpp -g -std=c++11 -o yorp -lpthread

bootstrap: bootstrap.cpp lData.cpp lightCurve.cpp dataPoint.cpp point.cpp
	g++ $(CFLAGS) -o $@ $^ $(LDFLAGS)

yorp_e: yorp_e.cpp
	g++ yorp_e.cpp -g -std=c++11 -o yorp_e -lpthread

minkowski: minkowski.f
	gfortran ./minkowski.f  -o minkowski

standardtri: standardtri.f
	gfortran ./standardtri.f -o standardtri

shape2obj: shape2obj.cpp
	g++ shape2obj.cpp -o shape2obj

%.o: %.c
	$(CC) $(CFLAGS) -c $<

clean:
	rm -f *.o convexinv yorp bootstrap yorp_e minkowski standardtri shape2obj 
