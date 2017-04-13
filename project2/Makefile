FC=gfortran
FFLAGS=-lm
CFLAGS=-std=c99 -O3 -MMD
LDFLAGS=-lm

all: serial

serial: serial.o serial_poisson.o common_poisson.o fst.o

-include *.d

.PHONY: clean

clean:
	$(RM) *.o *.d serial
