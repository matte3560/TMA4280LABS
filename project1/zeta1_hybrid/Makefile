CC=mpicc
CFLAGS=-O3 -std=gnu99 -Wall -Wpedantic -Wextra -fopenmp
LDFLAGS=-fopenmp
LDLIBS=-lm

SRC=main.c zeta.c test1.c
OBJ=$(SRC:.c=.o)
DEPS=$(OBJ:.o=.d)

TARGET=zeta
TEST1=test1

all: $(TARGET) $(TEST1)

$(TARGET): main.o zeta.o
	$(CC) $(LDFLAGS) $(LDLIBS) -o $@ $^

$(TEST1): test1.o zeta.o
	$(CC) $(LDFLAGS) $(LDLIBS) -o $@ $^

%.o: %.c
	$(CC) -c $(CFLAGS) -MMD -o $@ $<

-include $(DEPS)

.PHONY: clean vtest

clean:
	rm -f *.o *.d *.txt $(TARGET) $(TEST1)

vtest: $(TEST1)
	mpirun ./$<
