CC=mpicc
CFLAGS=-O3 -std=gnu99 -Wall -Wpedantic -Wextra
LDLIBS=-lm

SRC=main.c mach.c test1.c
OBJ=$(SRC:.c=.o)
DEPS=$(OBJ:.o=.d)

TARGET=mach
TEST1=test1

all: $(TARGET) $(TEST1)

$(TARGET): main.o mach.o
	$(CC) $(LDLIBS) -o $@ $^

$(TEST1): test1.o mach.o
	$(CC) $(LDLIBS) -o $@ $^

%.o: %.c
	$(CC) -c $(CFLAGS) -MMD -o $@ $<

-include $(DEPS)

.PHONY: clean vtest

clean:
	rm -f *.o *.d *.txt $(TARGET) $(TEST1)

vtest: $(TEST1)
	mpirun ./$<
