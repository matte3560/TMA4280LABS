CC=gcc
CFLAGS=-O3 -std=c99 -Wall -Wpedantic -Wextra
LDLIBS=-lm

HEADER=zeta.h
SOURCE=main.c utest.c zeta.c
COMMON_OBJ=zeta.o

TARGET=zeta
TEST=utest

all: $(TARGET)

$(TARGET): main.o $(COMMON_OBJ)
	$(CC) $(LDLIBS) -o $@ $^

$(TEST): utest.o $(COMMON_OBJ)
	$(CC) $(LDLIBS) -o $@ $^

$(SOURCE): $(HEADER)
	$(CC) $(CFLAGS) -c $<

.PHONY: clean

clean:
	rm -f *.o $(TARGET) $(TEST) *.txt
