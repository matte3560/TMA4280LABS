CC=gcc
CFLAGS=-O3 -std=gnu99 -Wall -Wpedantic -Wextra
LDLIBS=-lm

HEADER=zeta.h
SOURCE=main.c utest.c zeta.c
COMMON_OBJ=zeta.o

TARGET=zeta
TEST1=test1
TEST2=test2

all: $(TARGET) $(TEST1) $(TEST2)

$(TARGET): main.o $(COMMON_OBJ)
	$(CC) $(LDLIBS) -o $@ $^

$(TEST1): test1.o $(COMMON_OBJ)
	$(CC) $(LDLIBS) -o $@ $^

$(TEST2): test2.o $(COMMON_OBJ)
	$(CC) $(LDLIBS) -o $@ $^

$(SOURCE): $(HEADER)
	$(CC) $(CFLAGS) -c $<

.PHONY: clean utest vtest

utest: $(TEST1)
	./$(TEST1)

vtest: $(TEST2)
	./$(TEST2)

clean:
	rm -f *.o *.gch $(TARGET) $(TEST1) $(TEST2) *.txt
