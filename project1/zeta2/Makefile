CC=gcc
CFLAGS=-O3 -std=gnu99 -Wall -Wpedantic -Wextra -fopenmp
LDFLAGS=-fopenmp
LDLIBS=-lm

SRC=main.c test1.c test2.c zeta.c
OBJ=$(SRC:.c=.o)
COMMON_OBJ=zeta.o
DEPS=$(OBJ:.o=.d)

TARGET=zeta
TEST1=test1
TEST2=test2

all: $(TARGET) $(TEST1) $(TEST2)

$(TARGET): main.o $(COMMON_OBJ)
	$(CC) $(LDFLAGS) $(LDLIBS) -o $@ $^

$(TEST1): test1.o $(COMMON_OBJ)
	$(CC) $(LDFLAGS) $(LDLIBS) -o $@ $^

$(TEST2): test2.o $(COMMON_OBJ)
	$(CC) $(LDFLAGS) $(LDLIBS) -o $@ $^

%.o: %.c
	$(CC) -c $(CFLAGS) -MMD -o $@ $<

-include $(DEPS)

.PHONY: clean utest vtest

utest: $(TEST1)
	./$(TEST1)

vtest: $(TEST2)
	./$(TEST2)

clean:
	rm -f *.o *.d *.txt $(TARGET) $(TEST1) $(TEST2)
