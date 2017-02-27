CC=mpicc
CFLAGS=-O3 -std=gnu99 -Wall -Wpedantic -Wextra
LDLIBS=-lm

SRC=main.c zeta.c test1.c
OBJ=$(SRC:.c=.o)
DEPS=$(OBJ:.o=.d)

TARGET=zeta
TEST1=test1
TARGET_RD=zeta_rd
TEST1_RD=test1_rd

all: $(TARGET) $(TEST1) $(TARGET_RD) $(TEST1_RD)

$(TARGET): main.o zeta.o
	$(CC) $(LDLIBS) -o $@ $^

$(TEST1): test1.o zeta.o
	$(CC) $(LDLIBS) -o $@ $^

$(TARGET_RD): main.o zeta_rd.o
	$(CC) $(LDLIBS) -o $@ $^

$(TEST1_RD): test1.o zeta_rd.o
	$(CC) $(LDLIBS) -o $@ $^

%.o: %.c
	$(CC) -c $(CFLAGS) -MMD -o $@ $<

zeta_rd.o: zeta.c
	$(CC) -c $(CFLAGS) -MMD -DRECURSIVE_DOUBLING -o $@ $<

-include $(DEPS)

.PHONY: clean vtest vtest_rd

clean:
	rm -f *.o *.d *.txt $(TARGET) $(TEST1)

vtest: $(TEST1)
	mpirun ./$<

vtest_rd: $(TEST1_RD)
	mpirun ./$<
