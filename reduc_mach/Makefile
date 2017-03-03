CC=mpicc
CFLAGS=-O3 -std=gnu99 -Wall -Wpedantic -Wextra
LDLIBS=-lm

SRC=main.c mach.c test1.c
OBJ=$(SRC:.c=.o)
DEPS=$(OBJ:.o=.d)

TARGET=mach
TEST1=test1
TARGET_RD=mach_rd
TEST1_RD=test1_rd

all: $(TARGET) $(TEST1) $(TARGET_RD) $(TEST1_RD)

$(TARGET): main.o mach.o
	$(CC) $(LDLIBS) -o $@ $^

$(TEST1): test1.o mach.o
	$(CC) $(LDLIBS) -o $@ $^

$(TARGET_RD): main.o mach_rd.o
	$(CC) $(LDLIBS) -o $@ $^

$(TEST1_RD): test1.o mach_rd.o
	$(CC) $(LDLIBS) -o $@ $^

%.o: %.c
	$(CC) -c $(CFLAGS) -MMD -o $@ $<

mach_rd.o: mach.c
	$(CC) -c $(CFLAGS) -MMD -DRECURSIVE_DOUBLING -o $@ $<

-include $(DEPS)

.PHONY: clean vtest vtest_rd

clean:
	rm -f *.o *.d *.txt $(TARGET) $(TEST1) $(TARGET_RD) $(TEST1_RD)

vtest: $(TEST1)
	mpirun ./$<

vtest_rd: $(TEST1_RD)
	mpirun ./$<
