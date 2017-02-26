CC=gcc
CFLAGS=-O3 -std=gnu99 -Wall -Wpedantic -Wextra
LDLIBS=-lm

SRC=main.c zeta.c
OBJ=$(SRC:.c=.o)
DEPS=$(OBJ:.o=.d)

TARGET=zeta

all: $(TARGET)

$(TARGET): $(OBJ)
	$(CC) $(LDLIBS) -o $@ $^

%.o: %.c
	$(CC) -c $(CFLAGS) -MMD -o $@ $<

-include $(DEPS)

.PHONY: clean

clean:
	rm -f *.o *.d *.txt $(TARGET)
