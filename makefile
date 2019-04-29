# the compiler: gcc for C program
CC = gcc

# compiler flags:
#  -g    adds debugging information to the executable file
#  -Wall turns on most, but not all, compiler warnings
CFLAGS  = -g -Wall

# library flags:
#	-lgmp	gmp.h
#	-lm 	math.h
LIB = -lgmp -lm

# C files
FILES = test.c factor.c form.c arithmetic.c

# headers
HEADERS = project.h

# the build target executable:
TARGET = project


all: $(FILES) $(HEADERS)
	$(CC) $(CFLAGS) -o $(TARGET) $(FILES) $(LIB)

clean: 
	rm -f $(TARGET)

run: $(TARGET)
	./$(TARGET)