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
FILES = src/main.c src/cgm.c src/form.c src/auxiliary.c
FILES_FORM = src/main-form.c src/form.c src/auxiliary.c
# headers
HEADERS = headers/cgm.h

# the build target executable:
TARGET = cgm

# build test for form
FORM = test_form

all: $(FILES) $(HEADERS)
	$(CC) $(CFLAGS) -o $(TARGET) $(FILES) $(LIB)

form-only: $(FILES_FORM)
	$(CC) $(CFLAGS) -o $(FORM) $(FILES_FORM) $(LIB)
clean: 
	rm -f $(TARGET)

run: $(TARGET)
	./$(TARGET)
