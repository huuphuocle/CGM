# the compiler: gcc for C program
CC = gcc

# compiler flags:
#  -g    adds debugging information to the executable file
#  -Wall turns on most, but not all, compiler warnings
CFLAGS  = -g -Wall

# library flags:
#	-lgmp	gmp.h
#	-lm 	math.h
LIB = -lflint -lgmp -lm

# C files
FILES = main.c factor.c form.c arithmetic.c
FILES_FORM = test_form.c form.c arithmetic.c
# headers
HEADERS = cgm.h

# the build target executable:
TARGET = cgm

# build test for form
FORM = test_form

all: $(FILES) $(HEADERS)
	$(CC) $(CFLAGS) -o $(TARGET) $(FILES) $(LIB)

form: $(FILES_FORM)
	$(CC) $(CFLAGS) -o $(FORM) $(FILES_FORM) $(LIB)
clean: 
	rm -f $(TARGET)

run: $(TARGET)
	./$(TARGET)