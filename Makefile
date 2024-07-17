# Compiler
CC = gcc

# Compiler flags
CFLAGS = -ansi -Wall -Wextra -Werror -pedantic-errors

# Libraries to link
LIBS = -lm

# Source files
SRCS = symnmf.c utils.c

# Header files
HDRS = symnmf.h utils.h defs.h

# Object files
OBJS = $(SRCS:.c=.o)

# Output executable
OUT = symnmf

# Default target
all: $(OUT)

# Rule to link the executable
$(OUT): $(OBJS)
	$(CC) $(CFLAGS) -o $(OUT) $(OBJS) $(LIBS)

# Rule to compile source files into object files
%.o: %.c $(HDRS)
	$(CC) $(CFLAGS) -c $< -o $@

# Clean up the build
clean:
	rm -f $(OBJS) $(OUT)

.PHONY: all clean
