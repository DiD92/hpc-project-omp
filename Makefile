SRC=convolution.c
LIB=convolution.h
DST=convolution-omp
OPTS=-std=c99 -Wall -Wextra -pedantic-errors -O2
INC=-fopenmp
CC=gcc

conv-omp: $(SRC) $(LIB)
	@$(CC) $(OPTS) $(SRC) $(LIB) -o $(DST) $(INC)
	@echo Compilation complete!

clean: $(DST)
	@$(RM) $(DST)