SRC=convolution.c ppmparser.c
LIB=lib
DST=convolution-omp
OPTS=-std=c99 -Wall -Wextra -pedantic-errors -O2
INC=-fopenmp
CC=gcc

conv-omp: $(SRC) $(LIB)
	@$(CC) $(OPTS) $(SRC) -L $(LIB) -o $(DST) $(INC)
	@echo Compilation complete!

clean: $(DST)
	@$(RM) $(DST)