SRC=convolution.c ppmparser.c
LIB=lib
DST=convolution-omp
OPTS=-std=c99 -Wall -Wextra -pedantic-errors -O2 -g
INC=-fopenmp
CC=gcc

conv-omp: $(SRC) $(LIB)
	@$(CC) $(OPTS) $(SRC) -L $(LIB) -o $(DST) $(INC)
	@echo Compilation complete!

clean: $(DST)
	@$(RM) $(DST)

prof: $(DST)
	valgrind --tool=massif --time-unit=ms \
	--massif-out-file=mout/massif.out.%p \
	./convolution-omp ../images/NB.ppm ../kernels/kern3x3.txt \
	../result-images/NB-2.ppm 2
