#ifndef __CONVOLUTION_H__
#define __CONVOLUTION_H__
//--------------------------------------------------------------------------//
//
//  convolution.h
//
//  Created by Josep Lluis Lerida on 11/03/2015
//  Modified by Didac Semente Fernandez on 04/04/2016
//
// Type and method definitions for the convolution.c file.
//
//--------------------------------------------------------------------------//

//--------------------------------------------------------------------------//
// -- EXTERNAL LIBRARIES -------------------------------------------------- //
//--------------------------------------------------------------------------//

#include <stdio.h>

//--------------------------------------------------------------------------//

#include "ppmparser.h"

//--------------------------------------------------------------------------//
// -- METHOD DEFINITION --------------------------------------------------- //
//--------------------------------------------------------------------------//

int readImage(ImageData, FILE**, int, int, long int*);
int duplicateImageChunk(ImageData, ImageData);
KernelData readKernel(char*);
int initfilestore(ImageData, FILE**, char*, long*);
int savingChunk(ImageData, FILE**, int, int);
int convolve2D(int*, int*, int, int, float*, int, int);
void freeImagestructure(ImageData*);

//--------------------------------------------------------------------------//
#endif