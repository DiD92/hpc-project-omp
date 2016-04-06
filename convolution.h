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
// -- TYPE DEFINITIONS ---------------------------------------------------- //
//--------------------------------------------------------------------------//

struct imageppm {
    int height;
    int width;
    char *comment;
    int maxcolor;
    int P;
    int *R;
    int *G;
    int *B;
};

typedef struct imageppm* ImageData;

struct structkernel {
    int kernelX;
    int kernelY;
    float *vkern;
};

typedef struct structkernel* kernelData;

//--------------------------------------------------------------------------//
// -- METHOD DEFINITION --------------------------------------------------- //
//--------------------------------------------------------------------------//

ImageData initimage(char*, FILE**, int, int);
ImageData duplicateImageData(ImageData, int, int);

int readImage(ImageData, FILE**, int, int, long int*);
int duplicateImageChunk(ImageData, ImageData, int);
int initfilestore(ImageData, FILE**, char*, long*);
int savingChunk(ImageData, FILE**, int, int);
int convolve2D(int*, int*, int, int, float*, int, int);
void freeImagestructure(ImageData*);

//--------------------------------------------------------------------------//
