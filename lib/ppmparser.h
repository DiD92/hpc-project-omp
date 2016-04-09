#ifndef __PPM_PARSER_H__
#define __PPM_PARSER_H__
//--------------------------------------------------------------------------//
//
//  ppmparser.h
//
//  Created by Didac Semente Fernandez on 09/04/2016
//
// Parser used to read PPM files compliant with the PPM format standard,
// standard source: http://netpbm.sourceforge.net/doc/ppm.html
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

typedef struct structkernel* KernelData;

//--------------------------------------------------------------------------//
// -- METHOD DEFINITION --------------------------------------------------- //
//--------------------------------------------------------------------------//

ImageData parseFileHeader(char*, FILE**, int, int);
ImageData duplicateImageData(ImageData, int, int);

//--------------------------------------------------------------------------//
#endif
