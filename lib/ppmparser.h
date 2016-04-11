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
#include <stdint.h>
#include <sys/types.h>

//--------------------------------------------------------------------------//
// -- TYPE DEFINITIONS ---------------------------------------------------- //
//--------------------------------------------------------------------------//

struct imageppm {
    long headersize;
    intmax_t rastersize;
    int height;
    int width;
    char *comment;
    int maxcolor;
    int P;

    long rsize;
    long bsize;
    long gsize;

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

struct imgchunk {
    intmax_t start;
    intmax_t end;
};

typedef struct imgchunk* ImageChunk;

struct databucket {

    long bsize;
    long msize;
    int offset;
    int *data;
};

typedef struct databucket* DataBucket;

//--------------------------------------------------------------------------//
// -- METHOD DEFINITION --------------------------------------------------- //
//--------------------------------------------------------------------------//

FILE* openFile(char*, char*);
off_t fsize(const char*);
ImageChunk* calculateChunkSections(FILE**, ImageData, int);
DataBucket* initializeBuckets(int, long);
intmax_t getAdjustedPoint(FILE** f, intmax_t next);
ImageData parseFileHeader(char*, FILE**, int, int, double);
ImageData duplicateImageData(ImageData, int, int, double);

//--------------------------------------------------------------------------//
#endif
