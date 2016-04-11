//--------------------------------------------------------------------------//
//
//  ppmparser.h
//
//  Created by Didac Semente Fernandez on 09/04/2016
//
// Implementation of the ppmparser.h library.
//
//--------------------------------------------------------------------------//

//--------------------------------------------------------------------------//
// -- EXTERNAL LIBRARIES -------------------------------------------------- //
//--------------------------------------------------------------------------//

#include <string.h>
#include <stdlib.h>
#include <sys/stat.h>

//--------------------------------------------------------------------------//

#include "lib/ppmparser.h"

//--------------------------------------------------------------------------//
// -- MACRO DEFINITION -----------------------------------------------------//
//--------------------------------------------------------------------------//

#define MIN(X, Y) (((X) < (Y)) ? (X) : (Y))

//--------------------------------------------------------------------------//
// -- LIBRARY IMPLEMENTATION ---------------------------------------------- //
//--------------------------------------------------------------------------//

FILE* openFile(char* filename, char* mode) {
    return fopen(filename, mode);
}

off_t fsize(const char *filename) {
    struct stat st; 

    if (stat(filename, &st) == 0) {
        return st.st_size;
    }

    return -1; 
}

ImageChunk* calculateChunkSections(FILE** f, ImageData img, int partitions) {

    intmax_t partsize = (img->rastersize / (intmax_t) partitions);
    intmax_t sizeleft = img->rastersize;

    ImageChunk *chunkLst;

    chunkLst = (ImageChunk*) malloc(sizeof(ImageChunk) * partitions);

    if(chunkLst == NULL) {
        return NULL;
    }

    for(int i = 0; i < partitions; i++) {
        chunkLst[i] = (ImageChunk) malloc(sizeof(struct imgchunk));
        if(chunkLst[i] == NULL) {
            return NULL;
        }
    }

    chunkLst[0]->start = img->headersize;
    chunkLst[0]->end = getAdjustedPoint(f, img->headersize + partsize);

    sizeleft = sizeleft - (chunkLst[0]->end - chunkLst[0]->start);

    for(int i = 1; i < partitions; i++) {
        partsize = sizeleft / (partitions - i);
        chunkLst[i]->start = chunkLst[i-1]->end;
        chunkLst[i]->end = getAdjustedPoint(f, chunkLst[i]->start + partsize);
        sizeleft = sizeleft - (chunkLst[i]->end - chunkLst[i]->start);
    }

    return chunkLst;

}

DataBucket* initializeBuckets(int nbuckets, long bsize) {

    DataBucket *buckets;

    buckets = (DataBucket*) malloc(sizeof(DataBucket) * nbuckets);

    for(int i = 0; i < nbuckets; i++) {
        buckets[i] = malloc(sizeof(struct databucket));
        buckets[i]->bsize = 0;
        buckets[i]->msize = bsize;
        printf("%ld\n", buckets[i]->msize = bsize);
        buckets[i]->offset = 0;
        buckets[i]->data = calloc(sizeof(int), bsize);
    }

    return buckets;
}

intmax_t getAdjustedPoint(FILE** f, intmax_t next) {
    fseek(*f, next, SEEK_SET);
    char c;
    do {
        c = fgetc(*f);
    } while(c > 47 && c < 58);
    return (ftell(*f)-1);
}

// Open Image file and image struct initialization
ImageData parseFileHeader(char* nombre, FILE **fp, int partitions, int halo,
    double sizeMargin) {
    char c, comment[300];
    int i;
    ImageData img = NULL;
    long chunk = 0L;

    if ((*fp = openFile(nombre, "r")) == NULL) {
        perror("Error: ");
    } else {
        // Memory allocation
        img = (ImageData) malloc(sizeof(struct imageppm));

        // Storing file size
        img->rastersize = (intmax_t) fsize(nombre);
 
        // Reading the first line: Magical Number "P3"
        fscanf(*fp, "%c%d ", &c, &(img->P));
        
        // Reading the image comment
        for(i = 0; (c = fgetc(*fp)) != '\n'; i++) {
            comment[i] = c;
        }
        comment[i] = '\0';
        // Allocating information for the image comment
        img->comment = calloc(strlen(comment), sizeof(char));
        strcpy(img->comment, comment);
        // Reading image dimensions and color resolution
        fscanf(*fp, "%d %d %d\n", &img->width, &img->height, &img->maxcolor);

        img->headersize = ftell(*fp);
        img->rastersize = img->rastersize - img->headersize;

        chunk = img->width * (img->height / partitions);
        // We need to read halo extra rows.
        chunk = chunk + img->width * halo;
        chunk = (long) (chunk * sizeMargin);

        img->rsize = img->gsize = img->bsize = chunk;

        img->R = calloc(img->rsize, sizeof(int));
        img->G = calloc(img->gsize, sizeof(int));
        img->B = calloc(img->bsize, sizeof(int));
        if((img->R == NULL) || (img->G == NULL) || (img->B == NULL) ||
            (img->rastersize == -1)) {
            return NULL;
        }
    }

    return img;
}

// Duplicate the Image struct for the resulting image
ImageData duplicateImageData(ImageData src, int partitions, int halo, 
    double sizeMargin) {
    long chunk;
    // Struct memory allocation
    ImageData dst = (ImageData) malloc(sizeof(struct imageppm));

    // Copying the magic number
    dst->P = src->P;
    // Copying the string comment
    dst->comment = calloc(strlen(src->comment), sizeof(char));
    strcpy(dst->comment, src->comment);
    // Copying image dimensions and color resolution
    dst->width = src->width;
    dst->height = src->height;
    dst->maxcolor = src->maxcolor;
    chunk = dst->width * (dst->height / partitions);
    // We need to read an extra row.
    chunk = chunk + src->width * halo;
    chunk = (long) (chunk * sizeMargin);

    dst->rsize = src->rsize;
    dst->gsize = src->gsize;
    dst->bsize = src->bsize;

    dst->R = calloc(chunk, sizeof(int));
    dst->G = calloc(chunk, sizeof(int));
    dst->B = calloc(chunk, sizeof(int));
    if((dst->R == NULL) || (dst->G == NULL) || (dst->B == NULL)) {
        return NULL;
    }

    return dst;
}
