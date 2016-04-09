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

//--------------------------------------------------------------------------//

#include "lib/ppmparser.h"

//--------------------------------------------------------------------------//
// -- LIBRARY IMPLEMENTATION ---------------------------------------------- //
//--------------------------------------------------------------------------//

// Open Image file and image struct initialization
ImageData parseFileHeader(char* nombre, FILE **fp, int partitions, int halo) {
    char c, comment[300];
    int i, chunk = 0;
    ImageData img = NULL;

    if ((*fp = fopen(nombre, "r")) == NULL) {
        perror("Error: ");
    } else {
        // Memory allocation
        img = (ImageData) malloc(sizeof(struct imageppm));
 
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
        fscanf(*fp, "%d %d %d", &img->width, &img->height, &img->maxcolor);
        chunk = img->width * (img->height / partitions);
        // We need to read an extra row.
        chunk = chunk + img->width * halo;
        img->R = calloc(chunk, sizeof(int));
        img->G = calloc(chunk, sizeof(int));
        img->B = calloc(chunk, sizeof(int));
        if((img->R == NULL) || (img->G == NULL) || (img->B == NULL)) {
            return NULL;
        }
    }

    return img;
}

// Duplicate the Image struct for the resulting image
ImageData duplicateImageData(ImageData src, int partitions, int halo) {
    int chunk;
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
    dst->R = calloc(chunk, sizeof(int));
    dst->G = calloc(chunk, sizeof(int));
    dst->B = calloc(chunk, sizeof(int));
    if((dst->R == NULL) || (dst->G == NULL) || (dst->B == NULL)) {
        return NULL;
    }

    return dst;
}
