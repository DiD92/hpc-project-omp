//--------------------------------------------------------------------------//
//
//  convolution.c
//
//  Created by Josep Lluis Lerida on 11/03/2015
// Modified by Didac Semente Fernandez on 04/04/2016
//
// This program calculates the convolution for PPM images.
// The program accepts an PPM image file, a text definition of the kernel 
// matrix and the PPM file for storing the convolution results.
// The program allows to define image partitions for processing larger 
// images (>500MB).
// The 2D image is represented by 1D vector for chanel R, G and B. 
// The convolution is applied to each chanel separately.
//
//--------------------------------------------------------------------------//

//--------------------------------------------------------------------------//
// -- EXTERNAL LIBRARIES -------------------------------------------------- //
//--------------------------------------------------------------------------//

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/time.h>

//--------------------------------------------------------------------------//

#include "lib/convolution.h"

//--------------------------------------------------------------------------//
// -- MACRO DEFINITION -----------------------------------------------------//
//--------------------------------------------------------------------------//

#define F_MICROS_IN_SECOND 1000000.0

//--------------------------------------------------------------------------//
// -- AUXILIARY METHODS ----------------------------------------------------//
//--------------------------------------------------------------------------//

double toSeconds(suseconds_t);

//--------------------------------------------------------------------------//
// -- LIBRARY IMPLEMENTATION ---------------------------------------------- //
//--------------------------------------------------------------------------//

// Read the corresponding chunk from the source Image
int readImage(ImageData img, FILE **fp, int dim, int halosize, 
              long *position) {
    int haloposition;
    if (fseek(*fp, *position, SEEK_SET)) {
        perror("Error: ");
    }
    haloposition = dim - (img->width * (halosize * 2));
    for(int i = 0; i < dim; i++) {
        // When start reading the halo store the position in the image file
        if (halosize != 0 && i == haloposition) {
            *position=ftell(*fp);
        } 
        fscanf(*fp, "%d %d %d ", &img->R[i], &img->G[i], &img->B[i]);
    }

    return 0;
}

// Duplication of the  just readed source chunk 
// to the destiny image struct chunk
int duplicateImageChunk(ImageData src, ImageData dst) {
    
    memcpy((void*) dst, (void*) src, sizeof(dst));

    return 0;
}

// Open kernel file and reading kernel matrix. 
// The kernel matrix 2D is stored in 1D format.
KernelData readKernel(char* nombre) {
    FILE *fp;
    int i = 0, ksize = 0;
    KernelData kern = NULL;
    
    // Opening the kernel file
    fp = fopen(nombre, "r");
    if(!fp) {
        perror("Error: ");
    } else {
        // Memory allocation
        kern = (KernelData) malloc(sizeof(struct structkernel));
        
        // Reading kernel matrix dimensions
        fscanf(fp, "%d,%d,", &kern->kernelX, &kern->kernelY);
        ksize = (kern->kernelX * kern->kernelY);
        kern->vkern = (float*) malloc(ksize * sizeof(float));
        
        // Reading kernel matrix values
        for(i = 0; i < ksize; i++) {
            fscanf(fp, "%f,", &kern->vkern[i]);
        }

        fclose(fp);
    }

    return kern;
}

// Open the image file with the convolution results
int initfilestore(ImageData img, FILE** fp, char* nombre, long *position) {
    // File with the resulting image is created
    if ((*fp=fopen(nombre,"w")) == NULL) {
        perror("Error: ");
        return -1;
    }
    
    // Writing image header
    fprintf(*fp, "P%d\n%s\n%d %d\n%d\n", img->P, img->comment, img->width,
        img->height, img->maxcolor);
    *position = ftell(*fp);
    return 0;
}

// Writing the image partition to the resulting file. dim is 
// the exact size to write. offset is the displacement for avoid halos.
int savingChunk(ImageData img, FILE **fp, int dim, int offset){
    // Writing image partition
    int len = dim + offset;
    for(int i = offset; i < len; i++) {
        fprintf(*fp, "%d %d %d ", img->R[i], img->G[i], img->B[i]);
    }
    return 0;
}

// This function free the space allocated for the image structure.
void freeImagestructure(ImageData *src) {
    
    free((*src)->comment);
    free((*src)->R);
    free((*src)->G);
    free((*src)->B);
    free(*src);
}

//--------------------------------------------------------------------------//
// 2D convolution
// 2D data are usually stored in computer memory as contiguous 1D array.
// So, we are using 1D array for 2D data.
// 2D convolution assumes the kernel is center originated, which means, if
// kernel size 3 then, k[-1], k[0], k[1]. The middle of index is always 0.
// The following programming logics are somewhat complicated because of using
// pointer indexing in order to minimize the number of multiplications.
//
//
// signed integer (32bit) version:
//--------------------------------------------------------------------------//
int convolve2D(int* in, int* out, int dataSizeX, int dataSizeY,
               float* kernel, int kernelSizeX, int kernelSizeY) {
    int m, n;
    int *inPtr, *inPtr2, *outPtr;
    float *kPtr;
    int kCenterX, kCenterY;
    int rowMin, rowMax;                    // to check boundary of input array
    int colMin, colMax;                    //
    float sum;                             // temp accumulation buffer
    
    // Parameter validatin
    if(!in || !out || !kernel || dataSizeX <= 0 || kernelSizeX <= 0) { 
        return -1;
    }
    
    // Find centeral position of kernel (half of kernel size)
    kCenterX = (int) kernelSizeX / 2;
    kCenterY = (int) kernelSizeY / 2;
    
    // init working  pointers
    // note that  it is shifted (kCenterX, kCenterY),
    inPtr = inPtr2 = &in[dataSizeX * kCenterY + kCenterX];
    outPtr = out;
    kPtr = kernel;
    
    // start convolution
    // number of rows
    for(int i = 0; i < dataSizeY; ++i) {
        // compute the range of convolution, the current row of kernel 
        // should be between these
        rowMax = i + kCenterY;
        rowMin = i - dataSizeY + kCenterY;
        
        // number of columns
        for(int j = 0; j < dataSizeX; ++j) {
            // compute the range of convolution, the current column of kernel 
            // should be between these
            colMax = j + kCenterX;
            colMin = j - dataSizeX + kCenterX;
            
            sum = 0;                        // set to 0 before accumulate
            
            // flip the kernel and traverse all the kernel values
            // multiply each kernel value with underlying input data
            // kernel rows
            for(m = 0; m < kernelSizeY; ++m) {
                // check if the index is out of bound of input array
                if(m <= rowMax && m > rowMin) {
                    for(n = 0; n < kernelSizeX; ++n) {
                        // check the boundary of array
                        if(n <= colMax && n > colMin) {
                            sum += *(inPtr - n) * (*kPtr);
                        }
                        
                        ++kPtr;// next kernel
                    }
                } else {
                    // out of bound, move to next row of kernel
                    kPtr += kernelSizeX; 
                }
                
                // move input data 1 raw up
                inPtr -= dataSizeX;   
            }
            
            // convert integer number
            if(sum >= 0) { 
                *outPtr = (int)(sum + 0.5f);
            } else  { // For using with image editors like GIMP or others...
                *outPtr = (int)(sum - 0.5f);
            }
            
            kPtr = kernel; // reset kernel to (0,0)
            inPtr = ++inPtr2; // next input
            ++outPtr; // next output
        }
    }
    
    return 0;
}

//--------------------------------------------------------------------------//

double toSeconds(suseconds_t micros) {
    return (micros / F_MICROS_IN_SECOND);
}

//--------------------------------------------------------------------------//
// - MAIN METHOD -----------------------------------------------------------//
//--------------------------------------------------------------------------//
int main(int argc, char **argv) {
    int c, offset;
    int partitions, partsize, chunksize, halo, halosize;
    int iwidth, iheight;
    long position;
    double start, tstart, tend, tread, tcopy, tconv, tstore, treadk;
    struct timeval tim;
    FILE *fpsrc, *fpdst;
    ImageData source, output;
    KernelData kern;

    c = offset = 0;
    position = 0L;
    tstart = tend = tread = tcopy = tconv = tstore = treadk = 0.0;
    fpsrc = fpdst = NULL;
    source = output = NULL;
    kern = NULL;
    
    if(argc != 5) {
        printf("Usage: %s <image-file> <kernel-file> <result-file> "
            "<partitions>\n", argv[0]);
        
        printf("\n\nError, Missing parameters:\n");
        printf("format: ./serialconvolution image_file kernel_file "
            "result_file\n");
        printf("- image_file : source image path (*.ppm)\n");
        printf("- kernel_file: kernel path (text file with 1D "
            "kernel matrix)\n");
        printf("- result_file: result image path (*.ppm)\n");
        printf("- partitions : Image partitions\n\n");
        return -1;
    }
    
    // READING IMAGE HEADERS, KERNEL Matrix, DUPLICATE IMAGE DATA, 
    // OPEN RESULTING IMAGE FILE
    
    // Store number of partitions
    partitions = atoi(argv[4]);

    // Reading kernel matrix
    gettimeofday(&tim, NULL);
    start = tim.tv_sec + toSeconds(tim.tv_usec);
    tstart = start;
    if ( (kern = readKernel(argv[2])) == NULL) {
        return -1;
    }
    // The matrix kernel define the halo size to use with the image. 
    // The halo is zero when the image is not partitioned.
    if (partitions == 1) {
        halo = 0;
    } else { 
        halo = kern->kernelY;
    }
    gettimeofday(&tim, NULL);
    treadk = treadk + (tim.tv_sec + toSeconds(tim.tv_usec) - start);

    // Reading Image Header. Image properties: Magical number, comment, 
    // size and color resolution.
    gettimeofday(&tim, NULL);
    start = tim.tv_sec + toSeconds(tim.tv_usec);
    // Memory allocation based on number of partitions and halo size.
    source = parseFileHeader(argv[1], &fpsrc, partitions, halo);
    if (source == NULL) {
        return -1;
    }

    iwidth = source->width;
    iheight = source->height;

    gettimeofday(&tim, NULL);
    tread = tread + (tim.tv_sec + toSeconds(tim.tv_usec) - start);
    
    // Duplicate the image struct.
    gettimeofday(&tim, NULL);
    start = tim.tv_sec + toSeconds(tim.tv_usec);
    if ((output = duplicateImageData(source, partitions, halo)) == NULL) {
        return -1;
    }
    gettimeofday(&tim, NULL);
    tcopy = tcopy + (tim.tv_sec + toSeconds(tim.tv_usec) - start);
    
    // Initialize Image Storing file. Open the file and store the image header
    gettimeofday(&tim, NULL);
    start = tim.tv_sec + toSeconds(tim.tv_usec);
    if (initfilestore(output, &fpdst, argv[3], &position) != 0) {
        perror("Error: ");
        return -1;
    }
    gettimeofday(&tim, NULL);
    tstore = tstore + (tim.tv_sec + toSeconds(tim.tv_usec) - start);

    //----------------------------------------------------------------------//
    // CHUNK READING
    //----------------------------------------------------------------------//
    partsize  = (source->height * source->width) / partitions;

    while (c < partitions) {
        // Reading Next chunk.
        gettimeofday(&tim, NULL);
        start = tim.tv_sec + toSeconds(tim.tv_usec);
        if (c == 0) {
            halosize  = halo / 2;
            chunksize = partsize + (source->width * halosize);
            offset   = 0;
        } else if(c < (partitions - 1)) {
            halosize  = halo;
            chunksize = partsize + (source->width*halosize);
            offset    = (source->width * (halo / 2));
        } else {
            halosize  = halo / 2;
            chunksize = partsize + (source->width*halosize);
            offset    = (source->width * halosize);
        }
        
        if (readImage(source, &fpsrc, chunksize, halo / 2, &position)) {
            return -1;
        }
        gettimeofday(&tim, NULL);
        tread = tread + (tim.tv_sec + toSeconds(tim.tv_usec) - start);
        
        // Duplicate the image chunk
        gettimeofday(&tim, NULL);
        start = tim.tv_sec + toSeconds(tim.tv_usec);
        if ( duplicateImageChunk(source, output) ) {
            return -1;
        }
        gettimeofday(&tim, NULL);
        tcopy = tcopy + (tim.tv_sec+toSeconds(tim.tv_usec) - start);
        
        //------------------------------------------------------------------//
        // - CHUNK CONVOLUTION ---------------------------------------------//
        //------------------------------------------------------------------//
        gettimeofday(&tim, NULL);
        start = tim.tv_sec + toSeconds(tim.tv_usec);
        
        convolve2D(source->R, output->R, source->width, 
            (source->height/partitions) + halosize, kern->vkern, 
            kern->kernelX, kern->kernelY);
        convolve2D(source->G, output->G, source->width, 
            (source->height/partitions) + halosize, kern->vkern, 
            kern->kernelX, kern->kernelY);
        convolve2D(source->B, output->B, source->width, 
            (source->height/partitions) + halosize, kern->vkern, 
            kern->kernelX, kern->kernelY);
        
        gettimeofday(&tim, NULL);
        tconv = tconv + (tim.tv_sec + toSeconds(tim.tv_usec) - start);
        
        //------------------------------------------------------------------//
        // - CHUNK SAVING --------------------------------------------------//
        //------------------------------------------------------------------//
        gettimeofday(&tim, NULL);
        start = tim.tv_sec + toSeconds(tim.tv_usec);
        if (savingChunk(output, &fpdst, partsize, offset)) {
            perror("Error: ");
            return -1;
        }
        gettimeofday(&tim, NULL);
        tstore = tstore + (tim.tv_sec+toSeconds(tim.tv_usec) - start);
        c++;
    }

    fclose(fpsrc);
    fclose(fpdst);
    
    freeImagestructure(&source);
    freeImagestructure(&output);
    
    gettimeofday(&tim, NULL);
    tend = tim.tv_sec + toSeconds(tim.tv_usec);
    
    printf("-----------------------------------\n");
    printf("|          IMAGE INFO             |\n");
    printf("-----------------------------------\n");
    printf("Name: %s\n", argv[1]);
    printf("ISizeX : %d\n", iwidth);
    printf("ISizeY : %d\n", iheight);
    printf("kSizeX : %d\n", kern->kernelX);
    printf("kSizeY : %d\n", kern->kernelY);
    printf("-----------------------------------\n");
    printf("|         EXECUTION TIMES         |\n");
    printf("-----------------------------------\n");
    printf("%.6lfs elapsed in reading image file.\n", tread);
    printf("%.6lfs elapsed in copying image structure.\n", tcopy);
    printf("%.6lfs elapsed in reading kernel matrix.\n", treadk);
    printf("%.6lfs elapsed computing the convolution.\n", tconv);
    printf("%.6lfs elapsed in writing the resulting image.\n", tstore);
    printf("-----------------------------------\n");
    printf("%.6lfs elapsed in total.\n", tend-tstart);
    printf("-----------------------------------\n");
    return 0;
}

//--------------------------------------------------------------------------//
