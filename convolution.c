//--------------------------------------------------------------------------//
//
//  convolution.c
//
//  Created by Josep Lluis Lerida on 11/03/2015
//  Modified by Didac Semente Fernandez on 04/04/2016
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

#include <ctype.h>
#include <stdlib.h>
#include <string.h>
#include <sys/time.h>

//--------------------------------------------------------------------------//

#include "lib/convolution.h"

//--------------------------------------------------------------------------//
// -- MACRO DEFINITION -----------------------------------------------------//
//--------------------------------------------------------------------------//

#define F_MICROS_IN_SECOND 1000000.0
#define THREAD_NUM 1

#define TRUE 1
#define FALSE 0

#define REALLOC_MARGIN 10
#define INCREASE_FACTOR 100

//--------------------------------------------------------------------------//
// -- AUXILIARY METHODS ----------------------------------------------------//
//--------------------------------------------------------------------------//

double calculateExtraSize(int partitions);
double toSeconds(suseconds_t);
long checkForRealloc(void**, long, long, size_t, long);
long rebuildImage(ImageData, DataBucket*);

//--------------------------------------------------------------------------//
// -- LIBRARY IMPLEMENTATION ---------------------------------------------- //
//--------------------------------------------------------------------------//

// Read the corresponding chunk from the source Image
int readChunk(char* fileName, intmax_t *offset, intmax_t *limit, 
    DataBucket bucket) {
    intmax_t pos = *offset;
    int value = 0, mult = 10;
    int newValue = FALSE;
    int increase = INCREASE_FACTOR;
    long k = bucket->offset, bucketMemSize;
    char c;

    FILE *fp;
    int **temp = NULL;

    temp = (int**) malloc(sizeof(int*)); // Avoid breaking strict aliasing

    if((fp = openFile(fileName, "r")) == NULL) {
        perror("Error: ");
        return -1;
    }

    if(fseek(fp, pos, SEEK_SET)) {
        perror("Error: ");
        return -1;
    }

    while(pos < *limit) {
        pos = ftell(fp);
        c = fgetc(fp);
        if(isdigit(c)) {
            value = (value * mult) + (c - '0');
            newValue = TRUE;
        } else if(newValue) {
            bucket->data[k] = value;
            value = 0;
            newValue = FALSE;
            k++;
            // CHECKING IF WE ARE ABOUT TO FILL THE BUCKET
            *temp = bucket->data; 
            bucketMemSize = bucket->msize;
            bucket->msize = checkForRealloc((void**) temp, bucket->msize, 
                k + REALLOC_MARGIN, sizeof(bucket->data[0]), increase);
            bucket->data = *temp;
            if(bucketMemSize < bucket->msize) {
                increase *= 2;
            } else if(bucket->msize == -1) {
                perror("Error: ");
                return -1;
            }
        }
    }

    bucket->bsize = k;

    fclose(fp);
    free(temp);

    return 0;
}

// Duplication of the just readed source chunk 
// to the destiny image struct chunk
void* duplicateImageChunk(ImageData src, ImageData dst) {
    int** temp = NULL;
    temp = (int**) malloc(sizeof(int*)); // Avoid breaking strcit aliasing

    *temp = dst->R;
    dst->rsize = checkForRealloc((void**) temp, dst->rsize, src->rsize, 
        sizeof(dst->R[0]), src->rsize - dst->rsize);
    dst->R = *temp;

    *temp = dst->G;
    dst->gsize = checkForRealloc((void**) temp, dst->gsize, src->gsize, 
        sizeof(dst->G[0]), src->gsize - dst->gsize);
    dst->G = *temp;

    *temp = dst->B;
    dst->bsize = checkForRealloc((void**) temp, dst->bsize, src->bsize, 
        sizeof(dst->B[0]), src->bsize - dst->bsize);
    dst->B = *temp;

    free(temp);

    if(dst->rsize == -1 || dst->bsize == -1 || dst->gsize == -1) {
        return NULL;
    }
    
    if(memcpy((void*) dst->R, (void*) src->R, 
        dst->rsize * sizeof(dst->R[0])) == NULL) {
        return NULL;
    }

    if(memcpy((void*) dst->G, (void*) src->G, 
        dst->gsize * sizeof(dst->G[0])) == NULL) {
        return NULL;
    }

    return memcpy((void*) dst->B, (void*) src->B, 
        dst->bsize * sizeof(dst->B[0]));
}

// Open kernel file and reading kernel matrix. 
// The kernel matrix 2D is stored in 1D format.
KernelData readKernel(char* fileName) {
    FILE *fp;
    int ksize = 0;
    KernelData kern = NULL;
    
    // Opening the kernel file
    if((fp = openFile(fileName, "r")) == NULL) {
        perror("Error: ");
    } else {
        // Memory allocation
        kern = (KernelData) malloc(sizeof(struct structkernel));
        
        // Reading kernel matrix dimensions
        fscanf(fp, "%d,%d,", &kern->kernelX, &kern->kernelY);
        ksize = (kern->kernelX * kern->kernelY);
        kern->vkern = (float*) malloc(ksize * sizeof(float));
        
        // Reading kernel matrix values
        for(int i = 0; i < ksize; i++) {
            fscanf(fp, "%f,", &kern->vkern[i]);
        }

        fclose(fp);
    }

    return kern;
}

// Open the image file with the convolution results
int initfilestore(ImageData img, FILE** fp, char* fileName, long *position) {
    // File with the resulting image is created
    if((*fp = openFile(fileName, "w")) == NULL) {
        perror("Error: ");
        return -1;
    } 
    
    // Writing image header
    fprintf(*fp, "P%d\n%s\n%d %d\n%d\n", img->P, img->comment, img->width,
        img->height, img->maxcolor);
    *position = ftell(*fp);
    return 0;
}

// Writing the image chunk to the resulting file.
int savingChunk(ImageData img, FILE **fp, long *offset, long dataOffst, 
    long count){
    // Writing image partition
    fseek(*fp, *offset, SEEK_SET);
    
    for(long i = dataOffst; i < count; i++) {
        fprintf(*fp, "%d\n%d\n%d\n", img->R[i], img->G[i], img->B[i]);
    }
    *offset = ftell(*fp);
    return 0;
}

// This function frees the space allocated for the image structure.
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
    inPtr = inPtr2 = &in[(dataSizeX * kCenterY) + kCenterX];
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

double calculateExtraSize(int partitions) {
    double x = (double) partitions;
    return (x / (15 + 3*x)) - 0.058;
}

double toSeconds(suseconds_t micros) {
    return (micros / F_MICROS_IN_SECOND);
}

long checkForRealloc(void **ptr, long currSize, long margin, size_t posSize,
    long reallocIncrement) {

    long newSize = currSize;
    void *temp = NULL;

    if(newSize < margin) {
        newSize = newSize + reallocIncrement;
        if((temp = realloc(*ptr, newSize * posSize)) == NULL) {
            free(*ptr);
            return -1;
        } else {
            *ptr = temp;
        }
    }

    return newSize;
}

// Method used to fill the ImageData structure using the data found in the
// DataBucket list.
long rebuildImage(ImageData img, DataBucket *bucks) {
    long r, g, b, tsize;
    long rasterR, rasterG, rasterB;
    long increaseR, increaseG, increaseB;
    int flip, **temp;

    r = g = b = 0L;
    flip = 0;
    increaseR = increaseG = increaseB = INCREASE_FACTOR;
    temp = (int**) malloc(sizeof(int*)); // Avoid breaking strict aliasing

    for(int i = 0; i < 1; i++) {
        for(int j = 0; j < bucks[i]->bsize; j++) {
            switch(flip) {
                case 0:
                    img->R[r] = bucks[i]->data[j];
                    r++;
                    rasterR = img->rsize;
                    *temp = img->R;
                    img->rsize = checkForRealloc((void**) temp, img->rsize, 
                        (r + REALLOC_MARGIN), sizeof(img->R[0]), increaseR);
                    img->R = *temp;
                    if(rasterR < img->rsize) {
                        increaseR *= 2;
                    }
                    break;
                case 1:
                    img->G[g] = bucks[i]->data[j];
                    g++;
                    rasterG = img->gsize;
                    *temp = img->G;
                    img->gsize = checkForRealloc((void**) temp, img->gsize, 
                        (g + REALLOC_MARGIN), sizeof(img->G[0]), increaseG);
                    img->G = *temp;
                    if(rasterG < img->gsize) {
                        increaseG *= 2;
                    }
                    break;
                case 2:
                    img->B[b] = bucks[i]->data[j];
                    b++;
                    rasterB = img->bsize;
                    *temp = img->B;
                    img->bsize = checkForRealloc((void**) temp, img->bsize, 
                        (b + REALLOC_MARGIN), sizeof(img->B[0]), increaseB);
                    img->B = *temp;
                    if(rasterB < img->bsize) {
                        increaseB *= 2;
                    }
                    break;
            }
            flip = (flip + 1) % 3;
        }
        bucks[i]->offset = 0;
    }

    free(temp);

    tsize = (r + g + b);

    // Check for unaligned rasters
    // Either 1 Blue is missing from the image or
    // both 1 Green and 1 Blue.
    switch(tsize % 3) {
        case 0:
            break;
        case 2:
            bucks[0]->offset += 1;
            tsize -= 1;
        case 1:
            bucks[0]->offset += 1;
            tsize -= 1;
            break;
    }

    return (tsize / 3);
}

//--------------------------------------------------------------------------//
// - MAIN METHOD -----------------------------------------------------------//
//--------------------------------------------------------------------------//
int main(int argc, char **argv) {
    int c, offset;
    int partitions, halo, haloSize;
    int imgWidth, imgHeight;
    long position, chunkSize, iterSize, bucketSize;
    long convOffset, convSize;
    double start, tstart, tend, tread, tcopy, tconv, tstore, treadk;
    double extraSizeFactor;
    struct timeval tim;

    char *sourceFile, *outFile, *kernFile;

    FILE *fpsrc, *fpdst;
    ImageData source, output;
    KernelData kern;
    ImageChunk *chunkLst;
    DataBucket *buckets;

    c = offset = 0;
    position = 0L;
    tstart = tend = tread = tcopy = tconv = tstore = treadk = 0.0;
    sourceFile = outFile = kernFile = NULL;
    fpsrc = fpdst = NULL;
    source = output = NULL;
    kern = NULL;

    extraSizeFactor = 1.0;
    
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

    //Storing parameters
    sourceFile = argv[1];
    kernFile = argv[2];
    outFile = argv[3];
    partitions = atoi(argv[4]);
    
    // READING IMAGE HEADERS, KERNEL Matrix, DUPLICATE IMAGE DATA, 
    // OPEN RESULTING IMAGE FILE

    // Reading kernel matrix
    gettimeofday(&tim, NULL);
    start = tim.tv_sec + toSeconds(tim.tv_usec);
    tstart = start;
    if ((kern = readKernel(kernFile)) == NULL) {
        return -1;
    }
    // The matrix kernel defines the halo size to use with the image. 
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

    // Calculating extra size for memory assignment in order to avoid
    // calling realloc further in the execution

    extraSizeFactor = extraSizeFactor + calculateExtraSize(partitions);

    // Memory allocation based on number of partitions and halo size.
    if((source = parseFileHeader(sourceFile, &fpsrc, partitions, 
        halo, extraSizeFactor)) == NULL) {
        return -1;
    }

    imgWidth = source->width;
    imgHeight = source->height;

    gettimeofday(&tim, NULL);
    tread = tread + (tim.tv_sec + toSeconds(tim.tv_usec) - start);
    
    // Duplicate the image struct.
    gettimeofday(&tim, NULL);
    start = tim.tv_sec + toSeconds(tim.tv_usec);
    if ((output = duplicateImageData(source, partitions, halo, 
        extraSizeFactor)) == NULL) {
        return -1;
    }
    gettimeofday(&tim, NULL);
    tcopy = tcopy + (tim.tv_sec + toSeconds(tim.tv_usec) - start);
    
    // Initialize Image output file. Open the file and store the image header
    gettimeofday(&tim, NULL);
    start = tim.tv_sec + toSeconds(tim.tv_usec);
    if (initfilestore(output, &fpdst, outFile, &position) != 0) {
        perror("Error: ");
        return -1;
    }
    gettimeofday(&tim, NULL);
    tstore = tstore + (tim.tv_sec + toSeconds(tim.tv_usec) - start);

    bucketSize = (imgWidth * imgHeight * 3) / (partitions * 1);
    bucketSize = bucketSize + (imgWidth * halo);

    bucketSize = (long)((double) bucketSize * extraSizeFactor);

    chunkLst = calculateChunkSections(&fpsrc, source, partitions);

    if ((buckets = initializeBuckets(1, bucketSize)) == NULL) {
        perror("Error: ");
        return -1;
    }

    //----------------------------------------------------------------------//
    // CHUNK READING
    //----------------------------------------------------------------------//

    while (c < partitions) {

        // Reading chunk.
        gettimeofday(&tim, NULL);
        start = tim.tv_sec + toSeconds(tim.tv_usec);

        if (readChunk(sourceFile, &(chunkLst[c]->start), &(chunkLst[c]->end), 
            buckets[0])) {
             return -1;
        }

        iterSize = rebuildImage(source, buckets);

        gettimeofday(&tim, NULL);
        tread = tread + (tim.tv_sec + toSeconds(tim.tv_usec) - start);

        // Discarding incomplete row.
        convOffset = (iterSize % imgWidth);
        convSize = iterSize - convOffset;

        //Applying offset to bucket
        buckets[0]->offset += (convOffset * 3);

        // Rows to convolve needs to be bigger than kernel size, either way
        // there'll be problems in pixel alignment.

        haloSize = (halo / 2);

        if(c < (partitions - 1)) {
            chunkSize = (convSize / imgWidth) - haloSize;
            buckets[0]->offset += (imgWidth * (halo-1) * 3);
            if(c == 0) {
                offset = 0;
            } else {
                offset = haloSize;
            }

        } else {
            chunkSize = (convSize / imgWidth);
            offset = haloSize;
        }
        
        // Duplicate the image chunk
        gettimeofday(&tim, NULL);
        start = tim.tv_sec + toSeconds(tim.tv_usec);
        if (duplicateImageChunk(source, output) == NULL) {
            perror("Error: ");
            return -1;
        }
        gettimeofday(&tim, NULL);
        tcopy = tcopy + (tim.tv_sec + toSeconds(tim.tv_usec) - start);
        
        //------------------------------------------------------------------//
        // - CHUNK CONVOLUTION ---------------------------------------------//
        //------------------------------------------------------------------//
        gettimeofday(&tim, NULL);
        start = tim.tv_sec + toSeconds(tim.tv_usec);
        
        convolve2D(source->R, output->R, imgWidth, chunkSize, 
            kern->vkern, kern->kernelX, kern->kernelY);
        convolve2D(source->G, output->G, imgWidth, chunkSize, 
            kern->vkern, kern->kernelX, kern->kernelY);
        convolve2D(source->B, output->B, imgWidth, chunkSize, 
            kern->vkern, kern->kernelX, kern->kernelY);
        
        gettimeofday(&tim, NULL);
        tconv = tconv + (tim.tv_sec + toSeconds(tim.tv_usec) - start);
        
        //------------------------------------------------------------------//
        // - CHUNK SAVING --------------------------------------------------//
        //------------------------------------------------------------------//
        gettimeofday(&tim, NULL);

        start = tim.tv_sec + toSeconds(tim.tv_usec);
        if (savingChunk(output, &fpdst, &position, offset * imgWidth, 
                        chunkSize * imgWidth)) {
            perror("Error: ");
            return -1;
        }

        gettimeofday(&tim, NULL);
        tstore = tstore + (tim.tv_sec + toSeconds(tim.tv_usec) - start);

        adjustBucketContents(buckets, 1);

        c++;
    }

    fclose(fpsrc);
    fclose(fpdst);
    
    gettimeofday(&tim, NULL);
    tend = tim.tv_sec + toSeconds(tim.tv_usec);
    
    printf("-----------------------------------\n");
    printf("|          SYSTEM INFO            |\n");
    printf("-----------------------------------\n");
    printf("--------TYPE SIZES (BYTES)---------\n");
    printf("Size of short: ----> %ld\n", sizeof(short));
    printf("Size of int: ------> %ld\n", sizeof(int));
    printf("Size of long: -----> %ld\n", sizeof(long));
    printf("Size of intmax_t: -> %ld\n", sizeof(intmax_t));
    printf("Size of size_t: ---> %ld\n", sizeof(size_t));
    printf("Size of float: ----> %ld\n", sizeof(float));
    printf("Size of double: ---> %ld\n", sizeof(double));
    printf("-----------------------------------\n");
    printf("|          IMAGE INFO             |\n");
    printf("-----------------------------------\n");
    printf("Name: %s\n", sourceFile);
    printf("Header size (bytes): %ld\n", source->headersize);
    printf("Raster size (bytes): %jd\n", source->rastersize);
    printf("ISizeX : %d\n", imgWidth);
    printf("ISizeY : %d\n", imgHeight);
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

    //----------------------------------------------------------------------//
    // - MEMORY CLEANING  --------------------------------------------------//
    //----------------------------------------------------------------------//

    freeImagestructure(&source);
    freeImagestructure(&output);
    freeDataBuckets(buckets, 1);
    freeChunkList(chunkLst, partitions);
    free(kern->vkern);
    free(kern);

    //----------------------------------------------------------------------//

    return 0;
}

//--------------------------------------------------------------------------//
