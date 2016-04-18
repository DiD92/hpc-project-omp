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

#include <assert.h>
#include <omp.h>
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
#define THREAD_NUM 1

#define TRUE 1
#define FALSE 0

#define SIZE_MARGIN 1.10 // % of extra space reservation
#define REALLOC_MARGIN 10
#define INCREASE_FACTOR 100

//--------------------------------------------------------------------------//
// -- AUXILIARY METHODS ----------------------------------------------------//
//--------------------------------------------------------------------------//

double toSeconds(suseconds_t);
long checkForRealloc(void**, long, long, size_t, long);
long rebuildImage(ImageData, DataBucket*);


//--------------------------------------------------------------------------//
// -- LIBRARY IMPLEMENTATION ---------------------------------------------- //
//--------------------------------------------------------------------------//

// Read the corresponding chunk from the source Image
int readChunk(char* filen, intmax_t *offset, intmax_t *limit, DataBucket bucket) {
    intmax_t pos = *offset;
    int value = 0, mult = 10;
    int nvalue = FALSE;
    int increase = INCREASE_FACTOR;
    long k = bucket->offset, bucketms;
    char c;

    FILE *fp;

    if ((fp = openFile(filen, "r")) == NULL) {
        perror("Error: ");
        return -1;
    }

    if(fseek(fp, pos, SEEK_SET)) {
        perror("Error: ");
        return -1;
    }

    for(;pos < *limit;) {
        pos = ftell(fp);
        c = fgetc(fp);
        if(c > 47 && c < 58) {
            value = value * mult + (c - 48);
            nvalue = TRUE;
        } else if(nvalue) {
            bucket->data[k] = value;
            value = 0;
            nvalue = FALSE;
            k++;
            bucketms = bucket->msize;
            bucket->msize = checkForRealloc((void**) &(bucket->data), 
                bucket->msize, k + REALLOC_MARGIN, sizeof(bucket->data[0]),
                increase);
            if(bucketms < bucket->msize) {
                increase *= 2;
            }
            if(bucket->msize == -1) {
                perror("Error: ");
                return -1;
            }
            fflush(stdout);
        }
    }

    bucket->bsize = k;

    fclose(fp);

    return 0;
}
int readImage(ImageData img, FILE **fp, int chunksize, int halosize, 
    long *offset) {
    int haloposition;
    if (fseek(*fp, *offset, SEEK_SET)) {
        perror("Error: ");
    }
    haloposition = chunksize - (img->width * (halosize * 2));
    for(int i = 0; i < chunksize; i++) {
        // When start reading the halo store the position in the image file
        if (halosize != 0 && i == haloposition) {
            *offset = ftell(*fp);
        } 
        fscanf(*fp, "%d %d %d ", &img->R[i], &img->G[i], &img->B[i]);
    }

    return 0;
}

// Duplication of the  just readed source chunk 
// to the destiny image struct chunk
void* duplicateImageChunk(ImageData src, ImageData dst) {

    dst->rsize = checkForRealloc((void**) &(dst->R), dst->rsize, src->rsize, 
        sizeof(dst->R[0]), src->rsize - dst->rsize);

    dst->gsize = checkForRealloc((void**) &(dst->G), dst->gsize, src->gsize, 
        sizeof(dst->G[0]), src->gsize - dst->gsize);

    dst->bsize = checkForRealloc((void**) &(dst->B), dst->bsize, src->bsize, 
        sizeof(dst->B[0]), src->bsize - dst->bsize);

    if(dst->rsize == -1 || dst->bsize == -1 || dst->gsize == -1) {
        return NULL;
    }
    
    if(memcpy((void*) dst->R, (void*) src->R, dst->rsize * sizeof(dst->R[0])) == NULL) {
        return NULL;
    }

    if(memcpy((void*) dst->G, (void*) src->G, dst->gsize * sizeof(dst->G[0])) == NULL) {
        return NULL;
    }

    return memcpy((void*) dst->B, (void*) src->B, dst->bsize * sizeof(dst->B[0]));
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
    if((*fp = openFile(nombre, "w")) == NULL) {
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
int savingChunk(ImageData img, FILE **fp, long *offset, long dataOffst, 
    long count){
    printf("WRITING - %ld - %ld\n", dataOffst, count);
    // Writing image partition
    fseek(*fp, *offset, SEEK_SET);
    for(long i = dataOffst; i < count; i++) {
        fprintf(*fp, "%d\n%d\n%d\n", img->R[i], img->G[i], img->B[i]);
    }
    *offset = ftell(*fp);
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
int convolve2D(int* in, int* out, int dataOffst, int dataSizeX, int dataSizeY,
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

long checkForRealloc(void **ptr, long csize, long margin, size_t psize,
    long reallocInc) {
    long nsize = csize;
    void *temp = NULL, *temp2 = NULL;
    if(nsize < margin) {
        //printf("---NECESSARY REALLOC DETECTED---\n");
        temp2 = *ptr;
        nsize = nsize + reallocInc;
        if((temp = realloc(temp2, nsize * psize)) == NULL) {
            return -1;
        } else {
            *ptr = temp;
        }
        //printf("RESIZED: %ld -> %ld\n", csize, nsize);
        //printf("--------------------------------\n");
    }
    return nsize;
}

long rebuildImage(ImageData img, DataBucket *bucks) {
    long r, g, b, tsize;
    long rasterR, rasterG, rasterB;
    long increaseR, increaseG, increaseB;
    int flip;

    r = g = b = 0L;
    flip = 0;
    increaseR = increaseG = increaseB = INCREASE_FACTOR;

    for(int i = 0; i < 1; i++) {
        for(int j = 0; j < bucks[i]->bsize; j++) {
            switch(flip) {
                case 0:
                    img->R[r] = bucks[i]->data[j];
                    r++;
                    rasterR = img->rsize;
                    img->rsize = checkForRealloc((void**) &(img->R), 
                        img->rsize, (r + REALLOC_MARGIN), sizeof(img->R[0]),
                        increaseR);
                    if(rasterR < img->rsize) {
                        increaseR *= 2;
                    }
                    break;
                case 1:
                    img->G[g] = bucks[i]->data[j];
                    g++;
                    rasterG = img->gsize;
                    img->gsize = checkForRealloc((void**) &(img->G), 
                        img->gsize, (g + REALLOC_MARGIN), sizeof(img->G[0]),
                        increaseG);
                    if(rasterG < img->gsize) {
                        increaseG *= 2;
                    }
                    break;
                case 2:
                    img->B[b] = bucks[i]->data[j];
                    b++;
                    rasterB = img->bsize;
                    img->bsize = checkForRealloc((void**) &(img->B), 
                        img->bsize, (b + REALLOC_MARGIN), sizeof(img->B[0]),
                        increaseB);
                    if(rasterB < img->bsize) {
                        increaseB *= 2;
                    }
                    break;
            }
            flip = (flip + 1) % 3;
        }
        bucks[i]->offset = 0;
    }

    tsize = (r + g + b);

    printf("R: %ld G: %ld B: %ld\n", r, g, b);

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
    int c, offset, t;
    int partitions, partsize, halo, halosize;
    int iwidth, iheight;
    int threadId, threadError;
    long position, chunksize, itersize, bcksize;
    long convoffset, convsize;
    double start, tstart, tend, tread, tcopy, tconv, tstore, treadk;
    struct timeval tim;
    FILE *fpsrc, *fpdst;
    ImageData source, output;
    KernelData kern;
    ImageChunk *chunkLst;
    DataBucket *buckets;

    c = offset = t = 0;
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

    omp_set_dynamic(FALSE);
    omp_set_num_threads(THREAD_NUM);
    
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
    source = parseFileHeader(argv[1], &fpsrc, partitions, halo, SIZE_MARGIN);
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
    if ((output = duplicateImageData(source, partitions, halo, SIZE_MARGIN)) == NULL) {
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

    bcksize = (source->width * source->height * 3) / (partitions * THREAD_NUM);
    bcksize = bcksize + (source->width * halo);
    bcksize = (long) (bcksize * SIZE_MARGIN);

    printf("PSIZE: %ld\n", source->rsize * 3);
    printf("BSIZE: %ld\n", bcksize);
    printf("IRSIZE: %ld\n", source->rsize);

    chunkLst = calculateChunkSections(&fpsrc, source, partitions);

    if ((buckets = initializeBuckets(1, bcksize)) == NULL) {
        perror("Error: ");
        return -1;
    }

    //----------------------------------------------------------------------//
    // CHUNK READING
    //----------------------------------------------------------------------//

    while (c < partitions) {

        printf("-----------------------------------------------\n");
        // Reading Next chunk.
        gettimeofday(&tim, NULL);
        start = tim.tv_sec + toSeconds(tim.tv_usec);

        if (readChunk(argv[1], &(chunkLst[c]->start), &(chunkLst[c]->end), 
            buckets[0])) {
             return -1;
        }

        itersize = rebuildImage(source, buckets);

        // Discarding incomplete row.
        convoffset = (itersize % source->width);
        convsize = itersize - convoffset;

        //Applying offset to bucket
        buckets[0]->offset += (convoffset * 3);
        
        gettimeofday(&tim, NULL);
        tread = tread + (tim.tv_sec + toSeconds(tim.tv_usec) - start);

        // Rows to convolve needs to be bigger than kernel size, either way
        // there'll be problems in pixel alignment.
        if (c == 0) {
            halosize = halo / 2;
            chunksize = (convsize / source->width) - halosize;
            offset = 0; // Value is in rows
            //Applying offset to bucket
            buckets[0]->offset += (source->width * (halo-1) * 3);
        } else if(c < (partitions - 1)) {
            halosize = halo - 1;
            chunksize = (convsize / source->width) - (halo / 2);
            offset = (halo / 2);
            //Applying offset to bucket
            buckets[0]->offset += (source->width * (halo-1) * 3);
        } else {
            halosize = halo / 2;
            chunksize = (convsize / source->width);
            offset = halosize;
        }

        printf("Pixels aligned: %ld\nRows to convolve %ld\n", convsize, 
            chunksize);

        printf("TOTAL BUCKET OFFSET: %d\n", buckets[0]->offset);
        
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
        
        /*convolve2D(source->R, output->R, offset, source->width, chunksize, 
            kern->vkern, kern->kernelX, kern->kernelY);
        convolve2D(source->G, output->G, offset, source->width, chunksize, 
            kern->vkern, kern->kernelX, kern->kernelY);
        convolve2D(source->B, output->B, offset, source->width, chunksize, 
            kern->vkern, kern->kernelX, kern->kernelY);*/
        
        gettimeofday(&tim, NULL);
        tconv = tconv + (tim.tv_sec + toSeconds(tim.tv_usec) - start);
        
        //------------------------------------------------------------------//
        // - CHUNK SAVING --------------------------------------------------//
        //------------------------------------------------------------------//
        gettimeofday(&tim, NULL);

        start = tim.tv_sec + toSeconds(tim.tv_usec);
        if (savingChunk(output, &fpdst, &position, offset * source->width, 
                        chunksize * source->width)) {
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
    printf("|          IMAGE INFO             |\n");
    printf("-----------------------------------\n");
    printf("Name: %s\n", argv[1]);
    printf("Header size (bytes): %ld\n", source->headersize);
    printf("Raster size (bytes): %jd\n", source->rastersize);
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

    freeImagestructure(&source);
    freeImagestructure(&output);

    return 0;
}

//--------------------------------------------------------------------------//
