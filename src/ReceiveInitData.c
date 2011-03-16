// =============================================================================
// 
// Polonator G.007 Image Processing Software
//
// Church Lab, Harvard Medical School
// Written by Greg Porreca
//
// Release 1.0 -- 02-12-2008
//
// This software may be modified and re-distributed, but this header must appear
// at the top of the file.
//
// =============================================================================
//
/* ReceiveInitData.c
-
- Network code adapted from http://cs.baylor.edu/~donahoo/practical/CSockets
-
- Connects to server and receives data blocks; called by
- initialize_processor.c
-
- Written by Greg Porreca (Church Lab) 11-20-2007
- 
*/

#include <stdio.h>
#include <sys/socket.h>
#include <arpa/inet.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <errno.h>
#include <time.h>
#include "ProcessorParams.h"
#include "processor.h"

#define LEN 1000000

#ifdef SAVE_IMAGES
extern char CYCLE_NAME[];
#endif
void importFile(FILE*, char*);
void importImage(short unsigned int *, char*, int);
FILE* createFileList(char * filename); 

void ReceiveInitData(char *argv, int portnum,
            FILE *posfile,
            FILE *infofile,
            int curr_fcnum){

    int sock;                      // socket descriptor
    unsigned int data_size;        // length of data block to receive
    int bytesRcvd, totalBytesRcvd=0; // bytes read in single rcv() and total bytes read
    int imagesRcvd=0;
    char errorString[255];
    char sendBuffer[] = "1"; //CODE OF 1 ASKS SENDER FOR DATA

    short unsigned int *inputBuffer;//
    short unsigned int *image;
    short unsigned int *segmask_image;
    short unsigned int filenum;
    short unsigned int arraynum;
    short unsigned int fcnum;
    short unsigned int img_err_code;
    short unsigned int *beadpos_xcol;
    short unsigned int *beadpos_yrow;
    short unsigned int num_beads;
    int convert;
    short unsigned int sconvert;

    int i, j, k, l, z;

    int num_loops;

    int connection_isok=1;

    char log_string[255];

    #ifdef SAVE_IMAGES
    char outputfilename[255];
    FILE *imgfile;
    #endif

    //sock = GetSock(argv, portnum);

    data_size = (LEN+4) * sizeof(short unsigned int);
    if( (inputBuffer=(short unsigned int *)malloc(data_size)) == NULL){
        p_log_errorno((char*)"malloc() failed");
        exit(0);
    }
    if( (image=(short unsigned int *)malloc(data_size)) == NULL){
        p_log_errorno((char*)"malloc() failed");
        exit(0);
    }
    if( (segmask_image=(short unsigned int *)malloc(data_size)) == NULL){
        p_log_errorno((char*)"malloc() failed");
        exit(0);
    }
    if( (beadpos_xcol=(short unsigned int *)malloc(MAX_BEADS_PERFRAME * sizeof(short unsigned int))) == NULL){
        p_log_errorno((char*)"malloc() failed");
        exit(0);
    }
    if( (beadpos_yrow=(short unsigned int *)malloc(MAX_BEADS_PERFRAME * sizeof(short unsigned int))) == NULL){
        p_log_errorno((char*)"malloc() failed");
        exit(0);
    }

    /*
    //TELL SENDER WE WANT TO RECEIVE DATA
    if(send(sock, sendBuffer, 1, 0)!=1){
    if(errno = EPIPE){
      p_log_errorno((char*)"send() from client to request data failed because connection is broken");
      connection_isok=0;
    }
    else{
      p_log_errorno((char*)"send() from client to request data failed");
      exit(0);
    }
    }


    //CONSUME 'PRIMIMG' DATA SINCE FIRST IMAGE XFER
    //IS VERY SLOW
    if(connection_isok){
    bytesRcvd=0;
    totalBytesRcvd=0;
    while(totalBytesRcvd<data_size){
      if((bytesRcvd=recv(sock, ((char*)inputBuffer)+totalBytesRcvd, data_size-totalBytesRcvd, 0))<0){
    sprintf(errorString, "recv() failed or connection closed prematurely, %d bytes received (%d)", totalBytesRcvd, bytesRcvd);
    close(sock);
    p_log_errorno(errorString);
    exit(0);
      }
      totalBytesRcvd += bytesRcvd;
    }
    }
    */
#ifdef SAVE_IMAGES
    system("mkdir images");
    sprintf(outputfilename, "mkdir images/%s", CYCLE_NAME);
    system(outputfilename);
#endif

    sprintf(log_string, "raw_img/filelist_%s.txt", CYCLE_NAME);
    FILE* file_list = createFileList(log_string);
    char filenamebuffer[255];
    char filenamebuffer2[255];  
    char cwd[1024];

    // loop through all images, receiving blocks of data (1 image + header per block)
    if(1){
        for(i=curr_fcnum; i<curr_fcnum + 1; i++){
            for(j=0; j<ARRAYS_PER_FC; j++){
                for(k=0; k<IMGS_PER_ARRAY; k++){
                    sprintf(log_string, "Receiving data for frame %d %d %d", i, j, k);
                    p_log_simple(log_string);
                    bytesRcvd = 0;
                    totalBytesRcvd = 0;
  

                    importFile(file_list, filenamebuffer);
                    sprintf(filenamebuffer2, "%s/raw_img/0_%s/%s", getcwd(cwd, sizeof(cwd)), CYCLE_NAME,filenamebuffer);
                    importImage(inputBuffer+4, filenamebuffer2, 1000000);

#ifdef WHITE_BEADS
                    // invert for images where beads are lighter than background
                    // for(z=0; z<1000000; z++){
                    //  *(inputBuffer + z + 4) = 16383 - *(inputBuffer + z + 4);
                    //}
#endif
  
                    //flatten image and rescale
                    /*flatten_image(inputBuffer+4, image, 1);*/
                    image = inputBuffer + 4;
  
                    p_log((char*)"Segment image...");
                    //now segment image and get total number of beads;
                    find_objects(image,
                        &num_beads,
                        beadpos_xcol,
                        beadpos_yrow,
                        segmask_image);
  
                    p_log((char*)"Image segmented");
  
  
#ifdef SAVE_IMAGES
  
                    *(inputBuffer+1) = j;
                    *(inputBuffer+2) = k;

                    if(k%100==0)
                    {
                        sprintf(outputfilename, "images/%s/%02d_%04d.raw", CYCLE_NAME, *(inputBuffer+1), *(inputBuffer+2));
                        imgfile = fopen(outputfilename, "w");
                        fwrite(inputBuffer+4, sizeof(short unsigned int), 1000000, imgfile);
                        fwrite(image, sizeof(short unsigned int), 1000000, imgfile);
                        fwrite(segmask_image, sizeof(short unsigned int), 1000000, imgfile);
                        fclose(imgfile);
                    }
#endif
  
                    //now write data to files

                    //header in posfile for current image
                    //posfile is as follows:
                    //--4 bytes each:
                    //-1
                    //curr_fcnum
                    //curr_arraynum
                    //curr_imagenum
                    //num_beads in current image
                    //--2 bytes each:
                    //x for bead 0
                    //y for bead 0
                    //x for bead 1
                    //y for bead 1
                    //....
                    //--4 bytes:
                    //-1
                    convert=-1;
                    p_log((char*)"Write data to object table");
                    if(fwrite(&convert, sizeof(int), 1, posfile) < 1){
                        p_log_errorno((char*)"write of -1 to posfile failed");
                        exit(0);
                    }
                    if(fwrite(&i, sizeof(int), 1, posfile) < 1){
                        p_log_errorno((char*)"write of curr_fcnum to posfile failed");
                        exit(0);
                    }
                    if(fwrite(&j, sizeof(int), 1, posfile) < 1){
                        p_log_errorno((char*)"write of curr_arraynum to posfile failed");
                        exit(0);
                    }
                    if(fwrite(&k, sizeof(int), 1, posfile) < 1){
                        p_log_errorno((char*)"write of curr_imagenum to posfile failed");
                        exit(0);
                    }
                    convert = (int)num_beads;
                    if(fwrite(&convert, sizeof(int), 1, posfile) < 1){
                        p_log_errorno((char*)"write of num_beads to posfile failed");
                        exit(0);
                    }
                    for(l=0; l<num_beads; l++)
                    {
                        //position info to posfile
                        //don't write any values of 0 since this is the separator
                        //between images
                        if(fwrite(beadpos_xcol+l, sizeof(short unsigned int), 1, posfile) < 1){
                            p_log_errorno((char*)"write of beadpos_xcol to posfile failed");
                            exit(0);
                        }
                        if(fwrite(beadpos_yrow+l, sizeof(short unsigned int), 1, posfile) < 1){
                            p_log_errorno((char*)"write of beadpos_yrow to posfile failed");
                            exit(0);
                        }
                    } // end for l
                    sconvert=0;
                    if(fwrite(&sconvert, sizeof(short unsigned int), 1, posfile) < 1){
                        p_log_errorno((char*)"write of 0 to posfile failed");
                        exit(0);
                    }
  
  
                    //then, bead info to infofile;
                    //2 bytes per element, 8 bytes per row (image)
                    sconvert = i;
                    if(fwrite(&sconvert, sizeof(short unsigned int), 1, infofile) < 1){
                        p_log_errorno((char*)"write of curr_fcnum to infofile failed");
                        exit(0);
                    }
                    sconvert = j;
                    if(fwrite(&sconvert, sizeof(short unsigned int), 1, infofile) < 1){
                        p_log_errorno((char*)"write of curr_arraynum to infofile failed");
                        exit(0);
                    }
                    sconvert = k;
                    if(fwrite(&sconvert, sizeof(short unsigned int), 1, infofile) < 1){
                        p_log_errorno((char*)"write of curr_imagenum to infofile failed");
                        exit(0);
                    }
                    if(fwrite(&num_beads, sizeof(short unsigned int), 1, infofile) < 1){
                        p_log_errorno((char*)"write of num_beads to infofile failed");
                        exit(0);
                    }
                    p_log((char*)"Finished handling image");

                }//end for k
            }//end for j
        }//end for i
    }//end if connection_isok

    // tell acq we're finished receiving all data (we don't want it to close the connection
    // before the processor finishes receiving the last image)
    /*
    if(send(sock, sendBuffer, 1, 0)!=1){
    p_log_errorno((char*)"send() from client to signal all data received failed");
    exit(0);
    }


    if(shutdown(sock, 2)<0){
    p_log_errorno((char*)"ERROR on shutdown(sock)");
    exit(0);
    }//close(sock);
    */
} // end function

void ReceiveInitData_new(char *argv, int portnum,
        FILE *posfile,
        FILE *infofile,
        int curr_fcnum)
{

    int sock;                      // socket descriptor
    unsigned int data_size;        // length of data block to receive
    int bytesRcvd, totalBytesRcvd=0; // bytes read in single rcv() and total bytes read
    int imagesRcvd=0;
    char errorString[255];
    char sendBuffer[] = "1"; //CODE OF 1 ASKS SENDER FOR DATA

    short unsigned int *inputBuffer;//
    short unsigned int *image;
    short unsigned int *segmask_image;
    short unsigned int filenum;
    short unsigned int arraynum;
    short unsigned int fcnum;
    short unsigned int img_err_code;
    short unsigned int *beadpos_xcol;
    short unsigned int *beadpos_yrow;
    short unsigned int num_beads;
    int convert;
    short unsigned int sconvert;

    int i, j, k, l, z;

    int chao = 0;

    int num_loops;

    int connection_isok=1;

    char log_string[255];

    #ifdef SAVE_IMAGES
    char outputfilename[255];
    FILE *imgfile;16383,
    #endif

    sock = GetSock(argv, portnum);


    data_size = (LEN+4) * sizeof(short unsigned int);
    if( (inputBuffer=(short unsigned int *)malloc(data_size)) == NULL){
        p_log_errorno((char*)"malloc() failed");
        exit(0);
    }
    if( (image=(short unsigned int *)malloc(data_size)) == NULL){
        p_log_errorno((char*)"malloc() failed");
        exit(0);
    }
    if( (segmask_image=(short unsigned int *)malloc(data_size)) == NULL){
        p_log_errorno((char*)"malloc() failed");
        exit(0);
    }
    if( (beadpos_xcol=(short unsigned int *)malloc(MAX_BEADS_PERFRAME * sizeof(short unsigned int))) == NULL){
        p_log_errorno((char*)"malloc() failed");
        exit(0);
    }
    if( (beadpos_yrow=(short unsigned int *)malloc(MAX_BEADS_PERFRAME * sizeof(short unsigned int))) == NULL){
        p_log_errorno((char*)"malloc() failed");
        exit(0);
    }

    //TELL SENDER WE WANT TO RECEIVE DATA
    if(send(sock, sendBuffer, 1, 0)!=1){
        if(errno = EPIPE){
            p_log_errorno((char*)"send() from client to request data failed because connection is broken");
            connection_isok=0;
        }
        else{
            p_log_errorno((char*)"send() from client to request data failed");
            exit(0);
        }
    }

    /*
    //CONSUME 'PRIMIMG' DATA SINCE FIRST IMAGE XFER
    //IS VERY SLOW
    if(connection_isok){
    bytesRcvd=0;
    totalBytesRcvd=0;
    while(totalBytesRcvd<data_size){
      if((bytesRcvd=recv(sock, ((char*)inputBuffer)+totalBytesRcvd, data_size-totalBytesRcvd, 0))<0){
    sprintf(errorString, "recv() failed or connection closed prematurely, %d bytes received (%d)", totalBytesRcvd, bytesRcvd);
    close(sock);
    p_log_errorno(errorString);
    exit(0);
      }
      totalBytesRcvd += bytesRcvd;
    }
    }
    */

#ifdef SAVE_IMAGES
    system("mkdir autoexp_images");
    sprintf(outputfilename, "mkdir autoexp_images/%s", CYCLE_NAME);
    system(outputfilename);
#endif

    sprintf(outputfilename, "connection %d", connection_isok);
    p_log(outputfilename);

    // loop through all images, receiving blocks of data (1 image + header per block)
    if(connection_isok){
        for(i=curr_fcnum; i<1; i++){
            for(j=0; j<ARRAYS_PER_FC; j++){
                for(k=0; k<1; k++){
                    sprintf(log_string, "Receiving autoexposure images, total = %d;  current = %d",ARRAYS_PER_FC,j+1);
                    p_log(log_string);
                    bytesRcvd = 0;
                    totalBytesRcvd = 0;
  
                    // tell server we're ready for the next block
                    p_log((char*)"Request data");
                    if(send(sock, sendBuffer, 1, 0)!=1)
                    {
                        if(errno = EPIPE){
                            p_log_errorno((char*)"ERROR:\tsend() from client to request data failed because connection is broken");
                            connection_isok=0;
                            totalBytesRcvd=data_size;
                        }
                        else{
                            p_log_errorno((char*)"ERROR:\tsend() from client to request data failed");
                            exit(0);
                        }
                    }
  
                    // wait for full block to be received
                    while(totalBytesRcvd < data_size)
                    {

                        // was there an error during the recv?  if so, this is bad; crash
                        if((bytesRcvd = recv(sock, ((char*)inputBuffer) + totalBytesRcvd, data_size-totalBytesRcvd, MSG_WAITALL)) < 0)
                        {
                            /*if((bytesRcvd = recv(sock, ((char*)inputBuffer) + totalBytesRcvd, data_size-totalBytesRcvd, MSG_WAITALL)) <= 0){*/
                            /* if((bytesRcvd = recv(sock, ((char*)inputBuffer) + totalBytesRcvd, data_size-totalBytesRcvd, 0)) <= 0){*/
                            sprintf(errorString, "recv() failed or connection closed prematurely, %d bytes received", totalBytesRcvd);
                            close(sock);
                            p_log_errorno(errorString);
                            exit(0);
                        }

                        // was the connection broken?  if so, the acq software was probably stopped prematurely; recover gracefully
                        else if(bytesRcvd==0){
                            sprintf(log_string, "connection to acq appears to be broken while trying to receive image %d %d %d, received %d bytes so far", i, j, k, totalBytesRcvd);
                            p_log(log_string);
                            connection_isok=0;
                            totalBytesRcvd=data_size;
                        }

                        sprintf(log_string, "received %d bytes", bytesRcvd);
                        p_log(log_string);
                        totalBytesRcvd += bytesRcvd;
                    }
                    p_log((char*)"autoexposure images received");
                    if(!connection_isok){
                        p_log((char*)"ERROR:\tReceiveData: connection was broken; exiting");
                        break;
                    }

  
                    // NOW THAT WE HAVE THE IMAGE, DO SOMETHING WITH IT
                    //verify it is the image we think it is
                    /* if( (*(inputBuffer + 0) != i) || 
                        (*(inputBuffer + 1) != j) ||
                        (*(inputBuffer + 2) != k) || 
                        (*(inputBuffer + 3) == 1))
                    {
                        sprintf(log_string, "ERROR; expected image %d %d %d 0, got %d %d %d %d",
                            i, j, k,
                            *(inputBuffer),
                            *(inputBuffer+1),
                            *(inputBuffer+2),
                            *(inputBuffer+3));
                            p_log(log_string);
                            exit(0);
                    }*/
  
#ifdef WHITE_BEADS
                    //invert for images where beads are lighter than background
                    for(z=0; z<1000000; z++){
                        *(inputBuffer + z + 4) = 16383 - *(inputBuffer + z + 4);
                    }
#endif
  
                    //flatten image and rescale
                    /*flatten_image(inputBuffer+4, image, 1);*/
                    image = inputBuffer + 4;
  
                    p_log((char*)"Segment image...");
                    //now segment image and get total number of beads;
                    find_objects(image,
                        &num_beads,
                        beadpos_xcol,
                        beadpos_yrow,
                        segmask_image);
  
                    p_log((char*)"Image segmented");
  
  
#ifdef SAVE_IMAGES
                    sprintf(outputfilename, "autoexp_images/%02d_%04d.raw", *(inputBuffer+1), *(inputBuffer+2));
                    p_log(outputfilename);
                    imgfile = fopen(outputfilename, "w");
                    fwrite(inputBuffer+4, sizeof(short unsigned int), 1000000, imgfile);
                    fwrite(image, sizeof(short unsigned int), 1000000, imgfile);
                    fwrite(segmask_image, sizeof(short unsigned int), 1000000, imgfile);
                    fclose(imgfile);
#endif
  
                    //now write data to files

                    //header in posfile for current image
                    //posfile is as follows:
                    //--4 bytes each:
                    //-1
                    //curr_fcnum
                    //curr_arraynum
                    //curr_imagenum
                    //num_beads in current image
                    //--2 bytes each:
                    //x for bead 0
                    //y for bead 0
                    //x for bead 1
                    //y for bead 1
                    //....
                    //--4 bytes:
                    //-1
                    for(chao = 0;chao<15;chao++)
                    {
                        convert=-1;
                        p_log((char*)"Write data to object table");
                        if(fwrite(&convert, sizeof(int), 1, posfile) < 1){
                            p_log_errorno((char*)"write of -1 to posfile failed");
                            exit(0);
                        }
                        if(fwrite(&i, sizeof(int), 1, posfile) < 1){
                            p_log_errorno((char*)"write of curr_fcnum to posfile failed");
                            exit(0);
                        }
                        if(fwrite(&j, sizeof(int), 1, posfile) < 1){
                            p_log_errorno((char*)"write of curr_arraynum to posfile failed");
                            exit(0);
                        }
                        if(fwrite(&chao, sizeof(int), 1, posfile) < 1){
                            p_log_errorno((char*)"write of curr_imagenum to posfile failed");
                            exit(0);
                        }
                        convert = (int)num_beads;
                        if(fwrite(&convert, sizeof(int), 1, posfile) < 1){
                            p_log_errorno((char*)"write of num_beads to posfile failed");
                            exit(0);
                        }
                        for(l=0; l<num_beads; l++){
                            //position info to posfile
                            //don't write any values of 0 since this is the separator
                            //between images
                            if(fwrite(beadpos_xcol+l, sizeof(short unsigned int), 1, posfile) < 1){
                                p_log_errorno((char*)"write of beadpos_xcol to posfile failed");
                                exit(0);
                            }
                            if(fwrite(beadpos_yrow+l, sizeof(short unsigned int), 1, posfile) < 1){
                                p_log_errorno((char*)"write of beadpos_yrow to posfile failed");
                                exit(0);
                            }
                        } // end for l
                        sconvert=0;
                        if(fwrite(&sconvert, sizeof(short unsigned int), 1, posfile) < 1){
                            p_log_errorno((char*)"write of 0 to posfile failed");
                            exit(0);
                        }
      
      
                        //then, bead info to infofile;
                        //2 bytes per element, 8 bytes per row (image)
                        sconvert = i;
                        if(fwrite(&sconvert, sizeof(short unsigned int), 1, infofile) < 1){
                            p_log_errorno((char*)"write of curr_fcnum to infofile failed");
                            exit(0);
                        }
                        sconvert = j;
                        if(fwrite(&sconvert, sizeof(short unsigned int), 1, infofile) < 1){
                            p_log_errorno((char*)"write of curr_arraynum to infofile failed");
                            exit(0);
                        }
                        sconvert = chao;
                        if(fwrite(&sconvert, sizeof(short unsigned int), 1, infofile) < 1){
                            p_log_errorno((char*)"write of curr_imagenum to infofile failed");
                            exit(0);
                        }
                        if(fwrite(&num_beads, sizeof(short unsigned int), 1, infofile) < 1){
                            p_log_errorno((char*)"write of num_beads to infofile failed");
                            exit(0);
                        }
                        p_log((char*)"Finished handling image");
                    } // for chao
                }//end for k
            }//end for j
        }//end for i
    }//end if connection_isok

    // tell acq we're finished receiving all data (we don't want it to close the connection
    // before the processor finishes receiving the last image)
    if(send(sock, sendBuffer, 1, 0)!=1){
        p_log_errorno((char*)"send() from client to signal all data received failed");
        exit(0);
    }

    if(shutdown(sock, 2)<0){
        p_log_errorno((char*)"ERROR on shutdown(sock)");
        exit(0);
    }//close(sock);
} // end functio

void importImage(short unsigned int *buffer, char * filename, int len){
    FILE *thefile;
    unsigned long fileLen = len;

    //Open file
    thefile = fopen(filename, "rb");
    if (!thefile)
    {
        fprintf(stderr, "Unable to open file %s", filename);
        exit (1);
    }
    //fseek (thefile, 4000000, SEEK_SET);
    fread(buffer, sizeof(short unsigned int), fileLen, thefile);
    printf("tests: %s, %d, %d, %d\n", filename, buffer[0], buffer[1000], buffer[500000]);
    fclose(thefile);

}

void importFile(FILE * file_list, char* buffer) {
    char c;
    int i = 0;
    do {
        c = fgetc (file_list);
        buffer[i] = c;
        i++;
    } while (c != '\n');
    buffer[i-1] = 0;
}

FILE* createFileList(char * filename) {
    FILE * file_list;
    file_list = fopen( filename, "r");
    return file_list;
}

//
// ------------------------------------ END ReceiveInitData.c------------------------------------------
//
