// =============================================================================
// 
// Polonator G.007 Image Processing Software
//
// Church Lab, Harvard Medical School
// Written by Greg Porreca
//
// Release 1.0 -- 02-12-2008
// Release 1.1 -- 02-02-2009 modified to fix lost 'image send' command bug [GP]
//
// This software may be modified and re-distributed, but this header must appear
// at the top of the file.
//
// =============================================================================
//
/* receive_data.c
-
- Network code adapted from http://cs.baylor.edu/~donahoo/practical/CSockets
-
- Connects to server and receives data blocks; each is a raw image plus a
- header.  Passes the image and header info to ProcessImage for image
- image processing, which consists of cross-correlation registration, data
- extraction, and output to disk.  Called by processor.c
-
- Written by Greg Porreca (Church Lab) 12-14-2007
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
#include "Polonator_logger.h"
#include "processor.h"

#ifdef SAVE_FL_IMAGES
extern char CYCLE_NAME[];
#endif
void importFile(FILE*, char*);
void importImage(short unsigned int *, char*, int);
FILE* createFileList(char * filename); 

void ReceiveData(char *argv, int portnum,
         short unsigned int *reg_pix_list_xcols,
         short unsigned int *reg_pix_list_yrows,
         FILE *full_pixel_list,
         short unsigned int *beadvals,
         FILE *beadfile,
         FILE *sumfile,
         FILE *imgsumfile,
         FILE *reglogfile,
         int curr_fcnum){

    int sock;                      // socket descriptor
    unsigned int data_size;        // length of data block to receive
    int bytesRcvd, totalBytesRcvd=0; // bytes read in single rcv() and total bytes read
    int imagesRcvd=0;
    char errorString[255];
    char log_string[255];
    char sendBuffer[] = "18"; //CODE OF 1/8 ASKS SENDER FOR DATA

    short unsigned int *inputBuffer;//
    short unsigned int img_err_code;

    int i, j, k;

    int *offsetx_history, *offsety_history;

    clock_t time1, time2;
    int num_loops;

    int connection_isok=1;

    #ifdef SAVE_FL_IMAGES
    char outputfilename[255];
    FILE *imgfile;
    #endif
    sprintf(log_string, "raw_img/filelist_%s.txt", CYCLE_NAME);
    FILE* file_list = createFileList(log_string);
    char filenamebuffer[255];
    char filenamebuffer2[255];  
    char cwd[1024];

    //sock = GetSock(argv, portnum);


    data_size = ((NUM_XCOLS * NUM_YROWS)+4) * sizeof(short unsigned int);
    if( (inputBuffer=(short unsigned int *)malloc(data_size)) == NULL){
    p_log_errorno((char*)"malloc() failed");
    }

    //ALLOCATE MEMORY FOR OFFSET HISTORY ARRAYS (USED BY PROCESSIMAGE_REGISTER)
    if((offsetx_history=(int*)malloc(OFFSET_HISTORY_LENGTH * sizeof(int)))==NULL){
    p_log_errorno((char*)"malloc() offsetx_history failed");
    exit(0);
    }
    if((offsety_history=(int*)malloc(OFFSET_HISTORY_LENGTH * sizeof(int)))==NULL){
    p_log_errorno((char*)"malloc() offsety_history failed");
    exit(0);
    }
    for(i=0; i<OFFSET_HISTORY_LENGTH; i++){
    *(offsetx_history + i) = 0;
    *(offsety_history + i) = 0;
    }

    #ifdef SAVE_FL_IMAGES
    system("mkdir images");
    sprintf(outputfilename, "mkdir images/%s", CYCLE_NAME);
    system(outputfilename);
    #endif

    // loop through all images, receiving blocks of data (1 image + header per block)
    if(1){
    for(i=curr_fcnum; i<curr_fcnum + 1; i++){
      for(j=0; j<ARRAYS_PER_FC; j++){
    for(k=0; k<IMGS_PER_ARRAY; k++){
      bytesRcvd = 0;
      totalBytesRcvd = 0;
  
      p_log((char*)"STATUS:\tReceiveData: all data received for current image");
      imagesRcvd++;
      num_loops=0;

      // NOW THAT WE HAVE THE IMAGE, DO SOMETHING WITH IT
    #ifdef SAVE_FL_IMAGES
      *(inputBuffer) = 0;
      *(inputBuffer+1) = j;
      *(inputBuffer+2) = k;
  
      if(k%100==0){
        sprintf(outputfilename, "images/%s/%02d_%04d.raw", CYCLE_NAME, *(inputBuffer+1), *(inputBuffer+2));
        sprintf(log_string, "STATUS:\tReceiveData: write received data to image file %s", outputfilename);
        p_log(log_string);
        imgfile = fopen(outputfilename, "w");
        fwrite(inputBuffer+4, sizeof(short unsigned int), 1000000, imgfile);
        fclose(imgfile);
      }
    #endif
  
    #ifdef DEBUG1

      sprintf(log_string, "STATUS:\tReceiveData: calling ProcessImage for frame %d %d %d, expecting %d %d %d",
          *(inputBuffer),
          *(inputBuffer+1),
          *(inputBuffer+2),
          i, j, k);
      p_log(log_string);
  

      importFile(file_list, filenamebuffer);
      sprintf(filenamebuffer2, "%s/raw_img/%s/%s", getcwd(cwd, sizeof(cwd)), CYCLE_NAME,filenamebuffer);
      importImage(inputBuffer+4, filenamebuffer2, 1000000);
  
    #endif
      ProcessImage(reg_pix_list_xcols,
               reg_pix_list_yrows,
               full_pixel_list,
               inputBuffer + 4,
               inputBuffer,
               beadvals,
               beadfile,
               sumfile,
               imgsumfile,
               reglogfile,
               offsetx_history,
               offsety_history);
    } // end for k
      } // end for j
    } // end for i
    } // end if connection_isok

    free(offsetx_history);
    free(offsety_history);
    free(inputBuffer);

    }

    void ReceiveData_new(char *argv, int portnum,
         short unsigned int *reg_pix_list_xcols,
         short unsigned int *reg_pix_list_yrows,
         FILE *full_pixel_list,
         short unsigned int *beadvals,
         FILE *beadfile,
         FILE *sumfile,
         FILE *imgsumfile,
         FILE *reglogfile,
         FILE *bead_value_file,
         FILE *bead_sum_file,
         int curr_fcnum,
         fpos_t position){

    int sock;                      // socket descriptor
    unsigned int data_size;        // length of data block to receive
    int bytesRcvd, totalBytesRcvd=0; // bytes read in single rcv() and total bytes read
    int imagesRcvd=0;
    char errorString[255];
    char log_string[255];
    char sendBuffer[] = "18"; //CODE OF 1/8 ASKS SENDER FOR DATA
    char sendBuffer_new[] = "b"; //CODE OF 1/8 ASKS SENDER FOR DATA

    short unsigned int *inputBuffer;//
    short unsigned int img_err_code;

    int i, j, k;
    int h = 0;

    int *offsetx_history, *offsety_history;

    unsigned long long int array_value;

    clock_t time1, time2;
    int num_loops;
    unsigned long long int bead_result_value[ARRAYS_PER_FC][15];
    char send_buffer[] = "abcdefghijklmno";
    int connection_isok=1;

    #ifdef SAVE_FL_IMAGES
    char outputfilename[255];
    FILE *imgfile;
    #endif

    //sock = GetSock(argv, portnum);



    data_size = ((NUM_XCOLS * NUM_YROWS)+4) * sizeof(short unsigned int);
    if( (inputBuffer=(short unsigned int *)malloc(data_size)) == NULL){
        p_log_errorno((char*)"malloc() failed");
    }

    //ALLOCATE MEMORY FOR OFFSET HISTORY ARRAYS (USED BY PROCESSIMAGE_REGISTER)
    if((offsetx_history=(int*)malloc(OFFSET_HISTORY_LENGTH * sizeof(int)))==NULL){
        p_log_errorno((char*)"malloc() offsetx_history failed");
        exit(0);
    }
    if((offsety_history=(int*)malloc(OFFSET_HISTORY_LENGTH * sizeof(int)))==NULL){
        p_log_errorno((char*)"malloc() offsety_history failed");
        exit(0);
    }
    for(i=0; i<OFFSET_HISTORY_LENGTH; i++){
        *(offsetx_history + i) = 0;
        *(offsety_history + i) = 0;
    }


    //TELL SENDER WE WANT TO RECEIVE DATA
    if(send(sock, sendBuffer, 1, 0)!=1){
        if(errno = EPIPE){
          p_log_errorno((char*)"ERROR:\tReceivData: send() failed because connection is broken");
          connection_isok=0;
        }
        else{
            p_log_errorno((char*)"send() from client to request data failed");
            exit(0);
        }
    }
    if(send(sock, sendBuffer+1, 1, 0)!=1)
    {
        if(errno = EPIPE){
          p_log_errorno((char*)"ERROR:\tReceivData: send() failed because connection is broken");
          connection_isok=0;
        }
        else{
              p_log_errorno((char*)"send() from client to request data failed");
              exit(0);
        }
    } // end if

#ifdef SAVE_FL_IMAGES
    /* p_log((char*)"trying to create directory");
    h = system("mkdir autoexp_images");
    sprintf(outputfilename, "mkdir autoexp_images/%s", CYCLE_NAME);
    h = system(outputfilename);
    p_log(outputfilename);
    p_log((char*)"the above directory is created");
    if (h == -1)   p_log((char*)"ERROR in creating the above directory");*/
#endif

    // loop through all images, receiving blocks of data (1 image + header per block)
    if(connection_isok)
    {
        for(i=curr_fcnum; i<1; i++){
            for(j=0; j<ARRAYS_PER_FC; j++){
                for(k=0; k<15; k++)
                {
                    bytesRcvd = 0;
                    totalBytesRcvd = 0;
  
                    // tell server we're ready for the next block
    #ifdef DEBUG1
                    p_log((char*)"STATUS:\tReceiveData: request data from server");
                    p_log((char*)"STATUS:\tReceiveData: send first byte to request");
    #endif
                    if(send(sock, sendBuffer, 1, 0)!=1){
                        if(errno = EPIPE){
                            p_log_errorno((char*)"ERROR:\tsend() failed because connection is broken");
                            connection_isok=0;
                            totalBytesRcvd=data_size;
                        }
                        else{
                            p_log_errorno((char*)"ERROR:\tsend() from client to request data failed");
                            exit(0);
                        }
                    }
    #ifdef DEBUG1
                    p_log((char*)"STATUS:\tReceiveData: send second byte to request");
    #endif
                    if(send(sock, sendBuffer+1, 1, 0)!=1){
                        if(errno = EPIPE){
                            p_log_errorno((char*)"ERROR:\tsend() failed because connection is broken");
                            connection_isok=0;
                            totalBytesRcvd=data_size;
                        }
                        else{
                            p_log_errorno((char*)"ERROR:\tsend() from client to request data failed");
                            exit(0);
                        }
                    }

    #ifdef DEBUG1
                    p_log((char*)"STATUS:\tReceiveData: start receiving data");
                    p_log((char*)"STATUS:\tReceiveData: wait for data to arrive");
    #endif
                    // wait for full block to be received
                    while(totalBytesRcvd < data_size)
                    {
                        // was there an error during the recv?  if so, this is bad; crash
                        if((bytesRcvd = recv(sock, ((char*)inputBuffer) + totalBytesRcvd, data_size-totalBytesRcvd, 0)) < 0)
                        {
                            sprintf(errorString, "recv() failed, %d bytes received (%d)", totalBytesRcvd, bytesRcvd);
                            close(sock);
                            p_log_errorno(errorString);
                            // do a graceful restart here
                            connection_isok=0;
                            totalBytesRcvd=data_size;
                            bytesRcvd=0;
                        }
    
                        // was the connection broken?  if so, the acq software was probably stopped prematurely; recover gracefully
                        else if(bytesRcvd==0){
                            sprintf(log_string, "connection to acq appears to be broken while trying to receive image %d %d %d, received %d bytes so far", i, j, k, totalBytesRcvd);
                            p_log(log_string);
                            connection_isok=0;
                            totalBytesRcvd=data_size;
                        }
    /*      
    #ifdef DEBUG1
                        sprintf(log_string, "STATUS:\tReceiveData: received %d bytes of image data", bytesRcvd);
                        p_log(log_string);
    #endif
    */
                        totalBytesRcvd += bytesRcvd;
                        num_loops++;
                    }

                    if(!connection_isok)
                    {
                        i=curr_fcnum+1;
                        j=ARRAYS_PER_FC;    
                        k=IMGS_PER_ARRAY;
                        p_log((char*)"ERROR:\tReceiveData: connection was broken; exiting");
                        break;
                    }

                    p_log((char*)"STATUS:\tReceiveData: all data received for current image");
                    imagesRcvd++;
                    num_loops=0;

                    // NOW THAT WE HAVE THE IMAGE, DO SOMETHING WITH IT
    #ifdef SAVE_FL_IMAGES
                    if(k < 15){
                        sprintf(outputfilename, "autoexp_images/images/%s_%02d_%04d.raw", CYCLE_NAME, j, *(inputBuffer+2));
                        sprintf(log_string, "STATUS:\tReceiveData: write received data to image file %s", outputfilename);
                        p_log(log_string);
                        if((imgfile = fopen(outputfilename, "w")) == NULL){
                                sprintf(log_string, "ERROR opening output FL file %s", outputfilename);
                                p_log_errorno(log_string);
                                exit(0);
                        }
                        fwrite(inputBuffer+4, sizeof(short unsigned int), 1000000, imgfile);
                        fclose(imgfile);
                    }
    #endif
  
    #ifdef DEBUG1
                    sprintf(log_string, "STATUS:\tReceiveData: calling ProcessImage for frame %d %d %d, expecting %d %d %d",
                        *(inputBuffer),
                        *(inputBuffer+1),
                        *(inputBuffer+2),
                        i, j, k);
                    p_log(log_string);
    #endif

                    ProcessImage_new(reg_pix_list_xcols,
                        reg_pix_list_yrows,
                        full_pixel_list,
                        inputBuffer + 4,
                        inputBuffer,
                        beadvals,
                        beadfile,
                        sumfile,
                        imgsumfile,
                        reglogfile,
                        bead_value_file,
                        bead_sum_file,
                        offsetx_history,
                        offsety_history,
                        &array_value,
                        position);
  
                    bead_result_value[j][k] = array_value;
                } // end for k
            } // end for j
        } // end for i
    } // end if connection_isok

    for(j=0;j<ARRAYS_PER_FC;j++)
    {
        for(k=0; k<15; k++)
        {

            //sprintf(log_string, "STATUS:\tReceiveData: averaged FL bead value is %d when using gain %d", bead_result_value[j][k],(k*AUTO_EXPOSURE_STEP_SIZE+AUTO_EXPOSURE_STARTING_GAIN));
            sprintf(log_string, "STATUS:\tReceiveData: averaged FL bead value is %llu when using gain %d", bead_result_value[j][k],(k*AUTO_EXPOSURE_STEP_SIZE+AUTO_EXPOSURE_STARTING_GAIN));
            p_log_simple(log_string);

            if(bead_result_value[j][k]>AUTO_EXPOSURE_THRESHOLD) break;  
        } // end for k
        if (k==15) k=14;

        if(connection_isok)
        {
            if(k!= 14)
            {
                //sprintf(log_string, "STATUS:\tReceiveData: .............for lane %d, using autoexposure gain of %d , as its averaged FL bead value (over top %d \%) = %d, which is greater than the Threshold (%d)............", j,(k*AUTO_EXPOSURE_STEP_SIZE+AUTO_EXPOSURE_STARTING_GAIN),AUTO_EXPOSURE_TOP_PERCENT_VALUE,bead_result_value[j][k],AUTO_EXPOSURE_THRESHOLD);
                sprintf(log_string, "STATUS:\tReceiveData: .............for lane %d, using autoexposure gain of %d , as its averaged FL bead value (over top %d %%) = %llu, which is greater than the Threshold (%d)............", j,(k*AUTO_EXPOSURE_STEP_SIZE+AUTO_EXPOSURE_STARTING_GAIN),AUTO_EXPOSURE_TOP_PERCENT_VALUE,bead_result_value[j][k],AUTO_EXPOSURE_THRESHOLD);
                p_log_simple(log_string);
            }
            else{
                 //sprintf(log_string, "STATUS:\tReceiveData: .............for lane %d, using the maximum autoexposure gain of %d , as averaged FL bead values from previous lower gains (over top %d \%) are less than the Threshold (%d)............", j,(k*AUTO_EXPOSURE_STEP_SIZE+AUTO_EXPOSURE_STARTING_GAIN),AUTO_EXPOSURE_TOP_PERCENT_VALUE,AUTO_EXPOSURE_THRESHOLD);
                 sprintf(log_string, "STATUS:\tReceiveData: .............for lane %d, using the maximum autoexposure gain of %d , as averaged FL bead values from previous lower gains (over top %d %%) are less than the Threshold (%d)............", j,(k*AUTO_EXPOSURE_STEP_SIZE+AUTO_EXPOSURE_STARTING_GAIN),AUTO_EXPOSURE_TOP_PERCENT_VALUE,AUTO_EXPOSURE_THRESHOLD);
                 p_log_simple(log_string);
            }
            fprintf (bead_sum_file, "\n\n.............for lane %d, using autoexposure gain of %d .............", j,(k*AUTO_EXPOSURE_STEP_SIZE+AUTO_EXPOSURE_STARTING_GAIN));

            if(send(sock, &send_buffer[k], 1, 0)!=1){
                p_log_errorno((char*)"send() from client to signal all data received failed");
                //exit(0);
            }
            p_log_simple((char*)"STATUS:\tReceiveData: send gain for the current lane");
            sleep(1);
        } // end if
    } // end for j
    p_log_simple((char*)"STATUS:\t   .... Sending autoexposure gains finished ....   ");

    //sleep(5);
    // tell acq we're finished receiving all data (we don't want it to close the connection
    // before the processor finishes receiving the last image)
    if(connection_isok)
    {
        if(send(sock, sendBuffer, 1, 0)!=1){
          p_log_errorno((char*)"send() from client to signal all data received failed");
          //exit(0);
        }

        if(shutdown(sock, 2)<0){
            p_log_errorno((char*)"ERROR on shutdown(sock)");
            //exit(0);
        }//close(sock);
    } // end if

    free(offsetx_history);
    free(offsety_history);
    free(inputBuffer);

} 

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
// ------------------------------------ END receive_data.c------------------------------------------
//
