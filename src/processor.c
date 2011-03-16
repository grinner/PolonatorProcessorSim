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
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include "ProcessorParams.h"
#include "processor.h"
#include "Polonator_logger.h"
#include <sys/stat.h>
#ifdef SAVE_FL_IMAGES
char CYCLE_NAME[255];
#endif

int main(int argc, char *argv[]){
    FILE *full_pix_list;
    FILE *regpixfile;
    FILE *beadfile;
    FILE *sumfile;
    FILE *imgsumfile;
    FILE *reglogfile;
    FILE *bead_value_file;
    FILE *bead_sum_file;
    short unsigned int *reg_pix_list_xcols;
    short unsigned int *reg_pix_list_yrows;
    short unsigned int *beadvals;
    int fcnum;
    fpos_t pix_list_offset;
    int reg_X_pointer_offset;
    int reg_Y_pointer_offset;

    int portnum;// = atoi(argv[2]);
    char filename[255];
    char autoefilename[255];
    char autoefilename2[255];
    char rcvfilename[255];
    char posfilename[200];
    char posfilename2[200];
    char reglogfilename[200];
    char temp_string[255];
    char log_string[255];
    char hostname[255];

    int num_regs_skipped;
    int num_fcs;
    int loop_count = 0;
    int i;


    start_logger((char*)DEFAULT_PROCESSOR_LOGFILENAME, 1);

    if(argc !=4 ){
        sprintf(log_string, "Usage: %s <hostname> <Position filename> <Number of flowcells>", argv[0]);
        p_log(log_string);
        exit(0);
    }

    strcpy(posfilename, argv[2]);
    strcpy(posfilename2, posfilename);
    strcpy(hostname,argv[1]);
    portnum = PORTNUM;
    num_fcs = atoi(argv[3]);

    p_log((char*)"processor:\tAllocating memory for arrays ...");
    if( (reg_pix_list_xcols=(short unsigned int*)malloc(REG_PIXELS*IMGS_PER_ARRAY*ARRAYS_PER_FC*num_fcs*sizeof(short unsigned int))) == NULL){
        p_log((char*)"malloc() for reg_pixel_list_xcols failed...");
    }
    if( (reg_pix_list_yrows=(short unsigned int*)malloc(REG_PIXELS*IMGS_PER_ARRAY*ARRAYS_PER_FC*num_fcs*sizeof(short unsigned int))) == NULL){
        p_log((char*)"malloc() for reg_pixel_list_yrows failed...");
    }
    if((beadvals=(short unsigned int*)malloc(MAX_BEADS_PERFRAME*sizeof(short unsigned int))) == NULL){
        p_log((char*)"malloc() for beadvals failed...");
    }


    strcat(posfilename, ".reg");
    if( (regpixfile = fopen(posfilename, "r")) == NULL){
        sprintf(log_string, "ERROR opening reg pixel list file %s", posfilename);
        p_log(log_string);
        exit(0);
    }


    //LOAD REGISTRATION PIXEL LIST INTO MEMORY
    //IF A SET OF COORDS IS WITHIN THE BORDER
    //OF THE IMAGE (SIZE SPECIFIED BY SEARCH_XCOLS
    //AND SEARCH_YROWS, DO NOT USE THAT POINT -- INSTEAD,
    //USE THE PREVIOUS POINT
    num_regs_skipped=0;
    sprintf(log_string, "processor:\tLoading regpix file %s ...", posfilename);
    p_log(log_string);

    for(i=0; i<(REG_PIXELS*IMGS_PER_ARRAY*ARRAYS_PER_FC*num_fcs); i++)
    {
        if((fread(reg_pix_list_xcols+i, 2, 1, regpixfile)) < 1){
              sprintf(log_string, "ERROR reading from file regpix_list at data element %d", i);
              p_log_errorno(log_string);
              exit(0);
        }
        if((fread(reg_pix_list_yrows+i, 2, 1, regpixfile)) < 1){
            sprintf(log_string, "ERROR reading from file regpix_list at data element %d", i);
            p_log_errorno(log_string);
            exit(0);
        }

        if( (*(reg_pix_list_xcols+i) <= (SEARCH_XCOLS+MAX_ADDITIONAL_OFFSET+1)) ||
            (*(reg_pix_list_yrows+i) <= (SEARCH_YROWS+MAX_ADDITIONAL_OFFSET+1)) ||
            ((NUM_XCOLS - *(reg_pix_list_xcols+i)) <= (SEARCH_XCOLS+MAX_ADDITIONAL_OFFSET+1)) ||
            ((NUM_YROWS - *(reg_pix_list_yrows+i)) <= (SEARCH_YROWS+MAX_ADDITIONAL_OFFSET+1)))
        {
            num_regs_skipped++;
            if(i==0)
            {
                *(reg_pix_list_xcols + i) = (int)(NUM_XCOLS/2);
                *(reg_pix_list_yrows + i) = (int)(NUM_YROWS/2);
            }
            else{
                *(reg_pix_list_xcols + i) = *(reg_pix_list_xcols + i - 1);
                *(reg_pix_list_yrows + i) = *(reg_pix_list_yrows + i - 1);
            }
        }
   
    } // end for i
    fclose(regpixfile);
    /*  int percentage_skipped = (int)( ((float)(num_regs_skipped) / (REG_PIXELS*IMGS_PER_ARRAY*ARRAYS_PER_FC*num_fcs)) * 100);*/
    int percentage_skipped = (num_regs_skipped*100.0) / (REG_PIXELS*IMGS_PER_ARRAY*ARRAYS_PER_FC*num_fcs);
    int avg_regpix_perframe = REG_PIXELS - (int)( ((float)(percentage_skipped)/100) * REG_PIXELS);
    sprintf(log_string, "processor:\tSkipped %d reg coords (%d%%), using avg of %d/%d reg pixels at window size %d/%d ...",
        num_regs_skipped,
        percentage_skipped,
        avg_regpix_perframe,
        REG_PIXELS,
        SEARCH_XCOLS,
        SEARCH_YROWS);
    p_log(log_string);
    
    if(AUTO_EXPOSURE){
        p_log((char*)"trying to create the autoe cycle name directory");
        system("mkdir autoexp_images/beads");	
        system("mkdir autoexp_images/logs");
        system("mkdir autoexp_images/images");
    }

    //NOW WAIT FOR DATA TO ARRIVE, AND LOOP; FIRST GET FILENAME, THEN GET DATA
    while(loop_count<1)
    {
        memset(&pix_list_offset, 0, sizeof(pix_list_offset));
        if( (full_pix_list = fopen(posfilename2, "r")) == NULL){
            sprintf(log_string, "ERROR opening pixel list file %s", posfilename2);
            p_log_simple(log_string);
            perror(NULL);
            exit(0);
        }
        if(AUTO_EXPOSURE) run_autoexp(argv[1], atoi(argv[3]));

        p_log_simple((char*)"processor:\tRequesting filename from server ...");
        strcpy(filename, "beads/");
        if(loop_count==0) strcpy(rcvfilename, hostname);
        if(loop_count==1) strcpy(rcvfilename, "AM1a_G");
        if(loop_count==2) strcpy(rcvfilename, "AM1a_C");
        if(loop_count==3) strcpy(rcvfilename, "AM1a_T");
        //ReceiveFilename(argv[1], portnum, rcvfilename);

        sprintf(log_string, "processor:\tReceived filename %s from server ...", rcvfilename);
        p_log(log_string);
        sprintf(temp_string, "logs/%s.processorlog", rcvfilename);
        change_logfile(temp_string);
        p_log(log_string);

        p_log((char*)"processor:\tRequesting flowcell number from server ...");
        fcnum = 0;
        //fcnum = ReceiveFCNum(argv[1], portnum);
        sprintf(temp_string, "%d_", fcnum);
        strcat(filename, temp_string);
        strcat(filename, rcvfilename);
#ifdef SAVE_FL_IMAGES
        strcpy(CYCLE_NAME, temp_string);
        strcpy(CYCLE_NAME, rcvfilename);
#endif
        strcat(filename, ".beads");

        if( (beadfile = fopen(filename, "wb")) == NULL){
            sprintf(log_string, "ERROR opening output bead file %s", filename);
            p_log_errorno(log_string);
            exit(0);
        }
    

        sprintf(temp_string, "logs/%d_", fcnum);
        strcpy(reglogfilename, temp_string);
        strcat(reglogfilename, rcvfilename);
        strcat(reglogfilename, ".register-log");
        if( (reglogfile = fopen(reglogfilename, "wb")) == NULL){
            sprintf(log_string, "ERROR opening output reglog file %s", reglogfilename);
            p_log_errorno(log_string);
            exit(0);
        }

        sprintf(log_string, "processor:\tRequesting data from server for file %s ...", filename);
        p_log(log_string);

        strcat(filename, "ums");
        if( (sumfile = fopen(filename, "wb")) == NULL){
            sprintf(log_string, "ERROR opening output sum file %s", filename);
            p_log_errorno(log_string);
            exit(0);
        }

        strcat(filename, "_full");
        if( (imgsumfile = fopen(filename, "wb")) == NULL){
            sprintf(log_string, "ERROR opening output imgsum file %s", filename);
            p_log_errorno(log_string);
            exit(0);
        }

        // seek to correct place in object table, and advance pointer to
        // correct place in reg object table


        if(fcnum == 0){ // fcnum is either 0 or 1
            // seek to start of object table file
            if(fseek(full_pix_list, 0, SEEK_SET)!=0){
                p_log_errorno((char*)"ERROR ");
            }

        // set reg pointer offsets
            reg_X_pointer_offset = 0;
            reg_Y_pointer_offset = 0;
        }
        else{
            // seek to position in object table where the current flowcell starts
            // we know where to go because the last time through the loop when we 
            // were on FC0, we recorded the position after ReceiveData finished
            fsetpos(full_pix_list, &pix_list_offset);

            // set reg pointer offsets
            reg_X_pointer_offset = 0;
            reg_Y_pointer_offset = 0;
        }


        p_log((char*)"READY TO RECEIVE IMAGE DATA");
        ReceiveData(argv[1], portnum, reg_pix_list_xcols + reg_X_pointer_offset,
            reg_pix_list_yrows + reg_Y_pointer_offset, full_pix_list, beadvals, beadfile, sumfile, imgsumfile, reglogfile, fcnum);

            // remember position in object table so we can get back later
        if(fcnum == 0){
            fgetpos(full_pix_list, &pix_list_offset);
        }
  

        fclose(beadfile);
        fclose(sumfile);
        fclose(imgsumfile);
        fclose(reglogfile);

        p_log((char*)"ReceiveData exited");
        change_logfile((char*)DEFAULT_PROCESSOR_LOGFILENAME);

        p_log((char*)"processor:\tReceiveData exited; waiting before attempting to reconnect...");
        loop_count++;
        //sleep(10);
        fclose(full_pix_list);
    } // end while
} // end function

int run_autoexp(char *argv, int cell_num){
    FILE *full_pix_list;
    FILE *regpixfile;
    FILE *beadfile;
    FILE *sumfile;
    FILE *imgsumfile;
    FILE *reglogfile;
    FILE *bead_value_file;
    FILE *bead_sum_file;
    short unsigned int *reg_pix_list_xcols;
    short unsigned int *reg_pix_list_yrows;
    short unsigned int *beadvals;
    int fcnum =0;
    fpos_t pix_list_offset;
    int reg_X_pointer_offset;
    int reg_Y_pointer_offset;

    int portnum;// = atoi(argv[2]);
    char filename[255];
    char autoefilename[255];
    char autoefilename2[255];
    char rcvfilename[255];
    char posfilename[200];
    char reglogfilename[200];
    char temp_string[255];
    char log_string[255];
    char outputfilename[500];

    int num_regs_skipped;
    int num_fcs;

    int i;
    int h =0;


    start_logger((char*)DEFAULT_PROCESSOR_LOGFILENAME, 1);

    strcpy(posfilename, "autoexp_pos");
    portnum = PORTNUM;
    num_fcs = cell_num;

    p_log((char*)"processor:\tAllocating memory for arrays ...");
    if( (reg_pix_list_xcols=(short unsigned int*)malloc(REG_PIXELS*ARRAYS_PER_FC*15*sizeof(short unsigned int))) == NULL){
        p_log((char*)"malloc() for reg_pixel_list_xcols failed...");
    }
    if( (reg_pix_list_yrows=(short unsigned int*)malloc(REG_PIXELS*ARRAYS_PER_FC*15*sizeof(short unsigned int))) == NULL){
        p_log((char*)"malloc() for reg_pixel_list_yrows failed...");
    }
    if((beadvals=(short unsigned int*)malloc(MAX_BEADS_PERFRAME*sizeof(short unsigned int))) == NULL){
        p_log((char*)"malloc() for beadvals failed...");
    }
    memset(&pix_list_offset, 0, sizeof(pix_list_offset));

    if( (full_pix_list = fopen(posfilename, "r")) == NULL){
        sprintf(log_string, "ERROR opening pixel list file %s", posfilename);
        p_log(log_string);
        perror(NULL);
        return 1;
    }
    strcat(posfilename, ".reg");
    if( (regpixfile = fopen(posfilename, "r")) == NULL){
        sprintf(log_string, "ERROR opening reg pixel list file %s", posfilename);
        p_log(log_string);
        exit(0);
    }


    //LOAD REGISTRATION PIXEL LIST INTO MEMORY
    //IF A SET OF COORDS IS WITHIN THE BORDER
    //OF THE IMAGE (SIZE SPECIFIED BY SEARCH_XCOLS
    //AND SEARCH_YROWS, DO NOT USE THAT POINT -- INSTEAD,
    //USE THE PREVIOUS POINT
    num_regs_skipped=0;
    sprintf(log_string, "processor:\tLoading regpix file %s ...", posfilename);
    p_log(log_string);

    sprintf(log_string, "processor:\tautoe %s ...", posfilename);
    p_log(log_string);

    for(i=0; i<(REG_PIXELS*ARRAYS_PER_FC*15); i++){
        if((fread(reg_pix_list_xcols+i, 2, 1, regpixfile)) < 1){
            sprintf(log_string, "ERROR reading from file regpix_list at data element %d", i);
            p_log_errorno(log_string);
            exit(0);
        }
        if((fread(reg_pix_list_yrows+i, 2, 1, regpixfile)) < 1){
            sprintf(log_string, "ERROR reading from file regpix_list at data element %d", i);
            p_log_errorno(log_string);
            exit(0);
        }

        if( (*(reg_pix_list_xcols+i) <= (SEARCH_XCOLS+MAX_ADDITIONAL_OFFSET+1)) ||
            (*(reg_pix_list_yrows+i) <= (SEARCH_YROWS+MAX_ADDITIONAL_OFFSET+1)) ||
            ((NUM_XCOLS - *(reg_pix_list_xcols+i)) <= (SEARCH_XCOLS+MAX_ADDITIONAL_OFFSET+1)) ||
            ((NUM_YROWS - *(reg_pix_list_yrows+i)) <= (SEARCH_YROWS+MAX_ADDITIONAL_OFFSET+1)))
        {
            num_regs_skipped++;
            if(i==0){
                *(reg_pix_list_xcols + i) = (int)(NUM_XCOLS/2);
                *(reg_pix_list_yrows + i) = (int)(NUM_YROWS/2);
            }
            else{
                *(reg_pix_list_xcols + i) = *(reg_pix_list_xcols + i - 1);
                *(reg_pix_list_yrows + i) = *(reg_pix_list_yrows + i - 1);
            }
        } // end if
   
    } // end for
    fclose(regpixfile);
    /*  int percentage_skipped = (int)( ((float)(num_regs_skipped) / (REG_PIXELS)) * 100);*/
    int percentage_skipped = num_regs_skipped*100.0 / (15*ARRAYS_PER_FC*REG_PIXELS);
    int avg_regpix_perframe = REG_PIXELS - (int)( ((float)(percentage_skipped)/100) * REG_PIXELS);
    sprintf(log_string, "processor:\tSkipped %d reg coords (%d%%), using avg of %d/%d reg pixels at window size %d/%d ...",
        num_regs_skipped,
        percentage_skipped,
        avg_regpix_perframe,
        REG_PIXELS,
        SEARCH_XCOLS,
        SEARCH_YROWS);
        p_log(log_string);

    ReceiveFilename(argv, portnum, rcvfilename);

    /*    p_log((char*)"processor:\tRequesting flowcell number from server ...");
     fcnum = ReceiveFCNum(argv[1], portnum); */
    strcpy(filename, "autoexp_images/beads/");
    sprintf(temp_string, "%d_", fcnum);
    strcat(filename, temp_string);
    strcat(filename, rcvfilename);
#ifdef SAVE_FL_IMAGES
    strcpy(CYCLE_NAME, temp_string);
    strcpy(CYCLE_NAME, rcvfilename);
#endif

    sprintf(log_string, "processor:\tReceived filename %s from server ...", rcvfilename);
    p_log(log_string);

    sprintf(temp_string, "autoexp_images/logs/%s.processorlog", rcvfilename);
    change_logfile(temp_string);
    p_log(log_string);


    strcat(filename, ".beads");

    if( (beadfile = fopen(filename, "wb")) == NULL){
        sprintf(log_string, "ERROR opening output bead file %s", filename);
        p_log_errorno(log_string);
        exit(0);
    }

    strcpy(autoefilename, "autoexp_images/");
    strcat(autoefilename,CYCLE_NAME);
    strcat(autoefilename,"_bead_value_file.txt");

    if( (bead_value_file  = fopen(autoefilename, "w")) == NULL){
        sprintf(log_string, "ERROR opening output bead file %s", filename);
        p_log_errorno(log_string);
        exit(0);
    }

    strcpy(autoefilename2, "autoexp_images/");
    strcat(autoefilename2,CYCLE_NAME);
    strcat(autoefilename2,"_bead_sum_file.txt");

    if( (bead_sum_file = fopen(autoefilename2, "w")) == NULL){
        sprintf(log_string, "ERROR opening output bead file %s", filename);
        p_log_errorno(log_string);
        exit(0);
    } 


    sprintf(temp_string, "autoexp_images/logs/%d_", fcnum);
    strcpy(reglogfilename, temp_string);
    strcat(reglogfilename, rcvfilename);
    strcat(reglogfilename, ".register-log");
    if( (reglogfile = fopen(reglogfilename, "wb")) == NULL){
        sprintf(log_string, "ERROR opening output reglog file %s", reglogfilename);
        p_log_errorno(log_string);
        exit(0);
    }

    sprintf(log_string, "processor:\tRequesting data from server for file %s ...", filename);
    p_log(log_string);

    strcat(filename, "ums");
    if( (sumfile = fopen(filename, "wb")) == NULL){
        sprintf(log_string, "ERROR opening output sum file %s", filename);
        p_log_errorno(log_string);
        exit(0);
    }

    strcat(filename, "_full");
    if( (imgsumfile = fopen(filename, "wb")) == NULL){
        sprintf(log_string, "ERROR opening output imgsum file %s", filename);
        p_log_errorno(log_string);
        exit(0);
    }

    // seek to correct place in object table, and advance pointer to
    // correct place in reg object table
    if(fcnum == 0){ // fcnum is either 0 or 1
        // seek to start of object table file
        if(fseek(full_pix_list, 0, SEEK_SET)!=0){
            p_log_errorno((char*)"ERROR ");
        }

        // set reg pointer offsets
        reg_X_pointer_offset = 0;
        reg_Y_pointer_offset = 0;
    }
    else{
        // seek to position in object table where the current flowcell starts
        // we know where to go because the last time through the loop when we 
        // were on FC0, we recorded the position after ReceiveData finished
        fsetpos(full_pix_list, &pix_list_offset);

        // set reg pointer offsets
        reg_X_pointer_offset = 0;
        reg_Y_pointer_offset = 0;
    }

    p_log((char*)"READY TO RECEIVE IMAGE DATA");
    ReceiveData_new(argv, portnum, reg_pix_list_xcols + reg_X_pointer_offset,
        reg_pix_list_yrows + reg_Y_pointer_offset, full_pix_list, beadvals, beadfile, sumfile, imgsumfile, reglogfile,bead_value_file,bead_sum_file, fcnum, pix_list_offset);

    // remember position in object table so we can get back later
    if(fcnum == 0){
        fgetpos(full_pix_list, &pix_list_offset);
    }
  
    fclose(beadfile);
    fclose(sumfile);
    fclose(imgsumfile);
    fclose(reglogfile);
    fclose(bead_value_file); 
    fclose(bead_sum_file);

    change_logfile((char*)DEFAULT_PROCESSOR_LOGFILENAME);

    p_log_simple((char*)"processor:\tAutoexposure exited; waiting before attempting to reconnect...");
    sleep(10);

    fclose(full_pix_list);

} // end function
