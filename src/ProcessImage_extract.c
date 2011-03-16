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
#include "Polonator_logger.h"
#include "ProcessorParams.h"

void ProcessImage_extract(short unsigned int *image,
                            int offset_xcols,
                            int offset_yrows,
                            short unsigned int img_error,
                            FILE *full_pixel_list,
                            int curr_imagenum,
                            int curr_arraynum,
                            int curr_fcnum,
                            short unsigned int *num_beads,
                            short unsigned int *bead_val,
                            unsigned long long int *beadsum){

    short unsigned int curr_val;
    int curr_val2;
    int f_imagenum, f_arraynum, f_fcnum, f_numobjs; // current vals read from file
    short int curr_xcol, curr_yrow;
    int i;
    char log_string[255];

    if(img_error == 1){
        sprintf(log_string, "WARNING:\tProcessImage_extract: current image (%d %d %d) has error flag set",
            curr_imagenum,
            curr_arraynum,
            curr_fcnum);
        p_log(log_string);
    }

    // SEEK TO CORRECT PLACE IN full_pixel_list FILE
    // I.E. THE FIRST BEAD LOCATION FOR THIS
    // (imagenum, arraynum, fcnum)
    /*  fprintf(stdout, "ProcessImage_extract: seeking to header for correct frame %d %d %d...\n",
      curr_fcnum, curr_arraynum, curr_imagenum);*/
    #ifdef DEBUG1
    p_log((char*)"STATUS:\tProcessImage_extract: seeking to correct place in pixel list file");
    #endif

    if(fread(&curr_val2, 4, 1, full_pixel_list)<1){
        p_log_errorno((char*)"ERROR:\tProcessImage_extract read pixel list file failed(1)");
        exit(1);
    }
    while(curr_val2!=-1){
        sprintf(log_string, "curr_val tested != -1; %d", curr_val2);
        p_log(log_string);
        if(fread(&curr_val, 4, 1, full_pixel_list)<1){
            p_log_errorno((char*)"ERROR:\tProcessImage_extract: read pixel list file failed(2)");
            exit(1);
        }
    }

#ifdef DEBUG1
    p_log((char*)"STATUS:\tProcessImage_extract: read values from pixel list file");
#endif

    // READ IMAGENUM, ARRAYNUM, FCNUM AND COMPARE TO EXPECTED
    // REPORT AN ERROR BUT DO NOT ATTEMPT TO CORRECT
    if(fread(&f_fcnum, 4, 1, full_pixel_list)<1){
        p_log_errorno((char*)"ERROR:\tProcessImage_extract: read pixel list file failed(3)");
        exit(1);
    }
    if(fread(&f_arraynum, 4, 1, full_pixel_list)<1){
        p_log_errorno((char*)"ERROR:\tProcessImage_extract: read pixel list file failed(4)");
        exit(1);
    }
    if(fread(&f_imagenum, 4, 1, full_pixel_list)<1){
        p_log_errorno((char*)"ERROR:\tProcessImage_extract: read pixel list file failed(5)");
        exit(1);
    }
    if(fread(&f_numobjs, 4, 1, full_pixel_list)<1){
        p_log_errorno((char*)"ERROR:\tProcessImage_extract: read pixel list file failed(6)");
        exit(1);
    }

    if((f_imagenum != curr_imagenum)||(f_arraynum != curr_arraynum)||(f_fcnum != curr_fcnum)){
        sprintf(log_string, "ERROR in ProcessImage_extract; expected image %d %d %d does not match pixel list file: %d %d %d",
            curr_fcnum,
            curr_arraynum,
            curr_imagenum,
            f_fcnum,
            f_arraynum,     
            f_imagenum);
        p_log(log_string);
    }

    // THE FLUORESCENCE IMAGE HERE IS MISSING; ADVANCE PIXEL LIST FILE POINTER
    // THIS WONT WORK SINCE POSITIONS ARE 2 BYTES AND FLAG IS 4 BYTES; LAST FLAG 
    // MUST BE 2 BYTES
    if(img_error){
#ifdef DEBUG1
        p_log((char*)"WARNING:\tProcessImage_extract: fluorescence image is missing; advancing pixel list pointer");
#endif
        curr_val = 1;
        while(curr_val!=0)
        {
            if(fread(&curr_val, 2, 1, full_pixel_list)<1){
                p_log_errorno((char*)"ERROR:\tProcessImage_extract: read pixel list file failed(7)");
                exit(1);
            }
        }
    }
    // CURRENT IMAGE IS VALID; EXTRACT PIXEL VALS
    else
    {
#ifdef DEBUG1
        sprintf(log_string, "STATUS:\tProcessImage_extract: extracting %d pixel values", f_numobjs);
        p_log(log_string);
#endif
        *beadsum = 0;
        for(i=0; i<f_numobjs; i++){
            if(fread(&curr_xcol, 2, 1, full_pixel_list)<1){
                p_log_errorno((char*)"ERROR:\tProcessImage_extract: read pixel list file failed(8)");
                exit(1);
            }
            if(fread(&curr_yrow, 2, 1, full_pixel_list)<1){
                p_log_errorno((char*)"ERROR:\tProcessImage_extract: read pixel list file failed(9)");
                exit(1);
            }
  
            // DETERMINE OFFSET X AND Y VALS AND VERIFY THEY ARE WITHIN RANGE
            // IF NOT, SET PIXEL VALUES TO 0 AS A FLAG (NO BEAD CAN HAVE VALUE OF 0)
            if(((curr_xcol + offset_xcols) < 0) ||
                ((curr_yrow + offset_yrows) < 0) ||
                ((curr_xcol + offset_xcols) >= NUM_XCOLS) ||
                ((curr_yrow + offset_yrows) >= NUM_YROWS))
            {
                *(bead_val + i) = 0;
            }
            else{
        // ASSIGN CURRENT VAL IN BEAD LIST; IF EXTRACTED VALUE == 0, SET IT
        // TO 1 SINCE 0 IS NOT IN THE VALID PIXEL DATA RANGE OF THE BASECALLER
                *(bead_val + i) = *(image + (NUM_XCOLS * (curr_yrow + offset_yrows)) + (curr_xcol + offset_xcols));
                *beadsum+=*(bead_val+i);
                if(*(bead_val + i) == 0){
                    *(bead_val + i) = 1;
                    *(beadsum)++;
                }
            }
        } // end for i
    } // end else
    *num_beads = f_numobjs;



#ifdef DEBUG1
    p_log((char*)"STATUS:\tProcessImage_extract: advancing pixel list file pointer");
#endif
    if(fread(&curr_val, 2, 1, full_pixel_list)<1){
        p_log_errorno((char*)"ERROR:\tProcessImage_extract: read pixel list file failed(10)");
        exit(1);
    }
    while(curr_val!=0){
        sprintf(log_string, "ERROR: curr_val tested != 0; %d", curr_val);
        p_log(log_string);
        if(fread(&curr_val, 2, 1, full_pixel_list)<1){
            p_log_errorno((char*)"ERROR:\tProcessImage_extract: read pixel list file failed(11)");
            exit(1);
        }
    } // end while
} // end function

int compare (const void * a, const void * b)
{
  return ( *(short unsigned int*)a - *(short unsigned int*)b );
}

void ProcessImage_extract_new(short unsigned int *image,
                                int offset_xcols,
                                int offset_yrows,
                                short unsigned int img_error,
                                FILE *full_pixel_list,
                                int curr_imagenum,
                                int curr_arraynum,
                                int curr_fcnum,
                                short unsigned int *num_beads,
                                short unsigned int *bead_val,
                                unsigned long long int *beadsum,
                                FILE *bead_value_file,
                                  FILE *bead_sum_file,
                                short unsigned int *array,
                                unsigned long long int *array_index,
                                fpos_t position)
{

    short unsigned int curr_val;
    int curr_val2;
    int f_imagenum, f_arraynum, f_fcnum, f_numobjs; // current vals read from file
    short int curr_xcol, curr_yrow;
    int i;
    char log_string[255];
      fpos_t pix_list_offset;
    unsigned long long int bead1 = 0;
    unsigned long long int bead2 = 0;
    unsigned long long int bead3 = 0;


    if(img_error == 1){
        sprintf(log_string, "WARNING:\tProcessImage_extract: current image (%d %d %d) has error flag set",
            curr_imagenum,
            curr_arraynum,
            curr_fcnum);
        p_log(log_string);
    }

    // SEEK TO CORRECT PLACE IN full_pixel_list FILE
    // I.E. THE FIRST BEAD LOCATION FOR THIS
    // (imagenum, arraynum, fcnum)
    /*  fprintf(stdout, "ProcessImage_extract: seeking to header for correct frame %d %d %d...\n",
      curr_fcnum, curr_arraynum, curr_imagenum);*/
#ifdef DEBUG1
    p_log((char*)"STATUS:\tProcessImage_extract: seeking to correct place in pixel list file");
#endif
    /*  
    if((curr_imagenum == 0)&&(curr_arraynum>0)){
          p_log((char*)"STATUS:\tProcessImage_extract: update position in pixel list file");
          fgetpos(full_pixel_list, &position);
    }
    else if(curr_arraynum>0){
          p_log((char*)"STATUS:\tProcessImage_extract: set current position in pixel list file");
          fsetpos(full_pixel_list, &position);
    }
    else if(curr_arraynum==0){
      p_log((char*)"set position to 0 in the file pos file");
      if(fseek(full_pixel_list, 0, SEEK_SET)!=0){
    p_log_errorno((char*)"ERROR");
      }
    }
    */
    if(fread(&curr_val2, 4, 1, full_pixel_list)<1){
        p_log_errorno((char*)"ERROR:\tProcessImage_extract read pixel list file failed(1)");
        exit(1);
    }
    while(curr_val2!=-1){
        sprintf(log_string, "curr_val tested != -1; %d", curr_val2);
        p_log(log_string);
        if(fread(&curr_val, 4, 1, full_pixel_list)<1){
            p_log_errorno((char*)"ERROR:\tProcessImage_extract: read pixel list file failed(2)");
            exit(1);
        }
    }

#ifdef DEBUG1
    p_log((char*)"STATUS:\tProcessImage_extract: read values from pixel list file");
#endif

    // READ IMAGENUM, ARRAYNUM, FCNUM AND COMPARE TO EXPECTED
    // REPORT AN ERROR BUT DO NOT ATTEMPT TO CORRECT
    if(fread(&f_fcnum, 4, 1, full_pixel_list)<1){
        p_log_errorno((char*)"ERROR:\tProcessImage_extract: read pixel list file failed(3)");
        exit(1);
    }
    if(fread(&f_arraynum, 4, 1, full_pixel_list)<1){
        p_log_errorno((char*)"ERROR:\tProcessImage_extract: read pixel list file failed(4)");
        exit(1);
    }
    if(fread(&f_imagenum, 4, 1, full_pixel_list)<1){
        p_log_errorno((char*)"ERROR:\tProcessImage_extract: read pixel list file failed(5)");
        exit(1);
    }
    if(fread(&f_numobjs, 4, 1, full_pixel_list)<1){
        p_log_errorno((char*)"ERROR:\tProcessImage_extract: read pixel list file failed(6)");
        exit(1);
    }

    array = (short unsigned int*) malloc(f_numobjs*sizeof(short unsigned int));

    if((f_imagenum != curr_imagenum)||(f_arraynum != curr_arraynum)||(f_fcnum != curr_fcnum)){
        sprintf(log_string, "ERROR in ProcessImage_extract; expected image %d %d %d does not match pixel list file: %d %d %d",
            curr_fcnum,
            curr_arraynum,
            curr_imagenum,
            f_fcnum,
            f_arraynum,     
            f_imagenum);
        p_log(log_string);
    }

    // THE FLUORESCENCE IMAGE HERE IS MISSING; ADVANCE PIXEL LIST FILE POINTER
    // THIS WONT WORK SINCE POSITIONS ARE 2 BYTES AND FLAG IS 4 BYTES; LAST FLAG 
    // MUST BE 2 BYTES
    if(img_error){
#ifdef DEBUG1
        p_log((char*)"WARNING:\tProcessImage_extract: fluorescence image is missing; advancing pixel list pointer");
#endif
        curr_val = 1;
        while(curr_val!=0){
            if(fread(&curr_val, 2, 1, full_pixel_list)<1){
                p_log_errorno((char*)"ERROR:\tProcessImage_extract: read pixel list file failed(7)");
                exit(1);
            }
        }
    }
    // CURRENT IMAGE IS VALID; EXTRACT PIXEL VALS
    else{
#ifdef DEBUG1
        sprintf(log_string, "STATUS:\tProcessImage_extract: extracting %d pixel values", f_numobjs);
        p_log(log_string);
#endif
        *beadsum = 0;
        for(i=0; i<f_numobjs; i++)
        {
            if(fread(&curr_xcol, 2, 1, full_pixel_list)<1){
                p_log_errorno((char*)"ERROR:\tProcessImage_extract: read pixel list file failed(8)");
                exit(1);
            }
            if(fread(&curr_yrow, 2, 1, full_pixel_list)<1){
                p_log_errorno((char*)"ERROR:\tProcessImage_extract: read pixel list file failed(9)");
                exit(1);
            }
  
            // DETERMINE OFFSET X AND Y VALS AND VERIFY THEY ARE WITHIN RANGE
            // IF NOT, SET PIXEL VALUES TO 0 AS A FLAG (NO BEAD CAN HAVE VALUE OF 0)
            if(((curr_xcol + offset_xcols) < 0) ||
                ((curr_yrow + offset_yrows) < 0) ||
                ((curr_xcol + offset_xcols) >= NUM_XCOLS) ||
                ((curr_yrow + offset_yrows) >= NUM_YROWS))
            {
                *(bead_val + i) = 0;
                array[i] = 0;
            }
            else{
            // ASSIGN CURRENT VAL IN BEAD LIST; IF EXTRACTED VALUE == 0, SET IT
            // TO 1 SINCE 0 IS NOT IN THE VALID PIXEL DATA RANGE OF THE BASECALLER
                *(bead_val + i) = *(image + (NUM_XCOLS * (curr_yrow + offset_yrows)) + (curr_xcol + offset_xcols));
                array[i] = *(bead_val + i);
                *beadsum+=*(bead_val+i);
                if(*(bead_val + i) == 0){
                    *(bead_val + i) = 1;
                    array[i] = *(bead_val + i);
                    *(beadsum)++;
                }
            }
        /*      sprintf(log_string, "STATUS:\tProcessImage_extract: current bead is %d, current bead value is %d, array index value = %d, array size is %f", i, *(bead_val + i), array[i], sizeof(array)/sizeof(short unsigned int));
            p_log(log_string);*/
        } // end for
    } // end else
    *num_beads = f_numobjs;

    qsort(array, f_numobjs, sizeof(short unsigned int), compare);

    for(i=0; i<f_numobjs; i++){

        bead1 += array[f_numobjs-1-i];
        bead3 += array[i];
        /*      sprintf(log_string, "STATUS:\tProcessImage: index value = %d, count from the back is %d, bead1 = %d", array[i], array[f_numobjs-1-i],bead1);
        p_log(log_string);
        */
        fprintf (bead_value_file, "bead index = %d, bead value = %d\n",i,array[i]);          
        if(f_numobjs>99)
        {
            if((i == (int)(f_numobjs*1.0/100.0)) ||
                (i == (int)(f_numobjs*2.0/100.0)) ||
                (i == (int)(f_numobjs*5.0/100.0)) ||
                (i == (int)(f_numobjs*10.0/100.0)) ||
                (i == (int)(f_numobjs*20.0/100.0)) ||
                (i == (int)(f_numobjs*30.0/100.0)) ||
                (i == (int)(f_numobjs*40.0/100.0))|| 
                (i == (int)(f_numobjs*50.0/100.0)) ||
                (i == (int)(f_numobjs*60.0/100.0)) ||
                (i == (int)(f_numobjs*70.0/100.0)) ||
                (i == (int)(f_numobjs*80.0/100.0)) ||
                (i == (int)(f_numobjs*90.0/100.0)) ||
                (i == (f_numobjs-1)))
            {

                bead2 = (int)(i*1.0/(f_numobjs*1.0/100.0))+1;

                //fprintf (bead_sum_file, "top %d percent (average of top %d bead values) is %d\n", bead2, i, (bead1/i));
                fprintf (bead_sum_file, "top %llu percent (average of top %d bead values) is %llu\n", bead2, i, (bead1/i));

                //fprintf (bead_sum_file, "bottom %d percent (average of bottom %d bead values) is %d\n", bead2, i, (bead3/i));
                fprintf (bead_sum_file, "bottom %llu percent (average of bottom %d bead values) is %llu\n", bead2, i, (bead3/i));

                //fprintf (bead_sum_file, "Signal to Noise Ratio (for top/bottom %d percent) = %f\n\n", bead2, (bead1*1.0/(bead3*1.0)));
                fprintf (bead_sum_file, "Signal to Noise Ratio (for top/bottom %llu percent) = %f\n\n", bead2, (bead1*1.0/(bead3*1.0)));

    /*          sprintf(log_string, "STATUS:\tProcessImage_extract: top %d percent (average of top %d bead values) is %d\n", bead2, i, (bead1/i));
                p_log_simple(log_string);   

                sprintf(log_string, "STATUS:\tProcessImage_extract: bottom %d percent (average of bottom %d bead values) is %d\n", bead2, i, (bead3/i));
                p_log_simple(log_string);   

                sprintf(log_string, "STATUS:\tProcessImage_extract: Signal to Noise Ratio (for top/bottom %d percent) = %f\n", bead2, (bead1*1.0/(bead3*1.0)));
                p_log_simple(log_string);   */

                if(i == (int)(f_numobjs*AUTO_EXPOSURE_TOP_PERCENT_VALUE*1.0/100.0))
                {
                    // sprintf(log_string, "STATUS:\tProcessImage_extract: use average of top %d percent, which is %d\n", bead2, (bead1/i));
                    sprintf(log_string, "STATUS:\tProcessImage_extract: use average of top %llu percent, which is %llu\n", bead2, (bead1/i));
                    // fprintf (bead_sum_file, "....... use average of top %d percent bead values to set the auto_exposure gain, which is %d .......\n\n", bead2, (bead1/i));
                    fprintf (bead_sum_file, "....... use average of top %llu percent bead values to set the auto_exposure gain, which is %llu .......\n\n", bead2, (bead1/i));
                    p_log_simple(log_string);
                    *array_index = (bead1/i);
                } // end if
            } // end if
        } // end if
    }

    if(f_numobjs==0){
        sprintf(log_string, " .... WARNING!!! No bead found in the lane-alignment current image, autoexposure will use the minium gain %d .... \n", f_numobjs);
        p_log_errorno(log_string);
        *array_index = 0;
    } 
    else if((f_numobjs<100) && (f_numobjs>0)){
        sprintf(log_string, " .... WARNING!!! less than 100 bead found in the lane-alignment current image, using the entire beads population %d .... \n", f_numobjs);
        p_log_simple(log_string);
        *array_index = bead3/f_numobjs;
    }

#ifdef DEBUG1
    p_log((char*)"STATUS:\tProcessImage_extract: advancing pixel list file pointer");
#endif
    if(fread(&curr_val, 2, 1, full_pixel_list)<1){
        p_log_errorno((char*)"ERROR:\tProcessImage_extract: read pixel list file failed(10)");
        exit(1);
    }
    while(curr_val!=0){
        sprintf(log_string, "ERROR: curr_val tested != 0; %d", curr_val);
        p_log(log_string);
        if(fread(&curr_val, 2, 1, full_pixel_list)<1){
            p_log_errorno((char*)"ERROR:\tProcessImage_extract: read pixel list file failed(11)");
            exit(1);
        }   
    }
} // end function

