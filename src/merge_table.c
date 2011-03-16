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
#include <math.h>
#include <iostream>
#include <unistd.h>
#include "ProcessorParams.h"
#include "processor.h"
#include "Polonator_logger.h"

int main(int argc, char *argv[]){

  FILE *posfile;
  FILE *cy3_posfile;
  FILE *fam_posfile;
  FILE *txred_posfile;
  FILE *final_posfile;
  FILE *obj_logfile;
 
  FILE *infofile;

  fpos_t offsetval;

  char regfilename[200];
  char logfilename[200];
  char log_string[500];
  char info_filename[200];

  int i, j, k, l,m,n,o,convert,sconvert, object_count;
  int curr_fc, curr_array, curr_img;
  int num_fcs;
  int cy5_proceed, cy3_proceed, txred_proceed,fam_proceed;
  long long int found_object = 0;
  long long int total_object = 0;

  int curr_numobjs,cy3_curr_numobjs,txred_curr_numobjs,fam_curr_numobjs, max_curr_numobjs = 0;

  int **objects_array;
  int p,q;
  short unsigned int cy5_x,cy3_x,cy3_y,txred_x,txred_y,fam_x,fam_y,cy5_y;

  short unsigned int *curr_img_x;
  short unsigned int *curr_img_y;
  short unsigned int *regindices;
  short unsigned int temp;

  start_logger((char*)"initialize_processor-log", 1);

  
  system("rm 00_*.info");
  system("rm autoexp_pos.info");

  objects_array = (int**) malloc(sizeof(int) * 4000);
  for (i = 0; i < 2000; i++)
  {
    objects_array[i] = (int*) malloc(sizeof(int) * 2000);
  }

	for(p=0;p<2000;p++){
	   for(q=0;q<2000;q++){
		objects_array[p][q] = -1;
		}
	}

  if((curr_img_x = (short unsigned int*) malloc(MAX_BEADS_PERFRAME * sizeof(short unsigned int))) == NULL){
    p_log_errorno((char*)"ERROR allocating memory for curr_img_x");
    exit(0);
  }
  if((curr_img_y = (short unsigned int*) malloc(MAX_BEADS_PERFRAME * sizeof(short unsigned int))) == NULL){
    p_log_errorno((char*)"ERROR allocating memory for curr_img_y");
    exit(0);
  }
  
  if( (obj_logfile=fopen("merge.log", "w")) == NULL){
    sprintf(log_string, "ERROR opening reglog file %s for output", argv[2]);
    p_log_errorno(log_string);
    exit(0);
  }


  if((cy3_posfile=fopen("00_cy3", "rb")) == NULL){
    sprintf(log_string, "ERROR opening position file %s for input", "cy3");
    p_log_errorno(log_string);
    exit(0);
  }

  if( (txred_posfile=fopen("00_txr", "rb")) == NULL){
    sprintf(log_string, "ERROR opening position file %s for input", "txred");
    p_log_errorno(log_string);
    exit(0);
  }

  if( (fam_posfile=fopen("00_fam", "rb")) == NULL){
    sprintf(log_string, "ERROR opening position file %s for input", "fam");
    p_log_errorno(log_string);
    exit(0);
  }

  if( (posfile=fopen("00_cy5", "rb")) == NULL){
    sprintf(log_string, "ERROR opening position file %s for input", "cy5");
    p_log_errorno(log_string);
    exit(0);
  }

  if( (final_posfile=fopen(argv[1], "wb")) == NULL){
    sprintf(log_string, "ERROR opening position file %s for input", argv[1]);
    p_log_errorno(log_string);
    exit(0);
  }
  
  strcpy(info_filename,argv[1]);
  strcat(info_filename,".info");
  if( (infofile=fopen(info_filename, "wb")) == NULL){
    sprintf(log_string, "ERROR opening reg file %s for output", info_filename);
    p_log_errorno(log_string);
    exit(0);
  }


  // loop through FC, ARRAY, IMAGE:
  for(i=0; i<1; i++){
    for(j=0; j<ARRAYS_PER_FC; j++){
      for(k=0; k<IMGS_PER_ARRAY; k++){
	// load positions for current image
	// ------
	// each 'record' in the pos file is as follows:
	// 4 bytes:  -1
	// 4 bytes:  flowcell # ([0..num_fcs])
	// 4 bytes:  array #    ([0..ARRAYS_PER_FC])
  	// 4 bytes:  image #    ([0..IMGS_PER_ARRAY])
	// 4 bytes:  # of beads
	// 2 bytes:  beadpos_xcol
	// 2 bytes:  beadpos_yrow
	// ...
	// 2 bytes:  0
	//
 	// ------

	if( (fread(&curr_fc, sizeof(int), 1, posfile)) < 1){
	  p_log_errorno((char*)"ERROR reading curr_fc(1) from posfile");
	}
	if( (fread(&curr_fc, sizeof(int), 1, posfile)) < 1){
	  p_log_errorno((char*)"ERROR reading curr_fc(2) from posfile");
	}
	if( (fread(&curr_array, sizeof(int), 1, posfile)) < 1){
	  p_log_errorno((char*)"ERROR reading curr_array from posfile");
	}
	if( (fread(&curr_img, sizeof(int), 1, posfile)) < 1){
	  p_log_errorno((char*)"ERROR reading curr_img from posfile");
	}

	if( (fread(&curr_numobjs, sizeof(int), 1, posfile)) < 1){
	  p_log_errorno((char*)"ERROR reading curr_numobjs from posfile");
	}

/*	sprintf(log_string, "flowcell %d, curr_array %d, curr_img %d",curr_fc, curr_array, curr_img);
	p_log(log_string);
*/
	if( (fread(&curr_fc, sizeof(int), 1, cy3_posfile)) < 1){
	  p_log_errorno((char*)"ERROR reading curr_fc(1) from posfile");
	}
	if( (fread(&curr_fc, sizeof(int), 1, cy3_posfile)) < 1){
	  p_log_errorno((char*)"ERROR reading curr_fc(2) from posfile");
	}
	if( (fread(&curr_array, sizeof(int), 1, cy3_posfile)) < 1){
	  p_log_errorno((char*)"ERROR reading curr_array from posfile");
	}
	if( (fread(&curr_img, sizeof(int), 1, cy3_posfile)) < 1){
	  p_log_errorno((char*)"ERROR reading curr_img from posfile");
	}

	if( (fread(&curr_fc, sizeof(int), 1, txred_posfile)) < 1){
	  p_log_errorno((char*)"ERROR reading curr_fc(1) from posfile");
	}
	if( (fread(&curr_fc, sizeof(int), 1, txred_posfile)) < 1){
	  p_log_errorno((char*)"ERROR reading curr_fc(2) from posfile");
	}
	if( (fread(&curr_array, sizeof(int), 1, txred_posfile)) < 1){
	  p_log_errorno((char*)"ERROR reading curr_array from posfile");
	}
	if( (fread(&curr_img, sizeof(int), 1, txred_posfile)) < 1){
	  p_log_errorno((char*)"ERROR reading curr_img from posfile");
	}

	if( (fread(&curr_fc, sizeof(int), 1, fam_posfile)) < 1){
	  p_log_errorno((char*)"ERROR reading curr_fc(1) from posfile");
	}
	if( (fread(&curr_fc, sizeof(int), 1, fam_posfile)) < 1){
	  p_log_errorno((char*)"ERROR reading curr_fc(2) from posfile");
	}
	if( (fread(&curr_array, sizeof(int), 1, fam_posfile)) < 1){
	  p_log_errorno((char*)"ERROR reading curr_array from posfile");
	}
	if( (fread(&curr_img, sizeof(int), 1, fam_posfile)) < 1){
	  p_log_errorno((char*)"ERROR reading curr_img from posfile");
	}

	max_curr_numobjs = 0;
	if(max_curr_numobjs<curr_numobjs) max_curr_numobjs=curr_numobjs; 

	if( (fread(&cy3_curr_numobjs, sizeof(int), 1, cy3_posfile)) < 1){
	  p_log_errorno((char*)"ERROR reading curr_numobjs from posfile");
	}
	if(max_curr_numobjs<cy3_curr_numobjs) max_curr_numobjs=cy3_curr_numobjs;
	if( (fread(&txred_curr_numobjs, sizeof(int), 1, txred_posfile)) < 1){
	  p_log_errorno((char*)"ERROR reading curr_numobjs from posfile");
	}
	if(max_curr_numobjs<txred_curr_numobjs) max_curr_numobjs=txred_curr_numobjs;
	if( (fread(&fam_curr_numobjs, sizeof(int), 1, fam_posfile)) < 1){
	  p_log_errorno((char*)"ERROR reading curr_numobjs from posfile");
	}
	if(max_curr_numobjs<fam_curr_numobjs) max_curr_numobjs=fam_curr_numobjs;

        l = 0;
	m = 0;
        n = 0;
	o = 0;
        p = 0;
	q = 0;
	cy5_proceed = 1;
	cy3_proceed = 1;
	txred_proceed = 1;
	fam_proceed = 1;
	found_object = 0;
	total_object = 0;

	for(p=0;p<2000;p++){
	   for(q=0;q<2000;q++){
		objects_array[p][q] = -1;
		}
	}

	while(l<max_curr_numobjs)
	{

		if(l<curr_numobjs)
		{

			  if((fread(&cy5_x, sizeof(short unsigned int), 1, posfile)) < 1){
			    sprintf(log_string, "ERROR reading x from posfile at %d %d %d %d",
				    i, j, k, l);
			    p_log_errorno(log_string);
			  }
			  if((fread(&cy5_y, sizeof(short unsigned int), 1, posfile)) < 1){
			    sprintf(log_string, "ERROR reading y from posfile at %d %d %d %d",
				    i, j, k, l);
			    p_log_errorno(log_string);
			  }

			  if ((objects_array[(int)cy5_x][(int)cy5_y] == -1)&&((int)cy5_y!=1000 && (int)cy5_x!=1000)){
				objects_array[(int)cy5_x][(int)cy5_y] = 1;
				total_object++;
			  }

		}
		if(l<cy3_curr_numobjs)
		{
			  if((fread(&cy3_x, sizeof(short unsigned int), 1, cy3_posfile)) < 1){
			    sprintf(log_string, "ERROR reading x from posfile at %d %d %d %d",
				    i, j, k, l);
			    p_log_errorno(log_string);
			  }

			  if((fread(&cy3_y, sizeof(short unsigned int), 1, cy3_posfile)) < 1){
			    sprintf(log_string, "ERROR reading y from posfile at %d %d %d %d",
				    i, j, k, l);
			    p_log_errorno(log_string);
			  }

			  if ((objects_array[(int)cy3_x][(int)cy3_y] == -1)&&((int)cy3_y!=1000 && (int)cy3_x!=1000)){
				objects_array[(int)cy3_x][(int)cy3_y] = 1;
				total_object++;
			  }
		}

		if(l<txred_curr_numobjs)
		{
			  if((fread(&txred_x, sizeof(short unsigned int), 1, txred_posfile)) < 1){
			    sprintf(log_string, "ERROR reading x from posfile at %d %d %d %d",
				    i, j, k, l);
			    p_log_errorno(log_string);
			  }
			  if((fread(&txred_y, sizeof(short unsigned int), 1, txred_posfile)) < 1){
			    sprintf(log_string, "ERROR reading y from posfile at %d %d %d %d",
				    i, j, k, l);
			    p_log_errorno(log_string);
			  }
			  if ((objects_array[(int)txred_x][(int)txred_y] == -1)&&((int)txred_y!=1000 && (int)txred_x!=1000)){
				objects_array[(int)txred_x][(int)txred_y] = 1;
				total_object++;
			  }
		}

		if(l<fam_curr_numobjs)
		{
			  if((fread(&fam_x, sizeof(short unsigned int), 1, fam_posfile)) < 1){
			    sprintf(log_string, "ERROR reading x from posfile at %d %d %d %d",
				    i, j, k, l);
			    p_log_errorno(log_string);
			  }
			  if((fread(&fam_y, sizeof(short unsigned int), 1, fam_posfile)) < 1){
			    sprintf(log_string, "ERROR reading y from posfile at %d %d %d %d",
				    i, j, k, l);
			    p_log_errorno(log_string);
			  }
			  if ((objects_array[(int)fam_x][(int)fam_y] == -1)&&((int)fam_y!=1000 && (int)fam_x!=1000)){
				objects_array[(int)fam_x][(int)fam_y] = 1;
				total_object++;
			  }
		}
		l++;
	}

	convert=-1;

	// make sure we are where we think we are in the position file
	if( (fread(&temp, sizeof(short unsigned int), 1, posfile)) < 1){
	  p_log_errorno((char*)"ERROR reading curr_img from posfile");
	}
	if(temp!=0){
	  sprintf(log_string, "ERROR: 1expecting delimiter value of 0, read %d", temp);
	  p_log(log_string);
	  exit(0);
	}

	// make sure we are where we think we are in the position file
	if( (fread(&temp, sizeof(short unsigned int), 1, cy3_posfile)) < 1){
	  p_log_errorno((char*)"ERROR reading curr_img from posfile");
	}
	if(temp!=0){
	  sprintf(log_string, "ERROR: 2expecting delimiter value of 0, read %d", temp);
	  p_log(log_string);
	  exit(0);
	}

	// make sure we are where we think we are in the position file
	if( (fread(&temp, sizeof(short unsigned int), 1, txred_posfile)) < 1){
	  p_log_errorno((char*)"ERROR reading curr_img from posfile");
	}
	if(temp!=0){
	  sprintf(log_string, "ERROR: 3expecting delimiter value of 0, read %d", temp);
	  p_log(log_string);
	  exit(0);
	}

	// make sure we are where we think we are in the position file
	if( (fread(&temp, sizeof(short unsigned int), 1, fam_posfile)) < 1){
	  p_log_errorno((char*)"ERROR reading curr_img from posfile");
	}
	if(temp!=0){
	  sprintf(log_string, "ERROR: 4expecting delimiter value of 0, read %d", temp);
	  p_log(log_string);
	  exit(0);
	}


	//p_log((char*)"Write data to object table");
	if(fwrite(&convert, sizeof(int), 1, final_posfile) < 1){
	p_log_errorno((char*)"write of -1 to posfile failed");
	exit(0);
	}
	if(fwrite(&i, sizeof(int), 1, final_posfile) < 1){
	p_log_errorno((char*)"write of curr_fcnum to posfile failed");
	exit(0);
	}
	if(fwrite(&j, sizeof(int), 1, final_posfile) < 1){
	p_log_errorno((char*)"write of curr_arraynum to posfile failed");
	exit(0);
	}
	if(fwrite(&k, sizeof(int), 1, final_posfile) < 1){
	p_log_errorno((char*)"write of curr_imagenum to posfile failed");
	exit(0);
	}
	convert = (int)total_object;
	if(fwrite(&convert, sizeof(int), 1, final_posfile) < 1){
	p_log_errorno((char*)"write of num_beads to posfile failed");
	exit(0);
	}
        object_count=0;
	for(q=0;q<1000;q++){
	   for(p=0;p<1000;p++){
		if(objects_array[p][q] == 1){
			if(fwrite(&p, sizeof(short unsigned int), 1, final_posfile) < 1){
				p_log_errorno((char*)"write of beadpos_xcol to posfile failed");
				exit(0);
			}
			if(fwrite(&q, sizeof(short unsigned int), 1, final_posfile) < 1){
				p_log_errorno((char*)"write of beadpos_yrow to posfile failed");
				exit(0);
			}
			object_count++;
		}
	    }
	}

	//sprintf(log_string, "Position %d %d %d... %d cy5_objects, %d cy3_objects, %d txred_objects, %d fam_objects, total_object %d, total_object_count %d\n", i, j, k, curr_numobjs,cy3_curr_numobjs,txred_curr_numobjs,fam_curr_numobjs,total_object,object_count);
	sprintf(log_string, "Position %d %d %d... %d cy5_objects, %d cy3_objects, %d txred_objects, %d fam_objects, total_object %lld, total_object_count %d\n", i, j, k, curr_numobjs,cy3_curr_numobjs,txred_curr_numobjs,fam_curr_numobjs,total_object,object_count);
	p_log_simple(log_string);
	fprintf(obj_logfile,log_string);

	sconvert=0;
	if(fwrite(&sconvert, sizeof(short unsigned int), 1, final_posfile) < 1){
		p_log_errorno((char*)"write of 0 to posfile failed");
		exit(0);
	}
	//write to the info file
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
	sconvert = (int)total_object;
	if(fwrite(&sconvert, sizeof(short unsigned int), 1, infofile) < 1){
	p_log_errorno((char*)"write of num_beads to infofile failed");
	exit(0);
	}
      }
    }
  }
  fclose(posfile);
  fclose(cy3_posfile);
  fclose(fam_posfile);
  fclose(txred_posfile);
  fclose(final_posfile);
  fclose(infofile);
  fclose(obj_logfile);

//merge the autoexp object table

  if((cy3_posfile=fopen("autoexp_pos_cy3", "rb")) == NULL){
    sprintf(log_string, "ERROR opening position file %s for input", "autoexp_pos_cy3");
    p_log_errorno(log_string);
    exit(0);
  }

  if( (txred_posfile=fopen("autoexp_pos_txr", "rb")) == NULL){
    sprintf(log_string, "ERROR opening position file %s for input", "autoexp_pos_txred");
    p_log_errorno(log_string);
    exit(0);
  }

  if( (fam_posfile=fopen("autoexp_pos_fam", "rb")) == NULL){
    sprintf(log_string, "ERROR opening position file %s for input", "autoexp_pos_fam");
    p_log_errorno(log_string);
    exit(0);
  }

  if( (posfile=fopen("autoexp_pos_cy5", "rb")) == NULL){
    sprintf(log_string, "ERROR opening position file %s for input", "autoexp_pos_cy5");
    p_log_errorno(log_string);
    exit(0);
  }

  if( (final_posfile=fopen("autoexp_pos", "wb")) == NULL){
    sprintf(log_string, "ERROR opening position file %s for input", argv[1]);
    p_log_errorno(log_string);
    exit(0);
  }
  
  strcpy(info_filename,"autoexp_pos");
  strcat(info_filename,".info");
  if( (infofile=fopen(info_filename, "wb")) == NULL){
    sprintf(log_string, "ERROR opening reg file %s for output", info_filename);
    p_log_errorno(log_string);
    exit(0);
  }


        i = 0;
	j = 0;
        k = 0;
    for(j=0; j<ARRAYS_PER_FC; j++){
      for(k = 0;k<15;k++){
	if( (fread(&curr_fc, sizeof(int), 1, posfile)) < 1){
	  p_log_errorno((char*)"ERROR reading autoe curr_fc(1) from posfile");
	}
	if( (fread(&curr_fc, sizeof(int), 1, posfile)) < 1){
	  p_log_errorno((char*)"ERROR reading autoe curr_fc(2) from posfile");
	}
	if( (fread(&curr_array, sizeof(int), 1, posfile)) < 1){
	  p_log_errorno((char*)"ERROR reading autoe curr_array from posfile");
	}
	if( (fread(&curr_img, sizeof(int), 1, posfile)) < 1){
	  p_log_errorno((char*)"ERROR reading curr_img from posfile");
	}

	if( (fread(&curr_numobjs, sizeof(int), 1, posfile)) < 1){
	  p_log_errorno((char*)"ERROR reading curr_numobjs from posfile");
	}

	if( (fread(&curr_fc, sizeof(int), 1, cy3_posfile)) < 1){
	  p_log_errorno((char*)"ERROR reading curr_fc(1) from posfile");
	}
	if( (fread(&curr_fc, sizeof(int), 1, cy3_posfile)) < 1){
	  p_log_errorno((char*)"ERROR reading curr_fc(2) from posfile");
	}
	if( (fread(&curr_array, sizeof(int), 1, cy3_posfile)) < 1){
	  p_log_errorno((char*)"ERROR reading curr_array from posfile");
	}
	if( (fread(&curr_img, sizeof(int), 1, cy3_posfile)) < 1){
	  p_log_errorno((char*)"ERROR reading curr_img from posfile");
	}

	if( (fread(&curr_fc, sizeof(int), 1, txred_posfile)) < 1){
	  p_log_errorno((char*)"ERROR reading curr_fc(1) from posfile");
	}
	if( (fread(&curr_fc, sizeof(int), 1, txred_posfile)) < 1){
	  p_log_errorno((char*)"ERROR reading curr_fc(2) from posfile");
	}
	if( (fread(&curr_array, sizeof(int), 1, txred_posfile)) < 1){
	  p_log_errorno((char*)"ERROR reading curr_array from posfile");
	}
	if( (fread(&curr_img, sizeof(int), 1, txred_posfile)) < 1){
	  p_log_errorno((char*)"ERROR reading curr_img from posfile");
	}

	if( (fread(&curr_fc, sizeof(int), 1, fam_posfile)) < 1){
	  p_log_errorno((char*)"ERROR reading curr_fc(1) from posfile");
	}
	if( (fread(&curr_fc, sizeof(int), 1, fam_posfile)) < 1){
	  p_log_errorno((char*)"ERROR reading curr_fc(2) from posfile");
	}
	if( (fread(&curr_array, sizeof(int), 1, fam_posfile)) < 1){
	  p_log_errorno((char*)"ERROR reading curr_array from posfile");
	}
	if( (fread(&curr_img, sizeof(int), 1, fam_posfile)) < 1){
	  p_log_errorno((char*)"ERROR reading curr_img from posfile");
	}

	max_curr_numobjs = 0;
	if(max_curr_numobjs<curr_numobjs) max_curr_numobjs=curr_numobjs; 

	if( (fread(&cy3_curr_numobjs, sizeof(int), 1, cy3_posfile)) < 1){
	  p_log_errorno((char*)"ERROR reading curr_numobjs from posfile");
	}
	if(max_curr_numobjs<cy3_curr_numobjs) max_curr_numobjs=cy3_curr_numobjs;
	if( (fread(&txred_curr_numobjs, sizeof(int), 1, txred_posfile)) < 1){
	  p_log_errorno((char*)"ERROR reading curr_numobjs from posfile");
	}
	if(max_curr_numobjs<txred_curr_numobjs) max_curr_numobjs=txred_curr_numobjs;
	if( (fread(&fam_curr_numobjs, sizeof(int), 1, fam_posfile)) < 1){
	  p_log_errorno((char*)"ERROR reading curr_numobjs from posfile");
	}
	if(max_curr_numobjs<fam_curr_numobjs) max_curr_numobjs=fam_curr_numobjs;

        l = 0;
	m = 0;
        n = 0;
	o = 0;
        p = 0;
	q = 0;
	cy5_proceed = 1;
	cy3_proceed = 1;
	txred_proceed = 1;
	fam_proceed = 1;
	found_object = 0;
	total_object = 0;

	for(p=0;p<2000;p++){
	   for(q=0;q<2000;q++){
		objects_array[p][q] = -1;
		}
	}

	while(l<max_curr_numobjs)
	{

		if(l<curr_numobjs)
		{

			  if((fread(&cy5_x, sizeof(short unsigned int), 1, posfile)) < 1){
			    sprintf(log_string, "ERROR reading x from autoe posfile at %d %d %d %d",
				    i, j, k, l);
			    p_log_errorno(log_string);
			  }
			  if((fread(&cy5_y, sizeof(short unsigned int), 1, posfile)) < 1){
			    sprintf(log_string, "ERROR reading y from autoe posfile at %d %d %d %d",
				    i, j, k, l);
			    p_log_errorno(log_string);
			  }

			  if ((objects_array[(int)cy5_x][(int)cy5_y] == -1)&&((int)cy5_y!=1000 && (int)cy5_x!=1000)){
				objects_array[(int)cy5_x][(int)cy5_y] = 1;
				total_object++;
			  }

		}
		if(l<cy3_curr_numobjs)
		{
			  if((fread(&cy3_x, sizeof(short unsigned int), 1, cy3_posfile)) < 1){
			    sprintf(log_string, "ERROR reading x from autoe posfile at %d %d %d %d",
				    i, j, k, l);
			    p_log_errorno(log_string);
			  }

			  if((fread(&cy3_y, sizeof(short unsigned int), 1, cy3_posfile)) < 1){
			    sprintf(log_string, "ERROR reading y from autoe posfile at %d %d %d %d",
				    i, j, k, l);
			    p_log_errorno(log_string);
			  }

			  if ((objects_array[(int)cy3_x][(int)cy3_y] == -1)&&((int)cy3_y!=1000 && (int)cy3_x!=1000)){
				objects_array[(int)cy3_x][(int)cy3_y] = 1;
				total_object++;
			  }
		}

		if(l<txred_curr_numobjs)
		{
			  if((fread(&txred_x, sizeof(short unsigned int), 1, txred_posfile)) < 1){
			    sprintf(log_string, "ERROR reading x from autoe posfile at %d %d %d %d",
				    i, j, k, l);
			    p_log_errorno(log_string);
			  }
			  if((fread(&txred_y, sizeof(short unsigned int), 1, txred_posfile)) < 1){
			    sprintf(log_string, "ERROR reading y from autoe posfile at %d %d %d %d",
				    i, j, k, l);
			    p_log_errorno(log_string);
			  }
			  if ((objects_array[(int)txred_x][(int)txred_y] == -1)&&((int)txred_y!=1000 && (int)txred_x!=1000)){
				objects_array[(int)txred_x][(int)txred_y] = 1;
				total_object++;
			  }
		}

		if(l<fam_curr_numobjs)
		{
			  if((fread(&fam_x, sizeof(short unsigned int), 1, fam_posfile)) < 1){
			    sprintf(log_string, "ERROR reading x from autoe posfile at %d %d %d %d",
				    i, j, k, l);
			    p_log_errorno(log_string);
			  }
			  if((fread(&fam_y, sizeof(short unsigned int), 1, fam_posfile)) < 1){
			    sprintf(log_string, "ERROR reading y from autoe posfile at %d %d %d %d",
				    i, j, k, l);
			    p_log_errorno(log_string);
			  }
			  if ((objects_array[(int)fam_x][(int)fam_y] == -1)&&((int)fam_y!=1000 && (int)fam_x!=1000)){
				objects_array[(int)fam_x][(int)fam_y] = 1;
				total_object++;
			  }
		}
		l++;
	}

	// make sure we are where we think we are in the position file
	if( (fread(&temp, sizeof(short unsigned int), 1, posfile)) < 1){
	  p_log_errorno((char*)"ERROR reading curr_img from autoe posfile");
	}
	if(temp!=0){
	  sprintf(log_string, "ERROR: 1expecting delimiter value of 0, read %d", temp);
	  p_log(log_string);
	  exit(0);
	}

	// make sure we are where we think we are in the position file
	if( (fread(&temp, sizeof(short unsigned int), 1, cy3_posfile)) < 1){
	  p_log_errorno((char*)"ERROR reading curr_img from autoe posfile");
	}
	if(temp!=0){
	  sprintf(log_string, "ERROR: 2expecting delimiter value of 0, read %d", temp);
	  p_log(log_string);
	  exit(0);
	}

	// make sure we are where we think we are in the position file
	if( (fread(&temp, sizeof(short unsigned int), 1, txred_posfile)) < 1){
	  p_log_errorno((char*)"ERROR reading curr_img from posfile");
	}
	if(temp!=0){
	  sprintf(log_string, "ERROR: 3expecting delimiter value of 0, read %d", temp);
	  p_log(log_string);
	  exit(0);
	}

	// make sure we are where we think we are in the position file
	if( (fread(&temp, sizeof(short unsigned int), 1, fam_posfile)) < 1){
	  p_log_errorno((char*)"ERROR reading curr_img from posfile");
	}
	if(temp!=0){
	  sprintf(log_string, "ERROR: 4expecting delimiter value of 0, read %d", temp);
	  p_log(log_string);
	  exit(0);
	}


	convert=-1;
	//p_log((char*)"Write data to object table");
	if(fwrite(&convert, sizeof(int), 1, final_posfile) < 1){
	p_log_errorno((char*)"write of -1 to posfile failed");
	exit(0);
	}
	if(fwrite(&i, sizeof(int), 1, final_posfile) < 1){
	p_log_errorno((char*)"write of curr_fcnum to posfile failed");
	exit(0);
	}
	if(fwrite(&j, sizeof(int), 1, final_posfile) < 1){
	p_log_errorno((char*)"write of curr_arraynum to posfile failed");
	exit(0);
	}
	if(fwrite(&k, sizeof(int), 1, final_posfile) < 1){
	p_log_errorno((char*)"write of curr_imagenum to posfile failed");
	exit(0);
	}
	convert = (int)total_object;
	if(fwrite(&convert, sizeof(int), 1, final_posfile) < 1){
	p_log_errorno((char*)"write of num_beads to posfile failed");
	exit(0);
	}
        object_count=0;
	for(q=0;q<1000;q++){
	   for(p=0;p<1000;p++){
		if(objects_array[p][q] == 1){
			if(fwrite(&p, sizeof(short unsigned int), 1, final_posfile) < 1){
				p_log_errorno((char*)"write of beadpos_xcol to posfile failed");
				exit(0);
			}
			if(fwrite(&q, sizeof(short unsigned int), 1, final_posfile) < 1){
				p_log_errorno((char*)"write of beadpos_yrow to posfile failed");
				exit(0);
			}
			object_count++;
		}
	    }
	}

	sprintf(log_string, "autoe1 Position %d %d %d... %d cy5_objects, %d cy3_objects, %d txred_objects, %d fam_objects, total_object %lld, total_object_count %d\n", curr_fc, curr_array, curr_img, curr_numobjs,cy3_curr_numobjs,txred_curr_numobjs,fam_curr_numobjs,total_object,object_count);
	p_log_simple(log_string);

	sconvert=0;
	if(fwrite(&sconvert, sizeof(short unsigned int), 1, final_posfile) < 1){
		p_log_errorno((char*)"write of 0 to posfile failed");
		exit(0);
	}
	//write to the info file
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
	sconvert = (int)total_object;
	if(fwrite(&sconvert, sizeof(short unsigned int), 1, infofile) < 1){
	p_log_errorno((char*)"write of num_beads to infofile failed");
	exit(0);
	}
      }
    }
	
  fclose(posfile);
  fclose(cy3_posfile);
  fclose(fam_posfile);
  fclose(txred_posfile);
  fclose(final_posfile);
  fclose(infofile);
  //system("rm 00_*");
}
