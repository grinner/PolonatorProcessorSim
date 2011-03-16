      // tell server we're ready for the next block
      p_log((char*)"Request data");
      
      if(send(sock, sendBuffer, 1, 0)!=1){
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
      while(totalBytesRcvd < data_size){

        // was there an error during the recv?  if so, this is bad; crash
        if((bytesRcvd = recv(sock, ((char*)inputBuffer) + totalBytesRcvd, data_size-totalBytesRcvd, MSG_WAITALL)) < 0){
        /*if((bytesRcvd = recv(sock, ((char*)inputBuffer) + totalBytesRcvd, data_size-totalBytesRcvd, MSG_WAITALL)) <= 0){*/
          /*      if((bytesRcvd = recv(sock, ((char*)inputBuffer) + totalBytesRcvd, data_size-totalBytesRcvd, 0)) <= 0){*/
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
      p_log((char*)"Data received");
      if(!connection_isok){
        i=curr_fcnum+1;
        j=ARRAYS_PER_FC;
        k=IMGS_PER_ARRAY;
        p_log((char*)"ERROR:\tReceiveData: connection was broken; exiting");
        break;
      }

      
      // NOW THAT WE HAVE THE IMAGE, DO SOMETHING WITH IT
      //verify it is the image we think it is
      if( (*(inputBuffer + 0) != i) || 
          (*(inputBuffer + 1) != j) ||
          (*(inputBuffer + 2) != k) || 
          (*(inputBuffer + 3) == 1)){
        sprintf(log_string, "ERROR; expected image %d %d %d 0, got %d %d %d %d",
            i, j, k,
            *(inputBuffer),
            *(inputBuffer+1),
            *(inputBuffer+2),
            *(inputBuffer+3));
        p_log(log_string);
        exit(0);
      }