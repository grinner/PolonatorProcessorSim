#!/usr/bin/env python
# encoding: utf-8
"""
createImageDataSet.py

Created by Nicholas Conway on 2011-03-14.
Copyright (c) 2011 __Wyss_Instiute__. All rights reserved.

This script is for the polonator to generate an image dataset that 
represents the input to the system from a saved dataset.  Since the polonator
only saves 1 in every 100 images, we need to mulitply that image by 100.

Use: usage 
python createImageDataSet 
in the directory of your choosing and it will 
generate a txt file with a list of files names to load.

"""

import sys
import os
import glob

MULTIPLE = 100
MULTIPLE2 = 80

def main():
    multip = MULTIPLE
    path = os.getcwd()
    file_list = glob.glob('*.raw')
    file_list.sort()
    out_file = open('image_set.txt', 'w')
    for item in file_list:
        if item.find('2100') != -1:  
            multip = MULTIPLE2
        # end if
        else:
            multip = MULTIPLE
        # end else         
        for i in range(multip):
            out_file.write(item + '\n')
        # end for
    # end for
    out_file.close()
# end def


if __name__ == '__main__':
    main()

