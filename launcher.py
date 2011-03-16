#!/usr/bin/env python
# encoding: utf-8
"""
launcher.py

Created by Nicholas Conway on 2011-03-14.
Copyright (c) 2011 __Wyss_Instiute__. All rights reserved.

This script is for the polonator to run the image processing pipeline


"""

import os


os.system("rm 0*")
os.system("rm 1*")
os.system("./4color_initialize_processor.pl; ./processor.pl")

