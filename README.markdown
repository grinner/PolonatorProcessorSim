# Polonator G.007 Image Processing simulator
Nick Conway, Wyss Institute 5/11/2011

## Release Background
This is a simulator for the original Polonator processing pipeline.
It does the following

1. removes the original TCP/IP communication with the acquisition machine
2. Tidies up the code for readability and organization
3. uses a data file set from the Polonator G.007 for input to create the object
table and does a run of sequencing
4. Use `creatImageDataSet.py` to create the input text files by running it in 
each `*.raw` image directory to create an entire set of input images.

# Download

    git clone git@github.com:grinner/PolonatorProcessorSim.git

Get a pre-indexed group of test images at:
Directory formats are:
    `AM[N][n]_X`

where 

* N is a number signaling the read position of an objects sequence
* n is a letter (a,b,c ...) signaling the number of times that the above read position has been sequenced. (a=1st, b=2nd, c=3rd...)
* X is the base color being interrogated (A,C,G,T).  

Note: during object finding the emission fluorophore is mapped to a base by the following:

* A = cy3 
* C = txred, Texas red
* G = fam
* T = cy5

[Randomized Rolonies](http://dl.dropbox.com/u/15200497/raw_img.zip)
These images represent an actual sequence run.  The image set consists of one object finding cycle and one sequencing cycle.  Only 1 in 100 image positions are recorded

[Rolonies on a grid](http://dl.dropbox.com/u/15200497/images/rimage/grid_images_021011_rolony.tar.gz), Note: this data set needs further explanation on how to run.
These images represent an actual sequencing run.  
The image set consists of one object finding cycle and 6 sequencing cycles.
Additionally the 1st position is a also sequenced 3 times 
Only 1 in 100 image positions are recorded 
There is a README file in the `*.tar.gz` file to describe this further

## Usage

Move the above test images into a subdirectory of the cloned repository named `raw_img` as follows:

    raw_img/
        00_00_cy3/
        00_00_cy5/
        00_00_fam/
        00_00_txr/
        AM1a_A/
        AM1a_C/
        AM1a_G/
        AM1a_T/
    
    
In addition, make sure you create the following subdirectories

    autoexp_images/
    beads/
    bin/
    images/
    logs/
    qc/
    tetrahedra/
    
Run these files in this directory with: 

    python launcher.py 