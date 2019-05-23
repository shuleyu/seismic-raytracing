# Download and Install
'''
$ git clone --recursive https://github.com/shuleyu/seismic-raytracing.git
'''

If "--recursive" is not added, the dependencies will not be downloaded. In such case, do:

$ cd ./seismic-raytracing

$ git submodule update --init --recursive

# Parameter file (INFILE)
The WORKDIR parameter specify the output folder of this program. The program will create this folder if it doesn't exist.

$ vim ./INFILE

# Task file (LIST.sh)
Specify which task to run.

$ vim ./LIST.sh

# Execution file (Run.sh)
When INFILE and LIST.sh are properly set, run this script:

$ ./Run.sh

# Example
A ScS reflection on an ellipse-shaped low velocity structure on the core mantle boundary:

![alt text](https://github.com/shuleyu/raytracing/blob/master/SRC/example2.png)

Ray is coming from left hand side. Green lines for polarity change. Line width represent displacement amplitude.
