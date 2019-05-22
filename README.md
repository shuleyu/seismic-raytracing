# Install
$ git clone --recursive https://github.com/shuleyu/seismic-raytracing.git
If "--recursive" is not added, the dependencies will not be cloned. In such cases do:
$ cd ./seismic-raytracing
$ git submodule update --init --recursive

# Input file
# Run the code
seismology 2D ray tracing. Use polygons to represent 2D regions.

Edit INFILE, execute Run.sh. Dependencies should be here: https://github.com/shuleyu/CPP-Library.

For example, the ScS reflection on an ellipse-shaped low velocity structure on the core mantle boundary:

![alt text](https://github.com/shuleyu/raytracing/blob/master/example2.png)

Ray comes from left (changing slowness for each subplot). red if polarity is positive, green means negative polarity.
