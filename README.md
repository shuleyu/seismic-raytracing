# Download and Install
```
$ git clone --recursive https://github.com/shuleyu/seismic-raytracing.git
```

If `--recursive` is not added, the dependencies will not be downloaded. In such case, do:

```
$ cd ./seismic-raytracing
$ git submodule update --init --recursive
```


# Parameter file (INFILE)
The WORKDIR parameter specify the output folder of this program. The program will create this folder if it doesn't exist.

```
$ vim ./INFILE
```

# Task file (LIST.sh)
Specify which task(s) to execute.

```
$ vim ./LIST.sh
```

# Execution file (Run.sh)
When INFILE and LIST.sh are properly set, run this script:

```
$ bash ./Run.sh
```

Run.sh will try to compile the code and run the task listed in LIST.sh using the parameters in INFILE.

With a proper WORKDIR (and gmt4 installed), by default parameters in INFILE, it should calculate the following example case and produce the example figure.

# Example
A ScS reflection on an ellipse-shaped low velocity structure on the core mantle boundary; ray is coming from left hand side:

![alt text](https://github.com/shuleyu/raytracing/blob/master/SRC/example1.png)

Green lines for polarity change. Line width are comparable to displacement amplitude.

Multiple ScS reflection on this structure (varying takeoff angle):

![alt text](https://github.com/shuleyu/raytracing/blob/master/SRC/example2.png)

