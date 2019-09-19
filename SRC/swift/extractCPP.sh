#!/bin/bash


SRCDIR="./seismic-raytracing/SRC"

# 1. Some common header.
cat > cppLibrary.cpp << EOF
#include<vector>
#include<set>
#include<cmath>
#include<algorithm>
#include<complex>
#include<atomic>
#include<thread>
#include<unistd.h>
#include<string.h>
EOF

# 2. Ray.hpp

cat ${SRCDIR}/Ray.hpp | grep -v "#include<" | grep -v "std::runtime_error" | grep -v "std::cerr" | grep -v "std::endl">> cppLibrary.cpp

# 3. cat the Library Files.
while read hppFile
do
    cat ${SRCDIR}/CPP-Library/${hppFile} | grep -v "#include<" | grep -v "std::runtime_error" | grep -v "std::cerr" | grep -v "std::endl">> cppLibrary.cpp
done<< EOF
Lon2180.hpp
Lon2360.hpp
LL2TP.hpp
Sph2Cart.hpp
DotDist.hpp
CrossProduct.hpp
CreateGrid.hpp
LineJunction.hpp
LocDist.hpp
PREM.hpp
PlaneWaveCoefficients.hpp
PointOnSegment.hpp
PointInPolygon.hpp
RayPath.hpp
SegmentJunction.hpp
EOF


# 4. RayTracing.fun.cpp

cat ${SRCDIR}/RayTracing.fun.cpp | grep -v "#include<" | grep -v "runtime_error" | grep -v "cout" | grep -v "endl" | grep -v -w "printf" >> cppLibrary.cpp

# 5. Check

c++ -c cppLibrary.cpp


exit 0
