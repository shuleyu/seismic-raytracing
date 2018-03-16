#!/bin/bash

# ==============================================================
# This script calculate the 1D ray-tracing result.
#
# Shule Yu
# Spet 11 2014
# ==============================================================

echo ""
echo "--> `basename $0` is running."
mkdir -p ${WORKDIR}
cd ${WORKDIR}
trap "rm -f ${WORKDIR}/tmpfile*$$ ${WORKDIR}/*_${RunNumber}; exit 1" SIGINT

# ==============================================
#           ! Work Begin !
# ==============================================

# C++ code.
${EXECDIR}/ConstructRayPath.out 0 3 3 << EOF
${WORKDIR}/tmpfile_LayerSetting_${RunNumber}
${WORKDIR}/tmpfile_RayPathSetting_${RunNumber}
${WORKDIR}/${OutFileName}
${TakeOffAngle}
${EVDE}
${CriticalAngle}
EOF

if [ $? -ne 0 ]
then
    echo "C++ code Failed ..."
    rm -f tmpfile*$$
    exit 1
fi

# Clean up.
rm -f tmpfile*$$

cd ${WORKDIR}

exit 0
