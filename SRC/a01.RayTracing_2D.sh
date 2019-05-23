#!/bin/bash

# ==============================================================
# This script runs the 2D ray-tracing program.
#
# Shule Yu
# Spet 11 2014
# ==============================================================

echo ""
echo "--> `basename $0` is running."
rm -f ${WORKDIR}/${RayFilePrefix}* ${WORKDIR}/${PolygonFilePrefix}*
trap "rm -f ${WORKDIR}/tmpfile*$$ ${WORKDIR}/*_${RunNumber}; exit 1" SIGINT

# ==============================================
#           ! Work Begin !
# ==============================================

# C++ code.
${EXECDIR}/TraceIt.out 6 8 0 << EOF
${DebugInfo}
${TS}
${TD}
${RD}
${StopAtSurface}
${nThread}
${WORKDIR}/tmpfile_InputRays_${RunNumber}
${WORKDIR}/tmpfile_LayerSetting_${RunNumber}
${WORKDIR}/tmpfile_KeyDepths_${RunNumber}
${WORKDIR}/tmpfile_1DRef_${RunNumber}
${WORKDIR}/tmpfile_Polygons_${RunNumber}
${WORKDIR}/${RayFilePrefix}
${WORKDIR}/${ReceiverFileName}
${WORKDIR}/${PolygonFilePrefix}
EOF

[ $? -ne 0 ] && echo "C++ code Failed ..." && rm -f tmpfile*$$ && exit 1

# Clean up.
rm -f tmpfile*$$

cd ${WORKDIR}

exit 0
