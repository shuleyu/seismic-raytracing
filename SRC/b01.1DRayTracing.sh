#!/bin/bash

# ==============================================================
# This script plot the 1D ray-tracing result.
#
# Shule Yu
# Spet 11 2014
# ==============================================================

echo ""
echo "--> `basename $0` is running. "
mkdir -p ${PLOTDIR}/tmpdir_$$
cd ${PLOTDIR}/tmpdir_$$
trap "rm -rf ${PLOTDIR}/tmpdir_$$; exit 1" SIGINT EXIT

# Plot parameters.
gmtset PAPER_MEDIA = letter
gmtset ANNOT_FONT_SIZE_PRIMARY = 8p
gmtset LABEL_FONT_SIZE = 10p
gmtset LABEL_OFFSET = 0.05c

# ================================================
#         ! Work Begin !
# ================================================

# Check Calculation
if ! [ -s ${WORKDIR}/${OutFileName} ]
then
    echo "    !=> In `basename $0`: Run a01 first ..."
    exit 1
fi

# Plot.
OUTFILE=tmp.ps
RE="6371.0"
PROJ="-JPa${PLOTSIZE}i/60"
REG="-R0/360/0/${RE}"

# distance grid.
psbasemap ${PROJ} ${REG} -Ba10f2 -P -K -Xc -Yc > ${OUTFILE}

# plot mantle.
psxy ${PROJ} ${REG} -Sc${PLOTSIZE}i -Gorange -O -K >> ${OUTFILE} << EOF
0 0
EOF

# plot transition zone.
psxy ${PROJ} ${REG} -Sc`echo "${PLOTSIZE}/${RE}*(${RE}-660/2-410/2)"| bc -l`i -W`echo "${PLOTSIZE}/${RE}*(660-410)/2" |bc -l`i,olivedrab -O -K >> ${OUTFILE} << EOF
0 0
EOF

# plot outter core.
psxy ${PROJ} ${REG} -Sc`echo "${PLOTSIZE}/${RE}*3480.0"| bc -l`i -Gred -O -K >> ${OUTFILE} << EOF
0 0
EOF

# plot inner core.
psxy ${PROJ} ${REG} -Sc`echo "${PLOTSIZE}/${RE}*1221.5"| bc -l`i -Gyellow -O -K >> ${OUTFILE} << EOF
0 0
EOF

# plot ray path.
awk '{ if ($1!=">") print $1,$2 ; else print $1}' ${WORKDIR}/${OutFileName} | psxy ${PROJ} ${REG} -m -O -K >> ${OUTFILE}


# plot source.
psxy ${PROJ} ${REG} -Sa0.2i -Gyellow -N -O -K >> ${OUTFILE} << EOF
0 `echo "${RE}-${EVDE}" | bc -l`
EOF

# seal the plot.
psxy -J -R -O >> ${OUTFILE} << EOF
EOF

# Make PDF.
Title=`basename $0`
ps2pdf tmp.ps ${PLOTDIR}/${Title%.sh}.pdf
tomini ${PLOTDIR}/${Title%.sh}.pdf

exit 0
