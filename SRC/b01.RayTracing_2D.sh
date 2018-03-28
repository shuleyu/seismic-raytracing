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
ls ${WORKDIR}/${OutFilePrefix}* >/dev/null 2>&1
[ $? -ne 0 ] && echo "    !=> In `basename $0`: Run a01 first ..." && exit 1

# Plot.
Rotate=`echo ${CenterAt} | awk '{print $1}'`
OUTFILE=tmp.ps
RE="6371.0"
PROJ="-JPa${PLOTSIZE}i/${Rotate}"
REG="-R0/360/0/${RE}"

# Move to center.
psxy ${PROJ} ${REG} -K -Xc -Yc > ${OUTFILE} << EOF
EOF

# Move to CenterAt.
X=`echo ${PLOTSIZE} ${CenterAt} ${Rotate}| awk '{print $3*cos(((90-$4)+$2)/180*3.1415926)/6371.0/2*$1}'`
Y=`echo ${PLOTSIZE} ${CenterAt} ${Rotate}| awk '{print -$3*sin(((90-$4)+$2)/180*3.1415926)/6371.0/2*$1}'`
psxy -J -R -O -K -X${X}i -Y${Y}i >> ${OUTFILE} << EOF
EOF

# distance grid.
[ `echo "${PLOTSIZE}>10" | bc` -eq 1 ] && BAXI="-Ba1f1" || BAXI="-Ba10f1"
psbasemap ${PROJ} ${REG} ${BAXI} -O -K >> ${OUTFILE}

# plot mantle.
psxy ${PROJ} ${REG} -Sc${PLOTSIZE}i -G230/230/230 -O -K >> ${OUTFILE} << EOF
0 0
EOF

# plot transition zone.
psxy ${PROJ} ${REG} -Sc`echo "${PLOTSIZE}/${RE}*(${RE}-660/2-410/2)"| bc -l`i -W`echo "${PLOTSIZE}/${RE}*(660-410)/2" |bc -l`i,200/200/200 -O -K >> ${OUTFILE} << EOF
0 0
EOF

# plot outter core.
psxy ${PROJ} ${REG} -Sc`echo "${PLOTSIZE}/${RE}*3480.0"| bc -l`i -G170/170/170 -O -K >> ${OUTFILE} << EOF
0 0
EOF

# plot inner core.
psxy ${PROJ} ${REG} -Sc`echo "${PLOTSIZE}/${RE}*1221.5"| bc -l`i -G140/140/140 -O -K >> ${OUTFILE} << EOF
0 0
EOF

# plot ray path.
for file in `ls ${WORKDIR}/${OutFilePrefix}*`
do
    RayNumber=${file##*_}
    RayColor=`grep -w ${RayNumber} ${WORKDIR}/${OutInfoFile} | awk 'NR==1 {print $1}'`
    if [ ${RayColor} = "black" ]
    then
        [ `head -n 1 ${file} | awk '{print $2}'` = "S" ] && RayColor="red" || RayColor="blue"
    fi
    psxy ${file} ${PROJ} ${REG} -W0.5p,${RayColor} -m -O -K >> ${OUTFILE}
done


# plot source.
awk '{print $1,$2}' ${WORKDIR}/tmpfile_InputRays_${RunNumber} | sort -u > tmpfile_sources_$$
while read theta depth
do
    psxy ${PROJ} ${REG} -Sa0.2i -Gyellow -N -O -K >> ${OUTFILE} << EOF
${theta} `echo "${RE}-${depth}" | bc -l`
EOF
done < tmpfile_sources_$$

# plot velocity anomalies.
for file in `ls ${WORKDIR}/${PolygonOutPrefix}*`
do
    psxy ${PROJ} ${REG} ${file} -m -L -W1p,black -O -K  >> ${OUTFILE}
done

# plot scale at the CMB.
PROJ2=`echo "${PLOTSIZE} ${Rotate}" | awk '{print "-JPa"$1*3480/6371"i/"$2}'`
MOVE=`echo "${PLOTSIZE}" |  awk '{print $1/2*2891/6371}'`
gmtset TICK_LENGTH = -0.2c
gmtset ANNOT_OFFSET_PRIMARY = -0.2c
psbasemap ${PROJ2} ${REG} ${BAXI} -X${MOVE}i -Y${MOVE}i -O -K >> ${OUTFILE}

# seal the plot.
psxy -J -R -O >> ${OUTFILE} << EOF
EOF

# Make PDF.
Title=`basename $0`
ps2pdf tmp.ps ${PLOTDIR}/${Title%.sh}.pdf
tomini ${PLOTDIR}/${Title%.sh}.pdf

exit 0
