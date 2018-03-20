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
OUTFILE=tmp.ps
RE="6371.0"
PROJ="-JPa${PLOTSIZE}i/60"
REG="-R0/360/0/${RE}"

# Move to center.
psxy ${PROJ} ${REG} -P -K -Xc -Yc > ${OUTFILE} << EOF
EOF

# Move to CenterAt.
X=`echo ${PLOTSIZE} ${CenterAt}| awk '{print $3*cos((30+$2)/180*3.1415926)/6371.0/2*$1}'`
Y=`echo ${PLOTSIZE} ${CenterAt}| awk '{print -$3*sin((30+$2)/180*3.1415926)/6371.0/2*$1}'`
psxy -J -R -O -K -X${X}i -Y${Y}i >> ${OUTFILE} << EOF
EOF

# distance grid.
psbasemap ${PROJ} ${REG} -Ba10f2 -O -K >> ${OUTFILE}

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
	grep -n ">" ${file} | awk 'BEGIN {FS=":"} {print $1}' > lines1
	awk '{print $1-1}' lines1 | awk 'NR>1 {print $0}' > lines2
	wc -l < ${file} >> lines2
	paste lines1 lines2 > lines
	
	while read l1 l2
	do
		[ "`awk -v L1=${l1} 'NR==L1 {print $2}' ${file}`" = "P" ] && Pen="-W0.5p,blue" ||  Pen="-W0.5p,red"
		awk -v L1=${l1} -v L2=${l2} '{ if (NR>L1 && NR<=L2) print $0}' ${file} | psxy ${PROJ} ${REG} -m -O -K ${Pen} >> ${OUTFILE}
	done < lines
done


# plot source.
psxy ${PROJ} ${REG} -Sa0.2i -Gyellow -N -O -K >> ${OUTFILE} << EOF
0 `echo "${RE}-${EVDE}" | bc -l`
EOF

# plot velocity anomalies.
psxy ${PROJ} ${REG} ${WORKDIR}/tmpfile_Polygons_${RunNumber} -m -L -W1p,black -O -K  >> ${OUTFILE}

# seal the plot.
psxy -J -R -O >> ${OUTFILE} << EOF
EOF

# Make PDF.
Title=`basename $0`
ps2pdf tmp.ps ${PLOTDIR}/${Title%.sh}.pdf
tomini ${PLOTDIR}/${Title%.sh}.pdf

exit 0
