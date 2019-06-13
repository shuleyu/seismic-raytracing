#!/bin/bash

# ==============================================================
# This script plot the 2D ray-tracing result.
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

color[1]=darkred
color[2]=green
color[3]=lightblue
color[4]=purple
color[5]=lightgreen
color[6]=cyan
color[7]=darkblue
color[8]=gold
color[9]=yellow

# ================================================
#         ! Work Begin !
# ================================================

# Check Calculation
ls ${WORKDIR}/${RayFilePrefix}* >/dev/null 2>&1
[ $? -ne 0 ] && echo "    !=> In `basename $0`: Run a01 first ..." && exit 1
RE="6371.0"

# Plot.
Cnt=0
while read PLOTSIZE Theta Radius
do
    Cnt=$((Cnt+1))
    OUTFILE=${Cnt}.ps

    PROJ="-JPa${PLOTSIZE}i/${Theta}"
    REG="-R0/360/0/${RE}"

    # Move to center.
    psxy ${PROJ} ${REG} -K -Xc -Yc > ${OUTFILE} << EOF
EOF

    # Move to the plot point.
    Y=`echo ${PLOTSIZE} ${Theta} ${Radius} ${RE}| awk '{print ($3-$4)/$4/2*$1}'`
    psxy -J -R -O -K -Y${Y}i >> ${OUTFILE} << EOF
EOF

    # distance grid.
    [ `echo "${PLOTSIZE}>10" | bc` -eq 1 ] && BAXI="-Ba1f1" || BAXI="-Ba10f1"
    gmtset TICK_LENGTH = 0.2c
    gmtset ANNOT_OFFSET_PRIMARY = 0.2c
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
    for file in `ls ${WORKDIR}/${RayFilePrefix}*`
    do

        # Plot choice. Rays with segments larger than 90 degree is not plotted.
#         FinalDist=`tail -n 1 ${file} | awk '{if ($1>90) print 0; else print 1}'`
#         [ ${FinalDist} -eq 0 ] && continue

        RayColor=`head -n 1 ${file} | awk '{print $2}'`
        RayType=`head -n 1 ${file} | awk '{print $3}'`
        if [ ${RayColor} = "0" ]
        then
            [ ${RayType} = "P" ] && plotColor="blue" || plotColor="red"
        else 
            plotColor=${color[$RayColor]}
        fi

        Amp=`head -n 1 ${file} | awk '{if ($8<0) A=-1; else A=1; printf "%.6lf",A*sqrt(A*$8)*4}'`
        [ `echo "${Amp}<0" | bc` -eq 1 ] && plotColor="darkgreen"
        Amp=`echo ${Amp} | awk '{if ($1<0) printf "%.12f", -$1; else printf "%.12f", $1}'`
        [ `echo "${Amp}<0.1" | bc` -eq 1 ] && Amp=0.1
        [ ${LineThicknessUseAmp} -ne 1 ] && Amp=0.5

        psxy ${file} ${PROJ} ${REG} -W${Amp}p,${plotColor} -m -O -K >> ${OUTFILE}
    done

    for file in `ls ${WORKDIR}/${RayFilePrefix}*`
    do
        # Plot choice. Rays with segments larger than 90 degree is not plotted.
#         FinalDist=`tail -n 1 ${file} | awk '{if ($1>90) print 0; else print 1}'`
#         [ ${FinalDist} -eq 0 ] && continue

        BeginR=`awk 'NR==2 {print $2}' ${file}`
        EndR=`tail -n 1 ${file} | awk '{print $2}'`
        TargetR=`echo ${BeginR} ${EndR} | awk '{print $1+($2-$1)*0.1}'`
        Position=`awk -v T=${TargetR} 'NR>1 {if ($2>T) print $1,$2-T; else print $1,T-$2}' ${file} | sort -g -k 2,2 | head -n 1 | awk '{print $1}'`

        pstext ${PROJ} ${REG} -N -O -K >> ${OUTFILE} << EOF
${Position} ${TargetR} 5 0 0 CM ${file##*_}
EOF
    done

    # plot info at surface.

    if [ `echo "${PLOTSIZE}>10" | bc` -eq 1 ]
    then
        ${SRCDIR}/FindColumn.sh ${WORKDIR}/${ReceiverFileName} "<Dist> <TravelTime>" > tmpfile_$$
        while read dist tt
        do
            pstext ${PROJ} ${REG} -N -O -K >> ${OUTFILE} << EOF
${dist} `echo ${PLOTSIZE} ${RE} | awk '{print $2+$2*2/$1*0.5}'` 15 0 0 LB ${dist} deg. ${tt} sec.
EOF
        done < tmpfile_$$
    fi

    # plot source.
    awk '{print $1,$2}' ${WORKDIR}/tmpfile_InputRays_${RunNumber} | sort -u > tmpfile_sources_$$
    while read theta depth
    do
        psxy ${PROJ} ${REG} -Sa0.2i -Gyellow -N -O -K >> ${OUTFILE} << EOF
${theta} `echo "${RE}-${depth}" | bc -l`
EOF
    done < tmpfile_sources_$$

    # plot velocity anomalies.
    for file in `ls ${WORKDIR}/${PolygonFilePrefix}* 2>/dev/null`
    do
        psxy ${PROJ} ${REG} ${file} -m -L -W0.3p,black -O -K  >> ${OUTFILE}
    done

    # plot basemap at the CMB.
    PROJ2=`echo "${PLOTSIZE} ${Theta} ${RE}" | awk '{print "-JPa"$1*3480/$3"i/"$2}'`
    MOVE=`echo "${PLOTSIZE} ${RE}" |  awk '{print $1/2*2891/$2}'`
    gmtset TICK_LENGTH = -0.2c
    gmtset ANNOT_OFFSET_PRIMARY = -0.2c
    psbasemap ${PROJ2} ${REG} ${BAXI} -X${MOVE}i -Y${MOVE}i -O -K >> ${OUTFILE}

    # seal the plot.
    psxy -J -R -O >> ${OUTFILE} << EOF
EOF

done < ${WORKDIR}/tmpfile_PlotWhere_${RunNumber} # done plot position loop.

# Make PDF.
Title=`basename $0`
cat `ls -rt *ps` > ${PLOTDIR}/${Title%.sh}.ps
ps2pdf ${PLOTDIR}/${Title%.sh}.ps ${PLOTDIR}/${Title%.sh}.pdf

exit 0
