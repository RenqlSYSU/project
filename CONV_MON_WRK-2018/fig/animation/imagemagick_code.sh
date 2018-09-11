#!/bin/bash
#====================================================================
# The purpose of this program is processing pictures by ImageMagick
# The introduction of this software can be found in this website:
# https://doc.xuwenliang.com/docs/imagemagick/94
#========================================================================

fig_path=/home/ys17-19/renql/project/CONV_MON_WRK-2018/fig/animation/
#in_fig="$fig_path"preci_map_-60-60.0000
in_fig="$fig_path"wind_map_-90-90.0000

# 1. remove the withe edges
for((i=1;i<=9;i++))
do
    in_fig_name="$in_fig"0"$i".png
    echo $in_fig_name
    convert $in_fig_name -bordercolor white -trim $in_fig_name
done

for((i=10;i<=12;i++))
do
    in_fig_name="$in_fig""$i".png
    echo $in_fig_name
    convert $in_fig_name -bordercolor white -trim $in_fig_name
done

