#!/bin/bash

chara="pre_case = (/\"(a\",\"(b\"/)"
sed -i "17s/.*/${chara}/" 1902-draw_1p2x4_month_wind_w-shad.ncl
sed -i "s/smth49/smth0/g" 1902-draw_1p2x4_month_wind_w-shad.ncl
ncl 1902-draw_1p2x4_month_wind_w-shad.ncl

chara="pre_case = (/\"(c\",\"(d\"/)"
sed -i "17s/.*/${chara}/" 1902-draw_1p2x4_month_wind_w-shad.ncl
sed -i "s/smth0/smth25/g" 1902-draw_1p2x4_month_wind_w-shad.ncl
ncl 1902-draw_1p2x4_month_wind_w-shad.ncl

chara="pre_case = (/\"(e\",\"(f\"/)"
sed -i "17s/.*/${chara}/" 1902-draw_1p2x4_month_wind_w-shad.ncl
sed -i "s/smth25/smth49/g" 1902-draw_1p2x4_month_wind_w-shad.ncl
ncl 1902-draw_1p2x4_month_wind_w-shad.ncl

