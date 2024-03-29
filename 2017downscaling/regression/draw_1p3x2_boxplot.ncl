;*********************************************************************
;filename:      draw_boxplot.ncl
;input_file:   area_avrg_ts_winter.nc
;According rainfall to divide ts(var,time) into five levels which are
;(/"sunny","light_rain","moderate_rain","heavy_rain","rainstorm"/),
;then rank data from small to large and look for quartiles
;next draw the boxplot
;                                            by Ql_Ren
;                                           2018/03/14
;*********************************************************************************
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/shea_util.ncl"

begin
path      = "/home/ys17-19/renql/project/2017downscaling"
filename  = path + "/regression/mdata/" + "area_avrg_ts_winter.nc"
fig_name  = path + "/regression/fig/" + "1p3x2_boxplot_winter_delay_flux_divg"
fig_title  = "Relationship between rainfall and meteorological elements"
fig_out    = "pdf";"x11"
var_name   = "ts"
draw_var   = (/"700omega","850divg","300divg","850vort","flux"   ,"delay flux-divg"/) ;the first one of ts is rain
;draw_var   = (/"700omega","850divg","300divg","850vort","flux"   ,"flux-divg"/) ;the first one of ts is rain
magn       = (/    100   ,1000000  ,1000000  ,1000000  ,   1     , 10000     /)
MAG    = "10~S~" + (/"-2","-6"     ,"-6"     ,"-6"     ,"0"      ,"-4"       /) + "~N~"
unit       = (/  "Pa/s"  ,  "/s"   ,  "/s"   ,  "/s"   ,"g/(m*m)","g/(s*m*m)"/) ;the unit of rain is mm
rain_lvl   = (/-10,0,10,25,50,150/)
lvl_name   = (/"sunny","light","moderate","heavy","storm"/)

nvar  = dimsizes(draw_var)
nrain = dimsizes(lvl_name)
y     = new((/nvar ,nrain ,5/),float) ;the y axis
x     = (/1,2,3,4,5/)   ;the x axis
plot  = new((/nvar/),graphic)
mark  = new((/nvar,nrain/),graphic)

;==========================================================================
;read the ts(var,time) in winter and look for quartiles 
;===========================================================================
f  = addfile(filename,"r")
ts = f->$var_name$ ;get ts(var,time)

do nv = 1,nvar, 1
    ts(nv,:) = ts(nv,:)*magn(nv-1)
end do

do nr = 0, nrain-1, 1
    select_time := ind(ts(0,:).gt.rain_lvl(nr).and.ts(0,:).le.rain_lvl(nr+1))
    printVarSummary(select_time)

    dimt:= dimsizes(select_time)
    print(dimt)

    x25  = toint(floor(0.25*(dimt+1)))
    xx25 = 0.25*(dimt+1)-x25
    x75  = toint(floor(0.75*(dimt+1)))
    xx75 = 0.75*(dimt+1)-x75

    do nv = 0, nvar-1, 1
        if(nv.eq.5) then
        select_time = select_time -1
        select_time = where(select_time.lt.0,0,select_time)
        end if
        
        var := ts(nv+1,select_time)
        qsort(var)

        y(nv,nr,1) = (var(x25)-var(x25-1))*xx25+var(x25-1)
        y(nv,nr,2) = dim_median(var)
        y(nv,nr,3) = (var(x75)-var(x75-1))*xx75+var(x75-1)
        y(nv,nr,0) = min(var)
        y(nv,nr,4) = max(var)
    end do

end do
printVarSummary(y)

;==========================================================================
;draw the box plot
;===========================================================================
;define the thickness and height of the title,tick and label
   tick_font      = 0.02
   subtitle_font  = 0.02
   title_font     = 0.03
   label_font     = 0.015
   tick_thick     = 5.0
   subtitle_thick = 5.0  ;not have this attribute
   title_thick    = 5.0
   label_thick    = 1.5
   mp_thick       = 5.0

    wks = gsn_open_wks(fig_out,fig_name)

;-------set the attribute of xy-plot res----------
    res                  = True
    res@gsnMaximize      = True
    res@gsnPaperOrientation = "landscape" ;"portrait"
    res@tiXAxisOn           = False
    res@tiYAxisOn           = False
    res@tiXAxisFontHeightF  = label_font
    res@tiXAxisFontThicknessF = 1
    res@tiXAxisFontColor    = "black"
    res@tmXBLabels  = lvl_name

;---------set the attribute of the panel------------
    resp  = True    
    resp@gsnMaximize   = True
    resp@gsnPaperOrientation = "landscape" ;"portrait"
    ;resp@txString      = fig_title
    ;resp@txFontHeightF = title_font 
    ;resp@txFontThicknessF = title_thick

do nv = 0, nvar-1, 1
    res@tiMainString     = draw_var(nv) + " ( " + MAG(nv) + " " + unit(nv) + " )"
    plot(nv) = boxplot(wks,x,y(nv,:,:),False,res,False)
end do
    gsn_panel(wks,plot,(/2,3/),resp)
end
