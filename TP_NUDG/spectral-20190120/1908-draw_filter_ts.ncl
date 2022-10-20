;***************************************************************************
;filename:      spec.ncl
;read the 27 years U,V,T at TP areas in three cases (CTRL, NUDG6h, NUDG24h)
;then do the spectral analysis and plot the figure
;                                           by Ql_Ren
;                                           2019/01/20
;******************************************************************************
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/shea_util.ncl"

begin
path    = "/home/ys17-19/renql"
var_name = (/"U","V","T"/)
case_option = 0 ;0 means NUDG experiment; 1 means noTP experiment

if(case_option.eq.0) then
    case    = (/"CTRL","NG6h","NG24h"/)
    filein  = path +(/"/model/AMIP-CTRL/AMIP_C5PM.cam.h1.YEAR.",\
                      "/model/TP-NUDG-6h/AMIP_C5PM_TP_NUDG.cam.h1.YEAR.",\
                      "/model/TP-NUDG-24h/AMIP_C5PM_TP_NUDG.cam.h1.YEAR."/)
    years   = ispan(1979,1980,1)   ;the calculate time
else
    ;remove TP experiment, sample number is ten years 
    case1   = (/"TP_CTRL","TP_CR"/) ;TP_CTRL is 0021-0030; TP_CR is 0041-0050
    case    = (/"CTRL","noTP"/)
    filein  = path + "/model/" + case1 + "/pro/" + case1 + ".cam.h1.YEAR.00" ; + year + ".daily." + var_name + ".nc"
    years   = ispan(21,30,1)   ;the calculate time
end if

;area of TP 
lats = 25 
latn = 40
lonl = 70
lonr = 105
lev  = (/850,500,200/)

nv = 0
nl = 2
nc = 0

fig_name  = path + "/project/TP_NUDG/spectral-20190120/fig/" + case(nc) + "filter_ts_" + var_name(nv) + "_" + lev(nl)
;fig_title = case(nc) + " " + var_name(nv) + " " + lev(nl)

f   = addfile(filein(0)+years(0)+".daily."+var_name(nv)+".nc","r")
time        = cd_calendar(f->time,0)
select_time = ind(time(:,1).ge.1.and.time(:,1).le.12)
ntime       = dimsizes(select_time)
ncase = dimsizes(case)
nlev  = dimsizes(lev)
nyear = dimsizes(years)

;define array to calc the average of 27 years
var  = new((/nyear,ntime/),float)
ts   = new((/7,ntime/),float)
;plot = new((/ncase*nlev/),graphic)

;---------------------------------------------------------------
;read the data, then calculate the area and season average
;--------------------------------------------------------------
do nt = 0, nyear-1 ,1
    f   = addfile(filein(nc)+years(nt)+".daily."+var_name(nv)+".nc","r")
    var1 = f->$var_name(nv)$(select_time,{lev(nl)},{lats:latn},{lonl:lonr})
    var(nt,:) = wgt_areaave_Wrap(var1,1.0,1.0,0) ;get var(nvar,ncase,nyear,365,nlev)
end do

snt = 1
ts(0,:) = var(snt,:)
ts(1,:) = dim_avg_n(var,0)

cf = ezfftf(ts(1,:))
cf(:,4:) = 0
ts(2,:) = ezfftb(cf,cf@xbar)
printVarSummary(cf)

ca1  = 2.5
ca2  = 10
ca3  = 90
fca1 = 1.0/ca1
fca2 = 1.0/ca2
fca3 = 1.0/ca3
cf2  = ezfftf(ts(0,:))
fcf  = fspan(0, 0.5, ntime/2) 

ifcf := ind(fcf.gt.fca1.or.fcf.lt.fca2)
term = cf2 
term(:,ifcf) = 0.0  
;ts(3,:) = ezfftb(term,cf2@xbar)
ts(3,:) = ezfftb(term,0)

ifcf := ind(fcf.gt.fca2.or.fcf.lt.fca3)
term = cf2 
term(:,ifcf) = 0.0  
ts(4,:) = ezfftb(term,0)

opt = False
ts(5,:) = bw_bandpass_filter(ts(0,:),fca2,fca1,opt,0)
ts(6,:) = bw_bandpass_filter(ts(0,:),fca3,fca2,opt,0)

ts_label = (/tostring(years(snt)),"Clim","Filter_clim" \
            ,tostring(ca1)+"days-"+tostring(ca2)+"days_ez"\
            ,tostring(ca2)+"days-"+tostring(ca3)+"days_ez"\
            ,tostring(ca1)+"days-"+tostring(ca2)+"days_bw"\
            ,tostring(ca2)+"days-"+tostring(ca3)+"days_bw"/)

;if(nv.eq.2)
;    ts = ts - 273.15
;end if

;************************************************
; plotting
;************************************************
;define the thickness and height of the title,tick and label
   tick_font      = 0.015
   subtitle_font  = 0.018
   title_font     = 0.03
   label_font     = 0.015
   tick_thick     = 5.0
   subtitle_thick = 5.0  ;not have this attribute
   title_thick    = 5.0
   label_thick    = 1.5
   mp_thick       = 2.0
   font_index     = 22

wks  = gsn_open_wks("pdf",fig_name)             ; send graphics to PNG file

    resxy     = True
    resxy@vpWidthF  = 0.45
    resxy@vpHeightF = 0.15
    resxy@gsnFrame      = False
    resxy@gsnDraw       = False
    resxy@gsnMaximize   = True
    resxy@gsnPaperOrientation = "portrait"
   
    resxy@xyDashPatterns    = (/0,1,0,1/)
    resxy@xyLineThicknesses = (/2.5,2.5,2.5,2.5/)   ; make second line thicker
    resxy@xyLineColors      = (/"red","red","blue","blue"/)
    ;resxy@xyDashPatterns    = (/0,0/)
    ;resxy@xyLineThicknesses = (/2.5,2.5/)   ; make second line thicker
    ;resxy@xyLineColors      = (/"red","blue"/)
    resxy@tmYLAutoPrecision = False
    resxy@tmYLPrecision     = 2 
    ;resxy@tiXAxisOn = False
    ;resxy@tiYAxisOn = False
    resxy@tiXAxisString = "";"day"          ; xaxis
    resxy@tiYAxisString = "";"Variance/freq"           ; yaxis
    ;resxy@trYLog  = True
    ;resxy@trYMinF = 0
    ;resxy@trYMaxF = 500
    resxy@tmXBMode            = "Explicit"
    resxy@tmXBValues          = ispan(0,365,73)
    resxy@tmXBLabels          = ispan(0,365,73)
    
    resxy@gsnYRefLine            = 0.0 
    resxy@gsnYRefLineDashPatterns = (/0/) 
    resxy@gsnYRefLineThicknesses  = (/3.0/)
    resxy@gsnYRefLineColors       = (/"black"/)
    
    ;legend
    resxy@pmLegendDisplayMode = "Always"
    resxy@lgOrientation       = "Vertical"
    resxy@lgPerimOn           = False
    resxy@lgPerimFill         = "SolidFill"
    resxy@lgPerimFillColor    = "white"
    resxy@xyExplicitLegendLabels = ts_label(3:6)
    resxy@lgLabelFontColor       = "black"
    resxy@lgLabelFontThicknessF  = label_thick
    resxy@lgLabelFontHeightF     = label_font
    resxy@pmLegendWidthF   = 0.07
    resxy@pmLegendHeightF  = 0.15
    resxy@pmLegendZone  = 0   ;set the positional origin in the center of the plot
    resxy@pmLegendOrthogonalPosF  = 0   ; move ref vector along Y-axis
    resxy@pmLegendParallelPosF    = 0.7    ; move ref vector along X-axis
    
    resxy@tmBorderThicknessF  = mp_thick
    resxy@tmXBLabelFont         = font_index
    resxy@tmYLLabelFont         = font_index
    resxy@tmXBLabelFontHeightF     = tick_font 
    resxy@tmXBLabelFontThicknessF  = tick_thick
    resxy@tmYLLabelFontHeightF     = tick_font
    resxy@tmYLLabelFontThicknessF  = tick_thick

    ;resxy@tmXBMinorOn          = False ; Turn off x-axis (bottom) minor ticks
    ;resxy@tmXTMinorOn          = False ; Turn off x-axis (tottom)minor ticks
    ;resxy@tmXTOn               = False ; Turn off x-axis (top) major ticks
    ;resxy@tmYRMinorOn          = False ; Turn off y-axis (right) minor ticks
    ;resxy@tmYLMinorOn          = False ; Turn off y-axis (left) minor ticks
    ;resxy@tmYROn               = False ; Turn off y-axis (right) major ticks
    
    resxy@tmXBMajorThicknessF     = mp_thick
    resxy@tmYLMajorThicknessF     = mp_thick
    resxy@tmXBMinorThicknessF     = mp_thick
    resxy@tmYLMinorThicknessF     = mp_thick
    ;resxy@tmXBMajorLengthF = 0.01
    ;resxy@tmYRMajorLengthF = 0.01
    
    resxy@tmXMajorGrid                   = True
    resxy@tmXMajorGridLineDashPattern    = 2
    resxy@tmXMajorGridLineColor          = "gray"
    resxy@tmXMajorGridThicknessF         = 2.0
                             
    resxy@gsnStringFont         = font_index
    resxy@gsnStringFontHeightF = subtitle_font
    resxy@gsnCenterString  = "";lev(nl) + " hPa" 
    resxy@gsnLeftString    = case(nc)
    resxy@gsnRightString   = var_name(nv) + " " + lev(nl) + "hPa" 

   plot = gsn_csm_xy(wks,ispan(0,364,1),ts(3:6,:),resxy)     
    
    resp  = True    
    resp@gsnMaximize   = True
    resp@gsnPaperOrientation = "portrait" ;"landscape";
    ;resp@gsnPanelLabelBar = True
    ;resp@pmLabelBarWidthF  = 0.5
    ;resp@pmLabelBarHeightF = 0.05
    ;resp@lbLabelFontHeightF = tick_font
    ;resp@lbLabelFontThicknessF = tick_thick
    resp@txFontHeightF = title_font 
    resp@txFontThicknessF = title_thick

;    resp@txString      = fig_title
;    gsn_panel(wks,plot,(/nlev,ncase/),resp)

draw(wks)
frame(wks)
end
