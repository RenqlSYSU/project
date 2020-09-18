;***************************************************************************
;filename:      draw_12p3X1_27y_ts.ncl
;input_data:    Areas_preci_27y_ts.nc 
;read the seasonly U,V,T time series at three levs
;and then draw the xy plot,will have 12 panels. 
;                                            by Ql_Ren
;                                           2018/01/08
;*****************************************************************************
begin
;define the thickness and height of the title,tick and label
   tick_font      = 0.02
   subtitle_font  = 0.02
   title_font     = 0.03
   label_font     = 0.005
   tick_thick     = 5.0
   subtitle_thick = 5.0  ;not have this attribute
   title_thick    = 5.0
   label_thick    = 1.5
   mp_thick       = 5.0

;the path of the file
path   = "/home/ys17-19/renql/project/TP_NUDG/observation-20190127/"
case   = "ERA-Interim"
var_name = (/"uppEKE-2.5-10day","uppEKE-10-90day","uppEKE-2.5-10day","uppEKE-10-90day","hgt","uwnd","intQ1"/)   ;,"Q1"   the left title 
unit     = (/"%"               ,"%"              ,"10~S~4~N~J"      ,"10~S~4~N~J"     ,"m"  ,"m/s" ,"10W/m2"/)   ;,"W/kg"
draw_var = (/"STEhgh"          ,"STElow"         ,"TEhgh"           ,"TElow"          ,"NBR","Jet" ,"intQ1"/)   ;North Ridge; TPH
stand = False;True
nv = 0

;area = all, TP, up ,NBR1,NBR2
lats = (/25 ,25 ,25 ,35 ,35 /);-20
latn = (/50 ,40 ,40 ,50 ,50 /);60
lonl = (/0  ,70 ,35 ,70 ,0  /);20
lonr = (/360,105,70 ,95 ,360/);220

filename = path + "mdata/" + case + "_Clim_month_ave_" + var_name + ".nc"
if(stand) then 
fig_name = path + "fig/STE_indices_annual_ts_stand"
else
fig_name = path + "fig/STE_indices_annual_ts"
end if
fig_out  = "pdf"

nvar = dimsizes(var_name)
ts   = new((/nvar,12/),float)
plot = new((/nvar/),graphic)
;==============================================
;read data var(nvar,ncase,nyear)
;===================================================
do nv = 0,1,1
f    = addfile(filename(nv),"r")
term = f->month_ave(:,{lats(0):latn(0)},{lonl(0):lonr(0)})
term = term/1000.0
ts(nv,:) =  wgt_areaave(term(:,{lats(1):latn(1)},{lonl(1):lonr(1)}),1.0,1.0,0)\ 
          - wgt_areaave(term(:,{lats(2):latn(2)},{lonl(2):lonr(2)}),1.0,1.0,0)
ts(nv,:) = 100*ts(nv,:)/wgt_areaave(term(:,{lats(2):latn(2)},{lonl(2):lonr(2)}),1.0,1.0,0)
printMinMax(ts(nv,:),0)
end do

do nv = 2,3,1
f    = addfile(filename(nv),"r")
term := f->month_ave(:,{lats(1):latn(1)},{lonl(1):lonr(1)})
term = term/10000.0
ts(nv,:) =  wgt_areaave(term(:,{lats(1):latn(1)},{lonl(1):lonr(1)}),1.0,1.0,0)
printMinMax(ts(nv,:),0)
end do

nv = 4 ;NBR
f     = addfile(filename(nv),"r")
term := f->month_ave(:,{500},{lats(0):latn(0)},{lonl(0):lonr(0)})
term = term/10.0
;ts(nv,:) =  wgt_areaave(term(:,{lats(3):latn(3)},{lonl(3):lonr(3)}),1.0,1.0,0)
ts(nv,:) =  wgt_areaave(term(:,{lats(3):latn(3)},{lonl(3):lonr(3)}),1.0,1.0,0) \
          - wgt_areaave(term(:,{lats(4):latn(4)},{lonl(4):lonr(4)}),1.0,1.0,0)
printMinMax(ts(nv,:),0)

nv = 5 ;jet
f     = addfile(filename(nv),"r")
term := f->month_ave(:,{200},{lats(1):latn(1)},{lonl(1):lonr(1)})
ts(nv,:) = wgt_areaave(term,1.0,1.0,0) 
printMinMax(ts(nv,:),0)

nv = 6;intQ1
f     = addfile(filename(nv),"r")
term := f->month_ave(:,{lats(1):latn(1)},{lonl(1):lonr(1)})
term = term/10.0
ts(nv,:) = wgt_areaave(term,1.0,1.0,0)
printMinMax(ts(nv,:),0)
delete(term)

if(stand) then 
ts = dim_standardize_n(ts,1,1)
end if
ave  = dim_avg_n_Wrap(ts,1) ;ave(nvar,ncase,nseason,nlev)
printVarSummary(ts)
;===============================================
;draw the figure
;=============================================
wks = gsn_open_wks(fig_out,fig_name)

    resxy     = True
    resxy@vpWidthF  = 0.3
    resxy@vpHeightF = 0.25
    resxy@gsnFrame      = False
    resxy@gsnDraw       = False
    resxy@gsnMaximize   = True
    resxy@gsnPaperOrientation = "portrait"
    
    resxy@xyDashPatterns    = (/0  ,2  ,0  ,2  , 0 , 0 , 0 /)
    resxy@xyLineThicknesses = (/4.0,4.0,4.0,4.0,3.0,3.0,3.0/)   ; make second line thicker
    resxy@xyLineColors      = (/"black","black","lightcyan4","lightcyan4","red","blue","green"/) ; change line color
    resxy@tiXAxisOn = True 
    resxy@tiYAxisOn = False
    resxy@tiXAxisString = "month"
    ;resxy@trYMinF = -3.0
    ;resxy@trYMaxF = 15.0
    resxy@trXMinF    = 1
    resxy@trXMaxF    = 12
    resxy@tmXBMode   = "Explicit"
    resxy@tmXBValues = ispan(1,12,1)
    resxy@tmXBLabels = ispan(1,12,1)
    
    ;resxy@tmXBLabelFontHeightF     = tick_font-0.01
    ;resxy@tmXBLabelFontThicknessF  = tick_thick-1.0 
    ;resxy@tmYRLabelFontHeightF     = tick_font
    ;resxy@tmYRLabelFontThicknessF  = tick_thick

    resxy@tmXBMinorOn          = False ; Turn off x-axis (bottom) minor ticks
    resxy@tmXTMinorOn          = False ; Turn off x-axis (tottom)minor ticks
    resxy@tmXTOn               = False ; Turn off x-axis (top) major ticks
    resxy@tmYRMinorOn          = False ; Turn off y-axis (right) minor ticks
    resxy@tmYLMinorOn          = False ; Turn off y-axis (left) minor ticks
    resxy@tmYROn               = False ; Turn off y-axis (right) major ticks
    
    ;resxy@tmXBMajorLengthF = 0.01
    ;resxy@tmYRMajorLengthF = 0.01

    resxy@tmXMajorGrid                   = True
    resxy@tmXMajorGridLineDashPattern    = 2
    resxy@tmXMajorGridLineColor          = "gray"
    resxy@tmXMajorGridThicknessF         = 2.0
                       
    resxy@gsnYRefLine             = 0.0 
    resxy@gsnYRefLineDashPatterns = 2
    resxy@gsnYRefLineThicknesses  = 2.0
    resxy@gsnYRefLineColors       = "gray"
                             
    ;legend
    resxy@pmLegendDisplayMode = "Always"
    resxy@lgOrientation       = "Vertical"
    resxy@lgPerimOn           = False
    resxy@lgPerimFill         = "SolidFill"
    resxy@lgPerimFillColor    = "white"
    resxy@xyExplicitLegendLabels = draw_var + " (" + unit + ")" 
    resxy@lgLabelFontColor       = "black"
    resxy@lgLabelFontThicknessF  = label_thick
    resxy@lgLabelFontHeightF     = label_font
    resxy@pmLegendWidthF   = 0.04
    resxy@pmLegendHeightF  = 0.07
    resxy@pmLegendZone  = 0   ;set the positional origin in the center of the plot
    resxy@pmLegendOrthogonalPosF  = 0.35   ; move ref vector along Y-axis
    resxy@pmLegendParallelPosF    = 0.20    ; move ref vector along X-axis

    resp  = True    
    resp@gsnMaximize   = True
    resp@gsnPaperOrientation = "landscape";"portrait"
    ;resp@gsnPanelLabelBar = True
    ;resp@pmLabelBarWidthF  = 0.5
    ;resp@pmLabelBarHeightF = 0.05
    ;resp@lbLabelFontHeightF = tick_font
    ;resp@lbLabelFontThicknessF = tick_thick
    resp@txFontHeightF = title_font 
    resp@txFontThicknessF = title_thick
;   
;    resxy@gsnCenterString  = vari(nv,0,ns,nl)+"  NUDG6h:"+vari(nv,1,ns,nl)+"  NUDG24h:"+vari(nv,2,ns,nl)
    resxy@gsnLeftString    = case 
;    resxy@gsnRightString   = var_name(nv) + " (" + unit(nv) + ")"
;    if(stand) then 
;    resxy@gsnYRefLine      = (/-1.5,0.0,1.5/)
;    else
;    resxy@gsnYRefLine      = (/ave(0),ave(1),ave(2)/)
;    end if
    plot = gsn_csm_xy(wks,ispan(1,12,1),(/ts/),resxy)

    gsn_panel(wks,plot,(/1,1/),resp)
end

