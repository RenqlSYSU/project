;******************************************************************************
;filename:      draw_1p4x3_month_sst_EOFts.ncl
;input_file:    var_name(nv) + area + _EOFn3.nc
;----------------------------------------------------------------------
; Concepts illustrated:
;   - read EOF time series data
;   - Draw the first three EOF principal components
;   - each panel have four time series figure(have one curves)
;                                            by Ql_Ren
;                                           2018/12/04
;*********************************************************************************
begin
unit     = (/"~S~o~N~C","hPa","W/m2","W/m2"/)    ;this can be used as right string 
var_name = (/"TS","PSL","FLNT","FSNT"/)
case     = (/"MX_AQUAP_SCYC","SF_AQUAP_SCYC","NO_TOPO_SCYC","REAL_WORLD_SCYC"/);" NG24H",THIS CAN BE USED AS LEFT STRING
pre_case = (/"(a","(b","(c","(d"/)
nv = 0

;-----EOF area-------------- 
lats = -30
latn = 30
lonl = 90
lonr = 300
;lev  = 200
neval = 4
lat_spacing = 15
lon_spacing = 20

opt = True
opt@jopt = 0 ;0 mean use the covariance matrix to compute EOFs (default), 1 mean use the correlation matrix
eof_kind = (/"cov","corr"/)
tname    = (/"_","_DJF_"/)

;----path of the file and output file-----------------------------
path      = "/home/ys17-19/renql/project/CONV_MON_WRK-2018" 
filename  = path + "/mdata/" + var_name(nv) + lats + "S-" + latn + "N_" + lonl + "E-" + lonr + "E" + tname(0) + eof_kind(opt@jopt) + "EOF" + neval + ".nc"
fig_name  = path + "/fig/"   + var_name(nv) + lats + "S-" + latn + "N_" + lonl + "E-" + lonr + "E" + tname(0) + eof_kind(opt@jopt) + "EOFts" + neval
fig_out   = "pdf";"x11"

;================================================================
;read data
;====================================================================
f = addfile(filename,"r")
ts_var  = dim_standardize_n_Wrap(f->eof_ts,1,2)    ;(ncase,neval,nyear) 
right_title   = f->pcvar
center_title  = f->eval
printVarSummary(ts_var)
printVarSummary(right_title)

ncase = dimsizes(ts_var&case)
ntime = dimsizes(ts_var&time)
sig   = new((/ncase,neval/),logical) 
plot  = new((/ncase*(neval-1)/),graphic)

do nc = 0,ncase-1,1
sig(nc,:) = eofunc_north(right_title(nc,:),ntime,True)
end do
subtitle_color = where(sig,"red","black")

;==============================================================
;draw figure 
;==================================================================
;define the thickness and height of the title,tick and label
   tick_font      = 0.02
   subtitle_font  = 0.02
   title_font     = 0.02
   label_font     = 0.02
   tick_thick     = 5.0
   subtitle_thick = 5.0  ;not have this attribute
   title_thick    = 5.0
   label_thick    = 5.0
   mp_thick       = 3.0
   font_index     = 22

    wks = gsn_open_wks(fig_out,fig_name)
    gsn_define_colormap(wks,"GMT_polar") 
    ;blue and red, has 20 colors, 2-11 is bule, 12-21 is red
    
;------------set the resxy for time series figure-----------
    resxy     = True
    ;resxy@vpXF      = 0.1
    ;resxy@vpYF      = 0.4
    resxy@vpWidthF  = 0.9
    resxy@vpHeightF = 0.2
    resxy@gsnFrame      = False
    resxy@gsnDraw       = False
    resxy@gsnMaximize   = True
    
    resxy@xyDashPatterns     = (/0/)
    resxy@xyLineThicknesses = (/3.0/)   ; make second line thicker
    resxy@xyLineColors      = (/"black"/) ; change line color
    resxy@tiXAxisOn = False
    resxy@tiYAxisOn = False
    resxy@trYMinF = -3.0
    resxy@trYMaxF = 3.0
    
    resxy@tmBorderThicknessF  = mp_thick
    resxy@tmXBLabelFont         = font_index
    resxy@tmYLLabelFont         = font_index
    resxy@tmXBLabelFontHeightF     = label_font
    resxy@tmXBLabelFontThicknessF  = label_thick 
    resxy@tmYRLabelFontHeightF     = label_font
    resxy@tmYRLabelFontThicknessF  = label_thick

    resxy@tmXBMinorOn          = False ; Turn off x-axis (bottom) minor ticks
    resxy@tmXTMinorOn          = False ; Turn off x-axis (tottom)minor ticks
    resxy@tmXTOn               = False ; Turn off x-axis (top) major ticks
    resxy@tmYRMinorOn          = False ; Turn off y-axis (right) minor ticks
    resxy@tmYLMinorOn          = False ; Turn off y-axis (left) minor ticks
    resxy@tmYROn               = False ; Turn off y-axis (right) major ticks
    
    resxy@tmXBMajorThicknessF     = mp_thick
    resxy@tmYLMajorThicknessF     = mp_thick
    resxy@tmXBMinorThicknessF     = mp_thick
    resxy@tmYLMinorThicknessF     = mp_thick
    ;resxy@tmXBMajorLengthF = 0.01
    ;resxy@tmYRMajorLengthF = 0.01
    resxy@tmXBMode           = "Manual"
    resxy@tmYLMode           = "Manual"
    resxy@tmXBTickSpacingF   = ntime/10
    ;resxy@tmYLTickSpacingF   = 1.0

    resxy@tmXMajorGrid                   = True
    resxy@tmXMajorGridLineDashPattern    = 2
    resxy@tmXMajorGridLineColor          = "gray"
    resxy@tmXMajorGridThicknessF         = 2.0
                       
    resxy@gsnYRefLine             = 0
    resxy@gsnYRefLineDashPatterns = 0 ; (/12,12/)   
    resxy@gsnYRefLineThicknesses  = 3.0 ;(/3.0,3.0/)
    resxy@gsnYRefLineColors       = "black"  ;(/"blue","red"/)
                             
    resxy@gsnStringFont        = font_index
    resxy@gsnStringFontHeightF = subtitle_font
    
;------------------------------------------------------------------
;draw the subgraph and change the left and right titles
;-----------------------------------------------------------------------
    resp  = True    
    resp@gsnMaximize   = True
    resp@txFont = font_index
    resp@txFontHeightF = title_font 
    resp@txFontThicknessF = title_thick
    ;resp@gsnPanelRowSpec  = True
    ;resp@gsnPanelBottom   = 0.2
    resp@gsnPaperOrientation = "landscape" ;"portrait"

do nc = 0,ncase-1,1
do ne = 0,neval-2,1
    np = nc*(neval-1) + ne
    resxy@gsnLeftString    = case(nc)
    resxy@gsnCenterString  = "EOF " + (ne+1)
    resxy@gsnRightStringFontColor    = subtitle_color(nc,ne)
    resxy@gsnRightString   = center_title(nc,ne)+"  "+right_title(nc,ne)+"%"
    plot(np) = gsn_csm_xy(wks,ts_var&time,(/ts_var(nc,ne,:)/),resxy)
    print(max(ts_var(nc,ne,:)))
    print(min(ts_var(nc,ne,:)))
end do
end do 
    
    resp@txString      = "EOF of " + var_name(nv)
    gsn_panel(wks,plot,(/ncase,(neval-1)/),resp)
end

