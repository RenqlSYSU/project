;***************************************************************************
;filename:      draw_1p3x3_UVT_annual_series_MethodologyPaper.ncl
;plot the every year's annual cycle(356) of TP average UVT(27 curves) 
;then overlaying the climate annual cycle(black bold curve)
;and one panel has nine plots
;this fig is used to the paper.
;                                            by Ql_Ren
;                                           2017/04/11
;******************************************************************************
begin
path     = "/home/ys17-19/renql"
draw_var = (/"U","V","T"/)
unit     = (/"m/s","m/s","~S~o~N~C"/)   ;this can be used as right string 
pre_case = (/"(a","(b","(c", "(d","(e","(f", "(g","(h","(i"/)
nv = 2

path    = "/home/ys17-19/renql/project/TP_NUDG/concept_map"
filein  = path + "/mdata/F2000_CAM5.cam.h1.ESM.clim." + draw_var + ".nc"
gtopfile = "/home/ys17-19/renql/project/TP_NUDG/analysis/mdata/gtopo30_0.9x1.25.nc"
case    = (/"input2CESM"/)
fig_name = path + "/fig/" + case(0) + "_annual_cycle"
fig_out  = "pdf" ;"eps" ;"x,11"

;area of TP 
lats = (/25 /)
latn = (/40 /)
lonl = (/70 /)
lonr = (/105/)
lev  = (/850,500,300,200/)
na = 0

;define the thickness and height of the title,tick and label
   tick_font      = 0.01
   subtitle_font  = 0.013
   title_font     = 0.03
   label_font     = 0.015
   tick_thick     = 5.0
   subtitle_thick = 5.0  ;not have this attribute
   title_thick    = 5.0
   label_thick    = 1.5
   mp_thick       = 4.0
   font_index     = 22

;==========================================================
;create an array(nvar,ncase,nyear) to store the data
;============================================================
f = addfile(gtopfile,"r")
phis   = f->PHIS({lats(na):latn(na)},{lonl(na):lonr(na)})
phis   = phis/9.8   ;convert unit from m2/s2 to gpm
ncase = dimsizes(case)
nvar  = dimsizes(draw_var)
nlev  = dimsizes(lev)
nday  = 350
plot  = new((/ncase*nlev/),graphic)
plot2 = new((/ncase*nlev/),graphic)
plot3 = new((/ncase*nlev/),graphic)
                     
;---------------------------------------------------------------
;read the data, then calculate the area and season average
;--------------------------------------------------------------
wks = gsn_open_wks(fig_out,fig_name)
    gsn_define_colormap(wks,"matlab_jet") 
    ;has 64 colors blue green yellow,red

do nv = 0,nvar-1,1
f   = addfile(filein(nv),"r")
var = f->$draw_var(nv)$(0:(nday-1),{lev},{lats(na):latn(na)},{lonl(na):lonr(na)})
var&time = ispan(1,350,1)

if(nv.eq.2) then
var = var - 273.15
end if

printVarSummary(var)
nlat = dimsizes(var&lat)
nlon = dimsizes(var&lon)
ts = reshape(var(lev|:,lat|:,lon|:,time|:),(/nlev,nlat*nlon,nday/))

;===============================================
;draw the figure
;=============================================
    resxy     = True
    resxy@vpWidthF  = 0.3
    resxy@vpHeightF = 0.1
    resxy@gsnFrame      = False
    resxy@gsnDraw       = False
    resxy@gsnMaximize   = True
    resxy@gsnPaperOrientation = "portrait"
    
    resxy@xyDashPatterns    = 0 ;(/0,0,0/)
    resxy@xyLineThicknesses = 2.0;(/2.0, 2.0,2.0/)   ; make second line thicker
    resxy@xyLineColors      = ispan(4,62,2)
    resxy@tiXAxisOn = False
    resxy@tiYAxisOn = False
    ;resxy@trYMinF = floor(min(var))
    ;resxy@trYMaxF = ceil(max(var))
    
    resxy@tmBorderThicknessF  = mp_thick
    resxy@tmXBLabelFont         = font_index
    resxy@tmYLLabelFont         = font_index
    resxy@tmXBLabelFontHeightF     = tick_font
    resxy@tmXBLabelFontThicknessF  = tick_thick 
    resxy@tmYLLabelFontHeightF     = tick_font
    resxy@tmYLLabelFontThicknessF  = tick_thick

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
    ;resxy@tmXBMode           = "Manual"
    ;resxy@tmXBTickSpacingF   = 30
    ;resxy@tmXBTickStartF     = 0
    ;resxy@tmXBTickEndF       = 360
    resxy@tmXBMode            = "Explicit"
    resxy@tmXBValues          = ispan(0,360,30)
    resxy@tmXBLabels          = ispan(0,360,30)
    ;resxy@tmXBLabels   = (/"Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec"/)
    
    resxy@tmXMajorGrid                   = True
    resxy@tmXMajorGridLineDashPattern    = 2
    resxy@tmXMajorGridLineColor          = "gray"
    resxy@tmXMajorGridThicknessF         = 2.0
                             
    resxy@gsnStringFont         = font_index
    resxy@gsnStringFontHeightF   = subtitle_font
    resxy@gsnCenterString  = "";lev(nl) + " hPa" 
    resxy@gsnLeftString    = ""
    resxy@gsnRightString   = ""

    resxy2 = resxy
    resxy2@xyMonoLineColor  = True
    resxy2@xyMonoLineThickness = True
    resxy2@xyLineThicknessF = 4.0;(/2.0, 2.0,2.0/)   ; make second line thicker
    resxy2@xyLineColor      = "black"

    resxy3 = resxy2
    resxy2@xyLineThicknessF = 6.0;(/2.0, 2.0,2.0/)   ; make second line thicker
    
    resp  = True    
    resp@gsnMaximize   = True
    resp@gsnPaperOrientation = "landscape"

do nl = 0,nlev-1,1
    np = nl
    resxy@gsnLeftString   = pre_case(np) + ") " + case
    resxy@gsnRightString  = lev(nl) + "hPa " + draw_var(nv) + " (" + unit(nv) + ")"  
    resxy@trYMinF = floor(min(var(:,nl,:,:)))
    resxy@trYMaxF =  ceil(max(var(:,nl,:,:)))
    plot(np)  = gsn_csm_xy(wks,var&time,(/ts(nl,::10,:)/),resxy)
end do
    
   gsn_panel(wks,plot,(/nlev,ncase/),resp)
end do
end

