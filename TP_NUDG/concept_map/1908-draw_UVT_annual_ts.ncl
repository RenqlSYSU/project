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
case_option = 2 ;0 means NUDG experiment; 1 means noTP experiment
pre_case = (/"(a","(b","(c", "(d","(e","(f", "(g","(h","(i"/)
nv = 2

;the path of the file
if(case_option.eq.0) then
    case    = (/"CTRL","NG6h","NG24h"/)
    filein  = path +(/"/model/AMIP-CTRL/AMIP_C5PM.cam.h1.YEAR.",\
                      "/model/TP-NUDG-6h/AMIP_C5PM_TP_NUDG.cam.h1.YEAR.",\
                      "/model/TP-NUDG-24h/AMIP_C5PM_TP_NUDG.cam.h1.YEAR."/)
    years    = ispan(1979,2005,1)   ;the calculate time
    gtopfile= path + "/project/TP_NUDG/analysis/mdata/" + "gtopo30_0.9x1.25.nc"
end if  
if(case_option.eq.1) then
    ;remove TP experiment, sample number is ten years 
    case    = (/"TP_CTRL","TP_CR"/) ;TP_CTRL is 0021-0030; TP_CR is 0041-0050
    filein  = path + "/model/" + case + "/pro/" + case + ".cam.h1.YEAR.00" ; + year + ".daily." + var_name + ".nc"
    ;gtopfile = "/home/ys17-19/renql/model/TP_CR/USGS-gtopo30_1.9x2.5_noTP_c050602C.nc" 
    gtopfile = "/home/ys17-19/renql/model/TP_CTRL/USGS-gtopo30_1.9x2.5_remap_c050602.nc" 
    years   = ispan(21,30,1)   ;the calculate time
end if
if(case_option.eq.2) then
    case    = (/"F2000","NGTR"/) ;
    filein  = path +(/"/model/F2000_CAM5/F2000_CAM5.cam.h1.YEAR.00",\
                      "/model/F2000_CAM5_TP_NUDG/F2000_CAM5_TP_NUDG.cam.h1.YEAR.00"/)
                      ;"/model/F2000_CAM5_NUDG/F2000_CAM5.cam.h1.YEAR.00"/)
    gtopfile = path + "/project/TP_NUDG/analysis/mdata/" + "gtopo30_0.9x1.25.nc"
    years   = ispan(1,30,1)   ;the calculate time
end if
fig_name = path + "/project/TP_NUDG/concept_map/fig/" + case(0) + "_annual_cycle-" + draw_var(nv) 
;fig_out  = "eps" ;"eps" ;"x,11"
fig_out  = "pdf" ;"eps" ;"x,11"

;area of TP 
lats = (/25 /)
latn = (/40 /)
lonl = (/70 /)
lonr = (/105/)
lev  = (/300,200/);850,500,,200
na = 0

;define the thickness and height of the title,tick and label
   tick_font      = 0.01
   subtitle_font  = 0.012
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
nyear = dimsizes(years)
ncase = dimsizes(case)
nvar  = dimsizes(draw_var)
nlev  = dimsizes(lev)
nday1 = 365
nday2 = 350
var   = new((/ncase,nyear,nday1,nlev/),float)
ave   = new((/ncase,      nday1,nlev/),float)
plot  = new((/ncase*nlev/),graphic)
plot2 = new((/ncase*nlev/),graphic)
plot3 = new((/ncase*nlev/),graphic)
                     
;---------------------------------------------------------------
;read the data, then calculate the area and season average
;--------------------------------------------------------------
wks = gsn_open_wks(fig_out,fig_name)
    gsn_define_colormap(wks,"matlab_jet") 
    ;has 64 colors blue green yellow,red

do nv = 2,2,1;nvar-1,1
do nc = 0,ncase-1,1
files = systemfunc("ls " + filein(nc) + "*.daily."+draw_var(nv)+".nc")
do nt = 0, nyear-1 ,1
    print("handle with " + files(nt))
    f   = addfile(files(nt),"r")
    temp = f->$draw_var(nv)$(:,{lev},{lats(na):latn(na)},{lonl(na):lonr(na)})
    temp = mask(temp,phis.gt.3000,True) 
    var(nc,nt,:,:) = wgt_areaave_Wrap(temp,1.0,1.0,0) ;get var(nvar,ncase,nyear,365,nlev)
end do
end do

if(case_option.eq.2) then 
    cf = ezfftf_n(dim_avg_n(var,1),1)
    cf(:,:,4:,:)  = 0
    ave = ezfftb_n(cf,cf@xbar,2)
else
    ave = dim_avg_n_Wrap(var,1)
end if

if(nv.eq.2) then
var = var - 273.15   ;convert unit from k to c
ave = ave - 273.15
end if

var1 = var(:,:,0:nday2-1,:)
ave1 = ave(:,  0:nday2-1,:)
var1&time = ispan(1,nday2,1)

printVarSummary(var)
printVarSummary(ave)

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
    resxy@tiXAxisOn = True;False
    resxy@tiYAxisOn = True ;False
    resxy@tiXAxisString = "day"
    resxy@tiYAxisString = draw_var(nv) + "(" + unit(nv) + ")" 
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
    ;resxy@tmXBLabels          = ispan(0,360,30)
    resxy@tmXBLabels   = (/"Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec"/)
    
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

do nl = 0,0,1;nlev-1,1
do nc = 0, ncase-1,1
    np = ncase*nl + nc
;    np = nlev*nc + nl
    resxy@gsnLeftString   = pre_case(np) + ") " + case(nc)
    resxy@gsnRightString  = "" ;lev(nl) + "hPa " + draw_var(nv) + " (" + unit(nv) + ")"  
    resxy@trYMinF = floor(min(var1(:,:,:,nl)))
    resxy@trYMaxF =  ceil(max(var1(:,:,:,nl)))
    plot(np)  = gsn_csm_xy(wks,var1&time,(/var1(nc,:,:,nl)/),resxy)
    plot2(np) = gsn_csm_xy(wks,var1&time,(/ave1(0 ,  :,nl)/),resxy2)
    ;plot3(np) = gsn_csm_xy(wks,var&day,(/ave(nc,nv,  :,nl)/),resxy3)
    overlay(plot(np),plot2(np))
    ;overlay(plot(np),plot3(np))
end do
end do
    
   gsn_panel(wks,plot,(/nlev,ncase/),resp)
end do
end

