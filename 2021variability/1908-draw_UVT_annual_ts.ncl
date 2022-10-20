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
draw_var = (/"preci" ,"U","V","T"/)
unit     = (/"mm/day","m/s","m/s","~S~o~N~C"/)   ;this can be used as right string 
case_option = 0 ;0 means NUDG experiment; 1 means noTP experiment
pre_case = (/"(a","(b","(c", "(d","(e","(f", "(g","(h","(i"/)
nv = 0
month_s = 5 
month_e = 9

;the path of the file
if(case_option.eq.0) then
    case    = (/"CTRL","NG15S","NG4550","NGEN"/) ;,"NUDG24h"
    years   = ispan(1979,2005,1)   ;the calculate time
    path2   = "/home/ys17-19/renql"
    filein  = path +(/"/model/AMIP-CTRL/AMIP_C5PM.cam.h1.YEAR.",\
                      "/model/AMIP_NG15S/AMIP_CAM5_NG15S.cam.h1.YEAR.",\
                      "/model/AMIP_NG4550/AMIP_CAM5_NG4550.cam.h1.YEAR.",\
                      "/model/AMIP_NGEN/AMIP_CAM5_NGEN.cam.h1.YEAR."/)
    gtopfile= path + "/project/TP_NUDG/analysis/mdata/" + "gtopo30_0.9x1.25.nc"
end if  
fig_name = path + "/project/2021variability/fig/" + case(0) + "_annual_cycle-" + draw_var(nv) 
fig_out  = "pdf" ;"eps" ;"x,11"

area = (/"EA","SEA","NEA"/); "TP","STP","CTP", ;S mean south,N mean North,EA mean East Asia
lats = (/ 27 , 22  , 30  /)
latn = (/ 37 , 30  , 45  /)
lonl = (/105 ,107  ,115  /)
lonr = (/125 ,130  ,130  /)
lev  = (/300,200/);850,500,,200
na = 0

;define the thickness and height of the title,tick and label
   tick_font      = 0.01
   subtitle_font  = 0.01
   title_font     = 0.03
   label_font     = 0.01
   tick_thick     = 5.0
   subtitle_thick = 5.0  ;not have this attribute
   title_thick    = 5.0
   label_thick    = 1.5
   mp_thick       = 2.0
   font_index     = 22

;==========================================================
;create an array(nvar,ncase,nyear) to store the data
;============================================================
f = addfile(gtopfile,"r")
phis   = f->PHIS({lats(na):latn(na)},{lonl(na):lonr(na)})
phis   = phis/9.8   ;convert unit from m2/s2 to gpm

f    = addfile(filein(0)+years(0)+".daily.PRECC.nc","r")
time = cd_calendar(f->time,0)
if(month_s.gt.month_e) then 
select_time = ind(time(:,1).ge.month_s.or.time(:,1).le.month_e)
else
select_time = ind(time(:,1).ge.month_s.and.time(:,1).le.month_e)
end if
nday  = dimsizes(select_time)
nyear = dimsizes(years)
ncase = dimsizes(case)
nvar  = dimsizes(draw_var)
nlev  = dimsizes(lev)
var   = new((/ncase,nyear,nday/),float)
ave   = new((/ncase,      nday/),float)
plot  = new((/ncase*nlev/),graphic)
plot2 = new((/ncase*nlev/),graphic)
plot3 = new((/ncase*nlev/),graphic)
                     
;---------------------------------------------------------------
;read the data, then calculate the area and season average
;--------------------------------------------------------------
wks = gsn_open_wks(fig_out,fig_name)
    gsn_define_colormap(wks,"matlab_jet") 
    ;has 64 colors blue green yellow,red

if(nv.eq.0) then 
    do nc = 0,ncase-1,1
    files1 = systemfunc("ls " + filein(nc) + "*.daily.PRECC.nc")
    files2 = systemfunc("ls " + filein(nc) + "*.daily.PRECL.nc")
    do nt = 0, nyear-1 ,1
        print("handle with " + files1(nt))
        f    = addfile(files1(nt),"r")
        temp = f->PRECC(select_time,{lats(na):latn(na)},{lonl(na):lonr(na)})
        f    = addfile(files2(nt),"r")
        temp = temp + f->PRECL(select_time,{lats(na):latn(na)},{lonl(na):lonr(na)})
        var(nc,nt,:) = wgt_areaave_Wrap(temp,1.0,1.0,0) ;get var(nvar,ncase,nyear,365,nlev)
    end do
    end do
    var = var*1000*3600*24 ;convert m/s to mm/day
else
    do nc = 0,ncase-1,1
    files = systemfunc("ls " + filein(nc) + "*.daily."+draw_var(nv)+".nc")
    do nt = 0, nyear-1 ,1
        print("handle with " + files(nt))
        f   = addfile(files(nt),"r")
        temp = f->$draw_var(nv)$(:,{lev},{lats(na):latn(na)},{lonl(na):lonr(na)})
        ;temp = mask(temp,phis.gt.3000,True) 
        var(nc,nt,:) = wgt_areaave_Wrap(temp,1.0,1.0,0) ;get var(nvar,ncase,nyear,365,nlev)
    end do
    end do
end if

if(nv.ge.1) then 
    cf = ezfftf_n(dim_avg_n(var,1),1)
    cf(:,:,4:,:)  = 0
    ave = ezfftb_n(cf,cf@xbar,2)
else
    ave = dim_avg_n_Wrap(var,1)
end if

if(nv.eq.3) then
var = var - 273.15   ;convert unit from k to c
ave = ave - 273.15
end if

var&time = select_time + 1

printVarSummary(var)
printVarSummary(ave)

;===============================================
;draw the figure
;=============================================
    resxy     = True
    resxy@vpWidthF  = 0.6
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
    resxy@tiXAxisFont   = font_index
    resxy@tiYAxisFont   = font_index
    resxy@tiXAxisFontHeightF   = tick_font
    resxy@tiYAxisFontHeightF   = tick_font
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
    resxy@tmXBMode           = "Manual"
    ;resxy@tmXBTickSpacingF   = 30
    resxy@tmXBTickStartF     = select_time(0)
    resxy@tmXBTickEndF       = select_time(nday-1)
    ;resxy@tmXBMode           = "Explicit"
    ;resxy@tmXBValues         = ispan(0,360,30)
    ;resxy@tmXBLabels         = ispan(0,360,30)
    ;resxy@tmXBLabels   = (/"Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec"/)
    
    resxy@tmXMajorGrid                   = True
    resxy@tmXMajorGridLineDashPattern    = 2
    resxy@tmXMajorGridLineColor          = "gray"
    resxy@tmXMajorGridThicknessF         = 2.0
                             
    resxy@gsnStringFont         = font_index
    resxy@gsnStringFontHeightF  = subtitle_font
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
    resp@gsnPaperOrientation = "portrait" ;"landscape";

np = -1
do nc = 0, ncase-1,1
    np = np + 1 
    resxy@gsnLeftString   = pre_case(np) + ") " + case(nc) + "  " + draw_var(nv) + " (" + unit(nv) + ")" 
    resxy@gsnRightString  = "(" + lats(na) + "-" + latn(na) + "N," + lonl(na) + "-" + lonr(na) + "E)";lev(nl) + "hPa " + 
    resxy@trYMinF = floor(min(var))
    resxy@trYMaxF =  ceil(max(var))
    plot(np)  = gsn_csm_xy(wks,var&time,(/var(nc,:,:)/),resxy)
    plot2(np) = gsn_csm_xy(wks,var&time,(/ave(0 ,  :)/),resxy2)
    overlay(plot(np),plot2(np))
    ;plot3(np) = gsn_csm_xy(wks,var&day,(/ave(nc,nv,  :)/),resxy3)
    ;overlay(plot(np),plot3(np))
end do
    
   gsn_panel(wks,plot,(/ncase,1/),resp)
end

