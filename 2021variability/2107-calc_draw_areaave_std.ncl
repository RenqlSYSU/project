;***************************************************************************
;filename:  calc_areaave_std.ncl 
;1. read variable from multiple case, get var(ncase,nyear,ntime,nlat,nlon)
;2. use Four levels of nested loops to do space smooth
;3. time filter, calc variance, areaave, difference and plot
;                                           by Ql_Ren
;                                           2021/07/23
;******************************************************************************
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/shea_util.ncl"

begin
path = "/home/ys17-19/renql"
lats = (/15 ,10 /)
latn = (/50 ,60 /)
lonl = (/90 ,90 /)
lonr = (/150,150/)
lats1 = (/ 27, 32 , 25 , 40 /)
latn1 = (/ 37, 32 , 25 , 40 /)
lonl1 = (/110, 118, 115, 115/)
lonr1 = (/125, 118, 115, 115/)
var_name = (/"PRECC","PRECL"/)
na = 1
ndist = 10
diff  = False
perc  = False;True

ifrunave= 1 ;0 = not filter; 1 = runave; 2 = bandpass filter; 3 = Fourier filter
nave = 5 ;five days run ave
ca = 80.0;90.0 ;day
cb = 8.0 ;2.5  ;10.0 ;
month_s = 6 
month_e = 8

case    = (/"CTRL","NG15S","NG4550","NGEN"/) ;
filein  = path +(/"/model/AMIP-CTRL/AMIP_C5PM.cam.h1.YEAR.",\
                  "/model/AMIP_NG15S/AMIP_CAM5_NG15S.cam.h1.YEAR.",\
                  "/model/AMIP_NG4550/AMIP_CAM5_NG4550.cam.h1.YEAR.",\
                  "/model/AMIP_NGEN/AMIP_CAM5_NGEN.cam.h1.YEAR."/)
years   = ispan(1979,2005,1)   ;the calculate time
fig_name = path + "/project/2021variability/fig/smth_JJApreci_vari"
pre_case   = (/"(a","(b","(c","(d","(e","(f","(g","(h","(i","(j","(k","(l"/)
fig_out    ="pdf";"x11"

f    = addfile(filein(0)+years(0)+".daily.PRECC.nc","r")
time = cd_calendar(f->time,0)
if(month_s.gt.month_e) then 
select_time = ind(time(:,1).ge.month_s.or.time(:,1).le.month_e)
else
select_time = ind(time(:,1).ge.(month_s-1).and.time(:,1).le.(month_e+1))
end if
vars := f->PRECC(0,{lats(na):latn(na)},{lonl(na):lonr(na)})
ntime = dimsizes(select_time)
ncase = dimsizes(case)
nyear = dimsizes(years)
narea = dimsizes(lats1)
nlat  = dimsizes(vars&lat)
nlon  = dimsizes(vars&lon)
var  = new((/ncase,nyear,ntime,nlat,nlon/),float)
smth = new((/ncase,ndist,narea/),float)
plot = new((/6/),graphic)

;==========================================================
;read the data
;============================================================
do nf = 0,3,1
do nt = 0,nyear-1,1
    nv = 0
    files := systemfunc("ls "+filein(nf)+"*.daily."+var_name(nv)+".nc" )
    f  = addfile(files(nt),"r")
    print(" handle with " + files(nt))
    var(nf,nt,:,:,:) = f->$var_name(nv)$(select_time,{lats(na):latn(na)},{lonl(na):lonr(na)})

    nv = 1
    files := systemfunc("ls "+filein(nf)+"*.daily."+var_name(nv)+".nc" )
    f  = addfile(files(nt),"r")
    var(nf,nt,:,:,:) = var(nf,nt,:,:,:) + f->$var_name(nv)$(select_time,{lats(na):latn(na)},{lonl(na):lonr(na)})
end do
end do
var = linmsg_n(var,(/-1/),3)
var = var*3600*24*1000 ;mm/day
printVarSummary(var)

;==========================================================
;space smooth, time filter and variance
;============================================================
time := cd_calendar(f->time(select_time),0)
select_time := ind(time(:,1).ge.month_s.and.time(:,1).le.month_e)
lat = vars&lat
lon = vars&lon
do ndis = 0, ndist-1, 1
threshold = ndis
if(threshold.gt.0) then 
do nla = 1, nlat-2, 1
do nlo = 1, nlon-2, 1
    nsmth = 1
    do nla1 = 0, nlat-1, 1
    do nlo1 = 0, nlon-1, 1
    distance = sqrt((lat(nla1)-lat(nla))^2+(lon(nlo1)-lon(nlo))^2)
    if(distance.le.threshold.and.distance.gt.0) then
        var(:,:,:,nla,nlo) = var(:,:,:,nla,nlo) + var(:,:,:,nla1,nlo1)
        nsmth = nsmth + 1
    end if
    end do
    end do
    var(:,:,:,nla,nlo) = var(:,:,:,nla,nlo)/nsmth
    print(lat(nla)+"N, "+lon(nlo)+"E : threshold = "+threshold+", nsmth = "+nsmth)
end do
end do 
end if

dimh = 2
if(ifrunave.eq.1) then 
    opt1= -1 ;end-point option.  -1, utilize cyclic conditions; 0, set missing values; 1, utilize reflective (symmetric) conditions
    var = runave_n(var, nave, opt1, dimh)
    var = var - conform(var,dim_avg_n(var,1),(/0,2,3,4/))  ;remove annual cycle
end if
if(ifrunave.eq.2) then 
    fca = 1.0/ca
    fcb = 1.0/cb
    opt2= False
    if(rmv) then 
    var = var - conform(var,dim_avg_n(var,1),(/0,2,3,4/))  ;remove annual cycle
    end if
    var = bw_bandpass_filter(var,fca,fcb,opt2,dimh)
end if

do na1 = 0, narea-1, 1
    term = dim_avg_n(dim_avg_n(var(:,:,select_time,:,:)^2,2),1)
    copy_VarMeta(vars,term(0,:,:))
    smth(:,ndis,na1) = wgt_areaave(term(:,{lats1(na1):latn1(na1)},{lonl1(na1):lonr1(na1)}),1.0,1.0,0)
end do
end do
delete(var)

if(diff) then 
do nc = 1, ncase-1, 1
    if(perc) then 
    smth(nc,:,:) = (smth(nc,:,:)-smth(0,:,:))/smth(0,:,:)
    else
    smth(nc,:,:) = smth(nc,:,:)-smth(0,:,:)
    end if
end do
end if

;************************************************
; plotting
;************************************************
;define the thickness and height of the title,tick and label
wks  = gsn_open_wks(fig_out,fig_name)             ; send graphics to PNG file

    tick_font      = 0.015
    subtitle_font  = 0.018
    title_font     = 0.03
    label_font     = 0.018
    tick_thick     = 5.0
    subtitle_thick = 5.0  ;not have this attribute
    title_thick    = 5.0
    label_thick    = 1.5
    mp_thick       = 2.0
    font_index     = 22

    resxy     = True
    ;resxy@vpWidthF  = 0.45
    ;resxy@vpHeightF = 0.15
    resxy@gsnFrame      = False
    resxy@gsnDraw       = False
    resxy@gsnMaximize   = True
    resxy@gsnPaperOrientation = "portrait"
    
    resxy@tmYLAutoPrecision = False
    resxy@tmYLPrecision     = 2 
    ;resxy@tiXAxisOn = False
    ;resxy@tiYAxisOn = False
    resxy@tiXAxisString = "smooth radius"
    resxy@tiYAxisString = "preci variance"           ; yaxis
    ;resxy@trYMinF = 0
    ;resxy@trYMaxF = max(ave(:,0,:))
    resxy@trXMinF = 0
    resxy@trXMaxF = ndist-1
    resxy@xyDashPatterns    = 0 
    resxy@xyLineThicknesses = 3.0 
    resxy@xyLineColors      = (/"black","red","blue","grey"/)
    
    ;legend
    resxy@pmLegendDisplayMode = "Always"
    resxy@lgOrientation       = "Vertical"
    resxy@lgPerimOn           = False
    resxy@lgPerimFill         = "SolidFill"
    resxy@lgPerimFillColor    = "white"
    resxy@xyExplicitLegendLabels = case
    resxy@lgLabelFontColor       = "black"
    resxy@lgLabelFontThicknessF  = label_thick
    resxy@lgLabelFontHeightF     = label_font
    resxy@lgLabelFont            = font_index
    resxy@pmLegendWidthF   = 0.1
    resxy@pmLegendHeightF  = 0.25
    resxy@pmLegendZone  = 0   ;set the positional origin in the center of the plot
    resxy@pmLegendOrthogonalPosF  = - 0.28   ; move ref vector along Y-axis
    resxy@pmLegendParallelPosF    = 0.28    ; move ref vector along X-axis
    
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
                             
    resxy@gsnStringFont        = font_index
    resxy@gsnStringFontHeightF = subtitle_font
    resxy@gsnCenterString  = "";lev(nl) + " hPa"

    resp  = True    
    resp@gsnMaximize   = True
    resp@gsnPaperOrientation = "portrait" ;"landscape";
    resp@txFontHeightF = title_font 
    resp@txFontThicknessF = title_thick
    if(ifrunave.eq.0) then 
    resp@txString      = " no filter"
    end if
    if(ifrunave.eq.1) then 
    resp@txString      = " runave="+nave+"day"
    end if
    if(ifrunave.eq.2) then 
    resp@txString      = " bw_filter="+cb+"-"+ca+"day"
    end if
    if(ifrunave.eq.3) then 
    resp@txString      = " Fourier_filter="+cb+"-"+ca+"day"
    end if

np = 0
do na1 = 0, narea-1, 1
    resxy@gsnLeftString   = pre_case(np) + ") " + month_s + "-" + month_e 
    resxy@gsnRightString   = "(" + lats1(na1) + "-" + latn1(na1) + "N," + lonl1(na1) + "-" + lonr1(na1) + "E)"
    plot(np) = gsn_csm_xy(wks, ispan(0,ndist-1,1), smth(:,:,na1), resxy)    
    np = np + 1
end do 
gsn_panel(wks,plot,(/2,2/),resp)
end

