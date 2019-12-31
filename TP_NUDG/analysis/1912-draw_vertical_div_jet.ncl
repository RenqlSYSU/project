;*********************************************************************
;calc and draw the climatic vertical distribution of divergence(shaded),
;U (contour) and v&w(vector).
;                                                     by Ql Ren
;                                                    2019/12/12
;*********************************************************************************
begin
path   = "/home/ys17-19/renql"
case_option = 0 ;0 means NUDG experiment; 1 means noTP experiment
vector = True ;False
draw_var = (/"Div","Vor"/)  ;this is the variable used to draw, can used as center string
var_name = (/"U","V","OMEGA"/)
unit     = (/"10~S~-6~N~s~S~-1~N~","10~S~-6~N~s~S~-1~N~"/)    ;this can be used as right string 
diff = False;True
dnv = 0

season = True
if(season) then 
months   = (/"DJF","MAM","JJA","SON"/)
else
months = (/"JF","MA","MJ","JA"/)
;months  = (/"May","Jun","Jul","Aug"/)
end if
nmonth = dimsizes(months)
month_s = 1
month_e = 12
lev = (/1000,850,700,600,500, 400,300,250,200,150,100/)
;lev  = (/1000,925,850,700,600,500,400,300,200,100/)

;area: east Asia
area = (/"TP","EA","NH"/)
lats = (/15  ,15  , 0  /);20 ;-20
latn = (/55  ,55  , 89 /) ;55 ;60
lonl = (/70  ,110 , 10 /) ;30 ;20
lonr = (/105 ,130 ,350 /) ;150;220
lat_spacing = 15
lon_spacing = 30
na = 1

if(case_option.eq.0) then
    case1 = (/"CTRL","NG6h"/)
    if(diff) then
    case2 = (/"CTRL","NG6h-CTRL"/)
    else
    case2 = case1
    end if
    filein  = path +(/"/model/AMIP-CTRL/AMIP_C5PM.cam.h1.YEAR.",\
                      "/model/TP-NUDG-6h/AMIP_C5PM_TP_NUDG.cam.h1.YEAR."/)
;                      "/model/TP-NUDG-24h/AMIP_C5PM_TP_NUDG.cam.h1.YEAR."/)
    years    = ispan(1979,2005,1)   ;the calculate time
end if

if(case_option.eq.1) then
    ;remove TP experiment, sample number is ten years 
    case1 = (/"TP_CTRL","TP_CR"/) ;TP_CTRL is 0021-0030; TP_CR is 0041-0050
    if(diff) then
    case2 = (/"CTRL","noTP-CTRL"/)
    else
    case2 = (/"CTRL","noTP"/)
    end if
    filein  = path + "/model/" + case1 + "/pro/" + case1 + ".cam.h1.YEAR.00" ; + year + ".daily." + var_name + ".nc"
    years   = ispan(21,30,1)   ;the calculate time
end if

if(case_option.eq.2) then
    case1 = (/"F2000","NG58"/)
    if(diff) then
    case2 = (/"F2000","NG58-F2000"/)
    else
    case2 = case1
    end if
    filein  = path +(/"/model/F2000_CAM5/F2000_CAM5.cam.h1.YEAR.00",\
                      "/model/F2000_CAM5_NUDG5-8/F2000_CAM5.cam.h1.YEAR.00"/)
;                      "/model/F2000_CAM5_NUDG/F2000_CAM5.cam.h1.YEAR.00",\
    years   = ispan(1,29,1)   ;the calculate time
end if
fig_name = path + "/project/TP_NUDG/analysis/fig/" + case1(0) + "_vertical_div_jet"
pre_case = (/"(a","(b","(c","(d","(e","(f","(g","(h","(i","(j","(k","(l"/)
fig_out  = "pdf";"x11""eps";

files := systemfunc("ls "+filein(0)+"*.daily.U.nc" )
f      = addfile(files(0),"r")
vars   = f->U(0,{lev},{lats(na):latn(na)},{lonl(na):lonr(na)})
time   = cd_calendar(f->time,0)
select_time = ind(time(:,1).ge.month_s.and.time(:,1).le.month_e)
ntime = dimsizes(select_time)
nyear = dimsizes(years)
ncase = dimsizes(case1)
nlat  = dimsizes(vars&lat)
nlon  = dimsizes(vars&lon)
nlev  = dimsizes(lev)
u   = new((/ncase,nyear,ntime,nlev,nlat,nlon/),float)
v   = new((/ncase,nyear,ntime,nlev,nlat,nlon/),float)
u1  = new((/ncase,nmonth,nlev,nlat/),float)
v1  = new((/ncase,nmonth,nlev,nlat/),float)
var = new((/ncase,nmonth,nlev,nlat/),float)
if(vector) then 
w   = new((/ncase,nyear,ntime,nlev,nlat,nlon/),float)
w1  = new((/ncase,nmonth,nlev,nlat/),float)
end if

;the number of the panel plot
    plot   = new((/ncase*nmonth/),graphic)
    vect   = new((/ncase*nmonth/),graphic)
    line   = new((/ncase*nmonth/),graphic)

;define the thickness and height of the title,tick and label
   tick_font      = 0.03
   subtitle_font  = 0.03
   title_font     = 0.03
   vcref_font     = 0.02
   tick_thick     = 8.0
   subtitle_thick = 8.0  ;not have this attribute
   title_thick    = 5.0
   vc_thick       = 2.0
   mp_thick       = 3.0
   font_index     = 22

;================================================================
;read data
;====================================================================
do nf = 0 , ncase-1,1
do nt = 0,nyear-1,1
    nv = 0
    files := systemfunc("ls "+filein(nf)+"*.daily."+var_name(nv)+".nc" )
    f      = addfile(files(nt),"r")
    u(nf,nt,:,:,:,:) = f->$var_name(nv)$(select_time,{lev},{lats(na):latn(na)},{lonl(na):lonr(na)})
    nv = 1
    files := systemfunc("ls "+filein(nf)+"*.daily."+var_name(nv)+".nc" )
    f      = addfile(files(nt),"r")
    v(nf,nt,:,:,:,:) = f->$var_name(nv)$(select_time,{lev},{lats(na):latn(na)},{lonl(na):lonr(na)})
    if(vector) then 
    nv = 2
    files := systemfunc("ls "+filein(nf)+"*.daily."+var_name(nv)+".nc" )
    f      = addfile(files(nt),"r")
    w(nf,nt,:,:,:,:) = f->$var_name(nv)$(select_time,{lev},{lats(na):latn(na)},{lonl(na):lonr(na)})
    end if
end do 
end do

if(dnv.eq.0) then 
var0 = uv2dv_cfd ( u, v, u&lat, u&lon, 2 )
var0 = var0*10^6
end if

if(dnv.eq.1) then 
var0 = uv2vr_cfd ( u, v, u&lat, u&lon, 2 )
var0 = var0*10^6
end if

var0 := dim_avg_n_Wrap(var0,1)
u := dim_avg_n_Wrap(u,1)
v := dim_avg_n_Wrap(v,1)
var0 := dim_avg_n_Wrap(var0,4)
u := dim_avg_n_Wrap(u,4)
v := dim_avg_n_Wrap(v,4)
if(vector) then 
w := dim_avg_n_Wrap(w,1)
w := dim_avg_n_Wrap(w,4)
end if

time1 = cd_calendar(f->time(select_time),0)
do nm = 0, nmonth-1, 1
if(nm.eq.0.and.season) then 
mons = 12  
mone = 2
select_time1 := ind(time1(:,1).eq.12.or.time1(:,1).le.2)
else
mons = 3*nm  
mone = 3*nm + 2
select_time1 := ind(time1(:,1).ge.mons.and.time1(:,1).le.mone)
end if
var(:,nm,:,:) = dim_avg_n_Wrap(var0(:,select_time1,:,:),1)
u1(:,nm,:,:) = dim_avg_n_Wrap(u(:,select_time1,:,:),1)
v1(:,nm,:,:) = dim_avg_n_Wrap(v(:,select_time1,:,:),1)
if(vector) then 
w1(:,nm,:,:) = dim_avg_n_Wrap(w(:,select_time1,:,:),1)
end if
end do
printVarSummary(var)

if(diff) then 
var(1,:,:,:)= var(1,:,:,:) - var(0,:,:,:)   ;get the difference between CTRL and NUDG
u1(1,:,:,:)= u1(1,:,:,:) - u1(0,:,:,:)   ;get the difference between CTRL and NUDG
v1(1,:,:,:)= v1(1,:,:,:) - v1(0,:,:,:)   ;get the difference between CTRL and NUDG
if(vector) then 
    w1(1,:,:,:)= w1(1,:,:,:) - w1(0,:,:,:)   ;get the difference between CTRL and NUDG
    copy_VarMeta(vars(:,:,0),w1(0,0,:,:))
end if
u1(0,:,:,:) = mask(u1(0,:,:,:), u1(0,:,:,:).lt.0.0 , False)
else
u1 = mask(u1, u1.lt.0.0 , False)
end if
printVarSummary(var)
w1 = - w1*100
copy_VarMeta(vars(:,:,0),var(0,0,:,:))
copy_VarMeta(vars(:,:,0),u1(0,0,:,:))
copy_VarMeta(vars(:,:,0),v1(0,0,:,:))

;==============================================================
;draw figure in twelve months (one month is a fgure)
;==================================================================
    wks = gsn_open_wks(fig_out,fig_name)
    gsn_define_colormap(wks,"BlueDarkRed18")
    ;blue and red, has 18 colors,2-10 is blue
    ;gsn_define_colormap(wks,"GMT_polar") 
    ;blue and red, has 20 colors, 2-11 is bule, 12-21 is red
    ;gsn_define_colormap(wks,"matlab_jet") 
    ;has 64 colors blue green yellow,red

    res                  = True
    res@gsnFrame         = False
    res@gsnDraw          = False
    res@gsnMaximize      = True
    res@gsnAddCyclic     = False
    res@gsnPaperOrientation = "portrait"
    res@vpHeightF        = 0.6
    res@vpWidthF         = 0.8

;--------set the form of the coordinate-----------------
    res@tmBorderThicknessF  = mp_thick
    res@tmXBLabelFont   = font_index
    res@tmYLLabelFont   = font_index
    res@tmXBLabelFontHeightF     = tick_font
    ; res@tmXBLabelFontAspectF   = 1.6
    res@tmXBLabelFontThicknessF  = tick_thick
    res@tmYLLabelFontHeightF     = tick_font
    ; res@tmYLLabelFontAspectF   = 1.3
    res@tmYLLabelFontThicknessF  = tick_thick
    
    res@tmXBMajorThicknessF     = mp_thick
    res@tmYLMajorThicknessF     = mp_thick
    ;  res@tmXBMajorOutwardLengthF = -0.
    ;  res@tmYLMajorOutwardLengthF = -0.
    ;  res@tmYRMajorOutwardLengthF = -0.
      
    res@tmXBMinorThicknessF     = mp_thick
    res@tmYLMinorThicknessF     = mp_thick
    ;  res@tmXBMinorOutwardLengthF = -0.
    ;  res@tmYLMinorOutwardLengthF = -0.
    ;  res@tmYRMinorOutwardLengthF = -0.

;---------set the attribute of the contour-----------------
    res@cnFillOn         = True
    res@cnLinesOn        = False
    res@cnLineLabelsOn   = False
    res@cnInfoLabelOn    = False
;    res@cnLineThicknessF = 1.0
    
    ;res@cnLevelSelectionMode = "ExplicitLevels"   ; set manual contour levels
    ;res@cnLevels             = (/-32,-24,-16,-10,-5,-1,-0.5,-0.1, 0, 0.1, 0.5, 1, 5, 10, 16, 24, 32/)

    res@gsnSpreadColors      = False              ; turn of the full color map
    res@gsnSpreadColorStart  = 9    ;This resource only takes effect if gsnSpreadColors is set to True 
    res@gsnSpreadColorEnd    = 0                ; end at color 200
    
    res@lbLabelBarOn         = False ;True   ; the label bar
    res@lbOrientation        = "Horizontal"
    res@lbLabelFontHeightF   = tick_font              ; make labels smaller
    res@lbLabelFont          = font_index              ; make labels smaller
    res@lbLabelStride        = 1
    res@lbBoxLineThicknessF  = mp_thick-1

    res@pmLabelBarZone               = 0
    res@pmLabelBarOrthogonalPosF     = 0
    res@pmLabelBarParallelPosF       = 0.5

;------------------------------------------------------------------
;draw the subgraph and change the left and right titles
;-----------------------------------------------------------------------
    res@gsnCenterString      = ""
    res@gsnStringFont        = font_index
    res@gsnStringFontHeightF = subtitle_font
    
    res@gsnPresHgtHeightLabelOn = False
    res@tiYAxisOn = False
    res@tmYRLabelsOn = False
    ;res@tmYLLabelsOn = False
    res@tmXBMode = "Manual"
    res@tmXBTickSpacingF   = 10
    res@tmXBTickStartF     = lats(na)
    res@tmXBTickEndF       = latn(na)

;---------set the attribute of the contour-----------------
    res_line = res
    res_line@cnFillOn         = False
    res_line@cnLinesOn        = True
    res_line@cnLineLabelsOn   = True
    res_line@cnInfoLabelOn    = False
    res_line@cnLineThicknessF = 3.0
    res_line@cnLineColor      = "black"
    res_line@gsnContourZeroLineThicknessF = 0
    res_line@gsnContourNegLineDashPattern = 5
    res_line@gsnPresHgtHeightLabelOn = False
    res_line@tiYAxisOn = False
    res_line@tmYLLabelsOn = False
   
    res_line@cnLevelSelectionMode = "AutomaticLevels"   ; set manual contour levels
    res_line@cnLevelSpacingF       = 10        ; set contour interval
    res_line@gsnCenterString      = ""
    res_line@gsnLeftString    = " "
    res_line@gsnRightString   = " "

;----------set the attribute of the vector-------------------
    resvc                  = res
    resvc@cnFillOn         = False
    resvc@cnLinesOn        = False
    resvc@cnLineLabelsOn   = False
    resvc@cnInfoLabelOn    = False
    resvc@gsnPresHgtHeightLabelOn = False
    resvc@tiYAxisOn = False
    resvc@tmYLLabelsOn = False
    resvc@vcGlyphStyle               = "LineArrow"  ;"CurlyVector"
    resvc@vcVectorDrawOrder          = "PostDraw"
    resvc@vcLineArrowThicknessF      = vc_thick
    resvc@vcMinDistanceF             = 0.02
    resvc@vcLineArrowHeadMinSizeF   = 0.015
    resvc@vcLineArrowHeadMaxSizeF   = 0.015
    resvc@vcRefAnnoOn                = False
    resvc@vcRefLengthF               = 0.05           ; define length of vec ref
    resvc@vcRefAnnoZone              = 0     ;set the positional origin in the center of the plot
    ;resvc@vcRefAnnoSide              = "Bottom"          
    ;resvc@vcRefAnnoJust              = "TopRight"
    resvc@vcRefAnnoOrthogonalPosF    = -0.85           ; move ref vector along Y-axis
    resvc@vcRefAnnoParallelPosF      = 0.25            ; move ref vector along X-axis
    resvc@vcRefAnnoPerimOn           = False         ;close the border
    resvc@vcRefAnnoFont              = font_index
    resvc@vcRefAnnoFontHeightF       = vcref_font
    resvc@vcRefAnnoString1On         = False
    resvc@vcRefAnnoString2On         = True
    resvc@vcRefAnnoString2           = "8(m/s,10~S~-2~N~hPa/s)"
    resvc@vcRefMagnitudeF            = 8.0             ; define vector ref mag
    resvc@gsnCenterString      = ""
    resvc@gsnLeftString    = " "
    resvc@gsnRightString   = " "

    resp  = True    
    resp@gsnMaximize   = True
    resp@gsnPaperOrientation = "landscape";"portrait"
    resp@gsnPanelLabelBar = True
    resp@pmLabelBarWidthF  = 0.8
    resp@pmLabelBarHeightF = 0.04
    resp@lbLabelFontHeightF = 0.015
    resp@lbLabelFontThicknessF = tick_thick
    resp@txFontHeightF = title_font 
    resp@txFontThicknessF = title_thick

if(dnv.eq.0) then 
    minlvl = (/-7,-7/);
    maxlvl = (/ 7, 7/);
    spacig = (/ 1, 1/);
end if

if(dnv.eq.1) then 
    minlvl = (/-21,-21/);
    maxlvl = (/ 21, 21/);
    spacig = (/ 3 , 3 /);
end if

    res@cnFillColors  = (/2,3,4,5,6,7,8,9,12,13,14,15,16,17,18,19/)
;do nm = 0,nmonth-1,1
do nc = 0, ncase-1,1
;if((nc.gt.0).or.(nv.eq.1).or.(nv.eq.2).or.(nv.eq.3)) then
;    res@cnFillPalette = "BlueDarkRed18"
;   ;blue and red, has 18 colors,2-10 is blue
;else
;   res@cnFillPalette = "precip3_16lev"
;   ;blue and yellow and red, has 17 colors,2-18,2 is white
;end if 
    do nm = 0,nmonth-1,1
        np = nmonth*nc + nm
        ;np = ncase*nm + nc
        res@gsnLeftString   = pre_case(np) + ") " + case2(nc) + " " + months(nm) 
        res@gsnRightString  = lonl(na) + "-" + lonr(na) + "E " + draw_var(dnv) + " (" + unit(dnv) + ")"
        res@gsnSpreadColors       = False              ; turn of the full color map
        res@cnLevelSelectionMode  = "ManualLevels" ; set manual contour levels
        res@cnLevelSpacingF       := spacig(nc)        ; set contour interval
        res@cnMaxLevelValF        := maxlvl(nc) ;32    ; set max contour level
        res@cnMinLevelValF        := minlvl(nc)       ; set min contour level
        plot(np) = gsn_csm_pres_hgt(wks,var(nc,nm,:,:),res)
        if(diff.and.(nc.eq.1)) then
        res_line@cnLevelSpacingF       = 8        ; set contour interval
        else
        res_line@cnLevelSpacingF       = 8        ; set contour interval
        end if
        line(np) = gsn_csm_pres_hgt(wks,u1(nc,nm,:,:),res_line)
        overlay(plot(np),line(np))
        vect(np) = gsn_csm_pres_hgt_vector(wks,var(nc,nm,:,:),v1(nc,nm,:,:),w1(nc,nm,:,:),resvc)
        overlay(plot(np),vect(np))
        print(max(var(nc,nm,:,:)))
        print(min(var(nc,nm,:,:)))
        print(max(v1(nc,nm,:,:)))
        print(min(v1(nc,nm,:,:)))
        print(max(w1(nc,nm,:,:)))
        print(min(w1(nc,nm,:,:)))
    end do 
end do ;for the month 

    gsn_panel(wks,plot,(/ncase,nmonth/),resp)
    ;gsn_panel(wks,plot,(/nmonth,ncase/),resp)
;end do
end

