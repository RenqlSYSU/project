;*********************************************************************
;Draw the vertical section of monthly EKE (shaded)
;                                                     by Ql Ren
;                                                    2019/08/12
;*********************************************************************************
begin
draw_var = (/"Qd"   ,"Qeddy","A"/)  ;this is the variable used to draw, can used as center string
unit     = (/"K/Day","K/Day","10~S~-11~N~s~S~-2~N~"/)  ;this can be used as right string 
draw_var = (/"f_Qd" ,"f_Qeddy","A"/)  ;this is the variable used to draw, can used as center string
unit     = "10~S~-11~N~" + (/"s~S~-2~N~","s~S~-2~N~","s~S~-2~N~"/)  ;this can be used as right string 
;months   = (/"Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec"/)
months   = (/"MJ","JA"/)
lev  = (/1000,925,850,800,750,700,650,600,550,500,450,400,350,300,250,200,150,100/) ;hPa
dnc  = 0    ;chose which case to draw
diff = True ;whether draw the difference

;area: east Asia
area = (/"TP","EA"/)
lats = (/15 ,15 /);20 ;-20
latn = (/55 ,55 /) ;55 ;60
lonl = (/70 ,105 /) ;30 ;20
lonr = (/105,140/) ;150;220
month_s = 5
month_e = 8
na = 1

;the path of the file
path    = "/home/ys17-19/renql/project/TP_NUDG"
case_option = 0 ;0 means NUDG experiment; 1 means noTP experiment

if(case_option.eq.0) then
    case1   = (/"CTRL","NUDG6h","NUDG24h"/)
    case2   = (/" CTRL"," NG6h-CTRL", " NG24h-CTRL"/)
    gtopfile = path + "/analysis/mdata/" + "gtopo30_0.9x1.25.nc"
    nsample = 27
else
    case1   = (/"TP_CTRL","TP_CR"/)
    case2   = (/"CTRL"   ,"noTP"/)
    ;gtopfile = "/home/ys17-19/renql/model/TP_CR/USGS-gtopo30_1.9x2.5_noTP_c050602C.nc" 
    gtopfile = "/home/ys17-19/renql/model/TP_CTRL/USGS-gtopo30_1.9x2.5_remap_c050602.nc" 
    nsample = 10
end if
filename = path + "/z_tendency-20191022/mdata/" + case1 + "-Clim_daily_Qd_Qeddy_A.nc"
fig_name = path + "/analysis/fig/1p2x3_mj-ja_height_lat_3Q_" + area(na)
;fig_title  = "Seasonal " + draw_var(nv) + " variance in " ;+season
pre_case   = (/"(a","(b","(c","(d","(e","(f","(g","(h","(i","(j","(k","(l"/)
fig_out    ="pdf";"x11"

f = addfile(gtopfile,"r")
phis   = f->PHIS({lats(na):latn(na)},{lonl(na):lonr(na)})
phis   = phis/9.8   ;convert unit from m2/s2 to gpm
nmonth = dimsizes(months)
nlat   = dimsizes(phis&lat)
nlon   = dimsizes(phis&lon)
nlev   = dimsizes(lev)
nvar   = dimsizes(draw_var)
ncase  = dimsizes(case1)

f = addfile(filename(0),"r")
vars = f->Qd(0,{lev},{lats(na):latn(na)},{lonl(na):lonr(na)})
time = cd_calendar(f->time,0)
select_time = ind(time(:,1).ge.month_s.and.time(:,1).le.month_e)
ntime = dimsizes(select_time)
var   = new((/ncase,nvar,ntime ,nlev,nlat,nlon/),float)
var_m = new((/ncase,nvar,nmonth,nlev,nlat/),float)

;the number of the panel plot
    plot     = new((/nvar*nmonth/),graphic)
    gtop     = new((/nvar*nmonth/),graphic)
    line     = new((/nvar*nmonth/),graphic)

;define the thickness and height of the title,tick and label
   tick_font      = 0.03
   subtitle_font  = 0.03
   title_font     = 0.03
   tick_thick     = 8.0
   subtitle_thick = 8.0  ;not have this attribute
   title_thick    = 5.0
   mp_thick       = 3.0
   font_index     = 22

;================================================================
;read data
;====================================================================
do nc = 0, ncase-1, 1
    f = addfile(filename(nc),"r")
    do nv = 0, nvar-1,1
    var(nc,nv,:,:,:,:)  = f->$var_name(nv)$(select_time,{lev},{lats(na):latn(na)},{lonl(na):lonr(na)}) ;get an array (ncase,nmonth,nlat,nlon)
    end do
end if
var = dim_avg_n(var,5)
var(:,2,:,:,:) = var(:,2,:,:,:)*10^11 ;convert units of A from s^-2 to 10^-11 s^-2
printVarSummary(var)

;calc ealy(56) and late(78) stage of the summer 
select_time := ind(time(:,1).ge.5.and.time(:,1).le.6)
var_m(:,:,0,:,:) = dim_avg_n(var(:,:,select_time,:,:),2)
select_time := ind(time(:,1).ge.7.and.time(:,1).le.8)
var_m(:,:,1,:,:) = dim_avg_n(var(:,:,select_time,:,:),2)
copy_VarMeta(vars(:,:,0),var_m(0,0,0,:,:))

;==============================================================
;draw figure in twelve months (one month is a fgure)
;==================================================================
    wks = gsn_open_wks(fig_out,fig_name)
    ;gsn_define_colormap(wks,"BlueDarkRed18")
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

    res@gsnSpreadColors      = False              ; turn of the full color map
    res@gsnSpreadColorStart  = 9    ;This resource only takes effect if gsnSpreadColors is set to True 
    res@gsnSpreadColorEnd    = 0                ; end at color 200
    
    res@lbLabelBarOn         = True   ; the label bar
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
    res@tmXBMode = "Manual"
    res@tmXBTickSpacingF   = 15
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
   
    res_line@cnLevelSelectionMode = "AutomaticLevels"   ; set manual contour levels
    res_line@cnLevelSpacingF       = 10        ; set contour interval

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

    minlvl = (/-2.4,-2.4,-4 /);
    maxlvl = (/2.4 ,2.4 , 4 /);
    spacig = (/0.3 ,0.3 ,0.5/);

    res@cnFillPalette = "BlueDarkRed18"
   ;blue and red, has 18 colors,2-10 is blue
do nm = 0,nmonth-1,1
do nv = 0,nvar-1,1
    np = nm*nvar + nv
    res@gsnLeftString   = pre_case(np) + ") " + case2(fnc) + " " + months(nm) 
    res@gsnRightString  = lonl(na) "-" lonr(na) + "E " + draw_var(nv) + " (" + unit(nv) + ")"
    plot(np) = gsn_csm_pres_hgt(wks,var(fnc,nv,nm,:,:),res)
    print(max(var(fnc,nv,nm,:,:)))
    print(min(var(fnc,nv,nm,:,:)))
end do 

;    resp@txString      = fig_title
    gsn_panel(wks,plot,(/nmonth,nvar/),resp)
end

