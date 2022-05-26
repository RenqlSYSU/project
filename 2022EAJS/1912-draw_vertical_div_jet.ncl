;*********************************************************************
;calc and draw the climatic vertical distribution of divergence(shaded),
;U (contour) and v&w(vector).
;                                                     by Ql Ren
;                                                    2019/12/12
;*********************************************************************************
begin
path  = "/home/ys17-19/renql/project/2022EAJS"
case_option = 0
vector = True ;False
cont_var = (/"div"/)
shad_var = (/"dudt","U"/)
cnv = 0
snv = 1
draw_var = (/"U","V","OMEGA",cont_var(cnv),shad_var(snv)/)
unit     = (/"10~S~-6~N~s~S~-1~N~","10~S~-6~N~s~S~-1~N~","10~S~-5~N~m/s~S~2~N~"/)    ;this can be used as right string 
diff = True ;False ;
sigmask = True ;False;
siglvl = 0.05

season = False;True
if(season) then 
var_name = (/"season_ave","season_vari"/) ;this is the var_name in data file
months   = (/"DJF","MAM","JJA","SON"/)
imonth = ispan(0,3,1)
else
var_name = (/"month_ave","month_vari"/) ;this is the var_name in data file
months   = (/"Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec"/)
imonth = ispan(5,7,1)
end if
nmonth = dimsizes(imonth)
lev  = (/1000,925,850,700,500,400, 300,200,100/) ;,150,100hPa

;area: east Asia
area = (/"TP","EA","NH"/)
lats = (/15  ,15  , 0  /);20 ;-20
latn = (/55  ,55  , 89 /) ;55 ;60
lonl = (/70  ,120 , 10 /) ;30 ;20
lonr = (/105 ,120 ,350 /) ;150;220
lat_spacing = 15
lon_spacing = 30
na = 1

if(case_option.eq.0) then
    case1   = (/"F2000","JET_NUDG"/)
    if(diff) then
    case2   = (/"F2000","diff"/)
    else
    case2   = (/"F2000","JET_NUDG"/)
    end if
    years    = ispan(1,30,1)   ;the calculate time
    gtopfile= "/home/ys17-19/renql/project/TP_NUDG/analysis/mdata/" + "gtopo30_0.9x1.25.nc"
end if
filename = path + "/data/" + case1 + "-Clim_month_ave_"
fig_name = path + "/fig/" + case1(0) + "_vertical_div_jet"
pre_case = (/"(a","(b","(c","(d","(e","(f","(g","(h","(i","(j","(k","(l"/)
fig_out  = "pdf";"x11""eps";

f = addfile(gtopfile,"r")
phis  = f->PHIS({lats(na):latn(na)},{lonl(na):lonr(na)})
phis  = phis/9.8   ;convert unit from m2/s2 to gpm
nmonth= dimsizes(imonth)
nvar  = dimsizes(draw_var)
nyear = dimsizes(years)
ncase = dimsizes(case1)
nlat  = dimsizes(phis&lat)
nlon  = dimsizes(phis&lon)
nlev  = dimsizes(lev)
var   = new((/ncase,nvar,nmonth,nlev,nlat,nlon/),float)
vari  = new((/ncase,nvar,nmonth,nlev,nlat,nlon/),float)

;the number of the panel plot
    plot   = new((/ncase*nmonth/),graphic)
    vect   = new((/ncase*nmonth/),graphic)
    line   = new((/ncase*nmonth/),graphic)
    divl   = new((/ncase*nmonth/),graphic)
    sig    = new((/ncase*nmonth/),graphic)

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
do nf = 0, ncase-1, 1
do nv = 0, nvar-1, 1
    f = addfile(filename(nf) + draw_var(nv) + ".nc","r")
    var(nf,nv,:,:,:,:)  = f->$var_name(0)$(imonth,{lev},{lats(na):latn(na)},{lonl(na):lonr(na)})
    vari(nf,nv,:,:,:,:) = f->$var_name(1)$(imonth,{lev},{lats(na):latn(na)},{lonl(na):lonr(na)})
end do 
end do
var  := dim_avg_n_Wrap(var ,5)
vari := dim_avg_n_Wrap(vari,5)
;vari = where(abs(vari).le.0.000001,0.001,vari)
prob = ttest(var(0,:,:,:,:),vari(0,:,:,:,:),nyear,var(1,:,:,:,:),vari(1,:,:,:,:),nyear,True,False)
copy_VarMeta(var(0,0,0,:,:),prob(0,0,:,:))
printVarSummary(prob)

if(cnv.eq.0) then 
    var(:,3,:,:,:) = var(:,3,:,:,:)*1000000 ;div*10^6
end if
if(snv.eq.0) then 
    var(:,4,:,:,:) = var(:,4,:,:,:)*100000 ;dudt*10^5
end if
;var(:,0,:,:,:) = mask(var(:,0,:,:,:), var(:,0,:,:,:).ge.0, True);zonal wind
var(:,2,:,:,:) = -var(:,2,:,:,:)*100 ;omega*100
if(diff) then 
    var(1,1:4,:,:,:)= var(1,1:4,:,:,:) - var(0,1:4,:,:,:)   ;get the difference between CTRL and NUDG
end if
if(sigmask) then 
    ;var(1,3:4,:,:,:) = mask(var(1,3:4,:,:,:), prob(3:4,:,:,:).lt.siglvl, True);cont_var, shad_var
    var(1,3,:,:,:) = mask(var(1,3,:,:,:), prob(0,:,:,:).lt.siglvl.or.prob(1,:,:,:).lt.siglvl.or.prob(2,:,:,:).lt.siglvl.or.prob(3,:,:,:).lt.siglvl.or.prob(4,:,:,:).lt.siglvl, True);v
    var(1,4,:,:,:) = mask(var(1,4,:,:,:), prob(0,:,:,:).lt.siglvl.or.prob(1,:,:,:).lt.siglvl.or.prob(2,:,:,:).lt.siglvl.or.prob(3,:,:,:).lt.siglvl.or.prob(4,:,:,:).lt.siglvl, True);w
    v = var(1,1,:,:,:)
    w = var(1,2,:,:,:)
    ;var(1,1,:,:,:) = mask(var(1,1,:,:,:), prob(1,:,:,:).lt.siglvl.or.prob(2,:,:,:).lt.siglvl, True);prob(0,:,:,:).lt.siglvl.or.
    ;var(1,2,:,:,:) = mask(var(1,2,:,:,:), prob(1,:,:,:).lt.siglvl.or.prob(2,:,:,:).lt.siglvl, True);prob(0,:,:,:).lt.siglvl.or.
    var(1,1,:,:,:) = mask(var(1,1,:,:,:), prob(0,:,:,:).lt.siglvl.or.prob(1,:,:,:).lt.siglvl.or.prob(2,:,:,:).lt.siglvl, True);v
    var(1,2,:,:,:) = mask(var(1,2,:,:,:), prob(0,:,:,:).lt.siglvl.or.prob(1,:,:,:).lt.siglvl.or.prob(2,:,:,:).lt.siglvl, True);w
end if
printVarSummary(var)

;==============================================================
;draw figure in twelve months (one month is a fgure)
;==================================================================
    wks = gsn_open_wks(fig_out,fig_name)

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
    
    res@lbLabelBarOn         = True   ;False ; the label bar
    res@lbOrientation        = "Horizontal"
    res@pmLabelBarWidthF  = 0.7
    res@pmLabelBarHeightF = 0.08
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
    res_line@cnLineThicknessF = 4.0
    res_line@cnLineColor      = "black"
    res_line@gsnContourZeroLineThicknessF = 8
    res_line@gsnContourNegLineDashPattern = 5
    res_line@gsnPresHgtHeightLabelOn = False
    res_line@tiYAxisOn = False
    res_line@tmYLLabelsOn = False
   
    levels = new((/10/),integer)
    levels(0) = 0
    levels(1) = 5
    levels(2:9) = ispan(10,80,10)
    res_line@cnLevelSelectionMode = "ExplicitLevels"    ;set manual contour levels
    res_line@cnLevels             = levels       ; set contour interval
    ;res_line@cnLevelSelectionMode = "AutomaticLevels"   ; set manual contour levels
    ;res_line@cnLevelSpacingF       = 10        ; set contour interval
    res_line@gsnCenterString      = ""
    res_line@gsnLeftString    = " "
    res_line@gsnRightString   = " "

;---------set the attribute of the contour-----------------
    div_line = True
    div_line@gsnFrame         = False
    div_line@gsnDraw          = False
    div_line@cnFillOn         = False
    div_line@cnLinesOn        = True
    div_line@cnLineLabelsOn   = True
    div_line@cnInfoLabelOn    = False
    div_line@cnLineThicknessF = 3.5
    div_line@cnLineColor      = "blueviolet";"gray40"
    div_line@gsnContourZeroLineThicknessF = 0
    div_line@gsnContourNegLineDashPattern = 2
    div_line@gsnContourPosLineDashPattern = 0
    ;div_line@gsnPresHgtHeightLabelOn = False
    div_line@tiYAxisOn = False
    div_line@tmYLLabelsOn = False
   
    div_line@cnLevelSelectionMode = "AutomaticLevels"   ; set manual contour levels
    div_line@cnLevelSpacingF      = 1.0; 0.5;      ; set contour interval
    div_line@gsnCenterString   = ""
    div_line@gsnLeftString     = " "
    div_line@gsnRightString    = " "

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
    resvc@vcRefAnnoOn                = True
    resvc@vcRefLengthF               = 0.05           ; define length of vec ref
    resvc@vcRefAnnoZone              = 0     ;set the positional origin in the center of the plot
    ;resvc@vcRefAnnoSide              = "Bottom"          
    ;resvc@vcRefAnnoJust              = "TopRight"
    resvc@vcRefAnnoOrthogonalPosF    = 0.33           ; move ref vector along Y-axis
    resvc@vcRefAnnoParallelPosF      = 0.5            ; move ref vector along X-axis
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
    resp@gsnPanelLabelBar = False ;True
    resp@pmLabelBarWidthF  = 0.8
    resp@pmLabelBarHeightF = 0.04
    resp@lbLabelFontHeightF = 0.015
    resp@lbLabelFontThicknessF = tick_thick
    resp@txFontHeightF = title_font 
    resp@txFontThicknessF = title_thick

if(snv.eq.0) then 
    minlvl = (/-21,-7,-7/);
    maxlvl = (/ 21, 7, 7/);
    spacig = (/ 3 , 1, 1/);
end if
if(snv.eq.1) then 
    minlvl = (/-35,-14,-14/);
    maxlvl = (/ 35, 14, 14/);
    spacig = (/ 5 ,  2,  2/);
end if

np = -1
do nm = 0,nmonth-1,1
do nc = 0, ncase-1,1
if(snv.le.1.or.nc.gt.0) then
    res@cnFillPalette = "BlueDarkRed18"
    res@cnFillColors  := (/0,1,2,3,4,5,6,7,10,11,12,13,14,15,16,17/)
    ;blue and red, has 18 colors,2-10 is blue
else
    res@cnFillPalette = "precip3_16lev"
    res@cnFillColors  := (/0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16/)
    ;blue and yellow and red, has 17 colors,2-18,2 is white
end if 
;    do nm = 0,nmonth-1,1
        np = np+1
        res@gsnLeftString   = pre_case(np) + ") " + case2(nc) + " " + months(imonth(nm)) 
        res@gsnRightString  = lonl(na) + "E " + draw_var(snv) + " (" + unit(2) + ")";+ "-" + lonr(na) 
        res@gsnSpreadColors       = False              ; turn of the full color map
        res@cnLevelSelectionMode  = "ManualLevels" ; set manual contour levels
        res@cnLevelSpacingF       := spacig(nc)        ; set contour interval
        res@cnMaxLevelValF        := maxlvl(nc) ;32    ; set max contour level
        res@cnMinLevelValF        := minlvl(nc)       ; set min contour level
        plot(np) = gsn_csm_pres_hgt(wks,var(nc,4,nm,:,:),res)
        divl(np) = gsn_csm_contour(wks,var(nc,3,nm,:,:),div_line)
        overlay(plot(np),divl(np))
        
        if(diff.and.(nc.eq.1)) then
        resvc@vcRefMagnitudeF      = 4.0             ; define vector ref mag
        resvc@vcRefAnnoString2     = 4 
        else
        resvc@vcRefMagnitudeF      = 8.0             ; define vector ref mag
        resvc@vcRefAnnoString2     = 8 
        end if
        line(np) = gsn_csm_pres_hgt(wks,var(nc,0,nm,:,:),res_line)
        overlay(plot(np),line(np))
        vect(np) = gsn_csm_pres_hgt_vector(wks,var(nc,4,nm,:,:),var(nc,1,nm,:,:),var(nc,2,nm,:,:),resvc)
        overlay(plot(np),vect(np))
        do nv = 0, nvar-1, 1
        print(draw_var(nv))
        printMinMax(var(nc,nv,nm,:,:),0)
        end do
    end do 
end do ;for the month 

;    gsn_panel(wks,plot,(/1,2/),resp)
;    gsn_panel(wks,plot,(/ncase,nmonth/),resp)
    gsn_panel(wks,plot,(/nmonth,ncase/),resp)
;end do
end

