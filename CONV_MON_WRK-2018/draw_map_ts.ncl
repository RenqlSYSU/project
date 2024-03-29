;***************************************************************************
;filename:      draw_fsnt_flnt_preci_global_ts.ncl
;1 read data in (ncase,12nmonth*100years,nlat,nlon),FSNT+FLNT
;2 calc global weighted annual mean，then draw the ts plot
;                                            by Ql_Ren
;                                           2018/07/19
;******************************************************************************
load "$NCARG_ROOT/lib/ncarg/nclscripts/contrib/time_axis_labels.ncl"

begin
unit     = (/"~S~o~N~C","hPa","W/m2","W/m2"/)    ;this can be used as right string 
var_name = (/"TS","PSL","FLNT","FSNT"/)
case     = (/"MX_AQUAP_SCYC","SF_AQUAP_SCYC","NO_TOPO_SCYC","REAL_WORLD_SCYC"/);" NG24H",THIS CAN BE USED AS LEFT STRING
pre_case = (/"(a","(b","(c","(d"/)
nv = 2

;the path of the file
path     = "/home/ys17-19/renql/project/CONV_MON_WRK-2018" 
filename = path + "/data/" + case + "/pro/B20f19-" + (/"pure-aqua","surf-aqua","topo","realworld"/) + ".cam.h0." ; ; + var_name + ".1-100.nc" 
fig_name = path + "/fig/animation/" + var_name(nv) + "_global_map"
fig_out  = "pdf";"x11"

f    = addfile(filename(0) + var_name(0) + ".1-100.nc","r")
vars = f->$var_name(0)$(0,:,:)
nlat   = dimsizes(vars&lat) 
nlon   = dimsizes(vars&lon)
ncase  = dimsizes(case)
nyear  = 100 
nmonth = 12
months   = (/"Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec"/)
var   = new((/ncase,100*12,nlat,nlon/),float)
var1  = new((/ncase,nyear ,nlat,nlon/),float)
plot  = new((/ncase/),graphic)

;================================================================
;read data
;====================================================================
var = 0
do nc = 0,ncase-1,1
    print(filename(nc) + var_name(nv))
    f    = addfile(filename(nc) + var_name(nv) + ".1-100.nc","r")
    var(nc,:,:,:) = f->$var_name(nv)$
end do

if(nv.eq.0) then
var = var - 273.15 ;convert unit from K to C
end if
if(nv.eq.1) then
var = var/100 ;convert unit from Pa to hPa
end if
printVarSummary(var)

;do ns = 0,nyear-1,1
;    nms = ns*nmonth
;    nme = ns*nmonth + 11 
;    var1(:,ns,:,:) = dim_avg_n(var(:,nms:nme,:,:),1) ;calc the annual mean
;end do
;var1 := dim_avg_n(var1(:,69:99,:,:),1)
;copy_VarMeta(vars,var1(0,:,:))
;printVarSummary(var1)

;===============================================
;draw the figure
;=============================================
;define the thickness and height of the title,tick and label
   tick_font      = 0.015
   subtitle_font  = 0.018
   title_font     = 0.03
   label_font     = 0.015
   tick_thick     = 5.0
   subtitle_thick = 5.0  ;not have this attribute
   title_thick    = 5.0
   label_thick    = 1.5
   mp_thick       = 3.0
   font_index     = 22

wks = gsn_open_wks(fig_out,fig_name)

    ;gsn_define_colormap(wks,"BlueDarkRed18") 
    ;;blue and red, has 18 colors,2-10 is blue
    ;gsn_define_colormap(wks,"precip3_16lev") 
    ;blue and yellow and red, has 17 colors,2-18,2 is white
    ;gsn_define_colormap(wks,"GMT_polar") 
    ;blue and red, has 20 colors, 2-11 is bule, 12-21 is red
    gsn_define_colormap(wks,"matlab_jet") 
    ;has 64 colors blue green yellow,red

    res                  = True
    res@gsnFrame         = False
    res@gsnDraw          = False
    res@gsnMaximize      = True
    res@gsnAddCyclic     = False
    res@gsnPaperOrientation = "portrait"
    ;res@vpHeightF        = 0.25
    ;res@vpWidthF         = 0.8

;-------set the attribute of the map-------------------
    res@mpProjection      = "Mollweide"       ; choose projection
    res@mpGridAndLimbOn   = True              ; turn on lat/lon lines
    res@mpPerimOn         = False             ; turn off box around plot
    res@mpGridLatSpacingF = 30.               ; spacing for lat lines
    res@mpGridLonSpacingF = 30.               ; spacing for lon lines
    ;res@mpCenterLonF     = 180.
    ;res@mpMaxLatF        = latn
    ;res@mpMinLatF        = lats
    ;res@mpMaxLonF        = lonr
    ;res@mpMinLonF        = lonl
    ;res@gsnMajorLatSpacing = 15
    ;res@gsnMajorLonSpacing = 30
    res@mpFillOn         = False
    ;res@mpLandFillColor = "gray"
    ;res@mpFillDrawOrder = "PreDraw"
    res@mpGeophysicalLineThicknessF = mp_thick
    res@mpGeophysicalLineColor      = "black";"white"

;---------set the attribute of the contour-----------------
    res@cnFillOn         = True
    res@cnLinesOn        = False ;True
    res@cnLineLabelsOn   = False
    res@cnInfoLabelOn    = False
    ;res@cnLineThicknessF = 0
    ;res@gsnContourZeroLineThicknessF = 4.0
    ;res@cnLineColor      = "blue"
   
    res@cnLevelSelectionMode = "AutomaticLevels"   ; set manual contour levels
    ;res@cnLevelSelectionMode = "ExplicitLevels"   ; set manual contour levels
    ;res@cnFillColors         = ispan(2,65,2)
    ;res@cnLevelSelectionMode  = "ManualLevels" ; set manual contour levels
    ;res@cnMinLevelValF        = 0          ; set min contour level
    ;res@cnMaxLevelValF        = 30 ;32        ; set max contour level
    ;res@cnLevelSpacingF       = 1        ; set contour interval
    ;
    res@gsnSpreadColors      = False  ;True              ; turn of the full color map
    res@gsnSpreadColorStart  = 11                ; start at color 17
    res@gsnSpreadColorEnd    = 2                ; end at color 200
    
    res@lbLabelBarOn         = False   ; the label bar
    res@lbOrientation        = "Horizontal"
    res@lbLabelFontHeightF   = tick_font              ; make labels smaller
    res@lbLabelFont          = font_index              ; make labels smaller
    res@lbLabelStride        = 1
    res@lbBoxLineThicknessF  = mp_thick-2

    res@pmLabelBarZone               = 0
    res@pmLabelBarOrthogonalPosF     = 0.15
    res@pmLabelBarParallelPosF       = 0.5

;------------------------------------------------------------------
;draw the subgraph and change the left and right titles
;-----------------------------------------------------------------------
    res@gsnCenterString      = ""
    res@gsnStringFont         = font_index
    res@gsnStringFontHeightF  = subtitle_font

    resp  = True    
    resp@gsnMaximize   = True
    resp@gsnPaperOrientation = "portrait" ;"landscape"
    resp@gsnPanelLabelBar = True
    resp@pmLabelBarWidthF  = 0.1
    resp@pmLabelBarHeightF = 0.8
    resp@lbLabelFontHeightF = tick_font
    resp@lbLabelFontThicknessF = tick_thick
    resp@lbOrientation = "Vertical"
    resp@txFontHeightF = title_font 
    resp@txFontThicknessF = title_thick

    minlvl = (/ 0,935 /)
    maxlvl = (/30,1025/)
    spacig = (/ 1,3   /)

do nt = 0,11,1;100*12-1,1 
ny = floor((nt)/12)+1
print(ny)
nm = toint(nt-(ny-1)*12)
print(nm)
do nc = 0, ncase-1,1
    res@gsnLeftString   = pre_case(nc) + ")" + case(nc)
    res@gsnRightString  = var_name(nv) + " (" + unit(nv) + ") " + months(nm) ;lev(nl) + "hPa " +  
    ;res@cnLevelSelectionMode  = "ManualLevels" ; set manual contour levels
    ;res@cnLevelSpacingF       := spacig(nv)        ; set contour interval
    ;res@cnMaxLevelValF        := maxlvl(nv) ;32    ; set max contour level
    ;res@cnMinLevelValF        := minlvl(nv)       ; set min contour level
    plot(nc)   = gsn_csm_contour_map(wks,var(nc,nt,:,:),res)
    print(max(var(nc,nt,:,:)))
    print(min(var(nc,nt,:,:)))
end do
    resp@txString      = "Year" + ny
    gsn_panel(wks,plot,(/ncase,1/),resp)
end do
end

