;***************************************************************************
;read cmip data and draw contour figure (one figure in each page)
;                                            by Ql_Ren
;                                           2020/08/17
;******************************************************************************
begin
lats = (/10 /);-20
latn = (/50 /);60
lonl = (/40 /);20
lonr = (/130/);220
lat_spacing = 15
lon_spacing = 20
lev  = (/850,700,500,300,200/)
nl = 1
na = 0

path1 = "/home/ys17-19/data-model/CMIP5/"
path2 = "/home/ys17-19/renql/project/2020SouthTrough/"
var_name = (/"zg"/)
unit     = (/"m" /)
case  = "CMIP5"
nv = 0

months   = (/"DJF","MAM","JJA","SON"/)
fig_name = path2 + "/fig/" + case(0) + "-DJF_" + lev(nl) + "_" + var_name(nv)
pre_case = (/"(a","(b","(c","(d","(e","(f","(g","(h","(i","(j","(k","(l"/)
fig_out  = "pdf";"x11";"eps";
barshare = False;True ;
barvert  = False;True

;==================================================
;figure setting
;==================================================
   tick_font      = 0.02
   subtitle_font  = 0.02
   title_font     = 0.02
   vcref_font     = 0.02
   tick_thick     = 8.0
   subtitle_thick = 8.0  ;not have this attribute
   title_thick    = 5.0
   vc_thick       = 2.0
   mp_thick       = 2.5;3.0;
   font_index     = 22

    res                  = True
    res@gsnFrame         = False
    res@gsnDraw          = False
    res@gsnMaximize      = True
    res@gsnAddCyclic     = False;True;
    res@gsnPaperOrientation = "landscape" ;"portrait"
    ;res@vpHeightF        = 0.25
    ;res@vpWidthF         = 0.8

    res@mpCenterLonF     = 180.
    res@mpMaxLatF        = latn(na)
    res@mpMinLatF        = lats(na)
    res@mpMaxLonF        = lonr(na)
    res@mpMinLonF        = lonl(na)
    res@gsnMajorLatSpacing = lat_spacing
    res@gsnMajorLonSpacing = lon_spacing
    res@mpFillOn         = False
    ;res@mpLandFillColor = "gray"
    ;res@mpFillDrawOrder = "PreDraw"
    res@mpGeophysicalLineThicknessF = mp_thick
    res@mpGeophysicalLineColor      = "black";"white"

;--------set the form of the coordinate-----------------
    res@tmBorderThicknessF  = mp_thick
    res@tmXBLabelFont   = font_index
    res@tmYLLabelFont   = font_index
    res@tmXBLabelFontHeightF     = tick_font
    res@tmXBLabelFontThicknessF  = tick_thick
    res@tmYLLabelFontHeightF     = tick_font
    res@tmYLLabelFontThicknessF  = tick_thick
    
    res@tmXBMajorThicknessF     = mp_thick
    res@tmYLMajorThicknessF     = mp_thick
    res@tmXBMinorThicknessF     = mp_thick
    res@tmYLMinorThicknessF     = mp_thick

;---------set the attribute of the contour-----------------
    res@cnFillOn         = True
    res@cnLinesOn        = False ;True
    res@cnLineLabelsOn   = False
    res@cnInfoLabelOn    = False
    
    res@gsnSpreadColors      = False  ;True              ; turn of the full color map
    res@gsnSpreadColorStart  = 11                ; start at color 17
    res@gsnSpreadColorEnd    = 2                ; end at color 200

    res@gsnCenterString      = ""
    res@gsnStringFont         = font_index
    res@gsnStringFontHeightF = subtitle_font

;--------set the attribute of gtop--------------
    res_gtop                  = True
    res_gtop@gsnDraw          = False
    res_gtop@gsnFrame         = False
    res_gtop@cnFillOn         = False
    res_gtop@cnLinesOn        = True      ; turn off contour lines
    res_gtop@cnLineLabelsOn   = False      ; turn off contour line labels
    res_gtop@cnInfoLabelOn    = False
    res_gtop@cnLineThicknessF = mp_thick + 0.5
    res_gtop@cnLineColor      = "black"
    res_gtop@cnLevelSelectionMode  = "ExplicitLevels" ; set manual contour levels
    res_gtop@cnLevels         = (/1500/)
    res_gtop@gsnCenterString      = ""
    res_gtop@gsnLeftString    = " "
    res_gtop@gsnRightString   = " "

    resp  = True    
    resp@gsnMaximize   = True
    if(barshare) then; .and.nmonth.eq.4
    res@lbLabelBarOn         = False   ; the label bar
    resp@gsnPaperOrientation = "landscaoe";"portrait"
    resp@gsnPanelLabelBar = True
    if(barvert) then 
    resp@lbOrientation     = "Vertical"
    resp@pmLabelBarWidthF  = 0.05
    resp@pmLabelBarHeightF = 0.8
    else
    resp@lbOrientation     = "Horizontal" ;"Vertical"
    resp@pmLabelBarWidthF  = 0.8 ;0.05
    resp@pmLabelBarHeightF = 0.05;0.8
    end if
    else
    res@lbLabelBarOn         = True   ; the label bar
    res@lbLabelFontHeightF   = tick_font
    res@lbLabelFont          = font_index              ; make labels smaller
    res@lbLabelStride        = 1
    res@lbBoxLineThicknessF  = mp_thick-1
    res@pmLabelBarZone               = 0
    if(barvert) then
    res@lbOrientation        = "Vertical";
    res@pmLabelBarWidthF  = 0.03
    res@pmLabelBarHeightF = 0.17
    res@pmLabelBarOrthogonalPosF     = 0.01
    res@pmLabelBarParallelPosF       = 0.5
    else
    res@lbOrientation        = "Horizontal";"Vertical";
    res@pmLabelBarWidthF  = 0.7
    res@pmLabelBarHeightF = 0.1
    res@pmLabelBarOrthogonalPosF     = 0.25 ;y direction
    res@pmLabelBarParallelPosF       = 0.5
    end if
    resp@gsnPaperOrientation = "landscape" ;"portrait"
    resp@gsnPanelLabelBar = False; True
    end if
    resp@lbLabelFontHeightF = 0.015 
    resp@lbLabelFontThicknessF = tick_thick
    resp@txFontHeightF = title_font 
    resp@txFontThicknessF = title_thick

gtopfile = "~/renql/project/TP_NUDG/analysis/mdata/gtopo30_0.9x1.25.nc"
f = addfile(gtopfile,"r")
phis = f->PHIS({lats(na):latn(na)},{lonl(na):lonr(na)})
phis = phis/9.8   ;convert unit from m2/s2 to gpm
delete(f)

wks = gsn_open_wks(fig_out,fig_name)
;=================================================================
;read cmip data
;=====================================================================
filepath  = systemfunc("find " + path1 + " -name " + var_name(nv) + "_Amon_*_historical_r1i1p1_*.nc")
modelname = str_get_field(filepath,5,"/")
nfile     = dimsizes(filepath)
nfield    = str_fields_count(filepath,"_")
numbyear  = new((/nfile/),string)
do nf = 0, nfile-1, 1
    numbyear(nf) = str_get_field(filepath(nf),nfield(nf),"_")
end do
year = stringtoint(str_get_cols(numbyear,0,3))

numb   = 0
nmodel = 0
print("file number " + nfile)
print(" ")
do nf = 0, nfile-1, 1
if((nf.le.nfile-2).and.(modelname(nf).eq.modelname(nf+1))) then 
    numb = numb + 1
else
    numb = numb + 1
    nmodel = nmodel+1
    filepaths  := filepath((nf-numb+1):nf)
    yearindex  := dim_pqsort(year((nf-numb+1):nf),1)
    fs := addfiles(filepaths(yearindex),"r")
    ListSetType(fs,"cat")
    time := cd_calendar(fs[:]->time,0)
    timeindex := ind(time(:,1).le.2.or.time(:,1).eq.12)
    print("month number" + dimsizes(timeindex))
    var := dim_avg_n_Wrap(fs[:]->$var_name(nv)$(timeindex,{lev(nl)*100},{lats(na):latn(na)},{lonl(na):lonr(na)}),0)
    dnames := getfilevardims(fs[0],var_name(nv))  
    dsizes := getfilevardimsizes(fs[0],var_name(nv))
    delete([/fs,filepaths,yearindex,time/])

;    print("")
;    print("=============================================================")
    print( nmodel  + " model name : " + modelname(nf) + "  filenumber : " + numb)
;    print("  dimension names : " + str_join(dnames,","))
;    print("  dimension sizes : " + str_join(""+dsizes,","))
;    print("")
    
    if(nv.eq.0) then ;zg, geopotential height
        res@cnFillPalette = "precip3_16lev"
        res@cnFillColors  := (/0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16/)
        minlvl = (/-65,2750,-3.5/); -24,-24
        maxlvl = (/-50,3200, 3.5/);  24, 24
        spacig = (/1  , 30 ,0.5 /);   3,  3
    end if

;    res@gsnLeftString   = pre_case(0) + ") " + months(0) + "  " + modelname(nf) + " " + str_join(dnames(2:3),"x") + ":" + str_join(""+dsizes(2:3),"x") 
    res@gsnLeftString   = modelname(nf) + " " + str_join(dnames(2:3),"x") + ":" + str_join(""+dsizes(2:3),"x") 
    res@gsnRightString  = lev(nl) + "hPa " + var_name(nv) + " (" + unit(nv) + ") " 
    res@cnLevelSelectionMode  = "ManualLevels" ; set manual contour levels
    res@cnLevelSpacingF       := spacig(nl)        ; set contour interval
    res@cnMaxLevelValF        := maxlvl(nl) ;32    ; set max contour level
    res@cnMinLevelValF        := minlvl(nl)       ; set min contour level
    plot = gsn_csm_contour_map_ce(wks, var, res)
    gtop = gsn_csm_contour(wks, phis, res_gtop)
    overlay(plot,gtop)
    printMinMax(var,0)
    gsn_panel(wks,plot,(/1,1/),resp)

    numb  = 0   
end if
end do

end
    
