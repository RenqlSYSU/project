;***************************************************************************
;read cmip data and draw contour figure (one figure in each page)
;                                            by Ql_Ren
;                                           2020/08/17
;******************************************************************************
begin
lats = (/20 ,20/);-20
latn = (/27 ,27/);60
lonl = (/60 ,85/);20
lonr = (/120,95/);220
lev  = (/850,700,500,300,200/)
nl = 1

path1 = "/home/ys17-19/data-model/CMIP5/"
path2 = "/home/ys17-19/renql/project/2020SouthTrough/"
var_name = (/"zg"/)
unit     = (/"m" /)
case  = "CMIP5"
nv = 0

months   = (/"DJF","MAM","JJA","SON"/)
fig_name = path2 + "/fig/" + case(0) + "-DJF_ts_" + lev(nl) + "_" + var_name(nv)
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
   tick_thick     = 8.0
   subtitle_thick = 8.0  ;not have this attribute
   title_thick    = 5.0
   mp_thick       = 2.5;3.0;
   font_index     = 22

wks = gsn_open_wks(fig_out,fig_name)

    gsn_define_colormap(wks,"matlab_jet") 
    ;has 64 colors blue green yellow,red,2-65

    resxy     = True
    resxy@vpWidthF  = 0.8
    resxy@vpHeightF = 0.4
    resxy@gsnFrame      = False
    resxy@gsnDraw       = False
    resxy@gsnMaximize   = True
    resxy@gsnPaperOrientation = "portrait"
    
    resxy@xyDashPatterns    = 0 ;(/0,0,0/)
    resxy@xyLineThicknesses = 2.0;(/2.0, 2.0,2.0/)   ; make second line thicker
    resxy@xyLineColors      = ispan(4,41,1)
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
    resxy@tmXBValues          = ispan(1,12,1)
    resxy@tmXBLabels          = ispan(1,12,1)
    
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
    resxy2@xyLineThicknessF = 6.0;(/2.0, 2.0,2.0/)   ; make second line thicker
    resxy2@xyLineColor      = "black"

    resp  = True    
    resp@gsnMaximize   = True
    resp@gsnPaperOrientation = "landscape"

;=================================================================
;read cmip data
;=====================================================================
varm  = new((/38,12/),float)
varo  = new((/12/),float)
mname = new((/38/),string)

f     = addfile("~/renql/project/2020MiddleEastJet/data/ERA-Interim_hgt-year.nc","r")
term := dim_avg_n_Wrap(f->month_ave(:,:,{lev(nl)},{lats(0):latn(0)},{lonl(0):lonr(0)}),0)
varo  = wgt_areaave(term(:,:,{lonl(1):lonr(1)}),1.0,1.0,0) - wgt_areaave(term(:,:,{lonl(0):lonr(0)}),1.0,1.0,0)
varo  = varo/9.8
delete([/f,term/])
printMinMax(varo,0) 

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
nmodel = -1
print("file number " + nfile)
print(" ")
do nf = 0, nfile-1, 1
if((nf.le.nfile-2).and.(modelname(nf).eq.modelname(nf+1))) then 
    numb = numb + 1
else
    numb = numb + 1
    nmodel = nmodel+1
    print( nmodel  + " model name : " + modelname(nf) + "  filenumber : " + numb)
    filepaths  := filepath((nf-numb+1):nf)
    yearindex  := dim_pqsort(year((nf-numb+1):nf),1)
    fs := addfiles(filepaths(yearindex),"r")
    ListSetType(fs,"cat")
    time := cd_calendar(fs[:]->time,0)
    
    do nm = 0, 11, 1
    timeindex := ind(time(:,1).eq.(nm+1))
    term := dim_avg_n_Wrap(fs[:]->$var_name(nv)$(timeindex,{lev(nl)*100},{lats(0):latn(0)},{lonl(0):lonr(0)}),0)
    varm(nmodel,nm) = wgt_areaave(term(:,{lonl(1):lonr(1)}),1.0,1.0,0) - wgt_areaave(term(:,{lonl(0):lonr(0)}),1.0,1.0,0)
    end do

    dnames := getfilevardims(fs[0],var_name(nv))  
    dsizes := getfilevardimsizes(fs[0],var_name(nv))
    mname(nmodel) = modelname(nf) + " " + str_join(""+dsizes(2:3),"x")
    delete([/fs,filepaths,yearindex,time/])
    numb  = 0   
end if
end do
    printMinMax(varm,0) 

    resxy@gsnLeftString   = pre_case(0) + ") (" + lonl(1) + "-" + lonr(1) + "E," + lats(1) + "-" + latn(1) \
                                     + "N) - (" + lonl(0) + "-" + lonr(0) + "E," + lats(0) + "-" + latn(0) + "N)"
    resxy@gsnRightString  = lev(nl) + "hPa " + var_name(nv) + " (" + unit(nv) + ")"  
    resxy@trYMinF = floor(min(varm))
    resxy@trYMaxF =  ceil(max(varm))
    plot  = gsn_csm_xy(wks,ispan(1,12,1),(/varm/),resxy)
    plot2 = gsn_csm_xy(wks,ispan(1,12,1),(/varo/),resxy2)
    overlay(plot,plot2)
    gsn_panel(wks,plot,(/1,1/),resp)
end
    
