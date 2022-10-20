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

draw_var = (/"zg","vort"/)
var_name = (/"zg","ua"  /)
unit     = (/"m" ,"10~S~-5~N~s~S~-1~N~"/)
case  = (/"CMIP5","AMIP"/)
nv = 1
nc = 0

path3 = "~/renql/project/2020SouthTrough"
if(nc.eq.1) then 
modelnumb = 19
path1 = "/home/ys17-19/data-model/AMIP/ -name " + var_name(nv) + "_Amon_*_r1i1p1_*.nc"
end if
if(nc.eq.0) then 
modelnumb = 32
path1 = "/home/ys17-19/data-model/CMIP5/ -name " + var_name(nv) + "_Amon_*_historical_r1i1p1_*.nc"
path2 = "/home/ys17-19/data-model/CMIP5/ -name va_Amon_*_historical_r1i1p1_*.nc"
end if
months   = (/"DJF","MAM","JJA","SON"/)
fig_name = path3 + "/fig/" + case(nc) + "-DJF_ts_" + lev(nl) + "_" + draw_var(nv)
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

    ;gsn_define_colormap(wks,"matlab_jet") 
    ;has 64 colors blue green yellow,red,2-65
    gsn_define_colormap(wks,"precip3_16lev")
    ;blue and yellow and red, has 17 colors,2-18,2 is white
    ;gsn_define_colormap(wks,"BlueDarkRed18")
    ;blue and red, has 18 colors,2-10 is blue

    resxy     = True
    resxy@vpWidthF  = 0.8
    resxy@vpHeightF = 0.4
    resxy@gsnFrame      = False
    resxy@gsnDraw       = False
    resxy@gsnMaximize   = True
    resxy@gsnPaperOrientation = "portrait"
    
    resxy@xyMonoLineThickness = True
    resxy@xyDashPatterns    = 0 ;(/0,0,0/)
    resxy@xyLineThicknesses = 5.0;(/2.0, 2.0,2.0/)   ; make second line thicker
    ;resxy@xyLineColors      = ispan(4,41,1)
    resxy@tiXAxisOn = False
    resxy@tiYAxisOn = False
    ;resxy@trYMinF = floor(min(var))
    ;resxy@trYMaxF = ceil(max(var))
    resxy@trXMinF    = 1
    resxy@trXMaxF    = 12
    resxy@tmXBMode   = "Explicit"
    resxy@tmXBValues = ispan(1,12,1)
    resxy@tmXBLabels = ispan(1,12,1)
    
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
    resp@gsnFrame      = False  ; Don't advance frame.

;=================================================================
;read cmip data
;=====================================================================
varm  = new((/modelnumb,12/),float)
varo  = new((/12/),float)
mname  = new((/modelnumb/),string)
colors = new((/modelnumb/),string)
nlat   = new((/modelnumb/),string)
nlon   = new((/modelnumb/),string)

if(nv.eq.0) then 
f     = addfile("~/renql/project/2020MiddleEastJet/data/ERA-Interim_hgt-year.nc","r")
term := dim_avg_n_Wrap(f->month_ave(:,:,{lev(nl)},{lats(0):latn(0)},{lonl(0):lonr(0)}),0)
varo  = wgt_areaave(term(:,:,{lonl(1):lonr(1)}),1.0,1.0,0) - wgt_areaave(term(:,:,{lonl(0):lonr(0)}),1.0,1.0,0)
varo  = varo/9.8
end if
if(nv.eq.1) then 
f     = addfile("~/renql/project/2020MiddleEastJet/data/ERA-Interim_sfvp_global-year.nc","r")
term := dim_avg_n_Wrap(f->vort_m(:,:,{lev(nl)},{lats(0):latn(0)},{lonl(0):lonr(0)}),0)
varo  = wgt_areaave(term(:,:,{lonl(1):lonr(1)}),1.0,1.0,0)
varo  = varo*100000
end if
delete([/f,term/])
printMinMax(varo,0) 

filepath  = systemfunc("find " + path1)
filepath2 = systemfunc("find " + path2)
nfile     = dimsizes(filepath)
nfield    = str_fields_count(filepath,"/")
modelname = new((/nfile/),string)
year      = new((/nfile/),integer)
do nf = 0, nfile-1, 1
    filename = str_split(str_get_field(filepath(nf),nfield(nf),"/"),"_")
    modelname(nf) = filename(2)
    year(nf) = stringtoint(str_get_cols(filename(5),0,3))
end do
print(year(0:5))

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
    if(nv.eq.0) then 
    term := dim_avg_n_Wrap(fs[:]->$var_name(nv)$(timeindex,{lev(nl)*100},{lats(0):latn(0)},{lonl(0):lonr(0)}),0)
    varm(nmodel,nm) = wgt_areaave(term(:,{lonl(1):lonr(1)}),1.0,1.0,0) - wgt_areaave(term(:,{lonl(0):lonr(0)}),1.0,1.0,0)
    else
    filepaths2  := filepath2((nf-numb+1):nf)
    fs2 := addfiles(filepaths2(yearindex),"r")
    ListSetType(fs2,"cat")
    ua := dim_avg_n_Wrap( fs[:]->$var_name(nv)$(timeindex,{lev(nl)*100},{lats(0):latn(0)},{lonl(0):lonr(0)}),0)
    va := dim_avg_n_Wrap(fs2[:]->va(timeindex,{lev(nl)*100},{lats(0):latn(0)},{lonl(0):lonr(0)}),0)
    term  := uv2vr_cfd ( ua, va, ua&lat, ua&lon, 2 )*100000
    copy_VarMeta(ua,term)
    varm(nmodel,nm) = wgt_areaave(term(:,{lonl(1):lonr(1)}),1.0,1.0,0)
    end if
    end do

    dnames := getfilevardims(fs[0],var_name(nv))  
    dsizes := getfilevardimsizes(fs[0],var_name(nv))
    mname(nmodel) = modelname(nf) + " " + str_join(""+dsizes(2:3),"x")
    nlat(nmodel) = dsizes(2)
    nlon(nmodel) = dsizes(3)
    delete([/fs,filepaths,yearindex,time/])
    numb  = 0   
end if
end do

;dim_lat = (/60,64,90, 96,108,120, 128,143,144, 160,192,320, 720/) ;AMIP
;dim_lon = (/128,144,180, 192,256,288, 320,640,1152/) ;9 AMIP
dim_lat = (/48,60,64,73,89,90,96,108,120,128,143,144,160,192,320/) ;38 CMIP5 zg
dim_lon = (/96,128,144, 180,192,256, 288,320,640/) ;CMIP5
ndimlat = dimsizes(dim_lat)
ndimlon = dimsizes(dim_lon)
numblon = new((/ndimlon/),integer)
colorindex = (/4,7,8, 9,10,11, 15,17,18/) 
;colorindex = (/"blue","red","green","grey","purple", "Pink"/) 
do nla = 0, ndimlon-1, 1
    colors = where(nlon.eq.dim_lon(nla),colorindex(nla),colors)
    numblon(nla) = dimsizes(ind(nlon.eq.dim_lon(nla)))
end do
    lgres                    = True
    lgres@vpWidthF           = 0.1     ; width of legend (NDC)
    lgres@vpHeightF          = 0.09     ; height of one legend (NDC)
;    lgres@lgBoxMajorExtentF  = 0.02          ; puts space between color boxes
    lgres@lgPerimOn          = False                 ; turn off perimeter
;    lgres@lgLineColors       = colorindex 
    lgres@lgMonoDashIndex    = True
    lgres@lgDashIndexes      = 0 ; resxy@xyDashPatterns
    lgres@lgLineThicknesses  = 10.0
    lgres@lgLineLabelStrings = ""
    lgres@lgLabelFontHeightF = 0.2 
    lgres@lgLabelFont        = font_index

    printMinMax(varm,0) 
    resxy@xyLineColors    = colors 
    resxy@gsnLeftString   = pre_case(0) + ") (" + lonl(1) + "-" + lonr(1) + "E," + lats(1) + "-" + latn(1) \
                                     + "N) - (" + lonl(0) + "-" + lonr(0) + "E," + lats(0) + "-" + latn(0) + "N)"
    resxy@gsnRightString  = lev(nl) + "hPa " + var_name(nv) + " (" + unit(nv) + ")"  
    resxy@trYMinF = floor(min(varm))
    resxy@trYMaxF =  ceil(max(varm))
    plot  = gsn_csm_xy(wks,ispan(1,12,1),(/varm/),resxy)
    plot2 = gsn_csm_xy(wks,ispan(1,12,1),(/varo/),resxy2)
    overlay(plot,plot2)
    gsn_panel(wks,plot,(/1,1/),resp)
    printMinMax(varo,0) 

xpos = (/0.2, 0.4, 0.6/)
nnn = 3
do nn = 0,2,1
    nns = nnn*nn
    nne = nnn*nn+2
    lgres@lgLineColors = colorindex(nns:nne) 
    gsn_legend_ndc(wks, nnn, numblon(nns:nne)+" nlon="+dim_lon(nns:nne),xpos(nn),0.2,lgres)
end do
frame(wks)
end
    
