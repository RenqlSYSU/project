;----------------------------------------------------------------------
;calc_std_11slid_corr_JJApreci_DJFnino3.ncl  20180914
;----------------------------------------------------------------------
; Concepts illustrated:
;   - Reading an ASCII station file data with several columns of data
;   - Calculate the average of JJA rainfall and DJF Nino3 index
;   - Calculate the 11years sliding correlation between above data
;   - calculate standard deviation of the sliding correlation
;   - Interplote standard deviation to grid data and plot it 
;----------------------------------------------------------------------
; It uses a newer method (NCL V6.4.0 and later) for specifying the
; lat/lon information via special lat1d/lon1d attributes attached
; to the data to be plotted. 
;----------------------------------------------------------------------

begin
  path      = "/home/ys17-19/renql/project/2018homework/statistic"
  preci_f   = path + "/data/160preci/preci0" + (/"6","7","8"/) + ".txt"
  nino3_f   = path + "/data/ersst5.nino.mth.81-10.ascii"
  station_f = path + "/data/ncc_observation_readme.txt"
  fig_name  = path + "/std_11slid_corr_JJApreci_DJFnino3"
  fig_out   = "pdf" 
  fig_title = "Standard Deviation"

  nstation  = 160
  years     = ispan(1951,2017,1)
  nyear     = dimsizes(years)
  nmonth    = 3
  window    = 11  ; the slide window
  ncorr     = nyear-window+1

  preci     = new((/nmonth,nstation,nyear/),float)
  nino3     = new((/nyear ,nmonth/),float)
  lat       = new((/nstation/),float)
  lon       = new((/nstation/),float)
  slid_corr = new((/ncorr,nstation/),float)

;----- Read JJA precipitation from ascii file which are all numeric--------------------
; Preci data is stored in a file with two dimensions (station160,year67),from 1951 to 2017
do nm = 0,nmonth-1,1
  preci(nm,:,:) = tofloat(asciiread(preci_f(nm),(/nstation,nyear/),"integer"))
end do

;-----Read DJF Nino3 SST from ascii file which have character----------------------------
; Nino3 data is stored in 10 columns: year(1950-2017), month(12), Nino1+2, Anom, 
;              NINO3, ANOM, NINO4, ANOM, NINO3.4, ANOM
; For string data, each line of the file is read as a string.
  line  = asciiread(nino3_f, -1, "string")
  printVarSummary(line)
; Use "str_get_field" to indicate which fields to read in.
; Each field is separated by an arbitrary number of spaces.
do nm = 0,nmonth-1,1
  nino3(:,nm) = tofloat(str_get_field(line((12+nm):815:12),5," ")) ;only read DJF
end do

;-----Read the information of 160 stations from ascii file which have header lines-------
; Station information is stored in five columns: number,name,number,lat,lon
  line := asciiread(station_f, -1, "string")
  printVarSummary(line)
  lat = tofloat(str_get_field(line(21:),1," "))  ; Attach as lat1d, lon1d
  lon = tofloat(str_get_field(line(21:),2," "))  ; for plotting later
  print(lat(1))
; if there are number, name, number in file, can not separate each field successfully
; so remove these three columns

;-----calculate the sliding correlation and standard deviation-------------
  preci := dim_avg_n_Wrap(preci,0) ;JJA average, get(nstation, nyear)
  nino3 := dim_avg_n_Wrap(nino3,1) ;DJF average, get(nyear)
do nn = 0,ncorr-1,1
  slid_corr(nn,:) = escorc_n(preci(:,nn:(nn+10)), nino3(nn:(nn+10)), 1, 0)
end do
  slid_corr!0 = "year"
  slid_corr&year = years(5:(nyear-6))
  printVarSummary(slid_corr)
  std  = dim_stddev_n_Wrap(slid_corr({1956:2002},:),0) ;get (nstation)
  std@lat1d = lat
  std@lon1d = lon
  printVarSummary(std)

;==============================================================
;draw figure 
;==================================================================
;define the thickness and height of the title,tick and label
   tick_font      = 0.02
   subtitle_font  = 0.02
   title_font     = 0.03
   tick_thick     = 5.0
   subtitle_thick = 5.0  ;not have this attribute
   title_thick    = 5.0
   mp_thick       = 3.0
   font_index     = 22

;draw the figure
    wks = gsn_open_wks(fig_out,fig_name)

    cmap = read_colormap_file("GMT_gray") ; read color map
    cmap = cmap(::-1,:)                   ; reverse the color map
    ;gsn_define_colormap(wks,"BlueDarkRed18")
    ;blue and red, has 18 colors,2-10 is blue
    ;blue and red, has 20 colors, 2-11 is bule, 12-21 is red
    ;gsn_define_colormap(wks,"matlab_jet") 
    ;has 64 colors blue green yellow,red

    res                  = True
    ;res@gsnFrame         = False
    ;res@gsnDraw          = False
    res@gsnMaximize      = True
    res@gsnAddCyclic     = False
    res@gsnPaperOrientation = "landscape"
    ;res@vpHeightF        = 0.25
    ;res@vpWidthF         = 0.8

;-------set the attribute of the map-------------------
    res@mpCenterLonF     = 180.
    res@mpMaxLatF        = 55
    res@mpMinLatF        = 20
    res@mpMaxLonF        = 140
    res@mpMinLonF        = 70
    res@gsnMajorLatSpacing = 10
    res@gsnMajorLonSpacing = 20
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
    ; res@tmXBLabelFontAspectF   = 1.6
    res@tmXBLabelFontThicknessF  = tick_thick
    res@tmYLLabelFontHeightF     = tick_font
    ; res@tmYLLabelFontAspectF   = 1.3
    res@tmYLLabelFontThicknessF  = tick_thick
    
    ;  res@tmXBMajorOutwardLengthF = -0.
    ;  res@tmYLMajorOutwardLengthF = -0.
    ;  res@tmYRMajorOutwardLengthF = -0.
      
    ;  res@tmXBMinorOutwardLengthF = -0.
    ;  res@tmYLMinorOutwardLengthF = -0.
    ;  res@tmYRMinorOutwardLengthF = -0.

;---------set the attribute of the contour-----------------
    res@cnFillOn         = True
    res@cnLinesOn        = True
    res@cnLineLabelsOn   = True
    res@cnInfoLabelOn    = False
    res@cnLineThicknessF = 5.0
    ;res@gsnContourZeroLineThicknessF = 4.0
    res@cnLineColor      = "black"
   
    res@cnLevelSelectionMode  = "ManualLevels" ; set manual contour levels
    res@cnMinLevelValF        = 0         ; set min contour level
    res@cnMaxLevelValF        = 0.8        ; set max contour level
    res@cnLevelSpacingF       = 0.1       ; set contour interval
    res@cnFillPalette         = cmap

    res@gsnSpreadColors      = True              ; turn of the full color map
    res@gsnSpreadColorStart  = 2                ; start at color 17
    res@gsnSpreadColorEnd    = 21                ; end at color 200
    
    res@lbLabelBarOn         = True   ; the label bar
    res@lbOrientation        = "horizontal"
    res@lbLabelFontHeightF   = tick_font              ; make labels smaller
    res@lbLabelFont          = font_index              ; make labels smaller
    res@lbLabelStride        = 1
    res@lbBoxLineThicknessF  = mp_thick-1

;------------------------------------------------------------------
;draw the subgraph and change the left and right titles
;-----------------------------------------------------------------------
    res@tiMainString         = fig_title
    res@tiMainFontHeightF    = 0.04
    res@tiMainFontColor      = "black"
    res@gsnCenterString      = ""
    res@gsnStringFont        = font_index
    res@gsnStringFontHeightF = subtitle_font

    plot = gsn_csm_contour_map(wks,std,res)
    print(max(std))
    print(min(std))
end

