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
; Preci data is stored in a file with two dimensions (station160,year67),from 1951 to 2017
; Station information is stored in five columns: number,name,number,lat,lon
; Nino3 data is stored in 10 columns: year(1950-2017), month(12), Nino1+2, Anom, 
;              NINO3, ANOM, NINO4, ANOM, NINO3.4, ANOM
;
  path      = "/home/ys17-19/renql/project/2018homework/statistic/"
  preci_f   = "data/160preci/preci0" + (/"6","7","8"/) + ".txt"
  nino3_f   = "data/ersst5.nino.mth.81-10.ascii"
  station_f = "data/ncc_observation_readme.txt"

  nstation  = 160
  years     = ispan(1951,2017,1)
  nyear     = dimsizes(years)
  nmonth    = 3
  window    = 11  ; the slide window

  preci     = new((/nmonth,nstation,nyear/),float)
  nino3     = new((/nyear ,nmonth/),float)
  lat       = new((/nstation/),float)
  lon       = new((/nstation/),float)

;----- Read JJA precipitation from ascii file which are all numeric--------------------
do nm = 0,nmonth-1,1
  preci(nm,:,:) = tofloat(asciiread(preci_f(nm),(/nstation,nyear/),"integer"))
end do

;-----Read DJF Nino3 SST from ascii file which have character----------------------------
; For string data, each line of the file is read as a string.
  line  = asciiread(nino3_f, -1, "strinng")
  printVarSummary(line)
; Use "str_get_field" to indicate which fields to read in.
; Each field is separated by an arbitrary number of spaces.
do nm = 0,nmonth-1,1
  nino3(:,nm) = tofloat(str_get_field(line((12+n):816:12),5," "))
end do

;-----Read the information of 160 stations from ascii file which have header lines-------
  line := asciiread(station_f, -1, "strinng")
  lat = tofloat(str_get_field(lines(21:),1," "))  ; Attach as lat1d, lon1d
  lon = tofloat(str_get_field(lines(21:),2," "))  ; for plotting later
; if there are number, name, number in file, can not separate each field successfully
; so remove these three columns


  wks = gsn_open_wks("png","station")       ; send graphics to PNG file
  cmap = read_colormap_file("WhViBlGrYeOrRe") ; read color map
  cmap = cmap(::-1,:)                         ; reverse the color map

  res                             = True
  res@gsnMaximize                 = True
  res@gsnFrame                    = False   ; Want to draw markers later.
  res@gsnDraw                     = False

  res@cnLineLabelPlacementMode    = "Constant"
  res@cnLineLabelFontColor        = "Gray15"
  res@cnLineDashSegLenF           = 0.3

  res@cnLevelSelectionMode        = "ManualLevels"
  res@cnMinLevelValF              = 15         	; 15.25
  res@cnMaxLevelValF              = 50		; 49.75
  res@cnLevelSpacingF             = 0.25

  res@cnFillOn                    = True
  res@cnFillPalette               = cmap(2:88,:)
  res@cnFillOpacityF              = 0.75 
  res@cnLinesOn                   = True
  res@cnLineLabelsOn              = True
  res@cnLevelFlags		  = new(139,"string")
  res@cnLevelFlags(:)		  = "NoLine"
  res@cnLevelFlags(0::20)         = "LineAndLabel"

  res@lbOrientation               = "vertical"
  res@lbBoxLinesOn                = False

;---Zoom in on map area of interest
  res@mpMinLatF                   = min(pwv@lat1d)-.5
  res@mpMinLonF                   = min(pwv@lon1d)-.5
  res@mpMaxLatF                   = max(pwv@lat1d)+.5
  res@mpMaxLonF                   = max(pwv@lon1d)+.5

  res@mpFillOn                    = False
  res@mpOutlineDrawOrder          = "PostDraw"
  res@mpFillDrawOrder             = "PreDraw"
  res@mpOutlineBoundarySets       = "GeophysicalAndUSStates"
  res@mpUSStateLineColor          = "Gray10"
  res@mpUSStateLineDashPattern    = 2

  res@pmTickMarkDisplayMode       = "Always"   ; nicer map tickmarks

;---Create the plot; it won't get drawn because gsnFrame was set to False.
  res@tiMainString = "GPS PWV (18Z)"
  plot = gsn_csm_contour_map(wks,pwv,res)

;---Draw markers on the plot in the lat/lon locations.
  mkres               = True
  mkres@gsMarkerIndex = 16     ; Filled circle
  gsn_coordinates(wks,plot,pwv,mkres)
end 
