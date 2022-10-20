
;   Example script to produce standard plots for a WRF bwave run

load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/wrf/WRFUserARW.ncl"
;load "./WRFUserARW.ncl"

begin
;
; The WRF ARW input file.  
; This needs to have a ".nc" appended, so just do it.

  a = addfile("../IDEAL/wrfout_bwave.nc","r")

;  We generate plots, but what kind do we prefer?

   type = "x11"
;  type = "pdf"
;  type = "ps"
;  type = "ncgm"

  wks = gsn_open_wks(type,"plt_BWave")


; Set some basic resources
    res = True
    res@MainTitle     = "WRF BWAVE"     
    res@InitTime      = False
    res@TimePos       = "Left"
    res@Footer        = False
    res@vpWidthF      = .3           ; Overwite basic plot size
    res@vpHeightF     = .6
    res@lbTitleOn     = False
    res@pmLabelBarSide = "Right"
    res@lbOrientation = "Vertical"

    pltres = True

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;  What times and how many time steps are in the data set?

  times  = wrf_user_list_times(a)  ; get times in the file
  ntimes = dimsizes(times)         ; number of times in the file

; The specific height levels that we want the data interpolated to.

  height_levels = (/ 250., 2000./)   ; heigth levels to plot
  nlevels         = dimsizes(height_levels)     ; number of height levels


; This is the big loop over all of the time periods to process.

;  do it = 0,ntimes-1,2
  do it = 18,ntimes-1

    time = it
    res@TimeLabel                  = times(it)

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; First get the variables we will need        

  p   = wrf_user_getvar(a, "pressure",time)    ; pressure
  th  = wrf_user_getvar(a,"th",time)    ; get temperature (C)
  u   = wrf_user_getvar(a,"ua",time)    ; ua is u averaged to mass points
  v   = wrf_user_getvar(a,"va",time)    ; va is v averaged to mass points
  w   = wrf_user_getvar(a,"wa",time)    ; vertical velocity
  z   = wrf_user_getvar(a, "z",time)    ; grid point height
  ter = wrf_user_getvar(a,"HGT",time)   ; need terrain height sometimes

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

    do level = 0,nlevels-1

      height = height_levels(level) 

 ; Pressure
      p_plane  = wrf_user_intrp3d( p,z,"h",height,0.,False)
      opts_p   = res
      opts_p@ContourParameters     = (/ 2. /)
      opts_p@cnInfoLabelFontHeightF = 0.01
      contour_p     = wrf_contour(a,wks,p_plane,opts_p)

 ; Theta
      th_plane = wrf_user_intrp3d(th,z,"h",height,0.,False)
      opts_th  = res
      opts_th@cnFillOn             = True
      opts_th@gsnSpreadColorEnd    = -10   
      opts_th@pmLabelBarOrthogonalPosF = -0.1
      contour_th     = wrf_contour(a,wks,th_plane,opts_th)

 ; Vertical Velocity
      w_plane  = wrf_user_intrp3d( w,z,"h",height,0.,False)
      w_plane = 100.*w_plane
      w_plane@units  = "m/3"
      opts_w   = res
      opts_w@ContourParameters     = (/ 1. /)
      opts_w@cnFillOn              = True
      opts_w@gsnSpreadColorEnd     = -3
      opts_w@pmLabelBarOrthogonalPosF = 0.0
      contour_w = wrf_contour(a,wks, w_plane,opts_w)

 ; Wind Vectors
      u_plane  = wrf_user_intrp3d( u,z,"h",height,0.,False)
      v_plane  = wrf_user_intrp3d( v,z,"h",height,0.,False)
      u_plane@description = "Wind"
      v_plane@description = "Wind"
      opts_vct = res
      opts_vct@NumVectors          = 15 
      opts_vct@vcGlyphStyle        = "LineArrow"
      vector         =  wrf_vector(a,wks,u_plane, v_plane,opts_vct)

      plot = wrf_overlays(a,wks,(/contour_p, contour_th, vector/),pltres)
      plot = wrf_overlays(a,wks,(/contour_w, vector/),pltres)

    end do

  ;       ************************************************************

  end do  ;       end of the time loop

end
