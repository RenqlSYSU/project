	f1     =  addfile("reg_nh_8_20_rego+uv200.nc","r") 
	f2     =  addfile("reg_nh_8_20_rego+uv850.nc","r")
   rgn_na = (/"SA","EA","NAM","NAF","AU","SAM","SAF" /)
   regu1 = f1->regu8
   regu2 = f2->regu8
   lon = regu1&lon 
   
   rego8 = f1->rego8 
   o8_t  = f1->o8_t 
   regu8 = rego8
   regv8 = rego8
   u8_t = rego8
   v8_t = rego8
     

   do i = 0,3
   if(i.ne.3)then 
   regu8(i,:,:,:) = f1->regu8(i,:,:,:)  
   regv8(i,:,:,:) = f1->regv8(i,:,:,:)  
           
   u8_t(i,:,:,:)  = f1->u8_t(i,:,:,:)   
   v8_t(i,:,:,:)  = f1->v8_t(i,:,:,:) 
   else 
   regu8(i,:,:,:) = f2->regu8(i,:,:,:)  
   regv8(i,:,:,:) = f2->regv8(i,:,:,:)  
           
   u8_t(i,:,:,:)  = f2->u8_t(i,:,:,:)   
   v8_t(i,:,:,:)  = f2->v8_t(i,:,:,:) 
   end if 
   end do 
    
regu8@_FillValue = -9.96921e+36
regv8@_FillValue = -9.96921e+36
rego8@_FillValue = -9.96921e+36


  rego8  = where(abs(o8_t) .le.0.05,rego8 ,rego8@_FillValue) 
  confi_wind = where(u8_t.gt.v8_t, v8_t, u8_t) 
  regu8  = where(confi_wind .le.0.05,regu8 ,rego8@_FillValue) 
  regv8  = where(confi_wind .le.0.05,regv8 ,rego8@_FillValue) 

;;;;==========================================================================
  name = "figs3_reg8_20_nh"  
  wks= gsn_open_wks("eps",name ) 
  base   = new (20, graphic)
  mid_p  = new (20, graphic)
  top_p  = base
  plot   = base
  
  res                      = True             ; plot mods desired
  res@gsnDraw              = False
  res@gsnFrame             = False              
  ;res@gsnMaximize          = True

  res@gsnLeftString        = " "
  res@gsnRightString       = " "

  res@tmXBLabelsOn         = False              ; no bottom labels
  res@tmYRLabelsOn         = False              ; no right labels
  res@tmYLLabelsOn         = False              ; do not draw left labels
  res@tmXTLabelsOn         = False              ; do not draw top labels

;***************************************
  resc = res  
  resc@mpFillOn             = False
  resc@mpGeophysicalLineColor     = "lightsalmon3" ;"wheat4" ;"grey50" ;"transparent" 
  resc@mpGeophysicalLineThicknessF= 1.2
  ;resc@gsnAddCyclic         = True

  resc@mpShapeMode   = "FreeAspect"  
  resc@vpWidthF              = 0.65
  resc@vpHeightF             = 0.4
  resc@mpMinLatF            =  -10;SLAT
  resc@mpMaxLatF            =  60;NLAT
  ;resc@mpCenterLonF         =  0.;180
    
  ;res@mpGridAndLimbOn   = True              ; turn on lat/lon lines
  ;res@mpPerimOn         = False             ; turn off box around plot
  ;res@mpGridLatSpacingF = 60.               ; spacing for lat lines
  ;res@mpGridLonSpacingF = 180.               ; spacing for lon lines

  resc@cnFillOn                 = True
  resc@cnLinesOn                = False
  resc@cnLineLabelsOn           = False         
  resc@cnInfoLabelOn            = False
  resc@lbLabelBarOn             = False ; True
  resc@gsnLeftString        = " "
  resc@gsnRightString       = " "
  resc@tmXBLabelDeltaF      = -0.6 
  resc@tmYLLabelDeltaF      = -0.6
  resc@tmXBLabelFontHeightF = 0.028 
  resc@tmYLLabelFontHeightF = 0.028
  
  ;resc@tmXBMode = "Explicit"
  ;resc@tmXBValues = ispan(-180,120,60)
  ;resc@tmXBLabels = (/"180","120W","60W","0","60E","120E" /)
  resc@tmXBTickSpacingF = 60

  resc@tmXBMinorOn      = True
  resc@tmXBMinorValues  = ispan(-360,360,20)
  
  resc@lbLabelBarOn             = False ; True
  resc@cnFillPalette            = "BlueWhiteOrangeRed"  
  resc@cnLevelSelectionMode     = "ExplicitLevels"
  resc@cnLevels                 = ispan(-12,12,1) 

;resv****************************************** 
  refmag = 0.5
  resv = True
  resv@gsnDraw              = False
  resv@gsnFrame             = False              
  resv@gsnLeftString        = " "
  resv@gsnRightString       = " "
  resv@gsnAddCyclic         = True

  ; resv@vcPositionMode            = "ArrowTail"
  resv@vcGlyphStyle              = "CurlyVector"
  resv@vcLineArrowThicknessF     = 0.5;0.8      
  resv@vcLineArrowColor          =  "grey20" 
  resv@vcMinDistanceF            = 0.021        
  resv@vcMinMagnitudeF           = 0.3          

    ;****²Î¿¼¼ýÍ·****
    resv@vcRefAnnoOn               = False  
    resv@vcRefMagnitudeF           = refmag  
    resv@vcRefLengthF              = 0.1 
    resv@vcRefAnnoBackgroundColor  = "white" 
    resv@vcRefAnnoFontHeightF      = 0.024   
    resv@vcRefAnnoString2          = ""           
;;;========================================================
  
  j = 0
  resc@gsnLeftString = "(a) " + rgn_na(j)
  resc@mpMinLonF            =  40 
  resc@mpMaxLonF            =  160 
  resc@mpCenterLonF         =  (40+160)/2
  do i = 0,16,4
      if (i .eq. 16) then
        resc@tmXBLabelsOn = True            
      else
        resc@tmXBLabelsOn = False           
      end if
    resc@tmYLLabelsOn = True 
    base(i)  = gsn_csm_contour_map_ce(wks,rego8(j,:,:,i/4),resc)       
    resc@gsnLeftString = " "
    plot(i) = gsn_csm_vector(wks,regu8(j,:,:,i/4),regv8(j,:,:,i/4),resv)     
    overlay (base(i),plot(i))         
  end do 
 
  resc@tmYLLabelsOn = False    

  j = 1
  resc@gsnLeftString = "(b) " + rgn_na(j) 
  resc@mpMinLonF            =  60 
  resc@mpMaxLonF            =  179
  resc@mpCenterLonF         =  (60+179)/2
  do i = 1,17,4
      if (i .ge. 16) then
        resc@tmXBLabelsOn = True            
      else
        resc@tmXBLabelsOn = False           
      end if
    base(i)  = gsn_csm_contour_map_ce(wks,rego8(j,:,:,(i-1)/4),resc)       
    resc@gsnLeftString = " "
     plot(i) = gsn_csm_vector(wks,regu8(j,:,:,(i-1)/4),regv8(j,:,:,(i-1)/4),resv)     
    overlay (base(i),plot(i))         
  end do 

  resc@tmYLLabelsOn = False    

  j = 2
  resc@mpMinLonF            =  -120 
  resc@mpMaxLonF            =  0 
  resc@mpCenterLonF         =  (-120-0)/2

  resc@gsnLeftString = "(c) " + rgn_na(j)
  do i = 2,18,4
      if (i .ge. 16) then
        resc@tmXBLabelsOn = True            
      else
        resc@tmXBLabelsOn = False           
      end if
    base(i)  = gsn_csm_contour_map_ce(wks,rego8(j,:,:,(i-2)/4),resc)       
    resc@gsnLeftString = " "
    plot(i) = gsn_csm_vector(wks,regu8(j,:,:,(i-2)/4),regv8(j,:,:,(i-2)/4),resv)     
    overlay (base(i),plot(i))         
  end do 

     resc@tmYRLabelsOn = True    

  j = 3
  resc@mpMinLonF            =  -121 
  resc@mpMaxLonF            =  40 
  resc@mpCenterLonF         = (-121+40)/2
  resc@gsnLeftString = "(d) " + rgn_na(j) ;+"_hgt" ;"_OLR"
  do i = 3,19,4
      if (i .ge. 16) then
        resc@tmXBLabelsOn              = True 
        resv@vcRefAnnoOn               = True           ; 矢量 
        resv@vcRefAnnoParallelPosF     = 1
        resv@vcRefAnnoOrthogonalPosF   = 0.2
        resv@vcRefAnnoPerimOn          = False          ;参考矢量的边框
        resv@vcRefAnnoString1          = refmag+" m s~S~-1 "       ;reference vector
      else
        resc@tmXBLabelsOn = False           
      end if
    base(i)  = gsn_csm_contour_map_ce(wks,rego8(j,:,:,(i-3)/4),resc)       
    resc@gsnLeftString = " "
    plot(i) = gsn_csm_vector(wks,regu8(j,:,:,(i-3)/4),regv8(j,:,:,(i-3)/4),resv)     
    overlay (base(i),plot(i))         
  end do 
     nlon = dimsizes(lon)
     
  ;;equatorial   
     resp                  = True                      ; polyline mods desired
     resp@gsLineColor      = "grey50"                     ; color of lines
     resp@gsLineThicknessF = 1.0                       ; thickness of lines
     resp@gsLineDashPattern= 12

   dum1 = new(20,graphic)
   do i = 0,19
    dum1(i) = gsn_add_polyline(wks,base(i),(/0,180,360/),(/0.,0.,0./),resp)
   end do 
;;;;************************************************************
  resP                     = True               ; panel options
  resP@gsnPanelCenter      = False
  ;resP@gsnPanelDebug       = True
  resP@gsnMaximize         = True               ; maximize image
;---Scale plots based on one of the smaller plots, and not the first one.
  resP@gsnPanelScalePlotIndex     = 4
  resP@gsnPanelTop                = 0.90
  resP@gsnPanelBottom             = 0.10
  
  resP@gsnPanelLabelBar         = True
  ;resP@lbBoxEndCapStyle         = "TriangleBothEnds"
  resP@pmLabelBarWidthF         = 0.65 
  resP@pmLabelBarHeightF        = 0.045
 resP@pmLabelBarOrthogonalPosF = 0.04  
  resP@pmLabelBarParallelPosF   = 0.02 
  resP@lbLabelFontHeightF       = 0.01 
  resP@lbTitleOn                =  True                ; turn on title
  resP@lbTitleString            =  "W m~S~-2"  
  resP@lbTitlePosition          = "Right"              ; title position
  resP@lbTitleFontHeightF       = .012                 ; make title smaller
  resP@lbTitleOffsetF           = 0.012
  resP@lbTitleDirection         = "Across"             ; title direction

  resP@gsnPanelFigureStrings = (/"-8","","","","-6","","","","-4","","","",\
                               "-2","","","","0","","",""/)
  resP@gsnPanelFigureStringsFontHeightF = 0.01
  resP@amJust              = "TopLeft"
 
  a1 = 0.056727
  a2 = 0.286727
  a3 = 0.516727
  a4 = 0.746727
  resP@gsnPanelXF = (/a1,a2,a3,a4,a1,a2,a3,a4,a1,a2,a3,a4,a1,a2,a3,a4,a1,a2,a3,a4/)

  gsn_panel(wks,base,(/5,4/),resP)  

