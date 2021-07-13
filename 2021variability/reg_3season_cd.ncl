begin 
;;---------------------------------------------------
;;band filter 
  ihp     = 2             ; band pass
  sigma   = 1.0           ; Lanczos sigma 
  
  a = 20                           
  nWgt    = 2*a+1                       ; loose 100 each end                            
  fca     = 1./30                       ; start freq
  fcb     = 1./8                        ; last  freq  -999.;           
  wgt     = filwgts_lanczos (nWgt, ihp, fca, fcb, sigma )   
;;read file-----------------------------------------------------  
  f1  = addfile("u_0401_0930.nc","r")    
  f2  = addfile("v_0401_0930.nc","r")   
  f4  = addfile("olr_0401_0930.nc","r")
  f5  = addfile("WAF200_reg_stream.nc","r")
  u      = f1 ->u  (0:37,:,:,{0:70},{0:180})   
  v      = f2 ->v  (0:37,:,:,{0:70},{0:180})   
  fx    = f5 ->fx(:,{0:70},{0:180})
  fy    = f5 ->fy(:,{0:70},{0:180})
  olr = short2flt(f4 ->olr(0:37,:,{0:70},{0:180}))    ;4.30-8.31(10:133)

  ;;;PC------------------------------------------------------------------
   f6  = addfile("eof5-8.nc","r")
   ts_mjja = doubletofloat(f6 ->eof_ts) ;[evn | 3] x [day | 4674] 
   ts_mjja := reshape(ts_mjja,(/3,38,123/))
   sday = (/ 0,45, 76,0/)  ;0:44 early summer; 45:75 meiyu season; 76:122 late summer 
   eday = (/44,75,122,122/)
   iday = 0 ;2 
     ii = iday*8 ;16 
   pc1  = reshape(ts_mjja(0,:,sday(iday):eday(iday)),38*(eday(iday)-sday(iday)+1))
   pc2  = reshape(ts_mjja(1,:,sday(iday):eday(iday)),38*(eday(iday)-sday(iday)+1))
                                        ;45 31 47
  name = "4x2yr_reg_mj" 
;;;---------------------------------------------------------- 
  ;bandpass filter
  u1   = wgt_runave_n_Wrap ( u  , wgt, 0, 1 ) 
  v1   = wgt_runave_n_Wrap ( v  , wgt, 0, 1 ) 
  olr1 = wgt_runave_n_Wrap ( olr, wgt, 0, 1 ) 

  ;;regression
   olr2    = olr1(        lat|:,lon|:,year|:,time|:)             
   u2      =   u1(level|:,lat|:,lon|:,year|:,time|:)
   v2      =   v1(level|:,lat|:,lon|:,year|:,time|:)
  
  dio = dimsizes(olr2)
  diu = dimsizes(u2)
  rego = new((/dio(0),dio(1),4/),typeof(olr2))                                             

  regu = new((/diu(0),diu(1),diu(2),4/),typeof(u2))
  regv = new((/diu(0),diu(1),diu(2),4/),typeof(v2)) 
  rego1 = rego
  regu1 = regu
  regv1 = regv
  
  o_tval    = new((/dio(0)*dio(1),4/),typeof(olr2))
  o_nptxy   = new((/dio(0)*dio(1),4/),typeof(olr2))
  o1_tval   = new((/dio(0)*dio(1),4/),typeof(olr2))
  o1_nptxy  = new((/dio(0)*dio(1),4/),typeof(olr2))

  u_tval    = new((/diu(0)*diu(1)*diu(2),4/),typeof(u2))
  u_nptxy   = new((/diu(0)*diu(1)*diu(2),4/),typeof(u2))
  u1_tval   = new((/diu(0)*diu(1)*diu(2),4/),typeof(u2))
  u1_nptxy  = new((/diu(0)*diu(1)*diu(2),4/),typeof(u2))
  v_tval    = new((/diu(0)*diu(1)*diu(2),4/),typeof(u2))
  v_nptxy   = new((/diu(0)*diu(1)*diu(2),4/),typeof(u2))
  v1_tval   = new((/diu(0)*diu(1)*diu(2),4/),typeof(u2))
  v1_nptxy  = new((/diu(0)*diu(1)*diu(2),4/),typeof(u2))
  ;[year | 38] x [lat | 73] x [lon | 144] x [time | 183]
;;;regress================================================================================================
    do i = 0,6,2   ;30-74 75-105 106-152
   olr3    := reshape(olr2(:,:,:,sday(iday)+30-i:eday(iday)+30-i),(/dio(0),dio(1),dio(2)*(eday(iday)-sday(iday)+1)/))
   u3      := reshape(u2(:,:,:,:,sday(iday)+30-i:eday(iday)+30-i),(/diu(0),diu(1),diu(2),diu(3)*(eday(iday)-sday(iday)+1)/))
   v3      := reshape(v2(:,:,:,:,sday(iday)+30-i:eday(iday)+30-i),(/diu(0),diu(1),diu(2),diu(3)*(eday(iday)-sday(iday)+1)/))
      rego(:,:,i/2)    = regCoef(dim_standardize_n_Wrap(pc1,1,0),olr3)
      regu(:,:,:,i/2)  = regCoef(dim_standardize_n_Wrap(pc1,1,0),u3)         
      regv(:,:,:,i/2)  = regCoef(dim_standardize_n_Wrap(pc1,1,0),v3)  
      rego1(:,:,i/2)   = regCoef(dim_standardize_n_Wrap(pc2,1,0),olr3)
      regu1(:,:,:,i/2) = regCoef(dim_standardize_n_Wrap(pc2,1,0),u3)         
      regv1(:,:,:,i/2) = regCoef(dim_standardize_n_Wrap(pc2,1,0),v3)  
      
       o_tval  (:,i/2)= rego@tval
       o_nptxy (:,i/2)= rego@nptxy 
       o1_tval (:,i/2)= rego1@tval 
       o1_nptxy(:,i/2)= rego1@nptxy 
       u_tval  (:,i/2)= regu@tval
       u_nptxy (:,i/2)= regu@nptxy 
       u1_tval (:,i/2)= regu1@tval
       u1_nptxy(:,i/2)= regu1@nptxy
       v_tval  (:,i/2)= regv@tval
       v_nptxy (:,i/2)= regv@nptxy 
       v1_tval (:,i/2)= regv1@tval
       v1_nptxy(:,i/2)= regv1@nptxy
    end do 
      ;t-test 
       o_t     =  student_t(o_tval ,38*(eday(iday)-sday(iday)+1)/8)   ;( o_nptxy-2)/8 
       o1_t    =  student_t(o1_tval,38*(eday(iday)-sday(iday)+1)/8)   ;(o1_nptxy-2)/8
       u_t     =  student_t(u_tval ,38*(eday(iday)-sday(iday)+1)/8)   ;( u_nptxy-2)/8
       u1_t    =  student_t(u1_tval,38*(eday(iday)-sday(iday)+1)/8)   ;(u1_nptxy-2)/8
       v_t     =  student_t(v_tval ,38*(eday(iday)-sday(iday)+1)/8)   ;( v_nptxy-2)/8
       v1_t    =  student_t(v1_tval,38*(eday(iday)-sday(iday)+1)/8)   ;(v1_nptxy-2)/8
       
       o_t  := onedtond(ndtooned(o_t) ,(/dio(0),dio(1),4/))         
       o1_t := onedtond(ndtooned(o1_t),(/dio(0),dio(1),4/))         
       u_t  := onedtond(ndtooned(u_t) ,(/diu(0),diu(1),diu(2),4/))         
       u1_t := onedtond(ndtooned(u1_t),(/diu(0),diu(1),diu(2),4/))         
       v_t  := onedtond(ndtooned(v_t) ,(/diu(0),diu(1),diu(2),4/))         
       v1_t := onedtond(ndtooned(v1_t),(/diu(0),diu(1),diu(2),4/))       

   copy_VarCoords_1(olr2(:,:,0,0:3)  ,rego)     
   copy_VarCoords_1(u2  (:,:,:,0,0:3),regu)   
   copy_VarCoords_1(v2  (:,:,:,0,0:3),regv)   

   copy_VarCoords_1(olr2(:,:,0,0:3)  ,rego1)   
   copy_VarCoords_1(u2  (:,:,:,0,0:3),regu1)   
   copy_VarCoords_1(v2  (:,:,:,0,0:3),regv1)   
   
   copy_VarCoords_1(olr2(:,:,0,0:3)  ,o_t )   
   copy_VarCoords_1(olr2(:,:,0,0:3)  ,o1_t)   
   copy_VarCoords_1(u2  (:,:,:,0,0:3),u_t )   
   copy_VarCoords_1(u2  (:,:,:,0,0:3),u1_t)   
   copy_VarCoords_1(u2  (:,:,:,0,0:3),v_t )   
   copy_VarCoords_1(u2  (:,:,:,0,0:3),v1_t)   

;;significance 95%  
  wind_confi = where(u_t.lt.v_t, u_t, v_t)
  regu = where(wind_confi.le.0.05,regu,0)
  regv = where(wind_confi.le.0.05,regv,0)
  wind_confi1 = where(u1_t.lt.v1_t, u1_t, v1_t)
  regu1 = where(wind_confi1.le.0.05,regu1,0)
  regv1 = where(wind_confi1.le.0.05,regv1,0)
  
  rego  = where(o_t .le.0.05,rego ,0.)  
  rego1 = where(o1_t.le.0.05,rego1,0.)  
  
  f_c1 = new((/diu(1),diu(2),8/),typeof(regu))
  f_c1(:,:,0:3) =  wind_confi(9,:,:,:)
  f_c1(:,:,4:7) = wind_confi1(9,:,:,:)
  f_c1!2 = "time"
  f_c1!0 = "lat"
  f_c1!1 = "lon"
  f_c1 := wind_c1(time|:,lat|:,lon|:)
  fx(iday*8:iday*8+7,:,:) = where(f_c1.le.0.05,fx(iday*8:iday*8+7,:,:),0)
  fy(iday*8:iday*8+7,:,:) = where(f_c1.le.0.05,fy(iday*8:iday*8+7,:,:),0)
;;;;plot==========================================================================
  wks= gsn_open_wks("eps",name ) 
refmag = 1.
  base   = new (8, graphic)
  top_p  = new (8, graphic)
  mid_p  = new (8, graphic)

  res                      = True             ; plot mods desired
  res@gsnDraw              = False
  res@gsnFrame             = False              

  res@vpWidthF              = 0.65
  res@vpHeightF             = 0.4
  res@gsnLeftString        = " "
  res@gsnRightString       = " "

  res@tmXBLabelsOn         = False              ; no bottom labels
  res@tmYRLabelsOn         = False              ; no right labels
  res@tmYLLabelsOn         = False              ; do not draw left labels
  res@tmXTLabelsOn         = False              ; do not draw top labels
 ; olr **************************************
  resc = res  
  resc@mpFillOn             = True
  resc@mpGeophysicalLineColor     = "grey70"
  resc@mpGeophysicalLineThicknessF= 1.6

  resc@gsnAddCyclic         = False
  resc@mpShapeMode          = "FreeAspect"
  resc@mpMinLatF            = 0
  resc@mpMaxLatF            = 70
  resc@mpMinLonF            = 0
  resc@mpMaxLonF            = 170

  resc@cnFillOn                 = True
  resc@cnLinesOn                = False
  resc@cnLineLabelsOn           = False         
  resc@cnInfoLabelOn            = False
  resc@lbLabelBarOn             = False ; True
  resc@gsnLeftString        = " "
  resc@gsnRightString       = " "
  resc@tmXBLabelDeltaF      = -0.6  ;使其更靠近坐标轴
  resc@tmYLLabelDeltaF      = -0.6
  resc@tmXBLabelFontHeightF = 0.028 ;字体大小
  resc@tmYLLabelFontHeightF = 0.028

; u v ****************************************** 
  resv = True
  resv@gsnDraw              = False
  resv@gsnFrame             = False              
  resv@gsnLeftString        = " "
  resv@gsnRightString       = " "

  resv@vcPositionMode            = "ArrowTail"  ;格点位置对应箭头的尾部
  resv@vcGlyphStyle              = "CurlyVector";"LineArrow"  ;"Fillarrow"
  resv@vcLineArrowThicknessF     = 1.          ; 箭头边界线的粗细
  resv@vcLineArrowColor          = "Foreground"  ;grey30"mediumvioletred";"maroon3";"orangered";"darkorchid1";"deeppink3";"maroon1" 
  resv@vcMinDistanceF            = 0.017         ; 箭头之间的最小距离
  resv@vcMinMagnitudeF           = 0.3          ; 小于该数值则不绘制
    ;****参考箭头****
    resv@vcRefAnnoOn               = True  
    resv@vcRefMagnitudeF           = refmag  ;单位长度箭头所表示的实际物理量的大小
    resv@vcRefLengthF              = 0.04    ;单位长度箭头在图形中的长度（单位坐标系中）
    resv@vcRefAnnoBackgroundColor  = "white" ;背景颜色     
    resv@vcRefAnnoPerimOn          = True   ;关闭边框                                           
    resv@vcRefAnnoFontHeightF      = 0.015   ;参考箭头标签字体大小          
    resv@vcRefAnnoString1On        = False   ;不设定参考箭头上方的字符        
    resv@vcRefAnnoString2On        = True    ;设定参考箭头下方的字符
    resv@vcRefAnnoString2          = refmag           
    resv@vcRefAnnoSide             = "Top"   ; 参考箭头放置图形上方
    resv@vcRefAnnoOrthogonalPosF   = -0.17   ; 调整其垂直位置
    resv@vcRefAnnoParallelPosF     = 0.99   
    ;resv@vcGlyphOpacityF           = 0.8     ;透明度
resv1 = resv
    resv1@vcRefAnnoOrthogonalPosF   = -0.34 ; 调整其垂直位置
;;waf =========================
resv2 = resv 
  refmag1 = 10;20. ;5.
  resv2@vcRefMagnitudeF           = refmag1  ;单位长度箭头所表示的实际物理量的大小
  resv2@vcRefAnnoString2          = refmag1           
  resv2@vcGlyphStyle              = "LineArrow";"CurlyVector"  ;"Fillarrow"
  resv2@vcLineArrowThicknessF     = 2.          ; 箭头边界线的粗细
  resv2@vcLineArrowColor          = "mediumvioletred";"maroon3";"orangered";"darkorchid1";"deeppink3";"maroon1" 
  resv2@vcRefAnnoOrthogonalPosF   = -0.34   ; 调整其垂直位置
  resv2@vcRefAnnoParallelPosF     = 0.89   
  resv2@vcMinDistanceF            = 0.028;32         ; 箭头之间的最小距离
  resv2@vcRefLengthF              = 0.02    ;单位长度箭头在图形中的长度（单位坐标系中）
  resv2@vcMinMagnitudeF           = 2.          ; 小于该数值则不绘制

;panel 1***********************************************************************  
  resc@cnFillPalette            = "BlueWhiteOrangeRed" ;cmap1 ;NCV_blue_red ;MPL_bwr；MPL_BrBG
  resc@cnLevelSelectionMode     = "ExplicitLevels"
  resc@cnLevels                 = ispan(-10, 10, 1);(-24, 24, 2)

  resc@gsnLeftString = "(a) UM" 
  resc@lbLabelBarOn             = False ; True
  
  j = 3
  i = 0 
      if (i%2.eq. 0) then
        resc@tmYLLabelsOn = True
        resc@tmYLOn       = True
      else
        resc@tmYLLabelsOn = False
      end if

    base(i)  = gsn_csm_contour_map_ce(wks,rego(:,:,j),resc)       
    top_p(i) = gsn_csm_vector(wks,regu({200},:,:,j),regv({200},:,:,j),resv1)     
    mid_p(i) = gsn_csm_vector(wks,fx(j+ii,:,:),fy(j+ii,:,:),resv2)
    overlay (base(i),mid_p(i))     
    overlay (base(i),top_p(i))  
    resv@vcRefAnnoOn               = False  
    resv1@vcRefAnnoOn               = False  
    resv2@vcRefAnnoOn               = False  

  resc@gsnLeftString = " " 

  do i = 2,4,2
    j = j-1
      if (i%2.eq. 0) then
        resc@tmYLLabelsOn = True
        resc@tmYLOn       = True
      else
        resc@tmYLLabelsOn = False
      end if
    base(i)  = gsn_csm_contour_map_ce(wks,rego(:,:,j),resc)       
    top_p(i) = gsn_csm_vector(wks,regu({200},:,:,j),regv({200},:,:,j),resv1)     
    mid_p(i) = gsn_csm_vector(wks,fx(j+ii,:,:),fy(j+ii,:,:),resv2)
    overlay (base(i),mid_p(i))     
    overlay (base(i),top_p(i))  
  end do 

 j = 0
 i = 6
      if (i .ge. 6) then
        resc@tmXBLabelsOn = True            
        resc@tmXBOn       = True 
      else
        resc@tmXBLabelsOn = False           
      end if

      if (i%2 .eq. 0) then
        resc@tmYLLabelsOn = True
        resc@tmYLOn       = True
      else
        resc@tmYLLabelsOn = False
      end if

    base(i)  = gsn_csm_contour_map_ce(wks,rego(:,:,j),resc)       
    top_p(i) = gsn_csm_vector(wks,regu({200},:,:,j),regv({200},:,:,j),resv1)     
    mid_p(i) = gsn_csm_vector(wks,fx(j+ii,:,:),fy(j+ii,:,:),resv2)
    overlay (base(i),mid_p(i))     
    overlay (base(i),top_p(i))     
;;panel 2************************************************************* 
  resc@tmXBLabelsOn         = False              ; no bottom labels
  resc@tmYRLabelsOn         = False              ; no right labels
  resc@tmYLLabelsOn         = False              ; do not draw left labels
  resc@tmXTLabelsOn         = False              ; do not draw top labels

    resc@gsnLeftString = "(b) DM" 
    resv@vcRefAnnoOn               = True  
    resv1@vcRefAnnoOn                = True  
    resv2@vcRefAnnoOn                = True 
    refmag = 30
    resv2@vcRefMagnitudeF           = refmag  ;单位长度箭头所表示的实际物理量的大小
    resv2@vcRefAnnoString2          = refmag 
    resv2@vcMinDistanceF            = 0.032         ; 箭头之间的最小距离
    resv2@vcMinMagnitudeF           = 5.            ; 小于该数值则不绘制

i = 1
    j = 3 
    base(i)  = gsn_csm_contour_map_ce(wks,rego1(:,:,j),resc)       
    top_p(i) = gsn_csm_vector(wks,regu1({200},:,:,j),regv1({200},:,:,j),resv1)     
    mid_p(i) = gsn_csm_vector(wks,fx(7+ii,:,:),fy(7+ii,:,:),resv2)
    overlay (base(i),mid_p(i))     
    overlay (base(i),top_p(i))  

    resc@gsnLeftString = "" 
    resv@vcRefAnnoOn                = False  
    resv1@vcRefAnnoOn               = False  
    resv2@vcRefAnnoOn               = False  

  j = 3
  do i = 3,5,2
    j = j -1 
    base(i)  = gsn_csm_contour_map_ce(wks,rego1(:,:,j),resc)       
    eop_p(i) = gsn_csm_vector(wks,regu1({200},:,:,j),regv1({200},:,:,j),resv1)     
    mid_p(i) = gsn_csm_vector(wks,fx(4+j+ii,:,:),fy(4+j+ii,:,:),resv2)
    overlay (base(i),mid_p(i))     
    overlay (base(i),top_p(i))  
  resc@gsnLeftString = " " 
  end do    
  resc@gsnLeftString = "" 
    resv@vcRefAnnoOn                = False  
    resv1@vcRefAnnoOn               = False  
 j = 0
 i = 7
      if (i .ge. 6) then
        resc@tmXBLabelsOn = True            
        resc@tmXBOn       = True 
      else
        resc@tmXBLabelsOn = False           
      end if
    base(i)  = gsn_csm_contour_map_ce(wks,rego1(:,:,j),resc)       
    top_p(i) = gsn_csm_vector(wks,regu1({200},:,:,j),regv1({200},:,:,j),resv1)     
    mid_p(i) = gsn_csm_vector(wks,fx(4+ii,:,:),fy(4+ii,:,:),resv2)
    overlay (base(i),mid_p(i))     
    overlay (base(i),top_p(i)) 
;;;;************************************************************
  resP                     = True               ; panel options
  resP@gsnPanelCenter      = False
  resP@gsnPanelDebug       = True
  resP@gsnMaximize         = True               ; maximize image
  resP@gsnLeftString = " "
  resP@gsnRightString= " "
;---Scale plots based on one of the smaller plots, and not the first one.
  resP@gsnPanelScalePlotIndex     = 4
  resP@gsnPanelTop                = 0.90
  resP@gsnPanelBottom             = 0.10
  
  resP@gsnPanelLabelBar         = True
  ;resP@lbBoxEndCapStyle         = "TriangleBothEnds"
  resP@pmLabelBarWidthF         = 0.45                 ; 色标宽度default is shorter
  resP@pmLabelBarHeightF        = 0.05                 ; 高度default is taller
 ; resP@pmLabelBarOrthogonalPosF = 0.11                ;调整色标垂直方向距离
 resP@pmLabelBarParallelPosF   = 0.03                  ;调整色标水平方向距离
  resP@lbLabelFontHeightF       = 0.012                ; make labels smaller ( default=0.02 )
  resP@lbTitleOn                =  True                ; turn on title
  resP@lbTitleString            = "W m~S~-2"           ; title string    
  resP@lbTitlePosition          = "Right"              ; title position
  resP@lbTitleFontHeightF       = .012                 ; make title smaller
  resP@lbTitleOffsetF           = 0.012
  resP@lbTitleDirection         = "Across"             ; title direction

  resP@gsnPanelFigureStrings = (/"-6","-6","-4","-4","-2","-2","0","0"/)
  resP@gsnPanelFigureStringsFontHeightF = 0.01
  resP@amJust              = "TopLeft"

  a1 = 0.891923
  a2 = 0.691923
  a3 = 0.491923
  a4 = 0.291923
  resP@gsnPanelYF = (/a1,a1,a2,a2,a3,a3,a4,a4/)
  a1 = 0.212247
  a2 = 0.525082
  resP@gsnPanelXF = (/a1,a2,a1,a2,a1,a2,a1,a2/)

  gsn_panel(wks,base,(/4,2/),resP)  

end 
