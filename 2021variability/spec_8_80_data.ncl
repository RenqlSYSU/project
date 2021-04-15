undef("f_specx_yr")
function f_specx_yr(prc1)
begin
  dims = dimsizes(prc1)

  d   = 0
  sm  = 1         ; periodogram
  pct = 0.10
  ;   dof   = specx_anal(prcwgt1(0,5:),d,sm,pct)
  prcspc = new((/dims(0),dims(1)/2/),typeof(prc1))
  spcavg = new(dims(1)/2,typeof(prc1))
  spcavg = 0.0
  r1zsum = 0.0
  splt2 = new((/dims(0),4,dims(1)/2/),typeof(prc1))
  splt1 = new((/4,dims(1)/2/),typeof(prc1))
  splt1  = 0.0
  
  do yr = 0,dims(0)-1
   dof    = specx_anal(prc1(yr,:),d,sm,pct)  
   r1     = dof@xlag1                        ; extract segment lag-1
   r1zsum = r1zsum  + 0.5*(log((1+r1)/(1-r1)))
   spcavg = spcavg + dof@spcx*dof@frq
   spcavg1 = dof@spcx*dof@frq
  end do

  r1z  = r1zsum/dims(0)                           ; average r1z
  r1   = (exp(2*r1z)-1)/(exp(2*r1z)+1)            ; transform back
   
  spcavg = spcavg/dims(0)

  df = 2*dims(0);123./11.;   ;; deg of freedom dof 
  df@spcx = spcavg   ;<ARRAY of 56 elements>
  df@frq  = dof@frq  ;<ARRAY of 56 elements>
  df@xlag1= r1      
 
  p = 1/dof@frq
  ;printVarSummary(p) 
  ;splt = specx_ci(df, 0.05, 0.99)
  ;*******************************************************************
   twopi = 2.*3.14159
   pLow  = 0.01
   pHigh = 0.95 
   N     = dimsizes (spcavg)      ; number specx estimates
    ;; df    = (/ dof /)
  
   r     = dof@xlag1
   r2    = 2.*r
   rsq   = r*r
  
   temp  = r2*cos(twopi*dof@frq)    ; vector
                                     
   mkov  = 1./(1. + rsq - temp)      ; Markov Model
  
   sum1  = sum (mkov)                ; sum Markov elements
   sum2  = sum (spcavg)           ; sum spectral elements      sdof@spcx
   scale = sum2/sum1                 ; scaling factor
  
   xLow  = chiinv (pLow,  df)/df     ; lower confidence
   xHigh = chiinv (pHigh, df)/df     ; upper confidence
  
   splt      = new ( (/5,N/), typeof(dof)) ; 4 spec curves
   splt(0,:) = spcavg                ; input spectrum
   splt(1,:) = mkov*scale            ; Markov         
   splt(2,:) = splt(1,:)*xHigh       ; low  ci for Markov
   splt(3,:) = splt(1,:)*xLow        ; high ci for Markov
   splt(4,:dimsizes(p)-1) = 1/dof@frq
  return(splt)
end 

begin
   dir   = "/data3/home/liufeiclim2014/"
   f1     = addfile(dir+"oyy/data/pre_8-80day.sam.gt300.5-9.nc","r")        
   f2     = addfile(dir+"oyy/data/pre_8-80day.sam.gt300.11-3.nc","r")        
   f5     = addfile(dir+"oyy/nclcode/20201030/monsoon.domain.1979-2017.nc","r")
   area    = f5->area
   
   lat_st = (/0   ,0   ,0    ,0    ,-35 ,-40  ,-40   /)
   lat_ed = (/40  ,60  ,30   ,30   ,0   ,0    ,0     /)
   lon_st = (/60  ,100 ,250  ,-20  ,100 ,275  ,5     /)
   lon_ed = (/100 ,140 ,310  ,45   ,155 ,325  ,50    /)
   rgn_na = (/"SA","EA","NAM","NAF","AU","SAM","SAF" /)

   lat1 = (/25  ,29  ,5.5    ,5.5  ,-16 ,-14  ,-13   /)
   lon1 = (/90  ,117 ,293  ,9    ,134 ,309  ,28     /)
   lat2 = (/19  ,29  ,19    ,5.5  ,-16 ,-14  ,-13   /)
   lon2 = (/79.5  ,117 ,264  ,9    ,134 ,309  ,28     /)

   ;9.96921e+36
    ;printVarSummary(area)
     i = 0;0,6

     if(i.le.3)then 
      day = 153
     else 
      day = 151
     end if 
;;;area_average-----------------------------------------------------
     if(i.le.2) then  
       pr8  = f1->pre8_80(:,{lat_st(i):lat_ed(i)},{lon_st(i):lon_ed(i)})
       area1 = area({lat_st(i):lat_ed(i)},{lon_st(i):lon_ed(i)})
     else if(i.eq.3)then
       pr8_1  = f1->pre8_80(: ,{lat_st(i):lat_ed(i)},:)
       p8_80  = lonFlip(pr8_1)      
       pr8  = p8_80(:,:,{lon_st(i):lon_ed(i)})
       
       area := area({lat_st(i):lat_ed(i)},:)
       area  := lonFlip(area)
       area1 = area(:,{lon_st(i):lon_ed(i)})
     else
       pr8  = f2->pre8_80(: ,{lat_st(i):lat_ed(i)},{lon_st(i):lon_ed(i)})
       area1 = area({lat_st(i):lat_ed(i)},{lon_st(i):lon_ed(i)})
     end if 
     end if
     area_con = conform(pr8, area1, (/1,2/))   
     pr8  = where(area_con.gt.1,pr8 ,pr8@_FillValue)
     pr8_wgt = dim_avg_n(pr8,(/1,2/))
;;;station----------------------------------------------------------
   ; if(i.le.3.and.i.ge.0)then
   ;   pr8_wgt  = f1->pre8_80(:,{lat1(i)},{lon1(i)})
   ; else
   ;   pr8_wgt  = f2->pre8_80(:,{lat1(i)},{lon1(i)})
   ; end if
;;;-----------------------------------------------------------------

;printVarSummary(pr8)
;printVarSummary(area1)
    ;print(pr8_wgt(::100)) ;[time | 39] x [lat | 80] x [lon | 80]
 
   pr8_wgt      = linmsg_n(pr8_wgt,-1,0)    ;5967
   pr8_wgt@_FillValue = -9.96921e+36
  if(any(ismissing(pr8_wgt))) then
    print("Your data contains some missing values. Beware.")    
  end if    
;exit
   dim_p = dimsizes(pr8_wgt)
   yr_pr = reshape(pr8_wgt,(/dim_p/day,day/))
    
   spec  = f_specx_yr(yr_pr)
printVarSummary(spec)
  path_out=dir+"oyy/nclcode/20210310/spec8_80_"+rgn_na(i)+"_sam.gt300.nc"
  system("rm -f "+path_out)
  ncdf=addfile(path_out,"c")
  ncdf->spec  = spec 

end 
