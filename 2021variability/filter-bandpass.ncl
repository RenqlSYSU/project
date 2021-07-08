begin
   f1     = addfile("precip.day.1979-2018.nc","r")   ;40x360x720                            ; open netcdf file
   precip = f7->precip
;;=======================================================================
;;band filter 
  ihp     = 2             ; band pass

 ;ihp     = 1             ; high pass
 ;ihp     = 0             ; low_pass  
  sigma   = 1.0           ; Lanczos sigma 
  
  a = 25                           
  nWgt    = 2*a+1                       ; loose 25 each end                            
  fca     = 1./80                       ; start freq
  fcb     = 1./8                        ; last  freq  -999.;     

  wgt8_80      = filwgts_lanczos (nWgt, ihp, fca, fcb, sigma )   ;123 X 24 先小后大

  pre8_80   = wgt_runave_n_Wrap ( precip, wgt8_80  , 0, 0 ) 
  ;========output===
setfileoption("nc","format","netcdf4")
fon="cpc_8-80day.nc"
system("rm -f "+fon)
fo=addfile(fon,"c")
fo->pre8_80  =pre8_80 
end 