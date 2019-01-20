;***************************************************************************
;filename:      spec.ncl
;read the 27 years U,V,T at TP areas in three cases (CTRL, NUDG6h, NUDG24h) and then 
;calculate the season and area average in order to test the right of the 
;NUDG experiment.
;                                            by Ql_Ren
;                                           2019/01/20
;******************************************************************************
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/shea_util.ncl"

begin
path    = "/users/yangsong3/renql/project/TP_NUDG"
filein  = path +(/"/data/AMIP-CTRL/AMIP_C5PM.cam.h1.YEAR.",\
                  "/data/TP-NUDG-24h/AMIP_C5PM_TP_NUDG.cam.h1.YEAR.",\
                  "/data/TP-NUDG-6h/AMIP_C5PM_TP_NUDG.cam.h1.YEAR."/)
var_name = (/"U","V","T"/)
case     = (/"CTRL","NG24h","NG6h"/)
years    = ispan(1979,2005,1)   ;the calculate time

;area of TP 
lats = 25 
latn = 40
lonl = 70
lonr = 105
lev  = (/850,500,200/)

nv = 2
nl = 2
nc = 0

fig_name  = path + "/20190120-spectral/fig/" + case(nc) + " " + var_name(nv) + " " + lev(nl)
fig_title = case(nc) + " " + var_name(nv) + " " + lev(nl)

;---------------------------------------------------------------
;read the data, then calculate the area and season average
;--------------------------------------------------------------
files   = systemfunc("ls " + filein(nc) + "*.daily." + var_name(nv) + ".nc")
f       = addfiles(files,"r")
ListSetType(f,"cat")
var :=  f[:]->$var_name(nv)$(:,{lev(nl)},{lats:latn},{lonl:lonr})
var := wgt_areaave_Wrap(var,1.0,1.0,0) ;get var(nvar,ncase,nyear,365,nlev)

; set spectral function arguments
  d   = 0    ; detrending opt: 0=>remove mean 1=>remove mean + detrend
  sm  = 21   ; smooth: should be at least 3 and odd
  pct = 0.10 ; percent taper: (0.0 <= pct <= 1.0) 0.10 common. 

; calculate spectrum
  sdof = specx_anal(var,d,sm,pct)

; calculate confidence interval [here 5 and 95%]
; return 4 curves to be plotted
  splt = specx_ci (sdof, 0.05, 0.95)
  printVarSummary(splt)                          ; splt(4,702)

  f = sdof@frq                                 ; convenience
  p   = 1/f                                    ; *highly non-linear*
  p!0 = "f"
  p&f = f
  p@long_name = "period"
  p@units     = "month"

;************************************************
; plotting
;************************************************
   wks  = gsn_open_wks("png",fig_name)             ; send graphics to PNG file

   res = True                                      ; plot mods desired
   res@tiMainString = fig_title 
   res@tiXAxisString = "Frequency (cycles/month)"  ; xaxis
   res@tiYAxisString = "Variance/frq_interval"     ; yaxis

;***********************************************
; Generate log plot showing "red noise" confidence bounds
; (a) log scaling and (b) the Band Width
;***********************************************

   res@trYMinF             = 0.00                 ; manually set lower limit
   res@trYMaxF             = 20.0                 ;   "          upper
   res@gsnFrame            = False                ; do not advance frame

; PLOT abscissa as period (1/frequency)
 print("====")
 print(f+"   "+p+"  "+splt(0,:) )             ; print 
 print("====")

   res@tiXAxisString = "Period (month)"          ; xaxis
   res@tiYAxisString = "Variance/freq"           ; yaxis

   ip   = ind(p.le.240)                        ; all indices for "short" periods
   plot = gsn_csm_xy(wks,p(ip), splt(:,ip),res)     

   res@trXReverse = True
   plot = gsn_csm_xy(wks,p(ip), splt(:,ip),res)    
frame(wks)
end
