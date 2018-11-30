;************************************************
;advection_explicit_schemes_exercise6.ncl
;
;This is the sixth homework of Atmospheric Modeling
;FTCS, Lax, and Lax-Wendroff schemes are used to difference the 
;advection equation. Periodic boundary conditions are used.
;Then plots the result.
;
;                                       by Renql
;                                      2018/11/29
;************************************************

begin
;---Define parameters------------------
c  = 1.0    ;phase velocity
dx = 0.02   ;grid spacing
dt = 0.002  ;time spacing
lx = 1.0    ;the length of x
nx = 50     ;the number of grids
nstp = toint(lx/(c*dt)+1)    ;integration steps
drawt = (/0, 0.1, 0.5, 1.0/)    ;the time used to draw
title = (/"FTCS","Lax","Lax-Wendroff"/) + "_scheme"
opt   = 0 ;0 mean FTCS, 1 mean Lax, 2 mean Lax-Wendroff

;---Boundary and initial value--------------
k     = 2*3.14/dx
theta = 1
x0    = 0
x = fspan(-(lx-dx)/2, (lx-dx)/2, nx)

f  = new((/nstp,nx/),float)
f!0 = "time"
f&time = fspan(0, lx/c, nstp)
f(0,:) = cos(k*(x-x0))*exp(-(x-x0)^2/(2*theta^2))

f1 = new((/dimsizes(drawt),nx/),float)
do nt = 0, dimsizes(drawt)-1, 1
    f1(nt,:) = cos(k*(x-c*drawt(nt)-x0))*exp(-(x-c*drawt(nt)-x0)^2/(2*theta^2))
end do
print(x(0:10)) 
print(f(0,0:10))

;---Explicit schemes for the advection equation---------
if(opt.eq.0) then
    do istp = 1, nstp-1, 1
    do ix = 0, nx-1, 1
        if(ix.eq.0) then
            f(istp,ix) = f(istp-1,ix) - c*dt/2/dx*( f(istp-1,ix+1) - f(istp-1,nx-1))
        else if(ix.eq.(nx-1)) then
            f(istp,ix) = f(istp-1,ix) - c*dt/2/dx*( f(istp-1,0) - f(istp-1,ix-1))
        else
            f(istp,ix) = f(istp-1,ix) - c*dt/2/dx*( f(istp-1,ix+1) - f(istp-1,ix-1))
        end if
        end if
    end do
    end do

else if(opt.eq.1) then
    do istp = 1, nstp-1, 1
    do ix = 0, nx-1, 1
        if(ix.eq.0) then
            f(istp,ix) = (0.5 - c*dt/2/dx)*f(istp-1,ix+1) + (0.5 + c*dt/2/dx)*f(istp-1,nx-1)
        else if(ix.eq.(nx-1)) then
            f(istp,ix) = (0.5 - c*dt/2/dx)*f(istp-1,0) + (0.5 + c*dt/2/dx)*f(istp-1,ix-1)
        else
            f(istp,ix) = (0.5 - c*dt/2/dx)*f(istp-1,ix+1) + (0.5 + c*dt/2/dx)*f(istp-1,ix-1)
        end if
        end if
    end do
    end do

else
    do istp = 1, nstp-1, 1
    do ix = 0, nx-1, 1
        if(ix.eq.0) then
            f(istp,ix) = (1-(c*dt/dx)^2)*f(istp-1,ix) + c*dt/2/dx*( (c*dt/dx-1)*f(istp-1,ix+1) + (c*dt/dx+1)*f(istp-1,nx-1))
        else if(ix.eq.(nx-1)) then
            f(istp,ix) = (1-(c*dt/dx)^2)*f(istp-1,ix) + c*dt/2/dx*( (c*dt/dx-1)*f(istp-1,0) + (c*dt/dx+1)*f(istp-1,ix-1))
        else
            f(istp,ix) = (1-(c*dt/dx)^2)*f(istp-1,ix) + c*dt/2/dx*( (c*dt/dx-1)*f(istp-1,ix+1) + (c*dt/dx+1)*f(istp-1,ix-1))
        end if
        end if
    end do
    end do
end if
end if

;************************************************
; plotting
;************************************************
;define the thickness and height of the title,tick and label
   tick_font      = 0.015
   subtitle_font  = 0.018
   title_font     = 0.03
   label_font     = 0.015
   tick_thick     = 5.0
   subtitle_thick = 5.0  ;not have this attribute
   title_thick    = 5.0
   label_thick    = 1.5
   mp_thick       = 4.0
   font_index     = 22

figname = title(opt) + "_k" + k + "_theta" + theta + "_dt" + dt
plot = new((/dimsizes(drawt)/),graphic)
wks  = gsn_open_wks("png",figname)      ; send graphics to PNG file
  
  res               = True                        ; plot mods desired
  res@gsnFrame      = False
  res@gsnDraw       = False
  res@vpHeightF     = 0.2
  res@vpWidthF      = 0.8
  res@gsnMaximize   = True
  res@tiXAxisString = "x"  ; xaxis
  res@tiYAxisString = "f"  ; xaxis
  res@tiYAxisAngleF = 0.0 ;let the yaxis horizontal
  
  res@xyDashPatterns    = (/0, 0/)
  res@xyLineThicknesses = (/4.0, 4.0/)   ; make second line thicker
  res@xyLineColors      = (/"black","red"/)

;  res@gsnYRefLine             = 0.0
;  res@gsnYRefLineDashPatterns = 11   
;  res@gsnYRefLineThicknesses  = 4.0
;  res@gsnYRefLineColors       = "black"

  res@tmBorderThicknessF  = mp_thick
  res@tmXBLabelFont         = font_index
  res@tmYLLabelFont         = font_index
  res@tmXBLabelFontHeightF     = tick_font 
  res@tmXBLabelFontThicknessF  = tick_thick
  res@tmYLLabelFontHeightF     = tick_font
  res@tmYLLabelFontThicknessF  = tick_thick
  
  res@tmXBMajorThicknessF     = mp_thick
  res@tmYLMajorThicknessF     = mp_thick
  res@tmXBMinorThicknessF     = mp_thick
  res@tmYLMinorThicknessF     = mp_thick
  
  res@tmXMajorGrid                = True          ; implement x grid 
  res@tmXMajorGridThicknessF      = 2.0           ; 2.0 is default
  res@tmXMajorGridLineColor       = "gray"
  res@tmXMajorGridLineDashPattern = 2             ; select short dash lines
  
  res@tmYMajorGrid                = True          ; implement y grid
  res@tmYMajorGridThicknessF      = 2.0           ; 2.0 is default
  res@tmYMajorGridLineColor       = "gray"
  res@tmYMajorGridLineDashPattern = 2             ; select short dash lines
  
  res@pmLegendDisplayMode = "Always"
  res@lgOrientation       = "Vertical"
  res@lgPerimOn           = False
  res@lgPerimFill         = "SolidFill"
  res@lgPerimFillColor    = "white"
  res@xyExplicitLegendLabels = (/"True","Numerical"/)
  res@lgLabelFontColor       = "black"
  res@lgLabelFontThicknessF  = label_thick
  res@lgLabelFontHeightF     = label_font
  res@pmLegendWidthF   = 0.08
  res@pmLegendHeightF  = 0.07
  res@pmLegendZone  = 0   ;set the positional origin in the center of the plot
  res@pmLegendOrthogonalPosF  = 0.4   ; move ref vector along Y-axis
  res@pmLegendParallelPosF    = 0.3    ; move ref vector along X-axis

  res@gsnStringFont        = font_index
  res@gsnStringFontHeightF = subtitle_font

  resp  = True    
  resp@gsnMaximize   = True
  resp@gsnPaperOrientation = "portrait"
  resp@txFontHeightF = title_font 
  resp@txFontThicknessF = title_thick
  resp@txString  = title(opt)

  res@gsnLeftString  = "dt = " + dt + "s" 
do nt = 0, dimsizes(drawt)-1, 1
  res@gsnRightString = "t = " + drawt(nt) + "s"
  plot(nt) = gsn_csm_xy(wks,x,(/f1(nt,:),f(toint(drawt(nt)/dt),:)/),res)   ; create plot  
end do

  gsn_panel(wks,plot,(/dimsizes(drawt),1/),resp)
end


