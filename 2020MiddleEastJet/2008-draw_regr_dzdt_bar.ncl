;***************************************************************************
;Draw the regression pattern based on the first three eof_ts 
;   1.horizontal distribution of 200hPa level
;   2.horizontal distribution of vertical average 
;   3.vertical section of zonal mean 
;   4.latitude distribution of vertical average and zonal mean
;   5.vertical section of meridional mean 
;   6.longtitude distribution of vertical average and meridional mean
;                                            by Ql_Ren
;                                           2020/02/25
;******************************************************************************
begin
lev  = (/925,850,500,200,100/)
nl = 3
ave_option = 1

area = (/"ME","EA","NA"/);,"GB"
lats = (/ 10 , 10 , 10 /);, 10  ;-90 ,-90 ,-90 ,-90 -20,-20 ,-20
latn = (/ 60 , 60 , 60 /);, 60  ;90 , 90 , 90 , 90  60 , 60 , 60
lonl = (/ 0  ,105 , 250/);, 0   ;30 ;20  0  ,105 , 250
lonr = (/ 70 ,175 , 320/);, 360 ;150;220 70 ,175 , 320
narea = dimsizes(area)
siglvl = 0.05
years = ispan(1979,2015,1)
nyear = dimsizes(years)

case   = "ERA-Interim"
path   = "/home/ys17-19/renql/project/2020MiddleEastJet/"
filename1 = path + "data/" + case + "_DJF_dzdt_in-year.nc"
filename2 = path + "data/" + case + "_DJF_dzdt_in-year.nc"
fig_name = path + "fig/eof_regr_" + lev(nl) + "dzdt_bar"
fig_out    = "pdf" ;"eps";"x11"
pre_case   = (/"(a","(b","(c", "(d","(e","(f", "(g","(h","(i", "(j","(k","(l", "(m","(n","(o", "(p","(q","(r"/)
plot1   = new((/2/),graphic)
plot2   = new((/2/),graphic)
plot3   = new((/2/),graphic)
plot4   = new((/2/),graphic)
plot5   = new((/2/),graphic)
plot6   = new((/2/),graphic)

var_name := (/"dzdt","dudt" ,"dTdt","dbdt"/)  ;
unit     := (/"10~S~-3~N~m~S~2~N~/s~S~3~N~","10~S~-6~N~m/s~S~2~N~","10~S~-6~N~K/s","10~S~-12~N~K/(m*s)"/)  ;this can be used as right string 
nnv = 1 ;used to read dzdt induced by different forcing
draw_var := (/"Qd","vort","Tadv","lQte","lFte","hQte","hFte","ureg"/)  ;,"advect","b3y" this is the variable used to draw, can used as center string
select := (/0,1,2,1,3,0,2,0/) ;(/1,2,3/) ;0 is Qd, 1 is Qd_t
nvar = dimsizes(draw_var)

;------------------EOF time series------------------------
neval = 4
ts    = new((/narea,neval,nyear/),float)
eval  = new((/narea,neval/),float)
pcvar = new((/narea,neval/),float)
sig   = new((/narea,neval/),logical)
do na = 0, narea-1,1
    f  = addfile(path+"data/eof_u_6kinds-"+area(na)+".nc","r")
    eof0 := f->eof0
    eval(na,:)  = eof0@eval(0:(neval-1))
    pcvar(na,:) = eof0@pcvar(0:(neval-1))
    sig(na,:)   = eofunc_north(pcvar(na,:),nyear,False)
    ts(na,:,:)  = f->eof_ts(0, 0:(neval-1), 0:(nyear-1))
end do
ts = ts/conform(ts,sqrt(eval),(/0,1/))
subtitle_color = where(sig,"red","black") 
delete([/f,eof0/])

;------------------ regr uwind  and dudt ------------------------
p_var = new((/narea,neval,nvar/),float)
n_var = new((/narea,neval,nvar/),float)
var  = new((/narea,neval,nvar/),double)
f2  = addfile(path + "/data/ERA-Interim_DJF_uwnd-year.nc","r")
fin = addfile(path + "data/" + case + "_DJF_dzdt_in-year.nc","r") ;vort,fvort,Tadvect
fbv = addfile(path + "data/" + case + "_DJF_dzdt_in3-year.nc","r");Tadvectx,Tadvecty,vortx,vorty
fex = addfile(path + "data/" + case + "_DJF_dzdt_ex-year.nc","r") ;Qd,Qte,Fte
fte = addfile(path + "data/" + case + "_DJF_dzdt_te-year.nc","r") 
do na = 0, narea-1,1
    u := f2->season_ave(:,{200},{lats(na):latn(na)},{lonl(na):lonr(na)})
    nlat = dimsizes(u&lat)
    nlon = dimsizes(u&lon)
    ureg  := regCoef_n(ts(na,:,:),u,1,0)
    utval := onedtond(ureg@tval ,(/neval,nlat,nlon/))
    udf   := onedtond(ureg@nptxy,(/neval,nlat,nlon/))
    copy_VarCoords(u(0,:,:), ureg(0,:,:))
    b := utval 
    b = 0.5
    uprob  := betainc( udf/(udf+utval^2), udf/2.0, b )
    delete([/utval,udf,b/])
    
    vars := new((/nvar-1,nyear,nlat,nlon/),double)
    vars(3:6,:,:,:) = fte->$var_name(nnv)$(select(3:6),:,{lev(nl)},{lats(na):latn(na)},{lonl(na):lonr(na)}) 
    vars(0:2,:,:,:) = fin->$var_name(nnv)$(select(0:2),:,{lev(nl)},{lats(na):latn(na)},{lonl(na):lonr(na)}) 
    vars(1,:,:,:) = dim_sum_n_Wrap(vars(0:1,:,:,:),0) ;vort
    vars(2,:,:,:) = dim_sum_n_Wrap(fbv->$var_name(nnv)$(0:1,:,{lev(nl)},{lats(na):latn(na)},{lonl(na):lonr(na)}),0) ;Tadvect
    vars(0,:,:,:) = fex->$var_name(nnv)$(0,:,{lev(nl)},{lats(na):latn(na)},{lonl(na):lonr(na)}) ;Qd
    vars = vars*1000000

    varr := regCoef_n(ts(na,:,:),vars,1,1)
    tval := onedtond(varr@tval ,(/neval,nvar-1,nlat,nlon/))
    df   := onedtond(varr@nptxy,(/neval,nvar-1,nlat,nlon/))
    lon_m := conform(tval,vars&lon,3)
    copy_VarCoords(vars(0,0,:,:), varr(0,0,:,:))
    b := tval 
    b = 0.5
    prob  := betainc( df/(df+tval^2), df/2.0, b )
    prob   = mask(prob,lon_m.eq.0,False)
    delete([/tval,df,lon_m,b/])
    printVarSummary(varr)

    if(ave_option.eq.0) then 
    pureg = mask(ureg, (uprob.lt.siglvl).and.(ureg.gt.0.0), True)
    nureg = mask(ureg, (uprob.lt.siglvl).and.(ureg.lt.0.0), True)
    var(na,:,nvar-1) = wgt_areaave(pureg,cos(u&lat*4.0*atan(1.0)/180.0),1.0,0) \
                      -wgt_areaave(nureg,cos(u&lat*4.0*atan(1.0)/180.0),1.0,0)
    uprob := conform(varr,uprob,(/0,2,3/))
    ureg  := conform(varr,ureg ,(/0,2,3/))
    pureg := mask(varr, (uprob.lt.siglvl).and.(ureg.gt.0.0), True);.and.(prob.lt.siglvl)
    nureg := mask(varr, (uprob.lt.siglvl).and.(ureg.lt.0.0), True);.and.(prob.lt.siglvl)
    var(na,:,0:(nvar-2)) = wgt_areaave(pureg,cos(u&lat*4.0*atan(1.0)/180.0),1.0,0) \
                          -wgt_areaave(nureg,cos(u&lat*4.0*atan(1.0)/180.0),1.0,0)
    delete([/pureg,nureg,uprob,ureg,prob,varr/])
    printMinMax(var,0)
    printMinMax(var(:,:,0:(nvar-2)),0)
    end if

    if(ave_option.gt.0) then
    pureg = where((uprob.lt.siglvl).and.(ureg.gt.0.0), 1, 0)
    nureg = where((uprob.lt.siglvl).and.(ureg.lt.0.0),-1, 0)
    p_var(na,:,nvar-1) = dim_sum_n(pureg,(/1,2/))
    n_var(na,:,nvar-1) = dim_sum_n(nureg,(/1,2/))
    uprob := conform(varr,uprob,(/0,2,3/))
    ureg  := conform(varr,ureg ,(/0,2,3/))
    pureg := where((uprob.lt.siglvl).and.(ureg.gt.0.0).and.(prob.lt.siglvl).and.(varr.gt.0.0), 1, 0)
    nureg := where((uprob.lt.siglvl).and.(ureg.lt.0.0).and.(prob.lt.siglvl).and.(varr.lt.0.0),-1, 0)
    p_var(na,:,0:(nvar-2)) = dim_sum_n(pureg,(/2,3/))
    n_var(na,:,0:(nvar-2)) = dim_sum_n(nureg,(/2,3/))
    delete([/pureg,nureg,uprob,ureg,prob,varr/])
    printMinMax(p_var,0)
    printMinMax(n_var,0)
    var(na,:,:) = p_var(na,:,:) - n_var(na,:,:)
    var(na,:,:) = 100.0*var(na,:,:)/conform(var(na,:,:),var(na,:,nvar-1),0)
    printMinMax(var,0)
    p_var(na,:,:) = 100.0*p_var(na,:,:)/conform(p_var(na,:,:),p_var(na,:,nvar-1),0)
    n_var(na,:,:) = 100.0*n_var(na,:,:)/abs(conform(n_var(na,:,:),n_var(na,:,nvar-1),0))
    printMinMax(p_var,0)
    printMinMax(n_var,0)
    end if
end do 

;==================================================
;draw figure
;=====================================================
;define the thickness and height of the title,tick and label
   tick_font      = 0.02
   subtitle_font  = 0.02
   title_font     = 0.02
   tick_thick     = 8.0
   subtitle_thick = 8.0  ;not have this attribute
   title_thick    = 5.0
   mp_thick       = 3.0
   font_index     = 22

wks = gsn_open_wks(fig_out,fig_name)

    resxy = True
    resxy@gsnDraw  = False 
    resxy@gsnFrame = False
    resxy@gsnMaximize   = True
    resxy@vpWidthF  = 0.8
    resxy@vpHeightF = 0.4
    resxy@gsnXYBarChart = True
    resxy@gsnXYBarChartBarWidth = 0.2           ; change bar widths

    resxy@trXMinF = 0.5
    resxy@trXMaxF = 8.5
    ;resxy@trYMinF = floor(min(var)) 
    ;resxy@trYMaxF =  ceil(max(var))
    ;resxy@tiYAxisString = unit(nnv)
    resxy@tiYAxisFont   = font_index
    resxy@gsnYRefLine   = 0.0 
    resxy@tmXBMode   = "Explicit"         ; explicit labels
    resxy@tmXBValues = (/1,2,3,4,5,6,7,8/)
    resxy@tmXBLabels = draw_var 

    resxy@tmBorderThicknessF  = mp_thick
    resxy@tmXBLabelFont       = font_index
    resxy@tmYLLabelFont       = font_index
    resxy@tmXBLabelFontHeightF = tick_font
    resxy@tmXTLabelFontHeightF = tick_font
    resxy@tmYLLabelFontHeightF = tick_font

    resxy@tmXBMinorOn          = False ; Turn off x-axis (bottom) minor ticks
    resxy@tmXTMinorOn          = False ; Turn off x-axis (tottom)minor ticks
    resxy@tmXTOn               = False ; Turn off x-axis (top) major ticks
    resxy@tmYRMinorOn          = False ; Turn off y-axis (right) minor ticks
    resxy@tmYLMinorOn          = False ; Turn off y-axis (left) minor ticks
    resxy@tmYROn               = False ; Turn off y-axis (right) major ticks
    resxy@tmXMajorGrid                   = True
    resxy@tmXMajorGridLineDashPattern    = 2
    resxy@tmXMajorGridLineColor          = "gray"
    resxy@tmXMajorGridThicknessF         = 2.0
    
    resxy@gsnStringFont         = font_index
    resxy@gsnStringFontHeightF   = subtitle_font
    resxy@gsnCenterString  = "";lev(nl) + " hPa" 
    resxy@gsnLeftString    = ""
    resxy@gsnRightString   = ""

    lbres                    = True          ; labelbar only resources
    lbres@vpWidthF           = 0.25           ; labelbar width
    lbres@vpHeightF          = 0.1           ; labelbar height
    lbres@lbBoxMajorExtentF  = 0.2          ; puts space between color boxes
    ;lbres@lbFillColors       = (/"red","green","blue"/)
    lbres@lbMonoFillPattern  = True          ; Solid fill pattern
    lbres@lbLabelFontHeightF = tick_font         ; font eight. default is small
    lbres@lbLabelJust        = "CenterLeft"  ; left justify labels
    lbres@lbPerimOn          = False
    lbres@lgPerimColor      = "white"

    resp  = True    
    resp@gsnMaximize   = True
    resp@gsnFrame      = False  ; Don't advance frame.
    resp@gsnPaperOrientation = "landscape";"portrait" ;
    resp@gsnPanelLabelBar = False
    resp@gsnPanelBottom = 0.13   ; Leave room for labelbar
    resp@txFontHeightF = title_font 
    resp@txFontThicknessF = title_thick

do ne = 0,1,1
    np = ne
    nev = ne + 1
    resxy@gsnLeftString  = pre_case(np) + ") EOF" + nev

    if(ave_option.le.1) then 
    resxy@trYMinF = floor(min(var)) 
    resxy@trYMaxF =  ceil(max(var))
    resxy@tiYAxisString = unit(nnv)
    resxy@gsnXYBarChartColors = (/"red"/)
    xstart = resxy@trXMinF + (1-resxy@gsnXYBarChartBarWidth*3)/2
    print("xstart = " + xstart)
    plot1(np) = gsn_csm_xy(wks,fspan(xstart,xstart+7.0,8),var(0,ne,:),resxy)     ; draw each time series
    
    resxy@gsnXYBarChartColors = (/"green"/)                    ; seperately, not
    xstart = resxy@trXMinF + (1-resxy@gsnXYBarChartBarWidth*3)/2 + resxy@gsnXYBarChartBarWidth
    print("xstart = " + xstart)
    plot2(np) = gsn_csm_xy(wks,fspan(xstart,xstart+7.0,8),var(1,ne,:),resxy)     ; draw each time series
    overlay(plot1(np),plot2(np))
    
    resxy@gsnXYBarChartColors = (/"blue"/)                   ; but tweaking where
    xstart = resxy@trXMinF + (1-resxy@gsnXYBarChartBarWidth*3)/2 + resxy@gsnXYBarChartBarWidth*2
    print("xstart = " + xstart)
    plot3(np) = gsn_csm_xy(wks,fspan(xstart,xstart+7.0,8),var(2,ne,:),resxy)     ; draw each time series
    overlay(plot1(np),plot3(np))
    else
    resxy@trXMaxF = 7.5
    resxy@tmXBValues := (/1,2,3,4,5,6,7/)
    resxy@tmXBLabels := draw_var(0:(nvar-2)) 
    resxy@trYMinF = floor(min(n_var(:,ne,0:(nvar-2)))) 
    resxy@trYMaxF =  ceil(max(p_var(:,ne,0:(nvar-2))))
    resxy@tiYAxisString = "%" 
    resxy@gsnXYBarChartColors = (/"red"/)
    xstart = resxy@trXMinF + (1-resxy@gsnXYBarChartBarWidth*3)/2
    plot1(np) = gsn_csm_xy(wks,fspan(xstart,xstart+6.0,7),p_var(0,ne,0:(nvar-2)),resxy)     ; draw each time series
    plot4(np) = gsn_csm_xy(wks,fspan(xstart,xstart+6.0,7),n_var(0,ne,0:(nvar-2)),resxy)     ; draw each time series
    overlay(plot1(np),plot4(np))
    
    resxy@gsnXYBarChartColors = (/"green"/)                    ; seperately, not
    xstart = resxy@trXMinF + (1-resxy@gsnXYBarChartBarWidth*3)/2 + resxy@gsnXYBarChartBarWidth
    plot2(np) = gsn_csm_xy(wks,fspan(xstart,xstart+6.0,7),p_var(1,ne,0:(nvar-2)),resxy)     ; draw each time series
    plot5(np) = gsn_csm_xy(wks,fspan(xstart,xstart+6.0,7),n_var(0,ne,0:(nvar-2)),resxy)     ; draw each time series
    overlay(plot1(np),plot2(np))
    overlay(plot1(np),plot5(np))
    
    resxy@gsnXYBarChartColors = (/"blue"/)                   ; but tweaking where
    xstart = resxy@trXMinF + (1-resxy@gsnXYBarChartBarWidth*3)/2 + resxy@gsnXYBarChartBarWidth*2
    plot3(np) = gsn_csm_xy(wks,fspan(xstart,xstart+6.0,7),p_var(2,ne,0:(nvar-2)),resxy)     ; draw each time series
    plot6(np) = gsn_csm_xy(wks,fspan(xstart,xstart+6.0,7),n_var(0,ne,0:(nvar-2)),resxy)     ; draw each time series
    overlay(plot1(np),plot3(np))
    overlay(plot1(np),plot6(np))
    end if
end do

gsn_panel(wks,plot1,(/2,1/),resp)

colors = (/"red","green", "blue"/)
xpos   = (/0.2  ,    0.4,   0.6/)
do na = 0,narea-1,1
lbres@lbFillColors = colors(na)
gsn_labelbar_ndc(wks,1,area(na),xpos(na),0.18,lbres)   ; draw right labelbar column
end do

frame(wks)
end
