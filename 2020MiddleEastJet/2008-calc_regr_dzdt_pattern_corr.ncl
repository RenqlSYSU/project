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
lev  = (/1000, 925, 850, 700, 600, 500, 400, 350, 300, 250, 200, 175, 150, 125, 100/)
area = (/"ME","EA","NA"/);,"GB"
lats = (/ 10 , 10 , 10 /);, 10  ;-90 ,-90 ,-90 ,-90 -20,-20 ,-20
latn = (/ 60 , 60 , 60 /);, 60  ;90 , 90 , 90 , 90  60 , 60 , 60
lonl = (/ 0  ,105 , 250/);, 0   ;30 ;20  0  ,105 , 250
lonr = (/ 70 ,175 , 320/);, 360 ;150;220 70 ,175 , 320
narea = dimsizes(area)
siglvl = 0.05
years = ispan(1979,2015,1)
nyear = dimsizes(years)
nlev  = dimsizes(lev)

case   = "ERA-Interim"
path   = "/home/ys17-19/renql/project/2020MiddleEastJet/"

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
f2  = addfile(path + "/data/ERA-Interim_DJF_uwnd-year.nc","r")
fin = addfile(path + "data/" + case + "_DJF_dzdt_in-year.nc","r") ;vort,fvort,Tadvect
fbv = addfile(path + "data/" + case + "_DJF_dzdt_in3-year.nc","r");Tadvectx,Tadvecty,vortx,vorty
fex = addfile(path + "data/" + case + "_DJF_dzdt_ex-year.nc","r") ;Qd,Qte,Fte
fte = addfile(path + "data/" + case + "_DJF_dzdt_te-year.nc","r") 
do na = 0, narea-1,1
    u := f2->season_ave(:,{lev},{lats(na):latn(na)},{lonl(na):lonr(na)})
    nlat = dimsizes(u&lat)
    nlon = dimsizes(u&lon)
    ureg  := regCoef_n(ts(na,:,:),u,1,0)
    utval := onedtond(ureg@tval ,(/neval,nlev,nlat,nlon/))
    udf   := onedtond(ureg@nptxy,(/neval,nlev,nlat,nlon/))
    copy_VarCoords(u(0,:,:,:), ureg(0,:,:,:))
    b := utval 
    b = 0.5
    uprob  := betainc( udf/(udf+utval^2), udf/2.0, b )
    delete([/utval,udf,b/])
    
    vars := new((/nvar-1,nyear,nlev,nlat,nlon/),double)
    vars(3:6,:,:,:,:) = fte->$var_name(nnv)$(select(3:6),:,{lev},{lats(na):latn(na)},{lonl(na):lonr(na)}) 
    vars(0:2,:,:,:,:) = fin->$var_name(nnv)$(select(0:2),:,{lev},{lats(na):latn(na)},{lonl(na):lonr(na)}) 
    vars(1,:,:,:,:) = dim_sum_n_Wrap(vars(0:1,:,:,:,:),0) ;vort
    vars(2,:,:,:,:) = dim_sum_n_Wrap(fbv->$var_name(nnv)$(0:1,:,{lev},{lats(na):latn(na)},{lonl(na):lonr(na)}),0) ;Tadvect
    vars(0,:,:,:,:) = fex->$var_name(nnv)$(0,:,{lev},{lats(na):latn(na)},{lonl(na):lonr(na)}) ;Qd
    vars = vars*1000000

    varr := regCoef_n(ts(na,:,:),vars,1,1)
    tval := onedtond(varr@tval ,(/neval,nvar-1,nlev,nlat,nlon/))
    df   := onedtond(varr@nptxy,(/neval,nvar-1,nlev,nlat,nlon/))
    lon_m := conform(tval,vars&lon,4)
    copy_VarCoords(vars(0,0,:,:,:), varr(0,0,:,:,:))
    b := tval 
    b = 0.5
    prob  := betainc( df/(df+tval^2), df/2.0, b )
    prob   = mask(prob,lon_m.eq.0,False)
    delete([/tval,df,lon_m,b/])
    printVarSummary(varr)

    ;The centred (uncentred) measures the similarity of two patterns after (without) removal of the global mean.
    opt    := 1 ;0 = centered pattern correlation; 1 = uncentered pattern correlation.
    weight := conform(u(0,:,:,:), cos(u&lat*4.0*atan(1.0)/180.0),(/1/)) 
    uprobm := conform(tofloat(varr), uprob,(/0,2,3,4/))
    uregm  := conform(tofloat(varr), ureg ,(/0,2,3,4/))
    ;uregm  := mask(uregm, uprobm.lt.siglvl, True)
    ;varrm  := mask(varr , uprobm.lt.siglvl, True)
    varrm  := varr
    printVarSummary(uregm)
    uregm  := reshape(uregm ,(/neval,nvar-1,nlev*nlat,nlon/))
    varrm  := reshape(varrm ,(/neval,nvar-1,nlev*nlat,nlon/))
    weight := reshape(weight,(/nlev*nlat,nlon/))
    spcor  := pattern_cor(uregm,tofloat(varrm),weight,opt) 
    print(area(na))
    print(spcor)
    delete([/uregm,spcor/])
    delete([/varrm,uprobm,weight/])
end do

end 

