;***************************************************************************
;Draw the regression pattern based on the first three eof_ts 
;   1.horizontal distribution of 200hPa level
;                                            by Ql_Ren
;                                           2020/02/25
;******************************************************************************
begin
lev  = (/925,850,500,200,100/)
nl = 3

area = (/"ME","EA","NA"/);,"GB"
lats = (/ 15 , 20 , 25 /);, 10  ;-90 ,-90 ,-90 ,-90 
latn = (/ 40 , 45 , 50 /);, 60  ;90 , 90 , 90 , 90 
lonl = (/ 20 ,110 , 250/);, 0   ;30 ;20
lonr = (/ 60 ,180 , 320/);, 360 ;150;220
narea = dimsizes(area)
years = ispan(1979,2015,1)
nyear = dimsizes(years)

case   = "ERA-Interim"
path   = "/home/ys17-19/renql/project/2020MiddleEastJet/"

var_name := (/"dzdt","dudt" ,"dTdt","dbdt"/)  ;
unit     := (/"10~S~-3~N~m~S~2~N~/s~S~3~N~","10~S~-6~N~m/s~S~2~N~","10~S~-6~N~K/s","10~S~-12~N~K/(m*s)"/)  ;this can be used as right string 
nnv = 1 ;used to read dzdt induced by different forcing
draw_var := (/"Qd","vort","Tadv","lQte","lFte","hQte","hFte"/)  ;,"advect","b3y" this is the variable used to draw, can used as center string
select := (/0,1,2,1,3,0,2/) ;(/1,2,3/) ;0 is Qd, 1 is Qd_t
nvar = dimsizes(draw_var)

;------------------EOF time series------------------------
neval = 4
f2  = addfile(path + "/data/ERA-Interim_DJF_uwnd-year.nc","r")
fin = addfile(path + "data/" + case + "_DJF_dzdt_in-year.nc","r") ;vort,fvort,Tadvect
fbv = addfile(path + "data/" + case + "_DJF_dzdt_in3-year.nc","r");Tadvectx,Tadvecty,vortx,vorty
fex = addfile(path + "data/" + case + "_DJF_dzdt_ex-year.nc","r") ;Qd,Qte,Fte
fte = addfile(path + "data/" + case + "_DJF_dzdt_te-year.nc","r") 

do na = 0, narea-1,1
    f = addfile(path+"data/eof_u_6kinds-"+area(na)+".nc","r")
    eof0 := f->eof0
    ;printVarSummary(eof0)
    eval  = eof0@eval(0:(neval-1))
    ts = f->eof_ts(0, 0:(neval-1), 0:(nyear-1))
    ;ts = ts/conform(ts,sqrt(eval),(/0/))
    nlat = dimsizes(eof0&lat)
    nlon = dimsizes(eof0&lon)
    
    vars = new((/nvar,nyear,nlat,nlon/),double)
    vars(3:6,:,:,:) = fte->$var_name(nnv)$(select(3:6),:,{lev(nl)},{lats(na):latn(na)},{lonl(na):lonr(na)}) 
    vars(0:2,:,:,:) = fin->$var_name(nnv)$(select(0:2),:,{lev(nl)},{lats(na):latn(na)},{lonl(na):lonr(na)}) 
    vars(1,:,:,:) = dim_sum_n_Wrap(vars(0:1,:,:,:),0) ;vort
    vars(2,:,:,:) = dim_sum_n_Wrap(fbv->$var_name(nnv)$(0:1,:,{lev(nl)},{lats(na):latn(na)},{lonl(na):lonr(na)}),0) ;Tadvect
    vars(0,:,:,:) = fex->$var_name(nnv)$(0,:,{lev(nl)},{lats(na):latn(na)},{lonl(na):lonr(na)}) ;Qd
    vars = vars*1000000

    optETS = True
    optETS@jopt = 1; Use the standardized data matrix to compute the time series. The default is to use data and evec
    ;the results of these two option is a little different, but it didn't matter much 

    do ne = 0, 1, 1
    do nv = 0, nvar-1, 1
    var_ts = eofunc_ts_n(vars(nv,:,:,:), eof0, optETS, 0) ;get (nvar,neval,nyear)
    corr  = esccr(var_ts(ne,:), ts(ne,:), 0) ;the last number is mxlag
    prob  = rtest(corr,nyear,0) ;the last number is currently not used
    print(area(na) + ne + " " + draw_var(nv) + " corr : " + corr + " ; prob : " + prob)
    end do
    end do

delete([/f,eof0,vars/])
end do

end

