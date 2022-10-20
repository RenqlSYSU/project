;***************************************************************************
;Draw the regression pattern based on the first three eof_ts 
;   1.horizontal distribution of 200hPa level
;                                            by Ql_Ren
;                                           2020/02/25
;******************************************************************************
function dof_breth99(ts,dim)
begin
    dims = dimsizes(ts)
    n = dims(dim)
    lagr = esacr_n(ts,2,dim)
    ;print("lag autocorrelation:"+lagr)
    dof = n*(1-lagr(1)*lagr(2))/(1+lagr(1)*lagr(2))
    return dof
end

function dof_metz91(ts1,ts2)
begin
    n = dimsizes(ts1)
    lag = toint(n/2.0)-1
    lagr1 = esacr_n(ts1,lag,0)
    lagr2 = esacr_n(ts2,lag,0)
    term = 1.0
    do nlag = 1, lag, 1
        term = term + tofloat(2*(1-lag/n)*lagr1(lag)*lagr2(lag))
    end do
    
    if (term.gt.1) then
        dof = toint(n/term)
    else
        dof = n
    end if
    return dof
end

begin
lev  = (/925,850,500,200,100/)
nl = 3
filter = False;True ;

area = (/"ME","EA","NA"/);,"GB"
if (.not.isvar("lats")) then
lats = (/ 15 , 20 , 25 /);, 10  ;-90 ,-90 ,-90 ,-90 
latn = (/ 40 , 45 , 50 /);, 60  ;90 , 90 , 90 , 90 
lonl = (/ 20 ,105 , 250/);, 0   ;30 ;20
lonr = (/ 60 ,175 , 320/);, 360 ;150;220
end if
narea = dimsizes(area)
years = ispan(1979,2015,1)
nyear = dimsizes(years)

case   = "ERA-Interim"
path   = "/home/ys17-19/renql/project/2020MiddleEastJet/"

var_name := (/"dzdt","dudt" ,"dTdt","dbdt"/)  ;
unit     := (/"10~S~-3~N~m~S~2~N~/s~S~3~N~","10~S~-6~N~m/s~S~2~N~","10~S~-6~N~K/s","10~S~-12~N~K/(m*s)"/)  ;this can be used as right string 
nnv = 1 ;used to read dzdt induced by different forcing
draw_var := (/"Qd","lQte","lFte","hQte","hFte","Tadv","vort"/)  ;,"advect","b3y" this is the variable used to draw, can used as center string
select := (/0,1,3,0,2,0,1/) ;(/1,2,3/) ;0 is Qd, 1 is Qd_t
nvar = dimsizes(draw_var)

;------------------EOF time series------------------------
neval = 4
fin = addfile(path + "data/" + case + "_DJF_dzdt_b1xb1yb2-year.nc","r") ;vort,fvort,Tadvect
ftv = addfile(path + "data/" + case + "_DJF_dzdt_in3-year.nc","r");Tadvectx,Tadvecty,vortx,vorty
fex = addfile(path + "data/" + case + "_DJF_dzdt_ex-year.nc","r") ;Qd,Qte,Fte
fte = addfile(path + "data/" + case + "_DJF_dzdt_te-year.nc","r") 

corr2 = new((/nvar,6/),double)
prob2 = new((/nvar,6/),double)
do na = 0, narea-1,1
    f = addfile(path+"data/eof2_u_6kinds-DJF_"+area(na)+".nc","r")
    eof0 := f->eof0
    ;printVarSummary(eof0)
    eval  = eof0@eval(0:(neval-1))
    ts = f->eof_ts(0, 0:(neval-1), 0:(nyear-1))
    ts = ts/conform(ts,sqrt(eval),(/0/))
    nlat = dimsizes(eof0&lat)
    nlon = dimsizes(eof0&lon)
    
    vars = new((/nvar,nyear,nlat,nlon/),double)
    vars(1:4,:,:,:) = fte->$var_name(nnv)$(select(1:4),:,{lev(nl)},{lats(na):latn(na)},{lonl(na):lonr(na)}) 
    vars(0,:,:,:) = fex->$var_name(nnv)$(0,:,{lev(nl)},{lats(na):latn(na)},{lonl(na):lonr(na)}) ;Qd
    vars(5,:,:,:) = dim_sum_n_Wrap(ftv->$var_name(nnv)$(0:1,:,{lev(nl)},{lats(na):latn(na)},{lonl(na):lonr(na)}),0) ;Tadvect
    vars(6,:,:,:) = dim_sum_n_Wrap(fin->$var_name(nnv)$(:,:,{lev(nl)},{lats(na):latn(na)},{lonl(na):lonr(na)}),0) ;vort
    vars = vars*1000000
    if(filter) then
    remove_mean = True
    return_info = True
    var := dtrend_msg_n(years,vars,remove_mean,return_info,1)

    ca = 9.0 ;year
    fca = 1.0/ca
    fcf = fspan(0, 0.5, nyear/2) 
    low_cf = ind(fcf.lt.fca)

    cf = ezfftf_n(var,1)
    term = cf 
    term(:,:,low_cf,:,:) = 0.0  
    vars = ezfftb_n(term,0,2)
    delete([/term,cf/])
    else
    var := vars
    end if

    optETS = True
    optETS@jopt = 0; Use the standardized data matrix to compute the time series. The default is to use data and evec
    ;the results of these two option is a little different, but it didn't matter much 

    do nv = 0, nvar-1, 1
    var_ts = eofunc_ts_n(vars(nv,:,:,:), eof0, optETS, 0) ;get (nvar,neval,nyear)
    do ne = 0, 1, 1
        dof1 = dof_breth99(ts(ne,:),0) 
        dof2 = dof_breth99(var_ts(ne,:),0)
        edof = toint((dof1+dof2)/2.0)
        print("number="+nyear+", var_ts effective dof="+dof2+", eof_ts effective dof="+dof1+", average="+edof)
        ;edof = dof_metz91(var_ts(ne,:), ts(ne,:))
        ;print("dof metz91: "+edof)
        corr  = esccr(var_ts(ne,:), ts(ne,:), 0) ;the last number is mxlag
        prob  = rtest(corr,edof,0) ;the last number is currently not used
        corr2(nv,na*2+ne) = corr
        prob2(nv,na*2+ne) = prob
        if(prob.le.0.1) then
        print(area(na) + ne + " eval=" + eval(ne) + " " + draw_var(nv) + " corr : " + corr + " ; prob : " + prob + " significant")
        else
        print(area(na) + ne + " eval=" + eval(ne) + " " + draw_var(nv) + " corr : " + corr + " ; prob : " + prob)
        end if
    end do
    print(" ")
    end do

delete([/f,eof0,vars/])
print(" ")
end do

opt      = True
opt@fout = path+"project_dzdt_corr.txt"
opt@row  = False 
write_matrix(corr2,"6f8.3",opt)     ;data只能是二维的，fmtf是打印的格式如"9f12.6"
opt@fout = path+"project_dzdt_prob.txt"
write_matrix(prob2,"6f8.3",opt)     ;data只能是二维的，fmtf是打印的格式如"9f12.6"
end

