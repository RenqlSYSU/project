;***************************************************************************
;use the reanalysis data to calculate diabatic heating, transient eddy heating
;and Transient eddy vorticity forcing 
;                                            by Ql_Ren
;                                           2020/02/04
;******************************************************************************
begin
lats =  -30;-20
latn =  90 ;60
lonl =  0  ;20
lonr = 360 ;220
lev  = (/1000,925,850,700,600, 500,400,350,300,250, 200,175,150,125,100/)

path1 = "/home/ys17-19/data-model/CMIP5/"
path2 = "/home/ys17-19/renql/project/2020SouthTrough/"

filepath  = systemfunc("find " + path1 + " -name zg_Amon_*_historical_r1i1p1_*.nc")
modelname = str_get_field(filepath,5,"/")
nfile     = dimsizes(filepath)
;print(modelname(12))

numb  = 0
do nf = 0, nfile-2, 1
if(modelname(nf).eq.modelname(nf+1)) then 
    numb = numb + 1
else
    print("model name :" + modelname(nf))
    f = addfile
    numb  = 0

    



