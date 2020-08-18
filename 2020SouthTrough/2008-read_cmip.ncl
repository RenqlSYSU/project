;***************************************************************************
;Find file and know its resolution variable condition 
;
;                                            by Ql_Ren
;                                           2020/08/17
;******************************************************************************
begin
lats =  -30;-20
latn =  90 ;60
lonl =  0  ;20
lonr = 360 ;220
lev  = (/1000,925,850,700,600, 500,400,350,300,250, 200,175,150,125,100/)

path1 = "/home/ys17-19/data-model/AMIP/ -name zg_Amon_*_r1i1p1_*.nc"
;path1 = "/home/ys17-19/data-observation/cmip6/"
;path1 = "/home/ys17-19/data-model/CMIP5/ -name zg_Amon_*_historical_r1i1p1_*.nc"
path2 = "/home/ys17-19/renql/project/2020SouthTrough/"
fileinfor = True ; False
test      = False

filepath  = systemfunc("find " + path1) 
nfile     = dimsizes(filepath)
nfield    = str_fields_count(filepath,"/")

modelname = new((/nfile/),string)
year      = new((/nfile/),string)
do nf = 0, nfile-1, 1
    filename = str_split(str_get_field(filepath(nf),nfield(nf),"/"),"_")
    modelname(nf) = filename(2)
    year(nf) = stringtoint(str_get_cols(filename(5),0,3))
end do
print(year(0:10))
print(modelname(0:10))

if(test) then 
filepath1 = filepath(0:2)
yearindex = dim_pqsort(year(0:2),1)
fs = addfiles(filepath1(yearindex),"r")
ListSetType(fs,"cat")
var = fs[:]->zg
printVarSummary(var)
end if

if(fileinfor) then 
numb   = 0
nmodel = 0
print("file number " + nfile)
print(" ")
do nf = 0, nfile-1, 1
print(nf + " " + filepath(nf))
if((nf.le.nfile-2).and.(modelname(nf).eq.modelname(nf+1))) then 
    numb = numb + 1
else
    numb = numb + 1
    nmodel = nmodel+1
    print("=============================================================")
    print( nmodel  + " model name : " + modelname(nf) + "  filenumber : " + numb)
    f = addfile(filepath(nf),"r")
    vnames := getfilevarnames(f)           ; all variable names
    nvars  := dimsizes(vnames)
    do nv = 0,nvars-1
        vtype  = getfilevartypes(f,vnames(nv))
        vatts  := getfilevaratts(f,vnames(nv))     
        dnames := getfilevardims(f,vnames(nv))  
        dsizes := getfilevardimsizes(f,vnames(nv)) 
        if(dimsizes(dnames).ge.3) then
            print("Variable '" + vnames(nv) + "' ;type : " + vtype)
            if(.not.any(ismissing(dnames))) then
            print("  dimension names : " + str_join(dnames,","))
            else
            print("  dimension names : ")
            end if
            print("  dimension sizes : " + str_join(""+dsizes,","))
            print("  attributes      : " + str_join(vatts,","))

            var := f->$vnames(nv)$@units
            print(vnames(nv) + "@units : " + var)
            var := f->$dnames(1)$
            print(dnames(1) + "@units : " + var@units)
            print(dnames(1) + " : " + str_join((/var/),","))
        end if
    end do
    numb  = 0   
    print("")
    print("")
end if
end do
end if

end
    
