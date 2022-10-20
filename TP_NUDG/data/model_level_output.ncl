;***************************************************************************
;Mission:
;1. read variable related with level, include lev,ilev,hybm,hyam,hyai,hybi
;2. print these variable
;                                           by Ql_Ren
;                                           2019/10/22
;******************************************************************************
begin
filein   = "/users/yangsong3/renql/F/AMIP_C5PM/exe/AMIP_C5PM.cam.h1.1980-01-01-00000.nc"
var_name = (/"lev","hyam","hybm","ilev","hyai","hybi"/)
nvar = dimsizes(var_name)

f = addfile(filein,"r")
do nv = 0, nvar-1,1
    vars := f->$var_name(nv)$
    print(var_name(nv))
    print(vars)
end do 

end

