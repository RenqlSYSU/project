;******************************************************************************
;filename:     read_nc_write_dat.ncl 
;input_file:   ERSST/ersst.v5.*.nc
;because the nc file's time@calendar is 360_day, so the grads can not read it
;this purpose of this file is to read the all nc files and then convert it to 
;one dat file. Next write the ctl file manually.
;                                            by Ql_Ren
;                                           2018/04/16
;******************************************************************************
begin
;------the path of the file
path     = "/home/ys17-19/renql/project"
filein   = path + "/help_others/baizih/ERSST/ersst*"
fileout  = path + "/help_others/baizih/ERSST_dat/ersst.v5.195001-201612.dat"
var_name = (/"sst","ssta"/)

;-------read the data-----------
files1  = systemfunc("ls " + filein)
f       = addfiles(files1,"r")
ListSetType(f,"cat")

var1    = f[:]->$var_name(0)$
var2    = f[:]->$var_name(1)$
ntime   = dimsizes(var1&time)
print(ntime)
printVarSummary(var1)

system("rm -f " + fileout)
do nt=0,ntime-1
    fbindirwrite(fileout,var1(nt,0,:,:))
    fbindirwrite(fileout,var2(nt,0,:,:))
end do

end
