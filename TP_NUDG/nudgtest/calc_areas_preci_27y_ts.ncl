;***************************************************************************
;filename:      calc_areas_preci_27y_ts.ncl
;output_file:   Areas_preci_27y_ts.nc 
;read the 27 years JJA preci(PRECC+PRECL) in two case (CTRL, NUDG) and then 
;calculate the season and area average in order to test the right of the 
;NUDG experiment.
;                                            by Ql_Ren
;                                           2017/12/21
;******************************************************************************
begin
;the path of the file
filein   = (/"/users/yangsong3/Model_Team/F/AMIP-CTRL/pro/AMIP_C5PM.cam.h1.YEAR.",\ ;+year+".daily.nc"
             "/users/yangsong3/Model_Team/F/TP-NUDG/pro/AMIP_C5PM_TP_NUDG.cam.h1.YEAR."/) ;+year+".daily.nc"
fileout  = "/users/yangsong3/Model_Team/F/nudgtest/27years_time_series.nc"
var_name = (/"PRECC","PRECL"/)
case     = (/"AMIP_CTRL","TP_NUDG"/)
years    = ispan(1979,2005,1)   ;the calculate time

;area of index, India,Tibet,Philippines,Maritime_continent
area  = (/"India","Tibet","Philippines","Maritime_continent",\
          "North China","Central East China","South China"/)
lats = (/20 ,25  ,5   ,-10 ,34  ,28  ,20 /)
latn = (/35 ,35  ,20  ,7.5 ,43  ,34  ,28 /)
lonl = (/60 ,90  ,120 ,85  ,107 ,107 ,107/)
lonr = (/80 ,103 ,140 ,120 ,122 ,122 ,122/)

;==========================================================
;create an array(nvar,ncase,nyear) to store the data
;============================================================
f    = addfile(filein(0)+years(0)+".daily."+var_name(0)+".nc","r")
time = cd_calendar(f->time,0)
select_time = ind(time(:,1).ge.6.and.time(:,1).le.8)
ntime = dimsizes(select_time)
nyear = dimsizes(years)
ncase = dimsizes(case)
narea = dimsizes(area)
index = new((/ncase,narea,nyear/),float)

;---------------------------------------------------------------
;read the data, then calculate the area and season average
;--------------------------------------------------------------
do nf = 0,1,1 
do na = 0,narea-1,1
    do t = 0,nyear-1,1
        index(nf,na,t) = 0
        do v = 0, dimsizes(var_name)-1,1  ;the aim of this loop is to add the PRECC and PRECL 
            f  = addfile(filein(nf)+years(t)+".daily."+var_name(v)+".nc","r")
            vars := f->$var_name(v)$(select_time,{lats(na):latn(na)},{lonl(na):lonr(na)})
            vars := dim_avg_n_Wrap(vars,0)  ;calculate the JJA average
            index(nf,na,t) = index(nf,na,t) + wgt_areaave_Wrap(vars,1.0,1.0,0)*3600*24*1000 ;calculate the area average and get the index
        end do ;the loop of precc and precl
    end do  ;the loop of the years
end do   ;the loop of the areas
end do   ;the loop of the case
printVarSummary(index)

;====================================================================
;save the data
;=======================================================================
system("rm -f " + fileout)
ncdf = addfile(fileout,"c")
setfileoption(ncdf,"DefineMode",True)

;Set All field attribution
r
fileAtt =  True
fileAtt@creation_date = systemfunc("date")
fileAtt@discription   = "The index is JJA preci. The two case are CTRL,NUDG" + \
                    "the seven areas are India,Tibet,Philippines,Maritime_continent," +\
                    "North China,Central East China,South China"
fileattdef(ncdf,fileAtt)  

;Define all Coordinate
dimNames = (/"case","area""year"/)
dimSizes = (/ncase,narea,nyear/)
dimUnlim = (/False, False,False/)
filedimdef(ncdf, dimNames, dimSizes, dimUnlim)

;Define var_names, type(string) and dims_names(string)
filevardef(ncdf,"case","integer",(/"case"/))
filevardef(ncdf,"area","integer",(/"area"/))
filevardef(ncdf,"year","integer",(/"year"/))
filevardef(ncdf,"index",typeof(index),("case","area","year"/))

;Define the attribute of the variables
index@long_name = "the time series of two case and seven areas of 27 years"
index@units     = "mm/day"
filevarattdef(ncdf,"index",index)

ncdf->case   = (/1,2/)  ;the two cases are (/"CTRL","NUDG"/)
ncdf->area   = (/1,2,3,4,5,6,7/)  ;the two cases are (/"CTRL","NUDG"/)
ncdf->year   = (/years/)
ncdf->index  = (/index/)
end

