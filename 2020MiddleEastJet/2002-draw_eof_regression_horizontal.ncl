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
draw_var = (/"uwnd"/)  ;this is the variable used to draw, can used as center string
var_name = (/"season_ave"/)  ;use to read variables
unit     = (/"m/s"/)  ;this can be used as right string
nv = 0
lev  = (/850,500,200,100/)
nl = 2

area = (/"ME","EA","NA"/);,"GB"
lats = (/ 15 , 20 , 25 /);, 10  ;-90 ,-90 ,-90 ,-90 
latn = (/ 40 , 45 , 50 /);, 60  ;90 , 90 , 90 , 90 
lonl = (/ 0  ,105 , 250/);, 0   ;30 ;20
lonr = (/ 70 ,175 , 320/);, 360 ;150;220
narea = dimsizes(area)
neval = 4
years = ispan(1979,2015,1)

path   = "/home/ys17-19/renql/project/2020MiddleEastJet"
filename   = path + "/data/eof_u_6kinds-" + area + ".nc"
fig_name   = path + "/fig/eof_u_ts"
fig_out    = "pdf" ;"eps";"x11"
pre_case   = (/"(a","(b","(c", "(d","(e","(f", "(g","(h","(i"/)
plot   = new((/3/),graphic)

;define the thickness and height of the title,tick and label
   tick_font      = 0.01
   subtitle_font  = 0.013
   title_font     = 0.03
   label_font     = 0.015
   tick_thick     = 5.0
   subtitle_thick = 5.0  ;not have this attribute
   title_thick    = 5.0
   label_thick    = 1.5
   mp_thick       = 4.0
   font_index     = 22

f = addfile(filename(na),"r")
var = f->eof_ts

;===============================================
;draw the figure
;=============================================
wks = gsn_open_wks(fig_out,fig_name)

