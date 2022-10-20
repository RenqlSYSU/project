;*************************************************
; table_4.ncl
;
; Concepts illustrated:
;   - Drawing a table using gsn_table
;   - Filling table cells with a given color
;   - Adding tickmark labels to a table using gsn_csm_blank_plot
;   - Attaching a labelbar to a plot
;   - Explicitly setting tickmarks and labels on the bottom X axis
;   - Explicitly setting tickmarks and labels on the left Y axis
;   - Turning off the top and right tickmarks
;   - Reading an ASCII file with several columns of data
;   - Reading a specific field from an ASCII file
;   - Using "getvalues" to retrieve the size of a plot
;
;*************************************************
;
; These files are loaded by default in NCL V6.2.0 and newer
; load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
; load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
; load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"


;----------------------------------------------------------------------
; Procedure to attach a horizontal labelbar to the bottom of plot.
;----------------------------------------------------------------------
undef("add_labelbar")
procedure add_labelbar(wks,plot,colors,labels)
local vph, vpw, nboxes, lbres, lbid, amres, annoid
begin
  getvalues plot                      ; Get plot size for use in
    "vpHeightF" : vph                 ; creating labelbar.
    "vpWidthF"  : vpw
  end getvalues

  nboxes = dimsizes(colors)

  lbres                    = True          ; labelbar only resources
  lbres@lbAutoManage       = True          ; Necessary to control sizes
  lbres@vpWidthF           = 0.15 * vpw    ; labelbar width
  lbres@vpHeightF          = 0.95 * vph    ; labelbar height
  lbres@lbFillColors       = colors        ; labelbar colors
  lbres@lbMonoFillPattern  = True          ; Solid fill pattern
  lbres@lbLabelFontHeightF = 0.05          ; font height. default is small
  lbres@lbLabelAlignment   = "BoxCenters"  ; center of box
  lbres@lbOrientation      = "Vertical"
  lbres@lbPerimOn          = False

  lbid = gsn_create_labelbar(wks,nboxes,labels,lbres)
   
;
; Now, create some annotation resources indicating how we want to
; attach the labelbar to the plot. Here, we are using the top right
; corner of the labelbar as the point which we are going to position
; it, and then we use amParallelPosF and amOrthogonalPosF to indicate
; where we want to place it.
;
; amParallelPosF/amOrthogonalPosF                                       ;
;    0.0/ 0.0  -  annotation in dead center of plot
;    0.5/ 0.5  -  annotation at bottom right of plot
;    0.5/-0.5  -  annotation at top right of plot
;   -0.5/-0.5  -  annotation at top left of plot
;   -0.5/ 0.5  -  annotation at bottom left of plot
;                                                                 
  amres                  = True
  amres@amJust           = "CenterLeft"
  amres@amParallelPosF   =  0.52
  amres@amOrthogonalPosF =  0.0
  plot@annoid = gsn_add_annotation(plot,lbid,amres)
end

;----------------------------------------------------------------------
; Main code
;----------------------------------------------------------------------
begin
;
; Read some dummy data with labels. We know the first row contains labels
; and the first column contains labels. We don't know how many rows or 
; columns there are.
;
  asc_file = "table4.txt"
  data     = asciiread(asc_file,-1,"string")

; Count rows and columns of actual data.
  ncols    = numAsciiCol(asc_file)
  nrows    = numAsciiRow(asc_file)-1   ; Don't include first row which doesn't have numerical data

; Read the values
  values      = asciiread(asc_file,(/nrows,ncols/),"float")
  text_values = "" + values      ; Create a string version for the table

; Read the labels
  ncol_labels = str_split(data(0)," ")
  nrow_labels = str_get_field(data(1:),1," ")


  ncr    = dimsizes(values)
  col_width  = 1./ncols
  row_width  = 1./nrows
  col_width2 = col_width/2.
  row_width2 = row_width/2.

;
; Define colors to use for the table. Be sure to include white and black for background
; and foreground.
;
  colors = (/"White","Black","Navy","MediumBlue","MediumTurquoise","SeaGreen1",\
             "Green","Yellow","Orange","Firebrick1","Red","Brown"/)

  labels = new(dimsizes(colors)-1,string)    ; for labelbar

; Assign a color to a particular range of data.
  cols_for_vals = new(ncr,"string")
  cols_for_txt = new(ncr,"string")

  do i=2,11
    rng_beg = (i-2)*0.1
    rng_end = (i-1)*0.1
    cols_for_vals = where(rng_beg.le.values.and.values.lt.rng_end,colors(i),\
                          cols_for_vals)
;---Create formatted labels for the labelbar
    labels(i-2) = sprintf("%3.1f", rng_beg) + " - " + \
                  sprintf("%3.1f",rng_end)
  end do

  cols_for_txt = where(cols_for_vals.eq."Navy".or.cols_for_vals.eq."MediumBlue","white","black")

; Start the graphics.

  wks = gsn_open_wks("png","table")  ; send graphics to PNG file
  gsn_define_colormap(wks,colors)

;
; Create a blank plot so we can get some tickmarks.
;
; Do the blank plot before the table, so we can maximize
; the blank plot. Later we'll make sure the table goes
; in the same location as the blank plot.
;
  resb                         = True
  resb@gsnDraw                 = False
  resb@gsnFrame                = False

  resb@vpYF                    = 0.92
  resb@vpXF                    = 0.10
  resb@vpHeightF               = 0.85
  resb@vpWidthF                = 0.75

; Explicitly label X axis. The blank plot goes from 0 to 1, by default.
  resb@tmXBMode                = "Explicit"
  resb@tmXBValues              = fspan(col_width2,1.-col_width2,ncols)
  resb@tmXBLabels              = ncol_labels
  resb@tmXBLabelFontHeightF    = 0.015            ; Make font a little smaller

; Explicitly label Y axis.
  resb@tmYLMode                = "Explicit"
  resb@tmYLValues              = fspan(row_width2,1.-row_width2,nrows)
;
; The Y labels need to be reversed, because the label in the first
; row in the ASCII file is associated with the top of the plot.
;
  resb@tmYLLabels              = nrow_labels(::-1)
  resb@tmYLLabelFontHeightF    = 0.015            ; Make font a little smaller

; Rotate the X axis labels 90 degrees.
  resb@tmXBLabelAngleF         = 90.
  resb@tmXBLabelJust           = "CenterRight"

; Turn off right and top tickmarks
  resb@tmYROn                  = False
  resb@tmXTOn                  = False

; Main title
  resb@tiMainString            = "Table of dummy data"

  blank = gsn_csm_blank_plot(wks,resb)
  add_labelbar(wks,blank,colors(2:),labels)   ; Attach labelbar
;
; Get position and size of the blank plot so we can
; be sure to draw the table in same location.
;
  getvalues blank
    "vpXF"      : vpx
    "vpYF"      : vpy
    "vpWidthF"  : vpw
    "vpHeightF" : vph
  end getvalues

; Start and end coordinates for table.
  x = (/vpx,vpx+vpw/)
  y = (/vpy-vph,vpy/)

; Set up resources for table
  rest               = True

  rest@gsLineColor   = -1                ; No border lines
  rest@gsFillColor   = cols_for_vals     ; Each cell of table
  rest@txFontColor   = cols_for_txt

  rest@txFontHeightF = 0.015             ; Size of centered text

  gsn_table(wks,ncr,x,y,text_values,rest)    ; Draw table
  draw(blank)                                ; Draw tickmarks and labelbar

; Advance the frame.
  frame(wks)

end
