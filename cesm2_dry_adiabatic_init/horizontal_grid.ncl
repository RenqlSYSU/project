    ;===================================================================
    ;  Define horizontal grid
    ;===================================================================

    ;Modified to include SE options, 17th Sept 2017 - Isla Simpson   

       if (res .eq. "8x16") then
           nlat = 8 
           nlon = 16 
           if (se) then 
            print( (/"8x16 is not a valid resolution for the SE dycore"/) )
            exit
           end if
       end if

       if (res .eq. "64x128") then
           nlat = 64 
           nlon = 128 
           if (se) then 
            print( (/"64x128 is not a valid resolution for the SE dycore"/) )
            exit
           end if
       end if

       if (res .eq. "128x256") then
           nlat = 128
           nlon = 256 
          if (se) then
            print( (/"128x256 is not a valid resolution for the SE dycore"/) )
            exit
           end if
       end if

       if (res .eq. "256x512") then
           nlat = 256 
           nlon = 512 
           if (se) then
            print( (/"256x512 is not a valid resolution for the SE dycore"/) )
            exit
           end if
       end if

       if (res .eq. "512x1024") then
           nlat = 512
           nlon = 1024
           if (se) then
            print( (/"512x1024 is not a valid resolution for the SE dycore"/) )
            exit
           end if
       end if

       if (res .eq. "ne16") then 
          in=addfile("./grids/ne16_lonlat.nc","r")
          lon=in->lon
          lat=in->lat
          nlat=dimsizes(lat)
          ncol=dimsizes(lat)

          if (.not.se) then 
            print( (/"ne16 is only a valid resolution for the SE dycore"/) )
            exit
          end if
       end if

       if (res .eq. "ne30") then 
          in=addfile("./grids/ne30_lonlat.nc","r")
          lon=in->lon
          lat=in->lat
          nlat=dimsizes(lat)
          ncol=dimsizes(lat)
          if (.not.se) then 
            print( (/"ne30 is only a valid resolution for the SE dycore"/) ) 
            exit
          end if
       end if

       if (res .eq. "ne120") then 
          in=addfile("./grids/ne120_lonlat.nc","r")
          lon=in->lon
          lat=in->lat
          nlat=dimsizes(lat)
          ncol=dimsizes(lat)
          if (.not.se) then 
            print( (/"ne120 is only a valid resolution for the SE dycore"/) ) 
            exit
          end if
       end if



       if (.not.se) then 

         gaus_info = gaus(nlat/2) ;fspan(-90.0d0+90.0d0/nlat,90.0d0-90.d0/nlat,(/nlat/))

         lat = gaus_info(:,0)
         lat!0="lat"
         lat@long_name = "latitude"
         lat@units = "degrees_north"

         gw = gaus_info(:,1)
         gw!0 = "lat"
         gw@long_name = "gauss weights"
  
         lon = fspan(0.d0,360.-360./nlon,(/nlon/))
         lon!0 = "lon"
         lon@long_name = "longitude"
         lon@units = "degrees_east"

       end if 

          
