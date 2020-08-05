#!/user/bin/ebv python
from ecmwfapi import ECMWFDataServer
server = ECMWFDataServer()

var_id = [131,132,135,130,129]
var_na = ["uwnd","vwnd","omega","air","hgt"]

# for year in range(2017,2018):
year = 2018

for ivar in range(0,4):
    server.retrieve({
        "class"   : "ei",
        "dataset" : "interim",
        "date"    : str(year)+"-01-01/to/"+str(year)+"-12-31",
        "expver"  : "1",
        "grid"    : "1.5/1.5",
        "levelist": "1/2/3/5/7/10/20/30/50/70/100/125/150/175/200/225/250/300/350/400/450/500/550/600/650/700/750/775/800/825/850/875/900/925/950/975/1000",
        "levtype" : "pl",
        "param"   : str(var_id[ivar])+".128",
        "step"    : "0",
        "stream"  : "oper",
        "time"    : "00:00:00/06:00:00/12:00:00/18:00:00",
        "type"    : "an",
        "format"  : "netcdf",
        "target"  : "/home/ys17-19/data-observation/ERA-interim/pressure/"+var_na[ivar]+"/"+var_na[ivar]+".interim."+str(year)+".nc",
    })
