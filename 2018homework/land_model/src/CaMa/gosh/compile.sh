#!/bin/sh 

clean=$1       # 1th argument in script : if == "yes" then clean before compiling 

CAMADIR=`pwd`/..
BASE=$CAMADIR  # code location

libs="mod lib src map out"
for lib in $libs
do 
  cd $BASE/$lib
# if [ $clean = "yes" ];then
    make clean
# fi
  echo "*********** $lib **********"
  make all
done 

echo "CaMa compilation has been completed"
