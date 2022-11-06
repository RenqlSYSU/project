#!/bin/sh

LID=`date +%y%m%d`

echo "Mission1.1: add project file..."
cd ~/uor_track

find . -name "*.ncl" | xargs git add
find . -name "*.sh" | xargs git add
find . -name "*.f90" | xargs git add
find . -name "*.F90" | xargs git add
find . -name "*.F" | xargs git add
#find . -name "*.vbs" | xargs git add
find . -name "*.py" | xargs git add
find . -name "*.txt" | xargs git add
find . -name "*.csh" | xargs git add
find . -name "*.md" | xargs git add
find . -name "*.m" | xargs git add
find . -name "*.in" | xargs git add

#git add */script/*
#git add */SourceMods*

echo "Mission1.2: remove project's deleted file..."
git status | grep 'deleted' | cut -d ':' -f 2 | xargs git rm --cached

git commit -m "${LID}"
git push --force origin master

