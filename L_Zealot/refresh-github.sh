#!/bin/sh

LID=`date +%y%m%d`

#--------------------------
# Mission1: project
#--------------------------
echo "Mission1: project..."
cd ~/L_Zealot/project

find . -name "*.ncl" | xargs git add
find . -name "*.sh" | xargs git add
find . -name "*.f90" | xargs git add
find . -name "*.F" | xargs git add
find . -name "*.vbs" | xargs git add
find . -name "*.py" | xargs git add
find . -name "*.txt" | xargs git add
find . -name "*.csh" | xargs git add
find . -name "*.md" | xargs git add

git add */script/*
git add */SourceMods*

git commit -m "${LID}"
git push --force origin master

#--------------------------
# Mission2: paperhub
#--------------------------

echo "Mission2: paperhub..."
cd ~/L_Zealot/workspace/paperhub
find . -name "*.php" | xargs git add
find . -name "*.dat" | xargs git add
find . -name "*.md" | xargs git add
find . -name "*.py" | xargs git add
find . -name "*.sh" | xargs git add

git commit -m "${LID}"
git push origin master
