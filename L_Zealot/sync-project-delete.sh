#!/bin/sh

date

echo "======================="
echo "Neptune --> workstation"
echo "======================="
rsync -vzrtuopg -progress --delete -e 'ssh -p 80' /home/yangsong3/L_Zealot/project/ workstation@222.200.180.66:/home/workstation/L_Zealot/project/ 
echo "======================="
echo "workstation --> Neptune"
echo "======================="
rsync -vzrtuopg -progress --delete -e 'ssh -p 80' workstation@222.200.180.66:/home/workstation/L_Zealot/project/ /home/yangsong3/L_Zealot/project/ 

echo "======================="
echo "   Neptune --> HPC4"
echo "======================="
rsync -vzrtuopg -progress --delete -e 'ssh -p 22' /home/yangsong3/L_Zealot/project/ yangsong3@hpc4.sysu.edu.cn:/users/yangsong3/L_Zealot/project/  
echo "======================="
echo "   HPC4 --> Neptune"
echo "======================="
rsync -vzrtuopg -progress --delete -e 'ssh -p 22' yangsong3@hpc4.sysu.edu.cn:/users/yangsong3/L_Zealot/project/ /home/yangsong3/L_Zealot/project/  


echo "======================="
echo "   Neptune --> Tianhe2"
echo "======================="
rsync -vzrtuopg -progress --delete -e 'ssh -i /home/yangsong3/sysu_hjkx_ys.id' /home/yangsong3/L_Zealot/project/ sysu_hjkx_ys@172.16.22.11:/HOME/sysu_hjkx_ys/WORKSPACE/L_Zealot/project/ #>& ~/cron/sys_report/$LID-project-report-th2
echo "======================="
echo "   Tianhe2 --> Neptune"
echo "======================="
rsync -vzrtuopg -progress --delete -e 'ssh -i /home/yangsong3/sysu_hjkx_ys.id' sysu_hjkx_ys@172.16.22.11:/HOME/sysu_hjkx_ys/WORKSPACE/L_Zealot/project/ /home/yangsong3/L_Zealot/project/ #>& ~/cron/sys_report/$LID-project-report-th2

