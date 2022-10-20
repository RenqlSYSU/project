#!/bin/sh

LID=`date +%y%m%d`

#Workstation East Campus
rsync -vzrtuopg -e 'ssh -p 80' /home/yangsong3/L_Zealot/project/ workstation@222.200.180.66:/home/workstation/L_Zealot/project/ >& ~/cron/sys_report/$LID-project-report
rsync -vzrtuopg -e 'ssh -p 80' workstation@222.200.180.66:/home/workstation/L_Zealot/project/ /home/yangsong3/L_Zealot/project/ >> ~/cron/sys_report/$LID-project-report
echo "Workstation done."

#HPC4 East Campus
rsync -vzrtuopg  /home/yangsong3/L_Zealot/project/ yangsong3@hpc4.sysu.edu.cn:/users/yangsong3/L_Zealot/project/  >> ~/cron/sys_report/$LID-project-report 
rsync -vzrtuopg  yangsong3@hpc4.sysu.edu.cn:/users/yangsong3/L_Zealot/project/ /home/yangsong3/L_Zealot/project/  >> ~/cron/sys_report/$LID-project-report
echo "HPC4 done."

#Mac Pro South Campus
rsync -vzrtuopg  /home/yangsong3/L_Zealot/project/ yangsong3@202.116.70.22:/Users/yangsong3/L_Zealot/project/  >> ~/cron/sys_report/$LID-project-report 
rsync -vzrtuopg  yangsong3@202.116.70.22:/Users/yangsong3/L_Zealot/project/ /home/yangsong3/L_Zealot/project/  >> ~/cron/sys_report/$LID-project-report
echo "Mac Pro done."

#Mac Mini Cal
#rsync -vzrtuopg  /home/yangsong3/L_Zealot/project/ zhenningli@169.229.2.23:/Users/zhenningli/project/  >> ~/cron/sys_report/$LID-project-report 
#rsync -vzrtuopg  zhenningli@169.229.2.23:/Users/zhenningli/project/ /home/yangsong3/L_Zealot/project/  >> ~/cron/sys_report/$LID-project-report
echo "Mac Mini done."

sh /home/yangsong3/cron/refresh-github.sh >> ~/cron/sys_report/$LID-project-report
