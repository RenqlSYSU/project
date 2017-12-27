#!/bin/sh

LID=`date +%y%m%d`
# Lab Server
#rsync -vzrtuopg  -e 'ssh -p 66' eesael@222.200.180.66:/home/eesael/www/html/L_Zealot/paperhub   /home/yangsong3/L_Zealot/workspace/  >& ~/cron/sys_report/$LID-paperhub-report
#rsync -vzrtuopg  -e 'ssh -p 66' /home/yangsong3/L_Zealot/workspace/paperhub eesael@222.200.180.66:/home/eesael/www/html/L_Zealot/ >> ~/cron/sys_report/$LID-paperhub-report

# Cal Mini
#rsync -vzrtuopg zhenningli@169.229.2.23:/Users/zhenningli/workspace/paperhub    /home/yangsong3/L_Zealot/workspace/  >& ~/cron/sys_report/$LID-paperhub-report
#rsync -vzrtuopg /home/yangsong3/L_Zealot/workspace/paperhub zhenningli@169.229.2.23:/Users/zhenningli/workspace/ >> ~/cron/sys_report/$LID-paperhub-report
