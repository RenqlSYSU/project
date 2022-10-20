#!/bin/sh

# Update paperhub
LID=`date +%y%m%d`
rsync -vzrtulopg -progress -e 'ssh -p 66' eesael@222.200.180.66:/home/eesael/www/html/L_Zealot/paperhub   /home/yangsong3/L_Zealot/workspace/  >& ~/cron/sys_report/$LID-paperhub-report
rsync -vzrtulopg -progress -e 'ssh -p 66' /home/yangsong3/L_Zealot/workspace/paperhub eesael@222.200.180.66:/home/eesael/www/html/L_Zealot/ >> ~/cron/sys_report/$LID-paperhub-report

# ETL process, refine gold data
cd /home/yangsong3/L_Zealot/workspace/economy/
sh refine-gold-data.sh
