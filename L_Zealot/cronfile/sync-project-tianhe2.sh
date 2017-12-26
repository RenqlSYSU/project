#!/bin/sh

LID=`date +%y%m%d`

rsync -vzrtuopg -progress -e 'ssh -i /home/yangsong3/sysu_hjkx_ys.id' /home/yangsong3/L_Zealot/project/ sysu_hjkx_ys@172.16.22.11:/HOME/sysu_hjkx_ys/WORKSPACE/L_Zealot/project/ #>& ~/cron/sys_report/$LID-project-report-th2

rsync -vzrtuopg -progress -e 'ssh -i /home/yangsong3/sysu_hjkx_ys.id' sysu_hjkx_ys@172.16.22.11:/HOME/sysu_hjkx_ys/WORKSPACE/L_Zealot/project/ /home/yangsong3/L_Zealot/project/ #>& ~/cron/sys_report/$LID-project-report-th2

