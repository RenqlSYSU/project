#!/bin/sh

#HPC4 East Campus
# rsync -vzrtuopg  /home/ys17-19/renql/project/ yangsong3@hpc4.sysu.edu.cn:/users/yangsong3/renql/project/  >> ~/cron/sys_report/$LID-project-report
 rsync -vzrtuopg -progress --delete  yangsong3@hpc4.sysu.edu.cn:/users/yangsong3/renql/project/ /home/ys17-19/renql/project/  >> ~/cron/sys_report/$LID-project-report
