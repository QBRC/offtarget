R --vanilla --slave --args $1 $2 $3 $4 $5 $6 $7 < /home/galaxy/galaxy-dist/tools/offtarget/offtarget_v7.R > dump

if [ -f ${__tool_data_path__}CorrectedZScore.csv ]
 then
   mv dump log.txt
   zip temp ${__tool_data_path__}*.*
   mv ${__tool_data_path__}temp.zip $8 
 else 
   cat dump >&2
fi

