qui set haverdir "K:/DLX/DATA/"
if "`c(os)'" == "Unix" {
   global path "/san/RDS/Work/cmf/b1crs02/Projects/etp_personal/haver_download/"
   global filesep "/"
}
else if "`c(os)'" == "Windows" {
   global path = "\\rb\b1\NYRESAN\RDS\Work\cmf\b1crs02\Projects\etp_personal\haver_download"
   global filesep "\"
}

cd $path 

import haver FTBTR1ME@USECON, clear
gen time1 = dofm(time)
format time1 %td
drop time 
rename time1 date

outsheet * using "tb1m.csv", comma replace 