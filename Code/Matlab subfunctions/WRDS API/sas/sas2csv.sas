* Pipe into .zip;
filename %s pipe "zip -r > %s";

* Taken from https://communities.sas.com/message/185633#185633;
* Read dataset variable names;
proc sql noprint;
 select '"'||trim(name)||'"'
 into :names
 separated by "','"
 from dictionary.columns
 where libname eq "%s" and
 memname eq "%s"
 ;
quit;
* Write data;
data _null_;
 set %s;
 file %s mod dsd dlm=',' lrecl=1000000;
 if _n_ eq 1 then put &names.;
 put (_all_) (+0); 
run;