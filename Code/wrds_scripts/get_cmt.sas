
%global WRDS;
%let WRDS = ~/WRDS/Output;

PROC SQL;
CREATE TABLE CMT AS
SELECT KYTREASNOX as cmtcode, MCALDT as date, TMRETADJ as value FROM CRSP.TFZ_MTH_FT Order By KYTREASNOX, date; 
Run;
Quit;

Proc export data = cmt
	outfile = "&WRDS/cmt_data.csv" dbms = dlm replace;
	delimiter = ',';
Run;
