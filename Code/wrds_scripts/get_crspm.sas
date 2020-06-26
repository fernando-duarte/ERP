
%global WRDS;
%let WRDS = ~/WRDS/Output;

PROC SQL;
	CREATE TABLE CCM_Link as 
	select distinct a.permno, a.permco, gvkey, liid as iid, date, prc, ret, shrout, altprc, b.shrcd, b.exchcd
	from crsp.msf as a,
	crsp.msenames 
	( 
		where = (shrcd in (10,11))
	) as b,
	crsp.Ccmxpf_linktable
	(
		where = (
		linktype in ('LU' 'LC')
		and LINKPRIM in ('P' 'C')
		and USEDFLAG = 1 )
	) as c
	where a.permno = b.permno = c.lpermno
	and NAMEDT <= a.date <= NAMEENDT;
Quit;

Proc sql; 
	create table crsp_m as
	select date, permno, permco, shrcd, exchcd, prc, ret, shrout, altprc
	from work.ccm_link;
Quit;

Proc export data = crsp_m
	outfile = "&WRDS/crsp_m.txt" dbms = dlm replace;
	delimiter = ' ';
Run;
