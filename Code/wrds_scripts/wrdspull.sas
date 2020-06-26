
%global WRDS;
%let WRDS = ~/WRDS/Output;

*************** COMPUSTAT **********************************;
proc sql;
        create table sp500 as
            select distinct datadate, bkvlps, dvpsxm, epsx12, prccm
                from comp.idx_mth where(GVKEYX="000003");
quit;
proc export data=sp500 
outfile="&WRDS/sp500_book_to_market.xlsx" dbms = 
xlsx replace; 
run;
*************** IBES ***********************************;
proc sql;
  create table EPS_estimates1 as select 
  ticker, oftic, statpers, measure, fiscalp, fpi, estflag, meanest, 
fpedats, actual
  from ibes.statsum_epsus
  where oftic="SPX" and fiscalp="ANN"and measure="EPS";
quit;

Proc sql;
	create table EPS_estimates2 as select 
	ticker, oftic, statpers, measure, fiscalp, fpi, estflag, meanest, fpedats, actual
	from ibes.statsum_epsint
	where oftic="SPX" and fiscalp="ANN"and measure="EPS";
Quit; 

Proc sql;
	create table EPS_estimates as
	select * from EPS_estimates1
	union
	select * from EPS_estimates2;
quit;

Proc sql;
	create table EPS_estimates as 
	select *, Count(*) as CNT from EPS_estimates
	GROUP BY actual, meanest, staters, fpedats, fyi
	Having Count(*) = 1;
Quit;

proc export data=EPS_estimates 
outfile="&WRDS/EPS_estimates.csv" dbms = csv replace;
run;
