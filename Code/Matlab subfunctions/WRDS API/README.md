[![DOI](https://zenodo.org/badge/doi/10.5281/zenodo.19056.svg)](http://dx.doi.org/10.5281/zenodo.19056)

![interaction](https://raw.githubusercontent.com/okomarov/wrds/master/logos.png)

### Description
High level Matlab API that interacts with the Wharton Reasearch Data Services (WRDS) Unix server and its SAS data sets through SSH2.
    
### Requirements
* An account with WRDS of the type that admits SSH connections (PhD or above). 
  See WRDS's [account types](http://wrds-web.wharton.upenn.edu/wrds/support/Additional%20Support/Account%20Types.cfm) for details.
* [Java enabled](http://www.mathworks.co.uk/help/matlab/ref/usejava.html)

### Syntax

Connect to WRDS Unix servers

    w = wrds('username', 'password')

Execute commands

    w.cmd('cmdstring')
    
Check what specific datasets a (database) library has 

    w.getDatasetNames('libref')
    
Check info for a dataset you want to download

    w.getDatasetInfo('libref.datasetname')
    
and, if you get an empty result, then your institution might not have acces to that particular database. Check on  https://wrds-web.wharton.upenn.edu/wrds/mywrds/subscriptions.cfm? which databases you are subscribed on!
    
Download the dataset to a zipped .csv

    w.getDataset('libref.datasetname')

Switch off verbosity

    w.isVerbose = false;

###Examples:

Print `Hello World` in the Unix shell and print result into Matlab's cmd window (Verbose)

    w = wrds('olegkoma','forgiveMeIfIDontTellYou');
    w.cmd('echo "Hello World!"')


Get info about the [CRSPA/MSI dataset](http://wrds-web.wharton.upenn.edu/wrds/tools/variable.cfm?library_id=137&file_id=67079)

    w.getDatasetInfo('CRSPA.MSI')
    ans = 
    ...
                      Alphabetic List of Variables and Attributes                   
                                                                                
    # Variable Type Len Format    Informat Label                                  
                                                                                
    1 DATE     Num    4 YYMMDDN8. YYMMDD8. Date of Observation                    
    4 ewretd   Num    8 10.6      E13.6    Equal-Weighted Return-incl. dividends  
    5 ewretx   Num    8 10.6      E13.6    Equal-Weighted Return-excl. dividends  
    7 spindx   Num    8 8.2       8.2      Level of the S&P 500 Index             
    6 sprtrn   Num    8 10.6      E13.6    Return on the S&P 500 Index            
    9 totcnt   Num    8 5.        5.       Total Market Count                     
    8 totval   Num    8 14.2      E13.6    Total Market Value                     
    11 usdcnt   Num    8 5.        5.       Count of Securities Used               
    10 usdval   Num    8 14.2      E13.6    Market Value of Securities Used        
    2 vwretd   Num    8 10.6      E13.6    Value-Weighted Return-incl. dividends  
    3 vwretx   Num    8 10.6      E13.6    Value-Weighted Return-excl. dividends  
    ...
    
Download it    

    w.getDataset('CRSPA.MSI')
