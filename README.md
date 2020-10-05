# The Equity Risk Premia (ERP): A Review of Models 

## 1	Introduction
The equity risk premium — the expected return on stocks in excess of the risk-free rate — is a fundamental quantity in all of asset pricing, both for theoretical and practical reasons. It is a key measure of aggregate risk-aversion and an important determinant of the cost of capital for corporations, the savings decisions of individuals, and budgeting plans for governments. Recently, the equity risk premium (ERP) has also moved to the forefront as a leading indicator of the evolution of the economy, a potential explanation for jobless recoveries, and a gauge of financial stability. 

## 2	Software Dependencies
- MATLAB 2020a with Statistics and Machine Learning toolboxes
- STATA version 16 or greater
- HAVER privilege via DLX on the Fed SAN
- Anaconda for the latest Python version
- WRDS account, available [here](https://wrds-www.wharton.upenn.edu/)
- FRED account for an API key, found [here](https://research.stlouisfed.org/docs/api/api_key.html)
- System environment with 25 GB of memory

## 3	Repository Location
The repository for the Equity Risk Premium has been moved to GitHub [here](https://github.com/fernando-duarte/ERP). Clone the repository to your RAN machine and make all changes there. When naming the local Git repository cloned from the remote please follow the proper program naming routine (no spaces, no rare characters, etc.) so that the program does not crash due to invalid directory.

## 4	Code Structure
### 4.1 	Outline
The code base is mapped as a linear series of Matlab file, outlined below:

- `Code/project_scripts/download_data.m` downloads all available data sources to be used in the model creation process, these datasets include the following processes
   - Download WRDS data through SAS code, understand that this script may break, if this occurs the downloaded files must be done manually (refer to section 6)
   - Download assorted data from the internet, script begins with Python script and extend to various URLs ranging from the Federal Reserve to NYU Stern
   - Download FRED data, which include TIPS data, coproate bond yeilds,  short dated Treasury yeilds and recessions
     
- `Code/project_scripts/NY_Fed_CM_DAPM.m` this file recreates the NY Fed ERP measures as produced in the Adrian, Crump and Moench (2014) paper [here](https://papers.ssrn.com/sol3/papers.cfm?abstract_id=1837531)
	- Download data from the web to compute the Fama French measures, the Shiller data which is formated to fit the Matlab envrionment (subject to Shiller data formating)
	- Merge downloaded data from Fred tables and Shiller data into a large data table.
	- Clean data files for Fama French, WRDS cmt data, and saving the data to a final table, cleaning and merging to match the dates  
  
- `Code/project_scripts/Duke_CFO.m` load in Duke CFO data, performing data clearing (e.g. removing NaN) and selection before data exporting
  
- `Code/project_scripts/cay.m` this file recreates the CAY measure as produced by Lettau and Ludvigson (2001) paper [here](http://schwert.ssb.rochester.edu/f532/jf2001_ll.pdf)
	- Downloaded Fred data for various variable figures, (e..g current dollars, personal income receipts)  before mathcing dates to allign dataset. 
	- Compute the C, A and Y measures to run regression on 8 lags to compute the CAY measures,  before exporting the fiures to .csv file. 
	
- `Code/project_scripts/collect_data.m` this file reads all the original data and converts it into Matlab financial time series objects, these include the Damodaran data [here](http://people.stern.nyu.edu/adamodar/pdfiles/papers/ERP2012.pdf), the CFO data, bond yield data, Fama-French data, the Shiller data, Thomson Reuters data, FRED data, CompuStat  data, Cay data, and equity/debt issuance data.
	
- `Code/project_scripts/create_variables.m` this file homogenizes independent data measures, merging them all into a larger dataset to be used in constructing the ERP.
	- Load time series from all data providers (e.g. CFO, Fama-French) compiling them into one timetable of a common time vector along some axis. 
	- Generate Campbell-Thompson variables, a more detailied breakdown of the variable’s calculation can be found [here](https://dash.harvard.edu/bitstream/handle/1/2622619/Campbell_Predicting.pdf?sequence=2&isAllowed=y) 
  
- `Code/project_scripts/create_ERP_measures.m` this file computes several measures used for the ERP, across various horizons as a function of years (~0.08y, 0.25y, 0.5y, 1y, 2y, 3y, 4y, 5y). The code follows by averaging across various model types for computing the equity risk premium. These models are as follows: **_Historical mean models_**, **_Dividend discount mean models_**, **_Cross-sectional regression models_**, **_Time-series regression models_**, **_Survey_** (e.g. Duke survey of CFOs). Export the ERP measures with recession dummy variables added
 
- `Code/project_scripts/erp_principal_components.m` this file runs a principal component analysis on each of the ERP model estimates (e.g. PCA of DDM), before being written to a .csv file with the addition of a recession dummy.

- `Code/project_scripts/create_term_structure.m` this file constructs the term structure of the ERP, the expected excess returns over different investment horizons as outlined within the ERP paper. This includes the counter-factual measurements, which leaves expected stock returns unmodified but adjusts the risk free rates from their actual values to the average nominal bond yields over the period 1960 – 2013. Further details are highlighted within the paper.  
  
- `Code/project_scripts/test_predictability.m` this file runs regressions to test the predictability of the ERP models, as well as the principal components.
  - Begin by testing the predictability of the ERP models, across the various time horizons computed, and selecting the best out-of-sample R-squared
  - Extend the regression process outlined above on the principal components, selecting the best by observing the out-of-sample R-squared

- `Code/project_scripts/vintages.m` finally, we compare past vintages with our newly created figure, choosing to either extend or leave the vintage matrix unchanged before the values are exported to a .csv file for review

### 4.2 `/Code`
  This folder contains all the code to run the project. The code is sorted into a variety of folders. All the scripts that are called by the main file run_all.m are stored in Code\project_scripts. Functions used to compute model measures are stored in three separate folders, with separate sub-folder for additional language scripts.

### 4.3  `/Input`
  This folder contains all the raw data that is needed to run the ERP project, these include resources from FRED, WRDS as well as various online resources. 

### 4.4  `/Temp`
  This folder contains data at various stages of the project. Each folder maps to a script in the code\project_scripts folder. Each folder contains an input folder and an output folder. The input is the data that enters that script and the output is the data produced by the script.

### 4.5  `/Output`
  This folder contains all the output. The most important files are:

- `Original_data.csv` contains information on the input data, its sources, and the dates it spans. Use this to identify which datasets are updating and which are not.
- `Constructed_variables.csv` contains information on the constructed variables. Use this to identify which variables are updating and which aren’t.
- `ERP_1yr_ahead.csv` contains time series on the 20 constructed models. Use this to identify which models are updating and which aren’t.
- `ERP_principal_components.csv` contains the ERP estimates at different time horizons. This is the main output.
- `ERP_1yr_ahead.csv` contains the 1-year ahead ERP forecasts taken from the `ERP_principal_components.csv`.

## 5	Running Code
1. Open the `run_all.m` file and change the variable `root_dir` to reflect the proper directory where the project code is stored. 
    ```
    % %    e.g. directory of file is stored in /home/rcerxr21/Work/ERP
    root_dir = [filesep ‘home’ filesep ‘rcerxr21’ filesep ‘Work’ filesep ‘TIC’] 
    ```
2. After creating an account with FRED and receiving an API key, simply assign your API key to the variable `api` in the `main.m` file 

3. After creating a WRDS account, go to the WRDS web portal [here](https://wrds-cloud.wharton.upenn.edu/SASStudio/), you will be prompted to enter credentials for SAS, enter your user name and password for your WRDS account. Now create a folder structure that matches the following outline:
    ```
    Your primary folder (subject to user account)
    -> smitch15			     
       ---> sasuser.v94 
       The primary WRDS storage file
       ---> WRDS			    
          - Code		    Will be storing .sas files for downloading data
          - Output		    Will be storing .txt .csv and .xlsx files
    ```
4. After constructing the folder structure demonstrated above, go to the run_all.m file and update the variables `wrds_code_dir` and `wrds_out_dir` to match the location of those two folders. In addition, change `wrds_username` and `wrds_password` to match your WRDS username and password.

5. In the `run_all.m` file we are going to update the beaURL variable that reflects National Accounts (NIPA) data. Go to the Bureau of Economic Analysis website [here](https://apps.bea.gov/histdata/histChildLevels.cfm?HMI=7). Select (click) the most recent version of data (e.g. 2020, Q1) and then under the main tab right click on Section 2. Select “copy link address" and then paste that in place of the variable assignment.

6. Run the `haver_data.do` Stata file locally on the SAN, located in `Code\san_scripts\` this should produce a file named `tb1m.csv`. Now load up WinSCP on your machine, navigate to the folder where the data was downloaded on the SAN and copy that data to the `\Input` folder in the project

7. Go to the CFO survey website, located [here](https://cfosurvey.fuqua.duke.edu/release/). Go to the current year and quarter and select United States High-Level Results. Then find the part in the pdf where they report the expected average annual S&P500 return over the next 1 and 10 years. The precise questions are: *“Over the next 10 years, I expect the average annual S&P 500 return will be: Expected return"* and *“Over the next year, I expect the average annual S&P 500 return will be: Expected return."* Add the mean surveyed return and the number of total surveyed to the `\Input\cfo_data.csv` for those two questions

8. Run the `run_all.m` file in the RAN terminal, interactive environment or through the batch environment as follows:
    ```
    % %    e.g. running code via batch on the FRBNY RAN HPC Cluster
    $ matlab20a-batch-withemail 25 run_all.m 
    ```
10. Now after the ERP vintages are created check to make sure the output vintages are reasonable, by comparing them to previous versions stored in the `\Output` folder.

11. If everything matches up correctly, make a commit to your local branch of the ERP project and push to GitHub.
  
## 6	Potential Bugs
 - The NY_Fed_CM_DAPM.m file is subject to breaking on account of changes to the Shiller dataset, which has changed its data format in the past. In the event of such changes, one will be required to analyze the format and modify the code to fit the new formatting convention. 
 
- The WRDS download within the download_data.m file is subject to breaking when downloading the crsp_m.txt file. In the event of this occurrence, you will have to download the data manually as follows:
  1. Log onto your WRDS SAS Studio
  2. Find the file in your folder structure
  3. Right click on the file and download it
  4. Transfer it to the \Input folder
  
## 7	Potential Extensions
1. Automate the updates for the tb1m.csv (this may conflict with security privileges)
2. Automate the updates for the cfo_data.csv (this may require PDF parsing) 

## 8	Contributors
- Rajesh Rao (Research Analyst 22’)
- Charlie Smith (Research Analyst 21’)
- Eric Huang (Research Analyst 17’)
