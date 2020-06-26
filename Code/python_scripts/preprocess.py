import csv
import os
import re
import math
import scipy.io as sio
import numpy as np

os.chdir("/home/rcecrs02/ERP")

#find line breaks
#use first line after line break as title
#potentially should use re.split to catch double double spaces

def add_space(str, n):
	if n > 0:
		return str[0:n] + " " + str[n:]
	else:
		return str

def read_quarter(str):
	q = 0
	qmatch = re.search("(\w+) Quarter", str)
	try:
	    qstr = qmatch.group(1)
	    if qstr == "First" or qstr == "1st":
	        q = 1
	    elif qstr == "Second" or qstr == "2nd":
	        q = 2
	    elif qstr == "Third" or qstr == "3rd":
	        q = 3
	    elif qstr == "Fourth" or qstr == "4th":
	        q = 4
	except:
	    pass  

	seasonmatch = re.search("(Winter|Spring|Summer|Autumn|Fall|Year-End)", str)
	try:
	    seasonstr = seasonmatch.group(0)
	    if seasonstr == "Winter":
	        q = 1
	    elif seasonstr == "Spring":
	        q = 2
	    elif seasonstr == "Summer":
	        q = 3
	    elif seasonstr =="Autumn" or seasonstr == "Fall" or seasonstr == "Year-End":
	        q = 4
	except:
	    pass
	return q



def clean_CFO_file(surveyFileName, outfileName):
    read = 0
    block = ""
    surveyFile = open(surveyFileName)
    for f in surveyFile:
        # stop reading if you reach a new page
        if re.match("(.*)\\page(.*)", f):
            read = 0
        # keep reading otherwise by appending to block
        if read == 1:
            block = block + f
		# the header should read: "Duke CFO magazine Global Business Outlook survey" etc.	
        if re.match("(.*)Duke(.*)CFO(.*)", f):
            header = f
            nl = next(surveyFile)
            introline = next(surveyFile)					
        if re.match("(.*). On(.*)", introline) and re.match("(.*)yield(.*)", introline) and not (re.match("(.*)Weighted(.*)", introline)) and not (re.match("(.*)OUTLIERS(.*)", introline)) and not (re.match("(.*)Outliers(.*)", introline)):
            intro = introline
            read = 1
			# store quarter of survey
            ymatch = re.search("20\d{2}", header)
            year = ymatch.group(0)
            q = ad_quarter(header)
            dateq = year + 'Q' + str(q)
			
			# get treasury yield date and number from survey question                
            dmatch = re.search("[A-Za-z]+(\s+)[0-9]+[a-z,]*(\s*)\W[0-9]+", intro)
            datestr = dmatch.group(0).replace("nd","").replace("rd","").replace("th","")
            datestr = datestr.strip()

            yieldmatch = re.search("([0-9.]+)%", intro)
            tyield = yieldmatch.group(1)
		# when you hit the header of the question starting with the word 	 
		# "On" (which is supposed to be about yields), save it and 		
		# examine the next line	
        elif re.match("(.*). On(.*)", introline):
            global prev_intro
            prev_intro = introline
            introline = next(surveyFile)			
			# only read the non-weighted responses to the yield questions
            if re.match("(.*)yield(.*)", introline) and not (re.match("(.*)Weighted(.*)", prev_intro)) and not (re.match("(.*)OUTLIERS(.*)", prev_intro)) and not (re.match("(.*)Outliers(.*)", prev_intro)):
                # since each introline is divided into two lines in rtf,
                # you need to combine them together
                intro = prev_intro + introline
                read = 1	
                
                # store quarter of survey
                ymatch = re.search("20\d{2}", header)
                year = ymatch.group(0)
                q = read_quarter(header)
                dateq = year + 'Q' + str(q)
                
                # get treasury yield date and number from survey question                
                dmatch = re.search("[A-Z a-z]+(\s+)[0-9]+[a-z,]*(\s*)\W[0-9]+", intro)
                datestr = dmatch.group(0).replace("nd","").replace("rd","").replace("th","")
                datestr = datestr.strip()

                yieldmatch = re.search("([0-9.]+)%", intro)
                tyield = yieldmatch.group(1)
                
                #filename = "Input/Duke CFO Global Business Outlook Survey/" + dateq + ".csv"


    # create table from the block that was read
    if not(block == ""):
        with open('Data/Input/Duke CFO Global Business Outlook Survey/raw/tmp.txt', 'a+') as file:
		# this is a temporary txt file from which the excel
		# table will be built
            file.write(block)
            temp_table = list(csv.reader(open('Data/Input/Duke CFO Global Business Outlook Survey/raw/tmp.txt', 'rb'), delimiter='\t'))	
            new_table = []
            joinstr = ""
            exp_return_1y = ''
            exp_return_10y = ''
		# iterate through the block that was read
        for j in range(1,len(temp_table)):
            temp_table[j]
            a = "('%s')" %"".join(i)
			# the question sentences are divided into 2(or 3)
			# lines so we need to read both/all lines
            if re.match("(.*)next 10 years(.*)", a) or re.match("(.*)next year(.*)", a):
                b = "('%s')" %"".join(temp_table[j+1])
				# the sentence of this question is in 3 				# lines
				#if re.match("(.*)greater(.*)",b):
				#	b = b + "('%s')" %"".join(temp_table[j+2])
                c = a + b
				# delete some redundant stuff
                c = c.replace("(","")
                c = c.replace(")","")
                c = c.replace("'","")
                c = c.replace("\par","")
                if re.match("(.*)}d(.*)",c):
                    d = re.search("(.*)}d(.*)",c)
                    c = d.group(1)
                i = c.split("\\tab")
				
                i[0] = joinstr + ' ' + ''.join(i[0])
                joinstr = ""
				#new_table.append(i)
				#labelstr = ''.join(i[0])				
                if re.match("(.*)next 10 years(.*)", c) and not re.match("(.*)There is a 1-in(.*)", c):	
                    exp_return_10y = ''.join(i[1])
                elif re.match("(.*)next year(.*)", c) and not re.match("(.*)There is a 1-in(.*)", c):
                    exp_return_1y = ''.join(i[1])

        # create output file with dates and means
		#out_csv = csv.writer(open(filename, 'wb'))
		#out_csv.writerows(new_table)

        headlist = [dateq, datestr, tyield, exp_return_10y, exp_return_1y]
        with open(outfileName, 'ab') as mainfile:
            mainfile_csv = csv.writer(mainfile)
            duplicate_row = False
            for line in open(outfileName):
                if dateq in line:
                    duplicate_row = True		
                if not duplicate_row:
                    mainfile_csv.writerow(headlist)
		
        return dateq
			

def clean_French_file(f):
	filename = "Data/Input/Fama French/raw/" + f
	tables = []
	block = []
	#read space delimited line-by-line
	for line in open(filename):
		row = line.split()
		#if there is a line break, start a new table
		if len(row) == 0 and not block == []:
			# add table if last row first column is numeric, otherwise must be header/footer which we just ignore
			xval = ''.join(e for e in block[len(block)-1][0] if e.isalnum())
			if str.isdigit(xval):
				tables.append(block)
			#ntab = ntab + 1
			block = []
		#otherwise, append rows until there is a line break
		elif not len(row)==0:
			# code missing values as -999
			row[:] = [re.sub("\*\*\*\*\*\*\*", " -999 ", r) for r in row] 
			# stupid exception for when there is actually no space between two entries because they used the thousands place
			# e.g. 25_Portfolios_ME_Prior_60_13.txt
			row[:] = [add_space(r, max(m.start()-4 for m in re.finditer("[.]", r)) if len(re.findall("[.]", r))==2 else 0) for r in row]
			row = (' '.join(r for r in row)).split()
			block.append(row)

	#iterate through all tables
	for tab in tables:
		values = []
		header2 = []
		# any row that is numeric in its first element we consider a row of values
		# anything else goes into "header2"
		for row in tab[1:]:
			if str.isdigit(''.join(row[0])):
				row = map(float, row)
				values.append(row)
			else:
				header2.append(row)
				if "Lo" in row:
					lo_index = row.index("Lo")
					row[lo_index] = "Lo" + row[lo_index + 1]
					row.remove(row[lo_index + 1])
					hi_index = row.index("Hi")
					row[hi_index] = "Hi" + row[hi_index + 1]
					row.remove(row[hi_index + 1])

		# possible patterns are: 
		# Case 1: first row: title, second row: header
		# Case 2: first row: header, if there is no title because there
		# is only one table in the file
		
		# We name the file for the table the filename + table title
		
		# Case 2
		if str.isdigit(''.join(tab[1][0])):
			smallfile = "Data/Input/Fama French/"+f[:-4]
			header = [tab[0]]
		# Case 1
		else:
			title = ''.join(tab[0])
			title = ''.join(e for e in title if e.isalnum())
			smallfile = "Data/Input/Fama French/" + f[:-4] + "_" + title
			header = header2

		
		assert len(header) == 1 or len(header) == 2  
		# Sometimes there are 2-line headers, in this case
		# We apply the headings down, (e.g. North_America_6_Portfolios_ME_Prior_12_2.txt)
		newhead = ['']
		nvals = len(values[0])
		if len(header) == 2:
			n = len(header[0])
			m = (nvals - 1) / n
			for i in range(nvals-1):
				j = int(math.floor(i / m))
				newhead.append(header[0][j] + "_" + header[1][i])
		else:
			# sometimes headers are 2 words each
			if len(header[0]) == (nvals-1)*2:
				for i in range(nvals-1):
					newhead.append(header[0][2*i] + " " + header[0][2*i+1])
			else:
				newhead = newhead + header[0][:]

		#write to csv and mat
		out_csv = csv.writer(open(smallfile + ".csv", 'wb'))
		out_csv.writerows([newhead])
		out_csv.writerows(values)
		headarr = np.asarray(newhead, dtype=object)
		valmat = np.matrix(values)
		struct = {'header': headarr, 'values': valmat}
		sio.savemat(smallfile+".mat", {'data':struct})

def main():
	# all French files
#	for f in os.listdir("Input/Fama French/raw"):
#		if re.search(".txt", f):
#			try:
#				clean_French_file(f)
#			except IndexError:
#				print f, " Index Error"
#			except AssertionError:
#				print f, " Assertion Error"

    # Duke CFO files
	outfileName = "Data/Input/Duke CFO Global Business Outlook Survey/clean/duke_erp.csv"
	with open(outfileName, 'wb') as outfile:
        	outfilecsv = csv.writer(outfile)
        	outfilecsv.writerow(['dateq', 'tyield_date', 'tyield', '10-yr expected return', '1-year expected return'])
	for f in os.listdir("Data/Input/Duke CFO Global Business Outlook Survey/raw"):
		if re.match("(.*)US-Topline(.*)", f):
            		#try:
				fnm = "Data/Input/Duke CFO Global Business Outlook Survey/raw/" + f
				dateq = clean_CFO_file(fnm, outfileName)
				file_name = "Data/Input/Duke CFO Global Business Outlook Survey/raw/US-Topline_" + dateq + ".rtf"
				if f == 'US-Topline_download.rtf':
					orig_name = "Data/Input/Duke CFO Global Business Outlook Survey/raw/US-Topline_download.rtf"
					try:
						os.rename(orig_name,file_name)
					except WindowsError:
						os.remove(file_name)
						os.rename(orig_name,file_name)
					
				os.remove("Data/Input/Duke CFO Global Business Outlook Survey/raw/tmp.txt")
            		#except:
               			#print f, " Error", 
		

if __name__ == '__main__':
    main()
