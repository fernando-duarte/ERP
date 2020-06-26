import urllib.request, re, os, zipfile
from bs4 import BeautifulSoup
def dlfile(url, opener, dataPath):
    # Open the url
    try:
        f = opener.open(url)

        # Open our local file for writing
        locFile 		= 	dataPath+os.path.basename(url)
        with open(locFile, "wb") as local_file:
            local_file.write(f.read())

        zipCheck 		= 	re.compile('\.zip')
        if zipCheck.search(locFile) is not None:
        	z 				= 	zipfile.ZipFile(locFile)
        	z.extractall(dataPath)
        	z.close()
        	os.remove(locFile)
            
    #handle errors
    except:
        print("HTTP Error:" [url])

def grabAllData(source,urlRoot,indexPage,dataPath,opener):
	url 				= 	urlRoot + indexPage
	doc 				= 	opener.open(url)

	soup 				= 	BeautifulSoup(doc,'lxml')
	links 				= 	[link.get('href') for link in soup.find_all('a')]

	zips 				= 	re.compile('\.zip|\.txt|\.xls|\.xlsx|\.csv|\.rtf')
	pullList 			= 	[]
	for l in links:
		if l is not None:
			if zips.search(l) is not None:
				pullList.append(l.replace(urlRoot,''))

	dataPathSource 		= 	dataPath + source
	if not os.path.exists(dataPathSource):
   		os.makedirs(dataPathSource)

	print('\n1. Downloading data from ' + source + ':')
	[dlfile(urlRoot+pull, opener, dataPathSource) for pull in pullList]

def main():

	dataPath 			= 	'pyData/'
	
	urlList 			= 	[['','http://mba.tuck.dartmouth.edu/pages/faculty/ken.french/','data_library.html'], 
						 ['','http://www.econ.yale.edu/~shiller/','data.htm'],
						 ['','http://faculty.haas.berkeley.edu/lettau/','data_cay.html'],
						 ['','http://www.stern.nyu.edu/~adamodar/pc/datasets/',''],
						 ['','http://www.cfosurvey.org/','']
						 ]

	proxy_handler 		= 	urllib.request.ProxyHandler()
	proxy_auth_handler 	= 	urllib.request.HTTPBasicAuthHandler()
	opener 				= 	urllib.request.build_opener(proxy_handler,proxy_auth_handler)
	

	[grabAllData(group[0],group[1],group[2],dataPath,opener) for group in urlList]
	#[dlfile('http://www.federalreserve.gov/econresdata/researchdata/feds200628.xls',urllib.request.build_opener(urllib.request.ProxyHandler(),urllib.request.HTTPBasicAuthHandler()),'pyData/')]
	#[dlfile('http://us.spindices.com/documents/additional-material/sp-500-eps-est.xlsx',urllib.request.build_opener(urllib.request.ProxyHandler(),urllib.request.HTTPBasicAuthHandler()),'pyData/')]

if __name__ == '__main__':
    main()	
