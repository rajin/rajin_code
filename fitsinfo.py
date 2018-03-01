from astropy.io import fits
import glob

KEYWORDS = ['NAXIS2', 'LAMPID', 'GRATING', 'GR-ANGLE', 'OBJECT']
keywordValues = []

def get_fits_info(filename):
    header = fits.getheader(filename)
    for keyword in KEYWORDS:
    	try:
    		keywordValues.append(header[keyword])
    	except KeyError:
		print('Header keyword %s does not exist in %s' % (keyword, filename))
    		keywordValues.append('None')

    return keywordValues

files = glob.glob('*fits')
for file in files:
	fileInfo = get_fits_info(file)
	print(fileInfo)
