'''Main code that should be run for the spectral reduction of SALT longslit data
Please consult the readme file to understand how it works.

'''

import numpy as np
#import matplotlib as plt
from pyraf.iraf import noao
#from pyraf.iraf import stsdas
from pyraf import iraf
import pyfits as pft
import os,sys,glob,string,time
import math as mth


#------------------------------------------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------------------------------------------
#Function which applies the following iraf task on an ARC image: Identify, Reidentify, Fitcoords, transform

def identify(properties,index,irafhome):
   satis = 'n'
   while satis == 'n':
  properties = np.array(properties)
	print properties
	filename = properties[index,0]
	lampid = properties[index,2]
	path = irafhome+lampid+'.txt'
 	print path,filename
	iraf.noao.twodspec.longslit.identify(images=filename, section='middle line', databas='database',
			            	coordli=path, units='',nsum=10,match=-3., maxfeat=50, zwidth=100.,ftype='emission',
				    	fwidth=5., cradius=6., thresho=0., minsep=2.,functio='spline3',order=3,sample= '*',
				    	niterat=0, low_rej=3., high_re=3., grow=0.,autowri='no', graphic='stdgraph', cursor='',  	       			    	crval='',cdelt='',aidpars='',mode='ql')
	iraf.noao.twodspec.longslit.reidentify(referenc=filename,images=filename,interac='no',section='middle line',newaps='yes',
					overrid='no',refit='no',trace='no',step=10.,nsum=10.,shift=0.,search=5.,nlost=10.,
					cradius=5.,thresho=0.,addfeat='no',coordli=path,match=-3,maxfeat=50,minsep=2,
					databas='database',logfile='logfile',plotfil='',verbose='yes',graphic='stdgraph',
					cursor='',answer='yes',crval='',cdelt='',aidpars='',mode='ql')
	namesplit = string.split(filename,'.')
	iraf.noao.twodspec.longslit.fitcoords(images=str(namesplit[0]),fitname='',interac='yes',combine='no',
					databas='database',deletio='deletions.db',
					functio='chebyshev',xorder=6,yorder=6,logfile='STDOUT,logfile',plotfil='plotfile',
					graphic='stdgraph',cursor ='',mode='ql')
	trans = 't'+filename
	iraf.noao.twodspec.longslit.transform(input=filename,output=trans,minput='',moutput='',fitnames=namesplit[0],
					databas='database',
					interpt='spline3',x1='INDEF',x2='INDEF',dx='INDEF',nx='INDEF',xlog='no',y1='INDEF',
					y2='INDEF',dy='INDEF',ny='INDEF',ylog='no',flux='yes',blank='INDEF',
					logfile='STDOUT,logfile',mode='ql')
	os.system('ds9 %s &' % (trans))
	satis = str(raw_input("Are satisfied with the transformed spectra (y|n)?"))
	if (satis == 'n'):
	  askdel = str(raw_input("Delete the spectra and database (y|n)?"))
	  if (askdel == 'y'):
		os.system('rm -r database/ %s deletion.db' % (trans))
	os.system('mv %s %s history/' % (filename,trans))
	return
    
#------------------------------------------------------------------------------------------------------------------------------------

def ccdgap(name):
	fimg = pft.open(name)
	prihdr = fimg[0].header
	scidata = fimg[0].data
	
	n1 = prihdr['NAXIS1']
	n2 = prihdr['NAXIS2']

	a,b = 1023, 1075 #variables to store the ccd 1st gap coord
	c,d = 2100, 2147 #variables to store the ccd 2nd gap coord
	e = n2 #ccd height

	for i in range(0,e):
	   avside1 = avside2 = avside3 = avside4 = 0.0
	   grad1 = grad2 = c1 = c2 = 0.0
	   for j in range(1,6):
	      avside1 = avside1 + (scidata[i,a-(j+2)]/5.0)
	      avside2 = avside2 + (scidata[i,b+(j+2)]/5.0)	
	      avside3 = avside3 + (scidata[i,c-(j+2)]/5.0)
	      avside4 = avside4 + (scidata[i,d+(j+2)]/5.0)

	   grad1 = (avside2-avside1)/((b-a)+5.0)
	   grad2 = (avside4-avside3)/((d-c)+5.0)
	   c1 = scidata[i,a-5]
	   c2 = scidata[i,c-5]

	   for j in range(0,(b-a)+1):
	      scidata[i,a+(j-1)] = grad1*(j+5.)+avside1
	   for j in range(0,(d-c)):
	      scidata[i,c+(j-1)] = grad2*(j+5.)+avside3

	namec = "c"+name
	pft.writeto(namec,data=scidata,header=prihdr,clobber=True)
	fimg.close()
	os.system('mv %s history/' % (name))
	return

#------------------------------------------------------------------------------------------------------------------------------------

def lacosmic(name,irafhome):
	import time
	iraf.task(lacos_spec=irafhome+'lacos_spec.cl')
	outname = 'la'+name
	pl = 'mask'+name
	iraf.lacos_spec(input=name,output=outname,outmask=pl,gain=1.,readn=2.89,
				xorder=9,yorder=0,sigclip=4.5,sigfrac=0.5,objlim=1.,niter=7,verbose='yes',mode='al')
	old = time.time()
	os.system('ds9 %s %s &' % (pl,outname))
	os.system('mv %s history/' % (name))
	return

#------------------------------------------------------------------------------------------------------------------------------------

def science(sciname,filename):
	# taking the arc file name and using it as the reference input for transformation
	namesplit = string.split(filename,'.')
	trans = 't'+sciname
	iraf.noao.twodspec.longslit.transform(input=sciname,output=trans,minput='',moutput='',fitnames=namesplit[0],
					databas='database',
					interpt='linear',x1='INDEF',x2='INDEF',dx='INDEF',nx='INDEF',xlog='no',y1='INDEF',
					y2='INDEF',dy='INDEF',ny='INDEF',ylog='no',flux='yes',blank='INDEF',
					logfile='STDOUT,logfile',mode='ql')
	os.system('ds9 %s &' % (trans))
	backn = 'b'+trans
	iraf.noao.twodspec.longslit.background(input=trans,output=backn,axis=2,interac='yes',sample='*',naverag=1,
					functio='spline3',order=3,low_rej=2.,high_re=1.5,niterat=1,grow=0.,
					graphic='stdgraph',cursor='',mode='al')
	os.system('ds9 %s &' % (backn))
	os.system('mv %s history/' % (sciname))
	os.system('mv %s history/' % (trans))
	return

#------------------------------------------------------------------------------------------------------------------------------------

def apall(apname,jj,refnam):
   satis = 'n'
   while (satis == 'n'):
	outname = 'til'+apname
	refname=''
	interap ='yes'
	if (jj>1):
		apallans = str(raw_input("Do you want to apply the same solution as the previous frame (y|n)?"))
		if (apallans == 'y'):
			refname=refnam
			interap ='no'

	iraf.noao.twodspec.apextract.apall(input=apname,output=outname,apertur='',format='strip',referen=refname,profile='',
						interac=interap,find=interap,recente=interap,resize=interap,edit=interap,trace=interap,
						fittrac=interap,extract='yes',extras='no',review='no',line=500,nsum=30,
						lower=-5.,upper=12.,apidtab='',b_funct='chebyshev',b_order=1,b_sampl='-60:-40,40:60',
						b_naver=-3,b_niter=0,b_low_r=3.,b_high_=3.,b_grow=0.,width=5.,radius=10.,
						thresho=0.,nfind=1,minsep=35.,maxsep=100000.,order='increasing',aprecen='',
						npeaks='INDEF',shift='yes',llimit=-15,ulimit=15,ylevel=0.1,peak='yes',
						bkg='yes',r_grow=0.,avglimi='no',t_nsum=45,t_step=30,t_nlost=10,t_funct='spline3',
						t_order=3,t_sampl='*',t_naver=1,t_niter=0,t_low_r=3.,t_high_=3.,t_grow=0.,
						backgro='none',skybox=1,weights='none',pfit='fit1d',clean='no',saturat='INDEF',
						readnoi=0.,gain=1.,lsigma=4.,usigma=4.,nsubaps=1.,mode='ql')
	namesplit = string.split(outname,'.')
	os.system('ds9 %s &' % (namesplit[0]+'.0001.fits'))
	satis = str(raw_input("Are satisfied with the extracted spectra (y|n)?"))
	if (satis == 'n'):
		askdel = str(raw_input("Delete the spectra (y|n)?"))
		if (askdel == 'y'):
			os.system('rm %s' % (namesplit[0]+'.0001.fits'))
	else:
		os.system('mv %s history/' % (apname))
	return

#------------------------------------------------------------------------------------------------------------------------------------
def err(name,name1):
	ffts = pft.open(name)
	prihdr = ffts[0].header
	scidata = ffts[0].data
	ERROR = 'ERROR'
	prihdr.update('OBJECT', ERROR)
	scidata = scidata+(2.89**2.0)
	scidata = abs(scidata)
	scidata1 = scidata**0.5
	ffts1 = pft.writeto(name1,scidata1,header=prihdr,clobber=True)
	ffts.close()
	return
	

#------------------------------------------------------------------------------------------------------------------------------------

def transform(name,filename):
	# taking the arc file name and using it as the reference input for transformation
	namesplit = string.split(filename,'.')
	trans = 't'+name
	iraf.noao.twodspec.longslit.transform(input=name,output=trans,minput='',moutput='',fitnames=namesplit[0],
					databas='database',interpt='linear',x1='INDEF',x2='INDEF',dx='INDEF',nx='INDEF',
					xlog='no',y1='INDEF',y2='INDEF',dy='INDEF',ny='INDEF',ylog='no',flux='yes',blank='INDEF',
					logfile='STDOUT,logfile',mode='ql')
	os.system('mv %s history/' % (name))
	return

#------------------------------------------------------------------------------------------------------------------------------------

def divide(name1,name2,switch):
	ffts1 = pft.open(name1)
	prihdr1 = ffts1[0].header
	scidata1 = ffts1[0].data
	ffts2 = pft.open(name2)
	prihdr2 = ffts2[0].header
	scidata2 = ffts2[0].data
	x,y  = np.where(scidata2==0)
	for i in range(0,len(x)):
		scidata2[x[i],y[i]] = 1.
	scidata3 = abs(scidata1/scidata2)
	if switch == 1:
		new = 'd'+name1
	else:
		new = 'fl'+name1

	ffts3 = pft.writeto(new,scidata3,header=prihdr1,clobber=True)
	ffts1.close()
	ffts2.close()
	#os.system('mv %s history/' % (name1))
	#os.system('mv %s history/' % (name2))
	return

#------------------------------------------------------------------------------------------------------------------------------------
#function to rename fits file with name which makes much more sense. Also SALT fits files consists usually of 2 layers of header and 1
#layer of data. This function makes the files become single layer with a header which has all the important information
def rename(fits_name):
	try:
		fimg = pft.open(fits_name)
		prihdr = fimg[0].header
		prihdr1 = fimg[1].header
		scidata = fimg[1].data

		n = prihdr1['NAXIS']
		n1 = prihdr1['NAXIS1']
		n2 = prihdr1['NAXIS2']
		object_name = prihdr['OBJECT']
		split_name = object_name.split()
		new_name = ''
		for i in range(0,len(split_name)):
			if i == 0: new_name = new_name+split_name[i]
			else: new_name = new_name+'_'+split_name[i]

		index = fits_name.split('.')[0][-3:]
		if len(new_name)>=11: new_name = new_name[:11]+'_'+index+'.fits'
		else: new_name = new_name+'_'+index+'.fits'

		telalt = prihdr['TELALT']
		airmass = 1.0/mth.sin(mth.radians(telalt))
		prihdr.update('airmass', airmass)

		prihdr.update('NAXIS', n);prihdr.update('NAXIS1', n1);prihdr.update('NAXIS2', n2)
		pft.writeto(new_name,data=scidata,header=prihdr,clobber=True)
		out = 1
	except:
		fimg = pft.open(fits_name)
		prihdr = fimg[0].header
		scidata = fimg[0].data

		object_name = prihdr['OBJECT']
		split_name = object_name.split()
		new_name = ''
		for i in range(0,len(split_name)):
			if i == 0: new_name = new_name+split_name[i]
			else: new_name = new_name+'_'+split_name[i]
	
		index = fits_name.split('.')[0][-3:]
		if len(new_name)>=11: new_name = new_name[:11]+'_'+index+'.fits'
		else: new_name = new_name+'_'+index+'.fits'

		telalt = prihdr['TELALT']
		airmass = 1.0/mth.sin(mth.radians(telalt))
		prihdr.update('airmass', airmass)

		pft.writeto(new_name,data=scidata,header=prihdr,clobber=True)
		out = 2
	return out

#------------------------------------------------------------------------------------------------------------------------------------
#function to read the header of fits files and output a few critical info needed in the program
def info_fits(fits_name):
	dummy_list = []
	fimg = pft.open(fits_name)
	prihdr = fimg[0].header
	scidata = fimg[0].data

	naxis2 = prihdr['NAXIS2']
	try:lampid = prihdr['LAMPID']
	except:lampid = 'NONE'
	try:grating = prihdr['GRATING']
	except:grating = 'NONE'
	try:gr_angle = prihdr['GR-ANGLE']
	except:gr_angle = 'NONE'
	try:object_n = prihdr['OBJECT']
	except:object_n = 'NONE'	
	dummy_list = [fits_name,naxis2,lampid,grating,gr_angle,object_n]
	return dummy_list
#------------------------------------------------------------------------------------------------------------------------------------
def input_val(script,lim1,lim2):
	while True:
	   try:
	       value = int(raw_input('%s : ' % (script)))
	   except ValueError: # just catch the exceptions you know!
	       print 'That\'s not a number!'
	   else:
	       if lim1 <= value < lim2: # this is faster
	           break
	       else:
	           print 'Out of range. Try again'
	return value

#------------------------------------------------------------------------------------------------------------------------------------

def flat(flist):
	f = open('flatlist','w')
	for i in range(0,len(flist)):
		f.write('%s\n' % (flist[i]))
	f.close()
	f_list = '@flatlist'
	combine_flat = 'com_flat.fits'

	iraf.images.immatch.imcombine(input=f_list,output=combine_flat,headers='',bpmasks='',rejmask='',nrejmas='',expmask='',
					sigmas ='',imcmb='$I',logfile='STDOUT',combine='median',reject ='none',project='no',
					outtype='real',outlimi='',offsets='none',masktyp='none',maskval=0,blank=0.,scale='none',
					zero='none',weight='none',statsec='',expname='',lthresh='INDEF',hthresh='INDEF',nlow=1,nhigh=1,
					nkeep=1,mclip='yes',lsigma=3.,hsigma =3.,rdnoise=0.,gain=1.,snoise =0.,sigscal=0.1,
					pclip=-0.5,grow=0.,mode='ql')

	illum_flat = 'il'+combine_flat

	iraf.noao.imred.ccdred.ccdproc.noproc='no'
	iraf.noao.imred.ccdred.ccdproc.fixpix='no'
	iraf.noao.imred.ccdred.ccdproc.oversca='no'
	iraf.noao.imred.ccdred.ccdproc.trim='no'
	iraf.noao.imred.ccdred.ccdproc.zerocor='no'
	iraf.noao.imred.ccdred.ccdproc.darkcor='no'
	iraf.noao.imred.ccdred.ccdproc.flatcor='no'
	iraf.noao.imred.ccdred.ccdproc.illumco='yes'
	iraf.noao.imred.ccdred.ccdproc.fringec='no'
	iraf.noao.imred.ccdred.ccdproc.readcor='no'
	iraf.noao.imred.ccdred.ccdproc.scancor='no'

	iraf.noao.imred.ccdred.mkillumflat(input=combine_flat,output=illum_flat,ccdtype='',xboxmin=3.,xboxmax=0.25,yboxmin=3.,
					yboxmax=0.25,clip='yes',lowsigm=2.5,highsig=2.5,divbyze=1.,ccdproc='',mode='ql')
	ccdgap(illum_flat)
	cillum_flat = 'c'+illum_flat
	os.system('ds9 %s &' % (cillum_flat))
	fimg = pft.open(cillum_flat)
	prihdr = fimg[0].header
	scidata = fimg[0].data
	n1 = prihdr['NAXIS1']
	n2 = prihdr['NAXIS2']
	satis = 'n'
	while satis == 'n' or satis =='no':
		y0 = input_val('Enter the first y-coord from where the mean will be computed, y0',0,n2)
		y1 = input_val('Enter the second y-coord from where the mean will be computed, y1',y0,n2)
		x0 = input_val('Enter the first x-coord from where the mean will be computed, x0',0,n1)
		x1 = input_val('Enter the second x-coord from where the mean will be computed, x1',x0,n1)
		satis = str(raw_input('Are you satisfied with the coordinates input(y or n): '))
	mean_flux = np.mean(scidata[y0:y1,x0:x1])
	scidata1 = scidata/mean_flux
	new_name = 'master_flat.fits'
	pft.writeto(new_name,data=scidata1,header=prihdr,clobber=True)
	for i in range(0,len(flist)):
		os.system('mv %s history/' % (flist[i]))
	os.system('mv %s history/' % (combine_flat))
	os.system('mv %s history/' % (illum_flat))
	os.system('mv %s history/' % (cillum_flat))
	return
	
#------------------------------------------------------------------------------------------------------------------------------------
def trim(imname,x,y):
	fimg = pft.open(imname)
	prihdr = fimg[0].header
	scidata = fimg[0].data

	newscidata = scidata[y[0]:y[1],x[0]:x[1]]
	prihdr.update('NAXIS1', (x[1]-x[0]))
	prihdr.update('NAXIS2', (y[1]-y[0]))

	new_name = 'tr'+imname
	pft.writeto(new_name,data=newscidata,header=prihdr,clobber=True)
	os.system('mv %s history/' % (imname))

#------------------------------------------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------------------------------------------

irafhome = '/home/rajin/iraf/'

# opening a file
os.system('mkdir history')
files = glob.glob('*fits')
#renaming the files to more appropriate names!
for i in range(0,len(files)):
	out = rename(files[i])
	if out == 1:
		os.system('mv %s history/' % (files[i]))

#Applying trimming if necessary
files = glob.glob('*fits')
os.system('ds9 %s &' % (files[0]))
trim_ans = str(raw_input("Do you want to apply any trimming (y|n)?"))
if trim_ans == 'y':
	satis = 'n'
	while satis == 'n' or satis =='no':
		y0 = input_val('Enter the first y-coord, y0',0,2500)
		y1 = input_val('Enter the second y-coord, y1',y0,2500)
		x0 = input_val('Enter the first x-coord, x0',0,3180)
		x1 = input_val('Enter the second x-coord, x1',x0,3180)
		satis = str(raw_input('Are you satisfied with the coordinate input(y or n): '))
	x = [x0,x1] ; y = [y0,y1]
	for i in range(0,len(files)):
		trim(files[i],x,y)

#saving the properties of the fits files in an array
files = glob.glob('*fits')
properties = [] 
for i in range(0,len(files)):
	properties.append(info_fits(files[i]))
properties = np.array(properties) # array with all the properties of the data observed
	
flat_detec = 'no'
science_list = [] ; flat_list = []
#creating list for each types of files: science, flats.
for j in range(0,len(properties)):
	test = properties[j,5]
	if test == 'ARC':
		ident = j
		arc_file = properties[j,0]
		print 'The program identifies that there is an arc lamp image in the input files!'
	elif test == 'FLAT':
		flat_detec = 'yes'
		flat_list.append(properties[j,0])
	elif test != 'ERROR':
		science_list.append(properties[j,0])

print 'This program aims to do the spectral reduction for salt data'

#This section will use the function identify to go through the identification of the Arc along with other functions to transform the 
#arc. To transform the arc, a function is created and this same function will later be applied on to the other (science) frames.
arc_lamp = str(raw_input("Do you want to identify the arc lamp (y|n)?"))
if (arc_lamp == 'y'):
	identify(properties,ident,irafhome)

prefix =''
#This section interpolates across ccd gap to get better results later when the background is subtracted
ccdgp = str(raw_input("Do you want to fill the ccd gaps of the science frames with a gradient function(y|n)?"))
if (ccdgp == 'y'):
	for j in range(0,len(science_list)):
		ccdgpname = prefix+science_list[j]
		ccdgap(ccdgpname)
	prefix = 'c'+prefix

#This section will use lacosmic to remove cosmic ray only on the science frames. Check the settings in the lacomic
#function of this code to your liking
lacos = str(raw_input("Do you want to treat science frame for cosmic ray removal with lacosmics (y|n)?"))	
mask = []
if (lacos == 'y'):
	for j in range(0,len(science_list)):
		lacosname = prefix+science_list[j]
		lacosmic(lacosname,irafhome)
		mask.append(lacosname)
	prefix = 'la'+prefix

# Detecting if there are flat field inside the folder and then applying flat field correction to the science image (at the end we will have 
# 2 sets of science images 1 flat fielded and 1 none flat fielded
if flat_detec == 'yes':
	flat(flat_list)
	for j in range(0,len(science_list)):
		divide(prefix+science_list[j],'master_flat.fits',2)
	flatprefix = 'fl'+prefix

for i in range(0,len(mask)):
	os.system('mv %s history/' % ('mask'+mask[i]))

#Here we have the start of the calculation for the error frames
errorans = str(raw_input("Do you want to  compute error frames (y|n)?"))
error=[] ; errorflt = []
if (errorans == 'y'):
	for j in range(0,len(science_list)):
		name = prefix+science_list[j]
		namesplt = name.split('.')
		err_name = 'err_'+namesplt[0][-3:]+'.fits'
		err(name,err_name)
		error.append(err_name)
		if flat_detec == 'yes':
			flt_err_name = 'fl'+err_name
			flt_name = flatprefix+science_list[j]
			err(flt_name,flt_err_name)
			errorflt.append(flt_err_name)


print 'Applying wavelength calibration and background subtraction'
for j in range(0,len(science_list)):
	sciename = prefix+science_list[j]
	science(sciename,arc_file)
	if flat_detec == 'yes':
		flsciename = flatprefix+science_list[j]
		science(flsciename,arc_file)

prefix = 'bt'+prefix
if flat_detec == 'yes':
	flatprefix = 'bt'+flatprefix


prefixerr =''
if (len(error) != 0):
# transforming the error frames
	for j in range(0,len(error)):
		transform(error[j],arc_file)
	prefixerr ='t'+prefixerr
	i = 0
	for j in range(0,len(science_list)):
		sciename = prefix+science_list[j]
		errname = prefixerr+error[i]
	 	print(" %s / %s" % (sciename,errname))
		divide(errname,sciename,1)
		i = i+1
	prefixerr = 'd'+prefixerr

prefixerrflt =''
if (len(errorflt) != 0):
# transforming the flat fielded error frames
	for j in range(0,len(errorflt)):
		transform(errorflt[j],arc_file)
	prefixerrflt ='t'+prefixerrflt
	i = 0
	for j in range(0,len(science_list)):
		fltsciename = flatprefix+science_list[j]
		flterrname = prefixerrflt+errorflt[i]
		print(" %s / %s" % (fltsciename,flterrname))
		divide(flterrname,fltsciename,1)
		i = i+1
	prefixerrflt = 'd'+prefixerrflt
	

apname1=''
apallext = str(raw_input("Do you want to extract your 2D aperture and correct for tilt(y|n)?"))
if (apallext == 'y'):
	count = 0
	for j in range(0,len(science_list)):
		count = count+1
		apname = prefix+science_list[j]
		if (count==1): apname1 = apname
		apall(apname,count,apname1)
	prefix = 'til'+prefix

	if flat_detec == 'yes':
		for j in range(0,len(science_list)):
			fltapname = flatprefix+science_list[j]
			apall(fltapname,count,apname1)
  
#applying apall solution on the error frames
	if (len(error) != 0):
		for j in range(0,len(error)):
			apall(prefixerr+error[j],2,apname1)

	if (len(errorflt) != 0):
		for j in range(0,len(errorflt)):
			apall(prefixerrflt+errorflt[j],2,apname1)

