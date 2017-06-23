from __future__ import print_function
import flipper
from flipper import liteMap
from astropy.io import fits
import orphics.tools.io as io
from scipy.ndimage.interpolation import zoom
import numpy as np
from orphics.tools.stats import bin2D
from orphics.tools.io import Plotter,dictFromSection,listFromConfig
from orphics.analysis import flatMaps as fmaps
import numpy as np

#hdulist=fits.open('../ACTdata/ACTPol.fits')
catalog= None #hdulist[1].data
lmap=liteMap.liteMapFromFits('../ACTdata/ACTPol_148_D56_PA1_S2_1way_I.fits')
RAs=[]
DECs=[]
widthStampArcminute = 20.
pixScale = 0.5
Nrands = 10000

ra_range = [0.,360.]
dec_range = [-10.,6.5]

def profiles(lite_map,width_stamp_arcminute,pix_scale,catalog=None,n_random_points=None,rand_ra_range=None,rand_dec_range=None):
        width_stamp_degrees = width_stamp_arcminute /60.
        Np = np.int(width_stamp_arcminute/pix_scale+0.5)
        print ("Expected width in pixels = ", Np)

        lmap = lite_map
	stack=0
	N=0

        if catalog is not None:
                looprange = range(0,len(catalog))
                assert n_random_points is None
                random = False
        else:
                assert n_random_points is not None
                assert len(rand_ra_range)==2
                assert len(rand_dec_range)==2
                looprange = range(0,n_random_points)
                random = True
        
	for i in looprange:
		if random:
                        ra = np.random.uniform(*rand_ra_range)
                        dec =np.random.uniform(*rand_dec_range)
			#ra=RArand[i]
			#dec=DECrand[i]
		else:
                        ra=catalog[i][1]
			dec=catalog[i][2]
		ix, iy = lmap.skyToPix(ra,dec)
		if ix>=20 and ix<lmap.Nx-20 and iy>=20 and iy<lmap.Ny-20:
			print(ra,dec)
			smap = lmap.selectSubMap(ra-width_stamp_degrees/2.,ra+width_stamp_degrees/2.,dec-width_stamp_degrees/2.,dec+width_stamp_degrees/2.)
                	cutout = zoom(smap.data.copy(),zoom=(float(Np)/smap.data.shape[0],float(Np)/smap.data.shape[1]))
                	print (cutout.shape)
			stack = stack + cutout
			xMap,yMap,modRMap,xx,yy = fmaps.getRealAttributes(smap)
			N=N+1.
	stack=stack/N
	#print(stack.shape())
	#print(smap.data.shape)
	print(stack)
	print(N)
	io.quickPlot2d(stack,"stackrandom.png")

	dt = pix_scale
	arcmax = 20.
	thetaRange = np.arange(0.,arcmax,dt)
	breal = bin2D(modRMap*180.*60./np.pi,thetaRange)
	cents,recons = breal.bin(stack)
	pl = Plotter(labelX='Distance from Center (arcminutes)',labelY='Temperature Fluctuation ($\mu K$)', ftsize=10)
	pl.add(cents,recons)
	pl._ax.axhline(y=0.,ls="--",alpha=0.5)
	pl.done("randomprofiles.png")
        return stack, cents, recons

stack, cents, recons = profiles(lmap,widthStampArcminute,pixScale,n_random_points=Nrands,rand_ra_range=ra_range,rand_dec_range=dec_range)


