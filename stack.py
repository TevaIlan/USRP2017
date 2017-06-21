import flipper
from flipper import liteMap
from astropy.io import fits
from __future__ import print_function
import orphics.tools.io as io

hdulist=fits.open('ACTPol.fits')
catalog=hdulist[1].data
map=liteMap.liteMapFromFits('ACTPol_148_D56_PA1_S2_1way_I.fits')
RAs=[]
DECs=[]
widthStampArcminute = 20.
widthStampDegrees = widthStampArcminute /60
stack=0
N=0
for i in range(0,len(catalog)):
	ra=catalog[i][1]
	dec=catalog[i][2]
	ix, iy = map.skyToPix(ra,dec)
	# print ix, iy ; assert ix>=0 ; assert ix<map.Nx ; assert iy>=0 ; assert iy<map.Ny when this is uncommented we get an AssertionError which I guess means that some of the values are outside the map
	# if 53.983333<=ra and ra<=338.533333:
	if ix>=0 and ix<map.Nx and iy>=0 and iy<map.Ny:
		smap = map.selectSubMap(ra-widthStampDegrees,ra+widthStampDegrees,dec-widthStampDegrees,dec+widthStampDegrees)
		#SMaps.append(smap)
		stack = stack + smap.data.copy()
		N=N+1.
stack=stack/N
#print(stack.shape())
print(smap.data.shape)
print(stack)
print(N)
io.quickPlot2d(stack,"stack.png")


