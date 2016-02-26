#!/usr/bin/env python

import os, sys, numpy, pyfits, scipy

import scipy, scipy.ndimage


def subsky(fitsfile, outputfile):

    hdu = pyfits.open(fitsfile)

    data = hdu['SCI'].data

    #
    # Apply median filter in spatial direction
    # this should get rid of the continuum without affecting the skylines
    #
    
    
    data_filtered = scipy.ndimage.filters.median_filter(
        input=data, 
        size=(1,75), 
        footprint=None, 
        output=None, 
        mode='reflect', 
        cval=0.0, 
        origin=0)

    skylines = hdu['SCI'].data - data_filtered

    # now combine all skylines
    skyspec = numpy.mean(skylines, axis=0)

    print skylines.shape
    print skyspec.shape

    skyspec_2d = skyspec.reshape((1, -1))
    print skyspec_2d.shape

    contnorm = numpy.mean(data_filtered, axis=0).reshape((1, -1))
    print contnorm.shape

    #skysub = (hdu['SCI'].data - skyspec_2d) / contnorm
    skysub = (hdu['SCI'].data - skyspec_2d)

    hdu['SCI'].data = skysub

    # fitsfile[:-5]+".skysub2.fits"
    hdu.writeto(outputfile, clobber=True)


if __name__ == "__main__":

    for fitsfile in sys.argv[1:]:
        output = fitsfile[:-5]+".skysub2.fits"
        subsky(fitsfile, output)





    
