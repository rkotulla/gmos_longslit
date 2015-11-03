#!/usr/bin/env python

import os, sys, numpy, pyfits

if __name__ == "__main__":
    
    fitsfile = sys.argv[1]
    ap = [int(d) for d in sys.argv[2].split(":")]

    hdu = pyfits.open(fitsfile)

    data = hdu['SCI'].data

    cutout = data[ap[0]-1:ap[1], :]
    
    spec1d = numpy.mean(cutout, axis=0)

    print spec1d.shape

    idx = numpy.arange(spec1d.shape[0]) + 1
    hdr = hdu['SCI'].header

    wl = (idx - hdr['CRPIX1']) * hdr['CD1_1'] + hdr['CRVAL1']

    try:
        wlcorr = [float(d) for d in sys.argv[3].split("-")]
        d_wl = wlcorr[1] - wlcorr[0]
        wl += d_wl
    except:
        pass

    try:
        print "masking"
        blocks = sys.argv[4].split(",")
        for block in blocks:
            b12 = [int(d) for d in block.split(":")]
            if (b12[0] < 0):
                spec1d[(-b12[0]-1):(-b12[1])] = numpy.NaN
            else:
                wl_mask = (wl > b12[0]) & (wl < b12[1])
                spec1d[wl_mask] = numpy.NaN
    except:
        pass

    comb = numpy.append(wl.reshape((-1,1)),
                        spec1d.reshape((-1,1)), axis=1)
    print comb.shape

    # output as text file
    outfile = fitsfile[:-5]+".spec"
    with open(outfile, "w") as of:
        numpy.savetxt(of, comb)
        print >>of, "\n\n\n\n\n"

    #
    # also write as 1-d fits file
    #
    primhdu = pyfits.PrimaryHDU(
        data = spec1d,
        )
    primhdu.header['CD1_1'] = wl[1]-wl[0]
    primhdu.header['CRPIX1'] = 1.
    primhdu.header['CRVAL1'] = wl[0]
    primhdu.writeto(fitsfile[:-5]+".1d.fits", clobber=True)

    
