#!/usr/bin/env python

import os, sys, numpy, pyfits
import bottleneck

def clobberfile(fn):
    if os.path.isfile(fn):
        os.remove(fn)
    return 

if __name__ == "__main__":
    
    fitsfile = sys.argv[1]
    ap = [int(d) for d in sys.argv[2].split(":")]

    hdu = pyfits.open(fitsfile)

    data = hdu['SCI'].data.astype(numpy.float32)
    hdr = hdu['SCI'].header
    dq = hdu['DQ'].data
    data[dq > 0] = numpy.NaN

    cutout = data[ap[0]-1:ap[1], :]
    clobberfile("optextract.fits")
    pyfits.PrimaryHDU(data=cutout, header=hdr).writeto("optextract.fits", clobber=True)

    ypos = numpy.arange(ap[0]-1, ap[1])+1

    skyregion_select = numpy.zeros((ypos.shape[0]), dtype=numpy.bool)
    if (len(ap) > 4):
        # we have some user-specified sky-regions to use
        skyregion_select[:] = False
        for ix in range(4, len(ap), 2):
            _y1 = ap[ix] + ap[2]
            _y2 = ap[ix+1] + ap[2]
            skyregion_select[(ypos>=_y1)&(ypos<=_y2)] = True
    # print skyregion_select

    apsize = ap[3]
    
    #
    # compute median spectrum in spatial direction - 
    # this is going to be used for sky-subtraction
    #
    skyspec = bottleneck.nanmedian(cutout, axis=0)
    # print skyspec.shape

    cutout_bg = cutout - skyspec
    clobberfile("optextract0.fits")
    pyfits.PrimaryHDU(data=cutout_bg, header=hdr).writeto("optextract0.fits", clobber=True)


    cutout_sky = numpy.array(cutout)
    cutout_sky[~skyregion_select] = numpy.NaN
    clobberfile("optextract_skysel.fits")
    pyfits.PrimaryHDU(data=cutout_sky).writeto("optextract_skysel.fits", clobber=True)
    cutout_sky_1d = bottleneck.nanmean(cutout_sky, axis=0)
    numpy.savetxt("spec_simple_sky", cutout_sky_1d)

    cutout_skysub = numpy.array(cutout) - cutout_sky_1d
    clobberfile("optextract_skysub.fits")
    pyfits.PrimaryHDU(data=cutout_skysub, header=hdr).writeto("optextract_skysub.fits", clobber=True)
    
    #
    # Apply 1-d median filter in wavelength direction to filter out cosmics
    #
    profile = bottleneck.nanmedian(cutout, axis=1)
    profile2 = bottleneck.nanmean(cutout, axis=1)

    profile_combined = numpy.zeros((profile.shape[0], 3))
    profile_combined[:,0] = ypos[:]
    profile_combined[:,1] = profile[:]
    profile_combined[:,2] = profile2[:]
    numpy.savetxt("lineprofile", profile_combined)

    cutout2 = cutout / profile.reshape((-1,1))
    clobberfile("optextract1.fits")
    pyfits.PrimaryHDU(data=cutout2, header=hdr).writeto("optextract1.fits", clobber=True)

    cutout3 = cutout * profile.reshape((-1,1))
    clobberfile("optextract2.fits")
    pyfits.PrimaryHDU(data=cutout3, header=hdr).writeto("optextract2.fits", clobber=True)

    # scaled_sky = bottleneck.nanmedian(cutout2[skyregion_select], axis=0)
    # numpy.savetxt("scaledsky", scaled_sky)
    # cutout_f = (cutout2 - scaled_sky) * profile.reshape((-1,1))
    # pyfits.PrimaryHDU(data=cutout_f).writeto("optextract2f.fits", clobber=True)

    scaled_sky = bottleneck.nanmedian(cutout[skyregion_select], axis=0)
    numpy.savetxt("scaledsky", scaled_sky)
    cutout_f = (cutout - scaled_sky) #* profile.reshape((-1,1))
    clobberfile("optextract2f.fits")
    pyfits.PrimaryHDU(data=cutout_f, header=hdr).writeto("optextract2f.fits", clobber=True)
    
    #
    # find peak position - assume center position given is good to +/- 2 pixels
    #
    _close = (ypos>ap[2]-2) & (ypos<ap[2]+2)
    _peak_i = numpy.argmax(profile[_close])
    _peak = ypos[_close][_peak_i]
    # print _peak_i, _peak

    #
    # find dimensions of aperture - include everythin up to the point where the
    # profile turns up again
    #
    if (apsize < 0):
        dprof = numpy.diff(profile)
        numpy.savetxt("diffprofile", dprof)
        out_of_profile = (ypos[:-1] > _peak) & (dprof > 0)
        _max = numpy.min(ypos[:-1][out_of_profile])-1

        out_of_profile = (ypos[:-1] < _peak) & (dprof < 0)
        _min = numpy.max(ypos[:-1][out_of_profile])
    else:
        _min = _peak - ((apsize-1)/2)
        _max = _peak + ((apsize-1)/2)

    print _min, _max, _max-_min+1

    profile_sel = numpy.zeros((ypos.shape[0]), dtype=numpy.bool)
    profile_sel[((ypos>=_min)&(ypos<=_max))] = True
    cutout_sel = numpy.array(cutout_f)
    cutout_sel[~profile_sel] = numpy.NaN
    clobberfile("optextract3.fits")
    pyfits.PrimaryHDU(data=cutout_sel, header=hdr).writeto("optextract3.fits", clobber=True)
    #pyfits.PrimaryHDU(data=(cutout*profile_sel.reshape((-1,1)))).writeto("optextract3.fits", clobber=True)


    simplespec = bottleneck.nanmean(cutout[profile_sel], axis=0)
    numpy.savetxt("spec_simple", simplespec)



    profile_bg = numpy.ones((ypos.shape[0]), dtype=numpy.bool)
    profile_bg[((ypos>=_min)&(ypos<=_max))] = False #0 #_min:_max
    cutout_bg = numpy.array(cutout)
    cutout_bg[~profile_bg] = numpy.NaN
    clobberfile("optextract4.fits")
    pyfits.PrimaryHDU(data=cutout_bg, header=hdr).writeto("optextract4.fits", clobber=True)
    #pyfits.PrimaryHDU(data=(cutout*profile_bg.reshape((-1,1)))).writeto("optextract4.fits", clobber=True)

    #
    # Now compute the optimally extracted spectrum
    #
    weighted_spec = cutout_sel * profile.reshape((-1,1))
    clobberfile("optextract5.fits")
    pyfits.PrimaryHDU(data=weighted_spec, header=hdr).writeto("optextract5.fits", clobber=True)
    source_spec = bottleneck.nansum(weighted_spec, axis=0) / bottleneck.nansum(profile.reshape((-1,1))[profile_sel])

    #
    # So far the source spectrum is the weighted average across the aperture
    # multiply it with the pixel width to get a sum
    #
    _weighted = numpy.ones_like(weighted_spec)
    _weighted[~numpy.isfinite(weighted_spec)] = 0.
    _mean_to_sum = bottleneck.nansum(_weighted, axis=0)
    source_spec *= _mean_to_sum

    # print source_spec.shape
    numpy.savetxt("optspec", source_spec)
    numpy.savetxt("noptspec", bottleneck.nanmean(cutout_sel, axis=0))



    # # output as text file
    # outfile = fitsfile[:-5]+".spec"
    # with open(outfile, "w") as of:
    #     numpy.savetxt(of, comb)
    #     print >>of, "\n\n\n\n\n"

    # sys.exit(0)

    spec1d = source_spec

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
        # print "masking"
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
    # print comb.shape

    # output as text file
    outfile = fitsfile[:-5]+".optspec"
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
    clobberfile(fitsfile[:-5]+".opt1d.fits")
    primhdu.writeto(fitsfile[:-5]+".opt1d.fits", clobber=True)

    
