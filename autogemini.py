#!/usr/bin/env python

import os, sys
import pyfits

import pyraf
from pyraf import iraf

sys.stdout.write("Loading GEMINI IRAF package, please be patient ...")
sys.stdout.flush()
from pyraf.iraf import gemini
from pyraf.iraf import gmos
print "done!"


def clobberfile(fn):
    if (os.path.isfile(fn)):
        os.remove(fn)

if __name__ == "__main__":

    filelist = sys.argv[1:]

    bias_list = []
    flat_list = {}
    arc_list = []
    object_list = []
    std_list = []

    for filename in filelist:

        #print filename
        print

        hdu = pyfits.open(filename)
        obstype = hdu[0].header["OBSTYPE"]
        print filename, obstype

        if (obstype == "BIAS"):
            bias_list.append(filename)

        elif (obstype == "FLAT"):
            grtilt = hdu[0].header['GRTILT']
            if (grtilt not in flat_list):
                flat_list[grtilt] = []
            
            flat_list[grtilt].append(filename)

        elif (obstype == "ARC"):
            arc_list.append(filename)

        elif (obstype == "OBJECT"):
            if (hdu[0].header['OBSCLASS'] == "acq"):
                print "Skipping acquisition image",filename
                continue
            object_list.append(filename)


    print "done reading all files"

    print flat_list

    master_bias = bias_list[0]

    #
    # Create all flat-fields for each of the grtilt angles
    #
    for grtilt in flat_list:

        flats = flat_list[grtilt]
        flat_out = "masterflat__%.2f.fits" % (grtilt)

        #
        # Delete all the intermediate files created by gsflat/gsreduce
        # Just to be sure we don;t run into trouble, delete them if they 
        # happen to exist
        #
        for ffn in flats:
            for prefix in ['g', 'gs']:
                _, bn = os.path.split(ffn)
                tmpfile = "%s%s" % (prefix, bn)
                clobberfile(tmpfile)

        print "computing masterflat:", flat_out
        clobberfile(flat_out)
        iraf.gemini.gmos.gsflat(
            inflats=",".join(flats),
            specflat=flat_out,
            order=23,
            bias=master_bias,
            fl_over=False,
            )

        #
        # Delete all the intermediate files created by gsflat/gsreduce to 
        # clean things up in case we need to run it again
        #
        for ffn in flats:
            for prefix in ['g', 'gs']:
                _, bn = os.path.split(ffn)
                tmpfile = "%s%s" % (prefix, bn)
                clobberfile(tmpfile)
        


    #
    # Now run the actual science frames, ARCs first
    #
    arc_specs = {}
    arc_mjd = {}

    for arcfile in arc_list:
        print "ARC:", arcfile
        _, bn = os.path.split(arcfile)
        for prefix in ['g', 'gs']:
            tmpfile = "%s%s" % (prefix, bn)
            clobberfile(tmpfile)

        #
        # get some basic specs about the arc
        #
        hdu = pyfits.open(arcfile)
        grating = hdu[0].header['GRATING']
        grtilt = hdu[0].header['GRTILT']
        mjdobs = 0.0 #hdu[0].header['MJD-OBS']

        arc_mjd[bn] = mjdobs
        if (not grating in arc_specs):
            arc_specs[grating] = {}
        if (not grtilt in arc_specs[grating]):
            arc_specs[grating][grtilt] = []
        arc_specs[grating][grtilt].append(bn)

        arc_reduced = "arc__"+bn
        trans_arc = "trans__"+bn

        if (os.path.isfile(trans_arc)):
            continue

        #
        # Reduce ARC spectrum
        #
        clobberfile(arc_reduced)
        iraf.gemini.gmos.gsreduce(
            inimages=arcfile,
            outimages=arc_reduced,
            outpref="gs",
            fl_flat=False,
            bias=master_bias,
            fl_over=False,
            fl_fixpix=False
            )

        #
        # Find wavelength solution
        #
        iraf.gemini.gmos.gswavelength(
            inimages=arc_reduced,
            fl_inter=False,  # run non-interactively
            )

        #
        # Transform ARC spectrum just for checking
        #
        clobberfile(trans_arc)
        iraf.gemini.gmos.gstransform(
            inimages=arc_reduced,
            wavtran=arc_reduced[:-5],
            outimages=trans_arc,
            )


    #
    # Now reduce the actual science spectra
    # 
    for obj_file in object_list:
        _, bn = os.path.split(obj_file)

        reduced = bn[:-5]+".red.fits"
        trans = bn[:-5]+".trans.fits"
        skysub = bn[:-5]+".skysub.fits"

        clobberfile(reduced)
        clobberfile(trans)
        clobberfile(skysub)

        iraf.gemini.gmos.gsreduce(
            inimages=obj_file,
            outimages=reduced,
            outpref="gs",
            fl_flat=False, # Apply flat field correction
            bias=master_bias,
            fl_over=False, # Subtract overscan level (done via BIAS)
            fl_fixpix=False, # Interpolate across chip gaps if mosaicing
            #
            fl_gscr=True, # Clean images for cosmic rays
            fl_gmos=True, # Mosaic science extensions
            fl_vard=True, # Create variance and data quality frames
            )
        
        #
        # from all arc spectra, find the one closest in time with the
        # right combination of grating and grating angle
        #
        hdu = pyfits.open(obj_file)
        grating = hdu[0].header['GRATING']
        grtilt = hdu[0].header['GRTILT']
        mjdobs = 0.0 #hdu[0].header['MJD-OBS']

        if (not grating in arc_specs):
            print "No ARC found for this grating: %s" % (grating)
            continue
        if (not grtilt in arc_specs[grating]):
            print "No ARC found for this grating angle: %.2f" % (grtilt)
            continue

        good_arcs = arc_specs[grating][grtilt]

        # Now find the one closest in MJD to the observation
        best_arc = good_arcs[0]
        arc_db_name = "arc__"+(best_arc[:-5] if best_arc.endswith(".fits") else best_arc)

        iraf.gemini.gmos.gstransform(
            inimages=reduced,
            wavtran=arc_db_name,
            outimages=trans,
            )

        # 594:661,1000:1130
        #
        # Select region to be used as sky
        #
        skyfile = "skymask_%s.txt" % (bn[:-5])
        sky_string = None
        if (os.path.isfile(skyfile)):
            # load sky regions
            with open(skyfile, "r") as sf:
                sky_string = sf.readlines()[0].strip()
        else:
            # ask user for sky regions
            while (True):
                usr = raw_input("sky region y1:y2 (ENTER when done) : ")
                if (usr == ""):
                    break
                sky_string = usr if sky_string == None else "%s,%s" % (sky_string, usr) 
            with open(skyfile, "w") as sf:
                print >>sf, sky_string

        if (not sky_string == None):
            iraf.gemini.gmos.gsskysub(
                input=trans,
                fl_answer=False,
                output=skysub,
                fl_inter=False,
                long_sample=sky_string,
                )
                
                
