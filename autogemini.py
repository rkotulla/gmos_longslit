#!/usr/bin/env python

import os, sys
import pyfits

import argparse

import pyraf
from pyraf import iraf



def clobberfile(fn):
    if (os.path.isfile(fn)):
        os.remove(fn)

if __name__ == "__main__":
    
    #
    # Handle command line options
    #

    parser = argparse.ArgumentParser(description='Some arguments')
    parser.add_argument('--myskysub', dest='myskysub', action='store_true',
                        default=False,
                        help='Use the custom sky-subtraction method')
    parser.add_argument('--redo', dest='redo', action='store_true',
                        default=False,
                        help='Recreate and overwrite files that already exist')
    parser.add_argument('files', nargs='+', metavar='files',
                        help='raw files to be processed')

    args = parser.parse_args()

    #print "redo:", args.redo

    #
    # Do work
    #

    filelist = args.files

    bias_list = []
    flat_list = {}
    arc_list = []
    object_list = []
    std_list = []

    for filename in filelist:

        #print filename
        #print

        hdu = pyfits.open(filename)
        obstype = hdu[0].header["OBSTYPE"]
        target = hdu[0].header['OBJECT'] if 'OBJECT' in hdu[0].header else "???"

        #print filename, obstype

        if (obstype == "BIAS"):
            print "%s: Found BIAS frame" % (filename)
            bias_list.append(filename)

        elif (obstype == "FLAT"):
            grtilt = hdu[0].header['GRTILT']
            grating = hdu[0].header['GRATING']
            ccdsum = "x".join(hdu[1].header['CCDSUM'].split())
            naxis2 = hdu[1].header['NAXIS2']

            if (grating not in flat_list):
                flat_list[grating] = {}
            if (grtilt not in flat_list[grating]):
                flat_list[grating][grtilt] = {}
            if (ccdsum not in flat_list[grating][grtilt]):
                flat_list[grating][grtilt][ccdsum] = {}
            if (naxis2 not in flat_list[grating][grtilt][ccdsum]):
                flat_list[grating][grtilt][ccdsum][naxis2] = []
            
            flat_list[grating][grtilt][ccdsum][naxis2].append(filename)
            print "%s: Found FLAT" % (filename)

        elif (obstype == "ARC"):
            arc_list.append(filename)
            print "%s: Found ARC" % (filename)

        elif (obstype == "OBJECT"):
            obsclass = hdu[0].header['OBSCLASS']
            if (obsclass in ["acq", "acqCal"]):
                print "%s: Skipping acquisition image" % (filename)
                continue
            elif (obsclass == 'partnerCal'):
                print "%s: Found Standard star %s" % (filename, target)
                std_list.append(filename)
            elif (obsclass == "dayCal"):
                print "%s: Skipping daytime calibration (%s)" % (filename, target)
                continue
            else:
                print "%s: Found science target %s" % (filename, target)
                object_list.append(filename)


    print "done reading all files"

    sys.stdout.write("Loading GEMINI IRAF package, please be patient ...")
    sys.stdout.flush()
    from pyraf.iraf import gemini
    from pyraf.iraf import gmos
    print "done!"

    print flat_list

    master_bias = bias_list[0]

    #
    # Create all flat-fields for each of the grtilt angles
    #
    for grating in flat_list:
        for grtilt in flat_list[grating]:
            for ccdsum in flat_list[grating][grtilt]:
                for naxis2 in flat_list[grating][grtilt][ccdsum]:
                    
                    flats = flat_list[grating][grtilt][ccdsum][naxis2]
        
                    flat_out = "masterflat__%s__%.4f__%s__%04d.fits" % (
                        grating, grtilt, ccdsum, naxis2)

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

                    if (not os.path.isfile(flat_out) or args.redo):
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

        arc_reduced = "arc__"+bn
        trans_arc = "trans__"+bn

        #
        # get some basic specs about the arc
        #
        hdu = pyfits.open(arcfile)
        grating = hdu[0].header['GRATING']
        grtilt = hdu[0].header['GRTILT']
        ccdsum = "x".join(hdu[1].header['CCDSUM'].split())
        naxis2 = hdu[1].header['NAXIS2']
        mjdobs = 0.0 #hdu[0].header['MJD-OBS']

        arc_mjd[bn] = mjdobs
        if (not grating in arc_specs):
            arc_specs[grating] = {}
        if (not grtilt in arc_specs[grating]):
            arc_specs[grating][grtilt] = {}
        if (not ccdsum in arc_specs[grating][grtilt]):
            arc_specs[grating][grtilt][ccdsum] = {}
        if (not naxis2 in arc_specs[grating][grtilt][ccdsum]):
            arc_specs[grating][grtilt][ccdsum][naxis2] = []
            
        arc_specs[grating][grtilt][ccdsum][naxis2].append(bn)

        if (os.path.isfile(trans_arc) and not args.redo):
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
        try:
            iraf.gemini.gmos.gswavelength(
                inimages=arc_reduced,
                fl_inter=False,  # run non-interactively
            )
        except:
            iraf.gemini.gmos.gswavelength(
                inimages=arc_reduced,
                fl_inter='NO',  # run non-interactively
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
    for obj_file in object_list+std_list:
        _, bn = os.path.split(obj_file)

        reduced = bn[:-5]+".red.fits"
        trans = bn[:-5]+".trans.fits"
        skysub = bn[:-5]+".skysub.fits"

        print "\n\nWorking on %s\n\n" % (bn)

        # Find the correct flat-field to be used
        obj_hdu = pyfits.open(obj_file)
        grating = obj_hdu[0].header['GRATING']
        grtilt = obj_hdu[0].header['GRTILT']
        ccdsum = "x".join(obj_hdu[1].header['CCDSUM'].split())
        naxis2 = obj_hdu[1].header['NAXIS2']
        ff_name = "masterflat__%s__%.4f__%s__%04d.fits" % (
            grating, grtilt, ccdsum, naxis2)
        

        #
        # Reduce the spectrum (apply bias, flat, etc)
        #
        if (not os.path.isfile(reduced) or args.redo):
            clobberfile(reduced)


            
            print "Using %s for %s" % (ff_name, obj_file)
            print grating, grtilt, ccdsum, naxis2
            #continue
        
            iraf.gemini.gmos.gsreduce(
                inimages=obj_file,
                outimages=reduced,
                outpref="gs",
                fl_flat=os.path.isfile(ff_name), 
                # Apply flat field correction if the right file exists
                flatim=ff_name,
                bias=master_bias,
                fl_over=False, # Subtract overscan level (done via BIAS)
                fl_fixpix=False, # Interpolate across chip gaps if mosaicing
                #
                fl_gscr=False, #True, # Clean images for cosmic rays
                fl_gmos=True, # Mosaic science extensions
                fl_vard=True, # Create variance and data quality frames
                )
        

        #
        # Rectify the reduced spectrum
        #

        if (not os.path.isfile(trans) or args.redo):
            clobberfile(trans)

            #
            # from all arc spectra, find the one closest in time with the
            # right combination of grating and grating angle
            #
            #hdu = pyfits.open(obj_file)
            #grating = hdu[0].header['GRATING']
            #grtilt = hdu[0].header['GRTILT']
            mjdobs = 0.0 #hdu[0].header['MJD-OBS']

            
            if (not grating in arc_specs):
                print "No ARC found for this grating: %s" % (grating)
                continue
            elif (not grtilt in arc_specs[grating]):
                print "No ARC found for this grating angle: %.2f" % (grtilt)
                continue
            elif (not ccdsum in arc_specs[grating][grtilt]):
                print "No ARC found for this binning: %s" % (ccdsum)
                continue
            elif (not naxis2 in arc_specs[grating][grtilt][ccdsum]):
                print "No ARC found for this y-dimension: %s" % (naxis2)
                continue

            good_arcs = arc_specs[grating][grtilt][ccdsum][naxis2]

            # Now find the one closest in MJD to the observation
            best_arc = good_arcs[0]
            arc_db_name = "arc__"+(best_arc[:-5] if best_arc.endswith(".fits") else best_arc)

            iraf.gemini.gmos.gstransform(
                inimages=reduced,
                wavtran=arc_db_name,
                outimages=trans,
                )


        #
        # Apply sky-subtraction
        #

        if (not os.path.isfile(skysub) or args.redo):
            clobberfile(skysub)
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

                

    #
    # Now all frames are sky-subtracted, let's work on the standard star to 
    # find the flux-calibration
    #
    for std in std_list:
        print std
