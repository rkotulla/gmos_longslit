#!/usr/bin/env python

import os, sys
import pyfits
import math

import argparse

import pyraf
from pyraf import iraf

import skysubtract

def clobberfile(fn):
    if (os.path.isfile(fn)):
        os.remove(fn)

import datetime

def load_fits_header_edit(fn):

    header_edit = {}

    with open(fn, "r") as fhe:
        lines = fhe.readlines()
        for line in lines:
            items = line.split()
            
            fn = items[0]
            fitskey = items[1]
            # value = items[2]

            try:
                value = float(items[2])
            except:
                value = items[2]

            _, filebase = os.path.split(fn)
           
            if (not filebase in header_edit):
                header_edit[filebase] = []
            header_edit[filebase].append((fitskey, value))
            
    return header_edit

def apply_fits_header_edit(fn, hdulist, fhe):

    _, base = os.path.split(fn)
    if (base in fhe):
        for ext in hdulist:
            for (key, value) in fhe[base]:
                if (key in ext.header):
                    ext.header[key] = value

if __name__ == "__main__":
    
    #
    # Handle command line options
    #

    parser = argparse.ArgumentParser(description='Some arguments')
    parser.add_argument('--myskysub', dest='myskysub', action='store_true',
                        default=False,
                        help='Use the custom sky-subtraction method')
    # parser.add_argument('--redo', dest='redo', action='store_true',
    #                     default=False,
    #                     help='Recreate and overwrite files that already exist')

    parser.add_argument('--redo', dest='redo', action='store',
                        default='',
                        help='Recreate and overwrite files that already exist')

    parser.add_argument('--fitsedit', dest='fitsedit',
                        action='store', default='',
                        help='filename with FITS headers to edit on-the-fly')
    parser.add_argument('--noarcs', dest='noarcs',
                        action='store_true', default=False,
                        help='skip processing all ARC frames & re-use existing database')
    parser.add_argument('files', nargs='+', metavar='files',
                        help='raw files to be processed')

    args = parser.parse_args()

    fits_header_edit = {}
    if (os.path.isfile(args.fitsedit)):
        fits_header_edit = load_fits_header_edit(args.fitsedit)
        
    print fits_header_edit

    redo = None
    if (args.redo != ''):
        redo = [s.upper() for s in args.redo.split(",")]

    command_list = open("commands.cl", "w")

    # os._exit(0)

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

    timestamps = {}
    arc_specs_info = {}

    for filename in filelist:

        #print filename
        #print
        _, bn = os.path.split(filename)

        hdu = pyfits.open(filename)
        apply_fits_header_edit(filename, hdu, fits_header_edit)

        obstype = hdu[0].header["OBSTYPE"]
        target = hdu[0].header['OBJECT'] if 'OBJECT' in hdu[0].header else "???"

        date_time = ("%sT%s" % (hdu[0].header['DATE'], hdu[0].header['UT']))[:19]
        ts = (datetime.datetime.strptime(date_time, "%Y-%m-%dT%H:%M:%S") - datetime.datetime.min).total_seconds()
        timestamps[bn] = ts
        
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

            _arcspec = {'grating': grating,
                        'grtilt': grtilt,
                        'ccdsum': ccdsum,
                        'naxis2': naxis2,
                        }
            _, bn = os.path.split(filename)
            arc_specs_info[bn] = _arcspec

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
    sys.stdout.write(" done!\n")
    sys.stdout.flush()
    # print "done!"

    print flat_list

    master_bias = bias_list[0]

    #
    # Create all flat-fields for each of the grtilt angles
    #
    
    print "Creating flat-fields!"
    for grating in flat_list:
        for grtilt in flat_list[grating]:
            for ccdsum in flat_list[grating][grtilt]:
                for naxis2 in flat_list[grating][grtilt][ccdsum]:
                    
                    flats = flat_list[grating][grtilt][ccdsum][naxis2]
        
                    flat_out = "masterflat__%s__%.4f__%s__%04d.fits" % (
                        grating, grtilt, ccdsum, naxis2)

                    if (os.path.isfile(flat_out) and (not redo == [] or 'FLAT' in redo)):
                        print("Skipping generation of flat %s" % (flat_out))
                        continue

                    
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
                        _stdout = iraf.gemini.gmos.gsflat(
                            inflats=",".join(flats),
                            specflat=flat_out,
                            order=23,
                            bias=master_bias,
                            fl_over=False,
                            Stdout=1
                            )

                        print >>command_list, """
                        iraf.gemini.gmos.gsflat(
                            inflats=%s,
                            specflat=flat_out,
                            order=23,
                            bias=master_bias,
                            fl_over=False,
                            )""" % (",".join(flats))


                    #
                    # Delete all the intermediate files created by gsflat/gsreduce to 
                    # clean things up in case we need to run it again
                    #
                    for ffn in flats:
                        for prefix in ['g', 'gs']:
                            _, bn = os.path.split(ffn)
                            tmpfile = "%s%s" % (prefix, bn)
                            clobberfile(tmpfile)
        
    print ("done with flat-field generation!")

    #
    # Now run the actual science frames, ARCs first
    #
    arc_specs = {}
    arc_mjd = {}
    print("Working on ARCs next...")
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
        apply_fits_header_edit(arcfile, hdu, fits_header_edit)

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

        if (os.path.isfile(trans_arc) and 
            (not redo == [] or 'ARC' in redo or args.noarcs)):
            print "Skipping processing of ARC %s" % (arcfile)
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

        print >>command_list, """
        iraf.gemini.gmos.gsreduce(
            inimages=%s,
            outimages=%s,
            outpref="gs",
            fl_flat=False,
            bias=%s,
            fl_over='no',
            fl_fixpix='yes'
            )""" % (arcfile, arc_reduced, master_bias)

        #
        # Find wavelength solution
        #
        try:
            iraf.gemini.gmos.gswavelength(
                inimages=arc_reduced,
                fl_inter=True,  # run non-interactively
            )

        except:
            iraf.gemini.gmos.gswavelength(
                inimages=arc_reduced,
                fl_inter='NO',  # run non-interactively
            )

        print >>command_list, """
            iraf.gemini.gmos.gswavelength(
                inimages=%s,
                fl_inter='yes',  # run non-interactively
            )""" % (arc_reduced)


        # #
        # # Transform ARC spectrum just for checking
        # #
        # clobberfile(trans_arc)
        # iraf.gemini.gmos.gstransform(
        #     inimages=arc_reduced,
        #     wavtran=arc_reduced[:-5],
        #     outimages=trans_arc,
        #     )

    print arc_specs_info

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
        apply_fits_header_edit(obj_file, obj_hdu, fits_header_edit)

        grating = obj_hdu[0].header['GRATING']
        grtilt = obj_hdu[0].header['GRTILT']
        ccdsum = "x".join(obj_hdu[1].header['CCDSUM'].split())
        naxis2 = obj_hdu[1].header['NAXIS2']
        ff_name = "masterflat__%s__%.5f__%s__%04d.fits" % (
            grating, grtilt, ccdsum, naxis2)
        

        #
        # Reduce the spectrum (apply bias, flat, etc)
        #
        if (not os.path.isfile(reduced) or (redo != None and 'OBJ' in redo)):
            clobberfile(reduced)

            # print "Using %s for %s" % (ff_name, obj_file)
            # print grating, grtilt, ccdsum, naxis2
        
            # Don't fix cosmics in STD exposures as this seems to 
            # screw things up
            crj_fix = obj_file in object_list

            print("Reducing %s\n   bias: %s\n   flat: %s" % (obj_file, master_bias, ff_name))
            # Check if the bias dimensions match the frame dimensions
            _stdout = iraf.gemini.gmos.gsreduce(
                inimages=obj_file,
                outimages=reduced,
                outpref="gs",
                fl_flat=os.path.isfile(ff_name), 
                # Apply flat field correction if the right file exists
                flatim=ff_name,
                bias=master_bias,
                fl_over=True, #False, # Subtract overscan level (done via BIAS)
                fl_fixpix=False, # Interpolate across chip gaps if mosaicing
                #
                fl_gscr=crj_fix, #True, # Clean images for cosmic rays
                fl_gmos=True, # Mosaic science extensions
                fl_vard=True, # Create variance and data quality frames
                Stdout=1,
                )
            print "\n".join(_stdout)
        

        #
        # Rectify the reduced spectrum
        #

        if (not os.path.isfile(trans) or (redo != None and 'WLCAL' in redo)): #args.redo or True):
            clobberfile(trans)

            #
            # from all arc spectra, find the one closest in time with the
            # right combination of grating and grating angle
            #
            #hdu = pyfits.open(obj_file)
            #grating = hdu[0].header['GRATING']
            #grtilt = hdu[0].header['GRTILT']
            mjdobs = 0.0 #hdu[0].header['MJD-OBS']

            print "ARCSEL   Working on %s" % (obj_file)
            print "ARCSEL   grating: %s  //  GRTILT: %8.5f  //  NAXIS: %d  //  BIN: %s" % (
                grating, grtilt, naxis2, ccdsum)
            if (not grating in arc_specs):
                print "ARCSEL   No ARC found for this grating: %s" % (grating)
                continue
            elif (not grtilt in arc_specs[grating]):
                print "No ARC found for this grating angle: %.2f" % (grtilt)
                continue
            elif (not ccdsum in arc_specs[grating][grtilt]):
                print "ARCSEL   No ARC found for this binning: %s" % (ccdsum)
                continue
            elif (not naxis2 in arc_specs[grating][grtilt][ccdsum]):
                print "ARCSEL   No ARC found for this y-dimension: %s" % (naxis2)
                continue

            good_arcs = []
            for ccdsum in arc_specs[grating][grtilt]:
                for naxis2 in arc_specs[grating][grtilt][ccdsum]:
                    good_arcs.extend(arc_specs[grating][grtilt][ccdsum][naxis2])
            #good_arcs = arc_specs[grating][grtilt][ccdsum][naxis2]
            print "ARCSEL   Found these arcs for frame %s:\nARCSEL   -- %s" % (
                obj_file,
                "\nARCSEL   -- ".join(['%s: %s' % (
                    fn, "%(grating)s / %(grtilt).4f / %(ccdsum)s / %(naxis2)d" % arc_specs_info[fn]) for fn in good_arcs]),
#                "\nARCSEL   -- ".join(good_arcs),
            )

            # Now find the one closest in MJD to the observation
            min_delta_t = 1e99
            for arc in good_arcs:
                dt = math.fabs(timestamps[bn] - timestamps[arc])
                if (dt < min_delta_t):
                    min_delta_t = dt
                    best_arc = arc
            print "ARCSEL   +++>> %s" % (best_arc)

            #best_arc = good_arcs[0]
            arc_db_name = "arc__"+(best_arc[:-5] if best_arc.endswith(".fits") else best_arc)
            # continue

            print("Running GSTransform (%s --> %s)" % (reduced, trans))
            _stdout = iraf.gemini.gmos.gstransform(
                inimages=reduced,
                wavtran=arc_db_name,
                outimages=trans,
                Stdout=1
                )
            #print _stdout


        #
        # Apply sky-subtraction
        #

        if (not os.path.isfile(skysub) or (redo != None and 'SKYSUB' in redo)): #args.redo):
            clobberfile(skysub)
            # 594:661,1000:1130
            #
            # Select region to be used as sky
            #
            print "Subtracting sky background"
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

            if (args.myskysub):
                skysubtract.subsky(
                    fitsfile=trans, 
                    outputfile=skysub)
            else:
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
    print "list of standard stars:", std_list


    for std_file in std_list:
        print std_file

        _, bn = os.path.split(std_file)

        reduced = bn[:-5]+".red.fits"
        trans = bn[:-5]+".trans.fits"
        skysub = bn[:-5]+".skysub.fits"

        spec1d = bn[:-5]+".1d.fits"
        fluxfile = bn[:-5]+".flux.fits"
        sensfile = bn[:-5]+".sens.fits"

        std_hdu = pyfits.open(std_file)
        starname = std_hdu[0].header['OBJECT']

        if (not os.path.isfile(spec1d) or (redo != None and "STD" in redo)):
            clobberfile(spec1d)
            clobberfile(fluxfile)
            clobberfile(sensfile)

            #
            # Extract the spectrum
            #
            print("Extracting 1-d spectrum for standard star")
            iraf.gemini.gmos.gsextract(
                inimages=skysub,
                outimages=spec1d,
                apwidth=1.0,                    # 1 arcsec
                fl_inter=False,                 # Run interactively?
                find=True,                      # Define apertures automatically?
                recenter=True,                  # Recenter apertures?
                trace=True,                     # Trace apertures?
                tfunction = "chebyshev",        # Trace fitting function
                torder = 5,                     # Trace fitting function order
                tnsum = 20,                     # Number of dispersion lines to sum for trace
                tstep = 50,                     # Tracing step
                weights = "none",               # Extraction weights (none|variance)
                clean = False,                  # Detect and replace bad pixels?
                lsigma = 3.,                    # Lower rejection threshold for cleaning
                usigma = 3.,                    # Upper rejection threshold for cleaning
            )

            #
            # Establish the sensitivity function
            #
            print("Establish sensitivity response function")
            clobberfile("std")
            clobberfile("sens.fits")
            iraf.gemini.gmos.gsstandard(
                input = spec1d,
                sfile = "std",                  # Output flux file (used by SENSFUNC)
                sfunction = "sens.fits",             # Output root sensitivity function image name
                sci_ext = "SCI",                # Name or number of science extension
                var_ext = "VAR",                # Name or number of variance extension
                dq_ext = "DQ",                  # Name or number of data quality extension
                key_airmass = "AIRMASS",        # Header keyword for airmass
                key_exptime = "EXPTIME",        # Header keyword for exposure time
                fl_inter = False,               # Run the task interactively
                starname = starname,            # Standard star name(s) in calibration list
                samestar = True,                # Same star in all apertures
                apertures = "",                 # Aperture selection list
                beamswitch = False,             # Beam switch spectra
                bandwidth = 'INDEF',            # Bandpass width
                bandsep = 'INDEF',              # Bandpass separation
                fnuzero = 3.6800000000000E-20,  # Absolute flux zero point
                caldir = "onedstds$spec50cal/", # Directory containing calibration data
                observatory = "Gemini-North",   # Observatory
                mag = "",                       # Magnitude of stars
                magband = "",                   # Magnitude types/bands (U|B|V|R|I|J|H|K|L|Lprime
                teff = "",                      # Effective temperature of spectral types
                ignoreaps = True,               # Ignore apertures and make one sensitivity funct
                extinction = "",                # Extinction file
                out_extincti = "extinct.dat",   # Output revised extinction file
                function = "spline3",           # Fitting function
            )

#    sys.exit(0)

    #
    # Now apply flux-calibration to all object frames
    # 
    for obj_file in object_list+std_list:

        #
        # Apply flux-calibration to spectra
        # 
        # Note: We use a flux-calibration factor of 1e15, so all resulting
        #       fluxes are in units of 1e-15 erg/s/A/cm^2
        #

        _, bn = os.path.split(obj_file)
        skysub = bn[:-5]+".skysub.fits"
        fluxcal = bn[:-5]+".fluxcal.fits"

        if (not os.path.isfile(fluxcal) or (redo != None and "FLUXCAL" in redo)):
            clobberfile(fluxcal)

            _stdout=iraf.gemini.gmos.gscalibrate(
                input = skysub,               # Input spectra to calibrate
                output = fluxcal,             # Output calibrated spectra
                sfunction = "sens",           # Input image root name for sensitivity function
                sci_ext = "SCI",              # Name of science extension
                var_ext = "VAR",              # Name of variance extension
                dq_ext = "DQ",                # Name of data quality extension
                key_airmass = "AIRMASS",      # Airmass header keyword
                key_exptime = "EXPTIME",      # Exposure time header keyword
                fl_vardq = False,              # Propagate VAR/DQ planes
                fl_ext = False,               # Apply extinction correction to input spectra
                fl_flux = True,               # Apply flux calibration to input spectra
                fl_scale = True,              # Multiply output with fluxscale
                fluxscale = 1e15,             # Value of the flux scale (fl_scale=yes)
                ignoreaps = True,             # Ignore aperture numbers in flux calibration
                fl_fnu = False,               # Create spectra having units of FNU
                extinction = "",              # Extinction file
                observatory = "Gemini-North", # Observatory
                verbose = True,               # Verbose?
                Stdout=1,
                )
