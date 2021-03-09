#!/usr/bin/python3
#
# scripMaker.py v0.1
#
# This script takes as input a netCDF4 grid and a netCDF4 mask
# to produce an output netCDF4 file with the grid converted to
# the SCRIP format.
#
# Developed by Fabio Viola adapting Pier Giuseppe Fogli's Matlab
# code to our specific use case. For debug purposes, note that
# at the end of many lines is reported the line number of the
# original Matlab equivalent instruction(s).
#
# Released under GPL v3.0 license
# Contact info: Fabio Viola <fabio.viola@cmcc.it>


# global reqs
import pdb
import sys
import numpy
import getopt
import logging
import datetime
import traceback
from netCDF4 import Dataset


# main
if __name__ == "__main__":

    ########################################################
    #
    # create a logger
    #
    ########################################################
    
    logging.basicConfig()
    logger = logging.getLogger()
    logger.setLevel(level=logging.DEBUG)

    
    ########################################################
    #
    # read input params
    #
    ########################################################

    logger.debug("Reading input params...")
    gridFile = None
    maskFile = None
    outputFile = None
    options, remainder = getopt.getopt(sys.argv[1:], 'o:g:m:', ['output=', 'grid=', 'mask='])
    for opt, arg in options:
        if opt in ('-o', '--output'):
            outputFile = arg
            logger.debug("Output file set to: %s" % arg)
        elif opt in ('-g', '--grid'):
            gridFile = arg
            logger.debug("Grid file set to: %s" % arg)
        elif opt in ('-m', '--mask'):
            maskFile = arg
            logger.debug("Mask file set to: %s" % arg)
        else:
            logger.warning("Option %s not recognised" % opt)

    if not(gridFile) or not(outputFile):
        logger.error("You must pass the two input files and an output file name")
        sys.exit(1)

        
    ########################################################
    #
    # open input files
    #
    ########################################################
    
    logger.debug("Opening input file %s" % gridFile)
    try:
        iFile1 = Dataset(gridFile)
    except FileNotFoundError:
        logger.error("File not found %s" % gridFile)
        sys.exit(2)

    logger.debug("Opening input file %s" % maskFile)
    if maskFile:
        try:
            iFile2 = Dataset(maskFile)
        except FileNotFoundError:
            logger.error("File not found %s" % maskFile)
            sys.exit(3)

        
    ########################################################
    #
    # processing...
    #
    ########################################################

    ### Processing starts at line ~80
    
    # read T coordinates
    lon = numpy.squeeze(iFile1.variables["glamt"]).transpose() ### 305, 319
    lat = numpy.squeeze(iFile1.variables["gphit"]).transpose() ### 306, 320
    nlon = lon.shape[0] ### 327
    nlat = lat.shape[1] ### 328
    grid_size = nlat * nlon ### 456
    grid_dims = [nlon, nlat] ### 457
    grid_corners = 4 ### 458
    grid_rank = 2 ### 459

    # read F coordinates
    clon = numpy.squeeze(iFile1.variables["glamf"]).transpose() ### 310, 321
    clat = numpy.squeeze(iFile1.variables["gphif"]).transpose() ### 311, 322

    # load T grid metrics
    e1t = numpy.squeeze(iFile1.variables["e1t"]).transpose() ### 315, 323
    e2t = numpy.squeeze(iFile1.variables["e2t"]).transpose() ### 316, 324

    # load the land-sea mask (0=land, 1=sea)
    if maskFile:
        tmask = numpy.squeeze(iFile2.variables["tmask"]) ### 334
        if len(tmask.shape) == 3: ### 335
            tmask = numpy.squeeze(tmask[0,:,:]) ### 336
            
    mask = numpy.ones(lon.shape) ### 339
    if maskFile:
        mask = numpy.where(tmask>0.0, 1.0, 0.0) ### 340
    mask = mask.transpose() ### 446
    
    # compute T grid corners coordinates (F points)
    tarea = e1t * e2t  ### 468
    erad = 6371229; # NEMO Earth radius (m) ### 479
    tarea = tarea/(erad^2) # Normalise (probably useless) ### 487    
    
    # Translate longitude from range [-180:180] to [0:360]
    # NOTE: this is probably useless for our use case. Keep for future?
    lon = numpy.where(lon<0, lon+360, lon) ### 451
    clon = numpy.where(clon<0, clon+360, clon) ### 452

    # convert to radiants
    rlat = numpy.radians(lat) ### 490
    rlon = numpy.radians(lon) ### 491
    rclat = numpy.radians(clat) ### 492
    rclon = numpy.radians(clon) ### 493
    
    # determine fxt
    fxt = numpy.radians(numpy.floor(numpy.min(lat[0,:])))

    # initialise nemo_clo and nemo_cla
    nemo_clo = numpy.zeros((nlon, nlat, grid_corners)) ### 496
    nemo_cla = numpy.zeros((nlon, nlat, grid_corners)) ### 497

    # Fill the corners arrays: anti-clockwise starting from upper-right

    # Corner 1 lon and lat
    nemo_clo[:,:,0] = rclon; ### 500
    nemo_cla[:,:,0] = rclat; ### 513
    
    # Corner 2 lon and lat
    nemo_clo[0,:,1] = rclon[-1,:]    ### 503
    nemo_clo[1:,:,1] = rclon[0:-1,:] ### 502
    nemo_cla[0,:,1] = rclat[-3,:]    ### 516
    nemo_cla[1:,:,1] = rclat[0:-1,:] ### 515

    # Corner 3 lon and lat
    nemo_clo[1:,1:,2] = rclon[0:-1,0:-1]  ### 505
    nemo_clo[0,1:,2] = rclon[-1, 0:-1]    ### 506
    nemo_clo[1:,0,2] = rclon[0:-1,0]      ### 507
    nemo_clo[0,0,2] = rclon[-1,0]         ### 508
    nemo_cla[1:,1:,2] = rclat[0:-1,0:-1]  ### 518
    nemo_cla[0,1:,2] = rclat[-3,0:-1]     ### 519
    nemo_cla[:,0,2] = fxt                 ### 520

    # Corner 4 lon and lat
    nemo_clo[:,1:,3] = rclon[:,0:-1]    ### 510
    nemo_clo[:,0,3] = rclon[:,0]        ### 511
    nemo_cla[:,1:,3] = rclat[:,0:-1]    ### 522
    nemo_cla[:,0,3] = fxt               ### 523
        
    # Reshape to SCRIP convention
    rlon = numpy.matrix.flatten(rlon)   ### 693
    rlat = numpy.matrix.flatten(rlat)   ### 694
    mask = numpy.matrix.flatten(mask)   ### 695
    tarea = numpy.matrix.flatten(tarea) ### 696

    # reshape nemo_clo
    nemo_clo_final = nemo_clo.reshape((grid_size, grid_corners))
    
    # reshape nemo_cla
    nemo_cla_final = nemo_cla.reshape((grid_size, grid_corners))
    
    ### Processing ends at line 698
    ### After that, it's just production of the output NetCDF files
    
    
    ########################################################
    #
    # open and write output file
    #
    ########################################################
        
    logger.debug("Opening output file...")
    try:
        oFile = Dataset(outputFile, "w") ### 717
    except FileNotFoundError:
        logger.error("Impossible to create file %s" % outputFile)
        sys.exit(4)

    # write attributes
    oFile.title = "SCRIP grid created with scripMaker.py" ### 727
    oFile.Conventions = "CF-1.0" ### 726
    oFile.institution = "Euro-Mediterranean Centre for Climate Change - CMCC" ### 728
    oFile.source = "%s %s" % (gridFile, maskFile)  
    oFile.contact = "Fabio Viola <fabio.viola@cmcc.it>" ### 736
    oFile.creation_date = datetime.datetime.today().strftime("%Y/%m/%d %H:%M") ### 737

    # create dimensions 
    logger.debug("Creating dimension grid_size")
    oFile.createDimension("grid_size", grid_size) ### 741
    logger.debug("Creating dimension grid_rank")
    oFile.createDimension("grid_rank", grid_rank) ### 742
    logger.debug("Creating dimension grid_corners")
    oFile.createDimension("grid_corners", grid_corners) ### 743
    
    # create variable grid_dims
    logger.debug("Creating variable grid_dims")
    gridDimsVar = oFile.createVariable("grid_dims", numpy.dtype('int32').char, ("grid_rank")) ### 746, 747, 748, 749
    gridDimsVar[:] = grid_dims ### 818
    
    # create variable grid_center_lat
    logger.debug("Creating variable grid_center_lat")
    gridCenterLatVar = oFile.createVariable("grid_center_lat", numpy.dtype('double').char, ("grid_size")) ### 752, 753, 754, 761
    gridCenterLatVar.setncattr("units", "radians") ### 755, 756
    gridCenterLatVar.setncattr("long_name", "latitude") ### 757, 758
    gridCenterLatVar.setncattr("bounds", "grid_corner_lat") ### 759, 760
    gridCenterLatVar[:] = rlat ### 819

    # create variable grid_center_lon
    logger.debug("Creating variable grid_center_lon")
    gridCenterLonVar = oFile.createVariable("grid_center_lon", numpy.dtype('double').char, ("grid_size")) ### 764, 765, 767, 773
    gridCenterLonVar.setncattr("units", "radians") ### 767, 768
    gridCenterLonVar.setncattr("long_name", "lonitude") ### 769, 770
    gridCenterLonVar.setncattr("bounds", "grid_corner_lon") ### 771, 772
    gridCenterLonVar[:] = rlon ### 820

    # create variable grid_imask
    logger.debug("Creating variable grid_imask")
    imaskVar = oFile.createVariable("grid_imask", numpy.dtype('int32').char, ("grid_size")) ### 776, 777, 778, 785
    imaskVar.setncattr("units", "1") ### 779, 780
    imaskVar.setncattr("standard_name", "sea_binary_imask") ### 781, 782
    imaskVar.setncattr("long_name", "land-sea imask (1=sea, 0=land)") ### 783, 784
    try:
        imaskVar[:] = mask ### 821
    except:
        pdb.set_trace()

    # create variable grid_area
    logger.debug("Creating variable grid_area")
    areaVar = oFile.createVariable("grid_area", numpy.dtype('double').char, ("grid_size")) # 789, 790, 791, 798
    areaVar.setncattr("units", "sr") ### 792, 793
    areaVar.setncattr("standard_name", "cell_area") ### 794, 795
    areaVar.setncattr("long_name", "area of grid cells (in steradians)") ### 796, 797
    areaVar[:] = tarea ### 822

    # create variable grid_corner_lat
    logger.debug("Creating variable grid_corner_lat")
    gridCornerLatVar = oFile.createVariable("grid_corner_lat", numpy.dtype('double').char, ("grid_size", "grid_corners")) ### 801, 802, 803, 806
    gridCornerLatVar.setncattr("units", "radians") ### 804, 805
    gridCornerLatVar[:] = nemo_cla_final ### 823

    # create variable grid_corner_lon
    logger.debug("Creating variable grid_corner_lon")
    gridCornerLonVar = oFile.createVariable("grid_corner_lon", numpy.dtype('double').char, ("grid_size", "grid_corners")) ### 809, 810, 811, 814
    gridCornerLonVar.setncattr("units", "radians") ### 812, 813
    gridCornerLonVar[:] = nemo_clo_final ### 824
    
    # close output file
    oFile.close()
