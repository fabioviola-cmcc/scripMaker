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
    model = None
    gridFile = None
    maskFile = None
    outputFile = None
    options, remainder = getopt.getopt(sys.argv[1:], 'o:g:m:', ['output=', 'grid=', 'mask=', 'model='])
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
        elif opt == '--model':
            model = arg
            if model in ["atm", "ocn"]:
                logger.debug("Model set to: %s" % arg)
            else:
                logger.error("Model not supported! Use 'atm' or 'ocn'")
                sys.exit(5)
        else:
            logger.warning("Option %s not recognised" % opt)

    if not(gridFile) or not(outputFile):
        logger.error("You must pass the two input files and an output file name")
        sys.exit(1)

    if not(model):
        logger.error("Model must be specified with --model=ocn|atm")
        sys.exit(6)

        
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
    # set variable names depending on the model
    #    
    ########################################################

    if model == "atm":
        t_lon_var = 'XLONG_M'
        t_lat_var = 'XLAT_M'
        f_lon_var = 'XLONG_U'
        f_lat_var = 'XLAT_V'
    else:
        t_lon_var = 'glamt'
        t_lat_var = 'gphit'
        f_lon_var = 'glamf'
        f_lat_var = 'gphif'        
    
            
    ########################################################
    #
    # read data from the grids and mask (if any)
    #
    ########################################################
    
    # read coordinates on T grid
    t_lon = numpy.squeeze(iFile1.variables[t_lon_var]).transpose()
    t_lat = numpy.squeeze(iFile1.variables[t_lat_var]).transpose()

    # read coordinates on F grid
    if model == "ocn":
        f_lon = numpy.squeeze(iFile1.variables[f_lon_var]).transpose()
        f_lat = numpy.squeeze(iFile1.variables[f_lat_var]).transpose()
    else:        
        f_lon = numpy.squeeze(iFile1.variables[f_lon_var]).transpose()[1:,:]
        f_lat = numpy.squeeze(iFile1.variables[f_lat_var]).transpose()[:,1:]
        
    # read metrics fro T grid
    if model == "ocn":
        e1t = numpy.squeeze(iFile1.variables["e1t"]).transpose()
        e2t = numpy.squeeze(iFile1.variables["e2t"]).transpose()
        tarea = numpy.matrix.flatten(e1t * e2t)
    
    # determine other metrics for T grid
    t_lon_size = t_lon.shape[0]
    t_lat_size = t_lat.shape[1]
    t_grid_size = t_lat_size * t_lon_size
    t_grid_dims = [t_lon_size, t_lat_size]

    # determine general metrics
    grid_corners = 4
    grid_rank = 2

    # read the mask
    if maskFile:
        tmask = numpy.squeeze(iFile2.variables["tmask"])
        if len(tmask.shape) == 3:
            tmask = numpy.squeeze(tmask[0,:,:])

            
    ########################################################
    #
    # processing CENTER Lon and Lat
    #
    ########################################################

    # convert to radians
    rad_center_lon = numpy.radians(t_lon)
    rad_center_lat = numpy.radians(t_lat)

    # flatten
    rad_center_lon = numpy.matrix.flatten(rad_center_lon)
    rad_center_lat = numpy.matrix.flatten(rad_center_lat)

    
    ########################################################
    #
    # processing CORNER Lon and Lat
    #
    ########################################################

    # get the lat and lon cell size
    lon_cell_size = f_lon[1,1] - f_lon[0,1]
    lat_cell_size = f_lat[0,1] - f_lat[0,0]    
    
    # initialise nemo_clo and nemo_cla
    corner_lon = numpy.zeros((t_lon_size, t_lat_size, grid_corners))
    corner_lat = numpy.zeros((t_lon_size, t_lat_size, grid_corners))

    # Corner 0 lon and lat
    logger.debug("Setting corner 0...")
    corner_lon[:,:,0] = f_lon
    corner_lat[:,:,0] = f_lat

    # Corner 1 lon and lat
    # - moving to corner 1, we have problems with longitude when we get to the left border (i.e. first column)
    # - latitude should not represent a problem on the right column
    logger.debug("Setting corner 1...")
    corner_lon[0,:,1] = f_lon[0,:] - lon_cell_size
    corner_lon[1:,:,1] = f_lon[0:-1,:]
    corner_lat[:,:,1] = f_lat

    # Corner 2 lon and lat
    # - this is the most problematic corner, since in the bottom-left end of the grid
    #   we end having a void point both for latitude and longitude
    logger.debug("Setting corner 2...")
    corner_lon[0:,:,2] = f_lon[0:] - lon_cell_size
    corner_lon[1:,:,2] = f_lon[0:-1]
    corner_lat[:,0,2] = f_lat[:,0] - lat_cell_size
    corner_lat[:,1:,2] = f_lat[:,0:-1]

    # Corner 3 lon and lat
    # - moving to corner 3, longitude should not represent a problem
    # - latitude is instead a problem when we get to the bottom of the grid
    logger.debug("Setting corner 3...")
    corner_lon[:,:,3] = f_lon
    corner_lat[:,0,3] = f_lat[:,0] - lat_cell_size
    corner_lat[:,1:,3] = f_lat[:,0:-1]
    
    # convert to radiants
    rad_corner_lat = numpy.radians(corner_lat)
    rad_corner_lon = numpy.radians(corner_lon)

    # reshape nemo_clo
    rad_corner_lon = rad_corner_lon.reshape((t_grid_size, grid_corners))
    
    # reshape nemo_cla
    rad_corner_lat = rad_corner_lat.reshape((t_grid_size, grid_corners))


    ########################################################
    #
    # processing the mask
    #
    ########################################################

    # initialise a mask with all 1
    mask = numpy.ones(t_lon.shape)

    # if there's a mask (e.g. with NEMO, change the values where needed)
    if maskFile:
        mask = numpy.where(tmask>0.0, 1.0, 0.0)

    # transpose and flatten the matrix
    mask = mask.transpose()
    mask = numpy.matrix.flatten(mask)
    
    
    ########################################################
    #
    # open and write output file
    #
    ########################################################
        
    logger.debug("Opening output file...")
    try:
        oFile = Dataset(outputFile, "w")
    except FileNotFoundError:
        logger.error("Impossible to create file %s" % outputFile)
        sys.exit(4)

    # write attributes
    oFile.title = "SCRIP grid created with scripMaker.py"
    oFile.Conventions = "CF-1.0"
    oFile.institution = "Euro-Mediterranean Centre for Climate Change - CMCC"
    oFile.source = "%s %s" % (gridFile, maskFile)  
    oFile.contact = "Fabio Viola <fabio.viola@cmcc.it>"
    oFile.creation_date = datetime.datetime.today().strftime("%Y/%m/%d %H:%M")

    # create dimensions 
    logger.debug("Creating dimension grid_size")
    oFile.createDimension("grid_size", t_grid_size)
    logger.debug("Creating dimension grid_rank")
    oFile.createDimension("grid_rank", grid_rank)
    logger.debug("Creating dimension grid_corners")
    oFile.createDimension("grid_corners", grid_corners)
    
    # create variable grid_dims
    logger.debug("Creating variable grid_dims")
    gridDimsVar = oFile.createVariable("grid_dims", numpy.dtype('int32').char, ("grid_rank"))
    gridDimsVar[:] = t_grid_dims
    
    # create variable grid_center_lat
    logger.debug("Creating variable grid_center_lat")
    gridCenterLatVar = oFile.createVariable("grid_center_lat", numpy.dtype('double').char, ("grid_size"))
    gridCenterLatVar.setncattr("units", "rad")
    gridCenterLatVar.setncattr("long_name", "latitude")
    gridCenterLatVar.setncattr("bounds", "grid_corner_lat")
    gridCenterLatVar[:] = rad_center_lat

    # create variable grid_center_lon
    logger.debug("Creating variable grid_center_lon")
    gridCenterLonVar = oFile.createVariable("grid_center_lon", numpy.dtype('double').char, ("grid_size"))
    gridCenterLonVar.setncattr("units", "rad")
    gridCenterLonVar.setncattr("long_name", "lonitude")
    gridCenterLonVar.setncattr("bounds", "grid_corner_lon")
    gridCenterLonVar[:] = rad_center_lon

    # create variable grid_imask
    logger.debug("Creating variable grid_imask")
    imaskVar = oFile.createVariable("grid_imask", numpy.dtype('int32').char, ("grid_size"))
    imaskVar.setncattr("units", "1")
    imaskVar.setncattr("standard_name", "sea_binary_imask")
    imaskVar.setncattr("long_name", "land-sea imask (1=sea, 0=land)")
    imaskVar[:] = mask

    # create variable grid_area
    if model == "ocn":
        logger.debug("Creating variable grid_area")
        areaVar = oFile.createVariable("grid_area", numpy.dtype('double').char, ("grid_size"))
        areaVar.setncattr("units", "sr")
        areaVar.setncattr("standard_name", "cell_area")
        areaVar.setncattr("long_name", "area of grid cells (in steradians)")
        areaVar[:] = tarea

    # create variable grid_corner_lat
    logger.debug("Creating variable grid_corner_lat")
    gridCornerLatVar = oFile.createVariable("grid_corner_lat", numpy.dtype('double').char, ("grid_size", "grid_corners"))
    gridCornerLatVar.setncattr("units", "rad")
    gridCornerLatVar[:] = rad_corner_lat

    # create variable grid_corner_lon
    logger.debug("Creating variable grid_corner_lon")
    gridCornerLonVar = oFile.createVariable("grid_corner_lon", numpy.dtype('double').char, ("grid_size", "grid_corners"))
    gridCornerLonVar.setncattr("units", "rad")
    gridCornerLonVar[:] = rad_corner_lon
    
    # close output file
    oFile.close()
