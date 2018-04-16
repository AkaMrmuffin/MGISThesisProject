# Name:     Downwind Distance Calculator (DDC)
# Purpose:  The purpose of this is to Calculate the average downwind distance from O&G 
#           facilities to nearest road at Alberta.
# Input:    Alberta O&G Facility Feature Class, Alberta Public Road, and ERA-Interim 10-meater U&V Wind Components.   
# Output:   Downwind Distance form O&G Facility to nearest downwind road intersections 
# Author:   Mozhou Gao
# Project:  MGIS Final Proejct 
# Created:  20/02/2018
# Copyright:(c) mozhou.gao 2018

## Loading side Packages
import arcpy as ac
import netCDF4 as nc
import numpy as np
from arcpy.sa import *
import pandas as pd
import os
import fnmatch
import datetime
import gc
from os import listdir
from os.path import isfile, join

## Record start time 
T1 = datetime.datetime.now()
print T1

def windcal(u,v):
    # windcal can calculate the overall wind direction based on
    # Input u and v wind component
    ws = (u**2 + v**2)**0.5
    wd = np.arctan2(v,u)
    wd_ang = wd *180/np.pi
    wd_ang = wd_ang % 360

    return wd_ang,ws

# Read Wind Files
mypath = r'C:\Users\mozhou\Desktop\GaussianPlumeModel\WindData\Data'
files = [f for f in listdir(mypath) if isfile(join(mypath, f))]
for f in files:
    print "Wind Data Loaded"
    winduv = nc.Dataset(r'C:\Users\mozhou\Desktop\GaussianPlumeModel\WindData\Data\%s' % f,'r')


    # Extract Dimensions
    lat = winduv.variables['latitude'][:]
    lon = winduv.variables['longitude'][:]
    time = winduv.variables['time'][:]
    la = len(lat)
    lo = len(lon)
    t = len(time)

    # Calculate cell size
    cellsize = lat[1]-lat[0]
    # Calculate Raster lower left/rght angle
    X,Y = np.meshgrid(lon,lat)
    lowerleftlat = Y[la-1][0]
    lowerleftlong = X[la-1][0]

    # Iterate the day of each year 
    index = 0
    
    # Store the distance
    # Distance Array
    MD = []
    
    # Set the workspace for splitted feature class
    ac.env.workspace =r"C:\Users\mozhou\Desktop\TestGDB\TheSplit.gdb"
    fcs = ac.ListFeatureClasses()
    for pt in fcs:
        # Assign the O&G facilities
        point = pt

        while index<t:
            # Part I
            # read u and v wind component based on time index 
            u = winduv.variables['u10'][index,:,:]
            v = winduv.variables['v10'][index,:,:]

            # Calculate wind direction and wind speed
            WD, WS = windcal(u,v)
            print "Finish the Calculation of Resultant Wind Direction at Day: " + str(index + 1)
            # Convert wind direction/wind speed numpy array to raster object
            WdRaster = ac.NumPyArrayToRaster(WD, arcpy.Point(lowerleftlong,lowerleftlat),
                                            0.125,0.125,-9999)


            # Define the projection for raster file
            # NAD1983(CSRS)
            sr = ac.SpatialReference(4617)
            ac.DefineProjection_management(WdRaster, sr)

            # prject Ratser file, and save it in the scratch workspace
            ac.env.workspace = r"C:\Users\mozhou\Desktop\TestGDB\ScrachRasterSpace"
            prjrs_file = "WD" + str(index+1)
            # set new projection - North_America_Equidistant_Conic (Preserve The Distance)
            proj = ac.SpatialReference(102010)
            ac.ProjectRaster_management (WdRaster,prjrs_file, proj, "NEAREST")

            print "wind direction raster projected " + str(index+1)

            # Check out the ArcGIS Spatial Analyst extension license
            ac.CheckOutExtension("Spatial")

            inPointFeatures = point
            rasterpath = "C:\\Users\\mozhou\\Desktop\\TestGDB\\ScrachRasterSpace\\" + prjrs_file
            # Execute ExtractValuesToPoints
            ExtractMultiValuesToPoints(inPointFeatures, rasterpath, "BILINEAR")

            print "Wind Direction Raster Extracted"

            # Clear raster scratch workspace
            ac.env.workspace = r"C:\Users\mozhou\Desktop\TestGDB\ScrachRasterSpace"
            rasters = ac.ListRasters()
            for rs in rasters:
                ac.Delete_management(rs)

            print "The temporary raster files Number for zone  "+ str(pt) + "of day " + str(index+1)+ "was deleted"

            # Part II
            # Set the temporary bearline and interception path
            bearline =r'C:\Users\mozhou\Desktop\TestGDB\Scratch.gdb\bearline'
            InterPt = r'C:\Users\mozhou\Desktop\TestGDB\Scratch.gdb\Interpoints'
            points = pt
            Road =  r'C:\Users\mozhou\Desktop\TestGDB\Roads.gdb\AB_roads_prj'
            # Set the related attributes' name
            distfield = "MaxD"
            bearingfield = prjrs_file
            distunit = 'KILOMETERS'
            angelunit = 'DEGREES'
            linetype = 'GEODESIC'
            # BearingDistanceToLine
            ac.BearingDistanceToLine_management(points, bearline, 'X', 'Y', distfield,distunit,bearingfield, angelunit, linetype)

            # Create a intersection point
            ac.Intersect_analysis ([Road, bearline],InterPt, "", "", 'POINT')

            print"Bearline and intersection points were created for " + str(pt)

            #PartIII
            # Need to Update the Cursor()
            #Assign Fields

            F1 = ['OBJECTID','SHAPE@']
            F2 = ['FID_bearline','SHAPE@']

            # Create search Cursors
            # Facility points
            c1 =ac.da.SearchCursor(points,F1)

            # Create empty list for record distance
            D = []
            # Feature ID 
            D.append(int(pt[1:5]))
            # Time ID
            D.append(index+1)
            # Calculate Distance
            for row1 in c1:

                f_id = row1[0]
                f_geometry = row1[1]
                distlist = []
                # Reinitialise the intersection points
                c2 = ac.da.SearchCursor(InterPt,F2)
                for row2 in c2:
                    interid  = row2[0]
                    In_geometry = row2[1]

                    if interid  == f_id:
                        dis = f_geometry.distanceTo(In_geometry)
                        distlist.append(dis)

                if len(distlist) == 0:
                    fac_dist = 20000
                else:
                    fac_dist = min(distlist)

                D.append(fac_dist)

            MD.append(D)
            # Delete temporary bearline and interception points
            ac.env.workspace = r'C:\Users\mozhou\Desktop\TestGDB\Scratch.gdb'
            fcs = arcpy.ListFeatureClasses()
            for fc in fcs:
                arcpy.Delete_management(fc)

            print "The scratch GDB is clear Now!"

            index = index + 1
            print "Day: " + str(index)
        
        # Report the finished feature class 
        print str(pt) + " Done!"


    
    # Report the finished year
    print "Good Luck with this Year!"
    # Convert result to array
    MD_arr=np.array(MD)
    # Close the .nc file  
    winduv.close()
    # Save the result of current year
    dist_df = pd.DataFrame(MD_arr)
    dist_df.to_csv(r"C:\Users\mozhou\Desktop\TestGDB\dist2%s.csv" % f,sep=',')

    
print 'Calculation Finished!!!'
# report the overall calculating time
T2 = datetime.datetime.now()
print T2-T1

