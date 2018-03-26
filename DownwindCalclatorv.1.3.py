import arcpy as ac
import netCDF4 as nc
import numpy as np
from arcpy.sa import *
import pandas as pd
import gc

import datetime
T1 = datetime.datetime.now()
print T1


def windcal(u,v):
    # windcal can calculate the overall wind direction based on
    # input u and v wind component
    ws = (u**2 + v**2)**0.5
    wd = np.arctan2(v,u)
    wd_ang = wd *180/np.pi
    wd_ang = wd_ang % 360

    return wd_ang,ws

# Create Name list
ac.env.workspace = r"C:\Users\mozhou\Desktop\GaussianPlumeModel\Split.gdb"
fc = ac.ListFeatureClasses()
name = []
for pt in fc:
    name.append(str(pt))

print "Creating Name list..."
print "the length of name list is: " + str(len(name))

# Read Wind File
print "Loading Wind Data..."
winduv = nc.Dataset(r"C:\Users\mozhou\Desktop\GaussianPlumeModel\WindData\Wind2015.nc",'r')

# Extract Dimension
lat = winduv.variables['latitude'][:]
lon = winduv.variables['longitude'][:]
time = winduv.variables['time'][:]
la = len(lat)
lo = len(lon)
t = len(time)


# Calculate Raster lower left/rght angle
X,Y = np.meshgrid(lon,lat)
lowerleftlat = Y[la-1][0]
lowerleftlong = X[la-1][0]

# Iterate Wind data
index = 365

# Distance Array
MD = []

while index<t:
    # Part I
    # read u and v wind component
    u = winduv.variables['u10'][index,:,:]
    v = winduv.variables['v10'][index,:,:]

    # calculate wind direction and wind speed
    WD, WS = windcal(u,v)
    print "Finish the Calculation of Resultant Wind Direction at Day: " + str(index + 1)
    # convert wind direction/wind speed numpy array to raster object
    ac.env.workspace = r"C:\Users\mozhou\Desktop\GaussianPlumeModel\ScrachRasterSpace"
    WdRaster = ac.NumPyArrayToRaster(WD, ac.Point(lowerleftlong,lowerleftlat), 0.125,0.125,-9999)

    # define the projection for raster file
    # NAD1983(CSRS)
    sr = ac.SpatialReference(4617)
    ac.DefineProjection_management(WdRaster, sr)

    # prject Ratser file, and save it in the scratchworkspace
    prjrs_file = "wd" + str(index+1)
    # set new projection
    proj = ac.SpatialReference(102010)
    ac.ProjectRaster_management (WdRaster,prjrs_file, proj, "NEAREST")

    print "wind direction raster projected " + str(index+1)

    ac.CheckOutExtension("Spatial")

    ac.env.workspace = r'C:\Users\mozhou\Desktop\GaussianPlumeModel\Split.gdb'
    i = 0
    fc = ac.ListFeatureClasses()
    for n in name:
        inPointFeatures = 'C:/Users/mozhou/Desktop/GaussianPlumeModel/Split.gdb/' + n

        ac.env.workspace = r'C:\Users\mozhou\Desktop\GaussianPlumeModel\ScrachRasterSpace'
        inRaster = prjrs_file
        # Execute ExtractValuesToPoints
        print "Extracting WindDirection..."
        ExtractMultiValuesToPoints(inPointFeatures, inRaster, "NONE")

        # Part II
        # Set the temporary bearline and interception path
        bearline = r'C:\Users\mozhou\Desktop\GaussianPlumeModel\Scratch.gdb\bearline'
        InterPt = r'C:\Users\mozhou\Desktop\GaussianPlumeModel\Scratch.gdb\Interpoints'
        points = inPointFeatures
        Road = r'C:\Users\mozhou\Desktop\TestGDB\Roads.gdb\AB_roads_prj'
        # Set the related attributes' name
        distfield = "MaxD"
        bearingfield = prjrs_file
        distunit = 'KILOMETERS'
        angelunit = 'DEGREES'
        linetype = 'GEODESIC'

        # BearingDistanceToLine//!!!!!!add Field X and Y!!!!
        print"Creating Bearing Line..."
        ac.BearingDistanceToLine_management(points, bearline, 'X', 'Y', distfield, distunit, bearingfield, angelunit, linetype)
        # create a intersection point
        print"Finding the Intersection..."
        ac.Intersect_analysis([Road, bearline], InterPt, "", "", 'POINT')


        # PartIII
        # Need to Update the Cursor()
        # Assign Fields

        F1 = ['OBJECTID', 'SHAPE@']
        F2 = ['FID_bearline', 'SHAPE@']

        # create search Cursors
        # facility points
        c1 = ac.da.SearchCursor(points, F1)

        # Create empty list for record distance
        print "Calculating the Downwind-Distance..."
        D =[]
        D.append(index)
        D.append(int(n[1:]))
        # Calculate Distance
        for row1 in c1:
            f_id = row1[0]
            f_geometry = row1[1]
            distlist = []
            # reinitialise the intersection points
            c2 = ac.da.SearchCursor(InterPt, F2)
            for row2 in c2:
                interid = row2[0]
                In_geometry = row2[1]
                if interid == f_id:
                    dis = f_geometry.distanceTo(In_geometry)
                    distlist.append(dis)
            if len(distlist) == 0:
                fac_dist = 20000
            else:
                fac_dist = min(distlist)

            D.append(fac_dist)

        MD.append(D)

        # Delete temporary bearing line and interception points
        print "Clearing the Scratch GeoDatabase..."
        ac.env.workspace = r'C:\Users\mozhou\Desktop\GaussianPlumeModel\Scratch.gdb'
        fcs = ac.ListFeatureClasses()
        for fc in fcs:
            ac.Delete_management(fc)

        i += 1
        print "Feature No."+str(i)+" Finished"
        gc.collect()

    # clear raster scratch workspace
    ac.env.workspace = r"C:\Users\mozhou\Desktop\GaussianPlumeModel\ScrachRasterSpace"
    rasters = ac.ListRasters()
    for rs in rasters:
        ac.Delete_management(rs)
    print "The temporary raster files day "+ str(index+1) + " was deleted"

    print "Day: "+ str(index)

    index += 1

# End of the While Loop
print "Loop end Here, Day: " + str(index)
MD_arr = np.array(MD)

dist_df = pd.DataFrame(MD_arr)
print dist_df.head()

disttable = r"C:\Users\mozhou\Desktop\GaussianPlumeModel\PlumeModel\ABtest3.csv"

dist_df.to_csv(disttable, sep=',')
winduv.close()

T2 = datetime.datetime.now()
print T2-T1
