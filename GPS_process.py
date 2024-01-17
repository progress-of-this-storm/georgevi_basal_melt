# let's do this
import os
import fnmatch
import numpy
import rasterio as rio
import matplotlib
import matplotlib.pyplot as plt
import rioxarray
from rasterio.plot import show
import geopandas as gpd
import scipy.io
import datetime
import pandas as pd
import datetime
from math import sin, cos, asin, sqrt, pi
from datetime import date, timedelta, time

path = '/Users/louie.bell/Cambridge/mphil/Essay (Practical) 3/data/DeTided_DeIBEd_GPS_data/'
directory = os.listdir('/Users/louie.bell/Cambridge/mphil/Essay (Practical) 3/data/DeTided_DeIBEd_GPS_data')
matFiles = [path+fname for fname in directory if fnmatch.fnmatch(fname, '*.mat')]
matNames = []
for fname in directory:
    if fnmatch.fnmatch(fname, '*.mat'):
        name, ext = os.path.splitext(fname)
        matNames.append(name)

print(str(len(matFiles)) + ' files read, filenames:')
for i in range(len(matFiles)):
    print('\t',i, matFiles[i])   
    
class gpsData():
    
    def __init__(self, filepath):
        self.filepath = filepath
        
    def readMat(self): 
        matFile = scipy.io.loadmat(self.filepath)
        return matFile
          
    def getArray(self, matFile, arrayName):
        matArray = numpy.array(matFile[arrayName])
        return matArray
    
# function to compute velocity based on Haversine distance and days elapsed
def processGPSVel(df):
    def distance(lat1, lon1, lat2, lon2):
        r = 6371 # km
        p = pi / 180
        a = 0.5 - cos((lat2-lat1)*p)/2 + cos(lat1*p) * cos(lat2*p) * (1-cos((lon2-lon1)*p))/2
        return 2 * r * asin(sqrt(a))
    
    if len(df) == 0:
        print('\t Processing failed: GPS data has length zero')
    else:
        
        seasonStartLat = df.loc[0][3]
        seasonEndLat = df.loc[len(df)-1][3]
        seasonStartLon = df.loc[0][4]
        seasonEndLon = df.loc[len(df)-1][4]

        d0 = datetime.date(df.loc[0][0], df.loc[0][1], df.loc[0][2])
        d1 = datetime.date(df.loc[len(df)-1][0], df.loc[len(df)-1][1], df.loc[len(df)-1][2])
        delta = d1-d0

        d_km = distance(seasonStartLat, seasonStartLon, seasonEndLat, seasonEndLon)
        d_m = d_km*1000

        vel_d = d_m/delta.days
        print(delta.days)
        return vel_d

seasons = ['DJF', 'MAM', 'JJA', 'SON']
years = [2019, 2020, 2021, 2022]
yearSeasons = []
for year in years:
    for season in seasons:
        yearSeason = str(year) + season
        yearSeasons.append(yearSeason)
yearSeasons

seasonGeoms = {}
for i in matNames:
    seasonGeoms[i] = {}

velocities = pd.DataFrame(columns = ['GPS', 'Sat'], index = yearSeasons)
GPSvelocities = []
for i in range(len(matFiles)):
    GPSdat = gpsData(matFiles[i])
    gpsMat = GPSdat.readMat()
    arrayNames = [x for x in gpsMat]
    #siteName = arrayNames[12].split("_")[1]  # get unique part of each array
    siteName = matNames[i]
    print('Site name processing: ' + siteName)
    gpsLat = GPSdat.getArray(gpsMat, arrayNames[9]) # get lat
    gpsLon = GPSdat.getArray(gpsMat, arrayNames[10]) # get lon
    gpsLon = gpsLon-360 # convert to more readable longitude
    gpsTime = GPSdat.getArray(gpsMat, arrayNames[11]) # get time in decimal days since 0000

    # function to process numpy array into more workable list
    def processArray(array):
        list = []
        for i in range(len(array)):
            val = array[i][0]
            list.append(val)
        return list

    gpsLatVals = processArray(gpsLat)
    gpsLonVals = processArray(gpsLon)
    gpsTimeVals = processArray(gpsTime)

    ## create lookup list of all possible decimal days
    lookupDecDays = numpy.arange(737000, 740000,1) # create boundary for all possible lookup dates
    
    # 2017 to 2026 will be our deliberately broad temporal envelope (data spans from 2019 to 2022)
    
    # create lookup list of all possible datetimes
    startDate = datetime.date(2017, 11, 2)
    dates = []
    dates.append(startDate)
    delta = timedelta(days = 1)
    for i in range(len(lookupDecDays)):
        date = dates[i]
        newDate = date + delta
        dates.append(newDate)

    # create lookup as dictionary 
    day_lookup = dict(zip(lookupDecDays, dates))
    
    # create decimal lookup for raw data
    def processDayFraction(day):
        day_frac = day % 1 # get day fraction with mod 1
        hrs = day_frac*24 # get day value in hours
        hrs_frac = hrs % 1 # get fractional part with modulo 1
        mins = hrs_frac*60 # get remainder in mins
        timeCor = time(hour = int(hrs), minute = int(mins)) # get corrected time
        return timeCor
    
    # get list of all corrected times
    times = [processDayFraction(x) for x in gpsTimeVals]
    #for i in gpsTimeVals:
    #    timeFormat = processDayFraction(i)
    #    times.append(timeFormat)

    # get list of decimal days only
    decDaysOnly = [int(x) for x in gpsTimeVals]
    
    pd.set_option("display.precision", 10)
    mypd = pd.DataFrame(data = {'Lat':gpsLatVals,
                                'Lon':gpsLonVals,
                                'decDate':gpsTimeVals,
                                'decDay':decDaysOnly,
                                'Time':times})
    mypd['Date'] = mypd['decDay'].map(day_lookup)

    datetimes = []
    for i in range(len(gpsTimeVals)):
        dateTime = datetime.datetime.combine(mypd['Date'][i],mypd['Time'][i])
        datetimes.append(dateTime)

    mypd['datetime'] = datetimes # add combined datetimes to the DataFrame
    
    ## lookup lat/lon based on season
    
    seasonMap = {
        'DJF':[12, 1, 2],
        'MAM':[3,4,5],
        'JJA':[6,7,8],
        'SON':[9,10,11]
    }
    
    for gpsyear in numpy.arange(2019,2023,1):
        
        for gpsseason in seasonMap:
            
            seasonArray = []
            for i in range(len(mypd)):
                iMonth = mypd.loc[i]['Date'].month
                iYear = mypd.loc[i]['Date'].year
                iDay = mypd.loc[i]['datetime'].day

                if iYear == gpsyear and iMonth in seasonMap[gpsseason]:
                    row = [iYear, iMonth, iDay, mypd.loc[i]['Lat'], mypd.loc[i]['Lon']]
                    seasonArray.append(row)

            seasonData = pd.DataFrame(columns = ['Year', 'Month','Day', 'Lat', 'Lon'],
                                      index = list(numpy.arange(0,len(seasonArray))))

            for i in range(len(seasonArray)):
                seasonData.loc[i] = seasonArray[i]
                
            lenSeason = len(seasonData)
            print(f'Length of {gpsyear},{gpsseason} is: {lenSeason}')    

            gpsVel = processGPSVel(seasonData)
            
            # get GPS geometries for each season
            from shapely.geometry import Point
            geometrySeason = [Point(xy) for xy in zip(seasonData.Lon, seasonData.Lat)] # get geom points
            newdf = seasonData.drop(['Lon', 'Lat'], axis = 1) # new df without lat and lon static vals
            gdf = gpd.GeoDataFrame(newdf, crs = 'EPSG:4326', geometry = geometrySeason) # new gdf with geometry col
            gdf = gdf.to_crs('EPSG:3031') # convert to crs
            
            geoms=[]
            for i in gdf['geometry']:
                x = i.x
                y = i.y
                geom = (x,y) # get list of geometry tuples
                geoms.append(geom)
                
            ys = str(gpsyear) + gpsseason # get yearSeason
            
            seasonGeoms[siteName][ys] = geoms # append nested dictionary with yearSeason and list of geom tuples
            
            print('\tVelocity of site ' + siteName +':\n\t' + str(gpsyear) + ' ' + gpsseason +": " + str(gpsVel) + ' m/d')
            
            GPSvelocities.append(gpsVel)
            
