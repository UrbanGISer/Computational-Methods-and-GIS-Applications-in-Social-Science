import arcpy
import json
import os

try:
	import urllib2  # Python 2
except ImportError:
	import urllib.request as urllib2  # Python 3

API_Key = arcpy.GetParameterAsText(0)
Input_Origins = arcpy.GetParameterAsText(1)
Input_Destinations = arcpy.GetParameterAsText(2)
Travel_Mode = arcpy.GetParameterAsText(3)
Output_Table = arcpy.GetParameterAsText(4)


# the output distance will have units of miles
basestring = 'https://maps.googleapis.com/maps/api/distancematrix/json?language=en&units=imperial&origins={0}&destinations={1}&key={2}&mode={3}'

# define function to derive travel time and distance from Google Maps API
def fetch_google_OD(download_link):
    nt_time=''
    nt_dist=''
    try:
        req = urllib2.urlopen(download_link)
        jsonout = json.loads(req.read())
        #nt_time = jsonout['rows'][0]['elements'][0]['duration']['text']
        #nt_dist = jsonout['rows'][0]['elements'][0]['distance']['text']
        nt_time = jsonout['rows'][0]['elements'][0]['duration']['value'] #meters
        nt_dist = jsonout['rows'][0]['elements'][0]['distance']['value'] #seconds

        # transform seconds to minutes and meters to miles
        nt_time = round(nt_time / 60, 2)
        nt_dist = round(nt_dist * 0.000621371 , 2)

        '''
        # transform nt_time to minutes especially when there are hours
        nt_time_temp = nt_time.split()
        if "day" in nt_time and "hour" in nt_time:
            if len(nt_time_temp) == 6:
                nt_time = (int(nt_time_temp[0]) * 24 + int(nt_time_temp[2])) *60 + int(nt_time_temp[4])
            elif len(nt_time_temp) == 4:
                nt_time = (int(nt_time_temp[0]) * 24 + int(nt_time_temp[2])) *60
            else:
                pass
        elif "day" not in nt_time and "hour" in nt_time:
            if len(nt_time_temp) == 4:
                nt_time = int(nt_time_temp[0]) * 60 + int(nt_time_temp[2])
            elif len(nt_time_temp) == 2:
                nt_time = int(nt_time_temp[0]) * 60
            else:
                pass
        elif "day" not in nt_time and "hour" not in nt_time:
            nt_time = int(nt_time_temp[0])
        else:
            pass
        
        # transform nt_dist to miles and leave numbers
        nt_dist_temp = nt_dist.split()
        if len(nt_dist_temp) == 2 and "mi" in nt_dist_temp:
            nt_dist = nt_dist_temp[0]
        '''

    except Exception:
        arcpy.AddMessage("The request link is invalid!")  
    return [nt_time, nt_dist]


# write results
File_OutputTable = open(Output_Table, 'w')
File_OutputTable.write('OriginID, DestinationID, Time_Min, Dist_Mile\n')

# read Origins and Destinations, fetch Google travel times
if int(arcpy.GetCount_management(Input_Origins)[0]) == 0:
    arcpy.AddError("{0} has no features.".format(Input_Origins))
    raise arcpy.ExecuteError

if int(arcpy.GetCount_management(Input_Destinations)[0]) == 0:
    arcpy.AddError("{0} has no features.".format(Input_Destinations))
    raise arcpy.ExecuteError

# defined to get lon/lat coordinates
wgs_srf = arcpy.SpatialReference(4326)
Origins_srf = arcpy.Describe(Input_Origins).spatialReference
Destinations_srf = arcpy.Describe(Input_Destinations).spatialReference

# if spatial reference of input features are unknown, define as WGS84
if Origins_srf.name == "Unknown" or Origins_srf.name == "GCS_WGS_1984":
    #arcpy.DefineProjection_management(Input_Origins, wgs_srf) 
    arcpy.AddError("please define a projected coordinate system for the input features")    
if Destinations_srf.name == "Unknown" or Destinations_srf.name == "GCS_WGS_1984":
    #arcpy.DefineProjection_management(Input_Destinations, wgs_srf) 
    arcpy.AddError("please define a projected coordinate system for the input features")

# can process point or polygon features or both
for OriginsRow in arcpy.da.SearchCursor(Input_Origins, ["OID@", "SHAPE@"]):
    OriginsID = OriginsRow[0]
    Origins_Shape = OriginsRow[1].projectAs(wgs_srf)
    Origins_latlon = "{},{}".format(Origins_Shape.centroid.Y, Origins_Shape.centroid.X)

    for DestinationsRow in arcpy.da.SearchCursor(Input_Destinations, ["OID@", "SHAPE@"]):
        DestinationsID = DestinationsRow[0]
        Destinations_Shape = DestinationsRow[1].projectAs(wgs_srf)
        Destinations_latlon = "{},{}".format(Destinations_Shape.centroid.Y, Destinations_Shape.centroid.X)

        arcpy.AddMessage("Estimate Travel time from {0} to {1}...".format(OriginsID, DestinationsID))
        Google_Request_Link = basestring.format(Origins_latlon, Destinations_latlon, API_Key, Travel_Mode.lower())
        print(Google_Request_Link)
        Google_Travel_Cost = fetch_google_OD(Google_Request_Link)
        File_OutputTable.write("{0}, {1}, {2}, {3}\n".format(OriginsID, DestinationsID, Google_Travel_Cost[0], Google_Travel_Cost[1]))
        File_OutputTable.flush()
File_OutputTable.close()





