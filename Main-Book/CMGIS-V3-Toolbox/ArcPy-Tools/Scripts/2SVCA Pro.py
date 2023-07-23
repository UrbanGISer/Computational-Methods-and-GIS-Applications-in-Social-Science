# **************************************************************************************************************
# Date: December 3, 2022
# Author: Lingbo Liu
# Description: This script is used to calculate healthcare service accessibility based on two-step virtual
#              floating catchment method ( 2SVCA)
# InputFeature 1: (1)supply layer, and its (2) ID,(3) supply amount(doctors or beds),(4) broadband speed
# InputFeature 2: (1)demand layer, and its (2) ID,(3) deamnd amount(population),(4) broadband speed,
#                 (5) broadband subscription rate
# InputFeature 3: (1)OD matrix table, and its (2) supply ID,(3) demand ID,(4) distance value,
# 				  (5) distance threshod
# OutputFeature: accessibility table
# **************************************************************************************************************

# Import system modules
import arcpy
import numpy as np
import pandas as pd
import math
from arcgis.features import GeoAccessor, GeoSeriesAccessor

# Get input parameters
customerFL = arcpy.GetParameterAsText(0)
customerIDField = arcpy.GetParameterAsText(1)
customerSizeField = arcpy.GetParameterAsText(2)
# Add broadbandAcc 1/0 info for supply
customerBand = arcpy.GetParameterAsText(3)
customerBandAvail = arcpy.GetParameterAsText(4)

facilityFL = arcpy.GetParameterAsText(5)
facilityIDField = arcpy.GetParameterAsText(6)
facilitySizeField = arcpy.GetParameterAsText(7)
facilityBand = arcpy.GetParameterAsText(8)
facilityBandAvail = arcpy.GetParameterAsText(9)

distanceType = arcpy.GetParameterAsText(10)
Threshold = arcpy.GetParameterAsText(11)

matrixTable = arcpy.GetParameterAsText(12)
matrixDemandId = arcpy.GetParameterAsText(13)
matrixSupplyId = arcpy.GetParameterAsText(14)
matrixCost = arcpy.GetParameterAsText(15)

scaleFactor = arcpy.GetParameterAsText(16)
AccName = arcpy.GetParameterAsText(17)


# -------------Data Preparing--------------------#
distanceTable = "Temp_ODMatrix"
distanceCustomerID = "IN_FID"
distanceFacilityID = "NEAR_FID"
distanceValue = "NEAR_DIST"


# Delete a field if it exists in the input feature class (or layer).
def DelAField(theFC, theField):
    flds = arcpy.ListFields(theFC, theField)
    for fld in flds:
        if fld:
            arcpy.DeleteField_management(theFC, theField)


# Add a field if it does not exist in the input feature class (or layer).
def AddAField(theFC, theField, fldType):
    class GetOutOfDef(Exception):
        pass

    try:
        flds = arcpy.ListFields(theFC, theField)
        for fld in flds:
            if fld:
                raise GetOutOfDef
        arcpy.AddField_management(theFC, theField, fldType)
    except GetOutOfDef:
        pass


# Delete intermediate data if they already exist
if arcpy.Exists("Temp_ODMatrix"):
    arcpy.Delete_management("Temp_ODMatrix")

if arcpy.Exists("Temp_Customer_Pt"):
    arcpy.Delete_management("Temp_Customer_Pt")

if arcpy.Exists("Temp_Facility_Pt"):
    arcpy.Delete_management("Temp_Facility_Pt")

# Assign default distance decay coefficient
# one-parameter distance decay function

# -------------Check Demand and Supply Features--------------------#
if distanceType in ["Euclidean", "Geodesic"]:
    # Check feature layer type and convert it to centroid layer if necessary
    customer_pt = customerFL
    desc = arcpy.Describe(customerFL)
    if desc.shapeType == "Polygon":
        customer_pt = "Temp_Customer_Pt"
        arcpy.FeatureToPoint_management(customerFL, customer_pt, "INSIDE")
        # Need to change customer id to NEW_FID for shapefile, which is FID + 1
    # if the input features are in the geographic coordinate system, set it as geodesic
    customer_srf = desc.spatialReference
    if customer_srf.name == "Unknown" or customer_srf.type != "Projected":
        distanceType == "Geodesic"

    facility_pt = facilityFL
    desc = arcpy.Describe(facilityFL)
    if desc.shapeType == "Polygon":
        facility_pt = "Temp_Facility_Pt"
        arcpy.FeatureToPoint_management(facilityFL, facility_pt, "INSIDE")
        # Need to change facility id to NEW_FID for shapefile, which is FID + 1

mssage = distanceType + " used for Diatance "
arcpy.SetProgressorLabel(mssage)

# Calculate Euclidean distance from customer to facility
if distanceType in ["Euclidean", "Geodesic"]:
    if distanceType == "Euclidean":  # output unit is the identical the input features
        arcpy.GenerateNearTable_analysis(
            customer_pt, facility_pt, distanceTable, method="PLANAR", closest="ALL"
        )
    elif distanceType == "Geodesic":  # output unit is meter
        arcpy.GenerateNearTable_analysis(
            customer_pt, facility_pt, distanceTable, method="GEODESIC", closest="ALL"
        )
    desc = arcpy.Describe(distanceTable)
    print()
    gdbpath = desc.path
    distPath = gdbpath + "\\" + distanceTable
    matrix = pd.DataFrame.spatial.from_table(distPath)
    matrix1 = matrix[["IN_FID", "NEAR_FID", "NEAR_DIST"]].rename(
        columns={"IN_FID": "DID", "NEAR_FID": "SID", "NEAR_DIST": "Dist"}
    )
else:  # type is External Table
    matrix = pd.DataFrame.spatial.from_table(matrixTable)
    matrix1 = matrix[[matrixSupplyId, matrixDemandId, matrixCost]].rename(
        columns={matrixSupplyId: "SID", matrixDemandId: "DID", matrixCost: "Dist"}
    )

arcpy.SetProgressorLabel("OD Matrix Preparation Done")

# -------------Check Distance Table Source---------------------#

demand = pd.DataFrame.spatial.from_featureclass(customerFL)
supply = pd.DataFrame.spatial.from_featureclass(facilityFL)

supply1 = supply[[facilityIDField, facilitySizeField, facilityBand]].rename(
    columns={facilityIDField: "SID", facilitySizeField: "Pcp", facilityBand: "bi"}
)
# set band affordablity
if facilityBandAvail != "":
    supply1["bandRF"] = supply[facilityBandAvail].astype(float)
else:
    supply1["bandRF"] = 1.0

demand1 = demand[[customerIDField, customerSizeField, customerBand]].rename(
    columns={
        customerIDField: "DID",
        customerSizeField: "Pop",
        customerBand: "bk",
    }
)
# set band affordablity
if customerBandAvail != "":
    demand1["bandR"] = demand[customerBandAvail].astype(float)
else:
    demand1["bandR"] = 1.0


matrix1 = matrix1[matrix1.Dist <= float(Threshold)]
matrixf = pd.merge(matrix1, supply1, on="SID")
matrixf = pd.merge(matrixf, demand1, on="DID")
matrixf = matrixf.apply(pd.to_numeric)
arcpy.AddMessage("Data Merged...")

# -------------Calculation--------------------------------#
# matrixf["wpop"] = matrixf.Pop * matrixf.bk * matrixf.bandR
matrixf["wpop"] = matrixf.Pop * matrixf.bk * matrixf.bandR
sumpop = (
    matrixf[["SID", "wpop"]]
    .groupby(["SID"])
    .sum()
    .reset_index()
    .rename(columns={"index": "SID", "wpop": "Sumpop"})
)
matrixf = pd.merge(matrixf, sumpop, on="SID")
matrixf["docpopr"] = (
    float(scaleFactor) * matrixf.bi * matrixf.bandRF * matrixf.Pcp / matrixf.Sumpop
)
sumdocpopr = (
    matrixf[["DID", "docpopr"]]
    .groupby(["DID"])
    .sum()
    .reset_index()
    .rename(columns={"index": "DID"})
)
demand = pd.merge(
    demand, sumdocpopr, how="left", left_on=customerIDField, right_on="DID"
)
demand["docpopr"] = demand["docpopr"].fillna(0)
demand = demand.rename(columns={customerBandAvail: "bandR"})
demand["docpopr"] = demand.docpopr * demand.bandR
acc = demand[[customerIDField, "docpopr"]]
acc = acc.rename(columns={"docpopr": AccName})

arcpy.SetProgressorLabel("Calcualtion Done, Saving")
acc.spatial.to_table(location=distPath, sanitize_columns=False)
# distPath equals distanceTable
# DelAField(customerFL, AccName)
arcpy.JoinField_management(
    customerFL, customerIDField, distanceTable, customerIDField, [AccName]
)
arcpy.Delete_management(distanceTable)
