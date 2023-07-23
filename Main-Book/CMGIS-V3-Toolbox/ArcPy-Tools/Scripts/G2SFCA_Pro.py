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

facilityFL = arcpy.GetParameterAsText(3)
facilityIDField = arcpy.GetParameterAsText(4)
facilitySizeField = arcpy.GetParameterAsText(5)

distanceType = arcpy.GetParameterAsText(6)
unitConversion = arcpy.GetParameterAsText(7)

distanceDecayCategory = arcpy.GetParameterAsText(8)
distanceDecayFunc = arcpy.GetParameterAsText(9)
distanceDecayCoeffbeta = arcpy.GetParameterAsText(10)
distanceDecayCoeffalpha = arcpy.GetParameterAsText(11)
discretelist = arcpy.GetParameterAsText(12)

matrixTable = arcpy.GetParameterAsText(13)
matrixSupplyId = arcpy.GetParameterAsText(14)
matrixDemandId = arcpy.GetParameterAsText(15)
matrixCost = arcpy.GetParameterAsText(16)

scaleFactor = arcpy.GetParameterAsText(17)
AccName = arcpy.GetParameterAsText(18)
exportMatrix = arcpy.GetParameterAsText(19)
exportfile = arcpy.GetParameterAsText(20)
twosfca = arcpy.GetParameterAsText(21)


# -------------Data Preparing--------------------#
distanceTable = "Temp_ODMatrix"
distanceCustomerID = "IN_FID"
distanceFacilityID = "NEAR_FID"
distanceValue = "NEAR_DIST"

if twosfca == "I2SFCA":
    # switch the input of demand and supplu
    facilityFL = arcpy.GetParameterAsText(0)
    facilityIDField = arcpy.GetParameterAsText(1)
    facilitySizeField = arcpy.GetParameterAsText(2)

    customerFL = arcpy.GetParameterAsText(3)
    customerIDField = arcpy.GetParameterAsText(4)
    customerSizeField = arcpy.GetParameterAsText(5)

    if distanceType == "External Table":
        matrixSupplyId = arcpy.GetParameterAsText(15)
        matrixDemandId = arcpy.GetParameterAsText(14)
        matrixCost = arcpy.GetParameterAsText(16)

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
powerDecayCoeffbeta = "3.0"
exponentialDecayCoeffbeta = "0.03"
sqrtexponentialDecayCoeffbeta = "0.3"
gaussianDecayCoeffbeta = "80"
lognormalDecayCoeffbeta = "0.3"

# two-parameter distance decay function
loglogisticDecayCoeffbeta = "1"
loglogisticDecayCoeffalpha = "1"
comppowerexponentialDecayCoeffbeta = "1"
comppowerexponentialDecayCoeffalpha = "1"

# Check if user specified distance decay coefficient
if distanceDecayCategory != "Discrete":
    if distanceDecayCoeffbeta == "":
        if distanceDecayFunc == "Power":
            distanceDecayCoeffbeta = powerDecayCoeffbeta
        elif distanceDecayFunc == "Exponential":
            distanceDecayCoeffbeta = exponentialDecayCoeffbeta
        elif distanceDecayFunc == "Square-root exponential":
            distanceDecayCoeffbeta = sqrtexponentialDecayCoeffbeta
        elif distanceDecayFunc == "Gaussian":
            distanceDecayCoeffbeta = gaussianDecayCoeffbeta
        elif distanceDecayFunc == "Log-normal":
            distanceDecayCoeffbeta = lognormalDecayCoeffbeta
        elif distanceDecayFunc == "Log-logistic":
            distanceDecayCoeffbeta = loglogisticDecayCoeffbeta
            if distanceDecayCoeffalpha == "":
                distanceDecayCoeffalpha = loglogisticDecayCoeffalpha
                arcpy.AddWarning(
                    "There is no alpha specified for Log-logistic distance decay function and the default alpha = {} is used!".fomrat(
                        loglogisticDecayCoeffalpha
                    )
                )
        elif distanceDecayFunc == "Compound power-exponential":
            distanceDecayCoeffbeta = comppowerexponentialDecayCoeffbeta
            if distanceDecayCoeffalpha == "":
                distanceDecayCoeffalpha = comppowerexponentialDecayCoeffalpha
                arcpy.AddWarning(
                    "There is no alpha specified for Compound power-exponential distance decay function and the default alpha = {} is used!".fomrat(
                        comppowerexponentialDecayCoeffalpha
                    )
                )
            arcpy.AddMessage("There is no distance decay function listed!")
        arcpy.AddMessage("The default coefficients are used!")

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

# *****************************************************************#
arcpy.SetProgressorLabel("Step 1/6: Initiation Finished, Prepare OD Matrix...")
# *****************************************************************#


# Calculate Euclidean distance from customer to facility
if distanceType == "Euclidean":  # output unit is the identical the input features
    arcpy.GenerateNearTable_analysis(
        customer_pt, facility_pt, distanceTable, method="PLANAR", closest="ALL"
    )
elif distanceType == "Geodesic":  # output unit is meter
    arcpy.GenerateNearTable_analysis(
        customer_pt, facility_pt, distanceTable, method="GEODESIC", closest="ALL"
    )
else:  # type is External Table
    arcpy.TableSelect_analysis(matrixTable, distanceTable)

arcpy.SetProgressorLabel("OD Matrix Preparation Done")

# *****************************************************************#
arcpy.SetProgressorLabel("Step 2/6:OD Joining Data to OD Matrix...")
# *****************************************************************#
# -------------Check Distance Table Source---------------------#

desc = arcpy.Describe(distanceTable)
gdbpath = desc.path
distPath = gdbpath + "\\" + distanceTable
matrix = pd.DataFrame.spatial.from_table(distPath)
if distanceType == "External Table":
    matrix1 = matrix[[matrixSupplyId, matrixDemandId, matrixCost]].rename(
        columns={matrixSupplyId: "DID", matrixDemandId: "SID", matrixCost: "Dist"}
    )
else:
    matrix1 = matrix[["IN_FID", "NEAR_FID", "NEAR_DIST"]].rename(
        columns={"IN_FID": "DID", "NEAR_FID": "SID", "NEAR_DIST": "Dist"}
    )
demand = pd.DataFrame.spatial.from_featureclass(customerFL)
supply = pd.DataFrame.spatial.from_featureclass(facilityFL)

demand1 = demand[[customerIDField, customerSizeField]].rename(
    columns={customerIDField: "DID", customerSizeField: "Pop"}
)
supply1 = supply[[facilityIDField, facilitySizeField]].rename(
    columns={
        facilityIDField: "SID",
        facilitySizeField: "Pcp",
    }
)
matrixf = pd.merge(matrix1, supply1, on="SID")
matrixf = pd.merge(matrixf, demand1, on="DID")
matrixf["weight"] = 0
matrixf["Dist"] = matrixf.Dist * float(unitConversion)

mssage = "Apply Unit Conversion Factor by Multiplying " + unitConversion
arcpy.SetProgressorLabel(mssage)
mssage = distanceDecayCategory + " use for Diatance Decay Function"
arcpy.SetProgressorLabel(mssage)

# *****************************************************************#
arcpy.SetProgressorLabel("Step 3/6:Apply Distance Decay Function...")
# *****************************************************************#


def decayfuntbeta(x, beta, func):
    if func == "Power":
        if x == 0:
            return 0
        else:
            return x ** ((-1) * beta)
    if func == "Exponential":
        return math.exp((-1) * x * beta)
    if func == "Square-root exponential":
        return math.exp((-1) * math.sqrt(x) * beta)
    if func == "Gaussian":
        return 1 / (
            math.sqrt(2 * math.pi) * beta * math.exp((-0.5) * x**2 / beta**2)
        )
    if func == "Log-normal":
        if x == 0:
            return 0
        else:
            return math.exp((-1) * beta * (math.log(x)) ** 2)


def decayfuntalphabeta(x, alpha, beta, func):
    if func == "Log-logistic":
        return 1 / (1 + (x / alpha) ** beta)
    if func == "Compound power-exponential":
        return math.exp((-1) * alpha * x**beta)


def normalizelist(x):
    x1 = (x - x.min()) / (x.max() - x.min())
    return x1


if distanceDecayCategory == "Continuous":
    betav = float(distanceDecayCoeffbeta)
    if distanceDecayFunc in ["Log-logistic", "Compound power-exponential"]:
        alphav = float(distanceDecayCoeffalpha)
        matrixf["weight"] = matrixf["Dist"].apply(
            lambda x: decayfuntalphabeta(x, alphav, betav, distanceDecayFunc)
        )
    else:
        matrixf["weight"] = matrixf["Dist"].apply(
            lambda x: decayfuntbeta(x, betav, distanceDecayFunc)
        )
    matrixf["weight"] = normalizelist(matrixf["weight"])
elif distanceDecayCategory == "Discrete":
    Dilist = discretelist.split(";")
    t_list = np.array(Dilist, dtype=np.float32).tolist()
    t_list.sort(reverse=True)
    nlen = len(t_list) + 1
    nlist = []
    for i in range(0, nlen):
        v = 1 / len(t_list)
        v1 = v * (i)
        nlist.append(v1)
    nlist1 = nlist[1:]
    for i in range(len(nlist1)):
        matrixf.loc[(matrixf["Dist"] <= t_list[i]), ["weight"]] = nlist1[i]

else:  # Hybrid
    Dilist = discretelist.split(";")
    t_list = np.array(Dilist, dtype=np.float32).tolist()
    t_list.sort()
    if len(t_list) == 1:
        matrixf1 = matrixf.loc[(matrixf["Dist"] <= t_list[0])]
        matrixf1["weight"] = 1
        matrixf2 = matrixf.loc[(matrixf["Dist"] > t_list[0])]
        # for continuous part
        betav = float(distanceDecayCoeffbeta)
        if distanceDecayFunc in ["Log-logistic", "Compound power-exponential"]:
            alphav = float(distanceDecayCoeffalpha)
            matrixf2["weight"] = matrixf2["Dist"].apply(
                lambda x: decayfuntalphabeta(x, alphav, betav, distanceDecayFunc)
            )
        else:
            matrixf2["weight"] = matrixf2["Dist"].apply(
                lambda x: decayfuntbeta(x, betav, distanceDecayFunc)
            )
        matrixf2["weight"] = normalizelist(matrixf2["weight"])
        matrixf = pd.concat([matrixf1, matrixf2])

    if len(t_list) > 1:
        arcpy.AddMessage("First two threshold values were used for Hybrid Model ")
        t_list1 = t_list[0:2]
        matrixf1 = matrixf.loc[(matrixf["Dist"] <= t_list[0])]
        matrixf2 = matrixf.loc[(matrixf["Dist"] > t_list[1])]
        matrixf1["weight"] = 1
        matrixf2["weight"] = 0
        matrixf3 = matrixf.loc[
            (matrixf["Dist"] > t_list[0]) & (matrixf["Dist"] <= t_list[1])
        ]
        # for continuous part
        betav = float(distanceDecayCoeffbeta)
        if distanceDecayFunc in ["Log-logistic", "Compound power-exponential"]:
            alphav = float(distanceDecayCoeffalpha)
            matrixf3["weight"] = matrixf3["Dist"].apply(
                lambda x: decayfuntalphabeta(x, alphav, betav, distanceDecayFunc)
            )
        else:
            matrixf3["weight"] = matrixf3["Dist"].apply(
                lambda x: decayfuntbeta(x, betav, distanceDecayFunc)
            )
        # normalize
        matrixf3["weight"] = normalizelist(matrixf3["weight"])
        matrixf = pd.concat([matrixf1, matrixf2, matrixf3])

arcpy.SetProgressorLabel("Distance Decay Calculation Done")
# *****************************************************************#
arcpy.SetProgressorLabel("Step 4/6:Calcualte Accessibility...")
# *****************************************************************#
# -----------Calculate dis field with distance decay--------------#

matrixf["wpop"] = matrixf.Pop * matrixf.weight
sumpop = (
    matrixf[["SID", "wpop"]]
    .groupby(["SID"])
    .sum()
    .reset_index()
    .rename(columns={"index": "SID", "wpop": "Sumpop"})
)
matrixf = pd.merge(matrixf, sumpop, on="SID")
matrixf["docpopr"] = float(scaleFactor) * matrixf.Pcp * matrixf.weight / matrixf.Sumpop
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
demand[AccName] = demand["docpopr"]
acc = demand[["DID", AccName]].rename(columns={"DID": customerIDField})

# *****************************************************************#
arcpy.SetProgressorLabel("Step 5/6:Attaching Result to Input Data...")
# *****************************************************************#
# -----------Calculate dis field with distance decay--------------#

acc.spatial.to_table(location=distPath, sanitize_columns=False)
# distPath equals distanceTable

# DelAField(customerFL, AccName)
arcpy.JoinField_management(
    customerFL, customerIDField, distanceTable, customerIDField, [AccName]
)

if exportMatrix == "Yes":
    matrixf.spatial.to_table(location=exportfile, sanitize_columns=False)
# *****************************************************************#
arcpy.SetProgressorLabel("Step 6/6:Cleaning...")
# *****************************************************************#
arcpy.Delete_management(distanceTable)
