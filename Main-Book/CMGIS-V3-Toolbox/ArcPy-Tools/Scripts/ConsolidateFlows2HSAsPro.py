# **************************************************************************************************************
# Date: December 7, 2020
# Author: Changzhen Wang
# Description: 1) This script is used to consolidate OD flows from any fine levels to a coarse level, such as HSA 
# 			   2) The originID, destinationID in the input edge file should be consistent to the Unique Node ID
#                 in the Input Node-HSA Feature Layer.
#              3) The Unique HSA ID Field is identical to the HSAID in the Output HSA Node File, and is idential
#				  to O and D in the Output HSA Edge File.
#              4) The sum of Weight in the Output HSA Edge File is equal to the service volumes in Input Edge File (>0)
#			   5) This tool will be automatically computing the number of destination ZIP codes within each HSA
# ***************************************************************************************************************


# Import system modules
import arcpy
import numpy as np


# Get input parameters
inputEdgeFile = arcpy.GetParameterAsText(0)
originIDField = arcpy.GetParameterAsText(1)
destIDField = arcpy.GetParameterAsText(2)
serviceFlowField = arcpy.GetParameterAsText(3)
inputNodeHSAFL = arcpy.GetParameterAsText(4)
inputNodeIDField = arcpy.GetParameterAsText(5)
inputHSAIDField = arcpy.GetParameterAsText(6)
inputPopField = arcpy.GetParameterAsText(7)
outputHSAEdgeFile = arcpy.GetParameterAsText(8)

# Obtain the name and path of the input file
descinputEdgeFile = arcpy.Describe(inputEdgeFile)
inputEdgeFileName = descinputEdgeFile.name
inputEdgeFilepath = descinputEdgeFile.path

# Define temporary files to save input edge file and aggregated edge file
inputEdgeFileTemp = inputEdgeFilepath + "\\" + "Temp_Input_Edge_File" # save non-zero flows
aggEdgeFileTemp = inputEdgeFilepath + "\\" + "Temp_Agg_Edge_File" # save aggregated flows
numDZoneIDFileTemp = inputEdgeFilepath + "\\" + "Temp_Num_DZoneID_File" # save identical records of destID and destHSAID
numDZoneIDFileTemp2 = inputEdgeFilepath + "\\" + "Temp_Num_DZoneID_File2" # save summarized no. of destID within destHSAID


# Delete intermediate data if they already exist
if arcpy.Exists(inputEdgeFileTemp):
	arcpy.Delete_management(inputEdgeFileTemp)	
if arcpy.Exists(aggEdgeFileTemp):
	arcpy.Delete_management(aggEdgeFileTemp)
if arcpy.Exists(numDZoneIDFileTemp):
	arcpy.Delete_management(numDZoneIDFileTemp)
if arcpy.Exists(numDZoneIDFileTemp2):
	arcpy.Delete_management(numDZoneIDFileTemp2)
if arcpy.Exists(outputHSAEdgeFile):
	arcpy.Delete_management(outputHSAEdgeFile)


# Convert input file to .dbf file if not in the geodatabase
if not inputEdgeFilepath.endswith(".gdb"):
	inputEdgeFileTemp = inputEdgeFilepath + "\\" + "Temp_Input_Edge_File.dbf"
	aggEdgeFileTemp = inputEdgeFilepath + "\\" + "Temp_Agg_Edge_File.dbf"
	numDZoneIDFileTemp = inputEdgeFilepath + "\\" + "Temp_Num_DZoneID_File.dbf"
	numDZoneIDFileTemp2 = inputEdgeFilepath + "\\" + "Temp_Num_DZoneID_File2.dbf"

# Define Origin and Destination HSA ID
orginHSAIDField = inputHSAIDField
destHSAIDField = inputHSAIDField + "_1"

# Define field to save no. of dest nodes within HSA
NumDZoneIDField = "Num_DZoneID"

# Delete intermediate fields if they already exist
fieldList = arcpy.ListFields(inputEdgeFile)
for field in fieldList:
	if field.name == inputNodeIDField:
		arcpy.DeleteField_management(inputEdgeFile, [inputNodeIDField])
	if field.name == inputHSAIDField:
		arcpy.DeleteField_management(inputEdgeFile, [inputHSAIDField])
	if field.name == orginHSAIDField:
		arcpy.DeleteField_management(inputEdgeFile, [orginHSAIDField])
	if field.name == destHSAIDField:
		arcpy.DeleteField_management(inputEdgeFile, [destHSAIDField])
	if field.name == NumDZoneIDField:
		arcpy.DeleteField_management(inputEdgeFile, [NumDZoneIDField])


# Select non-zero flows of the input edge file and save in the temporary file
arcpy.TableSelect_analysis(inputEdgeFile, inputEdgeFileTemp, """{0} > 0""".format(serviceFlowField))

# Join HSA ID to the Origin and Destination Node ID in the Edge file
arcpy.JoinField_management(inputEdgeFileTemp, originIDField, inputNodeHSAFL, inputNodeIDField, [inputHSAIDField])
arcpy.JoinField_management(inputEdgeFileTemp, destIDField, inputNodeHSAFL, inputNodeIDField, [inputHSAIDField])

# Summarize the flows based on the output unique ID field
arcpy.Statistics_analysis(inputEdgeFileTemp, aggEdgeFileTemp, [[serviceFlowField, "SUM"]], [orginHSAIDField, destHSAIDField])

# Convert table to array and select the unique pair of (destIDField, destHSAIDField),
# Then convert the array to a new table
destID_HSAID_arr = arcpy.da.TableToNumPyArray(inputEdgeFileTemp, (destIDField, destHSAIDField))
destID_HSAID_arr_uq = np.unique(destID_HSAID_arr, axis = 0)
arcpy.da.NumPyArrayToTable(destID_HSAID_arr_uq, numDZoneIDFileTemp)

# Summarize the number of destination nodes within destination HSAID
arcpy.Statistics_analysis(numDZoneIDFileTemp, numDZoneIDFileTemp2, [[destIDField, "COUNT"]], destHSAIDField)
countDestIDField = "COUNT_" + destIDField
arcpy.AlterField_management(numDZoneIDFileTemp2, countDestIDField, NumDZoneIDField, NumDZoneIDField)

# Redefine the O and D in the new flow table
newOrgIDField = "O_" + inputHSAIDField
newDestIDField = "D_" + inputHSAIDField
newserviceFlowField = "SUM_" + serviceFlowField
newserviceFlowField2 = "N" + serviceFlowField
arcpy.TableSelect_analysis(aggEdgeFileTemp, outputHSAEdgeFile, """{0} > 0""".format(newserviceFlowField))
arcpy.JoinField_management(outputHSAEdgeFile, destHSAIDField, numDZoneIDFileTemp2, destHSAIDField, [NumDZoneIDField])

arcpy.AlterField_management(outputHSAEdgeFile, orginHSAIDField, newOrgIDField, newOrgIDField)
arcpy.AlterField_management(outputHSAEdgeFile, destHSAIDField, newDestIDField, newDestIDField)
arcpy.AlterField_management(outputHSAEdgeFile, newserviceFlowField, newserviceFlowField2, newserviceFlowField2)

# Cleanup intermediate data
arcpy.Delete_management(inputEdgeFileTemp)
arcpy.Delete_management(aggEdgeFileTemp)
arcpy.Delete_management(numDZoneIDFileTemp)
arcpy.Delete_management(numDZoneIDFileTemp2)
