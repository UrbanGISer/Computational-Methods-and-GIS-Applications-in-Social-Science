# **************************************************************************************************************
# Date: December 13, 2020
# Author: Changzhen Wang
# Description: 1) This script is used to calculate regionalization indices for pre-existing areas, such as 
#                 1993 Dartmouth HSAs and HRRs
# **************************************************************************************************************


# Import system modules
import arcpy
import os
import numpy as np
import networkx as nx
import DartmouthCommBySpAdjMtrxPro as dtcom


# Get input parameters
inputNHPointFL = arcpy.GetParameterAsText(0)
inputNIDField = arcpy.GetParameterAsText(1)
inputHIDField = arcpy.GetParameterAsText(2)
inputPopField = arcpy.GetParameterAsText(3)
inputEdgeFile = arcpy.GetParameterAsText(4)
edgeOrgIDField = arcpy.GetParameterAsText(5)
edgeDestIDField = arcpy.GetParameterAsText(6)
edgeFlowField = arcpy.GetParameterAsText(7)
edgeDistField = arcpy.GetParameterAsText(8)
inputHPolyFL = arcpy.GetParameterAsText(9)
inputHIDField2 = arcpy.GetParameterAsText(10)
outputHSAs = arcpy.GetParameterAsText(11)


# Join indices of each community (or HSA/HRR) to feature class by keyID and multiple valueIDs
def dict2featurecls2(part_dict, keyID, valueIDs, edgepath, featcls):
	part_list = [(key, *value) for key, value in part_dict.items()]
	dts = {'names': (keyID, valueIDs[0], valueIDs[1], valueIDs[2]), 'formats': (np.int16, np.float, np.float, np.int16)}
	partition_arr = np.rec.fromrecords(part_list, dtype = dts)
	partitionTemp = edgepath + "\\" + "Temp_Parition_File" # must have a path
	if arcpy.Exists(partitionTemp):
		arcpy.Delete_management(partitionTemp)
	arcpy.da.NumPyArrayToTable(partition_arr, partitionTemp)    
	fieldList = arcpy.ListFields(featcls)
	for field in fieldList:
		if field.name == valueIDs[0]:
			arcpy.DeleteField_management(featcls, [valueIDs[0]])
		elif field.name == valueIDs[1]:
			arcpy.DeleteField_management(featcls, [valueIDs[1]])
		elif field.name == valueIDs[2]:
			arcpy.DeleteField_management(featcls, [valueIDs[2]])
		elif field.name == valueIDs[3]:
			arcpy.DeleteField_management(featcls, [valueIDs[3]])
	arcpy.JoinField_management(featcls, keyID, partitionTemp, keyID, [valueIDs[0], valueIDs[1], valueIDs[2]])
	# Calculate the PAC==P/(3.54√A)
	arcpy.AddField_management(featcls, valueIDs[3], "FLOAT", field_alias= valueIDs[3])
	arcpy.CalculateField_management(featcls, valueIDs[3], "!shape.length!/ (3.54 *  math.sqrt( !shape.area! ))", "PYTHON3")
	arcpy.Delete_management(partitionTemp)
	return featcls


DG = nx.DiGraph()
sum_pop = 0
MergeNum = 0
init_partition = dict() # save inital partition from destination ZIP codes
dpartition = dict() # save final partition from inputNIDField and inputHIDField

# read the table in the geodatabase to construct the nodes of a graph
if inputPopField != "":
	with arcpy.da.SearchCursor(inputNHPointFL, [inputNIDField, inputPopField, inputHIDField]) as cursor:
		for row in cursor:
			node = int(row[0])
			if row[2] is not None:
				DG.add_node(node, id = node, pop = int(row[1]))
				sum_pop += int(row[1])
				init_partition[node] = -node
				dpartition[node] = int(row[2])
else:
	with arcpy.da.SearchCursor(inputNHPointFL, [inputNIDField, inputHIDField]) as cursor:
		for row in cursor:
			node = int(row[0])
			if row[1] is not None:
				DG.add_node(node, id = node, pop = 0)
				sum_pop += 0
				init_partition[node] = -node
				dpartition[node] = int(row[1])

# read the table in the geodatabase to construct the edges of a graph
if edgeDistField != "":
	with arcpy.da.SearchCursor(inputEdgeFile, [edgeOrgIDField, edgeDestIDField, edgeFlowField, edgeDistField], """{} > 0""".format(edgeFlowField)) as cursor:
		for row in cursor:
			if row[0] in DG.nodes() and row[1] in DG.nodes():
				DG.add_edge(int(row[0]), int(row[1]), weight=float(row[2]), estTime=float(row[3]))
				init_partition[int(row[1])] = int(row[1])
else:
	with arcpy.da.SearchCursor(inputEdgeFile, [edgeOrgIDField, edgeDestIDField, edgeFlowField], """{} > 0""".format(edgeFlowField)) as cursor:
		for row in cursor:
			if row[0] in DG.nodes() and row[1] in DG.nodes():
				DG.add_edge(int(row[0]), int(row[1]), weight=float(row[2]), estTime=0)
				init_partition[int(row[1])] = int(row[1])

init_coms = len(set(v for v in init_partition.values() if v >= 0)) # allow 0 to be coded as unique value for HSAs or HRRs
arcpy.AddMessage("Total number of nodes is {0}, and total population is {1}".format(DG.number_of_nodes(), sum_pop))             
arcpy.AddMessage("Total number of flows is {0}".format(DG.number_of_edges())) 
arcpy.AddMessage("Total number of service volumes is {0}".format(DG.size(weight = 'weight')))
arcpy.AddMessage("The number of destination ZIP codes is {0}".format(init_coms))

# a temporary file to save the popoluation for each HSA or HRR
aggPopFileTemp = "Temp_Agg_Pop_File"
if arcpy.Exists(aggPopFileTemp):
	arcpy.Delete_management(aggPopFileTemp)

# Delete outputHSAs if exists
arcpy.CopyFeatures_management(inputHPolyFL, outputHSAs)
fieldlist = arcpy.ListFields(outputHSAs)
for field in fieldlist:
	if field.name == inputPopField:
		arcpy.Delete_management(outputHSAs, inputPopField)

if inputPopField != "":
	arcpy.Statistics_analysis(inputNHPointFL, aggPopFileTemp, [[inputPopField, "SUM"], [inputNIDField, "COUNT"]], [inputHIDField])
	newCountNIDField = "COUNT_" + inputNIDField
	newPopField = "SUM_" + inputPopField
	arcpy.AlterField_management(aggPopFileTemp, newPopField, inputPopField, inputPopField)
	arcpy.JoinField_management(outputHSAs, inputHIDField2, aggPopFileTemp, inputHIDField, [inputPopField, newCountNIDField])

# Calculate indices: LI and avgTime
comindices = dtcom.calculateindices(DG, dpartition, init_partition)
inputEdgeFilePath = os.path.split(inputEdgeFile)[0]
outputHSAs = dict2featurecls2(comindices, inputHIDField2, ["LI", "EstTime", "Num_D" + inputNIDField, "Compactness"], inputEdgeFilePath, outputHSAs)
if edgeDistField == "": # if no travel time between nodes, delete the travel time field
	arcpy.DeleteField_management(outputHSAs, "EstTime")

arcpy.Delete_management(aggPopFileTemp)




