# **************************************************************************************************************
# Date: December 13, 2020
# Author: Changzhen Wang
# Description: This script is used to delineate HSAs or HRRs by the refined Dartmouth method
# **************************************************************************************************************


# Import system modules
import arcpy
import os
import copy
import numpy as np
import networkx as nx
import scipy
import DartmouthCommBySpAdjMtrxPro as dtcom


# Get input parameters
inputPolyFL = arcpy.GetParameterAsText(0)
inputIDField = arcpy.GetParameterAsText(1)
inputPopField = arcpy.GetParameterAsText(2)
inputEdgeFile = arcpy.GetParameterAsText(3)
edgeOrgIDField = arcpy.GetParameterAsText(4)
edgeDestIDField = arcpy.GetParameterAsText(5)
edgeFlowField = arcpy.GetParameterAsText(6)
edgeDistField = arcpy.GetParameterAsText(7)
delineateMethod = arcpy.GetParameterAsText(8)
inputSpAdjMatrix = arcpy.GetParameterAsText(9)
thresholdSize = arcpy.GetParameterAsText(10)
miniLocalIndex = arcpy.GetParameterAsText(11)
outputType = arcpy.GetParameterAsText(12)
outputHSAs = arcpy.GetParameterAsText(13)

if thresholdSize == "":
	thresholdSize = None
elif thresholdSize == "1,000":
	thresholdSize = 1000
elif thresholdSize == "120,000":
	thresholdSize = 120000
else:
	thresholdSize = int(thresholdSize)

if miniLocalIndex == "":
	miniLocalIndex = None
else:
	miniLocalIndex = float(miniLocalIndex)


# Join partition results to feature class by nodeID and valueID
def dict2featurecls(part_dict, keyID, valueID, edgepath, featcls):
	fieldList = arcpy.ListFields(featcls)
	for field in fieldList:
		if field.name == valueID:
			arcpy.DeleteField_management(featcls, [valueID])			
	part_dict1 = {key : value for key, value in part_dict.items()}
	partition_arr = np.array(list(part_dict1.items()), dtype=[(keyID, np.int64),(valueID, np.int64)])
	partitionTemp = edgepath + "\\" + "Temp_Parition_File"
	if arcpy.Exists(partitionTemp):
		arcpy.Delete_management(partitionTemp)
	arcpy.da.NumPyArrayToTable(partition_arr, partitionTemp)	
	arcpy.JoinField_management(featcls, keyID, partitionTemp, keyID, [valueID])
	arcpy.Delete_management(partitionTemp)
	return featcls


# Join indices of each community (or HSA/HRR) to feature class by keyID and multiple valueIDs
def dict2featurecls2(part_dict, keyID, valueIDs, edgepath, featcls):	   
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
	part_list = [(key, *value) for key, value in part_dict.items()]
	dts = {'names': (keyID, valueIDs[0], valueIDs[1], valueIDs[2]), 'formats': (np.int16, np.float, np.float, np.int16)}
	partition_arr = np.rec.fromrecords(part_list, dtype = dts)
	partitionTemp = edgepath + "\\" + "Temp_Parition_File"
	if arcpy.Exists(partitionTemp):
		arcpy.Delete_management(partitionTemp)
	arcpy.da.NumPyArrayToTable(partition_arr, partitionTemp)
	arcpy.JoinField_management(featcls, keyID, partitionTemp, keyID, [valueIDs[0], valueIDs[1], valueIDs[2]])
	# Calculate the PAC==P/(3.54√A)
	arcpy.AddField_management(featcls, valueIDs[3], "FLOAT", field_alias= valueIDs[3])
	arcpy.CalculateField_management(featcls, valueIDs[3], "!shape.length!/ (3.54 *  math.sqrt( !shape.area! ))", "PYTHON3")
	arcpy.Delete_management(partitionTemp)
	return featcls

DG = nx.DiGraph()
init_partition = {} # save inital partition from destination ZIP codes

# Convert preliminary results to dbf table with two fields, ZoneID and ComID
addComID = "HSAID"
type = "HSAs"
if inputIDField == addComID:
	addComID = "HRRID"
	type = "HRRs"

# read the table in the geodatabase to construct the nodes of a graph
if inputPopField != "":
	with arcpy.da.SearchCursor(inputPolyFL, [inputIDField, inputPopField]) as cursor:
		for row in cursor:
			DG.add_node(int(row[0]), id = int(row[0]), pop = int(row[1]))
			init_partition[int(row[0])] = -int(row[0])
else:
	with arcpy.da.SearchCursor(inputPolyFL, [inputIDField]) as cursor:
		for row in cursor:
			DG.add_node(int(row[0]), id = int(row[0]), pop = 0)
			init_partition[int(row[0])] = -int(row[0])				

# Save the number of destination ZIP codes within the HRR
HRR_NumDNode = {}
# Define field to save no. of dest nodes within HSA
NumDZoneIDField = "Num_DZoneID"
inputPolyFL2 = inputPolyFL + "_CPTemp" # avoid the incorrect computation of NumDZoneIDField

if arcpy.Exists(inputPolyFL2):
	arcpy.Delete_management(inputPolyFL2)

# read the table in the geodatabase to construct the edges of a graph
inputEdgePath = os.path.split(inputEdgeFile)[0]
inputEdgeFdNames= [field.name for field in arcpy.ListFields(inputEdgeFile)]
if edgeDistField != "":
	if NumDZoneIDField in inputEdgeFdNames and type == "HRRs":		
		with arcpy.da.SearchCursor(inputEdgeFile, [edgeOrgIDField, edgeDestIDField, edgeFlowField, edgeDistField, NumDZoneIDField], """{} > 0""".format(edgeFlowField)) as cursor:
			for row in cursor:
				orgIDValue = int(row[0])
				destIDValue = int(row[1])
				DG.add_edge(orgIDValue, destIDValue, weight=float(row[2]), estTime=float(row[3]))
				init_partition[destIDValue] = destIDValue

				if not bool(HRR_NumDNode.get(destIDValue)):
					HRR_NumDNode[destIDValue] = int(row[4])
	else:
		with arcpy.da.SearchCursor(inputEdgeFile, [edgeOrgIDField, edgeDestIDField, edgeFlowField, edgeDistField], """{} > 0""".format(edgeFlowField)) as cursor:
			for row in cursor:
				orgIDValue = int(row[0])
				destIDValue = int(row[1])
				DG.add_edge(orgIDValue, destIDValue, weight=float(row[2]), estTime=float(row[3]))
				init_partition[destIDValue] = destIDValue
else:
	if NumDZoneIDField in inputEdgeFdNames and type == "HRRs":
		with arcpy.da.SearchCursor(inputEdgeFile, [edgeOrgIDField, edgeDestIDField, edgeFlowField, NumDZoneIDField], """{} > 0""".format(edgeFlowField)) as cursor:
			for row in cursor:
				orgIDValue = int(row[0])
				destIDValue = int(row[1])	
				DG.add_edge(orgIDValue, destIDValue, weight=float(row[2]), estTime=0)
				init_partition[destIDValue] = destIDValue

				if not bool(HRR_NumDNode.get(destIDValue)):
					HRR_NumDNode[destIDValue] = int(row[3])
	else:
		with arcpy.da.SearchCursor(inputEdgeFile, [edgeOrgIDField, edgeDestIDField, edgeFlowField], """{} > 0""".format(edgeFlowField)) as cursor:
			for row in cursor:
				orgIDValue = int(row[0])
				destIDValue = int(row[1])	
				DG.add_edge(orgIDValue, destIDValue, weight=float(row[2]), estTime=0)
				init_partition[destIDValue] = destIDValue

# Extract a subnetwork (only perserving maximal outgoing flow volumes) for using Huff-Dartmouth method
subDG = nx.DiGraph()
if delineateMethod == "Dartmouth Method":
	subDG = copy.deepcopy(DG)
elif delineateMethod == "Huff-Dartmouth Method":
	for node in DG.nodes():
		subDG.add_node(node, id = node, pop = DG.nodes[node]['pop'])
		out_neighbor_weight = dict()
		out_neighbor_estTime = dict()
		out_neighbors = list(DG.successors(node))
		for out_neighbor in out_neighbors:
			out_neighbor_weight[out_neighbor] = DG.get_edge_data(node, out_neighbor, {"weight": 0}).get("weight", 0)
			out_neighbor_estTime[out_neighbor]= DG.get_edge_data(node, out_neighbor, {"estTime": 0}).get("estTime", 0)			
		if len(out_neighbor_weight) > 0:
			out_neighbor_sls = {key: value for (key, value) in out_neighbor_weight.items() if value == max(out_neighbor_weight.values())}
			for destnode, odweight in out_neighbor_sls.items():
				subDG.add_edge(node, destnode, weight = odweight, estTime = out_neighbor_estTime[destnode])
else:
	arcpy.AddError("Please ensure you select either Dartmouth Method or Huff-Dartmouth Method!")


# Show detailed information about the network in the Messages tab of the tool after completing the computation
init_coms = len(set(v for v in init_partition.values() if v >= 0))
subnode_pop = dict(subDG.nodes(data = 'pop'))
arcpy.AddMessage("Total number of nodes is {0}, and total population is {1}".format(subDG.number_of_nodes(), sum(subnode_pop.values())))             
arcpy.AddMessage("Total number of flows is {0}".format(subDG.number_of_edges())) 
arcpy.AddMessage("Total number of service volumes is {0}".format(subDG.size(weight = 'weight')))
arcpy.AddMessage("The number of destination ZIP codes is {0}".format(init_coms))

# Delete outputHSAs if exists
if arcpy.Exists(outputHSAs):
	arcpy.Delete_management(outputHSAs)

# Read the spatial adjacency matrix
adjmx = None
if inputSpAdjMatrix == "":
	arcpy.AddError("Please copy the path of .npz file to the Input Spatial Adjacency Matrix File!")
elif not inputSpAdjMatrix.endswith(".npz"):
	arcpy.AddError("Please ensure the Input Spatial Adjacency Matrix File ends with .npz format!")
else:
	adjmx = scipy.sparse.load_npz(inputSpAdjMatrix).todense().nonzero()

# Save the final partition result
dpartition = dict()

# Two enforcements both need population field
if thresholdSize is None or thresholdSize <= 0 :
	arcpy.AddMessage("The output {} do not impose any population constraint! ".format(type))
	if miniLocalIndex is None or miniLocalIndex <= 0:
		arcpy.AddMessage("And the ouput {} do not impose any LI constraint!".format(type))
		dpartition = dtcom.find_partition_Dartmouth(subDG, adjmx, init_partition, None, None)
	else:
		if inputPopField == "":
			arcpy.AddError("Please input the Population Field!")
		else:
			arcpy.AddMessage("But the ouput {} impose LI constraint!".format(type))
			dpartition = dtcom.find_partition_Dartmouth(subDG, adjmx, init_partition, None, miniLocalIndex)
else:
	arcpy.AddMessage("The output {} impose population constraint!".format(type))
	if inputPopField == "":
		arcpy.AddError("Please input Population Field!")
	else:
		if miniLocalIndex is None or miniLocalIndex <= 0:
			arcpy.AddMessage("But the ouput {} do not impose any LI constraint!".format(type))
			dpartition = dtcom.find_partition_Dartmouth(subDG, adjmx, init_partition, thresholdSize, None)
		else:
			arcpy.AddMessage("And the ouput {} impose LI constraint!".format(type))
			dpartition = dtcom.find_partition_Dartmouth(subDG, adjmx, init_partition, thresholdSize, miniLocalIndex)

arcpy.AddMessage("The {} delineates {} {}".format(delineateMethod, len(set(dpartition.values())), type))

inputPolyFL = dict2featurecls(dpartition, inputIDField, addComID, inputEdgePath, inputPolyFL) # Join Partition to Input Polygon Layer

arcpy.CopyFeatures_management(inputPolyFL, inputPolyFL2)
inputPolyFL2FdNames= [field.name for field in arcpy.ListFields(inputPolyFL2)]
if NumDZoneIDField in inputPolyFL2FdNames:
	arcpy.DeleteField_management(inputPolyFL2, NumDZoneIDField)

if len(HRR_NumDNode) > 0 and type == "HRRs":
	inputPolyFL2 = dict2featurecls(HRR_NumDNode, inputIDField, NumDZoneIDField, inputEdgePath, inputPolyFL2)

if outputType == "Dissolved":
	if inputPopField != "":
		if len(HRR_NumDNode) > 0 and type == "HRRs":
			arcpy.Dissolve_management(inputPolyFL2, outputHSAs, addComID, [[inputIDField,"COUNT"], [inputPopField, "SUM"], [NumDZoneIDField, "SUM"]])
			arcpy.AlterField_management(outputHSAs, "SUM_"+ NumDZoneIDField, NumDZoneIDField, NumDZoneIDField) # change the SUM_Num_DZoneID to Num_DZoneID
		else:
			arcpy.Dissolve_management(inputPolyFL2, outputHSAs, addComID, [[inputIDField,"COUNT"], [inputPopField, "SUM"]])
		arcpy.AlterField_management(outputHSAs, "SUM_"+ inputPopField, inputPopField, inputPopField) # change the SUM_POPU to POPU	
							
	else:
		if len(HRR_NumDNode) > 0 and type == "HRRs":
			arcpy.Dissolve_management(inputPolyFL2, outputHSAs, addComID, [[inputIDField,"COUNT"], [NumDZoneIDField, "SUM"]])
			arcpy.AlterField_management(outputHSAs, "SUM_"+ NumDZoneIDField, NumDZoneIDField, NumDZoneIDField) # change the SUM_Num_DZoneID to Num_DZoneID
		else:
			arcpy.Dissolve_management(inputPolyFL2, outputHSAs, addComID, [[inputIDField,"COUNT"]])
	
	# Calculate indices: LI, avgTime, the number of destination nodes
	comindices = dtcom.calculateindices(DG, dpartition, init_partition)
	outputHSAs = dict2featurecls2(comindices, addComID, ["LI", "EstTime", "Num_D" + inputIDField, "Compactness"], inputEdgePath, outputHSAs)

	if edgeDistField == "": # if no travel time between nodes, delete the travel time field
		arcpy.DeleteField_management(outputHSAs, "EstTime")
elif outputType == "Not dissolved":
	arcpy.CopyFeatures_management(inputPolyFL, outputHSAs) # only output the one with one more field HSAID or HRRID
else:
	arcpy.AddError("Please select a correct output Type!")

arcpy.DeleteField_management(inputPolyFL, addComID)
if arcpy.Exists(inputPolyFL2):
	arcpy.Delete_management(inputPolyFL2)