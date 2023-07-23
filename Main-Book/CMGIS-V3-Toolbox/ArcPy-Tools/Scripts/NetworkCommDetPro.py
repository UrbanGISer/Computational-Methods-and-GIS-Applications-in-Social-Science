# **************************************************************************************************************
# Date: December 8, 2020
# Author: Changzhen Wang
# Description: 1) This script is used to delineate HSAs or HRRs by (Sc)Louvain or (Sc)Leiden method
# **************************************************************************************************************


# Import system modules
import arcpy
import os
import igraph as ig
import leidenalg as la
import numpy as np
import networkx as nx
import scipy
import NetworkRefineCommBySpAdjMtrPro as netrfcom


# need to match row and input field name
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
inputResolution = arcpy.GetParameterAsText(9)
imposeSpAdj = arcpy.GetParameterAsText(10)
inputSpAdjMatrix = arcpy.GetParameterAsText(11)
thresholdSize = arcpy.GetParameterAsText(12)
outputHSAs = arcpy.GetParameterAsText(13)


if inputResolution == "":
	arcpy.AddError("please input a number (>0) for Input Resolution")
else:
	inputResolution = float(inputResolution)

if thresholdSize == "":
	thresholdSize = None
elif thresholdSize == "1,000":
	thresholdSize = 1000
elif thresholdSize == "120,000":
	thresholdSize = 120000
else:
	thresholdSize = int(thresholdSize)


# Join partition results to feature class by nodeID and valueID
def dict2featurecls(part_dict, keyID, valueID, edgepath, featcls):
	fieldList = arcpy.ListFields(featcls)
	for field in fieldList:
		if field.name == valueID:
			arcpy.DeleteField_management(featcls, [valueID])
	part_dict1 = {(key+1) : (value+1) for key, value in part_dict.items()}
	partition_arr = np.array(list(part_dict1.items()), dtype=[(keyID, np.int64),(valueID, np.int64)])
	partitionTemp = edgepath + "\\" + "Temp_Parition_File"
	if arcpy.Exists(partitionTemp):
		arcpy.Delete_management(partitionTemp)
	arcpy.da.NumPyArrayToTable(partition_arr, partitionTemp)	
	arcpy.JoinField_management(featcls, keyID, partitionTemp, keyID, [valueID])
	arcpy.Delete_management(partitionTemp)
	return featcls


# Join partition results to feature class by keyID and multiple valueIDs
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
	part_list = [((key + 1), *value) for (key, value) in part_dict.items()]
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
		

G = nx.Graph()
DG = nx.DiGraph()
edges = []
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
			G.add_node(int(row[0])-1, id = int(row[0])-1, pop = int(row[1])) # the first item
			DG.add_node(int(row[0])-1, id = int(row[0])-1, pop = int(row[1]))
			init_partition[int(row[0])-1] = -1
else:
	with arcpy.da.SearchCursor(inputPolyFL, [inputIDField]) as cursor:
		for row in cursor:
			G.add_node(int(row[0])-1, id = int(row[0])-1, pop = 0) # the first item
			DG.add_node(int(row[0])-1, id = int(row[0])-1, pop = 0)
			init_partition[int(row[0])-1] = -1

# Save the number of destination ZIP codes within the HRR
HRR_NumDNode = {}
# Define field to save no. of dest nodes within HSA
NumDZoneIDField = "Num_DZoneID"

# read the table in the geodatabase to construct the edges of a graph
inputEdgePath = os.path.split(inputEdgeFile)[0]
inputEdgeFdNames= [field.name for field in arcpy.ListFields(inputEdgeFile)]
if edgeDistField != "":
	if NumDZoneIDField in inputEdgeFdNames and type == "HRRs":
		with arcpy.da.SearchCursor(inputEdgeFile, [edgeOrgIDField, edgeDestIDField, edgeFlowField, edgeDistField, NumDZoneIDField], """{} > 0""".format(edgeFlowField)) as cursor:
			for row in cursor:
				if G.has_edge(int(row[0])-1, int(row[1])-1):
					G[int(row[0])-1][int(row[1])-1]['weight'] += float(row[2])
					G[int(row[0])-1][int(row[1])-1]['estTime'] += float(row[3])
				else:
					G.add_edge(int(row[0])-1, int(row[1])-1, weight=float(row[2]), estTime= float(row[3]))
				DG.add_edge(int(row[0])-1, int(row[1])-1, weight=float(row[2]), estTime=float(row[3]))
				edges.append((int(row[0])-1, int(row[1])-1, float(row[2]), float(row[3])))
				destID = int(row[1])-1
				init_partition[destID] = destID

				if not bool(HRR_NumDNode.get(destID)):
					HRR_NumDNode[destID] = int(row[4]) - 1
	else:
		with arcpy.da.SearchCursor(inputEdgeFile, [edgeOrgIDField, edgeDestIDField, edgeFlowField, edgeDistField], """{} > 0""".format(edgeFlowField)) as cursor:
			for row in cursor:
				if G.has_edge(int(row[0])-1, int(row[1])-1):
					G[int(row[0])-1][int(row[1])-1]['weight'] += float(row[2])
					G[int(row[0])-1][int(row[1])-1]['estTime'] += float(row[3])
				else:
					G.add_edge(int(row[0])-1, int(row[1])-1, weight=float(row[2]), estTime= float(row[3]))
				DG.add_edge(int(row[0])-1, int(row[1])-1, weight=float(row[2]), estTime=float(row[3]))
				edges.append((int(row[0])-1, int(row[1])-1, float(row[2]), float(row[3])))
				destID = int(row[1])-1
				init_partition[destID] = destID
else:
	if NumDZoneIDField in inputEdgeFdNames and type == "HRRs":
		with arcpy.da.SearchCursor(inputEdgeFile, [edgeOrgIDField, edgeDestIDField, edgeFlowField, NumDZoneIDField], """{} > 0""".format(edgeFlowField)) as cursor:
			for row in cursor:
				if G.has_edge(int(row[0])-1, int(row[1])-1):
					G[int(row[0])-1][int(row[1])-1]['weight'] += float(row[2])
				else:
					G.add_edge(int(row[0])-1, int(row[1])-1, weight=float(row[2]), estTime= 0)
				DG.add_edge(int(row[0])-1, int(row[1])-1, weight=float(row[2]), estTime=0)
				edges.append((int(row[0])-1, int(row[1])-1, float(row[2]), 0))
				destID = int(row[1])-1
				init_partition[destID] = destID

				if not bool(HRR_NumDNode.get(destID)):
					HRR_NumDNode[destID] = int(row[3]) - 1
	else:
		with arcpy.da.SearchCursor(inputEdgeFile, [edgeOrgIDField, edgeDestIDField, edgeFlowField], """{} > 0""".format(edgeFlowField)) as cursor:
			for row in cursor:
				if G.has_edge(int(row[0])-1, int(row[1])-1):
					G[int(row[0])-1][int(row[1])-1]['weight'] += float(row[2])
				else:
					G.add_edge(int(row[0])-1, int(row[1])-1, weight=float(row[2]), estTime= 0)
				DG.add_edge(int(row[0])-1, int(row[1])-1, weight=float(row[2]), estTime=0)
				edges.append((int(row[0])-1, int(row[1])-1, float(row[2]), 0))
				destID = int(row[1])-1
				init_partition[destID] = destID

init_coms = len(set(v for k, v in init_partition.items() if k == v))
node_pop = dict(DG.nodes(data = 'pop'))
total_pop = sum(node_pop.values())
arcpy.AddMessage("Total number of nodes is {0}".format(G.number_of_nodes()))            
arcpy.AddMessage("Total number of flows is {0}".format(DG.number_of_edges())) 
arcpy.AddMessage("Total number of service volumes is {0}".format(G.size(weight = 'weight')))
if total_pop > 0:
	arcpy.AddMessage("Total population is {0}".format(total_pop))
arcpy.AddMessage("The number of destination ZIP codes is {0}".format(init_coms))

G1 = ig.Graph(directed = False)
G1.add_vertices(list(set(G.nodes)))
G1.vs["name"] = list(set(G.nodes))
G1.add_edges([(x, y) for (x, y, z, w) in edges])
G1.es['weight'] = [z for (x, y, z, w) in edges]
G1.es['estTime'] = [w for (x, y, z, w) in edges]

# Generate preliminary delineation results
prepartition = None
if delineateMethod == "ScLeiden":
	prepartition = la.find_partition(G1, la.RBConfigurationVertexPartition, None, G1.es["weight"], 20, 0, 1, resolution_parameter = inputResolution)
elif delineateMethod == "ScLouvain":
	optimiser = la.Optimiser()
	prepartition = la.RBConfigurationVertexPartition(G1, None, G1.es["weight"], None, inputResolution)
	prepartition_agg = prepartition.aggregate_partition()
	while optimiser.move_nodes(prepartition_agg) > 0:
		prepartition.from_coarse_partition(prepartition_agg)
		prepartition_agg = prepartition_agg.aggregate_partition()
else:
	# to avoid users change the parameters in tool properties
	arcpy.AddError("Please ensure the selected delineation Method is ScLeiden or ScLouvain!")


# Save preliminary results to a dict
partition_dict = {}
for name, membership in zip(G1.vs["name"], prepartition.membership):
	if imposeSpAdj == "No" or imposeSpAdj == "Yes":
		partition_dict[int(name)] = int(membership)
	else:
		arcpy.AddError("Please ensure whether you wanted to impose spatial adjacency matrix or not!")


# Calculate modularity and the number of HSAs
numComs = len(set(prepartition.membership))  
premodularity = ig.Graph.modularity(G1, prepartition.membership, "weight")
arcpy.AddMessage("The {} method generates {} {} with the modularity of {:.4f} when resolution = {}". format(delineateMethod[2:], numComs, type, premodularity, inputResolution))


# Delete outputHSAs if exists
if arcpy.Exists(outputHSAs):
	arcpy.Delete_management(outputHSAs)

outputHSAPath = os.path.split(outputHSAs)[0]
outputHSAName = os.path.split(outputHSAs)[1]
outputHSAs2 = ""

# use to process donut HSAs and island HSAs in the final result
HSAs_neighHSAs = "HSAs_NeighHSAs"
HSAs_neighHSAs_Freq = "HSAs_NeighHSAs_Freq"
inputPolyFL2 = inputPolyFL + "_CPTemp"

if outputHSAPath.endswith(".gdb"):
	outputHSAs2 = outputHSAPath + "\\" + outputHSAName + "NotDis"
	HSAs_neighHSAs = outputHSAPath + "\\" + HSAs_neighHSAs
	HSAs_neighHSAs_Freq = outputHSAPath + "\\" + HSAs_neighHSAs_Freq
if outputHSAName.endswith(".shp"):
	outputHSAs2 = outputHSAPath + "\\" + outputHSAName[:-4] + "NotDis.shp"
	HSAs_neighHSAs = outputHSAPath + "\\" + HSAs_neighHSAs + ".dbf"
	HSAs_neighHSAs_Freq = outputHSAPath + "\\" + HSAs_neighHSAs_Freq + ".dbf"

if arcpy.Exists(outputHSAs2):
	arcpy.Delete_management(outputHSAs2) 
if arcpy.Exists(HSAs_neighHSAs):
	arcpy.Delete_management(HSAs_neighHSAs)
if arcpy.Exists(HSAs_neighHSAs_Freq):
	arcpy.Delete_management(HSAs_neighHSAs_Freq)
if arcpy.Exists(inputPolyFL2):
	arcpy.Delete_management(inputPolyFL2)

arcpy.CopyFeatures_management(inputPolyFL, inputPolyFL2)
inputPolyFL2FdNames= [field.name for field in arcpy.ListFields(inputPolyFL2)]
if NumDZoneIDField in inputPolyFL2FdNames:
	arcpy.DeleteField_management(inputPolyFL2, NumDZoneIDField)

# Whether impose spatial adjacency matrix or not
if imposeSpAdj == "No":	
	inputPolyFL2 = dict2featurecls(partition_dict, inputIDField, addComID, inputEdgePath, inputPolyFL2)
	if len(HRR_NumDNode) > 0 and type == "HRRs":		
		inputPolyFL2 = dict2featurecls(HRR_NumDNode, inputIDField, NumDZoneIDField, inputEdgePath, inputPolyFL2)

	arcpy.CopyFeatures_management(inputPolyFL2, outputHSAs2)

	if inputPopField != "":
		if len(HRR_NumDNode) > 0 and type == "HRRs":
			arcpy.Dissolve_management(inputPolyFL2, outputHSAs, addComID, [[inputIDField,"COUNT"], [inputPopField, "SUM"], [NumDZoneIDField, "SUM"]])
			arcpy.AlterField_management(outputHSAs, "SUM_"+ NumDZoneIDField, NumDZoneIDField, NumDZoneIDField) # change the SUM_Num_DZoneID to Num_DZoneID
		else:
			arcpy.Dissolve_management(inputPolyFL2, outputHSAs, addComID, [[inputIDField,"COUNT"], [inputPopField, "SUM"]])
		arcpy.AlterField_management(outputHSAs, "SUM_"+ inputPopField, inputPopField, inputPopField)
	else:
		if len(HRR_NumDNode) > 0 and type == "HRRs":
			arcpy.Dissolve_management(inputPolyFL2, outputHSAs, addComID, [[inputIDField,"COUNT"], [NumDZoneIDField, "SUM"]])
			arcpy.AlterField_management(outputHSAs, "SUM_"+ NumDZoneIDField, NumDZoneIDField, NumDZoneIDField) # change the SUM_Num_DZoneID to Num_DZoneID
		else:
			arcpy.Dissolve_management(inputPolyFL2, outputHSAs, addComID, [[inputIDField,"COUNT"]])

	if arcpy.Exists(inputPolyFL2):
		arcpy.Delete_management(inputPolyFL2)
	
	comindices = netrfcom.calculateindices(DG, partition_dict, init_partition)
	outputHSAs = dict2featurecls2(comindices, addComID, ["LI", "EstTime", "Num_D" + inputIDField, "Compactness"], inputEdgePath, outputHSAs)

	if edgeDistField == "": # if no travel time between nodes, delete the travel time field
		arcpy.DeleteField_management(outputHSAs, "EstTime")

elif imposeSpAdj == "Yes":
	# Load spatial adjacency matrix
	adjmx = None
	if inputSpAdjMatrix == "":
		arcpy.AddError("Please copy the path of .npz file to the Input Spatial Adjacency Matrix File.")
	elif not inputSpAdjMatrix.endswith(".npz"):
		arcpy.AddError("Please ensure the input path file is end with .npz format")
	else:
		adjmx = scipy.sparse.load_npz(inputSpAdjMatrix).todense().nonzero() #the row and column represents the unique IDs of input polygon layer
	rfpartition = netrfcom.refined_partition_network(G, adjmx, thresholdSize, inputResolution, partition_dict, init_partition)
	inputPolyFL2 = dict2featurecls(rfpartition, inputIDField, addComID, inputEdgePath, inputPolyFL2) # Join Partition to Input Polygon Layer
	
	if len(HRR_NumDNode) > 0 and type == "HRRs":		
		inputPolyFL2 = dict2featurecls(HRR_NumDNode, inputIDField, NumDZoneIDField, inputEdgePath, inputPolyFL2)
	# Output the non-dissolved HSAs /HRRs
	arcpy.CopyFeatures_management(inputPolyFL2, outputHSAs2)

	# Dissolve the ZIP codes to HSAs/ HRRs
	if inputPopField != "":
		if len(HRR_NumDNode) > 0 and type == "HRRs":
			arcpy.AddMessage(sum(HRR_NumDNode.values()))
			arcpy.Dissolve_management(inputPolyFL2, outputHSAs, addComID, [[inputIDField,"COUNT"], [inputPopField, "SUM"], [NumDZoneIDField, "SUM"]])
			arcpy.AlterField_management(outputHSAs, "SUM_"+ NumDZoneIDField, NumDZoneIDField, NumDZoneIDField) # change the SUM_Num_DZoneID to Num_DZoneID
		else:
			arcpy.Dissolve_management(inputPolyFL2, outputHSAs, addComID, [[inputIDField,"COUNT"], [inputPopField, "SUM"]])
		arcpy.AlterField_management(outputHSAs, "SUM_"+ inputPopField, inputPopField, inputPopField)
	else:
		if len(HRR_NumDNode) > 0 and type == "HRRs":
			arcpy.Dissolve_management(inputPolyFL2, outputHSAs, addComID, [[inputIDField,"COUNT"], [NumDZoneIDField, "SUM"]])
			arcpy.AlterField_management(outputHSAs, "SUM_"+ NumDZoneIDField, NumDZoneIDField, NumDZoneIDField) # change the SUM_Num_DZoneID to Num_DZoneID
		else:
			arcpy.Dissolve_management(inputPolyFL2, outputHSAs, addComID, [[inputIDField,"COUNT"]])

	#arcpy.DeleteField_management(inputPolyFL, addComID)

	arcpy.PolygonNeighbors_analysis(outputHSAs, HSAs_neighHSAs, addComID)
	arcpy.Statistics_analysis(HSAs_neighHSAs, HSAs_neighHSAs_Freq, [["src_{}".format(addComID), "COUNT"]], "src_{}".format(addComID))
	arcpy.JoinField_management(HSAs_neighHSAs, "src_{}".format(addComID), HSAs_neighHSAs_Freq, "src_{}".format(addComID), ["FREQUENCY"])

	ComID_2_ComID = dict()
	__MIN = 0.00001
	with arcpy.da.SearchCursor(HSAs_neighHSAs, ["src_{}".format(addComID), "nbr_{}".format(addComID), "LENGTH"], '"FREQUENCY" = 1') as cursor:
		for row in cursor:
			with arcpy.da.SearchCursor(outputHSAs, [addComID, "SHAPE@LENGTH"]) as cursor2:
				for row2 in cursor2:
					diff = abs(float(row[2]) - float(row2[1]))
					if diff <= __MIN:
						ComID_2_ComID[row[0] - 1]= row[1] - 1

	arcpy.Delete_management(HSAs_neighHSAs)
	arcpy.Delete_management(HSAs_neighHSAs_Freq)

	if len(ComID_2_ComID) > 0:
		rfpartition2 = dict()
		for key, value in rfpartition.items():
			if value in ComID_2_ComID.keys():
				rfpartition2[key] = ComID_2_ComID[value]
			else:
				rfpartition2[key] = value
		arcpy.Delete_management(outputHSAs)
		com2node = netrfcom.partition_com2node(rfpartition2)
		rfpartition, com2node, com_old_vs_com_new = netrfcom.__renumber(rfpartition2, com2node)

		arcpy.Delete_management(inputPolyFL2)
		inputPolyFL2 = dict2featurecls(rfpartition, inputIDField, addComID, inputEdgePath, inputPolyFL2)
		
		if len(HRR_NumDNode) > 0 and type == "HRRs":
			inputPolyFL2 = dict2featurecls(HRR_NumDNode, inputIDField, NumDZoneIDField, inputEdgePath, inputPolyFL2)
		
		arcpy.Delete_management(outputHSAs2)
		arcpy.CopyFeatures_management(inputPolyFL2, outputHSAs2)

		if inputPopField != "":
			if len(HRR_NumDNode) > 0 and type == "HRRs":
				arcpy.Dissolve_management(inputPolyFL2, outputHSAs, addComID, [[inputIDField,"COUNT"], [inputPopField, "SUM"], [NumDZoneIDField, "SUM"]])
				arcpy.AlterField_management(outputHSAs, "SUM_"+ NumDZoneIDField, NumDZoneIDField, NumDZoneIDField) # change the SUM_Num_DZoneID to Num_DZoneID
			else:
				arcpy.Dissolve_management(inputPolyFL2, outputHSAs, addComID, [[inputIDField,"COUNT"], [inputPopField, "SUM"]])
			arcpy.AlterField_management(outputHSAs, "SUM_"+ inputPopField, inputPopField, inputPopField)
		else:
			if len(HRR_NumDNode) > 0 and type == "HRRs":
				arcpy.Dissolve_management(inputPolyFL2, outputHSAs, addComID, [[inputIDField,"COUNT"], [NumDZoneIDField, "SUM"]])
				arcpy.AlterField_management(outputHSAs, "SUM_"+ NumDZoneIDField, NumDZoneIDField, NumDZoneIDField) # change the SUM_Num_DZoneID to Num_DZoneID
			else:
				arcpy.Dissolve_management(inputPolyFL2, outputHSAs, addComID, [[inputIDField,"COUNT"]])
		arcpy.DeleteField_management(inputPolyFL2, addComID)

	if arcpy.Exists(inputPolyFL2):
		arcpy.Delete_management(inputPolyFL2)

	rfmodularity = netrfcom.modularity(rfpartition, G)
	numComs = len(set(rfpartition.values()))
	arcpy.AddMessage("The {} method generates {} {} with the modularity of {:.4f} when resolution = {}". format(delineateMethod, numComs, type, rfmodularity, inputResolution))		

	# Calculate indices: LI and avgTime
	comindices = netrfcom.calculateindices(DG, rfpartition, init_partition)
	outputHSAs = dict2featurecls2(comindices, addComID, ["LI", "EstTime", "Num_D" + inputIDField, "Compactness"], inputEdgePath, outputHSAs)

	if edgeDistField == "": # if no travel time between nodes, delete the travel time field
		arcpy.DeleteField_management(outputHSAs, "EstTime")
else:
	arcpy.AddError("Please select Yes or No in the Impose Spatial Continguity!") 

