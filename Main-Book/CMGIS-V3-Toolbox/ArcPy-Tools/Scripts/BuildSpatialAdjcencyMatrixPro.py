# **************************************************************************************************************
# Date: December 6, 2020
# Author: Changzhen Wang
# Description: 1) This script is used to construct a spatial adjacency matrix by queen continguity. 
# 			   2) It also accounts for isolated islands that are separated from mainland by physical barriers, 
# 			      such as river, ocean, mountain, wildlife management areas, or parks, etc.
#              3) It can accommodate users' input adjacency matrix that is manually created by the physical road.
#				  The input adjusted adjacency matrix should have two columns (no headers) to represent originID and destinationID. 
# 				  The originID, destinationID must be identical to the Unique ID Field of the layer for HSAs delineation.
# ***************************************************************************************************************

# Import system modules
import arcpy
import numpy as np
import csv
from scipy.sparse import csr_matrix
from scipy.sparse.csgraph import reverse_cuthill_mckee
from scipy.spatial import distance
import scipy.sparse
from osgeo import ogr
import os

# Get input parameters
inputPolyFL = arcpy.GetParameterAsText(0)
polyIDField = arcpy.GetParameterAsText(1)
isAdjusted = arcpy.GetParameterAsText(2)
inputAdjustedFile = arcpy.GetParameterAsText(3)
outputSpAdjFile = arcpy.GetParameterAsText(4)

# It is a very important and efficient function to construct cluster/island
# Algorithm reference: http://raphael.candelier.fr/?blog=Adj2cluster
def adj2cluster(A):
	# symmetrize adjacency matrix
	S = A + A.transpose()
	#print(S)
	# transform dense matrix into sparse matrix
	G_sparse = csr_matrix(S)
	#print(G_sparse)
	# Reverse Cuthill-McKee ordering
	r = np.flip(reverse_cuthill_mckee(G_sparse))
	#print(r)
	# Get the clusters
	mlabel = dict()
	clusterNum = 0
	mlabel[clusterNum] = [r[0]]
	for i in range(1, len(r)):
		if S[mlabel[clusterNum], r[i]].any():
			c = np.append(mlabel[clusterNum],r[i])
			mlabel[clusterNum] = c
		else:
			clusterNum = clusterNum + 1
			mlabel[clusterNum] = [r[i]]
	return mlabel

# This function doesn't work as the .svm can't convert to table
def adjacency_matrix(inputfile, inputID, append_idsfile):
	node_dict =dict()
	n = arcpy.GetCount_management(inputfile)
	spWeightMtrx = "Temp_Adj_Mtr1.swm"
	outputTempTab = "Temp_Adj_Table1.dbf"
	if arcpy.Exists(spWeightMtrx):
		arcpy.Delete_management(spWeightMtrx)
	if arcpy.Exists(outputTempTab):
		arcpy.Delete_management(outputTempTab)	
	arcpy.GenerateSpatialWeightsMatrix_stats(inputfile, inputID, spWeightMtrx, "CONTIGUITY_EDGES_CORNERS", None, None, None, None, None, None, None, None, None, None)
	arcpy.ConvertSpatialWeightsMatrixtoTable_stats(spWeightMtrx, outputTempTab) # this function does not work because of it should be 32 bits
	geometries = arcpy.CopyFeatures_management(inputfile, arcpy.Geometry())

# Generate spatial adjacent matrix either by Shapefile or Geodatabase
def adjacency_matrix_2(filename, gdbname, ID, type, append_idsfile):
	node_dict = dict()
	if type == "Shapefile":
		shp_driver = ogr.GetDriverByName("ESRI Shapefile")
		vector = shp_driver.Open(filename, 0)  # 0-read, 1-writable
		layer = vector.GetLayer(0)  # shapefile always only has only 1 layer
		n = layer.GetFeatureCount()
		adj = np.zeros((n, n), dtype=int) # define an zero matrix
		for i in range(n):  # in shapefile, required id starts from 0
			feature1 = layer.GetFeature(i)
			ID_i = feature1.GetField(ID)
			geom1 = feature1.GetGeometryRef()
			node_dict[ID_i-1] = [geom1.Centroid().GetX(), geom1.Centroid().GetY()]
			layer.SetSpatialFilter(geom1)
			if n > 1:
				for feature2 in layer:
					ID_j = feature2.GetField(ID) # FID could not be identified! 12-6-2020
					if ID_i != ID_j:
						adj[ID_i-1][ID_j-1] = adj[ID_j-1][ID_i-1] = 1
			else:
				pass
	elif type =="OpenFileGDB":
		gdb_driver = ogr.GetDriverByName("OpenFileGDB")
		gdb_ds = gdb_driver.Open(gdbname, 0)
		layer = gdb_ds.GetLayerByName(filename)
		n = layer.GetFeatureCount()
		adj = np.zeros((n, n), dtype=int) # define an zero matrix
		for i in range(1, (n+1)):  # in feature class, required id starts from 1
			feature1 = layer.GetFeature(i)
			ID_i = feature1.GetField(ID)
			geom1 = feature1.GetGeometryRef()
			node_dict[ID_i-1] = [geom1.Centroid().GetX(),geom1.Centroid().GetY()]
			layer.SetSpatialFilter(geom1)
			if n > 1:
				for feature2 in layer:
					ID_j = feature2.GetField(ID)  # pay attention to field type of ID
					if ID_i != ID_j:
						adj[ID_i-1][ID_j-1] = adj[ID_j-1][ID_i-1] = 1
			else:
				pass
	mlabel = adj2cluster(adj)
	#mlabel= sorted(mlabel0.items(), key = lambda item: len(item[1]))
	while len(mlabel) > 1:
		node_values = dict()
		nml = len(mlabel.keys())
		for i in range(nml):
			node_values[i] = [node_dict[k] for k in mlabel[i] if k in node_dict]
		mtx1 = [[np.inf]*nml for i in range(nml)]
		mtx2 = [[np.inf]*nml for i in range(nml)]
		mtx3 = [[np.inf]*nml for i in range(nml)]
		for i in range(nml):
			nodelist = mlabel[i]
			for j in range(i):
				nbr_list = mlabel[j]
				dist_list = distance.cdist(node_values[i], node_values[j], 'euclidean')

				min_dist = np.min(dist_list)
				a = np.where(dist_list == min_dist)  
				mtx1[i][j] = mtx1[j][i]= min_dist
				mtx2[i][j] = mtx2[j][i] = nodelist[a[0][0]]
				mtx3[i][j] = mtx3[j][i] = nbr_list[a[1][0]]
			#b = np.where(mtx1 == np.min(mtx1))
		#col_minindex = np.argmin(mtx1, axis=0)[:-1] # axis =0 indicates column
		col_minindex = np.argmin(mtx1, axis=0)
		for k in range(nml-1):
			rowidx = col_minindex[k]
			ID1 = mtx2[rowidx][k]
			ID2 = mtx3[rowidx][k]
			adj[ID1][ID2] = adj[ID2][ID1] = 1
		mlabel = adj2cluster(adj)
		nml = len(mlabel.keys())                             
	#print (len(mlabel))
	# this part is used to add spatial adjacency matrix from physical road
	if append_idsfile != "":
		append_ids = csv.reader(open(append_idsfile, "r", newline=""))
		#next(append_ids, None) #skip the header
		for row in append_ids:
			O_ID = int(row[0])-1
			D_ID = int(row[1])-1
			if (adj[O_ID][D_ID] == 0) or (adj[D_ID][O_ID] == 0):
				adj[O_ID][D_ID] = adj[D_ID][O_ID] = 1
	return csr_matrix(adj)


def write_adjacency2csv(adjmx0, outputfile, ID1, ID2):
	adjmx = adjmx0.todense().nonzero()
	rowlsts = adjmx[0]
	collsts = adjmx[1]
	#print(len(rowlsts), len(collsts))
	writer = csv.writer(open(outputfile, 'w', newline=''))
	O_ID = "O_" + ID1
	D_ID = "D_" + ID2
	OD_W = "weight" 
	writer.writerow([O_ID, D_ID, OD_W])
	for i in range(len(rowlsts)):       
		writer.writerow([rowlsts[i] + 1, collsts[i] + 1 , 1]) 


# Check the input Polygon Layer is Polygon, avoid Change the Filter in the Tool Properties to other types
desc = arcpy.Describe(inputPolyFL)
if desc.shapeType != "Polygon":
	arcpy.AddMessage("Please input a polygon layer!")
else:
	spAdj = None
	filename = desc.name
	filepath = desc.path
	if filename.endswith(".shp"):
		arcpy.AddMessage("You are processing Shapefile {0} under {1}".format(filename, filepath))
		filefullpath = filepath + "\\" + filename
		if isAdjusted == "No":
			spAdj = adjacency_matrix_2(filefullpath, "", polyIDField, "Shapefile", "")
		else:
			if inputAdjustedFile != "":
				spAdj = adjacency_matrix_2(filefullpath, "", polyIDField, "Shapefile", inputAdjustedFile)
			else:
				arcpy.AddMessage("If you want to adjust the spatial adjacent matrix by physical road, please input the adjusted matrix file!")
	else:
		if filepath.endswith(".gdb"):
			arcpy.AddMessage("You are processing Feature class {0} under {1}".format(filename, filepath))
			if isAdjusted == "No":
				spAdj = adjacency_matrix_2(filename, filepath, polyIDField, "OpenFileGDB" , "")
			else:
				if inputAdjustedFile != "":
					spAdj = adjacency_matrix_2(filename, filepath, polyIDField, "OpenFileGDB", inputAdjustedFile)
				else:
					arcpy.AddMessage("If you want to adjust the spatial adjacent matrix by physical road, please input the adjusted matrix file!")
		else:
			arcpy.AddMessage("Please make sure your input Polygon layer is in Geodatabase!")
	arcpy.AddMessage("Output file path {0}".format(outputSpAdjFile))

	# If output is not end with .npz, add .npz to the end
	if not outputSpAdjFile.endswith(".npz"):
		outputSpAdjFile = outputSpAdjFile + ".npz"
		arcpy.AddMessage("Output file {0} is not end with .npz, it is automatically appended!".format(outputSpAdjFile))
		
	scipy.sparse.save_npz(outputSpAdjFile, spAdj)

	# output the spatial adjacency matrix along with the npz file for checking
	outputfilePath = os.path.split(outputSpAdjFile)[0]
	outputfileName = os.path.split(outputSpAdjFile)[1][:-4]
	outputfilePathName = outputfilePath + "\\" + outputfileName + ".csv"
	write_adjacency2csv(spAdj, outputfilePathName, polyIDField, polyIDField)
	arcpy.AddMessage("The corresponding csv file {0} is also output!".format(outputfilePathName))




