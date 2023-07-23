# Import system modules
import arcpy

# Get input parameters
customerFL = arcpy.GetParameterAsText(0)
customerIDField = arcpy.GetParameterAsText(1)
customerSizeField = arcpy.GetParameterAsText(2)
facilityFL = arcpy.GetParameterAsText(3)
facilityIDField = arcpy.GetParameterAsText(4)
facilitySizeField = arcpy.GetParameterAsText(5)
distanceType = arcpy.GetParameterAsText(6)
distanceDecayFunc = arcpy.GetParameterAsText(7)
distanceDecayCoeffbeta = arcpy.GetParameterAsText(8)
distanceDecayCoeffalpha = arcpy.GetParameterAsText(9)
unitConversion = arcpy.GetParameterAsText(10)
scaleFactor = arcpy.GetParameterAsText(11)
outputType = arcpy.GetParameterAsText(12)
outputTable = arcpy.GetParameterAsText(13)
outputHSAs = arcpy.GetParameterAsText(14)

distanceTable = "Temp_ODMatrix"
distanceCustomerID = "IN_FID"
distanceFacilityID = "NEAR_FID"
distanceValue = "NEAR_DIST"

# Delete intermediate data if they already exist
if arcpy.Exists("Temp_ODMatrix"):
	arcpy.Delete_management("Temp_ODMatrix")

if arcpy.Exists("Temp_Sum_Potent"):
	arcpy.Delete_management("Temp_Sum_Potent")

if arcpy.Exists("Temp_Facility_MaxPotent"):
	arcpy.Delete_management("Temp_Facility_MaxPotent")

if arcpy.Exists("Temp_Customer_Pt"):
	arcpy.Delete_management("Temp_Customer_Pt")

if arcpy.Exists("Temp_Facility_Pt"):
	arcpy.Delete_management("Temp_Facility_Pt")

if arcpy.Exists(outputTable):
	arcpy.Delete_management(outputTable)
	
if arcpy.Exists(outputHSAs):
	arcpy.Delete_management(outputHSAs)

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
if distanceDecayCoeffbeta != "":
	if distanceDecayFunc == "Power":
		powerDecayCoeffbeta = distanceDecayCoeffbeta
	elif distanceDecayFunc == "Exponential":
		exponentialDecayCoeffbeta = distanceDecayCoeffbeta
	elif distanceDecayFunc == "Square-root exponential":
		sqrtexponentialDecayCoeffbeta = distanceDecayCoeffbeta
	elif distanceDecayFunc == "Gaussian":
		gaussianDecayCoeffbeta = distanceDecayCoeffbeta
	elif distanceDecayFunc == "Log-normal":
		lognormalDecayCoeffbeta = distanceDecayCoeffbeta
	elif distanceDecayFunc == "Log-logistic":
		loglogisticDecayCoeffbeta = distanceDecayCoeffbeta
		if distanceDecayCoeffalpha != "":
			loglogisticDecayCoeffalpha = distanceDecayCoeffalpha
		else:
			arcpy.AddWarning("There is no alpha specified for Log-logistic distance decay function and the default alpha = {} is used!".fomrat(loglogisticDecayCoeffalpha))
	elif distanceDecayFunc == "Compound power-exponential":
		comppowerexponentialDecayCoeffbeta = distanceDecayCoeffbeta
		if distanceDecayCoeffalpha != "":
			comppowerexponentialDecayCoeffalpha = distanceDecayCoeffalpha
		else:
			arcpy.AddWarning("There is no alpha specified for Compound power-exponential distance decay function and the default alpha = {} is used!".fomrat(comppowerexponentialDecayCoeffalpha))
	else:
		arcpy.AddMessage("There is no distance decay function listed!")
else:
	arcpy.AddMessage("The default coefficients are used!")

# Check feature layer type and convert it to centroid layer if necessary
customer_pt = customerFL
desc = arcpy.Describe(customerFL)
if desc.shapeType == "Polygon":
	customer_pt = "Temp_Customer_Pt"
	arcpy.FeatureToPoint_management(customerFL, customer_pt, "INSIDE")
	# Need to change customer id to NEW_FID for shapefile, which is FID + 1
	'''
	if customerIDField == "OBJECTID":
		fieldList = arcpy.ListFields(customerFL)
		for field in fieldList:
			if field.name == "NEW_FID":
				arcpy.DeleteField_management(customerFL, ["NEW_FID"])
		arcpy.AddField_management(customerFL, "NEW_FID", "LONG")
		arcpy.CalculateField_management(customerFL, "NEW_FID", "!OBJECTID! + 1", "PYTHON3")
		customerIDField = "NEW_FID"
	'''
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
	'''
	if facilityIDField == "OBJECTID":
		fieldList = arcpy.ListFields(facilityFL)
		for field in fieldList:
			if field.name == "NEW_FID":
				arcpy.DeleteField_management(facilityFL, ["NEW_FID"])
		arcpy.AddField_management(facilityFL, "NEW_FID", "LONG")
		arcpy.CalculateField_management(facilityFL, "NEW_FID", "!OBJECTID! + 1", "PYTHON3")
		facilityIDField = "NEW_FID"
	'''

# Calculate Euclidean distance from customer to facility
# The PointDistance_analysis has been deprecated in ArcGIS Pro
#arcpy.PointDistance_analysis(customer_pt, facility_pt, distanceTable)
# output IN_FID (=ObjectID), NEAR_FID (=ObjectID), NEAR_DIST
if distanceType == "Euclidean": # output unit is the identical the input features
	arcpy.GenerateNearTable_analysis(customer_pt, facility_pt, distanceTable, method = "PLANAR", closest = "ALL")
if distanceType == "Geodesic": # output unit is meter
	arcpy.GenerateNearTable_analysis(customer_pt, facility_pt, distanceTable, method = "GEODESIC", closest = "ALL")

# Scale distance value with unit conversion factor
expression0 = "!" + distanceValue + "! * " + unitConversion
arcpy.CalculateField_management(distanceTable, distanceValue, expression0, "PYTHON3")

# Cleanup temporary data
if arcpy.Exists("Temp_Customer_Pt"):
	arcpy.Delete_management("Temp_Customer_Pt")
if arcpy.Exists("Temp_Facility_Pt"):
	arcpy.Delete_management("Temp_Facility_Pt")

# Check if customer layer already has the facility id field
fieldList = arcpy.ListFields(customerFL)
for field in fieldList:
	if field.name == distanceCustomerID:
		arcpy.DeleteField_management(customerFL, [distanceCustomerID])
	elif field.name == distanceFacilityID:
		arcpy.DeleteField_management(customerFL, [distanceFacilityID])
	elif field.name == distanceValue:
		arcpy.DeleteField_management(customerFL, [distanceValue])
	elif field.name == "Probabilit":
		arcpy.DeleteField_management(customerFL, ["Probabilit"])

# Join customer layer to distance table, append customer size field
arcpy.JoinField_management(distanceTable, distanceCustomerID, customerFL, customerIDField, [customerSizeField])

# Join facility layer to distance table, append facility size field
arcpy.JoinField_management(distanceTable, distanceFacilityID, facilityFL, facilityIDField, [facilitySizeField])

# Add a field for potential of each facility to every customer
arcpy.AddField_management(distanceTable, "Potential", "DOUBLE")

# Calculate facility potential field with distance decay
if distanceDecayFunc == "Power":
	# Check if any distance is 0
	arcpy.MakeTableView_management(distanceTable, "Temp_DistanceView", '"' + distanceValue + '" = 0')
	cnt = int(arcpy.GetCount_management("Temp_DistanceView").getOutput(0))
	if cnt != 0:
		arcpy.AddWarning("Warning: distance between {0} customer-facility pair(s) is 0!\nA weight of 0 is assigned to such pair!".format(cnt))
	# Add a field for distance-decay based weights, calculate this field with codeblock
	arcpy.AddField_management(distanceTable, "Weight", "DOUBLE")
	weightExpression = "calculateWeight(!" + distanceValue + "!)"
	codeblock = """def calculateWeight(distance):
		if distance == 0:
			return 0
		else:
			return distance ** ((-1) * """ + powerDecayCoeffbeta + """) * """ + scaleFactor
	arcpy.CalculateField_management(distanceTable, "Weight", weightExpression, "PYTHON3", codeblock)
	# Calculate facility potential by applying weight to facility size
	expression1 = "!" + facilitySizeField + "! * !Weight!"
	arcpy.CalculateField_management(distanceTable, "Potential", expression1, "PYTHON3")
	# Delete weight field from distance table
	arcpy.DeleteField_management(distanceTable, ["Weight"])
elif distanceDecayFunc == "Exponential":
	expression1 = scaleFactor + " * !" + facilitySizeField + "! * math.exp((-1) * !" + distanceValue + "! * " + exponentialDecayCoeffbeta + ")"
	arcpy.CalculateField_management(distanceTable, "Potential", expression1, "PYTHON3")
elif distanceDecayFunc == "Square-root exponential":
	expression1 = scaleFactor + " * !" + facilitySizeField + "! * math.exp((-1) * math.sqrt(!" + distanceValue + "!) * " + sqrtexponentialDecayCoeffbeta + ")"
	arcpy.CalculateField_management(distanceTable, "Potential", expression1, "PYTHON3")
elif distanceDecayFunc == "Gaussian":
	expression1 = scaleFactor + " * !" + facilitySizeField + "! / (math.sqrt(2 * math.pi) * " + gaussianDecayCoeffbeta + ") * math.exp((-0.5) * !" + distanceValue + "! ** 2 / " + gaussianDecayCoeffbeta + " ** 2)"
	arcpy.CalculateField_management(distanceTable, "Potential", expression1, "PYTHON3")
elif distanceDecayFunc == "Log-normal":	
	#expression1 = scaleFactor + " * !" + facilitySizeField + "! * math.exp((-1) * (math.log(!" + distanceValue + "!)) ** 2 * " + lognormalDecayCoeffbeta + ")"
	#arcpy.CalculateField_management(distanceTable, "Potential", expression1, "PYTHON3")
	# Check if any distance is 0
	arcpy.MakeTableView_management(distanceTable, "Temp_DistanceView", '"' + distanceValue + '" = 0')
	cnt = int(arcpy.GetCount_management("Temp_DistanceView").getOutput(0))
	if cnt != 0:
		arcpy.AddWarning("Warning: distance between {0} customer-facility pair(s) is 0!\nA weight of 0 is assigned to such pair!".format(cnt))
	# Add a field for distance-decay based weights, calculate this field with codeblock
	arcpy.AddField_management(distanceTable, "Weight", "DOUBLE")
	weightExpression = "calculateWeight(!" + distanceValue + "!)"
	codeblock = """def calculateWeight(distance):
		if distance == 0:
			return 0
		else:				
			return math.exp((-1) * """ + lognormalDecayCoeffbeta + """ * (math.log(distance)) ** 2) * """ + scaleFactor
	arcpy.CalculateField_management(distanceTable, "Weight", weightExpression, "PYTHON3", codeblock)
	# Calculate facility potential by applying weight to facility size
	expression1 = "!" + facilitySizeField + "! * !Weight!"
	arcpy.CalculateField_management(distanceTable, "Potential", expression1, "PYTHON3")
	# Delete weight field from distance table
	arcpy.DeleteField_management(distanceTable, ["Weight"])
elif distanceDecayFunc == "Log-logistic": # suppose r=1
	expression1 = scaleFactor + " * !" + facilitySizeField + "! / (1 + (!" + distanceValue + "! /" + loglogisticDecayCoeffalpha + ") ** " + loglogisticDecayCoeffbeta + ")" 
	arcpy.CalculateField_management(distanceTable, "Potential", expression1, "PYTHON3")
else: # Compound power-exponential
	expression1 = scaleFactor + " * !" + facilitySizeField + "! * math.exp((-1) * " + comppowerexponentialDecayCoeffalpha + "* !" + distanceValue + "! ** " + comppowerexponentialDecayCoeffbeta + ")" 
	arcpy.CalculateField_management(distanceTable, "Potential", expression1, "PYTHON3")


# Summarize maximum and total facility potentials by each customer
arcpy.Statistics_analysis(distanceTable, "Temp_Sum_Potent", [["Potential", "MAX"], ["Potential", "SUM"]], distanceCustomerID)

# Join summary table to distance table, append maximum and total potentials
arcpy.JoinField_management(distanceTable, distanceCustomerID,
						   "Temp_Sum_Potent", distanceCustomerID,
						   ["MAX_Potential", "SUM_Potential"])

# Calculate probability and estimated flows of customer visiting a facility
arcpy.AddField_management(distanceTable, "Probabilit", "DOUBLE")
arcpy.CalculateField_management(distanceTable, "Probabilit", "!Potential! / !SUM_Potential!", "PYTHON3")
arcpy.AddField_management(distanceTable, "EstFlow", "LONG")
arcpy.CalculateField_management(distanceTable, "EstFlow", "!Probabilit! * !" + customerSizeField + "!" , "PYTHON3")

# output estimated flows in a table
if outputType == "Predicted OD flows":
	arcpy.TableSelect_analysis(distanceTable, outputTable)
# output preliminary HSAs
if outputType == "Dissolved preliminary HSAs":
	# Select records from distance table with facility potential equal to maximum potential
	arcpy.TableSelect_analysis(distanceTable, "Temp_Facility_MaxPotent", '"Potential" = "MAX_Potential"')
	# Join maximum facility potential table to customer layer, append facility id
	arcpy.JoinField_management(customerFL, customerIDField, "Temp_Facility_MaxPotent", distanceCustomerID, [distanceFacilityID])
	# Dissolve the customer layer to HSAs by facility id
	arcpy.Dissolve_management(customerFL, outputHSAs, [distanceFacilityID], [[customerSizeField,"SUM"], [distanceCustomerID, "COUNT"]])
	# Delete joined facility id field from customer layer
	arcpy.DeleteField_management(customerFL, [distanceFacilityID])

# Cleanup intermediate data
arcpy.Delete_management("Temp_ODMatrix")
arcpy.Delete_management("Temp_Sum_Potent")
arcpy.Delete_management("Temp_Facility_MaxPotent")
