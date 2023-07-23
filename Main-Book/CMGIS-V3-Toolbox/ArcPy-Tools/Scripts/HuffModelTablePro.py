# Import system modules
import arcpy

# Get input parameters

customerFL = arcpy.GetParameterAsText(0)
customerIDField = arcpy.GetParameterAsText(1)
customerSizeField = arcpy.GetParameterAsText(2)
facilityFL = arcpy.GetParameterAsText(3)
facilityIDField = arcpy.GetParameterAsText(4)
facilitySizeField = arcpy.GetParameterAsText(5)
distanceTable = arcpy.GetParameterAsText(6)
distanceCustomerID = arcpy.GetParameterAsText(7)
distanceFacilityID = arcpy.GetParameterAsText(8)
distanceValue = arcpy.GetParameterAsText(9)
distanceDecayFunc = arcpy.GetParameterAsText(10)
distanceDecayCoeffbeta = arcpy.GetParameterAsText(11)
distanceDecayCoeffalpha = arcpy.GetParameterAsText(12)
scaleFactor = arcpy.GetParameterAsText(13)
outputType = arcpy.GetParameterAsText(14)
outputTable = arcpy.GetParameterAsText(15)
outputHSAs = arcpy.GetParameterAsText(16)

# Delete intermediate data if they already exist
if arcpy.Exists("Temp_Dist_Table"):
	arcpy.Delete_management("Temp_Dist_Table")

if arcpy.Exists("Temp_Sum_Potent"):
	arcpy.Delete_management("Temp_Sum_Potent")

if arcpy.Exists("Temp_Facility_MaxPotent"):
	arcpy.Delete_management("Temp_Facility_MaxPotent")

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


# Copy distance table into geodatabase
if ".dbf" in distanceTable:
	arcpy.CopyRows_management(distanceTable, "Temp_Dist_Table")
	distanceTable = "Temp_Dist_Table"

# Check if distance table already has the following fields:
# facility size, potential, probability, estflow
fieldList = arcpy.ListFields(distanceTable)
for field in fieldList:
	if field.name == customerSizeField:
		arcpy.DeleteField_management(distanceTable, [customerSizeField])
	elif field.name == facilitySizeField:
		arcpy.DeleteField_management(distanceTable, [facilitySizeField])
	elif field.name == "Potential":
		arcpy.DeleteField_management(distanceTable, ["Potential"])
	elif field.name == "Weight":
		arcpy.DeleteField_management(distanceTable, ["Weight"])
	elif field.name == "MAX_Potential":
		arcpy.DeleteField_management(distanceTable, ["MAX_Potential"])
	elif field.name == "SUM_Potential":
		arcpy.DeleteField_management(distanceTable, ["SUM_Potential"])
	elif field.name == "Probability":
		arcpy.DeleteField_management(distanceTable, ["Probability"])
	elif field.name == "EstFlow":
		arcpy.DeleteField_management(distanceTable, ["EstFlow"])


# Check if customer layer already has the following fields:
# facility id, probability, estflow
fieldList = arcpy.ListFields(customerFL)
for field in fieldList:
	if field.name == distanceFacilityID:
		arcpy.DeleteField_management(customerFL, [distanceFacilityID])
	elif field.name == "Probability":
		arcpy.DeleteField_management(customerFL, ["Probability"])
	elif field.name == "EstFlow":
		arcpy.DeleteField_management(customerFL, ["EstFlow"])
	elif field.name == facilitySizeField:
		arcpy.DeleteField_management(customerFL, facilitySizeField)

# Join customer layer to distance table, append customer size field
arcpy.JoinField_management(distanceTable, distanceCustomerID, customerFL, customerIDField, [customerSizeField])

# Join facility layer to distance table, append facility size field
arcpy.JoinField_management(distanceTable, distanceFacilityID, facilityFL, facilityIDField, [facilitySizeField])

# Add a field for potential of each facility to every customer
arcpy.AddField_management(distanceTable, "Potential", "DOUBLE")

# Calculate facility potential field with distance decay
if distanceDecayFunc == "Power":
	# Check if any distance is 0
	if ".mdb" in distanceTable:
		arcpy.MakeTableView_management(distanceTable, "Temp_DistanceView", "[" + distanceValue + "] = 0")
	else:
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

# Calculate probability of customer visiting a facility
arcpy.AddField_management(distanceTable, "Probability", "DOUBLE")
arcpy.CalculateField_management(distanceTable, "Probability", "!Potential! / !SUM_Potential!", "PYTHON_9.3")
arcpy.AddField_management(distanceTable, "EstFlow", "LONG")
arcpy.CalculateField_management(distanceTable, "EstFlow", "!Probability! * !" + customerSizeField + "!" , "PYTHON3")

# output estimated flows in a table
if outputType == "Predicted OD flows":
	arcpy.TableSelect_analysis(distanceTable, outputTable)
# output preliminary HSAs
if outputType == "Dissolved preliminary HSAs":
	# Select records from distance table with facility potential equal to maximum potential
	if ".mdb" in distanceTable:
		arcpy.TableSelect_analysis(distanceTable, "Temp_Facility_MaxPotent", "[Potential] = [MAX_Potential]")
	else:
		arcpy.TableSelect_analysis(distanceTable, "Temp_Facility_MaxPotent", '"Potential" = "MAX_Potential"')
	# Join maximum facility potential table to customer layer, append facility id
	arcpy.JoinField_management(customerFL, customerIDField, "Temp_Facility_MaxPotent", distanceCustomerID, [distanceFacilityID])
	#arcpy.JoinField_management(customerFL, customerIDField, facilityFL, facilityIDField, [facilitySizeField])
	# Dissolve the customer layer to HSAs by facility id
	arcpy.Dissolve_management(customerFL, outputHSAs, [distanceFacilityID], [[customerSizeField,"SUM"], [distanceCustomerID, "COUNT"]])
	# Delete joined facility id field from customer layer
	arcpy.DeleteField_management(customerFL, [distanceFacilityID])
		
# Cleanup intermediate data
arcpy.Delete_management("Temp_Dist_Table")
arcpy.Delete_management("Temp_Sum_Potent")
arcpy.Delete_management("Temp_Facility_MaxPotent")
