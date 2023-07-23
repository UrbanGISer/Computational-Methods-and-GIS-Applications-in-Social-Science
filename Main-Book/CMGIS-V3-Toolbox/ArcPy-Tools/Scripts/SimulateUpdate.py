'''----------------------------------------------------------------------------------
 Tool Name:     Simulate by Probability
 Source Name:   SimulateProbability.py
 Version:       ArcGIS Pro 3.0
 Author:        Lingbo Liu.
 Required Argumuments:  An input feature class 
                      An input filed for probability or frequency
                        
 Description:   Simulate trip list according the given probability list.
----------------------------------------------------------------------------------'''


import arcpy
import numpy as np
import pandas as pd
import collections
import os

###################### Read Data ####################################
# Get input parameters
sim_table= arcpy.GetParameterAsText(0)
tableID= arcpy.GetParameterAsText(1)
prob_field = arcpy.GetParameterAsText(2)
sim_num = arcpy.GetParameterAsText(3)
#output_table = arcpy.GetParameterAsText(4)


###################### Simulation ####################################
with arcpy.da.SearchCursor(sim_table, [prob_field]) as rows:
	df = pd.DataFrame.from_records(data=rows,columns=rows.fields)

print(df.head(3))
df = df.rename(columns={prob_field:'InProb'})

# Build a list of probability with field name 'prob'
df['prob']=df.InProb/sum(df.InProb)

#TypeError: cannot perform reduce with flexible type
#df = df.astype(float)
num=int(sim_num)
# Random simulation according to probability
simlist=np.random.choice(range(0,df.shape[0]), num, p=df.prob).tolist()

# group by  frequency
ct = collections.Counter(simlist)

ct_df = pd.DataFrame.from_dict(ct , orient='index').reset_index()
ct_df = ct_df.rename(columns={0:'SimN'})

mg=pd.merge(df.reset_index(),ct_df, on='index',how='left')
mg['SimN'].fillna(0, inplace = True)
mg['OxID']= range(1,(df.shape[0]+1))
mg1=mg[['OxID','SimN']]
###################### Save ####################################
#this code for update input data, it takes 3 minuts
# Check if customer layer already has the facility id field
fieldList = arcpy.ListFields(sim_table)
for field in fieldList:
	if field.name == "SimN":
		arcpy.DeleteField_management(sim_table, ["SimN"])
# build temp file path
inputFilePath = os.path.split(sim_table)[0]
fileTemp = inputFilePath + "\\" + "Temp_File"
if arcpy.Exists(fileTemp):
	arcpy.Delete_management(fileTemp)

x = mg1.reset_index(drop = True)
z = np.rec.fromrecords(x.values, names=x.columns.tolist())
arcpy.da.NumPyArrayToTable(z, fileTemp)

arcpy.JoinField_management(sim_table, tableID, fileTemp, 'OxID', ['SimN'])
arcpy.DeleteField_management(sim_table, ["OxID"])
arcpy.Delete_management(fileTemp)

#mg1.to_csv(output_table,index=False)