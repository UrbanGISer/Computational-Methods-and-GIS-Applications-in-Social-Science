'''----------------------------------------------------------------------------------
 Tool Name:     Repeat Row by Column Value
 Source Name:   RepeatRow.py
 Version:       ArcGIS Pro 3.0
 Author:        Lingbo Liu
 Required Argumuments:  An input feature class 
                        An input filed for repeat times
                        
 Description:   Expand O-D pairs in to single rows.
----------------------------------------------------------------------------------'''


import arcpy
import numpy as np
import pandas as pd

###################### Read Data ####################################
# Get input parameters
sim_table= arcpy.GetParameterAsText(0)
key_field = arcpy.GetParameterAsText(1)
count_field = arcpy.GetParameterAsText(2)
output_table = arcpy.GetParameterAsText(3)


###################### Simulation ####################################
allfiled=key_field.split(";")
allfiled.append(count_field)

with arcpy.da.SearchCursor(sim_table, allfiled) as rows:
	df = pd.DataFrame.from_records(data=rows,columns=rows.fields)


df = df.rename(columns={count_field:'times'})
#df["times"] = df["times"].astype(str).astype(int)
df["times"] = df["times"].astype(int)

dfnew=df.loc[df.index.repeat(df.times)]
dfnew["Count"]=1

###################### Save ####################################

#dfnew.to_csv(output_table,index=False)
#pandas dataframe to numpy to feature
x = dfnew.reset_index(drop = True)
z = np.rec.fromrecords(x.values, names=x.columns.tolist())
arcpy.da.NumPyArrayToTable(z, output_table)