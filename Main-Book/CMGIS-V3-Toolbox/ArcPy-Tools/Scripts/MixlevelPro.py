#*****************************************************************************************
# Embedded version, with functions of Func_ClusteringOrder.py, Func_OneLevelClustering.py, Func_TackleIsolation.py
# all added to the top of the code.
#
#This script take the input with two levels of geographic unit, e.g. tract and county,
#and assign mixed-level cluster type to each basic unit, e.g. tract based on multiple constraints
#such as population (>=20000) and total cancer count (>=15).
#The results are ClusType = 1,2, 3 or 4
#Where 1 = low level clusters (tract)
#      2 = mid level clusters (county)
#      3 = high level clusters (multi-county)
#      4 = mixed low and high level clusters (tract and county)
#
# Suggested citation for the mixed-level regionalization (MLR) method:
# Mu, L., Wang, F., Chen, V. W., & Wu, X. (2015). A place-oriented, mixed-level regionalization method for
# constructing geographic areas in health data dissemination and analysis. Annals of the Association of American Geographers,
# 105(1), 48-66. doi:10.1080/00045608.2014.968910.
#
#  Programmed by Mu, Lan (mulan@uga.edu)
#  copyright: 2012-
#  Revised ArcGIS Pro Version 2022: Lingbo Liu
#*****************************************************************************************

#*****************************************************************************************
#This script is modifed to use a funciton to perform all it used to do, so it can be called
# from another Python scrip.
#
#Description: 
#      This scripts implements the spatial order command from arcinfo workstation.
#      It asigns order values to features based on their x,y location. The net effect is
#      spatially sorting input features such that features that are close have similar
#      order values
#Inputs:
#      1)Input feature class
#      2)Input spatial order field name to create the new double field to store order values
#      
#Outputs:
#      1) Input feature class as derived output
#*****************************************************************************************

#prints a GP message and writes the same message to a log file
def PrintAndWriteMessages(msg,severity=0):
    if severity == 0:
        arcpy.AddMessage(msg)
    elif severity == 1:
        arcpy.AddWarning(msg)
    elif severity == 2:
        arcpy.AddError(msg)

#return the fractional part from a double number
def GetFractionalPart(dbl):
    return dbl - math.floor(dbl)

#return the peano curve coordinate for a given x,y value
def Peano(x,y,k):
    if (k == 0 or (x==1 and y==1)):
        return 0.5
    if x <= 0.5:       
        if y <= 0.5:
            quad=0
        else: 
            quad=1
    elif y <= 0.5:
        quad = 3
    else: 
        quad = 2
    subpos = Peano(2 * abs(x - 0.5), 2 * abs(y - 0.5), k-1)
    
    if(quad == 1 or quad == 3):
        subpos = 1 - subpos
    
    return GetFractionalPart((quad + subpos - 0.5)/4.0)

# Add a field if it does not exist in the input feature class (or layer).
def AddAField(theFC, theField, fldType):
    class GetOutOfDef( Exception ):
        pass  
    try:
        flds = arcpy.ListFields(theFC, theField)
        for fld in flds:
            if fld:
                raise GetOutOfDef
        arcpy.AddField_management(theFC, theField, fldType)
    except GetOutOfDef:
        pass

# Delete a field if it exists in the input feature class (or layer).
def DelAField(theFC, theField):
    flds = arcpy.ListFields(theFC, theField)
    for fld in flds:
        if fld:
            arcpy.DeleteField_management(theFC, theField)

def ClusteringOrder(inp_fc, PeanoOrder_fld, AttriOrder_fld, NAttriOrder_fld, Order_fld, Attri_list, Wt_list, Pct_Peano):
    
    #Import the modules and create the geoprocessor
    import arcpy, sys,os, traceback,math
    #gp = arcgisscripting.create(9.3)

    #Use the lowest available license
    for product in ['Engine','ArcView', 'ArcEditor', 'EngineGeoDB','ArcInfo', 'ArcServer']:
        if arcpy.CheckProduct(product).lower() == 'available':
            arcpy.SetProduct(product)
            break
    arcpy.env.overwriteOutput = True
    try:
        Pct_Attri = 100 - Pct_Peano

        #First Validate the field name
        PeanoOrder_fld = arcpy.ValidateFieldName(PeanoOrder_fld,os.path.dirname(inp_fc))
        AttriOrder_fld = arcpy.ValidateFieldName(AttriOrder_fld,os.path.dirname(inp_fc))
        NAttriOrder_fld = arcpy.ValidateFieldName(NAttriOrder_fld,os.path.dirname(inp_fc))
        Order_fld = arcpy.ValidateFieldName(Order_fld,os.path.dirname(inp_fc))  
        
        #Add the double field
        result = AddAField(inp_fc,PeanoOrder_fld,"DOUBLE")
        
        result = AddAField(inp_fc,AttriOrder_fld,"DOUBLE")

        result = AddAField(inp_fc,NAttriOrder_fld,"DOUBLE")
        
        result = AddAField(inp_fc,Order_fld,"DOUBLE")

        #Get the extent for the feature class
        desc = arcpy.Describe(inp_fc)
        arcpy.env.workspace = desc.path

        #Just in case, conver the input feature class layer to feature class, so the extent will be adjusted to the right values.
        #Reason:  a feature class layer's extent always points back to the reference feature class, in spite of the fact that the layer is
        # only a subset of all the features.
        theRealFC = "theRealFC"
        arcpy.FeatureClassToFeatureClass_conversion(inp_fc, arcpy.env.workspace, theRealFC)
        if desc.datasetType == 'ShapeFile':
            desc2 = arcpy.Describe(theRealFC + ".shp") 
        if desc.datasetType == 'FeatureClass':
            desc2 = arcpy.Describe(theRealFC)
        
        extent  = desc2.Extent
        xmin = extent.XMin
        ymin = extent.YMin
        xmax = extent.XMax
        ymax = extent.YMax
        
        #compute some constants to scale the coordinates to unit square before calling Peano
        dx = xmax - xmin
        dy = ymax - ymin
        if dx >= dy:
            offsetx = 0.0
            offsety = (1.0 - dy / dx)/ 2.0
            scale = dx
        else:
            offsetx = (1.0 - dx / dy)/ 2.0
            offsety = 0.0
            scale = dy
        
        #If the input features are lines or polygons get their centroids
        useCentroids = False
        if desc2.ShapeType.lower() in ['polyline','polygon']:
            useCentroids = True
            
        #Get each point and compute it's peano curve coordinate and store it back
        #Get an update cursor
        rows = arcpy.da.UpdateCursor(inp_fc,["SHAPE@",PeanoOrder_fld])
        #rows.reset()
        #row = next(rows)
        for row in rows:
            #Get the X,Y coordinate for each feature
            if useCentroids:
                pnt = row[0].centroid        
            else:
                pnt = row[0]

            unitx = (pnt.X - xmin) / scale + offsetx
            unity = (pnt.Y - ymin) / scale + offsety

            peanoPos = Peano(unitx, unity, 32)
            #row.setValue(PeanoOrder_fld,peanoPos)
            row[1] = peanoPos
            rows.updateRow(row)
            arcpy.SetProgressorPosition()
            #row = next(rows)
        del rows
        arcpy.Delete_management(theRealFC)

        #Calculate weighted aggregated attribute scores, scale the values to 0-1, then calcualte the order values. 
        #Calculate weighted aggregated attribute scores and find min and max values.
        AggrWt = 0
        for wt in Wt_list.split(';'):
            AggrWt = AggrWt + float(wt)
            
        curFlds = Attri_list.split(';')
        curFlds.append(AttriOrder_fld)
        rows = arcpy.da.UpdateCursor(inp_fc, curFlds)
        #rows.reset()
        j = 0
        #row = next(rows)
        for row in rows:
            AggrScore = 0
            i = 0
            for attri in Attri_list.split(';'):
                wt = Wt_list.split(';')[i]
                AggrScore = AggrScore + float(row[curFlds.index(attri)]) * float(wt)       
                i = i + 1

            AggrScore = AggrScore / AggrWt
            row[curFlds.index(AttriOrder_fld)] = AggrScore
            rows.updateRow(row)
            
            if j == 0:
                minScore = AggrScore
                maxScore = AggrScore
            else:
                if AggrScore < minScore:
                    minScore = AggrScore
                if AggrScore > maxScore:
                    maxScore = AggrScore

            arcpy.SetProgressorPosition()                    
            #row = next(rows)
            j = j + 1
        del rows

        #Normalize weighted aggregated attribute scores then calculate the final order values.
        curFlds = [AttriOrder_fld,NAttriOrder_fld, PeanoOrder_fld,Order_fld]
        rows = arcpy.da.UpdateCursor(inp_fc, curFlds)
        #rows.reset()
        #row = next(rows)
        for row in rows:
            NAggrScore = (row[0] - minScore) / (maxScore - minScore)
            row[1] = NAggrScore
            CombOrder = (row[2] * Pct_Peano + row[1] * Pct_Attri) / 100
            row[3] = CombOrder
            
            rows.updateRow(row)
            arcpy.SetProgressorPosition()
        del rows
        
        #Set the derived output when the script is done
        #Useful to show the output when the tool is used in model builder
        #arcpy.SetParameterAsText(8,inp_fc)
    except:
        # Return any python specific errors as well as any errors from the geoprocessor
        tb = sys.exc_info()[2]
        tbinfo = traceback.format_tb(tb)[0]
        pymsg = "PYTHON ERRORS:\nTraceback Info:\n" + tbinfo + "\nError Info:\n    " + \
                str(sys.exc_type)+ ": " + str(sys.exc_value) + "\n"
        PrintAndWriteMessages(pymsg,2)

        msgs = "GP ERRORS:\n" + arcpy.GetMessages(2) + "\n"
        PrintAndWriteMessages(msgs,2)
#End of funciton ClusteringOrder        

def OneLevelClustering(inp_fc, order_item, constraint_list, capacity_list, out_class_item, isolate):
 
    #Only run the function if there are more than 2 features.
    if int(arcpy.GetCount_management(inp_fc).getOutput(0)) >=2:
        
        try:           
            #parse the input parameter strings to lists
            ConstraintList = constraint_list.split(';')
            CapacityList = capacity_list.split(';')            
            #Get the OID field name
            desc = arcpy.Describe(inp_fc)
            oid_fld = desc.OIDFieldName
            # Start of new codes
            # Get the path and full name of the input feature class
            arcpy.env.workspace = desc.path
            fullInp_fc = desc.BaseName + "." + desc.Extension 
            out_class_item = arcpy.ValidateFieldName(out_class_item, os.path.dirname(inp_fc))
            
            #Add the new field. If the field already exists, AddField has a warning
            result = AddAField(inp_fc, out_class_item,"Double")
            #Validate the isolatte field name for the workspace
            isolate = arcpy.ValidateFieldName(isolate,os.path.dirname(inp_fc))
            #Add the new field. If the field already exists, AddField has a warning
            result = AddAField(inp_fc,isolate,"Short")
            arcpy.CalculateField_management(inp_fc,isolate, 0,"PYTHON_9.3","#")
            #Calculate the Class field
            #Get an update cursor such that values are in ascending order based on OrderItem field

            # start of new codes

        # new code, put all needed fields in a list then convert to a string
            fields = [order_item]
            for constraint in ConstraintList:
                fields.append(constraint)
            fields.append("TheID")
            fields.append("included")
            fieldsStr = ";".join(fields)
            # use a field to record whether a recorded in included in the cluster or not
            #Add the new field. If the field already exists, AddField has a warning
            result = AddAField(inp_fc,"included","Short")
            arcpy.CalculateField_management(inp_fc,"included",0,"PYTHON_9.3","#")
            where_clause = "included = 0" 
            # end of new codes 
            # Generate a spatial weight matrix
            # *.swm file cannot be saved in a geodatabase 
            # arcpy.env.overWriteOutput = True
            if desc.datasetType == 'FeatureClass':
                upDir = os.path.dirname(desc.path)
                swmTab = upDir + "\\" + "tmpswm.swm"
                dbfTab = upDir + "\\" + "tmpswm.dbf"
            else:
                swmTab = "tmpswm.swm"
                dbfTab = "tmpswm.dbf"
            if arcpy.Exists(swmTab):
                arcpy.Delete_management(swmTab)
            if arcpy.Exists(dbfTab):
                arcpy.Delete_management(dbfTab)
            #GenerateSpatialWeightsMatrix_stats requires an input of a feature class (inp_fc might be a layer)
            tmpFC = "tmpFC"
            fcLoc = arcpy.env.workspace
            arcpy.FeatureClassToFeatureClass_conversion(inp_fc, fcLoc, tmpFC)         
            #fullpathFC = arcpy.env.workspace + "\\" + tmpFC
            if desc.datasetType == 'ShapeFile':
                tmpFC = tmpFC + ".shp"
            arcpy.GenerateSpatialWeightsMatrix_stats(tmpFC, "TheID", swmTab, "CONTIGUITY_EDGES_ONLY", "EUCLIDEAN", "1", "#", "#","NO_STANDARDIZATION") 
            arcpy.ConvertSpatialWeightsMatrixtoTable_stats(swmTab, dbfTab)

            # end of new codes

            tot_rec = int(arcpy.GetCount_management(inp_fc).getOutput(0))
          
            #rows = arcpy.UpdateCursor(inp_fc,where_clause,"",fields,order_item + " A")               
            rows = arcpy.da.UpdateCursor(inp_fc,fields, where_clause, sql_clause=(None, 'ORDER BY '+ order_item +' ASC'))
            #rows.reset()
            row = next(rows)
            
        ##    acc_capacity1 = 0
        ##    acc_capacity2 = 0
            acc_capacity = []
            for i in range(len(ConstraintList)):
                acc_capacity.append(0)
                
            class_value = 1
            count = 0
            roundCount = 0
            newRoundCount = tot_rec
            orphanIteration = 0
            #Store the class value and the list of OIDs in that class in a dict.
            classValueDict = {class_value : []}

            ######################################
            ## Major round of clustering starts ##
            ######################################

            PrintAndWriteMessages("\n"+ "Total records = " + str(tot_rec), 1)
            
            while row:
            #for row in  rows:
                #Should use TheID instead of thisoid to search, because sometimes TheID != thisoid
                #thisoid = row.GetValue(oid_fld)
                thisTheID = row[fields.index("TheID")]
                    
                count += 1 #counter for records being processed
                roundCount += 1 #counter for this round of search
                cur_capacity = []
                
                # Get capacity values according to thisTheID
                for i in range(len(ConstraintList)):
                    #msg = "The current constraint is: " +  str(ConstraintList[i]) 
                    #PrintAndWriteMessages(msg,1)
                    cur_capacity.append(row[fields.index(ConstraintList[i])])
                
                ### Start of new codes
                if count == 1:
                    isRookNbr = True
                    row[fields.index("included")] = 1
                    rows.updateRow(row)                   
                        
                # Starting from the second polygon, make sure the new polygon is a rook neighbor (share boundary)
                # of one of the current cluster members. This will get rid of the "jumping" problem.
                if count > 1:
                    isRookNbr = False
                    
                    # the first member in each new cluster is always in
                    if len(classValueDict[class_value]) == 0:
                        row[fields.index("included")] = 1
                        rows.updateRow(row)
                        isRookNbr = True                            
                    else:
                        swmRows = arcpy.da.SearchCursor(dbfTab, ["TheID", "NID"], "TheID = " + str(thisTheID))
                        for swmRow in swmRows:

                            if str(swmRow[1]) in classValueDict[class_value]:

                                row[fields.index("included")] = 1
                                rows.updateRow(row)
                                isRookNbr = True
                                #jump out of for swmRows loop
                                break               
                            ##del swmRow
                            ##swmRow = swmnext(rows)
                        del swmRows
                ### End of new codes
                # new code, the "if included == 1" condition, to see if the caps are reached
                if row[fields.index("included")] == 1:
                    for i in range(len(ConstraintList)):
                        acc_capacity[i] += cur_capacity[i]           
                   
                    #start of new codes
                    #Accumulate the transfer item till it is equal or just greater than capacity    
                    #classValueDict[class_value].append(str(row[fields.index("TheID")]))
                    classValueDict[class_value].append(str(thisTheID))
                    # modified to take more constraints.
                    satisfyAll = 1
                    
                    for i in range(len(ConstraintList)):
                        if float(acc_capacity[i]) >= float(CapacityList[i]):
                            satisfyAll = satisfyAll * 1

                        else:
                            satisfyAll = satisfyAll * 0
                    if satisfyAll == 1 and count < tot_rec: 
                        for i in range(len(ConstraintList)):
                            acc_capacity[i] = 0
                        class_value += 1
                        classValueDict[class_value] = []
                    if count == tot_rec:
                        msg = ""      
                        PrintAndWriteMessages(msg,0)
                    #update cursor and the count for the next round of search
                    #rows = arcpy.UpdateCursor(inp_fc,where_clause,"",fields,order_item + " A")
                    rows = arcpy.da.UpdateCursor(inp_fc,fields, where_clause, sql_clause=(None, 'ORDER BY ' + order_item + ' ASC'))
                    #rows.reset()
                    
                    newRoundLyr = "newroundlyr"
                    #arcpy.MakeFeatureLayer_management(inp_fc,newRoundLyr,where_clause)
                    arcpy.MakeFeatureLayer_management(inp_fc,newRoundLyr,where_clause)
                    newRoundCount = int(arcpy.GetCount_management(newRoundLyr).getOutput(0))
                    arcpy.Delete_management(newRoundLyr)
                    roundCount = 0
                    
                # When included = 0
                else:
                    count -= 1


                # Now, deal with "orphans" - after the first sequential search by order value, isolated units left with included = 0
                # Decide how many iterations wanted in looking for "orphans".  Pull out all included = 0 and go back to for row in rows loop.
                
                #if thisoid == lastoid:
                if roundCount == newRoundCount:
                    orphanlyr = "orphanlyr"
                    #arcpy.MakeFeatureLayer_management(inp_fc,orphanlyr,where_clause)
                    arcpy.MakeFeatureLayer_management(inp_fc,orphanlyr,where_clause)
                    orphanCount = int(arcpy.GetCount_management(orphanlyr).getOutput(0))
                    arcpy.Delete_management(orphanlyr)         
                    #If the cursor is not empty, then start a new cluster
                    if orphanCount > 0:
                        orphanIteration += 1                        
                        #if the previous cluster is not empty, then start a new cluster.
                        if len(classValueDict[class_value]) > 0:  
                            for i in range(len(ConstraintList)):
                                acc_capacity[i] = 0                                
                            class_value += 1
                            classValueDict[class_value] = []                            
                        #rows = arcpy.UpdateCursor(inp_fc,where_clause,"",fields,order_item + " A")
                        rows = arcpy.da.UpdateCursor(inp_fc,fields, where_clause, sql_clause=(None, 'ORDER BY ' + order_item + ' ASC'))
                        try:
                            row = next(rows)
                        except:
                            break
                        newRoundLyr = "newroundlyr"
                        #arcpy.MakeFeatureLayer_management(inp_fc,newRoundLyr,where_clause)
                        arcpy.MakeFeatureLayer_management(inp_fc,newRoundLyr,where_clause)
                        newRoundCount = int(arcpy.GetCount_management(newRoundLyr).getOutput(0))
                        arcpy.Delete_management(newRoundLyr)
                        roundCount = 0
                    else:
                        #del row
                        try:
                            row = next(rows)
                        except:
                            break
                # keep moving if it's not the last OID
                else:    
                    #del row
                    try:
                        row = next(rows)
                    except:
                        break
                            
                    
        ##	#Update the progressor after processing 5% rows
           
            del rows
            ####################################
            ## Major round of clustering ends ##
            ####################################


                
            # Delete empty entries in the calue dictionary 
            for k in classValueDict.keys():
                if len(classValueDict[k]) == 0:
                    del classValueDict[k]   

            #Update the class values by first making a selection set and then using Calculate field.
            for k in classValueDict.keys():
                where_clause = '"TheID" IN (' + ','.join(classValueDict[k]) + ')'
                lyr = "inp_fc_lyr"
                arcpy.MakeFeatureLayer_management(inp_fc, lyr, where_clause)
                arcpy.CalculateField_management(lyr,out_class_item,k,"PYTHON_9.3","#")
                arcpy.Delete_management(lyr)
            ##        position += len(classValueDict[k]) /2.0
            ##        arcpy.SetProgressorPosition(position)

            # start of new codes
            # If a cluster has a subtotal < capacity, then merge it to an adjacent cluster to make it work.

            # Execute the Summary Statistics tool using the SUM option
            tmptb = "sum_tmp"
            if arcpy.Exists(tmptb):
                arcpy.Delete_management(tmptb)
            #modified to inrease flexibility of constraints.
            #StatFields = transfer_item1 + " SUM;" + transfer_item2 + " SUM"
            StatFields = ""
            for i in range(len(ConstraintList)):
                StatFields = StatFields + ConstraintList[i] + " SUM;"
            StatFields = StatFields.rstrip(';')
            arcpy.Statistics_analysis(inp_fc, tmptb, StatFields, out_class_item)

            # Get a list of fields from the new summary table.
            flds = arcpy.ListFields(tmptb)
            # Retrieve the field with the total sum value.
            fldcount = 0
            fldsum = []
            for fld in flds:
        ##        msg = "A field in the sum_tmp table is: " + fld.name
        ##        PrintAndWriteMessages(msg,0)
                if fld.name.__contains__("SUM_"):
                    #modified to increase flexibility of consrtaints.
                    fldsum.append(fld.name)
                    fldcount += 1

            # Open a Search Cursor using field name.
            #where_clause3 = fld1 + " < " + str(capacity1) + " or " + fld2 + " < " + str(capacity2)
            #modified to increase flexibility of consrtaints.
            where_clause3 = ""
            for i in range(fldcount):
                where_clause3 = where_clause3 + fldsum[i] + " < " + str(CapacityList[i]) + " or "
            where_clause3 = where_clause3.rstrip(' or ')

##            ### Monitoring intermediate results...
##            PrintAndWriteMessages("tmptb # of records = " + arcpy.GetCount_management(tmptb).getOutput(0),1)
##            PrintAndWriteMessages("tmptb fields: ", 1)
##            fldnames = []
##            for fld in flds:
##                fldnames.append(fld.name)
##            PrintAndWriteMessages(fldnames,1)
##
##            tmprows = arcpy.da.SearchCursor(tmptb,"*")
##            for tmprow in tmprows:
##                msg = str(tmprow[0]) + ", " + str(tmprow[1]) + ", " + str(tmprow[2]) + ", " + str(tmprow[3])+ ", " + str(tmprow[4])
##                PrintAndWriteMessages(msg, 1)
##            ### End of Monitoring...

            rowFlds = fldsum
            rowFlds.append(out_class_item)
            rows = arcpy.da.SearchCursor(tmptb, rowFlds, where_clause3)

            
            #Create an empty list to record whether a cluster needs to be merged or not.
            clusAdjList = []
            isoClus = []
            rowCount = 0
            
            for rowLocal in rows:
                rowCount += 1
                FoundIt = False
                foundCount = 0
                # which previous cluster is a rook neighbor to the unsatisfied cluster and has a min subtotal to reach capacities?
                currentClus = rowLocal[rowFlds.index(out_class_item)]
                theCluster = currentClus
                #PrintAndWriteMessages("currentClus=" + str(currentClus), 1)
                #PrintAndWriteMessages(classValueDict[currentClus], 1)
                where_clause4 = '"TheID" IN (' + ','.join(classValueDict[currentClus]) + ')'
                swmFlds = "NID"
                swmRows = arcpy.da.SearchCursor(dbfTab, swmFlds, where_clause4)
                #swmRows.reset()
                #swmRow = swmnext(rows)
                        
                for swmRow in swmRows:
                    for i in range(class_value, 0, -1):
                        #look for other clusters than itself
                        if i !=currentClus:             
                            if str(swmRow[swmFlds.index("NID")]) in classValueDict[i]:
                                # check to see if the combined capacities from the two clusters meet the criteria
                                # modified to increase constraints flexibility
                                capA = []
                                for j in range(len(ConstraintList)):
                                    capA.append(rowLocal[rowFlds.index(fldsum[j])])
                                    
                                foundFlds=fldsum
                                foundRows = arcpy.da.SearchCursor(tmptb, foundFlds, out_class_item + " = " + str(i))
                                #foundRows.reset()
                                foundRow = next(foundRows)
                                
                                # modified to increase constraints flexibility                               
                                capB = []                                
                                for j in range(len(ConstraintList)):
                                    capB.append(foundRow[foundFlds.index(fldsum[j])])
                                
                                del foundRow, foundRows

                                # modified to increase constraints flexibility
                                newcap = []
                                for j in range(len(ConstraintList)):
                                    newcap.append(capA[j] + capB[j])
                                    
                                newSatisfyAll = 1
                                for j in range(len(ConstraintList)):
                                    if float(newcap[j]) >= float(CapacityList[j]):
                                        newSatisfyAll = newSatisfyAll * 1
                                    else:
                                        newSatisfyAll = newSatisfyAll * 0

                                #if newcap1 >= capacity1 and newcap2 >= capacity2:
                                if newSatisfyAll == 1:        
                                    foundCount += 1
                                    if foundCount == 1:
                                        mincap = []
                                        for k in range(len(ConstraintList)):
                                            mincap.append(newcap[k])
                                        FoundIt = True
                                        theCluster = i
                                        break
                                    #elif foundCount > 1 and newcap1 < mincap1 and newcap2 < mincap2:
                                    elif foundCount > 1:
                                        minSatisfyAll = 1
                                        for j in range(len(ConstraintList)):
                                            if float(newcap[j]) < float(mincap[j]):
                                                minSatisfyAll = minSatisfyAll * 1
                                            else:
                                                minSatisfyAll = minSatisfyAll * 0

                                        if  minSatisfyAll == 1:    
                                            for k in range(len(ConstraintList)):
                                                mincap[k] = newcap[k]
                                            FoundIt = True
                                            theCluster = i
                                            break
                                    
                #del swmRow
                #swmRow = swmnext(rows)
                del swmRows
               
                #for every unmet cluster
        ##        msg = "The current searched cluster is " + str(currentClus) + ", and FoundIt = " + str(FoundIt)
        ##        PrintAndWriteMessages(msg,0)

                sumvalues = ""
                for i in range(len(ConstraintList)):
                    sumvalues = sumvalues + str(rowLocal[rowFlds.index(fldsum[i])]) + ", "
                sumvalues = sumvalues.rstrip(', ')
                    
                where_clause5 = out_class_item + " = " + str(currentClus)
                lyr2 = "unmetClus_lyr"
                arcpy.MakeFeatureLayer_management(inp_fc,lyr2,where_clause5)                    
                if FoundIt == True:              
                    #clusAdjust += 1
                    arcpy.CalculateField_management(lyr2,out_class_item, theCluster,"PYTHON_9.3","#")
                    arcpy.CalculateField_management(lyr2,isolate, 0,"PYTHON_9.3","#")

                else:
                    # previous approach
                    arcpy.CalculateField_management(lyr2,isolate, 1,"PYTHON_9.3","#")                    
                    # new approach to force the cluster membership by aggregating all isolated ones or
                    # by merging to the nearest met cluster.  Record it here and process it out side of the while loop.
                    isoClus.append(currentClus)                        
                arcpy.Delete_management(lyr2)    
                clusAdjList.append(FoundIt)

            del rows
            #after the merging,  adjust all involving cluster members.  Continue to use rows generated from where_clause3
            #where_clause3 = fld1 + " < " + str(capacity1) + " or " + fld2 + " < " + str(capacity2)
            rowFlds = out_class_item
            rows = arcpy.da.SearchCursor(tmptb, rowFlds,where_clause3)
            clusAdjust = 0
            listInd = 0

            
            for rowLocal in rows:
                
                if len(clusAdjList) > 0:
                    if clusAdjList[listInd] == True:
                        where_clause6 = out_class_item + " > " + str(rowLocal[rowFlds.index(out_class_item)] - clusAdjust)
                        lyr3 = "adjClus_lyr"
                        arcpy.MakeFeatureLayer_management(inp_fc,lyr3,where_clause6)
                        exp = "!" + out_class_item + "! - 1"
                        arcpy.CalculateField_management(lyr3,out_class_item, exp,"PYTHON_9.3","#")
                        arcpy.Delete_management(lyr3)
                        clusAdjust += 1
                        listInd += 1
            del rows
            
            arcpy.Delete_management(tmptb)
            arcpy.Delete_management(swmTab)
            arcpy.Delete_management(dbfTab)
            
            msg = "Clustering at this level is done. After adjustment, the number of subclusters is " \
                  + str(class_value - clusAdjust) + " including " + str(len(isoClus)) + \
                  " remained unsatisfied isolated clusters."
            arcpy.SetProgressorLabel(msg)
            PrintAndWriteMessages(msg,1)
            # end of new codes
            
            #Set the output to be same as input
            #Useful to show the output when the tool is used in model builder
            #arcpy.SetParameterAsText(6,inp_fc)
        except:
            # Return any python specific errors as well as any errors from the geoprocessor
            tb = sys.exc_info()[2]
            tbinfo = traceback.format_tb(tb)[0]
            pymsg = "PYTHON ERRORS:\nTraceback Info:\n" + tbinfo + "\nError Info:\n    " + \
                    str(sys.exc_type)+ ": " + str(sys.exc_value) + "\n"
            PrintAndWriteMessages(pymsg,2)

            msgs = "GP ERRORS:\n" + arcpy.GetMessages(2) + "\n"
            PrintAndWriteMessages(msgs,2)

    #If there are less than 2 features:
    else:
        msgs = "\nOne-level clustering function is not performed because the number of features is less than 2."
        PrintAndWriteMessages(msgs,1)
#End of funciton OneLevelClustering

###########################################################        
# A function  to tackle isolated cluster at different levels.   
# Add a new field of integrated cluster membership and dissolve all clusters at all levels.
# format: T.UUU.S, where T = cluster type (1,2,3), UUU = upper level code (e.g. county FIPS), and S = sub cluster ID 
# clusType = 1,2, 3 or 4
# Where 1 = low level clusters (tract)
#       2 = mid level clusters (county)
#       3 = high level clusters (multi-county)
#       4 = mixed low and high level clusters (tract and county)
###########################################################

def TackleIsolation(inp_fc, out_class_item, isolate, upper_ID, cluster_type, cluster, mixed_clusters, constraint_list, capacity_list):

    def MergeIsolated(tmpMixedClusFC, theClus, fldIso, revType):
        # General rules to tackle isolated clusters:
        # 1) First priority, the minimum-capacity cluster(e.g. population) from the adjacent non-isolated
        #    Don't use NEAR_ANALYSIS in this step, because it will randomly choose an adjacent one. 
        # 2) Second priority, if no adjacent cluster, find the nearest non-adjacent, non-isolated cluster

        lyr5 = "lyr5"; lyr6 = "lyr6"; lyr7 = "lyr7"
        
        # first priority
        where_clause = cluster + " = \'" + theClus + "\'"
        arcpy.MakeFeatureLayer_management(tmpMixedClusFC, lyr5, where_clause)
        where_clause = fldIso + " = 0"
        arcpy.MakeFeatureLayer_management(tmpMixedClusFC, lyr6, where_clause)
        arcpy.SelectLayerByLocation_management(lyr6, "SHARE_A_LINE_SEGMENT_WITH", lyr5)
        # lyr6 is non-isolated clusters adjacent to the isolated cluster

        rows2 = arcpy.da.SearchCursor(lyr6,[fldsCap[0],cluster])
        count = 0
        for row2 in rows2:
            count += 1
            if count == 1:
                minCap = row2[0]
                minClus = row2[1]
            elif row2[0] < minCap:
                minCap = row2[0]
                minClus = row2[1]
        del rows2

        # If find an adjecant cluster to merge, Go back to original file, change the previous cluster type and cluster membership
        if count >= 1:
            theNewClusID = minClus
            msg = "Find an adjacent min-capacity non-isolated cluster " + str(minClus) 
            PrintAndWriteMessages(msg, 1)
            
        # second priority
        else:
            # no adjacent cluster, find a nearest non isolated cluster
            # lyr5 has only one feature, corresponding to each row
            # lyr7 is non-isolated clusters
            where_clause = fldIso + " = 0"
            arcpy.MakeFeatureLayer_management(tmpMixedClusFC, lyr7, where_clause)

##            tmpRows1 = arcpy.da.SearchCursor(lyr5, oid_fld)
##            for tmpRow1 in tmpRows1:
##                PrintAndWriteMessages(tmpRow1[0],2)        
##            tmpRows2 = arcpy.da.SearchCursor(lyr7, oid_fld)
##            for tmpRow2 in tmpRows2:
##                PrintAndWriteMessages(tmpRow2[0],1)
            
            arcpy.Near_analysis(lyr5, lyr7)
            rows3 = arcpy.da.SearchCursor(lyr5, "NEAR_FID")
            row3 = next(rows3)
            nearFID = row3[0]
            where_clause = oid_fld + " = " + str(nearFID)
            #PrintAndWriteMessages(where_clause, 1)
            rows4 = arcpy.da.SearchCursor(lyr7, cluster, where_clause)
            row4 = next(rows4)
            theNewClusID = row4[0]
            del rows3, rows4
            msg = "Find a near non-adjacent non-isolated cluster " + str(theNewClusID) 
            PrintAndWriteMessages(msg, 1)
        
        # Update original file (before clustering). It will be re-dissolved later.
        # Keep the original membership, "out_class_item", and change the final membership, "cluster"
        where_clause = cluster + " = \'" + theClus + "\'"
        #PrintAndWriteMessages(where_clause, 1)
        arcpy.MakeFeatureLayer_management(inp_fc, lyr7, where_clause)
        arcpy.CalculateField_management(lyr7, cluster, theNewClusID, "PYTHON_9.3","#")
        
        PrintAndWriteMessages("Cluster " + str(theClus) + " is reassigned as " + str(theNewClusID) , 1)

        if revType != "":
            arcpy.CalculateField_management(lyr7, cluster_type, revType, "PYTHON_9.3","#")

        # End of "def MergeIsolated()"

    ##
    ##Main code starts here
    ##

    #Import the modules and create the geoprocessor
    import arcpy, sys,os, traceback,math

    #Use the lowest available license
    for product in ['Engine','ArcView', 'ArcEditor', 'EngineGeoDB','ArcInfo', 'ArcServer']:
        if arcpy.CheckProduct(product).lower() == 'available':
            arcpy.SetProduct(product)
            break
    arcpy.env.overwriteOutput = True

    # clean up previously created layers
    tmpLyrs = ["lyr","lyrIso", "lyrNonIso", "MixedClusLyr", "lyr1_2", "ly3", "lyr4", "lyr5", "lyr6", "lyr7"]  
    for tmpLyr in tmpLyrs:
        if arcpy.Exists(tmpLyr):
            arcpy.Delete_management(tmpLyr)
    try:
        #First Validate the field name
        out_class_item = arcpy.ValidateFieldName(out_class_item,os.path.dirname(inp_fc))
        isolate = arcpy.ValidateFieldName(isolate,os.path.dirname(inp_fc))
        upper_ID = arcpy.ValidateFieldName(upper_ID,os.path.dirname(inp_fc))
        cluster_type = arcpy.ValidateFieldName(cluster_type,os.path.dirname(inp_fc))
        cluster = arcpy.ValidateFieldName(cluster,os.path.dirname(inp_fc))
        mixed_clusters = arcpy.ValidateFieldName(mixed_clusters,os.path.dirname(inp_fc))

        #parse the input parameter strings to lists
        ConstraintList = constraint_list.split(';')
        CapacityList = capacity_list.split(';')

        isoCount1 = 0       # one level
        isoCount21 = 0      # mixed level, subsenario 1
        isoCount22 = 0      # mixed level, subsenario 2
        
        # Get the path and full name of the input feature class
        desc = arcpy.Describe(inp_fc)
        arcpy.env.workspace = desc.path
        fullInp_fc = desc.BaseName + "." + desc.Extension
       
        ### Additional codes to tackle isolated cluster at different levels.
        lyr = "lyr"
        arcpy.MakeFeatureLayer_management(inp_fc, lyr)
        AddAField(lyr, cluster, "Text")  

        # For mixed level clustering, Add a new field of integrated cluster membership.
        # format: T.UUU.S, where T = cluster type (1,2,3), UUU = upper level code (e.g. county FIPS), and S = sub cluster ID            
        if upper_ID != "" and cluster_type != "":
            lyr1_2 = "lyr1_2"
            where_clause = cluster_type + " <> 3"
            arcpy.MakeFeatureLayer_management(inp_fc, lyr1_2, where_clause)
            exp = "str(!" + cluster_type + "!) + '.'+ !" + upper_ID + "! + '.' + str(!" + out_class_item + "!)"
            arcpy.CalculateField_management(lyr1_2, cluster, exp, "PYTHON_9.3","#")
            arcpy.Delete_management(lyr1_2)

            # for multi-county (upper level unit) cluster, the upper_id in the cluster ID is assigned as special, e.g. 000 
            lyr3 = "lyr3"
            where_clause = cluster_type + " = 3"
            arcpy.MakeFeatureLayer_management(inp_fc, lyr3, where_clause)            
            if int(arcpy.GetCount_management(lyr3).getOutput(0)) > 0:
                rows = arcpy.da.SearchCursor(lyr3, upper_ID)
                count  = 0
                for row in rows:
                    count = count+1
                    firstUpperID = str(row[0])
                    if count == 1:
                        break
                del rows

                specialID = []
                for i in range(len(firstUpperID)):
                    specialID.append("0")
                specialID = "".join(specialID)

                #msg = "Multi-county's upper ID is assigned as " + specialID
                #PrintAndWriteMessages(msg, 1)
                exp = "str(!" + cluster_type + "!) + '.' + str('" + specialID + "') + '.' + str(!" + out_class_item + "!)"
                arcpy.CalculateField_management(lyr3, cluster, exp, "PYTHON_9.3","#")    
            arcpy.Delete_management(lyr3) 

        # one level clustering, no upper_ID, no cluster_type
        if upper_ID == "" or cluster_type == "":
            exp = "!" + out_class_item + "!"
            arcpy.CalculateField_management(lyr, cluster, exp , "PYTHON_9.3","#")           
            
        ### Dissolve all clusters at all levels according to the newly added integrated clusters.  
        isAMapLayer = False
        mxd = arcpy.mp.ArcGISProject("CURRENT")
        m = mxd.listMaps()[0]
        for mapLyr in m.listLayers():
            if mapLyr.name == inp_fc:
                isAMapLayer = True
                break

        #PrintAndWriteMessages("isAMapLayer = " + str(isAMapLayer) ,2)
        
        # The property of featureClass should not be used  if the input is NOT a layer in the map document 
        if isAMapLayer == False:
            if desc.dataType == 'ShapeFile':
                tmpMixedClusFC = "tmpMixedClusters.shp"
                MixedClusFC = mixed_clusters + ".shp"
            if desc.dataType == 'FeatureClass':
                tmpMixedClusFC = "tmpMixedClusters"
                MixedClusFC = mixed_clusters
            
        # The property of featureClass is required when describing a layer in ArcMap
        if isAMapLayer == True:
            if desc.datasetType == 'ShapeFile':
                tmpMixedClusFC = "tmpMixedClusters.shp"
                MixedClusFC = mixed_clusters + ".shp"
            if desc.datasetType == 'FeatureClass':
                tmpMixedClusFC = "tmpMixedClusters"
                MixedClusFC = mixed_clusters

        #PrintAndWriteMessages(tmpMixedClusFC, 2)

        if arcpy.Exists(tmpMixedClusFC):
            arcpy.Delete_management(tmpMixedClusFC)
            
        mixStatFlds = []
        for i in ConstraintList:
            mixStatFlds.append(i + " sum")
        if cluster_type != "":
            mixStatFlds.append(cluster_type + " min")
        if upper_ID != "":
            mixStatFlds.append(upper_ID + " first")
        mixStatFlds.append(isolate + " min")
        mixStatFlds = ";".join(mixStatFlds)
        arcpy.Dissolve_management(lyr, tmpMixedClusFC, cluster, mixStatFlds)
        
        flds = arcpy.ListFields(tmpMixedClusFC)
        fldClusType = ""; fldUpID = ""; fldIso = ""; fldsCap = []
        for fld in flds:
            #msg = "A field name in the temp mixed clusters is " + fld.name 
            #PrintAndWriteMessages(msg, 1)
            if fld.name.__contains__(cluster_type[:6]) or fld.name.__contains__(cluster_type[:6].upper()):
                fldClusType = fld.name
            if fld.name.__contains__("first_") or fld.name.__contains__("FIRST_"):
                fldUpID = fld.name
            if fld.name.__contains__(isolate[:6]) or fld.name.__contains__(isolate[:6].upper()):
                fldIso = fld.name
            for i in range(len(CapacityList)):
                if fld.name.__contains__(ConstraintList[i][:6]) or fld.name.__contains__(ConstraintList[i][:6].upper()):
                    fldsCap.append(fld.name)
        ##    msg = "In the mixed clusters, cluster type and upperID variables are: " + fldClusType + \
        ##         " and " + fldUpID
        ##    PrintAndWriteMessages(msg, 1)
            #msg = "In the mixed clusters, capacity variables are: " + ",".join(fldsCap)
            #PrintAndWriteMessages(msg, 1)
                        
        # mixed-level, replace all upper units to 000 for cluster type = 3
        if upper_ID != "" and cluster_type != "":
            arcpy.MakeFeatureLayer_management(tmpMixedClusFC, "tmplyr", fldClusType + " = 3")
            if int(arcpy.GetCount_management("tmplyr").getOutput(0)) > 0:
                arcpy.CalculateField_management("tmplyr", fldUpID, "str('" + specialID + "')", "PYTHON_9.3","#")
            arcpy.Delete_management("tmplyr")
                    
        # Near_analysis to calculate the nearest cluster of each dissolved cluster (mainly for the isolated ones)
        lyrIso = "lyrIso"; lyrNonIso = "lyrNonIso"
        where_clause = fldIso + " = 1"
        arcpy.MakeFeatureLayer_management(tmpMixedClusFC, lyrIso, where_clause)
        where_clause = fldIso + " = 0"
        arcpy.MakeFeatureLayer_management(tmpMixedClusFC, lyrNonIso, where_clause)
           
        oid_fld = arcpy.Describe(tmpMixedClusFC).OIDFieldName
        
        MixedClusLyr = "MixedClusLyr"
        arcpy.MakeFeatureLayer_management(tmpMixedClusFC, MixedClusLyr)
            
        lyr4 = "lyr4"; lyr5 = "lyr5"; lyr6 = "lyr6"; lyr7 = "lyr7"
        
        ### one level
        # No upper_ID, one-level clustering
        if upper_ID == "" or cluster_type == "":
              # New way: look for nearest non-isolated cluster to merge the isolated cluster
            # Remember, nearFID is based on the oid_fld of tmpMixedClusFC (and feature layers created from it)
            isoCount1 = int(arcpy.GetCount_management(lyrIso).getOutput(0))
            
            isoRows = arcpy.da.SearchCursor(lyrIso, cluster)   #all isolated dissolved clusters

            # Alternative to loop through all isolated clusters and call the merge function.
            # For some reason, the above old code does not work.
            isoClusters = []
            for isoRow in isoRows:
                isoClusters.append(isoRow[0])

            for i in range(isoCount1):
                #PrintAndWriteMessages("iteraiton # and cluster: " + str(i+1) + ", " + isoClusters[i], 2)
                # call the function to find and recalculate the isolated cluster
                MergeIsolated(tmpMixedClusFC, isoClusters[i], fldIso, "")
                #PrintAndWriteMessages(str(i+1) + " merge is done", 1)
                
            del isoRows
            # End of new way.
            
            if arcpy.Exists(lyr4):
                arcpy.Delete_management(lyr4)
            if arcpy.Exists(lyr5):
                arcpy.Delete_management(lyr5)
            if arcpy.Exists(lyr6):
                arcpy.Delete_management(lyr6)
            if arcpy.Exists(lyr7):
                arcpy.Delete_management(lyr7)                
        # End of one level: if upper_ID == "" or cluster_type == "":  

        ### mixed level
        if upper_ID != "" and cluster_type != "":
            # sub scenario 1: isolation at the lower level (tract), merge to a nearby same-level cluster, connected or not, e.g. St Martin
            # for each isolated cluster, select all sub clusters within the same upper unit.
            # There will be no adjacent cluster, simply find the nearest (old way is min) cluster.
            # Change cluster membership in original feature class (e.g. tract)
            where_clause = fldIso + " = 1 and " + fldClusType + " = 1"
            arcpy.MakeFeatureLayer_management(tmpMixedClusFC, lyr4, where_clause)
            isoCount21 = int(arcpy.GetCount_management(lyr4).getOutput(0))
            # New way: look for nearest cluster to merge the isolated cluster
            # Remember, nearFID is based on the oil_fld of tmpMixedClusFC (and feature layers created from it)

            where_clause = fldIso + " = 0 and " + fldClusType + " = 1"
            arcpy.MakeFeatureLayer_management(tmpMixedClusFC, lyr5, where_clause)
            arcpy.Near_analysis(lyr4, lyr5)

            searchFlds = [oid_fld, fldUpID, cluster, "NEAR_FID"]
            with arcpy.da.SearchCursor(lyr4, searchFlds) as cursor:  #all isolated dissolved clusters
                for row in cursor:
                    theFID = row[0]
                    theUpID = row[1]
                    theClus = row[2] 
                    nearFID = row[3]
                    
                    # Is the NEAR_FID within the same upper unit as theClus?
                    where_clause = oid_fld + " = " + str(nearFID)
                    rows2 = arcpy.da.SearchCursor(tmpMixedClusFC, [fldUpID, cluster], where_clause)
                    row2 = next(rows2)
                    nearUpID = row2[0]
                    # the original nearFID is in the same upper unit
                    if nearUpID == theUpID:
                        theNewClusID = row2[1]
                    else:
                        # the original nearFID is not in the same upper unit
                        # search for nearest within the same upper unit
                        where_clause = fldUpID + " = \'" + str(theUpID) + "\'"
                        arcpy.MakeFeatureLayer_management(tmpMixedClusFC, lyr6, where_clause)
                        
                        arcpy.MakeFeatureLayer_management(lyr6, "lyr6Iso", fldIso + " = 1")
                        arcpy.MakeFeatureLayer_management(lyr6, "lyr6NonIso", fldIso + " = 0")
                        arcpy.Near_analysis("lyr6Iso", "lyr6NonIso")

                        where_clause = oid_fld + " = " + str(theFID)
                        rows3 = arcpy.da.SearchCursor("lyr6Iso", "NEAR_FID", where_clause)
                        row3 = next(rows3)
                        nearFID = row3[0]
                        where_clause = oid_fld + " = " + str(nearFID)
                        rows4 = arcpy.da.SearchCursor(tmpMixedClusFC, cluster, where_clause)
                        row4 = next(rows4)
                        theNewClusID = row4[0]                    
                    del row2, rows2, row3, rows3, row4, rows4
                    
                    # Update original file (before clustering). It will be re-dissolved later.
                    # Keep the original membership, "out_class_item", and change the final membership, "cluster"
                    where_clause = cluster + " = \'" + str(theClus) + "\'"
                    arcpy.MakeFeatureLayer_management(inp_fc, lyr7, where_clause)
                    arcpy.CalculateField_management(lyr7, cluster, "\'"+ theNewClusID + "\'", "PYTHON_9.3","#")
            # End of new way.
                
            if arcpy.Exists(lyr4):
                arcpy.Delete_management(lyr4)
            if arcpy.Exists(lyr5):
                arcpy.Delete_management(lyr5)
            if arcpy.Exists(lyr6):
                arcpy.Delete_management(lyr6)
            if arcpy.Exists(lyr7):
                arcpy.Delete_management(lyr7)
                                
            # sub scenario 2: isolation at the upper level (county), merge to nearby clusters,
            # give the priority to upper level cluster.
            # Look for all adjacent clusters, give type = 2 or 3 first priority, since it keeps upper boundary. Find closest cluster
            # If there is no type2 or type3 cluster, look for closest type1 cluster and merge
            # If a mixed-level 1 & 3 is needed, define and count it
            count4 = 0
            where_clause = fldIso + " = 1 and " + fldClusType + " = 3"
            arcpy.MakeFeatureLayer_management(tmpMixedClusFC, lyr4, where_clause)
            isoCount22 = int(arcpy.GetCount_management(lyr4).getOutput(0))
        ##    msg = "The number of isolated type 3 clusters is " + str(isoCount2)  
        ##    PrintAndWriteMessages(msg, 1)
            
            # modified way, still works.  
            rows = arcpy.da.SearchCursor(lyr4,cluster)
            for row in rows:
                theClus = row[0]
                where_clause = cluster + " = \'" + str(theClus) + "\'"
                arcpy.MakeFeatureLayer_management(tmpMixedClusFC, lyr5, where_clause)
                arcpy.SelectLayerByLocation_management(MixedClusLyr, "SHARE_A_LINE_SEGMENT_WITH", lyr5)
                where_clause = cluster + " <> \'" + str(theClus) + "\'"
                arcpy.MakeFeatureLayer_management(MixedClusLyr, lyr6, where_clause)
        ##        msg = "The number of features in the selected layer is " + str(int(arcpy.GetCount(lyr6).GetOutput(0)))
        ##        PrintAndWriteMessages(msg, 1)

                where_clause = fldClusType + " = 2 or " + fldClusType + " = 3" 
                rows2 = arcpy.da.SearchCursor(lyr6,[fldsCap[0],cluster],where_clause)
                count = 0
                for row2 in rows2:
                    count += 1
                    if count == 1:
                        minCap = row2[0]
                        minClus = row2[1]
                    elif row2[0] < minCap:
                        minCap = row2[0]
                        minClus = row2[1]
                del rows2
                
                # If find a type2 cluster to merge, Go back to original file, change the previous cluster type and cluster membership
                if count >= 1:
                    where_clause = cluster + " = \'" + minClus + "\'"
                    arcpy.MakeFeatureLayer_management(inp_fc, lyr7, where_clause)
                    exp = "\'" + theClus + "\'"
                    arcpy.CalculateField_management(lyr7, cluster, exp, "PYTHON_9.3","#")       
                    arcpy.CalculateField_management(lyr7, cluster_type, "3", "PYTHON_9.3","#")
                    arcpy.Delete_management(lyr7)
                else:
                    # no adjacent type2 cluster, find a minCap type 1 cluster
                    rows3 = arcpy.da.SearchCursor(lyr6, [fldsCap[0],cluster],fldClusType + " = 1")
                    countType1 = 0
                    for row3 in rows3:
                        countType1 += 1
                        if countType1 == 1:
                            minCap = row3[0]
                            minClus = row3[1]
                        elif row3[0] < minCap:
                            minCap = row3[0]
                            minClus = row3[1]
                    del rows3

                    # If find a type1 cluster to merge, Go back to original file, change the previous clusters
                    # define and add count a new cluster type 4, means mixed level of 1 and 3
                    count4 += 1
                    where_clause = cluster + " = \'" + minClus + "\' or " + cluster + " = \'" + theClus + "\'"
                    arcpy.MakeFeatureLayer_management(inp_fc, lyr7, where_clause)
                    exp = "\'4." + specialID + "." + str(count4) + "\'"

                    #temp testing, will be removed
                    msg = "where_clause: " + where_clause
                    PrintAndWriteMessages(msg,1)
                    msg = "reassign cluster membership: " + cluster + " = " + exp
                    PrintAndWriteMessages(msg,1)
                    
                    arcpy.CalculateField_management(lyr7, cluster, exp, "PYTHON_9.3","#")        
                    arcpy.CalculateField_management(lyr7, cluster_type, "4", "PYTHON_9.3","#")
                    arcpy.Delete_management(lyr7)
                    
                #row = next(rows)
            del rows
            # end of old way  

            if arcpy.Exists(lyr4):
                arcpy.Delete_management(lyr4)
            if arcpy.Exists(lyr5):
                arcpy.Delete_management(lyr5)
            if arcpy.Exists(lyr6):
                arcpy.Delete_management(lyr6)
            if arcpy.Exists(MixedClusLyr):
                arcpy.Delete_management(MixedClusLyr)
            # End of mix level: if upper_ID != "" and cluster_type != "":    

        # re-dissolved mixed clusters after isolation-removal
        if isoCount1 >= 1 or isoCount21 >= 1 or isoCount22 >= 1:
            #oldClusters = desc.path + "\\" + tmpMixedClusFC
            #PrintAndWriteMessages(oldClusters, 1)
            arcpy.Delete_management(tmpMixedClusFC)
            PrintAndWriteMessages("Old mixed clusters have been deleted and new mixed clusters are created.", 1)
            if arcpy.Exists(MixedClusFC):
                arcpy.Delete_management(MixedClusFC)
            arcpy.Dissolve_management(lyr, MixedClusFC, cluster, mixStatFlds)
            
            # replace all upper units to 000 for cluster type = 3 or 4
            if upper_ID != "" and cluster_type != "":
                arcpy.MakeFeatureLayer_management(MixedClusFC, "tmplyr", fldClusType + " >= 3")
                if int(arcpy.GetCount_management("tmplyr").getOutput(0)) > 0:
                    arcpy.CalculateField_management("tmplyr", fldUpID, "str('" + specialID + "')", "PYTHON_9.3","#")                    
                arcpy.Delete_management("tmplyr")
        else:
            if arcpy.Exists(MixedClusFC):
                arcpy.Delete_management(MixedClusFC)

            arcpy.Rename_management(tmpMixedClusFC, MixedClusFC)
            PrintAndWriteMessages("No isolated clusters are found.", 1)
        
        # rename fields in dissolved clusters, to make them the same as the original file.
        flds = arcpy.ListFields(MixedClusFC)
            
        for i in range(len(CapacityList)):
            if ConstraintList[i] in flds:
                pass
            else:
                AddAField(MixedClusFC, ConstraintList[i], "Double")
                exp = "!" + fldsCap[i] + "!"
                arcpy.CalculateField_management(MixedClusFC, ConstraintList[i], exp, "PYTHON_9.3","#")
                arcpy.DeleteField_management(MixedClusFC, fldsCap[i])

        # mixed level only
        if upper_ID != "" and cluster_type != "":
            if cluster_type in flds:
                pass
            else:
                AddAField(MixedClusFC, cluster_type, "Short")
                exp = "!" + fldClusType + "!"
                arcpy.CalculateField_management(MixedClusFC, cluster_type, exp, "PYTHON_9.3","#")
                arcpy.DeleteField_management(MixedClusFC, fldClusType)

            if upper_ID in flds:
                pass
            else:
                AddAField(MixedClusFC, upper_ID, "Text")
                exp = "!" + fldUpID + "!"
                arcpy.CalculateField_management(MixedClusFC, upper_ID, str(exp), "PYTHON_9.3","#")
                arcpy.DeleteField_management(MixedClusFC, fldUpID)

        if isolate in flds:
            pass
        else:
            AddAField(MixedClusFC, isolate, "Short")
        exp = "!" + fldIso + "!"
        arcpy.CalculateField_management(MixedClusFC, isolate, exp, "PYTHON_9.3","#")
        arcpy.DeleteField_management(MixedClusFC, fldIso)
        
        arcpy.Delete_management(lyr)
        return MixedClusFC
        ### End of tackling isolated clusters.

    except:
        # Return any python specific errors as well as any errors from the geoprocessor
        tb = sys.exc_info()[2]
        tbinfo = traceback.format_tb(tb)[0]
        pymsg = "PYTHON ERRORS:\nTraceback Info:\n" + tbinfo + "\nError Info:\n    " + \
                str(sys.exc_type)+ ": " + str(sys.exc_value) + "\n"
        PrintAndWriteMessages(pymsg,2)

        msgs = "GP ERRORS:\n" + arcpy.GetMessages(2) + "\n"
        PrintAndWriteMessages(msgs,2)
#End of funciton TackleIsolation
        
### End of all three embedded funcitons

#-------------------------Prepare a summary report-------------------------------------#
def MLRreport(inp_fc):
    import pandas as pd
    from arcgis.features import GeoAccessor, GeoSeriesAccessor
    # Prepare a summary report
    # summarize each type's upper unit count in mixed clusters.
    sdf = pd.DataFrame.spatial.from_featureclass(inp_fc)
    CountOrig=sdf.shape[0]
    CountUpUnit=sdf[[upper_ID]].groupby(upper_ID).count().shape[0]
    CountMixed=sdf[[cluster]].groupby(cluster).count().shape[0]
    dictTypeFrqMix1=sdf[[cluster_type,upper_ID]].groupby(cluster_type).count().to_dict()
    dictTypeFrqMix=dictTypeFrqMix1[upper_ID]
    dictTypeFrq1=sdf[[upper_ID,cluster_type,cluster]].groupby(upper_ID).first().groupby(cluster_type).count().to_dict()
    dictTypeFrq=dictTypeFrq1[cluster]
    CountType1 = 0;CountType2 = 0;CountType3 = 0;CountType4 = 0
    CntUpperT1 = 0 ;CntUpperT2 = 0;CntUpperT3 = 0;CntUpperT4 = 0
    if 1 in dictTypeFrqMix:        
        CountType1 = dictTypeFrqMix[1]
        CntUpperT1 = dictTypeFrq[1]    
    if 2 in dictTypeFrqMix:
        CountType2 = dictTypeFrqMix[2]
        CntUpperT2 = dictTypeFrq[2]
    if 3 in dictTypeFrqMix:
        CountType3 = dictTypeFrqMix[3]
        CntUpperT3 = dictTypeFrq[3]
    if 4 in dictTypeFrqMix:
        CountType4 = dictTypeFrqMix[4]
        CntUpperT4 = dictTypeFrq[4]

    msg = "\n" + "Result summary: " + "\n" + \
        "The number of original lower-level units is " + str(CountOrig) + ".\n" \
        "The number of original upper-level units is " + str(CountUpUnit) + ".\n" \
        "The number of final mixed clusters is " + str(CountMixed) + ", including \n" + \
        " - " + str(CountType1) + " type1 clusters (single or multi lower-level units) " + \
                "from " + str(CntUpperT1) + " upper units.\n" + \
        " - " + str(CountType2) + " type2 clusters (single upper-level unit) " + \
                "from " + str(CntUpperT2) + " upper units.\n" + \
        " - " + str(CountType3) + " type3 clusters (multi upper-level units) " + \
                "from " + str(CntUpperT3) + " upper units.\n" + \
        " - " + str(CountType4) + " type4 clusters (mixed upper and lower level units) " + \
                "from " + str(CntUpperT4) + " upper units.\n"
    PrintAndWriteMessages(msg,1)

    keepPct = (float(CntUpperT1) + float(CntUpperT2)) / float(CountUpUnit) * 100
    msg = "\n" + "Mixed-level clustering is done! \n" + \
        "The boundaries of " + str(keepPct) + "% upper-level units are preserved."
    PrintAndWriteMessages(msg,1) 
  
#Import the modules
import arcpy, sys,os, traceback,math


#Use the lowest available license
for product in ['Engine','ArcView', 'ArcEditor', 'EngineGeoDB','ArcInfo', 'ArcServer']:
    if arcpy.CheckProduct(product).lower() == 'available':
        arcpy.SetProduct(product)
        break


arcpy.env.overwriteOutput  = True

try:
    #Create a log file with messages in the same location as the script

    #-------------------------Get the inputs

    inp_fc0 = arcpy.GetParameterAsText(0)
    PeanoOrder_fld = arcpy.GetParameterAsText(1)   # Spatial order by Peano Curve (0-1)
    AttriOrder_fld = arcpy.GetParameterAsText(2)   # Aggregated weighted attributive order
    NAttriOrder_fld = arcpy.GetParameterAsText(3)  # Normalized attribute order (0-1)
    Order_fld = arcpy.GetParameterAsText(4)        # Integrated order field to be created for each one-level clustering
    Attri_list = arcpy.GetParameterAsText(5)
    Wt_list = arcpy.GetParameterAsText(6)
    Pct_Peano = float(arcpy.GetParameterAsText(7))

    constraint_list = arcpy.GetParameterAsText(8)  # List of vairables with lower limits
    capacity_list = arcpy.GetParameterAsText(9)    # Values of the lower limit for each constraint variable.
    out_class_item = arcpy.GetParameterAsText(10)   # Sub cluster membership in each one-level clustering
    isolate = arcpy.GetParameterAsText(11)          # Whether a cluster is nusatisfied and isoloated (0 or 1)
    
    upper_ID = arcpy.GetParameterAsText(12)         # string or number, field to represent the upper level unit, e.g. county
    cluster_type = arcpy.GetParameterAsText(13)     # short integer, new variable added to table for cluster types, "ClusType" 

    DissolveWt = arcpy.GetParameterAsText(14)      # Weight (a field) of calculating dissolved units
    cluster = arcpy.GetParameterAsText(15)         # final cluster ID field
    mixed_clusters = arcpy.GetParameterAsText(16)   # final feature class of mixed clusters

    Pct_Attri = 100 - Pct_Peano

    #parse the input parameter strings to lists
    ConstraintList = constraint_list.split(';')
    CapacityList = capacity_list.split(';')
    AttriList = Attri_list.split(';')
    WtList = Wt_list.split(';')

    # Get the path and full name of the input feature class
    desc = arcpy.Describe(inp_fc0)
    arcpy.env.workspace = desc.path
    fullInp_fc = desc.BaseName + "." + desc.Extension

    ## Major changes for using da.UpdateCursor and da.SearchCursor
    ## da.UpdateCursor sql_clause 'ORDER BY' option only works with (geo)database, does not work with dbf or info table. 
    ##
    if desc.datasetType == 'ShapeFile':
        resultGDB = "resultGDB"
        if arcpy.Exists(arcpy.env.workspace + "\\" + resultGDB + ".gdb"):
            arcpy.Delete_management(resultGDB + ".gdb")
        
        arcpy.CreateFileGDB_management(arcpy.env.workspace, resultGDB)
        GDBLoc = arcpy.env.workspace + "\\" + "resultGDB.gdb"
        arcpy.FeatureClassToGeodatabase_conversion(inp_fc0, GDBLoc)
        inp_fc = GDBLoc + "\\" + desc.BaseName
        arcpy.env.workspace = GDBLoc
    if desc.datasetType == 'FeatureClass':
        inp_fc = inp_fc0
        
    #Validate the output field name for the workspace
    PeanoOrder_fld = arcpy.ValidateFieldName(PeanoOrder_fld,os.path.dirname(inp_fc))
    AttriOrder_fld = arcpy.ValidateFieldName(AttriOrder_fld,os.path.dirname(inp_fc))
    NAttriOrder_fld = arcpy.ValidateFieldName(NAttriOrder_fld,os.path.dirname(inp_fc))
    Order_fld = arcpy.ValidateFieldName(Order_fld,os.path.dirname(inp_fc))
    out_class_item = arcpy.ValidateFieldName(out_class_item,os.path.dirname(inp_fc))
    isolate = arcpy.ValidateFieldName(isolate,os.path.dirname(inp_fc))
    cluster_type = arcpy.ValidateFieldName(cluster_type,os.path.dirname(inp_fc))

    # Execute the Simmarize Statistics tool using the SUM option and case field (e.g. county)
    msg = "\n" + "All-in-one mixed-clustering script starts to run." + \
          "\n" + "Be patient, this could take a while." + \
          "\n\n" + "Aggegate the units to the upper level to determine cluster type at the upper level." + \
          "The mixed-cluster results will have cluster type of 1,2, 3 or 4, where" + \
          "\n" + " # 1 = lower level clusters (e.g.tract)" + \
          "\n" + " # 2 = single upper level clusters (e.g. a county)" + \
          "\n" + " # 3 = multi upper level clusters (e.g. multi-county)" + \
          "\n" + " # 4 = mixed lower and upper level clusters (e.g., tract and county)"
          
    PrintAndWriteMessages(msg,1)  
    
    tmpTbUpper = "sum_upper"
    if arcpy.Exists(tmpTbUpper):
        arcpy.Delete_management(tmpTbUpper)

    sumfldsDscp = []
    for constraint in ConstraintList:
        sumfldsDscp.append(constraint + " sum")
    sumfldsDscp = ";".join(sumfldsDscp)
    
    arcpy.Statistics_analysis(inp_fc, tmpTbUpper, sumfldsDscp, upper_ID)
    
    #Add the new field to the summary table. If the field already exists, AddField has a warning
    result = AddAField(tmpTbUpper, PeanoOrder_fld, "DOUBLE")
    result = AddAField(tmpTbUpper, AttriOrder_fld, "DOUBLE")
    result = AddAField(tmpTbUpper, NAttriOrder_fld, "DOUBLE")
    result = AddAField(tmpTbUpper, Order_fld, "DOUBLE")
    result = AddAField(tmpTbUpper, out_class_item, "Short")
    result = AddAField(tmpTbUpper, isolate, "Short")
    result = AddAField(tmpTbUpper, cluster_type, "Short")
           
    # Get the sum fields to be used later
    flds = arcpy.ListFields(tmpTbUpper)
    fldSum = []
    for fld in flds:
        if fld.name.__contains__("SUM_") or fld.name.__contains__("sum_"):
            fldSum.append(fld.name)
            
    # Go through every record to assign cluster type to each summarized upper unit (e.g. county).
    curFlds = fldSum
    curFlds.append(cluster_type)
    #PrintAndWriteMessages(curFlds,2)

    rows = arcpy.da.UpdateCursor(tmpTbUpper, curFlds)
    typeIndex = curFlds.index(cluster_type)
    #rows.reset()
    #row = next(rows)
    for row in rows:
        measures = []
        for fld in fldSum:
            if fld.__contains__("SUM_") or fld.__contains__("sum_"):
                measures.append(row[curFlds.index(fld)])

        # Cluster type = 1, to be disaggregated within one upper unit 
        met2All = 1; metAll = 1
        for i in range(len(measures)):            
            if float(measures[i]) >= 2 * float(CapacityList[i]):
                met2All = met2All * 1
            else:
                met2All = met2All * 0

            if float(measures[i]) >= float(CapacityList[i]):
                metAll = metAll * 1
            else:
                metAll = metAll * 0

        if met2All == 1:
            row[typeIndex] = 1
            
        # Cluster type = 2, good at upper unit level
        elif metAll == 1:
            row[typeIndex] = 2

        # Cluster type = 3, to be agrregated at upper unit level 
        else:
            row[typeIndex] = 3

        rows.updateRow(row)
        #del row
        #row = next(rows)
    del row, rows

    # If cluster_type field exists in the input feature class, delete it; otherwise the following joinfield will be conflicted.   
    #Call a function
    DelAField(inp_fc, cluster_type)

    #arcpy.env.overwriteOutput = True        
    #joinField to join the sum table back to the original data table, and permanently add the cluster_type field
    #Usage: JoinField <in_data> <in_field> <join_table> <join_field> {fields;fields...}
    arcpy.JoinField_management(inp_fc, upper_ID, tmpTbUpper, upper_ID, cluster_type)

    #Get the OID field name
    #oid_fld = desc.OIDFieldName
    oid_fld = arcpy.Describe(inp_fc).OIDFieldName
    
    # use a field to copy FID/OID to the new integer field, because FID/OID cannot be used in GenerateSpatialWeightsMatrix
    result = AddAField(inp_fc,"TheID","Short")
    
      
    # A pair of "!" signs for using existing variables (attribute table) in Python expression.  
    exp = "!" + oid_fld + "!"
    arcpy.CalculateField_management(inp_fc, "TheID", exp, "PYTHON_9.3","#")

    # creat a list of feature class layers for inputing into clustering tool (script)
    lyr = "tmplyr"
    #ready to call the functions in other scripts

    #------------------------ Scenario  1: upper unit need to be disaggregated----------------#
    msg = "\n" + "Now start processing scenario 1: an upper unit needs to be disaggregated." + \
          "\n" + "Clustering order accounting for both spatial (Peano Curve) and attributive measures " + \
          "will be calculated and then lower units will be clustered within each upper unit."       
    PrintAndWriteMessages(msg,0)
    count1 = 0
    where_clause = cluster_type + " = 1"
    rows = arcpy.da.SearchCursor(tmpTbUpper, upper_ID, where_clause)

    for row in rows:
    #while count1 <= 3:  # use a small sample for testing only
        count1 += 1
        where_clause2 = upper_ID + " = '" + row[0] + "'"
        
        arcpy.MakeFeatureLayer_management(inp_fc, lyr, where_clause2)        
        ClusteringOrder \
            (lyr, PeanoOrder_fld, AttriOrder_fld, NAttriOrder_fld, Order_fld, Attri_list, Wt_list, Pct_Peano)
        OneLevelClustering \
            (lyr, Order_fld, constraint_list, capacity_list, out_class_item, isolate)
        
        msg = where_clause2 + " is done. Its cluster type is 1."
        PrintAndWriteMessages(msg,1)        

        arcpy.Delete_management(lyr)

    del rows

    #------------------------ Scenario  2: upper unit is good enough to be a cluster
    msg = "\n" + "Now start processing scenario 2: an upper unit is good enough to be a cluster."
    PrintAndWriteMessages(msg,0)
    count2 = 0      
    where_clause = cluster_type + " = 2"
    rows = arcpy.da.SearchCursor(tmpTbUpper, upper_ID, where_clause)
    for row in rows:
        count2 += 1    
        where_clause2 = upper_ID + " = '" + row[0] + "'"
        arcpy.MakeFeatureLayer_management(inp_fc, lyr, where_clause2)

        result = AddAField(lyr,"included","Short")
        arcpy.SetProgressorPosition()
        arcpy.CalculateField_management(lyr,"included",1,"PYTHON_9.3","#")
        arcpy.CalculateField_management(lyr,out_class_item, count2, "PYTHON_9.3","#")
        
        msg = where_clause2 + " is done. Its cluster type is 2."
        PrintAndWriteMessages(msg,1)
        
        #del row
        arcpy.Delete_management(lyr)
        #row = next(rows)
    del rows

    #------------------------ Scenario  3: upper units need to be aggregated to larger unit to satisfy the clustering criteria.
    msg = "\n" + "Now start processing scenario 3: upper units need to be aggregated to " + \
          "larger units to satisfy the clustering criteria."
    PrintAndWriteMessages(msg,0)
    count3 = 0   
    where_clause = cluster_type + " = 3"
    arcpy.MakeFeatureLayer_management(inp_fc, lyr, where_clause)
    count3 = int(arcpy.GetCount_management(lyr).getOutput(0))
    if count3 > 0:
        dissolvedUpperFC = "dissolvedUpperFC"
        if arcpy.Exists(dissolvedUpperFC):
            arcpy.Delete_management(dissolvedUpperFC)

        # Use sumfldsDscp defined previously in this script
        # DissolveWt could be one of the contraints or NOT.
        if DissolveWt in ConstraintList:
            statFields = sumfldsDscp + "; TheID MIN"
            DissolveWtIn = 1
        else:
            statFields = sumfldsDscp + "; TheID MIN; " + DissolveWt + " SUM"
            DissolveWtIn = 0

        msg = "The field (weight) used to aggregate attributes to upper level is " + DissolveWt
        PrintAndWriteMessages(msg,1)
        
        arcpy.Dissolve_management(lyr, dissolvedUpperFC, upper_ID, statFields)
        AddAField(dissolvedUpperFC, cluster_type, "Short")
        arcpy.CalculateField_management(dissolvedUpperFC,cluster_type, "3", "PYTHON_9.3","#")

        # Get the fields of SUM or MIN to be used later
        flds = arcpy.ListFields(dissolvedUpperFC)
        fldSumDissolve = []; fldMin = ""; fldDissolveWt = ""
        for fld in flds:
            if fld.name.__contains__("SUM_") or fld.name.__contains__("sum_"):
                if DissolveWtIn == 1 or (DissolveWtIn == 0 and fld.name.__contains__(DissolveWt[:5]) == False):
                    fldSumDissolve.append(fld.name)
            if fld.name.__contains__("MIN_") or fld.name.__contains__("min_"):
                fldMin = fld.name
            if fld.name.__contains__(DissolveWt[:5]):
                fldDissolveWt = fld.name         
       
        # Calculate weighted dissolved attributes, e.g. Factors 1, 2, & 3 weighted by population
        #Usage: JoinField <in_data> <in_field> <join_table> <join_field> {fields;fields...}
        DelAField(lyr, fldDissolveWt)
        arcpy.JoinField_management(lyr, upper_ID, dissolvedUpperFC, upper_ID, fldDissolveWt)

        WtAttriList = []
        for i in range(len(AttriList)):
            WtAttriList.append("W" + AttriList[i])
            AddAField(lyr, WtAttriList[i], "DOUBLE")
        
        curFlds = []
        curFlds.append(DissolveWt)
        curFlds.append(fldDissolveWt)
        curFlds.extend(WtAttriList)
        curFlds.extend(AttriList)
        #PrintAndWriteMessages(curFlds,2)
        
        rows = arcpy.da.UpdateCursor(lyr, curFlds, fldDissolveWt + " > 0")
        #rows.reset()
        #row = next(rows)
        for row in rows:
            #theRatio = row.getValue(DissolveWt) / row.getValue(fldDissolveWt)
            #for i in range(len(WtAttriList)):
            #     row.setValue(WtAttriList[i], row.getValue(AttriList[i]) * theRatio)

            theRatio = row[curFlds.index(DissolveWt)] / row[curFlds.index(fldDissolveWt)]
            for i in range(len(WtAttriList)):
                row[curFlds.index(WtAttriList[i])]= row[curFlds.index(AttriList[i])]
                 
            rows.updateRow(row)
            #del row
            #row = next(rows)
        del rows

        #summarize weighted attributes by upperID, then integrate them to the dissolved table
        tmpTbSumWtAttri = "sum_attri"
        if arcpy.Exists(tmpTbSumWtAttri):
            arcpy.Delete_management(tmpTbSumWtAttri)

        sumAttrifldsDscp = []
        for Wtattri in WtAttriList:
            sumAttrifldsDscp.append(Wtattri + " sum")
        sumAttrifldsDscp = ";".join(sumAttrifldsDscp)

        arcpy.Statistics_analysis(lyr, tmpTbSumWtAttri, sumAttrifldsDscp, upper_ID)
        flds = arcpy.ListFields(tmpTbSumWtAttri)
        fldDissolveAttri = []
        for fld in flds:
            if fld.name.__contains__("SUM_") or fld.name.__contains__("sum_"):
                fldDissolveAttri.append(fld.name)
        DissolveAttri_list = ";".join(fldDissolveAttri)        

        arcpy.JoinField_management(dissolvedUpperFC, upper_ID, tmpTbSumWtAttri, upper_ID, DissolveAttri_list)
                          
        # add a field TheID to copy MIN_TheID to the new integer field.
        result = AddAField(dissolvedUpperFC,"TheID","Short")
    ##    if result.MaxSeverity == 1:
    ##        msg = result.GetMessages(1).split(':')[1].lstrip() + " For the dissolved feature class, existing values of TheID will be overwritten."
    ##        PrintAndWriteMessages(msg,1)    
        exp = "!" + str(fldMin) + "!"
        arcpy.CalculateField_management(dissolvedUpperFC,"TheID", exp, "PYTHON_9.3","#")

        dissolvedUpperFlyr = "dissolvedUpperFlyr"
        arcpy.MakeFeatureLayer_management(dissolvedUpperFC, dissolvedUpperFlyr)
        DissolveConstraint_list = ";".join(fldSumDissolve)       

        ClusteringOrder \
            (dissolvedUpperFlyr, PeanoOrder_fld, AttriOrder_fld, \
             NAttriOrder_fld, Order_fld, DissolveAttri_list, Wt_list, Pct_Peano)

        # In case there are too few features
        try :
            if int(arcpy.GetCount_management(dissolvedUpperFlyr).getOutput(0)) > 1:
                OneLevelClustering \
                    (dissolvedUpperFlyr, Order_fld, DissolveConstraint_list, \
                    capacity_list, out_class_item, isolate)    
        except:
            msg = "There are isolate features ."
            PrintAndWriteMessages(msg,1)
            print (msg)
            # arcpy.CalculateField_management(dissolvedUpperFlyr, Order_fld,  "!" +Order_fld + "!" , "PYTHON_9.3","#")
            arcpy.CalculateField_management(dissolvedUpperFlyr, out_class_item, "!MIN_TheID!" , "PYTHON_9.3","#")
            arcpy.CalculateField_management(dissolvedUpperFlyr, "included", 0 , "PYTHON_9.3","#")
            arcpy.CalculateField_management(dissolvedUpperFlyr, isolate, 1 , "PYTHON_9.3","#")



        # AddJoin to join the dissolved table back to the original data table, and calculate/change field values
        desc2 = arcpy.Describe(dissolvedUpperFlyr)
        joinName = desc2.BaseName
        inpName = desc.BaseName

        # Call a function to ddd a field if not added previously.
        AddAField(lyr, Order_fld, "DOUBLE")
        AddAField(lyr, out_class_item, "Short")
        AddAField(lyr, "included", "Short")
        AddAField(lyr, isolate, "Short")
        
        arcpy.AddJoin_management(lyr, upper_ID, dissolvedUpperFlyr, upper_ID)     # lyr is previously defined as cluster_type = 3
        arcpy.CalculateField_management(lyr, inpName + "." + Order_fld, "!" + joinName + "." + Order_fld + "!", "PYTHON_9.3","#")
        arcpy.CalculateField_management(lyr, inpName + "." + out_class_item, "!" + joinName + "."  + out_class_item + "!", "PYTHON_9.3","#")
        arcpy.CalculateField_management(lyr, inpName + "." + "included", "!" + joinName + "." + "included" + "!", "PYTHON_9.3","#")
        arcpy.CalculateField_management(lyr, inpName + "." + isolate, "!" + joinName + "."  + isolate + "!", "PYTHON_9.3","#")
        arcpy.RemoveJoin_management(lyr, joinName)

        arcpy.Delete_management(lyr)  
        arcpy.Delete_management(dissolvedUpperFlyr)
        arcpy.Delete_management(dissolvedUpperFC)
        msg = where_clause + " is done. Its cluster type is 3."
        PrintAndWriteMessages(msg,1)

    else:
        msg = "There is no cluster type 3 in this data."
        PrintAndWriteMessages(msg,1)
        
    ### Additional codes to tackle isolated cluster at different levels.   
    # Add a new field of integrated cluster membership and dissolve all clusters at all levels.
    # format: T.UUU.S, where T = cluster type (1,2,3), UUU = upper level code (e.g. county FIPS), and S = sub cluster ID 
    msg = "\n" + "Now start tackling isolated clusters at all levels." + \
          "\n" + "Some clusters' type and membership might be changed and a new cluster type of '4' " + \
          "might be created."
    PrintAndWriteMessages(msg,0)
    
    #-------------------------Func_TackleIsolation
    MixedClusFC =  TackleIsolation(inp_fc, out_class_item, isolate, upper_ID, cluster_type, \
                    cluster, mixed_clusters, constraint_list, capacity_list)

    #-------------------------new code 2022 update cluster_type --------------#
    #Keep the orginal Cluster Type(1,2,3,4) and build a Update Type for Input Feature
    exp="!"+cluster+"!"+"[0]"
    result = AddAField(inp_fc,"NewClusType","Short")  
    arcpy.CalculateField_management(inp_fc, "NewClusType", exp, "PYTHON_9.3","#")
    # Update Type for the result map  mixed_clusters
    arcpy.CalculateField_management(mixed_clusters, cluster_type, exp, "PYTHON_9.3","#")

    
    #-------------------------new code 2022 using Spatial Embeded DataFrame--------------#
    MLRreport(inp_fc)  
    #-------------------------new code 2022 Add results to the map document--------------#
    mxd = arcpy.mp.ArcGISProject("CURRENT")
    m = mxd.listMaps()[0]
    desc = arcpy.Describe(inp_fc)
    datapath = desc.path+"\\"+mixed_clusters
    m.addDataFromPath(datapath)               

    PrintAndWriteMessages("End of mixed-level clustering.", 1)   

    
except:
    # Return any python specific errors as well as any errors from the geoprocessor
    tb = sys.exc_info()[2]
    tbinfo = traceback.format_tb(tb)[0]
    pymsg = "PYTHON ERRORS:\nTraceback Info:\n" + tbinfo + "\nError Info:\n    " + \
            str(sys.exc_type)+ ": " + str(sys.exc_value) + "\n"
    PrintAndWriteMessages(pymsg,2)

    msgs = "arcpy ERRORS:\n" + arcpy.GetMessages(2) + "\n"
    PrintAndWriteMessages(msgs,2)

