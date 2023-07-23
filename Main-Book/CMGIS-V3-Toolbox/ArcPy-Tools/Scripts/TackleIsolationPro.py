#
# Embedded code, with the function embedded at the top of the code.
#
# Tackle isolated features/clusters at all levels, and merge them with other clusters.
# A function  to tackle isolated cluster at different levels.   
# Add a new field of integrated cluster membership and dissolve all clusters at all levels.
# format: T.UUU.S, where T = cluster type (1,2,3), UUU = upper level code (e.g. county FIPS), and S = sub cluster ID 
# clusType = 1,2, 3 or 4
# Where 1 = low level clusters (tract)
#       2 = mid level clusters (county)
#       3 = high level clusters (multi-county)
#       4 = mixed low and high level clusters (tract and county)
######################################################################
def TackleIsolation(inp_fc, out_class_item, isolate, upper_ID, cluster_type, cluster, mixed_clusters, constraint_list, capacity_list):
    #prints a arcpy message and writes the same message to a log file
    def PrintAndWriteMessages(msg,severity=0):
        if severity == 0:
            arcpy.AddMessage(msg)
        elif severity == 1:
            arcpy.AddWarning(msg)
        elif severity == 2:
            arcpy.AddError(msg)
        #logfile.write(msg + "\n")
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
            row3 = rows3.next()
            nearFID = row3[0]
            where_clause = oid_fld + " = " + str(nearFID)
            #PrintAndWriteMessages(where_clause, 1)
            rows4 = arcpy.da.SearchCursor(lyr7, cluster, where_clause)
            row4 = rows4.next()
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
        #mxd = arcpy.mapping.MapDocument("CURRENT")
        mxd = arcpy.mp.ArcGISProject("CURRENT")
        m = mxd.listMaps()[0]
        #df = mxd.activeDataFrame
        #lyrs = arcpy.mapping.ListLayers(mxd, "*", df)
        #lyrs = arcpy.mp.ListLayers(mxd, "*", df)
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
            if desc.featureClass.dataType == 'ShapeFile':
                tmpMixedClusFC = "tmpMixedClusters.shp"
                MixedClusFC = mixed_clusters + ".shp"
            if desc.featureClass.dataType == 'FeatureClass':
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
            #where_clause = fldIso + " = 1"
            #arcpy.MakeFeatureLayer_management(tmpMixedClusFC, lyr4, where_clause)
            #isoCount1 = long(arcpy.GetCount_management(lyr4).getOutput(0))
            
##            #old way: look for smallest-capacity cluster to merge the isolated cluster
##            rows = arcpy.da.SearchCursor(lyr4, cluster)
##            for row in rows: 
##                theClus = row[0]
##                where_clause = cluster + " <> \'" + theClus + "\'"
##                arcpy.MakeFeatureLayer_management(tmpMixedClusFC, lyr5, where_clause )
##
##                rows2 = arcpy.da.SearchCursor(lyr5, [fldsCap[0],cluster])
##                count = 0
##                for row2 in rows2:
##                    count += 1
##                    if count == 1:
##                        minCap = row2[0]
##                        minClus = row2[1]
##                    elif row2[0] < minCap:
##                        minCap = row2[0]
##                        minClus = row2[1]              
##                del rows2
##
##                # Update original file (before clustering).  It will be re-dissolved later.
##                where_clause = cluster + " = \'" + theClus + "\'"
##                arcpy.MakeFeatureLayer_management(inp_fc, lyr6, where_clause)
##                exp = "\'" + minClus + "\'"
##                arcpy.CalculateField_management(lyr6, cluster, exp, "PYTHON_9.3","#")
##                arcpy.CalculateField_management(lyr6, isolate, "0", "PYTHON_9.3","#")
##                theSubClusID = minClus.split(".")[-1]
##                arcpy.CalculateField_management(lyr6, out_class_item, theSubClusID, "PYTHON_9.3","#")
##            del rows
            # New way: look for nearest non-isolated cluster to merge the isolated cluster
            # Remember, nearFID is based on the oid_fld of tmpMixedClusFC (and feature layers created from it)
            isoCount1 = int(arcpy.GetCount_management(lyrIso).getOutput(0))
            
            isoRows = arcpy.da.SearchCursor(lyrIso, cluster)   #all isolated dissolved clusters
##            thisCount = 0
##            for isoRow in isoRows:
##                thisCount += 1
##                theClus = isoRow[0]
##                PrintAndWriteMessages("iteraiton # and cluster: " + str(thisCount) + ", " + theClus, 2)
##                # call the function to find and recalculate the isolated cluster
##                MergeIsolated(tmpMixedClusFC, theClus, fldIso, "")
##                PrintAndWriteMessages("One MergeIsolated() is done", 1)
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
##            #old way: look for smallest-capacity cluster to merge the isolated cluster
##            rows = arcpy.da.SearchCursor(lyr4, [fldUpID, cluster])
##            for row in rows: 
##                theUpID = row[0]
##                theClus = row[1]
##                where_clause = fldUpID + " = \'" + theUpID + "\' AND " + cluster + " <> \'" + theClus + "\'"
##                arcpy.MakeFeatureLayer_management(tmpMixedClusFC, lyr5, where_clause )
##
##                rows2 = arcpy.da.SearchCursor(lyr5, [fldsCap[0],cluster])
##                count = 0
##                for row2 in rows2:
##                    count += 1
##                    if count == 1:
##                        minCap = row2[0]
##                        minClus = row2[1]
##                    elif row2[0] < minCap:
##                        minCap = row2[0]
##                        minClus = row2[1]              
##                del rows2
##
##                # Go back to original files
##                where_clause = cluster + " = \'" + theClus + "\'"
##                arcpy.MakeFeatureLayer_management(inp_fc, lyr6, where_clause)
##                exp = "\'" + minClus + "\'"
##                arcpy.CalculateField_management(lyr6, cluster, exp, "PYTHON_9.3","#")
##                arcpy.CalculateField_management(lyr6, isolate, "0", "PYTHON_9.3","#")
##                theSubClusID = minClus.split(".")[-1]
##                arcpy.CalculateField_management(lyr6, out_class_item, theSubClusID, "PYTHON_9.3","#")
##            del rows
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
                    #row2 = rows2.next()
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
                        #row3 = rows3.next()
                        row3 = next(rows3)
                        nearFID = row3[0]
                        where_clause = oid_fld + " = " + str(nearFID)
                        rows4 = arcpy.da.SearchCursor(tmpMixedClusFC, cluster, where_clause)
                        #row4 = rows4.next()
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
        ##        msg = "The number of features in the selected layer is " + str(long(arcpy.GetCount(lyr6).GetOutput(0)))
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
                    
                #row = rows.next()
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
            # End of mix level: if upper_ID <> "" and cluster_type <> "":    
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
### Original TackelIsolation.py code with "import Func_TackleIsolation" removed.
import arcpy, sys, os, traceback
arcpy.env.overwriteOutput  = True
#Get the inputs
inp_fc = arcpy.GetParameterAsText(0)
out_class_item = arcpy.GetParameterAsText(1)   # Sub cluster membership in each one-level clustering
isolate = arcpy.GetParameterAsText(2)          # Whether a cluster is nusatisfied and isoloated (0 or 1)
upper_ID = arcpy.GetParameterAsText(3)         # string or number, field to represent the upper level unit, e.g. county
cluster_type = arcpy.GetParameterAsText(4)     # short integer, new variable added to table for cluster types, "ClusType" 
cluster = arcpy.GetParameterAsText(5)          # final cluster ID field in the original file
mixed_clusters = arcpy.GetParameterAsText(6)    # final feature class of mixed clusters
constraint_list = arcpy.GetParameterAsText(7)  # List of vairables with lower limits
capacity_list = arcpy.GetParameterAsText(8)    # Values of the lower limit for each constraint variable.

#ready to call the functions in other scripts
##import Func_TackleIsolation
##MixedClusFC = Func_TackleIsolation.TackleIsolation(inp_fc, out_class_item, isolate, upper_ID, cluster_type, \
##                                     cluster, mixed_clusters, constraint_list, capacity_list)
MixedClusFC = TackleIsolation(inp_fc, out_class_item, isolate, upper_ID, cluster_type, \
                                     cluster, mixed_clusters, constraint_list, capacity_list)
# Add results to the map document
# arcpy.CopyFeatures_management(MixedClusFC, outp_fc)
#mxd = arcpy.mapping.MapDocument("CURRENT")
mxd = arcpy.mp.ArcGISProject("CURRENT")
#df = mxd.activeDataFrame
#resultLyr = arcpy.mapping.Layer(MixedClusFC)
m = mxd.listMaps()[0]
#resultLyr = "MixClusters"
#arcpy.MakeFeatureLayer_management(MixedClusFC,resultLyr)
desc = arcpy.Describe(inp_fc)
datapath = desc.path+"\\"+mixed_clusters
# arcpy.CopyFeatures_management(MixedClusFC,datapath)

#resultLyr = arcpy.mp.LayerFile(MixedClusFC)
#resultLyr_file = "MixClusters.lyr"
#arcpy.SaveToLayerFile_management("study_quads_lyr", resultLyr_file, "ABSOLUTE")
#arcpy.mapping.AddLayer(df, resultLyr, "TOP")
#m.addLayer(resultLyr_file)
m.addDataFromPath(datapath)
#The End
