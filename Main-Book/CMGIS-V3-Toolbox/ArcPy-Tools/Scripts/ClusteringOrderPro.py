#*****************************************************************************************
# Embedded version, the function (used to be a standalone file) is added to the top of the code 
#
#This script is a stand alone version to calculate clustering order (not the function version) 
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

def ClusteringOrder(inp_fc, PeanoOrder_fld, AttriOrder_fld, NAttriOrder_fld, Order_fld, Attri_list, Wt_list, Pct_Peano):
    
    #prints a GP message and writes the same message to a log file
    def PrintAndWriteMessages(msg,severity=0):
        if severity == 0:
            arcpy.AddMessage(msg)
        elif severity == 1:
            arcpy.AddWarning(msg)
        elif severity == 2:
            arcpy.AddError(msg)
        #logfile.write(msg + "\n")

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
        #Create a log file with messages in the same location as the script
        #logfile = open(os.path.splitext(sys.argv[0])[0] + "_log.txt","w",0)

##This part is in the standalone python script but should be removed in a function.
##        #Get the inputs
##        inp_fc = arcpy.GetParameterAsText(0)
##        PeanoOrder_fld = arcpy.GetParameterAsText(1)
##        AttriOrder_fld = arcpy.GetParameterAsText(2)
##        NAttriOrder_fld = arcpy.GetParameterAsText(3)
##        Order_fld = arcpy.GetParameterAsText(4)
##        Attri_list = arcpy.GetParameterAsText(5)
##        Wt_list = arcpy.GetParameterAsText(6)
##        Pct_Peano = float(arcpy.GetParameterAsText(7))
        Pct_Attri = 100 - Pct_Peano
                    
    ##    #Get a feature count to set up the progressor
    ##    tot_features = long(arcpy.GetCount_management(inp_fc).GetOutput(0))
    ##    
    ##    #Create a progressor
    ##    arcpy.SetProgressor("step","Computing Spatial Order",0,tot_features + 1,1)
        
        #Add the new field if it does not exist, else overwrite existing values
##        msg = "Adding spatial order fields " + PeanoOrder_fld + ", " + AttriOrder_fld + ", " + NAttriOrder_fld+ ", " + Order_fld
##        arcpy.SetProgressorLabel(msg)
##        PrintAndWriteMessages(msg,0)

        #First Validate the field name
        PeanoOrder_fld = arcpy.ValidateFieldName(PeanoOrder_fld,os.path.dirname(inp_fc))
        AttriOrder_fld = arcpy.ValidateFieldName(AttriOrder_fld,os.path.dirname(inp_fc))
        NAttriOrder_fld = arcpy.ValidateFieldName(NAttriOrder_fld,os.path.dirname(inp_fc))
        Order_fld = arcpy.ValidateFieldName(Order_fld,os.path.dirname(inp_fc))
        

    ##Don't delete previous created order_fld, so that previous values stay.  
    ##    ##Get a list of all existing fields
    ##    allFlds = [fld.name.lower() for fld in arcpy.ListFields(inp_fc)]
    ##    if order_fld.lower() in allFlds:
    ##        msg = "Deleting the existing field " + order_fld 
    ##        PrintAndWriteMessages(msg,1)
    ##        arcpy.DeleteField(inp_fc,order_fld)
        
        #Add the double field
        result = AddAField(inp_fc,PeanoOrder_fld,"DOUBLE")
##        if result.MaxSeverity == 1:
##            msg = result.GetMessages(1).split(':')[1].lstrip() + "Existing values will be overwritten."
##            PrintAndWriteMessages(msg,1)
##        arcpy.SetProgressorPosition()
        
        result = AddAField(inp_fc,AttriOrder_fld,"DOUBLE")
##        if result.MaxSeverity == 1:
##            msg = result.GetMessages(1).split(':')[1].lstrip() + "Existing values will be overwritten."
##            PrintAndWriteMessages(msg,1)
##        arcpy.SetProgressorPosition()

        result = AddAField(inp_fc,NAttriOrder_fld,"DOUBLE")
##        if result.MaxSeverity == 1:
##            msg = result.GetMessages(1).split(':')[1].lstrip() + "Existing values will be overwritten."
##            PrintAndWriteMessages(msg,1)
##        arcpy.SetProgressorPosition()
        
        result = AddAField(inp_fc,Order_fld,"DOUBLE")
##        if result.MaxSeverity == 1:
##            msg = result.GetMessages(1).split(':')[1].lstrip() + "Existing values will be overwritten."
##            PrintAndWriteMessages(msg,1)
##        arcpy.SetProgressorPosition()

##        msg = "Computing the Peano curve coordinates for input features"
##        arcpy.SetProgressorLabel(msg)
##        PrintAndWriteMessages(msg,0)
        #Get the extent for the feature class
        desc = arcpy.Describe(inp_fc)
        arcpy.env.workspace = desc.path
    ##    #use the feature class extent even if the input is a feature layer
    ##    desc = arcpy.Describe(desc.CatalogPath)

        #Just in case, conver the input feature class layer to feature class, so the extent will be adjusted to the right values.
        #Reason:  a feature class layer's extent always points back to the reference feature class, in spite of the fact that the layer is
        # only a subset of all the features.
        theRealFC = "theRealFC"
        arcpy.FeatureClassToFeatureClass_conversion(inp_fc, arcpy.env.workspace, theRealFC)
        if desc.featureClass.dataType == 'ShapeFile':
            desc2 = arcpy.Describe(theRealFC + ".shp") 
        if desc.featureClass.dataType == 'FeatureClass':
            desc2 = arcpy.Describe(theRealFC)
        
        extent  = desc2.Extent
        xmin = extent.XMin
        ymin = extent.YMin
        xmax = extent.XMax
        ymax = extent.YMax

        #testing 
##        msg = "The current working extent of Peano curve calculation is " + str(xmin) + ", " + str(ymin) + ", " + str(xmax) + ", " + str(ymax) + ". "
##        PrintAndWriteMessages(msg,1)
        
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
        #row = rows.next()
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
            #row = rows.next()
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
        #row = rows.next()
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
            #row = rows.next()
            j = j + 1
        del rows

        #Normalize weighted aggregated attribute scores then calculate the final order values.
        curFlds = [AttriOrder_fld,NAttriOrder_fld, PeanoOrder_fld,Order_fld]
        rows = arcpy.da.UpdateCursor(inp_fc, curFlds)
        #rows.reset()
        #row = rows.next()
        for row in rows:
            #Scale the weighted aggregated atrtibute values to 0-1, the same scale as the Peano order.

            #NAggrScore = (row.getValue(AttriOrder_fld) - minScore) / (maxScore - minScore)
            #row.setValue(NAttriOrder_fld,NAggrScore)
            #CombOrder = (row.getValue(PeanoOrder_fld) * Pct_Peano + row.getValue(NAttriOrder_fld) * Pct_Attri) / 100
            #row.setValue(Order_fld, CombOrder)

            NAggrScore = (row[0] - minScore) / (maxScore - minScore)
            row[1] = NAggrScore
            CombOrder = (row[2] * Pct_Peano + row[1] * Pct_Attri) / 100
            row[3] = CombOrder
            
            rows.updateRow(row)
            arcpy.SetProgressorPosition()
            #row = rows.next()
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
### End of function



### Original ClusteringOder.py code ###
#Import the modules
import arcpy, sys,os, traceback
arcpy.env.overwriteOutput  = True

#This part is in the standalone python script but should be moved inside when using a function.
#Get the inputs
inp_fc = arcpy.GetParameterAsText(0)
PeanoOrder_fld = arcpy.GetParameterAsText(1)
AttriOrder_fld = arcpy.GetParameterAsText(2)
NAttriOrder_fld = arcpy.GetParameterAsText(3)
Order_fld = arcpy.GetParameterAsText(4)
Attri_list = arcpy.GetParameterAsText(5)
Wt_list = arcpy.GetParameterAsText(6)
Pct_Peano = float(arcpy.GetParameterAsText(7))
Pct_Attri = 100 - Pct_Peano

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
        
# use a field to copy FID/OID to the new integer field, because FID/OID cannot be used in GenerateSpatialWeightsMatrix
AddAField(inp_fc,"TheID","Short")
    
# A pair of "!" signs for using existing variables (attribute table) in Python expression.
oid_fld = arcpy.Describe(inp_fc).OIDFieldName
exp = "!" + oid_fld + "!"
arcpy.CalculateField_management(inp_fc, "TheID", exp, "PYTHON_9.3","#")

#ready to call the functions in other scripts
##import Func_ClusteringOrder
##Func_ClusteringOrder.ClusteringOrder \
##            (inp_fc, PeanoOrder_fld, AttriOrder_fld, NAttriOrder_fld, Order_fld, Attri_list, Wt_list, Pct_Peano)

ClusteringOrder \
            (inp_fc, PeanoOrder_fld, AttriOrder_fld, NAttriOrder_fld, Order_fld, Attri_list, Wt_list, Pct_Peano)

#The End


