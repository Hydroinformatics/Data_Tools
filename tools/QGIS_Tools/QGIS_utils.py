# -*- coding: utf-8 -*-

import re, os, numpy, shutil
from osgeo import gdal, ogr, osr

#%%
# Finds the top (first) and bottom (last) rows that have "Data" in a raster. 
# Assumes that the NoData cells in the raster have a value of 0.
def Raster_row_boundaries(raster):
    
    row_sum_index = numpy.where(numpy.sum(raster, axis=1) != 0)
    top_row = row_sum_index[0][0]
    last_row = row_sum_index[0][-1]
    
    if top_row < 0:
        top_row = 0
    if last_row < 0:
        last_row = 0
        
    return top_row, last_row

# Finds the first (left) and last (right) columns that have "Data" in a raster
# Assumes that the NoData cells in the raster have a value of 0.
def Raster_col_boundaries(raster):
    
    col_sum_index = numpy.where(numpy.sum(raster,axis=0) != 0)
    left_col = col_sum_index[0][0]
    right_col = col_sum_index[0][-1]
    
    if left_col < 0:
        left_col = 0
    if right_col < 0:
        right_col = 0
        
    return left_col, right_col

#%%
def Read_Raster(raster_file, raster_NoData = -999):
    
    raster_ds = gdal.Open(raster_file, gdal.GA_ReadOnly)
    raster_NoData = raster_ds.GetRasterBand(1).GetNoDataValue()
    raster = raster_ds.GetRasterBand(1).ReadAsArray()     
    
    return raster, raster_NoData, raster_ds


    def ReadBilFile(bil):
        
        gdal.GetDriverByName('EHdr').Register()
        img = gdal.Open(bil)
        band = img.GetRasterBand(1)
        data = band.ReadAsArray()
        
        return data


#%%
def Save_Raster(array, ds, old_raster, newRaster_name, raster_NoData):
    
    [cols, rows] = old_raster.shape
    driver = gdal.GetDriverByName('GTiff')
    outRaster = driver.Create(newRaster_name, rows, cols, 1, gdal.GDT_Float64)
    outRaster.SetGeoTransform(ds.GetGeoTransform())
    outRaster.SetProjection(ds.GetProjection())
    outRaster.GetRasterBand(1).WriteArray(array)
    outRaster.GetRasterBand(1).SetNoDataValue(raster_NoData)
    outRaster.FlushCache()
    outRaster = None
    ds=None
    
    return

#%%
def Clip_Raster(InputRaster, OutputImage, RefImage, resample_method):
    
    InputSrc = gdal.Open(InputRaster, gdal.GA_ReadOnly)
    gdalformat = 'GTiff'
    datatype = InputSrc.GetRasterBand(1).DataType
    NoData_value = InputSrc.GetRasterBand(1).GetNoDataValue()
    if NoData_value == None:
        NoData_value = -999999
    ##########################################################
    RasterSrc = gdal.Open(RefImage, gdal.GA_ReadOnly)
    upx, xres, xskew, upy, yskew, yres = RasterSrc.GetGeoTransform()
    cols = RasterSrc.RasterXSize
    rows = RasterSrc.RasterYSize
 
    ulx = upx + 0*xres + 0*xskew
    uly = upy + 0*yskew + 0*yres
 
    llx = upx + 0*xres + rows*xskew
    lry = upy + 0*yskew + rows*yres
 
    lrx = upx + cols*xres + rows*xskew
    ury = upy + cols*yskew + rows*yres
 
    urx = upx + cols*xres + 0*xskew
    lly = upy + cols*yskew + 0*yres
    ##########################################################
    #RasterSrcClip = gdal.Open(RefImage, gdal.GA_ReadOnly)
    Projection = RasterSrc.GetProjectionRef()
    
    if lly > ury:
        old_ury = ury
        ury = lly
        lly = old_ury
        
    warp_opts = gdal.WarpOptions(
            format = gdalformat,
            outputType = datatype, 
            outputBounds = [llx, lly, urx, ury], 
            xRes = xres, 
            yRes = yres, 
            dstSRS = Projection, 
            resampleAlg = resample_method)  
    # resampleAlg = gdal.GRA_NearestNeighbour)

    OutTile = gdal.Warp(OutputImage, InputRaster, options=warp_opts)
    OutTile = None # Close dataset

    return

#%%
def Convert_Vector_To_Tiff(InputVector, OutputImage, RefImage, attribute_str = '', burnVal = 1, RefInputVector = '', tiffExtent = '', datatype = gdal.GDT_Byte):

    gdalformat = 'GTiff'
    #burnVal = 1 #value for the output image pixels
    ##########################################################
    # Get projection info from reference image
    #print RefImage
    Image = gdal.Open(RefImage, gdal.GA_ReadOnly)
    tiffgeo_transform = Image.GetGeoTransform()
    pixel_width = tiffgeo_transform[1]
    
    # Open Shapefile
    Shapefile = ogr.Open(InputVector)
    Shapefile_layer = Shapefile.GetLayer()
    geo_transform = Shapefile_layer.GetExtent()
    
    if RefInputVector != '':
        temp_Shapefile = ogr.Open(RefInputVector)
        temp_Shapefile_layer = temp_Shapefile.GetLayer()
        geo_transform = temp_Shapefile_layer.GetExtent()

    if tiffExtent != '':        
        x_min = tiffgeo_transform[0]
        y_max = tiffgeo_transform[3]
        x_max = x_min + tiffgeo_transform[1] * Image.RasterXSize
        y_min = y_max + tiffgeo_transform[5] * Image.RasterYSize
        
        #x_max = tiffgeo_transform[1]
        #y_min = tiffgeo_transform[2]
    
        x_res = int(numpy.round((x_max - x_min)/pixel_width))
        y_res = int(numpy.round((y_max - y_min)/pixel_width))
        
    else:
        x_min = geo_transform[0]
        y_max = geo_transform[3]
        #x_max = x_min + geo_transform[1] * Image.RasterXSize
        #y_min = y_max + geo_transform[5] * Image.RasterYSize
        x_max = geo_transform[1]
        y_min = geo_transform[2]
    
        x_res = int(numpy.round((x_max - x_min)/pixel_width))
        y_res = int(numpy.round((y_max - y_min)/pixel_width))

    if tiffgeo_transform[0] > x_min or tiffgeo_transform[3] < y_max:
        msg = 'Spatial extent of DEM is smaller than subbasin shapefile. RefRaster: ' + str(OutputImage) + ' Xlims: ' + str(tiffgeo_transform[0]) + '>' + str(x_min) + 'Ylims: ' + str(tiffgeo_transform[3]) + '<' + str(y_max)
        print(msg)
        
    # Rasterise
    #print x_res, y_res
    Output = gdal.GetDriverByName(gdalformat).Create(OutputImage, x_res, y_res, 1, datatype, options=['COMPRESS=DEFLATE'])
    Output.SetProjection(Image.GetProjectionRef())
    Output.SetGeoTransform((x_min, pixel_width, 0, y_min, 0, pixel_width))

    NoData_value = -999

    # Write data to band 1
    Band = Output.GetRasterBand(1)
    Band.Fill(NoData_value)
    Band.SetNoDataValue(NoData_value)
    Band.FlushCache()
    
    if attribute_str == '':
        gdal.RasterizeLayer(Output, [1], Shapefile_layer, burn_values=[burnVal])
    else:
        gdal.RasterizeLayer(Output, [1], Shapefile_layer, options=["ATTRIBUTE=" + attribute_str])

    # Close datasets
    Band = None
    Output = None
    Image = None
    Shapefile = None

    return    
#%%
def Get_Centroids_Raster(raster):
    upx, xres, xskew, upy, yskew, yres = raster.GetGeoTransform()
    cols = raster.RasterXSize
    rows = raster.RasterYSize
    
    centroids_X = np.ones((1,cols))*-999
    centroids_Y = np.ones((rows,1))*-999
    
    for j in range(rows):
        if j == 0:
            centroids_Y[j,0] = upy - np.abs(yres)/2.0
        else:
            centroids_Y[j,0] = centroids_Y[j-1,0] - np.abs(yres) 
            
    for i in range(cols):
        if i == 0:
            centroids_X[0,i] = upx + np.abs(xres)/2.0 
        else:
            centroids_X[0,i] = centroids_X[0,i-1] + np.abs(xres)
    
    return centroids_X, centroids_Y


#%%
    
#%%
def Copy_Vector(basefile,copyfile):
    
    shutil.copy2(basefile, copyfile)
    
    exist = os.path.isfile(basefile[:-4] +'.cpg')
    if exist:
        shutil.copy2(basefile[:-4] +'.cpg', copyfile[:-4] + '.cpg')
    
    exist = os.path.isfile(basefile[:-4] +'.dbf')
    if exist:
        shutil.copy2(basefile[:-4] +'.dbf', copyfile[:-4] + '.dbf')
    
    exist = os.path.isfile(basefile[:-4] +'.prj')
    if exist:
        shutil.copy2(basefile[:-4] +'.prj', copyfile[:-4] + '.prj')
    
    exist = os.path.isfile(basefile[:-4] +'.shx')
    if exist:
        shutil.copy2(basefile[:-4] +'.shx', copyfile[:-4] + '.shx')

#%%
def Reproject_Vector(input_layer, output_name, spatial_ref):
    driver = ogr.GetDriverByName('ESRI Shapefile')
    # get the input layer
    dataSource = driver.Open(input_layer, 0)
    inLayer = dataSource.GetLayer()
    # input SpatialReference
    inSpatialRef = inLayer.GetSpatialRef()
    
    # output SpatialReference
    outSpatialRef = spatial_ref
    
    # create the CoordinateTransformation
    coordTrans = osr.CoordinateTransformation(inSpatialRef, outSpatialRef)
        
    # create the output layer
    if os.path.exists(output_name):
        driver.DeleteDataSource(output_name)
        #os.remove(outputShapefile)
        
    outDataSet = driver.CreateDataSource(output_name)
    
    # add fields
    inLayerDefn = inLayer.GetLayerDefn()
    outLayer = outDataSet.CreateLayer(inLayer.GetName(), geom_type=inLayerDefn.GetGeomType())
    for i in range(0, inLayerDefn.GetFieldCount()):
        fieldDefn = inLayerDefn.GetFieldDefn(i)
        outLayer.CreateField(fieldDefn)
    
    # get the output layer's feature definition
    outLayerDefn = outLayer.GetLayerDefn()
    
    # loop through the input features
    inFeature = inLayer.GetNextFeature()
    while inFeature:
        # get the input geometry
        geom = inFeature.GetGeometryRef()
        # reproject the geometry
        geom.Transform(coordTrans)
        # create a new feature
        outFeature = ogr.Feature(outLayerDefn)
        # set the geometry and attribute
        outFeature.SetGeometry(geom)
        for i in range(0, outLayerDefn.GetFieldCount()):
            outFeature.SetField(outLayerDefn.GetFieldDefn(i).GetNameRef(), inFeature.GetField(i))
        # add the feature to the shapefile
        outLayer.CreateFeature(outFeature)
        # dereference the features and get the next input feature
        outFeature = None
        inFeature = inLayer.GetNextFeature()
        
    # Create projection file
    proj_outFileName = os.path.dirname(output_name) + '/' + os.path.basename(output_name)[:-4] + '.prj'
    with open(proj_outFileName, 'w') as prj_file:
        prj_file.write(str(spatial_ref.ExportToWkt()))
    prj_file.close()
    
    # Save and close the shapefiles
    dataSource = None
    outDataSet = None        
    
def Dissolve_Polygons(input_layer_file, output_name):

        driver = ogr.GetDriverByName("ESRI Shapefile")
        
        dataSource = driver.Open(input_layer_file, 0)
        inlayer = dataSource.GetLayer()
        in_spatialRef = inlayer.GetSpatialRef()
        
        # create the output layer
        outputShapefile = output_name
        if os.path.exists(outputShapefile):
            driver.DeleteDataSource(outputShapefile)
    
        outDataSet = driver.CreateDataSource(outputShapefile)
    
        # add fields
        inLayerDefn = inlayer.GetLayerDefn()
        outLayer = outDataSet.CreateLayer(inlayer.GetName(), geom_type=inLayerDefn.GetGeomType())
        for i in range(0, inLayerDefn.GetFieldCount()):
            fieldDefn = inLayerDefn.GetFieldDefn(i)
            outLayer.CreateField(fieldDefn)
    
        # get the output layer's feature definition
        outLayerDefn = outLayer.GetLayerDefn()
    
        # loop through the input features
        #inFeature = bnd_layer.GetNextFeature()
        cpoly=1
        #while inFeature:
        for inFeature in inlayer:
            # get the input geometry
            geom = ogr.CreateGeometryFromWkt(str(inFeature.GetGeometryRef()))
            # reproject the geometry
            if cpoly == 1:
                poly = geom
            else:
                poly = poly.Union(geom)
            cpoly = cpoly + 1
            #inFeature = bnd_layer.GetNextFeature()
                
        # create a new feature
        outFeature = ogr.Feature(outLayerDefn)
        # set the geometry and attribute
        outFeature.SetGeometry(poly)
        for i in range(0, outLayerDefn.GetFieldCount()):
            outFeature.SetField(outLayerDefn.GetFieldDefn(i).GetNameRef(), inFeature.GetField(i))
        # add the feature to the shapefile
        outLayer.CreateFeature(outFeature)
        # dereference the features and get the next input feature
        outFeature = None
        
        proj_outFileName = os.path.dirname(output_name) + '/' + os.path.basename(output_name)[:-4] + '.prj'
            
        with open(proj_outFileName, 'w') as prj_file:
            prj_file.write(str(in_spatialRef.ExportToWkt()))
        prj_file.close()
        
        dataSource = None
        outDataSet = None
    

#def Buffer_Polygon(input_layer_file, output_name, bufferdist):
#    
#    driver = ogr.GetDriverByName("ESRI Shapefile")
#    dataSource = driver.Open(input_layer_file)
#    inlayer = dataSource.GetLayer()
#
#    if os.path.exists(output_name):
#        driver.DeleteDataSource(output_name)
#    outputBuffer = driver.CreateDataSource(output_name)
#    bufferlyr = outputBuffer.CreateLayer(output_name, geom_type=ogr.wkbPolygon)
#    featureDefn = bufferlyr.GetLayerDefn()
#
#    for feature in inlayer:
#        ingeom = feature.GetGeometryRef()
#        geomBuffer = ingeom.Buffer(bufferdist)
#
#        outFeature = ogr.Feature(featureDefn)
#        outFeature.SetGeometry(geomBuffer)
#        bufferlyr.CreateFeature(outFeature)
#    
#    dataSource = None
#    outputBuffer = None


def Clip_Vector(input_layer_file, bnd_layer_file, clip_name = None, bufferdist = None):
        
    driver = ogr.GetDriverByName("ESRI Shapefile")
    
    dataSource = driver.Open(bnd_layer_file, 0)
    bnd_layer = dataSource.GetLayer()
    bnd_spatialRef = bnd_layer.GetSpatialRef()
    
    if bnd_layer.GetFeatureCount() > 1:        
        # create the output layer
        outputShapefile = os.path.dirname(input_layer_file) + '/temp_bnd_layer.shp'
        Dissolve_Polygons(bnd_layer_file, outputShapefile)
        
        bnd_layer_file = outputShapefile 
    
        dataSource = driver.Open(bnd_layer_file, 0)
        bnd_layer = dataSource.GetLayer()
        bnd_spatialRef = bnd_layer.GetSpatialRef()
        
        
    for bndfeat in bnd_layer:
        bnd_geom = bndfeat.geometry()
    
    dataSource = driver.Open(input_layer_file, 0)
    input_layer = dataSource.GetLayer()
    input_spatialRef = input_layer.GetSpatialRef()
    
    if input_spatialRef != bnd_spatialRef:
        outShapefile = os.path.dirname(input_layer_file) + '/' + os.path.basename(input_layer_file)[:-4] + '_proj.shp'
        Reproject_Vector(input_layer_file, outShapefile, bnd_spatialRef)
        
        input_layer_file = outShapefile
        dataSource = driver.Open(input_layer_file, 0)
        input_layer = dataSource.GetLayer()
        input_spatialRef = input_layer.GetSpatialRef()
    
    if clip_name:
        outShapefile = os.path.dirname(input_layer_file) + '/' + clip_name
    else:    
        outShapefile = os.path.dirname(input_layer_file) + '/' + os.path.basename(input_layer_file)[:-4] + '_clip.shp'
   
    # Remove output shapefile if it already exists
    if os.path.exists(outShapefile):
        driver.DeleteDataSource(outShapefile)
        #os.remove(outShapefile)

    inLayerDefn = input_layer.GetLayerDefn()
    outDs = driver.CreateDataSource(outShapefile)
    outLayer = outDs.CreateLayer(input_layer.GetName(), geom_type=inLayerDefn.GetGeomType())
    print 'Clipping ' + str(outShapefile)
    
    #input_layer.Clip(bnd_layer, outLayer,options = ['SKIP_FAILURES = YES'])
    
    # add fields
    idField = ogr.FieldDefn("OBJECTID", ogr.OFTInteger)
    outLayer.CreateField(idField)
    for i in range(0, inLayerDefn.GetFieldCount()):
        fieldDefn = inLayerDefn.GetFieldDefn(i)
        outLayer.CreateField(fieldDefn)
    
    # get the output layer's feature definition
    outLayerDefn = outLayer.GetLayerDefn()
    
    # loop through the input features
    countf = 1
    for inFeature in input_layer:
        # get the input geometry
        geom = inFeature.GetGeometryRef()
                
        if geom.Intersects(bnd_geom) or geom.Within(bnd_geom):
            # create a new feature
            outFeature = ogr.Feature(outLayerDefn)
            # set the geometry and attribute
            outFeature.SetGeometry(geom)
            for i in range(1, outLayerDefn.GetFieldCount()):
                tfield_name = outLayerDefn.GetFieldDefn(i).GetNameRef()
                tfield_val = inFeature.GetField(tfield_name)
                outFeature.SetField(tfield_name,tfield_val)
                
            outFeature.SetField('OBJECTID',int(countf))
            # add the feature to the shapefile
            outLayer.CreateFeature(outFeature)
            # dereference the features and get the next input feature
            outFeature = None
            #inFeature = input_layer.GetNextFeature()
            countf = countf + 1
        
    if clip_name:
        proj_outFileName = os.path.dirname(input_layer_file) + '/' + clip_name[:-4] + '.prj'
    else:
        proj_outFileName = os.path.dirname(input_layer_file) + '/' + os.path.basename(input_layer_file)[:-4] + '_clip' + '.prj'
        
    with open(proj_outFileName, 'w') as prj_file:
        prj_file.write(str(input_spatialRef.ExportToWkt()))
    prj_file.close()
    
    dataSource = None
    
    return outShapefile