import gdal, ogr, os, osr
import subprocess

def createFolder(directory):
    try:
        if not os.path.exists(directory):
            os.makedirs(directory)
    except OSError:
        print ('Error: Creating directory. ' +  directory)
        

#to get intersection, first need to create polygons
def polygon (imagepath):
    bag = gdal.Open(imagepath)  # replace it with your file
     # raster is projected                                          
    bag_gtrn = bag.GetGeoTransform()
    bag_proj = bag.GetProjectionRef()
    bag_srs = osr.SpatialReference(bag_proj)
    geo_srs =bag_srs.CloneGeogCS()                 # new srs obj to go from x,y -> f,?
    transform = osr.CoordinateTransformation( bag_srs, geo_srs)

    bag_bbox_cells = (
            (0., 0.),
            (0, bag.RasterYSize),
            (bag.RasterXSize, bag.RasterYSize),
            (bag.RasterXSize, 0),
            )

    geo_pts = []
    for x, y in bag_bbox_cells:
        x2 = bag_gtrn[0] + bag_gtrn[1] * x + bag_gtrn[2] * y
        y2 = bag_gtrn[3] + bag_gtrn[4] * x + bag_gtrn[5] * y
        geo_pts.append(x2)
        geo_pts.append(y2)

    mspoly = ogr.Geometry(ogr.wkbLinearRing)
    mspoly.AddPoint(geo_pts[0], geo_pts[1])
    mspoly.AddPoint(geo_pts[2], geo_pts[3])
    mspoly.AddPoint(geo_pts[4], geo_pts[5])
    mspoly.AddPoint(geo_pts[6], geo_pts[7])
    mspoly.AddPoint(geo_pts[0], geo_pts[1])

    poli1 = ogr.Geometry(ogr.wkbPolygon)
    poli1.AddGeometry(mspoly)
    return poli1

#to get intersection point coordinates
def intersect(poli_ms,poli_pan):
    wkt1 = str(poli_ms) #get "poli_ms" from polygon func
    wkt2 = str(poli_pan)

    poly1 = ogr.CreateGeometryFromWkt(wkt1)
    poly2 = ogr.CreateGeometryFromWkt(wkt2)
    intersection = poly1.Intersection(poly2)
    ring = intersection.GetGeometryRef(0)
    count = ring.GetPointCount()
    countlist=[]
    for i in range(0,count):
        countlist+=[list(ring.GetPoint(i))]        
    xdif=countlist[0][0]-countlist[3][0]
    ydif=countlist[0][1]-countlist[1][1]
    return xdif,ydif,countlist


def jp2_to_tif (imagepath,output):
     
         os.system('gdal_translate -ot UInt16 -of GTiff '+imagepath+' '+output+'trans.tif')

#Writes rasters from ReadRaster()
def array2raster(newRasterfn,oX,oY,pixelWidth,pixelHeight,proj,nob,image,i,j,tilesize,gdt2): #new filename, originx, originy, pixelwidth, pixelheight, projection, number of bands, byte array from ReadRaster(), tile's upper left coordinates(x)(image coordinate system), and y, tilesize, GetDataTypeByName()

    cols = tilesize
    rows = tilesize
    originX = oX
    originY = oY

    driver = gdal.GetDriverByName('GTiff')
    outRaster = driver.Create(newRasterfn, cols, rows, nob, gdal.GDT_UInt16)

    outRaster.SetGeoTransform((originX, pixelWidth, 0, originY, 0, pixelHeight))
    outRaster.WriteRaster(0,0, tilesize, tilesize,image,tilesize,tilesize,gdt2,range(1,nob+1),0,0,0)
    outRasterSRS = osr.SpatialReference()
    outRasterSRS.ImportFromEPSG(4326)
    outRaster.SetProjection(proj)
    outRaster.FlushCache()

#to get tiles' upper left coordinates
def origin(dif,geotrans,tilesize,i,j):
    x = geotrans[0] + geotrans[1]*i
    y = geotrans[3] + geotrans[5]*j
    return x,y

#
def func(imagepath,outpath,tilesize): #image path that you divide tiles, the path you want to write the tiled images, tilesize

    image = gdal.Open(imagepath) #opens image
    proj = image.GetProjection()  #gets projection features 
    geotrans = image.GetGeoTransform()  #gets upper left coordinates, spatial resolution
    pixel_x=int((dif[2][1][0]-geotrans[0])/geotrans[1]) #calculates intersections' upper left coordinates according to image pixel coordinate system 
    pixel_y=(int((dif[2][1][1]-geotrans[3])/geotrans[5]))
    
   
    
    
    nob = image.RasterCount #counts number of bands
    createFolder(outpath+"\\tiled") #creates a folder for tiled images
    band = image.GetRasterBand(1) #reads raster band for getting datatype
    gdtt = gdal.GetDataTypeName(band.DataType) #gets GDAL data type
    gdt2 = gdal.GetDataTypeByName(gdtt)  #gets WKT Raster pixel type 
    xsize = abs(dif[0]/geotrans[1]) #calculates number of pixels x and y 
    ysize = abs(dif[1]/geotrans[1])
    xrem = xsize % tilesize 
    yrem = ysize % tilesize
    x =int( xsize - xrem)  #number of pixels x and y without remainders
    y =int( ysize - yrem)
    bandlist=range(1,nob+1) #creates a band list with nob for ReadRaster()

    u=0 #u and v are variables just for naming the file
    v=0
    for i  in range (pixel_x, x, tilesize): #starts with upper left x and increases as much as tilesize (image pixel coordinate system)
        u=u+1
        for j in range (pixel_y, y, tilesize): #starts with upper left y and increases as much as tilesize (image pixel coordinate system) 
            dataraster_1 = image.ReadRaster(i,j,tilesize,tilesize,tilesize,tilesize,gdal.GDT_UInt16,tuple(bandlist)) #xoff, yoff, xsize, ysize, buf_xsize, buf_ysize, buf_type, band_list
            oX,oY = origin(dif,geotrans,tilesize,i,j) #intersection of polygons, self.GetGeoTransform(), tilesize, i, j
            array2raster(outpath+"\\tiled\\newRaster"+str(1)+"-"+str(u)+"-"+str(v)+".tif",oX,oY,geotrans[1],geotrans[5],proj,nob,dataraster_1,i,j,tilesize,gdt2)
            v=v+1

path_pan = "C:\\Users\\USER-1\\Desktop\\g\\PHR1A_ISTANBUL_20151224_PAN_ORT\\PHR1A_ISTANBUL_20151224_PAN_ORT\\IMG_PHR1A_P_001\\IMG_PHR1A_P_201512240903304_ORT_PHR1A_20180702_12090717g0otitmh09b_1_R1C1.JP2"      
path_ms = "C:\\Users\\USER-1\\Desktop\\g\\PHR1A_ISTANBUL_20151224_XS_ORT\\PHR1A_ISTANBUL_20151224_XS_ORT\\IMG_PHR1A_MS_001\\IMG_PHR1A_MS_201512240903304_ORT_PHR1A_20180702_1209161ogdvi94vm2rg_1_R1C1.JP2"
outpath_pan ="C:\\Users\\USER-1\\Desktop\\Uhuzam\\LandsatGörüntü\\pan"
outpath_ms="C:\\Users\\USER-1\\Desktop\\Uhuzam\\LandsatGörüntü\\ms"
poli1=polygon(path_pan)
poli2=polygon(path_ms)
dif=intersect(poli2,poli1)

func (path_pan,outpath_pan,1024)
func (path_ms,outpath_ms,256)