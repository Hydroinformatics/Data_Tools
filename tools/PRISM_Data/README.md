# Python scripts used to pre-process **PRISM _Recent Years_** data.

### Dependencies:
  - osgeo, gdal, zipfile, numpy
 
### Data source: [PRISM Data](http://www.prism.oregonstate.edu/)

### Included scripts
 - PRISM_Data_Parser: Script can be used to summarize and convert raster PRISM data into a CSV file. Data from the PRISM zip raster files (i.e., .bil files), will be extracted for a given extent (i.e., region of interest) and saved to a CSV file. The script also generates a CSV with the centroid of each raster cell (i.e., raster cell ID, Lat, Long) for which data was obtained. 
 - PRISM_Data_Parser_Yearly_Totals: Creates a time series of the spatial yearly average of total precipitation for raster cells included within a watershed boundary. Similar to the PRISM_Data_Parser script, given a shapefile with the delineated watersheds, the script will produce a JSON file with a time series of the mean (spatial) yearly total precipitation for all raster cells included within each of the delineated watersheds. 
 
### Example use of PRISM_Data_Parser
#### Step 1: Download the Recent Data PRISM data
  Download the Recent Data PRISM data via FTP following the instructions detailed in http://www.prism.oregonstate.edu/documents/PRISM_downloads_FTP.pdf
  - Example input parameters for WinSCP: 
    - File protocol: FTP
    - Encryption: No encryption
    - Host: prism.nacse.org
    - Username: anonymous
    - Port: 21
    - Password: email address
    
  Save the downloaded data under the same directory.
  
#### Step 2: Define region of interest
  Clip one of the PRISM rasters to the extent of interest and save as a .tif file in the same directory of Step 1.
  
#### Step 3: Run script
  Call the PRISM_Data_Parser.py code from the command windows. Arguments include, in order:
  1. Path to data directory
  2. TIFF raster file clip to extent of interest
  3. Beginning year (e.g., 1981)
  4. End year (e.g., 2009)
  
#### Step 4: Review output files
  1. PRISM_PPT_Data.csv: Each row represents a time step (year, mon., day) and each column is associated with a raster cell (referenced by the cell ID).
  2. PRISM_Raster_Centroids.csv: Includes the raster cell IDs and the latitude and longitude of their centroid.  