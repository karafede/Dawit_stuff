

###### extracting the variables from the ECMWF nc filess
###### saving the stack rasters as a tiff file with their time varibale in csv files on each folder 

# for exapmle the AOD estimates for 2013 are saved in a  single tiff with the time vector in csv file in 2013 folder.


source("Z:/_SHARED_FOLDERS/Air Quality/Phase 2/PhD_DG/GWR_V1/ECMWF model files/nc_ECMWF_fun.R")
import_nc_ECMWF(directory="Z:/_SHARED_FOLDERS/Air Quality/Phase 2/PhD_DG/GWR_V1/ECMWF model files", year_input=2013 )
import_nc_ECMWF(directory="Z:/_SHARED_FOLDERS/Air Quality/Phase 2/PhD_DG/GWR_V1/ECMWF model files", year_input=2014 )
import_nc_ECMWF(directory="Z:/_SHARED_FOLDERS/Air Quality/Phase 2/PhD_DG/GWR_V1/ECMWF model files", year_input=2015 )
import_nc_ECMWF(directory="Z:/_SHARED_FOLDERS/Air Quality/Phase 2/PhD_DG/GWR_V1/ECMWF model files", year_input=2016 )

