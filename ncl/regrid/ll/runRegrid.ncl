;Using created ESMF weights files, write a netcdf file of a variable on the destination grid

;=============================================================================================
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/esmf/ESMF_regridding.ncl"
;=============================================================================================

begin

external mpas_mslp "/glade/u/home/szapiro/scripts/mpas/ncl/mpas_mslp.so"

; ---------------- user options --------------------
;wgtFile = "/glade/scratch/szapiro/mpasRuns/2006/compare/regrid/arctic_2_.5NH.nc"
wgtFile = "midLat_2_.5NH.nc"
diri = "/glade/scratch/szapiro/mpasRuns/2006/midLat/run/"
files = systemfunc("ls -1 " + diri + "output*.nc")
nfiles = dimsizes(files)

outName = "midLat_mslp.nc"
nlat = 181 ;the grid we're interpolating to. There's probably a better way to set this.
nlon = 720
; ---------------- end user options --------------------

nf = 0
print("processing file "+files(nf))

ff = addfile(files(nf),"r")
fieldin1 = ff->height_200hPa(:,:)
dims = dimsizes(fieldin1)
delete(fieldin1)

z0 = .5*(ff->zgrid(:,0)+ff->zgrid(:,1))

nTime = dims(0)
nCells = dims(1)
print("number of times in file "+nTime)

; ----------------------create netcdf file ----------------------
system("rm "+outName)   ; remove any pre-existing file
fout     = addfile(outName ,"c")  ; open output netCDF file

    ; explicitly declare file definition mode. Improve efficiency.
        setfileoption(fout,"DefineMode",True)

    ; create global attributes of the file
        fAtt               = True            ; assign file attributes
	fAtt@title         = "Regridded variable"  
	fAtt@source_file   =  files(nf)
	fAtt@Conventions   = "None"   
	fAtt@creation_date = systemfunc ("date")        
	fileattdef( fout, fAtt )            ; copy file attributes    

    ; predefine the coordinate variables and their dimensionality
    ; Note: to get an UNLIMITED record dimension, we set the dimensionality
    ; to -1 (or the actual size) and set the dimension name to True.
        dimNames = (/"time", "lat", "lon"/)
	dimSizes = (/ -1   ,  nlat,  nlon /) 
	dimUnlim = (/ True , False, False/)   
	filedimdef(fout,dimNames,dimSizes,dimUnlim)

    ;===================================================================
    ; predefine the the dimensionality of the variables to be written out
    ;===================================================================
    ; Here we are using NCL functions to facilitate defining 
    ; each variable's dimension name(s) and type. 
    ; The following could be replaced with explicit, user defined dimension 
    ; names different from those associated with the variable in memory. 
    ; Say, PS(time,lat,lon) in the NCL script. They could be redefined for the file via: 
    ; filevardef(fout, "PS"   ,typeof(PS) ,(/"TIME","latitude","longitude"/)) 
    ;===================================================================                         
       filevardef(fout, "mslp"    , "double"  ,(/"time", "lat", "lon"/))

    ;===================================================================
    ; explicitly exit file definition mode. **NOT REQUIRED**
    ;===================================================================
        setfileoption(fout,"DefineMode",False)

Opt = False
mslp = new(nCells,double)
z0 = .5*(ff->zgrid(:,0)+ff->zgrid(:,1))

do iTime=0,nTime-1
p0 = ff->pressure_p(iTime,:,0)+ff->pressure_base(iTime,:,0)
thetaSurf = ff->theta(iTime,:,0)

mpas_mslp::mpas_mslp(nCells, mslp, thetaSurf, z0, p0)
mslp_regrid = ESMF_regrid_with_weights(mslp, wgtFile, Opt)


; output variables directly
fout->mslp(iTime,:,:)  = mslp_regrid

end do ;iTime
delete(mslp)

end

