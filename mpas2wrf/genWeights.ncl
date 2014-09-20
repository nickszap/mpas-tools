load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/wrf/WRFUserARW.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/esmf/ESMF_regridding.ncl"

;======================================================================
; Generate MPAS to WRF regrid weights
;
;Following templates in:
; http://www.ncl.ucar.edu/Applications/Scripts/ESMF_all_10.ncl
; http://www.ncl.ucar.edu/Applications/Scripts/ESMF_wgts_29.ncl
;======================================================================
    
; User options -----------------------
    InterpMethod= "patch"
    ;existing files
    mpasPath = "mpas.nc"
    wrfPath = "wrfout_d01_2003-07-13_12:00:00.UV.nc"
    
    ;files we'll make
    srcGridName = "MPAS_ESMF.nc"
    dstGridName = "World_0.25deg_SCRIP.nc"
    wgtFile     = "MPAS_2_World.nc"

    ;set to True if you want to skip any of these steps
    SKIP_mpas_GEN  = False
    SKIP_wrf_GEN = False
    SKIP_WGT_GEN   = False
; end user options -----------------------

;----------------------------------------------------------------------
; Convert MPAS to unstructured ESMF file.
;----------------------------------------------------------------------

;---Retrieve cell centers
    sfile   = addfile(srcFileName,"r")
    lonCell = sfile->lonCell
    latCell = sfile->latCell

;---Convert to degrees from radians
    r2d     = 180.0d/(atan(1)*4.0d)  ; Radian to Degree
    lonCell = lonCell*r2d
    latCell = latCell*r2d

    if(.not.SKIP_ESMF_GEN) then
;---Set some Options
      Opt                = True
      Opt@ForceOverwrite = True
      Opt@PrintTimings   = True
      Opt@InputFileName  = srcFileName
      Opt@Debug          = True

      unstructured_to_ESMF(srcGridName,latCell,lonCell,Opt)

;---Clean up
      delete(Opt)
    end if

;---Retrieve either one level, or all levels. Use '-1' for all.
    sfile = addfile(srcFilePath,"r")

    Opt                = True
    Opt@SrcTitle       = "WRF grid"   ; optional
    Opt@WgtFileName    = "WRF_to_Rect.WgtFile_"+InterpMethod+".nc"

;---Generate the names of the files containing the 
;   source and destination grid descriptions

    Opt@SrcFileName    = "WRF.SCRIP_grid_description.nc"         ; Name of source and
    Opt@DstFileName    = "Rectilinear.SCRIP_grid_description.nc" ; destination files
    Opt@SrcRegional    = True

;---Get the variable to be regridded               
    UM   = wrf_user_getvar(sfile,"ua",0)     ; On mass grid 
    printVarSummary(UM)                      ; (bottom_top,south_north,west_east)

    ua   = UM(0,:,:)                         ; (south_north,west_east)
                                             ; only need 'horizontal' spatial points 
    printVarSummary(ua)                      ; (bottom_top,south_north,west_east)
    print("ua: min="+min(ua)+"   max="+max(ua))

    delete(UM)                               ; no longer needed

;---Get the source lat/lon grid (mass grid)
    lat2d = sfile->XLAT(0,:,:)               ; (south_north,west_east)
    lon2d = sfile->XLONG(0,:,:)

    ua@lat2d = lat2d                         ; This information will be used by
    ua@lon2d = lon2d                         ; ESMF_regrid for the source grid

    dims  = dimsizes(lat2d)
    nlat  = dims(0)
    nlon  = dims(1)
