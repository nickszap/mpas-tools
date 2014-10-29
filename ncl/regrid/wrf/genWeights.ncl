;----------------------------------------------------------------------
; The following is adapted from:
; http://www.ncl.ucar.edu/Applications/Templates/ESMF_unstruct_to_curv.ncl
;
; This is an NCL/ESMF template file for regridding from an 
; unstructured grid to a curvilinear grid. It uses ESMF_regrid
; to do the regridding.
; Both grids are assumed to be contained in separate NetCDF files.
;
; For more information on ESMF_regrid, see: 
;    http://www.ncl.ucar.edu/Document/Functions/ESMF/ESMF_regrid.shtml
;----------------------------------------------------------------------
; This example uses the ESMF application "ESMF_RegridWeightGen" to 
; generate the weights.
; For more information about ESMF:
;        http://www.earthsystemmodeling.org/
; This script uses built-in functions that are only available in 
; NCL V6.1.0-beta and later.
;----------------------------------------------------------------------
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/esmf/ESMF_regridding.ncl"

begin

;Begin user options --------------------
src_file = "x4.arctic.output.nc" ;mpas filepath
dst_file = "wrf_domain_36km.nc" ;wrf filepath
wgt_file = "mpas2wrf.nc" ;filename of mapping weights we'll make
interpMethod = "patch" ;"bilinear", "patch", or (untested) "conservative"
;End user options ----------------------

;---Data file containing source grid
    sfile    = addfile(src_file,"r")
;---Get mesh and variable to regrid
    var_name = "theta_pv"
    var      = sfile->$var_name$(0,:)
    r2d = 57.2957795 ; conversion from radians to degrees.
    src_lat  = sfile->latCell * r2d
    src_lon  = sfile->lonCell * r2d

;---Data file containing destination grid
    dfile    = addfile(dst_file,"r")
    dst_lat  = dfile->XLAT(0,:,:)
    dst_lon  = dfile->XLONG(0,:,:)

;---Set up regridding options
    Opt                   = True

;---"bilinear" is the default. "patch" and "conserve" are other options.
    Opt@InterpMethod      = interpMethod
    Opt@WgtFileName       = wgt_file

    Opt@SrcGridLat        = src_lat           ; source grid
    Opt@SrcGridLon        = src_lon
    Opt@SrcRegional       = True
    Opt@SrcInputFileName  = src_file          ; optional, but good idea

    Opt@DstGridLat        = dst_lat         ; destination grid
    Opt@DstGridLon        = dst_lon
    Opt@DstRegional       = True
    
;;---Change (likely)
;;  Only needed for "conserve" interpolation.
;;  Opt@DstGridCornerLat  = ndtooned( dfile->latt_bounds )
;;  Opt@DstGridCornerLon  = ndtooned( dfile->lont_bounds )
 
    Opt@ForceOverwrite    = True
    Opt@PrintTimings      = True
    Opt@Debug             = True

    var_regrid = ESMF_regrid(var,Opt)     ; Do the regridding

    printVarSummary(var_regrid)

;----------------------------------------------------------------------
; Plotting section
;
; This section creates filled contour plots of both the original
; data and the regridded data, and panels them.
;----------------------------------------------------------------------
    wks = gsn_open_wks("png","unstruct_to_curv")

    res                       = True
    res@mpProjection          = "Orthographic"
    res@mpDataBaseVersion   = "MediumRes"
    res@mpCenterLatF        = 90
    res@mpCenterLonF        = 0.


    res@gsnMaximize           = True

    res@gsnDraw               = False
    res@gsnFrame              = False

    res@cnFillOn              = True
    res@cnLinesOn             = False
    res@cnLineLabelsOn        = False
    res@cnFillMode            = "RasterFill"

    res@lbLabelBarOn          = False    ; Turn on later in panel

    res@mpMinLatF             = min(dst_lat)
    res@mpMaxLatF             = max(dst_lat)
    res@mpMinLonF             = min(dst_lon)
    res@mpMaxLonF             = max(dst_lon)

;;--Change (maybe)
    ;mnmxint = nice_mnmxintvl( min(var), max(var), 18, False)
    res@cnLevelSelectionMode = "ManualLevels"
    res@cnMinLevelValF       = 280. ; mnmxint(0)
    res@cnMaxLevelValF       = 360. ; mnmxint(1)
    res@cnLevelSpacingF      = 5. ; mnmxint(2)

;---Resources for plotting regridded data
    res@gsnAddCyclic = False            ;;---Change (maybe)
    res@tiMainString  = "Curvilinear grid (" + Opt@InterpMethod + ")"

    plot_regrid = gsn_csm_contour_map(wks,var_regrid,res)

;---Resources for plotting original data
    res@gsnAddCyclic       = False            ;;---Change (maybe)
    res@sfXArray           = src_lon
    res@sfYArray           = src_lat
    res@tiMainString = "Original unstructured grid (" + \
                        dimsizes(src_lon) + " cells)"

    plot_orig = gsn_csm_contour_map(wks,var,res)  

;---Compare the plots in a panel
    pres                   = True
    pres@gsnMaximize       = True
    pres@gsnPanelLabelBar  = True

    gsn_panel(wks,(/plot_orig,plot_regrid/),(/2,1/),pres)
end
