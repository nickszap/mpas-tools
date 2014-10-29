;======================================================================
; based off of ESMF_regrid_10.ncl
;
;Generate weights for mpas -> lat/lon mesh

load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/esmf/ESMF_regridding.ncl"

begin
;---Retrieve data
    srcDir = "/glade/scratch/szapiro/mpasRuns/2006/midLat/run/"
    srcFileName = srcDir+"output.2006-08-01_00.00.00.nc"
    print("Weights file for "+srcFileName)
    sfile = addfile(srcFileName,"r")
    sp    = sfile->theta(0,:,0)

    lonCell = sfile->lonCell    ; Cell centers
    latCell = sfile->latCell

    r2d     = 180.0d/(atan(1)*4.0d)  ; Radian to Degree
    lonCell = lonCell*r2d
    latCell = latCell*r2d

    Opt                 = True                     ; Regridding optioins

    Opt@SrcFileName     = "midLat_ESMF.nc"           ; output files
    Opt@DstFileName     = "World_0.5deg_SCRIP.nc"
    Opt@WgtFileName     = "midLat_2_.5NH.nc"
    Opt@ForceOverwrite  = True

    Opt@SrcGridLat      = latCell                  ; source grid
    Opt@SrcGridLon      = lonCell
    Opt@InputFileName   = srcFileName              ; optional

    Opt@DstGridType     = ".5x.5"                ; destination grid
    Opt@DstTitle        = "NH Grid .5 degree"
    Opt@DstLLCorner     = (/0.00d,   0.00d /)
    Opt@DstURCorner     = (/ 90.0d, 359.5d /)

    Opt@Debug           = True
    ;;Opt@PrintTimings  = True

    sp_regrid = ESMF_regrid(sp,Opt)
    printVarSummary(sp_regrid)
end

