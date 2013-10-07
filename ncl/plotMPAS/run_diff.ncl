;=============================================================================================
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
;=============================================================================================

begin

;SET RESOURCES:
res = True
res@gsnDraw             = True
res@gsnFrame            = False
res@gsnMaximize         = True
res@gsnSpreadColors     = True

;res@mpProjection        = "CylindricalEquidistant"
res@mpProjection        = "Orthographic"
res@mpDataBaseVersion   = "MediumRes" 
res@mpCenterLatF        = 90
res@mpCenterLonF        = 0.
res@cnFillOn            = False
;res@cnFillMode          = "AreaFill"
res@cnFillMode          = "RasterFill"
res@cnRasterSmoothingOn = True
res@cnFillOn            = True
res@cnLinesOn           = False

res@cnLevelSelectionMode = "ManualLevels"
res@cnMinLevelValF  =  -300.
res@cnMaxLevelValF  =  300.
res@cnLevelSpacingF =  50.

;; limited area plot domain
res@mpMinLatF =  60.
res@mpMaxLatF =  90.
res@mpMinLonF = -180.
res@mpMaxLonF =  180.

wks = gsn_open_wks("png","test_diff")
gsn_define_colormap(wks,"BlueYellowRed")

diri = "/arctic1/nick/cases/2007/"
files7 = systemfunc("ls -1 " + diri + "x4.t.output*")
diri = "/arctic1/nick/cases/vduda/x4/"
files6 = systemfunc("ls -1 " + diri + "x4.t.output.2006-08-01*")
nfiles = 1

set_default_fillvalue("double",  -99999)

stride = 1

ff = addfile(files6(0),"r")

r2d = 57.2957795 ; conversion from radians to degrees.
res@sfXArray = ff->lonCell(:) * r2d
latCell = ff->latCell(:)
res@sfYArray = latCell * r2d

delete(ff)

iFile = 0
ff6 = addfile(files6(iFile),"r")
ff7 = addfile(files7(iFile),"r")

nTimes = 29

avg6 = dim_sum_n(ff6->height_500hPa(0:nTimes-1,:),0)/nTimes
avg7 = dim_sum_n(ff7->height_500hPa(0:nTimes-1,:),0)/nTimes

field = avg6-avg7

res@lbTitleString="m"
res@lbTitlePosition= "Bottom"
res@lbTitleFontHeightF= 0.018
res@lbTitleDirection= "Across"

;field@long_name = "theta_2pvu "+ff->xtime(iTime,0:18)
;field@_FillValue = missingVal
;print(" plotting time "+ff->xtime(iTime,0:18))
;field@units = "[m]"
plot = gsn_csm_contour_map(wks,field,res)
draw(plot)
frame(wks)

delete(ff6)
delete(ff7)


end

;=============================================================================================
