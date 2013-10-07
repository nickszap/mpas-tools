;=============================================================================================
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
;=============================================================================================

external mpas_mslp "./mpas_mslp.so"

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
res@cnMinLevelValF  =  98000.
res@cnMaxLevelValF  =  102400.
res@cnLevelSpacingF =  400.

sres = True
sres@gsnDraw   = False
sres@gsnFrame  = False
sres@cnFillOn  = False
sres@cnLinesOn = True
sres@cnLineLabelsOn = False
sres@cnLineThicknessF = 2.0
sres@cnLineLabelFontHeightF = 0.009 
sres@cnInfoLabelFontHeightF = 0.01          
sres@cnLineLabelPlacementMode = "Computed"
sres@cnLineLabelDensityF = 1.2

sres@cnLevelSelectionMode = "ManualLevels"
sres@cnMinLevelValF  =  98000.
sres@cnMaxLevelValF  =  102400.
sres@cnLevelSpacingF =  400.

wks = gsn_open_wks("png","test_slp")
gsn_define_colormap(wks,"BlueYellowRed")

diri = "/arctic1/nick/cases/vduda/x4/"
files = systemfunc("ls -1 " + diri + "x4.t.output.2006-07-24*")
nfiles = dimsizes(files)

set_default_fillvalue("double",  -99999)

stride = 1

ff = addfile(files(0),"r")

r2d = 57.2957795 ; conversion from radians to degrees.
res@sfXArray = ff->lonCell(:) * r2d
res@sfYArray = ff->latCell(:) * r2d
sres@sfXArray = ff->lonCell(:) * r2d
sres@sfYArray = ff->latCell(:) * r2d

delete(ff)

;; loop over times

do nf = 0, nfiles(0)-1

print("processing file "+files(nf))

ff = addfile(files(nf),"r")
fieldin1 = ff->height_200hPa(:,:)
dims = dimsizes(fieldin1)
delete(fieldin1)

z0 = .5*(ff->zgrid(:,0)+ff->zgrid(:,1))

nTime = dims(0)
nCells = dims(1)

nTime=3
print("number of times in file "+nTime)

mslp = new(nCells,double)

do iTime=0,nTime-1
p0 = ff->pressure_p(iTime,:,0)+ff->pressure_base(iTime,:,0)
thetaSurf = ff->theta(iTime,:,0)

mpas_mslp::mpas_mslp(nCells, mslp, thetaSurf, z0, p0)

mslp@long_name = "MSLP "+ff->xtime(iTime,0:18)
print(" plotting time "+ff->xtime(iTime,0:18))
mslp@units = "[Pa]"
plot = gsn_csm_contour_map(wks,mslp,res)
plot_ov = gsn_csm_contour(wks,mslp,sres)
overlay(plot,plot_ov)
draw(plot)
;draw(plot_ov)
frame(wks)

end do

delete(ff)
delete(mslp)

end do

end

;=============================================================================================
