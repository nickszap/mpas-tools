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
res@mpCenterLatF        = 90.
res@mpCenterLonF        = 0.
res@cnFillOn            = False
;res@cnFillMode          = "AreaFill"
res@cnFillMode          = "RasterFill"
res@cnRasterSmoothingOn = True
res@cnFillOn            = True
res@cnLinesOn           = False

;; limited area plot domain
;;res@mpMinLatF = -90.
;;res@mpMaxLatF =  90
;;res@mpMinLonF = -180.
;;res@mpMaxLonF =  180.


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
sres@cnLineLabelDensityF = 1.4


res@cnLevelSelectionMode = "ManualLevels"
res@cnMinLevelValF  =  -0.1
res@cnMaxLevelValF  =  2.1
res@cnLevelSpacingF =  0.1

sres@cnLevelSelectionMode = "ManualLevels"
sres@cnMinLevelValF  =  97000.
sres@cnMaxLevelValF  =  103400.
sres@cnLevelSpacingF =  400.

wks = gsn_open_wks("png","land_ice_sea_mslp")
gsn_define_colormap(wks,"BlueYellowRed")

;files = systemfunc("find /arctic1/nick/cases/v1.0/x4/2week/kf/x4.kf.output*.nc -type f -print")

files = (/ \
"/arctic1/nick/cases/vduda/x4/x4.t.output.2006-07-24_12.00.00.nc"/)

nfiles = dimsizes(files)
set_default_fillvalue("double",  -99999)

stride = 1

;; pull in some static data from the first output file

ff = addfile(files(0),"r")
r2d = 57.2957795 ; conversion from radians to degrees.
res@sfXArray = ff->lonCell(0::stride) * r2d
res@sfYArray = ff->latCell(0::stride) * r2d
sres@sfXArray = ff->lonCell(0::stride) * r2d
sres@sfYArray = ff->latCell(0::stride) * r2d
terrain_hgt_in = ff->zgrid(:,0)
terrain_hgt = terrain_hgt_in(0::stride)
delete(ff)

;; loop over times

do nf = 0, nfiles(0)-1
;do nf = nfiles(0)-1, nfiles(0)-1

print("processing file "+files(nf))

ff = addfile(files(nf),"r")

fieldin = ff->xland(:,:)
dims = dimsizes(fieldin)

nTime = dims(0)
nTime = 3

nCells = dims(1)
z0 = .5*(ff->zgrid(:,0)+ff->zgrid(:,1))
mslp = new(nCells,double)

do iTime=0,nTime-1
;mslp
p0 = ff->pressure_p(iTime,:,0)+ff->pressure_base(iTime,:,0)
thetaSurf = ff->theta(iTime,:,0)
mpas_mslp::mpas_mslp(nCells, mslp, thetaSurf, z0, p0)

field = fieldin(iTime,:)-ff->xice(iTime,:)
sres@cnInfoLabelOn = False

field@long_name = "sfc type, MSLP "+ff->xtime(iTime,0:18)
print(" plotting time "+ff->xtime(iTime,0:18))
field@units = " "
plot = gsn_csm_contour_map(wks,field,res)
plot_ov = gsn_csm_contour(wks,mslp,sres)
overlay(plot,plot_ov)
draw(plot)
frame(wks)

end do

delete(ff)
delete(mslp)

end do

end

;=============================================================================================
