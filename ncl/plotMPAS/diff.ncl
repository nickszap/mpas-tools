;=============================================================================================
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
;=============================================================================================

external mpas_vInterp "./mpas_vInterp.so"

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

sres = True
sres@gsnDraw   = False
sres@gsnFrame  = False
sres@cnFillOn  = False
sres@cnLinesOn = True
sres@cnLineLabelsOn = False ;;True
sres@cnLineThicknessF = 3.0
sres@cnLineLabelFontHeightF = 0.009
sres@cnInfoLabelFontHeightF = 0.01
sres@cnLineLabelPlacementMode = "Computed"
sres@cnLineLabelDensityF = 1.2


res@cnLevelSelectionMode = "ManualLevels"
res@cnMinLevelValF  =  -10.
res@cnMaxLevelValF  =  10.
res@cnLevelSpacingF =  1.0


sres@cnLevelSelectionMode = "ManualLevels"
sres@cnMinLevelValF  =  280.
sres@cnMaxLevelValF  =  360.
sres@cnLevelSpacingF =  5.

wks = gsn_open_wks("png","mpasError")
gsn_define_colormap(wks,"BlueYellowRed")

diri = "/arctic1/nick/cases/vduda/x4/"
;files = systemfunc("ls -1 " + diri + "x4.t.output.*")
files = systemfunc("ls -1 " + diri + "x4.kf.output.2006-08-01*")
nfiles = dimsizes(files)

diri = "/arctic1/nick/cases/cfsr/output/"
files_anl = systemfunc("ls -1 "+ diri + "x4.cfsr.output.2006-08-08_00*")
ianl_start = 0 ; start index into files to compare with
iFile_anl = ianl_start
;level = 5

set_default_fillvalue("double",  -99999)

stride = 1

ff = addfile(files(0),"r")

r2d = 57.2957795 ; conversion from radians to degrees.
res@sfXArray = ff->lonCell(:) * r2d
res@sfYArray = ff->latCell(:) * r2d
sres@sfXArray = ff->lonCell(:) * r2d
sres@sfYArray = ff->latCell(:) * r2d
latCell = ff->latCell(:)
delete(ff)

;; loop over times

do nf = 0, nfiles(0)-1

print("processing file "+files(nf))

ff = addfile(files(nf),"r")

;;fieldin1 = ff->height_200hPa(:,:)
;;dims = dimsizes(fieldin1)
;;delete(fieldin1)

;nTime = dims(0) - 1
nTime=28
tStart = 28 ;; 0
print("number of times in file "+nTime)

do iTime=tStart,nTime
print("Reference file is "+files_anl(iFile_anl))
ffcfsr = addfile(files_anl(iFile_anl),"r")
iFile_anl = iFile_anl+1

;;reference field
cfsrTime = 0
epv = ffcfsr->ertel_pv(cfsrTime, :,:)
theta = ffcfsr->theta(cfsrTime,:,:)
dims = dimsizes(epv) ;how to access dimensions of file instead of reading variable?
nCells = dims(0)
nLevels = dims(1)

field_cfsr = new(nCells,double)

pvuVal = todouble(2.0)
missingVal = todouble(-99999.0)
mpas_vInterp::mpas_vinterp(nCells, nLevels, pvuVal, latCell, epv, theta, field_cfsr, missingVal)


;;forecast field
epv = ff->ertel_pv(iTime, :,:)
theta = ff->theta(iTime,:,:)
dims = dimsizes(epv) ;how to access dimensions of file instead of reading variable?
nCells = dims(0)
nLevels = dims(1)

field_mpas = new(nCells,double)

pvuVal = todouble(2.0)
missingVal = todouble(-99999.0)
mpas_vInterp::mpas_vinterp(nCells, nLevels, pvuVal, latCell, epv, theta, field_mpas, missingVal)

;;compute difference field
field = field_mpas-field_cfsr

sres@cnInfoLabelOn = False ;Turn off 'contour from x to y by z' that gets cut off

field@long_name = "mpas-cfsr theta_2pvu" + " "+ff->xtime(iTime,0:18)
;field@_FillValue = missingVal
print(" plotting time "+ff->xtime(iTime,0:18))
field@units = "[K]"
plot = gsn_csm_contour_map(wks,field,res)
plot_ov = gsn_csm_contour(wks, field_mpas, sres)
overlay(plot,plot_ov)
draw(plot)
frame(wks)

delete(ffcfsr)

end do

delete(ff)
delete(field_cfsr)
delete(field_mpas)
delete(field)
;delete(theta)
;delete(thetacfsr)

end do

end

;=============================================================================================
