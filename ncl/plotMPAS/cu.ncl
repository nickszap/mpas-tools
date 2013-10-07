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
res@cnMinLevelValF  =  0.
res@cnMaxLevelValF  =  0.03
res@cnLevelSpacingF =  .003

wks = gsn_open_wks("png","cuColumn_2006080718_kf_x7_day")
gsn_define_colormap(wks,"BlueYellowRed")

;diri = "/arctic1/nick/cases/v1.0/x4/longer/"
;diri = "/arctic1/nick/cases/cfsr/output/"
diri = "/arctic1/nick/cases/vduda/x7/"
;files = systemfunc("ls -1 " + diri + "x4.t.output.2006-07-24*")
;files = systemfunc("ls -1 " + diri + "x7.t.output.2006-08-03_12.*")
files = systemfunc("ls -1 " + diri + "x7.kf.output.2006-08-07_1*")
nfiles = dimsizes(files)

set_default_fillvalue("double",  -99999)

stride = 1

ff = addfile(files(0),"r")

r2d = 57.2957795 ; conversion from radians to degrees.
res@sfXArray = ff->lonCell(:) * r2d
latCell = ff->latCell(:)
res@sfYArray = latCell * r2d

delete(ff)

;; loop over times

do nf = 0, nfiles(0)-1

print("processing file "+files(nf))

ff = addfile(files(nf),"r")
fieldin1 = ff->height_200hPa(:,:)
dims = dimsizes(fieldin1)
delete(fieldin1)

nTime = dims(0)
nCells = dims(1)
field = new(nCells,double)
;nTime=1
print("number of times in file "+nTime)

startTime = 3
do iTime=startTime,nTime-1
;sum over a day if have enough to do so. assume 6 hourly.
nTimesBack = min((/iTime, startTime/))
rthcuten = ff->rthcuten(iTime-nTimesBack:iTime, :,:)

;field = new(nCells,double)
field1 = dim_sum(rthcuten) ; sum rightmost dim...levels
field = dim_sum_n(field1, 0) ; sum over times

field@long_name = "rthcuten column "+ff->xtime(iTime,0:18)
;field@_FillValue = missingVal
print(" plotting time "+ff->xtime(iTime,0:18))
field@units = "[K/s]"
plot = gsn_csm_contour_map(wks,field,res)
draw(plot)
frame(wks)

end do

delete(ff)
;delete(field)

end do
delete(field)
delete(latCell)
end

;=============================================================================================
