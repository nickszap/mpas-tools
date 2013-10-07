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

res@cnLevelSelectionMode = "ManualLevels"
res@cnMinLevelValF  =  280.
res@cnMaxLevelValF  =  360.
res@cnLevelSpacingF =  5.

wks = gsn_open_wks("png","theta_2pvu_20070801_t_x4")
;wks = gsn_open_wks("png","theta_2pvu_20060807_kf_o3_x7")
;wks = gsn_open_wks("png","test_cfsr")
gsn_define_colormap(wks,"BlueYellowRed")

;diri = "/arctic1/nick/cases/v1.0/x4/longer/"
diri = "/arctic1/nick/cases/2007/"
;diri = "/arctic1/nick/cases/vduda/"
;diri = "/arctic1/nick/cases/vduda/x4/"
;diri = "/arctic1/nick/cases/cfsr/output/"
;files = systemfunc("ls -1 " + diri + "x4.kf.output.*")
files = systemfunc("ls -1 " + diri + "x4.t.output.*")
;files = systemfunc("ls -1 " + diri + "x1.t.output.*")
;files = systemfunc("ls -1 " + diri + "x4.kf.output.*")
;files = systemfunc("ls -1 " + diri + "x7.kf.o3*")
;files = systemfunc("ls -1 " + diri + "x4.cfsr.output.2006-08-16*")
;files = systemfunc("ls -1 " + diri + "x4.t.output.2006-07-24*")
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
;nTime=1
print("number of times in file "+nTime)

do iTime=0,nTime-1
epv = ff->ertel_pv(iTime, :,:)
theta = ff->theta(iTime,:,:)
dims = dimsizes(epv) ;how to access dimensions of file instead of reading variable?
nCells = dims(0)
nLevels = dims(1)

field = new(nCells,double)

pvuVal = todouble(2.0)
missingVal = todouble(-99999.0)
mpas_vInterp::mpas_vinterp(nCells, nLevels, pvuVal, latCell, epv, theta, field, missingVal)
;fortran is not case sensitive so use all lower case

field@long_name = "theta_2pvu "+ff->xtime(iTime,0:18)
field@_FillValue = missingVal
print(" plotting time "+ff->xtime(iTime,0:18))
field@units = "[K]"
plot = gsn_csm_contour_map(wks,field,res)
draw(plot)
frame(wks)

end do

delete(ff)
delete(field)
delete(epv)
delete(theta)

end do

delete(latCell)
end

;=============================================================================================
