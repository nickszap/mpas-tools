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

res@mpProjection        = "CylindricalEquidistant"
;res@mpProjection        = "Orthographic"
res@mpDataBaseVersion   = "MediumRes" 
res@mpCenterLatF        = 0. ; 90.
res@mpCenterLonF        = 0.
res@cnFillOn            = False
;res@cnFillMode          = "AreaFill"
res@cnFillMode          = "RasterFill"
res@cnRasterSmoothingOn = True
res@cnFillOn            = True
res@cnLinesOn           = False
;res@lbLabelPosition  = "Bottom"             ; label position
;res@gsnRightString  = ""  ;units are placed in this property I think

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

wks = gsn_open_wks("png","mesh_sqrtarea")
gsn_define_colormap(wks,"BlueYellowRed")

;files = systemfunc("find /arctic1/nick/cases/v1.0/x4/2week/kf/x4.kf.output*.nc -type f -print")

;files = (/ \
;"/arctic1/nick/cases/v1.0/x4/2week/kf/x4.kf.output.2006-07-24_12.00.00.nc"/)

files = (/ \
"/arctic1/nick/real/mesh/x6.380000.grid.nc"/)
;"/home/nickszap/research/mpas/mesh/scvt/grid.x6.163842.nc"/)

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
delete(ff)

;; loop over times

do nf = 0, nfiles(0)-1
;do nf = nfiles(0)-1, nfiles(0)-1

print("processing file "+files(nf))

ff = addfile(files(nf),"r")

fieldin = ff->areaCell(:) ; in m^2
fieldin = fieldin(:)/1.e6 ; in km^2
if (fieldin(0) .lt. 1.) then
  ;assume we were on unit sphere
  rEarth = 6371. ;in km
  rUnit = .001; in km
  rRatio = rEarth/rUnit
  fieldin = fieldin(:)*rRatio*rRatio
end if
fieldin = sqrt(fieldin(:)) ;equivalent length for square cell

printMinMax (fieldin, True)

dims = dimsizes(fieldin)

;nTime = dims(0)
nTime = 1

do iTime=0,nTime-1

field = fieldin(0::stride)
dimsf = dimsizes(field)

field@long_name = " "
field@_FillValue = -99999
;print(" plotting time "+ff->xtime(iTime,0:18))

field@units = "[km]"
;res@gsnBottomString = "[km]"

plot = gsn_csm_contour_map(wks,field,res)
draw(plot)
frame(wks)

end do

delete(ff)

end do

end

;=============================================================================================
