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
res@mpCenterLatF        = 90.
res@mpCenterLonF        = 0.
res@cnFillOn            = False
res@cnFillMode          = "AreaFill"
;res@cnFillMode          = "RasterFill"
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
sres@cnLineLabelsOn = True
sres@cnLineThicknessF = 2.0
sres@cnLineLabelFontHeightF = 0.009 
sres@cnInfoLabelFontHeightF = 0.01          
sres@cnLineLabelPlacementMode = "Computed"
sres@cnLineLabelDensityF = 1.2


;; 200 hPa plots
;; temperature plot params
res@cnLevelSelectionMode = "ManualLevels"
res@cnMinLevelValF  =  280.
res@cnMaxLevelValF  =  360.
res@cnLevelSpacingF =  5.

wks = gsn_open_wks("png","cfsr_2pvu_2006")
gsn_define_colormap(wks,"BlueYellowRed")

;files = (/ \
;"/arctic1/nick/cases/cfsr/press/nomads/pgbhnl.gdas.2006081500.grb2.mb.nc"/)
;files = (/ \
;"/arctic1/nick/cases/cfsr/2006-anl/pgbhnl.gdas.2006081500.grb2.trop.nc"/)
diri = "/arctic1/nick/cases/cfsr/2006-anl/"
;files = systemfunc("ls -1 " + diri + "pgbhnl.gdas.*.nc")
;files = systemfunc("ls -1 " + diri + "*trop.nc")
files = systemfunc("ls -1 " + diri + "pgbhnl.gdas.2006090100.grb2.trop.nc")

nfiles = dimsizes(files)
set_default_fillvalue("double",  -99999)

stride = 1

;; pull in some static data from the first output file

ff = addfile(files(0),"r")
lon1d = ff->lon_0(0::stride)
lat1d = ff->lat_0(0::stride)
delete(ff)

nlon = dimsizes(lon1d)
nlat = dimsizes(lat1d)
print(" Sizes are "+nlon+" "+nlat)

;meshgrid the lat/lons
lonArray = new(nlon*nlat,float)
latArray = new(nlon*nlat,float)
do iLat=0,nlat-1
do iLon=0,nlon-1
lonArray(iLat*nlon+iLon) = lon1d(iLon)
latArray(iLat*nlon+iLon) = lat1d(iLat)
end do
end do

;make it for all points a la meshgrid
;latAll = conform_dims((/nlat*nlon/),lat1d,0)
;lonAll = conform_dims((/nlat*nlon/),lon1d,0)

res@sfXArray = lonArray
res@sfYArray = latArray
sres@sfXArray = lonArray
sres@sfYArray = latArray

field = new(nlat*nlon,float)

;; loop over times

do nf = 0, nfiles(0)-1
;do nf = nfiles(0)-1, nfiles(0)-1

print("processing file "+files(nf))

ff = addfile(files(nf),"r")

;p = ff->PRES_P0_L7_GLL0(:,:)
;T = ff->TMP_P0_L7_GLL0(:,:)
;The 2pvu fields have interpolation to -1.99999999e-06 and 1.99999999e-06
;Since we plot over the Arctic, don't need to learn how to deal with
;masking (yet?)
p = ff->PRES_P0_L109_GLL0(1,:,:)
T = ff->TMP_P0_L109_GLL0(1,:,:)

;theta = T*(p/1.e5)**(-287.0/1004.0) ; ncl has no power??!??
theta = T*exp((-287.0/1004.0)*log(p/1.e5))
;dims = dimsizes(fieldin)

;print some variable info
printVarSummary(theta)
printMinMax(theta, True)

;nTime = dims(0)
nTime = 1

do iTime=0,nTime-1 ; no time variable in cfsr files

do iLat=0,nlat-1
do iLon=0,nlon-1
field(iLat*nlon+iLon) = theta(iLat,iLon)
end do
end do

field@long_name = "theta_2pvu "+ p@initial_time + " f" + p@forecast_time
field@_FillValue = -99999
;print(" plotting time "+ff->xtime(iTime,0:18))
field@units = ""

res@lbTitleString="K"
res@lbTitlePosition= "Bottom"
res@lbTitleFontHeightF= 0.018
res@lbTitleDirection= "Across"

sres@cnInfoLabelOn = False ;Turn off 'contour from x to y by z' that gets cut off

plot = gsn_csm_contour_map(wks,field,res)
draw(plot)
frame(wks)

end do

delete(ff)
;delete(fieldin)
;delete(fieldin1)

end do

end

;=============================================================================================
