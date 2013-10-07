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
sres@cnLineLabelsOn = True
sres@cnLineThicknessF = 2.0
sres@cnLineLabelFontHeightF = 0.009 
sres@cnInfoLabelFontHeightF = 0.01          
sres@cnLineLabelPlacementMode = "Computed"
sres@cnLineLabelDensityF = 1.2


;; 500 hPa plots
;; temperature plot params
res@cnLevelSelectionMode = "ManualLevels"
res@cnMinLevelValF  =  225.
res@cnMaxLevelValF  =  275.
res@cnLevelSpacingF =  5.

;; height plot params
;sres@cnLevelSelectionMode = "ManualLevels"
sres@cnMinLevelValF  =  4500.
sres@cnMaxLevelValF  =  6500.
sres@cnLevelSpacingF =  50.

wks = gsn_open_wks("png","test1")
gsn_define_colormap(wks,"BlueYellowRed")

files = (/ \
"/home/nickszap/research/cases/cfsr/nomads/pgbhnl.gdas.2006081500.grb2.mb.nc"/)

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
field1 = new(nlat*nlon,float)

;; loop over times

do nf = 0, nfiles(0)-1
;do nf = nfiles(0)-1, nfiles(0)-1

print("processing file "+files(nf))

ff = addfile(files(nf),"r")

pInd = 2 ; look in lv_ISBL0 of grb2-to-nc file
fieldin1 = ff->HGT_P0_L100_GLL0(pInd,:,:)
fieldin = ff->TMP_P0_L100_GLL0(pInd,:,:)
dims = dimsizes(fieldin)

;unroll 2d fields into flat arrays
field = new(nlat*nlon,float)
field1 = new(nlat*nlon,float)

;nTime = dims(0)
nTime = 1

do iTime=0,nTime-1 ; no time variable in cfsr files

do iLat=0,nlat-1
do iLon=0,nlon-1
field(iLat*nlon+iLon) = fieldin(iLat,iLon)
field1(iLat*nlon+iLon) = fieldin1(iLat,iLon)
end do
end do

field@long_name = "."
field@_FillValue = -99999
field1@_FillValue = -99999
;print(" plotting time "+ff->xtime(iTime,0:18))
;field@units = "[C]"

res@lbTitleString="K"
res@lbTitlePosition= "Bottom"
res@lbTitleFontHeightF= 0.018
res@lbTitleDirection= "Across"

sres@cnInfoLabelOn = False ;Turn off 'contour from x to y by z' that gets cut off

plot = gsn_csm_contour_map(wks,field,res)
plot_ov = gsn_csm_contour(wks,field1,sres)
overlay(plot,plot_ov)
draw(plot)
frame(wks)

end do

delete(ff)
delete(fieldin)
delete(fieldin1)

end do

end

;=============================================================================================
