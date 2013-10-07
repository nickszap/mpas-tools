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
sres@cnLineLabelsOn = False
;sres@cnLineThicknessF = 2.6
sres@cnLineLabelFontHeightF = 0.009 
sres@cnInfoLabelFontHeightF = 0.01          
sres@cnLineLabelPlacementMode = "Computed"
;sres@cnLineLabelDensityF = 0.80


;; 500 hPa plots
;; theta plot params at level 20
res@cnLevelSelectionMode = "ManualLevels"
res@cnMinLevelValF  =  300.
res@cnMaxLevelValF  =  380.
res@cnLevelSpacingF =  5.

;; height plot params
sres@cnLevelSelectionMode = "ManualLevels"
sres@cnMinLevelValF  =  4400.
sres@cnMaxLevelValF  =  6200.
sres@cnLevelSpacingF =  50.

wks = gsn_open_wks("pdf","thetal20_h500Pa")
gsn_define_colormap(wks,"BlueYellowRed")

files = systemfunc("find /arctic1/nick/cases/v1.0/x4/2week/kf/x4.kf.output*.nc -type f -print")

;files = (/ \
;"/arctic1/nick/cases/v1.0/x4/2week/kf/x4.kf.output.2006-07-24_12.00.00.nc"/)

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

fieldin1 = ff->height_500hPa(:,:)
fieldin = ff->theta(:,:,20)
dims = dimsizes(fieldin)

nTime = dims(0)
;nTime = 4

do iTime=0,nTime-1

field = fieldin(iTime,0::stride)
field1 = fieldin1(iTime,0::stride)
dimsf = dimsizes(field)
hdiff = abs(field1-terrain_hgt)

do i=0,dimsf(0)-1
  if(hdiff(i) .lt. 1.) then
    field(i) = -99999
    field1(i) = -99999
  end if
end do

field@long_name = "theta_l20, z 500hPa "+ff->xtime(iTime,0:18)
field@_FillValue = -99999
field1@_FillValue = -99999
print(" plotting time "+ff->xtime(iTime,0:18))
field@units = "[K]"
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
