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
sres@cnLevelSelectionMode = "ManualLevels"
sres@cnMinLevelValF  =  4500.
sres@cnMaxLevelValF  =  6500.
sres@cnLevelSpacingF =  50.

;wks = gsn_open_wks("png","x1_t_2007080100_6hrs_thgt_500")
wks = gsn_open_wks("png", "x1_cfsr_2007080100_1day_thgt_500")
gsn_define_colormap(wks,"BlueYellowRed")

;diri = "/arctic1/nick/cases/v1.0/"
;diri = "/arctic1/nick/cases/v1.0/x4/2week/tiedtke/"
;diri = "/arctic1/nick/cases/v1.0/x4/august/tiedtke/"
;diri = "/arctic1/nick/cases/v1.0/x4/august/kf/"
;diri = "/arctic1/nick/cases/v1.0/2007/"
diri = "/arctic1/nick/cases/cfsr/2007/"
;files = systemfunc("ls -1 " + diri + "x1.tiedtke.output.*")
files = systemfunc("ls -1 " + diri + "x1.cfsr.output.*")

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
fieldin = ff->temperature_500hPa(:,:)
dims = dimsizes(fieldin)
nTime = dims(0)

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

;field@long_name = "500 hPa T and height "+ff->xtime(iTime,0:18)
field@_FillValue = -99999
field1@_FillValue = -99999
print(" plotting time "+ff->xtime(iTime,0:18))
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
