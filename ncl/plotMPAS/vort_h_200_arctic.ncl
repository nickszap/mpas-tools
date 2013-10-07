;=============================================================================================
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
;=============================================================================================

external mpas_vort_cell "./mpas_vort_cell.so"
external mpas_filter_cells "./mpas_filter_cells.so"

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
sres@cnLineLabelsOn = False
;sres@cnLineThicknessF = 2.6
sres@cnLineLabelFontHeightF = 0.009 
sres@cnInfoLabelFontHeightF = 0.01          
sres@cnLineLabelPlacementMode = "Computed"
;sres@cnLineLabelDensityF = 0.80


res@cnLevelSelectionMode = "ManualLevels"
res@cnMinLevelValF  =  -0.00015
res@cnMaxLevelValF  =  0.00015
res@cnLevelSpacingF =  0.00001

sres@cnLevelSelectionMode = "ManualLevels"
sres@cnMinLevelValF  =  10500.
sres@cnMaxLevelValF  =  12600.
sres@cnLevelSpacingF =  100.

wks = gsn_open_wks("pdf","vort_h_200_arctic")
gsn_define_colormap(wks,"BlueYellowRed")

files = ("/arctic1/nick/cases/v1.0/x4/august/x4.163842.tiedtke.output.2006-08-01_00.00.00.nc")

;files = (/ \
;"/arctic1/nick/cases/v1.0/x4/august/x4.163842.tiedtke.output.2006-08-01_00.00.00.nc",\
;"/arctic1/nick/cases/v1.0/x4/august/x4.163842.tiedtke.output.2006-08-15_00.00.00.nc",\
;"/arctic1/nick/cases/v1.0/x4/august/x4.163842.tiedtke.output.2006-08-29_00.00.00.nc" /)


nfiles = dimsizes(files)

set_default_fillvalue("double",  -99999)

stride = 1

ff = addfile(files(0),"r")

r2d = 57.2957795 ; conversion from radians to degrees.
res@sfXArray = ff->lonCell(0::stride) * r2d
res@sfYArray = ff->latCell(0::stride) * r2d
sres@sfXArray = ff->lonCell(0::stride) * r2d
sres@sfYArray = ff->latCell(0::stride) * r2d
terrain_hgt_in = ff->zgrid(:,0)
terrain_hgt = terrain_hgt_in(0::stride)

verticesOnCell = ff->verticesOnCell
cellsOnCell = ff->cellsOnCell
nEdgesOnCell = ff->nEdgesOnCell
dv = dimsizes(verticesOnCell)
maxEdges = dv(1)

delete(ff)

;; loop over times

do nf = 0, nfiles(0)-1
;do nf = nfiles(0)-2, nfiles(0)-1

print("processing file "+files(nf))

ff = addfile(files(nf),"r")
fieldin1 = ff->height_200hPa(:,:)
fieldin_tmp = ff->vorticity_200hPa(:,:)

dims = dimsizes(fieldin_tmp)
nv = dims(1)
dims = dimsizes(fieldin1)
nc = dims(1)
fieldin = fieldin1

;nTime = dims(0)
nTime=5
print("number of times in file "+nTime)

do iTime=0,nTime-1

mpas_vort_cell::mpas_vort_cell1(nEdgesOnCell, verticesOnCell, maxEdges, nc, nv, fieldin_tmp(iTime,:), fieldin(iTime,:) )
nfilt=2
mpas_filter_cells::mpas_filter_cells(nEdgesOnCell, cellsOnCell, maxEdges, nc, fieldin(iTime,:), nfilt )

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

field@long_name = "200 hPa hgt, rel vort "+ff->xtime(iTime,0:18)
field@_FillValue = -99999
field1@_FillValue = -99999
print(" plotting time "+ff->xtime(iTime,0:18))
field@units = "[C]"
plot = gsn_csm_contour_map(wks,field,res)
plot_ov = gsn_csm_contour(wks,field1,sres)
overlay(plot,plot_ov)
draw(plot)
frame(wks)

end do

delete(ff)
delete(fieldin)
delete(fieldin1)
delete(fieldin_tmp)

end do

end

;=============================================================================================
