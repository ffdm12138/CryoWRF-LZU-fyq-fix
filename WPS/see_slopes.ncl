; This script may be used to calculate the maximum slopes in x- and
;  y-direction in WRF input topography.

; - Add a file "a"
; - specify the domain "domain" and corresponding grid distances -> WRF
;   domains are usually called "d0*". Current grid distances are set to: d01->
;   1350m, d02->450m, d03->150m and d04-> 50m
; - Chose slope direction "slodir" -> "x" or "y"
; - Chose the angle above which slopes should be plotted "plot" -> "all" or
;   45 to plot slopes above 45°
; - Give the number of "smoothing" cycles (will be used to name the output files)

; (c) Franziska Gerber, 2018


load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/wrf/WRFUserARW.ncl" ; Add special wrf functions


begin

inputfolder = "/home/lig/RUN/Tibet_2019-2020/wpsout/";"smooth1/"
domain = "d02"
slodir = "x"
plot = "all" ;"all" or degrees above which it should be plotted
smoothing = 1
method = "_121"  ; "": smooth-desmooth_special, "_sds": smooth-desmooth, "_121": 1-2-1



; Read file
a = addfile(inputfolder+"geo_em."+domain+".nc","r")

; Create figure
type = "png"					  ; Chose output file type
wks = gsn_open_wks(type,"slopes"+plot+"_"+slodir+"-"+domain+"_smooth"+smoothing+method) ; Create a plot workstation
opts = True 					  ; Set some Basic Plot options
opts@MainTitle = "GEOGRID FIELDS"

; Set grid distances for the respective domains
if (domain .eq. "d01")
  griddist = 27000
end if
if (domain .eq. "d02")
  griddist = 9000
end if
;if (domain .eq. "d03")
  griddist = 3000
;end if
if(domain .eq. "d04")
  griddist = 1000
end if



ter = a->HGT_M(0,:,:)				  ; Read the terrain height
printVarSummary(ter)		  		  ; Print a summary of the terrain height

dims = dimsizes(ter)				  ; Get dimensions of ter

slox = new((/dims(0)-1,dims(1)/),float)		  ; Calculate slopes
sloy = new((/dims(0),dims(1)-1/),float)
terx = new((/dims(0)-1,dims(1)/),float)
tery = new((/dims(0),dims(1)-1/),float)


; Slopes in y-direction
do i = 0,dims(0)-1
  do j = 1,dims(1)-1
     sloy(i,j-1) = atan((ter(i,j)-ter(i,j-1))/griddist)*360/(2*3.14)
     tery(i,j-1) = ter(i,j)-(ter(i,j)-ter(i,j-1))/2
  end do
end do
print("Maximum slope y-direction: "+max(abs(sloy)))


; slope in x-direction
do j = 0,dims(1)-1
  do i = 1,dims(0)-1
     slox(i-1,j) = atan((ter(i,j)-ter(i-1,j))/griddist)*360/(2*3.14)
     terx(i-1,j) = ter(i,j)-(ter(i,j)-ter(i-1,j))/2
  end do
end do
print("Maximum slope x-direction: "+max(abs(slox)))





if (plot .eq. "all")				  ; Create figure with all slopes
  res = opts ; Use basic options for this field
  res@cnFillOn = True ; Create a color fill plot
  res@gsnDraw = False
  res@gsnFrame = False

  rescont = res
  rescont@cnFillOn = False
  rescont@cnLineLabelsOn = False
  rescont@cnHighLabelsOn = False
  rescont@cnLowLabelsOn = False
  ;rescont@ContourParameters = (/ 700., 3500., 200. /)
  rescont@gsnContourLineThicknessesScale = 0.5


  if (slodir .eq. "x")
    contour = wrf_contour(a,wks,abs(slox),res)
    contour_lines = wrf_contour(a,wks,terx,rescont)
  end if
  if (slodir .eq. "y")
    contour = wrf_contour(a,wks,abs(sloy),res)
    contour_lines = wrf_contour(a,wks,tery,rescont)
  end if
  pltres = True ; Set plot options
  mpres = True ; Set map options
  plotting = wrf_map_overlays(a,wks,(/contour,contour_lines/),pltres,mpres) ; Plot the data over a map background
end if







if (plot .eq. 45)				  ; Plot slopes > 45°
  res = opts ; Use basic options for this field
  res@cnFillOn = True ; Create a color fill plot
  res@cnMissingValPerimOn = True
  res@cnMissingValPerimColor = "black"
  res@gsnDraw = False
  res@gsnFrame = False

  rescont = res
  rescont@cnFillOn = False
  rescont@cnLineLabelsOn = False
  rescont@cnHighLabelsOn = False
  rescont@cnLowLabelsOn = False
  ;rescont@ContourParameters = (/ 700., 3500., 200. /)
  rescont@gsnContourLineThicknessesScale = 0.5

  if (slodir .eq. "x")
    slox = mask(slox,slox .lt. plot, False)
    contour = wrf_contour(a,wks,abs(slox),res)
    contour_lines = wrf_contour(a,wks,terx,rescont)
  end if
  if (slodir .eq. "y")
    sloy = mask(sloy,sloy .lt. plot, False)
    contour = wrf_contour(a,wks,abs(sloy),res)
    contour_lines = wrf_contour(a,wks,tery,rescont)
  end if

  pltres = True ; Set plot options
  mpres = True ; Set map options
  plotting = wrf_map_overlays(a,wks,(/contour, contour_lines/),pltres,mpres) ; Plot the data over a map background
end if



end
