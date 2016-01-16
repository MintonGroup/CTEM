pro bar
device, decomposed=1, retain=2
set_plot, 'PS'
DEVICE,/ENCAPSULATED,filename='bar.eps',/INCHES,xsize=2.0,ysize=3.0,font_size=18
thisDevice = !D.NAME
loadct, 33
TVLCT, 255, 255, 255, 254 ; White color
   !P.Color = '000000'xL
   !P.Background = 'FFFFFF'xL

	;Window, 0, xsize=150, ysize=1000
cgColorbar,/Vertical, COLOR='black',MAXRANGE=100,MINRANGE=-400, font=0 ,position=[0.60,0.1,0.95,0.95]
;dwin = tvrd(true=1)
;write_jpeg, 'bar.jpg', dwin, true=1, quality=90
end
