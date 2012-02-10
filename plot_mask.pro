pro plot_mask, ips

DX = GETENV('DX')
DY = GETENV('DY')
FILE = GETENV('FILE')

mask = dblarr(DX,DY)

openr, lun, FILE, /get_lun
readf, lun, mask

mask_reversed=alog(REVERSE(ABS(mask))+1e-3)

if (ips) then begin
    psfile = "field.ps"
    set_plot, 'ps'
    device, filename=psfile, bits_per_pixel=8,/color $
                                ;  ,xsize=18,ysize=18,yoffset=10
            ,xsize=20,ysize=20,yoffset=10

endif else begin
    set_plot, 'x'
    window, xsize=500,ysize=500, retain=2
    
    
 endelse

px=!x.window*!d.x_vsize         ; position of window in device pixels
py=!y.window*!d.y_vsize

plot,[0],/nodata,/noerase,/xstyle,/ystyle,xrange=[Dx,0],yrange=[0,Dy]
TVSCL, mask_reversed,  !x.window(0),!y.window(0),xsize=!x.window(1)-!x.window(0),ysize=!y.window(1)-!y.window(0),/norm

if (ips) then begin
    device, /close
    set_plot, 'x'
endif

end
