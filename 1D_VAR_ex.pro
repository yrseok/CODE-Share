; YRSEOK, 241125
;Levenberg-Marquardt Method

pname	= 'LM_method'
;********** Do not use FLTARR, Use DBLARR!!!!!!!! **********
pic_path = '/home/yrseok/ms/rm2/PIC/'
file_path= '/home/yrseok/ms/rm2/'

ch      = 8
level   = 101           ;layer #, # of us76standardout
data    = DBLARR(4,level)       ;us76standardout
layer   = DBLARR(level) ;altitude = layer
pres    = DBLARR(level) ;Pressure
temp    = DBLARR(level) ;temperature    ; for HW2-2
zeta    = DBLARR(level) ;vertical coordinate, altitude
pn      = DBLARR(ch)
wf      = DBLARR(level, ch)     ;weighting function
transmit= DBLARR(ch, level) ;transmittance

;====================================
;            READ DATA
;====================================

   file = file_path+'us76standardout.dat'  ;layer, pressure, temp, density?
   OPENR, 1, file
   READF, 1, data
   CLOSE, 1

   layer = data(0,*)
   pres = data(1,*)*0.001  & p0 = pres(0)  ;Surface Pressure

;=== < HW2-1 Weighting Function >===
;=========< Calculation >=========
z1 = 2.
dlnpn = 0.75 ;Spacing,  ;polynomial idx ; known
zeta = -ALOG(pres(*)/p0)                ; -ln(p/p0)

pn = exp(-(z1+dlnpn*(FINDGEN(ch))))*p0  ; pn is pressure at maximum wf at zeta
                                        ; cal using zeta equation

FOR i=0, ch-1 DO BEGIN                  ; WF=dTransmittance/dz
        transmit(i,*) = exp(-pres(*)/pn(i))
        wf(*,i) = pres(*)/pn(i)*transmit(i,*)
ENDFOR
;====================================
;               PLOT
;====================================

        ct=['dark violet','blue','light sky blue',$
        'aqua','lime','yellow','dark orange','red','black']

;       win = WINDOW(DIMENSIONS=[600,600])
;       pp = PLOT(wf(*,0),zeta(*),COLOR=ct[0],/CURRENT,$
;               XTITLE='Weighting Function', YTITLE='Altitude, -ln(p/p0)',$
;               XRANGE=[0,0.4], YRANGE=[0,10], FONT_SIZE=13, THICK=1.5)

;       FOR j=1, ch-1 DO BEGIN
;               pl = PLOT(wf(*,j),zeta(*), COLOR=ct[j],/OVERPLOT,THICK=1.5)
;       ENDFOR

;       WIN.SAVE, +pname+'.png'
;       SPAWN, 'mv -f ' + pname + '.png ' + pic_path
;       WIN.CLOSE

;======< HW2-2 Simulated radiances >======
temp            = data(2,*)
pl_temp         = DBLARR(level)         ; Retrieval Temp

Wj              = DBLARR(ch, level)     ; each channel
Gi              = DBLARR(ch, level)     ; each channel
Cij             = DBLARR(ch, ch)        ; Cij, square matix
Li              = DBLARR(ch)
;Li = integral 0~infinte, B(temp)*weighting function dz, i=1,2,...,m -> channel
dzeta           = DBLARR(level-1) & FOR i=0, level-2 DO dzeta[i]=zeta[i+1]-zeta[i]

;Add Noise
;Noise has then been added to the simulated measurements,
;with an r.m.s. valuE of 0.5 K, and the retrieval repeated.
err_Li          = DBLARR(ch)
err_temp        = DBLARR(level)
rms             = 0.5
;rms            = 0.05
err             = IMSL_RANDOM(ch) + rms

PRINT, "error: ", err
;========= < Step1. Cal radiance > ============
FOR ii=0, ch-1 DO BEGIN
        Li[ii] = TOTAL((temp(*)*wf(*,ii))*dzeta(*))
        err_Li[ii] = Li[ii] + err[ii]           ;Noise

;Li            T(zeta)      ki(z)      dz
ENDFOR

PRINT, "Calculated radiance -----> ", Li

;========= < Step2. Retrieve Temp & Noise > ============

        FOR i=0, ch-1 DO BEGIN  ; Wj = zeta^j-1
                FOR z=0, level-1 DO Wj(i,z) = zeta(z)^(i)
        ;       FOR z=0, level-1 DO Wj(i,z) = sin(2*!PI*(i+1)*zeta(z)/zeta(-1))
        ENDFOR

        FOR i=0, ch-1 DO BEGIN  ; Cij
        FOR j=0, ch-1 DO BEGIN
                Cij[i,j] = TOTAL(Wj(j,*)*(wf(*,i)*dzeta(*)))
                ;Cij       integral Wj(z)*Ki(z)*dz
        ENDFOR
        ENDFOR

                Cji     = INVERT(Cij)           ;Inverse matrix


        FOR i=0, ch-1 DO BEGIN          ;Contribution Function
        FOR j=0, level-1 DO BEGIN
        Gi[i,j] = TOTAL(Wj(*,j)*Cji(*,i))       ;contribution function
        ;Gi     = sigma(Wj*Cji)
        ;contribution function depends on C & W
        ;W is selected
        ;C --> depending on WF and dzeta
        ENDFOR
        ENDFOR


        FOR k=0, level-1 DO BEGIN
                pl_temp[k] = TOTAL(Gi[*,k]*Li(*))
                err_temp[k] = TOTAL(Gi[*,k]*err_Li(*))
                ;BT        = sigma(Gi*Li)
        ENDFOR
;====================================
;               PLOT
;====================================

;ct_tb = ['black', 'deep pink', 'orange']

;win = WINDOW(DIMENSIONS=[600,600])
;pl = PLOT(temp, zeta, COLOR=ct_tb[0],/CURRENT, $
;         XRANGE=[150,300], YRANGE=[0,10], THICK=2., FONT_SIZE=13,$
;        XTITLE='Retrieved quantity', YTITLE='Altitude, -ln(p/p0)',$
;        NAME='US profile temp')

;pl2 = PLOT(pl_temp, zeta, COLOR=ct_tb[1],$
;        /OVERPLOT, THICK=3., NAME='Retrieved temp', LINESTYLE='dot')
;pl3 = PLOT(err_temp, zeta, COLOR=ct_tb[2], /OVERPLOT, THICK=3., $
;        NAME='Retrieved temp with 0.5 K RMSE', LINESTYLE='dash')

;name = ['US profile temp', 'Retrieved temp', 'Retrieved temp with 0.5 K RMSE']
;txt  = TEXT([0.24,0.44,0.735], [0.91,0.91,0.91], name, FONT_SIZE=13, /NORMAL,$
;        ALIGNMENT=0.5, VERTICAL_ALIGNMENT=0.5, COLOR=ct_tb)

;        WIN.SAVE, +pname2+'.png'
;        SPAWN, 'mv -f ' + pname2 + '.png ' + pic_path
;        WIN.CLOSE

;========== < HW2-3 Contribution function > ==========
; Reproduce contribution function

;win = WINDOW(DIMENSIONS=[600,600])

;pp = PLOT(Gi[0,*],zeta,/NODATA,/CURRENT, YTITLE='Altitude, -ln(p/p0)',$
;        XTITLE='Contribution function',FONT_SIZE=13, $
;        XRANGE=[-20,20], YRANGE=[0,10])

;FOR zz=0, ch-1 DO BEGIN
;        gg = PLOT(Gi[zz,*],zeta,COLOR=ct[zz],/OVERPLOT,THICK=2)
;        hh = PLOT(GI[zz,*]*0.02,zeta,COLOR=ct[zz],/OVERPLOT,THICK=2,LINESTYLE='--')
;ENDFOR


;WIN.SAVE, +pname3+'.png'
;SPAWN, 'mv -f ' + pname3 + '.png ' + pic_path
;WIN.CLOSE


;========= < HW2-4 Error amplification factor > ==========

;EAF = DBLARR(level)
;FOR i=0, level-1 DO EAF[i] = SQRT(TOTAL((Gi[*,i])^2D))  ;;    SQRT[Gi(zeta)^2]

;win = WINDOW(DIMENSIONS=[600,600])

;plt = PLOT(EAF(*),zeta,/CURRENT,/XLOG,YTITLE='Altitude, -ln(p/p0)',THICK=2,$
;        XTITLE='Error factor',FONT_SIZE=13, XRANGE=[1,10000],YRANGE=[0,10])



;WIN.SAVE, +pname4+'.png'
;SPAWN, 'mv -f ' + pname4 + '.png ' + pic_path
;WIN.CLOSE

;=== < HW3, SVD > ===
wtfn = transpose(wf)
SVDC, wtfn, S, Vt, U
; 8, 8x101, 8x8 (singular value)
; weighting, Singular value, left, right

;loadct, 39
;  black = 0
;  white = 255
;  ctop  = 254
;  cbot  = 0
;  ncolors = (ctop-cbot)+1

;  ctt = [31, 51, 91, 103, 131, 189, 211, 237, 0]

; WINDOW, 0, XSIZE=1400, YSIZE=1000, RETAIN=2
; WINDOW, 0, XSIZE=900, YSIZE=3600, RETAIN=2
; DEVICE, decomposed=0
;  !P.BACKGROUND=white
;  !P.MULTI = [0, 4, 2, 0, 0]  ; 4¿­ 2ÇàÀ¸·Î À§¿¡¼­ ¾Æ·¡ÂÊÀ¸·Î ±×¸²
;  !P.MULTI = [0, 3, 3, 0, 0]  ; 4¿­ 2ÇàÀ¸·Î À§¿¡¼­ ¾Æ·¡ÂÊÀ¸·Î ±×¸²

;  !P.charsize=3.2
;  !P.thick=3.
;  !P.charthick=1.8

;  FOR i=0,ch-1 DO BEGIN
;          plot, Vt(i,*), zeta,/nodata, TITLE=string(S(i),format='(f7.4)'),$
;                   xrange=[-0.3,0.3],yrange=[0.,10.],xthick=2.0,ythick=2.0,$
;          xtitle='singular vector', ytitle='-ln(p/p0)',color=black

;          oplot,Vt(i,*),zeta,color=ctt[i]

;  ENDFOR

;    WRITE_PNG, pic_path+pname5+'_ch9.png',  TVRD(/TRUE)
;    WRITE_PNG, pic_path+pname5+'.png',  TVRD(/TRUE)


;================= < HW5. Levenberg-Marquardt method > =====================
Sa	= DBLARR(level, level) 
Se	= DBLARR(ch, ch)	
; a priori covariance
 Var_a_prior	= DBLARR(level)
 Var_a_prior(*) = 100.
 Regress_coe = 0.95
 p_cord = FINDGEN(101)*0.1
 Del_z = ABS(p_cord(1)-p_cord(0))            ; in unit of log(p)
;
 Length_scale = -Del_z/(2.*ALOG(Regress_coe))
 Vert_scale = Del_z/Length_scale
;
 For i=0, level-1 DO BEGIN
     For j=0, level-1 DO BEGIN
         index = ABS(i-j)
;        Cov_a_prior(i,j) = Var_a_prior(i)*(Regress_coe)^(2.*index)
         Sa(i,j) = Var_a_prior(i)*EXP(-index*Vert_scale)
     ENDFOR
 ENDFOR

;; Error Covariance Matrix
;
 Se = FLTARR(ch, ch)
 FOR i=0, ch-1 DO BEGIN
     FOR j=0, ch-1 DO BEGIN
         Se(i,j) = 0.0         ; K(2)
         IF(i EQ j) THEN BEGIN
            Se(i,j) = 0.25              ; K(2)
         ENDIF
     ENDFOR
 ENDFOR

K_matrix = DBLARR(ch, level)
Fxi	 = DBLARR(ch)
Sa_inv	= INVERT(Sa)
Se_inv	= INVERT(Se)
gmma	= 0.1
xa	= temp
xi	= pl_temp
yy	= Li
tol = 1e-6                     ; Convergence threshold

max_iter= 100
dzeta	= [0] + dzeta
FOR iter=0, max_iter-1 DO BEGIN
	For ix =0, ch-1 do begin
		Fxi[ix] = total((xi(*)*wf[*,ix])*dzeta(*))
		k_matrix[ix,*] = wf(*,ix)*dzeta(*)
	endfor
	kt_matrix= transpose(k_matrix)
	residual = yy-Fxi


    AA = INVERT((1 + gmma) * Sa_inv + KT_matrix # Se_inv # K_matrix)
    BB = (KT_matrix # Se_inv # residual) - (Sa_inv # (xi - xa))
    delta_xi = AA # BB

    xi = xi + delta_xi

    IF max(ABS(delta_xi)) LT tol THEN BEGIN
	    PRINT, 'converged at iteration: ', iter
	    BREAK
    ENDIF
ENDFOR
PRINT, 'Final retrieved temperature: ', xi

ct_tb = ['black', 'deep pink', 'orange']

win = WINDOW(DIMENSIONS=[600,600])
pl = PLOT(temp, zeta, COLOR=ct_tb[0],/CURRENT, $
         XRANGE=[150,300], YRANGE=[0,10], THICK=2., FONT_SIZE=13,$
        XTITLE='Retrieved quantity', YTITLE='Altitude, -ln(p/p0)',$
        NAME='US profile temp')

pl2 = PLOT(pl_temp, zeta, COLOR=ct_tb[1],$
        /OVERPLOT, THICK=3., NAME='Retrieved temp', LINESTYLE='dot')
pl3 = PLOT(xi, zeta, color=ct_tb[2], $
	/overplot, thick=3, name='Levenberg-Marquardt method')

name = ['US profile temp', 'Retrieved temp', 'Levenberg-Marquardt method']
txt  = TEXT([0.24,0.44,0.735], [0.91,0.91,0.91], name, FONT_SIZE=13, /NORMAL,$
        ALIGNMENT=0.5, VERTICAL_ALIGNMENT=0.5, COLOR=ct_tb)


 WIN.SAVE, +pname+'.png'
 SPAWN, 'mv -f ' + pname + '.png ' + pic_path
 WIN.CLOSE

 win = window(dimensions=[600,600])
 pl3= PLOT(xi, zeta, color=ct_tb[2], XRANGE=[600,1300], YRANGE=[0,10], $
	 THICK=2., FONT_SIZE=13, XTITLE='Retrieved quantity',$
	 YTITLE='Altitude, -ln(p/p0)',$
        /current, title='Levenberg-Marquardt method')
win.save, +pname+'_LM.png'
spawn, 'mv -f '+pname+'_LM.png'+pic_path
win.close

END
