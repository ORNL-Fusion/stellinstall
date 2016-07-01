interactive=0
if(!d.name eq 'x' or !d.name eq 'X') then interactive=1
ans=' '
oldchrsiz=!p.charsize
;!p.charsize=0
;if interactive then !p.charsize=2
!p.multi=0
m=read_rfun_nc("map_g118162.93030.nc")
g=readg(strtrim(string(m.eqdskname),2))
eqdsk=strtrim(string(m.eqdskname),2)
afile=strtrim(string(m.eqdskname),2)
basename=strmid(afile,1,strlen(strtrim(string(m.eqdskname),2))-1)
afile='a'+strmid(afile,1,strlen(afile)-1)
filename='plot'+strtrim(basename,2)+'.ps'
a=reada(afile)
cmg=cmgcolors()
mystyle={xstyle:1,ystyle:1,charsize:2}
mu0=double(4e-7*!pi)
twopi=double(2*!pi)
col=n_elements(m.xs(0,*))
row=n_elements(m.xs(*,0))
gpsi=twopi*(findgen(g.mw)/(g.mw-1)*(g.ssibry-g.ssimag)+g.ssimag); Wb !!!
jtor=m.xs*0.
s=(m.pflx-m.pflx(0))/(m.pflx(row-1)-m.pflx(0))
s1=findgen(n_elements(g.pprime))/(n_elements(g.pprime)-1)
!p.multi=[0,2,1]
plot,s,m.pprime,/nodata,title="p' vs psin",subtitle=eqdsk,_extra=mystyle
oplot,s,m.pprime,psym=1,color=2
oplot,s1,g.pprime*twopi,color=4
plot,s,mu0*m.ffprime,/nodata,title="u0 ff' vs psin",subtitle=eqdsk,_extra=mystyle
oplot,s,mu0*m.ffprime,psym=1,color=2
oplot,s1,g.ffprim*twopi/mu0,color=4
if interactive then read,'$(a1)',ans,prompt='<CR>'
; J terms on midplane
rmid=[reverse(m.xs(*,col/2)),m.xs(1:row-1,0)]
for j=0,row-1 do for k=0,col-1 do jtor(j,k)=$
m.xs(j,k)*m.pprime(j)
jtor=-jtor/twopi
jmid=[reverse(jtor(*,col/2)),jtor(1:row-1,0)]
plot,rmid,jmid,/nodata,title="midplane Rp'",_extra=mystyle
oplot,rmid,jmid,psym=1,color=3
for j=0,row-1 do for k=0,col-1 do jtor(j,k)=$
mu0*m.ffprime(j)/m.xs(j,k)
jtor=-jtor/twopi
jmid=[reverse(jtor(*,col/2)),jtor(1:row-1,0)]
plot,rmid,jmid,/nodata,title="midplane u0 ff'/R'",_extra=mystyle
oplot,rmid,jmid,psym=1,color=3
if interactive then read,'$(a1)',ans,prompt='<CR>'
jtor=0
jmid=0
rmid=0
s1=0
jtor=m.xs*0.
!p.multi=0
for j=0,row-1 do for k=0,col-1 do jtor(j,k)=$
m.xs(j,k)*m.pprime(j)+mu0*m.ffprime(j)/m.xs(j,k)
jtor=-jtor/twopi	;	This is the real deal
rmid=[reverse(m.xs(*,col/2)),m.xs(1:row-1,0)]
jmid=[reverse(jtor(*,col/2)),jtor(1:row-1,0)]
plot,rmid,jmid,/nodata, title="MIDPLANE: J(R)==(R p'+u0 ff' /R) vs R",$
subtitle="Local current density (A/m^3)	"+eqdsk,_extra=mystyle
oplot,rmid,jmid,psym=1,color=3
if interactive then read,'$(a1)',ans,prompt='<CR>'
dl=m.pflx*0.
for j=0,row-1 do dl(j)=max(m.arcsur(j,*)-shift(m.arcsur(j,*),1))
; arcsur are equal arc length
; volume elements
jbar=m.pflx*0.
dlob=jbar
for j=1,row-1 do for k=0,col-1 do $
dlob(j)=dlob(j)+dl(j)/m.bp(j,k)
dlob(0)=dlob(1)-(dlob(2)-dlob(1))
for j=1,row-1 do for k=0,col-1 do $
jbar(j)=jbar(j)+dl(j)*jtor(j,k)/m.bp(j,k)
jbar=jbar/dlob
jbar(0)=jtor(0,0)
jovr=m.pflx*0.
dlob=jovr
for j=1,row-1 do for k=0,col-1 do $
dlob(j)=dlob(j)+dl(j)/m.bp(j,k)
dlob(0)=dlob(1)-(dlob(2)-dlob(1))
for j=1,row-1 do for k=0,col-1 do $
jovr(j)=jovr(j)+dl(j)*jtor(j,k)/m.xs(j,k)/m.bp(j,k)
jovr=jovr/dlob
jovr(0)=jtor(0,0)/m.xs(0,0)
plot,m.tflx,jovr,/nodata,_extra=mystyle,$
title="<J/R>==INTEGRAL {j.dl/R/Bp}/INTEGRAL {dl/Bp} vs tflx",subtitle=eqdsk
oplot,m.tflx,jovr,color=4,thick=2
if interactive then read,'$(a1)',ans,prompt='<CR>'
s=m.tflx/m.tflx(row-1)
dv=(m.vprime +shift(m.vprime,1))/2*(m.pflx-shift(m.pflx,1))
dv(0)=0
v=dv*0.
for j=0,row-1 do v(j)=total(dv(0:j))
subtitle=$
"plasma volume = "+string(v(row-1),'$(f7.4)')+$
"; EFIT  volume = "+string(a.d.volume/1e6,'$(f7.4)') +"  "
plot,m.tflx,v,/nodata,title=" Volume==(dV/dpsi)*dpsi vs tflx",subtitle=$
subtitle+eqdsk,_extra=mystyle
oplot,m.tflx,v,color=4,thick=2
if interactive then read,'$(a1)',ans,prompt='<CR>'
itor=jbar*0.
for j=1,row-1 do itor(j)=itor(j-1)+dv(j)*jovr(j)/twopi
subtitle="plasma current [MA] = "+string(itor(row-1)/1e6,'$(f7.4)')+$
"; EFIT  current [MA] = "+string(a.d.ipmhd/1e6,'$(f7.4)')+"      "
plot,m.tflx,itor,/nodata,$
title="enclosed toroidal plasma current vs tflx",subtitle=$
subtitle+eqdsk,_extra=mystyle
oplot,m.tflx,itor,color=4,thick=2
if interactive then read,'$(a1)',ans,prompt='<CR>'
iprime=deriv(m.tflx,itor)
s=m.tflx/m.tflx(row-1)
iprime=deriv(s,itor)
singular=-1
ac=svdfit(s,iprime,11,/double,singular=singular)
fit=s*0.
for i=0, n_elements(fit)-1 do fit(i)=polyval(ac,s(i))
plot,s,iprime,/nodata,title="dI/ds vs s and it's svd	"+eqdsk,$
subtitle="SVD check: singular = "+string(singular),_extra=mystyle
oplot,s,iprime,color=4,thick=2
for i=0,n_elements(ac)-1 do xyouts,0.8,i*max(iprime)/n_elements(ac)-1,$
string(ac(i),'$(x,e11.4)')+"      "+eqdsk,color=6
print,"EQDSK	=	"+"'"+eqdsk+"'"
thin=5
oplot,s(0:row-1:thin),fit(0:row-1:thin),color=2,psym=4
print,"from a-file:Volume=	",a.d.volume/1e6,format='(a,e11.5)'
print,"Computed Volume=  	",v(row-1),format='(a,e11.5)'
print,'from a-file:Ip=  	',a.d.ipmhd,format='(a,e11.5)'
print,"Ienclosed=  		",itor(row-1),format='(a,e11.5)'

print,"AC array:"
print,"AC = [	$"
for i=0,n_elements(ac)-2 do $
  print,string(ac(i)/ac(0),'$(x,e11.4)')+'	,	$' 
for i=n_elements(ac)-1,n_elements(ac)-1 do $
  print,string(ac(i)/ac(0),'$(x,e11.4)')+']'
print,'!p.charsize=',!p.charsize
!p.charsize=oldchrsiz
end
; IDL> print,max(gpsi),max(m.pflx)
;       2.2288829       2.2287714 
;IDL> print,m.percenflux
;      0.99995000
; gpsi is linear in its index , alpsi=-1
; so dr_mid is near zero at rmin,rmax and rises to 0.0015 at the axis.
; row=512 col=257
;IDL> plot,gpsi,g.ffprim/twopi/mu0^2 is identical to  oplot,m.pflx,m.ffprime
;IDL> plot,gpsi,g.pprime*twopi is identical to  oplot,m.pflx,m.pprime 

