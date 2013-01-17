
      MODULE physical_constants
      USE precision
      implicit none
!----------------------------------------------------------------------
!  Mathematical constants
!----------------------------------------------------------------------
      REAL(rprec), PARAMETER :: pi=3.14159265358979323846264338328_rprec
      REAL(rprec), PARAMETER :: pio2=1.570796326794896619231321691_rprec
      REAL(rprec), PARAMETER :: twopi=6.28318530717958647692528677_rprec
      REAL(rprec), PARAMETER :: sqrt2=1.41421356237309504880168872_rprec
      REAL(rprec), PARAMETER :: degree = twopi / 360
      REAL(rprec), PARAMETER :: one = 1.0_rprec
      REAL(rprec), PARAMETER :: zero = 0.0_rprec
      REAL(rprec), PARAMETER :: two = 2.0_rprec
      REAL(rprec), PARAMETER :: half = 0.5_rprec
      REAL(rprec), PARAMETER :: three = 3.0_rprec
      REAL(rprec), PARAMETER :: four = 4.0_rprec
      REAL(rprec), PARAMETER :: six = 6.0_rprec
      REAL(rprec), PARAMETER :: twofour = 24.0_rprec
!----------------------------------------------------------------------
!  Physical constants
!------------------------------------------------------------------
      REAL(rprec), PARAMETER :: mu0 = 2 * twopi * 1.0e-7_rprec
      END MODULE physical_constants
      subroutine furplm(p,grsq,x,y,v,xp,yp,gsq,np,nx,ny,ntx,nw,nh,is,js)
c ***********************************************************************
c given an x-y grid with a variable p defined on the grid this routine
c finds the level p curves of value v.  that is, it finds np points
c (xp(k),yp(k)) such that p(xp(k),yp(k)) is close to v.
c id is an integer between 1 and 12 or 15 which identifies the option
c used to find the last point.
c if in=1 the curve is moving up and if in=0 the curve is moving down.
c if il=1 the curve is moving left and if il=0 the curve is moving right.
c if the last point was a grid point id=15 and igp is an integer between
c 1 and 12 which identifies the grid point option used to find the last
c point.
c if ig=1 the grid point was entered from the left and if ig=0 the grid
c point was entered from the right.
c iq is an indicator used in that part of the program which searches the
c grid squares surrounding a grid point.
c if ntx is less than zero search for the first point starts from the
c top of the grid and goes down.
c if np eq -1 the search starts in the center of the grid
c if ntx gt 0 the search starts at the bottom of the grid
c if you would like further details see joanne helton.
c ***********************************************************************
cray  lcm (p),(grsq),(x),(y),(xp),(yp),(gsq)
      USE physical_constants
      USE precision
      implicit none
      INTEGER :: js, is, ntx, ny, nx, np, nw, nh, nt, jst, nx2, ny2
      INTEGER :: kprt, j, i, k, iq, ig, in, id, il, jy, ix, igp, l
      INTEGER :: ixx, iyy, jdec
      REAL(rprec) :: v, a, b, c, f, yq, dx, dy , dp, e, d, ff, ws
      REAL(rprec) :: xp(*),yp(*), x1, y1
      REAL(rprec) :: grsq(nw,nh),gsq(*)
      REAL(rprec) :: p(nw,nh),x(nw),y(nh)
      a=zero
      b=zero
      c=zero
      f=zero
      yq=zero
      nt=ntx
      jst=js
      jdec=1
      if(np.ne.-1) go to 1999
      jst=ny/2+1
      jdec=1
      nt=iabs(nt)
 1999 if(nt.gt.0) go to 2004
      jst=ny-1
      jdec=-1
      nt=-nt
 2004 continue
      nx2=nx/2
      ny2=ny/2
      dx=x(2)-x(1)
      dy=y(2)-y(1)
      dp=abs(p(nx2,ny2))
      d=10._dbl**(-1*8)*dp
      e=10._dbl**(-1*4)*dx*three
      ff=10._dbl**(-1*6)*dx
      kprt=0
2001  continue
      j=jst
      i=is
      k=0
      in=0
      iq=0
c ---------------------
c we find a first point
c ---------------------
   13 if((v-p(i,j))*(p(i+1,j)-v).gt.zero) go to 30
      if(abs(p(i,j)-v).ge.d) go to 32
      in=1
      ig=1
      go to 15
   32 i=i+1
      if(i.eq.nx) go to 25
      go to 13
   25 j=j+jdec
      if(jdec.lt.0) go to 2005
      if(j.gt.ny) go to 10
      i=1
      go to 13
 2005 if(j.lt.1) go to 10
      i=1
      go to 13
   30 if(i.gt.(nx-2)) go to 74
      if(i.eq.1) go to 88
      xp(1)=x(i)
      call fit(2,x(i-1),x(i),x(i+1),x(i+2),p(i-1,j),p(i,j),p(i+1,j),
     1         p(i+2,j),xp(1),v,yq)
      ws=xp(1)
      call cubic(x(i-1),x(i),x(i+1),x(i+2),ws,a,b,c,f)
      gsq(1)=a*grsq(i-1,j)+b*grsq(i,j)+c*grsq(i+1,j)+f*grsq(i+2,j)
      go to 75
   88 xp(1)=x(2)
      call fit(2,x(1),x(2),x(3),x(4),p(1,j),p(2,j),p(3,j),p(4,j),
     1         xp(1),v,yq)
      ws=xp(1)
      call cubic(x(1),x(2),x(3),x(4),ws,a,b,c,f)
      gsq(1)=a*grsq(1,j)+b*grsq(2,j)+c*grsq(3,j)+f*grsq(4,j)
      go to 75
   74 xp(1)=x(nx-2)
      call fit(2,x(nx-3),x(nx-2),x(nx-1),x(nx),p(nx-3,j),p(nx-2,j),
     1         p(nx-1,j),p(nx,j),xp(1),v,yq)
      ws=xp(1)
      call cubic(x(nx-3),x(nx-2),x(nx-1),x(nx),ws,a,b,c,f)
      gsq(1)=a*grsq(nx-3,j)+b*grsq(nx-2,j)+c*grsq(nx-1,j)+f*grsq(nx,j)
   75 yp(1)=y(j)
      in=1
      if(xp(1).le.x(i)) xp(1)=x(i)+ff
      if(xp(1).ge.x(i+1)) xp(1)=x(i+1)-ff
      id=20
      il=1
      if(jdec.lt.0) il=0
      k=1
      if(kprt.eq.1) write(66,102)k,xp(k),yp(k),p(i,j),i,j,ix,jy,ig,in,id
     x      ,iq,il
  102 format(i5,3e13.6,10i5)
      go to 5
c ---------------------------------------------------------------
c in the block below we determine where the next point is located
c ---------------------------------------------------------------
c ---------------------------------------------------------------------
c the curve is entering the grid square from the right or from the left
c if il equals 1 the curve is moving left
c ---------------------------------------------------------------------
    4 if(il.eq.1) go to 21
c --------------------------------------------------
c the following options are for a curve moving right
c --------------------------------------------------
      if(i+1.gt.nx) go to 10
c --------------------------
c first we check grid points
c --------------------------
      if(abs(v-p(i+1,j+1)).ge.d) go to 45
      ig=1
      in=1
      i=i+1
      j=j+1
      igp=11
      go to 15
   45 if(abs(v-p(i+1,j)).ge.d) go to 33
      ig=1
      in=0
      i=i+1
      igp=12
      go to 15
c --------------------------------
c now we check between grid points
c --------------------------------
   33 if((v-p(i,j))*(p(i+1,j)-v).le.zero) go to 50
      ix=i
      jy=j
      id=1
      go to 1
   50 if((v-p(i,j+1))*(p(i+1,j+1)-v).le.zero) go to 51
      ix=i
      jy=j+1
      id=2
      go to 1
   51 if((v-p(i+1,j))*(p(i+1,j+1)-v).le.zero) go to 52
      ix=i+1
      jy=j
      id=3
      go to 7
   52 continue
c -------------------------------------------------
c the following options are for a curve moving left
c -------------------------------------------------
      if(i-1.lt.0) go to 10
c --------------------------
c first we check grid points
c --------------------------
   21 if(abs(v-p(i-1,j)).ge.d) go to 46
      in=0
      ig=0
      i=i-1
      igp=1
      go to 15
   46 if(abs(v-p(i-1,j+1)).ge.d) go to 34
      ig=0
      in=1
      i=i-1
      j=j+1
      igp=2
      go to 15
c --------------------------------
c now we check between grid points
c --------------------------------
   34 if((v-p(i,j+1))*(p(i-1,j+1)-v).le. zero) go to 53
      ix=i-1
      jy=j+1
      id=4
      go to 1
   53 if((v-p(i-1,j+1))*(p(i-1,j)-v).le. zero) go to 54
      ix=i-1
      jy=j
      id=5
      go to 7
   54 if((v-p(i,j))*(p(i-1,j)-v).le. zero) go to 55
      ix=i-1
      jy=j
      id=6
      go to 1
   55 continue
c ---------------------------------------------------------------------
c the curve is entering the grid square from the top or from the bottom
c or from a grid point
c if ig equals 1 the grid point was entered from the left
c if in equals 1 the curve is moving up
c ---------------------------------------------------------------------
    5 if(id.eq.15.and.il.eq.1) i=i-1
   72 iq=0
   62 if(in.eq.1) go to 12
c -------------------------------------------------
c the following options are for a curve moving down
c -------------------------------------------------
      if(j-1.lt.1) go to 10
   11 if(abs(v-p(i+1,j-1)).ge.d) go to 42
      ig=1
      in=0
      i=i+1
      j=j-1
      igp=3
      go to 15
   42 if(abs(v-p(i,j-1)).ge.d) go to 17
      ig=0
      in=0
      j=j-1
      igp=5
      go to 15
c --------------------------------
c now we check between grid points
c --------------------------------
   17 if((v-p(i+1,j-1))*(p(i+1,j)-v).le.zero) go to 56
      if(abs(v-p(i+1,j)).le.d) go to 56
      if(id.eq.15.and.igp.eq.7) go to 56
      ix=i+1
      jy=j-1
      id=7
      go to 7
   56 if((v-p(i,j-1))*(p(i+1,j-1)-v).le.zero) go to 57
      ix=i
      jy=j-1
      id=8
      go to 1
   57 if((v-p(i,j))*(p(i,j-1)-v).le.zero) go to 44
      if(abs(v-p(i,j)).le.d) go to 44
      ix=i
      jy=j-1
      id=9
      go to 7
   44 if(iq.ne.1) go to 48
      in=1
      if(il.eq.1) i=i+1
      if(il.eq.0) i=i-1
      go to 72
   48 iq=1
      in=1
      go to 62
c -----------------------------------------------
c the following options are for a curve moving up
c -----------------------------------------------
   12 if(j+1.gt.ny) go to 10
      if(i.eq.nx) go to 10
c --------------------------
c first we check grid points
c --------------------------
      if(abs(v-p(i+1,j+1)).ge.d) go to 40
      ig=1
      in=1
      i=i+1
      j=j+1
      igp=7
      go to 15
   40 if(abs(v-p(i,j+1)).ge.d) go to 61
      ig=0
      in=1
      j=j+1
      igp=9
      go to 15
c --------------------------------
c now we check between grid points
c --------------------------------
   61 if((v-p(i,j+1))*(p(i+1,j+1)-v).le. zero) go to 59
      ix=i
      jy=j+1
      id=10
      go to 1
   59 if((v-p(i+1,j))*(p(i+1,j+1)-v).le. zero) go to 60
      if(abs(v-p(i+1,j)).le.d) go to 60
      ix=i+1
      jy=j
      id=11
      go to 7
   60 if((v-p(i,j+1))*(p(i,j)-v).le. zero) go to 41
      if(abs(v-p(i,j)).le.d) go to 41
      ix=i
      jy=j
      id=12
      go to 7
   41 if(iq.ne.1) go to 49
      in=0
      if(il.eq.1) i=i+1
      if(il.eq.0) i=i-1
      go to 72
   49 iq=1
      in=0
      go to 62
c ---------------------------------------------------------------------
c in the block below the x and y values of the point found are computed
c ---------------------------------------------------------------------
c and stored and control is returned to the appropriate point in the
c block above
    1 k=k+1
      if(k.eq.nt) write(66,106)
  106 format(/,"***nt points found in furplm***",/)
      if(k.eq.nt) go to 2000
      if(ix.eq.1) go to 2
      if(ix.gt.(nx-2)) go to 3
      xp(k)=x(ix)
      call fit(2,x(ix-1),x(ix),x(ix+1),x(ix+2),p(ix-1,jy),p(ix,jy),
     1         p(ix+1,jy),p(ix+2,jy),xp(k),v,yq)
      ws=xp(k)
      call cubic(x(ix-1),x(ix),x(ix+1),x(ix+2),ws,a,b,c,f)
 2009 continue
      gsq(k)=a*grsq(ix-1,jy)+b*grsq(ix,jy)+c*grsq(ix+1,jy)
     1+f*grsq(ix+2,jy)
 2010 continue
      go to 6
    2 xp(k)=x(2)
      call fit(2,x(1),x(2),x(3),x(4),p(1,jy),p(2,jy),p(3,jy),p(4,jy),
     1         xp(k),v,yq)
      ws=xp(k)
      call cubic(x(1),x(2),x(3),x(4),ws,a,b,c,f)
      gsq(k)=a*grsq(1,jy)+b*grsq(2,jy)+c*grsq(3,jy)+f*grsq(4,jy)
      go to 6
    3 xp(k)=x(nx-2)
      call fit(2,x(nx-3),x(nx-2),x(nx-1),x(nx),p(nx-3,jy),p(nx-2,jy),
     1         p(nx-1,jy),p(nx,jy),xp(k),v,yq)
      ws=xp(k)
      call cubic(x(nx-3),x(nx-2),x(nx-1),x(nx),ws,a,b,c,f)
      gsq(k)=a*grsq(nx-3,jy)+b*grsq(nx-2,jy)+c*grsq(nx-1,jy)
     1+f*grsq(nx,jy)
    6 yp(k)=y(jy)
      if((id.eq.4).or.(id.eq.6)) i=i-1
      if(xp(k).le.x(ix)) xp(k)=x(ix)+ff
      if(xp(k).ge.x(ix+1)) xp(k)=x(ix+1)-ff
      if(k.le.3) go to 20
      if((abs(xp(k)-xp(1)).lt.e.and.abs(yp(k)-yp(1)).lt.e)
     1.and.k.ne.1) go to 10
   20 if(jy.eq.1.and.(((id.eq.1).or.(id.eq.6)).or.(id.eq.8))) go to 10
      if(jy.eq.ny.and.(((id.eq.2).or.(id.eq.4)).or.(id.eq.10))) go to 10
      if(((id.eq.2).or.(id.eq.4)).or.(id.eq.10)) j=j+1
      if(id.eq.8) j=j-1
      in=1
      if(yp(k).lt.yp(k-1)) in=0
      il=1
      if(xp(k).gt.xp(k-1)) il=0

      if(k.le.2) go to 5       !       trap against xp(k-2)

      if((abs(xp(k)-xp(k-2)).lt.e.and.abs(yp(k)-yp(k-2)).lt.e).and.k.ge.
     x      3) write(66,107)
  107 format(/,"***point same as point before last***",/)
      if((abs(xp(k)-xp(k-2)).lt.e.and.abs(yp(k)-yp(k-2)).lt.e)
     1.and.k.ge.3) go to 2000
      if(kprt.eq.1) write(66,102)k,xp(k),yp(k),p(i,j),i,j,ix,jy,ig,in,id
     x      ,iq,il
      go to 5
    7 k=k+1
      if(k.eq.nt) write(66,106)
      if(k.eq.nt) go to 2000
      xp(k)=x(ix)
      if(jy.eq.1) go to 8
      if(jy.gt.(ny-2)) go to 9
      yp(k)=y(jy)
      call fit(2,y(jy-1),y(jy),y(jy+1),y(jy+2),p(ix,jy-1),p(ix,jy),
     1         p(ix,jy+1),p(ix,jy+2),yp(k),v,yq)
      ws=yp(k)
      call cubic(y(jy-1),y(jy),y(jy+1),y(jy+2),ws,a,b,c,f)
      gsq(k)=a*grsq(ix,jy-1)+b*grsq(ix,jy)+c*grsq(ix,jy+1)
     1+f*grsq(ix,jy+2)
      go to 14
    8 yp(k)=y(2)
      call fit(2,y(1),y(2),y(3),y(4),p(ix,1),p(ix,2),p(ix,3),p(ix,4),
     1         yp(k),v,yq)
      ws=yp(k)
      call cubic(y(1),y(2),y(3),y(4),ws,a,b,c,f)
      gsq(k)=a*grsq(ix,1)+b*grsq(ix,2)+c*grsq(ix,3)+f*grsq(ix,4)
      go to 14
    9 yp(k)=y(ny-2)
      call fit(2,y(ny-3),y(ny-2),y(ny-1),y(ny),p(ix,ny-3),p(ix,ny-2),
     1         p(ix,ny-1),p(ix,ny),yp(k),v,yq)
      ws=yp(k)
      call cubic(y(ny-3),y(ny-2),y(ny-1),y(ny),ws,a,b,c,f)
      gsq(k)=a*grsq(ix,ny-3)+b*grsq(ix,ny-2)+c*grsq(ix,ny-1)
     1+f*grsq(ix,ny)
   14 if(yp(k).le.y(jy)) yp(k)=y(jy)+ff
      if(yp(k).ge.y(jy+1)) yp(k)=y(jy+1)-ff
      if(k.le.3) go to 19
      if((abs(xp(k)-xp(1)).lt.e.and.abs(yp(k)-yp(1)).lt.e)
     1.and.k.ne.1) go to 10
   19 if(ix.eq.1.and.(((id.eq.5).or.(id.eq.9)).or.(id.eq.12))) go to 10
      if(ix.eq.nx.and.(((id.eq.3).or.(id.eq.7)).or.(id.eq.11))) go to 10
      if(id.eq.5) i=i-1
      if((id.eq.7).or.(id.eq.9)) j=j-1
      in=1
      if(yp(k).lt.yp(k-1)) in=0
      il=1
      if(xp(k).gt.xp(k-1)) il=0
      if(((id.eq.3).or.(id.eq.7.and.il.eq.0)).or.(id.eq.11.and.il.eq.0))
     1i=i+1
      if(kprt.eq.1) write(66,102)k,xp(k),yp(k),p(i,j),i,j,ix,jy,ig,in,id
     x      ,iq,il

      if(k.le.2) go to 4       !       trap against xp(k-2)

      if((abs(xp(k)-xp(k-2)).lt.e.and.abs(yp(k)-yp(k-2)).lt.e).and.k.ge.
     x      3) write(66,107)
      if((abs(xp(k)-xp(k-2)).lt.e.and.abs(yp(k)-yp(k-2)).lt.e)
     1.and.k.ge.3) go to 2000
      go to 4
   15 k=k+1
      if(k.eq.nt) write(66,106)
      if(k.eq.nt) go to 2000
      xp(k)=x(i)
      yp(k)=y(j)
      if((abs(xp(k)-xp(k-2)).lt.e.and.abs(yp(k)-yp(k-2)).lt.e).and.k.ge.
     x      3) write(66,107)
      if((abs(xp(k)-xp(k-2)).lt.e.and.abs(yp(k)-yp(k-2)).lt.e)
     1.and.k.ge.3) go to 2000
      gsq(k)=grsq(i,j)
      id=15
      if(k.le.3) go to 18
      if((abs(xp(k)-xp(1)).lt.e.and.abs(yp(k)-yp(1)).lt.e)
     1.and.k.ne.1) go to 10
   18 if(((i.eq.1.or.i.eq.nx).and.(j.eq.1.or.j.eq.ny)).and.k.ne.1)
     1go to 10
      in=1
      if(yp(k).lt.yp(k-1)) in=0
      il=1
      if(xp(k).gt.xp(k-1)) il=0
      if(kprt.eq.1) write(66,102)k,xp(k),yp(k),p(i,j),i,j,ix,jy,ig,in,id
     x      ,iq,il,igp
      go to 5
   10 continue
      np=k
      return
 2000 continue
c test if looping on separatrix
      do 2020 i=1,k
      x1=xp(i)
      y1=yp(i)
      do 2008 j=i+1,k
      if(abs(xp(j)-x1).gt.e) go to 2008
      if(abs(yp(j)-y1).le.e) go to 2030
 2008 continue
 2020 continue
      go to 2100
 2030 continue
      k=j-i+1
      do  j=1,k
         l=i+j-1
         gsq(j)=gsq(l)
         xp(j)=xp(l)
         yp(j)=yp(l)
      end do
      go to 10
 2100 continue
      write(6,*)'error in furpl--fur1'
      stop
      write(66,2002) v,(ixx,x(ixx),ixx=1,nx),(iyy,y(iyy),iyy=1,ny)
2002  format("  furplm diagnostic center",/,
     .     "  value =",e14.6,/,"   x and y grids follow",/,
     .(i5,e14.6,i5,e14.6,i5,e14.6,i5,e14.6))
      write(66,2003)
2003  format(//,"  k  ","       xp          yp         p(i,j)   ",
     ."  i    j    ix   jy   ig   in   id   iq   il")
      kprt=1
      go to 2001
      end
      subroutine fit(k,x1,x2,x3,x4,y1,y2,y3,y4,x,y,yp)
c --------------------------------------------
c this routine is taken from an nyu program
c set k=1 to find y,yp and k=2 to find x,yp
c --------------------------------------------
cray  lcm (x1),(x2),(x3),(x4),(y1),(y2),(y3),(y4),(x)
      USE physical_constants
      USE precision
      implicit none
      INTEGER :: k, iturn, i
      REAL(rprec) :: x1, x2, x3, x4, y1, y2, y3, y4, x, y, yp
      REAL(rprec) :: c1, c2, c3, c4, d1, d2, d3, d4, crit, xa, ya
      REAL(rprec) :: d12, d13, d14, d23, d24, d34, f, xb, yb, dydx
      iturn=0
      c1=y1/((x1-x2)*(x1-x3)*(x1-x4))
      c2=y2/((x2-x1)*(x2-x3)*(x2-x4))
      c3=y3/((x3-x1)*(x3-x2)*(x3-x4))
      c4=y4/((x4-x1)*(x4-x2)*(x4-x3))
      if(k.eq.2) go to 2
   1  d1=x-x1
      d2=x-x2
      d3=x-x3
      d4=x-x4
      d12=d1*d2
      d13=d1*d3
      d14=d1*d4
      d23=d2*d3
      d24=d2*d4
      d34=d3*d4
      f=(c1*d23+c2*d13+c3*d12)*d4+c4*d12*d3
      yp=c1*(d23+d24+d34)+c2*(d13+d14+d34)
     +  +c3*(d12+d14+d24)+c4*(d12+d13+d23)
      if(k.eq.2) go to 3
      y=f
      return
    2 if(y.ge.min(y1,y2,y3,y4).and.y.le.max(y1,y2,y3,y4)) go to 4
   21 continue
      write(6,11)x1,x2,x3,x4,y1,y2,y3,y4,y
      x=zero
      write(6,500)
  500 format("error in sub fit at label 21")
      write(6,*)'error in furpl--fit1'
      stop
      return
  11  format(" fit",9e14.5)
   4  crit=(abs(y1)+abs(y2)+abs(y3)+abs(y4))*1.d-05
      i=0
      xa=x1
      xb=x2
      ya=y1
      yb=y2
      if((x-xa)*(x-xb).lt.zero) go to 10
      xa=x2
      xb=x3
      ya=y2
      yb=y3
      if((x-xa)*(x-xb).lt.zero) go to 10
      xa=x3
      xb=x4
      ya=y3
      yb=y4
      if((x-xa)*(x-xb).lt.zero) go to 10
  12  xa=x2
      ya=y2
      xb=x3
      yb=y3
      x=(xa+xb)/two
      if((y-ya)*(y-yb).lt.zero) go to 1
      xa=x1
      ya=y1
      xb=x2
      yb=y2
      x=(xa+xb)/two
      if((y-ya)*(y-yb).lt.zero) go to 1
      xa=x3
      ya=y3
      xb=x4
      yb=y4
      x=(xa+xb)/two
      if((y-ya)*(y-yb).lt.zero) go to 1
      if(y.ne.y1) go to 13
      x=x1
      go to 1
   13 if(y.ne.y2) go to 14
      x=x2
      go to 1
   14 if(y.ne.y3) go to 15
      x=x3
      go to 1
   15 if(y.ne.y4) go to 16
      x=x4
      go to 1
   16 continue
      write(6,11)x,y
      go to 21
3     if(abs(f-y).lt.crit) iturn=1
      if(i.eq.1) go to 7
      dydx=(yb-ya)/(xb-xa)
      if(abs(yp-dydx).lt..2*abs(yp)) go to 7
      if((f-y)*(ya-y).lt.zero) go to 5
      xa=x
      ya=f
      go to 6
   5  xb=x
      yb=f
   6  x=(xa+xb)/two
      i=1
      go to 1
    7 if((f-y)*(ya-y).lt.zero) go to 8
      xa=x
      ya=f
      go to 9
   8  xb=x
      yb=f
   9  dydx=(yb-ya)/(xb-xa)
      if(abs(yp-dydx).lt..2*abs(yp))dydx=yp
      x=x-(f-y)/dydx
      if(iturn.eq.1) return
      i=0
      go to 1
   10 if((y-ya)*(y-yb).lt.zero) go to 1
      go to 12
      end
      subroutine cubic(p1,p2,p3,p4,p,a,b,c,d)
      USE precision
      USE physical_constants
      implicit none
      INTEGER :: k12, k23, k34
      REAL(rprec) :: d, c, b, a, p4, p3, p2, p1
      REAL(rprec) :: p1p2, p2p3, p1p3, p1p4, p2p4, p3p4
      REAL(rprec) :: pp1, pp2, pp3, pp4, p
cray  lcm (p1),(p2),(p3),(p4)
      p1p2=p1-p2
      k12=sign(1.5_dbl,p1p2)
      if(abs(p1p2).lt.zero) go to 1
      p2p3=p2-p3
      k23=sign(1.5_dbl,p2p3)
      if(abs(p2p3).lt.zero) go to 2
      p1p3=p1-p3
      if(abs(p1p3).lt.zero) go to 1
      p1p4=p1-p4
      if (abs(p1p4).lt.zero) go to 3
      p2p4=p2-p4
      if (abs(p2p4).lt.zero) go to 4
      p3p4=p3-p4
      k34=sign(1.5_dbl,p3p4)
      if(abs(p3p4).lt.zero) go to 4
      if (iabs(k12+k23+k34).eq.3)  go to 5
      if (iabs(k23+k34).eq.2) go to 1
      if (iabs(k12+k23).eq.2) go to 4
5     pp1=p-p1
      pp2=p-p2
      pp3=p-p3
      pp4=p-p4
      a=pp2*pp3*pp4/(p1p2*p1p3*p1p4)
      b=-pp1*pp3*pp4/(p1p2*p2p3*p2p4)
      c=pp1*pp2*pp4/(p1p3*p2p3*p3p4)
      d=-pp1*pp2*pp3/(p1p4*p2p4*p3p4)
      return
1     a=zero
      call parab(p2,p3,p4,p,b,c,d)
      return
2     a=zero
      b=1.
      c=zero
      d=zero
      return
3     if (iabs(k12+k23).eq.2) go to 4
      go to 1
4     d=zero
      call parab(p1,p2,p3,p,a,b,c)
      return
      end
      subroutine parab(p1,p2,p3,p,a,b,c)
      USE precision
      USE physical_constants
      implicit none
      INTEGER :: k12, k23, k34
      REAL(rprec) :: d, c, b, a, p4, p3, p2, p1
      REAL(rprec) :: p1p2, p2p3, p1p3, p1p4, p2p4, p3p4
      REAL(rprec) :: pp1, pp2, pp3, pp4, p

cray  lcm (p1),(p2),(p3)
      p1p2=p1-p2
      p2p3=p2-p3
      p1p3=p1-p3
      pp1=p-p1
      pp2=p-p2
      pp3=p-p3
      k12=sign(1.5_dbl,p1p2)
      if(abs(p1p2).lt.zero) go to 1
      k23=sign(1.5_dbl,p2p3)
      if(abs(p2p3).lt.zero) go to 2
      if(abs(p1p3).lt.zero) go to 1
      if (iabs(k12+k23).ne.2) go to 1
      a=pp2*pp3/(p1p2*p1p3)
      b=-pp1*pp3/(p1p2*p2p3)
      c=pp1*pp2/(p1p3*p2p3)
      return
1     a=zero
      b=pp3/p2p3
      c=-pp2/p2p3
      return
2     a=pp3/p1p3
      b=zero
      c=-pp1/p1p3
      return
      end
