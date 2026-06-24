C     Last change:  SBR  20 May 2009    2:37 pm

C**********************************************************************      
c*****                                                 ****************
c*****  Program for showing the radial distribution    ****************
c*****                                                 ****************
c**********************************************************************      
c     ------- DEFINITIONS ---------------------------------------------
c
c     natoms = total number of atoms
c     rad(i,j) = distance   matrix
c     rv(xi,j) = atomic positions matrix
c     itype(i) = type of the atom number i 
c     xmed,ymed,zmed = Punto medio de la distribucion     
c     y(m,1) = Set of radial numbers
c     y(m,2) = Pair correlation function (Total)
c     y(m,3) = Pair correlation function for type 1 atoms
c     y(m,4) = Pair correlation function for type 2 atoms 
c     zenbaki1 = Number of the type 1 atom nearest to the mass-center
c     zenbaki2 = Number of the type 2 atom nearest to the mass-center       
c
c**********************************************************************

	parameter (nbin=100,natmax=3000,neimax=600,ndbg=1197)
      dimension y(420,4),t(3,3),eval(3),evec(3,3),y3(420),y4(420)
	dimension nbig(natmax),nsmall(natmax),d(natmax),ws(420),wc(420)

	character*80 header
	integer natoms,ntypes,zenbaki1,zenbaki2
      real alat,latty,erradmin1,erradmin2
      REAL*8 ti,temp
      dimension rv(6,natmax),itype(natmax)
     1          ,errad(natmax)

      dimension amass(10),ielement(10),nqn(10,3),nzeta(10),zs(10,3)
     1          ,as(10),ap(10),ad(10),alfatf(10),rad(natmax,natmax)


      open (22,file='distri.d',status='old')
      open (1,file='analis.out',status='unknown')
c      	open (2,file='correl.out',status='unknown')
      open (7,file='distri.out',status='unknown')
     	open (8,file='archivo.xyz',status='unknown')
c	open (9,file='gnup2.out',status='new')


      read(22,9501) header
9501  format(a80)
      write(1,9001)header
9001  format(/,3x,a80,/)
      read(22,*) ao
      write(1,*) 'ao=', ao
      read(22,5002)
5002  format (//)
      read(22,9502) natoms
9502  format(i5)

      read(22,5003)
5003  format (/)

c     Adaptado para CONTCAR con un s¢lo tipo de  tomos
      ntypes = 1
c************************************************************************
c  check that the number of atoms is not too large
c************************************************************************

      if (natoms.gt.natmax) then
         write(1,9901) natoms,natmax
9901     format(1x,'natoms=',i6,' is greater than natmax=',i6)
         stop
       end if

      write(1,9002)natoms,ntypes
9002  format(1x,i10,' atoms',i10,' particle types')
       
9503  format(3e25.16)
      
c     Datos para el Pt

      amass(1) = 195.084
      ielement(1) = 79
      rcut= 4.1
      rnn= 2.9

9504  format(e25.16,i10)


c***************************************************************************
c     Reading with the new format
c***************************************************************************
       do j=1,natoms
       itype(j) = 1
       end do
       
       read(22,9505) ((rv(i,j),i=1,3),j=1,natoms)
       write(1,3000)(itype(j),(rv(i,j),i=1,3),j=1,natoms)

3000   format(i10,3e25.16)

9505   format(3e21.16)

       write(1,*) 'coordenadas directas'
       do j = 1,natoms
       do i = 1,3
       rv(i,j) = rv(i,j)*ao
       end do
       end do
       
       write(1,3000)(itype(j),(rv(i,j),i=1,3),j=1,natoms)

c****************************************************************************
c     print out types
c****************************************************************************

1000  continue
      write(1,9102)
 9102 format('   type  element      amass  ',/,
     1       '   ----  -------    ---------')
      write(1,9103)(i,ielement(i),amass(i),i=1,ntypes)
 9103 format(1x,i4,i9,2x,g11.4)


c***************************************************************************
c     calculates distance matrix    
c***************************************************************************

        do i=1,natoms-1
        rad(i,i)=0.
        do j=1+1,natoms

        rad(i,j)=((rv(1,i)-rv(1,j))**2+(rv(2,i)-rv(2,j))**2
     #           +(rv(3,i)-rv(3,j))**2)**.5

        rad(j,i)=rad(i,j)
        end do
        end do
        rad(natoms,natoms)=0.
        
c***********************************************************************
c      C lculo de distancia promedio de enlace
c***********************************************************************
       d_prom = 0.0
       n_enl = 0
       
       do i=1,natoms
       do j= i+1,natoms
       
       if(rad(i,j).lt.rnn) then
       d_prom = d_prom + rad(i,j)
       n_enl = n_enl + 1
       end if
       
       end do
       end do
       
       d_prom = d_prom/float(n_enl)
       
       write(1,5000) n_enl,d_prom
5000   format(1x,'nro. de enlaces',i10,'distancia de enlace promedio=',
     1 f15.4)

c**************************************************************************
c      escritura de matriz de distancias    
c**************************************************************************

        write(1,1225)
1225    format(//,1x,'******* positions',//)
        write(1,1222)(i,rv(1,i),rv(2,i),rv(3,i),i=1,natoms)
1222    format(1x,i4,3f15.3)
        write(1,1200)
1200    format(//,1x,'******* matriz de distancias')
c       call pegleg(rad,natoms,natmax)

c************************************************************************
c       elipsoide de inercia - centro de masa   
c************************************************************************

        xmed=0.
        ymed=0.
        zmed=0.

        do i=1,natoms
        xmed=xmed+rv(1,i)
        ymed=ymed+rv(2,i)
        zmed=zmed+rv(3,i)

        t(1,1)=t(1,1)+rv(1,i)**2
        t(2,2)=t(2,2)+rv(2,i)**2
        t(3,3)=t(3,3)+rv(3,i)**2
        t(1,2)=t(1,2)+rv(1,i)*rv(2,1)
        t(1,3)=t(1,3)+rv(1,i)*rv(3,1)
        t(2,3)=t(2,3)+rv(2,i)*rv(3,1)
        end do
        xmed=xmed/float(natoms)
        ymed=ymed/float(natoms)
        zmed=zmed/float(natoms)
        t(1,1)=t(1,1)/float(natoms)
        t(2,2)=t(2,2)/float(natoms)
        t(3,3)=t(3,3)/float(natoms)
        t(1,2)=t(1,2)/float(natoms)
        t(1,3)=t(1,3)/float(natoms)
        t(2,3)=t(2,3)/float(natoms)

        t(2,1)=t(1,2)
        t(3,1)=t(1,3)
        t(3,2)=t(2,3)

        write(1,*)' '
        write(1,*)' '
        write(1,*) '******* centro de masa'
        write(1,2222)xmed,ymed,zmed
2222    format(//,'  <x>=',f10.5,'  <y>=',f10.5,'  <z>=',f10.5,//)

c        write(1,*) '******* momento de inercia'
c        call pegleg(t,3,3)
c        call evcsf(3,t,3,eval,evec,3)
c        write(1,*) '******* autovalores'
c        write(1,*)' '
c        write(1,1203) eval(1),eval(2),eval(3)
c1203    format(/,'  eval(1)=',f10.5,'  eval(2)=',f10.5,
c     #           '  eval(3)=',f10.5,///)
c        write(1,*) '******* matriz de rotacion'
cc*****       normaliza los autovectores ********
c        do i=1,3
c        onori=0.
c        do j=1,3
c        onori=onori+evec(j,i)**2
c        end do
c        do j=1,3
c        evec(j,i)=evec(j,i)/onori**.5
c        end do
c        end do
c        call pegleg(evec,3,3)

c*************************************************************************
c        calculo del atomo de tipo 1 mas proximo al centro de masas
c*************************************************************************

      erradmin1=40.0

      do 6000 i=1,natoms

      errad(i)=((rv(1,i)-xmed)**2+(rv(2,i)-ymed)**2
     1     +(rv(3,i)-zmed)**2)**.5

	if (errad(i).LT.erradmin1.and.itype(i).eq.1) then
         erradmin1=errad(i)
         zenbaki1=i
         else
         goto 6000
      endif
 6000 continue         
      write(1,6001)zenbaki1,erradmin1
 6001 format('Number of type 1 atom',I8,/,'Distance from the
     1 mass-center:',f12.6)


c*************************************************************************
c        calculo del atomo de tipo 2 mas proximo al centro de masas
c*************************************************************************

      erradmin2=40.0

      do 7000 i=1,natoms

      errad(i)=((rv(1,i)-xmed)**2+(rv(2,i)-ymed)**2
     1     +(rv(3,i)-zmed)**2)**.5

      if (errad(i).LT.erradmin2.and.itype(i).eq.2) then
         erradmin2=errad(i)
         zenbaki2=i
         else
         goto 7000
      endif
 7000 continue         
      write(1,7001)zenbaki2,erradmin2
 7001 format('Number of type 2 atom',I8,/,'Distance from the
     1 mass-center:',f12.6)

c************************************************************************
c     histograma de distancias - funcion de correlacion de pares  
c************************************************************************
	if (erradmin1.le.erradmin2) then
	   na=zenbaki1
	   else
	   na=zenbaki2
	endif
        sigma=.1
        c=((2.*3.14159)**(-.5))/sigma
        n1=zenbaki1
        n2=zenbaki2

        do i=1,420
        y(i,2)=0.
        y(i,3)=0.
        y(i,4)=0.
        end do

c*****Pair correlation function for all types of atoms*********************

        do 1 m=7,420
        y(m,1)=float(m-1)*37.8/420.
c        do 2 n=1,natoms
        do 3 j=1,natoms
10      format(1x,2i4,2f10.3)
        if(abs(rad(na,j)-y(m,1)).gt..4) go to 3
        y(m,2)=y(m,2)+c*exp(-1.*(y(m,1)-rad(na,j))**2/2.
     #         /sigma**2)
3       continue
c2       continue
	 y(m,2)=y(m,2)/(4.*3.1416*y(m,1)**2)
1       continue

C*****Pair correlation function for type 1 atoms****************************

        do 7100 m=7,420
           do 7101 j=1,natoms
           if(itype(j).ne.1) goto 7101
           if(abs(rad(na,j)-y(m,1)).gt..4) go to 7101
           y(m,3)=y(m,3)+c*exp(-1.*(y(m,1)-rad(na,j))**2/2.
     #         /sigma**2)
7101       continue
	   y(m,3)=y(m,3)/(4.*3.1416*y(m,1)**2)
7100    continue

C*****Pair correlation function for type 2 atoms****************************

        do 7200 m=7,420
           do 7201 j=1,natoms
           if(itype(j).ne.2) goto 7201
           if(abs(rad(na,j)-y(m,1)).gt..4) go to 7201
           y(m,4)=y(m,4)+c*exp(-1.*(y(m,1)-rad(na,j))**2/2.
     #         /sigma**2)
7201       continue
	   y(m,4)=y(m,4)/(4.*3.1416*y(m,1)**2)
7200    continue


c        write(2,987)(y(m,1),y(m,2),m=1,420)
        write(7,988)(y(m,1),y(m,2),y(m,3),y(m,4),m=1,420)
c        write(8,987)(y(m,1),y(m,3),m=1,420)
c        write(9,987)(y(m,1),y(m,3),m=1,420)
	write(1,1201)
987     format(2f15.5)
988     format(4f15.5)
1201    format(//,1x,'******* funcion de correlacion de pares
     #, plot 1',//)
        write(1,36)(y(m,1),y(m,2),m=1,420)
366     format(5(1x,f5.2,f7.1))
36      format(5(1x,f5.2,f7.1,' | '))

c***************************************************************************
c     factor de estructura   
c***************************************************************************

        write(1,1205)
1205    format(//,1x,'******* factor de estructura, plot 3')
        do i=1,420
        y3(i)=y(i,2)
        w0=w0+y3(i)
        end do


        write(1,1211)
1211    format(/,'       k               f(k)     ',/)
        do 300 n=1,420
        ws(n)=0.
        wc(n)=0.
        akn=FLOAT(n-1)*3.14/12./2.
        y4(n)=akn
        DO 400 J=1,420
        xj=float(j-1)*12./100.
        Ws(n)=Ws(n)+y3(J)*sin(akn*xj)
        Wc(n)=Wc(n)+y3(J)*cos(akn*xj)
400     CONTINUE
300     CONTINUE        
c        WRITE(4,987)(y4(i),(Ws(i)**2+wc(i)**2)**.5/w0,i=1,420)
        WRITE(1,336)(y4(i),(Ws(i)**2+wc(i)**2)**.5/w0,i=1,420)
3366    format(5(1x,f5.2,f8.4))
336     format(5(1x,f5.2,f8.4,' |'))

c*****************   coordinacion      ********************************

        do i=1,natoms
        do j=1,natoms
        if(rad(i,j).lt.rcut.and.rad(i,j).gt.0.)nbig(i)=nbig(i)+1
        if(rad(i,j).lt.rnn.and.rad(i,j).gt.0.)nsmall(i)=nsmall(i)+1
        end do
        end do
        write(1,1230)
1230    format(1x,'******* coordinacion',//)
        write(1,*)'     n         nsmall     nbig'
        write(1,*)' '
        write(1,1220)(k,nsmall(k),nbig(k),k=1,natoms)
1220    format(2(1x,3i10,5x))

c****************   histograma de coordinacion     **********************

        do i=1,100
        y(i,2)=0.
        end do
        do 11 m=1,15
        y(m,1)=float(m)
        do 22 n=1,natoms
        if(abs(float(nsmall(n))-y(m,1)).gt..5) go to 33
        y(m,2)=y(m,2)+1.
33      continue
22      continue
11      continue
c        write(3,987)(y(m,1),y(m,2),m=1,15)
        write(1,2201)
2201    format(//,1x,'******* histograma de nn, plot 2',//)
        write(1,136)(y(m,1),y(m,2),m=1,15)
136     format(2(1x,f15.5))

c****************  distancias respecto del CM    *********************

        do i=1,natoms
        d(i)=((rv(1,i)-xmed)**2+(rv(2,i)-ymed)**2+(rv(3,
     #              i)-zmed)**2)**.5
        end do  
c       atom closest to CM
        amin=100.
        do i=1,natoms
        if(d(i).le.amin)then
        amin=d(i)
        imin=i
        end if
        end do

        write(1,2301)
2301    format(//,1x,'******* distancias alcentro',//)
        write(1,2303)(i,rad(i,imin),i=1,natoms)
2303    format(5(1x,i4,f8.3,' |'))

        do i=1,420
        y(i,2)=0.
        end do

        do 51 m=1,420
        y(m,1)=float(m-1)*10./420.
        do 52 n=1,natoms
        if(abs(rad(n,imin)-y(m,1)).gt..4) go to 52
        y(m,2)=y(m,2)+c*exp(-1.*(y(m,1)-rad(n,imin))**2/2.
     #         /sigma**2)
52      continue
51      continue

c        write(7,987)(y(m,1),y(m,2),m=1,420)
        write(1,1301)
1301    format(//,1x,'******* histograma de distancias al centro
     #, plot 4',//)
        write(1,36)(y(m,1),y(m,2),m=1,420)

c        call plot3d(rv,natoms,43)
        end

C========================================================================
c        SUBROUTINE PEGLEG(A,N,NDIM)
C========================================================================


c        DIMENSION A(NDIM,NDIM)

c17      FORMAT(1H0/1H0/)
c18      FORMAT(I4,3X,8F8.4)
c19      FORMAT(1H0/1H0,3X,8(5X,I3))

c        KITE=0
c20      LOW=KITE+1
c        KITE=KITE+8
c        KITE=MIN0(KITE,N)

c        WRITE(1,19)(I,I=LOW,KITE)

c        DO 32 I=1,N
c32      WRITE(1,18)I,(A(I,J),J=LOW,KITE)

c        IF(N-KITE)40,40,20

c40      WRITE(1,17)
c        RETURN
c        END

