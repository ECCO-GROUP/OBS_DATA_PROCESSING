      parameter(Nx=360,Ny=180)
      parameter(Mrey=529,MTP=364)
      parameter(MTPstrt=12)

      real TPtme93(MTP),Reytme93(Mrey),x(2),y(2)
      real sst1(Nx,Ny),sst2(Nx,Ny),sst(Nx,Ny)

      open(20,file='REYv2_sst_93-03_wklypanom',form='unformatted',
     *access='direct',status='old',recl=Nx*Ny*4)
      open(30,file='REY_sstpanom_11-364',form='unformatted',
     *access='direct',status='unknown',recl=Nx*Ny*4)

! set up TP times   (XI)  days
      TPtme93(11)=5.25
      do irep=12,MTP
	TPtme93(irep)=TPtme93(irep-1)+9.9156
      enddo

! set up REYNOLDS times (X)  days
      Reytme93(1)=6.5
      do irec=2,Mrey
	Reytme93(irec)=Reytme93(irec-1)+7
      enddo

      IREC1=2
      DO NR=MTPstrt,MTP

        XI=TPtme93(NR)
	do irec=IREC1,Mrey
	write(6,*)NR,XI,Reytme93(irec-1),Reytme93(irec)
	if(XI.ge.Reytme93(irec-1).and.XI.lt.Reytme93(irec))then
	  x(1)=Reytme93(irec-1)
	  x(2)=Reytme93(irec)
	  IREC1=irec
	  go to 100
	endif
	enddo
	STOP1

  100  read(20,rec=irec-1)sst1
       read(20,rec=irec)sst2
       do j=1,Ny
       do i=1,Nx
	 if(sst1(i,j).eq.-9999..or.sst2(i,j).eq.-9999.)then
	   sst(i,j)=-9999.
	 else
	   y(1)=sst1(i,j)
	   y(2)=sst2(i,j)
	   call xinterp(x,y,XI,2,sst(i,j))
	 endif
       enddo
       enddo
       write(30,rec=NR)sst
	enddo

	write(6,*)' JOB COMPLETE'
	stop
	end
      subroutine xINTERP (X,Y,XI,N,yinterp)
c !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
c  MUST BE MONOTONIC INCREASING!!!!!!!!!!!!!!!!!!!
c  note :  real*4 version
c !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      REAL X(n),Y(n),xi,yinterp
C
      if(x(2).lt.x(1))then
	 write(6,*)' STOPPiNG XINTERP....not monotonic increasing'
	 stop
      endif
      IF (XI.LT.X(1)) GO TO 800
C
      DO 100 I=2,N
      IF (XI.LE.X(I)) GO TO 400
 100  CONTINUE
      yinterp=Y(N)+((Y(N)-Y(N-1))*(XI-X(N)))/(X(N)-X(N-1))
      write(6,*)' outside range ..XI > x(n)'
      RETURN
C
  400 yinterp=(Y(I-1)*(X(I)-XI)-Y(I)*(X(I-1)-XI))/(X(I)-X(I-1))
      RETURN
C
 800  yinterp=Y(1)+((Y(2)-Y(1))*(XI-X(1)))/(X(2)-X(1))
      write(6,*)' outside range ..XI < x(1)'
      RETURN
      END
