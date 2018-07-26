      parameter (Nx=360,Ny=180,Mrey=421,nweeks=52,nyrs=11)
      real sst(Nx*Ny),sstav(nweeks,Nx*Ny),sstanom(Nx*Ny)
      integer npts(nweeks,Nx*Ny)
      integer istart(nyrs)
      data istart/1,53,105,157,209,262,314,366,418,470,522/    ! 93-03   (97 has 53 weeks)
      nxy=Nx*Ny

      open(20,file='REYv2_sst_93-03',form='unformatted',access=
     *'direct',recl=Nx*Ny*4,status='old')


      DO iweek=1,52
        do iyr=	1,nyrs
	  irec=istart(iyr)+iweek-1
	    if(irec.le.Mrey)then
              read(20,rec=IREC)sst
              do ij=1,Nxy
		if(sst(ij).ne.-9999.)then
                  sstav(iweek,ij)=sstav(iweek,ij)+sst(ij)
	          npts(iweek,ij)=npts(iweek,ij)+1
		endif
              enddo
	    endif
	enddo
      ENDDO

      open(30,file='REYv2_sst_93-03_wklymean',form='unformatted',
     *     access='direct',recl=Nx*Ny*4,status='unknown')
      DO iweek=1,52
	do ij=1,Nxy
	  if(npts(iweek,ij).ne.0)then
	    sstav(iweek,ij)=sstav(iweek,ij)/npts(iweek,ij)
	  else
	    sstav(iweek,ij)=-9999.
	  endif
	enddo
	write(30,rec=iweek)(sstav(iweek,ij),ij=1,Nxy)
      ENDDO
      endfile(30)

      open(30,file='REYv2_sst_93-03_wklypanom',form='unformatted',
     *     access='direct',recl=Nx*Ny*4,status='unknown')
      DO iweek=1,52
        do iyr=	1,nyrs
	  irec=istart(iyr)+iweek-1
	    if(irec.le.Mrey)then
              read(20,rec=IREC)sst
              do ij=1,Nxy
		if(sst(ij).ne.-9999.)then
                  sstanom(ij)=sst(ij)-sstav(iweek,ij)
		else
		  sstanom(ij)=-9999.
		endif
              enddo
	    endif
	write(30,rec=irec)sstanom
	enddo
      ENDDO

      irec=261    ! 53rd week of 97
      read(30,rec=irec-1)sstanom
      read(30,rec=irec+1)sst
      do ij=1,Nxy
	if(sstanom(ij).ne.-9999..and.sst(ij).ne.-9999.)then
	  sstanom(ij)=(sstanom(ij)+sst(ij))/2.
	else
	  sstanom(ij)=-9999.
	endif
      enddo
      write(30,rec=irec)sstanom
      endfile(30)

      

      stop
      end
