      program forcing_trans

      REAL*4 TmpArr (    360,    160  )
      real*4 fieldout(360, 160)
      character*80 infile, outfile

C

      type *,'Start program'
      type *,'Give output file name :'
      read (*,'(a)') outfile

      outrec = 0
      uuniti=1
      uunito=2

        OPEN (uUnito,FILE=outfile,form='unformatted',
     &        access='direct',recl=360*160*4,status='new')

      type *,'Opened output file'

10    type *,'Give input file name :'
      read (*,'(a)') infile
      if(infile.eq.'END')STOP

      type *,'Give nrec:'
      read *,nrec1, nrec2

      if(infile(1:3).eq.'END') goto 100

        OPEN (uUniti,FILE=infile,form='unformatted',
     &        access='direct',recl=360*160*4,status='OLD')

      type *,'Opened input file'

      do fluxRec = nrec1, nrec2
      READ (uuniti,rec=fluxrec,end=100) TmpArr

      outrec = outrec + 1

      do i=1,360
      do j=1,160
c       TmpArr(i,j) = TmpArr(i,j)/10.
        if(TmpArr(i,j).lt.-9990.) TmpArr(i,j) = 0.
        fieldout(i,j)=TmpArr(i,j)
      enddo
      enddo

c     if(fluxRec.eq.20.or.fluxRec.eq.31.or.fluxRec.eq.41) then
c        do i=1,180
c        do j=1,80
c          fieldout(i,j) = -9999.
c        enddo
c        enddo
c        goto 50
c     endif
c     type *,fluxrec
c     call one2two_r4(TmpArr,fieldout,0.)
50    continue
c     write(uunito) 180, 80, -79., 79., 1., 359.
c     print *,fieldout
      write(uunito,rec=outrec) fieldout
      enddo

      goto 10

100   close(uuniti)
!
!     fluxrec = 120 -12
!     READ (uunito,rec=fluxrec,end=100) TmpArr
!     fluxrec = 120
!     write(uunito,rec=outrec) TmpArr

      close(uunito)

      stop
      END

      SUBROUTINE ONE2TWO_R4 (fieldin,fieldout,spval)
C
C     ======== Routine arguments ======================
      real*4 fieldin(360,160)
      real*4 fieldout(180,80)
      REAL spval
CEndofinterface
C     ========= Local variables =========================
      INTEGER I,J,k,l,ii,jj
      real av,nav

      k=0
      l=0
      do j=1,160,2
      l=l+1
      k=0
      do i=1,360,2
      k=k+1
        av = 0.
        nav = 0.
        do 100 jj=j,j+1
        do 100 ii=i,i+1
          if(fieldin(ii,jj).lt.spval+1.) fieldin(ii,jj)=spval
          if(fieldin(ii,jj).ne.spval) then
              av=av + fieldin(ii,jj)
              nav = nav+1.
          endif
 100    continue
      if(nav.gt.0) then
      fieldout(k,l) = av/nav
      else
      fieldout(k,l) = -9999.
      endif

      enddo
      enddo

C
      RETURN
      END

