c   program to plot color contours for Detlef.
c  This version plots from 90E  to 450E (puts Pacific Ocean in middle of page)
c        strtlon = first datapoint  after 90 E
c
      parameter (Miz=1600,Mjz=800)
      logical coast,solid,ifcolor,ifrmv,ifpaths,nobar,iffill
     *,iftimew,bar_reg
      character word*8,head1*100,head2*100,proj*2,fin*80,info*80
      character timeword*11
      integer strlol,endlol,dellol,strlal,endlal,dellal

c
c
      real g(Miz*Mjz),g1(Miz,Mjz),g2(Miz,Mjz)
      real xcal1(10000),ycal1(10000)
c
      common /plot/ step,coast,solid,plat,plon,iz,jz,xlonmi,xlonma,
     %              ylatmi,ylatma
      common /map/ strlat,endlat,strlon,endlon,zmax,zmin,zincr,nchop,
     $             scale
     *,strlol,endlol,dellol,strlal,endlal,dellal,nl
      common /mas/ nch1,nch2,mtot1,mtot2,eincr,emax,nbad,nobar
     *,iftimew,bar_reg
c
      data step,plat,plon /45.,0.,270./,
ck      data step,plat,plon /45.,0.,335./,
     %     zincr,zmax,zmin /10.,0.,0./,
ck     %     strlat,endlat,strlon,endlon /-90.,90.,90.,450./,
     %     strlat,endlat,strlon,endlon /-90.,90.,270.,400./,
     %     proj /'CE'/

!     mtot1=7473
!     open(55,file='/data37/king/A/ATOC_baseline/Wunsch.paths/atoc_all',
!    &form='formatted',status='old')
!     do i=1,mtot1
!        read(55,*)dist,ycal1(i),xcal1(i)
! if(xcal1(i).lt.90.)xcal1(i)=xcal1(i)+360
!     enddo
      mtot1=0

      mtot2=0
      nbad=0
      solid=.true.
      ifcolor=.true.
c
c
      read(5,*)strlat,endlat,strlon,endlon
      read(5,*)strlol,endlol,dellol
      read(5,*)strlal,endlal,dellal,nl
      step=dellol
      plat=(strlat+endlat)/2
      plon=(strlon+endlon)/2

      read (5,'(a)')fin
      write(6,'(a)')fin
      read (5,'(a)') head1
      read (5,'(a)')word 
      read (5,*)nskip    ! no of arrays to skip reading 
      write(6,*)' reading field number ',nskip+1
      read (5,*)ifcolor,ncolor,ifrmv,ifpaths,nobar,iffill,iftimew,
     *bar_reg
      write(6,*)' nobar = ',nobar
      read (5,*)gmin,gmax,gincr,scale
c
      if(.not.ifpaths)mtot1=0
      nch1=index(head1,'$')-1
c
      open(10,file=fin,form='unformatted',status='old',access=
     *'direct',recl=360*180*4)
c     if(nskip.ne.0)then
c       do nn=1,nskip
c         read(10,rec=nn)(g(ij),ij=1,iz*jz)
c       enddo
c     endif
c
      do j=1,mjz
      do i=1,miz
        g1(i,j)=-9999.
        g2(i,j)=-9999.
      enddo
      enddo
      do ij=1,miz*mjz
	g(ij)=-9999.
      enddo


c     read(10)iz,jz,ylatmi,ylatma,xlonmi,xlonma
      iz=360
      jz=180
      ylatmi=-89.5
      ylatma=89.5
      xlonmi=.5
      xlonma=359.5

      dy=(ylatma-ylatmi)/(jz-1)
      write(6,*)' orig ',iz,jz,ylatmi,ylatma,xlonmi,xlonma
      do i=1,1000
        ymin=ylatmi-i*dy
        if(ymin.lt.strlat)go to 4
      enddo
    4 ymin=ymin+dy
      do i=1,1000
        ymax=ylatma+i*dy
        if(ymax.gt.endlat)go to 5
       enddo
    5  ymax=ymax-dy
       j1=((ylatmi-ymin)/dy+1.01)
       j2=((ylatma-ymin)/dy+1.01)
       jz=((ymax-ymin)/dy+1.01)
       ylatma=ymax
       ylatmi=ymin
       write(6,*)' lat pos ',ylatmi,ylatma,j1,j2,jz
       read(10,rec=nskip+1)((g1(i,j),i=1,iz),j=j1,j2)
 
c  fill in last long. position if necessary
      dx=(xlonma-xlonmi)/(iz-1)
      write(6,*)' ddeg = ',dx
      iiz=(360./dx+1.1)
      if(iz+1.eq.iiz)then
       write(6,*)' changing iz ',iz,iiz
       iz=iiz
       xlonma=xlonma+dx
       do j=j1,j2
        g1(iz,j)=g1(1,j)
       enddo
      endif

CK   extra
c     do j=j1,j2
c     do i=1,iz
c      if(g1(i,j).ne.-9999.)then
c  g1(i,j)=30
c      else
c  g1(i,j)=-30
c      endif
c     enddo
c     enddo
CK   extra
 
c  look for data holes
      if(iffill)then
      do j=2,jz-1
      do i=2,iz-1
        if(g1(i,j).eq.-9999.)then
           if(g1(i-1,j).ne.-9999..and.g1(i+1,j).ne.-9999..and.
     *      g1(i,j-1).ne.-9999..and.g1(i,j+1).ne.-9999.)then
           g1(i,j)=(g1(i-1,j)+g1(i+1,j)+g1(i,j-1)+g1(i,j+1))/4
           endif
        endif
      enddo
      enddo
      endif
 
c  now check lat limits
      j1=1
      if(strlat.gt.ylatmi)then
	 j1=(strlat-ylatmi)/dx+1
      endif
      j2=jz
      if(endlat.lt.ylatma)then
	 j2=(endlat-ylatmi)/dx+1
      endif
      ylatma=ylatmi+(j2-1)*dx
      ylatmi=ylatmi+(j1-1)*dx
      jz=j2-j1+1
      write(6,*)' new y plot limits ',j1,j2,jz,ylatmi,ylatma

c  now setup lon limits
      ipos=((strlon-xlonmi)/dx)+1
      strtlon=(ipos-1)*dx +xlonmi
      if(strtlon.lt.strlon)strtlon=strtlon+dx
      write(6,*)' strtlon ',strtlon
      newiz=(endlon-strtlon)/dx+1
      izjz=0
      Do j=j1,j2
      do i=1,iz
        if(g1(i,j).ne.-9999.)g1(i,j)=g1(i,j)*scale
        xlon=(i-1)*dx+xlonmi
        if(xlon.lt.strtlon)xlon=xlon+360
        xpos=(xlon-strtlon)/dx
        ipos=(xpos+1.1)
        g2(ipos,j)=g1(i,j)
      enddo
      if(endlon-strtlon.ge.360.)g2(newiz,j)=g2(1,j)
      do i=1,newiz
       izjz=izjz+1
       g(izjz)=g2(i,j)
      enddo
      enddo
      xlonmi=strtlon
      iz=newiz
      xlonma=xlonmi+(iz-1)*dx
c
      gmean=0
      gvar=0
      ng=0
      do ij=1,iz*jz
       if(g(ij).ne.-9999.)then
        gmean=gmean+g(ij)
	gvar=gvar+g(ij)*g(ij)
        ng=ng+1
       endif
      enddo
      gmean=gmean/ng
      gvar=gvar/ng-gmean*gmean
      sd=sqrt(gvar)
      ggmean=gmean
      write(6,*)' data mean,sd = ',gmean,sd,ng
      if(ifrmv)then
       write(6,*)' mean removed.'
       head2='mean removed$'
c      head2=' $'
      else
       gmean=0
       head2=' $'
      endif
      nch2=index(head2,'$')-1
c
      zmin=1.e10
      zmax=-1.e10
      do ij=1,iz*jz
       if(g(ij).ne.-9999.)then
        g(ij)=g(ij)-gmean
        zmax=amax1(zmax,g(ij))
        zmin=amin1(zmin,g(ij))
       endif
      enddo
      write(info,123)ggmean,zmin,zmax,ng
  123 format('mean = ',f8.2,' min = ',f8.2,' max = ',f8.2,
     *  ' npts = ',i8)
      write(6,*)'data min,max= ',zmin,zmax
      if(gmin.eq.999.)then
        write(6,*)'  enter contour min,max,incr'
        read(5,*)zmin,zmax,zincr
      else
        zmin=gmin
        zmax=gmax
        zincr=gincr
      endif
      emat=1.e30

      call opngks
      call getdate(timeword)
      call domapk(g,head1,head2,word,xcal1,ycal1,proj,xcal2,ycal2,emat,
     *ifcolor,ncolor,fin,info,timeword)
      call clsgks
c
      stop
      end
c
      subroutine getdate(newword)
      character*1 time_word*24,timeword(24),newword(11)
      equivalence(time_word,timeword)
      call fdate(time_word)
 1     k=0
       do i=5,10
         k=k+1
         newword(k)=timeword(i)
       enddo
       do i=20,24
         k=k+1
         newword(k)=timeword(i)
       enddo
      return
      end
      SUBROUTINE DOMAPK (ZMAT,HEAD1,HEAD2,WORD,XCAL1,YCAL1,PROJ,XCAL2,
     $                  YCAL2,EMAT,ifcolor,ncolor,fin,info,timeword)
c
c this 	SUBROUTINE contour plots zmat data array using ncar software v3
c also does a continental global map with the continents solid filled
c and will stipple parts of the map where the error is above a value
c and will plot a symbol at locations of tide gauges
c
      parameter (nwrk=90000,nxcs=200000,niam=3000000,niai=150)
      character head1,head2,word*8,proj*2,help*8,fin*80,info*80
      character timeword*11
      logical solid,coast,ifcolor,nobar,iftimew,bar_reg
c                                                                       
      common /plot/ step,coast,solid,plat,plon,iz,jz,xlonmi,xlonma,
     %              ylatmi,ylatma
      common /map/ strlat,endlat,strlon,endlon,zmax,zmin,zincr,nchop,
     $             scale
     *,strlol,endlol,dellol,strlal,endlal,dellal,nl
      common /mas/ nch1,nch2,mtot1,mtot2,eincr,emax,nbad,nobar
     *,iftimew,bar_reg
c
      dimension zmat(*),pl1(2),pl2(2),pl3(2),pl4(2),xcal1(*),ycal1(*),
     $          rwrk(nwrk),iwrk(nwrk),xcal2(*),ycal2(*),emat(*)
      dimension iama(niam),xcs(nxcs),ycs(nxcs),iaia(niai),igia(niai),
     $          iasf(13)
      dimension u(10000),v(10000)
ck  label bar dimensions
      dimension lind(30)
      character*6 llbs(30)
      integer strlol,endlol,dellol,strlal,endlal,dellal
      character cxi3*3
c
      data pl1(2),pl2(2),pl3(2),pl4(2) /0.,0.,0.,0./, iasf /13*1/
c
c declare a routine to color the areas represented by the area map (colram)
c and one to draw grid lines masked by the area map (colrln)
c and one to shade the contours where the error is large
      external colramk,colrln,shader,colorz
      write(6,*)' in domapk...ifcolor = ',ifcolor
c
c
c turn off the clipping indicator
      call gsclip (0)
c
c set all aspect source flags to "individual"
      call gsasf (iasf)
c
c   force solid fill (fill area interior style) 
      call gsfais (1)
c
c define color indices
       if(ncolor.eq.1)call dfclrs1           ! positive run
       if(ncolor.eq.2)call dfclrs2           ! -/+ run
       if(ncolor.eq.3)call dfclrs3           ! gray scale
       if(ncolor.eq.4)call dfclrs4           ! positive run (5th floor)
       if(ncolor.eq.5)call dfclrs2a          ! -/+ anomaly run
       if(ncolor.eq.6)call dfclrs6          ! seasonal-/+ anomaly run
       if(ncolor.eq.7)call dfclrs2c          ! -/+ anomaly run
       if(ncolor.eq.8)call dfclrs1a          ! + run (took out 1 green)
       if(ncolor.eq.9)call dfclrs3a          ! other gray scale
c
c use ezmap and ezmapa to create a background
      call set (0.,1.,0.,1.,0.,1.,0.,1.,1)
c
c
c PLOTCHAR (FOR HEADER LABELS)
c use duplex character set for sans serif letters for better reproduction
      call pcseti ('CD - CHOOSE DUPLEX CHARACTER SET',1)
c
      call plchhq (.5,.75,head1(1:nch1),8.,0,0)
      call plchhq (.5,.73,head2(1:nch2),8.,0,0)
c     idate=date()
c     write (help,'(a8)') idate
c     call plchhq (.994,.994,help,8.,0,+1.)
      call plchhq (.994,.974,word,8.,0,+1.)
      if(iftimew)call plchhq ( .99,.05,timeword,6.,0,+1.)
      write (6,'(''0plot'',(a))') head1(1:nch1)
      write (6,'(1x,(a))') head2(1:nch2)
      write (6,'(1x,a8)') word
      write (6,'(1x,a8)') idate
c
c *********************************************************************
c MAPPING
c initialize the area map
      call arinam(iama,niam)
c
      pl1(1)=strlat
      pl2(1)=strlon
      pl3(1)=endlat
      pl4(1)=endlon
c
      call maproj (proj,plat,plon,0.)
ck      call mappos (0.05,0.95,0.1,.90)
      call mappos (0.05,0.95,0.2,.8)
      call mapset ('CO',pl1,pl2,pl3,pl4)
      call mapstr ('GR',step)
      call mapstc ('OU - OUTLINE DATASET','CO')
c use vertical stripes to reduce the # of points defining subareas
c a large number(~150) adds significantly to compute time, but saves
c memory.  VS>0 implies 2 groups of boundary lines
c using vs>0 means dividing the areas into subgroups and this makes
c an additional group of areas (group 1 is the group of areas created
c by continental boundaries and group 2 is the group of areas created
c by vertical stripping.  group 3 will be the group of contour lines).
c
      if(.not.ifcolor)call mapsti('VS',10)
      if(ifcolor)call mapsti('VS',20)
c
c set maplbl to draw perimeter of map and leave off meridian labels
      call mapsti('LA',0)
c
c initialize EZMAP
      call mapint
c
c   add edges to area map (retrieve, project and add boundary lines)
c adds 2 groups to area map: group 1 (G1=1) is perimeter plus set of
c projected boundary lines.  If VS>0, group 2 (G2=2) is added to the area
c map.  the purpose of group 2 is to split up areas into smaller areas
c (vertical stripes)
c
c  set up contours for color shading
      if(ifcolor)then
        call cpseti ('SET - DO-SET-CALL FLAG',0)
        call cpseti ('MAP - MAPPING FLAG',1)
        call cpsetr ('XC1 - X COORDINATE AT I=1',xlonmi)
        call cpsetr ('XCM - X COORDINATE AT I=iz',xlonma)
        call cpsetr ('YC1 - Y COORDINATE AT J=1',ylatmi)
        call cpsetr ('YCN - Y COORDINATE AT J=jz',ylatma)
        call cpsetr ('SPV - SPECIAL VALUE',-9999.)
c add contour lines to area map , then color
c set cls=1 to have conpack pick levels from zmin to zmax by zincr
        call cpseti ('CLS - CONTOUR LEVEL SELECTOR',1)
        if (zincr .ne. 0.)
     $   call cpsetr ('CIS - CONTOUR INTERVAL SPECIFIER',zincr)
        if (zmax .ne. 0. .or. zmin .ne. 0.) then
         call cpsetr ('CMN - CONTOUR LEVEL MINIMUM',zmin)
         call cpsetr ('CMX - CONTOUR LEVEL MAXIMUM',zmax)
        endif
        call cprect (zmat,iz,iz,jz,rwrk,nwrk,iwrk,nwrk)
        call cpclam(zmat,rwrk,iwrk,iama)
        call arscam(iama,xcs,ycs,nxcs,iaia,igia,niai,colorz)
        call arinam(iama,niam)
      endif
c
      call MAPBLA (iama)
c
c   preprocess the area map - 3 zeroes for no speeding up is appropriate for
c   continental outline maps.  3 ones should be ok for contour areas.
      call arpram (iama,0,0,0)
      isu=niam-(iama(6)-iama(5)-1)
      print *,' space used in area map is ',isu
c
c area map is now created, but not drawn, ARSCAM will draw and color it
c   color the map: ARSCAM scans an existing area map, extracting the 
c definitions of all the areas (calls COLRAM once for each scanned area)
      if (solid) call ARSCAM(iama,xcs,ycs,nxcs,iaia,igia,niai,colramk)
c
c reduce minimum vector length to include all points on boundary
      call mapsti ('MV',1)
c next 3 calls are equivalent to 1 call to mapdrw
c      call mapgrd
      call maplbl
c think i'll try mapgrm (masked grid lines) and i'll put coastline call
c after grid line call to see if dashed line setting remains in force
      call mapgrm (iama,xcs,ycs,nxcs,iaia,igia,niai,colrln)
      if (coast) call maplot
c
c restore minimum vector length to default
      call mapsti ('MV',4)
c
c******************************************************************
C DRAW POLYMARKERS TO MARK LOCATIONS OF TIDE GAUGES
      if (mtot1 .gt. 0) then
         do 15 m=1,mtot1
 15      call maptrn (ycal1(m),xcal1(m),u(m),v(m))
ck 15      call plchhq (u,v,':KGL:c',8.,0.,0.)
      end if
      call points(u,v,mtot1,0,0)

      if (mtot2 .gt. 0) then
         do 16 m=1,mtot2
         call maptrn (ycal2(m),xcal2(m),u,v)
 16      call plchhq (u,v,'*',8.,0.,0.)
      end if
c
C LABEL LAT AND LON GRID LINES WITH DEGREES
c  maplbm will label map correctly, but starts at strlon,strlat
      call MAPLBM (strlat,endlat,strlon,endlon,step,dellal,nl)
      go to 98765
c  the following code is for this special  case
c    *strlol,endlol,dellol,strlal,endlal,dellal
      call maptrn(strlat,xi3,u3,v3)
      s1=(endlat-strlat)*.03
      s2=(endlon-strlon)*.04
      do  i3=strlol,endlol,dellol
        j3=i3
c        j3=360-i3
        xi3=i3
        write(cxi3,'(i3)')j3
        call maptrn(strlat-s1,xi3,u3,v3)
        call plchhq(u3,v3,cxi3,.010,0.,0.)
        call maptrn(strlat-s1,xi3+s2,u3,v3)
        call plchhq(u3,v3,':KRL:.',.010,0.,1.)
      enddo
        call maptrn(strlat-s1,xi3+1.2*s2,u3,v3)
        call plchhq(u3,v3,'E',.010,0.,0.)
c
      do  i3=strlal,endlal,dellal
        write(cxi3,'(i3)')i3
        yi3=i3
        call maptrn(yi3,strlon-s2,u3,v3)
        call plchhq(u3,v3,cxi3,.010,0.,1.)
        call maptrn(yi3,strlon-.7*s2,u3,v3)
        call plchhq(u3,v3,':KRL:.',.010,0.,1.)
      enddo
        call maptrn(yi3,strlon-.3*s2,u3,v3)
        call plchhq(u3,v3,'N',.010,0.,1.)
c  add ticmarks
      do  i2=strlol,endlol,dellol
        xi2=i2
        call maptrn(strlat,xi2,u1,v1)
        call maptrn(strlat-1.,xi2,u2,v2)
        call line(u1,v1,u2,v2)
      enddo
      do  i3=strlal,endlal,dellal
        yi3=i3
        call maptrn(yi3,strlon,u1,v1)
        call maptrn(yi3,strlon-1.,u2,v2)
        call line(u1,v1,u2,v2)
        call maptrn(yi3,endlon,u1,v1)
        call maptrn(yi3,endlon+1.,u2,v2)
        call line(u1,v1,u2,v2)
      enddo
98765 continue
c
c *****************************************************************
c CONTOURING
c
c tell CONPACK not to do SET call (already done)
c use mapping functioon 1 (EZMAP) background, and
c what range of X and Y coordinates to send into the mapping function
c the X coordinate will be longitude from xlonmi to xlonma
c the Y coordinate will be latitudes from ylatmi to ylatma
c
      call cpseti ('SET - DO-SET-CALL FLAG',0)
      call cpseti ('MAP - MAPPING FLAG',1)
      call cpsetr ('XC1 - X COORDINATE AT I=1',xlonmi)
      call cpsetr ('XCM - X COORDINATE AT I=iz',xlonma)
      call cpsetr ('YC1 - Y COORDINATE AT J=1',ylatmi)
      call cpsetr ('YCN - Y COORDINATE AT J=jz',ylatma)
      call cpsetr ('HLS - HIGH/LOW LABEL SIZE',0.)            !.007)
      call cpsetr ('HLW - HIGH/LOW LABEL WHITE SPACE WIDTH',.004)
      call cpseti ('NSD - NUMBER OF SIGNIFICANT DIGITS',2)
      call cpsetr ('ILS - INFORMATIONAL LABEL SIZE',0.)        !.006)
      call cpsetr ('ILW - INFORMATIONAL LABEL WHITE SP WIDTH',.003)
      call cpsetr ('ILY - INFORMATIONAL LABEL Y-COORDINATE',-.5)
      if(.not.ifcolor)then
        call cpsetr ('LLS - LINE LABEL SIZE',.0075)
      else
        call cpsetr ('LLS - LINE LABEL SIZE',0.)
      endif
      call cpsetr ('LLW - LINE LABEL WHITE SPACE WIDTH',.004)
      call cpsetc ('HLT - HIGH/LOW LABEL TEXT STRINGS',
     $                    '$ZDV$''$ZDV$')
      call cpseti ('LLP - LINE LABEL POSITIONING',3)
      call cpsetr ('SPV - SPECIAL VALUE',-9999.)
c 
c set contour levels
c set cls=1 to have conpack pick levels from zmin to zmax by zincr
      call cpseti ('CLS - CONTOUR LEVEL SELECTOR',1)
      if (zincr .ne. 0.)
     $   call cpsetr ('CIS - CONTOUR INTERVAL SPECIFIER',zincr)
      if (zmax .ne. 0. .or. zmin .ne. 0.) then
         call cpsetr ('CMN - CONTOUR LEVEL MINIMUM',zmin)
         call cpsetr ('CMX - CONTOUR LEVEL MAXIMUM',zmax)
      endif
      call cpseti ('LIS - LABEL INTERVAL SPECIFIER',1)
c
c initialize the drawing of the contour plot
      call cprect (zmat,iz,iz,jz,rwrk,nwrk,iwrk,nwrk)
c
c now force conpack to pick levels so they can be modified
      call cppkcl (zmat,rwrk,iwrk)
c
c have settings for contour levels - modify negative contour levels
c from solid to dashed.
c
c get number of contour lines
      call cpgeti ('NCL - CONTOUR LEVEL SELECTION FLAG',nclv)
c
c need to modify 'CLD' from solid to dashed.  CLD is a character ARRAY,
c so first need to set the parameter array index for neg. contours -
c positive contours will use the default (solid)
      do 104 iclv=1,nclv
c set parameter array index to contour level iclv
      call cpseti ('PAI - PARAMETER ARRAY INDEX',iclv)
      call cpgetr ('CLV - CONTOUR LEVEL VALUE',clev)
      if (clev .lt. 0.)  call cpsetc 
     $  ('CLD - CONTOUR LEVEL DASH PATTERN','$$''$$''$$''$$''$$''')
c
c  set up color bar and labels
      if(ifcolor)then
	ilev=clev
        write(llbs(iclv),'(i4,2x)')ilev
	write(llbs(iclv),'(f6.0)')clev
	if(zincr.lt.1.)write(llbs(iclv),'(f6.1)')clev
	if(zincr.lt..1)write(llbs(iclv),'(f6.2)')clev
        lind(iclv)=iclv+1
      endif
 104    continue
       if(ifcolor) lnclv=nclv-1
c
c modify contour interval table which controls which contours will be
c labelled if conpack picks contour interval
      if (zincr .eq. 0.) go to 200
      do 103 i=1,10
c set parameter array index to i'th table element
      call cpseti ('PAI - PARAMETER ARRAY INDEX',i)
      call cpgetr ('CIT - CONTOUR INTERVAL TABLE',c)
      if (c .gt. 0.)  call cpseti ('LIT - LABEL INTERVAL TABLE',2)
 103  continue
c
c draw contour lines and labels
 200  continue
      if(.not.ifcolor)then
	 call cpcldr (zmat,rwrk,iwrk)
         call cplbdr (zmat,rwrk,iwrk)
         call cpgeti ('RWU',krwu)
         call cpgeti ('IWU',iwu)
         write (6,*) '0rwu=',krwu,'    iwu=',iwu
      endif
      if(ncolor.eq.3.or.ncolor.eq.9)then
        call gsclip(1)
        call cpcldr (zmat,rwrk,iwrk)
        call gsclip(0)
      endif
c
c *****************************************************************
c ERROR STIPPLING
      if (emat(1) .gt. 10.e19) go to 300
c
c initialize the area map
      call arinam(iama,niam)
c 
c set contour levels
c set cls=1 to have conpack pick levels from zmin to zmax by zincr
      call cpseti ('CLS - CONTOUR LEVEL SELECTOR',1)
      call cpsetr ('CIS - CONTOUR INTERVAL SPECIFIER',eincr)
      call cpsetr ('CMN - CONTOUR LEVEL MINIMUM',0.)
      call cpsetr ('CMX - CONTOUR LEVEL MAXIMUM',emax)
c
c initialize the drawing of the contour plot
      call cprect (emat,iz,iz,jz,rwrk,nwrk,iwrk,nwrk)
c
c now force conpack to pick levels so they can be modified
      call cppkcl (emat,rwrk,iwrk)
c
c have settings for contour levels - divide map into area above and below eincr
c
c get number of contour lines
      call cpgeti ('NCL - CONTOUR LEVEL SELECTION FLAG',nclv)
c
      do 102 iclv=1,nclv
c set parameter array index to contour level iclv
      call cpseti ('PAI - PARAMETER ARRAY INDEX',iclv)
      call cpgetr ('CLV - CONTOUR LEVEL VALUE',clev)
      if (clev .lt. (eincr-.01) .or. clev .gt. (eincr+.01))  then
         call cpseti ('AIA - AREA IDENTIFIER ABOVE LINE',0)
         call cpseti ('AIB - AREA IDENTIFIER BELOW LINE',0)
                                                   else
         call cpseti ('AIA - AREA IDENTIFIER ABOVE LINE',2)
         call cpseti ('AIB - AREA IDENTIFIER BELOW LINE',1)         
      endif
 102  continue
c
c draw contour lines
c     call cpcldr (emat,rwrk,iwrk)
c
c
c CPCLAM will put contour lines onto area map
c   and then ARSCAM will draw it
      call cpclam (emat,rwrk,iwrk,iama)
c
c scan the area map.  The routine SHADER will be called to shade the
c area above the 4.5 contour line
c
      call arscam (iama,xcs,ycs,nxcs,iaia,igia,niai,shader)
c
 300  continue
c
c  now draw the color bar
      write(6,*)' before color bar....nobar = ',nobar
      if(nobar)go to 305
      if(ifcolor)then
        call set(0.,1.,0.,1.,0.,1.,0.,1.,1)
        call LBSETR('WFL - width of fill lines',2.)
        call LBSETR('WLB - width of LABEL lines',2.)
        call lbseti('CBL - color of box lines',1)
	if(bar_reg)then
        call lblbar(0,0.15,.85,0.05,.13,lnclv+1,1.,.5,lind,
     *  0,llbs,lnclv,1)
	else
	call lblbar(0,0.05,.95,0.18,.26,lnclv+1,1.,.5,lind,
     *  0,llbs,lnclv,1)
        endif
      else
	call set (0.,1.,0.,1.,0.,1.,0.,1.,1)
	call plchhq (.35,.08,fin(1:80),8.,0,0)
	call plchhq (.35,.05,info(1:80),8.,0,0)
      endif

 305  call mapiq
      call frame
      return
      end
      SUBROUTINE COLRAMK(XCS,YCS,NCS,IAI,IAG,NAI)
c
      dimension xcs(*),ycs(*),iai(*),iag(*)
c
      iai1=-1
      do 101 i=1,nai
      if(iag(i).eq.1)iai1=iai(i)
  101 continue
      if(iai1.gt.0)then
c positive area index
	 if(mapaci(iai1).ne.1)then
            call gsfaci(1)
c fill area color index set to 0 (background) if color index.ne.1 (foreground)
	    call gfa(ncs-1,xcs,ycs)
c routine GFA draws a filled area
c draw on map
         endif
      endif
      return
      end
      SUBROUTINE COLRLN (XCS,YCS,NCS,IAI,IAG,NAI)
c
      dimension xcs(*),ycs(*),iai(*),iag(*)
c
c for each line segment, get a set of points, two group identifiers and
c their associated area identifiers.  If both of the area identifiers
c are zero or negative, the segment is not drawn; otherwise use MAPACI 
c to see if the segment is over water, and, if so, draw the segment.
c
      call gsln (3)
c gsln=2, sets a dashed line; gsln=3, sets a dotted line; (1=solid)
      if (iai(1) .ge. 0 .and. iai(2) .ge. 0) then
c left and right area identifier =>0 implies inside area; <0 =>outside
         itm=max0(iai(1),iai(2))
         if (mapaci (itm) .eq. 1) call gpl (ncs,xcs,ycs)
c if the color index for this polyline segment =1(foreground), draw it
c GPL draws polylines
      endif
      call gsln(1)
      return
      end
      SUBROUTINE SHADER (XCS,YCS,NCS,IAI,IAG,NAI)
c
c this version of SHADER shades the polygon whose edge is defined by the
c points ((xcs(i),ycs(i),i=1,ncs), if and only if, relative to edge 
c group 3,its area identifier is a 2.  SOFTFILL is used to do the shading.
c
      dimension xcs(*),ycs(*),iai(*),iag(*)
      dimension dst(10000),ind(10000)
c
c turn off shading
      ish=0
c
c if the area identifier for group 3 is a 2, turn on shading
      do 101 i=1,nai
         if (iai(i) .eq. 2) ish=1
         write (6,*) iag(i),iai(i)
c         if (iag(i) .eq. 3 .and. iai(i) .eq. 2) ish=1
 101  continue
c
c if shading is turned on, shade the area
      if (ish .ne. 0) then
         call sfsetr ('SP - SPACING OF FILL LINES',.003)
         call sfseti ('AN - ANGLE OF FILL LINES',0)
         call sfseti ('DO - DOT-FILL FLAG',1)
         call sfwrld (xcs,ycs,ncs-1,dst,10000,ind,10000)
      endif
      return
      end
      SUBROUTINE MAPLBM (STRLAT,ENDLAT,STRLON,ENDLON,STEP,steplat,nl)
c 
      character chrs*3,chlb*3,add*2,label*20
      write(6,*)' nl = ', nl
c
c label bottom meridians
      do 10 rlon=strlon,endlon+.05,step
      xlon=rlon
      if (rlon .lt. 0) xlon=rlon+360.
      if (rlon .gt. 360) xlon=rlon-360.
      nchs=0
      i=xlon
      write (chrs,'(i3)') i
      if (i .ge. 100) then
         nchs=nchs+1
         chlb(nchs:nchs)=chrs(1:1)
      endif
      if (i .ge. 10) then
         nchs=nchs+1
         chlb(nchs:nchs)=chrs(2:2)
      endif
      nchs=nchs+1
      chlb(nchs:nchs)=chrs(3:3)
      label=chlb(1:nchs)//':S:0:N:E'
      nchs=nchs+8
c
      call maptrn (strlat,rlon,u,v)
      call line (u,v,u,v-3.)
 10   if(nl.ge.2)call plchhq (u,v-7.,label(1:nchs),8.,0.,0.)
c
c label top meridians
      tstep=step
      x=strlon
      if(x.gt.180.)x=x-360
      do 20 rlon=strlon,endlon+.05,step
      if (x .gt. 179. .and. tstep .gt. 0.) tstep=-tstep
      ix=abs(x)
      write (chrs,'(i3)') ix
      nchs=0
      if (ix .ge. 100) then
         nchs=nchs+1
         chlb(nchs:nchs)=chrs(1:1)
      endif
      if (ix .ge. 10) then
         nchs=nchs+1
         chlb(nchs:nchs)=chrs(2:2)
      endif
      nchs=nchs+1
      chlb(nchs:nchs)=chrs(3:3)
c add superscript and e(ast) or w(est)
      add='E'
      if (x .lt. 0. .or. tstep .lt. 0.) add='W'
ck
      if(rlon.gt.360.)add='E'
ck
      label=chlb(1:nchs)//':S:0:N:'//add(1:1)
      nchs=nchs+8
c
      call maptrn (endlat,rlon,u,v)
      call line (u,v,u,v+3.)
      if(nl.eq.4)call plchhq (u,v+7.,label(1:nchs),8.,0.,0.)
 20   x=x+tstep
c
c label left and right zonals
      tstep=step
      if (strlat .lt. 0.) tstep=-step
ck      x=abs(strlat)
ck      do 30 rlat=strlat,endlat+.05,step
      do 30 rlat=-60.,60.,60.
      x=abs(rlat)
      if (rlat .ge. 0.) tstep=step
      ix=abs(x)
      write (chrs,'(i2)') ix
      nchs=0
      if (ix .ge. 10) then
         nchs=nchs+1
         chlb(nchs:nchs)=chrs(1:1)
      endif
      nchs=nchs+1
      chlb(nchs:nchs)=chrs(2:2)
c add superscripts and n(orth) or s(outh)
      add='S '
      if (rlat .ge. 0.) add='N '
      label='  '//chlb(1:nchs)//':S:0:N:'//add
      nchs=nchs+11
c
      call maptrn (rlat,strlon,u,v)
      call line (u,v,u-3.,v)
      if(nl.eq.1.or.nl.ge.3)
     *call plchhq (u,v,label(3:nchs+1),8.,0.,+1.)
      call maptrn (rlat,endlon,u,v)
      call line (u,v,u+3.,v)
      if(nl.eq.1.or.nl.eq.4)
     *call plchhq (u,v,label(1:nchs),8.,0.,-1.)
 30   continue
c
      return
      end
cccccccccccccccccccccccccccccccccccccccccccccc

      subroutine colorz(xcs,ycs,ncs,iai,iag,nai)
c
c  this version of shader shades the plygon whose edge is defined by
c  the points ((xcs(i),ycs(i),i=1,ncs) if and only, relative to edge
c  group 3, its area identifier is a 1.  The package softfill is used
c  to do the shading.
c
      dimension xcs(*),ycs(*),iai(*),iag(*)
c
c  polygon will be filled
c
      ifll=1
c
c  if any of the area identifiers is negative, don't fill polygon
      do 101 i=1,nai
        if(iai(i).lt.0)ifll=0
  101 continue
c
c  otherwise, fill the polygon in the color implied by its area
c  identifier relative to edge group 3 (the contour line group)
c
      if(ifll.ne.0)then
        ifll=0
      do 102 i=1,nai
         if(iag(i).eq.3) ifll=iai(i)+1 
  102 continue
c
c   color the area
c
      if (ifll.gt.0.and.ifll.le.19)then
        call gsfaci(ifll)
        call gfa(ncs-1,xcs,ycs)
                                   endif
                   endif
c
c   done
c
      return
      end
c
      subroutine dfclrs1
c  good positive run:
c  define a set of RGB color triples for colors 0 through 15
c
      dimension rgbv(3,19)
c
c  define the RGB color triple needed below.
c
c  white is index 0. black is index  1
c  This is for white background, black writing.
      data rgbv/  0.00 , 0.00 , 0.00 ,
c  blue, light green
     +            0.00 , 0.00 , 1.00 ,
     +            0.00 , 0.50 , 1.00 ,
     +            0.00 , 1.00 , 1.00 ,
     +            0.00 , 1.00 , 0.00 ,
     +            0.33 , 1.00 , 0.00 ,
     +            0.66 , 1.00 , 0.00 ,
c   yellow to red
     +            1.00 , 1.00 , 0.00 ,
c    +            1.00 , 0.88 , 0.00 ,
     +            1.00 , 0.75 , 0.00 ,
c    +            1.00 , 0.63 , 0.00 ,
     +            1.00 , 0.50 , 0.00 ,
c    +            1.00 , 0.38 , 0.00 ,
     +            1.00 , 0.25 , 0.00 ,
     +            1.00 , 0.0  , 0.00 ,
     +            0.75 , 0.0  , 0.00 ,
     +            0.63 , 0.00 , 0.00 ,
     +            0.50 , 0.00 , 0.00 ,
     +            0.38 , 0.00 , 0.00 ,
     +            0.25 , 0.00 , 0.00 ,
     +            0.00 , 0.00 , 0.00 ,
     +            0.00 , 0.00 , 0.00 /
c
c    white
       call gscr(1,0,1.0,1.0,1.0)
c
       do 101 i=1,19
       call gscr(1,i,rgbv(1,i),rgbv(2,i),rgbv(3,i))
  101 continue
      return
      end

      subroutine dfclrs2
c  good general run
c
c  define a set of RGB color triples for colors 0 through 15
c
      dimension rgbv(3,20)
c
c  define the RGB color triple needed below.
c
c  white is index 0. black is index  1
c  This is for white background, black writing.
      data rgbv/  0.00 , 0.00 , 0.00 ,
c  blue, green, yellow
     +            0.00 , 0.00 , 1.00 ,
     +            0.00 , 0.25 , 1.00 ,
     +            0.00 , 0.50 , 1.00 ,
     +            0.00 , 0.75 , 1.00 ,
     +            0.00 , 1.00 , 1.00 ,
c 2 
     +            0.00 , 1.00 , 0.66 ,
c
     +            0.00 , 1.00 , 0.00 ,
     +            0.33 , 1.00 , 0.00 ,
     +            0.66 , 1.00 , 0.00 ,
c  yellow
     +            1.00 , 1.00 , 0.00 ,
c   yellow to red 
     +            1.00 , 0.75 , 0.00 ,
     +            1.00 , 0.50 , 0.00 ,
c    +            1.00 , 0.35 , 0.00 ,
     +            1.00 , 0.00 , 0.00 ,
c   red to  blue
     +            0.75 , 0.00 , 0.00 ,
     +            0.63 , 0.00 , 0.00 ,
     +            0.50 , 0.00 , 0.00 ,
     +            0.38 , 0.00 , 0.00 ,
     +            0.25 , 0.00 , 0.00 ,
     +            0.00 , 0.00 , 0.00 /
c
c    white
       call gscr(1,0,1.0,1.0,1.0)
c
       do 101 i=1,20
       call gscr(1,i,rgbv(1,i),rgbv(2,i),rgbv(3,i))
  101 continue
      return
      end

      subroutine dfclrs3_dark
c  good general run
c
c  define a set of RGB color triples for colors 0 through 15
c
c
c  define the RGB color triple needed below.
c
c  white is index 0. black is index  1
c  This is for white background, black writing.
c
c    white
       call gscr(1,0,1.0,1.0,1.0)
       call gscr(1,1,0.0,0.0,0.0)
c
       del=1./11
       do 101 i=1,19
       rgbv=1-(i-1)*del
       if(rgbv.lt.0.)rgbv=0.
       if(rgbv.gt.1.)rgbv=1.
       call gscr(1,i+1,rgbv,rgbv,rgbv)
  101 continue
      return
      end

      subroutine dfclrs3_1
c  good general run
c
c  define a set of RGB color triples for colors 0 through 15
c
c
c  define the RGB color triple needed below.
c
c  white is index 0. black is index  1
c  This is for white background, black writing.
c
c    white
       call gscr(1,0,1.0,1.0,1.0)
       call gscr(1,1,0.0,0.0,0.0)
c
       del=1./18
       do 101 i=1,19
       rgbv=1-(i-1)*del
       if(rgbv.lt.0.)rgbv=0.
       if(rgbv.gt.1.)rgbv=1.
       call gscr(1,i+1,rgbv,rgbv,rgbv)
  101 continue
      return
      end

      subroutine dfclrs3a
c  good general run
c
c  define a set of RGB color triples for colors 0 through 15
c
c
c  define the RGB color triple needed below.
c
c  white is index 0. black is index  1
c  This is for white background, black writing.
c
c    white
       call gscr(1,0,1.0,1.0,1.0)
       call gscr(1,1,0.0,0.0,0.0)
c
       del=1./8
       k=1
       ic=0
       do 101 i=1,19,2
       rgbv=1-(k-1)*del
       k=k+1
       if(rgbv.lt.0.)rgbv=0.
       if(rgbv.gt.1.)rgbv=1.
       ic=ic+1
       call gscr(1,ic+1,rgbv,rgbv,rgbv)
       ic=ic+1
       if(i.eq.3)ic=ic-1
       if(i.eq.5)ic=ic-1
       call gscr(1,ic+1,rgbv,rgbv,rgbv)
  101 continue
      return
      end

      subroutine dfclrs3
c  good general run
c
c  define a set of RGB color triples for colors 0 through 15
c
c
c  define the RGB color triple needed below.
c
c  white is index 0. black is index  1
c  This is for white background, black writing.
c
c    white
       call gscr(1,0,1.0,1.0,1.0)
       call gscr(1,1,0.0,0.0,0.0)
c
       del=1./9
       k=1
       do 101 i=1,19,2
       rgbv=1-(k-1)*del
       k=k+1
       if(rgbv.lt.0.)rgbv=0.
       if(rgbv.gt.1.)rgbv=1.
       call gscr(1,i+1,rgbv,rgbv,rgbv)
       call gscr(1,i+2,rgbv,rgbv,rgbv)
  101 continue
      return
      end

      subroutine dfclrs4
c  good positive run:
c  define a set of RGB color triples for colors 0 through 15
c
      dimension rgbv(3,19)
c
c  define the RGB color triple needed below.
c
c  white is index 0. black is index  1
c  This is for white background, black writing.
      data rgbv/  0.00 , 0.00 , 0.00 ,
     +            0.50 , 0.00 , 0.00 ,   !dark red
     +            0.75 , 0.00 , 0.00 ,   !mid red
     +            1.00 , 0.00 , 0.00 ,   !red
     +            1.00 , 0.63 , 0.00 ,   !dark orange
     +            1.00 , 0.88 , 0.00 ,   !orange
     +            1.00 , 1.00 , 0.00 ,   !yellow
     +            0.88 , 1.00 , 0.00 ,   !light green
     +            0.00 , 0.88 , 0.00 ,   !dark green
     +            0.00 , 1.00 , 0.00 ,   !green

     +            0.00 , 1.00 , 1.00 ,   !aqua
     +            0.00 , 0.60 , 1.00 ,   !mid blue
     +            0.00 , 0.00 , 1.00 ,   !blue      *
     +            0.50 , 0.00 , 1.00 ,   !purple
     +            0.75 , 0.38 , 1.00 ,   !light purple
     +            0.85 , 0.31 , 1.00 ,   !magenta3
     +            1.00 , 0.25 , 1.00 ,   !magenta
     +            0.30 , 0.00 , 0.80 ,   !dark magenta
     +            0.20 , 0.00 , 0.40 /   !blackish
c
c    white
       call gscr(1,0,1.0,1.0,1.0)
c
       do 101 i=1,19
       call gscr(1,i,rgbv(1,i),rgbv(2,i),rgbv(3,i))
  101 continue
      return
      end

      subroutine dfclrs2a
c  good anomaly general run (2 center positions are same color)
c
c  define a set of RGB color triples for colors 0 through 15
c
      dimension rgbv(3,20)
c
c  define the RGB color triple needed below.
c
c  white is index 0. black is index  1
c  This is for white background, black writing.
      data rgbv/  0.00 , 0.00 , 0.00 ,
c  blue, green, yellow
     +            0.00 , 0.00 , 1.00 ,
     +            0.00 , 0.25 , 1.00 ,
     +            0.00 , 0.50 , 1.00 ,
     +            0.00 , 0.75 , 1.00 ,
     +            0.00 , 1.00 , 1.00 ,
c 2 
     +            0.00 , 1.00 , 0.66 ,
c
     +            0.00 , 1.00 , 0.00 ,
     +            0.33 , 1.00 , 0.00 ,
c    +            0.66 , 1.00 , 0.00 ,
c  yellow
     +            1.00 , 1.00 , 0.63 ,
     +            1.00 , 1.00 , 0.63 ,
c   yellow to red 
     +            1.00 , 0.75 , 0.00 ,
     +            1.00 , 0.50 , 0.00 ,
c    +            1.00 , 0.35 , 0.00 ,
     +            1.00 , 0.00 , 0.00 ,
c   red to  blue
     +            0.75 , 0.00 , 0.00 ,
     +            0.63 , 0.00 , 0.00 ,
     +            0.50 , 0.00 , 0.00 ,
     +            0.38 , 0.00 , 0.00 ,
     +            0.25 , 0.00 , 0.00 ,
     +            0.00 , 0.00 , 0.00 /
c
c    white
       call gscr(1,0,1.0,1.0,1.0)
c
       do 101 i=1,20
       call gscr(1,i,rgbv(1,i),rgbv(2,i),rgbv(3,i))
  101 continue
      return
      end


      subroutine dfclrs2b
c  good anomaly general run (2 center positions are same color)
c
c  define a set of RGB color triples for colors 0 through 15
c
      dimension rgbv(3,20)
c
c  define the RGB color triple needed below.
c
c  white is index 0. black is index  1
c  This is for white background, black writing.
      data rgbv/  0.00 , 0.00 , 0.00 ,
c  blue, green, yellow
     +            0.   , 0.   ,  .75 ,
     +            0.00 , 0.   , 1.   ,
     +            0.00 , 0.23 , 1.00 ,
     +            0.00 , 0.49 , 1.00 ,
     +            0.00 , 0.68 , 1.00 ,
     +            0.00 , 0.87 , 1.00 ,
     +            0.00 , 0.84 , 0.50 ,
     +            0.00 , 1.00 , 0.00 ,
c grey for middle
     +             .93 ,  .93  , .93  ,
     +             .93 ,  .93 ,  .93 ,
c  yellow
     +            1.00 , 1.00 , 0.00 ,
c   yellow to red 
     +            1.00 , 0.75 , 0.00 ,
     +            1.00 , 0.50 , 0.00 ,
     +            1.00 , 0.35 , 0.00 ,
     +            1.00 , 0.00 , 0.00 ,
c   red to  blue
     +            0.75 , 0.00 , 0.00 ,
     +            0.63 , 0.00 , 0.00 ,
     +            0.50 , 0.00 , 0.00 ,
     +            0.38 , 0.00 , 0.00 /
c
c    white
       call gscr(1,0,1.0,1.0,1.0)
c
       do 101 i=1,20
       call gscr(1,i,rgbv(1,i),rgbv(2,i),rgbv(3,i))
  101 continue
      return
      end

      subroutine dfclrs2c
c  good anomaly general run (2 center positions are same color)
c
c  define a set of RGB color triples for colors 0 through 15
c
      dimension rgbv(3,20)
c
c  define the RGB color triple needed below.
c
c  white is index 0. black is index  1
c  This is for white background, black writing.
      data rgbv/  0.00 , 0.00 , 0.00 ,
c  blue, green, yellow
     +            0.   , 0.   ,  .75 ,
     +            0.00 , 0.   , 1.   ,
     +            0.00 , 0.23 , 1.00 ,
     +            0.00 , 0.49 , 1.00 ,
     +            0.00 , 0.68 , 1.00 ,
     +            0.00 , 0.87 , 1.00 ,
     +            0.00 , 0.84 , 0.50 ,
     +            0.00 ,  .92 ,  .25 ,
     +            0.00 , 1.00 , 0.00 ,
c  yellow
     +            1.00 , 1.00 , 0.00 ,
     +            1.   ,  .87 , 0.   ,
c   yellow to red 
     +            1.00 , 0.75 , 0.00 ,
     +            1.00 , 0.50 , 0.00 ,
     +            1.00 , 0.35 , 0.00 ,
     +            1.00 , 0.00 , 0.00 ,
c   red to  blue
     +            0.75 , 0.00 , 0.00 ,
     +            0.63 , 0.00 , 0.00 ,
     +            0.50 , 0.00 , 0.00 ,
     +            0.38 , 0.00 , 0.00 /
c
c    white
       call gscr(1,0,1.0,1.0,1.0)
c
       do 101 i=1,20
       call gscr(1,i,rgbv(1,i),rgbv(2,i),rgbv(3,i))
  101 continue
      return
      end

      subroutine dfclrs1a
c  good positive run:
c  define a set of RGB color triples for colors 0 through 15
c
      dimension rgbv(3,19)
c
c  define the RGB color triple needed below.
c
c  white is index 0. black is index  1
c  This is for white background, black writing.
      data rgbv/  0.00 , 0.00 , 0.00 ,
c  blue, light green
     +            0.00 , 0.00 , 1.00 ,
     +            0.00 , 0.50 , 1.00 ,
     +            0.00 , 1.00 , 1.00 ,
     +            0.00 , 1.00 , 0.00 ,
     +            0.33 , 1.00 , 0.00 ,
ck     +            0.66 , 1.00 , 0.00 ,
c   yellow to red
     +            1.00 , 1.00 , 0.00 ,
c    +            1.00 , 0.88 , 0.00 ,
     +            1.00 , 0.75 , 0.00 ,
c    +            1.00 , 0.63 , 0.00 ,
     +            1.00 , 0.50 , 0.00 ,
c    +            1.00 , 0.38 , 0.00 ,
ck     +            1.00 , 0.25 , 0.00 ,
     +            1.00 , 0.0  , 0.00 ,
     +            0.75 , 0.0  , 0.00 ,
     +            0.63 , 0.00 , 0.00 ,
     +            0.50 , 0.00 , 0.00 ,
     +            0.38 , 0.00 , 0.00 ,
     +            0.25 , 0.00 , 0.00 ,
     +            0.00 , 0.00 , 0.00 ,
     +            0.00 , 0.00 , 0.00 ,
     +            0.00 , 0.00 , 0.00 ,
     +            0.00 , 0.00 , 0.00 /
c
c    white
       call gscr(1,0,1.0,1.0,1.0)
c
       do 101 i=1,19
       call gscr(1,i,rgbv(1,i),rgbv(2,i),rgbv(3,i))
  101 continue
      return
      end
      subroutine dfclrs6
c
c  define a set of RGB color triples for colors 0 through 15
c
      dimension rgbv(3,19)
c
c  define the RGB color triple needed below.
c
c  white is index 0. black is index  1
c  This is for white background, black writing.
      data rgbv/  0.00 , 0.00 , 0.00 ,
     +            1.00 , 0.00 , 1.00 ,
     +            0.80 , 0.00 , 1.00 ,
c purple
     +            0.60 , 0.00 , 1.00 ,
     +            0.38 , 0.00 , 1.00 ,
c blue
     +            0.00 , 0.00 , 1.00 ,
     +            0.00 , 0.50 , 1.00 ,
c
     +            0.00 , 0.88 , 1.00 ,
c green
c    +            0.00 , 1.00 , 0.25 ,
     +            0.50 , 1.00 , 0.00 ,
c yellow
     +            1.00 , 1.00 , 0.00 ,
c
     +            1.00 , 0.75 , 0.00 ,
     +            1.00 , 0.50 , 0.00 ,
     +            1.00 , 0.25 , 0.00 ,
c red
     +            1.00 , 0.00 , 0.00 ,
     +            0.90 , 0.00 , 0.00 ,
     +            0.75 , 0.00 , 0.00 ,
     +            0.63 , 0.00 , 0.00 ,
c    +            0.45 , 0.00 , 0.00 ,
     +            0.25 , 0.00 , 0.00 ,
     +            0.00 , 0.00 , 0.00 /
c
c    white
      call gscr(1,0,1.0,1.0,1.0)
c
      do 101 i=1,19
        call gscr(1,i,rgbv(1,i),rgbv(2,i),rgbv(3,i))
  101 continue
      return
      end

