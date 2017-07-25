      program mesh2d
c
c--------------------------------------------------------------------
c
c     This program performs tha automatic subdivision of arbitrary 2d
c     domains into an unstructured mesh of triangular straight sided
c     elements.
c     This program was developed by Jaime Peraire and Ken Morgan from
c     the University College of Swansea, U.K. and was completed in 
c     October 1985.
c
c--------------------------------------------------------------------
c    
      parameter (mxele=1500000,mxpoi=700000,mxbou=20000,mxseg=500)
      parameter (mxhlp=7*mxele+8*mxpoi)
      parameter (mxhli=mxele+mxseg*(mxbou+2)*2+20*mxseg)
      parameter (nmax = 50000)
      character filename*200,st*80,convrt*5,yn*1
      dimension coor (2,mxpoi) , strec (4,mxpoi) , coorg(2,mxpoi)  
      dimension iel  (3,mxele) , intmeg(3,mxele) , ieleg(3,mxele) 
      dimension iside(mxhli)   , icone(3*mxele)  , lcoor(mxpoi)  
      dimension lcore(mxpoi)   , lcoid(mxpoi)    , lboud(mxpoi)
      dimension lwher(mxpoi)   , lhowm(mxpoi)    , lposi(mxpoi)
      dimension unkng(4,mxpoi) , help (mxhlp)    , delta(4,mxpoi)
      dimension lep  (mxpoi)   , iseg (mxbou)    , ibsid(4,mxbou)
      dimension unkno(4,mxpoi) , lbk(nmax)       , lplay(mxpoi)
      dimension lelay( mxele ) , lafter(mxpoi)
      common /plt/ iplot
c     write(*,49)
c     read(*,'(a)')yn
      yn='n'
      iplot=0
      if(yn.eq.'n'.or.yn.eq.'N') iplot=1
      lbk(1) = nmax
      i1=1
      i2=i1+mxseg
      i3=i2+mxseg
      i4=i3+mxseg
      i5=i4+mxseg
      i6=i5+mxseg
      i7=i6+2*mxbou
      i8=i7+mxseg*mxbou
      i9=i8+mxseg*mxbou
      i10=i9+2*mxseg
      i11=i10+2*mxseg
      j1=1
      j2=j1+2*mxpoi
      j3=j2+2*mxpoi
      j4=j3+mxpoi
      j5=j4+mxpoi
      j6=j5+mxpoi
777   write (*, 100)
      read  (*, '(a)') filename
      open  (8, file=filename, status = 'old', err=777)      
778   write (*, 200)  
      read  (*, '(a)') filename
      open  (11, file=filename, status = 'old', err=778)
780   write (*, 400)
      read  (*, '(a)') filename
      open  (9, file=filename, status = 'unknown', err=780)
779   write (*,300)
      read  (*, * ,err=779)io
879   write (*,600)
      read  (*, * ,err=879)iquad
      if(iquad.eq.2) then
       open  (19, file='p_'//filename, status = 'unknown')
       open  (96, file='L_'//filename, status = 'unknown')
      end if
      nelem=0
c
c *** read the background   grid data
c
      call  input(neleg ,npoig ,ieleg ,coorg ,delta ,unkng ,io   ,
     *            strec ,help  ,coor  ,mxpoi ,mxele ,lhowm ,lwher,
     *            icone ,iside ,lbk,  lbk)
      iconv=0
c     read(11,*,err=5,end=5) iconv
   5  if(iconv.eq.0)
     &call convex(neleg ,npoig ,ieleg ,coorg ,iside(i7),iside(i8) ,
     *                  lcoor ,iside(i1) )
      call   plot(1     ,x     ,y     ,st    ,xmin  ,xmax  ,ymin ,
     *           ymax  ) 
c
c *** fill in iside
c
      call   side(neleg ,npoig ,nsidg ,ieleg ,iside, lwher ,lhowm , 
     *          icone , 0    ) 
c
c *** fill in intmeg
c
      call  elecon(nsidg ,ieleg ,iside ,intmeg )
c
c *** expand the background mesh to include all the domain
c
      call expand(neleg ,npoig ,coorg ,coor  ,ieleg ,intmeg )
c
      st='  generating  boundary segments      please wait  '
      call   plot(2     ,x     ,y     ,st    ,xmin  ,xmax  ,ymin  ,
     *            ymax  )
      call rfiliv(lplay, mxpoi , 0)
      call rfiliv(lelay, mxele , 0)
      call genmsh(mbcs  ,iside(i1)    ,iside(i8)    ,npoig ,neleg , 
     *            ieleg ,coorg ,intmeg,nboun ,node  ,help(j1)     ,
     *            coor  ,delta ,iside(i2)    ,lcoor ,iside(i6)    ,
     *            lboud ,strec ,mxseg ,mxpoi ,mxbou ,
     *            xmin  ,xmax  ,ymin  ,ymax  ,iside(i4)    ,
     *            iside(i3)    ,iside(i7)    ,help(j2)     ,lposi ,
     *            lcoid        ,iel          ,nelem ,
     *            nel   ,       npl          ,iside(i9)    ,
     *            deltmi,help(j3)     ,help(j4) ,lcore ,help(j6)  ,
     *            lep   ,help(j5)     ,lwher    ,lhowm ,icone     ,
     *            mxele ,unkng ,unkno,io ,lbk ,lbk, lplay , lelay ,
     *            nlay,lafter)
      call   plot(6     ,x     ,y     ,st    ,xmin  ,xmax  ,ymin  ,
     *            ymax  )
      st='              generating  finished    '
      call   plot(2     ,x     ,y     ,st    ,xmin  ,xmax  ,ymin  ,   
     *            ymax  )
c
19999 continue
      call getopt (iwhat )
      if (iwhat.eq.1) then
c
c     find out boundary points
c
        call  bound(node  ,nelem ,lcoor ,iel   )
c
c     find out real and ideal number of conectivities
c
        call conere(node  ,nelem ,iel   ,lcore )
        call coneid(node  ,nelem ,iel   ,lcoor ,lcoid ,coor  ,strec )
c
c     fill in iside
c
        call   side(nelem ,node  ,nside ,iel   ,iside ,lwher ,lhowm ,  
     *              icone ,0     )
c 
        st='              eating 3s....please wait    '
        call   plot(2     ,x     ,y     ,st    ,xmin  ,xmax  ,ymin  , 
     *              ymax  )
        call   eat3(node  ,nelem ,nside ,iel   ,iside ,lcoor ,lposi ,
     *              lwher ,lhowm ,lcore ,lcoid ,icone ,coor  ,nel   ,
     *              npl   ,lep   ,strec ,lboud ,unkno )
        st='         swapping diagonals....please wait  '
        call   plot(2     ,x     ,y     ,st    ,xmin  ,xmax  ,ymin  ,  
     *              ymax  )
        call coneid(node  ,nelem ,iel   ,lcoor ,lcoid ,coor  ,strec )
        call swapdi(node  ,nelem ,nside ,iside ,iel   ,lcore ,lcoid , 
     *              lcoor ,coor  ,nel   ,npl   )
        nsmoo = 2
        st='             smoothing....please wait      '
        call   plot(2     ,x     ,y     ,st    ,xmin  ,xmax  ,ymin  ,  
     *              ymax  )
        call   side(nelem ,node  ,nside ,iel   ,iside ,lwher ,lhowm ,  
     *              icone ,0     )
        call smooth(2     ,3     ,nelem ,node  ,nsmoo ,lcoor ,lcore ,  
     *              iel   ,coor  ,help  ,nel   ,npl   ,lwher ,lhowm ,
     *              icone )
        st='        swapping diagonals....please wait  '
        call   plot(2     ,x     ,y     ,st    ,xmin  ,xmax  ,ymin  ,  
     *              ymax  )
        call swapdi(node  ,nelem ,nside ,iside ,iel   ,lcore ,lcoid ,   
     *              lcoor ,coor  ,nel ,npl )
        nsmoo = 2
        st='             smoothing....please wait      '
        call   plot(2     ,x     ,y     ,st    ,xmin  ,xmax  ,ymin  ,  
     *              ymax  )
        call   side(nelem ,node  ,nside ,iel   ,iside ,lwher ,lhowm ,  
     *              icone ,0     )
        call smooth(2     ,3     ,nelem ,node  ,nsmoo ,lcoor ,lcore ,  
     *              iel   ,coor  ,help  ,nel   ,npl  ,lwher ,lhowm ,
     *              icone )
        st='             cosmetics  finished    '   
        call   plot(2     ,x     ,y     ,st    ,xmin  ,xmax  ,ymin  ,  
     *              ymax  )
        goto 19999
      else if (iwhat.eq.2) then
        st='             checking mesh....please wait'
        call   plot(2     ,x     ,y     ,st    ,xmin  ,xmax  ,ymin  ,  
     *              ymax  )
        call areach(node  ,nelem ,iel   ,coor  )
        goto 19999
      else if (iwhat.eq.3) then
        st='              saving     data   '
        call   plot(2     ,x     ,y     ,st    ,xmin  ,xmax  ,ymin  ,  
     *              ymax  )
        call   side(nelem ,node  ,nside ,iel   ,iside ,lwher ,lhowm ,  
     *               icone ,0    )
        call  elecon(nside ,iel ,iside ,intmeg )
        call oustar(node  ,nelem ,coor  ,iel   ,lboud ,iside ,nside ,
     *              iseg  ,ibsid ,unkno ,strec ,intmeg,lplay ,lelay ,
     *              nlay  ,iquad ,lafter)
        st='              saving   fininshed    '
        call   plot(2     ,x     ,y     ,st    ,xmin  ,xmax  ,ymin  ,  
     *              ymax  )
        goto 19999
      else if(iwhat.eq.4) then
        nne=nelem
c       if(iq.eq.1) nne=nelem-nel/2
        st='        npoin='//convrt(node)//'     nelem='//convrt(nne)
        call   plot(2     ,x     ,y     ,st    ,xmin  ,xmax  ,ymin  ,  
     *              ymax  )
        call   plot(3     ,x     ,y     ,st    ,xmin  ,xmax  ,ymin  ,  
     *              ymax  )
        call gmplot(node  ,coor  ,nelem ,iel   ,iq    ,nel   )
        call   plot(6     ,x     ,y     ,st    ,xmin  ,xmax  ,ymin  ,  
     *              ymax  )
        goto 19999
      else if(iwhat.eq.5) then
        call   plot(7     ,x     ,y     ,st    ,xmin  ,xmax  ,ymin  ,  
     *              ymax  )
        stop 
      else
        goto 19999
      end if
1     format(a80)
 49   format(' Do you want to plot (Y/N)  ?', $)
100   format(' input  geometry  filename  :', $)
200   format(' input background filename  :', $)
400   format(' input   result   filename  :', $)
500   format(' input   unknown  filename  :', $)
300   format(' remeshing ? (1 yes) ,or', 
     *       ' generating initial mesh ? (0 yes)  : ',$)
600   format(' Triangles in viscous regions ? (1) ,or', 
     *       ' quadrilaterals in viscous regions? (2)  : ',$)
      end
c
c--------------------------------------------------------------------
c
c     subroutine to read the background grid
c
c--------------------------------------------------------------------
c
      subroutine input(neleg ,npoig ,ieleg ,coorg ,delta ,unkng ,      
     *                 io    ,help0 ,help1 ,help2 ,mxpoi ,mxele ,
     *                 lhowm ,lwher ,icone ,iside ,lbk, abk )
      dimension coorg(2,*),delta(4,*),unkng(4,*),help0(*)
      dimension help1(*)  ,help2(*)  ,ieleg(3,*),lhowm(*)
      dimension lwher(*)  ,icone(*)  ,iside(4,*)
      dimension lbk(*), abk(*)
      character *80 text
      if(io.eq.0)then
      read(11,'(a)',err=10) text
      read(11,*,err=10) npoig,neleg
      if( npoig .gt. mxpoi ) then
        stop 'error 30 '
      end if
      if( neleg .gt. mxele ) then
        stop 'error 40 '
      end if
      do 117 in=1,npoig
        read(11,*,err=10) j,(coorg(i,j),i=1,2),(unkng(i,j),i=1,4),
     *                    d1,d2,an1,an2
      delta(1,j)=d2
      delta(2,j)=d1/d2
      an=sqrt(an1*an1+an2+an2)
      delta(3,j)=an1/an
      delta(4,j)=an2/an
117   continue
      do 116 ie=1,neleg
        read(11,*,err=10) j,(ieleg(i,j),i=1,3)
116   continue
      NMAX = lbk(1)
      read(11,*)
      read(11,*) n0,n1
      nma    = 6+5*n0+10*n1
      if(nma.gt.NMAX) stop 'increase NMAX'
      lbk(1) = n0
      lbk(2) = n1
      lbk(3) = 6
      lbk(4) = lbk(3)+ 5*n0
      lbk(5) = lbk(4)+10*n1
      read(11,*)
      do 400 ip=1,n0
      read(11,*)
      ish   = lbk(3)+(ip-1)*5
      read(11,*) abk(ish   ),abk(ish+ 1),abk(ish+ 2),
     -            abk(ish+ 3),abk(ish+ 4)
  400 continue
      read(11,*)
      do 500 ip=1,n1
      read(11,*)
      ish   = lbk(4)+(ip-1)*10
      read(11,*) abk(ish   ),abk(ish+ 1),abk(ish+ 2),
     -            abk(ish+ 3),abk(ish+ 4)
      read(11,*) abk(ish+ 5),abk(ish+ 6),abk(ish+ 7),
     -            abk(ish+ 8),abk(ish+9)
  500 continue
      else
      lbk(1) = 0
      lbk(2) = 0
      lbk(3) = 0
      read(11,*,err=10) nl   
      do 819 il = 1 , nl
      read(11,'(a)',err=10) text
 819  continue
      read(11,'(a)',err=10) text
      read(11,*,err=10)neleg, npoig
      if( npoig .gt. mxpoi ) then
        stop 'error 30 '
      end if
      if( neleg .gt. mxele ) then
        stop 'error 40 '
      end if
      read(11,'(a)',err=10) text
      do 196 ie=1,neleg
        read(11,*,err=10) j,(ieleg(i,j),i=1,3)
196   continue
      read(11,'(a)',err=10) text
      do 118 in=1,npoig
        read(11,*,err=10) j,(coorg(i,j),i=1,2)
c    *                  ,unkng(2,j),unkng(3,j),unkng(1,j),unkng(4,j),a,a
118   continue
      read(11,'(a)',err=10) text
      do 119 in=1,npoig
        read(11,*,err=10) j,(unkng(i,j),i=1,4)
119   continue
      endif
c
c *** connectivities
c **** get the value of the refinement indicators
c
      i1= 1
      i2=7*neleg+1
      i3=i2+npoig
      i4=i3+4*npoig
      i5=i4+npoig
      i6=i5+npoig
c
      if(io.eq.1) then
        call   side(neleg ,npoig ,nsidg ,ieleg ,iside, lwher ,lhowm , 
     *            icone ,0 ) 
        call ref2(npoig ,neleg ,coorg ,unkng ,ieleg, delta,
     *                       help1(i1) ,help1(i2) ,help1(i3) ,help0 ,
     *                       help2,lhowm,lwher,icone,help1(i4),
     *                       help1(i5),help1(i6) )
      end if
      call refin (neleg ,npoig ,coorg ,delta )
      return
  10  stop 'error 10'
      end
c
c--------------------------------------------------------------------
c
c     this subroutine is  called  once  for each boundary segment
c     it interpolates additional boundary points according to the
c     spacing specified by the background grid
c
c--------------------------------------------------------------------
c
      subroutine interp(npoig,neleg,npoin,coorg,ieleg,intmeg,
     *                  delta,coord,coorn,nbno,nnn,lcoor,
     *                  lbou,nbou,ibs,ilast,lboud,
     *                  mxbou,mxpoi,mxseg,nn,ln,icond,
     *                  lelch ,deltmi,ibl ,lbk,abk)
      parameter ( mxs =50000 )
      dimension coorg(2,*) , coord(2,*)  , delta(4,*)
      dimension ieleg(3,*) ,intmeg(3,*)  , lboud(*)
      dimension coorn(2,*) ,nbno(mxseg,*), ln(mxseg,*)
      dimension nn(*)      ,icond(*)
      dimension x(2,mxs),xx(2,mxs),tsp(2,mxs)
      dimension lcoor(*),nnn(*)
      dimension lbou(2,*),lelch(*)
      dimension lbk(*), abk(*)
c
c *** store the intersection points
c
      if(nbou.eq.0) then
        nbou=2
        npoin=2
        lbou(1,1)=ln(ibs,1)
        lbou(1,2)=ln(ibs,nn(ibs))
        lbou(2,1)=1
        lbou(2,2)=2
        do 108 j=1,2
          coord(j,1)=coorn(j,ln(ibs,1))
          coord(j,2)=coorn(j,ln(ibs,nn(ibs)))
108     continue
      else
        kfir=ln(ibs,1)
        ksec=ln(ibs,nn(ibs))
        ifir=1
        isec=1
        do 100 ibou=1,nbou
          if(lbou(1,ibou).eq.kfir) ifir=0
          if(lbou(1,ibou).eq.ksec) isec=0
100     continue
        if(ifir.eq.0) goto 110
        nbou=nbou+1
        npoin=npoin+1
        coord(1,npoin)=coorn(1,kfir)
        coord(2,npoin)=coorn(2,kfir)
        lbou(1,nbou)=kfir
        lbou(2,nbou)=npoin
110     if(isec.eq.0) goto 120
        nbou=nbou+1
        npoin=npoin+1
        coord(1,npoin)=coorn(1,ksec)
        coord(2,npoin)=coorn(2,ksec)
        lbou(1,nbou)=ksec
        lbou(2,nbou)=npoin
120     continue
      endif

c
c *** transfer
c
      do 200 i=1,nn(ibs)
        do 201 j=1,2
          x(j,i)=coorn(j,ln(ibs,i))
201     continue
200   continue
c
      call split(2,nn(ibs),np,x,tsp,xx,npoig,neleg,deltmi,
     *           coorg ,ieleg ,intmeg,delta ,lelch ,ilast,ibl,
     *           lbk,abk)
c
      nnn(ibs)=np
      if(np.gt.mxbou) stop 'error 80'
      nbno(ibs,1)=ln(ibs,1)
      nbno(ibs,np)=ln(ibs,nn(ibs))
c
c *** now store all points but the end ones
c
      if(np.eq.2) goto 9000
      do 9001 ip=2,np-1
        npoin=npoin+1
        if ( npoin .gt. mxpoi) stop 'error 90'
        nbno(ibs,ip)=npoin
        do 9002 j=1,2
          coord(j,npoin)=xx(j,ip)
 9002   continue
        lboud(npoin)=icond(ibs)
9001  continue
9000  continue
      return
  20  stop ' error 20'
      end
c
c-------------------------------------------------------------------
c
c     this subroutine fill the closed domain by triangular elements
c
c--------------------------------------------------------------------
c
      subroutine triang(nonf,npfrt,nqfrt,nregi,nonr,iel,
     *                    coor,nelem,node,neleg,npoig,ieleg,intmeg,
     *                    coorg,delta,strec,toler,ilast,qfr,near,
     *                    near1,howf,howf1,ncheck,mxele,mxpoi,lelch,
     *                    lep , ibl , lmb,unkng,unkno,lmp,
     *                    lbk,abk)
      dimension nregi(*),lep(*) , lmb(*) , lmp(*)
      dimension npfrt(*),nqfrt(*),unkng(4,*),unkno(4,*)
      dimension coor(2,*),strec(4,*),qfr(7,*)
      dimension iel(3,*),near(*),near1(*),howf(*)
      dimension howf1(*),ar(3),lelch(*)
      dimension ncheck(*),intmeg(3,*)
      dimension ieleg(3,*),coorg(2,*),delta(4,*)
      dimension lbk(*),abk(*),xr(2)
      common /plt/ iplot
      character *80 st
      character *5 convrt
      xtr(xq,yq,ax,ay,alph)=ax*alph*(ax*xq+ay*yq)+ay*(ay*xq-ax*yq)
      ytr(xq,yq,ax,ay,alph)=ay*alph*(ax*xq+ay*yq)-ax*(ay*xq-ax*yq)
      xba(xq,yq,ax,ay,alph)=(ax*(ax*xq+ay*yq)/alph)+ay*(ay*xq-ax*yq)
      yba(xq,yq,ax,ay,alph)=(ay*(ax*xq+ay*yq)/alph)-ax*(ay*xq-ax*yq)
      itedg = 0
      aw = 1.1
      do 19 ip = 1,mxpoi
       lmp(ip) = 0
 19   continue
        do 570 ik=1,node
          ncheck(ik)=-100
570     continue
        do 297 nl=1,nonf
        kn=nqfrt(nl)
        kn1=npfrt(nl)
        xn=coor(1,kn)
        yn=coor(2,kn)
        xn1=coor(1,kn1)
        yn1=coor(2,kn1)
        alph=min(strec(2,kn1),strec(2,kn))
c       if(ibl.eq.1) alph = 1.
        anx=0.5*(strec(3,kn1)+strec(3,kn))
        any=0.5*(strec(4,kn1)+strec(4,kn))
        anm=sqrt(anx*anx+any*any)
        anx=anx/anm
        any=any/anm
        qfr(1,nl)=anx
        qfr(2,nl)=any
        qfr(3,nl)=xba(xn,yn,anx,any,alph)
        qfr(4,nl)=yba(xn,yn,anx,any,alph)
        qfr(5,nl)=xba(xn1,yn1,anx,any,alph)
        qfr(6,nl)=yba(xn1,yn1,anx,any,alph)
        x12f=qfr(3,nl)-qfr(5,nl)
        y12f=qfr(4,nl)-qfr(6,nl)
        qfr(7,nl)=sqrt(x12f*x12f+y12f*y12f)
        lmp(kn) = 1
        lmp(kn1) = 1
 297    continue
c
c     set up ncheck values for region
c
        do 580 ik=1,nonf
          k1=npfrt(ik)
          ncheck(k1)=2
580     continue
        disw=0.0
        ibs = 1
160     continue
        nl=nonf
161     continue
        kn=nqfrt(nl)
        kn1=npfrt(nl)
        xn=coor(1,kn)
        yn=coor(2,kn)
        xn1=coor(1,kn1)
        yn1=coor(2,kn1)
        alph=min(strec(2,kn1),strec(2,kn))
c       if(ibl.eq.1) alph =1.
        anx=qfr(1,nl)
        any=qfr(2,nl)
        xnf=qfr(3,nl)
        ynf=qfr(4,nl)
        xn1f=qfr(5,nl)
        yn1f=qfr(6,nl)
        alen1f=qfr(7,nl)
        if(ibl.eq.1.and.ibs.eq.1) then
         if(lmb(kn).gt.0.and.lmb(kn1).gt.0) goto 162
        else
         if(ibl.eq.1.and.lmb(kn).le.0.and.lmb(kn1).le.0) goto 159
         if(alen1f.lt.aw*disw) goto 162
        end if
159     disw1=0.0
        call order(nl,npfrt,nqfrt,coor,npoig,neleg,coorg,
     *             ieleg,delta,strec,disw,disw1,qfr,ibl,lmb,ibs,
     *             lmp)
        goto 161
162     continue
        tole1=0.00001*alen1f
c
c      find out average spacing
c
        average=min(strec(1,kn1),strec(1,kn))
        x12=xn-xn1
        y12=yn-yn1
        alen1=sqrt(x12*x12+y12*y12)
c
c      create a new node
c
      csafe=1.0
      dside=csafe*average
      if(ibl.eq.1)then
       tolr = 2.00
      else
       tolr = 1.4
      endif
      if(dside.gt.tolr*alen1f) dside=tolr*alen1f
      if(dside.lt.0.55*alen1f) dside=0.55*alen1f
      twod2=2.0*dside*dside
      xbarf=0.5*(xnf+xn1f)
      ybarf=0.5*(ynf+yn1f)
      xdiff=xnf-xn1f
      ydiff=ynf-yn1f
      dkarg=xdiff*xdiff+ydiff*ydiff
      d12=sqrt(dkarg)
      hkarg=0.5*twod2-0.25*d12*d12
      hk=sqrt(hkarg)
      xcbf=-hk*ydiff/d12
      ycbf= hk*xdiff/d12
      xtempf=xbarf+xcbf
      ytempf=ybarf+ycbf
      xtemp=xtr(xtempf,ytempf,anx,any,alph)
      ytemp=ytr(xtempf,ytempf,anx,any,alph)
c
c     loop over possible nodes - find closest neighbours
c
      a=yn1-yn
      b=xn-xn1
      c=(xn1-xn)*yn1+(yn-yn1)*xn1
      radius=dside
      h1=0.8*radius
      inum=0
      do 110 kp=1,nonr
        ken=nregi(kp)
        if(ken.eq.kn.or.ken.eq.kn1) go to 110
        xken=coor(1,ken)
        yken=coor(2,ken)
        xkenf=xba(xken,yken,anx,any,alph)
        ykenf=yba(xken,yken,anx,any,alph)
        xdiff1=xkenf-xtempf
        ydiff1=ykenf-ytempf
        distf=sqrt(xdiff1*xdiff1+ydiff1*ydiff1)
        if(distf.gt.h1) go to 110
        if(a*xken+b*yken+c.le.0.0) goto 110

        inum=inum+1
        howf1(inum)=distf
        near1(inum)=ken
110   continue
c
c     decide which of the nodes is chosen
c
       if((lmb(kn).eq.100.or.lmb(kn1).eq.100).and.itedg.le.3)then
         inum=1
         near(inum)=0
         howf(inum)=0.0
         goto 812
        end if
       if(inum.eq.0) then
       if(ibl.eq.0.or.(ibl.eq.1.and.lmb(kn).gt.0.and.lmb(kn1).gt.0))then
         inum=1
         near(1)=0
         howf(1)=0.0
       end if
       else
c
c     order them
c
         do 599 i=1,inum
           comp=1.e+6
           do 598 j=1,inum
             if(near1(j).eq.0) goto 598
             if(howf1(j).gt.comp) goto 598
             is=j
             comp=howf1(j)
598        continue
           near(i)=near1(is)
           howf(i)=howf1(is)
           near1(is)=0
599      continue
c
c      add the new point to the list
c
        if(ibl.ne.1.or.(lmb(kn).ge.1.and.lmb(kn1).ge.1)) then
         inum=inum+1
         near(inum)=0
         howf(inum)=0.0
        end if
       endif
c
c      select ---> start by the closest
c
 812   do 601 i=1,inum
       kp=near(i)
         if(kp.eq.0) then
           xp=xtemp
           yp=ytemp
         else
           if(ibl.eq.1.and.lmb(kn).ge.1.and.lmb(kn1).ge.1.and.
     &                     lmb(kp).ge.1) goto 601
           xp=coor(1,kp)
           yp=coor(2,kp)
         endif
c
c      see if this connection is possible
c
       call possib(kn1,kn,kp,xn1,yn1,xn,yn,xp,yp,nonf,
     *               nonr,npfrt,nqfrt,nregi,coor,iyon)

       if(iyon.ne.0) goto 603
  601  continue
c
c      we are in trouble !!!!!
c      find the 30 existing nodes that give maximum angle
c
       inum=0
       ang1=0.0
       do 210 kp=1,nonr
       ken=nregi(kp)
       if(ken.eq.kn.or.ken.eq.kn1) go to 210
       xken=coor(1,ken)
       yken=coor(2,ken)
       if(a*xken+b*yken+c.le.0.0) goto 210
c
c      see if this connection is possible
c
       call possib(kn1,kn,ken,xn1,yn1,xn,yn,xken,yken,
     *               nonf,nonr,npfrt,nqfrt,nregi,coor,iyon)
       if(iyon.eq.0) goto 210
       xkenf=xba(xken,yken,anx,any,alph)
       ykenf=yba(xken,yken,anx,any,alph)
       xdiff1=xkenf-xn1f
       ydiff1=ykenf-yn1f
       xdiff2=xkenf-xnf
       ydiff2=ykenf-ynf
       distf1=sqrt(xdiff1*xdiff1+ydiff1*ydiff1)
       distf2=sqrt(xdiff2*xdiff2+ydiff2*ydiff2)
       cosa=(xdiff1*xdiff2+ydiff1*ydiff2)/(distf1*distf2)
       if(cosa.gt.1.0) cosa=1.0
       if(cosa.lt.-1.0) cosa=-1.0
       angl=acos(cosa)
       if(angl.lt.ang1) go to 210
       howf(inum+1)=angl
       near(inum+1)=ken
       if(inum.eq.0) goto 311
       do 310 i=1,inum
       k=i
       angi=howf(i)
       if(angi.gt.angl) goto 310
       do 410 j=1,inum-k
       l=inum-j
       near(l+1)=near(l)
       howf(l+1)=howf(l)
410    continue
       near(k)=ken
       howf(k)=angl
       goto 311
310    continue
311    continue
       if(inum.lt.30) inum=inum+1
       if(inum.eq.30) ang1=howf(30)
c      if(ibl.eq.1) goto 926
  210  continue
c
c      select ---> start by the closest
c
 926   if(inum.eq.0) goto 703
       wfar=1.e+6
       do 701 i=1,inum
       kp=near(i)
       angi=howf(i)
       xp=coor(1,kp)
       yp=coor(2,kp)
c
c      check the size
c
       xpf=xba(xp,yp,anx,any,alph)
       ypf=yba(xp,yp,anx,any,alph)
       d1=sqrt((xpf-xn1f)**2+(ypf-yn1f)**2)
       d2=sqrt((xpf-xnf)**2+(ypf-ynf)**2)
       di=max(d1,d2)
       if(di.lt.wfar) kp1=kp
       if(di.lt.wfar) wfar=di
701    continue
       kp=kp1
       xp=coor(1,kp)
       yp=coor(2,kp)
       goto 603
703    continue
       if(ibl.eq.1) then
         kp= 0
         xp=xtemp
         yp=ytemp
c
c      see if this connection is possible
c
       call possib(kn1,kn,kp,xn1,yn1,xn,yn,xp,yp,nonf,
     *               nonr,npfrt,nqfrt,nregi,coor,iyon)

       if(iyon.ne.0) goto 603
       end if
       print*,' cannot find the connectivity'
c
c      try another side
c
       lmp(kn) = 0
       lmp(kn1) = 0
       disw1=1.01*alen1f
       call order(nl,npfrt,nqfrt,coor,npoig,neleg,coorg,
     *            ieleg,delta,strec,disw,disw1,qfr,ibl,lmb,ibs,
     *            lmp)
       goto 161
c
c      nothing wrong with it
c
  603  continue
       indic=1
       knear=kp
       if(knear.ne.0) indic=0

      if(indic.eq.0) go to 620
      node=node+1
      if(node.gt.mxpoi) stop 'error 90'
      coor(1,node)=xp
      coor(2,node)=yp
      lmb ( node )=1
      if(ibl.eq.1) lmb ( node )=0
c
c     find out interpolated variables
c
      inorm=1
      ilast = lep(kn)
      call findel(npoig,neleg,coorg,ieleg,intmeg,xp,yp,ilast,
     *            ar,i1,i2,i3,ienr,lelch)
      call getval(npoig,i1,i2,i3,ar,delta,dis,alp,anx,any,inorm)
       xr(1) = xp
       xr(2) = yp
      if(lbk(1).ne.0) then
        do 117 ip = 1,lbk(1)
        i1   = lbk(3)+(ip-1)*5
        sr   = spapt(abk(i1),xr)
        dis = min(dis,sr)    
  117   continue
      endif    
c
      if(lbk(2).ne.0) then
        do 127 ip = 1,lbk(2)
        i1   = lbk(4)+(ip-1)*10
        i2   = i1+5
        sr   = spaln(abk(i1),abk(i2),xr)
        dis = min(dis,sr)    
  127   continue
      endif    
c
      strec(1,node)=dis
      strec(2,node)=alp
      strec(3,node)=anx
      strec(4,node)=any
      lmp(node) = 1
      lmp(kn)  = 1
      lmp(kn1) = 1
      call getval(npoig,i1,i2,i3,ar,unkng,dis,alp,anx,any,0)
      unkno(1,node)=dis
      unkno(2,node)=alp
      unkno(3,node)=anx
      unkno(4,node)=any
      lep  ( node )=ienr
      nonr=nonr+1
      nuno=nonr
      nregi(nuno)=node
      knear=node
      ncheck(knear)=-100
620   continue
c
c       form element
c
      nelem=nelem+1
      if(nelem.gt.mxele) stop 'error 100'
      iel(1,nelem)=kn1
      iel(2,nelem)=kn
      iel(3,nelem)=knear
c
c     plot element boundaries
c
      if(iplot.eq.0) then
      xplot=coor(1,knear)
      yplot=coor(2,knear)
      call plot(4,xplot,yplot,st,xmin,xmax,ymin,ymax)
      do 6 ipoin=1,3
      nodeno=iel(ipoin,nelem)
      xplot=coor(1,nodeno)
      yplot=coor(2,nodeno)
      call plot(5,xplot,yplot,st,xmin,xmax,ymin,ymax)
6     continue
      else
      if(mod(nelem,100).eq.0) then
        st='  nelem ='//convrt(nelem)//'  npoin = '//convrt(node)
     &                               //'  nofrt = '//convrt(nonf)
        call plot(2,x,y,st,xmin,xmax,ymin,ymax)
      end if
      end if
c
c      update front and active nodes
c
      if(lmb(kn).eq.100)itedg = itedg + 1
      if(lmb(kn1).eq.100)itedg = itedg + 1
      nl2=nl+1
      nqfrt(nl2)=kn
      nqfrt(nl)=knear
      npfrt(nl2)=knear
      if(ncheck(knear).lt.0) ncheck(knear)=0
      ncheck(knear)=ncheck(knear)+2
      nonf=nl2
        ks=nqfrt(nl)
        ks1=npfrt(nl)
        xs=coor(1,ks)
        ys=coor(2,ks)
        xs1=coor(1,ks1)
        ys1=coor(2,ks1)
        alph=min(strec(2,ks1),strec(2,ks))
c       if(ibl.eq.1) alph = 1.
        anx=0.5*(strec(3,ks1)+strec(3,ks))
        any=0.5*(strec(4,ks1)+strec(4,ks))
        anm=sqrt(anx*anx+any*any)
        anx=anx/anm
        any=any/anm
        qfr(1,nl)=anx
        qfr(2,nl)=any
        qfr(3,nl)=xba(xs,ys,anx,any,alph)
        qfr(4,nl)=yba(xs,ys,anx,any,alph)
        qfr(5,nl)=xba(xs1,ys1,anx,any,alph)
        qfr(6,nl)=yba(xs1,ys1,anx,any,alph)
        x12f=qfr(3,nl)-qfr(5,nl)
        y12f=qfr(4,nl)-qfr(6,nl)
        qfr(7,nl)=sqrt(x12f*x12f+y12f*y12f)
        ks=nqfrt(nl2)
        ks1=npfrt(nl2)
        xs=coor(1,ks)
        ys=coor(2,ks)
        xs1=coor(1,ks1)
        ys1=coor(2,ks1)
        alph=min(strec(2,ks1),strec(2,ks))
c       if(ibl.eq.1) alph = 1.
        anx=0.5*(strec(3,ks1)+strec(3,ks))
        any=0.5*(strec(4,ks1)+strec(4,ks))
        anm=sqrt(anx*anx+any*any)
        anx=anx/anm
        any=any/anm
        qfr(1,nl2)=anx
        qfr(2,nl2)=any
        qfr(3,nl2)=xba(xs,ys,anx,any,alph)
        qfr(4,nl2)=yba(xs,ys,anx,any,alph)
        qfr(5,nl2)=xba(xs1,ys1,anx,any,alph)
        qfr(6,nl2)=yba(xs1,ys1,anx,any,alph)
        x12f=qfr(3,nl2)-qfr(5,nl2)
        y12f=qfr(4,nl2)-qfr(6,nl2)
        qfr(7,nl2)=sqrt(x12f*x12f+y12f*y12f)
c
c      delete sides from active list
c
      ncht=0
      npass=1
      ntop=kn1
      nbot=knear
      marker=0
      knon=nonf-2
360   do 300 kp=1,knon
      ktop=npfrt(kp)
      kbot=nqfrt(kp)
      if(marker.gt.0) go to 330
      if((ktop.eq.nbot).and.(kbot.eq.ntop)) go to 320
      go to 300
320   marker=1
      nonf=nonf-2
      ncheck(ktop)=ncheck(ktop)-2
      ncheck(kbot)=ncheck(kbot)-2
      ncht=1
      go to 300
330   kp1=kp-1
      npfrt(kp1)=ktop
      nqfrt(kp1)=kbot
      do 397 iz=1,7
        qfr(iz,kp1)=qfr(iz,kp)
 397  continue
300   continue
      npass=npass+1
      if(npass.gt.2) go to 340
      if(marker.eq.0) go to 350
      npfrt(nonf)=knear
      nqfrt(nonf)=kn
      do 497 iz=1,7
        qfr(iz,nonf)=qfr(iz,nl2)
 497  continue
      marker=0
      knon=knon-1
      go to 400
350   knon=knon+1
400   ntop=knear
      nbot=kn
      go to 360
340   continue
c
      knon = 0
      if(ibl.eq.1) then
       if(lmb(kn).eq.2) lmb(kn) = -2  
       if(lmb(kn1).eq.2) lmb(kn1) = -2  
       if(lmb(kn).eq.4) lmb(kn) = -3  
       if(lmb(kn1).eq.4) lmb(kn1) = -3  
       if(lmb(kn).eq.3) lmb(kn) = 4  
       if(lmb(kn1).eq.3) lmb(kn1) = 4  
       do 632 i = 1,nonf
        np= npfrt(i)
        nq= nqfrt(i)
        if(lmb(np).gt.0.or.lmb(nq).gt.0) knon = knon+1
 632   continue
       if(knon.eq.0) return
      end if
c
c     remove nodes from active list
c
      if(ncht.eq.0) go to 370
      ired=0
      knon=nonr
      do 380 kp=1,knon
      kkn=nregi(kp)
      kch=ncheck(kkn)
      if(kch.ne.0) go to 390
      ired=ired+1
      go to 380
390   kp1=kp-ired
      nregi(kp1)=nregi(kp)
380   continue
      nonr=knon-ired
370   continue
      if(nonf.gt.0) go to 160
      return
5000  format(4i5)
      end
c
c-------------------------------------------------------------------
c
c      this subroutine finds out whether connection whith point kp
c      is possible iyon=1 or not iyon=0
c
c-------------------------------------------------------------------
c
       subroutine possib(kn1,kn,kp,xn1,yn1,xn,yn,xp,yp,
     *                     nonf,nonr,npfrt,nqfrt,nregi,
     *                     coor,iyon)
       dimension nregi(*)
       dimension npfrt(*),nqfrt(*)
       dimension coor(2,*)
c      deter(p1,q1,p2,q2,p3,q3)=p2*q3-p3*q2-p1*q3+p3*q1+p1*q2-p2*q1
       deter(p1,q1,p2,q2,p3,q3)=(p2-p1)*(q3-q1)-(p3-p1)*(q2-q1)
       iyon=1
c
c      loop over the front nodes
c
       xmin=min(xn1,xn,xp) - 0.00001
       xmax=max(xn1,xn,xp) + 0.00001
       ymin=min(yn1,yn,yp) - 0.00001
       ymax=max(yn1,yn,yp) + 0.00001
       do 1000 it=1,nonr
         kj=nregi(it)
         if(kj.eq.kn1.or.kj.eq.kn.or.kj.eq.kp) goto 1000
         xt=coor(1,kj)
         yt=coor(2,kj)
c
c      check if the point is interior
c
         if(min(xmin,xt).eq.xt.or.min(ymin,yt).eq.yt.or.
     &      max(xmax,xt).eq.xt.or.max(ymax,yt).eq.yt) goto 1000
         iin=1
         area2=deter(xn1,yn1,xn,yn,xp,yp)
         if(area2.lt.1.e-7) goto 900
         a1=deter(xt,yt,xn,yn,xp,yp)/area2
         a2=deter(xn1,yn1,xt,yt,xp,yp)/area2
         a3=deter(xn1,yn1,xn,yn,xt,yt)/area2
         wcomp=min(a1,a2,a3)
         if(wcomp.lt.-0.00001) iin=0
  900    if(iin.eq.0) goto 1000
         iyon=0
         return
 1000  continue
c
c      equation of the mid-base : kp line
c
       xmb=0.5*(xn1+xn)
       ymb=0.5*(yn1+yn)
       as=ymb-yp
       bs=xp-xmb
       cs=(xmb-xp)*ymb+(yp-ymb)*xmb
c
c      loop over the front sides : check for intersection
c
       do 2000 ir=1,nonf
       knt1=npfrt(ir)
       if(knt1.eq.0) goto 2000
       knt=nqfrt(ir)
       if(knt1.eq.kn1.and.knt.eq.kn) goto 2000
       if(knt1.eq.kp.or.knt.eq.kp) goto 2000
       xnt1=coor(1,knt1)
       ynt1=coor(2,knt1)
       xnt=coor(1,knt)
       ynt=coor(2,knt)
       if(max(xmin,xnt1,xnt).eq.xmin.or.max(ymin,ynt1,ynt).eq.ymin.or.
     & min(xmax,xnt1,xnt).eq.xmax.or.min(ymax,ynt1,ynt).eq.ymax)goto2000
       at=ynt1-ynt
       bt=xnt-xnt1
       ct=(xnt1-xnt)*ynt1+(ynt-ynt1)*xnt1
       s1=at*xmb+bt*ymb+ct
       s2=at*xp+bt*yp+ct
       s3=as*xnt1+bs*ynt1+cs
       s4=as*xnt+bs*ynt+cs
       sig1=s1*s2
       sig2=s3*s4
       if(sig1.gt.0.0.or.sig2.gt.0.0) goto 2000
       iyon=0
       return
2000   continue
       return
       end
c
c-------------------------------------------------------------------
c
c     this subroutine order the front in decending order
c
c--------------------------------------------------------------------
c
      subroutine order(nl,npfrt,nqfrt,coor,npoig,
     *                 neleg,coorg,ieleg,delta,strec,disw,
     *                 disw1,qfr,ibl,lmb,ibs,lmp)
      dimension npfrt(*),nqfrt(*),lmb(*),lmp(*)
      dimension delta(4,*),coor(2,*),strec(4,*)
      dimension coorg(2,*),ieleg(3,*),qfr(7,*)
      dimension t(7)
      disw=1.e+6
      ibs = 0
      iwh = 0
      if(ibl.eq.1) then
        do 1000 il=1,nl
         kn = npfrt(il)
         kn1= nqfrt(il)
         if(lmb(kn).le.0.and.lmb(kn1).le.0) goto 1000
         if(lmb(kn).ge.1.and.lmb(kn1).ge.1) then
           dis = qfr(7,il)
           if(dis.lt.disw1) goto 1000
           if(dis.gt.disw) goto 1000
           if(lmp(npfrt(il)).eq.0.and.lmp(nqfrt(il)).eq.0)goto 1000
           iwh=il
           ibs=1
           disw = dis
         end if
 1000   continue
        if(iwh.ne.0) goto 200
      end if
      do 2 il=1,nl
       kn = npfrt(il)
       kn1= nqfrt(il)
       if(ibl.eq.1.and.lmb(kn).le.0.and.lmb(kn1).le.0) goto 2
       if(lmb(kn).eq.100.and.lmb(kn1).le.0)goto 2
       dis = qfr(7,il)
       if(dis.lt.disw1) goto 2
       if(dis.gt.disw) goto 2
       if(lmp(npfrt(il)).eq.0.and.lmp(nqfrt(il)).eq.0)goto 2
       iwh=il
       disw = dis
 2    continue
      if(iwh.eq.0) then
      do 2000 il=1,nl
      dis=qfr(7,il)
      if(dis.lt.disw1) goto 2000
      if(dis.gt.disw) goto 2000
      iwh=il
      disw=dis
2000  continue
      end if
c
c     swap values
c
 200  ip=npfrt(iwh)
      iq=nqfrt(iwh)
      do 297 iz=1,7
        t(iz)=qfr(iz,iwh)
 297  continue
      npfrt(iwh)=npfrt(nl)
      nqfrt(iwh)=nqfrt(nl)
      do 298 iz=1,7
        qfr(iz,iwh)=qfr(iz,nl)
 298  continue
      npfrt(nl)=ip
      nqfrt(nl)=iq
      do 299 iz=1,7
        qfr(iz,nl)=t(iz)
 299  continue
      return
      end
c
c ------------------------------------------------------------------
c
c   moving points
c
c -----------------------------------------------------------------
c
      subroutine mshmov( intma, coord, lwher, lhowm, icone, lmb,
     &                   deltm,  nlay,  ilay,    n1,    n2, hmin,
     &                   unkno,    io,     h, index, m1, m2 ,
     &                   tn1, tn2 )
c
      dimension intma(3,*)  , coord(2,*) , lmb(*) , lwher(*)
      dimension lhowm( * )  , icone( * ) , unkno(4,*)
      dimension net(3) , nxi(3) ,  x(3) , y(3) , geo(6)
      save ho,ro,uo,po,so,to
      data  nxi/-1.0 , 1.0 , 0.0 /
      data  net/-1.0 , 0.0 , 1.0 /
c
      print*,'minimum spacing on this layer =', deltm
c
      h = hmin*exp((ilay-1)*log(deltm/0.7/hmin)/max(1,(nlay-1)))
      read(1,*) h
c
      if(io.eq.0) goto 9
      dr = 0.0
      du = 0.0
      dp = 0.0
      ds = 0.0
      dt = 0.0
      do 8 ke = n1,n2
       i1 = intma(1,ke)
       i2 = intma(2,ke)
       i3 = intma(3,ke)
       if((lmb(i1).eq.0.and.lmb(i2).gt.0).or.
     &    (lmb(i2).eq.0.and.lmb(i1).gt.0)) then
         d = ((coord(1,i1)-coord(1,i2))**2+
     &        (coord(2,i2)-coord(2,i1))**2)**0.5
         r1 = unkno(1,i1)
         r2 = unkno(1,i2)
         dr = max(dr,abs((r1-r2)/d),0.000001)
         u1 = ((unkno(2,i1)/r1)**2+(unkno(3,i1)/r1)**2)**0.5
         u2 = ((unkno(2,i2)/r2)**2+(unkno(3,i2)/r2)**2)**0.5
         du = max(du,abs((u1-u2)/d),0.000001)
         p1 = 0.4*(unkno(4,i1)-0.5*r1*u1**2)
         p2 = 0.4*(unkno(4,i2)-0.5*r2*u2**2)
         dp = max(dp,abs((p1-p2)/d),0.000001)
         s1 = u1*sqrt(r1/(1.4*p1))
         s2 = u2*sqrt(r2/(1.4*p2))
         ds = max(ds,abs((s1-s2)/d),0.000001)
         t1 = 1.4*(unkno(4,i1)/r1-0.5*u1**2)
         t2 = 1.4*(unkno(4,i2)/r2-0.5*u2**2)
         dt = max(dt,abs((t1-t2)/d),0.000001)
       end if
       if((lmb(i1).eq.0.and.lmb(i3).gt.0).or.
     &    (lmb(i3).eq.0.and.lmb(i1).gt.0)) then
         d = ((coord(1,i1)-coord(1,i3))**2+
     &        (coord(2,i3)-coord(2,i1))**2)**0.5
         r1 = unkno(1,i1)
         r2 = unkno(1,i3)
         dr = max(dr,abs((r1-r2)/d),0.000001)
         u1 = ((unkno(2,i1)/r1)**2+(unkno(3,i1)/r1)**2)**0.5
         u2 = ((unkno(2,i3)/r2)**2+(unkno(3,i3)/r2)**2)**0.5
         du = max(du,abs((u1-u2)/d),0.000001)
         p1 = 0.4*(unkno(4,i1)-0.5*r1*u1**2)
         p2 = 0.4*(unkno(4,i3)-0.5*r2*u2**2)
         dp = max(dp,abs((p1-p2)/d),0.000001)
         s1 = u1*sqrt(r1/(1.4*p1))
         s2 = u2*sqrt(r2/(1.4*p2))
         ds = max(ds,abs((s1-s2)/d),0.000001)
         t1 = 1.4*(unkno(4,i1)/r1-0.5*u1**2)
         t2 = 1.4*(unkno(4,i3)/r2-0.5*u2**2)
         dt = max(dt,abs((t1-t2)/d),0.000001)
       end if
       if((lmb(i3).eq.0.and.lmb(i2).gt.0).or.
     &    (lmb(i2).eq.0.and.lmb(i3).gt.0)) then
         d = ((coord(1,i2)-coord(1,i3))**2+
     &        (coord(2,i3)-coord(2,i2))**2)**0.5
         r1 = unkno(1,i3)
         r2 = unkno(1,i2)
         dr = max(dr,abs((r1-r2)/d),0.000001)
         u1 = ((unkno(2,i3)/r1)**2+(unkno(3,i3)/r1)**2)**0.5
         u2 = ((unkno(2,i2)/r2)**2+(unkno(3,i2)/r2)**2)**0.5
         du = max(du,abs((u1-u2)/d),0.000001)
         p1 = 0.4*(unkno(4,i3)-0.5*r1*u1**2)
         p2 = 0.4*(unkno(4,i2)-0.5*r2*u2**2)
         dp = max(dp,abs((p1-p2)/d),0.000001)
         s1 = u1*sqrt(r1/(1.4*p1))
         s2 = u2*sqrt(r2/(1.4*p2))
         ds = max(ds,abs((s1-s2)/d),0.000001)
         t1 = 1.4*(unkno(4,i3)/r1-0.5*u1**2)
         t2 = 1.4*(unkno(4,i2)/r2-0.5*u2**2)
         dt = max(dt,abs((t1-t2)/d),0.000001)
       end if
   8  continue
      if(ilay.eq.1) then
        ro = dr
        uo = du
        po = dp
        so = ds
        to = dt
        ho = hmin
      end if
c     h  = max(min(ro/dr,uo/du,po/dp,so/ds)*hmin,ho)
c     h = max(min((ro/dr)**2*hmin,2.5*ho),ho)
      if(index.eq.1) then
       tp = ro
       bt = dr
       iz = 1
      end if
      if(index.eq.2) then
       tp = uo
       bt = du
       iz = 2
      end if
      if(index.eq.3) then
       tp = po
       bt = dp
      end if
      if(index.eq.4) then
       tp = so
       bt = ds
       iz = 2
      end if
      if(index.eq.5) then
       tp = to
       bt = dt
       iz = 2
      end if
      ro = max (ro,dr)
      uo = max (uo,du)
      po = max (po,dp)
      so = max (so,ds)
      to = max (to,dt)
c     h = max(min((tp/bt)**iz*hmin,2.5*ho),ho)
      print *,'hmin for layer',ilay,'   is   ',h
      print *,'d for layer',ilay,'   is   ',bt
      ho = h
   9  do 10 ke = n1,n2
       ie = ke
       icount = 0
       j1 = 0
       i1 = intma(1,ie)
       i2 = intma(2,ie)
       i3 = intma(3,ie)
       if(lmb(i1).ge.1.and.lmb(i2).ge.1) j1= 3
       if(lmb(i2).ge.1.and.lmb(i3).ge.1) j1= 1
       if(lmb(i3).ge.1.and.lmb(i1).ge.1) j1= 2
       j2 = j1 + 1
       if(j2.gt.3) j2 = 1
       j3 = j2 + 1
       if(j3.gt.3) j3 = 1
       ind = 0
       if(j1.eq.0) then
         if(lmb(i1).ge.2.and.lmb(i1).ne.100) j1 = 2
         if(lmb(i2).ge.2.and.lmb(i2).ne.100) j1 = 3
         if(lmb(i3).ge.2.and.lmb(i3).ne.100) j1 = 1
         j2 = j1 - 1
         if(j2.lt.1) j2 = 3
         ind = 1
       end if
       if(j1.eq.0)then
c
c ***   trailing edge
c
         if(lmb(i1).eq.100.or.lmb(i2).eq.100.or.lmb(i3).eq.100)then
           if(lmb(i1).eq.100)kp = i1
           if(lmb(i2).eq.100)kp = i2 
           if(lmb(i3).eq.100)kp = i3 
           nst = lwher(kp)
           nhm = lhowm(kp)
           do 667 iel1 = nst+1, nst+nhm
            ielem = icone(iel1)
            if(ielem.lt.n1.or.ielem.gt.n2)goto 667
            kp1 = intma(1,ielem)
            kp2 = intma(2,ielem)
            kp3 = intma(3,ielem)
            kkp = 0
            if(lmb(kp1).eq.1.and.lmb(kp2).eq.100)kkp = kp3
            if(lmb(kp2).eq.1.and.lmb(kp3).eq.100)kkp = kp1
            if(lmb(kp3).eq.1.and.lmb(kp1).eq.100)kkp = kp2
            if(kkp.eq.0)goto 667
            do 668 iel2 = nst+1,nst+nhm
             iell = icone(iel2)
             if(iell.lt.n1.or.iell.gt.n2)goto 668 
             kkp1 = intma(1,iell)
             kkp2 = intma(2,iell)
             kkp3 = intma(3,iell)
             if(kkp1.eq.kkp.and.kkp2.eq.kp)j1 = 3
             if(kkp2.eq.kkp.and.kkp3.eq.kp)j1 = 1
             if(kkp3.eq.kkp.and.kkp1.eq.kp)j1 = 2
             if(j1.eq.0)goto 668
             j3 = j1 + 1
             if(j3.gt.3) j3 = 1
             j2 = j3 + 1
             if(j2.gt.3) j2 = 1
             ind = 0
             ie = iell
             goto 111
 668        continue
 667       continue
         else
           goto 10
         endif
       endif 
 111   ip = intma(j1,ie)
       if(lmb(ip).eq.-10) goto 10
c      if(h.ge.deltm) then
c       lmb(ip) = -10
c       goto 10
c      end if
       if(lhowm(intma(j1,ie)).gt.3) then
         nst = lwher(ip)
         nhm = lhowm(ip)
         als = 10000000.
         ie = 0
         do 610 iee = nst+1,nst+nhm
          iel = icone(iee)
          i1 = intma(1,iel)
          i2 = intma(2,iel)
          i3 = intma(3,iel)
          if(i1.eq.ip) then
            ic = i1
            il = i3
            ir = i2
          end if
          if(i2.eq.ip) then
            ic = i2
            il = i1
            ir = i3
          end if
          if(i3.eq.ip) then
            ic = i3
            il = i2
            ir = i1
          end if
          if(lmb(ir).le.0.or.lmb(il).le.0)goto 610
          d = ((coord(1,ic)-coord(1,ir))**2+
     &         (coord(2,ic)-coord(2,ir))**2)**0.5
          dd = ((coord(1,il)-coord(1,ir))**2+
     &         (coord(2,il)-coord(2,ir))**2)**0.5
          costh = (coord(1,ic)-coord(1,ir))*(coord(1,il)-coord(1,ir))
     &          + (coord(2,ic)-coord(2,ir))*(coord(2,il)-coord(2,ir))
          costh = costh/(d*dd)
          theta = acos(costh)
          hnew = h/sin(theta)
          if((d-hnew).gt.0) then
            u1 = hnew/(d-hnew)
            u2 = 1.
          else
            u1 = 1.0
            u2 = 0.0
          endif
          xb = u1+u2
          yb = u1+u2
          xt = u1*coord(1,ic)+u2*coord(1,ir)
          yt = u1*coord(2,ic)+u2*coord(2,ir)
          xc = xt/xb
          yc = yt/yb
          xo = coord(1,ic)
          yo = coord(2,ic)
          coord(1,ic) = xc
          coord(2,ic) = yc
          jcount = 0
 617      alm = 0.
          do 711 jee = nst+1,nst+nhm
           jel = icone(jee)
           j1 = intma(1,jel)
           j2 = intma(2,jel)
           j3 = intma(3,jel)
           x(1)=coord(1,j1)
           y(1)=coord(2,j1)
           x(2)=coord(1,j2)
           y(2)=coord(2,j2)
           x(3)=coord(1,j3)
           y(3)=coord(2,j3)
           x21=x(2)-x(1)
           x31=x(3)-x(1)
           y21=y(2)-y(1)
           y31=y(3)-y(1)
           rj =x21*y31-x31*y21
           if(rj.lt.0.0) then
             if(jcount.gt.5) then
               coord(1,ic) = xo
               coord(2,ic) = yo
               print *,' jcount=5 for point',ip
               goto  610
             else
               jcount = jcount + 1
               coord(1,ic) = (xo-xc)/6.*jcount+xc
               coord(2,ic) = (yo-yc)/6.*jcount+yc
               goto 617
             end if
           end if
           if(lmb(j1).le.0.and.lmb(j2).le.0) goto 711
           if(lmb(j1).le.0.and.lmb(j3).le.0) goto 711
           if(lmb(j3).le.0.and.lmb(j2).le.0) goto 711
           if(j1.eq.ip) then
             jc = j1
             jl = j3
             jr = j2
           end if
           if(j2.eq.ip) then
             jc = j2
             jl = j1
             jr = j3
           end if
           if(j3.eq.ip) then
             jc = j3
             jl = j2
             jr = j1
           end if
           alr = ((coord(1,jc)-coord(1,jr))**2+
     &            (coord(2,jc)-coord(2,jr))**2)**0.5
           all = ((coord(1,jc)-coord(1,jl))**2+
     &            (coord(2,jc)-coord(2,jl))**2)**0.5
           alm = max(alm,alr,all)
 711      continue
          if(alm.lt.als) then
            als = alm
            ie  = iel
          end if
          coord(1,ic) = xo
          coord(2,ic) = yo
 610     continue   
         if(ie.eq.0) then
           print *,' can not move point ',ip
           print *,(icone(j), j=nst+1,nst+nhm)
           lmb(ip) = -10
           goto 10
         end if
         i1 = intma(1,ie)
         i2 = intma(2,ie)
         i3 = intma(3,ie)
         if(lmb(i1).ge.1.and.lmb(i2).ge.1) j1= 3
         if(lmb(i2).ge.1.and.lmb(i3).ge.1) j1= 1
         if(lmb(i3).ge.1.and.lmb(i1).ge.1) j1= 2
         j2 = j1 + 1
         if(j2.gt.3) j2 = 1
         j3 = j2 + 1
         if(j3.gt.3) j3 = 1
       end if
       if(ind.eq.0) then
         dd = ((coord(1,intma(j3,ie))-coord(1,intma(j2,ie)))**2+
     &        (coord(2,intma(j3,ie))-coord(2,intma(j2,ie)))**2)**0.5
         costh = (coord(1,intma(j1,ie)) - coord(1,intma(j2,ie)))
     &          *(coord(1,intma(j3,ie)) - coord(1,intma(j2,ie)))
     &         + (coord(2,intma(j1,ie)) - coord(2,intma(j2,ie)))
     &          *(coord(2,intma(j3,ie)) - coord(2,intma(j2,ie)))
       else
        costh = 0.
        dd    = 1.
       end if
       d = ((coord(1,ip)-coord(1,intma(j2,ie)))**2+
     &      (coord(2,ip)-coord(2,intma(j2,ie)))**2)**0.5
       costh = costh /(d*dd)
       theta = acos(costh)
       hnew = h/sin(theta)
       if((d-hnew).gt.0) then
         u1 = hnew/(d-hnew)
         u2 = 1.
       else
         u1 = 1.0
         u2 = 0.0
       endif
       xb = u1+u2
       yb = u1+u2
       xt = u1*coord(1,ip)+u2*coord(1,intma(j2,ie))
       yt = u1*coord(2,ip)+u2*coord(2,intma(j2,ie))
       xo = coord(1,ip)
       yo = coord(2,ip)
       xn = xt/xb
       yn = yt/yb
       coord(1,ip) = xn
       coord(2,ip) = yn
       lmb  ( ip ) =-10
c
  607  nst = lwher(ip)
       nhm = lhowm(ip)
       do 620 iel = nst+1,nst+nhm
        ielem = icone(iel)
        in1=intma(1,ielem)
        in2=intma(2,ielem)
        in3=intma(3,ielem)
        x(1)=coord(1,in1)
        y(1)=coord(2,in1)
        x(2)=coord(1,in2)
        y(2)=coord(2,in2)
        x(3)=coord(1,in3)
        y(3)=coord(2,in3)
c
        x21=x(2)-x(1)
        x31=x(3)-x(1)
        y21=y(2)-y(1)
        y31=y(3)-y(1)
c
        rj =x21*y31-x31*y21
        rj1=1./rj
c
        xix= y31*rj1
        xiy=-x31*rj1
        etx=-y21*rj1
        ety= x21*rj1
c
c *** form n,x & n,y
c
        rnxi=nxi(1)
        rnet=net(1)
        geo(1  )=xix*rnxi+etx*rnet
        geo(1+3)=xiy*rnxi+ety*rnet
        rnxi=nxi(2)
        rnet=net(2)
        geo(2  )=xix*rnxi+etx*rnet
        geo(2+3)=xiy*rnxi+ety*rnet
        rnxi=nxi(3)
        rnet=net(3)
        geo(3  )=xix*rnxi+etx*rnet
        geo(3+3)=xiy*rnxi+ety*rnet
c
c *** jacobian
c
        anx=geo(1)
        any=geo(1+3)
        dist1=1.0/sqrt(anx*anx+any*any)
        anx=geo(2)
        any=geo(2+3)
        dist2=1.0/sqrt(anx*anx+any*any)
        anx=geo(3)
        any=geo(3+3)
        dist3=1.0/sqrt(anx*anx+any*any)
        rel=min(dist1,dist2,dist3)
        if(rj.le.0) then
         if(icount.gt.5) then
           coord(1,ip) = xo
           coord(2,ip) = yo
           goto  10
         else
           icount = icount + 1
           print *,' creating negative volume by moving point',ip
           coord(1,ip) = (xo-xn)/6.*icount+xn
           coord(2,ip) = (yo-yn)/6.*icount+yn
           goto 607
         end if
        end if
 620   continue
  10  continue
      do 112 ip = m1,m2
       if(lmb(ip).ne.0)goto 112
       nst = lwher(ip)
       nhm = lhowm(ip)
       do 113 iel = nst+1,nst+nhm
        ie = icone(iel)
        if(ie.lt.n1.or.ie.gt.n2)goto 113
        i1 = intma(1,ie)
        i2 = intma(2,ie)
        i3 = intma(3,ie)
        isum = lmb(i1)+lmb(i2)+lmb(i3)
        if(isum.ge.90)then
          if(lmb(i1).eq.100)it = i1
          if(lmb(i2).eq.100)it = i2
          if(lmb(i3).eq.100)it = i3
          hlen = (coord(1,ip)-coord(1,it))*tn1 +
     &           (coord(2,ip)-coord(2,it))*tn2
c         coord(1,ip) = coord(1,it) + 0.5*hlen*tn1
c         coord(2,ip) = coord(2,it) + 0.5*hlen*tn2
          coord(1,ip) = coord(1,it) + 1.0*hlen*tn1
          coord(2,ip) = coord(2,it) + 1.0*hlen*tn2
          lmb(ip)= -100
        endif
 113   continue
 112  continue 
c
      return
      end
c 
c--------------------------------------------------------------------
c
c     plot generated mesh
c
c--------------------------------------------------------------------
c
      subroutine gmplot( node,coor,nelem,iel,iq,nel)
      dimension coor(2,*)
      dimension iel(3,*)
      character*80 st
c
c       plot element boundaries
c
      if(iq.eq.1) then 
      do 20 ie=1,nel,2
        nodeno=iel(3,ie+1)
        xplot=coor(1,nodeno)
        yplot=coor(2,nodeno)
        call plot(4,xplot,yplot,st,xmin,xmax,ymin,ymax)
        do 30 ipoin=1,3
          nodeno=iel(ipoin,ie)
          xplot=coor(1,nodeno)
          yplot=coor(2,nodeno)
          call plot(5,xplot,yplot,st,xmin,xmax,ymin,ymax)
   30   continue
          nodeno=iel(3,ie+1)
          xplot=coor(1,nodeno)
          yplot=coor(2,nodeno)
          call plot(5,xplot,yplot,st,xmin,xmax,ymin,ymax)
   20 continue
      end if
      nst=nel
      if(iq.eq.0) nst=0
      do 25 ie=nst+1,nelem
        nodeno=iel(3,ie)
        xplot=coor(1,nodeno)
        yplot=coor(2,nodeno)
        call plot(4,xplot,yplot,st,xmin,xmax,ymin,ymax)
        do 35 ipoin=1,3
          nodeno=iel(ipoin,ie)
          xplot=coor(1,nodeno)
          yplot=coor(2,nodeno)
          call plot(5,xplot,yplot,st,xmin,xmax,ymin,ymax)
   35   continue
   25 continue
      return
      end
c
c--------------------------------------------------------------------
c
c    this  subroutine   finds   out  in  which element does the point
c    (x,y) lie, and calculates the area coordinates for interpolation
c
c--------------------------------------------------------------------
c
      subroutine findel(npoig,neleg,coorg,ieleg,intmeg,x,y,ilast,
     1                  ar,i1,i2,i3,ie,lelch)
c
c     this  subroutine   finds   out  in  which element does the point
c     (x,y) lie, and calculates the area coordinates for interpolation
c
      dimension coorg(2,*),ieleg(3,*)
      dimension intmeg(3,*)
      dimension ar(3),ilo(3),or(3)
      dimension lelch(*)
c
c     determinant function
c
c     deter(p1,q1,p2,q2,p3,q3)=p2*q3-p3*q2-p1*q3+p3*q1+p1*q2-p2*q1
      deter(p1,q1,p2,q2,p3,q3)=(p2-p1)*(q3-q1)-(p3-p1)*(q2-q1)
c *** initialise the list of elements which will be tried
c
c
      do 5 i=1,neleg
        lelch(i)=0
 5    continue
c
    9 continue
 1000 ie=ilast
      lelch(ie)=lelch(ie)+1
      if(lelch(ie).gt.3)go to 2
      i1=ieleg(1,ie)
      i2=ieleg(2,ie)
      i3=ieleg(3,ie)
      x1=coorg(1,i1)
      x2=coorg(1,i2)
      x3=coorg(1,i3)
      y1=coorg(2,i1)
      y2=coorg(2,i2)
      y3=coorg(2,i3)
c     area coordinates
c
      area2=deter(x1,y1,x2,y2,x3,y3)
      ar(1)=deter(x ,y ,x2,y2,x3,y3)/area2
      ar(2)=deter(x1,y1,x ,y ,x3,y3)/area2
      ar(3)=deter(x1,y1,x2,y2,x ,y )/area2
c
c     order them
c
      do 101 i=1,3
      or(i)=ar(i)
      ilo(i)=i
  101 continue
c
      do 102 i=1,3
      j1=i+1
      do 102 j=j1,3
      if(or(i)-or(j)) 102,102,103
  103 temp=or(i)
      itemp=ilo(i)
      or(i)=or(j)
      ilo(i)=ilo(j)
      or(j)=temp
      ilo(j)=itemp
  102 continue
c
      if(or(1).ge.-1.e-3) goto 2000
c
c     get the next element
c
      inext=intmeg(ilo(1),ie)
      if(inext.ne.0) then
            ilast=inext
            go to 1000
       else
            inext=intmeg(ilo(2),ie)
            if(or(2).lt.0.0.and.inext.ne.0)then
              ilast=inext
              goto 1000
            else
             go to 2
            endif
         endif
c
c    look for the closest nodal point
c
 2          adis=1.e+6
            do 1010 ip=1,npoig
            xp=coorg(1,ip)
            yp=coorg(2,ip)
            adi1=(x-xp)**2+(y-yp)**2
            if(adi1.gt.adis) goto 1010
            adis=adi1
            kpo=ip
 1010       continue
c
c      now get an element containing this point
c
2001           iex=0
  20           iex=iex+1
 1119          in=1
1020           ipo=ieleg(in,iex)
                if(ipo.eq.kpo)then
                    ilast=iex
                    i1 = kpo
                    i2 = 1
                    i3 = 1
                    ar(1) = 1.
                    ar(2) = 0.
                    ar(3) = 0.
                    return
                endif
            if(in.lt.3) then
                 in=in+1
                 goto 1020
            end if
            if (iex.lt.neleg) then
                 iex=iex+1
                 goto 1119
            end if
c
c *** find the next closest point
       go to 2
c
 2000 continue
c
      return
      end
c
c--------------------------------------------------------------------
c
c     interpolate
c
c--------------------------------------------------------------------
c
      subroutine getval(npoig,i1,i2,i3,ar,delta,dis,alp,anx,any,
     *                    inorm)
      dimension delta(4,*),ar(3)
      if(ar(1).lt.0) ar(1) = 0.
      if(ar(2).lt.0) ar(2) = 0.
      if(ar(3).lt.0) ar(3) = 0.
      a = ar(1)+ar(2)+ar(3)
      ar(1) = ar(1)/a
      ar(2) = ar(2)/a
      ar(3) = ar(3)/a
      dis=ar(1)*delta(1,i1)+ar(2)*delta(1,i2)+ar(3)*delta(1,i3)
      alp=ar(1)*delta(2,i1)+ar(2)*delta(2,i2)+ar(3)*delta(2,i3)
      anx=ar(1)*delta(3,i1)+ar(2)*delta(3,i2)+ar(3)*delta(3,i3)
      any=ar(1)*delta(4,i1)+ar(2)*delta(4,i2)+ar(3)*delta(4,i3)
c
c     normalize if required
c
      if(inorm.eq.0) return
      anm=sqrt(anx*anx+any*any)
      anx=anx/anm
      any=any/anm
      return
      end
c
c-------------------------------------------------------------------
c
c     smooth out the grid in nsmoo steps
c
c-------------------------------------------------------------------
c

       subroutine smooth(ndimn,nnode,nelem,npoin,nsmoo,
     &                   lcoor,lcore,intmat,coord,coor0,
     &                   nel,npl,lwher,lhowm,icone)
c
       dimension    coord(2,*),coor0(2,*)
       dimension    rdivn(30),x(6,2)
       dimension lcoor(*),lcore(*)
       dimension intmat(3,*),node(6)
       dimension lhowm(*),lwher(*),icone(*)
c
c      deter(p1,q1,p2,q2,p3,q3)=p2*q3-p3*q2-p1*q3+p3*q1+p1*q2-p2*q1
       deter(p1,q1,p2,q2,p3,q3)=(p2-p1)*(q3-q1)-(p3-p1)*(q2-q1)
       ndivn=30
c
       do 500 idivn=1,ndivn
       rdivn(idivn)=1./real(idivn)
  500  continue
c
c     -----smooth out the grid in nsmoo steps
c
       if(nsmoo.eq.0) goto 10001
c
       do 10000 ismoo=1,nsmoo
c
c     ----set coor0=0
c
       do 1200 ip=1,npoin
       do 1201 id=1,2
       coor0(id,ip)=0.0
1201    continue
1200    continue
c
c     ----loop over the elements
c
       do 2000 ielem=1,nelem
        do 2100 ic=1,nnode
        in=intmat(ic,ielem)
        do 2101 id=1,ndimn
        x(      ic,id)=coord(id,in)
        x(nnode+ic,id)=coord(id,in)
 2101  continue
        node(ic)=in
 2100  continue
c
        do 2200 ic=1,nnode
        in=node(ic)
        if(lcoor(in).ne.0) goto 2199
         do 2300 jc=1,nnode-1
         do 2301 id=1,ndimn
         coor0(id,in)=coor0(id,in)+x(ic+jc,id)
 2301   continue
 2300   continue
 2199  continue
 2200  continue
c
 2000 continue
c
       do 3000 icoor=npl+1,npoin
       if(lcoor(icoor).ne.0) goto 3000
       is=2*lcore(icoor)
       cn=rdivn(is)
        xold=coord(1,icoor)
        yold=coord(2,icoor)
        do 3100 id=1,ndimn
        coord(id,icoor)=cn*coor0(id,icoor)
 3100  continue
        nelsur=lhowm(icoor)
        nstart=lwher(icoor)+1
        do 64 ncone=nstart,nstart+nelsur-1
          ie=icone(ncone)
          i1=intmat(1,ie)
          i2=intmat(2,ie)
          i3=intmat(3,ie)
c
          x1=coord(1,i1)
          x2=coord(1,i2)
          x3=coord(1,i3)
          y1=coord(2,i1)
          y2=coord(2,i2)
          y3=coord(2,i3)
c
      if(deter(x1,y1,x2,y2,x3,y3).gt.0.0) goto 64
      coord(1,icoor)=xold
      coord(2,icoor)=yold
      go to 3000
 64   continue
 3000 continue
c
10000 continue
10001 continue
c
c
      return
      end
c
c--------------------------------------------------------------------
c
c     this subroutine swaps the diagonals to improve the mesh   
c
c-------------------------------------------------------------------
c
      subroutine swapdi(npoin,nelem,nside,iside,intmat,lcore,
     1                  lcoid,lcoor,coord,nel,npl)
      dimension iside(4,*),intmat(3,*),lcoor(*)
      dimension lcore(*),lcoid(*),coord(2,*)
      character*80 st
      character*5 convrt, stnum
c
c     ij=1   :  obtain a more even distribution per node
c     ij=0   :  maximize the minimum angle per element
c   
c     determinant function
c
c     deter(p1,q1,p2,q2,p3,q3)=p2*q3-p3*q2-p1*q3+p3*q1+p1*q2-p2*q1
       deter(p1,q1,p2,q2,p3,q3)=(p2-p1)*(q3-q1)-(p3-p1)*(q2-q1)
c
c     angle function
c
      co(x1,y1,x2,y2,x3,y3)=((x3-x2)*(x1-x2)+(y3-y2)*(y1-y2))/
     *   sqrt(((x3-x2)**2+(y3-y2)**2)*((x1-x2)**2+(y1-y2)**2))
c
      ij=0
 1000 ichan=0
c
c     loop over the sides
c
      print *,nel,npl
      do 2000 is=1,nside
      i1=iside(1,is)
      i2=iside(2,is)
      if(i1.le.npl.or.i2.le.npl) goto 2000
      ie1=iside(3,is)
      ie2=iside(4,is)
      if(ie1.le.nel.or.ie2.le.nel) goto 2000
c
c     check for boundary sides
c
      if(ie2.eq.0) goto 2000
c
c     determine i3 & i4
c
      do 5000 in=1,3
      in1=in+1
      if(in1.gt.3) in1=in1-3
      in2=in+2
      if(in2.gt.3) in2=in2-3
      ip11=intmat(in,ie1)
      ip12=intmat(in1,ie1)
      ip13=intmat(in2,ie1)
      if(ip11.eq.i1.and.ip12.eq.i2) i3=ip13
      ip21=intmat(in,ie2)
      ip22=intmat(in1,ie2)
      ip23=intmat(in2,ie2)
      if(ip21.eq.i2.and.ip22.eq.i1) i4=ip23
 5000 continue
c
c     find 'deficit' or 'superhavit' of connectivities
c
      if(ij.eq.1)then
        ih1=abs(lcore(i1)-lcoid(i1))
        ih2=abs(lcore(i2)-lcoid(i2))
        ih3=abs(lcore(i3)-lcoid(i3))
        ih4=abs(lcore(i4)-lcoid(i4))
        ihf1=abs(lcore(i1)-1-lcoid(i1))
        ihf2=abs(lcore(i2)-1-lcoid(i2))
        ihf3=abs(lcore(i3)+1-lcoid(i3))
        ihf4=abs(lcore(i4)+1-lcoid(i4))
c
c     check if it is worth to swap
c
        iswap=0
        iactu=ih1+ih2+ih3+ih4
        ifutu=ihf1+ihf2+ihf3+ihf4
        if(iactu.gt.ifutu) iswap=1
        iam=max(ih1,ih2,ih3,ih4)
        ifm=max(ihf1,ihf2,ihf3,ihf4)
        if(iactu.eq.ifutu.and.iam.gt.(ifm+1)) iswap=1
        if(iswap.eq.0) goto 2000
      else
c
c     find angles
c
        cosa1=co(coord(1,i1),coord(2,i1),coord(1,i2),
     *       coord(2,i2),coord(1,i3),coord(2,i3))
        cosa2=co(coord(1,i1),coord(2,i1),coord(1,i2),
     *       coord(2,i2),coord(1,i4),coord(2,i4))
        cosa3=co(coord(1,i2),coord(2,i2),coord(1,i1),
     *       coord(2,i1),coord(1,i4),coord(2,i4))
        cosa4=co(coord(1,i2),coord(2,i2),coord(1,i1),
     *       coord(2,i1),coord(1,i3),coord(2,i3))
        alph1=cosa1
        if(cosa2.gt.alph1)alph1=cosa2
        if(cosa3.gt.alph1)alph1=cosa3
        if(cosa4.gt.alph1)alph1=cosa4
 
        cosa1=co(coord(1,i2),coord(2,i2),coord(1,i4),
     *       coord(2,i4),coord(1,i3),coord(2,i3))
        cosa2=co(coord(1,i1),coord(2,i1),coord(1,i4),
     *       coord(2,i4),coord(1,i3),coord(2,i3))
        cosa3=co(coord(1,i1),coord(2,i1),coord(1,i3),
     *       coord(2,i3),coord(1,i4),coord(2,i4))
        cosa4=co(coord(1,i2),coord(2,i2),coord(1,i3),
     *       coord(2,i3),coord(1,i4),coord(2,i4))
        alph2=cosa1
        if(cosa2.gt.alph2)alph2=cosa2
        if(cosa3.gt.alph2)alph2=cosa3
        if(cosa4.gt.alph2)alph2=cosa4
        iswap=0
        if(alph1.gt.alph2)iswap=1
        if(iswap.eq.0)go to 2000
      endif
c
c     check area
c
      x1=coord(1,i1)
      x2=coord(1,i4)
      x3=coord(1,i3)
      y1=coord(2,i1)
      y2=coord(2,i4)
      y3=coord(2,i3)
      if(deter(x1,y1,x2,y2,x3,y3).le.0.0000001) goto 2000
      x1=coord(1,i3)
      x2=coord(1,i4)
      x3=coord(1,i2)
      y1=coord(2,i3)
      y2=coord(2,i4)
      y3=coord(2,i2)
      if(deter(x1,y1,x2,y2,x3,y3).le.0.0000001) goto 2000
c
c     swap
c
      ichan=ichan+1
      intmat(1,ie1)=i1
      intmat(2,ie1)=i4
      intmat(3,ie1)=i3
      intmat(1,ie2)=i3
      intmat(2,ie2)=i4
      intmat(3,ie2)=i2
      lcore(i1)=lcore(i1)-1
      lcore(i2)=lcore(i2)-1
      lcore(i3)=lcore(i3)+1
      lcore(i4)=lcore(i4)+1
c
c     detect the sides i1-i4 and i3-i2
c
      ist1=0
      ist2=0
      do 4000 ist=1,nside
      i1t=iside(1,ist)
      i2t=iside(2,ist)
      if((i1t.eq.i1.and.i2t.eq.i4).or.(i1t.eq.i4.and.i2t.eq.i1))
     *ist1=ist
      if((i1t.eq.i3.and.i2t.eq.i2).or.(i1t.eq.i2.and.i2t.eq.i3))
     *ist2=ist
      if(ist1.ne.0.and.ist2.ne.0) goto 4001
 4000 continue
 4001 continue
      if(iside(3,ist1).eq.ie2) iside(3,ist1)=ie1
      if(iside(4,ist1).eq.ie2) iside(4,ist1)=ie1
      if(iside(3,ist2).eq.ie1) iside(3,ist2)=ie2
      if(iside(4,ist2).eq.ie1) iside(4,ist2)=ie2
c
c     update iside
c
      iside(1,is)=i4
      iside(2,is)=i3

 2000 continue
      if (ichan.eq.1) then
        st = '           side has been swapped.'
      else
        st = '          sides have been swapped.'
      end if
      stnum = convrt (ichan)
      do 10, i3=1, 5
        st (i3:i3) = stnum (i3:i3)
10    continue
      call plot(2,x,y,st,xmin,xmax,ymin,ymax)
      if(ichan.ge.1) goto 1000
      return
      end
c
c--------------------------------------------------------------------
c
c       this subroutine removes the points where only three
c       elements coincide
c
c--------------------------------------------------------------------
c
        subroutine eat3(npoin,nelem,nside,intmat,iside,lcoor,lposi,
     *                  lwher,lhowm,lcore,lcoid,icone,coord,nel,npl,
     *                  lep  ,strec, lboud,unkno)
        character*5 convrt, stnum
        character*80 st
        dimension  intmat(3,*),iside(4,*),lcoor(*)
        dimension  lposi(*),lwher(*),lhowm(*),lboud(*)
        dimension  lcore(*),lcoid(*),icone(*)
        dimension coord(2,*),lep(*),strec(4,*),unkno(4,*)
        ntres=0
        kpoin=0
c
c       loop over number of nodes
c
        do 1000 ip=1,npoin
          kpoin=kpoin+1
          lposi(ip)=kpoin
c
c         if boundary point leave it
c
          if(lcoor(ip).ne.0) goto 1000
c
c       check number of elements
c
          if(lcore(ip).ne.3) goto 1000
c
c       a point with only three elements
c
          ntres=ntres+1
          kpoin=kpoin-1
          lposi(ip)=0
c
c       get the elements from icone
c
          iloca=lwher(ip)
          ie1=icone(iloca+1)
          ie2=icone(iloca+2)
          ie3=icone(iloca+3)
c
c       get new conectivity point for ie1 from ie2
c
          ip1=intmat(1,ie1)
          ip2=intmat(2,ie1)
          ip3=intmat(3,ie1)
          do 1002 in=1,3
            ipt=intmat(in,ie2)
            if(ipt.ne.ip1.and.ipt.ne.ip2.and.ipt.ne.ip3) ino=ipt
1002      continue
c
c       replace conectivity
c
          do 1003 in=1,3
            if(intmat(in,ie1).eq.ip) intmat(in,ie1)=ino
            lcore(intmat(in,ie1))=lcore(intmat(in,ie1))-1
1003      continue
          do 1004 in=1,3
            intmat(in,ie2)=0
            intmat(in,ie3)=0
1004      continue
c
c       end loop over points
c
1000    continue
c
c       transfer
c
        do 2000 ip=1,npoin
          il=lposi(ip)
          if(il.eq.0) goto 2000
          lcoor(il)=lcoor(ip)
          lcore(il)=lcore(ip)
          lcoid(il)=lcoid(ip)
          lep(il)  = lep(ip)
          lboud(il)  = lboud(ip)
          strec(1,il) = strec(1,ip)
          strec(2,il) = strec(2,ip)
          strec(3,il) = strec(3,ip)
          strec(4,il) = strec(4,ip)
          unkno(1,il) = unkno(1,ip)
          unkno(2,il) = unkno(2,ip)
          unkno(3,il) = unkno(3,ip)
          unkno(4,il) = unkno(4,ip)
          do 2001 id=1,2
            coord(id,il)=coord(id,ip)
2001      continue
2000    continue
        je=0
        do 3000 ie=1,nelem
          if(intmat(1,ie).eq.0) goto 3000
          je=je+1
          do 3001 in=1,3
            iold=intmat(in,ie)
            inew=lposi(iold)
            intmat(in,je)=inew
3001      continue
3000    continue
c
c       get npoin and nelem
c
        npoin=npoin-ntres
        nelem=nelem-2*ntres
c
c       output nr. of eaten points
c
        st = '              number of 3s removed = '
        stnum = convrt (ntres)
        do 10, i3=1, 5
          st (37+i3:37+i3) = stnum (i3:i3)
10      continue 
        call plot(2,x,y,st,xmin,xmax,ymin,ymax)
c
c       fill in iside
c
        call side(nelem,npoin,nside,intmat,iside,lwher,lhowm,icone,0)
        return
        end
c
c-------------------------------------------------------------------
c
c            this subroutine is an output-for-input routine
c
c--------------------------------------------------------------------
c
      subroutine oustar(npoin,nelem,coord,intmat,lpoin,
     *                  iside,nside,iseg ,ibsid ,unkno,strec,
     *                  ieq,lplay,lelay,nlay,iquad,lafter)
c
      dimension ibsid(4,*),iseg(*),unkno(4,*),strec(4,*)
      dimension coord(2,*),intmat(3,*),lpoin(*),iside(4,*)
      dimension lplay(*)  , lelay(*) , ieq(4,*) , lafter(*)
      dimension list(100)
c
      nb = 0
      do 2500 is=1,nside
       if(iside(4,is).ne.0) goto 2500
       nb =nb + 1
       i1=iside(1,is)
       i2=iside(2,is)
       ie=iside(3,is)
       ib=lpoin(i1)
       if(ib.eq.00) ib=lpoin(i2)
       ibsid(1,nb) = i1
       ibsid(2,nb) = i2
       ibsid(3,nb) = ie
       ibsid(4,nb) = ib
       iseg ( nb ) = ib
2500  continue
      do i = 1 , npoin
       lpoin(i) = 0
      end do

      npimp = 0
      do i = 1 , nb
       ist = ibsid(1,i)
       is  = ibsid(1,i)
       if(lafter(is).ne.0) npimp = npimp + 1
 32    if(lafter(is).ne.0) then
         write(196,*) is,coord(1,is),coord(2,is)
         lpoin(ist) = lpoin(ist) + 1
         is = lafter(is)
         goto 32
       end if
      end do
      write(96,*) npimp
      do i = 1 , nb
       is = ibsid(1,i)
       if(lpoin(is).ne.0) then
         list(1) = lpoin(is)
         list(2) = is
         nlst = 2
 34      if(lafter(is).ne.0) then
           is = lafter(is)
           nlst = nlst + 1
           list(nlst) = is
           goto 34
         end if
         write(96,'(52I7)')(list(j),j=1,nlst)
       end if
      end do
c
      if(iquad.eq.2) then
      neq = 0
      do il = 1 , nlay
       do i = 1 , nside
        i1 = iside(1,i)
        i2 = iside(2,i)
        ie1= iside(3,i)
        ie2= iside(4,i)
        if(ie1.eq.0.or.ie2.eq.0) goto 11
        if(lelay(ie1).le.0.or.lelay(ie2).le.0) goto 11
        if(lelay(ie1).ne.lelay(ie2)) goto 11
        i3 = intmat(1,ie1)+intmat(2,ie1)+intmat(3,ie1)-i1-i2
        i4 = intmat(1,ie2)+intmat(2,ie2)+intmat(3,ie2)-i1-i2           
        ism = 100000
        if(lplay(i1).gt.0.and.ism.gt.lplay(i1)) ism = lplay(i1)
        if(lplay(i2).gt.0.and.ism.gt.lplay(i2)) ism = lplay(i3)
        if(lplay(i3).gt.0.and.ism.gt.lplay(i3)) ism = lplay(i3)
        if(lplay(i4).gt.0.and.ism.gt.lplay(i4)) ism = lplay(i4)
        ing = 0
        if(lplay(i1).lt.0) ing = ing + 1
        if(lplay(i2).lt.0) ing = ing + 1
        if(lplay(i3).lt.0) ing = ing + 1
        if(lplay(i4).lt.0) ing = ing + 1
c       ism= min(lplay(i1),lplay(i2),lplay(i3),lplay(i4))
        if(ism.ne.il.or.ing.eq.2) goto 11
        neq = neq + 1
        lelay(ie1) = -1
        lelay(ie2) = -1
        ieq(1,neq) = i1
        ieq(2,neq) = i4
        ieq(3,neq) = i2
        ieq(4,neq) = i3
  11    continue
       end do
       print *, 'layer=',il,' nquads=',neq
      end do
        
      write(19,*) 1
      write(19,'(a)') 'title'
      write(19,'(a)') 'nelem   npoin   nboun'
      write(19,*) neq*4+(nelem-2*neq),npoin,nb
      write(19,'(a)') '  connectivities '
      ns = 0 
      do 2000 i=1,neq
       do k = 1 , 4
        i1 = ieq(k,i)
        kk = k + 1
        if(kk.gt.4) kk = 1
        i2 = ieq(kk,i)
        i3 = i1
        ns = ns + 1
        write(19,3)ns,i1,i2,i3,0
       end do
 2000 continue
      do 2001 i=1,nelem
       if(lelay(i).ne.-1) then
         ns = ns + 1
         write(19,3)ns,(intmat(j,i),j=1,3),0
       end if
 2001 continue
      else
        neq = 0
        do  i = 1 , nelem
         lelay(i) = 0
        end do
      end if
c
      write(9,*) 1
      write(9,'(a)') 'title'
      write(9,'(a)') 'nelem   npoin   nboun'
      if(iquad.eq.2) then
        write(9,*) neq,nelem-2*neq,npoin,nb
      else
        write(9,*) nelem,npoin,nb
      end if
      write(9,'(a)') '  connectivities '
      do 2002 i=1,neq
       write(9,3)i,(ieq(j,i),j=1,4)
 2002 continue
      do 2003 i=1,nelem
       if(lelay(i).ne.-1) then
         neq = neq + 1
         write(9,3)neq,(intmat(j,i),j=1,3),0
       end if
 2003 continue
c
      write(9,'(a)') '  coordinates '
      if(iquad.eq.2) write(19,'(a)') '  coordinates '
      do 1008 i=1,npoin
       sp = strec(1,i)
       st = strec(2,i)
       write(9,5)i,(coord(j,i),j=1,2),sp,st
       if(iquad.eq.2) write(19,5)i,(coord(j,i),j=1,2),sp,st
1008  continue
      write(9,'(a)') '  unknown '
      if(iquad.eq.2) write(19,'(a)') '  unknown '
      do 1005 i=1,npoin
       ro=unkno(1,i)
       uv=unkno(2,i)
       vv=unkno(3,i)
       en=unkno(4,i)
1202   write(9,5)i,ro,uv,vv,en
       if(iquad.eq.2)write(19,5)i,ro,uv,vv,en
1005  continue
      write(9,'(a)') '  boundary sides '
      if(iquad.eq.2) write(19,'(a)') '  boundary sides '
      do 1009 i = 1 , nb
       write(9,3)(ibsid(j,i),j=1,4),iseg(i)
       if(iquad.eq.2)write(19,3)(ibsid(j,i),j=1,4),iseg(i)
1009  continue
      return
1     format(5i8)
3     format(i10,i10,10i10)
5     format(i10,6e17.8)
      end
c
c-------------------------------------------------------------------
c
c      subroutine to reduce the spacing of the background grid
c
c ------------------------------------------------------------------
c
      subroutine refin(neleg ,npoig ,coorg ,delta )
c
      dimension coorg(2,*) ,delta(4,*)
      character yn
c
 10   write(*,100)
      read(*,'(a)',err=10)yn      
      if(yn.eq.'n'.or.yn.eq.'N')return
c
      xmax=-1.0e+6
      xmin= 1.0e+6
      ymax=-1.0e+6
      ymin= 1.0e+6
c
      do 20 ip=1,npoig
        xl=coorg(1,ip)
        yl=coorg(2,ip)
        xmin = min (xmin, xl)
        xmax = max (xmax, xl)
        ymin = min (ymin, yl)
        ymax = max (ymax, yl)
 20   continue
      write (*,200)xmin,ymin,xmax,ymax
 30   write (*,300)
      write (*,350)
      read  (*,*,err=30)x1,y1
      write (*,400)
      read  (*,*,err=30)x2,y2
 40   write (*,500)
      read  (*,*,err=40)fact1,fact2
      do 60 ip=1,npoig 
        if(coorg(1,ip).ge.x1.and.coorg(1,ip).le.x2.and.
     *     coorg(2,ip).ge.y1.and.coorg(2,ip).le.y2)then
          d2=delta(1,ip)
          d1=delta(2,ip)*delta(1,ip)
          d1=d1/fact1
          d2=d2/fact2
          delta(1,ip)=d2
          delta(2,ip)=d1/d2
        endif
 60   continue
 70   write (*,700)
      read  (*,'(a)',err=70)yn
      if(yn.eq.'n'.or.yn.eq.'N')then
        return
      else
        go to 30
      endif
c
100   format('  do you wish to alter the spacing in the ', 
     *       '  background grid (Y/N) ? ',$)
200   format('  the background grid is contained in the window',/
     *       ' xmin =  ',e15.5,',   ymin =  ',e15.5,/
     *       ' xmax =  ',e15.5,',   ymax =  ',e15.5,/)
300   format('  define the window where you would like ', 
     *       '  to alter the spacing :',$)
350   format('  input xmin,ymin  :  ',$)
400   format('  input xmax,ymax  :  ',$)
500   format('  by what factor would you like to reduce the spacing ?',
     *       /,'  first direction,   second direction  : ',$)
700   format('  would you like to define another window  (Y/N) ?', $)
c
      return
      end
c
c -----------------------------------------------------------------------
c
c     this sub. generates points on a ferguson spline
c     defined by npi points according to the spacing defined in
c     the information points xip. sip stores 1./spacing
c
c --------------------------------------------------------------------------
c
      subroutine split(ndim,nps,npn,xsp,tsp,xnp,npbg,nebg,spmin,
     *                 coorg,ieleg,intmeg,delta ,lelch ,ilast,ibl,
     *                  lbk,abk)
c
c
      parameter (naux = 1000000 ,ndmx = 2)
      dimension xsp(ndim,*),tsp(ndim,*),xnp(ndim,*)
      dimension r1(2*ndmx),r2(2*ndmx),r(2*ndmx)
      dimension aux(naux),al(naux),coef(5,naux)
      dimension xne(naux),xl(naux),sip(naux)
      dimension tip(ndmx,naux),xip(ndmx,naux)
      dimension coorg(2,*),ieleg(3,*),intmeg(3,*)
      dimension delta(4,*),lelch(*),ar(3)
      dimension lbk(*),abk(*),xxr(2)
c
      if(nps.gt.naux) stop ' error 120 '
c
      eps = 1.e-05
      tlim = 25.0
      tol = 1.e-06
c
      call rfillv(xl ,naux,0.0)
c
c *** interpolates splines for the support points.
c
      do 30 id=1,ndim
      do 10 ip=1,nps
      xl(ip) = xsp(id,ip)
   10 continue
      call spline(2,nps,xl,xne,aux)
      do 20 ip=1,nps
      tsp(id,ip) = xne(ip)
   20 continue
   30 continue
c
c *** calculates the length of each segment.
c
      u1 = 0.0
      u2 = 1.0
      ssp = 0.0
      smin= 1.e+30
      rumin = 1.e+30
      ns = nps-1
c
      do 50 is=1,ns
      do 40 id=1,ndim
      id1 = id+ndim
      r1(id ) = xsp(id,is)
      r1(id1) = tsp(id,is)
      is1 = is+1
      r2(id ) = xsp(id,is1)
      r2(id1) = tsp(id,is1)
   40 continue
      call coeff(ndim,r1,r2,a1,a2,a3,a4,a5,rum)
      coef(1,is) = a1
      coef(2,is) = a2
      coef(3,is) = a3
      coef(4,is) = a4
      coef(5,is) = a5
      rumin = min(rumin,rum)
      call lengt(0,a1,a2,a3,a4,a5,u1,u2,s,eps)
      al(is) = s
      smin= min(smin,s)
      ssp = ssp+s
   50 continue
      call findel(npbg,nebg,coorg,ieleg,intmeg,xsp(1,1),xsp(2,1),
     *             ilast,ar,i1,i2,i3,ienr,lelch)
      smin=min(smin,delta(1,ieleg(1,ienr)),delta(1,ieleg(2,ienr)),
     &        delta(1,ieleg(3,ienr)))
      xxr(1) = xsp(1,1)
      xxr(2) = xsp(2,1)
      if(lbk(1).ne.0) then
        do ip = 1,lbk(1)
         i1   = lbk(3)+(ip-1)*5
         sr   = spapt(abk(i1),xxr)
         smin = min(smin,sr,abk(i1+2))
        end do
      endif
c
      if(lbk(2).ne.0) then
        do ip = 1,lbk(2)
         i1   = lbk(4)+(ip-1)*10
         i2   = i1+5
         sr   = spaln(abk(i1),abk(i2),xxr)
         smin = min(smin,sr,abk(i1+2),abk(i2+2))
        end do
      endif
      do 56 ip=2,nps
       call findel(npbg,nebg,coorg,ieleg,intmeg,xsp(1,ip),xsp(2,ip),
     *             ilast,ar,i1,i2,i3,ienr,lelch)
       xxr(1) = xsp(1,ip)
       xxr(2) = xsp(2,ip)
       if(lbk(1).ne.0) then
         do ic = 1,lbk(1)
          i1   = lbk(3)+(ic-1)*5
          sr   = spapt(abk(i1),xxr)
          smin = min(smin,sr)
         end do
       endif
c
       if(lbk(2).ne.0) then
         do ic = 1,lbk(2)
          i1   = lbk(4)+(ic-1)*10
          i2   = i1+5
          sr   = spaln(abk(i1),abk(i2),xxr)
          smin = min(smin,sr)
         end do
       endif

       do 57 ie =1, nebg
        if(lelch(ie).ne.0)smin=min(smin,delta(1,ieleg(1,ie)),
     &             delta(1,ieleg(2,ie)),delta(1,ieleg(3,ie)))
  57   continue
  56  continue
      spmin=smin
c
c *** determine the number of information points (min = 3) from
c     the minimum spacing spmin in the B.G.
c
      nin = max(int(ssp/spmin+0.5),2)
      npi = nin+1
      if(npi.gt.naux) stop ' error 130 '
      xr = ssp/float(nin)
c
c *** places the points.
c
      ik = 1
      np = 1
      do 60 id = 1,ndim
      xip(id,np) = xsp(id,1)
      tip(id,np) = tsp(id,1)
   60 continue
c
      s1=xr
      do 100 is=1,ns
      u1 = 0.0
      u0 = is-1
      a1 = coef(1,is)
      a2 = coef(2,is)
      a3 = coef(3,is)
      a4 = coef(4,is)
      a5 = coef(5,is)
      do 70 id=1,ndim
      id1 = id+ndim
      r1(id ) = xsp(id,is)
      r1(id1) = tsp(id,is)
      is1 = is+1
      r2(id ) = xsp(id,is1)
      r2(id1) = tsp(id,is1)
   70 continue
      sis = al(is)
   80 continue
c
      if(s1.gt.sis)then
       s1 = s1-sis
       goto 100
      else
       call markp(a1,a2,a3,a4,a5,rumin,u1,u2,s1)
       call fgcurv(1,ndim,r1,r2,u2,r)
       np = np+1
       do 90 id=1,ndim
       xip(id,np) = r(id)
       tip(id,np) = r(id+ndim)
   90  continue
       ik = ik+1
       s1 = s1+xr
       if(s1.gt.eps) goto 80
      endif
  100 continue
      if(abs(npi-np).gt.1) stop ' error 140 '
      do 110 id=1,ndim
      xip(id,npi) = xsp(id,nps)
      tip(id,npi) = tsp(id,nps)
  110 continue
c
c *** computes the vector that contains the interpolated spacings.
c
c ...                              first point (k=0)
      tx = tip(1,1)
      ty = tip(2,1)
c
      ilast=1
      inorm=1
      call findel(npbg,nebg,coorg,ieleg,intmeg,xip(1,1),xip(2,1),
     *             ilast,ar,i1,i2,i3,ienr,lelch)
      call getval(npbg,i1,i2,i3,ar,delta,del,alp,anx,any,inorm)
       xxr(1) = xip(1,1)
       xxr(2) = xip(2,1)
      if(lbk(1).ne.0) then
        do 117 ip = 1,lbk(1)
        i1   = lbk(3)+(ip-1)*5
        sr   = spapt(abk(i1),xxr)
        del = min(del,sr)    
  117   continue
      endif    
c
      if(lbk(2).ne.0) then
        do 127 ip = 1,lbk(2)
        i1   = lbk(4)+(ip-1)*10
        i2   = i1+5
        sr   = spaln(abk(i1),abk(i2),xxr)
        del = min(del,sr)    
  127   continue
      endif    
c
c     if(ibl.eq.1.and.alp.gt.1.5) alp = 1.5
      call trans(del,alp,anx,any,a11,a12,a21,a22)
c
      ttx = a11*tx+a12*ty
      tty = a21*tx+a22*ty
      tm  = tx*tx+ty*ty
      ttm = ttx*ttx+tty*tty
      spa = sqrt(tm/ttm)
      sip(1) = 1./spa
c ...                            last point (k=0)
      tx = tip(1,npi)
      ty = tip(2,npi)
c
      call findel(npbg,nebg,coorg,ieleg,intmeg,xip(1,npi),xip(2,npi),
     *             ilast,ar,i1,i2,i3,ienr,lelch)
      call getval(npbg,i1,i2,i3,ar,delta,del,alp,anx,any,inorm)
       xxr(1) = xip(1,npi)
       xxr(2) = xip(2,npi)
      if(lbk(1).ne.0) then
        do 217 ip = 1,lbk(1)
        i1   = lbk(3)+(ip-1)*5
        sr   = spapt(abk(i1),xxr)
        del = min(del,sr)    
  217   continue
      endif    
c
      if(lbk(2).ne.0) then
        do 227 ip = 1,lbk(2)
        i1   = lbk(4)+(ip-1)*10
        i2   = i1+5
        sr   = spaln(abk(i1),abk(i2),xxr)
        del = min(del,sr)    
  227   continue
      endif    
c
c     if(ibl.eq.1.and.alp.gt.1.5) alp = 1.5
      call trans(del,alp,anx,any,a11,a12,a21,a22)
c
      ttx = a11*tx+a12*ty
      tty = a21*tx+a22*ty
      tm  = tx*tx+ty*ty
      ttm = ttx*ttx+tty*tty
      spa = sqrt(tm/ttm)
      sip(npi) = 1./spa
c ...                          intermediate points.
      do 120 ip=2,npi-1
      ip1 = ip+1
      x1 = xip(1,ip)
      y1 = xip(2,ip)
      tx1 = tip(1,ip)
      ty1 = tip(2,ip)
      x2 = xip(1,ip1)
      y2 = xip(2,ip1)
      tx2 = tip(1,ip1)
      ty2 = tip(2,ip1)
c
      d2x = 6.*(x2-x1)-2.*(tx2+2.*tx1)
      d2y = 6.*(y2-y1)-2.*(ty2+2.*ty1)
      vpz = tx1*d2y-ty1*d2x
      vpm = sqrt(vpz*vpz)
      vtm = sqrt(tx1*tx1+ty1*ty1)
      cur = vpm/(vtm*vtm*vtm)
c
      call findel(npbg,nebg,coorg,ieleg,intmeg,x1,y1,
     *            ilast,ar,i1,i2,i3,ienr,lelch)
      call getval(npbg,i1,i2,i3,ar,delta,del,alp,anx,any,inorm)
       xxr(1) = x1
       xxr(2) = y1
      if(lbk(1).ne.0) then
        do 317 ipp = 1,lbk(1)
        i1   = lbk(3)+(ipp-1)*5
        sr   = spapt(abk(i1),xxr)
        del = min(del,sr)    
  317   continue
      endif    
c
      if(lbk(2).ne.0) then
        do 327 ipp = 1,lbk(2)
        i1   = lbk(4)+(ipp-1)*10
        i2   = i1+5
        sr   = spaln(abk(i1),abk(i2),xxr)
        del = min(del,sr)    
  327   continue
      endif    
c
c     if(ibl.eq.1.and.alp.gt.1.5) alp = 1.5
      call trans(del,alp,anx,any,a11,a12,a21,a22)
c
      ttx = a11*tx+a12*ty
      tty = a21*tx+a22*ty
      ttz = a31*tx+a32*ty
      tm  = tx*tx+ty*ty
      ttm = ttx*ttx+tty*tty
      spa = sqrt(tm/ttm)
      if(cur*spa.gt.tlim) spa=tlim/cur
      sip(ip) = 1./spa
  120 continue
c
c *** calculates the number of elements.
c
      xne(1) = 0.0
      ane = 0.0
      do 130 is=1,npi-1
      ane = ane+0.5*(sip(is)+sip(is+1))*xr
      xne(is+1) = ane
  130 continue
      nel = max(int(ane+0.5),2)
      sc1 = ane/real(nel)
      npn = nel+1
c
c *** calculate spacings xl.
c
      ik = 0
      x1 = 0.0
      ik = ik+1
      xl(ik) = 0.0
      an = real(ik)*sc1
      do 160 is=1,npi-1
      an1 = xne(is)
      an2 = xne(is+1)
      x2 = x1+xr
      if(an.gt.an2) goto 150
      d1 = sip(is)
      d2 = sip(is+1)
      a = 0.5*(d2-d1)
      b = d1*x2-d2*x1
      f = 0.5*(x1*x1*(d1+d2)-2.*d1*x2*x1)
  140 continue
      ane = an-an1
      c = f-xr*ane
      if(abs(a/d2).lt.tol) then
       x = -c/b
      else
       d = sqrt(b*b-4.*a*c)
       x4 = 0.5*(-b+d)/a
       x5 = 0.5*(-b-d)/a
       x = x4
       if(x5.ge.x1.and.x5.le.x2) x = x5
      endif
      ik = ik+1
      xl(ik) = x
      an = real(ik)*sc1
      if(an.le.an2) goto 140
  150 continue
      x1 = x2
  160 continue
      xl(npn) = ssp
c
c *** places the points.
c
      ik = 1
      s1 = xl(ik+1)-xl(ik)
      np = 1
      do 170 id = 1,ndim
      xnp(id,np) = xsp(id,1)
  170 continue
c
      do 210 is=1,ns
      u1 = 0.0
      u0 = is-1
      a1 = coef(1,is)
      a2 = coef(2,is)
      a3 = coef(3,is)
      a4 = coef(4,is)
      a5 = coef(5,is)
      do 180 id=1,ndim
      id1 = id+ndim
      r1(id ) = xsp(id,is)
      r1(id1) = tsp(id,is)
      is1 = is+1
      r2(id ) = xsp(id,is1)
      r2(id1) = tsp(id,is1)
  180 continue
      sis = al(is)
  190 continue
c
      if(s1.gt.sis) then
       s1 = s1-sis
       goto 210
      else
       call markp(a1,a2,a3,a4,a5,rumin,u1,u2,s1)
       call fgcurv(0,ndim,r1,r2,u2,r)
       np = np+1
       do 200 id=1,ndim
       xnp(id,np) = r(id)
  200  continue
       ik = ik+1
       s1 = s1+xl(ik+1)-xl(ik)
       if(s1.gt.eps) goto 190
      endif
  210 continue
      if(abs(npn-np).gt.1) stop ' error 150 '
      do 220 id=1,ndim
      xnp(id,npn) = xsp(id,nps)
  220 continue
c
      return
      end
c
c ---------------------------------------------------------------------------
c
c *** this sub. finds the tangents in the points defining a ferguson splines.
c     ib is an indicator of the geometrical constraints at the ends.
c
c ----------------------------------------------------------------------------
c
      subroutine spline(ib,n,r,t,aux)
c             ib  =  1 ............ specified tangents: t(1),t(n).
c             ib  =  2 ............ zero second derivatives.
c
c     note: the number of points is n. the tridiagonal system of n-2
c           equations is solved by gauss elimination & backsubstitution.
c
      dimension aux(*),r(*),t(*)
c
c *** first row i = 2.   ib = 1  -> [ 4, 1] ;  ib = 2  -> [ 3.5, 1]
c
      if(n.le.1) then
       pause ' spline: n=1 wrong number of points'
      else if(n.eq.2) then
       t(1) = r(2)-r(1)
       t(2) = t(1)
      else if(n.eq.3) then
       if(ib.eq.1) then
        t(2) = 0.25*(3.*(r(3)-r(1))-t(1)-t(3))
       else
        t(1) = -1.25*r(1)+1.50*r(2)-0.25*r(3)
        t(2) = -0.50*r(1)          +0.50*r(3)
        t(3) =  0.25*r(1)-1.50*r(2)+1.25*r(3)
       endif
      else
       rv = 3.*(r(3)-r(1))
       if(ib.eq.1) then
        bet = 4.0
        rv = rv-t(1)
       else
        bet = 3.5
        rv = rv-1.5*(r(2)-r(1))
       endif
       t(2) = rv/bet
c
c *** rows of the type  [ 1,  4,  1 ]
c
       do 100 j=3,n-2
       aux(j) = 1./bet
       bet = 4.-aux(j)
       if(bet.eq.0.) pause ' error in spline: zero pivot'
       rv = 3.*(r(j+1)-r(j-1))
       t(j) = (rv-t(j-1))/bet
  100  continue
c
c *** last row i = n-1   ib = 1  -> [ 1, 4] ;  ib = 2  -> [ 1, 3.5 ].
c
       aux(n-1) = 1./bet
       rv = 3.*(r(n)-r(n-2))
       if(ib.eq.1) then
        bet = 4.
        rv = rv-t(n)
       else
        bet = 3.5
        rv = rv-1.5*(r(n)-r(n-1))
       endif
       bet = bet-aux(n-1)
       if(bet.eq.0.) pause ' error in spline: zero pivot'
       t(n-1) = (rv-t(n-2))/bet
c
c *** backsubstitution.
c
       do 200 j=n-2,2,-1
       t(j) = t(j)-aux(j+1)*t(j+1)
  200  continue
c
c *** end values when ib = 2.
c
       if(ib.eq.2) then
        t(1) = 1.5*(r(2)-r(1  ))-0.5*t(2)
        t(n) = 1.5*(r(n)-r(n-1))-0.5*t(n-1)
       endif
      endif
c
      return
      end
c
c ----------------------------------------------------------------------
c
c *** this sub. computes the coeficients of the polynomial |r'|**2
c
c ------------------------------------------------------------------------
c
      subroutine coeff(ndim,r1,r2,a1,a2,a3,a4,a5,rumin)
c
      dimension r1(2*ndim),r2(2*ndim)
c
      a1 = 0.0
      a2 = 0.0
      a3 = 0.0
      a4 = 0.0
      a5 = 0.0
      am1 = 0.
      am2 = 0.
c
      do 1000 id=1,ndim
      id1 = id+ndim
      r12 = r1(id)-r2(id)
      p = r1(id1)
      q = r2(id1)
      pp = p*p
      qq = q*q
      s = 3.*( 2.*r12+   p+q)
      t = 2.*(-3.*r12-2.*p-q)
      a1 = a1+pp
      a2 = a2+2.*p*t
      a3 = a3+t*t+2.*p*s
      a4 = a4+2.*t*s
      a5 = a5+s*s
      am1 = am1+pp
      am2 = am2+qq
 1000 continue
c
      rumin = min(am1,am2)
      rumin = sqrt(rumin)
c
      return
      end
c
c -----------------------------------------------------------------------------
c
c *** this sub. computes the length of a segment of a cubic.
c
c ------------------------------------------------------------------------------
c
      subroutine lengt(in,a1,a2,a3,a4,a5,u1,u2,s,eps)
c
      parameter(nit=40)
c
      eps1 = eps
      os = -1.e+30
c
      f1 = sqrt(a1+u1*(a2+u1*(a3+u1*(a4+u1*a5))))
      f2 = sqrt(a1+u2*(a2+u2*(a3+u2*(a4+u2*a5))))
c
      u21 = u2-u1
      st = 0.5*u21*(f1+f2)
      ost = st
      kt = 1
      do 200 it=1,nit
      tnm = kt
      del = u21/tnm
      x = u1+0.5*del
      sum = 0.0
      do 100 jt=1,kt
      sum = sum+sqrt(a1+x*(a2+x*(a3+x*(a4+x*a5))))
      x = x+del
  100 continue
      st = 0.5*(st+u21*sum/tnm)
      kt = kt*2
      s = (4.*st-ost)/3.
      if(in.eq.0) eps1 = eps*abs(os)
      if(abs(s-os).le.eps1) goto 300
      os = s
      ost = st
  200 continue
      stop ' error 160 '
  300 return
      end
c
c -----------------------------------------------------------------------------
c
c *** this sub. calculates the position u2 of a point on a f.s. such that
c     the length of the cubic segment u1,u2 is sl.
c
c ------------------------------------------------------------------------------
c
      subroutine markp(a1,a2,a3,a4,a5,rumin,u1,u2,sl)
c
      parameter(nit=20)
c
      eps1 = 1.e-03
      eps2 = max(rumin*eps1*1.e-02,1.e-05)
      ux = u1
      u2 = u1
      ss = sl
      do 100 it=1,nit
      ff = sqrt(a1+u2*(a2+u2*(a3+u2*(a4+u2*a5))))
      u2 = u2+ss/ff
      if(abs(u2-ux).le.eps1) goto 200
      call lengt(1,a1,a2,a3,a4,a5,ux,u2,vi,eps2)
      ss = ss-vi
      ux = u2
  100 continue
      stop ' error 160 '
  200 return
      end
c
c ----------------------------------------------------------------------------
c
c *** to evaluate the transformation matrix from an ellipse to an
c     unstretched circle of unit radius
c
c ------------------------------------------------------------------------------
c
      subroutine trans(d,s,a1,a2,t11,t12,t21,t22)
c
      an=sqrt(a1*a1+a2*a2)
      an1=a2/an
      an2=-a1/an
c
      t11=(a1*a1/s+an1*an1)/d
      t12=(a1*a2/s+an1*an2)/d
      t21=t12
      t22=(a2*a2/s+an2*an2)/d
c
      return
      end
c 
c --------------------------------------------------------------------------
c
c *** expression of a ferguson curve segment.
c
c ----------------------------------------------------------------------------
c
      subroutine fgcurv(ider,ndim,r1,r2,u,r)
c
      dimension r1(3*ndim),r2(3*ndim),r(3*ndim)
c
      do 1000 id=1,ndim
c
      id1 = id+ndim
      id2 = id1+ndim
      r12 = r2(id)-r1(id)
      a1 = r1(id)
      a2 = r1(id1)
      a3 =  3.*r12-2.*r1(id1)-r2(id1)
      a4 = -2.*r12+   r1(id1)+r2(id1)
      r(id) = a1+u*(a2+u*(a3+u*a4))
      if(ider.eq.0) goto 1000
      r(id1) = a2+u*(2.*a3+3.*u*a4)
      if(ider.eq.1) goto 1000
      r(id2) = 2.*a3+6.*u*a4
c
 1000 continue
c
      return
      end
c
c-------------------------------------------------------------------
c
c     this subroutine finds out the boundary points
c     lcoor(ipoin).eq.0 --> interior
c     lcoor(ipoin).ne.0 --> boundary
c
c--------------------------------------------------------------------
c
      subroutine bound(npoin,nelem,lcoor,intmat)
      dimension lcoor(*),intmat(3,*)
      do 1000 ip=1,npoin
        lcoor(ip)=0
 1000 continue
c
c     loop over the elements
c
      do 2000 ie=1,nelem
        do 2001 in=1,3
          in1=in+1
          if(in1.gt.3) in1=in1-3
          in2=in+2
          if(in2.gt.3) in2=in2-3
          ip=intmat(in,ie)
          lcoor(ip)=lcoor(ip)+intmat(in2,ie)-intmat(in1,ie)
2001    continue
2000  continue
      return
      end
c
c--------------------------------------------------------------------
c
c
c--------------------------------------------------------------------
c
      subroutine conere(npoin,nelem,intmat,lcore)
      dimension intmat(3,*),lcore(*)
      do 1000 ip=1,npoin
      lcore(ip)=0
 1000 continue
      do 2000 ie=1,nelem
      do 2001 in=1,3
      ip=intmat(in,ie)
      lcore(ip)=lcore(ip)+1
 2001 continue
 2000 continue
      return
      end
c
c--------------------------------------------------------------------
c
c     this subr. finds out the optimal number of conectivities
c     for each node
c
c--------------------------------------------------------------------
c
      subroutine coneid(npoin,nelem,intmat,lcoor,lcoid,coord,strec)
      dimension intmat(3,*),lcoor(*),lcoid(*)
      dimension coord(2,*),strec(4,*)
      xba(xq,yq,ax,ay,alph)=(ax*(ax*xq+ay*yq)/alph)+ay*(ay*xq-ax*yq)
      yba(xq,yq,ax,ay,alph)=(ay*(ax*xq+ay*yq)/alph)-ax*(ay*xq-ax*yq)
      pi=3.1415927
      do 1000 ip=1,npoin
      lcoid(ip)=6
 1000 continue
c
c     search for a side to start
c
      do 2000 ip=1,npoin
      ic=ip
      if(lcoor(ip).ne.0) goto 2001
 2000 continue
 2001 continue
      ia=0
      ib=0
      do 3000 ie=1,nelem
      do 3001 in=1,3
      ip=intmat(in,ie)
      if(ip.ne.ic) goto 3001
      inb=in-1
      if(inb.lt.1) inb=inb+3
      if(lcoor(intmat(inb,ie)).ne.0) ib=intmat(inb,ie)
      if(ib.ne.0) goto 3002
 3001 continue
 3000 continue
 3002 continue
      imemo=ic
 4000 ia=ib-lcoor(ic)
      alph=strec(2,ic)
      if(alph.le.0.00001) alph=1
      anx=strec(3,ic)
      any=strec(4,ic)
      x1r=coord(1,ib)-coord(1,ic)
      y1r=coord(2,ib)-coord(2,ic)
      x2r=coord(1,ia)-coord(1,ic)
      y2r=coord(2,ia)-coord(2,ic)
      x1=xba(x1r,y1r,anx,any,alph)
      y1=yba(x1r,y1r,anx,any,alph)
      x2=xba(x2r,y2r,anx,any,alph)
      y2=yba(x2r,y2r,anx,any,alph)
      bot=(sqrt(x1*x1+y1*y1)*sqrt(x2*x2+y2*y2))
      if(bot.le..000001) bot=1
      cosa=(x1*x2+y1*y2)/bot
      if(cosa.gt.1.0) cosa=1.0
      if(cosa.lt.-1.0) cosa=-1.0
      theta=acos(cosa)
      if((x2*y1-x1*y2).lt.0.0) theta=2.*pi-theta
      divi=(theta+(pi/6.))/(pi/3.)
      lcoid(ic)=divi
      if(lcoid(ic).lt.1) lcoid(ic)=1
      ib=ic
      ic=ia
      if(ic.ne.imemo) goto 4000
c
      return
      end
c
c--------------------------------------------------------------------
c
c
c--------------------------------------------------------------------
c
      subroutine side(nelem,npoin,iloca,intmat,iside,lwher,lhowm,
     1                icone,ind)
      dimension intmat(3,*),iside(4,*)
      dimension lwher(*),lhowm(*),icone(*)
c
c     fill in lhowm : nr. of elements per node
c
      do 1490 ip=1,npoin
      lhowm(ip)=0
 1490 continue
      do 1500 ie=1,nelem
      do 1500 in=1,3
      ip=intmat(in,ie)
      lhowm(ip)=lhowm(ip)+1
 1500 continue
 1501 continue
c
c      fill in lwher : location of each node inside icone
c
      lwher(1)=0
      do 1600 ip=2,npoin
      lwher(ip)=lwher(ip-1)+lhowm(ip-1)
 1600 continue
c
c      fill in icone : elements in each node
c
      do 1690 ip=1,npoin
      lhowm(ip)=0
 1690 continue
      do 1700 ie=1,nelem
      do 1701 in=1,3
      ip=intmat(in,ie)
      lhowm(ip)=lhowm(ip)+1
      jloca=lwher(ip)+lhowm(ip)
      icone(jloca)=ie
 1701 continue
 1700 continue
      if(ind.eq.1) return
c
c     loop over the nodes
c
      iloca=0
      do 3000 ip=1,npoin
      iloc1=iloca
      iele=lhowm(ip)
      if(iele.eq.0) goto 3000
c
c     initialize iside ----> important for boundary sides
c
      do 3001 is=1,iele+2
      iside(3,is+iloc1)=0
      iside(4,is+iloc1)=0
 3001 continue
      iwher=lwher(ip)
c
c     loop over elements surrounding the point ip
c
      ip1=ip
      do 3090 iel=1,iele
      ie=icone(iwher+iel)
c
c     find out position of ip in the coneivity matrix
c
      do 3091 in=1,3
      in1=in
      ipt=intmat(in,ie)
      if(ipt.eq.ip) goto 3092
 3091 continue
 3092 continue
      do 3100 j=1,2
      in2=in1+j
      if(in2.gt.3) in2=in2-3
      ip2=intmat(in2,ie)
      if(ip2.lt.ip1) goto 3100
c
c     check the side ----->  new or old
c
      if(iloca.eq.iloc1) goto 7304
      do 5600 is=iloc1+1,iloca
      jloca=is
      if(iside(2,is).eq.ip2) goto 7303
 5600 continue
 7304 continue
c
c     new side
c
      iloca=iloca+1
      iside(1,iloca)=ip1
      iside(2,iloca)=ip2
      iside(2+j,iloca)=ie
      goto 3012
c
c     old side
c
 7303 continue
      iside(2+j,jloca)=ie
 3012 continue
 3100 continue
c
c     end loop over elements surrounding point ip
c
 3090 continue
      do 8000 is=iloc1+1,iloca
      if(iside(3,is).ne.0) goto 8000
      iside(3,is)=iside(4,is)
      iside(4,is)=0
      iside(1,is)=iside(2,is)
      iside(2,is)=ip1
 8000 continue
c
c     end loop over points
c
 3000 continue
      return
      end
c
c--------------------------------------------------------------------
c
c      this subroutine checks the validity of the mesh
c
c--------------------------------------------------------------------
c
      subroutine areach(npoin,nelem,intmat,coord)   
      dimension coord(2,*),intmat(3,*)
      character*80 st
c
c     determinant function
c
c     deter(p1,q1,p2,q2,p3,q3)=p2*q3-p3*q2-p1*q3+p3*q1+p1*q2-p2*q1 
      deter(p1,q1,p2,q2,p3,q3)=(p2-p1)*(q3-q1)-(p3-p1)*(q2-q1)
      isucc=0
      do 1000 ie=1,nelem
        i1=intmat(1,ie)
        i2=intmat(2,ie)
        i3=intmat(3,ie)
 
        x1=coord(1,i1)
        x2=coord(1,i2)
        x3=coord(1,i3)
        y1=coord(2,i1)
        y2=coord(2,i2)
        y3=coord(2,i3)
        if(deter(x1,y1,x2,y2,x3,y3).gt.0.0) then
          goto 1000
        else
          write (2, 110)ie,i1,x1,y1,i2,x2,y2,i3,x3,y3
          isucc=1
        end if
1000  continue
      if(isucc.eq.0) then
        st='           the checking has been successful.'
        call plot(2,x,y,st ,xmin,xmax,ymin,ymax)
       else
        st='      the checking has not been successful.'
        call plot(2,x,y,st ,xmin,xmax,ymin,ymax)
      end if 
      return
110   format('  area error in element ',i5,/,
     *       '  node 1 =',i5,'  x1 =',f10.5,'  y1 =',f10.5,/,
     *       '  node 2 =',i5,'  x2 =',f10.5,'  y2 =',f10.5,/,
     *       '  node 3 =',i5,'  x3 =',f10.5,'  y3 =',f10.5,/)
      end
c
c--------------------------------------------------------------------
c
c     this subroutine fills in the element conectivity matrix needed
c     for the fast searching algorithm
c
c-------------------------------------------------------------------
c
      subroutine elecon(nside,intmat,iside,intmel)
      dimension intmat(3,*),iside(4,*),intmel(3,*)
c
c     loop over the sides
c
      do 1000 is=1,nside
        ip1=iside(1,is)
        ie1=iside(3,is)
        ie2=iside(4,is)
c
c     first ie1
c
        i1=intmat(1,ie1)
        i2=intmat(2,ie1)
        i3=intmat(3,ie1)
        if(ip1.eq.i1) ipos=3
        if(ip1.eq.i2) ipos=1
        if(ip1.eq.i3) ipos=2
c
c     go into intmel
c
        intmel(ipos,ie1)=ie2
c
c     now ie2 if ie2.ne.0
c
        if(ie2.eq.0) goto 1000
        i1=intmat(1,ie2)
        i2=intmat(2,ie2)
        i3=intmat(3,ie2)
        if(ip1.eq.i1) ipos=2
        if(ip1.eq.i2) ipos=3
        if(ip1.eq.i3) ipos=1
c
c     go into intmel
c
        intmel(ipos,ie2)=ie1
c 
c     end loop over the sides
c
1000  continue
      return
      end
c
c--------------------------------------------------------------------
c
c     this subroutine adds triangles to the background grid in order
c     to make a convex region
c
c--------------------------------------------------------------------
c
      subroutine convex(nelem ,npoin ,intmat ,coord ,npfrt ,nqfrt ,
     *                  lcoor ,nregi )
      dimension coord(2,*) ,intmat(3,*) ,lcoor(*)
      dimension npfrt(*)   ,nqfrt(*)    ,nregi(*)
c
c *** first find out the boundary points
c
      call  bound(npoin ,nelem ,lcoor ,intmat )
c
c *** store boundary points
c
      nside=0
      do 1000 ip=1,npoin
        if(lcoor(ip).eq.0) goto 1000
        nside=nside+1
        nregi(nside)=ip
1000  continue
c
c *** now set up the front
c *** identify a side
c
      iside=0
      iel=0
1999  iel=iel+1
      ihow=0
      do 2001 in=1,3
        ip=intmat(in,iel)
        if(lcoor(ip).ne.0) ihow=ihow+1
2001  continue
      if(ihow.ge.2) goto 2002
2000  goto 1999
2002  continue
c
c *** get the first side
c
      do 3000 in=1,3
        ip=intmat(in,iel)
        if(lcoor(ip).eq.0) goto 3000
        in1=in+1
        if(in1.eq.4) in1=1
        ip1=intmat(in1,iel)
        if(lcoor(ip1).eq.0) goto 3000
c
c *** check for the next four sides
c
        ip2=ip1+lcoor(ip)
        if(ip2.le.0) goto 3000
        if(lcoor(ip2).eq.0) goto 3000
        ip3=ip+lcoor(ip2)
        if(ip3.le.0) goto 3000
        if(lcoor(ip3).eq.0) goto 3000
        ip4=ip2+lcoor(ip3)
        if(ip4.le.0) goto 3000
        if(lcoor(ip4).eq.0) goto 3000
        ip5=ip3+lcoor(ip4)
        if(ip5.le.0) goto 3000
        if(lcoor(ip5).eq.0) goto 3000
        goto 3001
3000  continue
      goto 1999
3001  continue
      iold=ip1
c
c *** the first side is ip1-ip
c     we can start filling the front
c
4000  iside=iside+1
      npfrt(iside)=ip1
      nqfrt(iside)=ip
      ikeep=ip
      if(lcoor(ip).eq.0) goto 4010
      ip=ip1+lcoor(ip)
      if(ip.le.0) goto 4010
      goto 4011
 4010 stop 'error 400'
 4011 continue
      lcoor(ikeep)=0
      ip1=ikeep
      if(ip1.ne.iold) goto 4000
c
c *** check if iside.eq.nside
c
      if(iside.eq.nside) goto 4500
c
c *** there are more holes, find them
c
      goto 1999
 4500 continue
c
c *** loop over the number of sides
c
      npbou=nside
      iposi=0
 6000 iposi=iposi+1
      kn1=npfrt(iposi)
      if(kn1.eq.0) goto 5900
      kn=nqfrt(iposi)
      xn=coord(1,kn)
      yn=coord(2,kn)
      xn1=coord(1,kn1)
      yn1=coord(2,kn1)
      a=yn1-yn
      b=xn-xn1
      c=(xn1-xn)*yn1+(yn-yn1)*xn1
      inum=0
      ang1=0.0
      do 6500 is=1,npbou
      ken=nregi(is)
      if(ken.eq.kn.or.ken.eq.kn1) goto 6500
      xken=coord(1,ken)
      yken=coord(2,ken)
      if(a*xken+b*yken+c.le.0.0) goto 6500
      xdif1=xken-xn1
      ydif1=yken-yn1
      xdif2=xken-xn
      ydif2=yken-yn
      dist1=sqrt(xdif1*xdif1+ydif1*ydif1)
      dist2=sqrt(xdif2*xdif2+ydif2*ydif2)
      cosa=(xdif1*xdif2+ydif1*ydif2)/(dist1*dist2)
      if(cosa.gt.1.0) cosa=1.0
      if(cosa.lt.-1.0) cosa=-1.0
      angl=acos(cosa)
      if(angl.lt.0.1) goto 6500
c
c *** see if it is a possible connectivity
c
      call possib(kn1   ,kn    ,ken   ,xn1   ,yn1   ,xn    ,yn    ,     
     *            xken  ,yken  ,nside ,npbou ,npfrt ,nqfrt ,nregi ,     
     *            coord ,iyon  )
c
      if(iyon.eq.0) goto 6500
      inum=1
      if(angl.lt.ang1) goto 6500
      ang1=angl
      kpo=ken
 6500 continue
c
      if(inum.eq.0) goto 5900
c
c *** create a new element
c
      nelem=nelem+1
      intmat(1,nelem)=kn1
      intmat(2,nelem)=kn
      intmat(3,nelem)=kpo
c
c *** update the front
c
      npfrt(iposi)=0
      nqfrt(iposi)=0
      ind=0
      do 5700 is=1,nside
        if(npfrt(is).ne.kpo.or.nqfrt(is).ne.kn1) goto 5700
        ind=1
        npfrt(is)=0
        nqfrt(is)=0
        goto 5701
 5700 continue
 5701 continue
      if(ind.eq.1) goto 5790
      nside=nside+1
      npfrt(nside)=kn1
      nqfrt(nside)=kpo
 5790 ind=0
      do 5800 is=1,nside
        if(npfrt(is).ne.kn.or.nqfrt(is).ne.kpo) goto 5800
        ind=1
        npfrt(is)=0
        nqfrt(is)=0
        goto 5801
 5800 continue 
 5801 continue
      if(ind.eq.1) goto 5900
      nside=nside+1
      npfrt(nside)=kpo
      nqfrt(nside)=kn
 5900 if(iposi.lt.nside) goto 6000
      return
      end
c
c--------------------------------------------------------------------
c
c     this subroutine modifies the coordinates of the background
c     grid boundary  points  in  order  to  completely cover the
c     domain of interest
c
c--------------------------------------------------------------------
c
      subroutine expand(neleg,npoig,coorg,coor0,ieleg,intmeg)
      dimension coorg(2,*),coor0(2,*)
      dimension  ieleg(3,*),intmeg(3,*)
      data expa/0.1/
      call rfillv(coor0,2*npoig,0.0)
c
c *** loop over the elements
c
      ae=1.e+6
      do 1000 ie=1,neleg
        do 1001 in=1,3
          if(intmeg(in,ie).ne.0) goto 1001
c
c *** move nodes
c
          in1=in+1
          in2=in+2
          if(in1.gt.3) in1=in1-3
          if(in2.gt.3) in2=in2-3
          ip1=ieleg(in1,ie)
          ip2=ieleg(in2,ie)
          ax=coorg(2,ip2)-coorg(2,ip1)
          ay=coorg(1,ip1)-coorg(1,ip2)
          at=sqrt(ax*ax+ay*ay)
          if(at.lt.ae) ae=at
1001    continue
1000  continue

      exp=expa*ae
      do 2000 ie=1,neleg
        do 2001 in=1,3
          if(intmeg(in,ie).ne.0) goto 2001
c
c *** move nodes
c
          in1=in+1
          in2=in+2
          if(in1.gt.3) in1=in1-3
          if(in2.gt.3) in2=in2-3
          ip1=ieleg(in1,ie)
          ip2=ieleg(in2,ie)
          ax=coorg(2,ip2)-coorg(2,ip1)
          ay=coorg(1,ip1)-coorg(1,ip2)
          at=sqrt(ax*ax+ay*ay)
          ax=ax*exp/at
          ay=ay*exp/at
          coor0(1,ip1)=coor0(1,ip1)+ax
          coor0(2,ip1)=coor0(2,ip1)+ay
          coor0(1,ip2)=coor0(1,ip2)+ax
          coor0(2,ip2)=coor0(2,ip2)+ay
2001    continue
2000  continue
c
c *** add
c
      do 3000 i=1,npoig
        do 3001 j=1,2
          coorg(j,i)=coorg(j,i)+coor0(j,i)
3001    continue
3000  continue
      return
      end
c
c----------------------------------------------------------------------
c
c     this subroutine generates the required parameters for re-meshing 
c     from a previously obtained solution
c
c----------------------------------------------------------------------
c
      subroutine ref2(npoig,neleg,coorg,unkng,intmag,strec,
     &                geomg,mmatg,deri2,help,derip,lhowm,lwher,icone,
     &                factor,unkey,fc)
       real mmatg
       dimension coorg(2,npoig),unkng(4,npoig)
       dimension geomg(7,neleg),mmatg(npoig),derip(2,npoig)
       dimension deri2(4,npoig),strec(4,npoig),help(4,npoig)
       dimension  intmag(3,neleg),lhowm(npoig),lwher(npoig)
       dimension icone(3*neleg),factor(npoig) ,unkey(npoig)
       dimension fc(*)
c
c      obtain globally the geometrical values needed
c
       ndimn=2
       namat=4
       nnode=3
       ngeom=7
       call getgeo(neleg,npoig,nnode,ndimn,ngeom,
     &             intmag,coorg,geomg,nstop)
c
c     initialize
c
       do 1702 ipoin=1,npoig
       factor(ipoin)=0.
       strec(1,ipoin)=1.e+6
       strec(2,ipoin)=1.e+6
 1702  continue
c
c     -----inquire for maximum and minimum spacing
c
       write(*,100)
 100   format(' enter maximum   spacing  desired  : ', $ )
       read(*,*) deltma
       write(*,200) 
 200   format(' enter minimum   spacing  desired  : ', $ )
       read(*,*) deltmi
       write(*,250) 
 250   format(' enter      scaling       desired  : ', $ )
       read(*,*) deltas
       write(*,300)
 300   format(' enter maximum stretching desired  : ', $ )
       read(*,*) str
c
c     -----iquire for refinement criteria
c
       print *,' what do you want to use as key variable ?'
       print *,' 1 - density'
       print *,' 2 - velocity (modulus)'
       print *,' 3 - density + mach'
       print *,' 4 - entropy'
       print *,' 5 - density + entropy'
       print *,' 6 - mach number'
       print *,' 7 - pressure'
       write(*,400)
 400   format(' enter option chosen   : ' ,$ )
       read(*,*) icrit
       print *,' what do you want to use as eps ?'
       read(*,*) eps
       print *,' what do you want to use as emax ?'
       read(*,*) emax
       iback=0 
       ivari=1
       if(icrit.eq.2) ivari=2
       if(icrit.eq.4.or.icrit.eq.5) ivari=3
       if(icrit.eq.6.or.icrit.eq.3) ivari=4
       if(icrit.eq.7) ivari=5
       if(icrit.eq.3.or.icrit.eq.5) iback=1
       goto 1210
 1209  if(iback.eq.1) then
       ivari=1
       iback=0
       else
       endif
 1210  continue
c
c     -----obtain the gradient at the nodes
c
       call getvar(nnode,namat,ngeom,2,neleg,ndimn,npoig,
     &             coorg,intmag,derip,deri2,geomg,
     &             mmatg,unkng,ivari)
c
c      find out principal directions and variations
c
       do 1601 ipoin=1,npoig
       amean=0.5*(deri2(1,ipoin)+deri2(4,ipoin))
       adevi=0.5*(deri2(1,ipoin)-deri2(4,ipoin))
       tauxy=deri2(2,ipoin)
       si1=amean+sqrt(adevi**2+tauxy**2)
       si2=amean-sqrt(adevi**2+tauxy**2)
       if((abs(si1)+abs(si2)).lt.1.e-8) then
         help(1,ipoin)=1.e-9
         help(2,ipoin)=1.e-9
         help(3,ipoin)=1.0
         help(4,ipoin)=0.0
       else
         if(amean.ge.0.0) then
           help(1,ipoin)=abs(si1)
           help(2,ipoin)=abs(si2)
           alpha=atan2(tauxy,adevi)/2.
           help(3,ipoin)=-sin(alpha)
           help(4,ipoin)= cos(alpha)
         else
           help(1,ipoin)=abs(si2)
           help(2,ipoin)=abs(si1)
           alpha=atan2(tauxy,adevi)/2.
           help(3,ipoin)= cos(alpha)
           help(4,ipoin)= sin(alpha)
         endif
       endif
 1601  continue
c
c      search for the maximum spacing and stretching
c
       print *,'input first point'
       read(*,*) npstr 
       amax=0.0
       do 1201 ipoin=npstr,npoig
       if(help(1,ipoin).gt.amax) amax=help(1,ipoin)
 1201  continue
       amaxsq=sqrt(amax)
       do 1202 ipoin=1,npoig
       delta=deltas*amaxsq/sqrt(max(1.e-9,help(1,ipoin)))
       delt2=deltas*amaxsq/sqrt(max(1.e-9,help(2,ipoin)))
       if(delta.gt.deltma) delta=deltma
       if(delt2.gt.deltma) delt2=deltma
       if(delta.lt.deltmi) delta=deltmi
       if(delt2.lt.deltmi) delt2=deltmi
       help(1,ipoin)=delta
       help(2,ipoin)=delt2
 1202  continue
       do 1302 ipoin=1,npoig
       delta=help(1,ipoin)
       strec(2,ipoin)=min(strec(2,ipoin),help(2,ipoin))
       if(delta.ge.strec(1,ipoin)) goto 1302
       strec(1,ipoin)=delta
       strec(3,ipoin)=help(3,ipoin)
       strec(4,ipoin)=help(4,ipoin)
 1302  continue
       call error(npoig,neleg,unkng,ivari,lhowm,lwher,icone,
     1           geomg,ngeom,nnode,intmag,fc,eps,unkey)               
       do 7252 ip = 1, npoig 
        factor(ip)=max(fc(ip),factor(ip))
 7252  continue
       if(iback.eq.1) goto 1209
       do 1402 ipoin=1,npoig
       delta=strec(1,ipoin)
       delt2=strec(2,ipoin)
       stre=delt2/delta
       if(stre.gt.str) stre=str
       if(stre.lt.1.0) stre=1.0
       strec(2,ipoin)=stre
 1402  continue
 3420  continue
c
       if(emax.eq.1.or.eps.eq.1) return
       do 56 ip=1,npoig                                                         
            if(factor(ip).ge.emax)then                                          
              strec(1,ip)=deltmi
              strec(2,ip)=str
            else
              strec(1,ip)= (deltmi-deltma)*factor(ip)/emax+deltma
              strec(2,ip)=-(1.-str)*factor(ip)/emax+1.
            endif     
  56   continue                                                                 
c
      return
      end
c
      subroutine error(npoin,nelem,unkno,ivari,lhowm,lwher,icone,
     &                 geome,ngeom,nnode,intmat,err,eps,unkey)
c
c *** calculates the dimensionless error according to lohner
c
      dimension unkno(4,npoin),lhowm(npoin),lwher(npoin)
      dimension icone(3*nelem),geome(7,nelem),intmat(3,nelem)
      dimension err(npoin),unkey(npoin)
c
c *** find the key variable
c
      gamma=1.4
      do 500 ip=1,npoin
        ro=unkno(1,ip)
        u1=unkno(2,ip)
        u2=unkno(3,ip)
        en=unkno(4,ip)
        vsq=u1*u1+u2*u2
        pre=(gamma-1.)*ro*(en-0.5*vsq)
        if(ivari.eq.1)unkey(ip)=ro
        if(ivari.eq.2)unkey(ip)=sqrt(vsq)
        if(ivari.eq.3)unkey(ip)=pre/(ro*gamma)
        if(ivari.eq.4)unkey(ip)=sqrt((vsq*ro)/(gamma*pre))
        if(ivari.eq.5)unkey(ip)=pre 
c
500   continue
c *** loop over the points
c
      do 7000 ip=1,npoin
       topa=0.0
       topb=0.0
       topc=0.0
       topd=0.0
       bota=0.0
       botb=0.0
       botc=0.0
       botd=0.0
c *** loop over the elements surrounding ip
c
       nelsur=lhowm(ip)
       nstart=lwher(ip)+1
       do 6500 ncone=nstart,nstart+nelsur-1
        etopa=0.0
        etopb=0.0
        ebote1=0.0
        ebote2=0.0
        iel=icone(ncone)
        rj=geome(ngeom,iel)*0.5
        do 6000 in=1,nnode
         etopa=etopa+unkey(intmat(in,iel))*geome(in,iel)
         etopb=etopb+unkey(intmat(in,iel))*
     &            geome(in+nnode,iel)
         ebote1=ebote1+abs(unkey(intmat(in,iel))*
     &             geome(in,iel))
         ebote2=ebote2+abs(unkey(intmat(in,iel))*
     &              geome(in+nnode,iel))
         if(intmat(in,iel).eq.ip)ix=in
6000    continue
       ebota=abs(etopa)
       ebotb=abs(etopb)
       etopc=etopa*geome(ix+nnode,iel)*rj
       etopd=etopb*geome(ix+nnode,iel)*rj
       etopa=etopa*geome(ix,iel)*rj
       etopb=etopb*geome(ix,iel)*rj
c
       ebota=ebota+eps*ebote1
       ebotc=ebota*abs(geome(ix+nnode,iel))*rj
       ebota=ebota*abs(geome(ix,iel))*rj
       ebotb=ebotb+eps*ebote2
       ebotd=ebotb*abs(geome(ix+nnode,iel))*rj
       ebotb=ebotb*abs(geome(ix,iel))*rj
       topa=topa+etopa**2
       topb=topb+etopb**2
       topc=topc+etopc**2
       topd=topd+etopd**2
       bota=bota+ebota**2
       botb=botb+ebotb**2
       botc=botc+ebotc**2
       botd=botd+ebotd**2
6500  continue
c
      topa=sqrt(topa+topb+topc+topd)
      bota=sqrt(bota+botb+botc+botd)
c
      err(ip)=(topa/bota)
7000   continue
c
      return      
      end
c
c---------------------------------------------------------------------
c
c     this subroutine evaluates n,x & n,y for each element
c     shape function (linear triangles) and the jacobian (2=2*area)
c
c---------------------------------------------------------------------
c
       subroutine getgeo(nelem,npoin,nnode,ndimn,ngeom,
     &                   intmat,coord,geome,nstop)
       dimension    coord(ndimn,npoin),geome(ngeom,nelem)
       dimension    x(3),y(3),nxi(3),net(3)
       dimension  intmat(nnode,nelem)
       data  nxi/-1.0 , 1.0 , 0.0 /
       data  net/-1.0 , 0.0 , 1.0 /
       do 1000 ielem=1,nelem
c
c      pick up the values needed
c
       do 1001 inode=1,nnode
       in=intmat(inode,ielem)
       x(inode)=coord(1,in)
       y(inode)=coord(2,in)
 1001 continue
c
c      evaluate the geometrical quantities needed
c
       x21=x(2)-x(1)
       x31=x(3)-x(1)
       y21=y(2)-y(1)
       y31=y(3)-y(1)
       rj =x21*y31-x31*y21
       if(rj.gt.0.0) goto 2030
       write(6,12) ielem
 12    format('  some element have got negative area ',i5)
       nstop=1
       go to 1000
 2030  continue
       rj1=1./rj
       xix= y31*rj1
       xiy=-x31*rj1
       etx=-y21*rj1
       ety= x21*rj1
c
c      form n,x & n,y
c
       do 1002 in=1,3
       rnxi=nxi(in)
       rnet=net(in)
       geome(in  ,ielem)=xix*rnxi+etx*rnet
       geome(in+3,ielem)=xiy*rnxi+ety*rnet
 1002 continue
       geome(7,ielem)=rj
c
c      end of loop over the elements
c
 1000 continue
      return
      end
c
c---------------------------------------------------------------------
c
c     this subroutine evaluate the first and second derivatives
c     using variational recovery 
c
c---------------------------------------------------------------------
c
       subroutine getvar(nnode,namat,ngeom,nderv,nelem,
     &                   ndimn,npoin,coord,intmat,derip,
     &                   deri2,geome,mmat,unkno,
     &                   ivari)
       real mmat
       dimension    unkno(namat,npoin),mmat(npoin)
       dimension    geome(ngeom,nelem),coord(ndimn,npoin)
       dimension    derip(nderv,npoin),deri2(4,npoin)
       dimension  intmat(nnode,nelem),node(3)
       dimension  nx(4),ny(4)
c
c     initialize
c
       call rfillv(derip,nderv*npoin,0.0)
       call rfillv(deri2,4*npoin,0.0)
c
c     obtain mmat (already inverted)
c
       call gtmmat(nelem,npoin,nnode,ngeom,intmat,
     &             geome,mmat)
       gamma=1.4
c
c      loop over the elements
c
       do 7000 ielem=1,nelem
c
       roelx=0.0
       roely=0.0
c
c     obtain the values needed
c
       do 7100 inode=1,nnode
       ipoin=intmat(inode,ielem)
       ro=unkno(1,ipoin)
       u1=unkno(2,ipoin)
       u2=unkno(3,ipoin)
       en=unkno(4,ipoin)
       vsq=u1*u1+u2*u2
       pre=(gamma-1.)*ro*(en-0.5*vsq)
       if(ivari.eq.1) ropoi=ro
       if(ivari.eq.2) ropoi=sqrt(vsq)
       if(ivari.eq.3) ropoi=pre/(ro**gamma)
       if(ivari.eq.4) ropoi=sqrt((vsq*ro)/(gamma*pre))
       if(ivari.eq.5) ropoi=pre
       rnx  =geome(inode      ,ielem)
       rny  =geome(inode+nnode,ielem)
       roelx=roelx+rnx*ropoi
       roely=roely+rny*ropoi
       node(inode)=ipoin
         nx(inode)=rnx
         ny(inode)=rny
 7100 continue
c
c     element area
c
      rj=geome(ngeom,ielem)*0.5/3.0
c
c     assembly
c
       do 1101 inode=1,nnode
       ipoin=node(inode)
       derip(1,ipoin)=derip(1,ipoin)+rj*roelx
       derip(2,ipoin)=derip(2,ipoin)+rj*roely
 1101  continue
c
c      end of loop over the elements
c
 7000 continue
c
c      multiply by mmat(inverted) and divide by the nr. of el./per node
c
       do 7001 ipoin=1,npoin
       do 7002 j=1,2
       derip(j,ipoin)=derip(j,ipoin)*mmat(ipoin)
 7002  continue
 7001  continue
c
c     obtain the second variations
c
       do 8000 ielem=1,nelem
       roexx=0.0
       roexy=0.0
       roeyx=0.0
       roeyy=0.0
c
c      obtain the values needed
c
       do 8100 inode=1,nnode
       ipoin=intmat(inode,ielem)
       ropox=derip(1,ipoin)
       ropoy=derip(2,ipoin)
       rnx  =geome(inode      ,ielem)
       rny  =geome(inode+nnode,ielem)
       roexx=roexx+rnx*ropox
       roexy=roexy+rny*ropox
       roeyx=roeyx+rnx*ropoy
       roeyy=roeyy+rny*ropoy
       node(inode)=ipoin
         nx(inode)=rnx
         ny(inode)=rny
 8100 continue
c
c     element area
c
      rj=geome(ngeom,ielem)*0.5/3.0
c
c     assembly
c
       do 8101 inode=1,nnode
       ipoin=node(inode)
       deri2(1,ipoin)=deri2(1,ipoin)+rj*roexx
       deri2(2,ipoin)=deri2(2,ipoin)+rj*roexy
       deri2(3,ipoin)=deri2(3,ipoin)+rj*roeyx
       deri2(4,ipoin)=deri2(4,ipoin)+rj*roeyy
 8101  continue
c
c     end of loop over the elements
c
 8000 continue
c
c     multiply by mmat(inverted) and divide by the nr. of el./per node
c
       do 8001 ipoin=1,npoin
       do 8002 j=1,4
       deri2(j,ipoin)=deri2(j,ipoin)*mmat(ipoin)
 8002  continue
 8001  continue
c
c     symmetrize !!!!!! (equate cross derivatives)
c
       do 8060 ipoin=1,npoin
       tauxy=(deri2(2,ipoin)+deri2(3,ipoin))/2.
       deri2(2,ipoin)=tauxy
       deri2(3,ipoin)=tauxy
 8060  continue
c
      return
      end
c
c---------------------------------------------------------------------
c   
c      evaluate the lamped mass matrix 
c
c----------------------------------------------------------------------
c
       subroutine gtmmat(nelem,npoin,nnode,ngeom,
     *                   intmat,geome,mmat)
       real mmat
       dimension geome(ngeom,*),mmat(*)
       dimension  intmat(nnode,*)
       c6=1./6.
       call rfillv(mmat,npoin,0.0)
       do 1000 ielem=1,nelem
c
c      jacobian of the element
c
         rj =geome(7,ielem)
         rj6=rj*c6
c
c      add to mmat
c
         do 1020 inode=1,nnode
           in=intmat(inode,ielem)
           mmat(in)=mmat(in)+rj6
1020     continue
1000   continue
       do 2000 i=1,npoin
         mmat(i)=1./mmat(i)
2000   continue
       return
       end
c
c---------------------------------------------------------------------
c
c     this subroutine fills a vector with integer values
c
c---------------------------------------------------------------------
c
      subroutine rfiliv(ia,n,ivalue)
      dimension  ia(n)
      do 1000 i=1,n
        ia(i)=ivalue
1000  continue
      return
      end
c
c---------------------------------------------------------------------
c
c     this subroutine fills a vector with real values
c
c---------------------------------------------------------------------
c
       subroutine rfillv(a,na,c)
       dimension a(na)
       do 1000 i=1,na
         a(i)=c
 1000  continue
       return
       end
c
c---------------------------------------------------------------------
c
c      function to convert from number to character
c
c----------------------------------------------------------------------
c
       character*5 function convrt (number)
       character*1 ch 
       dimension ch(5)
       do 5, n=1, 5
         ch(n)=' '
5      continue
       k=0
       i1=number
10     j= mod (i1, 10)
       k=k+1
       ch(k) = char(ichar('0')+j)
       i1 = i1/10
       if (i1.ne.0) goto 10
       convrt = ch(5)//ch(4)//ch(3)//ch(2)//ch(1)
       return
       end
c
c---------------------------------------------------------------------
c
c      plotting subroutine
c
c----------------------------------------------------------------------
c
       subroutine plot(iop,x,y,st,xmin,xmax,ymin,ymax)
c      include '/usr/include/f77/usercore77.h'
       integer  white, black
       common /plt/ iplot
       parameter (white=254, black=255,nil=-1 )
       character* 80 st,textar(10)
       dimension xarray(10),yarray(10)
c      integer vsurf(vwsurfsize)
       goto (10,20,30,40,50,60,70) iop
  10   if (iplot.eq.0) then
c      call getviewsurface(vsurf)
c      call initializecore (basic, noinput, twod)
c      call initializevwsurf(vsurf, true)
c      call selectvwsurf(vsurf)
c      call setndcspace2(1.,0.781)
c      call defcolorindices (vsurf, white, white, 1.0, 1.0, 1.0)
c      call defcolorindices (vsurf, black, black, 0.0, 0.0, 0.0)
c      call setviewport2 (0.,1.,0.0,.781)
c      call setwindow (0.,35.,0.,5.)
c      call createtempseg ()
c      call fillbox (black,0.,0.,35.,5.)
c      call closetempseg()
c      call setviewport2 (0.,1.,.66,.781)
c      call setwindow (0.,35.,0.,5.)
c      call createtempseg ()
c      call settextindex(white)
c      call setlineindex(white)
c      call fillbox (black,0.,0.,35.,5.)
c      call setcharprecision(character)
c      call setfont (roman)
c      call setcharsize (.5,.5)
c      call moveabs2(4.,2.5)
c      st='   swan2d   finite element mesh generator  '
c      call text(st)
c      call closetempseg()
       else
       st='   swan2d   finite element mesh generator  '
       print *,st
       end if
       return
  20   if(iplot.eq.0) then
c      call setviewport2 (0.,1.,0.,.11)
c      call setwindow (0.,35.,0.,5.)
c      call createtempseg ()
c      call settextindex(white)
c      call setlineindex(white)
c      call fillbox (black,0.,0.,35.,5.)
c      call setcharprecision(character)
c      call setfont (roman)
c      call setcharsize (.5,.5)
c      call moveabs2(4.,2.5)
c      call text(st)
c      call closetempseg()
       else
       print *,st
       end if
       return
  30   if(iplot.eq.0) then
c      call setviewport2 (0. ,1. , 0.11  , .66)
c      call setwindow (xmin,xmax,ymin,ymax)
c      call createtempseg ()
c      call settextindex(white)
c      call setlineindex(white)
c      call fillbox (black,xmin,ymin,xmax-xmin,ymax-ymin)
c      call closetempseg()
c      call setviewport2 (0.25,.75 , 0.135  , .635)
c      xa=xmax-xmin
c      ya=ymax-ymin
c      xb=(xmin+xmax)/2.
c      yb=(ymin+ymax)/2.
c      if(xa.gt.ya)then
c         xl=xb-0.5*xa
c         xr=xb+0.5*xa
c         yl=yb-0.5*xa
c         yu=yb+0.5*xa
c       else
c         xl=xb-0.5*ya
c         xr=xb+0.5*ya
c         yl=yb-0.5*ya
c         yu=yb+0.5*ya
c       end if
c      call setwindow (xl,xr,yl,yu)
c      call createtempseg ()
       end if
       return
  40   if(iplot.eq.0)then
c        call moveabs2(x,y)
       end if
       return
  50   if(iplot.eq.0)then
c        call lineabs2(x,y)
       end if
       return
  60   if(iplot.eq.0) then
c        call closetempseg()
       end if
       return
  70   if(iplot.eq.0) then
c      call deselectvwsurf(vsurf)
c      call terminatecore()
       end if
       return
       entry getopt(iop)
       if(iplot.eq.0) then
c      textar(1) = 'quit.'
c      textar(2) = 'plot the mesh .'
c      textar(3) = 'save data. '
c      textar(4) = 'check the mesh. '
c      textar(5) = 'mesh cosmetics.'
c      stepsize = 80./(2*5+1)
c      do 150 i=1,5
c        yarray(i) = (2*i-1)*stepsize
c        xarray(i) = 10.
c150   continue
c      xarray(5+1)=nil
c      call setviewport2 (0.74 ,0.99 , 0.3  , .5)
c      call setwindow (0.,80.,0.,80.)
c      call createtempseg ()
c      call settextindex(black)
c      call setlineindex(black)
c      call fillbox (white,0.,0.,80.,80.)
c      call setcharprecision(character)
c      call setfont (roman)
c      call setcharsize (2.,2.)
c      i=1
c190   if (xarray(i).ne.nil) then
c        call fillbox (white,xarray(i), yarray(i), 60., 
c    1                 stepsize)
c        call moveabs2 (xarray(i)+60./10.0, 
c    1                        yarray(i)+stepsize/3.0)
c         call text (textar(i))
c        i=i+1
c        goto 190
c      end if
c      call initializedevice (locator, 1)
c      call setechosurface (locator, 1, vsurf)
c      call setecho (locator, 1, 1)
c      call setechoposition (locator, 1, 0.5 ,0.5)
c      call initializedevice (button, 1)
c      call setechosurface (button, 1, vsurf)
c  55  call awtbuttongetloc2 (0.0,1,num, xd,yd)
c      if(num.eq.0) goto 55
c      xstep = (80.-0.)/(0.99-0.74)
c      ystep = (80.-0.)/(0.5-0.3)
c      viewx = (xd-0.74)*xstep+0.
c      viewy = (yd-0.3)*ystep+0.
c      n = int(viewy/stepsize)+1
c      if ((n.gt.0).and.mod(n,2).eq.0) then
c        index = n/2
c      else
c        goto 55
c      end if
c      call terminatedevice (locator, 1)
c      call terminatedevice (button, 1)
c      iop=5-index+1
c      call fillbox (black,0.,0.,80.,80.)
c      call closetempseg()
       else
  43   print *, '                      '
       print * ,'  1.  mesh cosmetics  '
       print * ,'  2.  check the mesh  '
       print * ,'  3.  save data file  '
       print * ,'  4.  plot the  mesh  '
       print * ,'  5.       quit       '
       print *, '                      '
       write(*,49) 
  49   format('    input  option     :',$)
       read(*,*,err=43) iop
       if(iop.lt.1.or.iop.gt.5) goto 43
       end if
       return
       end
c      subroutine fillbox (index, x, y, dx, dy)
c      integer fillindex
c      dimension x1(4), y1(4)
c      call inqfillindex (fillindex)
c      call setfillindex (index)
c      x1(1) = x
c      x1(2) = x
c      x1(3) = x+dx
c      x1(4) = x+dx
c      y1(1) = y
c      y1(2) = y+dy
c      y1(3) = y+dy
c      y1(4) = y  
c      call polygonabs2 (x1, y1, 4)
c      call setfillindex (fillindex)
c      call moveabs2 (x, y)
c      call linerel2 (0, dy)
c      call linerel2 (dx, 0)
c      call linerel2 (0, -dy)
c      call linerel2 (-dx, 0)       
c      return   
c      end    
c
c*    [spaln] computes the spacing at a point from a source line               *
c*-----------------------------------------------------------------------------*
      real function spaln(xsrc1,xsrc2,xr)
      real   xsrc1(*),xsrc2(*),xr(*)
      real   xsrc3(5)
c
      tolg = 1.e-05
c     ... calculates the projection
      xsrc3(1) = xsrc2(1)-xsrc1(1)
      xsrc3(2) = xsrc2(2)-xsrc1(2)
      al       = sqrt(xsrc3(1)**2+xsrc3(2)**2)
      if(al.lt.tolg) stop ' error in defining line source'
      al1      = 1./al
      xsrc3(1) = xsrc3(1)*al1
      xsrc3(2) = xsrc3(2)*al1
      sca = (xr(1)-xsrc1(1))*xsrc3(1)+
     -      (xr(2)-xsrc1(2))*xsrc3(2)
      if(sca.le.0.0) then
        spaln  = spapt(xsrc1,xr) 
      else if(sca.ge.al) then
        spaln  = spapt(xsrc2,xr) 
      else
        w2       = sca*al1
        w1       = 1.-w2
        xsrc3(1) = w1*xsrc1(1)+w2*xsrc2(1)
        xsrc3(2) = w1*xsrc1(2)+w2*xsrc2(2)
        xsrc3(3) = w1*xsrc1(3)+w2*xsrc2(3)
        xsrc3(4) = w1*xsrc1(4)+w2*xsrc2(4)
        xsrc3(5) = w1*xsrc1(5)+w2*xsrc2(5)
        spaln    = spapt(xsrc3,xr) 
      endif
c
      return
      end
c*-----------------------------------------------------------------------------*
c*    [spapt] computes the spacing at a point from a source point              *
c*-----------------------------------------------------------------------------*
      real function spapt(xsrc,xr)
      real xsrc(*),xr(*)
c
      BIG = 500.
      CLG2 = log(2.0)
c
      r   = sqrt((xsrc(1)-xr(1))**2+
     -           (xsrc(2)-xr(2))**2)
      if(r.le.xsrc(4)) then
        spapt = xsrc(3)
      else
        r     = r-xsrc(4)
        rr    = xsrc(5)-xsrc(4)
        if(rr.le.0.) stop 'incorrect source distances'
        rr    = CLG2/rr
        ae    = min(r*rr,BIG)
        spapt = xsrc(3)*exp(ae)
      endif
c
      return
      end
c
c -------------------------------------------------------------------
c
c      subroutine to generate boundary points
c
c -------------------------------------------------------------------
c
      subroutine genmsh(nbcs  ,ibcs  ,nbno  ,npoig ,neleg ,ieleg ,    
     *                  coorg ,intmeg,nboun ,node  ,coorn ,coor  ,    
     *                  delta ,nnn   ,lcoor ,lbou  ,
     *                  lboud ,strec ,mxseg ,mxpoi ,mxbou ,xmin  ,     
     *                  xmax  ,ymin  ,ymax  ,icond ,nn    ,ln    ,
     *                  coorl ,lpr   ,lboul ,iel   ,nelem ,nel   ,
     *                  npl   ,lelch ,deltmi,npfrt ,nqfrt ,
     *                  nregi ,helpa ,lep   ,lmb   ,lwher ,lhowm ,
     *                  icone ,mxele ,unkng , unkno ,io , lbk , abk,
     *                  lplay ,lelay ,nlay,lafter)
      parameter (  mxb = 50000 )
      dimension nbno(mxseg,*) ,ieleg(3,*) ,coorg(2,*) ,coorn(2,*) 
      dimension intmeg(3,*)   ,coor(2,*)  ,delta(4,*) ,strec(4,*)
      dimension nnn(*)        ,lcoor(*)   ,lbou(2,*)  ,ibcs(*)          
      dimension lboud(*)      ,ln(mxseg,*),nn(*)      ,helpa(*)    
      dimension icond(*)      ,coorl(*) ,lpr(*)     ,lboul(*)
      dimension iel(3,*)      ,lelch(*)   ,npfrt(*)   ,nqfrt(*)
      dimension nregi(*)      ,lep(*)     ,ar(3)      ,lmb(*)
      dimension lwher(*)      ,lhowm(*)   ,icone(*)  , mn(2)
      dimension unkng(4,*)    ,unkno(4,*) ,  d1(mxb) , d2(mxb)
      dimension lafter(*)     ,lstak( mxb)
      dimension kode(5), lcent(0:mxb),nback(mxb),nhead(mxb),lgen(mxb)
      dimension lbk(*), abk(*),xr(2),hlay(1000), lplay(*),lelay(*)
      common /plt/ iplot
      character *80 text,st
      character *5 convrt
      data pi/3.1415927/
      data kode/1,2,3,1,2/
c
      ilast= 1
      mcomp = 0
      nend = 0
      xmax=-1.0e+6
      xmin= 1.0e+6
      ymax=-1.0e+6
      ymin= 1.0e+6
c
c *** read basic data corresponding to the new mesh
c
      read(8,'(a)',err=24) text
      read(8,*,err=24) nfn,nbcs,nbvis,nlay,hmin
      lheight = 1
      if(nlay.lt.0) lheight = -1
      nlay = abs(nlay)
      if(nfn.gt.mxpoi) stop 'error 50'
      if(nbcs.gt.mxseg) stop 'error 60'
c
c *** specify coordinates of fixed nodes
c *** find out max and min coordinates for plotting
c
      do 10 in=1,nfn
        read(8,*,err=24) ip,(coorn(j,ip),j=1,2)
        xl=coorn(1,ip)
        yl=coorn(2,ip)
        xmin = min (xmin, xl)
        xmax = max (xmax, xl)
        ymin = min (ymin, yl)
        ymax = max (ymax, yl)
10    continue
      do 15 ibs=1,nbcs
        read(8,*,err=24) ib,nib,ico
        il=abs(ib)
        nn(il)=nib
        icond(il)=ico
        ibcs(il)=ib
        if(nn(il).gt.mxbou) stop 'error 70'
        read(8,*,err=24) (ln(il,i),i=1,nn(il))
  15  continue
c
c ***  for viscous adaption select variable for B.L
c
      if(nbvis.ne.0.and.io.ne.0) then
      print *,'input index for b.l:(1=r , 2=u , 3=p , 4=m , 5=t)'
      read(*,*) idn
      end if
c
c *** find minimum grid spacing
c
      deltmi=1.e+6
      do 11 ip=1,npoig
      deltmi=min(deltmi,delta(1,ip))
 11   continue
c
c *** specify region boundary segments
c
  18  node =0
      nbou =0
      ne   =0
      lnode=0
      call rfiliv(lboud,mxpoi,0)
      call rfiliv(lboul,mxpoi,0)
      call rfiliv(lafter,mxpoi,0)
      call rfiliv(lmb,mxpoi,0)
      call rfiliv(lgen,mxb,0)
      call rfiliv(lcent,mxb,0)
      call rfiliv(nback,mxb,0)
      call rfiliv(nhead,mxb,0)
      lcent(0) = 0
      if(nbvis.eq.0) goto 999
c
      ncoun = 0
      do 20 ibs=1,nbvis
        ncoun= node
        call interp(npoig ,neleg ,node  ,coorg ,ieleg ,intmeg,delta , 
     *              coor  ,coorn ,nbno  ,nnn   ,lcoor ,lbou  ,nbou  ,  
     *              ibs   ,ilast ,lboud ,mxbou ,mxpoi ,mxseg ,nn    ,
     *              ln    ,icond ,lelch ,deltmi,     0, lbk,abk)
        npsg = node - ncoun
        st=' segment'//convrt(ibs)//' has finished ,  Np = '
     *              //convrt(npsg)
        call   plot(2     ,x     ,y     ,st    ,xmin  ,xmax  ,ymin  ,
     *              ymax  )
   20 continue
      nboun=node
      do 31 ibs=1,nbvis
        do 34 i=1,2
          j=1
          if(i.eq.2) j=nnn(ibs)
          kpoi=nbno(ibs,j)
          do 32 ibou=1,nbou
            knew=lbou(2,ibou)
            ktry=lbou(1,ibou)
            if(ktry.eq.kpoi) goto 33
32        continue
33        nbno(ibs,j)=knew
34      continue
31    continue
c
c *** identify trailing edge
c
      do 766 ib = 1, nbcs
       ip1 = ln(ib,1)
       if(lmb(ip1).eq.100)goto 766
       do 767 ibb = 1,nbcs
        ip2 = ln(ibb,nn(ibb))
        if(ip1.eq.ip2.and.(ib.le.nbvis.or.ibb.le.nbvis))then
          ip3 = ln(ib,2)
          ip4 = ln(ibb,nn(ibb)-1)
          cosalp = (coorn(1,ip3)-coorn(1,ip1))*
     &             (coorn(1,ip4)-coorn(1,ip2))
     &           + (coorn(2,ip3)-coorn(2,ip1))*
     &             (coorn(2,ip4)-coorn(2,ip2))
c
          tlen1 = (coorn(1,ip3)-coorn(1,ip1))**2
     &          + (coorn(2,ip3)-coorn(2,ip1))**2
          tlen2 = (coorn(1,ip4)-coorn(1,ip2))**2 
     &          + (coorn(2,ip4)-coorn(2,ip2))**2
          cosalp = cosalp/sqrt(tlen1*tlen2)
          tnx1 = (coorn(2,ip3) - coorn(2,ip1))/sqrt(tlen1)
          tny1 = (coorn(1,ip1) - coorn(1,ip3))/sqrt(tlen1)
          tnx2 = (coorn(2,ip2) - coorn(2,ip4))/sqrt(tlen2)
          tny2 = (coorn(1,ip4) - coorn(1,ip2))/sqrt(tlen2)
          tlen = (tnx1+tnx2)**2 + (tny1+tny2)**2
          tn1  = (tnx1 + tnx2)/sqrt(tlen)
          tn2  = (tny1 + tny2)/sqrt(tlen)
          if(cosalp.gt.0.8.and.ibb.le.nbvis.and.ib.le.nbvis)then
            lmb(nbno(ib,1)) = 100
            lplay(nbno(ib,1)) = -1
            print *,nbno(ib,1),' is a trailing edge'
            goto 766
          endif
        endif
 767   continue
 766  continue
c
      spacing = 10000.
      in=0
      lnode=0
      do 520 iseg=1,nbvis
       noseg=abs(ibcs(iseg))
       np1=nbno(noseg,1)
       nk=nnn(noseg)
       np2=nbno(noseg,nk)
       nl1=nk-1
       nl0= 1
       nst=1
       na =0
       if(ibcs(iseg).lt.0) then
         nl0 = nl1
         nl1 = 1
         nst = -1
         na  =  1
       end if
       do 550 kn= nl0,nl1,nst
        lnode=lnode+1
        npfrt(lnode)=nbno(noseg,kn+na)
        nqfrt(lnode)=nbno(noseg,kn+1-na)
        i = nbno(noseg,kn+na)
        x=coor(1,i)
        y=coor(2,i)
        inorm=1
        call findel(npoig ,neleg ,coorg ,ieleg ,intmeg,x    ,y    ,
     *              ilast ,ar    ,i1    ,i2    ,i3    ,ienr  ,lelch)
        call getval(npoig ,i1    ,i2    ,i3    ,ar    ,delta,dis  ,
     *              alp   ,anx   ,any   ,inorm )
        xr(1) = x
        xr(2) = y
        if(lbk(1).ne.0) then
          do 117 ip = 1,lbk(1)
          j1   = lbk(3)+(ip-1)*5
          sr   = spapt(abk(j1),xr)
          dis = min(dis,sr)    
  117     continue
        endif    
c
        if(lbk(2).ne.0) then
          do 127 ip = 1,lbk(2)
          j1   = lbk(4)+(ip-1)*10
          j2   = j1+5
          sr   = spaln(abk(j1),abk(j2),xr)
          dis = min(dis,sr)    
  127     continue
        endif    
c
        strec(1,i)=dis
        strec(2,i)=alp
        strec(3,i)=anx
        strec(4,i)=any
        call getval(npoig ,i1    ,i2    ,i3    ,ar    ,unkng,dis  ,
     *             alp   ,anx   ,any   ,0 )
        unkno(1,i)=dis
        unkno(2,i)=alp
        unkno(3,i)=anx
        unkno(4,i)=any
        lep(i) = ienr
        if(lmb(i).eq.100) spacing = min(spacing,strec(1,i))
  550  continue
       i = nbno(noseg,kn+na)
       x=coorn(1,i)
       y=coorn(2,i)
       inorm=1
       call findel(npoig ,neleg ,coorg ,ieleg ,intmeg,x    ,y    ,
     *            ilast ,ar    ,i1    ,i2    ,i3    ,ienr  ,lelch)
       call getval(npoig ,i1    ,i2    ,i3    ,ar    ,delta,dis  ,
     *             alp   ,anx   ,any   ,inorm )
       xr(1) = x
       xr(2) = y
       if(lbk(1).ne.0) then
         do 217 ip = 1,lbk(1)
          j1   = lbk(3)+(ip-1)*5
          sr   = spapt(abk(j1),xr)
          dis = min(dis,sr)    
  217   continue
       endif    
c
       if(lbk(2).ne.0) then
         do 227 ip = 1,lbk(2)
          j1   = lbk(4)+(ip-1)*10
          j2   = j1+5
          sr   = spaln(abk(j1),abk(j2),xr)
          dis = min(dis,sr)    
  227    continue
       endif    
c
       strec(1,i)=dis
       strec(2,i)=alp
       strec(3,i)=anx
       strec(4,i)=any
       call getval(npoig ,i1    ,i2    ,i3    ,ar    ,unkng,dis  ,
     *             alp   ,anx   ,any   ,0 )
       unkno(1,i)=dis
       unkno(2,i)=alp
       unkno(3,i)=anx
       unkno(4,i)=any
       lep(i) = ienr
       if(lmb(i).eq.100) spacing = min(spacing,strec(1,i))
  520 continue
c
      do 311 i = 1 , lnode
       i1 = npfrt(i)
       i2 = nqfrt(i)
       do 312 j = 1 , lnode
        if(j.eq.i) goto 312
        if(nqfrt(j).eq.i1) then
          nback(i) = j
          nhead(j) = i
          lcent(i) = i1
          if(lmb(i1) .eq. 0) lmb(i1) = 1
          goto 311
        end if
  312  continue
       lcent(i)= i1
  311 continue
c
      do 313 i = 1 , lnode
       i0 = lcent(nback(i))
       i1 = lcent(i)
       i2 = lcent(nhead(i))
       if(lmb(i1).ne.0) then
         dx1=coor(1,i2)-coor(1,i1)
         dy1=coor(2,i2)-coor(2,i1)
         dx2=coor(1,i1)-coor(1,i0)
         dy2=coor(2,i1)-coor(2,i0)
         d1(i)= dy1/sqrt(dx1*dx1+dy1*dy1)+
     &            dy2/sqrt(dx2*dx2+dy2*dy2)
         d2(i)=-dx1/sqrt(dx1*dx1+dy1*dy1)-
     &            dx2/sqrt(dx2*dx2+dy2*dy2)
         al   = sqrt(d1(i)**2+d2(i)**2)
         d1(i)= d1(i)/al
         d2(i)= d2(i)/al
       end if
       if(nback(i).eq.0) lgen(i) = -1
       if(nhead(i).eq.0) lgen(i) = -1
  313 continue
c
      ne    = 0
      read(8,*) (hlay(il),il=1,nlay)
      if(lheight.lt.0) then
       do i = nlay,2,-1
        hlay(i)=hlay(i)-hlay(i-1)
       end do
      end if
      il    = 0
  17  il    = il + 1
      knode = lnode
      w     = float(il)/float(nlay)
      do 314 i = 1 , knode
       ip = lcent(i)
       if(lmb(ip).eq.100.or.lmb(lcent(nhead(i))).eq.100.or.
     &    lmb(lcent(nback(i))).eq.100) goto 491 
       dx1= d1(nback(i))
       dy1= d2(nback(i))
       dx2= d1(nhead(i))
       dy2= d2(nhead(i))
       dd = sqrt((dx1+dx2)**2+(dy1+dy2)**2)
       d1(i) = (1-w)*d1(i) + w*(dx1+dx2)/dd
       d2(i) = (1-w)*d2(i) + w*(dy1+dy2)/dd
       dd    = sqrt(d1(i)**2+d2(i)**2)
       d1(i) = d1(i)/dd
       d2(i) = d2(i)/dd
  491  if(lmb(ip) .ne. 100.or.lgen(i).lt.il-1) goto 314
       sp = strec(1,ip)
       xne= coor(1,ip) - d1(i) * spacing
       yne= coor(2,ip) - d2(i) * spacing
       node = node + 1
       coor(1,node) = xne
       coor(2,node) = yne
       do 319 j = 1 , lnode
        if(j.eq.i) goto 319
        mn(1) = lcent(nback(j))
        mn(2) = lcent(j)
        if(mn(1).eq.ip.or.mn(2).eq.ip) goto 319
        xx = coor(1,mn(2))
        yy = coor(2,mn(2))
        dd = sqrt((xx-xne)**2+(yy-yne)**2)
        if(dd.lt.hlay(il)) then
          node = node -1
          goto 314
        end if
        call int2d(mn,ip,node,coor,ifnd)
        if(ifnd.eq.1) then
          node = node -1
          goto 314
        end if
  319  continue
       lnode = lnode+ 1
       npfrt(lnode) = node
       nqfrt(lnode) = ip
       lcent(lnode) = node
       nback(lnode) = lnode+1
       nhead(lnode) = i
       lnode = lnode + 1
       npfrt(lnode) = ip
       nqfrt(lnode) = node
       nhead(nback(i)) = lnode
       nback(lnode) = nback(i)
       nback(i)     = lnode -1
       nhead(lnode) = lnode -1
       lcent(lnode) = ip
       lmb(ip)      = 1
       lmb(node)    = 100
       lplay(node)  = -1
       lgen(lnode-1)= il
       lgen(lnode)  = il
       d1(lnode-1)  = d1(i)
       d2(lnode-1)  = d2(i)
       i0           = lcent(nback(i))
       i1           = lcent(i)
       i2           = lcent(nhead(i))
       dx1=coor(1,i2)-coor(1,i1)
       dy1=coor(2,i2)-coor(2,i1)
       dx2=coor(1,i1)-coor(1,i0)
       dy2=coor(2,i1)-coor(2,i0)
       d1(i)= dy1/sqrt(dx1*dx1+dy1*dy1)+
     &          dy2/sqrt(dx2*dx2+dy2*dy2)
       d2(i)=-dx1/sqrt(dx1*dx1+dy1*dy1)-
     &          dx2/sqrt(dx2*dx2+dy2*dy2)
       al   = sqrt(d1(i)**2+d2(i)**2)
       d1(i)= d1(i)/al
       d2(i)= d2(i)/al
       i0           = lcent(nback(lnode))
       i1           = lcent(lnode)
       i2           = lcent(nhead(lnode))
       dx1=coor(1,i2)-coor(1,i1)
       dy1=coor(2,i2)-coor(2,i1)
       dx2=coor(1,i1)-coor(1,i0)
       dy2=coor(2,i1)-coor(2,i0)
       d1(lnode)= dy1/sqrt(dx1*dx1+dy1*dy1)+
     &          dy2/sqrt(dx2*dx2+dy2*dy2)
       d2(lnode)=-dx1/sqrt(dx1*dx1+dy1*dy1)-
     &          dx2/sqrt(dx2*dx2+dy2*dy2)
       al   = sqrt(d1(lnode)**2+d2(lnode)**2)
       d1(lnode)= d1(lnode)/al
       d2(lnode)= d2(lnode)/al
       x=coor(1,node)
       y=coor(2,node)
       inorm=1
       call findel(npoig ,neleg ,coorg ,ieleg ,intmeg,x    ,y    ,
     *             ilast ,ar    ,i1    ,i2    ,i3    ,ienr  ,lelch)
       call getval(npoig ,i1    ,i2    ,i3    ,ar    ,delta,dis  ,
     *             alp   ,anx   ,any   ,inorm )
       xr(1) = x
       xr(2) = y
       if(lbk(1).ne.0) then
         do 617 ip = 1,lbk(1)
          j1   = lbk(3)+(ip-1)*5
          sr   = spapt(abk(j1),xr)
          dis = min(dis,sr)
  617    continue
       endif
c
       if(lbk(2).ne.0) then
         do 627 ip = 1,lbk(2)
          j1   = lbk(4)+(ip-1)*10
          j2   = j1+5
          sr   = spaln(abk(j1),abk(j2),xr)
          dis = min(dis,sr)
  627    continue
       endif
c
       strec(1,node)=dis
       strec(2,node)=alp
       strec(3,node)=anx
       strec(4,node)=any
       call getval(npoig ,i1    ,i2    ,i3    ,ar    ,unkng,dis  ,
     *             alp   ,anx   ,any   ,0 )
       unkno(1,node)=dis
       unkno(2,node)=alp
       unkno(3,node)=anx
       unkno(4,node)=any
       lep(node) = ienr
  314 continue
c
      ist = 0
      do i = 1 , lnode
       if(lplay(lcent(i)).ne.-1) lplay(lcent(i)) = il
c      if(lmb(lcent(i)).eq.100) then
c      ist = ist + 1
c      lstak(ist) = nback(i)
c      ist = ist + 1
c      lstak(ist) = nhead(i)
c      write(99,*) i,nback(i),nhead(i),ist
c      end if
      end do
c     jst = 1
c
c     do kk = 1 , 2
      do 315 i = 1 , lnode
c213   i = lstak(jst)
c      jst = jst + 1
       i0 = lcent(nback(i))
       i1 = lcent(i)
       i2 = lcent(nhead(i))
c      do kst = 1, ist
c       if(lstak(kst).eq.nback(i)) goto 214
c      end do
c      ist = ist + 1
c      lstak(ist) = nback(i)
c214   do kst = 1, ist
c       if(lstak(kst).eq.nhead(i)) goto 215
c      end do
c      ist = ist + 1
c      lstak(ist) = nhead(i)
c215   continue
c      write(99,*) i,nback(i),nhead(i),ist
c      if(((lmb(i0).eq.100.or.lmb(i2).eq.100).and.kk.eq.1).or.
c    &    ((lmb(i0).ne.100.and.lmb(i2).ne.100).and.kk.eq.2)) then
       if(lgen(i).lt.il-1) then
c	   write(99,*) ' lgen'
c	   write(99,*) il,i,lgen(i),nback(i),nhead(i),
c    &               i1,i0,i2,coor(1,i1),coor(2,i1)
c	   write(99,*) '-----------------------------------'
	   goto 315
	 end if
       if(lmb(i1).eq.100) then
c	   write(99,*) ' lmb'
c	   write(99,*) il,i,lmb(i1),i1,i0,i2,coor(1,i1),coor(2,i1)
c	   write(99,*) '-----------------------------------'
	   goto 315
	 end if
       if(lgen(nback(i)).lt.il-1.or.lgen(nhead(i)).lt.il-1) then
c	   write(99,*) ' lgen b h'
c	   write(99,*) il,i,nback(i),nhead(i),lgen(nback(i)),
c    &               lgen(nhead(i)),i1,i0,i2,
c    &               coor(1,i1),coor(2,i1)
c	   write(99,*) '-----------------------------------'
	   goto 315
	 end if
       if(hlay(il).ge.strec(1,i1)) then
c	   write(99,*) ' strec'
c	   write(99,*) il,i,i1,coor(1,i1),coor(2,i1),strec(1,i1)
c	   write(99,*) il,i,i0,coor(1,i0),coor(2,i0),strec(1,i0)
c	   write(99,*) il,i,i2,coor(1,i2),coor(2,i2),strec(1,i2)
c	   write(99,*) '-----------------------------------'
	   goto 315
	 end if
       xn = coor(1,i1) - d1(i) * hlay(il)
       yn = coor(2,i1) - d2(i) * hlay(il)
c      dd0= sqrt((coor(1,i0)-xn)**2+(coor(2,i0)-yn)**2)
c      dd2= sqrt((coor(1,i2)-xn)**2+(coor(2,i2)-yn)**2)
c      if(dd0.lt.hlay(il).or.dd2.lt.hlay(il)) goto 315
       node = node +1
       coor(1,node) = xn
       coor(2,node) = yn
c
c  check intersection
c
       do 316 j = 1 , lnode
        if(j.eq.i) goto 316
        mn(1) = lcent(nback(j))
        mn(2) = lcent(j)
        if(i0.eq.0) goto 327
        if(mn(2).ne.i1.and.mn(2).ne.i0.and.mn(2).ne.i2) then
          xx = coor(1,mn(2))
          yy = coor(2,mn(2))
          dd = sqrt((xx-xn)**2+(yy-yn)**2)
          if(dd.lt.hlay(il)) then
            node = node -1
c	      write(99,*),' spacing'
c	      write(99,*) il,i,nback(i),nhead(i),i1,i0,i2
c	      write(99,*) '-----------------------------------'
            goto 315
          end if
        end if
        if(mn(1).eq.i0.or.mn(2).eq.i0) goto 316
        call int2d(mn,i0,node,coor,ifnd)
        if(ifnd.eq.1) then
          node = node -1
c	      write(99,*),' intersection'
c	      write(99,*) il,i,nback(i),nhead(i),i0,node,mn(1),mn(2)
c	      write(99,*) i0,coor(1,i0),coor(2,i0)
c	      write(99,*) node,coor(1,node),coor(2,node)
c	      write(99,*) mn(1),coor(1,mn(1)),coor(2,mn(1))
c	      write(99,*) mn(2),coor(1,mn(2)),coor(2,mn(2))
c	      write(99,*) '-----------------------------------'
          goto 315
        end if
  327   if(mn(1).eq.i2.or.mn(2).eq.i2) goto 316
        call int2d(mn,node,i2,coor,ifnd)
        if(ifnd.eq.1) then
          node = node -1
c	      write(99,*),' intersection'
c	      write(99,*) il,i,nback(i),nhead(i),node,i2,mn(1),mn(2)
c	      write(99,*) node,coor(1,node),coor(2,node)
c	      write(99,*) i2,coor(1,i2),coor(2,i2)
c	      write(99,*) mn(1),coor(1,mn(1)),coor(2,mn(1))
c	      write(99,*) mn(2),coor(1,mn(2)),coor(2,mn(2))
c	      write(99,*) '-----------------------------------'
          goto 315
        end if
  316  continue
c
       x=coor(1,node)
       y=coor(2,node)
       inorm=1
       call findel(npoig ,neleg ,coorg ,ieleg ,intmeg,x    ,y    ,
     *             ilast ,ar    ,j1    ,j2    ,j3    ,ienr  ,lelch)
       call getval(npoig ,j1    ,j2    ,j3    ,ar    ,delta,dis  ,
     *             alp   ,anx   ,any   ,inorm )
       xr(1) = x
       xr(2) = y
       if(lbk(1).ne.0) then
         do 417 ip = 1,lbk(1)
          k1   = lbk(3)+(ip-1)*5
          sr   = spapt(abk(k1),xr)
          dis = min(dis,sr)
  417    continue
       endif
c
       if(lbk(2).ne.0) then
         do 427 ip = 1,lbk(2)
          k1   = lbk(4)+(ip-1)*10
          k2   = k1+5
          sr   = spaln(abk(k1),abk(k2),xr)
          dis = min(dis,sr)
  427    continue
       endif
c
       strec(1,node)=dis
       strec(2,node)=alp
       strec(3,node)=anx
       strec(4,node)=any
       call getval(npoig ,j1    ,j2    ,j3    ,ar    ,unkng,dis  ,
     *             alp   ,anx   ,any   ,0 )
       unkno(1,node)=dis
       unkno(2,node)=alp
       unkno(3,node)=anx
       unkno(4,node)=any
       lep(node) = ienr
c
       if(i0.ne.0) then
         ne = ne  + 1
         iel(1,ne) = i0
         iel(2,ne) = i1
         iel(3,ne) = node
         lelay(ne) = nback(i)
       end if
       if(i2.ne.0) then
         ne = ne  + 1
         iel(1,ne) = i1
         iel(2,ne) = i2
         iel(3,ne) = node
         lelay(ne) = i
       else
c
c  *** check for open end
c
       end if
       lgen( i)  = il
       lcent(i)  = node
       lafter(i1) = node
c
c     end if
  315 continue
c     end do
c     if(jst.le.ist) goto 213
      print*,' layer=',il,' nelem=',ne
c
      if(il.lt.nlay) goto 17
      do i = 1 , node
       if(lmb(i).eq.100) lafter(i) = 0
      end do
c
c       set up initial set of nodes for region
c
 999  continue
      npl = node
      nel = ne
      nelem = ne
      do 324 i = 1 , lnode
       lplay(lcent(i)) = il + 1
       npfrt(i) = lcent(i)
       nqfrt(i) = lcent(nhead(i))
 324  continue
      print *,' do u want to return ??? (0= yes)  '
      read(*,*) irt
      if(irt.eq.0) return
       print *,npl,nel,nelem,node
      ncoun = node
      do 900 ibs=nbvis+1,nbcs
        call interp(npoig ,neleg ,node  ,coorg ,ieleg ,intmeg,delta , 
     *              coor  ,coorn ,nbno  ,nnn   ,lcoor ,lbou  ,nbou  ,  
     *              ibs   ,ilast ,lboud ,mxbou ,mxpoi ,mxseg ,nn    ,
     *              ln    ,icond ,lelch ,deltmi,     0,lbk,abk)
        npsg = node - ncoun
        ncoun= node
        st=' segment'//convrt(ibs)//' has finished ,  Np = '
     *              //convrt(npsg)
        call   plot(2     ,x     ,y     ,st    ,xmin  ,xmax  ,ymin  ,
     *              ymax  )
 900  continue
      do 61 ibs=nbvis+1,nbcs
        do 64 i=1,2
          j=1
          if(i.eq.2) j=nnn(ibs)
          kpoi=nbno(ibs,j)
          do 62 ibou=1,nbou
            knew=lbou(2,ibou)
            ktry=lbou(1,ibou)
            if(ktry.eq.kpoi) goto 63
  62      continue
  63      nbno(ibs,j)=knew
  64    continue
  61  continue
c
      do 920 iseg=nbvis+1,nbcs
      noseg=abs(ibcs(iseg))
      np1=nbno(noseg,1)
      nk=nnn(noseg)
      np2=nbno(noseg,nk)
      nl1=nk-1
      if(ibcs(iseg).lt.0) go to 940
      do 950 kn=1,nl1
      lnode=lnode+1
      npfrt(lnode)=nbno(noseg,kn)
      nqfrt(lnode)=nbno(noseg,kn+1)
  950 continue
      go to 920
  940 do 970 kn=1,nl1
      in=nl1+1-kn
      lnode=lnode+1
      npfrt(lnode)=nbno(noseg,in+1)
      nqfrt(lnode)=nbno(noseg,in)
  970 continue
  920 continue
c
      do 307 knode=1,lnode
        nregi(knode)=npfrt(knode)
307   continue
c
      nofrt=lnode
      nonr =lnode
c
c *** set up strec for the already existing nodes
c
      do 1070 j=1,nonr
          i=nregi(j)
          x=coor(1,i)
          y=coor(2,i)
          inorm=1
          call findel(npoig ,neleg ,coorg ,ieleg ,intmeg,x    ,y    ,
     *                ilast ,ar    ,i1    ,i2    ,i3    ,ienr  ,lelch)
          call getval(npoig ,i1    ,i2    ,i3    ,ar    ,delta,dis  ,
     *                alp   ,anx   ,any   ,inorm )
       xr(1) = x
       xr(2) = y
      if(lbk(1).ne.0) then
        do 717 ip = 1,lbk(1)
        i1   = lbk(3)+(ip-1)*5
        sr   = spapt(abk(i1),xr)
        dis = min(dis,sr)    
  717   continue
      endif    
c
      if(lbk(2).ne.0) then
        do 727 ip = 1,lbk(2)
        i1   = lbk(4)+(ip-1)*10
        i2   = i1+5
        sr   = spaln(abk(i1),abk(i2),xr)
        dis = min(dis,sr)    
  727   continue
      endif    
c
          strec(1,i)=dis
          strec(2,i)=alp
          strec(3,i)=anx
          strec(4,i)=any
          call getval(npoig ,i1    ,i2    ,i3    ,ar    ,unkng,dis  ,
     *                alp   ,anx   ,any   ,0 )
          unkno(1,i)=dis
          unkno(2,i)=alp
          unkno(3,i)=anx
          unkno(4,i)=any
          lep(i) = ienr
1070  continue
c
      st='              generating    please wait  '
      call   plot(2     ,x     ,y     ,st    ,xmin  ,xmax  ,ymin  ,
     *            ymax  )
       call triang(nofrt ,npfrt ,nqfrt ,nregi ,nonr  ,iel   ,coor  ,
     &             nelem ,node  ,neleg ,npoig ,ieleg ,intmeg,coorg ,
     &             delta ,strec ,toler ,ilast ,helpa ,lcoor ,lboul ,
     &             coorl(1)   ,coorl(mxpoi), lpr ,mxele ,mxpoi ,
     &             lelch ,lep  ,  0 ,lmb,unkng,unkno,lwher,lbk,abk)
c
      call   plot(3     ,x     ,y     ,st    ,xmin  ,xmax  ,ymin  ,
     *            ymax  )
      call gmplot(node,coor,ne,iel,iq,ne)
      return
  24  stop ' error 20 '
      end
c --------------------------------------------------------------------
c
      subroutine int2d(mn,m1,m2,coord,ifnd)
c
      dimension  coord(2,*) ,  mn(2)
c
      err                   = 1.0e-13
      ifnd                  = 0
      x1                    = coord(1,m1)
      y1                    = coord(2,m1)
      x2                    = coord(1,m2)
      y2                    = coord(2,m2)
      x3                    = coord(1,mn(1))
      y3                    = coord(2,mn(1))
      x4                    = coord(1,mn(2))
      y4                    = coord(2,mn(2))
      if(min(x1,x2).lt.max(x3,x4).or.max(x1,x2).gt.min(x3,x4).or.
     &   min(y1,y2).lt.max(y3,y4).or.max(y1,y2).gt.min(y3,y4)) return
      if(abs(x2-x1).le.err) then
        if(abs(x3-x4).le.err) return
        a2                  = (y4-y3)/(x4-x3)
        b2                  = y3 - a2 * x3
        xc                  = x1
        yc                  = a2 * xc + b2
      else if(abs(x3-x4).le.err) then
        a1                  = (y2-y1)/(x2-x1)
        b1                  = y1 - a1 * x1
        xc                  = x3
        yc                  = a1 * xc + b1
      else
        a2                  = (y4-y3)/(x4-x3)
        b2                  = y3 - a2 * x3
        a1                  = (y2-y1)/(x2-x1)
        b1                  = y1 - a1 * x1
        if(abs(a1-a2).le.err) return
        xc                  = (b2-b1)/(a1-a2)
        yc                  = a1 * xc + b1
      end if
      if(abs(yc-y3).lt.err.and.abs(xc-x3).lt.err) then
         print*,1,m1,m2,mn
         print *,x1,y1,x2,y2
         print *,x3,y3,x4,y4
         print *,xc,yc,a1,a2,b1,b2
        ifnd                = 1
      else if(abs(yc-y4).lt.err.and.abs(xc-x4).lt.err) then
         print*,2,m1,m2,mn
         print *,x1,y1,x2,y2
         print *,x3,y3,x4,y4
         print *,xc,yc,a1,a2,b1,b2
        ifnd                = 1
      else if(yc.le.max(y3,y4).and.yc.ge.min(y3,y4).and.
     &        xc.le.max(x3,x4).and.xc.ge.min(x3,x4).and.
     &        yc.le.max(y1,y2).and.yc.ge.min(y1,y2).and.
     &        xc.le.max(x1,x2).and.xc.ge.min(x1,x2)) then
         print*,3,m1,m2,mn
         print *,x1,y1,x2,y2
         print *,x3,y3,x4,y4
         print *,xc,yc,a1,a2,b1,b2
        ifnd                = 1
      end if
      return
      end


