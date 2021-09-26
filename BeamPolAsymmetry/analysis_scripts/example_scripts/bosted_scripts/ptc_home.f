      program ptc
c TODO 

c need HMS Cer cut on sp18 runs? save sh dist. 
c just like we do for cointime

c fine tune the three sets of cuts

c apply effic. cuts to iacc=1,2,3

c analyze al/d ratios

c add few more bins in pt

c apply cointime spectra for pi to get eff., contamination

c get the sig0's working

c Study ep elastic

c 

c check model for pi- Delta++

c try to add RF timing cut for kaons

c study beam cosphi, cos2phi terms

c add fits to ptb.f

c include R in sidis model

c pion absorbtion corrections

c det. efficiency corections

c try Harut binning in (x,q2)
c The range in 
c$x$ with $x_{min}-x_{max}$ (0.005-0.99), and
crange in
c $Q^2$ with $Q^2_{min}-Q^2_{max}$(1-2500) were
cdivided to
c $N_x$=60, and $N_q$=100 
cbins in log scale, respectively, with 
c$x_i=x_{min}e^{i\Delta x}$ and 
c$Q^2_i=Q^2_{min}e^{i\Delta Q^2}$, 

c Checks done
c fix simc excl. etc. on ptcrat rat2 plots
c HMS mispointing help phie distribution compared SIMC/
c pion MM at 0.94 for spring18 ok
c omega MM peak is around 0.77 (should be 0.78)
c look at omega peak at high pdelta. Look goodc
c runs 5600-17 had P=2.58 but should have been 2.53
c see if cer time cut can reduce proron contamination. Answer: no
c do kaon-lt runs and spring19 csv runs
c why big L / R assymetry in randoms in sp18?? 
c due to trig. inefficiency
c check of detector calibrations, adc windows, etc.
c seems to be big dependence on current spring18
c  see 3433-37, 3469-72. Fixed with large corr2.
c problems with runs 3421 and 23 compared to 25...
c fix problems with current_cut/current, dt runs 4914-37
c careful treatment of kaon decays in simc.
c add pi+pi0 channel
c fix problems with excl. model at lower W
c fix huge fluctuations in excl, rho, dummy 3696-99 for example
c rerun simc with 4x more sidis
c fix neg simc xsections.
c changed to DDS frag. fun. in SIMC for pi, K. Much better for K
c see if HKDS works better than DSS frag. fun.
c yes, changed to HKNS for pi, kept DDS for K
c fixed reason pion mm peak was a bit low?
c fixed slope in phicm at pi by adding HMS phi offset
c don;t add identical runs together first for cnt files
c KLT why neg. ppi<-12 not working? Now it works.
c are runs 3700 etc 2% higher than runs 3400 (two kin.) No.




      implicit none
      integer i,j,k,kk,kkk,irun,it,itlast,itrig,icc
      integer imps,ineg,ipos,icoin,itprev,sumi(20)
      integer accep(20,30,70,30),loedge(20),hiedge(20)
      integer lo_simc_hms(30,30),hi_simc_hms(30,30),irunlast
      real accepc(20,30,70,30),avy1,avy2,rat
      real accepz(20,30,70,30)
      integer okacc(2,30,70,30),totdiff(20),im,aerotray
      character*80 fname
      character*100 sstring
      real*8 weight,normfac,sum1,sum2,corr,tcut,e0_e,e0_p
      real*8 dpe,dpp,dppcorr,dthe,dthp,dphie,dphip,pe,the,pp
      real*8 thp,chrg,tefe,tefp,dt,ctpinew,thp_rad,th,bcm1
      real*8 dpp_init, xptar_init, yptar_init,yt_orig,yrast
      real*8 zk,phicmk,cthcmk,ptk,sumacc(10,2),r1,r2,r3,r4,r5
      real*8 sum1m(20),sum2m(20),sratem(20),mmpi2_cent_h
      real*8 mmpi2_cent_s,rex,rexer,tcutk,decdist,m2final
      real*8 sratemer(20),sratemx(20,6),sratemxer(20,6)
      real*8 sum1z(20),sum2z(20),xmaxaero,ymaxaero,aeromin
      real*8 sratez(20),sratezer(20),fstate(8,22,4)
      real*8 sratezx(20,6),sratezxer(20,6)
      real*8 m_mu/.105/,m_pi/0.139/,m_k/0.495/
      integer ipart,idist
      real*8 avnfac(8),sum3tot(8)
      real*8 ctpi,ctk,ctp,cere,cerng,ceraero,cerhg
      real*8 pre,prp,she,shp,avrter(20),eprot,aerot,aeront,evtype
      real*8 rate1,pipeak,kpeak,sum3,pav(2),thav(2),fcoin,avrt(20)
      real*8 savrt,savrter,timecorr,corr2,pgdsceff,hgdsceff,corr4
      real*8 toff(3000:7000),offset,kin(7000,5),savrtnox,savrtnoxer
      real*8 rate(7000),rateer(7000),fry,totchi2(10),totdf(10)
      real*8 ratem(20),ratemer(20),pirm(20),piaccm(20),corr3
      real*8 ratez(20),ratezer(20),pirz(20),piaccz(20),bcmcorr
      real*8 avrate(8,6,15), avrateer(8,6,15),rt(95,20)
      real*8 rter(95,20),chi,df,w,w2
      real*8 savrate(8,6,15), savrateer(8,6,15)
      real*8 savratenox(8,6,15), savratenoxer(8,6,15)
      real*8 srt(95),srter(95),aeroth(40,16)
      real*8 avkin(2,16,15,20,8)
      real*8 xexit,yexit,crad,voffset,hwid
      integer pexit,hexit,mmpi2hs(20),icase
      real*8 avratexz(4,4,4,6),avratexzer(4,4,4,6),gdsc
      real*8 srate(7000),srateer(7000),curmin,curmax
      real*8 sratex(7000,6),sratexer(7000,6),normfacx
      real*8 current,hbeta,pbeta,hztar,pztar,hxfp,hyfp,hdpfp,hypfp,sum4
      real*8 pxfp,pyfp,pxpfp,pypfp,ppprev,thpprev,sigsv(10),ptkeff
      real*8 ep0prev,the0prev,bcmhms(8000),piaccH,piaccL,x,q2,the0r
      integer ctpih(40),ctkh(40),cteh(40),itsv(7000),ix,iy,ctkhf(40)
      integer ctpihfa(40),ctphf(40),ctpihw(50,10),yth(19,16),ii,jj
      real*8 pir, piacc, ratec(7000), ratecer(7000), thplast,ep0,the0
      real*8 fact,rateck(6,20),ratecker(6,20),rateak(6,20)
      real*8 ratesk(6,20),ratesker(6,20),xaero,yaero
      real*8 pirk(6,20),piacck(6,20),pirks(6,20),piaccks(6,20)
      integer dpph(100),dpeh(100),dthph(100),dtheh(100),doit,nrtc
      integer dpphistrun(15),ythistrun(15)
      integer dphieh(100),dphiph(100),ctpihf(40),nrt,ikin,itp,ith
      integer hgh(30,50),ctrfrun(21,2),ctrfrn(10,20,2),xyaeroh(4,21)
      integer hghp(30,50),ngh(16)
c counts in bins of pt,phi,z,xq2,pol
      integer ixq2, ipt,iphi,iz,ihel,irr,iq2bin,ixbin
      integer cnts(2,16,15,20,8),iselec
      integer cntsct(20,40,4)
      integer cntsmmpi(50,8),cntsz(50,8),izz
      integer cntse(2,16,15,2),cntseh(2,16,15,8,2)
      integer cntsemch(2,16,15,8)
      real*8 cntsemc(2,16,15,3)
      real cntsmc(2,16,15,20,16),cntsmcmmpi(50,16),cntsmcz(50,16)
      integer ispireal, ispiacc, iskreal, iskacc, ispiaccL,ispiaccH
      real*8 sumcnt(3,2),sumcntsf(3,2),avsinphi,avsinphier
      real*8 rtc(99),rtcer(99),avrtc,avrtcer,ztnew,zdnew
      real*8 pplast,Empi,Emk,Emp,phi,mmpi2,mmk2,mmp2,hmean,pmean
      real*8 amp/0.9383/, ampi/0.1396/, amk/0.494/
      real*8 e0,ep,ppi,p_x(2),p_y(2),p_z(2),thpi,sigcc,sigcm
      real*8 v1,v2,v3,elreal
      integer mmpih(20,12),phih(20),meeh(20),mmph(20,14),eprh(14,21)
      integer meeht(20,2),mmpht(20),mmphs(20),mmpihs(20),mmpihsa(20)
      integer mmphf(400),mmphfs(400),mmpheta(400),ipm,betah(2,20)
      integer cereh(12,16),cerepih(12,16),ps1,ps2,ps3,ps4,ps5,p6
      integer shh(16),sheh(16),sh2h(16),prh(16,2),ctev(60,6,2)
      integer aeroh(16),aerohp(16),hnrf,pnrf,i1,i2,i3,i4,i5
      integer cthist(20,4,40),ctrfhist(20,4,40,2),ip,imax
      integer goodhodoh,goodhodop,mmpihfs(50)
      real*8 heffh,heffp,fptimeh,fptimep,okhms,okshms
      real*8 epv(4),p1vp(4),phicm,cthcm,pt,pt2,nu,zpi,mee,zh,zp,zs,zd
      real*8 hgtry,pgtry,prf,hrf,pfpt,hfpt,pfptc,pfptcr,pfptcc(10)
      real*8 hdcx,hdcxp,hdcy,hdcyp,pdcx,pdcxp,pdcy,pdcyp
      real*8 hfptc,hfptcr,hfptcc(20),dpshmslo,dpshmshi
      real*8 dphmslo,dphmshi,dthhms,dthshms,dphihms,dphishms
      real*8 e0last,ep0last,the0last,dtav,chrgtot
      real*8 tefeav,tefpav,corrav,corr2av,corr3av
      character*900 string
      character*80 title
      integer irate_e,irate_p,ielclean,ielcleanp
c 2 hodo planes, 100 x or y bins 3 y or x bins
      integer ctx(4,100,3,10),ctrf(5000,2)
      integer cty(4,3,100,10)
      real*8 xpad,ypad,rff,rfff,tmp
      integer npt,testing,usenewacc,ilo(4),ihi(4)
      real*8 ycoeff(10)
      real*8 del(100000),yfp(100000),ypfp(100000),
     >  ytarg(100000)
      common/rastdatacmn/del,yfp,ypfp,ytarg,npt
      integer simcevents,doaccstudy,accv(70)
      real*8 pi,phi0e,phi0p,xoff,xpoff
      real*8 dpz, dpn, dthz, dthn, dphiz, dphin, yn, yz
c target endcaps/walls in mm
      integer pxdist(100),pydist(100,7),sumk(20)

      real endcap_entr(2)/0.150, 0.130/
      real endcap_exit(2)/0.191, 0.188/
      real walls(2)/0.219, 0.184/
c For LH2, that's between 71.5 - 72.87 kg/m^3 (~1.2% variation), and for 
c LD2 it's between 165.1 - 168.7 kg/m^3 (~2.2% variation).
      real targden(2)/0.0723, 0.167/
      real tlen(2)/ 10.00, 10.00/
c density of Al used
      real alden/2.81/
c dummy foils in gm/cm2
      real dummy(2)/0.1816, 0.1815/
c so dummy subtraction factor is 
c h: .034 * 2.8 / 0.363 = .262
c d  .032 * 2.8 / 0.363 = .260

c if want to make code run fast, skkip simc events
      simcevents=1000000000

c cuts. Mark recommends -12 to 22 for SHMS (that is where matrix
c elements were fitted).
      dphmslo = -10.
      dphmshi =  10.
      dpshmslo = -12.
      dpshmshi =  22.
      dthhms = 0.070
      dthshms = 0.070
      dphihms = 0.025
      dphishms = 0.030
c a bit larger limits
      dphmslo = -12.
      dphmshi =  12.
      dpshmslo = -16.
      dpshmshi =  24.
      dthhms = 0.070
      dthshms = 0.070
      dphihms = 0.030
      dphishms = 0.030

! define spectrometer angles to match SIMC
      pi=3.141592653589793
      phi0e = 1.5 * pi
      phi0p = 0.5 * pi

      testing=0
      if(testing.ne.1) then
       open(unit=7,file='ptc.chist')
       open(unit=8,file='ptc.chistk')
       open(unit=20,file='ptc.ctpk')
       open(unit=14,file='ptc.phi')
       open(unit=15,file='ptc.mmpi')
       open(unit=16,file='ptc.out')
       open(unit=17,file='ptc.chistf')
       open(unit=18,file='ptc.rate')
       open(unit=19,file='ptc.chiste')
       open(unit=31,file='ptc.mmp')
       open(unit=32,file='ptc.mee')
       open(unit=33,file='ptc.chistkf')
       open(unit=36,file='ptcctb.top')
       open(unit=38,file='ptc.ctxy')
       open(unit=43,file='ptc.rf')
       open(unit=44,file='ptc.rff')
       open(unit=45,file='ptc.ztar')
       open(unit=46,file='ptc.beta')
       open(unit=47,file='ptc.aero')
       open(unit=48,file='ptc.sh')
       open(unit=52,file='ptc.she')
       open(unit=53,file='ptc.sh2')
       open(unit=59,file='ptc.pr')
       open(unit=54,file='ptc.ng')
       open(unit=49,file='ptcjpsi.txt')
       open(unit=51,file='ptc.cere')
       open(unit=55,file='ptc.mrate')
       open(unit=56,file='ptc.zrate')
       open(unit=81,file='ptcdpphist.txt')
       open(unit=82,file='ptcythist.txt')
       open(unit=85,file='ptc.eprh')
       open(unit=91,file='ptc.cmpsimc')
       open(unit=181,file='summedrunlist.txt')
      endif

      curmin=1000
      curmax=0.


       do kk=1,40
        do icc=1,16
         aeroth(kk,icc)=0.
        enddo
       enddo

c read in read in new acceptance limits from 4/20/20 based
c on fall18 runs (output ptc.accep)
       usenewacc = 1
       open(unit=5,file='ptc.goodaccep')
       do ipt=1,30
         do iphi=1,30
          read(5,'(10i4)') i1,i2,ilo,ihi
          do ith=1,70
           do irr=1,2
            i3 = okacc(irr,ipt,ith,iphi)
            okacc(irr,ipt,ith,iphi)=0
            if(ith.gt.ilo(2*irr-1) .and.
     >         ith.gt.ilo(2*irr) .and.
     >         ith.lt.ihi(2*irr-1).and.
     >         ith.lt.ihi(2*irr)) then
              okacc(irr,ipt,ith,iphi)=1
            endif
c            write(6,'(2i2,20i3)') i3,okacc(irr,ipt,ith,iphi),
c     >       ipt,i1,iphi,i2,irr,ith,ilo,ihi
           enddo
          enddo
         enddo
       enddo


c Time offsets
c no longer used
c      open(unit=5,file='ptc.toff')
c      do irun=3001,7000
c       toff(irun)= -1.0 
c      enddo
c      do i=1,10000
c       read(5,*,end=2) irun,offset
c       toff(irun)=offset
c      enddo
c 2    close(unit=5)

c BCM1 with current corrections from hms.f
c don't use this anymore
c      open(unit=5,file='ptc.hmsbcm')
c      do i=1,10000
c       read(5,*,end=3) irun,tmp
c       bcmhms(irun)=tmp
c      enddo
c 3    close(unit=5)

      nrt = 0
      nrtc = 0
      do kk=1,20
       avrt(kk) = 0.
       avrter(kk) = 0.
      enddo
      savrt = 0.
      savrtnox = 0.
      avrtc =0.
      savrter = 0.
      savrtnoxer = 0.
      avrtcer = 0.
      open(unit=5,file='corrrunlistnew.txt')
      read(5,'(a)') title
      do i=1,2509
c      do i=1,1720
c      do i=1,831
        read(5,'(i4,i2,10f7.3,2i8,2i10,f7.2,3f7.3,f7.1,6i6)') 
     >   irun,it,e0,ep0,the0,pp,thp,chrg,dt,
     >   tefe,tefp,current,irate_e,irate_p,
     >   ielclean,ielcleanp,bcm1,v1,v2,v3,elreal,
     >   ps1,ps2,ps3,ps4,ps5,p6

        kin(i,1)= abs(ep0)
        kin(1,2)= the0
        the0r = the0
        the0 = the0 * pi/180.
        ep0 = abs(ep0)
        thp_rad = thp * pi/180.

 
       if(it.ne.itlast.or.
cc     >   pp.ne.12345. .or.
     >   abs(pp - pplast).gt.0.02 .or.
     >   abs(thp-thplast).gt.0.1) then
         if(nrt.gt.0) then
          ikin = 0
          if(irun.le.6008 .or. irun.ge. 6554) then
           if(abs(abs(pplast) - 1.925).lt.0.01) ikin = 1
           if(abs(abs(pplast) - 2.485).lt.0.01) ikin = 2
           if(abs(abs(pplast) - 2.602).lt.0.01) ikin = 3
           if(abs(abs(pplast) - 3.339).lt.0.01) ikin = 4
           if(abs(abs(pplast) - 2.006).lt.0.01) ikin = 5
           if(abs(abs(pplast) - 2.575).lt.0.01) ikin = 6
           if(abs(abs(pplast) - 3.43).lt.0.01) ikin = 7
           if(abs(abs(pplast) - 4.79).lt.0.01 .and.
     >      irun.ge.5391) ikin = 8
          endif
          itp = itlast
          if(pplast.lt.0.) itp = itp + 3
          ith = int((thplast-5.9)/2.) + 1
c special case
          if(ikin.eq.5) ith = int((thplast - 12.9)/3.) + 1
          do kk=1,10
           avrt(kk) = avrt(kk) / avrter(kk)
           avrter(kk) = 1./sqrt(avrter(kk))
          enddo
          savrt = savrt / savrter
          savrter = 1./sqrt(savrter)
          savrtnox = savrtnox / savrtnoxer
          savrtnoxer = 1./sqrt(savrtnoxer)
          doit = 0
c          if(irun.lt.4212) doit=1
          if(irun.lt.4233) doit=1
          if(irun.gt.4275 .and. irun.lt.4313 ) doit=1
          if(irun.gt.4323.and.irun.lt.9000) doit=1
          if(doit.gt.0) then
           if(ikin.gt.0) then
            avrate(ikin,itp,ith) = avrt(1)
            avrateer(ikin,itp,ith) = avrter(1)
            savrate(ikin,itp,ith) = savrt
            savrateer(ikin,itp,ith) = savrter
            savratenox(ikin,itp,ith) = savrtnox
            savratenoxer(ikin,itp,ith) = savrtnoxer
           endif
           if(curmax/curmin.gt.1.4) write(16,'(''current scan'')')
           do kk=1,10
            chi = 0.
            df = 0.
            do k=1,nrt
             chi = chi + (rt(k,kk) - avrt(kk))**2 / rter(k,kk)**2
             df = df + 1.
            enddo
            totchi2(kk) = totchi2(kk) + chi
            totdf(kk) = totdf(kk) + df
            if(curmax/curmin.gt.1.4 .or. kk.eq.1) then
             write(16,'(8x,3i5,3f8.3,f4.0)') ikin,itp,ith,
     >        avrt(kk),avrter(kk),chi/df,df
             if(kk.eq.1 .and. chi .gt. df + 6.) write(16,
     >        '(''warning, big chisq'')')
            endif
           enddo ! loop over kk
           write(16,'(8x,3i5,4f8.3)') ikin,itp,ith,
     >       savrt,savrter,savrtnox,savrtnoxer
          endif
         endif
         write(16,'()')
        endif
c write out all individual files: no summing over runs
c start longg
        if( pp.ne.12345) then

c write out cnts
        close(unit=22)
       write(fname,
     >    '(''/group/c-sidis/bosted/Cntfiles/Cntdecay'',
     >     i4,''.txt'')') irunlast
c        write(fname,
c     >    '(''/group/c-sidis/bosted/Cnt0/Cntdecay'',
c     >     i4,''.txt'')') irunlast
        open(unit=22,file=fname)
        do icase=1,7
         sum3 = 0.
         do idist=1,22
          sum3 = sum3 + fstate(icase,idist,4)
         enddo
         do idist=1,22
          do ipart=1,3
           if(fstate(icase,idist,4).ne.0.) then
            fstate(icase,idist,ipart) =
     >       fstate(icase,idist,ipart) /
     >       fstate(icase,idist,4)
           endif
          enddo
          write(22,'(2i3,4f8.3)') icase,idist,
     >      fstate(icase,idist,4)/sum3,
     >      (fstate(icase,idist,ipart),ipart=1,3)
         enddo
        enddo
        do icase=1,8
         irr = icase*2
         normfac = avnfac(icase) / chrgtot / 
     >    max(1.,sum3tot(icase))
         do im=1,50
          cntsmcmmpi(im,irr) = cntsmcmmpi(im,irr) * normfac 
         enddo
         do izz=1,50
          cntsmcz(izz,irr) = cntsmcz(izz,irr) * normfac
         enddo
         do ipt =1,16
          do iphi=1,15
           do iz=1,20
             do iq2bin=1,2
              cntsmc(iq2bin,ipt,iphi,iz,irr) = 
     >         cntsmc(iq2bin,ipt,iphi,iz,irr)  * normfac
             enddo
           enddo
          enddo
         enddo
        if(icase.eq.2) then
         do ipt =1,16
          do iphi=1,15
           do iq2bin=1,2
            cntsemc(iq2bin,ipt,iphi,2) = 
     >       cntsemc(iq2bin,ipt,iphi,2)  * normfac
           enddo
          enddo
         enddo
        endif
        enddo
        close(unit=22)
       write(fname,
     >    '(''/group/c-sidis/bosted/Cntfiles/Cntdata'',
     >     i4,''.txt'')') irunlast
c        write(fname,
c     >    '(''/group/c-sidis/bosted/Cnt0/Cntdata'',
c     >     i4,''.txt'')') irunlast
        open(unit=22,file=fname)
        do im =1,50
         write(22,'(i3,8i5,7(f6.0,e12.4),2e12.4)') im, 
     >     (cntsmmpi(im,irr),irr=1,8),
     >     (cntsmcmmpi(im,irr),irr=1,16)
        enddo
        do izz=1,50
         write(22,'(i3,8i5,7(f6.0,e12.4),2e12.4)') izz,
     >     (cntsz(izz,irr),irr=1,8),
     >     (cntsmcz(izz,irr),irr=1,16)
        enddo

        do iq2bin=1,2
        do iz=1,20
         do iphi=1,15
          do ipt=1,16
           sum3 = 0.
           do irr=1,2
            sum3 = sum3 + cnts(iq2bin,ipt,iphi,iz,irr)
            sum3 = sum3 + cnts(iq2bin,ipt,iphi,iz,irr+6)
           enddo
           sum3 = sum3 + cntsmc(iq2bin,ipt,iphi,iz,1)
           sum3 = sum3 + cntsmc(iq2bin,ipt,iphi,iz,3)
           if(sum3.gt.0) then
            write(22,'(i2,3i3,8i5,f6.0,e12.4,f6.0,e12.4,
     >       3f7.1,f7.3,3f7.1,f7.3,
     >       5(f6.0,e12.4),2e12.4)') 
     >       iq2bin,ipt,iphi,iz,
     >       (cnts(iq2bin,ipt,iphi,iz,irr),irr=1,8),
     >       (cntsmc(iq2bin,ipt,iphi,iz,irr),irr=1,4),
     >       (avkin(iq2bin,ipt,iphi,iz,irr)/
     >        max(1.,cntsmc(iq2bin,ipt,iphi,iz,1)),irr=1,8),
     >       (cntsmc(iq2bin,ipt,iphi,iz,irr),irr=5,16)
           endif
          enddo
         enddo
        enddo
        enddo
c coincidence spectra
        close(unit=22)
        write(fname,
     >    '(''/group/c-sidis/bosted/Cntfiles/Cntct'',
     >     i4,''.txt'')') irunlast
c        write(fname,
c     >    '(''/group/c-sidis/bosted/Cnt0/Cntct'',
c     >     i4,''.txt'')') irunlast
        open(unit=22,file=fname)
        do iz=1,20
         do kk=1,40
          write(22,'(4i6)') (cntsct(iz,kk,im),im=1,4)
         enddo
        enddo
c exclusive pion counts
        close(unit=22)
        write(fname,
     >    '(''/group/c-sidis/bosted/Cntfiles/Cnte'',
     >     i4,''.txt'')') irunlast
c        write(fname,
c     >    '(''/group/c-sidis/bosted/Cnt0/Cnte'',
c     >     i4,''.txt'')') irunlast
        open(unit=22,file=fname)
        do iq2bin=1,2
         do iphi=1,7
          do ipt=1,4
           pir = cntse(iq2bin,ipt,iphi,1)
           piacc = cntse(iq2bin,ipt,iphi,2)
           rex = (pir - piacc/4.) / chrg / dt / 
     >      tefe / tefp / corr3 / corr / corr2
           rexer = sqrt(pir + piacc/16.)/chrg/dt/
     >      tefe / tefp / corr3 / corr / corr2
           if(cntsemc(iq2bin,ipt,iphi,2).gt.0.) then
            rex =  rex/cntsemc(iq2bin,ipt,iphi,2)
            rexer = rexer/cntsemc(iq2bin,ipt,iphi,2)
           endif
           write(22,'(i2,2i3,2i5,f6.0,2e12.4,2f8.3)') 
     >       iq2bin,ipt,iphi,
     >       (cntse(iq2bin,ipt,iphi,irr),irr=1,2),
     >       (cntsemc(iq2bin,ipt,iphi,irr),irr=1,3),
     >       rex, rexer
          enddo
         enddo
        enddo
        do iq2bin=1,2
         do iphi=1,7
          do ipt=1,4
           write(22,'(i2,2i3,8i5/8x,8i5/8x,8i5)') 
     >       iq2bin,ipt,iphi,
     >       (cntseh(iq2bin,ipt,iphi,im,1),im=1,8),
     >       (cntseh(iq2bin,ipt,iphi,im,2),im=1,8),
     >       (cntsemch(iq2bin,ipt,iphi,im),im=1,8)
          enddo
         enddo
        enddo
         DO IPT=1,16
         DO IPHI=1,15
          DO IQ2BIN=1,2
           CNTSE(IQ2BIN,IPT,IPHI,1)=0.
           CNTSE(IQ2BIN,IPT,IPHI,2)=0.
           DO IM=1,8
            CNTSEH(IQ2BIN,IPT,IPHI,IM,1)=0.
            CNTSEH(IQ2BIN,IPT,IPHI,IM,2)=0.
            CNTSEMCH(IQ2BIN,IPT,IPHI,IM)=0.
           ENDDO
           DO IRR=1,3
            CNTSEMC(IQ2BIN,IPT,IPHI,IRR)=0.
           ENDDO
          ENDDO
         ENDDO
        ENDDO
        do iz=1,20
         do kk=1,40
          do jj=1,4
            cntsct(iz,kk,jj) = 0
          ENDDO
         ENDDO
        ENDDO
        DO IPT =1,16
         DO IPHI=1,15
          DO IZ=1,20
           DO IQ2BIN=1,2
            DO IRR=1,8
             CNTS(IQ2BIN,IPT,IPHI,IZ,IRR)=0
            ENDDO
            DO IRR=1,8
             AVKIN(IQ2BIN,IPT,IPHI,IZ,IRR)=0.
            ENDDO
            DO IRR=1,16
             CNTSMC(IQ2BIN,IPT,IPHI,IZ,IRR)=0
            ENDDO
           ENDDO
          ENDDO
         ENDDO
        ENDDO
        DO IPT=1,50
         DO IRR=1,8
          CNTSMMPI(IPT,IRR)=0
          CNTSZ(IPT,IRR)=0
         ENDDO
         DO IRR=1,16
          CNTSMCMMPI(IPT,IRR)=0
          CNTSMCZ(IPT,IRR)=0
         ENDDO
        ENDDO
        write(181,182) irunlast,itlast,e0last,ep0last,
     >   the0last,pplast,thplast,chrgtot,
     >   dtav/chrgtot,
     >   tefeav/chrgtot,
     >   tefpav/chrgtot,
     >   corrav/chrgtot,
     >   corr2av/chrgtot,
     >   corr3av/chrgtot
 182    format(i4,i2,5f7.3,f8.2,6f7.3)
        chrgtot = 0.
        tefeav = 0.
        tefpav = 0.
        dtav = 0.
        corr3av = 0.
        corrav = 0.
        corr2av = 0.
        do irr=1,8
         avnfac(irr) = 0.
         sum3tot(irr) = 0.
        enddo

         write(16,'(1x,''starting it, pp, thp='',i2,f7.2,f7.1)')
     >    it,pp,thp
        write(16,'(''  rn  cur  rate   err simm/d  w/ex'',
     >     ''   lt tefe tefe  r/r cor3  cor cor2 hsef psef'')')
         itlast = it
         thplast = thp
         pplast = pp
         e0last = e0
         ep0last = ep0
         the0last = the0r
         nrt = 0
         do kk=1,10
           avrt(kk) = 0.
           avrter(kk) = 0.
         enddo
         savrt = 0.
         savrter = 0.
         savrtnox = 0.
         savrtnoxer = 0.
         curmin=1000.
         curmax=0.
! end longg section on starting new kin.
        endif

c tried this, but not good
c      dpshmslo = -12.
c      if(irun.gt.4400 .and. irun.le.5334 ) dpshmslo = -19.5
c      if(irun.ge.7871 .and. irun.lt.8300 ) dpshmslo = -19.5

c apply BCM4A correction (from code hmsf)
c        corr = 0.995 + 0.035 * (log(60.) - log(current)) /
c     >                         (log(60.) - log(2.))
c        chrg = chrg * corr
c use bcm1 for spring18 csv runs
c        if(irun.ge.4860 .and. irun.le.5350) chrg = bcm1
c changed to use bcm1 for all runs now
c fix stange gain change
        if(irun.ge.3651.and.irun.lt.4400) bcm1 = bcm1 * 0.94
        if(current.le.60.) then
          bcmcorr = 1.00 + 0.045 * (log(60.) - log(current)) /
     >                          (log(60.) - log(2.))
        else
         bcmcorr = 1. + 0.010 * (current - 60.) / 25.
        endif
        bcm1 = bcm1 * bcmcorr
c override bcm4a with bcm1
        chrg = bcm1

c change dt from percent to fraction
       dt = dt/100.

c aerogel tray. 2= index 1.015, 1= index 1.011
       aerotray = 2
       if(irun.ge.4965.and.irun.le.5378) aerotray=1
       if(irun.ge.7940.and.irun.le.8356) aerotray=1
       xmaxaero = 100.
       ymaxaero = 100.
       aeromin = 2.5
       if(aerotray.eq.1) then
        xmaxaero = 40.
        ymaxaero = 25.
        aeromin = 2.5
       endif
c get central mmpi
! get kinematic vectors
! positive yptar makes theta smaller in HMS
        ep = ep0
        dphie = 0.
        dthe = 0.
        ppi = abs(pp) 
        dthp = 0.
        dphip = 0.
	 call physics_angles(the0, phi0e, dthe, dphie,
     >     ep,p_x(1),p_y(1),p_z(1))
	 call physics_angles(thp_rad, phi0p, dthp, dphip, 
     >     ppi,p_x(2),p_y(2),p_z(2))
        Empi = e0 + amp - ppi - sqrt(ampi**2 + ep**2)
        mmpi2_cent_h = 
     >     (Empi)**2 - 
     >     (p_x(1) + p_x(2))**2 -
     >     (p_y(1) + p_y(2))**2 -
     >     (p_z(1) + p_z(2)- e0)**2
c this is for e- in shms, pi- in hms
        Empi = e0 + amp - ep - sqrt(ampi**2 + ppi**2)
        mmpi2_cent_s = 
     >     (Empi)**2 - 
     >     (p_x(1) + p_x(2))**2 -
     >     (p_y(1) + p_y(2))**2 -
     >     (p_z(1) + p_z(2)- e0)**2
c        write(6,'(''mmpi2='',i5,2f8.2)') irun,mmpi2_cent_h,mmpi2_cent_s

! rf spacing in nsec. Only for Fall (not hooked up in Spring)
! Right now, the master oscillator frequency is 499001553.45 Hz. We get every other
! bunch, so the bunch space should be 4.008 ns like you found below.
! see fort.44 for check: seems good for whole run
! MOFC1FREQ - this gives the absolute frequency in Hz (EPICS)
! MOFC1DELTA - this give the deviation from exactly 499 MHz

       rff = 4.000 * 1.0020

       kin(i,3)=pp
       kin(i,4)=thp
       itsv(i)=it

c runs with low HMS elclean / chrg (copied from hms.f)
       ikin=1
 
c Skip some problem runs
        if(ikin.eq.1 .and. 
c    >    ((irun.ge.8180 .and. irun.le. 8190) .or. 
c    >     (irun.ge.5240 .and. irun.le. 5250) .or.
c    >     (irun.ge.5290 .and. irun.le. 5310) .or. 
c    >     (irun.ge.8040 .and. irun.le. 8050) .or. 
c    >     (irun.ge.6490 .and. irun.le. 6520))  .and. 
c     >     irun.ge.3300.and.irun.le.3450.and.
c     >     irun.ge.6080.and.irun.le.6085.and.
c     >     irun.eq.3716 .and.
c     >     thp.gt.18.0 .and.
c     >     it.eq.1 .and.pp.gt.0.0 .and.
c     >   mmpi2_cent_h.lt.1.5 .and.
c    >    10*(irun/10).eq.irun .and. irun.lt.4400 .and.
c    >      (irun.eq.3478.or.irun.eq.5603.or.irun.gt.8226).and.
c current scans
c    >  ((irun.ge.3470 .and.  irun.le. 3472) .or.
c    >   (irun.ge.3517 .and.  irun.le. 3519) .or.
c    >   (irun.ge.3982 .and.  irun.le. 3984) .or.
c    >   (irun.ge.3435 .and.  irun.le. 3437) .or.
c    >   (irun.ge.4051 .and.  irun.le. 4054) .or.
c    >   (irun.ge.5600 .and.  irun.le. 5617) .or.
c    >   (irun.ge.3460 .and.  irun.le. 3478) .or.
c    >   (irun.ge.4882 .and.  irun.le. 4890) .or.
c    >   (irun.ge.4910 .and.  irun.le. 4912) .or.
c    >   (irun.ge.4965 .and.  irun.le. 4975) .or.
c    >   (irun.ge.4982 .and.  irun.le. 4998) .or.
c    >   (irun.ge.5040 .and.  irun.le. 5066) .or.
c    >   (irun.ge.4965 .and.  irun.le. 4975) .or.
c    >   (irun.ge.5128 .and.  irun.le. 5220) .or.
c    >   (irun.ge.5367 .and.  irun.le. 5378) .or.
c    >   (irun.ge.5415 .and.  irun.le. 5436) .or.
c    >   (irun.ge.5437 .and.  irun.le. 5440) .or.
c    >   (irun.ge.5634 .and.  irun.le. 5636) .or.
c    >   (irun.ge.5671 .and.  irun.le. 5727) .or.
c    >   (irun.ge.5751 .and.  irun.le. 5754) .or.
c    >   (irun.ge.5963 .and.  irun.le. 5965) .or.
c    >   (irun.ge.6129 .and.  irun.le. 6135) .or.
c    >   (irun.ge.6429 .and.  irun.le. 6433) .or.
c    >   (irun.ge.6459 .and.  irun.le. 6464) .or.
c    >   (irun.ge.8039 .and.  irun.le. 8084) .or.
c    >   (irun.ge.5687 .and.  irun.le. 5689) .or.
c    >   (irun.gt.8226)) .and.
c     >   irun.ge.3500  .and.irun.le.3520.and.
c Kaon 10.6 runs
c     >   irun.ge.4865  .and.irun.le.5360.and.
c     >   irun.ge.5460  .and.irun.le.5472.and.
c     >   irun.ge.5775  .and.irun.le.5490.and.
c     >   irun.gt.5400 .and. irun.le.9000 .and.
c     >   irun.ge.7500  .and.irun.le.7900.and.
c     >  (irun.eq.3424 .or. irun.eq.3442 .or.
c     >   irun.eq.3472 .or.irun.eq.3842 .or.
c     >   irun.eq.3929 .or.irun.eq.3947 .or.
c     >   irun.eq.4070 .or.irun.eq.5406 .or.
c     >   irun.eq.5411 .or.irun.eq.5507.or.
c     >    (irun.gt.7600 .and. irun.lt.7850) .or.
c     >   irun.eq.5993.or. irun.eq.7604) .and.
c     >  pp.gt.4.0.and.irun.lt.6100.and.
c     >    irun.lt.3440 .and.
c     >    (irun.eq.3900.or.irun.eq.6025) .and.
c     >   ((irun.ge.3800 .and.irun.le.3830) .or. 
c     >    irun.ge.5400 .and. irun.le.6000 .and.
c     >    pp.lt.4.2 .and. irate_p.lt.150000 .and.
c     >   irun.ge.3800 .and.irun.le.3830.and.
c     >   irun.ge.5820 .and.irun.le.5830.and.
c     >   it.eq.1 .and.
c     >   abs(pp).lt. 3.0.and.
c     >   ((irun.ge.6192 .and.irun.le.6293) .or.
c     >    (irun.ge.3421 .and.irun.le.3425 .and.it.eq.1)).and.
c short run: low rate
     >    irun.ne.4188 .and.
c low rate for no apparant reason
     >    irun.ne.4324 .and.
c looks like BCM problem KLT10
     >    irun.ne.4892 .and.irun.ne.4893.and.
c bad runs in csv 2019
     >    irun.ne.7599 .and.
     >    irun.ne.7600 .and.
     >    irun.ne.7678 .and.
c short low rate run klt
     >    irun.ne.7914 .and.
     >    irun.ne.7915 .and.
     >    irun.ne.7917 .and.
     >    irun.ne.7980 .and.
c positron runs
     >   (irun.lt.4212.or.irun.gt.4231)
c elastic runs
     >   .and.(irun.lt.6009.or.irun.gt.6017)
     >   ) then
c Original skim files
c        write(fname,
c     >    '(''/work/hallc/sane/bosted/ptc/Skim'',i4,''.txt'')') irun
c new ones
        write(fname,
     >    '(''/work/hallc/c-sidis18/bosted/Skimfiles/Skim'',
     >     i4,''.txt'')') irun
        OPEN(UNIT=9,FILE=FNAME)

        DO K=1,15
         DPPHISTRUN(K)=0.
         YTHISTRUN(K)=0.
        ENDDO
        DO K=1,50
         DO KK=1,10
          CTPIHW(K,KK)=0
         ENDDO
        ENDDO
        PIR = 0.
        PIACC = 0.
        PIACCL = 0.
        PIACCH = 0.
        DO IM=1,20
         PIRM(IM) = 0
         PIACCM(IM) = 0
         PIRZ(IM) = 0
         PIACCZ(IM) = 0
        ENDDO
        DO K=1,6
         DO KK=1,20
          PIRK(K,KK)=0.
          PIACCK(K,KK)=0.
          PIRKS(K,KK)=0.
          PIACCKS(K,KK)=0.
         ENDDO
        ENDDO
        DO K=1,40
         CTPIH(K)=0
         CTPIHF(K)=0
         CTPIHFA(K)=0
         CTKH(K)=0
         CTPHF(K)=0
         CTKHF(K)=0
         CTEH(K)=0
        ENDDO
        DO K=1,14
         DO KK=1,21
          EPRH(K,KK)=0
         ENDDO
        ENDDO
        DO K=1,12
         DO KK=1,16
          CEREH(K,KK)=0
          CEREPIH(K,KK)=0
         ENDDO
        ENDDO
        DO K=1,16
         NGH(K)=0
         SHH(K)=0
         prh(k,1)=0
         prh(k,2)=0
         SHEH(K)=0
         SH2H(K)=0
         AEROH(K)=0
         AEROHp(K)=0
         DO KK=1,19
          YTH(KK,K)=0
         ENDDO
        ENDDO
        DO K=1,20
         BETAH(1,K)=0
         BETAH(2,K)=0
        enddo
        DO K=1,21
         CTRFRUN(K,1)=0
         CTRFRUN(K,2)=0
        enddo
        do k=1,20
         DO KK=1,10
          CTRFRN(KK,K,1)=0
          CTRFRN(KK,K,2)=0
         ENDDO
         DO KK=1,4
          XYAEROH(KK,K)=0.
          XYAEROH(KK,21)=0.
         ENDDO
         DO KK=1,12
          MMPIH(K,KK)=0
          MMPH(K,KK)=0
         ENDDO
         PHIH(K)=0
         MEEH(K)=0
        ENDDO
        do k=1,400
         MMPHF(K)=0
         MMPHETA(K)=0
        enddo
        DO ii =1,8
         DO j=1,22
          DO k=1,4
           fstate(ii,j,k)=0.
          ENDDO
         ENDDO
        ENDDO

        DO IHEL=1,3
         do j=1,2
          SUMCNT(IHEL,j)=0
          SUMCNTSF(IHEL,j)=0
         enddo
        ENDDO

        DO J=1,7654321
C        DO J=1,10000
C         READ(9a)',end=10,err=10) string
         read(9,'(a)',end=10) string
         if(j.lt.2) write(6,'(a)') string(1:60)
         read(string,*) ctpi,ctk,ctp
         if(ctpi.eq.0.0 .and.ctk.eq.0. .and. ctp.eq.0.) then
c            write(6,'(''finished'',i5,i7)') irun,j
            read(9,*) sum1,sum2,sum3,sum4
c override track effeciency SHMS with new one
c            tefp = sum4 / sum2
            tefp = sum4 / sum2
            if(tefp.lt.0.80) write(16,'(''big bad tefp='',f8.3)') tefp
            read(9,*) sum1,sum2,sum3,sum4
c override track effeciency HMS with new one
c            tefe = sum4 / sum2
c changed to sum3 because pbeta is mostly 0. for some run periods (kaonlt)
            tefe = sum4 / sum2
c skkip 200 lines
            do jj=1,200
             read(9,'(a)',end=10,err=10) string
            enddo
c get percentage of cointimeraw within +/-100 nsec
            read(9,*) sum1,sum2,sum3
c            timecorr = 1. - sum3
! Get efficiency of pgdsc
            read(9,*) i1,i2,pgdsceff
            goto 10
         endif ! end check on end of file
! read in an event
! cahnged order of dthe, dphie to match PeterB.C
! dth is xprar, dphi is yptar
         if(irun.lt.4400) then
          read(string,*,end=999) ctpi,ctk,ctp,imps,ineg,ipos,
     >     dpe,dthe,dphie,dpp,dthp,dphip,cere,pre,she,
     >     cerng,ceraero,cerhg,prp,shp,hgtry,pgtry,pbeta,
     >     hdcx,hdcxp,hdcy,hdcyp,pdcx,pdcxp,pdcy,pdcyp,
     <     prf,hrf,pfpt,hfpt,pztar,hztar,gdsc
     >     ,goodhodoh,goodhodop,hbeta,ctpinew,aerot,aeront,evtype
          i1 = min(15,max(1,int((ctpi+30.)/4.)+1))
          i2 = min(6,max(1,int(float(irate_p)/1.e5)+1))
          if(evtype.lt.3) then
           ctev(i1,i2,1) = ctev(i1,i2,1) + ps4
          else
           ctev(i1,i2,2) = ctev(i1,i2,2) + 1
          endif

c beam helicity
          if(irun.lt.5360) then
           ihel = 3 ! undefined
           if(ipos.eq.1 .and. imps.eq.ineg) then
            ihel = max(1,imps) ! 0->1 2->2
           endif
          else
           if(imps.eq.1) ihel=3
           if(imps.eq.0) ihel=1
           if(imps.eq.2) ihel=2
           if(imps.lt.0 .or. imps.gt.2) then
            ihel=3
            write(6,'(''error hel'',3i3)') imps,ipos,ineg
           endif
          endif

c corresction for HMS set too low by 0.16% in sp18 for 5.27 setting
c and also an angle correction
c add the 0.4%
c Now new values in run list
c          if(ep0.gt.5.2) then 
c           dpe = dpe - 0.56
c          endif
c          dpp = dpp - 0.4

c vertical angle offset for all of spring 18
c Carlos found 2.85 mr, so my value is similar
c (comes from making cos(phi) dist. even)
           dthe = dthe + 0.0027

c this is for fall18, spring19
         else
          evtype=4
          read(string,*,end=999) ctpi,ctk,ctp,imps,ineg,ipos,
     >     dpe,dthe,dphie,dpp,dthp,dphip,cere,pre,she,
     >     cerng,ceraero,cerhg,prp,shp,hgtry,pgtry,pbeta,
     >     hdcx,hdcxp,hdcy,hdcyp,pdcx,pdcxp,pdcy,pdcyp,
     <     prf,hrf,pfpt,hfpt,pztar,hztar,gdsc
     >     ,goodhodoh,goodhodop,hbeta,ctpinew,aerot,aeront

c also shift for other runs
           dthe = dthe + 0.0027

c beam helicity
          if(irun.lt.5360) then
           ihel = 3 ! undefined
           if(ipos.eq.1 .and. imps.eq.ineg) then
            ihel = max(1,imps) ! 0->1 2->2
           endif
          else
           if(imps.eq.1) ihel=3
           if(imps.eq.0) ihel=1
           if(imps.eq.2) ihel=2
           if(imps.lt.0 .or. imps.gt.2) then
            ihel=3
            write(6,'(''error hel'',3i3)') imps,ipos,ineg
           endif
          endif
         endif
c         goodhodoh = 1.
c         goodhodop = 1.
c         hbeta = 1.

c change pr to not be normalized by pp
c         prp= prp * abs(pp)

c correction to hms xptar. Data were processed with
c an offset of -5 mr. Carlos found we should use
c +3 mr on average with new DC. In reality, this
c correction should depend on hdelta as
c it is probably due to DC alignment, but for
c now lets just do a constant offset of 5 mr
c reran with zero offset for spring18, but -5 by mistake for
c fall18

         k = int((dpe+12.)/24.*12.)+1
         k = max(1,min(k,12))
         kk = max(1,min(15,int(cere/1.5)+1))
         if(she.gt.0.7) then
          cereh(k,kk) = cereh(k,kk)+1
          cereh(k,16) = cereh(k,16)+1
         else
          cerepih(k,kk) = cerepih(k,kk)+1
          cerepih(k,16) = cerepih(k,16)+1
         endif

         kk = max(1,min(15,int(she*10.)+1))
         sheh(kk) = sheh(kk) + 1
         sheh(16) = sheh(16) + 1

c SHMS dipole exit
         pexit=1
         xexit = pdcx -307. * pdcxp
         yexit = pdcy -307. * pdcyp
         crad = 23.81
         voffset = crad - 24.035
         hwid = 11.549/2.
         if(abs(yexit) .lt. hwid) then
           if(abs(xexit) .gt.  (crad + voffset)) pexit=0 
         else
           if ( yexit .ge. hwid) then
            if ( (xexit - voffset)**2 + 
     >           (yexit - hwid)**2 .gt.
     >            crad**2) pexit=0;
           endif
           if ( yexit .le. -1.*hwid) then
            if ( (xexit - voffset)**2 + 
     >           (yexit + hwid)**2 .gt.
     >            crad**2) pexit=0;
           endif
         endif
c         if(j.lt.10000.and.pexit.eq.0) write(6,'(i2,6f8.3)') 
c     >    pexit,xexit,yexit,pdcx,pdcxp,pdcy,pdcyp

c HMS dipole exit
         xexit = hdcx - 148. * hdcxp
         yexit = hdcy - 148. * hdcyp
         hexit=1
         if ( (xexit - 2.8)**2 + 
     >        (yexit)**2 .gt. 46.607**2) hexit=0;
c         if(j.gt.10000.and.hexit.eq.0) write(6,'(i3,6f8.3)') 
c     >    hexit,xexit,yexit,hdcx,hdcxp,hdcy,hdcyp


c option to  using corrected dpp. Doesn't help at least
c for positive dpp, so no need. Haven't checked dpp<-12
c because cannot find suitable run.
c         call shmscorr(pdcx,pdcxp,pdcy,pdcyp,dpp,dppcorr)
c         dpp = dppcorr

c fill in acceptance array and fill in ok flags
         okhms = 0
         if(abs(dpe).lt.13 .and.
     >    abs(dthe).lt.0.100 .and. 
     >    abs(dphie).lt.0.030) then
          ith = int((dthe + 0.100) / 0.200 * 70)+1
          iphi = int((dphie + 0.030) / 0.060 * 30.)+1
          ipt = int((dpe + 13.) / 26. * 30.)+1
          accep(1,ipt,ith,iphi) = accep(1,ipt,ith,iphi) + 1
          okhms = okacc(1,ipt,ith,iphi)
          if(okhms.eq.1) accep(5,ipt,ith,iphi) = 
     >     accep(5,ipt,ith,iphi) + 1
! try making adjustments to DC 
! this takes a long time, so limit number of event per run
          if(j.lt.2) then
           fry = 0.001
	   call mc_hms_recon (hdcx,hdcxp,hdcy,hdcyp,fry,
     >      dpz, dphiz, dthz, yz)
           if(j.lt.1) write(6,'(''hms z'',8f8.4)') 
     >      dpe, dpz, dthe, dthz, dphie, dphiz
           do k=1,3
            do kk=1,3
c             xoff = -0.5 + 0.5 * (k-1)
             xoff = 0.
             fry = -0.101 + 0.1*(k-1)
             xpoff = -0.001 + 0.0010 * (kk-1)
	     call mc_hms_recon (hdcx+xoff,hdcxp+xpoff,hdcy,hdcyp,
     >         fry,dpn, dphin, dthn, yn)
             if(j.lt.1)  write(6,'(2i3,8f7.3)') 
     >         k,kk,dpz, dpn, dthz, dthn, dphiz, dphin, yz, yn
             if(abs(dpn).lt.13 .and.
     >        abs(dthn).lt.0.100 .and. 
     >        abs(dphin).lt.0.030) then
              ith = int((dthn + 0.100) / 0.200 * 70)+1
              iphi = int((dphin + 0.030) / 0.060 * 30.)+1
              ipt = int((dpn + 13.) / 26. * 30.)+1
              kkk = 10 + (k-1) * 3 + kk
              accep(kkk,ipt,ith,iphi) = 
     >          accep(kkk,ipt,ith,iphi) + 1
             endif
            enddo
           enddo
          endif ! j
         endif ! dpe, ...

         okshms = 0
         if(abs(dpp-5.).lt.17 .and.
     >    abs(dthp).lt.0.100 .and. 
     >    abs(dphip).lt.0.030) then
          ith = int((dthp + 0.100) / 0.200 * 70)+1
          iphi = int((dphip + 0.030) / 0.060 * 30.)+1
          ipt = int((dpp + 12.) / 34. * 30.)+1
          accep(3,ipt,ith,iphi) = accep(3,ipt,ith,iphi) + 1
          okshms = okacc(2,ipt,ith,iphi)
          if(okshms.eq.1) accep(7,ipt,ith,iphi) = 
     >     accep(7,ipt,ith,iphi) + 1
         endif
c use same th and phi range for dpp<-12 as for -12
         if(dpp.lt.-12 .and.
     >    abs(dthp).lt.0.100 .and. 
     >    abs(dphip).lt.0.030) then
          ith = int((dthp + 0.100) / 0.200 * 70)+1
          iphi = int((dphip + 0.030) / 0.060 * 30.)+1
          ipt = 1
          okshms = okacc(2,ipt,ith,iphi)
         endif

c correct to make pion peak at 0 nsec (bin 10.5)
         if(irun.le.3418 ) ctpi = ctpi + 0.5
         if(irun.gt.3418.and.irun.le.3428 ) ctpi = ctpi - 0.30
         if(irun.ge.3429.and.irun.lt.3550 ) ctpi = ctpi - 0.12
         if(irun.ge.3550 .and. irun.lt.4253) ctpi = ctpi -0.05
         if(irun.ge.4253 .and. irun.lt.4400) ctpi = ctpi +0.12 
         if(irun.ge.4400 .and. irun.lt.5304) ctpi = ctpi -0.6 
         if(irun.ge.5304 .and. irun.lt.5360) ctpi = ctpi -0.10 
         if(irun.ge.5360 .and. irun.lt.6156) ctpi = ctpi - 0.0 
         if(irun.ge.6156 .and. irun.lt.7592 ) ctpi = ctpi- 0.8 
         if(irun.ge.7592 .and. irun.lt.7665 ) ctpi = ctpi- 0.20 
         if(irun.gt.7665) ctpi = ctpi- 0.00

c try making corr. a bit bigger to account for 2 m difference 
c between fp and center of hodo
         ctk = ctk * 1.05 + ctpi
         ctp = ctp * 1.05 + ctpi

! histograms of p and ytarg by run
         if(dpp.gt.-25. .and. dpp.lt. 45.) then
          k = int((dpp+25.) / 70. * 14.)+1
          dpphistrun(k)=dpphistrun(k)+1
         endif
         dpphistrun(15)=dpphistrun(15)+1

         if(pgtry.gt.-20. .and. pgtry.lt. 20.) then
          k = int((pgtry+20.) / 40. * 14.)+1
          ythistrun(k)=ythistrun(k)+1
         endif
         ythistrun(15)=ythistrun(15)+1

! histograms of p, th, phi
         if(dpp.gt.-50. .and. dpp.lt. 50.) then
          k = int(dpp+50.)+1
          dpph(k)=dpph(k)+1
         endif
         if(dpp.gt.-20. .and. dpp.lt. 20.) then
          k = int((dpe+20.) * 2.5 )+1
          dpeh(k)=dpeh(k)+1
         endif
         if(dthp.gt.-0.05 .and. dthp.lt. 0.0499) then
          k = int(dthp*1000. + 50.)+1
          dthph(k)=dthph(k)+1
         endif
         if(dthe.gt.-0.05 .and. dthe.lt. 0.0499) then
          k = int(dthe*1000. + 50.)+1
          dtheh(k)=dtheh(k)+1
         endif
         if(dphip.gt.-0.05 .and. dphip.lt. 0.0499) then
          k = int(dphip*1000. + 50.)+1
          dphiph(k)=dphiph(k)+1
         endif
         if(dphie.gt.-0.10 .and. dphie.lt. 0.0999) then
          k = int(dphie*500. + 50.)+1
          dphieh(k)=dphieh(k)+1
         endif
! beta histograms (no cuts)
         k = min(20,max(1,int((hbeta-0.9)/0.01)))
         betah(1,k) = betah(1,k)+1
         k = min(19,max(1,int((pbeta-0.8)/0.4*20)))
         betah(2,k) = betah(2,k)+1
         betah(2,20) = betah(2,20)+1

! get kinematic vectors
! positive yptar makes theta smaller in HMS
        the = the0 - dthe
        ep = ep0 * (1. + dpe/100.)
c wrong
c        p_x(1) =  ep * sin(the) * cos(dphie)
c        p_y(1) = -ep * sin(the) * sin(dphie)
c        p_z(1) =  ep * cos(the) 

! positive yptar makes theta bigger in SHMS
        ppi = abs(pp) * (1. + dpp/100.)
        thpi = thp * pi/180. + dthp
c wrong
c        p_x(2) = -ppi * sin(thpi) * cos(dphip)
c        p_y(2) = -ppi * sin(thpi) * sin(dphip)
c        p_z(2) =  ppi * cos(thpi) 
! above is wrong: this is right:

	 call physics_angles(the0, phi0e, dthe, dphie,
     >     ep,p_x(1),p_y(1),p_z(1))

c xxx for test
c         dthp = dthp + 0.005

	 call physics_angles(thp_rad, phi0p, dthp, dphip, 
     >     ppi,p_x(2),p_y(2),p_z(2))

c for test. No this is wrong
c	 call physics_angles(the0, phi0e, dphie, dthe,
c     >     ep,p_x(1),p_y(1),p_z(1))
c	 call physics_angles(thp_rad, phi0p, dphip, dthp,  
c     >     ppi,p_x(2),p_y(2),p_z(2))

        Empi = e0 + amp - ep - sqrt(ampi**2 + ppi**2)
        Emk = e0 + amp - ep - sqrt(amk**2 + ppi**2)
        Emp = e0 + amp - ep - sqrt(amp**2 + ppi**2)

        w2=(e0 + amp - ep)**2 - p_x(1)**2 - p_y(1)**2 - (p_z(1)-e0)**2
        w=0.
        if(w2.gt.0.) w = sqrt(w2)

        mmpi2 = 
     >     (Empi)**2 - 
     >     (p_x(1) + p_x(2))**2 -
     >     (p_y(1) + p_y(2))**2 -
     >     (p_z(1) + p_z(2)- e0)**2

        if(mmpi2.gt.0.8 .and. mmpi2 .lt. 1.0) then
         kk = min(20, max(1, 
     >    int( 20. * (p_x(1) + p_x(2) + 0.5)) + 1))
         pxdist(kk) = pxdist(kk)+1
         do kkk=1,7
          sum1 = ppi * 0.005 * (kkk-4)  
          kk = min(20, max(1, 
     >     int( 20. * (p_y(1) + p_y(2) + sum1 + 0.5)) + 1))
          pydist(kk,kkk) = pydist(kk,kkk)+1
         enddo
        endif
c look for radiated ep elastic events (proton in shms)
         if(ceraero.lt.0.5 .and.
     >      abs(ctp).lt.2.0 .and.
     >      goodhodoh .gt. 0. and.
     >      goodhodop .gt. 0. .and.
     >      cere.gt.2. .and. she.gt.0.7) then
          pt2 = (p_x(1) + p_x(2))**2 +
     >          (p_y(1) + p_y(2))**2
          e0_e = abs(ep)/(1. - 2.*abs(ep)*sin(the/2.)**2/amp)
          eprot = sqrt(amp**2 + ppi**2)
cx this is wrong! (cant use thpi!)
          e0_p  = amp * (Eprot - amp)  / (amp - Eprot +  ppi*cos(thpi))
          if(pt2.lt.0.001 .and. abs(e0_e - e0_p).lt.1.) then
c            write(6,'(6f7.3)') dpp,pt2, e0_e, e0_p, e0_e - e0_p
            k = min(14,max(1,int((dpp+25)/5.) + 1))
            kk = min(20,max(1,int((e0_e - e0_p + 1.)/0.1)+1))
            eprh(k,kk) = eprh(k,kk) + 1
            eprh(k,21) = eprh(k,21) + 1
           endif
         endif

! aero and sh histos
         xaero = pdcx + 231. * pdcxp
         yaero = pdcy + 231. * pdcyp

! Cuts on delta regions to be used and also electron PID
         if(dpe.gt.dphmslo .and. dpe.lt.dphmshi .and.
     >      dpp.gt.dpshmslo .and. dpp.lt.dpshmshi .and.
     >      abs(dthe).lt.dthhms .and.
     >      abs(dthp).lt.dthshms .and.
     >      abs(dphie).lt.dphihms .and.
     >      abs(dphip).lt.dphishms .and.
c     >       okhms.eq.1 .and. okshms.eq.1 .and.
     >       pexit.eq.1 .and. hexit.eq.1.and. 
c add cuts at aerogel xy
     >      abs(xaero).lt.xmaxaero .and. abs(yaero).lt.ymaxaero.and.
c requiring pexti and hexit only lowers rate by 0.1%
c don't require cere for spring18
     >      (cere.gt.1.1 .or. irun.lt.4400) .and.
     >       she.gt.0.65) then
c         write(6,'(i4,i6,f8.4)') irun,j,dthe
         
          if(ceraero.gt.2.5) then
           k = min(15,max(1,int(shp/0.1)+1))
           shh(k) = shh(k)+1
           shh(16) = shh(16)+1
           if(cerhg.gt.4.) then
            sh2h(k) = sh2h(k)+1 
            sh2h(16) = sh2h(16)+1
           endif
           k = min(15,max(1,int(prp/0.01)+1))
           if(shp.lt.0.6.and.abs(ctpi).lt.2.0) then
            prh(k,1) = prh(k,1)+1
            prh(16,1) = prh(16,1)+1
            if(shp.gt.0.05) then
             prh(k,2) = prh(k,2)+1
             prh(16,2) = prh(16,2)+1
            endif
           endif
          endif

          if(prp.gt.0.02 .and. abs(ctp).lt.1) then
           k = min(15,max(1,int(ceraero/1.)+1))
           aerohp(k) = aerohp(k)+1
           aerohp(16) = aerohp(16)+1
           if(j.lt.1) write(6,'("aero p",5f7.2)')
     >       ceraero,aerot,aeront
          endif
          if(prp.gt.0.02 .and. abs(ctpi).lt.1) then
           k = min(15,max(1,int(ceraero/1.)+1))
           aeroh(k) = aeroh(k)+1
           aeroh(16) = aeroh(16)+1
           if(j.gt.1000.and.j.lt.1000) 
     >       write(6,'("aero pi",5f7.2)')
     >       ceraero,aerot,aeront
           k = min(20,max(1,int((xaero+50.)/100.*20.)+1))
           xyaeroh(1,k) = xyaeroh(1,k)+1
           xyaeroh(1,21) = xyaeroh(1,21)+1
           if(ceraero.gt.2.5) then
            xyaeroh(2,k) = xyaeroh(2,k)+1
            xyaeroh(2,21) = xyaeroh(2,21)+1
           endif
           k = min(20,max(1,int((yaero+50.)/100.*20.)+1))
           xyaeroh(3,k) = xyaeroh(3,k)+1
           xyaeroh(3,21) = xyaeroh(3,21)+1
           if(ceraero.gt.2.5) then
            xyaeroh(4,k) = xyaeroh(4,k)+1
            xyaeroh(4,21) = xyaeroh(4,21)+1
           endif
          endif

! Coin time pion histos at SHMS hodo planes
c         ppi = abs(pp) * (1. + dpp/100.)
          if(ctpi.gt.-3.0 .and. ctpi.lt.1. .and.
     >      ceraero.gt.2.5 .and.
     >      (ppi.lt.3.0 .or. cerhg.gt.2.0)) then
           k = max(1,min(10,int((ctpi + 3.0) / 0.4 ) + 1)) 
           xpad = pdcx + 56. * pdcxp
           ypad = pdcy + 56. * pdcyp
           if(abs(xpad).lt.60. .and.abs(ypad).lt.40.) then
            ix = int((xpad+60.)/120.*100.)+1
            iy = int((ypad+40.)/80.*3.)+1
            ctx(1,ix,iy,k) = ctx(1,ix,iy,k)+1
            ix = int((xpad+60.)/120.*3.)+1
            iy = int((ypad+40.)/80.*100.)+1
            cty(1,ix,iy,k) = cty(1,ix,iy,k)+1
           endif
           xpad = pdcx + 266. * pdcxp
           ypad = pdcy + 266. * pdcyp
           if(abs(xpad).lt.60. .and.abs(ypad).lt.40.) then
            ix = int((xpad+60.)/120.*100.)+1
            iy = int((ypad+40.)/80.*3.)+1
            ctx(2,ix,iy,k) = ctx(2,ix,iy,k)+1
            ix = int((xpad+60.)/120.*3.)+1
            iy = int((ypad+40.)/80.*100.)+1
            cty(2,ix,iy,k) = cty(2,ix,iy,k)+1
           endif
           xpad = hdcx + 90. * hdcxp
           ypad = hdcy + 90. * hdcyp
           if(abs(xpad).lt.60. .and.abs(ypad).lt.40.) then
            ix = int((xpad+60.)/120.*100.)+1
            iy = int((ypad+40.)/80.*3.)+1
            ctx(3,ix,iy,k) = ctx(3,ix,iy,k)+1
            ix = int((xpad+60.)/120.*3.)+1
            iy = int((ypad+40.)/80.*100.)+1
            cty(3,ix,iy,k) = cty(3,ix,iy,k)+1
           endif
           xpad = hdcx + 310. * hdcxp
           ypad = hdcy + 310. * hdcyp
           if(abs(xpad).lt.60. .and.abs(ypad).lt.40.) then
            ix = int((xpad+60.)/120.*100.)+1
            iy = int((ypad+40.)/80.*3.)+1
            ctx(4,ix,iy,k) = ctx(4,ix,iy,k)+1
            ix = int((xpad+60.)/120.*3.)+1
            iy = int((ypad+40.)/80.*100.)+1
            cty(4,ix,iy,k) = cty(4,ix,iy,k)+1
           endif
          endif

c get SHMS and  hodo time relative to rf. Sometimes need to correct
c for alternate bucket.
          hfptcr = hfpt - hrf + 401.
c run 5038 is first one where rf timing works
          if(irun.ge.5038.and.irun.le.5046) hfptcr = hfptcr- rff/2.
          if(irun.ge.5183.and.irun.le.5185) hfptcr = hfptcr- rff/2.
          if(irun.ge.5244.and.irun.le.5377) hfptcr = hfptcr- rff/2.
          if(irun.ge.5602.and.irun.le.5662) hfptcr = hfptcr- rff/2.
          if(irun.ge.5840.and.irun.le.5941) hfptcr = hfptcr- rff/2.
          if(irun.ge.6049.and.irun.le.6068) hfptcr = hfptcr- rff/2.
          if(irun.ge.6080.and.irun.le.6129) hfptcr = hfptcr- rff/2.
          if(irun.ge.6219.and.irun.le.6267) hfptcr = hfptcr- rff/2.
          if(irun.ge.6312.and.irun.le.6511) hfptcr = hfptcr- rff/2.
          if(irun.ge.6518.and.irun.le.6600) hfptcr = hfptcr- rff/2.
          if(irun.ge.7702.and.irun.le.7736) hfptcr = hfptcr- rff/2.
          if(irun.ge.7978.and.irun.le.7989) hfptcr = hfptcr- rff/2.
c          if(irun.ge.7702.and.irun.le.7736) hfptcr = hfptcr- rff/2.
c          if(irun.ge.8050.and.irun.le.8180) hfptcr = hfptcr- rff/2.
          hfptc = hfptcr - rff*int(hfptcr/rff)
          hnrf = int(hfptcr/rff)
          pfptcr = pfpt - prf + 401.
c adjust for even/odd buckets
          if(irun.ge.5038.and.irun.le.5046) pfptcr = pfptcr- rff/2.
          if(irun.ge.5183.and.irun.le.5185) pfptcr = pfptcr- rff/2.
          if(irun.ge.5244.and.irun.le.5377) pfptcr = pfptcr- rff/2.
          if(irun.ge.5602.and.irun.le.5662) pfptcr = pfptcr- rff/2.
          if(irun.ge.5840.and.irun.le.5941) pfptcr = pfptcr- rff/2.
          if(irun.ge.6049.and.irun.le.6068) pfptcr = pfptcr- rff/2.
          if(irun.ge.6080.and.irun.le.6129) pfptcr = pfptcr- rff/2.
          if(irun.ge.6219.and.irun.le.6267) pfptcr = pfptcr- rff/2.
          if(irun.ge.6312.and.irun.le.6511) pfptcr = pfptcr- rff/2.
          if(irun.ge.6518.and.irun.le.6600) pfptcr = pfptcr- rff/2.
          if(irun.ge.7702.and.irun.le.7736) pfptcr = pfptcr- rff/2.
          if(irun.ge.7978.and.irun.le.7989) pfptcr = pfptcr- rff/2.
          if(irun.ge.8050.and.irun.le.8070) pfptcr = pfptcr- rff/2.
          if(irun.ge.8073.and.irun.le.8180) pfptcr = pfptcr- rff/2.

          if(irun.lt.5360) pfptcr = pfptcr + 0.10

          if(irun.gt.7590) pfptcr = pfptcr- 1.0

c adjustment for lower beam energies
c 6 gev
          if(irun.gt.7870.and.irun.lt.7978) pfptcr = pfptcr - 0.45
c 8gev
          if(irun.ge.7978.and.irun.lt.7990) pfptcr = pfptcr - 0.25
          if(irun.ge.7990.and.irun.lt.9999) pfptcr = pfptcr - 0.22

          pfptc = pfptcr - rff*int(pfptcr/rff)

          pnrf = int(pfptcr/rff)
c         if((j/100)*100.eq.j) write(6,'(''rf'', 3i5,6f8.2)') 
c     >    hnrf, pnrf, hnrf - pnrf -10,
c     >    ctpi,hrf,hfpt,hfptcr,hfptc,pfptc

! Don't Apply target position correction (ignoring cos(thp) term)
c (it is already taken into account by using rf timing!)
c         pfptc  = pfptc  - hztar / 30. 


! coin time pion in bins of p and four cerenkov combinations, for 
! runs with positive SHMS polairy
! Also, fp - rf time
          if(ctpi.gt.-2.0 .and. ctpi.lt.6. .and.
     >     shp .lt. 0.7 .and. ppi.gt. 1.8) then
           ipm=1
           if(pp.lt.0.) ipm=2
           if(ppi.lt. 3.8) then
            ip = int((ppi-1.8)*10)+1
            k = max(1,min(40,int((ctpi + 2.0) / 0.2 ) + 1)) 
            if(ceraero.gt.2.5 .and.cerhg.lt.1.5) 
     >       cthist(ip,1,k) = cthist(ip,1,k) + 1
            if(ceraero.gt.2.5 .and.cerhg.ge.1.5) 
     >       cthist(ip,2,k) = cthist(ip,2,k) + 1
            if(ceraero.le.2.5 .and.cerhg.lt.1.5) 
     >       cthist(ip,3,k) = cthist(ip,3,k) + 1
            if(ceraero.le.2.5 .and.cerhg.ge.1.5) 
     >       cthist(ip,4,k) = cthist(ip,4,k) + 1
            k = max(1,min(40,int((pfptc) / 0.1 ) + 1)) 
            if(ctpi.lt.2.and. irun.gt.5400) then
             if(ceraero.gt.2.5 .and.cerhg.lt.0.5) 
     >        ctrfhist(ip,2,k,ipm) = ctrfhist(ip,2,k,ipm) + 1
             if(ceraero.gt.2.5 .and.cerhg.ge.1.5) 
     >        ctrfhist(ip,3,k,ipm) = ctrfhist(ip,3,k,ipm) + 1
             if(ceraero.le.0.5 .and.cerhg.lt.0.5) 
     >        ctrfhist(ip,4,k,ipm) = ctrfhist(ip,4,k,ipm) + 1
             ctrfhist(ip,1,k,ipm) = ctrfhist(ip,1,k,ipm) + 1
            endif ! tighter ctpi cut
           endif ! check on pion momentum<3.8
! RF timing 
           if(ceraero.gt.2.5 .and. hfptcr.gt.0. .and.
     >      hfptcr.lt.1000.) then
            k = int(hfptcr*5.)
            ctrf(k,2) = ctrf(k,2)+1
            k = max(1,min(20,int((hfptc) / 0.1 ) + 1)) 
            ctrfrun(k,2) = ctrfrun(k,2) + 1
            ctrfrun(21,2) = ctrfrun(21,2) + 1
            do kk=1,10
             k = max(1,min(20,int((hfptcc(kk)) / 0.1 ) + 1)) 
             ctrfrn(kk,k,2) = ctrfrn(kk,k,2) + 1
            enddo
           endif
           if(ceraero.gt.2.5 .and. pfptcr.gt.0. .and.
     >      pfptcr.lt.1000.) then
            k = int(pfptcr*5.)
            ctrf(k,1) = ctrf(k,1)+1
            k = max(1,min(20,int((pfptc) / 0.2 ) + 1)) 
            ctrfrun(k,1) = ctrfrun(k,1) + 1
            ctrfrun(21,1) = ctrfrun(21,1) + 1
            do kk=1,10
             k = max(1,min(20,int((pfptcc(kk)) / 0.1 ) + 1)) 
             ctrfrn(kk,k,1) = ctrfrn(kk,k,1) + 1
            enddo
c           if((j/100)*100.eq.j) write(6,
c     >      '(4f8.3)') ppi,ctpi,hgtry,pgtry
           endif
          endif

! ytarg histo
          zh = hztar
          zp = pztar
          call newz(dpp,pdcy,pdcyp,ztnew)
          ztnew = ztnew / sin(thp/57.3)
          zd = zh - zp
          zdnew = zh - ztnew 

          if(abs(ctpi).lt.2. .and. ceraero.gt.2.5) then
           kk = min(16,max(1,int(zh+9.)))
           if(abs(dpe).lt.8.) then
            yth(1,kk) = yth(1,kk)+1
           else
            yth(2,kk) = yth(2,kk)+1
           endif

           kk = min(16,max(1,int(zp+9.)))
           if(dpp.gt.-10.and.dpp.lt.20.) then
            yth(3,kk) = yth(3,kk)+1
           else
            yth(4,kk) = yth(4,kk)+1
           endif

           kk = min(16,max(1,int(zd+9.)))
           if(dpp.gt.-10.and.dpp.lt.20.) then
            yth(5,kk) = yth(5,kk)+1
           else
            yth(6,kk) = yth(6,kk)+1
           endif
          endif

! all events
          if(abs(ctpi).gt.0. .and. ceraero.gt.2.5) then
           k=7
           if(dpp.gt. 10.and.dpp.lt. 20.) k=8 
           if(dpp.gt.  0.and.dpp.lt. 10.) k=9  
           if(dpp.gt. -5.and.dpp.lt.  0.) k=10
           if(dpp.gt.-10.and.dpp.lt. -5.) k=11
           if(dpp.gt.-15.and.dpp.lt.-10.) k=12
           if(dpp.lt.-15.) k=13
           kk = min(16,max(1,int(zp+9.)))
           yth(k ,kk) = yth( k,kk)+1
          endif

! coin time histos in 1 nsec bins
c from this study, looks like optimum ceraero cut is 2.5
c         write(6,'(/i8,8f8.3)') j,ctpi,ceraero,shp
c also look at kaons
         if(ctpi.gt.-18.0 .and. ctpi.lt.28. .and.
     >     ctk.gt.-18. .and. ctk .lt.28. .and.
     >     shp .gt. 0.02 .and.
     >     shp .lt. 0.7) then
           k = max(1,min(50,int((ctpi + 18.) / 1. ) + 1)) 
           ctpihw(k,1) = ctpihw(k,1)+1
           if(ceraero.gt.1.) ctpihw(k,2) = ctpihw(k,2)+1
           if(ceraero.gt.2.) ctpihw(k,3) = ctpihw(k,3)+1
           if(ceraero.gt.3.) ctpihw(k,4) = ctpihw(k,4)+1
           if(ceraero.gt.4.) ctpihw(k,5) = ctpihw(k,5)+1
           k = max(1,min(50,int((ctk + 18.) / 1. ) + 1)) 
           ctpihw(k,6) = ctpihw(k,6)+1
           if(ceraero.lt.3.)  ctpihw(k,7) = ctpihw(k,7)+1
           if(ceraero.lt.2.)  ctpihw(k,8) = ctpihw(k,8)+1
           if(ceraero.lt.1.)  ctpihw(k,9) = ctpihw(k,9)+1
           if(ceraero.lt.0.1) ctpihw(k,10) = ctpihw(k,10)+1
         endif


! Get flags for real or accidental pions, kaons
! reals
         ispireal = 0
         ispiacc = 0
         ispiaccL = 0
         ispiaccH = 0
         iskreal = 0
         iskacc = 0

c spring18 and 19 have narrow ctpi cuts
c use narrower cut for kaons
         if(irun.lt.4400 .or. irun.gt.7590) then
c this is too tight: cuts about 5% of events
c           tcut = 0.75
           tcut = 0.9
           tcutk = 0.6
           if(abs(ctpi - 0.0).lt.tcut) ispireal = 1
           if(abs(ctk - 0.0).lt.tcutk) iskreal = 1
! accidentals. Avoid ctpi <10 region, has protons
           if(abs(ctpi  -8.0).lt.tcut) ispiacc = 1
           if(abs(ctpi -12.0).lt.tcut) ispiacc = 1
           if(abs(ctpi  +4.0).lt.tcut) ispiacc = 1
           if(abs(ctpi  +8.0).lt.tcut) ispiacc = 1
           if(abs(ctpi  -8.0).lt.tcut) ispiaccL = 1
           if(abs(ctpi -12.0).lt.tcut) ispiaccL = 1
           if(abs(ctpi  +4.0).lt.tcut) ispiaccH = 1
           if(abs(ctpi  +8.0).lt.tcut) ispiaccH = 1
           if(abs(ctk   -8.0).lt.tcutk) iskacc = 1
           if(abs(ctk  -12.0).lt.tcutk) iskacc = 1
           if(abs(ctk   +4.0).lt.tcutk) iskacc = 1
           if(abs(ctk   +8.0).lt.tcutk) iskacc = 1
         else
c Fall18 runs
           tcut=2.6
           if(abs(ctpi - 0.0).lt.tcut) ispireal = 1
           if(abs(ctk - 0.0).lt.tcut) iskreal = 1
! accidentals
           if(abs(ctpi - 8.0).lt.tcut) ispiacc = 1
           if(abs(ctpi -16.0).lt.tcut) ispiacc = 1
           if(abs(ctpi + 8.0).lt.tcut) ispiacc = 1
           if(abs(ctpi +16.0).lt.tcut) ispiacc = 1
           if(abs(ctpi - 8.0).lt.tcut) ispiaccL = 1
           if(abs(ctpi -16.0).lt.tcut) ispiaccL = 1
           if(abs(ctpi + 8.0).lt.tcut) ispiaccH = 1
           if(abs(ctpi +16.0).lt.tcut) ispiaccH = 1
           if(abs(ctk  -24.0).lt.tcut) iskacc = 1
           if(abs(ctk  -16.0).lt.tcut) iskacc = 1
           if(abs(ctk  + 8.0).lt.tcut) iskacc = 1
           if(abs(ctk  +16.0).lt.tcut) iskacc = 1
          endif

c is this very likely an electron in SHMS?
          iselec = 1
          if(shp.lt.0.80) iselec = 0
          if(cerng .lt. 1 .and. irun.lt.4400) iselec = 0 
          if(cerhg .lt. 1) iselec = 0
          
! Get pion pt and phicm
         nu = e0 - ep
         q2 = 2. * e0 * ep * (1. - p_z(1)/ep)
         x = q2 / 2. / amp / nu
         epv(1)=p_x(1)
         epv(2)=p_y(1)
         epv(3)=p_z(1)
         epv(4)=ep
         p1vp(1)=p_x(2)
         p1vp(2)=p_y(2)
         p1vp(3)=p_z(2)

         p1vp(4)= sqrt(ppi**2 + ampi**2)
         zpi = p1vp(4) / nu
!xxx test
c        call getphi(e0,epv,p1vp,phicm)
c        write(6,'(/8f7.4)') 
c    >    p_x(1), p1vp(1), p_x(1) + p1vp(1),
c    >    p_y(1), p1vp(2), p_y(1) + p1vp(2),phicm
c
c        p1vp(2) = p_y(2) + ppi * 0.020
c        call getphi(e0,epv,p1vp,phicm)
c        write(6,'(8f7.4)') 
c    >    p_x(1), p1vp(1), p_x(1) + p1vp(1),
c    >    p_y(1), p1vp(2), p_y(1) + p1vp(2),phicm
c
c        p1vp(1) = p_x(2) + ppi * 0.020
c        call getphi(e0,epv,p1vp,phicm)
c        write(6,'(8f7.4)') 
c    >    p_x(1), p1vp(1), p_x(1) + p1vp(1),
c    >    p_y(1), p1vp(2), p_y(1) + p1vp(2),phicm

c        p1vp(1)=p_x(2)
c        p1vp(2)=p_y(2)
cxxx end test

         call getphi(e0,epv,p1vp,phicm)
         call getcos(e0,epv,p1vp,ampi,cthcm,pt)

         p1vp(4)= sqrt(ppi**2 + amk**2)
         zk = p1vp(4) / nu
         call getphi(e0,epv,p1vp,phicmk)
         call getcos(e0,epv,p1vp,amk,cthcmk,ptk)

c aerogel time hitograms
         if(ceraero.gt.0.0 .and. aerot.ne.0)  then
          icc=int((aerot+50.)/5)+1
          icc=min(15,max(1,icc))
          kk = int(ceraero)+1
          kk = min(20,max(1,kk))
          if(abs(ctp).lt.0.5) then
           aeroth(kk,icc) = aeroth(kk,icc)+1
           aeroth(kk,16) = aeroth(kk,16)+1
          endif
          if(abs(ctpi).lt.0.5) then
           aeroth(20+kk,icc) = aeroth(20+kk,icc)+1
           aeroth(20+kk,16) = aeroth(20+kk,16)+1
          endif
         endif

c cointime / rf in bins of z for 4 aerogel ranges
         if(ctpi.gt.-30.0 .and. ctpi.lt.30. 
     >    .and. (ppi .lt. 3.0 .or. cerhg .gt. 0.5)
     >    .and. shp.gt.0.02
     >     ) then
          iz = min(20,int(zpi*20.)+1)
          if(abs(ctpi).lt.2.0) then
           kk = int((ctpi+2.)/0.1)+1
           if(irun.gt.4400) then
            kk = min(40,int(pfptc*10.)+1)
           endif
           if(ceraero.le.0.) then
            cntsct(iz,kk,1) = cntsct(iz,kk,1) + 1 
           endif
           if(ceraero.gt.0.0 .and. ceraero.le.2.5) then
            cntsct(iz,kk,2) = cntsct(iz,kk,2) + 1 
           endif
           if(ceraero.gt.2.5 .and. ceraero.le.4.0) then
            cntsct(iz,kk,3) = cntsct(iz,kk,3) + 1 
           endif
           if(ceraero.gt.4.0) then
            cntsct(iz,kk,4) = cntsct(iz,kk,4) + 1 
           endif
          endif
         endif

c This is the main definition of good pions
c aerogel efficiency is about 99% for 2.5, and 97% for 3.5
c TODO change for aerotray=1??
         if(ctpi.gt.-30.0 .and. ctpi.lt.30. 
c put this below
     >     .and. ceraero.gt.2.5 
     >    .and. (ppi .lt. 3.0 .or. cerhg .gt. 0.5)
     >    .and. shp.gt.0.02
     >     ) then

c simple counters
          if(ispireal.eq.1) pir = pir + 1
          if(ispiacc.eq.1) piacc = piacc + 1
          if(ispiaccL.eq.1) piaccL = piaccL + 1
          if(ispiaccH.eq.1) piaccH = piaccH + 1
          im = max(1,min(20,int(mmpi2/0.5)+1))
          if(ispireal.eq.1) pirm(im) = pirm(im) + 1
          if(ispiacc.eq.1) piaccm(im) = piaccm(im) + 1
          im = min(20,int(zpi/0.05)+1)
          if(ispireal.eq.1) pirz(im) = pirz(im) + 1
          if(ispiacc.eq.1) piaccz(im) = piaccz(im) + 1

c counters versus kinematic
          k = max(1,min(20,int((dpe+15.)/30. * 20.)+1))
          if(ispireal.eq.1) pirk(1,k) = pirk(1,k) + 1
          if(ispiacc.eq.1) piacck(1,k) = piacck(1,k) + 1

          k = max(1,min(20,int((dthe*1000.+70.)/140.*20.)+1))
          if(ispireal.eq.1) pirk(2,k) = pirk(2,k) + 1
          if(ispiacc.eq.1) piacck(2,k) = piacck(2,k) + 1

          k = max(1,min(20,int((dphie*1000.+35.)/ 70.*20.)+1))
           if(ispireal.eq.1) pirk(3,k) = pirk(3,k) + 1
          if(ispiacc.eq.1) piacck(3,k) = piacck(3,k) + 1

          k = max(1,min(20,int((dpp+25.)/70. * 20.)+1))
          if(ispireal.eq.1) pirk(4,k) = pirk(4,k) + 1
          if(ispiacc.eq.1) piacck(4,k) = piacck(4,k) + 1

          k = max(1,min(20,int((dthp*1000.+70.)/140.*20.)+1))
          if(ispireal.eq.1) pirk(5,k) = pirk(5,k) + 1
          if(ispiacc.eq.1) piacck(5,k) = piacck(5,k) + 1

          k = max(1,min(20,int((dphip*1000.+ 35.)/ 70.*20.)+1))
          if(ispireal.eq.1) pirk(6,k) = pirk(6,k) + 1
          if(ispiacc.eq.1) piacck(6,k) = piacck(6,k) + 1
          
c big array
          pt2 = pt**2
          irr = 0
          if(ispireal.eq.1) irr=1
          if(ispiacc.eq.1) irr=2
          if(dphie.gt.0.) iq2bin=1
          if(dphie.le.0.) iq2bin=2 
c just using one q2bin now
c           iq2bin=1


          if(irr.gt.0 .and. zpi.lt.1.0
     >     .and.sqrt(pt2).lt.1.00) then
           iz = int(zpi*20.)+1
           izz = int(zpi*50.)+1
           iphi = max(1,min(15,int(phicm/6.29*15 ) + 1)) 
c          if(iz.eq.10) write(6,'(6f8.3)') ppi,nu,p1vp(4),zpi
c           ipt = int(pt2 / 0.5 * 12) + 1
           ipt = int(sqrt(pt2) / 1.00 * 16.) + 1
           im = max(1,min(50,int(mmpi2/0.2)+1))
           phih(iphi) = phih(iphi)+1
           phih(20) = phih(20)+1
c pions
            cnts(iq2bin,ipt,iphi,iz,irr) = 
     >        cnts(iq2bin,ipt,iphi,iz,irr) + 1
            sumcnt(ihel,irr) = sumcnt(ihel,irr) + 1
            sumcntsf(ihel,irr) = sumcntsf(ihel,irr) + 
     >        sin(phicm)
            cntsmmpi(im,irr) = cntsmmpi(im,irr) + 1
            cntsz(izz,irr) = cntsz(izz,irr) + 1
            if(cere.lt.-1000.) write(6,'(''cer'',4f8.3)')
     >       cere,dpe,dphie,dthe
c version with tighter cuts
           if((ppi .lt. 2.8 .or. cerhg .gt. 0.5).and.
     >        (ppi .lt. 3.0 .or. cerhg .gt. 1.5).and.
     >        ceraero.gt.5.0 .and.
     >        iselec .eq. 0 .and. 
     >        prp.gt.0.005.and.
     >        shp.gt.0.05) then
             cnts(iq2bin,ipt,iphi,iz,irr+2) = 
     >        cnts(iq2bin,ipt,iphi,iz,irr+2) + 1
             cntsmmpi(im,irr+2) = cntsmmpi(im,irr+2) + 1
             cntsz(izz,irr+2) = cntsz(izz,irr+2) + 1
c require good rf time also except when doesnt exist
c this assumes peak has been set a 1 nsec. 
             if((irun.lt.5040 .or. abs(pfptc-1.00).lt.0.5)) then
              cnts(iq2bin,ipt,iphi,iz,irr+4) = 
     >         cnts(iq2bin,ipt,iphi,iz,irr+4) + 1
              cntsmmpi(im,irr+4) = cntsmmpi(im,irr+4) + 1
              cntsz(izz,irr+4) = cntsz(izz,irr+4) + 1
             endif
            else
c             write(6,'(i6,3i3)') iselec,pexit,hexit
            endif
          endif
c exclusive pion counts
          if(irr.gt.0 .and. mmpi2.gt.0.76 .and.mmpi2.lt.1.0 .and.
     >      iselec .eq. 0 .and. 
     >      (ppi .lt. 3.0 .or. cerhg .gt. 0.5)
     >     .and.sqrt(pt2).lt.0.3) then
           iphi = max(1,min(15,int(phicm/6.29*7 ) + 1)) 
           ipt = int(sqrt(pt2) / 0.3 * 4.) + 1
           im = max(1,min(8,int((mmpi2-0.76)/0.24*8.)+1))
           cntse(1,ipt,iphi,irr) = 
     >       cntse(1,ipt,iphi,irr) + 1 
           cntseh(1,ipt,iphi,im,irr) = 
     >        cntseh(1,ipt,iphi,im,irr) + 1
          endif

c histogram of cointime
          if(ctpi.gt.-10.0 .and. ctpi.lt.30.0 .and.
     >      ceraero.gt.2.5) then
            k = max(1,min(40,int((ctpi + 10.) / 1. ) + 1)) 
            ctpih(k) = ctpih(k)+1
          endif

         endif ! end main definition good pions

c This is the main definition of good kaons
         if(ctk.gt.-30.0 .and. ctk .lt.30. .and.
     >     ceraero.lt.3.5 .and.
c changed: only keep kaons for p<2.8 (also done for SIMC)
c took out
c     >     ppi.lt.2.8 .and.
     >     iselec .eq. 0 .and. 
     >     cerhg .lt. 0.5 .and.
     >     shp.gt.0.02) then

c big array
         pt2 = ptk**2
         irr = 0
         if(iskreal.eq.1) irr=7
         if(iskacc.eq.1) irr=8
         if(dphie.gt.0.) iq2bin=1
         if(dphie.le.0.) iq2bin=2
c just one bin
c         iq2bin=1

         if(irr.gt.0 .and. zk.lt.1.0 .and.sqrt(pt2).lt.1.00) then
           Emk = e0 + amp - ep - sqrt(amk**2 + ppi**2)
           mmk2 = 
     >      (Emk)**2 - 
     >      (p_x(1) + p_x(2))**2 -
     >      (p_y(1) + p_y(2))**2 -
     >      (p_z(1) + p_z(2)- e0)**2
           im = max(1,min(50,int(mmk2/0.2)+1))
          iz = int(zk*20.)+1
          izz = int(zk*50.)+1
          iphi = max(1,min(15,int(phicmk/6.29*15 ) + 1)) 
c          ipt = int(pt2 / 0.5 * 12) + 1
          ipt = int(sqrt(pt2) / 1.00 * 16.) + 1
c kaons in irr=7,8
            cnts(iq2bin,ipt,iphi,iz,irr) = 
     >       cnts(iq2bin,ipt,iphi,iz,irr) + 1
            cntsmmpi(im,irr) = cntsmmpi(im,irr) + 1
            cntsz(izz,irr) = cntsz(izz,irr) + 1
          endif

         endif ! end main definition good kaons

         if(ctpi.gt.-2.0 .and. ctpi.lt.6. .and.
     >      shp .lt. 0.7) then
            k = max(1,min(40,int((ctpi + 2.0) / 0.2 ) + 1)) 
            ctpihfa(k) = ctpihfa(k)+1
            if(ceraero.gt.2.5) then
             ctpihf(k) = ctpihf(k)+1
            endif
          endif
         if(ctk.gt.-2.0 .and. ctk.lt.6. .and.
     >     shp .lt. 0.7.and.ceraero.eq.0.) then
           k = max(1,min(40,int((ctk + 2.0) / 0.2 ) + 1)) 
           ctkhf(k) = ctkhf(k)+1
         endif
         if(ctp.gt.-2.0 .and. ctp.lt.6. .and.
     >     shp .lt. 0.7.and.ceraero.eq.0.) then
           k = max(1,min(40,int((ctp + 2.0) / 0.2 ) + 1)) 
           ctphf(k) = ctphf(k)+1
         endif
         if(ctk.gt.-10.and.ctk.lt.30.0 .and. 
     >     ceraero.lt.0.5 .and.
     >     shp .lt. 0.7.and.abs(pp).lt.3.0) then
           k = max(1,min(40,int((ctk + 9.5) / 1. ) + 1)) 
           ctkh(k) = ctkh(k)+1
         endif
c         if(ctk.gt.-2.and.ctk.lt.6.0 .and. 
c     >     ceraero.lt.0.5 .and.
c     >     shp .lt. 0.7.and.abs(pp).lt.3.0) then
c           k = max(1,min(40,int((ctk + 2.0) / 0.2 ) + 1)) 
c           ctkhf(k) = ctkhf(k)+1
c         endif
         if(ctpi.gt.-10.0 .and. ctpi.lt.30.0 .and.
     >     ceraero.gt.2.5 .and.
     >     shp .gt. 0.8.and.cerhg.gt.1.0.and.
     >     (cerng.gt.1.0.or.irun.gt.4400)) then
           k = max(1,min(40,int((ctpi + 9.5) / 1. ) + 1)) 
           cteh(k) = cteh(k)+1
         endif
        endif ! cuts on dpp, dpe


        if(ceraero.gt.2.5 .and. shp .lt. 0.7 .and.
     >     dpe.gt.-12.0 .and. dpe.lt.12. .and.
     >     she.gt.0.7 .and. cere.gt.0.5 .and.
     >     mmpi2.gt.0.85**2 .and. mmpi2.lt.1.05**2 .and.
     >     dpp.gt.-20.0 .and. dpp.lt.50.) then
         k = max(1,min(20,int((sqrt(mmpi2)-0.85)/0.010 ) + 1)) 
         kk = max(1,min(12,int((dpp + 20.)/5.) + 1))
         if(abs(ctpi).lt.2.0) then
           mmpih(k,kk) = mmpih(k,kk)+1
         endif
         if(abs(ctpi).lt.2.0 .and. it.eq.1 .and. pp.gt.0.0
     >    .and. dpp.gt.-12.0 .and. dpp.lt. 22) 
     >      mmpihs(k) = mmpihs(k)+1
         if(abs(ctpi).gt.2.0 .and. it.eq.1 .and. pp.gt.0.0
     >    .and. dpp.gt.-12.0 .and. dpp.lt. 22) 
     >      mmpihsa(k) = mmpihsa(k)+1
        endif

        mmk2 = 
     >     (Emk)**2 - 
     >     (p_x(1) + p_x(2))**2 -
     >     (p_y(1) + p_y(2))**2 -
     >     (p_z(1) + p_z(2)- e0)**2

        mmp2 = 
     >     (Emp)**2 - 
     >     (p_x(1) + p_x(2))**2 -
     >     (p_y(1) + p_y(2))**2 -
     >     (p_z(1) + p_z(2)- e0)**2
        if(ceraero.lt.0.5 .and. shp .lt. 0.7 .and.
     >     cerhg .lt. 0.5 .and.
     >     dpe.gt.-12.0 .and. dpe.lt.12. .and.
     >     she.gt.0.7 .and. cere.gt.1.5 .and.
c     >     mmp2.gt.0.7**2.and.mmp2.lt.1.1**2 .and.
     >     mmp2.gt.0.0**2.and.
     >     dpp.gt.-25.0 .and. dpp.lt.45.) then
c           k = max(1,min(20,int((sqrt(mmp2)-0.7)/0.02 ) + 1)) 
         k = max(1,min(20,int((sqrt(mmp2)-0.0)/0.1 ) + 1)) 
         kk = max(1,min(14,int((dpp+25.)/5. ) + 1)) 
         if(abs(ctp).lt.2.0) mmph(k,kk) = mmph(k,kk)+1
         if(abs(ctp).lt.2.0 .and. it.eq.1 .and. pp.gt.0.0
     >    .and.dpp.gt.-12.0 .and. dpp.lt. 22.) 
     >    mmphs(k) = mmphs(k)+1
         k = max(1,min(400,int((sqrt(mmp2))/0.005 ) + 1)) 
         if(abs(ctp).lt.2.0) mmphf(k) = mmphf(k)+1
         if(abs(ctp).lt.2.0 .and. it.eq.1.and. pp.gt.0.0
     >    .and. dpp.gt.-12.0 .and. dpp.lt. 22.)
     >    mmphfs(k) = mmphfs(k)+1
c         write(6,'(''w'',3f8.3)') w,mmp2,ctp
         if(abs(w-1.5).lt.0.1 .and. it.eq.1 .and. pp.gt. 0.) then
          if(abs(ctp).lt.2.0) mmpheta(k) = mmpheta(k)+1
         endif
        endif

! e,e mass
        mee = sqrt(
     >     (ep + ppi)**2 - 
     >     (p_x(1) + p_x(2))**2 -
     >     (p_y(1) + p_y(2))**2 -
     >     (p_z(1) + p_z(2))**2)
        if(ceraero.gt.2.5 .and. shp .gt. 0.85 .and.
     >     cerhg .gt. 4.0 .and.
     >     (cerng .gt. 4.0 .or. irun.gt. 4400) .and.
     >     she.gt.0.7 .and. cere.gt.1.5 .and.
     >     mee.gt.3.00.and.mee.lt.3.20 .and.
     >     abs(ctpi).lt.1.0) then
c         write(6,'(''mee'',4f8.3)') mee,ceraero,shp,ctpi
         write(49,'(9f7.3)') mee,shp,p_x(1),p_y(1),p_z(1),
     >    p_x(2),p_y(2),p_z(2),pp
         k = max(1,min(20,int((mee-3.00)/0.010 ) + 1)) 
         meeh(k) = meeh(k)+1
         if(irun.lt.9990) then
          if(pp.lt.0.) then
           meeht(k,1) = meeht(k,1)+1
          else
           meeht(k,2) = meeht(k,2)+1
          endif
         endif
        endif

! particle ID plots
 
! noble gas
        if(ceraero.gt.2.5 .and. shp .gt. 0.7 .and.
     >     cerhg.gt.1.0) then
         kk = max(1,min(15,int(cerng/2.)+1))
         ngh(kk) = ngh(kk) + 1
         ngh(16) = ngh(16) + 1
        endif

! hg spectra versus delta. 
        if(ceraero.gt.2.5 .and. shp .lt. 0.8 .and.
     >     dpe.gt.-10.0 .and. dpe.lt.10. .and.
     >     abs(ctpi).lt.1.0 .and.
     >     dpp.gt.-5.0 .and. dpp.lt.5.0 .and.
     >     ppi.gt.3.3) then
           k =  int((dpp + 5.)*3.)+1
           kk = min(50,int(cerhg)+1)
         if(k.gt.0.and.k.lt.31.and.kk.gt.0.and.kk.lt.51) then
          hgh(k,kk) = hgh(k,kk) + 1
         endif
        endif
! hg spectra versus ppi 
        if(ceraero.gt.2.5 .and. shp .lt. 0.8 .and.
     >     dpe.gt.-10.0 .and. dpe.lt.10. .and.
     >     abs(ctpi).lt.1.0 .and.
c     >     dpp.gt.-5.0 .and. dpp.lt.5.0 .and.
     >     ppi.gt.2.0) then
           k =  int((ppi - 2.0)/0.1)+1
           kk = min(50,int(cerhg * 2.0)+1)
         if(k.gt.0.and.k.lt.31.and.kk.gt.0.and.kk.lt.51) then
          hghp(k,kk) = hghp(k,kk) + 1
         endif
        endif

        enddo ! LOOP OVER events
 10     continue

! Read in Hodo and track efficiencies
! these override anything done above
        write(fname,
     >    '(''/work/hallc/c-sidis18/bosted/Skimfiles/Skimeff'',
     >     i4,''.txt'')') irun
        open(unit=9,file=fname)
        read(9,*) i1,i2,heffh
        read(9,*) i1,i2,heffp
c took this out
c        corr = heffh * heffp

c corr now takes into account low trig eff for sp18
        corr = 1.0
        if(irun.lt.4400) corr = 0.95

c i4/i5 is fraction of events passing ptdcmult>4 cut
c implemented 4/30/2020 in PeterB.C
        read(9,*) i1,i2,sum1,i3,i4,sum2,i4,i5
        corr3 = float(i4)/float(i5)

        read(9,*) i1,i2,i3,i4,i5
        tefe = float(i5) / float(i2)
        read(9,*) i1,i2,i3,i4,i5
        tefp = float(i5) / float(i2)
        read(9,*) i1,i2,hgdsceff
        if(i2.eq.0) hgdsceff=1.0
        read(9,*) i1,i2,pgdsceff
c these are really just 3/4 or 4/4 eff (depending on
c ptracking parameter at bottom. 
c        corr = hgdsceff * pgdsceff

c pass over five lines
        read(9,*) sum1
        sum1 = 0.
        sum2 = 0.
        sum3 = 0.
        sum4 = 0.
c get new way of doing tefp and timecorr
        do kk=1,100
         read(9,*) i1,i2,i3,r1,r2,r3,r4
         if(kk.gt.20.and.kk.lt.80.and.r4.lt.0.1.and.
     >    r4.gt.0.) then
          sum1 = sum1 + r3/r4**2
          sum2 = sum2 + 1./r4**2
c          write(6,'(3i5,4f7.3)') i1,i2,i3,r1,r2,r3,r4
         endif
         if(kk.ge.60.and.kk.lt.76) sum3 = sum3 + i1
         if(kk.ge.76.and.kk.lt.92) sum4 = sum4 + i1
        enddo
        if(sum2.gt.0.) then
c this doesn't change anything much
c         write(6,'(''new trk eff p'',2f7.3)') tefp,sum1/sum2
c         tefp = sum1 / sum2
c this doesn't work
c         write(6,'(''new timecorr '',2f7.3)') timecorr,sum3/sum4
c         timecorr = sum3 / sum4
        endif

! Read in SIMC file. Normfac is in simc/outfiles/*.his
! SIDIS pions
        srate(i)=0.
        srateer(i)=0.
        write(fname,
     >    '(''/group/c-sidis/bosted/simc/outfiles/simc_'',
     >     i4,''_pi_rad.hist'')') irun
        open(unit=9,file=fname)
        do jj=1,110
         read(9,'(a)',err=14,end=14) sstring
        enddo
c        write(6,'(a)') sstring
        read(sstring,'(30x,e15.4)') normfac
c Correction to make cross section per nucleon
c        write(6,'(''normfac='',e12.4)') normfac
c aluminum factor 27 times attenuation factor
        if(it.eq.3) normfac = normfac * 27. * 0.8
        write(fname,
     >    '(''/work/hallc/c-sidis18/bosted/simctxt/simc_'',
     >     i4,''_pi_rad'')') irun
        open(unit=9,file=fname)
        sum1 = 0.
        sum2 = 0.
        sum3 = 0.
        do im=1,20
         sum1m(im)=0.
         sum2m(im)=0.
         sratem(im)=0.
         sratemer(im)=0.
         sum1z(im)=0.
         sum2z(im)=0.
         sratez(im)=0.
         sratezer(im)=0.
        enddo

! order written is yptar, then xptar
 188	format(e12.4,21f9.4,2e12.4,2f8.3)
        do jj=1,simcevents
         read(9,188,end=13) weight,
     >    dpe,dphie,dthe,hztar,hdcx,hdcxp,hdcy,hdcyp,
     >    dpp,dphip,dthp,pztar,pdcx,pdcxp,pdcy,pdcyp,
     >    dpp_init, xptar_init, yptar_init,yt_orig,yrast,
     >    sigcc,sigcm,decdist,m2final

         if(jj.lt.0  ) then
          write(6,'(''simc 1'',2f8.3)') decdist,m2final
c          write(6,'(''simc chk'',8f7.2)') 
c     >     dpp, dpp_init, dthp*10., xptar_init,
c     >     dphip*10., yptar_init, hztar, 
c     >     yt_orig/sin(thp * pi/180.)
         endif
         sum3 = sum3 + 1.
         dthe = dthe / 100.
c need to change sign!
c         dphie = -1. * dphie
         dphie = dphie / 100.
         dthp = dthp / 100.
         dphip = dphip / 100.
         pdcxp = pdcxp / 100.
         pdcyp = pdcyp / 100.
         hdcxp = hdcxp / 100.
         hdcyp = hdcyp / 100.
         xaero = pdcx + 231. * pdcxp
         yaero = pdcy + 231. * pdcyp

c fill in acceptance array
         okhms = 0
         okshms = 0
         if(abs(dpe).lt.13 .and.
     >    abs(dthe).lt.0.100 .and. 
     >    abs(dphie).lt.0.030) then
          ith = int((dthe + 0.100) / 0.200 * 70)+1
          iphi = int((dphie + 0.030) / 0.060 * 30.)+1
          ipt = int((dpe + 13.) / 26. * 30.)+1
          accep(2,ipt,ith,iphi) = accep(2,ipt,ith,iphi) + 1
          okhms = okacc(1,ipt,ith,iphi)
          if(okhms.eq.1) accep(6,ipt,ith,iphi) = 
     >     accep(6,ipt,ith,iphi) + 1
         endif
         if(abs(dpp-5.).lt.17 .and.
     >    abs(dthp).lt.0.100 .and. 
     >    abs(dphip).lt.0.030) then
          ith = int((dthp + 0.100) / 0.200 * 70)+1
          iphi = int((dphip + 0.030) / 0.060 * 30.)+1
          ipt = int((dpp + 12.) / 34. * 30.)+1
          accep(4,ipt,ith,iphi) = accep(4,ipt,ith,iphi) + 1
          okshms = okacc(2,ipt,ith,iphi)
          if(okshms.eq.1) accep(8,ipt,ith,iphi) = 
     >     accep(8,ipt,ith,iphi) + 1
         endif
c use same th and phi range for dpp<-12 as for -12
         if(dpp.lt.-12 .and.
     >    abs(dthp).lt.0.100 .and. 
     >    abs(dphip).lt.0.030) then
          ith = int((dthp + 0.100) / 0.200 * 70)+1
          iphi = int((dphip + 0.030) / 0.060 * 30.)+1
          ipt = 1
          okshms = okacc(2,ipt,ith,iphi)
         endif

         pexit=1
         xexit = pdcx -307. * pdcxp
         yexit = pdcy -307. * pdcyp
         crad = 23.81
         voffset = crad - 24.035
         hwid = 11.549/2.
         if(abs(yexit) .lt. hwid) then
           if(abs(xexit) .gt.  (crad + voffset)) pexit=0 
         else
           if ( yexit .ge. hwid) then
            if ( (xexit - voffset)**2 + 
     >           (yexit - hwid)**2 .gt.
     >            crad**2) pexit=0;
           endif
           if ( yexit .le. -1.*hwid) then
            if ( (xexit - voffset)**2 + 
     >           (yexit + hwid)**2 .gt.
     >            crad**2) pexit=0;
           endif
         endif
c         if(pexit.eq.0) 
c     >    write(6,'(''pexit'',i7,i2,6f8.3)') jj,
c     >      pexit,xexit,yexit,pdcx,pdcxp,pdcy,pdcyp

c HMS dipole exit
         xexit = hdcx - 148. * hdcxp
         yexit = hdcy - 148. * hdcyp
         hexit=1
         if ( (xexit - 2.8)**2 + 
     >        (yexit)**2 .gt. 46.607**2) hexit=0;
c         if(hexit.eq.0) 
c     >    write(6,'(''hexit'',i7,i2,6f8.3)') jj,
c     >      pexit,xexit,yexit,pdcx,pdcxp,pdcy,pdcyp

         if(dphie.gt.0.) iq2bin=1
         if(dphie.le.0.) iq2bin=2
c just one bin now
c        iq2bin=1

         the = the0 - dthe
         ep = ep0 * (1. + dpe/100.)
c         p_x(1) =  ep * sin(the) * cos(dphie)
c         p_y(1) = -ep * sin(the) * sin(dphie)
c         p_z(1) =  ep * cos(the) 
         ppi = abs(pp) * (1. + dpp/100.)
         thpi = thp * 3.1415928/180. + dthp
c         p_x(2) = -ppi * sin(thpi) * cos(dphip)
c         p_y(2) = -ppi * sin(thpi) * sin(dphip)
c         p_z(2) =  ppi * cos(thpi) 
	 call physics_angles(the0, phi0e, dthe, dphie, 
     >     ep,p_x(1),p_y(1),p_z(1))
	 call physics_angles(thp_rad, phi0p, dthp, dphip,
     >     ppi,p_x(2),p_y(2),p_z(2))
         nu = e0 - ep
         epv(1)=p_x(1)
         epv(2)=p_y(1)
         epv(3)=p_z(1)
         epv(4)=ep
         p1vp(1)=p_x(2)
         p1vp(2)=p_y(2)
         p1vp(3)=p_z(2)
         p1vp(4)= sqrt(ppi**2 + ampi**2)
         zpi = p1vp(4) / nu
         call getphi(e0,epv,p1vp,phicm)
         call getcos(e0,epv,p1vp,ampi,cthcm,pt)

         Empi = e0 + amp - ep - sqrt(ampi**2 + ppi**2)
         Emk = e0 + amp - ep - sqrt(amk**2 + ppi**2)
         Emp = e0 + amp - ep - sqrt(amp**2 + ppi**2)

         w2=(e0 + amp - ep)**2 - p_x(1)**2 - p_y(1)**2 - (p_z(1)-e0)**2
         w=0.
         if(w2.gt.0.) w = sqrt(w2)

         mmpi2 = 
     >     (Empi)**2 - 
     >     (p_x(1) + p_x(2))**2 -
     >     (p_y(1) + p_y(2))**2 -
     >     (p_z(1) + p_z(2)- e0)**2

         if(dpe.gt.dphmslo .and. dpe.lt.dphmshi .and.
     >      abs(dthe).lt.dthhms .and.
     >      abs(dthp).lt.dthshms .and.
     >      abs(dphie).lt.dphihms .and.
     >      abs(dphip).lt.dphishms .and.
c     >       okhms.eq.1 .and. okshms.eq.1 .and.
     >       pexit.eq.1 .and. hexit.eq.1.and. 
c add cuts at aerogel xy
     >      abs(xaero).lt.xmaxaero .and. abs(yaero).lt.ymaxaero.and.
     >      dpp.gt.dpshmslo .and. dpp.lt.dpshmshi) then
          sum1 = sum1 + weight
          sum2 = sum2 + 1.
          im = min(20,int(mmpi2/0.5)+1)
          sum1m(im) = sum1m(im) + weight
          sum2m(im) = sum2m(im) + 1
          im = min(20,int(zpi/0.05)+1)
          sum1z(im) = sum1z(im) + weight
          sum2z(im) = sum2z(im) + 1
c final state counter
          ipart = 0
          if(abs(m2final-m_mu).lt.0.01) ipart=1
          if(abs(m2final-m_pi).lt.0.01) ipart=2
          if(abs(m2final-m_k).lt.0.05) ipart=3
          idist=min(22,max(1,int((decdist/0.1)+1)))
          if(ipart.ne.0) fstate(1,idist,ipart) = 
     >      fstate(1,idist,ipart) +weight
          fstate(1,idist,4) = fstate(1,idist,4) +
     >     weight
          if(jj.lt.0 .and.ipart.ne.2) then
           write(6,'("f 1",2i3,f6.2,2e12.2)') ipart,
     >      idist,m2final,fstate(1,idist,ipart),
     >      fstate(1,idist,4)
          endif
c counters versus kinematic
          k = max(1,min(20,int((dpe+15.)/30. * 20.)+1))
          pirks(1,k) = pirks(1,k) + 1
          piaccks(1,k) = piaccks(1,k) + weight

          k = max(1,min(20,int((dthe*1000.+70.)/140.*20.)+1))
          pirks(2,k) = pirks(2,k) + 1
          piaccks(2,k) = piaccks(2,k) + weight

          k = max(1,min(20,int((dphie*1000.+ 35.)/ 70.*20.)+1))
          pirks(3,k) = pirks(3,k) + 1
          piaccks(3,k) = piaccks(3,k) + weight

          k = max(1,min(20,int((dpp+25.)/70. * 20.)+1))
          pirks(4,k) = pirks(4,k) + 1
          piaccks(4,k) = piaccks(4,k) + weight

          k = max(1,min(20,int((dthp*1000.+70.)/140.*20.)+1))
          pirks(5,k) = pirks(5,k) + 1
          piaccks(5,k) = piaccks(5,k) + weight

          k = max(1,min(20,int((dphip*1000.+ 35.)/ 70.*20.)+1))
          pirks(6,k) = pirks(6,k) + 1
          piaccks(6,k) = piaccks(6,k) + weight

          pt2 = pt**2
          if(zpi.lt.1.0 .and. sqrt(pt2) .lt. 1.00) then
           iz = int(zpi*20.)+1
           izz = int(zpi*50.)+1
           im = max(1,min(50,int(mmpi2/0.2)+1))
           iphi = max(1,min(15,int(phicm/6.29*15 ) + 1)) 
c           ipt = int(pt2 / 0.5 * 12) + 1
           ipt = int(sqrt(pt2) / 1.00 * 16.) + 1
           cntsmc(iq2bin,ipt,iphi,iz,1) = 
     >       cntsmc(iq2bin,ipt,iphi,iz,1) + 1
           cntsmc(iq2bin,ipt,iphi,iz,2) = 
     >       cntsmc(iq2bin,ipt,iphi,iz,2) + weight
           cntsmcmmpi(im,1) = cntsmcmmpi(im,1) + 1
           cntsmcmmpi(im,2) = cntsmcmmpi(im,2) + weight
           cntsmcz(izz,1) = cntsmcz(izz,1) + 1
           cntsmcz(izz,2) = cntsmcz(izz,2) + weight
c  with cut on mmpi2
c          if(mmpi2.gt.2.0) then
c           cntsmc(iq2bin,ipt,iphi,iz,5) = 
c    >       cntsmc(iq2bin,ipt,iphi,iz,5) + 1
c           cntsmc(iq2bin,ipt,iphi,iz,6) = 
c    >       cntsmc(iq2bin,ipt,iphi,iz,6) + weight
c          endif
             avkin(iq2bin,ipt,iphi,iz,1) = 
     >        avkin(iq2bin,ipt,iphi,iz,1) + dpe 
             avkin(iq2bin,ipt,iphi,iz,2) = 
     >        avkin(iq2bin,ipt,iphi,iz,2) + dthe*1000.
             avkin(iq2bin,ipt,iphi,iz,3) = 
     >        avkin(iq2bin,ipt,iphi,iz,3) + dphie*1000.
c changed from hztar to z
             avkin(iq2bin,ipt,iphi,iz,4) = 
     >        avkin(iq2bin,ipt,iphi,iz,4) + zpi
             avkin(iq2bin,ipt,iphi,iz,5 ) = 
     >        avkin(iq2bin,ipt,iphi,iz,5) + dpp
             avkin(iq2bin,ipt,iphi,iz,6) = 
     >        avkin(iq2bin,ipt,iphi,iz,6) + dthp*1000.
             avkin(iq2bin,ipt,iphi,iz,7) = 
     >        avkin(iq2bin,ipt,iphi,iz,7) + dphip*1000.
             avkin(iq2bin,ipt,iphi,iz,8) = 
     >        avkin(iq2bin,ipt,iphi,iz,8) + pt2 
          endif
         endif
        enddo
 13     srate(i) = sum1 * normfac / sum3
        srateer(i) = srate(i) / sqrt(sum2)
        do im=1,20
         sratem(im) = sum1m(im) * normfac / sum3
         sratemer(im) = sratem(im) / sqrt(max(1.,sum2m(im)))
         sratez(im) = sum1z(im) * normfac / sum3
         sratezer(im) = sratez(im) / sqrt(max(1.,sum2z(im)))
        enddo
        write(6,'(''normfac='',5e10.3)') 
     >    normfac,srate(i),sum1,sum2,sum3
        avnfac(1) = avnfac(1) + normfac * chrg
        sum3tot(1) = sum3tot(1) + sum3
        do k=1,6
         do kk=1,20
          ratesk(k,kk) = piaccks(k,kk)  * normfac / sum3
          ratesker(k,kk) = ratesk(k,kk)/sqrt(max(1.,pirks(k,kk)))
         enddo
        enddo

 14     continue

! Exclusive  pions, endcaps, and rho, no_rad, K, Knorad
        do icase=1,6
        sratex(i,icase)=0.
        sratexer(i,icase)=0.
        normfacx = 0.
        sum1 = 0.
        sum2 = 0.
        sum3 = 0.
        sum4 = 0.
        do im=1,20
         sum1m(im)=0.
         sum2m(im)=0.
         sum1z(im)=0.
         sum2z(im)=0.
         sratemx(im,icase)=0.
         sratemxer(im,icase)=0.
         sratezx(im,icase)=0.
         sratezxer(im,icase)=0.
        enddo
        if(icase.eq.1) write(fname,
     >    '(''/group/c-sidis/bosted/simc/outfiles/simc_'',
     >     i4,''_pi_excl.hist'')') irun
        if(icase.eq.2) write(fname,
     >    '(''/group/c-sidis/bosted/simc/outfiles/simc_'',
     >     i4,''_pi_rho.hist'')') irun
        if(icase.eq.3) write(fname,
     >    '(''/group/c-sidis/bosted/simc/outfiles/simc_'',
     >     i4,''_endcap.hist'')') irun
        if(icase.eq.4) write(fname,
     >    '(''/group/c-sidis/bosted/simc/outfiles/simc_'',
     >     i4,''_pi_norad.hist'')') irun
        if(icase.eq.5) write(fname,
     >    '(''/group/c-sidis/bosted/simc/outfiles/simc_'',
     >     i4,''_k_rad.hist'')') irun
        if(icase.eq.6) write(fname,
     >    '(''/group/c-sidis/bosted/simc/outfiles/simc_'',
     >     i4,''_k_norad.hist'')') irun
        open(unit=9,file=fname)
! sometimes on line 111, sometimes on 110
        do jj=1,112
         read(9,'(a)',err=16,end=16) sstring
c         write(6,'(i4,a)') jj,sstring(20:26)
         if(sstring(20:26).eq."normfac") then
c           write(6,'(a)') sstring
          read(sstring(30:60),*) normfacx
         endif
        enddo
        if(normfacx.eq.0.) write(6,
     >    '(''error normfac_x='',e12.4)') normfacx
c Correction to make cross section per nucleon
c not needed for exclusie 
c if(it.eq.3) normfacx = normfacx * 27.
c normalize endcaps down from Dummy, multiply by 27
c to get per nucleon for SIDIS, and put 0.8 for attentuation
        if(icase.eq.3) then
         if(it.eq.1) normfacx = normfacx * 0.262 * 27. * 0.8
         if(it.eq.2) normfacx = normfacx * 0.260 * 27. * 0.8
         if(it.eq.3) normfacx = normfacx * 0.000
        endif
c make rho twice as big for deuteron, dummy because
c only calculates for protons
        if(icase.eq.2) then
         if(it.eq.2) normfacx = normfacx * 2.
         if(it.eq.3) normfacx = normfacx * 2.
        endif
        if(icase.eq.1) write(fname,
     >    '(''/work/hallc/c-sidis18/bosted/simctxt/simc_'',
     >     i4,''_pi_excl'')') irun
        if(icase.eq.2) write(fname,
     >    '(''/work/hallc/c-sidis18/bosted/simctxt/simc_'',
     >     i4,''_pi_rho'')') irun
        if(icase.eq.3) write(fname,
     >    '(''/work/hallc/c-sidis18/bosted/simctxt/simc_'',
     >     i4,''_endcap'')') irun
        if(icase.eq.4) write(fname,
     >    '(''/work/hallc/c-sidis18/bosted/simctxt/simc_'',
     >     i4,''_pi_norad'')') irun
        if(icase.eq.5) write(fname,
     >    '(''/work/hallc/c-sidis18/bosted/simctxt/simc_'',
     >     i4,''_k_rad'')') irun
        if(icase.eq.6) write(fname,
     >    '(''/work/hallc/c-sidis18/bosted/simctxt/simc_'',
     >     i4,''_k_norad'')') irun
        open(unit=9,file=fname)

        do jj=1,simcevents
         read(9,'(a)',end=15,err=999) string
c fix stars
         do k=13, 13+189, 9
          if(string(k:k+8).eq.'*********') then
           write(string(k:k+8),'("  999.999")')
          endif
         enddo

c         if(icase.eq.1) then
          read(string,'(e12.4,21f9.4,2e12.4,2f8.3)',end=15,err=999) 
     >     weight,
     >     dpe,dphie,dthe,hztar,hdcx,hdcxp,hdcy,hdcyp,
     >     dpp,dphip,dthp,pztar,pdcx,pdcxp,pdcy,pdcyp,
     >     dpp_init, xptar_init, yptar_init,yt_orig,yrast
     >    ,sigcc,sigcm,decdist,m2final
         if(jj.lt.0.and.icase.eq.5) 
     >     write(6,'(''simc 5'',2f8.3)') decdist,m2final
c         else
c order of yptar, then xptar fixed 4/25/2020
c          read(string,'(e12.4,21f9.4,2e12.4)',end=15,err=999) weight,
c     >     dpe,dphie,dthe,hztar,hdcx,hdcxp,hdcy,hdcyp,
c     >     dpp,dphip,dthp,pztar,pdcx,pdcxp,pdcy,pdcyp,
c     >     dpp_init, xptar_init, yptar_init,yt_orig,yrast
c         endif
         if(weight.lt.0.) then
          write(6,'(''error, neg. wieight'',i2,e10.2)') icase,
     >     weight
          write(6,'(a)') string(1:50)
          write(6,'(a)') string(51:100)
          write(6,'(a)') string(101:150)
          write(6,'(a)') string(151:200)
          write(6,'(a)') string(201:250)
          weight=0.
         endif
         if(weight.gt.1.0.and.icase.eq.2) then
          write(6,'(''error, big giantic wieight'',i2,e10.2)') icase,
     >     weight
          weight=0.
         endif
         sum3 = sum3 + 1.
         sum4 = sum4 + sigcc
         dthe = dthe / 100.
c need to change sign!
c         dphie = -1. * dphie
         dphie = dphie / 100.
         dthp = dthp / 100.
         dphip = dphip / 100.
         pdcxp = pdcxp / 100.
         pdcyp = pdcyp / 100.
         hdcxp = hdcxp / 100.
         hdcyp = hdcyp / 100.
         xaero = pdcx + 231. * pdcxp
         yaero = pdcy + 231. * pdcyp

         okhms = 0
         okshms = 0
         if(abs(dpe).lt.13 .and.
     >    abs(dthe).lt.0.100 .and. 
     >    abs(dphie).lt.0.030) then
          ith = int((dthe + 0.100) / 0.200 * 70)+1
          iphi = int((dphie + 0.030) / 0.060 * 30.)+1
          ipt = int((dpe + 13.) / 26. * 30.)+1
          if(icase.eq.1) then
c took out adding excl to acc list
c           accep(2,ipt,ith,iphi) = accep(2,ipt,ith,iphi) + 1
          endif
          okhms = okacc(1,ipt,ith,iphi)
         endif
         if(abs(dpp-5.).lt.17 .and.
     >    abs(dthp).lt.0.100 .and. 
     >    abs(dphip).lt.0.030) then
          ith = int((dthp + 0.100) / 0.200 * 70)+1
          iphi = int((dphip + 0.030) / 0.060 * 30.)+1
          ipt = int((dpp + 12.) / 34. * 30.)+1
          if(icase.eq.1) then
           accep(4,ipt,ith,iphi) = accep(4,ipt,ith,iphi) + 1
          endif
          okshms = okacc(2,ipt,ith,iphi)
         endif
c use same th and phi range for dpp<-12 as for -12
         if(dpp.lt.-12 .and.
     >    abs(dthp).lt.0.100 .and. 
     >    abs(dphip).lt.0.030) then
          ith = int((dthp + 0.100) / 0.200 * 70)+1
          iphi = int((dphip + 0.030) / 0.060 * 30.)+1
          ipt = 1
          okshms = okacc(2,ipt,ith,iphi)
         endif

         if(dphie.gt.0.) iq2bin=1
         if(dphie.le.0.) iq2bin=2
c just one bin now
c        iq2bin=1

         the = the0 - dthe
         ep = ep0 * (1. + dpe/100.)
c wrong way
c         p_x(1) =  ep * sin(the) * cos(dphie)
c         p_y(1) = -ep * sin(the) * sin(dphie)
c         p_z(1) =  ep * cos(the) 
         ppi = abs(pp) * (1. + dpp/100.)
         thpi = thp * 3.1415928/180. + dthp
c         p_x(2) = -ppi * sin(thpi) * cos(dphip)
c         p_y(2) = -ppi * sin(thpi) * sin(dphip)
c         p_z(2) =  ppi * cos(thpi) 
c correct way 
	 call physics_angles(the0, phi0e, dthe, dphie, 
     >     ep,p_x(1),p_y(1),p_z(1))
	 call physics_angles(thp_rad, phi0p, dthp, dphip,
     >     ppi,p_x(2),p_y(2),p_z(2))

         nu = e0 - ep
         epv(1)=p_x(1)
         epv(2)=p_y(1)
         epv(3)=p_z(1)
         epv(4)=ep
         p1vp(1)=p_x(2)
         p1vp(2)=p_y(2)
         p1vp(3)=p_z(2)
         p1vp(4)= sqrt(ppi**2 + ampi**2)
         zpi = p1vp(4) / nu
         call getphi(e0,epv,p1vp,phicm)
         call getcos(e0,epv,p1vp,ampi,cthcm,pt)
! check
         Empi = e0 + amp - ep - sqrt(ampi**2 + ppi**2)
         mmpi2 = 
     >     (Empi)**2 - 
     >     (p_x(1) + p_x(2))**2 -
     >     (p_y(1) + p_y(2))**2 -
     >     (p_z(1) + p_z(2)- e0)**2
         if(jj.lt.0.and.icase.eq.1)
     >     write(6,'(7f8.3)') p_x(1),p_y(1),p_z(1),
     >               p_x(2),p_y(2),p_z(2),mmpi2
         if(mmpi2.gt.0.87.and.mmpi2.lt.0.919.and.icase.eq.1) then
          k = int((mmpi2-0.87)/0.05 * 20)+1
          if(k.gt.20) write(6,'(''error k mmpi2'',i3,f8.3)') k,mmpi2
          mmpi2hs(k) = mmpi2hs(k)+1
         endif

! acceptance cuts
         if(dpe.gt.dphmslo .and. dpe.lt.dphmshi .and.
     >      abs(dthe).lt.dthhms .and.
     >      abs(dthp).lt.dthshms .and.
     >      abs(dphie).lt.dphihms .and.
     >      abs(dphip).lt.dphishms .and.
c     >       okhms.eq.1 .and. okshms.eq.1 .and.
     >       pexit.eq.1 .and. hexit.eq.1.and. 
c add cuts at aerogel xy
     >      abs(xaero).lt.xmaxaero .and. abs(yaero).lt.ymaxaero.and.
c for kaons, put momentum cut
c     >      (icase.lt.5 .or. ppi.lt.2.8) .and.
     >      dpp.gt.dpshmslo .and. dpp.lt.dpshmshi) then

c final state counter
          ipart = 0
          if(abs(m2final-m_mu).lt..01) ipart=1
          if(abs(m2final-m_pi).lt..01) ipart=2
          if(abs(m2final-m_k).lt..11) ipart=3
          idist=min(22,max(1,int((decdist/0.1)+1)))
          if(ipart.ne.0) fstate(icase+1,idist,ipart) = 
     >      fstate(icase+1,idist,ipart) +weight
          fstate(icase+1,idist,4) = 
     >     fstate(icase+1,idist,4) + weight
          if(jj.lt.0 .and.icase.eq.5) then
           write(6,'("f 5",3i3,f6.2,2e12.2)') ipart,
     >      idist,icase,m2final,
     >      fstate(icase+1,idist,ipart),
     >      fstate(icase+1,idist,4)
          endif

          sum1 = sum1 + weight
          sum2 = sum2 + 1.
          im = max(1,min(20,int(mmpi2/0.5)+1))
          sum1m(im) = sum1m(im) + weight
          sum2m(im) = sum2m(im) + 1
          im = max(1,min(20,int(zpi/0.05)+1))
          sum1z(im) = sum1z(im) + weight
          sum2z(im) = sum2z(im) + 1
          pt2 = pt**2
          if(zpi.lt.1.0 .and. sqrt(pt2) .lt. 1.00) then
           iz = int(zpi*20.)+1
           izz = int(zpi*50.)+1
           im = max(1,min(50,int(mmpi2/0.2)+1))
           iphi = max(1,min(15,int(phicm/6.29*15 ) + 1)) 
c           ipt = int(pt2 / 0.5 * 12) + 1
           ipt = int(sqrt(pt2) / 1.00 * 16.) + 1
           i3 = 1 + 2*icase
           i4 = 2 + 2*icase
           cntsmc(iq2bin,ipt,iphi,iz,i3) = 
     >       cntsmc(iq2bin,ipt,iphi,iz,i3) + 1
           cntsmc(iq2bin,ipt,iphi,iz,i4) = 
     >       cntsmc(iq2bin,ipt,iphi,iz,i4) + weight
           cntsmcmmpi(im,i3) = cntsmcmmpi(im,i3) + 1
           cntsmcmmpi(im,i4) = cntsmcmmpi(im,i4) + weight
           cntsmcz(izz,i3) = cntsmcz(izz,i3) + 1
           cntsmcz(izz,i4) = cntsmcz(izz,i4) + weight
c get sig for pi_norad, k_norad
           if(icase.eq.4) then
            cntsmc(iq2bin,ipt,iphi,iz,15) = 
     >        cntsmc(iq2bin,ipt,iphi,iz,15) + sigcc
            cntsmcmmpi(im,15) = cntsmcmmpi(im,15)+sigcc
            cntsmcz(izz,15) = cntsmcz(izz,15) + sigcc
           endif
           if(icase.eq.6) then
            cntsmc(iq2bin,ipt,iphi,iz,16) = 
     >        cntsmc(iq2bin,ipt,iphi,iz,16) + sigcc
            cntsmcmmpi(im,16) = cntsmcmmpi(im,16)+sigcc
            cntsmcz(izz,16) = cntsmcz(izz,16) + sigcc
           endif
c  with mmpi2 cut
c           if(mmpi2.gt.2.0) then
c            cntsmc(iq2bin,ipt,iphi,iz,7) = 
c     >       cntsmc(iq2bin,ipt,iphi,iz,7) + 1
c            cntsmc(iq2bin,ipt,iphi,iz,8) = 
c     >       cntsmc(iq2bin,ipt,iphi,iz,8) + weight
c           endif
          endif
c exclusive pion counts
          if(icase.eq.1 .and. mmpi2.gt.0.76 .and.
     >     mmpi2.lt.1.0 .and.
     >     sqrt(pt2).lt.0.3) then
           iphi = max(1,min(15,int(phicm/6.29*7 ) + 1)) 
           ipt = int(sqrt(pt2) / 0.3 * 4.) + 1
           cntsemc(1,ipt,iphi,1) = 
     >       cntsemc(1,ipt,iphi,1) + 1 
           cntsemc(1,ipt,iphi,2) = 
     >       cntsemc(1,ipt,iphi,2) + weight
           cntsemc(1,ipt,iphi,3) = 
     >       cntsemc(1,ipt,iphi,3) + sigcm
           im = max(1,min(8,int((mmpi2-0.76)/0.24*8.)+1))
           cntsemch(1,ipt,iphi,im) = 
     >        cntsemch(1,ipt,iphi,im) + 1
          endif
         endif
        enddo
 15     sratex(i,icase) = sum1 * normfacx / max(1.,sum3)
        sratexer(i,icase) = sratex(i,icase) / sqrt(max(1.,sum2))
        do im=1,20
         sratemx(im,icase) = sum1m(im) * normfacx / sum3
         sratemxer(im,icase) = sratemx(im,icase) / 
     >     sqrt(max(1.,sum2m(im)))
         sratezx(im,icase) = sum1z(im) * normfacx / sum3
         sratezxer(im,icase) = sratezx(im,icase) / 
     >     sqrt(max(1.,sum2z(im)))
        enddo
        write(6,'(''normfac_x='',i2,5e10.3)') icase,
     >    normfacx,sratex(i,icase),sum1,sum2,sum3
        i4 = 2 + 2*icase
        avnfac(icase+1) = avnfac(icase+1) + normfacx * chrg
        sum3tot(icase+1) = sum3tot(icase+1) + sum3
 16     continue
        enddo ! loop over icase

! write out diagnostic histograms
        write(43,'(i4,20i3)') irun,(int(100.*ctrfrun(k,1)/
     >   float(ctrfrun(21,1))),k=1,20)
c        write(43,'( 20i4)') (ctrfrun(k,2)/10,k=2,18)
c this is to find position of cointime peak
        do ii=1,2
        write(44,'(i4)') irun
        do kk=1,10
         sum1=0.
         sum2=0.
         do k=1,20
          sum1 = sum1 + k * ctrfrn(kk,k,ii)
          sum2 = sum2 +     ctrfrn(kk,k,ii)
         enddo
         sum1 = sum1 / sum2
         sum3 = 0.
         do k=1,20
          sum3 = sum3 + (k-sum1)**2 * ctrfrn(kk,k,ii)
         enddo
c         write(6,'(2f8.3)') sum1,sqrt(sum3/sum2)
         sigsv(kk) = sqrt(sum3/sum2)
        enddo
        sum1=1000.
        do kk=1,10
         if(sigsv(kk).lt.sum1) then
          sum1 = sigsv(kk)
          sum2 = 4.002 + 0.001*kk
         endif
        enddo
        write(44,'(i4,f6.3,11f6.2)') irun,sum2,pp,
     >    (sigsv(kk),kk=1,10)
        enddo

        write(45,'(/i4,4f7.2)') irun,ep0,the0*57.,pp,thp
        do k=1,13
         write(45,'(16i4)') (yth(k,kk),kk=1,16)
        enddo

        write(81,'(i4,14f5.0)') irun,
     >   (1000.*float(dpphistrun(k))/
     >    float(dpphistrun(15)),k=1,14)
        write(82,'(i4,14f5.0)') irun,
     >   (1000.*float(ythistrun(k))/
     >    float(ythistrun(15)),k=1,14)

        write(46,'(i5)') irun
c hbeta is not read in
c        write(46,'(20i4)') (betah(1,k)/10,k=1,20)
        write(46,'(i5,19i4)') betah(2,20),
     >    (int(100.*betah(2,k)/betah(2,20)),k=1,19)

        write(47,'(/i5,i2,f7.1,f6.2)') 
     >    irun,it,irate_p/1000.,pp
        write(47,'(i5,15i4)') aeroh(16),
     >    (int(1000.*aeroh(k)/aeroh(16)),k=1,15)
        write(47,'(i5,15i4)') aerohp(16),
     >    (int(1000.*aerohp(k)/aerohp(16)),k=1,15)
c        do kk=1,4
c         write(47,'(20i4)') (int(1000.*xyaeroh(kk,k)/xyaeroh(kk,21)),
c     >    k=1,20)
c        enddo
        write(48,'(i5)') irun
        write(48,'(20i4)') (int(1000.*shh(k)/shh(16)),k=1,15)
        write(52,'(i5)') irun
        write(52,'(20i4)') (int(1000.*sheh(k)/sheh(16)),k=1,15)
        write(53,'(i5)') irun
        write(53,'(20i4)') (int(1000.*sh2h(k)/sh2h(16)),k=1,15)
        write(59,'(/i5,f6.2)') irun,pp
        write(59,'(20i4)') (int(1000.*prh(k,1)/prh(16,1)),k=1,15)
c note: with shp cut still divided by 16,1
        write(59,'(20i4)') (int(1000.*prh(k,2)/prh(16,1)),k=1,15)
        write(54,'(i5)') irun
        write(54,'(16i4)') ngh(16),
     >   (int(1000.*ngh(k)/ngh(16)),k=1,15)

        write(51,'(/i5)') irun
c        do k=1,12
        do k=7,7
         write(51,'(i6,15i4)') cereh(k,16),
     >   (int(1000.*float(cereh(k,kk))/
     >    float(cereh(k,16))),kk=1,15)
        enddo
c        do k=1,12
        do k=7,7
c         write(51,'(i6,15i4)') cerepih(k,16),
c     >   (int(1000.*float(cerepih(k,kk))/
c     >    float(cerepih(k,16))),kk=1,15)
        enddo

        sum3=0.
        do k=1,39
         sum1 = (ctpih(k) * k +
     >      ctpih(k+1) * (k+1) ) 
         sum2 = (ctpih(k)  +
     >      ctpih(k+1) ) 
         if(sum2.gt.sum3) then
          sum3 = sum2
          pipeak = sum1 / sum2
         endif
        enddo
        write(6,'(i4,7i6)') irun,(ctpih(kk),kk=17,23)

c Electronic dead time not measured by EDTM
c casued by having elclean in ref. time signal 150 nsec
c after the 3/4 signal
c correction is bigger for sp18 due to 100 nsec long hodo
c pulses compared to 35 (?) after that. 
c This is based on fitting the current scans.
        corr2 = 1.0
        if(irun.lt.4400) then
         corr2 = 1.0 - 0.280 * (float(irate_p)/1.e6 + elreal/1000.)
        endif
        if(irun.gt.4400) then
         corr2 = 1.0 - 0.050 * (float(irate_p)/1.e6 + elreal/1000.)
        endif

! timecorr currently using nevent32
! Get rate averaged over run
        timecorr = 1. 
        ratec(i) = (pir - piacc/4.) / chrg / dt / 
     >   tefe / tefp / corr2 / corr3 / corr
        ratecer(i) = sqrt(pir + piacc/16.)/chrg/dt/
     >    tefe / tefp / corr2 / corr3 / corr
        chrgtot = chrgtot + chrg
        dtav = dtav + dt * chrg
        tefeav = tefeav + tefe * chrg
        tefpav = tefpav + tefp * chrg
        corrav = corrav + corr * chrg
        corr2av = corr2av + corr2 * chrg
        corr3av = corr3av + corr3 * chrg

        do im=1,20
         ratem(im) = (pirm(im) - piaccm(im)/4.) / chrg / dt / 
     >    tefe / tefp / corr3 / corr / corr2
         ratemer(im ) = sqrt(pirm(im) + piaccm(im)/16.)/chrg/dt/
     >    tefe / tefp / corr3 / corr / corr2
         ratez(im) = (pirz(im) - piaccz(im)/4.) / chrg / dt / 
     >    tefe / tefp / corr3 / corr / corr2
         ratezer(im ) = sqrt(pirz(im) + piaccz(im)/16.)/chrg/dt/
     >    tefe / tefp / corr3 / corr / corr2
        enddo
        fact = 1. / chrg / dt / 
     >   tefe / tefp / corr3 / corr / corr2
        write(91,'(''irun='',i5)') irun
        do k=1,6
         do kk=1,20
          rateck(k,kk) = (pirk(k,kk) - piacck(k,kk)/4.) * fact
          rateak(k,kk) =               piacck(k,kk)/4.  * fact
          ratecker(k,kk) = fact * sqrt(pirk(k,kk) + 
     >     piacck(k,kk)/16.)
         enddo
         write(91,'(21f5.0)') ratec(i),(rateck(k,kk),kk=1,20)
         write(91,'(20i4)') (int(1000.*rateck(k,kk)/ratec(i)),kk=1,20)
         write(91,'(20i4)') (int(1000.*rateak(k,kk)/ratec(i)),kk=1,20)
         write(91,'(20i4)') (int(1000.*ratesk(k,kk)/srate(i)),kk=1,20)
         write(91,'(/)')
        enddo
        if(it.ne.itprev .or. 
     >   abs(thp - thpprev).gt.0.1 .or.
     >   abs(pp - ppprev).gt.0.1 .or.
     >   abs(ep0 - ep0prev).gt.0.1.or.
     >   abs(the0*57.3 -the0prev).gt.0.1) then
         if(nrtc.gt.0) then
          avrtc = avrtc / avrtcer
          avrtcer = 1./sqrt(avrtcer)
          chi = 0.
          df = 0.
          do k=1,nrtc
           chi = chi + (rtc(k) - avrtc)**2 / rtcer(k)**2
           df = df + 1.
          enddo
          write(18,'(6x,38x,2f6.1,f6.1,f4.0)')
     >      avrtc,avrtcer,chi/df,df
         endif
         avrtc = 0.
         avrtcer = 0.
         nrtc = 0
         itprev = it
         ppprev = pp
         thpprev = thp
         ep0prev = ep0
         the0prev = the0*57.3
         write(18,'()')
        endif ! check if changing kinematics
         nrtc = nrtc + 1
c took out * corr2 not sure why it was there!
         rtc(nrtc) = ratec(i) 
         rtcer(nrtc) = ratecer(i)
         avrtc = avrtc + rtc(nrtc)/rtcer(nrtc)**2
         avrtcer = avrtcer + 1./rtcer(nrtc)**2
        write(18,
     >  '(i4,i2,f4.1,f4.0,f4.1,f4.0,2f6.0,2f5.1,2f6.1,f5.2)') irun,it,
     >    ep0,the0*57.,pp,thp,pir,piacc/4.,chrg,
     >    current,ratec(i),ratecer(i),tefp
        write(7,'(i5,i2,f6.2,40i4)') irun,it,pipeak,
     >    (ctpih(k)/10,k=1,40)
c     >    (ctpirh(k),k=1,20)
        write(77,'(/i5)') irun
        do kk=1,10
         write(77,'(20i4)') (ctpihw(k,kk)/30,k=12,31)
        enddo
        do k=1,30
           write(36,'(i3,i7)') k,ctpih(k)
        enddo
        write(36,'(''hist'')')
        write(20,'(i5,f8.2)') irun,pipeak-10.5
        if(abs(pipeak-10.5).gt.0.25) writE(20,'(''WARNING'')')
        sum3=0.
        do k=2,19
         sum1 = (ctpihf(k) * k +
     >      ctpihf(k+1) * (k+1) ) 
         sum2 = (ctpihf(k)  +
     >      ctpihf(k+1) ) 
         if(sum2.gt.sum3) then
          sum3 = sum2
          pipeak = sum1 / sum2
         endif
        enddo
        sum3 = 0.
        do k=6,15
         sum3 = sum3 + ctpihfa(k)
        enddo
        write(17,'(/''p ='',f5.2,i5,i3,f6.2)') pp,irun,it,pipeak
        write(17,'(i5,38i3)') int(sum3),
     >    (int(400*ctpihfa(k)/sum3),k=3,40)
        sum3=0.
        do k=6,15
         sum3 = sum3 + ctpihf(k)
        enddo
        write(17,'(i5, 38i3)') int(sum3),
     >    (int(400*ctpihf(k)/sum3),k=3,40)

        sum3 = 0.
        do k=6,15
         sum3 = sum3 + ctkhf(k)
        enddo
        write(17,'(i5,38i3)') int(sum3),
     >    (int(400*ctkhf(k)/sum3),k=3,40)

        sum3 = 0.
        do k=6,15
         sum3 = sum3 + ctphf(k)
        enddo
        write(17,'(i5,38i3)') int(sum3), 
     >    (int(400*ctphf(k)/sum3),k=3,40)


        write(8,'(i5,20i4)') irun,(ctkh(k),k=1,20)
        write(33,'(i5,20i4)') irun,(ctkhf(k),k=1,20)
        write(19,'(i5,20i4)') irun,(cteh(k),k=1,20)
        if(it.eq.1) then
         write(31,'(''irun,it='',i5,i2)') irun,it
         do kk=1,12
          write(31,'(i2,20i4)') kk,(mmph(k,kk),k=1,20)
         enddo
        endif

        write(32,'(i5,20i4)') irun,(meeh(k),k=1,20)
c old way
        rate1 = icoin / chrg / dt
        rate(i)=rate1
        fcoin = icoin
        rateer(i)=sqrt(fcoin)/chrg/dt
c new way
        rate(i) = ratec(i)
        rateer(i) = ratecer(i)
cxxx test
        irunlast = irun
        write(16,'(i4,f5.1,f6.1,f5.1,2f6.2,3f5.2,
     >   6f5.2,f6.0,2f6.0,2f7.3,2f8.0)') 
     >   irun,current,
     >   ratec(i),ratecer(i),srate(i)/ratec(i),
c     >   (srate(i) + sratex(i,1))/ratec(i),
c     >   (srate(i) + sratex(i,1) + sratex(i,2))/ratec(i),
c     >   (srate(i) + sratex(i,1) + sratex(i,2) + sratex(i,3))/ratec(i),
c just excl. and endcap, no rho
     >   (srate(i) + sratex(i,1) + sratex(i,3))/ratec(i),
     >   dt,tefe,tefp,piacc/(pir-piacc/4.)/current*10.,
     >   corr3,corr,corr2,heffh,heffp,
     >   float(irate_p)/1000.,pir,piacc/4.,
     >   piaccL/(pir-piaccL/2.)/current*100.,
     >   piaccH/(pir-piaccH/2.)/current*100.,
     >   float(ielclean)/chrg/100., float(ielcleanp)/chrg/100.
        nrt = nrt + 1
        do kk=1,10
c does use irate_p  + elreal work better?
        corr4 = 1.0 + 0.025 * float(kk-1) * 
     >    (float(irate_p) / 1.E6 + elreal / 1000.)
         rt(nrt,kk) = ratec(i) * corr4
         rter(nrt,kk) = ratecer(i) * corr4
         avrt(kk) = avrt(kk) + ratec(i)/ratecer(i)**2 / corr4
         avrter(kk) = avrter(kk) + 1./ratecer(i)**2 / corr4**2
        enddo
        savrtnox = savrtnox + srate(i)/srateer(i)**2
        savrtnoxer = savrtnoxer + 1./srateer(i)**2
! here the four SIMC results are combined!
! changed to exclude the rho
c        srate(i) = srate(i) + sratex(i,1) + sratex(i,2) + sratex(i,3) 
c        srateer(i) = sqrt(srateer(i)**2 + sratexer(i,1)**2 + 
c     >   sratexer(i,2)**2 + sratexer(i,3)**2)
        srate(i) = srate(i) + sratex(i,1) + sratex(i,3) 
        srateer(i) = sqrt(srateer(i)**2 + sratexer(i,1)**2 + 
     >   sratexer(i,3)**2)
        srt(nrt) = srate(i)
        srter(nrt) = srateer(i)
        savrt = savrt + srate(i)/srateer(i)**2
        savrter = savrter + 1./srateer(i)**2
        if(current.gt.curmax) curmax = current
        if(current.lt.curmin) curmin = current
! Special case: last run
        if(irun.eq.6559) then
         avrate(8,3,11) = ratec(i)
         avrateer(8,3,11) = ratecer(i)
        endif

! compare to simc
c        write(6,
c     >   '(i4,e10.2,f6.2,e10.2,f6.2,f8.3)')
c     >   i,rate(i),rateer(i)/rate(i),
c     >   srate(i),srateer(i)/srate(i),
c     >   srate(i)/rate(i)

c hist of rad. ep elastic events
        write(85,'(''run, it '',i4,i3)') irun, it
        do k=1,14 
         write(85,'(i2,i5,20i3)') k,eprh(k,21),
     >    ((100*eprh(k,kk))/max(1,eprh(k,21)),kk=1,20)
        enddo

c hist. of mmpi
        if(mmpi2_cent_h.lt.1.8) then
c         write(15,'(''run, it, pp '',i4,i3,f6.2)') irun, it,pp
         do k=1,20
          sumk(k)=0
         enddo
         do kk=1,12 
c          write(15,'(20i4)') (mmpih(k,kk)/10,k=1,12)
          do k=1,20
           sumk(k) = sumk(k)+mmpih(k,kk)
          enddo
         enddo
         if(irun.lt.3422) write(15,'(4x,12f5.2)') 
     >     (10.*(0.85 + 0.010*(k-0.5)),k=4,15)
         write(15,'(i4,20i5)') irun,(sumk(k)/10,k=4,15)
        endif

c hist. of phicm
        write(14,'(/''run, it '',i4,i3,2f8.3)') irun, it,pp,thp
        write(14,'(i5,15i3)') phih(20),
     >    (int(100.*phih(k)/phih(20)),k=1,15)
        write(144,'(i4,2f8.3,6f6.0)') irun,
     >   (sumcntsf(1,1) - sumcntsf(2,1))/
     >   (sumcnt(1,1)+sumcnt(2,1)),
     >   1./sqrt(sumcnt(1,1)+sumcnt(2,1)),
     >   sumcnt(1,1), sumcnt(2,1), sumcnt(3,1),
     >   sumcnt(1,2), sumcnt(2,2), sumcnt(3,2)


       write(55,'(/i5,i3,2f7.2)') irun,it,pp,thp
       write(55,'(20f4.1)') (min(9.9, max(0.,
     >   ratem(im)/max(1.e-20,sratem(im)))),im=1,20)
       write(55,'(20f4.1)') (min(9.9, max(0.,
     >   ratem(im)/max(1.e-20,
     >   (sratem(im)+sratemx(im,1))))),im=1,20)
       write(55,'(20f4.1)') (min(9.9, max(0.,
     >   ratemer(im)/max(1.e-20,
     >   (sratem(im)+sratemx(im,1))))),im=1,20)
       do im=1,20
        write(55,'(i2,9f7.2)') im,
     >   ratem(im),ratemer(im),
     >   sratem(im),sratemer(im),
     >   sratemx(im,1),sratemxer(im,1),
     >   min(9.9, max(0.,
     >   ratem(im)/max(1.e-20,
     >   sratem(im)))),
     >   min(9.9, max(0.,
     >   ratem(im)/max(1.e-20,
     >   (sratem(im)+sratemx(im,1))))),
     >   min(9.9, max(0.,
     >   ratemer(im)/max(1.e-20,
     >   (sratem(im)+sratemx(im,1)))))
       enddo
       write(56,'(/i5,i3,2f7.2)') irun,it,pp,thp
       do im=5,20
        write(56,'(i2,6f6.1,3f6.2)') im,
     >   ratez(im),ratezer(im),
     >   sratez(im),sratezer(im),
     >   sratezx(im,1),sratezxer(im,1),
     >   min(99.9, max(0.,
     >   ratez(im)/max(1.e-20,
     >   sratez(im)))),
     >   min(99.9, max(0.,
     >   ratez(im)/max(1.e-20,
     >   (sratez(im)+sratezx(im,1))))),
     >   min(99.9, max(0.,
     >   ratezer(im)/max(1.e-20,
     >   (sratez(im)+sratezx(im,1)))))
       enddo

       endif ! check on run range to analyze
      enddo ! loop over runs
      write(6,'(''total chi2 df'')')
      do kk=1,10
       write(6,'(i3,2f10.1,f8.3)') kk,totchi2(kk),totdf(kk),
     >   totchi2(kk)/totdf(kk)
      enddo
      write(32,'(i5,20i4)') irun,(meeht(k,1),k=1,20)
      write(32,'(i5,20i4)') irun,(meeht(k,2),k=1,20)
      do k=1,20
         write(32,'(f7.3,i5,f6.1)') 3.0  + 0.010 * (k-0.5),
     >    meeht(k,2),sqrt(float(meeht(k,2)))
      enddo
      do k=1,20
         write(32,'(f7.3,i5,f5.2)') 2.80 + 0.020 * (k-0.5),
     >    meeht(k,1),sqrt(float(meeht(k,1)))
      enddo


c Plot of coin time pi for 20 momenta, four conditions
      close(unit=22)
      open(unit=22,file='ptcct.top')
      write(22,'(''set device postscript'')')
      do ip=1,20
       ix = int((ip+4)/5)
       iy = ip - 5*(ix-1)
       imax = 0
       do k=1,40
        do j=1,4
         imax = max(imax,cthist(ip,j,k))
        enddo
       enddo
       write(22,122) ix,iy,imax,1.75+0.1*ip
 122   format(1x,'set window x ',i1,' of 4 y ',i1,' of 5'/
     >   1x,'set limits y 0. ',i8/
     >   1x,'title top ',1h','P=',f5.2,1h')
       do k=1,40
        write(22,'(1x,f7.2,i8)') -2.+0.2*(k-0.5), cthist(ip,1,k)
       enddo
       write(22,'(''join'')')
       do k=1,40
        write(22,'(1x,f7.2,i8)') -2.+0.2*(k-0.5), cthist(ip,2,k)
       enddo
       write(22,'(''set pattern .05 .05 .05 .05 ; join pattern'')')
       do k=1,40
        write(22,'(1x,f7.2,i8)') -2.+0.2*(k-0.5), cthist(ip,3,k)
       enddo
       write(22,'(''set pattern .02 .02 .02 .02 ; join pattern'')')
       do k=1,40
        write(22,'(1x,f7.2,i8)') -2.+0.2*(k-0.5), cthist(ip,4,k)
       enddo
       write(22,'(''set pattern .02 .02 .06 .06 ; join pattern'')')
      enddo

      do ipm=1,2
      close(unit=22)
      if(ipm.eq.1) open(unit=22,file='ptcctrf.top')
      if(ipm.eq.2) open(unit=22,file='ptcctrfneg.top')
      write(22,'(''set device postscript'')')
      do ip=1,20
       ix = int((ip+4)/5)
       iy = ip - 5*(ix-1)
       imax = 0
       do k=1,40
        do j=1,4
         imax = max(imax,ctrfhist(ip,j,k,ipm))
        enddo
       enddo
       write(22,122) ix,iy,imax,1.75+0.1*ip
       do k=1,40
        write(22,'(1x,f7.2,i8)') 0.1*(k-0.5), ctrfhist(ip,1,k,ipm)
       enddo
       write(22,'(''join'')')
       do k=1,40
        write(22,'(1x,f7.2,i8)') 0.1*(k-0.5), ctrfhist(ip,2,k,ipm)
       enddo
       write(22,'(''set pattern .05 .05 .05 .05 ; join pattern'')')
       do k=1,40
        write(22,'(1x,f7.2,i8)') 0.1*(k-0.5), ctrfhist(ip,3,k,ipm)
       enddo
       write(22,'(''set pattern .02 .02 .02 .02 ; join pattern'')')
       do k=1,40
        write(22,'(1x,f7.2,i8)') 0.1*(k-0.5), ctrfhist(ip,4,k,ipm)
       enddo
       write(22,'(''set pattern .02 .02 .06 .06 ; join pattern'')')
      enddo
      enddo ! ipm

c      write(31,'(i5,i2,20i4)') irun,it,(mmphs(k),k=1,20)
      do k=1,400
         write(31,'(f7.3,i7,f8.2,i7,f8.2)')  0.005 * (k-0.5),
     >    mmphfs(k),sqrt(float(mmphfs(k))),
     >    mmpheta(k),sqrt(float(mmpheta(k)))
      enddo

      do k=1,20
         write(15,'(f7.3,i7,f8.2,i7,f8.2)') 
     >    0.85 + 0.010 * (k-0.5),
     >    mmpihs(k),sqrt(float(mmpihs(k))),
     >    mmpihsa(k),sqrt(float(mmpihsa(k)))
      enddo

c plots of on-line rates
      open(unit=22,file='ptcr1.top')
      write(22,'(''set device postscript'')')
      do ix=1,2
       do iy=1,2
        write(22,101) ix,iy
 101    format(1x,
     >   1x,'set window x ',i1,' of 2 y ',
     >    i1,' of 2'/ 
     >   1x,'set font duplex'/
     >   1x,'set intensity 4'/
     >   1x,'set bar size 0.'/
     >   1x,'set order x y dy SYM'/
     >   1x,'set sym 9O size 1.0'/
     >   1x,'set scale y log'/
     >   1x,'title left ',1h',' ',1h'/
     >   1x,'title ',1h','ratec',1h'/
     >   1x,'title bottom ',1h','theta (deg)',1h'/
     >   1x,'set ticks size 0.05')
        if(ix.eq.1.and.iy.eq.2) write(22,102)
 102    format(1x,'title top',1h',
     >    'x=0.3 Q2=3.0 P= 2.53',1h')
        if(ix.eq.1.and.iy.eq.1) write(22,103)
 103    format(1x,'title top',1h',
     >    'x=0.3 Q2=3.0 P= -2.53',1h')
        if(ix.eq.2.and.iy.eq.2) write(22,104)
 104    format(1x,'title top',1h',
     >    'x=0.3 Q2=3.0 P= 1.96',1h')
        if(ix.eq.2.and.iy.eq.1) write(22,105)
 105    format(1x,'title top',1h',
     >    'x=0.3 Q2=3.0 P= -1.96',1h')
        do i=1,5000
         if(abs(kin(i,1)-5.27).lt.0.1.and.
     >      kin(i,3).gt.2.2.and.
     >      kin(i,3).lt.2.6 .and.
     >      ix.eq.1 .and.iy.eq.2) then
          write(22,'(3f10.4,i4,''O'')') 
     >      kin(i,4) -0.03+ 0.006*( i - int(i/10) * 10),
     >      ratec(i),ratecer(i),itsv(i)
         endif
         if(abs(kin(i,1)-5.27).lt.0.1.and.
     >      kin(i,3).lt.-2.2.and.
     >      kin(i,3).gt.-2.6 .and.
     >      ix.eq.1 .and.iy.eq.1) then
          write(22,'(3f10.4,i4,''O'')') 
     >      kin(i,4) -0.03+ 0.006*( i - int(i/10) * 10),
     >      ratec(i),ratecer(i),itsv(i)
         endif
         if(abs(kin(i,1)-5.27).lt.0.1.and.
     >      kin(i,3).gt.1.9.and.
     >      kin(i,3).lt.2.1 .and.
     >      ix.eq.2 .and.iy.eq.2) then
          write(22,'(3f10.4,i4,''O'')') 
     >      kin(i,4) -0.03+ 0.006*( i - int(i/10) * 10),
     >      ratec(i),ratecer(i),itsv(i)
         endif
         if(abs(kin(i,1)-5.27).lt.0.1.and.
     >      kin(i,3).lt.-1.9.and.
     >      kin(i,3).gt.-2.1 .and.
     >      ix.eq.2 .and.iy.eq.1) then
          write(22,'(3f10.4,i4,''O'')') 
     >      kin(i,4) -0.03+ 0.006*( i - int(i/10) * 10),
     >      ratec(i),ratecer(i),itsv(i)
         endif
        enddo
        write(22,'(''plot'')')
       enddo
      enddo
      close(unit=22)
      open(unit=22,file='ptcr2.top')
      write(22,'(''set device postscript'')')
c      do ix=1,2
      do ix=1,2
       do iy=1,2
        write(22,111) ix,iy
 111    format(1x,
     >   1x,'set window x ',i1,' of 2 y ',
     >    i1,' of 2'/ 
     >   1x,'set font duplex'/
     >   1x,'set intensity 4'/
     >   1x,'set bar size 0.'/
     >   1x,'set order x y dy SYM'/
     >   1x,'set sym 9O size 1.0'/
     >   1x,'set scale y log'/
     >   1x,'title left ',1h',' ',1h'/
     >   1x,'title ',1h','ratec',1h'/
     >   1x,'title bottom ',1h','theta (deg)',1h'/
     >   1x,'set ticks size 0.05')
        if(ix.eq.1.and.iy.eq.2) write(22,112)
 112    format(1x,'title top',1h',
     >    'x=0.3 Q2=4.1 P= 3.40',1h')
        if(ix.eq.1.and.iy.eq.1) write(22,113)
 113    format(1x,'title top',1h',
     >    'x=0.3 Q2=4.1 P= -3.40',1h')
        if(ix.eq.2.and.iy.eq.2) write(22,114)
 114    format(1x,'title top',1h',
     >    'x=0.3 Q2=4.1 P= 2.65',1h')
        if(ix.eq.2.and.iy.eq.1) write(22,115)
 115    format(1x,'title top',1h',
     >    'x=0.3 Q2=4.1 P= -2.65',1h')
        do i=1,5000
         if(abs(kin(i,1)-3.32).lt.0.1.and.
     >      kin(i,3).gt.3.2.and.
     >      kin(i,3).lt.3.6 .and.
     >      ix.eq.1 .and.iy.eq.2) then
          write(22,'(3f10.4,i4,''O'')') 
     >      kin(i,4) -0.03+ 0.006*( i - int(i/10) * 10),
     >      ratec(i),ratecer(i),itsv(i)
         endif
         if(abs(kin(i,1)-3.32).lt.0.1.and.
     >      kin(i,3).lt.-3.2.and.
     >      kin(i,3).gt.-3.6 .and.
     >      ix.eq.1 .and.iy.eq.1) then
          write(22,'(3f10.4,i4,''O'')') 
     >      kin(i,4) -0.03+ 0.006*( i - int(i/10) * 10),
     >      ratec(i),ratecer(i),itsv(i)
         endif
         if(abs(kin(i,1)-3.32).lt.0.1.and.
     >      kin(i,3).gt.2.5.and.
     >      kin(i,3).lt.2.8 .and.
     >      ix.eq.2 .and.iy.eq.2) then
          write(22,'(3f10.4,i4,''O'')') 
     >      kin(i,4) -0.03+ 0.006*( i - int(i/10) * 10),
     >      ratec(i),ratecer(i),itsv(i)
         endif
         if(abs(kin(i,1)-3.32).lt.0.1.and.
     >      kin(i,3).lt.-2.5.and.
     >      kin(i,3).gt.-2.8 .and.
     >      ix.eq.2 .and. iy.eq.1) then
          write(22,'(3f10.4,i4,''O'')') 
     >      kin(i,4) -0.03+ 0.006*( i - int(i/10) * 10),
     >      ratec(i),ratecer(i),itsv(i)
         endif
        enddo
        write(22,'(''plot'')')
       enddo
      enddo

! Do dummy subtraction
      do ikin=1,8
       do ith=1,15
        if(avrate(ikin,1,ith).ne.0. .and.
     >     avrate(ikin,3,ith).ne.0.) then
         avrate(ikin,1,ith) =
     >    avrate(ikin,1,ith) - 
     >    0.262 * avrate(ikin,3,ith)
         avrateer(ikin,1,ith) = sqrt(
     >    avrateer(ikin,1,ith)**2 + 
     >    (0.262 * avrateer(ikin,3,ith))**2)
        endif
        if(avrate(ikin,2,ith).ne.0. .and.
     >     avrate(ikin,3,ith).ne.0.) then
         avrate(ikin,2,ith) =
     >    avrate(ikin,2,ith) - 
     >    0.262 * avrate(ikin,3,ith)
         avrateer(ikin,2,ith) = sqrt(
     >    avrateer(ikin,2,ith)**2 + 
     >    (0.260 * avrateer(ikin,3,ith))**2)
        endif
        if(avrate(ikin,4,ith).ne.0. .and.
     >     avrate(ikin,6,ith).ne.0.) then
         avrate(ikin,4,ith) =
     >    avrate(ikin,4,ith) - 
     >    0.262 * avrate(ikin,6,ith)
         avrateer(ikin,4,ith) = sqrt(
     >    avrateer(ikin,4,ith)**2 + 
     >    (0.262 * avrateer(ikin,6,ith))**2)
        endif
        if(avrate(ikin,5,ith).ne.0. .and.
     >     avrate(ikin,6,ith).ne.0.) then
         avrate(ikin,5,ith) =
     >    avrate(ikin,5,ith) - 
     >    0.262 * avrate(ikin,6,ith)
         avrateer(ikin,5,ith) = sqrt(
     >    avrateer(ikin,5,ith)**2 + 
     >    (0.260 * avrateer(ikin,6,ith))**2)
        endif
       enddo
      enddo

! Results
      do ikin=1,8
       do ith=1,15
        if(avrate(ikin,1,ith).ne.0.) then
         write(16,'(2i4,10f6.2)') ikin,ith,
     >    (avrate(ikin,itp,ith) /
     >     avrate(ikin,  1,ith),
     >     avrateer(ikin,itp,ith) /
     >     avrate(ikin,  1,ith),itp=2,6)
        endif
       enddo
      enddo

! Plots of p, th, phi
      open(unit=22,file='ptckin.top')
      write(22,'(''set device postscript'')')
      iy=1
      ix=1
      write(22,131) ix,iy
 131  format(1x,
     >   1x,'set window x ',i1,' of 3 y ',
     >    i1,' of 2'/ 
     >   1x,'set font duplex'/
     >   1x,'set intensity 4'/
     >   1x,'set bar size 0.'/
     >   1x,'set order x y dy SYM'/
     >   1x,'set sym 9O size 1.0'/
     >   1x,'set ticks size 0.05')
      write(22,132)
 132  format(
     >   1x,'title bottom ',1h','dpp',1h')
      do k=1,100
       write(22,'(f8.4,i10)') -50.+1.0*(k-0.5),
     >  dpph(k)/10
      enddo
      write(22,'(''hist'')')
      ix=2
      write(22,131) ix,iy
      write(22,133)
 133  format(
     >   1x,'title bottom ',1h','dthp',1h')
      do k=1,100
       write(22,'(f8.4,i10)') -50.+1.0*(k-0.5),
     >  dthph(k)/10
      enddo
      write(22,'(''hist'')')
      ix=3
      write(22,131) ix,iy
      write(22,134)
 134  format(
     >   1x,'title bottom ',1h','dphip',1h')
      do k=1,100
       write(22,'(f8.4,i10)') -50.+1.0*(k-0.5),
     >  dphiph(k)/10
      enddo
      write(22,'(''hist'')')
      iy=2
      ix=1
      write(22,131) ix,iy
      write(22,142)
 142  format(
     >   1x,'title bottom ',1h','dpe',1h')
      do k=1,100
       write(22,'(f8.4,i10)') -20.+0.4*(k-0.5),
     >  dpeh(k)/10
      enddo
      write(22,'(''hist'')')
      ix=2
      write(22,131) ix,iy
      write(22,143)
 143  format(
     >   1x,'title bottom ',1h','dthe',1h')
      do k=1,100
       write(22,'(f8.4,i10)') -50.+1.0*(k-0.5),
     >  dtheh(k)/10
      enddo
      write(22,'(''hist'')')
      ix=3
      write(22,131) ix,iy
      write(22,144)
 144  format(
     >   1x,'title bottom ',1h','dphie',1h')
      do k=1,100
       write(22,'(f8.4,i10)') -100.+2.0*(k-0.5),
     >  dphieh(k)/10
      enddo
      write(22,'(''hist'')')


c plots of combined rate ratios versus th_shms
      open(unit=22,file='ptcrat.top')
      write(22,'(''set device postscript'')')
      do ix=1,2
       do iy=1,2
        write(22,201) ix,iy
 201    format(1x,
     >   1x,'set window x ',i1,' of 2 y ',
     >    i1,' of 2'/ 
     >   1x,'set font duplex'/
     >   1x,'set intensity 4'/
     >   1x,'set bar size 0.'/
     >   1x,'set order x y dy SYM'/
     >   1x,'set sym 9O size 1.0'/
     >   1x,'set limits y 0 2.'/
     >   1x,'title left ',1h',' ',1h'/
     >   1x,'title ',1h','ratio',1h'/
     >   1x,'title bottom ',1h','theta (deg)',1h'/
     >   1x,'set ticks size 0.05')
        if(ix.eq.1.and.iy.eq.2) write(22,202)
 202    format(1x,'title top',1h',
     >    'x=0.3 Q2=3.0 z=0.45',1h')
        if(ix.eq.1.and.iy.eq.1) write(22,203)
 203    format(1x,'title top',1h',
     >    'x=0.3 Q2=3.0 z=0.33',1h')
        if(ix.eq.2.and.iy.eq.2) write(22,204)
 204    format(1x,'title top',1h',
     >    'x=0.3 Q2=4.1 z=0.45',1h')
        if(ix.eq.2.and.iy.eq.1) write(22,205)
 205    format(1x,'title top',1h',
     >    'x=0.3 Q2=4.1 z=0.33',1h')
        ikin = 2*(ix-1)+iy
        do itp=2,6
         if(itp.eq.2) write(22,'(''set symbol 9O size 1.2'')')
         if(itp.eq.3) write(22,'(''set symbol 2O size 1.2'')')
         if(itp.eq.4) write(22,'(''set symbol 1O size 1.2'')')
         if(itp.eq.5) write(22,'(''set symbol 9O size 1.2'')')
         if(itp.eq.6) write(22,'(''set symbol 2O size 1.2'')')
         do ith=1,15
          if(avrate(ikin,1,ith).ne.0.) then
           if(avrate(ikin,itp,ith).gt.0.) 
     >      write(22,'(3f10.4)') ith*2. + 4.,
     >      avrate(ikin,itp,ith) /
     >      avrate(ikin,  1,ith),
     >      avrateer(ikin,itp,ith) /
     >      avrate(ikin,  1,ith)
          endif
         enddo
         write(22,'(''plot'')')
         if(itp.ge.5) then
          write(22,'(''set sym size 0.2 ; plot'')')
          write(22,'(''set sym size 0.4 ; plot'')')
          write(22,'(''set sym size 0.6 ; plot'')')
          write(22,'(''set sym size 0.8 ; plot'')')
          write(22,'(''set sym size 1.0 ; plot'')')
         endif
c simc ratios
         do ith=1,15
          if(savrate(ikin,1,ith).ne.0.) then
           fact=1.0
           if(savrate(ikin,itp,ith).gt.0.) 
     >      write(22,'(3f10.4)') ith*2. + 4.,
     >      fact*savrate(ikin,itp,ith) /
     >      savrate(ikin,  1,ith)
          endif
         enddo
         write(22,'(''join 1'')')
c simc ratios
         do ith=1,15
          if(savratenox(ikin,1,ith).ne.0.) then
           fact=1.0
           if(savratenox(ikin,itp,ith).gt.0.) 
     >      write(22,'(3f10.4)') ith*2. + 4.,
     >      fact*savratenox(ikin,itp,ith) /
     >      savratenox(ikin,  1,ith)
          endif
         enddo
         write(22,'(''set pattern 0.05 0.05 0.05 0.05'')')
         write(22,'(''join pattern 1'')')
        enddo
       enddo
      enddo
      close(unit=22)

c plots of on-line combined rate ratios page 2
      open(unit=22,file='ptcrat2.top')
      write(22,'(''set device postscript'')')
      do ix=1,2
       do iy=1,2
        write(22,201) ix,iy
        if(ix.eq.1.and.iy.eq.1) write(22,232)
 232    format(1x,'title top',1h',
     >    'x=0.45 Q2=4.5 z=0.35',1h')
        if(ix.eq.1.and.iy.eq.2) write(22,233)
 233    format(1x,'title top',1h',
     >    'x=0.45 Q2=4.5 z=0.45',1h')
        if(ix.eq.2.and.iy.eq.1) write(22,234)
 234    format(1x,'title top',1h',
     >    'x=0.45 Q2=4.5 z=0.60',1h')
        if(ix.eq.2.and.iy.eq.2) write(22,235)
 235    format(1x,'title top',1h',
     >    'x=0.45 Q2=4.5 z=0.90',1h')
        ikin = 2*(ix-1)+iy + 4
        do itp=2,6
         if(itp.eq.2) write(22,'(''set symbol 9O size 1.2'')')
         if(itp.eq.3) write(22,'(''set symbol 2O size 1.2'')')
         if(itp.eq.4) write(22,'(''set symbol 1O size 1.2'')')
         if(itp.eq.5) write(22,'(''set symbol 9O size 1.2'')')
         if(itp.eq.6) write(22,'(''set symbol 2O size 1.2'')')
         do ith=1,15
          if(avrate(ikin,1,ith).ne.0.0 .and.
     >      (ikin.ne.8.or.ith.ne.7)) then
           th = ith*2. + 4.
           if(ix.eq.1 .and.iy.eq.1) th = 10. + 3.*ith
           if(avrate(ikin,itp,ith).gt.0.) 
     >      write(22,'(3f10.4)') th,
     >      avrate(ikin,itp,ith) /
     >      avrate(ikin,  1,ith),
     >      avrateer(ikin,itp,ith) /
     >      avrate(ikin,  1,ith)
          endif
         enddo
         write(22,'(''plot'')')
         if(itp.ge.5) then
          write(22,'(''set sym size 0.2 ; plot'')')
          write(22,'(''set sym size 0.4 ; plot'')')
          write(22,'(''set sym size 0.6 ; plot'')')
          write(22,'(''set sym size 0.8 ; plot'')')
          write(22,'(''set sym size 1.0 ; plot'')')
         endif
c simc ratios
         do ith=1,15
          if(savrate(ikin,1,ith).ne.0.) then
           th = ith*2. + 4.
           if(ix.eq.1 .and.iy.eq.1) th = 10. + 3.*ith
           fact=1.0
c now done above
c           if(itp.eq.3.or.itp.eq.6) fact=0.75

           if(savrate(ikin,itp,ith).gt.0.) 
     >      write(22,'(3f10.4)') th,
     >      fact * savrate(ikin,itp,ith) /
     >      savrate(ikin,  1,ith)
          endif
         enddo
         write(22,'(''join 1'')')
c simc ratios
         do ith=1,15
          if(savratenox(ikin,1,ith).ne.0.) then
           fact=1.0
           if(savratenox(ikin,itp,ith).gt.0.) 
     >      write(22,'(3f10.4)') ith*2. + 4.,
     >      fact*savratenox(ikin,itp,ith) /
     >      savratenox(ikin,  1,ith)
          endif
         enddo
         write(22,'(''set pattern 0.05 0.05 0.05 0.05'')')
         write(22,'(''join pattern 1'')')
        enddo
       enddo
      enddo
      close(unit=22)

c plots of LH2/SIMC ratios versus th_shms
      open(unit=22,file='ptcrats.top')
      write(22,'(''set device postscript'')')
      do ix=1,2
       do iy=1,2
        write(22,201) ix,iy
        if(ix.eq.1.and.iy.eq.2) write(22,202)
        if(ix.eq.1.and.iy.eq.1) write(22,203)
        if(ix.eq.2.and.iy.eq.2) write(22,204)
        if(ix.eq.2.and.iy.eq.1) write(22,205)
        ikin = 2*(ix-1)+iy
        write(22,'(''set symbol 9O size 1.2'')')
        do ith=1,15
         if(avrate(ikin,1,ith).ne.0.) then
          write(22,'(3f10.4)') ith*2. + 4.,
     >      avrate(ikin,1,ith) /
     >      savrate(ikin,1,ith),
     >      avrateer(ikin,1,ith) /
     >      savrate(ikin,1,ith)
         endif
        enddo
        write(22,'(''plot'')')
        write(22,'(''set sym size 0.2 ; plot'')')
        write(22,'(''set sym size 0.4 ; plot'')')
        write(22,'(''set sym size 0.6 ; plot'')')
        write(22,'(''set sym size 0.8 ; plot'')')
        write(22,'(''set sym size 1.0 ; plot'')')
       enddo
      enddo
      close(unit=22)

c plots of LH2/SIMC ratios versus th_shms p. 2
      open(unit=22,file='ptcrat2s.top')
      write(22,'(''set device postscript'')')
      do ix=1,2
       do iy=1,2
        write(22,201) ix,iy
        if(ix.eq.1.and.iy.eq.2) write(22,232)
        if(ix.eq.1.and.iy.eq.1) write(22,233)
        if(ix.eq.2.and.iy.eq.2) write(22,234)
        if(ix.eq.2.and.iy.eq.1) write(22,235)
        ikin = 2*(ix-1)+iy + 4
        write(22,'(''set symbol 9O size 1.2'')')
        do ith=1,15
         if(avrate(ikin,1,ith).ne.0.) then
          write(22,'(3f10.4)') ith*2. + 4.,
     >      avrate(ikin,1,ith) /
     >      savrate(ikin,1,ith),
     >      avrateer(ikin,1,ith) /
     >      savrate(ikin,1,ith)
         endif
        enddo
        write(22,'(''plot'')')
        write(22,'(''set sym size 0.2 ; plot'')')
        write(22,'(''set sym size 0.4 ; plot'')')
        write(22,'(''set sym size 0.6 ; plot'')')
        write(22,'(''set sym size 0.8 ; plot'')')
        write(22,'(''set sym size 1.0 ; plot'')')
       enddo
      enddo
      close(unit=22)

c hg spectrum for pions versus delta
      open(unit=69,file='ptc.hgcer')
      write(69,'('' dpp hgnphe=1, 2, 3....18'')')
      do k=1,30
       write(69,'(f4.1,25i4)') 2. + .1*(float(k)-0.5)
     >  ,(hghp(k,j)/10,j=1,18)
      enddo

c hg spectrum for pions versus delta
      open(unit=69,file='hgnphepe.txt')
      write(69,'('' dpp hgnphe=1, 2, 3....18'')')
      do k=1,30
       write(69,'(f4.1,25i4)') -5.+(float(k)-0.5)/3
     >  ,(hgh(k,j)/10,j=1,18)
      enddo


      do i=1,4
       do iy=1,3
        do ix=1,100
         sum3=0.
         do k=1,9
          sum1 = (ctx(i,ix,iy,k) * k +
     >      ctx(i,ix,iy,k+1) * (k+1) ) 
          sum2 = (ctx(i,ix,iy,k)  +
     >      ctx(i,ix,iy,k+1) ) 
          if(sum2.gt.sum3) then
           sum3 = sum2
           pipeak = sum1 / sum2
          endif
         enddo
         write(38,'(i1,i2,i4,i7,f5.2,20i4)')
     >    i,iy,ix,int(sum3),pipeak,
     >    (ctx(i,ix,iy,k)/10,k=1,10)
        enddo
       enddo
      enddo
      do i=1,4
       do ix=1,3
        do iy=1,100
         sum3=0.
         do k=1,9
          sum1 = (cty(i,ix,iy,k) * k +
     >      cty(i,ix,iy,k+1) * (k+1) ) 
          sum2 = (cty(i,ix,iy,k)  +
     >      cty(i,ix,iy,k+1) ) 
          if(sum2.gt.sum3) then
           sum3 = sum2
           pipeak = sum1 / sum2
          endif
         enddo
         write(38,'(i1,i2,i4,i7,f5.2,20i4)')
     >    i,iy,ix,int(sum3),pipeak,
     >    (cty(i,ix,iy,k)/10,k=1,10)
        enddo
       enddo
      enddo
      do k=1,5000
c       write(16,'(f8.1,2i10)') 0.2*(k-0.5),ctrf(k,1),ctrf(k,2)
      enddo

      write(15,'(''mc'',20i5)') (mmpi2hs(k),k=1,20)

c find fiducial acceptance region
c phi is yptar (horizontal)
      open(unit=22,file='ptc.accep')
      do ipt=1,30
       do iphi=1,30
        do i=1,4
         sum1=0.
         do ith=1,70
          sum1 = sum1 + accep(i,ipt,ith,iphi)
         enddo
         sum2 = 0.
         do ith=1,70 
          sum2 = sum2 + accep(i,ipt,ith,iphi)
          if(sum2.lt.0.05*sum1) loedge(i) = ith
         enddo
         sum2 = 0.
         do ith=70,1,-1
          sum2 = sum2 + accep(i,ipt,ith,iphi)
          if(sum2.lt.0.05*sum1) hiedge(i) = ith
         enddo
        enddo
        write(22,'(10i4)') ipt,iphi,(loedge(i),i=1,4),
     >   (hiedge(i),i=1,4)
! save simc hms result here
        lo_simc_hms(ipt,iphi)=loedge(2)
        hi_simc_hms(ipt,iphi)=hiedge(2)
       enddo
      enddo

c find which x, xp offset in HMS gives best agreement
      do i=1,20
       totdiff(i)=0
      enddo
      do ipt=6,25
       do iphi=6,25
        do i=11,19
         sum1=0.
         do ith=1,70
          sum1 = sum1 + accep(i,ipt,ith,iphi)
         enddo
         sum2 = 0.
         do ith=1,70 
          sum2 = sum2 + accep(i,ipt,ith,iphi)
          if(sum2.lt.0.05*sum1) loedge(i) = ith
         enddo
         sum2 = 0.
         do ith=70,1,-1
          sum2 = sum2 + accep(i,ipt,ith,iphi)
          if(sum2.lt.0.05*sum1) hiedge(i) = ith
         enddo
         totdiff(i) = totdiff(i) + 
     >     abs(loedge(i) - lo_simc_hms(ipt,iphi)) +
     >     abs(hiedge(i) - hi_simc_hms(ipt,iphi))
c         write(6,'(3i3,4i4)') ipt,iphi,i,
c     >     loedge(i),lo_simc_hms(ipt,iphi), 
c     >     hiedge(i),hi_simc_hms(ipt,iphi)
        enddo
       enddo
      enddo
      do i=11,19
c       write(6,'(''totdiff'',i3,i10)') i,totdiff(i)
      enddo
      open(unit=22,file='ptc.avacc')
      do ipt=4,27
       do ith=1,70
        do i=1,20
         sumi(i)=0.
         do iphi=1,30
          sumi(i) = sumi(i) + accep(i,ipt,ith,iphi)
         enddo
        enddo
        write(22,'(2i3,2i6,9i5)') ipt,ith,
     >   sumi(1)/10,sumi(2)/10,
     >   (sumi(i),i=11,19)
       enddo
      enddo

c aerogel timing. No useful cut here
c      do kk=1,40
c       write(6,'(i2,f7.0,15f4.0)') kk, 
c     >  aeroth(kk,16),(100.*(aeroth(kk,icc)/
c     >  aeroth(kk,16)),icc=1,15)
c      enddo
      sum1=0.
      do kk=21,40
       sum1 = sum1 + aeroth(kk,16)
      enddo
c      do kk=21,40
c       write(6,'(i4,f7.3)') kk,aeroth(kk,16)/sum1
c      enddo

c ratios evtype = 2 to 4 
      do i2=1,6
       write(16,'(''evtype 2 to 4 ratio'',i2)') i2
       do i1=1,15
        write(16,'(f6.1,2i8,f10.3)'), -30. + 4.*(i1-0.5),
     >   ctev(i1,i2,1),ctev(i1,i2,2),
     >   float(ctev(i1,i2,1))/max(1.,float(ctev(i1,i2,2)))
       enddo
      enddo

c ratios data/MC 
      do i=1,7,2
       do ith=1,70
        do ipt=1,30,3
         ip = (ipt+2)/3
         sumacc(ip,1) = 0.
         sumacc(ip,2) = 0.
         do iphi=1,30
          sumacc(ip,1) = sumacc(ip,1) + accep(i,ipt,ith,iphi)
     >       + accep(i,ipt+1,ith,iphi) + accep(i,ipt+2,ith,iphi)
          sumacc(ip,2) = sumacc(ip,2) + accep(i+1,ipt,ith,iphi)
     >       + accep(i+1,ipt+1,ith,iphi) + accep(i+1,ipt+2,ith,iphi)
         enddo
        enddo
        write(16,'(i1,i3,10f7.2)') i,ith,
     >   (sumacc(ip,1)/max(1.,sumacc(ip,2)),ip=1,10)
       enddo
      enddo

      do i=5,7,2
       do ith=10,60
        do iphi=1,30,3
         ip = (iphi+2)/3
         sumacc(ip,1) = 0.
         sumacc(ip,2) = 0.
         do ipt=1,30
          sumacc(ip,1) = sumacc(ip,1) + accep(i,ipt,ith,iphi)
     >       + accep(i,ipt,ith,iphi+1) + accep(i,ipt,ith,iphi+2)
          sumacc(ip,2) = sumacc(ip,2) + accep(i+1,ipt,ith,iphi)
     >       + accep(i+1,ipt,ith,iphi+1) + accep(i+1,ipt,ith,iphi+2)
         enddo
        enddo
        write(16,'(i1,i3,10f7.2)') i,ith,
     >   (sumacc(ip,1)/max(1.,sumacc(ip,2)),ip=1,10)
       enddo
      enddo
 999  write(6,'(i3,i3,a/a)') icase,i,fname,string(1:80)

      do kk=1,20
       write(16,'(i3,8i6)') kk,pxdist(kk),
     >   (pydist(kk,j),j=1,7)
      enddo

      stop
      end

      subroutine getphi(ebeam,elf4,had,phi)
      implicit none
c
c input ebeam   -> the beam energy   
c       elf4    -> scattered electron px,py,pz
c       had     -> hadron px,py,pz
c
c  output  the angle phi in the photon frame (phi)
c          assuming target is a proton

       real*8 qiu4(4),el04(4),elf4(4),tnorm(4),had(4)
       real*8 vangle,vdotm,phi,pitwo,ebeam,y
c
        pitwo=2*acos(-1.0)   ! 2*pi
c
c define the initial  el0 4-momentum
c
        el04(1)=0
        el04(2)=0
        el04(3)=ebeam
        qiu4(1) = el04(1) - elf4(1)
        qiu4(2) = el04(2) - elf4(2)
        qiu4(3) = el04(3) - elf4(3)
        call crossm(qiu4,el04,tnorm)
        phi=vangle(qiu4,el04,qiu4,had)
        y=vdotm(tnorm,had,3)
        if (y.lt.0) phi = pitwo - phi
        return
        end

       subroutine crossm(a,b,c)
       implicit none
       real*8 a(4),b(4),c(4)
       c(1)=a(2)*b(3)-a(3)*b(2)
       c(2)=a(3)*b(1)-a(1)*b(3)
       c(3)=a(1)*b(2)-a(2)*b(1)
       return
       end
c
       real*8 function vdotm(a,b,n)
       implicit none
       real*8 a(4),b(4),s
       integer i,n
       s=0.0
       do i=1,3
         s = s + a(i)*b(i)
       enddo
       if(n.eq.4) s=s-a(n)*b(n)
       vdotm=s
       return
       end
c
       real*8 function vangle(a,b,c,d)
       implicit none
       real*8 a(4),b(4),c(4),d(4),xm,ym,vcos
       real*8 x(4),y(4),pi,vdotm
       pi=acos(-1.0)
       call crossm(a,b,x)
       call crossm(c,d,y)
       xm=vdotm(x,x,3)
       ym=vdotm(y,y,3)
       if(xm.gt.0.0 .and. ym.gt.0.0) then
         vcos=vdotm(x,y,3)/sqrt(xm)/sqrt(ym)
         if(abs(vcos).lt.1.0) then
            vangle=acos(vcos)
         else
c           if(abs(vcos).gt.1.0) write(6,'(1x,''error, vcos'',
c     >        7f8.3)') x,y,vcos
            if(vcos.ge.1.0)  vangle=0
            if(vcos.le.-1.0)  vangle=pi
         endif 
       else
         write(6,'(1x,''xm,ym='',10f8.3)') xm,ym,c,d,y
         vangle=0
       endif
       return
       end
 
      subroutine getcos(e0,z_e,z_p,mp,costh,pt)
      implicit none
      real*8  e0,z_e(4),z_p(4),pe,pp,q2,z_q(4),
     >  z_w(4),pw, pt
      real*8 gamma,beta,cost,ecm,pcm,pl,mp,costh,pq,w

      pe = z_e(4)
      Pp = sqrt(z_p(1)**2 + z_p(2)**2 +z_p(3)**2)

      Z_Q(1)   = -Z_E(1) 
      Z_Q(2)   = -Z_E(2) 
      Z_Q(3)   =  E0 - Z_E(3) 
      Z_Q(4)   =  E0 - Z_E(4)
      q2 = z_e(1)**2 + z_e(2)**2 + (E0 - z_e(3))**2 - 
     >   (E0 - z_e(4))**2
      w = sqrt(0.9383**2 + 2.*0.9383*(E0-z_e(4))-q2)
      PQ       =  SQRT(Z_Q(4)**2+Q2)

      Z_W(1)   = Z_Q(1) 
      Z_W(2)   = Z_Q(2) 
      Z_W(3)   = Z_Q(3)
      Z_W(4)   = SQRT(PQ**2+W**2)
      PW       = PQ

      GAMMA = Z_W(4)/W
      BETA  = PW/Z_W(4)
      COST  =(Z_W(1)*Z_P(1)+Z_W(2)*Z_P(2)+Z_W(3)*Z_P(3))/(PW*PP)
      ECM   = GAMMA*(Z_P(4)-BETA*PP*COST)
      PCM   = SQRT(ECM**2-MP**2)
      PL    = GAMMA*(PP*COST-BETA*Z_P(4))
      pt = sqrt(pcm**2 - pl**2)
      costh = -1.1
      if(pcm.ne.0.) costh = pl/pcm
      return
      end


! improved ytarg for SHMS
      subroutine  newz(del,yfp,ypfp,y)
      implicit none
      real*8 del,yfp,ypfp,y
      real*8 xval(10)/
     >    0.319E+00,
     >   -0.774E+02,
     >    0.273E-01,
     >    0.228E+01,
     >   -0.444E-02,
     >   -0.178E+03,
     >   -0.177E-01,
     >    0.221E+02,
     >    0.255E-02,
     >    0.521E-01/

        y = xval(1) * yfp +
     >      xval(2) * ypfp +
     >      xval(3) * del +
     >      xval(4) * yfp  * ypfp +
     >      xval(5) * yfp **2 +
     >      xval(6) * ypfp **2 +
     >      xval(7) * yfp  * del +
     >      xval(8) * ypfp * del +
     >      xval(9) * del**2  +
     >      xval(10)             
        return
        end

c apply optics correction to SHMS delta
      subroutine shmscorr(pdcx,pdcxp,pdcy,pdcyp,dpp,dppcorr)
      implicit none
      real*8 pdcx,pdcxp,pdcy,pdcyp,dpp,dppcorr,x(4),y,xval(35)
      integer i

      real*8 delcorI(35)/
     >   0.3080E+01,
     >   0.3467E-01,
     >   0.2849E-01,
     >  -0.2840E-01,
     >   0.1959E-01,
     >  -0.5584E-03,
     >   0.2295E-03,
     >  -0.5946E-02,
     >  -0.1155E-03,
     >   0.1136E-02,
     >   0.1675E-02,
     >  -0.1800E-02,
     >  -0.3024E-03,
     >  -0.8011E-03,
     >   0.6904E-04,
     >  -0.5475E-05,
     >   0.1552E-03,
     >   0.8019E-04,
     >  -0.1393E-06,
     >  -0.8338E-04,
     >   0.7756E-05,
     >   0.4899E-05,
     >  -0.6119E-04,
     >  -0.4821E-04,
     >  -0.9406E-05,
     >   0.2010E-03,
     >   0.4101E-05,
     >  -0.1131E-03,
     >   0.1229E-04,
     >   0.1584E-04,
     >  -0.1634E-03,
     >   0.2213E-04,
     >   0.2146E-04,
     >  -0.2847E-04,
     >  -0.1118E-04/
      real*8 delcorH(35)/
     >   0.6903E+00,
     >   0.1539E-01,
     >   0.7762E-02,
     >  -0.1224E-01,
     >   0.1751E-03,
     >  -0.8011E-03,
     >  -0.2192E-03,
     >  -0.3854E-02,
     >   0.2729E-03,
     >   0.3245E-03,
     >   0.1623E-02,
     >  -0.7236E-03,
     >  -0.1781E-03,
     >  -0.3349E-04,
     >  -0.5477E-03,
     >   0.6722E-05,
     >  -0.1099E-03,
     >   0.6942E-04,
     >  -0.7247E-06,
     >  -0.2062E-04,
     >   0.9084E-05,
     >   0.6717E-05,
     >   0.4034E-04,
     >  -0.3140E-04,
     >  -0.9142E-05,
     >   0.3675E-04,
     >   0.1783E-04,
     >   0.1110E-03,
     >   0.4201E-05,
     >  -0.2906E-05,
     >  -0.2211E-04,
     >   0.1502E-04,
     >   0.9743E-05,
     >  -0.4351E-04,
     >  -0.6663E-05/
      real*8 delcorG(35)/
     >   0.1662E+01,
     >  -0.5426E+00,
     >  -0.4593E-01,
     >   0.1601E+00,
     >   0.2593E+00,
     >  -0.4485E-03,
     >   0.6722E-02,
     >  -0.9802E-02,
     >   0.1921E-02,
     >   0.5744E-03,
     >  -0.4498E-02,
     >   0.1154E-01,
     >   0.9370E-04,
     >  -0.3938E-02,
     >  -0.4350E-02,
     >   0.2679E-04,
     >  -0.5061E-03,
     >   0.8402E-04,
     >  -0.8631E-05,
     >  -0.2249E-04,
     >   0.3740E-04,
     >   0.1649E-04,
     >   0.1416E-03,
     >   0.8376E-04,
     >  -0.1529E-04,
     >   0.7915E-04,
     >   0.2418E-04,
     >   0.2875E-03,
     >   0.2703E-05,
     >   0.4173E-05,
     >  -0.9455E-04,
     >  -0.1481E-04,
     >  -0.6484E-04,
     >   0.5366E-04,
     >   0.2219E-04/
      real*8 delcorF(35)/
     >   0.1166E-01,
     >  -0.1073E-01,
     >   0.7001E-02,
     >   0.3576E-02,
     >  -0.1036E-02,
     >  -0.3685E-02,
     >   0.4576E-03,
     >   0.9232E-03,
     >   0.1646E-02,
     >  -0.1507E-03,
     >  -0.2188E-03,
     >  -0.9253E-03,
     >   0.1678E-03,
     >  -0.4011E-04,
     >  -0.1184E-04,
     >   0.1231E-03,
     >  -0.1795E-03,
     >  -0.1573E-03,
     >   0.7723E-06,
     >  -0.2015E-04,
     >   0.7892E-06,
     >  -0.1099E-04,
     >   0.5480E-04,
     >   0.2222E-04,
     >  -0.9220E-05,
     >   0.1281E-04,
     >  -0.8403E-05,
     >   0.1979E-03,
     >   0.3188E-05,
     >  -0.5569E-05,
     >  -0.1660E-05,
     >  -0.4443E-04,
     >  -0.2993E-04,
     >   0.7957E-04,
     >   0.1026E-04/
      real*8 delcorE(35)/
     >   0.9015E+00,
     >  -0.2850E+00,
     >  -0.3250E+00,
     >   0.1784E+00,
     >   0.2237E+00,
     >  -0.5578E-01,
     >   0.2462E-01,
     >   0.3267E-01,
     >   0.1328E-01,
     >   0.4542E-01,
     >  -0.6385E-02,
     >  -0.3683E-01,
     >  -0.3034E-01,
     >  -0.2734E-02,
     >  -0.1668E-01,
     >  -0.2275E-03,
     >  -0.6309E-03,
     >  -0.1844E-03,
     >  -0.8500E-03,
     >   0.9546E-04,
     >  -0.1648E-03,
     >   0.6685E-04,
     >   0.4655E-03,
     >  -0.2727E-03,
     >   0.1882E-03,
     >  -0.1345E-03,
     >  -0.2861E-03,
     >   0.5717E-04,
     >   0.1550E-04,
     >  -0.8202E-04,
     >   0.4235E-04,
     >  -0.5958E-03,
     >   0.3998E-03,
     >   0.6468E-03,
     >   0.9820E-03/
      real*8 delcorD(35)/
     >   0.5146E+02,
     >   0.4719E+01,
     >   0.9932E+00,
     >   0.1218E+00,
     >  -0.6629E+00,
     >  -0.1476E-01,
     >   0.1144E+00,
     >   0.4807E-02,
     >   0.4108E-01,
     >   0.1449E-01,
     >  -0.2482E-01,
     >   0.4377E-01,
     >   0.3198E-02,
     >  -0.2844E-01,
     >  -0.1213E-01,
     >  -0.4178E-04,
     >  -0.4434E-03,
     >  -0.6943E-03,
     >  -0.4957E-03,
     >  -0.3643E-03,
     >  -0.2264E-03,
     >   0.7206E-03,
     >   0.9498E-03,
     >  -0.5239E-03,
     >   0.1280E-03,
     >  -0.2438E-03,
     >  -0.2339E-04,
     >   0.1065E-02,
     >   0.4963E-03,
     >  -0.2865E-03,
     >   0.1783E-03,
     >   0.3073E-03,
     >   0.2590E-03,
     >   0.4419E-03,
     >   0.3631E-03/
      real*8 delcorC(35)/
     >   0.1803E+03,
     >   0.1549E+02,
     >   0.2286E+01,
     >   0.8640E+00,
     >  -0.1914E+01,
     >  -0.1195E+00,
     >   0.4285E+00,
     >   0.2150E-01,
     >   0.1511E+00,
     >   0.1511E+00,
     >  -0.8258E-01,
     >   0.1506E+00,
     >  -0.6568E-01,
     >  -0.1181E+00,
     >  -0.7646E-01,
     >  -0.9379E-03,
     >  -0.2577E-02,
     >  -0.3421E-02,
     >  -0.2506E-02,
     >  -0.1478E-02,
     >  -0.9571E-03,
     >   0.2377E-02,
     >   0.2921E-02,
     >  -0.1253E-02,
     >   0.3808E-03,
     >   0.3541E-03,
     >  -0.6542E-03,
     >   0.6187E-02,
     >   0.1086E-02,
     >  -0.9090E-03,
     >   0.1059E-03,
     >   0.1310E-02,
     >  -0.2179E-04,
     >   0.2405E-02,
     >   0.2701E-02/
      real*8 delcorB(35)/
     >   0.7240E+02,
     >   0.3687E+01,
     >   0.4572E+01,
     >  -0.1353E+01,
     >  -0.7066E+00,
     >   0.1112E+00,
     >   0.1906E-01,
     >   0.1462E+00,
     >   0.5796E-01,
     >   0.5180E-01,
     >  -0.8295E-01,
     >  -0.1349E+00,
     >   0.7098E-01,
     >   0.9002E-03,
     >  -0.8145E-01,
     >   0.3144E-02,
     >  -0.2182E-02,
     >  -0.2191E-02,
     >  -0.1544E-02,
     >  -0.7824E-03,
     >  -0.2311E-03,
     >   0.7003E-02,
     >  -0.8734E-03,
     >  -0.3305E-02,
     >  -0.1321E-02,
     >  -0.6612E-03,
     >   0.1292E-02,
     >   0.7102E-03,
     >   0.9856E-03,
     >   0.6751E-03,
     >   0.5225E-03,
     >   0.2192E-02,
     >  -0.8836E-03,
     >   0.1737E-02,
     >  -0.5982E-03/
      real*8 delcorA(35)/
     >   0.1552E+03,
     >   0.8702E+01,
     >  -0.6154E+01,
     >  -0.3239E+01,
     >   0.2221E+01,
     >  -0.3908E+00,
     >   0.1345E-01,
     >   0.1308E+00,
     >   0.2316E+00,
     >   0.1661E+00,
     >  -0.1678E+00,
     >  -0.2654E+00,
     >   0.4527E-01,
     >   0.2899E-01,
     >  -0.1749E+00,
     >   0.2176E-01,
     >  -0.2916E-02,
     >  -0.1664E-01,
     >  -0.1370E-02,
     >   0.3799E-02,
     >   0.4724E-02,
     >   0.1361E-01,
     >  -0.9711E-03,
     >  -0.5640E-02,
     >  -0.4430E-02,
     >   0.1135E-02,
     >   0.4153E-02,
     >  -0.7566E-03,
     >  -0.3884E-02,
     >   0.1059E-02,
     >  -0.1347E-02,
     >  -0.6841E-03,
     >  -0.4961E-02,
     >  -0.3170E-02,
     >  -0.3723E-02/

      x(1) = pdcx
      x(2) = pdcy
      x(3) = pdcxp *1000.
      x(4) = pdcyp *1000.

      do i=1,35
       if(dpp.lt.-20) xval(i) = delcorA(i)
       if(dpp.ge.-20.and.dpp.lt.-19) xval(i) = delcorB(i)
       if(dpp.ge.-19.and.dpp.lt.-17) xval(i) = delcorC(i)
       if(dpp.ge.-17.and.dpp.lt.-15) xval(i) = delcorD(i)
       if(dpp.ge.-15.and.dpp.lt.-10) xval(i) = delcorE(i)
       if(dpp.ge.-10.and.dpp.lt. 18) xval(i) = delcorF(i)
       if(dpp.ge. 18.and.dpp.lt. 24) xval(i) = delcorG(i)
       if(dpp.ge. 24.and.dpp.lt. 30) xval(i) = delcorH(i)
       if(dpp.ge. 30) xval(i) = delcorI(i)
      enddo

c all 0
      y = 
     >    xval( 1) * x(1)**0 * x(2)**0 * x(3)**0 * x(4)**0 +
c power total 1
     >    xval( 2) * x(1)**1 * x(2)**0 * x(3)**0 * x(4)**0 +
     >    xval( 3) * x(1)**0 * x(2)**1 * x(3)**0 * x(4)**0 +
     >    xval( 4) * x(1)**0 * x(2)**0 * x(3)**1 * x(4)**0 +
     >    xval( 5) * x(1)**0 * x(2)**0 * x(3)**0 * x(4)**1 +
c power total 2
     >    xval( 6) * x(1)**1 * x(2)**1 * x(3)**0 * x(4)**0 +
     >    xval( 7) * x(1)**1 * x(2)**0 * x(3)**1 * x(4)**0 +
     >    xval( 8) * x(1)**1 * x(2)**0 * x(3)**0 * x(4)**1 +
     >    xval( 9) * x(1)**0 * x(2)**1 * x(3)**1 * x(4)**0 +
     >    xval(10) * x(1)**0 * x(2)**1 * x(3)**0 * x(4)**1 +
     >    xval(11) * x(1)**0 * x(2)**0 * x(3)**1 * x(4)**1 +
     >    xval(12) * x(1)**2 * x(2)**0 * x(3)**0 * x(4)**0 +
     >    xval(13) * x(1)**0 * x(2)**2 * x(3)**0 * x(4)**0 +
     >    xval(14) * x(1)**0 * x(2)**0 * x(3)**2 * x(4)**0 +
     >    xval(15) * x(1)**0 * x(2)**0 * x(3)**0 * x(4)**2 +
c power total 3
     >    xval(16) * x(1)**2 * x(2)**1 * x(3)**0 * x(4)**0 +
     >    xval(17) * x(1)**2 * x(2)**0 * x(3)**1 * x(4)**0 +
     >    xval(18) * x(1)**2 * x(2)**0 * x(3)**0 * x(4)**1 +
     >    xval(19) * x(1)**0 * x(2)**2 * x(3)**1 * x(4)**0 +
     >    xval(20) * x(1)**0 * x(2)**2 * x(3)**0 * x(4)**1 +
     >    xval(21) * x(1)**0 * x(2)**0 * x(3)**2 * x(4)**1 +
     >    xval(22) * x(1)**1 * x(2)**2 * x(3)**0 * x(4)**0 +
     >    xval(23) * x(1)**1 * x(2)**0 * x(3)**2 * x(4)**0 +
     >    xval(24) * x(1)**1 * x(2)**0 * x(3)**0 * x(4)**2 +
     >    xval(25) * x(1)**0 * x(2)**1 * x(3)**2 * x(4)**0 +
     >    xval(26) * x(1)**0 * x(2)**1 * x(3)**0 * x(4)**2 +
     >    xval(27) * x(1)**0 * x(2)**0 * x(3)**1 * x(4)**2 +
     >    xval(28) * x(1)**3 * x(2)**0 * x(3)**0 * x(4)**0 +
     >    xval(29) * x(1)**0 * x(2)**3 * x(3)**0 * x(4)**0 +
     >    xval(30) * x(1)**0 * x(2)**0 * x(3)**3 * x(4)**0 +
     >    xval(31) * x(1)**0 * x(2)**0 * x(3)**0 * x(4)**3 +
     >    xval(32) * x(1)**1 * x(2)**1 * x(3)**1 * x(4)**0 +
     >    xval(33) * x(1)**1 * x(2)**1 * x(3)**0 * x(4)**1 +
     >    xval(34) * x(1)**1 * x(2)**0 * x(3)**1 * x(4)**1 +
     >    xval(35) * x(1)**0 * x(2)**1 * x(3)**1 * x(4)**1  

      dppcorr = dpp + y

c      write(6,'(7f8.3)') pdcx,pdcxp,pdcy,pdcyp,dpp,dppcorr,y

      return
      end

	subroutine physics_angles(theta0,phi0,dx,dy,ep,px,py,pz )

!Generate physics angles in lab frame.  Theta is angle from beamline.
!phi is projection on x-y plane (so phi=0 is down, sos=pi/2, hms=3*pi/2.
!
!theta=acos( (cos(theta0)-dy*sin(theta0)*sin(phi0))/ sqrt(1+dx**2+dy**2) )
!phi=atan( (dy*cos(theta0)+sin(theta0)*sin(phi0)) / (sin(theta0)*cos(phi0)+dx) )
!
! Note that these formulae assume phi0=pi/2 or 3*pi/2 (thus the sin(phi0)
! gives a -/+ sign for the HMS/SOS).  Thus, we set the cos(phi0) term to zero.

	real*8 dx,dy		!dx/dy (xptar/yptar) for event.
	real*8 theta0,phi0	!central physics angles of spectrometer.
	real*8 theta,phi	!physics angles for event.
	real*8 r,sinth,costh,sinph,cosph	!intermediate variables.
	real*8 tmp,pi,px,py,pz,ep

	pi=3.141592653589793

	costh = cos(theta0)
	sinth = sin(theta0)
	sinph = sin(phi0)
	cosph = cos(phi0)
	r = sqrt( 1. + dx**2 + dy**2 )

	if (abs(cosph).gt.0.0001) then	!phi not at +/- pi/2
	  write(6,*) 'theta,phi bad'
	  write(6,*) 'phi0=',phi0,'=',phi0*180/pi,'degrees'
	endif

	tmp = (costh - dy * sinth * sinph) / r
	if (abs(tmp).gt.1) write(6,*) 'tmp=',tmp
	theta = acos( (costh - dy * sinth * sinph) / r )
	if (dx.ne.0.0) then
	  phi = atan( (dy*costh + sinth*sinph) / dx )	!gives -90 to 90 deg.
	  if (phi.le.0) phi=phi+pi			!make 0 to 180 deg.
	  if (sinph.lt.0.) phi=phi+pi		!add pi to phi for HMS
	else
	  phi = phi0
	endif
        px =  ep * sin(theta) * cos(phi)
        py =  ep * sin(theta) * sin(phi)
        pz =  ep * cos(theta) 

	return
	end

	subroutine mc_hms_recon (xfp,xpfp,yfp,ypfp,fry,
     >   delta_p,delta_t,delta_phi,y_tgt)
C+______________________________________________________________________________
!
! MC_HMS_RECON : Reconstruct target quantities from tracks.
!		   This subroutine is part of the MC_HMS program.
!
! Right-handed coordinates are assumed: X=down, Z=downstream, Y = (Z cross X)
!
! Inputs are from common block in track_*.inc (except for fry):
!  xs, ys, fry  are in cm.
!  dxdzs, dydzs are unitless slopes (we say "radians" to contrast to "mr").
! Matrix Elements want meters and radians, so have to convert cm-->m.
! Output:
!  delta_p is deviation from central momentum(unitless). Convert to % for outpu+
!  delta_t, delta_phi are in "radians".
!  y_tgt is in meters, convert to cm for output.
!
! Author: D. Potterveld, ANL, 18-Mar-1993
!
! Modification History:
!
! 2-August-1995 (RMM) Hardcoded in an extra index on the coeff, expon, and
!		      n_terms variables to specify HMS. (HMS = 1)
!
!  19-AUG-1993	(DHP) Modified to use COSY INFINITY reconstruction coefficients.
C-______________________________________________________________________________
     
	implicit none

	integer*4 specnum
	parameter (specnum = 1)			!this is the HMS routine

C Argument definitions.

        real*8 xfp,xpfp,yfp,ypfp
	real*8	delta_p,delta_t,delta_phi,y_tgt
	real*8	fry			!vertical position at target (+y=down)

C Cosy reconstruction matrix elements.

	integer*4 max_elements
	parameter (max_elements = 1000)

	real*8		coeff(1,4,max_elements)
	integer*2	expon(1,5,max_elements)
	integer*4	n_terms(1),max_order
	real*8		sum(4),hut(5),term

C Misc. variables.

	integer*4	i,j
	integer*4	chan

	logical		firsttime	/.true./
	character*132	line

C No amnesia, please...

	save

C ============================= Executable Code ================================

C First time through, read in coefficients from data file.

	if (firsttime) then
c	  if (.not.locforunt(chan)) stop 'MC_HMS_RECON: No I/O channels!'
c	  open (unit=chan,status='old',readonly,file='hms/recon_cosy.dat')
          chan=89
	  open (unit=chan,status='old',file='hmsrecon_cosy.dat')

! Skkip past header.

	  line = '!'
	  do while (line(1:1).eq.'!')
	    read (chan,1001) line
	  enddo

! Read in coefficients and exponents.

	  n_terms(specnum) = 0
	  max_order = 0
	  do while (line(1:4).ne.' ---')
	    n_terms(specnum) = n_terms(specnum) + 1
	    if (n_terms(specnum).gt.max_elements)
     >		stop 'WCRECON: too many COSY terms!'
	    read (line,1200) (coeff(specnum,i,n_terms(specnum)),i=1,4),
     >				(expon(specnum,j,n_terms(specnum)),j=1,5)
	    read (chan,1001) line
	    max_order = max(max_order, expon(specnum,1,n_terms(specnum)) + 
     >			expon(specnum,2,n_terms(specnum)) +
     >			expon(specnum,3,n_terms(specnum)) +
     >			expon(specnum,4,n_terms(specnum)) +
     >			expon(specnum,5,n_terms(specnum)))
	  enddo
c	  write(6,*) 'HMS: N_TERMS, MAX_ORDER = ',n_terms(specnum),max_order
	  close (unit=chan)
	  firsttime = .false.
	endif

C Reset COSY sums.

	do i = 1,4
	  sum(i) = 0.
	enddo

C Convert hut quantities to right-handed coordinates, in meters and "radians".
C Make sure hut(5) is non-zero, to avoid taking 0.0**0 (which crashes)
c make this an input now.
c        fry = 0.001
	hut(1) = xfp/100.		!cm --> m
	hut(2) = xpfp			!slope ("radians")
	hut(3) = yfp/100.		!cm --> m
	hut(4) = ypfp			!slope ("radians")
	hut(5) = fry/100.		!vert. position at target(cm-->m)
	if (abs(hut(5)).le.1.e-30) hut(5)=1.e-30

C Compute COSY sums.

	do i = 1,n_terms(specnum)
	  term =  hut(1)**expon(specnum,1,i) * hut(2)**expon(specnum,2,i)
     >		* hut(3)**expon(specnum,3,i) * hut(4)**expon(specnum,4,i)
     >		* hut(5)**expon(specnum,5,i)
	  sum(1) = sum(1) + term*coeff(specnum,1,i)
	  sum(2) = sum(2) + term*coeff(specnum,2,i)
	  sum(3) = sum(3) + term*coeff(specnum,3,i)
	  sum(4) = sum(4) + term*coeff(specnum,4,i)
	enddo
     
C Load output values.

	delta_phi = sum(1)		!slope ("radians")
	y_tgt	  = sum(2)*100.		!m --> cm
	delta_t   = sum(3)		!slope ("radians")
	delta_p   = sum(4)*100.		!percent deviation
     
      return

C ============================ Format Statements ===============================

1001	format(a)
1200	format(1x,4g16.9,1x,5i1)

      END
