      program mpi_xy8
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!     Monte Carlo simulation of site-diluted (2+1)d XY model
!       TEST
!     periodic boundary conditions, sequential update of sites
!     uses KISS05 random number generator
!     serial + parallel via conditional compilation
!     -------------------------------------------------------------------
!     History:
!
!     mpi_xy1:     04 June 2015  first version, based on mpi_hw5
!     mpi_xy2:     05 June 2015  performance optimatizations in wolff
!     mpi_xy3:     09 June 2015  fixed bug in error of xis and xit and missing initialization of en2
!     mpi_xy4:     15 June 2015  hot and cold starts
!     mpi_xy5:     15 June 2015  split NMESS for improved estimators
!     mpi_xy6:     18 June 2015  also measure dm/dT directly
!     mpi_xy7:     26 Jan  2016  loop over Lt
!     mpi_xy8:     03 June 2016	 fix missing initialization of conf2xis and conf2xit, stop recalculation of cos and sin
!                  23 May  2019  change impurity loop (294) to classical physics vs quantum physics.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Preprocessor directives
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!#define PARALLEL
#define VERSION 'mpi_xy8'

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Data types
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      implicit none
      integer,parameter         :: r8b= SELECTED_REAL_KIND(P=14,R=99)    ! 8-byte reals !precision = 14 decimals, exponent range = 99
      integer,parameter         :: i4b= SELECTED_INT_KIND(8)             ! 4-byte integers !-10^8 to 10^8
      integer,parameter         :: ilog=kind(.true.)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Simulation parameters
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      integer(i4b),parameter    :: L=32                                ! L=linear system size
      integer(i4b),parameter    :: NLT=1                                ! number of L
      integer(i4b), parameter   :: LTARRAY(NLT)=(/32/) ! system sizes
      integer(i4b),parameter    :: LTMAX=LTARRAY(NLT), L3MAX=L*L*LTMAX  !L3MAX = 32*32*32 = 32768
      real(r8b),   parameter    :: TMAX=2.5D0, TMIN=1.0D0             ! max and min temperatures
      real(r8b),   parameter    :: DT0=0.05                              ! temp step, must be positive
      integer(i4b), parameter   :: COLDSTART = -1                       ! set to 1 for cold start and to -1 for hot start

      real(r8b),   parameter    :: IMPCONC=0.000000D0
      integer(i4b),parameter    :: NCONF=100                         ! number of disorder configs
      logical(ilog), parameter  :: CANON_DIS=.false.

      integer(i4b),parameter    :: NEQ=100,NMESS=500                    ! Monte Carlo sweeps, must be even

      integer(i4b),parameter    :: IRINIT=1                             ! LFSR seed


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Internal constants - do not touch !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      real(r8b), parameter      :: pi=3.141592653589793D0
      integer(i4b),parameter    :: NTEMP= 1+NINT((TMAX-TMIN)/DT0)        ! number of temperatures
      integer(i4b),parameter    :: TDSIZE=17                            ! size of MPI data transfer array


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Variable declarations
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      real (r8b)      :: sx(0:L3MAX-1),sy(0:L3MAX-1)             ! XY spin
      logical(ilog)   :: occu(0:L3MAX-1)                      ! occupation of site with spin
      real (r8b)      :: sxnew, synew, slen
      real (r8b)      :: nsx,nsy, dE                   ! sum over neighboring spins

      real(r8b)       :: mx,my,sweepmag,sweepen                  ! magnetization vector
      real(r8b)       :: sweepmagqt,sweepmagq2, sweepmagq3,magqt,magqs     ! mag(qtime),mag(qspace),
      real(r8b)       :: Gt,Gs,Gtcon,Gscon                         ! G(qspace), G(qtime), connected versions
      real(r8b)       :: mag, mag2, mag4, bin, susc,en, en2,sph     ! magnetization, its square, energy
      real(r8b)       :: enmag, dmdT                                ! <e m>
      real(r8b)       :: mag1half,mag2half,en1half,en2half
      real(r8b)       :: xit,xis,xitcon,xiscon                      ! correlation lengths in space an time, connected versions
      real(r8b)       :: glxit,glxis,glxitcon,glxiscon              ! global correlation lengths in space an time, connected versions

      real(r8b)       :: confmag(NTEMP),confmag2(NTEMP)
      real(r8b)       :: conf2mag(NTEMP), confmag4(NTEMP)              ! configuration averages
      real(r8b)       :: conflogmag(NTEMP)
      real(r8b)       :: confsusc(NTEMP), confbin(NTEMP)
      real(r8b)       :: conf2susc(NTEMP), conf2bin(NTEMP)              ! conf. av. of squared observables
      real(r8b)       :: confen(NTEMP),confsph(NTEMP)
      real(r8b)       :: conf2en(NTEMP),conf2sph(NTEMP)
      real(r8b)       :: confGt(NTEMP),confGs(NTEMP),confGtcon(NTEMP),confGscon(NTEMP)
      real(r8b)       :: confxit(NTEMP),confxis(NTEMP)
      real(r8b)       :: conf2xit(NTEMP),conf2xis(NTEMP)
      real(r8b)       :: confxitcon(NTEMP),confxiscon(NTEMP)
      real(r8b)       :: confdmdT(NTEMP),conf2dmdT(NTEMP)
      real(r8b)       :: confdlnmdT(NTEMP),conf2dlnmdT(NTEMP)

      integer(i4b)    :: m1(0:L3MAX-1)            ! neighbor table
      integer(i4b)    :: m2(0:L3MAX-1)
      integer(i4b)    :: m3(0:L3MAX-1)
      integer(i4b)    :: m4(0:L3MAX-1)
      integer(i4b)    :: m5(0:L3MAX-1)
      integer(i4b)    :: m6(0:L3MAX-1)

      integer(i4b)    :: Lt, iLT, L3                   ! LT, counter, volume

      real(r8b)       :: qspace,qtime          ! minimum q values for correlation length
	  real(r8b)		  :: cosspace(0:L-1)
	  real(r8b)		  :: sinspace(0:L-1)
	  real(r8b)		  :: costime(0:L3MAX-1)
	  real(r8b)		  :: sintime(0:L3MAX-1)

      integer(i4b)    :: iconf,init            ! current disorder config
      integer(i4b)    :: totconf               ! total number of confs run so far
      integer(i4b)    :: N_IMPSITE,iimp        ! number of impurity sites

      integer(i4b)    :: i1,i2,i3             ! coordinates
      integer(i4b)    :: is                 ! site index
      integer(i4b)    :: ispace,itime

      integer(i4b)    :: isweep                ! Monte Carlo sweep
      real(r8b)       :: T, dT, beta               ! Monte Carlo temperature
      integer(i4b)    :: itemp

      integer(i4b)    :: nclsweep,nspsweep        ! number of clusters, spins flipped in single sweep
      real(r8b)       :: totncl,totnsp            ! total numbers of clusters, spins flipped
      integer(i4b)    :: ncluster
      real(r8b)       :: avclsize

      real(r8b),external         :: rkiss05
      external                      kissinit

      character avenfile*15,avmafile*15,avcofile*15,avdmfile*15

! Now the MPI stuff !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#ifdef PARALLEL
      include 'mpif.h'
      integer(i4b)              :: ierr
      integer(i4b)              :: id,myid                  ! process index
      integer(i4b)              :: numprocs              ! total number of processes
      integer(i4b)              :: status(MPI_STATUS_SIZE)
      real(r8b)                 :: transdata(TDSIZE),auxdata(TDSIZE)     ! MPI transfer array
#endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Start of main program
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


! Set up MPI !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#ifdef PARALLEL
      call MPI_INIT(ierr)
      call MPI_COMM_RANK( MPI_COMM_WORLD, myid, ierr )
      call MPI_COMM_SIZE( MPI_COMM_WORLD, numprocs, ierr )

      if (myid==0) then
         print *,'Program ',VERSION,' running on', numprocs, ' processes'
         print *,'--------------------------------------------------'
         print *,'L= ', L
         print *,'MC steps: ', NEQ, ' + ', NMESS
      endif ! of if (myid==0)
#else
      print *,'Program ',VERSION,' running on single processor'
      print *,'--------------------------------------------------'
      print *,'L= ', L
      print *,'MC steps: ', NEQ, ' + ', NMESS
#endif

	  qspace=2*pi/L
	  do i1=0, L-1
		cosspace(i1)=cos(qspace*i1) !i1 = x
		sinspace(i1)=sin(qspace*i1) !i1 = x
	  enddo

! Loop over Lt !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   LT_loop: do iLT=1,NLT
      Lt = LTARRAY(iLT)
      L3 = L*L*Lt

      qtime=2*pi/Lt                 ! minimum q values for correlation length
	  do i1=0, Lt-1
		costime(i1)=cos(qtime*i1)
		sintime(i1)=sin(qtime*i1)
	  enddo

      avenfile='aven0000000.dat'
      avmafile='avma0000000.dat'
      avcofile='avco0000000.dat'
      avdmfile='avdm0000000.dat'
      write(avenfile(5:7),'(I3.3)') L
      write(avmafile(5:7),'(I3.3)') L
      write(avcofile(5:7),'(I3.3)') L
      write(avdmfile(5:7),'(I3.3)') L
      write(avenfile(8:11),'(I4.4)') Lt
      write(avmafile(8:11),'(I4.4)') Lt
      write(avcofile(8:11),'(I4.4)') Lt
      write(avdmfile(8:11),'(I4.4)') Lt

! Set up neigbor table

    do i1=0, Lt-1
     do i2=0, L-1
      do i3=0, L-1
          is = L*(L*i1 + i2) + i3

          if (i1.eq.Lt-1) then
              m1(is)=is - L*L*(Lt-1)
          else
              m1(is)=is + L*L
          endif

          if (i1.eq.0) then
              m2(is)=is + L*L*(Lt-1)
          else
              m2(is)=is - L*L
          endif

          if (i2.eq.L-1) then
              m3(is)= is - L*(L-1)
          else
              m3(is)= is + L
          endif

          if (i2.eq.0) then
              m4(is)=is + L*(L-1)
          else
              m4(is)=is - L
          endif

          if (i3.eq.L-1) then
              m5(is)= is - (L-1)
          else
              m5(is)= is + 1
          endif

          if (i3.eq.0) then
              m6(is)= is + (L-1)
          else
              m6(is)= is - 1
          endif
      enddo
      enddo
      enddo


      confmag(:)   = 0.D0
      conf2mag(:)  = 0.D0
      confmag2(:)  = 0.D0
      confmag4(:)  = 0.D0
      conflogmag(:)= 0.D0
      confsusc(:)  = 0.D0
      confbin(:)   = 0.D0
      conf2susc(:) = 0.D0
      conf2bin(:)  = 0.D0
      confen(:)    = 0.D0
      confsph(:)   = 0.D0
      conf2en(:)   = 0.D0
      conf2sph(:)  = 0.D0
      confGt(:)    = 0.D0
      confGs(:)    = 0.D0
      confGtcon(:) = 0.D0
      confGscon(:) = 0.D0
      confxit(:)   = 0.D0
      confxis(:)   = 0.D0
	   conf2xit(:)  = 0.D0 ! fix v8
      conf2xis(:)  = 0.D0 ! fix v8
      confxitcon(:)= 0.D0
      confxiscon(:)= 0.D0
      confdmdT(:)  = 0.D0
      conf2dmdT(:) = 0.D0
      confdlnmdT(:)= 0.D0
      conf2dlnmdT(:)=0.D0


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  Loop over disorder configurations
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#ifdef PARALLEL
      disorder_loop: do iconf=myid+1,(NCONF/numprocs)*numprocs,numprocs
         if (myid==0) print *, 'dis. conf.', iconf
#else
      disorder_loop: do iconf=1,NCONF
         print *, 'dis. conf.', iconf
#endif

      N_IMPSITE=nint(L*L*impconc)

! Initialize random number generator, IRINIT must be positive
      init=irinit+iconf-1
      call kissinit(init)

! Initialize impurities
      occu(:)=.true. !no impurities
      iimp=0 !for first loop
      IF (CANON_DIS) THEN
      do while (iimp.lt.N_IMPSITE) !number of impurities < num_total_impurities
         ispace = int(rkiss05()*L3MAX)
         if (occu(ispace)) then
             occu(ispace) = .false.
!             do itime=0,Lt-1
!                occu(ispace+L*L*itime)=.false.
!                print *,ispace+L*L*itime
!             enddo
            iimp=iimp+1
         endif
      enddo
      ELSE
      do ispace=0,L3MAX-1
         if (rkiss05()<IMPCONC) then
            occu(ispace) = .false.
            ! do itime=0,Lt-1
            !    occu(ispace+L*L*itime)=.false.
            ! enddo
            iimp=iimp+1
         endif
      enddo
      ENDIF


! Initialize spins
      if (COLDSTART==1) then
        do is=0, L3-1
        if (occu(is)) then
           sx(is)=1.D0
           sy(is)=0.D0
        endif
        enddo
      else
        do is=0, L3-1
        if (occu(is)) then
          slen=1.D0
          do while (slen>0.5D0)
             sx(is)= rkiss05()-0.5D0
             sy(is)= rkiss05()-0.5D0
             slen=sqrt(sx(is)*sx(is) + sy(is)*sy(is) )
          enddo
          slen=1.D0/slen		! changed v8
          sx(is)=sx(is)*slen
          sy(is)=sy(is)*slen
        endif
        enddo
      endif

! Loop over temperatures
         if (COLDSTART==1) then
            T=TMIN
            dT=DT0
         else
            T=TMAX
            dT=-DT0
         endif
      temperature: do itemp=1,NTEMP
      beta=1.D0/T
!      print *,'T= ',real(T)

!     Equilibration, carry out NEQ full MC steps

      do isweep=1,NEQ/2
         call metro_sweep
         call wolff_sweep(0,nclsweep,nspsweep)
      enddo
      totncl=0.D0
      totnsp=0.D0
      do isweep=1, NEQ/2
         call metro_sweep
         call wolff_sweep(0,nclsweep,nspsweep)
         totncl=totncl+nclsweep
         totnsp=totnsp+nspsweep
      enddo
      avclsize=totnsp/totncl
      ncluster=L3/avclsize+1

! Measuring loop, carry out NMESS full MC sweeps

      mag= 0.D0
      mag2=0.D0
      mag4=0.D0
      en=  0.D0
      en2= 0.D0
      magqt=0.D0
      magqs=0.D0
      Gt=  0.D0
      Gs=  0.D0
      enmag=0.D0

      totncl=0.D0
      totnsp=0.D0
      do isweep=1, NMESS/2
         call metro_sweep
         call wolff_sweep(ncluster,nclsweep,nspsweep)
         totncl=totncl+nclsweep
         totnsp=totnsp+nspsweep

         call measurement
         call corr_func
      enddo     ! of do isweep ...
      mag1half=mag
      en1half=en

      call metro_sweep                                    ! separate the two halfs
      call wolff_sweep(ncluster,nclsweep,nspsweep)

      do isweep=1, NMESS/2
         call metro_sweep
         call wolff_sweep(ncluster,nclsweep,nspsweep)
         totncl=totncl+nclsweep
         totnsp=totnsp+nspsweep

         call measurement
         call corr_func
      enddo     ! of do isweep ...
      mag2half=mag-mag1half
      en2half=en-en1half
      avclsize=totnsp/totncl


      mag=mag/NMESS
      mag2=mag2/NMESS
      mag4=mag4/NMESS
      en=en/NMESS
      en2=en2/NMESS
      mag1half=2.D0*mag1half/NMESS
      mag2half=2.D0*mag2half/NMESS
      en1half=2.D0*en1half/NMESS
      en2half=2.D0*en2half/NMESS
      enmag=enmag/NMESS
      magqt=magqt/NMESS
      magqs=magqs/NMESS
      Gt=Gt/NMESS
      Gs=Gs/NMESS

      susc=(mag2 - mag1half*mag2half)*L3*beta
      bin=1-mag4/(3*mag2**2)
      sph=(en2-en1half*en2half)*L3*beta**2
      dmdT=(enmag-0.5*(en1half*mag2half+en2half*mag1half))*L3*beta**2

      Gtcon=Gt-magqt**2
      Gscon=Gs-magqs**2


      xit= (mag2 - Gt)/ (Gt*qtime*qtime)
      xit=sqrt(abs(xit))
      xis= (mag2 - Gs)/ (Gs*qspace*qspace)
      xis=sqrt(abs(xis))
      xitcon= ((mag2 - mag**2) - Gtcon)/ (Gtcon*qtime*qtime)
      xitcon=sqrt(abs(xitcon))
      xiscon= ((mag2 - mag**2) - Gscon)/ (Gscon*qspace*qspace)
      xiscon=sqrt(abs(xiscon))


#ifdef PARALLEL
!! Package data for transmission
      transdata(1)=mag
      transdata(2)=mag2
      transdata(3)=mag4
      transdata(4)=en
      transdata(5)=en2
      transdata(6)=susc
      transdata(7)=bin
      transdata(8)=sph
      transdata(9)=Gt
      transdata(10)=Gs
      transdata(11)=Gtcon
      transdata(12)=Gscon
      transdata(13)=xit
      transdata(14)=xis
      transdata(15)=xitcon
      transdata(16)=xiscon
      transdata(17)=dmdT


      if(myid.ne.0) then                                                  ! Send data
         call MPI_SEND(transdata,TDSIZE,MPI_DOUBLE_PRECISION,0,myid,MPI_COMM_WORLD,ierr)
      else                                                                 ! Receive data
         do id=0,numprocs-1
            if (id==0) then
               auxdata(:)=transdata(:)
            else
               call MPI_RECV(auxdata,TDSIZE,MPI_DOUBLE_PRECISION,id,id,MPI_COMM_WORLD,status,ierr)
            endif
            confmag(itemp)=confmag(itemp)       +auxdata(1)
            conf2mag(itemp)=conf2mag(itemp)     +(auxdata(1))**2
            confmag2(itemp)=confmag2(itemp)     +auxdata(2)
            confmag4(itemp)=confmag4(itemp)     +auxdata(3)
            conflogmag(itemp)=conflogmag(itemp) +log(auxdata(1))
            confsusc(itemp)=confsusc(itemp)     +auxdata(6)
            confbin(itemp)=confbin(itemp)       +auxdata(7)
            conf2susc(itemp)=conf2susc(itemp)   +(auxdata(6))**2
            conf2bin(itemp)=conf2bin(itemp)     +(auxdata(7))**2
            confen(itemp)=confen(itemp)         +auxdata(4)
            confsph(itemp)=confsph(itemp)       +auxdata(8)
            conf2en(itemp)=conf2en(itemp)       +(auxdata(4))**2
            conf2sph(itemp)=conf2sph(itemp)     +(auxdata(8))**2
            confGt(itemp)=confGt(itemp)         +auxdata(9)
            confGs(itemp)=confGs(itemp)         +auxdata(10)
            confGtcon(itemp)=confGtcon(itemp)   +auxdata(11)
            confGscon(itemp)=confGscon(itemp)   +auxdata(12)
            confxit(itemp)=confxit(itemp)       +auxdata(13)
            confxis(itemp)=confxis(itemp)       +auxdata(14)
            conf2xit(itemp)=conf2xit(itemp)     +(auxdata(13))**2
            conf2xis(itemp)=conf2xis(itemp)     +(auxdata(14))**2
            confxitcon(itemp)=confxitcon(itemp) +auxdata(15)
            confxiscon(itemp)=confxiscon(itemp) +auxdata(16)
            confdmdT(itemp)=confdmdT(itemp)     +auxdata(17)
            conf2dmdT(itemp)=conf2dmdT(itemp)   +(auxdata(17))**2
            confdlnmdT(itemp)=confdlnmdT(itemp) +auxdata(17)/auxdata(1)
            conf2dlnmdT(itemp)=conf2dlnmdT(itemp)+(auxdata(17)/auxdata(1))**2

         enddo
      endif

#else
      confmag(itemp)=confmag(itemp)      +mag
      conf2mag(itemp)=conf2mag(itemp)    +mag**2
      confmag2(itemp)=confmag2(itemp)    +mag2
      confmag4(itemp)=confmag4(itemp)    +mag4
      conflogmag(itemp)=conflogmag(itemp)+log(mag)
      confsusc(itemp)=confsusc(itemp)    +susc
      confbin(itemp)=confbin(itemp)      +bin
      conf2susc(itemp)=conf2susc(itemp)  +susc**2
      conf2bin(itemp)=conf2bin(itemp)    +bin**2
      confen(itemp)=confen(itemp)        +en
      confsph(itemp)=confsph(itemp)      +sph
      conf2en(itemp)=conf2en(itemp)      +en**2
      conf2sph(itemp)=conf2sph(itemp)    +sph**2
      confGt(itemp)=confGt(itemp)        +Gt
      confGs(itemp)=confGs(itemp)        +Gs
      confGtcon(itemp)=confGtcon(itemp)  +Gtcon
      confGscon(itemp)=confGscon(itemp)  +Gscon
      confxit(itemp)=confxit(itemp)       +xit
      confxis(itemp)=confxis(itemp)       +xis
      conf2xit(itemp)=conf2xit(itemp)     +xit**2
      conf2xis(itemp)=conf2xis(itemp)     +xis**2
      confxitcon(itemp)=confxitcon(itemp) +xitcon
      confxiscon(itemp)=confxiscon(itemp) +xiscon
      confdmdT(itemp)=confdmdT(itemp)     +dmdT
      conf2dmdT(itemp)=conf2dmdT(itemp)   +dmdT**2
      confdlnmdT(itemp)=confdlnmdT(itemp) +dmdT/mag
      conf2dlnmdT(itemp)=conf2dlnmdT(itemp)+(dmdT/mag)**2

#endif

      T=T+dT
      enddo temperature

!!! output !!!!!!!!!!!!!!!!!!!!!!!!

#ifdef PARALLEL
      totconf=iconf+numprocs-1
      if (myid==0) then
#else
      totconf=iconf
#endif


      open(7,file=avenfile,status='replace')
      rewind(7)
      write(7,*) 'program ', VERSION
      write(7,*) 'spatial + temporal system size', L, Lt
      write(7,*) 'equilibration steps     ',  NEQ
      write(7,*) 'measurement steps   ', NMESS
      write(7,*) 'impurity conc.  ', real(impconc),',  number of imps. ',N_IMPSITE
      write(7,*) 'disorder configurations ', NCONF
      write(7,*) 'disorder configurations processed ',totconf,' of  ',NCONF
      write(7,*) 'CANON_DIS ',CANON_DIS
      write(7,*) 'COLDSTART ',COLDSTART
      write(7,*) 'LFSR-Init        ', IRINIT
      write(7,*) '-----------------'
      write(7,*) '  T    [<energy>]  std.en/sqrt(totconf)    [<spec.heat>]    std.spec.heat/sqrt(totconf)'
      if (COLDSTART==1) then
         T=TMIN
         dT=DT0
      else
         T=TMAX
         dT=-DT0
      endif
      do itemp=1,NTEMP
        write(7,'(1x,25(e13.7,1x))') t,confen(itemp)/totconf,&
                                     (sqrt(conf2en(itemp)/totconf-(confen(itemp)/totconf)**2))/sqrt(1.D0*totconf),&
                                     confsph(itemp)/totconf,&
                                     (sqrt(conf2sph(itemp)/totconf-(confsph(itemp)/totconf)**2))/sqrt(1.D0*totconf)
        T=T+dT
      enddo
      close(7)
      open(7,file=avmafile,status='replace')
      rewind(7)
      write(7,*) 'program ', VERSION
      write(7,*) 'spatial + temporal system size', L, Lt
      write(7,*) 'equilibration steps     ',  NEQ
      write(7,*) 'measurement steps   ', NMESS
      write(7,*) 'impurity conc.  ', real(impconc),',  number of imps. ',N_IMPSITE
      write(7,*) 'disorder configurations ', NCONF
      write(7,*) 'disorder configurations processed ',totconf,' of  ',NCONF
      write(7,*) 'CANON_DIS ',CANON_DIS
      write(7,*) 'COLDSTART ',COLDSTART
      write(7,*) 'LFSR-Init        ', IRINIT
      write(7,*) '-----------------'
      write(7,*) ' T    [<mag>]  [<mag>^2]  std.mag/sqrt(totconf)  [<susc>] std.susc/sqrt(iconf)  [<Binder>] ',&
              'std.bin/sqrt(totconf) [log <mag>]  Global.Binder'
      if (COLDSTART==1) then
         T=TMIN
         dT=DT0
      else
         T=TMAX
         dT=-DT0
      endif
      do itemp=1,NTEMP
         write(7,'(1x,25(e13.7,1x))') t,confmag(itemp)/totconf,conf2mag(itemp)/totconf,&
              (sqrt(conf2mag(itemp)/totconf-(confmag(itemp)/totconf)**2))/sqrt(1.D0*totconf),confsusc(itemp)/totconf,&
              (sqrt(conf2susc(itemp)/totconf-(confsusc(itemp)/totconf)**2))/sqrt(1.D0*totconf), confbin(itemp)/totconf,&
              (sqrt(conf2bin(itemp)/totconf-(confbin(itemp)/totconf)**2))/sqrt(1.D0*totconf), conflogmag(itemp)/totconf,&
              1-confmag4(itemp)/totconf/(3*(confmag2(itemp)/totconf)**2)
         T=T+dT
      enddo
      close(7)
      open(7,file=avcofile,status='replace')
      rewind(7)
      write(7,*) 'program ', VERSION
      write(7,*) 'spatial + temporal system size', L, Lt
      write(7,*) 'equilibration steps     ',  NEQ
      write(7,*) 'measurement steps   ', NMESS
      write(7,*) 'impurity conc.  ', real(impconc),',  number of imps. ',N_IMPSITE
      write(7,*) 'disorder configurations ', NCONF
      write(7,*) 'disorder configurations processed ',totconf,' of  ',NCONF
      write(7,*) 'CANON_DIS ',CANON_DIS
      write(7,*) 'COLDSTART ',COLDSTART
      write(7,*) 'LFSR-Init        ', IRINIT
      write(7,*) '-----------------'
      write(7,*) ' T  [xit]  [xis] [xit]/Lt  [xis]/L [xitcon]  [xiscon] [xitcon]/Lt  [xiscon]/L', &
				' glxit glxis  glxit/Lt glxis/L glxitcon glxiscon glxitcon/Lt glxiscon/L',&
				' std.dev.(xit/Lt)  std.dev.(xis/L)'
      if (COLDSTART==1) then
         T=TMIN
         dT=DT0
      else
         T=TMAX
         dT=-DT0
      endif
      do itemp=1,NTEMP
         glxit= (confmag2(itemp) - confGt(itemp))/ (confGt(itemp)*qtime*qtime)
         glxit=sqrt(abs(glxit))
         glxis= (confmag2(itemp) - confGs(itemp))/ (confGs(itemp)*qspace*qspace)
         glxis=sqrt(abs(glxis))
         glxitcon= ((confmag2(itemp) - (confmag(itemp)**2)/totconf) - confGtcon(itemp))/ (confGtcon(itemp)*qtime*qtime)
         glxitcon=sqrt(abs(glxitcon))
         glxiscon= ((confmag2(itemp) - (confmag(itemp)**2)/totconf) - confGscon(itemp))/ (confGscon(itemp)*qspace*qspace)
         glxiscon=sqrt(abs(glxiscon))

         write(7,'(1x,25(e13.7,1x))') t,confxit(itemp)/totconf, confxis(itemp)/totconf,&
              confxit(itemp)/totconf/Lt, confxis(itemp)/totconf/L,&
              confxitcon(itemp)/totconf,confxiscon(itemp)/totconf,confxitcon(itemp)/totconf/Lt,confxiscon(itemp)/totconf/L,&
              glxit,glxis,glxit/Lt,glxis/L,glxitcon,glxiscon,glxitcon/Lt,glxiscon/L,&
              (sqrt(conf2xit(itemp)/totconf-(confxit(itemp)/totconf)**2))/sqrt(1.D0*totconf)/Lt,&
              (sqrt(conf2xis(itemp)/totconf-(confxis(itemp)/totconf)**2))/sqrt(1.D0*totconf)/L
         T=T+DT
      enddo
      close(7)
      open(7,file=avdmfile,status='replace')
      rewind(7)
      write(7,*) 'program ', VERSION
      write(7,*) 'spatial + temporal system size', L, Lt
      write(7,*) 'equilibration steps     ',  NEQ
      write(7,*) 'measurement steps   ', NMESS
      write(7,*) 'impurity conc.  ', real(impconc),',  number of imps. ',N_IMPSITE
      write(7,*) 'disorder configurations ', NCONF
      write(7,*) 'disorder configurations processed ',totconf,' of  ',NCONF
      write(7,*) 'CANON_DIS ',CANON_DIS
      write(7,*) 'COLDSTART ',COLDSTART
      write(7,*) 'LFSR-Init        ', IRINIT
      write(7,*) '-----------------'
      write(7,*) '  T    [<m>]  std.m/sqrt(totconf)   [<dmdT>]   std.dmdT/sqrt(totconf)   [<dlnmdT>]'//&
                                                                'std.dlnmdT/sqrt(totconf)  [<dmdT>]/[<m>]'
      if (COLDSTART==1) then
         T=TMIN
         dT=DT0
      else
         T=TMAX
         dT=-DT0
      endif
      do itemp=1,NTEMP
        write(7,'(1x,25(e13.7,1x))') t,confmag(itemp)/totconf,&
              (sqrt(conf2mag(itemp)/totconf-(confmag(itemp)/totconf)**2))/sqrt(1.D0*totconf), confdmdT(itemp)/totconf,&
              (sqrt(conf2dmdT(itemp)/totconf-(confdmdT(itemp)/totconf)**2))/sqrt(1.D0*totconf),confdlnmdT(itemp)/totconf,&
              (sqrt(conf2dlnmdT(itemp)/totconf-(confdlnmdT(itemp)/totconf)**2))/sqrt(1.D0*totconf),&
              confdmdT(itemp)/confmag(itemp)
        T=T+dT
      enddo
      close(7)

#ifdef PARALLEL
      endif ! of if (myid==0)
#endif

      enddo disorder_loop

      enddo LT_loop

#ifdef PARALLEL
      call MPI_FINALIZE(ierr)
#endif

      stop
      contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine Wolff_sweep(ntobeflipped,nclflipped,nspflipped)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!     Performs a Wolff sweep consisting of ncluster single cluster flips
!
!     If ntobeflipped <= 0, subroutine determines length of sweep by counting
!     flipped spins rather than clusters (usefull for equilibration when
!     cluster size not known)
!
!     Do NOT use ntobeflipped <= 0 for measurement because the last cluster
!     flipped to reach L3 flipped spins is biased towards large clusters
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      integer(i4b)    :: ntobeflipped          ! INPUT: number of spins to be flipped
      integer(i4b)    :: nclflipped            ! OUTPUT: number of clusters flipped
      integer(i4b)    :: nspflipped            ! OUTPUT: number of spins flipped

      integer(i4b)    :: stack(0:L3-1)         ! stack for cluster construction
      integer(i4b)    :: sp                    ! stackpointer
      integer(i4b)    :: current,neighbor      ! current site in cluster construction, its neighbor
      integer(i4b)    :: isize                 ! size of current cluster
      real(r8b)       :: nx,ny                 ! reflection direction
      real(r8b)       :: scalar1,scalar2, snsn          ! scalar products in addition probability
      real(r8b)       :: padd

      nclflipped=0
      nspflipped=0


      cluster_loop: do

         is = int(L3*rkiss05())             ! seed site for Wolff cluster
         if (occu(is)) then                  ! is site occupied?

         stack(0)=is
         sp=1

         slen=1.D0
         do while (slen>0.5D0)
            nx= rkiss05()-0.5D0
            ny= rkiss05()-0.5D0
            slen= sqrt(nx*nx+ny*ny)
         enddo
         slen=1.D0/slen		! changed v8
         nx= nx*slen
         ny= ny*slen

         nclflipped=nclflipped+1
         isize=1
         scalar2= 2.D0*(nx*sx(is)+ny*sy(is))               ! scalar product for p_add

         sx(is)=sx(is)-nx*scalar2                               ! now flip seed spin
         sy(is)=sy(is)-ny*scalar2

         do while(sp.gt.0)                   ! now build the cluster
           sp=sp-1
           current = stack(sp)              ! get site from stack

               scalar1=- (nx*sx(current)+ny*sy(current))                ! scalar product for p_add

               neighbor=m6(current)
               if (occu(neighbor)) then
               scalar2=  2.D0*(nx*sx(neighbor)+ny*sy(neighbor) )   ! scalar product for p_add
               snsn=scalar1*scalar2
               if (snsn>0) then
                  padd=1.D0-exp(-beta*snsn)                                                 ! check whether parallel
                  if(rkiss05().lt.padd) then
                     sx(neighbor)=sx(neighbor)-nx*scalar2                               ! now flip current spin
                     sy(neighbor)=sy(neighbor)-ny*scalar2
                     stack(sp)=neighbor
                     sp=sp+1
                     isize=isize+1
                  endif
               endif
               endif

               neighbor=m5(current)
               if (occu(neighbor)) then
               scalar2=  2.D0*(nx*sx(neighbor)+ny*sy(neighbor) )   ! scalar product for p_add
               snsn=scalar1*scalar2
               if (snsn>0) then
                  padd=1.D0-exp(-beta*snsn)                                                 ! check whether parallel
                  if(rkiss05().lt.padd) then
                     sx(neighbor)=sx(neighbor)-nx*scalar2                               ! now flip current spin
                     sy(neighbor)=sy(neighbor)-ny*scalar2
                     stack(sp)=neighbor
                     sp=sp+1
                     isize=isize+1
                  endif
               endif
               endif

               neighbor=m4(current)
               if (occu(neighbor)) then
               scalar2=  2.D0*(nx*sx(neighbor)+ny*sy(neighbor) )   ! scalar product for p_add
               snsn=scalar1*scalar2
               if (snsn>0) then
                  padd=1.D0-exp(-beta*snsn)                                                 ! check whether parallel
                  if(rkiss05().lt.padd) then
                     sx(neighbor)=sx(neighbor)-nx*scalar2                               ! now flip current spin
                     sy(neighbor)=sy(neighbor)-ny*scalar2
                     stack(sp)=neighbor
                     sp=sp+1
                     isize=isize+1
                  endif
               endif
               endif

               neighbor=m3(current)
               if (occu(neighbor)) then
               scalar2=  2.D0*(nx*sx(neighbor)+ny*sy(neighbor) )   ! scalar product for p_add
               snsn=scalar1*scalar2
               if (snsn>0) then
                  padd=1.D0-exp(-beta*snsn)                                                 ! check whether parallel
                  if(rkiss05().lt.padd) then
                     sx(neighbor)=sx(neighbor)-nx*scalar2                               ! now flip current spin
                     sy(neighbor)=sy(neighbor)-ny*scalar2
                     stack(sp)=neighbor
                     sp=sp+1
                     isize=isize+1
                  endif
               endif
               endif

               neighbor=m2(current)
               if (occu(neighbor)) then
               scalar2=  2.D0*(nx*sx(neighbor)+ny*sy(neighbor) )   ! scalar product for p_add
               snsn=scalar1*scalar2
               if (snsn>0) then
                  padd=1.D0-exp(-beta*snsn)                                                 ! check whether parallel
                  if(rkiss05().lt.padd) then
                     sx(neighbor)=sx(neighbor)-nx*scalar2                               ! now flip current spin
                     sy(neighbor)=sy(neighbor)-ny*scalar2
                     stack(sp)=neighbor
                     sp=sp+1
                     isize=isize+1
                  endif
               endif
               endif

               neighbor=m1(current)
               if (occu(neighbor)) then
               scalar2=  2.D0*(nx*sx(neighbor)+ny*sy(neighbor) )   ! scalar product for p_add
               snsn=scalar1*scalar2
               if (snsn>0) then
                  padd=1.D0-exp(-beta*snsn)                                                 ! check whether parallel
                  if(rkiss05().lt.padd) then
                     sx(neighbor)=sx(neighbor)-nx*scalar2                               ! now flip current spin
                     sy(neighbor)=sy(neighbor)-ny*scalar2
                     stack(sp)=neighbor
                     sp=sp+1
                     isize=isize+1
                  endif
               endif
               endif


            enddo       ! of cluster building and flipping

         nspflipped   = nspflipped   +isize

         endif          ! of if(occu(is))

         if (ntobeflipped>0) then
            if (nclflipped.ge.ntobeflipped) exit cluster_loop
         else
            if (nspflipped.ge.L3) exit cluster_loop
         endif
      enddo  cluster_loop

      end subroutine wolff_sweep


     subroutine metro_sweep
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!    carries out one Metropolis sweep
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      real(r8b)    :: angle

      do is=0, L3-1
      if (occu(is)) then
          nsx=0.D0
          nsy=0.D0
          if (occu(m1(is))) then
             nsx=nsx+sx(m1(is))
             nsy=nsy+sy(m1(is))
          endif
          if (occu(m2(is))) then
             nsx=nsx+sx(m2(is))
             nsy=nsy+sy(m2(is))
          endif
          if (occu(m3(is))) then
             nsx=nsx+sx(m3(is))
             nsy=nsy+sy(m3(is))
          endif
          if (occu(m4(is))) then
             nsx=nsx+sx(m4(is))
             nsy=nsy+sy(m4(is))
          endif
          if (occu(m5(is))) then
             nsx=nsx+sx(m5(is))
             nsy=nsy+sy(m5(is))
          endif
          if (occu(m6(is))) then
             nsx=nsx+sx(m6(is))
             nsy=nsy+sy(m6(is))
          endif

          slen=1.D0
          do while (slen>0.5D0)
             sxnew= rkiss05()-0.5D0
             synew= rkiss05()-0.5D0
             slen=sqrt(sxnew*sxnew + synew*synew)
          enddo
          slen=1.D0/slen		! changed v8
          sxnew=sxnew*slen
          synew=synew*slen

          dE= nsx*(sx(is)-sxnew)+ nsy*(sy(is)-synew)
          if (dE<0.or.(exp(-dE*beta)>rkiss05())) then
             sx(is)=sxnew
             sy(is)=synew
          endif
      endif     ! of if(occu(is))
      enddo     ! of do is ...


      end subroutine metro_sweep


      subroutine measurement
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! measures energy and magnetization
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      sweepen=0.D0
      mx=0.D0
      my=0.D0
      do is=0,L3-1
         if (occu(is)) then
            mx=mx+sx(is)
            my=my+sy(is)
            nsx=0.D0
            nsy=0.D0
            if (occu(m1(is))) then
               nsx=nsx+sx(m1(is))
               nsy=nsy+sy(m1(is))
            endif
            if (occu(m3(is))) then
               nsx=nsx+sx(m3(is))
               nsy=nsy+sy(m3(is))
            endif
            if (occu(m5(is))) then
               nsx=nsx+sx(m5(is))
               nsy=nsy+sy(m5(is))
            endif
            sweepen=sweepen-sx(is)*nsx
            sweepen=sweepen-sy(is)*nsy
         endif
      enddo
      sweepen=sweepen/L3
      en=en+sweepen
      en2=en2+sweepen**2

      sweepmag=   sqrt(mx*mx+my*my)/L3
      mag= mag +  sweepmag
      mag2=mag2+  sweepmag**2
      mag4=mag4+  sweepmag**4

      enmag=enmag+ sweepen*sweepmag

      end subroutine measurement

      subroutine corr_func
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! measures correlation function
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      real(r8b)   :: Remqtx,Remqty,Immqtx,Immqty
      real(r8b)   :: Remq2x,Remq2y,Immq2x,Immq2y
      real(r8b)   :: Remq3x,Remq3y,Immq3x,Immq3y

      Remqtx=0
      Remqty=0
      Immqtx=0
      Immqty=0

      Remq2x=0
      Remq2y=0
      Immq2x=0
      Immq2y=0

      Remq3x=0
      Remq3y=0
      Immq3x=0
      Immq3y=0
	  is=-1
      do i1=0, Lt-1
      do i2=0, L-1
      do i3=0, L-1
		is = L*(L*i1 + i2) + i3 				! changed v8
         if (occu(is)) then
            Remqtx=Remqtx+sx(is)*costime(i1)	! changed v8
            Remqty=Remqty+sy(is)*costime(i1)
            Immqtx=Immqtx+sx(is)*sintime(i1)
            Immqty=Immqty+sy(is)*sintime(i1)

            Remq2x=Remq2x+sx(is)*cosspace(i2)
            Remq2y=Remq2y+sy(is)*cosspace(i2)
            Immq2x=Immq2x+sx(is)*sinspace(i2)
            Immq2y=Immq2y+sy(is)*sinspace(i2)

            Remq3x=Remq3x+sx(is)*cosspace(i3)
            Remq3y=Remq3y+sy(is)*cosspace(i3)
            Immq3x=Immq3x+sx(is)*sinspace(i3)
            Immq3y=Immq3y+sy(is)*sinspace(i3)
         endif
      enddo
      enddo
      enddo
      sweepmagqt=Remqtx**2+Remqty**2 + Immqtx**2+Immqty**2
      sweepmagqt=sqrt(sweepmagqt)/L3
      magqt=magqt+sweepmagqt
      Gt=Gt+sweepmagqt**2

      sweepmagq2=(Remq2x**2+Remq2y**2 + Immq2x**2+Immq2y**2)
      sweepmagq3=(Remq3x**2+Remq3y**2 + Immq3x**2+Immq3y**2)
      sweepmagq2=sqrt(sweepmagq2)/L3
      sweepmagq3=sqrt(sweepmagq3)/L3
      magqs=magqs+0.5*(sweepmagq2+sweepmagq3)
      Gs=Gs+0.5*(sweepmagq2**2+sweepmagq3**2)

      end subroutine corr_func

      end




!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Random number generator KISS05 after a suggestion by George Marsaglia
! in "Random numbers for C: The END?" posted on sci.crypt.random-numbers
! in 1999
!
! version as in "double precision RNGs" in  sci.math.num-analysis
! http://sci.tech-archive.net/Archive/sci.math.num-analysis/2005-11/msg00352.html
!
! The  KISS (Keep It Simple Stupid) random number generator. Combines:
! (1) The congruential generator x(n)=69069*x(n-1)+1327217885, period 2^32.
! (2) A 3-shift shift-register generator, period 2^32-1,
! (3) Two 16-bit multiply-with-carry generators, period 597273182964842497>2^59
! Overall period > 2^123
!
!
! A   call to rkiss05() gives one random real in the interval [0,1),
! i.e., 0 <= rkiss05 < 1
!
! Before using rkiss05   call kissinit(seed) to initialize
! the generator by random integers produced by Park/Millers
! minimal standard LCG.
! Seed should be any positive integer.
!
! FORTRAN implementation by Thomas Vojta, vojta@mst.edu
! built on a module found at www.fortran.com
!
! History:
!        v0.9     Dec 11, 2010    first implementation
!        V0.91    Dec 11, 2010    inlined internal function for the SR component
!        v0.92    Dec 13, 2010    extra shuffle of seed in kissinit
!        v0.93    Aug 13, 2012    changed integer representation test to avoid data statements
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


      FUNCTION rkiss05()
      implicit none
      integer,parameter      :: r8b= SELECTED_REAL_KIND(P=14,R=99)   ! 8-byte reals
      integer,parameter      :: i4b= SELECTED_INT_KIND(8)            ! 4-byte integers

     real(r8b),parameter    :: am=4.656612873077392578d-10       ! multiplier 1/2^31

      real(r8b)             :: rkiss05
      integer(i4b)          :: kiss
      integer(i4b)          :: x,y,z,w              ! working variables for the four generators
      common /kisscom/x,y,z,w

      x = 69069 * x + 1327217885
      y= ieor (y, ishft (y, 13)); y= ieor (y, ishft (y, -17)); y= ieor (y, ishft (y, 5))
      z = 18000 * iand (z, 65535) + ishft (z, - 16)
      w = 30903 * iand (w, 65535) + ishft (w, - 16)
      kiss = ishft(x + y + ishft (z, 16) + w , -1)
      rkiss05=kiss*am
      END FUNCTION rkiss05


      SUBROUTINE kissinit(iinit)
      implicit none
      integer,parameter      :: r8b= SELECTED_REAL_KIND(P=14,R=99)   ! 8-byte reals
      integer,parameter      :: i4b= SELECTED_INT_KIND(8)            ! 4-byte integers

      integer(i4b) idum,ia,im,iq,ir,iinit
      integer(i4b) k,x,y,z,w,c1
      real(r8b)    rkiss05,rdum
      parameter (ia=16807,im=2147483647,iq=127773,ir=2836)
      common /kisscom/x,y,z,w

      !!! Test integer representation !!!
      c1=-8
      c1=ishftc(c1,-3)
!     print *,c1
      if (c1.ne.536870911) then
         print *,'Nonstandard integer representation. Stoped.'
         stop
      endif

      idum=iinit
      idum= abs(1099087573 * idum)               ! 32-bit LCG to shuffle seeds
      if (idum.eq.0) idum=1
      if (idum.ge.IM) idum=IM-1

      k=(idum)/IQ
      idum=IA*(idum-k*IQ)-IR*k
      if (idum.lt.0) idum = idum + IM
      if (idum.lt.1) then
         x=idum+1
      else
         x=idum
      endif
      k=(idum)/IQ
      idum=IA*(idum-k*IQ)-IR*k
      if (idum.lt.0) idum = idum + IM
      if (idum.lt.1) then
         y=idum+1
      else
         y=idum
      endif
      k=(idum)/IQ
      idum=IA*(idum-k*IQ)-IR*k
      if (idum.lt.0) idum = idum + IM
      if (idum.lt.1) then
         z=idum+1
      else
         z=idum
      endif
      k=(idum)/IQ
      idum=IA*(idum-k*IQ)-IR*k
      if (idum.lt.0) idum = idum + IM
      if (idum.lt.1) then
         w=idum+1
      else
         w=idum
      endif

      rdum=rkiss05()

      return
      end subroutine kissinit
