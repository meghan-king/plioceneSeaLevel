! Pseudo-spectral ice age sea level model, using algorithm of Kendall et al, 2005: On post-glacial sea level - II...
! Written by Sam Goldberg, Harvard University EPS '16 for senior thesis
! Advised by Jerry Mitrovica
! 12/2015-1/2016
! CHANGES (prior to tracking with Git on 2018-01-15):
!  2017-11-22:
!  > Replaced "form = 'binary'" with "form = 'unformatted', access = 'stream', status = '<old>/<replace>', where the
!     last option depends on specific read/write needs." This is to bring the code into compliance with GNU compilers.
!  > Added "status = 'old'" to READ statement for planetary model Love number files (line 174 at the time of this
!     writing).
!  > Ensure everything fits within 120 columns.
!  > Cleaned up the indentations for do...enddo, if...endif, etc.
!  by Erik N.-H. Chan
!  2018-01-12:
!  > Added if..else block around lines 519-525 (at the time of writing) to avoid NaN (= 0 / 0 or dividing by 0).
!  by Erik N.-H. Chan

! For all variables except SL, the prefix d- is equivalent to δ- in Kendall, i.e. incremental change over one timestep.
! The prefix delta- is equivalent to Δ, i.e. total change since first timestep.
! SL is different - dSL refers to spatially heterogeneous change only (the script SL in Kendall), and deltaSL is the
!  total change including eustatic - both since first timestep
! The suffix -xy denotes a quantity in the spatial domain; -lm denotes a quantity in the spectral domain.
! The prefix old- generally refers to the previous timestep, used to calculate increments, although sometimes
!  to previous iteration to check convergence.
! The suffix -0 generally refers to the value at first timestep, used to calculate total changes.
! nouter refers to the iteration of the outer loop (full time series).
! ninner refers to the inner loop within each timestep to converge the ocean load.
! n refers to the timestep, between 1 and nsteps.
! Other than the above, I have endeavored to be faithful to the notation of Kendall et al.

! The INPUT directory should contain the following files:
!  - Ice: named 'ice1','ice2','ice3',...,'ice(nsteps)', containing ice thickness in metres on a Gauss-Legendre grid.
!  - Times: named 'times', containing the times, in years, of each ice files. Present is 0, past is negative, etc.
!  - Present-day bedrock topography: named 'topobr', with values metres on a Gauss-Legendre grid,
!  - Love numbers in JXM's maxwell.f output format: name is specified in the "planetmodel" variable in user_specs_mod.
! The OUTPUT directory will contain the following files:
!  - 'SL1','SL2',...: Relative sea level, in metres, on a Gauss-Legendre grid,
!  - 'times':  Times in years of each output file,
!  - 'R1','G1','R2','G2',...: Radial and Geoid components of RSL change (only if parameter calcRG is enabled).
! NOTE: All of the grid files (i.e., ice, topography, sea-level, radial, geoid) could be in ASCII text or flat binary
!  formats. For these files, the order of values in read and write operations are as follows: Along lines of equal
!  latitude, from 0 degrees towards 360 degrees of longitude; each of these lines are then ordered from 0 colatitude
!  towards 180 (i.e., from North to South).

include 'spharmt.f90' ! Spherical harmonic transform module

!================================================================================================USER SPECIFICATIONS===!
module user_specs_mod
!______________________________________________________________________________________________________________________!

   ! Input directory
   character(*), parameter :: inputfolder  = 'Pliocene Ice/' !MEK Added new input folder

   ! Planetary model directory
   ! The filename of the desired model (i.e., the Love numbers) given by the planetmodel variable below. This is now
   !  incorporated into the planets_mod module below, since automation limits the freedom in naming, which could be
   !  complex and highly customized.
   character(*), parameter :: planetfolder = 'Love numbers/' !MEK Added new folder

   ! Output directory
   character(*), parameter :: outputfolder = 'Outputs/' !MEK Edited Output folder

   ! Common file extensions
   ! Since the code does not explicitly specify the precision (single or double) of the input/output files, an
   !  extension could be used to designate the type. Alternatively, the files could also be delimited text files (to
   !  be implemented later). Specify the common extension here (applicable to the names of ALL input/output files):
   character(4), parameter :: ext = ''    ! '.sgl' | '.dbl' | '.txt' | ''
   ! ... and their file type:
   character(*), parameter :: fType = 'binary'  ! 'binary' | 'text' !MEK changed to binary

   ! Planet selection
   character(*), parameter :: whichplanet = 'earth'                  ! 'earth', 'Mars', etc.
   character(*), parameter :: planetmodel = 'prem.l96C.ump5.lm5' ! MEK Added new folder

   ! Model parameters
   integer, parameter :: norder = 256           ! Max spherical harmonic degree/order
   integer, parameter :: npam = 500             ! Max relaxation modes
   integer, parameter :: nglv = 256             ! Number of GL points in latitude
   integer, parameter :: nsteps = 661           ! Number of timesteps/ice files
   logical, parameter :: checkmarine = .true.   ! T to check for floating marine-based ice
                                                ! F to assume all ice is grounded
   logical, parameter :: tpw = .true.           ! T to incorporate rotational feedback, F for non-_rotating planet
   real, parameter :: epsilon1 = 1.0E-5         ! Inner loop convergence criterion
   real, parameter :: epsilon2 = 1.0E-5         ! Outer loop convergence criterion
                                                !  (if doing a convergence check for outer loop, see below)
   integer, parameter :: numouter = 3           ! Number of outer loops (enter 99 for convergence check) !MEK changed from '1' to '3'
   logical, parameter :: calcRG = .false.       ! T to calculate the radial and geoid displacements, F only to do RSL -
                                                !  only works for fixed number of outer loops (no convergence check)!
end module user_specs_mod

!=================================================================================PHYSICSAL & MATHEMATICAL CONSTANTS===!
module constants_mod
!______________________________________________________________________________________________________________________!
   real, parameter :: pi = 3.1415926535898      ! Pi
   complex, parameter :: ii=(0.0,1.0)           ! Square root of -1
   real, parameter :: gravConst = 6.67408E-11   ! Gravitational constant (m^3/kg/s^2)
end module constants_mod

!========================================================================================================PLANETS_MOD===!
module planets_mod
!______________________________________________________________________________________________________________________!
   use constants_mod

   real :: radius
   real :: mass
   real :: rhoi
   real :: rhow
   real :: gacc
   real :: omega
   real :: acoef, ccoef
   real :: moiA, moiC
   real :: kf

   contains

   !=========================================================================================PLANETS_MOD: EARTH_INIT===!
   subroutine earth_init
   !___________________________________________________________________________________________________________________!
      radius = 6.371E6              ! Radius of the Earth (m)
      mass = 5.976E24               ! Mass of the Earth (kg)
      rhoi = 920.0                  ! Density of ice (kg/m^3)
      rhow = 1000.0                 ! Density of fresh water (kg/m^3)
      gacc = 9.80665                ! Acceleration due to gravity at the Earth's surface (m/s^2)
      omega = 7.292e-5              ! Rotation rate of the Earth (rad/s)
      moiA=0.3296145*mass*radius**2 ! Principal moment of inertia of the Earth
      moiC=0.3307007*mass*radius**2 ! Principal moment of inertia of the Earth
      kf = 0.9342+0.008             ! Fluid (Tidal) Love number

   end subroutine earth_init

   !==========================================================================================PLANETS_MOD: MARS_INIT===!
   subroutine mars_init
   !___________________________________________________________________________________________________________________!
      radius = 3.3899E6             ! Radius of Mars (m)
      mass = 6.4185E23              ! Mass of Mars (kg)
      rhoi = 1220.0                 ! Density of ice mix on Mars (kg/m^3)
      rhow = 1000.0                 ! Density of fresh water (kg/m^3)
      gacc = 3.713                  ! Acceleration due to gravity at Mars's surface (m/s^2)
      omega = 7.08819118E-5         ! Rotation rate of Mars (rad/s)
      moiA=0.363914*mass*radius**2  ! Principal moments of inertia of Mars
      moiC=0.365905*mass*radius**2  ! Principal moments of inertia of Mars
      !----------------------------------------------------------------------------------!
      ! Lithospheric thickness and corresponding kf (based on Zharkov and Gudkova, 2005)
      !    15    |    48    |    58    |    84    |    110   |    164    |    200    [km]
      ! 1.203959 | 1.127384 | 1.110566 | 1.065899 | 1.023186 | 0.9458358 | 0.8986673
      !----------------------------------------------------------------------------------!
      kf = 0.899                    ! Fluid (Tidal) Love number

   end subroutine mars_init

end module planets_mod

!======================================================================================================================!
!                                                      MAIN BLOCK                                                      !
!______________________________________________________________________________________________________________________!

!=======================================================================================================MAIN PROGRAM===!
program sl_model
!______________________________________________________________________________________________________________________!
use spharmt
use user_specs_mod
use planets_mod
implicit none

!==========================================================
!                        VARIABLES
!___________________(Edit with caution)____________________

! Inputs
real, dimension(nglv,2*nglv,nsteps) :: icexy    ! Spatial inputs of ice model
real, dimension(nglv,2*nglv) :: truetopo        ! Present-day topography
real, dimension(nsteps) :: times                ! Timesteps of ice model (years)
real, dimension(npam,norder) :: rprime,r,s      ! Love numbers
real, dimension(norder) :: ke,he                ! Love numbers
real, dimension(npam,norder) :: rprimeT,rT      ! Love numbers (tidal)
real, dimension(norder) :: kTE, hTE             ! Love numbers (tidal)

! Model calculations
real, dimension(nglv,2*nglv,nsteps) :: sl, topoxy, cxy   ! Big arrays of sea level change, topography, ocean function
real, dimension(nglv,2*nglv) :: deltaslxy, dslxy         ! Total sea level change,
                                                         !  total spatially heterogeneous sea level change
real, dimension(nglv,2*nglv) :: icestarxy                ! Grounded ice thickness
real, dimension(nglv,2*nglv) :: beta, cstarxy, cstar0    ! Grounded ice mask, ice-free ocean function
real, dimension(nglv,2*nglv) :: tOxy, rOxy, tTxy         ! Projections used to calculate loads and shoreline migration
complex, dimension(0:norder,0:norder) :: cstarlm,oldcstarlm,tOlm,rOlm,dSlm,olddSlm,&
                                         icestarlm,dicestarlm,deltaicestarlm,oldicestarlm,icestar0, &
                                         t0lm,oldt0lm,tTlm,oldtTlm,dsllm,deltasllm  ! Above, in spectral domain
complex, dimension(0:norder,0:norder,nsteps) :: dicestar,dS,deltaicestar,deltaS     ! Big arrays of changes in loads,
                                                                                    !  used in Love number viscous
                                                                                    !  response
real :: conserv                                          ! Uniform geoid shift (ΔΦ/g)
real :: ttl, ekhl                                        ! Used in Love number calculations
real, dimension(nsteps,norder) :: lovebeta               ! Used in Love number calculations
complex, dimension(0:norder,0:norder) :: viscous         ! Used in Love number calculations
real :: xi, zeta                                         ! Convergence checks

! For calculating R and G separately
real, dimension(nglv,2*nglv) :: rrxy, drrxy_computed
real, dimension(nglv,2*nglv,nsteps) :: rr,gg             ! R and G (radial displacement and geoid change)
complex, dimension(0:norder,0:norder) :: rrlm, dgglm, drrlm_computed
real, dimension(nsteps,norder) :: lovebetarr
complex :: viscousrr

! Rotation calculations
! real, parameter :: --------------------------------->     ! Fluid Love number is defined in the planetary modules
complex, dimension(0:2,0:2) :: lambda, oldlambda, lambda0   ! Rotational driving
complex, dimension(0:2,0:2,nsteps) :: dlambda, deltalambda  ! Big arrays of changes in rotational driving
real, dimension(3) :: mm, sum_m, oldm                       ! Rotational perturbation vector
real, dimension(3,nsteps) :: dm                             ! Big array of changes in m
real, dimension(3,3) :: il, oldil, sum_il                   ! Load perturbations to MOI
real, dimension(3,3,nsteps) :: dil                          ! Big array of changes in IL
complex, dimension(0:2,0:2) :: dsl_rot, rr_rot              ! Sea level change from rotational feedbacks
real :: ekhTE                                               ! Used in Love number calculations
real, dimension(nsteps) :: lovebetatt, lovebetattrr         ! Used in Love number calculations
real :: betatt, betattprime                                 ! Used in Love number calculations
complex :: viscoustt,viscousttrr                            ! Used in Love number calculations

! Miscellaneous variables
integer :: ninner,nouter                                    ! Iteration of inner and outer loops
integer :: i,j,k,l,m,n,nn                                   ! Do-loop indices
character(3) :: numstr, iterstr                             ! String for timestep number for reading/writing files
integer :: counti, countf,counti1,countf1,countrate         ! Computation timing
type(sphere) :: spheredat                                   ! SH transform data to be passed to subroutines

! For Jerry's code to read in Love numbers
integer :: legord(norder),nmod(norder),nmodes(norder),ll,nm,np
real :: xn
real, dimension(3,norder) :: elast,asymv,telast,tasymv
real, dimension(npam,norder) :: resh,resl,resk,tresh,tresl,tresk
real :: taurr,taurt,dmx,dmy

! Planetary values
if (whichplanet == 'earth' .or. whichplanet == 'Earth' .or. whichplanet == 'EARTH') then
   call earth_init
elseif (whichplanet == 'mars' .or. whichplanet == 'Mars' .or. whichplanet == 'MARS') then
   call mars_init
else
   write(*,*) 'The parameters for the planet you entered are not built in.'
   write(*,*) 'Please check your spelling for the variable whichplanet in the user_specs_mod module.'
   write(*,*) 'If you preferred, you could create a new subroutine for the new planet in the planets_mod module.'
   write(*,*) 'Terminating: program sl_model'
   stop
endif

!===========================================================
!                      BEGIN EXECUTION
!___________________________________________________________

call system_clock(count=counti, count_rate = countrate)  ! Total computation time
call spharmt_init(spheredat,2*nglv,nglv,norder,radius)   ! Initialize spheredat (for SH transform subroutines)

! Read in input files
write(*,'(A)') 'Reading input files'


! Ice files
do n = 1,nsteps
   write(numstr,'(I3)') n !MEK changed from -1 to nothing
   numstr = trim(adjustl(numstr))
   if (fType == 'text') then        ! If the files are space-delimited text files
      open(unit = 1, file = inputfolder//'ice'//trim(numstr)//ext, form = 'formatted', access = 'sequential', &
      & status = 'old')
      read(1,*) icexy(:,:,n)
   elseif (fType == 'binary') then  ! If the files are "flat" binary
      open(unit = 1, file = inputfolder//'ice'//trim(numstr)//ext, form = 'unformatted', access = 'stream', &
      & status = 'old')
      read(1) icexy(:,:,n)
   endif
   close(1)
enddo


! Timesteps
if (fType == 'text') then        ! If the files are space-delimited text files
   open(unit = 1, file = inputfolder//'times'//ext, form = 'formatted', access = 'sequential', status = 'old')
   open(unit = 2, file = outputfolder//'times'//ext, form = 'formatted', access = 'sequential', status = 'replace')
   read(1,*) times
   ! Write to output folder for plotting results, etc.
   write(2,'(ES14.4E2)') times
elseif (fType == 'binary') then  ! If the files are "flat" binary
   open(unit = 1, file = inputfolder//'times'//ext,  form = 'unformatted', access = 'stream', status = 'old')
   open(unit = 2, file = outputfolder//'times'//ext, form = 'unformatted', access = 'stream', status = 'replace')
   read(1) times
   ! Write to output folder for plotting results, etc.
   write(2) times
endif
close(1)
close(2)

! Present topography
open(unit = 1, file = planetfolder//'topobr'//ext, form = 'formatted', access = 'sequential', status = 'old')
read(1,*) truetopo
!if (fType == 'text') then        ! If the files are space-delimited text files
   !open(unit = 1, file = planetfolder//'topobr'//ext, form = 'formatted', access = 'sequential', status = 'old')
   !read(1,*) truetopo
!elseif (fType == 'binary') then  ! If the files are "flat" binary
   !open(unit = 1, file = planetfolder//'topobr'//ext, form = 'unformatted', access = 'stream', status = 'old')
   !read(1) truetopo
!endif
close(1)
! truetopo(:,:)=truetopo(:,:)-icexy(:,:,nsteps) ! If using ice topography

! Read in Love numbers (Jerry's output from 'maxwell.f')
open(unit = 2, file = planetfolder//planetmodel, status = 'old')
! Following code borrowed from Jerry
read(2,*)
do j = 1,norder
   read(2,*) legord(j),nmodes(j)
   nm=nmodes(j)
   xn=real(legord(j))
   ll=legord(j)
   nmod(ll)=nm
   read(2,*) (s(i,ll),i=1,nm)
   read(2,*) elast(1,ll),elast(2,ll),taurr,taurt,elast(3,ll)
   read(2,*) asymv(1,ll),asymv(2,ll),taurr,taurt,asymv(3,ll)
   read(2,*) (resh(i,ll),i=1,nm)
   read(2,*) (resl(i,ll),i=1,nm)
   read(2,*) (resk(i,ll),i=1,nm)
   if(xn .lt. 1.5) cycle
   read(2,*) telast(1,ll),telast(2,ll),dmx,dmy,telast(3,ll)
   read(2,*) tasymv(1,ll),tasymv(2,ll),dmx,dmy,tasymv(3,ll)
   read(2,*) (tresh(i,ll),i=1,nm)
   read(2,*) (tresl(i,ll),i=1,nm)
   read(2,*) (tresk(i,ll),i=1,nm)
enddo
! 1001  format(2i5)
! 1002  format(20a4)
! 1003  format(i10,1x,i5)
! 1020  format(5e16.8)
close(2)

! Divide by l (numbers are multiplied by l in Jerry's output)
do l = 1,norder
   resl(:,l)=resl(:,l)/real(l)
   resk(:,l)=resk(:,l)/real(l)
   tresl(:,l)=tresl(:,l)/real(l)
   tresk(:,l)=tresk(:,l)/real(l)
enddo

! Assign Love number inputs to the letters used in Kendall et al
do l = 1,norder
   he(l) = elast(1,l)
   ke(l) = elast(3,l)/real(l)
   hTE(l) = telast(1,l)
   kTE(l) = telast(3,l)/real(l)
enddo
rprime(:,:)=resk(:,:)
r(:,:)=resh(:,:)
rprimeT(:,:)=tresk(:,:)
rT(:,:)=tresh(:,:)

!===========================================================
!                       CALCULATIONS
!___________________________________________________________

! Initialize topography to be present topography (STEP 1)
do n=1,nsteps
   topoxy(:,:,n)=truetopo(:,:) ! (eq. 48)
enddo
! Calculate ocean function from present topography as first guess
do j = 1,2*nglv
   do i = 1,nglv
      if (truetopo(i,j)<0) then
         cxy(i,j,:) = 1
      else
         cxy(i,j,:) = 0
      endif
   enddo
enddo

! Decompose initial topography (STEP 2) (used to check convergence of outer loop)
call spat2spec(topoxy(:,:,1),t0lm,spheredat)

! BEGIN OUTER LOOP
do nouter = 1,numouter

   write(*,'(A,I2)') 'Outer loop iteration ', nouter
   call system_clock(counti1) ! Time for each outer loop

   do n = 1,nsteps ! Time iteration
      write(*,'(A,I4,A,EN15.4E2,A)') '  Timestep',n,', ',times(n),' years'

      ! Calculate icestar (STEP 3) (eq.43)
      if (checkmarine) then
         do j = 1,2*nglv
            do i = 1,nglv
                if (topoxy(i,j,n)>0) then ! If not marine
                   icestarxy(i,j) = icexy(i,j,n)
                elseif (icexy(i,j,n)>(abs(topoxy(i,j,n))*rhow/rhoi)) then ! If marine, but thick enough to be grounded
                   icestarxy(i,j) = icexy(i,j,n)
                else ! If floating ice
                   icestarxy(i,j) = 0
                endif
            enddo
         enddo
      else ! If not checking for floating ice
         icestarxy(:,:)=icexy(:,:,n)
      endif

      write(numstr,'(I3)') n
      numstr = trim(adjustl(numstr))
      if (fType == 'text') then        ! If the files are space-delimited text files
         open(unit = 1, file = outputfolder//'ICE'//numstr, form = 'formatted', access = 'sequential', &
         & status = 'replace')
         write(1,'(ES14.4E2)') icestarxy(:,:) ! MEK Removed topoxy(:,:,n) + icexy(:,:,n)
      elseif (fType == 'binary') then  ! If the files are "flat" binary
         open(unit = 1, file = outputfolder//'ICE'//numstr, form = 'unformatted', access = 'stream', &
         & status = 'replace')
         write(1) icestarxy(:,:) ! MEK Removed topoxy(:,:,n) + icexy(:,:,n)
      endif
      close(1)

      ! Calculate beta (STEP 3) (eq. 44)
      do j = 1,2*nglv
         do i = 1,nglv
            if (icestarxy(i,j)==0) then
               beta(i,j)=1
            else
               beta(i,j)=0
            endif
         enddo
      enddo

      call spat2spec(icestarxy,icestarlm,spheredat) ! Decompose ice field

      if (n==1) then
         dicestarlm(:,:)=0.0 ! No change at first timestep
         icestar0(:,:)=icestarlm(:,:) ! Save initial to calculate total change
      else
         dicestarlm(:,:)=icestarlm(:,:)-oldicestarlm(:,:) ! Incremental change
      endif
      oldicestarlm(:,:)=icestarlm(:,:) ! Save to calculate increment on next timestep
      deltaicestarlm(:,:)=icestarlm(:,:)-icestar0(:,:) ! Total change since time0
      dicestar(:,:,n)=dicestarlm(:,:) ! Save into big matrix
      deltaicestar(:,:,n)=deltaicestarlm(:,:) ! Save into big matrix

      ! Calculate cstar (STEP 4)
      cstarxy(:,:)=cxy(:,:,n)*beta(:,:) ! (eq. 65)
      if (n==1) then
         cstar0 = cstarxy ! Value at time0, used in tO calculation below for shoreline migration
      endif
      oldcstarlm(:,:) = cstarlm(:,:) ! Save previous to use in eq. 79
      call spat2spec(cstarxy,cstarlm,spheredat) ! Decompose cstar
      if (n==1) then
         oldcstarlm(:,:) = cstarlm(:,:) ! So that there will be no change at time0
      endif
      ! Calculate tO (STEP 4)
      tOxy(:,:)=topoxy(:,:,1)*(cstarxy(:,:)-cstar0(:,:)) ! (eq. 70)
      call spat2spec(tOxy,tOlm,spheredat) ! Decompose tO

      ! Calculate initial dS guess for first iteration (STEP 5)
      if (nouter==1) then
         if (n==1) then
            dSlm(:,:) = (0.0,0.0) ! No difference on first timestep
            tTxy(:,:)=truetopo(:,:)*cstarxy(:,:) ! (eq. 80) - tT is the script T in Kendall
            call spat2spec(tTxy,tTlm,spheredat) ! Decompose tT
            oldtTlm(:,:) = tTlm(:,:) ! Save to calculate difference on next timestep
         else
            tTxy(:,:)=truetopo(:,:)*cstarxy(:,:) ! (eq. 80) - tT is the script T in Kendall
            call spat2spec(tTxy,tTlm,spheredat) ! Decompose tT
            dSlm(:,:) = (oldcstarlm(:,:)/oldcstarlm(0,0)) &
                         * ((-rhoi/rhow)*dicestarlm(0,0) + tTlm(0,0) - oldtTlm(0,0)) &
                         -(tTlm(:,:)-oldtTlm(:,:))  ! (eq. 79)
            oldtTlm(:,:) = tTlm(:,:) ! Save to calculate difference on next timestep
         endif
      else
         dSlm(:,:) = dS(:,:,n) ! Use converged load from last outer loop (eq. 82)
      endif

    ! BEGIN INNER LOOP
    ninner = 1
    200 CONTINUE

      !------\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/    Rotation    \/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/------!

      if (tpw) then ! If incorporating rotation

         ! MOI perturbations (Matsuyama et al, 2006) (only need (1,3),(2,3), and (3,3))
         !il(1,1) = 4*pi*radius**4*( 1/(3*sqrt(5.0))*(rhow*deltaS(2,0,n)+rhoi*deltaicestar(2,0,n)) &
         !           - sqrt(2.0/15.0) * real(rhow*deltaS(2,2,n) + rhoi*deltaicestar(2,2,n)) )
         !il(2,2) = 4*pi*radius**4*( 1/(3*sqrt(5.0)) * (rhow*deltaS(2,0,n)+rhoi*deltaicestar(2,0,n)) &
         !           + sqrt(2.0/15.0) * real(rhow*deltaS(2,2,n) + rhoi*deltaicestar(2,2,n)) )
         il(3,3) = (-8*pi*radius**4/(3*sqrt(5.0))) * (rhow*deltaS(2,0,n)+rhoi*deltaicestar(2,0,n))
         write(*,*) 'deltaS(2,0,n)', deltaS(2,0,n)
         !il(1,2) = (8*pi*radius**4/(sqrt(30.0))) * real(-1.0*ii*(rhow*deltaS(2,2,n)+rhoi*deltaicestar(2,2,n)))
         il(1,3) = (8*pi*radius**4/(sqrt(30.0))) * real(rhow*deltaS(2,1,n)+rhoi*deltaicestar(2,1,n))
         il(2,3) = (-8*pi*radius**4/(sqrt(30.0))) * real(-1.0*ii*(rhow*deltaS(2,1,n)+rhoi*deltaicestar(2,1,n)))
         !il(2,1) = il(1,2)
         !il(3,1) = il(1,3)

         if (n==1) then
            dil(:,:,n)=il(:,:)
         else
            dil(:,:,n)=il(:,:)-oldil(:,:)
         endif
         oldil(:,:)=il(:,:)

         ! Calculate m (rotation vector perturbation)
         !  - From Jerry's written notes, also Mitrovica, Wahr, Matsuyama, and Paulson (2005) eq. 2
         sum_il(:,:) = 0.0
         sum_m(:) = 0.0
         do nn = 1,n-1 ! Sum over all previous timesteps
            betatt = 0.0
            betattprime = 0.0
            do k = 1,nmod(2) ! Sum over k=1,K
               betatt = betatt+(rprime(k,2)/s(k,2))*(1.0-exp(-1.0*s(k,2)*(times(n)-times(nn))/1000.0))
               betattprime = betattprime+(rprimeT(k,2)/s(k,2))*(1.0-exp(-1.0*s(k,2)*(times(n)-times(nn))/1000.0))
            enddo
            sum_il(:,:) = sum_il(:,:)+dil(:,:,nn)*betatt
            sum_m(:) = sum_m(:)+dm(:,nn)*betattprime
         enddo
         do i = 1,2
            mm(i)=(1/(1-kTE(2)/kf))*((1/(moiC-moiA))*((1+kE(2))*il(i,3)+sum_il(i,3))+(1/kf)*sum_m(i))
         enddo
         mm(3) = (1.0/moiC)*((1+kE(2))*il(3,3)+sum_il(3,3))

         if (n==1) then
            dm(:,n)=mm(:)
         else
            dm(:,n)=mm(:)-oldm(:)
         endif
         oldm(:)=mm(:)

         ! Calculate lambda (rotational driving) from m (Milne and Mitrovica 1998)
         lambda(0,0) = (radius**2*omega**2/3.0)*((mm(1)**2+mm(2)**2+mm(3)**2)+2.0*mm(3))
         lambda(2,0) = (radius**2*omega**2/(6.0*sqrt(5.0)))*(mm(1)**2+mm(2)**2-2.0*mm(3)**2-4.0*mm(3))
         lambda(2,1) = (radius**2*omega**2/sqrt(30.0))*((mm(1)-ii*mm(2))*(1.0+mm(3)))
         lambda(2,2) = (radius**2*omega**2/sqrt(5.0*24.0))*(mm(2)**2-mm(1)**2+2.0*ii*mm(1)*mm(2))

         ! Simplified linear form from Mitrovica, Milne and Davis (2001) - makes negligible difference
         !lambda(0,0) = (2.0*radius**2*omega**2/3.0)*mm(3)
         !lambda(2,0) = (-2.0*radius**2*omega**2/(3*sqrt(5.0)))*mm(3)
         !lambda(2,1) = (radius**2*omega**2/sqrt(30.0))*(mm(1)-ii*mm(2))
         !lambda(2,2) = 0.0

         if (n==1) then
           lambda0(:,:)=lambda(:,:)
           dlambda(:,:,n)=(0.0,0.0)
           deltalambda(:,:,n)=(0.0,0.0)
         else
           dlambda(:,:,n)=lambda(:,:)-oldlambda(:,:)
           deltalambda(:,:,n)=lambda(:,:)-lambda0(:,:)
         endif
         oldlambda(:,:) = lambda(:,:)

         ! Calculate effect on sea level (Love numbers) (Kendall)
         dsl_rot(:,:)=(0.0,0.0)
         if (n/=1) then
            ekhTE = 1+kTE(2)-hTE(2)
            do nn = 1,n-1 ! Sum over all previous timesteps
               lovebetatt(nn) = 0.0
               do k = 1,nmod(2) ! Sum over k=1,K
                  lovebetatt(nn) = lovebetatt(nn) + ((rprimeT(k,2)-rT(k,2))/s(k,2)) &
                                   * (1-exp(-1.0*s(k,2)*(times(n)-times(nn))/1000.0)) ! (eq. B27)
               enddo
            enddo
            do m = 0,2
               viscoustt = (0.0,0.0)
               do nn = 1,n-1 ! Sum the loads over all previous timesteps to get the sum on the fourth line of eq. B28
                  viscoustt = viscoustt + lovebetatt(nn)*dlambda(2,m,nn)
               enddo
               dsl_rot(2,m) = (1/gacc)*(ekhTE*(deltalambda(2,m,n-1)+dlambda(2,m,n))+viscoustt) ! (eq. B28/B25)
            enddo
         endif

         if (nouter==numouter .and. calcRG) then ! For R calculations
            rr_rot(:,:) = (0.0,0.0)
            if (n/=1) then
               do nn = 1,n-1 ! Sum over all previous timesteps
                  lovebetattrr(nn) = 0.0
                  do k = 1,nmod(2) ! Sum over k=1,K
                     lovebetattrr(nn) = lovebetattrr(nn) &
                                       + ((rT(k,2))/s(k,2)) &
                                       * (1-exp(-1.0*s(k,2)*(times(n)-times(nn))/1000.0)) ! (eq. B27)
                  enddo
               enddo
               do m = 0,2
                  viscousttrr = (0.0,0.0)
                  do nn = 1,n-1 ! Sum the loads over all previous timesteps to get the sum on the fourth line of eq. B28
                     viscousttrr = viscousttrr + lovebetattrr(nn)*dlambda(2,m,nn)
                  enddo
                  rr_rot(2,m) = (1/gacc)*(hTE(2)*(deltalambda(2,m,n-1)+dlambda(2,m,n))+viscousttrr)
               enddo
            endif
         endif

      else ! If not incorporating rotation
         dsl_rot(:,:) = (0.0,0.0)
         rr_rot(:,:) = (0.0,0.0)
      endif ! tpw

      !------/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\    End Rotation    /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\------!

      ! Calculate beta (eq. B14) - viscous response factor
      do l = 1,norder
         do nn = 1,n-1 ! Sum over all previous timesteps
            lovebeta(nn,l) = 0.0
            do k = 1,nmod(l) ! Sum over k=1,K
               lovebeta(nn,l) = lovebeta(nn,l) + ((rprime(k,l)-r(k,l))/s(k,l)) &
                                 * (1-exp(-1.0*s(k,l)*(times(n)-times(nn))/1000.0)) ! (eq. B14)
            enddo
         enddo
      enddo

      if (nouter==numouter .and. calcRG) then ! For R calculations
         do l = 1,norder
            do nn = 1,n-1 ! Sum over all previous timesteps
               lovebetarr(nn,l) = 0.0
               do k = 1,nmod(l) ! Sum over k=1,K (modes)
                  lovebetarr(nn,l) = lovebetarr(nn,l) &
                                     + ((r(k,l))/s(k,l)) * (1-exp(-1.0*s(k,l)*(times(n)-times(nn))/1000.0)) ! (eq. B14)
               enddo
            enddo
         enddo
      endif

      ! Compute viscous response outside inner loop, since it doesn't depend on current timestep
      viscous(:,:) = (0.0,0.0)
      do l = 1,norder
         do m = 0,l
            do nn = 1,n-1 ! Sum the loads over all previous timesteps to get the sum on the second line of eq. B18
               viscous(l,m) = viscous(l,m) + lovebeta(nn,l)*(rhoi*dicestar(l,m,nn)+rhow*dS(l,m,nn))
            enddo
         enddo
      enddo



      ! Calculate change in sea level from change in load (STEP 6)
      if (n==1) then
         dsllm(:,:) = (0.0,0.0) ! No change on first timestep
      else
         do l = 1,norder
            ttl = (4.0*pi*(radius**3))/(mass*(2*real(l)+1.0)) ! (eq. B16)
            ekhl = 1 + ke(l) - he(l) ! (eq. B13)
            do m = 0,l
               ! Total  = _____________________________elastic_____________________________ + _viscous_
               dsllm(l,m) = ttl * ekhl * (rhoi*deltaicestar(l,m,n) + rhow*deltaS(l,m,n-1) + rhow*dSlm(l,m)) &
                             + ttl * viscous(l,m) ! (eq. B18)
            enddo
         enddo
      endif

      ! Add rotational effects (calculated above)
      dsllm(0:2,0:2)=dsllm(0:2,0:2)+dsl_rot(0:2,0:2)

      ! Convert dSL (total spatially heterogeneous change since time0) to spatial domain
      call spec2spat(dslxy,dsllm,spheredat)

      ! Compute r0 (STEP 7)
      rOxy(:,:) = dslxy(:,:)*cstarxy(:,:) ! (eq. 68)
      call spat2spec(rOxy,rOlm,spheredat) ! Decompose r0

      ! Compute conservation term (STEP 8)
      conserv = (1/cstarlm(0,0))*(-1.0*(rhoi/rhow)*deltaicestar(0,0,n)-rOlm(0,0)+tOlm(0,0)) ! (eq. 78)

      ! Compute change in ocean load again (STEP 9)
      olddSlm = dSlm ! Save previous-iterate load to check convergence
      if (n==1) then
         dSlm(:,:)=(0.0,0.0)     ! No change at first timestep
      else
         dSlm(:,:)=-1.0*deltaS(:,:,n-1)+rOlm(:,:)+conserv*cstarlm(:,:)-tOlm(:,:) ! (eq. 73)
      endif

      ! Calculate convergence criterion for inner loop
      if ((sum(abs(dSlm))-sum(abs(olddSlm)))==0.0 .and. sum(abs(olddSlm))==0.0) then
         xi = 0 ! Otherwise xi = 0 / 0 = NaN at the first loop of the first timestep.
      elseif (sum(abs(olddSlm))==0.0) then
         xi = abs((sum(abs(dSlm))-sum(abs(olddSlm)))/1E-8) ! Avoid dividing by 0
      else
         xi = abs((sum(abs(dSlm))-sum(abs(olddSlm)))/sum(abs(olddSlm))) ! (eq. 83)
      endif
      if (xi>epsilon1) then ! If not converged
         ninner = ninner+1

      ! Save total load changes, used in Love number calculation
      if (n==1) then
        deltaS(:,:,n)=dSlm(:,:)
        else
        deltaS(:,:,n)=deltaS(:,:,n-1)+dSlm(:,:)
      endif

         GOTO 200 ! Inner loop
      endif
      write(*,'(A,I3,A)') '  ',ninner,' iterations' ! End inner loop

      if (nouter==numouter .and. calcRG) then ! For R calculations
         if (n==1) then
            rrlm(:,:) = (0.0,0.0)! No change on first timestep
         else
            do l = 1,norder
               ttl = (4.0*pi*(radius**3))/(mass*(2.0*real(l)+1.0)) ! (eq. B16)
               do m = 0,l
                  ! Viscous response
                  viscousrr = (0.0,0.0)
                  do nn = 1,n-1 ! Sum the loads over all previous timesteps to get the sum on the second line of eq. B18
                     viscousrr = viscousrr + lovebetarr(nn,l)*(rhoi*dicestar(l,m,nn) + rhow*dS(l,m,nn))
                  enddo
                  rrlm(l,m) = ttl*he(l)*(rhoi*deltaicestar(l,m,n)+rhow*deltaS(l,m,n-1)+rhow*dSlm(l,m))+ttl*viscousrr
               enddo
            enddo
         endif
         rrlm(0:2,0:2) =rrlm(0:2,0:2)+rr_rot(0:2,0:2)
         call spec2spat(rrxy,rrlm,spheredat)
         rr(:,:,n) = rrxy(:,:)
      endif

      dS(:,:,n)=dSlm(:,:) ! Save converged load into big matrix

      ! Compute total sea level change (STEP 10) (eq. 85)
      deltasllm(:,:) = dsllm(:,:) ! Spatially heterogeneous component
      deltasllm(0,0) = deltasllm(0,0) + conserv ! Add uniform conservation term to (0,0)
      call spec2spat(deltaslxy,deltasllm,spheredat) ! Synthesize deltasl

      ! Save into big array
      sl(:,:,n) = deltaslxy(:,:)


   enddo ! n = 1,nsteps -> Time iteration

   ! Update topography fieldS, ocean functions (STEP 11)
   oldt0lm = t0lm ! Save old initial topography to check convergence
   do n = 1,nsteps
      topoxy(:,:,n) = truetopo(:,:)+sl(:,:,nsteps)-sl(:,:,n) ! (eq. 39)
      ! Update ocean function
      do j = 1,2*nglv
         do i = 1,nglv
            if (topoxy(i,j,n)>=0.0) then
               cxy(i,j,n) = 0.0
            else
               cxy(i,j,n) = 1.0
            endif
         enddo
      enddo
   enddo

   ! Decompose initial topography to check convergence (STEP 12)
   call spat2spec(topoxy(:,:,1),t0lm,spheredat)

   call system_clock(countf1) ! Time for each loop
   write(*,'(A,F6.2,A)') '  Loop time ',real(countf1-counti1)/real(countrate),' seconds'

   ! Calculate convergence criterion (if using convergence check for outer loop)
   if (numouter==99) then
      zeta = abs((SUM(abs(t0lm))-SUM(abs(oldt0lm)))/SUM(abs(oldt0lm))) ! (eq. 86)
      if (zeta<epsilon2) then ! If converged
         EXIT ! Get out of outer loop
      endif
   endif

enddo ! nouter = 1,numouter (Outer loop)

300 continue
if (numouter==99) then
   write(*,*)
   write(*,*) '  **************'
   write(*,*) '  * Converged! *'
   write(*,*) '  **************'
   write(*,*)
   write(*,'(I3,A)') nouter,' topography iterations'
   write(*,*)
endif


!===========================================================
!                          OUTPUT
!___________________________________________________________

write(*,'(A)') 'Writing output files...'

do n = 1,nsteps
   write(numstr,'(I3)') n
   numstr = trim(adjustl(numstr))
   if (fType == 'text') then        ! If the files are space-delimited text files
      open(unit = 1, file = outputfolder//'SL'//trim(numstr)//ext, form = 'formatted', access = 'sequential', &
      & status = 'replace')
      write(1,'(ES14.4E2)') sl(:,:,n) - sl(:,:,1) ! Subtract present to get RSL !MEK changed '1' to 'nsteps'
                                                !NG: for future, present is at the start of the run and the field
!                                               !"sl(:,:,1)" is zero.  for ice age code, subtract sl(:,:,nsteps).
   elseif (fType == 'binary') then  ! If the files are "flat" binary
      open(unit = 1, file = outputfolder//'SL'//trim(numstr)//ext, form = 'unformatted', access = 'stream', &
      & status = 'replace')
      write(1) sl(:,:,n) - sl(:,:,1) ! Subtract present to get RSL !MEK changed '1' to 'nsteps'
   endif
   close(1)

!   if (calcRG) then
!
!      if (fType == 'text') then        ! If the files are space-delimited text files
!         open(unit = 1, file = outputfolder//'R'//trim(numstr)//ext, form = 'formatted', access = 'sequential', &
!         & status = 'replace')
!         write(1,'(ES14.4E2)') rr(:,:,n)
!      elseif (fType == 'binary') then  ! If the files are "flat" binary
!         open(unit = 1, file = outputfolder//'R'//trim(numstr)//ext, form = 'unformatted', access = 'stream', &
!         & status = 'replace')
!         write(1) rr(:,:,n)
!      endif
!      close(1)
!
!      ! Compute geoid displacement
!      gg(:,:,n) = sl(:,:,n)+rr(:,:,n)

!      if (fType == 'text') then        ! If the files are space-delimited text files
!         open(unit = 1, file = outputfolder//'G'//trim(numstr)//ext, form = 'formatted', access = 'sequential', &
!         & status = 'replace')
!         write(1,'(ES14.4E2)') gg(:,:,n)
!      elseif (fType == 'binary') then  ! If the files are "flat" binary
!         open(unit = 1, file = outputfolder//'G'//trim(numstr)//ext, form = 'unformatted', access = 'stream', &
!         & status = 'replace')
!         write(1) gg(:,:,n)
!      endif
!      close(1)
!
!   endif

!   if (fType == 'text') then        ! If the files are space-delimited text files
!      open(unit = 1, file = outputfolder//'ICE'//numstr, form = 'formatted', access = 'sequential', &
!      & status = 'replace')
!      write(1,'(ES14.4E2)') icestarxy(:,:) ! MEK Removed topoxy(:,:,n) + icexy(:,:,n)
!   elseif (fType == 'binary') then  ! If the files are "flat" binary
!      open(unit = 1, file = outputfolder//'ICE'//numstr, form = 'unformatted', access = 'stream', &
!      & status = 'replace')
!      write(1) icestarxy(:,:) ! MEK Removed topoxy(:,:,n) + icexy(:,:,n)
!   endif
!   close(1)

!   if (fType == 'text') then        ! If the files are space-delimited text files
!      open(unit = 1, file = outputfolder//'TOPO'//numstr, form = 'formatted', access = 'sequential', &
!      & status = 'replace')
!      write(1,'(ES14.4E2)') topoxy(:,:,n)
!   elseif (fType == 'binary') then  ! If the files are "flat" binary
!      open(unit = 1, file = outputfolder//'TOPO'//numstr, form = 'unformatted', access = 'stream', &
!      & status = 'replace')
!      write(1) topoxy(:,:,n)
!   endif
!   close(1)
enddo

call system_clock(countf) ! Total time
write(*,'(A,F7.2,A)') 'Done! Total time ',real(countf-counti)/real(countrate),' seconds'

end program sl_model
