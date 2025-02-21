program oxmodel
   implicit none
   
   integer maxzanz, maxtimeanz, maxlatanz   !needed for the dimensions of the arrays below
   parameter(maxlatanz=13)
   parameter(maxzanz=20)         !maximum number of levels in the model
   parameter(maxtimeanz=500000)  !maximum number of integration time steps
   
   !meteorology (these variables contain the meteorological parameters
   !of the standard atmosphere. They are set in the routine readatmosphaere):
   real p(maxzanz), t(maxzanz), air(maxzanz), o2(maxzanz)
         
   !species (These are the values that are calculated by the model):
   real o(maxtimeanz,maxzanz)      !concentration of O in molecules per cm^3
   real o3(maxtimeanz,maxzanz)     !concentration of O3 in molecules per cm^3
   real ox(maxtimeanz,maxzanz)     !concentration of Ox in molecules per cm^3
   real topo3(maxtimeanz,maxzanz)  !overhead ozone in DU

   !additional things that are calculated during the model run and need to be saved:
   integer year(maxtimeanz),month(maxtimeanz),day(maxtimeanz) !current date (greogorian calendar, UT)
   real time(maxtimeanz) !current time (UT)
   real sza(maxtimeanz) !current solar zenith angle (degree)

   !reaction constants and photolysis frequencies (have to be calculated):
   real KO_O2   !reaction constant for O+O2+M->O3+M
   real KO_O3   !reaction constant for O+O3->O2+O2
   real JO2   !reaction constant for O2+hv->O+O
   real JO3   !reaction constant for O3+hv->O2+O

   !variables that are needed for the functions make_JO2 und make_JO3
   !They are initialised in the section "Program initialisation".
   real JO2log_array(24, 30, 40, 5)   !lookup table f�r log(JO2)array
   real JO3log_array(24, 30, 40, 5)   !lookup table f�r log(JO3)
   real jgridz(24),jgridsza(30),jgridtopo3(40),jgridalbedo(5) !Definition of the grid of the lookup tables
   integer jgridzanz,jgridszaanz,jgridtopo3anz,jgridalbedoanz !Definition of the grid of the lookup tables (size of the lookup table in all dimensions)
   
   !model level, model parameters and start date and time:
   real z(maxzanz) !model level (km)
   real lat(maxlatanz)
   real lon !longitude and latitude of the  model run (degree)
   real albedo(maxlatanz) !albedo at the location of the model run
   integer startyear,startmonth,startday !start date of the model run (gregorian datum in UT)
   real starttime !start time of the model run (UT)
   
   !variables for the integration of the model
   real integrationtime !length of the model run in hours
   real timestep !timestep of the integration in hours

   !indizes for the loops in the model
   integer zi !index of the loop over the model levels
   integer zanz !number of model levels (maximal value for zi)
   integer*8 ti !index for the timesteps
   integer tianz !number of timesteps (can be calculated from integration time and timestep)
   integer li !index of the loop over the model latitudes
   integer latanz

   !output
   character*4 filename
   character*8 filename_lat
   character*8 filespec_albedo

   !functions
   real make_JO2, make_JO3
   real make_sza
   real make_KO_O2, make_KO_O3

   !********
   ! Setup:
   !********

   !some basic parameters:
   albedo=0.25           !fill in number (real, albedo at the location of the model run)
   lon=0.              !fill in number (real, longitude of the model run in degrees)
   
   !latitude grid.
   latanz=13
   lat(1)=90.
   lat(2)=75.
   lat(3)=60.
   lat(4)=45.
   lat(5)=30.
   lat(6)=15.
   lat(7)=0.
   lat(8)=-15.
   lat(9)=-30.
   lat(10)=-45.
   lat(11)=-60.
   lat(12)=-75.
   lat(13)=-90.

   albedo(1)=0.58
   albedo(2)=0.45
   albedo(3)=0.26
   albedo(4)=0.18
   albedo(5)=0.18
   albedo(6)=0.15
   albedo(7)=0.11
   albedo(8)=0.12
   albedo(9)=0.13
   albedo(10)=0.13
   albedo(11)=0.15
   albedo(12)=0.675
   albedo(13)=0.75

   !Start date und -time.
   startyear=2000     !fill in number (integer, year)
   startmonth=1       !fill in number (integer, month)
   startday=1         !fill in number (integer, day)
   starttime=0.       !fill in number (real, time in UT)
         
   !integration time and integration timestep
   integrationtime=438300.   !fill in number (real, length of the modellrun in hours)
   timestep=1.               !fill in number (timestep of the integration in hours)
   tianz = int(integrationtime/timestep)
   
   !vertical grid:
   zanz=18
   z(1)=1.
   z(2)=4.
   z(3)=7.
   z(4)=10.
   z(5)=13.
   z(6)=16.
   z(7)=19.
   z(8)=22.
   z(9)=25.
   z(10)=30.
   z(11)=35.
   z(12)=40.
   z(13)=45.
   z(14)=50.
   z(15)=55.
   z(16)=60.
   z(17)=65.
   z(18)=70.

   !*********
   ! initialise the program run:
   !*********

   !Initialise lookup tables for the photolysis frequencies. The routinen make_JO2 und make_JO3 interpolate the photolysis frequencies
   !from precalculated 4-dimensional lookup tables. These tables initialised here once and are then available throughout the program run.

   !intialise the grid of these tables:
   call makejgrid(jgridz,jgridzanz,jgridsza,jgridszaanz,jgridtopo3,jgridtopo3anz,jgridalbedo,jgridalbedoanz)

   !read JO2_array
   open(10,file='inputs/J-O2.dat',status='old')
   read(10,*) JO2log_array
   close(10)
   
   !read JO3_array
   open(10,file='inputs/J-O3.dat',status='old')
   read(10,*) JO3log_array
   close(10)

   !Read standard atmosphere. Standard values for temperature (t), pressure (p), O2-concentration (o2) and concentration of the 
   !air molecules (air) are read and provided on the vertical grid ( z(1...zanz) ) that has been defined above.
   call readatmosphaere(z, zanz, maxzanz, t, p, o2, air)
   
   !*****************************************
   ! Insert the main body of the program here
   !*****************************************
   do li=1, latanz
      write(*,*) "Latitude:", lat(li), "(Albedo: ", albedo(li), ")"
      ! Set initial values for O, Ox, O3, overhead ozone, datetime
      o(1,:) = 0
      o3(1,:) = 0
      ox(1,:) = 0
      topo3(1,:) = 0
      
      year(1) = startyear
      month(1) = startmonth
      day(1) = startday
      time(1) = starttime
      
      ! Loop over time
      do ti=2, int(integrationtime)
         !write(*,*) "ti:", ti
         ! Calculate current date and time from timestep
         call newdate(year(ti-1),month(ti-1),day(ti-1),time(ti-1),timestep,year(ti),month(ti),day(ti),time(ti))
   
         ! Calculate sza
         sza(ti) = make_sza(lon,lat(li),year(ti),month(ti),day(ti),time(ti))
         ! Loop over level
         do zi=1, zanz
            ! Calculate reaction constants
            KO_O2 = make_KO_O2(t(zi),air(zi))
            KO_O3 = make_KO_O3(t(zi))
   
            ! Calculate photolysis freq
            JO2 = make_JO2(jgridz, jgridzanz, jgridsza, jgridszaanz, jgridtopo3, &
               jgridtopo3anz, jgridalbedo, jgridalbedoanz, JO2log_array, &
               z(zi), sza(ti), topo3(ti-1, zi), albedo(li))
            JO3 = make_JO3(jgridz, jgridzanz, jgridsza, jgridszaanz, jgridtopo3, &
               jgridtopo3anz, jgridalbedo, jgridalbedoanz, JO3log_array, &
               z(zi), sza(ti), topo3(ti-1, zi), albedo(li))

            ! Calculate new Ox conc. at level zi
            ox(ti,zi) = ox(ti-1,zi) + (2*JO2*o2(zi) - 2*KO_O3*o(ti-1,zi)*o3(ti-1,zi)) * timestep*60*60
   
            ! At altitude zi calculate O, O3 from Ox
            o3(ti,zi) = KO_O2 * ox(ti,zi) * o2(zi)/(JO3 + KO_O2 * o2(zi))
            o(ti,zi) = JO3*ox(ti,zi)/(KO_O2*o2(zi) + JO3)
         enddo ! finish level loop
         ! Calculate profile of overhead ozone
         call make_topo3(z,zanz,maxzanz,ti,maxtimeanz,o3,topo3)
      enddo ! finish time loop
   
      !**********
      ! write program output
      !**********

      !one file per vertical level: annually on 1st January at 12:00 UT
   do zi=1,zanz
      write(filename,fmt='(i4.4)') int(z(zi))
      write(filespec_albedo,fmt='(f0.2)') albedo(li)
      if (lat(li).lt.0) then
        write(filename_lat,fmt='(i2.2)') abs(int(lat(li)))
        open(10,file='results/lat'//filename_lat(1:2)//'n_km'//filename(3:4)//'.dat')
        write(10,fmt='(a14)') '-'//filename_lat(1:2)//'° latitude'
      else
        write(filename_lat,fmt='(i2.2)') int(lat(li))
        open(10,file='results/lat'//filename_lat(1:2)//'_km'//filename(3:4)//'.dat')
        write(10,fmt='(a13)') filename_lat(1:2)//'° latitude'
      endif
      write(10,fmt='(a15)') filespec_albedo//' albedo'
      write(10,fmt='(a13)') filename(3:4)//'km altitude'
      write(10,fmt='(a195)') '  alt.[km]     year     month       day      hour  SZA[deg]         '//&
       &'o2-vmr         ox-vmr          o-vmr         o3-vmr    '//&
       &'[o2][cm^-3]   [ox][1/cm^3]    [o][1/cm^3]   [o3][1/cm^3]      topo3[DU]'
      do ti=1,tianz
         if ((month(ti).eq.1).and.(day(ti).eq.1).and.(abs(time(ti)-12.).le.0.1)) then
            write(10,fmt='(f9.0,1x,i9,1x,i9,1x,i9,1x,f9.2,1x,f9.2,9(1x,E14.4))')&
             &z(zi),year(ti),month(ti),day(ti),time(ti),sza(ti),&
             &o2(zi)/air(zi),ox(ti,zi)/air(zi),o(ti,zi)/air(zi),o3(ti,zi)/air(zi),&
             &o2(zi),ox(ti,zi),o(ti,zi),o3(ti,zi),topo3(ti,zi)
         endif
      enddo
   enddo
   close(10)

   !one vertical profile per year in separate files: annually on 1st January at 12:00 UT
   do ti=1,tianz
       if ((month(ti).eq.1).and.(day(ti).eq.1).and.(abs(time(ti)-12.).le.0.1)) then
         write(filename,fmt='(i4.4)') year(ti)
         write(filespec_albedo,fmt='(f0.2)') albedo(li)
         if (lat(li).lt.0) then
            write(filename_lat,fmt='(i2.2)') abs(int(lat(li)))
            open(10,file='results/lat'//filename_lat(1:2)//'n_year'//filename//'.dat')
            write(10,fmt='(a14)') '-'//filename_lat(1:2)//'° latitude'
         else
            write(filename_lat,fmt='(i2.2)') int(lat(li))
            open(10,file='results/lat'//filename_lat(1:2)//'_year'//filename//'.dat')
            write(10,fmt='(a13)') filename_lat(1:2)//'° latitude'
         endif
         write(10,fmt='(a15)') filespec_albedo//' albedo'
         write(10,fmt='(a9)') 'Year '//filename
         write(10,fmt='(a195)') '  alt.[km]     year     month       day      hour  SZA[deg]         '//&
          &'o2-vmr         ox-vmr          o-vmr         o3-vmr    '//&
          &'[o2][cm^-3]   [ox][1/cm^3]    [o][1/cm^3]   [o3][1/cm^3]      topo3[DU]'
         do zi=1,zanz
            write(10,fmt='(f9.0,1x,i9,1x,i9,1x,i9,1x,f9.2,1x,f9.2,9(1x,E14.4))')&
             &z(zi),year(ti),month(ti),day(ti),time(ti),sza(ti),&
             &o2(zi)/air(zi),ox(ti,zi)/air(zi),o(ti,zi)/air(zi),o3(ti,zi)/air(zi),&
             &o2(zi),ox(ti,zi),o(ti,zi),o3(ti,zi),topo3(ti,zi)
         enddo
       endif
   enddo
   close(10)

   !for the last year for all levels: a file with data on the first of each month at 12:00UT  
   do zi=1,zanz
      write(filename,fmt='(i4.4)') int(z(zi))
      write(filespec_albedo,fmt='(f0.2)') albedo(li)
      if (lat(li).lt.0) then
         write(filename_lat,fmt='(i2.2)') abs(int(lat(li)))
         open(10,file='results/lat'//filename_lat(1:2)//'n_annual_cycle'//filename(3:4)//'km.dat')
         write(10,fmt='(a14)') '-'//filename_lat(1:2)//'° latitude'
      else
         write(filename_lat,fmt='(i2.2)') int(lat(li))
         open(10,file='results/lat'//filename_lat(1:2)//'_annual_cycle'//filename(3:4)//'km.dat')
         write(10,fmt='(a13)') filename_lat(1:2)//'° latitude'
      endif
      write(10,fmt='(a15)') filespec_albedo//' albedo'
      write(10,fmt='(a13)') filename(3:4)//'km altitude'
      write(10,fmt='(a195)') '  alt.[km]     year     month       day      hour  SZA[deg]         '//&
       &'o2-vmr         ox-vmr          o-vmr         o3-vmr    '//&
       &'[o2][cm^-3]   [ox][1/cm^3]    [o][1/cm^3]   [o3][1/cm^3]      topo3[DU]'
      do ti=1,tianz
         if ((year(ti).eq.year(tianz-int(24*365/timestep))).and.(day(ti).eq.1).and.(abs(time(ti)-12.).le.0.1)) then
            write(10,fmt='(f9.0,1x,i9,1x,i9,1x,i9,1x,f9.2,1x,f9.2,9(1x,E14.4))')&
             &z(zi),year(ti),month(ti),day(ti),time(ti),sza(ti),&
             &o2(zi)/air(zi),ox(ti,zi)/air(zi),o(ti,zi)/air(zi),o3(ti,zi)/air(zi),&
             &o2(zi),ox(ti,zi),o(ti,zi),o3(ti,zi),topo3(ti,zi)
         endif
      enddo
   enddo
   close(10)

   !on the mid-summer day of the last year: hourly data (one file per level)
   do zi=1,zanz
      write(filename,fmt='(i4.4)') int(z(zi))
      write(filespec_albedo,fmt='(f0.2)') albedo(li)
      if (lat(li).lt.0) then
         write(filename_lat,fmt='(i2.2)') abs(int(lat(li)))
         open(10,file='results/lat'//filename_lat(1:2)//'n_diurnal_cycle'//filename(3:4)//'km.dat')
         write(10,fmt='(a14)') '-'//filename_lat(1:2)//'° latitude'
      else
         write(filename_lat,fmt='(i2.2)') int(lat(li))
         open(10,file='results/lat'//filename_lat(1:2)//'_diurnal_cycle'//filename(3:4)//'km.dat')
         write(10,fmt='(a13)') filename_lat(1:2)//'° latitude'
      endif
      write(10,fmt='(a15)') filespec_albedo//' albedo'
      write(10,fmt='(a13)') filename(3:4)//'km altitude'
      write(10,fmt='(a195)') '  alt.[km]     year     month       day      hour  SZA[deg]         '//&
       &'o2-vmr         ox-vmr          o-vmr         o3-vmr    '//&
       &'[o2][cm^-3]   [ox][1/cm^3]    [o][1/cm^3]   [o3][1/cm^3]      topo3[DU]'
      do ti=1,tianz
         if ((year(ti).eq.year(tianz-int(24*365/timestep))).and.(month(ti).eq.6).and.(day(ti).eq.21)) then
            write(10,fmt='(f9.0,1x,i9,1x,i9,1x,i9,1x,f9.2,1x,f9.2,9(1x,E14.4))')&
            & z(zi),year(ti),month(ti),day(ti),time(ti),sza(ti),&
            &o2(zi)/air(zi),ox(ti,zi)/air(zi),o(ti,zi)/air(zi),o3(ti,zi)/air(zi),&
            &o2(zi),ox(ti,zi),o(ti,zi),o3(ti,zi),topo3(ti,zi)
         endif
      enddo
   enddo
   close(10)

   enddo ! finish latitude loop
   write(*,*) "Finished!"
end ! finish program
               
!***********************************************************************
function make_KO_O2(t,air)
!***********************************************************************
   implicit none
   real t, air, make_KO_O2
   
   !*****************************************
   ! Insert code to calculate KO_O2
   !*****************************************
   make_KO_O2 = 6.0e-34*((t/300.)**(-2.4)) * air
end
    
!***********************************************************************
function make_KO_O3(t)
!***********************************************************************
   implicit none
   real t, make_KO_O3
   
   !*****************************************
   ! Insert code to calculate KO_O3
   !*****************************************
    make_KO_O3 = 8.0e-12*exp(-2060./t)
end

!***********************************************************************
subroutine readatmosphaere(z, zanz, maxzanz, t, p, o2, air)
!***********************************************************************
!   z ...........Height      [km]
!   zscal .......Scale heiht [cm]
!   t ...........Temperature [K]
!   p ...........Pressure    [hPa]
!   o2 ..........O2    concentration [molecules / cm**3]
!   n2 ..........N2    concentration [molecules / cm**3]
!   air .........Air   concentration [molecules / cm**3]
!   o ...........O     concentration [molecules / cm**3]
!   o1d .........O(1D) concentration [molecules / cm**3]
!   o3 ..........O3    concentration [molecules / cm**3]
!   topo3 .......overhead O3 column  [DU]

   implicit none
   real stz(200), stzscal(200), stt(200)
   real stp(200), sto2(200), stn2(200)
   real stair(200), sto(200)
   real sto1d(200), sto3(200), sttopo3(200)
   integer sti, stanz

   real zscal(maxzanz), t(maxzanz)
   real p(maxzanz), o2(maxzanz), n2(maxzanz)
   real air(maxzanz), o(maxzanz)
   real o1d(maxzanz), o3(maxzanz), topo3(maxzanz)

   real z(maxzanz)
   integer zi, zanz, maxzanz

   stanz=120

   open(10,file='inputs/standard_atmos.dat', status='old')
   read(10,*)
   read(10,*)
   do sti = 1,stanz
      read(10,*) stz(sti), stzscal(sti), stt(sti), stp(sti), &
     &    sto2(sti), stn2(sti), stair(sti), sto(sti), &
     &    sto1d(sti), sto3(sti)     
   enddo
   close(10)
   
   sttopo3(stanz)=sto3(stanz)*&
     & stzscal(stanz)*28.97/48.
   do sti = stanz-1,1,-1
      sttopo3(sti) = &
     &    sttopo3(sti+1) + &
     &     (stz(sti+1)-stz(sti))&
     &     *0.5*(sto3(sti)+sto3(sti+1))&
     &     *1.e+5
   enddo
   do sti = 1, stanz
      sttopo3(sti)=sttopo3(sti)/2.69e+16
   enddo

   !Auf die Modellevel umsetzen      
   sti=0
   do zi=1,zanz
100      continue
         sti=sti+1
         if (sti.gt.stanz) then
            write(*,*) 'Modelllevel ist nicht in der '//&
     &          'Standardatmosphaere vorhanden !'
            stop
         endif
      if (abs(stz(sti)-z(zi)).gt.0.01) goto 100
      zscal(zi)=stzscal(sti)
      t(zi)=stt(sti)
      p(zi)=stp(sti)
      o2(zi)=sto2(sti)
      n2(zi)=stn2(sti)
      air(zi)=stair(sti)
      o(zi)=sto(sti)
      o1d(zi)=sto1d(sti)
      o3(zi)=sto3(sti)
      topo3(zi)=sttopo3(sti)
   enddo
   
end

!***********************************************************************
function make_JO2(jgridz,jgridzanz,jgridsza,jgridszaanz,jgridtopo3,jgridtopo3anz,jgridalbedo,jgridalbedoanz,&
 &JO2log_array,z,sza,topo3,albedo)
!***********************************************************************
   implicit none
   
   real make_JO2
   real JO2log_array(24, 30, 40, 5), JO2log
   real jgridz(24),jgridsza(30),jgridtopo3(40),jgridalbedo(5)
   integer jgridzanz,jgridszaanz,jgridtopo3anz,jgridalbedoanz
   real z,sza,topo3,albedo

   call read_from_matrix_4d_24_30_40_5&
     & (jgridz,jgridzanz,jgridsza,jgridszaanz,&
     & jgridtopo3,jgridtopo3anz,jgridalbedo,jgridalbedoanz,&
     & JO2log_array,z,sza,topo3,albedo,JO2log)


   if (JO2log.gt.-20.) then 
      make_JO2=10.**JO2log
   else
      make_JO2=0.
   endif
   
end

!***********************************************************************
function make_JO3(jgridz,jgridzanz,jgridsza,jgridszaanz,jgridtopo3,jgridtopo3anz,jgridalbedo,jgridalbedoanz,&
 &JO3log_array,z,sza,topo3,albedo)
!***********************************************************************
   implicit none
   
   real make_JO3
   real JO3log_array(24, 30, 40, 5), JO3log
   real jgridz(24),jgridsza(30),jgridtopo3(40),jgridalbedo(5)
   integer jgridzanz,jgridszaanz,jgridtopo3anz,jgridalbedoanz
   real z,sza,topo3,albedo
 
   call read_from_matrix_4d_24_30_40_5&
     & (jgridz,jgridzanz,jgridsza,jgridszaanz,&
     & jgridtopo3,jgridtopo3anz,jgridalbedo,jgridalbedoanz,&
     & JO3log_array,z,sza,topo3,albedo,JO3log)

   if (JO3log.gt.-20.) then 
      make_JO3   = 10.**JO3log
   else
      make_JO3   = 0.
   endif      

end

!***********************************************************************
function make_sza(lon,lat,year,month,day,time)
!***********************************************************************
   implicit none

   real make_sza
   real lat,lon
   real beta, lambda
   integer year,month,day
   real time
   integer it(5)
   real*8 tis
   real*8 also,deso
   real sonndekl, sonnrekta
   real posdekl, posrekta
   real*8 thetg,tu
   real greenwichrekta
   real great_circle_distance
   real twopi,pi,rad,erdr
   real winkelkorrrad,winkelkorrdeg

   twopi =  2.*4.*atan(1.)
   pi = 4.*atan(1.)
   rad   =  twopi / 360.
   erdr=6371.

   beta = lat * rad
   lambda  = lon * rad

   it(1)=year
   it(2)=month
   it(3)=day
   it(4)=int(time)
   it(5)=int((time-real(int(time)))*60.+0.5)
   tis=0.0d+00
   
   call zjulia (it, tis, thetg, tu)
   greenwichrekta=winkelkorrrad(sngl(thetg))
   call sunpos (tu, also, deso)

   sonnrekta=winkelkorrdeg(sngl(also)/rad)
   sonndekl=winkelkorrdeg(sngl(deso)/rad)

   posrekta=winkelkorrdeg((lambda+greenwichrekta)/rad)
   posdekl=lat
   
   if (sonndekl.ge.270.) then
      sonndekl=sonndekl-360.
   else if (sonndekl.ge.90.) then
      sonnrekta=winkelkorrrad(sonnrekta+180.)
      sonndekl=180.-sonndekl
   end if 

   make_sza=great_circle_distance(posdekl,posrekta,sonndekl,sonnrekta)/erdr/rad

end

!***********************************************************************
subroutine make_topo3(z,zanz,maxzanz,ti,maxtimeanz,o3,topo3)
!***********************************************************************
   implicit none
   
   integer maxzanz,zanz, zi
   integer*8 ti
   integer maxtimeanz
   real z(maxzanz), o3(maxtimeanz,maxzanz)
   real zscal
   real topo3(maxtimeanz,maxzanz)

   !Skalenhoehe (in cm) in ~100km Hoehe
   zscal=6.5e+5
   
   topo3(ti,zanz)=o3(ti,zanz)*zscal
   !*28.97/48.
   do zi = zanz-1,1,-1
      topo3(ti,zi)=topo3(ti,zi+1)+&
     &    (z(zi+1)-z(zi))*0.5*(o3(ti,zi)+o3(ti,zi+1))*1.e+5
   enddo

   do zi = 1,zanz
      topo3(ti,zi)=topo3(ti,zi)/2.69e+16
   enddo
   
end

!***********************************************************************
subroutine newdate(oldyear,oldmonth,oldday,oldtime,increment,year,month,day,time)
!***********************************************************************
   implicit none

   integer oldyear,oldmonth,oldday
   real oldtime
   integer year,month,day
   real time
   real increment,realyear
   integer febdayanz

   year=oldyear
   realyear=real(year)
   month=oldmonth
   day=oldday
   time=oldtime+increment

10    continue
   if(time.ge.24.) then
      time=time-24.
      day=day+1
  
      febdayanz=28
      if ( mod(realyear,4.).lt.0.0001 ) febdayanz=29
      if ( mod(realyear,100.).lt.0.0001 ) febdayanz=28
      if ( mod(realyear,1000.).lt.0.0001 ) febdayanz=29
   
      if ((month.eq.1).or.(month.eq.3).or.(month.eq.5)&
     &    .or.(month.eq.7).or.(month.eq.8).or.(month.eq.10)) then
         if (day.gt.31) then
            month=month+1
            day=1
         endif
      elseif ((month.eq.4).or.(month.eq.6).or.(month.eq.9).or.&
     &   (month.eq.11)) then
         if (day.gt.30) then
            month=month+1
            day=1
         endif
      elseif (month.eq.2) then
         if (day.gt.febdayanz) then
            month=month+1
            day=1
         endif
      elseif (month.eq.12) then
         if (day.gt.31) then
            year=year+1
            month=1
            day=1
         endif
      endif
   endif
   if (time.ge.24.) goto 10
end


!***********************************************************************
!***********************************************************************
!***********************************************************************
!***********************************************************************
!***********************************************************************
! From here on the code contains additional subroutines, which are
! called by the subroutines above or which are needed for
! them. You don't need to deal with the following code at all.
!***********************************************************************
!***********************************************************************
!***********************************************************************
!***********************************************************************
!***********************************************************************

!***********************************************************************
subroutine makejgrid(jgridz,jgridzanz,jgridsza,jgridszaanz,jgridtopo3,jgridtopo3anz,jgridalbedo,jgridalbedoanz)
!***********************************************************************
   implicit none
   real jgridz(24),jgridsza(30),jgridtopo3(40),jgridalbedo(5)
   integer jgridzanz,jgridszaanz,jgridtopo3anz,jgridalbedoanz
   real jgridtopo3tmp
   integer jgridtopo3i
 
   !jgridz=geom. Hoehe in km
   jgridz(1)=1.
   jgridz(2)=4.
   jgridz(3)=7.
   jgridz(4)=10.
   jgridz(5)=13.
   jgridz(6)=16.
   jgridz(7)=19.
   jgridz(8)=22.
   jgridz(9)=25.
   jgridz(10)=30.
   jgridz(11)=35.
   jgridz(12)=40.
   jgridz(13)=45.
   jgridz(14)=50.
   jgridz(15)=55.
   jgridz(16)=60.
   jgridz(17)=65.
   jgridz(18)=70.
   jgridz(19)=75.
   jgridz(20)=80.
   jgridz(21)=85.
   jgridz(22)=90.
   jgridz(23)=95.
   jgridz(24)=100.
   jgridzanz=24
 
   !jgridsza
   jgridsza(1)=0.
   jgridsza(2)=20.
   jgridsza(3)=30.
   jgridsza(4)=40.
   jgridsza(5)=45.
   jgridsza(6)=50.
   jgridsza(7)=55.
   jgridsza(8)=60.
   jgridsza(9)=65.
   jgridsza(10)=67.5
   jgridsza(11)=70.
   jgridsza(12)=72.5
   jgridsza(13)=75.
   jgridsza(14)=77.5
   jgridsza(15)=80.
   jgridsza(16)=82.
   jgridsza(17)=84.
   jgridsza(18)=85.
   jgridsza(19)=86.
   jgridsza(20)=87.
   jgridsza(21)=88.
   jgridsza(22)=89.
   jgridsza(23)=90.
   jgridsza(24)=91.
   jgridsza(25)=92.
   jgridsza(26)=93.
   jgridsza(27)=94.
   jgridsza(28)=95.
   jgridsza(29)=96.
   jgridsza(30)=180.
   jgridszaanz=30

 
   !jgridtopo3=overhead ozone in DU
   jgridtopo3i=0

   jgridtopo3tmp = 0.0

   do while ( jgridtopo3tmp .le.1.0 ) 
      jgridtopo3i=jgridtopo3i+1
      jgridtopo3(jgridtopo3i)=jgridtopo3tmp
      jgridtopo3tmp=jgridtopo3tmp + 0.5
   enddo

   jgridtopo3tmp = 2.0
   do while ( jgridtopo3tmp .le.10.0 ) 
      jgridtopo3i=jgridtopo3i+1
      jgridtopo3(jgridtopo3i)=jgridtopo3tmp
      jgridtopo3tmp=jgridtopo3tmp + 1.0
   enddo

   jgridtopo3tmp = 15.0
   do while ( jgridtopo3tmp .le.50.0 ) 
      jgridtopo3i=jgridtopo3i+1
      jgridtopo3(jgridtopo3i)=jgridtopo3tmp
      jgridtopo3tmp=jgridtopo3tmp + 5.0
   enddo
   jgridtopo3tmp = 60.0
   do while ( jgridtopo3tmp .le.100.0 ) 
      jgridtopo3i=jgridtopo3i+1
      jgridtopo3(jgridtopo3i)=jgridtopo3tmp
      jgridtopo3tmp=jgridtopo3tmp + 10.0
   enddo
   jgridtopo3tmp = 120.0
   do while ( jgridtopo3tmp .le.200.0 ) 
      jgridtopo3i=jgridtopo3i+1
      jgridtopo3(jgridtopo3i)=jgridtopo3tmp
      jgridtopo3tmp=jgridtopo3tmp + 20.0
   enddo
   jgridtopo3tmp = 225.0
   do while ( jgridtopo3tmp .le.300.0 ) 
      jgridtopo3i=jgridtopo3i+1
      jgridtopo3(jgridtopo3i)=jgridtopo3tmp
      jgridtopo3tmp=jgridtopo3tmp + 25.0
   enddo
   jgridtopo3tmp = 350.0
   do while ( jgridtopo3tmp .le.600.0 ) 
      jgridtopo3i=jgridtopo3i+1
      jgridtopo3(jgridtopo3i)=jgridtopo3tmp
      jgridtopo3tmp=jgridtopo3tmp + 50.0
   enddo
   jgridtopo3anz=jgridtopo3i
   
   !jgridalbedo=Albedo
   jgridalbedo(1)=.05
   jgridalbedo(2)= .1
   jgridalbedo(3)= .2
   jgridalbedo(4)= .5
   jgridalbedo(5)= .75
   jgridalbedoanz=5
   
end

!***********************************************************************
subroutine zjulia (it,tis,thetg,tu)
!***********************************************************************
   !calculates the julian date and the and die right ascension of the Greenwich meridian
   !                                                                 
   !input:                                                           
   !      it(1-5) universal time                                     
   !      it(1)   year          (integer)                            
   !      it(2)   month         (integer)                            
   !      it(3)   day           (integer)                            
   !      it(4)   hours         (integer)                            
   !      it(5)   minutes       (integer)                            
   !      tis     seconds       (real)                               
   !                                                                 
   ! output:                                                         
   !      thetg   position of the Grenwich meridian at the input date/time [rad]     
   !      tu      julianian century                           
   !***********************************************************************

   implicit real*8 (a-h,o-z)
   integer it(5)
   real ti(5)
   integer i


   twopi =  2.*4.*atan(1.)
   rad   =  twopi / 360.

   ! it(1)-it(5),tis is the date in year,month,day,hour,minute,seconds

   do 80 i = 1, 5
	 ti(i) = float(it(i))
80    continue

   ! rjdays: calculate the julian date at the beginning of the year (1st Jan 00:00 UT)

   if (it(2) .le. 2) then
	 factor = 365.d0 * ti(1) + ti(3) + 31.d0 * (ti(2) - 1.d0)&
     &   + dint((ti(1) - 1.d0)/4.d0)&
     &   - dint(0.75d0 * (dint((ti(1) - 1.d0)/100.d0) + 1.d0))

   else if (it(2) .ge. 3) then
	 factor = 365.d0 * ti(1) + ti(3) + 31.d0 * (ti(2) - 1.d0)&
     &   - dint(0.4d0 * ti(2) + 2.3d0) + dint(ti(1)/4.d0)&
     &   - dint(0.75d0 * (dint(ti(1)/100.d0) + 1.d0))
   end if
   
   rjdays = 2415020.5d+00 - 693961.d0 + factor

   tu = (rjdays - 2415020.0d0)/36525.d0


   ! right ascension of the greenwich meridian at 00:00 UT

   thetg0 = 99.6909833d0 + 36000.7689d0 * tu + 0.00038708d0 * tu**2
   thetg0 = (thetg0/360.d0 - dint(thetg0/360.d0)) * 360.d0

   ! right ascension of the greenwich meridian at the current time

   thetg = thetg0 + (ti(4)*60.d0 + ti(5) + tis/60.d0)*0.25068447d0

   ! thetg in radian

   thetg = thetg * rad

   ! julianian date
   rjday0 = rjdays + ti(4)/24.d0 + ti(5)/1440.d0 + tis/86371.d0
   ! julian century
   tu = (rjday0 - 2415020.0d0)/36525.d0
   
end
!***********************************************************************
subroutine sunpos (tu, also, deso)
!***********************************************************************
   ! calculates the current position of the sun based on the julian century
   ! the position is returned as right ascension, declination and distance earth-sun
   !
   ! input:                                                               
   !    tu     julian century                         [-]              
   !                                                                 
   ! output:                                                              
   !    also   right acension of the sun              [rad]            
   !    deso   declination of the sun                 [rad]            
   !   (rse    distance earth-sun                     [m])             
   !                                                                 
   ! local parameters:                                                     
   !   exzene   exzentricityof earth's orbit          [-]              
   !   rsem     average distance earth-sun            [km]             
   !   epso     inclination of the ecliptic versus the equatorial plane [grad]           
   !   siepso   sin(epso)                                               
   !   coepso   cos(epso)                                               
!***********************************************************************

! calculation of the right ascension and the declination of the sun and the distance earth-sun

   implicit real*8 (a-h, o-z)
   common /basic/ rad
   pi =  4.*datan(1.0d+00)
   twopi =  2.*4.*datan(1.0d+00)
   rad   =  twopi / 360.0d+00
   
   data exzene/1.675104d-02/, rsem/149.597893d+06/
   data epso/23.452294d0/

! calculation of siepso and coepso
   siepso = dsin (epso*rad)
   coepso = dcos (epso*rad)

   rmlso = 4.88162797d0 + 628.3319509d0 * tu&
     &        + 5.27962099d-06 * tu * tu

   rmaso = 6.25658357d0 + 628.3019457d0 * tu - 2.61799388d-06 * tu * tu - 5.75958653d-08 * tu * tu * tu

   sinma = dsin(rmaso)
   rmaso2 = 2.0d0 * rmaso
   sinma2 = dsin(rmaso2)
   rmaso3 = 3.0d0 * rmaso
   sinma3 = dsin(rmaso3)

   dmlso1 = 3.35007223d-02 - 8.35838179d-5 * tu - 2.44346095d-07 * tu * tu
   dmlso2 = 3.5070646d-04 - 1.74532925d-06 * tu

   dmlso = dmlso1 * sinma + dmlso2 * sinma2 + 5.11381471d-06 * sinma3

   rlaso = rmlso + dmlso
   rlaso = rlaso - twopi * dint(rlaso/twopi)

! right ascension of the sun [rad]

   crlaso = dcos(rlaso)
   talaso = dtan(rlaso)

   also1 = talaso * coepso
   also = datan(also1)

   if (crlaso .lt. 0.0d0) also = also + pi
   if (also .lt. 0.0d0) also = also + twopi

! declinateion of the sun [rad]

   srlaso = dsin(rlaso)
   deso1 = srlaso * siepso

   deso = dasin(deso1)

! distance earth-sun [m]

   rse1 = 1.0000002d0 * (1.0d0 - exzene * exzene)
   rse2 = rmaso + dmlso
   rse3 = 1.0d0 + exzene * dcos(rse2)

   rse = rsem * rse1/rse3 * 1000.0d0

end

!***********************************************************************
subroutine read_from_matrix_4d_24_30_40_5(tabt,tabtanz,tabx,tabxanz,taby,tabyanz,tabz,tabzanz,tabelle,t,x,y,z,werttxyz)
!***********************************************************************
!  no extrapolation => t, x, y, z must be within the range of the matrix!

   implicit none

   real tabelle(24,30,40,5)
   real tabt(24),tabx(30),taby(40),tabz(5)
   integer tabtanz,tabxanz,tabyanz,tabzanz
   real t,x,y,z
   real wert000z,wert010z,wert011z,wert001z,wert0x0z,wert0x1z
   real wert0xyz
   real wert100z,wert110z,wert111z,wert101z,wert1x0z,wert1x1z
   real wert1xyz
   real werttxyz
   
   integer it,ix,iy,iz
   integer it0,it1,ix0,ix1,iy0,iy1,iz0,iz1

   real interpol 
   
   if ( (t.lt.tabt(1)).or.(t.gt.tabt(tabtanz)).or.&
     &     (x.lt.tabx(1)).or.(x.gt.tabx(tabxanz)).or.&
     &     (y.lt.taby(1)).or.(y.gt.taby(tabyanz)).or.&
     &     (z.lt.tabz(1)).or.(z.gt.tabz(tabzanz)) ) then
       write(*,*) t, tabt(1), tabt(tabtanz)
       write(*,*) x, tabx(1), tabx(tabxanz)
       write(*,*) y, taby(1), taby(tabyanz)
       write(*,*) z, tabz(1), tabz(tabzanz)

       write(*,*) 'read_from_matrix_4d: out of range of lookup table!'
      stop
   end if
   
   if (t.eq.tabt(1)) then
      it0=1
      it1=2
   else
      it=0
5        continue
         it=it+1
      if (t.gt.tabt(it)) goto 5
      it0=it-1
      it1=it
   end if 

   if (x.eq.tabx(1)) then
      ix0=1
      ix1=2
   else
      ix=0
10       continue
         ix=ix+1
      if (x.gt.tabx(ix)) goto 10
      ix0=ix-1
      ix1=ix
   end if

   if (y.eq.taby(1)) then
      iy0=1
      iy1=2
   else
      iy=0
15       continue
         iy=iy+1
      if (y.gt.taby(iy)) goto 15
      iy0=iy-1
      iy1=iy
   end if

   if (z.eq.tabz(1)) then
      iz0=1
      iz1=2
   else
      iz=0
20       continue
         iz=iz+1
      if (z.gt.tabz(iz)) goto 20
      iz0=iz-1
      iz1=iz
   end if
   
   ! write(*,*) 'pos:', it0, ix0, iy0, iz0

   wert000z=interpol(real(tabelle(it0,ix0,iy0,iz0)),&
     & real(tabelle(it0,ix0,iy0,iz1)),tabz(iz0),tabz(iz1),z)
   wert010z=interpol(real(tabelle(it0,ix1,iy0,iz0)),&
     & real(tabelle(it0,ix1,iy0,iz1)),tabz(iz0),tabz(iz1),z)
   wert011z=interpol(real(tabelle(it0,ix1,iy1,iz0)),&
     & real(tabelle(it0,ix1,iy1,iz1)),tabz(iz0),tabz(iz1),z)
   wert001z=interpol(real(tabelle(it0,ix0,iy1,iz0)),&
     & real(tabelle(it0,ix0,iy1,iz1)),tabz(iz0),tabz(iz1),z)

   wert100z=interpol(real(tabelle(it1,ix0,iy0,iz0)),&
     & real(tabelle(it1,ix0,iy0,iz1)),tabz(iz0),tabz(iz1),z)
   wert110z=interpol(real(tabelle(it1,ix1,iy0,iz0)),&
     & real(tabelle(it1,ix1,iy0,iz1)),tabz(iz0),tabz(iz1),z)
   wert111z=interpol(real(tabelle(it1,ix1,iy1,iz0)),&
     & real(tabelle(it1,ix1,iy1,iz1)),tabz(iz0),tabz(iz1),z)
   wert101z=interpol(real(tabelle(it1,ix0,iy1,iz0)),&
     & real(tabelle(it1,ix0,iy1,iz1)),tabz(iz0),tabz(iz1),z)
     
   ! write(*,*) 'werte0:', wert000z, wert010z, wert011z, wert001z
   ! write(*,*) 'werte1:', wert100z, wert110z, wert111z, wert101z

   wert0x0z=interpol(wert000z,wert010z,tabx(ix0),tabx(ix1),x)
   wert0x1z=interpol(wert001z,wert011z,tabx(ix0),tabx(ix1),x)

   wert1x0z=interpol(wert100z,wert110z,tabx(ix0),tabx(ix1),x)
   wert1x1z=interpol(wert101z,wert111z,tabx(ix0),tabx(ix1),x)


   wert0xyz=interpol(wert0x0z,wert0x1z,taby(iy0),taby(iy1),y)               
   wert1xyz=interpol(wert1x0z,wert1x1z,taby(iy0),taby(iy1),y)  
   
   
   werttxyz=interpol(wert0xyz,wert1xyz,tabt(it0),tabt(it1),t)
         
          
end

!***********************************************************************
function interpol(werta,wertb,xa,xb,xinter)
!***********************************************************************

   real interpol,werta,wertb,xa,xb,xinter

   ! if ((werta.lt.-900.).or.(wertb.lt.-900.).
   !c   or.(xa.lt.-900.).or.(xb.lt.-900.)) then
   !     interpol=-999.999 
   !     return
   ! end if

   interpol=werta + (wertb-werta)/(xb-xa)*(xinter-xa)

end 

!***********************************************************************
function winkelkorrrad(winkel)
!***********************************************************************
   ! moves angles in rad into the interval [0,2*pi)

    real winkelkorrrad,winkel,pi
    pi=3.1416

    winkelkorrrad=winkel

10     if (winkelkorrrad.lt.0.) then 
       winkelkorrrad=winkelkorrrad+2*pi
       goto 10
    end if

20     if (winkelkorrrad.ge.(2*pi)) then
       winkelkorrrad=winkelkorrrad-2*pi
       goto 20
    end if

end
!***********************************************************************
function winkelkorrdeg(winkel)
!***********************************************************************
   ! moves angles in deg into the interval [0,360)

    real winkelkorrdeg,winkel

    winkelkorrdeg=winkel

10     if (winkelkorrdeg.lt.0.) then 
       winkelkorrdeg=winkelkorrdeg+360.
       goto 10
    end if

20     if (winkelkorrdeg.ge.360.) then
       winkelkorrdeg=winkelkorrdeg-360.
       goto 20
    end if

end

!***********************************************************************
function great_circle_distance(lat1,lon1,lat2,lon2)
!***********************************************************************
   real great_circle_distance
   real lat1
   real lon1
   real lat2
   real lon2
   real hilfwinkel

   real winkelkorrrad,winkelkorrdeg

   real cospsi,erdr,psi
   real pi,rr

   pi=3.1416
   rr=2.*pi/360.
   erdr=6371.
  
   hilfwinkel=winkelkorrdeg(lon1-lon2)
   cospsi=(cos(rr*lat1)*cos(rr*lat2)*cos(rr*hilfwinkel)+&
     &    sin(rr*lat1)*sin(rr*lat2))

   cospsi=min(cospsi,1.)
   cospsi=max(cospsi,-1.)
   psi=acos(cospsi)
   great_circle_distance=winkelkorrrad(psi)*erdr

end
