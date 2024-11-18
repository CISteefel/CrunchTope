      subroutine rdpar(adh,adhh,adhv,aphi,bdh,bdhh,bdhv,bdot,bdoth,
     $ bdotv,cco2,dadhh,dadhv,dbdhh,dbdhv,dbdth,dbdtv,dhfe,dvfe,
     $ ipch,ipchmx,ipcv,ipcvmx,itgenf,nacdpr,narxmx,narxt,ndat0s,
     $ ndbmax,ndbptg,ndbptl,nerr,noutpt,ntprmx,ntprt,nttyo,nwarn,
     $ prehw,presg,q500fl,tdamax,tdamin,tempc,uakey,udbfmt,udbval,
     $ xdbval,xhfe,xlke,xvfe)
c
c     This subroutine reads the data describing the temperature grid
c     and data for miscellaneous parameters which are represented
c     on this grid from the DATA0 file. Such parameters include the
c     pressure, the Debye-Huckel parameters A(gamma,10), B(gamma),
c     and A(phi), the B-dot parameter of Helgeson (1969), and the
c     equilibrium constant for the "Eh" reaction. Depending on the
c     values of the ipch and ipcv flags, pressure derivatives of
c     various order of these parameters may also be read. If pressure
c     corrections are available, a grid for the half-width of the
c     pressure correction band is also read. This subroutine also
c     reads the nominal temperature limits for the data file. These
c     limits are written on the DATA1 and DATA1F files. This suboutine
c     also reads the coefficients needed to calculate the acivity
c     coefficient of CO2(aq) using the Drummond (1981) equation.
c
c     The counter nerr is incremented by one for each error encountered.
c     The counter nwarn is similarly incremented for each warning.
c
c     This suboutine is called by:
c
c       EQPT/eqpt.f
c
c-----------------------------------------------------------------------
c
c     Principal input:
c
c       ipch   = enthalpy functions data grid flag:
c                  -1 = no enthalpy grids
c                   0 = enthalpy grids present
c                   1 = grids present for the enthalpy and its first
c                         pressure derivative
c                   2 = grids present for the enthalpy and its first
c                         and second pressure derivatives
c       ipcv   = volume functions data flag:
c                  -1 = no volume grids
c                   0 = volume grids present
c                   1 = grids present for the volume and its first
c                         pressure derivative
c                   2 = grids present for the volume and its first
c                         and second pressure derivatives
c       ndat0s = unit number of the stripped DATA0 file
c       narxt  = array of numbers of coefficients in the temperature
c                  ranges
c       ntprt  = the number of temperature ranges on the standard
c                  temperature grid
c       uakey  = string specifying the type of data file ("SEDH" or
c                  "Pitzer") being processed
c
c     Principal output:
c
c       adh    = array of A(gamma,10) values on the standard
c                  temperature grid
c       adhh   = array of A(H) values on the standard temperature
c                  grid; A(H) contains but is not equal to
c                  dA(gamma,10)/dT
c       adhv   = array of A(V) values on the standard temperature
c                  grid; A(V) contains but is not equal to
c                  dA(gamma,10)/dP
c       aphi   = array of A(phi) values on the standard temperature
c                  grid; A(phi) = A(gamma,e)/3 = 2.303A(gamma,10)/3,
c                  so the pressure derivatives of A(phi) are obtained
c                  from A(V) and dnA(V)/dPn.
c       bdh    = array of B(gamma) values on the standard temperature
c                  grid
c       bdhh   = array of B(H) values on the standard temperature grid
c       bdhv   = array of B(V) values on the standard temperature grid
c       bdot   = array of B-dot values on the standard temperature
c                  grid
c       bdoth  = array of B-dot(H)values on the standard temperature
c                  grid
c       bdotv  = array of B-dot(V) values on the standard temperature
c                  grid
c       cco2   = array of coefficients for the Drummond (1981)
c                  equation
c       dadhh  = array of dnA(H)/dPn values on the standard temperature
c                  grid
c       dadhv  = array of dnA(V)/dPn values on the standard temperature
c                  grid
c       dbdhh  = array of dnB(H)/dPn values on the standard temperature
c                  grid
c       dbdhv  = array of dnA(V)/dPn values on the standard temperature
c                  grid
c       dbdth  = array of dnB-dot(H)/dPn values on the standard
c                  temperature grid
c       dbdtv  = array of dnB-dot(V)/dPn values on the standard
c                  temperature grid
c       dhfe   = array of pressure derivatives of the enthalpy of
c                  reaction for the "Eh" reaction on the standard
c                  temperature grid
c       dvfe   = array of pressure derivatives of the volume of
c                  reaction for the "Eh" reaction on the standard
c                  temperature grid
c       ndbptg = the number of distinct points on the standard
c                  temperature grid.
c       nerr   = cumulative error counter
c       nwarn  = cumulative warning counter
c       presg  = array of standard pressures on the standard
c                  temperature grid
c       prehw  = array of half-width values for the standard pressure
c                  envelope on the standard temperature grid
c       tdamax = the nominal maximum temperature (C) for which the
c                  data file is valid
c       tdamin = the nominal minimum temperature (C) for which the
c                  data file is valid
c       tempc  = array of temperatures defining the standard
c                  temperature grid
c       udbfmt = the format used when reading a line of data on the
c                  standard temperature grid
c       xlke   = array of log K values for the "Eh" reaction
c                  on the standard temperature grid
c       xhfe   = array of enthalpy of reaction values for the "Eh"
c                  reaction hon the standard temperature grid
c       xvfe   = array of volume of reaction values for the "Eh"
c                  reaction on the standard temperature grid
c
c-----------------------------------------------------------------------
c
      implicit none
c
c-----------------------------------------------------------------------
c
c     Calling sequence variable declarations.
c
      integer ipchmx,ipcvmx,narxmx,ndbmax,ntprmx
c
      integer ndat0s,noutpt,nttyo
c
      integer nacdpr(ntprmx),narxt(ntprmx)
c
      integer ipch,ipcv,itgenf,ndbptg,ndbptl,nerr,ntprt,nwarn
c
      logical q500fl
c
      character(len=16) udbval(ndbmax)
c
      character(len=16) udbfmt
      character(len=8) uakey
c
      real(8) tdamax,tdamin
c
      real(8) cco2(5)
c
      real(8) adh(narxmx,ntprmx),adhh(narxmx,ntprmx),
     $ adhv(narxmx,ntprmx),aphi(narxmx,ntprmx),bdh(narxmx,ntprmx),
     $ bdhh(narxmx,ntprmx),bdhv(narxmx,ntprmx),bdot(narxmx,ntprmx),
     $ bdoth(narxmx,ntprmx),bdotv(narxmx,ntprmx),
     $ dadhh(narxmx,ntprmx,ipchmx),dadhv(narxmx,ntprmx,ipcvmx),
     $ dbdhh(narxmx,ntprmx,ipchmx),dbdhv(narxmx,ntprmx,ipcvmx),
     $ dbdth(narxmx,ntprmx,ipchmx),dbdtv(narxmx,ntprmx,ipcvmx),
     $ dhfe(narxmx,ntprmx,ipchmx),dvfe(narxmx,ntprmx,ipcvmx),
     $ prehw(narxmx,ntprmx),presg(narxmx,ntprmx),tempc(narxmx,ntprmx),
     $ xdbval(ndbmax),xhfe(narxmx,ntprmx),xlke(narxmx,ntprmx),
     $ xvfe(narxmx,ntprmx)
c
c-----------------------------------------------------------------------
c
c     Local variable declarations.
c
      integer i,ipc,j2,j3,n,ndbloc,ntpr
c
      integer ilnobl
c
      logical qend,qerr,q500nd
c
      character(len=24) ux24,ux24lc
      character(len=80) ulbufa,ulbufb,uline,ustr80
c
c-----------------------------------------------------------------------
c
c     Skip to the 'Temperature limits' line.
c
      ux24 = 'Temperature limits (degC'
      ux24lc = ux24
      call locase(ux24lc)
  100 read (ndat0s,1000,end=110,err=995) ustr80
 1000 format(a)
      call lejust(ustr80)
      call locase(ustr80)
      if (ustr80(1:24) .ne. ux24lc(1:24)) go to 100
      go to 120
c
  110 j2 = ilnobl(ux24)
      write (noutpt,1010) ux24(1:j2)
      write (nttyo,1010) ux24(1:j2)
 1010 format(/' * Error - (EQPT/rdpar) Unexpectedly encountered',
     $ ' end-of-file',/7x,'while searching the DATA0 file for the',
     $ ' find the line containing',/7x,'"',a,'".')
      stop
  120 continue
c
c     Read in the parameters. First, read the nominal temperature
c     limits.
c
      read (ndat0s,1002,end=990,err=995) tdamin,tdamax
 1002 format(5x,2f10.4)
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Read the standard temperature grid. The format of this grid will
c     be analyzed and will serve as a template for that of subsequent
c     data presented on this grid. The original EQ3/6 format was
c     four points per line. Up to six points per line are now possible.
c     The format will be determined by counting the number of data
c     points on the first line.
c
c       ndbptl = maximum number of data points per line (1-6)
c       ndbptg = number of data points in the entire grid
c       ndbloc = number of lines containing the grid
c
      ndbptg = narxt(1)
      do ntpr = 2,ntprt
        ndbptg = ndbptg + narxt(ntpr) - 1
      enddo
c
      ux24 = 'Temperatures            '
      ux24lc = ux24
      call locase(ux24lc)
      read (ndat0s,1000,end=990,err=995) ustr80
      call lejust(ustr80)
      call locase(ustr80)
      if (ustr80(1:24) .ne. ux24lc(1:24)) then
        ux24 = 'Temperature grid (degC) '
        ux24lc = ux24
        call locase(ux24lc)
      endif
      if (ustr80(1:24) .ne. ux24lc(1:24)) then
        ux24 = 'Temperature grid (C)    '
        ux24lc = ux24
        call locase(ux24lc)
      endif
c
      if (ustr80(1:24) .ne. ux24lc(1:24)) then
        j2 = ilnobl(ux24)
        j3 = ilnobl(ustr80)
        j3 = min(j3,24)
        write (noutpt,1020) ustr80(1:j3),ux24(1:j2)
        write (nttyo,1020) ustr80(1:j3),ux24(1:j2)
 1020   format(/' * Error - (EQPT/rdpar) Found the string "',a,'"',
     $  /7x,'on the DATA0 file where the string "',a,'"',
     $  /7x,'was expected.')
        stop
      endif
c
c     Analyze the grid format.
c
      read (ndat0s,1030,end=990,err=995) uline
 1030 format(a80)
      backspace(ndat0s)
c
      ulbufa = uline
      ulbufb = ' '
      ndbptl = 0
      do n = 1,6
        call lejust(ulbufa)
        i = index(ulbufa,' ')
        if (i .le. 1) go to 130
        ndbptl = ndbptl + 1
        ulbufb = ulbufa(i:80)
        ulbufa = ulbufb
      enddo
  130 continue
c
      ndbloc = ndbptg/ndbptl
      if (mod(ndbptg,ndbptl) .gt. 0) ndbloc = ndbloc + 1
      udbfmt = '( (5x,6f10.4) )'
      write (udbfmt(7:7),'(i1)') ndbptl
c
c     The format analysis is complete.
c     Read the temperatures on the log K temperature grid.
c     Return the data in the xdbval holding array.
c
      q500nd = .false.
      call rdgrid(ndat0s,ndbmax,ndbptg,ndbptl,qend,qerr,q500nd,
     $ udbval,xdbval)
      if (qend) go to 990
      if (qerr) go to 995
c     read (ndat0s,udbfmt,end=990,err=995) (xdbval(i), i = 1,ndbptg)
c
c     Load the data into the tempc array.
c
c     Calling sequence substitutions:
c       tempc for zdbval
c
      call ldbar2(nacdpr,narxmx,narxt,ndbmax,ntprmx,ntprt,
     $ xdbval,tempc)
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Read the standard pressure grid.
c
      ux24 = 'Pressures               '
      ux24lc = ux24
      call locase(ux24lc)
      read (ndat0s,1000,end=990,err=995) ustr80
      call lejust(ustr80)
      call locase(ustr80)
      if (ustr80(1:24) .ne. ux24lc(1:24)) then
        ux24 = 'Pressure grid (bars)    '
        ux24lc = ux24
        call locase(ux24lc)
      endif
c
      if (ustr80(1:24) .ne. ux24lc(1:24)) then
        j2 = ilnobl(ux24)
        j3 = ilnobl(ustr80)
        j3 = min(j3,24)
        write (noutpt,1020) ustr80(1:j3),ux24(1:j2)
        write (nttyo,1020) ustr80(1:j3),ux24(1:j2)
        stop
      endif
c
c     Read the pressures on the log K temperature grid.
c     Return the data in the xdbval holding array.
c
      q500nd = .false.
      call rdgrid(ndat0s,ndbmax,ndbptg,ndbptl,qend,qerr,q500nd,
     $ udbval,xdbval)
      if (qend) go to 990
      if (qerr) go to 995
c
c     Load the data into the presg array.
c
c     Calling sequence substitutions:
c       presg for zdbval
c
      call ldbar2(nacdpr,narxmx,narxt,ndbmax,ntprmx,ntprt,
     $ xdbval,presg)
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
      if (ipcv .ge. 0) then
c
c       Read the standard pressure envelope half-width.
c
        ux24 = 'Pressure envelope half-w'
        ux24lc = ux24
        call locase(ux24lc)
        read (ndat0s,1000,end=990,err=995) ustr80
        call lejust(ustr80)
        call locase(ustr80)
c
        if (ustr80(1:24) .ne. ux24lc(1:24)) then
          j2 = ilnobl(ux24)
          j3 = ilnobl(ustr80)
          j3 = min(j3,24)
          write (noutpt,1020) ustr80(1:j3),ux24(1:j2)
          write (nttyo,1020) ustr80(1:j3),ux24(1:j2)
          stop
        endif
c
c       Read the pressure half-widths on the log K temperature grid.
c       Return the data in the xdbval holding array.
c
        q500nd = .false.
        call rdgrid(ndat0s,ndbmax,ndbptg,ndbptl,qend,qerr,q500nd,
     $  udbval,xdbval)
        if (qend) go to 990
        if (qerr) go to 995
c
c       Load the data into the prehw array.
c
c       Calling sequence substitutions:
c         prehw for zdbval
c
        call ldbar2(nacdpr,narxmx,narxt,ndbmax,ntprmx,ntprt,
     $  xdbval,prehw)
      endif
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
      if (uakey(1:8) .eq .'SEDH    ') then
c
c       Read the grid for the Debye-Huckel A(gamma,10) parameter and
c       the grids for related parameters.
c
        ux24 = 'Debye Huckel A (ADH)    '
        ux24lc = ux24
        call locase(ux24lc)
        read (ndat0s,1000,end=990,err=995) ustr80
        call lejust(ustr80)
        call locase(ustr80)
        if (ustr80(1:24) .ne. ux24lc(1:24)) then
          ux24 = 'Debye-Huckel A_gamma (kg'
          ux24lc = ux24
          call locase(ux24lc)
        endif
        if (ustr80(1:24) .ne. ux24lc(1:24)) then
          ux24 = 'Debye-Huckel A_gamma    '
          ux24lc = ux24
          call locase(ux24lc)
        endif
c
        if (ustr80(1:24) .ne. ux24lc(1:24)) then
          j2 = ilnobl(ux24)
          j3 = ilnobl(ustr80)
          j3 = min(j3,24)
          write (noutpt,1020) ustr80(1:j3),ux24(1:j2)
          write (nttyo,1020) ustr80(1:j3),ux24(1:j2)
          stop
        endif
c
c       Read the Debye-Huckel A(gamma,10) values on the log K
c       temperature grid. Return the data in the xdbval holding array.
c
        q500nd = q500fl
        call rdgrid(ndat0s,ndbmax,ndbptg,ndbptl,qend,qerr,q500nd,
     $  udbval,xdbval)
        if (qend) go to 990
        if (qerr) go to 995
c
c       Load the data into the adh array.
c
c       Calling sequence substitutions:
c         adh for zdbval
c
        call ldbar2(nacdpr,narxmx,narxt,ndbmax,ntprmx,ntprt,
     $  xdbval,adh)
c
        if (ipch .ge. 0) then
c
c         Read the grid for the Debye-Huckel A(H) parameter.
c
          ux24 = 'Debye-Huckel A_H (kcal*k'
          ux24lc = ux24
          call locase(ux24lc)
          read (ndat0s,1000,end=990,err=995) ustr80
          call lejust(ustr80)
          call locase(ustr80)
          if (ustr80(1:24) .ne. ux24lc(1:24)) then
            ux24 = 'Debye-Huckel A_H        '
            ux24lc = ux24
            call locase(ux24lc)
          endif
c
          if (ustr80(1:24) .ne. ux24lc(1:24)) then
            j2 = ilnobl(ux24)
            j3 = ilnobl(ustr80)
            j3 = min(j3,24)
            write (noutpt,1020) ustr80(1:j3),ux24(1:j2)
            write (nttyo,1020) ustr80(1:j3),ux24(1:j2)
            stop
          endif
c
c         Read the Debye-Huckel A(H) values on the log K temperature
c         grid. Return the data in the xdbval holding array.
c
          q500nd = q500fl
          call rdgrid(ndat0s,ndbmax,ndbptg,ndbptl,qend,qerr,q500nd,
     $    udbval,xdbval)
          if (qend) go to 990
          if (qerr) go to 995
c
c         Load the data into the adhh array.
c
c         Calling sequence substitutions:
c           adhh for zdbval
c
          call ldbar2(nacdpr,narxmx,narxt,ndbmax,ntprmx,ntprt,
     $    xdbval,adhh)
c
          do ipc = 1,ipch
c
            read (ndat0s,1000,end=990,err=995) ustr80
c
c           Read the Debye-Huckel A(H) pressure derivatives on the
c           log K temperature grid. Return the data in the xdbval
c           holding array.
c
            q500nd = q500fl
            call rdgrid(ndat0s,ndbmax,ndbptg,ndbptl,qend,qerr,q500nd,
     $      udbval,xdbval)
            if (qend) go to 990
            if (qerr) go to 995
c
c           Load the data into the dadhh array.
c
c           Calling sequence substitutions:
c             dadhh for zdbval
c             ipchmx for ipcmax
c
            call ldbar3(ipc,ipchmx,nacdpr,narxmx,narxt,ndbmax,
     $      ntprmx,ntprt,xdbval,dadhh)
          enddo
        endif
c
        if (ipcv .ge. 0) then
c
c         Read the grid for the Debye-Huckel A(V) parameter.
c
          ux24 = 'Debye-Huckel A_V (cm**3*'
          ux24lc = ux24
          call locase(ux24lc)
          read (ndat0s,1000,end=990,err=995) ustr80
          call lejust(ustr80)
          call locase(ustr80)
          if (ustr80(1:24) .ne. ux24lc(1:24)) then
            ux24 = 'Debye-Huckel A_V        '
            ux24lc = ux24
            call locase(ux24lc)
          endif
c
          if (ustr80(1:24) .ne. ux24lc(1:24)) then
            j2 = ilnobl(ux24)
            j3 = ilnobl(ustr80)
            j3 = min(j3,24)
            write (noutpt,1020) ustr80(1:j3),ux24(1:j2)
            write (nttyo,1020) ustr80(1:j3),ux24(1:j2)
            stop
          endif
c
c         Read the Debye-Huckel A(V) values on the log K temperature
c         grid. Return the data in the xdbval holding array.
c
          q500nd = q500fl
          call rdgrid(ndat0s,ndbmax,ndbptg,ndbptl,qend,qerr,q500nd,
     $    udbval,xdbval)
          if (qend) go to 990
          if (qerr) go to 995
c
c         Load the data into the adhv array.
c
c         Calling sequence substitutions:
c           adhv for zdbval
c
          call ldbar2(nacdpr,narxmx,narxt,ndbmax,ntprmx,ntprt,
     $    xdbval,adhv)
c
          do ipc = 1,ipcv
c
            read (ndat0s,1000,end=990,err=995) ustr80
c
c           Read the Debye-Huckel A(V) pressure derivatives on the
c           log K temperature grid. Return the data in the xdbval
c           holding array.
c
            q500nd = q500fl
            call rdgrid(ndat0s,ndbmax,ndbptg,ndbptl,qend,qerr,q500nd,
     $      udbval,xdbval)
            if (qend) go to 990
            if (qerr) go to 995
c
c           Load the data into the dadhv array.
c
c           Calling sequence substitutions:
c             dadhv for zdbval
c             ipcvmx for ipcmax
c
            call ldbar3(ipc,ipcvmx,nacdpr,narxmx,narxt,ndbmax,
     $      ntprmx,ntprt,xdbval,dadhv)
          enddo
        endif
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c       Read the grid for the Debye-Huckel B(gamma) parameter and
c
        ux24 = 'Debye Huckel B (BDH)    '
        ux24lc = ux24
        call locase(ux24lc)
        read (ndat0s,1000,end=990,err=995) ustr80
        call lejust(ustr80)
        call locase(ustr80)
        if (ustr80(1:24) .ne. ux24lc(1:24)) then
          ux24 = 'Debye-Huckel B_gamma (kg'
          ux24lc = ux24
          call locase(ux24lc)
        endif
        if (ustr80(1:24) .ne. ux24lc(1:24)) then
          ux24 = 'Debye-Huckel B_gamma    '
          ux24lc = ux24
          call locase(ux24lc)
        endif
c
        if (ustr80(1:24) .ne. ux24lc(1:24)) then
          j2 = ilnobl(ux24)
          j3 = ilnobl(ustr80)
          j3 = min(j3,24)
          write (noutpt,1020) ustr80(1:j3),ux24(1:j2)
          write (nttyo,1020) ustr80(1:j3),ux24(1:j2)
          stop
        endif
c
c       Read the Debye-Huckel B(gamma) values on the log K temperature
c       grid. Return the data in the xdbval holding array.
c
        q500nd = q500fl
        call rdgrid(ndat0s,ndbmax,ndbptg,ndbptl,qend,qerr,q500nd,
     $  udbval,xdbval)
        if (qend) go to 990
        if (qerr) go to 995
c
c       Load the data into the bdh array.
c
c       Calling sequence substitutions:
c         bdh for zdbval
c
        call ldbar2(nacdpr,narxmx,narxt,ndbmax,ntprmx,ntprt,
     $  xdbval,bdh)
c
        if (ipch .ge. 0) then
c
c         Read the grid for the Debye-Huckel B(H) parameter.
c
          ux24 = 'Debye-Huckel B_H (cal*kg'
          ux24lc = ux24
          call locase(ux24lc)
          read (ndat0s,1000,end=990,err=995) ustr80
          call lejust(ustr80)
          call locase(ustr80)
          if (ustr80(1:24) .ne. ux24lc(1:24)) then
            ux24 = 'Debye-Huckel B_H        '
            ux24lc = ux24
            call locase(ux24lc)
          endif
c
          if (ustr80(1:24) .ne. ux24lc(1:24)) then
            j2 = ilnobl(ux24)
            j3 = ilnobl(ustr80)
            j3 = min(j3,24)
            write (noutpt,1020) ustr80(1:j3),ux24(1:j2)
            write (nttyo,1020) ustr80(1:j3),ux24(1:j2)
            stop
          endif
c
c         Read the Debye-Huckel B(H) values on the log K temperature
c         grid. Return the data in the xdbval holding array.
c
          q500nd = q500fl
          call rdgrid(ndat0s,ndbmax,ndbptg,ndbptl,qend,qerr,q500nd,
     $    udbval,xdbval)
          if (qend) go to 990
          if (qerr) go to 995
c
c         Load the data into the bdhh array.
c
c         Calling sequence substitutions:
c           bdhh for zdbval
c
          call ldbar2(nacdpr,narxmx,narxt,ndbmax,ntprmx,ntprt,
     $    xdbval,bdhh)
c
          do ipc = 1,ipch
c
            read (ndat0s,1000,end=990,err=995) ustr80
c
c           Read the Debye-Huckel B(H) pressure derivatives on the
c           log K temperature grid. Return the data in the xdbval
c           holding array.
c
            q500nd = q500fl
            call rdgrid(ndat0s,ndbmax,ndbptg,ndbptl,qend,qerr,q500nd,
     $      udbval,xdbval)
            if (qend) go to 990
            if (qerr) go to 995
c
c           Load the data into the dbdhh array.
c
c           Calling sequence substitutions:
c             dbdhh for zdbval
c             ipchmx for ipcmax
c
            call ldbar3(ipc,ipchmx,nacdpr,narxmx,narxt,ndbmax,
     $      ntprmx,ntprt,xdbval,dbdhh)
          enddo
        endif
c
        if (ipcv .ge. 0) then
c
c         Read the grid for the Debye-Huckel B(V) parameter.
c
          ux24 = 'Debye-Huckel B_V (cm**2*'
          ux24lc = ux24
          call locase(ux24lc)
          read (ndat0s,1000,end=990,err=995) ustr80
          call lejust(ustr80)
          call locase(ustr80)
          if (ustr80(1:24) .ne. ux24lc(1:24)) then
            ux24 = 'Debye-Huckel B_V        '
            ux24lc = ux24
            call locase(ux24lc)
          endif
c
          if (ustr80(1:24) .ne. ux24lc(1:24)) then
            j2 = ilnobl(ux24)
            j3 = ilnobl(ustr80)
            j3 = min(j3,24)
            write (noutpt,1020) ustr80(1:j3),ux24(1:j2)
            write (nttyo,1020) ustr80(1:j3),ux24(1:j2)
            stop
          endif
c
c         Read the Debye-Huckel B(V) values on the log K temperature
c         grid. Return the data in the xdbval holding array.
c
          q500nd = q500fl
          call rdgrid(ndat0s,ndbmax,ndbptg,ndbptl,qend,qerr,q500nd,
     $    udbval,xdbval)
          if (qend) go to 990
          if (qerr) go to 995
c
c         Load the data into the bdhv array.
c
c         Calling sequence substitutions:
c           bdhv for zdbval
c
          call ldbar2(nacdpr,narxmx,narxt,ndbmax,ntprmx,ntprt,
     $    xdbval,bdhv)
c
          do ipc = 1,ipcv
c
            read (ndat0s,1000,end=990,err=995) ustr80
c
c           Read the Debye-Huckel B(V) pressure derivatives on the
c           log K temperature grid. Return the data in the xdbval
c           holding array.
c
            q500nd = q500fl
            call rdgrid(ndat0s,ndbmax,ndbptg,ndbptl,qend,qerr,q500nd,
     $      udbval,xdbval)
            if (qend) go to 990
            if (qerr) go to 995
c
c           Load the data into the dbdhv array.
c
c           Calling sequence substitutions:
c             dbdhv for zdbval
c             ipcvmx for ipcmax
c
            call ldbar3(ipc,ipcvmx,nacdpr,narxmx,narxt,ndbmax,
     $      ntprmx,ntprt,xdbval,dbdhv)
          enddo
        endif
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c       Read the grid for the B-dot parameter and the grids for
c       related parameters.
c
        ux24 = 'Bdot                    '
        ux24lc = ux24
        call locase(ux24lc)
        read (ndat0s,1000,end=990,err=995) ustr80
        call lejust(ustr80)
        call locase(ustr80)
        if (ustr80(1:24) .ne. ux24lc(1:24)) then
          ux24 = 'B-dot                   '
          ux24lc = ux24
          call locase(ux24lc)
        endif
c
        if (ustr80(1:24) .ne. ux24lc(1:24)) then
          j2 = ilnobl(ux24)
          j3 = ilnobl(ustr80)
          j3 = min(j3,24)
          write (noutpt,1020) ustr80(1:j3),ux24(1:j2)
          write (nttyo,1020) ustr80(1:j3),ux24(1:j2)
          stop
        endif
c
c       Read the B-dot values on the log K temperature grid.
c       Return the data in the xdbval holding array.
c
        q500nd = q500fl
        call rdgrid(ndat0s,ndbmax,ndbptg,ndbptl,qend,qerr,q500nd,
     $  udbval,xdbval)
        if (qend) go to 990
        if (qerr) go to 995
c
c       Load the data into the bdot array.
c
c       Calling sequence substitutions:
c         bdot for zdbval
c
        call ldbar2(nacdpr,narxmx,narxt,ndbmax,ntprmx,ntprt,
     $  xdbval,bdot)
c
        if (ipch .ge. 0) then
c
c         Read the grid for the B-dot(H) parameter.
c
          ux24 = 'B-dot_H                 '
          ux24lc = ux24
          call locase(ux24lc)
          read (ndat0s,1000,end=990,err=995) ustr80
          call lejust(ustr80)
          call locase(ustr80)
c
          if (ustr80(1:24) .ne. ux24lc(1:24)) then
            j2 = ilnobl(ux24)
            j3 = ilnobl(ustr80)
            j3 = min(j3,24)
            write (noutpt,1020) ustr80(1:j3),ux24(1:j2)
            write (nttyo,1020) ustr80(1:j3),ux24(1:j2)
            stop
          endif
c
c         Read the B-dot(H) values on the log K temperature grid
c         Return the data in the xdbval holding array.
c
          q500nd = q500fl
          call rdgrid(ndat0s,ndbmax,ndbptg,ndbptl,qend,qerr,q500nd,
     $    udbval,xdbval)
          if (qend) go to 990
          if (qerr) go to 995
c
c         Load the data into the bdoth array.
c
c         Calling sequence substitutions:
c           bdoth for zdbval
c
          call ldbar2(nacdpr,narxmx,narxt,ndbmax,ntprmx,ntprt,
     $    xdbval,bdoth)
c
          do ipc = 1,ipch
c
            read (ndat0s,1000,end=990,err=995) ustr80
c
c           Read the B-dot(H) pressure derivatives on the
c           log K temperature grid. Return the data in the xdbval
c           holding array.
c
            q500nd = q500fl
            call rdgrid(ndat0s,ndbmax,ndbptg,ndbptl,qend,qerr,q500nd,
     $      udbval,xdbval)
            if (qend) go to 990
            if (qerr) go to 995
c
c           Load the data into the dbdth array.
c
c           Calling sequence substitutions:
c             dbdth for zdbval
c             ipchmx for ipcmax
c
            call ldbar3(ipc,ipchmx,nacdpr,narxmx,narxt,ndbmax,
     $      ntprmx,ntprt,xdbval,dbdth)
          enddo
        endif
c
        if (ipcv .ge. 0) then
c
c         Read the grid for the B-dot(V) parameter.
c
          ux24 = 'B-dot_V                 '
          ux24lc = ux24
          call locase(ux24lc)
          read (ndat0s,1000,end=990,err=995) ustr80
          call lejust(ustr80)
          call locase(ustr80)
c
          if (ustr80(1:24) .ne. ux24lc(1:24)) then
            j2 = ilnobl(ux24)
            j3 = ilnobl(ustr80)
            j3 = min(j3,24)
            write (noutpt,1020) ustr80(1:j3),ux24(1:j2)
            write (nttyo,1020) ustr80(1:j3),ux24(1:j2)
            stop
          endif
c
c         Read the B-dot(V) values on the log K temperature
c         grid. Return the data in the xdbval holding array.
c
          q500nd = q500fl
          call rdgrid(ndat0s,ndbmax,ndbptg,ndbptl,qend,qerr,q500nd,
     $    udbval,xdbval)
          if (qend) go to 990
          if (qerr) go to 995
c
c         Load the data into the bdotv array.
c
c         Calling sequence substitutions:
c           bdotv for zdbval
c
          call ldbar2(nacdpr,narxmx,narxt,ndbmax,ntprmx,ntprt,
     $    xdbval,bdotv)
c
          do ipc = 1,ipcv
c
            read (ndat0s,1000,end=990,err=995) ustr80
c
c           Read the B-dot(V) pressure derivatives on the
c           log K temperature grid. Return the data in the xdbval
c           holding array.
c
            q500nd = q500fl
            call rdgrid(ndat0s,ndbmax,ndbptg,ndbptl,qend,qerr,q500nd,
     $      udbval,xdbval)
            if (qend) go to 990
            if (qerr) go to 995
c
c           Load the data into the dbdtv array.
c
c           Calling sequence substitutions:
c             dbdtv for zdbval
c             ipcvmx for ipcmax
c
            call ldbar3(ipc,ipcvmx,nacdpr,narxmx,narxt,ndbmax,
     $      ntprmx,ntprt,xdbval,dbdtv)
          enddo
        endif
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c       Read the coefficients for computing the activity coefficient
c       of CO2(aq) from the Drummond (1981) equation.
c
        ux24 = 'Cco2   (coefficients for'
        ux24lc = ux24
        call locase(ux24lc)
        read (ndat0s,1000,end=990,err=995) ustr80
        call lejust(ustr80)
        call locase(ustr80)
c
        if (ustr80(1:24) .ne. ux24lc(1:24)) then
          ux24 = 'Cco2                    '
          ux24lc = ux24
          call locase(ux24lc)
        endif
c
        if (ustr80(1:24) .ne. ux24lc(1:24)) then
          j2 = ilnobl(ux24)
          j3 = ilnobl(ustr80)
          j3 = min(j3,24)
          write (noutpt,1020) ustr80(1:j3),ux24(1:j2)
          write (nttyo,1020) ustr80(1:j3),ux24(1:j2)
          stop
        endif
c
        read (ndat0s,1060,end=990,err=995) (cco2(n), n = 1,5)
 1060   format(5x,f10.4,11x,f12.7,/10x,f5.1,11x,f12.4,/5x,f10.6)
      endif
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
      if (uakey(1:8) .eq .'Pitzer  ') then
c
c       Read the grid for the Debye-Huckel A(phi) parameter and
c       the grids for related parameters.
c
        ux24 = 'Debye Huckel Aphi       '
        ux24lc = ux24
        call locase(ux24lc)
        read (ndat0s,1000,end=990,err=995) ustr80
        call lejust(ustr80)
        call locase(ustr80)
        if (ustr80(1:24) .ne. ux24lc(1:24)) then
          ux24 = 'Debye_Huckel A_phi      '
          ux24lc = ux24
          call locase(ux24lc)
        endif
c
        if (ustr80(1:24) .ne. ux24lc(1:24)) then
          j2 = ilnobl(ux24)
          j3 = ilnobl(ustr80)
          j3 = min(j3,24)
          write (noutpt,1020) ustr80(1:j3),ux24(1:j2)
          write (nttyo,1020) ustr80(1:j3),ux24(1:j2)
          stop
        endif
c
c       Read the Debye-Huckel A(phi) values on the log K temperature
c       grid. Return the data in the xdbval holding array.
c
        q500nd = q500fl
        call rdgrid(ndat0s,ndbmax,ndbptg,ndbptl,qend,qerr,q500nd,
     $  udbval,xdbval)
        if (qend) go to 990
        if (qerr) go to 995
c
c       Load the data into the aphi array.
c
c       Calling sequence substitutions:
c         aphi for zdbval
c
        call ldbar2(nacdpr,narxmx,narxt,ndbmax,ntprmx,ntprt,
     $  xdbval,aphi)
c
        if (ipch .ge. 0) then
c
c         Read the grid for the Debye-Huckel A(H) parameter.
c
          ux24 = 'Debye-Huckel A_H (kcal*k'
          ux24lc = ux24
          call locase(ux24lc)
          read (ndat0s,1000,end=990,err=995) ustr80
          call lejust(ustr80)
          call locase(ustr80)
          if (ustr80(1:24) .ne. ux24lc(1:24)) then
            ux24 = 'Debye-Huckel A_H        '
            ux24lc = ux24
            call locase(ux24lc)
          endif
c
          if (ustr80(1:24) .ne. ux24lc(1:24)) then
            j2 = ilnobl(ux24)
            j3 = ilnobl(ustr80)
            j3 = min(j3,24)
            write (noutpt,1020) ustr80(1:j3),ux24(1:j2)
            write (nttyo,1020) ustr80(1:j3),ux24(1:j2)
            stop
          endif
c
c         Read the Debye-Huckel A(H) values on the log K temperature
c         grid. Return the data in the xdbval holding array.
c
          q500nd = q500fl
          call rdgrid(ndat0s,ndbmax,ndbptg,ndbptl,qend,qerr,q500nd,
     $    udbval,xdbval)
          if (qend) go to 990
          if (qerr) go to 995
c
c         Load the data into the adhh array.
c
c         Calling sequence substitutions:
c           adhh for zdbval
c
          call ldbar2(nacdpr,narxmx,narxt,ndbmax,ntprmx,ntprt,
     $    xdbval,adhh)
c
          do ipc = 1,ipch
c
            read (ndat0s,1000,end=990,err=995) ustr80
c
c           Read the Debye-Huckel A(H) pressure derivatives on the
c           log K temperature grid. Return the data in the xdbval
c           holding array.
c
            q500nd = q500fl
            call rdgrid(ndat0s,ndbmax,ndbptg,ndbptl,qend,qerr,q500nd,
     $      udbval,xdbval)
            if (qend) go to 990
            if (qerr) go to 995
c
c           Load the data into the dadhh array.
c
c           Calling sequence substitutions:
c             dadhh for zdbval
c             ipchmx for ipcmax
c
            call ldbar3(ipc,ipchmx,nacdpr,narxmx,narxt,ndbmax,
     $      ntprmx,ntprt,xdbval,dadhh)
          enddo
        endif
c
        if (ipcv .ge. 0) then
c
c         Read the grid for the Debye-Huckel A(V) parameter.
c
          ux24 = 'Debye-Huckel A_V (cm**3*'
          ux24lc = ux24
          call locase(ux24lc)
          read (ndat0s,1000,end=990,err=995) ustr80
          call lejust(ustr80)
          call locase(ustr80)
          if (ustr80(1:24) .ne. ux24lc(1:24)) then
            ux24 = 'Debye-Huckel A_V        '
            ux24lc = ux24
            call locase(ux24lc)
          endif
c
          if (ustr80(1:24) .ne. ux24lc(1:24)) then
            j2 = ilnobl(ux24)
            j3 = ilnobl(ustr80)
            j3 = min(j3,24)
            write (noutpt,1020) ustr80(1:j3),ux24(1:j2)
            write (nttyo,1020) ustr80(1:j3),ux24(1:j2)
            stop
          endif
c
c         Read the Debye-Huckel A(V) values on the log K temperature
c         grid. Return the data in the xdbval holding array.
c
          q500nd = q500fl
          call rdgrid(ndat0s,ndbmax,ndbptg,ndbptl,qend,qerr,q500nd,
     $    udbval,xdbval)
          if (qend) go to 990
          if (qerr) go to 995
c
c         Load the data into the adhv array.
c
c         Calling sequence substitutions:
c           adhv for zdbval
c
          call ldbar2(nacdpr,narxmx,narxt,ndbmax,ntprmx,ntprt,
     $    xdbval,adhv)
c
          do ipc = 1,ipcv
c
            read (ndat0s,1000,end=990,err=995) ustr80
c
c           Read the Debye-Huckel A(V) pressure derivatives on the
c           log K temperature grid. Return the data in the xdbval
c           holding array.
c
            q500nd = q500fl
            call rdgrid(ndat0s,ndbmax,ndbptg,ndbptl,qend,qerr,q500nd,
     $      udbval,xdbval)
            if (qend) go to 990
            if (qerr) go to 995
c
c           Load the data into the dadhv array.
c
c           Calling sequence substitutions:
c             dadhv for zdbval
c             ipcvmx for ipcmax
c
            call ldbar3(ipc,ipcvmx,nacdpr,narxmx,narxt,ndbmax,
     $      ntprmx,ntprt,xdbval,dadhv)
          enddo
        endif
      endif
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Read the grid for the log K of the "Eh" reaction and the grids
c     for related parameters.
c
      ux24 = 'log k for eh reaction   '
      ux24lc = ux24
      call locase(ux24lc)
      read (ndat0s,1000,end=990,err=995) ustr80
      call lejust(ustr80)
      call locase(ustr80)
      if (ustr80(1:24) .ne. ux24lc(1:24)) then
        ux24 = 'Eh reaction: logKr (2H2O'
        ux24lc = ux24
        call locase(ux24lc)
      endif
      if (ustr80(1:24) .ne. ux24lc(1:24)) then
        ux24 = 'Eh reaction: logKr      '
        ux24lc = ux24
        call locase(ux24lc)
      endif
c
      if (ustr80(1:24) .ne. ux24lc(1:24)) then
        j2 = ilnobl(ux24)
        j3 = ilnobl(ustr80)
        j3 = min(j3,24)
        write (noutpt,1020) ustr80(1:j3),ux24(1:j2)
        write (nttyo,1020) ustr80(1:j3),ux24(1:j2)
        stop
      endif
c
c     Read the "Eh" reaction log K values on the log K temperature
c     grid. Return the data in the xdbval holding array.
c
      q500nd = q500fl
      call rdgrid(ndat0s,ndbmax,ndbptg,ndbptl,qend,qerr,q500nd,
     $ udbval,xdbval)
      if (qend) go to 990
      if (qerr) go to 995
c
c     Load the data into the xlke array.
c
c     Calling sequence substitutions:
c       xlke for zdbval
c
      call ldbar2(nacdpr,narxmx,narxt,ndbmax,ntprmx,ntprt,
     $ xdbval,xlke)
c
      if (ipch .ge. 0) then
c
c       Read the grid for the enthalpy of reaction of the "Eh"
c       reaction.
c
        ux24 = 'Eh reaction: delh0r (kca'
        ux24lc = ux24
        call locase(ux24lc)
        read (ndat0s,1000,end=990,err=995) ustr80
        call lejust(ustr80)
        call locase(ustr80)
        if (ustr80(1:24) .ne. ux24lc(1:24)) then
          ux24 = 'Eh reaction: delH0r     '
          ux24lc = ux24
          call locase(ux24lc)
        endif
c
        if (ustr80(1:24) .ne. ux24lc(1:24)) then
          j2 = ilnobl(ux24)
          j3 = ilnobl(ustr80)
          j3 = min(j3,24)
          write (noutpt,1020) ustr80(1:j3),ux24(1:j2)
          write (nttyo,1020) ustr80(1:j3),ux24(1:j2)
          stop
        endif
c
c       Read the "Eh" reaction enthalpy function values on the
c       log K temperature grid. Return the data in the xdbval
c       holding array.
c
        q500nd = q500fl
        call rdgrid(ndat0s,ndbmax,ndbptg,ndbptl,qend,qerr,q500nd,
     $  udbval,xdbval)
        if (qend) go to 990
        if (qerr) go to 995
c
c       Load the data into the xhfe array.
c
c       Calling sequence substitutions:
c         xhfe for zdbval
c
        call ldbar2(nacdpr,narxmx,narxt,ndbmax,ntprmx,ntprt,
     $  xdbval,xhfe)
c
        do ipc = 1,ipch
c
          read (ndat0s,1000,end=990,err=995) ustr80
c
c         Read the "Eh" reaction enthalpy function derivatives on
c         the log K temperature grid. Return the data in the xdbval
c         holding array.
c
          q500nd = q500fl
          call rdgrid(ndat0s,ndbmax,ndbptg,ndbptl,qend,qerr,q500nd,
     $    udbval,xdbval)
          if (qend) go to 990
          if (qerr) go to 995
c
c         Load the data into the dhfe array.
c
c         Calling sequence substitutions:
c           dhfe for zdbval
c           ipchmx for ipcmax
c
          call ldbar3(ipc,ipchmx,nacdpr,narxmx,narxt,ndbmax,
     $    ntprmx,ntprt,xdbval,dhfe)
        enddo
      endif
c
      if (ipcv .ge. 0) then
c
c       Read the grid for the volume of reaction of the "Eh"
c       reaction.
c
        ux24 = 'Eh reaction: delv0r (cm*'
        ux24lc = ux24
        call locase(ux24lc)
        read (ndat0s,1000,end=990,err=995) ustr80
        call lejust(ustr80)
        call locase(ustr80)
        if (ustr80(1:24) .ne. ux24lc(1:24)) then
          ux24 = 'Eh reaction: delV0r     '
          ux24lc = ux24
          call locase(ux24lc)
        endif
c
        if (ustr80(1:24) .ne. ux24lc(1:24)) then
          j2 = ilnobl(ux24)
          j3 = ilnobl(ustr80)
          j3 = min(j3,24)
          write (noutpt,1020) ustr80(1:j3),ux24(1:j2)
          write (nttyo,1020) ustr80(1:j3),ux24(1:j2)
          stop
        endif
c
c       Read the "Eh" reaction volume function values on the
c       log K temperature grid. Return the data in the xdbval
c       holding array.
c
        q500nd = q500fl
        call rdgrid(ndat0s,ndbmax,ndbptg,ndbptl,qend,qerr,q500nd,
     $  udbval,xdbval)
        if (qend) go to 990
        if (qerr) go to 995
c
c       Load the data into the xvfe array.
c
c       Calling sequence substitutions:
c         xvfe for zdbval
c
        call ldbar2(nacdpr,narxmx,narxt,ndbmax,ntprmx,ntprt,
     $  xdbval,xvfe)
c
        do ipc = 1,ipcv
c
          read (ndat0s,1000,end=990,err=995) ustr80
c
c         Read the "Eh" reaction volume function derivatives on
c         the log K temperature grid. Return the data in the xdbval
c         holding array.
c
          q500nd = q500fl
          call rdgrid(ndat0s,ndbmax,ndbptg,ndbptl,qend,qerr,q500nd,
     $    udbval,xdbval)
          if (qend) go to 990
          if (qerr) go to 995
c
c         Load the data into the dvfe array.
c
c         Calling sequence substitutions:
c           dvfe for zdbval
c           ipcvmx for ipcmax
c
          call ldbar3(ipc,ipcvmx,nacdpr,narxmx,narxt,ndbmax,ntprmx,
     $    ntprt,xdbval,dvfe)
        enddo
      endif
c
      read (ndat0s,1000,end=990,err=995) ustr80
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
      go to 999
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
  990 write (noutpt,2000)
      write (nttyo,2000)
 2000 format(/' * Error - (EQPT/rdpar) Unexpectedly encountered',
     $ /7x,'end-of-file while reading the DATA0 file.')
      stop
c
  995 write (noutpt,2010)
      write (nttyo,2010)
 2010 format(/' * Error - (EQPT/rdpar) Encountered a read format',
     $ /7x,'error while reading the DATA0 file.')
      stop
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
  999 continue
      end
