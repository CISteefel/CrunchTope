      subroutine rdgrid(ndat0s,ndbmax,ndbptg,ndbptl,qend,qerr,q500nd,
     $ udbval,xdbval)
c
c     This subroutine reads data from the standard log K temperature
c     grid. Other kinds of data are read from such grids.
c
c     This suboutine is called by:
c
c       EQPT/rdpar.f
c       EQPT/pcraq.f
c       EQPT/pcrsg.f
c
c-----------------------------------------------------------------------
c
c     Principal input:
c
c       ndbptg = the number of points on the grid
c       ndbptl = the maximum number of points on a single line
c
c     Principal output:
c
c       qend   = logical flag, = .true. if end-of-file was encountered
c       qerr   = logical flag, = .true. if a read format error was
c                  encountered
c       udbval = string array, the equivalent of xdbval
c       xdbval = array containing the data read from the grid
c
c-----------------------------------------------------------------------
c
      implicit none
c
c-----------------------------------------------------------------------
c
c     Calling sequence variable declarations.
c
      integer ndbmax
c
      integer ndat0s
c
      integer ndbptg,ndbptl
c
      logical qend,qerr,q500nd
c
      character(len=16) udbval(ndbmax)
c
      real(8) xdbval(ndbmax)
c
c-----------------------------------------------------------------------
c
c     Local variable declarations.
c
      integer i,ii,j,jj,k,n
c
      character(len=80) uline,ulbufa,ulbufb
      character(len=16) ux16
c
c-----------------------------------------------------------------------
c
      qend = .false.
      qerr = .false.
c
c     Note: the following coding generally assumes that a blank
c     space trails all but the last number on a line. If a blank
c     space is not found, it is assumed that a number ends after
c     four decimal places.
c
      n = 0
      do i = 1,ndbptg/ndbptl
        read (ndat0s,'(a)',end=990,err=995) uline
        ulbufa = uline
        do k = 1,ndbptl
          call lejust(ulbufa)
          ii = index(ulbufa,' ')
          jj = index(ulbufa,'.')
          if (jj .eq. 0) then
            jj = 80
          else
            jj = jj + 5
          endif
          ii = min(ii,jj)
          if (ii .gt. 1) then
            ux16 = ulbufa(1:ii - 1)
            call locase(ux16)
            udbval(n + k) = ux16
            if (ux16(1:7) .eq. 'no_data') then
              xdbval(n + k) = 9999999.
            else
              read (ux16,'(f10.4)',err=995) xdbval(n + k)
            endif
          else
            udbval(n + k) = ' '
            xdbval(n + k) = 0.
          endif
          if (k .lt. ndbptl) then
            ulbufb = ulbufa(ii:80)
            ulbufa = ulbufb
          endif
        enddo
c       read (uline,udbfmt,err=995) (xdbval(n + k), k = 1,ndbptl)
c       Note: udbfmt = '( (5x,6f10.4) )' with the "6" repeat
c       count replaced by ndbptl.
        n = n + ndbptl
      enddo
c
      j = mod(ndbptg,ndbptl)
      if (j .gt. 0) then
        read (ndat0s,'(a)',end=990,err=995) uline
        ulbufa = uline
        do k = 1,j
          call lejust(ulbufa)
          ii = index(ulbufa,' ')
          jj = index(ulbufa,'.')
          if (jj .eq. 0) then
            jj = 80
          else
            jj = jj + 5
          endif
          ii = min(ii,jj)
          if (ii .gt. 1) then
            ux16 = ulbufa(1:ii - 1)
            call locase(ux16)
            udbval(n + k) = ux16
            if (ux16(1:7) .eq. 'no_data') then
              xdbval(n + k) = 9999999.
            else
              read (ux16,'(f10.4)',err=995) xdbval(n + k)
            endif
          else
            udbval(n + k) = ' '
            xdbval(n + k) = 0.
          endif
          if (k .lt. j) then
            ulbufb = ulbufa(ii:80)
            ulbufa = ulbufb
          endif
        enddo
c       read (uline,udbfmt,err=995) (xdbval(n + k), k = 1,j)
        n = n + j
      endif
c
      if (q500nd) then
c
c       Apply the following filter. Treat values of "500.0000" as
c       "no data". This accommodates older EQ3/6 data files which
c       use "500.0000" to mean no data.
c
        do i = 1,ndbptg
          if (xdbval(i).gt.499.9999 .and. xdbval(i).lt.500.0001) then
            udbval(i) = 'no_data'
            xdbval(i) = 9999999.
          endif
        enddo
      endif
c
      go to 999
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
  990 qend = .true.
      go to 999
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
  995 qerr = .true.
      go to 999
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
  999 continue
      end
