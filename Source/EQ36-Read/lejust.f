      subroutine lejust(ustr)
c
c     This subroutine left-justifies the non-blank portion of the string
c     ustr.
c
c     This subroutine is called by:
c
c       Any
c
c-----------------------------------------------------------------------
c
c     Input:
c
c       ustr   = the input string variable
c
c     Output:
c
c       ustr   = the output string variable
c
c-----------------------------------------------------------------------
c
      implicit none
c
c-----------------------------------------------------------------------
c
c     Calling sequence variable declarations.
c
      character*(*) ustr
c
c-----------------------------------------------------------------------
c
c     Local variable declarations.
c
      integer j,jj,jbl,j1,nchars
c
      integer ifnobl
c
c-----------------------------------------------------------------------
c
c     Get the length of the string variable.
c
      nchars = len(ustr)
c
c     Get the position of the first non-blank character and the number
c     of blanks on the left-hand-side.
c
      j1 = ifnobl(ustr)
      jbl = j1 - 1
c
      if (jbl .gt. 0) then
        do jj = j1,nchars
          j = jj - jbl
          ustr(j:j) = ustr(jj:jj)
        enddo
        do j = nchars - jbl + 1,nchars
          ustr(j:j) = ' '
        enddo
      endif
c
      end
