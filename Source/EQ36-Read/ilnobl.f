      integer function ilnobl(ustr)
c
c     This subroutine finds the position of the last non-blank character
c     in the string ustr.
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
c       ilnobl = the position of the first non-blank character
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
      integer j,nchars
c
c-----------------------------------------------------------------------
c
c     Get the length of the string variable.
c
      nchars = len(ustr)
c
c     Find the first non-blank character.
c
      ilnobl = 0
      do j = nchars,1,-1
        if (ustr(j:j) .ne. ' ') then
          ilnobl = j
          go to 999
        endif
      enddo
c
  999 continue
      end
