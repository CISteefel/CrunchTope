      subroutine locase(ustr)
c
c     This subroutine converts a string from upper case to lower case.
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
      integer idel,j,nchars
c
      character*1 u1
c
c-----------------------------------------------------------------------
c
c
      idel = ichar('A') - ichar('a')
      if (idel .ne. 0) then
        nchars = len(ustr)
        do j = 1,nchars
          u1 = ustr(j:j)
          if (u1.ge.'A' .and. u1.le.'Z')
     $    ustr(j:j) = char(ichar(u1) - idel)
        enddo
      endif
c
      end
