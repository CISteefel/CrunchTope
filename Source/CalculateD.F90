      subroutine dhconst(dhad,dhbd,tempk,tconv)

      implicit real*8 (a-h,o-z)
 
!!  compute temperature in [deg C]
 
      tempc = tempk - tconv
    
!!  this subroutine which computes the Debye-Huckel
!!  constants as a function of temperature was taken directly
!!  from WATEQ2 (Ball et al. 1979). The mathematics is documented
!!  in Truesdell and Jones (1974).
 
      s1 = 374.11d0-tempc
      s2 = s1**(1.0d0/3.0d0)
      bob1 = (1.0d0+0.1342489d0*s2-3.946263d-03*s1)
      bob2 = (3.1975d0-0.3151548d0*s2-1.203374d-3*s1)
      bob3 = (7.48908d-13*s1**4.0d0)
      bob4 = bob2+bob3
      s3 = dsqrt(bob1/bob4)
      if (tempk.lt.373.16d0) then
         cc1 = 87.74d0-tempc*(tempc*(1.41d-6*tempc-9.398d-4)+0.4008d0)
      else
         cc1 = 5321d0/tempk+233.76d0-tempk*(tempk*(8.29d-7*tempk - 1.417d-3)+.9297d0)
      endif
      cc1 = dsqrt(cc1*tempk)
      dhad = 18246d2*s3/cc1**3.0d0
      dhbd = 50.29d0*s3/cc1
      return
      end subroutine dhconst