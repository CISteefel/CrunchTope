        NOTE:  Inner loop is read first, so JX first (could change this to cycle over JY, but this is the way it is now)
		
		DO jy = 0,ny
          DO jx= 1,nx
            READ(23,*,END=1020) qy(jx,jy,jz)
          END DO
        END DO

Velocity in Y coordinate direction
1,0
2,0
3,0
4,0

1,1
2,1
3,1
4,1

1,2
2,2
3,2
4,2

1,3
2,3
3,3
4,3

1,4
2,4
3,4
4,4


1.0
2.0
3.0
4.0

1.0
2.0
3.0
4.0

1.1
2.1
3.1
4.1

1.2
2.2
3.2
4.2

1.3
2.3
3.3
4.3

1.4
2.4
3.4
4.4
