program hbond_main
  implicit none
  integer :: ic, nconf, npart, nwat, iw, jw
  integer :: ib, ix, iy, iz, ncheck, it, ip
  integer :: ih, jh
  integer :: nc, nhc, list(10), list_ix(10), list_iy(10), list_iz(10)
  integer :: nhbond(10), nhbond5(10)
  
  real*8  :: x(125), y(125), z(125), ri, pi
  real*8  :: hx1, hy1, hz1, hx2, hy2, hz2, roh
  real*8  :: hx3, hy3, hz3, hx4, hy4, hz4
  real*8  :: ohx1(125), ohy1(125), ohz1(125)
  real*8  :: ohx2(125), ohy2(125), ohz2(125)
  real*8  :: xij, yij, zij, r, box
  real*8  :: noh1(50), noh2(50), noh3(50), noh4(50), ncoord(50)
  real*8  :: xi, yi, zi, xj, yj, zj, theta, sum, ave
  character(len=2) :: atnm

!  print *, 'NCONF '
!  read *, nconf

  nconf  = 1250

  ncoord = 0.0d0
  noh1 = 0.0d0
  noh2 = 0.0d0
  noh3 = 0.0d0
  noh4 = 0.0d0
  pi   = 3.14159265358979323846264338d0
  box   = 9.858d0
  open (15, file='trajectory.xyz', status='old')
  open (25, file='hbond_info.dat', status='unknown')
  open (40, file='snapshot.xyz',   status='unknown')

  nhbond = 0
  nhbond5 = 0
  
  do ic = 1, nconf
     read (15, *) npart
     nwat = npart/3
     ncheck = 0
     read (15, *)
    
     do iw = 1, nwat
        read (15, *) atnm, x(iw), y(iw), z(iw)
        read (15, *) atnm, hx1, hy1, hz1
        read (15, *) atnm, hx2, hy2, hz2

        ohx1(iw) = hx1 - x(iw)
        ohy1(iw) = hy1 - y(iw)
        ohz1(iw) = hz1 - z(iw)
        ohx2(iw) = hx2 - x(iw)
        ohy2(iw) = hy2 - y(iw)
        ohz2(iw) = hz2 - z(iw)

        roh   = sqrt(ohx1(iw)**2 + ohy1(iw)**2 + ohz1(iw)**2)
        ohx1(iw) = ohx1(iw)/roh
        ohy1(iw) = ohy1(iw)/roh
        ohz1(iw) = ohz1(iw)/roh
        
        roh   = sqrt(ohx2(iw)**2 + ohy2(iw)**2 + ohz2(iw)**2)
        ohx2(iw) = ohx2(iw)/roh
        ohy2(iw) = ohy2(iw)/roh
        ohz2(iw) = ohz2(iw)/roh
     end do

     ! -- O -- O --
     do iw = 1, nwat
        xi = x(iw)
        yi = y(iw)
        zi = z(iw)
        
        hx1 = ohx1(iw)
        hy1 = ohy1(iw)
        hz1 = ohz1(iw)
        
        hx2 = ohx2(iw)
        hy2 = ohy2(iw)
        hz2 = ohz2(iw)

        nc  = 0
        nhc = 0
        list = 0
        
        do jw = 1, nwat
           if (jw .eq. iw) cycle
           xj = x(jw)
           yj = y(jw)
           zj = z(jw)

           hx3 = ohx1(jw)
           hy3 = ohy1(jw)
           hz3 = ohz1(jw)
           
           hx4 = ohx2(jw)
           hy4 = ohy2(jw)
           hz4 = ohz2(jw)
        
           do ix = -1, 1
              do iy = -1, 1
                 do iz = -1, 1
                    xij = (xj + ix*box) - xi
                    yij = (yj + iy*box) - yi
                    zij = (zj + iz*box) - zi
        
                    r  = sqrt(xij*xij + yij*yij + zij*zij)
                    
                    if (r .lt. 3.36d0) then
                       xij = xij/r
                       yij = yij/r
                       zij = zij/r
                       nc = nc + 1
                       list(nc) = jw
                       list_ix(nc) = ix
                       list_iy(nc) = iy
                       list_iz(nc) = iz
                       
                       ncoord(iw) = ncoord(iw) + 1
                       ! iw - jw
                       ! calculate angle of H1-O-O
                       theta = acos(xij*hx1 + yij*hy1 + zij*hz1)*180.0d0/pi
                       if (theta .le. 40.0d0) then
                          nhc = nhc + 1
                          noh1(iw) = noh1(iw) + 1.0d0
                       end if
                       ! calculate angle of H2-O-O
                       theta = acos(xij*hx2 + yij*hy2 + zij*hz2)*180.0d0/pi
                       if (theta .le. 40.0d0) then
                          nhc = nhc + 1
                          noh2(iw) = noh2(iw) + 1.0d0
                       end if

                       ! jw-iw
                       ! calculate angle of H3-O-O
                       theta = acos(-xij*hx3 - yij*hy3 - zij*hz3)*180.0d0/pi
                       if (theta .le. 40.0d0) then
                          nhc = nhc + 1
                          noh3(iw) = noh3(iw) + 1.0d0
                       end if
                       ! calculate angle of H4-O-O
                       theta = acos(-xij*hx4 - yij*hy4 - zij*hz4)*180.0d0/pi
                       if (theta .le. 40.0d0) then
                          nhc = nhc + 1
                          noh3(iw) = noh3(iw) + 1.0d0
                       end if
                       
                    end if
! DEBUG
                 end do
              end do
           end do

           
           
        end do

        nhbond(nhc) = nhbond(nhc) + 1
        if (nc .eq. 5) then
           nhbond5(nhc) = nhbond5(nhc) + 1
        end if
        if (nc .ne. nhc) then
           write (40, *) (nc+1)*3
           write (40, *) nc, nhc
           !xi = x(iw)
           !yi = y(iw)
           !zi = z(iw)
           write (40, *) ' O ', xi, yi, zi
           write (40, *) ' H ', xi+hx1, yi+hy1, zi+hz1
           write (40, *) ' H ', xi+hx2, yi+hy2, zi+hz2
           do it = 1, nc
              jw = list(it)
              ix = list_ix(it)
              iy = list_iy(it)
              iz = list_iz(it)
              xi = x(jw) + ix*box 
              yi = y(jw) + iy*box
              zi = z(jw) + iz*box
              write (40, *) ' O ', xi, yi, zi
              write (40, *) ' H ', xi+ohx1(jw), yi+ohy1(jw), zi+ohz1(jw)
              write (40, *) ' H ', xi+ohx2(jw), yi+ohy2(jw), zi+ohz2(jw)
           end do
        end if
           
     end do
     !stop
     
  end do

  ncoord = ncoord /(nconf)
  noh1   = noh1   /(nconf)
  noh2   = noh2   /(nconf)
  noh3   = noh3   /(nconf)
  !noh4   = noh4   /(nconf)
  
  do iw = 1, nwat
     sum = noh1(iw) + noh2(iw) + noh3(iw) !+ noh4(iw)
     write (25, *) ncoord(iw), sum, noh1(iw), noh2(iw), noh3(iw)!, noh4(iw)
  end do

  sum = 0.0d0
  
  do nhc = 1, 10
     sum = sum + nhbond(nhc)
  end do

  ave = 0.0d0
  do nhc = 1, 10
!     write (35, *) nhc, dble(nhbond(nhc))/sum
     ave = ave + dble(nhc)*dble(nhbond(nhc))/sum
  end do

  print *, 'AVE', ave

 
  
!  sum = 0.0d0
!  do nhc = 1, 10
!     sum = sum + nhbond5(nhc)
!  end do

!  do nhc = 1, 10
!     write (36, *) nhc, dble(nhbond5(nhc))/sum
!  end do

end program hbond_main
