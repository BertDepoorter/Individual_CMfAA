module mod_usr
  use mod_hd

  implicit none

contains

  subroutine usr_init()

    usr_init_one_grid => rm1d_init_one_grid

    call set_coordinate_system("Cartesian")
    call hd_activate()

  end subroutine usr_init

  ! Initialize one grid
  subroutine rm1d_init_one_grid(ixGmin1,ixGmax1,ixmin1,ixmax1,w,x)
    integer, intent(in) :: ixGmin1,ixGmax1, ixmin1,ixmax1
    double precision, intent(in) :: x(ixGmin1:ixGmax1,1:ndim)
    double precision, intent(inout) :: w(ixGmin1:ixGmax1,1:nw)

    ! iprob==1 rarefaction wave & shock
    ! test one of lecture notes
    if (iprob==1) then
        where (abs(x(ixmin1:ixmax1,1))<0.5d0)
           w(ixmin1:ixmax1,rho_)   = 1d0
           w(ixmin1:ixmax1,mom(1)) = 0.0d0
          
        elsewhere
           w(ixmin1:ixmax1,rho_)   = 0.125d0
           w(ixmin1:ixmax1,mom(1)) = 0.0d0
          
        end where

    ! iprob==2  shock & shock
    ! test 2 of the lecture notes
    else if (iprob==2) then
        where (abs(x(ixmin1:ixmax1,1))<0.5d0)
           w(ixmin1:ixmax1,rho_)   = 0.445d0
           w(ixmin1:ixmax1,mom(1)) = 0.31061d0
         
        elsewhere
           w(ixmin1:ixmax1,rho_)   = 0.5d0
           w(ixmin1:ixmax1,mom(1)) = 0d0
          
        end where
    ! iprob==3  rarefaction wave
    ! test 3 of the lecture notes
    else if (iprob==3) then
        where (abs(x(ixmin1:ixmax1,1))<0.5d0)
           w(ixmin1:ixmax1,rho_)   = 0.5d0
           w(ixmin1:ixmax1,mom(1)) = 0d0
         
        elsewhere
           w(ixmin1:ixmax1,rho_)   = 0.445d0
           w(ixmin1:ixmax1,mom(1)) = 0.31061d0
       
        end where
    else
        call mpistop("iprob not available!")
    end if

    call hd_to_conserved(ixGmin1,ixGmax1,ixmin1,ixmax1,w,x)

  end subroutine rm1d_init_one_grid

end module mod_usr
