! Ensemble Kalman Square Root Filter 
! Copyright (C) 2022 J.J.D. Hooghiem 

! This program is free software: you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.

! This program is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
! GNU General Public License for more details.

! You should have received a copy of the GNU General Public License
! along with this program.  If not, see <https://www.gnu.org/licenses/>.

module enkf_core
  
  implicit none (type,external)

  ! Reference and definition of the BLAS double dot function
  real(kind=8),external :: ddot

  ! Reference and definition of the BLAS routines required 
  external :: dgemv,dger,daxpy

  private 

  public :: enksrf, dckcmm

contains

  subroutine enksrf(nobs,nmembers,nparams,obs,mrej,rej_thr,hx,hxp,xp,r, hphr,x,rejected,assimilate)
    !
    ! This subroutine implements the Ensemble Kalman Square Root Filter 
    ! algorithm as described in 
    !    "Ensemble Data Assimilation without Perturbed Observations"
    ! by Whitaker and Hamill 2002. DOI: 10.1175/1520-0493(2002)130<1913:EDAWPO>2.0.CO;2
    ! and  
    !    "An ensemble data assimilation system to estimate CO2 surface 
    !    "fluxes from atmospheric trace gas observations" 
    ! by Peters et al. 2005. DOI: 10.1029/2005JD006157
    ! 

    ! Input variables, not changed during routine
    integer*8,intent(in)                              :: nobs     ! the amount of observations
    integer*8,intent(in)                              :: nmembers ! the amount of ensemble members used to represent the covariance
                                                                  ! strucure (sometimes called particle number) 
    integer*8,intent(in)                              :: nparams  ! number of statevector elements to be estimated
    logical,dimension(nobs), intent(in)               :: mrej     ! may reject wether an observations may be rejected
    logical,dimension(nobs), intent(in)               :: assimilate ! wether we want to assimilate this observation 
    real*8, dimension(nobs), intent(in)               :: obs      ! observation values
    real*8, dimension(nobs), intent(in)               :: R        ! observation error  
    real*8, dimension(nobs), intent(in)               :: rej_thr  ! rejection threshold for observation that may be rejected

    ! local variables
    integer*8                                         :: i,j       !
    real*8                                            :: n_fac     ! for the repeated factor of  1/(N-1), n=nmembers
    real*8                                            :: res,alpha ! residual, alpha factor (see Whitaker and Hamill 2002)
    real*8,dimension(nparams)                         :: PHt       ! Estimated product of PH-transpose
    real*8,dimension(nparams)                         :: KG        ! Kalman Gain 
    real*8,dimension(nobs)                            :: fac       ! factor used for updates
    real*8,dimension(nmembers)                        :: HXp_n     ! temporary row of deviations  

    !
    ! Output (updated in the loop)
    ! 
    real*8, dimension(nobs,nmembers),intent(inout)    :: HXp     ! observed deviations
    real*8, dimension(nobs),intent(inout)             :: Hx,HPHR ! simulationed mean value of the observations, HPH+R
    real*8, dimension(nparams,nmembers),intent(inout) :: Xp      ! statevecter deviations
    real*8, dimension(nparams),intent(inout)          :: x       ! statevector
    integer*8, dimension(nobs),intent(inout)          :: rejected ! returns True if observation is rejected

    ! initialize variables
    PHt=0.0 
    KG=0.0 
    fac=0.0 
    n_fac=1.0/(float(nmembers)-1.0)
       
    ! start processing observations one at a time
    do i=1,nobs
        if (assimilate(i).eqv..False.) then
                cycle 
        endif
        ! compute difference between forecast and observed:
        res           = obs(i) - Hx(i)

        ! We should be checking if we reject observations here by comparing to res(idual)
        if (mrej(i).eqv..True.) then 
                if ( abs(res) > rej_thr(i) * sqrt(R(i) )) then 
                        rejected(i)=1
                        cycle
                endif
        endif 

        ! Local copy of the deviations row HXp(i,:) will be updated later
        HXp_n=HXp(i,:)

        ! calculate PHt using blas/lapac DGEMV 
        ! Computes PHt = 1/(N-1) HXp_n . Xp 
        call DGEMV('N',nparams,nmembers,n_fac,Xp,nparams,HXp_n,1,0.0,PHt,1)

        ! estimates HPH+R=1/(N-1) HXp_n . HXp_n + R  , using DDOT. 
        HPHR(i) = n_fac * DDOT(nmembers,HXp_n,1,HXp_n,1) + R(i)

        ! Compute the Kalman Gain
        KG = PHt/HPHR(i)

        ! alpha factor
        alpha=1.0 / (1.0 + sqrt( R(i) / HPHR(i) ))

        ! We now have the info to compute the updated version of
        ! x, Xp, HXp, and Hx  
        ! some of this stuf is not final yet 

        ! could it be more effiecient to do this and the part of fac in order to update Hx?
        ! call DGEMV('N',nobs,nmembers,res*n_fac*(1/HPHR(i)),HXp,nmembers,HXp(i,:),1,1.0,Hx,1)

        ! Approximate HK, and store in fac
        call DGEMV( 'N' , nobs , nmembers , n_fac * (1.0 / HPHR(i)),HXp,nobs,HXp_n,1,0.0,fac,1)

        ! update Hx = Hx + HK * residual 
        call daxpy(nobs,res,fac,1,Hx,1)

        ! update x = x  + res * KG 
        call daxpy(nparams,res,KG,1,x,1)

        ! update deviations HX = HX - HK HXp(i) using outer product
        ! call DGER(nobs,nmembers,-1.0*alpha,fac,1,HXp_n,1,HXp,nobs)

        ! update Xp 
        ! call DGER(nparams,nmembers,-1.0*alpha,KG,1,HXp_n,1,Xp,nparams)

        ! update deviations HX = HX - HK HXp(i) using outer product in loop
        do j=1,nmembers
                ! update Xp
                call daxpy(nparams,-1.0*alpha*HXp_n(j),KG,1,Xp(:,j),1)
                ! update HXp for column j
                call daxpy(nobs,-1.0*alpha*HXp_n(j) ,fac,1,HXp(:,j),1)
        enddo ! loop over nmembers 

    ! end loop over observations
    enddo 

  end subroutine  enksrf
  subroutine dckcmm(N,M,K,NM,A,B,C,R)
         !
         ! Computes the 
         !       matmul( kron(A,B), C ) 
         !       where A(N,N) and B(M,M) are cholesky decompositions
         !       C is a matrix with dimension (NM,K)
         !       and stores the result in R in R(NM,K)
         !         
         integer*8,intent(in) :: N, M, K, NM
         real*8,dimension(N,N),intent(in) :: A
         real*8,dimension(M,M),intent(in) :: B
         real*8,dimension(NM,K),intent(in) :: C
         real*8,dimension(NM,K),intent(inout) :: R 
         integer*8            :: i,j,h,l,g
         real*8 :: aa
         do i=1,N
          do j=1,i
           do l=1,M
            do g=1,K
             aa=A(i,j)*C(l+(j-1)*M,g)
             do h=l,M
              R(h+(i-1)*M,g) =R(h+(i-1)*M,g)+  aa* B(h,l) 
             enddo
            enddo
           enddo
          enddo
         enddo

  end subroutine dckcmm

end module enkf_core
