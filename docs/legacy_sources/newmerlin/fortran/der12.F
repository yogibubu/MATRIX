!############################################################
!#                  SUBROUTINE: GRAD_X_TO_S                 #
!############################################################
subroutine grad_x_to_s(gsvec, gxvec, bpmat, pmat, nat3, nvib)
    
    !subroutine grad_x_to_s
    !
    !DESCRIPTION
    !  Converts the Cartesian gradient to the internal
    !  coordinates space.
    !
    !ARGUMENTS
    !  gsvec  >>>>  gradient in non-redundant internal
    !               coordinates
    !  gxvec  >>>>  gradient in Cartesian coordinates
    !  bpmat  >>>>  pseudoinverse of the non-redundant Wilson 
    !               B matrix
    !  pmat   >>>>  projection matrix (B trasposta per B)
    !  nat3   >>>>  number of Cartesian coordinates
    !  nvib   >>>>  number of normal modes
    !
    
    implicit none 
    integer :: nat3, nvib
    real(kind=dp) :: gxvec(nat3), pgvec(nat3), gsvec(nvib), &
        & pmat(nat3,nat3), bpmat(nat3,nvib)
    
    !Project out rotational and translational motions from
    !the Cartesian gradient of the moments of inertia:
    pgvec = matmul(pmat, gxvec)
    
    !Convert the projected gradient:
    gsvec = matmul(transpose(bpmat), pgvec)
    
    return
end subroutine grad_x_to_s
    
    
!############################################################
!#                  SUBROUTINE: HESS_X_TO_S                 #
!############################################################
subroutine hess_x_to_s(hsmat, hxmat, bpmat, gsvec, dbnmat, pmat, nat3, nvib, hxproj) 
    
    !SUBROUTINE hess_x_to_s
    !
    !DESCRIPTION
    !  Converts the Cartesian Hessian matrix to the internal
    !  coordinates space.
    !
    !ARGUMENTS
    !  hsmat  >>>>  hessian matrix in non-redundant internal
    !               coordinates
    !  hxmat  >>>>  Cartesian hessian matrix
    !  bpmat  >>>>  pseudoinverse of the non-redundant Wilson 
    !               B matrix
    !  gsvec  >>>>  gradient in non-redundant internal
    !               coordinates
    !  dbnmat >>>>  first-order Cartesian derivatives of
    !               the non-redundant Wilson B matrix
    !  pmat   >>>>  projection matrix 
    !  nat3   >>>>  number of Cartesian coordinates
    !  nvib   >>>>  number of normal modes
    !
    
    implicit none 
    integer :: nat3, nvib, i, j, k 
    real(kind=dp) :: gsvec(nvib), pmat(nat3,nat3), bpmat(nat3,nvib), gtmat(nat3,nat3), &
        & dbnmat(nvib,nat3,nat3), hxmat(nat3,nat3), hpmat(nat3,nat3), hsmat(nvib,nvib), &
        & htmat(nat3,nat3)
    real(kind=dp), optional :: hxproj(nat3,nat3)

    !Compute gradient term (correction to the Cartesian 
    !Hessian matrix):
    gtmat = zero
    do i = 1, nat3
        do j = 1, nat3
            do k = 1, nvib
                gtmat(i,j) = gtmat(i,j) + gsvec(k)*dbnmat(k,i,j)
            end do
        end do
    end do
    
    !Correct Hessian matrix:
    htmat = hxmat - gtmat 
    
    !Project out rotational and translational motions from
    !the Cartesian Hessian matrix:
    hpmat = matmul(pmat, matmul(htmat, pmat))
    
    if (present(hxproj)) hxproj = hpmat
    
    !Convert the Cartesian Hessian matrix to the internal 
    !coordinates space:
    hsmat = matmul(matmul(transpose(bpmat), hpmat), bpmat)
    
    return
end subroutine hess_x_to_s

