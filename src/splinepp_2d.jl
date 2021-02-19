"""

Convert 2d spline in B form to spline in pp form with periodic boundary 
conditions. This is a special case of the procedure spline_pp_b_to_pp_2d 
for the double periodic case to avoid the select case statements
- self : arbitrary degree 2d spline 
- n_cells : number of gridcells
- b_coeffs(n_cells(1)*n_cells(2)) : coefficients of spline in B-form
- pp_coeffs((self%spline1%degree+1)*(self%spline2%degree+1),n_cells(1)*n_cells(2)) : coefficients of spline in pp-form
"""
function spline_pp_b_to_pp_2d_periodic( self, n_cells, b_coeffs, pp_coeffs)

    degree1 = self[1].degree
    degree2 = self[2].degree

    for j=1:n_cells(2), i=1:n_cells(1)
        spline_pp_b_to_pp_2d_cell(self, n_cells, b_coeffs, pp_coeffs, i,j)
    end

end 

#=


  !> Convert 2d spline in B form to spline in pp form   
  subroutine sll_s_spline_pp_b_to_pp_2d( self, n_cells, b_coeffs, pp_coeffs)
    type( sll_t_spline_pp_2d), intent(inout)::  self !< arbitrary degree 2d spline 
    sll_int32,  intent(in) ::  n_cells(2) !< number of gridcells
    sll_real64, intent(in) :: b_coeffs(:)!(n_cells(1)*n_cells(2))   !< coefficients of spline in B-form
    sll_real64, intent(out):: pp_coeffs((self%spline1%degree+1)*(self%spline2%degree+1),n_cells(1)*n_cells(2))  !< coefficients of spline in pp-form
    sll_int32 :: i,j,l1,l2
    sll_int32 :: degree1,degree2
    sll_real64 :: pp_coeffs_local(self%spline1%degree+1,self%spline2%degree+1)
    sll_int32 :: offset1, offset2, cell, n_coeffs(2), index, upper(2)


    degree1= self%spline1%degree
    degree2= self%spline2%degree
    n_coeffs(1) = self%spline1%n_coeffs
    n_coeffs(2) = self%spline2%n_coeffs

    select case ( self%spline1%boundary_conditions )
    case (sll_p_boundary_periodic )
       offset1 = -degree1
    case (sll_p_boundary_clampeddiri:sll_p_boundary_clampeddiri_clamped)
       offset1 = -1
    case default
       offset1 = 0
    end select

    upper(1) = n_cells(1)-degree1+1
    upper(2) = n_cells(2)-degree2+1
    do j=1, n_cells(2)
       do i=1, n_cells(1)
          cell = (j-1)*n_cells(1)+i
          do l2=0,degree2
             select case ( self%spline2%boundary_conditions )
             case( sll_p_boundary_periodic )
                offset2 = modulo(j-degree2+l2-1,n_coeffs(2))*n_coeffs(1) ! periodic
             case( sll_p_boundary_clampeddiri )
                offset2 = (j+l2-2)*n_coeffs(1)
                if ( j==1 .and. l2 == 0 ) then
                   pp_coeffs_local(:,1) = 0.0_f64
                   !print*, 'a'
                   !continue
                   !goto 100
                   cycle
                elseif ( j == n_cells(2) .and. l2 == degree2 ) then
                   pp_coeffs_local(:,degree2+1) = 0.0_f64
                   !continue
                   !goto 100
                   cycle
                end if
             case( sll_p_boundary_clamped_clampeddiri )
                offset2 = (j+l2-1)*n_coeffs(1)
                if ( j == n_cells(2) .and. l2 == degree2 ) then
                   pp_coeffs_local(:,degree2+1) = 0.0_f64
                   !continue
                   !goto 100
                   cycle
                end if
             case( sll_p_boundary_clampeddiri_clamped )
                offset2 = (j+l2-2)*n_coeffs(1)
                if ( j==1 .and. l2 == 0 ) then
                   pp_coeffs_local(:,1) = 0.0_f64
                   !continue
                   !goto 100
                   cycle
                end if
                !case ( sll_p_boundary_clampeddiri:sll_p_boundary_clampeddiri_clamped)
                !   offset2 = (j+l2-2)*n_coeffs(1)
             case default
                offset2 = (j+l2-1)*n_coeffs(1)
             end select

             do l1=0,degree1
                select case( self%spline1%boundary_conditions )
                case ( sll_p_boundary_periodic )
                   index = modulo(i+offset1+l1-1,n_coeffs(1))+1
                   !print*, i, j, l1, l2, offset2, index
                   self%spline1%scratch_b(l1+1) = b_coeffs( offset2+index )
                case ( sll_p_boundary_clamped )
                   index = i + offset1 + l1
                   self%spline1%scratch_b(l1+1) = b_coeffs( offset2+index )
                   ! Set to zero for Dirichlet if we are at the boundary
                case ( sll_p_boundary_clamped_clampeddiri )
                   if ( i == n_cells(1) .and. l1 == degree1 ) then
                      self%spline1%scratch_b(l1+1) = 0.0_f64
                   else
                      index = i + offset1 + l1
                      self%spline1%scratch_b(l1+1) = b_coeffs( offset2+index )
                   end if
                case ( sll_p_boundary_clampeddiri_clamped )
                   if ( i == 1 .and. l1 == 0 ) then
                      self%spline1%scratch_b(l1+1) = 0.0_f64
                   else
                      index = i + offset1 + l1
                      self%spline1%scratch_b(l1+1) = b_coeffs( offset2+index )
                   end if
                case ( sll_p_boundary_clampeddiri )
                   if ( i == 1 .and. l1 == 0 ) then
                      self%spline1%scratch_b(l1+1) = 0.0_f64
                   elseif ( i == n_cells(1) .and. l1 == degree1 ) then
                      self%spline1%scratch_b(l1+1) = 0.0_f64
                   else
                      index = i + offset1 + l1
                      self%spline1%scratch_b(l1+1) = b_coeffs( offset2+index )
                   end if
                end select
             end do
             ! For clamped splines, we need to use the boundary coefficients
             select case( self%spline1%boundary_conditions )
             case ( sll_p_boundary_periodic )
                call sll_s_spline_pp_b_to_pp_1d_cella(degree1, self%spline1%poly_coeffs, &
                     self%spline1%scratch_b,pp_coeffs_local(:,l2+1))
             case default
                if ( i > degree1-1 ) then
                   if ( i < upper(1)+1 ) then
                      call sll_s_spline_pp_b_to_pp_1d_cella(degree1, self%spline1%poly_coeffs, &
                           self%spline1%scratch_b,pp_coeffs_local(:,l2+1))
                   else
                      call sll_s_spline_pp_b_to_pp_1d_cella(degree1, &
                           self%spline1%poly_coeffs_boundary_right(:,:,i-upper(1)), &
                           self%spline1%scratch_b,pp_coeffs_local(:,l2+1))
                   end if
                else

                   call sll_s_spline_pp_b_to_pp_1d_cella(degree1, &
                        self%spline1%poly_coeffs_boundary_left(:,:,i), &
                        self%spline1%scratch_b,pp_coeffs_local(:,l2+1))
                end if


             end select
             !call sll_s_spline_pp_b_to_pp_1d_cell(self%spline1, &
             !     self%spline1%scratch_b,pp_coeffs_local(:,l2+1))
             !100          continue             
          end do
          do l1=0,degree1
             self%spline2%scratch_b =  pp_coeffs_local(l1+1,:)
             select case( self%spline2%boundary_conditions )
             case ( sll_p_boundary_periodic )
                call sll_s_spline_pp_b_to_pp_1d_cella(degree2, self%spline2%poly_coeffs, &
                     self%spline2%scratch_b,self%spline2%scratch_p)
             case default
                if ( j > degree2-1 ) then
                   if ( j < upper(2)+1 ) then
                      call sll_s_spline_pp_b_to_pp_1d_cella(degree2, &
                           self%spline2%poly_coeffs, &
                           self%spline2%scratch_b,self%spline2%scratch_p)
                   else
                      call sll_s_spline_pp_b_to_pp_1d_cella(degree2, &
                           self%spline2%poly_coeffs_boundary_right(:,:,j-upper(2)), &
                           self%spline2%scratch_b,self%spline2%scratch_p)
                   end if
                else

                   call sll_s_spline_pp_b_to_pp_1d_cella(degree2, &
                        self%spline2%poly_coeffs_boundary_left(:,:,j), &
                        self%spline2%scratch_b,self%spline2%scratch_p)
                end if


             end select


             !call sll_s_spline_pp_b_to_pp_1d_cell(self%spline2, &
             !     self%spline2%scratch_b,self%spline2%scratch_p)
             do l2=0,degree2
                pp_coeffs( l2*(degree1+1)+l1+1 , cell) = self%spline2%scratch_p(l2+1)
             end do
          end do
       end do
    end do

  end subroutine sll_s_spline_pp_b_to_pp_2d



!!$  !> Convert 2d spline in B form in a cell to spline in pp form with periodic boundary conditions 
!!$  subroutine sll_s_spline_pp_b_to_pp_2d_cell_clamped(spline1,spline2,n_cells, b_coeffs, pp_coeffs,i,j)
!!$    type( sll_t_spline_pp_1d), intent(inout)::  spline1 !< arbitrary degree 1d spline
!!$    type( sll_t_spline_pp_1d), intent(inout)::  spline2 !< arbitrary degree 1d spline
!!$    sll_int32, intent(in)    :: n_cells(2) !< number of gridcells
!!$    sll_real64,intent(in)    :: b_coeffs(n_cells(1)*n_cells(2))   !< coefficients of spline in B-form
!!$    sll_real64,intent(inout) :: pp_coeffs((spline1%degree+1)*(spline2%degree+1),n_cells(1)*n_cells(2))  !< coefficients of spline in pp-form
!!$
!!$    !> convert b-coefficients to pp-coefficients in first dimension
!!$    do j=1, n_cells(2)
!!$       do i=1, n_cells(1)
!!$          do l=0, degree(2)
!!$             spline1%scratch_b=b_coeffs(i-degree1+(j-degp2+l)*n_cells(1):i+(j-degp2+l)*n_cells(1))
!!$
!!$  end subroutine sll_s_spline_pp_b_to_pp_2d_cell_clamped


  !> Convert 2d spline in B form in a cell to spline in pp form with periodic boundary conditions 
  subroutine sll_s_spline_pp_b_to_pp_2d_cell(spline1,spline2,n_cells, b_coeffs, pp_coeffs,i,j)
    type( sll_t_spline_pp_1d), intent(inout)::  spline1 !< arbitrary degree 1d spline
    type( sll_t_spline_pp_1d), intent(inout)::  spline2 !< arbitrary degree 1d spline
    sll_int32, intent(in)    :: n_cells(2) !< number of gridcells
    sll_real64,intent(in)    :: b_coeffs(n_cells(1)*n_cells(2))   !< coefficients of spline in B-form
    sll_real64,intent(inout) :: pp_coeffs((spline1%degree+1)*(spline2%degree+1),n_cells(1)*n_cells(2))  !< coefficients of spline in pp-form
    sll_int32, intent(in)    :: i,j !< indices 
    sll_int32 :: k,l
    sll_int32 :: degree1,degree2,degp1,degp2
    degree1= spline1%degree
    degree2= spline2%degree
    degp1=degree1+1
    degp2=degree2+1
    !> convert b-coefficients in pp-coefficients in first dimension
    if (i>degree1) then
       if(j>degree2) then
          do l=0,degree2
             spline1%scratch_b=b_coeffs(i-degree1+(j-degp2+l)*n_cells(1):i+(j-degp2+l)*n_cells(1))
             call sll_s_spline_pp_b_to_pp_1d_cell(spline1, spline1%scratch_b, pp_coeffs(1+l*degp1:degp1*(l+1),i+n_cells(1)*(j-1)))
          end do
       else 
          !> use of modulo for boundary cells in second dimension 
          do l=0,degree2
             spline1%scratch_b=b_coeffs(i-degree1+modulo(j-degp2+l,n_cells(2))*n_cells(1):i+modulo(j-degp2+l,n_cells(2))*n_cells(1))
             call sll_s_spline_pp_b_to_pp_1d_cell(spline1, spline1%scratch_b, pp_coeffs(1+l*degp1:degp1*(l+1),i+n_cells(1)*(j-1)))
          end do
       end if
    else 
       !> use of modulo for boundary cells in both dimensions 
       do l=0,degree2
          do k=0, degree1
             spline1%scratch_b(k+1)=b_coeffs(modulo(i-degp1+k,n_cells(1))+1+modulo(j-degp2+l,n_cells(2))*n_cells(1))
          end do

          call sll_s_spline_pp_b_to_pp_1d_cell(spline1, spline1%scratch_b, pp_coeffs(1+l*degp1:degp1*(l+1),i+n_cells(1)*(j-1)))
       end do
    end if

    !> convert b-coefficients in pp_coefficients in second dimension
    do k=1, degp1
       do l=1,degp2
          spline2%scratch_b(l)=pp_coeffs(k+degp1*(l-1),i+n_cells(1)*(j-1))
       end do
       call sll_s_spline_pp_b_to_pp_1d_cell(spline2, spline2%scratch_b, spline2%scratch_p)
       do l=1,degp2
          pp_coeffs(k+degp1*(l-1),i+n_cells(1)*(j-1))=spline2%scratch_p(l)
       end do
    end do
  end subroutine sll_s_spline_pp_b_to_pp_2d_cell


  !> Perform two times a 1d hornerschema on the poly_coeffs
  subroutine sll_s_spline_pp_horner_m_2d(self, val, degree, x)
    type( sll_t_spline_pp_2d), intent(in)::  self !< arbitrary degree 2d spline 
    sll_real64, intent(out):: val(:,:) !< array of values
    sll_int32, intent(in)  :: degree(2) !< degree of the spline
    sll_real64, intent(in) :: x(2) !< point at which we evaluate our spline
    sll_int32 :: i

    do i=1, degree(1)+1!size(val,1)
       val(i,1)=sll_f_spline_pp_horner_1d(degree(1), self%spline1%poly_coeffs, x(1), i)
    end do
    do i=1, degree(2)+1!size(val,2)
       val(i,2)=sll_f_spline_pp_horner_1d(degree(2), self%spline2%poly_coeffs, x(2), i)
    end do
  end subroutine sll_s_spline_pp_horner_m_2d

  !> Perform a 2d hornerschema on the pp_coeffs at the indices
  function sll_f_spline_pp_horner_2d(degree, pp_coeffs, x, indices, n_cells) result(res)
    sll_int32,  intent(in) :: degree(2) !< degree of the spline
    sll_real64, intent(in) :: pp_coeffs(:,:)  !< coefficients of spline in pp-form
    sll_real64, intent(in) :: x(2) !< point at which we evaluate our spline
    sll_int32,  intent(in) :: indices(2) !< indices of cell in which is x
    sll_int32,  intent(in) :: n_cells(2) !< number of gridcells
    sll_real64 :: res !< value of the splinefunction at point x
    sll_real64 :: pp_coeffs_1d(degree(2)+1,1)
    sll_int32  :: i
    !> Perform a 1d hornerschema in the first dimension
    do i=0,degree(2)
       pp_coeffs_1d(i+1,1)=sll_f_spline_pp_horner_1d(degree(1), pp_coeffs(1+i*(degree(1)+1):(degree(1)+1)*(i+1),1+(indices(2)-1)*n_cells(1):n_cells(1)*indices(2)), x(1),indices(1))
    end do
    !> Perform a 1d hornerschema in the second dimension
    res=sll_f_spline_pp_horner_1d(degree(2), pp_coeffs_1d, x(2),1)
  end function sll_f_spline_pp_horner_2d

=#
