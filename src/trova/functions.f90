!> Computes the sum of tensor values within specified longitude and latitude bounds.
!!
!! This subroutine iterates over the tensor and sums the values in the third dimension
!! for each grid cell defined by the longitude and latitude arrays.
!!
!! @param result (real(8), intent(out)) The output array containing the summed values.
!! @param tensor (real(8), intent(in)) The input tensor containing the data to be summed.
!! @param lon (real(8), intent(in)) The array of longitude values defining the grid cells.
!! @param lat (real(8), intent(in)) The array of latitude values defining the grid cells.
!! @param numPdY (integer, intent(in)) The number of grid points in the Y direction.
!! @param numPdX (integer, intent(in)) The number of grid points in the X direction.
!! @param nlen (integer, intent(in)) The length of the first dimension of the tensor.
!! @param npart (integer, intent(in)) The length of the second dimension of the tensor.

subroutine K_dq(result, tensor, lon, lat, numPdY, numPdX, nlen, npart)
    implicit none
    integer, intent(in) :: numPdY, numPdX, nlen, npart
    real(8), intent(in) :: lat(numPdY+1, numPdX+1)
    real(8), intent(in) :: lon(numPdY+1, numPdX+1)
    real(8), intent(in) :: tensor(nlen, npart, 3)
    real(8), intent(out) :: result(numPdY, numPdX)
    integer :: i, j, k, n

    result = 0.0

    do i = 1, nlen
        do j = 1, npart
            if (tensor(i, j, 1) /= -999.9 .and. tensor(i, j, 2) /= -999.9) then
                do k = 1, numPdY
                    do n = 1, numPdX
                        if (tensor(i, j, 1) > lon(k, n) .and. tensor(i, j, 1) < lon(k+1, n+1) .and. &
                            tensor(i, j, 2) > lat(k, n) .and. tensor(i, j, 2) < lat(k+1, n+1)) then
                            result(k, n) = result(k, n) + tensor(i, j, 3)
                            exit 
                        endif
                    enddo
                enddo
            endif
        enddo
    enddo
end subroutine K_dq

!> Computes the sum of tensor values within specified longitude, latitude, and vertical layer bounds.
!!
!! This subroutine iterates over the tensor and sums the values in the third dimension
!! for each grid cell defined by the longitude and latitude arrays, and within the specified vertical layer bounds.
!!
!! @param result (real(8), intent(out)) The output array containing the summed values.
!! @param tensor (real(8), intent(in)) The input tensor containing the data to be summed.
!! @param z0 (real(8), intent(in)) The lower bound of the vertical layer.
!! @param z1 (real(8), intent(in)) The upper bound of the vertical layer.
!! @param lon (real(8), intent(in)) The array of longitude values defining the grid cells.
!! @param lat (real(8), intent(in)) The array of latitude values defining the grid cells.
!! @param numPdY (integer, intent(in)) The number of grid points in the Y direction.
!! @param numPdX (integer, intent(in)) The number of grid points in the X direction.
!! @param nlen (integer, intent(in)) The length of the first dimension of the tensor.
!! @param npart (integer, intent(in)) The length of the second dimension of the tensor.

subroutine K_dq_layers(result, tensor, z0, z1, lon, lat, numPdY, numPdX, nlen, npart)
    implicit none
    integer, intent(in) :: numPdY, numPdX, nlen, npart
    real(8), intent(in) :: lat(numPdY+1, numPdX+1)
    real(8), intent(in) :: lon(numPdY+1, numPdX+1)
    real(8), intent(in) :: tensor(nlen, npart, 4)
    real(8), intent(in) :: z0, z1
    real(8), intent(out) :: result(numPdY, numPdX)
    integer :: i, j, k, n

    result = 0.0

    do i = 1, nlen
        do j = 1, npart
            if (tensor(i, j, 1) /= -999.9 .and. tensor(i, j, 2) /= -999.9) then
                if (tensor(i, j, 4) > z0 .and. tensor(i, j, 4) <= z1) then
                    do k = 1, numPdY
                        do n = 1, numPdX
                            if (tensor(i, j, 1) > lon(k, n) .and. tensor(i, j, 1) < lon(k+1, n+1) .and. &
                                tensor(i, j, 2) > lat(k, n) .and. tensor(i, j, 2) < lat(k+1, n+1)) then
                                result(k, n) = result(k, n) + tensor(i, j, 3)
                                exit
                            endif
                        enddo
                    enddo
                endif
            endif
        enddo
    enddo
end subroutine K_dq_layers

!> Computes the average of tensor values within specified longitude and latitude bounds.
!!
!! This subroutine iterates over the tensor and computes the average of the values in the third dimension
!! for each grid cell defined by the longitude and latitude arrays.
!!
!! @param result (real(8), intent(out)) The output array containing the averaged values.
!! @param tensor (real(8), intent(in)) The input tensor containing the data to be averaged.
!! @param lon (real(8), intent(in)) The array of longitude values defining the grid cells.
!! @param lat (real(8), intent(in)) The array of latitude values defining the grid cells.
!! @param numPdY (integer, intent(in)) The number of grid points in the Y direction.
!! @param numPdX (integer, intent(in)) The number of grid points in the X direction.
!! @param nlen (integer, intent(in)) The length of the first dimension of the tensor.
!! @param npart (integer, intent(in)) The length of the second dimension of the tensor.

subroutine K_dq_por(result, tensor, lon, lat, numPdY, numPdX, nlen, npart)

    implicit none
    integer, intent(in) :: numPdY, numPdX, nlen, npart
    real(8), intent(in) :: lat(numPdY+1, numPdX+1)
    real(8), intent(in) :: lon(numPdY+1, numPdX+1)
    real(8), intent(in) :: tensor(nlen, npart, 3)
    real(8), intent(out) :: result(numPdY, numPdX)
    real(8) :: count_part(numPdY, numPdX)
    integer :: i, j, k, n

    result = 0.0
    count_part = 0.0

    do i = 1, nlen
        do j = 1, npart
            if (tensor(i, j, 1) /= -999.9 .and. tensor(i, j, 2) /= -999.9) then
                do k = 1, numPdY
                    do n = 1, numPdX
                        if (tensor(i, j, 1) > lon(k, n) .and. tensor(i, j, 1) < lon(k+1, n+1) .and. &
                            tensor(i, j, 2) > lat(k, n) .and. tensor(i, j, 2) < lat(k+1, n+1)) then
                            result(k, n) = result(k, n) + tensor(i, j, 3)
                            count_part(k, n) = count_part(k, n) + 1
                            exit
                        endif
                    enddo
                enddo
            endif
        enddo
    enddo

    do k = 1, numPdY
        do n = 1, numPdX
            if (count_part(k, n) /= 0.0) then
                result(k, n) = result(k, n) / count_part(k, n)
            else
                result(k, n) = 0.0
            endif
        enddo
    enddo
end subroutine K_dq_por

!> Computes the sum of tensor values within specified thresholds.
!!
!! This subroutine iterates over the tensor and computes the sum of the values in the third dimension
!! for each grid cell defined by the matrix and matrix_ind arrays, applying the specified threshold.
!!
!! @param result (real(8), intent(out)) The output array containing the summed values.
!! @param matrix (real(8), intent(in)) The input tensor containing the data to be summed.
!! @param matrix_ind (real(8), intent(in)) The input tensor containing the indicator values.
!! @param threshold (real) The threshold value for filtering.
!! @param npart (integer) The number of parts in the tensor.
!! @param ntime (integer) The number of time steps in the tensor.

subroutine K_dq_So(result,matrix,matrix_ind,threshold,npart, ntime)

    real *8, intent(in) :: matrix(npart, ntime)
    real *8, intent(in) :: matrix_ind(npart, ntime)
    real *8 :: matrix_(npart, ntime)
    real *8, intent(out) :: result(npart, ntime)
    integer :: i,j,k,npart
    real threshold 
    
    result = 0.0
      
    do i=1,npart
        do j=1,ntime
            if (int(matrix(i,j))==int(-999.9)) then
                matrix_(i,j)=0.0
            else
               matrix_(i,j)=matrix(i,j)
            endif
        enddo
    enddo
    
    do i=1,npart
        if (matrix_(i,1).lt.0) then
            result(i,1)=0.0
        else
            result(i,1)=matrix_(i,1)
        endif 
         
        do j= 2, ntime
            if (matrix_(i,j).ge.0) then
               result(i,j)=matrix_(i,j)
            else if (matrix_(i,j) .lt. 0 .and.  matrix_ind(i,j) .eq. 1) then
               result(i,j)=0.0
            else if (threshold .lt. matrix_(i,j) .and.   matrix_(i,j).lt. 0) then
               result(i,j)=0.0
            else
               suma=sum(result(i,1:j))
                do k=1,j
                    if (suma/=0) then
                       aux=result(i,k)-(result(i,k)/suma)*abs(matrix_(i,j))
                    else
                       aux=0
                    endif
                
                    if (aux .lt. 0) then
                       result(i,k)=0
                    else
                       result(i,k)=aux
                    endif
                enddo
            endif
        enddo
    enddo
end subroutine K_dq_So

!> Filters particles based on a specified threshold.
!!
!! This subroutine filters particles based on a specified threshold and step value (paso).
!! It updates the output array with the filtered particles and counts the number of filtered particles.
!!
!! @param output (real, intent(out)) The output array containing the filtered particles.
!! @param count_part (integer, intent(out)) The count of filtered particles.
!! @param matrix (real, intent(in)) The input matrix containing particle data.
!! @param matrix_ref (real, intent(in)) The reference matrix for filtering.
!! @param paso (integer, intent(in)) The step value indicating the direction of the simulation (-1 for backward, 1 for forward).
!! @param threshold (real, intent(in)) The threshold value for filtering.
!! @param numP (integer, intent(in)) The number of particles.

subroutine filter_part2(output,count_part,matrix,matrix_ref, paso, threshold, numP)
  
    real, intent(in) :: matrix(numP,11)
    real, intent(in) :: matrix_ref(numP, 11)
    integer, intent(in) :: paso
    real, intent(in) :: threshold 
    real, intent(out) :: output(numP, 11)
    integer, intent(out) :: count_part
    integer :: i, numP
   
    output(:,:)=-999.9
    count_part=0
    if (paso==-1)then
        do i=1, numP
            if (matrix_ref(i,4) .le. threshold .and. matrix_ref(i,4) /= -999.9) then
               output(i,:)=matrix(i,:)
               count_part=count_part+1
            endif
        enddo
    endif
    
    if (paso==1)then 
        do i=1, numP
            if (matrix_ref(i,3) .ge. threshold) then
               output(i,:)=matrix(i,:)
            endif
        enddo
    endif
end subroutine filter_part2

!> Filters particles based on height within specified layers.
!!
!! This subroutine filters particles based on a specified height range (lowerlayer to upperlayer) and step value (paso).
!! It updates the output array with the filtered particles and counts the number of filtered particles.
!!
!! @param output (real, intent(out)) The output array containing the filtered particles.
!! @param count_part (integer, intent(out)) The count of filtered particles.
!! @param matrix (real, intent(in)) The input matrix containing particle data.
!! @param matrix_ref (real, intent(in)) The reference matrix for filtering.
!! @param paso (integer, intent(in)) The step value indicating the direction of the simulation (-1 for backward, 1 for forward).
!! @param lowerlayer (real, intent(in)) The lower bound of the height range for filtering.
!! @param upperlayer (real, intent(in)) The upper bound of the height range for filtering.
!! @param numP (integer, intent(in)) The number of particles.

subroutine filter_part_by_height(output,count_part,matrix,matrix_ref, paso, lowerlayer, upperlayer, numP)
  
    real, intent(in) :: matrix(numP,11)
    real, intent(in) :: matrix_ref(numP, 11)
    integer, intent(in) :: paso
    real, intent(in) :: lowerlayer, upperlayer 
    real, intent(out) :: output(numP, 11)
    integer, intent(out) :: count_part
    integer :: i, numP
      
    output(:,:)=-999.9
    count_part=0
    
    if (paso==-1)then
        do i=1, numP
            if (matrix_ref(i,5) .ge. lowerlayer .and.  matrix_ref(i,5) .le. upperlayer .and. matrix_ref(i,5) /= -999.9) then
               output(i,:)=matrix(i,:)
               count_part=count_part+1
            endif
        enddo
    endif
    
    if (paso==1)then 
        do i=1, numP
            if (matrix_ref(i,5) .ge.  lowerlayer .and.  matrix_ref(i,5) .le. upperlayer .and. matrix_ref(i,5) /= -999.9 ) then
               output(i,:)=matrix(i,:)
               count_part=count_part+1
            endif
        enddo
    endif
end subroutine filter_part_by_height

!> Reads a FLEXPART binary output file and filters particles by domain and/or height.
!!
!! This subroutine reads particle data from a FLEXPART unformatted binary file and stores
!! the output in a matrix. Optionally, it filters the particles based on:
!! - a geographical domain (longitude and latitude bounds),
!! - a vertical height threshold.
!!
!! Particles that do not meet the specified filtering criteria are skipped.
!! The output array is initialized with -999 and filled with valid particle data.
!!
!! @param output_ (real*8, intent(out)) A 2D array of shape (nparts, 11) where each row contains:
!!   - 1: Particle ID (npoint)
!!   - 2: Longitude (xlon)
!!   - 3: Latitude (ylat)
!!   - 4: Specific humidity (qvi)
!!   - 5: Vertical height (ztra1)
!!   - 6: Topography (topo)
!!   - 7: Density (rhoi)
!!   - 8: PBL height (hmixi)
!!   - 9: Tropopause height (tri)
!!   - 10: Temperature (tti)
!!   - 11: Particle mass (xmass)
!! @param filename (character, intent(in)) Path to the FLEXPART binary file to read.
!! @param nparts (integer(kind=4), intent(in)) Maximum number of particles to read.
!! @param x_l (real, intent(in)) Left (minimum) longitude of the geographic domain.
!! @param y_l (real, intent(in)) Bottom (minimum) latitude of the geographic domain.
!! @param x_r (real, intent(in)) Right (maximum) longitude of the geographic domain.
!! @param y_r (real, intent(in)) Top (maximum) latitude of the geographic domain.
!! @param limit_domain (integer, intent(in)) Flag to apply domain filter (1 = apply, 0 = ignore).
!! @param limit_height (integer, intent(in)) Flag to apply height filter (1 = apply, 0 = ignore).
!! @param value_height (real, intent(in)) Maximum allowed height if limit_height = 1.

subroutine read_binary_file_flexpart_height(output_, filename, nparts, x_l, y_l, x_r, y_r, limit_domain, limit_height, value_height)

    implicit none
    integer, parameter :: rk = kind(1.0)
    integer(kind=4), intent(in) :: nparts
    real*8, intent(out) :: output_(nparts, 11)
    character(*), intent(in) :: filename
    real, intent(in) :: x_l, y_l, x_r, y_r, value_height
    integer, intent(in) :: limit_domain, limit_height

    integer :: i, cant, ios
    integer :: unitpartout, itime
    integer :: npoint, itramem
    real :: xlon, ylat, ztra1, topo, pvi, qvi, rhoi, hmixi, tri, tti, xmass

    cant = 1
    unitpartout = 1

    output_(:,:) = -999

    open(unitpartout, file=filename, form='unformatted', IOSTAT=ios)
    if (ios /= 0) then
        print *, 'Error openning the file:', filename
        return
    endif

    read(unitpartout) itime

    do i = 1, nparts
        read(unitpartout, IOSTAT=ios) npoint, xlon, ylat, ztra1, itramem, topo, pvi, qvi, rhoi, hmixi, tri, tti, xmass
        if (ios /= 0) exit

        if (limit_height == 1) then
            if (ztra1 > value_height) cycle
        endif

        ! Filtro por dominio si se activa
        if (limit_domain == 1) then
            if (xlon < x_l .or. xlon > x_r) cycle
            if (ylat < y_l .or. ylat > y_r) cycle
        endif

        output_(cant,1) = npoint   ! parcel id
        output_(cant,2) = xlon     ! longitude
        output_(cant,3) = ylat     ! latitude
        output_(cant,4) = qvi      ! specific humidity
        output_(cant,5) = ztra1    ! vertical height
        output_(cant,6) = topo     ! topography
        output_(cant,7) = rhoi     ! density
        output_(cant,8) = hmixi    ! PBL height
        output_(cant,9) = tri      ! tropopause height
        output_(cant,10) = tti     ! temperature
        output_(cant,11) = xmass   ! mass

        cant = cant + 1
    end do

    close(unitpartout)

end subroutine read_binary_file_flexpart_height