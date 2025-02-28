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
end subroutine

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
    
end subroutine

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

end subroutine

!> Reads a binary file and extracts relevant data within specified domain limits.
!!
!! This subroutine reads a binary file and extracts relevant data, storing it in the output array.
!! If domain limits are specified, only data within those limits is included in the output.
!!
!! @param output_ (real(8), intent(out)) The output array containing the extracted data.
!! @param filename (character(500), intent(in)) The name of the binary file to read.
!! @param nparts (integer, intent(in)) The number of parts in the binary file.
!! @param x_l (real, intent(in)) The lower X-coordinate of the domain.
!! @param y_l (real, intent(in)) The lower Y-coordinate of the domain.
!! @param x_r (real, intent(in)) The upper X-coordinate of the domain.
!! @param y_r (real, intent(in)) The upper Y-coordinate of the domain.
!! @param limit_domian (integer, intent(in)) Whether to limit the domain (1 for yes, 0 for no).

subroutine read_binary_file(output_,filename, nparts,x_l,y_l, x_r,y_r, limit_domian)

    integer,parameter :: rk=kind(1.0)
    real(rk)         :: b1
    integer          :: id1,i,k,idx,index_id, ty
    integer(kind=4)  :: n1,nparts, cant
    integer bytes,aux_bytes
    real *8  :: matrix(nparts,13)
    real *8 :: output(nparts,11)
    real *8, intent(out) :: output_(nparts,11)
    character(500) :: filename
    real, intent(in) :: x_l,y_l, x_r,y_r
    
    open(1,file=filename,access="stream",form="unformatted")
    bytes=(nparts+1)*60 +12
    read(1, pos=1)ty

    if (ty==4)then
       i=17
       aux_bytes=69
    else   
       i=25
       aux_bytes=77
    endif
   
    index_id=1
    do while (i<bytes-60)
	    read(1,pos=i) id1 
	    matrix(index_id,1)=id1
	    idx=2
        do k=i+4,aux_bytes-1,4
		    read(1,pos=k) n1
		    read(1,pos=k) b1 
            matrix(index_id,idx)=b1
		    idx=idx+1
        end do
	    i=i+60
	    aux_bytes=aux_bytes+60
	    index_id=index_id+1
    end do
    
    output(:,1)=matrix(:,1)
    output(:,2)=matrix(:,2)
    output(:,3)=matrix(:,3)
    output(:,4)=matrix(:,8)
    output(:,5)=matrix(:,4)
    output(:,6)=matrix(:,6)
    output(:,7)=matrix(:,9)
    output(:,8)=matrix(:,10)
    output(:,9)=matrix(:,11)
    output(:,10)=matrix(:,12)
    output(:,11)=matrix(:,13)
 
    output_(:,:)=-999.
    
   if (limit_domian == 1) then
        cant=1
        do j=1, nparts-1
           if (output(j,2) .ge. x_l .and. output(j,2) .le. x_r .and. &
              output(j,3) .ge. y_l .and. output(j,3) .le. y_r) then
              output_(cant,:)=output(j,:)
              cant=cant+1
           endif
        enddo
    else
        output_(:,:)=output(:,:)
    endif
    close(1)
return
end subroutine read_binary_file

!> Determines the size of a file in bytes.
!!
!! This subroutine determines the size of a specified file in bytes and returns the size.
!!
!! @param bytes (integer, intent(out)) The size of the file in bytes.
!! @param filename (character(500), intent(in)) The name of the file to check.

subroutine len_file(bytes, filename)

    character(500) :: filename
    character(len=256) :: message
    integer :: ios,x
    logical :: foundit
    integer, intent(out) :: bytes
    inquire(file=filename,exist=foundit,size=x,iostat=ios,iomsg=message)
    bytes=x
end subroutine

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
end subroutine

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

subroutine filter_part(output,count_part,matrix,matrix_ref, paso, threshold, numP)
  
    real, intent(in) :: matrix(numP,4)
    real, intent(in) :: matrix_ref(numP, 4)
    integer, intent(in) :: paso
    real, intent(in) :: threshold 
    real, intent(out) :: output(numP, 4)
    integer, intent(out) :: count_part
    integer :: i, numP
   
    output(:,:)=-999.9
    count_part=0
    if (paso==-1)then
        do i=1, numP
            if (matrix_ref(i,3) .le. threshold .and. matrix_ref(i,3) /= -999.9) then
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
end subroutine

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
end subroutine

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
end subroutine

!> Searches for rows in the matrix that match the values in the list.
!!
!! This subroutine searches for rows in the matrix that match the values in the list.
!! It updates the output array with the matching rows.
!!
!! @param output (real, intent(out)) The output array containing the matching rows.
!! @param matrix (real, intent(in)) The input matrix containing data.
!! @param lista (real, intent(in)) The list of values to search for.
!! @param len_lista (integer, intent(in)) The length of the list.
!! @param numP (integer, intent(in)) The number of rows in the matrix.

subroutine search_row(output,matrix,lista,len_lista,numP)
   
    real, intent(in) :: matrix(numP,11)
    real, intent(in) :: lista(len_lista)
    !integer, intent(in) :: lista(len_lista)
    real, intent(out) :: output(len_lista, 11)
    integer :: i,j, numP
    
    output(:,:)=-999.9
   
    do i=1, len_lista
        do j=1, numP
            if (int(matrix(j,1))==int(lista(i))) then
              output(i,:) = matrix(j,:)
           endif 
        enddo
    enddo
        
end subroutine

!> Determines the indices of elements in value_mascara that match value_mask.
!!
!! This subroutine searches through the value_mascara array and sets the corresponding
!! indices in the vector array where the elements match the value_mask.
!!
!! @param vector (integer, intent(out)) The output array containing the indices of matching elements.
!! @param value_mascara (integer, intent(in)) The input array to search through.
!! @param value_mask (integer, intent(in)) The value to match in the value_mascara array.
!! @param len_value_mascara (integer, intent(in)) The length of the value_mascara array.

subroutine determined_id(vector, value_mascara,value_mask,len_value_mascara)

    integer, intent(in) :: value_mascara(len_value_mascara)
    integer, intent(out) :: vector(len_value_mascara)
    integer :: i,j, len_value_mascara, value_mask
    
    vector = -999

    do j=1,len_value_mascara
        if (value_mascara(j)==value_mask) then
           vector(j)=j-1
        endif
    enddo
end subroutine

!> Computes the difference between two matrices based on a specified step value (paso).
!!
!! This subroutine computes the difference between two matrices (matrix1 and matrix2) based on a specified step value (paso).
!! It updates the output array with the computed differences.
!!
!! @param output (real(8), intent(out)) The output array containing the computed differences.
!! @param matrix1 (real(8), intent(in)) The first input matrix.
!! @param matrix2 (real(8), intent(in)) The second input matrix.
!! @param paso (real, intent(in)) The step value indicating the direction of the simulation (-1 for backward, 1 for forward).
!! @param dx (integer, intent(in)) The number of rows in the matrices.
!! @param dy (integer, intent(in)) The number of columns in the matrices.

subroutine Kdif(output, matrix1, matrix2, paso, dx, dy)

    integer :: dx,dy, i
    real *8, intent(in):: matrix1(dx,dy)
    real *8, intent(in):: matrix2(dx,dy)
    real, intent(in) :: paso
    real *8, intent(out):: output(dx,dy-1)
    
    output(:,:)=-999.9

    if (paso == -1.) then
        do i=1,dx
            if (int(matrix1(i,1)) /= int(-999.9) .and. int(matrix2(i,1)) /= int(-999.9)) then 
               output(i,3)=matrix2(i,4)-matrix1(i,4)
               output(i,2)=matrix1(i,3)
               output(i,1)=matrix1(i,2)
               output(i,4)=matrix1(i,5)
            endif
        enddo
    endif 

    if (paso == 1.) then
        do i=1,dx
            if (int(matrix2(i,1)) /= int(-999.9) .and. int(matrix1(i,1)) /= int(-999.9)) then 
                output(i,3)=matrix2(i,4)-matrix1(i,4)
                output(i,2)=matrix2(i,3)
                output(i,1)=matrix2(i,2)
                output(i,4)=matrix2(i,5)
            endif
        enddo
    endif 
end subroutine