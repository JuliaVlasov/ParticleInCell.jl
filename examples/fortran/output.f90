module output_m

use hdf5
use mesh_fields_m

implicit none

contains

subroutine write_data(istep, fields)

integer,           intent(in) :: istep
type(fields_3d_t), intent(in) :: fields

integer             :: error
integer(hid_t)      :: file_id, dataset_id, dataspace_id
integer(hsize_t)    :: data_dims(3)
character(len=4)    :: cplot
integer             :: nx, ny, nz


!Initialize FORTRAN interface.
CALL h5open_f (error)
!Create a new file using default properties.
call int2string(istep, cplot)
call H5Fcreate_f("fields-"//cplot//".h5", H5F_ACC_TRUNC_F, file_id, error)

!Write the data.
nx = fields%mesh%nx
ny = fields%mesh%ny
nz = fields%mesh%nz
data_dims = [nx+1, ny+1, nz+1]

call H5Screate_simple_f(3, data_dims, dataspace_id, error);

call H5Dcreate_f(file_id, "/ex", H5T_NATIVE_DOUBLE, &
		 dataspace_id, dataset_id, error)
call H5Dwrite_f(dataset_id, H5T_NATIVE_DOUBLE, fields%e(1,:,:,:), data_dims, error)
call H5Dclose_f(dataset_id,error)
call H5Dcreate_f(file_id, "/ey", H5T_NATIVE_DOUBLE, &
		 dataspace_id, dataset_id, error)
call H5Dwrite_f(dataset_id, H5T_NATIVE_DOUBLE, fields%e(2,:,:,:), data_dims, error)
call H5Dclose_f(dataset_id,error)
call H5Dcreate_f(file_id, "/ez", H5T_NATIVE_DOUBLE, &
		 dataspace_id, dataset_id, error)
call H5Dwrite_f(dataset_id, H5T_NATIVE_DOUBLE, fields%e(3,:,:,:), data_dims, error)
call H5Dclose_f(dataset_id,error)
call H5Dcreate_f(file_id, "/rho", H5T_NATIVE_DOUBLE, &
		 dataspace_id, dataset_id, error)
call H5Dwrite_f(dataset_id, H5T_NATIVE_DOUBLE, fields%rho, data_dims, error)
call H5Dclose_f(dataset_id,error)

call H5Sclose_f(dataspace_id,error);

!Terminate access to the file.
CALL h5fclose_f(file_id, error)

!Close FORTRAN interface.
CALL h5close_f(error)


open(newunit=file_id,file="fields-"//cplot//".xmf")
write(file_id,"(a)")"<?xml version='1.0' ?>"
write(file_id,"(a)")"<!DOCTYPE Xdmf SYSTEM 'Xdmf.dtd' []>"
write(file_id,"(a)")"<Xdmf xmlns:xi='http://www.w3.org/2003/XInclude' Version='2.2'>"
write(file_id,"(a)")"<Domain>"
write(file_id,"(a)")"<Grid Name='mesh' GridType='Uniform'>"
write(file_id,"(a,3i5,a)")"<Topology TopologyType='3DCoRectMesh' NumberOfElements='", &
                           nz+1,ny+1,nx+1,"'/>"
write(file_id,"(a)")"<Geometry GeometryType='ORIGIN_DXDYDZ'>"
write(file_id,"(a)")"<DataItem Dimensions='3' NumberType='Float' Format='XML'>"
write(file_id,"(3f12.5)") fields%mesh%xmin, fields%mesh%ymin, fields%mesh%zmin
write(file_id,"(a)")"</DataItem>"
write(file_id,"(a)")"<DataItem Dimensions='3' NumberType='Float' Format='XML'>"
write(file_id,"(3f12.5)") fields%mesh%dx, fields%mesh%dy, fields%mesh%dz
write(file_id,"(a)")"</DataItem>"
write(file_id,"(a)")"</Geometry>"

write(file_id,"(a)")"<Attribute Name='ex' AttributeType='Scalar' Center='Node'>"
write(file_id,"(a,3i5,a)")"<DataItem Dimensions='",nz+1,ny+1,nx+1, &
                             "' NumberType='Float' Precision='8' Format='HDF'>"
write(file_id,"(a)")"fields-"//cplot//".h5:/ex"
write(file_id,"(a)")"</DataItem>"
write(file_id,"(a)")"</Attribute>"
write(file_id,"(a)")"<Attribute Name='ey' AttributeType='Scalar' Center='Node'>"
write(file_id,"(a,3i5,a)")"<DataItem Dimensions='",nz+1,ny+1,nx+1, &
                             "' NumberType='Float' Precision='8' Format='HDF'>"
write(file_id,"(a)")"fields-"//cplot//".h5:/ey"
write(file_id,"(a)")"</DataItem>"
write(file_id,"(a)")"</Attribute>"
write(file_id,"(a)")"<Attribute Name='ez' AttributeType='Scalar' Center='Node'>"
write(file_id,"(a,3i5,a)")"<DataItem Dimensions='",nz+1,ny+1,nx+1, &
                             "' NumberType='Float' Precision='8' Format='HDF'>"
write(file_id,"(a)")"fields-"//cplot//".h5:/ez"
write(file_id,"(a)")"</DataItem>"
write(file_id,"(a)")"</Attribute>"

write(file_id,"(a)")"<Attribute Name='rho' AttributeType='Scalar' Center='Node'>"
write(file_id,"(a,3i5,a)")"<DataItem Dimensions='",nz+1,ny+1,nx+1, &
                             "' NumberType='Float' Precision='8' Format='HDF'>"
write(file_id,"(a)")"fields-"//cplot//".h5:/rho"
write(file_id,"(a)")"</DataItem>"
write(file_id,"(a)")"</Attribute>"

write(file_id,"(a)")"</Grid>"
write(file_id,"(a)")"</Domain>"
write(file_id,"(a)")"</Xdmf>"
close(file_id)

end subroutine write_data

!> Convert an integer < 9999 to a 4 characters string
subroutine int2string( istep, cstep )

integer         , intent(in ) :: istep   !< input integer
character(len=*), intent(out) :: cstep   !< output string

integer          :: l
character(len=8) :: str_fmt

l = len(cstep)

if ( istep >= 0  .and. istep < 10**l) then
   str_fmt="(I0"//char(l+48)//"."//char(l+48)//")"
   write(cstep,str_fmt) istep
else
   print*, 'int2string', 'index is negative or too big'
   print*, 'index =', istep, ' cstep length = ', l
   cstep = 'xxxx'
end if

end subroutine int2string


end module output_m
