struct Fields

   nx
   ny
   ex
   ey
   bz
   jx
   jy
   r0
   r1

   function Fields(nx, ny)

       ex = OffsetArray{Float64}(undef, 0:nx-1,0:ny)
       ey = OffsetArray{Float64}(undef, 0:nx,0:ny-1)
       bz = OffsetArray{Float64}(undef, 0:nx-1,0:ny-1)
       jx = OffsetArray{Float64}(undef, 0:nx-1,0:ny)
       jy = OffsetArray{Float64}(undef, 0:nx,0:ny-1)
       r0 = OffsetArray{Float64}(undef, 0:nx,0:ny) 
       r1 = OffsetArray{Float64}(undef, 0:nx,0:ny)

       new( nx, ny, ex, ey, bz, jx, jy, r0, r1)

   end

end
