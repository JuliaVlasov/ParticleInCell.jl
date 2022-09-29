export read_particles

function read_particles( filename, mesh :: Mesh )

    nbpart = countlines(filename)
    println( " nbpart   : ", nbpart )
    println( " filename : ", filename )

    nx, ny = mesh.nx, mesh.ny
    dx, dy = mesh.dx, mesh.dy
    dimx   = mesh.xmax - mesh.xmin
    dimy   = mesh.ymax - mesh.ymin

    weight    = (dimx * dimy) / nbpart
    particles = Particles( nbpart, weight)

    open(filename) do f

        for (k,line) in enumerate(eachline(f))
            str_values = split(line)
            ix = parse(Int32,   str_values[1])
            iy = parse(Int32,   str_values[2])
            dpx = parse(Float64, str_values[3])
            dpy = parse(Float64, str_values[4])
            particles.v[1,k] = parse(Float64, str_values[5])
            particles.v[2,k] = parse(Float64, str_values[6])
            particles.x[1,k] = (dpx+ix) * dx
            particles.x[2,k] = (dpy+iy) * dy
        end

    end

    particles

end

