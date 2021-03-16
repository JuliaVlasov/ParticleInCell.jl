import Statistics: mean
import Base.Threads: @threads, @sync, @spawn, nthreads, threadid

export CloudInCell

struct CloudInCell end

export compute_rho

function compute_rho(p, kernel::CloudInCell, m :: TwoDGrid)

    nx, ny = m.nx, m.ny
    dx, dy = m.dx, m.dy
    rho = zeros(nx+1, ny+1)
    nbpart = size(p.array)[2]

    @inbounds for ipart = 1:nbpart
        xp = p.array[1, ipart] / dx
        yp = p.array[2, ipart] / dy

        i = trunc(Int, xp) + 1
        j = trunc(Int, yp) + 1

        dxp = xp - i + 1
        dyp = yp - j + 1

        a1 = (1 - dxp) * (1 - dyp)
        a2 = dxp * (1 - dyp)
        a3 = dxp * dyp
        a4 = (1 - dxp) * dyp

        ip1 = mod1(i + 1, nx)
        jp1 = mod1(j + 1, ny)

        w = p.array[5, ipart]

        rho[i, j] += a1 * w
        rho[ip1, j] += a2 * w
        rho[ip1, jp1] += a3 * w
        rho[i, jp1] += a4 * w

    end

    for i = 1:nx
       rho[i,1] += rho[i,ny+1]
    end
    for j = 1:ny
       rho[1,j] += rho[nx+1,j]
    end
    
    rho ./= ( dx * dy )

    return rho[1:nx,1:ny] .- mean(rho[1:nx,1:ny])


end

export compute_current!

function compute_current!(jx, jy, m::TwoDGrid, kernel::CloudInCell, p)

    nbpart = size(p.array)[2]
    nx, ny = m.nx, m.ny
    dx, dy = m.dx, m.dy

    fill!(jx, 0)
    fill!(jy, 0)

    scaling = 1 / ( dx * dy )

    ntid = nthreads()
    jxloc = [zero(jx) for _ in 1:ntid]
    jyloc = [zero(jy) for _ in 1:ntid]

    chunks = Iterators.partition(1:nbpart, nbpart÷ntid)

    @sync for chunk in chunks
        @spawn begin
            tid = threadid()
            @inbounds for ipart in chunk
                xp = p.array[1, ipart] / dx
                yp = p.array[2, ipart] / dy

                i = floor(Int, xp) + 1
                j = floor(Int, yp) + 1

                dxp = xp - i + 1
                dyp = yp - j + 1

                a1 = (1 - dxp) * (1 - dyp)
                a2 = dxp * (1 - dyp)
                a3 = dxp * dyp
                a4 = (1 - dxp) * dyp

                w = p.array[5, ipart]

                w1 = p.array[3, ipart] * scaling * w
                w2 = p.array[4, ipart] * scaling * w

                jxloc[tid][i, j] += a1 * w1
                jyloc[tid][i, j] += a1 * w2

                jxloc[tid][i+1, j] += a2 * w1
                jyloc[tid][i+1, j] += a2 * w2

                jxloc[tid][i+1, j+1] += a3 * w1
                jyloc[tid][i+1, j+1] += a3 * w2

                jxloc[tid][i, j+1] += a4 * w1
                jyloc[tid][i, j+1] += a4 * w2
            end
        end
    end

    jx .= reduce(+,jxloc)
    jy .= reduce(+,jyloc)

    for i = 1:nx+1
        jx[i, 1] += jx[i, ny+1]
        jx[i, ny+1] = jx[i, 1]
        jy[i, 1] += jy[i, ny+1]
        jy[i, ny+1] = jy[i, 1]
    end
    for j = 1:ny+1
        jx[1, j] += jx[nx+1, j]
        jx[nx+1, j] = jx[1, j]
        jy[1, j] += jy[nx+1, j]
        jy[nx+1, j] = jy[1, j]
    end

end
