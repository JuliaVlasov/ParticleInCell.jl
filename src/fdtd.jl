export FDTD

struct FDTD

    ex::Array{Float64,2}
    ey::Array{Float64,2}
    bz::Array{Float64,2}

    function FDTD(mesh)

        nx, ny = mesh.nx, mesh.ny
        ex = zeros(nx, ny + 1)
        ey = zeros(nx + 1, ny)
        bz = zeros(nx, ny)

        new(ex, ey, bz)

    end

end

export faraday!

function faraday!(fdtd::FDTD, m::TwoDGrid, dt)

    nx, ny = m.nx, m.ny
    dx, dy = m.dx, m.dy

    for i = 1:nx, j = 1:ny
        dex_dy = (fdtd.ex[i, mod1(j + 1, ny)] - fdtd.ex[i, j]) / dy
        dey_dx = (fdtd.ey[mod1(i + 1, nx), j] - fdtd.ey[i, j]) / dx
        fdtd.bz[i, j] += dt * (dex_dy - dey_dx)
    end

end

export ampere_maxwell!

function ampere_maxwell!(fdtd::FDTD, m::TwoDGrid, jx, jy, dt)

    nx, ny = m.nx, m.ny
    dx, dy = m.dx, m.dy

    for i = 1:nx, j = 1:ny+1
        dbz_dy = (fdtd.bz[i, mod1(j, ny)] - fdtd.bz[i, mod1(j - 1, ny)]) / dy
        fdtd.ex[i, j] += dt * dbz_dy - dt * 0.5 * (jx[i, j] + jx[mod1(i + 1, nx), j])
    end

    for i = 1:nx+1, j = 1:ny
        dbz_dx = (fdtd.bz[mod1(i, nx), mod1(j, ny)] - fdtd.bz[mod1(i - 1, nx), j]) / dx
        fdtd.ey[i, j] -= dt * dbz_dx - dt * 0.5 * (jy[i, j] + jy[i, mod1(j + 1, ny)])
    end

end

export update_fields!

function update_fields!(ex, ey, bz, m::TwoDGrid, fdtd::FDTD)

    nx, ny = m.nx, m.ny

    for i = 1:nx+1, j = 1:ny+1
        ex[i, j] = 0.5 * (fdtd.ex[mod1(i - 1, nx), j] + fdtd.ex[mod1(i, nx), j])
        ey[i, j] = 0.5 * (fdtd.ey[i, mod1(j - 1, ny)] + fdtd.ey[i, mod1(j, ny)])
        bz[i, j] =
            0.25 * (
                fdtd.bz[mod1(i - 1, nx), mod1(j - 1, ny)] +
                fdtd.bz[mod1(i, nx), mod1(j - 1, ny)] +
                fdtd.bz[mod1(i - 1, nx), mod1(j, ny)] +
                fdtd.bz[mod1(i, nx), mod1(j, ny)]
            )
    end

end




export compute_energy

compute_energy(fdtd::FDTD, m::TwoDGrid) = 0.5 * log(sum(fdtd.ex .^ 2) * m.dx * m.dy)
