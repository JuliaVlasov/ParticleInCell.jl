export push_v!

function push_v!(p, dt)

    for ipart = 1:p.nbpart

        v1 = p.vit[ipart, 1]
        v2 = p.vit[ipart, 2]
        e1 = p.epx[ipart]
        e2 = p.epy[ipart]

        v1 = v1 + 0.5dt * e1
        v2 = v2 + 0.5dt * e2

        tantheta = 0.5dt * p.bpz[ipart]
        sintheta = 2 * tantheta / (1 + tantheta * tantheta)

        v1 = v1 + v2 * tantheta
        v2 = v2 - v1 * sintheta
        v1 = v1 + v2 * tantheta

        p.vit[ipart, 1] = v1 + 0.5dt * e1
        p.vit[ipart, 2] = v2 + 0.5dt * e2

    end

end

export push_x!

function push_x!(p, mesh, dt)

    p.pos .+= p.vit .* dt
    update_cells!(p, mesh)

end
