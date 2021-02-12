export push_v!

function push_v!(p, dt)

    for ipart = 1:p.nbpart

        v1 = p.vit[1, ipart]
        v2 = p.vit[2, ipart]
        e1 = p.ebp[1, ipart]
        e2 = p.ebp[2, ipart]

        v1 += 0.5dt * e1
        v2 += 0.5dt * e2

        tantheta = 0.5dt * p.ebp[3, ipart]
        sintheta = 2 * tantheta / (1 + tantheta * tantheta)

        v1 += v2 * tantheta
        v2 += - v1 * sintheta
        v1 += v2 * tantheta

        p.vit[1, ipart] = v1 + 0.5dt * e1
        p.vit[2, ipart] = v2 + 0.5dt * e2

    end

end

export push_x!

function push_x!(p, mesh, dt)

    p.pos .+= p.vit .* dt

    update_cells!(p, mesh)

end
