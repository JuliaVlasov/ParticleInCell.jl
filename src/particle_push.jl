export push_v!

function push_v!(p, dt)

    for ipart = 1:p.nbpart

        e1 = p.epx[ipart]
        e2 = p.epy[ipart]

        p.vit[1, ipart] += 0.5dt * e1
        p.vit[2, ipart] += 0.5dt * e2

        tantheta = 0.5dt * p.bpz[ipart]
        sintheta = 2 * tantheta / (1 + tantheta * tantheta)

        p.vit[1, ipart] += p.vit[2, ipart] * tantheta
        p.vit[2, ipart] += - p.vit[1, ipart] * sintheta
        p.vit[1, ipart] += p.vit[2, ipart] * tantheta

        p.vit[1, ipart] += 0.5dt * e1
        p.vit[2, ipart] += 0.5dt * e2

    end

end

export push_x!

function push_x!(p, mesh, dt)

    p.pos .+= p.vit .* dt

    update_cells!(p, mesh)

end
