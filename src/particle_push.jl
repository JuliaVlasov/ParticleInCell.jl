export push_v!

function push_v!(p, dt)

    for ipart = 1:p.nbpart

        e1 = p.epx[ipart]
        e2 = p.epy[ipart]

        p.vit[ipart, 1] += 0.5dt * e1
        p.vit[ipart, 2] += 0.5dt * e2

        tantheta = 0.5dt * p.bpz[ipart]
        sintheta = 2 * tantheta / (1 + tantheta * tantheta)

        p.vit[ipart, 1] += p.vit[ipart,2] * tantheta
        p.vit[ipart, 2] += - p.vit[ipart,1] * sintheta
        p.vit[ipart, 1] += p.vit[ipart,2] * tantheta

        p.vit[ipart, 1] += 0.5dt * e1
        p.vit[ipart, 2] += 0.5dt * e2

    end

end

export push_x!

function push_x!(p, mesh, dt)

    p.pos .+= p.vit .* dt
    update_cells!(p, mesh)

end
