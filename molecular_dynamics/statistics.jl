function potential_energy(particles, box)
    num_atoms = length(particles)
    e_pot = 0
    for m in 1:num_atoms, n in (m+1):num_atoms
        pos_m = particles[m].pos
        pos_n = particles[n].pos
        dr = pos_m - pos_n - round.((pos_m - pos_n)./ box, digits=0) .* box # PBC
        r2 = sum(dr.*dr) / particles[n].σ^2
        r6 = r2^3
        r12 = r6^2
        e_pair = particles[n].ϵ *(1.0/r12 - 1.0/r6)
        e_pot += e_pair
    end
    return e_pot
end

function kinetic_energy(particles)::Float64
    e_kin = 0.0
    for n in 1:length(particles)
        vsq = sum(particles[n].vel .* particles[n].vel)
        e_kin += 0.5 * particles[n].μ * vsq
    end
    return e_kin
end

function pair_correlation(particles)
    r_step = particles[1].σ /20
    r_min = 0.75*particles[1].σ
    r_max = 3*particles[1].σ
    r_border = r_min : r_step : r_max
    bins = zeros(length(r_border))
    for m in 1:num_atoms, n in (m+1):num_atoms
        pos_m = particles[m].pos
        pos_n = particles[n].pos
        dr = abs(pos_m - pos_n)
        dr += - round.((pos_m - pos_n)./ box, digits=0)
        dist = sqrt(sum(dr.*dr))
        nb = max(1, Int16(round((dist - r_min) / r_step, digits=0)))
        if nb <= length(bins)
            bins[nb] += 1
        end
    end
    n = 0
    g_r = zeros(length(bins))
    for r in r_border
        n += 1
        g_r[n] = bins[n] / r^2 / (4 * pi) 
    end
    return g_r, r_border
end
