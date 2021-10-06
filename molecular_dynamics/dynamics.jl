using Distributions

function force(particles, box)
    num_atoms = length(particles)
    f_partial = zeros(num_atoms, num_atoms, 3)

    for m in 1:num_atoms, n in (m+1):num_atoms
        dr = particles[m].pos - particles[n].pos
        dr += - round.(dr./box, digits=0) .* box # PBC
        r2 = sum(dr.*dr) / particles[n].σ^2
        r4 = r2 * r2
        r8 = r4 * r4
        r14 = r8 * r4 * r2
        f_partial[m, n, :] = - particles[n].ϵ * (6.0/r8 - 12/r14) * dr
        f_partial[n, m, :] = - f_partial[m, n, :]
    end

    f_total = sum(f_partial, dims=2)
    f_total = dropdims(f_total, dims=2)
    return f_total
end

function velocity_verlet(particles, dt)
    num_atoms = length(particles)
    dt2 = dt/2
    dtsq2 = dt^2 / 2
    for n in 1:num_atoms
        particles[n].pos += dt*particles[n].vel + dtsq2*particles[n].acc
        particles[n].vel += dt2*particles[n].acc
    end
    force_update = force(particles, box)
    for n in 1:num_atoms
        particles[n].acc = force_update[n, :] / particles[n].μ
        particles[n].vel += dt2*particles[n].acc
    end
    return particles
end

function rescale_velocities(particles, temperature)
    num_atoms = length(particles)
    k_B = 1  # in units of epsilon, temperature in units of 125
    for n in 1:num_atoms
        boltzmann_factor = sqrt(particles[n].μ / (2 * k_B * temperature))
        particles[n].vel = rand(Normal(0, 1), 3) * boltzmann_factor
    end
    return particles
end

function rebox_particles(particles, box)
    num_atoms = length(particles)
    for n in 1:num_atoms
        dr = particles[n].pos
        dr += - floor.(dr ./ box) .* box
        particles[n].pos = dr
    end
    return particles
end

function distance(pos1, pos2, box)
    dr = pos2 -  pos1
    dr += - floor.(dr ./ box) .* box
    dist = sqrt(sum(dr .* dr))
end
