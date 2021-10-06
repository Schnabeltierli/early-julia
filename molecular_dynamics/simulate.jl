using Plots
gr()
include("types.jl")
include("dynamics.jl")
include("statistics.jl")

# LJ parameters for liquid Argon: ϵ/k_B = 125.7, σ=0.3345nm
# Triple point	83.8058 K, 69 kPa
# Critical point 150.87 K, 4.898 MPa
# Melting point	83.80 K
# Boiling point	87.30 K

box = [5, 5, 5]
num_atoms = 125
particles_generator = []
for n=1:num_atoms
    nx = (n-1) % 5 + 0.5
    ny = round(floor((n-1) / 5), digits=0) + 0.5
    nz = round(floor((n-1) / 25), digits=0) + 0.5
    pos = [nx, ny, nz]
    particle = Particle(pos, zeros(3), zeros(3), 1.0, 1.0, 1.0)
    push!(particles_generator, particle)
end
particles = particles_generator
length(particles)

e_pot = round(potential_energy(particles, box), digits=7)
e_kin = round(kinetic_energy(particles), digits=7)
e_tot = e_kin + e_pot
println("initial energy: total $e_tot, potential $e_pot, kinetic $e_kin")
ep_list = zeros(0); ek_list = zeros(0); t_list = zeros(0)
dt = 0.02
temperature = 1.0  # 1.0 is approx. T = 125 K
particles = rescale_velocities(particles, temperature)
dt_steps = 10000
num_collect = 50 # number of data collection steps
intermediate = ceil(dt_steps / num_collect)
g_r0, r_border = pair_correlation(particles)
g_r = zeros(length(g_r0))

for t in 1:dt_steps
    particles = velocity_verlet(particles, dt)
    if t % intermediate == 0
        e_pot = round(potential_energy(particles, box), digits = 7)
        e_kin = round(kinetic_energy(particles), digits = 7)
        e_tot = round(e_kin + e_pot, digits = 7)
        println("energy after $t steps: total $e_tot, potential $e_pot, kinetic $e_kin")
        append!(ek_list, e_kin)
        append!(ep_list, e_pot)
        append!(t_list, t)
        g_r0, r_border = pair_correlation(particles)
        for n in 1:length(r_border)
            g_r[n] += g_r0[n]
        end
        particles = rescale_velocities(particles, temperature)
        particles = rebox_particles(particles, box)
    end
end

# energy
plot(t_list, ek_list,
    xlabel="t",
    ylabel="energy",
    title="Energy",
    label="kinetic energy",
    legend = true)
plot!(t_list, ep_list, label="potential energy")
# pair correlation
g_r1 = g_r / num_collect
plot(r_border, g_r1,
    xlabel="r",
    ylabel="g(r)",
    title="Pair Correlation",
    legend= false)
# simulation box
x = zeros(0); y = zeros(0); z = zeros(0)
particles = rebox_particles(particles, box)
for n in 1:num_atoms
    append!(x, particles[n].pos[1])
    append!(y, particles[n].pos[2])
    append!(z, particles[n].pos[3])
end
scatter(x, y, z,
    title = "Simulation Box",
    legend=false)
