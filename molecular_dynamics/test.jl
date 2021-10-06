using Plots
gr()
include("types.jl")
include("dynamics.jl")
include("statistics.jl")

# test dynamics of a pair of particles

box = [10, 10, 10]
num_atoms = 2
particles_generator = []
for n=1:num_atoms
    pos = [n*1.03, 0, 0]
    particle = Particle(pos, zeros(3), zeros(3), 1.0, 1.0, 1.0)
    push!(particles_generator, particle)
end
particles = particles_generator
length(particles)

e_pot = round(potential_energy(particles, box), digits=7)
e_kin = round(kinetic_energy(particles), digits=7)
e_tot = e_kin + e_pot
println("initial energy: total $e_tot, potential $e_pot, kinetic $e_kin")
r_list = zeros(0)
v_list = zeros(0)
a_list = zeros(0); t_list = zeros(0)
dt = 0.01
dt_steps = 2000
intermediate = ceil(dt_steps / 100)

for t in 1:dt_steps
    particles = velocity_verlet(particles, dt)
    if t % intermediate == 0
        e_pot = round(potential_energy(particles, box), digits = 7)
        e_kin = round(kinetic_energy(particles), digits = 7)
        e_tot = round(e_kin + e_pot, digits = 7)
        println("energy after $t steps: total $e_tot, potential $e_pot, kinetic $e_kin")
        dr = particles[1].pos - particles[2].pos
        dist = sqrt(sum(dr.*dr))
        dvx = particles[1].vel[1] - particles[2].vel[1]
        dax = particles[1].acc[1] - particles[2].acc[1]
        append!(r_list, dist)
        append!(v_list, dvx)
        append!(a_list, dax)
        append!(t_list, t)
    end
end

plot(t_list, r_list,
    xlabel="t",
    ylabel="dist",
    title="Pair Distance",
    legend = false)
