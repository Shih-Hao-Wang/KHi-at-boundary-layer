using Oceananigans, Printf, JLD2, Random

const Re₀ = 50.0
const Ri₀ = 1/80
const Pr  = 1.0
const tᵣ  = 100.0
const R₀  = 1.0
const h   = 3.0
const Lx  = 64.0
const Lz  = 48.0
const Nx, Nz = 768, 576
const Cd  = 2e-3

wallbc = :noslip                            
fname = "khiwall2d_Re50_h3_noslip_v2fine"   

arch = CPU()
grid = RectilinearGrid(arch; size=(Nx, Nz), x=(0, Lx), z=(0, Lz),
                       topology=(Periodic, Flat, Bounded), halo=(4, 4))

U★(z) = tanh(z - h)
B★(z) = Ri₀ * tanh(R₀ * (z - h))
relax_u(x, z, t, u) = (U★(z) - u) / tᵣ
relax_b(x, z, t, b) = (B★(z) - b) / tᵣ
Fu = Forcing(relax_u, field_dependencies = :u)
Fb = Forcing(relax_b, field_dependencies = :b)

if wallbc == :noslip
    u_bcs = FieldBoundaryConditions(bottom = ValueBoundaryCondition(0.0))
else
    drag_u(x, t, u) = -Cd * abs(u) * u
    u_bcs = FieldBoundaryConditions(bottom = FluxBoundaryCondition(drag_u, field_dependencies = :u))
end

model = NonhydrostaticModel(grid;
    timestepper = :RungeKutta3,
    advection   = Centered(order = 4),     
    closure     = ScalarDiffusivity(ν = 1/Re₀, κ = 1/(Re₀*Pr)),
    buoyancy    = BuoyancyTracer(),
    tracers     = :b,
    forcing     = (; u = Fu, b = Fb),
    boundary_conditions = (; u = u_bcs))

Random.seed!(1234)                          
u₀(x, z) = U★(z)
b₀(x, z) = B★(z) + Ri₀ * 1e-3 * (2rand() - 1)
set!(model, u = u₀, b = b₀)

simulation = Simulation(model; Δt = 0.02, stop_time = 800.0)
wizard = TimeStepWizard(cfl = 0.5, max_change = 1.1, max_Δt = 0.1)  
simulation.callbacks[:wizard] = Callback(wizard, IterationInterval(10))

u, v, w = model.velocities
b  = model.tracers.b
Ū  = Field(Average(u, dims = 1))
B̄  = Field(Average(b, dims = 1))
ζ  = Field(∂z(u) - ∂x(w))

progress(sim) = @printf("t=%7.1f  Δt=%.3f  max|w|=%.3e\n",
                        time(sim), sim.Δt, maximum(abs, w))
simulation.callbacks[:progress] = Callback(progress, TimeInterval(10.0))

simulation.output_writers[:slices] = JLD2Writer(model, (; b, ζ, u, w);   # ← v2:加存 w
    filename = fname * "_slices.jld2",
    schedule = TimeInterval(5.0), overwrite_existing = true)
simulation.output_writers[:profiles] = JLD2Writer(model, (; Ū, B̄);
    filename = fname * "_profiles.jld2",
    schedule = TimeInterval(2.0), overwrite_existing = true)
simulation.output_writers[:checkpoint] = Checkpointer(model;
    schedule = TimeInterval(500.0), prefix = fname * "_chk",
    overwrite_existing = true)

run!(simulation)
@info "done → $(fname)"
