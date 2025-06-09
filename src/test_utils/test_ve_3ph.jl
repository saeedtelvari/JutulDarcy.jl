using JutulDarcy, Jutul
using Test

G = CartesianMesh((2, 2), (2.0, 2.0))
G = reservoir_domain(G)


nc = number_of_cells(G)
pv = copy(pore_volume(G))
timesteps = [1.0, 2.0]*3600*24


context = DefaultContext()

@assert isa(context, JutulContext)
parameters = nothing




inj = 1
prod = nc
pv[inj] *= 1000
pv[prod] *= 1000

bar = 1e5
p0 = 100*bar # 100 bar

mu = 1e-3    # 1 cP
cl = 1e-5/bar
pRef = 100*bar
rhoLS = 1000
A = AqueousPhase()
L = LiquidPhase()
V = VaporPhase()
sys = ImmiscibleSystem([A, L, V])
model = SimulationModel(G, sys, context = context)

kr = BrooksCoreyRelativePermeabilities(sys, [2, 2, 2])
s = model.secondary_variables
s[:RelativePermeabilities] = kr
s[:PhaseMassDensities] = ConstantCompressibilityDensities(sys, pRef, rhoLS, cl)

tot_time = sum(timesteps)
forces = setup_forces(model)

p_init = repeat([p0], nc)
p_init[inj] = 2*p0
p_init[prod] = p0/2

s_inj = [0.5, 0.3, 0.2]
s_res = [0.3, 0.3, 0.4]

s_inj = [0.1, 0.0, 0.9]
s_res = [0.9, 0.0, 0.1]

s_inj = [0.1, 0.9, 0.0]
s_res = [0.9, 0.1, 0.0]

s_init = repeat(s_inj, 1, nc)
s_init[:, inj] .= s_inj

# State is dict with pressure in each cell
init = Dict(:Pressure => p_init, :Saturations => s_init)


linear_solver = :auto
linear_solver = reservoir_linsolve(model, linear_solver = linear_solver)

if isnothing(parameters)
    parameters = setup_parameters(model)
end
state0 = setup_state(model, init)

sim = Simulator(model, state0 = state0, parameters = parameters)

arg = NamedTuple()
arg = (linear_solver = linear_solver, )
cfg = simulator_config(sim; info_level = -1, debug_level = 1, arg...)
simulate(sim, timesteps, forces = forces, config = cfg)


 # simulate_reservoir(state0, model, dt)









