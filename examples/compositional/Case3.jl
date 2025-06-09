# # A more complex compositional model
# This example sets up a more complex compositional simulation with five
# different components. Other than that, the example is similar to the others
# that include wells and is therefore not commented in great detail.
using MultiComponentFlash

h2o = MolecularProperty(0.018015268, 22.064e6, 647.096, 5.595e-05, 0.3442920843)
co2 = MolecularProperty(0.0440098, 7.3773e6, 304.1282, 9.412e-05, 0.22394)
c1 = MolecularProperty(0.0160, 4.60e6, 190.6, 9.863e-5, 0.011)

bic = zeros(3, 3)
bic[2,1] = 1 # CO2-H2O interaction
bic[1,2] = 1 # H2O-CO2 interaction
bic[1,3] = 0.5 # CO2-CO2 interaction
bic[3,1] = 0.5 # CH4-H2O interaction
mixture = MultiComponentMixture([h2o, co2, c1], A_ij = bic, names = ["H2O", "CO2", "CH4"])
eos = GenericCubicEOS(mixture, PengRobinson())

using Jutul, JutulDarcy, GLMakie
Darcy, bar, kg, meter, Kelvin, day = si_units(:darcy, :bar, :kilogram, :meter, :Kelvin, :day)
nx = 20
ny = 1
nz = 10

dims = (nx, ny, nz)
g = CartesianMesh(dims, (4000.0, 1.0, 100.0))
nc = number_of_cells(g)
K = repeat([0.5*Darcy], 1, nc)
res = reservoir_domain(g, porosity = 0.25, permeability = K, temperature = 333.15*Kelvin)
# Set up a vertical well in the first corner, perforated in all layers
prod = setup_vertical_well(g, K, nx, ny, name = :Producer)
# Set up an injector in the opposite corner, perforated in all layers
inj = setup_vertical_well(g, K, 1, 1, name = :Injector)

rhoLS = 1000.0*kg/meter^3
rhoVS = 100.0*kg/meter^3

rhoS = [rhoLS, rhoVS]
L, V = LiquidPhase(), VaporPhase()
# Define system and realize on grid
sys = MultiPhaseCompositionalSystemLV(eos, (L, V))
model, parameters = setup_reservoir_model(res, sys, wells = [inj, prod], block_backend = true);
kr = BrooksCoreyRelativePermeabilities(sys, 2.0, 0.0, 1.0)
model = replace_variables!(model, RelativePermeabilities = kr)

push!(model[:Reservoir].output_variables, :Saturations)

# Zc = zeros(nc)'
Zc = ones(nc)'
Zm = ones(nc)'

# println("Number of cells: ", Zc)
Zc[1:div(nc,2)] .= 0.0 # Set first half cells to be CO2
# println("Number of cells: ", Zm)
Zm[div(nc,2)+1:nc] .= 0.0 # Set second half cells to be CH4


# println("Number of cells: ", Zc)
Zw = zeros(nc)'
Z = vcat(Zw, Zc, Zm)
P = ones(nc).*100*bar
state0 = setup_reservoir_state(model, Pressure =P, OverallMoleFractions = Z);

dt = repeat([5.0]*day, 365)
rate_target = TotalRateTarget(0.00)
I_ctrl = InjectorControl(rate_target, [0, 1, 0], density = rhoVS)
bhp_target = BottomHolePressureTarget(100*bar)
P_ctrl = ProducerControl(bhp_target)

controls = Dict()
controls[:Injector] = I_ctrl
controls[:Producer] = P_ctrl
forces = setup_reservoir_forces(model, control = controls)
ws, states = simulate_reservoir(state0, model, dt, parameters = parameters, forces = forces);
# ## Once the simulation is done, we can plot the states
# ### CO2 mole fraction
sg = states[end][:OverallMoleFractions][2, :]
fig, ax, p = plot_cell_data(g, sg)
fig
# ### Gas saturation
sg = states[end][:Saturations][2, :]
fig, ax, p = plot_cell_data(g, sg)
fig
# ### Pressure
p = states[end][:Pressure]
fig, ax, p = plot_cell_data(g, p)
fig

plot_reservoir(model, states, step = length(dt))

# Create a new figure for the line plots
f = GLMakie.Figure()
ax = GLMakie.Axes(f[1, 1], xlabel = "Cell index", ylabel = "CO2 mass")

# Plot all lines in the same axes
for i in 1:10
       co2_reshaped = reshape(states[i][:TotalMasses][2,:], (nx, nz))
       GLMakie.lines!(f, co2_reshaped[:,25], label = "Step $i")
end

# Add legend
# GLMakie.axislegend(ax)

# Display the figure
display(f)
