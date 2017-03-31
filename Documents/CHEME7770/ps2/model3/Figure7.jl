# Script to estimate the acetate and celmass in W3110
# Varma A, Palsson BO (1994) Stoichiometric flux balance models quantitatively predict growth and metabolic by-product
# secretion in wild-type Escherichia coli W3110. Appl Environ Microbiol 60: 3724-31.

# include -
include("include.jl")

# setup the time-scale -
time_start = 0.0
time_stop = 10.0
time_step = 0.1
time_array = collect(time_start:time_step:time_stop)
number_of_timesteps = length(time_array)

# Fire up the max cellmass -
data_dictionary = maximize_cellmass_data_dictionary(time_start,time_stop,time_step)

# Problem specific kinetic parameters -
vmax_glucose_uptake = 10.5    # this was given in the abstract, make sure to cite that
K_glucose_uptake = 18.0 # found on bionumbers

vmax_acetate_uptake = 11.3 #this is also from the Palsson paper!
K_acetate_uptake = 5.0 #found on bionumbers


# initialize the problem -
number_of_external_states = 3
state_array = zeros(number_of_timesteps,number_of_external_states)

# set the ic -
state_array[1,1] = 11.11   # 1 glucose
state_array[1,2] = 0.5     # 2 acetate
state_array[1,3] = 0.005   # 3 cellmass

# capture the exit flags -
exit_flag_array = Int[]

# main loop -
for time_step_index = 1:number_of_timesteps-1

  # make a deepcopy of the data_dictionary -
  copy_data_dictionary = deepcopy(data_dictionary)

  # grab the state -
  glucose = state_array[time_step_index,1]
  acetate = state_array[time_step_index,2]
  cellmass = state_array[time_step_index,3]

  # calculate glucose uptake -
  qGlc = vmax_glucose_uptake*(glucose)/(K_glucose_uptake+glucose)


  @show qGlc

  # setup the species bounds -
  species_bounds_array = copy_data_dictionary["species_bounds_array"]
  species_bounds_array[82,1] = -qGlc
  species_bounds_array[82,2] = -0.99*qGlc

  #getting it to utilize acetate once glucose runs out
  qA = 0.0
  if state_array[time_step_index,1] < 2.0

    # calculate glucose uptake -
    qA = vmax_acetate_uptake*(acetate)/(K_acetate_uptake+acetate)

    @show qA

    species_bounds_array[74,1] = -qA
    species_bounds_array[74,2] = 0.0
  end

  # calculate the fluxes using the LP -
  (objective_value, flux_array, dual_array, uptake_array, exit_flag) = FluxDriver(copy_data_dictionary)

  # grab the fluxes from the flux_array -
  mu = 1.23*flux_array[24]

  qAcetate_production = 0.0
  if state_array[time_step_index,1] < 2.0
    qAcetate_production = -qA
  else
    qAcetate_production = flux_array[35]
  end

  #what is 1.23

  @show mu

  # update the external state -
  state_array[time_step_index+1,1] = glucose -1.0*qGlc*cellmass
  state_array[time_step_index+1,2] = acetate + qAcetate_production*cellmass
  state_array[time_step_index+1,3] = cellmass + (mu)*cellmass

  # correct negatives -
  idx_nz = find(state_array[time_step_index+1,:].<0)
  state_array[time_step_index+1,idx_nz] = 0.0

  # capture the exit flag -
  push!(exit_flag_array,exit_flag)

end

# ok so now we have the growth/production/degradation rates, now use them to find product amount
C = state_array[:,3]
G = state_array[:,1]
A = state_array[:,2]

compare_A_t = [2.0;4.5;7.0;7.5;8.0]
compare_A_v = [0.3;1.0;3.0;3.8;2.5]
compare_G_t = [1.0;3.0;4.0;6.0;7.0;7.5]
compare_G_v = [11.0;10.0;9.5;7.0;3.0;0.0]
compare_C_t = [2.0;4.0;6.0;7.0;8.0]
compare_C_v = [0.01;0.05;0.2;0.5;0.75]

# ok now lets make a graph out of this

clf()

plot(time_array,G,color="orange")
plot(time_array,A,color="blue")
plot(compare_A_t,compare_A_v,color = "blue","--")
plot(compare_G_t,compare_G_v,color = "orange","--")
xlabel("Time (hr)")
ylabel("Acetate (B), Glucose (O), [mM]")
title("Aerobic Batch Culture Analysis")
grid("on")

savefig("/Users/karenna/Documents/CHEME7770/ps2/figs/Figure7met.pdf")

clf()

plot(time_array,C,color="purple")
plot(compare_C_t,compare_C_v,color = "purple","--")
xlabel("Time (hr)")
ylabel("X [g/L]")
title("Aerobic Batch Culture Analysis")
grid("on")

savefig("/Users/karenna/Documents/CHEME7770/ps2/figs/Figure7cell.pdf")
#legend()
#is this making separate plots or no
