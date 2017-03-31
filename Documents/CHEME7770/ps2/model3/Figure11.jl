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
data_dictionary = maximize_anaerobic_cellmass_data_dictionary(time_start,time_stop,time_step)

# Problem specific kinetic parameters -
vmax_glucose_uptake = 18.5    # this was given in the abstract, make sure to cite that
K_glucose_uptake = 25.0

#this is not a real value I just made it up
#OK I need to determine this with Bionumbers
#mess with this

# initialize the problem -
number_of_external_states = 5
state_array = zeros(number_of_timesteps,number_of_external_states)

# set the ic -
state_array[1,1] = 11.0   # 1 glucose
state_array[1,2] = 0.001    # 2 acetate
state_array[1,3] = 0.01   # 3 cellmass
state_array[1,4] = 0.001   # 4 formate
state_array[1,5] = 0.001   # 5 ethanol


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
  formate = state_array[time_step_index,4]
  ethanol = state_array[time_step_index,5]

  # calculate glucose uptake -
  qGlc = vmax_glucose_uptake*(glucose)/(K_glucose_uptake+glucose)

  @show qGlc

  # setup the species bounds -
  species_bounds_array = copy_data_dictionary["species_bounds_array"]
  species_bounds_array[82,1] = -qGlc
  species_bounds_array[82,2] = -0.99*qGlc


  # calculate the fluxes using the LP -
  (objective_value, flux_array, dual_array, uptake_array, exit_flag) = FluxDriver(copy_data_dictionary)

  # grab the fluxes from the flux_array -
  mu = 1.23*flux_array[24]
  qAcetate_production = flux_array[35]
  qFormate_production = flux_array[41]
  qEthanol_production = flux_array[40]

  #what is 1.23

  @show mu

  # update the external state -
  state_array[time_step_index+1,1] = glucose -1.0*qGlc*cellmass
  state_array[time_step_index+1,2] = acetate + qAcetate_production*cellmass
  state_array[time_step_index+1,3] = cellmass + mu*cellmass
  state_array[time_step_index+1,4] = formate + qFormate_production*cellmass
  state_array[time_step_index+1,5] = ethanol + qEthanol_production*cellmass

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
F = state_array[:,4]
E = state_array[:,5]

compare_A_t = [2.0;4.0;6.0;8.0;9.0]
compare_A_v = [0.0;0.75;1.5;4.5;7.5]
compare_G_t = [2.0;4.0;6.0;8.0;9.0]
compare_G_v = [11.0;10.0;8.5;5.0;3.0]
compare_C_t = [2.0;4.0;6.0;8.0;9.0]
compare_C_v = [0.01;0.025;0.05;0.125;0.25]
compare_F_t = [2.0;4.0;6.0;8.0;9.0]
compare_F_v = [0.0;1.0;4.0;10.0;14.0]
compare_E_t = [2.0;4.0;6.0;8.0;9.0]
compare_E_v = [0.0;0.5;1.5;3.5;7.0]


# ok now lets make a graph out of this
clf()

plot(time_array,G,color="orange")
plot(compare_G_t,compare_G_v,color = "orange", "--")
plot(time_array,A,color="blue")
plot(compare_A_t,compare_A_v,color = "blue", "--")
plot(time_array,F,color="green")
plot(compare_F_t,compare_F_v,color = "green", "--")
plot(time_array,E,color="pink")
plot(compare_E_t,compare_E_v,color = "pink", "--")
xlabel("Time (hr)")
ylabel("Acetate (B), Glucose (O), Formate (G), Ethanol (P) [mM]")
title("Anaerobic Batch Culture Analysis")
grid("on")

savefig("/Users/karenna/Documents/CHEME7770/ps2/figs/Figure11met.pdf")

clf()

plot(time_array,C,color="purple")
plot(compare_C_t,compare_C_v,color = "purple", "--")
xlabel("Time (hr)")
ylabel("X [g/L]")
title("Anaerobic Batch Culture Analysis")
grid("on")


savefig("/Users/karenna/Documents/CHEME7770/ps2/figs/Figure11cell.pdf")
