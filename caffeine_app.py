import matplotlib.pyplot as plt
import numpy as np
import streamlit as st

# ... your other imports and functions ...

def plot_caffeine_over_time(mass, caffeine_amount, nM_adenosine, gender, adenosine_Kd_A1, adenosine_Kd_A2A, caffeine_Kd_A1, caffeine_Kd_A2A):
    # Define constants
    half_life = 4.0  # half-life of caffeine in hours
    time = np.linspace(0, 24, 300)  # 24 hours, sampled at 300 points
    caffeine = np.zeros_like(time)

    # For the first hour, linearly increase caffeine to its max
    absorption_mask = time < 1
    caffeine[absorption_mask] = caffeine_amount * time[absorption_mask]

    # After absorption, decrease based on half-life
    decay_mask = time >= 1
    decay_times = time[decay_mask] - 1  # Subtract the absorption hour
    caffeine[decay_mask] = caffeine_amount * (0.5 ** (decay_times / half_life))

    # Plot
    # Create a figure with two y-axes
    fig, ax1 = plt.subplots()

    # Plot the caffeine levels on the first y-axis
    ax1.plot(time, caffeine, label='Caffeine (mg)', color='tab:blue')
    ax1.set_xlabel('Time (s)')
    ax1.set_ylabel('Caffeine (mg)', color='tab:blue')
    ax1.tick_params(axis='y', labelcolor='tab:blue')

    # Create a second y-axis
    ax2 = ax1.twinx()

    # Calculate fC_simultaneous_A1 based on current inputs
    INPUTS = [[mass, caff, nM_adenosine, gender, adenosine_Kd_A1, adenosine_Kd_A2A, caffeine_Kd_A1, caffeine_Kd_A2A] for caff in caffeine]
    fC_simultaneous_A1 = [receptor_occupancy(*L)["A1_occupancy"][0] for L in INPUTS]

    # Plot fC_simultaneous_A1 on the second y-axis
    ax2.plot(time, fC_simultaneous_A1, label='fC_simultaneous_A1', color='tab:red')
    ax2.set_ylabel('fC_simultaneous_A1', color='tab:red')
    ax2.tick_params(axis='y', labelcolor='tab:red')

    # Calculate fC_simultaneous_A1 based on current inputs
    INPUTS = [[mass, caff, nM_adenosine, gender, adenosine_Kd_A1, adenosine_Kd_A2A, caffeine_Kd_A1, caffeine_Kd_A2A] for caff in caffeine]
    fC_simultaneous_A2 = [receptor_occupancy(*L)["A2A_occupancy"][0] for L in INPUTS]
    # Plot fC_simultaneous_A1 on the second y-axis
    ax2.plot(time, fC_simultaneous_A2, label='fC_simultaneous_A1', color='tab:orange')
    ax2.set_ylabel('fC_simultaneous_A2', color='tab:red')

    # Show legend
    ax1.legend(loc='upper left')
    ax2.legend(loc='upper right')

    # Display the plot in Streamlit
    st.pyplot(fig)


def receptor_occupancy(mass, mg_caffeine, nM_adenosine, gender, adenosine_Kd_A1, adenosine_Kd_A2A, caffeine_Kd_A1, caffeine_Kd_A2A):
    # toggle:
    # man vs woman for vd_per_kg
    # using 50 mL for blood in brain volume but the blood volume in the brain would be around 42-56 ml.
    # todo:
    # fix label overlap with slider
    # give other slider inputs to receptor_occupancy

    # Given constants
    if gender == 'Male':
        vd_per_kg = .733  # MEN # Volume of distribution per kg in mL/kg
    elif gender == 'Female':
        vd_per_kg = .591  # WOMEN # Volume of distribution per kg in mL/kg

    molecular_weight_adenosine = 267.2413  # g/mol
    molecular_weight_caffeine = 194.19  # g/mol
    
    fix = 10 # applied one order of magnitude change to Kd to fix it
    Kd_adenosine_A1 = adenosine_Kd_A1 * 1e-3 * fix  # 900 μM in Molar for A1 receptor # .009 works
    Kd_caffeine_A1 = caffeine_Kd_A1  # 40 μM in Molar for A1 receptor
    Kd_adenosine_A2A = adenosine_Kd_A2A * 1e-3  # Average value for A2A receptor
    Kd_caffeine_A2A = caffeine_Kd_A2A  # Average value for A2A receptor
    
    # Calculate volume of distribution for the given mass
    vd_total = vd_per_kg * mass # (mL)
    
    # Calculate caffeine concentration in the "theoretical distribution volume"
    caffeine_density = mg_caffeine / vd_total # (mg / mL) 
    caffeine_mass = caffeine_density * .001 * 1.120 # (g) # the blood volume in the brain would be around 42-56 ml.
    ## ## ## ##
    # caffeine_mass = .4275
    ## ## ## ##


    # Convert nM adenosine to M and calculate its concentration in the brain (nM to g)
    adenosine_mass = nM_adenosine * 1e-9 * molecular_weight_adenosine
    
    # Calculate the mass ratio of caffeine to adenosine in the brain
    ratio = caffeine_mass / adenosine_mass
    
    μM_caffeine = caffeine_mass / molecular_weight_caffeine * 1e6 #(μM)
    μM_adenosine = nM_adenosine / 1000
    # Calculate the individual probabilities of binding to A1 receptors
    f_caffeine_A1 = μM_caffeine / (μM_caffeine + Kd_caffeine_A1)
    f_adenosine_A1 = μM_adenosine / (μM_adenosine + Kd_adenosine_A1)
    
    # Calculate the individual probabilities of binding to A2A receptors
    f_caffeine_A2A = μM_caffeine / (μM_caffeine + Kd_caffeine_A2A)
    f_adenosine_A2A = μM_adenosine / (μM_adenosine + Kd_adenosine_A2A)
    
    # Calculate the percentage of A1 receptors occupied by each molecule in the presence of competition
    fA_simultaneous_A1 = μM_adenosine / (μM_adenosine + Kd_adenosine_A1 * (1 + μM_caffeine/Kd_caffeine_A1))
    fC_simultaneous_A1 = μM_caffeine / (μM_caffeine + Kd_caffeine_A1 * (1 + μM_adenosine/Kd_adenosine_A1))
    
    # Calculate the percentage of A2A receptors occupied by each molecule in the presence of competition
    fA_simultaneous_A2A = μM_adenosine / (μM_adenosine + Kd_adenosine_A2A * (1 + μM_caffeine/Kd_caffeine_A2A))
    fC_simultaneous_A2A = μM_caffeine / (μM_caffeine + Kd_caffeine_A2A * (1 + μM_adenosine/Kd_adenosine_A2A))
    
    return {
        "caffeine_density": caffeine_density,
        "caffeine_mass": caffeine_mass,
        "μM_adenosine": μM_adenosine,
        "adenosine_mass": adenosine_mass,
        "ratio_caffeine_to_adenosine": ratio,
        "binding_probability_A1": (f_caffeine_A1 * 100 , f_adenosine_A1 * 100),
        "binding_probability_A2A": (f_caffeine_A2A * 100, f_adenosine_A2A * 100),
        "A1_occupancy": (fC_simultaneous_A1 * 100, fA_simultaneous_A1 * 100),
        "A2A_occupancy": (fC_simultaneous_A2A * 100, fA_simultaneous_A2A * 100)
    }




# Main function for Streamlit app
def main():
    st.title("Caffeine Receptor Occupancy Calculator")

    # Sliders and selectors
    mass = st.slider("Body Mass (kg):", 30, 120, 80)
    mg_caffeine = st.slider("Caffeine (mg):", 0, 1000, 800)
    nM_adenosine = st.slider("Adenosine (nM):", 0, 500, 30)
    gender = st.radio("Gender:", ['Male', 'Female'])

    # Advanced settings (collapsed by default)
    with st.expander("Advanced Settings"):
        adenosine_Kd_A1 = st.slider("Adenosine Kd A1 (nM):", 0.3, 1.5, 0.9)
        adenosine_Kd_A2A = st.slider("Adenosine Kd A2A (nM):", 10.0, 50.0, 30.0)
        caffeine_Kd_A1 = st.slider("Caffeine Kd A1 (μM):", 20.0, 60.0, 40.0)
        caffeine_Kd_A2A = st.slider("Caffeine Kd A2A (μM):", 10.0, 50.0, 45.0)

    # Compute results
    result = receptor_occupancy(mass, mg_caffeine, nM_adenosine, gender, adenosine_Kd_A1, adenosine_Kd_A2A, caffeine_Kd_A1, caffeine_Kd_A2A)
    
    # Display results
    st.subheader("Results")
    for key, val in result.items():
        st.write(f"{key}: {val}")

    # Plot the caffeine levels
    plot_caffeine_over_time(mass, mg_caffeine, nM_adenosine, gender, adenosine_Kd_A1, adenosine_Kd_A2A, caffeine_Kd_A1, caffeine_Kd_A2A)

# Run the app
if __name__ == "__main__":
    main()