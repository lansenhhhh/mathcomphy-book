"""

"""
import numpy as np
import matplotlib.pyplot as pyplot
from particle3D import Particle3D
import sys

def read_data():
    '''
    create particle3d objects from input file.
    
    return a list of particle3d objects
    '''
    
    filename = sys.argv[1]
    n = []
    particles=[]
    with open(filename, 'r') as file:
        for line in file:
            n.append(line)
    for i in range(len(n)):
         particles.append(Particle3D.read_line(n[i]))
   
    return particles


def compute_separations(particles):
    separations = np.zeros((len(particles),len(particles), 3)) 
    '''
    reshape the 3D array in the form of [n,n,3]
    n is number of particles
    particles are in form [p1, p2, p3.......]
    so len(particles) is equal to n
    '''
    
    for i in range (len(particles)):
        for j in range (i+1, len(particles)):
            '''
            the aim of choosing range of j from i to n is to avoid overcounting           
            '''
            
            separation_x = particles[i].position[0] - particles[j].position[0]
            separation_y = particles[i].position[1] - particles[j].position[1]
            separation_z = particles[i].position[2] - particles[j].position[2]
            
            '''
            particle[i] means the ith component of particle
            for example, position[0] is the x-compoonent
            '''
            
            separations[i,j] = [separation_x, separation_y,separation_z]
            separations[j,i] = [-1*separation_x, -1*separation_y,-1*separation_z]
            
            '''
            the vector from a to b is the negative of the vector from b to a.
            
            '''
            
    return separations
            
                        
     
def compute_forces_potential(particles, separations):
    forces = np.zeros((len(particles), 3))
    potential = 0.0

    for i in range(len(particles)):
        for j in range(len(particles)):
            if (j > i) and particles[j].label != "Mercury":
                separation = np.linalg.norm(separations[i, j])
                Force = (8.887692593e-10 * particles[j].mass * particles[i].mass * separations[i,j])/(separation)**3
                #work out the force between the two particles
                forces[j] += Force
                forces[i] -= Force
                
                potential -= (8.887692593e-10 * particles[j].mass * particles[i].mass) / (separation)    
    return forces, potential


def compute_energy_deviation(energy, initial_energy): # compute energy deviation
    E_min=np.min(energy)
    E_max=np.max(energy)
    deviation = abs(E_min-E_max)/initial_energy
    return deviation

def compute_distance(p1,p2): 
    return np.linalg.norm(p1.position - p2.position)
    

def main():
    # notice when gives incorrect inputs.
    if len(sys.argv)!=5:
        sys.exit("Please input at console: %run Simulation.py <input_file> <output_file> <number_of_timesteps> <timestep_size>")
    try:
        dt = float(sys.argv[4])
        nstep = int(sys.argv[3])
    except ValueError:
        sys.exit("Input ValueError: Please input at console: %run Simulation.py <input_file> <output_file> <number_of_timesteps> <timestep_size>")
        

    # Initial conditions of the system
    particles = read_data()
    time = 0.0

    
    # reconstruct particles into solar objects in order, spectialize Sun, Moon, and Earth
    solar_particles = [None, None, None]
    others = []
    for particle in particles:
        if particle.label == "Sun":
            solar_particles[0] = particle
        elif particle.label == "Moon":
            solar_particles[1] = particle
        elif particle.label == "Earth":
            solar_particles[2] = particle
        else:
            others.append(particle)
    solar_particles.extend(others)
    
        
    # subtract the centre-of-mass velocity
    com_velocity = Particle3D.com_velocity(solar_particles)
    for i in range(len(solar_particles)):
        solar_particles[i].velocity -= com_velocity
        
        
    # Initialise arrays that we will store results in
    n = len(solar_particles)
    times = np.zeros(nstep)
    energys = np.zeros(nstep)
    energy_0 = 0
    positions = np.zeros((n, nstep, 3))
    forces = np.zeros((n, nstep, 3)) #added
    solar_particles_velocity_x = np.zeros((n,nstep))
    
    # store initial velocities for all particles for orbit period measurement
    init_velocities = np.zeros((n,3))
    for j in range(n):
        init_velocities[j] = solar_particles[j].velocity


    # compute initial forces for first loop iteration
    separations = compute_separations(solar_particles)
    forces, _ = compute_forces_potential(solar_particles, separations)
    
    # apsides
    distances_to_sun = np.zeros((n,nstep)) #initialize arrays for 
    distances_to_earth = np.zeros((nstep))
    
    
    # Main time integration loop
    lines = [] #empty list for output file
    lines_2 = []
    
    for i in range(nstep):
        times[i] = time
        time += dt

        # update all particle positions     
        for j in range(n):
            solar_particles[j].update_position_2nd(dt, forces[j])
            
        # store particle positions in array
            if j != 1:
                positions[j, i] = solar_particles[j].position
            
            # planets distances to sun
                distances_to_sun[j, i] = compute_distance(solar_particles[0], solar_particles[j])
            
            
        # moon distance to sun
        distances_to_earth[i] = compute_distance(solar_particles[1], solar_particles[2])
        
        # get new separations and new forces on all particles, and the potential
        separations_new = compute_separations(solar_particles)
        
        force_new, potential = compute_forces_potential(solar_particles, separations_new)
            
        
        # update all particle velocities
        for j in range(n):
            solar_particles[j].update_velocity(dt, (forces[j] + force_new[j]) /2)
            
            solar_particles_velocity_x[j,i] = solar_particles[j].velocity[0]
            
            #check for complete period by checking x_compnent_velocity's value with significance of 5%
            
        # replace forces with new forces for next iteration
        forces = force_new
        # compute the kinetic energy and save the total energy
        KE = Particle3D.total_kinetic_energy(solar_particles)
        
        Total_energy = KE + potential
        energys[i] = Total_energy
        if i == 0:
            energy_0 = Total_energy
        
        # write output file with correct form in each timestep.
        lines.append(n)
        lines.append("point = " + str(i+1))
        for j in range(n):
            l = str(solar_particles[j])
            lines.append(l)
            
    # find energy deviation
    deviation = compute_energy_deviation(energys, energy_0)
    
    # some initilizations
    time_a = np.zeros(n)
    time_b = np.zeros(n)
    periods = np.zeros(n)
    finish_orbit_index = np.zeros(n)
    half_periods = np.zeros(n)
    #loop except sun and moon
    for j in range(2,n):
        for i in range(nstep-1): 
        #if particle's x_velocity changes signs at a timestep means it passes through the y-axis of its orbit at this timestep.
            if (solar_particles_velocity_x[j,i] * solar_particles_velocity_x[j,i+1]) < 0:
                time_b[j] = time_a[j]
                time_a[j] = i
                half_periods[j] += (time_a[j] - time_b[j]) * (dt) # time spends in the "half period"
                finish_orbit_index[j] += 1 
                #index determine how many times it goes across the y-axis, happens two times in a period.
            else:
                pass
        if finish_orbit_index[j] != 0: #avoid dividing by 0
            periods[j] =  2 * half_periods[j]/finish_orbit_index[j]
        
        
    # find planets sun apsides
    max_PS_distances = np.zeros(n)
    min_PS_distances = np.full(n, np.inf)
    for j in range(2,n):
        max_PS_distances[j] = np.max(distances_to_sun[j,:])
        min_PS_distances[j] = np.min(distances_to_sun[j,:])
        
    # find earth moon apsides
    max_EM_distance = 0
    min_EM_distance = np.inf
    max_EM_distance = np.max(distances_to_earth)
    min_EM_distance = np.min(distances_to_earth)

    # output for periods, apsides...
    
    # timestep and duration
    lines_2.append("timestep is: " + str(dt) + " days , Simulation lasts for " + str(nstep * dt) + " days") 
    lines_2.append(" ")
    # energy deviation
    lines_2.append("Energy deviation is:" + str(deviation) + " M_earth AU^2 / day^2'")
    lines_2.append("which is " + str(deviation/Total_energy) +" of Total energy")
    lines_2.append(" ")
    
    # Moon apsides
    lines_2.append("Moon's aphelion: " + str(max_EM_distance) + "AU")
    lines_2.append("Moon's perihelion: " + str(min_EM_distance) + "AU")
    lines_2.append(" ")
    
    # planets upsides
    for j in range(2,n):
        if finish_orbit_index[j] >= 2:
            lines_2.append(str(solar_particles[j].label) + "'s aphelion: "+ str(max_PS_distances[j]) + "AU")
            lines_2.append(str(solar_particles[j].label) + "'s perihelion: "+ str(min_PS_distances[j]) + "AU")
        else:
            lines_2.append(str(solar_particles[j].label) + " hasn't finish a complete orbit to determine apsides")
    lines_2.append(" ")
    
    # periods
    for j in range(2,n):
        if finish_orbit_index[j] < 4:
            lines_2.append(str(solar_particles[j].label) + " needs longer simulation duration to determine its peirod")
        else:
            lines_2.append(str(solar_particles[j].label) + " has period of: " + str(periods[j]) + " days, index is: " + str(finish_orbit_index[j]))
    lines_2.append(" ")
    lines_2.append("The larger index the better accuracy")
    lines_2.append(" ")
    
    
    # output in a file as name set in console.
    with open(sys.argv[2], 'w') as file:
        for line in lines_2:
            file.write(str(line) + "\n")
        for line in lines:
            file.write(str(line) + '\n')
            
    
    pyplot.title('Halleys Comet Trajectory')
    pyplot.xlabel('x / AU')
    pyplot.ylabel('y / AU')
    pyplot.plot(positions[-1, :, 0], positions[-1, :, 1], label='Halleys Comet')
    pyplot.legend()
    pyplot.show()


   
# Add this code after the `main()` function
if __name__ == "__main__":
    main()
    

