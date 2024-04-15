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
            if 'Earth' in line or 'Sun' in line or 'Moon' or 'Mercury' in line:  # 只保留地球、太阳和月球的数据
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
    total_potential = 0 
    '''
    initialise the value 
    '''
    
    
    for i in range(len(particles)):
        for j in range(i+1, len(particles)):
            separation = np.linalg.norm(separations[i, j]) 

            Gravity = ((8.887692593e-10 * particles[i].mass * particles[j].mass) / (separation**2)) * (separations[i, j] / separation)  

            forces[j] += Gravity
            forces[i] -= Gravity
            '''
            This is Newton's 3rd law: forces are in same magnitude but opposite direction
            '''
                
            potential = -(8.887692593e-10 * particles[j].mass * particles[i].mass) / separation
            total_potential += potential
    
    return forces, total_potential


def main():

    if len(sys.argv)!=5:
        print ("Please input at console: %run unit3.py <input> <output> <number_of_timesteps> <timestep_size(dt)>")
        sys.exit()
    '''
    sys.argv[] is used to pass arguments from an external input to a python file when running the file
    in this case, sys.argv[1] is <input> etc.
    !=5 used to check if the length equals to 5 or not
    '''
    try:
        nstep = int(sys.argv[3])
        dt = float(sys.argv[4])
    except ValueError:
        print ("Input ValueError: Please input at console: %run unit3.py <input> <output> <number_of_timesteps> <timestep_size(dt)>")
        sys.exit()
        
    '''
    # Initial conditions of the system
    '''
    particles = read_data()
    time = 0.0
    '''
    # subtract the centre-of-mass velocity
    '''
    com_velocity = Particle3D.com_velocity(particles)
    for i in range(len(particles)):
        particles[i].velocity -= com_velocity
        '''
        com_velocity is equal to total momentum over total mass
        '''
        
    # Initialise arrays 
    n = len(particles)
    times = np.zeros(nstep)
    energy = np.zeros(nstep)
    positions = np.zeros((n, nstep, 3))
    
    forces = np.zeros((n, nstep, 3))
    '''
    initialise the value
    '''
    
    separations = compute_separations(particles)
    forces, _ = compute_forces_potential(particles, separations)
    lines = []
    
    '''
    i want to use the force not force & potential, so (forces, _ means) just take the first one (force)
    '''
    
    # Main time integration loop
    for i in range(nstep):
        times[i] = time
        time += dt

        # update all particle positions     
        for j in range(n):
            particles[j].update_position_2nd(dt, forces[j])
            '''
            the update position 2nd is in particle_3D file
            ft means the force so here I use forces[j]
            '''
            
        # store particle positions in array
            positions[j, i] = particles[j].position
            '''
            position is an array
            '''
            
        #  get new separations and new forces on all particles, and the potential
        separations_new = compute_separations(particles)
        
        force_new, potential = compute_forces_potential(particles, separations_new)
        '''
        update the force and potential
        '''
            
        #  update all particle velocities
        for j in range(n):
            particles[j].update_velocity(dt, (0.5*(forces[j] + force_new[j])))
            '''
            this is the Verlet method
            '''
            
        #  replace forces with new forces for next iteration
        forces = force_new

        #  compute the kinetic energy and save the total energy
        KE = Particle3D.total_kinetic_energy(particles)        
        Total_energy = KE + potential
        energy[i] = Total_energy
        
        # write output file with correct form in each timestep.
        lines.append(n)
        lines.append("point = " + str(i+1))
        for j in range(n):
            l = str(particles[j])
            lines.append(l)
    # output file as name set in console.
    with open(sys.argv[2], 'w') as file:
        for line in lines:
            file.write(str(line) + '\n')
            '''
            \n indicates a newline
            '''
    # Make two plots, of the Mercury - Sun x distance,
    # and of the trajectory x-vs-y of the Earth's orbit.

    pyplot.title('Mercury-Sun Location')
    pyplot.xlabel('time / days')
    pyplot.ylabel('x / AU')
    pyplot.plot(times, positions[1, :, 0] - positions[0, :, 0])
    pyplot.show()

    pyplot.title('Earth Trajectory')
    pyplot.xlabel('x / AU')
    pyplot.ylabel('y / AU')
    pyplot.plot(positions[2, :, 0],  positions[2, :, 1])
    pyplot.show()

    pyplot.title('Total Energy')
    pyplot.xlabel('x / days')
    pyplot.ylabel('Energy / M_earth AU^2 / day^2')
    pyplot.plot(times, energy)
    pyplot.show()
    
    '''
    [1, :, 0] means get the first element of all columns in the second row
    positions[1, :, 0] - positions[0, :, 0] means one element minus the corresponding one
    '''


# This python standard code makes it so that the "main"
# function is only called when you run the file directly,
# not when you just import it from another python file.
if __name__ == "__main__":
    main()

