"""
CompMod Ex2: Particle3D, a class to describe point particles in 3D space

An instance describes a particle in Euclidean 3D space: 
velocity and position are [3] arrays

Author: Jiufu Chen
Number: FILL IN YOUR MATRICULATION NUMBER HERE,
        THE ONES STARTING WITH AN S2196890.

"""
import numpy as np


class Particle3D(object):
    """
    Class to describe point-particles in 3D space

    Attributes
    ----------
    label: name of the particle
    mass: mass of the particle
    position: position of the particle
    velocity: velocity of the particle

    Methods
    -------
    __init__
    __str__
    kinetic_energy: computes the kinetic energy
    momentum: computes the linear momentum
    update_position_1st: updates the position to 1st order
    update_position_2nd: updates the position to 2nd order
    update_velocity: updates the velocity

    Static Methods
    --------------
    read_file: initializes a P3D instance from a file handle
    total_kinetic_energy: computes total K.E. of a list of particles
    com_velocity: computes centre-of-mass velocity of a list of particles
    """

    def __init__(self, label, mass, position, velocity):
        self.label=str(label)
        self.mass=mass
        self.position=np.array(position)
        self.velocity=np.array(velocity)
        """
        Initialises a particle in 3D space.

        Parameters
        ----------
        label: str
            name of the particle
        mass: float
            mass of the particle
        position: [3] float array
            position vector
        velocity: [3] float array
            velocity vector
        """
        ...

    def __str__(self):
        """
        Return an XYZ-format string. The format is
        label    x  y  z

        Returns
        -------
        str
        """
        x = self.position[0]
        y = self.position[1]
        z = self.position[2]
        Lb = str(self.label) + ' ' + str(x) + ' ' + str(y) + ' ' + str(z)
        return Lb
        ...
        return (f"{self.label} {self.position[0]} {self.position[1]} {self.position[2]}")

    def kinetic_energy(self):
        v=np.linalg.norm(self.velocity)
        ke=self.mass*v**2/2
        
        """
        Returns the kinetic energy of a Particle3D instance

        Returns
        -------
        ke: float
            1/2 m v**2
        """
        ...
        return ke  # ...

    def momentum(self):
        p = self.mass*self.velocity

        """
        set p=[] as one blank series
        i in range of self.velocity init series
        p will contain every datas of self.velocity series multiple of self.mass
        in order of p=m1v1+m2v2+m3v3
        """
        
        ...
        return p  

    def update_position_1st(self, dt):
        self.position+=dt*self.velocity

        return self.position
        """
        By equation r(t + dt) = r(t) + dt · v(t) in initial funtion
        """
        ...

    def update_position_2nd(self, dt, f):
        """
         Given a time-step and force, second-order update of the particle's position based on its velocity
         
         Parameter
         -------
         dt: float
         time-step
         f: numpy array
         3D force vector
         """
        self.position+=dt*self.velocity+(dt**2)*f/(2*self.mass)
        return self.position    
  
    '''
    def update_position_2nd(self, dt,force):
        """
        return the renew position by two times
        returns
        ---
        renewone = position0 + dt *v + dr**2 *froce/(2mass)
        """
        position0 = self.position
        v = self.velocity
        m = self.mass
        new_position = position0 + dt * v + dt**2 * force/ (2*m)
        self.position = new_position
        return new_position
        '''

    def update_velocity(self,dt,f):
        self.velocity+=(f*dt)/(self.mass)
        return self.velocity
        """
        By equation v(t + dt) = v(t) + dt · f(t)/m in initial funtion
        """
        ...

    @staticmethod
    def read_line(line):
        n=line.split()
        label=str(n[0])
        mass=float(n[1])
        position= [float(n[2]),float(n[3]),float(n[4])]
        velocity= [float(n[5]),float(n[6]),float(n[7])]
        
        """
        Creates a Particle3D instance given a line of text.

        The input line should be in the format:
        label   <mass>  <x> <y> <z>    <vx> <vy> <vz>

        Parameters
        ----------
        filename: str
            Readable file handle in the above format

        Returns
        -------
        p: Particle3D
        """
        ...
        return Particle3D(label, mass, position, velocity) # ...

    @staticmethod
    def total_kinetic_energy(particles):
        ke=0
        for i in particles:
            ke+=i.kinetic_energy()
        """
       Total KE as a series from 1 order correlating particles funtion range to i, then add
        every order of particles.kinetic_energy until i order.
        """
        ...
        return ke  # ...

    @staticmethod
    def com_velocity(particles):
        m=0
        p=0
        for i in particles:
            m+=i.mass
            p+=i.momentum()
        com_v=p/m
        """
        Computes the CoM velocity of a list of P3D's

        Parameters
        ----------
        particles: list
            A list of Particle3D instances

        Returns
        -------
        com_vel: array
            Centre-of-mass velocity
        """
        ...
        return com_v # ...
