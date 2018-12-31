# By Kipling Crossing
# kip.crossing@gmail.com
import numpy as np
from Euler import rotation_matrix


class Coil(object):
    # Put some of these into init
    coil_I = 1
    frequency = 10000

    coil_position = [0,0,0]
    coil_axis = [0,0,1]

    coil_radius_in = 0.005
    coil_radius_out = 0.015
    coil_radius_ave = (coil_radius_in + coil_radius_out)/2.0
    coil_length = 0.065
    coil_depth = coil_radius_out - coil_radius_in
    wire_thickness = 0.0005
    coil_turns = (coil_length/wire_thickness) \
    *(coil_depth/wire_thickness)

    coil_inductance = ((0.8/25.4)*(coil_radius_ave**2)*(coil_turns**2)) \
    /(6*coil_radius_ave + 9*coil_length + 10*coil_depth)

    X_L = 2*np.pi*frequency*coil_inductance

    coil_resistance = 12.0

    coil_impedance = np.sqrt(coil_resistance**2 + X_L**2)

    coil_current = 1/coil_impedance
    F = coil_current*coil_turns

    _winding_sides = 128
    winding_radius = coil_radius_in

    def __init__(self):
        pass

    @property
    def permittivity(self):
        """str: permittivity of free space"""
        return 8.854187817*(10**-12)   #E0

    @property
    def permeability(self):
        """str: permeability of free space"""
        return np.pi*4*(10**-7)      #u0


    def Small_Flux_Density(self, element_location, element_length, output_location):
        element_output_vector = np.subtract(output_location, element_location)

        element_output_unit_vector = element_output_vector / np.linalg.norm(element_output_vector)
        # Doesn't work for z = 0!!! So don't consider this case.

        # Magnetic flux density
        dB = np.dot(np.dot(1.0 / (np.linalg.norm(element_output_vector) ** 2)\
        , np.cross(element_length, element_output_unit_vector))\
        ,(self.permeability*self.coil_I/(4*np.pi)))

        return dB



    def Flux_Dencity(self,flux_location):

        small_length = 2 * self.coil_radius * np.tan(np.pi / self._coil_sides)

        #Create an class that inherits numpy and includes a methof that creates the unit vector for any vector in 3D
        unit_coil_axis = self.coil_axis / np.linalg.norm(self.coil_axis)

        c = np.cross(self.coil_axis, [1, 0, 0])

        start_radius_vector = (c / np.linalg.norm(c)) * self.coil_radius

        total_flux = [0, 0, 0]
        for i in range(self._coil_sides):
            new_radius_vector = np.dot(rotation_matrix(self.coil_axis, i * (2 * np.pi / self._coil_sides)), start_radius_vector)

            small_length_location = np.add(self.coil_position, new_radius_vector)

            unit_new_radius_vector = new_radius_vector / np.linalg.norm(new_radius_vector)
            unit_small_length_vector = np.cross(unit_new_radius_vector, unit_coil_axis)

            small_length_vector = np.dot(small_length, unit_small_length_vector)

            small_flux = self.Small_Flux_Density(small_length_location, small_length_vector, flux_location)

            total_flux = np.add(total_flux, small_flux)

        return total_flux


class Mesh(object):
    vector_mesh = {}
    _zero_vector = [0,0,0]
    def __init__(self, spacing, x_start, x_finish, y_start, y_finish, z_start, z_finish):
        self.spacing = spacing
        self.x_start = x_start
        self.x_finish = x_finish
        self.y_start = y_start
        self.y_finish = y_finish
        self.z_start = z_start
        self.z_finish = z_finish


    def Generate_Mesh(self):
        for i in range(int((1/self.spacing)*self.x_start),int((1/self.spacing)*self.x_finish)):
            for j in range(int((1/self.spacing)*self.y_start),int((1/self.spacing)*self.y_finish)):
                for k in range(int((1/self.spacing)*self.y_start),int((1/self.spacing)*self.y_finish)):
                    location = (i*self.spacing,j*self.spacing,k*self.spacing)
                    self.vector_mesh[location] = self._zero_vector



T_Coil = Coil()
print(T_Coil.coil_turns)
r = T_Coil.F
print(r)
