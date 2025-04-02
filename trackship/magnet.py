import numpy as np

class Magnet():

    def __init__(self, cp_core, cp_yoke, b_goal, mag_type):
        '''Default constructor
        
        :param cp_core:
            The corner points of the core.

        :param cp_yoke:
            The corner points of the yoke.

        :param B_goal:
            The desired field.

        :param type:
            The magnet type.
        '''

        self.cp_core = cp_core
        self.cp_yoke = cp_yoke
        self.b_goal = b_goal
        self.type = mag_type

        if (mag_type != 1) and (mag_type != 3):
            print(f"Magnet type {mag_type} unknown!")

    def compute_core_surface(self):
        '''Compute the surface of the core.
        
        :return:
            The surface area.
        '''

        s = np.sum(self.cp_core[:-1, 0]*self.cp_core[1:, 2]
                   - self.cp_core[1:, 0]*self.cp_core[:-1, 2])
        s += self.cp_core[-1, 0]*self.cp_core[0, 2]
        s -=  self.cp_core[0, 0]*self.cp_core[-1, 2]

        return 0.5*abs(s)

    def compute_yoke_surface(self):
        '''Compute the surface of the core.
        
        :return:
            The surface area.
        '''

        s = np.sum(self.cp_yoke[:-1, 0]*self.cp_yoke[1:, 2]
                   - self.cp_yoke[1:, 0]*self.cp_yoke[:-1, 2])
        s += self.cp_yoke[-1, 0]*self.cp_yoke[0, 2]
        s -=  self.cp_yoke[0, 0]*self.cp_yoke[-1, 2]

        return 0.5*abs(s)
    
    def compute_B_core(self):
        '''Comoute the core field
        
        :return:
            The value of the core field.
        
        '''
        if self.type == 1:
            return self.b_goal
        else:
            S_core = self.compute_core_surface()
            S_yoke = self.compute_yoke_surface()
            return self.b_goal*S_yoke/S_core
            
    def compute_B_yoke(self):
        '''Comoute the yoke field
        
        :return:
            The value of the yoke field.
        
        '''
        if self.type == 1:
            S_core = self.compute_core_surface()
            S_yoke = self.compute_yoke_surface()
            return self.b_goal*S_core/S_yoke
        else:
            return self.b_goal
        
    def plot(self, ax, B_max=1.9):
        '''Plot the configuration.
        
        :param ax:
            The matplotlib axes object.
            
        :return:
            None.
        '''

        B_core = self.compute_B_core()
        B_yoke = self.compute_B_yoke()

        if B_core/B_max > 1.0:
            o_core = 1.0
        else:
            o_core = B_core/B_max

        if B_yoke/B_max > 1.0:
            o_yoke = 1.0
        else:
            o_yoke = B_yoke/B_max

        if self.type == 1:

            ax.fill(self.cp_core[:, 0], self.cp_core[:, 2], color=[1.0, 0., 0., o_core])
            ax.fill(self.cp_yoke[:, 0], self.cp_yoke[:, 2], color=[0.0, 0., 1.0, o_yoke])
        elif self.type == 3:
            ax.fill(self.cp_core[:, 0], self.cp_core[:, 2], color=[0.0, 0., 1.0, o_core])
            ax.fill(self.cp_yoke[:, 0], self.cp_yoke[:, 2], color=[1.0, 0., 0., o_yoke])

        return None
    
    def to_snoopy(self, snoopy_dict, y_void, coil_material="copper_water_cooled.json",
                  max_turns=1, current_density=-1, NI = 10000, coil_diam=9.0,
                  insulation=0.5, yoke_spacer=5.0, material="aisi1010.json",
                  field_density=5, delta_x=0.5,	delta_y=0.5, delta_z=0.5):
        '''Translate the magnet parameters into the
        snoopy parameters.

        :return:
            A disctionary with the snoopy variables.
        '''

        if (self.cp_core.shape[0] != 4) or (self.cp_yoke.shape[0] != 4):
            print('for snoopy magnets we need 4 keypoints!')
            return None

        # sainity checks
        if abs(self.cp_yoke[0, 2] - self.cp_yoke[1, 2]) > 1e-7:
            print('The shape is not compatible for snoo.py!')

        if abs(self.cp_yoke[2, 2] - self.cp_yoke[3, 2]) > 1e-7:
            print('The shape is not compatible for snoo.py!')

        if abs(self.cp_core[0, 2] - self.cp_core[1, 2]) > 1e-7:
            print('The shape is not compatible for snoo.py!')

        if abs(self.cp_core[2, 2] - self.cp_core[3, 2]) > 1e-7:
            print('The shape is not compatible for snoo.py!')

        if self.type == 1:
            snoopy_dict["yoke_type"].append("Mag1")
        elif self.type == 3:
            snoopy_dict["yoke_type"].append("Mag3")
        else:
            print(f"Yoke type {self.type} unknown!")
            return None

        snoopy_dict["coil_material"].append(coil_material)
        snoopy_dict["max_turns"].append(max_turns)
        snoopy_dict["J_tar(A/mm2)"].append(current_density)
        snoopy_dict["NI(A)"].append(NI)
        snoopy_dict["coil_diam(mm)"].append(coil_diam)
        snoopy_dict["insulation(mm)"].append(insulation)
        snoopy_dict["yoke_spacer(mm)"].append(yoke_spacer)
        snoopy_dict["material"].append(material)
        snoopy_dict["field_density"].append(field_density)
        snoopy_dict["delta_x(m)"].append(delta_x)
        snoopy_dict["delta_y(m)"].append(delta_y)
        snoopy_dict["delta_z(m)"].append(delta_z)

        snoopy_dict["Z_pos(m)"].append(self.cp_core[0, 2])
        snoopy_dict["Z_len(m)"].append(self.cp_core[3, 2] - self.cp_core[0, 2])

        snoopy_dict["Xmgap1(m)"].append(self.cp_core[0, 0])
        snoopy_dict["Xmgap2(m)"].append(self.cp_core[3, 0])

        snoopy_dict["Xcore1(m)"].append(self.cp_core[1, 0])
        snoopy_dict["Xcore2(m)"].append(self.cp_core[2, 0])

        snoopy_dict["Xvoid1(m)"].append(self.cp_yoke[0, 0])
        snoopy_dict["Xvoid2(m)"].append(self.cp_yoke[3, 0])

        snoopy_dict["Xyoke1(m)"].append(self.cp_yoke[1, 0])
        snoopy_dict["Xyoke2(m)"].append(self.cp_yoke[2, 0])

        snoopy_dict["Ycore1(m)"].append(y_void)
        snoopy_dict["Ycore2(m)"].append(y_void)

        snoopy_dict["Yvoid1(m)"].append(y_void)
        snoopy_dict["Yvoid2(m)"].append(y_void)

        snoopy_dict["Yyoke1(m)"].append(y_void + snoopy_dict["Xcore1(m)"][-1] - snoopy_dict["Xmgap1(m)"][-1])
        snoopy_dict["Yyoke2(m)"].append(y_void + snoopy_dict["Xcore2(m)"][-1] - snoopy_dict["Xmgap2(m)"][-1])

        return