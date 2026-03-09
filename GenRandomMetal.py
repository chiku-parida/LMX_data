import random
from ase import Atoms
import ase
from ase.io import read, write
from ase.visualize import view
#from nglview import show_ase
#import nglview

from clease import NewStructures, Evaluate
from clease.settings import Concentration, CECrystal, CEBulk
from clease.tools import reconfigure
import numpy as np
import logging
import clease.plot_post_process as pp
import matplotlib.pyplot as plt
from ase.db import connect
logging.basicConfig(level=logging.INFO)
from clease.corr_func import CorrFunction
import os 



class CrysGen:
    def __init__(self, dbname, template, metals, halogens, concs, struct_per_gen=100, mixing_metals=True, mixing_halogens=True):
        """
        Initialize the CrysGen: CrystalStructureGenerator.

        Parameters:
        - dbname: The name of the database to store the generated structures
        - template: The template structure to use for the crystal structure generation
        - metals: List of metal elements for the cataion sites
        - halogens: List of halogen elements for the anion sites
        - concs: List of concentrations for the metals and halogens
        - mixing_metals: Boolean flag to cataion sites (default is False)
        - mixing_halogens: Boolean flag to anion sites (default is False)
        """
        self.dbname = dbname
        self.template = template
        self.metals = metals
        self.halogens = halogens
        self.concs = concs
        self.mixing_metals = mixing_metals
        self.mixing_halogens = mixing_halogens
        self.struct_per_gen = struct_per_gen
        self.structure = None

    def F11(self):
        """
        Generate a crystal structure using the specified templates, number of atoms,
        lattice constant,.
        """
        f11 = read('/home/charles/phd_project/task1/All Structure/LiInF.cif')
        if self.template == 'F11':
            if self.mixing_halogens and self.mixing_metals:
                M1 = self.metals[0]
                M2 = self.metals[1]
                X1 = self.halogens[0]
                X2 = self.halogens[1]
                C1 = self.concs[0]
                db_name = self.dbname
                os.system(f"rm {db_name}")

                #cell information

                sg = 'P 1 21/m 1'
                sgn = 11
                a, b, c, alpha, beta, gamma = f11.get_cell_lengths_and_angles()


                #wyckoff positions

                Li1 = [0.993, -0.0729, 0.8836] #--> 4f 1
                Li2 = [0.502,  -0.0729, 0.3845] #--> 4f  1
                Li3 = [1.000, 0.115, 0.6581] #--> 4f  0.5
                Li4 = [0.507, 0.114, 0.1590] #--> 4f 0.5
                In1 = [0.49583, 0.25, 0.88389] #--> 2e 1
                In2 = [1.02027, 0.25, 0.38288] #--> 2e 1
                F1 = [1.3004, 0.25, 0.7032] #--> 2e 1
                F2 = [0.7321, 0.25, 1.0509] #--> 2e 1
                F3 = [0.7247, 0.0783, 0.7956] #--> 4f 1
                F4 = [1.2401, 0.0819, 0.9579] #--> 4f 1
                F5 = [0.7911, 0.25, 0.5492] #--> 2e 1
                F6 = [0.2206, 0.25, 0.2051] #--> 2e 1
                F7 = [0.7842, 0.0802, 0.2939] #--> 4f 1
                F8 = [0.2690, 0.0794, 0.4595] #--> 4f 1






                basis_elements=[['Li'], #4  --> 4
                                ['Li'], #4 --> 4
                                ['X', 'Li'], #2  --> 4
                                ['X', 'Li'], #2  --> 4
                                [M1, M2],#['In'], #2
                                [M1, M2],#['In'], #2
                                #['F'], ['F'], ['F'], ['F'], ['F'], ['F'], ['F'], ['F']] #2
                                #[X1], [X1], [X1], [X1], 
                                [X1, X2], [X1, X2], [X1, X2], [X1, X2], [X1, X2], [X1, X2], [X1, X2], [X1, X2]] #2
                                

                grouped_basis=[[0, 1], [2, 3], [4, 5], [6, 7, 8, 9, 10, 11, 12, 13]]



                A_eq = [[1, 0, 0, 0, 0, 0, 0],
                        [0, 1, 1, 0, 0, 0, 0],
                        [0, 0, 0, 1, 1, 0, 0],
                        [0, 0, 0, 0, 0, 1, 1],
                        [0, 0, 4, 0, 0, 0, 0],
                        [0, 0, 0, 0, 0, 24, 0],
                    ]


                b_eq = [1, 1, 1, 1, 2, 24*C1]

                con = Concentration(basis_elements=basis_elements, grouped_basis=grouped_basis, A_eq=A_eq, b_eq=b_eq) 






                setting = CECrystal(
                            cellpar=[a, b, c, alpha, beta, gamma],
                            basis=[Li1, Li2, Li3, Li4, In1, In2, F1, F2, F3, F4, F5, F6, F7, F8],
                            concentration=con,
                            spacegroup=sgn,
                            size=[(3, 0, 0), (0, 2, 0), (0, 0, 1)],
                            db_name=db_name,
                            max_cluster_dia=[5.0, 5.0],
                        );
                setting.get_prim_cell_id(write_if_missing=True)

                ns = NewStructures(settings=setting, generation_number=0, struct_per_gen=self.struct_per_gen);

                ns.generate_random_structures()
            else:
                raise ValueError("Both mixing_metals and mixing_halogens must be True.")
        else:
            raise ValueError("Invalid templaTE. Use the template name 'F11'.")
        
    def Cl12(self):
        """
        Generate a crystal structure using the specified templates, number of atoms,
        lattice constant,.
        """
        cl12 = read('/home/charles/phd_project/task1/All Structure/LiScCl.cif')
        if self.template == 'Cl12':
            if self.mixing_halogens and self.mixing_metals:
                M1 = self.metals[0]
                M2 = self.metals[1]
                X1 = self.halogens[0]
                X2 = self.halogens[1]
                C1 = self.concs[0]
                db_name = self.dbname
                os.system(f"rm {db_name}")

                #cell information

                sg = 'C 1 2/m 1'
                sgn = 12
                z = 2
                a, b, c, alpha, beta, gamma = cl12.get_cell_lengths_and_angles()

                #wyckoff positions

                Li1 = [0.5, 0.8347, 0] #--> 4g 1
                Li2 = [0, 0.169, 0.5] #--> 4h 0.365
                Li3 = [0.5, 0, 0.5] #--> 2d 0.27
                Sc1 = [0, 0, 0] #--> 2a 1
                Cl1 = [0.7579, 0, 0.2313] #--> 4i 1
                Cl2 = [0.2387, 0.8391, 0.2371] #--> 8j 1


                basis_elements=[['Li'], #4  --> 4
                                ['X', 'Li'], #1.46  --> 4
                                ['X', 'Li'], #0.54  --> 2
                                [M1, M2], #2
                                [X1, X2], [X1, X2]] #12 --> 12

                grouped_basis=[[0], [1, 2], [3], [4, 5]]


                A_eq = [[1, 0, 0, 0, 0, 0, 0],
                        [0, 1, 1, 0, 0, 0, 0],
                        [0, 0, 0, 1, 1, 0, 0],
                        [0, 0, 0, 0, 0, 1, 1],
                        [0, 0, 6, 0, 0, 0, 0],
                        [0, 0, 0, 0, 0, 12, 0],
                        ]


                b_eq = [1, 1, 1, 1, 2, 12*C1]

                con = Concentration(basis_elements=basis_elements, grouped_basis=grouped_basis, A_eq=A_eq, b_eq=b_eq) 
                setting = CECrystal(
                                    cellpar=[a, b, c, alpha, beta, gamma],
                                    basis=[Li1, Li2, Li3, Sc1, Cl1, Cl2],
                                    concentration=con,
                                    spacegroup=sgn,
                                    size=[(2, 0, 0), (0, 2, 0), (0, 0, 2)],
                                    db_name=db_name,
                                    max_cluster_dia=[5.0, 5.0],
                                );  
                setting.get_prim_cell_id(write_if_missing=True)

                ns = NewStructures(settings=setting, generation_number=0, struct_per_gen=self.struct_per_gen);

                ns.generate_random_structures()
            else:
                raise ValueError("Both mixing_metals and mixing_halogens must be True.")
        else:
            raise ValueError("Invalid templaTE. Use the template name 'Cl12'.")
        
    def Br12(self):
        """
        Generate a crystal structure using the specified templates, number of atoms,
        lattice constant,.
        """
        br12 = read('/home/charles/phd_project/task1/All Structure/LiErBr.cif')
        if self.template == 'Br12':
            if self.mixing_halogens and self.mixing_metals:
                M1 = self.metals[0]
                M2 = self.metals[1]
                X1 = self.halogens[0]
                X2 = self.halogens[1]
                C1 = self.concs[0]
                db_name = self.dbname
                os.system(f"rm {db_name}")

                #cell information

                sg = 'C 1 2/m 1'
                sgn = 12
                a, b, c, alpha, beta, gamma = br12.get_cell_lengths_and_angles()

                #wyckoff positions

                Li1 = [0.5, 0.8349, 0] #--> 4g 1
                Li2 = [0, 0.155, 0.5] #--> 4h 0.5
                Er1 = [0, 0, 0] #--> 2a 1
                Br1 = [0.7492, 0, 0.2409] #--> 4i 1
                Br2 = [0.2464, 0.8338, 0.2438] #--> 8j 1


                basis_elements=[['Li'], #4  --> 4
                                ['X', 'Li'], #2  --> 4
                                [M1, M2], #2
                                [X1, X2], [X1, X2]] #12 --> 12

                grouped_basis=[[0], [1], [2], [3, 4]]

                A_eq = [[1, 0, 0, 0, 0, 0, 0],
                        [0, 1, 1, 0, 0, 0, 0],
                        [0, 0, 0, 1, 1, 0, 0],
                        [0, 0, 0, 0, 0, 1, 1],
                        [0, 0, 4, 0, 0, 0, 0],
                        [0, 0, 0, 0, 0, 12, 0],
                        ]


                b_eq = [1, 1, 1, 1, 2,12*C1]




                con = Concentration(basis_elements=basis_elements, grouped_basis=grouped_basis, A_eq=A_eq, b_eq=b_eq) 

                setting = CECrystal(
                                    cellpar=[a, b, c, alpha, beta, gamma],
                                    basis=[Li1, Li2, Er1, Br1, Br2],
                                    concentration=con,
                                    spacegroup=sgn,
                                    size=[(2, 0, 0), (0, 2, 0), (0, 0, 2)],
                                    db_name=db_name,
                                    max_cluster_dia=[5.0, 5.0],
                                    );
                setting.get_prim_cell_id(write_if_missing=True)

                ns = NewStructures(settings=setting, generation_number=0, struct_per_gen=self.struct_per_gen);

                ns.generate_random_structures()
            else:
                raise ValueError("Both mixing_metals and mixing_halogens must be True.")
        else:
            raise ValueError("Invalid templaTE. Use the template name 'Br12'.")
        
    def F15(self):
        """
        Generate a crystal structure using the specified templates, number of atoms,
        lattice constant,.
        """
        f15 = read('/home/charles/phd_project/task1/All Structure/LiGaF.cif')
        if self.template == 'F15':
            if self.mixing_halogens and self.mixing_metals:
                M1 = self.metals[0]
                M2 = self.metals[1]
                X1 = self.halogens[0]
                X2 = self.halogens[1]
                C1 = self.concs[0]
                db_name = self.dbname
                os.system(f"rm {db_name}")

                #cell information

                sg = 'C 1 2/c 1'
                sgn = 15
                a, b, c, alpha, beta, gamma = f15.get_cell_lengths_and_angles()

                #wyckoff positions

                Li1 = [0.0221, 0.3500, 0.5373] # 8 f 1.
                Li2 = [0.2028, -0.1410, 0.9317] # 8 f 1.
                Li3 = [0, 0.1885, 0.25] # 4 e 1.
                Li4 = [0.1408, -0.1477, 0.2681] # 8 f 1.
                Li5 = [0.3321, 0.2101, 0.9481] # 8 f 1.
                Ga1 = [0, 0, 0] # 4 a 1.
                Ga2 = [0.338090, 0.00276, 0.70527] # 8 f 1.
                F1 = [0.09362, 0.3407, 0.1960] # 8 f 1.
                F2 = [0.23064, 0.3408, 0.8881] # 8 f 1.
                F3 = [0.10043, 0.3425, 0.6933] # 8 f 1.
                F4 = [0.24871, 0.0205, 0.83718] # 8 f 1.
                F5 = [0.07634, 0.1616, -0.0629] # 8 f 1.
                F6 = [-0.07322, 0.1634, 0.0658] # 8 f 1.
                F7 = [0.42093, -0.0071, 0.56809] # 8 f 1.
                F8 = [0.26220, 0.1538, 0.60803] # 8 f 1.
                F9 = [0.07206, 0.0051, 0.16521] # 8 f 1.




                basis_elements=[['Li'], ['Li'], ['Li'], ['Li'], ['Li'], #36
                                [M1, M2],[M1, M2], #12
                                #['Na', 'K'], ['Na', 'K'], ['Na', 'K'], ['Na', 'K'],
                                [X1, X2], [X1, X2], [X1, X2], [X1, X2], [X1, X2], [X1, X2], [X1, X2], [X1, X2], [X1, X2]] 
                                #[X2, X1], [X2, X1], [X2, X1], [X2, X1], [X2, X1], [X2, X1]] #72

                grouped_basis=[[0, 1, 2, 3, 4], [5, 6], [7, 8, 9, 10, 11, 12, 13, 14, 15]]


                A_eq = [[1, 0, 0, 0, 0],
                        [0, 1, 1, 0, 0],
                        [0, 0, 0, 1, 1],
                        [0, 0, 0, 72, 0],
                        ]


                b_eq = [1, 1, 1, 72*C1]

                con = Concentration(basis_elements=basis_elements, grouped_basis=grouped_basis, A_eq=A_eq, b_eq=b_eq) 

                setting = CECrystal(
                                    cellpar=[a, b, c, alpha, beta, gamma],
                                    basis=[Li1, Li2, Li3, Li4, Li5, Ga1, Ga2, F1, F2, F3, F4, F5, F6, F7, F8, F9],
                                    concentration=con,
                                    spacegroup=sgn,
                                    size=[(2, 0, 0), (0, 2, 0), (0, 0, 1)],
                                    db_name=db_name,
                                    max_cluster_dia=[5.0, 5.0],
                                    );

                setting.get_prim_cell_id(write_if_missing=True)
                ns = NewStructures(settings=setting, generation_number=0, struct_per_gen=self.struct_per_gen);

                ns.generate_random_structures()
            else:
                raise ValueError("Both mixing_metals and mixing_halogens must be True.")
        else:
            raise ValueError("Invalid templaTE. Use the template name 'F15'.")
        
    def Cl62(self):
        """
        Generate a crystal structure using the specified templates, number of atoms,
        lattice constant,.
        """
        cl62 = read('/home/charles/phd_project/task1/All Structure/LiYbCl.cif')
        if self.template == 'Cl62':
            if self.mixing_halogens and self.mixing_metals:
                M1 = self.metals[0]
                M2 = self.metals[1]
                X1 = self.halogens[0]
                X2 = self.halogens[1]
                C1 = self.concs[0]
                db_name = self.dbname
                os.system(f"rm {db_name}")

                #cell information

                sg = 'P n m a'
                sgn = 62
                z = 4
                a, b, c, alpha, beta, gamma = cl62.get_cell_lengths_and_angles()

                #wyckoff positions

                Li1 = [0.8632, 0.5914, 0.4939] # 8 d 0.705
                Li2 = [0.6168, 0.4234, 0.4943] # 8 d 0.795
                Yb1 = [0.3756, 0.25, 0.0107] # 4 c 1.
                Cl1 = [0.7908, 0.4186, 0.7272] # 8 d 1.
                Cl2 = [0.4568, 0.75, 0.2287] # 4 c 1.
                Cl3 = [0.4582, 0.4183, 0.2503] # 8 d 1.
                Cl4 = [0.2048, 0.25, 0.2396] # 4 c 1.



                basis_elements=[['X', 'Li'], #5.64  --> 8
                                ['X', 'Li'], #6.36  --> 8
                                [M1, M2], #4 --> 4
                                [X1, X2], [X1, X2], [X1, X2], [X1, X2]] #12 --> 12

                grouped_basis=[[0,1], [2], [3, 4, 5, 6]]


                A_eq = [[1, 1, 0, 0, 0, 0],
                        [0, 0, 1, 1, 0, 0],
                        [0, 0, 0, 0, 1, 1],
                        [0, 16, 0, 0, 0, 0],
                        [0, 0, 0, 0, 12, 0],
                        ]


                b_eq = [1, 1, 1, 12, 12*C1]

                con = Concentration(basis_elements=basis_elements, grouped_basis=grouped_basis, A_eq=A_eq, b_eq=b_eq) 


                setting = CECrystal(
                                    cellpar=[a, b, c, alpha, beta, gamma],
                                    basis=[Li1, Li2, Yb1, Cl1, Cl2, Cl3, Cl4],
                                    concentration=con,
                                    spacegroup=sgn,
                                    size=[(1, 0, 0), (0, 1, 0), (0, 0, 2)],
                                    db_name=db_name,
                                    max_cluster_dia=[5.0, 5.0],
                                    );

                setting.get_prim_cell_id(write_if_missing=True)
                ns = NewStructures(settings=setting, generation_number=0, struct_per_gen=self.struct_per_gen);

                ns.generate_random_structures()
            else:
                raise ValueError("Both mixing_metals and mixing_halogens must be True.")
        else:
            raise ValueError("Invalid templaTE. Use the template name 'Cl62'.")
        
    def Cl164(self):
        """
        Generate a crystal structure using the specified templates, number of atoms,
        lattice constant,.
        """
        cl164 = read('/home/charles/phd_project/task1/All Structure/LiErCl.cif')
        if self.template == 'Cl164':
            if self.mixing_halogens and self.mixing_metals:
                M1 = self.metals[0]
                M2 = self.metals[1]
                X1 = self.halogens[0]
                X2 = self.halogens[1]
                C1 = self.concs[0]
                db_name = self.dbname
                os.system(f"rm {db_name}")

                #cell information

                sg = 'P -3 m 1'
                sgn = 164
                z = 3
                a, b, c, alpha, beta, gamma = cl164.get_cell_lengths_and_angles()

                #wyckoff positions

                Li1 = [0.3397, 0.3397, 0] # 6 g 1.
                Li2 = [0.2884, 0, 0.5] # 6 h 0.5
                Er1 = [0, 0, 0] # 1 a 1.
                Er2 = [0.3333, 0.6667, 0.5100] # 2 d 0.97
                Er3 = [0.3333, 0.6667, 0.0506] # 2 d 0.03
                Cl1 = [0.1135, 0.8865, 0.7681] # 6 i 1.
                Cl2 = [0.2215, 0.7785, 0.2676] # 6 i 1.
                Cl3 = [0.4454, 0.5546, 0.7559] # 6 i 1.


                basis_elements=[['Li'], #6  --> 6
                                ['X', 'Li'], #3  --> 6
                                [M1, M2], #1 --> 1
                                [M2, M1, 'X'], #1.94   -> 2
                                [M2, M1, 'X'], #0.06   -> 2
                                [X1, X2], [X1, X2], [X1, X2]] #18 --> 18

                grouped_basis=[[0], [1], [2], [3, 4], [5, 6, 7]]


                A_eq = [[1, 0, 0, 0, 0, 0, 0, 0, 0],
                        [0, 1, 1, 0, 0, 0, 0, 0, 0],
                        [0, 0, 0, 1, 0, 0, 0, 0, 0],
                        [0, 0, 0, 0, 1, 1, 1, 0, 0],
                        [0, 0, 0, 0, 0, 0, 0, 1, 1],
                        [0, 0, 6, 0, 0, 0, 0, 0, 0],
                        [0, 0, 0, 0, 0, 0, 4, 0, 0],
                        [0, 0, 0, 0, 0, 0, 0, 18, 0],
                        ]


                b_eq = [1, 1, 1, 1, 1, 3, 2, 18*C1]
        
                con = Concentration(basis_elements=basis_elements, grouped_basis=grouped_basis, A_eq=A_eq, b_eq=b_eq) 


                setting = CECrystal(
                                    cellpar=[a, b, c, alpha, beta, gamma],
                                    basis=[Li1, Li2, Er1, Er2, Er3, Cl1, Cl2, Cl3],
                                    concentration=con,
                                    spacegroup=sgn,
                                    size=[(1, 0, 0), (0, 1, 0), (0, 0, 2)],
                                    db_name=db_name,
                                    max_cluster_dia=[5.0, 5.0],
                                    );
                setting.get_prim_cell_id(write_if_missing=True)

                ns = NewStructures(settings=setting, generation_number=0, struct_per_gen=self.struct_per_gen);

                ns.generate_random_structures()
            else:
                raise ValueError("Both mixing_metals and mixing_halogens must be True.")
        else:
            raise ValueError("Invalid templaTE. Use the template name 'Cl164'.")
        
    def F165(self):
        """
        Generate a crystal structure using the specified templates, number of atoms,
        lattice constant,.
        """
        f165 = read('/home/charles/phd_project/task1/All Structure/LiScF.cif')
        if self.template == 'F165':
            if self.mixing_halogens and self.mixing_metals:
                M1 = self.metals[0]
                M2 = self.metals[1]
                X1 = self.halogens[0]
                X2 = self.halogens[1]
                C1 = self.concs[0]
                db_name = self.dbname
                os.system(f"rm {db_name}")

                #cell information

                sg = 'P -3 c 1'
                sgn = 165
                z = 6
                a, b, c, alpha, beta, gamma = f165.get_cell_lengths_and_angles()

                #wyckoff positions

                Li1 = [0, 0.7000, 0.25] # 6 f 1
                Li2 = [0.6372, 0.6574, 0.4629] # 12 g 1
                Sc1 = [0, 0, 0] # 2 b 1
                Sc2 = [0.6667, 0.3333, 0.26391] # 4 d 1
                F1 = [0.7846, 0.5535, 0.14500] # 12 g 1
                F2 = [0.5416, 0.4298, 0.38168] # 12 g 1
                F3 = [0.1179, 0.2267, 0.10763] # 12 g 1


                basis_elements=[['Li'], ['Li'], #18  --> 18
                                [M1, M2], #6 --> 6
                                [M1, M2],
                                [X1, X2], [X1, X2], [X1, X2]] #36 --> 36

                grouped_basis=[[0, 1], [2,3], [4, 5, 6]]


                A_eq = [[1, 0, 0, 0, 0],
                        [0, 1, 1, 0, 0],
                        [0, 0, 0, 1, 1],
                        [0, 0, 0, 36, 0]
                        ]


                b_eq = [1, 1, 1, 36*C1]
        
                con = Concentration(basis_elements=basis_elements, grouped_basis=grouped_basis, A_eq=A_eq, b_eq=b_eq) 


                setting = CECrystal(
                                    cellpar=[a, b, c, alpha, beta, gamma],
                                    basis=[Li1, Li2, Sc1, Sc2, F1, F2, F3],
                                    concentration=con,
                                    spacegroup=sgn,
                                    size=[(2, 0, 0), (0, 2, 0), (0, 0, 1)],
                                    db_name=db_name,
                                    max_cluster_dia=[5.0, 5.0],
                                    );

                setting.get_prim_cell_id(write_if_missing=True)
                ns = NewStructures(settings=setting, generation_number=0, struct_per_gen=self.struct_per_gen);

                ns.generate_random_structures()
            else:
                raise ValueError("Both mixing_metals and mixing_halogens must be True.")
        else:
            raise ValueError("Invalid templaTE. Use the template name 'F165'.")

