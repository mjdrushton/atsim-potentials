import collections

Element_Data = collections.namedtuple('Element_Data', ['atomic_number', 'atomic_mass','covalent_radius'])

reference_data = {
  'Ru' : Element_Data(atomic_number = 44, atomic_mass = 101.06999999999999, covalent_radius = 1.46),
  'Re' : Element_Data(atomic_number = 75, atomic_mass = 186.20699999999999, covalent_radius = 1.51),
  'Ra' : Element_Data(atomic_number = 88, atomic_mass = 226.0, covalent_radius              = 2.21),
  'Rb' : Element_Data(atomic_number = 37, atomic_mass = 85.467799999999997, covalent_radius = 2.2000000000000002),
  'Rn' : Element_Data(atomic_number = 86, atomic_mass = 222.0, covalent_radius              = 1.5),
  'Rh' : Element_Data(atomic_number = 45, atomic_mass = 102.9055, covalent_radius           = 1.4199999999999999),
  'Be' : Element_Data(atomic_number = 4, atomic_mass  = 9.0121819999999992, covalent_radius = 0.95999999999999996),
  'Ba' : Element_Data(atomic_number = 56, atomic_mass = 137.327, covalent_radius            = 2.1499999999999999),
  'Bi' : Element_Data(atomic_number = 83, atomic_mass = 208.98038, covalent_radius          = 1.48),
  'Br' : Element_Data(atomic_number = 35, atomic_mass = 79.903999999999996, covalent_radius = 1.2),
  'H'  : Element_Data(atomic_number = 1, atomic_mass  = 1.0079400000000001, covalent_radius = 0.31),
  'P'  : Element_Data(atomic_number = 15, atomic_mass = 30.973761, covalent_radius          = 1.0700000000000001),
  'Os' : Element_Data(atomic_number = 76, atomic_mass = 190.22999999999999, covalent_radius = 1.4399999999999999),
  'Ge' : Element_Data(atomic_number = 32, atomic_mass = 72.640000000000001, covalent_radius = 1.2),
  'Gd' : Element_Data(atomic_number = 64, atomic_mass = 157.25, covalent_radius             = 1.96),
  'Ga' : Element_Data(atomic_number = 31, atomic_mass = 69.722999999999999, covalent_radius = 1.22),
  'Pr' : Element_Data(atomic_number = 59, atomic_mass = 140.90764999999999, covalent_radius = 2.0299999999999998),
  'Pt' : Element_Data(atomic_number = 78, atomic_mass = 195.078, covalent_radius            = 1.3600000000000001),
  'Pu' : Element_Data(atomic_number = 94, atomic_mass = 244.0, covalent_radius              = 1.8700000000000001),
  'C'  : Element_Data(atomic_number = 6, atomic_mass  = 12.0107, covalent_radius            = 0.76000000000000001),
  'Pb' : Element_Data(atomic_number = 82, atomic_mass = 207.19999999999999, covalent_radius = 1.46),
  'Pa' : Element_Data(atomic_number = 91, atomic_mass = 231.03587999999999, covalent_radius = 2.0),
  'Pd' : Element_Data(atomic_number = 46, atomic_mass = 106.42, covalent_radius             = 1.3899999999999999),
  'Cd' : Element_Data(atomic_number = 48, atomic_mass = 112.411, covalent_radius            = 1.4399999999999999),
  'Po' : Element_Data(atomic_number = 84, atomic_mass = 209.0, covalent_radius              = 1.3999999999999999),
  'Pm' : Element_Data(atomic_number = 61, atomic_mass = 145.0, covalent_radius              = 1.99),
  'Ho' : Element_Data(atomic_number = 67, atomic_mass = 164.93031999999999, covalent_radius = 1.9199999999999999),
  'Hf' : Element_Data(atomic_number = 72, atomic_mass = 178.49000000000001, covalent_radius = 1.75),
  'Hg' : Element_Data(atomic_number = 80, atomic_mass = 200.59, covalent_radius             = 1.3200000000000001),
  'He' : Element_Data(atomic_number = 2, atomic_mass  = 4.0026020000000004, covalent_radius = 0.28000000000000003),
  'Mg' : Element_Data(atomic_number = 12, atomic_mass = 24.305, covalent_radius             = 1.4099999999999999),
  'K'  : Element_Data(atomic_number = 19, atomic_mass = 39.098300000000002, covalent_radius = 2.0299999999999998),
  'Mn' : Element_Data(atomic_number = 25, atomic_mass = 54.938048999999999, covalent_radius = 1.6100000000000001),
  'O'  : Element_Data(atomic_number = 8, atomic_mass  = 15.9994, covalent_radius            = 0.66000000000000003),
  'S'  : Element_Data(atomic_number = 16, atomic_mass = 32.064999999999998, covalent_radius = 1.05),
  'W'  : Element_Data(atomic_number = 74, atomic_mass = 183.84, covalent_radius             = 1.6200000000000001),
  'Zn' : Element_Data(atomic_number = 30, atomic_mass = 65.409000000000006, covalent_radius = 1.22),
  'Eu' : Element_Data(atomic_number = 63, atomic_mass = 151.964, covalent_radius            = 1.98),
  'Zr' : Element_Data(atomic_number = 40, atomic_mass = 91.224000000000004, covalent_radius = 1.75),
  'Er' : Element_Data(atomic_number = 68, atomic_mass = 167.25899999999999, covalent_radius = 1.8899999999999999),
  'Ni' : Element_Data(atomic_number = 28, atomic_mass = 58.693399999999997, covalent_radius = 1.24),
  'Na' : Element_Data(atomic_number = 11, atomic_mass = 22.98977, covalent_radius           = 1.6599999999999999),
  'Nb' : Element_Data(atomic_number = 41, atomic_mass = 92.906379999999999, covalent_radius = 1.6399999999999999),
  'Nd' : Element_Data(atomic_number = 60, atomic_mass = 144.24000000000001, covalent_radius = 2.0099999999999998),
  'Ne' : Element_Data(atomic_number = 10, atomic_mass = 20.1797, covalent_radius            = 0.57999999999999996),
  'Np' : Element_Data(atomic_number = 93, atomic_mass = 237.0, covalent_radius              = 1.8999999999999999),
  'Fr' : Element_Data(atomic_number = 87, atomic_mass = 223.0, covalent_radius              = 2.6000000000000001),
  'Fe' : Element_Data(atomic_number = 26, atomic_mass = 55.844999999999999, covalent_radius = 1.52),
  'B'  : Element_Data(atomic_number = 5, atomic_mass  = 10.811, covalent_radius             = 0.83999999999999997),
  'F'  : Element_Data(atomic_number = 9, atomic_mass  = 18.998403199999998, covalent_radius = 0.56999999999999995),
  'Sr' : Element_Data(atomic_number = 38, atomic_mass = 87.620000000000005, covalent_radius = 1.95),
  'N'  : Element_Data(atomic_number = 7, atomic_mass  = 14.0067, covalent_radius            = 0.70999999999999996),
  'Kr' : Element_Data(atomic_number = 36, atomic_mass = 83.798000000000002, covalent_radius = 1.1599999999999999),
  'Si' : Element_Data(atomic_number = 14, atomic_mass = 28.0855, covalent_radius            = 1.1100000000000001),
  'Sn' : Element_Data(atomic_number = 50, atomic_mass = 118.70999999999999, covalent_radius = 1.3899999999999999),
  'Sm' : Element_Data(atomic_number = 62, atomic_mass = 150.36000000000001, covalent_radius = 1.98),
  'V'  : Element_Data(atomic_number = 23, atomic_mass = 50.941499999999998, covalent_radius = 1.53),
  'Sc' : Element_Data(atomic_number = 21, atomic_mass = 44.955910000000003, covalent_radius = 1.7),
  'Sb' : Element_Data(atomic_number = 51, atomic_mass = 121.76000000000001, covalent_radius = 1.3899999999999999),
  'Se' : Element_Data(atomic_number = 34, atomic_mass = 78.959999999999994, covalent_radius = 1.2),
  'Co' : Element_Data(atomic_number = 27, atomic_mass = 58.933199999999999, covalent_radius = 1.5),
  'Cm' : Element_Data(atomic_number = 96, atomic_mass = 247.0, covalent_radius              = 1.6899999999999999),
  'Cl' : Element_Data(atomic_number = 17, atomic_mass = 35.453000000000003, covalent_radius = 1.02),
  'Ca' : Element_Data(atomic_number = 20, atomic_mass = 40.078000000000003, covalent_radius = 1.76),
  'Ce' : Element_Data(atomic_number = 58, atomic_mass = 140.11600000000001, covalent_radius = 2.04),
  'Xe' : Element_Data(atomic_number = 54, atomic_mass = 131.29300000000001, covalent_radius = 1.3999999999999999),
  'Tm' : Element_Data(atomic_number = 69, atomic_mass = 168.93421000000001, covalent_radius = 1.8999999999999999),
  'Cs' : Element_Data(atomic_number = 55, atomic_mass = 132.90545, covalent_radius          = 2.4399999999999999),
  'Cr' : Element_Data(atomic_number = 24, atomic_mass = 51.996099999999998, covalent_radius = 1.3899999999999999),
  'Cu' : Element_Data(atomic_number = 29, atomic_mass = 63.545999999999999, covalent_radius = 1.3200000000000001),
  'La' : Element_Data(atomic_number = 57, atomic_mass = 138.90549999999999, covalent_radius = 2.0699999999999998),
  'Li' : Element_Data(atomic_number = 3, atomic_mass  = 6.9409999999999998, covalent_radius = 1.28),
  'Tl' : Element_Data(atomic_number = 81, atomic_mass = 204.38329999999999, covalent_radius = 1.45),
  'Lu' : Element_Data(atomic_number = 71, atomic_mass = 174.96700000000001, covalent_radius = 1.8700000000000001),
  'Th' : Element_Data(atomic_number = 90, atomic_mass = 232.03809999999999, covalent_radius = 2.0600000000000001),
  'Ti' : Element_Data(atomic_number = 22, atomic_mass = 47.866999999999997, covalent_radius = 1.6000000000000001),
  'Te' : Element_Data(atomic_number = 52, atomic_mass = 127.59999999999999, covalent_radius = 1.3799999999999999),
  'Tb' : Element_Data(atomic_number = 65, atomic_mass = 158.92534000000001, covalent_radius = 1.9399999999999999),
  'Tc' : Element_Data(atomic_number = 43, atomic_mass = 98.0, covalent_radius               = 1.47),
  'Ta' : Element_Data(atomic_number = 73, atomic_mass = 180.9479, covalent_radius           = 1.7),
  'Yb' : Element_Data(atomic_number = 70, atomic_mass = 173.03999999999999, covalent_radius = 1.8700000000000001),
  'Dy' : Element_Data(atomic_number = 66, atomic_mass = 162.5, covalent_radius              = 1.9199999999999999),
  'I'  : Element_Data(atomic_number = 53, atomic_mass = 126.90447, covalent_radius          = 1.3899999999999999),
  'U'  : Element_Data(atomic_number = 92, atomic_mass = 238.02891, covalent_radius          = 1.96),
  'Y'  : Element_Data(atomic_number = 39, atomic_mass = 88.905850000000001, covalent_radius = 1.8999999999999999),
  'Ac' : Element_Data(atomic_number = 89, atomic_mass = 227.0, covalent_radius              = 2.1499999999999999),
  'Ag' : Element_Data(atomic_number = 47, atomic_mass = 107.8682, covalent_radius           = 1.45),
  'Ir' : Element_Data(atomic_number = 77, atomic_mass = 192.21700000000001, covalent_radius = 1.4099999999999999),
  'Am' : Element_Data(atomic_number = 95, atomic_mass = 243.0, covalent_radius              = 1.8),
  'Al' : Element_Data(atomic_number = 13, atomic_mass = 26.981538, covalent_radius          = 1.21),
  'As' : Element_Data(atomic_number = 33, atomic_mass = 74.921599999999998, covalent_radius = 1.1899999999999999),
  'Ar' : Element_Data(atomic_number = 18, atomic_mass = 39.948, covalent_radius             = 1.0600000000000001),
  'Au' : Element_Data(atomic_number = 79, atomic_mass = 196.96655000000001, covalent_radius = 1.3600000000000001),
  'At' : Element_Data(atomic_number = 85, atomic_mass = 210.0, covalent_radius              = 1.5),
  'In' : Element_Data(atomic_number = 49, atomic_mass = 114.818, covalent_radius            = 1.4199999999999999),
  'Mo' : Element_Data(atomic_number = 42, atomic_mass = 95.939999999999998, covalent_radius = 1.54)
  }
