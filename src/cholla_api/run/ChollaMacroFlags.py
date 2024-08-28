class ChollaMacroFlags:
    '''
    Cholla Macro Flags object
        Holds information related to the compiler flags used for a Cholla
            executable. 
        Initialized with:
        - DualEnergy (bool): whether Gas Internal Energy saved with Dual Energy
            Formalism
        - Gravity (bool): whether gravity was solved on GPU
        - Particles (bool): whether particles were used & tracked on GPU
        - Cosmology (bool): whether equations of moation was solved in 
            comoving frame
        - Cooling_GRACKLE (bool): whether GRACKLE was used for cooling & chemistry
        - Cooling_GPU (bool): whether cooling & chemistry solved onn GPU
        - OTF_Analysis (bool): whether on-the-fly cosmological analysis was 
            performed
    '''
    def __init__(self, DualEnergy, Gravity, Particles, Cosmology, Cooling_GRACKLE, Cooling_GPU, OTF_Analysis):
        self.hydro = True
        self.DualEnergy = DualEnergy
        self.Gravity = Gravity # assume solved on GPU
        self.Particles = Particles # assume solved on GPU
        self.Cosmology = Cosmology
        self.Cooling_GRACKLE = Cooling_GRACKLE
        self.Cooling_GPU = Cooling_GPU
        self.OTF_Analysis = OTF_Analysis
        if self.Cooling_GRACKLE and self.Cooling_GPU:
            print("--- Error (!) Cooling + Chemistry can only be solved on GPU ~OR~ GRACKLE. ---")

