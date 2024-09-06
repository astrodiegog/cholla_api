import numpy as np
import matplotlib.pyplot as plt

from cholla_api.data.ChollaBox import ChollaBox
from cholla_api.analysis.ChollaHydroCalculator import ChollaHydroCalculator


class ChollaSnapVizHydroDiagnostic:
    '''
    Cholla Snapshot Vizualization Diagnostic object
        This object will create snapshot-wide diagnostic plots for hydro values

        Initialized with:
        - ChollaSnap (ChollaSnap): ChollaSnap object
        - ChollaGrid (ChollaGrid): ChollaGrid object
        - ImgDir (str): path of directory to place plots
        - chMacroFlags (ChollaMacroFlags): ChollaMacroFlags, holding macro
                compiling information
        - dpi (int): (optional) dots-per-inch to specify plot resolution
        - dtype (datatype): (optional) precision to calculate data 
    '''

    def __init__(self, chSnap, chGrid, chMacroFlags, ImgDir, dpi=256, precision=np.float32):
        self.snap = chSnap

        self.boxHeads = chGrid.get_BoxHeads()
        self.totnumcells = chGrid.nx_global * chGrid.ny_global * chGrid.nz_global
        self.ImgDir = ImgDir

        self.n_savedvals = 5
        _ = '''self.DE_flag = DE_flag
        self.CHEM_flag = CHEM_flag
        if self.DE_flag:
            self.n_savedvals += 1 # gas energy also saved
        if self.CHEM_flag:
            self.n_savedvals += 7 # HI, HII, HeI, HeII, HeIII, e-, temp also saved

        if self.CHEM_flag and not self.DE_flag:
            print("EXITING, chem flag needs DE flag")
            exit()
        '''
        self.MacroFlags = chMacroFlags
        if self.MacroFlags.DualEnergy:
            self.n_savedvals += 1
        if self.MacroFlags.Cooling_GPU:
            self.n_savedvals += 7 # HI, HII, HeI, HeII, HeIII, e-, temp also saved

        self.ImgDir = ImgDir
        self.precision = precision
        self.dpi = dpi

        # number of histogram bins
        self.nhist_bins = 100

    def get_alldistrfPath(self, img_ext='png'):
        '''
        Returns the file path for all distributions plot

        Args:
            img_ext (str): (optional) image file extension
        Returns:
            (str): image directory suffixed with file name
        '''

        return self.ImgDir + "/" + f"alldistr_{self.snap.SnapHead.nSnap:.0f}" + "." + img_ext

    def get_alldistr_CHEMFLAG_fPath(self, img_ext='png'):
        '''
        Returns the file path for all distributions plot

        Args:
            img_ext (str): (optional) image file extension
        Returns:
            (str): image directory suffixed with file name
        '''

        flags_str = "_"
        if self.MacroFlags.Cooling_GPU:
            flags_str += "CHEMGPU_ON"
        else:
            flags_str += "CHEMGPU_OFF"
        

        return self.ImgDir + "/" + f"alldistr_{self.snap.SnapHead.nSnap:.0f}" + flags_str + "." + img_ext


    def plot_alldistr(self, img_ext='png'):
        '''
        Plot the histograms of all hydro values for this snapshot.
        Create a (self.n_savedvals + 1 x 1) plot. From top-to-bottom, order of
            histograms are: density, x-momentum, y-momentum, z-momentum, 
            energy, gas-energy (if DE_flag was true in initialization),
            gas-energy calculation

        Currently assumes cosmological comoving units: density saved in [h2 Msun / kpc3] 
            units, momentum saved in [h2 Msun / kpc3 km / s] units, energy saved in
            [h2 Msun / kpc3 km2 / s2] units

        Args:
            save_plot (bool): (optional) whether to save plot
            img_ext (str): (optional) image file extension
        Returns:
            ...
        '''

        # set histogram limits
        l_density_histmin, l_density_histmax = -1.0, 3.0
        l_mom_histmin, l_mom_histmax = -1.0, 6.0
        l_energy_histmin, l_energy_histmax = 1.0, 8.0
        if self.MacroFlags.Cooling_GPU:
            l_gasenergy_histmin, l_gasenergy_histmax = 1.0, 8.0
            l_chemdensity_histmin, l_chemdensity_histmax = -5.0, 3.0
            l_edensity_histmin, l_edensity_histmax = -2.0, 3.0
            l_temp_histmin, l_temp_histmax = 2.0, 8.0
        else:
            l_gasenergy_histmin, l_gasenergy_histmax = -6.0, 1.0

        # set the histogram bins
        l_density_bins = np.linspace(l_density_histmin, l_density_histmax, self.nhist_bins)
        l_mom_bins = np.linspace(l_mom_histmin, l_mom_histmax, self.nhist_bins)
        l_energy_bins = np.linspace(l_energy_histmin, l_energy_histmax, self.nhist_bins)
        l_gasenergy_bins = np.linspace(l_gasenergy_histmin, l_gasenergy_histmax, self.nhist_bins)
        if self.MacroFlags.Cooling_GPU:
            l_chemdensity_bins = np.linspace(l_chemdensity_histmin, l_chemdensity_histmax, self.nhist_bins)
            l_edensity_bins = np.linspace(l_edensity_histmin, l_edensity_histmax, self.nhist_bins)
            l_temp_bins = np.linspace(l_temp_histmin, l_temp_histmax, self.nhist_bins)


        # set histogram centers + widths
        density_centers = (l_density_bins[:-1] + l_density_bins[1:]) / 2.
        density_width = np.diff(l_density_bins)
        mom_centers = (l_mom_bins[:-1] + l_mom_bins[1:]) / 2.
        mom_width = np.diff(l_mom_bins)
        energy_centers = (l_energy_bins[:-1] + l_energy_bins[1:]) / 2.
        energy_width = np.diff(l_energy_bins)
        gasenergy_centers = (l_gasenergy_bins[:-1] + l_gasenergy_bins[1:]) / 2.
        gasenergy_width = np.diff(l_gasenergy_bins)
        if self.MacroFlags.Cooling_GPU:
            chemdensity_centers = (l_chemdensity_bins[:-1] + l_chemdensity_bins[1:]) / 2.
            chemdensity_width = np.diff(l_chemdensity_bins)
            edensity_centers = (l_edensity_bins[:-1] + l_edensity_bins[1:]) / 2.
            edensity_width = np.diff(l_edensity_bins)
            temp_centers = (l_temp_bins[:-1] + l_temp_bins[1:]) / 2.
            temp_width = np.diff(l_temp_bins)


        # initialize histogram arrays
        l_density_snap = np.zeros(self.nhist_bins - 1)
        l_momx_snap = np.zeros(self.nhist_bins - 1)
        l_momy_snap = np.zeros(self.nhist_bins - 1)
        l_momz_snap = np.zeros(self.nhist_bins - 1)
        l_energy_snap = np.zeros(self.nhist_bins - 1)
        if self.MacroFlags.DualEnergy:
            l_gasenergy_snap = np.zeros(self.nhist_bins - 1)
        l_gasenergycalc_snap = np.zeros(self.nhist_bins - 1)
        if self.MacroFlags.Cooling_GPU:
            l_HIdensity_snap = np.zeros(self.nhist_bins - 1)
            l_HIIdensity_snap = np.zeros(self.nhist_bins - 1)
            l_HeIdensity_snap = np.zeros(self.nhist_bins - 1)
            l_HeIIdensity_snap = np.zeros(self.nhist_bins - 1)
            l_HeIIIdensity_snap = np.zeros(self.nhist_bins - 1)
            l_edensity_snap = np.zeros(self.nhist_bins - 1)
            l_temp_snap = np.zeros(self.nhist_bins - 1)

        for boxhead in self.boxHeads:
            # create Box object
            box = ChollaBox(self.snap.SnapPath, boxhead, self.MacroFlags)
            boxhydrocalc = ChollaHydroCalculator(box.BoxHead.local_dims)

            # ensure snap has redshift for title
            if(boxhead.nBox == 0):
                self.snap.SnapHead.set_cosmoinfo(box)

            # grab density, momxyz, energy, gasenergy
            density = box.get_hydrodata(box.density_str)
            mom_x = box.get_hydrodata(box.momx_str)
            mom_y = box.get_hydrodata(box.momy_str)
            mom_z = box.get_hydrodata(box.momz_str)
            energy = box.get_hydrodata(box.energy_str)
            if self.MacroFlags.DualEnergy:
                gasenergy = box.get_hydrodata(box.gasenergy_str)
            gasenergycalc = boxhydrocalc.int_energy(energy, density, mom_x, mom_y, mom_z)
            
            # grab chem info
            if self.MacroFlags.Cooling_GPU:
                HIdensity = box.get_hydrodata(box.HIdensity_str)
                HIIdensity = box.get_hydrodata(box.HIIdensity_str)
                HeIdensity = box.get_hydrodata(box.HeIdensity_str)
                HeIIdensity = box.get_hydrodata(box.HeIIdensity_str)
                HeIIIdensity = box.get_hydrodata(box.HeIIIdensity_str)
                edensity = box.get_hydrodata(box.edensity_str)
                temp = box.get_hydrodata(box.temp_str)


            # histogram da log
            l_density_box, _ = np.histogram(np.log10(density.flatten()), bins=l_density_bins)
            l_momx_box, _ = np.histogram(np.log10(np.abs(mom_x).flatten()), bins=l_mom_bins)
            l_momy_box, _ = np.histogram(np.log10(np.abs(mom_y).flatten()), bins=l_mom_bins)
            l_momz_box, _ = np.histogram(np.log10(np.abs(mom_z).flatten()), bins=l_mom_bins)
            l_energy_box, _ = np.histogram(np.log10(energy.flatten()), bins=l_energy_bins)
            if self.MacroFlags.DualEnergy:
                l_gasenergy_box, _ = np.histogram(np.log10(gasenergy.flatten()), bins=l_gasenergy_bins)
            l_gasenergycalc_box, _ = np.histogram(np.log10(gasenergycalc.flatten()), bins=l_gasenergy_bins)
            if self.MacroFlags.Cooling_GPU:
                l_HIdensity_box, _ = np.histogram(np.log10(HIdensity.flatten()), bins=l_chemdensity_bins)
                l_HIIdensity_box, _ = np.histogram(np.log10(HIIdensity.flatten()), bins=l_chemdensity_bins)
                l_HeIdensity_box, _ = np.histogram(np.log10(HeIdensity.flatten()), bins=l_chemdensity_bins)
                l_HeIIdensity_box, _ = np.histogram(np.log10(HeIIdensity.flatten()), bins=l_chemdensity_bins)
                l_HeIIIdensity_box, _ = np.histogram(np.log10(HeIIIdensity.flatten()), bins=l_chemdensity_bins)
                l_edensity_box, _ = np.histogram(np.log10(edensity.flatten()), bins=l_edensity_bins)
                l_temp_box, _ = np.histogram(np.log10(temp.flatten()), bins=l_temp_bins)

            # place onto global array
            l_density_snap += l_density_box
            l_momx_snap += l_momx_box
            l_momy_snap += l_momy_box
            l_momz_snap += l_momz_box
            l_energy_snap += l_energy_box
            if self.MacroFlags.DualEnergy:
                l_gasenergy_snap += l_gasenergy_box
            l_gasenergycalc_snap += l_gasenergycalc_box
            if self.MacroFlags.Cooling_GPU:
                l_HIdensity_snap += l_HIdensity_box
                l_HIIdensity_snap += l_HIIdensity_box
                l_HeIdensity_snap += l_HeIdensity_box
                l_HeIIdensity_snap += l_HeIIdensity_box
                l_HeIIIdensity_snap += l_HeIIIdensity_box
                l_edensity_snap += l_edensity_box
                l_temp_snap += l_temp_box

            print(f"Done with box {box.BoxHead.nBox:.0f}")

        # normalize all histograms to become a PDF
        l_density_snap /= self.totnumcells
        l_momx_snap /= self.totnumcells
        l_momy_snap /= self.totnumcells
        l_momz_snap /= self.totnumcells
        l_energy_snap /= self.totnumcells
        if self.MacroFlags.DualEnergy:
            l_gasenergy_snap /= self.totnumcells
        l_gasenergycalc_snap /= self.totnumcells
        if self.MacroFlags.Cooling_GPU:
            l_HIdensity_snap /= self.totnumcells
            l_HIIdensity_snap /= self.totnumcells
            l_HeIdensity_snap /= self.totnumcells
            l_HeIIdensity_snap /= self.totnumcells
            l_HeIIIdensity_snap /= self.totnumcells
            l_edensity_snap /= self.totnumcells
            l_temp_snap /= self.totnumcells


        fig, ax = plt.subplots(nrows = 6, ncols = 1, figsize=(10,10))
        ax_density, ax_momx, ax_momy, ax_momz, ax_energy, ax_gasenergycalc = ax
        
        if self.MacroFlags.DualEnergy:
            fig, ax = plt.subplots(nrows = 7, ncols = 1, figsize=(10,12))
            ax_density, ax_momx, ax_momy, ax_momz, ax_energy, ax_gasenergy, ax_gasenergycalc = ax
        if self.MacroFlags.Cooling_GPU:
            fig, ax = plt.subplots(nrows = 7, ncols = 2, figsize=(20,12))
            ax_density, ax_momx, ax_momy, ax_momz, ax_energy, ax_gasenergy, ax_gasenergycalc = ax[:,0]
            ax_HIdensity, ax_HIIdensity, ax_HeIdensity, ax_HeIIdensity, ax_HeIIIdensity, ax_edensity, ax_temp = ax[:,1]


        # do da plottin

        # plot density
        im0 = ax_density.bar(density_centers, l_density_snap, width=density_width, 
                             align='center', facecolor='none', edgecolor='k')

        # plot momx
        im1 = ax_momx.bar(mom_centers, l_momx_snap, width=mom_width, align='center', 
                          facecolor='none', edgecolor='k')
        
        # plot momy
        im2 = ax_momy.bar(mom_centers, l_momy_snap, width=mom_width, align='center', 
                          facecolor='none', edgecolor='k')

        # plot momz
        im3 = ax_momz.bar(mom_centers, l_momz_snap, width=mom_width, align='center',
                          facecolor='none', edgecolor='k')

        # plot energy
        im4 = ax_energy.bar(energy_centers, l_energy_snap, width=energy_width,
                            align='center', facecolor='none', edgecolor='k')

        if self.MacroFlags.DualEnergy:
            # plot gasenergy
            im6 = ax_gasenergy.bar(gasenergy_centers, l_gasenergy_snap, 
                                   width=gasenergy_width, align='center', facecolor='none', 
                                   edgecolor='k')
        
        # plot calculated gasenergy
        im5 = ax_gasenergycalc.bar(gasenergy_centers, l_gasenergycalc_snap,
                                   width=gasenergy_width, align='center', facecolor='none', 
                                   edgecolor='k')

        if self.MacroFlags.Cooling_GPU:
            # plot chem species
            im10 = ax_HIdensity.bar(chemdensity_centers, l_HIdensity_snap, width=chemdensity_width,
                                    align='center', facecolor='none', edgecolor='k')
            
            im11 = ax_HIIdensity.bar(chemdensity_centers, l_HIIdensity_snap, width=chemdensity_width,
                                     align='center', facecolor='none', edgecolor='k')

            im12 = ax_HeIdensity.bar(chemdensity_centers, l_HeIdensity_snap, width=chemdensity_width,
                                     align='center', facecolor='none', edgecolor='k')

            im13 = ax_HeIIdensity.bar(chemdensity_centers, l_HeIIdensity_snap, width=chemdensity_width,
                                      align='center', facecolor='none', edgecolor='k')

            im14 = ax_HeIIIdensity.bar(chemdensity_centers, l_HeIIIdensity_snap, width=chemdensity_width,
                                        align='center', facecolor='none', edgecolor='k')

            im15 = ax_edensity.bar(edensity_centers, l_edensity_snap, width=edensity_width,
                                   align='center', facecolor='none', edgecolor='k')

            im16 = ax_temp.bar(temp_centers, l_temp_snap, width=temp_width,
                               align='center', facecolor='none', edgecolor='k')


        # define x and y string labels
        density_str = r'$\rho$'
        densityunit_str = r"$\log_{10} [h^2 \rm{M}_{\odot} \rm{kpc}^{-3}]$"
        density_label = density_str + r'$\ -\ $' + densityunit_str
        density_prob = r"$\rm{P}($" + density_str + r"$)$"

        momx_str = r'$p_x$'
        momy_str = r'$p_y$'
        momz_str = r'$p_z$'
        momunit_str = r"$\log_{10} [h^2 \rm{M}_{\odot} \rm{kpc}^{-3} \rm{km} \rm{s}^{-1}]$"
        momx_label = momx_str + r'$\ -\ $' + momunit_str
        momy_label = momy_str + r'$\ -\ $' + momunit_str
        momz_label = momz_str + r'$\ -\ $' + momunit_str
        momx_prob = r"$\rm{P}($" + momx_str + r"$)$"
        momy_prob = r"$\rm{P}($" + momy_str + r"$)$"
        momz_prob = r"$\rm{P}($" + momz_str + r"$)$"

        energy_str = r'$E$'
        if self.MacroFlags.DualEnergy:
            gasenergy_str = r'$GE$'
        gasenergycalc_str = r'$GE_{\rm{calc}}$'
        energyunit_str = r"$\log_{10} [h^2 \rm{M}_{\odot} \rm{kpc}^{-3} \rm{km}^{2} \rm{s}^{-2}]$"
        energy_label = energy_str + r'$\ -\ $' + energyunit_str
        if self.MacroFlags.DualEnergy:
            gasenergy_label = gasenergy_str + r'$\ -\ $' + energyunit_str
        gasenergycalc_label = gasenergycalc_str + r'$\ -\ $' + energyunit_str
        energy_prob = r"$\rm{P}($" + energy_str + r"$)$"
        gasenergy_prob = r"$\rm{P}($" + gasenergy_str + r"$)$"
        gasenergycalc_prob = r"$\rm{P}($" + gasenergycalc_str + r"$)$"

        if self.MacroFlags.Cooling_GPU:
            HIdensity_str = r'$\rho_{\rm{HI}}$'
            HIIdensity_str = r'$\rho_{\rm{HII}}$'
            HeIdensity_str = r'$\rho_{\rm{HeI}}$'
            HeIIdensity_str = r'$\rho_{\rm{HeII}}$'
            HeIIIdensity_str = r'$\rho_{\rm{HeIII}}$'
            edensity_str = r'$\rho_{\rm{e}}$'
            tempunit_str = r"$\log_{10} [K]$"
            temp_str = r'$\rm{T}$'
            
            HIdensity_label = HIdensity_str + r'$\ -\ $' + densityunit_str
            HIIdensity_label = HIIdensity_str + r'$\ -\ $' + densityunit_str
            HeIdensity_label = HeIdensity_str + r'$\ -\ $' + densityunit_str
            HeIIdensity_label = HeIIdensity_str + r'$\ -\ $' + densityunit_str
            HeIIIdensity_label = HeIIIdensity_str + r'$\ -\ $' + densityunit_str
            edensity_label = edensity_str + r'$\ -\ $' + densityunit_str
            temp_label = temp_str + r'$\ -\ $' + tempunit_str

            HIdensity_prob = r"$\rm{P}($" + HIdensity_str + r"$)$"
            HIIdensity_prob = r"$\rm{P}($" + HIIdensity_str + r"$)$"
            HeIdensity_prob = r"$\rm{P}($" + HeIdensity_str + r"$)$"
            HeIIdensity_prob = r"$\rm{P}($" + HeIIdensity_str + r"$)$"
            HeIIIdensity_prob = r"$\rm{P}($" + HeIIIdensity_str + r"$)$"
            edensity_prob = r"$\rm{P}($" + edensity_str + r"$)$"
            temp_prob = r"$\rm{P}($" + temp_str + r"$)$"


        # apply x-labels
        _ = ax_density.set_xlabel(density_label)
        _ = ax_momx.set_xlabel(momx_label)
        _ = ax_momy.set_xlabel(momy_label)
        _ = ax_momz.set_xlabel(momz_label)
        _ = ax_energy.set_xlabel(energy_label)
        if self.MacroFlags.DualEnergy:
            _ = ax_gasenergy.set_xlabel(gasenergy_label)
        _ = ax_gasenergycalc.set_xlabel(gasenergycalc_label)
        if self.MacroFlags.Cooling_GPU:
            _ = ax_HIdensity.set_xlabel(HIdensity_label)
            _ = ax_HIIdensity.set_xlabel(HIIdensity_label)
            _ = ax_HeIdensity.set_xlabel(HeIdensity_label)
            _ = ax_HeIIdensity.set_xlabel(HeIIdensity_label)
            _ = ax_HeIIIdensity.set_xlabel(HeIIIdensity_label)
            _ = ax_edensity.set_xlabel(edensity_label)
            _ = ax_temp.set_xlabel(temp_label)

        # apply y-labels
        _ = ax_density.set_ylabel(density_prob)
        _ = ax_momx.set_ylabel(momx_prob)
        _ = ax_momy.set_ylabel(momy_prob)
        _ = ax_momz.set_ylabel(momz_prob)
        _ = ax_energy.set_ylabel(energy_prob)
        if self.MacroFlags.DualEnergy:
            _ = ax_gasenergy.set_ylabel(gasenergy_prob)
        _ = ax_gasenergycalc.set_ylabel(gasenergycalc_prob)
        if self.MacroFlags.Cooling_GPU:
            _ = ax_HIdensity.set_ylabel(HIdensity_prob)
            _ = ax_HIIdensity.set_ylabel(HIIdensity_prob)
            _ = ax_HeIdensity.set_ylabel(HeIdensity_prob)
            _ = ax_HeIIdensity.set_ylabel(HeIIdensity_prob)
            _ = ax_HeIIIdensity.set_ylabel(HeIIIdensity_prob)
            _ = ax_edensity.set_ylabel(edensity_prob)
            _ = ax_temp.set_ylabel(temp_prob)

        # add background grids
        _ = ax_density.grid(which='both', axis='both', alpha=0.3)
        _ = ax_momx.grid(which='both', axis='both', alpha=0.3)
        _ = ax_momy.grid(which='both', axis='both', alpha=0.3)
        _ = ax_momz.grid(which='both', axis='both', alpha=0.3)
        _ = ax_energy.grid(which='both', axis='both', alpha=0.3)
        if self.MacroFlags.DualEnergy:
            _ = ax_gasenergy.grid(which='both', axis='both', alpha=0.3)
        _ = ax_gasenergycalc.grid(which='both', axis='both', alpha=0.3)
        if self.MacroFlags.Cooling_GPU:
            _ = ax_HIdensity.grid(which='both', axis='both', alpha=0.3)
            _ = ax_HIIdensity.grid(which='both', axis='both', alpha=0.3)
            _ = ax_HeIdensity.grid(which='both', axis='both', alpha=0.3)
            _ = ax_HeIIdensity.grid(which='both', axis='both', alpha=0.3)
            _ = ax_HeIIIdensity.grid(which='both', axis='both', alpha=0.3)
            _ = ax_edensity.grid(which='both', axis='both', alpha=0.3)
            _ = ax_temp.grid(which='both', axis='both', alpha=0.3)


        # title shtuff
        title_str = r"$\rm{Hydro \ Info} -  $" + rf"$z = {self.snap.SnapHead.z:.4f}$"
        fig.suptitle(title_str, y=0.95)

        fig.tight_layout()

        # savefig
        #plt.savefig(self.get_alldistrfPath(img_ext), dpi=self.dpi)
        plt.savefig(self.get_alldistr_CHEMFLAG_fPath(img_ext), dpi=self.dpi)


