# Insert the correct code sequence for plotting in this cell.
# %load ../Snippets/rev_GrowthPlot.py 

Time, Biomass = my_data[:,0], my_data[:,1:]
DataFile = os.path.join('..','Data','TempGrowthExperiment.csv')
LnBiomass = np.log(Biomass)
[plt.scatter(Time, X, label=Exp) for Exp,X in enumerate(LnBiomass.T)]
plt.legend([r'{}:{}$^\circ$C'.format(Idx, T) for Idx, T in enumerate(temperatures)], bbox_to_anchor=(1.05, 1), loc='upper left'); plt.xlabel('time, h'); plt.ylabel('ln(Biomass)')
my_data = np.genfromtxt(DataFile, delimiter=',', skip_header=1)