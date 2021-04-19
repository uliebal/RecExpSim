DataFile = 'Strain_characterization_1.csv'
my_data = np.array([np.genfromtxt(DataFile, delimiter=',', skip_header=1)])
Time, Biomass = my_data[:,0], my_data[:,1:]
LnBiomass = np.log(Biomass)
[plt.scatter(Time, X, label=Exp) for Exp,X in enumerate(LnBiomass.T)]
plt.legend()
# %load Snippets/rev_GrowthPlot.py 
