Idx_optT, Linear_optT = None, None 
MB = np.mean(Biomass[Linear_optT:,Idx_optT])
GR = np.polyfit(Time[:Linear_optT],LnBiomass[:Linear_optT,Idx_optT],1)
print('max biomass: {:.0f}\nmax growth rate: {:.2f}'.format(MB,GR[0]))
# %load Snippets/rev_GrowthPars.py