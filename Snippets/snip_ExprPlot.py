DataFile = 'Production_Experiments.csv'
my_data = np.genfromtxt(DataFile, delimiter=',', skip_header=1).reshape(-1,7)
GCcont, Express = my_data[:,2], my_data[:,6]
plt.plot(GCcont,Express, linestyle = '--', marker = 'x', color = 'grey')
plt.gca().set(xlabel='GC-cont', ylabel='rel. expression', xlim=(.4,.8), ylim=(0,1))
plt.savefig('RelExpress_Vs_GCcont_allProm.png', format='png')
# %load Snippets/rev_ExprPlot.py 
