# %load Snippets/BioMotifSim_XNOR_ori.py
# Get x and y values for plotting

x_log = np.logspace(-2, 2, 200)
y_log = np.logspace(-2, 2, 200)
x = np.linspace(0, 2, 200)
y = np.linspace(0, 2, 200)
xx, yy = np.meshgrid(x, y)
xx_log, yy_log = np.meshgrid(x_log, y_log)

# Parameters (steep Hill functions)
nx = 10
ny = 10

def myresponse(x, y, nx, ny):
    # Add equation for the interaction of transcription factors A and B based on the input x and y.
    # Use the Hill function with OR logic, single occupancy from above.
    A = None # Add your code here
    B = None # Add your code here
    return None # Add your code here

# Generate plots
p_and = xyz_im_plot(
    xx,
    yy,
    myresponse(xx, yy, nx, ny),
    xx_log,
    yy_log,
    myresponse(xx_log, yy_log, nx, ny),
    title="XNOR logic gate",
)


bokeh.io.show(
    bokeh.layouts.column(
        bokeh.layouts.row(p_and, bokeh.models.Spacer(width=30)),
    )
)
# %load Snippets/BioMotifSim_XNOR_sol.py