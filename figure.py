import seaborn as sns

colors=  sns.color_palette("magma",10)
# %matplotlib inline
# %config InBackend.figure_format = 'retina'
def figure_adjust(x):
  for ax in plt.gcf().get_axes():
    ax.tick_params(labelsize=x,direction='out',length=8,width=1.5,pad=8)
    ax.spines['top'].set_linewidth(standard)
    ax.spines['bottom'].set_linewidth(standard)
    ax.spines['right'].set_linewidth(standard)
    ax.spines['left'].set_linewidth(standard)
    labels = ax.get_xticklabels() + ax.get_yticklabels()
    [label.set_fontname('Times new roman') for label in labels]
    [label.set_weight('normal') for label in labels]
standard = 2
font1 = {'family': 'Times new roman','weight': 'normal','size': 24}
font2 = {'family': 'Times new roman','weight': 'normal','size': 20}
font3 = {'family': 'Times new roman','weight': 'normal','size': 16}
font4 = {'family': 'Times new roman','weight': 'normal','size': 14}
def adjust_ax(ax):
    ax.tick_params(axis='x', direction='in')
    ax.spines['top'].set_linewidth(standard)
    ax.spines['bottom'].set_linewidth(standard)
    ax.spines['right'].set_linewidth(standard)
    ax.spines['left'].set_linewidth(standard)
    ax.tick_params(axis='y', labelsize =20, direction='in', length=10, width=1.5,  right=True)
    ax.tick_params(axis='x', labelsize =20, direction='in', length=10, width=1.5,  bottom=True)
    labels = ax.get_yticklabels() + ax.get_xticklabels()
    [label.set_fontname('Times new roman') for label in labels]
    [label.set_weight('normal') for label in labels]
    ax.minorticks_on( )
    ax.tick_params(axis='y', which='minor', length=5, width= 1.5, direction='in', right=True)
    ax.tick_params(axis='x', which='minor', length=5, width= 1.5, direction='in', bottom =True)