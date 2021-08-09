import matplotlib, os
package_directory = os.path.dirname(os.path.abspath(__file__))
matplotlib.style.use(package_directory+"/matplotlib_style_file")

def grid_selfmade(ax,minor=True,color='white',linewidth=0.85):
    ax.grid(False)
    xlim = ax.get_xlim()
    ylim = ax.get_ylim()
    ax.set_xlim(xlim)
    ax.set_ylim(ylim)
    y_axis=ax.yaxis.get_view_interval()
    for xx in ax.xaxis.get_ticklocs():
        ax.plot([xx,xx],y_axis,linestyle='-',color=color,linewidth=linewidth,zorder=0)
    if minor==True:
        for xx in ax.xaxis.get_ticklocs(minor=True):
            ax.plot([xx,xx],y_axis,linestyle='-',color='#fafafa',linewidth=0.38,zorder=-1)
    x_axis=ax.xaxis.get_view_interval()
    for yy in ax.yaxis.get_ticklocs():
        ax.plot(x_axis,[yy,yy],linestyle='-',color=color,linewidth=linewidth,zorder=0)
    if minor==True:
        for yy in ax.yaxis.get_ticklocs(minor=True):
            ax.plot(x_axis,[yy,yy],linestyle='-',color='#fafafa',linewidth=0.38,zorder=-1)
            
    ax.set_xlim(xlim)
    ax.set_ylim(ylim)
    
