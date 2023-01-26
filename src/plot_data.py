import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
from constants import ASSETS_FOLDER_NAME
from os import path
from utils import filepath_or_default

def plot_data(Rftinf, Rft0, S, fobj_last_value, fitted_data, t_list=None, Q_estimated=None, Tho_estimated=None, Tco_estimated=None, Ud_estimated=None, dirpath=None, filename=None):
    tf_was_given = t_list is not None
    number_of_plots = 5 if tf_was_given else 1   
    
    fig = plt.figure(figsize=(10, 20 if tf_was_given else 10))
    
    ax = fig.add_subplot(number_of_plots, 1, 1)
    ax.plot(fitted_data.index, fitted_data["Rft"], linewidth=1.2, label="Given")
    ax.plot(fitted_data.index, fitted_data["Rft_fitted"], alpha = 0.5, label = 
            r"Fitted: $R_{ft}^{\infty}$ = " + f"{Rftinf:.3E}" +
            r", $R_{ft}^0$ = " + f"{Rft0:.3E}" + 
            r", $S$ = " + f"{S:.3E}")
    ax.set_ylabel(r"$R_{ft} \, (m^2K / W)$")
    ax.set_xlabel("Time (days)")
    ticks_loc = ax.get_yticks().tolist()
    ax.yaxis.set_major_locator(mticker.FixedLocator(ticks_loc))
    ax.set_yticklabels([f"{x:.1E}" for x in ticks_loc])
    ax.set_title(r"Total fouling factor exponential fit by minimizing the sum of normalized squared residuals:"+
                 "\n"+
                 r"$\sum_{i=1}^{\# \, of \, given \, days} \left[ \frac{\left(R_{ft, given, i} - R_{ft, fitted, i}\right)^2}{R_{ft, given, i}} \right] =$" +
                 f"${fobj_last_value:.1E}$",
                 pad = 20)
    ax.legend()
    
    if tf_was_given:
        ax2 = fig.add_subplot(number_of_plots, 1, 2)
        ax2.plot(t_list, Tho_estimated, "r")
        ax2.set_ylabel("Temperature (K)")
        ax2.set_xlabel("Time (days)")
        ax2.set_title("Estimated output hot temperature over time")
        
        ax3 = fig.add_subplot(number_of_plots, 1, 3)
        ax3.plot(t_list, Tco_estimated, "b")
        ax3.set_ylabel("Temperature (K)")
        ax3.set_xlabel("Time (days)")
        ax3.set_title("Estimated output cold temperature over time")
        
        ax4 = fig.add_subplot(number_of_plots, 1, 4)
        ax4.plot(t_list, Ud_estimated, "k")
        ax4.set_ylabel(r"U $(W/m^2K)$")
        ax4.set_xlabel("Time (days)")
        ax4.set_title("Estimated global heat transfer coefficient over time")
        
        ax5 = fig.add_subplot(number_of_plots, 1, 5)
        ax5.plot(t_list, Q_estimated, "orange")
        ax5.set_ylabel("Q (W)")
        ax5.set_xlabel("Time (days)")
        ax5.set_title("Estimated heat transfer rate over time")
        
    fig.tight_layout()
    
    filepath = filepath_or_default(filename, dirpath)
    
    fig.savefig(filepath, dpi=800)