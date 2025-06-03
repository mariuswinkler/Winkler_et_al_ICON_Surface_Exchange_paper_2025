#%%
import functions as fn
import numpy as np
import matplotlib.pylab as plt
import matplotlib.gridspec as gridspec

def compute_exchange_coefficients_for_threshold(min_wind_threshold, mwind_values):
    dz = 12.5
    pqm1 = 0.00908387     ###c.f. ERA5-tropics[975hPa] : 0.01628391 ###c.f. ngc4008-tropics[ml90]: 0.00908387
    qsat_sfc = 0.00908104 ###c.f. ERA5-tropics[1000hPa]: 0.01644187 ###c.f. ngc4008-tropics[2m]  : 0.00908104
    thetam1 = 289 
    theta_sfc = 290.0
    rough_m = 0.01 # Roughness length in meters

    arrays = {
        'cD': [], 'cH': [], 'cD_neutral': [], 'cH_neutral': [],
        'RIB': [], 'mwind_out': [], 'stab_mom': [], 'stab_heat': []
    }

    for mwind in mwind_values:
        cD, cH, cD_neutral, cH_neutral, RIB, mwind_out, stab_mom, stab_heat = fn.sfc_exchange_coefficients(
            dz, pqm1, thetam1, mwind, rough_m, theta_sfc, qsat_sfc, min_wind_threshold
        )
        arrays['cD'].append(cD)
        arrays['cH'].append(cH)
        arrays['cD_neutral'].append(cD_neutral)
        arrays['cH_neutral'].append(cH_neutral)
        arrays['RIB'].append(RIB)
        arrays['mwind_out'].append(mwind_out)
        arrays['stab_mom'].append(stab_mom)
        arrays['stab_heat'].append(stab_heat)

    # Convert lists to arrays
    for key in arrays:
        arrays[key] = np.array(arrays[key])

    # Add derived variables
    arrays['km'] = arrays['cD'] * dz * np.sqrt(mwind_values**2)
    arrays['kh'] = arrays['cH'] * dz * np.sqrt(mwind_values**2)
    arrays['pr'] = arrays['cD'] / arrays['cH']
    arrays['tau'] = arrays['cD'] * mwind_values**2

    return arrays

#%%
mwind_values = np.arange(0.01, 10, 0.01)

thresholds = {
    "1MIN": 1.0,
    "4MIN": 4.0
}

results = {key: compute_exchange_coefficients_for_threshold(thresh, mwind_values)
           for key, thresh in thresholds.items()}

#%%
## Plotting
SIZE = 15
plt.rcParams["axes.labelsize"] = SIZE
plt.rcParams["legend.fontsize"] = SIZE
plt.rcParams["xtick.labelsize"] = SIZE
plt.rcParams["ytick.labelsize"] = SIZE
plt.rcParams["font.size"] = SIZE
fig = plt.figure(figsize=(15, 18), facecolor="w", edgecolor="k")
fig.suptitle(f'Forward:\n ICON "sfc_exchange_coefficients" Routine', y=0.94)

G = gridspec.GridSpec(4, 2, hspace=0.5, wspace=0.3)

######
######

ax1 = plt.subplot(G[0, 0])

ax1.plot(mwind_values, results['1MIN']['mwind_out'], color='red', label='1MIN')
ax1.plot(mwind_values, results['4MIN']['mwind_out'], color='blue', label='4MIN')

ax1.set_title(r"Wind Speed", fontsize=SIZE)
ax1.set_ylabel(r"Processed Wind Speed / ms$^{-1}$")
ax1.set_xlabel(r"")
ax1.set_xlim(left=0, right=max(mwind_values))
ax1.set_xticks([0,5,10])
#ax1.set_xticklabels([])
ax1.set_ylim(bottom=0, top=10)
ax1.set_yticks([0,5,10])

ax1.spines[["left", "bottom"]].set_position(("outward", 20))
ax1.spines[["right", "top"]].set_visible(False)

ax1.annotate('A', xy=(0, 1), xycoords='axes fraction', fontsize=20, fontweight='bold', ha='center', va='center', xytext=(-20, 20), textcoords='offset points')
#ax1.legend(loc="upper center", bbox_to_anchor=(0.5, -0.35), fancybox=True, shadow=False, ncol=2, fontsize=SIZE)

######
######

ax2 = plt.subplot(G[0, 1])

ax2.plot(mwind_values, results['1MIN']['RIB'], color='red', label='1MIN')
ax2.plot(mwind_values, results['4MIN']['RIB'], color='blue', label='4MIN')
ax2.axhline(y=0.0, color='grey', ls='dashed', lw=1)

ax2.set_title("Richardson Number (RIB)", fontsize=SIZE)
ax2.set_ylabel(r"$/ - $")
ax2.set_xlabel(r"")
ax2.set_xlim(left=0, right=10)
ax2.set_xticks([0,5,10])
#ax2.set_xticklabels([])
ax2.set_ylim(bottom=-0.5,top=0)
#ax2.set_yticks([YLOW, 0, YTOP])

ax2.spines[["left", "bottom"]].set_position(("outward", 20))
ax2.spines[["right", "top"]].set_visible(False)

ax2.annotate('B', xy=(0, 1), xycoords='axes fraction', fontsize=20, fontweight='bold', ha='center', va='center', xytext=(-20, 20), textcoords='offset points')
#ax2.legend(loc="upper center", bbox_to_anchor=(0.5, -0.35), fancybox=True, shadow=False, ncol=2, fontsize=SIZE)

######
######

ax3 = plt.subplot(G[1, 0])

ax3.plot(mwind_values, results['1MIN']['stab_mom'], color='red', ls='solid', label='1MIN')
ax3.plot(mwind_values, results['4MIN']['stab_mom'], color='blue', ls='solid', label='4MIN')
ax3.axhline(y=1.0, color='grey', ls='dashed', lw=1)

ax3.set_title(r"Stability Function Momentum", fontsize=SIZE)
ax3.set_ylabel(r"")
ax3.set_xlabel(r"")
ax3.set_xlim(left=0, right=10)
ax3.set_xticks([0,5,10])
#ax3.set_xticklabels([])
ax3.set_ylim(bottom=0,top=2)
#ax3.set_yticks([0,0.0005,0.001,0.0015])

ax3.spines[["left", "bottom"]].set_position(("outward", 20))
ax3.spines[["right", "top"]].set_visible(False)

ax3.annotate('C', xy=(0, 1), xycoords='axes fraction', fontsize=20, fontweight='bold', ha='center', va='center', xytext=(-20, 20), textcoords='offset points')
#ax3.legend(loc="upper center", bbox_to_anchor=(0.5, -0.35), fancybox=True, shadow=False, ncol=2, fontsize=SIZE)

######
######

ax4 = plt.subplot(G[1, 1])

ax4.plot(mwind_values, results['1MIN']['cD'], color='red', label='1MIN')
ax4.plot(mwind_values, results['4MIN']['cD'], color='blue', label='4MIN')

ax4.set_title(r"$c_D$ (called km in ICON Routine)", fontsize=SIZE)
ax4.set_ylabel(r"$/ - $")
ax4.set_xlabel(r"")
ax4.set_xlim(left=0, right=10)
ax4.set_xticks([0,5,10])
#ax4.set_xticklabels([])
ax4.set_ylim(bottom=0, top=0.006)
#ax4.set_yticks([YLOW, 0, YTOP])

ax4.spines[["left", "bottom"]].set_position(("outward", 20))
ax4.spines[["right", "top"]].set_visible(False)

ax4.annotate('D', xy=(0, 1), xycoords='axes fraction', fontsize=20, fontweight='bold', ha='center', va='center', xytext=(-20, 20), textcoords='offset points')
#ax4.legend(loc="upper center", bbox_to_anchor=(-0.15, -0.2), fancybox=True, shadow=False, ncol=2, fontsize=SIZE)

######
######

ax5 = plt.subplot(G[2, 0])

ax5.plot(mwind_values, results['1MIN']['pr'], color='red', label='1MIN')
ax5.plot(mwind_values, results['4MIN']['pr'], color='blue', label='4MIN')

ax5.set_title(r"Prandtl Number: $\mathrm{Pr} = c_D / c_H$", fontsize=SIZE)
ax5.set_ylabel(r"/ $-$")
ax5.set_xlabel(r"")
ax5.set_xlim(left=0, right=10)
ax5.set_xticks([0,5,10])
#ax5.set_xticklabels([])
ax5.set_ylim(bottom=0, top=1)
#ax5.set_yticks([YLOW, 0, YTOP])

ax5.spines[["left", "bottom"]].set_position(("outward", 20))
ax5.spines[["right", "top"]].set_visible(False)

ax5.annotate('E', xy=(0, 1), xycoords='axes fraction', fontsize=20, fontweight='bold', ha='center', va='center', xytext=(-20, 20), textcoords='offset points')
#ax5.legend(loc="upper center", bbox_to_anchor=(0.5, -0.35), fancybox=True, shadow=False, ncol=2, fontsize=SIZE)

######
######

ax6 = plt.subplot(G[2, 1])

ax6.plot(mwind_values, results['1MIN']['km'], color='red', label='1MIN')
ax6.plot(mwind_values, results['4MIN']['km'], color='blue', label='4MIN')

ax6.set_title(r"km = $c_D \Delta z V_{10}$", fontsize=SIZE)
ax6.set_ylabel(r"/ m$^{2}$s$^{-2}$")
ax6.set_xlabel(r"")
ax6.set_xlim(left=0, right=10)
ax6.set_xticks([0,5,10])
#ax6.set_xticklabels([])
ax6.set_ylim(bottom=0, top=0.4)
#ax6.set_yticks([YLOW, 0, YTOP])

ax6.spines[["left", "bottom"]].set_position(("outward", 20))
ax6.spines[["right", "top"]].set_visible(False)

ax6.annotate('F', xy=(0, 1), xycoords='axes fraction', fontsize=20, fontweight='bold', ha='center', va='center', xytext=(-20, 20), textcoords='offset points')
#ax6.legend(loc="upper center", bbox_to_anchor=(0.5, -0.35), fancybox=True, shadow=False, ncol=2, fontsize=SIZE)

######
######

ax7 = plt.subplot(G[3, 0])

ax7.plot(mwind_values, results['1MIN']['tau'], color='red', label='1MIN')
ax7.plot(mwind_values, results['4MIN']['tau'], color='blue', label='4MIN')

ax7.set_title(r"$\tau = c_D \cdot V_{10}^2$", fontsize=SIZE)
ax7.set_ylabel(r"/ Nm$^{-2}$")
ax7.set_xlabel(r"Input Wind Speed $V_{10}$ / ms$^{-1}$")
ax7.set_xlim(left=0, right=10)
ax7.set_xticks([0,5,10])
#ax7.set_xticklabels([])
ax7.set_ylim(bottom=-0.001)#, top=0.05)
#ax7.set_yticks([YLOW, 0, YTOP])

ax7.spines[["left", "bottom"]].set_position(("outward", 20))
ax7.spines[["right", "top"]].set_visible(False)

ax7.annotate('G', xy=(0, 1), xycoords='axes fraction', fontsize=20, fontweight='bold', ha='center', va='center', xytext=(-20, 20), textcoords='offset points')
#ax7.legend(loc="upper center", bbox_to_anchor=(0.5, -0.35), fancybox=True, shadow=False, ncol=2, fontsize=SIZE)

######
######

ax8 = plt.subplot(G[3, 1])

ax8.plot(mwind_values, results['1MIN']['tau'], color='red', label='1MIN')
ax8.plot(mwind_values, results['4MIN']['tau'], color='blue', label='4MIN')

ax8.set_title(r"Zoom: $\tau = c_D \cdot V_{10}^2$", fontsize=SIZE)
ax8.set_ylabel(r"/ Nm$^{-2}$")
ax8.set_xlabel(r"Input Wind Speed $V_{10}$ / ms$^{-1}$")
ax8.set_xlim(left=0, right=5)
ax8.set_xticks([0,5])
#ax8.set_xticklabels([])
ax8.set_ylim(bottom=-0.001, top=0.05)
#ax8.set_yticks([YLOW, 0, YTOP])

ax8.spines[["left", "bottom"]].set_position(("outward", 20))
ax8.spines[["right", "top"]].set_visible(False)

ax8.annotate('H', xy=(0, 1), xycoords='axes fraction', fontsize=20, fontweight='bold', ha='center', va='center', xytext=(-20, 20), textcoords='offset points')
ax8.legend(loc="upper center", bbox_to_anchor=(-0.15, -0.4), fancybox=True, shadow=False, ncol=2, fontsize=SIZE)

#plt.tight_layout()

filename = f'00_components_of_subroutine.png'
plt.savefig('figs/'+filename, facecolor='white', bbox_inches='tight')

plt.show()


#%%
cD_diag_1MIN = (results['1MIN']['tau'])**(0.5) / mwind_values #mwind_out_array
cD_diag_4MIN = (results['4MIN']['tau'])**(0.5) / mwind_values #mwind_out_array_4MIN
fig = plt.figure(figsize=(15, 10), facecolor="w", edgecolor="k")
fig.suptitle(f'Forward:\n ICON "sfc_exchange_coefficients" Routine', y=1)

G = gridspec.GridSpec(2, 2, hspace=0.5, wspace=0.3)

######
######

ax = plt.subplot(G[0, 0])

ax.plot(mwind_values, results['1MIN']['cD'], color='red', label='1MIN')
ax.plot(mwind_values, results['4MIN']['cD'], color='blue', label='4MIN')

ax.set_title(r"$c_D$ from Subroutine", fontsize=SIZE)
ax.set_ylabel(r"/ $-$ ")
ax.set_xlabel(r"")
ax.set_xlim(left=0, right=10)
ax.set_xticks([0,5,10])
#ax.set_xticklabels([])
ax.set_ylim(bottom=0)#, top=0.006)
#ax.set_yticks([YLOW, 0, YTOP])

ax.spines[["left", "bottom"]].set_position(("outward", 20))
ax.spines[["right", "top"]].set_visible(False)

ax.annotate('1', xy=(0, 1), xycoords='axes fraction', fontsize=20, fontweight='bold', ha='center', va='center', xytext=(-20, 20), textcoords='offset points')

######
######

ax = plt.subplot(G[0, 1])

ax.plot(mwind_values, results['1MIN']['tau'], color='red', label='1MIN')
ax.plot(mwind_values, results['4MIN']['tau'], color='blue', label='4MIN')
#ax.plot(mwind_values, tau, color='lightsteelblue', ls='dashed', label='1MIN & 4MIN')

ax.set_title(r"$\tau_{\mathrm{diag}}$ diagnosed from $c_D$", fontsize=SIZE)
ax.set_ylabel(r"/ Nm$^{-2}$")
ax.set_xlabel(r"")
ax.set_xlim(left=0, right=10)
ax.set_xticks([0,5,10])
#ax.set_xticklabels([])
ax.set_ylim(bottom=-0.001)#, top=0.05)
#ax.set_yticks([YLOW, 0, YTOP])

ax.spines[["left", "bottom"]].set_position(("outward", 20))
ax.spines[["right", "top"]].set_visible(False)

ax.annotate('2', xy=(0, 1), xycoords='axes fraction', fontsize=20, fontweight='bold', ha='center', va='center', xytext=(-20, 20), textcoords='offset points')

######
######

ax = plt.subplot(G[1, 0])

ax.plot(mwind_values, cD_diag_1MIN, color='red', label='1MIN')
ax.plot(mwind_values, cD_diag_4MIN, color='blue', label='4MIN')
#ax.plot(mwind_values, tau, color='lightsteelblue', ls='dashed', label='1MIN & 4MIN')

ax.set_title(r"$c_{D-\mathrm{diag}}$ diagnosed via: $c_{D-\mathrm{diag}} = \tau_{\mathrm{diag}}^2 / V_{10}$", fontsize=SIZE)
ax.set_ylabel(r"/ $-$ ")
ax.set_xlabel(r"Input Wind Speed $V_{10}$ / ms$^{-1}$")
ax.set_xlim(left=0, right=10)
ax.set_xticks([0,5,10])
#ax.set_xticklabels([])
ax.set_ylim(bottom=0)#, top=0.006)
#ax.set_yticks([YLOW, 0, YTOP])

ax.spines[["left", "bottom"]].set_position(("outward", 20))
ax.spines[["right", "top"]].set_visible(False)

ax.annotate('3', xy=(0, 1), xycoords='axes fraction', fontsize=20, fontweight='bold', ha='center', va='center', xytext=(-20, 20), textcoords='offset points')

ax.legend(loc="upper center", bbox_to_anchor=(1.2, -0.35), fancybox=True, shadow=False, ncol=2, fontsize=SIZE)

######
######

ax = plt.subplot(G[1, 1])

ax.plot(mwind_values, results['1MIN']['tau'], color='red', label='1MIN')
ax.plot(mwind_values, results['4MIN']['tau'], color='blue', label='4MIN')
#ax.plot(mwind_values, tau, color='lightsteelblue', ls='dashed', label='1MIN & 4MIN')

ax.set_title(r"Zoom of (2)", fontsize=SIZE)
ax.set_ylabel(r"/ Nm$^{-2}$")
ax.set_xlabel(r"Input Wind Speed $V_{10}$ / ms$^{-1}$")
ax.set_xlim(left=0, right=5)
ax.set_xticks([0,5])
#ax.set_xticklabels([])
ax.set_ylim(bottom=-0.001, top=0.05)
#ax.set_yticks([YLOW, 0, YTOP])

ax.spines[["left", "bottom"]].set_position(("outward", 20))
ax.spines[["right", "top"]].set_visible(False)

ax.annotate('4', xy=(0, 1), xycoords='axes fraction', fontsize=20, fontweight='bold', ha='center', va='center', xytext=(-20, 20), textcoords='offset points')

#plt.tight_layout()

filename = f'01_cD_vs_tau.png'
plt.savefig('figs/'+filename, facecolor='white', bbox_inches='tight')

plt.show()


# %%
