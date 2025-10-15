import numpy as np
import scipy as sp
from scipy.constants import gravitational_constant as G
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D



def parameters(filename):
    tle = open(filename,"r")
    lines = tle.read().split('\n')
    # name = lines[0]
    e = float(lines[2][26:33])*10e-8
    mean_motion = float(lines[2][52:63])*((2*np.pi)/(24*60*60))
    mean_anom = float(lines[2][43:51])*(np.pi/180)
    inclination = float(lines[2][8:16])*(np.pi/180)
    arg_peri = float(lines[2][34:42])*(np.pi/180)
    right_ascending = float(lines[2][17:25])*(np.pi/180)

    return e,mean_motion, mean_anom, inclination, arg_peri, right_ascending


def compute_anomalies(t_sec, ma, e, n=0):
    # Compute true, eccentric and mean anomalies from mean anomaly and eccentricity
    # Propagate mean anomaly using mean motion (vector over time)
    ma_vec = n*t_sec +ma

    ## Using Newton's method of rooting finding    
    # Define Kepler's equation as a local function (for optimisation)
    def kepler_equation(E):
        return E - e*np.sin(E)  #- ma_vec[x]
    
    # Function for the derivative of Kepler's Equation
    def dkepler_equation(E):
        return 1-e*np.cos(E)
    
    # Calculating initial guess based on given mean anomaly
    if ma < np.pi:
        init_guess = ma - e/2
    elif ma > np.pi:
        init_guess = ma + e/2
    else:
        init_guess = np.pi


    ## Calc eccentric anomaly 
    E = []      # creating empty list for eccentric anomaly

    # looping through each element in mean anomaly vector
    for a in ma_vec:
        ratio = 1e-6
        tol = 1e-6
        while abs(ratio) >= tol:
            fE = kepler_equation(init_guess) - a
            dfE = dkepler_equation(init_guess)
            ratio = (fE/dfE)
            if abs(ratio) <tol:
                break
            init_guess = init_guess - ratio 
        E.append(init_guess)

    E = np.array(E) 

    ## Calc true anomaly
    tanton2 = np.sqrt((1+e)/(1-e))*np.tan(E/2)
    ta = 2*np.arctan(tanton2)


    # Wrap anomalies to -pi:pi
    def wrap_angles(angles):
        for i in range(0,len(angles)):
            x = angles[i] % (2*np.pi)
            if x > np.pi:
                x = x - 2*np.pi
            angles[i] = x
            
            
        return angles


    ta = wrap_angles(ta)
    ea = wrap_angles(E)
    ma_vec = wrap_angles(ma_vec)

    return ta, ea, ma_vec


def compute_orbital_velocity(a, e, ta, mu):

    h = np.sqrt(a*mu*(1-e**2))
    v_r = (mu/h)*e*np.sin(ta)  # Radial velocity
    v_n = (mu/h)*(1+e*np.cos(ta)) # Normal/Tangential velocity

    return v_r, v_n


def elements_to_perifocal(ta, a, e, mu, h):
    # Calc perifocal distance in m
    r = h**2/(mu*(1+e*np.cos(ta)))

    # Compute perifocal coordinates
    p = r*np.cos(ta)
    q = r*np.sin(ta)
    w = np.zeros_like(p)

    # Compute perifocal velocities
    dp = (-mu*np.sin(ta))/h
    dq = (mu*(e+np.cos(ta)))/h
    dw = np.zeros_like(p)

    return p, q, w, dp, dq, dw


def perifocal_to_eci(p, q, w, dp, dq, dw, i, raan, argp):


    # Transform coordinates to ECI frame
    x = np.zeros_like(p)
    y = np.zeros_like(p)
    z = np.zeros_like(p)

    dx = np.zeros_like(p)
    dy = np.zeros_like(p)
    dz = np.zeros_like(p)

    trans_matrix = perifocal_to_eci_matrix(i, raan, argp)
    peri_matrix = np.array([p,q,w])
    vel_peri_matrix = np.array([dp,dq,dw])
    pos_matrix = trans_matrix @ peri_matrix
    vel_matrix = trans_matrix @ vel_peri_matrix
    # Compute coordinates
    x = pos_matrix[0]
    y = pos_matrix[1]
    z = pos_matrix[2]

    dx = vel_matrix[0]
    dy = vel_matrix[1]
    dz = vel_matrix[2]
    return x, y, z, dx, dy, dz


def perifocal_to_eci_matrix(i, raan, argp):

    # Calculate transformation matrix from perifocal to ECI frame
    sraan = np.sin(raan)
    craan = np.cos(raan)
    si = np.sin(i)
    ci = np.cos(i)
    sargp = np.sin(argp)
    cargp = np.cos(argp)

    p_to_e = np.array([[craan*cargp-sraan*ci*sargp, -craan*sargp-sraan*ci*cargp, sraan*si],
                    [sraan*cargp+craan*ci*sargp, -sraan*sargp+craan*ci*cargp, -craan*si],
                    [si*sargp, si*cargp, ci]])

    return p_to_e

def magnitude_of_radius (p,q,w,t_sec):
    rmag_list = []
    for i in range(0,len(t_sec)):
        rmag = np.sqrt(p[i]**2 + q[i]**2 + w[i]**2)
        rmag_list.append(rmag)
    return rmag_list

def startup_plotting(font_size=16,line_width=1.5,output_dpi=300,tex_backend=False):

    if tex_backend:
        try:
            plt.rcParams.update({
                    "text.usetex": True,
                    "font.family": "serif",
                    "font.serif": ["Computer Modern Roman"],
                        })
        except:
            print("WARNING: LaTeX backend not configured properly. Not using.")
            plt.rcParams.update({"font.family": "serif",
                    "font.serif": ["Computer Modern Roman"],
                        })
    
    # TODO: Colour-scheme.

    # Turn on axes grids.
    plt.rcParams.update({"axes.grid" : True, 
                        "legend.framealpha": 1,
                        "legend.edgecolor": [1,1,1],
                        "lines.linewidth": line_width,
                        "savefig.dpi": output_dpi,
                        "savefig.format": 'pdf'})

    # Change default font sizes.
    plt.rc('font', size=font_size) #controls default text size
    plt.rc('axes', titlesize=font_size) #fontsize of the title
    plt.rc('axes', labelsize=font_size) #fontsize of the x and y labels
    plt.rc('xtick', labelsize=font_size) #fontsize of the x tick labels
    plt.rc('ytick', labelsize=font_size) #fontsize of the y tick labels
    plt.rc('legend', fontsize=0.85*font_size) #fontsize of the legend 