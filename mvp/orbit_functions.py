import numpy as np
import math
from dataclasses import dataclass


# SI
MU_EARTH = 3.986e14   # [m^3/s^2]
R_EARTH = 6378e3      # [m]



def stumpf_S(z):
    """
    Stumpf function S(z).

    TODO: figure out what tolerance for |z|>tol I actually should use
    """
    if z > 1e-12:
        sqrt_z = np.sqrt(z)
        return (sqrt_z - np.sin(sqrt_z)) / (sqrt_z**3)
    elif z < -1e-12:
        sqrt_neg_z = np.sqrt(-z)
        return (np.sinh(sqrt_neg_z) - sqrt_neg_z) / (sqrt_neg_z**3)
    else:
        return 1/6

def stumpf_C(z):
    """
    Stumpf function C(z).

    TODO: figure out what tolerance for |z|>tol I actually should use
    """
    if z > 1e-12:
        sqrt_z = np.sqrt(z)
        return (1 - np.cos(sqrt_z)) / z
    elif z < -1e-12:
        sqrt_neg_z = np.sqrt(-z)
        return (np.cosh(sqrt_neg_z) - 1) / (-z)
    else:
        return 1/2
    
def stumpf_dSdz(z):
    """
    Stumpf function S(z) derivative with respect to z.
    """
    if np.abs(z) < 1e-2:
        # use power series expansion for small z
        return -1/math.factorial(5) + 2*z/math.factorial(7) - 3*z**2/math.factorial(9) + 4*z**3/math.factorial(11)
    
    S = stumpf_S(z)
    C = stumpf_C(z)
    
    return (C - 3*S) / (2*z)

def stumpf_dCdz(z):
    """
    Stumpf function C(z) derivative with respect to z.
    """
    if np.abs(z) < 1e-2:
        # use power series expansion for small z
        return -1/math.factorial(4) + 2*z/math.factorial(6) - 3*z**2/math.factorial(8) + 4*z**3/math.factorial(10)

    S = stumpf_S(z)
    C = stumpf_C(z)

    return (1 - z*S - 2*C) / (2*z)




@dataclass
class GaussSolution:
    r1: np.ndarray
    r2: np.ndarray
    dt: float
    v1: np.ndarray
    v2: np.ndarray

def gauss_uv(r1_vec, r2_vec, dt, mu=MU_EARTH):
    # standarize inputs
    r1_vec = np.array(r1_vec, dtype=float)
    r2_vec = np.array(r2_vec, dtype=float)

    r1 = np.linalg.norm(r1_vec)
    r2 = np.linalg.norm(r2_vec)
    

    # step 0: compute nu, angle between two vectors
    #   a dot b     = |a|*|b|*cos(nu)
    #   |a cross b| = |a|*|b|*sin(nu)
    #   tan(nu) =  |a cross b| / a dot b
    nu_short = math.atan2(
        np.linalg.norm(np.cross(r1_vec, r2_vec)), 
        r1_vec.dot(r2_vec)
    )
    nu_long = 2*np.pi - nu_short

    # short way
    (f, g, gdot) = gauss_uv_fg_solver(r1, r2, dt, mu, nu_short)
    v1_short = (r2_vec - f*r1_vec) / g
    v2_short = (gdot*r2_vec - r1_vec) / g

    # long way
    (f, g, gdot) = gauss_uv_fg_solver(r1, r2, dt, mu, nu_long)
    v1_long = (r2_vec - f*r1_vec) / g
    v2_long = (gdot*r2_vec - r1_vec) / g

    # TODO: get perigee values and eliminate ones that cross through earth... or maybe do this by consumer of this function?

    sol_short = GaussSolution(
        r1_vec,
        r2_vec,
        dt,
        v1_short,
        v2_short,
    )

    sol_long = GaussSolution(
        r1_vec,
        r2_vec,
        dt,
        v1_long,
        v2_long,
    )

    return sol_short, sol_long


def gauss_uv_fg_solver(
    r1: float, 
    r2: float, 
    dt: float, 
    mu: float, 
    nu: float
) -> tuple[float, float, float]:
    # direction of motion
    DM = np.sign(np.pi - nu)

    # step 1: from r1 and r2 and "direction of motion", evaluate the constant (eq. 5-15, eq. 5-37 BMW)
    A = DM * math.sqrt(r1 * r2 * (1 + math.cos(nu)))

    # step 2: pick a trial value for z
    #         z = deltaE^2 for elliptical
    #         z = deltaF^2 for hyperolic
    #         good initial guess is z = 0
    z = 0
    dz = 0

    # newton raphson parameters
    dt_tol = 1e-4
    nr_max_iter = 300
    nr_counter = 0

    # bisection search parameters
    bs_max_iter = 100
    bs_counter = 0

    while(True):
        if nr_counter >= nr_max_iter:
            print("Warning: max newton-raphson iterations reached when solving for z")
            break
        if bs_counter >= bs_max_iter:
            print("Warning: max bisection-search iterations reached when solving for z")
            break
        
        bs_counter += 1

        z_trial = z + dz

        # DEBUG
        # print(z_trial)

        # step 3: evaluate S and C for selected z (eq. 4-37 and eq. 4-38 BMW)
        S = stumpf_S(z_trial)
        C = stumpf_C(z_trial)

        # step 4: determine aux variable y (eq. 5-17 BMW)
        y = r1 + r2 - A * (1 - z_trial*S) / math.sqrt(C)

        # step 4.5: check if y is negative (only can happen with short way trajectories)
        #           keep halving step size until y is no longer negative
        if y < 0:
            dz = 0.5*dz
            continue
        
        # z_trial is accepted
        z = z_trial

        nr_counter += 1
        bs_counter = 0

        # step 5: determine x (eq. 5-18 BMW)
        x = math.sqrt(y / C)

        # step 6: check trial value of z by computing dt_trial (eq. 5-20 BMW)
        #         then compare to true dt 
        #         then newton-raphson iterate

        dt_trial = (S*x**3 + A*math.sqrt(y)) / math.sqrt(mu)
        dt_error = dt_trial - dt

        # DEBUG
        # print(f"time error: {dt_error}")

        if np.abs(dt_error) < dt_tol:
            break
        
        Sprime = stumpf_dSdz(z)
        Cprime = stumpf_dCdz(z)

        dtdz = 1/(math.sqrt(mu)) * ( x**3 * (Sprime - 3*S*Cprime/(2*C)) + A/8*(3*S*math.sqrt(y)/C + A/x) )
        
        # clip to prevent huge jumps
        dz =  np.clip(-dt_error / dtdz, -3, 3)

    # step 7: evaluate f, g, gdot (eq. 5-21, 5-22, 5-23 BMW)
    #         then compute v1 and v2 (eq. 5-24, 5-25 BMW)

    S = stumpf_S(z)
    C = stumpf_C(z)
    y = r1 + r2 - A * (1 - z*S) / math.sqrt(C)

    f = 1 - y/r1
    g = A*math.sqrt(y/mu)
    gdot = 1 - y/r2

    return (f, g, gdot)



def propagate_orbit(r0_vec, v0_vec, dt, mu):
    """
    Propagates an orbit using the universal variable formulation for Keplers problem.
    Fundamentals of Astrodynamics (BMWS), Sec 4.4.4.
    
    Args:
        r0_vec (array-like): Initial position vector (km)
        v0_vec (array-like): Initial velocity vector (km/s)
        dt (float): Time of flight (s)
        mu (float): Gravitational parameter (km^3/s^2)
        
    Returns:
        r_vec (np.array): Final position vector
        v_vec (np.array): Final velocity vector
    """
    # edge case
    if dt == 0:
        return np.array(r0_vec, dtype=float), np.array(v0_vec, dtype=float)

    # precompute some things that are used frequently
    r0_vec = np.array(r0_vec, dtype=float)
    v0_vec = np.array(v0_vec, dtype=float)
    
    r0 = np.linalg.norm(r0_vec)
    v0 = np.linalg.norm(v0_vec)
    rv0 = np.dot(r0_vec, v0_vec)

    # alpha = 1/a (a is semi-major axis)
    alpha = -v0**2/mu + 2/r0 
    
    # initial guess for chi (universal anomaly)
    if alpha > 1e-6:
        # ellipse (BMW eq. 4.75)
        chi = np.sqrt(mu) * dt * alpha
        # print("DEBUG: ellipse: ", chi)

    elif alpha < -1e-6:
        # hyperbola (BMW eq. 4.76)
        # TODO: hyperbolas are broken right now, need to debug this further 
        #       NVM, they seem to be working but need more testing
        chi = (
            np.sign(dt) * 
            np.sqrt(-1/alpha) * 
            np.log( 
                -2 * mu * dt * alpha / 
                (rv0 + np.sign(dt) * np.sqrt(-mu*alpha) * (1 - r0*alpha))
            )
        )
        # print("DEBUG: hyperbola: ", chi)
        
    else:
        # parabola 
        #   1/6 chi^3 + sigma0/2 chi^2 + r0 chi = sqrt(mu) dt
        #   assume small step (linear term dominates)
        #   TODO: make a better guess
        chi = np.sqrt(mu) * dt / r0
        # print("DEBUG: parabola: ", chi)

    # newton-raphson
    tol = 1e-8
    max_iter = 100
    counter = 0
    while(counter < max_iter):
        counter += 1

        # compute z from chi
        z = alpha * chi**2
        S = stumpf_S(z)
        C = stumpf_C(z)
        
        # universal kepler equation (BMW eq. 4.39 / 4.41)
        #   sqrt(mu)*dt = rv0/sqrt(mu) * chi^2 * C + (1 - r0*alpha) * chi^3 * S + r0 * chi
        #   f = sqrt(mu) * (tn - t)
        f = (
            (rv0 / np.sqrt(mu)) * chi**2 * C  
            + (1 - r0 * alpha) * chi**3 * S 
            + r0 * chi 
            - np.sqrt(mu) * dt
        )
        
        # universal variable formulation (BMW eq. 4.40 / 4.44)
        #   sqrt(mu) * dt/dx = r
        #                    = chi^2 * C + rv0/sqrt(mu) * chi * (1 - zS) + r0 * (1 - zC)
        fp = C * chi**2 + (rv0 / np.sqrt(mu)) * chi * (1 - z*S) + r0 * (1 - z*C) # (1 - alpha * r0) * chi**2 * C + r0
        

        # update chi (BMW eq. 4.42)
        dchi = f/fp
        chi = chi - dchi
        
        if np.abs(dchi) < tol:
            break
    
    if counter == max_iter:
        print("Warning: max iterations reached when solving for x (universal anomaly)")

    # evaluate f and g functions, then compute r_vec
    f = 1 - (chi**2 / r0) * C            # (BMW eq. 4.58)
    g = dt - (chi**3 / np.sqrt(mu)) * S  # (BMW eq. 4.61)
    r_vec = f * r0_vec + g * v0_vec      # (BMW eq. 4.45)
    r = np.linalg.norm(r_vec)
    
    # evaluate fdot and gdot functions, then compute v_vec
    fdot = (np.sqrt(mu) / (r * r0)) * (alpha * chi**3 * S - chi)  # (BMW eq. 4.63)
    gdot = 1 - (chi**2 / r) * C                                   # (BMW eq. 4.62)
    v_vec = fdot * r0_vec + gdot * v0_vec                         # (BMW eq. 4.46)
    
    return r_vec, v_vec