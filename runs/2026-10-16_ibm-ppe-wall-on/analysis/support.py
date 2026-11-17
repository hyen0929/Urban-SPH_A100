import numpy as np

def W(q):
    q = np.asarray(q)
    out = np.zeros_like(q, dtype=float)
    m = q < 1.0
    t = 1.0 - q[m]
    out[m] = t**8 * (1 + 8*q[m] + 25*q[m]**2 + 32*q[m]**3)
    return out                      # A dropped: all results are ratios

def dWdr(q):
    q = np.asarray(q)
    out = np.zeros_like(q, dtype=float)
    m = q < 1.0
    t = 1.0 - q[m]
    out[m] = -22*q[m]*t**7*(16*q[m]**2 + 7*q[m] + 1)
    return out

def analyse(name, h, kappa, shell, d0):
    R = kappa*h
    N = 1400
    # cylindrical coords: z along wall normal (wall at z = -d0), s = radial
    z = np.linspace(-R, R, 2*N+1)
    s = np.linspace(0, R, N+1)
    Z, S = np.meshgrid(z, s, indexing='ij')
    r = np.sqrt(Z**2 + S**2)
    q = r/R
    w  = W(q)
    dw = dWdr(q)
    gz = np.where(r > 0, dw*(Z/np.maximum(r, 1e-30)), 0.0)   # z-component of grad W
    jac = 2*np.pi*S
    def integ(f, mask):
        g = np.where(mask, f*jac, 0.0)
        return np.trapezoid(np.trapezoid(g, s, axis=1), z)
    full  = r <= R
    solid = (Z < -d0) & full                       # missing side
    slab  = (Z < -d0) & (Z > -d0-shell) & full     # what one marker layer fills
    Wfull, Wsolid, Wslab = integ(w,full), integ(w,solid), integ(w,slab)
    Gsolid, Gslab = integ(gz,solid), integ(gz,slab)
    print(f"{name}  h={h} kappa={kappa} R={R} m, marker shell={shell} m, fluid at d0={d0} m")
    print(f"   Shepard  : solid side holds {Wsolid/Wfull*100:5.1f}% of the support; "
          f"markers recover {Wslab/Wsolid*100:5.1f}% of it "
          f"-> deficit {Wsolid/Wfull*100:4.1f}% -> {(Wsolid-Wslab)/Wfull*100:4.1f}%")
    print(f"   sum gradW: markers cancel {Gslab/Gsolid*100:5.1f}% of the one-sided normal term")
    print()

for d0 in (1.25, 2.50):
    analyse("Leipzig", 3.75, 2.0, 2.5, d0)
analyse("Kori   ", 7.50, 2.0, 2.5, 2.50)
