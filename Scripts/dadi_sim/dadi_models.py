import numpy as np
import dadi
import matplotlib.pyplot as plt

# --- 1. Load your REAL data SFS ---
# This is the SFS for your two populations, e.g., Utah and Mogollon Rim.
# Replace 'path/to/your/sfs.fs' with the actual file path.
# data_sfs = dadi.Spectrum.from_file('path/to/your/sfs.fs')
ns = data_sfs.sample_sizes # e.g., (20, 24) for 10 and 12 diploid individuals

# --- 2. Set up the grid for dadi's calculations ---
# The grid points should be larger than the sample sizes. A good rule of thumb
# is to try [n1+10, n2+10, n3+10] or larger for the list.
pts_l = [ns[0] + 20, ns[0] + 30, ns[0] + 40]

print(f"Data SFS loaded for sample sizes: {ns}")
print(f"Grid points for calculation: {pts_l}")

# Let's project down the SFS for speed and to remove singletons if needed.
# This is a common practice, but optional. Here we project down to 18,18.
# You should choose projection sizes appropriate for your data.
#proj_ns = (18, 18)
#data_sfs = data_sfs.project(proj_ns)
#ns = data_sfs.sample_sizes # Update sample sizes after projection
#print(f"SFS projected down to sample sizes: {ns}")


# ===========================================================================
# Model 1: Strict Isolation (SI)
# ===========================================================================
def strict_isolation(params, ns, pts):
    """
    A simple strict isolation model.
    params = (nu1, nu2, T_split)
    nu1: Size of population 1 after split.
    nu2: Size of population 2 after split.
    T_split: Time of split in units of 2*Nref generations.
    """
    nu1, nu2, T_split = params
    xx = dadi.Numerics.default_grid(pts)
    phi = dadi.PhiManip.phi_1D(xx)
    phi = dadi.PhiManip.split_1_to_2(phi, ns[0], ns[1])
    # Integrate forward in time with no migration
    phi.integrate([nu1, nu2], T_split, m=[[0, 0], [0, 0]])
    fs = dadi.Spectrum.from_phi(phi, ns, (xx, xx))
    return fs

# ===========================================================================
# Model 2: Isolation with Migration (IM)
# ===========================================================================
def isolation_with_migration(params, ns, pts):
    """
    An isolation-with-migration model with continuous, asymmetric gene flow.
    params = (nu1, nu2, T_split, m12, m21)
    nu1, nu2: Population sizes.
    T_split: Time of split.
    m12: Migration rate from pop 2 into pop 1 (2*Nref*m12).
    m21: Migration rate from pop 1 into pop 2 (2*Nref*m21).
    """
    nu1, nu2, T_split, m12, m21 = params
    xx = dadi.Numerics.default_grid(pts)
    phi = dadi.PhiManip.phi_1D(xx)
    phi = dadi.PhiManip.split_1_to_2(phi, ns[0], ns[1])
    # Define the migration matrix
    migration_matrix = [[0, m12], [m21, 0]]
    phi.integrate([nu1, nu2], T_split, m=migration_matrix)
    fs = dadi.Spectrum.from_phi(phi, ns, (xx, xx))
    return fs

# ===========================================================================
# Model 3: Secondary Contact (SC)
# ===========================================================================
def secondary_contact(params, ns, pts):
    """
    A secondary contact model with a period of isolation followed by gene flow.
    params = (nu1, nu2, T_iso, T_sc, m12, m21)
    nu1, nu2: Population sizes (assumed constant).
    T_iso: Time of strict isolation period.
    T_sc: Time of secondary contact period (total split time = T_iso + T_sc).
    m12, m21: Migration rates during secondary contact.
    """
    nu1, nu2, T_iso, T_sc, m12, m21 = params
    xx = dadi.Numerics.default_grid(pts)
    phi = dadi.PhiManip.phi_1D(xx)
    phi = dadi.PhiManip.split_1_to_2(phi, ns[0], ns[1])
    # Isolation period
    phi.integrate([nu1, nu2], T_iso, m=[[0, 0], [0, 0]])
    # Secondary contact period
    migration_matrix = [[0, m12], [m21, 0]]
    phi.integrate([nu1, nu2], T_sc, m=migration_matrix)
    fs = dadi.Spectrum.from_phi(phi, ns, (xx, xx))
    return fs

# ===========================================================================
# Model 4: Strict Isolation with a Founder Bottleneck (SI+B)
# ===========================================================================
def si_with_bottleneck(params, ns, pts):
    """
    A strict isolation model where population 1 experiences a founder bottleneck.
    params = (nuB, nu1, nu2, T_split)
    nuB: Size of population 1 during the bottleneck (instantaneous).
    nu1: Final size of population 1 after recovery.
    nu2: Size of population 2.
    T_split: Total time since the split.
    """
    nuB, nu1, nu2, T_split = params
    xx = dadi.Numerics.default_grid(pts)
    phi = dadi.PhiManip.phi_1D(xx)
    # To model an instantaneous bottleneck, we split and integrate with
    # the final sizes, but pass nuB as the initial size for pop 1.
    phi = dadi.PhiManip.split_1_to_2(phi, ns[0], ns[1])
    # This integration function allows for exponential growth/decline.
    # We are going from nuB to nu1 in pop 1 over T_split.
    # Pop 2 is constant at nu2.
    nu1_func = lambda t: nuB * (nu1/nuB)**(t/T_split)
    nu_func = (nu1_func, lambda t: nu2) # nu2 is constant
    phi.integrate(nu_func, T_split, m=[[0, 0], [0, 0]])
    fs = dadi.Spectrum.from_phi(phi, ns, (xx, xx))
    return fs


