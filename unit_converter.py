# code unit in cgs
l_u = 1.5e13  # [cm]
m_u = 2e33  # [g]
G_ = 6.67e-8
v_u = (G_*m_u/l_u)**0.5  # [cm/s]

t_u = l_u/v_u
dm_u = m_u/t_u
den_u = m_u/l_u**3
p_u = m_u/l_u/t_u**2
B_u = l_u**-0.5 * m_u**0.5 * t_u**-1
L_u = m_u * l_u**2 / t_u**3


unit_list = [
    'time: '+str(t_u),
    'length: '+str(l_u),
    'mass: '+str(m_u),
    'density: '+str(den_u),
    'mass loss: '+str(dm_u),
    'velocity: '+str(v_u),
    'pressure: '+str(p_u),
    'magnetic field: '+str(B_u),
    'luminosity: '+str(L_u)
]
print(unit_list)


def convert_l(_l): return _l/l_u


_a = 3.5e12
a_binary = convert_l(_a)

_R_s = 7e10  # r_sun (stellar radius)
R_s = convert_l(_R_s)

_R_w_s = 7e11  # r_sun (wind-launchiung surface)
R_w_s = convert_l(_R_w_s)

_R_NS = 1e5  # 10 KM
R_NS = convert_l(_R_NS)

_R_w_NS = 7e11
R_w_NS = convert_l(_R_w_NS)

# MASS [g]


def convert_m(_m): return _m/m_u


_m_star = 3e34
m_star = convert_m(_m_star)
_m_NS = 2.8e33
m_NS = convert_m(_m_NS)

# MASS LOSS


def convert_dm(_dm): return _dm/dm_u


_Mdot_s = 6.3e19
Mdot_s = convert_dm(_Mdot_s)

# VELOCITY


def convert_v(_v): return _v/v_u


_v_w_s = 3e8  # 0.01c
v_w_s = convert_v(_v_w_s)  # stellar wind speed at the surface of the star

# LUMINOSITY


def convert_L(_L): return _L/L_u


_L_w_NS = 1e34
L_w_NS = convert_L(_L_w_NS)

# MAGNETIC FIELD


def convert_B(_B): return _B/B_u


_B_star = 10 * (_R_s/_R_w_s)**2
B_star = convert_B(_B_star)
_B_NS = 1e14 * (_R_NS/_R_w_NS)**2
B_NS = convert_B(_B_NS)

# PRESSURE


def convert_p(_p): return _p/p_u


unit_list = [
    'a: '+str(a_binary),
    'R_s: '+str(R_s),
    'R_NS: '+str(R_NS),
    'R_w_s: '+str(R_w_s),
    'R_w_NS: '+str(R_w_NS),
    'm_star: '+str(m_star),
    'm_NS: '+str(m_NS),
    'Mdot_s: '+str(Mdot_s),
    'v_w_s: '+str(v_w_s),
    'L_w_NS: '+str(L_w_NS),
    'B_star: '+str(B_star),
    'B_NS: '+str(B_NS)
]

print(unit_list)
