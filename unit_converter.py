_l = 7e10 # [cm]
_m = 2e33 # [g]
_G = 6.67e-8 # [cm^3/g/s^2]

l = 1
G = 1 

t = 1 # s (2e-8)

m = l**3/t**2/_G # g
dm = m/t # g/s
v = l/t # cm/s
p = m/l/t**2 # g/cm/s^2
B = l**-0.5 * m**0.5 * t**-1 # cm^-1/2 g^-1/2 s^-1
L = m * l**2 / t**3 # g cm^2 / s^3

print(l**3/m/t**2)

unit_list = [
    'time: '+str(t), 
    'length: '+str(l), 
    'mass: '+str(m), 
    'mass loss: '+str(dm),
    'velocity: '+str(v), 
    'pressure: '+str(p), 
    'magnetic field: '+str(B), 
    'luminosity: '+str(L)
    ]
print(unit_list)

# Units to find

# LENGTH [cm]
def convert_l(_l): return _l/l

_a = 3 
a_binary = convert_l(_a)

_R_w_s = 1 # r_sun (wind-launchiung surface)
R_w_s = convert_l(_R_w_s)

_R_w_NS = 0.5
R_w_NS = convert_l(_R_w_NS)

# MASS [g]
def convert_m(_m): return _m/m

_m_star = 3e34
m_star = convert_m(_m_star)
_m_NS = 2.8e33
m_NS = convert_m(_m_NS)

# MASS LOSS
def convert_dm(_dm): return _dm/dm

_Mdot_s = 6.3e19
Mdot_s = convert_dm(_Mdot_s)

# VELOCITY
def convert_v(_v): return _v/v
_v_w_s = 3e8
v_w_s = convert_v(_v_w_s)  # stellar wind speed at the surface of the star

# LUMINOSITY
def convert_L(_L): return _L/L
_L_w_NS = 1e34
L_w_NS = convert_L(_L_w_NS)

# MAGNETIC FIELD
def convert_B(_B): return _B/B
_B_star = 10
B_star = convert_B(_B_star)
_B_NS = 1e14
B_NS = convert_B(_B_NS)

# PRESSURE
def convert_p(_p): return _p/p
_p_w_s = 1e-5 
p_w_s = convert_p(_p_w_s)
_p_w_NS = 1e-5 
p_w_NS = convert_p(_p_w_NS)

unit_list = [
    'a: '+str(a_binary), 
    'R_w_s: '+str(R_w_s), 
    'R_w_NS: '+str(R_w_NS), 
    'm_star: '+str(m_star),
    'm_NS: '+str(m_NS),
    'Mdot_s: '+str(Mdot_s), 
    'v_w_s: '+str(v_w_s),
    'L_w_NS: '+str(L_w_NS), 
    'B_star: '+str(B_star), 
    'B_NS: '+str(B_NS), 
    'p_w_s: '+str(p_w_s),
    'p_w_NS: '+str(p_w_NS)
]

print(unit_list)