//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
// C headers
// C++ headers
#include <math.h>

#include <algorithm>
#include <atomic>
#include <cmath>
#include <cstdio>   // fopen(), fprintf(), freopen()
#include <cstring>  // strcmp()
#include <sstream>
#include <stdexcept>
#include <string>
// Athena++ headers
#include "../athena.hpp"
#include "../athena_arrays.hpp"
#include "../coordinates/coordinates.hpp"
#include "../eos/eos.hpp"
#include "../field/field.hpp"
#include "../globals.hpp"
#include "../hydro/hydro.hpp"
#include "../mesh/mesh.hpp"
#include "../parameter_input.hpp"

void BinarySrc(MeshBlock *pmb, const Real time, const Real dt, const AthenaArray<Real> &prim,
               const AthenaArray<Real> &prim_scalar, const AthenaArray<Real> &bcc, AthenaArray<Real> &cons,
               AthenaArray<Real> &cons_scalar);

double gamma_eos = 5.0 / 3;
double rho_amb = 1e-10;
double cs_amb = 1e-13;
double p_amb = 0.0;

// orbital parameters
double m_star = 0.0;
double m_NS = 0.0;
double a_binary = 0.0;
double e_binary = 0.0;

// solar wind parameters
double rho_w_s = 0;
double R_w_s = 0.0;
double v_w_s = 0;
double p_w_s = 0.0;
double cs_w_s = 0.0;
double B_star = 0.0;
double Mdot_s = 0.0;

// pulsar wind parameters
bool pulsar_wind = false;
double rho_w_NS = 0;
double R_w_NS = 0.0;
double v_w_NS = 0.0;
double p_w_NS = 0.0;
double cs_w_NS = 0.0;
double B_NS = 0;
double sigma_NS = 0;

double pi = 3.14159265358979311599796346854;

template <typename Fun>
auto root_bisection(Fun f, decltype(f(0)) low, decltype(f(0)) high) -> decltype(f(0)) {
    using Scalar = decltype(f(0));

    for (; fabs((high - low)) > fabs(high) * 1e-6;) {
        Scalar mid = 0.5 * (high + low);
        if (f(mid) > 0)
            high = mid;
        else
            low = mid;
    }
    return 0.5 * (high + low);
}

double mean_anomaly2true_anomaly(double M_anomaly, double e) {
    double E_anomaly = 0;

    if (fabs(M_anomaly) > 1e-6) {
        E_anomaly = root_bisection([=](double x) -> double { return (x - e * sin(x) - M_anomaly); }, 0, 2 * pi);
    }
    return 2 * atan2(sqrt(1 + e) * sin(E_anomaly * 0.5), sqrt(1 - e) * cos(0.5 * E_anomaly));
}

double time2true_anomaly(double t, double e, double a, double m_tot) {
    double M = sqrt(m_tot / a / a / a) * t;
    M = (M / (2 * pi) - floor(M / (2 * pi))) * 2 * pi;
    return mean_anomaly2true_anomaly(M, e);
}

double elliptic_orbit_r(double mu, double e, double a) { return a * (1 - e * e) / (1 + e * cos(mu)); }

double elliptic_orbit_v(double m_tot, double e, double a) { return sqrt(m_tot / (a * (1 - e * e))); }

double v_to_gamma(double v) { return 1 / sqrt(1 - v * v / 3e10 / 3e10); }

double gamma_to_v(double gamma) { return 3e10 * sqrt(1 - (1 / gamma) * (1 / gamma)); }

static Real A1(const Real x1, const Real x2, const Real x3, double x_offset, double y_offset, double R, double B);

static Real A2(const Real x1, const Real x2, const Real x3, double x_offset, double y_offset, double R, double B);

static Real A3(const Real x1, const Real x2, const Real x3);

void SetBfield(MeshBlock *pmb, Coordinates *pco, FaceField &b, int is, int ie, int js, int je, int ks, int ke,
               double x_offset, double y_offset, double R, double B, int ngh);

double print_wind_info(std::string name, double ei, double ek, double eB, double p_g, double S, double v, double rho,
                       double cs, double B, double R_w) {
    double p_mag = eB;
    double F = S * v;
    double Lt = (ei + p_g) * F;
    double Lk = ek * F;
    double Lm = (eB + p_mag) * F;
    double Ltot = Lt + Lk + Lm;
    std::cout << name << " :\n";
    std::cout << " -Launch Radius r_w: " << R_w << "\n";
    std::cout << "  --Density(r_w): " << rho << std::endl;
    std::cout << "  --B(r_w): " << B << std::endl;
    std::cout << "  --Sound speed(r_w): " << cs << std::endl;
    std::cout << "  --Velocity(r_w): " << v << std::endl;
    std::cout << " -Total Luminosity: " << Ltot << std::endl;
    std::cout << "   --Thermal Luminosity: " << Lt << ' ' << 100 * Lt / Ltot << "%" << std::endl;
    std::cout << "   --Kinetic Luminosity: " << Lk << ' ' << 100 * Lk / Ltot << "%" << std::endl;
    std::cout << "   --Mangetic Luminosity: " << Lm << ' ' << 100 * Lm / Ltot << "%" << std::endl;
    std::cout << "   --sigma: " << Lm / (Lt + Lk) << std::endl;
    std::cout << std::endl;

    return Ltot;
}

void Mesh::InitUserMeshData(ParameterInput *pin) {
    // read hydro parameters from the input file
    gamma_eos = pin->GetReal("hydro", "gamma");
    rho_amb = pin->GetReal("problem", "rho_amb");  // density of the ambient medium
    cs_amb = pin->GetReal("problem", "cs_amb");    // pressure of the ambient medium
    p_amb = cs_amb * cs_amb * rho_amb / gamma_eos;

    // read binary orbit parameters from the input file
    a_binary = pin->GetReal("problem", "a_binary");  // semi-major axis of the binary system
    e_binary = pin->GetReal("problem", "e_binary");  // eccentricity of the binary system
    m_star = pin->GetReal("problem", "m_star");      // mass of the star
    m_NS = pin->GetReal("problem", "m_NS");          // mass of the NS

    // read wind parameters from the input file

    // solar wind parameters

    R_w_s = pin->GetReal("problem", "R_w_s");    // winding launch radius
    v_w_s = pin->GetReal("problem", "v_w_s");    // velocity of the wind at wind launch surface
    cs_w_s = pin->GetReal("problem", "cs_w_s");  // sound speed of the wind at the wind launch surface
    Mdot_s = pin->GetReal("problem", "Mdot_s");  //
    B_star = pin->GetReal("problem", "B_star");  // magnetic field of the star at the wind launch surface R_w_s

    // pulsar wind parameters
    pulsar_wind = pin->GetBoolean("problem", "pulsar_wind");  // if the pulsar wind is enabled
    R_w_NS = pin->GetReal("problem", "R_w_NS");               // radius of the NS wind launch surface
    v_w_NS = pin->GetReal("problem", "v_w_NS");               // velocity of the wind from the NS
    cs_w_NS = pin->GetReal("problem", "cs_w_NS");             // sound speed of the wind from the NS
    sigma_NS = pin->GetReal("problem", "sigma_NS");
    B_NS = pin->GetReal("problem", "B_NS");  // magnetic field of the NS at the wind launch surface R_w_NS

    // calculate solar wind parameters
    rho_w_s = Mdot_s / (4 * PI * R_w_s * R_w_s * v_w_s);  //
    p_w_s = cs_w_s * cs_w_s * rho_w_s / gamma_eos;

    // calculate the pulsar wind parameters

    // sigma = (B^2/2 + P_mag)/(e_i+e_k+p_g) = B^2/(p_g/(gamma_eos-1) + 0.5*rho*v^2 + p_g)
    rho_w_NS = B_NS * B_NS / sigma_NS / (cs_w_NS * cs_w_NS / (gamma_eos - 1) + 0.5 * v_w_NS * v_w_NS);
    p_w_NS = cs_w_NS * cs_w_NS * rho_w_NS / gamma_eos;

    // calculate & print the wind info
    std::cout << "Printing info in code unit\n";
    double v_orb = elliptic_orbit_v(m_star + m_NS, e_binary, a_binary);
    std::cout << " -Orbital velocity: " << v_orb << "\n";
    std::cout << " -Ambient density:  " << rho_amb << "\n";
    std::cout << " -Ambient sound speed: " << cs_amb << "\n\n";

    double ei_w_s = p_w_s / (gamma_eos - 1);        // internal energy density of the wind from the star
    double ek_w_s = 0.5 * rho_w_s * v_w_s * v_w_s;  // kinetic energy density of the wind from the star
    double eB_w_s = B_star * B_star / 2;            // magnetic energy density of the wind from the star
    double S_w_s = 4 * pi * R_w_s * R_w_s;          //  launch surface area of the wind from the star

    double ei_w_NS = p_w_NS / (gamma_eos - 1);          // internal energy of the wind from the NS
    double ek_w_NS = 0.5 * rho_w_NS * v_w_NS * v_w_NS;  // kinetic energy of the wind from the NS
    double eB_w_NS = B_NS * B_NS / 2;                   // magnetic energy of the wind from the NS
    double S_w_NS = 4 * pi * R_w_NS * R_w_NS;           //  launch surface area of the wind from the NS

    double L_s =
        print_wind_info("Stellar Wind", ei_w_s, ek_w_s, eB_w_s, p_w_s, S_w_s, v_w_s, rho_w_s, cs_w_s, B_star, R_w_s);
    if (pulsar_wind) {
        double L_NS = print_wind_info("Pulsar Wind", ei_w_NS, ek_w_NS, eB_w_NS, p_w_NS, S_w_NS, v_w_NS, rho_w_NS,
                                      cs_w_NS, B_NS, R_w_NS);
        std::cout << " -Wind Luminosity Ratio (pulsar/stellar): " << L_NS / L_s << std::endl;
    }
    // set the source function
    EnrollUserExplicitSourceFunction(BinarySrc);
    return;
}

void get_orbit_info(double time, double &x_s, double &y_s, double &vx_s, double &vy_s, double &x_NS, double &y_NS,
                    double &vx_NS, double &vy_NS) {
    const double m_tot = m_star + m_NS;
    // calculate the current position and velocity of the star and the NS
    double mu = time2true_anomaly(time, e_binary, a_binary, m_tot);
    double r_orb = elliptic_orbit_r(mu, e_binary, a_binary);
    double v_orb = elliptic_orbit_v(m_tot, e_binary, a_binary);
    double rx = r_orb * cos(mu);
    double ry = r_orb * sin(mu);
    double vx = -v_orb * sin(mu);
    double vy = v_orb * (cos(mu) + e_binary);

    x_s = m_NS / m_tot * rx;
    y_s = m_NS / m_tot * ry;
    vx_s = m_NS / m_tot * vx;
    vy_s = m_NS / m_tot * vy;

    x_NS = -m_star / m_tot * rx;
    y_NS = -m_star / m_tot * ry;
    vx_NS = -m_star / m_tot * vx;
    vy_NS = -m_star / m_tot * vy;
}

void UpdateSrc(MeshBlock *pmb, const Real time, const Real dt, AthenaArray<Real> &cons, const AthenaArray<Real> &bcc) {
    const double eps_s = R_w_s / 10;    // softening length
    const double eps_NS = R_w_NS / 10;  // softening length
    double x_s = 0, y_s = 0, vx_s = 0, vy_s = 0, x_NS = 0, y_NS = 0, vx_NS = 0, vy_NS = 0;
    get_orbit_info(time, x_s, y_s, vx_s, vy_s, x_NS, y_NS, vx_NS, vy_NS);

    for (int k = pmb->ks; k <= pmb->ke; ++k) {
        for (int j = pmb->js; j <= pmb->je; ++j) {
            for (int i = pmb->is; i <= pmb->ie; ++i) {
                Real x = pmb->pcoord->x1v(i);
                Real y = pmb->pcoord->x2v(j);
                Real z = pmb->pcoord->x3v(k);

                double r2s = sqrt((x - x_s) * (x - x_s) + (y - y_s) * (y - y_s) + z * z + eps_s * eps_s);
                double r2NS = sqrt((x - x_NS) * (x - x_NS) + (y - y_NS) * (y - y_NS) + z * z + eps_NS * eps_NS);

                if (r2s < R_w_s) {
                    double v_wind = 0;      // core region
                    if (r2s > R_w_s / 2) {  // winding region
                        v_wind = v_w_s;
                    }
                    double vx_ = (v_wind * (x - x_s) / r2s + vx_s);
                    double vy_ = (v_wind * (y - y_s) / r2s + vy_s);
                    double vz_ = (v_wind * z / r2s);

                    cons(IDN, k, j, i) = rho_w_s;        // density
                    cons(IM1, k, j, i) = rho_w_s * vx_;  // momentum x
                    cons(IM2, k, j, i) = rho_w_s * vy_;  // momentum y
                    cons(IM3, k, j, i) = rho_w_s * vz_;  // momentum z
                    cons(IEN, k, j, i) =
                        p_w_s / (gamma_eos - 1) + 0.5 * rho_w_s * (vx_ * vx_ + vy_ * vy_ + vz_ * vz_);  // Energy
                }

                if (pulsar_wind && r2NS < R_w_NS) {
                    double v_wind = 0;        // core region
                    if (r2NS > R_w_NS / 2) {  // winding region
                        v_wind = v_w_NS;
                    }
                    double vx_ = (v_wind * (x - x_NS) / r2NS + vx_NS);
                    double vy_ = (v_wind * (y - y_NS) / r2NS + vy_NS);
                    double vz_ = (v_wind * z / r2NS);

                    cons(IDN, k, j, i) = rho_w_NS;        // density
                    cons(IM1, k, j, i) = rho_w_NS * vx_;  // momentum x
                    cons(IM2, k, j, i) = rho_w_NS * vy_;  // momentum y
                    cons(IM3, k, j, i) = rho_w_NS * vz_;  // momentum z
                    cons(IEN, k, j, i) =
                        p_w_NS / (gamma_eos - 1) + 0.5 * rho_w_NS * (vx_ * vx_ + vy_ * vy_ + vz_ * vz_);  // Energy
                }

                // gravity of the star
                cons(IM1, k, j, i) -= dt * cons(IDN, k, j, i) * m_star * (x - x_s) / (r2s * r2s * r2s);  // momentum x
                cons(IM2, k, j, i) -= dt * cons(IDN, k, j, i) * m_star * (y - y_s) / (r2s * r2s * r2s);  // momentum y
                cons(IM3, k, j, i) -= dt * cons(IDN, k, j, i) * m_star * z / (r2s * r2s * r2s);

                // gravity of the NS
                cons(IM1, k, j, i) -= dt * cons(IDN, k, j, i) * m_NS * (x - x_NS) / (r2NS * r2NS * r2NS);  // momentum x
                cons(IM2, k, j, i) -= dt * cons(IDN, k, j, i) * m_NS * (y - y_NS) / (r2NS * r2NS * r2NS);  // momentum y
                cons(IM3, k, j, i) -= dt * cons(IDN, k, j, i) * m_NS * z / (r2NS * r2NS * r2NS);

                // std::cout << "Dens: " << cons(IDN, k, j, i) << std::endl;
                // std::cout << "Momentum(x): " << cons(IM1, k, j, i) << std::endl;
            }
        }
    }
    return;
}

void MeshBlock::ProblemGenerator(ParameterInput *pin) {
    // set the initial condition of the ambiant medium
    for (int k = ks; k <= ke; k++) {
        for (int j = js; j <= je; j++) {
            for (int i = is; i <= ie; i++) {
                phydro->u(IDN, k, j, i) = rho_amb;
                phydro->u(IM1, k, j, i) = 0.0;
                phydro->u(IM2, k, j, i) = 0.0;
                phydro->u(IM3, k, j, i) = 0.0;
                phydro->u(IEN, k, j, i) = p_amb / (gamma_eos - 1);
            }
        }
    }

    double x_s = 0, y_s = 0, vx_s = 0, vy_s = 0, x_NS = 0, y_NS = 0, vx_NS = 0, vy_NS = 0;
    get_orbit_info(0, x_s, y_s, vx_s, vy_s, x_NS, y_NS, vx_NS, vy_NS);

    if (MAGNETIC_FIELDS_ENABLED) {
        // Star magnetic field
        SetBfield(this, pcoord, pfield->b, is, ie, js, je, ks, ke, x_s, y_s, R_w_s, B_star, NGHOST);
        //  NS magnetic field
        if (pulsar_wind) {
            SetBfield(this, pcoord, pfield->b, is, ie, js, je, ks, ke, x_NS, y_NS, R_w_NS, B_NS, NGHOST);
        }
        pfield->CalculateCellCenteredField(pfield->b, pfield->bcc, pcoord, is, ie, js, je, ks, ke);
    }

    // set the initial condition of the wind
    UpdateSrc(this, 0, 0, phydro->u, pfield->bcc);

    // add magnetic field contribution to the energy // previously we did this in the source function
    // if you forget to do this, the primitive variables derived from the conserved variables will be wrong
    for (int k = ks; k <= ke; k++) {
        for (int j = js; j <= je; j++) {
            for (int i = is; i <= ie; i++) {
                double Bx = pfield->b.x1f(k, j, i);
                double By = pfield->b.x2f(k, j, i);
                double Bz = pfield->b.x3f(k, j, i);
                phydro->u(IEN, k, j, i) += 0.5 * (Bx * Bx + By * By + Bz * Bz);
            }
        }
    }
}

void BinarySrc(MeshBlock *pmb, const Real time, const Real dt, const AthenaArray<Real> &prim,
               const AthenaArray<Real> &prim_scalar, const AthenaArray<Real> &bcc, AthenaArray<Real> &cons,
               AthenaArray<Real> &cons_scalar) {
    UpdateSrc(pmb, time, dt, cons, bcc);
}

void SetBfield(MeshBlock *pmb, Coordinates *pco, FaceField &b, int is, int ie, int js, int je, int ks, int ke,
               double x_offset, double y_offset, double R, double B, int ngh) {
    BoundaryValues *pbval = pmb->pbval;
    AthenaArray<Real> a1, a2, a3;
    int nx1 = pmb->block_size.nx1 + 2 * ngh + 1;
    int nx2 = pmb->block_size.nx2 + 2 * ngh + 1;
    int nx3 = pmb->block_size.nx3 + 2 * ngh + 1;
    a1.NewAthenaArray(nx3, nx2, nx1);
    a2.NewAthenaArray(nx3, nx2, nx1);
    a3.NewAthenaArray(nx3, nx2, nx1);

    int level = pmb->loc.level;
    for (int k = ks; k <= ke + 1; k++) {
        for (int j = js; j <= je + 1; j++) {
            for (int i = is; i <= ie + 1; i++) {
                a1(k, j, i) = A1(pco->x1v(i), pco->x2f(j), pco->x3f(k), x_offset, y_offset, R, B);
                a2(k, j, i) = A2(pco->x1f(i), pco->x2v(j), pco->x3f(k), x_offset, y_offset, R, B);
                a3(k, j, i) = A3(pco->x1f(i), pco->x2f(j), pco->x3v(k));
            }
        }
    }

    AthenaArray<Real> area, len, len_p1;
    area.NewAthenaArray(nx1);
    len.NewAthenaArray(nx1);
    len_p1.NewAthenaArray(nx1);
    Real real_min = std::numeric_limits<Real>::min();

    // for 1,2,3-D
    for (int k = ks; k <= ke; ++k) {
        int jl = js;
        int ju = je + 1;
        for (int j = jl; j <= ju; ++j) {
            pco->Face2Area(k, j, is, ie, area);
            pco->Edge3Length(k, j, is, ie + 1, len);
            for (int i = is; i <= ie; ++i) {
                area(i) = (area(i) < real_min) ? real_min : area(i);
                b.x2f(k, j, i) += -1.0 * (len(i + 1) * a3(k, j, i + 1) - len(i) * a3(k, j, i)) / area(i);
            }
        }
    }
    for (int k = ks; k <= ke + 1; ++k) {
        for (int j = js; j <= je; ++j) {
            pco->Face3Area(k, j, is, ie, area);
            pco->Edge2Length(k, j, is, ie + 1, len);
            for (int i = is; i <= ie; ++i) {
                b.x3f(k, j, i) += (len(i + 1) * a2(k, j, i + 1) - len(i) * a2(k, j, i)) / area(i);
            }
        }
    }

    // for 2D and 3D
    if (pmb->block_size.nx2 > 1) {
        for (int k = ks; k <= ke; ++k) {
            for (int j = js; j <= je; ++j) {
                pco->Face1Area(k, j, is, ie + 1, area);
                pco->Edge3Length(k, j, is, ie + 1, len);
                pco->Edge3Length(k, j + 1, is, ie + 1, len_p1);
                for (int i = is; i <= ie + 1; ++i) {
                    b.x1f(k, j, i) += (len_p1(i) * a3(k, j + 1, i) - len(i) * a3(k, j, i)) / area(i);
                }
            }
        }
        for (int k = ks; k <= ke + 1; ++k) {
            for (int j = js; j <= je; ++j) {
                pco->Face3Area(k, j, is, ie, area);
                pco->Edge1Length(k, j, is, ie, len);
                pco->Edge1Length(k, j + 1, is, ie, len_p1);
                for (int i = is; i <= ie; ++i) {
                    b.x3f(k, j, i) -= (len_p1(i) * a1(k, j + 1, i) - len(i) * a1(k, j, i)) / area(i);
                }
            }
        }
    }
    // for 3D only
    if (pmb->block_size.nx3 > 1) {
        for (int k = ks; k <= ke; ++k) {
            for (int j = js; j <= je; ++j) {
                pco->Face1Area(k, j, is, ie + 1, area);
                pco->Edge2Length(k, j, is, ie + 1, len);
                pco->Edge2Length(k + 1, j, is, ie + 1, len_p1);
                for (int i = is; i <= ie + 1; ++i) {
                    b.x1f(k, j, i) -= (len_p1(i) * a2(k + 1, j, i) - len(i) * a2(k, j, i)) / area(i);
                }
            }
        }
        for (int k = ks; k <= ke; ++k) {
            int jl = js;
            int ju = je + 1;
            for (int j = jl; j <= ju; ++j) {
                pco->Face2Area(k, j, is, ie, area);
                pco->Edge1Length(k, j, is, ie, len);
                pco->Edge1Length(k + 1, j, is, ie, len_p1);
                for (int i = is; i <= ie; ++i) {
                    area(i) = (area(i) < real_min) ? real_min : area(i);
                    b.x2f(k, j, i) += (len_p1(i) * a1(k + 1, j, i) - len(i) * a1(k, j, i)) / area(i);
                }
            }
        }
    }
    return;
}

static Real A1(const Real x1, const Real x2, const Real x3, double x_offset, double y_offset, double R, double B) {
    Real a1 = 0.0;
    double x1c = x1 - x_offset;
    double x2c = x2 - y_offset;
    double x3c = x3 - 0;
    Real rc = std::max(sqrt(x1c * x1c + x2c * x2c + x3c * x3c), R / 2.);
    a1 = B * R * R / rc / rc / rc * (-1. * x2c);
    return (a1);
}

static Real A2(const Real x1, const Real x2, const Real x3, double x_offset, double y_offset, double R, double B) {
    Real a2 = 0.0;
    double x1c = x1 - x_offset;
    double x2c = x2 - y_offset;
    double x3c = x3 - 0;
    Real rc = std::max(sqrt(x1c * x1c + x2c * x2c + x3c * x3c), R / 2.);
    a2 = B * R * R / rc / rc / rc * (1. * x1c);
    return (a2);
}

static Real A3(const Real x1, const Real x2, const Real x3) {
    Real a3 = 0.0;
    return (a3);
}