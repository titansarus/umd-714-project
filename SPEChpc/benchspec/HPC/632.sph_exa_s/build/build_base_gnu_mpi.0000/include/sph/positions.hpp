#pragma once

#include <vector>
#include <cmath>

namespace sphexa
{
namespace sph
{
template <typename T, class Dataset>
void computePositions(const std::vector<int> &l, Dataset &d)
{
    const int n = l.size();
    const int *clist = l.data();

    const T *grad_P_x = d.grad_P_x.data();
    const T *grad_P_y = d.grad_P_y.data();
    const T *grad_P_z = d.grad_P_z.data();
    const T *dt = d.dt.data();
    const T *du = d.du.data();
    T *x = d.x.data();
    T *y = d.y.data();
    T *z = d.z.data();
    T *vx = d.vx.data();
    T *vy = d.vy.data();
    T *vz = d.vz.data();
    T *x_m1 = d.x_m1.data();
    T *y_m1 = d.y_m1.data();
    T *z_m1 = d.z_m1.data();
    T *u = d.u.data();
    T *du_m1 = d.du_m1.data();
    T *dt_m1 = d.dt_m1.data();

    const BBox<T> bbox = d.bbox;
#if defined(SPEC_OPENMP) || defined(SPEC_OPENMP_TARGET)
#pragma omp parallel for
#endif
    for (int pi = 0; pi < n; pi++)
    {
        int i = clist[pi];

        // ADD COMPONENT DUE TO THE GRAVITY HERE
        T ax = -(grad_P_x[i]); //-G * fx
        T ay = -(grad_P_y[i]); //-G * fy
        T az = -(grad_P_z[i]); //-G * fz

#ifndef NDEBUG
        if (std::isnan(ax) || std::isnan(ay) || std::isnan(az))
            printf("ERROR::UpdateQuantities(%d) acceleration: (%f %f %f)\n", i, ax, ay, az);
#endif

        // Update positions according to Press (2nd order)
        T deltaA = dt[i] + 0.5 * dt_m1[i];
        T deltaB = 0.5 * (dt[i] + dt_m1[i]);

        T valx = (x[i] - x_m1[i]) / dt_m1[i];
        T valy = (y[i] - y_m1[i]) / dt_m1[i];
        T valz = (z[i] - z_m1[i]) / dt_m1[i];

        vx[i] = valx + ax * deltaA;
        vy[i] = valy + ay * deltaA;
        vz[i] = valz + az * deltaA;

        x_m1[i] = x[i];
        y_m1[i] = y[i];
        z_m1[i] = z[i];

        // x[i] = x + dt[i] * valx + ax * dt[i] * deltaB;
        x[i] += dt[i] * valx + (vx[i] - valx) * dt[i] * deltaB / deltaA;
        y[i] += dt[i] * valy + (vy[i] - valy) * dt[i] * deltaB / deltaA;
        z[i] += dt[i] * valz + (vz[i] - valz) * dt[i] * deltaB / deltaA;

        if (bbox.PBCx && x[i] < bbox.xmin)
            x[i] += (bbox.xmax - bbox.xmin);
        else if (bbox.PBCx && x[i] > bbox.xmax)
            x[i] -= (bbox.xmax - bbox.xmin);
        if (bbox.PBCy && y[i] < bbox.ymin)
            y[i] += (bbox.ymax - bbox.ymin);
        else if (bbox.PBCy && y[i] > bbox.ymax)
            y[i] -= (bbox.ymax - bbox.ymin);
        if (bbox.PBCz && z[i] < bbox.zmin)
            z[i] += (bbox.zmax - bbox.zmin);
        else if (bbox.PBCz && z[i] > bbox.zmax)
            z[i] -= (bbox.zmax - bbox.zmin);

        // Update the energy according to Adams-Bashforth (2nd order)
        deltaA = 0.5 * dt[i] * dt[i] / dt_m1[i];
        deltaB = dt[i] + deltaA;

        u[i] += du[i] * deltaB - du_m1[i] * deltaA;

#ifndef NDEBUG
        if (std::isnan(u[i]))
            printf("ERROR::UpdateQuantities(%d) internal energy: u %f du %f dB %f du_m1 %f dA %f\n", i, u[i], du[i], deltaB, du_m1[i],
                   deltaA);
#endif

        du_m1[i] = du[i];
        dt_m1[i] = dt[i];
    }
}
} // namespace sph
} // namespace sphexa
