namespace stream_cpu
{
    double dx1(unsigned int l, double* arr) {
        return (arr[l + 1] - arr[l - 1]) / (2.0 * config.hx);
    };

    double dy1(unsigned int l, double* arr) {
        return (arr[l + config.offset] - arr[l - config.offset]) / (2.0 * config.hy);
    };

    double dx2(unsigned int l, double* arr) {
        return (arr[l + 1] - 2.0 * arr[l] + arr[l - 1]) / (config.hx * config.hx);
    };

    double dy2(unsigned int l, double* arr) {
        return (arr[l + config.offset] - 2.0 * arr[l] + arr[l - config.offset]) / (config.hy * config.hy);
    };

    double dx1_eq_0_forward(unsigned int l, double* f) {
        return (4.0 * f[l + 1] - f[l + 2]) / 3.0;
    }
    double dx1_eq_0_back(unsigned int l, double* f) {
        return (4.0 * f[l - 1] - f[l - 2]) / 3.0;
    }
    double dy1_eq_0_up(unsigned int l, double* f) {
        return (4.0 * f[l + config.offset] - f[l + 2 * config.offset]) / 3.0;
    }
    double dy1_eq_0_down(unsigned int l, double* f) {
        return (4.0 * f[l - config.offset] - f[l - 2 * config.offset]) / 3.0;
    }

    double dx1_forward(unsigned int l, double *f) {
        return -0.5 * (3.0 * f[l] - 4.0 * f[l + 1] + f[l + 2]) / config.hx;
    }
    double dx1_back(unsigned int l, double *f) {
        return  0.5 * (3.0 * f[l] - 4.0 * f[l - 1] + f[l - 2]) / config.hx;
    }
    double dy1_up(unsigned int l, double *f) {
        return  -0.5 * (3.0 * f[l] - 4.0 * f[l + config.offset] + f[l + 2 * config.offset]) / config.hx;
    }
    double dy1_down(unsigned int l, double *f) {
        return  0.5 * (3.0 * f[l] - 4.0 * f[l - config.offset] + f[l - 2 * config.offset]) / config.hx;
    }


    #define VX_ dy1(l, ksi)
    #define VY_ -dx1(l, ksi)

    void vorticity(double* omega_new, double* omega, double* ksi, double* T, double* C)
    {
        auto InnerComputing = [&](unsigned int l) {
            return omega[l] + config.tau * (
                (dx1(l, ksi) * dy1(l, omega) - dy1(l, ksi) * dx1(l, omega)) //nonlinear term
                + (dx2(l, omega) + dy2(l, omega)) /** config.Pr*/

                + config.grav_y * config.Ra / config.Pr * (dx1(l, T) - config.density_x)
                - config.grav_x * config.Ra / config.Pr * (dy1(l, T) - config.density_y)

                + config.grav_y * config.Ra / config.Pr * config.K * (dx1(l, C) - config.density_x)
                - config.grav_x * config.Ra / config.Pr * config.K * (dy1(l, C) - config.density_y)
                );
        };

        for (unsigned int j = 0; j <= config.ny; ++j) {
            for (unsigned int i = 0; i <= config.nx; ++i) {
                unsigned int l = i + config.offset * j;

                if (l < config.N) {
                    /* INNER */
                    if (i > 0 && i < config.nx && j > 0 && j < config.ny) {
                        omega_new[l] = InnerComputing(l);
                    }
                    else {
                        if (j == 0 && (i > 0 && i < config.nx)) {
                            omega_new[l] = -0.5 / (config.hy * config.hy) * (8.0 * ksi[l + config.offset] - ksi[l + config.offset * 2]);
                            continue;
                        }
                        else if (j == config.ny && (i > 0 && i < config.nx)) {
                            omega_new[l] = -0.5 / (config.hy * config.hy) * (8.0 * ksi[l - config.offset] - ksi[l - config.offset * 2]);
                            continue;
                        }

                        if (config.xbc == 0) { // closed
                            if (i == 0 && (j > 0 && j < config.ny))
                                omega_new[l] = -0.5 / (config.hx * config.hx) * (8.0 * ksi[l + 1] - ksi[l + 2]);
                            if (i == config.nx && (j > 0 && j < config.ny))
                                omega_new[l] = -0.5 / (config.hx * config.hx) * (8.0 * ksi[l - 1] - ksi[l - 2]);
                            continue;
                        }
                        else if (config.xbc == 1) { // periodic
                            if (i == 0 && (j > 0 && j < config.ny)) {
                                unsigned int ll = config.nx - 1 + config.offset * j;
                                omega_new[l] = InnerComputing(ll);
                                continue;
                            }
                            if (i == config.nx && (j > 0 && j < config.ny)) {
                                unsigned int ll = 1 + config.offset * j;
                                omega_new[l] = InnerComputing(ll);
                                continue;
                            }
                        }
                        {
                            omega_new[l] = 0;
                            omega_new[l] = 0;
                        }
                    }
                }
            }
        }
    }

    void vorticity_Soret(double *omega_new, double *omega, double *ksi, double *T, double *C)
    {
        auto InnerComputing = [&](unsigned int l) {
            return omega[l] + config.tau * (
                (dx1(l, ksi) * dy1(l, omega) - dy1(l, ksi) * dx1(l, omega)) //nonlinear term
                + (dx2(l, omega) + dy2(l, omega)) /** config.Pr*/

                + config.grav_y * config.Ra / config.Pr * (dx1(l, T) /*- config.density_x*/)
                - config.grav_x * config.Ra / config.Pr * (dy1(l, T) /*- config.density_y*/)

                + config.grav_y * config.Ra / config.Pr * (dx1(l, C) /*- config.density_x*/)
                - config.grav_x * config.Ra / config.Pr * (dy1(l, C) /*- config.density_y*/)
                );
        };

        for (unsigned int j = 0; j <= config.ny; ++j) {
            for (unsigned int i = 0; i <= config.nx; ++i) {
                unsigned int l = i + config.offset * j;

                if (l < config.N) {
                    /* INNER */
                    if (i > 0 && i < config.nx && j > 0 && j < config.ny) {
                        omega_new[l] = InnerComputing(l);
                    }
                    else {
                        if (j == 0 && (i > 0 && i < config.nx)) {
                            omega_new[l] = -0.5 / (config.hy * config.hy) * (8.0 * ksi[l + config.offset] - ksi[l + config.offset * 2]);
                            continue;
                        }
                        else if (j == config.ny && (i > 0 && i < config.nx)) {
                            omega_new[l] = -0.5 / (config.hy * config.hy) * (8.0 * ksi[l - config.offset] - ksi[l - config.offset * 2]);
                            continue;
                        }

                        if (config.xbc == 0) { // closed
                            if (i == 0 && (j > 0 && j < config.ny))
                                omega_new[l] = -0.5 / (config.hx * config.hx) * (8.0 * ksi[l + 1] - ksi[l + 2]);
                            if (i == config.nx && (j > 0 && j < config.ny))
                                omega_new[l] = -0.5 / (config.hx * config.hx) * (8.0 * ksi[l - 1] - ksi[l - 2]);
                            continue;
                        }
                        else if (config.xbc == 1) { // periodic
                            if (i == 0 && (j > 0 && j < config.ny)) {
                                unsigned int ll = config.nx - 1 + config.offset * j;
                                omega_new[l] = InnerComputing(ll);
                                continue;
                            }
                            if (i == config.nx && (j > 0 && j < config.ny)) {
                                unsigned int ll = 1 + config.offset * j;
                                omega_new[l] = InnerComputing(ll);
                                continue;
                            }
                        }
                        {
                            omega_new[l] = 0;
                            omega_new[l] = 0;
                        }
                    }
                }
            }
        }
    }

    void vorticity_quadratic_temperature(double *omega_new, double *omega, double *ksi, double *T)
    {
        auto InnerComputing = [&](unsigned int l) {
            return omega[l] + config.tau * (
                (dx1(l, ksi) * dy1(l, omega) - dy1(l, ksi) * dx1(l, omega)) //nonlinear term
                + (dx2(l, omega) + dy2(l, omega))

                + config.grav_y * config.Ra / config.Pr * (dx1(l, T)) * T[l]
                - config.grav_x * config.Ra / config.Pr * (dy1(l, T)) * T[l]
                );
        };

        for (unsigned int j = 0; j <= config.ny; ++j) {
            for (unsigned int i = 0; i <= config.nx; ++i) {
                unsigned int l = i + config.offset * j;

                if (l < config.N) {
                    /* INNER */
                    if (i > 0 && i < config.nx && j > 0 && j < config.ny) {
                        omega_new[l] = InnerComputing(l);
                    }
                    else {
                        if (j == 0 && (i > 0 && i < config.nx)) {
                            omega_new[l] = -0.5 / (config.hy * config.hy) * (8.0 * ksi[l + config.offset] - ksi[l + config.offset * 2]);
                            continue;
                        }
                        else if (j == config.ny && (i > 0 && i < config.nx)) {
                            omega_new[l] = -0.5 / (config.hy * config.hy) * (8.0 * ksi[l - config.offset] - ksi[l - config.offset * 2]);
                            continue;
                        }

                        if (config.xbc == 0) { // closed
                            if (i == 0 && (j > 0 && j < config.ny))
                                omega_new[l] = -0.5 / (config.hx * config.hx) * (8.0 * ksi[l + 1] - ksi[l + 2]);
                            if (i == config.nx && (j > 0 && j < config.ny))
                                omega_new[l] = -0.5 / (config.hx * config.hx) * (8.0 * ksi[l - 1] - ksi[l - 2]);
                            continue;
                        }
                        else if (config.xbc == 1) { // periodic
                            if (i == 0 && (j > 0 && j < config.ny)) {
                                unsigned int ll = config.nx - 1 + config.offset * j;
                                omega_new[l] = InnerComputing(ll);
                                continue;
                            }
                            if (i == config.nx && (j > 0 && j < config.ny)) {
                                unsigned int ll = 1 + config.offset * j;
                                omega_new[l] = InnerComputing(ll);
                                continue;
                            }
                        }
                        {
                            omega_new[l] = 0;
                            omega_new[l] = 0;
                        }
                    }
                }
            }
        }
    }


    void poisson_stream(double* ksi_new, double* ksi, double* omega)
    {
        double tau = 0.1 * config.hx * config.hy; // Вычисляем tau

        for (unsigned int j = 0; j <= config.ny; ++j) {
            for (unsigned int i = 0; i <= config.nx; ++i) {
                unsigned int l = i + config.offset * j;

                if (l < config.N) {
                    /* INNER */
                    if (i > 0 && i < config.nx && j > 0 && j < config.ny) {
                        ksi_new[l] = ksi[l] + tau * (dx2(l, ksi) + dy2(l, ksi) + omega[l]);
                    }
                    else {
                        if (j == 0 /*&& (i > 0 && i < config.nx)*/) {
                            ksi_new[l] = 0.0;
                            continue;
                        }
                        else if (j == config.ny /*&& (i > 0 && i < config.nx)*/) {
                            ksi_new[l] = 0.0;
                            continue;
                        }

                        if (config.xbc == 0) { // closed
                            if (i == 0 && (j > 0 && j < config.ny))
                                ksi_new[l] = 0.0;
                            if (i == config.nx && (j > 0 && j < config.ny))
                                ksi_new[l] = 0.0;
                            continue;
                        }
                        else if (config.xbc == 1) { // periodic
                            if (i == 0 && (j > 0 && j < config.ny)) {
                                unsigned int ll = config.nx - 1 + config.offset * j;
                                ksi_new[l] = ksi[ll] + tau * (dx2(ll, ksi) + dy2(ll, ksi) + omega[ll]);
                                continue;
                            }
                            if (i == config.nx && (j > 0 && j < config.ny)) {
                                unsigned int ll = 1 + config.offset * j;
                                ksi_new[l] = ksi[ll] + tau * (dx2(ll, ksi) + dy2(ll, ksi) + omega[ll]);
                                continue;
                            }
                        }
                        {
                            ksi_new[l] = 0;
                        }
                    }
                }
            }
        }
    }

    void poisson_stream_v2(double* ksi_new, double* ksi, double* omega)
    {
        double tau = 0.25 * config.hx * config.hy; // Вычисляем tau

        for (unsigned int j = 0; j <= config.ny; ++j) {
            for (unsigned int i = 0; i <= config.nx; ++i) {
                unsigned int l = i + config.offset * j;

                auto inner = [&](unsigned int l)
                {
                    return 0.25 * (ksi[l - 1] + ksi[l + 1] + ksi[l + config.offset] + ksi[l - config.offset]) +  0.25 * config.hx * config.hy * omega[l];
                };
                //if (l < config.N) 
                {
                    /* INNER */
                    if (i > 0 && i < config.nx && j > 0 && j < config.ny) {
                        ksi_new[l] = inner(l);
                    }
                    else {
                        if (j == 0 /*&& (i > 0 && i < config.nx)*/) {
                            ksi_new[l] = 0.0;
                            continue;
                        }
                        else if (j == config.ny /*&& (i > 0 && i < config.nx)*/) {
                            ksi_new[l] = 0.0;
                            continue;
                        }

                        if (config.xbc == 0) { // closed
                            if (i == 0 && (j > 0 && j < config.ny))
                                ksi_new[l] = 0.0;
                            if (i == config.nx && (j > 0 && j < config.ny))
                                ksi_new[l] = 0.0;
                            continue;
                        }
                        else if (config.xbc == 1) { // periodic
                            if (i == 0 && (j > 0 && j < config.ny)) {
                                unsigned int ll = config.nx - 1 + config.offset * j;
                                ksi_new[l] = inner(ll);
                                continue;
                            }
                            if (i == config.nx && (j > 0 && j < config.ny)) {
                                unsigned int ll = 1 + config.offset * j;
                                ksi_new[l] = inner(ll);
                                continue;
                            }
                        }
                        {
                            ksi_new[l] = 0;
                        }
                    }
                }
            }
        }
    }


    void temperature_2d(double* T, double* T0, double* ksi)
    {
        for (unsigned int j = 0; j <= config.ny; ++j) {
            for (unsigned int i = 0; i <= config.nx; ++i) {
                unsigned int l = i + config.offset * j;

                if (l < config.N) {
                    /* INNER */
                    if (i > 0 && i < config.nx && j > 0 && j < config.ny) {
                        T[l] = T0[l]
                            + config.tau * (
                                -dy1(l, ksi) * dx1(l, T0) + dx1(l, ksi) * dy1(l, T0)
                                + (VX_ * config.density_x + VY_ * config.density_y)
                                + (dx2(l, T0) + dy2(l, T0)) / config.Pr
                                );
                        continue;
                    }
                    else {
                        if (j == 0) {
                            T[l] = 0.0;
                            continue;
                        }
                        else if (j == config.ny) {
                            T[l] = 0.0;
                            continue;
                        }

                        if (config.xbc == 0) { // closed
                            if (i == 0 && (j > 0 && j < config.ny)) {
                                T[l] = dx1_eq_0_forward(l, T0);
                                continue;
                            }
                            if (i == config.nx && (j > 0 && j < config.ny)) {
                                T[l] = dx1_eq_0_back(l, T0);
                                continue;
                            }
                        }
                        else if (config.xbc == 1) { // periodic
                            if (i == 0 && (j > 0 && j < config.ny)) {
                                unsigned int ll = config.nx - 1 + config.offset * j;
                                T[l] = T0[ll];
                                continue;
                            }
                            if (i == config.nx && (j > 0 && j < config.ny)) {
                                unsigned int ll = 1 + config.offset * j;
                                T[l] = T0[ll];
                                continue;
                            }
                        }

                        T[l] = 0;
                    }
                }
            }
        }
    }
    void temperature_2d_full(double* T, double* T0, double* ksi)
    {
        for (unsigned int j = 0; j <= config.ny; ++j) {
            for (unsigned int i = 0; i <= config.nx; ++i) {
                unsigned int l = i + config.offset * j;

                if (l < config.N) {
                    /* INNER */
                    if (i > 0 && i < config.nx && j > 0 && j < config.ny) {
                        T[l] = T0[l]
                            + config.tau * (
                                -dy1(l, ksi) * dx1(l, T0) + dx1(l, ksi) * dy1(l, T0)
                                + (dx2(l, T0) + dy2(l, T0)) / config.Pr
                                );
                            continue;
                    }
                    else {
                        if (j == 0) {
                            T[l] = 1.0;
                            continue;
                        }
                        else if (j == config.ny) {
                            T[l] = 0.0;
                            continue;
                        }

                        if (config.xbc == 0) { // closed
                            if (i == 0 && (j > 0 && j < config.ny)) {
                                T[l] = dx1_eq_0_forward(l, T0);
                                continue;
                            }
                            if (i == config.nx && (j > 0 && j < config.ny)) {
                                T[l] = dx1_eq_0_back(l, T0);
                                continue;
                            }
                        }
                        else if (config.xbc == 1) { // periodic
                            if (i == 0 && (j > 0 && j < config.ny)) {
                                unsigned int ll = config.nx - 1 + config.offset * j;
                                T[l] = T0[ll];
                                continue;
                            }
                            if (i == config.nx && (j > 0 && j < config.ny)) {
                                unsigned int ll = 1 + config.offset * j;
                                T[l] = T0[ll];
                                continue;
                            }
                        }

                        T[l] = 0;
                    }
                }
            }
        }
    }


    void temperature_2d_flux(double* T, double* T0, double* ksi)
    {
        for (unsigned int j = 0; j <= config.ny; ++j) {
            for (unsigned int i = 0; i <= config.nx; ++i) {
                unsigned int l = i + config.offset * j;

                if (l < config.N) {
                    /* INNER */
                    if (i > 0 && i < config.nx && j > 0 && j < config.ny) {
                        T[l] = T0[l]
                            + config.tau * (
                                -dy1(l, ksi) * dx1(l, T0) + dx1(l, ksi) * dy1(l, T0)
                                + (VX_ * config.density_x + VY_ * config.density_y)
                                + (dx2(l, T0) + dy2(l, T0)) / config.Pr
                                );
                            continue;
                    }
                    else {
                        if (j == 0) {
                            T[l] = dy1_eq_0_up(l, T0);
                            continue;
                        }
                        else if (j == config.ny) {
                            T[l] = dy1_eq_0_down(l, T0);
                            continue;
                        }

                        if (config.xbc == 0) { // closed
                            if (i == 0 && (j > 0 && j < config.ny)) {
                                T[l] = dx1_eq_0_forward(l, T0);
                                continue;
                            }
                            if (i == config.nx && (j > 0 && j < config.ny)) {
                                T[l] = dx1_eq_0_back(l, T0);
                                continue;
                            }
                        }
                        else if (config.xbc == 1) { // periodic
                            if (i == 0 && (j > 0 && j < config.ny)) {
                                unsigned int ll = config.nx - 1 + config.offset * j;
                                T[l] = T0[ll];
                                continue;
                            }
                            if (i == config.nx && (j > 0 && j < config.ny)) {
                                unsigned int ll = 1 + config.offset * j;
                                T[l] = T0[ll];
                                continue;
                            }
                        }

                        T[l] = 0;
                    }
                }
            }
        }
    }
    void temperature_2d_flux_full(double* T, double* T0, double* ksi)
    {
        for (unsigned int j = 0; j <= config.ny; ++j) {
            for (unsigned int i = 0; i <= config.nx; ++i) {
                unsigned int l = i + config.offset * j;

                if (l < config.N) {
                    /* INNER */
                    if (i > 0 && i < config.nx && j > 0 && j < config.ny) {
                        T[l] = T0[l]
                            + config.tau * (
                                -dy1(l, ksi) * dx1(l, T0) + dx1(l, ksi) * dy1(l, T0)
                                + (dx2(l, T0) + dy2(l, T0)) / config.Pr
                                );
                            continue;
                    }
                    else {
                        if (j == 0) {
                            T[l] = dy1_eq_0_up(l, T0) - (-1 * config.hy * 2.0 / 3.0);
                            continue;
                        }
                        else if (j == config.ny) {
                            T[l] = dy1_eq_0_down(l, T0) + (-1 * config.hy * 2.0 / 3.0);
                            continue;
                        }

                        if (config.xbc == 0) { // closed
                            if (i == 0 && (j > 0 && j < config.ny)) {
                                T[l] = dx1_eq_0_forward(l, T0);
                                continue;
                            }
                            if (i == config.nx && (j > 0 && j < config.ny)) {
                                T[l] = dx1_eq_0_back(l, T0);
                                continue;
                            }
                        }
                        else if (config.xbc == 1) { // periodic
                            if (i == 0 && (j > 0 && j < config.ny)) {
                                unsigned int ll = config.nx - 1 + config.offset * j;
                                T[l] = T0[ll];
                                continue;
                            }
                            if (i == config.nx && (j > 0 && j < config.ny)) {
                                unsigned int ll = 1 + config.offset * j;
                                T[l] = T0[ll];
                                continue;
                            }
                        }

                        T[l] = 0;
                    }
                }
            }
        }
    }


    void concentration_2d(double* C, double* C0, double* ksi)
    {
        for (unsigned int j = 0; j <= config.ny; ++j) {
            for (unsigned int i = 0; i <= config.nx; ++i) {
                unsigned int l = i + config.offset * j;

                if (l < config.N) {
                    /* INNER */
                    if (i > 0 && i < config.nx && j > 0 && j < config.ny) {
                        C[l] = C0[l]
                            + config.tau * (
                                -dy1(l, ksi) * dx1(l, C0) + dx1(l, ksi) * dy1(l, C0)
                                + (VX_ * config.density_x + VY_ * config.density_y)
                                + (dx2(l, C0) + dy2(l, C0)) / (config.Le * config.Pr)
                                );
                        continue;
                    }
                    else {
                        if (j == 0) {
                            C[l] = dy1_eq_0_up(l, C0);
                            continue;
                        }
                        else if (j == config.ny) {
                            C[l] = dy1_eq_0_down(l, C0);
                            continue;
                        }

                        if (config.xbc == 0) { // closed
                            if (i == 0 && (j > 0 && j < config.ny)) {
                                C[l] = dx1_eq_0_forward(l, C0);
                                continue;
                            }
                            if (i == config.nx && (j > 0 && j < config.ny)) {
                                C[l] = dx1_eq_0_back(l, C0);
                                continue;
                            }
                        }
                        else if (config.xbc == 1) { // periodic
                            if (i == 0 && (j > 0 && j < config.ny)) {
                                unsigned int ll = config.nx - 1 + config.offset * j;
                                C[l] = C0[ll];
                                continue;
                            }
                            if (i == config.nx && (j > 0 && j < config.ny)) {
                                unsigned int ll = 1 + config.offset * j;
                                C[l] = C0[ll];
                                continue;
                            }
                        }

                        C[l] = 0;
                    }
                }
            }
        }
    }
    void concentration_2d_full(double* C, double* C0, double* ksi)
    {
        for (unsigned int j = 0; j <= config.ny; ++j) {
            for (unsigned int i = 0; i <= config.nx; ++i) {
                unsigned int l = i + config.offset * j;

                if (l < config.N) {
                    /* INNER */
                    if (i > 0 && i < config.nx && j > 0 && j < config.ny) {
                        C[l] = C0[l]
                            + config.tau * (
                                -dy1(l, ksi) * dx1(l, C0) + dx1(l, ksi) * dy1(l, C0)
                                + (dx2(l, C0) + dy2(l, C0)) / (config.Le * config.Pr)
                                );
                            continue;
                    }
                    else {
                        if (j == 0) {
                            C[l] = dy1_eq_0_up(l, C0) - (-1 * config.hy * 2.0 / 3.0);
                            continue;
                        }
                        else if (j == config.ny) {
                            C[l] = dy1_eq_0_down(l, C0) + (-1 * config.hy * 2.0 / 3.0);
                            continue;
                        }

                        if (config.xbc == 0) { // closed
                            if (i == 0 && (j > 0 && j < config.ny)) {
                                C[l] = dx1_eq_0_forward(l, C0);
                                continue;
                            }
                            if (i == config.nx && (j > 0 && j < config.ny)) {
                                C[l] = dx1_eq_0_back(l, C0);
                                continue;
                            }
                        }
                        else if (config.xbc == 1) { // periodic
                            if (i == 0 && (j > 0 && j < config.ny)) {
                                unsigned int ll = config.nx - 1 + config.offset * j;
                                C[l] = C0[ll];
                                continue;
                            }
                            if (i == config.nx && (j > 0 && j < config.ny)) {
                                unsigned int ll = 1 + config.offset * j;
                                C[l] = C0[ll];
                                continue;
                            }
                        }

                        C[l] = 0;
                    }
                }
            }
        }
    }

    void concentration_2d_full_Soret(double *C, double *C0, double *T0, double *ksi)
    {
        for (unsigned int j = 0; j <= config.ny; ++j) {
            for (unsigned int i = 0; i <= config.nx; ++i) {
                unsigned int l = i + config.offset * j;

                if (l < config.N) {
                    /* INNER */
                    if (i > 0 && i < config.nx && j > 0 && j < config.ny) {
                        C[l] = C0[l]
                            + config.tau * (
                                -dy1(l, ksi) * dx1(l, C0) + dx1(l, ksi) * dy1(l, C0)
                                + ((dx2(l, C0) + dy2(l, C0))
                                + (dx2(l, T0) + dy2(l, T0)) * config.psi) / config.Sc
                                );
                            continue;
                    }
                    else {
                        if (j == 0) {
                            C[l] = dy1_eq_0_up(l, C0) - config.psi * (dy1_up(l, T0) * config.hy * 2.0 / 3.0);
                            continue;
                        }
                        else if (j == config.ny) {
                            C[l] = dy1_eq_0_down(l, C0) + config.psi * (dy1_down(l, T0) * config.hy * 2.0 / 3.0);
                            continue;
                        }

                        if (config.xbc == 0) { // closed
                            if (i == 0 && (j > 0 && j < config.ny)) {
                                C[l] = dx1_eq_0_forward(l, C0) - config.psi * (dx1_forward(l, T0) * config.hx * 2.0 / 3.0);
                                continue;
                            }
                            if (i == config.nx && (j > 0 && j < config.ny)) {
                                C[l] = dx1_eq_0_back(l, C0) + config.psi * (dx1_back(l, T0) * config.hx * 2.0 / 3.0);
                                continue;
                            }
                        }
                        else if (config.xbc == 1) { // periodic
                            if (i == 0 && (j > 0 && j < config.ny)) {
                                unsigned int ll = config.nx - 1 + config.offset * j;
                                C[l] = C0[ll];
                                continue;
                            }
                            if (i == config.nx && (j > 0 && j < config.ny)) {
                                unsigned int ll = 1 + config.offset * j;
                                C[l] = C0[ll];
                                continue;
                            }
                        }

                        C[l] = 0;
                    }
                }
            }
        }
    }


    void swap_one(double* f_old, double* f_new)
    {
        for (unsigned int l = 0; l < config.N; l++)
        {
            f_old[l] = f_new[l];
        }
    }

    void swap_three(double* f_old, double* f_new, double* f2_old, double* f2_new, double* f3_old, double* f3_new)
    {
        for (unsigned int l = 0; l < config.N; l++)
        {
            f_old[l] = f_new[l];
            f2_old[l] = f2_new[l];
            f3_old[l] = f3_new[l];
        }
    }


    struct CuPoisson
    {
        unsigned int k = 0;
        double eps = 0, res = 0, res0 = 0;
        double eps_iter = 1e-5;

        void solve(double *ksi, double *ksi0, double *omega)
        {
            k = 0;
            eps = 1.0;
            res = 0.0;
            res0 = 0.0;

            auto reduce = [this](double* f)
            {
                double s = 0;
                for (unsigned int l = 0; l < config.N; l++)
                    s += abs(f[l]);
                return s;
            };


            for (k = 1; k < 100000; k++)
            {
                poisson_stream(ksi, ksi0, omega);
                res = reduce(ksi);
                eps = abs(res - res0) / (res0 + 1e-5);
                res0 = res;

                std::swap(ksi, ksi0);

                if (eps < eps_iter)	break;
                if (k % 1000 == 0) std::cout << "k = " << k << ", eps = " << eps << std::endl;
            }
            //if (k > 100) std::cout << "device k = " << k << ", eps = " << eps << std::endl;
        }
    };
}