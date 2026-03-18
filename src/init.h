#pragma once

//folders
std::string folder_fields;
std::string folder_final;
std::string folder_lines;

void init_parameters(Configuration& c)
{
	addFolder("fields", folder_fields);
	addFolder("final", folder_final);
	addFolder("lines", folder_lines);

	ReadingFile par("parameters.txt");

	int bc;
	par.reading<int>(bc, "xbc", 0); c.xbc = static_cast<bc_type>(bc);
	par.reading<int>(bc, "ybc", 0); c.ybc = static_cast<bc_type>(bc);

	par.reading<unsigned int>(c.heatflux, "heatflux", 1);

	par.reading<double>(c.Ly, "Ly", 1.0);
	par.reading<double>(c.Lx, "Lx", 1.0);
	par.reading<unsigned int>(c.nx, "nx", 20);
	par.reading<unsigned int>(c.ny, "ny", 20);
	c.nx = (unsigned int)c.Lx * c.nx;
	c.ny = (unsigned int)c.Ly * c.ny;
	c.nz = (unsigned int)c.Ly * c.nz;
	par.reading<double>(c.alpha, "alpha", 0);
	par.reading<double>(c.Ra, "Ra", 5000);
	par.reading<double>(c.Rad, "Rad", 5000);
	par.reading<double>(c.Pr, "Pr", 10);
	par.reading<double>(c.Le, "Le", 130);
	par.reading<double>(c.Sc, "Sc", 1000);
	par.reading<double>(c.psi, "psi", -0.2);
	par.reading<double>(c.K, "K", 0);
	par.reading<double>(c.tau, "tau", 0.0001);
	par.reading<double>(c.grav_y, "grav_y", 1.0);
	par.reading<double>(c.incr_parameter, "increment", -20);


	if (c.nx > 0) c.dim = 1;
	if (c.ny > 0) c.dim = 2;
	if (c.nz > 0) c.dim = 3;

	c.hx = c.Lx / c.nx;
	c.hy = c.Ly / c.ny;
	c.hz = c.Lz / c.nz;

	unsigned int q = 1; //1 = FD, 0 = FV
	if (c.dim == 1)
	{
		c.ny = c.nz = 1 - q;
		c.hy = c.hz = 0;
		c.Sx = 1;
		c.Sy = c.Sz = 0;
		c.dV = c.hx;
		c.N = (c.nx + q);
	}

	if (c.dim == 2)
	{
		c.nz = 1 - q;
		c.hz = 0;
		c.Sx = c.hy;
		c.Sy = c.hx;
		c.Sz = 0;
		c.dV = c.hx * c.hy;
		c.N = (c.nx + q) * (c.ny + q);
		c.offset = c.nx + q;
		c.offset2 = 0;
	}

	if (c.dim == 3)
	{
		c.Sx = c.hy * c.hz;
		c.Sy = c.hx * c.hz;
		c.Sz = c.hx * c.hy;
		c.dV = c.hx * c.hy * c.hz;
		c.N = (c.nx + q) * (c.ny + q) * (c.nz + q);
		c.offset = c.nx + q;
		c.offset2 = (c.nx + q) * (c.ny + q);
	}

	c.Nbytes = c.N * sizeof(double);

	double pi = 3.1415926535897932384626433832795;
	auto make_vector = [&pi](double angle, double (*func)(double))
	{	return std::floor(func(angle * pi / 180.0) * 1e+10) / 1e+10;	};
	auto set_angles = [&make_vector, &c](double a, double b)
	{
		c.grav_y = make_vector(a, cos);
		c.grav_x = make_vector(a, sin);

		c.vibr_y = make_vector(b, cos);
		c.vibr_x = make_vector(b, sin);

		c.density_y = make_vector(a, cos);
		c.density_x = make_vector(a, sin);
	};

	set_angles(c.alpha, c.beta);

}


void init_fields(Configuration& c, Arrays& h)
{
	double PI = acos(-1);
	//for (unsigned int k = 0; k <= c.nz; k++) {
		for (unsigned int j = 0; j <= c.ny; j++) {
			for (unsigned int i = 0; i <= c.nx; i++)
			{
				unsigned int q = i + config.offset * j;
				double y = c.hy * j;
				double x = c.hx * i;
				h.C[q] = 0; // 1 - y;
				h.T[q] = 0; // 1 - y;

				//if (y > 0.5) h.C[q] += 0.5;
				//else h.C[q] += -0.5;
				h.ksi[q] = (1e-6) * sin(2 * PI * x / c.Lx) * sin(PI * y);
				//h.C0[q] = h.C[q] = (1e-5) * sin(2 * PI * x / c.Lx) * (1 - y) * c.psi;
			}
		}
	//}

}