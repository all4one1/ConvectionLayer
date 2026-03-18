#pragma once

#include "FromOuterSparse/SparseMatrix.h"
#include "types_project.h"

namespace stream_cpu
{
	double vF2(unsigned int l, double* v, double* f, unsigned int shift, double h)
	{
		double FR, FL;

		double VR = (v[l + shift] + v[l]) * 0.5;
		double VL = (v[l] + v[l - shift]) * 0.5;

		if (VR > 0) FR = f[l];
		else if (VR < 0) FR = f[l + shift];

		if (VL > 0) FL = f[l - shift];
		else if (VL < 0) FL = f[l];
		return (VR * FR - VL * FL) / h;
	};

	void solveJacobi(double* f, double* f0, double* bb, int NN, SparseMatrix& M)
	{
		unsigned int k = 0;
		double eps_iter = 1e-5;

		for (k = 1; k < 100000; k++)
		{
			double s = 0;
			for (int j = 0; j < NN; j++)
			{
				s = M.line(j, f0);
				f[j] = f0[j] + (bb[j] - s) / M[j][j];
			}


			double max = 0;
			double dif;
			for (int i = 0; i < NN; i++)
			{
				dif = abs(f0[i] - f[i]);
				if (dif > max)
					max = dif;
			}
			for (int j = 0; j < NN; j++)
				f0[j] = f[j];

			if (max < eps_iter)	break;

			if (k % 1000 == 0) cout << "config k = " << k << ", eps = " << max << endl;
		}

	}
	void transform_to_velocity(Configuration& c, double* ksi, double* vx, double* vy)
	{
		unsigned int l, ll;
		for (unsigned int j = 0; j <= c.ny; j++) {
			for (unsigned int i = 0; i <= c.nx; i++) {
				l = i + c.offset * j;

				auto inner_x = [&](unsigned int q) {return (ksi[q + c.offset] - ksi[q - c.offset]) / (2.0 * c.hy);	};
				auto inner_y = [&](unsigned int q) {return -(ksi[q + 1] - ksi[q - 1]) / (2.0 * c.hx);	};

				if (i > 0 && i < c.nx && j > 0 && j < c.ny)
				{
					vx[l] = inner_x(l);
					vy[l] = inner_y(l);
					continue;
				}
				else
				{

					if (j == 0 || j == c.ny)
					{
						vx[l] = vy[l] = 0;
						continue;
					}

					if (c.xbc == closed)
					{
						if (i == 0 || i == c.nx)
						{
							vx[l] = vy[l] = 0;
							continue;
						}
					}
					else if (c.xbc == periodic)
					{
						if (i == 0)
						{
							ll = (c.nx - 1) + c.offset * j;
							vx[l] = inner_x(ll);
							vy[l] = inner_y(ll);
							continue;
						}
						if (i == c.nx)
						{
							ll = 1 + c.offset * j;
							vx[l] = inner_x(ll);
							vy[l] = inner_y(ll);
							continue;
						}
					}
				}
			}
		}
	}

	void form_matrix_with_2d_laplace_FD(SparseMatrix& M, Configuration& c, bc_type bcx, bc_type bcy, double coef = 1.0, double tau = 1e+10)
	{
		M.resize(c.N);

		double ax = -1. / (c.hx * c.hx) * coef;
		double ay = -1. / (c.hy * c.hy) * coef;

		double at = 1.0 / tau;

		for (unsigned int j = 0; j <= c.ny; j++) {
			for (unsigned int i = 0; i <= c.nx; i++) {
				unsigned int l = i + c.offset * j;
				std::map<int, double> m;
				/* INNER */
				if (i > 0 && i < c.nx && j > 0 && j < c.ny)
				{
					m[l - c.offset] = ay;
					m[l - 1] = ax;
					m[l] = (-2.0 * ax - 2.0 * ay) + at;
					m[l + 1] = ax;
					m[l + c.offset] = ay;

					M.add_line_with_map(m, l);
					continue;
				}
				else {
					if (bcy == Dirichlet)
					{
						if (j == 0) {
							m[l] = 1.0;
							M.add_line_with_map(m, l);
							continue;
						}
						else if (j == c.ny) {
							m[l] = 1.0;
							M.add_line_with_map(m, l);
							continue;
						}
					}
					else if (bcy == Neumann)
					{
						if (j == 0) {
							m[l] = 1.0;
							m[l + c.offset] = -1.0;
							M.add_line_with_map(m, l);
							continue;
						}
						else if (j == c.ny) {
							m[l] = 1.0;
							m[l - c.offset] = -1.0;
							M.add_line_with_map(m, l);
							continue;
						}
					}

					if (bcx == Dirichlet) { 
						if (i == 0 && (j > 0 && j < c.ny))
						{
							m[l] = 1.0;
							M.add_line_with_map(m, l);
						}
						if (i == c.nx && (j > 0 && j < c.ny))
						{
							m[l] = 1.0;
							M.add_line_with_map(m, l);
						}
						continue;
					}
					else if (bcx == periodic) { 
						if (i == 0 && (j > 0 && j < c.ny)) {
							unsigned int ll = c.nx - 1 + c.offset * j;
							m[l] = 1;
							m[ll] = -1;
							M.add_line_with_map(m, l);
							continue;
						}
						if (i == c.nx && (j > 0 && j < c.ny)) {
							unsigned int ll = 1 + c.offset * j;
							m[l] = 1;
							m[ll] = -1;
							M.add_line_with_map(m, l);
							continue;
						}
					}
					else if (bcx == Neumann) {
						if (i == 0 && (j > 0 && j < c.ny)) {
							m[l] = 1;
							m[l + 1] = -1;
							//m[l] = 3;	m[l + 1] = -4; m[l + 2] = 1;
							M.add_line_with_map(m, l);
							continue;
						}
						if (i == c.nx && (j > 0 && j < c.ny)) {
							m[l] = 1;
							m[l - 1] = -1;
							//m[l] = 3;	m[l - 1] = -4;	m[l - 2] = 1;
							M.add_line_with_map(m, l);
							continue;
						}
					}

				}
			}
		}
	}
	void form_matrix_with_2d_laplace_FV(SparseMatrix& M, Configuration& c, bc_type bcx, bc_type bcy, double coef = 1.0, double tau = 1e+10)
	{
		M.resize(c.N);

		double ax = c.Sx / c.hx * coef;
		double ay = c.Sy / c.hy * coef;
		double at = c.dV / tau;

		auto index = [&](Side side, unsigned int l, bool periodic = false)
		{
			if (side == Side::center)
				return l;

			else if (!periodic)
			{
				switch (side)
				{
				//case Side::center:
				//	return l;
				//	break;
				case Side::west:
					return l - 1;
					break;
				case Side::east:
					return l + 1;
					break;
				case Side::south:
					return l - c.offset;
					break;
				case Side::north:
					return l + c.offset;
					break;
				case Side::front:
					return l - c.offset2;
					break;
				case Side::back:
					return l + c.offset2;
					break;
				}
			}
			else if (periodic)
			{
				switch (side)
				{
				//case Side::center:
				//	return l;
				//	break;
				case Side::west:
					return l + (c.nx - 1);
					break;
				case Side::east:
					return l - (c.nx - 1);
					break;
				case Side::south:
					return l + (c.ny - 1) * c.offset;
					break;
				case Side::north:
					return l - (c.ny - 1) * c.offset;
					break;
				case Side::front:
					return l + (c.nz - 1) * c.offset2;
					break;
				case Side::back:
					return l - (c.nz - 1) * c.offset2;
					break;
				}
			}
			return UINT_MAX;
		};

		for (unsigned int j = 0; j < c.ny; j++) {
			for (unsigned int i = 0; i < c.nx; i++) {
				unsigned int l = i + c.offset * j;
				std::map<int, double> m;

				auto flux = [&](Side side, double a, bool border, bc_type bc)
				{
					if (border)
					{
						if (bc == Neumann)
						{
							//
						}
						if (bc == Dirichlet)
						{
							m[index(Side::center, l)] += a / 0.5;
						}
						if (bc == periodic)
						{
							m[index(Side::center, l)] += a;
							m[index(side, l, true)] = -a;
						}
					}
					else
					{
						m[index(Side::center, l)] += a;
						m[index(side, l)] = -a;
					}
				};

				flux(Side::west, ax, i == 0, bcx);
				flux(Side::east, ax, i == c.nx - 1, bcx);

				if (c.dim > 1)
				{
					flux(Side::south, ay, j == 0, bcy);
					flux(Side::north, ay, j == c.ny - 1, bcy);
				}

				m[index(Side::center, l)] += at;

				M.add_line_with_map(m, l);
			}
		}
	}

	void form_rhs_for_heat_equation(double* b, bool reset, SparseMatrix& M, Configuration& c, bc_type bcx, bc_type bcy, double coef = 1.0, double tau = 1e+10)
	{
		//auto vF = [this](ScalarVariable& f, int i, int j, int k, double Sx, double Sy, double Sz)
		//{
		//	double vxF = Sx * (f.get_at_side(Side::east, i, j, k) * vx.get_for_centered_cell(Side::east, i, j, k)
		//		- f.get_at_side(Side::west, i, j, k) * vx.get_for_centered_cell(Side::west, i, j, k));
		//	double vyF = Sy * (f.get_at_side(Side::north, i, j, k) * vy.get_for_centered_cell(Side::north, i, j, k)
		//		- f.get_at_side(Side::south, i, j, k) * vy.get_for_centered_cell(Side::south, i, j, k));
		//	double vzF = Sz * (f.get_at_side(Side::back, i, j, k) * vz.get_for_centered_cell(Side::back, i, j, k)
		//		- f.get_at_side(Side::front, i, j, k) * vz.get_for_centered_cell(Side::front, i, j, k));

		//	return (vxF + vyF + vzF);
		//};


		//for (unsigned int j = 0; j < c.ny; j++) {
		//	for (unsigned int i = 0; i < c.nx; i++) {
		//		unsigned int l = i + c.offset * j;
		//		std::map<int, double> m;
		//		if (reset) b[l] = 0.0;

		//		auto flux = [&](Side side, double a, bool border, bc_type bc)
		//		{
		//			if (border)
		//			{
		//				if (bc == Neumann)
		//				{
		//					//
		//				}
		//				if (bc == Dirichlet)
		//				{
		//					m[l] += F.boundary(side) * a / 0.5;
		//				}
		//				if (bc == periodic)
		//				{

		//				}
		//			}
		//			else
		//			{
		//				m[index(Side::center, l)] += a;
		//				m[index(side, l)] = -a;
		//			}
		//		};



		//		b[l] = F0(i, j, k) * DV / tau;
		//		b[l] += -vF(F0, i, j, k, SX * rx, SY * ry, SZ * rz);



		//		auto bc = [&F, b, i, j, k, l](Side side, double S, double h, double coef, bool border)
		//		{
		//			if (border)
		//			{
		//				if (F.boundary.type(side) == MathBoundary::Neumann)
		//				{
		//					b[l] += F.boundary.normal_deriv_oriented(side) * S * coef;
		//				}
		//				else if (F.boundary.type(side) == MathBoundary::Dirichlet)
		//				{
		//					b[l] += F.boundary(side) * coef * S / (0.5 * h);
		//				}
		//				else if (F.boundary.type(side) == MathBoundary::Periodic)
		//				{
		//					//nothing to do
		//				}
		//			}
		//		};


		//		bc(Side::west, SX, hx, PhysCoef, i == 0);
		//		bc(Side::east, SX, hx, PhysCoef, i == nx - 1);


		//		if (dim > 1)
		//		{
		//			bc(Side::south, SY, hy, PhysCoef, j == 0);
		//			bc(Side::north, SY, hy, PhysCoef, j == ny - 1);
		//		}

		//	}
		//}

	}

	void form_rhs_temperature(Configuration& c, double* rhs, double* f, double *vx, double *vy, bool reset = true, double tau = 1e+10)
	{
		auto vF = [&](double* f, double* vx, double* vy, unsigned int l)
		{
			double dvx = 0, dvy = 0;
			dvx = 0.5 * (f[l + 1] - f[l - 1]) / c.hx * vx[l];
			dvy = 0.5 * (f[l + c.offset] - f[l - c.offset]) / c.hy * vy[l];

			//dvx = vF2(l, vx, f, 1, c.hx);
			//dvy = vF2(l, vy, f, c.offset, c.hy);
			return dvx + dvy - (vx[l] * config.density_x + vy[l] * config.density_y);
		};

		unsigned int l; // ll;
		for (unsigned int j = 0; j <= c.ny; j++) {
			for (unsigned int i = 0; i <= c.nx; i++) {
				l = i + c.offset * j;
				if (reset) rhs[l] = 0.0;

				auto inner = [&](unsigned int q)
				{
					return f[q] / tau -vF(f, vx, vy, q);
				};


				if (i > 0 && i < c.nx && j > 0 && j < c.ny)
				{
					rhs[l] = inner(l);
				}
				else 
				{

					if (j == 0)
					{
						rhs[l] = 0.0;
						continue;
					}
					else if (j == c.ny)
					{
						rhs[l] = 0.0;
						continue;
					}
					
					if (c.xbc == closed)
					{
						if (i == 0)
							rhs[l] = 0;
						if (i == c.nx)
							rhs[l] = 0;
						continue;
					}
					else if (c.xbc == periodic)
					{
						//if (i == 0) ll = (c.nx - 1) + c.offset * j;
						//if (i == c.nx) ll = (1) + c.offset * j;

						rhs[l] = 0; // inner(ll);
					}
				}
			}
		}
	}
	void form_rhs_concentration(Configuration& c, double* rhs, double* f, double* vx, double* vy, bool reset = true, double tau = 1e+10)
	{
		auto vF = [&](double* f, double* vx, double* vy, unsigned int l)
		{
			double dvx = 0, dvy = 0;;
			dvx = 0.5 * (f[l + 1] - f[l - 1]) / c.hx * vx[l];
			dvy = 0.5 * (f[l + c.offset] - f[l - c.offset]) / c.hy * vy[l];

			//dvx = vF2(l, vx, f, 1, c.hx);
			//dvy = vF2(l, vy, f, c.offset, c.hy);
			return dvx + dvy - (vx[l] * config.density_x + vy[l] * config.density_y);
		};

		unsigned int l;  //, ll;
		for (unsigned int j = 0; j <= c.ny; j++) {
			for (unsigned int i = 0; i <= c.nx; i++) {
				l = i + c.offset * j;
				if (reset) rhs[l] = 0.0;

				auto inner = [&](unsigned int q)
				{
					return f[q] / tau - vF(f, vx, vy, q);
				};


				if (i > 0 && i < c.nx && j > 0 && j < c.ny)
				{
					rhs[l] = inner(l);
				}
				else
				{

					if (j == 0)
					{
						rhs[l] = 0.0;
						continue;
					}
					else if (j == c.ny)
					{
						rhs[l] = 0.0;
						continue;
					}

					if (c.xbc == closed)
					{
						if (i == 0)
							rhs[l] = 0;
						if (i == c.nx)
							rhs[l] = 0;
						continue;
					}
					else if (c.xbc == periodic)
					{
						//if (i == 0) ll = (c.nx - 1) + c.offset * j;
						//if (i == c.nx) ll = (1) + c.offset * j;

						rhs[l] = 0; // inner(ll);
					}
				}
			}
		}
	}
	void form_rhs_vorticity(Configuration& c, double* rhs, double *T, double *C, double *ksi, double* omega, double* vx, double* vy, bool reset = true, double tau = 1e+10)
	{
		auto vF = [&](double* f, double* vx, double* vy, unsigned int l)
		{
			double dvx = 0, dvy = 0;
			dvx = 0.5 * (f[l + 1] - f[l - 1]) / c.hx * vx[l];
			dvy = 0.5 * (f[l + c.offset] - f[l - c.offset]) / c.hy * vy[l];
			return dvx + dvy;
		};
		unsigned int l; // ll;
		for (unsigned int j = 0; j <= c.ny; j++) {
			for (unsigned int i = 0; i <= c.nx; i++) {
				l = i + c.offset * j;
				if (reset) rhs[l] = 0.0;

				auto inner = [&](unsigned int q)
				{
					return omega[l] / tau
						//+ (-vy[l] * dy1(l, omega) - vx[l] * dx1(l, omega))
						- vF(omega, vx, vy, l)
						+ c.grav_y * c.Ra / c.Pr * (dx1(l, T) - c.density_x)
						- c.grav_x * c.Ra / c.Pr * (dy1(l, T) - c.density_y)

						+ c.grav_y * c.Ra / c.Pr * c.K * (dx1(l, C) - c.density_x)
						- c.grav_x * c.Ra / c.Pr * c.K * (dy1(l, C) - c.density_y);
				};


				if (i > 0 && i < c.nx && j > 0 && j < c.ny)
				{
					rhs[l] = inner(l);
				}
				else
				{
					if (j == 0)
					{
						rhs[l] = -0.5 / (c.hy * c.hy) * (8.0 * ksi[l + c.offset] - ksi[l + c.offset * 2]);
						continue;
					}
					else if (j == c.ny)
					{
						rhs[l] = -0.5 / (c.hy * c.hy) * (8.0 * ksi[l - c.offset] - ksi[l - c.offset * 2]);
						continue;
					}

					if (c.xbc == closed)
					{
						if (i == 0)
						{
							rhs[l] = -0.5 / (c.hx * c.hx) * (8.0 * ksi[l + 1] - ksi[l + 2]);
							continue;
						}
						if (i == c.nx) {
							rhs[l] = -0.5 / (c.hx * c.hx) * (8.0 * ksi[l - 1] - ksi[l - 2]);
							continue;
						}
					}
					else if (c.xbc == periodic)
					{
						//if (i == 0) ll = (c.nx - 1) + c.offset * j;
						//if (i == c.nx) ll = (1) + c.offset * j;

						rhs[l] = 0; // inner(ll);
					}
				}
			}
		}
	}
	void form_rhs_stream(Configuration& c, double* rhs, double* omega, bool reset = true)
	{
		unsigned int l; // ll;
		for (unsigned int j = 0; j <= c.ny; j++) {
			for (unsigned int i = 0; i <= c.nx; i++) {
				l = i + c.offset * j;
				if (reset) rhs[l] = 0.0;

				auto inner = [&](unsigned int q)
				{
					return omega[l];
				};

				if (i > 0 && i < c.nx && j > 0 && j < c.ny)
				{
					rhs[l] = inner(l);
				}
				else
				{
					if (j == 0)
					{
						rhs[l] = 0;
						continue;
					}
					else if (j == c.ny)
					{
						rhs[l] = 0;
						continue;
					}

					if (c.xbc == closed)
					{
						if (i == 0)
						{
							rhs[l] = 0;
							continue;
						}
						if (i == c.nx) {
							rhs[l] = 0;
							continue;
						}
					}
					else if (c.xbc == periodic)
					{
						//if (i == 0) ll = (c.nx - 1) + c.offset * j;
						//if (i == c.nx) ll = (1) + c.offset * j;

						rhs[l] = 0; // inner(ll);
					}
				}
			}
		}
	}

	

	struct ImplicitStream
	{
		unsigned int N;
		SparseMatrix SMT, SMC, SMO, SMP;

		ImplicitStream(unsigned int N_, Configuration& c) : N(N_)
		{
			form_matrix_with_2d_laplace_FD(SMT, c, c.xbc == periodic ? periodic : Neumann, Dirichlet, (1.0 / c.Pr), c.tau);
			form_matrix_with_2d_laplace_FD(SMC, c, c.xbc == periodic ? periodic : Neumann, Neumann, (1.0 / (c.Le * c.Pr)), c.tau);
			form_matrix_with_2d_laplace_FD(SMO, c, c.xbc == periodic ? periodic : Dirichlet, Dirichlet, 1.0, c.tau);
			form_matrix_with_2d_laplace_FD(SMP, c, c.xbc == periodic ? periodic : Dirichlet, Dirichlet, 1.0);
		}
	};
}