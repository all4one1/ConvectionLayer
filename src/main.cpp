#include <vector>
#include <iostream>
#include "Extras.h"
#include "types_project.h"
Configuration config;
#include "init.h"
#include "CPUSolvers/cpu_stream.h"
#include "CPUSolvers/cpu_stream_impl.h"
using std::cout;
using std::endl;




RUN_STATE run;

int main(int argc, char** argv)
{
	ReadingFile par("parameters.txt");
	par.reading<int>(run.read_recovery, "continue", 0); // 0 - full start, 1 - full continue, 2 - read fields, state is clear
	par.reading<int>(run.cpu_mode, "cpu", 1);
	par.reading<double>(run.timeq_limit, "time_limit", 20000);
	par.reading<double>(run.timeq_minimal, "time_minimal", 0);
	par.reading<int>(run.call_max, "call_max", 100);
	par.reading_string(run.note, "note", "");

	init_parameters(config);
	unsigned int &N = config.N;

	// init objects
	Arrays array;
	FuncTimer ftimer;
	StatValues stat;
	Backup backup("recovery", true);
	std::ofstream w_final, w_temporal;

	// for the stream function - vorticity method
	allocate_host_arrays({ &array.T, &array.T0, &array.C, &array.C0, 
		&array.omega, &array.omega0, &array.ksi, &array.ksi0, 
		&array.buffer, &array.buffer2, &array.vx, &array.vy, &array.rhs }, config.N);

	stream_cpu::CuPoisson cpuPoisson;
	stream_cpu::ImplicitStream IS(N, config);

	//reading state from file
	if (run.read_recovery)
	{
		int error = backup.read(run.iter, run.timeq, run.call_i, config, { array.ksi, array.omega, array.T, array.C }, run.read_recovery == 1);

		if (error == 0)
		{
			if (run.cpu_mode)
			{
				for (unsigned int l = 0; l < config.N; l++)
				{
					array.ksi0[l] = array.ksi[l];
					array.omega0[l] = array.omega[l];
					array.T0[l] = array.T[l];
					array.C0[l] = array.C[l];
				}
			}
			auto open_type = run.read_recovery == 1 ? std::ofstream::app : std::ofstream::out;
			w_final.open("w_final.dat", open_type);
			w_temporal.open("w_temporal.dat", open_type);
			if (run.read_recovery == 2) deleteFilesInDirectory(L"fields");
		}
		if (error == 1)
		{
			run.read_recovery = 0;
		}
	}
	if (!run.read_recovery)
	{
		w_final.open("w_final.dat");
		w_temporal.open("w_temporal.dat");
		deleteFilesInDirectory(L"fields");
		init_fields(config, array);
	}


reset: //we go here when a steady state has been reached and we start with a new value of Ra
	if (run.call_i >= run.call_max) return 0;
	if (run.call_i > 0 && run.iter == 0)	config.Ra += config.incr_parameter;
	
	//objects to check if steady state has been established
	Checker check_ksi(&stat.ksi_sum, &run.timeq, Checker::ExitType::Relative, "ksi", 1e-6);
	Checker check_omega(&stat.omega_sum, &run.timeq, Checker::ExitType::Relative, "omega");
	Checker check_C(&stat.C_sum, &run.timeq, Checker::ExitType::Relative, "C");

	ftimer.start("main");

	while (run.stop_signal == 0)
	{
		run.iter++;
		run.timeq += config.tau;
		if (run.timeq > run.timeq_limit) run.stop_signal = 2;


		ftimer.start("calc");

		if (run.cpu_mode == 1) // explicit
		{
			stream_cpu::vorticity_quadratic_temperature(array.omega, array.omega0, array.ksi, array.T);
			stream_cpu::temperature_2d_flux_full(array.T, array.T0, array.ksi);

			stream_cpu::swap_three(array.omega0, array.omega, array.T0, array.T, array.C0, array.C);
			cpuPoisson.solve(array.ksi, array.ksi0, array.omega);
		}
		if (run.cpu_mode == 2) //implicit
		{
			//implemented for old tasks
		}


		ftimer.end("calc");
		
		// OUTPUT
		if (run.iter == 1 
			|| (run.every_time(config.tau, 1.0)  && run.timeq < 1000) 
			|| (run.every_time(config.tau, 10.0) && run.timeq >= 1000))
		{

			//compute statistics
			Nu_y_for_fixed_flux(config, array.T, stat.Nu);
			stat.ksi_max = absmax(array.ksi, N);
			stat.ksi_sum = sum_abs(array.ksi, N);
			stat.omega_sum = sum_abs(array.omega, N);
			stat.C_sum_signed = sum_signed(array.C, N);
			check_ksi.update();
			check_omega.update();
			stream_cpu::transform_to_velocity(config, array.ksi, array.vx, array.vy);
			velocity_stats(config, array.vx, array.vy, stat);

			//print on screen
			cout << endl << "Ra = " << config.Ra << ", t= " << run.timeq << ", " << run.iter << endl;
			cout << "ksi= " << stat.ksi_max << ", Vmax= " << stat.Vmax << ", Nu=" << stat.Nu << endl;
			if (!run.note.empty()) cout << "note: " << run.note << endl;

			// write to time-dependent values in file
			if (run.iter == 1) w_temporal << "t, time(sec), time(sec)v2, max_ksi, omega_sum, ksi_point, T_point, Vmax, Ek, Nu" << " Ra=" << config.Ra << endl;
			w_temporal << run.timeq << " " << ftimer.update_and_get("main") << " " << ftimer.get("calc")
				<< " " << stat.ksi_max << " " << stat.omega_sum << " " << array.ksi[INDEX(10, 10, 0)] << " " << array.T[INDEX(10, 10, 0)] 
				<< " " << stat.Vmax << " " << stat.Ek << " " << stat.Nu
				<< endl;

			// check if steady state
			if (run.timeq >= run.timeq_minimal)
			{
				cout << "ksi check   = " << check_ksi.dif << " " << check_ksi.dif_rel << endl;
				cout << "omega check = " << check_omega.dif << " " << check_omega.dif_rel << endl;
				if (check_ksi.ready_to_exit) run.stop_signal = 1;
			}

			// write fields 
			if (run.every_time(config.tau, 20))
			{
				write_fields2d(_path("fields"), _str(config.Ra) + " " + _str(run.timeq), config,
					{ array.ksi, array.omega, array.T, array.C, array.vx, array.vy },
					"ksi, omega, T, C, vx, vy"
				);
			}

			//write recovery file to continue 
			if (run.every_time(config.tau, 20))
			{
				backup.save(run.iter, run.timeq, run.call_i, config, { array.ksi, array.omega, array.T, array.C }, "ksi, omega, T, C");
			}
		}

		//write final state when a steady state has been reached
		if (run.stop_signal > 0)
		{
			if (run.call_i == 0)	w_final << "Ra, ksi_max, ksi_max2, omega_sum, Vmax, Ek, time(sec), t, Nu" << endl;

			w_final << config.Ra << " " << stat.ksi_max << " " << pow(stat.ksi_max, 2) << " " << stat.omega_sum << " " << stat.Vmax << " " << stat.Ek
				<< " " << ftimer.update_and_get("main") << " " << run.timeq
				<< " " << stat.Nu
				<< endl;
			backup.save(run.iter, run.timeq, run.call_i, config,	{ array.ksi, array.omega, array.T, array.C },	"ksi, omega, T, C");
			
			write_fields2d(_path("final"), _str(config.Ra) + " " + _str(run.timeq), config,
				{ array.ksi, array.omega, array.T, array.C, array.vx, array.vy },
				"ksi, omega, T, C, vx, vy");

			if (run.stop_signal == 1)
			{
				run.iter = 0; run.timeq = 0; run.stop_signal = 0; 
				run.call_i++;
				goto reset;
			}
			if (run.stop_signal == 2) break;
		}

		if (run.stop_signal == -1)
		{
			break;
		}
	} cout << "run.stop: " << run.stop_signal << endl;


	return 0;
}
