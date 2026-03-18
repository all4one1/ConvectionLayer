#pragma once

#include <string>
#include <windows.h>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <fstream>
#include <algorithm>
#include <vector>
#include <map>
#include "types_project.h"


using std::string;
using std::ifstream;
using std::stringstream;
using std::istringstream;
using std::ostringstream;
using std::cout;
using std::endl;


struct RUN_STATE
{
	int stop_signal = 0, read_recovery = 0, call_i = 0, call_max = 1000, cpu_mode = 0;
	double timeq = 0, time_sec = 0, timeq_limit = 10000, timeq_minimal = 0;
	size_t iter = 0;
	std::string note = "";
	int test = 0;
	bool every(int n) { return iter % n == 0; }
	bool every_time(double tau, double t) { int tt = (int)(ceil(1.0 / tau) * t); return iter % (tt) == 0; }
	void reset()
	{
		iter = 0;
		timeq = 0;
		stop_signal = 0;
	}
};

struct ReadingFile
{
private:
	ifstream read;
	string str, substr, buffer;
	string file_name;
	stringstream ss;
	istringstream iss;
	ostringstream oss;
	int stat = 0, pos = 0;

public:
	ReadingFile(string name)
	{
		file_name = name;
		open_file(file_name);
		stat = 0;
	}
	ReadingFile()
	{
		stat = 0;
	}

	void open_file(string file_name) {
		read.open(file_name.c_str());
		if (read.good()) {
			cout << endl << "the parameter file \"" << file_name << "\" has been read " << endl << endl;
			#ifdef LoggingON
			Log << endl << "the parameter file \"" << file_name << "\" has been read " << endl << endl;
			#endif // LoggingON

			oss << read.rdbuf();
			buffer = oss.str();
			iss.str(buffer);
		}
		else {
			cout << "the parameter file has been not found, default parameters will be initialized " << endl;
			#ifdef LoggingON
			Log << "the parameter file has been not found, default parameters will be initialized " << endl;
			#endif // LoggingON
			buffer = "";
			iss.str(buffer);
		}
	}

	//template <typename T>
	//int reading(std::map<T, string> &m, T def_var = 0, T min = 0, T max = 0)
	//{
	//	return reading<T>(m.first, m.second, dev_var, min, max);
	//}

	template <typename T>
	int reading(T& var, string parameter_name, T def_var = 0, T min = 0, T max = 0) {
		int ret = 0;
		stat = 0;
		transform(parameter_name.begin(), parameter_name.end(), parameter_name.begin(), ::tolower);
		iss.clear();
		iss.seekg(0);

		while (getline(iss, str))
		{
			//substr.clear();
			ss.str("");	ss.clear();	ss << str;	ss >> substr;
			transform(substr.begin(), substr.end(), substr.begin(), ::tolower);
			if (substr == parameter_name) {
				ret = 1;
				pos = (int)ss.tellg();
				while (ss >> substr) {
					if (substr == "=")
					{
						ss >> var;
						stat = 1;
						break;
					}
				}

				if (stat == 0) {
					ss.clear();
					ss.seekg(pos);
					ss >> var;
				}
				break;
			}
		}
		if (iss.fail())
		{
			var = def_var;
		}
		if (min != max && (min + max) != 0) {
			if (var > max || var < min)
			{
				cout << "Warning: \"" + parameter_name + "\" should not be within this range" << endl;
				var = def_var;
			}
		}


		#ifdef LoggingON
		Log << parameter_name << "= " << var << endl;
		#endif // LoggingON

		return ret; //return 1 if read
	}

	void reading_string(string& var, string parameter_name, string def_var) {
		stat = 0;
		transform(parameter_name.begin(), parameter_name.end(), parameter_name.begin(), ::tolower);
		iss.clear();
		iss.seekg(0);

		while (getline(iss, str))
		{
			//substr.clear();
			ss.str("");	ss.clear();	ss << str;	ss >> substr;
			transform(substr.begin(), substr.end(), substr.begin(), ::tolower);
			if (substr == parameter_name) {
				pos = (int)ss.tellg();
				while (ss >> substr) {
					if (substr == "=")
					{
						ss >> var;
						stat = 1;
						break;
					}
				}

				if (stat == 0) {
					ss.clear();
					ss.seekg(pos);
					ss >> var;
				}
				break;
			}
		}
		if (iss.fail())
		{
			var = def_var;
		}
		#ifdef LoggingON
		Log << parameter_name << "= " << var << endl;
		#endif // LoggingON
	}

	bool reading_string(string parameter_name) {
		stat = 0;
		transform(parameter_name.begin(), parameter_name.end(), parameter_name.begin(), ::tolower);
		iss.clear();
		iss.seekg(0);

		while (getline(iss, str))
		{
			//substr.clear();
			ss.str("");	ss.clear();	ss << str;	ss >> substr;
			transform(substr.begin(), substr.end(), substr.begin(), ::tolower);
			if (substr == parameter_name) 
			{
				return true;
			}
		}
		#ifdef LoggingON
		Log << parameter_name << "= " << var << endl;
		#endif // LoggingON
		return false;
	}
};

struct FuncTimer
{
private:
	std::map <std::string, double> timer, t1, t2;
	std::string active;
	double dt, total;
public:
	FuncTimer()
	{
		dt = total = 0;
	}
	void start(std::string s)
	{
		active = s;
		t1[s] = clock();
	}
	void end(std::string s)
	{
		if (t1.find(s) == t1.end())
		{
			std::cout << s + " trigger not started" << std::endl;
			return;
		}
		else
		{
			t2[s] = clock();
			dt = (t2[s] - t1[s]) / CLOCKS_PER_SEC;
			//if (s == active)
			{
				timer[s] += dt;
			}
		}
	}
	double get(std::string s)
	{
		if (t1.find(s) == t1.end())
		{
			std::cout << s + " trigger not started" << std::endl;
			return 0.0;
		}
		else
		{
			return timer[s];
		}
	}

	double update_and_get(std::string s)
	{
		if (t1.find(s) == t1.end())
		{
			std::cout << s + " trigger not started" << std::endl;
			return 0.0;
		}
		else
		{
			t2[s] = clock();
			dt = (t2[s] - t1[s]) / CLOCKS_PER_SEC;
			return timer[s] + dt;
		}
	}
	double get_last_diff(std::string s)
	{
		if (t1.find(s) == t1.end())
		{
			std::cout << s + " trigger not started" << std::endl;
			return 0.0;
		}
		else
		{
			t2[s] = clock();
			dt = (t2[s] - t1[s]) / CLOCKS_PER_SEC;
			return dt = (t2[s] - t1[s]) / CLOCKS_PER_SEC;
		}
	}

	std::string get_info()
	{
		int n = int(timer.size());
		std::ostringstream oss;
		oss << "Calculation time in seconds. Number of cases: " << n << ".\n";

		for (auto& it : timer)
		{
			oss << it.first << ": " << it.second << std::endl;
		}
		return oss.str();
	}
	void show_info()
	{
		std::cout << get_info() << std::endl;
	}
	void write_info(std::string path = "timer_report.dat")
	{
		std::ofstream file(path);
		file << get_info() << std::endl;
	}

};


size_t allocate_host_arrays(std::vector<double**> v, size_t N)
{
	size_t Nbytes = sizeof(double) * N;
	size_t total_bytes = 0;
	for (auto it : v)
	{
		*it = new double[N];
		memset(*it, 0, Nbytes);
		total_bytes += Nbytes;
	}
	return total_bytes;
}

struct Checker
{
	enum class ExitType { Relative, Absolute } type;
	std::string name;
	std::vector<double> v, t;
	double* binded_value = nullptr, * binded_time = nullptr;
	double t_estimated = 0.0;
	double last_time = 0;
	size_t N = 0, i2 = 0, i1 = 0, i0 = 0, iter = 0;
	double eps_default = 1e-6;
	double relative = 1, absolute = 1;
	double dif = 1, dif_rel = 1, dif2 = 0, dif_rel2 = 0;
	double deriv = 0, deriv2 = 0;
	size_t i_check = 0, check_limit = 20;
	bool ready_to_exit = false;

	Checker(double* ptr, double* t_ptr, ExitType type_ = ExitType::Relative, std::string name_ = "Checker", double eps = -1)
	{
		iter = 0;
		N = 3;
		binded_value = ptr;
		binded_time = t_ptr;
		v.resize(N);
		std::fill(v.begin(), v.end(), 0);

		t.resize(N);
		std::fill(t.begin(), t.end(), 0);

		type = type_;
		name = name_;

		if (eps > 0) eps_default = eps;
	}

	void update()
	{
		if (last_time == binded_time[0]) return;
		i2 = iter % N;
		i1 = (iter - 1) % N;
		i0 = (iter - 2) % N;

		v[i2] = binded_value[0];
		t[i2] = binded_time[0];
		iter++;
		last_time = binded_time[0];
		if (iter < N) return;

		dif = v[i2] - v[i1];
		dif_rel = (v[i2] - v[i1]) / v[i1];
		//dif2 = (v[i2] - v[i1]) - (v[i1] - v[i0]);
		//dif_rel2 = ((v[i2] - v[i1]) - (v[i1] - v[i0])) / (v[i2] - v[i1]);
		double dt = t[i2] - t[i1];
		deriv = (v[i2] - v[i1]) / dt;
		deriv2 = (v[i2] - 2.0 * v[i1] + v[i0]) / (dt * dt);


		absolute = abs(dif);
		relative = abs(dif / v[i1]);

		check();
		//cout << iter << ": " << i2 << " " << i1 << " " << i0 << endl;
		//cout << iter << ": " << v[i2] << " " << v[i1] << " " << v[i0] << endl;
	}

	void update(double& val, double& t)
	{
		binded_value = &val;
		binded_time = &t;
		this->update();
	}

	bool check(double eps = -1)
	{
		update();
		if (eps < 0) eps = eps_default;
		if (relative < eps && type == ExitType::Relative
			|| absolute < eps && type == ExitType::Absolute)
		{
			i_check++;
			cout << name + ",i= " << i_check << endl;
			if (i_check >= check_limit)
			{
				ready_to_exit = true;
				return true;
			}
		}
		else
		{
			i_check = 0;
		}
		ready_to_exit = false;
		return false;
	}

	double estimate_finish(double eps = -1)
	{
		update();
		if (eps < 0) eps = eps_default;

		double e2 = (v[i2] - v[i1]) / v[i1];
		double e1 = (v[i1] - v[i0]) / v[i0];
		double dt = t[i2] - t[i1];

		return t[i2] + dt * (eps - e2) / (e2 - e1);
	}

	bool operator()()
	{
		return ready_to_exit;
	}

	//bool check_long_term();
};


#define INDEX(i, j, k) i + config.offset * j + config.offset2 * k

void pauseWithLineInfo(const std::string& file, int line)
{
	std::cout << "Paused at file: " << file << ", line: " << line << std::endl;
	std::cout << "Press Enter to continue...";
	std::cin.ignore();
}
#define pause pauseWithLineInfo(__FILE__, __LINE__);


void write_fields2d(std::string folder, std::string name, Configuration& c, std::vector<double*> v = {}, 
	std::string head = "T, vx, vy, vz, C, Psi")
{
	string str = folder + name + ".txt";

	std::ofstream all(str.c_str());

	all << "x, y, " + head << endl;
	
	//all << setprecision(16) << fixed;
	for (unsigned int j = 0; j <= c.ny; j++) {
		for (unsigned int i = 0; i <= c.nx; i++) {
			unsigned int l = i + c.offset * j;
			all << i * c.hx << " " << j * c.hy << " ";
			for (auto& it : v)
			{
				all << it[l] << " ";
			}
			all << endl;
		}
	}
};

void write_x_line(unsigned int j, std::string folder, std::string name, Configuration& c, std::vector<double*> v = {}, std::string head = "T, vx, vy, vz, C, Psi")
{
	string str = folder + name + ".txt";

	std::ofstream w(str.c_str());

	w << "x, " + head << endl;
	for (unsigned int i = 0; i <= c.nx; i++) 
	{
		unsigned int l = i + c.offset * j;
		w << i * c.hx << " " << j * c.hy << " ";
		for (auto& it : v)
		{
			w << it[l] << " ";
		}
		w << endl;
	}
}

//std::string section_x_max(unsigned int j, Configuration &c, std::vector<double *> v)
//{
//	size_t n = v.size();
//	std::vector<double> maxv(n);
//
//	for (unsigned int i = 0; i <= c.nx; i++)
//	{
//		for (size_t i = 0; i < n; i++)
//		{
//			maxv[i] = v[i][0];
//		}
//
//		unsigned int l = i + c.offset * j;
//		for (size_t i = 0; i < n; i++)
//		{
//			if (v[i][l] > maxv[i])
//				maxv[i] = v[i][l];
//		}
//	}
//}

void deleteFilesInDirectory(const std::wstring& directoryPath) {
	WIN32_FIND_DATAW findFileData;
	HANDLE hFind = INVALID_HANDLE_VALUE;

	// Create a search path for all files in the directory
	std::wstring searchPath = directoryPath + L"\\*";

	// Find the first file in the directory
	hFind = FindFirstFileW(searchPath.c_str(), &findFileData);
	if (hFind == INVALID_HANDLE_VALUE) {
		std::wcerr << L"Error opening directory: " << directoryPath << std::endl;
		return;
	}

	do {
		// Skip the special directories "." and ".."
		if (wcscmp(findFileData.cFileName, L".") == 0 || wcscmp(findFileData.cFileName, L"..") == 0) {
			continue;
		}

		// Construct the full file path
		std::wstring filePath = directoryPath + L"\\" + findFileData.cFileName;

		// Check if the entry is a regular file
		if (!(findFileData.dwFileAttributes & FILE_ATTRIBUTE_DIRECTORY)) {
			// Delete the file
			if (DeleteFileW(filePath.c_str()))
			{
				//std::wcout << L"Deleted file: " << filePath << std::endl;
			}
			else {
				std::wcerr << L"Error deleting file: " << filePath << std::endl;
			}
		}
	} while (FindNextFileW(hFind, &findFileData) != 0);

	// Close the search handle
	FindClose(hFind);
}
std::streamoff fileSize(const std::string& filePath) { 
	std::ifstream file(filePath, std::ios::binary | std::ios::ate); 
	if (!file.is_open()) { std::cerr << "Error!" << std::endl;        return false; }
	std::streampos bytes = file.tellg();
	return bytes.operator std::streamoff(); 	//file size in bytes
}

void addFolder(std::string fname, std::string &global)
{
	#ifdef __linux__ 
		string str = "mkdir -p " + fname +"/";
		global = fname + "/";
		system(str.c_str());
	#endif
	
	#ifdef _WIN32
		CreateDirectoryA(fname.c_str(), NULL);
		global = fname + "//";
	#endif
}

std::string _path(std::string fname)
{
	#ifdef __linux__ 
	return fname + "/";
	#endif

	#ifdef _WIN32
	return fname + "//";
	#endif
}



double sum_signed(double* f, unsigned int N)
{
	double s = 0;

	for (unsigned int i = 0; i < N; i++)
		s += f[i];
	return s;
}
double sum_abs(double* f, unsigned int N)
{
	double s = 0;

	for (unsigned int i = 0; i < N; i++)
		s += abs(f[i]);
	return s;
}
double signedmax(double* f, unsigned int N)
{
	double m = f[0];

	for (unsigned int i = 0; i < N; i++)
	{
		if (f[i] > m) m = f[i];
	}
	return m;
}
double signedmin(double* f, unsigned int N)
{
	double m = f[0];

	for (unsigned int i = 0; i < N; i++)
	{
		if (f[i] < m) m = f[i];
	}
	return m;
}
double absmax(double *f, unsigned int N)
{
	double m = f[0];

	for (unsigned int i = 0; i < N; i++)
	{
		if (abs(f[i]) > m) m = abs(f[i]);
	}
	return m;
}


std::string _str(double n, int pres = 5)
{
	std::stringstream ss;
	ss << std::setprecision(pres) << n;
	return ss.str();
}
std::string _str(int n, int pres = 5)
{
	std::stringstream ss;
	ss << std::setprecision(pres) << n;
	return ss.str();
}
std::string _str(size_t n, int pres = 5)
{
	std::stringstream ss;
	ss << std::setprecision(pres) << n;
	return ss.str();
}

void velocity_stats(Configuration &c, double *vx, double *vy, StatValues &stat) {
	double V = 0, VX = 0, VY = 0;
	stat.Ek = stat.Vmax = stat.Vx = stat.Vy = 0.0;

	unsigned int l = 0;
	for (unsigned int j = 0; j <= c.ny; j++) {
		for (unsigned int i = 0; i <= c.nx; i++) {
			l = i + c.offset * j;
			VX = abs(vx[l]);
			VY = abs(vy[l]);

			V = vx[l] * vx[l] + vy[l] * vy[l];
			stat.Ek += V;
			if (sqrt(V) > stat.Vmax) stat.Vmax = sqrt(V);
			if (VX > stat.Vx) stat.Vx = VX;
			if (VY > stat.Vy) stat.Vy = VY;
		}
	}
	stat.Ek = stat.Ek / 2.0 * c.hx * c.hy;

	stat.Cu = stat.Vx * c.tau / c.hx + stat.Vy * c.tau / c.hy;

	std::vector<double> Pe = { (stat.Vx * c.hx + stat.Vy * c.hy) * c.Pr, (stat.Vx * c.hx + stat.Vy * c.hy) * 1, (stat.Vx * c.hx + stat.Vy * c.hy) * (c.Le * c.Pr) };

	stat.Pe = *(std::max_element(Pe.begin(), Pe.end()));
}

void make_full(Configuration &c, double* full, double* add)
{
	for (unsigned int j = 0; j <= c.ny; j++) {
		for (unsigned int i = 0; i <= c.nx; i++)
		{
			int l = i + c.offset * j;
			double x = i * c.hx;
			double y = j * c.hy;
			full[l] = -c.grav_x * x - c.grav_y * y + add[l];
		}
	}
};

void Nu_y(Configuration &c, double* f, double &Nu_top, double &Nu_down, double plus_const = 0) {
	Nu_top = 0;
	Nu_down = 0;

	auto dy1 = [&](unsigned int l)
	{
		return abs((f[l + c.offset] - f[l])) / c.hy;
	};
	auto dy1_ = [&](unsigned int l)
	{
		return abs((f[l] - f[l - c.offset])) / c.hy;
	};

	unsigned int l;
	for (unsigned int i = 0; i <= c.nx - 1; i++)
	{
		l = i;
		Nu_down += (dy1(l) + dy1(l + 1)) / 2.0 * c.hx;
		l = i + c.offset * (c.ny - 1);
		Nu_top  += (dy1_(l) + dy1_(l + 1)) / 2.0 * c.hx;
	}

	Nu_down = Nu_down / c.Lx + plus_const;
	Nu_top = Nu_top / c.Lx + plus_const;

}

void Nu_y_for_fixed_flux(Configuration &c, double *f, double &Nu_) {
	Nu_ = 0;
	unsigned int l;
	double s_top = 0, s_bottom = 0;
	for (unsigned int i = 1; i <= c.nx - 1; i++)
	{
		l = i + c.offset * (c.ny - 1);
		s_top += f[l];
		
		l = i;
		s_bottom += f[l];
	}
	s_top = s_top / (c.nx - 1);
	s_bottom = s_bottom / (c.nx - 1);

	Nu_ = 1.0 / (s_bottom - s_top);
}



struct Backup
{
	std::string backup_name[2];
	bool second_backup;
	Backup(std::string name = "recovery", bool second = false)
	{
		backup_name[0] = name + ".txt";
		backup_name[1] = name + "2.txt";
		second_backup = second;
	}

	void save(size_t iter, double time, int call_i, Configuration& c, std::vector<double*> v = {}, std::string head = "")
	{
		stringstream ss, ss2; string str;
		ss.str(""); ss.clear();
		//ss << setprecision(15);
		ss << iter * c.tau;


		std::ofstream write(backup_name[0]);
		write << "x, y, " + head << endl;
		write << "iter,time,call_i,Ra= " << iter << " " << time << " " << call_i << " " << c.Ra << endl;
		write << std::setprecision(10) << std::fixed;
		for (unsigned int j = 0; j <= c.ny; j++) {
			for (unsigned int i = 0; i <= c.nx; i++) {
				unsigned int l = i + c.offset * j;
				write << i * c.hx << " " << j * c.hy << " ";
				for (auto& it : v)
				{
					write << it[l] << " ";
				}
				write << endl;
			}
		}

		if (second_backup)
		{
			std::ofstream write(backup_name[1]);
			write << "x, y, " + head << endl;
			write << "iter,time,call_i,Ra= " << iter << " " << time << " " << call_i << " " << c.Ra << endl;
			write << std::setprecision(10) << std::fixed;
			for (unsigned int j = 0; j <= c.ny; j++) {
				for (unsigned int i = 0; i <= c.nx; i++) {
					unsigned int l = i + c.offset * j;
					write << i * c.hx << " " << j * c.hy << " ";
					for (auto& it : v)
					{
						write << it[l] << " ";
					}
					write << endl;
				}
			}
		}
	};

	int read(size_t &iter, double &time, int &call_i, Configuration& c, std::vector<double*> v = {}, bool read_state = true)
	{
		string str;
		string substr;
		stringstream ss;

		std::ifstream read; 
		if (check2files(backup_name[0], backup_name[1]))
			read.open(backup_name[0]);
		else 
			read.open(backup_name[1]);

		if (!read.is_open())
		{
			std::cout << "no recovery file" << endl;
			return 1;
		}

		getline(read, str); //skip header
		getline(read, str); //read iter, time, Ra

		if (read_state) 
		{
			ss << str;
			ss >> substr; //skip name
			ss >> substr; iter = atoi(substr.c_str());
			ss >> substr; time = atof(substr.c_str());
			ss >> substr; call_i = atoi(substr.c_str());
			ss >> substr; c.Ra = atof(substr.c_str());
		}

		for (unsigned int j = 0; j <= c.ny; j++) {
			for (unsigned int i = 0; i <= c.nx; i++) {
				unsigned int l = i + c.offset * j;

				ss.str(""); 
				ss.clear(); 
				getline(read, str); 
				ss << str;

				ss >> substr; ss >> substr;  //skip reading x,y
				for (auto& it : v)
				{
					ss >> substr;
					it[l] = atof(substr.c_str());
				}
			}
		}
		return 0;
	}

	void read_fields_only(Configuration& c, std::vector<double*> v = {}, unsigned int skip_lines = 0)
	{
		string str;
		string substr;
		stringstream ss;

		std::ifstream read;
		if (check2files(backup_name[0], backup_name[1]))	read.open(backup_name[0]);
		else												read.open(backup_name[1]);

		for (unsigned int q = 0; q < skip_lines; q++)
			getline(read, str);

		for (unsigned int j = 0; j <= c.ny; j++) {
			for (unsigned int i = 0; i <= c.nx; i++) {
				unsigned int l = i + c.offset * j;

				ss.str("");
				ss.clear();
				getline(read, str);
				ss << str;

				ss >> substr; ss >> substr;  //skip reading x,y
				for (auto& it : v)
				{
					ss >> substr;
					it[l] = atof(substr.c_str());
				}
			}
		}
	}


	bool check2files(const std::string& filePath1, const std::string& filePath2) {
		std::ifstream file1(filePath1, std::ios::binary | std::ios::ate); 
		std::ifstream file2(filePath2, std::ios::binary | std::ios::ate);
		if (!file2.is_open()) 
		{
			return true;
		}
		std::streampos size1 = file1.tellg();
		std::streampos size2 = file2.tellg();
		return size1 >= size2; 
	}

};


struct Trajectory
{
	std::vector<double> x, y, z;
	Configuration& c;
	size_t n;
	std::ofstream w;
	Trajectory(Configuration &c_, size_t n_, bool append = false) : c(c_), n(n_)
	{
		x.resize(n);
		y.resize(n);
		z.resize(n);


		if (append) 
		//if (1)
		{
			string str, substr;
			stringstream ss;
			
			str = get_last_line("trajectory.dat"); 
			ss << str;

			for (size_t i = 0; i < n; i++)
			{
				ss >> substr; //skip time
				ss >> substr; x[i] = atof(substr.c_str());
				ss >> substr; y[i] = atof(substr.c_str());
				ss >> substr; z[i] = atof(substr.c_str());
			}
			w.open("trajectory.dat", std::ofstream::app);
		}
		else
		{
			w.open("trajectory.dat");
		}
	}

	std::string get_last_line(const std::string &filename)
	{
		std::ifstream file(filename, std::ios::binary | std::ios::ate);
		if (!file)
			return "";

		std::streamoff pos = file.tellg();
		if (pos == 0)
			return "";

		char ch;
		std::string line;

		pos--;

		// пропустить \n в конце файла
		while (pos >= 0) {
			file.seekg(pos);
			file.get(ch);

			if (ch != '\n' && ch != '\r')
				break;

			pos--;
		}

		while (pos >= 0) {
			file.seekg(pos);
			file.get(ch);

			if (ch == '\n')
				break;

			line.push_back(ch);
			pos--;
		}

		std::reverse(line.begin(), line.end());

		return line;
	}

	void trace_all(double t, double *vx, double *vy, double *vz = nullptr)
	{
		for (size_t p = 0; p < n; p++)
		{
			trace(x[p], y[p], z[p], vx, vy, vz, c);
			w << std::setprecision(12);
			w << t << " " << x[p] << " " << y[p] << " ";
			if (vz != nullptr) w << z[p] << " ";
		}
		w << endl;
	}
	
	void trace(double& x, double& y, double& z, double* vx, double* vy, double* vz, Configuration& c)
	{
		unsigned int i1, j1, k1, i2, j2, k2;
		double VX, VY, VZ, x1, x2, y1, y2, z1, z2, h3;
		double inp;
		h3 = c.dV;
		#define zero 1e-10

		//i1, i2
		if (modf((x / c.hx), &inp) < zero) {
			i1 = unsigned int(x / c.hx);
			i2 = i1 + 1;
			if (i1 == c.nx) {
				i2 = c.nx;
				i1 = i2 - 1;
			}
		}
		else {
			i1 = unsigned int(floor(x / c.hx));
			i2 = unsigned int(ceil(x / c.hx));
		}

		//j1,j2
		if (modf((y / c.hy), &inp) < zero) {
			j1 = unsigned int(y / c.hy);
			j2 = j1 + 1;
			if (j1 == c.ny) {
				j2 = c.ny;
				j1 = j2 - 1;
			}
		}
		else {
			j1 = unsigned int(floor(y / c.hy));
			j2 = unsigned int(ceil(y / c.hy));
		}

		//k1,k2
		if (modf((z / c.hz), &inp) < zero) {
			k1 = unsigned int(z / c.hz);
			k2 = k1 + 1;
			if (k1 == c.nz) {
				k2 = c.nz;
				k1 = k2 - 1;
			}
		}
		else {
			k1 = unsigned int(floor(z / c.hz));
			k2 = unsigned int(ceil(z / c.hz));
		}


		x1 = c.hx * i1; x2 = c.hx * i2;
		y1 = c.hy * j1; y2 = c.hy * j2;
		z1 = c.hz * k1; z2 = c.hz * k2;
		if (vz == nullptr)
		{
			z2 = 1;
			z1 = 0;
		}

		VX = vx[i1 + c.offset * j1 + c.offset2 * k1] * (x2 - x) * (y2 - y) * (z2 - z) +
			vx[i1 + c.offset * j1 + c.offset2 * k2] * (x2 - x) * (y2 - y) * (-z1 + z) +
			vx[i1 + c.offset * j2 + c.offset2 * k1] * (x2 - x) * (-y1 + y) * (z2 - z) +
			vx[i1 + c.offset * j2 + c.offset2 * k2] * (x2 - x) * (-y1 + y) * (-z1 + z) +
			vx[i2 + c.offset * j1 + c.offset2 * k1] * (-x1 + x) * (y2 - y) * (z2 - z) +
			vx[i2 + c.offset * j1 + c.offset2 * k2] * (-x1 + x) * (y2 - y) * (-z1 + z) +
			vx[i2 + c.offset * j2 + c.offset2 * k1] * (-x1 + x) * (-y1 + y) * (z2 - z) +
			vx[i2 + c.offset * j2 + c.offset2 * k2] * (-x1 + x) * (-y1 + y) * (-z1 + z);
		VX = VX / h3;

		VY = vy[i1 + c.offset * j1 + c.offset2 * k1] * (x2 - x) * (y2 - y) * (z2 - z) +
			vy[i1 + c.offset * j1 + c.offset2 * k2] * (x2 - x) * (y2 - y) * (-z1 + z) +
			vy[i1 + c.offset * j2 + c.offset2 * k1] * (x2 - x) * (-y1 + y) * (z2 - z) +
			vy[i1 + c.offset * j2 + c.offset2 * k2] * (x2 - x) * (-y1 + y) * (-z1 + z) +
			vy[i2 + c.offset * j1 + c.offset2 * k1] * (-x1 + x) * (y2 - y) * (z2 - z) +
			vy[i2 + c.offset * j1 + c.offset2 * k2] * (-x1 + x) * (y2 - y) * (-z1 + z) +
			vy[i2 + c.offset * j2 + c.offset2 * k1] * (-x1 + x) * (-y1 + y) * (z2 - z) +
			vy[i2 + c.offset * j2 + c.offset2 * k2] * (-x1 + x) * (-y1 + y) * (-z1 + z);
		VY = VY / h3;
		//cout << VX * c.tau << " " << VY * c.tau << endl;
			 
		x = x + VX * c.tau;
		y = y + VY * c.tau;

		if (vz != nullptr)
		{
			VZ = vz[i1 + c.offset * j1 + c.offset2 * k1] * (x2 - x) * (y2 - y) * (z2 - z) +
				vz[i1 + c.offset * j1 + c.offset2 * k2] * (x2 - x) * (y2 - y) * (-z1 + z) +
				vz[i1 + c.offset * j2 + c.offset2 * k1] * (x2 - x) * (-y1 + y) * (z2 - z) +
				vz[i1 + c.offset * j2 + c.offset2 * k2] * (x2 - x) * (-y1 + y) * (-z1 + z) +
				vz[i2 + c.offset * j1 + c.offset2 * k1] * (-x1 + x) * (y2 - y) * (z2 - z) +
				vz[i2 + c.offset * j1 + c.offset2 * k2] * (-x1 + x) * (y2 - y) * (-z1 + z) +
				vz[i2 + c.offset * j2 + c.offset2 * k1] * (-x1 + x) * (-y1 + y) * (z2 - z) +
				vz[i2 + c.offset * j2 + c.offset2 * k2] * (-x1 + x) * (-y1 + y) * (-z1 + z);
			VZ = VZ / h3;
			z = z + VZ * c.tau;
		}

		if (x < 0 + c.hx) x = x + c.Lx;
		if (x > c.Lx - c.hx) x = x - c.Lx;
		if (y < 0) y = y + c.Ly;
		if (y > c.hy * c.ny) y = y - c.Ly;
		if (z < 0 + c.hz) z = z + c.Lz;
		if (z > c.Lz - c.hz) z = z - c.Lz;
	}
};


