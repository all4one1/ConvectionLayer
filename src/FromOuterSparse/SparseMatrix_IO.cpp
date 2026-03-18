// this source file is for output, print, debug and other auxiliary stuff

#include "SparseMatrix.h"
#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream> 
using std::cout;
using std::endl;
using std::ofstream;



void SparseMatrix::save_compressed_matrix(std::string filename)
{
	std::ofstream w(filename);
	w << Nfull << " " << nval << " " << nrow << endl;

	size_t n = val.size();
	for (size_t i = 0; i < n; i++)
	{
		w << val[i];
		if (i != (n - 1)) w << " ";
		else w << endl;
	}
	n = col.size();
	for (size_t i = 0; i < col.size(); i++)
	{
		w << col[i];
		if (i != (n - 1)) w << " ";
		else w << endl;
	}
	n = row.size();
	for (size_t i = 0; i < row.size(); i++)
	{
		w << row[i];
		if (i != (n - 1)) w << " ";
		else w << endl;
	}
}
void SparseMatrix::save_compressed_matrix_with_rhs(double* b, std::string filename)
{
	save_compressed_matrix(filename);
	std::ofstream w(filename, std::ios_base::app);
	for (int i = 0; i < Nfull; i++)
	{
		w << b[i];
		if (i != (Nfull - 1)) w << " ";
		else w << endl;
	}
}
void SparseMatrix::save_compressed_matrix_binary(std::string filename)
{
	std::ofstream w(filename, std::ios_base::app);
	w.write(reinterpret_cast<const char*>(&Nfull), sizeof(int));
	w.write(reinterpret_cast<const char*>(&nval), sizeof(int));
	w.write(reinterpret_cast<const char*>(&nrow), sizeof(int));

	w.write(reinterpret_cast<const char*>(val.data()), sizeof(double) * nval);
	w.write(reinterpret_cast<const char*>(col.data()), sizeof(int) * nval);
	w.write(reinterpret_cast<const char*>(row.data()), sizeof(int) * nrow);
}
void SparseMatrix::read_compressed_matrix(std::string filename)
{
	std::ifstream read(filename);
	//std::istringstream iss;
	if (!read.good())
	{
		cout << "Error: no such a file '" << filename << "'" << endl;
		return;
	}
	else
	{
		//std::ostringstream oss;
		//oss << read.rdbuf();
		//iss.str(oss.str());
	}
	std::string line;
	//std::string symbol;
	double f;
	int i;
	//getline(iss, line);
	getline(read, line);
	std::stringstream ss;
	ss << line;
	ss >> Nfull; 	resize(Nfull);
	ss >> nval; 
	ss >> nrow;

	//values
	getline(read, line);
	ss.clear();
	ss << line;
	while (ss >> f)
		val.push_back(f);

	//columns
	getline(read, line);
	ss.clear();
	ss << line;
	while (ss >> i)
		col.push_back(i);

	//rows
	getline(read, line);
	ss.clear();
	ss << line;
	//because nrow is already allocated
	for (int r = 0; r < nrow; r++)
	{
		ss >> i;
		row[r] = i;
	}
}
void SparseMatrix::read_compressed_matrix_with_rhs(double* b, std::string filename)
{
	read_compressed_matrix(filename);

	std::ifstream read(filename);
	std::string line;

	getline(read, line);
	getline(read, line);
	getline(read, line);
	getline(read, line);
	getline(read, line);

	std::stringstream ss;
	ss << line;

	for (int i = 0; i < Nfull; i++)
	{
		ss >> b[i];
	}
}
void SparseMatrix::read_compressed_matrix_binary(std::string filename)
{
	std::ifstream r(filename, std::ios::binary);
	r.read(reinterpret_cast<char*>(&Nfull), sizeof(int));
	resize(Nfull);

	r.read(reinterpret_cast<char*>(&nval), sizeof(int));
	r.read(reinterpret_cast<char*>(&nrow), sizeof(int));

	val.resize(nval);
	col.resize(nval);

	r.read(reinterpret_cast<char*>(val.data()), sizeof(double) * nval);
	r.read(reinterpret_cast<char*>(col.data()), sizeof(int) * nval);
	r.read(reinterpret_cast<char*>(row.data()), sizeof(int) * nrow);
}




void SparseMatrix::save_full_matrix(int precision, std::string filename)
{
	ofstream out(filename);
	//ofstream out2("info.dat");
	out << std::fixed << std::setprecision(precision);
	for (unsigned int k = 0; k < row.size() - 1; k++)
	{
		std::vector <double> line(Nfull);
		std::vector <int> t(Nfull);
		for (int i = 0; i < Nfull; i++)
		{
			line[i] = 0.0;
			//	t[i] = -1;
		}

		for (int j = row[k]; j < row[k + 1]; j++)
		{
			line[col[j]] = val[j];
			//	t[col[j]] = type[j];
		}

		for (int i = 0; i < Nfull; i++)
		{
			out << line[i] << " ";
			//	if (t[i] == -1) out2 << "_" << " ";
			//	else out2 << t[i] << " ";
		}
		out << endl;
		//out2 << endl;
	}
}
void SparseMatrix::save_full_matrix_v2(int precision, std::string filename)
{
	ofstream out(filename);
	out << std::fixed << std::setprecision(precision);
	for (int i = 0; i < Nfull; i++) {
		for (int j = 0; j < Nfull; j++)
		{
			out << get_element(i, j) << " ";
		}
		out << endl;
	}
}
void SparseMatrix::save_full_matrix_with_rhs(int precision, double *b, std::string filename)
{
	ofstream out(filename);
	out << std::fixed << std::setprecision(precision);
	for (unsigned int k = 0; k < row.size() - 1; k++)
	{
		std::vector <double> line(Nfull);
		std::vector <int> t(Nfull);
		for (int i = 0; i < Nfull; i++)
		{
			line[i] = 0.0;
		}

		for (int j = row[k]; j < row[k + 1]; j++)
		{
			line[col[j]] = val[j];
		}

		for (int i = 0; i < Nfull; i++)
		{
			out << line[i] << " ";
		}
		out << "	 " << b[k];
		out << endl;
	}
}
void SparseMatrix::read_full_matrix_with_rhs(int N, double* b, std::string filename)
{
	Nfull = N;
	resize(Nfull);
	std::ifstream read(filename);
	if (!read.good())
	{
		cout << "Error: no such a file '" << filename << "'" << endl;
		return;
	}

	for (int i = 0; i < N; i++)
	{
		std::string line;
		std::stringstream ss;
		getline(read, line);
		ss << line;
		double f = 0.0;
		std::map<int, double> m;
		for (int j = 0; j < N; j++)
		{
			ss >> f;
			if (f != 0.0)
			{
				m[j] = f;
			}
		}
		add_line_with_map(m, i);
		if (b != nullptr)
		{
			ss >> b[i];
		}
	}


}



void SparseMatrix::print_storage()
{
	double S = double(val.capacity() * 8 + col.capacity() * 4 + row.capacity() * 4 + diag.capacity() * 8 + type.capacity() * 4);
	int num = int(val.capacity());
	cout << num << " elements, " << S / 1024 / 1024 << " MB approx. matrix memory usage" << " \n\n";
}
void SparseMatrix::print_index_ij(int l)
{
	if (!(l < nval)) { cout << " out of range " << endl; return; }

	int j = col[l];
	int i;
	for (int q = 0; q < Nfull; q++)
	{
		if (l < row[q + 1])
		{
			i = q;
			break;
		}
	}

	cout << "[" << i << "][" << j << "]" << endl;
}
void SparseMatrix::save_type()
{
	ofstream m("type.dat");

	for (int i = 0; i < Nfull; i++) {
		for (int j = 0; j < Nfull; j++)
		{
			std::string str = "_____";
			int t = get_type(i, j);
			if (t == center) str = "centr";
			if (t == south) str = "south";
			if (t == north) str = "north";
			if (t == west) str = "west_";
			if (t == east) str = "east_";

			m << str << " ";
		}
		m << endl;
	}
}
void SparseMatrix::print_compressed_matrix()
{
	//ofstream seq("seq.dat");
	cout << "Sequential (compressed) format: " << endl;
	cout << "val: ";
	for (int i = 0; i < nval; i++)
		cout << val[i] << " ";
	cout << endl;

	cout << "col: ";
	for (int i = 0; i < nval; i++)
		cout << col[i] << " ";
	cout << endl;

	cout << "row: ";
	for (int i = 0; i < Nfull + 1; i++)
		cout << row[i] << " ";
	cout << endl;
}
void SparseMatrix::print_full_matrix()
{
	for (int i = 0; i < Nfull; i++) {
		for (int j = 0; j < Nfull; j++)
		{
			std::cout << get_element(i, j) << " ";
		}
		std::cout << endl;
	}
	cout << "Number of (non-zero) elements: " << nval << endl;
	[this]() {
		int count_zero = 0;
		for (auto it : val)
			if (abs(it) < 1e-15)
				count_zero++;	
		if (count_zero > 0)
			cout << "Number of zero elements: " << count_zero << endl;
	} ();
}

void SparseMatrix::add_submatrix_rightward(SparseMatrix& B, int stride_x, int stride_y)
{
	for (int i = 0; i < B.Nfull; i++)
	{
		int global_i = i + stride_y;
		//int l1 = row[global_i];
		int l2 = row[global_i + 1];

		vecpair v = B.get_a_row(i);

		for (auto& it : v.second)
			it += stride_x;

		insert_many_in_one_row(l2, global_i, v.first, v.second);
	}
}

vecpair SparseMatrix::get_a_row(int ii)
{
	std::vector<int> vi;
	std::vector<double> vd;

	for (int j = row[ii]; j < row[ii + 1]; j++)
	{
		vi.push_back(col[j]);
		vd.push_back(val[j]);
	}
	return std::make_pair(vd, vi);
}




