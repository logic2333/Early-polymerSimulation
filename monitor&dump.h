// inline

#pragma once
#include "chain.h"
#include <fstream>
#include <omp.h>
#include <iomanip>
#include <bitset>
using namespace std;


struct Cmonitor {
	static vector<unsigned> movage;
	static vector<double> distances;
	static void getDistance() {
		unsigned sum = 0;
		const int chainCnt = chains_act::chain_heads.size();
#pragma omp parallel for reduction(+:sum)
		for (int i = 0; i < chainCnt; i++)
			sum += (chains_act::chain_heads[i] - chains_act::chain_tails[i]).len();
		distances.push_back(double(sum) / chainCnt);
	}
	static vector<double> getMovage() {
		vector<double> res(5);	// 0-A, 1-B, 2-C, 3-A+B, 4-A+B+C
		for (unsigned i = 0; i < 3; i++)
			if (chains_act::molCnt[i])
				res[i] = double(movage[i]) / chains_act::molCnt[i];
		res[3] = double(movage[0] + movage[1]) / (chains_act::molCnt[0] + chains_act::molCnt[1]);
		res[4] = double(movage[0] + movage[1] + movage[2]) / (chains_act::molCnt[0] + chains_act::molCnt[1] + chains_act::molCnt[2]);
		return res;
	}
};

vector<unsigned> Cmonitor::movage;
vector<double> Cmonitor::distances;

struct dump {
	ofstream data_file, distance_file;
	dump() noexcept {}
	dump(const string& datafileName, const string& distancefileName) {
		data_file.open(datafileName.c_str(), ios::binary);
		distance_file.open(distancefileName.c_str());
	}
	void write_cubic() {
		if (data_file) {
			const int a = chains_act::all.size();
			data_file.write((char*)&coord_3D::_r, sizeof(unsigned short));
			auto data = new bitset<2 * MAXSPACE>();
#pragma omp parallel for
			for (int i = 0; i < a; i++)
				switch (chains_act::all[i].occup) {
				case 'A': (*data)[2 * i] = 1; (*data)[2 * i + 1] = 1; break;
				case 'B': (*data)[2 * i + 1] = 1; break;
				case 'C': (*data)[2 * i] = 1; break;
				case 'D': break;
				}
			data_file.write((char*)data, ceil(a / 4.0));
			delete data;
		}
	}
	void write_surroundings() {
		if (data_file) {
			const int a = chains_act::all.size();
			unsigned short* data = new unsigned short[a];
#pragma omp parallel for
			for (int i = 0; i < a; i++)
				data[i] = chains_act::all[i].getSurroundings(st, true);
			data_file.write((char*)data, sizeof(unsigned short) * a);
			delete[] data;
		}
	}
	void write_energy() {
		if (data_file) {
			const int a = chains_act::all.size();
			double* data = new double[a];
#pragma omp parallel for
			for (int i = 0; i < a; i++) {
				vector<double> surroundings;
				chains_act::all[i].getSurroundings(st, surroundings);
				data[i] = calcEnergy(chains_act::all[i].occup, surroundings);
			}
			data_file.write((char*)data, sizeof(double) * a);
			delete[] data;
		}
	}
	void write_distance() {
		if (distance_file) {
			distance_file << "立方体边长：" << coord_3D::_r << ' '
				<< "链数：" << chains_act::chain_heads.size() << ' '
				<< "链长：" << Input::chain_len << ' '
				<< "刚度：" << Input::chain_rigid << ' '
				<< "记录数：" << Cmonitor::distances.size() << ' '
				<< "采样率：" << (Input::EC_MCSMax < 10000 ? 1 : Input::EC_MCSMax / 10000) << endl;
			distance_file << fixed << setprecision(2);
			for (auto& i : Cmonitor::distances)
				distance_file << i << endl;
		}
	}
	void end() {
		if (data_file.is_open()) data_file.close();
		if (distance_file.is_open()) distance_file.close();
	}
};

void save() {
	dump dumper(Input::data_file_name, Input::distance_file_name);
	dumper.write_cubic();
	dumper.write_surroundings();
	dumper.write_energy();
	dumper.write_distance();
	dumper.end();
}