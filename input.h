// inline

#pragma once
#include <string>
#include <graphics.h>
#include <vector>
#include <bitset>
#include "SimpleIni.h"
#include <Windows.h>
using namespace std;

struct Input {
	static string working_path;
	static bitset<5> display_elements;		// show_bond, show_D, show_C, show_B, show_A
	static unsigned canvas_size;
	static unsigned EC_MCSMax, EC_w, EC_m, EC_n;
	static unsigned blockCnt, blockConcent, C_Concent, chain_len, chain_rigid, C_distr;
	static string chain_distr;
	static string gifName;
	static unsigned gifRate, gifStop, shotRate;
	static string data_file_name, distance_file_name;
	static bool interactions_continuous, gif;
	static bool screen;
	static unsigned add_time_stamp, add_amount;
	Input() noexcept {}
	Input(const char* iniName) {
		CSimpleIniCaseA ini;
		ini.LoadFile(iniName);
		set_working_path(ini.GetValue("path", "working_path"));
		CopyFile(iniName, (working_path + iniName).c_str(), false);
#pragma omp sections
		{
#pragma omp section
			{
				set_polymer(
					ini.GetLongValue("polymer", "block_count"),
					ini.GetLongValue("polymer", "block_concentration"),
					ini.GetLongValue("polymer", "chain_length"),
					ini.GetValue("polymer", "block_distribution"),
					ini.GetLongValue("polymer", "rigidity")
				);
			}
#pragma omp section
			{
				set_C(
					ini.GetLongValue("solution", "C_concentration"),
					ini.GetLongValue("solution", "C_distribution")
				);
				set_interactions(
					ini.GetBoolValue("interaction", "is_continuous")
				);
				set_dump(
					ini.GetValue("dump", "data_file_name"),
					ini.GetValue("dump", "distance_file_name")
				);
			}
#pragma omp section
			{
				set_add_C(
					ini.GetLongValue("add_C", "time_stamp"),
					ini.GetLongValue("add_C", "amount")
				);
			}
#pragma omp section
			{
				set_EC(
					ini.GetLongValue("end_condition", "MCSMax"),
					ini.GetLongValue("end_condition", "which_one"),
					ini.GetLongValue("end_condition", "mol_percent"),
					ini.GetLongValue("end_condition", "no_move_MCS")
				);
			}
#pragma omp section
			{
				set_pic(
					ini.GetValue("pic", "gif_name"),
					ini.GetLongValue("pic", "gif_rate"),
					ini.GetLongValue("pic", "gif_stoppage"),
					ini.GetLongValue("pic", "shot_bmp_rate")
				);
			}
#pragma omp section
			{
				set_display(
					ini.GetBoolValue("display", "put_screen"),
					ini.GetValue("display", "display_elements"),
					ini.GetLongValue("display", "canvas_size")
				);
			}
		}
	}
	void set_working_path(const char* ipt) {
		working_path = ipt;
	}
	// block_count, block_concentration, chain_length; block_distribution, rigidity
	void set_polymer(unsigned ipt_1, unsigned ipt_2, unsigned ipt_3, const char* ipt_4, unsigned ipt_5) {
		blockCnt = ipt_1; blockConcent = ipt_2; chain_len = ipt_3; chain_distr = ipt_4; chain_rigid = ipt_5;
	}
	void set_add_C(unsigned ipt_1, unsigned ipt_2) {
		add_time_stamp = ipt_1; add_amount = ipt_2;
	}
	// C_concentration, C_distribution_mode: 0 - UP 1 - DOWN 2 - RANDOM
	void set_C(unsigned ipt_1, unsigned ipt_2) noexcept {
		C_Concent = ipt_1; C_distr = ipt_2;
	}
	// partial_silent, display_elements, canvas_size
	void set_display(bool ipt_1, const char* ipt_2, unsigned ipt_3) {
		screen = ipt_1;
		display_elements = bitset<5>(ipt_2);
		canvas_size = ipt_3;
	}
	// Interactions: file_name, is_continuous
	void set_interactions(bool ipt_1) {
		interactions_continuous = ipt_1;
		CopyFile("interaction.csv", (working_path + "interaction.csv").c_str(), false);
	}
	// end_condition: MCSMax, which_one, mol_percent, no_move_MCS
	void set_EC(unsigned ipt_1, unsigned ipt_2, unsigned ipt_3, unsigned ipt_4) noexcept {
		EC_MCSMax = ipt_1; EC_w = ipt_2; EC_m = ipt_3; EC_n = ipt_4;
	}
	// gif_name, gif_rate(frame per second), gif_stoppage(frames), shot_bmp_rate(by MCS)
	void set_pic(const string& ipt_1, unsigned ipt_2, unsigned ipt_3, unsigned ipt_4) {
		gif = (ipt_1 != "");
		gifName = working_path + ipt_1; gifRate = ipt_2; gifStop = ipt_3; shotRate = ipt_4;
	}
	// data_file_name, distance_file_name
	void set_dump(const string& ipt_1, const string& ipt_2) {
		data_file_name = working_path + ipt_1;
		distance_file_name = working_path + ipt_2;
	}
};