// Project.cpp: 定义控制台应用程序的入口点。
//

#include "stdafx.h"
#include "monitor&dump.h"
#include "drawing.h"
#include "end_condition.h"
#include "interaction.h"
#include "chain.cpp"
#include <conio.h>

string Input::working_path, Input::data_file_name, Input::chain_distr, Input::distance_file_name, Input::gifName;
bitset<5> Input::display_elements;
unsigned Input::canvas_size, Input::EC_MCSMax, Input::EC_w, Input::EC_m, Input::EC_n;
unsigned Input::blockCnt, Input::blockConcent, Input::C_Concent, Input::chain_len, Input::chain_rigid, Input::C_distr;
bool Input::gif, Input::interactions_continuous, Input::screen;
unsigned Input::gifRate, Input::gifStop, Input::shotRate;
unsigned Input::add_time_stamp, Input::add_amount;

typedef enum { ESC = 1, MCS, SAVE, CONTINUE, NIL } CMD;

CMD getCommand() noexcept {
	if (GetAsyncKeyState(VK_ESCAPE) < 0) return ESC;
	if (GetAsyncKeyState(VK_SPACE) < 0) return MCS;
	if (GetAsyncKeyState('S') < 0) return SAVE;
	if (GetAsyncKeyState('C') < 0) return CONTINUE;
	return NIL;
}


int main() {
	printf("------ 高分子相行为模拟 V1.0 By Logic ------\n");
	printf("按任意键加载配置（simulation.ini）...");
	_getch();
	Input ipt("simulation.ini");
	chains_act CA;
	CA.init(ipt);
	sort(chains_act::chains.begin(), chains_act::chains.end());
	Cmonitor monitor; monitor.getDistance();
	Ccanvas canvas(ipt.canvas_size);						
	endCondition end_condition(ipt.EC_MCSMax, ipt.EC_w, ipt.EC_m, ipt.EC_n);	
	Interaction inter("interaction.csv", ipt.interactions_continuous);
	if (canvas.size) {
		canvas.init();
#pragma omp sections
		{
#pragma omp section
			{
				canvas.drawChains(ipt.display_elements);
				canvas.drawSol(ipt.display_elements);
			}
#pragma omp section
			{
				if (Input::gif) canvas.initGIF(ipt.gifName.c_str(), ipt.gifRate, ipt.gifStop);
			}
#pragma omp section
			{
				if (Input::screen) {
					canvas.openScreen();
					canvas.putScreen();
				}
			}
		}
	}
	printf("加载完毕！\n");
	printf("按任意键开始模拟...");
	_getch();
	printf("\n开始模拟！命令说明如下：\n");
	printf("按END     暂停并查看当前帧数（如果不响应请试长按），***暂停后***：\n");
	printf("按s       截图并保存当前状态和数据到文件\n");
	printf("按空格键  运行一个MCS（产生下一帧动画，即动一下）\n");
	printf("按c       继续连续模拟\n");
	printf("按ESC     截图保存并退出\n");
	printf(">");
	unsigned MCSCnt = 0; bool simulation_end = false;
	while (!simulation_end) {
		while (GetAsyncKeyState(VK_END) >= 0)
			if (one_mcs(canvas, end_condition, inter, MCSCnt)) {
				simulation_end = true; break;
			}
		if (!simulation_end) {
			printf("暂停...当前帧数：%d\n>", MCSCnt);
			while (const CMD cmd = getCommand()) {
				if (cmd == CONTINUE) {
					printf("继续...\n>");  break;
				}
				else if (cmd == SAVE) {
					save();
					if (canvas.size) canvas.shotImage(MCSCnt, ipt.working_path);
					printf("保存成功！\n>");
				}
				else if (cmd == MCS) {
					if (one_mcs(canvas, end_condition, inter, MCSCnt)) {
						simulation_end = true; break;
					}
					printf("动一下\n>");
				}
				else if (cmd == ESC) {
					simulation_end = true; break;
				}
			}
		}
			
	}
	printf("保存并退出，保存中...\n");
#pragma omp sections
	{
#pragma omp section
		{
			save();
		}
#pragma omp section
		{
			if (canvas.size) {
				canvas.shotImage(MCSCnt, ipt.working_path);
				if (Input::gif) canvas.endGIF();
				if (Input::screen) canvas.closeScreen();
			}
		}
	}
	return 0;
}