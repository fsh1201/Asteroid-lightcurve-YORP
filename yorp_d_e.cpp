#include <iostream>
#include <fstream>
#include <sstream>

using namespace std;

int main(int argc, char** argv)
{
	if (argc < 4)
	{
		cout << "yorp_d_e path lcpath parameter" << endl;
		cout << "path: the path where you want to save the results" << endl;
		cout << "lcpath: the path of light curve file" << endl;
		cout << "parameter: the path of input parameters" << endl;
		exit(1);
	}
	string path(argv[1]);
	string lcpath(argv[2]);
	string parpath(argv[3]);
	int nbs, nth, dplc;

	ifstream par(parpath, ios::in);
	string line;
	for (int i = 0; i < 17; i++)
		getline(par, line);
	par >> nbs;
	getline(par, line);
	par >> nth;
	getline(par, line);
	par >> dplc;
	getline(par, line);
	par.close();

	string cmd;

	cmd = "cp " + lcpath + " " + path + "/";
	cout << cmd << endl;
	system(cmd.c_str());

	if (dplc == 0 || dplc == 1)
	{
		cmd = "rm -rf " + path + "/inverse " + path + "/bootstrap" + to_string(dplc);
		cout << cmd << endl;
		system(cmd.c_str());

		cmd = "mkdir " + path + "/inverse " + path + "/bootstrap" + to_string(dplc) + " " + path + "/inverse/in " + path + "/inverse/out "
			+ path + "/bootstrap" + to_string(dplc) + "/in " + path + "/bootstrap" + to_string(dplc) + "/out "
			+ path + "/bootstrap" + to_string(dplc) + "/lcbs " + path + "/inverse/out/best ";
		cout << cmd << endl;
		system(cmd.c_str());
	}
	if (dplc == 2)
	{
		cmd = "rm -rf " + path + "/inverse " + path + "/bootstrap0 " + path + "/bootstrap1";
		cout << cmd << endl;
		system(cmd.c_str());

		cmd = "mkdir " + path + "/inverse " + path + "/bootstrap0 " + path + "/bootstrap1 " + path + "/inverse/in " + path + "/inverse/out "
			+ path + "/bootstrap0/in " + path + "/bootstrap0/out " + path + "/bootstrap1/in " + path + "/bootstrap1/out "
			+ path + "/bootstrap0/lcbs " + path + "/bootstrap1/lcbs " + path + "/inverse/out/best ";
		cout << cmd << endl;
		system(cmd.c_str());
	}

	cmd = "cp " + parpath + " " + path + "/inverse/in/in";
	cout << cmd << endl;
	system(cmd.c_str());

	cmd = "./yorp " + lcpath + " " + path + "/inverse/in/in " + path + "/inverse/out/out_area " + path + "/inverse/out/out_par "
		+ path + "/inverse/out/out_lcs " + path + "/inverse/out/yorp_chi2.txt";
	cout << cmd << endl;
	system(cmd.c_str());

	if (dplc == 0)
	{
		cmd = "python ./ychi2bsin.py " + path + "/inverse/in/in " + path + "/inverse/out/yorp_chi2.txt.sort.txt "
			+ path + "/bootstrap0/in/in";
		cout << cmd << endl;
		system(cmd.c_str());

		cmd = "cp " + lcpath + " " + path + "/bootstrap0/lcbs/lc.txt";
		cout << cmd << endl;
		system(cmd.c_str());

		cmd = "./yorp_e " + path + "/bootstrap0/lcbs/lc.txt" + " " + path + "/bootstrap0/in/in " + path + "/bootstrap0/out/out_area "
			+ path + "/bootstrap0/out/out_par " + path + "/bootstrap0/out/out_lcs " + path + "/bootstrap0/out/yorp_chi2.txt "
			+ to_string(nbs) + " " + to_string(nth) + " " + to_string(dplc);
		cout << cmd << endl;
		system(cmd.c_str());

		cmd = "python ./fig1st.py " + path + "/bootstrap0/results " + path + "/bootstrap0/out/yorp_chi2.txt";
		cout << cmd << endl;
		system(cmd.c_str());
	}

	if (dplc == 1)
	{
		cmd = "python ./ychi2bsin.py " + path + "/inverse/in/in " + path + "/inverse/out/yorp_chi2.txt.sort.txt "
			+ path + "/bootstrap1/in/in";
		cout << cmd << endl;
		system(cmd.c_str());

		cmd = "cp " + lcpath + " " + path + "/bootstrap1/lcbs/lc.txt";
		cout << cmd << endl;
		system(cmd.c_str());

		cmd = "./yorp_e " + path + "/bootstrap1/lcbs/lc.txt" + " " + path + "/bootstrap1/in/in " + path + "/bootstrap1/out/out_area "
			+ path + "/bootstrap1/out/out_par " + path + "/bootstrap1/out/out_lcs " + path + "/bootstrap1/out/yorp_chi2.txt "
			+ to_string(nbs) + " " + to_string(nth) + " " + to_string(dplc);
		cout << cmd << endl;
		system(cmd.c_str());

		cmd = "python ./fig1st.py " + path + "/bootstrap1/results " + path + "/bootstrap1/out/yorp_chi2.txt";
		cout << cmd << endl;
		system(cmd.c_str());
	}

	if (dplc == 2)
	{
		cmd = "python ./ychi2bsin.py " + path + "/inverse/in/in " + path + "/inverse/out/yorp_chi2.txt.sort.txt "
			+ path + "/bootstrap0/in/in";
		cout << cmd << endl;
		system(cmd.c_str());

		cmd = "cp " + lcpath + " " + path + "/bootstrap0/lcbs/lc.txt";
		cout << cmd << endl;
		system(cmd.c_str());

		cmd = "./yorp_e " + path + "/bootstrap0/lcbs/lc.txt" + " " + path + "/bootstrap0/in/in " + path + "/bootstrap0/out/out_area "
			+ path + "/bootstrap0/out/out_par " + path + "/bootstrap0/out/out_lcs " + path + "/bootstrap0/out/yorp_chi2.txt "
			+ to_string(nbs) + " " + to_string(nth) + " 0";
		cout << cmd << endl;
		system(cmd.c_str());

		cmd = "python ./fig1st.py " + path + "/bootstrap0/results " + path + "/bootstrap0/out/yorp_chi2.txt";
		cout << cmd << endl;
		system(cmd.c_str());


		cmd = "python ./ychi2bsin.py " + path + "/inverse/in/in " + path + "/inverse/out/yorp_chi2.txt.sort.txt "
			+ path + "/bootstrap1/in/in";
		cout << cmd << endl;
		system(cmd.c_str());

		cmd = "cp " + lcpath + " " + path + "/bootstrap1/lcbs/lc.txt";
		cout << cmd << endl;
		system(cmd.c_str());

		cmd = "./yorp_e " + path + "/bootstrap1/lcbs/lc.txt" + " " + path + "/bootstrap1/in/in " + path + "/bootstrap1/out/out_area "
			+ path + "/bootstrap1/out/out_par " + path + "/bootstrap1/out/out_lcs " + path + "/bootstrap1/out/yorp_chi2.txt "
			+ to_string(nbs) + " " + to_string(nth) + " 1";
		cout << cmd << endl;
		system(cmd.c_str());

		cmd = "python ./fig1st.py " + path + "/bootstrap1/results " + path + "/bootstrap1/out/yorp_chi2.txt";
		cout << cmd << endl;
		system(cmd.c_str());
	}

	ifstream result(path + "/inverse/out/yorp_chi2.txt.sort.txt", ios::in);
	getline(result, line);
	getline(result, line);
	stringstream ssline(line);
	double t;
	int nt, i;
	for (int j = 0; j < 12; j++)
	{
		ssline >> t;
	}
	ssline >> nt >> i;
	cmd = "cp " + path + "/inverse/out/out_area_" + to_string(nt) + "_" + to_string(i) + ".txt " + path + "/inverse/out/best/";
	cout << cmd << endl;
	system(cmd.c_str());
	cmd = "cp " + path + "/inverse/out/out_par_" + to_string(nt) + "_" + to_string(i) + ".txt " + path + "/inverse/out/best/";
	cout << cmd << endl;
	system(cmd.c_str());
	cmd = "cp " + path + "/inverse/out/out_lcs_" + to_string(nt) + "_" + to_string(i) + ".txt " + path + "/inverse/out/best/";
	cout << cmd << endl;
	system(cmd.c_str());
	cmd = "cat " + path + "/inverse/out/best/out_area_" + to_string(nt) + "_" + to_string(i)
		+ ".txt | ./minkowski | ./standardtri > " + path + "/inverse/out/best/shape.txt";
	cout << cmd << endl;
	system(cmd.c_str());
	cmd = "cat " + path + "/inverse/out/best/shape.txt | ./shape2obj > " + path + "/inverse/out/best/shape.obj";
	cout << cmd << endl;
	system(cmd.c_str());

	return 0;
}