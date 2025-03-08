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

	cmd = "rm -rf " + path + "/inverse " + path + "/bootstrap";
	cout << cmd << endl;
	system(cmd.c_str());

	cmd = "mkdir " + path + "/inverse " + path + "/bootstrap " + path + "/inverse/in " + path + "/inverse/out "
		+ path + "/bootstrap/in " + path + "/bootstrap/out " + path + "/bootstrap/lcbs "
		+ path + "/inverse/out/best ";
	cout << cmd << endl;
	system(cmd.c_str());

	cmd = "cp " + parpath + " " + path + "/inverse/in/in";
	cout << cmd << endl;
	system(cmd.c_str());

	cmd = "./yorp " + lcpath + " " + path + "/inverse/in/in " + path + "/inverse/out/out_area " + path + "/inverse/out/out_par "
		+ path + "/inverse/out/out_lcs " + path + "/inverse/out/yorp_chi2.txt";
	cout << cmd << endl;
	system(cmd.c_str());

	cmd = "python ./ychi2bsin.py " + path + "/inverse/in/in " + path + "/inverse/out/yorp_chi2.txt.sort.txt "
		+ path + "/bootstrap/in/in";
	cout << cmd << endl;
	system(cmd.c_str());

	cmd = "cp " + lcpath + " " + path + "/bootstrap/lcbs/lc.txt";
	cout << cmd << endl;
	system(cmd.c_str());

	cmd = "./yorp_e " + path + "/bootstrap/lcbs/lc.txt" + " " + path + "/bootstrap/in/in " + path + "/bootstrap/out/out_area "
		+ path + "/bootstrap/out/out_par " + path + "/bootstrap/out/out_lcs " + path + "/bootstrap/out/yorp_chi2.txt "
		+ to_string(nbs) + " " + to_string(nth) + " " + to_string(dplc);
	cout << cmd << endl;
	system(cmd.c_str());

	cmd = "python ./fig1st.py " + path + "/bootstrap/results " + path + "/bootstrap/out/yorp_chi2.txt";
	cout << cmd << endl;
	system(cmd.c_str());

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