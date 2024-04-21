#include "Header.h"
#include "Force.h"
#include "Restraint.h"
#include <fstream>

double det2(double a11, double a12, double a21, double a22)
{
	return a11 * a22 - a12 * a21;
}
double det3(double a11, double a12, double a13, double a21, double a22, double a23, double a31, double a32, double a33)
{
	return a11 * det2(a22, a23, a32, a33) - a12 * det2(a21, a23, a31, a33) + a13 * det2(a21, a22, a31, a32);
}

int Read(vector<double3_>& pt_list, vector<vc>& hex_list, Force& F, vector<Restraint>& R, double& E, double& Nu, int& number_of_nodes_of_elem, vector<vector<int>>& list_nodes_with_elem_nums)
{
	string filename = "C:/Users/Igor Volov/Desktop/FC files/Igor7k.fc";
	//string filename = "C:/Users/kolychev.SALDLAB/Desktop/proga Igor/FC files/fidesys11.fc";

	std::ifstream fc_file(filename, std::ios::in);

	if (!fc_file)
	{
		throw std::runtime_error("cannot open fc file: " + filename);
	}
	auto _root = nlohmann::json::parse(fc_file);

	fc_file.close();

	const auto mesh = _root["mesh"];
	const int mesh_node_count = mesh["nodes_count"];
	const int mesh_elems_count = mesh["elems_count"];

	const string tmp1 = _root["mesh"]["elem_types"];
	//cout << tmp1 << endl;
	string num_of_nodes_of_elem0;
	base64_decode(tmp1, num_of_nodes_of_elem0);
	//const string num_of_nodes_of_elem = num_of_nodes_of_elem0;
	const int num_of_nodes_of_elem = *reinterpret_cast<const int*>(num_of_nodes_of_elem0.c_str());
	//cout << "--" << num_of_nodes_of_elem << endl;	//triangles 168430090 //squares 12
	int num_nodes = COUNT_OF_NODES;
	//cout << "num_nodes = " << num_nodes << endl;
	const string tmp = _root["mesh"]["nids"];
	string mesh_nids;
	base64_decode(tmp, mesh_nids);

	const string tmp3 = _root["mesh"]["nodes"];
	string mesh_nodes_tmp;
	base64_decode(tmp3, mesh_nodes_tmp);
	const string mesh_nodes = mesh_nodes_tmp;

	pt_list.resize(mesh_node_count);
	list_nodes_with_elem_nums.resize(mesh_node_count);

	//! Read Cubit numeration
	std::map<int, int> _map_node_numeration;
	for (int node_ID = 0; node_ID < mesh_node_count; node_ID++) {
		const int mesh_node_ID = *reinterpret_cast<const int*>(mesh_nids.c_str() + (node_ID) * sizeof(int));
		_map_node_numeration.insert(std::pair<int, int>(mesh_node_ID, node_ID));
		pt_list[node_ID].x = *reinterpret_cast<const double*>(mesh_nodes.c_str() + (node_ID * 3 + 0) * sizeof(double));
		pt_list[node_ID].y = *reinterpret_cast<const double*>(mesh_nodes.c_str() + (node_ID * 3 + 1) * sizeof(double));
		pt_list[node_ID].z = *reinterpret_cast<const double*>(mesh_nodes.c_str() + (node_ID * 3 + 2) * sizeof(double));
	}

	const string tmp2 = _root["mesh"]["elems"];
	string mesh_elems;
	base64_decode(tmp2, mesh_elems);
	hex_list.resize(mesh_elems_count);

	for (int element_ID = 0; element_ID < mesh_elems_count; element_ID++) {
		for (int j = 0; j < num_nodes; j++)
		{
			const int node_number = *reinterpret_cast<const int*>(mesh_elems.c_str() + (element_ID * COUNT_OF_NODES + j) * sizeof(int)); // [elem][j].asInt();
			map<int, int>::iterator map_iterator;
			map_iterator = _map_node_numeration.find(node_number);
			//const int glob_node = FCFindInMap(_map_node_numeration, node_number, "mesh.elems");
			hex_list[element_ID].n[j] = map_iterator->second;
		}
	}

	for (int j = 0; j < mesh_elems_count; j++)
	{
		for (int k = 0; k < num_nodes; k++)
		{
			list_nodes_with_elem_nums[hex_list[j].n[k]].push_back(j);
		}
	}

	const string Jung0 = _root["materials"][0]["elasticity"][0]["constants"][0];
	string Jung1;
	base64_decode(Jung0, Jung1);
	double Jung = *reinterpret_cast<const double*>(Jung1.c_str());
	const string Poison0 = _root["materials"][0]["elasticity"][0]["constants"][1];
	string Poison1;
	base64_decode(Poison0, Poison1);
	double Poison = *reinterpret_cast<const double*>(Poison1.c_str());
	E = Jung; Nu = Poison;

	int len = 6;
	
	const int apply_to_size = _root["loads"][0]["apply_to_size"];

	vector<int> v(apply_to_size);
	const string Force0 = _root["loads"][0]["apply_to"];
	string Force1;
	base64_decode(Force0, Force1);
	for (int i = 0; i < apply_to_size; i++)
	{
		v[i] = *reinterpret_cast<const int*>(Force1.c_str() + i * sizeof(int));
	}
	//pt_list[node_ID].x = *reinterpret_cast<const double*>(mesh_nodes.c_str() + (node_ID * 3 + 0) * sizeof(double));


	vector<string> data0(len), data1(len);
	vector<double> data(len);
	for (int i = 0; i < data0.size(); i++)
	{
		data0[i] = _root["loads"][0]["data"][i];
		base64_decode(data0[i], data1[i]);
		data[i] = *reinterpret_cast<const double*>(data1[i].c_str());
		//cout << "data[" << i << "] = " << data[i] << endl;
	}
	double3_ f; f.x = data[0]; f.y = data[1]; f.z = data[2];

	F = Force(f, v);
	
	//vector<string> data0(len), data1(len);
	//vector<double> data(len);
	for (int k = 0; k < _root["restraints"].size(); k++)
	{
		const int apply_to_size_restraints = _root["restraints"][k]["apply_to_size"];
		vector<int> flag;
		for (int i = 0; i < _root["restraints"][k]["flag"].size(); i++)
		{
			flag.push_back(_root["restraints"][k]["flag"][i]);
		}
		//for (int i = 0; i < flag.size(); i++)
		//{
			//cout << flag[i] << " ";
		//}
		//cout << endl;

		vector<int> w(apply_to_size_restraints);
		const string Restraints0 = _root["restraints"][k]["apply_to"];
		string Restraints1;
		base64_decode(Restraints0, Restraints1);
		for (int i = 0; i < apply_to_size_restraints; i++)
		{
			w[i] = *reinterpret_cast<const int*>(Restraints1.c_str() + i * sizeof(int));
			//cout << "w[" << i << "] = " << w[i] << endl;
		}

		data0.clear(); data1.clear(); data.clear();
		data0 = data1 = vector<string>(len);
		data = vector<double>(len);
		for (int i = 0; i < data0.size(); i++)
		{
			data0[i] = _root["restraints"][k]["data"][i];
			base64_decode(data0[i], data1[i]);
			data[i] = *reinterpret_cast<const double*>(data1[i].c_str());
			//cout << "data[" << i << "] = " << data[i] << endl;
		}

		double3_ r; r.x = data[0]; r.y = data[1]; r.z = data[2];
		R.push_back(Restraint(r, w, flag));
	}

	return 0;
}