#include "mpi.h"
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <cstdlib>
#include <map>
#include <iomanip>
#include <cmath>
#include <algorithm>
#include <set>
#include "mstate.h"
using namespace std;

// store conf id to conf name and number of e and H


void help_message();
int load_ms(vector <mstate> &all_states);
int write_ms(vector <mstate> &states);
int merge_ms(vector <mstate> &red, vector <mstate> &all);
int merge_crg(vector <mstate> &red, vector <mstate> &all);
int output_ms_out(vector <mstate> &states);
int output_ms_crg(vector <mstate> &states, const int &argc, char *argv[]);
int load_sub_binary(vector <mstate> &states, char *sub_name);
int write_sub_binary(vector <mstate> &states, char *sub_name);
int only_key_res(vector <mstate> &states, char *file_name);

void out_id(map < unsigned short, pair < string, pair <int, int> > > dic);
bool sortState(const mstate &lhs, const mstate &rhs) {
	return lhs.count < rhs.count;
}
int main(int argc, char *argv[])
{
	time_t time_start = time(NULL);
	// read "head3.lst"
	mstate::load_id2conf();

	//out_id(mstate::id2conf);
	vector <mstate> all_states;
	load_ms(all_states);
	int n_states = all_states.size();

	int n_samp = 0;

	for (int i=0; i<argc; i++) {
		string flag(argv[i]);
		if (flag == "-h") {
			help_message();
			return 0;
		}
		else if (flag == "-f") {
			only_key_res(all_states, argv[i+1]);
		}
	}



	int (*merge_f)(vector <mstate> &red, vector <mstate> &all);
	if (argc == 2 && string(argv[1]) == "-c")
	{
		merge_f = merge_crg;
	}
	else {
		merge_f = merge_ms;
	}


	cout << n_states << " micro states loaded" << endl;

	set<mstate> myset;
	for (int i=0; i<n_states; i++) {
		set<mstate>::iterator setIt = myset.find(all_states[i]);
		if (setIt != myset.end()) {
			mstate newState(all_states[i]);
			newState.count += setIt->count;
			newState.H += setIt->H;
			newState.Hsq += setIt->Hsq;
			myset.erase(setIt);
			myset.insert(newState);
		}
		else {
			myset.insert(all_states[i]);
		}
	}
	cout << "Number of unique states: " << myset.size() << endl;

	vector<mstate> red_states;
	red_states.insert(red_states.begin(), myset.begin(), myset.end());
    sort(red_states.begin(), red_states.end(), sortState);
        
    for (int i=0; i<argc; i++) {
    	string flag(argv[i]);
    	if (flag == "-s") {
    		output_ms_out(red_states);
    	}
    }

	write_ms(red_states);
	output_ms_crg(red_states, argc, argv);
	cout << "Total time: " << time(NULL) - time_start << endl;

	return 0;
}



int load_ms(vector <mstate> &all_states)
{
    int n_spe;
    unsigned short conf;
    char buffer[9];

    vector <string> resName;
    vector <unsigned short> iconf;

    mstate samp_state;
      
    ifstream ifile("ms.dat", ios::in | ios::binary);
    ifile.read((char *) &n_spe, 4);
    for (int i_spe=0; i_spe<n_spe; i_spe++)
    {
        ifile.read(buffer, 8); 
        buffer[8] = '\0';
        resName.push_back(string(buffer));
    }
    mstate::n_spe = n_spe;
    mstate::head_spe = resName;


    while (ifile.read((char *) &conf, 2))
    {
    	iconf.push_back(conf);
        for (int i_spe=1; i_spe<n_spe; i_spe++)
        {
            ifile.read((char *) &conf, 2);
            iconf.push_back(conf);
        }
        ifile.read((char *) &samp_state.H, 8);
        ifile.read((char *) &samp_state.Hsq, 8);
        ifile.read((char *) &samp_state.count, 4);

        samp_state.iconf = iconf;
        iconf.clear();



        all_states.push_back(samp_state);
    }
    
    return 0;
}

int merge_ms(vector <mstate> &red, vector <mstate> &all)
{
    for (vector <mstate>::iterator i_all=all.begin(); i_all<all.end(); i_all++)
    {
        vector <mstate>::iterator i_red;
        for (i_red=red.begin(); i_red<red.end(); i_red++)
        {
            if (*i_red == *i_all)
            {
                (*i_red).merge(*i_all);
                break;
            }
        }
        if (i_red == red.end())
        {
            red.push_back(*i_all);
        }
    }

    return 0;
}
                  
int merge_crg(vector <mstate> &red, vector <mstate> &all)
{
    for (vector <mstate>::iterator i_all=all.begin(); i_all<all.end(); i_all++)
    {
        vector <mstate>::iterator i_red;
        for (i_red=red.begin(); i_red<red.end(); i_red++)
        {
            if ((*i_red).same_crg(*i_all))
            {
                (*i_red).merge(*i_all);
                break;
            }
        }
        if (i_red == red.end())
        {
            red.push_back(*i_all);
        }
    }

    return 0;
}

int output_ms_out(vector <mstate> &states)
{
    ofstream ofile("ms_out", ios::out);

    ofile << setiosflags(ios::right | ios::fixed);
    for (int i_spe=0; i_spe<mstate::n_spe; i_spe++)
    {
       ofile << setw(15) << mstate::head_spe[i_spe];
    }
    ofile << setw(6) << "H"
          << setw(4) << "e"
          << setw(5) << "crg"
          << setw(10) << "H_avg"
          << setw(10) << "H_std"
          << setw(10) << "occ" 
          << setw(10) << "count"; 
    ofile << endl;

    //cout << states.size() << " uniqe states of the specified residues" << endl;

    int n_samp=0;
    for (size_t j=0; j<states.size(); j++)
    {
        n_samp += states[j].count;
    }

    //out_id(mstate::id2conf);
    for (size_t j=0; j<states.size(); j++)
    {
        states[j].stat();
        for (int i_spe=0; i_spe<mstate::n_spe; i_spe++)
        {
         //   cout << mstate::id2conf[1].first << '\t' << mstate::id2conf[states[j].iconf[i_spe]].first << endl;
            ofile << setw(15) << mstate::id2conf[states[j].iconf[i_spe]].first;
        }
        ofile << setw(6) << states[j].p
              << setw(4) << states[j].e
              << setw(5) << states[j].crg
              << setw(10) << setprecision(2) << states[j].H_avg
              << setw(10) << setprecision(2) << states[j].H_std
              << setw(9) << setprecision(2) << (float) states[j].count * 100 /n_samp << "%"
              << setw(10) << states[j].count; 
        ofile << endl;
    }
   

    ofile.close();

    return 0;
} 
    
int output_ms_crg(vector <mstate> &states, const int &argc, char *argv[])
{
    vector <mstate> crg_states;
    bool same = true;

    if (!(argc == 2 && string(argv[1]) == "-c")) {
        merge_crg(crg_states, states);
    }
    else {
        crg_states = states;
    }

    ofstream ofile("ms_crg", ios::out);

    ofile << setiosflags(ios::right | ios::fixed);
    for (int i_spe=0; i_spe<mstate::n_spe; i_spe++)
    {
       ofile << setw(9) << mstate::head_spe[i_spe];
    }
    ofile << setw(6) << "H"
          << setw(4) << "e"
          << setw(5) << "crg"
          << setw(10) << "H_avg"
          << setw(10) << "H_std"
          << setw(10) << "occ"
          << setw(10) << "count";
    ofile << endl;

    //cout << crg_states.size() << " uniqe charge states loaded" << endl;

    int n_samp=0;
    for (size_t j=0; j<states.size(); j++)
    {
        n_samp += states[j].count;
    }

    for (size_t j=0; j<crg_states.size(); j++)
    {
        crg_states[j].stat();
        for (int i_spe=0; i_spe<mstate::n_spe; i_spe++)
        {

            ofile << setw(9) 
                  << mstate::id2conf[crg_states[j].iconf[i_spe]].second.second
                   - mstate::id2conf[crg_states[j].iconf[i_spe]].second.first;
        }
        ofile << setw(6) << crg_states[j].p
              << setw(4) << crg_states[j].e
              << setw(5) << crg_states[j].crg
              << setw(10) << setprecision(2) << crg_states[j].H_avg
              << setw(10) << setprecision(2) << crg_states[j].H_std
              << setw(9) << setprecision(2) << (float)crg_states[j].count * 100 /n_samp << "%"
              << setw(10) << crg_states[j].count;
        ofile << endl;
    }


    ofile.close();

    return 0;
}

int write_sub_binary(vector <mstate> &sub_states, char *sub_name)
{
    ofstream ofile(sub_name, ios::out | ios::binary);

    unsigned int s_size=sub_states.size();
    ofile.write((char *) &s_size, sizeof(unsigned int));

    if (s_size > 0)
    {
        for (size_t i=0; i<sub_states.size(); i++)
        {
            sub_states[i].write_binary(ofile);
        }
     }

     ofile.close();

    return 0;
}

int load_sub_binary(vector <mstate> &sub_states, char *sub_name)
{
    ifstream ifile(sub_name, ios::in | ios::binary);
    unsigned int s_size;

    ifile.read((char *) &s_size, sizeof(unsigned int));

    if (s_size > 0)
    {
        for (size_t i=0; i<s_size; i++)
        {
            mstate buff_state;
            buff_state.load_binary(ifile);
            sub_states.push_back(buff_state);
        }
    }

    ifile.close();

    return 0;
}

int write_ms(vector <mstate> &states)
{
    ofstream ofile("bak_ms.dat", ios::out | ios::binary);

    int n_spe = mstate::n_spe;
    vector <string> head(mstate::head_spe);

    ofile.write((char *) &n_spe, sizeof(int));
    for (int i_spe=0; i_spe<n_spe; i_spe++) {
        ofile.write((char *) head[i_spe].c_str(), 8);
    }
    
    for (int i=0; i<states.size(); i++) {
        for (int i_spe=0; i_spe<n_spe; i_spe++) {
            ofile.write((char *) &states[i].iconf[i_spe], 2);
        }
        ofile.write((char *) &states[i].H, 8);
        ofile.write((char *) &states[i].Hsq, 8);
        ofile.write((char *) &states[i].count, 4);
    }
    ofile.close();
    
    return 0;
} 

void help_message()
{
    //cout << "syntax: mpirun -n 6 ./te [-c] [-f input]" 
         //<< endl;
}


vector <string> load_key_res(char *file_name)
{
    string line;
    vector <string> key_res;

    ifstream res_file(file_name);
    if (res_file.is_open()) {
        while (!res_file.eof()) {
            getline(res_file, line);
            if (line[0] == '#') continue;
            key_res.push_back(line);
        }
        key_res.pop_back();
        res_file.close();
    }
    //printf("load in %d key residues\n", (int) key_res.size());
    for (int i=0; i<key_res.size(); i++) 
        //cout << key_res[i] << endl;

    return key_res;
}

int only_key_res(vector <mstate> &states, char *file_name)
{
    vector <string> key_res = load_key_res(file_name);
    vector <string> all_res = states[0].head_spe;
    int n_key = key_res.size();
    int n_all = all_res.size();
    states[0].head_spe = key_res;
    states[0].n_spe    = n_key;

    vector <int> key_index;
    
    for (int i=0; i<n_key; i++) {
        for (int j=0; j<n_all; j++) {
            if (all_res[j] == key_res[i]) {
                key_index.push_back(j);
                break;
            }
        }
    }

    for (int i=0; i<states.size(); i++) {
        mstate tmp_state = states[i];
        states[i].iconf.clear();
        for (int j=0; j<n_key; j++) {
            states[i].iconf.push_back(tmp_state.iconf[key_index[j]]);
        }
    }
        
    return 0;
}

void out_id(map < unsigned short, pair < string, pair <int, int> > > dic)
{
    //cout << "size of the dic is:  " << dic.size() << endl;
    for (int i=0; i<dic.size(); i++) {
        //printf("id: %d, name: %s, ne: %d, nh: %d\n", i, dic[i].first.c_str(), dic[i].second.first, dic[i].second.second);
    }
}
