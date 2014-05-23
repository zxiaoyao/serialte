#include "mstate.h"

map<unsigned short, pair<string, pair<int, int> > > mstate::id2conf;
int mstate::n_spe = 0;
vector <string> mstate::head_spe;


int mstate::load_id2conf()
{
    string line;

    ifstream ifile("head3.lst", ios::in);
    getline(ifile, line);
    while (getline(ifile, line))
    {
        int iConf = atoi(line.substr(0, 5).c_str());
        int ne    = atoi(line.substr(46, 3).c_str());
        int nH    = atoi(line.substr(49, 3).c_str());

        string confName = line.substr(6, 14);
        id2conf[iConf-1] = pair<string, pair<int, int> > (confName, pair<int, int>(ne, nH));
    }

    return 0;
}

bool mstate::operator==(mstate &source)
{
	bool same = true;
	for (size_t i=0; i<iconf.size(); i++) {
		if (iconf[i] != source.iconf[i]) {
			same = false;
			break;
		}
	}

	return same;
}

bool mstate::operator <(const mstate &rhs) const {
	return iconf < rhs.iconf;
}

bool mstate::same_crg(mstate &source)
{
    bool same = true;
    for (size_t i=0; i<iconf.size(); i++)
    {
        if (id2conf[iconf[i]].second != id2conf[source.iconf[i]].second)
        {
            same = false;
            break;
        }
    }

    return same;
}

mstate& mstate::merge(mstate &source)
{
    count += source.count;
    H     += source.H;
    Hsq   += source.Hsq;

    return *this;
}

mstate& mstate::stat()
{
    e = p = crg = 0;
    H_avg = H_std = 0.0;
    for (int i=0; i<n_spe; i++) {
        e += id2conf[iconf[i]].second.first;
        p += id2conf[iconf[i]].second.second;
    }
    crg = p - e;

    H_avg = H/(double) count;
    if (Hsq - H * H / (double) count < 0.0) {
        H_std = 0.0;
    }
    else {
		H_std = sqrt((Hsq - H * H / (double) count)/(double) count);
    }

    return *this;
}

void mstate::load_binary(ifstream &ifile)
{
    iconf.resize(n_spe);;
    for (int i=0; i<n_spe; i++) {
        ifile.read((char *) &iconf[i], sizeof(unsigned short));
    }
    ifile.read((char *) &H, sizeof(double));
    ifile.read((char *) &Hsq, sizeof(double));
    ifile.read((char *) &count, sizeof(int));
}

void mstate::write_binary(ofstream &ofile)
{
    for (int i=0; i<n_spe; i++) {
        ofile.write((char *) &iconf[i], sizeof(unsigned short));
    }
    ofile.write((char *) &H, sizeof(double));
    ofile.write((char *) &Hsq, sizeof(double));
    ofile.write((char *) &count, sizeof(int));
}
