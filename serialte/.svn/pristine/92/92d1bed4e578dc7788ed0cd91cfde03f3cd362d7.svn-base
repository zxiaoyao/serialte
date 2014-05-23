#ifndef _MSTATE_H_
#define _MSTATE_H_

#include <iostream>
#include <fstream>

#include <string>
#include <vector>
#include <map>

#include <cstdlib>
#include <cmath>

using namespace std;
class mstate
{
    public:
        mstate() : H(0.0), Hsq(0.0), count(0), e(0), p(0), crg(0), H_avg(0.0), H_std(0.0) {};
        mstate(const mstate &rhs): H(rhs.H), Hsq(rhs.Hsq), count(rhs.count), e(rhs.e),
        		p(rhs.p), crg(rhs.crg), H_avg(rhs.H_avg), H_std(rhs.H_std), iconf(rhs.iconf) {};
	
	    // check whether have the same microstate with source state
        bool operator== (mstate &source);

        bool operator<(const mstate &rhs) const;
		
		// check whether have the same charge state with source microstate
        bool same_crg(mstate &source);
	
		// merge another same microstate "source"
        mstate& merge(mstate &source);
	
		// some statistics of this microstate
        mstate& stat();
		
		// load one microstate from a binary file
        void load_binary(ifstream &ifile);
		// write one microstate to a binary file
        void write_binary(ofstream &ofile);
	
	
	    // read from head3.lst
        static int load_id2conf(void);

	
	/* data members of a microstate:
	    H:			enthalpy of the state
		Hsq:		square of the enthalpy
	    count:		times of appearence of this state
	    hashvalue:	hash value of a microstate, used to compare two states
		e:			number of electron
		p:			number of proton
		crg:		charge of the state
		H_avg:		average enthalpy of this same microstate (probably we don't select all the residues to be the special ones)
		H_std:		std of a series of this same microstate
		iconf:		the conformer sequence of a microstate
	 */

        double H;
        double Hsq;
        int count;
//        int hashvalue;

        int e;
        int p;
        int crg;
        double H_avg;
        double H_std;
        vector <unsigned short> iconf;
	/* static members of the microstate class
		head_spe:		a sequence of special residue names for all the microstates
		n_spe:			number of special residue names
		id2conf:		a map between ID and conformer name  (id, (name, (e, p)))
	 */
        static vector <string> head_spe;
        static int n_spe;
        static map < unsigned short, pair < string, pair <int, int> > > id2conf;
};

#endif
