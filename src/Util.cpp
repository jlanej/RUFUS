#include <unistd.h>
#include <ios>
#include <bitset>
#include <string>
#include <vector>
#include <iostream>
#include <fstream>
#include <sstream>
#include <math.h>
#include <sys/resource.h>
#include <stdio.h>
#include <stdlib.h>


#include "Util.h"

using namespace std;

bool Util::fncomp(char lhs, char rhs) 
{
  return lhs<rhs;
}

 const vector<string> Util::Split(const string& line, const char delim)
{
  vector<string> tokens;
  stringstream lineStream(line);
  string token;
  while(getline(lineStream, token, delim)){
    tokens.push_back(token);
  }
  return tokens;
}

string Util::trim(string s)
{
  string newS = "";
  cout << s << endl;
  for (int i = 0; i < s.size(); i++)
    {
      if ((int)s.c_str()[i]>=33)
	{
	  newS = newS+s.c_str()[i];
	  cout << s.c_str()[i] << endl;
	}
    }
  cout << newS << endl;
  return newS;
}

unsigned long Util::HashToLong(const string& hash)
{
  unsigned long result = 0;
  for(int i=0; i<hash.length();i++)
    {
      char c = hash[i];
      unsigned long bits;
      switch(c) {
        case 'A': bits = 0; break;
        case 'C': bits = 2; break;
        case 'G': bits = 1; break;
        case 'T': bits = 3; break;
        default:
          cout << "ERROR, invalid character - " << c << endl;
          bits = 0;
          break;
      }
      result |= (bits << (i * 2));
    }
  return result;
}

unsigned long HashToLongTB(const string& hash, const string& calledby)
{
  unsigned long result = 0;
  int bitcout = 0;
  for(int i=0; i<hash.size();i++)
    {
      char c = hash[i];
      unsigned long bits;
      switch(c) {
        case 'A': bits = 0; break;
        case 'C': bits = 2; break;
        case 'G': bits = 1; break;
        case 'T': bits = 3; break;
        default:
          cout << "ERROR, invalid character - " << c << ", in hash " << hash << ", called by " << calledby <<  endl;
          bits = 0;
          break;
      }
      result |= (bits << bitcout);
      bitcout += 2;
    }
  return result;
}

string Util::LongToHash(unsigned long LongHash, int HashSize)
{
  string value;
  value.reserve(HashSize);
  for (int i = 0; i < HashSize; i++)
    {
      unsigned long bits = (LongHash >> (i * 2)) & 3;
      switch(bits) {
        case 0: value += 'A'; break;
        case 2: value += 'C'; break;
        case 1: value += 'G'; break;
        case 3: value += 'T'; break;
      }
    }
  return value;
}

/*string Util::RevComp(string Sequence)
{
  string NewString = "";
  for(int i = Sequence.length(); i > 0; --i)
    {
      switch(char(Sequence.c_str()[i]))
	{
	case 'A':
	  NewString += 'T';
	case 'C':
	  NewString += 'G';
	case 'G':
	  NewString += 'C';
	case 'T':
	  NewString += 'A';
	case 'N':
	  NewString += 'N';
	default: 
	  cout << "ERROR IN RevComp - \n" << Sequence.c_str()[i] << endl;
	}
    }
  return NewString;
  }*/

  string Util::RevComp (const string& Sequence)
{
  int len = Sequence.size();
  string NewString;
  NewString.reserve(len);
  for(int i = len-1; i>=0; i--)
    {
      char C = Sequence[i];
      switch(C) {
        case 'A': NewString += 'T'; break;
        case 'C': NewString += 'G'; break;
        case 'G': NewString += 'C'; break;
        case 'T': NewString += 'A'; break;
        case 'N': NewString += 'N'; break;
        default: cout << "\nERROR IN RevComp - " << C << "\n"; break;
      }
    }
  return NewString;
}

string Util::RevQual(const string& Sequence)
{
  int len = Sequence.size();
  string NewString;
  NewString.reserve(len);
  for(int i = len-1; i>=0; i--)
    {
      unsigned char C = Sequence[i];
      if (C != '\0')
	{NewString += C;}
    }
  return NewString;
}

void Util::process_mem_usage(double& vm_usage, double& resident_set, double& MAXvm, double& MAXrss)
{
  using std::ios_base;
  using std::ifstream;
  using std::string;

  vm_usage     = 0.0;
  resident_set = 0.0;

  // 'file' stat seems to give the most reliable results                                                                                                               

  ifstream stat_stream("/proc/self/stat",ios_base::in);

  // dummy vars for leading entries in stat that we don't care about                                                                                                   
  //                                                                                                                                                                   
  string pid, comm, state, ppid, pgrp, session, tty_nr;
  string tpgid, flags, minflt, cminflt, majflt, cmajflt;
  string utime, stime, cutime, cstime, priority, nice;
  string O, itrealvalue, starttime;

  // the two fields we want                                                                                                                                            

  unsigned long vsize;
  long rss;

  stat_stream >> pid >> comm >> state >> ppid >> pgrp >> session >> tty_nr
              >> tpgid >> flags >> minflt >> cminflt >> majflt >> cmajflt
              >> utime >> stime >> cutime >> cstime >> priority >> nice
              >> O >> itrealvalue >> starttime >> vsize >> rss; // don't care about the rest                                                                           

  stat_stream.close();

  long page_size_kb = sysconf(_SC_PAGE_SIZE) / 1024; // in case x86-64 is configured to use 2MB pages                                                                  
  vm_usage     = vsize / 1024.0;
  resident_set = rss * page_size_kb;
  if (vm_usage > MAXvm){MAXvm = vm_usage;}
  if (resident_set > MAXrss){MAXrss = resident_set;}
}

