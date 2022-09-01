#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include <set>
#include <map>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cmath>
#include <algorithm>

using namespace std;

typedef unsigned long long mask_t;

int popcount(mask_t x)
{
  int c = 0;
  for (; x > 0; x &= x -1) c++;
  return c;
}

vector < vector <double> > make_f2_submat(const vector < vector <double> > &f2,
					  const vector <int> &unmixed_inds, mask_t mask) {
  int num_unmixed = unmixed_inds.size();
  int n = popcount(mask);
  vector < vector <double> > d(n, vector <double> (n));
  int isquash = 0;
  for (int i = 0; i < num_unmixed; i++)
    if ((mask>>i)&1) {
      int jsquash = 0;
      for (int j = 0; j < num_unmixed; j++)
	if ((mask>>j)&1) {
	  d[isquash][jsquash] = f2[unmixed_inds[i]][unmixed_inds[j]];
	  jsquash++;
	}
      isquash++;
    }
  return d;
}

double compute_tree_err(const vector < vector <double> > &d_in) {
  int num_unmixed = d_in.size();
  if (num_unmixed < 3) return 0;
  vector < vector <double> > d(2*num_unmixed, vector <double> (2*num_unmixed));
  for (int i = 0; i < num_unmixed; i++)
    for (int j = 0; j < num_unmixed; j++)
      d[i][j] = d_in[i][j];
  vector < vector <double> > dtree(2*num_unmixed, vector <double> (2*num_unmixed, 1e9));
  vector <bool> used(num_unmixed);

  int num_nodes = num_unmixed; // will increase

  double r[num_unmixed*2];
  int nodes_left = num_unmixed;
  // main neighbor joining loop
  for (; nodes_left > 2; nodes_left--) {
    for (int i = 0; i < num_nodes; i++)
      if (!used[i]) {
	r[i] = 0;
	for (int k = 0; k < num_nodes; k++)
	  if (!used[k])
	    r[i] += d[i][k];
	r[i] /= nodes_left-2;
      }
    double min_D = 1e9;
    int best_i = 0, best_j = 0;
    for (int i = 0; i < num_nodes; i++)
      if (!used[i])
	for (int j = i+1; j < num_nodes; j++)
	  if (!used[j]) {
	    double D = d[i][j] - r[i] - r[j];
	    if (D < min_D) {
	      min_D = D;
	      best_i = i; best_j = j;
	    }
	  }
    int k = num_nodes;
    for (int m = 0; m < num_nodes; m++)
      d[k][m] = d[m][k] = 0.5 * (d[best_i][m]+d[best_j][m]-d[best_i][best_j]);
    d[best_i][k] = d[k][best_i] = 0.5*(d[best_i][best_j]+r[best_i]-r[best_j]);
    d[best_j][k] = d[k][best_j] = 0.5*(d[best_i][best_j]+r[best_j]-r[best_i]);
    dtree[best_i][k] = dtree[k][best_i] = d[k][best_i];
    dtree[best_j][k] = dtree[k][best_j] = d[k][best_j];
    used[best_i] = used[best_j] = true;
    used.push_back(false);
    num_nodes++;
  }
  // last neighbor joining step
  vector <int> last2;
  for (int i = 0; i < num_nodes; i++)
    if (!used[i])
      last2.push_back(i);
  dtree[last2[0]][last2[1]] = dtree[last2[1]][last2[0]] = d[last2[0]][last2[1]];

  // compute distances along tree using Floyd-Warshall
  for (int k = 0; k < num_nodes; k++)
    for (int i = 0; i < num_nodes; i++)
      if (dtree[i][k] < 1e9)
	for (int j = 0; j < num_nodes; j++)
	  dtree[i][j] = min(dtree[i][j], dtree[i][k] + dtree[k][j]);

  // compute and return errors
  double max_err = 0;
  for (int i = 0; i < num_unmixed; i++)
    for (int j = i+1; j < num_unmixed; j++) {
      double err = dtree[i][j]-d[i][j];
      max_err = max(max_err, abs(err));
    }
  return max_err;
}

void output_bests(const set < pair <double, mask_t> > &bests,
		  const vector <string> &unmixed_pops) {
  const int max_to_output = 25;
  int num_unmixed = unmixed_pops.size();
  int iter = 0;
  for (set < pair <double, mask_t> >::iterator it = bests.begin();
       it != bests.end() && iter < max_to_output; it++, iter++) {
    double err = it->first; mask_t mask = it->second;
    int n = popcount(mask);
    cout << "size = " << n << ", max err = " << err << endl;
    cout << "present:";
    for (int i = 0; i < num_unmixed; i++)
      if ((mask>>i)&1)
	cout << " " << unmixed_pops[i];
    cout << endl;
    cout << "absent:";
    for (int i = 0; i < num_unmixed; i++)
      if (!((mask>>i)&1))
	cout << " " << unmixed_pops[i];
    cout << endl;
    cout << endl;
  }
  cout << "----------------------------------------------------------------------" << endl;
  cout << endl;
}

int main(int argc, char *argv[]) {
  
  if (!(argc == 4 || argc == 5)) {
    cerr << "usage:" << endl;
    cerr << "- arg1 = .f2.tab file (from output of compute_moment_stats)" << endl;
    cerr << "- arg2 = max number of subsets of each size to analyze (e.g., 10000;" << endl;
    cerr << "             increase if you wish to perform a more exhaustive search)" << endl;
    cerr << "- arg3 = file containing list of populations to choose from" << endl;
    cerr << "- (optional) arg4 = file containing populations required to be in the tree" << endl;
    exit(1);
  }

  // input f2 tab matrix
  ifstream f2_in(argv[1]);
  string str;
  getline(f2_in, str); // header
  vector <string> pops;
  istringstream iss(str);
  while (iss >> str) pops.push_back(str);
  int num_pops = pops.size(); cerr << "num pops: " << num_pops << endl;
  vector < vector <double> > f2(num_pops, vector <double> (num_pops));
  for (int i = 0; i < num_pops; i++) {
    f2_in >> str;
    for (int j = 0; j < num_pops; j++)
      f2_in >> f2[i][j];
  }
  f2_in.close();

  int max_subsets; sscanf(argv[2], "%d", &max_subsets);

  // input unmixed pops
  ifstream unmixed_in(argv[3]);
  vector <int> unmixed_inds;
  vector <string> unmixed_pops;
  while (unmixed_in >> str) {
    unmixed_pops.push_back(str);
    bool found = false;
    for (int i = 0; i < num_pops; i++)
      if (pops[i] == str) {
	found = true;
	unmixed_inds.push_back(i);
	break;
      }
    if (!found) {
      cerr << "error: unmixed pop not found: " << str << endl;
      exit(1);
    }
  }
  int num_unmixed = unmixed_inds.size();
  if (num_unmixed > 64) {
    cerr << "error: maximum of 64 unmixed pops allowed" << endl;
    exit(1);
  }

  // input required pops
  vector <int> required_inds_in_unmixed;
  int min_leaves; set < pair <double, mask_t> > cur_bests;
  if (argc == 5) {
    ifstream required_in(argv[4]);
    mask_t mask = 0;
    while (required_in >> str) {
      bool found = false;
      for (int i = 0; i < num_unmixed; i++)
	if (pops[unmixed_inds[i]] == str) {
	  found = true;
	  required_inds_in_unmixed.push_back(i);
	  cerr << "requiring pop" << i << ": " << pops[unmixed_inds[i]] << endl;
	  mask |= ((mask_t) 1)<<i;
	  break;
	}
      if (!found) {
	cerr << "error: required pop not found: " << str << endl;
	exit(1);
      }
    }
    min_leaves = required_inds_in_unmixed.size();
    cur_bests.insert(make_pair(compute_tree_err(make_f2_submat(f2, unmixed_inds, mask)), mask));
    if (min_leaves >= 4) output_bests(cur_bests, unmixed_pops);
  }
  else {
    min_leaves = 0;
    cur_bests.insert(make_pair(0.0, 0));
  }

  // iterate through subsets of expanding sizes
  for (int leaves = min_leaves+1; leaves <= num_unmixed; leaves++) {
    set < pair <double, mask_t> > next_bests;
    set <mask_t> next_masks_seen; // record subsets seen to avoid repeating work
    
    int iter = 0;
    for (set < pair <double, mask_t> >::iterator it = cur_bests.begin(); it != cur_bests.end();
	 it++) {
      mask_t cur_mask = it->second;
      for (int i = 0; i < num_unmixed; i++)
	if (!((cur_mask>>i)&1)) {
	  mask_t mask = cur_mask | (((mask_t) 1)<<i);
	  if (next_masks_seen.find(mask) == next_masks_seen.end()) {
	    next_masks_seen.insert(mask);
	    next_bests.insert(make_pair(compute_tree_err(make_f2_submat(f2, unmixed_inds, mask)),
					mask));
	  }
	}
      iter++;
      if (iter == max_subsets)
	break;
    }

    if (leaves >= 4) output_bests(next_bests, unmixed_pops);
    cur_bests = next_bests;
  }
}
