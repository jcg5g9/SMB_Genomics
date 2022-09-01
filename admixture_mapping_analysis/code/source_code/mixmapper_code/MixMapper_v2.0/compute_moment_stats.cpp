/*
  NOTE:
  - in the bootstrap version, the FIRST rep is the one using all the data (not bootstrapped)
  - in the jackknife version, the LAST rep is the one using all the data (no chroms left out)

  We do not recommend jackknifing for MixMapper analysis.  Bootstrapping is more suitable...
  - jackknifing makes the branch choices appear overly stable: each rep uses 21/22 of the data
    then, after multiplying in the jackknife std err, confidence intervals exceed branch lengths
    => hard to interpret
  - also, not clear how to combine jackknifing with bootstrapping over indivs
 */

#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include <set>
#include <map>
#include <numeric>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cmath>

#include "GenoReader.cpp"

using namespace std;

double sq(double x) { return x*x; }

double std_dev(const vector <double> &x) {
  int n = x.size();
  double xbar = accumulate(x.begin(), x.end(), 0.0) / n;
  double sum = 0;
  for (int i = 0; i < n; i++)
    sum += sq(x[i]-xbar);
  return sqrt(sum / (n-1));
}

int line_count(const char *filename) {
  FILE *fptr = fopen(filename, "r");
  if (fptr == NULL) {
    cerr << "Error opening file " << string(filename) << endl;
    exit(1);
  }
  char line[1024];
  int ans = 0;
  while (fgets(line, 1024, fptr))
    ans++;
  fclose(fptr);
  return ans;
}

int main(int argc, char *argv[]) {

  if (!(argc == 8)) {
    cerr << "usage:" << endl;
    cerr << "- arg1 = .ind file" << endl;
    cerr << "- arg2 = .snp file" << endl;
    cerr << "- arg3 = .geno file" << endl;
    cerr << "- arg4 = output prefix (-> .pops.txt, .f2_boots.txt, .h_boots.txt, .f2.tab," << endl;
    cerr << "                           .neg_f3.txt)" << endl;
    cerr << "- arg5 = resample individuals in each population? (y/n)" << endl;
    cerr << "- arg6 = number of bootstrap replicates" << endl; /*(or j to jackknife over chroms)*/
    cerr << "- arg7 = number of snp blocks (e.g., 50)" << endl; /*(or j to jackknife over chroms)*/
    /*
    cerr << "Note: We do not recommend jackknifing for MixMapper analysis." << endl;
    cerr << "      Bootstrapping is more suitable; see the README for details." << endl;
    */
    exit(1);
  }

  const char *ind_filename = argv[1];
  const char *snp_filename = argv[2];
  const char *geno_filename = argv[3];
  string output_prefix(argv[4]);
  char resample_indivs; sscanf(argv[5], "%c", &resample_indivs);
  if (!(resample_indivs == 'y' || resample_indivs == 'n')) {
    cerr << "need y/n for arg4 (resample individuals in each population?); got "
	 << resample_indivs << endl;
    exit(1);
  }
  bool jackknife_chroms = false;
  int num_resamples; sscanf(argv[6], "%d", &num_resamples);
  int num_blocks; sscanf(argv[7], "%d", &num_blocks);
  if (argv[6][0] == 'j' || argv[7][0] == 'j') jackknife_chroms = true;
  cerr << "params set:" << endl;
  cerr << "- resample individuals in each population? " << resample_indivs << endl;
  string boots_or_jacks;
  if (jackknife_chroms) {
    cerr << "- jackknifing; one rep per chrom followed by one with none removed" << endl;
    boots_or_jacks = "jacks";
  }
  else {
    cerr << "- number of bootstrap resamples: " << num_resamples << endl;
    cerr << "- number of snp blocks: " << num_blocks << endl;
    boots_or_jacks = "boots";
  }

  // generate list of populations using ind file
  map <string, int> pop_name_ind;
  vector <string> pop_names;
  vector <int> pop_inds; // index of population that each individual belongs to; -1 if Ignore
  map <string, vector <int> > pop_samples_map;
  vector < vector <int> > pop_samples;
  ifstream ind_in(ind_filename);
  string indiv_id, gender, pop_name;
  while (ind_in >> indiv_id >> gender >> pop_name) {
    if (pop_name.find("Ignore") != string::npos ||
	pop_name.find("ignore") != string::npos) {
      cerr << "ignoring line: " << indiv_id << " " << gender << " " << pop_name << endl;
      pop_inds.push_back(-1);
    }
    else {
      pop_samples_map[pop_name].push_back(pop_inds.size());
      if (pop_name_ind.find(pop_name) == pop_name_ind.end()) {
	pop_name_ind[pop_name] = pop_names.size();
	pop_names.push_back(pop_name);
      }
      pop_inds.push_back(pop_name_ind[pop_name]);
    }
  }
  ind_in.close();
  int numindivs = pop_inds.size();
  cerr << "number of samples (including Ignore): " << numindivs << endl;
  ofstream fout_pops((output_prefix + ".pops.txt").c_str());
  for (int p = 0; p < (int) pop_names.size(); p++) {
    pop_samples.push_back(pop_samples_map[pop_names[p]]);
    cerr << pop_names[p] << "\t" << pop_samples[p].size() << endl;
    fout_pops << pop_names[p] << endl;
  }
  fout_pops.close();
  int num_pops = pop_names.size();
  cerr << "number of pops: " << num_pops << endl;

  int num_snps = line_count(snp_filename);
  cerr << "number of snps: " << num_snps << endl;

  vector <int> snp_chroms(num_snps);
  vector <bool> snp_ignore(num_snps);
  
  // read snp file to determine chrom of each snp
  set <int> chroms;
  FILE *snp_file = fopen(snp_filename, "r");
  const int buf_size = 1024;
  char snpline[buf_size], ID[buf_size]; int chrom;
  int tot_ignore = 0, tot_non_ignore = 0;
  for (int s = 0; s < num_snps; s++) {
    if (fgets(snpline, buf_size, snp_file) == NULL) {
      cerr << "error: snp file ended at line " << s << "; expected " << num_snps << endl;
      exit(1);
    }
    sscanf(snpline, "%s%d", ID, &chrom); // ignore genpos, physpos in 3rd and 4th cols
    if (!(1 <= chrom && chrom <= 22)) {
      snp_ignore[s] = true;
      tot_ignore++;
    }
    else {
      chroms.insert(chrom);
      tot_non_ignore++;
    }
    snp_chroms[s] = chrom;
  }
  cout << "number of ignored snps: " << tot_ignore << endl;
  cout << "number of non-ignored snps: " << tot_non_ignore << endl;
  vector <int> uniq_chrom_ids(chroms.begin(), chroms.end());
  fclose(snp_file);

  if (jackknife_chroms) num_resamples = chroms.size() + 1;

  vector < vector <double> > h_boots;
  vector < vector < vector <double> > > f2_boots;
  char line[numindivs+10];

  for (int r = 0; r < num_resamples; r++) {

    cerr << "rep " << r;

    int tot_snps = 0;
    vector <double> h_r(num_pops);
    vector < vector <double> > f2_r(num_pops, vector <double> (num_pops));
    
    /***** resample snp blocks *****/
    vector <int> bootstrap_snp_mults(num_snps, 1);
    if (jackknife_chroms) {
      if (r < num_resamples-1) { // not the last, which uses all chromosomes
	for (int s = 0; s < num_snps; s++)
	  bootstrap_snp_mults[s] = snp_chroms[s] != uniq_chrom_ids[r];
      }
    }
    else {
      if (r) { // not the first, which doesn't resample
	vector <int> bootstrap_block_mults(num_blocks);
	for (int k = 0; k < num_blocks; k++)
	  bootstrap_block_mults[rand()%num_blocks]++;
	int num_per_block = (num_snps+num_blocks-1)/num_blocks;
	for (int s = 0; s < num_snps; s++)
	  bootstrap_snp_mults[s] = bootstrap_block_mults[s/num_per_block];
      }
    }

    /***** resample each population *****/
    vector <int> bootstrap_sample_mults(numindivs, 1);
    if (resample_indivs == 'y') {
      if (r)  // not the first, which doesn't resample
	for (int x = 0; x < num_pops; x++) {
	  int pop_size = pop_samples[x].size();
	  for (int k = 0; k < pop_size; k++) {
	    bootstrap_sample_mults[pop_samples[x][k]]--;
	    bootstrap_sample_mults[pop_samples[x][rand()%pop_size]]++;
	  }
	}
    }

    /***** read geno file and compute f2 statistics *****/
    try {
      IOUtils::GenoReader genoReader(geno_filename, numindivs, num_snps);

      for (int s = 0; s < num_snps; s++) {
	if ((s & 0x3fff) == 0)
	  cerr << "." << flush;
	genoReader.read_line(line);
	if (snp_ignore[s] || bootstrap_snp_mults[s] == 0) continue;

	/***** compute multipliers for bias term *****/
	// ... and compute population allele totals a[pop_ind], b[pop_ind]
	vector <double> a(num_pops), b(num_pops), bias_mult(num_pops);
	vector <int> n(num_pops), num_sq(num_pops);
	for (int i = 0; i < numindivs; i++) {
	  int x = pop_inds[i];
	  if (x == -1) continue;
	  int gtype = line[i]-'0';	
	  if (gtype != 9) {
	    int sample_mult = bootstrap_sample_mults[i];
	    a[x] += sample_mult * gtype;
	    b[x] += sample_mult * (2-gtype);
	    n[x] += 2 * sample_mult; // diploid: 2 * ...
	    num_sq[x] += 2 * sample_mult * sample_mult;
	  }
	}
	for (int x = 0; x < num_pops; x++) {
	  bias_mult[x] = num_sq[x] / sq(n[x]);
	  bias_mult[x] /= 1-bias_mult[x];
	}
      
	// bad_snp if any pop has <= 1 good sample
	bool bad_snp = false;
	for (int x = 0; x < num_pops; x++)
	  if (n[x] == 0 || n[x]*n[x] == num_sq[x])
	    bad_snp = true;
	if (bad_snp) continue;
      
	int snp_mult = bootstrap_snp_mults[s];
	tot_snps += snp_mult;
      
	// augment f2 and h totals
	for (int x = 0; x < num_pops; x++) {

	  double p_x = a[x]/(double) (a[x]+b[x]);
	  double bias_x = bias_mult[x] * a[x]*b[x] / sq(a[x]+b[x]);

	  h_r[x] += 2 * (p_x*(1-p_x) + bias_x) * snp_mult;

	  for (int y = x+1; y < num_pops; y++) {
	
	    double p_y = a[y]/(double) (a[y]+b[y]);
	    double bias_y = bias_mult[y] * a[y]*b[y] / sq(a[y]+b[y]);
	
	    f2_r[x][y] += ((p_x-p_y)*(p_x-p_y) - bias_x - bias_y) * snp_mult;
	  }
	}
      }
      if (!genoReader.check_eof()) {
	printf("ERROR: expected EOF after %d snps, but file still has data\n", num_snps);
	exit(1);
      }
      else
	cerr << " done (" << num_snps << ")" << endl;
      genoReader.close();
    }
    catch (string errstr) {
      cout << endl << "ERROR: " << errstr << endl; exit(1);
    }
    cout << "number of snps actually used in bootstrap rep (inc. resamples): " << tot_snps << endl;
    // create data matrices for this rep
    for (int x = 0; x < num_pops; x++) {
      h_r[x] /= tot_snps;
      for (int y = x+1; y < num_pops; y++)
	f2_r[x][y] = f2_r[y][x] = f2_r[x][y] / tot_snps;
    }
    h_boots.push_back(h_r);
    f2_boots.push_back(f2_r);
  }


  /***** write output *****/

  // .f2.tab file (tab-delimited f2 matrix with row and column headers)
  ofstream fout_f2_tab((output_prefix + ".f2.tab").c_str());
  int r_all = jackknife_chroms ? num_resamples-1 : 0; // index of all-data rep
  for (int x = 0; x < num_pops; x++)
    fout_f2_tab << '\t' << pop_names[x];
  fout_f2_tab << endl;
  for (int x = 0; x < num_pops; x++) {
    fout_f2_tab << pop_names[x];
    for (int y = 0; y < num_pops; y++)
      fout_f2_tab << '\t' << f2_boots[r_all][x][y];
    fout_f2_tab << endl;
  }
  fout_f2_tab.close();

  // .f2_boots or .f2_jacks file (f2 values for all bootstrap/jackknife reps)
  ofstream fout_f2((output_prefix + ".f2_" + boots_or_jacks + ".txt").c_str());
  for (int x = 0; x < num_pops; x++) {
    for (int y = 0; y < num_pops; y++) {
      for (int r = 0; r < num_resamples; r++) {
	if (r) fout_f2 << ' ';
	fout_f2 << f2_boots[r][x][y];
      }
      fout_f2 << endl;
    }
  }
  fout_f2.close();

  // .h_boots or .h_jacks file (heterozygosity values for all bootstrap/jackknife reps)
  ofstream fout_h((output_prefix + ".h_" + boots_or_jacks + ".txt").c_str());
  for (int x = 0; x < num_pops; x++) {
    for (int r = 0; r < num_resamples; r++) {
      if (r) fout_h << ' ';
      fout_h << h_boots[r][x];
    }
    fout_h << endl;
  }
  fout_h.close();

  // .neg_f3.txt (triples with negative f3 values, indicating admixture)
  FILE *fptr_neg_f3 = fopen((output_prefix + ".neg_f3.txt").c_str(), "w");
  int n = num_resamples-1;
  vector <double> smallest_f3_zscore(num_pops, 1000.0);
  for (int c = 0; c < num_pops; c++)
    for (int a = 0; a < num_pops; a++)
      for (int b = a+1; b < num_pops; b++)
	if (c != a && c != b) {
	  double f3_all_data, f3_std;
	  if (jackknife_chroms) {
	    // for jackknife, last one is with none left out
	    f3_all_data = (f2_boots[n][a][c] + f2_boots[n][b][c] - f2_boots[n][a][b]) / 2;
	    vector <double> f3(n);
	    for (int i = 0; i < n; i++)
	      f3[i] = (f2_boots[i][a][c] + f2_boots[i][b][c] - f2_boots[i][a][b]) / 2;
	    f3_std = std_dev(f3) * (n-1) / sqrt(n);
	  }
	  else {
	    // for bootstrap, first one is with none left out
	    f3_all_data = (f2_boots[0][a][c] + f2_boots[0][b][c] - f2_boots[0][a][b]) / 2;
	    vector <double> f3(n);
	    for (int i = 1; i < num_resamples; i++)
	      f3[i-1] = (f2_boots[i][a][c] + f2_boots[i][b][c] - f2_boots[i][a][b]) / 2;
	    f3_std = std_dev(f3);
	  }
	  smallest_f3_zscore[c] = min(smallest_f3_zscore[c], f3_all_data / f3_std);
	  if (f3_all_data < 0)
	    fprintf(fptr_neg_f3, "%-20s%-20s%-20s%10.5f%10.5f%10.2f\n",
		    pop_names[c].c_str(), pop_names[a].c_str(),
		    pop_names[b].c_str(), f3_all_data, f3_std, f3_all_data / f3_std);
	}
  fclose(fptr_neg_f3);

  cerr << "smallest f3 z-score for each population (negative indicates admixture):" << endl;
  for (int c = 0; c < num_pops; c++)
    fprintf(stderr, "%-20s%10.2f\n", pop_names[c].c_str(), smallest_f3_zscore[c]);

  cerr << "printing filtered pop list (eliminating pops with an f3 z-score < -2) to stdout"
       << endl;
  for (int c = 0; c < num_pops; c++)
    if (smallest_f3_zscore[c] >= -2)
      cout << pop_names[c] << endl;
}
