#include "GenoReader.hpp"

#include <string>
#include <cstring>

namespace IOUtils {

  using std::string;
  using std::max;

/*
  constructor:
  - open file
  - check that formatting is correct (either 0129 matching #indivs or GENO %d %d %x %x)
  - set is_packed flag for packed or not
  - set rlen
  - move file pointer to start of actual data

  method: read line into preallocated buffer (convert to 0129)

  destructor: close file
 */

  GenoReader::GenoReader(const char *filename, int _numindivs, int numsnps) {
    numindivs = _numindivs;
    fin = fopen(filename, "rb");
    if (fin == NULL) throw string("unable to open geno file\n");

    // check if it's packed; set is_packed accordingly
    rlen = max((numindivs+3)/4, MIN_RLEN); // must be consistent with mcio.c!
    fread(buf, 1, rlen, fin);
    is_packed = false;
    if (strncmp(buf, "GENO", 4) == 0) {
      int xnind, xnsnp, xihash, xshash;
      sscanf(buf,"GENO %d %d %x %x", &xnind, &xnsnp, &xihash, &xshash);
      if (xnind != numindivs)
	throw string("packed geno file has wrong number of indivs according to header");
      if (xnsnp != numsnps)
	throw string("packed geno file has wrong number of snps according to header");
      is_packed = true;
    }

    if (!is_packed) {
      fclose(fin);
      fin = fopen(filename, "r"); // open as text file
    }
  }
  
  void GenoReader::read_line(char *line) {
    const char lookup[] = "0129";
    if (is_packed) { // read the packed line and convert to 0129
      int num_read = fread(buf, 1, rlen, fin);
      if (num_read != rlen)
	throw string("error reading packed geno file (probably premature EOF)");
      for (int num = 0; num < numindivs; num++) {
	int wnum = num >> 2 ;
	int wplace = num & 3 ; 

	int rshft = (3-wplace) << 1 ;
	int b = buf[wnum] >> rshft ;

	b = b & 3 ;
	line[num] = lookup[b];
      }
    }
    else {
      if (fgets(buf, BUF_LEN, fin) == NULL)
	throw string("premature EOF (fewer lines of data than SNPs)");
      if ((int) strlen(buf) != numindivs+1) {
	char errbuf[1000];
	sprintf(errbuf, "geno file line has wrong length:\n"
		"- expected %d (one digit per indiv)\n"
		"- got %d\n"
		"check data files...\n"
		"- does .geno file match .ind and .snp files?\n"
		"- is .geno file format is actually eigenstrat?\n"
		"http://genepath.med.harvard.edu/~reich/InputFileFormats.htm",
		numindivs, (int) strlen(buf)-1);
	throw string(errbuf);
      }
      memcpy(line, buf, numindivs);
    }
  }
  
  bool GenoReader::check_eof(void) {
    int c = fgetc(fin);
    if (c != EOF) {
      ungetc(c, fin);
      return false;
    }
    return true;
  }

  void GenoReader::close(void) {
    fclose(fin);
  }

}
