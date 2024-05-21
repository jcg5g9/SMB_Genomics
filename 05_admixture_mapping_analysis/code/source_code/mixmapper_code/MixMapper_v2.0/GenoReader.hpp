#ifndef GENOREADER_HPP
#define GENOREADER_HPP

#include <cstdio>

namespace IOUtils {

  class GenoReader {

    static const int MIN_RLEN = 48; // must be consistent with mcio.c!

    FILE *fin;
    int numindivs;
    bool is_packed;
    int rlen; // number of bytes to read per SNP, if packed format
    
    static const int BUF_LEN = 1<<20; // 1 MB
    char buf[BUF_LEN];

  public:

    GenoReader(const char *filename, int _numindivs, int numsnps);
    void read_line(char *buf);
    bool check_eof(void);
    void close(void);
    
  };

}

#endif
