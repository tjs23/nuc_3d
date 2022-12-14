cdef extern from "string.h":
  ctypedef int size_t
  void *memcpy(void *dst,void *src,size_t len)
  void *memmove(void *dst,void *src,size_t len)
  void *memset(void *b,int c,size_t len)

cdef extern from "stdlib.h":
  void free(void *)
  void *malloc(size_t)
  void *calloc(size_t,size_t)
  void *realloc(void *,size_t)
  int c_abs "abs" (int)
  void qsort(void *base, size_t nmemb, size_t size,
             int (*compar)(void *,void *))

cdef extern from "math.h":
   double sqrt(double x)

cdef extern from "stdio.h":
  ctypedef struct FILE:
    pass
  FILE *fopen(char *,char *)
  FILE *freopen(char *path, char *mode, FILE *stream)
  int fileno(FILE *stream)
  int dup2(int oldfd, int newfd)
  int fflush(FILE *stream)

  FILE * stderr
  FILE * stdout
  int fclose(FILE *)
  int sscanf(char *str,char *fmt,...)
  int printf(char *fmt,...)
  int sprintf(char *str,char *fmt,...)
  int fprintf(FILE *ifile,char *fmt,...)
  char *fgets(char *str,int size,FILE *ifile)

cdef extern from "ctype.h":
  int toupper(int c)
  int tolower(int c)
  
cdef extern from "unistd.h":
  char *ttyname(int fd)
  int isatty(int fd)  

cdef extern from "string.h":
  int strcmp(char *s1, char *s2)
  int strncmp(char *s1,char *s2,size_t len)
  char *strcpy(char *dest,char *src)
  char *strncpy(char *dest,char *src, size_t len)
  char *strdup(char *)
  char *strcat(char *,char *)
  size_t strlen(char *s)
  int memcmp( void * s1, void *s2, size_t len )

cdef extern from "Python.h":
   long _Py_HashPointer(void*)
   FILE* PyFile_AsFile(object)


cdef extern from "stdint.h":
  ctypedef int int8_t
  ctypedef int int16_t
  ctypedef int int32_t
  ctypedef int int64_t
  ctypedef int uint8_t
  ctypedef int uint16_t
  ctypedef int uint32_t
  ctypedef int uint64_t

cdef extern from "zlib.h":
  ctypedef void * gzFile
  ctypedef int64_t z_off_t

  int gzclose(gzFile fp)
  int gzread(gzFile fp, void *buf, unsigned int n)
  char *gzerror(gzFile fp, int *errnum)

  gzFile gzopen( char *path, char *mode)
  gzFile gzdopen (int fd, char *mode)
  char * gzgets(gzFile file, char *buf, int len)
  int gzeof( gzFile file )

cdef extern from "samtools/bam.h":

  # constants
  int BAM_DEF_MASK
  # IF _IOLIB=2, bamFile = BGZF, see bgzf.h
  # samtools uses KNETFILE, check how this works

  ctypedef struct tamFile:
      pass

  ctypedef struct bamFile:
      pass

  ctypedef struct bam1_core_t:
      int32_t tid 
      int32_t pos
      uint32_t bin
      uint32_t qual
      uint32_t l_qname
      uint32_t flag
      uint32_t n_cigar
      int32_t l_qseq
      int32_t mtid 
      int32_t mpos 
      int32_t isize

  ctypedef struct bam1_t:
    bam1_core_t core
    int l_aux
    int data_len
    int m_data
    uint8_t *data

  ctypedef struct bam_pileup1_t:
      bam1_t *b 
      int32_t qpos 
      int indel
      int level
      uint32_t is_del
      uint32_t is_head
      uint32_t is_tail

  ctypedef int (*bam_pileup_f)(uint32_t tid, uint32_t pos, int n, bam_pileup1_t *pl, void *data)

  ctypedef int (*bam_fetch_f)(bam1_t *b, void *data)

  ctypedef struct bam_header_t:
     int32_t n_targets
     char **target_name
     uint32_t *target_len
     void *hash
     void *rg2lib
     int l_text
     char *text

  ctypedef struct bam_index_t:
      int32_t n
      uint64_t n_no_coor

  ctypedef struct bam_plbuf_t:
      pass

  ctypedef struct pair64_t:
      uint64_t u, v
      
  ctypedef struct bam_iter_t:
      int from_first
      int tid, beg, end, n_off, i, finished
      uint64_t curr_off
      pair64_t *off

  # ctypedef __bam_iter_t * bam_iter_t

  bam1_t * bam_init1()
  void bam_destroy1(bam1_t *)

  bamFile razf_dopen(int data_fd, char *mode)

  int64_t bam_seek( bamFile fp, uint64_t voffset, int where)
  int64_t bam_tell( bamFile fp )

  # void bam_init_header_hash(bam_header_t *header)

  ###############################################
  # stand-ins for samtools macros
  uint32_t * bam1_cigar( bam1_t * b)
  char * bam1_qname( bam1_t * b)
  uint8_t * bam1_seq( bam1_t * b)
  uint8_t * bam1_qual( bam1_t * b)
  uint8_t * bam1_aux( bam1_t * b)

  ###############################################
  # bam iterator interface
  bam_iter_t bam_iter_query( bam_index_t *idx, int tid, int beg, int end)

  int bam_iter_read(bamFile fp, bam_iter_t iter, bam1_t *b)

  void bam_iter_destroy(bam_iter_t iter)

  ###############################################

  bam1_t * bam_dup1( bam1_t *src ) 
  
  bam1_t * bam_copy1(bam1_t *bdst, bam1_t *bsrc)
  bam_index_t *bam_index_load(char *f )

  void bam_index_destroy(bam_index_t *idx)

  int bam_parse_region(bam_header_t *header, char *str, int *ref_id, int *begin, int *end)

  ###############################################
  bam_plbuf_t *bam_plbuf_init(bam_pileup_f func, void *data)

  int bam_fetch(bamFile fp, bam_index_t *idx, int tid, int beg, int end, void *data, bam_fetch_f func)

  int bam_plbuf_push(bam1_t *b, bam_plbuf_t *buf)

  void bam_plbuf_destroy(bam_plbuf_t *buf)
  ########################################
  # pileup iterator interface
  ctypedef struct bam_plp_t:
      pass

  ctypedef bam_pileup1_t * const_bam_pileup1_t_ptr "const bam_pileup1_t *"

  ctypedef int (*bam_plp_auto_f)(void *data, bam1_t *b)

  bam_plp_t bam_plp_init( bam_plp_auto_f func, void *data)
  int bam_plp_push( bam_plp_t iter,  bam1_t *b)
  bam_pileup1_t * bam_plp_next( bam_plp_t iter, int *_tid, int *_pos, int *_n_plp)
  bam_pileup1_t * bam_plp_auto( bam_plp_t iter, int *_tid, int *_pos, int *_n_plp)
  void bam_plp_set_mask(bam_plp_t iter, int mask)
  void bam_plp_reset(bam_plp_t iter)
  void bam_plp_destroy(bam_plp_t iter)
  void bam_plp_set_maxcnt(bam_plp_t iter, int maxcnt)

  ##################################################

  int bam_read1( bamFile fp, bam1_t *b)
  int bam_validate1( bam_header_t *header, bam1_t *b)
  int bam_write1( bamFile fp, bam1_t *b)

  bam_header_t *bam_header_init()

  int bam_header_write( bamFile fp, bam_header_t *header)

  bam_header_t *bam_header_read( bamFile fp )

  void bam_header_destroy(bam_header_t *header)

  bam1_t * bam_dup1( bam1_t *src ) 
  
  bam1_t * bam_copy1(bam1_t *bdst, bam1_t *bsrc)

  # functions for dealing with the auxillary string
  uint8_t *bam_aux_get(bam1_t *b,  char tag[2])

  void bam_aux_append(bam1_t *b, char tag[2], char type, int len, uint8_t *data)

  int bam_aux_del(bam1_t *b, uint8_t *s)

  # type conversion functions
  int32_t bam_aux2i(uint8_t *s)
  float bam_aux2f(uint8_t *s)
  double bam_aux2d(uint8_t *s)
  char bam_aux2A( uint8_t *s)
  char *bam_aux2Z( uint8_t *s)
  int bam_aux_type2size( int )
  
  # determine indexing bin for a region
  int bam_reg2bin(uint32_t beg, uint32_t end)

  # calculate alignment end position from a cigar string
  uint32_t bam_calend(bam1_core_t *c, uint32_t *cigar)

cdef extern from *:
    ctypedef char* const_char_ptr "const char*"

cdef extern from "samtools/sam.h":

  ctypedef struct samfile_t_un:
    tamFile tamr
    bamFile bam
    FILE *tamw
    
  ctypedef struct samfile_t:
     int type
     samfile_t_un x
     bam_header_t *header

  samfile_t *samopen( const_char_ptr fn, char * mode, void *aux)

  int sampileup( samfile_t *fp, int mask, bam_pileup_f func, void *data)

  void samclose(samfile_t *fp)

  int samread(samfile_t *fp, bam1_t *b)

  int samwrite(samfile_t *fp, bam1_t *b)

  # functions not declared in sam.h but available as extern
  int bam_prob_realn(bam1_t *b, char *ref)
  int bam_cap_mapQ(bam1_t *b, char *ref, int thres)

