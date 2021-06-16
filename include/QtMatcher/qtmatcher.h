/****************************************************************************
 *                                                                          *
 * Copyright (C) 2015 Neutrino International Inc.                           *
 *                                                                          *
 * Author : Brian Lin <lin.foxman@gmail.com>, Skype: wolfram_lin            *
 *                                                                          *
 ****************************************************************************/

#ifndef QT_MATCHER_H
#define QT_MATCHER_H

#include <QtCore>
#include <QtNetwork>
#include <QtSql>
#include <QtScript>
#include <Essentials>

QT_BEGIN_NAMESPACE

#ifndef QT_STATIC
#    if defined(QT_BUILD_QTMATCHER_LIB)
#      define Q_MATCHER_EXPORT Q_DECL_EXPORT
#    else
#      define Q_MATCHER_EXPORT Q_DECL_IMPORT
#    endif
#else
#    define Q_MATCHER_EXPORT
#endif

namespace N
{

/*****************************************************************************
 *                                                                           *
 *                        Approximate String Matching                        *
 *                                                                           *
 *****************************************************************************/

namespace Pattern
{

enum                                  {
  NDB_VERSION            = 150        , // Version number.
  NDB_NUM_TABLES         = 256        , // The number of hash tables.
  ENDIAN_BYTEORDER_CHECK = 0x62445371 , // A constant for byte-order checking.
  exact                  = 0          , // Exact match.
  dice                   = 1          , // Approximate string matching with dice coefficient.
  cosine                 = 2          , // Approximate string matching with cosine coefficient.
  jaccard                = 3          , // Approximate string matching with jaccard coefficient.
  overlap                = 4          , // Approximate string matching with overlap coefficient.
}                                     ;

class  Q_MATCHER_EXPORT Exact          ;
class  Q_MATCHER_EXPORT Dice           ;
class  Q_MATCHER_EXPORT Cosine         ;
class  Q_MATCHER_EXPORT Jaccard        ;
class  Q_MATCHER_EXPORT Overlap        ;
class  Q_MATCHER_EXPORT NGramGenerator ;
class  Q_MATCHER_EXPORT MurmurHash2    ;
struct Q_MATCHER_EXPORT ndbTableRef    ;
class  Q_MATCHER_EXPORT ndbException   ;

/* This class implements the traits of exact matching. */
class Q_MATCHER_EXPORT Exact
{
  public:

    int min_size  (int qsize,double alpha) ;
    int max_size  (int qsize,double alpha) ;
    int min_match (int qsize,int rsize,double alpha) ;

  protected:

  private:

};

/* This class implements the traits of dice coefficient. */
class Q_MATCHER_EXPORT Dice
{
  public:

    int min_size(int qsize, double alpha) ;
    int max_size(int qsize, double alpha) ;
    int min_match(int qsize, int rsize, double alpha) ;

  protected:

  private:

};

/* This class implements the traits of cosine coefficient. */
class Q_MATCHER_EXPORT Cosine
{
  public:

    int min_size(int qsize,double alpha) ;
    int max_size(int qsize,double alpha) ;
    int min_match(int qsize,int rsize,double alpha) ;

  protected:

  private:

};

/* This class implements the traits of Jaccard coefficient. */
class Q_MATCHER_EXPORT Jaccard
{
  public:

    int min_size  (int qsize,double alpha) ;
    int max_size  (int qsize,double alpha) ;
    int min_match (int qsize,int rsize,double alpha) ;

  protected:

  private:

};

/* This class implements the traits of overlap coefficient. */
class Q_MATCHER_EXPORT Overlap
{
  public:

    int min_size  (int qsize,double alpha) ;
    int max_size  (int qsize,double alpha) ;
    int min_match (int qsize,int rsize,double alpha) ;

  protected:

  private:

};

/**
 * Obtain a set of letter n-grams in a string.
 *  @param  str     The string.
 *  @param  ins     The insert iterator that receives the set of n-grams.
 *  @param  n       The unit of n-grams.
 *  @param  be      \c true to generate n-grams that encode begin and end of
 *                  a string.
 */
template < class string_type , class insert_iterator >
static void NGrams ( const string_type & str         ,
                     insert_iterator     ins         ,
                     int                 n           ,
                     bool                be          )
{
  typedef typename string_type::value_type   char_type                      ;
  typedef std::basic_stringstream<char_type> stringstream_type              ;
  typedef std::map<string_type, int>         ngram_stat_type                ;
  const   char_type                          mark = (char_type)0x01         ;
  string_type                                src                            ;

  if (be)                                                                   {
    // Append marks for begin/end of the string.
    for (int i = 0 ; i < n-1 ; ++i ) src += mark                            ;
    src += str                                                              ;
    for (int i = 0 ; i < n-1 ; ++i ) src += mark                            ;
  } else
  if ((int)str.length() < n)                                                {
    // Pad marks when the string is shorter than n.
     src = str                                                              ;
     for (int i = 0;i < n - (int)str.length();++i) src += mark              ;
  } else src = str                                                          ;
  // Count n-grams in the string.
  ngram_stat_type stat                                                      ;
  for (typename string_type::size_type i = 0 ; i < src.length()-n+1 ; ++i ) {
    string_type ngram = src.substr(i, n)                                    ;
    ++stat[ngram]                                                           ;
  }                                                                         ;
  // Convert the n-gram stat into a set.
  typename ngram_stat_type::const_iterator it                               ;
  for (it = stat.begin();it != stat.end();++it)                             {
    *ins = it->first                                                        ;
    // Append numbers if the same n-gram occurs more than once.
    for (int i = 2;i <= it->second;++i)                                     {
      stringstream_type ss                                                  ;
      ss << it->first << i                                                  ;
      *ins = ss.str()                                                       ;
    }                                                                       ;
  }                                                                         ;
}

/* This class generates n-grams for a string. */
class Q_MATCHER_EXPORT NGramGenerator
{
  public:

    /* Constructs an instance as a tri-gram generator. */
    NGramGenerator (void) : m_n  ( 3     )
                          , m_be ( false )
    {
    }

    /**
     * Constructs an instance as an n-gram generator.
     *  @param  n       The unit of n-grams.
     *  @param  be      \c true to generate n-grams that encode begin and
     *                  end of a string.
     */
    NGramGenerator (int n,bool be = false) : m_n  ( n  )
                                           , m_be ( be )
    {
    }

    /**
     * Sets the parameters for n-gram generation.
     *  @param  n       The unit of n-grams.
     *  @param  be      \c true to generate n-grams that encode begin and
     *                  end of a string.
     */
    void set(int n,bool be = false)
    {
      m_n  = n  ;
      m_be = be ;
    }

    /**
     * Gets the unit of n-grams.
     *  @return int     The unit of n-grams.
     */
    int get_n(void) const
    {
      return m_n ;
    }

    /**
     * Gets the flag for representing a begin/end of letters.
     *  @return bool    \c true if n-grams encoding the begin and end of a
     *                  string are generated.
     */
    bool get_be(void) const
    {
      return m_be ;
    }

    /**
     * Obtain a set of letter n-grams in a string.
     *  @param  str     The string.
     *  @param  ins     The insert iterator that receives the set of n-grams.
     */
    template < class string_type , class insert_iterator >
    void operator()(const string_type & str , insert_iterator ins) const
    {
      NGrams ( str , ins , m_n , m_be ) ;
    }

  protected:

    int  m_n  ; // The unit of n-grams.
    bool m_be ; // The flag for begin/end of tokens.

  private:

};

/**
 * Murmur Hash 2
 *
 *  This code makes the following assumption about how your machine behaves
 *      - We can read a 4-byte value from any address without crashing.
 *
 *  It also has a few limitations:
 *      - It will not work incrementally.
 *      - It will not produce the same results on little-endian and big-endian
 *        machines.
 */

class Q_MATCHER_EXPORT MurmurHash2 : public std::binary_function<const void *,size_t,unsigned int>
{
  public:

    inline unsigned int operator() (const void * key,size_t size) const
    {
      // 'm' and 'r' are mixing constants generated offline.
      // They're not really 'magic', they just happen to work well.
      const unsigned int m    = 0x5bd1e995        ;
      const int          r    = 24                ;
      // Initialize the hash to a 'random' value
      const unsigned int seed = 0x87654321        ;
            unsigned int h    = seed ^ size       ;
      // Mix 4 bytes at a time into the hash
      const char       * data = (const char *)key ;

      while ( size >= 4 )                         {
        unsigned int k = get32bits(data)          ;
        k    *= m                                 ;
        k    ^= k >> r                            ;
        k    *= m                                 ;
        h    *= m                                 ;
        h    ^= k                                 ;
        data += 4                                 ;
        size -= 4                                 ;
      }
      // Handle the last few bytes of the input array
      switch ( size )                             {
        case 3: h ^= data[2] << 16                ;
        case 2: h ^= data[1] <<  8                ;
        case 1: h ^= data[0]                      ;
                h *= m                            ;
      }                                           ;
      // Do a few final mixes of the hash to ensure the last few
      // bytes are well-incorporated.
      h ^= h >> 13                                ;
      h *= m                                      ;
      h ^= h >> 15                                ;
      return h                                    ;
    }

  protected:

    inline static unsigned int get32bits(const char * d)
    {
      return *reinterpret_cast<const unsigned int *>(d) ;
    }

  private:

};

struct ndbTableRef
{
  unsigned int offset ;
  unsigned int num    ;
};

Q_MATCHER_EXPORT unsigned int TableBegin(void) ;

/* Exception class for the NDB builder. */
class Q_MATCHER_EXPORT ndbException : public std::invalid_argument
{
  public:

    ndbException(const std::string & msg) : std::invalid_argument ( msg )
    {
    }

};

/* NDB builder */
template <typename hash_function>
class ndbBuilder
{
  protected:

    // A bucket structure.
    struct bucket {
      unsigned int hash   ; // Hash value of the record.
      unsigned int offset ; // Offset address to the actual record.

      bucket(void) : hash   ( 0 )
                   , offset ( 0 )
      {
      }

      bucket(unsigned int h,unsigned int o) : hash   ( h )
                                            , offset ( o )
      {
      }

    };

    // A hash table is a vector of buckets.
    typedef std::vector<bucket> hashtable    ;

    std::ofstream &  m_os                    ; // Output stream.
    unsigned int     m_begin                 ;
    unsigned int     m_cur                   ;
    hashtable        m_ht [ NDB_NUM_TABLES ] ; // Hash tables.

  public:

    /**
     * Constructs an object.
     *  @param  os          The output stream to which this class write the
     *                      database. This stream must be opened in the
     *                      binary mode (\c std::ios_base::binary).
     */
    ndbBuilder ( std::ofstream & os ) : m_os ( os )
    {
      m_begin = m_os.tellp (                 ) ;
      m_cur   = TableBegin (                 ) ;
      m_os.seekp           ( m_begin + m_cur ) ;
    }

    /**
     * Destructs an object.
     */
    virtual ~ndbBuilder(void)
    {
      this -> close ( ) ;
    }

    /**
     * Inserts a pair of key and value to the database.
     *  Any key in the database should be unique, but this library does not
     *  check duplicated keys.
     *  @param  key         The pointer to the key.
     *  @param  ksize       The size of the key.
     *  @param  value       The pointer to the value.
     *  @param  vsize       The size of the value.
     */
    template <class key_t, class value_t>
    void put(const key_t * key, size_t ksize, const value_t * value, size_t vsize)
    {
      // Write out the current record.
      write_uint32 ( (unsigned int)ksize                           ) ;
      m_os.write   ( reinterpret_cast<const char *>(key  ) , ksize ) ;
      write_uint32 ( (unsigned int)vsize                           ) ;
      m_os.write   ( reinterpret_cast<const char *>(value) , vsize ) ;
      // Compute the hash value and choose a hash table.
      unsigned int hv = hash_function()(static_cast<const void *>(key),ksize) ;
      hashtable  & ht = m_ht [ hv % NDB_NUM_TABLES ]                          ;
      // Store the hash value and offset to the hash table.
      ht.push_back(bucket(hv, m_cur)) ;
      // Increment the current position.
      m_cur += sizeof(unsigned int) + ksize + sizeof(unsigned int) + vsize    ;
    }

  protected:

    void close(void)
    {
      // Check the consistency of the stream offset.
      if (m_begin + m_cur != (unsigned int)m_os.tellp()) {
        throw ndbException("Inconsistent stream offset") ;
      }                                                  ;
      // Store the hash tables. At this moment, the file pointer refers to
      // the offset succeeding the last key/value pair.
      for (size_t i = 0 ; i < NDB_NUM_TABLES ; ++i ) {
        hashtable & ht = m_ht [ i ]                  ;
        // Do not write an empty hash table.
        if (!ht.empty())                             {
          // An actual table will have the double size; half elements
          // in the table are kept empty.
          int n = ht.size() * 2                      ;
          // Allocate the actual table.
          bucket * dst = new bucket[n]               ;
          // Put hash elements to the table with the open-address method.
          typename hashtable::const_iterator it      ;
          for (it = ht.begin();it != ht.end();++it)  {
            int k = (it->hash >> 8) % n              ;
            // Find a vacant element.
            while (dst[k].offset != 0)               {
              k = (k+1) % n                          ;
            }                                        ;
            // Store the hash element.
            dst[k].hash   = it->hash                 ;
            dst[k].offset = it->offset               ;
          }
          // Write out the new table.
          for (int k = 0;k < n;++k)                  {
            write_uint32 ( dst[k] . hash   )         ;
            write_uint32 ( dst[k] . offset )         ;
          }
          // Free the table.
          delete[] dst                               ;
        }
      }
      // Store the current position.
      unsigned int offset = (unsigned int)m_os.tellp() ;
      // Rewind the stream position to the beginning.
      m_os.seekp   ( m_begin)                            ;
      // Write the file header.
      char chunkid[4] = {'N','C','D','B'}            ;
      m_os.write   ( chunkid , 4            )        ;
      write_uint32 ( offset  - m_begin      )        ;
      write_uint32 ( NDB_VERSION            )        ;
      write_uint32 ( ENDIAN_BYTEORDER_CHECK )        ;
      // Write references to hash tables. At this moment, dbw->cur points
      // to the offset succeeding the last key/data pair.
      for (size_t i = 0;i < NDB_NUM_TABLES;++i)      {
        // Offset to the hash table (or zero for non-existent tables).
        write_uint32 ( m_ht[i].empty() ? 0 : m_cur ) ;
        // Bucket size is double to the number of elements.
        write_uint32 ( m_ht[i].size () * 2         ) ;
        // Advance the offset counter.
        m_cur += sizeof(unsigned int) * 2 * m_ht[i].size() * 2 ;
      }
      // Seek to the last position.
      m_os.seekp ( offset )                          ;
    }

    inline void write_uint32(unsigned int value)
    {
      m_os.write(reinterpret_cast<const char *>(&value), sizeof(value)) ;
    }

};

/* NDB++ reader */
template <typename hash_function>
class ndbReader
{
  protected:

    struct bucket_t
    {
      unsigned int hash   ; // Hash value of the record.
      unsigned int offset ; // Offset address to the actual record.
    } ;

    struct hashtable_t
    {
      unsigned int     num     ; // Number of elements in the table.
      const bucket_t * buckets ; // Buckets (array of bucket).
    } ;

    const unsigned char * m_buffer                ; // Pointer to the memory block.
    size_t                m_size                  ; // Size of the memory block.
    bool                  m_own                   ; //
    hashtable_t           m_ht [ NDB_NUM_TABLES ] ; // Hash tables.
    size_t                m_n                     ;

  public:

    /* Constructs an object. */
    ndbReader(void) : m_buffer ( NULL  )
                    , m_size   ( 0     )
                    , m_own    ( false )
                    , m_n      ( 0     )
    {
    }

    /**
     * Constructs an object by opening a database on memory.
     *  @param  buffer      The pointer to the memory image of the database.
     *  @param  size        The size of the memory image.
     *  @param  own         If this is set to \c true, this library will call
     *                      delete[] when the database is closed.
     */
    ndbReader(const void * buffer,size_t size,bool own) : m_buffer ( NULL  )
                                                        , m_size   ( 0     )
                                                        , m_own    ( false )
                                                        , m_n      ( 0     )
    {
      this -> open ( buffer , size , own ) ;
    }

    /**
     * Constructs an object by opening a database from an input stream.
     *  @param  ifs         The input stream from which this library reads
     *                      a database.
     */
    ndbReader ( std::ifstream & ifs ) : m_buffer ( NULL  )
                                      , m_size   ( 0     )
                                      , m_own    ( false )
                                      , m_n      ( 0     )
    {
      this -> open ( ifs ) ;
    }

    /**
     * Destructs the object.
     */
    virtual ~ndbReader(void)
    {
      this -> close ( ) ;
    }

    /**
     * Tests if the database is opened.
     *  @return bool        \c true if the database is opened,
     *                      \c false otherwise.
     */
    bool is_open(void) const
    {
      return ( m_buffer != NULL ) ;
    }

    /**
     * Obtains the number of elements in the database.
     *  @return size_t          The number of elements.
     */
    size_t size(void) const
    {
      return m_n ;
    }

    /**
     * Tests if the database is empty.
     *  @return bool        \c true if the number of records is zero,
     *                      \c false otherwise.
     */
    bool empty(void) const
    {
      return ( m_n == 0 ) ;
    }

    /**
     * Opens the database from an input stream.
     *  @param  ifs         The input stream from which this library reads
     *                      a database.
     */
    size_t open(std::ifstream & ifs)
    {
      char                   chunk [ 4 ]           ;
      char                   size  [ 4 ]           ;
      std::istream::pos_type offset = ifs.tellg()  ;
      do                                           {
        // Read a chunk identifier.
        ifs . read ( chunk , 4 )                   ;
        if (ifs.fail()) break                      ;
        // Check the chunk identifier.
        if (std::strncmp(chunk, "NCDB", 4) != 0)   {
          break                                    ;
        }                                          ;
        // Read the size of the chunk.
        ifs . read ( size , 4 )                    ;
        if (ifs.fail()) break                      ;
        // Allocate a memory block for the chunk.
        unsigned int    chunk_size = read_uint32(reinterpret_cast<unsigned char *>(size)) ;
        unsigned char * block      = new unsigned char [ chunk_size ] ;
        // Read the memory image from the stream.
        ifs . seekg ( 0 , std::ios_base::beg )     ;
        if (ifs.fail()) break                      ;
        ifs.read(reinterpret_cast<char *>(block),chunk_size) ;
        if (ifs.fail()) break                      ;
        return this->open(block, chunk_size, true) ;
      } while ( 0 )                                ;
      ifs . seekg ( offset , std::ios::beg )       ;
      return 0                                     ;
    }

    /**
     * Opens the database from a memory image.
     *  @param  buffer      The pointer to the memory image of the database.
     *  @param  size        The size of the memory image.
     *  @param  own         If this is set to \c true, this library will call
     *                      delete[] when the database is closed.
     */
    size_t open(const void * buffer, size_t size, bool own = false)
    {
      const unsigned char * p = reinterpret_cast<const unsigned char *>(buffer) ;
      // Make sure that the size of the chunk is larger than the minimum size.
      if (size < TableBegin())                       {
        throw ndbException("The memory image is smaller than a chunk header.") ;
      }
      // Check the chunk identifier.
      if (memcmp ( p , "NCDB" , 4 ) != 0)            {
        throw ndbException("Incorrect chunk header") ;
      }
      p += 4                                         ;
      // Read the chunk header.
      unsigned int csize     = read_uint32(p)        ;
      p += sizeof(unsigned int)                      ;
      unsigned int version   = read_uint32(p)        ;
      p += sizeof(unsigned int)                      ;
      unsigned int byteorder = read_uint32(p)        ;
      p += sizeof(unsigned int)                      ;
      // Check the byte-order consistency.
      if ( byteorder != ENDIAN_BYTEORDER_CHECK )     {
        throw ndbException("Inconsistent byte order") ;
      }
      // Check the version number.
      if ( version != NDB_VERSION )                  {
        throw ndbException("Incompatible CDB++ versions");
      }
      // Check the chunk size.
      if ( size < csize )                            {
        throw ndbException("The memory image is smaller than a chunk size.");
      }
      // Set memory block and size.
      m_buffer = reinterpret_cast<const unsigned char *>(buffer) ;
      m_size   = size                                ;
      m_own    = own                                 ;
      // Set pointers to the hash tables.
      m_n      = 0                                   ;
      const ndbTableRef * ref = reinterpret_cast<const ndbTableRef *>(p) ;
      for ( size_t i = 0 ; i < NDB_NUM_TABLES ; ++i ) {
        if ( ref[i] . offset )                       {
          // Set the buckets.
          m_ht[i].buckets = reinterpret_cast<const bucket_t*>(m_buffer + ref[i].offset);
          m_ht[i].num     = ref[i].num               ;
        } else                                       {
          // An empty hash table.
          m_ht[i].buckets = NULL                     ;
          m_ht[i].num     = 0                        ;
        }
        // The number of records is the half of the table size.
        m_n += ( ref[i] . num / 2 )                  ;
      }
      return (size_t)csize ;
    }

    /* Closes the database. */
    void close(void)
    {
      if (m_own && m_buffer != NULL) {
        delete[] m_buffer            ;
      }
      m_buffer = NULL                ;
      m_size   = 0                   ;
      m_n      = 0                   ;
    }

    /**
     * Finds the key in the database.
     *  @param  key         The pointer to the key.
     *  @param  ksize       The size of the key.
     *  @param  vsize       The pointer of a variable to which the size of the
     *                      value returned. This parameter can be \c NULL.
     *  @return const void* The pointer to the value.
     */
    const void * get(const void * key, size_t ksize, size_t * vsize) const
    {
      unsigned int        hv = hash_function()(key, ksize)   ;
      const hashtable_t * ht = &m_ht [ hv % NDB_NUM_TABLES ] ;
      if ( ht->num && ht->buckets != NULL )                  {
        int              n = ht->num                         ;
        int              k = (hv >> 8) % n                   ;
        const bucket_t * p = NULL                            ;
        while ( p = &ht->buckets[k] , p->offset )            { // ?? lamba expression
          if ( p->hash == hv )                               {
            const unsigned char * q = m_buffer + p->offset   ;
            if ( read_uint32(q) == ksize &&
                 memcmp(key, q + sizeof(unsigned int), ksize) == 0) {
              q += sizeof(unsigned int) + ksize              ;
              if (vsize != NULL) *vsize = read_uint32(q)     ;
              return q + sizeof(unsigned int)                ;
            }                                                ;
          }                                                  ;
          k = (k+1) % n                                      ;
        }                                                    ;
      }                                                      ;
      if (vsize != NULL) *vsize = 0                          ;
      return NULL                                            ;
    }

  protected:

    inline unsigned int read_uint32(const unsigned char * p) const
    {
      return *reinterpret_cast<const unsigned int *>(p) ;
    }

};

typedef ndbBuilder < MurmurHash2 > Builder ;
typedef ndbReader  < MurmurHash2 > Reader  ;

/**
 * A writer for an n-gram database.
 *  This template class builds an n-gram database. The first template
 *  argument (string_tmpl) specifies the type of a key (string), the second
 *  template argument (value_tmpl) specifies the type of a value associated
 *  with a key, and the third template argument (ngram_generator_tmpl)
 *  customizes generation of feature sets (n-grams) from keys.
 *
 *  This class is inherited by writer_base, which adds the functionality of
 *  managing a master string table (list of strings).
 *
 *  @param  string_tmpl             The type of a string.
 *  @param  value_tmpl              The value type.
 *                                  This is required to be an integer type.
 *  @param  ngram_generator_tmpl    The type of an n-gram generator.
 */
template < class string_tmpl , class value_tmpl , class ngram_generator_tmpl >
class ngramdb_writer_base
{
  public:

    typedef string_tmpl                      string_type          ; // The type representing a string.
    typedef value_tmpl                       value_type           ; // The type of values associated with key strings.
    typedef ngram_generator_tmpl             ngram_generator_type ; // The function type for generating n-grams from a key string.
    typedef typename string_type::value_type char_type            ; // The type representing a character.

  protected:

    typedef std::vector<string_type>          ngrams_type  ; // The type of an array of n-grams.
    typedef std::vector<value_type>           values_type  ; // The vector type of values associated with an n-gram.
    typedef std::map<string_type,values_type> hashdb_type  ; // The type implementing an index (associations from n-grams to values).
    typedef std::vector<hashdb_type>          indices_type ; // The vector of indices for different n-gram sizes.

    indices_type                 m_indices ; // The vector of indices.
    const ngram_generator_type & m_gen     ; // The n-gram generator.
    std::stringstream            m_error   ; // The error message.

  public:

    /**
     * Constructs an object.
     *  @param  gen             The n-gram generator.
     */
    ngramdb_writer_base(const ngram_generator_type & gen) : m_gen ( gen )
    {
    }

    /* Destructs an object. */
    virtual ~ngramdb_writer_base(void)
    {
    }

    /* Clears the database. */
    void clear(void)
    {
      m_indices . clear (    ) ;
      m_error   . str   ( "" ) ;
    }

    /**
     * Checks whether the database is empty.
     *  @return bool    \c true if the database is empty, \c false otherwise.
     */
    bool empty(void)
    {
      return m_indices . empty ( ) ;
    }

    /**
     * Returns the maximum length of keys in the n-gram database.
     *  @return int     The maximum length of keys.
     */
    int max_size(void) const
    {
      return (int)m_indices . size ( ) ;
    }

    /**
     * Checks whether an error has occurred.
     *  @return bool    \c true if an error has occurred.
     */
    bool fail(void) const
    {
      return ! m_error . str ( ) . empty ( ) ;
    }

    /**
     * Returns an error message.
     *  @return std::string The string of the error message.
     */
    std::string error(void) const
    {
      return m_error . str ( ) ;
    }

    /**
     * Inserts a string to the n-gram database.
     *  @param  key         The key string.
     *  @param  value       The value associated with the string.
     */
    bool insert(const string_type & key,const value_type & value)
    {
      // Generate n-grams from the key string.
      ngrams_type ngrams                     ;
      m_gen(key, std::back_inserter(ngrams)) ;
      if (ngrams.empty()) return false       ;
      // Resize the index array for the number of the n-grams;
      // we build an index for each n-gram number.
      if (m_indices.size() < ngrams.size())  {
        m_indices.resize(ngrams.size())      ;
      }
      hashdb_type & index = m_indices[ngrams.size()-1] ;
      // Store the associations from the n-grams to the value.
      typename ngrams_type::const_iterator it ;
      for (it = ngrams.begin();it != ngrams.end();++it) {
        const string_type & ngram = *it      ;
        typename hashdb_type::iterator iti = index.find(ngram);
        if (iti == index.end()) {
          // Create a new posting array.
          values_type v(1);
          v[0] = value;
          index.insert(typename hashdb_type::value_type(ngram,v));
        } else {
          // Append the value to the existing posting array.
          iti->second.push_back(value);
        }
      }
      return true ;
    }

    /**
     * Stores the n-gram database to files.
     *  @param  name        The prefix of file names.
     *  @return bool        \c true if the database is successfully stored,
     *                      \c false otherwise.
     */
    bool store(const std::string & base)
    {
      // Write out all the indices to files.
      for (int i = 0;i < (int)m_indices.size();++i)    {
        if (!m_indices[i].empty())                     {
          std::stringstream ss                         ;
          ss << base << '.' << i+1 << ".cdb"           ;
          bool b = this->store(ss.str(), m_indices[i]) ;
          if (!b) return false                         ;
        }                                              ;
      }                                                ;
      return true                                      ;
    }

  protected:

    bool store(const std::string & name, const hashdb_type & index)
    {
      // Open the database file with binary mode.
      std::ofstream ofs(name.c_str(), std::ios::binary)          ;
      if (ofs.fail())                                            {
        m_error << "Failed to open a file for writing: " << name ;
        return false                                             ;
      }

      try                                                        {
        // Open a NDB writer
        Builder dbw(ofs)                                         ;
        // Put associations: n-gram -> values.
        typename hashdb_type::const_iterator it                  ;
        for (it = index.begin();it != index.end();++it)          {
          // Put an association from an n-gram to its values.
          dbw . put                                              (
            it->first.c_str()                                    ,
            sizeof(char_type) * it->first.length()               ,
            &it->second[0]                                       ,
            sizeof(it->second[0]) * it->second.size()          ) ;
        }
      } catch (const ndbException & e)                           {
        m_error << "NDB error: " << e.what()                     ;
        return false                                             ;
      }
      return true                                                ;
    }

};

/**
 * A SimString database writer.
 *  This template class builds a SimString database. The first template
 *  argument (string_tmpl) specifies the type of a character, and the second
 *  template argument (ngram_generator_tmpl) customizes generation of feature
 *  sets (n-grams) from strings.
 *
 *  Inheriting the base class ngramdb_writer_base that builds indices from
 *  n-grams to string IDs, this class maintains associations between strings
 *  and string IDs.
 *
 *  @param  string_tmpl             The type of a string.
 *  @param  ngram_generator_tmpl    The type of an n-gram generator.
 */
template < class string_tmpl , class ngram_generator_tmpl = NGramGenerator >
class writer_base : public ngramdb_writer_base < string_tmpl, unsigned int , ngram_generator_tmpl >
{
  public:

    typedef string_tmpl string_type                               ; // The type representing a string.
    typedef unsigned int value_type                               ; // The type of values associated with key strings.
    typedef ngram_generator_tmpl             ngram_generator_type ; // The function type for generating n-grams from a key string.
    typedef typename string_type::value_type char_type            ; // The type representing a character.
    typedef ngramdb_writer_base<string_tmpl,unsigned int,ngram_generator_tmpl> base_type; // The type of the base class.

  protected:

    std::string   m_name        ; // The base name of the database.
    std::ofstream m_ofs         ; // The output stream for the string collection.
    int           m_num_entries ; // The number of strings in the database.

  public:

    /**
     * Constructs a writer object.
     *  @param  gen         The n-gram generator used by this writer.
     */
    writer_base(const ngram_generator_type & gen) : base_type     ( gen )
                                                  , m_num_entries ( 0   )
    {
    }

    /**
     * Constructs a writer object by opening a database.
     *  @param  gen         The n-gram generator used by this writer.
     *  @param  name        The name of the database.
     */
    writer_base(const ngram_generator_type & gen,const std::string & name)
              : base_type     ( gen )
              , m_num_entries ( 0   )
    {
      this -> open ( name ) ;
    }

    /* Destructs a writer object. */
    virtual ~writer_base(void)
    {
      this -> close ( ) ;
    }

    /**
     * Opens a database.
     *  @param  name        The name of the database.
     *  @return bool        \c true if the database is successfully opened,
     *                      \c false otherwise.
     */
    bool open(const std::string & name)
    {
      m_num_entries = 0;
      // Open the master file for writing.
      m_ofs.open(name.c_str(), std::ios::binary) ;
      if (m_ofs.fail()) {
        this->m_error << "Failed to open a file for writing: " << name;
        return false;
      }
      // Reserve the region for a file header.
      if (!this->write_header(m_ofs)) {
        m_ofs.close();
        return false;
      }
      m_name = name;
      return true;
    }

    /**
     * Closes the database.
     *  @param  name        The name of the database.
     *  @return bool        \c true if the database is successfully opened,
     *                      \c false otherwise.
     */
    bool close(void)
    {
      bool b = true;
      // Write the n-gram database to files.
      if (!m_name.empty()) {
        b &= this->store(m_name);
      }
      // Finalize the file header, and close the file.
      if (m_ofs.is_open()) {
        b &= this->write_header(m_ofs);
        m_ofs.close();
      }
      // Initialize the members.
      m_name.clear();
      m_num_entries = 0;
      return b;
    }

    /**
     * Inserts a string to the database.
     *  @param  str         The string to be inserted.
     *  @return bool        \c true if the string is successfully inserted,
     *                      \c false otherwise.
     */
    bool insert(const string_type & str)
    {
      // This will be the offset address to access the key string.
      value_type off = (value_type)(std::streamoff)m_ofs.tellp() ;
      // Write the key string to the master file.
      m_ofs.write(reinterpret_cast<const char*>(str.c_str()), sizeof(char_type) * (str.length()+1));
      if (m_ofs.fail()) {
        this->m_error << "Failed to write a string to the master file." ;
        return false;
      }
      ++m_num_entries;
      // Insert the n-grams of the key string to the database.
      return base_type::insert(str,off);
    }

  protected:

    bool write_header(std::ofstream & ofs)
    {
      unsigned int num_entries = m_num_entries;
      unsigned int max_size    = (unsigned int)this->max_size();
      unsigned int size        = (unsigned int)m_ofs.tellp();
      // Seek to the beginning of the master file, to which the file header
      // is to be written.
      ofs.seekp(0);
      if (ofs.fail()) {
        this->m_error << "Failed to seek the file pointer for the master file.";
        return false;
      }
      // Write the file header.
      m_ofs.write  ("CIDB", 4) ;
      write_uint32 (ENDIAN_BYTEORDER_CHECK);
      write_uint32 (NDB_VERSION);
      write_uint32 (size);
      write_uint32 (sizeof(char_type));
      write_uint32 (this->m_gen.get_n());
      write_uint32 (static_cast<int>(this->m_gen.get_be()));
      write_uint32 (num_entries);
      write_uint32 (max_size);
      if (ofs.fail()) {
        this->m_error << "Failed to write a file header to the master file.";
        return false ;
      }
      return true;
    }

    inline void write_uint32(unsigned int value)
    {
      m_ofs.write(reinterpret_cast<const char *>(&value), sizeof(value)) ;
    }

};

#ifdef XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

#include <iostream>
#include <string>

class memory_mapped_file_base
{
public:
    typedef size_t size_type;

    memory_mapped_file_base() {}
    virtual ~memory_mapped_file_base() {}

    void open(const std::string& path, std::ios_base::openmode mode) {}
    bool is_open() const {return false; }
    void close() {}
    void resize(size_type size) {}
    size_type size() const {return 0; }
    char* data() const {return NULL; }
    const char* const_data() const {return NULL; }
    static int alignment() {return 0; }
};

#if     defined(_WIN32)
#include "memory_mapped_file_win32.h"
#define memory_mapped_file memory_mapped_file_win32

#else
#include "memory_mapped_file_posix.h"
#define memory_mapped_file memory_mapped_file_posix

#endif

#define __MEMORY_MAPPED_FILE_POSIX_H__

#include <fcntl.h>
#include <memory.h>
#include <stdio.h>
#include <stdlib.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <sys/mman.h>

class memory_mapped_file_posix :
    public memory_mapped_file_base
{
public:
    typedef size_t size_type;

protected:
    int                     m_fd;
    std::ios_base::openmode m_mode;
    void*                   m_data;
    size_type               m_size;

public:
    memory_mapped_file_posix()
    {
        m_fd = -1;
        m_mode = std::ios_base::in;
        m_data = NULL;
        m_size = 0;
    }

    virtual ~memory_mapped_file_posix()
    {
        close();
    }

    void open(const std::string& path, std::ios_base::openmode mode)
    {
        int flags = 0;
        struct stat buf;

        if (mode & std::ios_base::in) {
            flags = O_RDONLY;
        }
        if (mode & std::ios_base::out) {
            flags = O_RDWR | O_CREAT;
        }
        if (mode & std::ios_base::trunc) {
            flags |= (O_RDWR | O_TRUNC);
        }

        m_fd = ::open(path.c_str(), flags, S_IRUSR | S_IWUSR | S_IRGRP | S_IROTH);
        if (m_fd != -1) {
            if (::fstat(m_fd, &buf) == 0) {
                m_mode = mode;
                this->resize((size_type)buf.st_size);
            } else {
                ::close(m_fd);
                m_fd = -1;
            }
        }
    }

    bool is_open() const
    {
        return (m_fd != -1);
    }

    void close()
    {
        this->free();
        if (m_fd != -1) {
            ::close(m_fd);
            m_fd = -1;
        }
    }

    bool resize(size_type size)
    {
        if (size == 0) {
            this->free();
            return true;
        }

        if (m_fd == -1) {
            return false;
        }

        this->free();

        if ((m_mode & std::ios_base::out) && m_size < size) {
            /* Try to expand the file to the specified size. */
            if (::lseek(m_size, size, SEEK_SET) >= 0) {
                char c;
                if (read(m_fd, &c, sizeof(char)) == -1) {
                    c = 0;
                }
                if (write(m_fd, &c, sizeof(char)) == -1) {
                    return false;   // Failed to write the last position.
                }
            } else {
                return false;       // Failed to expand the file.
            }
        }

        /* Map the file into process memory. */
        m_data = ::mmap(
            NULL,
            size,
            (m_mode & std::ios_base::out) ? (PROT_READ | PROT_WRITE) : PROT_READ,
            MAP_SHARED,
            m_fd,
            0);

        m_size = size;
        return true;
    }

    void free()
    {
        if (m_data != NULL) {
            ::munmap(m_data, m_size);
            m_data = NULL;
        }
        m_size = 0;
    }

    size_type size() const
    {
        return m_size;
    }

    char* data() const
    {
        return reinterpret_cast<char*>(m_data);
    }

    const char* const_data() const
    {
        return reinterpret_cast<const char*>(m_data);
    }

    static int alignment()
    {
        return 0;
    }
};

/*__MEMORY_MAPPED_FILE_POSIX_H__*/

// __MEMORY_MAPPED_FILE_WIN32_H__

#define NOMINMAX    // To fix min/max conflicts with STL.

#include <memory.h>
#include <windows.h>

class memory_mapped_file_win32 :
    public memory_mapped_file_base
{
public:
    typedef size_t size_type;

protected:
        HANDLE	                m_hFile;
        HANDLE	                m_hMapping;
    std::ios_base::openmode m_mode;
        char*	                m_data;
        size_type               m_size;

public:
    memory_mapped_file_win32()
    {
        m_hFile = INVALID_HANDLE_VALUE;
        m_hMapping = INVALID_HANDLE_VALUE;
        m_mode = 0;
        m_data = NULL;
        m_size = 0;
    }

    virtual ~memory_mapped_file_win32()
    {
        close();
    }

    void open(const std::string& path, std::ios_base::openmode mode)
    {
                DWORD dwDesiredAccess = 0;
                DWORD dwCreationDisposition = 0;

        if (mode & std::ios_base::in) {
                        dwDesiredAccess |= GENERIC_READ;
                        dwCreationDisposition = OPEN_EXISTING;
                }
        if (mode & std::ios_base::out) {
                        dwDesiredAccess |= GENERIC_WRITE;
                        dwCreationDisposition = CREATE_NEW;
                }
                if (mode & std::ios_base::trunc) {
                        dwDesiredAccess = (GENERIC_READ | GENERIC_WRITE);
                        dwCreationDisposition = CREATE_ALWAYS;
                }

                m_hFile = CreateFileA(
                        path.c_str(),
                        dwDesiredAccess,
                        0,
                        NULL,
                        dwCreationDisposition,
                        FILE_ATTRIBUTE_NORMAL,
                        NULL
                        );

                if (m_hFile != INVALID_HANDLE_VALUE) {
            m_mode = mode;
            this->resize((size_type)GetFileSize(m_hFile, NULL));
        }
        }

    bool is_open() const
    {
        return (m_hFile != INVALID_HANDLE_VALUE);
    }

    void close()
    {
        this->free();
        if (m_hFile != INVALID_HANDLE_VALUE) {
                    CloseHandle(m_hFile);
                    m_hFile = INVALID_HANDLE_VALUE;
        }
        }

    bool resize(size_type size)
    {
            if (size == 0) {
            this->free();
                    return true;
            }

        if (m_hFile == INVALID_HANDLE_VALUE) {
            return false;
        }

        this->free();
        DWORD flProtect = (m_mode & std::ios_base::out) ? PAGE_READWRITE : PAGE_READONLY;
        m_hMapping = CreateFileMappingA(
                    m_hFile,
                    NULL,
            flProtect,
                    0,
                    (DWORD)size,
                    NULL
                    );

            if (m_hMapping == NULL) {
                    CloseHandle(m_hFile);
                    m_hFile = NULL;
            return false;
            }

        DWORD dwDesiredAccess = (m_mode & std::ios_base::out) ? FILE_MAP_ALL_ACCESS : FILE_MAP_READ;
        m_data = (char*)MapViewOfFile(
                    m_hMapping,
                    dwDesiredAccess,
                    0,
                    0,
                    0
                    );

            if (m_data == NULL) {
                    CloseHandle(m_hMapping);
                    m_hMapping = NULL;
                    CloseHandle(m_hFile);
                    m_hFile = NULL;
            return false;
            }

            m_size = size;
        return true;
    }

    void free()
    {
            if (m_data != NULL) {
                    UnmapViewOfFile(m_data);
                    m_data = NULL;
            }
            if (m_hMapping != INVALID_HANDLE_VALUE) {
                    CloseHandle(m_hMapping);
                    m_hMapping = NULL;
            }
            m_size = 0;
    }

    size_type size() const
    {
        return m_size;
    }

    char* data() const
    {
        return m_data;
    }

    const char* const_data() const
    {
        return m_data;
    }

    static int alignment()
    {
        return 0;
    }
};

/*__MEMORY_MAPPED_FILE_WIN32_H__*/

/**
 * A reader for an n-gram database.
 *  @param  value_tmpl              The value type.
 *                                  This is required to be an integer type.
 */
template <
    class value_tmpl
>
class ngramdb_reader_base
{
public:
    /// The type of a value.
    typedef value_tmpl value_type;

protected:
    // An inverted list of SIDs.
    struct inverted_list_type
    {
        int num;
        const value_type* values;

        friend bool operator<(
            const inverted_list_type& x,
            const inverted_list_type& y
            )
        {
            return (x.num < y.num);
        }
    };
    // An array of inverted lists.
    typedef std::vector<inverted_list_type> inverted_lists_type;

    // A hash table that retrieves SIDs from n-grams.
    typedef cdbpp::cdbpp hashtbl_type;

    // An index containing strings of the same size.
    struct index_type
    {
        // The memory image of the database.
        memory_mapped_file  image;
        // The index.
        hashtbl_type        table;
    };

    // Indices with different sizes of strings.
    typedef std::vector<index_type> indices_type;

    // A candidate string of retrieved results.
    struct candidate_type
    {
        // The SID.
        value_type  value;
        // The overlap count (frequency of the SID in the inverted lists).
        int         num;

        candidate_type(value_type v, int n)
            : value(v), num(n)
        {
        }
    };

    // An array of candidates.
    typedef std::vector<candidate_type> candidates_type;

    // An array of SIDs retrieved.
    typedef std::vector<value_type> results_type;

protected:
    // The array of the indices.
    indices_type m_indices;
    // The maximum size of strings in the database.
    int m_max_size;
    // The database name (base name of indices).
    std::string m_name;
    // The error message.
    std::stringstream m_error;


public:
    /**
     * Constructs an object.
     */
    ngramdb_reader_base()
    {
    }

    /**
     * Destructs an object.
     */
    virtual ~ngramdb_reader_base()
    {
    }

    /**
     * Checks whether an error has occurred.
     *  @return bool    \c true if an error has occurred.
     */
    bool fail() const
    {
        return !m_error.str().empty();
    }

    /**
     * Returns an error message.
     *  @return std::string The string of the error message.
     */
    std::string error() const
    {
        return m_error.str();
    }

    /**
     * Opens an n-gram database.
     *  @param  name        The name of the database.
     *  @param  max_size    The maximum size of the strings.
     */
    void open(const std::string& name, int max_size)
    {
        m_name = name;
        m_max_size = max_size;
        // The maximum size corresponds to the number of indices in the database.
        m_indices.resize(max_size);
    }

    /**
     * Closes an n-gram database.
     */
    void close()
    {
        m_name.clear();
        m_indices.clear();
        m_error.str("");
    }

    /**
     * Performs an overlap join on inverted lists retrieved for the query.
     *  @param  query       The query object that stores query n-grams,
     *                      threshold, and conditions for the similarity
     *                      measure.
     *  @param  results     The SIDs that satisfies the overlap join.
     */
    template <class measure_type, class query_type>
    void overlapjoin(const query_type& query, double alpha, results_type& results)
    {
        int i;
        const int qsize = query.size();

        // Allocate a vector of postings corresponding to n-gram queries.
        inverted_lists_type posts(qsize);

        // Compute the range of n-gram lengths for the candidate strings;
        // in other words, we do not have to search for strings whose n-gram
        // lengths are out of this range.
        const int xmin = std::max(measure_type::min_size(query.size(), alpha), 1);
        const int xmax = std::min(measure_type::max_size(query.size(), alpha), m_max_size);

        // Loop for each length in the range.
        for (int xsize = xmin;xsize <= xmax;++xsize) {
            // Access to the n-gram index for the length.
            hashtbl_type& tbl = open_index(m_name, xsize);
            if (!tbl.is_open()) {
                // Ignore an empty index.
                continue;
            }

            // Search for string entries that match to each query n-gram.
            // Note that we do not traverse each entry here, but only obtain
            // the number of and the pointer to the entries.
            typename query_type::const_iterator it;
            for (it = query.begin(), i = 0;it != query.end();++it, ++i) {
                size_t vsize;
                const void *values = tbl.get(
                    it->c_str(),
                    sizeof(it->at(0)) * it->length(),
                    &vsize
                    );
                posts[i].num = (int)(vsize / sizeof(value_type));
                posts[i].values = reinterpret_cast<const value_type*>(values);
            }

            // Sort the query n-grams by ascending order of their frequencies.
            // This reduces the number of initial candidates.
            std::sort(posts.begin(), posts.end());

            // The minimum number of n-gram matches required for the query.
            const int mmin = measure_type::min_match(qsize, xsize, alpha);
            // A candidate must match to one of n-grams in these queries.
            const int min_queries = qsize - mmin + 1;

            // Step 1: collect candidates that match to the initial queries.
            candidates_type cands;
            for (i = 0;i < min_queries;++i) {
                candidates_type tmp;
                typename candidates_type::const_iterator itc = cands.begin();
                const value_type* p = posts[i].values;
                const value_type* last = posts[i].values + posts[i].num;

                while (itc != cands.end() || p != last) {
                    if (itc == cands.end() || (p != last && itc->value > *p)) {
                        tmp.push_back(candidate_type(*p, 1));
                        ++p;
                    } else if (p == last || (itc != cands.end() && itc->value < *p)) {
                        tmp.push_back(candidate_type(itc->value, itc->num));
                        ++itc;
                    } else {
                        tmp.push_back(candidate_type(itc->value, itc->num+1));
                        ++itc;
                        ++p;
                    }
                }
                std::swap(cands, tmp);
            }

            // No initial candidate is found.
            if (cands.empty()) {
                continue;
            }

            // Step 2: count the number of matches with remaining queries.
            for (;i < qsize;++i) {
                candidates_type tmp;
                typename candidates_type::const_iterator itc;
                const value_type* first = posts[i].values;
                const value_type* last = posts[i].values + posts[i].num;

                // For each active candidate.
                for (itc = cands.begin();itc != cands.end();++itc) {
                    int num = itc->num;
                    if (std::binary_search(first, last, itc->value)) {
                        ++num;
                    }

                    if (mmin <= num) {
                        // This candidate has sufficient matches.
                        results.push_back(itc->value);
                    } else if (num + (qsize - i - 1) >= mmin) {
                        // This candidate still has the chance.
                        tmp.push_back(candidate_type(itc->value, num));
                    }
                }
                std::swap(cands, tmp);

                // Exit the loop if all candidates are pruned.
                if (cands.empty()) {
                    break;
                }
            }

            if (!cands.empty()) {
                // Step 2 was not performed.
                typename candidates_type::const_iterator itc;
                for (itc = cands.begin();itc != cands.end();++itc) {
                    if (mmin <= itc->num) {
                        results.push_back(itc->value);
                    }
                }
            }
        }
    }

protected:
    /**
     * Open the index storing strings of the specific size.
     *  @param  base            The base name of the indices.
     *  @param  size            The size of strings.
     *  @return hashtbl_type&   The hash table of the index.
     */
    hashtbl_type& open_index(const std::string& base, int size)
    {
        index_type& index = m_indices[size-1];
        if (!index.table.is_open()) {
            std::stringstream ss;
            ss << base << '.' << size << ".cdb";
            index.image.open(ss.str().c_str(), std::ios::in);
            if (index.image.is_open()) {
                index.table.open(index.image.data(), index.image.size());
            }
        }

        return index.table;
    }
};



/**
 * A SimString database reader.
 *  This template class retrieves string from a SimString database.
 *
 *  Inheriting the base class ngramdb_reader_base that retrieves string IDs
 *  from a query feature set, this class manages the master string table,
 *  which maintains associations between strings and string IDs.
 */
class reader
    : public ngramdb_reader_base<uint32_t>
{
public:
    /// The type of an n-gram generator.
    typedef ngram_generator ngram_generator_type;
    /// The type of the base class.
    typedef ngramdb_reader_base<uint32_t> base_type;

protected:
    int m_ngram_unit;
    bool m_be;
    int m_char_size;

    /// The content of the master file.
    std::vector<char> m_strings;

public:
    /**
     * Constructs an object.
     */
    reader()
    {
    }

    /**
     * Destructs an object.
     */
    virtual ~reader()
    {
        close();
    }

    /**
     * Opens a SimString database.
     *  @param  name        The name of the SimString database.
     *  @return bool        \c true if the database is successfully opened,
     *                      \c false otherwise.
     */
    bool open(const std::string& name)
    {
        uint32_t num_entries, max_size;

        // Open the master file.
        std::ifstream ifs(name.c_str(), std::ios_base::in | std::ios_base::binary);
        if (ifs.fail()) {
            this->m_error << "Failed to open the master file: " << name;
            return false;
        }

        // Obtain the size of the master file.
        ifs.seekg(0, std::ios_base::end);
        size_t size = (size_t)ifs.tellg();
        ifs.seekg(0, std::ios_base::beg);

        // Read the image of the master file.
        m_strings.resize(size);
        ifs.read(&m_strings[0], size);
        ifs.close();

        // Check the file header.
        const char* p = &m_strings[0];
        if (size < 36 || std::strncmp(p, "SSDB", 4) != 0) {
            this->m_error << "Incorrect file format";
            return false;
        }
        p += 4;

        // Check the byte order.
        if (BYTEORDER_CHECK != read_uint32(p)) {
            this->m_error << "Incompatible byte order";
            return false;
        }
        p += 4;

        // Check the version.
        if (SIMSTRING_STREAM_VERSION != read_uint32(p)) {
            this->m_error << "Incompatible stream version";
            return false;
        }
        p += 4;

        // Check the chunk size.
        if (size != read_uint32(p)) {
            this->m_error << "Inconsistent chunk size";
            return false;
        }
        p += 4;

        // Read the unit of n-grams, begin/end flag.
        m_char_size = (int)read_uint32(p);
        p += 4;
        m_ngram_unit = (int)read_uint32(p);
        p += 4;
        m_be = (read_uint32(p) != 0);
        p += 4;

        // Read the number of enties.
        num_entries = read_uint32(p);
        p += 4;

        // Read the maximum size of strings in the database.
        max_size = read_uint32(p);

        base_type::open(name, (int)max_size);
        return true;
    }

    /**
     * Closes the database.
     */
    void close()
    {
        base_type::close();
    }

    int char_size() const
    {
        return m_char_size;
    }

    /**
     * Retrieves strings that are similar to the query.
     *  @param  query           The query string.
     *  @param  measure         The similarity measure.
     *  @param  alpha           The threshold for approximate string matching.
     *  @param  ins             The insert iterator that receives retrieved
     *                          strings.
     *  @see    ::simstring::exact, ::simstring::dice, ::simstring::cosine,
     *          ::simstring::jaccard, ::simstring::overlap
     */
    template <class string_type, class insert_iterator>
    void retrieve(
        const string_type& query,
        int measure,
        double alpha,
        insert_iterator ins
        )
    {
        switch (measure) {
        case exact:
            this->retrieve<simstring::measure::exact>(query, alpha, ins);
            break;
        case dice:
            this->retrieve<simstring::measure::dice>(query, alpha, ins);
            break;
        case cosine:
            this->retrieve<simstring::measure::cosine>(query, alpha, ins);
            break;
        case jaccard:
            this->retrieve<simstring::measure::jaccard>(query, alpha, ins);
            break;
        case overlap:
            this->retrieve<simstring::measure::overlap>(query, alpha, ins);
            break;
        }
    }

    /**
     * Retrieves strings that are similar to the query.
     *  @param  measure_type    The similarity measure.
     *  @param  query           The query string.
     *  @param  alpha           The threshold for approximate string matching.
     *  @param  ins             The insert iterator that receives retrieved
     *                          strings.
     *  @see    ::simstring::measure::exact, ::simstring::measure::dice,
     *          ::simstring::measure::cosine, ::simstring::measure::jaccard,
     *          ::simstring::measure::overlap
     */
    template <class measure_type, class string_type, class insert_iterator>
    void retrieve(
        const string_type& query,
        double alpha,
        insert_iterator ins
        )
    {
        typedef std::vector<string_type> ngrams_type;
        typedef typename string_type::value_type char_type;

        ngram_generator_type gen(m_ngram_unit, m_be);
        ngrams_type ngrams;
        gen(query, std::back_inserter(ngrams));

        typename base_type::results_type results;
        base_type::overlapjoin<measure_type>(ngrams, alpha, results);

        typename base_type::results_type::const_iterator it;
        const char* strings = &m_strings[0];
        for (it = results.begin();it != results.end();++it) {
            const char_type* xstr = reinterpret_cast<const char_type*>(strings + *it);
            *ins = xstr;
        }
    }


protected:
    inline uint32_t read_uint32(const char* p) const
    {
        return *reinterpret_cast<const uint32_t*>(p);
    }
};


#endif

}

}

QT_END_NAMESPACE

#endif
