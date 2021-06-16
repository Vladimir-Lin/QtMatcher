#include <qtmatcher.h>

unsigned int N::Pattern::TableBegin(void)
{
  return ( 16 + sizeof(N::Pattern::ndbTableRef) * N::Pattern::NDB_NUM_TABLES ) ;
}
