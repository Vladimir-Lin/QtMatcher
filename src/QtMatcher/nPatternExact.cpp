#include <qtmatcher.h>

int N::Pattern::Exact::min_size(int qsize, double alpha)
{ Q_UNUSED ( alpha ) ;
  return qsize ;
}

int N::Pattern::Exact::max_size(int qsize, double alpha)
{ Q_UNUSED ( alpha ) ;
  return qsize ;
}

int N::Pattern::Exact::min_match(int qsize, int rsize, double alpha)
{ Q_UNUSED ( rsize ) ;
  Q_UNUSED ( alpha ) ;
  return qsize ;
}
