#include <qtmatcher.h>

int N::Pattern::Overlap::min_size(int qsize,double alpha)
{ Q_UNUSED ( qsize ) ;
  Q_UNUSED ( alpha ) ;
  return 1 ;
}

int N::Pattern::Overlap::max_size(int qsize,double alpha)
{ Q_UNUSED ( qsize ) ;
  Q_UNUSED ( alpha ) ;
  return (int)INT_MAX ;
}

int N::Pattern::Overlap::min_match(int qsize,int rsize,double alpha)
{
  return (int)std::ceil ( alpha * min ( qsize , rsize ) ) ;
}
