#include <qtmatcher.h>

int N::Pattern::Jaccard::min_size(int qsize,double alpha)
{
  return (int)std::ceil(alpha * qsize) ;
}

int N::Pattern::Jaccard::max_size(int qsize,double alpha)
{
  return (int)std::floor(qsize / alpha) ;
}

int N::Pattern::Jaccard::min_match(int qsize,int rsize,double alpha)
{
  return (int)std::ceil(alpha * (qsize + rsize) / (1 + alpha)) ;
}
