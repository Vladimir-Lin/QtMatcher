#include <qtmatcher.h>

int N::Pattern::Cosine::min_size(int qsize,double alpha)
{
  return (int)std::ceil(alpha * alpha * qsize) ;
}

int N::Pattern::Cosine::max_size(int qsize,double alpha)
{
  return (int)std::floor(qsize / (alpha * alpha)) ;
}

int N::Pattern::Cosine::min_match(int qsize,int rsize,double alpha)
{
  return (int)std::ceil(alpha * std::sqrt((double)qsize * rsize)) ;
}
