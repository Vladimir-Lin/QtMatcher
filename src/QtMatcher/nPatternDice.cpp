#include <qtmatcher.h>

int N::Pattern::Dice::min_size(int qsize, double alpha)
{
  return (int)std::ceil(alpha * qsize / (2. - qsize)) ;
}

int N::Pattern::Dice::max_size(int qsize, double alpha)
{
  return (int)std::floor((2. - alpha) * qsize / alpha) ;
}

int N::Pattern::Dice::min_match(int qsize, int rsize, double alpha)
{
  return (int)std::ceil(0.5 * alpha * (qsize + rsize)) ;
}
