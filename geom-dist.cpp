// geometric_distribution
#include <iostream>
#include <random>

int main()
{
  const int nrolls = 10000; // number of experiments
  const int nstars = 100;   // maximum number of stars to distribute
  const int n = 10;
  const double pr = 0.15;

  std::default_random_engine generator;
  std::geometric_distribution<int> distribution(pr);

  int p[n] = {};

  for (int i = 0; i < nrolls; ++i)
  {
    int number = distribution(generator);
    if (number < n)
      ++p[number];
  }

  std::cout << "geometric_distribution (" << pr << "):" << std::endl;
  for (int i = 0; i < n; ++i)
    std::cout << i << ": " << std::string(p[i] * nstars / nrolls, '*') << std::endl;

  return 0;
}