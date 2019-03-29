#ifndef PRIORITYQUEUE_H
#define PRIORITYQUEUE_H

#include <queue>
#include <limits>

using namespace std;

struct paire
{
    paire() : dist(std::numeric_limits<float>::max()), fib1(std::numeric_limits<unsigned int>::max()), fib2(std::numeric_limits<unsigned int>::max()) {}
    float dist;
    unsigned int fib1;
    unsigned int fib2;
};

class comparisonPaires
{
public:
    bool operator() (const paire p1, const paire p2) const
    {
        return p1.dist > p2.dist;
    }
};

class oppositeComparisonPaires
{
public:
    bool operator() (const paire p1, const paire p2) const
    {
        return p1.dist < p2.dist;
    }
};

template <class Q>
void clearQueue(Q & q) {
    q = Q();
}

#endif // PRIORITYQUEUE_H
