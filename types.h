#ifndef TYPE_H
#define TYPE_H

int NUM_BCELLS = 100;
int NUM_ANTIBODIES = 10;
int NUM_FITTEST = 50;
int NUM_TRIALS = 100;
int NUM_GENERATIONS = 100;
char ACIDS[20] = {1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20};

typedef struct rule
{
  int ops[4];
} rule;

typedef struct epitope
{
  char acids[8];
} epitope;

typedef struct antibody
{
  int score;
  char acids[8];
  char affinity[8];
  char affinity_prev[8];
  char delta_affinity[8];
} antibody;

typedef struct Bcell
{
  int fitness;
  rule genome[8];
  antibody* antibodies;
} Bcell;

#endif
