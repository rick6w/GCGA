#include "debug.h"
#include "types.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

rule* initRuleset()
{
  int n=0;
  rule* ruleset = (rule*)malloc(256 * sizeof(rule));
  rule* rule = &(ruleset[n]);
  for(int i=0; i<4; i++) {
    for(int j=0; j<4; j++) {
      for(int k=0; k<4; k++) {
        for(int l=0; l<4; l++) {
          rule->ops[0] = i;
          rule->ops[1] = j;
          rule->ops[2] = k;
          rule->ops[3] = l;
          rule = &(ruleset[n++]);
        }
      }
    }
  }
  return ruleset;
}

epitope* initEpitopes(int n)
{
  int r;
  epitope* epitopes = (epitope*) malloc(n * sizeof(epitope));
  for(int i=0; i<n; i++)
  {
    for(int j=0; j<8; j++)
    {
      r=rand()%20;
      epitopes[i].acids[j] = ACIDS[r];
    }
  }
  return epitopes;
}

antibody* initAntibodies(int n)
{
  int r;
  antibody* antibodies = (antibody*) malloc(n * sizeof(antibody));
  for(int i=0; i<n; i++)
  {
    for(int j=0; j<8; j++)
    {
      r=rand()%20;
      antibodies[i].acids[j] = ACIDS[r];
      antibodies[i].score = 0;
    }
  }
  return antibodies;
}

Bcell* initBcells(int numBcells, int numAntibodies, rule* ruleset)
{
  Bcell* bcells = (Bcell*) malloc(numBcells * sizeof(Bcell));
  for(int i=0; i<numBcells; i++)
  {
    for(int j=0; j<8; j++) bcells[i].genome[j] = ruleset[rand()%256];
    bcells[i].antibodies = initAntibodies(numAntibodies);
    bcells[i].fitness = 0;
  }
  return bcells;
}

void applyOp(char op, char* acid)
{
  switch(op)
  {
    case 0:
      if(*acid<20) (*acid)++;
      else (*acid)=1;
      break;
    case 1:
      if(*acid>1) (*acid)--;
      else (*acid)=20;
      break;
    case 2:
      break;
    case 3:
      *acid = ACIDS[rand()%20];
      break;
    default:
      break;
  }
}

void modifyAntibody(antibody* antibody, rule* genome)
{
  char delta_affinity;
  for(int i=0; i<8; i++)
  {
    delta_affinity = antibody->delta_affinity[i];
    applyOp(genome[i].ops[delta_affinity], &(antibody->acids[i]));
  }
}

void computeDeltaAffinity(antibody* antibody)
{
  char affinity, affinity_prev, max=0x1F;
  for(int i=0; i<8; i++)
  {
    affinity = antibody->affinity[i];
    affinity_prev = antibody->affinity_prev[i];
    if((affinity==max) && (affinity==affinity_prev))
      antibody->delta_affinity[i] = 2;
    else if(affinity > affinity_prev)
      antibody->delta_affinity[i] = 0;
    else if(affinity < affinity_prev)
      antibody->delta_affinity[i] = 1;
    else antibody->delta_affinity[i] = 3;
  }
}

void computeAffinities(Bcell* bcell, epitope* epitopes, int n)
{
  epitope* epitope;
  antibody* antibody;
  for(int i=0; i<n; i++)
  {
    antibody = &(bcell->antibodies[i]);
    memcpy(antibody->affinity_prev,
           antibody->affinity,
           sizeof(antibody->affinity));
    for(int j=0; j<n; j++)
    {
      epitope = &(epitopes[j]);
      for(int k=0; k<8; k++)
      {
        antibody->affinity[k] = ((antibody->acids[k])
                              ^ ~(epitope->acids[k]))
                              & 0x1f; //zero 3 msbs
      }
    }
    computeDeltaAffinity(antibody);
  }
}

void computeFitness(Bcell* bcells, int numBcells, int numAntibodies)
{
  int sum;
  antibody* antibody;
  for(int i=0; i<numBcells; i++)
  {
    sum = 0;
    for(int j=0; j<numAntibodies; j++)
    {
      antibody = &(bcells[i].antibodies[j]);
      for(int k=0; k<8; k++) antibody->score += (antibody->affinity[k]);
      sum += (antibody->score);
    }
    bcells[i].fitness = sum;
  }
}

//comparison function for qsort()
int cmp(const void* a, const void* b)
{
  Bcell* a1 = (Bcell*)a;
  Bcell* b1 = (Bcell*)b;
  return (a1->fitness < b1->fitness)-(a1->fitness > b1->fitness);
}

Bcell* selectFittest(Bcell* bcells, int numBcells, int numFittest)
{
  Bcell* fittest = (Bcell*) malloc(numFittest * sizeof(Bcell));
  qsort(bcells, numBcells, sizeof(Bcell), cmp);
  for(int i=0; i<numFittest; i++) fittest[i] = bcells[i];
  return fittest;
}

void mutate(){}

void breedGenome(Bcell* parent1, Bcell* parent2, Bcell* child)
{
  int split = rand()%8;
  for(int i=0; i<split; i++) child->genome[i] = parent1->genome[i];
  for(int i=split; i<8; i++) child->genome[i] = parent2->genome[i];
}

Bcell* breedBcells(Bcell* bcells,
                   int numBcells,
                   int numNewBcells,
                   int numAntibodies)
{
  Bcell *parent1, *parent2;
  Bcell* newBcells = (Bcell*) malloc(numNewBcells * sizeof(Bcell));
  for(int i=0; i<numNewBcells; i++)
  {
    parent1 = &(bcells[rand()%numBcells]);
    parent2 = &(bcells[rand()%numBcells]);
    breedGenome(parent1, parent2, &(newBcells[i]));
    newBcells[i].antibodies = initAntibodies(numAntibodies);
    newBcells[i].fitness = 0;
    //printf("*********%d fitness=%d\n", i, newBcells[i].fitness);
  }
  return newBcells;
}

void evolve(Bcell* bcells,
            int numBcells,
            int numAntibodies,
            int numFittest,
            int numTrials,
            int numGenerations)
{
  epitope* epitopes;
  for(int i=0; i<numGenerations; i++)
  {
    epitopes = initEpitopes(numAntibodies);
    for(int trial=0; trial<numTrials; trial++)
    {
      for(int j=0; j<numBcells; j++)
      {
        computeAffinities(&(bcells[j]), epitopes, numAntibodies);
        for(int k=0; k<numAntibodies; k++)
        {
          modifyAntibody(&(bcells[j].antibodies[k]), bcells[j].genome);
        }
      }
    }
    computeFitness(bcells, numBcells, numAntibodies);
    //printf("Stats after gen%d:\n", i);
    //printStats(bcells, numBcells);
    //printBcells(bcells, numBcells, numAntibodies);
    Bcell* fittest = selectFittest(bcells, numBcells, numFittest);
    free(bcells);
    free(epitopes);
    bcells = breedBcells(fittest, numFittest, numBcells, numAntibodies);
    free(fittest);
  }
  printFittestRule(bcells, numBcells);
  free(bcells);
}

void run(int numBcells,
         int numAntibodies,
         int numFittest,
         int numTrials,
         int numGenerations,
         rule* ruleset)
{
  Bcell* bcells = initBcells(numBcells, numAntibodies, ruleset);
  evolve(bcells,
         numBcells,
         numAntibodies,
         numFittest,
         numTrials,
         numGenerations);
}

void main()
{
  srand(12345);
  rule* RULESET = initRuleset();

  for (int i=0; i<50; i++)
  {
    run(NUM_BCELLS,
        NUM_ANTIBODIES,
        NUM_FITTEST,
        NUM_TRIALS,
        NUM_GENERATIONS,
        RULESET);
  }

  //printRules(RULESET);
}

