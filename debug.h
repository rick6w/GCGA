//#ifndef DEBUG_H
//#define DEBUG_H


#include <stdio.h>
#include <stdlib.h>
#include "types.h"

#ifndef DEBUG_H
#define DEBUG_H

void printEpitope(epitope* epitope)
{
  for(int i=0; i<8; i++) printf("%d ", epitope->acids[i]);
  printf("\n");
}

void printEpitopes(epitope* epitopes, int numEpitopes)
{
  printf("Epitopes:\n");
  for(int i=0; i<numEpitopes; i++)
  {
    printf("#%d: ", i);
    printEpitope(&(epitopes[i]));
  }
  printf("-----\n");
}

void printAntibodies_2rows(antibody* antibodies)
{
  int i, j, n=NUM_ANTIBODIES/2, m=NUM_ANTIBODIES;
  antibody* a;
  printf("ANTIBODIES:\n----------\n");
  for(i=0; i<n; i++) printf("Antibody #%2d                  ", i);
  printf("\n");
  for(i=0; i<n; i++) printf("Score: %7d                ", antibodies[i].score);
  printf("\n");
  for(i=0; i<n; i++) printf("Amino acids:                  ");
  printf("\n");
  for(i=0; i<n; i++) {
    a = &(antibodies[i]);
    for(j=0; j<8; j++) printf("%2d ", a->acids[j]);
    printf("      ");
  }
  printf("\n");
  for(i=0; i<n; i++) printf("Previous affinities:          ");
  printf("\n");
  for(i=0; i<n; i++) {
    a = &(antibodies[i]);
    for(j=0; j<8; j++) printf("%2X ", a->affinity_prev[j]);
    printf("      ");
  }
  printf("\n");
  for(i=0; i<n; i++) printf("Current affinities:           ");
  printf("\n");
  for(i=0; i<n; i++) {
    a = &(antibodies[i]);
    for(j=0; j<8; j++) printf("%2X ", a->affinity[j]);
    printf("      ");
  }
  printf("\n");
  for(i=0; i<n; i++) printf("Delta affinity values:        ");
  printf("\n");
  for(i=0; i<n; i++) {
    a = &(antibodies[i]);
    for(j=0; j<8; j++) printf("%2d ", (a->delta_affinity[j])+1);
    printf("      ");
  }
  printf("\n\n");

  for(i=n; i<m; i++) printf("Antibody #%2d                  ", i);
  printf("\n");
  for(i=n; i<m; i++) printf("Score: %7d                ", antibodies[i].score);
  printf("\n");
  for(i=n; i<m; i++) printf("Amino acids:                  ");
  printf("\n");
  for(i=n; i<m; i++) {
    a = &(antibodies[i]);
    for(j=0; j<8; j++) printf("%2d ", a->acids[j]);
    printf("      ");
  }
  printf("\n");
  for(i=n; i<m; i++) printf("Previous affinities:          ");
  printf("\n");
  for(i=n; i<m; i++) {
    a = &(antibodies[i]);
    for(j=0; j<8; j++) printf("%2X ", a->affinity_prev[j]);
    printf("      ");
  }
  printf("\n");
  for(i=n; i<m; i++) printf("Current affinities:           ");
  printf("\n");
  for(i=n; i<m; i++) {
    a = &(antibodies[i]);
    for(j=0; j<8; j++) printf("%2X ", a->affinity[j]);
    printf("      ");
  }
  printf("\n");
  for(i=n; i<m; i++) printf("Delta affinity values:        ");
  printf("\n");
  for(i=n; i<m; i++) {
    a = &(antibodies[i]);
    for(j=0; j<8; j++) printf("%2d ", (a->delta_affinity[j])+1);
    printf("      ");
  }
  printf("\n\n");
}

void printAntibody(antibody* antibody)
{
  printf("Score: %d\n", antibody->score);
  printf("Amino acids\n");
  for(int i=0; i<8; i++) printf("%d ", antibody->acids[i]);
  printf("\nPrevious affinities:\n");
  for(int i=0; i<8; i++) printf("%X ", antibody->affinity_prev[i]);
  printf("\nCurrent affinities:\n");
  for(int i=0; i<8; i++) printf("%X ", antibody->affinity[i]);
  printf("\nDelta affinity values:\n");
  for(int i=0; i<8; i++) printf("%d ", antibody->delta_affinity[i]);
  printf("\n\n");
}

void printAntibodies(antibody* antibodies, int numAntibodies)
{
  printf("Antibodies:\n");
  for(int i=0; i<numAntibodies; i++)
  {
    printf("#%d:\n", i);
    printAntibody(&(antibodies[i]));
  }
  printf("-----\n");
}

void printRule(rule* rule)
{
  printf("[");
  for(int i=0; i<4; i++) printf("%d ", (rule->ops[i])+1);
  printf("]");
}

void printGenome(Bcell* bcell)
{
  for(int i=0; i<8; i++) printRule(&(bcell->genome[i]));
  printf("\n");
}

void printBcell(Bcell* bcell, int numAntibodies)
{
  printf("Fitness: %d\n", bcell->fitness);
  printf("Rule Genome:\n");
  printGenome(bcell);
  printAntibodies_2rows(bcell->antibodies);
  printf("-----\n");
}

void printBcells(Bcell* bcells, int numBcells, int numAntibodies)
{
  printf("B CELLS:\n########\n");
  for(int i=0; i<numBcells; i++)
  {
    printf("B Cell #%d\n", i);
    printBcell(&(bcells[i]), numAntibodies);
  }
  printf("########\n\n");
}

void printStats(Bcell* bcells, int numBcells)
{
  Bcell* bcell;
  for(int i=0; i<numBcells; i++)
  {
    bcell = &(bcells[i]);
    printf("Rule Genome: ");
    printGenome(bcell);
    printf("Fitness: %d\n\n", bcell->fitness);
  }
}

void printAverage(Bcell* bcells, int numBcells)
{
  double avg=0;
  Bcell* bcell;
  for(int i=0; i<numBcells; i++) avg += bcells[i].fitness;
  avg /= numBcells;
  printf("%f ", avg);
}

//comparison function for qsort()
int compare(const void* a, const void* b)
{
  Bcell* a1 = (Bcell*)a;
  Bcell* b1 = (Bcell*)b;
  return (a1->fitness < b1->fitness)-(a1->fitness > b1->fitness);
}

void printFittestRule(Bcell* bcells, int numBcells)
{
  qsort(bcells, numBcells, sizeof(Bcell), compare);
  printGenome(&(bcells[0]));
}

void printRules(rule* rules)
{
  for(int i=0; i<256; i++)
  {
    printf("Rule #%d\n", i);
    for(int j=0; j<4; j++) printf("%d: %d\n", j, rules[i].ops[j]);
    printf("-----\n\n");
  }
}

#endif
