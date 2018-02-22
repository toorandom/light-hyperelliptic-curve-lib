
/* Campo Zp */


#ifndef _HYPER_H
#define _HYPER_H 1
extern int p; 

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>

typedef struct sq_matrix
{
  int **c;
  int n;
} sq_matrix_t;

/* Polinomio en K[x] con t terminos y coeficientes en K almacenados en *c
 * el elemento mem es para uso interno al momento de limpiar y liberar 
 */
typedef struct poly
{
  unsigned int t;
  int *c; 
  unsigned int mem;
} poly_t;


/* point (x,y) en KxK */
typedef struct point
{
  int x, y;
} point_t;

/* divisores en forma de mumford */

typedef struct mumford {
poly_t *u;
poly_t *v;
} mumford_t;

/* prototipos */
poly_t *
poly_alloc (int);

mumford_t *
mumford_alloc (int);

point_t *
point_alloc (void);

mumford_t *
mumford_build_g2 (point_t *, point_t *);

#endif /* _HYPER_H */
