#include "hyper.h"

/* aloja memoria para un punto */
point_t *
point_alloc (void)
{
  point_t *r;
  r = calloc (1, sizeof (struct point));
  return r;
}

/* libera memoria de un punto */
void
point_free (point_t * x)
{
  free (x);
}

/* aloja memoria para n puntos */
point_t **
pointset_alloc (unsigned int n)
{
  int i;
  point_t **r = calloc (n, sizeof (point_t *));
  for (i = 0; i < n; i++)
    r[i] = point_alloc ();
  return r;
}

/* libera memoria de n puntos en x */
void
pointset_free (point_t ** x, unsigned int n)
{
  int i;
  for (i = 0; i < n; i++)
    point_free (x[i]);
  free (x);
}

/* asigna las coordenadas x,y a un punto q */
void
point_assign (point_t * q, unsigned int x, unsigned int y)
{
  q->x = mod_p (x);
  q->y = mod_p (y);
}
