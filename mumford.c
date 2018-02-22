#include "hyper.h"

void
divisorform_g2_print(mumford_t *D) {
int x[2],y[2],r;

if(poly_iszero(D->u)) {
        printf("(0,0)-2inf\n");
        return;
}
                

r=poly_getroots(D->u,&x);
if(r == 2) {
y[0]=poly_eval(D->v,x[0]);
y[1]=poly_eval(D->v,x[1]);
printf("(%d,%d)+(%d,%d)-2inf\n",x[0],y[0],x[1],y[1]); 
}
if(r == 1) {
y[0]=poly_eval(D->v,x[0]);
printf("(%d,%d)+(I,J)-2inf\n",x[0],y[0]); 
}
if(r == 0)
printf("(I,J)+(K,L)-2inf\n");
return;
}


/* Aloja un divisor de mumford de genero g */
mumford_t *
mumford_alloc (int g)
{
  mumford_t *r;
  r = calloc (1, sizeof (struct mumford));
  r->u = poly_alloc (g);
  r->v = poly_alloc (g - 1);
  return r;
}


mumford_t *
mumford_build_g2 (point_t * a, point_t * b)
{
  mumford_t *r;
  int m;
  poly_t *t1, *t2;
  t1 = poly_alloc (1);
  t2 = poly_alloc (1);
  t1->c[0] = mod_p (-1 * a->x);
  t2->c[0] = mod_p (-1 * b->x);
  t1->c[1] = t2->c[1] = 1;
  r = mumford_alloc (2);
  poly_mul (t1, t2, r->u);
//  printf ("t1=");
//  poly_print (t1);
//  printf ("t2=");
//  poly_print (t2);
//  printf ("t1*t2=");
//  poly_print (r->u);

  m = mod_p ((a->y - b->y) * mod_inv (mod_p (a->x - b->x), p));
//  printf ("%d-%d = %d, M=%d\n", a->y, b->y, a->y - b->y, m);
  r->v->c[0] = mod_p ((-1 * m * a->x) + a->y);
  r->v->c[1] = m;
  r->v->t = 1;
  poly_free (t1);
  poly_free (t2);
  return r;
}

int 
mumford_equal(mumford_t *a, mumford_t *b) {
if(!poly_equal(a->u,b->u)) 
        return 0;
if(!poly_equal(a->v,b->v)) 
        return 0;
return 1;
}

void
mumford_free (mumford_t * d)
{
  poly_free (d->u);
  poly_free (d->v);
  free (d);
}
