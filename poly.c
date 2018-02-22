/*
 * Biblioteca para manejar polinomios en formato de espacio vectorial, 
 * i.e. a_{n} x^(n) + a_{n-1} x^(n-1) + ... + a_{1} x^(1) + a_{0} x^(0) = (a_{0},a_{1}, ... , a_{n-1}, a_{n})
 * es decir el coeficiente de grado menor es el primero  
 * esta biblioteca implementa el anillo (dominio entero) <Zp[x],+,*>
 "* el point de esto es implementar el calculo del polinomio de lagrange sobre un campo finito de caracteristica p!=2
 * 3/Ago/2011
 * rduarte@ciencias.unam.mx
 */

#include "hyper.h"


/* llena de ceros los coeficientes del polinomio  */
void
poly_clean (poly_t * f)
{
  int i;
  if ((int) f->t <= 0)
    {
      f->t = 0;
      f->c[0] = 0;
      return;
    }
  for (i = 0; i <= f->mem; i++)
    f->c[i] = 0;
}

int
poly_equal(poly_t *f1, poly_t *f2) {
        int i;
        if(f1->t != f2->t) 
                return 0;
        for(i=0;i<f1->t ;i++) 
                if(f1->c[i] != f2->c[i]) 
                        return 0;
       
        return 1;

}

/* sobreescribe el polinomio b con el polinomio a */
void
poly_copy (poly_t * b, poly_t * a)
{
  int i;
  poly_clean (b);
  if ((int) a->t <= 0)
    {
      b->c[0] = a->c[0];
      return;
    }
  for (i = 0; i <= a->t; i++)
    b->c[i] = a->c[i];
  b->t = a->t;
}

/* obtiene raices de Zp (manera chafa, solo con fines demostrativos) */

int
poly_getroots(poly_t *f, int *r) {
int i,j=0;
for(i=0;i<p;i++)
        if(poly_eval(f,i) == 0)  {
                r[j] = i;
                j++;
        }
return j;
}


/* aloja t terminos para un nuevo polinomio f, este polinomio sera de grado t (puede ser 0) */
poly_t *
poly_alloc (int t)
{
  poly_t *f;
  f = calloc (1, sizeof (struct poly));
  f->c = calloc (t + 2, sizeof (int));
  f->t = t;
  f->mem = t + 1;
  return f;
}

/* libera memoria usada por los coeficientes de f */
void
poly_free (poly_t * f)
{
  free (f->c);
  free (f);
}


/* pone los coeficientes negativos del polinomio f como terminos de Zp y quita los ceros arriba del grado en caso de existir */
void
poly_normalize (poly_t * f)
{
  int i;
  if ((signed int) f->t < 0)
    {
      poly_clean (f);
      f->t = 0;
      return;
    }

  for (i = 0; i < f->t; i++)
    f->c[i] = mod_p (f->c[i]);


  i = f->t;
  if (f->c[f->t] == 0)
    while ((f->c[i] == 0) & ((int) f->t > 0))
      {
	i--;
	f->t -= 1;
      }

  if (f->t < 0)
    {
      poly_clean (f);
      f->t = 0;
      return;
    }
}


/* evalua el polinomio f en x mod p */
int
poly_eval (poly_t * f, int x)
{
  int i, r = f->c[0];
  for (i = 1; i <= f->t; i++)
    {
      r += (f->c[i] * mod_exp (x, i));
      r %= p;
    }
  return r;
}

/* suma polinomios a+b=c se divide en casos para que sea mas rapido */
void
poly_add (poly_t * a, poly_t * b, poly_t * c)
{
  int i;

  if ((a->t == 0) & (b->t == 0))
    {
      c->c[0] = a->c[0] + b->c[0];
      c->t = 0;
      return;
    }
  poly_clean (c);
  c->t = max (a->t, b->t);

  if (poly_iszero (a) & poly_iszero (b))
    {
      poly_clean (c);
      c->t = 0;
      return;
    }
  if (poly_iszero (a))
    {
      poly_copy (c, b);
      return;
    }

  if (poly_iszero (b))
    {
      poly_copy (c, a);
      return;
    }

  for (i = 0; i <= min (a->t, b->t); i++)
    {
      c->c[i] = a->c[i] + b->c[i];
      c->c[i] %= p;
    }
  for (i = min (a->t, b->t) + 1; i <= max (a->t, b->t); i++)
    {
      if (max (a->t, b->t) == a->t)
	c->c[i] = a->c[i] % p;
      else
	c->c[i] = b->c[i] % p;
    }
  poly_normalize (c);
}

/* calcula g  tal que g = kf con k en Zp */
void
poly_scalar_mul (int k, poly_t * f)
{
  int i;
  poly_normalize (f);
  if (k == 0)
    {
      poly_clean (f);
      f->t = 0;
      return;
    }

  for (i = 0; i <= f->t; i++)
    {
      f->c[i] *= k;
      f->c[i] %= p;
    }
  poly_normalize (f);
}




/* multiplica polinomios, tambien se usan varios casos para agilizarla */
void
poly_mul (poly_t * a, poly_t * b, poly_t * c)
{
  int i, j;
  poly_clean (c);
  poly_normalize (a);
  poly_normalize (b);
  if (((a->t == 0) & (a->c[0] == 0)) || ((b->t == 0) & (b->c[0] == 0)))
    {
      c->t = 0;
      return;
    }


  if ((a->t == 0) & (b->t == 0))
    {
      c->t = 0;
      c->c[0] = a->c[0] * b->c[0];
      c->c[0] %= p;
      return;
    }
  if (a->t == 0)
    {
      poly_copy (c, b);
      poly_scalar_mul (a->c[0], c);
      return;
    }
  if (b->t == 0)
    {
      poly_copy (c, a);
      poly_scalar_mul (b->c[0], c);
      return;
    }

  poly_normalize (a);
  poly_normalize (b);
  c->t = a->t + b->t;
  poly_clean (c);
  for (i = 0; i <= a->t; i++)
    for (j = 0; j <= b->t; j++)
      {
	c->c[i + j] += (a->c[i] * b->c[j]);
	c->c[i + j] %= p;
      }
  poly_normalize (c);
}


/* imprime el polinomio f */
/*
void
poly_print (poly_t * f)
{
  int i;
  poly_normalize (f);

  if ((signed int) f->t <= 0)
    {
      printf ("(%d)\n", f->c[0]);
      return;
    }



  printf ("(");
  for (i = 0; i <= f->t; i++)
    {
      if (i != f->t)
	printf ("%d,", f->c[i]);
      else
	printf ("%d", f->c[i]);
    }
  printf (")\n");
}
*/

void
poly_print (poly_t * f)
{
  int i;
  poly_normalize (f);
  if ((f->t <= 0) & (f->c[0] == 0))
    {
      printf ("(0)\n");
      return;
    }

  if (f->t == 0)
    {
      printf ("(%d)\n", mod_p (f->c[0]));
      return;
    }
  printf ("(");
  for (i = 0; i <= f->t; i++)
    {
      if (i != 0)
	printf ("+");

      printf ("%d", mod_p (f->c[i]));

      if (i > 0)
	printf ("*x");

      if (i > 1)
	printf ("^%d", i);
    }
  printf (")\n");

}

/* regresa verdadero si f es el polinomio 0 */
int
poly_iszero (poly_t * f)
{
  int i;
  for (i = 0; i < f->t; i++)
    if (f->c[i] != 0)
      return 0;
  return 1;
}

/* mueve n lugares a la derecha los coeficientes del polinomio, o multiplica por x^n como lo quieras ver */
void
poly_rshift (poly_t * a, unsigned int n)
{
  int i;
  poly_t *t;
  if (!n)
    return;
  t = poly_alloc (a->t);
  poly_copy (t, a);
  for (i = 0; i <= a->t; i++)
    a->c[i + n] = t->c[i];
  for (i = 0; i < n; i++)
    a->c[i] = 0;
  a->t += n;
  poly_free (t);
  return;
}


/* mueve n lugares a la izquierda los coeficientes del polinomio, o divide entre x^n y pone ceros cuando el grado del dividendo es mayor que el del divisor */

void
poly_lshift (poly_t * a, unsigned int n)
{
  int i;
  poly_t *t;
  if (!n)
    return;
  t = poly_alloc (a->t);
  poly_copy (t, a);
  for (i = 0; i <= a->t; i++)
    a->c[i] = t->c[i + n];
  for (i = t->t - n + 1; i <= n + t->t; i++)
    a->c[i] = 0;
  a->t -= n;
  poly_free (t);
  return;
}


/* calcula q y r tal que q=a/b y bq+r=a , nota que q y r ya deben tener memoria lista */
void
poly_div (poly_t * a, poly_t * b, poly_t * q, poly_t * r)
{
  poly_t *d, *dneg, *N, *D, *Nt;
  N = poly_alloc (a->t);
  Nt = poly_alloc (a->t );
  D = poly_alloc (b->t);
  d = poly_alloc (b->t + a->t );
  dneg = poly_alloc (b->t + a->t );
  poly_copy (N, a);
  poly_copy (D, b);
  poly_copy (d, D);
  q->t = a->t - b->t;
  if (N->t >= D->t)
    {
      poly_clean (q);
      while (N->t >= D->t)
	{
	  poly_copy (d, D);
	  poly_rshift (d, N->t - D->t);
	  q->c[N->t - D->t] = mod_p (N->c[N->t] * mod_inv (d->c[d->t], p));
	  poly_normalize (q);
	  poly_scalar_mul (q->c[N->t - D->t], d);
	  poly_copy (dneg, d);
	  poly_scalar_mul (-1, dneg);
	  poly_add (N, dneg, Nt);
	  poly_copy (N, Nt);
	  poly_normalize (N);
	  poly_clean (Nt);
	  poly_clean (dneg);
	}
      poly_copy (r, N);
    }
  else
    {
      poly_clean (q);
      poly_copy (r, N);
    }
  poly_free (dneg);
  poly_free (Nt);
  poly_free (N);
  poly_free (D);
  poly_free (d);
  return;
}


/*
 * Dados s apuntadores a puntos (**q) te regresa el polinomio f de grado s-1 que pasa por los s puntos sobre ZpxZp
 * Nota: este aloja memoria para el apuntador que te regresa, hay que hacerle poly_free(f) despues
 *
 * Nota2: Ya funciona pero falla cuando alguno de los puntos tiene coordenada y=0 (ha de ser una tonteria)
 */

poly_t *
lagrange (point_t ** q, int s)
{
  poly_t **l, *h, *temp, *L;
  int i, j, m, inv_m;
  l = calloc (s, sizeof (struct poly));
  h = poly_alloc (1);
  L = poly_alloc (s);
  poly_clean (L);
  L->t = s;
  temp = poly_alloc (s);
  poly_clean (temp);
  h->c[1] = 1;
  h->t = 1;
  for (j = 0; j < s; j++)
    {
      l[j] = poly_alloc (s);
      l[j]->t = 0;
      l[j]->c[0] = 1;

      for (i = 0; i < s; i++)
	{
	  if (i != j)
	    {
	      m = mod_p (q[j]->x - q[i]->x);
	      inv_m = mod_inv (m, p);
	      h->c[0] = mod_p (q[i]->x * -1);
	      poly_scalar_mul (inv_m, h);
	      poly_mul (h, l[j], temp);
	      h->c[1] = 1;
	      poly_copy (l[j], temp);
	      poly_clean (temp);
	    }
	}
      poly_scalar_mul (q[j]->y, l[j]);
      poly_add (L, l[j], temp);
      poly_copy (L, temp);
      poly_free (l[j]);
    }
  poly_free (h);
  poly_free (temp);
  free (l);
  return L;
}
