#include "hyper.h"

int p = 31;
/* calcula x mod p aunque x sea negativo porque la instruccion mod del procesador no sabe poner en terminos de Zp */
unsigned int 
mod_p (int x)
{
  int f;
  if (x < 0)
    {   
      f = ((x / p) * -1) + 1;
      return ((f * p) + x)%p;
    }   
  return x % p;

}


/*
* Esta funcion calcula la combinacion lineal de x,y que da como resultado el gcd(x,y) 
* o sea:
*
* x*r[0] + y*r[1] = gcd(x,y)
*
* recordemos que r[1] nos sirve como el inverso multiplicativo modulo x si x,y son primos relativos
* esta funcion es usada por la funcion mod_inv();
*/

void
egcd (int x, int y, int *rv)
{
  int t, respaldoy = y;
  rv[0] = x;
  rv[1] = 1;
  rv[2] = 0;

  int c, d;
  c = 0;
  d = 1;

  int q;
  int t1, t2;

  while (y != 0)
    {
      q = rv[0] / y;
      t1 = y;
      y = rv[0] - y * q;
      rv[0] = t1;

      t1 = c;
      t2 = d;

      c = rv[1] - c * q;
      d = rv[2] - d * q;
      rv[1] = t1;
      rv[2] = t2;

      if (rv[1] < 0)
        {
          t = (respaldoy + rv[1]) % respaldoy;
          rv[1] = t;
        }
    }
}
/*
 * Calcula el inverso de x modulo y
 * usando el algoritmo de euclides extendido
 */
int
mod_inv (int x, int y)
{
  int r[3];
  egcd (mod_p (x), y, (int *) &r);
  r[1] = mod_p (r[1]);
  return r[1];
}


/*
 * Calcula el maximo comun divisor de x,y
 */
int
gcd (int x, int y)
{
  int g;
  if (x < 0)
    x = -x;
  if (y < 0)
    y = -y;
  if (x + y == 0)
    exit (-1);
  g = y;
  while (x > 0)
    {
      g = x;
      x = y % x;
      y = g;
    }
  return g;
}


/* funciones que regresan maximo y minimo de 2 elementos */
int
max (int x, int y)
{
  return (x >= y) ? x : y;
}

int
min (int x, int y)
{
  return (x <= y) ? x : y;
}
/* exponenciacion modulo p muy rapida descomponiendo en bits el exponente e incrementando la base cuando el bit es 1 */
int
mod_exp (int b, int e)
{
  int r = 1;
  while (e > 0)
    {
      if (e & 1)
        r = (r * b) % p;
      e = e >> 1;
      b = (b * b) % p;
    }
  return r;
}

