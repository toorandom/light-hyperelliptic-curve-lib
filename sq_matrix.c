#include "hyper.h"


/* Cambia la k-esima columna de 'a' (k=1...n) por la columna c */
void
sq_matrix_change_column (sq_matrix_t * a, int *c, int k)
{
  int i;
  for (i = 0; i < a->n; i++)
    a->c[i][k] = c[i];
}

void
sq_matrix_print (sq_matrix_t * a, int *c)
{

  int i, j;
  for (i = 0; i < a->n; i++)
    {
      if (i > 0)
	printf ("\n");
      for (j = 0; j < a->n; j++)
	{printf ("%d\t", a->c[i][j]);}
      if(c != NULL)
      printf("|%d",c[i]);
	//printf ("%d\t", mod_p (a->c[i][j]));
   }
  printf("\n\n");
}

sq_matrix_t *
sq_matrix_alloc (unsigned int n)
{
  sq_matrix_t *r;
  int i;
  r = calloc (1, sizeof (struct sq_matrix));
  r->n = n;
  r->c = calloc ((n*n)+1, sizeof (int));
  for (i = 0; i < n; i++)
    r->c[i] = calloc (n, sizeof (int));
  return r;
}

void
sq_matrix_free(sq_matrix_t *m) {
 int i;
 for(i=0;i<m->n;i++)
    free(m->c[i]);
  m->n=0;
  free(m->c);
}

int
sq_matrix_determinant2 (sq_matrix_t * a)
{

  int i, j, i_c, j_c, c;
  sq_matrix_t *menor;
  int det = 0;
  menor = sq_matrix_alloc (a->n -1);

  if (a->n == 1)
    return a->c[0][0];
  if (a->n == 2)
    return (a->c[0][0] * a->c[1][1] - a->c[0][1] * a->c[1][0]);

  for (c = 0; c < a->n; c++)
    {
      i_c = 0;
      for (i = 1; i < a->n; i++)
	{
	  j_c = 0;
	  for (j = 0; j < a->n; j++)
	    {
	      if (j == c)
		continue;
	      menor->c[i_c][j_c] = a->c[i][j];
	      j_c++;
	    }
	  i_c++;
	}
      det +=
	((c % 2 == 0) ? 1 : -1) * a->c[0][c] * sq_matrix_determinant (menor);
     
    }
  return mod_p (det);
  return det;

}





void
sq_matrix_copy (sq_matrix_t * a, sq_matrix_t * b)
{
  int i, j;
  if (a->n != b->n)
    {
      fprintf (stderr, "Cannot copy matrices of different dimension\n");
      exit (EXIT_FAILURE);
    }
  for (i = 0; i < b->n; i++)
    for (j = 0; j < a->n; j++)
      a->c[i][j] = b->c[i][j];

  a->n = b->n;
}





int sq_matrix_determinant(sq_matrix_t *a)
{
   int deta11dp,deta11do,deta11,deta12dp,deta12do,deta12,deta13dp,deta13do,deta13,deta14dp,deta14do,deta14,det;
   /*
  printf("Determinante print m\n");
   sq_matrix_print(a,NULL);
  printf("end determinante print m\n");
  */

   deta11dp=(a->c[1][1]*a->c[2][2]*a->c[3][3])+(a->c[1][2]*a->c[2][3]*a->c[3][1])+(a->c[1][3]*a->c[2][1]*a->c[3][2]);

   deta11do=(a->c[1][2]*a->c[2][1]*a->c[3][3])+(a->c[1][1]*a->c[2][3]*a->c[3][2])+(a->c[1][3]*a->c[2][2]*a->c[3][1]);

   deta11=deta11dp-deta11do;


   deta12dp=(a->c[1][0]*a->c[2][2]*a->c[3][3])+(a->c[1][2]*a->c[2][3]*a->c[3][0])+(a->c[1][3]*a->c[2][0]*a->c[3][2]);

   deta12do=(a->c[1][2]*a->c[2][0]*a->c[3][3])+(a->c[1][0]*a->c[2][3]*a->c[3][2])+(a->c[1][3]*a->c[2][2]*a->c[3][0]);

   deta12=deta12dp-deta12do;


   deta13dp=(a->c[1][0]*a->c[2][1]*a->c[3][3])+(a->c[1][1]*a->c[2][3]*a->c[3][0])+(a->c[1][3]*a->c[2][0]*a->c[3][1]);

   deta13do=(a->c[1][1]*a->c[2][0]*a->c[3][3])+(a->c[1][0]*a->c[2][3]*a->c[3][1])+(a->c[1][3]*a->c[2][1]*a->c[3][0]);

   deta13=deta13dp-deta13do;


   deta14dp=(a->c[1][0]*a->c[2][1]*a->c[3][2])+(a->c[1][1]*a->c[2][2]*a->c[3][0])+(a->c[1][2]*a->c[2][0]*a->c[3][1]);

   deta14do=(a->c[1][1]*a->c[2][0]*a->c[3][2])+(a->c[1][0]*a->c[2][2]*a->c[3][1])+(a->c[1][2]*a->c[2][1]*a->c[3][0]);

   deta14=deta14dp-deta14do;

   

   det=(a->c[0][0]*(-1*-1)*deta11)+(a->c[0][1]*(-1*-1*-1)*deta12)+(a->c[0][2]*(-1*-1*-1*-1)*deta13)+(a->c[0][3]*(-1*-1*-1*-1*-1)*deta14);
   
   return det;
}


void mumford_monic(mumford_t *D) {
        int c = mod_inv(D->u->c[2],p);
        if(c != 1){
        D->u->c[0] = mod_p(D->u->c[0]*c);
        D->u->c[1] = mod_p(D->u->c[1]*c);
        D->u->c[2] = mod_p(D->u->c[2]*c);
        }
        return;
}




void
system_solve (sq_matrix_t * a, int *coef, int *res)
{
  int *det_sys = calloc (a->n, sizeof (int));
  sq_matrix_t *t;
  int det,inv_det, i;
  /*
printf("### SISTEMA DE EQ ###\n");
  sq_matrix_print(a,coef);
printf("### END SISTEMA DE EQ ###\n");
*/

//  printf("constantes=%d,%d,%d,%d\n",coef[0],coef[1],coef[2],coef[3]);
  det = sq_matrix_determinant(a);
//  sq_matrix_print(a);
  inv_det = mod_inv(det,p);
//  printf("Determinante inverso = %d\n",inv_det);
  if (det == 0)
    {
      fprintf (stderr, "The system has no solutions\n");
      sq_matrix_print(a,NULL);
      exit(1);
    }

  t = sq_matrix_alloc (a->n);

  sq_matrix_copy (t, a);
//  printf ("Matriz original\n");
//  sq_matrix_print (t);
//  printf ("su determinante es %d\n", det);
  for (i = 0; i < a->n; i++)
    {
      sq_matrix_change_column (t, coef, i);
      det_sys[i] = sq_matrix_determinant (t);
//      sq_matrix_print (t);
//      printf ("su determinante es %d\n", det_sys[i]);
      sq_matrix_copy (t, a);
      res[i] = mod_p(det_sys[i] * inv_det);
    }
  /*
  printf("Soluciones:\n");
  for(i=0;i<a->n;i++)
          printf("%d ",res[i]);
          */
  //printf("\n");
return;
}
void mumford_print(mumford_t *D) {

printf("u=");
poly_print(D->u);
printf("v=");
poly_print(D->v);
return;

}
mumford_t *hyper_g2_double_sup_2(mumford_t *D) {
sq_matrix_t *m;
mumford_t *DplusD = mumford_alloc(2);
int a,b,c,d,coefs[4],sols[4];

/* Curva hipereliptica, parte derecha de  y^2 = x^5 - 1 */
poly_t *hyperC = poly_alloc(5);
/* aqui almacenaremos el polinomio de interpolacion de grado 3 */
poly_t *L = poly_alloc(3);
/* aqui almacenaremos u*u */
poly_t *u2 = poly_alloc((2*D->u->t)+1);
/* aqui almacenaremos L^2 = L*L - C donde C es la curva hipereliptica en este caso deg(L)=3 por eso almacenamos 6 y deg(C)=5 */
poly_t *L2  = poly_alloc(6);
/* variable Temporal */
poly_t *L2_t  = poly_alloc(6);
/* residuo temporal  */
poly_t *r = poly_alloc(2);
/* vars temporales */
poly_t *q_t = poly_alloc(3);
poly_t *r_t = poly_alloc(3);
/* obtenemos u*u */
poly_mul(D->u,D->u,u2);
/* almacenamos memoria para el sistema de 4 variables por 4 ecuaciones*/
m = sq_matrix_alloc(4);

/* rellenamos la parte derecha de la curva hipereliptica g=2 para manipularla x^5 -1 mod p y la multiplicamos por -1 (la restaremos a L*L)  */
hyperC->c[5]=mod_p(-1*1);
hyperC->c[4]=0;
hyperC->c[3]=0;                        //CORRECCION
hyperC->c[2]=0;
hyperC->c[1]=0;
hyperC->c[0]=1; //mod_p(-1*30);




a = D->u->c[1];
b = D->u->c[0];
c = D->v->c[1];
d = D->v->c[0];

/* la matriz resultante para duplicar un divisor:
 
a^2-b   	-a      	1	0	|c
ab      	-b      	0	1	|d
6a^2c-6ad-6bc   -4ac+4d         2c      0       |10ab-5a^3
6abc-6bd        -4bc            2d      0       |5b^2-5ba^2

*/
/* primer renglon con matriz aumentada */
m->c[0][0] = mod_p((a*a)-b);
m->c[0][1] = mod_p(-1*a);
m->c[0][2] = 1;
m->c[0][3] = 0;
coefs[0] = mod_p(c);

/* segundo renglon con matriz aumentada */
m->c[1][0] = mod_p(a*b);
m->c[1][1] = mod_p(-1*b);
m->c[1][2] = 0;
m->c[1][3] = 1;
coefs[1] = mod_p(d);

/* tercer renglon con matriz aumentada */
m->c[2][0] = mod_p((6*a*a*c)-(6*a*d)-(6*b*c));
m->c[2][1] = mod_p((4*d)-(4*a*c));
m->c[2][2] = mod_p(2*c);
m->c[2][3] = 0;
coefs[2] = mod_p((10*a*b)-(5*a*a*a));

/* cuarto renglon con matriz aumentada */
m->c[3][0] = mod_p((6*a*b*c)-(6*b*d));
m->c[3][1] = mod_p(-1*4*b*c);
m->c[3][2] = mod_p(2*d);
m->c[3][3] = 0;
coefs[3] = mod_p((5*b*b)-(5*b*a*a));
int i;
/*
for(i=0;i<4;i++)
        printf("%d ",coefs[i]);
printf("\n");
*/
system_solve(m,(int *)&coefs,(int *)&sols);
/*
printf("Matriz 2D\n");
sq_matrix_print(m);
*/
L->c[0] = sols[3];
L->c[1] = sols[2];
L->c[2] = sols[1];
L->c[3] = sols[0]; //CORRECCION: Transpuesto!!!! 
L->t = 3;
/*
printf("en 2D L=");
poly_print(L);
*/

/* Multiplicamos L*L */
poly_mul(L,L,L2);
/* Le restamos L2-C */
poly_add(L2,hyperC,L2_t);
/* Dividimos entre u*u para obtener u'(x), r deberia de ser el polinomio cero */
poly_div(L2_t,u2,DplusD->u,r_t);

/*
printf("r=");
poly_print(r);
*/
/* Calculamos v'(x) */
poly_div(L,DplusD->u,q_t,DplusD->v);
/*
printf("q_t=");
poly_print(q_t);
*/
/* proyectamos con infinito */
DplusD->v->c[0] = mod_p(-1*DplusD->v->c[0]);
DplusD->v->c[1] = mod_p(-1*DplusD->v->c[1]);
mumford_monic(DplusD);
return DplusD;

}


mumford_t *hyper_g2_add_disjoint_supp(mumford_t *D1, mumford_t *D2) {
sq_matrix_t *m;
mumford_t *D3 = mumford_alloc(2);
int a,b,c,d,A,B,C,D,coefs[4],sols[4];
/* Curva hipereliptica, parte derecha de  y^2 = x^5 - 1 */
poly_t *hyperC = poly_alloc(5);
/* aqui almacenaremos el polinomio de interpolacion de grado 3 */
poly_t *L = poly_alloc(3);
/* aqui almacenaremos u*u' */
poly_t *u1u2 = poly_alloc(D1->u->t+D2->u->t+1);
/* aqui almacenaremos L^2 = L*L - C donde C es la curva hipereliptica en este caso deg(L)=3 por eso almacenamos 6 y deg(C)=5 */
poly_t *L2  = poly_alloc(6);
/* variable Temporal */
poly_t *L2_t  = poly_alloc(6);
/* residuo temporal  */
poly_t *r = poly_alloc(2);
/* vars temporales */
poly_t *q_t = poly_alloc(3);
poly_t *r_t = poly_alloc(3);
/* obtenemos u*u' */
poly_mul(D1->u,D2->u,u1u2);
/* almacenamos memoria para el sistema de 4 variables por 4 ecuaciones*/
m = sq_matrix_alloc(4);

/* rellenamos la parte derecha de la curva hipereliptica g=2 para manipularla x^5 -1 mod p y la multiplicamos por -1 (la restaremos a L*L)  */
hyperC->c[5]=mod_p(-1*1);
hyperC->c[4]=0;
hyperC->c[3]=0;                        //CORRECCION
hyperC->c[2]=0;
hyperC->c[1]=0;
hyperC->c[0]=1; //mod_p(-1*30);

/* Solo para ayudarnos con la notacion, nombramos las variables a,b,c,d,A,B,C,D tal que D1=(x^2+ax+b,cx+d) D2=(x^2+Ax+B,Cx+D) */
a = D1->u->c[1];
b = D1->u->c[0];
c = D1->v->c[1];
d = D1->v->c[0];

A = D2->u->c[1];
B = D2->u->c[0];
C = D2->v->c[1];
D = D2->v->c[0];

/* la matriz resultante es:

a^2-b	-a	1	0	|c
ab	-b	0	1	|d
A^2-B	-A	1	0	|C
AB	-B	0	1	|D

*/


/* primer renglon con matriz aumentada */
m->c[0][0] = mod_p((a*a)-b);
m->c[0][1] = mod_p(-1*a);
m->c[0][2] = 1;
m->c[0][3] = 0;
coefs[0] = mod_p(c);

/* segundo renglon con matriz aumentada */
m->c[1][0] = mod_p(a*b);
m->c[1][1] = mod_p(-1*b);
m->c[1][2] = 0;
m->c[1][3] = 1;
coefs[1] = mod_p(d);

/* tercer renglon con matriz aumentada */
m->c[2][0] = mod_p((A*A)-B);
m->c[2][1] = mod_p(-1*A);
m->c[2][2] = 1;
m->c[2][3] = 0;
coefs[2] = mod_p(C);

/* cuarto renglon con matriz aumentada */
m->c[3][0] = mod_p(A*B);
m->c[3][1] = mod_p(-1*B);
m->c[3][2] = 0;
m->c[3][3] = 1;
coefs[3] = mod_p(D);

system_solve(m,(int *)&coefs,(int *)&sols);
/*
printf("Matriz D1+D2\n");
sq_matrix_print(m,coefs);
*/

L->c[0] = sols[3];
L->c[1] = sols[2];
L->c[2] = sols[1];
L->c[3] = sols[0]; //CORRECCION: Transpuesto!!!! 
L->t = 3;

//printf("L=");
//poly_print(L);

/* Multiplicamos L*L */
poly_mul(L,L,L2);
/* Le restamos L2-C */
poly_add(L2,hyperC,L2_t);
/* Dividimos entre u*u' para obtener u''(x), r deberia de ser el polinomio cero */
poly_div(L2_t,u1u2,D3->u,r_t);
/***
printf("r=");
poly_print(r);
printf("u''(x)=");
poly_print(D3->u);
****/
/* Calculamos v''(x) */
poly_div(L,D3->u,q_t,D3->v);
/**
printf("-v''(x)=");
poly_print(D3->v);
**/
/* proyectamos con infinito */
D3->v->c[0] = mod_p(-1*D3->v->c[0]);
D3->v->c[1] = mod_p(-1*D3->v->c[1]);
mumford_monic(D3);
return D3;
}

void div_copy(mumford_t *Dd, mumford_t *Ds) {
        int uts,vts;
        uts = Ds->u->t;
        vts = Ds->v->t;
        memset(Dd->u,0,sizeof(int)*uts);
        memset(Dd->v,0,sizeof(int)*vts);
        poly_copy(Dd->u,Ds->u);
        poly_copy(Dd->v,Ds->v);
        memcpy(&Dd->u->t,&Ds->u->t,sizeof(int));
        memcpy(&Dd->v->t,&Ds->v->t,sizeof(int));
        return;
}

int mumford_zero(mumford_t *D){
        //if(poly_iszero(D->u) && poly_iszero(D->v))
          //      return 1;
        return 0;
}

mumford_t *hyper_g2_add(mumford_t *D1, mumford_t *D2) {
        mumford_t *Df;
        if(mumford_zero(D1)) 
                return D2;
        if(mumford_zero(D2))
                return D1;

        if(mumford_equal(D1,D2)) {
                Df = hyper_g2_double_sup_2(D1);
                mumford_monic(Df);
        }
        else { 
                Df = hyper_g2_add_disjoint_supp(D1,D2);
                mumford_monic(Df);
        }
        return Df;
}

mumford_t *hyper_g2_nadd(unsigned int n,mumford_t *D) {
        mumford_t *T1,*T2,*Df;
        int i;
        T1 = mumford_alloc(2);
        T2 = mumford_alloc(2);
        Df = mumford_alloc(2);
        if(n == 1) 
                return D;
        if(n == 0)
                return T1;
        div_copy(T1,D);
        for(i=0;i<n-1;i++) {
                /*
                printf("Iteracion %d\n",i);
                printf("D=");
                divisorform_g2_print(D);
                printf("T1=");
                divisorform_g2_print(T1);
                */
                
                T2 = hyper_g2_add(D,T1);
                /*
                printf("D+T1=");
                divisorform_g2_print(T2);
                divisorform_g2_print(T2);
                */
                
                div_copy(T1,T2);
                //printf("\n\n");
        }
        return T2;

}

int
main (int argc, char **argv)
{
mumford_t *Db,*Apub,*Bpub,*SA,*SB;
int apriv,bpriv;


/* los 2 puntos de la base */
point_t *a1,*b1;

/* alojamos memoria */
a1 = point_alloc();
b1 = point_alloc();

/* puntos K-racionales de la curva
 * (1,0),(2,0),(3,5),(3,26),(4,0),(6,5),(6,26),(7,2),(7,29),(8,0),(11,6),(11,25),(12,5),(12,26),(13,6),(13,25)
 * (14,2),(14,29),(16,0),(17,5),(17,26),(19,2),(19,29),(21,6),(21,25),(22,6),(22,25),(24,5),(24,26),(25,2)
 * (25,29),(26,6),(26,25),(28,2),(28,29)
 */

/*
D   =(7,2)+(28,29)-2inf
T1  =(3,5)+(21,25)-2inf
*/



/* punto 1 para D1 es (13,6) en la curva */
a1->x = 12;
a1->y = 5;

b1->x = 28;
b1->y = 2;
/*
point_free(a1);
point_free(b1);
point_free(a2);
point_free(b2);
*/

/* Construimos los divisores con los puntos racionales */
Db = mumford_build_g2(a1,b1) ;
/*
mumford_print(D1);
mumford_print(D2);


divisorform_g2_print(D1);
divisorform_g2_print(D2);
*/


/* apriv = 16 , bpriv=50-16 */

apriv=16;
bpriv=34;

Apub = hyper_g2_nadd(apriv,Db);
printf("Llave privada de A=%d y publica:\n",apriv);
mumford_print(Apub);
divisorform_g2_print(Apub);


Bpub = hyper_g2_nadd(bpriv,Db);
printf("Llave privada de B=%d y publica:\n",bpriv);
mumford_print(Bpub);
divisorform_g2_print(Bpub);

/* A y B intercambian llaves publicas y calculan */

/* A calcula de la llave publica de B (Bpub) */
printf("Secreto obtenido por A es\n######## BEGIN SECRETO A ########\n");
SA = hyper_g2_nadd(apriv,Bpub);
mumford_print(SA);
divisorform_g2_print(SA);
printf("###### END SECRETO A ##########\n");

/* B calcula de la llave publica de A (Apub) */
printf("Secreto obtenido por B es:\n####### BEGIN SECRETO B #########\n");
SB = hyper_g2_nadd(bpriv,Apub);
mumford_print(SB);
divisorform_g2_print(SB);
printf("###### END SECRETO B ##########\n");







/*
D4 = hyper_g2_add_disjoint_supp(D2,D2);
mumford_print(D4);
divisorform_g2_print(D4);
*/


/*
printf("Double normal\n");
D1 = hyper_g2_double_sup_2(D2);
mumford_print(D1);
divisorform_g2_print(D1);
*/


/* sumamos :) */
/*
D3 =  hyper_g2_double_sup_2(D1);
printf("2[D1]:\n");
mumford_print(D3);
divisorform_g2_print(D3);
D2 =  hyper_g2_add_disjoint_supp(D3,D2);
printf("3[D1]:\n");
mumford_print(D2);
divisorform_g2_print(D2);
*/
return 0;
}
