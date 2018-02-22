#define h(x,y) ((((y*y)%31) - ((x*x*x*x*x)%31) + 1)+(31*1000))%31
int main(void) {
        int i,j;
        for(i=0;i<31;i++)
                for(j=0;j<31;j++)
                        if(h(i,j) == 0)
                                printf("(%d,%d)\n",i,j);
        return 0;

}
