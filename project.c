/********************************************************************
Title: Program on Splitting or Fractional Step Method
*********************************************************************/

#include<stdio.h>
#include<conio.h>
#include<math.h>
main()
{

int i,j,k,n,p,q;
float T0[27][1], T1[27][1], T2[27][1], T[27][1], C1[29][29], C2[29][29], C3[29][29];
float C[5][5][5], B[27][27], temp;

for(p=0; p<29; p++){
	for(q=0; q<29; q++){
		C1[p][q] = 0;
		C2[p][q] = 0;
		C3[p][q] = 0;
			   }
		   }
clrscr();

for(n=0; n<27; n++){
       T0[n][1] = 300; //Initial Condition
       printf("T0[%d] = %f\n",n,T0[n]);
       }

p=1;

 /******************************************************************
		    Step 1 Begins
       Contructing Coefficient Matrix for Step 1
 *******************************************************************/

for(i=0; i<5; i++){
	for(j=0; j<5; j++){
		for(k=0; k<5; k++){
			C[i][j][k] = 0;
				  }
			  }
		  }

for(i=1; i<4; i++)
	{ for(j=1; j<4; j++)
		{ for(k=1; k<4; k++)
			{ if(i==1){
				    C[1][j][k] = 1;
				    C[2][j][k] = -0.001067;
				  }
			  else{
			  if(i==3){
				    C[i][j][k] = 1;
				  }
			  else{
				C[i][j][k] = 1;
				C[i-1][j][k] = -0.000533;
				C[i+1][j][k] = -0.000533;
			      }
			      }
			  C1[p][p] = C[i][j][k];
			  C1[p][p+1] = C[i+1][j][k];
			  C1[p][p-1] = C[i-1][j][k];
			  p++;
			}
		}
	}
/*********** Finding Inverse of the Coefficient Matrix Step 1***************/

for (i=0; i<27; i++)
    {
	for (j=0; j<27; j++)
	{
	    if ( i == j)
		B[i][j] = 1;
	    else
		B[i][j] = 0;
	}
    }

    if (C1[1][1]==0)
	{
	for(i=1;i<29;i++)
	    {
		temp=C1[1][i];
		C1[1][i]=C1[2][i];
		C1[2][i]=temp;

		temp=B[0][i-1];
		B[0][i-1]=B[1][i-1];
		B[1][i-1]=temp;
	    }
	}

    temp = C1[1][1];

    for(i=1;i<29;i++)
    {
	for(j=1;j<29;j++)
	{
	    C1[i][j]/=temp;
	    B[i-1][j-1]/=temp;
	}
    }

    for (i=2;i<29;i++)
    {
	temp = C1[i][0];
	for(j=1;j<29;j++)
	{
	    C1[i][j]-=C1[1][j]*temp;
	    B[i-1][j-1]-=B[0][j-1]*temp;
	}
    }
 /************** Matrix Multiplication for Step 1 ***************/

for(i=0;i<27;i++)
        {
            for(j=0;j<1;j++)
            {
                T1[i][j]=0;
                for(k=0;k<27;k++)
                {
                    T1[i][j]+=B[i][k]*T0[k][j];

                }
                printf("%d\t",T1[i][j]);
            }
            printf("\n");
	}

 /*********************************************************************
		    Step 2 Begins
       Contructing Coefficient Matrix for Step 2
 *********************************************************************/

for(i=0; i<5; i++){
	for(j=0; j<5; j++){
		for(k=0; k<5; k++){
			C[i][j][k] = 0;
				  }
			  }
		  }

for(i=1; i<4; i++)
	{ for(j=1; j<4; j++)
		{ for(k=1; k<4; k++)
			{ if(j==1){
				    C[i][1][k] = 1;
				    C[i][2][k] = -0.001067;
				  }
			  else{
			  if(j==3){
				    C[i][j][k] = 1;
				  }
			  else{
				C[i][j][k] = 1;
				C[i][j-1][k] = -0.000533;
				C[i][j+1][k] = -0.000533;
			      }
			      }
			  C2[p][p] = C[i][j][k];
			  C2[p][p+1] = C[i][j+1][k];
			  C2[p][p-1] = C[i][j+1][k];
			  p++;
			}
		}
	}


/*********** Finding Inverse of the Coefficient Matrix Step 2***************/

for (i=0; i<27; i++)
    {
	for (j=0; j<27; j++)
	{
	    if ( i == j)
		B[i][j] = 1;
	    else
		B[i][j] = 0;
	}
    }

    if (C2[1][1]==0)
	{
	for(i=1;i<29;i++)
	    {
		temp=C2[1][i];
		C2[1][i]=C2[2][i];
		C2[2][i]=temp;

		temp=B[0][i-1];
		B[0][i-1]=B[1][i-1];
		B[1][i-1]=temp;
	    }
	}

    temp = C2[1][1];

    for(i=1;i<29;i++)
    {
	for(j=1;j<29;j++)
	{
	    C2[i][j]/=temp;
	    B[i-1][j-1]/=temp;
	}
    }

    for (i=2;i<29;i++)
    {
	temp = C2[i][0];
	for(j=1;j<29;j++)
	{
	    C1[i][j]-=C1[1][j]*temp;
	    B[i-1][j-1]-=B[0][j-1]*temp;
	}
    }
/************** Matrix Multiplication for Step 2 ***************/

for(i=0;i<27;i++)
        {
            for(j=0;j<1;j++)
            {
                T2[i][j]=0;
                for(k=0;k<27;k++)
                {
                    T2[i][j]+=B[i][k]*T1[k][j];

                }
                printf("%d\t",T2[i][j]);
            }
            printf("\n");
	}
 
 /********************************************************************
		    Step 3 Begins
       Contructing Coefficient Matrix for Step 3
 *********************************************************************/

for(i=0; i<5; i++){
	for(j=0; j<5; j++){
		for(k=0; k<5; k++){
			C[i][j][k] = 0;
				  }
			  }
		  }

for(i=1; i<4; i++)
	{ for(j=1; j<4; j++)
		{ for(k=1; k<4; k++)
			{ if(k==1){
				    C[i][j][1] = 1;
				    C[k][j][2] = -0.001067;
				  }
			  else{
			  if(k==3){
				    C[i][j][k] = 1;
				  }
			  else{
				C[i][j][k] = 1;
				C[i][j][k-1] = -0.000533;
				C[i][j][k+1] = -0.000533;
			      }
			      }
			  C3[p][p] = C[i][j][k];
			  C3[p][p+1] = C[i][j][k+1];
			  C3[p][p-1] = C[i][j][k-1];
			  p++;
			}
		}
	}


/*********** Finding Inverse of the Coefficient Matrix Step 3***************/

for (i=0; i<27; i++)
    {
	for (j=0; j<27; j++)
	{
	    if ( i == j)
		B[i][j] = 1;
	    else
		B[i][j] = 0;
	}
    }

    if (C3[1][1]==0)
	{
	for(i=1;i<29;i++)
	    {
		temp=C3[1][i];
		C3[1][i]=C3[2][i];
		C3[2][i]=temp;

		temp=B[0][i-1];
		B[0][i-1]=B[1][i-1];
		B[1][i-1]=temp;
	    }
	}

    temp = C3[1][1];

    for(i=1;i<29;i++)
    {
	for(j=1;j<29;j++)
	{
	    C3[i][j]/=temp;
	    B[i-1][j-1]/=temp;
	}
    }

    for (i=2;i<29;i++)
    {
	temp = C3[i][0];
	for(j=1;j<29;j++)
	{
	    C3[i][j]-=C3[1][j]*temp;
	    B[i-1][j-1]-=B[0][j-1]*temp;
	}
    }
/************** Matrix Multiplication for Step 3 ***************/

for(i=0;i<27;i++)
	{
	    for(j=0;j<1;j++)
	    {
		T[i][j]=0;
		for(k=0;k<27;k++)
		{
		    T[i][j]+=B[i][k]*T2[k][j];

		}
/************************ Output of the Solution *********************/

		printf("%f\t",T[i][j]);
	    }
	    printf("\n");
	}

getch();
return(0);
}
