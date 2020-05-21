#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#define MAX 100 
int i,j,k,n=3,R=3,C=3,M=3,N=4,sum,order;

float gaussian(float A[10][10]){
	float c,x[10],sum=0.0;
    printf("\nEnter the order of matrix: ");
    scanf("%d",&n);
    printf("\nEnter the elements of augmented matrix R-wise:\n\n");
    for(i=1; i<=n; i++)
    {
        for(j=1; j<=(n+1); j++)
        {
            printf("A[%d][%d] : ", i,j);
            scanf("%f",&A[i][j]);
        }
    }
    for(j=1; j<=n; j++) /* loop for the generation of upper triangular matrix*/
    {
        for(i=1; i<=n; i++)
        {
            if(i>j)
            {
                c=A[i][j]/A[j][j];
                for(k=1; k<=n+1; k++)
                {
                    A[i][k]=A[i][k]-c*A[j][k];
                }
            }
        }
    }
    x[n]=A[n][n+1]/A[n][n];
    /* this loop is for backward substitution*/
    for(i=n-1; i>=1; i--)
    {
        sum=0;
        for(j=i+1; j<=n; j++)
        {
            sum=sum+A[i][j]*x[j];
        }
        x[i]=(A[i][n+1]-sum)/A[i][i];
    }
    printf("\nThe solution is: \n");
    for(i=1; i<=n; i++)
    {
        printf("\nx%d=%f\t",i,x[i]); /* x1, x2, x3 are the required solutions*/
    }
    return(0);
}	

float gaussjordan(float matrix[10][10]){
    float ratio,a;
    printf("Enter order of matrix: ");
    scanf("%d", &n);
    printf("Enter the matrix: \n");
    for(i = 0; i < n; i++){
        for(j = 0; j < n; j++){
            scanf("%f", &matrix[i][j]);
        }
    }
    for(i = 0; i < n; i++){
        for(j = n; j < 2*n; j++){
            if(i==(j-n))
                matrix[i][j] = 1.0;
            else
                matrix[i][j] = 0.0;
        }
    }
    for(i = 0; i < n; i++){
        for(j = 0; j < n; j++){
            if(i!=j){
                ratio = matrix[j][i]/matrix[i][i];
                for(k = 0; k < 2*n; k++){
                    matrix[j][k] -= ratio * matrix[i][k];
                }
            }
        }
    }
    for(i = 0; i < n; i++){
        a = matrix[i][i];
        for(j = 0; j < 2*n; j++){
            matrix[i][j] /= a;
        }
    }
    printf("The inverse matrix is: \n");
    for(i = 0; i < n; i++){
        for(j = n; j < 2*n; j++){
            printf("%.2f", matrix[i][j]);
            printf("\t");
        }
        printf("\n");
    }
    return 0;
}

float power(float A[40][40]){
float x[40],z[40],e[40],zmax,emax;
    printf("\nEnter the order of matrix:");
    scanf("%d",&n);
    printf("\nEnter matrix elements R-wise\n");
    for(i=1; i<=n; i++)
    {
        for(j=1; j<=n; j++)
        {
            printf("A[%d][%d]=", i,j);
            scanf("%f",&A[i][j]);
        }
    }
    printf("\nEnter the Cumn vector\n");
    for(i=1; i<=n; i++)
    {
        printf("X[%d]=",i);
        scanf("%f",&x[i]);
    }
    do
    {
        for(i=1; i<=n; i++)
        {
            z[i]=0;
            for(j=1; j<=n; j++)
            {
                z[i]=z[i]+A[i][j]*x[j];
            }
        }
        zmax=fabs(z[1]);
        for(i=2; i<=n; i++)
        {
            if((fabs(z[i]))>zmax)
                zmax=fabs(z[i]);
        }
        for(i=1; i<=n; i++)
        {
            z[i]=z[i]/zmax;
        }
        for(i=1; i<=n; i++)
        {
            e[i]=0;
            e[i]=fabs((fabs(z[i]))-(fabs(x[i])));
        }
        emax=e[1];
        for(i=2; i<=n; i++)
        {
            if(e[i]>emax)
                emax=e[i];
        }
        for(i=1; i<=n; i++)
        {
            x[i]=z[i];
        }
    }
    while(emax>0.001);
    printf("\n The required eigen value is %f",zmax);
    printf("\n\nThe required eigen vector is :\n");
    for(i=1; i<=n; i++)
    {
        printf("%f\t",z[i]);
    }
	printf("\n");
    //getch();
}

void luDecomposition(int mat[][MAX], int n) { 
	int lower[n][n], upper[n][n];
	for(int i = 0; i < n; i++) 
    for( int j = 0; j < n; j++) 
    { 
        lower[i][j]= 0; 
        upper[i][j]=0; 
    } 

	// Decomposing matrix into Upper and Lower 
	// triangular matrix 
	for (int i = 0; i < n; i++) { 

		// Upper Triangular 
		for (int k = i; k < n; k++) { 

			// Summation of L(i, j) * U(j, k) 
			int sum = 0; 
			for (int j = 0; j < i; j++) 
				sum += (lower[i][j] * upper[j][k]); 

			// Evaluating U(i, k) 
			upper[i][k] = mat[i][k] - sum; 
		} 

		// Lower Triangular 
		for (int k = i; k < n; k++) { 
			if (i == k) 
				lower[i][i] = 1; // Diagonal as 1 
			else { 

				// Summation of L(k, j) * U(j, i) 
				int sum = 0; 
				for (int j = 0; j < i; j++) 
					sum += (lower[k][j] * upper[j][i]); 

				// Evaluating L(k, i) 
				lower[k][i] = (mat[k][i] - sum) / upper[i][i]; 
			} 
		} 
	}
 
	printf("Lower Triangular\t\t");
		printf("Upper Triangular\n"); 

	// Displaying the result : 
	for (int i = 0; i < n; i++) 
	{ 
		// Lower 
		for (int j = 0; j < n; j++) 
			printf("%d\t",lower[i][j]); 
		printf("\t"); 

		// Upper 
		for (int j = 0; j < n; j++) 
			printf("%d\t",upper[i][j]);  
		printf("\n"); 
	} 
}

void swap(int mat[R][C], int R1, int R2, int C) { 
	for (int i = 0; i < C; i++) 
	{ 
		int temp = mat[R1][i]; 
		mat[R1][i] = mat[R2][i]; 
		mat[R2][i] = temp; 
	} 
} 

void display(int mat[R][C], int R, int C) { 
	for (int i = 0; i < R; i++) 
	{ 
		for (int j = 0; j < C; j++) 
			printf(" %d", mat[i][j]); 
		printf("\n"); 
	} 
} 

int rankOfMatrix(int mat[R][C]) { 
	int rank = C; 

	for (int row = 0; row < rank; row++) 
	{ 
		// Before we visit current row 'row', we make 
		// sure that mat[row][0],....mat[row][row-1] 
		// are 0. 

		// Diagonal element is not zero 
		if (mat[row][row]) 
		{ 
		for (int col = 0; col < R; col++) 
		{ 
			if (col != row) 
			{ 
				// This makes all entries of current 
				// column as 0 except entry 'mat[row][row]' 
				double mult = (double)mat[col][row] / 
									mat[row][row]; 
				for (int i = 0; i < rank; i++) 
				mat[col][i] -= mult * mat[row][i]; 
			} 
		} 
		} 

		// Diagonal element is already zero. Two cases 
		// arise: 
		// 1) If there is a row below it with non-zero 
		// entry, then swap this row with that row 
		// and process that row 
		// 2) If all elements in current column below 
		// mat[r][row] are 0, then remvoe this column 
		// by swapping it with last column and 
		// reducing number of columns by 1. 
		else
		{ 
			int reduce = 1;//true; 

			/* Find the non-zero element in current 
				column */
			for (int i = row + 1; i < R; i++) 
			{ 
				// Swap the row with non-zero element 
				// with this row. 
				if (mat[i][row]) 
				{ 
					swap(mat, row, i, rank); 
					reduce = -1; 
					break ; 
				} 
			} 

			// If we did not find any row with non-zero 
			// element in current columnm, then all 
			// values in this column are 0. 
			if (reduce) 
			{ 
				// Reduce number of columns 
				rank--; 

				// Copy the last column here 
				for (int i = 0; i < R; i ++) 
					mat[i][row] = mat[i][rank]; 
			} 

			// Process this row again 
			row--; 
		} 

	// Uncomment these lines to see intermediate results 
	printf("Intermediate Steps\n");
	display(mat, R, C); 
	printf("\n"); 
	} 
	return rank; 
} 

void transpose ( int mat[3][3], int m[3][3] ){
    int i, j ;
    for ( i = 0 ; i < MAX ; i++ )
    {
        for ( j = 0 ; j < MAX ; j++ )
             m[i][j] = mat[j][i] ;
    }
}
 
int echelon(){
	printf("Enter order of the matrix:");
	int order;scanf("%d",&order);
	int row=order;int col=order;
	float mat[row][col], ratio,  det = 1, factor = 1, temp;
	int i, j, k, m, row_intrchng_count = 0,shift = 0;
	char c ='r';
while(c != 'q')
{             
   //printf("Enter order of the matrix:");
   //scanf("%d",&order);
   //row=order;col=order;
   printf("Enter the matrix of order %dx%d\n", row, col );
   for(i = 0; i<row; i++)
                for(j = 0; j<col; j++)
                scanf("%f", &mat[i][j]);


   for(i = 0; i<row-1; i++)//for loop for gettting row echelon form of the matrix
   {
                for(j = i; j<col; j++)
      {
         if(j>i)
         {
                if(mat[i][i+shift] != 0)//mat[i][i+shift] is pivot element of the matrix
            {
                if(mat[j][i+shift] != 0)//check whether the value is zero or needed to be made zero
               {


                               temp = mat[j][i+shift];
                               factor *= mat[i][i+shift];
                               for(k =0; k<col; k++)//elementary row operation of this loop does not change the value of the determinant
                               {
                                    mat[j][k] *= mat[i][i+shift];
                               }
                                for(k =0; k<col; k++)//elementary row operation of this loop does not change the value of the determinant
                               {
                                    mat[j][k] -= mat[i][k]*temp;
                               }
               }
                                }
            else
            {
                for(k = i+1; k<row; k++)
               if(mat[k][i] != 0)	//satisfaction of this condition means interchange of two rows and so determinant is multiplied by -1 per swapping of the two rows
                  {
                                row_intrchng_count++;
                     break;
                  }
               for(m = 0; m<col; m++)//this for loop interchange the two rows of the matrix
               {
               mat[i][m] = mat[i][m]+mat[k][m];
                  mat[k][m] = mat[i][m]-mat[k][m];
                  mat[i][m] = mat[i][m]-mat[k][m];
               }
               j--; //to cancel the effect of the increment that take place after the end } of j loop
               if(k == row)    //this condition satisfaction means that matrix is singular
               {

                  mat[row-1][row-1] = 0;
                  i =row;
                  break;
               }

            }
         }
      }
   }

    if(mat[row-1][row-1] != 0)
    {
                for(i = 0; i<row; i++)
                det *= mat[i][i];
                det = pow(-1, row_intrchng_count)*det/factor;
    }
    else
                det = 0;
    printf("the determinant of the matrix = %f\n", det);

   for(i = 0; i<row; i++)
                for(j = 0; j<col; j++)
      {
                printf("%5.3f\t", mat[i][j]);
                if(j == (col-1))
                                printf("\n");
      }

   scanf("%c", &c);
}
   getch();
}

int isOrthogonal(int a[][MAX], int m, int n) { 
if (m != n) 
	return -1; 

// Find transpose 
int trans[n][n]; 
for (int i = 0; i < n; i++) 
	for (int j = 0; j < n; j++) 
	trans[i][j] = a[j][i]; 

// Find product of a[][] 
// and its transpose 
int prod[n][n]; 
for (int i = 0; i < n; i++) 
{ 
	for (int j = 0; j < n; j++) 
	{ 

	int sum = 0; 
	for (int k = 0; k < n; k++) 
	{ 

		// Since we are multiplying with 
		// transpose of itself. We use 
		sum = sum + (a[i][k] * a[j][k]); 
	} 

	prod[i][j] = sum; 
	} 
} 

// Check if product is identity matrix 
for (int i = 0; i < n; i++) 
{ 
	for (int j = 0; j < n; j++) 
	{ 
	if (i != j && prod[i][j] != 0) 
		return -1; 
	if (i == j && prod[i][j] != 1) 
		return -1; 
	} 
} 

return 1; 
} 
 
int determinant( int mat[3][3] ){
    sum=0;int p ;
    j = 1 ; k = MAX - 1 ;

    for ( i = 0 ; i < MAX ; i++ )
    {
        p = pow ( -1, i ) ;

        if ( i == MAX - 1 )
            k = 1 ;
        sum = sum + ( mat[0][i] * ( mat[1][j] *
                                    mat[2][k] - mat[2][j] *
                                    mat[1][k] ) ) * p ;
        j = 0 ;
    }

    return sum ;
}
 
int main()
{
	int c;
	while(c!=10)
	{
		float matrix[n][n],A[n][n],b[n];int trans[N][N];
		printf("Please select your choice:\n");
		printf("1:Singular/Non Singular\n2:Gaussian Elimination\n3:Transpose of a matrix\n4:Orthogonality of Matrix\n5:Echelon form of Matrix\n6:LU Decomposition\n7:Rank of a Matrix\n8:Rayleigh Power Method\n9:Gauss Jordan Method\n10.Exit\n");
		scanf("%d",&c);
		switch(c){
		case 1: printf("Singular/Non Singular:\n");
			{
				//printf()
				int a[3][3], i, j;
				long determinant;
				printf("Enter the 9 elements of matrix: ");
				for(i = 0 ;i < 3;i++)
				  for(j = 0;j < 3;j++)
					   scanf("%d", &a[i][j]);
				printf("\nThe matrix is\n");
				for(i = 0;i < 3; i++){
					printf("\n");
					for(j = 0;j < 3; j++)
						printf("%d\t", a[i][j]);
			}
			determinant = a[0][0] * ((a[1][1]*a[2][2]) - (a[2][1]*a[1][2])) -a[0][1] * (a[1][0]* a[2][2] - a[2][0] * a[1][2]) + a[0][2] * (a[1][0] * a[2][1] - a[2][0] * a[1][1]);
			if ( determinant == 0)
				printf("\nMatrix is Singular");
			else
				printf("\nThe given matrix is not Singular");
		return 0;
			}
			break;
		case 2: printf("Gaussian Elimination:");
				gaussian(A);
				break;
		case 3: printf("Transpose:\n");
		{
			int m, n, c, d, matrix[10][10], transpose[10][10];
			printf("Enter the number of rows and columns of matrix\n");
			scanf("%d%d", &m, &n);
			printf("Enter elements of the matrix\n");
			for (c = 0; c < m; c++)
				for(d = 0; d < n; d++)
					scanf("%d", &matrix[c][d]);
 
		    for (c = 0; c < m; c++)
			    for( d = 0 ; d < n ; d++ )
					transpose[d][c] = matrix[c][d];
 
			printf("Transpose of the matrix:\n");
 
		    for (c = 0; c < n; c++) {
			    for (d = 0; d < m; d++)
					printf("%d\t", transpose[c][d]);
			printf("\n");
		}
		return 0;
		}
		break;
		case 4: printf("Orthogonality of Matrix:\n");
				{ 

					int a[][MAX] = {{1, 0, 0}, 
									{0, 1, 0}, 
									{0, 0, 1}}; 
				
					if (isOrthogonal(a, 3, 3)) 
					printf("Yes"); 
					else
					printf("No"); 
					return 0; 
				}
				break;
		case 5: printf("Echelon form of Matrix:\n");
				echelon();
				break;
		case 6: printf("LU Decomposition:\n");
				{ 
				int mat[][MAX] = { { 2, -1, -2 }, 
						   { -4, 6, 3 }, 
						   { -4, -2, 8 } }; 
	  
				luDecomposition(mat, 3); 
				return 0; 
				}  	
				break;
		case 7: printf("Rank of a Matrix\n");
				{ 
					int mat[][3] = {{1,   2,   1}, 
					{0,  3,   1}, 
					{-2,   1,   4}}; 
					printf("Rank of the matrix is : %d", 
					rankOfMatrix(mat)); 
					return 0; 
				}
				break;
		case 8: printf("Rayleigh Power Method");
				power(A);
				break;
		case 9: printf("Gauss Jordan Method\n");
				gaussjordan(matrix);
				break;
		}
	}
}