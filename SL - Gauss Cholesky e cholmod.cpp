/**
    - Eliminação de Gauss
    - Clolesky
    - Cholesky modificado
*/

#include <stdio.h>
#include <math.h>

#define epsilon 0.0001
/** METODO DE ELIMINACAO DE GAUSS
    Fonte: Wikipedia ( adaptamos a nossas necessidades! )
**/

int forwardSubstitution(int n, double **a) {
	int i, j, k, max;
	double t;
	for (i = 0; i < n; ++i) {
		max = i;
		for (j = i + 1; j < n; ++j)
			if (a[j][i] > a[max][i])
				max = j;
		for (j = 0; j < n + 1; ++j) {
			t = a[max][j];
			a[max][j] = a[i][j];
			a[i][j] = t;
		}
		for (j = n; j >= i; --j){
			for (k = i + 1; k < n; ++k){
			    if (a[i][i]< 10e-8){
                    return (-1);
			    }
				a[k][j] -= a[k][i]/a[i][i] * a[i][j];
			}
		}
	}
	return (0);
}

int reverseElimination(int n, double **a) {
	int i, j;
	for (i = n - 1; i >= 0; --i) {
	    if (a[i][i] < 10e-8){
            return (0);
	    }
		a[i][n] = a[i][n] / a[i][i];
		a[i][i] = 1;
		for (j = i - 1; j >= 0; --j) {
			a[j][n] -= a[j][i] * a[i][n];
			a[j][i] = 0;
		}
	}
	return (0);
}

void gauss(int n, double **Efx, double *d) {
	int i, j;

    double **a = new double*[n+1];
    for (int i = 0 ; i < n ; i++) {
        a[i] = new double[n+1];
    }

    for (int i = 0 ; i < n ; i++) {
        for (int j = 0 ; j < n ; j++) {
            a[i][j] = Efx[i][j];
        }
    }

	for (int j = 0 ; j < n ; j++) {
        a[j][n] = d[j];
    }

	while (forwardSubstitution(n,a) == -1){

        for (int i = 0 ; i < n ; i++){
            a[i][i]+=10;
        }
	}
	while (reverseElimination(n,a) == -1){

        for (int i = 0 ; i < n ; i++){
            a[i][i]+=10;
        }
	}

    for (int j = 0 ; j < n ; j++) {
        d[j] = a[j][n];
    }

    for (int i = 0 ; i < n ; i++) {
        delete []a[i];
    }
    delete []a;
}

/// FIM DA ELIMINACAO DE GAUSS


int cholesky (double **A, int n){
    double **G = new double*[n];
    for (int i = 0 ; i < n ; i++){
        G[i] = new double[n];
        for (int j = 0 ; j < n ; j++){
            G[i][j] = 0;
        }
    }
    double sq;
    for (int i = 0 ; i < n ; i++){

        if (i == 0){
            G[0][0] = sqrt (A[0][0]);
            for (int j = 1 ; j < n ; j++){
                G[j][0] = A[j][0]/G[0][0];
            }
        }

        else {
            for (int j = 1 ; j < i ; j++){
                sq = 0;
                for (int k = 0 ; k<j ; k++){
                    sq = sq + G[i][k]*G[j][k];
                }
                G[i][j] = (A[i][j] - sq)/G[j][j];
            }

            sq = 0;
            for (int k = 0 ; k < i ; k++){
                sq = sq + G[i][k]*G[i][k];
            }
            sq = A[i][i] - sq;
            G[i][i] = sqrt(sq);
        }
    }


    for (int i = 0 ; i < n ; i++){
        for (int j = 0 ; j < n ; j++){
            A[i][j] = G[i][j];
        }
    }
    return (0);
}

int chol_mod (double **A, int n){
    double **G = new double*[n];
    for (int i = 0 ; i < n ; i++){
        G[i] = new double[n];
        for (int j = 0 ; j < n ; j++){
            G[i][j] = 0;
        }
    }
    double sq;
    for (int i = 0 ; i < n ; i++){

        if (i == 0){
            G[0][0] = sqrt (fabs(A[0][0]));

            for (int j = 1 ; j < n ; j++){
                if (fabs(G[0][0])<10e-8){
                    if (G[0][0]>=0)
                        G[0][0] = epsilon;
                    else
                        G[0][0] = -epsilon;
                }
                G[j][0] = A[j][0]/G[0][0];
            }
        }

        else {
            for (int j = 1 ; j < i ; j++){
                sq = 0;
                for (int k = 0 ; k<j ; k++){
                    sq = sq + G[i][k]*G[j][k];
                }
                if (fabs(G[j][j])<10e-8){
                    if (G[j][j]>=0)
                        G[j][j] = epsilon;
                    else
                        G[j][j] = -epsilon;
                }
                G[i][j] = (A[i][j] - sq)/G[j][j];
            }

            sq = 0;
            for (int k = 0 ; k < i ; k++){
                sq = sq + G[i][k]*G[i][k];
            }
            sq = A[i][i] - sq;
            G[i][i] = sqrt(fabs (sq));

        }
    }
    for (int i = 0 ; i < n ; i++){
        for (int j = 0 ; j < n ; j++){
            A[i][j] = G[i][j];
        }
    }
    return (0);
}

int getd (double **A, double *d, int n){
    double aux;
    for (int i = 0 ; i < n ; i++){
        aux = 0.0;
        for (int j = 0 ; j < i ; j++){
            aux = aux + A[i][j]*d[j];
        }
        if (fabs(A[i][i]) < 10e-8){
            A[i][i] = epsilon;
        }
        d[i] = (d[i]-aux)/A[i][i];
    }

    for (int i = n-1 ; i >=0 ; i--){
        aux = 0.0;
        for (int j = i+1 ; j < n ; j++){
            aux = aux + A[j][i]*d[j];
        }
        if (fabs(A[i][i]) < 10e-8){
            A[i][i] = epsilon;
        }
        d[i] = (d[i]-aux)/A[i][i];
    }

    return 0;
}


/// TESTE DE UNIDADE
//int main (){
//    int n = 3;
//    double **A = new double*[n];
//    for (int i = 0 ; i < n ; i++) {
//        A[i] = new double[n];
//    }
//    double *d = new double [n];
//    A[0][0] = 4;
//    A[1][0] = 2;
//    A[0][1] = 2;
//    A[1][1] = 10;
//    A[0][2] = -4;
//    A[2][0] = -4;
//    A[2][1] = 4;
//    A[1][2] = 4;
//    A[2][2] = 9;
//
//    d[0] = 0;
//    d[1] = 6;
//    d[2] = 5;

//    chol_mod (A, n);
//    getd (A, d, n);

//    for (int i = 0 ; i < n ; i++){
//        printf ("\n");
//        for (int j = 0 ; j < n ; j++){
//            printf ("%lf, ",A[i][j]);
//        }
//    }
//
//    printf ("\n");printf ("\n");
//    for (int i = 0 ; i < n ; i++){
//        printf ("%lf, ",d[i]);
//    }
//}
