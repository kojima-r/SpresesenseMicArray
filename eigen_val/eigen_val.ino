#include <cstdio>
#include <cmath>
#include <cstring>

struct Complex {
    float real;
    float imag;

    Complex() : real(0), imag(0) {}
    Complex(float r, float i) : real(r), imag(i) {}

    Complex operator+(const Complex& other) const {
        return Complex(real + other.real, imag + other.imag);
    }

    Complex operator-(const Complex& other) const {
        return Complex(real - other.real, imag - other.imag);
    }

    Complex operator*(const float& other) const {
    	return Complex(real*other, imag*other);
    }
    Complex operator*(const Complex& other) const {
        return Complex(real * other.real - imag * other.imag, real * other.imag + imag * other.real);
    }

    Complex operator/(const float& scalar) const {
        return Complex(real / scalar, imag / scalar);
    }
};

inline Complex conj(const Complex& a) {
    return Complex(a.real,-a.imag);
}
// Function to print a 2D array
void printArray(Complex* arr, int rows, int cols) {
    for (int i = 0; i < rows; ++i) {
        for (int j = 0; j < cols; ++j) {
            printf("%2.2f + %2.2fi ", arr[i * cols + j].real, arr[i * cols + j].imag);
        }
        printf("\n");
    }
}

int run() {
    const int N = 3; // Assuming N is 3, adjust accordingly

    Complex* a = new Complex[N * N];
    Complex* b = new Complex[N * N];
    Complex* A = new Complex[N * N];
    Complex* X = new Complex[N * N];
    Complex* An = new Complex[N * N];
    Complex* C = new Complex[N * N];

    // Initialize a matrix
    a[0] = {1.0, 0.0}; a[1] = {5.0, 0.0}; a[2] = {4.0, 0.0};
    a[3] = {2.0, 0.0}; a[4] = {4.0, 0.0}; a[5] = {-7.0, 0.0};
    a[6] = {2.0, 0.0}; a[7] = {7.0, 0.0}; a[8] = {14.0, 0.0};

    // Initialize b matrix
    for (int i = 0; i < N * N; ++i) {
        b[i] = a[i] + Complex(0, a[i].real); // a+1j*a
    }

    // Initialize A and X matrices
    for (int i = 0; i < N; ++i) {
        for (int j = 0; j < N; ++j) {
            A[i*N+j]=b[i*N+j]+conj(b[j*N+i]);
	}
        X[i*N+i] = Complex(1.0, 0.0); // Identity matrix
    }
    
    delete[] a;
    delete[] b;
    //printArray(A,N,N);
    //printArray(X,N,N);
    //printf("===\n");
    for (int iter = 0; iter < 10; ++iter) {
        float max_A = 0.0;
        int p = 0, q = 0;

        // Find p, q
        for (int i = 0; i < N; ++i) {
            for (int j = i + 1; j < N; ++j) {
                if (i != j) {
                    float a = A[i * N + j].real * A[i * N + j].real + A[i * N + j].imag * A[i * N + j].imag;
                    if (a > max_A) {
                        max_A = a;
                        p = i;
                        q = j;
                    }
                }
            }
        }

        printf("max:%f (%d,%d)\n", max_A , p ,q);
	if(max_A<1.0e-10){
            break;
	}
        float theta = M_PI / 4.0 - 0.5 * std::atan2((A[p * N + p].real - A[q * N + q].real), 2.0 * std::sqrt(max_A));
        float phi = -std::atan2(A[p * N + q].imag, A[p * N + q].real);
        printf("%f %f\n", theta, phi);

        float cos_theta = std::cos(theta);
        float sin_theta = std::sin(theta);
        Complex e_neg_phi(std::cos(phi),-std::sin(phi));
        Complex e_phi(std::cos(phi),std::sin(phi));

        // Update An matrix
        for (int i = 0; i < N; ++i) {
            for (int j = i; j < N; ++j) {
                if (i == p && j != q && j != p) {
                    An[p * N + j] = A[p * N + j] * cos_theta + A[q * N + j] * sin_theta * e_neg_phi;
                    An[j * N + p] = conj(An[p * N + j]);
                } else if (i != q && i != p && j == p) {
                    An[p * N + i] = A[p * N + i] * cos_theta + A[q * N + i] * sin_theta * e_neg_phi;
                    An[i * N + p] = conj(An[p * N + i]);

                } else if (i == q && j != q) {
                    An[q * N + j] = A[p * N + j] * sin_theta*e_phi - A[q * N + j] * cos_theta;
                    An[j * N + q] = conj(An[q * N + j]);
                } else if (i != q && i != p && j == q) {
                    An[q * N + i] = A[p * N + i] * sin_theta*e_phi - A[q * N + i] * cos_theta;
                    An[i * N + q] = conj(An[q * N + i]);

                } else if (i == p && j == p) {
                    An[p * N + p] = A[p * N + p] * cos_theta * cos_theta + A[q * N + q] * sin_theta * sin_theta +
                                    A[p * N + q] * sin_theta * cos_theta * e_phi + A[q * N + p] * sin_theta * cos_theta *e_neg_phi;
                } else if (i == q && j == q) {
                    An[q * N + q] = A[p * N + p] * sin_theta * sin_theta + A[q * N + q] * cos_theta * cos_theta -
                                    A[p * N + q] * sin_theta * cos_theta * e_phi - A[q * N + p] * sin_theta * cos_theta * e_neg_phi;
                } else if (i == p && j == q) {
                    An[p * N + q] = A[p * N + p] * sin_theta * cos_theta - A[q * N + q] * sin_theta * cos_theta -
                                    A[p * N + q] * cos_theta * cos_theta * e_phi + A[q * N + p] * sin_theta * sin_theta * e_neg_phi;
                    An[p * N + q] = e_neg_phi*An[p * N + q];
                    An[q * N + p] = conj(An[p * N + q]);
                } else {
                    An[i * N + j] = A[i * N + j];
                    An[j * N + i] = A[j * N + i];
                }
            }
        }

        // Update X matrix
        for (int i = 0; i < N; ++i) {
            Complex bip = X[i * N + p];
            Complex biq = X[i * N + q];
            X[i * N + p] = bip * cos_theta + biq * sin_theta*e_phi;
            X[i * N + q] = bip * sin_theta*e_neg_phi - biq * cos_theta;
        }

        // Display An matrix
        printArray(An, N, N);
	memcpy(A,An,sizeof(Complex)*N*N);
        // Display diagonal elements of An
        for (int i = 0; i < N; ++i) {
            printf(">> %f %f\n", An[i * N + i].real,An[i * N + i].imag);
	}
    }

    printArray(X,N,N);
    delete[] A;
    delete[] X;
    delete[] An;

    return 0;
}

void setup() {
  // put your setup code here, to run once:
  Serial.begin(115200);
  while (!Serial);

  Serial.println("Init Audio Library");
  
  run();
}

void loop() {
  // put your main code here, to run repeatedly:

}
