#include <cstdio>
#include <cmath>

struct Complex {
    double real;
    double imag;

    Complex(double r, double i) : real(r), imag(i) {}

    Complex operator+(const Complex& other) const {
        return Complex(real + other.real, imag + other.imag);
    }

    Complex operator-(const Complex& other) const {
        return Complex(real - other.real, imag - other.imag);
    }

    Complex operator*(const Complex& other) const {
        return Complex(real * other.real - imag * other.imag, real * other.imag + imag * other.real);
    }

    Complex operator/(const double& scalar) const {
        return Complex(real / scalar, imag / scalar);
    }
};

// Function to print a 2D array
void printArray(Complex* arr, int rows, int cols) {
    for (int i = 0; i < rows; ++i) {
        for (int j = 0; j < cols; ++j) {
            printf("%f + %fi ", arr[i * cols + j].real, arr[i * cols + j].imag);
        }
        printf("\n");
    }
}

int main() {
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
        b[i] = a[i] + Complex(a[i].real, -a[i].imag); // conjugate of a[i]
    }

    // Initialize A and X matrices
    for (int i = 0; i < N * N; ++i) {
        A[i] = b[i] + Complex(b[i].real, -b[i].imag); // conjugate of b[i]
        X[i] = (i % (N + 1) == 0) ? Complex(1.0, 0.0) : Complex(0.0, 0.0); // Identity matrix
    }

    for (int iter = 0; iter < 10; ++iter) {
        double max_A = 0.0;
        int p = 0, q = 0;

        // Find p, q
        for (int i = 0; i < N; ++i) {
            for (int j = i + 1; j < N; ++j) {
                if (i != j) {
                    double a = std::sqrt(A[i * N + j].real * A[i * N + j].real + A[i * N + j].imag * A[i * N + j].imag);
                    if (a > max_A) {
                        max_A = a;
                        p = i;
                        q = j;
                    }
                }
            }
        }

        printf("%f\n", max_A);

        double theta = M_PI / 4.0 - 0.5 * std::atan2((A[p * N + p].real - A[q * N + q].real), 2.0 * max_A);
        double phi = -std::atan2(A[p * N + q].imag, A[p * N + q].real);
        printf("%f %f\n", theta, phi);

        double cos_theta = std::cos(theta);
        double sin_theta = std::sin(theta);

        // Update An matrix
        for (int i = 0; i < N; ++i) {
            for (int j = i; j < N; ++j) {
                if (i == p && j != q && j != p) {
                    An[p * N + j] = A[p * N + j] * cos_theta + A[q * N + j] * sin_theta;
                    An[j * N + p] = An[p * N + j];
                } else if (i != q && i != p && j == p) {
                    An[p * N + i] = A[p * N + i] * cos_theta + A[q * N + i] * sin_theta;
                    An[i * N + p] = An[p * N + i];
                } else if (i == q && j != q) {
                    An[q * N + j] = A[p * N + j] * sin_theta - A[q * N + j] * cos_theta;
                    An[j * N + q] = An[q * N + j];
                } else if (i != q && i != p && j == q) {
                    An[q * N + i] = A[p * N + i] * sin_theta - A[q * N + i] * cos_theta;
                    An[i * N + q] = An[q * N + i];
                } else if (i == p && j == p) {
                    An[p * N + p] = A[p * N + p] * cos_theta * cos_theta + A[q * N + q] * sin_theta * sin_theta +
                                    A[p * N + q] * sin_theta * cos_theta + A[q * N + p] * sin_theta * cos_theta;
                } else if (i == q && j == q) {
                    An[q * N + q] = A[p * N + p] * sin_theta * sin_theta + A[q * N + q] * cos_theta * cos_theta -
                                    A[p * N + q] * sin_theta * cos_theta - A[q * N + p] * sin_theta * cos_theta;
                } else if (i == p && j == q) {
                    An[p * N + q] = A[p * N + p] * sin_theta * cos_theta - A[q * N + q] * sin_theta * cos_theta -
                                    A[p * N + q] * cos_theta * cos_theta - A[q * N + p] * sin_theta * sin_theta;
                    An[q * N + p] = An[p * N + q];
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
            X[i * N + p] = bip * cos_theta + biq * sin_theta;
            X[i * N + q] = bip * sin_theta - biq * cos_theta;
        }

        // Display An matrix
        printArray(An, N, N);

        // Compute C matrix
        for (int i = 0; i < N; ++i) {
            for (int j = 0; j < N; ++j) {
                Complex sum = {0.0, 0.0};
                for (int k = 0; k < N; ++k) {
                    sum = sum + Complex(An[k * N + i].real * An[k * N + j].real - An[k * N + i].imag * An[k * N + j].imag,
                                        An[k * N + i].real * An[k * N + j].imag + An[k * N + i].imag * An[k * N + j].real);
                }
                C[i * N + j] = sum;
            }
        }

        // Threshold small values in C
        for (int i = 0; i < N * N; ++i) {
            if (std::sqrt(C[i].real * C[i].real + C[i].imag * C[i].imag) < 1.0e-5) {
                C[i] = {0.0, 0.0};
            }
        }

        // Display C matrix
        printArray(C, N, N);
        delete[] A;
        A = An;

        // Display diagonal elements of An
        printf(">> %f %f\n", An[p * N + p].real, An[q * N + q].real);
    }

    delete[] a;
    delete[] b;
    delete[] A;
    delete[] X;
    delete[] An;
    delete[] C;

    return 0;
}
