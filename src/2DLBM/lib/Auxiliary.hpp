#ifndef AUXILIARY_HPP
#define AUXILIARY_HPP

#include <vector>
#include <stdexcept>
#include <initializer_list>

#include <numeric> // per std::inner_product

template <typename T>
T dotProduct(const std::vector<T>& a, const std::vector<T>& b)
{
    if (a.size() != b.size()) {
        throw std::invalid_argument("Vector sizes must match.");
    }
    return std::inner_product(a.begin(), a.end(), b.begin(), T{0});
}


/**
 * @brief A generic N-dimensional matrix class with simple (i, j, k...) element access.
 * 
 * Internally stored as a 1D vector for performance, but accessed using intuitive multi-index syntax.
 */
template <class T>
class Matrix {
private:
    std::vector<T> data;             // Flat storage
    std::vector<int> dimensions;     // Sizes of each dimension
    std::vector<int> strides;        // Strides for flat index calculation

    // Compute strides for fast index flattening
    void computeStrides() {
        strides.resize(dimensions.size());
        int stride = 1;
        for (int i = dimensions.size() - 1; i >= 0; --i) {
            strides[i] = stride;
            stride *= dimensions[i];
        }
    }

    // Flatten multi-index to flat index
    int getFlatIndex(std::initializer_list<int> indices) const {
        if (indices.size() != dimensions.size()) {
            throw std::runtime_error("Incorrect number of indices");
        }

        int flatIndex = 0;
        int i = 0;
        for (int idx : indices) {
            if (idx < 0 || idx >= dimensions[i]) {
                throw std::runtime_error("Index out of bounds");
            }
            flatIndex += idx * strides[i];
            ++i;
        }
        return flatIndex;
    }

public:
    // Constructor with dimension sizes
    Matrix(const std::vector<int>& dims) : dimensions(dims) {
        int totalSize = 1;
        for (int d : dims) totalSize *= d;
        data.resize(totalSize);
        computeStrides();
    }

    // Default constructor
    Matrix() = default;

    // Access operator for 2D matrix: mat(i, j)
    T& operator()(int i, int j) {
        return data.at(i * strides[0] + j * strides[1]);
    }

    const T& operator()(int i, int j) const {
        return data.at(i * strides[0] + j * strides[1]);
    }

    // Getters
    const std::vector<int>& shape() const { return dimensions; }
    int totalSize() const { return data.size(); }

    // Set value at (i, j, ...)
    void set(std::initializer_list<int> indices, const T& value) {
        data.at(getFlatIndex(indices)) = value;
    }

    // Get flat index of a given (i, j, ...) coordinate
    int flatIndex(std::initializer_list<int> indices) const {
        return getFlatIndex(indices);
    }

    T getCopy(int i, int j) const {
        return data.at(i * strides[0] + j * strides[1]);
    }
};

std::vector <int> evaluateBoundary(const std::vector<int>& indices, const Matrix<bool> &obstacleMatrix)
{
    const int &col = indices.at(0);
    const int &row = indices.at(1);
    const int &NX = obstacleMatrix.shape().at(0);
    const int &NY = obstacleMatrix.shape().at(1);
    std::vector<int> boundary_here(4);

    // Evaluate boundary array, considering near objects
    // Horizontal 
    if( col>0 && obstacleMatrix(col-1,row))
    {
        boundary_here.at(0) = -1;
    }
    else if(col<(NX-1) && obstacleMatrix(col+1,row)) 
    {
        boundary_here.at(0) = 1;
    }
    else 
    {
        boundary_here.at(0) = 0;
    }

    // Vertical 
    if (row>0 && obstacleMatrix(col,row-1))
    {
        boundary_here.at(1) = -1;
    }
    else if (row<(NY-1) && obstacleMatrix(col,row+1))
    {
        boundary_here.at(1) = 1;
    }
    else 
    {
        boundary_here.at(1) = 0;
    }

    // First Diag
    if (col>0 && row>0 && obstacleMatrix(col-1,row-1))
    {
        boundary_here.at(2) = -1;
    }
    else if (col<(NX-1) && row<(NY-1) && obstacleMatrix(col+1, row+1))
    {
        boundary_here.at(2) = 1;
    }
    else 
    {
        boundary_here.at(2) = 0;
    }

    // Second Diag
    if (col>0 && row<(NY-1) && obstacleMatrix(col-1,row+1))
    {
        boundary_here.at(3) = 1;
    }
    else if (col<(NX-1) && row>0 && obstacleMatrix(col+1, row-1))
    {
        boundary_here.at(3) = -1;
    }
    else 
    {
        boundary_here.at(3) = 0;
    }  
    return boundary_here;
}

#endif // AUXILIARY_HPP

/*
EXAMPLE 

Matrix<double> mat({100, 50});  // 2D Matrix

mat(0, 0) = 1.0;
double x = mat(10, 20);

Matrix<double> cube({10, 10, 10});
cube(1, 2, 3) = 5.5;

*/