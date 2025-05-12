#ifndef AUXILIARY_CPP
#define AUXILIARY_CPP

#include <vector>
#include <stdexcept>
#include <initializer_list>

#include <numeric>

template <typename T>
T dotProduct(const std::vector<T>& a, const std::vector<T>& b)
{
    if (a.size() != b.size()) {
        throw std::invalid_argument("Vector sizes must match.");
    }
    return std::inner_product(a.begin(), a.end(), b.begin(), T{0});
}

template <class T>
class Matrix
{
private:
    std::vector<T> data;
    std::vector<int> dimensions;

    int calculateFlatIndex(const std::vector<int> &indices) const
    {
        int flatIndex = 0;
        for (int i = indices.size() - 1; i >= 0; --i)
        {
            flatIndex = flatIndex * dimensions[i] + indices[i];
        }

        return flatIndex;
    }

public:
    Matrix(const std::vector<int> &dims)
    {
        dimensions = dims;
        int totalSize = 1;
        for (int dim : dims)
        {
            totalSize *= dim;
        }
        data.resize(totalSize);
    }

    Matrix() = default;

    T &operator()(const std::vector<int> &indices)
    {
        return data.at(calculateFlatIndex(indices));
    }

    const T operator()(int x, int y) const
    {
        return data.at(x + y * dimensions[0]);
    }

    T &operator()(int x, int y)
    {
        return data.at(x + y * dimensions[0]);
    }

    T &operator()(int x, int y, int z)
    {
        return data.at(x + y * dimensions[0] + z * dimensions[0] * dimensions[1]);
    }

    const std::vector<int> &shape() const
    {
        return dimensions;
    }

    int size() const
    {
        return data.size();
    }

    const std::vector<int> indicesAt(int flatIndex) const
    {
        if (flatIndex >= static_cast<int>(data.size()))
        {
            throw std::runtime_error("Invalid flat index");
        }

        std::vector<int> indices(dimensions.size());
        for (int i = 0; i < static_cast<int>(dimensions.size()); ++i)
        {
            indices[i] = flatIndex % dimensions[i];
            flatIndex /= dimensions[i];
        }
        return indices;
    }

    T &at(int flatIndex)
    {
        return const_cast<T &>(static_cast<const Matrix *>(this)->at(flatIndex));
    }

    const T &at(int flatIndex) const
    {
        return data.at(flatIndex);
    }

    void set(const std::vector<int> &indices, const T &value)
    {
        data.at(calculateFlatIndex(indices)) = value;
    }

    void setAt(int flatIndex, const T &value)
    {
        data.at(flatIndex) = value;
    }

    T getCopy(int x, int y) const
    {
        return data.at(x + y * dimensions[0]);
    }
};

#endif // AUXILIARY_HPP

/*
EXAMPLE 

Matrix<double> mat({100, 50});  // 2D Matrix

mat(0, 0) = 1.0;
double x = mat(10, 20);

Matrix<double> cube({10, 10, 10});
cube(1, 2, 3) = 5.5;

*/