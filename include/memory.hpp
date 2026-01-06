#ifndef MEMORY_HPP
#define MEMORY_HPP

#include <cstdlib>
#include <cstddef>
#include <new>

/**
 * Aligned memory allocation utilities for HPC.
 * Uses 64-byte alignment for AVX-512 compatibility.
 */

constexpr size_t CACHE_LINE_SIZE = 64;  // AVX-512 alignment

/**
 * Allocate aligned memory.
 * @param count Number of elements to allocate
 * @param alignment Memory alignment (default: 64 bytes for AVX-512)
 * @return Pointer to aligned memory, or nullptr on failure
 */
template<typename T>
T* aligned_alloc_array(size_t count, size_t alignment = CACHE_LINE_SIZE) {
    if (count == 0) return nullptr;
    
    size_t size = count * sizeof(T);
    
    // Ensure size is a multiple of alignment for aligned_alloc
    size_t aligned_size = ((size + alignment - 1) / alignment) * alignment;
    
    void* ptr = std::aligned_alloc(alignment, aligned_size);
    if (!ptr) {
        throw std::bad_alloc();
    }
    
    return static_cast<T*>(ptr);
}

/**
 * Free aligned memory.
 * @param ptr Pointer to memory allocated with aligned_alloc_array
 */
template<typename T>
void aligned_free_array(T* ptr) {
    std::free(ptr);
}

/**
 * RAII wrapper for aligned arrays.
 * Automatically frees memory on destruction.
 */
template<typename T>
class AlignedArray {
public:
    AlignedArray() : data_(nullptr), size_(0) {}
    
    explicit AlignedArray(size_t count) 
        : data_(aligned_alloc_array<T>(count)), size_(count) {
        // Zero-initialize
        for (size_t i = 0; i < size_; ++i) {
            data_[i] = T{};
        }
    }
    
    ~AlignedArray() {
        if (data_) {
            aligned_free_array(data_);
        }
    }
    
    // Move constructor
    AlignedArray(AlignedArray&& other) noexcept 
        : data_(other.data_), size_(other.size_) {
        other.data_ = nullptr;
        other.size_ = 0;
    }
    
    // Move assignment
    AlignedArray& operator=(AlignedArray&& other) noexcept {
        if (this != &other) {
            if (data_) {
                aligned_free_array(data_);
            }
            data_ = other.data_;
            size_ = other.size_;
            other.data_ = nullptr;
            other.size_ = 0;
        }
        return *this;
    }
    
    // No copy
    AlignedArray(const AlignedArray&) = delete;
    AlignedArray& operator=(const AlignedArray&) = delete;
    
    // Access
    T& operator[](size_t index) { return data_[index]; }
    const T& operator[](size_t index) const { return data_[index]; }
    
    T* data() { return data_; }
    const T* data() const { return data_; }
    
    size_t size() const { return size_; }
    
    T* begin() { return data_; }
    T* end() { return data_ + size_; }
    const T* begin() const { return data_; }
    const T* end() const { return data_ + size_; }

private:
    T* data_;
    size_t size_;
};

#endif // MEMORY_HPP

