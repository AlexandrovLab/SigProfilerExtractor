/**
 * @file cuda_nmf_single.cu
 *
 * @brief Run NMF's multiplicative update with CUDA
 *
 * @authors Mark Barnes, Jason Dai, Marcos Díaz-Gay, Johannes Blaschke
 * Contact: mdbarnes@ucsd.edu
 *
 */
#include <cfloat>
#include <string>
#include <cmath>
#include <iostream>
#include <sstream>
#include <cuda.h>
#include <cuda_runtime.h>
#include <cuda_runtime_api.h>
#include <cublas_v2.h>
#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <pybind11/stl.h>


namespace py = pybind11;

typedef float dtype;

// number of blocks per 1 dimension of the grid
const int block_size = 32;
const int ew_block_size = 128;

// number of total blocks in the grid
const int blockSize = block_size * block_size;



__global__ void sumCommSingleBlock(const dtype *a, double *out, const int arraySize) {
    int idx = threadIdx.x;

    double sum = 0;
    for (int i = idx; i < arraySize; i += blockSize)
        sum += a[i];
    __shared__ double r[blockSize];
    r[idx] = sum;
    __syncthreads();
    for (int size = blockSize/2; size>0; size/=2) { //uniform
        if (idx<size)
            r[idx] += r[idx+size];
        __syncthreads();
    }
    if (idx == 0)
        *out = r[0];
}


void sum_matrix_elements_launch(
        cudaStream_t * stream,
        const dtype * d_C,
        double * d_sum,
        int nr_rows, int nr_cols,
        int gpu_id
    ) {

    cudaSetDevice(gpu_id);
    int arraySize = nr_rows * nr_cols;
    sumCommSingleBlock<<<1, blockSize, 0, * stream>>>(d_C, d_sum, arraySize);
}


// matrix logarithm element wise (naive) kernel: log(C) (element wise)
__global__ void ew_log_kernel(dtype *C, int dx, int dy) {

    int idx = threadIdx.x+blockDim.x*blockIdx.x; // create thread x index
    int idy = threadIdx.y+blockDim.y*blockIdx.y; // create thread y index

    if ((idx < dx) && (idy < dy)){
        C[idx*dy+idy] = log(C[idx*dy+idy]);
    }
}

void ew_log_launch(
        cudaStream_t * stream,
        dtype * d_C,
        int nr_rows, int nr_cols,
        int gpu_id
    ) {

    cudaSetDevice(gpu_id);
    dim3 block(block_size, block_size);  // dim3 variable holds 3 dimensions
    dim3 grid((nr_rows+block.x-1)/block.x, (nr_cols+block.y-1)/block.y);
    ew_log_kernel<<<grid, block, 0, * stream>>>(d_C, nr_rows, nr_cols);
}


// matrix multiply element wise (naive) kernel: C = A * B (element wise)
__global__ void ew_multiplication_kernel(const dtype *A, const dtype *B, dtype *C, int dx, int dy) {

    int idx = threadIdx.x+blockDim.x*blockIdx.x; // create thread x index
    int idy = threadIdx.y+blockDim.y*blockIdx.y; // create thread y index

    if ((idx < dx) && (idy < dy)){
        C[idx*dy+idy] = A[idx*dy+idy] * B[idx*dy+idy];
    }
}


void ew_multiplication_launch(
        cudaStream_t * stream,
        const dtype * d_A, const dtype * d_B, dtype * d_C,
        int nr_rows, int nr_cols,
        int gpu_id
    ) {

    cudaSetDevice(gpu_id);
    dim3 block(block_size, block_size);  // dim3 variable holds 3 dimensions
    dim3 grid((nr_rows+block.x-1)/block.x, (nr_cols+block.y-1)/block.y);
    ew_multiplication_kernel<<<grid, block, 0, * stream>>>(d_A, d_B, d_C, nr_rows, nr_cols);
    // cudaDeviceSynchronize();

}


// element wise matrix division (naive) kernel: C = A / B (element wise)
__global__ void ew_division_kernel(const dtype *A, const dtype *B, dtype *C, int dx, int dy) {
    int threadId = blockIdx.x*blockDim.x + threadIdx.x;
    if(threadId<(dx*dy)) {
        C[threadId] = A[threadId] / B[threadId];
    }
}


void ew_division_launch(
        cudaStream_t * stream, 
        const dtype * d_A, const dtype * d_B, dtype * d_C,
        int nr_rows, int nr_cols,
        int gpu_id
    ) {

    cudaSetDevice(gpu_id);
    dim3 block(ew_block_size);
    dim3 grid((nr_rows*nr_cols-1)/ew_block_size+1);
    ew_division_kernel<<<grid, block, 0, * stream>>>(d_A, d_B, d_C, nr_rows, nr_cols);
    // cudaDeviceSynchronize();
}


// Multiply matrices  A with tranpose  B on GPU and save the result in C
// C(m,n) = A(m,k) * B^T(k,n)
// C-array style based on: https://stackoverflow.com/questions/56043539/cublassgemm-row-major-multiplication
void gpu_blas_mmul_rightT(
        cublasHandle_t * handle,
        dtype * C, const dtype * A, const dtype * B,
        const int m, const int k, const int n, int gpu_id
    ) {

    // m is rowsA, k is colsA, n is colsB
    // lda=k <= swap row <-> column as fortran order is expected
    int lda=k, ldb=k, ldc=n;
    const dtype alpha = 1;
    const dtype beta  = 0;

    cudaSetDevice(gpu_id);

    // cuBLAS single-precision matrix-matrix multiply
    // http://docs.nvidia.com/cuda/cublas/index.html#cublas-lt-t-gt-gemm
    cublasSgemm(* handle,
            CUBLAS_OP_T, CUBLAS_OP_N,
            n, m, k,
            & alpha,
            B, ldb, // expects fortran array-oder => swap A and B
            A, lda,
            & beta,
            C, ldc);
}


// Multiply matrices transpose A with B on GPU and save the result in C
// C(m,n) = A^T(m,k) * B(k,n)
// C-array style based on: https://stackoverflow.com/questions/56043539/cublassgemm-row-major-multiplication
void gpu_blas_mmul_leftT(
        cublasHandle_t * handle,
        dtype * C, const dtype * A, const dtype * B,
        const int m, const int k, const int n, int gpu_id
    ) {
    // m is rowsA, k is colsA, n is colsB
    // lda=k <= swap row <-> column as fortran order is expected
    int lda=m, ldb=n, ldc=n;
    const dtype alpha = 1;
    const dtype beta  = 0;

    cudaSetDevice(gpu_id);

    // cuBLAS single-precision matrix-matrix multiply
    // http://docs.nvidia.com/cuda/cublas/index.html#cublas-lt-t-gt-gemm
    cublasSgemm(* handle,
            CUBLAS_OP_N, CUBLAS_OP_T,
            n, m, k,
            & alpha,
            B, ldb, // expects fortran array-oder => swap A and B
            A, lda,
            & beta,
            C, ldc);
    


}


// Multiply the arrays A and B on GPU and save the result in C
// C(m,n) = A(m,k) * B(k,n)
// C-array style based on: https://stackoverflow.com/questions/56043539/cublassgemm-row-major-multiplication
void gpu_blas_mmul(
        cublasHandle_t * handle,
        dtype * C, const dtype * A, const dtype * B, 
        const int m, const int k, const int n, int gpu_id
    ) {

    // lda=k <= swap row <-> column as fortran order is expected
    int lda=k, ldb=n, ldc=n;
    const dtype alpha = 1;
    const dtype beta  = 0;

    cudaSetDevice(gpu_id);

    // cuBLAS single-precision matrix-matrix multiply
    // http://docs.nvidia.com/cuda/cublas/index.html#cublas-lt-t-gt-gemm
    cublasSgemm(* handle,
            CUBLAS_OP_N, CUBLAS_OP_N,
            n, m, k,
            & alpha,
            B, ldb, // expects fortran array-oder => swap A and B
            A, lda,
            & beta,
            C, ldc);



}


typedef struct NMFHandle {

    dtype * d_V, * d_W, * d_H, * d_ones;

    // d_wh = W@H
    dtype * d_wh;
    // d_vwh = V/(W@H)
    dtype * d_vwh;
    // d_numerator_w
    dtype * d_numerator_w;
    // d_denominator_w
    dtype * d_denominator_w;
    // d_NDW = d_numerator_w/d_denominator_w
    dtype * d_NDW;
    // d_numerator_h
    dtype * d_numerator_h;
    // d_denominator_h
    dtype * d_denominator_h;
    // d_NDH = d_numerator_w/d_denominator_w
    dtype * d_NDH;
    // d_sum
    double * d_sum;
    // d_Vsum
    double * d_Vsum;
    // d_WHsum
    double * d_WHsum;

    cublasHandle_t handle_a;
    cublasHandle_t handle_b;
    cublasHandle_t handle_c;
    cublasHandle_t handle_d;
    cublasHandle_t handle_e;
    cublasHandle_t handle_f;
    cublasHandle_t handle_g;
    cudaStream_t stream;
    cudaGraph_t graph;
    cudaGraph_t kl_loss_graph;
    cudaGraphExec_t instance;
    cudaGraphExec_t kl_loss_instance;

    bool graph_created = false;
} NMFHandle;



template <class T> class ptr_wrapper {
    public:
        ptr_wrapper() : ptr(nullptr) {}
        ptr_wrapper(T * ptr) : ptr(ptr) {}
        ptr_wrapper(const ptr_wrapper& other) : ptr(other.ptr) {}
        T & operator* () const { return * ptr; }
        T * operator->() const { return   ptr; }
        T * get() const { return ptr; }
        void destroy() { delete ptr; }
        T& operator[](std::size_t idx) const { return ptr[idx]; }
    private:
        T* ptr;
};



ptr_wrapper<NMFHandle> NMFHandleNew() {
    return new NMFHandle;
}



void NMFHandleCreate(
        py::array_t<dtype> V, 
        py::array_t<dtype> W,
        py::array_t<dtype> H,
        py::array_t<dtype> ones,
        ptr_wrapper<NMFHandle> nmf_handle,
        int gpu_id
    ) {

    // Define shape of objects
    int nr_rows_V = V.shape(0);
    int nr_cols_V = V.shape(1);
    int nr_rows_W = W.shape(0);
    int nr_cols_W = W.shape(1);
    int nr_rows_H = H.shape(0);
    int nr_cols_H = H.shape(1);
    int nr_rows_ones = ones.shape(0);
    int nr_cols_ones = ones.shape(1);

    cudaSetDevice(gpu_id);
    cudaStreamCreate(& nmf_handle->stream);

    // Declare and allocate device variables
    cudaMalloc((void **) & nmf_handle->d_V, nr_rows_V*nr_cols_V*sizeof(dtype));
    cudaMalloc((void **) & nmf_handle->d_W, nr_rows_W*nr_cols_W*sizeof(dtype));
    cudaMalloc((void **) & nmf_handle->d_H, nr_rows_H*nr_cols_H*sizeof(dtype));
    cudaMalloc((void **) & nmf_handle->d_ones, nr_rows_ones*nr_cols_ones*sizeof(dtype));

    //__________________________________________________________________________
    // Local temporary storate -- this is free'd at the end of this function
    //

    // These can be re-used (they just need to be re-computed between W and H)
    // d_wh = W@H
    cudaMalloc((void **) & nmf_handle->d_wh, nr_rows_W*nr_cols_H*sizeof(dtype));
    // d_vwh = V/(W@H)
    cudaMalloc((void **) & nmf_handle->d_vwh, nr_rows_V*nr_cols_V*sizeof(dtype));


    // Local storage for the W update
    // d_numerator_w
    cudaMalloc((void **) & nmf_handle->d_numerator_w, nr_rows_W*nr_cols_W*sizeof(dtype));
    // d_denominator_w
    cudaMalloc((void **) & nmf_handle->d_denominator_w, nr_rows_ones*nr_cols_H*sizeof(dtype));
    // d_NDW = d_numerator_w/d_denominator_w
    cudaMalloc((void **) & nmf_handle->d_NDW, nr_rows_W*nr_cols_W*sizeof(dtype));

    // Local storage for the H update
    // d_numerator_h
    cudaMalloc((void **) & nmf_handle->d_numerator_h, nr_rows_H*nr_cols_H*sizeof(dtype));
    // d_denominator_h
    cudaMalloc((void **) & nmf_handle->d_denominator_h, nr_rows_H*nr_cols_H*sizeof(dtype));
    // d_NDH = d_numerator_w/d_denominator_w
    cudaMalloc((void **) & nmf_handle->d_NDH, nr_rows_H*nr_cols_H*sizeof(dtype));

    // Local storage for dtypes
    // d_sum
    cudaMalloc((void **) & nmf_handle->d_sum, sizeof(double));
    // d_Vsum
    cudaMalloc((void **) & nmf_handle->d_Vsum, sizeof(double));
    // d_WHsum
    cudaMalloc((void **) & nmf_handle->d_WHsum, sizeof(double));

    //--------------------------------------------------------------------------

    cublasCreate(& nmf_handle->handle_a);
    cublasSetStream(nmf_handle->handle_a, nmf_handle->stream);
    cublasCreate(& nmf_handle->handle_b);
    cublasSetStream(nmf_handle->handle_b, nmf_handle->stream);
    cublasCreate(& nmf_handle->handle_c);
    cublasSetStream(nmf_handle->handle_c, nmf_handle->stream);
    cublasCreate(& nmf_handle->handle_d);
    cublasSetStream(nmf_handle->handle_d, nmf_handle->stream);
    cublasCreate(& nmf_handle->handle_e);
    cublasSetStream(nmf_handle->handle_e, nmf_handle->stream);
    cublasCreate(& nmf_handle->handle_f);
    cublasSetStream(nmf_handle->handle_f, nmf_handle->stream);
    cublasCreate(& nmf_handle->handle_g);
    cublasSetStream(nmf_handle->handle_g, nmf_handle->stream);

    nmf_handle->graph_created = false;
}



void NMFHandleDestroy(
        ptr_wrapper<NMFHandle> nmf_handle
    ) {

    cublasDestroy(nmf_handle->handle_a);
    cublasDestroy(nmf_handle->handle_b);
    cublasDestroy(nmf_handle->handle_c);
    cublasDestroy(nmf_handle->handle_d);
    cublasDestroy(nmf_handle->handle_e);
    cublasDestroy(nmf_handle->handle_f);
    cublasDestroy(nmf_handle->handle_g);

    // Free memory
    cudaFree(nmf_handle->d_V);
    cudaFree(nmf_handle->d_W);
    cudaFree(nmf_handle->d_H);
    cudaFree(nmf_handle->d_ones);
    cudaFree(nmf_handle->d_wh);
    cudaFree(nmf_handle->d_vwh);
    cudaFree(nmf_handle->d_numerator_w);
    cudaFree(nmf_handle->d_denominator_w);
    cudaFree(nmf_handle->d_NDW);
    cudaFree(nmf_handle->d_numerator_h);
    cudaFree(nmf_handle->d_denominator_h);
    cudaFree(nmf_handle->d_NDH);
    cudaFree(nmf_handle->d_sum);
    cudaFree(nmf_handle->d_Vsum);
    cudaFree(nmf_handle->d_WHsum);

    cudaStreamDestroy(nmf_handle->stream);

    delete nmf_handle.get();
}



int fit_loop(
        py::array_t<dtype> V, 
        py::array_t<dtype> W,
        py::array_t<dtype> H,
        py::array_t<dtype> ones,
        ptr_wrapper<NMFHandle> nmf_handle,
        int test_conv,
        int max_iterations,
        int min_iterations,
        dtype tolerance,
        int gpu_id
    ) {

    cudaSetDevice(gpu_id);
    // Define shape of objects
    int nr_rows_V = V.shape(0);
    int nr_cols_V = V.shape(1);
    int nr_rows_W = W.shape(0);
    int nr_cols_W = W.shape(1);
    int nr_rows_H = H.shape(0);
    int nr_cols_H = H.shape(1);
    int nr_rows_ones = ones.shape(0);
    int nr_cols_ones = ones.shape(1);
    double loss_init = 0;
    double kl_loss = 0;
    double sum = 0;
    double Vsum = 0;
    double WHsum = 0;
    double prev_loss = 0;
    bool loss_converged = false;


    // Get numpy buffers
    py::buffer_info V_buf = V.request();
    py::buffer_info W_buf = W.request();
    py::buffer_info H_buf = H.request();
    py::buffer_info ones_buf = ones.request();

    // Obtain numpy data pointer
    dtype * V_ptr = reinterpret_cast<dtype*>(V_buf.ptr);
    dtype * W_ptr = reinterpret_cast<dtype*>(W_buf.ptr);
    dtype * H_ptr = reinterpret_cast<dtype*>(H_buf.ptr);
    dtype * ones_ptr = reinterpret_cast<dtype*>(ones_buf.ptr);


    // Copy host to device
    cudaMemcpy(
        nmf_handle->d_V, V_ptr, nr_rows_V*nr_cols_V*sizeof(dtype), cudaMemcpyHostToDevice
    );
    cudaMemcpy(
        nmf_handle->d_W, W_ptr, nr_rows_W*nr_cols_W*sizeof(dtype), cudaMemcpyHostToDevice
    );
    cudaMemcpy(
        nmf_handle->d_H, H_ptr, nr_rows_H*nr_cols_H*sizeof(dtype), cudaMemcpyHostToDevice
    );
    cudaMemcpy(
        nmf_handle->d_ones, ones_ptr, nr_rows_ones*nr_cols_ones*sizeof(dtype), cudaMemcpyHostToDevice
    );
    
    


    for(int iter=0; iter < max_iterations / test_conv; iter++) {
        for(int i=0; i < test_conv; i++) {
            if (! nmf_handle->graph_created) {
                cudaStreamBeginCapture(nmf_handle->stream, cudaStreamCaptureModeGlobal);

                //__________________________________________________________________________
                // W - Update
                //

                // w @ h
                gpu_blas_mmul(& nmf_handle->handle_c, nmf_handle->d_wh, nmf_handle->d_W, nmf_handle->d_H, nr_rows_W, nr_cols_W, nr_cols_H, gpu_id);
                // v / (w@h)
                ew_division_launch(& nmf_handle->stream, nmf_handle->d_V, nmf_handle->d_wh, nmf_handle->d_vwh, nr_rows_V, nr_cols_V, gpu_id);
                // numerator_w = v / (w@h) @ transpose(h)
                gpu_blas_mmul_rightT(& nmf_handle->handle_a, nmf_handle->d_numerator_w, nmf_handle->d_vwh, nmf_handle->d_H, nr_rows_V, nr_cols_V, nr_rows_H, gpu_id);
                // denominator_w = ones @ transpose(h)
                gpu_blas_mmul_rightT(& nmf_handle->handle_b, nmf_handle->d_denominator_w, nmf_handle->d_ones, nmf_handle->d_H, nr_rows_ones, nr_cols_ones, nr_rows_H, gpu_id);
                // numerator_w / denominator_w
                ew_division_launch(& nmf_handle->stream, nmf_handle->d_numerator_w, nmf_handle->d_denominator_w, nmf_handle->d_NDW, nr_rows_W, nr_cols_W, gpu_id);
                // update d_W
                ew_multiplication_launch(& nmf_handle->stream, nmf_handle->d_W, nmf_handle->d_NDW, nmf_handle->d_W, nr_rows_W, nr_cols_W, gpu_id);

                //--------------------------------------------------------------------------


                // Re-calulate V/(W@H) after W has been updated
                gpu_blas_mmul(& nmf_handle->handle_d, nmf_handle->d_wh, nmf_handle->d_W, nmf_handle->d_H, nr_rows_W, nr_cols_W, nr_cols_H, gpu_id);
                ew_division_launch(& nmf_handle->stream, nmf_handle->d_V, nmf_handle->d_wh, nmf_handle->d_vwh, nr_rows_V, nr_cols_V, gpu_id);


                //__________________________________________________________________________
                // H - Update
                //

                // numerator_h = transpose(w) @ (v / (w@h))
                gpu_blas_mmul_leftT(& nmf_handle->handle_e, nmf_handle->d_numerator_h, nmf_handle->d_W, nmf_handle->d_vwh, nr_cols_W, nr_rows_W, nr_cols_V, gpu_id);
                // denominator_h = transpose(w) @ ones
                gpu_blas_mmul_leftT(& nmf_handle->handle_f, nmf_handle->d_denominator_h, nmf_handle->d_W, nmf_handle->d_ones, nr_cols_W, nr_rows_W, nr_cols_ones, gpu_id);
                // numerator_h / denominator_h
                ew_division_launch(& nmf_handle->stream, nmf_handle->d_numerator_h, nmf_handle->d_denominator_h, nmf_handle->d_NDH, nr_rows_H, nr_cols_H, gpu_id);
                // update d_H
                ew_multiplication_launch(& nmf_handle->stream, nmf_handle->d_H, nmf_handle->d_NDH, nmf_handle->d_H, nr_rows_H, nr_cols_H, gpu_id);

                cudaStreamEndCapture(nmf_handle->stream, & nmf_handle->graph);
                cudaGraphInstantiate(& nmf_handle->instance, nmf_handle->graph, NULL, NULL, 0); 
                nmf_handle->graph_created = true;
            }
            cudaGraphLaunch(nmf_handle->instance, nmf_handle->stream);

        }
        cudaDeviceSynchronize();

        //__________________________________________________________________________
        // loss_converged
        //

        // calculate kl_loss
        if(iter == 0) {
            cudaStreamBeginCapture(nmf_handle->stream, cudaStreamCaptureModeGlobal);

            // create reconstruction: d_w @ d_h = d_w @ d_h
            gpu_blas_mmul(& nmf_handle->handle_g, nmf_handle->d_wh, nmf_handle->d_W, nmf_handle->d_H, nr_rows_W, nr_cols_W, nr_cols_H, gpu_id);

            // d_v / d_w @ d_h
            ew_division_launch(& nmf_handle->stream, nmf_handle->d_V, nmf_handle->d_wh, nmf_handle->d_vwh, nr_rows_V, nr_cols_V, gpu_id);
            
            //log(d_v / d_w @ d_h)
            ew_log_launch(& nmf_handle->stream, nmf_handle->d_vwh, nr_rows_V, nr_cols_V, gpu_id);

            //ew multiply v * log(d_v / d_w @ d_h)
            ew_multiplication_launch(& nmf_handle->stream, nmf_handle->d_V, nmf_handle->d_vwh, nmf_handle->d_vwh, nr_rows_V, nr_cols_V, gpu_id);

            // sum of v*log()
            sum_matrix_elements_launch(& nmf_handle->stream, nmf_handle->d_vwh, nmf_handle->d_sum, nr_rows_V, nr_cols_V, gpu_id);

            // sum of v
            sum_matrix_elements_launch(& nmf_handle->stream, nmf_handle->d_V, nmf_handle->d_Vsum, nr_rows_V, nr_cols_V, gpu_id);
            
            // sum of reconstruction (WH)
            sum_matrix_elements_launch(& nmf_handle->stream, nmf_handle->d_wh, nmf_handle->d_WHsum, nr_rows_V, nr_cols_V, gpu_id);

            //cudagraph
            cudaStreamEndCapture(nmf_handle->stream, & nmf_handle->kl_loss_graph);
            cudaGraphInstantiate(& nmf_handle->kl_loss_instance, nmf_handle->kl_loss_graph, NULL, NULL, 0); 
            
        }
        // launch cuda graph
        cudaGraphLaunch(nmf_handle->kl_loss_instance, nmf_handle->stream);
        cudaDeviceSynchronize();

        // copy sum to host
        cudaMemcpy((&sum), nmf_handle->d_sum, sizeof(double), cudaMemcpyDeviceToHost);
        // copy Vsum to host
        cudaMemcpy((&Vsum), nmf_handle->d_Vsum, sizeof(double), cudaMemcpyDeviceToHost);
        // copy WHsum to host
        cudaMemcpy((&WHsum), nmf_handle->d_WHsum, sizeof(double), cudaMemcpyDeviceToHost);

        // calculate kl_loss
        kl_loss = sum - Vsum + WHsum;

        loss_converged = false;
        if(iter == 0) {
            loss_init = kl_loss;
        }
        else if( ((prev_loss - kl_loss) / loss_init) < tolerance ) {
            loss_converged = true;
        }
        prev_loss = kl_loss;
        
        //--------------------------------------------------------------------------

        if(((iter + 1) >= (min_iterations / test_conv)) && loss_converged) {
            // Copy result into z -- returns result
            cudaMemcpy(
                W_ptr, nmf_handle->d_W, nr_rows_W*nr_cols_W*sizeof(dtype), cudaMemcpyDeviceToHost
            );

            // Copy result into z -- returns result
            cudaMemcpy(
                H_ptr, nmf_handle->d_H, nr_rows_H*nr_cols_H*sizeof(dtype), cudaMemcpyDeviceToHost
            );
            return (iter+1) * test_conv;
        } // check loss converged
    } // outer loop end

    // if it does not converge copy over results and print out max iterations
    // Copy result into z -- returns result
    cudaMemcpy(
        W_ptr, nmf_handle->d_W, nr_rows_W*nr_cols_W*sizeof(dtype), cudaMemcpyDeviceToHost
    );

    // Copy result into z -- returns result
    cudaMemcpy(
        H_ptr, nmf_handle->d_H, nr_rows_H*nr_cols_H*sizeof(dtype), cudaMemcpyDeviceToHost
    );
    return max_iterations;
} // fit_loop end



PYBIND11_MODULE(cuda_nmf_single, cudaBind) {
    cudaBind.doc() = "Performs a cuBLAS Matrix-Matrix multiply.";
    cudaBind.def("fit_loop", fit_loop);

    py::class_<ptr_wrapper<NMFHandle>>(cudaBind, "NMFHandle_t");
    cudaBind.def("NMFHandleNew", NMFHandleNew);
    cudaBind.def("NMFHandleCreate", NMFHandleCreate);
    cudaBind.def("NMFHandleDestroy", NMFHandleDestroy);
}