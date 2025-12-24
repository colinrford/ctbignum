import std;
import lam.ctbignum;

template<std::size_t Na, std::size_t Nb,typename T>
void mul_wrapper__internal(T* result, T* a, T* b)
{

  using lam::cbn::big_int;

  big_int<Na,T> A;
  big_int<Nb,T> B;

  std::memcpy(&(A[0]), a, Na * sizeof(T));
  std::memcpy(&(B[0]), b, Nb * sizeof(T));

  auto r = lam::cbn::mul(A,B);

  std::memcpy(result, &(r[0]), r.size() * sizeof(T));

}

// macro to define functions 
// (declared as extern "C", to be able to mangle names in our own simpler way)
#define MUL_WRAPPER(T_NAME, T_TYPE, N1, N2) extern "C" void mul_wrapper_##T_NAME##_##N1##_##N2(T_TYPE* r, T_TYPE* a, T_TYPE* b) \
{ mul_wrapper__internal<N1,N2>(r,a,b); }

MUL_WRAPPER(uint8_t, std::uint8_t, 1, 1); // the name of the function will be: mul_wrapper_uint8_t_1_1
MUL_WRAPPER(uint8_t, std::uint8_t, 2, 1); // etc
MUL_WRAPPER(uint16_t, std::uint16_t, 2, 2);

