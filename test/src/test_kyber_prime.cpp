
import std;
import lam.ctbignum;

using namespace lam::cbn;
using namespace lam::cbn::literals;

// Kyber parameter q = 3329
using KyberField = ZqElement<std::uint32_t, 3329>;

int main() {
    std::println("Testing Kyber Prime Field (q=3329) construction...");

    // 1. Literal construction
    constexpr auto a = KyberField{10_Z32};
    constexpr auto b = KyberField{3330_Z32}; // Should reduce to 1

    static_assert(a == KyberField{10_Z32});
    static_assert(b == KyberField{1_Z32});
    
    std::println("a = {}", a); 
    std::println("b = {}", b); 

    // 2. Arithmetic
    constexpr auto c = a + b; // 11
    constexpr auto d = a * a; // 100
    constexpr auto e = KyberField{1665_Z32} * KyberField{2_Z32}; // 3330 -> 1

    static_assert(c == KyberField{11_Z32});
    static_assert(d == KyberField{100_Z32});
    static_assert(e == KyberField{1_Z32});

    std::println("10 + 1 = {}", c);
    std::println("10 * 10 = {}", d);
    std::println("1665 * 2 = {}", e);

    // 3. Inverse
    constexpr auto two = KyberField{2_Z32};
    constexpr auto inv2 = KyberField{1_Z32} / two;
    
    std::println("1 / 2 = {}", inv2);
    static_assert(inv2 == KyberField{1665_Z32});

    std::println("Kyber Prime construction successful!");
    return 0;
}
