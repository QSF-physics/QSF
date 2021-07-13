#pragma once
#define BIT(n) 1 << (n)
#define SET(x,y) x |= (BIT(y))
#define CLEAR(x,y) x &= ~(BIT(y))
// #define READ(x,y) ((0u == (x & (BIT(y))))?0u:1u)
#define READ(x,y) (0u == (x & (BIT(y))))
#define TOGGLE(x,y) (x ^= (BIT(y)))
#define CLEAR_LOWER(x) (x & ~(x/2))

#define STRINGIFY_(x) #x
#define STRINGIFY(x) STRINGIFY_(x)

#if defined(_MSC_VER)
#define FORCE_INLINE __forceinline
#else
#define FORCE_INLINE __attribute__((always_inline))
#endif