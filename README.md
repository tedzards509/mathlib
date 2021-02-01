A (for now) header-only math lib for 64-bit(x86_64, aarch64,
risc-v NOT TESTED)

It is not designed to be light weight or easy to develop. It
should be fast. It is header-only making it easy to
integrate into different projects, build processes and
platforms.

please don't contribute to the bindings until the design is
thought through

A vector can represent

- A math vector for linear algebra

- A collection of values to iterate over like pixels

A matrix can represent

- A math matrix

Since matrices and complex numbers are always used
mathematically, they don't have integer counterparts.

it is designed for speed and usage of SSE, AVX, AVX512,
webSIMD, CUDA and NEON thus it makes sense to group pixels
into a vec64 since AVX512/VEC64 = 8-bit. Maybe I'll add an
Array Type that sorts this stuff automatically. Every type
should not have memory overhead. I have some ideas
concerning the direction to guide this library, but it
should grow as I find useful usages. This library is a c++
one, but I plan to add rust, Lua and Python bindings. It
supports multiple precisions like

- signed integer i8, i16, i32, i64
- unsigned integer u8, u16, u32, u64
- float at least 32, 64 bit

The feature level for each type(a combination of precision
and an amount like vec4f32) can differ, I don't use
templates since I'm not comfortable and don't know how to
detect the right vectorisation method.

"main.cpp" exists for testing purposes only.

Integrate the aml_lua_binding.cpp into your c++ application
if you want to. For simple or trivial functions(like
addition) the aml version is slower.

This library does NOT guarantee correctness of the result.
It ignores IEEE 754 edge cases(overflow, INF, NaN, div 0,
sqrt(-x), consistency, ...) in favor of speed. In the future
this library may add a safe math flag. Do not use this
library in safety critical mechanisms(encryption, aerospace,
automotive). This library has undefined behavior(e.g. an
internal division by 0) and it's not documented.