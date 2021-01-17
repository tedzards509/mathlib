A (for now) header-only math lib

please don't contribute for now

A vector can represent

- A math vector for linear algebra

- A collection of values to iterate over like pixels

it is designed for speed and usage of SSE, AVX, AVX512,
webSIMD, CUDA and NEON thus it makes sense to group pixels
into a vec64 since AVX512/VEC64 = 8-bit. Every type should
not have memory overhead. I have some ideas concerning the
direction to guide this library, but it should grow by
finding useful usages. This library is a c++ one, but I plan
to add rust, Lua and Python bindings. It supports multiple
precisions like

- signed integer i8, i16, i32, i64
- unsigned integer u8, u16, u32, u64
- float at least 32, 64, 80, 128 bit

The feature level for each type(a combination of precision
and an amount like vec4f32) can differ, I don't use
templates since I'm not comfortable and don't know how to
detect the right vectorisation method.

"main.cpp" exists for testing purposes only.

Integrate the aml_lua_binding.cpp into your c++ application
if you want to. For simple or trivial functions the aml
version is slower.