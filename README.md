# micromat

Eigen-like library for block operations

## Requirements

- _Build and installation:_
    - Builds with C++03 and above.
    - Also builds with some microcontroller compilers.
    - No library dependency beside `libstdc++`
- _Runtime:_
    - There is no runtime requirement, not even an operating system.

## Features

This is a library for block based math.
It is written in a way, that allows it to be used on a PC or a microcontroller alike.
That's the reason why this library will only support types with a fixed maximum size.

Only a portion of Eigen's API is implemented.
Especially the const / mutable code duplication is only performed when needed.

Implemented API, see Eigen for meaning and usage:  
(e.t.) = lazy evaluated expression template

  * All Block objects:
    - `RowsAtCompileTime`, `ColsAtCompileTime`
    - `Constant`, `Ones`, `Zero` (e.t.)
    - `(row, col)`
    - `=`, `+=`, `-=`, `*=`, `/=` Also usable with Eigen expressions on the
    right hand side.
    - `+`, `-`, `*`, `/` (e.t.)
    - `setConstant`
    - `sum`
    - `abs`, `abs2` (e.t.)
    - `block` (e.t.)
    - `replicate` (e.t.)
  * `Array`
    - `MaxRowsAtCompileTime`, `MaxColsAtCompileTime`
  * `Map`

## Restrictions

The intention for this library is to be able to write code that runs on PC (with Eigen or micromat) and on microcontroller (with micromat) alike.
This imposes some requirements on the code written:

* the code may not depend on new, malloc or any other dynamic memory scheme.
* The API should be a subset of Eigen. There should not be additional features.

## Internals

* _curiously recurring template pattern_
  is a central design pattern used by micromat.
* _class BlockMixins_:
  This is a central class providing block operations.
  Deriving classes only need to implement access to indivitual members and size information.
* todo: Currently `+,-,*,/` only operate element-wise.
  This has to be reworked to support `Matrix` types.
* (todo: Expand this section when working at the code again.)

## Testing

The main test `TesttMatImpl.h` runs some basic operations and checks their results.
It has compile time switches to run either with `Eigen` or with `micromat`.
It runs either with `Eigen` or with `micromat` depending on a preprocessor macro.
This macro is set by `TestEigen.cpp` and `TestMicromat.cpp`.
With this trick we can prove API compatibility.
