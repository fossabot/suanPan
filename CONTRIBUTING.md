# Contribution

## Code Format

Please use the file `Resharper.DotSettings` to format the code. It is the configuration that can be imported into [ReSharper C++](https://www.jetbrains.com/resharper-cpp/).

The `.clang-format` generated can also be used to override format settings in various IDEs.

## Naming Convention

1. Please use `snake_case` instead of `CamelCase`.
2. Please use meaningful variable names instead of abstract names such as `a`, `tt`, etc.
3. Please avoid abbreviations, except for common ones such as `ptr` and `tmp`. For variables related to the corresponding theories, they can be spelled out. It does not hurt to define `epsilon_vol` instead of `ev` to represent volumetric strain variable. This is for the ease of readability and maintenance.

## Some Tips

1. Do not use variadic static variables.
2. Do not use raw pointers, instead, please use smart pointers.
3. Please always provide an example model that tests the code.
