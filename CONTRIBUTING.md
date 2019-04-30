
- Pluto now uses C99; it originally used ANSI C / C89 and so a lot of the 
  existing code is still in that form. The old code is typically updated when it 
  is touched.

- The LLVM coding style is used. Automatic formatting for this can be done
  using clang-format:

  $ clang-format -style=LLVM

  A .clang-format corresponding to this config exists in its top-level directory.

- Functions, structure members, and non-obvious variables should have comments.

- Contributions should be submitted as pull requests on the Github repo:
  https://github.com/bondhugula/pluto
