- Make sure your code is indented the same way as surrounding code

- Use spaces instead of tabs for indenting (for vim, set expandtab). All of
  Pluto uses a 4-space indent, try using the same. For vim, it should
  be  "set ts=4 sw=4 expandtab"
- For for/while loops, place the opening brace on the same line.
``` 
for (i=0; i<N; i++) {

    }
```
- For function definitions, place opening brace on the next line: </br>
    `void foo(void) `</br>
    `{`</br></br>
    `}` 

- For if/else that have their blocks on a separate line, use braces (even if the
  block is a single line).
```
    if (...) {
        code;
    }else{
        code;
    }
```
- No need of braces if the then/else block appears on the same line</br>
    `if (...) code;`</br>
    `else more code;`

- Leave spaces in assignments for better readability. For eg.,  </br>
`   sol[j] = cst->val[i][j];`</br>
      instead of </br>
`   sol[j]=cst->val[i][j];`

- Leave space between arguments</br>
  `foo(a, b, c); `</br>
  `printf("%s\n", s);`</br>
      instead of</br>
  `foo(a,b,c);`</br>
  `printf("%s\n",s);`

- Any new functions added should at least have comments on what they do
