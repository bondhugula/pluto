#pragma scop
for (i = 0; i < N; ++i) {
  a->b = c;
}
c->d = a->b;
toto = a->d(); 
// Getter/Setter are not supported, only method calls.
// a->d() = bla bla is NOT ok
#pragma endscop
