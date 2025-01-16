To run `read-fcc-higgs-v2.cpp` as a ROOT script, do

```
root "read-fcc-higgs-v2.cpp(\"FILENAME\", \"OUTPUT\")"
```

Note the quote symbols! `\"` is simply an escape character.

~~The boolean at the end is there for future purposes, but has no use right now. Still, this is required to prevent the potential error from ROOT.~~ Removed the third argument of the macro, 16 Jan 2025.