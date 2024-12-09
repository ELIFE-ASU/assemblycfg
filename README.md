# Context-Free Grammar (CFG) string assembly

Directed string assembly index calculator using smallest grammar algorithm re-pair. This will find a short valid path, but there is no guarantee that it will find the shortest possible valid path.

CFG/cfg_ai.py has two useful functions: ai_upper and ai_upper_with_pathways. Both return the same path length, but ai_upper_with_pathways prints the joins along the way. 

# Installation
Use pip to install this package.

`pip install git+https://github.com/ELIFE-ASU/CFG.git`

# Examples

See the examples folder for examples of how to use the package.

## Example 1: abracadabra.
```console
Processing abracadabra
START SYMBOLS: a,b,r,c,d
JOINS: 
a + b = ab
ab + r = abr
abr + a = abra
abra + c = abrac
abrac + a = abraca
abraca + d = abracad
abracad + abra = abracadabra
Path Length: 7
```
This example demonstrates its ability to find and reuse paths properly. 

## Example 2: aaaaaaa (7a's)
```console
Processing aaaaaaa
START SYMBOLS: a
JOINS: 
a + a = aa
aa + aa = aaaa
aaaa + aa = aaaaaa
aaaaaa + a = aaaaaaa
Path Length: 4
```
This is the classic example of LZ & breakage failing because it reuses implicitly built strings to create 7a in 3 steps. Since this algorithm guarantees valid assembly pathways, it will never do that. 
It properly adheres to the rules of assembly joining.

## Example 3: ababcdcd
```console
Processing ababcdcd
START SYMBOLS: a,b,c,d
JOINS: 
a + b = ab
c + d = cd
ab + ab = abab
abab + cd = ababcd
ababcd + cd = ababcdcd
Path Length: 5
```

This example demonstrates that this algorithm is not limited to depth = 1 as it can build 'ab' and 'cd' independently before eventually joining them. 
It is not restricted to building strings from left to right and is not limited by assembly width/depth.
