---
title: Know your PRNGs
author: gzt
date: '2019-03-09'
slug: know-your-prngs
categories:
  - C
tags:
  - R
  - statistics
  - computing
authors: []
---

TL;DR version: if all you're using the standalone `R` math library in C for is generating uniform
random numbers, I have a [little C program](https://github.com/gzt/replaceR) to remove 
that dependency as long as you don't mind your seeds not having the same output as `R`. I might
even fix that later.

I fell into a rabbit hole recently - I do a fair amount of work in C using the `R` standalone math
library for pseudo-random number generation (PRNGs) and sometimes some of its special functions
(rather than, say, working in C so it will be called from `R` or trying to implement a PRNG 
myself or a special function myself - leave RNGs and numerical analysis to experts). The usage is 
fairly simple - you call the standalone library, initialize your PRNG with `set_seed(seed1, seed2)`
 (this deserves some discussion, maybe later, it is an important topic), 
 then go ahead and call `runif()` or whatever other function you're interested in. 
 There are a few pitfalls but that's beyond the scope of this 
 post. Please feel free to reach out if you ever need help, though!

Anyway, I've been working on a cluster recently where it just seems to be a pain to get anything
installed on there. Rather than sort all that out, I looked at my code and realized all I was 
doing in this set of programs with the `R` math library was calling `runif()`. So, what would it
take to replace it? I was simply using the default and I didn't need perfect replication - 
that is to say, I didn't need to get the same results as `R` with the same seeds, I just needed
to be able to reproduce results if desired. I also want to be able to go back to the "real"
library if I need something more than just the PRNG, so I only want to have to change a line 
or two in a header to switch.

So, what do I need to do?

1. Find out what PRNG is used by default for `runif()` and how it is implemented. Hopefully
the C code in the the library is not too heavily `SEXP`ed up or not too different from the 
reference implementation of the method. I'm an okay programmer but these are things that you
don't want to mess with.
2. Find out how `set_seed()` and initialization work and write a replacement function.
3. Put it into my own code.


## So what PRNG are we using anyway

The helpful `R` documentation says the default method in `R` is the 
[Mersenne Twister](http://www.math.sci.hiroshima-u.ac.jp/~m-mat/MT/MT2002/emt19937ar.html) 
algorithm, which seems to be a *de facto* standard. The reference implementation is in 
[pretty clean C code](http://www.math.sci.hiroshima-u.ac.jp/~m-mat/MT/MT2002/CODES/mt19937ar.c) and seems fairly straightforward. Neat! There's a little warning on the site that
old methods of initialization aren't so good, so we'll have to see what `R` is doing and
be careful when deciding what to do ourselves. In the future I might want to replicate the
other methods that `R` offers but that adds some overhead.

The default method for the `R` standalone library called from C, though is the Marsaglia Multicarry.
I could implement that to completely duplicate the `Rmath.h` experience or I could duplicate the
actual `R` experience. I'm leaning toward the latter, as Marsaglia seems a little too easy!

Now to check the `R` code to [see how they implement this](https://svn.r-project.org/R/trunk/src/main/RNG.c).
They have a helpful comment about what they're doing.

```

/* ===================  Mersenne Twister ========================== */
/* From http://www.math.keio.ac.jp/~matumoto/emt.html */
/* New URL (accessed 2018-11-08):
   http://www.math.sci.hiroshima-u.ac.jp/~m-mat/eindex.html

   The initialization method in the 1998 code and paper had a minor
   issue that was addressed with new initialization approaches in an
   update in 2002.  R has always used a different initialization
   approach and is not affected by that issue.
*/

/* A C-program for MT19937: Real number version([0,1)-interval)
   (1999/10/28)
     genrand() generates one pseudorandom real number (double)
   which is uniformly distributed on [0,1)-interval, for each
   call. sgenrand(seed) sets initial values to the working area
   of 624 words. Before genrand(), sgenrand(seed) must be
   called once. (seed is any 32-bit integer.)
   Integer generator is obtained by modifying two lines.
     Coded by Takuji Nishimura, considering the suggestions by
   Topher Cooper and Marc Rieffel in July-Aug. 1997.

   Copyright (C) 1997, 1999 Makoto Matsumoto and Takuji Nishimura.
   When you use this, send an email to: matumoto@math.keio.ac.jp
   with an appropriate reference to your work.

   REFERENCE
   M. Matsumoto and T. Nishimura,
   "Mersenne Twister: A 623-Dimensionally Equidistributed Uniform
   Pseudo-Random Number Generator",
   ACM Transactions on Modeling and Computer Simulation,
   Vol. 8, No. 1, January 1998, pp 3--30.
*/
```

Nice, we don't need to worry about their initialization being out of date. 
Looking at their code, it looks like a simple port of the standard implementation. 
This is a relief, as I do not want to have to do any real coding. One nice thing:
They change a division into a multiplication.

```
    return ( (double)y * 2.3283064365386963e-10 ); /* reals: [0,1)-interval */
```
Maybe it doesn't matter in the end. #1 is done: the library's algorithm isn't really 
different from the reference implementation.

An important note here: a typical PRNG will generate a unsigned 32-bit `int` with approximately
uniform distribution in its range. `R` takes that `int` and turns it into a uniform `\([0,1)\)` by 
dividing by `\(2^32\)`. It uses these uniform values to generate other random variables and perform
sampling. Frankly, you might want to uniformly sample from a large set of integers more 
efficiently than that or generate uniform `double`s with more precision, so we might look into 
that. `R` has recently [begun fixing their default for sampling](http://r.789695.n4.nabble.com/Bias-in-R-s-random-integers-td4752563.html) 
from large sets of integers, so that's good.

## Initialization in R 

Pseudo-random number generators have some stored *state* and then they have two parts: 

1. Generating the next state.
2. Generating the next random number from the state.

Many generators put all their work in one or the other of these steps. If you know the state,
you know the entire future of the PRNG. If there is an easy function for inverting step 2,
you can find out the state of the PRNG by observing enough random numbers.
The Mersenne Twister has a state that is 623-dimensional and you obtain the next random number 
by reading the next entry of the state and doing a simple, invertible function of it. 
When you have gone through all 623 entries, it 
generates a new state (ie 623 new random numbers). The other PRNGs in `R` have smaller states.
However, you do not provide `R` with 623 32-bit integers, you provide it with 1 in `R` 
itself or 2 when calling it from C (note Marsaglia Multicarry has a state of size 2 and an identity function 
for step 2 as well). So, what are we doing?

In looking at the code, they aren't using the outdated initialization which the authors
of the MT algorithm warn against. It *looks* like they initialize by using the good ol' 
LCG (with a multiplier of `69069`, nice) with whatever (one) seed is provided to fill out 
the state space and then generate a new state based on that. 
`R` from `set_seed()` in C uses both values for initializing `runif()` with Marsaglia Multicarry,
and perhaps we should want to emulate this behavior rather than consider how MT is initialized.
For users of MT in R, this could be an issue for some, as, like, 
this means there are only `2^32` possible initializations, but it's not a big problem.
However, what are the odds we'll be providing an actual 32 bits of entropy to the generator?
It would be nice to be able to add in more entropy! The reference implementation also provides
a method for initializing with an arbitrary array of values. This is good.

Setting the seed can be a real issue if you want a reproducible simulation and/or are running a 
lot of simulations (say, you're operating in parallel). If you don't care about reproducing
your work and are only running one simulation, you can let `R` do its thing if you're working
interactively or, if you're in C, you can use some tried-and-true methods like generating 
something based on the time or your process ID or some hardware source or entropy. If you're 
doing something just for demonstration purposes and want to also provide some kind of time stamp,
I like to make a seed based on the date (`20190310` for today, for instance). When running
things in parallel, whether you want it to be completely reproducible or not, initialization is
going to require some care - basing something on the time is usually okay but if you're starting
a bunch of processes in near proximity to each other, there's a possibility of collisions and
you have to consider how the nearness of seeds is related to the results of your PRNG even if
they don't collide (spoiler: for some, including MT, this is a big problem). 

All this to say: I made this take both provided seeds in `set_seed()` into account as well as
provided a method for using the `init_by_array()` function. I also included a function that 
will initialize based on a hash of the arguments provided to call the program. This is handy
for me because I will often run several copies of the same program with only minor changes
in the arguments but I also want the results to be reproducible (at least while testing). 
The hash makes it so that each program will have a rather different seed even with only a 
minor change in the arguments but will have the same seed if I run it again with the exact 
same arguments.

## Using This

Include the header for `mt19337` and `replaceR` and you can use `set_seed` and `runif` just as
you did before. I provide a set of 
[reference results](https://github.com/gzt/replaceR/blob/master/ReplaceResults.txt) to check
whether your installation is working. 

You can compile and run the example file here:

```
gcc mt19937ar.c replaceR.c testreplaceR.c -o testreplaceR
./testreplaceR
```
To initialize with a hash of the arguments, you can call `hash_init_rand(argc, argv)` or
if you are combining it with a set of seeds, `hash_set_seed(argc, argv, seed1, seed2)`. 
As written, these are like the reference implementation and the library in that they use 
global variables for state. You can get around this.

## Future thoughts

The MT random number generator has a few difficulties: it fails statistical tests, has a 
large state, has problems when it comes to seeds that differ very little (problem for parallelism), 
and is slow for what it gives you. Also, observing a few hundred draws will
tell you everything you need to know about the state, which you may or may not care about.

There are other PRNGs that may or may not have these problems. In the future, I'll talk about
implementing them in `R` or what difficulties there are in doing that (I've been looking at
[PCG](https://www.pcg-random.org), for one). 

## Footnote: Marsaglia Multicarry

It's pretty easy code, it fits on three lines:

```
  /* 0177777(octal) == 65535(decimal)*/
      I1= 36969*(I1 & 0177777) + (I1>>16);
      I2= 18000*(I2 & 0177777) + (I2>>16);
      return fixup(((I1 << 16)^(I2 & 0177777)) * i2_32m1); /* in [0,1) */
```

`I1` is one seed, `I2` is another seed, `i2_32m1` is `\(2^{32}-1\)`. Figuring out what's going on there is an 
exercise, I guess.
