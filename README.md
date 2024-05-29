# DMSolitonFinder

## Installation

**Method 1:** One can directly import the package from the web every time, using
```
Import["https://raw.githubusercontent.com/hongyi18/DMSolitonFinder/main/DMSolitonFinder.wl"]
```

**Method 2:** Download the package to the Applications subdirectory of $UserBaseDirectory, using
```
URLDownload["https://raw.githubusercontent.com/hongyi18/DMSolitonFinder/main/DMSolitonFinder.wl", FileNameJoin[{$UserBaseDirectory, "Applications", "DMSolitonFinder.wl"}]]
```
To load the package, simply run
```
<<DMSolitonFinder`
```

## Tutorial
Please follow the Mathematica notebook [Get started on DMSolitonFinder.nb](https://github.com/hongyi18/DMSolitonFinder/blob/main/Get%20started%20on%20DMSolitonFinder.nb) for a quick tutorial. For definitions and equations used in the package, refer to the paper [arXiv:123456]().

## Common questions in finding optimal soliton profiles
**1. How can I find out what options are available for each function?** <br>
The usage message of each function has listed all available options. It can be called with "?FunctionName". For example, run `?ShootFields`.

**2. How can I find soliton profiles with better convergence at large radii?** <br>
The convergence at large radii is controlled by the option `AmpTolerance` in function `ShootFields`. It is the reqiested ratio of the typical amplitudes of $f[r]$ at large radii, defined as the region between `OptionValue[AmpCheckRangeRatio]*rf` and `rf`, to the central amplitude. The default value of `AmpTolerance` is $10^{-4}$. You can try smaller values for this option for better convergence.

**3. The function `ShoofFields` is stopped at very few loops or even at the beginning. What happens?** <br>
To start the shooting method, we always need an initial guess for boundary values. It could be that the initial guess for the value of $\Psi[r]$ at the inner boundary is too big. The default value of this guess is `InitialValuePsi->-1`. One can try taking smaller values than $-1$, or slightly changing the value of `InitialStepPsi`.
