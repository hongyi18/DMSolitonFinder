# DMSolitonFinder

DMSolitonFinder is a Mathematica package devoted to solving dark matter soliton profiles in an automated way. It will dynamically adjust the boundary conditions and the
spatial range of solutions until localized solutions are found with the requested precision and accuracy. This package is [featured in the Staff Picks and Publication Materials columns on Wolfram Community](https://community.wolfram.com/groups/-/m/t/3203564).

## Installation/update

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

**Note:** If the installation/update fails, please check your internet connection by following the guide [Troubleshooting Internet Connectivity Problems](https://reference.wolfram.com/language/tutorial/TroubleshootingInternetConnectivity.html).

## Tutorial
Please follow the Mathematica notebook [Get started on DMSolitonFinder.nb](https://github.com/hongyi18/DMSolitonFinder/blob/main/Get%20started%20on%20DMSolitonFinder.nb) for a quick tutorial. 
For a more advanced example demonstrating the mass-radius relation for solitons with self-interactions, refer to [this post](https://community.wolfram.com/groups/-/m/t/3203564) on Wolfram Community. For definitions and equations used in the package, refer to the paper [arXiv:2406.05031](https://arxiv.org/abs/2406.05031).

<a href="https://github.com/hongyi18/DMSolitonFinder/blob/main/Get%20started%20on%20DMSolitonFinder.nb" download>Get started on DMSolitonFinder.nb</a>

## Common questions in finding optimal soliton profiles
**1. How can I find out what options are available for each function?** <br>
The usage message of each function has listed all available options. It can be called with "?FunctionName". For example, run 
```
?ShootFields
```

**2. How can I find soliton profiles with better convergence at large radii?** <br>
The convergence at large radii is controlled by the option `AmpTolerance` in function `ShootFields`. It is the reqiested ratio of the typical amplitudes of $f[r]$ at large radii, defined as the region between `OptionValue[AmpCheckRangeRatio]*rf` and `rf`, to the central amplitude. The default value of `AmpTolerance` is $10^{-4}$. You can try smaller values for this option for better convergence.

**3. The function `ShoofFields` has been stopped at very few loops or even at the beginning. What happened?** <br>
To start the shooting method, we always need an initial guess for boundary values. It could be that the initial guess for the value of $\Psi[r]$ at the inner boundary is too big. The default value of this guess is `InitialValuePsi->-1`. One can try taking smaller values than $-1$, or slightly changing the value of `InitialStepPsi`.

**4. What can I do if ``ShootFields`` fails to find soliton profiles?** <br>
One thing you can try is slightly changing the values of the options ``FirstMinimumLocationRatio`` or ``AmpCheckRangeRatio``. It could also be that Mathematica fails to solve the field equations accurately, so improving ``WorkingPrecision``, ``AccuracyGoal``, and ``PrecisionGoal`` might help.

**5. What can I do if soliton profiles are not accurate at large radii?** <br>
In addition to increasing ``WorkingPrecision``, ``AccuracyGoal``, and ``PrecisionGoal``, you can also try slightly changing the values of the options ``InitialMaxRadius`` and ``InitialStepRadius`` in ``ShootFields``, which control the initial range of solutions and the initial step used to increase the radius.

**6. What can I do if desired solutions are not found after the requested number of loops?** <br>
Try increasing `MaxIteration`, whose default value is $10^{3}$.

**7. How can I see the intermediate values of solutions at each loop, e.g., for debugging purpose?** <br>
This can be done by setting ``IfPrintProgress->1``. In this case, one can monitor the intermediate values of $\\{\Psi \text{i}, \text{d} \Psi, \text{rf}, \text{drf}\\}$ at each loop, which refers to the current value of $\Psi[r]$ at the inner boundary, current step for increasing $\Psi \text{i}$, current range of the solutions, and current step for increasing $\text{rf}$.

## Contact
If you find any bugs, please report them to my email address displayed on [my website](https://hongyi18.github.io/).
