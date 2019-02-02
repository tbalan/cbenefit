# cbenefit 0.2.3

* Fixed an issue where row names could interfere with the algorithm
* Fixed an issue that arises when the number of treated and controls are not equal. A warning is issued if matching is done without replacement in this case.
* Added the `replace` argument, that overrides the other stuff in `matchit_args`.
* The call of the function is now also returned.
* Added a `NEWS.md` file to track changes to the package.
