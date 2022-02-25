## Introduction
\
This module provides python functions for soil-property estimation and other soil-related processing, including runoff and erosion.

These functions are embedded within high-complexity watershed models like SWAT or APEX, but sometimes researchers may need these functions as stand-alone building blocks.



## Soil Structure

This sub-module contains functions for estimating several soil-structural paramters.
Most importantly, it gives several methods for calculating the USLE soil-erodibility factor (K). 



## Soil Water

This sub-module gives several methods for different soil water parameters. 
These functions can be combined to estimate a range of important parameters based on soil texture, bulk density and soil organic matter


## Runoff
Contains implementations of the NRCS curve-number method, including methods to estimate peak runoff rates.
Additionally, it contains several different implementations of the green and ampt method for runoff estimation.


## CN_GA_conversions
Contains a number of functions for putting the results of a curve-number estimate in terms of Green-Ampt and vice-versa.


##Erosion

Contains a few functions for estimating erosion, including the USLE and the MUSLE.


## Installation
\
To install, download this repository and run  

```pip freeze -r requirements.txt ```





## References:

Many functions here were adapted from the technical documentation and source code of the Soil and Water Assessment Tool (SWAT).

To better understand these functions, the SWAT technical guide is a useful resource: https://swat.tamu.edu/media/99192/swat2009-theory.pdf 

Where relevant, page numbers are included in function documentations.





Copyright ⓒ 2022 Benjamin Dube.
MIT License

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.


