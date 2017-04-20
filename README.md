# venice

Author: Jean Coupon.

## Description

**[NEW]: fits support added for input files (faster reading, names for input columns, and [CFITSIO filters](http://heasarc.gsfc.nasa.gov/docs/software/fitsio/filters.html)).**

**WARNING: the input and output catalogue format is `fits` if read from STDIN or STDOUT, respectively (the format is set from the file extension, otherwise).**


`venice` is a mask utility program that reads a mask file (DS9 or fits type) and a catalogue of objects to:

1. create a pixelized mask,
2. find objects inside/outside a mask,
3. or generate a random catalogue of objects inside/outside a mask.

The program reads the mask file and checks if a point, giving its coordinates, is inside or outside the mask, i.e. inside or outside at least one polygon of the mask.

If the input mask is in DS9 `.reg` format, the following region types are accepted:
- `polygon`
- `box`
- `circle`
- `ellipse`

**The convention used is `0` when the object is INSIDE the mask and `1` when the object is OUTSIDE.**

## Installation

To compile it, you first need to install the gsl and cfitsio libraries (http://www.gnu.org/software/gsl/, http://heasarc.gsfc.nasa.gov/fitsio/).

Then run:
```shell
$ make
```
or
```shell
$ make PREFIX_GSL=DIRECTORY_NAME PREFIX_CFITSIO=DIRECTORY_NAME
```
if gsl and cfisio libraries are installed in different directories than `/usr/local`.


If you want to use a different compiler than gcc,
simply edit Makefile or type (only tested with gcc and icc):
```shell
$ make CC=my_favorite_compiler
```

**WARNING: venice without fits library is no longer supported (you must install the cfitsio library).**

## Usage

```
Usage: venice -m mask.[reg,fits]               [OPTIONS] -> binary mask for visualization
    or venice -m mask.[reg,fits] -cat file.cat [OPTIONS] -> objects in/out of mask
    or venice -m mask.[reg,fits] -r            [OPTIONS] -> random catalogue

Options:
    -r                       create random catalogue
    -cat FILE                input catalogue file name, default:stdin
    -o FILE                  output file name, default:stdout
    -catfmt [ascii,fits]     input catalogue format, default:fits if stdin
    -ofmt [ascii,fits]       output file format, default:fits if stdout
    -f [outside,inside,all]  output format, default:outside
    -[x,y]col N              column id for x and y (starts at 1)
    -coord [cart,spher]      coordinate type, default:cart
    -[x,y]min value          lower limit for x and y
    -[x,y]max value          upper limit for x and y
    -nz file_nz.in           redshift distribution for random objects
    -z zmin,zmax             redshift range for random objects (if volume limited)
    -seed  N                 random seed
    -npart N                 number of random objects
    -cd                      multiply npart by the mask area (for constant density)
    -flagName NAME           name of the flag colum for fits files. default: flag
    -h, --help               this message
```


IMPORTANT NOTICES:
- the general convention for venice is `0`:INSIDE the mask, `1`:OUTSIDE the mask,
- the gsl library is used to generate an improved random catalogue. In order to initialize venice with a different seed, set `-seed number`,
- for `.reg` masks, only ds9-types: `polygon`,`box`,`circle`,`ellipse` are supported,
- for `.fits` masks, the input catalogue should be given in image coordinates (x,y) and without the RA/DEC option; all points are returned with the pixel value added at the end of the line,
- if no output file is provided (`-o OUTFILE`), the result is prompted to the standard output (i.e. the terminal),
- the format for `file_nz.in` should be: `z n(z)` in histogram form. GSL convention: `bin[i]` corresponds to `range[i] <= x < range[i+1]`; the upper limit for the last bin is set to 100.
- if the input catalogue is given in fits format, the output catalogue also has to be in fits

## Options

### 1. Create a pixelized mask

```shell
$ venice -m mask[.reg,.fits] [OPTIONS]
```

For this task, don't add a catalogue in input (if you add a catalogue in input, e.g. `-cat file.cat`, the program will automatically switch to the task #2, see below for details). If the mask file is a DS9-type mask file (with the extension `.reg`), only ds9-types: `polygon`,`box`,`circle`,`ellipse` will be taken into account, see DS9 help for more details). If the mask is a `.fits` file, the convention is to consider object with value =0 as outside the mask, but the output will be in any case a "pixelized mask" in ASCII format with:
- 0 when the center of the pixel is inside the mask.
- 1 when the center of the pixel is outside the mask.

Command-line options:
- `-nx number`: number of pixels in the x direction. Default = 512.
- `-ny number`: number of pixels in the y direction. Default = 512.
- `-o outputfile`: output file where the pixelized mask is written [Default = mask.out].
- `-xmin value`: the minimum coordinate in the x direction.
- `-xmax value`: the maximum coordinate in the x direction.
- `-ymin value`: the minimum coordinate in the y direction.
- `-ymax value`: the maximum coordinate in the y direction. The default value for the coordinates limits are definied by the mask limits.

Example:
How to create a pixelized (10 X 10 pixels) mask with a mask file named mask[.reg,.fits] and put the results in pixel_mask.out:

```shell
$ venice -m mask[.reg,.fits] -nx 10 -ny 10 -o pixel_mask.out
```

The result in pixel_mask.out will look like this:

```
0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0
0.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  0.0
0.0  1.0  1.0  1.0  1.0  1.0  1.0  0.0  1.0  0.0
0.0  1.0  1.0  1.0  0.0  0.0  0.0  0.0  1.0  0.0
0.0  1.0  1.0  1.0  0.0  0.0  0.0  0.0  1.0  0.0
0.0  1.0  1.0  1.0  0.0  0.0  0.0  0.0  1.0  0.0
0.0  1.0  0.0  0.0  1.0  1.0  1.0  1.0  1.0  0.0
0.0  1.0  0.0  0.0  1.0  1.0  1.0  1.0  1.0  0.0
0.0  1.0  0.0  0.0  1.0  1.0  1.0  1.0  1.0  0.0
0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0
```

### 2. Find objects inside/outside a mask

For ds9-type masks, `venice` determines if a point is inside a polygon by drawing a line between a point and a second point (outside the polygon) and count how many times the line crosses the sides of the polygon. If the number is odd, the object is inside, if the number is even, the point is outside (Press et al 2007, Numerical recipes in c++).

In order to improve the speed, the process is made in 2 steps. A first quick check is performed to know if the point is inside or outside the square drawn by the extrema of the polygon (less time consuming) and then the above test is made.

This is fast. The only limitation is that all mask entries (`polygon`,`box`,`circle`,`ellipse`) are converted into polygons and stored in memory as a binary tree. A very big mask file (> 10^5 entries) will take a significant amount of memory.

The typical time for a CFHTLS wide catalogue and its terapix mask is about 5 seconds (200 000 objects in the catalogue and 1000 polygons in the mask).

One must provide a catalogue of objects with the same coordinate system as the mask file:

```shell
$ venice -m mask[.reg,.fits] -cat file.cat [OPTIONS]
```

The program checks if the objects in the catalogue are inside or outside the mask.

If `-f all` is set, the output file contains the initial catalogue with an additionnal column for the mask flag: the flag is 0 when the object is inside the mask and 1 when outside. In case of a fits file, the `-f` option is automatically set to `all` and the value of the pixel is given at the end of the line. Objects outside the mask limits are set to -99. A simple way to select objects with the desired pixel value (here "0") is to write:

```shell
$ venice -m mask.fits -cat old.cat | awk '$NF == 0 {print $0}' > new.cat
```

Options:
- `-f  [outside,inside,all]`: format of the catalogue. Default = outside. `outside` means only objects OUTSIDE the mask are written in the output file with the format: x y. `inside` means only objects INSIDE the mask are written. `all` means ALL objects are written in the catalogue with a flag 0 or 1 (resp. inside or outside) with the format: `col1 col2 ... flag`.
- `-xcol value`: column number of the x coordinate. Default=1.
- `-ycol value`: column number of the y coordinate. Default=2.
- `-o outputfile`: output file where the catalogue and the flag are written. Default = mask.out.
- `-xmin value`: the minimum coordinate in the x direction.
- `-xmax value`: the maximum coordinate in the x direction.
- `-ymin value`: the minimum coordinate in the y direction.
- `-ymax value`: the maximum coordinate in the y direction.
- `-flagName NAME`: for output fits file, the name of the flag column

Example:
How to create a new catalogue newcat.cat with a mask file named mask[.reg,.fits]
and a catalogue oldcat.cat for which the coordinate column numbers are
4 (x) and 5 (y):

```shell
$ venice -m mask[.reg,.fits] -cat oldcat.cat -xcol 4 -ycol 5 -o newcat.cat
```

### 3. Generate a random catalogue of objects inside/outside a mask

Given a mask file, the program generates a random catalogue and flag the objects if inside or outside the mask. The coordinates are drawn from a uniform distribution.

```shell
$ venice -m mask[.reg,.fits] -r [OPTIONS]
```

If `-f all` is set, the output file contains the initial catalogue with an additionnal column for the mask flag: the flag is 0 when the object is inside the mask and 1 when outside. In case of a fits file, the `-f` option is automatically set to `all` and the value of the pixel is given at the end of the line. Objects outside the mask limits are set to -99. A simple way to select objects with the desired pixel value (here "0") is to write:

```shell
$ venice -m mask.fits -r | awk '$NF == 0 {print $0}' > random.cat
```

Options:
- `-coord`: cart or spher. spher allows to draw a uniform distribution on a sphere. The coordinates must be in degrees, ra=[0:360],dec=[-90.0:90.0]. Default = cart.
- `-o outputfile`: output file where the catalogue is written. Default = mask.out.
- `-npart number`: number of objects. Default = 1 000 000.
- `-f  [outside,inside,all]`, format of the catalogue. Default = outside. `outside` means only objects OUTSIDE the mask are written in the output file with the format: x y. `inside` means only objects INSIDE the mask are written. `all` means ALL objects are written in the catalogue with a flag 0 or 1 (resp. inside or outside) with the format: `col1 col2 ... flag`.
- `-xcol value`: column number of the x coordinate. Default=1.
- `-ycol value`: column number of the y coordinate. Default=2.
- `-o outputfile`: output file where the catalogue and the flag are written. Default = mask.out.
- `-xmin value`: the minimum coordinate in the x direction.
- `-xmax value`: the maximum coordinate in the x direction.
- `-ymin value`: the minimum coordinate in the y direction.
- `-ymax value`: the maximum coordinate in the y direction.
- `-seed value`: (must be > 0), the seed for the GSL random generator (default: 20091982)
- `-nz file_nz`: to provide a file with the redshift distribution from which the random objects will be drawn. Note: if the binning is too small, this will "kill" the large scale power along the line of sight direction
- `-z zmin,zmax`: to have the random point number follow the volume size between zmin and zmax. This garantee a constant density as function of redshift. If the data sample is volume limited, this is the right option to use (instead of `-nz`).

The default value for the coordinates limits are definied by the mask limits. If you don't provide a mask (so if you only want a random catalogue with no mask), you have to define all of these values.

IMPORTANT: `npart` is the number of DRAWN objects, then if the mask is not empty the number of objects "outside" will be < npart. Tip: the ratio n_outside/npart is the unmasked area of your field (=1 if no mask).

Examples:

How to create a random catalogue random.cat (mask file mask[.reg,.fits]) that contains only objects outside the mask.

```shell
$ venice -m mask[.reg,.fits] -r -f outside -npart 10000000 -o random.cat -seed 1234567
```

How to create a random catalogue between 0.0 and 1.0 in both coordinates
with no mask, following a redshift distribution given in file_nz.in.

```shell
$ venice -r -xmin 0.0 -xmax 1.0 -ymin 0.0 -ymax 1.0 -nz file_nz.in
```
