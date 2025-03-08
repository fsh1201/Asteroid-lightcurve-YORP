# Asteroid-lightcurve-YORP

## Requirements

gcc, g++, gfortran

## Usage

### Build

    $ make

### YORP Detection

    $ yorp lc input_parameter out_area out_par out_lcs yorp_chi2
**lc**: light curve file, the same format as DAMIT, see an [example lc](https://astro.troja.mff.cuni.cz/projects/damit/LightCurves/exportAllForAsteroid/208/plaintext/A208.lc.txt), and the [doc](https://github.com/fsh1201/Asteroid-lightcurve-YORP/blob/main/convexinv_doc.pdf)

**input_parameter**: the range of initial guesses of the spin state and YORP strength, degree and order of laplace series, scattering parameters, iteration stop condition. See the file [in](https://github.com/fsh1201/Asteroid-lightcurve-YORP/blob/main/in) for an example

**out_area**: areas and normal vectors

**out_par**: spin parameters obtained by iteration

**out_lcs**: output lightcurves, the same format as the **lc**

**yorp_chi2**: lambda(deg) beta(deg) p(hour) yorp(rad/d^2) rchisq jd0(JD) phi0(deg) a d k b c

### Uncertainty Estimation

    $ yorp_e lc input_parameter out_area out_par out_lcs yorp_chi2 nbootstrap nth

**nbootstrap**: bootstrap time, 8000+ would be great

**nth**: the number of processors to do this calculation

### YORP Detection & Uncertainty Estimation

    $ yorp_d_e path lcpath parameter
    
**path**: the path where you want to save the results

**lcpath**: the path of light curve file

**parameter**: the path of input parameters, an example see [in](https://github.com/fsh1201/Asteroid-lightcurve-YORP/blob/main/in)

### Areas to Shape

    $ cat out_area | minkowski | standardtri > shape.txt

    $ cat shape.txt | shape2obj > shape.obj

or

    $ cat out_area | minkowski | standardtri | shape2obj > shape.obj
